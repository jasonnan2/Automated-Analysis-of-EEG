classdef ERPanalysis_struct 
    properties
       DATA
       info
       analysisResults
    end
    
    
    methods
        
        function obj = ERPanalysis_struct(myStruct)
            disp('Initializing object')
            % obj is the formatted data structure with properties of DATA
            % and info
            % Initializing the properties to the object
            %-------------------------------------------------------------------%
            
            obj.DATA=myStruct.DATA;
            obj.info=myStruct.info;
            
        end
        
        function obj = cleanDatasets(obj)
            disp('Removing missing data')
            % Getting rid of the nan's and missing subjects
            for i = 1:length(obj.info.groupNames)
                obj.DATA.(obj.info.groupNames{i}) = cleanDataset(obj,obj.DATA.(obj.info.groupNames{i}));
            end
        end
        
        function data = cleanDataset(obj,data)
            % Get field names from the data structure
            fieldNames = obj.info.variables;

            % Loop over each field of the data structure
            for v = 1:numel(fieldNames)
                % Get the field data
                fieldData = data.(fieldNames{v});
                % Find missing data in the field data
                missingidx = find(isnan(squeeze(sum(sum(sum(fieldData, 1), 2), 3))));
                % Remove missing data from the field data
                fieldData(:,:,:,missingidx) = [];
                % Assign the cleaned field data back to the data structure
                data.(fieldNames{v}) = fieldData;
                data.missingSubs{v}=[fieldNames{v} data.subList(missingidx)];
            end
        end
        

        function obj = standardPipeline(obj,baselineTime)
            % baseline time | vector with start and end time of baseline in mS 
            [~,baselineIDX]=min(abs(baselineTime'-obj.info.timeAxis)');
            obj.info.baselineIDX=baselineIDX;
            properties = obj.info.variables;
            disp('---------------------------------')
            for j = 1:length(properties)
                property = properties{j};
                disp("Processing "+property)
                combined = combine_groups(obj,property);
                cleaned = rej5SD(combined);
                baselineCorrected = baselineCorrection(cleaned,obj.info.baselineIDX);
                obj = split_combined(obj, baselineCorrected,property);
                disp('---------------------------------')
            end
        end
        
        function combined = combine_groups(obj,property)
            disp('Combining groups')
            % stacking groups together for 5SD
            data = cell(1, length(obj.info.groupNames)); % empty cell array
            for i = 1:length(data)
                data{i} = obj.DATA.(obj.info.groupNames{i}).(property);
            end
            combined =cat(4,data{:});
        end

        function obj = split_combined(obj, groupedData,property)
            disp('Splitting groups')
            % Reordering the data back into groups in obj
            for i=1:length(obj.info.groupNames)
                obj.DATA.(obj.info.groupNames{i}).(property)=groupedData(:,:,:,1:size(obj.DATA.(obj.info.groupNames{i}).(property),4));
                groupedData(:,:,:,1:size(obj.DATA.(obj.info.groupNames{i}).(property),4))=[];
            end
        end
        
        function obj=plotERP(obj,timeRange,properties)
            obj.info.timeRange=timeRange;
            timeNames = fieldnames(timeRange);
            N = numel(fieldnames(obj.DATA));
            combinations = nchoosek(1:N,2); % This will give you a matrix where each row is a combination of two groups

            for p=1:length(properties)
                property=properties{p};
                for comb = 1:size(combinations, 1)
                    firstSub=combinations(comb, 1);
                    secondSub=combinations(comb, 2);

                    group1 = obj.DATA.(obj.info.groupNames{firstSub}).(property);
                    group2 = obj.DATA.(obj.info.groupNames{secondSub}).(property);

                    figure
                    figCount=0;
                    for t=1:length(timeNames)
                        [~, startIdx] = min(abs(obj.info.timeAxis - timeRange.(timeNames{t})(1)));
                        [~, endIdx] = min(abs(obj.info.timeAxis - timeRange.(timeNames{t})(2)));
                        
                        obj.info.timeIDX.(timeNames{t})=[startIdx endIdx];
         
                        for f=1:length(obj.info.freq_list)
                            figCount = figCount + 1;
                            subplot(3,4,figCount)
                            s1=squeeze(nanmean(group1(f,:,startIdx:endIdx,:),3)); % get to chan x sub
                            s2=squeeze(nanmean(group2(f,:,startIdx:endIdx,:),3)); % get to chan x sub
                            
                            chanlocs=topoSignificance(s1,s2,obj.info.chanlocs,{'.','+'});
                            plotSigTopo(s1,s2,chanlocs,{'.','+'})
                            
                            if t==1
                                title(obj.info.freq_list{f})
                            end
                            if f==1
                                text(-1,-.2, timeNames{t},'fontweight','bold','fontsize',16,'rotation',90);         
                            end
                        end
                    end
                    sgtitle(property+" "+obj.info.groupNames(secondSub)+"-"+obj.info.groupNames(firstSub))
                end
            end
        end
    end
end

function combined = rej5SD(combined)
    disp('5 SD Outlier')
    for f = 1:size(combined, 1)
        for chan = 1:size(combined, 2)
            for t = 1:size(combined, 3)
                zval = zscore(squeeze(combined(f, chan, t, :)));
                combined(f, chan, t, zval > 5) = nan;
            end
        end
    end
end

function baselineCorrected = baselineCorrection(cleaned,baselineIDX)
    disp('Baseline Correction')
    baselineCorrected=nan(size(cleaned));
    for f=1:size(cleaned,1)
        for chan=1:size(cleaned,2)
            for sub=1:size(cleaned,4)
                data=squeeze(cleaned(f,chan,:,sub));
                baselineData=nanmean(data(baselineIDX(1):baselineIDX(end)));
                baselineCorrected(f,chan,:,sub)=data-baselineData;
            end
        end
    end
end

%% Plotting scalpmap with significant electrodes
function new_chan = plotSigTopo(s1,s2,chanlocs,key)
    % s1 and s2 are chan x sub size
    plotdata=nanmean(s2,2)-nanmean(s1,2);
    topoplot(plotdata,chanlocs,'headrad','rim','electrodes','labels','efontsize' ,16); 
    colorbar(); 
    caxis([min(plotdata),max(plotdata)])
    colormap(parula(128))
    brighten(0.5)
end

%% Function to calculate which electrodes are significant btw two groups
function chanlocs=topoSignificance(s1,s2,chanlocs,key)

    % pre and post are data matrices size N subjects x C channels
    % key is a two element cell array with markers for non-sig and sig
    % new_chan is the chan_locs
    % Perform t-tests
    
    s1=squeeze(s1);
    s2=squeeze(s2);
    for chan = 1:length(chanlocs)
        % Get non-NaN indices
        [~,p(chan)] = ttest2(s1(chan,:), s2(chan,:));
        if isnan(p(chan))
            disp('asdfawrsasdbaweasbas')
            waitforbuttonpress
        end
        
    end
    %[~,p]=fdr(reshape(p,[1,numel(p)]),0.05);
    mask=key((p<0.05)+1);
    [chanlocs.labels]=mask{:};
end
        







    