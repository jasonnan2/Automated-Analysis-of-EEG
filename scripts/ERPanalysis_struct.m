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
            obj.analysisResults=struct;
            
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
        
        function obj=plotScalpmap(obj,timeRange,properties)
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
                            s1=squeeze(nanmean(group1(f,:,startIdx:endIdx,:),3)); % get to chan x sub of first subject
                            s2=squeeze(nanmean(group2(f,:,startIdx:endIdx,:),3)); % get to chan x sub of second subject
                            
                            pvals=topoSignificance(s1,s2); % get 1x n channel of p values 
                            chanlabels={obj.info.chanlocs.labels}; % just the cell array of labels
                            
                            chans2add = chanlabels(pvals<0.05); % significant channels to add to list
                
                            % Check if the field already exists
                            if isfield(obj.analysisResults, 'sigElectrodes') && isfield(obj.analysisResults.sigElectrodes, property) && ...
                                isfield(obj.analysisResults.sigElectrodes.(property), timeNames{t}) && isfield(obj.analysisResults.sigElectrodes.(property).(timeNames{t}), obj.info.freq_list{f})
                            
                                % Get the existing cell array
                                existingArray = obj.analysisResults.sigElectrodes.(property).(timeNames{t}).(obj.info.freq_list{f});
                                % Check if the new letter is different from the existing ones
                                chans2add = setdiff(chans2add, existingArray);
                                obj.analysisResults.sigElectrodes.(property).(timeNames{t}).(obj.info.freq_list{f}) = [existingArray, chans2add];
                            elseif ~isempty(chans2add)
                                % Create a new cell array with the new letter
                                obj.analysisResults.sigElectrodes.(property).(timeNames{t}).(obj.info.freq_list{f}) = chans2add;
                            end

                            plotSigTopo(s1,s2,obj.info.chanlocs,pvals,{'.','+'})
                            
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
        
        function plotERPs(obj,vars2plot,freq2plot,times2plot)
            N=numel(fieldnames(obj.DATA));
            timeNames = fieldnames(obj.info.timeRange);
            sigElectrodes=obj.analysisResults.sigElectrodes;
            for p=1:length(vars2plot)
                property=vars2plot{p};
                time_list = intersect(fieldnames(sigElectrodes.(property)),times2plot);

                for t=1:length(time_list)

                    time=time_list{t};
                    sig=sigElectrodes.(property).(time);
                    [~,~,ib] = intersect(fieldnames(sig),freq2plot);
                    freq_list=freq2plot(sort(ib));

                    % get longest electrode map
                    maxLength = 0;
                    for i = 1:length(freq_list)
                        currentLength = length(sig.(freq_list{i}));  % Get the length of the current cell array
                        if currentLength > maxLength
                            maxLength = currentLength;  % Update the maximum length
                        end
                    end

                    if ~isempty(freq_list)
                        figure
                        for f=1:length(freq_list)

                            for c=1:length(sig.(freq_list{f}))

                                subplot(length(freq_list),maxLength,(f-1)*maxLength+c)

                                freq=freq_list{f};
                                chan=sig.(freq){c};
                                freqIdx = find(strcmp(obj.info.freq_list, freq));
                                timeIdx = obj.info.timeIDX.(time);
                                elecIdxs = find(strcmp({obj.info.chanlocs.labels},chan ));
                                data=[];
                                for n=1:N
                                    data(:,n) = squeeze(nanmean(obj.DATA.(obj.info.groupNames{n}).(property)(freqIdx,elecIdxs,timeIdx(1):timeIdx(2),:),4));
                                end

                                plot(obj.info.timeAxis(timeIdx(1):timeIdx(2)),data)
                                title(chan)
                                if c==1
                                    ylabel(freq,'fontweight','bold','fontsize',16)
                                end
                            end
                        end
                        Lgnd = legend(obj.info.groupNames);
                        Lgnd.Position(1) = 0.7;
                        Lgnd.Position(2) = 0.9;
                        sgtitle(strcat(property,'-',time))
                    end
                end
            end

        end
        
    end
end
%%
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
%%
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
function new_chan = plotSigTopo(s1,s2,chanlocs,pvals, key)

    %[~,p]=fdr(reshape(p,[1,numel(p)]),0.05);
    mask=key((pvals<0.05)+1);
    [chanlocs.labels]=mask{:};
    
    
    % s1 and s2 are chan x sub size
    plotdata=nanmean(s2,2)-nanmean(s1,2);
    topoplot(plotdata,chanlocs,'headrad','rim','electrodes','labels','efontsize' ,16); 
    colorbar(); 
    caxis([min(plotdata),max(plotdata)])
    colormap(parula(128))
    brighten(0.5)
end



%% Function to calculate which electrodes are significant btw two groups
function pvals=topoSignificance(s1,s2)

    % s1 and s2 are data matrices size C channels x N subjects
    % key is a two element cell array with markers for non-sig and sig
    % new_chan is the chan_locs
    % Perform t-tests
    
    s1=squeeze(s1);
    s2=squeeze(s2);
    for chan = 1:size(s1,1)
        % Get non-NaN indices
        [~,pvals(chan)] = ttest2(s1(chan,:), s2(chan,:));
        if isnan(pvals(chan))
            disp('asdfawrsasdbaweasbas')
            waitforbuttonpress
        end
    end

end
        





    