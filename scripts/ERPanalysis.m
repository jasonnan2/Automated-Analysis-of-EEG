classdef ERPanalysis
    properties
       scalpData
       sourceData
       info
       analysisResults
    end
    
    
    methods
        
        function obj = ERPanalysis(myStruct)
            disp('Initializing object')
            % obj is the formatted data structure with properties of DATA
            % and info
            % Initializing the properties to the object
            %-------------------------------------------------------------------%
            
            obj.scalpData=myStruct.scalpData;
            
            if isfield(myStruct,'sourceData')
                obj.sourceData=myStruct.sourceData;
            end
            obj.info=myStruct.info;
            
            %obj.analysisResults=struct;
            
        end
        
        function obj = cleanDatasets(obj,runType)
            disp('Removing missing data')
            % Getting rid of the nan's and missing subjects
            for i = 1:length(obj.info.groupNames)
                obj.(runType).(obj.info.groupNames{i}) = cleanDataset(obj,obj.(runType).(obj.info.groupNames{i}));
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
        

        function obj = standardPipeline(obj,baselineTime,runType)
            % baseline time | vector with start and end time of baseline in mS 
            [~,baselineIDX]=min(abs(baselineTime'-obj.info.timeAxis)');
            obj.info.baselineIDX=baselineIDX;
            properties = obj.info.variables;
            disp('---------------------------------')
            for j = 1:length(properties)
                property = properties{j};
                disp("Processing "+property)
                combined = combine_groups(obj,runType,property);
                cleaned = rej5SD(combined);
                baselineCorrected = baselineCorrection(cleaned,obj.info.baselineIDX);
                obj = split_combined(obj, baselineCorrected,property,runType);
                disp('---------------------------------')
            end
        end
        
        function combined = combine_groups(obj,runType,property)
            disp('Combining groups')
            % stacking groups together for 5SD
            data = cell(1, length(obj.info.groupNames)); % empty cell array
            for i = 1:length(data)
                data{i} = obj.(runType).(obj.info.groupNames{i}).(property);
            end
            combined =cat(4,data{:});
        end

        function obj = split_combined(obj, groupedData,property,runType)
            disp('Splitting groups')
            % Reordering the data back into groups in obj
            for i=1:length(obj.info.groupNames)
                obj.(runType).(obj.info.groupNames{i}).(property)=groupedData(:,:,:,1:size(obj.(runType).(obj.info.groupNames{i}).(property),4));
                groupedData(:,:,:,1:size(obj.(runType).(obj.info.groupNames{i}).(property),4))=[];
            end
        end
        
        function obj=plotScalpmap(obj,timeRange,properties)
            obj.info.timeRange=timeRange;
            timeNames = fieldnames(timeRange);
            N = numel(fieldnames(obj.scalpData));
            combinations = nchoosek(1:N,2); % This will give you a matrix where each row is a combination of two groups

            for p=1:length(properties)
                property=properties{p};
                for comb = 1:size(combinations, 1)
                    firstSub=combinations(comb, 1);
                    secondSub=combinations(comb, 2);

                    group1 = obj.scalpData.(obj.info.groupNames{firstSub}).(property);
                    group2 = obj.scalpData.(obj.info.groupNames{secondSub}).(property);

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
            N=numel(fieldnames(obj.scalpData));
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
                                    data(:,n) = squeeze(nanmean(obj.scalpData.(obj.info.groupNames{n}).(property)(freqIdx,elecIdxs,timeIdx(1):timeIdx(2),:),4));
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
        function plotNetwork(obj,netwrk,timeRange,properties)
            obj.info.timeRange=timeRange;
            timeNames = fieldnames(obj.info.timeRange);
            N = length(obj.info.groupNames);
            combinations = nchoosek(1:N,2);
            for p=1:length(properties)
                property=properties{p};
                figure
                count=0;
                for t=1:length(timeNames)
                    time=timeNames{t};
                    [~, startIdx] = min(abs(obj.info.timeAxis - timeRange.(time)(1)));
                    [~, endIdx] = min(abs(obj.info.timeAxis - timeRange.(time)(2)));
                    obj.info.timeIDX.(time)=[startIdx endIdx];
                    timeIdx = obj.info.timeIDX.(time);
                    for f = 1:length(obj.info.freq_list)
                        % Generating plot data
                        data=[];sem=[];
                        for n=1:N
                            for i=1:length(netwrk)
                                subdata = squeeze(nanmean(nanmean(obj.sourceData.(obj.info.groupNames{n}).(property)(f,netwrk(i).roi,timeIdx(1):timeIdx(2),:),2),3));
                                data(i,n) = nanmean(subdata,'all');
                                sem(i,n) = std(subdata)/sqrt(length(subdata)); % Standard error of mean across N=51 subjects
                            end
                        end

                        count=count+1;
                        subplot(length(timeNames),length(obj.info.freq_list),count)
                        b=bar(data,'grouped'); hold on
                        % Find the number of groups and the number of bars in each group
                        [ngroups, nbars] = size(data);
                        % Calculate the width for each bar group
                        groupwidth = min(0.8, nbars/(nbars + 1.5));
                        % Set the position of each error bar in the centre of the main bar
                        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
                        for i = 1:nbars
                            % Calculate center of each bar
                            x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
                            errorbar(x, data(:,i), sem(:,i), 'k', 'linestyle', 'none');
                        end

                        for comb = 1:size(combinations, 1)
                            firstSub=combinations(comb, 1);
                            secondSub=combinations(comb, 2);

                            for i=1:length(netwrk)
                                group1 = squeeze(nanmean(nanmean(obj.sourceData.(obj.info.groupNames{firstSub}).(property)(f,netwrk(i).roi,timeIdx(1):timeIdx(2),:),2),3));
                                group2 = squeeze(nanmean(nanmean(obj.sourceData.(obj.info.groupNames{secondSub}).(property)(f,netwrk(i).roi,timeIdx(1):timeIdx(2),:),2),3));

                                [~,pvals(comb,i)] = ttest2(group1, group2);

                            end
                            group1Location = (1:ngroups) - groupwidth/2 + (2*firstSub-1) * groupwidth / (2*nbars);
                            group2Location = (1:ngroups) - groupwidth/2 + (2*secondSub-1) * groupwidth / (2*nbars);
                            A=[group1Location;group2Location]';
                            groupingKey = mat2cell(A, ones(1, size(A, 1)), size(A, 2));
                            sigstar(groupingKey,pvals(comb,:))
                        end

                        set(gca,'xticklabel',{netwrk.name},'fontweight','bold','fontsize',12,'Xticklabelrotation',90)%This line should replace the numbers with your categorical label%errorbar(CVmall,CVsemall);
                        if t==1
                            title(obj.info.freq_list{f})
                        end
                        if f==1
                           ylabel(timeNames{t},'fontweight','bold','fontsize',16,'rotation',90);         
                        end

                     end
                end
                Lgnd = legend(obj.info.groupNames);
                Lgnd.Position(1) = 0.7;
                Lgnd.Position(2) = 0.9;
                sgtitle(property)
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
        





    