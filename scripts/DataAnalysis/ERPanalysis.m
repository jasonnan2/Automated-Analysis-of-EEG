classdef ERPanalysis
    properties
       scalpData
       sourceData
       info
       analysisResults
    end
    
    methods
        
        function obj = ERPanalysis(myStruct,baselineTime, timeRange)
            % obj           | formatted data structure with properties of DATA
            % baseline time | vector with start and end time of baseline in mS 
            disp('Initializing object')
            
            % and info
            % Initializing the properties to the object
            %-------------------------------------------------------------------%
            
            obj.scalpData=myStruct.scalpData;
            if isfield(myStruct,'sourceData')
                obj.sourceData=myStruct.sourceData;
            end
            obj.info=myStruct.info;
            [~,baselineIDX]=min(abs(baselineTime'-obj.info.timeAxis)');
            obj.info.baselineIDX=baselineIDX;
            obj.info.timeRange = timeRange;
            timeNames = fieldnames(obj.info.timeRange);
            for t = 1:length(timeNames)
                [startIdx, endIdx] = obj.getTimeIndices(timeRange, timeNames{t});
                obj.info.timeIDX.(timeNames{t}) = [startIdx endIdx];
            end
        end
        
        function [startIdx, endIdx] = getTimeIndices(obj, timeRange, timeName)
            [~, startIdx] = min(abs(obj.info.timeAxis - timeRange.(timeName)(1)));
            [~, endIdx] = min(abs(obj.info.timeAxis - timeRange.(timeName)(2)));
        end
        
        function obj = cleanDatasets(obj,runType)
            disp('Removing missing data')
            % Getting rid of the nan's and missing subjects
            for i = 1:length(obj.info.groupNames)
                obj.(runType).(obj.info.groupNames{i}) = obj.cleanDataset(obj.(runType).(obj.info.groupNames{i}));
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
        

        function obj = standardPipeline(obj,runType)
            properties = obj.info.variables;
            disp('---------------------------------')
            for j = 1:length(properties)
                property = properties{j};
                disp("Processing "+property)
                combined = obj.combine_groups(runType,property);
                cleaned = rej5SD(combined);
                baselineCorrected = baselineCorrection(cleaned,obj.info.baselineIDX);
                obj = obj.split_combined( baselineCorrected,property,runType);
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
        
        function obj=plotScalpmap(obj,properties)
            
            timeNames = fieldnames(obj.info.timeIDX);
            N = numel(fieldnames(obj.scalpData));
            combinations = nchoosek(1:N,2); % This will give you a matrix where each row is a combination of two groups

            for p=1:length(properties)
                property=properties{p}; % get property name
                for comb = 1:size(combinations, 1)

                    group1=obj.info.groupNames{combinations(comb, 1)};
                    group2=obj.info.groupNames{combinations(comb, 2)};

                    figure
                    figCount=0;
                    for t=1:length(timeNames)
                        timeName=timeNames{t};
                        for f=1:length(obj.info.freq_list)
                            freq=obj.info.freq_list{f};

                            figCount = figCount + 1;
                            subplot(3,4,figCount)

                            s1 = squeeze(nanmean(obj.getGroupData(group1,property,freq,timeName),3));
                            s2 = squeeze(nanmean(obj.getGroupData(group2,property,freq,timeName),3));

                            pvals=calGroupSig(s1,s2); % get 1x n channel of p values 
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
                    sgtitle(property+" "+group2+"-"+group1)
                end
            end
        end
        
        function plotERPs(obj,vars2plot,freq2plot,times2plot)
            color_list={'g','b','r'};
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
                                hold on
                                h = gobjects(N, 1); % Preallocate an array for the line handles
                                for n=1:N
                                    d=obj.scalpData.(obj.info.groupNames{n}).(property)(freqIdx,elecIdxs,timeIdx(1):timeIdx(2),:);
                                    data(:,n) = squeeze(nanmean(d,4));
                                    eb = shadedErrorBar(obj.info.timeAxis(timeIdx(1):timeIdx(2)),data(:,n),std(squeeze(d)')/sqrt(size(d,4)), color_list{n});
                                    h(n) = eb.mainLine;
                                    set(get(get(eb.patch, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
                                end

                                hold off
                                %plot(obj.info.timeAxis(timeIdx(1):timeIdx(2)),data)

                                title(chan)
                                if c==1
                                    ylabel(freq,'fontweight','bold','fontsize',16)
                                end
                            end
                        end
                        Lgnd = legend(h,obj.info.groupNames);
                        Lgnd.Position(1) = 0.7;
                        Lgnd.Position(2) = 0.9;
                        sgtitle(strcat(property,'-',time))
                    end
                end
            end
        end
        function plotNetwork(obj,netwrk,properties)
            timeNames = fieldnames(obj.info.timeRange);
            N = length(obj.info.groupNames);
            combinations = nchoosek(1:N,2);
            for p=1:length(properties)
                property=properties{p};
                figure
                count=0;
                for t=1:length(timeNames)
                    timeName=timeNames{t};
                    timeIdx = obj.info.timeIDX.(timeName);
                    for f = 1:length(obj.info.freq_list)
                        freq=obj.info.freq_list{f};
                        % Generating plot data
                        data=[];sem=[];
                        for n=1:N
                            for i=1:length(netwrk)

                                subdata = squeeze(nanmean(nanmean(getGroupData(obj,obj.info.groupNames{n},property,freq,timeName,netwrk(i).roi),2),3));
                                data(i,n) = nanmean(subdata,'all');
                                sem(i,n) = std(subdata)/sqrt(length(subdata));
                            end
                        end
                        count=count+1;
                        subplot(length(timeNames),length(obj.info.freq_list),count)
                        b=bar(data,'grouped'); hold on
                        % Find the number of groups and the number of bars in each group
                        [ngroups, nbars] = size(data);
                        % Calculate the width for each bar group
                        groupwidth = min(0.8, nbars/(nbars + 1.5));
                        for i = 1:nbars
                            % Calculate center of each bar
                            x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
                            errorbar(x, data(:,i), sem(:,i), 'k', 'linestyle', 'none');
                        end
                        for comb = 1:size(combinations, 1)
                            group1=obj.info.groupNames{combinations(comb, 1)};
                            group2=obj.info.groupNames{combinations(comb, 2)};
                            s1=[];s2=[];
                            for i=1:length(netwrk)
                                % collapsing networks
                                s1(i,:) = squeeze(nanmean(nanmean(getGroupData(obj,group1,property,freq,timeName,netwrk(i).roi),2),3));
                                s2(i,:) = squeeze(nanmean(nanmean(getGroupData(obj,group2,property,freq,timeName,netwrk(i).roi),2),3));
                            end
                            %[~,pvals(comb,:)] = ttest2(s1, s2);
                            pvals(comb,:)=calGroupSig(s1,s2);

                            group1Location = (1:ngroups) - groupwidth/2 + (2*combinations(comb, 1)-1) * groupwidth / (2*nbars);
                            group2Location = (1:ngroups) - groupwidth/2 + (2*combinations(comb, 2)-1) * groupwidth / (2*nbars);
                            A=[group1Location;group2Location]';
                            groupingKey = mat2cell(A, ones(1, size(A, 1)), size(A, 2));
                            sigstar(groupingKey,pvals(comb,:))
                        end
                        set(gca,'xticklabel',{netwrk.name},'fontweight','bold','fontsize',12,'Xticklabelrotation',90)%This line should replace the numbers with your categorical label%errorbar(CVmall,CVsemall);
                        if t==1
                            title(freq)
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
        
%         function s = getGroupData(obj,group,property,freq,timeName,chans)
%             % group | name of the group
%             % property | string of variable
%             % freq | string of frequency band
%             % chans | vector of chans to include, if 'all' or empty indicates use all channels
%             % timeName | string of timeName
% 
%             if nargin < 6
%                 chans = 'all';
%             end
% 
%             freqIdx = find(strcmp(obj.info.freq_list, freq));
%             timeIdx = obj.info.timeIDX.(timeName);
% 
%             if isempty(chans) || (ischar(chans) && strcmp(chans, 'all'))
%                 chans = 1:size(obj.scalpData.(group).(property), 2); % assuming channels are the 2nd dimension
%             end
%             s = obj.scalpData.(group).(property)(freqIdx,chans,timeIdx(1):timeIdx(2),:);
%         end
        
        
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
function pvals=calGroupSig(s1,s2)

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
        





    