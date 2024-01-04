classdef ScalpAnalysis < DataAnalysis
    properties
        DATA
        scalpResults
    end
    
    methods
        function obj = ScalpAnalysis(DATA, info, baselineTime, timeRange)
            % Call superclass constructor
            obj = obj@DataAnalysis(info, baselineTime, timeRange);
            obj.DATA=DATA;
        end
        
        function obj=plotScalpmap(obj,properties)
            
            timeNames = fieldnames(obj.info.timeIDX);
            N = numel(fieldnames(obj.DATA));
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

                            pvals=obj.calGroupSig(s1,s2); % get 1x n channel of p values 
                            chanlabels={obj.info.chanlocs.labels}; % just the cell array of labels

                            chans2add = chanlabels(pvals<0.05); % significant channels to add to list

                            % Check if the field already exists
                            if isfield(obj.scalpResults, 'sigElectrodes') && isfield(obj.scalpResults.sigElectrodes, property) && ...
                                isfield(obj.scalpResults.sigElectrodes.(property), timeNames{t}) && isfield(obj.scalpResults.sigElectrodes.(property).(timeNames{t}), obj.info.freq_list{f})

                                % Get the existing cell array
                                existingArray = obj.scalpResults.sigElectrodes.(property).(timeNames{t}).(obj.info.freq_list{f});
                                % Check if the new letter is different from the existing ones
                                chans2add = setdiff(chans2add, existingArray);
                                obj.scalpResults.sigElectrodes.(property).(timeNames{t}).(obj.info.freq_list{f}) = [existingArray, chans2add];
                            elseif ~isempty(chans2add)
                                % Create a new cell array with the new letter
                                obj.scalpResults.sigElectrodes.(property).(timeNames{t}).(obj.info.freq_list{f}) = chans2add;
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
        
        function plotERPs(obj,vars2plot,freq2plot,times2plot,errorType)
            if nargin <5 || isempty(errorType)
                errorType='none';
            end
            
            if ~ismember(errorType, {'none', 'sem', '95CI'})
                error('Invalid value for errorType. Allowed values are ''none'', ''sem'', and ''95CI''.');
            end
            
            color_list={'g','b','r'};
            N=numel(fieldnames(obj.DATA));
            timeNames = fieldnames(obj.info.timeRange);
            sigElectrodes=obj.scalpResults.sigElectrodes;
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

                                %subplot(length(freq_list),maxLength,(f-1)*maxLength+c)
                                subplot(maxLength,length(freq_list),(c-1) * length(freq_list) + f)

                                freq=freq_list{f};
                                chan=sig.(freq){c};
                                freqIdx = find(strcmp(obj.info.freq_list, freq));
                                timeIdx = obj.info.timeIDX.(time);
                                elecIdxs = find(strcmp({obj.info.chanlocs.labels},chan ));
                                data=[];
                                hold on
                                h = gobjects(N, 1); % Preallocate an array for the line handles
                                for n=1:N
                                    d=obj.DATA.(obj.info.groupNames{n}).(property)(freqIdx,elecIdxs,timeIdx(1):timeIdx(2),:);
                                    data(:,n) = squeeze(nanmean(d,4));
                                    CI=[];
                                    sem=std(squeeze(d),0,2)/sqrt(size(d,4));
                                    
                                    ts = tinv([0.025  0.975],size(d,4)-1);   % T-Score
                                    CI(:,1) = ts(2)*sem; %upperCI
                                    CI(:,2) = ts(1)*sem; %lowerCI
                                    
                                    
                                    if strcmp(errorType,'sem')
                                        eb = shadedErrorBar(obj.info.timeAxis(timeIdx(1):timeIdx(2)),data(:,n),[sem,-sem], color_list{n},1);
                                    elseif strcmp(errorType,'95CI')
                                        eb = shadedErrorBar(obj.info.timeAxis(timeIdx(1):timeIdx(2)),data(:,n)',CI', color_list{n},1);
                                    end
                                    
                                    if strcmp(errorType,'sem') || strcmp(errorType,'95CI')
                                        h(n) = eb.mainLine;
                                        set(get(get(eb.patch, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
                                    end
                                    
                                end
                                hold off
                                if strcmp(errorType,'none')
                                    h=plot(obj.info.timeAxis(timeIdx(1):timeIdx(2)),data);
                                    arrayfun(@(x) set(h(x), 'Color', color_list{x}), 1:N);
                                end

                                ylabel(chan)
                                if c==1
                                    title(freq,'fontweight','bold','fontsize',16)
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
        function plotScalpBar(obj,vars2plot,freq2plot,times2plot)
            N=numel(fieldnames(obj.DATA));
            combinations = nchoosek(1:N,2);
            timeNames = fieldnames(obj.info.timeRange);
            sigElectrodes=obj.scalpResults.sigElectrodes;
            for p=1:length(vars2plot)
                property=vars2plot{p};
                time_list = intersect(fieldnames(sigElectrodes.(property)),times2plot);

                for t=1:length(time_list)

                    timeName=time_list{t};
                    sig=sigElectrodes.(property).(timeName);
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

                                %subplot(length(freq_list),maxLength,(f-1)*maxLength+c)
                                subplot(maxLength,length(freq_list),(c-1) * length(freq_list) + f)

                                freq=freq_list{f};
                                chan=sig.(freq){c};
                                freqIdx = find(strcmp(obj.info.freq_list, freq));
                                timeIdx = obj.info.timeIDX.(timeName);
                                elecIdxs = find(strcmp({obj.info.chanlocs.labels},chan ));
                                data=[];sem=[];
                                hold on
                                for n=1:N
                                    subdata = squeeze(nanmean(obj.getGroupData(obj.info.groupNames{n},property,freq,timeName,elecIdxs),3));
                                    data(n) = nanmean(subdata);
                                    sem(n) = std(subdata)/sqrt(length(subdata));
                                end
                                obj.plotErrBar(data,sem)
                                [ngroups, nbars] = size(data);
                                groupwidth = min(0.8, nbars/(nbars + 1.5));

                                for comb = 1:size(combinations, 1)
                                    group1=obj.info.groupNames{combinations(comb, 1)};
                                    group2=obj.info.groupNames{combinations(comb, 2)};
                                    s1=[];s2=[];
                                    s1 = squeeze(nanmean(obj.getGroupData(group1,property,freq,timeName,elecIdxs),3));
                                    s2 = squeeze(nanmean(obj.getGroupData(group2,property,freq,timeName,elecIdxs),3));
                                    pvals(comb,:)=obj.calGroupSig(s1,s2);
                                    group1Location = (1:ngroups) - groupwidth/2 + (2*combinations(comb, 1)-1) * groupwidth / (2*nbars);
                                    group2Location = (1:ngroups) - groupwidth/2 + (2*combinations(comb, 2)-1) * groupwidth / (2*nbars);
                                    A=[group1Location;group2Location]';
                                    groupingKey = mat2cell(A, ones(1, size(A, 1)), size(A, 2));
                                    sigstar(groupingKey,pvals(comb,:))
                                end
                                hold off
                                xlabel(chan)
                                xticks([])
                                if c==1
                                    title(freq,'fontweight','bold','fontsize',16)
                                end
                            end
                        end
                        Lgnd = legend(obj.info.groupNames);
                        Lgnd.Position(1) = 0.7;
                        Lgnd.Position(2) = 0.9;
                        sgtitle(strcat(property,'-',timeName))
                    end
                end
            end
        end
        %--------------END OF CLASS METHODS-----------------------%
    end
    %--------------END OF CLASS-----------------------%
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

