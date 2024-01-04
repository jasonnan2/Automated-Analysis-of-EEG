classdef SourceAnalysis < DataAnalysis
    properties
        DATA
        sourceResults
    end
    methods
        function obj = SourceAnalysis(DATA, info, baselineTime, timeRange)
            % Call superclass constructor
            obj = obj@DataAnalysis( info, baselineTime, timeRange);
            obj.DATA=DATA;
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
                            pvals(comb,:)=obj.calGroupSig(s1,s2);

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
    end
end
