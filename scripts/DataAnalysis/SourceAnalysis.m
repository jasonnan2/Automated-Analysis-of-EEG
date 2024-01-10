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
                        
                        obj.plotErrBar(data,sem)
                        [ngroups, nbars] = size(data);
                        groupwidth = min(0.8, nbars/(nbars + 1.5));

                        for comb = 1:size(combinations, 1)
                            group1=obj.info.groupNames{combinations(comb, 1)};
                            group2=obj.info.groupNames{combinations(comb, 2)};
                            s1=[];s2=[];
                            for i=1:length(netwrk)
                                % collapsing networks
                                s1(i,:) = squeeze(nanmean(nanmean(obj.getGroupData(group1,property,freq,timeName,netwrk(i).roi),2),3));
                                s2(i,:) = squeeze(nanmean(nanmean(obj.getGroupData(group2,property,freq,timeName,netwrk(i).roi),2),3));
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
        function obj=plotBrainmap(obj,properties,combinations)
            hm = headModel.loadFromFile('headModel_templateFile_newPEBplus.mat');
            T = hm.indices4Structure(hm.atlas.label);
            T = double(T)';
            timeNames = fieldnames(obj.info.timeIDX);
            N = numel(fieldnames(obj.DATA));
            
            if ~exist('combinations')
                combinations = nchoosek(1:N,2); % This will give you a matrix where each row is a combination of two groups
            end
            
            for p=1:length(properties)
                property=properties{p}; % get property name
                for comb = 1:size(combinations, 1)
                    group1=obj.info.groupNames{combinations(comb, 1)};
                    group2=obj.info.groupNames{combinations(comb, 2)};    
                    for t=1:length(timeNames)
                        timeName=timeNames{t};
                        plotdata=[];
                        for f=1:length(obj.info.freq_list)
                            freq=obj.info.freq_list{f};
                            % average in time dimension
                            s1 = squeeze(nanmean(obj.getGroupData(group1,property,freq,timeName),3));
                            s2 = squeeze(nanmean(obj.getGroupData(group2,property,freq,timeName),3));
                            pvals(:,f)=obj.calGroupSig(s1,s2); % get 1x n channel of p values 
                            plotdata(:,f)=nanmean(s2,2)-nanmean(s1,2); % average across subjects and subtract
                            roi2add = obj.info.roi(pvals(:,f)<=0.05);
                            %%% adding ROI's to results
                            % Check if the field already exists
                            if isfield(obj.sourceResults, 'sigROIs') && isfield(obj.sourceResults.sigROIs, property) && ...
                                isfield(obj.sourceResults.sigROIs.(property), timeNames{t}) && isfield(obj.sourceResults.sigROIs.(property).(timeNames{t}), obj.info.freq_list{f})
                                % Get the existing cell array
                                existingArray = obj.sourceResults.sigROIs.(property).(timeNames{t}).(obj.info.freq_list{f});
                                % Check if the new letter is different from the existing ones
                                roi2add = setdiff(roi2add, existingArray);
                                obj.sourceResults.sigROIs.(property).(timeNames{t}).(obj.info.freq_list{f}) = [existingArray; roi2add];
                            elseif ~isempty(roi2add)
                                % Create a new cell array with the new letter
                                obj.sourceResults.sigROIs.(property).(timeNames{t}).(obj.info.freq_list{f}) = roi2add;
                            end
                        end
                        
                        
                        %plotdata(pvals>0.05)=0; % Uncorrected p-value thresholding only for PostPre difference condition
                        
                        tohide = T'*(pvals>0.05);
                        X=T'*plotdata;
                        X(tohide==1)=0;
                        plot68roi(hm,X , 1,obj.info.freq_list)
                        hAx = axes('Position', [0, 0, 1, 1], 'Visible', 'off');
                        text(0.5, 1, property+" "+timeName+" "+group2+"-"+group1, 'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
                    end
                end
            end
        end
        %--------------END OF CLASS METHODS-----------------------%
    end
    %--------------END OF CLASS-----------------------%
end
