classdef ScalpObject < DataAnalysis
    properties
        DATA
        scalpResults
    end
    
    methods
        function obj = ScalpObject(DATA, info, baselineTime, timeRange)
            % Call superclass constructor
            obj = obj@DataAnalysis(info, baselineTime, timeRange);
            obj.DATA=DATA;
        end
        
        function obj=ScalpAnalysis(obj)
            % obj=ScalpAnalysis(obj)
            % runs all comparisons between groups in scalp space

            % Outputs:
            %
            %   Results are saved into the following fields of the input object `obj`:
            %
            %   - obj.scalpResults.sigElectrodesP.(property).(timeName).(group1_group2).(freq)
            %       Raw p-values for group comparisons at each electrode.
            %       These are not FDR-corrected.
            %
            %   - obj.scalpResults.sigElectrodes.(property).(timeName).(freq)
            %       Significance results across electrodes 
            %       These are not FDR-corrected


            properties=obj.info.variables;
            timeNames = fieldnames(obj.info.timeIDX);
            N = numel(fieldnames(obj.DATA));
            combinations = nchoosek(1:N,2); % This will give you a matrix where each row is a combination of two groups

            for p=1:length(properties)
                property=properties{p}; % get property name
                for comb = 1:size(combinations, 1)
                    group1=obj.info.groupNames{combinations(comb, 1)};
                    group2=obj.info.groupNames{combinations(comb, 2)};
                    for t=1:length(timeNames)
                        timeName=timeNames{t};
                        for f=1:length(obj.info.freq_list)
                            freq=obj.info.freq_list{f};

                            s1 = squeeze(nanmean(obj.getGroupData(group1,property,freq,timeName),3));
                            s2 = squeeze(nanmean(obj.getGroupData(group2,property,freq,timeName),3));
                            [pvals, stats] =obj.calGroupSig(s1,s2,obj.info.experimentalDesign); % get 1x n channel of p values
                            obj.scalpResults.sigElectrodesP.(property).(timeNames{t}).(append(group1,"_",group2)).(obj.info.freq_list{f})=pvals;
                            %obj.scalpResults.sigElectrodesStats.(property).(timeNames{t}).(append(group1,"_",group2)).(obj.info.freq_list{f})=stats;

                            chanlabels={obj.info.chanlocs.labels}; % just the cell array of labels
%                             fdrP = fdr(pvals);
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
                        end
                    end
                end
            end
        end
        
        function plotScalpMap(obj,varargin)
            
            %   [...] = plotERPs(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies additional
            %   parameters and their values. Will only plot conditions that satisfy 
            %   all conditions specified. Default run will plot all combinations
            %   Valid parameters are the following:
            %
            %        Parameter  Value
            %         'vars2plot'    cell array of variables in the object
            %                        to plot. Default is to plot all
            %                        variables
            %                  
            %         'freq2plot'    cell array of frequencies in the object
            %                        to plot. Default is to plot all
            %                        frequencies
            %                   
            %         'times2plot'   cell array of time ranges in the object
            %                        to plot. Default is to plot all
            %                        time ranges
            %         'combinations' 2D matrix defining which groups to run
            %                        comparisons on.
            %                        ex: [1,2] will do groups 1 vs 2
            %                            [1,2;1,3] is 1 vs 2 and 1 vs 3
            %         'FDRflag'      1/0 flag for plotting fdr corrected p-values correcting
            %                        within each scalp map (# electrodes)
            %                        Default is 1 (plot corrected values).
            
            N = length(obj.info.groupNames);
            pnames = {'vars2plot','freq2plot','times2plot','combinations','FDRflag'};
            dflts  = {obj.info.variables,obj.info.freq_list, fieldnames(obj.info.timeRange),nchoosek(1:N,2),1 };
            [vars2plot,freq_list,timeNames,combinations, FDRflag] = parseArgs(pnames,dflts,varargin{:});
            

            for p=1:length(vars2plot)
                property=vars2plot{p}; % get property name
                for comb = 1:size(combinations, 1)

                    group1=obj.info.groupNames{combinations(comb, 1)};
                    group2=obj.info.groupNames{combinations(comb, 2)};

                    figure
                    figCount=0;
                    for t=1:length(timeNames)
                        timeName=timeNames{t};
                        
                        for f=1:length(freq_list)
                            
                            freq=freq_list{f};
                            figCount = figCount + 1;
                            subplot(length(timeNames),length(freq_list),figCount)
                            s1 = squeeze(nanmean(obj.getGroupData(group1,property,freq,timeName),3));
                            s2 = squeeze(nanmean(obj.getGroupData(group2,property,freq,timeName),3));
                            pvals=obj.calGroupSig(s1,s2,obj.info.experimentalDesign); % get 1x n channel of p values
                            if FDRflag
                                pvals = fdr(pvals);
                            end
                            plotSigTopo(s1,s2,obj.info.chanlocs,pvals,{'.','+'})
                            if t==1
                                title(freq_list{f})
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
        
        function plotERPs(obj,varargin)
            %   [...] = plotERPs(...,'PARAM1',VAL1,'PARAM2',VAL2,...)
            %   Plots overlayed ERPs between groups. Non-FDR corrected mask
            %   Specifies additional parameters and their values.
            %   Valid parameters are the following:
            %
            %        Parameter  Value
            %         'vars2plot'    cell array of variables in the object
            %                        to plot. Default is to plot all
            %                        variables
            %                  
            %         'freq2plot'    cell array of frequencies in the object
            %                        to plot. Default is to plot all
            %                        frequencies
            %                   
            %         'times2plot'   cell array of time ranges in the object
            %                        to plot. Default is to plot all
            %                        time ranges
            %         'chans2plot'   cell array of channels to plot.
            %                        Default is all significant ones
            %         'FDRflag'      1/0 flag for plotting fdr corrected p-values correcting
            %                        within each scalp map (# electrodes)
            %                        Default is 0 (plot uncorrected values).
            %
            %         'errorType'    string value of the type of error to
            %                        plot. choices are 'none','sem', and
            %                        '95CI'. Default is 'none'
            %         'color_list'   cell array to determine colors
            %                        default is {'g','b','r'}


            N=length(obj.info.groupNames);
            pnames = {'vars2plot','freq2plot','times2plot','chans2plot','FDRflag','errorType','color_list'};
            dflts  = {obj.info.variables,obj.info.freq_list, fieldnames(obj.info.timeRange),'all',0,'none',{'g','b','r'}};
            [vars2plot,freq2plot,times2plot,chans2plot, FDRflag,errorType,color_list] = parseArgs(pnames,dflts,varargin{:});
            if ~ismember(errorType, {'none', 'sem', '95CI'})
                error('Invalid value for errorType. Allowed values are ''none'', ''sem'', and ''95CI''.');
            end
           
            sigElectrodes=obj.scalpResults.sigElectrodes;
            for p=1:length(vars2plot)
                property=vars2plot{p};
                time_list = intersect(fieldnames(sigElectrodes.(property)),times2plot);

                for t=1:length(time_list)
                    
                    time=time_list{t};
                    sig=sigElectrodes.(property).(time);
                    [~,~,ib] = intersect(fieldnames(sig),freq2plot);
                    freq_list=freq2plot(sort(ib));

                    if isequal(chans2plot,'all')
                        
                        % Get max length of all frequencies in one view
                        maxLength = 0;
                        for i = 1:length(freq_list)
                            currentLength = length(sig.(freq_list{i}));  % Get the length of the current cell array
                            if currentLength > maxLength
                                maxLength = currentLength;  % Update the maximum length
                            end
                        end
                    else
                        chan_list = chans2plot;
                        maxLength=length(chan_list);
                    end

                    % get longest electrode map
                    

                    if ~isempty(freq_list)
%                         figure
                        for f=1:length(freq_list)
                            chan_list = sig.(freq_list{f});
                            for c=1:length(chan_list)

                                %subplot(length(freq_list),maxLength,(f-1)*maxLength+c)
                                subplot(maxLength,length(freq_list),(c-1) * length(freq_list) + f)

                                freq=freq_list{f};
                                chan=chan_list{c};
                                freqIdx = find(strcmp(obj.info.freq_list, freq));
                                timeIdx = obj.info.timeIDX.(time);
                                elecIdxs = find(strcmp({obj.info.chanlocs.labels},chan ));
                                data=[];
                                hold on
                                h = gobjects(N, 1); % Preallocate an array for the line handles
                                for n=1:N

                                    obj.getGroupData(group1,property,freq,timeName)
                                    
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
        function plotScalpBar(obj,varargin)
            %   [...] = plotScalpBar(...,'PARAM1',VAL1,'PARAM2',VAL2,...) 
            %   Plots grouped bar plots which are significantly different. 
            %   Non-FDR corrected pvalues. 
            %   specifies additional parameters and their values.  
            %   Valid parameters are the following:
            %
            %        Parameter  Value
            %         'vars2plot'    cell array of variables in the object
            %                        to plot. Default is to plot all
            %                        variables
            %                  
            %         'freq2plot'    cell array of frequencies in the object
            %                        to plot. Default is to plot all
            %                        frequencies
            %                   
            %         'times2plot'   cell array of time ranges in the object
            %                        to plot. Default is to plot all
            %                        time ranges
            %         'chans2plot'   cell array of channels in the object
            %                        to plot. Default is to plot
            %                        significant channels
            %         'groups'       cell array of groups in the object
            %                        to plot. Default is to plot all
            %                        groups. Can also be used to rearrange
            %                        bar order
            pnames = {'vars2plot','freq2plot','times2plot','chans2plot','groups'};
            dflts  = {obj.info.variables,obj.info.freq_list, fieldnames(obj.info.timeRange),{obj.info.chanlocs.labels},fieldnames(obj.DATA)};
            [vars2plot,freq2plot,times2plot,chans2plot,groups] = parseArgs(pnames,dflts,varargin{:});
            N=numel(groups);
            combinations = nchoosek(1:N,2);

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
                            chan_list = intersect(chans2plot,sig.(freq_list{f}));
%                             if length(chan_list)<length(chans2plot)
%                                 chan_list=chans2plot;
%                             end

                            for c=1:length(chan_list)
                                
                                if length(freq_list)==1
                                    maxLength=length(chan_list);
                                end
                                
                                %subplot(length(freq_list),maxLength,(f-1)*maxLength+c)
                                subplot(maxLength,length(freq_list),(c-1) * length(freq_list) + f);

                                freq=freq_list{f};
                                chan=chan_list{c};
                                freqIdx = find(strcmp(obj.info.freq_list, freq));
                                timeIdx = obj.info.timeIDX.(timeName);
                                elecIdxs = find(strcmp({obj.info.chanlocs.labels},chan ));
                                data=[];sem=[];
                                hold on
                                % Get all the group data
                                for n=1:N
                                    subdata = squeeze(nanmean(obj.getGroupData(groups{n},property,freq,timeName,elecIdxs),3));
                                    data(n) = nanmean(subdata);
                                    sem(n) = std(subdata)/sqrt(length(subdata));
                                end
                                obj.plotErrBar(data,sem); % plot actual bars
                                ylim([-0.025 0.045])
                                % getting significance bars
                                [ngroups, nbars] = size(data);
                                groupwidth = min(0.8, nbars/(nbars + 1.5));
                                for comb = 1:size(combinations, 1)
                                    group1=groups{combinations(comb, 1)};
                                    group2=groups{combinations(comb, 2)};
                                    s1=[];s2=[];
                                    s1 = squeeze(nanmean(obj.getGroupData(group1,property,freq,timeName,elecIdxs),3));
                                    s2 = squeeze(nanmean(obj.getGroupData(group2,property,freq,timeName,elecIdxs),3));
                                    pvals(comb,:)=obj.calGroupSig(s1,s2,obj.info.experimentalDesign);
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
                        Lgnd = legend(groups,'fontweight','bold','fontsize',9);
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
function plotSigTopo(s1,s2,chanlocs,pvals, key)
% plotSigTopo  Plot a topographic map of EEG channel differences with significant electrodes
%
%   plotSigTopo(s1, s2, chanlocs, pvals, key) computes and plots the
%   difference in EEG signal between two conditions (s2 - s1) across channels.
%   Electrodes with p < 0.05 are labeled using the provided key.
%
%   Inputs:
%       s1       - (channels x subjects) data matrix for condition 1
%       s2       - (channels x subjects) data matrix for condition 2
%       chanlocs - EEGLAB-style structure of channel locations
%       pvals    - Vector of p-values for each channel
%       key      - Cell array of markers to indicate sig vs no sig. ie:
%                  {'*','.')
%
%   This function:
%       - Computes average difference between s2 and s1
%       - Masks and labels significant electrodes (p < 0.05)
%       - Plots scalp topography using EEGLAB's topoplot
%
%   Example:
%       plotSigTopo(data1, data2, EEG.chanlocs, pvals, EEG.chanlabels)


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

