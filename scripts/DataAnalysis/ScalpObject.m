classdef ScalpObject < DataAnalysis
    properties
        DATA
        scalpResults
    end
    
    methods
        function obj = ScalpObject(DATA, info, baselineTime, timeRange,cfg)
            % ScalpObject  Constructor for scalp-level EEG analysis object
            %
            %   obj = ScalpObject(DATA, info, baselineTime, timeRange, cfg)
            %
            %   Initializes a ScalpObject instance, inheriting from the DataAnalysis
            %   superclass. This object is designed to process scalp-level EEG data 
            %   (e.g., channel-space time series) for group-level comparisons and 
            %   statistical analysis.
            %
            %   Inputs:
            %     DATA         - Struct containing scalp EEG data for each group
            %     info         - Struct with experiment metadata (group names, variables, frequencies, etc.)
            %     baselineTime - 1x2 vector specifying baseline window [start, end] in ms
            %     timeRange    - Struct defining named time windows (e.g., time1 = [0 500])
            %     cfg          - (Optional) Struct of user-defined analysis settings
            %
            %   Output:
            %     obj          - Instantiated ScalpObject with all properties initialized
            %
            %   Example:
            %     obj = ScalpObject(project.scalpData, project.info, [-250 -50], timeRange, cfg);
            %
            if nargin<5 || isempty(cfg)
                cfg=struct();
            end
            % Call superclass constructor
            obj = obj@DataAnalysis(info, baselineTime, timeRange,cfg);
            obj.DATA=DATA;
        end
        
        function obj=calChanData(obj)
            % calChanData  Compute scalp-level channel data across time windows
            %
            %   obj = calChanData(obj) averages the scalp EEG data across the time
            %   dimension for each group, frequency band, and time window specified in
            %   the object's metadata. The resulting channel-level data is stored in 
            %   `obj.scalpResults.chanData`.
            %
            %   Inputs:
            %     obj - ScalpObject instance containing:
            %              • obj.DATA: EEG data organized by group
            %              • obj.info: struct with fields:
            %                  - groupNames: group labels
            %                  - variables: names of EEG variables (e.g., 'ERP', 'Power')
            %                  - freq_list: frequency bands
            %                  - timeIDX: named time window indices
            %
            %   Outputs:
            %     obj - The same object with the following field populated:
            %              • obj.scalpResults.chanData.(group).(variable).(time).(freq)
            %                where data is averaged across time.
            %
            %   Example:
            %     scalpObj = scalpObj.calChanData();
            disp('-----------Calulating channel data and saving in obj.scalpResults.chanData-----------')
            properties=obj.info.variables;
            timeNames = fieldnames(obj.info.timeIDX);
            N = numel(fieldnames(obj.DATA));

            for p=1:length(properties)
                property=properties{p}; % get property name
                for n = 1:N
                    group=obj.info.groupNames{n};
                    for t=1:length(timeNames)
                        timeName=timeNames{t};
                        for f=1:length(obj.info.freq_list)
                            freq=obj.info.freq_list{f};
                            data = squeeze(nanmean(obj.getGroupData(group,property,freq,timeName),3));
                            obj.scalpResults.chanData.(obj.info.groupNames{n}).(property).(timeNames{t}).(obj.info.freq_list{f}) = data;
                        end
                    end
                end
            end
        end
        
        function obj=analyzeScalp(obj,varargin)
            
            %   [...] = analyzeScalp(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies additional
            %   parameters and their values. Will only analyze conditions that satisfy 
            %   all conditions specified. Default run will do all
            %   combinations. Option to plot scalp topo plots
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
            %         'FDRflag'      1/0 flag for fdr corrected p-values 
            %                        within each scalp map (# electrodes)
            %                        Default is 1 (fdr corrected values).
            %         'toPlot'       1/0 flag for plotting scalp topo
            %                        plots. Default is 1
            %         'isnormal'    'auto' 1 0 |
            %                        indicates distribution of the data for
            %                        significance testing. normal (1 flag) will use
            %                        ttest/ttest2. not normal (0 flag) will
            %                        use signrank/ranksum test. 'auto' will
            %                        determine if the data is normal with
            %                        adtest. 'auto' is default.
            %   Outputs:
            %     obj - The same object with the following field populated:
            %              • obj.scalpResults.sigElectrodes
            %              • obj.scalpResults.sigElectrodesP

            
            cfg = parseCfgOrArgs(obj, varargin{:});

            % Use fields directly
            vars2plot  = cfg.vars2plot;
            freq2plot  = cfg.freq2plot;
            times2plot  = cfg.times2plot;
            combinations = cfg.combinations;
            FDRflag    = cfg.FDRflag;
            toPlot     = cfg.toPlot;
            isnormal   = cfg.isnormal;
            disp('-----------Calculating group differences and saving in obj.scalpObject.scalpResults------------')


            % Initialize results
            if ~isfield(obj.scalpResults,'sigElectrodes')
                obj.scalpResults.sigElectrodes=struct();
                obj.scalpResults.sigElectrodesP=struct();
            end

            % normality test
            if strcmp(isnormal,'auto')
                % Get all Data
                all_data=[];
                for p=1:length(vars2plot)
                    property=vars2plot{p};
                    
                    for t=1:length(times2plot)
                        timeName=times2plot{t};
                        for f = 1:length(freq2plot)
                            freq=freq2plot{f};
                            for comb = unique(combinations)
                                groupName=obj.info.groupNames{comb};
                                tempData = obj.scalpResults.chanData.(groupName).(property).(timeName).(freq);
                                all_data = [all_data reshape(tempData,1,[])];
                            end
                        end
                    end
                end
                isnormal = ~adtest(all_data);
            end
            obj.scalpResults.isnormal=isnormal;
            for p=1:length(vars2plot)
                property=vars2plot{p}; % get property name
                for comb = 1:size(combinations, 1)

                    group1=obj.info.groupNames{combinations(comb, 1)};
                    group2=obj.info.groupNames{combinations(comb, 2)};
                    if toPlot
                        figure
                        figCount=0;
                    end
                    for t=1:length(times2plot)
                        timeName=times2plot{t};
                        
                        for f=1:length(freq2plot)
                            freq=freq2plot{f};
                            s1 = squeeze(nanmean(obj.getGroupData(group1,property,freq,timeName),3));
                            s2 = squeeze(nanmean(obj.getGroupData(group2,property,freq,timeName),3));
                            [pvals]=obj.calGroupSig(s1,s2,obj.info.experimentalDesign,isnormal); % get 1x n channel of p values
                            if FDRflag
                                pvals = fdr(pvals);
                            end

                            obj = addScalpSig(obj,property, timeName, append(group2,"_",group1), freq, pvals); 
                            if toPlot
                                figCount = figCount + 1;
                                subplot(length(times2plot),length(freq2plot),figCount)
                                plotSigTopo(s1,s2,obj.info.chanlocs,pvals,{'.','+'})
                                if t==1
                                    title(freq2plot{f})
                                end
                                if f==1
                                    text(-1,-.2, times2plot{t},'fontweight','bold','fontsize',16,'rotation',90);         
                                end
                            end
                        end
                    end
                    if toPlot
                        sgtitle(property+" "+group2+"-"+group1)
                    end
                end
            end
        end
        
        function plotERPs(obj,varargin)
            % plotERPs  Plot overlayed ERPs across groups for selected channels, frequencies, and time windows
            %
            %   plotERPs(obj, 'PARAM1', val1, 'PARAM2', val2, ...)
            %
            %   This method visualizes averaged ERP waveforms for different groups
            %   using the scalp-level data in the object. It supports plotting across
            %   multiple frequencies, channels, and time windows with optional error
            %   bars.
            %
            %   Parameters (Name-Value pairs) Default values are in cfg:
            %      
            %     'vars2plot'   - Cell array of variable names to plot 
            %     'freq2plot'   - Cell array of frequency bands to plot 
            %     'times2plot'  - Cell array of time window names to plot 
            %     'chans2plot'  - Cell array of channel labels to plot 
            %     'groups2plot  - Numeric array of groups to plot
            %     'errorType'   - Type of error bar to include:
            %                       • 'none'  - No error bars (default)
            %                       • 'sem'   - Standard error of the mean
            %                       • '95CI'  - 95% confidence intervals
            %     'color_list'  - Cell array of color codes for plotting each group
            %                     
            %
            %   Behavior:
            %     - If 'chans2plot' is set to 'all', only significant electrodes (from obj.scalpResults.sigElectrodes)
            %       are included for each frequency/time combination.
            %     - Uses `shadedErrorBar` for optional error visualizations.
            %
            %   Example:
            %     plotERPs(obj, 'vars2plot', {'ERP'}, 'freq2plot', {'theta'}, ...
            %                  'times2plot', {'time1'}, 'chans2plot', {'Fz','Cz'}, ...
            %                  'errorType', 'sem', 'color_list', {'r','b','g'})
            %

            cfg = parseCfgOrArgs(obj, varargin{:});

            % Use fields directly
            vars2plot  = cfg.vars2plot;
            freq2plot  = cfg.freq2plot;
            times2plot  = cfg.times2plot;
            chans2plot = cfg.chans2plot;
            errorType = cfg.errorType;
            color_list = cfg.color_list;
            groups2plot = cfg.groups2plot;
            N=length(groups2plot);

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
                        figure
                        for f=1:length(freq_list)
                            if isequal(chans2plot,'all')
                                chan_list = sig.(freq_list{f});
                            end
                            for c=1:length(chan_list)

%                                 subplot(length(freq_list),maxLength,(f-1)*maxLength+c)
                                subplot(maxLength,length(freq_list),(c-1) * length(freq_list) + f)

                                freq=freq_list{f};
                                chan=chan_list{c};
                                timeIdx = obj.info.timeIDX.(time);
                                elecIdxs = find(strcmp({obj.info.chanlocs.labels},chan ));
                                data=[];
                                hold on
                                h = gobjects(N, 1); % Preallocate an array for the line handles
                                for n=1:N
                                    d=obj.getGroupData(obj.info.groupNames{n},property,freq,time,elecIdxs);
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
                        Lgnd = legend(h,obj.info.groupNames(groups2plot));
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
            %         'groups2plot'  vector array of groups in the object
            %                        to plot. Default is to plot all
            %                        groups. Can also be used to rearrange
            %                        bar order
            %         'color_list'   cell array to determine colors
            %                        default is {'r','b','g','m','k','c','y'}


            cfg = parseCfgOrArgs(obj, varargin{:});
            % Use fields directly
            vars2plot  = cfg.vars2plot;
            freq2plot  = cfg.freq2plot;
            times2plot  = cfg.times2plot;
            chans2plot = cfg.chans2plot;
            color_list = cfg.color_list;
            groups2plot= cfg.groups2plot;

            groups = obj.info.groupNames(groups2plot);
            combinations = nchoosek(groups2plot,2);
            N=length(groups);

            sigElectrodes=obj.scalpResults.sigElectrodes;
            for p=1:length(vars2plot)
                property=vars2plot{p};
                time_list = intersect(fieldnames(sigElectrodes.(property)),times2plot);
                
                for t=1:length(time_list)

                    timeName=time_list{t};
                    sig=sigElectrodes.(property).(timeName);
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

                    if ~isempty(freq_list)
                        figure
                        for f=1:length(freq_list)
                            
                            if isequal(chans2plot,'all')
                                chan_list = sig.(freq_list{f});
                            end

                            for c=1:length(chan_list)
                                
                                if length(freq_list)==1
                                    maxLength=length(chan_list);
                                end
                                
                                %subplot(length(freq_list),maxLength,(f-1)*maxLength+c)
                                subplot(maxLength,length(freq_list),(c-1) * length(freq_list) + f);

                                freq=freq_list{f};
                                chan=chan_list{c};
                                elecIdxs = find(strcmp({obj.info.chanlocs.labels},chan ));
                                data=[];sem=[];
                                hold on
                                % Get all the group data
                                for n=1:N
                                    subdata = squeeze(nanmean(obj.getGroupData(groups{n},property,freq,timeName,elecIdxs),3));
                                    data(n) = nanmean(subdata);
                                    sem(n) = std(subdata)/sqrt(length(subdata));
                                end
                                obj.plotErrBar(data,sem,color_list); % plot actual bars
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

        function obj = addScalpSig(obj,property, timeName, groupDiff, freq, pvals)
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
            obj.scalpResults.sigElectrodesP.(property).(timeName).(groupDiff).(freq)=pvals;
    
            chanlabels={obj.info.chanlocs.labels}; % just the cell array of labels
            chans2add = chanlabels(pvals<0.05); % significant channels to add to list
    
            % Check if the field already exists


            if isfield(obj.scalpResults, 'sigElectrodes') && isfield(obj.scalpResults.sigElectrodes, property) && ...
                isfield(obj.scalpResults.sigElectrodes.(property), timeName) && isfield(obj.scalpResults.sigElectrodes.(property).(timeName), freq)
                % Get the existing cell array
                existingArray = obj.scalpResults.sigElectrodes.(property).(timeName).(freq);
                % Check if the new letter is different from the existing ones
                chans2add = setdiff(chans2add, existingArray);
                obj.scalpResults.sigElectrodes.(property).(timeName).(freq) = [existingArray, chans2add];
            elseif ~isempty(chans2add)
                % Create a new cell array with the new letter
                obj.scalpResults.sigElectrodes.(property).(timeName).(freq) = chans2add;
              
            end

            if ~isfield(obj.scalpResults.sigElectrodes, property)
                obj.scalpResults.sigElectrodes.(property)=struct();
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

