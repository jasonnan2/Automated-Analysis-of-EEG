classdef SourceObject < DataAnalysis
    properties
        DATA
        sourceResults
        netConnectivity
    end
    methods
        function obj = SourceObject(DATA, info, baselineTime, timeRange,cfg)
            % SourceObject  Constructor for source-level analysis object
            %
            %   obj = SourceObject(DATA, info, baselineTime, timeRange, cfg)
            %
            %   Initializes a SourceObject instance, which inherits from the
            %   DataAnalysis superclass. This object manages source-level neural
            %   data and configuration for ROI/network-based processing.
            %
            %   Inputs:
            %     DATA         - Struct containing source-level data for each group
            %     info         - Struct with metadata (group names, frequencies, variables, etc.)
            %     baselineTime - 1x2 vector specifying baseline window [start, end] in ms
            %     timeRange    - Struct defining named time windows (e.g., time1 = [0 500])
            %     cfg          - (Optional) Struct of configuration options (e.g., frequencies, plots)
            %
            %   Output:
            %     obj          - Initialized SourceObject with inherited DataAnalysis properties
            %
            %   Example:
            %     obj = SourceObject(project.sourceData, project.info, [-250 -50], timeRange, cfg);

            if nargin<5 || isempty(cfg)
                cfg=struct();
            end
            % Call superclass constructor
            obj = obj@DataAnalysis( info, baselineTime, timeRange,cfg);
            obj.DATA=DATA;
        end
        
        function plotNetwork(obj,varargin)
            %   [...] = plotNetwork(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies additional
            %   parameters and their values. Will only plot conditions that satisfy 
            %   all conditions specified. Default run will plot all
            %   combinations. Plots shown are masked by fdr correction
            %   across brain region. 
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
            %         'groups2plot'  array defining which groups show
            %                        comparisons on.
            %                        ex: [1,2] will do groups 1 vs 2
            %                            [1,2,3] will show gruops 1,2,3
            %         'FDRflag'      1/0 flag for plotting fdr corrected p-values correcting
            %                        within each scalp map (# networks)
            %                        Default is 1 (plot corrected values).
            %         'isnormal'    'auto' 1 0 |
            %                        indicates distribution of the data for
            %                        significance testing. normal (1 flag) will use
            %                        ttest/ttest2. not normal (0 flag) will
            %                        use signrank/ranksum test. 'auto' will
            %                        determine if the data is normal with
            %                        adtest. 'auto' is default.
            %
            %         'color_list'   cell array to determine colors
            %                        default is {'r','b','g','m','k','c','y'}
            
            cfg = parseCfgOrArgs(obj, varargin{:});
            % Use fields directly
            vars2plot  = cfg.vars2plot;
            freq2plot  = cfg.freq2plot;
            times2plot  = cfg.times2plot;
            groups2plot= cfg.groups2plot;
            FDRflag = cfg.FDRflag;
            isnormal=cfg.isnormal;
            color_list=cfg.color_list;

            netwrk=obj.info.netwrk;
            combinations = nchoosek(groups2plot,2);


            % normality test
            if strcmp(isnormal,'auto')
                % Get all Data
                all_data=[];
                for p=1:length(vars2plot)
                    property=vars2plot{p};
                    for t=1:length(times2plot)
                        timeName=times2plot{t};
                        for f = 1:length(obj.info.freq_list)
                            freq=obj.info.freq_list{f};
                            for comb = unique(combinations)
                                groupName=obj.info.groupNames{comb};
                                tempData = obj.sourceResults.netData.(groupName).(property).(timeName).(freq);
                                all_data = [all_data reshape(tempData,1,[])];
                            end
                        end
                    end
                end
                isnormal = ~adtest(all_data);
            end
            obj.sourceResults.isnormalNet=isnormal;

            for p=1:length(vars2plot)
                property=vars2plot{p};
                figure
                count=0;
                for t=1:length(times2plot)
                    timeName=times2plot{t};
                    for f = 1:length(freq2plot)
                        freq=freq2plot{f};
                        % Generating plot data
                        data=[];sem=[];
                        
                        for comb = 1:length(groups2plot)
                            group=obj.info.groupNames{groups2plot(comb)};
                            tempData = obj.sourceResults.netData.(group).(property).(timeName).(freq); % size net x groups
                            data(:,comb) = nanmean(tempData,2);
                            sem(:,comb) = std(tempData,[],2)/sqrt(size(tempData,2));
                        end
                        
                        count=count+1;
                        subplot(length(times2plot),length(freq2plot),count)
                        
                        obj.plotErrBar(data,sem,color_list)
                        [ngroups, nbars] = size(data);
                        groupwidth = min(0.8, nbars/(nbars + 1.5));
                        
                        combinations = nchoosek(groups2plot,2);
                        for comb = 1:size(combinations, 1)
                            group1=obj.info.groupNames{combinations(comb, 1)};
                            group2=obj.info.groupNames{combinations(comb, 2)};
                            s1=obj.sourceResults.netData.(group1).(property).(timeName).(freq);
                            s2=obj.sourceResults.netData.(group2).(property).(timeName).(freq);

                            pvals(comb,:)=obj.calGroupSig(s1,s2,obj.info.experimentalDesign,isnormal);
                            if FDRflag
                                pvals(comb,:)=fdr(pvals(comb,:));
                            end

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
                           ylabel(times2plot{t},'fontweight','bold','fontsize',16,'rotation',90);         
                        end
                     end
                end
                Lgnd = legend(obj.info.groupNames(groups2plot));
                Lgnd.Position(1) = 0.7;
                Lgnd.Position(2) = 0.9;
                sgtitle(property)
            end
        end

        function obj=calRoiData(obj)
            % calRoiData  Computes ROI-averaged source activity for each group, frequency, and time window
            %
            %   obj = calRoiData(obj) averages the source-level neural data across the
            %   time dimension within each defined time window and stores the resulting
            %   ROI-level data in the `obj.sourceResults.roiData` field. The function
            %   processes all groups, frequency bands, and time windows defined in the
            %   `obj.info` structure.
            %
            %   Input:
            %     obj  - SourceObject instance from running
            %            SourceObject(project.sourceData, project.info, baselineTime, timeRange, cfg);
            %             
            %
            %   Output:
            %     obj  - The input object with the following field populated:
            %              • obj.sourceResults.roiData: ROI-level averages indexed by
            %                group, variable, time window, and frequency.
            %
            %   Example usage:
            %     sourceObject = sourceObject.calRoiData();

            disp('-----------Calulating ROI data and saving in obj.sourceResults.roiData-----------')

            timeNames = fieldnames(obj.info.timeIDX);
            N = numel(fieldnames(obj.DATA));
            combinations = nchoosek(1:N,2); % This will give you a matrix where each row is a combination of two groups
            properties=obj.info.variables;
            for p=1:length(properties)
                property=properties{p}; % get property name
                for n = 1:N
                    group=obj.info.groupNames{n};
                    for t=1:length(timeNames)
                        timeName=timeNames{t};
                        for f=1:length(obj.info.freq_list)
                            freq=obj.info.freq_list{f};
                            % average in time dimension
                            data = squeeze(nanmean(obj.getGroupData(group,property,freq,timeName),3));
                            obj.sourceResults.roiData.(obj.info.groupNames{n}).(property).(timeName).(freq) = data;
                        end
                    end
                end
            end
        end
        function obj = calNetData(obj, netwrk)
        % calNetData Calculate network-level data from ROI-level source results
        %
        %   obj = calNetData(obj, netwrk)
        %
        %   This method processes ROI-level source data and aggregates it based on
        %   predefined brain network definitions provided in the `netwrk` structure.
        %
        %   Inputs:
        %     obj     - Class object containing source-level analysis results
        %     netwrk  - Struct array defining brain networks. Each element must contain:
        %                 .name : Name of the network (e.g., 'FPN', 'DMN')
        %                 .roi  : Vector of ROI indices that define the network
        %   Outputs:
        %     obj - Updated object containing network-level averaged data
        %
        %     Results are saved into:
        %
        %     - obj.sourceResults.netData.(group).(property).(timeName).(freq)
        %         A matrix of size [Nnetworks x Nsubjects] containing network-averaged
        %         values for each subject, computed by averaging ROI-level data within
        %         each predefined network.
        % Validate netwrk input
        
            if nargin<1
                error('Require netwrk definition')
            end
            if ~isstruct(netwrk)
                error('Input "netwrk" must be a struct array.');
            end
            
            requiredFields = {'name', 'roi'};
            for i = 1:numel(netwrk)
                for f = 1:numel(requiredFields)
                    if ~isfield(netwrk(i), requiredFields{f})
                        error('Each element in "netwrk" must contain a "%s" field.', requiredFields{f});
                    end
                end
                if ~ischar(netwrk(i).name) && ~isstring(netwrk(i).name)
                    error('netwrk(%d).name must be a string or character vector.', i);
                end
                if ~isnumeric(netwrk(i).roi) || isempty(netwrk(i).roi)
                    error('netwrk(%d).roi must be a non-empty numeric vector.', i);
                end
            end
            obj.info.netwrk = netwrk;
            disp('-----------Calulating network data and saving in obj.sourceResults.netData-----------')

            timeNames = fieldnames(obj.info.timeIDX);
            N = numel(fieldnames(obj.DATA));
            properties=obj.info.variables;
            for p=1:length(properties)
                property=properties{p}; % get property name
                for n = 1:N
                    group=obj.info.groupNames{n};
                    for t=1:length(timeNames)
                        timeName=timeNames{t};
                        for f=1:length(obj.info.freq_list)
                            freq=obj.info.freq_list{f};
                            % average in time dimension
                            data = squeeze(nanmean(obj.getGroupData(group,property,freq,timeName),3)); % size roi x Nsubs
                            netData =[];
                            for net=1:length(netwrk)
                                netData(net,:) = nanmean(data(netwrk(net).roi,:),1);
                            end
                            obj.sourceResults.netData.(group).(property).(timeName).(freq) = netData;
                        end
                    end
                end
            end
        end

        function obj=analyzeRoi(obj,varargin)
            
            %   [...] = analyzeRoi(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies additional
            %   parameters and their values. Will only plot conditions that satisfy 
            %   all conditions specified. Default run will plot all
            %   combinations. Plots shown are masked by fdr correction
            %   across brain region. 
            %   Valid parameters are the following:
            %
            %        Parameter  Value
            %         'hmFile'       filepath to a specific headmodel file
            %                        to use in plotting. From
            %                        github.com/aojeda/dsi. If none is
            %                        specified, the default will be used. 
            %         'rois2plot'    numeric array of variables in the object
            %                        to plot. Default is to plot all
            %                        ROIs
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
            %              • obj.sourceResults.sigROIs
            %              • obj.scalpResults.sigROIsP
            
            cfg = parseCfgOrArgs(obj, varargin{:});
            % Use fields directly
            hmFile = cfg.hmFile;
            rois2plot = cfg.rois2plot;
            combinations = cfg.combinations;
            vars2plot  = cfg.vars2plot;
            freq2plot  = cfg.freq2plot;
            times2plot  = cfg.times2plot;
            FDRflag = cfg.FDRflag;
            toPlot = cfg.toPlot;
            isnormal=cfg.isnormal;

            % Initialize results structs
            if ~isfield(obj.sourceResults,'sigROIs')
                obj.sourceResults.sigROIs=struct();
                obj.sourceResults.sigROIsP=struct();
            end
            disp('-----------Calculating group differences and saving in obj.sourceResults------------')
            

            % check if headmodel and plotting functions exist
            if toPlot
                % Try to get default head model path
                try
                    defaultHMfile = headModel.getDefaultTemplateFilename();
                catch
                    error('Cannot find head model for plotting. Ensure github.com/aojeda/dsi is installed.');
                end
            
                % Validate head model path
                if ~exist(hmFile, 'file')
                    if ~isempty(defaultHMfile) && exist(defaultHMfile, 'file')
                        hmFile = defaultHMfile;
                        disp('Using default head model for plotting.');
                    else
                        error('Can not find any head model files.');
                    end
                end
            end



            % normality test
            if strcmp(isnormal,'auto')
                % Get all Data
                all_data=[];
                for p=1:length(vars2plot)
                    property=vars2plot{p};
                    for t=1:length(times2plot)
                        timeName=times2plot{t};
                        for f = 1:length(obj.info.freq_list)
                            freq=obj.info.freq_list{f};
                            for comb = unique(combinations)
                                groupName=obj.info.groupNames{comb};
                                tempData = obj.sourceResults.roiData.(groupName).(property).(timeName).(freq);
                                all_data = [all_data reshape(tempData(rois2plot,:),1,[])];
                            end
                        end
                    end
                end
                isnormal = ~adtest(all_data);
            end
            obj.sourceResults.isnormal=isnormal;
            if toPlot
                hm = headModel.loadFromFile(hmFile);
                T = hm.indices4Structure(hm.atlas.label);
                T = double(T)';
            end
            for p=1:length(vars2plot)
                property=vars2plot{p}; % get property name
                for comb = 1:size(combinations, 1)
                    group1=obj.info.groupNames{combinations(comb, 1)};
                    group2=obj.info.groupNames{combinations(comb, 2)};    
                    for t=1:length(times2plot)
                        timeName=times2plot{t};
                        plotdata=[];pvals=[];
                        for f=1:length(freq2plot)
                            freq=freq2plot{f};
                            % average in time dimension
                            s1 = squeeze(nanmean(obj.getGroupData(group1,property,freq,timeName),3));
                            s2 = squeeze(nanmean(obj.getGroupData(group2,property,freq,timeName),3));
                            pvals(:,f)=obj.calGroupSig(s1,s2,obj.info.experimentalDesign,isnormal); % get 1x n channel of p values

                            if FDRflag
                                tempP = fdr(pvals(rois2plot,f)); % only fdr correct ROIs of interest
                                pvals(rois2plot,f) = tempP;
                            end
                            unusedROI = setdiff(1:size(pvals,1),rois2plot);
                            pvals(unusedROI,f)=nan; % remove p values of roi not used
                            obj=addSourceSig(obj,property, timeName, append(group2,"_",group1), freq, pvals(:,f));
                            plotdata(:,f)=nanmean(s2,2)-nanmean(s1,2); % average across subjects and subtract
                            plotdata(unusedROI,f)=0;
                        end
                        if toPlot
                            tohide = T'*(pvals>0.05);
                            X=T'*plotdata;
                            X(tohide==1)=0;
                            plot68roi(hm,X , 1,freq2plot)
                            hAx = axes('Position', [0, 0, 1, 1], 'Visible', 'off');
                            text(0.5, 1, property+" "+timeName+" "+group2+"-"+group1, 'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
                        end
                    end
                end
            end
        end
        function obj = calConnectivity(obj,netwrk,varargin)
            % calculateConnectivity Compute functional connectivity between brain regions or networks
            %
            %   connMatrix = calculateConnectivity()
            %
            %   This function computes connectivity using one of two approaches:
            %     1. 'roi' method:
            %         - Compute ROI-to-ROI connectivity matrix
            %         - Apply Fisher Z-transform
            %         - Average connectivity within each network
            %         - Convert back to connectivity values
            %     2. 'net' method:
            %         - Average ROI time series within each network
            %         - Compute network-to-network connectivity directly
            %         - Intra-network connectivity defaults to ROI-to-ROI method
            %
            %   Inputs:
            %     obj         - SourceObject
            %     netwrk      - Struct array defining brain networks. Each element must contain:
            %                 .name : Name of the network (e.g., 'FPN', 'DMN')
            %                 .roi  : Vector of ROI indices that define the network
            %   Optional: 
            %     method      'net' or 'roi'
            %
            %   Notes:
            %     - Intra-network connectivity is always computed using ROI-level method
            %     - Assumes data has been preprocessed (e.g., detrended, filtered)
            
            properties=obj.info.variables;

            cfg = parseCfgOrArgs(obj, varargin{:});
            % Use fields directly
            method  = cfg.netConMethod;

            timeNames = fieldnames(obj.info.timeRange);
            N = length(obj.info.groupNames);            
            if ~isfield(obj.info,'netwrk')
                obj.info.netwrk=netwrk;
            end
            disp('-----------Calulating network connectivity and saving in obj.netConnectivity-----------')


            for p=1:length(properties)
                property=properties{p};
                
                for t=1:length(timeNames)
                    timeName=timeNames{t};
                    for f = 1:length(obj.info.freq_list)
                        freq=obj.info.freq_list{f};
                        % Generating plot data
                        data=[];
                        for n=1:N % iterates through groups
                            groupdata = squeeze(getGroupData(obj,obj.info.groupNames{n},property,freq,timeName)); % size roi x time x sub
                            all_sub_net=[];
                            for sub=1:size(groupdata,3)
                                subdata = squeeze(groupdata(:,:,sub)); % size roi x time
                                net_mat=[];
                                for a =  1:length(netwrk)
                                    for b = 1:length(netwrk)
               
                                        if a==b | strcmp(method,'roi') % intra network connectivity by averaging roi
                                            netrois=netwrk(a).roi;
                                            intraR_mat=corr(subdata(netrois,:)','rows','complete','Type','Spearman'); 
                                            intraFisherz=atanh(intraR_mat);
                                            intraFisherz(intraFisherz==inf)=nan;
                                            net_mat(a,b) = tanh(nanmean(intraFisherz,'all'));
                                        elseif strcmp(method,'net')
                                            net1=netwrk(a).roi;
                                            net2=netwrk(b).roi;
                                            net1Data = nanmean(subdata(net1,:));
                                            net2Data = nanmean(subdata(net2,:));
                                            net_mat(a,b)=corr(net1Data',net2Data','rows','complete','Type','Spearman');
                                        else 
                                            error("Connectivity Method must be 'roi' or 'net'")
                                        end
                                    end
                                end
                                all_sub_net(:,:,sub)=net_mat;
                            end
                            obj.netConnectivity.(obj.info.groupNames{n}).(property).(timeNames{t}).(obj.info.freq_list{f}) = all_sub_net;
                        end
                    end
                end
            end
        end
        function plotNetConnectivity(obj,varargin)
            %   [...] = analyzeNetConnectivity(...,'PARAM1',VAL1,'PARAM2',VAL2,...)
            %   Plots and anlyzes network connectivity data. Must run calConnectivity()
            %   first. Each dataset will generate 3 plots, group 1, group 2, and
            %   group1-group2. Individual groups will test connectivity wrt 0.
            %   Group1-group2 will test differences between groups.
            %   Normality test will determine if ttest/ttest2 or signrank/ranksum
            %
            %   Specifies additional parameters and their values. 
            %   Will only plot conditions that satisfy all conditions specified. 
            %   Default run will plot all combinations
            
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
            %         'isnormal'    'auto' 1 0 |
            %                        indicates distribution of the data for
            %                        significance testing. normal (1 flag) will use
            %                        ttest/ttest2. not normal (0 flag) will
            %                        use signrank/ranksum test. 'auto' will
            %                        determine if the data is normal with
            %                        adtest. 'auto' is default.

            cfg = parseCfgOrArgs(obj, varargin{:});
            % Use fields directly
            combinations = cfg.combinations;
            vars2plot  = cfg.vars2plot;
            freq2plot  = cfg.freq2plot;
            times2plot  = cfg.times2plot;
            FDRflag = cfg.FDRflag;
            isnormal=cfg.isnormal;

            netwrk=obj.info.netwrk;
            
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
                                tempData = obj.netConnectivity.(groupName).(property).(timeName).(freq);
                                all_data = [all_data reshape(tempData,1,[])];
                            end
                        end
                    end
                end
                isnormal = ~adtest(atanh(all_data));
            end
            
            for p=1:length(vars2plot)
                property=vars2plot{p};
            
                for t=1:length(times2plot)
                    timeName=times2plot{t};
                    for comb = 1:size(combinations, 1)
                        figure
                        count=1;
                        subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.1], [0.1 0.1], [0.1 0.1]);
                        group1=obj.info.groupNames{combinations(comb, 1)};
                        group2=obj.info.groupNames{combinations(comb, 2)};
                        for f = 1:length(freq2plot)
                            freq=freq2plot{f};
                            s1 = squeeze(obj.netConnectivity.(group1).(property).(timeName).(freq));
                            s2 = squeeze(obj.netConnectivity.(group2).(property).(timeName).(freq));
                            pvalsDiff=[]; pvals1=[];pvals2=[];
                            for a=1:size(s1,1)
                                pvalsDiff(:,a)=obj.calGroupSig(squeeze(atanh(s1(a,:,:))),squeeze(atanh(s2(a,:,:))),obj.info.experimentalDesign,isnormal);
                                pvals1(:,a) = obj.calGroupSig(squeeze(atanh(s1(a,:,:))),zeros(size(s1,[2,3])),'paired',isnormal);
                                pvals2(:,a) = obj.calGroupSig(squeeze(atanh(s2(a,:,:))),zeros(size(s2,[2,3])),'paired',isnormal);
                            end
                            
                            if FDRflag==1
                                pvalsDiff = fdr_matCorrect(pvalsDiff);
                                pvals1 = fdr_matCorrect(pvals1);
                                pvals2 = fdr_matCorrect(pvals2);
                            end
            
                            subplot(length(freq2plot),3,count)
                            plotNetConn(tanh(nanmean(atanh(s1),3)), pvals1,netwrk)
                            
                            count=count+1;
            
                            if f==1
                                title(group1)
                            end
            
                            ylabel(freq,'fontweight','bold','fontsize',12)
            
                            subplot(length(freq2plot),3,count)
                            plotNetConn(tanh(nanmean(atanh(s2),3)), pvals1,netwrk)
                            count=count+1;
                            if f==1
                                title(group2)
                            end
            
                            subplot(length(freq2plot),3,count)
                            plotNetConn(tanh(nanmean(atanh(s2)-atanh(s1),3)), pvalsDiff,netwrk)
                            count=count+1;
                            if f==1
                                title(group2 + "-" + group1)
                            end
                        end
                    end
                    sgtitle(property+" "+timeName)
                end
            end
        end
        function obj = addSourceSig(obj,property, timeName, groupDiff, freq, pvals)
            obj.sourceResults.sigROIsP.(property).(timeName).(groupDiff).(freq)=pvals;
            roi2add = obj.info.roi(pvals<=0.05);
            %%% adding ROI's to results
            % Check if the field already exists
            if isfield(obj.sourceResults, 'sigROIs') && isfield(obj.sourceResults.sigROIs, property) && ...
                isfield(obj.sourceResults.sigROIs.(property), timeName) && isfield(obj.sourceResults.sigROIs.(property).(timeName), freq)
                % Get the existing cell array
                existingArray = obj.sourceResults.sigROIs.(property).(timeName).(freq);
                % Check if the new letter is different from the existing ones
                roi2add = setdiff(roi2add, existingArray);
                obj.sourceResults.sigROIs.(property).(timeName).(freq) = [existingArray; roi2add];
            elseif ~isempty(roi2add)
                % Create a new cell array with the new letter
                obj.sourceResults.sigROIs.(property).(timeName).(freq) = roi2add;
            end

            if ~isfield(obj.sourceResults.sigROIs, property)
                obj.sourceResults.sigROIs.(property) = struct();
            end
        end


        %--------------END OF CLASS METHODS-----------------------%
    end
    %--------------END OF CLASS-----------------------%
end
