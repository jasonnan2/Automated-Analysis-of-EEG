classdef SourceObject < DataAnalysis
    properties
        DATA
        sourceResults
        netConnectivity
    end
    methods
        function obj = SourceObject(DATA, info, baselineTime, timeRange)
            % Call superclass constructor
            obj = obj@DataAnalysis( info, baselineTime, timeRange);
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
            %         'groups2plot' array defining which groups show
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
            
            pnames = {'vars2plot','freq2plot','times2plot','groups2plot','FDRflag','isnormal'};
            dflts  = {obj.info.variables,obj.info.freq_list, fieldnames(obj.info.timeRange),1:length(obj.info.groupNames),1,'auto'};
            [vars2plot,freq_list,timeNames,groups2plot,FDRflag,isnormal] = parseArgs(pnames,dflts,varargin{:});
            netwrk=obj.info.netwrk;

            % normality test
            if strcmp(isnormal,'auto')
                % Get all Data
                all_data=[];
                for p=1:length(vars2plot)
                    property=vars2plot{p};
                    for t=1:length(timeNames)
                        timeName=timeNames{t};
                        for f = 1:length(obj.info.freq_list)
                            freq=obj.info.freq_list{f};
                            for comb = unique(combinations)
                                groupName=obj.info.groupNames{comb};
                                tempData = obj.sourceResults.netData.(groupName).(property).(timeName).(freq);
                                all_data = [all_data reshape(tempData(rois2plot,:),1,[])];
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
                for t=1:length(timeNames)
                    timeName=timeNames{t};
                    for f = 1:length(obj.info.freq_list)
                        freq=freq_list{f};
                        % Generating plot data
                        data=[];sem=[];
                        
                        for comb = 1:length(groups2plot)
                            group=obj.info.groupNames{groups2plot(comb)};
                            tempData = obj.sourceResults.netData.(group).(property).(timeName).(freq); % size net x groups
                            data(:,comb) = nanmean(tempData,2);
                            sem(:,comb) = std(tempData,[],2)/sqrt(size(tempData,2));
                        end
                        
                        count=count+1;
                        subplot(length(timeNames),length(freq_list),count)
                        
                        obj.plotErrBar(data,sem)
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

        function obj=calRoiData(obj)

            % Outputs:
            %
            %   Results are saved into the following fields of the input object `obj`:
            %


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
            
            N = length(obj.info.groupNames);
            pnames = {'rois2plot','vars2plot','freq2plot','times2plot','combinations','FDRflag','toPlot','isnormal'};
            dflts  = {length(obj.info.roi),obj.info.variables,obj.info.freq_list, fieldnames(obj.info.timeRange),nchoosek(1:N,2),1,1,'auto'};
            [rois2plot,vars2plot,freq_list,timeNames,combinations,FDRflag,toPlot,isnormal] = parseArgs(pnames,dflts,varargin{:});
            % normality test
            if strcmp(isnormal,'auto')
                % Get all Data
                all_data=[];
                for p=1:length(vars2plot)
                    property=vars2plot{p};
                    for t=1:length(timeNames)
                        timeName=timeNames{t};
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
                hm = headModel.loadFromFile('headModel_templateFile_newPEBplus.mat');
                T = hm.indices4Structure(hm.atlas.label);
                T = double(T)';
            end
            for p=1:length(vars2plot)
                property=vars2plot{p}; % get property name
                for comb = 1:size(combinations, 1)
                    group1=obj.info.groupNames{combinations(comb, 1)};
                    group2=obj.info.groupNames{combinations(comb, 2)};    
                    for t=1:length(timeNames)
                        timeName=timeNames{t};
                        plotdata=[];pvals=[];
                        for f=1:length(freq_list)
                            freq=freq_list{f};
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
                            obj=addSourceSig(obj,property, timeName, groupDiff, freq, pvals(:,f));
                            plotdata(:,f)=nanmean(s2,2)-nanmean(s1,2); % average across subjects and subtract
                            plotdata(unusedROI,f)=0;
                        end
                        if toPlot
                            tohide = T'*(pvals>0.05);
                            X=T'*plotdata;
                            X(tohide==1)=0;
                            plot68roi(hm,X , 1,freq_list)
                            hAx = axes('Position', [0, 0, 1, 1], 'Visible', 'off');
                            text(0.5, 1, property+" "+timeName+" "+group2+"-"+group1, 'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
                        end
                    end
                end
            end
        end
        function obj = calConnectivity(obj,netwrk,method)
            % calculateConnectivity Compute functional connectivity between brain regions or networks
            %
            %   connMatrix = calculateConnectivity(data, method, roiLabels, networkLabels)
            %
            %   This function computes connectivity using one of two approaches:
            %     1. 'avg_roi' method:
            %         - Compute ROI-to-ROI connectivity matrix
            %         - Apply Fisher Z-transform
            %         - Average connectivity within each network
            %         - Convert back to connectivity values
            %     2. 'network' method:
            %         - Average ROI time series within each network
            %         - Compute network-to-network connectivity directly
            %         - Intra-network connectivity defaults to ROI-to-ROI method
            %
            %   Inputs:
            %     obj            - SourceObject
            %     method         - 'roi' or 'network'
            %
            %   Notes:
            %     - Intra-network connectivity is always computed using ROI-level method
            %     - Assumes data has been preprocessed (e.g., detrended, filtered)
            
            properties=obj.info.variables;
            timeNames = fieldnames(obj.info.timeRange);
            N = length(obj.info.groupNames);
            combinations = nchoosek(1:N,2);
            
            if ~isfield(obj.info,'netwrk')
                obj.info.netwrk=netwrk;
            end
            obj.info.netConnectivity_Method=method;
            
            for p=1:length(properties)
                property=properties{p};
                
                for t=1:length(timeNames)
                    timeName=timeNames{t};
                    timeIdx = obj.info.timeIDX.(timeName);
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
            
            N = length(obj.info.groupNames);
            pnames = {'vars2plot','freq2plot','times2plot','combinations','FDRflag'};
            dflts  = {obj.info.variables,obj.info.freq_list, fieldnames(obj.info.timeRange),nchoosek(1:N,2),1};
            [vars2plot,freq_list,timeNames,combinations, FDRflag] = parseArgs(pnames,dflts,varargin{:});
            netwrk=obj.info.netwrk;
            % Get all Data
            all_data=[];
            for p=1:length(vars2plot)
                property=vars2plot{p};
                for t=1:length(timeNames)
                    timeName=timeNames{t};
                    for f = 1:length(obj.info.freq_list)
                        freq=obj.info.freq_list{f};
                        for comb = unique(combinations)
                            groupName=obj.info.groupNames{comb};
                            tempData = obj.netConnectivity.(groupName).(property).(timeName).(freq);
                            all_data = [all_data reshape(tempData,1,[])];
                        end
                    end
                end
            end
            
            isnormal = ~adtest(atanh(all_data));
            
            for p=1:length(vars2plot)
                property=vars2plot{p};
            
                for t=1:length(timeNames)
                    timeName=timeNames{t};
                    figure
                    count=1;
                    subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.1], [0.1 0.1], [0.1 0.1]);
            
            
                    for comb = 1:size(combinations, 1)
                            group1=obj.info.groupNames{combinations(comb, 1)};
                            group2=obj.info.groupNames{combinations(comb, 2)};
                            
                        for f = 1:length(obj.info.freq_list)
                            freq=obj.info.freq_list{f};
                            s1 = squeeze(obj.netConnectivity.(group1).(property).(timeName).(freq));
                            s2 = squeeze(obj.netConnectivity.(group2).(property).(timeName).(freq));
                            pvalsDiff=[];
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
            
                            subplot(length(obj.info.freq_list),3,count)
                            plotNetConn(tanh(nanmean(atanh(s1),3)), pvals1,netwrk)
                            
                            count=count+1;
            
                            if f==1
                                title(group1)
                            end
            
                            ylabel(freq,'fontweight','bold','fontsize',12)
            
                            subplot(length(obj.info.freq_list),3,count)
                            plotNetConn(tanh(nanmean(atanh(s2),3)), pvals1,netwrk)
                            count=count+1;
                            if f==1
                                title(group2)
                            end
            
                            subplot(length(obj.info.freq_list),3,count)
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
        end


        %--------------END OF CLASS METHODS-----------------------%
    end
    %--------------END OF CLASS-----------------------%
end
