classdef DataAnalysis
    properties
        info
        neuralBehMdl
    end
    methods
        function obj = DataAnalysis(info,baselineTime, timeRange,cfg)
            % obj           | formatted data structure with properties of DATA
            % baseline time | vector with start and end time of baseline in mS 
            % cfg           | configuration struct for definitions. If left
            %                 empty defaults are used

            disp('Initializing object')
           
            % Initializing the properties to the object
            %-------------------------------------------------------------------%
            obj.info=info;
            [~,baselineIDX]=min(abs(baselineTime'-obj.info.timeAxis)');
            obj.info.baselineIDX=baselineIDX;
            obj.info.timeRange = timeRange;
            timeNames = fieldnames(obj.info.timeRange);
            for t = 1:length(timeNames)
                [startIdx, endIdx] = obj.getTimeIndices(timeRange, timeNames{t});
                timeIDX.(timeNames{t}) = [startIdx endIdx];
            end
            obj.info.timeIDX=timeIDX;
            obj = obj.setDefaultCfg(cfg);
        end

        function obj = setDefaultCfg(obj,cfg)
            N=length(obj.info.groupNames);
            if ~isfield(cfg, 'freq2plot');     cfg.freq2plot = obj.info.freq_list; end
            if ~isfield(cfg, 'rois2plot');     cfg.rois2plot = length(obj.info.roi); end
            if ~isfield(cfg, 'times2plot');    cfg.times2plot = fieldnames(obj.info.timeRange); end
            if ~isfield(cfg, 'vars2plot');     cfg.vars2plot = obj.info.variables; end
            if ~isfield(cfg, 'combinations');  cfg.combinations = nchoosek(1:N,2); end
            if ~isfield(cfg, 'errorType');     cfg.errorType = 'none'; end
            if ~isfield(cfg, 'FDRflag');       cfg.FDRflag = 1; end
            if ~isfield(cfg, 'toPlot');        cfg.toPlot = 1; end
            if ~isfield(cfg, 'isnormal');      cfg.isnormal = 'auto'; end
            if ~isfield(cfg, 'chans2plot');    cfg.chans2plot = 'all'; end
            if ~isfield(cfg, 'color_list');    cfg.color_list = {'r','b','g','m','k','c','y'}; end
            if ~isfield(cfg, 'groups2plot');   cfg.groups2plot = 1:length(obj.info.groupNames); end
            if ~isfield(cfg, 'hmFile');        cfg.hmFile = ''; end
            if ~isfield(cfg, 'netConMethod');  cfg.hmFile = 'net'; end
            obj.info.cfg = cfg;
        end
        function [cfgOut] = parseCfgOrArgs(obj, varargin)
            % parseCfgOrArgs  Unified parser for cfg struct or name-value pairs
            % 
            % Usage:
            %   cfg = parseCfgOrArgs(obj, varargin)
            %
            % Automatically fills missing fields from obj.info defaults
            % Default values from obj.info
            defaults = obj.info.cfg;
        
            % If a single struct is passed (cfg), convert to name-value pairs
            if numel(varargin) == 1 && isstruct(defaults)
                cfgStruct = defaults;
                varargin = reshape([fieldnames(cfgStruct), struct2cell(cfgStruct)]', 1, []);
            end
        
            % Parameter names and defaults in cell form for parseArgs
            pnames = fieldnames(defaults);
            dflts = struct2cell(defaults);
            
            % Safe unpacking from varargout-based parseArgs
            parsed = cell(1, numel(pnames));
            [parsed{:}] = parseArgs(pnames, dflts, varargin{:});
            
            cfgOut = cell2struct(parsed(:), pnames(:), 1);
        end
        
        function [startIdx, endIdx] = getTimeIndices(obj, timeRange, timeName)
            [~, startIdx] = min(abs(obj.info.timeAxis - timeRange.(timeName)(1)));
            [~, endIdx] = min(abs(obj.info.timeAxis - timeRange.(timeName)(2)));
        end
        
        function obj = cleanDatasets(obj)
            disp('Removing missing data')
            % Getting rid of the nan's and missing subjects
            for i = 1:length(obj.info.groupNames)
                obj.DATA.(obj.info.groupNames{i}) = obj.cleanDataset(obj.DATA.(obj.info.groupNames{i}));
            end
        end
        
        function data = cleanDataset(obj,data)
            % Get field names from the data structure
            fieldNames = obj.info.variables;
            allmissing={};
            % Loop over each field of the data structure
            for v = 1:numel(fieldNames)
                % Get the field data
                fieldData = data.(fieldNames{v});
                % Find missing data in the field data, treats all 0 as
                % missing data too
                missingidx = find(isnan(squeeze(sum(sum(sum(fieldData, 1), 2), 3))) | squeeze(sum(sum(sum(abs(fieldData), 1), 2), 3))==0);
                % Remove missing data from the field data
                fieldData(:,:,:,missingidx) = [];
                % Assign the cleaned field data back to the data structure
                data.(fieldNames{v}) = fieldData;
                data.missingSubs{v}=[fieldNames{v} reshape(data.subList(missingidx), 1, [])];
                %data.subList(missingidx)=[];
            end
        end
        

        function obj = standardProcessing(obj,varargin)
            %        Parameter  Value
            %         'outlier'    String to define either 5SD or 3MAD
            %                      between subjects outlier removal
            %                      Default is 5SD
            %         'calBaseline' 1/0 binary if you want to calculate
            %                       baseline. Default is true
            
            pnames = {'outlier','calBaseline'};
            dflts  = {'5SD',1};
            [outlier,calBaseline] = parseArgs(pnames,dflts,varargin{:});

            properties = obj.info.variables;
            disp('---------------------------------')
            for j = 1:length(properties)
                property = properties{j};
                disp("Processing "+property)
                combined = obj.combine_groups(property);
                
                if strcmp(outlier,'5SD')
                    cleaned = obj.rej5SD(combined);
                elseif strcmp(outlier,'3MAD')
                    cleaned = obj.rej3MAD(combined);
                elseif strcmp(outlier,'none')
                    cleaned=combined;
                end
                
                baselineCorrected = obj.baselineCorrection(cleaned,obj.info.baselineIDX);
                if calBaseline
                    obj = obj.split_combined( baselineCorrected,property);
                else
                    obj = obj.split_combined( cleaned,property);
                end
                
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
        function s = getGroupData(obj,group,property,freq,timeName,chans)
            % group | name of the group
            % property | string of variable
            % freq | string of frequency band. if 'all', gets all freq
            % chans | vector of chans to include, if 'all' or empty indicates use all channels
            % timeName | string of timeName

            if nargin < 6
                chans = 'all';
            end
            
            if strcmp(freq,'all')
                freqIdx=1:size(obj.DATA.(group).(property),1);
            else
                freqIdx = find(strcmp(obj.info.freq_list, freq));
            end
            timeIdx = obj.info.timeIDX.(timeName);

            if isempty(chans) || (ischar(chans) && strcmp(chans, 'all'))
                chans = 1:size(obj.DATA.(group).(property), 2); % assuming channels are the 2nd dimension
            end
            s = obj.DATA.(group).(property)(freqIdx,chans,timeIdx(1):timeIdx(2),:);
        end
        
        function obj=applyFunc(obj,func)
            % func is function handel in the format of @(x)function(x,optional arguments)
            for p=1:length(obj.info.variables)
                property=obj.info.variables{p};
                for g=1:length(obj.info.groupNames)
                    group=obj.info.groupNames{g};
                    for f=1:length(obj.info.freq_list)
                        data = squeeze(obj.DATA.(group).(property)(f,:,:,:));
                        for c=1:size(data,1)
                            x =func(squeeze(data(c,:,:))); % size time x N subs
                            obj.DATA.(group).(property)(f,c,:,:) = x;
                            x=[];
                        end
                    end
                end
            end
        end
        
        function obj = calSigTbl(obj,chanTypes)
            % Extracts significantly different variables and saves into 
            % obj.results.sigValues.(property).(time)
            %
            % chanTypes | 'all' |'sig' string value indicating what channels to save as
            %             table. default is to save only significant channels 'sig' 
            
            if nargin<2
                chanTypes='sig';
            else
                chanTypes=lower(chanTypes);
            end
            if ~(strcmp(chanTypes,'all') | strcmp(chanTypes,'sig'))
                error("chanTypes must be specified as 'all' or 'sig'"+chanTypes+": was gotten instead")
            end
            
            fnames = fieldnames(obj);
            results = fnames{contains(fnames,'Results')};
            obj.(results).neuralBehMdl = struct();
            N=numel(fieldnames(obj.DATA));
            timeNames = fieldnames(obj.info.timeRange);
            
            if isfield(obj.(results),'sigElectrodes')
                sigChans=obj.(results).sigElectrodes;
                chanSet={obj.info.chanlocs.labels};
            elseif isfield(obj.(results),'sigROIs')
                sigChans=obj.(results).sigROIs;
                chanSet=obj.info.roi;
            end

            varNames=fieldnames(sigChans);
            times2plot = fieldnames(obj.info.timeRange);
            for p=1:length(varNames)
                property=varNames{p};
                p_idx= find(strcmp(obj.info.variables, property));
                time_list = intersect(fieldnames(sigChans.(property)),times2plot);

                for t=1:length(time_list)
                    timeName=time_list{t};
                    sig=sigChans.(property).(timeName);
                    
                    if strcmp(chanTypes,'sig')
                        freq_list=fieldnames(sig);
                    elseif strcmp(chanTypes,'all')
                        freq_list=obj.info.freq_list;
                    end
                    
                    % get longest electrode map
                    labels=[]; tbldata=[];
                    
                    for f=1:length(freq_list)
                        
                        if strcmp(chanTypes,'sig')
                            allChans=sig.(freq_list{f});
                        elseif strcmp(chanTypes,'all')
                            allChans=chanSet;
                        else
                            error("chanTypes must be specified as 'all' or 'sig'"+chanTypes+": was gotten instead")
                        end
                        
                        for c=1:length(allChans)
                            freq=freq_list{f};
                            chan=allChans{c};
                            elecIdxs = find(strcmp({obj.info.chanlocs.labels},chan ));
                            if isempty(elecIdxs)
                                elecIdxs = find(strcmp(obj.info.roi,chan ));
                            end
                            data=[];allsubs=[];groupname=[];
                            for n=1:N % iterate groups
                               subdata = [squeeze(nanmean(obj.getGroupData(obj.info.groupNames{n},property,freq,timeName,elecIdxs),3))];
                               data =[data; subdata];
                               subsOrig = obj.DATA.(obj.info.groupNames{n}).subList;
                               missing = obj.DATA.(obj.info.groupNames{n}).missingSubs{p_idx}(2:end);
                               subsOrig(contains(subsOrig,missing))=[];
                               allsubs=[allsubs,reshape(subsOrig, 1, [])];
                               groupname=[groupname repmat(obj.info.groupNames(n),1,length(subdata))];
                            end
                            tbldata = [tbldata data];
                            labels=[labels strcat(property,"_",timeName,"_",freq,"_",chan)];
                        end
                    end
                    tbl = [table(allsubs', groupname', 'VariableNames', {'subID', 'group'}) array2table(tbldata, 'VariableNames', labels)];
                    obj.(results).sigValues.(property).(timeName)=tbl;
                end
            end
        end
        function obj=NeurBehMdl(obj,neuralVar,behTbl,keyColumnName,baseModel,modelname)
            % neuralVar     | cell array of neural variable to test against
            %                 base model
            % behTbl        | table with behavior/redcap data which you want to correlate
            %                 with neural data. Subject names must be included and in  
            %                 the exact same format as the subList in neural data
            % keyColumnName | string input denoting column name for subjectID
            % baseModel     | model definition for fitlm, neural data will be appened on 
            %                 ex: basemodel='target ~ var1+' ->  basemodel='target ~ +var1+neural'
            %                 ex: basemodel='target ~ var1*' ->  basemodel='target ~ var1*neural'

            fnames = fieldnames(obj);
            results = fnames{contains(fnames,'Results')};
            timeNames = fieldnames(obj.info.timeRange);
            for p=1:length(neuralVar)
                property = neuralVar{p};
                if nargin<5
                    modelname=makeValidFieldName(append(baseModel,"_",property));
                else
                    modelname=append(modelname,"_",property);
                end
                sigStruct=obj.(results).sigValues.(property);
                timeNames=fieldnames(sigStruct);
                NBtbl=[];
                for t=1:length(timeNames)
                    neuralTbl = sigStruct.(timeNames{t});
                    neuralTbl.Properties.VariableNames = strrep(neuralTbl.Properties.VariableNames, ' ', ''); % removes any spaces
                    varNames = setdiff(neuralTbl.Properties.VariableNames,{'subID','group', behTbl.Properties.VariableNames{:}});
                    % last cleaning up of subject names
                    neuralTbl.subID=lower(neuralTbl.subID);
                    behTbl.(keyColumnName)=lower(behTbl.(keyColumnName));
                    neuralTbl.subID = cellfun(@(x) erase(x, ['_' extractAfter(x, '_')]), neuralTbl.subID, 'UniformOutput', false);
                    % Join tables
                    tbldata = innerjoin(neuralTbl, behTbl, 'LeftKeys', 'subID', 'RightKeys', keyColumnName,'RightVariables', ...
                                setdiff(behTbl.Properties.VariableNames,neuralTbl.Properties.VariableNames));
                
                    if height(tbldata)~=height(neuralTbl)
                        warning('New Table is missing some subjects, please make sure behTbl has all the subejcts present')
                    end
                    obj.(results).sigValues.(property).(timeNames{t}) = tbldata;
                    NBtbl =[ NBtbl;obj.calNeurBehMdl(baseModel,tbldata,varNames)];
                end

                count=0;
                while isfield(obj.neuralBehMdl,modelname)
                    modelname = modelname+string(count);
                    count=count+1;
                end

                obj.neuralBehMdl.(modelname)= NBtbl;
            end
        end
        %--------------END OF CLASS METHODS-----------------------%
    end
    
    methods (Static)
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
        
        function combined = rej3MAD(combined)
            disp('3 MAD Outlier')
            for f = 1:size(combined, 1)
                for chan = 1:size(combined, 2)
                    data=squeeze(combined(f, chan, :));
                    out = isoutlier(data);
                    data(out)=nan;
                    combined(f, chan, :) = data;
%                     for t = 1:size(combined, 3)
%                         out = isoutlier(squeeze(combined(f, chan, t, :)));
%                         combined(f, chan, t, out) = nan;
%                     end
                end
            end
        end

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

        % Function to calculate which electrodes are significant btw two groups
        function [pvals,stats]=calGroupSig(s1,s2,experimentalDesign, isnormal)
            % s1 and s2 are data matrices size C channels x N subjects
            % experimentalDesign is string for 'paired' or 'twoSample'
            % Perform t-tests
            if nargin<4
                isnormal=1; % default to ttest
            end

            s1=squeeze(s1);
            s2=squeeze(s2);
            
            if iscolumn(s1)
                s1=s1';
            end
            if iscolumn(s2)
                s2=s2';
            end
            
            for chan = 1:size(s1,1)
                % Get non-NaN indices
                if isnormal
                    if strcmp(experimentalDesign,'paired')
                        [~,pvals(chan),~,stats(chan)] = ttest(s1(chan,:), s2(chan,:));
                    elseif strcmp(experimentalDesign,'twoSample')
                        [~,pvals(chan),~,stats(chan)] = ttest2(s1(chan,:), s2(chan,:));
                    end
                elseif ~isnormal
                    if strcmp(experimentalDesign,'paired')
                        [pvals(chan),~,stats(chan)] = signrank(s1(chan,:), s2(chan,:));
                    elseif strcmp(experimentalDesign,'twoSample')
                        [pvals(chan),~,stats(chan)] = ranksum(s1(chan,:), s2(chan,:));
                    end
                end

                if isnan(pvals(chan))
                    disp('asdfawrsasdbaweasbas')
                    waitforbuttonpress
                end
            end
        end
        function plotErrBar(data,sem,color_list)
            if nargin <3 
                color_list={'r','b','g','m','k','c','y'};
            end
            if isrow(data) | iscolumn(data)
                data = vertcat(data,nan(size(data)));
                sem = vertcat(sem,nan(size(sem)));
                xlim([0.5 1.5])
            end
            b=bar(data,'grouped'); hold on % data size N networks x G groups
            
            for i=1:length(b)
                b(i).FaceColor=color_list{i};
            end
            
            % Find the number of groups and the number of bars in each group
            [ngroups, nbars] = size(data);
            % Calculate the width for each bar group
            groupwidth = min(0.8, nbars/(nbars + 1.5));
            for i = 1:nbars
                % Calculate center of each bar
                x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
                errorbar(x, data(:,i), sem(:,i), 'k', 'linestyle', 'none');
            end
        end

        function neural_beh_tbl=calNeurBehMdl(baseModel,tbldata,varNames)
            % Returns robust fitlm
            % baseModel | model definition for fitlm, neural data will be appened on 
            %            ex: basemodel='target ~ var1+' ->  basemodel='target ~ +var1+neural'
            %            ex: basemodel='target ~ var1*' ->  basemodel='target ~ var1*neural'
            % tbldata | data for fitlm, can be found after running calSigTbl
            %           must include all variables in baseModel
            % varNames | neural data names from tbldata
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Output | table with p values only
            parts = strsplit(varNames{1},'_');
            neural_beh_tbl=table();
            pvals=[];
            for i=1:length(varNames)
                modelDef =  baseModel+varNames{i};
                model = fitlm(tbldata,modelDef,'RobustOpts','on'); 
                models{i}=model;
                neuralidx = find(strcmp(varNames{i},model.CoefficientNames));
                pvals(i) = model.Coefficients.pValue(neuralidx); % Check here to make sure its always the second
            end
            neural_beh_tbl.(baseModel)=varNames';
            %neural_beh_tbl.neuralP=pvals';
            neural_beh_tbl.model=models';
            %neural_beh_tbl(neural_beh_tbl.p>0.05,:)=[];
        end



        %--------------END OF STATIC METHODS-----------------------%
    end 
    %------------------END OF CLASS-------------------------%
end

