classdef DataAnalysis
    %%% TODO %%%
    % Add in option for paired or 2 sample ttest
    % Add in option to choose what pairings to run 
    properties
        info
    end
    methods
        function obj = DataAnalysis(info,baselineTime, timeRange)
            % obj           | formatted data structure with properties of DATA
            % baseline time | vector with start and end time of baseline in mS 
            disp('Initializing object')
            
            % and info
            % Initializing the properties to the object
            %-------------------------------------------------------------------%

            obj.info=info;
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
        

        function obj = standardPipeline(obj)
            properties = obj.info.variables;
            disp('---------------------------------')
            for j = 1:length(properties)
                property = properties{j};
                disp("Processing "+property)
                combined = obj.combine_groups(property);
                cleaned = obj.rej5SD(combined);
                baselineCorrected = obj.baselineCorrection(cleaned,obj.info.baselineIDX);
                obj = obj.split_combined( baselineCorrected,property);
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
            % freq | string of frequency band
            % chans | vector of chans to include, if 'all' or empty indicates use all channels
            % timeName | string of timeName

            if nargin < 6
                chans = 'all';
            end

            freqIdx = find(strcmp(obj.info.freq_list, freq));
            timeIdx = obj.info.timeIDX.(timeName);

            if isempty(chans) || (ischar(chans) && strcmp(chans, 'all'))
                chans = 1:size(obj.DATA.(group).(property), 2); % assuming channels are the 2nd dimension
            end
            s = obj.DATA.(group).(property)(freqIdx,chans,timeIdx(1):timeIdx(2),:);
        end
        %--------------END OF STATIC METHODS-----------------------%
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
        %--------------END OF STATIC METHODS-----------------------%
    end 
    %------------------END OF CLASS-------------------------%
end

