classdef ERPanalysis
    properties
        groups;
        variables;
        baselineIDX;
    end
    methods
        function obj = ERPanalysis(groupPaths,variables,baselineIDX)
            disp('Initializing object')
            % groupPaths is a cell array of file paths
            % Initializing the properties to the object
            for i = 1:length(groupPaths)
                obj.groups{i} = load(groupPaths{i});
            end
            obj.variables=variables;
            obj.baselineIDX=baselineIDX;
        end
        function obj = cleanDatasets(obj)
            disp('Removing missing data')
            % Getting rid of the nan's and missing subjects
            for i = 1:length(obj.groups)
                obj.groups{i} = cleanDataset(obj.groups{i});
            end
        end
        function combined = combine_groups(obj,property)
            disp('Combining groups')
            % stacking groups together for 5SD
            data = cell(1, length(obj.groups));
            for i = 1:length(obj.groups)
                data{i} = obj.groups{i}.(property);
            end
            combined =cat(4,data{:});
        end
 
        function obj = split_combined(obj, groupedData,property)
            disp('Splitting groups')
            % Reordering the data back into groups in obj
            for i=1:length(obj.groups)
                obj.groups{i}.(property)=groupedData(:,:,:,1:length(obj.groups{i}.subjectcoll));
                groupedData(:,:,:,1:length(obj.groups{i}.subjectcoll))=[];
            end
        end
        
        function obj = standardPipeline(obj)
            properties = obj.variables;
            disp('---------------------------------')
            for j = 1:length(properties)
                property = properties{j};
                disp("Processing "+property)
                combined = combine_groups(obj,property);
                cleaned = rej5SD(combined);
                baselineCorrected = baselineCorrection(cleaned,obj.baselineIDX);
                obj = split_combined(obj, baselineCorrected,property);
                disp('---------------------------------')
            end
        end
    end
end

function baselineCorrected = baselineCorrection(cleaned,baselineIDX)
    disp('Baseline Correction')
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

function group = cleanDataset(group)
    missingidx = find(isnan(squeeze(sum(sum(sum(group.EVgain, 1), 2), 3))));
    group.EVgain(:,:,:,missingidx) = [];
    group.EVloss(:,:,:,missingidx) = [];
    group.GLbias(:,:,:,missingidx) = [];
    group.subjectcoll(missingidx) = [];
end