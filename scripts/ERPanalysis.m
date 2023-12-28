classdef ERPanalysis
    properties
        groups;
        groupNames;
        freq_list;
        variables;
        baselineIDX;
        timeAxis;
    end
    methods
        function obj = ERPanalysis(groupPaths,variables,baselineTime)
            disp('Initializing object')
            % groupPaths is a cell array of file paths
            % Initializing the properties to the object
            
            for i = 1:length(groupPaths)
                obj.groups{i} = load(groupPaths{i});
            end
            obj.variables=variables;
            obj.timeAxis=obj.groups{1}.timeAxis;
            [~,baselineIDX]=min(abs(baselineTime'-obj.timeAxis)');
            obj.baselineIDX=baselineIDX;
            obj.groupNames=cellfun(@(x) x.groupname, obj.groups,'uniformoutput',false);
            obj.freq_list=obj.groups{1}.freq_list;
            
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

%% Plotting scalpmap with significant electrodes
function plotSigTopo(s1,s2,chanlocs,key)
    % s1 and s2 are chan x sub size
    
    plotdata=nanmean(s2,2)-nanmean(s1,2);
    
    new_chan=topoSignificance(s1,s2,chanlocs,{'.','+'});
    topoplot(plotdata,new_chan,'headrad','rim','electrodes','labels'); 
    colorbar(); 
    caxis([min(plotdata),max(plotdata)])
end

%% Function to calculate which electrodes are significant btw two groups
function chanlocs=topoSignificance(s1,s2,chanlocs,key)

    % pre and post are data matrices size N subjects x C channels
    % key is a two element cell array with markers for non-sig and sig
    % new_chan is the chan_locs
    % Perform t-tests
    
    s1=squeeze(s1);
    s2=squeeze(s2);
    for chan = 1:length(chanlocs)
        % Get non-NaN indices
        [~,p(chan)] = ttest2(s1(chan,:), s2(chan,:));
        if isnan(p(chan))
            disp('asdfawrsasdbaweasbas')
            waitforbuttonpress
        end
        
    end
    %[~,p]=fdr(reshape(p,[1,numel(p)]),0.05);
    mask=key((p<0.05)+1);
    [chanlocs.labels]=mask{:};
end
        
    