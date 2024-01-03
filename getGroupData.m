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
        chans = 1:size(obj.sourceData.(group).(property), 2); % assuming channels are the 2nd dimension
    end
    s = obj.sourceData.(group).(property)(freqIdx,chans,timeIdx(1):timeIdx(2),:);
end