


for n=1:N
    subdata = squeeze(nanmean(obj.getGroupData(obj.info.groupNames{n},property,freq,timeName,elecIdxs),3));
    data(n) = nanmean(subdata);
    sem(n) = std(subdata)/sqrt(length(subdata));
end
obj.plotErrBar(data,sem); % plot actual bars

% getting significance bars
[ngroups, nbars] = size(data);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for comb = 1:size(combinations, 1)
    group1=obj.info.groupNames{combinations(comb, 1)};
    group2=obj.info.groupNames{combinations(comb, 2)};
    s1=[];s2=[];
    s1 = squeeze(nanmean(obj.getGroupData(group1,property,freq,timeName,elecIdxs),3));
    s2 = squeeze(nanmean(obj.getGroupData(group2,property,freq,timeName,elecIdxs),3));
    pvals(comb,:)=obj.calGroupSig(s1,s2);
    group1Location = (1:ngroups) - groupwidth/2 + (2*combinations(comb, 1)-1) * groupwidth / (2*nbars);
    group2Location = (1:ngroups) - groupwidth/2 + (2*combinations(comb, 2)-1) * groupwidth / (2*nbars);
    A=[group1Location;group2Location]';
    groupingKey = mat2cell(A, ones(1, size(A, 1)), size(A, 2));
    sigstar(groupingKey,pvals(comb,:))
end