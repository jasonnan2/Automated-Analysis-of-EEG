%% Climate luckydoor Main analysis
% Jason Nan
% 12/19/2023
clear all;close all;clc
cd('/media/owner/data3/Jason/Active/climateLD')
addpath('scripts')
%%
control=load('/media/owner/data3/Jason/Active/climateLD/analysis_results/climateLDcontrol_amp_rmv.mat')
witnessed=load('/media/owner/data3/Jason/Active/climateLD/analysis_results/climateLDwitnessed_amp_rmv.mat')
happenedto=load('/media/owner/data3/Jason/Active/climateLD/analysis_results/climateLDhappenedto_amp_rmv.mat')
%%
control=cleanDataset(control);
witnessed=cleanDataset(witnessed);
happenedto=cleanDataset(happenedto);

c=control.EVgain;
w=witnessed.EVgain;
h=happenedto.EVgain;


%% Baseline -500-0


%% 5SD rejection
function combined=rej5SD(c,w,h)
combined=cat(4,c,w,h);
for f=1:size(combined,1)
    for chan=1:size(combined,2)
        for t=1:size(combined,3)
            zval = zscore(squeeze(combined(f,chan,t,:)));
            combined(f,chan,t,zval>5)=nan;
        end
    end
end
end

%% get rid of missing subjects
function group = cleanDataset(group)
    missingidx=find(isnan(squeeze(sum(sum(sum(group.EVgain,1),2),3))));
    group.EVgain(:,:,:,missingidx)=[];
    group.EVloss(:,:,:,missingidx)=[];
    group.GLbias(:,:,:,missingidx)=[];
    group.subjectcoll(missingidx)=[];
end



