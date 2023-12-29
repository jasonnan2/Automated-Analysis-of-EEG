%% Climate luckydoor Main analysis
% Jason Nan
% 12/19/2023
clear all;close all;clc
cd('C:\Users\jason\OneDrive\Desktop\MS Beng\Neatlabs\Climate_LuckyDoor')
addpath('scripts')
addpath('A:\eeglab_OldLSL_DataAna04072023')
eeglab
%% Main Analysis with Class
% Create an object of the class
dataPath='A:/ClimateLD/analysis_results/ClimateLD_allgroups_amp_rmv.mat';
load(dataPath)
variables = {'EVgain', 'EVloss', 'GLbias'};
baselineTime=[-250 -50];

analysis = ERPanalysis_struct(CLIMATELD); % set up the class, need filepath of group data, variables to run and baseline correction range
analysis = analysis.cleanDatasets(); % remove any missing subjects
analysis = analysis.standardPipeline(baselineTime); % standard processing pipeline, 5SD outlier, baseline correction

%%

% choice 0-500
% immediate reward 500-1000
% cum reward 1000-1500

timeRange.choice=[0 500];
timeRange.imReward=[500 1000];
timeRange.cumReward=[1000 1500];
Vars2Plot={'EVgain'};
analysis = analysis.plotERP(timeRange,Vars2Plot);
%%


topoplot(zeros(1,24),chanlocs_str,'headrad','rim','electrodes','labels'); 
%%

% 'EVgain-alpha-imReward-CPz-C3-F7'
% 'EVgain-alpha-cumReward-CPz'
% 'EVgain-alpha-choice-F7-P7-O1'
% 
% 'EVgain-beta-imReward-C3-FC3-FP2'
% 'EVgain-beta-cumReward-C3-FC3-FP2'
% 'EVgain-beta-choice-C3'
% 
% 'EVgain-theta-imreward-P3'
% 'EVgain-theta-choice-FC3'

EVgain_alpha = struct('freq', 'alpha', 'metric', 'EVgain', 'times', {'imReward', 'cumReward', 'choice'}, 'electrodes', {{'CPz', 'C3', 'F7'}, {'CPz'}, {'F7', 'P7', 'O1'}});
EVgain_beta = struct('freq', 'beta', 'metric', 'EVgain', 'times', {'imReward', 'cumReward', 'choice'}, 'electrodes', {{'C3', 'FC3', 'FP2'}, {'C3', 'FC3', 'FP2'}, {'C3'}});
EVgain_theta = struct('freq', 'theta', 'metric', 'EVgain', 'times', {'imReward', 'choice'}, 'electrodes', {{'P3'}, {'FC3'}});

%% Plotting Line plots for EVgain

obj=analysis;
p=1;
N=numel(fieldnames(obj.DATA));
property=Vars2Plot{p};
timeNames = fieldnames(timeRange);
N = numel(fieldnames(obj.DATA));
combinations = nchoosek(1:N,2); % This will give you a matrix where each row is a combination of two groups

for p=1:length(properties)
    property=properties{p};
    
    for n=1:N
        data(:,:,:,n) = nanmean(obj.DATA.(obj.info.groupNames{n}).(property),4);
    end

    figure
    figCount=0;
    for t=1:length(timeNames)
        for f=1:length(obj.info.freq_list)
            figCount = figCount + 1;
            subplot(3,4,figCount)

            [~, startIdx] = min(abs(obj.info.timeAxis - timeRange.(timeNames{t})(1)));
            [~, endIdx] = min(abs(obj.info.timeAxis - timeRange.(timeNames{t})(2)));

            s1=squeeze(nanmean(group1(f,:,startIdx:endIdx,:),3)); % get to chan x sub
            s2=squeeze(nanmean(group2(f,:,startIdx:endIdx,:),3)); % get to chan x sub
            plotSigTopo(s1,s2,obj.info.chanlocs,{'.','+'})
            
            
        end
    end
end

            
            
%%
varSets = {
'EVgain-alpha-imReward-CPz-C3-F7',
'EVgain-alpha-cumReward-CPz',
'EVgain-alpha-choice-F7-P7-O1',
'EVgain-beta-imReward-C3-FC3-FP2',
'EVgain-beta-cumReward-C3-FC3-FP2',
'EVgain-beta-choice-C3',
'EVgain-theta-imreward-P3',
'EVgain-theta-choice-FC3'}

varString=varSets{1}
[metric,freqIdx,timeIdx,elecIdxs]=splitConditions(obj,varString)
%%
function [metric,freqIdx,timeIdx,elecIdxs]=splitConditions(obj,varString)
    % Loop over each variable set
    % Split the variable set into its components
    components = split(varString, '-');
    metric = components{1};
    freq = components{2};
    time = components{3};
    elecs = components(4:end);

    % Find the indices of the components in the dataset
    freqIdx = find(strcmp(obj.info.freq_list, freq));
    timeIdx = obj.info.timeIDX.(time);
    elecIdxs = find(ismember({obj.info.chanlocs.labels}, elecs));
end  
            
            
            
            
            
            
            
            
            
            
            
            
            
            