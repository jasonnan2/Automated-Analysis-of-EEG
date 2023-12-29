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
vars2plot={'EVgain'};
analysis = analysis.plotScalpmap(timeRange,vars2plot);

vars2plot={'EVgain'};
freq2plot={'theta','alpha','beta'};
times2plot={'choice','imReward','cumReward'};

analysis.plotERPs(vars2plot,freq2plot,times2plot)

    
            