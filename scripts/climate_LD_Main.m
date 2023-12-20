%% Climate luckydoor Main analysis
% Jason Nan
% 12/19/2023
clear all;close all;clc
cd('C:\Users\jason\OneDrive\Desktop\MS Beng\Neatlabs\Climate_LuckyDoor')
addpath('scripts')

%%
% Create an object of the class
groupPaths = {'A:/ClimateLD/analysis_results/climateLDcontrol_amp_rmv.mat', 'A:/ClimateLD/analysis_results/climateLDwitnessed_amp_rmv.mat','A:/ClimateLD/analysis_results/climateLDhappenedto_amp_rmv.mat'};
variables = {'EVgain', 'EVloss', 'GLbias'};
baselineIDX=[1,126];


analysis = ERPanalysis(groupPaths,variables,baselineIDX);
analysis = analysis.cleanDatasets();
finishedAnalysis = analysis.standardPipeline();





