%% Climate luckydoor Main analysis - with subclasses
% Jason Nan
% 12/19/2023
clear all;close all;clc
cd('C:\Users\jason\OneDrive\Desktop\MS Beng\Neatlabs\Climate_LuckyDoor')
addpath('scripts')
addpath('scripts/DataAnalysis')
addpath('A:\eeglab_OldLSL_DataAna04072023')
eeglab
%% Loading in data and general preprocessing
% loading in data
dataPath='A:/ClimateLD/analysis_results/ClimateLD_allgroups_amp_rmv.mat';
load(dataPath)
% Defining time in mS for baseline correction
baselineTime=[-250 -50];
% Defining time ranges of interest
timeRange.choice=[0 500];
timeRange.imReward=[500 1000];
timeRange.cumReward=[1000 1500];

% defining network
fpn=[5 6 55 66 59 60];
con=[3 4 19 20 37 38 39 40 41 42 57 58 67 68];
admn=[11 12 25 26 29 30 53 54];
pdmn=[15 16 21 22 51 52];
mtldmn=[9 10 17 18 31 32 35 36 61 62 65 66 1];
vis=[7 8 13 14 23 24 27 28 43 44];
sm=[33 34 49 50 45 46];
van=[2 47 48 63 64];

netwrk(1).name='FPN';
netwrk(1).roi=fpn;
netwrk(2).name='CON';
netwrk(2).roi=con;
netwrk(3).name='aDMN';
netwrk(3).roi=admn;
netwrk(4).name='pDMN';
netwrk(4).roi=pdmn;
netwrk(5).name='mtlDMN';
netwrk(5).roi=mtldmn;
netwrk(6).name='Visual';
netwrk(6).roi=vis;
netwrk(7).name='SM';
netwrk(7).roi=sm;
netwrk(8).name='VAN';
netwrk(8).roi=van;


%% Scalp ERP Analysis

scalpObject=ScalpAnalysis(CLIMATELD.scalpData, CLIMATELD.info, baselineTime, timeRange);
% remove any missing subjects if any from Scalp data
scalpObject = scalpObject.cleanDatasets(); 
% standard processing pipeline, 5SD outlier, baseline correction
scalpObject = scalpObject.standardPipeline(); 

% Plotting scalp topo plots for variables in vars2plot and each timeRange
% speficied
vars2plot={'EVgain'};
scalpObject = scalpObject.plotScalpmap(vars2plot);
% Plotting significant channels from topo plots as line plots with all
% groups overlayed. Only selecting conditions which are in the specified
% vectors for var, frequency, and time ranges. 
vars2plot={'EVgain'};
freq2plot={'alpha'};
times2plot={'choice','imReward','cumReward'};
scalpObject.plotERPs(vars2plot,freq2plot,times2plot)

%% brainNetwork Analysis
sourceObject=SourceAnalysis(CLIMATELD.sourceData, CLIMATELD.info, baselineTime, timeRange);
sourceObject = sourceObject.cleanDatasets(); 
sourceObject = sourceObject.standardPipeline(); 

sourceObject.plotNetwork(netwrk,vars2plot)
















