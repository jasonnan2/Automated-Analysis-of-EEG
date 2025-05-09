%% Source Analysis
clear all;close all;clc
addpath('scripts')
addpath('scripts/DataAnalysis')
addpath('scripts/DataAnalysis/functions')
eeglab % Required

%% Loading in data and definitions
dataPath = 'sample.mat'; % random sample data
load(dataPath)

% Defining time in mS for baseline correction
baselineTime=[-250 -50];
% Defining time ranges of interest - replace 'time1' with whatever you want
timeRange=struct();
timeRange.time1=[0 500];
timeRange.time2=[500 1000];
timeRange.time3=[1000 1500];

% Defining optional parameters to test for
cfg = struct();
cfg.vars2plot={'NeuralVarName'};
cfg.freq2plot={'theta','alpha','beta'}; 
cfg.times2plot={'time1'};
cfg.groups2plot=[1:3]; % for network plots
cfg.combinations = [1,2; 1,3];
cfg.rois2plot = [42,43,47,48]; 
cfg.FDRflag=1; 
cfg.toPlot=1;
cfg.isnormal='auto';
cfg.hmFile = 'headModel_templateFile_newPEBplus.mat';
cfg.netConMethod = 'net';

% defining network for source localization - https://github.com/aojeda/dsi 
fpn=[5 6 55 66 59 60];netwrk(1).name='FPN';netwrk(1).roi=fpn;
con=[3 4 19 20 37 38 39 40 41 42 57 58 67 68];netwrk(2).name='CON';netwrk(2).roi=con;
admn=[11 12 25 26 29 30 53 54];netwrk(3).name='aDMN';netwrk(3).roi=admn;
pdmn=[15 16 21 22 51 52];netwrk(4).name='pDMN';netwrk(4).roi=pdmn;
mtldmn=[9 10 17 18 31 32 35 36 61 62 65 66 1];netwrk(5).name='mtlDMN';netwrk(5).roi=mtldmn;
vis=[7 8 13 14 23 24 27 28 43 44];netwrk(6).name='Visual';netwrk(6).roi=vis;
sm=[33 34 49 50 45 46];netwrk(7).name='SM';netwrk(7).roi=sm;
van=[2 47 48 63 64];netwrk(8).name='VAN';netwrk(8).roi=van;

%% repeat general processing steps
sourceObject=SourceObject(project.sourceData, project.info, baselineTime, timeRange,cfg);
sourceObject = sourceObject.cleanDatasets(); 
sourceObject = sourceObject.standardProcessing();
sourceObject = sourceObject.calRoiData();
% sourceObject = sourceObject.analyzeRoi('FDRflag',0);
% sourceObject = sourceObject.calSigTbl();
sourceObject = sourceObject.calNetData(netwrk);
sourceObject = sourceObject.calConnectivity(netwrk);

%% Source Plotting
close all
sourceObject.plotNetwork(); % grouped network bar plots for specified measure
sourceObject.plotNetConnectivity('FDRflag',1)

%% Source Behavior Model - same format as scalp
% Formatting behavior table for fitlm 
behTbl=table();
rng(1234)
behTbl.Subject=project.scalpData.Group1.subList';
behTbl.behVar1=rand(1,10)';

neuralVar={'NeuralVarName'};
baseModel="behVar1 ~ 1+ "; % define model for fitlm - automatically appends neuralVar
keyColumnName='Subject';
sourceObject=sourceObject.NeurBehMdl(neuralVar,behTbl,keyColumnName,baseModel,'modelName');
sourceObject.neuralBehMdl.modelName_NeuralVarName = extractp(sourceObject.neuralBehMdl.modelName_NeuralVarName,"NeuralVarName",0);
