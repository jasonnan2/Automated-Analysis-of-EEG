%% Source Analysis
clear all;close all;clc
addpath('scripts')
addpath('scripts/DataAnalysis')
addpath('scripts/DataAnalysis/functions')
addpath('samplePlots')
eeglab % Required

%% Loading in data and definitions
dataPath = './samplePlots/sample.mat'; % random sample data
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
cfg.vars2plot={'NeuralVarName'};        % Based on study definition
cfg.freq2plot={'theta','alpha','beta'}; % Based on study definition
cfg.times2plot={'time1'};               % Based on study definition
cfg.groups2plot=[1:2];                  % Based on study definition
cfg.combinations = [1,2];          % Based on study definition
cfg.rois2plot = [42,43,47,48];          % Based on ROI order and list
cfg.FDRflag=1;                          % 1/0 to calculate FDR or not
cfg.toPlot=1;                           % For brain plots during analyzeRoi()
cfg.isnormal='auto';                    % 'auto', 1, 0
cfg.hmFile = 'headModel_templateFile_newPEBplus.mat'; 
cfg.netConMethod = 'net';               % 'net' or 'roi'

% defining network for source localization - https://github.com/aojeda/dsi 
fpn=[5 6 55 66 59 60];netwrk(1).name='FPN';netwrk(1).roi=fpn;
con=[3 4 19 20 37 38 39 40 41 42 57 58 67 68];netwrk(2).name='CON';netwrk(2).roi=con;
admn=[11 12 25 26 29 30 53 54];netwrk(3).name='aDMN';netwrk(3).roi=admn;
pdmn=[15 16 21 22 51 52];netwrk(4).name='pDMN';netwrk(4).roi=pdmn;
mtldmn=[9 10 17 18 31 32 35 36 61 62 65 66 1];netwrk(5).name='mtlDMN';netwrk(5).roi=mtldmn;
vis=[7 8 13 14 23 24 27 28 43 44];netwrk(6).name='Visual';netwrk(6).roi=vis;
sm=[33 34 49 50 45 46];netwrk(7).name='SM';netwrk(7).roi=sm;
van=[2 47 48 63 64];netwrk(8).name='VAN';netwrk(8).roi=van;

%% Create SourceObject and run each step of the source analysis pipeline
sourceObject = SourceObject(project.sourceData, project.info, baselineTime, timeRange, cfg);

sourceObject = sourceObject.cleanDatasets();        % Clean or remove bad datasets
sourceObject = sourceObject.standardProcessing();   % Apply baseline correction and artifact rejection
sourceObject = sourceObject.calRoiData();           % Calculates ROI-level time series
sourceObject = sourceObject.analyzeRoi();           % Run statistical tests based on cfg, generate ROI plots
sourceObject = sourceObject.calSigTbl();            % Compile significant findings into table
sourceObject = sourceObject.calNetData(netwrk);     % Compute average signal per network
sourceObject = sourceObject.calConnectivity(netwrk);% Compute connectivity between networks

%% Plot Source Results
close all
sourceObject.plotNetwork();                         % Plot bar graphs of network-level signal
sourceObject.plotNetConnectivity('FDRflag',1);      % Plot network connectivity matrices

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
