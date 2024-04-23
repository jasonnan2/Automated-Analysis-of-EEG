clear all;close all;clc
addpath('scripts')
addpath('scripts/DataAnalysis')
addpath('scripts/DataAnalysis/functions')
%% Loading in data and definitions
dataPath = 'sample.mat'; % random Data
load(dataPath)

% Formatting behavior table for fitlm 
%baseModel=""; % sepcify model for fitlm
%behTbl=readtable(); % read behavior table for neural behavior analysis
behTbl=table();
behTbl.Subject=project.scalpData.Group1.subList';
behTbl.behVar1=rand(1,10)';

% Defining time in mS for baseline correction
baselineTime=[-250 -50];
% Defining time ranges of interest - replace 'time1' with whatever you want
timeRange=struct();
timeRange.time1=[0 500];
timeRange.time2=[500 1000];
timeRange.time3=[1000 1500];

% defining network for source localization - https://github.com/aojeda/dsi 
fpn=[5 6 55 66 59 60];netwrk(1).name='FPN';netwrk(1).roi=fpn;
con=[3 4 19 20 37 38 39 40 41 42 57 58 67 68];netwrk(2).name='CON';netwrk(2).roi=con;
admn=[11 12 25 26 29 30 53 54];netwrk(3).name='aDMN';netwrk(3).roi=admn;
pdmn=[15 16 21 22 51 52];netwrk(4).name='pDMN';netwrk(4).roi=pdmn;
mtldmn=[9 10 17 18 31 32 35 36 61 62 65 66 1];netwrk(5).name='mtlDMN';netwrk(5).roi=mtldmn;
vis=[7 8 13 14 23 24 27 28 43 44];netwrk(6).name='Visual';netwrk(6).roi=vis;
sm=[33 34 49 50 45 46];netwrk(7).name='SM';netwrk(7).roi=sm;
van=[2 47 48 63 64];netwrk(8).name='VAN';netwrk(8).roi=van;

%% Scalp Analysis
scalpObject=ScalpObject(project.scalpData, project.info, baselineTime, timeRange); % Initialize scalp processing
scalpObject = scalpObject.cleanDatasets(); % remove any missing subjects if any from Scalp data
scalpObject = scalpObject.standardProcessing(); % standard processing pipeline, 5SD outlier, baseline correction
scalpObject = scalpObject.ScalpAnalysis(); % calcualtes all significant electrodes between groups across all conditions
%%% Transforming data 
%func=@(x) abs(hilbert(x));
%scalpObject = scalpObject.applyFunc(func);
scalpObject = scalpObject.calSigTbl('sig'); % creating table of significant neural attributes

%% Scalp Plotting
close all
vars2plot={'NeuralVarName'}; freq2plot={'theta','alpha','beta'}; times2plot={'time1'}; errorType='sem'; chans2plot={'Pz'};
scalpObject.plotScalpMap('vars2plot',vars2plot,'freq2plot',freq2plot,'times2plot',times2plot,'combinations',[1,2]);
scalpObject.plotERPs('vars2plot',vars2plot,'freq2plot',freq2plot,'times2plot',times2plot)
scalpObject.plotScalpBar('vars2plot',vars2plot,'freq2plot',freq2plot,'times2plot',times2plot,'chans2plot',chans2plot)

%% Scalp Behavior Models
% scalpObject = scalpObject.calSigTbl('sig'); % creating table of significant neural attributes
% time_list=fieldnames(scalpObject.info.timeRange);
% for t=1:length(time_list)
%     tbldata = scalpObject.scalpResults.sigValues.NeuralVarName.(time_list{t}); % get table data
%     [~,tbldata] = calculate_group_means(tbldata, groups, append(time_list(t),'_',groupnames)); % only use newtable. If testing all, then use first output
%     scalpObject.scalpResults.sigValues.exRareg.(time_list{t}) = tbldata;
% end

% behavior and neural analysis
neuralVar={'NeuralVarName'}; % Specify here
baseModel="behVar1 ~ 1+ "; % define model for fitlm - automatically appends neuralVar
keyColumnName='Subject'; % column name with subjectIds
scalpObject=scalpObject.NeurBehMdl(neuralVar,behTbl,keyColumnName,baseModel,'modelName');
scalpObject.neuralBehMdl.modelName_NeuralVarName = extractp(scalpObject.neuralBehMdl.modelName_NeuralVarName,"NeuralVarName",1); % extract all p-values from fitlm


%% Source Analysis
% repeat general processing steps
sourceObject=SourceObject(project.sourceData, project.info, baselineTime, timeRange);
sourceObject = sourceObject.cleanDatasets(); 
sourceObject = sourceObject.standardProcessing();
sourceObject = sourceObject.SourceAnalysis();
sourceObject = sourceObject.calSigTbl();

%% Source Plotting
sourceObject.plotNetwork(netwrk,vars2plot); % grouped network bar plots for specified measure
sourceObject.plotBrainmap('vars2plot',vars2plot,'freq2plot',freq2plot,'times2plot',times2plot,'combinations',[1,2]); % full roi plot for specified measure

%% Source Behavior Model - same format as scalp
neuralVar={'NeuralVarName'};
baseModel="behVar1 ~ 1+ "; % define model for fitlm - automatically appends neuralVar
keyColumnName='Subject';
sourceObject=sourceObject.NeurBehMdl(neuralVar,behTbl,keyColumnName,baseModel,'modelName');
sourceObject.neuralBehMdl.modelName_NeuralVarName = extractp(sourceObject.neuralBehMdl.modelName_NeuralVarName,"NeuralVarName",0);
