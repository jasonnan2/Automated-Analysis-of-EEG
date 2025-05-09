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

% Defining config structure for focused analysis
cfg=struct();
cfg.vars2plot={'NeuralVarName'};
cfg.freq2plot={'theta','alpha','beta'}; 
cfg.times2plot={'time1'};
cfg.combinations = [1,2];
cfg.groups2plot = [1,2];
cfg.errorType='none';
cfg.FDRflag=1; 
cfg.toPlot=1;
cfg.isnormal='auto';

%% Scalp Analysis
close all
scalpObject=ScalpObject(project.scalpData, project.info, baselineTime, timeRange,cfg); % Initialize scalp processing
scalpObject = scalpObject.cleanDatasets(); % remove any missing subjects if any from Scalp data
scalpObject = scalpObject.standardProcessing(); % standard processing pipeline, 5SD outlier, baseline correction
scalpObject = scalpObject.calChanData(); % calcualtes all significant electrodes between groups across all conditions

%% Opional name value args to bypass config settings

close all
scalpObject = scalpObject.analyzeScalp('toPlot',0);
scalpObject.plotERPs('freq2plot',{'theta'})
scalpObject.plotScalpBar('freq2plot',{'theta'},'chans2plot',{'Pz'}); 

%% Scalp Behavior Models
% Formatting behavior table for fitlm 
behTbl=table();
rng(1234)
behTbl.Subject=project.scalpData.Group1.subList';
behTbl.behVar1=rand(1,10)';
scalpObject = scalpObject.calSigTbl('sig'); 

% behavior and neural analysis
neuralVar={'NeuralVarName'}; % Specify here
baseModel="behVar1 ~ 1+ "; % define model for fitlm - automatically appends neuralVar
keyColumnName='Subject'; % column name with subjectIds
scalpObject=scalpObject.NeurBehMdl(neuralVar,behTbl,keyColumnName,baseModel,'modelName');
scalpObject.neuralBehMdl.modelName_NeuralVarName = extractp(scalpObject.neuralBehMdl.modelName_NeuralVarName,"NeuralVarName",1); % extract all p-values from fitlm
scalpObject.neuralBehMdl.modelName_NeuralVarName

%% Extra Functionality and operations
%% Transforming data - operates on time series data
func=@(x) abs(hilbert(x));
scalpObject = scalpObject.applyFunc(func);
%% Modify Variables for neural behavioral models

time_list=fieldnames(scalpObject.info.timeRange);
for t=1:length(time_list)
    tbldata = scalpObject.scalpResults.sigValues.NeuralVarName.(time_list{t}); % get original table data
    [~,tbldata] = calculate_group_means(tbldata, groups, append(time_list(t),'_',groupnames)); % only use newtable. If testing all, then use first output
    scalpObject.scalpResults.sigValues.NeuralVarName.(time_list{t}) = tbldata;
end










