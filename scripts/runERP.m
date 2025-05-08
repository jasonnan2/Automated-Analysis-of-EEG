clear all;close all;clc
addpath('scripts')
addpath('scripts/DataAnalysis')
addpath('scripts/DataAnalysis/functions')
eeglab % Required
%% Loading in data and definitions
dataPath = 'sample.mat'; % random sample data
load(dataPath)

% Formatting behavior table for fitlm 
%baseModel=""; % sepcify model for fitlm
%behTbl=readtable(); % read behavior table for neural behavior analysis
behTbl=table();
rng(1234)
behTbl.Subject=project.scalpData.Group1.subList';
behTbl.behVar1=rand(1,10)'

% Defining time in mS for baseline correction
baselineTime=[-250 -50];
% Defining time ranges of interest - replace 'time1' with whatever you want
timeRange=struct();
timeRange.time1=[0 500];
timeRange.time2=[500 1000];
timeRange.time3=[1000 1500];

% Defining optional parameters to test for
vars2plot={'NeuralVarName'};
freq2plot={'theta','alpha','beta'}; 
times2plot={'time1'};
combinations = [1,2];
errorType='sem';
chans2plot={''}; 
FDRflag=0; 
toPlot=1;
isnormal='auto';

%% Scalp Analysis
scalpObject=ScalpObject(project.scalpData, project.info, baselineTime, timeRange); % Initialize scalp processing
scalpObject = scalpObject.cleanDatasets(); % remove any missing subjects if any from Scalp data
scalpObject = scalpObject.standardProcessing(); % standard processing pipeline, 5SD outlier, baseline correction
%%% Transforming data - operates on time series data
%func=@(x) abs(hilbert(x));
%scalpObject = scalpObject.applyFunc(func);
scalpObject = scalpObject.calChanData(); % calcualtes all significant electrodes between groups across all conditions

scalpObject = scalpObject.analyzeScalp('vars2plot',vars2plot,'freq2plot',freq2plot, ...
    'times2plot',times2plot,'combinations',combinations,'FDRflag',FDRflag,'toPlot',toPlot,'isnormal',isnormal);

%% Scalp Plotting ERP and Bar graph
close all
scalpObject.plotERPs('vars2plot',vars2plot,'freq2plot',freq2plot,'times2plot',times2plot)
scalpObject.plotScalpBar('vars2plot',vars2plot,'freq2plot',freq2plot,'times2plot',times2plot)

%% Scalp Behavior Models

scalpObject = scalpObject.calSigTbl('sig'); % creating table of significant neural attributes saved in scalpObject.scalpresults.sigValues.(property).(time)
% Modify Variables for neural behavioral models
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
scalpObject.neuralBehMdl.modelName_NeuralVarName

