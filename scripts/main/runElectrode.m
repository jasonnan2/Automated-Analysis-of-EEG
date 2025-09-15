clear all;close all;clc
addpath('scripts')
addpath('scripts/DataAnalysis')
addpath('scripts/DataAnalysis/functions')
addpath('samplePlots')
% addpath(genpath('Automated-Analysis-of-EEG'))
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

% Defining config structure for focused analysis
cfg=struct();
cfg.vars2plot={'NeuralVarName'};        % Based on study definition
cfg.freq2plot={'theta','alpha','beta'}; % Based on study definition
cfg.times2plot={'time1'};      % Based on study definition
cfg.combinations = [1,2];               % Based on study definition
cfg.groups2plot = [1,2];               % Based on study definition
cfg.chans2plot = 'all';                 % 'all' - all significant or cell array of any electrodes to plot
cfg.errorType='none';                   % 'none', 'sem','95CI'
cfg.FDRflag=1;                          % 1/0 to calculate FDR or not
cfg.toPlot=1;                           % For electrode topo plots during analyzeElectrode()
cfg.isnormal='auto';                    % 'auto', 1, 0

%% Electrode Analysis
close all
electrodeObject=ElectrodeObject(project.electrodeData, project.info, baselineTime, timeRange,cfg); % Initialize electrode processing
electrodeObject = electrodeObject.cleanDatasets(); % remove any missing subjects if any from Electrode data
electrodeObject = electrodeObject.standardProcessing(); % standard processing pipeline, 5SD outlier, baseline correction
electrodeObject = electrodeObject.calChanData(); % calcualtes all significant electrodes between groups across all conditions

%% Plotting Electrode Analysis
close all
electrodeObject = electrodeObject.analyzeElectrode(); % Calculates electrode differences based on cfg file and plots electrode topo plots
electrodeObject.plotERPs() % plotting ERP after analyze Electrode. Must run analyzeElectrode() before plotting
electrodeObject.plotElectrodeBar(); % plot grouped bar plots after analyze Electrode. Must run analyzeElectrode() before plotting

%% Electrode Behavior Models
% Formatting behavior table for fitlm 
behTbl=table();
rng(1234)
behTbl.Subject=project.electrodeData.Group1.subList';
behTbl.behVar1=rand(1,10)';
electrodeObject = electrodeObject.calSigTbl('sig'); 

% behavior and neural analysis
neuralVar={'NeuralVarName'}; % Specify here
baseModel="behVar1 ~ 1+ "; % define model for fitlm - automatically appends neuralVar
keyColumnName='Subject'; % column name with subjectIds
electrodeObject=electrodeObject.NeurBehMdl(neuralVar,behTbl,keyColumnName,baseModel,'modelName');
electrodeObject.neuralBehMdl.modelName_NeuralVarName = extractp(electrodeObject.neuralBehMdl.modelName_NeuralVarName,"NeuralVarName",1); % extract all p-values from fitlm
electrodeObject.neuralBehMdl.modelName_NeuralVarName

%% Extra Functionality and operations
%% Transforming data - operates on time series data
func=@(x) abs(hilbert(x));
electrodeObject = electrodeObject.applyFunc(func);
%% Modify Variables for neural behavioral models

time_list=fieldnames(electrodeObject.info.timeRange);
for t=1:length(time_list)
    tbldata = electrodeObject.electrodeResults.sigValues.NeuralVarName.(time_list{t}); % get original table data
    [~,tbldata] = calculate_group_means(tbldata, groups, append(time_list(t),'_',groupnames)); % only use newtable. If testing all, then use first output
    electrodeObject.electrodeResults.sigValues.NeuralVarName.(time_list{t}) = tbldata;
end