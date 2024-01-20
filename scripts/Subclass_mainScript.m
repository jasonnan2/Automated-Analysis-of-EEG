%% Climate luckydoor Main analysis - with subclasses
% Jason Nan
% 12/19/2023 - 
clear all;close all;clc
cd('C:\Users\jason\OneDrive\Desktop\MS Beng\Neatlabs\Climate_LuckyDoor')
addpath('scripts')
addpath('scripts/DataAnalysis')
addpath('A:\eeglab_OldLSL_DataAna04072023')
eeglab
%% Loading in data and definitions
%dataPath='A:/ClimateLD/analysis_results/formattedData/ClimateLD_allgroups_amp_rmv.mat';
%dataPath='A:/ClimateLD/analysis_results/formattedData/ClimateLD_combinedGroup_amp_rmv.mat';
%dataPath='A:\ClimateLD\analysis_results\formattedData\ClimateLD_traumaGroup_amp_rmv.mat';
dataPath='A:/ClimateLD/analysis_results/formattedData/ClimateLD_combinedGroup_exRareG_amp_rmv.mat';
load(dataPath)
CLIMATELD.info.experimentalDesign='twoSample';

% Formatting behavior table for fitlm 
baseModel="ex_WS_RareG ~ 1 + group + age +";
tbl=readtable("A:\ClimateLD\ClimateStudy_LD_Beh2.xlsx",'sheet','All');
behTbl = tbl(:, {'Subject','age','ex_WS_RareG'});
behTbl.Subject = cellfun(@(x) erase(x, ['_' extractAfter(x, '_')]), behTbl.Subject, 'UniformOutput', false);
behTbl.Subject = cellfun(@(x) replace(x, ' ',''), behTbl.Subject, 'UniformOutput', false);

% Defining time in mS for baseline correction
baselineTime=[-250 -50];
% Defining time ranges of interest
timeRange=struct();
timeRange.choice=[0 500];
timeRange.imReward=[500 1000];
timeRange.cumReward=[1000 1500];
%timeRange.all=[0 1500];

% defining network
fpn=[5 6 55 66 59 60];netwrk(1).name='FPN';netwrk(1).roi=fpn;
con=[3 4 19 20 37 38 39 40 41 42 57 58 67 68];netwrk(2).name='CON';netwrk(2).roi=con;
admn=[11 12 25 26 29 30 53 54];netwrk(3).name='aDMN';netwrk(3).roi=admn;
pdmn=[15 16 21 22 51 52];netwrk(4).name='pDMN';netwrk(4).roi=pdmn;
mtldmn=[9 10 17 18 31 32 35 36 61 62 65 66 1];netwrk(5).name='mtlDMN';netwrk(5).roi=mtldmn;
vis=[7 8 13 14 23 24 27 28 43 44];netwrk(6).name='Visual';netwrk(6).roi=vis;
sm=[33 34 49 50 45 46];netwrk(7).name='SM';netwrk(7).roi=sm;
van=[2 47 48 63 64];netwrk(8).name='VAN';netwrk(8).roi=van;

%% Scalp Analysis
scalpObject=ScalpAnalysis(CLIMATELD.scalpData, CLIMATELD.info, baselineTime, timeRange); % Initialize scalp processing
scalpObject = scalpObject.cleanDatasets(); % remove any missing subjects if any from Scalp data
scalpObject = scalpObject.standardPipeline(); % standard processing pipeline, 5SD outlier, baseline correction
%func=@(x) abs(hilbert(x));
%scalpObject = scalpObject.applyFunc(func);
vars2plot={'exRareg'};
scalpObject = scalpObject.plotScalpmap(vars2plot);% scalp topo plots for specified measure
scalpObject = scalpObject.calSigTbl(); % creating table of significant neural attributes

%%% adding in two more models
time_list=fieldnames(scalpObject.info.timeRange);
for t=1:length(time_list)
    tbldata = scalpObject.scalpResults.sigValues.exRareg.(time_list{t}); % get table data
    %varNames = setdiff(tbldata.Properties.VariableNames,{'age','ex_WS_RareG','subID','group'});
    cols = contains(tbldata.Properties.VariableNames, 'alpha') & contains(tbldata.Properties.VariableNames, 'F');
    tbldata.("exRareg_"+time_list{t}+"_alpha_frontal") = mean(tbldata{:, cols},2);
    cols = contains(tbldata.Properties.VariableNames, 'alpha') & contains(tbldata.Properties.VariableNames, 'P');
    tbldata.("exRareg_"+time_list{t}+"_alpha_parietal") = mean(tbldata{:, cols},2);
    scalpObject.scalpResults.sigValues.exRareg.(time_list{t})=tbldata;
end
%fdr(scalpObject.scalpResults.sigElectrodesP.exRareg.all.alpha)
% behavior and neural analysis
baseModel="ex_WS_RareG ~ 1  + age+ group+"; keyColumnName='Subject';
scalpObject=scalpObject.NeurBehMdl(behTbl,keyColumnName,baseModel,'ex_WS_RareG_linear');
baseModel="ex_WS_RareG ~ 1  + age+ group*"; keyColumnName='Subject';
scalpObject=scalpObject.NeurBehMdl(behTbl,keyColumnName,baseModel,'ex_WS_RareG_combo');

scalpObject.scalpResults.neuralBehMdl.ex_WS_RareG_linear = extractp(scalpObject.scalpResults.neuralBehMdl.ex_WS_RareG_linear,'',0);
scalpObject.scalpResults.neuralBehMdl.ex_WS_RareG_combo = extractp(scalpObject.scalpResults.neuralBehMdl.ex_WS_RareG_combo,"group_Control:",0);
%%
freq2plot={'alpha','broadband'};
times2plot={'imReward','cumReward','choice'};
close all
errorType='95CI';
% Plotting significant channels from topo plots as line plots and bargraph 
% Only selecting conditions which are specified above
scalpObject.plotERPs(vars2plot,freq2plot,times2plot,errorType)
scalpObject.plotScalpBar(vars2plot,freq2plot,times2plot)
%% brainNetwork Analysis
% repeat general processing steps
sourceObject=SourceAnalysis(CLIMATELD.sourceData, CLIMATELD.info, baselineTime, timeRange);
sourceObject = sourceObject.cleanDatasets(); 
sourceObject = sourceObject.standardPipeline();
sourceObject.plotNetwork(netwrk,vars2plot); % grouped network bar plots for specified measure
%combinations=[1 3;2 3];
combinations=[1,2];
sourceObject = sourceObject.plotBrainmap(vars2plot,combinations); % full roi plot for specified measure
sourceObject = sourceObject.calSigTbl();
baseModel="ex_WS_RareG ~ 1 + age +group+"; keyColumnName='Subject';
sourceObject=sourceObject.NeurBehMdl(behTbl,keyColumnName,baseModel,'ex_WS_RareG_linear')
baseModel="ex_WS_RareG ~ 1 + age +group*"; keyColumnName='Subject';
sourceObject=sourceObject.NeurBehMdl(behTbl,keyColumnName,baseModel,'ex_WS_RareG_combo')
%fdr(sourceObject.sourceResults.sigROIsP.exRareg.all.alpha(:,end))
sourceObject.sourceResults.neuralBehMdl.ex_WS_RareG_linear = extractp(sourceObject.sourceResults.neuralBehMdl.ex_WS_RareG_linear,'',0)
sourceObject.sourceResults.neuralBehMdl.ex_WS_RareG_combo = extractp(sourceObject.sourceResults.neuralBehMdl.ex_WS_RareG_combo,"group_Control:",0)

%%

sourceChoicetbl = sourceObject.sourceResults.sigValues.exRareg.choice;
scalpChoicetbl = scalpObject.scalpResults.sigValues.exRareg.choice;

%%

figure
scatterplot(





















%% Custom Neural Beh models
% 
% varNames={'exRareg_all_alpha_frontal','exRareg_all_alpha_parietal'};
% baseModel="ex_WS_RareG ~ 1 +  age +group+"; keyColumnName='Subject';
% 
% neural_beh_tbl=table();
% pvals=[];
% for i=1:length(varNames)
%     modelDef =  baseModel+varNames{i};
%     model = fitlm(tbldata,modelDef,'RobustOpts','on')
%     pvals(i) = model.Coefficients.pValue(end); % Check here to make sure its always the second
% end
% neural_beh_tbl.(baseModel)=varNames';
% neural_beh_tbl.p=pvals';
% neural_beh_tbl(neural_beh_tbl.p>0.05,:)=[];
% 
