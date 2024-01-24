%% Climate luckydoor Main analysis - with subclasses
% Jason Nan
% 12/19/2023 - 
clear all;close all;clc
cd('C:\Users\jason\OneDrive\Desktop\MS Beng\Neatlabs\Climate_LuckyDoor')
addpath('scripts')
addpath('scripts/DataAnalysis')
addpath('scripts/DataAnalysis/functions')
addpath('A:\eeglab_OldLSL_DataAna04072023')
eeglab
%% Loading in data and definitions
dataPath='A:/ClimateLD/analysis_results/formattedData/ClimateLD_allgroups_amp_rmv.mat';
%dataPath='A:/ClimateLD/analysis_results/formattedData/ClimateLD_combinedGroup_amp_rmv.mat';
%dataPath='A:\ClimateLD\analysis_results\formattedData\ClimateLD_traumaGroup_amp_rmv.mat';
%dataPath='A:/ClimateLD/analysis_results/formattedData/ClimateLD_combinedGroup_exRareG_amp_rmv.mat';
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
scalpObject=ScalpObject(CLIMATELD.scalpData, CLIMATELD.info, baselineTime, timeRange); % Initialize scalp processing
scalpObject = scalpObject.cleanDatasets(); % remove any missing subjects if any from Scalp data
scalpObject = scalpObject.standardProcessing(); % standard processing pipeline, 5SD outlier, baseline correction
scalpObject = scalpObject.ScalpAnalysis(); % calcualtes all significant electrodes between groups across all conditions
%func=@(x) abs(hilbert(x));
%scalpObject = scalpObject.applyFunc(func);
scalpObject = scalpObject.calSigTbl(); % creating table of significant neural attributes

vars2plot={'EVgain'};
freq2plot={'broadband'};
times2plot={'choice'};

scalpObject = scalpObject.plotScalpMap('vars2plot',vars2plot,'freq2plot',freq2plot,'times2plot',times2plot,'combinations',[1,3]);% scalp topo plots for specified measure

%%% adding in two more models into sigValues table to run NeurBehMdls
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
neuralVar={'EVgain'};
baseModel="ex_WS_RareG ~ 1  + age+ group+"; keyColumnName='Subject';
scalpObject=scalpObject.NeurBehMdl(neuralVar,behTbl,keyColumnName,baseModel,'ex_WS_RareG_linear');
baseModel="ex_WS_RareG ~ 1  + age+ group*"; keyColumnName='Subject';
scalpObject=scalpObject.NeurBehMdl(neuralVar,behTbl,keyColumnName,baseModel,'ex_WS_RareG_combo');

scalpObject.neuralBehMdl.ex_WS_RareG_linear = extractp(scalpObject.neuralBehMdl.ex_WS_RareG_linear_EVgain,'',1);
scalpObject.neuralBehMdl.ex_WS_RareG_combo = extractp(scalpObject.neuralBehMdl.ex_WS_RareG_combo_EVgain,"group_Control:",0);
%%
close all
errorType='sem';
% Plotting significant channels from topo plots as line plots and bargraph 
% Only selecting conditions which are specified above
scalpObject.plotERPs('vars2plot',vars2plot,'freq2plot',freq2plot,'times2plot',times2plot)
scalpObject.plotScalpBar('vars2plot',vars2plot,'freq2plot',freq2plot,'times2plot',times2plot)
%% Source Analysis
% repeat general processing steps
sourceObject=SourceObject(CLIMATELD.sourceData, CLIMATELD.info, baselineTime, timeRange);
sourceObject = sourceObject.cleanDatasets(); 
sourceObject = sourceObject.standardProcessing();
sourceObject = sourceObject.SourceAnalysis();
sourceObject = sourceObject.calSigTbl();

%sourceObject.plotNetwork(netwrk,vars2plot); % grouped network bar plots for specified measure
sourceObject = sourceObject.plotBrainmap('vars2plot',vars2plot,'freq2plot',freq2plot,'times2plot',times2plot,'combinations',[1,3]); % full roi plot for specified measure

neuralVar={'EVgain'};
baseModel="ex_WS_RareG ~ 1 + age +group+"; keyColumnName='Subject';
sourceObject=sourceObject.NeurBehMdl(neuralVar,behTbl,keyColumnName,baseModel,'ex_WS_RareG_linear');

baseModel="ex_WS_RareG ~ 1 + age +group*"; keyColumnName='Subject';
sourceObject=sourceObject.NeurBehMdl(neuralVar,behTbl,keyColumnName,baseModel,'ex_WS_RareG_combo');

%fdr(sourceObject.sourceResults.sigROIsP.exRareg.all.alpha(:,end))
sourceObject.neuralBehMdl.ex_WS_RareG_linear = extractp(sourceObject.neuralBehMdl.ex_WS_RareG_linear,'',0);
sourceObject.neuralBehMdl.ex_WS_RareG_combo = extractp(sourceObject.neuralBehMdl.ex_WS_RareG_combo,"group_Control:");
sourceObject.neuralBehMdl.ex_WS_RareG_combo = extractp(sourceObject.neuralBehMdl.ex_WS_RareG_combo,"group_HappenedTo:");

%% 
sourceChoicetbl = sourceObject.sigValues.exRareg.choice;
source_mdl=sourceObject.neuralBehMdl.ex_WS_RareG_linear.model{1};

scalpChoicetbl = scalpObject.sigValues.exRareg.choice;
scalp_mdl=scalpObject.neuralBehMdl.ex_WS_RareG_combo.model{1};

%%
close all
aoctool(sourceChoicetbl.exRareg_choice_alpha_posteriorcingulateR,sourceChoicetbl.ex_WS_RareG,sourceChoicetbl.group,0.05,'pccR','ex_WS_rareG')
aoctool(scalpChoicetbl.exRareg_choice_alpha_Pz,scalpChoicetbl.ex_WS_RareG,scalpChoicetbl.group,0.05,'Pz','ex_WS_rareG')

%% Final Figures
close all
var1=scalpChoicetbl.exRareg_choice_alpha_Pz;
var2=scalpChoicetbl.ex_WS_RareG;
P=polyfit(var1,var2,1);
Bfit=polyval(P,var1);
hold on
plot(var1,var2,'b.')
plot(var1,Bfit,'r-')
hold off
axis square
title('ex WS RareG','fontweight','bold')
xlabel('Average alpha Pz during Choice period')
[r,p]=corr(var1,var2,'Type','Spearman','Rows','Complete')

%%
sourceChoicetbl;
mdlCon=fitlm(scalpChoicetbl(strcmp(sourceChoicetbl.group,'Control'),:),"ex_WS_RareG ~ 1 + age +exRareg_choice_alpha_Pz",'RobustOpts','on');
mdlHap=fitlm(scalpChoicetbl(strcmp(sourceChoicetbl.group,'HappenedTo'),:),"ex_WS_RareG ~ 1 + age +exRareg_choice_alpha_Pz",'RobustOpts','on');
close all
h1 = plot(mdlCon); % h1 contains the lines for model 1
set(h1(1),'Color','r'); % make data line for model 1 green
set(h1(2:end),'Color','r'); % make fit and confidence lines for model 1 cyan
hold on
h2 = plot(mdlHap); % h2 contains the lines for model 2
set(h2(1),'Color','b'); % make data line for model 2 orange
set(h2(2:end),'Color','b'); % make fit and confidence lines for model 2 magenta
legend([h1(1); h2(1)],{'Control','Happenedto'});
title('scalp Pz fitlm for groups seperatly')
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
