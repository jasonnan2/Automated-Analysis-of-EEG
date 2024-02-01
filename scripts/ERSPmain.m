clear all;close all;clc
cd('C:\Users\jason\OneDrive\Desktop\MS Beng\Neatlabs\Climate_LuckyDoor')
addpath('scripts')
addpath('scripts/DataAnalysis')
addpath('scripts/DataAnalysis/functions')
load('C:\Users\jason\OneDrive\Desktop\MS Beng\Neatlabs\Climate_LuckyDoor\scripts\New Folder\subjectcoll.mat')
%% Loading in data and definitions
dataPath='A:/ClimateLD/analysis_results/formattedData/ERSP_sample.mat'
load(dataPath)

%% ERSP Analysis
baselineTime=[-750 -550];
% Defining time ranges of interest
timeRange=struct();
timeRange.all=[-948 1704];
timeRange.short=[-200 600]

freqRanges.theta=[3 8];
freqRanges.alpha=[8 13];
freqRanges.beta=[13,30];

erspObject=ERSPObject(ESRP.DATA, ESRP.info, baselineTime, timeRange); % Initialize scalp processing
erspObject = erspObject.cleanDatasets(); % remove any missing subjects if any from Scalp data
erspObject = erspObject.standardProcessing('outlier','3MAD'); % standard processing pipeline, 5SD outlier, baseline correction
scalpObj = erspObject.createERP(freqRanges); % convert ersp data to ERP data by averaging freq. 

erspObject.plotERSP('chans2plot',{'Oz'},'times2plot',{'short'})

%%
timeRange=struct();
timeRange.t=[244 577];
timeRange.a=[211 550];
timeRange.b=[314 592];

scalpObject=ScalpObject(scalpObj.DATA, scalpObj.info, baselineTime, timeRange); % Initialize scalp processing
scalpObject = scalpObject.standardProcessing('outlier','none','calBaseline',0); % standard processing pipeline, 5SD outlier, baseline correction
scalpObject = scalpObject.ScalpAnalysis(); % calcualtes all significant electrodes between groups across all conditions
scalpObject = scalpObject.calSigTbl(); % creating table of significant neural attributes

% Scalp Plotting
close all
vars2plot={'emotion'}; errorType='sem';
scalpObject.plotScalpMap('vars2plot',vars2plot,'freq2plot',{'alpha'},'times2plot',{'a'});
caxis([-5 5])
%scalpObject.plotScalpMap('vars2plot',vars2plot,'freq2plot',{'theta'},'times2plot',{'theta'});
%scalpObject.plotScalpMap('vars2plot',vars2plot,'freq2plot',{'beta'},'times2plot',{'beta'});
%scalpObject.plotERPs('vars2plot',vars2plot,'freq2plot',freq2plot,'times2plot',times2plot)
%scalpObject.plotScalpBar('vars2plot',vars2plot,'freq2plot',freq2plot,'times2plot',times2plot)
