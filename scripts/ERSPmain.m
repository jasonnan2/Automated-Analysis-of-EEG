clear all;close all;clc
cd('C:\Users\jason\OneDrive\Desktop\MS Beng\Neatlabs\Climate_LuckyDoor')
addpath('scripts')
addpath('scripts/DataAnalysis')
addpath('scripts/DataAnalysis/functions')
addpath 'A:\eeglab_OldLSL_DataAna04072023'
eeglab
%% Loading in data and definitions
dataPath='A:\ClimateLD\analysis_results\formattedData\ClimateLD_ersp_allgroups_exRareG_amp_rmv.mat';
load(dataPath)
CLIMATELD.info.experimentalDesign='twoSample';
%% ERSP Analysis
baselineTime=[-250 -50];
% Defining time ranges of interest
timeRange=struct();
timeRange.base=[-250 500];
timeRange.choice=[0 500];
freqRanges.theta=[3 8];
freqRanges.alpha=[8 13];
freqRanges.beta=[13,30];

erspObject=ERSPObject(CLIMATELD.erspData, CLIMATELD.info, baselineTime, timeRange); % Initialize scalp processing
erspObject = erspObject.cleanDatasets(); % remove any missing subjects if any from Scalp data
erspObject = erspObject.standardProcessing('outlier','3MAD'); % standard processing pipeline, 5SD outlier, baseline correction
scalpObj = erspObject.createERP(freqRanges); % convert ersp data to ERP data by averaging freq. 

%%
close all
erspObject.plotERSPgroups('chans2plot',{'Fz','Pz'},'groups2plot',{'HappenedTo','Witnessed','Control'},'cRange',[-0.5 1])

han=axes(gcf,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylab=ylabel(han,'Freq (Hz)','fontweight','bold','fontsize',30);
xlab=xlabel(han,'Time (mSec)','fontweight','bold','fontsize',16);

ylab.Position(1) =ylab.Position(1)-0.04; % change horizontal position of ylabel
ylab.Position(2) = ylab.Position(2); 

xlab.Position(1) =xlab.Position(1); % change horizontal position of xlabel
xlab.Position(2) = xlab.Position(2)-0.01; 
sgtitle('')
%%
timeRange=struct();
timeRange.choice=[0 500];

scalpObject=ScalpObject(scalpObj.DATA, scalpObj.info, baselineTime, timeRange); % Initialize scalp processing
%scalpObject = scalpObject.standardProcessing('outlier','none','calBaseline',0); % standard processing pipeline, 5SD outlier, baseline correction
scalpObject = scalpObject.ScalpAnalysis(); % calcualtes all significant electrodes between groups across all conditions
scalpObject = scalpObject.calSigTbl(); % creating table of significant neural attributes

% Scalp Plotting
close all
errorType='sem';
scalpObject.plotScalpMap('freq2plot',{'alpha'});
%caxis([-5 5])
%scalpObject.plotScalpMap('vars2plot',vars2plot,'freq2plot',{'theta'},'times2plot',{'theta'});
%scalpObject.plotScalpMap('vars2plot',vars2plot,'freq2plot',{'beta'},'times2plot',{'beta'});
%scalpObject.plotERPs('vars2plot',vars2plot,'freq2plot',freq2plot,'times2plot',times2plot)
%scalpObject.plotScalpBar('vars2plot',vars2plot,'freq2plot',freq2plot,'times2plot',times2plot)
