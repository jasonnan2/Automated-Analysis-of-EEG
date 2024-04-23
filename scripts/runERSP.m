clear all;close all;clc
addpath('scripts')
addpath('scripts/DataAnalysis')
addpath('scripts/DataAnalysis/functions')
%% Loading in data and definitions
dataPath='';
load(dataPath)
CLIMATELD.info.experimentalDesign='twoSample';
%% ERSP Analysis
baselineTime=[-250 -50];
% Defining time ranges of interest
timeRange=struct();
timeRange.base=[-250 500];
%timeRange.choice=[0 500];
freqRanges.theta=[3 8];
freqRanges.alpha=[8 13];
freqRanges.beta=[13,30];

erspObject=ERSPObject(CLIMATELD.erspData, CLIMATELD.info, baselineTime, timeRange); % Initialize scalp processing
erspObject = erspObject.cleanDatasets(); % remove any missing subjects if any from Scalp data
erspObject = erspObject.standardProcessing('outlier','3MAD'); % standard processing pipeline, 5SD outlier, baseline correction
% scalpObj = erspObject.createERP(freqRanges); % convert ersp data to ERP data by averaging freq. 

%%
close all
erspObject.plotERSPgroups('chans2plot',{'Fz','Pz'},'groups2plot',{'HappenedTo','Control'},'cRange',[-0.5 1])

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
