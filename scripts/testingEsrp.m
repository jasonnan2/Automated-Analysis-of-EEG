clear all;close all;clc
cd('C:\Users\jason\OneDrive\Desktop\MS Beng\Neatlabs\Climate_LuckyDoor')
addpath('scripts')
addpath('scripts/DataAnalysis')
addpath('scripts/DataAnalysis/functions')
%% Loading in data and definitions
clear all;clc
dataPath='A:\ClimateLD\analysis_results\formattedData\ClimateLD_allgroups_exRareG_amp_rmv.mat';
load(dataPath)
baselineTime=[-250 -50];
timeRange.a=[-250 500];
scalpObject=ScalpObject(CLIMATELD.scalpData, CLIMATELD.info, baselineTime, timeRange); % Initialize scalp processing
scalpObject = scalpObject.cleanDatasets(); 
%% 

p=scalpObject.info.variables;
groups=scalpObject.info.groupNames;

for g=1:length(groups)
    data=double(squeeze(scalpObject.DATA.(groups{g}).(p{1})(4,:,:,:)));
    EEG=struct();
    for n=1:size(data,3)
        EEG.data=data(:,:,n);
        EEG.trials=size(EEG.data,3);
        EEG.nbchan=size(EEG.data,1);
        EEG.times=-500:4:1496;
        EEG.srate=250;
        try
            [ERSP_Pxx, ~, EEGtimes, freq, ~] = timeFrequencyDecomposition_jasonTest(EEG, [], [1 45]);
        catch
            ERSP_Pxx=nan(24,30,1,474);
        end
        erspData.DATA.(groups{g}).(p{1})(:,:,:,n)=squeeze(ERSP_Pxx);
        
    end
end


%%
dataPath='A:\ClimateLD\analysis_results\formattedData\ClimateLD_ersp_allgroups_exRareG_amp_rmv.mat';
load(dataPath)
CLIMATELD.info.experimentalDesign='twoSample';
CLIMATELD.erspData.Control.exRareg=erspData.DATA.Control.exRareg;
CLIMATELD.erspData.Witnessed.exRareg=erspData.DATA.Witnessed.exRareg;
CLIMATELD.erspData.HappenedTo.exRareg=erspData.DATA.HappenedTo.exRareg;

%%
baselineTime=[-250 -50];
% Defining time ranges of interest
timeRange=struct();
timeRange.base=[-250 500];

freqRanges.theta=[3 8];
freqRanges.alpha=[8 13];
freqRanges.beta=[13,30];

erspObject=ERSPObject(CLIMATELD.erspData, CLIMATELD.info, baselineTime, timeRange); % Initialize scalp processing
erspObject = erspObject.cleanDatasets(); % remove any missing subjects if any from Scalp data
erspObject = erspObject.standardProcessing('outlier','3MAD'); % standard processing pipeline, 5SD outlier, baseline correction
scalpObj = erspObject.createERP(freqRanges); % convert ersp data to ERP data by averaging freq. 

%%
close all
erspObject.plotERSPgroups('chans2plot',{'Fz','Pz'},'groups2plot',{'HappenedTo','Witnessed','Control'},'cRange',[-0.1 .1])

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





















