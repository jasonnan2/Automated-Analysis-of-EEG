dataPath='A:\ClimateLD\analysis_results\formattedData\ClimateLD_rest.mat';
%dataPath='A:/ClimateLD/analysis_results/formattedData/ClimateLD_allgroups_exRareG_amp_rmv.mat';
load(dataPath)
CLIMATELD.info.experimentalDesign='twoSample';
baselineTime=[nan];
% Defining time ranges of interest
timeRange=struct();
timeRange.all=[0 3996];

scalpObject=ScalpObject(CLIMATELD.scalpData, CLIMATELD.info, baselineTime, timeRange); % Initialize scalp processing
scalpObject = scalpObject.cleanDatasets(); % remove any missing subjects if any from Scalp data
scalpObject = scalpObject.standardProcessing('calBaseline',0); % standard processing pipeline, 5SD outlier, baseline correction
scalpObject = scalpObject.ScalpAnalysis(); % calcualtes all significant electrodes between groups across all conditions
%func=@(x) abs(hilbert(x));
%scalpObject = scalpObject.applyFunc(func);
scalpObject = scalpObject.calSigTbl(); % creating table of significant neural attributes

%% Scalp Plotting
close all
vars2plot={'rest'}; freq2plot={'alpha'}; times2plot={'all'}; errorType='sem';
scalpObject.plotScalpMap('vars2plot',vars2plot,'freq2plot',freq2plot,'times2plot',times2plot,'combinations',[1,2]);
%scalpObject.plotERPs('vars2plot',vars2plot,'freq2plot',freq2plot,'times2plot',times2plot)
%scalpObject.plotScalpBar('vars2plot',vars2plot,'freq2plot',freq2plot,'times2plot',times2plot)

