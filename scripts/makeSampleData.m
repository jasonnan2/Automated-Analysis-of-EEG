%% Generate Project structure (This will be format to use package) 
load('sample.mat','chanlocs','roi') % loading in chanlocs and roi samples

project=struct();
project.info=struct();
project.info.freq_list = {'theta'  'alpha'  'beta'  'broadband'};
project.info.timeAxis = -500:4:1496;
project.info.groupNames = {'Group1','Group2','Group3'};
project.info.variables = {'NeuralVarName'};
project.info.chanlocs = chanlocs;
project.info.roi = roi;
project.info.experimentalDesign = 'twoSample';

%% Generating the sample Scalp Data
rng(42)
% Parameters
nFreq = 4; nChan = 24; nTime = 500; nSubj = 10;
subList = arrayfun(@(x) sprintf('sub%d', x), 1:nSubj, 'UniformOutput', false);
fs = 250; % Sampling rate
timeVec = linspace(-0.5, 1.5, nTime); % in seconds

% Initialize data with Gaussian noise
sharedNoise = randn(nFreq, nChan, nTime, nSubj)/4;
Group1 = sharedNoise+randn(nFreq, nChan, nTime, nSubj)/20; 
Group2 = sharedNoise+randn(nFreq, nChan, nTime, nSubj)/20; 
Group3 = sharedNoise+randn(nFreq, nChan, nTime, nSubj)/20; 

% Insert ERP into Group1: theta (1), channels 10 & 14
erpTime = (timeVec >= 0 & timeVec <= 0.5); % 0–500ms
erp = exp(-((timeVec - 0.2).^2) / (2 * 0.02^2)); % peak at 200ms, sigma=50ms
erp = erp .* erpTime; % zero outside 0–500ms
Group1(1,10,:,:) = squeeze(Group1(1,10,:,:)) + repmat(erp', 1,nSubj);
Group1(1,14,:,:) = squeeze(Group1(1,14,:,:)) + repmat(erp', 1,nSubj);

% Insert ERP into Group2: alpha (2), channel 19
erpTime = (timeVec >= 0 & timeVec <= 0.5); % 0–500ms
erp = exp(-((timeVec - 0.3).^2) / (2 * 0.02^2)); % peak at 300ms, sigma=50ms
erp = erp .* erpTime; % zero outside 0–500ms
Group2(2,19,:,:) = squeeze(Group2(2,19,:,:)) + repmat(erp', 1,nSubj);

% Insert ERP into Group3: alpha (2), channel 15
erpTime = (timeVec >= 0 & timeVec <= 0.5); % 0–500ms
erp = exp(-((timeVec - 0.25).^2) / (2 * 0.02^2)); % peak at 250ms, sigma=50ms
erp = erp .* erpTime; % zero outside 0–500ms
Group3(2,19,:,:) = squeeze(Group3(2,19,:,:)) + repmat(erp', 1,nSubj);

% Assign to project struct
project.scalpData.Group1.NeuralVarName = Group1;
project.scalpData.Group2.NeuralVarName = Group2;
project.scalpData.Group3.NeuralVarName = Group3;

% Replace with actual subIDs
project.scalpData.Group1.subList = subList;
project.scalpData.Group2.subList = subList;
project.scalpData.Group3.subList = subList;


%% Generating sample Source Data

nROI=68;
% Initialize data with Gaussian noise
sharedNoise = randn(nFreq, nROI, nTime, nSubj)/4;
Group1 = sharedNoise+randn(nFreq, nROI, nTime, nSubj)/20; 
Group2 = sharedNoise+randn(nFreq, nROI, nTime, nSubj)/20; 
Group3 = sharedNoise+randn(nFreq, nROI, nTime, nSubj)/20; 

sharedAlphaSource = exp(-((timeVec - 0.3).^2) / (2 * 0.03^2)) .* sin(2*pi*10*timeVec/1000);
sharedAlphaSource = 500*repmat(sharedAlphaSource, nSubj, 1)';  

% Insert source into Group1 and Group2 for sample theta connectivity
for roi=[15 16 21 22 51 52] % pDMN rois
    Group1(1,roi,:,:) = squeeze(Group1(1,roi,:,:)) + sharedAlphaSource;
end
for roi=[2 47 48 63 64] % VAN rois
    Group1(1,roi,:,:) = squeeze(Group1(1,roi,:,:)) + sharedAlphaSource;
end


% Insert source for ROI differences 
extraSource = exp(-((timeVec - 0.1).^2) / (2 * 0.03^2)) .* sin(2*pi*10*timeVec/1000);
extraSource = repmat(extraSource, nSubj, 1)';
Group3(2,42,:,:) = squeeze(Group3(2,42,:,:)) + 500*extraSource;
Group3(2,43,:,:) = squeeze(Group3(2,43,:,:)) + 500*extraSource;

% Assign to project struct
project.sourceData.Group1.NeuralVarName = Group1;
project.sourceData.Group2.NeuralVarName = Group2;
project.sourceData.Group3.NeuralVarName = Group3;

% Replace with actual subIDs
project.sourceData.Group1.subList = subList;
project.sourceData.Group2.subList = subList;
project.sourceData.Group3.subList = subList;


save('./sample.mat','project','-append')

