%% Generate Project structure (This will be format to use package) 
load('sample.mat','chanlocs','roi') % Load example channel locations and ROI definitions

% Initialize the main project structure
project = struct();
project.info = struct();

% Define frequency bands, time axis, group labels, and variable names
project.info.freq_list = {'theta', 'alpha', 'beta', 'broadband'};
project.info.timeAxis = -500:4:1496; % Time axis in milliseconds
project.info.groupNames = {'Group1','Group2','Group3'};
project.info.variables = {'NeuralVarName'}; % Name of neural data field
project.info.chanlocs = chanlocs; % Channel location structure
project.info.roi = roi;           % Region of interest structure
project.info.experimentalDesign = 'twoSample'; % Experimental design type
%% Replace with actual data - sample data generated below
% % Store generated scalp data into the project structure
% project.scalpData.Group1.NeuralVarName = Group1EEG;
% project.scalpData.Group2.NeuralVarName = Group2EEG;
% project.scalpData.Group3.NeuralVarName = Group3EEG;
% 
% % Assign subject IDs
% project.scalpData.Group1.subList = subList1;
% project.scalpData.Group2.subList = subList2;
% project.scalpData.Group3.subList = subList3;
% 
% % Store generated source data into the project structure
% project.sourceData.Group1.NeuralVarName = Group1Source;
% project.sourceData.Group2.NeuralVarName = Group2Source;
% project.sourceData.Group3.NeuralVarName = Group3Source;
% 
% % Assign subject IDs
% project.sourceData.Group1.subList = subList1;
% project.sourceData.Group2.subList = subList2;
% project.sourceData.Group3.subList = subList3;

%% Generating the sample Scalp Data
rng(42) % For reproducibility

% Define data dimensions
nFreq = 4; nChan = 24; nTime = 500; nSubj = 10;
subList = arrayfun(@(x) sprintf('sub%d', x), 1:nSubj, 'UniformOutput', false);
fs = 250; % Sampling rate in Hz
timeVec = linspace(-0.5, 1.5, nTime); % Time vector in seconds

% Create shared Gaussian noise base for all groups
sharedNoise = randn(nFreq, nChan, nTime, nSubj)/4;

% Create scalp EEG data for each group with small variations
Group1 = sharedNoise + randn(nFreq, nChan, nTime, nSubj)/20; 
Group2 = sharedNoise + randn(nFreq, nChan, nTime, nSubj)/20; 
Group3 = sharedNoise + randn(nFreq, nChan, nTime, nSubj)/20; 

% Insert ERP-like signal into Group1: theta band (1), channels 10 & 14
erpTime = (timeVec >= 0 & timeVec <= 0.5);
erp = exp(-((timeVec - 0.2).^2) / (2 * 0.02^2)) .* erpTime;
Group1(1,10,:,:) = squeeze(Group1(1,10,:,:)) + repmat(erp', 1, nSubj);
Group1(1,14,:,:) = squeeze(Group1(1,14,:,:)) + repmat(erp', 1, nSubj);

% Insert ERP into Group2: alpha band (2), channel 19
erp = exp(-((timeVec - 0.3).^2) / (2 * 0.02^2)) .* erpTime;
Group2(2,19,:,:) = squeeze(Group2(2,19,:,:)) + repmat(erp', 1, nSubj);

% Insert ERP into Group3: alpha band (2), channel 19
erp = exp(-((timeVec - 0.25).^2) / (2 * 0.02^2)) .* erpTime;
Group3(2,19,:,:) = squeeze(Group3(2,19,:,:)) + repmat(erp', 1, nSubj);

% Store generated scalp data into the project structure
project.scalpData.Group1.NeuralVarName = Group1;
project.scalpData.Group2.NeuralVarName = Group2;
project.scalpData.Group3.NeuralVarName = Group3;

% Assign subject IDs
project.scalpData.Group1.subList = subList;
project.scalpData.Group2.subList = subList;
project.scalpData.Group3.subList = subList;

%% Generating sample Source Data

nROI = 68; % Number of source ROIs

% Create shared noise and group-specific data for source signals
sharedNoise = randn(nFreq, nROI, nTime, nSubj)/4;
Group1 = sharedNoise + randn(nFreq, nROI, nTime, nSubj)/20; 
Group2 = sharedNoise + randn(nFreq, nROI, nTime, nSubj)/20; 
Group3 = sharedNoise + randn(nFreq, nROI, nTime, nSubj)/20; 

% Generate alpha oscillation burst for insertion
sharedAlphaSource = exp(-((timeVec - 0.3).^2) / (2 * 0.03^2)) .* sin(2*pi*10*timeVec/1000);
sharedAlphaSource = 500 * repmat(sharedAlphaSource, nSubj, 1)';  % Scale & repeat for each subject

% Insert alpha oscillation into specific ROIs for Group1 (e.g., pDMN and VAN)
for roi = [15 16 21 22 51 52] % pDMN
    Group1(1,roi,:,:) = squeeze(Group1(1,roi,:,:)) + sharedAlphaSource;
end
for roi = [2 47 48 63 64] % VAN
    Group1(1,roi,:,:) = squeeze(Group1(1,roi,:,:)) + sharedAlphaSource;
end

% Insert additional signal into Group3: alpha band (2), ROIs 42 and 43
extraSource = exp(-((timeVec - 0.1).^2) / (2 * 0.03^2)) .* sin(2*pi*10*timeVec/1000);
extraSource = repmat(extraSource, nSubj, 1)';
Group3(2,42,:,:) = squeeze(Group3(2,42,:,:)) + 500 * extraSource;
Group3(2,43,:,:) = squeeze(Group3(2,43,:,:)) + 500 * extraSource;

% Store generated source data into the project structure
project.sourceData.Group1.NeuralVarName = Group1;
project.sourceData.Group2.NeuralVarName = Group2;
project.sourceData.Group3.NeuralVarName = Group3;

% Assign subject IDs
project.sourceData.Group1.subList = subList;
project.sourceData.Group2.subList = subList;
project.sourceData.Group3.subList = subList;

% Append generated project data back to the sample file
save('./sample.mat','project','-append')
