%% Neatlabs Interns EEGProcessing 
cd('/media/owner/data3/Jason/Active/climateLD')
addpath('scripts')
addpath ('/media/owner/data3/eeglab')
eeglab
%% Set data and result folder paths
clear
clc

%%% Change here
rawData =  '/media/owner/Seagate Desktop Drive/climate study/happenedto';
results = '/media/owner/data3/Jason/Active/climateLD/processed_data/happenedto';
%%% 

files1 = pickfiles(rawData,{'.xdf'});
files = char(files1);
Nsubjects = size(files,1);
if ~exist(results,'dir')
    mkdir(results);
end

%% Set parameters
subjectIndices = 1:Nsubjects;
runOnChannels = true;
runOnSources = false;
taskNames = {'go_green','middle_fish','lost_star','lucky_door','face_off','two_tap'};
Ntasks = length(taskNames);
%---------
% PEB+ parameters
windowSize = (40/1000)*500; % 40 ms
overlaping = 25;
solverType = 'bsbl';
saveFull = false;
account4artifacts = true;
postprocCallback = [];
src2roiReductionType = 'm2power';

%templatefile = headModel.getDefaultTemplateFilename();
templatefile = '/media/owner/data3/eeglab/plugins/BrainEAnalysis/headModel_templateFile_32chan.mat';
conductivity = [0.33, 0.022, 0.33];
orientation = true;
%---------
files = files(subjectIndices,:);
fid = fopen(fullfile(results,char("ERROR_log_"+char(datetime('today')))),'w');
c = onCleanup(@()fclose(fid));

%% Run pipeline for every subject
for subject= Nsubjects
    
    rawFile = deblank(files(subject,:));
    if ~isempty(rawFile)
        [pathName, rawFileName] = fileparts(rawFile);
        
        % Create subject folder if needed
        subjectFolder = fullfile(results,rawFileName);
        if ~exist(subjectFolder,'dir')
            mkdir(subjectFolder);
        end
        
        disp(['Processing subject ' rawFileName ' ' num2str(subject) '/' num2str(Nsubjects)]);
        for task= 4
            disp(['Processing task' num2str(task)]);
            try
                
                % Create task folder if needed
                subjectTaskFolder = fullfile(results,rawFileName,taskNames{task});
%                                 %remove directory
%                                 if exist(subjectTaskFolder,'dir')
%                                     rmdir(subjectTaskFolder,'s');
%                                 end
                %create a folder directory
                if ~exist(subjectTaskFolder,'dir')
                    mkdir(subjectTaskFolder);
                end
                
                % Epoching
                disp('Epoching');
                %                 epochFiles = pickfiles(subjectTaskFolder,{'epoch_','.set'});
                if task == 1
                    epochFiles = pickfiles(subjectTaskFolder,{'processvars_nogo1_','.mat'});
                elseif task == 4
                    epochFiles = pickfiles(subjectTaskFolder,{'processvars_ch_','.mat'});
                else
                    epochFiles = pickfiles(subjectTaskFolder,{'processvars_','.mat'});
                end
                if isempty(epochFiles)
                    behFile = pickfiles(pathName,{'.xlsx'},{'.xlsx'},{'lock'});
                    if isempty(behFile)
                        behFile = pickfiles(pathName,'trial_summary.csv');
                    end
                    switch taskNames{task}
                        case 'go_green'
                            if size(behFile,1) ~= 1
                                behFileind = find(contains(cellstr(behFile),'gowait'));
                            else
                                behFileind = 1;
                            end
                            [EEGset, ProcessVars] = go_wait.preProcessing(rawFile, deblank(behFile(behFileind,:)));
                        case 'middle_fish'
                            if size(behFile,1) ~= 1
                                behFileind = find(contains(cellstr(behFile),'middlefish'));
                            else
                                behFileind = 1;
                            end
                            [EEGset, ProcessVars] = middle_fish.preProcessing(rawFile, deblank(behFile(behFileind,:)));
                        case 'lost_star'
                            if size(behFile,1) ~= 1
                                behFileind = find(contains(cellstr(behFile),'loststar'));
                            else
                                behFileind = 1;
                            end
                            [EEGset, ProcessVars] = lost_star.preProcessing(rawFile, deblank(behFile(behFileind,:)));
                        case 'lucky_door'
                            if size(behFile,1) ~= 1
                                behFileind = find(contains(cellstr(behFile),'luckydoor'));
                            else
                                behFileind = 1;
                            end
                            [EEGset, ProcessVars] = lucky_door.preProcessing(rawFile, deblank(behFile(behFileind,:)));
                        case 'face_off'
                            if size(behFile,1) ~= 1
                                behFileind = find(contains(cellstr(behFile),'faceoff'));
                            else
                                behFileind = 1;
                            end
                            [EEGset, ProcessVars] = face_off.preProcessing(rawFile, deblank(behFile(behFileind,:)));
                        case 'two_tap'
                            if size(behFile,1) ~= 1
                                behFileind = find(contains(cellstr(behFile),'twotap'));
                            else
                                behFileind = 1;
                            end
                            [EEGset, ProcessVars] = two_tap.preProcessing(rawFile, deblank(behFile(behFileind,:)));
                    end
                    for k=1:length(EEGset)
                        pipelineSettings = ProcessVars(k);
                        if task == 1
                            EEGset{k}.etc.pipelineSettingsFile = fullfile(subjectTaskFolder,['processvars_nogo1_' taskNames{task} '_' rawFileName '_' EEGset{k}.setname,'.mat']);
                            pop_saveset(EEGset{k},'filepath',subjectTaskFolder,'filename',['epoch_nogo1_' taskNames{task} '_' rawFileName '_' EEGset{k}.setname],'savemode','onefile');
                        elseif task == 4
                            EEGset{k}.etc.pipelineSettingsFile = fullfile(subjectTaskFolder,['processvars_ch_' taskNames{task} '_' rawFileName '_' EEGset{k}.setname,'.mat']);
                            pop_saveset(EEGset{k},'filepath',subjectTaskFolder,'filename',['epoch_ch_' taskNames{task} '_' rawFileName '_' EEGset{k}.setname],'savemode','onefile');
                            
                        else
                            EEGset{k}.etc.pipelineSettingsFile = fullfile(subjectTaskFolder,['processvars_' taskNames{task} '_' rawFileName '_' EEGset{k}.setname,'.mat']);
                            pop_saveset(EEGset{k},'filepath',subjectTaskFolder,'filename',['epoch_' taskNames{task} '_' rawFileName '_' EEGset{k}.setname],'savemode','onefile');
                            
                        end
                        save(EEGset{k}.etc.pipelineSettingsFile,'pipelineSettings');
                    end
                end
                
                % PEB+
                disp('PEB+');
                if task == 1
                    epochFiles = pickfiles(subjectTaskFolder,{'epoch_nogo1_','.set'},{'epoch_nogo1_','.set'},{'peb'});
                    nepochFiles = size(epochFiles,1);
                    pebFiles = pickfiles(subjectTaskFolder,{'peb_cleaned_nogo1_','.set'});
                    npebFiles = size(pebFiles,1);
                elseif task == 4
                    epochFiles = pickfiles(subjectTaskFolder,{'epoch_ch_','.set'},{'epoch_ch_','.set'},{'peb'});
                    nepochFiles = size(epochFiles,1);
                    pebFiles = pickfiles(subjectTaskFolder,{'peb_cleaned_ch_','.set'});
                    npebFiles = size(pebFiles,1);
                else
                    epochFiles = pickfiles(subjectTaskFolder,{'epoch_','.set'},{'epoch_','.set'},{'peb'});
                    nepochFiles = size(epochFiles,1);
                    pebFiles = pickfiles(subjectTaskFolder,{'peb_cleaned_','.set'});
                    npebFiles = size(pebFiles,1);
                end
                if npebFiles < nepochFiles
                    for k=1:size(epochFiles,1)
                        [~, fileName] = fileparts(deblank(epochFiles(k,:)));
                        if isempty(cell2mat(regexp(string(pebFiles),fileName(end-3:end))))
                            EEG = pop_loadset(deblank(epochFiles(k,:)));
                            windowSize = (40/1000)*EEG.srate; % 40 ms
                            
                            % Run PEB and save cleaned data
                            EEG = pop_forwardModel(EEG, templatefile, conductivity, orientation);
                            EEG = pop_pebp(EEG, saveFull, account4artifacts, src2roiReductionType);
                            if task == 1
                                pop_saveset(EEG,'filepath',subjectTaskFolder,'filename',['peb_cleaned_nogo1_' fileName],'savemode','onefile');
                                
                                EEG = moveSource2DataField(EEG, ' ');                        % Save ROI data
                                pop_saveset(EEG,'filepath',subjectTaskFolder,'filename',['peb_roi_nogo1_' fileName],'savemode','onefile');
                            elseif task == 4
                                pop_saveset(EEG,'filepath',subjectTaskFolder,'filename',['peb_cleaned_ch_' fileName],'savemode','onefile');
                                
                                EEG = moveSource2DataField(EEG, ' ');                        % Save ROI data
                                pop_saveset(EEG,'filepath',subjectTaskFolder,'filename',['peb_roi_ch_' fileName],'savemode','onefile');
                            else
                                pop_saveset(EEG,'filepath',subjectTaskFolder,'filename',['peb_cleaned_' fileName],'savemode','onefile');
                                
                                EEG = moveSource2DataField(EEG, ' ');                        % Save ROI data
                                pop_saveset(EEG,'filepath',subjectTaskFolder,'filename',['peb_roi_' fileName],'savemode','onefile');
                            end
                        end
                    end
                end
                
                
                % Time-frequency decomposition (channels)
                disp('TF channels');
                if runOnChannels
                    if task == 1
                        pebFiles = pickfiles(subjectTaskFolder,{'peb_cleaned_nogo1_','.set'});
                        npebFiles = size(pebFiles,1);
                        tfFiles = pickfiles(subjectTaskFolder,{'power_channel_nogo1_','.mat'},{'power_channel_nogo1_','.mat'},{'_v'});
                        ntfFiles = size(tfFiles,1);
                    elseif task == 4
                        pebFiles = pickfiles(subjectTaskFolder,{'peb_cleaned_ch_','.set'});
                        npebFiles = size(pebFiles,1);
                        tfFiles = pickfiles(subjectTaskFolder,{'power_channel_ch_','.mat'},{'power_channel_ch_','.mat'},{'_v'});
                        ntfFiles = size(tfFiles,1);
                    else
                        pebFiles = pickfiles(subjectTaskFolder,{'peb_cleaned_','.set'});
                        npebFiles = size(pebFiles,1);
                        tfFiles = pickfiles(subjectTaskFolder,{'power_channel_','.mat'},{'power_channel_','.mat'},{'_v'});
                        ntfFiles = size(tfFiles,1);
                    end
                    if ntfFiles < npebFiles
                        if task == 1
                            baselineFile = pickfiles(subjectTaskFolder,{'peb_cleaned_nogo1_','stim.set'});
                            pebFiles = pickfiles(subjectTaskFolder,{'peb_cleaned_nogo1_','.set'});
                        elseif task == 4
                            baselineFile = pickfiles(subjectTaskFolder,{'peb_cleaned_ch_','stim.set'});
                            pebFiles = pickfiles(subjectTaskFolder,{'peb_cleaned_ch_','.set'});
                        else
                            baselineFile = pickfiles(subjectTaskFolder,{'peb_cleaned_','stim.set'});
                            pebFiles = pickfiles(subjectTaskFolder,{'peb_cleaned_','.set'});
                        end
                        pebFiles = setdiff(pebFiles,baselineFile,'rows');
                        if ~isempty(baselineFile) && ~isempty(pebFiles)
                            pebFiles = char(baselineFile,pebFiles);
                        elseif isempty(pebFiles)
                            pebFiles = baselineFile;
                        end
                        for k=1:size(pebFiles,1)
                            [~, fileName] = fileparts(deblank(pebFiles(k,:)));
                            if isempty(cell2mat(regexp(string(tfFiles),fileName(end-3:end))))
                                EEG = pop_loadset(deblank(pebFiles(k,:)));
                                %                                 load(EEG.etc.pipelineSettingsFile,'pipelineSettings');
                                if task == 1
                                    load(fullfile(subjectTaskFolder,['processvars_nogo1_' taskNames{task} '_' rawFileName '_' EEG.setname,'.mat']),'pipelineSettings');
                                elseif task == 4
                                    load(fullfile(subjectTaskFolder,['processvars_ch_' taskNames{task} '_' rawFileName '_' EEG.setname,'.mat']),'pipelineSettings');
                                else
                                    load(fullfile(subjectTaskFolder,['processvars_' taskNames{task} '_' rawFileName '_' EEG.setname,'.mat']),'pipelineSettings');
                                end
                                label = {EEG.chanlocs.labels};
                                if k==1
                                    [TFPower,TFPhase, timeAxis, freqAxis, TFBaseline] = timeFrequencyDecomposition(EEG, pipelineSettings.baseline, pipelineSettings.cutoffFreq);
                                else
                                    [TFPower,TFPhase, timeAxis, freqAxis] = timeFrequencyDecomposition(EEG, TFBaseline, pipelineSettings.cutoffFreq);
                                end
                                if task == 1
                                    save(fullfile(subjectTaskFolder,['power_channel_nogo1_' fileName '.mat']),'TFPower','label','timeAxis', 'freqAxis','pipelineSettings');
                                    save(fullfile(subjectTaskFolder,['phase_channel_nogo1_' fileName '.mat']),'TFPhase','label','timeAxis', 'freqAxis','pipelineSettings');
                                elseif task == 4
                                    save(fullfile(subjectTaskFolder,['power_channel_ch_' fileName '.mat']),'TFPower','label','timeAxis', 'freqAxis','pipelineSettings');
                                    save(fullfile(subjectTaskFolder,['phase_channel_ch_' fileName '.mat']),'TFPhase','label','timeAxis', 'freqAxis','pipelineSettings');
                                    
                                else
                                    save(fullfile(subjectTaskFolder,['power_channel_' fileName '.mat']),'TFPower','label','timeAxis', 'freqAxis','pipelineSettings');
                                    save(fullfile(subjectTaskFolder,['phase_channel_' fileName '.mat']),'TFPhase','label','timeAxis', 'freqAxis','pipelineSettings');
                                    
                                end
                            end
                        end
                    end
                end
                
                
            catch ME
                ME.message
 %                                               keyboard
                fprintf(fid,'%s\t%s\t%s\n',rawFileName,taskNames{task},ME.message);
            end
        end
    end
end
