%% lucky door 
%% generating brain map differences

%load the peb cleaned file!

% rawData = getdir;
% results = getdir;

rawData = uigetdir('/media/owner/data3/Jason/Active/climateMF/climate_study/control');
results = uigetdir('/media/owner/Seagate Desktop Drive/TT results');

files11 = pickfiles(rawData,{'csucres.xdf'});
% files21 = pickfiles(rawData,{'braine.xdf'});
files21p1 = pickfiles(rawData,{'braine.xdf'});
files22 = pickfiles(rawData,{'BrainE.xdf'});

files = char(files21p1, files22);
Nsubjects = size(files,1);

taskNames = {'go_green','middle_fish','lost_star','lucky_door','face_off','two_tap','rest'};
Ntasks = length(taskNames);
block = 1;


cndnNames = {'stim'};
cndn = 1;

clustypeNames = {'thetaClus','alphaClus','betaClus'};

%src localization
solverType = 'bsbl';
saveFull = true;
account4artifacts = true;

src2roiReductionType = 'mpower';
templatefile = 'headModel_templateFile_newPEBplus.mat'; %headModel.getDefaultTemplateFilename();
%         templatefile = headModel.getDefaultTemplateFilename();
conductivity = [0.33, 0.022, 0.33];
orientation = true;
windowSize = (40/1000)*250; %EEG.srate;

hm = headModel.loadFromFile(templatefile);
T = hm.indices4Structure(hm.atlas.label);
T = double(T)';
P = sparse(bsxfun(@rdivide,T, sum(T,2)))';

%labelnames = {'selective attention','response inhibition','conflict processing','working memory','emotion processing','internal attention','global cognition'};
EEGdata = [];
EEGtimes = [];
EEGbaselinedata = [];
EEGbasetimes = [];

tabledata3 = [];
tabledata4 = [];

% redcapfilename = '/data2/RedCap/BrainE_RedCap_032020.xlsx';
% behdatacol = struct;
% 
% %RedCap Table
% redcaptable = readtable (redcapfilename,'Sheet',1);
% rc = struct;
% rc.age = redcaptable.age;
% rc.pilotid = redcaptable.pilotid;
% rc.pilotid = str2double(erase(string(rc.pilotid),"Pilotii")); % Get only id number for comparison
% rc.depression = redcaptable.depressionscore;
% rc.loneliness = redcaptable.loneliness;
% rc.wisdom = redcaptable.wisdom;

theta2 = [];
beta2 = [];
alpha2 = [];
beta1 = [];
alpha1 = [];

%%
for subject= 1:Nsubjects
    
    tabledata2 = [];
    
    rawFile = deblank(files(subject,:));
    [pathName, rawFileName] = fileparts(rawFile);
    disp(['Processing subject ' rawFileName ' ' num2str(subject) '/' num2str(Nsubjects)]);
    
    try
        pilotname = 'ucsd';
        [~,endIndex] = regexp(rawFileName,'ucsd');
        [startIndex,~] = regexp(rawFileName,'_brainemood');
        subjectindex = str2num(rawFileName(endIndex+1:startIndex-1));
    end
    
    try
        pilotname = 'vadr';
        [~,endIndex] = regexp(rawFileName,'vadr');
        [startIndex,~] = regexp(rawFileName,'_braineii');
        subjectindex = str2num(rawFileName(endIndex+1:startIndex-1));
    end
    
    try
        pilotname = 'vastudy';
        [~,endIndex] = regexp(rawFileName,'vastudy');
        [startIndex,~] = regexp(rawFileName,'_braineii');
        subjectindex = str2num(rawFileName(endIndex+1:startIndex-1));
    end
    
    if isempty(subjectindex)
        pilotname = 'Pilotii';
        [~,endIndex] = regexp(rawFileName,'pilotii');
        [startIndex,~] = regexp(rawFileName,'_braineii');
        subjectindex = str2num(rawFileName(endIndex+1:startIndex-1));
    end
    
    if isempty(subjectindex)
        pilotname = 'Pilotii';
        [~,endIndex] = regexp(rawFileName,'pilotii');
        [startIndex,~] = regexp(rawFileName,'_2_braineii');
        subjectindex = str2num(rawFileName(endIndex+1:startIndex-1));
    end
    
    if isempty(subjectindex)
        pilotname = 'pilot';
        [~,endIndex] = regexp(rawFileName,'pilot');
        [startIndex,~] = regexp(rawFileName,'_braine');
        subjectindex = str2num(rawFileName(endIndex+1:startIndex-1));
    end
    
    if isempty(subjectindex)
        pilotname = 'Pilot';
        [~,endIndex] = regexp(rawFileName,'Pilot');
        [startIndex,~] = regexp(rawFileName,'_BrainE');
        subjectindex = str2num(rawFileName(endIndex+1:startIndex-1));
    end
    
    if isempty(subjectindex)
        pilotname = 'pilot';
        [~,endIndex] = regexp(rawFileName,'pilot');
        [startIndex,~] = regexp(rawFileName,'_BrainE');
        subjectindex = str2num(rawFileName(endIndex+1:startIndex-1));
    end
    
    if isempty(subjectindex)
        pilotname = 'tms';
        [~,endIndex] = regexp(rawFileName,'tms');
        [startIndex,~] = regexp(rawFileName,'_tmse_1_before');
        subjectindex = str2num(rawFileName(endIndex+1:startIndex-1));
    end
    
     if isempty(subjectindex)
        pilotname = 'psyc';
        [~,endIndex] = regexp(rawFileName,'psyc');
        [startIndex,~] = regexp(rawFileName,'_csucres');
        subjectindex = str2num(rawFileName(endIndex+1:startIndex-1));
     end
    
     if isempty(subjectindex)
        pilotname = 'basicn';
        [~,endIndex] = regexp(rawFileName,'basicn');
        [startIndex,~] = regexp(rawFileName,'_csucres');
        subjectindex = str2num(rawFileName(endIndex+1:startIndex-1));
     end
     
    
    meantscoreAct = zeros(68,7);
    taskset_invalid = [];
    %     theta1 = [];
    %     beta1 = [];
    %     alpha1 = [];
    
    for taski= 4
        
        tabledata = zeros(68,499);
        theta = zeros(68,1);
        beta = zeros(68,1);
        alpha = zeros(68,1);
        
        try
            tscoreAct = zeros(3,68);
            
            task = taski;
            block = 1;
            
            %if taski == 7, task = 1; block = 2; end
            
            %3-7 = theta(1), 8-12 = beta(2), 13-30 = alpha(3)
            
            %Focusing on ERS and ERD
            if task == 1
                pcfile.tfinfo{1} = [-1000 1000 3 7];
                pcfile.tfinfo{2} = [-1000 1000 8 12];
                pcfile.tfinfo{3} = [-1000 1000 13 30];
                
            elseif task == 6
                pcfile.tfinfo{1} = [-1000 1000 3 7];
                pcfile.tfinfo{2} = [-1000 1000 8 12];
                pcfile.tfinfo{3} = [-1000 1000 13 30];
                
            elseif task == 3
                pcfile.tfinfo{1} = [-1000 1000 3 7];
                pcfile.tfinfo{2} = [-1000 1000 8 12];
                pcfile.tfinfo{3} = [-1000 1000 13 30];
                
            elseif task == 5
                pcfile.tfinfo{1} = [-1000 1000 3 7];
                pcfile.tfinfo{2} = [-1000 1000 8 12];
                pcfile.tfinfo{3} = [-1000 1000 13 30];
                
            elseif task == 2
                pcfile.tfinfo{1} = [-1000 1000 3 7];
                pcfile.tfinfo{2} = [-1000 1000 8 12];
                pcfile.tfinfo{3} = [-1000 1000 13 30];
                
            elseif task == 4
                pcfile.tfinfo{1} = [-1000 1500 3 30];
                pcfile.tfinfo{2} = [-1000 1500 8 12];
                pcfile.tfinfo{3} = [-1000 1500 13 30];
            end
            
            disp(['Processing task' num2str(task)]);
            
            
            % Create task folder if needed
            subjectTaskFolder = fullfile(results,rawFileName,taskNames{task});
            
            
            EEGmain = struct;
            EEGmain_baseline = struct;
            
            specName = 'tmapStr_newpeak_flank';
            if task == 1
                specName1 = 'nogo1_extd';
            elseif task == 5
                specName1 = '.';
            elseif task == 4
                specName1 = '_ch_';
            else
                specName1 = '.';
            end
            
            %pipeline settings, block info is present
            %                             processvarFiles = pickfiles(subjectTaskFolder, {'processvars', specName1, taskNames{task}, rawFileName, '_stim.mat'});
            processvarFiles = pickfiles(subjectTaskFolder, {'processvars', specName1, taskNames{task}, rawFileName, '_stim.mat'});
            if isempty(processvarFiles)
                if task == 1
                    specName1 = [];
                    specName1 = '_nogo1_';
                    processvarFiles = pickfiles(subjectTaskFolder, {'processvars', specName1, taskNames{task}, rawFileName, '_stim.mat'});
                end
            else
            end
            load(processvarFiles,'pipelineSettings');
            %                             Tblock = pipelineSettings.trialsetind{blocki}(pipelineSettings.blockind{blocki});
            trialsetind = pipelineSettings.trialsetind{block};
            blockind = (pipelineSettings.blockind{block});
            
            if task == 1 
            trialsetind = [pipelineSettings.trialsetind{1}];
            blockind = [pipelineSettings.blockind{1};pipelineSettings.blockind{2}];
            end
            
            
            if task == 4
                bl = unique(pipelineSettings.behtabledata.Block(abs(pipelineSettings.behtabledata.ChoiceValue) == 30));
                ex = unique(pipelineSettings.behtabledata.Block(abs(pipelineSettings.behtabledata.ChoiceValue) == 20));
                choice_b1_rarel = contains(pipelineSettings.behtabledata.ChoiceSeries,'RareL') & (pipelineSettings.behtabledata.Block == bl(end));
                choice_ex_rarel = contains(pipelineSettings.behtabledata.ChoiceSeries,'RareL') & (pipelineSettings.behtabledata.Block == ex(end));
                
                %baseline (sameEV) RareL
                trialsind{1} = (trialsetind(blockind(choice_b1_rarel(trialsetind(blockind)))));
                
                %diffEV RareL
                trialsind{2} = (trialsetind(blockind(choice_ex_rarel(trialsetind(blockind)))));
                
                choice_b1_rareg = contains(pipelineSettings.behtabledata.ChoiceSeries,'RareG') & (pipelineSettings.behtabledata.Block == bl(end));
                choice_ex_rareg = contains(pipelineSettings.behtabledata.ChoiceSeries,'RareG') & (pipelineSettings.behtabledata.Block == ex(end));
                
                %sameEV RareG
                trialsind{3} = (trialsetind(blockind(choice_b1_rareg(trialsetind(blockind)))));
                
                %diffEV RareG
                trialsind{4} = (trialsetind(blockind(choice_ex_rareg(trialsetind(blockind)))));
            end
            
            %                             DesignMatrix = pipelineSettings.Xglm{block};
            RTall = pipelineSettings.behtabledata.ResponseTime;
            %         valueall = pipelineSettings.behtabledata.ChoiceValue;
            %accall = pipelineSettings.behtabledata.Accuracy;
            
            
            pebFiles = pickfiles(subjectTaskFolder,{'peb_cleaned_ch_epoch_ch', taskNames{task},rawFileName, 'stim.set'});
            
            EEGmain = pop_loadset(deblank(pebFiles));
            
            %when running OA only! 
            
%             if subject == 26 
%                 EEGmain.data = EEGmain.data([1:7 9:13 15:22 24], :, :);
%                 EEGmain.chanlocs(8) = [];
%                 EEGmain.chanlocs(14) = [];
%                 EEGmain.chanlocs(23) = [];
%             end 
%             
%             if subject == 13 
%                 EEGmain.data = EEGmain.data([1:11 13:24], :, :);
%                 EEGmain.chanlocs(12) = [];
%             end
%             
%             if subject == 38 
%                 EEGmain.data = EEGmain.data([1:4 6:24], :, :);
%                 EEGmain.chanlocs(5) = [];
%             end

            %match group 
%             if subject == 43 
%                 EEGmain.data = EEGmain.data([1 3 5:12 14:24],:,:);
%                 EEGmain.chanlocs(2) = [];
%                 EEGmain.chanlocs(4) = [];
%                 EEGmain.chanlocs(13) = [];
%             end
            %suicidality group 
%             if subject == 3
%                 EEGmain.data = EEGmain.data([1:13 15:24],:,:);
%                 EEGmain.chanlocs(14) = [];
%             end
%             if subject == 4 
%                 EEGmain.data = EEGmain.data([1:5 7:14 16:24],:,:);
%                 EEGmain.chanlocs(6) = [];
%                 EEGmain.chanlocs(15) = [];
%             end
            
            %checking channel config
            load('chanlocs_file.mat');
            if size(EEGmain.data,1) < 24
                EEGmain_dummy = EEGmain;
                
                [a,b] = intersect({chanlocs_str.labels}, {EEGmain_dummy.chanlocs.labels});
                chinterp = setdiff(1:24,b);
                
                for chinterp1 = chinterp
                EEGmain_dummy.chanlocs(chinterp1+1:end+1) = EEGmain_dummy.chanlocs(chinterp1:end);
                EEGmain_dummy.chanlocs(chinterp1) = chanlocs_str(chinterp1);
                EEGmain_dummy.nbchan = length(EEGmain_dummy.chanlocs);
                
                EEGmain_dummy.data(chinterp1+1:end+1,:,:)  = EEGmain_dummy.data(chinterp1:end,:,:);
                EEGmain_dummy.data(chinterp1,:,:)  = zeros(1,size(EEGmain_dummy.data,2),size(EEGmain_dummy.data,3));
                end
                %spherical interpolation
                EEGmain = eeg_interp(EEGmain_dummy,chinterp);
            end
            
            EEGmainref = EEGmain;
            
            
            %now iterate through all the freq bands, filter, and save the
            %activity
            
            trialsind{1} = trialsetind(blockind);
            trialsfac = 1; %factor to be multiplied with every trial
            faci_tilda = 1;
            
            
            
            filename = ['brainmap_']; % taskNames{task} '_block_' num2str(blocki)]; % '_'  clustypeNames{clustype}];
            
            
            for clustype = 1 %:3 %pos, neg
                
                
                
                EEG = EEGmain;
                
                
                tfinfoarr = pcfile.tfinfo{clustype};
                
                if isempty(tfinfoarr) == 0
                    timesind1 = find(EEG.times > tfinfoarr(1));
                    timesind1 = timesind1(1);
                    timesind2 = find(EEG.times < tfinfoarr(2));
                    timesind2 = timesind2(end);
                    
                    %adding the following to save the pop
                    %inverse solution
                    if (timesind2 - timesind1) < 80
                        timesind2 = timesind1 + 80;
                    end
                    
                    freqind1 = tfinfoarr(3); freqind2 = tfinfoarr(4);
                    
                    %                             filename = ['EEGsrc_full_'  faciName{faci} '_'  num2str(facii) '_PC_' num2str(pci) '_' clustypeNames{clustype} '_' num2str(clusti)];
                    
                    fileEEGsrc =  pickfiles(results, {filename});
                    %                 if isempty(fileEEGsrc)
                    
                     %EEG = pop_eegfiltnew(EEGmain, freqind1, freqind2, 826,0,[],0);
                    
                     for tri = 1:4
                         if tri == 1
                             
                             trialtype = trialsind{1};
                             
                             EEG.data = permute([permute(EEG.data(:,timesind1:timesind2,trialtype(trialtype < size(EEG.data,3))),[3 1 2]).*trialsfac],[2 3 1]);
                             EEG.data = nanmean(EEG.data,3);
                             
                             EEG.times = EEGmain.times(timesind1:timesind2);
                             EEG.trials = size(EEG.data,3); %length(trialsind{facii});
                             EEG.pnts = length(EEG.times);
                             
                             EEG.epoch = [];
                             EEG.etc.subjectNames = subjectindex;
                             
                             EEG = pop_forwardModel(EEG, templatefile, conductivity, orientation);
                             EEG = pop_pebp(EEG,saveFull, account4artifacts, src2roiReductionType);
                             
                             %need flanks to avoid edge
                             %artifacts while taking
                             %envelope (each side flank
                             %is of 500 msec in
                             %length(125 samples) in 256
                             %sampling rate
                             flanklen = 125;
                             %                                                     flankact = reshape(randsample(EEG.etc.src.act(:),size(EEG.etc.src.act,1)*(flanklen)),[size(EEG.etc.src.act,1), flanklen]);
                             %                                                     flankactFull = reshape(randsample(EEG.etc.src.actFull(:),size(EEG.etc.src.actFull,1)*(flanklen)),[size(EEG.etc.src.actFull,1), flanklen]);
                             %                                                 flanki_s = randi(size(EEG.etc.src.act,2)*size(EEG.etc.src.act,3)-flanklen-1,[68,1]); flankact = EEG.etc.src.act(:,flanki_s:flanki_s+flanklen-1);
                             
                             EEG1 = EEG; EEG2 = [];
                             flanki_s = randi(size(EEG.etc.src.act,2)*size(EEG.etc.src.act,3)-flanklen-1,[68,1]); flankact = EEG.etc.src.act(:,flanki_s:flanki_s+flanklen-1);
                             
                             %                         if subi == 1
                             EEG1.etc.src.act(:,1:end+flanklen*2,1) = envelope(cat(2,flankact,EEG.etc.src.act(:,:,1),flankact)',1,'peak')';
                             %                         else
                             %                             EEG1.etc.src.act(:,:,subi) = envelope(cat(2,flankact,EEG.etc.src.act(:,1:end,subi),flankact)',1,'peak')';
                             %                         end
                             %restoring the time
                             %length
                             tempvar_1 = EEG1.etc.src.act(:,flanklen+1:end-flanklen,1);
                             EEG2.etc.src.act(:,:,1) = tempvar_1;
                             
                             EEG.etc.src.act = EEG2.etc.src.act; clear EEG1 EEG2;
                             
                         elseif tri == 2
                             
                             EEG = EEGmain;
                             trialtype = trialsind{2};
                             
                             EEG.data = permute([permute(EEG.data(:,timesind1:timesind2,trialtype(trialtype < size(EEG.data,3))),[3 1 2]).*trialsfac],[2 3 1]);
                             EEG.data = nanmean(EEG.data,3);
                             
                             EEG.times = EEGmain.times(timesind1:timesind2);
                             EEG.trials = size(EEG.data,3); %length(trialsind{facii});
                             EEG.pnts = length(EEG.times);
                             
                             EEG.epoch = [];
                             EEG.etc.subjectNames = subjectindex;
                             
                             EEG = pop_forwardModel(EEG, templatefile, conductivity, orientation);
                             EEG = pop_pebp(EEG,saveFull, account4artifacts, src2roiReductionType);
                             
                             %need flanks to avoid edge
                             %artifacts while taking
                             %envelope (each side flank
                             %is of 500 msec in
                             %length(125 samples) in 256
                             %sampling rate
                             flanklen = 125;
                             %                                                     flankact = reshape(randsample(EEG.etc.src.act(:),size(EEG.etc.src.act,1)*(flanklen)),[size(EEG.etc.src.act,1), flanklen]);
                             %                                                     flankactFull = reshape(randsample(EEG.etc.src.actFull(:),size(EEG.etc.src.actFull,1)*(flanklen)),[size(EEG.etc.src.actFull,1), flanklen]);
                             %                                                 flanki_s = randi(size(EEG.etc.src.act,2)*size(EEG.etc.src.act,3)-flanklen-1,[68,1]); flankact = EEG.etc.src.act(:,flanki_s:flanki_s+flanklen-1);
                             
                             EEG1 = EEG; EEG2 = [];
                             flanki_s = randi(size(EEG.etc.src.act,2)*size(EEG.etc.src.act,3)-flanklen-1,[68,1]); flankact = EEG.etc.src.act(:,flanki_s:flanki_s+flanklen-1);
                             
                             %                         if subi == 1
                             EEG1.etc.src.act(:,1:end+flanklen*2,1) = envelope(cat(2,flankact,EEG.etc.src.act(:,:,1),flankact)',1,'peak')';
                             %                         else
                             %                             EEG1.etc.src.act(:,:,subi) = envelope(cat(2,flankact,EEG.etc.src.act(:,1:end,subi),flankact)',1,'peak')';
                             %                         end
                             %restoring the time
                             %length
                             tempvar_2 = EEG1.etc.src.act(:,flanklen+1:end-flanklen,1);
                             EEG2.etc.src.act(:,:,1) = tempvar_2;
                             
                             EEG.etc.src.act = EEG2.etc.src.act; clear EEG1 EEG2;
                         elseif tri == 3 
                             
                             EEG = EEGmain;
                             trialtype = trialsind{3};
                             
                             EEG.data = permute([permute(EEG.data(:,timesind1:timesind2,trialtype(trialtype < size(EEG.data,3))),[3 1 2]).*trialsfac],[2 3 1]);
                             EEG.data = nanmean(EEG.data,3);
                             
                             EEG.times = EEGmain.times(timesind1:timesind2);
                             EEG.trials = size(EEG.data,3); %length(trialsind{facii});
                             EEG.pnts = length(EEG.times);
                             
                             EEG.epoch = [];
                             EEG.etc.subjectNames = subjectindex;
                             
                             EEG = pop_forwardModel(EEG, templatefile, conductivity, orientation);
                             EEG = pop_pebp(EEG,saveFull, account4artifacts, src2roiReductionType);
                             
                             %need flanks to avoid edge
                             %artifacts while taking
                             %envelope (each side flank
                             %is of 500 msec in
                             %length(125 samples) in 256
                             %sampling rate
                             flanklen = 125;
                             %                                                     flankact = reshape(randsample(EEG.etc.src.act(:),size(EEG.etc.src.act,1)*(flanklen)),[size(EEG.etc.src.act,1), flanklen]);
                             %                                                     flankactFull = reshape(randsample(EEG.etc.src.actFull(:),size(EEG.etc.src.actFull,1)*(flanklen)),[size(EEG.etc.src.actFull,1), flanklen]);
                             %                                                 flanki_s = randi(size(EEG.etc.src.act,2)*size(EEG.etc.src.act,3)-flanklen-1,[68,1]); flankact = EEG.etc.src.act(:,flanki_s:flanki_s+flanklen-1);
                             
                             EEG1 = EEG; EEG2 = [];
                             flanki_s = randi(size(EEG.etc.src.act,2)*size(EEG.etc.src.act,3)-flanklen-1,[68,1]); flankact = EEG.etc.src.act(:,flanki_s:flanki_s+flanklen-1);
                             
                             %                         if subi == 1
                             EEG1.etc.src.act(:,1:end+flanklen*2,1) = envelope(cat(2,flankact,EEG.etc.src.act(:,:,1),flankact)',1,'peak')';
                             %                         else
                             %                             EEG1.etc.src.act(:,:,subi) = envelope(cat(2,flankact,EEG.etc.src.act(:,1:end,subi),flankact)',1,'peak')';
                             %                         end
                             %restoring the time
                             %length
                             tempvar_3 = EEG1.etc.src.act(:,flanklen+1:end-flanklen,1);
                             EEG2.etc.src.act(:,:,1) = tempvar_3;
                             
                             EEG.etc.src.act = EEG2.etc.src.act; clear EEG1 EEG2;
                         else
                             EEG = EEGmain;
                             
                             trialtype = trialsind{4};
                             
                             EEG.data = permute([permute(EEG.data(:,timesind1:timesind2,trialtype(trialtype < size(EEG.data,3))),[3 1 2]).*trialsfac],[2 3 1]);
                             EEG.data = nanmean(EEG.data,3);
                             
                             EEG.times = EEGmain.times(timesind1:timesind2);
                             EEG.trials = size(EEG.data,3); %length(trialsind{facii});
                             EEG.pnts = length(EEG.times);
                             
                             EEG.epoch = [];
                             EEG.etc.subjectNames = subjectindex;
                             
                             EEG = pop_forwardModel(EEG, templatefile, conductivity, orientation);
                             EEG = pop_pebp(EEG,saveFull, account4artifacts, src2roiReductionType);
                             
                             %need flanks to avoid edge
                             %artifacts while taking
                             %envelope (each side flank
                             %is of 500 msec in
                             %length(125 samples) in 256
                             %sampling rate
                             flanklen = 125;
                             %                                                     flankact = reshape(randsample(EEG.etc.src.act(:),size(EEG.etc.src.act,1)*(flanklen)),[size(EEG.etc.src.act,1), flanklen]);
                             %                                                     flankactFull = reshape(randsample(EEG.etc.src.actFull(:),size(EEG.etc.src.actFull,1)*(flanklen)),[size(EEG.etc.src.actFull,1), flanklen]);
                             %                                                 flanki_s = randi(size(EEG.etc.src.act,2)*size(EEG.etc.src.act,3)-flanklen-1,[68,1]); flankact = EEG.etc.src.act(:,flanki_s:flanki_s+flanklen-1);
                             
                             EEG1 = EEG; EEG2 = [];
                             flanki_s = randi(size(EEG.etc.src.act,2)*size(EEG.etc.src.act,3)-flanklen-1,[68,1]); flankact = EEG.etc.src.act(:,flanki_s:flanki_s+flanklen-1);
                             
                             %                         if subi == 1
                             EEG1.etc.src.act(:,1:end+flanklen*2,1) = envelope(cat(2,flankact,EEG.etc.src.act(:,:,1),flankact)',1,'peak')';
                             %                         else
                             %                             EEG1.etc.src.act(:,:,subi) = envelope(cat(2,flankact,EEG.etc.src.act(:,1:end,subi),flankact)',1,'peak')';
                             %                         end
                             %restoring the time
                             %length
                             tempvar_4 = EEG1.etc.src.act(:,flanklen+1:end-flanklen,1);
                             EEG2.etc.src.act(:,:,1) = tempvar_4;
                             
                             EEG.etc.src.act = EEG2.etc.src.act; clear EEG1 EEG2;
                             
                             
                         end
                     end
                     

%                     %%%%%%%%%%%%%%%%plotting
%                     a = [find(EEG.times > 200)];
%                     b = [find(EEG.times > 500)];
%                     peakt = [a(1) a(1) b(1)];
%                     lenvals = [25 25 25]; %make alpha 2 peaks [25 13 7]
%                     lengthpeak = 1;
%                     basetimeinds = 63:112;
%                     
%                     Actvalue{subject,taski,clustype} = tempvar;
%                     tempvar1 = nanmean(tempvar(:,peakt(clustype)-lenvals(clustype)*lengthpeak:peakt(clustype)+lenvals(clustype)*lengthpeak),2)-nanmean(tempvar(:,basetimeinds),2);
%                     tempvar1 = tempvar; %(:,peakt(clustype)-lenvals(clustype)*lengthpeak:peakt(clustype)+lenvals(clustype)*lengthpeak)-nanmean(tempvar(:,basetimeinds),2);

                    GLbias1 = tempvar_1 - tempvar_3;
                    
                    EVloss1 = tempvar_2 - tempvar_1;
                    
                    EVgain1 = tempvar_4 - tempvar_3;
                    %record all ROI for each clustype

                    
                    
                    %taskinew = [1 3 4 0 5 0 2];
                    
                    %                     meantemplate = nanmean(templateAct{taskinew(taski),clustype},2);
                    %                     stdtemplate = nanstd(templateAct{taskinew(taski),clustype},[],2);
                    %
                    %                     zscoreAct = (tempvar1 - meantemplate)./stdtemplate;
                    %                     tscoreAct(clustype,:) = (zscoreAct * 10) + 50;
                    %
                    
                    
                    
                end
            end
            
            
            
        end
        
        %         meantscoreAct(:,taski) = squeeze(nanmean(tscoreAct,1));
        %         if ~any(find(tscoreAct(:)))
        %             taskset_invalid = [taskset_invalid taski];
        %         end

        
    end
    
    subjectid = erase(rawFileName,"_braineii");
    subjectnumber = str2double(erase(subjectid,"pilotii"));
    
    %index = find(rc.pilotid == subjectnumber);
    
%     if  ~isempty (index)
%         behdatacol(subject).redcap(1) = (rc.pilotid(index));
%         behdatacol(subject).redcap(2) = str2double(char(rc.depression(index)));
%         behdatacol(subject).redcap(3) = str2double(char(rc.loneliness(index)));
%         behdatacol(subject).redcap(4) = str2double(char(rc.wisdom(index)));
%         
%     else
%         behdatacol(subject).redcap(1) = nan;
%         behdatacol(subject).redcap(2) = nan;
%         behdatacol(subject).redcap(3) = nan;
%         behdatacol(subject).redcap(4) = nan;
%     end
%     
%     temp = cell2mat({behdatacol.redcap}');
%     
    
%     theta1(:,subject) = theta;
%     beta1(:,subject) = beta;
%     alpha1(:,subject) = alpha;

       GLbias_match(:,:,subject) = GLbias1;
                    
       EVloss_match(:,:,subject) = EVloss1;
                    
       EVgain_match(:,:,subject) = EVgain1;
    
    if isempty(subjectindex)
        subjectindex = nan;
    end
    subjectcoll(subject) = subjectindex;
    
    %tabledata3 = [tabledata2(subject,:), 'VariableNames', {'theta','alpha','beta',[subject]}];
    
    
end