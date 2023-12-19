%% Climate luckydoor analysis EEG 
% Jason Nan
% 12/19/2023
clear all;close all;clc
cd('/media/owner/data3/Jason/Active/climateLD')
addpath ('/media/owner/data3/eeglab')
addpath('scripts')
eeglab
%% lucky door 
%% generating brain map differences
fid = fopen(fullfile(char("ERROR_log_"+char(datetime('today')))),'w');
c = onCleanup(@()fclose(fid));
%load the peb cleaned file!

% rawData = getdir;
% results = getdir;

rawData = uigetdir('/media/owner/data3/Jason/Active/climateMF/climate_study/control');
results = uigetdir('/media/owner/Seagate Desktop Drive/TT results');

files11 = pickfiles(rawData,{'csucres.xdf'});
% files21 = pickfiles(rawData,{'braine.xdf'});
files21p1 = pickfiles(rawData,{'braine.xdf'});
files22 = pickfiles(rawData,{'BrainE.xdf'});

files = char(files21p1, files22); % for Control
%files = char(files11); % for witnessed and happenedto list only
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
% size freq x chan x time  x Subjects
GLbias=nan(4,24,500,Nsubjects);
EVloss=nan(4,24,500,Nsubjects); 
EVgain=nan(4,24,500,Nsubjects);
%results=uigetdir('/media/owner/data3/Jason/Active/climateLD/processed_data_DONTNEED/happenedto'); % run for [13 15 16 17 18 19] in happened to
for subject= 1:Nsubjects%[13 15 16 17 18 19]%
    
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
     
    try

        task = 4;
        block = 1;
        pcfile.tfinfo{1} = [-1000 1500 3 7];
        pcfile.tfinfo{2} = [-1000 1500 8 12];
        pcfile.tfinfo{3} = [-1000 1500 13 30];
        pcfile.tfinfo{4} = [-1000 1500 1 45];



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
        processvarFiles = pickfiles(subjectTaskFolder, {'processvars', specName1, taskNames{task}, rawFileName, '_stim.mat'});
        

        load(processvarFiles,'pipelineSettings');
        trialsetind = pipelineSettings.trialsetind{block};
        blockind = (pipelineSettings.blockind{block});

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
 

        RTall = pipelineSettings.behtabledata.ResponseTime;
        %         valueall = pipelineSettings.behtabledata.ChoiceValue;
        %accall = pipelineSettings.behtabledata.Accuracy;

        pebFiles = pickfiles(subjectTaskFolder,{'peb_cleaned_ch_epoch_ch', taskNames{task},rawFileName, 'stim.set'});

        EEGmain = pop_loadset(deblank(pebFiles));

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

        EEGmainref = EEGmain; % store clean EEGmain

        %now iterate through all the freq bands, filter, and save the
        %activity

        trialsind{1} = trialsetind(blockind);
        trialsfac = 1; %factor to be multiplied with every trial
        faci_tilda = 1;

        for f=1:4 % iterate frequency
            EEGmain=EEGmainref;
            tfinfoarr = pcfile.tfinfo{f};
            timesind1 = find(EEGmain.times > tfinfoarr(1));
            timesind1 = timesind1(1);
            timesind2 = find(EEGmain.times < tfinfoarr(2));
            timesind2 = timesind2(end);

            %adding the following to save the pop
            %inverse solution
            if (timesind2 - timesind1) < 80
                timesind2 = timesind1 + 80;
            end
            freqind1 = tfinfoarr(3); freqind2 = tfinfoarr(4);
            EEGfilt = pop_eegfiltnew(EEGmain, freqind1, freqind2, 826,0,[],0);

            for tri = 1:4

                EEG=EEGfilt;
                trialtype = trialsind{tri};

                EEG.data = permute([permute(EEG.data(:,timesind1:timesind2,trialtype(trialtype < size(EEG.data,3))),[3 1 2]).*trialsfac],[2 3 1]);
                broadband = permute([permute(EEGmain.data(:,timesind1:timesind2,trialtype(trialtype < size(EEGmain.data,3))),[3 1 2]).*trialsfac],[2 3 1]);
                amp_rmv = find(squeeze(max(max(abs(broadband))))>100);
                EEG.data(:,:,amp_rmv)=[];
                

                EEG.data = nanmean(EEG.data,3);

                EEG.times = EEGmain.times(timesind1:timesind2);
                EEG.trials = size(EEG.data,3); %length(trialsind{facii});
                EEG.pnts = length(EEG.times);

                EEG.epoch = [];
                EEG.etc.subjectNames = subjectindex;

                tempData(f,:,:,tri)=EEG.data; clear EEG

            end
        end
            GLbias(:,:,:,subject) = tempData(:,:,:,1) - tempData(:,:,:,3);      
            EVloss(:,:,:,subject) = tempData(:,:,:,2) - tempData(:,:,:,1);         
            EVgain(:,:,:,subject) = tempData(:,:,:,4) - tempData(:,:,:,3);
            clear tempData
    catch ME
        ME.message
        fprintf(fid,'%i\t%s\t%s\n',subject,rawFileName,ME.message);
    end
    
    subjectcoll{subject} = rawFileName;
end


