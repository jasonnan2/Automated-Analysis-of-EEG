%% Climate luckydoor Main analysis
% Jason Nan
% 12/19/2023
clear all;close all;clc
cd('C:\Users\jason\OneDrive\Desktop\MS Beng\Neatlabs\Climate_LuckyDoor')
addpath('scripts')
addpath('A:\eeglab_OldLSL_DataAna04072023')
eeglab
load('chanlocs_file.mat') % chanlocs_str
%% Main Analysis with Class
% Create an object of the class
timeAxis=-500:4:1496;

groupPaths = {'A:/ClimateLD/analysis_results/climateLDcontrol_amp_rmv.mat', 'A:/ClimateLD/analysis_results/climateLDwitnessed_amp_rmv.mat','A:/ClimateLD/analysis_results/climateLDhappenedto_amp_rmv.mat'};
variables = {'EVgain', 'EVloss', 'GLbias'};

baselineTime=[-250 -50];

analysis = ERPanalysis(groupPaths,variables,baselineTime); % set up the class, need filepath of group data, variables to run and baseline correction range
analysis = analysis.cleanDatasets(); % remove any missing subjects
finishedAnalysis = analysis.standardPipeline(); % standard processing pipeline, 5SD outlier, baseline correction

%%

% choice 0-500
% immediate reward 500-1000
% cum reward 1000-1500

timeRange.choice=[0 500];
timeRange.imReward=[500 1000];
timeRange.cumReward=[1000 1500];
timeNames = fieldnames(timeRange);

obj=finishedAnalysis;
properties = obj.variables;
groups = obj.groups;

groupNames=cellfun(@(x) x.groupname, groups,'uniformoutput',false);
freq_list=groups{1}.freq_list;
timeAxis=groups{1}.timeAxis;

% Assuming groups is your 1xN cell array
N = length(groups);
combinations = nchoosek(1:N,2); % This will give you a matrix where each row is a combination of two groups

for p=1:length(properties)
    property=properties{p};
    for comb = 1:size(combinations, 1)
        firstSub=combinations(comb, 1);
        secondSub=combinations(comb, 2);
        group1 = groups{firstSub}.(property);
        group2 = groups{secondSub}.(property);
        
        figure
        figCount=0;
        for t=1:length(timeNames)
            for f=1:length(freq_list)
                figCount = figCount + 1;
                subplot(3,4,figCount)
                
                [~, startIdx] = min(abs(timeAxis - timeRange.(timeNames{t})(1)));
                [~, endIdx] = min(abs(timeAxis - timeRange.(timeNames{t})(2)));
                
                s1=squeeze(nanmean(group1(f,:,startIdx:endIdx,:),3)); % get to chan x sub
                s2=squeeze(nanmean(group2(f,:,startIdx:endIdx,:),3)); % get to chan x sub
                plotSigTopo(s1,s2,chanlocs_str,{'.','+'})
                if t==1
                    title(freq_list{f})
                end
                if f==1
                    text(-1,-.2, timeNames{t},'fontweight','bold','fontsize',16,'rotation',90);         
                end
            end
        end
        sgtitle(property+" "+groupNames(secondSub)+"-"+groupNames(firstSub))
    end
end

%%
function plotSigTopo(s1,s2,chanlocs,key)
    % s1 and s2 are chan x sub size
    plotdata=nanmean(s2,2)-nanmean(s1,2);
    new_chan=topoSignificance(s1,s2,chanlocs,{'.','+'});
    topoplot(plotdata,new_chan,'headrad','rim','electrodes','labels','efontsize' ,16); 
    colorbar(); 
    caxis([min(plotdata),max(plotdata)])
    colormap(parula(128))
    brighten(0.5)
end

%%
function chanlocs=topoSignificance(s1,s2,chanlocs,key)

    % pre and post are data matrices size N subjects x C channels
    % key is a two element cell array with markers for non-sig and sig
    % new_chan is the chan_locs
    % Perform t-tests
    
    s1=squeeze(s1);
    s2=squeeze(s2);
    for chan = 1:length(chanlocs)
        % Get non-NaN indices
        [~,p(chan)] = ttest2(s1(chan,:), s2(chan,:));
        if isnan(p(chan))
            disp('asdfawrsasdbaweasbas')
            waitforbuttonpress
        end
        
    end
    %[~,p]=fdr(reshape(p,[1,numel(p)]),0.05);
    mask=key((p<0.05)+1);
    [chanlocs.labels]=mask{:};
end
        
    