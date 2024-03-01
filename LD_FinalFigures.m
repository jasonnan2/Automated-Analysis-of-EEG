%% Final Figures
clear all;close all;clc
cd('C:\Users\jason\OneDrive\Desktop\MS Beng\Neatlabs\Climate_LuckyDoor')
addpath('scripts')
addpath('scripts/DataAnalysis')
addpath('scripts/DataAnalysis/functions')
addpath('A:\eeglab_OldLSL_DataAna04072023')
eeglab
load('A:\ClimateLD\Final\scalpObject.mat')
%% Scalp map alpha band choice
close all
scalpObject.plotScalpMap('freq2plot',{'alpha'},'times2plot',{'choice'});
ylabel('')
title('Directly Exposed - Other','fontsize',16)
sgtitle('')
%% Bar plot
close all
scalpObject.plotScalpBar('freq2plot',{'alpha'},'times2plot',{'choice'},'chans2plot',{'Pz'},'groups',{'HappenedTo','Control'});
xlabel('')
title('Avg Pz Alpha Activity','fontsize',16)
sgtitle('')
legend('Directly Exposed','Other','fontsize',9)
legend('boxoff')

%% Neural Behavior

tbl=readtable("A:\ClimateLD\Final\scalpTables.xlsx");
baseModel="ex_WS_RareG ~ 1 + age +group*exRareg_choice_alpha_Pz";
model = fitlm(tbl,baseModel,'RobustOpts','on');

control_neural=tbl.exRareg_choice_alpha_Pz(strcmp(tbl.group,'Control'));
happened_neural=tbl.exRareg_choice_alpha_Pz(strcmp(tbl.group,'HappenedTo'));
control_beh=tbl.ex_WS_RareG(strcmp(tbl.group,'Control'));
happened_beh=tbl.ex_WS_RareG(strcmp(tbl.group,'HappenedTo'));

figure
hold on
[r,p] = plotSpearmans(happened_neural,happened_beh,'r')
[r,p] = plotSpearmans(control_neural,control_beh,'b')
xlabel('Avg Pz Activity','fontweight','bold','fontsize',12);
ylabel(['Win Stay on High',newline, 'Expected Value Choices'],'fontweight','bold','fontsize',12);
[~, hobj, ~, ~] =legend('','Directly Exposed','','Other','fontweight','bold','fontsize',9);
legend('boxoff')
hl = findobj(hobj,'type','line');
set(hl,'LineWidth', 2);
%% single scalp
%load('A:\ClimateLD\Final\allgroup_scalpObject.mat')
%close all
subplot(1,3,1)
s = nanmean(squeeze(nanmean(scalpObject.getGroupData('HappenedTo','exRareg','alpha','choice'),3)),2);
topoplot(s,scalpObject.info.chanlocs,'headrad','rim','electrodes','labels','efontsize',9);
colorbar()
caxis([-1.25,1.25])
colormap(parula(128))
brighten(0.5)
title('Directly Exposed','fontsize',16)

subplot(1,3,2)
s = nanmean(squeeze(nanmean(scalpObject.getGroupData('Witnessed','exRareg','alpha','choice'),3)),2);
topoplot(s,scalpObject.info.chanlocs,'headrad','rim','electrodes','labels','efontsize',9);
colorbar()
caxis([-1.25,1.25])
colormap(parula(128))
brighten(0.5)
title('Indirectly Exposed','fontsize',16)

subplot(1,3,3)
s = nanmean(squeeze(nanmean(scalpObject.getGroupData('Control','exRareg','alpha','choice'),3)),2);
topoplot(s,scalpObject.info.chanlocs,'headrad','rim','electrodes','labels','efontsize',9);
colorbar()
caxis([-1.25,1.25])
colormap(parula(128))
brighten(0.5)
title('Non-exposed','fontsize',16)
%%
function [r,p] = plotSpearmans(var1,var2,c)
    P=polyfit(var1,var2,1);
    Bfit=polyval(P,var1);
    hold on
    plot(var1,var2,c+".",'markersize',12)
    plot(var1,Bfit,c+"-")
    hold off
    axis square
    [r,p]=corr(var1,var2,'Type','Spearman','Rows','Complete');
end
