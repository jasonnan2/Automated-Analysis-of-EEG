
%%
% >>Also see if the frontal activity shows significant differences when individuals are collapsed by trauma (see final column in the 3rd spreadsheet in attached xlsx, 2= yes trauma, 1- no trauma), 
% and/or if it shows correlation with EVgain metrics.

% 2.	Two groups split by trauma with current processing pipeline

%%
% >>For GLbias can you confirm if the parallel behavior model (as for EVgain below, robust fitlm) doesn’t show any behavior effect or does it?
% GLbias ~ 1 + ExpGroup + age + gender + ethnicity + sesscore + anxiety + depression

tbl=readtable("A:\ClimateLD\climateTrauma_beh_redcap.xlsx",'sheet','LD_beh_results')
tbl.gender = categorical(tbl.gender);
tbl.ethnicity = categorical(tbl.ethnicity);

var_name = 'Glbias ~ 1 + group + age + gender + ethnicity + sesscore + anxiety + depression';

model = fitlm(tbl,var_name,'RobustOpts','on') 

