clear
tbl=readtable("A:\ClimateLD\ClimateStudy_LD_Beh2.xlsx",'sheet','All');
%tbl.gender = categorical(tbl.gender);
%tbl.ethnicity = categorical(tbl.ethnicity);
tbl.TwoGroup = categorical(tbl.TwoGroup);
tbl.traumaBin = categorical(tbl.recenttrauma>=2);
var_name = 'ex_WS_RareG ~ 1 + Group  + age + gender + ethnicity + sesscore + anxiety + depression';
model = fitlm(tbl,var_name,'RobustOpts','on') 

%model = fitlm(tbl,"ex_WS_RareG ~1+age+TwoGroup",'RobustOpts','on') 

