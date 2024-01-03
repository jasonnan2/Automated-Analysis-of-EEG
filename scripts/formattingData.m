CLIMATELD=struct();
% Scalp Data
load('A:\ClimateLD\analysis_results\climateLDcontrol_amp_rmv.mat')
CLIMATELD.scalpData.Control.groupname='Control';
CLIMATELD.scalpData.Control.EVgain=EVgain;
CLIMATELD.scalpData.Control.EVloss=EVloss;
CLIMATELD.scalpData.Control.GLbias=GLbias;
CLIMATELD.scalpData.Control.subList=subjectcoll;

load('A:\ClimateLD\analysis_results\climateLDwitnessed_amp_rmv.mat')
CLIMATELD.scalpData.Witnessed.groupname='Witnessed';
CLIMATELD.scalpData.Witnessed.EVgain=EVgain;
CLIMATELD.scalpData.Witnessed.EVloss=EVloss;
CLIMATELD.scalpData.Witnessed.GLbias=GLbias;
CLIMATELD.scalpData.Witnessed.subList=subjectcoll;

load('A:\ClimateLD\analysis_results\climateLDhappenedto_amp_rmv.mat')
CLIMATELD.scalpData.HappenedTo.groupname='HappenedTo';
CLIMATELD.scalpData.HappenedTo.EVgain=EVgain;
CLIMATELD.scalpData.HappenedTo.EVloss=EVloss;
CLIMATELD.scalpData.HappenedTo.GLbias=GLbias;
CLIMATELD.scalpData.HappenedTo.subList=subjectcoll;

% source Data
load('A:\ClimateLD\analysis_results\climateLDcontrol_source_amp_rmv.mat')
CLIMATELD.sourceData.Control.groupname='Control';
CLIMATELD.sourceData.Control.EVgain=EVgain;
CLIMATELD.sourceData.Control.EVloss=EVloss;
CLIMATELD.sourceData.Control.GLbias=GLbias;
CLIMATELD.sourceData.Control.subList=subjectcoll;

load('A:\ClimateLD\analysis_results\climateLDwitnessed_source_amp_rmv.mat')
CLIMATELD.sourceData.Witnessed.groupname='Witnessed';
CLIMATELD.sourceData.Witnessed.EVgain=EVgain;
CLIMATELD.sourceData.Witnessed.EVloss=EVloss;
CLIMATELD.sourceData.Witnessed.GLbias=GLbias;
CLIMATELD.sourceData.Witnessed.subList=subjectcoll;

load('A:\ClimateLD\analysis_results\climateLDhappenedto_source_amp_rmv.mat')
CLIMATELD.sourceData.HappenedTo.groupname='HappenedTo';
CLIMATELD.sourceData.HappenedTo.EVgain=EVgain;
CLIMATELD.sourceData.HappenedTo.EVloss=EVloss;
CLIMATELD.sourceData.HappenedTo.GLbias=GLbias;
CLIMATELD.sourceData.HappenedTo.subList=subjectcoll;

CLIMATELD.info.freq_list=freq_list;
CLIMATELD.info.timeAxis=timeAxis;
CLIMATELD.info.groupNames=fieldnames(CLIMATELD.scalpData);

CLIMATELD.info.variables=setdiff(fieldnames(CLIMATELD.scalpData.Control),{'groupname','subList'});

load('chanlocs_file.mat') % chanlocs_str
CLIMATELD.info.chanlocs=chanlocs_str;

load('68roi.mat')
CLIMATELD.info.roi=roi;


save('A:\ClimateLD\analysis_results\ClimateLD_allgroups_amp_rmv.mat','CLIMATELD')



