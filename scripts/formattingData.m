CLIMATELD=struct();
load('A:\ClimateLD\analysis_results\climateLDcontrol.mat')
CLIMATELD.DATA.Control.groupname='Control';
CLIMATELD.DATA.Control.EVgain=EVgain;
CLIMATELD.DATA.Control.EVloss=EVloss;
CLIMATELD.DATA.Control.GLbias=GLbias;
CLIMATELD.DATA.Control.subList=subjectcoll;

load('A:\ClimateLD\analysis_results\climateLDwitnessed.mat')
CLIMATELD.DATA.Witnessed.groupname='Witnessed';
CLIMATELD.DATA.Witnessed.EVgain=EVgain;
CLIMATELD.DATA.Witnessed.EVloss=EVloss;
CLIMATELD.DATA.Witnessed.GLbias=GLbias;
CLIMATELD.DATA.Witnessed.subList=subjectcoll;

load('A:\ClimateLD\analysis_results\climateLDhappenedto.mat')
CLIMATELD.DATA.HappenedTo.groupname='HappenedTo';
CLIMATELD.DATA.HappenedTo.EVgain=EVgain;
CLIMATELD.DATA.HappenedTo.EVloss=EVloss;
CLIMATELD.DATA.HappenedTo.GLbias=GLbias;
CLIMATELD.DATA.HappenedTo.subList=subjectcoll;

CLIMATELD.info.freq_list=freq_list;
CLIMATELD.info.timeAxis=timeAxis;
CLIMATELD.info.groupNames=fieldnames(CLIMATELD.DATA);

CLIMATELD.info.variables=setdiff(fieldnames(CLIMATELD.DATA.Control),{'groupname','subList'});


load('chanlocs_file.mat') % chanlocs_str
CLIMATELD.info.chanlocs=chanlocs_str;
save('A:\ClimateLD\analysis_results\ClimateLD_allgroups.mat','CLIMATELD')


