**Code Files**

scripts/DataAnalysis - source code for objects and required functions

scripts/runERP.m - sample code for running different types of ERP and source analysis

scripts/runERSP.m - sample code for running ERSP analysis

Input Data Structure

**Info               | Structure with metadata about the dataset**

  -experimentalDesign | string of either ‘paired’ or ‘twoSample’ denoting which ttest to preform
  
  -freq_list          | cell array of frequency names in same order as dataset ex {‘theta’,’alpha’]
  
  -timeAxis           | 1xT vector of timestamps which the datapoints correspond to in mS
  
  -groupNames         | cell array of the names of the groups in dataset. Ex {‘pre’,’post’}
  
  -variables          | cell array of variables you have
  
  -Chanlocs           | chan locs of your EEG from eeglab
  
  -Roi                | only for source processing, cell array of the names of the roi’s from BSBL
  
  
**DATA               | Structure of data, each field is a different group**

  -EX, DATA.pre would be a structure for the pre group, DATA.post would be for post group
  
  -Each group structure contains:
  
  -Groupname     | string of group name
  
  -Variable name | matrix of size F frequencies, C channels, T time points, N subjects
  
    Ex: DATA.pre.Var1 ; DATA.pre.Var2 = rand(F_freq, C_chan, T_time, N_sub)
    
  -subList       | cell array of subject names corresponding to the 4th dimension of variable name
  
    Ex: DATA.subList = {'sub1','sub2','sub3','sub4'}
    
  -timeRange     | structure with 2 element vectors for different timings of interest in mS.
  
    Ex: timeRange.choice=[100 200] would indicate the ‘choice’ period is in the 100-200mS range
    

**Sample Structure found in sample.mat**




