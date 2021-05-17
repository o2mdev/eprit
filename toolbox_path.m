 % restoredefaultpath;
 
 epr_toolbox_path = 'D:\EPR-IT';
  
 % Common routines
 addpath([epr_toolbox_path, filesep, 'common'], '-begin')

 % Report generators
  addpath([epr_toolbox_path, filesep, 'Reports'], '-begin')
  
 % Forward projection
  addpath([epr_toolbox_path, filesep, 'radon'], '-begin')
  
  % FBP image reconstruction 32 and 64 bit
  addpath([epr_toolbox_path, filesep, 'iradon'], '-begin')

  % multy-stage reconstruction
  addpath([epr_toolbox_path, filesep, 'iradon_mstage'], '-begin')
  
  % EPR imaging
  addpath([epr_toolbox_path, filesep, 'epri'], '-begin')
  addpath([epr_toolbox_path, filesep, 'eprfit'], '-begin')

  % Project viewer 
  addpath([epr_toolbox_path, filesep, 'pviewer'], '-begin')
   
  % Volume visualization
  % addpath([epr_toolbox_path, filesep, '3dparty', filesep, 'sliceomatic2'], '-end')

  % Registration software
  addpath([epr_toolbox_path, filesep, 'Arbuz2.0'], '-begin')

  % Image visualization software
  addpath([epr_toolbox_path, filesep, 'ibGUI'], '-begin')

