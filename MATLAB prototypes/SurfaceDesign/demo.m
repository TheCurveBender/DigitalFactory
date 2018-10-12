%Demonstrate

addpath('nurbs_toolbox')
addpath('Etc')

opt=MultiOpt();
opt.nBlocks = 5; %Number of blocks
opt.yvals  = [0 2.75 5.5]; %y-values for design curves
[curves] = makeSurfaceDesign(opt);
