%*************************************************
% Execute user-defined paths *********************
REPOBASEPATH = '/home/sbuczko1/git/';
% $$$ REPOBASEPATH = '/asl/packages/';

addpath ~/git/pull_stats_DEV/util
% rtp_prod2_DEV
addpath(sprintf('%s/rtp_prod2_DEV/util', REPOBASEPATH));
addpath(sprintf('%s/rtp_prod2_DEV/grib', REPOBASEPATH));
addpath(sprintf('%s/rtp_prod2_DEV/emis', REPOBASEPATH));
addpath(genpath(sprintf('%s/rtp_prod2_DEV/airs', REPOBASEPATH)));

% swutils (will move under matlib soon)
addpath(sprintf('%s/swutils', REPOBASEPATH));

% matlib
addpath(sprintf('%s/matlib/clouds/sarta', REPOBASEPATH));  % driver_cloudy_sarta
addpath(sprintf('%s/matlib/rtptools', REPOBASEPATH));   % for cat_rtp

addpath /asl/matlib/aslutil
addpath /asl/matlib/time
%*************************************************
