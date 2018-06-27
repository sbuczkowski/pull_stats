function config = config()

% Input file information
config.rtpFileDir = '/asl/data/rtp_airibrad_v5/';
config.rtpFileGlob = 'era_airibrad*_random.rtp';
config.asType = 'random';

% Output file information
config.psOutPutDir = '/home/sbuczko1/testoutput/2015/airs/random/';

% run klayers?
config.runKlayers = false;

% latitude binning
config.nbins = 20; % generate 2*nbins bins in total

% robs filter specification
config.sFilter = 'find(p.iudef(4,:) == 68 & p.landfrac == 1)';  % night + land
config.sDescriptor = '_night_land';