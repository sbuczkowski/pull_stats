function pull_stats_airxbcal_site(year,siteID);

%*******************************************************************************
% sequentially read in a year's worth of site rtp files. For each
% day read in, collect and bin rtp data by site ID from
% prof.iudef(2,:) and concatenate to rolling tally. At end of year,
% output separate files for each site
%*******************************************************************************

% $$$ addpath /asl/matlib/h4tools
% $$$ addpath /asl/rtp_prod/airs/utils
% $$$ addpath /asl/packages/rtp_prod2/util
% $$$ addpath /home/sergio/MATLABCODE/PLOTTER  %
% $$$                                          % equal_area_spherical_bands

% record run start datetime in output stats file for tracking
trace.RunDate = datetime('now','TimeZone','local','Format', ...
                         'd-MMM-y HH:mm:ss Z');
trace.Reason = 'Normal site pull_stats runs';
trace.klayers = false;
trace.droplayers = false;

cstr =[ 'bits1-4=NEdT[0.08 0.12 0.15 0.20 0.25 0.30 0.35 0.4 0.5 0.6 0.7' ...
  ' 0.8 1.0 2.0 4.0 nan]; bit5=Aside[0=off,1=on]; bit6=Bside[0=off,1=on];' ...
  ' bits7-8=calflag&calchansummary[0=OK, 1=DCR, 2=moon, 3=other]' ];

basedir = ['/asl/rtp/rtp_airxbcal_v5/' int2str(year) '/site'];
dayfiles = dir(fullfile(basedir, 'era_airxbcal_day*_site.rtp'));
fprintf(1,'>>> numfiles = %d\n', length(dayfiles));

% $$$ nsites = 20;   % use the 20 canonical AIRS site
nchans = 2378; % AIRS channel space
nlevels = 60;    % ERA

% initialize accumulators
accum_robs = [];  
accum_rcalc = [];
accum_rtime = [];
accum_rlat = [];
accum_rlon = [];
accum_satzen = [];
accum_solzen = [];
accum_spres = [];
accum_stemp = [];
accum_ptemp = [];
accum_gas_1 = [];
accum_gas_3 = [];
accum_iudef1 = [];
accum_iudef2 = [];
accum_iudef4 = [];
accum_nlevs = [];
accum_plevs = [];

iday = 1;
% $$$ for giday = 1:100:length(dayfiles)
for giday = 1:length(dayfiles)
   fprintf(1, '>>> year = %d  :: giday = %d\n', year, giday);
   a = dir(fullfile(basedir,dayfiles(giday).name));
   a.bytes;
   if a.bytes > 100000
      [h,ha,p,pa] = rtpread(fullfile(basedir,dayfiles(giday).name));
      f = h.vchan;  % AIRS proper frequencies

   % get count of obs at each site
   % selection reason is stored in p.iudef(1,:)
   % site ID is stored in p.iudef(2,:)
% $$$    edges = 0.5:1.0:20.5;
% $$$    [nsiteobs, edges, siteID] = histcounts(p.iudef(2,:), edges);
% $$$    maxobs = max(nsiteobs);  % find site with most obs this day

   siteind = find(p.iudef(2,:) == siteID);
   nsiteobs = length(siteind);

   if ~nsiteobs
       continue;
   end
   
% $$$    % allocate space to store this day's obs
% $$$    robs = NaN(nchans, nsiteobs, 'single');
% $$$    rcalc = NaN(nchans, nsiteobs, 'single');
% $$$    rtimes = NaN(nsiteobs, 'double');
% $$$    rlat = NaN(nsiteobs, 'single');
% $$$    rlon = NaN(nsiteobs, 'single');
% $$$    satzen = NaN(nsiteobs, 'single');
% $$$    solzen = NaN(nsiteobs, 'single');
% $$$    spres = NaN(nsiteobs, 'single');
% $$$    stemp = NaN(nsiteobs, 'single');
% $$$    ptemp = NaN(nlevels, nsiteobs, 'single');
% $$$    gas_1 = NaN(nlevels, nsiteobs, 'single');
% $$$    gas_3 = NaN(nlevels, nsiteobs, 'single');
% $$$    iudef1 = NaN(nsiteobs, 'single');
% $$$    iudef2 = NaN(nsiteobs, 'single');
% $$$    iudef4 = NaN(nsiteobs, 'single');
% $$$    nlevs = NaN(nsiteobs, 'single');
% $$$    plevs = NaN(nlevels, nsiteobs, 'single');

   % Look for bad channels and initialize counts
   [nedt,ab,ical] = calnum_to_data(p.calflag,cstr);
   n = length(p.rlat);
   count_all = ones(2378,n);
   for i=1:2378
       % Find bad channels
       k = find( p.robs1(i,:) == -9999 | ical(i,:) ~= 0 | nedt(i,:) > 1);
       %          % These are the good channels
       %          kg = setdiff(1:n,k);
       % NaN's for bad channels
       p.robs1(i,k) = NaN;
       count_all(i,k) = 0;
   end   % for i=1:2378
    
   % accumulate into the annual tally
   accum_rtime = cat(2, accum_rtime, p.rtime(siteind));  
   accum_robs = cat(2, accum_robs, p.robs1(:,siteind));
   accum_rcalc = cat(2, accum_rcalc, p.rcalc(:,siteind));
   accum_rlat = cat(2, accum_rlat, p.rlat(siteind));   
   accum_rlon = cat(2, accum_rlon, p.rlon(siteind));   
   accum_satzen = cat(2, accum_satzen, p.satzen(siteind)); 
   accum_solzen = cat(2, accum_solzen, p.solzen(siteind)); 
   accum_spres = cat(2, accum_spres, p.spres(siteind));  
   accum_stemp = cat(2, accum_stemp, p.stemp(siteind));  
   accum_ptemp = cat(2, accum_ptemp, p.ptemp(:,siteind));
   accum_gas_1 = cat(2, accum_gas_1, p.gas_1(:,siteind));
   accum_gas_3 = cat(2, accum_gas_3, p.gas_3(:,siteind));
   accum_iudef1 = cat(2, accum_iudef1, p.iudef(1,siteind));
   accum_iudef2 = cat(2, accum_iudef2, p.iudef(2,siteind));
   accum_iudef4 = cat(2, accum_iudef4, p.iudef(4,siteind));
   accum_nlevs = cat(2, accum_nlevs, p.nlevs(siteind));  
   accum_plevs = cat(2, accum_plevs, p.plevs(:,siteind));
   
   iday = iday + 1
   end   % if a.bytes > 1000000

end  % giday

% $$$ outfileDir='/asl/data/stats/airs';
outfileDir='/home/sbuczko1/WorkingFiles/data/stats/airs';

outfileName = sprintf('airxbcal_site-%d_%d_era_rad_stats', siteID, year);
outfilePath = fullfile(outfileDir, outfileName)
eval_str = sprintf('save %s accum_*', outfilePath);
eval(eval_str);
end

