function pull_stats_airxbcal_site(year);

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

nsites = 20;   % use the 20 canonical AIRS site
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
   edges = 0.5:1.0:20.5;
   [nsiteobs, edges, siteID] = histcounts(p.iudef(2,:), edges);
   maxobs = max(nsiteobs);  % find site with most obs this day

   % allocate space to store this day's obs
   robs = NaN(nsites, nchans, maxobs, 'single');
   rcalc = NaN(nsites, nchans, maxobs, 'single');
   rtimes = NaN(nsites, maxobs, 'double');
   rlat = NaN(nsites, maxobs, 'single');
   rlon = NaN(nsites, maxobs, 'single');
   satzen = NaN(nsites, maxobs, 'single');
   solzen = NaN(nsites, maxobs, 'single');
   spres = NaN(nsites, maxobs, 'single');
   stemp = NaN(nsites, maxobs, 'single');
   ptemp = NaN(nsites, nlevels, maxobs, 'single');
   gas_1 = NaN(nsites, nlevels, maxobs, 'single');
   gas_3 = NaN(nsites, nlevels, maxobs, 'single');
   iudef1 = NaN(nsites, maxobs, 'single');
   iudef2 = NaN(nsites, maxobs, 'single');
   iudef4 = NaN(nsites, maxobs, 'single');
   nlevs = NaN(nsites, maxobs, 'single');
   plevs = NaN(nsites, nlevels, maxobs, 'single');

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

   % loop over sites
   for s = 1:nsites
       if ~nsiteobs   % no obs for this site. move to next
           continue;  
       end   % if ~nsiteobs
       
       siteind = find(siteID == s);
       rtimes(s,1:nsiteobs(s)) = p.rtime(siteind);
       robs(s,:,1:nsiteobs(s)) = p.robs1(:,siteind);
       rcalc(s,:,1:nsiteobs(s)) = p.rcalc(:,siteind);
       rlat(s,1:nsiteobs(s)) = p.rlat(siteind);
       rlon(s,1:nsiteobs(s)) = p.rlon(siteind);
       satzen(s,1:nsiteobs(s)) = p.satzen(siteind);
       solzen(s,1:nsiteobs(s)) = p.solzen(siteind); 
       spres(s,1:nsiteobs(s)) = p.spres(siteind);
       stemp(s,1:nsiteobs(s)) = p.stemp(siteind);
       ptemp(s,:,1:nsiteobs(s)) = p.ptemp(:,siteind);
       gas_1(s,:,1:nsiteobs(s)) = p.gas_1(:,siteind);
       gas_3(s,:,1:nsiteobs(s)) = p.gas_3(:,siteind);
       iudef1(s,1:nsiteobs(s)) = p.iudef(1,siteind); 
       iudef2(s,1:nsiteobs(s)) = p.iudef(2,siteind);
       iudef4(s,1:nsiteobs(s)) = p.iudef(4,siteind); 
       nlevs(s,1:nsiteobs(s)) = p.nlevs(siteind);
       plevs(s,:,1:nsiteobs(s)) = p.plevs(:,siteind);
   
   end   % for s = 1:20

   % accumulate into the annual tally
   accum_rtime = cat(2, accum_rtime, rtimes);
   accum_robs = cat(3, accum_robs, robs);
   accum_rcalc = cat(3, accum_rcalc, rcalc);
   accum_rlat = cat(2, accum_rlat, rlat);
   accum_rlon = cat(2, accum_rlon, rlon);
   accum_satzen = cat(2, accum_satzen, satzen);
   accum_solzen = cat(2, accum_solzen, solzen);
   accum_spres = cat(2, accum_spres, spres);
   accum_stemp = cat(2, accum_stemp, stemp);
   accum_ptemp = cat(3, accum_ptemp, ptemp);
   accum_gas_1 = cat(3, accum_gas_1, gas_1);
   accum_gas_3 = cat(3, accum_gas_3, gas_3);
   accum_iudef1 = cat(2, accum_iudef1, iudef1);
   accum_iudef2 = cat(2, accum_iudef2, iudef2);
   accum_iudef4 = cat(2, accum_iudef4, iudef4);
   accum_nlevs = cat(2, accum_nlevs, nlevs);
   accum_plevs = cat(3, accum_plevs, plevs);
   
   iday = iday + 1
   end   % if a.bytes > 1000000

end  % giday

% clear variables for re-use (data /should/ just overwrite but,
% just to be safe)
clear rtimes robs rcalc rlat rlon satzen solzen spres stemp ptemp ...
    gas_1 gas_3 iudef1 iudef2 iudef4 nlevs plevs

% $$$ outfileDir='/asl/data/stats/airs';
outfileDir='/home/sbuczko1/WorkingFiles/data/stats/airs';

% loop over sites and break out individual tallies to separate
% files
for s=1:nsites
    % pull current site out of accumulated arrays
    rtimes = squeeze(accum_rtime(s,:));
    robs = squeeze(accum_robs(s,:,:));
    rcalc = squeeze(accum_rcalc(s,:,:));
    rlat = squeeze(accum_rlat(s,:));
    rlon = squeeze(accum_rlon(s,:));
    satzen = squeeze(accum_satzen(s,:));
    solzen = squeeze(accum_solzen(s,:));
    spres = squeeze(accum_spres(s,:));
    stemp = squeeze(accum_stemp(s,:));
    ptemp = squeeze(accum_ptemp(s,:,:));
    gas_1 = squeeze(accum_gas_1(s,:,:));
    gas_3 = squeeze(accum_gas_3(s,:,:));
    iudef1 = squeeze(accum_iudef1(s,:));
    iudef2 = squeeze(accum_iudef2(s,:));
    iudef4 = squeeze(accum_iudef4(s,:));
    nlevs = squeeze(accum_nlevs(s,:));
    plevs = squeeze(accum_plevs(s,:,:));

    % find NaN entries in obs dimension
    removeIndices = isnan(rtimes);

    % remove NaN obs from arrays
    rtimes(removeIndices) = [];
    robs(:,removeIndices) = [];
    rcalc(:,removeIndices) = [];
    rlat(removeIndices) = [];
    rlon(removeIndices) = [];
    satzen(removeIndices) = []; 
    solzen(removeIndices) = [];
    spres(removeIndices) = [];
    stemp(removeIndices) = [];
    ptemp(:,removeIndices) = [];
    gas_1(:,removeIndices) = [];
    gas_3(:,removeIndices) = [];
    iudef1(removeIndices) = []; 
    iudef2(removeIndices) = []; 
    iudef4(removeIndices) = []; 
    nlevs(removeIndices) = [];
    plevs(:,removeIndices) = [];

    outfileName = sprintf('airxbcal_site-%d_%d_era_rad_stats', s, year);
    outfilePath = fullfile(outfileDir, outfileName)
    eval_str = sprintf(['save %s rtimes robs rcalc rlat rlon *zen ' ...
                        'gas_* iudef* *levs *temp spres'], outfilePath);
    eval(eval_str);
end

