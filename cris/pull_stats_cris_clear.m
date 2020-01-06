function pull_stats_cris_clear(year, filter);

%**************************************************
% need to make this work on daily concat files: look for loop over
% granules, this will need to be removed. Also break out by fov
% (add loop and index over p.ifov)
%
% following the 'file in; file out' standard, this routine will
% read in ONE daily concatenated rtp file and calculate its
% stats. There will be a driver function above this that will feed
% rtp file paths to this routine and can provide for slurm indexing
% and chunking
%**************************************************

%year = 2014;

addpath /asl/matlib/h4tools
addpath /asl/rtp_prod/airs/utils
addpath /asl/packages/ccast/motmsc/utils/
addpath /home/sbuczko1/git/pull_stats_DEV/util
addpath /asl/rtp_prod/cris/unapod
addpath /home/sbuczko1/git/pull_stats_DEV/util  %
                                         % equal_area_spherical_bands
addpath /home/sbuczko1/git/pull_stats_DEV/cris/util  % cris_lowres_chans
addpath /asl/matlib/rtptools  % mmwater_rtp.m

[sID, sTempPath] = genscratchpath();

% check for existence of configuration struct
bCheckConfig = false;
if nargin == 3
    bCheckConfig = true;
end

% record run start datetime in output stats file for tracking
trace.RunDate = datetime('now','TimeZone','local','Format', ...
                         'd-MMM-y HH:mm:ss Z');
trace.Reason = 'Normal pull_stats runs';
if bCheckConfig & isfield(cfg, 'Reason')
    trace.Reason = cfg.Reason;
end

bRunKlayers = true;
klayers_exec = ['/asl/packages/klayersV205/BinV201/' ...
                'klayers_airs_wetwater'];
if bCheckConfig & isfield(cfg, 'klayers') & cfg.klayers == false
    bRunKlayers = false;
end
trace.klayers = bRunKlayers;

trace.droplayers = false;
if bCheckConfig & isfield(cfg, 'droplayers') & cfg.droplayers == true
    trace.droplayers = true;
end

rtpdir = '/asl/rtp/rtp_cris_ccast_lowres/';
if bCheckConfig & isfield(cfg, 'rtpdir')
    rtpdir = cfg.rtpdir;
end

statsdir = '/asl/data/stats/cris';
if bCheckConfig & isfield(cfg, 'statsdir')
    statsdir = cfg.statsdir;
end

basedir = fullfile(rtpdir, 'clear', int2str(year));
dayfiles = dir(fullfile(basedir, 'cris_era_csarta_clear_d*.rtp'));
ndays = length(dayfiles);
fprintf(1,'>>> numfiles = %d\n', ndays);

% Get proper frequencies for these data
[n1,n2,n3,userLW,userMW,userSW, ichan] = cris_lowres_chans();
f = cris_vchan(2, userLW, userMW, userSW);
nchans = length(f);

% calculate latitude bins
nbins=20; % gives 2N+1 element array of lat bin boundaries
latbin_edges = equal_area_spherical_bands(nbins);
nlatbins = length(latbin_edges)-1;

nlevs = 101;  % klayers output
nfovs = 9;    % FOVs/FOR

% allocate final accumulator arrays
robs = zeros(ndays, nlatbins, nfovs, nchans);
rclr = zeros(ndays, nlatbins, nfovs, nchans);
rclrbias_std = zeros(ndays, nlatbins, nfovs, nchans);

lat_mean = zeros(ndays, nlatbins, nfovs);
lon_mean = zeros(ndays, nlatbins, nfovs);
solzen_mean = zeros(ndays, nlatbins, nfovs);
rtime_mean = zeros(ndays, nlatbins, nfovs); 
count = zeros(ndays, nlatbins, nfovs, nchans);
tcc_mean = zeros(ndays, nlatbins, nfovs);
stemp_mean = zeros(ndays, nlatbins, nfovs);
ptemp_mean = zeros(ndays, nlatbins, nfovs, nlevs);
gas1_mean = zeros(ndays, nlatbins, nfovs, nlevs);
gas3_mean = zeros(ndays, nlatbins, nfovs, nlevs);
spres_mean = zeros(ndays, nlatbins, nfovs);
nlevs_mean = zeros(ndays, nlatbins, nfovs);
iudef4_mean = zeros(ndays, nlatbins, nfovs);
mmwater_mean = zeros(ndays, nlatbins, nfovs);
satzen_mean = zeros(ndays, nlatbins, nfovs);
plevs_mean = zeros(ndays, nlatbins, nfovs, nlevs);

iday = 1;
% for giday = 1:50:length(dayfiles)
for giday = 1:length(dayfiles)
   fprintf(1, '>>> year = %d  :: giday = %d\n', year, giday);
   a = dir(fullfile(basedir,dayfiles(giday).name));
   if a.bytes > 100000
      [h,ha,p,pa] = rtpread(fullfile(basedir,dayfiles(giday).name));

      % cris rtps lost valid p.iudef(4,:) somewhere in this last
      % re-processing
      p.iudef(4,:) = p.solzen>90.0;  % set p.iudef(4) == 1 for
                                     % descending/night, 0 for ascending/day

      switch filter
        case 1
          k = find(p.iudef(4,:) == 1); % descending node (night)
          sDescriptor='_desc';
        case 2
          k = find(p.iudef(4,:) == 1 & p.landfrac == 0); % descending node
                                                         % (night), ocean
          sDescriptor='_desc_ocean';
        case 3
          k = find(p.iudef(4,:) == 1 & p.landfrac == 1); % descending node
                                                        % (night), land
          sDescriptor='_desc_land';
        case 4
          k = find(p.iudef(4,:) == 0); % ascending node (day)
          sDescriptor='_asc';
        case 5
          k = find(p.iudef(4,:) == 0 & p.landfrac == 0); % ascending node
                                                         % (day), ocean
          sDescriptor='_asc_ocean';
        case 6
          k = find(p.iudef(4,:) == 0 & p.landfrac == 1); % ascending node
                                                        % (day), land
          sDescriptor='_asc_land';
      end

      pp = rtp_sub_prof(p, k);

      % run klayers on the rtp data to convert levels -> layers
      % save calcs as the re-run of klayers wipes them out
      tmp_rclr = pp.rclr;
      tmp_tcc = pp.tcc;
      fprintf(1, '>>> running klayers... ');
      fn_rtp1 = fullfile(sTempPath, ['cris_' sID '_1.rtp']);
      rtpwrite(fn_rtp1, h,ha,pp,pa);
      clear pp;
      fn_rtp2 = fullfile(sTempPath, ['cris_' sID '_2.rtp']);
      klayers_run = [klayers_exec ' fin=' fn_rtp1 ' fout=' fn_rtp2 ...
                     ' > ' sTempPath '/kout.txt'];
      unix(klayers_run);
      fprintf(1, 'Done\n');

      % Read klayers output into local rtp variables
      [h,ha,pp,pa] = rtpread(fn_rtp2);
      % restore rclr
      pp.rclr = tmp_rclr;
      pp.tcc = tmp_tcc;
      clear tmp_rclr tmp_rcld tmp_tcc;
      
      % get column water
      mmwater = mmwater_rtp(h, pp);

      % Check for obs with layer profiles that go lower than
      % topography. Need to check nlevs and NaN out any layers
      % at or below this level

      % ** Any layers-sensitive variables added in averaging code below must
      % ** be checked here first.
      for i=1:length(pp.nlevs)
          badlayers = pp.nlevs(i) : 101;
          pp.plevs(badlayers, i) = NaN;
          pp.gas_1(badlayers, i) = NaN;
          pp.gas_3(badlayers, i) = NaN;
          pp.ptemp(badlayers, i) = NaN;
      end


      % initialize counts and look for bad channels (what should
      % the iasi bad channel test look like?)
      [nchans, nobs] = size(pp.robs1);
      
      % loop over latitude bins
      for ilat = 1:nlatbins-1
          % subset based on latitude bin
          inbin = find(pp.rlat > latbin_edges(ilat) & pp.rlat <= ...
                       latbin_edges(ilat+1));
          p = rtp_sub_prof(pp,inbin);
          
          for z = 1:9  % loop over FOVs to further sub-select
              infov = find(p.ifov == z);
              p2 = rtp_sub_prof(p, infov);
              
              count_infov = ones(length(infov), nchans);
              binwater = mmwater(infov);
              % Loop over obs in day
              % Radiance mean and std
              
              r  = p2.robs1;
              rc = p2.rclr;
              
              % Convert r to rham
              r = box_to_ham(r);  % assumes r in freq order!!  Needed
                                  % for lowres
              robs(iday,ilat,z,:) = nanmean(r,2);
              rclr(iday,ilat,z,:) = nanmean(rc,2);
              rclrbias_std(iday, ilat, z, :) = nanstd(r-rc,0,2);
              
              lat_mean(iday,ilat,z) = nanmean(p2.rlat);
              lon_mean(iday,ilat,z) = nanmean(p2.rlon);
              solzen_mean(iday,ilat,z) = nanmean(p2.solzen);
              rtime_mean(iday,ilat,z)  = nanmean(p2.rtime);
              count(iday,ilat,z,:) = sum(count_infov);
              stemp_mean(iday,ilat,z) = nanmean(p2.stemp);
              iudef4_mean(iday,ilat,z) = nanmean(p2.iudef(4,:));
              ptemp_mean(iday,ilat,z,:) = nanmean(p2.ptemp,2);
              gas1_mean(iday,ilat,z,:) = nanmean(p2.gas_1,2);
              gas3_mean(iday,ilat,z,:) = nanmean(p2.gas_3,2);
              spres_mean(iday,ilat,z) = nanmean(p2.spres);
              nlevs_mean(iday,ilat,z) = nanmean(p2.nlevs);
              satzen_mean(iday,ilat,z) = nanmean(p2.satzen);
              plevs_mean(iday,ilat,z,:) = nanmean(p2.plevs,2);
              mmwater_mean(iday,ilat) = nanmean(binwater);
          end  % infov (z)
      end  % end loop over ilat
          
          iday = iday + 1
   end % if a.bytes > 1000000
end  % giday
eval_str = ['save /asl/data/stats/cris/rtp_cris_lowres_rad_'  int2str(year) ...
            '_clear' sDescriptor ' robs rclr *_std *_mean count '];
eval(eval_str);
