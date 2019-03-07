function pull_stats_airicrad_clear_scanangle_latbin(year, filter, ...
                                                  ilat, cfg);
% PULL_STATS_AIRICRAD_RAND Create stats accumulations from rtp
%

addpath /asl/matlib/h4tools
addpath /asl/rtp_prod/airs/utils
addpath /asl/packages/rtp_prod2/util
addpath /home/sergio/MATLABCODE/PLOTTER  %
                                         % equal_area_spherical_bands
addpath /asl/matlib/rtptools  % mmwater_rtp.m

[sID, sTempPath] = genscratchpath();

% check for existence of configuration struct
bCheckConfig = false;
if nargin == 4
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

rtpdir = '/asl/rtp/rtp_airicrad_v6/';
if bCheckConfig & isfield(cfg, 'rtpdir')
    rtpdir = cfg.rtpdir;
end

statsdir = '/asl/data/stats/airs/clear';
if bCheckConfig & isfield(cfg, 'statsdir')
    statsdir = cfg.statsdir;
end

basedir = fullfile(rtpdir, 'clear', int2str(year));
dayfiles = dir(fullfile(basedir, 'era_airicrad_day*_clear.rtp'));
ndays = length(dayfiles);
% $$$ ndays = 32;
fprintf(1,'>>> numfiles = %d\n', ndays);

% calculate latitude bins
nbins=20; % gives 2N+1 element array of lat bin boundaries
latbinedges = equal_area_spherical_bands(nbins);
% $$$ nlatbins = length(latbinedges)-1;
nlatbins = 1;

nchans = 2645;  % AIRICRAD/L1C channel space
nlevs = 101;  % klayers output
nfovs = 90; % FOVs

% allocate final accumulator arrays
robs = zeros(nlatbins, nfovs, nchans);
rclr = zeros(nlatbins, nfovs, nchans);

rclrbias_std = zeros(nlatbins, nfovs, nchans);

l1cproc_mean = zeros(nlatbins, nfovs, nchans);
l1csreason_mean = zeros(nlatbins, nfovs, nchans);
dbtun_mean = zeros(nlatbins, nfovs);

lat_mean = zeros(nlatbins, nfovs);
lon_mean = zeros(nlatbins, nfovs);
solzen_mean = zeros(nlatbins, nfovs);
rtime_mean = zeros(nlatbins, nfovs); 
count = zeros(nlatbins, nfovs, nchans);
tcc_mean = zeros(nlatbins, nfovs);
stemp_mean = zeros(nlatbins, nfovs);
ptemp_mean = zeros(nlatbins, nfovs, nlevs);
gas1_mean = zeros(nlatbins, nfovs, nlevs);
gas3_mean = zeros(nlatbins, nfovs, nlevs);
spres_mean = zeros(nlatbins, nfovs);
nlevs_mean = zeros(nlatbins, nfovs);
iudef4_mean = zeros(nlatbins, nfovs);
mmwater_mean = zeros(nlatbins, nfovs);
satzen_mean = zeros(nlatbins, nfovs);
satazi_mean = zeros(nlatbins, nfovs);
plevs_mean = zeros(nlatbins, nfovs, nlevs);

iday = 1;
% $$$ for giday = 60:91
for giday = 1:length(dayfiles)
   fprintf(1, '>>> year = %d  :: giday = %d\n', year, giday);
   a = dir(fullfile(basedir,dayfiles(giday).name));
   if a.bytes < 100000
        fprintf(2, '**>> ERROR: short input rtp file %s\n', dayfiles(giday).name); 
        continue;
   end
       
      [h,ha,p,pa] = rtpread(fullfile(basedir,dayfiles(giday).name));
      f = h.vchan;  % AIRS proper frequencies
      
      % sanity check on p.robs1 as read in. (There have been
      % instances where this array is short on the spectral
      % dimension which fails in rad2bt. We trap for this here)
      obs = size(p.robs1);
      chans = size(f);
      if obs(1) ~= chans(1)
          fprintf(2, ['**>> ERROR: obs/vchan spectral channel ' ...
                      'mismatch in %s. Bypassing day.\n'], dayfiles(giday).name);
          continue;
      end
            
      switch filter
        case 1
          k = find(p.iudef(4,:) == 1); % descending node (night)
          sDescriptor='desc';
        case 2
          k = find(p.iudef(4,:) == 1 & p.landfrac == 0); % descending
                                                     % node (night) ocean
          sDescriptor='desc_ocean';
        case 3
          k = find(p.iudef(4,:) == 1 & p.landfrac == 1); % descending node
                                                        % (night), land
          sDescriptor='desc_land';
        case 4
          k = find(p.iudef(4,:) == 0); % ascending node (day)
          sDescriptor='asc';
        case 5
          k = find(p.iudef(4,:) == 0 & p.landfrac == 0); % ascending node
                                                         % (day), ocean
          sDescriptor='asc_ocean';
        case 6
          k = find(p.iudef(4,:) == 0 & p.landfrac == 1); % ascending node
                                                        % (day), land
          sDescriptor='asc_land';
      end

      % also find the required latbin
      inbin = find(p.rlat > latbinedges(ilat) & p.rlat <= ...
                   latbinedges(ilat+1));

      filteredobs = intersect(k, inbin);
      if length(filteredobs) == 0
          continue;
      end
      prof = rtp_sub_prof(p, filteredobs);
      clear p;

      % check for empty profile struct after subsetting
      if length(prof.rtime) == 0
          fprintf(2, ['>>> No obs found after filter subset for %s. ' ...
                      'Continuing to next day.\n'], dayfiles(giday).name);
          continue;  % jump to next day
      end

      % concatenate rtp structs
      if ~exist('pp')
          pp = prof;
      else
          % build structure descriptors
          fnames = fieldnames(prof);
          nnames = length(fnames);
          % now the concatenation
          for j=1:nnames
              pp.(fnames{j}) = [pp.(fnames{j}) prof.(fnames{j})];
          end

      end
      clear prof;
end

      if bRunKlayers
          % klayers kills previous sarta in the rtp structures so
          % we need to save values and re-insert after klayers
          % finishes
          tmp_rclr = pp.rclr;
          tmp_tcc = pp.tcc;
          tmp_l1cproc = pp.l1cproc;
          tmp_l1csreason = pp.l1csreason;
          tmp_dbtun = pp.dbtun;
          
          % run klayers on the rtp data to convert levels -> layers
          fprintf(1, '>>> running klayers... ');
          fn_rtp1 = fullfile(sTempPath, ['airs_' sID '_1.rtp']);
          rtpwrite(fn_rtp1, h,ha,pp,pa);
          clear pp;
          fn_rtp2 = fullfile(sTempPath, ['airs_' sID '_2.rtp']);
          klayers_run = [klayers_exec ' fin=' fn_rtp1 ' fout=' fn_rtp2 ...
                         ' > ' sTempPath '/kout.txt'];
          unix(klayers_run);
          fprintf(1, 'Done\n');

          % Read klayers output into local rtp variables
          [h,ha,pp,pa] = rtpread(fn_rtp2);
          
          f = h.vchan;  % AIRS proper frequencies

          % restore sarta values
          pp.rclr = tmp_rclr;
          pp.tcc = tmp_tcc;
          pp.l1cproc = tmp_l1cproc;
          pp.l1csreason = tmp_l1csreason;
          pp.dbtun = tmp_dbtun;
          clear tmp_rclr tmp_tcc tmp_l1cproc tmp_l1csreason tmp_dbtun;
          
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
              pp.palts(badlayers, i) = NaN;
              pp.gas_1(badlayers, i) = NaN;
              pp.gas_2(badlayers, i) = NaN;
              pp.gas_3(badlayers, i) = NaN;
              pp.gas_4(badlayers, i) = NaN;
              pp.gas_5(badlayers, i) = NaN;
              pp.gas_6(badlayers, i) = NaN;
              pp.gas_12(badlayers, i) = NaN;          
              pp.ptemp(badlayers, i) = NaN;
          end
          
      end

      % Initialize counts
      n = length(pp.rlat);
      count_all = ones(nchans,n);
      for i=1:nchans
         % Find bad channels (l1c incorporates l1b calnum by
         % inserting interpolated values for any place where l1b
         % would be 'bad'. i.e. there should be no 'bad' channels)
         k = find( pp.robs1(i,:) == -9999);
%          % These are the good channels
%          kg = setdiff(1:n,k);
% NaN's for bad channels
         pp.robs1(i,k) = NaN;
         pp.rclr(i,k) = NaN;
         count_all(i,k) = 0;
      end

      % loop over FOVs
      for ifov = 1:90
          infov = find(pp.xtrack == ifov);
          if length(infov) == 0
              continue;
          end
          bincount = count_all(:,infov); 
          binwater = mmwater(infov);

          % Radiance mean and std
          r  = pp.robs1(:,infov);
          clr_calc = pp.rclr(:,infov);
          
          % spectral
          robs(1,ifov,:) = nanmean(r,2);
          rclr(1,ifov,:) = nanmean(clr_calc,2);
          rclrbias_std(1,ifov,:) = nanstd(r-clr_calc,0,2);

          l1cproc_mean(1,ifov, :) = nanmean(pp.l1cproc(:,infov), 2);
          l1csreason_mean(1,ifov, :) = nanmean(pp.l1csreason(:,infov), ...
                                                   2);
          dbtun_mean(1,ifov) = nanmean(pp.dbtun(infov));
          
          lat_mean(1,ifov) = nanmean(pp.rlat(infov));
          lon_mean(1,ifov) = nanmean(pp.rlon(infov));
          solzen_mean(1,ifov) = nanmean(pp.solzen(infov));
          rtime_mean(1,ifov)  = nanmean(pp.rtime(infov));
          count(1,ifov,:) = sum(bincount,2)';
          tcc_mean(iday, 1) = nanmean(pp.tcc(infov));
          stemp_mean(1,ifov) = nanmean(pp.stemp(infov));
          ptemp_mean(1,ifov,:) = nanmean(pp.ptemp(:,infov),2);
          gas1_mean(1,ifov,:) = nanmean(pp.gas_1(:,infov),2);
          gas3_mean(1,ifov,:) = nanmean(pp.gas_3(:,infov),2);
          spres_mean(1,ifov) = nanmean(pp.spres(infov));
          nlevs_mean(1,ifov) = nanmean(pp.nlevs(infov));
          iudef4_mean(1,ifov) = nanmean(pp.iudef(4,infov));
          mmwater_mean(1,ifov) = nanmean(binwater);
          satzen_mean(1,ifov) = nanmean(pp.satzen(infov));
          satazi_mean(1,ifov) = nanmean(pp.satazi(infov));          
          plevs_mean(1,ifov,:) = nanmean(pp.plevs(:,infov),2);
      end  % end loop over ifov

outfile = fullfile(statsdir, sprintf('rtp_airicrad_era_rad_kl_1year_scanangle_lbin_%d_%s_clear_%s', ...
           ilat, int2str(year), sDescriptor));
eval_str = ['save ' outfile [' robs rcl* *_mean count latbinedges ' ...
                    'trace']];
fprintf(1,'>> Executing save command: %s\n', eval_str);
eval(eval_str);
