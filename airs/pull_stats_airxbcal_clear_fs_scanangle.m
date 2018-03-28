function pull_stats_airxbcal_clear_fs_scanangle(year, filter, cfg);

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

addpath /asl/matlib/h4tools
addpath /asl/rtp_prod/airs/utils
addpath /asl/packages/rtp_prod2/util
addpath /home/sergio/MATLABCODE/PLOTTER  %
                                         % equal_area_spherical_bands
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

rtpdir = '/asl/rtp/rtp_airxbcal_v5/';
if bCheckConfig & isfield(cfg, 'rtpdir')
    rtpdir = cfg.rtpdir;
end

statsdir = '/asl/data/stats/airs';
if bCheckConfig & isfield(cfg, 'statsdir')
    statsdir = cfg.statsdir;
end

basedir = fullfile(rtpdir, int2str(year), 'clear');
dayfiles = dir(fullfile(basedir, 'era_airxbcal_day*_clear.rtp'));
fprintf(1,'>>> numfiles = %d\n', length(dayfiles));

iday = 1;
% $$$ for giday = 1:length(dayfiles)
for giday = 18:20
    fprintf(1, '>>> year = %d  :: giday = %d\n', year, giday);
   a = dir(fullfile(basedir,dayfiles(giday).name));
   if a.bytes > 100000
      [h,ha,p,pa] = rtpread(fullfile(basedir,dayfiles(giday).name));
      f = h.vchan;  % AIRS proper frequencies
      
      % sanity check on p.robs1 as read in. (There have been
      % instances where this array is short on the spectral
      % dimension which fails in rad2bt. We trap for this here)
      [ochans, nobs] = size(p.robs1);
      [nchans,~] = size(f);
      if ochans ~= nchans
          fprintf(2, ['**>> ERROR: obs/vchan spectral channel ' ...
                      'mismatch in %s. Bypassing day.\n'], dayfiles(giday).name);
          continue;
      end
            
      switch filter
        case 1
          k = find(p.iudef(4,:) == 68); % descending node (night)
          sDescriptor='desc';
        case 2
          k = find(p.iudef(4,:) == 68 & p.landfrac == 0); % descending
                                                     % node (night) ocean
          sDescriptor='desc_ocean';
        case 3
          k = find(p.iudef(4,:) == 68 & p.landfrac == 1); % descending node
                                                        % (night), land
          sDescriptor='desc_land';
        case 4
          k = find(p.iudef(4,:) == 65); % ascending node (day)
          sDescriptor='asc';
        case 5
          k = find(p.iudef(4,:) == 65 & p.landfrac == 0); % ascending node
                                                         % (day), ocean
          sDescriptor='asc_ocean';
        case 6
          k = find(p.iudef(4,:) == 65 & p.landfrac == 1); % ascending node
                                                        % (day), land
          sDescriptor='asc_land';
      end

      pp = rtp_sub_prof(p, k);
      p = pp;
      k = find(abs(p.rlat) < 60);
      pp = rtp_sub_prof(p, k);
      clear p;
      
      if bRunKlayers
          % klayers kills previous sarta in the rtp structures so
          % we need to save values and re-insert after klayers
          % finishes
          rcalc = pp.rcalc;
          tcc = pp.tcc;
          
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
          pp.rcalc = rcalc;
          pp.tcc = tcc;
          clear rcalc tcc;
          
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
         pp.rcalc(i,k) = NaN;
         count_all(i,k) = 0;
      end

      % Loop over latitude bins
      nfors = 90; 
      for ifor = 1:nfors
          % subset based on latitude bin
          inbin = find(pp.xtrack == ifor);
          p = rtp_sub_prof(pp,inbin);
          bincount = count_all(:,inbin); 
          binwater = mmwater(inbin);
          
          % Radiance mean and std
          r  = p.robs1;
          clr_calc = p.rcalc;
          
          % spectral
          robs(iday,ifor,:) = nanmean(r,2);
          rclr(iday,ifor,:) = nanmean(clr_calc,2);
          rclrbias_std(iday,ifor,:) = nanstd(r-clr_calc,0,2);
          
          lat_mean(iday,ifor) = nanmean(p.rlat);
          lon_mean(iday,ifor) = nanmean(p.rlon);
          scanang_mean(iday,ifor) = nanmean(p.scanang);
          solzen_mean(iday,ifor) = nanmean(p.solzen);
          rtime_mean(iday,ifor)  = nanmean(p.rtime);
          count(iday,ifor,:) = sum(bincount,2)';
          tcc_mean(iday, ifor) = nanmean(p.tcc);
          stemp_mean(iday,ifor) = nanmean(p.stemp);
          ptemp_mean(iday,ifor,:) = nanmean(p.ptemp,2);
          gas1_mean(iday,ifor,:) = nanmean(p.gas_1,2);
          gas3_mean(iday,ifor,:) = nanmean(p.gas_3,2);
          spres_mean(iday,ifor) = nanmean(p.spres);
          nlevs_mean(iday,ifor) = nanmean(p.nlevs);
          iudef4_mean(iday,ifor) = nanmean(p.iudef(4,:));
          mmwater_mean(iday,ifor) = nanmean(binwater);
          satzen_mean(iday,ifor) = nanmean(p.satzen);
          plevs_mean(iday,ifor,:) = nanmean(p.plevs,2);
      end  % end loop over latitudes
          iday = iday + 1
   end % if a.bytes > 1000000
end  % giday

outfile = fullfile(statsdir, sprintf('rtp_airxbcal_era_rad_kl_018-020_%s_random_fs_scanang_%s', ...
           int2str(year), sDescriptor));
eval_str = ['save ' outfile [' robs rcl* *_mean count* ' ...
                    'trace']];
eval(eval_str);
