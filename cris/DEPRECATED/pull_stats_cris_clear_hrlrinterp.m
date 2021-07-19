function pull_stats_cris_clear_hrlrinterp(year, filter);

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
addpath ~/git/rtp_prod2/util
addpath /asl/rtp_prod/cris/unapod
addpath /home/sbuczko1/git/pull_stats/util  %
                                         % equal_area_spherical_bands
addpath /home/sbuczko1/git/rtp_prod2/cris  % cris_lowres_chans
addpath /asl/matlib/rtptools  % mmwater_rtp.m

[sID, sTempPath] = genscratchpath();
sTempPath = '/home/sbuczko1/Work/scratch';

klayers_exec = ['/asl/packages/klayersV205/BinV201/' ...
                'klayers_airs_wetwater'];

% Get proper frequencies for these data
[n1,n2,n3,userLW,userMW,userSW, ichan] = cris_lowres_chans();
f = cris_vchan(2, userLW, userMW, userSW);
nchans = length(f);

basedir = fullfile('/asl/rtp/rtp_cris_ccast_hires/clear', ...
                   int2str(year));
% cris_era_csarta_clear_d20180422.rtp
dayfiles = dir(fullfile(basedir, 'cris_era_csarta_clear_d*.rtp'));
ndays = length(dayfiles);

% calculate latitude bins
nbins=20; % gives 2N+1 element array of lat bin boundaries
latbins = equal_area_spherical_bands(nbins);
nlatbins = length(latbins);

nfovs = 9;
nlevs = 101;

% allocate final accumulator arrays
robs = nan(ndays, nlatbins, nfovs, nchans);
rclr = nan(ndays, nlatbins, nfovs, nchans);
rbias_std = nan(ndays, nlatbins, nfovs, nchans);

lat_mean = nan(ndays, nlatbins, nfovs);
lon_mean = nan(ndays, nlatbins, nfovs);
solzen_mean = nan(ndays, nlatbins, nfovs);
rtime_mean = nan(ndays, nlatbins, nfovs); 
count = nan(ndays, nlatbins, nfovs, nchans);
stemp_mean = nan(ndays, nlatbins, nfovs);
ptemp_mean = nan(ndays, nlatbins, nfovs, nlevs);
gas1_mean = nan(ndays, nlatbins, nfovs, nlevs);
gas3_mean = nan(ndays, nlatbins, nfovs, nlevs);
spres_mean = nan(ndays, nlatbins, nfovs);
nlevs_mean = nan(ndays, nlatbins, nfovs);
iudef4_mean = nan(ndays, nlatbins, nfovs);
mmwater_mean = nan(ndays, nlatbins, nfovs);
satzen_mean = nan(ndays, nlatbins, nfovs);
satazi_mean = nan(ndays, nlatbins, nfovs);
plevs_mean = nan(ndays, nlatbins, nfovs, nlevs);


iday = 1;
% for giday = 1:50:length(dayfiles)
for giday = 1:length(dayfiles)
   fprintf(1, '>>> year = %d  :: giday = %d\n', year, giday);
   a = dir(fullfile(basedir,dayfiles(giday).name));
   if a.bytes > 100000
      [h,ha,p,pa] = rtpread(fullfile(basedir,dayfiles(giday).name));

      %%%%%%% interp hires robs/rclr down to lowres
      addpath /asl/packages/ccast/source
      
      % specify the target res
      opt1 = struct;
      opt1.user_res = 'lowres';  % 'lowres' or 'midres'
      opt1.inst_res = 'hires3';  % nominal inst res
      wlaser = 773.1301;         % nominal wlaser

      % first robs
      r1LW = p.robs1(1:717,:);
      v1LW = h.vchan(1:717);
      
      r1MW = p.robs1(718:1586,:);
      v1MW = h.vchan(718:1586);
      % interpolate the MW user grid
      [instMW, userMW] = inst_params('MW', wlaser, opt1);
      [r2MW, v2MW] = finterp(r1MW(:,:), v1MW, userMW.dv);
      ix = find(userMW.v1 <= v2MW & v2MW <= userMW.v2);
      v2MW = v2MW(ix);
      r2MW = r2MW(ix, :); 
      if userMW.v1 ~= v2MW(1) | userMW.v2 ~= v2MW(end)
          error('MW grid mismatch')
      end


      r1SW = p.robs1(1587:2223,:);
      v1SW = h.vchan(1587:2223);
      % interpolate the SW user grid
      [instSW, userSW] = inst_params('SW', wlaser, opt1);
      [r2SW, v2SW] = finterp(r1SW(:,:), v1SW, userSW.dv);
      ix = find(userSW.v1 <= v2SW & v2SW <= userSW.v2);
      v2SW = v2SW(ix);
      r2SW = r2SW(ix, :); 
      if userSW.v1 ~= v2SW(1) | userSW.v2 ~= v2SW(end)
          error('SW grid mismatch')
      end
      
      p.robs1  = cat(1, r1LW,r2MW,r2SW);

      % now rclr
      r1LW = p.rclr(1:717,:);
      v1LW = h.vchan(1:717);
      
      r1MW = p.rclr(718:1586,:);
      v1MW = h.vchan(718:1586);
      % interpolate the MW user grid
      [instMW, userMW] = inst_params('MW', wlaser, opt1);
      [r2MW, v2MW] = finterp(r1MW(:,:), v1MW, userMW.dv);
      ix = find(userMW.v1 <= v2MW & v2MW <= userMW.v2);
      v2MW = v2MW(ix);
      r2MW = r2MW(ix, :); 
      if userMW.v1 ~= v2MW(1) | userMW.v2 ~= v2MW(end)
          error('MW grid mismatch')
      end


      r1SW = p.rclr(1587:2223,:);
      v1SW = h.vchan(1587:2223);
      % interpolate the SW user grid
      [instSW, userSW] = inst_params('SW', wlaser, opt1);
      [r2SW, v2SW] = finterp(r1SW(:,:), v1SW, userSW.dv);
      ix = find(userSW.v1 <= v2SW & v2SW <= userSW.v2);
      v2SW = v2SW(ix);
      r2SW = r2SW(ix, :); 
      if userSW.v1 ~= v2SW(1) | userSW.v2 ~= v2SW(end)
          error('SW grid mismatch')
      end
      
      p.rclr  = cat(1, r1LW,r2MW,r2SW);

      % clean up 
      clear r1LW r1MW r1SW v1LW v1MW v1SW r2LW r2MW r2SW v2LW v2MW ...
          v2SW ix userSW userMW opt1 wlaser
      %%%%%%%

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

      if h.ptype == 0
          % run klayers on the rtp data to convert levels -> layers
          fprintf(1, '>>> running klayers... ');
          % first remove rcalc field and save it for later restore
          tmp_rclr = pp.rclr;
          pp = rmfield(pp, 'rclr');

          fn_rtp1 = fullfile(sTempPath, ['cris_' sID '_1.rtp']);
          fprintf(1, '>>>> writing temp file %s\n', fn_rtp1);
          rtpwrite(fn_rtp1, h,ha,pp,pa);
          clear pp;
          fn_rtp2 = fullfile(sTempPath, ['cris_' sID '_2.rtp']);
          fprintf(1, '>>>> klayers output to %s\n', fn_rtp2);
          klayers_run = [klayers_exec ' fin=' fn_rtp1 ' fout=' fn_rtp2 ...
                         ' > ' sTempPath '/kout.txt'];
          unix(klayers_run);

          % Read klayers output into local rtp variables
          [h,ha,pp,pa] = rtpread(fn_rtp2);

          % restore rcalc
          pp.rclr = tmp_rclr;
          clear tmp_rclr;

          fprintf(1, 'Done\n');
      end
      
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
      count_all = int8(ones(nchans, nobs, nfovs));
      
      % loop over latitude bins
      for ilat = 1:nlatbins-1
          % subset based on latitude bin
          inbin = find(pp.rlat > latbins(ilat) & pp.rlat <= ...
                       latbins(ilat+1));
          p = rtp_sub_prof(pp,inbin);
          
          for z = 1:9  % loop over FOVs to further sub-select
              ifov = find(p.ifov == z);
              p2 = rtp_sub_prof(p, ifov);
              
              bincount = count_all(:,inbin,z); 
              binwater = mmwater(inbin);
              % Loop over obs in day
              % Radiance mean and std
              
              r  = p2.robs1;
              rc = p2.rclr;
              
              % Convert r to rham
              r = box_to_ham(r);  % assumes r in freq order!!  Needed
                                  % for lowres
              robs(iday,ilat,z,:) = nanmean(r,2);
              rcal(iday,ilat,z,:) = nanmean(rc,2);
              rbias_std(iday, ilat,z,:) = nanstd(r-rc,0,2);
              
              lat_mean(iday,ilat,z) = nanmean(p2.rlat);
              lon_mean(iday,ilat,z) = nanmean(p2.rlon);
              solzen_mean(iday,ilat,z) = nanmean(p2.solzen);
              rtime_mean(iday,ilat,z)  = nanmean(p2.rtime);
              count(iday,ilat,z) = sum(bincount(1,:))';
              stemp_mean(iday,ilat,z) = nanmean(p2.stemp);
              iudef4_mean(iday,ilat,z) = nanmean(p2.iudef(4,:));
              ptemp_mean(iday,ilat,z,:) = nanmean(p2.ptemp,2);
              gas1_mean(iday,ilat,z,:) = nanmean(p2.gas_1,2);
              gas3_mean(iday,ilat,z,:) = nanmean(p2.gas_3,2);
              spres_mean(iday,ilat,z) = nanmean(p2.spres);
              nlevs_mean(iday,ilat,z) = nanmean(p2.nlevs);
              satzen_mean(iday,ilat,z) = nanmean(p2.satzen);
              satazi_mean(iday,ilat,z) = nanmean(p2.satazi);
              plevs_mean(iday,ilat,z,:) = nanmean(p2.plevs,2);
              mmwater_mean(iday,ilat) = nanmean(binwater);
% $$$               scanang_mean(iday,ilat,z) = nanmean(p.scanang);
          end  % ifov (z)
      end  % end loop over ilat
          
          iday = iday + 1
   end % if a.bytes > 1000000
end  % giday
eval_str = ['save /home/sbuczko1/Work/cris-interp-testing/rtp_cris_lowres_rad_'  int2str(year) ...
            '_clear' sDescriptor ' robs rcal *_std *_mean count '];

eval(eval_str);
