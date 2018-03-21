function pull_stats_airicrad_rand_fs_scanangle(year, filter, cfg);

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

nchans = 2645;  % AIRICRAD/L1C channel space

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

rtpdir = '/asl/rtp/rtp_airicrad_v6/';
if bCheckConfig & isfield(cfg, 'rtpdir')
    rtpdir = cfg.rtpdir;
end

statsdir = '/asl/data/stats/airs';
if bCheckConfig & isfield(cfg, 'statsdir')
    statsdir = cfg.statsdir;
end

llat = 60;
if bCheckConfig & isfield(cfg, 'latlimit')
    llat = cfg.latlimit;
end

basedir = fullfile(rtpdir, int2str(year), 'random_fs');
dayfiles = dir(fullfile(basedir, 'era_airicrad_day*_random_fs.rtp'));
fprintf(1,'>>> numfiles = %d\n', length(dayfiles));



% calculate latitude bins
% $$$ nbins=20; % gives 2N+1 element array of lat bin boundaries
% $$$ latbins = equal_area_spherical_bands(nbins);
% $$$ nlatbins = length(latbins);

iday = 1;
% $$$ for giday = 1:100:length(dayfiles)
for giday = 1:length(dayfiles)
   fprintf(1, '>>> year = %d  :: giday = %d\n', year, giday);
   a = dir(fullfile(basedir,dayfiles(giday).name));
   if a.bytes > 100000
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
      k = find(abs(p.rlat) <= llat);
      pp = rtp_sub_prof(p, k);
      clear p;
      
      if bRunKlayers
          % klayers kills previous sarta in the rtp structures so
          % we need to save values and re-insert after klayers
          % finishes
          rcalc = pp.rcalc;
          sarta_rclearcalc = pp.sarta_rclearcalc;
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
          pp.sarta_rclearcalc = sarta_rclearcalc;
          pp.tcc = tcc;
          clear rcalc sarta_rclearcalc tcc;
          
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

          % Remove 'clouds' that have only partial defining
          % characteristics (effects 0.2% of obs in test rtp data)
          % cfrac=cngwat=0 but cpsize/cprtop/cprbot ~= 0
          a = find(p.cfrac==0 & p.cngwat==0 & p.cprtop>0);
          p.ctype(a)=-1;
          a = find(p.cfrac2==0 & p.cngwat2==0 & p.cprtop2>0);
          p.ctype2(a)=-1;

          % replace 'bad data' placeholders in cprtop/cprbot (-9999
          % -> NaN, esesntially except for 'clouds' from above)
          a = find(p.ctype==-1);
          p.cprtop(a)=NaN; p.cprbot(a)=NaN;
          a = find(p.ctype2==-1);
          p.cprtop2(a)=NaN; p.cprbot2(a)=NaN;

          % massage situations where we have two of same cloud type
          % in a single obs i.e. p.ctype=p.ctype2 (average to a
          % single cloud)
          dblwat = find(p.ctype==101 & p.ctype2==101);
          p.cfrac2(dblwat) = p.cfrac(dblwat) + p.cfrac2(dblwat);
          p.cpsize2(dblwat) = (p.cpsize(dblwat) + ...
                               p.cpsize2(dblwat))/2;
          p.cprtop2(dblwat) = (p.cprtop(dblwat) + ...
                               p.cprtop2(dblwat))/2;
          p.cprbot2(dblwat) = (p.cprbot(dblwat) + ...
                               p.cprbot2(dblwat))/2;
          p.cngwat2(dblwat) = (p.cfrac(dblwat).*p.cngwat(dblwat) + ...
                               p.cfrac2(dblwat).*p.cngwat2(dblwat))./ ...
                               p.cfrac2(dblwat);
          p.cfrac(dblwat)=0; p.cpsize(dblwat)=0;
          p.cngwat(dblwat)= 0; p.cprtop(dblwat)=NaN;
          p.cprbot(dblwat)=NaN; p.ctype(dblwat)=-1;
          
          dblice = find(p.ctype==201 & p.ctype2==201);
          p.cfrac(dblice) = p.cfrac(dblice) + p.cfrac2(dblice);
          p.cpsize(dblice) = (p.cpsize(dblice) + ...
                               p.cpsize2(dblice))/2;
          p.cprtop(dblice) = (p.cprtop(dblice) + ...
                               p.cprtop2(dblice))/2;
          p.cprbot(dblice) = (p.cprbot(dblice) + ...
                               p.cprbot2(dblice))/2;
          p.cngwat(dblice) = (p.cfrac(dblice).*p.cngwat(dblice) + ...
                               p.cfrac2(dblice).*p.cngwat2(dblice))./ ...
                               p.cfrac2(dblice);
          p.cfrac2(dblice)=0; p.cpsize2(dblice)=0;
          p.cngwat2(dblice)= 0; p.cprtop2(dblice)=NaN;
          p.cprbot2(dblice)=NaN; p.ctype2(dblice)=-1;
          
          % at this point cloud fields should now be set so that
          % ctype/cpsize/cngwat and similar contain only ice clouds
          % (or none) and ctype2/cpsize2/cngwat2 and similar
          % contain only water clouds (or none)
          trace.NOTE_ON_CLOUD_VARS = ['ctype_mean and similar are average ' ...
                              'ICE cloud parameters, ctype2_mean and similar ' ...
                              'are average WATER cloud parameters'];

          % pull out cloud parameters for averaging
          % build logical array for each cloud type and pull out
          % only values associated with a cloud for averaging
          lIce = (p.ctype==201 & p.cfrac>0 & p.cngwat>0 & p.cpsize>0 ...
                   & p.cprbot>0 & p.cprtop>0);

          ice_ind = find(lIce);
          count_ice(iday, ifor) = length(ice_ind);
          ctype_mean(iday, ifor) = nanmean(p.ctype(ice_ind));
          cngwat_mean(iday, ifor) = nanmean(p.cngwat(ice_ind));
          cpsize_mean(iday, ifor) = nanmean(p.cpsize(ice_ind));
          cprbot_mean(iday, ifor) = nanmean(p.cprbot(ice_ind));
          cprtop_mean(iday, ifor) = nanmean(p.cprtop(ice_ind));
          
          lWat = (p.ctype2==101 & p.cfrac2>0 & p.cngwat2>0 & p.cpsize2>0 ...
                   & p.cprbot2>0 & p.cprtop2>0);

          wat_ind = find(lWat);
          count_water(iday, ifor) = length(wat_ind); 
          ctype2_mean(iday, ifor) = nanmean(p.ctype2(wat_ind));
          cngwat2_mean(iday, ifor) = nanmean(p.cngwat2(wat_ind));
          cpsize2_mean(iday, ifor) = nanmean(p.cpsize2(wat_ind));
          cprbot2_mean(iday, ifor) = nanmean(p.cprbot2(wat_ind));
          cprtop2_mean(iday, ifor) = nanmean(p.cprtop2(wat_ind));

          % cloud fraction gets averaged over ALL obs
          cfrac_mean(iday, ifor) = nanmean(p.cfrac);
          cfrac2_mean(iday, ifor) = nanmean(p.cfrac2);
          cfrac12_mean(iday, ifor) = nanmean(p.cfrac12);
          
          % Radiance mean and std
          r  = p.robs1;
          cldy_calc = p.rcalc;
          clr_calc = p.sarta_rclearcalc;
          
          % spectral
          robs(iday,ifor,:) = nanmean(r,2);
          rcldy(iday,ifor,:) = nanmean(cldy_calc,2);
          rcldybias_std(iday,ifor,:) = nanstd(r-cldy_calc,0,2);
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

outfile = fullfile(statsdir, sprintf('rtp_airicrad_era_rad_kl_l%d_%s_random_fs_scanang_%s', ...
           llat, int2str(year), sDescriptor));
eval_str = ['save ' outfile [' robs rcl* *_mean count* ' ...
                    'trace']];
eval(eval_str);
