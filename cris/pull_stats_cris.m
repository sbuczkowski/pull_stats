function pull_stats_cris(year, filter, cfg);

[sID, sTempPath] = genscratchpath();
if isfield(cfg, 'stemppath')
    sTempPath = cfg.sTempPath;
    sID = '99999';
end
fprintf(1, '>> TempPath: %s\n', sTempPath)

% record run start datetime in output stats file for tracking
trace.RunDate = datetime('now','TimeZone','local','Format', ...
                         'd-MMM-y HH:mm:ss Z');
trace.Reason = 'Normal pull_stats runs';
if isfield(cfg, 'reason')
    trace.Reason = cfg.reason;
end

bRunKlayers = true;
klayers_exec = ['/asl/packages/klayersV205/BinV201/' ...
                'klayers_airs_wetwater'];
if isfield(cfg, 'klayers') & cfg.klayers == false
    bRunKlayers = false;
end
trace.klayers = bRunKlayers;
trace.klayers_exec = klayers_exec;

trace.droplayers = false;
if isfield(cfg, 'droplayers') & cfg.droplayers == true
    trace.droplayers = true;
end

model = 'era';
if isfield(cfg, 'model')
    model = cfg.model;
end

instname = 'npp';   % default is NPP
if isfield(cfg, 'instname') 
    instname = cfg.instname;
end

bHiRes = false;
fscale = 'lowres';
if isfield(cfg, 'fscale') & strcmp(cfg.fscale, 'hires')
    bHiRes = true;
    fscale = 'hires';
end

rtpsrcdir = sprintf('/asl/rtp/cris/%s_ccast_%s', instname, fscale);
if isfield(cfg, 'rtpsrcdir')
    rtpsrcdir = cfg.rtpsrcdir;
end

descriptor = 'random';
if isfield(cfg, 'descriptor')
    descriptor = cfg.descriptor;
end
needscloudfields = (~strcmp(descriptor, 'clear') & ~strcmp(descriptor, ...
                                                  'dcc'));

statsdir = sprintf('/asl/stats/cris/%s/%s', instname, descriptor);
if isfield(cfg, 'statsdir')
    statsdir = cfg.statsdir;
end
        if exist(statsdir) == 0
            mkdir(statsdir);
        end

rta = 'csarta';
if isfield(cfg, 'rta')
    rta = cfg.rta;
end

basedir = fullfile(rtpsrcdir, sprintf('%s', descriptor), ...
                   int2str(year));
fprintf(1, '>> Basedir: %s\n', basedir);
% $$$ namefilter = sprintf('%s_rtp_d*_%s.rtp', model, descriptor);
namefilter = cfg.namefilter;
fprintf(1, '>> namefilter: %s\n', namefilter)
dayfiles = dir(fullfile(basedir, namefilter));
                                         
                                         
ndays = length(dayfiles);
fprintf(1,'>>> numfiles = %d\n', ndays);

% Get proper frequencies for these data
if true == bHiRes
    [n1,n2,n3,userLW,userMW,userSW, ichan] = cris_hires_chans();
else
    [n1,n2,n3,userLW,userMW,userSW, ichan] = cris_lowres_chans();
end
f = cris_vchan(2, userLW, userMW, userSW);
nchans = length(f);

% calculate latitude bins
nbins=20; % gives 2N+1 element array of lat bin boundaries
if isfield(cfg, 'nbins')
    nbins = cfg.nbins;
end
latbin_edges = equal_area_spherical_bands(nbins);
nlatbins = length(latbin_edges)-1;

nlevs = 101;  % klayers output
nfovs = 9;    % FOVs/FOR

% allocate final accumulator arrays
robs = nan(ndays, nlatbins, nfovs, nchans);
rclr = nan(ndays, nlatbins, nfovs, nchans);
rclrbias_std = nan(ndays, nlatbins, nfovs, nchans);

lat_mean = nan(ndays, nlatbins, nfovs);
lon_mean = nan(ndays, nlatbins, nfovs);
solzen_mean = nan(ndays, nlatbins, nfovs);
rtime_mean = nan(ndays, nlatbins, nfovs); 
count = nan(ndays, nlatbins, nfovs, nchans);
% $$$ tcc_mean = nan(ndays, nlatbins, nfovs);
stemp_mean = nan(ndays, nlatbins, nfovs);
ptemp_mean = nan(ndays, nlatbins, nfovs, nlevs);
gas1_mean = nan(ndays, nlatbins, nfovs, nlevs);
gas3_mean = nan(ndays, nlatbins, nfovs, nlevs);
spres_mean = nan(ndays, nlatbins, nfovs);
nlevs_mean = nan(ndays, nlatbins, nfovs);
iudef4_mean = nan(ndays, nlatbins, nfovs);
mmwater_mean = nan(ndays, nlatbins, nfovs);
satzen_mean = nan(ndays, nlatbins, nfovs);
plevs_mean = nan(ndays, nlatbins, nfovs, nlevs);

if needscloudfields
    rcld = nan(ndays, nlatbins, nfovs, nchans);
    rcldbias_std = nan(ndays, nlatbins, nfovs, nchans);

    count_ice = nan(ndays, nlatbins, nfovs);
    ctype_mean = nan(ndays, nlatbins, nfovs);
    cngwat_mean = nan(ndays, nlatbins, nfovs);
    cpsize_mean = nan(ndays, nlatbins, nfovs);
    cprbot_mean = nan(ndays, nlatbins, nfovs);
    cprtop_mean = nan(ndays, nlatbins, nfovs);

    count_water = nan(ndays, nlatbins, nfovs);
    ctype2_mean = nan(ndays, nlatbins, nfovs);
    cngwat2_mean = nan(ndays, nlatbins, nfovs);
    cpsize2_mean = nan(ndays, nlatbins, nfovs);
    cprbot2_mean = nan(ndays, nlatbins, nfovs);
    cprtop2_mean = nan(ndays, nlatbins, nfovs);

    % cloud fraction gets averaged over ALL obs
    cfrac_mean = nan(ndays, nlatbins, nfovs);
    cfrac2_mean = nan(ndays, nlatbins, nfovs);
    cfrac12_mean = nan(ndays, nlatbins, nfovs);

end


iday = 1;
for giday = 1:ndays
   fprintf(1, '>>> year = %d  :: giday = %d\n', year, giday);
   a = dir(fullfile(basedir,dayfiles(giday).name));
   if a.bytes > 100000
       day_str = sprintf('giday = %d, iday = %d',giday, iday);
       
       [h,ha,p_rtp,pa] = rtpread(fullfile(basedir,dayfiles(giday).name));
      

       nightobs = find(p_rtp.iudef(4,:) == 1);
       dayobs = find(p_rtp.iudef(4,:) == 0);
       landobs = find(p_rtp.landfrac == 1);
       oceanobs = find(p_rtp.landfrac == 0);

      switch filter
        case 1 % descending node (night), global
          k = nightobs;
          sDescriptor='desc';
        case 2 % descending node (night), ocean
          k = intersect(nightobs, oceanobs);
          sDescriptor='desc_ocean';
        case 3 % descending node (night), land
          k = intersect(nightobs, landobs);
          sDescriptor='desc_land';
        case 4 % ascending node (day)
          k = dayobs;
          sDescriptor='asc';
        case 5 % ascending node (day), ocean
          k = intersect(dayobs, oceanobs);
          sDescriptor='asc_ocean';
        case 6 % ascending node (day), land
          k = intersect(dayobs, landobs);
          sDescriptor='asc_land';
      end
      if length(k) == 0
          % no obs survive filter so on to next loop iteration
          continue
      end

      p_filt = rtp_sub_prof(p_rtp, k);
      clear p_rtp
      % run klayers on the rtp data to convert levels -> layers
      % save calcs as the re-run of klayers wipes them out
      tmp_rclr = p_filt.rclr;
% $$$       tmp_tcc = p_filt.tcc;
      if needscloudfields
          tmp_rcld = p_filt.rcld;
      end

      % rtps are in levels but h.ptype is set to 1 for layers. need
      % to reset this to 0 before running klayers (because we have
      % started trying to preserve attributes from the
      % klayers/sarta chain. ultimately, this needs to be fixed in
      % rtp generation)
% $$$       load('/home/sbuczko1/Work/scratch/cris_header')
% $$$       h = hr;
% $$$       clear hr
      h.ptype=0;
      h.pfields=5;
      h.ngas=2;
      h.glist=[1;3];
      h.gunit=[21; 21];
      fprintf(1, '>>> running klayers... ');
      fn_rtp1 = fullfile(sTempPath, ['cris_' sID '_1.rtp']);
      rtpwrite(fn_rtp1, h,ha,p_filt,pa);
      clear p_filt;
      fn_rtp2 = fullfile(sTempPath, ['cris_' sID '_2.rtp']);
      klayers_run = [klayers_exec ' fin=' fn_rtp1 ' fout=' fn_rtp2 ...
                     ' > /home/sbuczko1/Work/scratch/kout.txt'];
      unix(klayers_run);
      fprintf(1, 'Done\n');

      % Read klayers output into local rtp variables
      [h,ha,p_filt,pa] = rtpread(fn_rtp2);
      % restore rclr
      p_filt.rclr = tmp_rclr;
% $$$       p_filt.tcc = tmp_tcc;
      clear tmp_rclr 
      if needscloudfields
          p_filt.rcld = tmp_rcld;
          clear tmp_rcld
      end
      
      % get column water
      mmwater = mmwater_rtp(h, p_filt);

      % Check for obs with layer profiles that go lower than
      % topography. Need to check nlevs and NaN out any layers
      % at or below this level

      % ** Any layers-sensitive variables added in averaging code below must
      % ** be checked here first.
      for i=1:length(p_filt.nlevs)
          badlayers = p_filt.nlevs(i) : 101;
          p_filt.plevs(badlayers, i) = NaN;
          p_filt.gas_1(badlayers, i) = NaN;
          p_filt.gas_3(badlayers, i) = NaN;
          p_filt.ptemp(badlayers, i) = NaN;
      end

      % initialize counts and look for bad channels (what should
      % the iasi bad channel test look like?)
      [nchans, nobs] = size(p_filt.robs1);
      
      % loop over latitude bins
      for ilat = 1:nlatbins-1
          % subset based on latitude bin
          inbin = find(p_filt.rlat > latbin_edges(ilat) & p_filt.rlat <= ...
                       latbin_edges(ilat+1));
% $$$           p_inbin = rtp_sub_prof(p_filt,inbin);

          daylat_str = sprintf('%s, ilat = %d', day_str, ilat);
          
          for z = 1:9  % loop over FOVs to further sub-select
              infov = find(p_filt.ifov == z);
              inbin_fov = intersect(inbin, infov);
              p_infov = rtp_sub_prof(p_filt, inbin_fov);

              daylatfov_str = sprintf('%s, z = %d', daylat_str, z);
              
              count_infov = ones(length(inbin_fov), nchans);
              % QA/QC checks for bad chans, etc go here and set
              % elements of count_infov to zero

              binwater = mmwater(inbin_fov);        

              if needscloudfields
                  % Remove 'clouds' that have only partial defining
                  % characteristics (effects 0.2% of obs in test rtp data)
                  % cfrac=cngwat=0 but cpsize/cprtop/cprbot ~= 0
                  a = find(p_infov.cfrac==0 & p_infov.cngwat==0 & p_infov.cprtop>0);
                  p_infov.ctype(a)=-1;
                  a = find(p_infov.cfrac2==0 & p_infov.cngwat2==0 & p_infov.cprtop2>0);
                  p_infov.ctype2(a)=-1;

                  % replace 'bad data' placeholders in cprtop/cprbot (-9999
                  % -> NaN, esesntially except for 'clouds' from above)
                  a = find(p_infov.ctype==-1);
                  p_infov.cprtop(a)=NaN; p_infov.cprbot(a)=NaN;
                  a = find(p_infov.ctype2==-1);
                  p_infov.cprtop2(a)=NaN; p_infov.cprbot2(a)=NaN;

                  % massage situations where we have two of same cloud type
                  % in a single obs i.e. p_infov.ctype=p_infov.ctype2 (average to a
                  % single cloud)
                  dblwat = find(p_infov.ctype==101 & p_infov.ctype2==101);
                  p_infov.cfrac2(dblwat) = p_infov.cfrac(dblwat) + p_infov.cfrac2(dblwat);
                  p_infov.cpsize2(dblwat) = (p_infov.cpsize(dblwat) + ...
                                        p_infov.cpsize2(dblwat))/2;
                  p_infov.cprtop2(dblwat) = (p_infov.cprtop(dblwat) + ...
                                        p_infov.cprtop2(dblwat))/2;
                  p_infov.cprbot2(dblwat) = (p_infov.cprbot(dblwat) + ...
                                        p_infov.cprbot2(dblwat))/2;
                  p_infov.cngwat2(dblwat) = (p_infov.cfrac(dblwat).*p_infov.cngwat(dblwat) + ...
                                        p_infov.cfrac2(dblwat).*p_infov.cngwat2(dblwat))./ ...
                      p_infov.cfrac2(dblwat);
                  p_infov.cfrac(dblwat)=0; p_infov.cpsize(dblwat)=0;
                  p_infov.cngwat(dblwat)= 0; p_infov.cprtop(dblwat)=NaN;
                  p_infov.cprbot(dblwat)=NaN; p_infov.ctype(dblwat)=-1;
                  
                  dblice = find(p_infov.ctype==201 & p_infov.ctype2==201);
                  p_infov.cfrac(dblice) = p_infov.cfrac(dblice) + p_infov.cfrac2(dblice);
                  p_infov.cpsize(dblice) = (p_infov.cpsize(dblice) + ...
                                       p_infov.cpsize2(dblice))/2;
                  p_infov.cprtop(dblice) = (p_infov.cprtop(dblice) + ...
                                       p_infov.cprtop2(dblice))/2;
                  p_infov.cprbot(dblice) = (p_infov.cprbot(dblice) + ...
                                       p_infov.cprbot2(dblice))/2;
                  p_infov.cngwat(dblice) = (p_infov.cfrac(dblice).*p_infov.cngwat(dblice) + ...
                                       p_infov.cfrac2(dblice).*p_infov.cngwat2(dblice))./ ...
                      p_infov.cfrac2(dblice);
                  p_infov.cfrac2(dblice)=0; p_infov.cpsize2(dblice)=0;
                  p_infov.cngwat2(dblice)= 0; p_infov.cprtop2(dblice)=NaN;
                  p_infov.cprbot2(dblice)=NaN; p_infov.ctype2(dblice)=-1;
                  
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
                  lIce = (p_infov.ctype==201 & p_infov.cfrac>0 & p_infov.cngwat>0 & p_infov.cpsize>0 ...
                          & p_infov.cprbot>0 & p_infov.cprtop>0);

                  ice_ind = find(lIce);
                  count_ice(iday, ilat,z) = length(ice_ind);
                  ctype_mean(iday, ilat,z) = nanmean(p_infov.ctype(ice_ind));
                  cngwat_mean(iday, ilat,z) = nanmean(p_infov.cngwat(ice_ind));
                  cpsize_mean(iday, ilat,z) = nanmean(p_infov.cpsize(ice_ind));
                  cprbot_mean(iday, ilat,z) = nanmean(p_infov.cprbot(ice_ind));
                  cprtop_mean(iday, ilat,z) = nanmean(p_infov.cprtop(ice_ind));
                  
                  lWat = (p_infov.ctype2==101 & p_infov.cfrac2>0 & p_infov.cngwat2>0 & p_infov.cpsize2>0 ...
                          & p_infov.cprbot2>0 & p_infov.cprtop2>0);

                  wat_ind = find(lWat);
                  count_water(iday, ilat,z) = length(wat_ind); 
                  ctype2_mean(iday, ilat,z) = nanmean(p_infov.ctype2(wat_ind));
                  cngwat2_mean(iday, ilat,z) = nanmean(p_infov.cngwat2(wat_ind));
                  cpsize2_mean(iday, ilat,z) = nanmean(p_infov.cpsize2(wat_ind));
                  cprbot2_mean(iday, ilat,z) = nanmean(p_infov.cprbot2(wat_ind));
                  cprtop2_mean(iday, ilat,z) = nanmean(p_infov.cprtop2(wat_ind));

                  % cloud fraction gets averaged over ALL obs
                  cfrac_mean(iday, ilat,z) = nanmean(p_infov.cfrac);
                  cfrac2_mean(iday, ilat,z) = nanmean(p_infov.cfrac2);
                  cfrac12_mean(iday, ilat,z) = nanmean(p_infov.cfrac12);
              end
              
              % Convert robs to Hamming           
              r = box_to_ham(p_infov.robs1);  % assumes r in freq order!!  Needed
                                  % for lowres

              robs(iday,ilat,z,:) = nanmean(r,2);
              rclr(iday,ilat,z,:) = nanmean(p_infov.rclr,2);
              rclrbias_std(iday, ilat,z,:) = nanstd(r-p_infov.rclr,0,2);

              if needscloudfields
                  rcld(iday,ilat,z,:) = nanmean(p_infov.rcld,2);
                  rcldbias_std(iday, ilat,z,:) = nanstd(r-p_infov.rcld,0, ...
                                                      2);
              end
              
              lat_mean(iday,ilat,z) = nanmean(p_infov.rlat);
              lon_mean(iday,ilat,z) = nanmean(p_infov.rlon);
              solzen_mean(iday,ilat,z) = nanmean(p_infov.solzen);
              rtime_mean(iday,ilat,z)  = nanmean(p_infov.rtime);
              count(iday,ilat,z,:) = sum(count_infov);
              stemp_mean(iday,ilat,z) = nanmean(p_infov.stemp);
              iudef4_mean(iday,ilat,z) = nanmean(p_infov.iudef(4,:));
              ptemp_mean(iday,ilat,z,:) = nanmean(p_infov.ptemp,2)';
              gas1_mean(iday,ilat,z,:) = nanmean(p_infov.gas_1,2)';
              gas3_mean(iday,ilat,z,:) = nanmean(p_infov.gas_3,2)';
              spres_mean(iday,ilat,z) = nanmean(p_infov.spres);
              nlevs_mean(iday,ilat,z) = nanmean(p_infov.nlevs);
              satzen_mean(iday,ilat,z) = nanmean(p_infov.satzen);
              plevs_mean(iday,ilat,z,:) = nanmean(p_infov.plevs,2)';
              mmwater_mean(iday,ilat) = nanmean(binwater);
% $$$               scanang_mean(iday,ilat,z) = nanmean(p.scanang);
          end  % ifov (z)
      end  % end loop over ilat
          
          iday = iday + 1
   end % if a.bytes > 1000000
end  % giday
outfile = fullfile(statsdir, sprintf('cris_%s_%s_rad_scanangle_%4d_%s_%s', ...
           instname, fscale, year, descriptor, ...
                                     sDescriptor));

eval_str = sprintf('save %s trace robs rclr *_std *_mean count -v7.3', ...
                   outfile);
if needscloudfields
    eval_str = sprintf('save %s trace robs rclr rcld *_std *_mean count -v7.3', ...
                       outfile);
end
fprintf(1, '>> Executing save command: \n\t%s\n', eval_str)
eval(eval_str);

fprintf(1, 'Done\n')
