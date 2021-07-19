function pull_stats_cris_hires_2(year, filter);

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
addpath /home/sergio/MATLABCODE/PLOTTER  %
                                         % equal_area_spherical_bands
addpath /asl/matlib/rtptools  % mmwater_rtp.m
addpath /home/sbuczko1/git/rtp_prod2/cris  % cris_hires_chans
klayers_exec = '/asl/packages/klayersV205/BinV201/klayers_airs_wetwater';

[sID, sTempPath] = genscratchpath();

% calculate latitude bins
nbins=20; % gives 2N+1 element array of lat bin boundaries
latbins = equal_area_spherical_bands(nbins);
nlatbins = length(latbins);

% Get proper frequencies for these data
[n1,n2,n3,userLW,userMW,userSW, ichan] = cris_hires_chans();
f = cris_vchan(2, userLW, userMW, userSW);

basedir = '/asl/rtp/rtp_cris_ccast_hires/clear_daily/2015';
sSubset = 'clear';
dayfiles = dir(fullfile(basedir, ['rtp*_' sSubset '.rtp']));

iday = 1;
% for giday = 1:50:length(dayfiles)
for giday = 1:length(dayfiles)
   fprintf(1, '>>> year = %d  :: giday = %d\n', year, giday);
   a = dir(fullfile(basedir,dayfiles(giday).name));
   if a.bytes > 100000
      [h,ha,p,pa] = rtpread(fullfile(basedir,dayfiles(giday).name));

% $$$       % run klayers on the rtp data (Sergio is asking for this to
% $$$       % convert levels to layers for his processing?)
% $$$ 
% $$$       % first remove rcalc field and save it for later restore
% $$$       rcalc = p.rcalc;
% $$$       p = rmfield(p, 'rcalc');
% $$$       
% $$$       fprintf(1, '>>> running klayers... ');
% $$$       fn_rtp1 = fullfile(sTempPath, ['airs_' sID '_1.rtp']);
% $$$       rtpwrite(fn_rtp1, h,ha,p,pa)
% $$$       fn_rtp2 = fullfile(sTempPath, ['airs_' sID '_2.rtp']);
% $$$       klayers_run = [klayers_exec ' fin=' fn_rtp1 ' fout=' fn_rtp2 ...
% $$$                      ' > ' sTempPath '/kout.txt'];
% $$$       unix(klayers_run);
% $$$       [h,ha,p,pa] = rtpread(fn_rtp2);
% $$$       % restore rcalc
% $$$       p.rcalc = rcalc;
% $$$       clear rcalc;
% $$$       
% $$$       fprintf(1, 'Done\n');
% $$$ 
% $$$       % get column water
% $$$       mmwater = mmwater_rtp(h, p);

      switch filter
        case 1
          k = find(p.solzen > 90); % descending node (night)
          sDescriptor='_desc';
        case 2
          k = find(p.solzen > 90 & p.landfrac == 0); % descending
                                                     % node (night) ocean
          sDescriptor='_desc_ocean';
        case 3
          k = find(p.solzen > 90 & p.landfrac == 1); % descending node
                                                        % (night), land
          sDescriptor='_desc_land';
        case 4
          k = find(p.solzen < 90); % ascending node (day)
          sDescriptor='_asc';
        case 5
          k = find(p.solzen < 90 & p.landfrac == 0); % ascending node
                                                         % (day), ocean
          sDescriptor='_asc_ocean';
        case 6
          k = find(p.solzen < 90 & p.landfrac == 1); % ascending node
                                                        % (day), land
          sDescriptor='_asc_land';
      end

      pp = rtp_sub_prof(p, k);
% $$$       mmwater = mmwater(k);

      % initialize counts and look for bad channels (what should
      % the iasi bad channel test look like?)
      [nchans, nobs] = size(pp.robs1);
      nfovs = 9;
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
          
% Radiance mean and std
         r  = p.robs1;
         rc = p.rcalc;

         % Convert r to rham
% $$$          r = box_to_ham(r);  % assumes r in freq order!!

         % B(T) bias mean and std
% $$$          bto = real(rad2bt(f,r));
% $$$          btc = real(rad2bt(f,rc));
% $$$          btobs(iday,ilat,:) = nanmean(bto,2);
% $$$          btcal(iday,ilat,:) = nanmean(btc,2);
% $$$          bias(iday,ilat,:)  = nanmean(bto-btc,2);
% $$$          bias_std(iday,ilat,:) = nanstd(bto-btc,0,2);
              robs(iday,ilat,:,z) = nanmean(r,2);
              rcal(iday,ilat,:,z) = nanmean(rc,2);
              rbias_std(iday, ilat,:) = nanstd(r-rc,0,2);
         lat_mean(iday,ilat,z) = nanmean(p.rlat);
         lon_mean(iday,ilat,z) = nanmean(p.rlon);
         solzen_mean(iday,ilat,z) = nanmean(p.solzen);
         rtime_mean(iday,ilat,z)  = nanmean(p.rtime);
         count(iday,ilat,:,z) = sum(bincount,2)';
         stemp_mean(iday,ilat,z) = nanmean(p.stemp);
         ptemp_mean(iday,ilat,:,z) = nanmean(p.ptemp,2);
         gas1_mean(iday,ilat,:,z) = nanmean(p.gas_1,2);
         gas3_mean(iday,ilat,:,z) = nanmean(p.gas_3,2);
         spres_mean(iday,ilat,z) = nanmean(p.spres);
         nlevs_mean(iday,ilat,z) = nanmean(p.nlevs);
         iudef4_mean(iday,ilat,z) = nanmean(p.iudef(4,:));
         satzen_mean(iday,ilat,z) = nanmean(p.satzen);
         plevs_mean(iday,ilat,:,z) = nanmean(p.plevs,2);
         end % ifov(z)
      end  % end loop over latitudes
      iday = iday + 1
   end % if a.bytes > 1000000
end  % giday
eval_str = ['save /asl/data/stats/cris/rtp_cris_hires_era_sinc_rad_'  int2str(year) ...
            '_clear' sDescriptor ' robs rcal rbias_std *_mean count '];

eval(eval_str);
