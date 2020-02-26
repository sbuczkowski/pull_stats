function pull_stats_cris_hires_2_secantbin(year, filter);

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
klayers_exec = '/asl/packages/klayersV205/BinV201/klayers_airs_wetwater';

[sID, sTempPath] = genscratchpath();

% calculate latitude bins
% $$$ nbins=20; % gives 2N+1 element array of lat bin boundaries
% $$$ latbins = equal_area_spherical_bands(nbins);
% $$$ nlatbins = length(latbins);
nbins = 1;
latbins = [-50 50];
description = ['cris ccast hires binned according to sec(satzen). ' ...
               'Lat range restricted (-50, 50)'];

% Get proper frequencies for these data
[n1,n2,n3,userLW,userMW,userSW, ichan] = cris_hires_chans();
f = cris_vchan(2, userLW, userMW, userSW);

% $$$ basedir = '/home/sbuczko1/WorkingFiles/rtp_cris_ccast_hires/clear/2016';
% $$$ dayfiles = dir(fullfile(basedir, 'cris_hr_merra_*_clear.rtp'));
basedir = '/asl/rtp/rtp_cris_ccast_hires/clear/2016';
dayfiles = dir(fullfile(basedir, 'cris_hr_era_csarta_*_clear.rtp'));

iday = 1;
% for giday = 1:50:length(dayfiles)
for giday = 1:length(dayfiles)
   fprintf(1, '>>> year = %d  :: giday = %d\n', year, giday);
   a = dir(fullfile(basedir,dayfiles(giday).name));
   if a.bytes > 100000
      [h,ha,p,pa] = rtpread(fullfile(basedir,dayfiles(giday).name));

      latfiltered = find(p.rlat >= -50 & p.rlat <= 50);
      
      switch filter
        case 1
          k = find(p.iudef(4,:) == 1); % descending node (night)
          sDescriptor='_night';
        case 2
          k = find(p.iudef(4,:) == 1 & p.landfrac == 0); % descending
                                                     % node (night) ocean
          sDescriptor='_night_ocean';
        case 3
          k = find(p.iudef(4,:) == 1 & p.landfrac == 1); % descending node
                                                        % (night), land
          sDescriptor='_night_land';
        case 4
          k = find(p.iudef(4,:) == 0); % ascending node (day)
          sDescriptor='_day';
        case 5
          k = find(p.iudef(4,:) == 0 & p.landfrac == 0); % ascending node
                                                         % (day), ocean
          sDescriptor='_day_ocean';
        case 6
          k = find(p.iudef(4,:) == 0 & p.landfrac == 1); % ascending node
                                                        % (day), land
          sDescriptor='_day_land';
      end

      pp = rtp_sub_prof(p, intersect(k, latfiltered));
% $$$       mmwater = mmwater(k);

      % initialize counts and look for bad channels (what should
      % the iasi bad channel test look like?)
      [nchans, nobs] = size(pp.robs1);
      nfovs = 9;
      count_all = int8(ones(nchans, nobs, nfovs));

      % calculate 10 equal secant bins for secant(satzen)
      [N,edges, bind] = histcounts(secd(pp.satzen),10);
      
      % loop over secant bins
      for isec = 1:10
          % subset based on sec(satzen) bin
          inbin = find(bind == isec);
          p = rtp_sub_prof(pp,inbin);
          
          for z = 1:9  % loop over FOVs to further sub-select
              ifov = find(p.ifov == z);
              p2 = rtp_sub_prof(p, ifov);
              
              bincount = count_all(:,inbin,z); 
          
% Radiance mean and std
         r  = p.robs1;
         rc = p.rcalc;

         % Convert r to rham
         r = box_to_ham(r);  % assumes r in freq order!!

         % B(T) bias mean and std
% $$$          bto = real(rad2bt(f,r));
% $$$          btc = real(rad2bt(f,rc));
% $$$          btobs(iday,isec,:) = nanmean(bto,2);
% $$$          btcal(iday,isec,:) = nanmean(btc,2);
% $$$          bias(iday,isec,:)  = nanmean(bto-btc,2);
% $$$          bias_std(iday,isec,:) = nanstd(bto-btc,0,2);
              robs(iday,isec,:,z) = nanmean(r,2);
              rcal(iday,isec,:,z) = nanmean(rc,2);
              rbias_std(iday, isec,:) = nanstd(r-rc,0,2);
         lat_mean(iday,isec,z) = nanmean(p.rlat);
         lon_mean(iday,isec,z) = nanmean(p.rlon);
         solzen_mean(iday,isec,z) = nanmean(p.solzen);
         rtime_mean(iday,isec,z)  = nanmean(p.rtime);
         count(iday,isec,:,z) = sum(bincount,2)';
         stemp_mean(iday,isec,z) = nanmean(p.stemp);
         ptemp_mean(iday,isec,:,z) = nanmean(p.ptemp,2);
         gas1_mean(iday,isec,:,z) = nanmean(p.gas_1,2);
         gas3_mean(iday,isec,:,z) = nanmean(p.gas_3,2);
         spres_mean(iday,isec,z) = nanmean(p.spres);
         nlevs_mean(iday,isec,z) = nanmean(p.nlevs);
         iudef4_mean(iday,isec,z) = nanmean(p.iudef(4,:));
         satzen_mean(iday,isec,z) = nanmean(p.satzen);
         plevs_mean(iday,isec,:,z) = nanmean(p.plevs,2);
         end % ifov(z)
      end  % end loop over latitudes
      iday = iday + 1
   end % if a.bytes > 1000000
end  % giday
eval_str = ['save /asl/data/stats/cris/cris_hr_era_rad_csarta_secbin_'...
            int2str(year) '_clear' sDescriptor ' robs rcal rbias_std *_mean edges count description'];

eval(eval_str);
