function pull_stats_cris_dcc(year, filter);

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

% Get proper frequencies for these data
[n1,n2,n3,userLW,userMW,userSW, ichan] = cris_lowres_chans();
f = cris_vchan(2, userLW, userMW, userSW);

basedir = fullfile('/asl/rtp/rtp_cris_ccast_lowres/dcc', ...
                   int2str(year));
dayfiles = dir(fullfile(basedir, 'era_rtp_d*_dcc.rtp'));

% calculate latitude bins
nbins=20; % gives 2N+1 element array of lat bin boundaries
latbins = equal_area_spherical_bands(nbins);
nlatbins = length(latbins);

iday = 1;
% for giday = 1:50:length(dayfiles)
for giday = 1:length(dayfiles)
   fprintf(1, '>>> year = %d  :: giday = %d\n', year, giday);
   a = dir(fullfile(basedir,dayfiles(giday).name));
   if a.bytes > 100000
      [h,ha,p,pa] = rtpread(fullfile(basedir,dayfiles(giday).name));

      switch filter
        case 1
          k = find(p.iudef(4,:) == 1); % descending node (night)
          sDescriptor='_desc';
% $$$         case 2
% $$$           k = find(p.iudef(4,:) == 1 & p.landfrac == 0); % descending node
% $$$                                                          % (night), ocean
% $$$           sDescriptor='_desc_ocean';
% $$$         case 3
% $$$           k = find(p.iudef(4,:) == 1 & p.landfrac == 1); % descending node
% $$$                                                         % (night), land
% $$$           sDescriptor='_desc_land';
        case 4
          k = find(p.iudef(4,:) == 0); % ascending node (day)
          sDescriptor='_asc';
% $$$         case 5
% $$$           k = find(p.iudef(4,:) == 0 & p.landfrac == 0); % ascending node
% $$$                                                          % (day), ocean
% $$$           sDescriptor='_asc_ocean';
% $$$         case 6
% $$$           k = find(p.iudef(4,:) == 0 & p.landfrac == 1); % ascending node
% $$$                                                         % (day), land
% $$$           sDescriptor='_asc_land';
      end

      pp = rtp_sub_prof(p, k);

      % initialize counts and look for bad channels (what should
      % the iasi bad channel test look like?)
      [nchans, nobs] = size(pp.robs1);
      nfovs = 9;
      count_all = int8(ones(nchans, nobs, nfovs));
      
      % grab any obs within +-60 deg lat
      
          % subset based on latitude bin
          inbin = find(pp.rlat >= -60 & pp.rlat <= 60);
          p = rtp_sub_prof(pp,inbin);
          
          for z = 1:9  % loop over FOVs to further sub-select
              ifov = find(p.ifov == z);
              p2 = rtp_sub_prof(p, ifov);
              
              bincount = count_all(:,inbin,z); 
              
              % Loop over obs in day
              % Radiance mean and std
              
              r  = p2.robs1;
              rc = p2.rcalc;
              
              % Convert r to rham
              r = box_to_ham(r);  % assumes r in freq order!!  Needed
                                  % for lowres
              robs(iday,:,z) = nanmean(r,2);
              rcal(iday,:,z) = nanmean(rc,2);
              rbias_std(iday,:,z) = nanstd(r-rc,0,2);
% $$$               bto = real(rad2bt(f,r));
% $$$               btc = real(rad2bt(f,rc));
% $$$               btobs(iday,:,z) = nanmean(bto,2);
% $$$               btcal(iday,:,z) = nanmean(btc,2);
% $$$               bias(iday,:,z)  = nanmean(bto-btc,2);
% $$$               bias_std(iday,:,z) = nanstd(bto-btc,0,2);
              lat_mean(iday,z) = nanmean(p2.rlat);
              lon_mean(iday,z) = nanmean(p2.rlon);
              solzen_mean(iday,z) = nanmean(p2.solzen);
              rtime_mean(iday,z)  = nanmean(p2.rtime);
              count(iday,z) = sum(bincount(1,:))';
              stemp_mean(iday,z) = nanmean(p2.stemp);
              iudef4_mean(iday,z) = nanmean(p2.iudef(4,:));
              ptemp_mean(iday,:,z) = nanmean(p.ptemp,2);
              gas1_mean(iday,:,z) = nanmean(p.gas_1,2);
% $$$               gas2_mean(iday,:,z) = nanmean(p.gas_2,2);                
              gas3_mean(iday,:,z) = nanmean(p.gas_3,2);
% $$$               gas4_mean(iday,:,z) = nanmean(p.gas_4,2);
% $$$               gas5_mean(iday,:,z) = nanmean(p.gas_5,2);
% $$$               gas6_mean(iday,:,z) = nanmean(p.gas_6,2);
% $$$               gas9_mean(iday,:,z) = nanmean(p.gas_9,2);
% $$$               gas12_mean(iday,:,z) = nanmean(p.gas_12,2);                    
% $$$               spres_mean(iday,z) = nanmean(p.spres);
              nlevs_mean(iday,z) = nanmean(p.nlevs);
              satzen_mean(iday,z) = nanmean(p.satzen);
              plevs_mean(iday,:,z) = nanmean(p.plevs,2);
% $$$               scanang_mean(iday,`z) = nanmean(p.scanang);
          end  % ifov (z)
          
          iday = iday + 1
   end % if a.bytes > 1000000
end  % giday
eval_str = ['save /asl/data/stats/cris/rtp_cris_lowres_rad_'  int2str(year) ...
            '_dcc' sDescriptor ' robs  *_mean count '];
% $$$             '_clear' sDescriptor '  btobs btcal bias bias_std *_mean count '];
eval(eval_str);
