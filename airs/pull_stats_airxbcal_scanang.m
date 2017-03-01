function pull_stats_airxbcal_scanang(year, filter);

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
addpath /home/sergio/MATLABCODE/PLOTTER  % equal_area_spherical_bands

trace.RunDate = datetime('now','TimeZone','local','Format', ...
                         'd-MMM-y HH:mm:ss Z');

cstr =[ 'bits1-4=NEdT[0.08 0.12 0.15 0.20 0.25 0.30 0.35 0.4 0.5 0.6 0.7' ...
  ' 0.8 1.0 2.0 4.0 nan]; bit5=Aside[0=off,1=on]; bit6=Bside[0=off,1=on];' ...
  ' bits7-8=calflag&calchansummary[0=OK, 1=DCR, 2=moon, 3=other]' ];

basedir = ['/asl/data/rtp_airxbcal_v5/' int2str(year) '/clear'];
dayfiles = dir(fullfile(basedir, 'ecmwf_airxbcal*.rtp'));
nfiles = length(dayfiles);
fprintf(1,'>>> numfiles = %d\n', nfiles);
if nfiles == 0
    fprintf(2, '>>> Error: No input files found. Exiting\n');
    return;
end

% calculate bins
% use p.xtrack as scanang proxy and bin every three FOVs
binedges = [1 [1:30]*3];   
nbins = length(binedges);

iday = 1;
% for giday = 1:50:length(dayfiles)
for giday = 1:length(dayfiles)
   fprintf(1, '>>> year = %d  :: giday = %d\n', year, giday);
   a = dir(fullfile(basedir,dayfiles(giday).name));
   a.bytes;
   if a.bytes > 100000
      [h,ha,p,pa] = rtpread(fullfile(basedir,dayfiles(giday).name));
      f = h.vchan;  % AIRS proper frequencies
      
      switch filter
        case 1         % descending node (night), tropics
          k = find(p.iudef(4,:) == 68 & abs(p.rlat) <= 30);
          sDescriptor='_desc_trop';
        case 2          % descending node (night), mid-lat
          k = find(p.iudef(4,:) == 68 & abs(p.rlat) > 30 & abs(p.rlat) <= 50);
          sDescriptor='_desc_midlat';
      end

      pp = rtp_sub_prof(p, k);

      % Look for bad channels and initialize counts
      [nedt,ab,ical] = calnum_to_data(pp.calflag,cstr);
      n = length(pp.rlat);
      count_all = ones(2378,n);
      for i=1:2378
         % Find bad channels
         k = find( pp.robs1(i,:) == -9999 | ical(i,:) ~= 0 | nedt(i,:) > 1);
%          % These are the good channels
%          kg = setdiff(1:n,k);
% NaN's for bad channels
         pp.robs1(i,k) = NaN;
         pp.rcalc(i,k) = NaN;
         count_all(i,k) = 0;
      end

      % Loop over latitude bins
      for ibin = 1:nbins-1
           % subset based on latitude bin
           inbin = find(pp.xtrack > binedges(ibin) & pp.xtrack <= ...
                        binedges(ibin+1));
           p = rtp_sub_prof(pp,inbin);
           bincount = count_all(:,inbin);
           
           % Radiance mean and std
           r  = p.robs1;
           rc = p.rcalc;

           % B(T) bias mean and std
           bto = real(rad2bt(f,r));
           btc = real(rad2bt(f,rc));
           btobs(iday,ibin,:) = nanmean(bto,2);
           btcal(iday,ibin,:) = nanmean(btc,2);
           bias_std(iday,ibin,:) = nanstd(bto-btc,0,2);
           lat_mean(iday,ibin) = nanmean(p.rlat);
           lon_mean(iday,ibin) = nanmean(p.rlon);
           solzen_mean(iday,ibin) = nanmean(p.solzen);
           rtime_mean(iday,ibin)  = nanmean(p.rtime);
           count(iday,ibin,:) = sum(bincount,2)';
           stemp_mean(iday,ibin) = nanmean(p.stemp);
           ptemp_mean(iday,ibin,:) = nanmean(p.ptemp,2);
           gas1_mean(iday,ibin,:) = nanmean(p.gas_1,2);
           gas3_mean(iday,ibin,:) = nanmean(p.gas_3,2);
           spres_mean(iday,ibin) = nanmean(p.spres);
           nlevs_mean(iday,ibin) = nanmean(p.nlevs);
           iudef4_mean(iday,ibin) = nanmean(p.iudef(4,:));
           satzen_mean(iday,ibin) = nanmean(p.satzen);
           plevs_mean(iday,ibin,:) = nanmean(p.plevs,2);
       end  % end loop over latitudes

      iday = iday + 1
   end % if a.bytes > 1000000
end  % giday
eval_str = ['save ~/testoutput/2015/airs/clear/rtp_airxbcal_ecmwf_scanang_'  int2str(year) ...
            '_clear' sDescriptor ' btobs btcal bias_std *_mean count trace'];
eval(eval_str);
