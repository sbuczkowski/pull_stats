function pull_stats_airxbcal(year, filter, cfg);

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

cstr =[ 'bits1-4=NEdT[0.08 0.12 0.15 0.20 0.25 0.30 0.35 0.4 0.5 0.6 0.7' ...
  ' 0.8 1.0 2.0 4.0 nan]; bit5=Aside[0=off,1=on]; bit6=Bside[0=off,1=on];' ...
  ' bits7-8=calflag&calchansummary[0=OK, 1=DCR, 2=moon, 3=other]' ];

basedir = ['/asl/rtp/rtp_airxbcal_v5/' int2str(year) '/clear'];
dayfiles = dir(fullfile(basedir, 'era_airxbcal_day*_clear.rtp'));
fprintf(1,'>>> numfiles = %d\n', length(dayfiles));

% calculate latitude bins
nbins=20; % gives 2N+1 element array of lat bin boundaries
latbins = equal_area_spherical_bands(nbins);
nlatbins = length(latbins);

iday = 1;
%for giday = 1:10:length(dayfiles)
for giday = 1:length(dayfiles)
   fprintf(1, '>>> year = %d  :: giday = %d\n', year, giday);
   a = dir(fullfile(basedir,dayfiles(giday).name));
   a.bytes;
   if a.bytes > 100000
      [h,ha,p,pa] = rtpread(fullfile(basedir,dayfiles(giday).name));
      f = h.vchan;  % AIRS proper frequencies
      
      switch filter
        case 1
          k = find(p.iudef(4,:) == 68); % descending node (night)
          sDescriptor='_desc';
        case 2
          k = find(p.iudef(4,:) == 68 & p.landfrac == 0); % descending node
                                                         % (night), ocean
          sDescriptor='_desc_ocean';
        case 3
          k = find(p.iudef(4,:) == 68 & p.landfrac == 1); % descending node
                                                        % (night), land
          sDescriptor='_desc_land';
        case 4
          k = find(p.iudef(4,:) == 65); % ascending node (night)
          sDescriptor='_asc';
        case 5
          k = find(p.iudef(4,:) == 65 & p.landfrac == 0); % ascending node
                                                         % (night), ocean
          sDescriptor='_asc_ocean';
        case 6
          k = find(p.iudef(4,:) == 65 & p.landfrac == 1); % ascending node
                                                        % (night), land
          sDescriptor='_asc_land';
        case 7
          k = find(abs(p.rlat) < 30 & p.landfrac == 0); % tropical
                                                        % ocean
          sDescriptor='_tropocean';
      end

      pp = rtp_sub_prof(p, k);
      clear p;
      
      if bRunKlayers
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
          
          % get column water
          mmwater = mmwater_rtp(h, pp);
      end
      
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
         count_all(i,k) = 0;
      end

      % Loop over latitude bins
      for ilat = 1:nlatbins-1
           % subset based on latitude bin
           inbin = find(pp.rlat > latbins(ilat) & pp.rlat <= ...
                        latbins(ilat+1));
           p = rtp_sub_prof(pp,inbin);
           bincount = count_all(:,inbin);
           if bRunKlayers
               binwater = mmwater(inbin);
           end
           % Radiance mean and std
           r  = p.robs1;
           rc = p.rcalc;

           robs(iday,ilat,:) = nanmean(r,2);
           rcal(iday,ilat,:) = nanmean(rc,2);
           rbias_std(iday,ilat,:) = nanstd(r-rc,0,2);
           
           lat_mean(iday,ilat) = nanmean(p.rlat);
           lon_mean(iday,ilat) = nanmean(p.rlon);
           solzen_mean(iday,ilat) = nanmean(p.solzen);
           rtime_mean(iday,ilat)  = nanmean(p.rtime);
           count(iday,ilat,:) = sum(bincount,2)';
           stemp_mean(iday,ilat) = nanmean(p.stemp);
           ptemp_mean(iday,ilat,:) = nanmean(p.ptemp,2);
           gas1_mean(iday,ilat,:) = nanmean(p.gas_1,2);
           gas3_mean(iday,ilat,:) = nanmean(p.gas_3,2);
           spres_mean(iday,ilat) = nanmean(p.spres);
           nlevs_mean(iday,ilat) = nanmean(p.nlevs);
           iudef4_mean(iday,ilat) = nanmean(p.iudef(4,:));
           satzen_mean(iday,ilat) = nanmean(p.satzen);
           plevs_mean(iday,ilat,:) = nanmean(p.plevs,2);
           if bRunKlayers
               mmwater_mean(iday,ilat) = nanmean(binwater);
           end
       end  % end loop over latitudes

      iday = iday + 1
   end % if a.bytes > 1000000
end  % giday
% $$$ startdir='/asl/data/stats/airs';
startdir='/home/sbuczko1/WorkingFiles';
eval_str = ['save ' startdir '/rtp_airxbcal_era_rad_kl_'  int2str(year) ...
            '_clear' sDescriptor ' robs rcal rbias_std *_mean count trace '];
eval(eval_str);
