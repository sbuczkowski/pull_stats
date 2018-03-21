function pull_stats_cris_clear_fs_scanangle(year, filter, cfg);

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
addpath /home/sbuczko1/git/rtp_prod2/cris  % cris_lowres_chans
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

rtpdir = '/asl/rtp/rtp_cris_ccast_lowres/';
if bCheckConfig & isfield(cfg, 'rtpdir')
    rtpdir = cfg.rtpdir;
end

statsdir = '/asl/data/stats/cris';
if bCheckConfig & isfield(cfg, 'statsdir')
    statsdir = cfg.statsdir;
end

basedir = fullfile(rtpdir, 'clear', int2str(year));
dayfiles = dir(fullfile(basedir, 'cris_lr_era_d*_clear.rtp'));
fprintf(1,'>>> numfiles = %d\n', length(dayfiles));

% Get proper frequencies for these data
[n1,n2,n3,userLW,userMW,userSW, ichan] = cris_lowres_chans();
f = cris_vchan(2, userLW, userMW, userSW);

nfors = 30;  % number of FORs cross-track
nfovs = 9;   % number of FOVs per FOR

iday = 1;
for giday = 1:length(dayfiles)
% $$$ for giday = 18:20
    fprintf(1, '>>> year = %d  :: giday = %d\n', year, giday);
    a = dir(fullfile(basedir,dayfiles(giday).name));
    if a.bytes > 100000
        [h,ha,p,pa] = rtpread(fullfile(basedir,dayfiles(giday).name));

      switch filter
        case 1
          k = find(p.solzen > 90); % descending node (night)
          sDescriptor='desc';
        case 2
          k = find(p.solzen > 90 & p.landfrac == 0); % descending node
                                                         % (night), ocean
          sDescriptor='desc_ocean';
        case 3
          k = find(p.solzen > 90 & p.landfrac == 1); % descending node
                                                        % (night), land
          sDescriptor='desc_land';
        case 4
          k = find(p.solzen <= 90); % ascending node (day)
          sDescriptor='asc';
        case 5
          k = find(p.solzen <= 90 & p.landfrac == 0); % ascending node
                                                         % (day), ocean
          sDescriptor='asc_ocean';
        case 6
          k = find(p.solzen <= 90 & p.landfrac == 1); % ascending node
                                                        % (day), land
          sDescriptor='asc_land';
      end

      pp = rtp_sub_prof(p, k);

      % limit obs collected to lat +-60
      p = pp;
      k = find(abs(p.rlat) < 60);
      pp = rtp_sub_prof(p, k);
      
      % run klayers on the rtp data to convert levels -> layers
      % save calcs as the re-run of klayers wipes them out
      rclr = pp.rcalc;
      fprintf(1, '>>> running klayers... ');
      fn_rtp1 = fullfile(sTempPath, ['cris_' sID '_1.rtp']);
      rtpwrite(fn_rtp1, h,ha,pp,pa);
      clear pp;
      fn_rtp2 = fullfile(sTempPath, ['cris_' sID '_2.rtp']);
      klayers_run = [klayers_exec ' fin=' fn_rtp1 ' fout=' fn_rtp2 ...
                     ' > ' sTempPath '/kout.txt'];
      unix(klayers_run);
      fprintf(1, 'Done\n');

      % Read klayers output into local rtp variables
      [h,ha,pp,pa] = rtpread(fn_rtp2);
      % restore sarta_rclearcalc
      pp.rcalc = rclr;
      clear rclr
      
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
      for ifor = 1:nfors
          % subset based on latitude bin
          inbin = find(pp.xtrack == ifor);
          p = rtp_sub_prof(pp,inbin);
          
          for z = 1:9  % loop over FOVs to further sub-select
              ifov = find(p.ifov == z);
              p2 = rtp_sub_prof(p, ifov);
              
              bincount = count_all(:,inbin,z); 
              binwater = mmwater(inbin);        
              % Loop over obs in day
              % Radiance mean and std
              
              r  = p2.robs1;
              rc = p2.rcalc;
              
              % leave as sinc for test
              % Convert r to rham
              r = box_to_ham(r);  % assumes r in freq order!!  Needed
                                  % for lowres
              
              robs(iday,ifor,:,z) = nanmean(r,2);
              rcal(iday,ifor,:,z) = nanmean(rc,2);
              rbias_std(iday, ifor,:) = nanstd(r-rc,0,2);
              
              lat_mean(iday,ifor,z) = nanmean(p2.rlat);
              lon_mean(iday,ifor,z) = nanmean(p2.rlon);
              solzen_mean(iday,ifor,z) = nanmean(p2.solzen);
              rtime_mean(iday,ifor,z)  = nanmean(p2.rtime);
              count(iday,ifor,z) = sum(bincount(1,:))';
              stemp_mean(iday,ifor,z) = nanmean(p2.stemp);
              iudef4_mean(iday,ifor,z) = nanmean(p2.iudef(4,:));
              ptemp_mean(iday,ifor,:,z) = nanmean(p2.ptemp,2);
              gas1_mean(iday,ifor,:,z) = nanmean(p2.gas_1,2);
              gas3_mean(iday,ifor,:,z) = nanmean(p2.gas_3,2);
              spres_mean(iday,ifor,z) = nanmean(p2.spres);
              nlevs_mean(iday,ifor,z) = nanmean(p2.nlevs);
              satzen_mean(iday,ifor,z) = nanmean(p2.satzen);
              plevs_mean(iday,ifor,:,z) = nanmean(p2.plevs,2);
              mmwater_mean(iday,ifor) = nanmean(binwater);
% $$$               scanang_mean(iday,ifor,z) = nanmean(p2.scanang);
          end  % ifov (z)
      end  % end loop over ifor
          
          iday = iday + 1
   end % if a.bytes > 1000000
end  % giday
outfile = fullfile(statsdir, sprintf('rtp_cris_lowres_era_rad_kl_%s_clear_scanang_%s', ...
           int2str(year), sDescriptor));
eval_str = ['save ' outfile ' robs rcal *_std *_mean count '];
eval(eval_str);
