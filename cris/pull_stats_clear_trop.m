function pull_stats_clear_trop(year);

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
addpath /asl/matlib/rtptools  % mmwater_rtp.m

klayers_exec = ['/asl/packages/klayersV205/BinV201/' ...
                'klayers_airs_wetwater'];

[sID, sTempPath] = genscratchpath();

basedir = fullfile('/asl/rtp/rtp_cris_UW_hires/clear', ...
                   int2str(year));
dayfiles = dir(fullfile(basedir, 'uwcris_ecmwf_20180121_clear.rtp'));

iday = 1;
% for giday = 1:50:length(dayfiles)
for giday = 1:length(dayfiles)
    fprintf(1, '>>> year = %d  :: giday = %d\n', year, giday);
    a = dir(fullfile(basedir,dayfiles(giday).name));
    if a.bytes > 100000
        [h,ha,p,pa] = rtpread(fullfile(basedir,dayfiles(giday).name));

        % make initial selection for night ocean observations
% $$$         k = find(p.iudef(4,:) == 1 & p.landfrac == 0); % descending node
% $$$                                                        % (night), ocean
        k = find(p.solzen > 90 & p.landfrac == 0); % descending node
                                                       % (night), ocean
        sDescriptor='_desc_ocean';
        pp = rtp_sub_prof(p, k);

      % run klayers on the rtp data (Sergio is asking for this to
      % convert levels to layers for his processing?)

      % first remove rcalc field and save it for later restore
      rcalc = pp.rclr;
      pp = rmfield(pp, 'rclr');
      
      fprintf(1, '>>> running klayers... ');
      fn_rtp1 = fullfile(sTempPath, ['airs_' sID '_1.rtp']);
      rtpwrite(fn_rtp1, h,ha,pp,pa)
      fn_rtp2 = fullfile(sTempPath, ['airs_' sID '_2.rtp']);
      klayers_run = [klayers_exec ' fin=' fn_rtp1 ' fout=' fn_rtp2 ...
                     ' > ' sTempPath '/kout.txt'];
      unix(klayers_run);
      [h,ha,pp,pa] = rtpread(fn_rtp2);
      % restore rcalc
      pp.rclr = rcalc;
      clear rcalc;
      
      fprintf(1, 'Done\n');

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
% $$$         nfovs = 9;
% $$$         count_all = int8(ones(nchans, nobs, nfovs));
        nfors = 30;
        
        % subset based on latitude bin for tropical (+- 30)
        inbin = find(pp.rlat > -30 & pp.rlat < 30);
        p = rtp_sub_prof(pp,inbin);
        sDescriptor = [sDescriptor '_trop'];
        
        for z = 1:9  % loop over FOVs to further sub-select
            ifov = find(p.ifov == z);
            p2 = rtp_sub_prof(p, ifov);
% $$$         for z = 1:30  % loop over FORs to further sub-select
% $$$             ifor = find(p.xtrack == z);
% $$$             p2 = rtp_sub_prof(p, ifor);
            
% $$$             bincount = int8(ones(nchans, length(ifor)));
            bincount = int8(ones(nchans, length(ifov)));                            
            
            % Loop over obs in day
            % Radiance mean and std
            
            r  = p2.robs1;
            rc = p2.rclr;
            
            % Convert r to rham
% $$$             r = box_to_ham(r);  % assumes r in freq order!!  Needed
% $$$                                 % for lowres
            bias_std(iday,:,z) = nanstd(r-rc,0,2);
            robs(iday,:,z) = nanmean(r,2);
            rcal(iday,:,z) = nanmean(rc,2);
            lat_mean(iday,z) = nanmean(p2.rlat);
            lon_mean(iday,z) = nanmean(p2.rlon);
            solzen_mean(iday,z) = nanmean(p2.solzen);
            rtime_mean(iday,z)  = nanmean(p2.rtime);
            count(iday,z) = sum(bincount(1,:))';
            stemp_mean(iday,z) = nanmean(p2.stemp);
            iudef4_mean(iday,z) = nanmean(p2.iudef(4,:));
            ptemp_mean(iday,:,z) = nanmean(p2.ptemp,2);
            gas1_mean(iday,:,z) = nanmean(p2.gas_1,2);
            gas2_mean(iday,:,z) = nanmean(p2.gas_2,2);                
            gas3_mean(iday,:,z) = nanmean(p2.gas_3,2);
            gas4_mean(iday,:,z) = nanmean(p2.gas_4,2);
            gas5_mean(iday,:,z) = nanmean(p2.gas_5,2);
            gas6_mean(iday,:,z) = nanmean(p2.gas_6,2);
            gas9_mean(iday,:,z) = nanmean(p2.gas_9,2);
            gas12_mean(iday,:,z) = nanmean(p2.gas_12,2);                    
            spres_mean(iday,z) = nanmean(p2.spres);
            nlevs_mean(iday,z) = nanmean(p2.nlevs);
            satzen_mean(iday,z) = nanmean(p2.satzen);
            plevs_mean(iday,:,z) = nanmean(p2.plevs,2);
            scanang_mean(iday,z) = nanmean(p2.scanang);
        end  % ifov (z)
        
        iday = iday + 1
    end % if a.bytes > 1000000
end  % giday
eval_str = ['save /asl/data/stats/cris2/rtp_UWcris_hires_021_binFOV_'  int2str(year) ...
            '_clear' sDescriptor ' robs rcal bias_std *_mean count '];

eval(eval_str);
