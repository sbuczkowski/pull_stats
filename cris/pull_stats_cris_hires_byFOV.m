function pull_stats_cris_hires_byFOV(year, ifov, cfg);

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
func_name = 'pull_stats_cris_hires'
%year = 2014;

addpath /asl/matlib/h4tools
addpath /asl/rtp_prod/airs/utils
% $$$ addpath /asl/packages/ccast/motmsc/utils/
addpath /home/sbuczko1/git/rtp_prod2/util
addpath /home/sbuczko1/git/pull_stats/util  % cat_rtp_basic
addpath /asl/rtp_prod/cris/unapod
addpath /asl/matlib/aslutil   % int2bits
addpath /home/sergio/MATLABCODE/PLOTTER  %
                                         % equal_area_spherical_bands
addpath /home/sbuczko1/git/rtp_prod2/cris % cris_lowres_chans
addpath /asl/matlib/rtptools  % mmwater_rtp.m
addpath('/home/sbuczko1/git/swutils');  % githash


[sID, sTempPath] = genscratchpath();

% check for existence of configuration struct
bCheckConfig = false;
if nargin == 3
    bCheckConfig = true;
end

% record run start datetime in output stats file for tracking
trace.githash = githash(func_name);
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

rtpdir = '/asl/rtp/rtp_cris_ccast_hires/';
if bCheckConfig & isfield(cfg, 'rtpdir')
    rtpdir = cfg.rtpdir;
end

statsdir = '/asl/data/stats/cris';
if bCheckConfig & isfield(cfg, 'statsdir')
    statsdir = cfg.statsdir;
end

sPattern='rtp_cris_hires_%4d_FOV-%d_%s_%s';
if bCheckConfig & isfield(cfg, 'spattern')
    sPattern = cfg.spattern;
end

svars = ' robs rclr rbias_std *_mean count trace latbinedges';
if bCheckConfig & isfield(cfg, 'svars')
    svars = cfg.svars;
end

filter = 1
if bCheckConfig & isfield(cfg, 'filter')
    filter = cfg.filter;
end

sSubset = 'allfov';
fprintf(1, '>> Running allfov by FOV stats\n');
basedir = fullfile(rtpdir, sSubset, int2str(year));
gdays = dir(basedir);
gdays(1:2) = [];  % dir picks up . and .. directories. Remove them
                  % from array (assumes dir always puts them
                  % first. DANGEROUS)
ndays = length(gdays);
fprintf(1, '>> Found %d days to process\n', ndays);

% Get proper frequencies for these data
[n1,n2,n3,userLW,userMW,userSW, ichan] = cris_hires_chans();
f = cris_vchan(2, userLW, userMW, userSW);
nchans = length(f);

% calculate latitude bins
% $$$ nbins=20; % gives 2N+1 element array of lat bin boundaries
% $$$ latbinedges = equal_area_spherical_bands(nbins);

% single +-50 lat bin
latbinedges = [-50 50];

nlatbins = length(latbinedges)-1;

nlevs = 101;  % klayers output
nfors = 45;   % FORs 
nfovs = 9;    % FOVs/FOR

% $$$ for ifov = 1:nfovs  % loop over FOVs 

    % allocate final accumulator arrays
    robs = zeros(ndays, nfors, nlatbins, nchans);
    rclr = zeros(ndays, nfors, nlatbins, nchans);
    rcld = zeros(ndays, nfors, nlatbins, nchans);
    rcldbias_std = zeros(ndays, nfors, nlatbins, nchans);
    rclrbias_std = zeros(ndays, nfors, nlatbins, nchans);

    lat_mean = zeros(ndays, nfors, nlatbins);
    lon_mean = zeros(ndays, nfors, nlatbins);
    solzen_mean = zeros(ndays, nfors, nlatbins);
    rtime_mean = zeros(ndays, nfors, nlatbins); 
    count = zeros(ndays, nfors, nlatbins, nchans);
    tcc_mean = zeros(ndays, nfors, nlatbins);
    stemp_mean = zeros(ndays, nfors, nlatbins);
    ptemp_mean = zeros(ndays, nfors, nlatbins, nlevs);
    gas1_mean = zeros(ndays, nfors, nlatbins, nlevs);
    gas3_mean = zeros(ndays, nfors, nlatbins, nlevs);
    spres_mean = zeros(ndays, nfors, nlatbins);
    nlevs_mean = zeros(ndays, nfors, nlatbins);
    iudef4_mean = zeros(ndays, nfors, nlatbins);
    mmwater_mean = zeros(ndays, nfors, nlatbins);
    satzen_mean = zeros(ndays, nfors, nlatbins);
    plevs_mean = zeros(ndays, nfors, nlatbins, nlevs);

    iday = 1;
    for giday = 1:ndays
        fprintf(1, '>>> year = %d  :: giday = %d\n', year, giday);
        gbasedir = fullfile(basedir, gdays(giday).name);
        fprintf(1, '>> looking for input granule files in %s\n', ...
                gbasedir);
        grans = dir(sprintf('%s/%s', gbasedir, '*.rtp'));
        ngrans = length(grans);
        fprintf(1,'>>> number of granules found = %d\n', ngrans);
        
        bFirstGranRead = false;
        for igran = 1:ngrans
            gfile = fullfile(gbasedir,grans(igran).name);
            fprintf(1, '>> Reading granule %d : %s\n', igran, gfile)
            [head,hattr,prof,pattr] = rtpread(gfile);
            f = head.vchan;

            % sanity check on p.robs1 as read in. (There have been
            % instances where this array is short on the spectral
            % dimension which fails in rad2bt. We trap for this here)
            obs = size(prof.robs1);
            chans = size(f);
            if obs(1) ~= chans(1)
                fprintf(2, ['**>> ERROR: obs/vchan spectral channel ' ...
                            'mismatch in %s. Bypassing day.\n'], dayfiles(giday).name);
                continue;
            end

            switch filter
              case 1
                k = find(prof.iudef(4,:) == 1); % descending node (night)
                sDescriptor='desc';
              case 2
                k = find(prof.iudef(4,:) == 1 & prof.landfrac == 0); % descending node
                                                                     % (night), ocean
                sDescriptor='desc_ocean';
              case 3
                k = find(prof.iudef(4,:) == 1 & prof.landfrac == 1); % descending node
                                                                     % (night), land
                sDescriptor='desc_land';
              case 4
                k = find(prof.iudef(4,:) == 0); % ascending node (day)
                sDescriptor='asc';
              case 5
                k = find(prof.iudef(4,:) == 0 & prof.landfrac == 0); % ascending node
                                                                     % (day), ocean
                sDescriptor='asc_ocean';
              case 6
                k = find(prof.iudef(4,:) == 0 & prof.landfrac == 1); % ascending node
                                                                     % (day), land
                sDescriptor='asc_land';
            end

            infov = find(prof.ifov == ifov);
            keep_inds = intersect(k, infov);        

            pp = rtp_sub_prof(prof, keep_inds);
            clear prof

            fprintf(1, '>>> Found %d filtered obs in granule %d\n', ...
                    length(pp.rtime), igran);
            if length(pp.rtime) == 0
                % no obs in current day for current filter. jump to next
                % day
                fprintf(2, ['>>> No obs in gran %d for current filter. NEXT ' ...
                            'GRANULE\n'], igran);
                continue;
            end

            % run klayers on the rtp data (Sergio is asking for this to
            % convert levels to layers for his processing?)

            % first remove rcalc field and save it for later restore
            tmp_rclr = pp.rclr;
            pp = rmfield(pp, 'rclr');
            
            fprintf(1, '>>> running klayers... ');
            fn_rtp1 = fullfile(sTempPath, ['airs_' sID '_1.rtp']);
            rtpwrite(fn_rtp1, head,hattr,pp,pattr)
            fn_rtp2 = fullfile(sTempPath, ['airs_' sID '_2.rtp']);
            klayers_run = [klayers_exec ' fin=' fn_rtp1 ' fout=' fn_rtp2 ...
                           ' > ' sTempPath '/kout.txt'];
            unix(klayers_run);
            [head,hattr,pp,pattr] = rtpread(fn_rtp2);
            % restore rcalc
            pp.rclr = tmp_rclr;
            clear tmp_rclr;
            
            fprintf(1, 'Done\n');

            if ~bFirstGranRead
                h = head;
                ha = hattr;
                p = pp;
                pa = pattr;

                bFirstGranRead = true;  % first granule has been read
            else
                [h, p] = cat_rtp_basic(h, p, head, pp);
            end
            
        end % end loop over granules
            clear prof

            % get column water
            mmwater = mmwater_rtp(h, p);

            % Check for obs with layer profiles that go lower than
            % topography. Need to check nlevs and NaN out any layers
            % at or below this level

            % ** Any layers-sensitive variables added in averaging code below must
            % ** be checked here first.
            for i=1:length(p.nlevs)
                badlayers = p.nlevs(i) : 101;
                p.plevs(badlayers, i) = NaN;
                p.gas_1(badlayers, i) = NaN;
                p.gas_3(badlayers, i) = NaN;
                p.ptemp(badlayers, i) = NaN;
            end

            % initialize counts and look for bad channels (what should
            % the iasi bad channel test look like?)
            [nchans, nobs] = size(p.robs1);

            % loop over FORs
            for ifor = 1:nfors
                infor = find(p.atrack == ifor);
                % loop over latitude bins
                for ilat = 1:nlatbins
                    % subset based on latitude bin
                    inbin = find(p.rlat > latbinedges(ilat) & p.rlat <= ...
                                 latbinedges(ilat+1));
                    keep_ind = intersect(infor, inbin);
                    pp = rtp_sub_prof(p,keep_ind);

                    count_infov = ones(length(pp.rlat), nchans);
                    binwater = mmwater(inbin);

                    % Convert r to rham
                    r = box_to_ham(pp.robs1);  % assumes r in freq order!!

                    robs(iday,ifor,ilat,:) = nanmean(r,2);
                    rclr(iday,ifor,ilat,:) = nanmean(pp.rclr,2);
                    rbias_std(iday,ifor,ilat,:) = nanstd(r-pp.rclr,0,2);
                    
                    lat_mean(iday,ifor,ilat) = nanmean(pp.rlat);
                    lon_mean(iday,ifor,ilat) = nanmean(pp.rlon);
                    solzen_mean(iday,ifor,ilat) = nanmean(pp.solzen);
                    rtime_mean(iday,ifor,ilat)  = nanmean(pp.rtime);
                    count(iday,ifor,ilat,:) = sum(count_infov);
                    stemp_mean(iday,ifor,ilat) = nanmean(pp.stemp);
                    iudef4_mean(iday,ifor,ilat) = nanmean(pp.iudef(4,:));
                    ptemp_mean(iday,ifor,ilat,:) = nanmean(pp.ptemp,2);
                    gas1_mean(iday,ifor,ilat,:) = nanmean(pp.gas_1,2);
                    gas3_mean(iday,ifor,ilat,:) = nanmean(pp.gas_3,2);
                    spres_mean(iday,ifor,ilat) = nanmean(pp.spres);
                    nlevs_mean(iday,ifor,ilat) = nanmean(pp.nlevs);
                    satzen_mean(iday,ifor,ilat) = nanmean(pp.satzen);
                    plevs_mean(iday,ifor,ilat,:) = nanmean(pp.plevs,2);
                    mmwater_mean(iday,ifor,ilat) = nanmean(binwater);
                end  % latbins
            end  % FORs
            iday = iday + 1

    end  % iday
    
    % save all days to single yearly file per fov
    sname = sprintf(sPattern, year,ifov,sSubset,sDescriptor);

    eval_str = sprintf('save %s/%s %s', statsdir, sname, svars);
    fprintf(1, '>> Executing save command: %s\n', eval_str)
    eval(eval_str);

% $$$ end  % ifov


