function pull_stats_iasi_dcc(year, filter)

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

addpath /asl/matlab2012/rtptools  % rtpread_12.m
addpath /asl/matlib/h4tools
addpath /asl/rtp_prod/airs/utils
addpath /asl/packages/ccast/motmsc/utils/
addpath ~/git/rtp_prod2/util
addpath ~/git/pull_stats/util
addpath /asl/matlib/aslutil  % mktemp

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

rtpdir = '/asl/rtp/rtp_iasi1/';
if bCheckConfig & isfield(cfg, 'rtpdir')
    rtpdir = cfg.rtpdir;
end

statsdir = '/asl/data/stats/iasi/dcc';
if bCheckConfig & isfield(cfg, 'statsdir')
    statsdir = cfg.statsdir;
end

% collect some system parameters to log
[~, hostname] = system('hostname');
slurm_job_id = getenv('SLURM_JOB_ID');
slurm_array_job_id = getenv('SLURM_ARRAY_JOB_ID')
fprintf(1, '*** Hostname: %s\tJobID: %s\tArray JobID: %s\n', hostname, ...
        slurm_job_id, slurm_array_job_id);
slurm_job_partition = getenv('SLURM_JOB_PARTITION');
slurm_restart_count = getenv('SLURM_RESTART_COUNT');
fprintf(1, '*** Partition: %s\tRestart Count: %s\n', slurm_job_partition, ...
        slurm_restart_count);
slurm_submit_host = getenv('SLURM_SUBMIT_HOST');
slurm_submit_dir = getenv('SLURM_SUBMIT_DIR');
fprintf(1, '*** Submit host: %s\tSubmit dir: %s\n', slurm_submit_host, ...
        slurm_submit_dir);
[sID, sTempPath] = genscratchpath();
fprintf(1, '*** Temp path: %s\tTemp sID: %s\n', sTempPath, sID);
fprintf(1, '*** Task run start %s\n', char(datetime('now')));
fprintf(1, '*** Running filter %d on year %d\n', filter, year);


% record run start datetime in output stats file for tracking
trace.RunDate = datetime('now','TimeZone','local','Format', ...
                         'd-MMM-y HH:mm:ss Z');
trace.Content = 'radiance';
trace.SlurmJobID = slurm_job_id;

basedir = fullfile(rtpdir, 'dcc', int2str(year));
dayfiles = dir(fullfile(basedir, 'iasi1_era_d*_dcc.rtp_1'));
ndays = length(dayfiles);
fprintf(1,'>>> numfiles = %d\n', ndays);

% calculate latitude bins
nbins=20; % gives 2N+1 element array of lat bin boundaries
latbinedges = equal_area_spherical_bands(nbins);
nlatbins = length(latbinedges)-1;

nchans = 8461;  % IASI channel space
nlevs = 101; % klayers output
nfovs = 4; % IASI FOV count

% allocate final accumulator arrays
robs = nan(ndays, nlatbins, nchans, nfovs);
rcal = nan(ndays, nlatbins, nchans, nfovs);

lat_mean = nan(ndays, nlatbins, nfovs);
lon_mean = nan(ndays, nlatbins, nfovs);
solzen_mean = nan(ndays, nlatbins, nfovs);
rtime_mean = nan(ndays, nlatbins, nfovs); 
count = nan(ndays, nlatbins, nchans, nfovs);
stemp_mean = nan(ndays, nlatbins, nfovs);
ptemp_mean = nan(ndays, nlatbins, nlevs, nfovs);
gas1_mean = nan(ndays, nlatbins, nlevs, nfovs);
gas3_mean = nan(ndays, nlatbins, nlevs, nfovs);
spres_mean = nan(ndays, nlatbins, nfovs);
nlevs_mean = nan(ndays, nlatbins, nfovs);
iudef4_mean = nan(ndays, nlatbins, nfovs);
mmwater_mean = nan(ndays, nlatbins, nfovs);
satzen_mean = nan(ndays, nlatbins, nfovs);
satazi_mean = nan(ndays, nlatbins, nfovs);
plevs_mean = nan(ndays, nlatbins, nlevs, nfovs);

iday = 1;
for giday = 1:length(dayfiles)
    fprintf(1, '%s >>> year = %d  :: giday = %d\n', char(datetime('now')), ...
                                                      year, giday);
    if dayfiles(giday).bytes > 100000
        try
            [h,ha,p,pa] = rtpread_12(fullfile(basedir, dayfiles(giday).name));
        catch
            fprintf(2, '%s >>>> ERROR: File issue giday= %d with %s\n', ...
                    char(datetime('now')), giday, dayfiles(giday).name);
            continue;
        end

        % check for proper model levels and move to next day if the
        % wrong model is mixed in
        [levs, obs] = size(p.ptemp);
        if levs ~= 60  % not era in the rtp file
            fprintf(2, '%s >>>> ERROR: model issue giday= %d with %s\n', ...
                    char(datetime('now')), giday, dayfiles(giday).name);
            continue;
        end
        
        f = h.vchan;

        ptype = h.ptype; % grab input ptype for traceability
        
        switch filter
          case 1
            k = find(p.iudef(4,:) == 1); % day
            sDescriptor='day';
          case 2
            k = find(p.iudef(4,:) == 1 & p.landfrac == 0); % day, ocean
            sDescriptor='day_ocean';
          case 3
            k = find(p.iudef(4,:) == 1 & p.landfrac == 1); % day, land
            sDescriptor='day_land';
          case 4
            k = find(p.solzen > 90); % night
            sDescriptor='night';
          case 5
            fprintf(1, '>>> Using solzen for day/night selection <<<\n');
            k = find(p.solzen > 90 & p.landfrac == 0); % night, ocean
            sDescriptor='night_ocean';
          case 6
            k = find(p.solzen > 90 & p.landfrac == 1); % night, land
            sDescriptor='night_land';
        end

        pp = rtp_sub_prof(p, k);
        clear p

        if bRunKlayers
            % klayers kills previous sarta in the rtp structures so
            % we need to save values and re-insert after klayers
            % finishes
            tmp_robs1 = pp.robs1;
            tmp_rcalc = pp.rcalc;

            % run klayers on the rtp data to convert levels ->
            % layers
            fprintf(1, '>>> running klayers... ');
            fn_rtp1 = fullfile(sTempPath, ['iasi_' sID '_1.rtp']);
            outfiles = rtpwrite_12(fn_rtp1, h,ha,pp,pa)
            % run klayers on first half of spectrum
            fprintf(1, '>>> Running klayers on first half of spectrum.\n');
            fbase = ['iasi_' sID '_2.rtp'];
            fn_rtp2 = fullfile(sTempPath, [fbase '_1']);
            klayers_run = [klayers_exec ' fin=' outfiles{1} ' fout=' fn_rtp2 ...
                           ' > ' sTempPath '/kout.txt'];
            unix(klayers_run);
            fprintf(1, '>>>>> Done\n');
            % run klayers on second half of spectrum
            fprintf(1, '>>> Running klayers on second half of spectrum.\n');
            fn_rtp2 = fullfile(sTempPath, [fbase '_2']);
            klayers_run = [klayers_exec ' fin=' outfiles{2} ' fout=' fn_rtp2 ...
                           ' > ' sTempPath '/kout.txt'];
            unix(klayers_run);
            fprintf(1, '>>>>> Done\n');
            fprintf(1, '>>> Reading in klayers output.\n');
            [h,ha,pp,pa] = rtpread_12(fullfile(sTempPath, [fbase '_1']));

            % restore sarta values
            pp.robs1 = tmp_robs1;
            pp.rcalc = tmp_rcalc;
            clear tmp_robs1 tmp_rcalc;
            
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

        % initialize counts and look for bad channels (what should
        % the iasi bad channel test look like?)
        [nchans, nobs] = size(pp.robs1);
        nfovs = 4;
        count_all = int8(ones(nchans, nobs, nfovs));

        % loop over latitude bins
        for ilat = 1:nlatbins
            % subset based on latitude bin
            inbin = find(pp.rlat > latbinedges(ilat) & pp.rlat <= ...
                         latbinedges(ilat+1));
            if isempty(inbin)
                % no obs in latbin. Onward to next
                continue;
            end
            
            p = rtp_sub_prof(pp,inbin);

            for z = 1:4  % loop over FOVs to further sub-select
                ifov = find(p.ifov == z);
                p2 = rtp_sub_prof(p, ifov);

                bincount = count_all(:,inbin,z); 

                % Loop over obs in day
                % Radiance mean and std

                r  = p2.robs1;
                rc = p2.rcalc;

                robs(iday,ilat,:,z) = nanmean(r,2);
                rcal(iday,ilat,:,z) = nanmean(rc,2);
                rbias_std(iday,ilat,:,z) = nanstd(r-rc,0,2);

                lat_mean(iday,ilat,z) = nanmean(p2.rlat);
                lon_mean(iday,ilat,z) = nanmean(p2.rlon);
                solzen_mean(iday,ilat,z) = nanmean(p2.solzen);
                rtime_mean(iday,ilat,z)  = nanmean(p2.rtime);
                count(iday,ilat,z) = sum(bincount(1,:))';
                stemp_mean(iday,ilat,z) = nanmean(p2.stemp);
                iudef4_mean(iday,ilat,z) = nanmean(p2.iudef(4,:));
% $$$                 ptemp_mean(iday,ilat,:,z) = nanmean(p2.ptemp,2);
% $$$                 gas1_mean(iday,ilat,:,z) = nanmean(p2.gas_1,2);
% $$$                 gas3_mean(iday,ilat,:,z) = nanmean(p2.gas_3,2);
                spres_mean(iday,ilat,z) = nanmean(p2.spres);
                nlevs_mean(iday,ilat,z) = nanmean(p2.nlevs);
                satzen_mean(iday,ilat,z) = nanmean(p2.satzen);
% $$$                 plevs_mean(iday,ilat,:,z) = nanmean(p2.plevs,2);

            end  % ifov (z)
        end  % end loop over ilat
            
            iday = iday + 1;
    end % if a.bytes > 1000000
end  % giday

savefile = sprintf('/asl/data/stats/iasi/dcc/rtp_iasi_era_%d_rad_dcc_%s', year, sDescriptor);
save(savefile, 'robs','rcal', 'rbias_std', '*_mean','count', ...
     'latbinedges', 'trace', '-v7.3')

fprintf(1, '*** Task end time: %s\n', char(datetime('now')));
