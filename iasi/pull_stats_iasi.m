function pull_stats_iasi(year, filter)

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
addpath /home/sergio/MATLABCODE/PLOTTER  %
                                         % equal_area_spherical_bands
addpath /asl/matlib/aslutil  % mktemp

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

klayers_exec = '/asl/packages/klayersV205/BinV201/klayers_airs_wetwater';

basedir = fullfile('/asl/rtp/rtp_iasi1/clear', int2str(year));

% record run start datetime in output stats file for tracking
trace.RunDate = datetime('now','TimeZone','local','Format', ...
                         'd-MMM-y HH:mm:ss Z');
trace.Content = 'radiance';
trace.SlurmJobID = slurm_job_id;

% Chris built the IASI rtp generation such that output rtp files go
% into a YYYY/MM/DD/*.rtp_[12] tree. Each leaf branch constitutes
% one day and has two rtp files, each containing about half the
% IASI spectral space. For a given year, then, we'll grab a list of
% MM directories to loop over. Within each of those, we will do the
% same to grab a list of DD directories. I will leave this
% abstracted out separately so that changes will be minimal if we
% reorganize the rtp directories to conform to airs/cris YYYY/DDD
% style

monthdirs = dir(fullfile(basedir, '*'));
i=1;  % day counter
dayfiles=monthdirs(1);  % initialize structure

% the results of dir derived this way include the directories '.'
% and '..' in elements 1&2. interesting (valid) month and day
% directories will start at index 3
fprintf(1, '%s Starting file search...\n', char(datetime('now')));
for imonth=3:numel(monthdirs)
    mbasedir = fullfile(basedir, monthdirs(imonth).name);
    daydirs = dir(fullfile(mbasedir, '*'));
    for iday = 3:numel(daydirs)
        dbasedir = fullfile(mbasedir, daydirs(iday).name);
        sdayfiles = dir(fullfile(dbasedir, '*.rtp_*'));
        % check for existence of daily files before using them for
        % processing
        switch numel(sdayfiles)
          case 2
            % two files as expected. continue
            dayfiles(i) = sdayfiles(1);
            % dayfiles is collecting files info from disparate
            % paths so, capture full path info in the name field
            dayfiles(i).name = fullfile(dbasedir, sdayfiles(1).name);
            i=i+1;
          case 1
            % only one of two files. report and go to next day
            fprintf(2, ['%s >>> ERROR: %s contains only one rtp file. ' ...
                        'Moving on to next day.\n'], ...
                    char(datetime('now')), dbasedir);
            continue;
          case 0
            % no files found for day. report and go to next day
            fprintf(2, ['%s >>> ERROR: %s contains no rtp files. Moving ' ...
                        'on to next day.\n'], char(datetime('now')), ...
                dbasedir);
            continue;
        end
        
    end
end
fprintf(1, '%s File search complete.\n', char(datetime('now')));

% calculate latitude bins
nbins=20; % gives 2N+1 element array of lat bin boundaries
latbins = equal_area_spherical_bands(nbins);
nlatbins = length(latbins);

iday = 1;

for giday = 1:length(dayfiles)
    fprintf(1, '%s >>> year = %d  :: giday = %d\n', char(datetime('now')), ...
                                                      year, giday);
    if dayfiles(giday).bytes > 100000
        try
            [h,ha,p,pa] = rtpread_12(dayfiles(giday).name);
        catch
            fprintf(2, '%s >>>> ERROR: File issue giday= %d with %s\n', ...
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

        % run klayers on the rtp data (Sergio is asking for this to
        % convert levels to layers for his processing?)
        % *** Actually, Chris is keeping klayers results in his rtp
        % output files so, I can just use the values already there ***
% $$$ 
% $$$         % first remove rcalc field and save it for later restore
% $$$         rcalc = p.rcalc;
% $$$         p = rmfield(p, 'rcalc');
% $$$         
% $$$         h.ptype = 0; % force ptype to LEVPRO so klayers will run (even if it
% $$$                      % was run in rtp generation)
% $$$ 
% $$$         fprintf(1, '>>> running klayers... ');
% $$$         fn_rtp1 = fullfile(sTempPath, ['iasi_' sID '_1.rtp']);
% $$$         outfiles = rtpwrite_12(fn_rtp1, h,ha,p,pa)
% $$$         % run klayers on first half of spectrum
% $$$         fprintf(1, '>>> Running klayers on first half of spectrum.\n');
% $$$         fbase = ['iasi_' sID '_2.rtp'];
% $$$         fn_rtp2 = fullfile(sTempPath, [fbase '_1']);
% $$$         klayers_run = [klayers_exec ' fin=' outfiles{1} ' fout=' fn_rtp2 ...
% $$$                        ' > ' sTempPath '/kout.txt'];
% $$$         unix(klayers_run);
% $$$         fprintf(1, '>>>>> Done\n');
% $$$         % run klayers on second half of spectrum
% $$$         fprintf(1, '>>> Running klayers on second half of spectrum.\n');
% $$$         fn_rtp2 = fullfile(sTempPath, [fbase '_2']);
% $$$         klayers_run = [klayers_exec ' fin=' outfiles{2} ' fout=' fn_rtp2 ...
% $$$                        ' > ' sTempPath '/kout.txt'];
% $$$         unix(klayers_run);
% $$$         fprintf(1, '>>>>> Done\n');
% $$$         fprintf(1, '>>> Reading in klayers output.\n');
% $$$         [h,ha,p,pa] = rtpread_12(fullfile(sTempPath, [fbase '_1']));
% $$$         % restore rcalc
% $$$         p.rcalc = rcalc;
% $$$         clear rcalc;
% $$$
% $$$         ptype = 1; % reset ptype to LAYPRO since we've re-run klayers
% $$$         
% $$$         fprintf(1, 'Done\n');

        % initialize counts and look for bad channels (what should
        % the iasi bad channel test look like?)
        [nchans, nobs] = size(pp.robs1);
        nfovs = 4;
        count_all = int8(ones(nchans, nobs, nfovs));

        % grab levels/layer indication for traceability from most
        % current head struct
        trace.ptype = ptype;

        % loop over latitude bins
        for ilat = 1:nlatbins-1
            % subset based on latitude bin
            inbin = find(pp.rlat > latbins(ilat) & pp.rlat <= ...
                         latbins(ilat+1));
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
                ptemp_mean(iday,ilat,:,z) = nanmean(p.ptemp,2);
                gas1_mean(iday,ilat,:,z) = nanmean(p.gas_1,2);
                gas2_mean(iday,ilat,:,z) = nanmean(p.gas_2,2);                
                gas3_mean(iday,ilat,:,z) = nanmean(p.gas_3,2);
                gas4_mean(iday,ilat,:,z) = nanmean(p.gas_4,2);
                gas5_mean(iday,ilat,:,z) = nanmean(p.gas_5,2);
                gas6_mean(iday,ilat,:,z) = nanmean(p.gas_6,2);
                gas9_mean(iday,ilat,:,z) = nanmean(p.gas_9,2);
                gas12_mean(iday,ilat,:,z) = nanmean(p.gas_12,2);                    
                spres_mean(iday,ilat,z) = nanmean(p.spres);
                nlevs_mean(iday,ilat,z) = nanmean(p.nlevs);
                satzen_mean(iday,ilat,z) = nanmean(p.satzen);
                plevs_mean(iday,ilat,:,z) = nanmean(p.plevs,2);
                scanang_mean(iday,ilat,z) = nanmean(p.scanang);
            end  % ifov (z)
        end  % end loop over ilat
            
            iday = iday + 1;
    end % if a.bytes > 1000000
end  % giday

savefile = sprintf('/asl/data/stats/iasi/rtp_iasi_%d_rad_clear_%s', year, sDescriptor);
save(savefile, 'robs','rcal', 'rbias_std', '*_mean','count', 'trace')

fprintf(1, '*** Task end time: %s\n', char(datetime('now')));