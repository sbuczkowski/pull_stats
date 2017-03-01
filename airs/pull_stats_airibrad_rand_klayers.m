function pull_stats_airibrad_rand_klayers(year, filter);

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
addpath ~/git/rtp_prod2/util
addpath /home/sergio/MATLABCODE/PLOTTER  %
                                         % equal_area_spherical_bands
addpath /asl/matlib/rtptools  % mmwater_rtp.m

% record run start datetime in output stats file for tracking
trace.RunDate = datetime('now','TimeZone','local','Format', ...
                         'd-MMM-y HH:mm:ss Z');

klayers_exec = ['/asl/packages/klayersV205/BinV201/' ...
                'klayers_airs_wetwater'];
cstr =[ 'bits1-4=NEdT[0.08 0.12 0.15 0.20 0.25 0.30 0.35 0.4 0.5 0.6 0.7' ...
        ' 0.8 1.0 2.0 4.0 nan]; bit5=Aside[0=off,1=on]; bit6=Bside[0=off,1=on];' ...
        ' bits7-8=calflag&calchansummary[0=OK, 1=DCR, 2=moon, 3=other]' ];

[sID, sTempPath] = genscratchpath();

basedir = fullfile('/asl/rtp/rtp_airibrad_v5/', ...
                   int2str(year), 'random');
dayfiles = dir(fullfile(basedir, 'era_airibrad*_random.rtp'));
ndayfiles = length(dayfiles);
fprintf(1,'>>> numfiles = %d\n', ndayfiles);

% calculate latitude bins
nbins=20; % gives 2N+1 element array of lat bin boundaries
latbins = equal_area_spherical_bands(nbins);
nlatbins = length(latbins);

% initialize cloud summary arrays with -9999 for Sergio
cfrac_mean = ones(ndayfiles,nlatbins) * -9999;
cprtop_mean = cfrac_mean;
cprbot_mean = cfrac_mean;
cngwat_mean = cfrac_mean;
cpsize_mean = cfrac_mean;
cfrac2_mean = cfrac_mean;
cprtop2_mean = cfrac_mean;
cprbot2_mean = cfrac_mean;
cngwat2_mean = cfrac_mean;
cpsize2_mean = cfrac_mean;
cfrac12_mean = cfrac_mean;

iday = 1;
% for giday = 1:50:ndayfiles
for giday = 1:ndayfiles
    fprintf(1, '>>> year = %d  :: giday = %d\n', year, giday);
    a = dir(fullfile(basedir,dayfiles(giday).name));
    if a.bytes > 100000
        [h,ha,pp,pa] = rtpread(fullfile(basedir, ...
                                       dayfiles(giday).name));

        % subset for basic type and reduce load going in to klayers
        switch filter
          case 1
            k = find(pp.iudef(4,:) == 68); % descending node (night)
            sDescriptor='_desc';
          case 2
            k = find(pp.iudef(4,:) == 68 & pp.landfrac == 0); % descending
                                                            % node (night) ocean
            sDescriptor='_desc_ocean';
          case 3
            k = find(pp.iudef(4,:) == 68 & pp.landfrac == 1); % descending node
                                                            % (night), land
            sDescriptor='_desc_land';
          case 4
            k = find(pp.iudef(4,:) == 65); % ascending node (day)
            sDescriptor='_asc';
          case 5
            k = find(pp.iudef(4,:) == 65 & pp.landfrac == 0); % ascending node
                                                            % (day), ocean
            sDescriptor='_asc_ocean';
          case 6
            k = find(pp.iudef(4,:) == 65 & pp.landfrac == 1); % ascending node
                                                            % (day), land
            sDescriptor='_asc_land';
        end

        p = rtp_sub_prof(pp, k);
        clear pp;
        
        % save some sarta output from rtp generation as following
        % klayers run will wipe it out
        sarta_rclearcalc = p.sarta_rclearcalc;
        sarta_rcldycalc = p.rcalc;
        
        f = h.vchan;  % AIRS proper frequencies

        % sanity check on p.robs1 as read in. (There have been
        % instances where this array is short on the spectral
        % dimension which fails in rad2bt. We trap for this here)
        obs = size(p.robs1);
        chans = size(f);
        if obs(1) ~= chans(1)
            fprintf(2, ['**>> ERROR: obs/vchan spectral channel ' ...
                        'mismatch in %s. Bypassing day.\n'], dayfiles(giday).name);
            continue;
        end
        
        % run klayers on the rtp data (Sergio is asking for this to
        % convert levels to layers for his processing?)
        fprintf(1, '>>> running klayers... ');
        fn_rtp1 = fullfile(sTempPath, ['airs_' sID '_1.rtp']);
        rtpwrite(fn_rtp1, h,ha,p,pa);
        clear p;
        fn_rtp2 = fullfile(sTempPath, ['airs_' sID '_2.rtp']);
        klayers_run = [klayers_exec ' fin=' fn_rtp1 ' fout=' fn_rtp2 ...
                       ' > ' sTempPath '/kout.txt'];
        unix(klayers_run);
        [h,ha,pp,pa] = rtpread(fn_rtp2);
        fprintf(1, 'Done\n');

        % restore relevant sarta output from original rtp structs
        pp.rcalc = sarta_rclearcalc;
        pp.rcldycalc = sarta_rcldycalc;
        clear sarta_rclearcalc sarta_rcldycalc;
        
        % get column water
        mmwater = mmwater_rtp(h, pp);
        
        % Initialize counts
        n = length(pp.rlat);
        count_all = ones(2378,n);

        % Find bad channels and clear them from spectra
        [nedt,ab,ical] = calnum_to_data(pp.calflag,cstr);
        for i=1:2378
            k = find( pp.robs1(i,:) == -9999 | ical(i,:) ~= 0 | nedt(i,:) > 1);
            %          % These are the good channels
            %          kg = setdiff(1:n,k);
            % NaN's for bad channels
            pp.robs1(i,k) = NaN;
            pp.rcalc(i,k) = NaN;
            pp.rcldycalc(i,k) = NaN;
            count_all(i,k) = 0;
        end

        % Loop over latitude bins
        for ilat = 1:nlatbins-1
            % subset based on latitude bin
            inbin = find(pp.rlat > latbins(ilat) & pp.rlat <= ...
                         latbins(ilat+1));
            p = rtp_sub_prof(pp,inbin);
            bincount = count_all(:,inbin); 
            binwater = mmwater(inbin);
            
            % Collect output stats for the bin
            r  = p.robs1;
            rc = p.rcalc;
            rcldycalc = p.rcldycalc;

            robs(iday,ilat,:) = nanmean(r,2);
            rcldy(iday,ilat,:) = nanmean(rcldycalc,2);
            rclr(iday,ilat,:) = nanmean(rc,2);
            bias_std(iday,ilat,:) = nanstd(r-rc,0,2);
            count(iday,ilat,:) = sum(bincount,2)';
            lat_mean(iday,ilat) = nanmean(p.rlat);
            lon_mean(iday,ilat) = nanmean(p.rlon);
            solzen_mean(iday,ilat) = nanmean(p.solzen);
            rtime_mean(iday,ilat)  = nanmean(p.rtime);
            stemp_mean(iday,ilat) = nanmean(p.stemp);
            ptemp_mean(iday,ilat,:) = nanmean(p.ptemp,2);
            gas1_mean(iday,ilat,:) = nanmean(p.gas_1,2);
            gas2_mean(iday,ilat,:) = nanmean(p.gas_2,2);
            gas3_mean(iday,ilat,:) = nanmean(p.gas_3,2);
            gas4_mean(iday,ilat,:) = nanmean(p.gas_4,2);
            gas5_mean(iday,ilat,:) = nanmean(p.gas_5,2);
            gas6_mean(iday,ilat,:) = nanmean(p.gas_6,2);
            gas9_mean(iday,ilat,:) = nanmean(p.gas_9,2);
            gas12_mean(iday,ilat,:) = nanmean(p.gas_12,2);
            spres_mean(iday,ilat) = nanmean(p.spres);
            nlevs_mean(iday,ilat) = nanmean(p.nlevs);
            iudef4_mean(iday,ilat) = nanmean(p.iudef(4,:));
            satzen_mean(iday,ilat) = nanmean(p.satzen);
            plevs_mean(iday,ilat,:) = nanmean(p.plevs,2);
            mmwater_mean(iday,ilat) = nanmean(binwater);

            % average and store cloud variables
            ind1 = find((p.ctype == 201 & p.cfrac > 0 & p.cngwat > 0 & p.cprtop > 0) |  ...
                        (p.ctype2 == 201 & p.cfrac2 > 0 & p.cngwat2 ...
                         > 0 & p.cprtop2 > 0));
            cfrac_mean(iday,ilat) = nanmean(p.cfrac(ind1));
            cprtop_mean(iday,ilat) = nanmean(p.cprtop(ind1));
            cprbot_mean(iday,ilat) = nanmean(p.cprbot(ind1));
            cngwat_mean(iday,ilat) = nanmean(p.cngwat(ind1));
            cpsize_mean(iday,ilat) = nanmean(p.cpsize(ind1));

            ind2 = find((p.ctype == 101 & p.cfrac > 0 & p.cngwat > 0 & p.cprtop > 0) |  ...
                        (p.ctype2 == 101 & p.cfrac2 > 0 & p.cngwat2 ...
                         > 0 & p.cprtop2 > 0));
            cfrac2_mean(iday,ilat) = nanmean(p.cfrac2(ind2));
            cprtop2_mean(iday,ilat) = nanmean(p.cprtop2(ind2));
            cprbot2_mean(iday,ilat) = nanmean(p.cprbot2(ind2));
            cngwat2_mean(iday,ilat) = nanmean(p.cngwat2(ind2));
            cpsize2_mean(iday,ilat) = nanmean(p.cpsize2(ind2));

            cfrac12_mean(iday,ilat) = nanmean(p.cfrac12);
            
        end  % end loop over latitudes
            iday = iday + 1
    end % if a.bytes > 1000000
end  % giday
eval_str = ['save ~/WorkingFiles/data/stats/airs/rtp_airibrad_era_rad_'  int2str(year) ...
            '_random_kl' sDescriptor ' robs rcldy rclr bias_std *_mean count latbins'];
eval(eval_str);
