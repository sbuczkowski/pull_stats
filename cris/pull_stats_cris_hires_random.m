function pull_stats_cris_hires_random(year, filter, cfg);

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

rtpdir = '/asl/rtp/rtp_cris_ccast_hires/';
if bCheckConfig & isfield(cfg, 'rtpdir')
    rtpdir = cfg.rtpdir;
end

statsdir = '/asl/data/stats/cris';
if bCheckConfig & isfield(cfg, 'statsdir')
    statsdir = cfg.statsdir;
end

basedir = fullfile(rtpdir, 'random', int2str(year));
dayfiles = dir(fullfile(basedir, 'cris2_ecmwf_csarta_random_d*.rtp'));
fprintf(1,'>>> numfiles = %d\n', length(dayfiles));
ndays = length(dayfiles);

% Get proper frequencies for these data
[n1,n2,n3,userLW,userMW,userSW, ichan] = cris_hires_chans();
f = cris_vchan(2, userLW, userMW, userSW);
nchans = length(f);

% calculate latitude bins
nbins=20; % gives 2N+1 element array of lat bin boundaries
latbin_edges = equal_area_spherical_bands(nbins);
nlatbins = length(latbin_edges) - 1;

nlevs = 101;  % klayers output

% allocate final accumulator arrays
robs = zeros(ndays, nlatbins, nchans);
rclr = zeros(ndays, nlatbins, nchans);
rcld = zeros(ndays, nlatbins, nchans);
rcldbias_std = zeros(ndays, nlatbins, nchans);
rclrbias_std = zeros(ndays, nlatbins, nchans);

lat_mean = zeros(ndays, nlatbins);
lon_mean = zeros(ndays, nlatbins);
solzen_mean = zeros(ndays, nlatbins);
rtime_mean = zeros(ndays, nlatbins); 
count = zeros(ndays, nlatbins, nchans);
tcc_mean = zeros(ndays, nlatbins);
stemp_mean = zeros(ndays, nlatbins);
ptemp_mean = zeros(ndays, nlatbins, nlevs);
gas1_mean = zeros(ndays, nlatbins, nlevs);
gas3_mean = zeros(ndays, nlatbins, nlevs);
spres_mean = zeros(ndays, nlatbins);
nlevs_mean = zeros(ndays, nlatbins);
iudef4_mean = zeros(ndays, nlatbins);
mmwater_mean = zeros(ndays, nlatbins);
satzen_mean = zeros(ndays, nlatbins);
plevs_mean = zeros(ndays, nlatbins, nlevs);

iday = 1;
% for giday = 1:50:length(dayfiles)
for giday = 1:length(dayfiles)
    fprintf(1, '>>> year = %d  :: giday = %d\n', year, giday);
    a = dir(fullfile(basedir,dayfiles(giday).name));
    if a.bytes < 100000
        fprintf(2, '**>> ERROR: short input rtp file %s\n', dayfiles(giday).name); 
        continue;
    end
    
    [h,ha,prof,pa] = rtpread(fullfile(basedir,dayfiles(giday).name));
    f = h.vchan;  % CrIS proper frequencies

    % nadir sub-select (**for 06132018 testing purposes**)
    nadir_inds = find(prof.xtrack == 15 | prof.xtrack == 16);
    p = rtp_sub_prof(prof, nadir_inds);
    fprintf(1, 'Pre/post nadir selection obs count: %d / %d\n', ...
            length(prof.robs1), length(p.robs1));
    clear prof
    trace.PROCESSING_NOTE = ['Full-swath subset to nadir: xtrack ' ...
                        '= 15,16'];
    
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
    clear p;

    if bRunKlayers
        % run klayers on the rtp data to convert levels -> layers
        % save calcs as the re-run of klayers wipes them out
        tmp_rclr = pp.rclr;
        tmp_rcld = pp.rcld;
        tmp_tcc = pp.tcc;
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
        % restore rclr
        pp.rclr = tmp_rclr;
        pp.rcld = tmp_rcld;
        pp.tcc = tmp_tcc;
        clear tmp_rclr tmp_rcld tmp_tcc;
        
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
    end % end if bRunKlayers
        
        % initialize counts and look for bad channels (what should
        % the iasi bad channel test look like?)
        [nchans, nobs] = size(pp.robs1);
        nfovs = 9;
        count_all = int8(ones(nchans, nobs, nfovs));
        
        % loop over latitude bins
        for ilat = 1:nlatbins
            % subset based on latitude bin
            inbin = find(pp.rlat > latbin_edges(ilat) & pp.rlat <= ...
                         latbin_edges(ilat+1));
            p = rtp_sub_prof(pp,inbin);

            for z = 1:9  % loop over FOVs to further sub-select
                ifov = find(p.ifov == z);
                p2 = rtp_sub_prof(p, ifov);
                
                bincount = count_all(:,inbin,z); 
                binwater = mmwater(inbin);

                % Remove 'clouds' that have only partial defining
                % characteristics (effects 0.2% of obs in test rtp data)
                % cfrac=cngwat=0 but cpsize/cprtop/cprbot ~= 0
                a = find(p2.cfrac==0 & p2.cngwat==0 & p2.cprtop>0);
                p2.ctype(a)=-1;
                a = find(p2.cfrac2==0 & p2.cngwat2==0 & p2.cprtop2>0);
                p2.ctype2(a)=-1;

                % replace 'bad data' placeholders in cprtop/cprbot (-9999
                % -> NaN, esesntially except for 'clouds' from above)
                a = find(p2.ctype==-1);
                p2.cprtop(a)=NaN; p2.cprbot(a)=NaN;
                a = find(p2.ctype2==-1);
                p2.cprtop2(a)=NaN; p2.cprbot2(a)=NaN;

                % massage situations where we have two of same cloud type
                % in a single obs i.e. p2.ctype=p2.ctype2 (average to a
                % single cloud)
                dblwat = find(p2.ctype==101 & p2.ctype2==101);
                p2.cfrac2(dblwat) = p2.cfrac(dblwat) + p2.cfrac2(dblwat);
                p2.cpsize2(dblwat) = (p2.cpsize(dblwat) + ...
                                      p2.cpsize2(dblwat))/2;
                p2.cprtop2(dblwat) = (p2.cprtop(dblwat) + ...
                                      p2.cprtop2(dblwat))/2;
                p2.cprbot2(dblwat) = (p2.cprbot(dblwat) + ...
                                      p2.cprbot2(dblwat))/2;
                p2.cngwat2(dblwat) = (p2.cfrac(dblwat).*p2.cngwat(dblwat) + ...
                                      p2.cfrac2(dblwat).*p2.cngwat2(dblwat))./ ...
                    p2.cfrac2(dblwat);
                p2.cfrac(dblwat)=0; p2.cpsize(dblwat)=0;
                p2.cngwat(dblwat)= 0; p2.cprtop(dblwat)=NaN;
                p2.cprbot(dblwat)=NaN; p2.ctype(dblwat)=-1;
                
                dblice = find(p2.ctype==201 & p2.ctype2==201);
                p2.cfrac(dblice) = p2.cfrac(dblice) + p2.cfrac2(dblice);
                p2.cpsize(dblice) = (p2.cpsize(dblice) + ...
                                     p2.cpsize2(dblice))/2;
                p2.cprtop(dblice) = (p2.cprtop(dblice) + ...
                                     p2.cprtop2(dblice))/2;
                p2.cprbot(dblice) = (p2.cprbot(dblice) + ...
                                     p2.cprbot2(dblice))/2;
                p2.cngwat(dblice) = (p2.cfrac(dblice).*p2.cngwat(dblice) + ...
                                     p2.cfrac2(dblice).*p2.cngwat2(dblice))./ ...
                    p2.cfrac2(dblice);
                p2.cfrac2(dblice)=0; p2.cpsize2(dblice)=0;
                p2.cngwat2(dblice)= 0; p2.cprtop2(dblice)=NaN;
                p2.cprbot2(dblice)=NaN; p2.ctype2(dblice)=-1;
                
                % at this point cloud fields should now be set so that
                % ctype/cpsize/cngwat and similar contain only ice clouds
                % (or none) and ctype2/cpsize2/cngwat2 and similar
                % contain only water clouds (or none)
                trace.NOTE_ON_CLOUD_VARS = ['ctype_mean and similar are average ' ...
                                    'ICE cloud parameters, ctype2_mean and similar ' ...
                                    'are average WATER cloud parameters'];

                % pull out cloud parameters for averaging
                % build logical array for each cloud type and pull out
                % only values associated with a cloud for averaging
                lIce = (p2.ctype==201 & p2.cfrac>0 & p2.cngwat>0 & p2.cpsize>0 ...
                        & p2.cprbot>0 & p2.cprtop>0);

                ice_ind = find(lIce);
                count_ice(iday, ilat) = length(ice_ind);
                ctype_mean(iday, ilat) = nanmean(p2.ctype(ice_ind));
                cngwat_mean(iday, ilat) = nanmean(p2.cngwat(ice_ind));
                cpsize_mean(iday, ilat) = nanmean(p2.cpsize(ice_ind));
                cprbot_mean(iday, ilat) = nanmean(p2.cprbot(ice_ind));
                cprtop_mean(iday, ilat) = nanmean(p2.cprtop(ice_ind));
                
                lWat = (p2.ctype2==101 & p2.cfrac2>0 & p2.cngwat2>0 & p2.cpsize2>0 ...
                        & p2.cprbot2>0 & p2.cprtop2>0);

                wat_ind = find(lWat);
                count_water(iday, ilat) = length(wat_ind); 
                ctype2_mean(iday, ilat) = nanmean(p2.ctype2(wat_ind));
                cngwat2_mean(iday, ilat) = nanmean(p2.cngwat2(wat_ind));
                cpsize2_mean(iday, ilat) = nanmean(p2.cpsize2(wat_ind));
                cprbot2_mean(iday, ilat) = nanmean(p2.cprbot2(wat_ind));
                cprtop2_mean(iday, ilat) = nanmean(p2.cprtop2(wat_ind));

                % cloud fraction gets averaged over ALL obs
                cfrac_mean(iday, ilat) = nanmean(p2.cfrac);
                cfrac2_mean(iday, ilat) = nanmean(p2.cfrac2);
                cfrac12_mean(iday, ilat) = nanmean(p2.cfrac12);

                % Radiance mean and std
                
                r  = p2.robs1;
                clr_calc = p2.rclr;
                cldy_calc = p2.rcld;
                
                % leave as sinc for test
                % Convert r to rham
                r = box_to_ham(r);  % assumes r in freq order!!  Needed
                                    % for lowres
                
                robs(iday,ilat,:,z) = nanmean(r,2);
                rclr(iday,ilat,:,z) = nanmean(clr_calc,2);
                rcld(iday,ilat,:,z) = nanmean(cldy_calc,2);
                rclrbias_std(iday, ilat,:) = nanstd(r-clr_calc,0,2);
                rcldbias_std(iday, ilat,:) = nanstd(r-cldy_calc,0,2);
                
                lat_mean(iday,ilat,z) = nanmean(p2.rlat);
                lon_mean(iday,ilat,z) = nanmean(p2.rlon);
                solzen_mean(iday,ilat,z) = nanmean(p2.solzen);
                rtime_mean(iday,ilat,z)  = nanmean(p2.rtime);
                count(iday,ilat,z) = sum(bincount(1,:))';
                stemp_mean(iday,ilat,z) = nanmean(p2.stemp);
                iudef4_mean(iday,ilat,z) = nanmean(p2.iudef(4,:));
                ptemp_mean(iday,ilat,:,z) = nanmean(p.ptemp,2);
                gas1_mean(iday,ilat,:,z) = nanmean(p.gas_1,2);
                gas3_mean(iday,ilat,:,z) = nanmean(p.gas_3,2);
                spres_mean(iday,ilat,z) = nanmean(p.spres);
                nlevs_mean(iday,ilat,z) = nanmean(p.nlevs);
                satzen_mean(iday,ilat,z) = nanmean(p.satzen);
                plevs_mean(iday,ilat,:,z) = nanmean(p.plevs,2);
                mmwater_mean(iday,ilat) = nanmean(binwater);

            end  % ifov (z)
        end  % end loop over ilat
            
            iday = iday + 1
end  % giday
outfile = fullfile(statsdir, sprintf('rtp_cris2_hires_ecmwf_rad_kl_%s_random_nadir_%s', ...
                                     int2str(year), sDescriptor));
eval_str = ['save ' outfile ' robs rclr rcld *_std *_mean count trace'];
eval(eval_str);
