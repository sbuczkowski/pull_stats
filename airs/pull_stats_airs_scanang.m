function pull_stats_airs_scanang(year, filter, cfg);

% Built around AIRICRAD AIRS data product

[sID, sTempPath] = genscratchpath();
if isfield(cfg, 'stemppath')
    sTempPath = cfg.stemppath;
    rng('shuffle');
    sID = sprintf('%03d', randi(999));
end
fprintf(1, '>> TempPath: %s\n', sTempPath)

% record run start datetime in output stats file for tracking
trace.RunDate = datetime('now','TimeZone','local','Format', ...
                         'd-MMM-y HH:mm:ss Z');
assert(isfield(cfg, 'reason'))
trace.Reason = cfg.reason;

bRunKlayers = true;
klayers_exec = ['/asl/packages/klayersV205/BinV201/' ...
                'klayers_airs_wetwater'];
if isfield(cfg, 'klayers') & cfg.klayers == false
    bRunKlayers = false;
end
trace.klayers = bRunKlayers;
trace.klayers_exec = klayers_exec;

trace.droplayers = false;
if isfield(cfg, 'droplayers') & cfg.droplayers == true
    trace.droplayers = true;
end

model = 'era';
if isfield(cfg, 'model')
    model = cfg.model;
end

instName = cfg.instname; 

rtpsrcdir = sprintf('/asl/rtp/rtp_airicrad_v6');
if isfield(cfg, 'rtpsrcdir')
    rtpsrcdir = cfg.rtpsrcdir;
end

descriptor = 'random';
if isfield(cfg, 'descriptor')
    descriptor = cfg.descriptor;
end

statsdir = sprintf('/asl/data/stats/%s/%s', instName,descriptor);
if isfield(cfg, 'statsdir')
    statsdir = cfg.statsdir;
end

rta = 'sarta';
if isfield(cfg, 'rta')
    rta = cfg.rta;
end

bIncludeCalcs = true;
if isfield(cfg, 'includecalcs') & cfg.includecalcs == 0
    bIncludeCalcs = false;
end

basedir = fullfile(rtpsrcdir, sprintf('%s', descriptor), ...
                   int2str(year));
fprintf(1, '>> Basedir: %s\n', basedir);
% era_airicrad_day241_clear.rtp
namefilter = sprintf('%s_airicrad_d*_%s.rtp', model, descriptor);
fprintf(1, '>> namefilter: %s\n', namefilter)
dayfiles = dir(fullfile(basedir, namefilter));


ndays = length(dayfiles);
fprintf(1,'>>> numfiles = %d\n', ndays);

% Get proper frequencies for these data
nchans = 2645;  %AIRICRAD/L1C channel space

nfors = 90;

nlevs = 101;  % klayers output

% allocate final accumulator arrays for common variables
robs_mean = nan(ndays, nfors, nchans);

if bIncludeCalcs
    rclr_mean = nan(ndays, nfors, nchans);
    rclrbias_std = nan(ndays, nfors, nchans);
end

lat_mean = nan(ndays, nfors);
lon_mean = nan(ndays, nfors);
solzen_mean = nan(ndays, nfors);
rtime_mean = nan(ndays, nfors); 
count = nan(ndays, nfors, nchans);
tcc_mean = nan(ndays, nfors);
stemp_mean = nan(ndays, nfors);
ptemp_mean = nan(ndays, nfors, nlevs);
gas1_mean = nan(ndays, nfors, nlevs);
gas3_mean = nan(ndays, nfors, nlevs);
spres_mean = nan(ndays, nfors);
nlevs_mean = nan(ndays, nfors);
iudef4_mean = nan(ndays, nfors);
mmwater_mean = nan(ndays, nfors);
scanang_mean = nan(ndays, nfors);
satzen_mean = nan(ndays, nfors);
satazi_mean = nan(ndays, nfors);
plevs_mean = nan(ndays, nfors, nlevs);
dbtun_mean = nan(ndays, nfors);

% AIRS-specific variables (for all subsets)
l1cproc_mean = nan(ndays, nfors, nchans);
l1csreason_mean = nan(ndays, nfors, nchans);


if strcmp(descriptor, 'clear')
    % add in clear subset specific variables
    dbtun = nan(ndays, nfors);
else

    % add in random/site/dcc specific variables (basically stuff that
    % tags along with sarta_cloudy runs)
    if bIncludeCalcs
        rcld_mean = nan(ndays, nfors, nchans);
        rcldbias_std = nan(ndays, nfors, nchans);
    end

    count_ice = nan(ndays, nfors);
    ctype_mean = nan(ndays, nfors);
    cngwat_mean = nan(ndays, nfors);
    cpsize_mean = nan(ndays, nfors);
    cprbot_mean = nan(ndays, nfors);
    cprtop_mean = nan(ndays, nfors);

    count_water = nan(ndays, nfors);
    ctype2_mean = nan(ndays, nfors);
    cngwat2_mean = nan(ndays, nfors);
    cpsize2_mean = nan(ndays, nfors);
    cprbot2_mean = nan(ndays, nfors);
    cprtop2_mean = nan(ndays, nfors);

    cfrac_mean = nan(ndays, nfors);
    cfrac2_mean = nan(ndays, nfors);
    cfrac12_mean = nan(ndays, nfors);

end

iday = 1;
for giday = 1:ndays
    fprintf(1, '>>> year = %d  :: giday = %d\n', year, giday);
    a = dir(fullfile(basedir,dayfiles(giday).name));
    if a.bytes < 100000
        continue
    end
    day_str = sprintf('giday = %d, iday = %d',giday, iday);
    
    [h,ha,p_rtp,pa] = rtpread(fullfile(basedir,dayfiles(giday).name));
    

    ascobs = find(p_rtp.iudef(4,:) == 65); % day
    descobs = find(p_rtp.iudef(4,:) == 68); % night
    landobs = find(p_rtp.landfrac == 1);
    oceanobs = find(p_rtp.landfrac == 0);

    switch filter
      case 1 % descending node (night), global
        k = descobs;
        sDescriptor='desc';
      case 2 % descending node (night), ocean
        k = intersect(descobs, oceanobs);
        sDescriptor='desc_ocean';
      case 3 % descending node (night), land
        k = intersect(descobs, landobs);
        sDescriptor='desc_land';
      case 4 % ascending node (day)
        k = ascobs;
        sDescriptor='asc';
      case 5 % ascending node (day), ocean
        k = intersect(ascobs, oceanobs);
        sDescriptor='asc_ocean';
      case 6 % ascending node (day), land
        k = intersect(ascobs, landobs);
        sDescriptor='asc_land';
    end

    if length(k) == 0
        continue;
    end
    p_filt = rtp_sub_prof(p_rtp, k);
    clear p_rtp
    % run klayers on the rtp data to convert levels -> layers
    % save calcs as the re-run of klayers wipes them out
    if bIncludeCalcs
        tmp_rclr = p_filt.rclr;
        if ~strcmp(descriptor, 'clear')
            tmp_rcld = p_filt.rcld;
        end
    end

    tmp_tcc = p_filt.tcc;
    tmp_l1cproc = p_filt.l1cproc;
    tmp_l1csreason = p_filt.l1csreason;

    if strcmp(descriptor, 'clear')
        tmp_dbtun = p_filt.dbtun;
    end

    fn_rtp1 = fullfile(sTempPath, ['airs_' sID '_1.rtp']);
    fn_rtp2 = fullfile(sTempPath, ['airs_' sID '_2.rtp']);
    fprintf(1, '>>> running klayers : writing output to %s...', fn_rtp2);
    rtpwrite(fn_rtp1, h,ha,p_filt,pa);
    clear p_filt;
    klayers_run = [klayers_exec ' fin=' fn_rtp1 ' fout=' fn_rtp2 ...
                   ' > ' sTempPath '/kout.txt'];
    unix(klayers_run);
    fprintf(1, 'Done\n');

    % Read klayers output into local rtp variables
    [h,ha,p_filt,pa] = rtpread(fn_rtp2);
    % restore rclr
    if bIncludeCalcs
        p_filt.rclr = tmp_rclr;
        clear tmp_rclr
        if ~strcmp(descriptor, 'clear')
            p_filt.rcld = tmp_rcld;
            clear tmp_rcld
        end
    end
    p_filt.tcc = tmp_tcc;
    p_filt.l1cproc = tmp_l1cproc;
    p_filt.l1csreason = tmp_l1csreason;
    clear tmp_tcc tmp_l1cproc tmp_l1csreason
    if strcmp(descriptor, 'clear')
        p_filt.dbtun = tmp_dbtun;
        clear tmp_dbtun
    end
    
    % get column water
    mmwater = mmwater_rtp(h, p_filt);

    % Check for obs with layer profiles that go lower than
    % topography. Need to check nlevs and NaN out any layers
    % at or below this level

    % ** Any layers-sensitive variables added in averaging code below must
    % ** be checked here first.
    for i=1:length(p_filt.nlevs)
        badlayers = p_filt.nlevs(i) : 101;
        p_filt.plevs(badlayers, i) = NaN;
        p_filt.gas_1(badlayers, i) = NaN;
        p_filt.gas_3(badlayers, i) = NaN;
        p_filt.ptemp(badlayers, i) = NaN;
    end

    % initialize counts and look for bad channels (what should
    % the iasi bad channel test look like?)
    [nchans, nobs] = size(p_filt.robs1);
    
    % loop over latitude bins
    for ifor = 1:nfors
        % subset based on latitude bin
        inbin = find(p_filt.xtrack == ifor)
        p_inbin = rtp_sub_prof(p_filt,inbin);

        daylat_str = sprintf('%s, ifor = %d', day_str, ifor);
        
        
        count_inbin = ones(length(inbin), nchans);
        % QA/QC checks for bad chans, etc go here and set
        % elements of count_infov to zero

        binwater = mmwater(inbin);        

        if ~strcmp(descriptor, 'clear')
            % Remove 'clouds' that have only partial defining
            % characteristics (effects 0.2% of obs in test rtp data)
            % cfrac=cngwat=0 but cpsize/cprtop/cprbot ~= 0
            a = find(p_inbin.cfrac==0 & p_inbin.cngwat==0 & p_inbin.cprtop>0);
            p_inbin.ctype(a)=-1;
            a = find(p_inbin.cfrac2==0 & p_inbin.cngwat2==0 & p_inbin.cprtop>0);
            p_inbin.ctype2(a)=-1;

            % replace 'bad data' placeholders in cprtop/cprbot (-9999
            % -> NaN, esesntially except for 'clouds' from above)
            a = find(p_inbin.ctype==-1);
            p_inbin.cprtop(a)=NaN; p_inbin.cprbot(a)=NaN;
            a = find(p_inbin.ctype2==-1);
            p_inbin.cprtop(a)=NaN; p_inbin.cprbot2(a)=NaN;

            % massage situations where we have two of same cloud type
            % in a single obs i.e. p_inbin.ctype=p_inbin.ctype2 (average to a
            % single cloud)
            dblwat = find(p_inbin.ctype==101 & p_inbin.ctype2==101);
            p_inbin.cfrac2(dblwat) = p_inbin.cfrac(dblwat) + p_inbin.cfrac2(dblwat);
            p_inbin.cpsize2(dblwat) = (p_inbin.cpsize(dblwat) + ...
                                       p_inbin.cpsize2(dblwat))/2;
            p_inbin.cprtop(dblwat) = (p_inbin.cprtop(dblwat) + ...
                                      p_inbin.cprtop(dblwat))/2;
            p_inbin.cprbot2(dblwat) = (p_inbin.cprbot(dblwat) + ...
                                       p_inbin.cprbot2(dblwat))/2;
            p_inbin.cngwat2(dblwat) = (p_inbin.cfrac(dblwat).*p_inbin.cngwat(dblwat) + ...
                                       p_inbin.cfrac2(dblwat).*p_inbin.cngwat2(dblwat))./ ...
                p_inbin.cfrac2(dblwat);
            p_inbin.cfrac(dblwat)=0; p_inbin.cpsize(dblwat)=0;
            p_inbin.cngwat(dblwat)= 0; p_inbin.cprtop(dblwat)=NaN;
            p_inbin.cprbot(dblwat)=NaN; p_inbin.ctype(dblwat)=-1;
            
            dblice = find(p_inbin.ctype==201 & p_inbin.ctype2==201);
            p_inbin.cfrac(dblice) = p_inbin.cfrac(dblice) + p_inbin.cfrac2(dblice);
            p_inbin.cpsize(dblice) = (p_inbin.cpsize(dblice) + ...
                                      p_inbin.cpsize2(dblice))/2;
            p_inbin.cprtop(dblice) = (p_inbin.cprtop(dblice) + ...
                                      p_inbin.cprtop(dblice))/2;
            p_inbin.cprbot(dblice) = (p_inbin.cprbot(dblice) + ...
                                      p_inbin.cprbot2(dblice))/2;
            p_inbin.cngwat(dblice) = (p_inbin.cfrac(dblice).*p_inbin.cngwat(dblice) + ...
                                      p_inbin.cfrac2(dblice).*p_inbin.cngwat2(dblice))./ ...
                p_inbin.cfrac2(dblice);
            p_inbin.cfrac2(dblice)=0; p_inbin.cpsize2(dblice)=0;
            p_inbin.cngwat2(dblice)= 0; p_inbin.cprtop(dblice)=NaN;
            p_inbin.cprbot2(dblice)=NaN; p_inbin.ctype2(dblice)=-1;
            
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
            lIce = (p_inbin.ctype==201 & p_inbin.cfrac>0 & p_inbin.cngwat>0 & p_inbin.cpsize>0 ...
                    & p_inbin.cprbot>0 & p_inbin.cprtop>0);

            ice_ind = find(lIce);
            count_ice(iday, ifor) = length(ice_ind);
            ctype_mean(iday, ifor) = nanmean(p_inbin.ctype(ice_ind));
            cngwat_mean(iday, ifor) = nanmean(p_inbin.cngwat(ice_ind));
            cpsize_mean(iday, ifor) = nanmean(p_inbin.cpsize(ice_ind));
            cprbot_mean(iday, ifor) = nanmean(p_inbin.cprbot(ice_ind));
            cprtop_mean(iday, ifor) = nanmean(p_inbin.cprtop(ice_ind));
            
            lWat = (p_inbin.ctype2==101 & p_inbin.cfrac2>0 & p_inbin.cngwat2>0 & p_inbin.cpsize2>0 ...
                    & p_inbin.cprbot2>0 & p_inbin.cprtop>0);

            wat_ind = find(lWat);
            count_water(iday, ifor) = length(wat_ind); 
            ctype2_mean(iday, ifor) = nanmean(p_inbin.ctype2(wat_ind));
            cngwat2_mean(iday, ifor) = nanmean(p_inbin.cngwat2(wat_ind));
            cpsize2_mean(iday, ifor) = nanmean(p_inbin.cpsize2(wat_ind));
            cprbot2_mean(iday, ifor) = nanmean(p_inbin.cprbot2(wat_ind));
            cprtop_mean(iday, ifor) = nanmean(p_inbin.cprtop(wat_ind));

            % cloud fraction gets averaged over ALL obs
            cfrac_mean(iday, ifor) = nanmean(p_inbin.cfrac);
            cfrac2_mean(iday, ifor) = nanmean(p_inbin.cfrac2);
            cfrac12_mean(iday, ifor) = nanmean(p_inbin.cfrac12);
        end
        
        r = p_inbin.robs1;
        robs_mean(iday,ifor,:) = nanmean(r,2);
        if bIncludeCalcs
            rclr_mean(iday,ifor,:) = nanmean(p_inbin.rclr,2);
            rclrbias_std(iday, ifor,:) = nanstd(r-p_inbin.rclr,0, ...
                                                2);
            if ~strcmp(descriptor, 'clear')
                rcld_mean(iday,ifor,:) = nanmean(p_inbin.rcld,2);
                rcldbias_std(iday, ifor,:) = nanstd(r-p_inbin.rcld,0, ...
                                                    2);
            end
        end

        if strcmp(descriptor, 'clear')
            dbtun_mean(iday,ifor) = nanmean(p_inbin.dbtun);
        end
        
        lat_mean(iday,ifor) = nanmean(p_inbin.rlat);
        lon_mean(iday,ifor) = nanmean(p_inbin.rlon);
        solzen_mean(iday,ifor) = nanmean(p_inbin.solzen);
        rtime_mean(iday,ifor)  = nanmean(p_inbin.rtime);
        count(iday,ifor,:) = sum(count_inbin);
        stemp_mean(iday,ifor) = nanmean(p_inbin.stemp);
        iudef4_mean(iday,ifor) = nanmean(p_inbin.iudef(4,:));
        ptemp_mean(iday,ifor,:) = nanmean(p_inbin.ptemp,2)';
        gas1_mean(iday,ifor,:) = nanmean(p_inbin.gas_1,2)';
        gas3_mean(iday,ifor,:) = nanmean(p_inbin.gas_3,2)';
        spres_mean(iday,ifor) = nanmean(p_inbin.spres);
        nlevs_mean(iday,ifor) = nanmean(p_inbin.nlevs);
        satzen_mean(iday,ifor) = nanmean(p_inbin.satzen);
        satazi_mean(iday,ifor) = nanmean(p_inbin.satazi);
        plevs_mean(iday,ifor,:) = nanmean(p_inbin.plevs,2)';
        mmwater_mean(iday,ifor) = nanmean(binwater);
        scanang_mean(iday,ifor) = nanmean(p_inbin.scanang);
        l1cproc_mean(iday,ifor,:) = nanmean(p_inbin.l1cproc,2)';
        l1csreason_mean(iday,ifor,:) = nanmean(p_inbin.l1csreason,2)';
        tcc_mean(iday, ifor) = nanmean(p_inbin.tcc);
        
    end  % end loop over ifor
    
        iday = iday + 1
end  % giday
outfile = fullfile(statsdir, sprintf('rtp_airicrad_rad_kl_scanangle_%4d_%s_%s', ...
                                     year, descriptor, ...
                                     sDescriptor));

if bIncludeCalcs
    eval_str = sprintf('save %s *_std *_mean count trace dayfiles -v7.3', ...
                   outfile);
    if ~strcmp(descriptor, 'clear')
        eval_str = sprintf('save %s *_std *_mean count trace dayfiles -v7.3', ...
                           outfile);
    end
else
    eval_str = sprintf('save %s *_mean count trace dayfiles -v7.3', ...
                   outfile);
end

fprintf(1, '>> Executing save command: \n\t%s\n', eval_str)
eval(eval_str);

fprintf(1, 'Done\n')
