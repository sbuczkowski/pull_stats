function pull_stats_random_template(year);

addpath /asl/matlib/h4tools
addpath /asl/rtp_prod/airs/utils
addpath ~/git/rtp_prod2/util
addpath /home/sergio/MATLABCODE/PLOTTER  % equal_area_spherical_bands

% record run start datetime in output stats file for tracking
trace.RunDate = datetime('now','TimeZone','local','Format', ...
                         'd-MMM-y HH:mm:ss Z');

%**********************************************************************
% changes to input/output locations, filenames and binning are made
% here

nBins = 20; % number of lat bins to separate obs into (generates
            % 2*nBins actual bins via a (2*nBins)+1 element array
            % of bin edges

% source dir for input rtp files
sRtpFileDir = '/asl/data/rtp_airibrad_v5';
% regexp to find input rtp file names (can vary era/ecmwf,
% airibrad/airxbcal and, random/clear/site/dcc. Currently one such
% rtp file per day
sRtpFileGlob = 'era_airibrad*_random.rtp';
% type: random/clear/dcc/site
asType = 'random';
% output directory
sOutputDir = '/home/sbuczko1/testoutput/2015/airs/random/';
%**********************************************************************

% needed for calnum_to_data bad data indexing (copied from
% elsewhere. True provenance unknown)
cstr =[ 'bits1-4=NEdT[0.08 0.12 0.15 0.20 0.25 0.30 0.35 0.4 0.5 0.6 0.7' ...
        ' 0.8 1.0 2.0 4.0 nan]; bit5=Aside[0=off,1=on]; bit6=Bside[0=off,1=on];' ...
        ' bits7-8=calflag&calchansummary[0=OK, 1=DCR, 2=moon, 3=other]' ];

% locate input files
sBaseDir = fullfile(sRtpFileDir, int2str(year), asType);
dayfiles = dir(fullfile(sBaseDir, sRtpFileGlob));
fprintf(1,'>>> numfiles = %d\n', length(dayfiles));

% calculate latitude bins
latbins = equal_area_spherical_bands(nBins); % gives 2N+1 element array of lat bin boundaries
nlatbins = length(latbins);

% loop over days, computing stats individually for each day
iday = 1;
% for giday = 1:50:length(dayfiles)
for giday = 1:length(dayfiles)
    fprintf(1, '>>> year = %d  :: giday = %d\n', year, giday);
    a = dir(fullfile(sBaseDir,dayfiles(giday).name));
    if a.bytes > 100000
        [h,ha,p,pa] = rtpread(fullfile(sBaseDir,dayfiles(giday).name));
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
        
        %**************************************************
        % filter robs 
        k = find(p.iudef(4,:) == 68 & p.landfrac == 0); % night +
                                                        % ocean
        sDescriptor = 'night_ocean';
        %**************************************************

        % subset according to robs filter
        pp = rtp_sub_prof(p, k);

        % Initialize counts
        [nedt,ab,ical] = calnum_to_data(pp.calflag,cstr);
        n = length(pp.rlat);
        count_all = ones(2378,n);
        for i=1:2378
            % Find bad channels
            k = find( pp.robs1(i,:) == -9999 | ical(i,:) ~= 0 | nedt(i,:) > 1);

            % NaN's for bad channels
            pp.robs1(i,k) = NaN;
            pp.rcalc(i,k) = NaN;
            pp.sarta_rclearcalc(i,k) = NaN;
            count_all(i,k) = 0;
        end

        % Loop over latitude bins
        for ilat = 1:nlatbins-1
            % subset based on latitude bin
            inbin = find(pp.rlat > latbins(ilat) & pp.rlat <= ...
                         latbins(ilat+1));
            p = rtp_sub_prof(pp,inbin);
            bincount = count_all(:,inbin); 
            
            % Radiance mean and std
            r  = p.robs1;
            rc = p.rcalc;
            rclear = p.sarta_rclearcalc;
            
            robs(iday,ilat,:) = nanmean(r,2);
            rcal(iday,ilat,:) = nanmean(rc,2);
            rclrcal(iday,ilat,:) = nanmean(rclear,2);
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
        end  % end loop over latitudes
            iday = iday + 1
    end % if a.bytes > 1000000
end  % giday

% Output aggregate stats as mat file 
sOutputFile = [sOutputDir '/rtp_airibrad_rad_' int2str(year) '_random_' ...
               sDescriptor];
eval_str = ['save ' sOutputFile ...
            ' robs rcal rclrcal *_mean count latbins trace'];
eval(eval_str);
