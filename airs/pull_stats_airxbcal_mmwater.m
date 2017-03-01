function pull_stats_airxbcal_mmwater(year, filter)

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
addpath /asl/packages/rtp_prod2/util
addpath /home/sergio/MATLABCODE/PLOTTER  % equal_area_spherical_bands
addpath /asl/matlib/rtptools  % mmwater_rtp.m

trace.RunDate = datetime('now','TimeZone','local','Format', ...
                         'd-MMM-y HH:mm:ss Z');
klayers_exec = ['/asl/packages/klayersV205/BinV201/' ...
                'klayers_airs_wetwater'];

[sID, sTempPath] = genscratchpath();

cstr =[ 'bits1-4=NEdT[0.08 0.12 0.15 0.20 0.25 0.30 0.35 0.4 0.5 0.6 0.7' ...
  ' 0.8 1.0 2.0 4.0 nan]; bit5=Aside[0=off,1=on]; bit6=Bside[0=off,1=on];' ...
  ' bits7-8=calflag&calchansummary[0=OK, 1=DCR, 2=moon, 3=other]' ];

basedir = ['/asl/data/rtp_airxbcal_v5/' int2str(year) '/clear'];
dayfiles = dir(fullfile(basedir, 'ecmwf_airxbcal*.rtp'));
nfiles = length(dayfiles);
fprintf(1,'>>> numfiles = %d\n', nfiles);
if nfiles == 0
    fprintf(2, '>>> Error: No input files found. Exiting\n');
    return;
end

% calculate bins
% bin mmwater in 5mm bins out to 100mm
binedges = [0:5:100];   
nbins = length(binedges);

iday = 1;
% for giday = 1:50:length(dayfiles)
for giday = 1:length(dayfiles)
   fprintf(1, '>>> year = %d  :: giday = %d\n', year, giday);
   a = dir(fullfile(basedir,dayfiles(giday).name));
   a.bytes;
   if a.bytes > 100000
       % read in rtp file
       [h,ha,p,pa] = rtpread(fullfile(basedir,dayfiles(giday).name));
       
       % pre-filter rtp data to minimize load going into klayers
       switch filter
         case 1         % satzen +-10
           k = find(abs(p.satzen) < 10);
           sDescriptor='_satzen-10';
% $$$          case 2          % descending node (night), mid-lat
% $$$            k = find(p.iudef(4,:) == 68 & abs(p.rlat) > 30 & abs(p.rlat) <= 50);
% $$$            sDescriptor='_desc_midlat';
       end
       
       pp = rtp_sub_prof(p, k);

       % write out subset data for input to klayers
       fn_rtp1 = fullfile(sTempPath, ['airs_' sID '_1.rtp']);
       rtpwrite(fn_rtp1, h,ha,pp,pa);

       % run klayers on the rtp data (Sergio is asking for this to
       % convert levels to layers for his processing?)
       fprintf(1, '>>> running klayers... ');
       fn_rtp2 = fullfile(sTempPath, ['airs_' sID '_2.rtp']);
       klayers_run = [klayers_exec ' fin=' fn_rtp1 ' fout=' fn_rtp2 ...
                      ' > ' sTempPath '/kout.txt'];
       unix(klayers_run);
       fprintf(1, 'Done\n');

       % Read klayers output into local rtp variables
       [h,ha,pp,pa] = rtpread(fn_rtp2);

       % get column water
       mmwater = mmwater_rtp(h, pp);

       f = h.vchan;  % AIRS proper frequencies

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
         pp.rcalc(i,k) = NaN;
         count_all(i,k) = 0;
      end

      % Loop over latitude bins
      for ibin = 1:nbins-1
           % subset based on latitude bin
           inbin = find(mmwater > binedges(ibin) & mmwater <= ...
                        binedges(ibin+1));
           p = rtp_sub_prof(pp,inbin);
           bincount = count_all(:,inbin);
           
           % Radiance mean and std
           r  = p.robs1;
           rc = p.rcalc;

           % B(T) bias mean and std
           bto = real(rad2bt(f,r));
           btc = real(rad2bt(f,rc));
           btobs(iday,ibin,:) = nanmean(bto,2);
           btcal(iday,ibin,:) = nanmean(btc,2);
           bias_std(iday,ibin,:) = nanstd(bto-btc,0,2);
           lat_mean(iday,ibin) = nanmean(p.rlat);
           lon_mean(iday,ibin) = nanmean(p.rlon);
           solzen_mean(iday,ibin) = nanmean(p.solzen);
           rtime_mean(iday,ibin)  = nanmean(p.rtime);
           count(iday,ibin,:) = sum(bincount,2)';
           stemp_mean(iday,ibin) = nanmean(p.stemp);
           ptemp_mean(iday,ibin,:) = nanmean(p.ptemp,2);
           gas1_mean(iday,ibin,:) = nanmean(p.gas_1,2);
           gas3_mean(iday,ibin,:) = nanmean(p.gas_3,2);
           spres_mean(iday,ibin) = nanmean(p.spres);
           nlevs_mean(iday,ibin) = nanmean(p.nlevs);
           iudef4_mean(iday,ibin) = nanmean(p.iudef(4,:));
           satzen_mean(iday,ibin) = nanmean(p.satzen);
           plevs_mean(iday,ibin,:) = nanmean(p.plevs,2);
       end  % end loop over latitudes

      iday = iday + 1
   end % if a.bytes > 1000000
end  % giday
eval_str = ['save ~/testoutput/2015/airs/clear/rtp_airxbcal_ecmwf_mmwater_'  int2str(year) ...
            '_clear' sDescriptor ' btobs btcal bias_std *_mean count ' ...
                    'trace binedges'];
eval(eval_str);
