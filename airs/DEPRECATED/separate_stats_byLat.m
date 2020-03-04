% separate an existing stats files and create files by latitude bin

clear all
load /asl/data/stats/airs/clear/rtp_airicrad_era_rad_kl_16day_2018_clear_desc.mat

% latbinedges and trace are resaved without modification

base_count = count;
clear count;
base_lat_mean = lat_mean;
clear lat_mean;
base_rclr = rclr;
clear rclr;
base_spres_mean = spres_mean;
clear spres_mean;
base_dbtun_mean = dbtun_mean;
clear dbtun_mean;
base_rclrbias_std = rclrbias_std;
clear rclrbias_std;
base_stemp_mean = stemp_mean;
clear stemp_mean;
base_gas1_mean = gas1_mean;
clear gas1_mean;
base_lon_mean = lon_mean;
clear lon_mean;
base_robs = robs;
clear robs;
base_tcc_mean = tcc_mean;
clear tcc_mean;
base_gas3_mean = gas3_mean;
clear gas3_mean;
base_mmwater_mean = mmwater_mean;
clear mmwater_mean;
base_rtime_mean = rtime_mean;
clear rtime_mean;
base_iudef4_mean = iudef4_mean;
clear iudef4_mean;
base_nlevs_mean = nlevs_mean;
clear nlevs_mean;
base_satazi_mean = satazi_mean;
clear satazi_mean;
base_l1cproc_mean = l1cproc_mean;
clear l1cproc_mean;
base_plevs_mean = plevs_mean;
clear plevs_mean;
base_satzen_mean = satzen_mean;
clear satzen_mean;
base_l1csreason_mean = l1csreason_mean;
clear l1csreason_mean;
base_ptemp_mean = ptemp_mean;
clear ptemp_mean;
base_solzen_mean = solzen_mean;
clear solzen_mean;

for ilat = 2:2:38
    count = squeeze(base_count(:,ilat,:));
    rclr = squeeze(base_rclr(:,ilat,:));
    robs = squeeze(base_robs(:,ilat,:));
    rclrbias_std = squeeze(base_rclrbias_std(:,ilat,:));
    l1cproc_mean = squeeze(base_l1cproc_mean(:,ilat,:));
    l1csreason_mean = squeeze(base_l1csreason_mean(:,ilat,:));
    plevs_mean = squeeze(base_plevs_mean(:,ilat,:));
    ptemp_mean = squeeze(base_ptemp_mean(:,ilat,:));
    gas1_mean = squeeze(base_gas1_mean(:,ilat,:));
    gas3_mean = squeeze(base_gas3_mean(:,ilat,:));
    
    lat_mean = squeeze(base_lat_mean(:,ilat));
    lon_mean = squeeze(base_lon_mean(:,ilat));
    spres_mean = squeeze(base_spres_mean(:,ilat));
    dbtun_mean = squeeze(base_dbtun_mean(:,ilat));
    stemp_mean = squeeze(base_stemp_mean(:,ilat));
    tcc_mean = squeeze(base_tcc_mean(:,ilat));
    mmwater_mean = squeeze(base_mmwater_mean(:,ilat));
    rtime_mean = squeeze(base_rtime_mean(:,ilat));
    iudef4_mean = squeeze(base_iudef4_mean(:,ilat));
    nlevs_mean = squeeze(base_nlevs_mean(:,ilat));
    satazi_mean = squeeze(base_satazi_mean(:,ilat));
    satzen_mean = squeeze(base_satzen_mean(:,ilat));
    solzen_mean = squeeze(base_solzen_mean(:,ilat));

    outfile = sprintf(['/asl/data/stats/airs/clear/' ...
                       'rtp_airicrad_era_rad_kl_16day_lb%d_2018_clear_desc.mat'], ilat);;
    save(outfile, 'count','rclr','robs','rclrbias_std','l1cproc_mean','l1csreason_mean', ...
        'plevs_mean','ptemp_mean','gas1_mean','gas3_mean','lat_mean','lon_mean', ...
        'spres_mean','dbtun_mean','stemp_mean','tcc_mean','mmwater_mean', ...
        'rtime_mean','iudef4_mean','nlevs_mean','satazi_mean','satzen_mean', ...
         'solzen_mean','latbinedges','trace')

end

    

