function run_pull_stats_airxbcal_site(year)
set_process_dirs;
addpath(genpath(rtp_sw_dir));

% grab the slurm array index for this process
slurmindex = str2num(getenv('SLURM_ARRAY_TASK_ID'));

siteID = 1+slurmindex;
disp(year)
disp(siteID)
pull_stats_airxbcal_site(year,siteID);
