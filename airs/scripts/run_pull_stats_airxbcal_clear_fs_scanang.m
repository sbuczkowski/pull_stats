function run_pull_stats_airxbcal_clear_fs_scanang(filter)
set_process_dirs;
addpath(genpath(rtp_sw_dir));

% grab the slurm array index for this process
slurmindex = str2num(getenv('SLURM_ARRAY_TASK_ID'));

cfg.klayers = true;

year = 2002+slurmindex;
disp(year)
disp(filter)
pull_stats_airxbcal_clear_fs_scanangle(year, filter, cfg);
