function run_pull_stats_airxbcal_dcc(filter)
set_process_dirs;
addpath(genpath(rtp_sw_dir));

% grab the slurm array index for this process
slurmindex = str2num(getenv('SLURM_ARRAY_TASK_ID'));

year = 2002+slurmindex;
disp(year)
disp(filter)
pull_stats_airxbcal_dcc(year,filter);
