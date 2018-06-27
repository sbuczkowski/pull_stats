function run_pull_stats_random(filter)
set_process_dirs;
addpath(genpath(rtp_sw_dir));
addpath /home/sbuczko1/git/pull_stats/cris

% grab the slurm array index for this process
slurmindex = str2num(getenv('SLURM_ARRAY_TASK_ID'));

year = 2012+slurmindex;
disp(year)
pull_stats_cris_random(year,filter);
