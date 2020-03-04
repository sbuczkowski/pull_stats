function run_pull_stats_random(filter)
% $$$ set_process_dirs;
% $$$ addpath(genpath(rtp_sw_dir));
addpath('/asl/packages/rtp_prod2/util')
addpath('~/git/pull_stats/iasi')

% grab the slurm array index for this process
slurmindex = str2num(getenv('SLURM_ARRAY_TASK_ID'));

year = 2007+slurmindex;
disp(year)
pull_stats_iasi_rand(year,filter);
