function run_pull_stats_byFOR(year, filter)
addpath /home/sbuczko1/git/pull_stats/cris

% grab the slurm array index for this process
slurmindex = str2num(getenv('SLURM_ARRAY_TASK_ID'));

ifor = slurmindex;
disp(ifor)
pull_stats_cris_hires_byFOR(year, filter, ifor);
