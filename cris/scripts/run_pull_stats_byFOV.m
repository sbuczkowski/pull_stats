function run_pull_stats_byFOV(ifov, cfg)
addpath /home/sbuczko1/git/pull_stats/cris

% grab the slurm array index for this process
slurmindex = str2num(getenv('SLURM_ARRAY_TASK_ID'));

year = 2012+slurmindex;
disp(year)
pull_stats_cris_hires_byFOV(year,ifov,cfg);
