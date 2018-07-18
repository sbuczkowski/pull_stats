function run_pull_stats_airixcal_clear(filter)
addpath ../

% grab the slurm array index for this process
slurmindex = str2num(getenv('SLURM_ARRAY_TASK_ID'));

cfg.klayers = true;

year = 2002+slurmindex;
disp(year)
disp(filter)
pull_stats_airixcal_clear(year, filter, cfg);
