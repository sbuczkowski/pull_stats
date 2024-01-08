function run_pull_stats_scanang(filter, cfg)

% grab the slurm array index for this process
slurmindex = str2num(getenv('SLURM_ARRAY_TASK_ID'));

year = 2002+slurmindex;
disp(year)
pull_stats_airs_scanang(year,filter,cfg);
