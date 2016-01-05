function run_pull_stats_random_template()

% grab the slurm array index for this process
slurmindex = str2num(getenv('SLURM_ARRAY_TASK_ID'));

year = 2002+slurmindex;
disp(year)
pull_stats_airibrad_rand(year);
