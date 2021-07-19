function run_pull_stats_clear(filter)
addpath /home/sbuczko1/git/pull_stats_DEV/cris

% grab the slurm array index for this process
slurmindex = str2num(getenv('SLURM_ARRAY_TASK_ID'));
if isempty(slurmindex)
    slurmindex=0;
end

year = 2012+slurmindex;
disp(year)
pull_stats_cris_clear(year,filter);
