function run_pull_stats(filter, wKLayers)
set_process_dirs;
addpath(genpath(rtp_sw_dir));

% grab the slurm array index for this process
slurmindex = str2num(getenv('SLURM_ARRAY_TASK_ID'));

year = 2002+slurmindex;
disp(year)
disp(filter)
if wKLayers == true
    pull_stats_airibrad_rand_klayers(year,filter);
else
    pull_stats_airibrad_rand(year,filter);
end