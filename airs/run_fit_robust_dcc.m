function run_fit_robust_dcc()
set_process_dirs;
addpath(genpath(rtp_sw_dir));

% grab the slurm array index for this process
slurmindex = str2num(getenv('SLURM_ARRAY_TASK_ID'));

lati = slurmindex + 1;
disp(lati)
pause(slurmindex)
fit_robust_series_dcc_plat(lati);
