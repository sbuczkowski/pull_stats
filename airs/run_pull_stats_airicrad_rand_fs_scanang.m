function run_pull_stats_airicrad_rand_fs_scanang(filter)
set_process_dirs;
addpath(genpath(rtp_sw_dir));

% grab the slurm array index for this process
slurmindex = str2num(getenv('SLURM_ARRAY_TASK_ID'));

%cfg.rtpdir = '/home/sbuczko1/WorkingFiles/rtp_airicrad_v6';
%cfg.statsdir = '/home/sbuczko1/WorkingFiles/data/stats/airs';
cfg.klayers = true;
cfg.latlimit = 50; % +-50

year = 2002+slurmindex;
disp(year)
disp(filter)
pull_stats_airicrad_rand_fs_scanangle(year, filter, cfg);
