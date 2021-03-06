function run_pull_stats_airicrad_clear(filter)
% grab the slurm array index for this process
slurmindex = str2num(getenv('SLURM_ARRAY_TASK_ID'));

cfg.klayers = true;
cfg.rtpdir = '/asl/rtp/rtp_airicrad_v6';
cfg.statsdir = '/asl/data/stats/airs/clear';

year = 2002+slurmindex;
disp(year)
disp(filter)
pull_stats_airicrad_clear(year, filter, cfg);
