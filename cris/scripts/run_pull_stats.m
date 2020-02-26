function run_pull_stats(filter, cfg)

include_addpaths;

% grab the slurm array index for this process
slurmindex = str2num(getenv('SLURM_ARRAY_TASK_ID'));

switch cfg.instName
  case 'cris'
    baseyear = 2012;
  case 'cris2'
    baseyear = 2018;
  otherwise
    error('Improper instrument name in configuration')
end

year = baseyear+slurmindex;
disp(year)
pull_stats_cris(year,filter,cfg);
