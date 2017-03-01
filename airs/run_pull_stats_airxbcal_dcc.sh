#!/bin/bash
#
# usage: sbatch run_pull_stats filter
#
# where filter = {1..6}
# 1 = descending (night), land & ocean
# 2 = desc, ocean
# 3 = desc, land
# 4 = ascending (day, land & ocean
# 5 = asc, ocean
# 6 = asc, land
# 7 = tropical, ocean

# sbatch options
#SBATCH --job-name=RUN_AIRS_PULL_STATS_DCC
# partition = dev/batch
#SBATCH --partition=batch
# qos = short/normal/medium/long/long_contrib
#SBATCH --qos=medium
#SBATCH --account=pi_strow
#SBATCH -N1
#SBATCH --mem-per-cpu=18000
#SBATCH --cpus-per-task 1
#SBATCH --time=02:00:00
# airxbcal has data from 2002 to present: 14 years
##SBATCH --array=13

# matlab options
MATLAB=/usr/cluster/matlab/current/bin/matlab
MATOPT=' -nojvm -nodisplay -nosplash'

echo "Executing srun of run_pull_stats"
$MATLAB $MATOPT -r "addpath('~/git/pull_stats/airs/'); run_pull_stats_airxbcal_dcc($1); exit"
    
echo "Finished with srun of run_pull_stats"



