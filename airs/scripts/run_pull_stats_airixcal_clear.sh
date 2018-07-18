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
#SBATCH --job-name=RUN_AIRS_PULL_STATS
# partition = dev/batch
#SBATCH --partition=batch
# qos = short/normal/medium/long/long_contrib
#SBATCH --qos=medium
##SBATCH --account=pi_strow
#SBATCH -N1
#SBATCH --mem=18000
#SBATCH --cpus-per-task 1
#SBATCH --time=10:00:00
# airxbcal has data from 2002 to present: 14 years
#SBATCH --array=0-16

#SBATCH --mail-user=sbuczko1@umbc.edu
##SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
##SBATCH --mail-type=TIME_LIMIT_50

#SBATCH -o /home/sbuczko1/logs/sbatch/run_airs_pull_stats-%A_%a.out
#SBATCH -e /home/sbuczko1/logs/sbatch/run_airs_pull_stats-%A_%a.err

# matlab options
MATLAB=/usr/cluster/matlab/current/bin/matlab
MATOPT=' -nojvm -nodisplay -nosplash'

echo "Executing srun of run_airs_pull_stats"
$MATLAB $MATOPT -r "addpath('~/git/pull_stats/airs/'); run_pull_stats_airixcal_clear($1); exit"
    
echo "Finished with srun of run_airs_pull_stats"



