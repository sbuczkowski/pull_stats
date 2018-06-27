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

# sbatch options
#SBATCH --job-name=RUN_IASI_PULL_STATS_DCC
# partition = dev/batch
#SBATCH --partition=batch
# qos = short/normal/medium/long/long_contrib
#SBATCH --qos=medium
#SBATCH --account=pi_strow
#SBATCH -N1
#SBATCH --mem=18000
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
##SBATCH --array=0-9

#SBATCH --mail-user=sbuczko1@umbc.edu
##SBATCH --mail-type=BEGIN
##SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=TIME_LIMIT_50

#SBATCH -o /home/sbuczko1/logs/sbatch/pull_stats_iasi_dcc-%A_%a.out
#SBATCH -e /home/sbuczko1/logs/sbatch/pull_stats_iasi_dcc-%A_%a.err

# matlab options
MATLAB=/usr/cluster/matlab/current/bin/matlab
MATOPT=' -nojvm -nodisplay -nosplash'

JOBSTEP=0

echo "Executing srun of run_pull_stats_dcc"
$MATLAB $MATOPT -r "addpath('~/git/pull_stats/iasi'); run_pull_stats_dcc($1); exit"
    
echo "Finished with srun of run_pull_stats_dcc"


