#!/bin/bash
#
# usage: sbatch run_pull_stats_random_template

# sbatch options
#SBATCH --job-name=RUN_AIRS_PULL_STATS
# partition = dev/batch
#SBATCH --partition=batch
# qos = short/normal/medium/long/long_contrib
#SBATCH --qos=medium
#SBATCH --account=pi_strow
#SBATCH -N1
#SBATCH --mem=18000
#SBATCH --cpus-per-task 1
#SBATCH --time=04:00:00
# airxbcal has data from 2002 to present: 14 years
##SBATCH --array=0-13

#SBATCH --mail-user=sbuczko1@umbc.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=TIME_LIMIT_50

#SBATCH -o /home/sbuczko1/logs/sbatch/run_airs_pull_stats-%A_%a.out
#SBATCH -e /home/sbuczko1/logs/sbatch/run_airs_pull_stats-%A_%a.err

# matlab options
MATLAB=/usr/cluster/matlab/current/bin/matlab
MATOPT=' -nojvm -nodisplay -nosplash'

echo "Executing sbatch run of run_pull_stats_random_template"
$MATLAB $MATOPT -r "loadAddPaths_template; run_pull_stats_random_template(); exit"
    
echo "Finished with sbatch run of run_pull_stats_random_template"



