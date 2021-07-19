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
#SBATCH --job-name=RUN_LR_PULL_STATS
# partition = dev/batch
#SBATCH --partition=high_mem
# qos = short/normal/medium/long/long_contrib
#SBATCH --qos=normal+
#SBATCH --account=pi_strow
#SBATCH -N1
#SBATCH --mem=18000
#SBATCH --cpus-per-task 1
#SBATCH --time=01:30:00
# low res has data from 2012 to present: 4 years
#SBATCH --array=0-7
#SBATCH --requeue

#SBATCH -o /home/sbuczko1/LOGS/sbatch/run_cris_pull_stats-%A_%a.out
#SBATCH -e /home/sbuczko1/LOGS/sbatch/run_cris_pull_stats-%A_%a.err

# matlab options
MATLAB=matlab
MATOPT=' -nojvm -nodisplay -nosplash'

echo "Executing srun of run_pull_stats"
$MATLAB $MATOPT -r "include_addpaths_DEV; cfg=ini2struct('$1'); run_pull_stats($2,cfg); exit"
    
echo "Finished with srun of run_pull_stats"



