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
#SBATCH --partition=high_mem
# qos = short/normal/medium/long/long_contrib
#SBATCH --qos=normal+
#SBATCH --account=pi_strow
#SBATCH -N1
#SBATCH --mem=18000
#SBATCH --cpus-per-task 1
#SBATCH --time=03:59:00
# airicrad has data from 2002 to present: 14 years
#SBATCH --array=0

#SBATCH --mail-user=sbuczko1@umbc.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=TIME_LIMIT_50

#SBATCH -o /home/sbuczko1/LOGS/sbatch/run_airs_pull_stats-%A_%a.out
#SBATCH -e /home/sbuczko1/LOGS/sbatch/run_airs_pull_stats-%A_%a.err

# matlab options
MATLAB=matlab
MATOPT=' -nojvm -nodisplay -nosplash'
which matlab

echo "Executing srun of run_airs_pull_stats"
$MATLAB $MATOPT -r "addpath('/home/sbuczko1/git/pull_stats_DEV/airs');\
                    addpath('/home/sbuczko1/git/pull_stats_DEV/util/');\
                    addpath('/home/sbuczko1/git/swutils');\
                    cfg=ini2struct('$2');\
                    run_pull_stats($1, cfg);\
                    exit"
    
echo "Finished with srun of run_airs_pull_stats"



