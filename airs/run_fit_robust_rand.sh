#!/bin/bash
#

# sbatch options
#SBATCH --job-name=RUN_FIT_ROBUST_RAND
# partition = dev/batch
#SBATCH --partition=batch
# qos = short/normal/medium/long/long_contrib
#SBATCH --qos=short
#SBATCH --account=pi_strow
#SBATCH -N1
#SBATCH --mem-per-cpu=18000
#SBATCH --cpus-per-task 1
#SBATCH --time=00:20:00
# airxbcal has data from 2002 to present: 14 years

# matlab options
MATLAB=/usr/cluster/matlab/current/bin/matlab
MATOPT=' -nojvm -nodisplay -nosplash'

echo "Executing srun of run_fit_robust"
$MATLAB $MATOPT -r "addpath('~/git/pull_stats/airs/'); run_fit_robust_rand(); exit"
    
echo "Finished with srun of run_fit_robust"



