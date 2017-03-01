#!/bin/bash

# makes sbatch calls to run run_pull_stats for all extant filter combinations
#
# currently this makes sense for random subset only

for i in {1..6}; do
    sbatch --comment="MERRA first of month timeseries test" --array=0-14 run_pull_stats_airibrad_merraTseries.sh $i
done




    
    
