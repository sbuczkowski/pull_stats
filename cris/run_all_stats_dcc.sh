#!/bin/bash

# makes sbatch calls to run run_pull_stats for all extant filter combinations
#
# currently this makes sense for random subset only

for i in {1,4}; do
    sbatch run_pull_stats_dcc.sh $i
done




    
    
