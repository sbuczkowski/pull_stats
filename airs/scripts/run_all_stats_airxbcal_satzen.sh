#!/bin/bash

# makes sbatch calls to run run_pull_stats for all extant filter combinations
#
# currently this makes sense for random subset only
wKLayers=$1

for i in {1..2}; do
    sbatch run_pull_stats_airxbcal_satzen.sh $i $wKLayers
done




    
    
