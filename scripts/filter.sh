#!/bin/bash

mkdir -p /outputfolder/filter_results
resultfolder=/outputfolder/filter_results/run_"$(date +"%Y-%m-%d_%H-%M-%S")"
mkdir -p "$resultfolder"

# Run the filtering
Rscript /scripts/Rscripts/workflow.R "$resultfolder" "$debug_flag" "$internal_BL" "$threads"
#----------------------
