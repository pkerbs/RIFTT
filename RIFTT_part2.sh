#!/bin/bash

# REQUIRED PARAMETERS
	# General
		anno=""			# Genome annotation file (GENCODE)
		outputfolder=""	# Output folder of the detection pipeline (Part I)
		clintable=""	# Excel table with columns: cohort, sample, Karyotype, otherCyto
	# Optional
		internal_BL=1	# Whether to use internal blacklist of fusion genes
		user_BL=""		# Path to own fusion blacklist (xlsx file, first column with fusion labels)
#---------------------------------------------------------

#(DO NOT EDIT THIS SECTION)
# Run the filtering
	cmd="singularity exec"
	cmd+=" --bind $outputfolder:/outputfolder"
	cmd+=",$anno:/anno.gtf"
	cmd+=",$clintable:/clintable.xlsx"
	if [[ -f "$user_BL" ]]; then
		cmd+=",$user_BL:/blacklist_user.xlsx"
	fi
	cmd+=" --env debug_flag=$debug_flag,internal_BL=$internal_BL,threads=$threads"
	cmd+=" RIFTT.sif"
	cmd+=" /scripts/filter.sh"
	eval $cmd
#----------------------
