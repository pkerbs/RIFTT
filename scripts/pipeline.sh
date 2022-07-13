#!/bin/bash

run="$(date +"%Y-%m-%d_%H-%M-%S")"

# SET PARAMS
	source "/scripts/params.sh"
	mkdir -p "$outputfolder/pipeline_logs"
	steps=($tasks)
	if [[ ${steps[0]} -eq 1 && ${steps[2]} -eq 1 ]];then export threads=$((threads/2));fi
#----------------------

# FUSIONCATCHER
	if [[ ${steps[0]} -eq 1 ]]; then
		logfile="$outputfolder/pipeline_logs/$sample_name"_"$run"_fusioncatcher.log
		/scripts/fusioncatcher.sh "/read1.fastq.gz" "/read2.fastq.gz" >"$logfile" 2>&1 &
	fi
#----------------------

# FASTP
	if [[ ${steps[1]} -eq 1 ]]; then
		logfile="$outputfolder/pipeline_logs/$sample_name"_"$run"_fastp.log
		/scripts/fastp.sh "/read1.fastq.gz" "/read2.fastq.gz" >"$logfile" 2>&1
	fi
#----------------------

# STAR mapping
	if [[ ${steps[2]} -eq 1 ]]; then
		
		# If trimming by FastP was performed, take trimmed fastq files and
		# get read length from fastp report for --sjdbOverhang parameter
		if [[ ${steps[1]} -eq 1 ]]; then
			read1="$trimbase/$sample_name"_R1.fastq.gz
			read2="$trimbase/$sample_name"_R2.fastq.gz
			overhang=$(grep "mean length after filtering" "$trimbase"/"$sample_name"_fastp_report.html)
			overhang=${overhang/*, /""}
			overhang=${overhang/bp<*>/""}
			
		# Else take untrimmed input and get longest read for --sjdbOverhang parameter
		else
			read1="/read1.fastq.gz"
			read2="/read2.fastq.gz"
			overhang=$(zcat "$read1" "$read2" | awk '{if(NR%4==2 && length>max){max = length} } END {print max}')
		fi
		
		export overhang=$((overhang - 1))
		logfile="$outputfolder/pipeline_logs/$sample_name"_"$run"_star.log
		/scripts/star.sh "$read1" "$read2" >"$logfile" 2>&1
	fi
#----------------------

# ARRIBA
	if [[ ${steps[3]} -eq 1 ]]; then
		logfile="$outputfolder/pipeline_logs/$sample_name"_"$run"_arriba.log
		/scripts/arriba.sh >"$logfile" 2>&1 &
	fi
#----------------------

# FEATURECOUNTS
	if [[ ${steps[4]} -eq 1 ]]; then
		logfile="$outputfolder/pipeline_logs/$sample_name"_"$run"_featurecounts.log
		/scripts/featurecounts.sh >"$logfile" 2>&1 &
	fi
#----------------------

# INSERT SIZE ESTIMATION
	if [[ ${steps[5]} -eq 1 ]]; then
		logfile="$outputfolder/pipeline_logs/$sample_name"_"$run"_insertsize.log
		/scripts/insertsize.sh >"$logfile" 2>&1
	fi
#----------------------

wait # wait for everything to finish