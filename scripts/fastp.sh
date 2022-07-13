#!/bin/bash

# SET PARAMS
	read1="$1"
	read2="$2"
#----------------------

mkdir -p "$trimbase"

"$fastp" -i "$read1" -I "$read2" \
-o "$trimbase"/"$sample_name"_R1.fastq.gz -O "$trimbase"/"$sample_name"_R2.fastq.gz \
--unpaired1 "$trimbase"/"$sample_name"_R1_unpaired.fastq.gz --unpaired2 "$trimbase"/"$sample_name"_R2_unpaired.fastq.gz \
-5 -3 -M 13 -W 1 -l 35 -p \
-R "$sample_name report" -h "$trimbase"/"$sample_name"_fastp_report.html -j "$trimbase"/"$sample_name"_fastp_report.json
