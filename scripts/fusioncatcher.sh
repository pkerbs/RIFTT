#!/bin/bash

# SET PARAMS
	read1="$1"
	read2="$2"
#----------------------

mkdir -p "$fcbase"
mkdir -p "$fcbase"/"$sample_name"

fusioncatcher -d "$fcdata" \
-i "$read1","$read2" \
-p "$threads" \
-o "$fcbase/$sample_name" \
--aligners blat,star,bowtie2
