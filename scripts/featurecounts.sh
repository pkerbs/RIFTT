#!/bin/bash

mkdir -p "$featCbase"
mkdir -p "$featCbase/raw"
mkdir -p "$featCbase/reformatted"

featureCounts -a "$anno" \
-s "$strandness" \
-p \
--countReadPairs \
-Q 0 \
-C \
-M -O \
--fraction \
-t exon \
-g gene_id \
-T "$threads" \
-o "$featCbase/raw/$sample_name.fc" \
"$mapbase/$sample_name/$sample_name.bam"

# REFORMAT COUNTSFILE	
	awk 'NR!=1 && NR!=2' "$featCbase"/raw/"$sample_name".fc | awk -F '\t' 'BEGIN {OFS = FS} {print $1,$7}' > "$featCbase"/reformatted/"$sample_name".fc
