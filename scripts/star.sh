#!/bin/bash

# SET PARAMS
	read1="$1"
	read2="$2"
#----------------------

mkdir -p "$mapbase"
mkdir -p "$mapbase"/"$sample_name"

STAR \
--runThreadN "$threads" \
--genomeDir "$star_index" --genomeLoad NoSharedMemory --sjdbGTFfile "$anno" --sjdbOverhang "$overhang" \
--readFilesIn "$read1" "$read2" --readFilesCommand zcat \
--outSAMtype BAM Unsorted --outSAMunmapped Within \
--outFileNamePrefix "$mapbase/$sample_name/" --outTmpDir "$mapbase/$sample_name/"STARtemp \
--outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 \
--chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50

samtools view -b -f 4 -o "$mapbase/$sample_name/$sample_name"_unmapped.bam "$mapbase/$sample_name/"Aligned.out.bam &
samtools view -b -F 4 -o "$mapbase/$sample_name/$sample_name"_mapped.bam "$mapbase/$sample_name/"Aligned.out.bam
samtools sort -@ "$threads" -m 1G -T "$mapbase/$sample_name/$sample_name".tmp -O bam -o "$mapbase/$sample_name/$sample_name".bam "$mapbase/$sample_name/$sample_name"_mapped.bam
samtools index "$mapbase/$sample_name/$sample_name".bam
rm -f "$mapbase/$sample_name/$sample_name"_mapped.bam "$mapbase/$sample_name/"Aligned.out.bam
