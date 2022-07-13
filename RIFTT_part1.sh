#!/bin/bash

# REQUIRED PARAMETERS
	# General
		threads=16			# Number of threads to run the tools of the pipeline
		outputfolder=""		# Common folder for the output of the pipeline
		genomebuild="hg38"	# ["hg19", "hg38"]
	
	# Sample information
		sample_name="$1"										# Name of sample
		fastq_folder="$2"										# Folder of fastq files. Forward/Reverse read files will be searched by sample_name
		read1=`find "$fastq_folder" -name "$sample_name"*R1*`	# First read (automatic detection)
		read2=`find "$fastq_folder" -name "$sample_name"*R2*`	# Second read (automatic detection)
		strandness=0											# [0 -> unstranded, 1 -> stranded, 2 -> reversely stranded]

	# Reference files
		ref=""			# Genome reference file (Fasta)
		anno=""			# Genome annotation file (GENCODE)
		starindex=""	# Index folder for STAR
		fcdata=""		# Data folder for FusionCatcher
	
	# Steps to perform
		FusionCatcher=1		# Fusion calling by FusionCatcher
		FastP=1				# Read trimming before STAR mapping
		STAR=1				# Mapping by STAR
		Arriba=1			# Fusion calling by Arriba
		FeatureCounts=1		# Read counting
		Picard=1 			# Insert size estimation
#---------------------------------------------------------

#(DO NOT EDIT THIS SECTION)
tasks="$FusionCatcher $FastP $STAR $Arriba $FeatureCounts $Picard"
cmd="singularity run"
cmd+=" --bind $ref:/genome.fa,$anno:/anno.gtf,$starindex:/star_index,$fcdata:/fcdata"
cmd+=",$outputfolder:/outputfolder,$read1:/read1.fastq.gz,$read2:/read2.fastq.gz"
cmd+=" --env sample_name=$sample_name,strandness=$strandness"
cmd+=",threads=$threads,genomebuild=$genomebuild,tasks='$tasks'"
cmd+=" RIFTT.sif"
eval $cmd
#---------------------------------------------------------
