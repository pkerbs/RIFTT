#!/bin/bash

# DEFINE GLOBAL PARAMETERS
	export ref="/genome.fa"
	export anno="/anno.gtf"
	export star_index="/star_index"
	export fcdata="/fcdata"
			
	# Make folders accessible through environment variables
		outputfolder="/outputfolder"
		export fcbase="$outputfolder/fusioncatcher"
		export trimbase="$outputfolder/trimmedfastq"
		export mapbase="$outputfolder/mapping"
		export arbase="$outputfolder/arriba"
		export featCbase="$outputfolder/featurecounts"
		export isbase="$outputfolder/insertsizes"
		
	# Make tools available as environment variables
		export PATH="/tools/arriba_v2.2.1/:$PATH"
		export PATH="/tools/fusioncatcher-1.33/bin/:$PATH"
		export PATH="/tools/subread-2.0.3-Linux-x86_64/bin/:$PATH"
		export PATH="/tools/samtools-1.14/:$PATH"
		export PATH="/tools/STAR-2.7.10a/bin/Linux_x86_64/:$PATH"
		export picard="/tools/picard/picard.jar"
		export fastp="/tools/fastp.0.23.1"
		
	# Add tools for FusionCatcher to PATH
		export PATH="/tools/blat/:$PATH"
		export PATH="/tools/bbmap/:$PATH"
		export PATH="/tools/bowtie-1.2.3-linux-x86_64/:$PATH"
		export PATH="/tools/bowtie2-2.3.5.1-linux-x86_64/:$PATH"
		export PATH="/tools/picard/:$PATH"
		export PATH="/tools/seqtk-1.2-r101c/:$PATH"
		export PATH="/tools/fastqtk-0.27/:$PATH"
			
	# Add parameters for Arriba to PATH
		case $genomebuild in
			hg19)
				genome="hg19_hs37d5_GRCh37"
			;;
			hg38)
				genome="hg38_GRCh38"
			;;
		esac
		export blacklist="/tools/arriba_v2.2.1/database/blacklist_"$genome"_v2.2.1.tsv.gz"
		export known_fusions="/tools/arriba_v2.2.1/database/known_fusions_"$genome"_v2.2.1.tsv.gz"
		export prot_domains="/tools/arriba_v2.2.1/database/protein_domains_"$genome"_v2.2.1.gff3"	
#-----------------------------
