#!/bin/bash

mkdir -p "$arbase"
mkdir -p "$arbase"/"$sample_name"

arriba -x "$mapbase"/"$sample_name"/"$sample_name".bam \
-o "$arbase/$sample_name"/fusions.tsv -O "$arbase/$sample_name"/fusions.discarded.tsv \
-a "$ref" -g "$anno" -b "$blacklist" -k "$known_fusions" -t "$known_fusions" -p "$prot_domains" \
#-d structural_variants_from_WGS.tsv
