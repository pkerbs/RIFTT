#!/bin/bash
# RIFTT Part 1 - fusion detection for a single sample.
# Usage:  ./RIFTT_part1.sh <sample_name> <fastq_folder>
# Configuration is read from config/params.conf (copy config/params.conf.example).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${RIFTT_CONFIG:-$HERE/config/params.conf}"
SIF="${RIFTT_SIF:-$HERE/RIFTT.sif}"

[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG"; echo "  cp config/params.conf.example config/params.conf  and edit it."; exit 1; }
[[ -f "$SIF" ]]    || { echo "ERROR: container not found: $SIF (build it or set RIFTT_SIF)"; exit 1; }
# shellcheck disable=SC1090
source "$CONFIG"

sample_name="${1:?usage: RIFTT_part1.sh <sample_name> <fastq_folder>}"
fastq_folder="${2:?usage: RIFTT_part1.sh <sample_name> <fastq_folder>}"
read1=$(find "$fastq_folder" -maxdepth 1 -name "${sample_name}*R1.fastq.gz" | head -1)
read2=$(find "$fastq_folder" -maxdepth 1 -name "${sample_name}*R2.fastq.gz" | head -1)
[[ -f "$read1" && -f "$read2" ]] || { echo "ERROR: could not find ${sample_name}*R1/R2.fastq.gz in $fastq_folder"; exit 1; }

for v in REF ANNO STARINDEX FCDATA OUTPUTFOLDER; do
  [[ -n "${!v:-}" ]] || { echo "ERROR: $v is not set in $CONFIG"; exit 1; }
done
mkdir -p "$OUTPUTFOLDER"

tasks="$STEP_FUSIONCATCHER $STEP_FASTP $STEP_STAR $STEP_ARRIBA $STEP_FEATURECOUNTS $STEP_PICARD"

singularity run \
  --bind "$REF:/genome.fa,$ANNO:/anno.gtf,$STARINDEX:/star_index,$FCDATA:/fcdata,$OUTPUTFOLDER:/outputfolder,$read1:/read1.fastq.gz,$read2:/read2.fastq.gz" \
  --env "sample_name=$sample_name,strandness=$STRANDNESS,threads=$THREADS,genomebuild=$GENOMEBUILD" \
  --env "tasks=$tasks" \
  "$SIF"
