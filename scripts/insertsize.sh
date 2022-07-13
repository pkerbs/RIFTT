#!/bin/bash

mkdir -p "$isbase"

java -jar "$picard" \
CollectInsertSizeMetrics \
I="$mapbase"/"$sample_name"/"$sample_name".bam \
O="$isbase"/"$sample_name"_is_metrics.txt \
H="$isbase"/"$sample_name"_is_histogram.pdf \
M=0.5