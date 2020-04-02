#!/usr/bin/env bash

technology="$1"
shift #remove first arg

if [ "$technology" = "nanopore" ]
then
    bash scripts/downsample_nanopore_reads.sh "$@"
elif [ "$technology" = "illumina_PE" ]
then
    bash scripts/downsample_pe_illumina_reads.sh "$@"
elif [ "$technology" = "illumina_SE" ]
then
    bash scripts/downsample_se_illumina_reads.sh "$@"
else
    echo "Invalid technology given: $technology"
    exit 1
fi