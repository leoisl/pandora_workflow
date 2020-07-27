#!/usr/bin/env bash
set -eux
JOB_NAME="snakemake_master_process."$(date --iso-8601='minutes')
LOG_DIR=logs/
MEMORY=4000
PROFILE="lsf"

mkdir -p "$LOG_DIR"
bsub -R "select[mem>$MEMORY] rusage[mem=$MEMORY] span[hosts=1]" \
    -M "$MEMORY" \
    -o "$LOG_DIR"/"$JOB_NAME".o \
    -e "$LOG_DIR"/"$JOB_NAME".e \
    -J "$JOB_NAME" \
      snakemake --snakefile Snakefile_compare_nodenovo --profile "$PROFILE" --stats "$LOG_DIR"/Snakefile_compare_nodenovo.stats "$@"

