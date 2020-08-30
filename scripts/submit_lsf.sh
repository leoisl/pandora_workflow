#!/usr/bin/env bash
set -eux

MEMORY=16000
LOCAL_CORES=16
LOG_DIR=logs/
JOB_NAME="snakemake_master_process."$(date --iso-8601='minutes')

mkdir -p "$LOG_DIR"

bsub -R "select[mem>$MEMORY] rusage[mem=$MEMORY] span[hosts=1]" \
    -M "$MEMORY" \
    -n "$LOCAL_CORES" \
    -o "$LOG_DIR"/"$JOB_NAME".o \
    -e "$LOG_DIR"/"$JOB_NAME".e \
    -J "$JOB_NAME" \
      bash scripts/run_pipeline_lsf.sh "$@"

