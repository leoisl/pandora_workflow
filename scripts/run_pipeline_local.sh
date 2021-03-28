#!/usr/bin/env bash
set -eux

LOG_DIR=logs/
mkdir -p "$LOG_DIR"
snakemake --use-singularity "$@"
exit 0
