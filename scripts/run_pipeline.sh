#!/usr/bin/env bash
set -eux

PROFILE="lsf"
LOG_DIR=logs/
LOCAL_CORES=16

mkdir -p "$LOG_DIR"
snakemake --snakefile Snakefile_map_with_denovo --profile "$PROFILE" --stats "$LOG_DIR"/Snakefile_map_with_denovo.stats "$@"
sleep 60 # avoid locked directory issues
snakemake --snakefile Snakefile_get_denovo_updated_prg --local-cores "$LOCAL_CORES" --profile "$PROFILE" --stats "$LOG_DIR"/Snakefile_get_denovo_updated_prg.stats "$@"
sleep 60 # avoid locked directory issues
snakemake --snakefile Snakefile_compare --profile "$PROFILE" --stats "$LOG_DIR"/Snakefile_compare.stats "$@"
