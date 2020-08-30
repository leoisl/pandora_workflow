#!/usr/bin/env bash
set -eux

PROFILE="lsf"
LOG_DIR=logs/
LOCAL_CORES=16

snakemake --snakefile Snakefile_map_with_denovo        --profile "$PROFILE" --local-cores "$LOCAL_CORES" --stats "$LOG_DIR"/Snakefile_map_with_denovo.stats "$@"
sleep 60 # avoid locked directory issues
snakemake --snakefile Snakefile_get_denovo_updated_prg --profile "$PROFILE" --local-cores "$LOCAL_CORES" --stats "$LOG_DIR"/Snakefile_get_denovo_updated_prg.stats --keep-going "$@"
sleep 60 # avoid locked directory issues
snakemake --snakefile Snakefile_compare                --profile "$PROFILE" --local-cores "$LOCAL_CORES" --stats "$LOG_DIR"/Snakefile_compare.stats "$@"
exit 0
