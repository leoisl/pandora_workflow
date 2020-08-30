#!/usr/bin/env bash
set -eux

LOG_DIR=logs/
mkdir -p "$LOG_DIR"

snakemake --snakefile Snakefile_map_with_denovo        --use-singularity --stats "$LOG_DIR"/Snakefile_map_with_denovo.stats                     "$@" >"$LOG_DIR"/Snakefile_map_with_denovo.out 2>"$LOG_DIR"/Snakefile_map_with_denovo.err
sleep 60 # avoid locked directory issues
snakemake --snakefile Snakefile_get_denovo_updated_prg --use-singularity --stats "$LOG_DIR"/Snakefile_get_denovo_updated_prg.stats --keep-going "$@" >"$LOG_DIR"/Snakefile_get_denovo_updated_prg.out 2>"$LOG_DIR"/Snakefile_get_denovo_updated_prg.err
sleep 60 # avoid locked directory issues
snakemake --snakefile Snakefile_compare                --use-singularity --stats "$LOG_DIR"/Snakefile_compare.stats                             "$@" >"$LOG_DIR"/Snakefile_compare.out 2>"$LOG_DIR"/Snakefile_compare.err
exit 0
