#!/usr/bin/env bash
set -eux

LOG_DIR=logs/
mkdir -p "$LOG_DIR"

snakemake --snakefile Snakefile_map_with_denovo        --use-singularity "$@"
sleep 60 # avoid locked directory issues
snakemake --snakefile Snakefile_get_denovo_updated_prg --use-singularity "$@"
sleep 60 # avoid locked directory issues
snakemake --snakefile Snakefile_compare                --use-singularity "$@"
exit 0
