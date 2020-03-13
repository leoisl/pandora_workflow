#!/usr/bin/env bash
set -eux

PROFILE="lsf"
LOG_DIR=logs/

mkdir -p "$LOG_DIR"
snakemake --snakefile Snakefile_map_with_denovo --profile "$PROFILE" --verbose --stats "$LOG_DIR"/Snakefile_map_with_denovo.stats "$@"
snakemake --snakefile Snakefile_get_denovo_updated_prg --profile "$PROFILE" --verbose --stats "$LOG_DIR"/Snakefile_get_denovo_updated_prg.stats "$@"
snakemake --snakefile Snakefile_compare --profile "$PROFILE" --verbose --stats "$LOG_DIR"/Snakefile_compare.stats "$@"
