#!/usr/bin/env bash
set -eux

CORES=112
LOG_DIR=logs/


mkdir -p "$LOG_DIR"
snakemake -j "$CORES" --snakefile Snakefile_map_with_denovo --stats "$LOG_DIR"/Snakefile_map_with_denovo.stats --report "$LOG_DIR"/Snakefile_map_with_denovo.report.html  --use-singularity >"$LOG_DIR"/Snakefile_map_with_denovo.out 2>"$LOG_DIR"/Snakefile_map_with_denovo.err
sleep 60 # avoid locked directory issues
snakemake -j "$CORES" --snakefile Snakefile_get_denovo_updated_prg --stats "$LOG_DIR"/Snakefile_get_denovo_updated_prg.stats --report "$LOG_DIR"/Snakefile_get_denovo_updated_prg.report.html --use-singularity >"$LOG_DIR"/Snakefile_get_denovo_updated_prg.out 2>"$LOG_DIR"/Snakefile_get_denovo_updated_prg.err
sleep 60 # avoid locked directory issues
snakemake -j "$CORES" --snakefile Snakefile_compare --stats "$LOG_DIR"/Snakefile_compare.stats --report "$LOG_DIR"/Snakefile_compare.report.html --use-singularity >"$LOG_DIR"/Snakefile_compare.out 2>"$LOG_DIR"/Snakefile_compare.err
