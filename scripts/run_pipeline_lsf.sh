#!/usr/bin/env bash
set -eux

PROFILE="lsf"
snakemake --profile "$PROFILE" "$@"