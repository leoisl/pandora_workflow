#!/usr/bin/env bash
if [[ "$#" -lt 4 ]]; then
    echo "Error: Illegal number of parameters"
    echo -e "usage:\n$(basename "$0") <reads> <ref> <covg> <outname>"
    exit 1
fi

set -euv

reads="$1"
ref="$2"
covg="$3"
outname="$4"

genome_size=$(grep -v '^>' "$ref" | wc | awk '{print $3-$1}')
num_bases_to_keep=$((genome_size * covg))


python downsample_illumina_reads/downsample_illumina_se_reads.py --reads $reads \
--number_of_bases $num_bases_to_keep --out_reads $outname