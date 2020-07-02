from pathlib import Path
import itertools
from snakemake.utils import min_version

min_version("5.4.0")

# ======================================================
# Config files
# ======================================================
configfile: "config.yaml"

# ======================================================
# Functions
# ======================================================


# ======================================================
# Variables
# ======================================================
sample_data_dir = config["sample_data_dir"]
msas_dir = config["msas_dir"]
prgs_dir = config["prgs_dir"]
analysis_output_dir = config["analysis_output_dir"]
msas_csv = config["msas_csv"]
samples = config["samples"]
make_prg_timeout_in_second = int(config["make_prg_timeout_in_second"])


output_files = []
for technology, sample, coverage, strategy in itertools.product(
                                                    config["technologies"],
                                                    config["samples"],
                                                    config["coverages"],
                                                    config["subsample"]["strategies"]):
    output_files.extend([
        f"{analysis_output_dir}/{technology}/{coverage}x/{strategy}/{sample}/map_with_discovery",
    ])
    output_files.extend([
        f"{analysis_output_dir}/{technology}/{coverage}x/{strategy}/prgs/denovo_updated.prg.fa",
    ])
    output_files.extend([
        f"{analysis_output_dir}/{technology}/{coverage}x/{strategy}/compare_nodenovo_local_genotyping/pandora_multisample_genotyped_local.vcf",
        f"{analysis_output_dir}/{technology}/{coverage}x/{strategy}/compare_nodenovo_global_genotyping/pandora_multisample_genotyped_global.vcf",
        f"{analysis_output_dir}/{technology}/{coverage}x/{strategy}/compare_withdenovo_local_genotyping/pandora_multisample_genotyped_local.vcf",
        f"{analysis_output_dir}/{technology}/{coverage}x/{strategy}/compare_withdenovo_global_genotyping/pandora_multisample_genotyped_global.vcf",
    ])


# deduplicate
output_files = list(set(output_files))

# ======================================================
# Rules
# ======================================================
rule all:
    input:
         output_files

include: "rules/map_with_denovo.smk"
include: "rules/get_denovo_updated_prg.smk"
include: "rules/compare.smk"