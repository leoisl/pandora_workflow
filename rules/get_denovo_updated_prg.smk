from pathlib import Path
import pandas as pd
import fileinput
from scripts.utils import *


def get_genes_with_denovo_paths(analysis_output_dir, technology, coverage, sub_strategy, samples):
    denovo_dirs = [f"{analysis_output_dir}/{technology}/{coverage}x/{sub_strategy}/{sample}/map_with_discovery/denovo_paths" for sample in samples]
    denovo_dirs = [Path(denovo_dir) for denovo_dir in denovo_dirs]
    genes = set()
    for denovo_dir in denovo_dirs:
        for file in denovo_dir.glob("*.fa"):
            gene = file.name.split(".")[0]
            genes.add(gene)
    return genes


def get_genes_without_denovo_paths(genes_with_denovo_paths, msas_csv):
    msa_paths_as_str = pd.read_csv(msas_csv)["msas_absolute_paths"]
    msa_paths = [Path(msa_path_as_str) for msa_path_as_str in msa_paths_as_str]
    all_genes = {p.name.replace(".fa", "") for p in msa_paths}
    assert len(msa_paths) == len(all_genes)
    genes_without_denovo_paths = all_genes - genes_with_denovo_paths
    return genes_without_denovo_paths


def aggregate_prgs_with_denovo_path_input(wildcards):
    genes_with_denovo_paths = get_genes_with_denovo_paths(analysis_output_dir, wildcards.technology, wildcards.coverage,
                                                          wildcards.sub_strategy, samples)
    input_files = []
    for gene in genes_with_denovo_paths:
        tool = "custom"
        input_files.append(
            f"{analysis_output_dir}/{wildcards.technology}/{wildcards.coverage}x/{wildcards.sub_strategy}/prgs/{tool}/{gene}.prg.fa"
        )
    return input_files


rule aggregate_prgs_without_denovo_path:
    input:
        map_with_discovery_dirs = expand(analysis_output_dir+"/{{technology}}/{{coverage}}x/{{sub_strategy}}/{sample}/map_with_discovery", sample=config["samples"])
    output:
        prgs_without_denovo_paths = analysis_output_dir+"/{technology}/{coverage}x/{sub_strategy}/prgs/denovo_updated.prgs_without_denovo_paths.fa",
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    params:
        original_prg = config["original_prg"]
    log:
        "logs/aggregate_prgs_without_denovo_path/{technology}/{coverage}x/{sub_strategy}/.log"
    run:
        genes_with_denovo_paths = get_genes_with_denovo_paths(analysis_output_dir, wildcards.technology, wildcards.coverage,
                                                          wildcards.sub_strategy, samples)
        genes_without_denovo_paths = get_genes_without_denovo_paths(genes_with_denovo_paths, msas_csv)

        get_PRGs_from_original_PRG_restricted_to_list_of_genes(params.original_prg, output.prgs_without_denovo_paths,
                                                               genes_without_denovo_paths)



rule run_clustalo_after_adding_MSA_path:
    input:
        map_with_discovery_dirs = expand(analysis_output_dir+"/{{technology}}/{{coverage}}x/{{sub_strategy}}/{sample}/map_with_discovery", sample=config["samples"]),
        msa = msas_dir + "/{clustering_tool}/{gene}.fa"
    output:
        updated_msa = analysis_output_dir+"/{technology}/{coverage}x/{sub_strategy}/msas/{clustering_tool}/{gene}.clustalo.fa",
        appended_msa = analysis_output_dir+"/{technology}/{coverage}x/{sub_strategy}/msas/{clustering_tool}/{gene}.fa",
    threads: 8
    shadow: "shallow"
    resources:
        mem_mb = lambda wildcards, attempt: {1: 4000, 2: 16000, 3: 32000}.get(attempt, 64000)
    params:
        denovo_dirs = lambda wildcards, input: [map_with_discovery_dir+"/denovo_paths"
                                                for map_with_discovery_dir in input.map_with_discovery_dirs]
    singularity: config["make_prg_dependencies_img"]
    log:
        "logs/run_clustalo_after_adding_MSA_path/{technology}/{coverage}x/{sub_strategy}/{clustering_tool}/{gene}.log"
    script:
        "../scripts/run_clustalo_after_adding_MSA_path.py"


rule run_make_prg:
    input:
        updated_msa = rules.run_clustalo_after_adding_MSA_path.output.updated_msa
    output:
        prg = analysis_output_dir+"/{technology}/{coverage}x/{sub_strategy}/prgs/{clustering_tool}/{gene}.prg.fa"
    threads: 1
    shadow: "shallow"
    resources:
        mem_mb = lambda wildcards, attempt: {1: 4000, 2: 16000, 3: 32000}.get(attempt, 64000)
    params:
        log_level = "DEBUG",
        make_prg_script = "scripts/make_prg_from_msa.py",
        max_nesting_lvl = config.get("max_nesting_lvl", 5),
        prefix = lambda wildcards, output: output.prg.replace("".join(Path(output.prg).suffixes), ""),
        original_prg = config["original_prg"],
    singularity: config["make_prg_dependencies_img"]
    log:
        "logs/run_make_prg/{technology}/{coverage}x/{sub_strategy}/{clustering_tool}/{gene}.log"
    script:
        "../scripts/run_make_prg.py"


def concatenate_several_prgs_into_one(input_prgs, output_prg):
    with open(output_prg, "w") as fout, fileinput.input(input_prgs) as fin:
        for line in fin:
            if is_header(line):
                fout.write(line)
            else:
                prg_sequence = get_PRG_sequence(line)
                fout.write(prg_sequence + "\n")


rule aggregate_prgs_with_denovo_path:
    input:
        map_with_discovery_dirs = expand(analysis_output_dir+"/{{technology}}/{{coverage}}x/{{sub_strategy}}/{sample}/map_with_discovery", sample=config["samples"]),
        prgs = aggregate_prgs_with_denovo_path_input,
    output:
        prgs_with_denovo_paths = analysis_output_dir+"/{technology}/{coverage}x/{sub_strategy}/prgs/denovo_updated.prgs_with_denovo_paths.fa",
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/aggregate_prgs_with_denovo_path/{technology}/{coverage}x/{sub_strategy}/.log"
    run:
        concatenate_several_prgs_into_one(input.prgs, output.prgs_with_denovo_paths)


rule aggregate_prgs:
    input:
        prgs_with_denovo_paths = rules.aggregate_prgs_with_denovo_path.output.prgs_with_denovo_paths,
        prgs_without_denovo_paths = rules.aggregate_prgs_without_denovo_path.output.prgs_without_denovo_paths
    output:
        prg = analysis_output_dir+"/{technology}/{coverage}x/{sub_strategy}/prgs/denovo_updated.prg.fa",
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/aggregate_prgs/{technology}/{coverage}x/{sub_strategy}/.log"
    params:
        original_prg = config["original_prg"]
    run:
        concatenate_several_prgs_into_one([input.prgs_with_denovo_paths, input.prgs_without_denovo_paths],
                                          output.prg)

        # check original prg and new prg have the same number of sequences
        prgs_in_original = 0
        with open(params.original_prg) as fh:
            for line in fh:
                if line.startswith(">"):
                    prgs_in_original += 1

        prgs_in_new = 0
        with open(output.prg) as fh:
            for line in fh:
                if line.startswith(">"):
                    prgs_in_new += 1

        assert prgs_in_original == prgs_in_new, f"Original PRG ({prgs_in_original}) and new PRG ({prgs_in_new}) dont have the same number of entries!"
