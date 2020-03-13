from pathlib import Path
from rules.utils import get_technology_param

rule index_prg_updated_with_denovo_paths:
    input:
        get_denovo_updated_prg(analysis_output_dir+"/{technology}/{coverage}x/{sub_strategy}/prgs/denovo_updated.prg.fa"),
    output:
        index=analysis_output_dir+"/{technology}/{coverage}x/{sub_strategy}/prgs/denovo_updated.prg.fa.k15.w14.idx",
        kmer_prgs=directory(analysis_output_dir+"/{technology}/{coverage}x/{sub_strategy}/prgs/kmer_prgs"),
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000
    params:
        pandora=config["pandora_executable"],
    log:
        "logs/index_prg_updated_with_denovo_paths/{technology}/{coverage}x/{sub_strategy}.log"
    shell:
        "{params.pandora} index -t {threads} {input} > {log} 2>&1"



rule create_tsv_for_reads:
    input:
        [map_with_denovo(file for file in
         expand(sample_data_dir + "/{sample}/{sample}.{{coverage}}x.{{sub_strategy}}.{{technology}}.fastq", sample=config["samples"]))]
    output:
        tsv = sample_data_dir + "/samples.{coverage}x.{sub_strategy}.{technology}.tsv"
    threads: 1
    resources:
        mem_mb=200
    log:
        "logs/create_tsv_for_reads/{coverage}x/{sub_strategy}/{technology}.log"
    shell:
        """
        for path in {input}
        do
            filename=$(basename $path)
            sample_name=${{filename/.fastq}}
            echo -e \"$sample_name\t$(realpath $path)\" >> {output.tsv} 2>> {log}
        done
        """

rule compare_with_denovo:
    input:
        read_index=rules.create_tsv_for_reads.output.tsv,
        prg=get_denovo_updated_prg(analysis_output_dir+"/{technology}/{coverage}x/{sub_strategy}/prgs/denovo_updated.prg.fa"),
        prg_index=rules.index_prg_updated_with_denovo_paths.output.index,
    output:
        vcf=analysis_output_dir+"/{technology}/{coverage}x/{sub_strategy}/compare_with_denovo/pandora_multisample_genotyped.vcf",
        vcf_ref=analysis_output_dir+"/{technology}/{coverage}x/{sub_strategy}/compare_with_denovo/pandora_multisample.vcf_ref.fa",
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 30000
    params:
        pandora=config["pandora_executable"],
        log_level="debug",
        outdir=lambda wildcards, output: str(Path(output.vcf).parent),
        technology_param = lambda wildcards: get_technology_param(wildcards)
    log:
        "logs/compare_with_denovo/{technology}/{coverage}x/{sub_strategy}.log"
    shell:
        """
        {params.pandora} compare --prg_file {input.prg} \
            --read_index {input.read_index} \
            --outdir {params.outdir} \
            -t {threads} \
            --genotype \
            --max_covg 100000 \
            {params.technology_param} \
            --log_level {params.log_level} > {log} 2>&1
        """

rule compare_no_denovo:
    input:
        read_index=rules.create_tsv_for_reads.output.tsv,
        prg=config["original_prg"],
        prg_index=map_with_denovo(rules.index_original_prg.output.index),
    output:
        vcf=analysis_output_dir+"/{technology}/{coverage}x/{sub_strategy}/compare_no_denovo/pandora_multisample_genotyped.vcf",
        vcf_ref=analysis_output_dir+"/{technology}/{coverage}x/{sub_strategy}/compare_no_denovo/pandora_multisample.vcf_ref.fa",
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 30000
    params:
        pandora=config["pandora_executable"],
        log_level="debug",
        outdir=lambda wildcards, output: str(Path(output.vcf).parent),
        technology_param = lambda wildcards: get_technology_param(wildcards)
    log:
        "logs/compare_no_denovo/{technology}/{coverage}x/{sub_strategy}.log"
    shell:
        """
        {params.pandora} compare \
            --prg_file {input.prg} \
            --read_index {input.read_index} \
            --outdir {params.outdir} \
            -t {threads} \
            --genotype \
            --max_covg 100000 \
            {params.technology_param} \
            --log_level {params.log_level} > {log} 2>&1
        """
