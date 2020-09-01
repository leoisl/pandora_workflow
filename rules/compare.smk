from pathlib import Path
from scripts.utils import get_technology_param

rule index_prg_updated_with_denovo_paths:
    input:
        output_folder+"/{technology}/{coverage}x/{sub_strategy}/prgs/denovo_updated.prg.fa",
    output:
        index=output_folder+"/{technology}/{coverage}x/{sub_strategy}/prgs/denovo_updated.prg.fa.k15.w14.idx",
        kmer_prgs=directory(output_folder+"/{technology}/{coverage}x/{sub_strategy}/prgs/kmer_prgs"),
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000
    singularity: pandora_container
    log:
        "logs/index_prg_updated_with_denovo_paths/{technology}/{coverage}x/{sub_strategy}.log"
    shell:
        "pandora index -t {threads} {input} > {log} 2>&1"


rule create_tsv_for_reads:
    input:
        reads = lambda wildcards: [get_reads(subsampled_reads, wildcards.technology, sample, wildcards.coverage, wildcards.sub_strategy) for sample in samples],
    output:
        tsv = output_folder+"/{technology}/{coverage}x/{sub_strategy}/reads.tsv"
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


rule compare_withdenovo:
    input:
        read_index=rules.create_tsv_for_reads.output.tsv,
        prg=output_folder+"/{technology}/{coverage}x/{sub_strategy}/prgs/denovo_updated.prg.fa",
        prg_index=rules.index_prg_updated_with_denovo_paths.output.index,
    output:
        vcf=    output_folder+"/{technology}/{coverage}x/{sub_strategy}/compare_withdenovo_{genotyping_mode}_genotyping/pandora_multisample_genotyped_{genotyping_mode}.vcf",
        vcf_ref=output_folder+"/{technology}/{coverage}x/{sub_strategy}/compare_withdenovo_{genotyping_mode}_genotyping/pandora_multisample.vcf_ref.fa",
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 30000
    params:
        log_level="info",
        outdir=lambda wildcards, output: str(Path(output.vcf).parent),
        technology_param = lambda wildcards: get_technology_param(wildcards)
    log:
        "logs/compare_withdenovo/{technology}/{coverage}x/{sub_strategy}/{genotyping_mode}.log"
    singularity: pandora_container
    shell:
        """
            pandora compare --prg_file {input.prg} \
            --read_index {input.read_index} \
            --outdir {params.outdir} \
            -t {threads} \
            --genotype {wildcards.genotyping_mode} \
            --max_covg 100000 \
            {params.technology_param} \
            --log_level {params.log_level} > {log} 2>&1
        """


rule compare_nodenovo:
    input:
        read_index = rules.create_tsv_for_reads.output.tsv,
        prg = output_folder + "/prgs/prg.fa",
        prg_index = output_folder + "/prgs/prg.fa.k15.w14.idx",
    output:
        vcf=    output_folder+"/{technology}/{coverage}x/{sub_strategy}/compare_nodenovo_{genotyping_mode}_genotyping/pandora_multisample_genotyped_{genotyping_mode}.vcf",
        vcf_ref=output_folder+"/{technology}/{coverage}x/{sub_strategy}/compare_nodenovo_{genotyping_mode}_genotyping/pandora_multisample.vcf_ref.fa",
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 30000
    params:
        log_level="info",
        outdir=lambda wildcards, output: str(Path(output.vcf).parent),
        technology_param = lambda wildcards: get_technology_param(wildcards)
    log:
        "logs/compare_nodenovo/{technology}/{coverage}x/{sub_strategy}/{genotyping_mode}.log"
    singularity: pandora_container
    shell:
        """
            pandora compare \
            --prg_file {input.prg} \
            --read_index {input.read_index} \
            --outdir {params.outdir} \
            -t {threads} \
            --genotype {wildcards.genotyping_mode} \
            --max_covg 100000 \
            {params.technology_param} \
            --log_level {params.log_level} > {log} 2>&1
        """
