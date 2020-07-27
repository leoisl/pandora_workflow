from pathlib import Path
from scripts.utils import get_technology_param


rule create_tsv_for_reads:
    input:
        expand(sample_data_dir + "/{sample}/{sample}.{{coverage}}x.{{sub_strategy}}.{{technology}}.fastq", sample=config["samples"])
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


rule compare_nodenovo:
    input:
        read_index=rules.create_tsv_for_reads.output.tsv,
        prg=config["original_prg"],
        prg_index=config["original_prg"] + ".k15.w14.idx",
    output:
        vcf=    analysis_output_dir+"/{technology}/{coverage}x/{sub_strategy}/compare_nodenovo_{genotyping_mode}_genotyping/pandora_multisample_genotyped_{genotyping_mode}.vcf",
        vcf_ref=analysis_output_dir+"/{technology}/{coverage}x/{sub_strategy}/compare_nodenovo_{genotyping_mode}_genotyping/pandora_multisample.vcf_ref.fa",
    threads: 54
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 50000
    params:
        pandora=config["pandora_executable"],
        log_level="info",
        outdir=lambda wildcards, output: str(Path(output.vcf).parent),
        technology_param = lambda wildcards: get_technology_param(wildcards)
    log:
        "logs/compare_nodenovo/{technology}/{coverage}x/{sub_strategy}/{genotyping_mode}.log"
    shell:
        """
        {params.pandora} compare \
            --prg_file {input.prg} \
            --read_index {input.read_index} \
            --outdir {params.outdir} \
            -t {threads} \
            --genotype {wildcards.genotyping_mode} \
            --max_covg 100000 \
            {params.technology_param} \
            --log_level {params.log_level} > {log} 2>&1
        """
