rule subsample_nanopore:
    input:
        reads=sample_data_dir+"/{sample}/{sample}.nanopore.fastq.gz",
        ref=sample_data_dir+"/{sample}/{sample}.ref.fa",
    output:
        subsampled_reads = sample_data_dir+"/{sample}/{sample}.{coverage}x.{sub_strategy}.nanopore.fastq"
    params:
        mean_q_weight=config['subsample']['mean_q_weight'],
        min_length=config['subsample']['min_length'],
        seed=config['subsample']['seed'],
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 1000 * attempt
    singularity: config["container"]
    log:
        "logs/subsample/{sub_strategy}/nanopore.{sample}.{coverage}x.log"
    shell:
        """
        bash scripts/downsample_reads.sh \
            nanopore \
            {input.reads} \
            {input.ref} \
            {wildcards.coverage} \
            {output.subsampled_reads} \
            {wildcards.sub_strategy} \
            {params.min_length} \
            {params.mean_q_weight} \
            {params.seed} 2> {log}  
        """


rule subsample_PE_illumina:
    input:
        reads_1=sample_data_dir+"/{sample}/{sample}.illumina_1.fastq.gz",
        reads_2=sample_data_dir+"/{sample}/{sample}.illumina_2.fastq.gz",
        ref=sample_data_dir+"/{sample}/{sample}.ref.fa",
    output:
        subsampled_reads_1 = sample_data_dir+"/{sample}/{sample}.{coverage}x.random.illumina.1.fastq",
        subsampled_reads_2 = sample_data_dir+"/{sample}/{sample}.{coverage}x.random.illumina.2.fastq"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 1000 * attempt
    log:
        "logs/subsample_PE_illumina/illumina.{sample}.{coverage}x.log"
    shell:
        """
        bash scripts/downsample_reads.sh \
            illumina_PE \
            {input.reads_1} \
            {input.reads_2} \
            {input.ref} \
            {wildcards.coverage} \
            {output.subsampled_reads_1} \
            {output.subsampled_reads_2} 
        """


rule concat_both_subsampled_PE_illumina_reads:
    input:
         subsampled_reads_1 = rules.subsample_PE_illumina.output.subsampled_reads_1,
         subsampled_reads_2 = rules.subsample_PE_illumina.output.subsampled_reads_2
    output:
         subsampled_reads = sample_data_dir+"/{sample}/{sample}.{coverage}x.random.illumina.fastq"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 1000 * attempt
    log:
        "logs/concat_both_subsampled_PE_illumina_reads/illumina.{sample}.{coverage}x.log"
    shell:
        "cat {input.subsampled_reads_1} {input.subsampled_reads_2} > {output.subsampled_reads}"


rule subsample_SE_illumina:
    input:
        reads=sample_data_dir+"/{sample}/{sample}.illumina.fastq.gz",
        ref=sample_data_dir+"/{sample}/{sample}.ref.fa",
    output:
        subsampled_reads = sample_data_dir+"/{sample}/{sample}.{coverage}x.random.illumina.fastq",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 1000 * attempt
    log:
        "logs/subsample_SE_illumina/illumina.{sample}.{coverage}x.log"
    shell:
        """
        bash scripts/downsample_reads.sh \
            illumina_SE \
            {input.reads} \
            {input.ref} \
            {wildcards.coverage} \
            {output.subsampled_reads} 
        """

# try to subsample paired before of single
ruleorder: subsample_PE_illumina > concat_both_subsampled_PE_illumina_reads > subsample_SE_illumina

rule index_original_prg:
    input:
        config["original_prg"]
    output:
        index=config["original_prg"] + ".k15.w14.idx",
        kmer_prgs=directory(prgs_dir + "/kmer_prgs")
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000
    log:
        "logs/index_original_prg.log"
    singularity: config["container"]
    shell:
        "pandora index -t {threads} {input} > {log} 2>&1"


rule map_with_discovery:
    input:
        prg=config["original_prg"],
        index=rules.index_original_prg.output.index,
        reads=sample_data_dir+"/{sample}/{sample}.{coverage}x.{sub_strategy}.{technology}.fastq",
        ref=sample_data_dir+"/{sample}/{sample}.ref.fa",
    output:
        outdir=directory(analysis_output_dir + "/{technology}/{coverage}x/{sub_strategy}/{sample}/map_with_discovery")
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 30000,
    params:
        log_level="info",
        use_discover=True,
    log:
        "logs/map_with_discovery/{technology}/{coverage}x/{sub_strategy}/{sample}.log"
    singularity: config["container"]
    shell:
        """
        bash scripts/pandora_map.sh pandora {input.prg} \
            {input.reads} \
            {output.outdir} \
            {threads} \
            {wildcards.technology} \
            {params.log_level} \
            {log} \
            {params.use_discover} \
            {input.ref}
        """
