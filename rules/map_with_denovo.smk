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
