rule index_original_prg:
    input:
        config["original_prg"]
    output:
        index=config["original_prg"] + ".k15.w14.idx",
        kmer_prgs=directory(prgs_dir + "/kmer_prgs")
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000
    params:
        pandora=config["pandora_executable"],
    log:
        "logs/index_original_prg.log"
    shell:
        "{params.pandora} index -t {threads} {input} > {log} 2>&1"

rule index_prg_updated_with_denovo_paths:
    input:
        analysis_output_dir+"/{technology}/{coverage}x/{sub_strategy}/prgs/denovo_updated.prg.fa",
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
