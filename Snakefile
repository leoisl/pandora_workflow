import itertools
from pathlib import Path
import pandas as pd

configfile: "config.yaml"

########################################################################################################################
# helper functions
########################################################################################################################
def update_to_absolute_path_core(path_series):
    return path_series.apply(lambda path: str(Path(path).absolute()))
def update_to_absolute_path(df, columns):
    for column in columns:
        df[column] = update_to_absolute_path_core(df[column])
    return df

def get_reads(subsampled_reads, technology, sample, coverage, sub_strategy):
    assert sample in subsampled_reads.sample_id.to_list()
    sample_path = subsampled_reads[subsampled_reads.sample_id == sample]["subsampled_reads_dir"].tolist()[0]
    return f"{sample_path}/{sample}.{coverage}x.{sub_strategy}.{technology}.fastq"

def get_technology_param(technology):
    assert technology in ["illumina", "nanopore"]
    if technology == "illumina":
        return "--illumina"
    else:
        return ""
########################################################################################################################
########################################################################################################################


########################################################################################################################
# setup global vars
########################################################################################################################
output_folder = config['output_folder']
coverages = config["coverages"]
subsamplings = config["subsamplings"]
technologies = config["technologies"]

samples_csv = config["samples"]
samples_df = pd.read_csv(samples_csv)
samples_df = update_to_absolute_path(samples_df, ["sample_path"])
samples = samples_df.sample_id.to_list()

subsampled_reads_dir = config["subsampled_reads_dir"]
subsampled_reads = pd.read_csv(subsampled_reads_dir)
subsampled_reads = update_to_absolute_path(subsampled_reads, ["subsampled_reads_dir"])

msas_dir = config["msas_dir"]
pandora_URL = config["pandora_URL"]
make_prg_URL = config["make_prg_URL"]
mafft_container = config["containers"]["mafft"]
########################################################################################################################
########################################################################################################################


########################################################################################################################
# setup output files
########################################################################################################################
output_files = []
for technology, sample, coverage, subsampling in itertools.product(
                                                    technologies,
                                                    samples,
                                                    coverages,
                                                    subsamplings):
    output_files.extend([
        f"{output_folder}/{technology}/{coverage}x/{subsampling}/compare_nodenovo/pandora_multisample_genotyped.vcf",
        f"{output_folder}/{technology}/{coverage}x/{subsampling}/compare_withdenovo/pandora_multisample_genotyped.vcf",
    ])

# deduplicate
output_files = list(set(output_files))
########################################################################################################################
########################################################################################################################


# ======================================================
# Rules
# ======================================================
rule all:
    input:
         output_files

rule download_pandora:
    output:
        pandora_exec = output_folder + "/tools/pandora"
    params:
        outdir = lambda wildcards, output: str(Path(output.pandora_exec).parent),
    threads: 1
    resources:
        mem_mb=200
    log:
        "logs/download_pandora.log"
    shell:
        """
        (cd {params.outdir} && \
        wget {pandora_URL} && \
        gunzip pandora* && \
        mv pandora* pandora && \
        chmod +x ./pandora) >{log} 2>&1
        """


rule download_make_prg:
    output:
        make_prg_exec = output_folder + "/tools/make_prg"
    params:
        outdir = lambda wildcards, output: str(Path(output.make_prg_exec).parent),
    threads: 1
    resources:
        mem_mb=200
    log:
        "logs/download_make_prg.log"
    shell:
        """
        (cd {params.outdir} && \
        wget {make_prg_URL} && \
        mv make_prg* make_prg && \
        chmod +x ./make_prg) >{log} 2>&1
        """


rule create_read_index:
    input:
        reads = lambda wildcards: [get_reads(subsampled_reads, wildcards.technology, sample, wildcards.coverage, wildcards.sub_strategy) for sample in samples],
    output:
        read_index = output_folder+"/{technology}/{coverage}x/{sub_strategy}/reads.tsv"
    threads: 1
    resources:
        mem_mb=200
    log:
        "logs/create_read_index/{coverage}x/{sub_strategy}/{technology}.log"
    shell:
        """
        for path in {input}
        do
            filename=$(basename $path)
            sample_name=${{filename/.fastq}}
            echo -e \"$sample_name\t$(realpath $path)\" >> {output.read_index} 2>> {log}
        done
        """


rule make_prg_from_msa:
    input:
        msas_dir = msas_dir,
        make_prg_exec = rules.download_make_prg.output.make_prg_exec
    output:
        prg_file = output_folder+"/prgs/ecoli_pangenome_PRG.prg.fa",
        prg_folder = directory(output_folder+"/prgs/ecoli_pangenome_PRG_prgs"),
        update_DS = output_folder+"/prgs/ecoli_pangenome_PRG.update_DS",
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000
    params:
        output_prefix = output_folder+"/prgs/ecoli_pangenome_PRG"
    log:
        "logs/make_prg_from_msa.log"
    shell:
        "{input.make_prg_exec} from_msa --input {input.msas_dir} --output_prefix {params.output_prefix} -t {threads} >{log} 2>&1"


rule index_original_prg:
    input:
        prg = rules.make_prg_from_msa.output.prg_file,
        pandora_exec = rules.download_pandora.output.pandora_exec
    output:
        index = output_folder+"/prgs/ecoli_pangenome_PRG.prg.fa.k15.w14.idx",
        kmer_prgs = directory(output_folder+"/prgs/kmer_prgs")
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20000
    log:
        "logs/index_original_prg.log"
    shell:
        "{input.pandora_exec} index -t {threads} {input.prg} >{log} 2>&1"


rule pandora_discover:
    input:
        prg = rules.make_prg_from_msa.output.prg_file,
        index = rules.index_original_prg.output.index,
        reads_index = rules.create_read_index.output.read_index,
        pandora_exec = rules.download_pandora.output.pandora_exec
    output:
        outdir=directory(output_folder + "/{technology}/{coverage}x/{sub_strategy}/pandora_discover_out")
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20000,
    params:
        technology_param = lambda wildcards: get_technology_param(wildcards.technology)
    log:
        "logs/pandora_discover/{technology}/{coverage}x/{sub_strategy}/pandora_discover.log"
    shell:
        "{input.pandora_exec} discover --outdir {output.outdir} -t {threads} --max-covg 100000 {params.technology_param} "
        "{input.prg} {input.reads_index} >{log} 2>&1"


rule update_prg:
    input:
        update_DS = rules.make_prg_from_msa.output.update_DS,
        pandora_discover_out = rules.pandora_discover.output.outdir,
        make_prg_exec = rules.download_make_prg.output.make_prg_exec
    output:
        prg_file = output_folder+ "/{technology}/{coverage}x/{sub_strategy}/prgs_updated/ecoli_pangenome_PRG.prg.fa"
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 20000
    params:
        output_prefix = output_folder+ "/{technology}/{coverage}x/{sub_strategy}/prgs_updated/ecoli_pangenome_PRG"
    container: mafft_container
    log:
        "logs/{technology}/{coverage}x/{sub_strategy}/update_prg.log"
    shell:
        "{input.make_prg_exec} update --update_DS {input.update_DS} --denovo_paths {input.pandora_discover_out}/denovo_paths.txt "
        "--output_prefix {params.output_prefix} -t {threads} >{log} 2>&1"


rule index_updated_prg:
    input:
        prg = rules.update_prg.output.prg_file,
        pandora_exec = rules.download_pandora.output.pandora_exec
    output:
        index = output_folder+"/{technology}/{coverage}x/{sub_strategy}/prgs_updated/ecoli_pangenome_PRG.prg.fa.k15.w14.idx",
        kmer_prgs = directory(output_folder+"/{technology}/{coverage}x/{sub_strategy}/prgs_updated/kmer_prgs")
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20000
    log:
        "logs/{technology}/{coverage}x/{sub_strategy}/index_updated_prg.log"
    shell:
        "{input.pandora_exec} index -t {threads} {input.prg} >{log} 2>&1"


rule compare_withdenovo:
    input:
        read_index=rules.create_read_index.output.read_index,
        prg=rules.index_updated_prg.input.prg,
        prg_index=rules.index_updated_prg.output.index,
        pandora_exec = rules.download_pandora.output.pandora_exec
    output:
        vcf=    output_folder+"/{technology}/{coverage}x/{sub_strategy}/compare_withdenovo/pandora_multisample_genotyped.vcf",
        vcf_ref=output_folder+"/{technology}/{coverage}x/{sub_strategy}/compare_withdenovo/pandora_multisample.vcf_ref.fa",
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20000
    params:
        outdir=lambda wildcards, output: str(Path(output.vcf).parent),
        technology_param = lambda wildcards: get_technology_param(wildcards.technology)
    log:
        "logs/compare_withdenovo/{technology}/{coverage}x/{sub_strategy}/pandora_compare_withdenovo.log"
    shell:
        """
            {input.pandora_exec} compare  \
            --outdir {params.outdir} \
            --genotype \
            --max-covg 100000 \
            -t {threads} \
            {params.technology_param} \
             {input.prg} \
             {input.read_index} > {log} 2>&1
        """


rule compare_nodenovo:
    input:
        read_index=rules.create_read_index.output.read_index,
        prg=rules.index_original_prg.input.prg,
        prg_index=rules.index_original_prg.output.index,
        pandora_exec = rules.download_pandora.output.pandora_exec
    output:
        vcf=    output_folder+"/{technology}/{coverage}x/{sub_strategy}/compare_nodenovo/pandora_multisample_genotyped.vcf",
        vcf_ref=output_folder+"/{technology}/{coverage}x/{sub_strategy}/compare_nodenovo/pandora_multisample.vcf_ref.fa",
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20000
    params:
        outdir=lambda wildcards, output: str(Path(output.vcf).parent),
        technology_param = lambda wildcards: get_technology_param(wildcards.technology)
    log:
        "logs/compare_nodenovo/{technology}/{coverage}x/{sub_strategy}/pandora_compare_nodenovo.log"
    shell:
        """
            {input.pandora_exec} compare  \
            --outdir {params.outdir} \
            --genotype \
            --max-covg 100000 \
            -t {threads} \
            {params.technology_param} \
             {input.prg} \
             {input.read_index} > {log} 2>&1
        """
