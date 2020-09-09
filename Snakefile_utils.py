from pathlib import Path
import pandas as pd

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

def get_nanopore_reads(subsampled_reads, sample_name, subsampling, coverage):
    assert sample_name in subsampled_reads.sample_id.to_list()
    sample_path = subsampled_reads[subsampled_reads.sample_id == sample_name]["sample_path"].tolist()[0]
    return f"{sample_path}/{sample_name}.{coverage}x.{subsampling}.nanopore.fastq"

def get_uncompressed_reference(references, reference_id):
    assert reference_id in references.reference_id.to_list()
    uncompressed_reference_path = references[references.reference_id == reference_id]["uncompressed_file"].tolist()[0]
    return uncompressed_reference_path

def get_assembly(samples_df, sample_name):
    assert sample_name in samples_df.sample_id.to_list()
    sample_path = samples_df[samples_df.sample_id == sample_name]["sample_path"].tolist()[0]
    return f"{sample_path}/{sample_name}.ref.fa"

def get_config_vars(config):
    output_folder = config['output_folder']
    original_prg = config['original_prg']
    original_prg = str(Path(original_prg).absolute())

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
    msas_csv = config["msas_csv"]
    msas = pd.read_csv(msas_csv)
    msas = update_to_absolute_path(msas, ["msa"])

    make_prg_timeout_in_second = int(config["make_prg_timeout_in_second"])
    make_prg_memory_limit = int(config["make_prg_memory_limit"])
    clustalo_timeout_in_second = int(config["clustalo_timeout_in_second"])

    pandora_container = config["pandora_container"]

    ignore_adding_denovo_paths_for_these_genes = config["ignore_adding_denovo_paths_for_these_genes"]

    return output_folder, original_prg, coverages, subsamplings, technologies, samples_df, samples, pandora_container, \
            make_prg_timeout_in_second, make_prg_memory_limit, clustalo_timeout_in_second, msas_dir, msas, subsampled_reads, ignore_adding_denovo_paths_for_these_genes