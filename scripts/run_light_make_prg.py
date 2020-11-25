from pathlib import Path
import sys
sys.path.append(str(Path().absolute()))
import os
from multiprocessing import Pool

processes = snakemake.threads
MSAs = snakemake.input.MSAs
output_folder = snakemake.output.prgs
max_nesting_lvl = snakemake.params.max_nesting_lvl
min_match_length = snakemake.params.min_match_length

def infer_gene_from_msa(msa):
    MSA_path = Path(msa)
    MSA_path = MSA_path.with_suffix("")
    return MSA_path.name


def run_make_prg(msa):
    gene = infer_gene_from_msa(msa)
    os.system(f"make_prg from_msa --max_nesting {max_nesting_lvl} --min_match_length {min_match_length} "
              f"--prefix {gene} {msa}")
    os.makedirs(output_folder, exist_ok=True)
    os.rename(f"{gene}.max_nest{max_nesting_lvl}.min_match{min_match_length}.prg", f"{output_folder}/{gene}.prg.fa")


with Pool(processes=processes) as pool:
    pool.starmap(run_make_prg, [MSAs])
