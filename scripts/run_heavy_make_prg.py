from pathlib import Path
import sys
sys.path.append(str(Path().absolute()))
import os
import time
import logging
logging.basicConfig(
    format="%(asctime)s [%(levelname)s]: %(message)s",
    level=logging.INFO,
    filename=str(snakemake.log),
    filemode="w",
)

msa = snakemake.input.updated_msa
gene = snakemake.wildcards.gene
prg = snakemake.output.prg
max_nesting_lvl = snakemake.params.max_nesting_lvl
min_match_length = snakemake.params.min_match_length

start = time.time()
logging.info(f"Running make_prg for {gene}")
os.system(f"make_prg from_msa --max_nesting {max_nesting_lvl} --min_match_length {min_match_length} "
          f"--prefix {gene} {msa}")
os.rename(f"{gene}.max_nest{max_nesting_lvl}.min_match{min_match_length}.prg", prg)
stop = time.time()
runtime = stop - start
logging.info(f"Finished running make_prg for {gene}")
logging.info(f"make_prg runtime for heavy {gene} in seconds: {runtime:.3f}")


