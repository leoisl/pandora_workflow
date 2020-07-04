from pathlib import Path
import sys
sys.path.append(str(Path().absolute()))
from snakemake.utils import read_job_properties
from typing import TextIO

import subprocess
from utils import *
import logging
log_level = snakemake.params.log_level
logging.basicConfig(
    filename=str(snakemake.log),
    filemode="w",
    level=log_level,
    format="[%(asctime)s]:%(levelname)s: %(message)s",
    datefmt="%d/%m/%Y %I:%M:%S %p",
)



def build_prg_after_adding_denovo_paths(
    make_prg_script: str,
    max_nesting_lvl: int,
    prefix: str,
    msa: str,
    prg: TextIO,
    gene: str,
    original_prg: Path,
    timeout_in_seconds: int,
    mem_mb_given_to_job: int,
    max_mb_allowed: int,
    run_status_fh: TextIO
):
    try:
        logging.info("Building PRG for MSA.")
        subprocess.check_call(
            f"python3 {make_prg_script} -v --max_nesting {max_nesting_lvl} --prefix {prefix} {msa}",
            shell=True,
            timeout=timeout_in_seconds
        )
        logging.info("Finished building PRG.")

        #  Adding header info to PRG and renaming to correct filepath
        prg.write(f">{gene}\n")
        tmp_prg = Path(f"{prefix}.max_nest{max_nesting_lvl}.min_match7.prg")
        prg.write(tmp_prg.read_text() + "\n")
        run_status_fh.write("SUCCESS")

    except subprocess.TimeoutExpired:
        # we timed-out, let's just copy the previous PRG
        logging.info("[WARNING]: Timeout reached, copying previous PRG.")
        with original_prg.open() as original_prg_fh:
            get_PRGs_from_original_PRG_restricted_to_list_of_genes(original_prg_fh, prg, [gene])
        run_status_fh.write("FAIL : TimeoutExpired")
    except KeyboardInterrupt: # memory limit reached
        if mem_mb_given_to_job >= max_mb_allowed:
            logging.info("[WARNING]: Memory limit exceeded, copying previous PRG.")
            with original_prg.open() as original_prg_fh:
                get_PRGs_from_original_PRG_restricted_to_list_of_genes(original_prg_fh, prg, [gene])
            run_status_fh.write("FAIL : Memory Limit Exceeded")

def get_mem_mb_given_to_job():
    jobscript = sys.argv[1]
    job_properties = read_job_properties(jobscript)
    return int(job_properties["resources"]["mem_mb"])


def main():
    updated_msa = Path(snakemake.input.updated_msa)
    prg = Path(snakemake.output.prg)
    run_status = Path(snakemake.output.run_status)
    original_prg = Path(snakemake.params.original_prg)
    mem_mb_given_to_job = get_mem_mb_given_to_job()
    max_mb_allowed = int(snakemake.params.make_prg_memory_limit)


    with prg.open("w") as prg_fh, run_status.open("w") as run_status_fh:
        build_prg_after_adding_denovo_paths(
            snakemake.params.make_prg_script,
            snakemake.params.max_nesting_lvl,
            snakemake.params.prefix,
            str(updated_msa),
            prg_fh,
            snakemake.wildcards.gene,
            original_prg,
            snakemake.params.make_prg_timeout_in_second,
            mem_mb_given_to_job,
            max_mb_allowed,
            run_status_fh
        )


main()
