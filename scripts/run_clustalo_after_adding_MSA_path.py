from pathlib import Path
import sys
sys.path.append(str(Path().absolute()))

from utils import *
import subprocess
from typing import TextIO, List
from snakemake import shell
import logging
log_level = snakemake.params.log_level
logging.basicConfig(
    filename=str(snakemake.log),
    filemode="w",
    level=log_level,
    format="[%(asctime)s]:%(levelname)s: %(message)s",
    datefmt="%d/%m/%Y %I:%M:%S %p",
)


def run_msa_after_adding_denovo_paths(
    appended_msa: str, updated_msa: str, threads: int, timeout_in_seconds: int, run_status_fh: TextIO,
        gene:str, ignore_adding_denovo_paths_for_these_genes: List[str]):

    if gene in ignore_adding_denovo_paths_for_these_genes:
        # ignore
        logging.info(f"Ignoring running multiple sequence alignment on gene {gene}.")

        # let's just create and empty MSA to tell snakemake things are fine, they were "done"
        shell(f"echo '' > {updated_msa}")

        # log this
        run_status_fh.write("FAIL : Ignored\n")
    else:
        try:
            logging.info("Running multiple sequence alignment.")
            logging.info(f"Number of sequences: {get_number_of_sequences_in_fasta(appended_msa)}")

            time_spent = run_command_and_time_it(
                ["clustalo", "--dealign", "--threads", str(threads), "--in", appended_msa, "--out", updated_msa],
                timeout_in_seconds)
            logging.info(f"Time in seconds: {time_spent}")
            run_status_fh.write("SUCCESS\n")
            logging.info("Multiple sequence alignment finished.")
        except subprocess.TimeoutExpired:
            # we timed-out
            logging.info("[WARNING]: Timeout reached, making MSA empty.")

            # let's just create and empty MSA to tell snakemake things are fine, they were "done"
            shell(f"echo '' > {updated_msa}")

            # log this
            run_status_fh.write("FAIL : TimeoutExpired\n")


def main():
    appended_msa = Path(snakemake.input.appended_msa)
    updated_msa = Path(snakemake.output.updated_msa)
    timeout_in_seconds = snakemake.params.clustalo_timeout_in_second
    run_status = Path(snakemake.output.run_status)
    gene = snakemake.wildcards.gene
    ignore_adding_denovo_paths_for_these_genes = snakemake.params.ignore_adding_denovo_paths_for_these_genes

    with run_status.open("w") as run_status_fh:
        run_msa_after_adding_denovo_paths(
            str(appended_msa), str(updated_msa), snakemake.threads, timeout_in_seconds, run_status_fh,
            gene, ignore_adding_denovo_paths_for_these_genes)

main()
