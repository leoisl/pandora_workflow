from pathlib import Path
import sys
sys.path.append(str(Path().absolute()))

from utils import *
import subprocess
from typing import List, TextIO
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


def get_denovo_path_filepaths(denovo_dirs: List[str]) -> List[Path]:
    logging.debug("Getting filepaths for all de novo paths.")
    denovo_paths = []
    for denovo_dir in denovo_dirs:
        denovo_paths.extend(
            list(Path(denovo_dir).glob(f"{snakemake.wildcards.gene}.*.fa"))
        )

    logging.debug(f"Got {len(denovo_paths)} de novo paths")
    return denovo_paths


def append_denovo_paths_to_msa(
    denovo_paths: List[Path], fh_out: TextIO, old_msa: Path
) -> None:
    logging.info("Appending de novo paths to old MSA.")
    fh_out.write(old_msa.read_text())

    for p in denovo_paths:
        read_counter = 1
        sample = p.parts[-4]
        name = p.with_suffix("").name

        with p.open() as fasta:
            for line in fasta:
                if line.startswith(">"):
                    fh_out.write(
                        f"{line.rstrip()}_sample={sample}_{name}_path{read_counter}\n"
                    )
                    read_counter += 1
                else:
                    fh_out.write(line)

    logging.info("De novo paths appended to old MSA.")


def run_msa_after_adding_denovo_paths(
    appended_msa: str, updated_msa: str, threads: int, gene: str, original_prg: Path, timeout_in_seconds: int, msa_run_status_filename: str
):
    with open(msa_run_status_filename, "w") as msa_run_status_fh:
        try:
            logging.info("Running multiple sequence alignment.")
            subprocess.check_call(
                f"clustalo --dealign --threads {threads} --in {appended_msa} --out {updated_msa}",
                shell=True,
                timeout=timeout_in_seconds
            )
            msa_run_status_fh.write("SUCCESS\n")
            logging.info("Multiple sequence alignment finished.")
        except subprocess.TimeoutExpired:
            # we timed-out
            logging.info("[WARNING]: Timeout reached, copying previous PRG.")

            # let's just touch the new MSA to tell snakemake things are fine, they were "done"
            shell(f"touch {updated_msa}")

            # log this
            msa_run_status_fh.write("FAIL : TimeoutExpired\n")

            # now we make files to skip make_prg
            # let's just copy the previous PRG
            prg_filename = updated_msa.replace(".clustalo.fa", ".prg.fa")
            with original_prg.open() as original_prg_fh, open(prg_filename, "w") as prg:
                get_PRGs_from_original_PRG_restricted_to_list_of_genes(original_prg_fh, prg, [gene])

            # log this also
            prg_run_status_filename = msa_run_status_filename.replace("msas_run_status", "prgs_run_status")
            with open(prg_run_status_filename, "w") as prg_run_status_fh:
                prg_run_status_fh.write("FAIL : Clustalo\n")




def main():
    old_msa = Path(snakemake.input.msa)
    denovo_dirs = snakemake.params.denovo_dirs
    updated_msa = Path(snakemake.output.updated_msa)
    appended_msa = Path(snakemake.output.appended_msa)
    gene = snakemake.wildcards.gene
    original_prg = Path(snakemake.params.original_prg)
    timeout_in_seconds = snakemake.params.clustalo_timeout_in_second
    run_status = snakemake.output.run_status

    if not appended_msa.parent.is_dir():
        logging.debug("Creating parent directory for new MSA.")
        appended_msa.parent.mkdir(parents=True, exist_ok=True)

    denovo_paths = get_denovo_path_filepaths(denovo_dirs)

    if not denovo_paths:
        raise Exception("Error: no denovo paths found.")
    else:
        with appended_msa.open("w") as fh_out:
            append_denovo_paths_to_msa(denovo_paths, fh_out, old_msa)

        run_msa_after_adding_denovo_paths(
            str(appended_msa), str(updated_msa), snakemake.threads, gene, original_prg, timeout_in_seconds, run_status
        )

main()
