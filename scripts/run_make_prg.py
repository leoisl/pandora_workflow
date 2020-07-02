import subprocess
from pathlib import Path
import logging
log_level = snakemake.params.log_level
logging.basicConfig(
    filename=str(snakemake.log),
    filemode="w",
    level=log_level,
    format="[%(asctime)s]:%(levelname)s: %(message)s",
    datefmt="%d/%m/%Y %I:%M:%S %p",
)
from scripts.utils import *


def build_prg_after_adding_denovo_paths(
    make_prg_script: str,
    max_nesting_lvl: int,
    prefix: str,
    msa: str,
    prg: TextIO,
    gene: str,
    original_prg: Path,
    timeout_in_seconds: int
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

    except subprocess.TimeoutExpired:
        # we timed-out, let's just copy the previous PRG
        logging.info("[WARNING]: Timeout reached, copying previous PRG.")
        with original_prg.open() as original_prg_fh:
            get_PRGs_from_original_PRG_restricted_to_list_of_genes(original_prg_fh, prg, [gene])




def main():
    updated_msa = Path(snakemake.input.updated_msa)
    prg = Path(snakemake.output.prg)
    original_prg = Path(snakemake.params.original_prg)

    with prg.open("w") as prg_fh:
        build_prg_after_adding_denovo_paths(
            snakemake.params.make_prg_script,
            snakemake.params.max_nesting_lvl,
            snakemake.params.prefix,
            str(updated_msa),
            prg_fh,
            snakemake.wildcards.gene,
            original_prg,
            snakemake.params.make_prg_timeout_in_second,
        )


main()
