from snakemake.shell import shell
from typing import List, TextIO
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


def build_prg_after_adding_denovo_paths(
    make_prg_script: str,
    max_nesting_lvl: int,
    prefix: str,
    msa: str,
    prg: TextIO,
    gene: str,
):
    logging.info("Building PRG for MSA.")
    shell(
        f"python3 {make_prg_script} -v --max_nesting {max_nesting_lvl} --prefix {prefix} {msa}"
    )
    logging.info("Finished building PRG.")
    logging.debug("Adding header info to PRG and renaming to correct filepath.")
    prg.write(f">{gene}\n")
    tmp_prg = Path(f"{prefix}.max_nest{max_nesting_lvl}.min_match7.prg")
    prg.write(tmp_prg.read_text() + "\n")


def main():
    updated_msa = Path(snakemake.input.updated_msa)
    prg = Path(snakemake.output.prg)

    with prg.open("w") as prg_fh:
        build_prg_after_adding_denovo_paths(
            snakemake.params.make_prg_script,
            snakemake.params.max_nesting_lvl,
            snakemake.params.prefix,
            str(updated_msa),
            prg_fh,
            snakemake.wildcards.gene,
        )


main()
