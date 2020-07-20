from pathlib import Path
import sys
sys.path.append(str(Path().absolute()))

from utils import *
from typing import List, TextIO
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



def main():
    old_msa = Path(snakemake.input.msa)
    denovo_dirs = snakemake.params.denovo_dirs
    appended_msa = Path(snakemake.output.appended_msa)

    if not appended_msa.parent.is_dir():
        logging.debug("Creating parent directory for new MSA.")
        appended_msa.parent.mkdir(parents=True, exist_ok=True)

    denovo_paths = get_denovo_path_filepaths(denovo_dirs)

    if not denovo_paths:
        raise Exception("Error: no denovo paths found.")
    else:
        with appended_msa.open("w") as fh_out:
            append_denovo_paths_to_msa(denovo_paths, fh_out, old_msa)

main()
