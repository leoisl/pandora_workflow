from snakemake import shell
from pathlib import Path

def main():
    previous_gene = Path(snakemake.input.previous_gene)
    msa = Path(snakemake.output.msa)
    run_MSA = snakemake.params.run_MSA
    threads = snakemake.threads

    if not run_MSA:
        shell(f"cp {previous_gene} {msa}")
    else:
        shell(f"clustalo --dealign --threads {threads} --in {previous_gene} --out {msa}")

main()