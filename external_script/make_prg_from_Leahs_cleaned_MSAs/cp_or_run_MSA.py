from snakemake import shell
from pathlib import Path

def there_is_only_one_sequence(fasta_file):
    nb_of_seqs = 0
    with open(fasta_file) as fasta_file_fh:
        for line in fasta_file_fh:
            if line.startswith(">"):
                nb_of_seqs += 1
    return nb_of_seqs == 1


def main():
    previous_gene = Path(snakemake.input.previous_gene)
    msa = Path(snakemake.output.msa)
    run_MSA = snakemake.params.run_MSA
    threads = snakemake.threads

    if not run_MSA:
        shell(f"cp {previous_gene} {msa}")
    else:
        if there_is_only_one_sequence(previous_gene):
            shell(f"cp {previous_gene} {msa}")
        else:
            shell(f"clustalo --dealign --threads {threads} --in {previous_gene} --out {msa}")

main()