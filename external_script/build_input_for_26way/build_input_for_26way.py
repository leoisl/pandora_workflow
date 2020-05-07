from glob import glob
from pathlib import Path
import os
import argparse

def create_symlink_if_it_does_not_already_exist(src, dest):
    if dest.exists():
        print(f"Skipping creating symlink for {dest}: already exists")
    else:
        print(f"Creating symlink {src} -> {dest}")
        os.symlink(src, dest)


def create_genomes(input_genomes_dir, output_dir):
    output_dir=Path(output_dir).absolute()
    output_dir.mkdir(parents=True, exist_ok=True)
    for file in glob(f"{input_genomes_dir}/*.fasta"):
        file = Path(file).absolute()
        genome = file.with_suffix("").name
        genome_dir = output_dir / genome
        genome_dir.mkdir(exist_ok=True)
        create_symlink_if_it_does_not_already_exist(file, genome_dir / f"{genome}.ref.fa")

def create_illumina_data(input_illumina_dir, output_dir):
    output_dir = Path(output_dir).absolute()
    output_dir.mkdir(parents=True, exist_ok=True)
    for file in glob(f"{input_illumina_dir}/*.fastq.gz"):
        file = Path(file).absolute()
        filename = file.with_suffix("").with_suffix("").name
        genome = filename[:-2]
        illumina_pair = filename[-1]
        genome_dir = output_dir / genome
        genome_dir.mkdir(exist_ok=True)
        create_symlink_if_it_does_not_already_exist(file, genome_dir / f"{genome}.illumina_{illumina_pair}.fastq.gz")

def create_nanopore_data(input_nanopore_dir, output_dir):
    output_dir = Path(output_dir).absolute()
    output_dir.mkdir(parents=True, exist_ok=True)
    for file in glob(f"{input_nanopore_dir}/*.fastq.gz"):
        file = Path(file).absolute()
        filename = file.with_suffix("").with_suffix("").name
        genome = filename
        genome_dir = output_dir / genome
        genome_dir.mkdir(exist_ok=True)
        create_symlink_if_it_does_not_already_exist(file, genome_dir / f"{genome}.nanopore.fastq.gz")


def get_args():
    parser = argparse.ArgumentParser(description='Build the input for the 26-way.')
    parser.add_argument('--genomes_dir', type=str, required=True)
    parser.add_argument('--illumina_dir', type=str, required=True)
    parser.add_argument('--nanopore_dir', type=str, required=True)
    parser.add_argument('--output_dir', type=str, required=True)
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    create_genomes(args.genomes_dir, args.output_dir)
    create_illumina_data(args.illumina_dir, args.output_dir)
    create_nanopore_data(args.nanopore_dir, args.output_dir)

if __name__ == "__main__":
    main()