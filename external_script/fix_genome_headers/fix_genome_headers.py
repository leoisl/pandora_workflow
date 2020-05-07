from glob import glob
from pathlib import Path
import os
import argparse


def fix_genome_headers(input_genomes_dir, output_genomes_dir):
    output_genomes_dir=Path(output_genomes_dir).absolute()
    output_genomes_dir.mkdir(parents=True, exist_ok=True)

    for file in glob(f"{input_genomes_dir}/*.fasta"):
        print(f"Processing {file}...")
        input_file = Path(file).absolute()
        output_file = output_genomes_dir / input_file.name
        with open(input_file) as input_file_fh, open(output_file, "w") as output_file_fh:
            header_index = 0
            for line in input_file_fh:
                if line.startswith(">"):
                    line = f">{header_index} {line[1:]}"
                    header_index+=1

                output_file_fh.write(line)


def get_args():
    parser = argparse.ArgumentParser(description='Fix genome headers, putting an id before the header itself.')
    parser.add_argument('--genomes_dir', type=str, required=True)
    parser.add_argument('--output_dir', type=str, required=True)
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    fix_genome_headers(args.genomes_dir, args.output_dir)

if __name__ == "__main__":
    main()