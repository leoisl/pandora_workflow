import argparse
import pysam
import pandas as pd


def get_number_of_bases(file):
    with pysam.FastxFile(file) as fastx_file:
        number_of_bases = sum([len(entry.sequence) for entry in fastx_file])
    return number_of_bases


def get_args():
    parser = argparse.ArgumentParser(description='Get read coverage given assembly and reads.')
    parser.add_argument('--assembly_and_reads_csv', type=str, required=True)
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    assembly_and_reads_csv = pd.read_csv(args.assembly_and_reads_csv)

    for assembly, reads in zip(assembly_and_reads_csv["assembly"], assembly_and_reads_csv["reads"]):
        number_of_bases_in_assembly = get_number_of_bases(assembly)
        number_of_bases_in_read = get_number_of_bases(reads)
        coverage = number_of_bases_in_read / number_of_bases_in_assembly
        print(f"Coverage of {reads} ({number_of_bases_in_read} bps) on {assembly} ({number_of_bases_in_assembly} bps):")
        print(coverage)

if __name__ == "__main__":
    main()