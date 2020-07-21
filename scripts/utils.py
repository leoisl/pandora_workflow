from typing import TextIO, List
from datetime import datetime
import subprocess

def get_technology_param(wildcards):
    if wildcards.technology=="illumina":
        return "--illumina"
    else:
        return ""


def is_header(line):
    return line.startswith(">")


def get_gene(line):
    stripped_line = line.rstrip()
    gene = stripped_line[1:]
    return gene


def get_PRG_sequence(line):
    prg_sequence = line.rstrip()
    line_ends_digit = prg_sequence[-1].isdigit()
    if line_ends_digit:
        prg_sequence += " "
    return prg_sequence


def get_PRGs_from_original_PRG_restricted_to_list_of_genes(original_prg: TextIO, new_prg: TextIO, genes: List):
    for line in original_prg:
        if is_header(line):
            gene = get_gene(line)
        else:
            if gene in genes:
                prg_sequence = get_PRG_sequence(line)
                new_prg.write(">" + gene + "\n")
                new_prg.write(prg_sequence + "\n")


def get_number_of_sequences_in_fasta(fasta_filepath: str) -> int:
    with open(fasta_filepath) as fasta_filehandler:
        number_of_sequences_in_fasta = sum([1 for line in fasta_filehandler if line.startswith(">")])
        return number_of_sequences_in_fasta


def run_command_and_time_it(list_of_command_and_args, timeout_in_seconds=None):
    before = datetime.now()
    subprocess.check_call(list_of_command_and_args, timeout=timeout_in_seconds)
    after = datetime.now()
    time_spent_in_command = after - before
    return time_spent_in_command.total_seconds()