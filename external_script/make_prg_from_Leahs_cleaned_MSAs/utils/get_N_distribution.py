# configs
input_dir = "/hps/nobackup2/iqbal/leandro/gal_msa/10perc_renamed_with_piggy"

from Bio import SeqIO
import os
import sys

print("locus,allele,number_of_Ns,number_of_bases")
for file in os.listdir(input_dir):
    print(f"Processing {file}...", file=sys.stderr)
    source_file_full_path = os.path.join(input_dir, file)
    for record in SeqIO.parse(source_file_full_path, "fasta"):
        allele = record.id
        number_of_Ns = record.seq.count("N")
        number_of_bases = len(record.seq)
        print(f"{file},{allele},{number_of_Ns},{number_of_bases}")