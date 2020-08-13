gene_list = []
with open("/hps/nobackup2/iqbal/leandro/gal_msa/panx_vs_gal/panx_unmapped_genes.txt") as fin:
    for line in fin:
        gene_list.append(line.strip())

PRG_id_to_sequence = {}
id = None
sequence = None
with open("/hps/nobackup/iqbal/leandro/pandora_analysis_pipeline/data_24_way/data_for_pipeline/prgs/ecoli_pangenome_PRG_210619.fa") as fin:
    for line in fin:
        if id is None:
            id = line[1:].strip()
        else:
            sequence = line[:-1]
            PRG_id_to_sequence[id] = sequence
            id = None

for gene in gene_list:
    print(f">{gene}\n{PRG_id_to_sequence[gene]}")
