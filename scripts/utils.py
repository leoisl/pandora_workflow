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


def get_PRGs_from_original_PRG_restricted_to_list_of_genes(original_prg_filename, new_prg_filename, genes):
    with open(original_prg_filename) as original_prg, open(new_prg_filename, "w") as new_prg:
        for line in original_prg:
            if is_header(line):
                gene = get_gene(line)
            else:
                if gene in genes:
                    prg_sequence = get_PRG_sequence(line)
                    new_prg.write(">" + gene + "\n")
                    new_prg.write(prg_sequence + "\n")
