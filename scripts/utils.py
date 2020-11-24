from pathlib import Path
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


class LocusComplexity:
    # caching/memoization
    __locus_complexity_records = None

    def __init__(self, gene, nb_of_seqs):
        self.gene = gene
        self.nb_of_seqs = nb_of_seqs

    @staticmethod
    def get_locus_complexity_records(updated_msas_dir: Path) -> List["LocusComplexity"]:
        # implements locus_complexity_records caching/memoization
        if LocusComplexity.__locus_complexity_records is None:
            LocusComplexity.__locus_complexity_records = []
            file_list = [file for file in Path(updated_msas_dir).iterdir() if file.is_file()]
            for file in file_list:
                gene = file.with_suffix("").name
                nb_of_seqs = get_number_of_sequences_in_fasta(str(file))
                LocusComplexity.__locus_complexity_records.append(LocusComplexity(gene=gene, nb_of_seqs=nb_of_seqs))

        return LocusComplexity.__locus_complexity_records

    def __str__(self):
        return f"{self.gene}: nb_of_seqs = {self.nb_of_seqs}"
    def __repr__(self):
        return str(self)
    def __lt__(self, other):
        return self.gene < other.gene
    def __eq__(self, other):
        return (self.gene, self.nb_of_seqs) == (other.gene, other.nb_of_seqs)


def get_light_MSAs(updated_msas_dir: Path, complex_MSA_sequence_threshold: int) -> List[str]:
    locus_complexity_records = LocusComplexity.get_locus_complexity_records(updated_msas_dir)
    light_records = filter(lambda locus_complexity_record: locus_complexity_record.nb_of_seqs <= complex_MSA_sequence_threshold,
                           locus_complexity_records)
    light_MSAs = [updated_msas_dir / f"{light_record.gene}.fa" for light_record in light_records]
    return [str(light_MSA) for light_MSA in light_MSAs]

def get_heavy_MSAs(updated_msas_dir: Path, complex_MSA_sequence_threshold: int) -> List[str]:
    locus_complexity_records = LocusComplexity.get_locus_complexity_records(updated_msas_dir)
    heavy_records = filter(lambda locus_complexity_record: locus_complexity_record.nb_of_seqs > complex_MSA_sequence_threshold,
                           locus_complexity_records)
    heavy_MSAs = [updated_msas_dir / f"{heavy_record.gene}.fa" for heavy_record in heavy_records]
    return [str(heavy_MSA) for heavy_MSA in heavy_MSAs]

def get_light_PRGs(updated_msas_dir: Path,
                   light_prgs_dir: Path,
                   complex_MSA_sequence_threshold: int) -> List[str]:
    locus_complexity_records = LocusComplexity.get_locus_complexity_records(updated_msas_dir)
    light_records = filter(
        lambda locus_complexity_record: locus_complexity_record.nb_of_seqs > complex_MSA_sequence_threshold,
        locus_complexity_records)
    light_PRGs = [light_prgs_dir / f"{light_record.gene}.prg.fa" for light_record in light_records]
    return [str(light_PRG) for light_PRG in light_PRGs]

def get_heavy_PRGs(updated_msas_dir: Path,
                   heavy_prgs_dir: Path,
                   complex_MSA_sequence_threshold: int) -> List[str]:
    locus_complexity_records = LocusComplexity.get_locus_complexity_records(updated_msas_dir)
    heavy_records = filter(
        lambda locus_complexity_record: locus_complexity_record.nb_of_seqs > complex_MSA_sequence_threshold,
        locus_complexity_records)
    heavy_PRGs = [heavy_prgs_dir / f"{heavy_record.gene}.prg.fa" for heavy_record in heavy_records]
    return [str(heavy_PRG) for heavy_PRG in heavy_PRGs]


# TODO: promote these to unit tests
MSA_path = Path("/home/leandro/git/myforks/pandora_analysis_pipeline/sample_data/msas/custom")

def test_get_locus_complexity_records():
    locus_complexity_records = LocusComplexity.get_locus_complexity_records(MSA_path)
    assert sorted(locus_complexity_records) == sorted([
        LocusComplexity("GC00000085_1", 16),
        LocusComplexity("GC00001725", 307),
        LocusComplexity("GC00001920", 307),
        LocusComplexity("GC00003727_5", 18),
        LocusComplexity("GC00005143", 29),
        LocusComplexity("GC00006430_1", 1),
        LocusComplexity("GC00007254", 6),
        LocusComplexity("GC00007867", 4),
        LocusComplexity("GC00008389_1", 1),
        LocusComplexity("GC00009588_1", 1),
    ])

def test_get_light_MSAs():
    light_MSAs = get_light_MSAs(MSA_path, 10)
    assert sorted(light_MSAs) == sorted([
        str(MSA_path/"GC00006430_1.fa"),
        str(MSA_path/"GC00007254.fa"),
        str(MSA_path/"GC00007867.fa"),
        str(MSA_path/"GC00008389_1.fa"),
        str(MSA_path/"GC00009588_1.fa"),
    ])

def test_get_heavy_MSAs():
    heavy_MSAs = get_heavy_MSAs(MSA_path, 10)
    assert sorted(heavy_MSAs) == sorted([
        str(MSA_path/"GC00000085_1.fa"),
        str(MSA_path/"GC00001725.fa"),
        str(MSA_path/"GC00001920.fa"),
        str(MSA_path/"GC00003727_5.fa"),
        str(MSA_path/"GC00005143.fa"),
    ])


test_get_locus_complexity_records()
test_get_light_MSAs()
test_get_heavy_MSAs()
