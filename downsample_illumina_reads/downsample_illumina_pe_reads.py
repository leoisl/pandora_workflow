import argparse
import random
import pysam
from abc import ABC, abstractmethod


class DownsampleIlluminaReads(ABC):
    def __init__(self, number_of_bases):
        self._number_of_bases = number_of_bases

    @property
    def number_of_bases(self):
        return self._number_of_bases

    @abstractmethod
    def _get_read_id_to_number_of_bases(self):
        pass

    @abstractmethod
    def _output_reads(self, read_pair_ids_to_output):
        pass

    def downsample_illumina_reads(self):
        read_id_to_number_of_bases = self._get_read_id_to_number_of_bases()

        number_of_reads = len(read_id_to_number_of_bases)
        random_order_of_reads = self._get_list_with_random_order_of_read_ids(number_of_reads)

        read_ids_to_output = self._get_reads_until_bases_are_saturated(read_id_to_number_of_bases, random_order_of_reads)

        self._output_reads(read_ids_to_output)


    def _get_read_id_to_number_of_bases_core (self, reads_iterator):
        read_id_to_number_of_bases = []
        for record in reads_iterator:
            number_of_bases = len(record.sequence)
            read_id_to_number_of_bases.append(number_of_bases)
        return read_id_to_number_of_bases

    def _shuffle_list(self, list):
        random.shuffle(list)
        return list

    def _get_list_with_random_order_of_read_ids(self, number_of_reads):
        list_with_random_order_of_read_ids = list(range(number_of_reads))
        return self._shuffle_list(list_with_random_order_of_read_ids)

    def _get_reads_until_bases_are_saturated(self, read_id_to_number_of_bases, random_order_of_reads):
        read_ids_to_output = set()
        number_of_bases_output = 0

        for read_index in random_order_of_reads:
            if number_of_bases_output >= self.number_of_bases:
                break
            number_of_bases_in_read = read_id_to_number_of_bases[read_index]
            number_of_bases_output += number_of_bases_in_read
            read_ids_to_output.add(read_index)

        return read_ids_to_output



class DownsampleIlluminaPEReads(DownsampleIlluminaReads):
    def __init__(self, reads1, reads2, number_of_bases, out_reads1, out_reads2):
        super().__init__(number_of_bases)
        self._reads1 = reads1
        self._reads2 = reads2
        self._out_reads1 = out_reads1
        self._out_reads2 = out_reads2

    @property
    def reads1(self):
        return self._reads1
    @property
    def reads2(self):
        return self._reads2
    @property
    def out_reads1(self):
        return self._out_reads1
    @property
    def out_reads2(self):
        return self._out_reads2

    def _get_read_id_to_number_of_bases(self):
        return self._get_read_pair_id_to_number_of_bases()

    def _get_read_pair_id_to_number_of_bases(self):
        with pysam.FastxFile(self.reads1) as reads1_fastx_file, pysam.FastxFile(self.reads2) as reads2_fastx_file:
            self._get_read_pair_id_to_number_of_bases_core(reads1_fastx_file, reads2_fastx_file)

    def _get_read_pair_id_to_number_of_bases_core(self, reads1_fastx_file, reads2_fastx_file):
        left_read_id_to_number_of_bases = self._get_read_id_to_number_of_bases_core(reads1_fastx_file)
        right_read_id_to_number_of_bases = self._get_read_id_to_number_of_bases_core(reads2_fastx_file)
        read_pair_id_to_number_of_bases = [left_bases + right_bases for left_bases, right_bases in
                                           zip(left_read_id_to_number_of_bases, right_read_id_to_number_of_bases)]
        return read_pair_id_to_number_of_bases

    def _output_reads(self, read_pair_ids_to_output):
        with pysam.FastxFile(self.reads1) as reads1_fastx_file, pysam.FastxFile(self.reads2) as reads2_fastx_file, \
             open(self.out_reads1, "w") as out_reads1_file, open(self.out_reads2, "w") as out_reads2_file:
            self._output_reads_core(read_pair_ids_to_output, reads1_fastx_file, reads2_fastx_file,
                               out_reads1_file, out_reads2_file)

    def _output_reads_core(self, read_pair_ids_to_output, reads1_fastx_file, reads2_fastx_file,
                            out_reads1_file, out_reads2_file):
        read_pair_id = 0
        for record1, record2 in zip(reads1_fastx_file, reads2_fastx_file):
            if read_pair_id in read_pair_ids_to_output:
                out_reads1_file.write(f"{str(record1)}\n")
                out_reads2_file.write(f"{str(record2)}\n")
            read_pair_id += 1


class DownsampleIlluminaSEReads(DownsampleIlluminaReads):
    def __init__(self, reads, number_of_bases, out_reads):
        super().__init__(number_of_bases)
        self._reads = reads
        self._out_reads = out_reads

    @property
    def reads(self):
        return self._reads
    @property
    def out_reads(self):
        return self._out_reads

    def _get_read_id_to_number_of_bases(self):
        with pysam.FastxFile(self.reads) as reads_fastx_file:
            return self._get_read_id_to_number_of_bases_core(reads_fastx_file)

    def _output_reads(self, read_ids_to_output):
        with pysam.FastxFile(self.reads) as reads_fastx_file, open(self.out_reads, "w") as out_reads_file:
            self._output_reads_core(read_ids_to_output, reads_fastx_file, out_reads_file)

    def _output_reads_core(self, read_ids_to_output, reads_fastx_file, out_reads_file):
        read_id = 0
        for record in reads_fastx_file:
            if read_id in read_ids_to_output:
                out_reads_file.write(f"{str(record)}\n")
            read_id += 1





# untested functions below
def get_args():
    parser = argparse.ArgumentParser(description='Downsample illumina paired-end reads.')
    parser.add_argument('--reads1')
    parser.add_argument('--reads2')
    parser.add_argument("--number_of_bases", type=int, help="Number of bases to have in the output files")
    parser.add_argument('--out_reads1')
    parser.add_argument('--out_reads2')

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = get_args()
    downsampler = DownsampleIlluminaPEReads(args.reads1, args.reads2, args.number_of_bases, args.out_reads1, args.out_reads2)
    downsampler.downsample_illumina_pe_reads()