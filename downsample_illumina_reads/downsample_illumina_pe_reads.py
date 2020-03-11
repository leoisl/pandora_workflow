import argparse
import pysam
from downsample_illumina_reads.downsample_illumina_reads import DownsampleIlluminaReads


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
    downsampler.downsample_illumina_reads()