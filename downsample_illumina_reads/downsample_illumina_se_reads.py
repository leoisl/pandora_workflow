import argparse
import pysam
from downsample_illumina_reads.downsample_illumina_reads import DownsampleIlluminaReads


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
    parser = argparse.ArgumentParser(description='Downsample illumina single-end reads.')
    parser.add_argument('--reads')
    parser.add_argument("--number_of_bases", type=int, help="Number of bases to have in the output files")
    parser.add_argument('--out_reads')

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = get_args()
    downsampler = DownsampleIlluminaSEReads(args.reads, args.number_of_bases, args.out_reads)
    downsampler.downsample_illumina_reads()