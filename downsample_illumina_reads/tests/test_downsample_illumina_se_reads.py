from downsample_illumina_reads.downsample_illumina_se_reads import DownsampleIlluminaSEReads
from unittest.mock import Mock, patch, call
from unittest import TestCase


class TestDownsampleIlluminaSEReads(TestCase):
    def setUp(self) -> None:
        self.dummy_downsampler = DownsampleIlluminaSEReads(None, None, None)
        self.reads_fastx_file = ["read_1", "read_2"]
        self.out_reads_file_write_mock = Mock()
        self.out_reads_file = Mock(write=self.out_reads_file_write_mock)

    def test_____init__(self):
        downsampler = DownsampleIlluminaSEReads("reads", "number_of_bases", "out_reads")
        self.assertEqual(downsampler.reads, "reads")
        self.assertEqual(downsampler.number_of_bases, "number_of_bases")
        self.assertEqual(downsampler.out_reads, "out_reads")

    def test___output_reads_core___read_pair_ids_to_output_is_empty___no_reads_are_output(self):
        read_pair_ids_to_output = []
        self.dummy_downsampler._output_reads_core(read_pair_ids_to_output, self.reads_fastx_file, self.out_reads_file)
        self.out_reads_file_write_mock.assert_not_called()


    def test___output_reads_core___read_pair_ids_contains_one_read___one_read_is_output(self):
        read_pair_ids_to_output = [1]
        self.dummy_downsampler._output_reads_core(read_pair_ids_to_output, self.reads_fastx_file, self.out_reads_file)
        self.out_reads_file_write_mock.assert_called_once_with("read_2\n")

    def test___output_reads_core___read_pair_ids_contains_two_reads___two_reads_are_output(self):
        read_pair_ids_to_output = [0, 1]
        self.dummy_downsampler._output_reads_core(read_pair_ids_to_output, self.reads_fastx_file, self.out_reads_file)
        self.assertEqual(self.out_reads_file_write_mock.call_count, 2)
        self.out_reads_file_write_mock.has_calls(call("read_1\n"), call("read_2\n"), any_order=False)
