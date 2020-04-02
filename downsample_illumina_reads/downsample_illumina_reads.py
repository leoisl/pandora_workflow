import random
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
