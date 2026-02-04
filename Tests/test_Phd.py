# Revisions copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Phd module."""

import unittest
import pytest

from Bio import SeqIO
from Bio.Sequencing import Phd


class PhdTestOne(unittest.TestCase):
    def setUp(self):
        self.handle = open("Phd/phd1")

    def tearDown(self):
        self.handle.close()

    def test_check_SeqIO(self):
        """Test phd1 using parser via SeqIO."""
        records = SeqIO.parse(self.handle, "phd")
        # Contig 1
        record = next(records)
        assert record.id == "34_222_(80-A03-19).b.ab1"
        assert record.name == "34_222_(80-A03-19).b.ab1"
        assert record.description == "34_222_(80-A03-19).b.ab1"
        assert record.seq.startswith("ctccgtcggaacatcatcggatcctatcaca")
        assert record.seq.endswith("ctctcctctccctccctccgactccaaagcgtg")
        assert record.letter_annotations["phred_quality"][:10] == [9, 9, 10, 19, 22, 37, 28, 28, 24, 22]
        assert record[:10].format("fasta") == ">34_222_(80-A03-19).b.ab1\nctccgtcgga\n"
        assert record[:10].format("qual") == ">34_222_(80-A03-19).b.ab1\n9 9 10 19 22 37 28 28 24 22\n"
        assert record[:10].format("fastq") == "@34_222_(80-A03-19).b.ab1\nctccgtcgga\n+\n**+47F==97\n"
        assert record[:10].format("fastq-illumina") == "@34_222_(80-A03-19).b.ab1\nctccgtcgga\n+\nIIJSVe\\\\XV\n"
        # Contig 2
        record = next(records)
        assert record.id == "425_103_(81-A03-19).g.ab1"
        assert record.name == "425_103_(81-A03-19).g.ab1"
        assert record.letter_annotations["phred_quality"][:10] == [14, 17, 22, 10, 10, 10, 15, 8, 8, 9]
        # Contig 3
        record = next(records)
        assert record.id == "425_7_(71-A03-19).b.ab1"
        assert record.name == "425_7_(71-A03-19).b.ab1"
        assert record.letter_annotations["phred_quality"][:10] == [10, 10, 10, 10, 8, 8, 6, 6, 6, 6]
        # Make sure that no further records are found
        with pytest.raises(StopIteration):
            next(records)

    def test_check_record_parser(self):
        """Test phd1 file in detail."""
        records = Phd.parse(self.handle)
        # Record 1
        record = next(records)
        assert record.file_name == "34_222_(80-A03-19).b.ab1"
        assert record.comments["abi_thumbprint"] == 0
        assert record.comments["call_method"] == "phred"
        assert record.comments["chem"] == "term"
        assert record.comments["chromat_file"] == "34_222_(80-A03-19).b.ab1"
        assert record.comments["dye"] == "big"
        assert record.comments["phred_version"] == "0.020425.c"
        assert record.comments["quality_levels"] == 99
        assert record.comments["time"] == "Fri Feb 13 09:16:11 2004"
        assert record.comments["trace_array_max_index"] == 10867
        assert record.comments["trace_array_min_index"] == 0
        assert record.comments["trace_peak_area_ratio"] == pytest.approx(0.1467, abs=5e-8)
        assert record.comments["trim"][0] == 3
        assert record.comments["trim"][1] == 391
        assert record.comments["trim"][2] == pytest.approx(0.05, abs=5e-8)
        center = len(record.sites) // 2
        assert record.sites[0] == ("c", "9", "6")
        assert record.sites[1] == ("t", "9", "18")
        assert record.sites[2] == ("c", "10", "26")
        assert record.sites[3] == ("c", "19", "38")
        assert record.sites[4] == ("g", "22", "49")
        assert record.sites[5] == ("t", "37", "65")
        assert record.sites[6] == ("c", "28", "76")
        assert record.sites[7] == ("g", "28", "87")
        assert record.sites[8] == ("g", "24", "100")
        assert record.sites[9] == ("a", "22", "108")
        assert record.sites[center - 5] == ("c", "11", "5259")
        assert record.sites[center - 4] == ("c", "11", "5273")
        assert record.sites[center - 3] == ("t", "9", "5286")
        assert record.sites[center - 2] == ("g", "10", "5300")
        assert record.sites[center - 1] == ("a", "10", "5316")
        assert record.sites[center] == ("t", "8", "5323")
        assert record.sites[center + 1] == ("c", "8", "5343")
        assert record.sites[center + 2] == ("g", "8", "5352")
        assert record.sites[center + 3] == ("c", "8", "5366")
        assert record.sites[center + 4] == ("c", "8", "5378")
        assert record.sites[-10] == ("c", "8", "10756")
        assert record.sites[-9] == ("c", "8", "10764")
        assert record.sites[-8] == ("a", "8", "10769")
        assert record.sites[-7] == ("a", "8", "10788")
        assert record.sites[-6] == ("a", "8", "10803")
        assert record.sites[-5] == ("g", "10", "10816")
        assert record.sites[-4] == ("c", "11", "10826")
        assert record.sites[-3] == ("g", "11", "10840")
        assert record.sites[-2] == ("t", "11", "10855")
        assert record.sites[-1] == ("g", "11", "10864")
        assert record.seq[:10] == "ctccgtcgga"
        assert record.seq[-10:] == "ccaaagcgtg"
        assert record.seq_trimmed[:10] == "cgtcggaaca"
        assert record.seq_trimmed[-10:] == "tatttcggag"
        # Record 2
        record = next(records)
        center = len(record.sites) // 2
        assert record.file_name == "425_103_(81-A03-19).g.ab1"
        assert record.comments["abi_thumbprint"] == 0
        assert record.comments["call_method"] == "phred"
        assert record.comments["chem"] == "term"
        assert record.comments["chromat_file"] == "425_103_(81-A03-19).g.ab1"
        assert record.comments["dye"] == "big"
        assert record.comments["phred_version"] == "0.020425.c"
        assert record.comments["quality_levels"] == 99
        assert record.comments["time"] == "Tue Feb 17 10:31:15 2004"
        assert record.comments["trace_array_max_index"] == 10606
        assert record.comments["trace_array_min_index"] == 0
        assert record.comments["trace_peak_area_ratio"] == pytest.approx(0.0226, abs=5e-8)
        assert record.comments["trim"][0] == 10
        assert record.comments["trim"][1] == 432
        assert record.comments["trim"][2] == pytest.approx(0.05, abs=5e-8)
        assert record.sites[0] == ("c", "14", "3")
        assert record.sites[1] == ("g", "17", "11")
        assert record.sites[2] == ("g", "22", "23")
        assert record.sites[3] == ("g", "10", "35")
        assert record.sites[4] == ("a", "10", "53")
        assert record.sites[5] == ("t", "10", "68")
        assert record.sites[6] == ("c", "15", "75")
        assert record.sites[7] == ("c", "8", "85")
        assert record.sites[8] == ("c", "8", "94")
        assert record.sites[9] == ("a", "9", "115")
        assert record.sites[center - 5] == ("c", "33", "5140")
        assert record.sites[center - 4] == ("c", "28", "5156")
        assert record.sites[center - 3] == ("g", "25", "5167")
        assert record.sites[center - 2] == ("c", "28", "5178")
        assert record.sites[center - 1] == ("c", "18", "5193")
        assert record.sites[center] == ("a", "16", "5204")
        assert record.sites[center + 1] == ("a", "15", "5213")
        assert record.sites[center + 2] == ("a", "10", "5230")
        assert record.sites[center + 3] == ("a", "10", "5242")
        assert record.sites[center + 4] == ("t", "8", "5249")
        assert record.sites[-10] == ("c", "8", "10489")
        assert record.sites[-9] == ("c", "8", "10503")
        assert record.sites[-8] == ("c", "8", "10514")
        assert record.sites[-7] == ("a", "8", "10516")
        assert record.sites[-6] == ("g", "8", "10530")
        assert record.sites[-5] == ("c", "8", "10550")
        assert record.sites[-4] == ("c", "10", "10566")
        assert record.sites[-3] == ("a", "8", "10574")
        assert record.sites[-2] == ("a", "7", "10584")
        assert record.sites[-1] == ("g", "7", "10599")
        assert record.seq[:10] == "cgggatccca"
        assert record.seq[-10:] == "cccagccaag"
        assert record.seq_trimmed[:10] == "cctgatccga"
        assert record.seq_trimmed[-10:] == "ggggccgcca"
        # Record 3
        record = next(records)
        center = len(record.sites) // 2
        assert record.file_name == "425_7_(71-A03-19).b.ab1"
        assert record.comments["abi_thumbprint"] == 0
        assert record.comments["call_method"] == "phred"
        assert record.comments["chem"] == "term"
        assert record.comments["chromat_file"] == "425_7_(71-A03-19).b.ab1"
        assert record.comments["dye"] == "big"
        assert record.comments["phred_version"] == "0.020425.c"
        assert record.comments["quality_levels"] == 99
        assert record.comments["time"] == "Thu Jan 29 11:46:14 2004"
        assert record.comments["trace_array_max_index"] == 9513
        assert record.comments["trace_array_min_index"] == 0
        assert record.comments["trace_peak_area_ratio"] == pytest.approx(100.0, abs=5e-8)
        assert record.comments["trim"][0] == -1
        assert record.comments["trim"][1] == -1
        assert record.comments["trim"][2] == 0.05
        assert record.sites[0] == ("a", "10", "7")
        assert record.sites[1] == ("c", "10", "13")
        assert record.sites[2] == ("a", "10", "21")
        assert record.sites[3] == ("t", "10", "28")
        assert record.sites[4] == ("a", "8", "33")
        assert record.sites[5] == ("a", "8", "40")
        assert record.sites[6] == ("a", "6", "50")
        assert record.sites[7] == ("t", "6", "53")
        assert record.sites[8] == ("c", "6", "66")
        assert record.sites[9] == ("a", "6", "68")
        assert record.sites[center - 5] == ("a", "6", "4728")
        assert record.sites[center - 4] == ("t", "10", "4737")
        assert record.sites[center - 3] == ("a", "10", "4746")
        assert record.sites[center - 2] == ("a", "8", "4756")
        assert record.sites[center - 1] == ("t", "8", "4759")
        assert record.sites[center] == ("t", "8", "4768")
        assert record.sites[center + 1] == ("a", "8", "4775")
        assert record.sites[center + 2] == ("g", "10", "4783")
        assert record.sites[center + 3] == ("t", "8", "4788")
        assert record.sites[center + 4] == ("g", "8", "4794")
        assert record.sites[-10] == ("a", "8", "9445")
        assert record.sites[-9] == ("t", "6", "9453")
        assert record.sites[-8] == ("c", "6", "9462")
        assert record.sites[-7] == ("t", "6", "9465")
        assert record.sites[-6] == ("g", "6", "9478")
        assert record.sites[-5] == ("c", "6", "9483")
        assert record.sites[-4] == ("t", "6", "9485")
        assert record.sites[-3] == ("t", "8", "9495")
        assert record.sites[-2] == ("t", "3", "9504")
        assert record.sites[-1] == ("n", "0", "9511")
        assert record.seq[:10] == "acataaatca"
        assert record.seq[-10:] == "atctgctttn"
        # Make sure that no further records are found
        with pytest.raises(StopIteration):
            next(records)


class PhdTestTwo(unittest.TestCase):
    def setUp(self):
        self.handle = open("Phd/phd2")

    def tearDown(self):
        self.handle.close()

    def test_check_SeqIO(self):
        """Test phd2 using parser via SeqIO."""
        records = SeqIO.parse(self.handle, "phd")
        # Contig 1
        record = next(records)
        assert record.id == "ML4924R"
        assert record.name == "ML4924R"
        assert record.description == "ML4924R"
        assert record.seq.startswith("actttggtcgcctgcaggtaccggtccgnga")
        assert record.seq.endswith("agaagctcgttctcaacatctccgttggtgaga")
        assert record.letter_annotations["phred_quality"][:10] == [6, 6, 6, 8, 8, 12, 18, 16, 14, 11]
        assert record[:10].format("fasta") == ">ML4924R\nactttggtcg\n"
        assert record[:10].format("qual") == ">ML4924R\n6 6 6 8 8 12 18 16 14 11\n"
        assert record[:10].format("fastq") == "@ML4924R\nactttggtcg\n+\n'''))-31/,\n"
        assert record[:10].format("fastq-illumina") == "@ML4924R\nactttggtcg\n+\nFFFHHLRPNK\n"
        # Make sure that no further records are found
        with pytest.raises(StopIteration):
            next(records)


class PhdTest454(unittest.TestCase):
    def setUp(self):
        self.handle = open("Phd/phd_454")

    def tearDown(self):
        self.handle.close()

    def test_check_SeqIO(self):
        """Test phd_454 using parser via SeqIO."""
        records = SeqIO.parse(self.handle, "phd")
        # Contig 1
        record = next(records)
        assert record.id == "EBE03TV04IHLTF.77-243"
        assert record.name == "EBE03TV04IHLTF.77-243"
        assert record.description == "EBE03TV04IHLTF.77-243 1"
        assert record.seq == "ggggatgaaagggatctcggtggtaggtga"
        assert record.letter_annotations["phred_quality"][:10] == [37, 37, 37, 37, 37, 37, 37, 37, 37, 37]
        assert record.format("fasta") == ">EBE03TV04IHLTF.77-243 1\nggggatgaaagggatctcggtggtaggtga\n"
        assert (record.format("qual") == ">EBE03TV04IHLTF.77-243 1\n"
            "37 37 37 37 37 37 37 37 37 37 "
            "37 37 37 26 26 26 30 33 33 33\n"
            "33 33 36 36 33 33 33 36 26 22\n")
        assert (record.format("fastq") == "@EBE03TV04IHLTF.77-243 1\n"
            "ggggatgaaagggatctcggtggtaggtga\n"
            "+\n"
            "FFFFFFFFFFFFF;;;?BBBBBEEBBBE;7\n")
        assert record[:10].format("fastq-illumina") == "@EBE03TV04IHLTF.77-243 1\nggggatgaaa\n+\neeeeeeeeee\n"
        # Make sure that no further records are found
        with pytest.raises(StopIteration):
            next(records)


class PhdTestSolexa(unittest.TestCase):
    def setUp(self):
        self.handle = open("Phd/phd_solexa")

    def tearDown(self):
        self.handle.close()

    def test_check_SeqIO(self):
        """Test phd2 using parser via SeqIO."""
        records = SeqIO.parse(self.handle, "phd")
        # Contig 1
        record = next(records)
        assert record.id == "HWI-EAS94_4_1_1_537_446"
        assert record.name == "HWI-EAS94_4_1_1_537_446"
        assert record.description == "HWI-EAS94_4_1_1_537_446 1"
        assert record.seq == "gccaatcaggtttctctgcaagcccctttagcagctgagc"
        assert record.letter_annotations["phred_quality"] == [
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                28,
                23,
                30,
                30,
                30,
                30,
                30,
                30,
                28,
                22,
                8,
                22,
                7,
                15,
                15,
                15,
                10,
                10,
                11,
                15,
            ]
        assert record.format("fasta") == ">HWI-EAS94_4_1_1_537_446 1\ngccaatcaggtttctctgcaagcccctttagcagctgagc\n"
        assert (record.format("qual") == ">HWI-EAS94_4_1_1_537_446 1\n"
            "30 30 30 30 30 30 30 30 30 30 "
            "30 30 30 30 30 30 30 30 30 30\n"
            "28 23 30 30 30 30 30 30 28 22 "
            "8 22 7 15 15 15 10 10 11 15\n")
        assert (record.format("fastq") == "@HWI-EAS94_4_1_1_537_446 1\n"
            "gccaatcaggtttctctgcaagcccctttagcagctgagc\n"
            "+\n"
            "????????????????????=8??????=7)7(000++,0\n")
        assert (record.format("fastq-illumina") == "@HWI-EAS94_4_1_1_537_446 1\n"
            "gccaatcaggtttctctgcaagcccctttagcagctgagc\n"
            "+\n"
            "^^^^^^^^^^^^^^^^^^^^\\W^^^^^^\\VHVGOOOJJKO\n")
        # Contig 2
        record = next(records)
        assert record.id == "HWI-EAS94_4_1_1_602_99"
        assert record.name == "HWI-EAS94_4_1_1_602_99"
        assert record.description == "HWI-EAS94_4_1_1_602_99 1"
        assert record.seq == "gccatggcacatatatgaaggtcagaggacaacttgctgt"
        assert record.letter_annotations["phred_quality"] == [
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                30,
                16,
                30,
                28,
                22,
                22,
                22,
                14,
                15,
                15,
                5,
                10,
                15,
                10,
                5,
            ]
        assert record.format("fasta") == ">HWI-EAS94_4_1_1_602_99 1\ngccatggcacatatatgaaggtcagaggacaacttgctgt\n"
        assert (record.format("qual") == ">HWI-EAS94_4_1_1_602_99 1\n"
            "30 30 30 30 30 30 30 30 30 30 "
            "30 30 30 30 30 30 30 30 30 30\n"
            "30 30 30 30 30 30 16 30 28 22 "
            "22 22 14 15 15 5 10 15 10 5\n")
        assert (record.format("fastq") == "@HWI-EAS94_4_1_1_602_99 1\n"
            "gccatggcacatatatgaaggtcagaggacaacttgctgt\n"
            "+\n"
            "??????????????????????????1?=777/00&+0+&\n")
        assert (record.format("fastq-illumina") == "@HWI-EAS94_4_1_1_602_99 1\n"
            "gccatggcacatatatgaaggtcagaggacaacttgctgt\n"
            "+\n"
            "^^^^^^^^^^^^^^^^^^^^^^^^^^P^\\VVVNOOEJOJE\n")
        # Make sure that no further records are found
        with pytest.raises(StopIteration):
            next(records)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
