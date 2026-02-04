"""Tests for SeqIO GFA module."""

import unittest
import pytest

from Bio import BiopythonWarning
from Bio import SeqIO


class TestRead(unittest.TestCase):
    def test_read_GFA1(self):
        """Test parsing valid GFA 1.x files."""
        records = list(SeqIO.parse("GFA/seq.gfa", "gfa1"))
        assert len(records) == 8
        assert records[6].id == "MTh13014"
        assert records[6].seq == "TTAGGTCTCCACCCCTGACTCCCCTCAGCCATAGAAGGCCCCACCCCAGTCTCAGCCCTACTCCACTCAAGCACTATAGTTGTAGCAGGAATCTTCTTACTCATCCGCTTCCACCCCCTAGCAGAAAATAGCCCACTAATCCAAACTCTAACACTATGCTTAGGCGCTATCACCACTCTGTTCGCAGCAGTCTGCGCCCTTACACAAAATGACATCAAAAAAATCGTAGCCTTCTCCACTTCAAGTCAACTAGGACTCATAATAGTTACAATCGGCATCAACCAACCACACCTAGCATTCCTGCACATCTGTACCCACGCCTTCTTCAAAGCCATACTATTTATGTGCTCCGGGTCCATCATCCACAACCTTAACAATGAACAAGATATTCGAAAAATAGGAGGACTACTCAAAACCATACCTCTCACTTCAACCTCCCTCACCATTGGCAGCCTAGCATTAGCAGGAATACCTTTCCTCACAGGTTTCTACTCCAAAGACC"
        assert records[0].annotations["SN"] == ("Z", "MT_human")
        assert records[0].annotations["SO"] == ("i", "0")

        records = list(SeqIO.parse("GFA/seq_with_len.gfa", "gfa1"))
        assert len(records) == 9
        assert records[8].seq == "GAAAAATTGCCCTTGGTTTTCGCTTCGCTCAAACTCTATTGAACTTCGCTTTCGCTCAGTTCGTCGGGGCAATTTTTTGGTTAATACTT"

        records = list(SeqIO.parse("GFA/fake_with_checksum.gfa", "gfa1"))
        assert len(records) == 1
        assert records[0].seq == "AAA"

        records = list(SeqIO.parse("GFA/no_seq.gfa", "gfa1"))
        assert len(records) == 9
        assert len(records[0]) == 528

    def test_read_GFA2(self):
        """Test parsing valid GFA 2.0 files."""
        records = list(SeqIO.parse("GFA/fake_gfa2.gfa", "gfa2"))
        assert len(records) == 1
        assert records[0].seq == "AAA"


class TestCorrupt(unittest.TestCase):
    def test_corrupt_gfa2(self):
        """Check a GFA 1.x file does not parse in GFA 2."""
        with pytest.raises(ValueError):
            list(SeqIO.parse("GFA/seq.gfa", "gfa2"))

    def test_corrupt_segment_fields(self):
        """Check a GFA file with invalid fields on a segment line."""
        with pytest.raises(ValueError):
            list(SeqIO.parse("GFA/corrupt_segment_fields.gfa", "gfa1"))

    def test_corrupt_len(self):
        """Check a GFA file with an incorrect length."""
        with pytest.warns(BiopythonWarning):
            list(SeqIO.parse("GFA/corrupt_len.gfa", "gfa1"))

    def test_corrupt_checksum(self):
        """Check a GFA file with an incorrect checksum."""
        with pytest.warns(BiopythonWarning):
            list(SeqIO.parse("GFA/corrupt_checksum.gfa", "gfa1"))

    def test_corrupt_tag_name(self):
        """Check a GFA file with an invalid tag name."""
        with pytest.warns(BiopythonWarning):
            list(SeqIO.parse("GFA/corrupt_tag_name.gfa", "gfa1"))

    def test_corrupt_tag_type(self):
        """Check a GFA file with an incorrect tag type."""
        with pytest.warns(BiopythonWarning):
            list(SeqIO.parse("GFA/corrupt_tag_type.gfa", "gfa1"))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
