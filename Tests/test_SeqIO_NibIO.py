"""Tests for SeqIO NibIO module."""

import unittest
import pytest
from io import BytesIO

from Bio import SeqIO


class TestNibReaderWriter(unittest.TestCase):
    def test_read_even(self):
        with open("Nib/test_even.fa") as handle:
            record = SeqIO.read(handle, "fasta")
        sequence = record.seq
        with open("Nib/test_even_bigendian.nib", "rb") as handle:
            record = SeqIO.read(handle, "nib")
        assert sequence == record.seq
        with open("Nib/test_even_littleendian.nib", "rb") as handle:
            record = SeqIO.read(handle, "nib")
        assert sequence == record.seq

    def test_read_odd(self):
        with open("Nib/test_odd.fa") as handle:
            record = SeqIO.read(handle, "fasta")
        sequence = record.seq
        with open("Nib/test_odd_bigendian.nib", "rb") as handle:
            record = SeqIO.read(handle, "nib")
        assert sequence == record.seq
        with open("Nib/test_odd_littleendian.nib", "rb") as handle:
            record = SeqIO.read(handle, "nib")
        assert sequence == record.seq

    def test_write_even(self):
        with open("Nib/test_even.fa") as handle:
            record = SeqIO.read(handle, "fasta")
        sequence = record.seq
        handle = BytesIO()
        n = SeqIO.write(record, handle, "nib")
        assert n == 1
        handle.flush()
        handle.seek(0)
        record = SeqIO.read(handle, "nib")
        handle.close()
        assert sequence == record.seq

    def test_write_odd(self):
        with open("Nib/test_odd.fa") as handle:
            record = SeqIO.read(handle, "fasta")
        sequence = record.seq
        handle = BytesIO()
        n = SeqIO.write(record, handle, "nib")
        assert n == 1
        handle.flush()
        handle.seek(0)
        record = SeqIO.read(handle, "nib")
        handle.close()
        assert sequence == record.seq


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
