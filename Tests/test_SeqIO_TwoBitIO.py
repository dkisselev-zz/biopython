"""Tests for SeqIO TwoBitIO module."""

import unittest
import pytest

from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
from Bio.Seq import UndefinedSequenceError
from Bio.SeqRecord import SeqRecord


class Parsing(unittest.TestCase):
    """Test parsing 2bit files."""

    def setUp(self):
        path = "TwoBit/sequence.fa"
        records = SeqIO.parse(path, "fasta")
        self.records = list(records)

    def test_littleendian(self, step=5):
        path = "TwoBit/sequence.littleendian.2bit"
        with open(path, "rb") as stream:
            records = SeqIO.parse(stream, "twobit")
            assert records.byteorder == "little"
            assert len(self.records) == len(records)
            for record1, record2 in zip(self.records, records):
                assert record1.id == record2.id
                seq1 = record1.seq
                seq2 = record2.seq
                assert seq1 == seq2
                n = len(seq1)
                for i in range(0, n, step):
                    for j in range(i, n, step):
                        assert seq1[i:j] == seq2[i:j]
                        assert repr(seq1[i:j]) == repr(seq2[i:j])

    def test_bigendian(self, step=5):
        path = "TwoBit/sequence.bigendian.2bit"
        with open(path, "rb") as stream:
            records = SeqIO.parse(stream, "twobit")
            assert len(records) == 6
            assert records.byteorder == "big"
            for record1, record2 in zip(self.records, records):
                assert record1.id == record2.id
                seq1 = record1.seq
                seq2 = record2.seq
                assert seq1 == seq2
                n = len(seq1)
                for i in range(0, n, step):
                    for j in range(i, n, step):
                        assert seq1[i:j] == seq2[i:j]
                        assert repr(seq1[i:j]) == repr(seq2[i:j])

    def test_sequence_long(self):
        path = "TwoBit/sequence.long.2bit"
        with open(path, "rb") as stream:
            with pytest.raises(ValueError) as cm:
                SeqIO.parse(stream, "twobit")
            assert str(cm.value) == "version-1 twoBit files with 64-bit offsets for index are currently not supported"


class TestComparisons(unittest.TestCase):
    """Test comparisons of sequences read from 2bit files to Seq and other objects."""

    def setUp(self):
        path = "TwoBit/sequence.bigendian.2bit"
        self.stream = open(path, "rb")
        records = SeqIO.parse(self.stream, "twobit")
        record1 = next(records)
        record2 = next(records)
        self.seq1a = record1.seq
        self.seq2a = record2.seq
        path = "TwoBit/sequence.fa"
        records = SeqIO.parse(path, "fasta")
        record1 = next(records)
        record2 = next(records)
        self.seq1b = record1.seq
        self.seq2b = record2.seq

    def tearDown(self):
        self.stream.close()

    def test_eq(self):
        seq1a = self.seq1a
        seq2a = self.seq2a
        seq1b = self.seq1b
        seq2b = self.seq2b
        assert seq1a == seq1b
        assert seq2a == seq2b
        assert seq1a == seq1a
        assert seq2a == seq2a
        with pytest.raises(UndefinedSequenceError):
            seq1a == Seq(None, len(seq1a))
        with pytest.raises(UndefinedSequenceError):
            seq2a == Seq(None, len(seq2a))
        with pytest.raises(UndefinedSequenceError):
            seq1a == Seq(None, 10)
        with pytest.raises(UndefinedSequenceError):
            seq2a == Seq(None, 10)
        with pytest.raises(UndefinedSequenceError):
            Seq(None, len(seq1a)) == seq1a
        with pytest.raises(UndefinedSequenceError):
            Seq(None, len(seq2a)) == seq2a
        with pytest.raises(UndefinedSequenceError):
            Seq(None, 10) == seq1a
        with pytest.raises(UndefinedSequenceError):
            Seq(None, 10) == seq2a

    def test_ne(self):
        seq1a = self.seq1a
        seq2a = self.seq2a
        seq1b = self.seq1b
        seq2b = self.seq2b
        assert seq1a != seq2a
        assert seq1a != seq2b
        with pytest.raises(UndefinedSequenceError):
            seq1a != Seq(None, len(seq1a))
        with pytest.raises(UndefinedSequenceError):
            seq2a != Seq(None, len(seq2a))
        with pytest.raises(UndefinedSequenceError):
            seq1a != Seq(None, 10)
        with pytest.raises(UndefinedSequenceError):
            seq2a != Seq(None, 10)
        with pytest.raises(UndefinedSequenceError):
            Seq(None, len(seq1a)) != seq1a
        with pytest.raises(UndefinedSequenceError):
            Seq(None, len(seq2a)) != seq2a
        with pytest.raises(UndefinedSequenceError):
            Seq(None, 10) != seq1a
        with pytest.raises(UndefinedSequenceError):
            Seq(None, 10) != seq2a

    def test_lt(self):
        seq1 = self.seq1a
        seq2 = self.seq2a
        assert seq1 < seq2
        assert "AA" < seq1
        assert seq1 < "TT"
        assert "AA" < seq2
        assert seq2 < "TTT"
        assert b"AA" < seq1
        assert seq1 < b"TT"
        assert b"AA" < seq2
        assert seq2 < b"TTT"
        assert Seq("AA") < seq1
        assert seq1 < Seq("TT")
        assert Seq("AA") < seq2
        assert seq2 < Seq("TTT")
        assert MutableSeq("AA") < seq1
        assert seq1 < MutableSeq("TT")
        assert MutableSeq("AA") < seq2
        assert seq2 < MutableSeq("TTT")
        with pytest.raises(UndefinedSequenceError):
            seq1 < Seq(None, len(seq1))
        with pytest.raises(UndefinedSequenceError):
            seq2 < Seq(None, len(seq2))
        with pytest.raises(UndefinedSequenceError):
            seq1 < Seq(None, 10)
        with pytest.raises(UndefinedSequenceError):
            seq2 < Seq(None, 10)

    def test_le(self):
        seq1 = self.seq1a
        seq2 = self.seq2a
        assert seq1 <= seq2
        assert seq1 <= "TT"
        assert "TT" <= seq2
        assert seq1 <= b"TT"
        assert "TT" <= seq2
        assert seq1 <= Seq("TT")
        assert "TT" <= seq2
        assert seq1 <= MutableSeq("TT")
        assert MutableSeq("TT") <= seq2
        assert "AA" <= seq1
        assert "AA" <= seq2
        assert b"AA" <= seq1
        assert b"AA" <= seq2
        assert Seq("AA") <= seq1
        assert Seq("AA") <= seq2
        assert MutableSeq("AA") <= seq1
        assert MutableSeq("AA") <= seq2
        assert "GC" <= seq1
        assert "GC" <= seq2
        assert b"GC" <= seq1
        assert b"GC" <= seq2
        assert Seq("GC") <= seq1
        assert Seq("GC") <= seq2
        assert MutableSeq("GC") <= seq1
        assert MutableSeq("GC") <= seq2
        with pytest.raises(UndefinedSequenceError):
            seq1 <= Seq(None, len(seq1))
        with pytest.raises(UndefinedSequenceError):
            seq2 <= Seq(None, len(seq2))
        with pytest.raises(UndefinedSequenceError):
            seq1 <= Seq(None, 10)
        with pytest.raises(UndefinedSequenceError):
            seq2 <= Seq(None, 10)

    def test_gt(self):
        seq1 = self.seq1a
        seq2 = self.seq2a
        assert seq2 > seq1
        assert "TT" > seq1
        assert seq2 > "TT"
        assert b"TT" > seq1
        assert seq2 > b"TT"
        assert Seq("TT") > seq1
        assert seq2 > Seq("TT")
        assert MutableSeq("TT") > seq1
        assert seq2 > MutableSeq("TT")
        assert seq1 > "AA"
        assert seq2 > "AA"
        assert seq1 > b"AA"
        assert seq2 > b"AA"
        assert seq1 > Seq("AA")
        assert seq2 > Seq("AA")
        assert seq1 > MutableSeq("AA")
        assert seq2 > MutableSeq("AA")
        assert seq1 > "GC"
        assert seq2 > "GC"
        assert seq1 > b"GC"
        assert seq2 > b"GC"
        assert seq1 > Seq("GC")
        assert seq2 > Seq("GC")
        assert seq1 > MutableSeq("GC")
        assert seq2 > MutableSeq("GC")
        with pytest.raises(UndefinedSequenceError):
            seq1 > Seq(None, len(seq1))
        with pytest.raises(UndefinedSequenceError):
            seq2 > Seq(None, len(seq2))
        with pytest.raises(UndefinedSequenceError):
            seq1 > Seq(None, 10)
        with pytest.raises(UndefinedSequenceError):
            seq2 > Seq(None, 10)

    def test_ge(self):
        seq1 = self.seq1a
        seq2 = self.seq2a
        assert seq2 >= seq1
        assert "TT" >= seq1
        assert seq2 >= "TT"
        assert b"TT" >= seq1
        assert seq2 >= b"TT"
        assert Seq("TT") >= seq1
        assert seq2 >= Seq("TT")
        assert MutableSeq("TT") >= seq1
        assert seq2 >= MutableSeq("TT")
        assert seq1 >= "AA"
        assert seq2 >= "AA"
        assert seq1 >= b"AA"
        assert seq2 >= b"AA"
        assert seq1 >= Seq("AA")
        assert seq2 >= Seq("AA")
        assert seq1 >= MutableSeq("AA")
        assert seq2 >= MutableSeq("AA")
        assert seq1 >= "GC"
        assert seq2 >= "GC"
        assert seq1 >= b"GC"
        assert seq2 >= b"GC"
        assert seq1 >= Seq("GC")
        assert seq2 >= Seq("GC")
        assert seq1 >= MutableSeq("GC")
        assert seq2 >= MutableSeq("GC")
        with pytest.raises(UndefinedSequenceError):
            seq1 >= Seq(None, len(seq1))
        with pytest.raises(UndefinedSequenceError):
            seq2 >= Seq(None, len(seq2))
        with pytest.raises(UndefinedSequenceError):
            seq1 >= Seq(None, 10)
        with pytest.raises(UndefinedSequenceError):
            seq2 >= Seq(None, 10)


class TestBaseClassMethods(unittest.TestCase):
    """Test if methods from the base class are called correctly."""

    def setUp(self):
        path = "TwoBit/sequence.bigendian.2bit"
        self.stream = open(path, "rb")
        records = SeqIO.parse(self.stream, "twobit")
        self.record1_twobit = next(records)
        self.seq1_twobit = self.record1_twobit.seq
        self.record2_twobit = next(records)
        self.seq2_twobit = self.record2_twobit.seq
        path = "TwoBit/sequence.fa"
        records = SeqIO.parse(path, "fasta")
        self.record1_fasta = next(records)
        self.seq1_fasta = self.record1_fasta.seq
        self.record2_fasta = next(records)
        self.seq2_fasta = self.record2_fasta.seq

    def tearDown(self):
        self.stream.close()

    def test_getitem(self):
        assert self.seq1_twobit == self.seq1_fasta
        assert self.seq2_twobit == self.seq2_fasta
        assert self.seq1_twobit[:] == self.seq1_fasta[:]
        assert self.seq2_twobit[:] == self.seq2_fasta[:]
        assert self.seq1_twobit[30] == self.seq1_fasta[30]
        assert self.seq2_twobit[30] == self.seq2_fasta[30]
        assert self.seq1_twobit[-30] == self.seq1_fasta[-30]
        assert self.seq2_twobit[-30] == self.seq2_fasta[-30]
        assert self.record1_twobit.seq == self.record1_fasta.seq
        assert self.record2_twobit.seq == self.record2_fasta.seq
        assert self.record1_twobit[:].seq == self.record1_fasta[:].seq
        assert self.record2_twobit[:].seq == self.record2_fasta[:].seq
        assert self.record1_twobit[30] == self.record1_fasta[30]
        assert self.record2_twobit[30] == self.record2_fasta[30]
        assert self.record1_twobit[-30] == self.record1_fasta[-30]
        assert self.record2_twobit[-30] == self.record2_fasta[-30]

    def test_bytes(self):
        b = bytes(self.seq1_twobit)
        assert isinstance(b, bytes)
        assert len(b) == 480
        assert b == bytes(self.seq1_fasta)
        b = bytes(self.seq1_twobit[:10])
        assert len(b) == 10
        assert isinstance(b, bytes)
        assert b == b"GTATACCCCT"

    def test_hash(self):
        assert hash(self.seq1_twobit) == hash(self.seq1_fasta)

    def test_add(self):
        assert isinstance(self.seq1_twobit + "ABCD", Seq)
        assert self.seq1_twobit + "ABCD" == self.seq1_fasta + "ABCD"
        assert isinstance(self.record1_twobit + "ABCD", SeqRecord)
        record1_twobit = self.record1_twobit + "ABCD"
        record1_fasta = self.record1_fasta + "ABCD"
        assert self.record1_twobit.seq == self.record1_fasta.seq

    def test_radd(self):
        assert isinstance("ABCD" + self.seq1_twobit, Seq)
        assert "ABCD" + self.seq1_twobit == "ABCD" + self.seq1_fasta
        assert isinstance("ABCD" + self.record1_twobit, SeqRecord)
        record1_twobit = "ABCD" + self.record1_twobit
        record1_fasta = "ABCD" + self.record1_fasta
        assert self.record1_twobit.seq == self.record1_fasta.seq

    def test_mul(self):
        assert isinstance(2 * self.seq1_twobit, Seq)
        assert 2 * self.seq1_twobit == 2 * self.seq1_fasta
        assert isinstance(self.seq1_twobit * 2, Seq)
        assert self.seq1_twobit * 2 == self.seq1_fasta * 2

    def test_contains(self):
        for seq in (
            self.seq1_twobit,
            self.seq1_fasta,
            self.record1_twobit,
            self.record1_fasta,
        ):
            assert "ACCCCT" in seq
            assert "ACGTACGT" not in seq

    def test_repr(self):
        assert isinstance(repr(self.seq1_twobit), str)
        assert repr(self.seq1_twobit) == repr(self.seq1_fasta)

    def test_str(self):
        assert isinstance(str(self.seq1_twobit), str)
        assert str(self.seq1_twobit) == str(self.seq1_fasta)

    def test_count(self):
        assert self.seq1_twobit.count("CT") == self.seq1_fasta.count("CT")
        assert self.seq1_twobit.count("CT", 75) == self.seq1_fasta.count("CT", 75)
        assert self.seq1_twobit.count("CT", 125, 250) == self.seq1_fasta.count("CT", 125, 250)
        assert self.record1_twobit.count("CT") == self.record1_fasta.count("CT")
        assert self.record1_twobit.count("CT", 75) == self.record1_fasta.count("CT", 75)
        assert self.record1_twobit.count("CT", 125, 250) == self.record1_fasta.count("CT", 125, 250)

    def test_find(self):
        assert self.seq1_twobit.find("CT") == self.seq1_fasta.find("CT")
        assert self.seq1_twobit.find("CT", 75) == self.seq1_fasta.find("CT", 75)
        assert self.seq1_twobit.find("CT", 75, 100) == self.seq1_fasta.find("CT", 75, 100)
        assert self.seq1_twobit.find("CT", None, 100) == self.seq1_fasta.find("CT", None, 100)

    def test_rfind(self):
        assert self.seq1_twobit.rfind("CT") == self.seq1_fasta.rfind("CT")
        assert self.seq1_twobit.rfind("CT", 450) == self.seq1_fasta.rfind("CT", 450)
        assert self.seq1_twobit.rfind("CT", None, 100) == self.seq1_fasta.rfind("CT", None, 100)
        assert self.seq1_twobit.rfind("CT", 75, 100) == self.seq1_fasta.rfind("CT", 75, 100)

    def test_index(self):
        assert self.seq1_twobit.index("CT") == self.seq1_fasta.index("CT")
        assert self.seq1_twobit.index("CT", 75) == self.seq1_fasta.index("CT", 75)
        assert self.seq1_twobit.index("CT", None, 100) == self.seq1_fasta.index("CT", None, 100)
        for seq in (self.seq1_twobit, self.seq1_fasta):
            with pytest.raises(ValueError):
                seq.index("CT", 75, 100)

    def test_rindex(self):
        assert self.seq1_twobit.rindex("CT") == self.seq1_fasta.rindex("CT")
        assert self.seq1_twobit.rindex("CT", None, 100) == self.seq1_fasta.rindex("CT", None, 100)
        for seq in (self.seq1_twobit, self.seq1_fasta):
            with pytest.raises(ValueError):
                seq.rindex("CT", 450)
            with pytest.raises(ValueError):
                seq.rindex("CT", 75, 100)

    def test_startswith(self):
        for seq in (self.seq1_twobit, self.seq1_fasta):
            assert seq.startswith("GTAT")
            assert seq.startswith("TGGG", start=10)
            assert seq.startswith("TGGG", start=10, end=14)
            assert not seq.startswith("TGGG", start=10, end=12)

    def test_endswith(self):
        for seq in (self.seq1_twobit, self.seq1_fasta):
            assert seq.endswith("ACCG")
            assert seq.endswith("ACCG", 476)
            assert seq.endswith("GCAC", 472, 478)
            assert not seq.endswith("GCAC", 476, 478)

    def test_split(self):
        assert self.seq1_twobit.split() == self.seq1_fasta.split()
        assert self.seq1_twobit.split("C") == self.seq1_fasta.split("C")
        assert self.seq1_twobit.split("C", 1) == self.seq1_fasta.split("C", 1)

    def test_rsplit(self):
        assert self.seq1_twobit.rsplit() == self.seq1_fasta.rsplit()
        assert self.seq1_twobit.rsplit("C") == self.seq1_fasta.rsplit("C")
        assert self.seq1_twobit.rsplit("C", 1) == self.seq1_fasta.rsplit("C", 1)

    def test_strip(self):
        assert self.seq1_twobit.strip("G") == self.seq1_fasta.strip("G")

    def test_lstrip(self, chars=None):
        assert self.seq1_twobit.lstrip("G") == self.seq1_fasta.lstrip("G")

    def test_rstrip(self, chars=None):
        assert self.seq1_twobit.rstrip("G") == self.seq1_fasta.rstrip("G")

    def test_upper(self):
        seq1_twobit_upper = self.seq1_twobit.upper()
        seq1_fasta_upper = self.seq1_fasta.upper()
        assert seq1_twobit_upper == seq1_fasta_upper
        assert seq1_twobit_upper[140:210] == seq1_fasta_upper[140:210]
        seq2_twobit_upper = self.seq2_twobit.upper()
        seq2_fasta_upper = self.seq2_fasta.upper()
        assert seq2_twobit_upper == seq2_fasta_upper
        assert seq2_twobit_upper[140:210] == seq2_fasta_upper[140:210]

    def test_lower(self):
        seq1_twobit_lower = self.seq1_twobit.lower()
        seq1_fasta_lower = self.seq1_fasta.lower()
        assert seq1_twobit_lower == seq1_fasta_lower
        assert seq1_twobit_lower[140:210] == seq1_fasta_lower[140:210]
        seq2_twobit_lower = self.seq2_twobit.lower()
        seq2_fasta_lower = self.seq2_fasta.lower()
        assert seq2_twobit_lower == seq2_fasta_lower
        assert seq2_twobit_lower[140:210] == seq2_fasta_lower[140:210]

    def test_isupper(self):
        assert self.seq1_twobit.isupper() == self.seq1_fasta.isupper()
        assert self.seq2_twobit.isupper() == self.seq2_fasta.isupper()
        assert self.record1_twobit.isupper() == self.record1_fasta.isupper()
        assert self.record2_twobit.isupper() == self.record2_fasta.isupper()

    def test_islower(self):
        assert self.seq1_twobit.islower() == self.seq1_fasta.islower()
        assert self.seq2_twobit.islower() == self.seq2_fasta.islower()
        assert self.record1_twobit.islower() == self.record1_fasta.islower()
        assert self.record2_twobit.islower() == self.record2_fasta.islower()

    def test_replace(self):
        # seq.transcribe uses seq._data.replace
        assert self.seq1_twobit.transcribe() == self.seq1_fasta.transcribe()

    def test_translate(self):
        # seq.reverse_complement uses seq._data.translate
        assert self.seq1_twobit.reverse_complement() == self.seq1_fasta.reverse_complement()
        record1_twobit = self.record1_twobit.reverse_complement()
        record1_fasta = self.record1_fasta.reverse_complement()
        assert record1_twobit.seq == record1_fasta.seq

    def test_defined(self):
        assert self.seq1_twobit.defined
        assert self.seq2_twobit.defined
        assert self.seq1_twobit.defined_ranges == ((0, len(self.seq1_twobit)),)
        assert self.seq2_twobit.defined_ranges == ((0, len(self.seq2_twobit)),)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
