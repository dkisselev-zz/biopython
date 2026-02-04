# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for seq module."""

import copy
import unittest
import pytest
import warnings

try:
    import numpy as np
except ImportError:
    np = None

from Bio import BiopythonWarning
from Bio import Seq
from Bio.Data.CodonTable import standard_dna_table
from Bio.Data.CodonTable import TranslationError
from Bio.Data.IUPACData import ambiguous_dna_complement
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Data.IUPACData import ambiguous_rna_complement
from Bio.Data.IUPACData import ambiguous_rna_values

test_seqs = [
    Seq.Seq("TCAAAAGGATGCATCATG"),
    Seq.Seq("T"),
    Seq.Seq("ATGAAACTG"),
    Seq.Seq("ATGAARCTG"),
    Seq.Seq("AWGAARCKG"),  # Note no U or T
    Seq.Seq("".join(ambiguous_rna_values)),
    Seq.Seq("".join(ambiguous_dna_values)),
    Seq.Seq("AWGAARCKG"),
    Seq.Seq("AUGAAACUG"),
    Seq.Seq("ATGAAA-CTG"),
    Seq.Seq("ATGAAACTGWN"),
    Seq.Seq("AUGAAA==CUG"),
    Seq.Seq("AUGAAACUGWN"),
    Seq.Seq("AUGAAACTG"),  # U and T
    Seq.MutableSeq("ATGAAACTG"),
    Seq.MutableSeq("AUGaaaCUG"),
    Seq.Seq("ACTGTCGTCT"),
]
protein_seqs = [
    Seq.Seq("ATCGPK"),
    Seq.Seq("T.CGPK"),
    Seq.Seq("T-CGPK"),
    Seq.Seq("MEDG-KRXR*"),
    Seq.MutableSeq("ME-K-DRXR*XU"),
    Seq.Seq("MEDG-KRXR@"),
    Seq.Seq("ME-KR@"),
    Seq.Seq("MEDG.KRXR@"),
]


class TestSeq(unittest.TestCase):
    def setUp(self):
        self.s = Seq.Seq("TCAAAAGGATGCATCATG")

    def test_as_string(self):
        """Test converting Seq to string."""
        assert "TCAAAAGGATGCATCATG" == self.s

    def test_seq_construction(self):
        """Test Seq object initialization."""
        sequence = bytes(self.s)
        s = Seq.Seq(sequence)
        assert isinstance(s, Seq.Seq)
        assert s == self.s
        s = Seq.Seq(bytearray(sequence))
        assert isinstance(s, Seq.Seq)
        assert s == self.s
        s = Seq.Seq(sequence.decode("ASCII"))
        assert isinstance(s, Seq.Seq)
        assert s == self.s
        s = Seq.Seq(self.s)
        assert isinstance(s, Seq.Seq)
        assert s == self.s
        s = Seq.Seq(Seq.MutableSeq(sequence))
        assert isinstance(s, Seq.Seq)
        assert s == self.s
        with pytest.raises(UnicodeEncodeError):
            Seq.Seq("ÄþÇÐ")  # All are Latin-1 characters
        with pytest.raises(UnicodeEncodeError):
            Seq.Seq("あいうえお")  # These are not

    def test_repr(self):
        """Test representation of Seq object."""
        assert "Seq('TCAAAAGGATGCATCATG')" == repr(self.s)

    def test_truncated_repr(self):
        seq = "TCAAAAGGATGCATCATGTCAAAAGGATGCATCATGTCAAAAGGATGCATCATGTCAAAAGGA"
        expected = "Seq('TCAAAAGGATGCATCATGTCAAAAGGATGCATCATGTCAAAAGGATGCATCATG...GGA')"
        assert expected == repr(Seq.Seq(seq))

    def test_length(self):
        """Test len method on Seq object."""
        assert 18 == len(self.s)

    def test_first_nucleotide(self):
        """Test getting first nucleotide of Seq."""
        assert "T" == self.s[0]

    def test_last_nucleotide(self):
        """Test getting last nucleotide of Seq."""
        assert "G" == self.s[-1]

    def test_slicing(self):
        """Test slicing of Seq."""
        assert "AA" == self.s[3:5]

    def test_reverse(self):
        """Test reverse using -1 stride."""
        assert "GTACTACGTAGGAAAACT" == self.s[::-1]

    def test_extract_third_nucleotide(self):
        """Test extracting every third nucleotide (slicing with stride 3)."""
        assert "TAGTAA" == self.s[0::3]
        assert "CAGGTT" == self.s[1::3]
        assert "AAACCG" == self.s[2::3]

    def test_concatenation_of_seq(self):
        t = Seq.Seq("T")
        u = self.s + t
        assert str(self.s) + "T" == u
        assert self.s + Seq.Seq("T") == "TCAAAAGGATGCATCATGT"

    def test_replace(self):
        assert "ATCCCA" == Seq.Seq("ATC-CCA").replace("-", "")

    def test_cast_to_list(self):
        assert list("ATC") == list(Seq.Seq("ATC"))
        assert list("ATC") == list(Seq.MutableSeq("ATC"))
        assert list("") == list(Seq.MutableSeq(""))
        assert list("") == list(Seq.Seq(""))
        with pytest.raises(Seq.UndefinedSequenceError):
            list(Seq.Seq(None, length=3))
        with pytest.raises(Seq.UndefinedSequenceError):
            list(Seq.Seq({3: "ACGT"}, length=10))


class TestSeqStringMethods(unittest.TestCase):
    def setUp(self):
        self.s = Seq.Seq("TCAAAAGGATGCATCATG")
        self.dna = [
            Seq.Seq("ATCG"),
            Seq.Seq("gtca"),
            Seq.MutableSeq("GGTCA"),
            Seq.Seq("CTG-CA"),
        ]
        self.rna = [
            Seq.Seq("AUUUCG"),
            Seq.MutableSeq("AUUCG"),
            Seq.Seq("uCAg"),
            Seq.MutableSeq("UC-AG"),
            Seq.Seq("U.CAG"),
        ]
        self.nuc = [Seq.Seq("ATCG")]
        self.protein = [
            Seq.Seq("ATCGPK"),
            Seq.Seq("atcGPK"),
            Seq.Seq("T.CGPK"),
            Seq.Seq("T-CGPK"),
            Seq.Seq("MEDG-KRXR*"),
            Seq.MutableSeq("ME-K-DRXR*XU"),
            Seq.Seq("MEDG-KRXR@"),
            Seq.Seq("ME-KR@"),
            Seq.Seq("MEDG.KRXR@"),
        ]
        self.test_chars = ["-", Seq.Seq("-"), Seq.Seq("*"), "-X@"]

    def test_string_methods(self):
        for a in self.dna + self.rna + self.nuc + self.protein:
            assert a.lower() == str(a).lower()
            assert a.upper() == str(a).upper()
            assert a.islower() == str(a).islower()
            assert a.isupper() == str(a).isupper()
            assert a.strip() == str(a).strip()
            assert a.lstrip() == str(a).lstrip()
            assert a.rstrip() == str(a).rstrip()

    def test_mutableseq_upper_lower(self):
        seq = Seq.MutableSeq("ACgt")
        lseq = seq.lower()
        assert lseq == "acgt"
        assert seq == "ACgt"
        assert lseq.islower()
        assert not seq.islower()
        lseq = seq.lower(inplace=False)
        assert lseq == "acgt"
        assert seq == "ACgt"
        assert lseq.islower()
        assert not seq.islower()
        lseq = seq.lower(inplace=True)
        assert lseq == "acgt"
        assert lseq is seq
        assert lseq.islower()
        assert lseq.islower()
        seq = Seq.MutableSeq("ACgt")
        useq = seq.upper()
        assert useq == "ACGT"
        assert seq == "ACgt"
        assert useq.isupper()
        assert not seq.isupper()
        useq = seq.upper(inplace=False)
        assert useq == "ACGT"
        assert seq == "ACgt"
        assert useq.isupper()
        assert not seq.isupper()
        useq = seq.upper(inplace=True)
        assert useq == "ACGT"
        assert useq is seq
        assert useq.isupper()
        assert seq.isupper()

    def test_hash(self):
        with warnings.catch_warnings(record=True):
            hash(self.s)

    def test_not_equal_comparsion(self):
        """Test __ne__ comparison method."""
        assert Seq.Seq("TCAAA") != Seq.Seq("TCAAAA")

    def test_less_than_comparison(self):
        """Test __lt__ comparison method."""
        assert self.s[:-1] < self.s

    def test_less_than_comparison_of_incompatible_types(self):
        """Test incompatible types __lt__ comparison method."""
        with pytest.raises(TypeError):
            self.s < 1

    def test_less_than_or_equal_comparison(self):
        """Test __le__ comparison method."""
        assert self.s <= self.s

    def test_less_than_or_equal_comparison_of_incompatible_types(self):
        """Test incompatible types __le__ comparison method."""
        with pytest.raises(TypeError):
            self.s <= 1

    def test_greater_than_comparison(self):
        """Test __gt__ comparison method."""
        assert self.s > self.s[:-1]

    def test_greater_than_comparison_of_incompatible_types(self):
        """Test incompatible types __gt__ comparison method."""
        with pytest.raises(TypeError):
            self.s > 1

    def test_greater_than_or_equal_comparison(self):
        """Test __ge__ comparison method."""
        assert self.s >= self.s

    def test_greater_than_or_equal_comparison_of_incompatible_types(self):
        """Test incompatible types __ge__ comparison method."""
        with pytest.raises(TypeError):
            self.s >= 1

    def test_add_method_using_wrong_object(self):
        with pytest.raises(TypeError):
            self.s + {}

    def test_radd_method_using_wrong_object(self):
        assert self.s.__radd__({}) == NotImplemented

    def test_contains_method(self):
        assert "AAAA" in self.s

    def test_startswith(self):
        assert self.s.startswith("TCA")
        assert self.s.startswith(("CAA", "CTA"), 1)

    def test_endswith(self):
        assert self.s.endswith("ATG")
        assert self.s.endswith(("ATG", "CTA"))

    def test_append_nucleotides(self):
        self.test_chars.append(Seq.Seq("A"))
        assert 5 == len(self.test_chars)

    def test_append_proteins(self):
        self.test_chars.append(Seq.Seq("K"))
        self.test_chars.append(Seq.Seq("K-"))
        self.test_chars.append(Seq.Seq("K@"))

        assert 7 == len(self.test_chars)

    def test_stripping_characters(self):
        for a in self.dna + self.rna + self.nuc + self.protein:
            for char in self.test_chars:
                str_char = str(char)
                assert a.strip(char) == str(a).strip(str_char)
                assert a.lstrip(char) == str(a).lstrip(str_char)
                assert a.rstrip(char) == str(a).rstrip(str_char)
                try:
                    removeprefix = str(a).removeprefix(str_char)
                    removesuffix = str(a).removesuffix(str_char)
                except AttributeError:
                    if str(a).startswith(str_char):
                        removeprefix = str(a)[len(str_char) :]
                    else:
                        removeprefix = str(a)
                    if str_char and str(a).endswith(str_char):
                        removesuffix = str(a)[: -len(str_char)]
                    else:
                        removesuffix = str(a)

                assert a.removeprefix(char) == removeprefix
                assert a.removesuffix(char) == removesuffix

    def test_finding_characters(self):
        for a in self.dna + self.rna + self.nuc + self.protein:
            for char in self.test_chars:
                str_char = str(char)
                assert a.find(char) == str(a).find(str_char)
                assert a.find(char, 2, -2) == str(a).find(str_char, 2, -2)
                assert a.rfind(char) == str(a).rfind(str_char)
                assert a.rfind(char, 2, -2) == str(a).rfind(str_char, 2, -2)

    def test_counting_characters(self):
        from Bio.SeqRecord import SeqRecord

        for a in self.dna + self.rna + self.nuc + self.protein:
            r = SeqRecord(a)
            for char in self.test_chars:
                str_char = str(char)
                n = str(a).count(str_char)
                assert a.count(char) == n
                assert r.count(char) == n
                n = str(a).count(str_char, 2, -2)
                assert a.count(char, 2, -2) == n
                assert r.count(char, 2, -2) == n

    def test_splits(self):
        for a in self.dna + self.rna + self.nuc + self.protein:
            for char in self.test_chars:
                str_char = str(char)
                assert a.split(char) == str(a).split(str_char)
                assert a.rsplit(char) == str(a).rsplit(str_char)

                for max_sep in [0, 1, 2, 999]:
                    assert a.split(char, max_sep) == str(a).split(str_char, max_sep)


class TestSeqAddition(unittest.TestCase):
    def setUp(self):
        self.dna = [
            Seq.Seq("ATCG"),
            Seq.Seq("gtca"),
            Seq.MutableSeq("GGTCA"),
            Seq.Seq("CTG-CA"),
            "TGGTCA",
        ]
        self.rna = [
            Seq.Seq("AUUUCG"),
            Seq.MutableSeq("AUUCG"),
            Seq.Seq("uCAg"),
            Seq.MutableSeq("UC-AG"),
            Seq.Seq("U.CAG"),
            "UGCAU",
        ]
        self.nuc = [Seq.Seq("ATCG"), "UUUTTTACG"]
        self.protein = [
            Seq.Seq("ATCGPK"),
            Seq.Seq("atcGPK"),
            Seq.Seq("T.CGPK"),
            Seq.Seq("T-CGPK"),
            Seq.Seq("MEDG-KRXR*"),
            Seq.MutableSeq("ME-K-DRXR*XU"),
            "TEDDF",
        ]

    def test_addition_dna_rna_with_generic_nucleotides(self):
        for a in self.dna + self.rna:
            for b in self.nuc:
                c = a + b
                assert c == str(a) + str(b)

    def test_addition_dna_rna_with_generic_nucleotides_inplace(self):
        for a in self.dna + self.rna:
            for b in self.nuc:
                c = b + a
                b += a  # can't change 'a' as need value next iteration
                assert c == b

    def test_addition_rna_with_rna(self):
        self.rna.pop(3)
        for a in self.rna:
            for b in self.rna:
                c = a + b
                assert c == str(a) + str(b)

    def test_addition_rna_with_rna_inplace(self):
        self.rna.pop(3)
        for a in self.rna:
            for b in self.rna:
                c = b + a
                b += a
                assert c == b

    def test_addition_dna_with_dna(self):
        for a in self.dna:
            for b in self.dna:
                c = a + b
                assert c == str(a) + str(b)

    def test_addition_dna_with_dna_inplace(self):
        for a in self.dna:
            for b in self.dna:
                c = b + a
                b += a
                assert c == b

    def test_addition_dna_with_rna(self):
        self.dna.pop(4)
        self.rna.pop(5)
        for a in self.dna:
            for b in self.rna:
                assert str(a) + str(b) == a + b
                assert str(b) + str(a) == b + a
                # Check in place works
                c = a
                c += b
                assert c == str(a) + str(b)
                c = b
                c += a
                assert c == str(b) + str(a)

    def test_addition_proteins(self):
        self.protein.pop(2)
        for a in self.protein:
            for b in self.protein:
                c = a + b
                assert c == str(a) + str(b)

    def test_addition_proteins_inplace(self):
        self.protein.pop(2)
        for a in self.protein:
            for b in self.protein:
                c = b + a
                b += a
                assert c == b

    def test_adding_protein_with_nucleotides(self):
        for a in self.protein[0:5]:
            for b in self.dna[0:3] + self.rna[0:4]:
                assert str(a) + str(b) == a + b
                a += b

    def test_adding_generic_nucleotide_with_other_nucleotides(self):
        for a in self.nuc:
            for b in self.dna + self.rna + self.nuc:
                c = a + b
                assert c == str(a) + str(b)

    def test_adding_generic_nucleotide_with_other_nucleotides_inplace(self):
        for a in self.nuc:
            for b in self.dna + self.rna + self.nuc:
                c = b + a
                b += a
                assert c == b


class TestSeqMultiplication(unittest.TestCase):
    def test_mul_method(self):
        """Test mul method; relies on addition method."""
        for seq in test_seqs + protein_seqs:
            assert seq * 3 == seq + seq + seq
        if np is not None:
            factor = np.intc(3)  # numpy integer
            for seq in test_seqs + protein_seqs:
                assert seq * factor == seq + seq + seq

    def test_mul_method_exceptions(self):
        """Test mul method exceptions."""
        for seq in test_seqs + protein_seqs:
            with pytest.raises(TypeError):
                seq * 3.0
            with pytest.raises(TypeError):
                seq * ""

    def test_rmul_method(self):
        """Test rmul method; relies on addition method."""
        for seq in test_seqs + protein_seqs:
            assert 3 * seq == seq + seq + seq
        if np is not None:
            factor = np.intc(3)  # numpy integer
            for seq in test_seqs + protein_seqs:
                assert factor * seq == seq + seq + seq

    def test_rmul_method_exceptions(self):
        """Test rmul method exceptions."""
        for seq in test_seqs + protein_seqs:
            with pytest.raises(TypeError):
                3.0 * seq
            with pytest.raises(TypeError):
                "" * seq

    def test_imul_method(self):
        """Test imul method; relies on addition and mull methods."""
        for seq in test_seqs + protein_seqs:
            original_seq = seq * 1  # make a copy
            seq *= 3
            assert seq == original_seq + original_seq + original_seq
        if np is not None:
            factor = np.intc(3)  # numpy integer
            for seq in test_seqs + protein_seqs:
                original_seq = seq * 1  # make a copy
                seq *= factor
                assert seq == original_seq + original_seq + original_seq

    def test_imul_method_exceptions(self):
        """Test imul method exceptions."""
        for seq in test_seqs + protein_seqs:
            with pytest.raises(TypeError):
                seq *= 3.0
            with pytest.raises(TypeError):
                seq *= ""


class TestMutableSeq(unittest.TestCase):
    def setUp(self):
        sequence = b"TCAAAAGGATGCATCATG"
        self.s = Seq.Seq(sequence)
        self.mutable_s = Seq.MutableSeq(sequence)

    def test_mutableseq_construction(self):
        """Test MutableSeq object initialization."""
        sequence = bytes(self.s)
        mutable_s = Seq.MutableSeq(sequence)
        assert isinstance(mutable_s, Seq.MutableSeq)
        assert mutable_s == self.s
        mutable_s = Seq.MutableSeq(bytearray(sequence))
        assert isinstance(mutable_s, Seq.MutableSeq)
        assert mutable_s == self.s
        mutable_s = Seq.MutableSeq(sequence.decode("ASCII"))
        assert isinstance(mutable_s, Seq.MutableSeq)
        assert mutable_s == self.s
        mutable_s = Seq.MutableSeq(self.s)
        assert isinstance(mutable_s, Seq.MutableSeq)
        assert mutable_s == self.s
        mutable_s = Seq.MutableSeq(Seq.MutableSeq(sequence))
        assert mutable_s == self.s
        assert isinstance(mutable_s, Seq.MutableSeq)
        with pytest.raises(UnicodeEncodeError):
            Seq.MutableSeq("ÄþÇÐ")  # All are Latin-1 characters
        with pytest.raises(UnicodeEncodeError):
            Seq.MutableSeq("あいうえお")  # These are not

    def test_repr(self):
        assert "MutableSeq('TCAAAAGGATGCATCATG')" == repr(self.mutable_s)

    def test_truncated_repr(self):
        seq = "TCAAAAGGATGCATCATGTCAAAAGGATGCATCATGTCAAAAGGATGCATCATGTCAAAAGGA"
        expected = (
            "MutableSeq('TCAAAAGGATGCATCATGTCAAAAGGATGCATCATGTCAAAAGGATGCATCATG...GGA')"
        )
        assert expected == repr(Seq.MutableSeq(seq))

    def test_equal_comparison(self):
        """Test __eq__ comparison method."""
        assert self.mutable_s == "TCAAAAGGATGCATCATG"

    def test_not_equal_comparison(self):
        """Test __ne__ comparison method."""
        assert self.mutable_s != "other thing"

    def test_less_than_comparison(self):
        """Test __lt__ comparison method."""
        assert self.mutable_s[:-1] < self.mutable_s

    def test_less_than_comparison_of_incompatible_types(self):
        with pytest.raises(TypeError):
            self.mutable_s < 1

    def test_less_than_comparison_with_str(self):
        assert self.mutable_s[:-1] <= "TCAAAAGGATGCATCATG"

    def test_less_than_or_equal_comparison(self):
        """Test __le__ comparison method."""
        assert self.mutable_s[:-1] <= self.mutable_s

    def test_less_than_or_equal_comparison_of_incompatible_types(self):
        with pytest.raises(TypeError):
            self.mutable_s <= 1

    def test_less_than_or_equal_comparison_with_str(self):
        assert self.mutable_s[:-1] <= "TCAAAAGGATGCATCATG"

    def test_greater_than_comparison(self):
        """Test __gt__ comparison method."""
        assert self.mutable_s > self.mutable_s[:-1]

    def test_greater_than_comparison_of_incompatible_types(self):
        with pytest.raises(TypeError):
            self.mutable_s > 1

    def test_greater_than_comparison_with_str(self):
        assert self.mutable_s > "TCAAAAGGATGCATCAT"

    def test_greater_than_or_equal_comparison(self):
        """Test __ge__ comparison method."""
        assert self.mutable_s >= self.mutable_s

    def test_greater_than_or_equal_comparison_of_incompatible_types(self):
        with pytest.raises(TypeError):
            self.mutable_s >= 1

    def test_greater_than_or_equal_comparison_with_str(self):
        assert self.mutable_s >= "TCAAAAGGATGCATCATG"

    def test_add_method(self):
        """Test adding wrong type to MutableSeq."""
        with pytest.raises(TypeError):
            self.mutable_s + 1234

    def test_radd_method_wrong_type(self):
        assert self.mutable_s.__radd__(1234) == NotImplemented

    def test_contains_method(self):
        assert "AAAA" in self.mutable_s

    def test_startswith(self):
        assert self.mutable_s.startswith("TCA")
        assert self.mutable_s.startswith(("CAA", "CTA"), 1)

    def test_endswith(self):
        assert self.mutable_s.endswith("ATG")
        assert self.mutable_s.endswith(("ATG", "CTA"))

    def test_as_string(self):
        assert "TCAAAAGGATGCATCATG" == self.mutable_s

    def test_length(self):
        assert 18 == len(self.mutable_s)

    def test_converting_to_immutable(self):
        assert isinstance(Seq.Seq(self.mutable_s), Seq.Seq)

    def test_first_nucleotide(self):
        assert "T" == self.mutable_s[0]

    def test_setting_slices(self):
        assert Seq.MutableSeq("CAAA") == self.mutable_s[1:5], "Slice mutable seq"

        self.mutable_s[1:3] = "GAT"
        assert Seq.MutableSeq("TGATAAAGGATGCATCATG") == self.mutable_s, "Set slice with string and adding extra nucleotide"

        self.mutable_s[1:3] = self.mutable_s[5:7]
        assert Seq.MutableSeq("TAATAAAGGATGCATCATG") == self.mutable_s, "Set slice with MutableSeq"
        if np is not None:
            one, three, five, seven = np.array([1, 3, 5, 7])  # numpy integers
            assert Seq.MutableSeq("AATA") == self.mutable_s[one:five], "Slice mutable seq"

            self.mutable_s[one:three] = "GAT"
            assert Seq.MutableSeq("TGATTAAAGGATGCATCATG") == self.mutable_s, "Set slice with string and adding extra nucleotide"

            self.mutable_s[one:three] = self.mutable_s[five:seven]
            assert Seq.MutableSeq("TAATTAAAGGATGCATCATG") == self.mutable_s, "Set slice with MutableSeq"

    def test_setting_item(self):
        self.mutable_s[3] = "G"
        assert Seq.MutableSeq("TCAGAAGGATGCATCATG") == self.mutable_s
        if np is not None:
            i = np.intc(3)
            self.mutable_s[i] = "X"
            assert Seq.MutableSeq("TCAXAAGGATGCATCATG") == self.mutable_s

    def test_deleting_slice(self):
        del self.mutable_s[4:5]
        assert Seq.MutableSeq("TCAAAGGATGCATCATG") == self.mutable_s

    def test_deleting_item(self):
        del self.mutable_s[3]
        assert Seq.MutableSeq("TCAAAGGATGCATCATG") == self.mutable_s

    def test_appending(self):
        self.mutable_s.append("C")
        assert Seq.MutableSeq("TCAAAAGGATGCATCATGC") == self.mutable_s

    def test_inserting(self):
        self.mutable_s.insert(4, "G")
        assert Seq.MutableSeq("TCAAGAAGGATGCATCATG") == self.mutable_s

    def test_popping_last_item(self):
        assert "G" == self.mutable_s.pop()

    def test_remove_items(self):
        self.mutable_s.remove("G")
        assert Seq.MutableSeq("TCAAAAGATGCATCATG") == self.mutable_s, "Remove first G"

        with pytest.raises(ValueError):
            self.mutable_s.remove("Z")

    def test_count(self):
        assert 7 == self.mutable_s.count("A")
        assert 2 == self.mutable_s.count("AA")

    def test_index(self):
        assert 2 == self.mutable_s.index("A")
        with pytest.raises(ValueError):
            self.mutable_s.index("8888")

    def test_reverse(self):
        """Test using reverse method."""
        self.mutable_s.reverse()
        assert Seq.MutableSeq("GTACTACGTAGGAAAACT") == self.mutable_s

    def test_reverse_with_stride(self):
        """Test reverse using -1 stride."""
        assert Seq.MutableSeq("GTACTACGTAGGAAAACT") == self.mutable_s[::-1]

    def test_complement(self):
        self.mutable_s.complement(inplace=True)
        assert "AGTTTTCCTACGTAGTAC" == self.mutable_s

    def test_complement_rna(self):
        m = self.mutable_s.complement_rna()
        assert self.mutable_s == "TCAAAAGGATGCATCATG"
        assert isinstance(m, Seq.MutableSeq)
        assert m == "AGUUUUCCUACGUAGUAC"
        m = self.mutable_s.complement_rna(inplace=True)
        assert self.mutable_s == "AGUUUUCCUACGUAGUAC"
        assert isinstance(m, Seq.MutableSeq)
        assert m == "AGUUUUCCUACGUAGUAC"

    def test_reverse_complement_rna(self):
        m = self.mutable_s.reverse_complement_rna()
        assert self.mutable_s == "TCAAAAGGATGCATCATG"
        assert isinstance(m, Seq.MutableSeq)
        assert m == "CAUGAUGCAUCCUUUUGA"
        m = self.mutable_s.reverse_complement_rna(inplace=True)
        assert self.mutable_s == "CAUGAUGCAUCCUUUUGA"
        assert isinstance(m, Seq.MutableSeq)
        assert m == "CAUGAUGCAUCCUUUUGA"

    def test_transcribe(self):
        r = self.mutable_s.transcribe()
        assert self.mutable_s == "TCAAAAGGATGCATCATG"
        assert isinstance(r, Seq.MutableSeq)
        assert r == "UCAAAAGGAUGCAUCAUG"
        r = self.mutable_s.transcribe(inplace=True)
        assert self.mutable_s == "UCAAAAGGAUGCAUCAUG"
        assert isinstance(r, Seq.MutableSeq)
        assert r == "UCAAAAGGAUGCAUCAUG"
        d = self.mutable_s.back_transcribe()
        assert self.mutable_s == "UCAAAAGGAUGCAUCAUG"
        assert isinstance(d, Seq.MutableSeq)
        assert d == "TCAAAAGGATGCATCATG"
        d = self.mutable_s.back_transcribe(inplace=True)
        assert self.mutable_s == "TCAAAAGGATGCATCATG"
        assert isinstance(d, Seq.MutableSeq)
        assert d == "TCAAAAGGATGCATCATG"

    def test_complement_mixed_aphabets(self):
        seq = Seq.MutableSeq("AUGaaaCTG")
        seq.complement_rna(inplace=True)
        assert "UACuuuGAC" == seq

    def test_complement_rna_string(self):
        seq = Seq.MutableSeq("AUGaaaCUG")
        seq.complement_rna(inplace=True)
        assert "UACuuuGAC" == seq

    def test_complement_dna_string(self):
        seq = Seq.MutableSeq("ATGaaaCTG")
        seq.complement(inplace=True)
        assert "TACtttGAC" == seq

    def test_reverse_complement(self):
        self.mutable_s.reverse_complement(inplace=True)
        assert "CATGATGCATCCTTTTGA" == self.mutable_s

    def test_extend_method(self):
        self.mutable_s.extend("GAT")
        assert Seq.MutableSeq("TCAAAAGGATGCATCATGGAT") == self.mutable_s

    def test_extend_with_mutable_seq(self):
        self.mutable_s.extend(Seq.MutableSeq("TTT"))
        assert Seq.MutableSeq("TCAAAAGGATGCATCATGTTT") == self.mutable_s

    def test_delete_stride_slice(self):
        del self.mutable_s[4 : 6 - 1]
        assert Seq.MutableSeq("TCAAAGGATGCATCATG") == self.mutable_s

    def test_extract_third_nucleotide(self):
        """Test extracting every third nucleotide (slicing with stride 3)."""
        assert Seq.MutableSeq("TAGTAA") == self.mutable_s[0::3]
        assert Seq.MutableSeq("CAGGTT") == self.mutable_s[1::3]
        assert Seq.MutableSeq("AAACCG") == self.mutable_s[2::3]

    def test_set_wobble_codon_to_n(self):
        """Test setting wobble codon to N (set slice with stride 3)."""
        self.mutable_s[2::3] = "N" * len(self.mutable_s[2::3])
        assert Seq.MutableSeq("TCNAANGGNTGNATNATN") == self.mutable_s
        if np is not None:
            start, step = np.array([2, 3])  # numpy integers
            self.mutable_s[start::step] = "X" * len(self.mutable_s[2::3])
            assert Seq.MutableSeq("TCXAAXGGXTGXATXATX") == self.mutable_s


class TestAmbiguousComplements(unittest.TestCase):
    def test_ambiguous_values(self):
        """Test that other tests do not introduce characters to our values."""
        assert "-" not in ambiguous_dna_values
        assert "?" not in ambiguous_dna_values


class TestComplement(unittest.TestCase):
    def test_complement_ambiguous_dna_values(self):
        for ambig_char, values in sorted(ambiguous_dna_values.items()):
            compl_values = Seq.Seq(values).complement()
            ambig_values = ambiguous_dna_values[ambiguous_dna_complement[ambig_char]]
            assert sorted(compl_values) == sorted(ambig_values)

    def test_complement_ambiguous_rna_values(self):
        for ambig_char, values in sorted(ambiguous_rna_values.items()):
            # Will default to DNA if neither T nor U found...
            if "u" in values or "U" in values:
                compl_values = Seq.Seq(values).complement_rna().transcribe()
            else:
                compl_values = Seq.Seq(values).complement().transcribe()
            ambig_values = ambiguous_rna_values[ambiguous_rna_complement[ambig_char]]
            assert sorted(compl_values) == sorted(ambig_values)

    def test_complement_incompatible_letters(self):
        seq = Seq.Seq("CAGGTU")
        dna = seq.complement()
        assert "GTCCAA" == dna
        rna = seq.complement_rna()
        assert "GUCCAA" == rna

    def test_complement_of_mixed_dna_rna(self):
        seq = "AUGAAACTG"  # U and T
        dna = Seq.complement(seq)
        assert "TACTTTGAC" == dna
        rna = Seq.complement_rna(seq)
        assert "UACUUUGAC" == rna

    def test_complement_of_rna(self):
        seq = "AUGAAACUG"
        rna = Seq.complement_rna(seq)
        assert "UACUUUGAC" == rna

    def test_complement_of_dna(self):
        seq = "ATGAAACTG"
        assert "TACTTTGAC" == Seq.complement(seq)

    def test_immutable(self):
        from Bio.SeqRecord import SeqRecord

        r = SeqRecord(Seq.Seq("ACGT"))
        with pytest.raises(TypeError) as cm:
            Seq.complement(r, inplace=True)
        assert str(cm.value) == "SeqRecords are immutable"
        with pytest.raises(TypeError) as cm:
            Seq.complement("ACGT", inplace=True)
        assert str(cm.value) == "strings are immutable"
        with pytest.raises(TypeError) as cm:
            Seq.complement_rna(r, inplace=True)
        assert str(cm.value) == "SeqRecords are immutable"
        with pytest.raises(TypeError) as cm:
            Seq.complement_rna("ACGT", inplace=True)
        assert str(cm.value) == "strings are immutable"


class TestReverseComplement(unittest.TestCase):
    def test_reverse_complement(self):
        test_seqs_copy = copy.copy(test_seqs)
        test_seqs_copy.pop(13)

        for nucleotide_seq in test_seqs_copy:
            if not isinstance(nucleotide_seq, Seq.Seq):
                continue
            if "u" in nucleotide_seq or "U" in nucleotide_seq:
                expected = Seq.reverse_complement_rna(nucleotide_seq)
                assert repr(expected) == repr(nucleotide_seq.reverse_complement_rna())
                assert repr(expected[::-1]) == repr(nucleotide_seq.complement_rna())
                assert nucleotide_seq.complement_rna() == Seq.reverse_complement_rna(nucleotide_seq)[::-1]
                assert nucleotide_seq.reverse_complement_rna() == Seq.reverse_complement_rna(nucleotide_seq)
            else:
                expected = Seq.reverse_complement(nucleotide_seq)
                assert repr(expected) == repr(nucleotide_seq.reverse_complement())
                assert repr(expected[::-1]) == repr(nucleotide_seq.complement())
                assert nucleotide_seq.complement() == Seq.reverse_complement(nucleotide_seq)[::-1]
                assert nucleotide_seq.reverse_complement() == Seq.reverse_complement(nucleotide_seq)

    def test_reverse_complement_of_mixed_dna_rna(self):
        seq = "AUGAAACTG"  # U and T
        dna = Seq.reverse_complement(seq)
        assert "CAGTTTCAT" == dna
        rna = Seq.reverse_complement_rna(seq)
        assert "CAGUUUCAU" == rna

    def test_reverse_complement_of_rna(self):
        seq = "AUGAAACUG"
        dna = Seq.reverse_complement(seq)
        assert "CAGTTTCAT" == dna

    def test_reverse_complement_of_dna(self):
        seq = "ATGAAACTG"
        assert "CAGTTTCAT" == Seq.reverse_complement(seq)

    def test_immutable(self):
        from Bio.SeqRecord import SeqRecord

        r = SeqRecord(Seq.Seq("ACGT"))
        with pytest.raises(TypeError) as cm:
            Seq.reverse_complement(r, inplace=True)
        assert str(cm.value) == "SeqRecords are immutable"
        with pytest.raises(TypeError) as cm:
            Seq.reverse_complement("ACGT", inplace=True)
        assert str(cm.value) == "strings are immutable"
        with pytest.raises(TypeError) as cm:
            Seq.reverse_complement_rna(r, inplace=True)
        assert str(cm.value) == "SeqRecords are immutable"
        with pytest.raises(TypeError) as cm:
            Seq.reverse_complement_rna("ACGT", inplace=True)
        assert str(cm.value) == "strings are immutable"


class TestDoubleReverseComplement(unittest.TestCase):
    def test_reverse_complements(self):
        """Test double reverse complement preserves the sequence."""
        sorted_amb_rna = sorted(ambiguous_rna_values)
        sorted_amb_dna = sorted(ambiguous_dna_values)
        for sequence in [
            Seq.Seq("".join(sorted_amb_dna)),
            Seq.Seq("".join(sorted_amb_dna).replace("X", "")),
            Seq.Seq("AWGAARCKG"),  # Note no U or T
        ]:
            reversed_sequence = sequence.reverse_complement()
            assert sequence == reversed_sequence.reverse_complement()
        for sequence in [
            Seq.Seq("".join(sorted_amb_rna)),
            Seq.Seq("".join(sorted_amb_rna).replace("X", "")),
            Seq.Seq("AWGAARCKG"),  # Note no U or T
        ]:
            reversed_sequence = sequence.reverse_complement_rna()
            assert sequence == reversed_sequence.reverse_complement_rna()


class TestTranscription(unittest.TestCase):
    def test_transcription_dna_into_rna(self):
        for nucleotide_seq in test_seqs:
            expected = Seq.transcribe(nucleotide_seq)
            assert str(nucleotide_seq).replace("t", "u").replace("T", "U") == expected

    def test_transcription_dna_string_into_rna(self):
        seq = "ATGAAACTG"
        assert "AUGAAACUG" == Seq.transcribe(seq)

    def test_seq_object_transcription_method(self):
        for nucleotide_seq in test_seqs:
            if isinstance(nucleotide_seq, Seq.Seq):
                assert repr(Seq.transcribe(nucleotide_seq)) == repr(nucleotide_seq.transcribe())

    def test_back_transcribe_rna_into_dna(self):
        for nucleotide_seq in test_seqs:
            expected = Seq.back_transcribe(nucleotide_seq)
            assert str(nucleotide_seq).replace("u", "t").replace("U", "T") == expected

    def test_back_transcribe_rna_string_into_dna(self):
        seq = "AUGAAACUG"
        assert "ATGAAACTG" == Seq.back_transcribe(seq)

    def test_seq_object_back_transcription_method(self):
        for nucleotide_seq in test_seqs:
            if isinstance(nucleotide_seq, Seq.Seq):
                expected = Seq.back_transcribe(nucleotide_seq)
                assert repr(nucleotide_seq.back_transcribe()) == repr(expected)


class TestTranslating(unittest.TestCase):
    def setUp(self):
        self.test_seqs = [
            Seq.Seq("TCAAAAGGATGCATCATG"),
            Seq.Seq("ATGAAACTG"),
            Seq.Seq("ATGAARCTG"),
            Seq.Seq("AWGAARCKG"),  # Note no U or T
            Seq.Seq("".join(ambiguous_rna_values)),
            Seq.Seq("".join(ambiguous_dna_values)),
            Seq.Seq("AUGAAACUG"),
            Seq.Seq("ATGAAACTGWN"),
            Seq.Seq("AUGAAACUGWN"),
            Seq.MutableSeq("ATGAAACTG"),
            Seq.MutableSeq("AUGaaaCUG"),
        ]

    def test_translation(self):
        for nucleotide_seq in self.test_seqs:
            nucleotide_seq = nucleotide_seq[: 3 * (len(nucleotide_seq) // 3)]
            if "X" not in nucleotide_seq:
                expected = Seq.translate(nucleotide_seq)
                assert expected == nucleotide_seq.translate()

    def test_gapped_seq_with_gap_char_given(self):
        seq = Seq.Seq("ATG---AAACTG")
        assert "M-KL" == seq.translate(gap="-")
        with pytest.raises(TranslationError):
            seq.translate(gap="~")

        seq = Seq.Seq("GTG---GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
        assert "V-AIVMGR*KGAR*" == seq.translate(gap="-")
        with pytest.raises(TranslationError):
            seq.translate(gap=None)

        seq = Seq.Seq("ATG~~~AAACTG")
        with pytest.raises(TranslationError):
            seq.translate(gap="-")

        seq = Seq.Seq("ATG---AAACTGTAG")
        assert "M-KL*" == seq.translate(gap="-")
        assert "M-KL@" == seq.translate(gap="-", stop_symbol="@")
        with pytest.raises(TranslationError):
            seq.translate(gap="~")

        seq = Seq.Seq("ATG~~~AAACTGTAG")
        with pytest.raises(TranslationError):
            seq.translate(gap="-")

    def test_gapped_seq_no_gap_char_given(self):
        seq = Seq.Seq("ATG---AAACTG")
        with pytest.raises(TranslationError):
            seq.translate(gap=None)

    def test_translation_wrong_type(self):
        """Test translation table cannot be CodonTable."""
        seq = Seq.Seq("ATCGTA")
        with pytest.raises(ValueError):
            seq.translate(table=ambiguous_dna_complement)

    def test_translation_of_string(self):
        seq = "GTGGCCATTGTAATGGGCCGC"
        assert "VAIVMGR" == Seq.translate(seq)

    def test_translation_of_gapped_string_with_gap_char_given(self):
        seq = "GTG---GCCATTGTAATGGGCCGC"
        expected = "V-AIVMGR"
        assert expected == Seq.translate(seq, gap="-")
        with pytest.raises(TypeError):
            Seq.translate(seq, gap=[])
        with pytest.raises(ValueError):
            Seq.translate(seq, gap="-*")

    def test_translation_of_gapped_string_no_gap_char_given(self):
        seq = "GTG---GCCATTGTAATGGGCCGC"
        with pytest.raises(TranslationError):
            Seq.translate(seq)

    def test_translation_to_stop(self):
        for nucleotide_seq in self.test_seqs:
            nucleotide_seq = nucleotide_seq[: 3 * (len(nucleotide_seq) // 3)]
            if "X" not in nucleotide_seq:
                short = Seq.translate(nucleotide_seq, to_stop=True)
                assert short == Seq.translate(nucleotide_seq).split("*")[0]

        seq = "GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
        assert "VAIVMGRWKGAR" == Seq.translate(seq, table=2, to_stop=True)

    def test_translation_on_proteins(self):
        """Check translation fails on a protein."""
        for s in protein_seqs:
            if len(s) % 3 != 0:
                with pytest.warns(BiopythonWarning):
                    with pytest.raises(TranslationError):
                        Seq.translate(s)

                with pytest.warns(BiopythonWarning):
                    with pytest.raises(TranslationError):
                        s.translate()
            else:
                with pytest.raises(TranslationError):
                    Seq.translate(s)

                with pytest.raises(TranslationError):
                    s.translate()

    def test_translation_of_invalid_codon(self):
        for codon in ["TA?", "N-N", "AC_", "Ac_"]:
            with pytest.raises(TranslationError):
                Seq.translate(codon)

    def test_translation_of_glutamine(self):
        for codon in ["SAR", "SAG", "SAA"]:
            assert "Z" == Seq.translate(codon)

    def test_translation_of_asparagine(self):
        for codon in ["RAY", "RAT", "RAC"]:
            assert "B" == Seq.translate(codon)

    def test_translation_of_leucine(self):
        for codon in ["WTA", "MTY", "MTT", "MTW", "MTM", "MTH", "MTA", "MTC", "HTA"]:
            assert "J" == Seq.translate(codon)

    def test_translation_with_bad_table_argument(self):
        table = {}
        with pytest.raises(ValueError) as cm:
            Seq.translate("GTGGCCATTGTAATGGGCCGC", table=table)
        assert str(cm.value) == "Bad table argument"
        table = b"0x"
        with pytest.raises(TypeError) as cm:
            Seq.translate("GTGGCCATTGTAATGGGCCGC", table=table)
        assert str(cm.value) == "table argument must be integer or string"

    def test_translation_with_codon_table_as_table_argument(self):
        table = standard_dna_table
        assert "VAIVMGR" == Seq.translate("GTGGCCATTGTAATGGGCCGC", table=table)

    def test_translation_incomplete_codon(self):
        with pytest.warns(BiopythonWarning):
            Seq.translate("GTGGCCATTGTAATGGGCCG")

    def test_translation_extra_stop_codon(self):
        seq = "GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGTAG"
        with pytest.raises(TranslationError):
            Seq.translate(seq, table=2, cds=True)

    def test_translation_using_cds(self):
        seq = "GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
        assert "MAIVMGRWKGAR" == Seq.translate(seq, table=2, cds=True)

        seq = "GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCG"  # not multiple of three
        with pytest.raises(TranslationError):
            Seq.translate(seq, table=2, cds=True)

        seq = "GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA"  # no stop codon
        with pytest.raises(TranslationError):
            Seq.translate(seq, table=2, cds=True)

        seq = "GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"  # no start codon
        with pytest.raises(TranslationError):
            Seq.translate(seq, table=2, cds=True)

    def test_translation_using_tables_with_ambiguous_stop_codons(self):
        """Check for error and warning messages.

        Here, 'ambiguous stop codons' means codons of unambiguous sequence
        but with a context sensitive encoding as STOP or an amino acid.
        Thus, these codons appear within the codon table in the forward
        table as well as in the list of stop codons.
        """
        seq = "ATGGGCTGA"
        with pytest.raises(ValueError):
            Seq.translate(seq, table=28, to_stop=True)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            Seq.translate(seq, table=28)
            message = str(w[-1].message)
            assert message.startswith("This table contains")
            assert message.endswith("be translated as amino acid.")


class TestStopCodons(unittest.TestCase):
    def setUp(self):
        self.misc_stops = "TAATAGTGAAGAAGG"

    def test_stops(self):
        for nucleotide_seq in [self.misc_stops, Seq.Seq(self.misc_stops)]:
            assert "***RR" == Seq.translate(nucleotide_seq)
            assert "***RR" == Seq.translate(nucleotide_seq, table=1)
            assert "***RR" == Seq.translate(nucleotide_seq, table="SGC0")
            assert "**W**" == Seq.translate(nucleotide_seq, table=2)
            assert "**WRR" == Seq.translate(nucleotide_seq, table="Yeast Mitochondrial")
            assert "**WSS" == Seq.translate(nucleotide_seq, table=5)
            assert "**WSS" == Seq.translate(nucleotide_seq, table=9)
            assert "**CRR" == Seq.translate(nucleotide_seq, table="Euplotid Nuclear")
            assert "***RR" == Seq.translate(nucleotide_seq, table=11)
            assert "***RR" == Seq.translate(nucleotide_seq, table="Bacterial")

    def test_translation_of_stops(self):
        assert Seq.translate("TAT") == "Y"
        assert Seq.translate("TAR") == "*"
        assert Seq.translate("TAN") == "X"
        assert Seq.translate("NNN") == "X"

        assert Seq.translate("TAt") == "Y"
        assert Seq.translate("TaR") == "*"
        assert Seq.translate("TaN") == "X"
        assert Seq.translate("nnN") == "X"

        assert Seq.translate("tat") == "Y"
        assert Seq.translate("tar") == "*"
        assert Seq.translate("tan") == "X"
        assert Seq.translate("nnn") == "X"


class TestAttributes(unittest.TestCase):
    def test_seq(self):
        s = Seq.Seq("ACGT")
        with pytest.raises(AttributeError):
            s.dog
        s.dog = "woof"
        assert "dog" in dir(s)
        assert s.dog == "woof"
        del s.dog
        with pytest.raises(AttributeError):
            s.dog
        assert "dog" not in dir(s)
        with pytest.raises(AttributeError):
            s.cat
        s.dog = "woof"
        s.cat = "meow"
        assert "dog" in dir(s)
        assert "cat" in dir(s)
        assert s.dog == "woof"
        assert s.cat == "meow"
        del s.dog
        with pytest.raises(AttributeError):
            s.dog
        assert "dog" not in dir(s)
        assert "cat" in dir(s)
        assert s.cat == "meow"
        del s.cat
        with pytest.raises(AttributeError):
            s.cat
        assert "cat" not in dir(s)
        s.dog = "woof"
        s.dog = "bark"
        assert "dog" in dir(s)
        assert s.dog == "bark"
        del s.dog
        with pytest.raises(AttributeError):
            s.dog
        assert "dog" not in dir(s)

    def test_mutable_seq(self):
        s = Seq.MutableSeq("ACGT")
        with pytest.raises(AttributeError):
            s.dog
        s.dog = "woof"
        assert "dog" in dir(s)
        assert s.dog == "woof"
        del s.dog
        with pytest.raises(AttributeError):
            s.dog
        assert "dog" not in dir(s)
        with pytest.raises(AttributeError):
            s.cat
        s.dog = "woof"
        s.cat = "meow"
        assert "dog" in dir(s)
        assert "cat" in dir(s)
        assert s.dog == "woof"
        assert s.cat == "meow"
        del s.dog
        with pytest.raises(AttributeError):
            s.dog
        assert "dog" not in dir(s)
        assert "cat" in dir(s)
        assert s.cat == "meow"
        del s.cat
        with pytest.raises(AttributeError):
            s.cat
        assert "cat" not in dir(s)
        s.dog = "woof"
        s.dog = "bark"
        assert "dog" in dir(s)
        assert s.dog == "bark"
        del s.dog
        with pytest.raises(AttributeError):
            s.dog
        assert "dog" not in dir(s)


class TestSeqDefined(unittest.TestCase):
    def test_zero_length(self):
        zero_length_seqs = [
            Seq.Seq(""),
            Seq.Seq(None, length=0),
            Seq.Seq({}, length=0),
            Seq.MutableSeq(""),
        ]

        for seq in zero_length_seqs:
            assert seq.defined, repr(seq)
            assert seq.defined_ranges == (), repr(seq)

    def test_undefined(self):
        seq = Seq.Seq(None, length=1)
        assert not seq.defined
        assert seq.defined_ranges == ()
        seq = Seq.Seq({3: "ACGT"}, length=10)
        assert not seq.defined
        assert seq.defined_ranges == ((3, 7),)

    def test_defined(self):
        seqs = [
            Seq.Seq("T"),
            Seq.Seq({0: "A"}, length=1),
            Seq.Seq({0: "A", 1: "C"}, length=2),
        ]

        for seq in seqs:
            assert seq.defined, repr(seq)
            assert seq.defined_ranges == ((0, len(seq)),), repr(seq)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
