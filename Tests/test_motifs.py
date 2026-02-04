# Copyright 2008 by Bartek Wilczynski.  All rights reserved.
# Revisions copyright 2019 by Victor Lin.
# Adapted from test_Mymodule.py by Jeff Chang.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Tests for motifs module."""

import math
import tempfile
import unittest
import pytest

try:
    import numpy as np
except ImportError:

    pytest.skip("Install numpy if you want to use Bio.motifs.", allow_module_level=True)

from Bio import motifs
from Bio.Seq import Seq


class TestBasic(unittest.TestCase):
    """Basic motif tests."""

    def test_format(self):
        m = motifs.create([Seq("ATATA")])
        m.name = "Foo"
        s1 = format(m, "pfm")
        expected_pfm = """  1.00   0.00   1.00   0.00  1.00
  0.00   0.00   0.00   0.00  0.00
  0.00   0.00   0.00   0.00  0.00
  0.00   1.00   0.00   1.00  0.00
"""
        s2 = format(m, "jaspar")
        expected_jaspar = """>None Foo
A [  1.00   0.00   1.00   0.00   1.00]
C [  0.00   0.00   0.00   0.00   0.00]
G [  0.00   0.00   0.00   0.00   0.00]
T [  0.00   1.00   0.00   1.00   0.00]
"""
        assert s2 == expected_jaspar
        s3 = format(m, "transfac")
        expected_transfac = """P0      A      C      G      T
01      1      0      0      0      A
02      0      0      0      1      T
03      1      0      0      0      A
04      0      0      0      1      T
05      1      0      0      0      A
XX
//
"""
        assert s3 == expected_transfac
        with pytest.raises(ValueError):
            format(m, "foo_bar")

    def test_relative_entropy(self):
        m = motifs.create([Seq("ATATA"), Seq("ATCTA"), Seq("TTGTA")])
        assert len(m.alignment) == 3
        assert m.background == {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}
        assert m.pseudocounts == {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0}
        assert np.allclose(
                m.relative_entropy,
                np.array([1.0817041659455104, 2.0, 0.4150374992788437, 2.0, 2.0]),
            )
        m.background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
        assert np.allclose(
                m.relative_entropy,
                np.array(
                    [
                        0.8186697601117167,
                        1.7369655941662063,
                        0.5419780939258206,
                        1.7369655941662063,
                        1.7369655941662063,
                    ]
                ),
            )
        m.background = None
        assert m.background == {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}
        pseudocounts = math.sqrt(len(m.alignment))
        m.pseudocounts = {
            letter: m.background[letter] * pseudocounts for letter in "ACGT"
        }
        assert np.allclose(
                m.relative_entropy,
                np.array(
                    [
                        0.3532586861097656,
                        0.7170228827697498,
                        0.11859369972847714,
                        0.7170228827697498,
                        0.7170228827697499,
                    ]
                ),
            )
        m.background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
        assert np.allclose(
                m.relative_entropy,
                np.array(
                    [
                        0.19727984803857979,
                        0.561044044698564,
                        0.20984910512125132,
                        0.561044044698564,
                        0.5610440446985638,
                    ]
                ),
            )

    def test_reverse_complement(self):
        """Test if motifs can be reverse-complemented."""
        background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
        pseudocounts = 0.5
        m = motifs.create([Seq("ATATA")])
        m.background = background
        m.pseudocounts = pseudocounts
        received_forward = format(m, "transfac")
        expected_forward = """\
P0      A      C      G      T
01      1      0      0      0      A
02      0      0      0      1      T
03      1      0      0      0      A
04      0      0      0      1      T
05      1      0      0      0      A
XX
//
"""
        assert received_forward == expected_forward
        expected_forward_pwm = """\
        0      1      2      3      4
A:   0.50   0.17   0.50   0.17   0.50
C:   0.17   0.17   0.17   0.17   0.17
G:   0.17   0.17   0.17   0.17   0.17
T:   0.17   0.50   0.17   0.50   0.17
"""
        assert str(m.pwm) == expected_forward_pwm
        m = m.reverse_complement()
        received_reverse = format(m, "transfac")
        expected_reverse = """\
P0      A      C      G      T
01      0      0      0      1      T
02      1      0      0      0      A
03      0      0      0      1      T
04      1      0      0      0      A
05      0      0      0      1      T
XX
//
"""
        assert received_reverse == expected_reverse
        expected_reverse_pwm = """\
        0      1      2      3      4
A:   0.17   0.50   0.17   0.50   0.17
C:   0.17   0.17   0.17   0.17   0.17
G:   0.17   0.17   0.17   0.17   0.17
T:   0.50   0.17   0.50   0.17   0.50
"""
        assert str(m.pwm) == expected_reverse_pwm
        # Same but for RNA motif.
        background_rna = {"A": 0.3, "C": 0.2, "G": 0.2, "U": 0.3}
        pseudocounts = 0.5
        m_rna = motifs.create([Seq("AUAUA")], alphabet="ACGU")
        m_rna.background = background_rna
        m_rna.pseudocounts = pseudocounts
        expected_forward_rna_counts = """\
        0      1      2      3      4
A:   1.00   0.00   1.00   0.00   1.00
C:   0.00   0.00   0.00   0.00   0.00
G:   0.00   0.00   0.00   0.00   0.00
U:   0.00   1.00   0.00   1.00   0.00
"""
        assert str(m_rna.counts) == expected_forward_rna_counts
        expected_forward_rna_pwm = """\
        0      1      2      3      4
A:   0.50   0.17   0.50   0.17   0.50
C:   0.17   0.17   0.17   0.17   0.17
G:   0.17   0.17   0.17   0.17   0.17
U:   0.17   0.50   0.17   0.50   0.17
"""
        assert str(m_rna.pwm) == expected_forward_rna_pwm
        expected_reverse_rna_counts = """\
        0      1      2      3      4
A:   0.00   1.00   0.00   1.00   0.00
C:   0.00   0.00   0.00   0.00   0.00
G:   0.00   0.00   0.00   0.00   0.00
U:   1.00   0.00   1.00   0.00   1.00
"""
        assert str(m_rna.reverse_complement().counts) == expected_reverse_rna_counts
        expected_reverse_rna_pwm = """\
        0      1      2      3      4
A:   0.17   0.50   0.17   0.50   0.17
C:   0.17   0.17   0.17   0.17   0.17
G:   0.17   0.17   0.17   0.17   0.17
U:   0.50   0.17   0.50   0.17   0.50
"""
        assert str(m_rna.reverse_complement().pwm) == expected_reverse_rna_pwm
        # Same thing, but now start with a motif calculated from a count matrix
        m = motifs.create([Seq("ATATA")])
        counts = m.counts
        m = motifs.Motif(counts=counts)
        m.background = background
        m.pseudocounts = pseudocounts
        received_forward = format(m, "transfac")
        assert received_forward == expected_forward
        assert str(m.pwm) == expected_forward_pwm
        m = m.reverse_complement()
        received_reverse = format(m, "transfac")
        assert received_reverse == expected_reverse
        assert str(m.pwm) == expected_reverse_pwm
        # Same, but for RNA count matrix
        m_rna = motifs.create([Seq("AUAUA")], alphabet="ACGU")
        counts = m_rna.counts
        m_rna = motifs.Motif(counts=counts, alphabet="ACGU")
        m_rna.background = background_rna
        m_rna.pseudocounts = pseudocounts
        assert str(m_rna.counts) == expected_forward_rna_counts
        assert str(m_rna.pwm) == expected_forward_rna_pwm
        assert str(m_rna.reverse_complement().counts) == expected_reverse_rna_counts
        assert str(m_rna.reverse_complement().pwm) == expected_reverse_rna_pwm


class TestAlignAce(unittest.TestCase):
    """Testing parsing AlignAce output files."""

    def test_alignace_parsing(self):
        """Test if Bio.motifs can parse AlignAce output files."""
        with open("motifs/alignace.out") as stream:
            record = motifs.parse(stream, "AlignAce")
        assert record.version == "AlignACE 4.0 05/13/04"
        assert record.command == "./AlignACE -i test.fa"
        assert len(record.parameters) == 7
        assert record.parameters["expect"] == "10"
        assert record.parameters["gcback"] == "0.38"
        assert record.parameters["minpass"] == "200"
        assert record.parameters["seed"] == "1227623309"
        assert record.parameters["numcols"] == "10"
        assert record.parameters["undersample"] == "1"
        assert record.parameters["oversample"] == "1"
        assert len(record.sequences) == 10
        assert record.sequences[0] == "SEQ1; M: CTCAATCGTAGA at 52"
        assert record.sequences[1] == "SEQ2; M: CTCAATCGTAGA at 172"
        assert record.sequences[2] == "SEQ3; M: CTCAATCGTAGA at 112"
        assert record.sequences[3] == "SEQ4; M: CTCAATCGTAGA at 173"
        assert record.sequences[4] == "SEQ5; M: CTCAATCGTAGA at 185"
        assert record.sequences[5] == "SEQ6; M: CTCAATCGTAGA at 105"
        assert record.sequences[6] == "SEQ7; M: CTCAATCGTAGA at 177"
        assert record.sequences[7] == "SEQ8; M: CTCAATCGTAGA at 172"
        assert record.sequences[8] == "SEQ9; M: CTCAATCGTAGA at 93"
        assert record.sequences[9] == "SEQ10; M: CTCAATCGTAGA at 3"
        assert len(record) == 16
        assert record[0].alphabet == "ACGT"
        assert len(record[0].alignment.sequences) == 11
        assert record[0].alignment.sequences[0] == "TCTACGATTGAG"
        assert record[0].alignment.sequences[1] == "TCTACGATTGAG"
        assert record[0].alignment.sequences[2] == "TCTACGATTGAG"
        assert record[0].alignment.sequences[3] == "TCTACGATTGAG"
        assert record[0].alignment.sequences[4] == "TCTACGATTGAG"
        assert record[0].alignment.sequences[5] == "TCTACGATTGAG"
        assert record[0].alignment.sequences[6] == "TCTACGATTGAG"
        assert record[0].alignment.sequences[7] == "TCTACGATTGAG"
        assert record[0].alignment.sequences[8] == "TCTACGATTGAG"
        assert record[0].alignment.sequences[9] == "TCAAAGATAGAG"
        assert record[0].alignment.sequences[10] == "TCTACGATTGAG"
        assert record[0].mask == (1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1)
        assert record[0].score == pytest.approx(57.9079, abs=5e-8)
        assert str(record[0]) == """\
TCTACGATTGAG
TCTACGATTGAG
TCTACGATTGAG
TCTACGATTGAG
TCTACGATTGAG
TCTACGATTGAG
TCTACGATTGAG
TCTACGATTGAG
TCTACGATTGAG
TCAAAGATAGAG
TCTACGATTGAG"""
        motif = record[0][2:-1]
        assert motif.alphabet == "ACGT"
        assert str(motif) == """\
TACGATTGA
TACGATTGA
TACGATTGA
TACGATTGA
TACGATTGA
TACGATTGA
TACGATTGA
TACGATTGA
TACGATTGA
AAAGATAGA
TACGATTGA"""
        assert motif.mask == (0, 1, 1, 1, 1, 1, 0, 1, 1)
        assert record[1].alphabet == "ACGT"
        assert len(record[1].alignment.sequences) == 22
        assert record[1].alignment.sequences[0] == "GCGAAGGAAGCAGCGCGTGTG"
        assert record[1].alignment.sequences[1] == "GGCACCGCCTCTACGATTGAG"
        assert record[1].alignment.sequences[2] == "CAGAGCTTAGCATTGAACGCG"
        assert record[1].alignment.sequences[3] == "CTAATGAAAGCAATGAGAGTG"
        assert record[1].alignment.sequences[4] == "CTTGTGCCCTCTAAGCGTCCG"
        assert record[1].alignment.sequences[5] == "GAGCACGACGCTTTGTACCTG"
        assert record[1].alignment.sequences[6] == "CGGCACTTAGCAGCGTATCGT"
        assert record[1].alignment.sequences[7] == "CTGGTTTCATCTACGATTGAG"
        assert record[1].alignment.sequences[8] == "GGGCCAATAGCGGCGCCGGAG"
        assert record[1].alignment.sequences[9] == "GTGGAGTTATCTTAGTGCGCG"
        assert record[1].alignment.sequences[10] == "GAGAGGTTATCTACGATTGAG"
        assert record[1].alignment.sequences[11] == "CTGCTCCCCGCATACAGCGCG"
        assert record[1].alignment.sequences[12] == "CAGAACCGAGGTCCGGTACGG"
        assert record[1].alignment.sequences[13] == "GTGCCCCAAGCTTACCCAGGG"
        assert record[1].alignment.sequences[14] == "CGCCTCTGATCTACGATTGAG"
        assert record[1].alignment.sequences[15] == "GTGCTCATAGGGACGTCGCGG"
        assert record[1].alignment.sequences[16] == "CTGCCCCCCGCATAGTAGGGG"
        assert record[1].alignment.sequences[17] == "GTAAAGAAATCGATGTGCCAG"
        assert record[1].alignment.sequences[18] == "CACCTGCAATTGCTGGCAGCG"
        assert record[1].alignment.sequences[19] == "GGCGGGCCATCCCTGTATGAA"
        assert record[1].alignment.sequences[20] == "CTCCAGGTCGCATGGAGAGAG"
        assert record[1].alignment.sequences[21] == "CCTCGGATCGCTTGGGAAGAG"
        assert record[1].mask == (1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1)
        assert record[1].score == pytest.approx(19.6235, abs=5e-8)
        assert str(record[1]) == """\
GCGAAGGAAGCAGCGCGTGTG
GGCACCGCCTCTACGATTGAG
CAGAGCTTAGCATTGAACGCG
CTAATGAAAGCAATGAGAGTG
CTTGTGCCCTCTAAGCGTCCG
GAGCACGACGCTTTGTACCTG
CGGCACTTAGCAGCGTATCGT
CTGGTTTCATCTACGATTGAG
GGGCCAATAGCGGCGCCGGAG
GTGGAGTTATCTTAGTGCGCG
GAGAGGTTATCTACGATTGAG
CTGCTCCCCGCATACAGCGCG
CAGAACCGAGGTCCGGTACGG
GTGCCCCAAGCTTACCCAGGG
CGCCTCTGATCTACGATTGAG
GTGCTCATAGGGACGTCGCGG
CTGCCCCCCGCATAGTAGGGG
GTAAAGAAATCGATGTGCCAG
CACCTGCAATTGCTGGCAGCG
GGCGGGCCATCCCTGTATGAA
CTCCAGGTCGCATGGAGAGAG
CCTCGGATCGCTTGGGAAGAG"""
        motif = record[1][2:-1]
        assert motif.alphabet == "ACGT"
        assert str(motif) == """\
GAAGGAAGCAGCGCGTGT
CACCGCCTCTACGATTGA
GAGCTTAGCATTGAACGC
AATGAAAGCAATGAGAGT
TGTGCCCTCTAAGCGTCC
GCACGACGCTTTGTACCT
GCACTTAGCAGCGTATCG
GGTTTCATCTACGATTGA
GCCAATAGCGGCGCCGGA
GGAGTTATCTTAGTGCGC
GAGGTTATCTACGATTGA
GCTCCCCGCATACAGCGC
GAACCGAGGTCCGGTACG
GCCCCAAGCTTACCCAGG
CCTCTGATCTACGATTGA
GCTCATAGGGACGTCGCG
GCCCCCCGCATAGTAGGG
AAAGAAATCGATGTGCCA
CCTGCAATTGCTGGCAGC
CGGGCCATCCCTGTATGA
CCAGGTCGCATGGAGAGA
TCGGATCGCTTGGGAAGA"""

        assert record[2].alphabet == "ACGT"
        assert len(record[2].alignment.sequences) == 18
        assert record[2].alignment.sequences[0] == "GTGCGCGAAGGAAGCAGCGCG"
        assert record[2].alignment.sequences[1] == "CAGAGCTTAGCATTGAACGCG"
        assert record[2].alignment.sequences[2] == "GTGCCCGATGACCACCCGTCG"
        assert record[2].alignment.sequences[3] == "GCCCTCTAAGCGTCCGCGGAT"
        assert record[2].alignment.sequences[4] == "GAGCACGACGCTTTGTACCTG"
        assert record[2].alignment.sequences[5] == "CGGCACTTAGCAGCGTATCGT"
        assert record[2].alignment.sequences[6] == "GGGCCAATAGCGGCGCCGGAG"
        assert record[2].alignment.sequences[7] == "GCGCACTAAGATAACTCCACG"
        assert record[2].alignment.sequences[8] == "CGGCCCGTTGTCCAGCAGACG"
        assert record[2].alignment.sequences[9] == "CTGCTCCCCGCATACAGCGCG"
        assert record[2].alignment.sequences[10] == "GTGCCCCAAGCTTACCCAGGG"
        assert record[2].alignment.sequences[11] == "GTGCTCATAGGGACGTCGCGG"
        assert record[2].alignment.sequences[12] == "CTGCCCCCCGCATAGTAGGGG"
        assert record[2].alignment.sequences[13] == "CGCCGCCATGCGACGCAGAGG"
        assert record[2].alignment.sequences[14] == "AACCTCTAAGCATACTCTACG"
        assert record[2].alignment.sequences[15] == "GACCTGGAGGCTTAGACTTGG"
        assert record[2].alignment.sequences[16] == "GCGCTCTTCCCAAGCGATCCG"
        assert record[2].alignment.sequences[17] == "GGGCCGTCAGCTCTCAAGTCT"
        assert record[2].mask == (1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1)
        assert record[2].score == pytest.approx(19.1804, abs=5e-8)
        assert str(record[2]) == """\
GTGCGCGAAGGAAGCAGCGCG
CAGAGCTTAGCATTGAACGCG
GTGCCCGATGACCACCCGTCG
GCCCTCTAAGCGTCCGCGGAT
GAGCACGACGCTTTGTACCTG
CGGCACTTAGCAGCGTATCGT
GGGCCAATAGCGGCGCCGGAG
GCGCACTAAGATAACTCCACG
CGGCCCGTTGTCCAGCAGACG
CTGCTCCCCGCATACAGCGCG
GTGCCCCAAGCTTACCCAGGG
GTGCTCATAGGGACGTCGCGG
CTGCCCCCCGCATAGTAGGGG
CGCCGCCATGCGACGCAGAGG
AACCTCTAAGCATACTCTACG
GACCTGGAGGCTTAGACTTGG
GCGCTCTTCCCAAGCGATCCG
GGGCCGTCAGCTCTCAAGTCT"""
        motif = record[2][2:-1]
        assert motif.alphabet == "ACGT"
        assert str(motif) == """\
GCGCGAAGGAAGCAGCGC
GAGCTTAGCATTGAACGC
GCCCGATGACCACCCGTC
CCTCTAAGCGTCCGCGGA
GCACGACGCTTTGTACCT
GCACTTAGCAGCGTATCG
GCCAATAGCGGCGCCGGA
GCACTAAGATAACTCCAC
GCCCGTTGTCCAGCAGAC
GCTCCCCGCATACAGCGC
GCCCCAAGCTTACCCAGG
GCTCATAGGGACGTCGCG
GCCCCCCGCATAGTAGGG
CCGCCATGCGACGCAGAG
CCTCTAAGCATACTCTAC
CCTGGAGGCTTAGACTTG
GCTCTTCCCAAGCGATCC
GCCGTCAGCTCTCAAGTC"""

        assert record[3].alphabet == "ACGT"
        assert len(record[3].alignment.sequences) == 16
        assert record[3].alignment.sequences[0] == "GCCCCAAGCTTACCCAGGGAC"
        assert record[3].alignment.sequences[1] == "GCCGTCTGCTGGACAACGGGC"
        assert record[3].alignment.sequences[2] == "GCCGACGGGTGGTCATCGGGC"
        assert record[3].alignment.sequences[3] == "GCCAATAGCGGCGCCGGAGTC"
        assert record[3].alignment.sequences[4] == "GCCCCCCGCATAGTAGGGGGA"
        assert record[3].alignment.sequences[5] == "GCCCGTACCGGACCTCGGTTC"
        assert record[3].alignment.sequences[6] == "GCCTCATGTACCGGAAGGGAC"
        assert record[3].alignment.sequences[7] == "GACACGCGCCTGGGAGGGTTC"
        assert record[3].alignment.sequences[8] == "GCCTTTGGCCTTGGATGAGAA"
        assert record[3].alignment.sequences[9] == "GGCCCTCGGATCGCTTGGGAA"
        assert record[3].alignment.sequences[10] == "GCATGTTGGGAATCCGCGGAC"
        assert record[3].alignment.sequences[11] == "GACACGCGCTGTATGCGGGGA"
        assert record[3].alignment.sequences[12] == "GCCAGGTACAAAGCGTCGTGC"
        assert record[3].alignment.sequences[13] == "GCGATCAGCTTGTGGGCGTGC"
        assert record[3].alignment.sequences[14] == "GACAAATCGGATACTGGGGCA"
        assert record[3].alignment.sequences[15] == "GCACTTAGCAGCGTATCGTTA"
        assert record[3].mask == (1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1)
        assert record[3].score == pytest.approx(18.0097, abs=5e-8)
        assert str(record[3]) == """\
GCCCCAAGCTTACCCAGGGAC
GCCGTCTGCTGGACAACGGGC
GCCGACGGGTGGTCATCGGGC
GCCAATAGCGGCGCCGGAGTC
GCCCCCCGCATAGTAGGGGGA
GCCCGTACCGGACCTCGGTTC
GCCTCATGTACCGGAAGGGAC
GACACGCGCCTGGGAGGGTTC
GCCTTTGGCCTTGGATGAGAA
GGCCCTCGGATCGCTTGGGAA
GCATGTTGGGAATCCGCGGAC
GACACGCGCTGTATGCGGGGA
GCCAGGTACAAAGCGTCGTGC
GCGATCAGCTTGTGGGCGTGC
GACAAATCGGATACTGGGGCA
GCACTTAGCAGCGTATCGTTA"""
        motif = record[3][2:-1]
        assert motif.alphabet == "ACGT"
        assert str(motif) == """\
CCCAAGCTTACCCAGGGA
CGTCTGCTGGACAACGGG
CGACGGGTGGTCATCGGG
CAATAGCGGCGCCGGAGT
CCCCCGCATAGTAGGGGG
CCGTACCGGACCTCGGTT
CTCATGTACCGGAAGGGA
CACGCGCCTGGGAGGGTT
CTTTGGCCTTGGATGAGA
CCCTCGGATCGCTTGGGA
ATGTTGGGAATCCGCGGA
CACGCGCTGTATGCGGGG
CAGGTACAAAGCGTCGTG
GATCAGCTTGTGGGCGTG
CAAATCGGATACTGGGGC
ACTTAGCAGCGTATCGTT"""
        assert record[4].alphabet == "ACGT"
        assert len(record[4].alignment.sequences) == 15
        assert record[4].alignment.sequences[0] == "CGGCACAGAGCTT"
        assert record[4].alignment.sequences[1] == "ATCCGCGGACGCT"
        assert record[4].alignment.sequences[2] == "CGCCTGGGAGGGT"
        assert record[4].alignment.sequences[3] == "CGGAAGGGACGTT"
        assert record[4].alignment.sequences[4] == "ACACACAGACGGT"
        assert record[4].alignment.sequences[5] == "TGCCAGAGAGGTT"
        assert record[4].alignment.sequences[6] == "AGACTGAGACGTT"
        assert record[4].alignment.sequences[7] == "AATCGTAGAGGAT"
        assert record[4].alignment.sequences[8] == "CGTCTCGTAGGGT"
        assert record[4].alignment.sequences[9] == "CGTCGCGGAGGAT"
        assert record[4].alignment.sequences[10] == "CTTCTTAGACGCT"
        assert record[4].alignment.sequences[11] == "CGACGCAGAGGAT"
        assert record[4].alignment.sequences[12] == "ATGCTTAGAGGTT"
        assert record[4].alignment.sequences[13] == "AGACTTGGGCGAT"
        assert record[4].alignment.sequences[14] == "CGACCTGGAGGCT"
        assert record[4].mask == (1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1)
        assert record[4].score == pytest.approx(16.8287, abs=5e-8)
        assert str(record[4]) == """\
CGGCACAGAGCTT
ATCCGCGGACGCT
CGCCTGGGAGGGT
CGGAAGGGACGTT
ACACACAGACGGT
TGCCAGAGAGGTT
AGACTGAGACGTT
AATCGTAGAGGAT
CGTCTCGTAGGGT
CGTCGCGGAGGAT
CTTCTTAGACGCT
CGACGCAGAGGAT
ATGCTTAGAGGTT
AGACTTGGGCGAT
CGACCTGGAGGCT"""
        motif = record[4][2:-1]
        assert motif.alphabet == "ACGT"
        assert str(motif) == """\
GCACAGAGCT
CCGCGGACGC
CCTGGGAGGG
GAAGGGACGT
ACACAGACGG
CCAGAGAGGT
ACTGAGACGT
TCGTAGAGGA
TCTCGTAGGG
TCGCGGAGGA
TCTTAGACGC
ACGCAGAGGA
GCTTAGAGGT
ACTTGGGCGA
ACCTGGAGGC"""
        assert record[5].alphabet == "ACGT"
        assert len(record[5].alignment.sequences) == 18
        assert record[5].alignment.sequences[0] == "GTGCGCGAAGGAAGCAGCGCGTG"
        assert record[5].alignment.sequences[1] == "TTGAGCCGAGTAAAGGGCTGGTG"
        assert record[5].alignment.sequences[2] == "CAATGCTAAGCTCTGTGCCGACG"
        assert record[5].alignment.sequences[3] == "CAACTCTCTATGTAGTGCCCGAG"
        assert record[5].alignment.sequences[4] == "CGACGCTTTGTACCTGGCTTGCG"
        assert record[5].alignment.sequences[5] == "CGAGTCAATGACACGCGCCTGGG"
        assert record[5].alignment.sequences[6] == "CGATACGCTGCTAAGTGCCGTCC"
        assert record[5].alignment.sequences[7] == "CCGGGCCAATAGCGGCGCCGGAG"
        assert record[5].alignment.sequences[8] == "CCACGCTTCGACACGTGGTATAG"
        assert record[5].alignment.sequences[9] == "CCGAGCCTCATGTACCGGAAGGG"
        assert record[5].alignment.sequences[10] == "CTGCTCCCCGCATACAGCGCGTG"
        assert record[5].alignment.sequences[11] == "CCGAGGTCCGGTACGGGCAAGCC"
        assert record[5].alignment.sequences[12] == "GTGCTCATAGGGACGTCGCGGAG"
        assert record[5].alignment.sequences[13] == "CCCTACTATGCGGGGGGCAGGTC"
        assert record[5].alignment.sequences[14] == "GCCAGCAATTGCAGGTGGTCGTG"
        assert record[5].alignment.sequences[15] == "CTCTGCGTCGCATGGCGGCGTGG"
        assert record[5].alignment.sequences[16] == "GGAGGCTTAGACTTGGGCGATAC"
        assert record[5].alignment.sequences[17] == "GCATGGAGAGAGATCCGGAGGAG"
        assert record[5].mask == (1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1)
        assert record[5].score == pytest.approx(15.0441, abs=5e-8)
        assert str(record[5]) == """\
GTGCGCGAAGGAAGCAGCGCGTG
TTGAGCCGAGTAAAGGGCTGGTG
CAATGCTAAGCTCTGTGCCGACG
CAACTCTCTATGTAGTGCCCGAG
CGACGCTTTGTACCTGGCTTGCG
CGAGTCAATGACACGCGCCTGGG
CGATACGCTGCTAAGTGCCGTCC
CCGGGCCAATAGCGGCGCCGGAG
CCACGCTTCGACACGTGGTATAG
CCGAGCCTCATGTACCGGAAGGG
CTGCTCCCCGCATACAGCGCGTG
CCGAGGTCCGGTACGGGCAAGCC
GTGCTCATAGGGACGTCGCGGAG
CCCTACTATGCGGGGGGCAGGTC
GCCAGCAATTGCAGGTGGTCGTG
CTCTGCGTCGCATGGCGGCGTGG
GGAGGCTTAGACTTGGGCGATAC
GCATGGAGAGAGATCCGGAGGAG"""
        motif = record[5][2:-1]
        assert motif.alphabet == "ACGT"
        assert str(motif) == """\
GCGCGAAGGAAGCAGCGCGT
GAGCCGAGTAAAGGGCTGGT
ATGCTAAGCTCTGTGCCGAC
ACTCTCTATGTAGTGCCCGA
ACGCTTTGTACCTGGCTTGC
AGTCAATGACACGCGCCTGG
ATACGCTGCTAAGTGCCGTC
GGGCCAATAGCGGCGCCGGA
ACGCTTCGACACGTGGTATA
GAGCCTCATGTACCGGAAGG
GCTCCCCGCATACAGCGCGT
GAGGTCCGGTACGGGCAAGC
GCTCATAGGGACGTCGCGGA
CTACTATGCGGGGGGCAGGT
CAGCAATTGCAGGTGGTCGT
CTGCGTCGCATGGCGGCGTG
AGGCTTAGACTTGGGCGATA
ATGGAGAGAGATCCGGAGGA"""
        assert record[6].alphabet == "ACGT"
        assert len(record[6].alignment.sequences) == 20
        assert record[6].alignment.sequences[0] == "GCGCGTGTGTGTAAC"
        assert record[6].alignment.sequences[1] == "GCACAGAGCTTAGCA"
        assert record[6].alignment.sequences[2] == "GGTGGTCATCGGGCA"
        assert record[6].alignment.sequences[3] == "GCGCGTGTCATTGAC"
        assert record[6].alignment.sequences[4] == "GGACGGCACTTAGCA"
        assert record[6].alignment.sequences[5] == "GCGCGTCCCGGGCCA"
        assert record[6].alignment.sequences[6] == "GCTCGGCCCGTTGTC"
        assert record[6].alignment.sequences[7] == "GCGCGTGTCCTTTAA"
        assert record[6].alignment.sequences[8] == "GCTGATCGCTGCTCC"
        assert record[6].alignment.sequences[9] == "GCCCGTACCGGACCT"
        assert record[6].alignment.sequences[10] == "GGACGTCGCGGAGGA"
        assert record[6].alignment.sequences[11] == "GCGGGGGGCAGGTCA"
        assert record[6].alignment.sequences[12] == "GGACGTACTGGCACA"
        assert record[6].alignment.sequences[13] == "GCAGGTGGTCGTGCA"
        assert record[6].alignment.sequences[14] == "GCGCATACCTTAACA"
        assert record[6].alignment.sequences[15] == "GCACGGGACTTCAAC"
        assert record[6].alignment.sequences[16] == "GCACGTAGCTGGTAA"
        assert record[6].alignment.sequences[17] == "GCTCGTCTATGGTCA"
        assert record[6].alignment.sequences[18] == "GCGCATGCTGGATCC"
        assert record[6].alignment.sequences[19] == "GGCCGTCAGCTCTCA"
        assert record[6].mask == (1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1)
        assert record[6].score == pytest.approx(13.3145, abs=5e-8)
        assert str(record[6]) == """\
GCGCGTGTGTGTAAC
GCACAGAGCTTAGCA
GGTGGTCATCGGGCA
GCGCGTGTCATTGAC
GGACGGCACTTAGCA
GCGCGTCCCGGGCCA
GCTCGGCCCGTTGTC
GCGCGTGTCCTTTAA
GCTGATCGCTGCTCC
GCCCGTACCGGACCT
GGACGTCGCGGAGGA
GCGGGGGGCAGGTCA
GGACGTACTGGCACA
GCAGGTGGTCGTGCA
GCGCATACCTTAACA
GCACGGGACTTCAAC
GCACGTAGCTGGTAA
GCTCGTCTATGGTCA
GCGCATGCTGGATCC
GGCCGTCAGCTCTCA"""
        motif = record[6][2:-1]
        assert motif.alphabet == "ACGT"
        assert str(motif) == """\
GCGTGTGTGTAA
ACAGAGCTTAGC
TGGTCATCGGGC
GCGTGTCATTGA
ACGGCACTTAGC
GCGTCCCGGGCC
TCGGCCCGTTGT
GCGTGTCCTTTA
TGATCGCTGCTC
CCGTACCGGACC
ACGTCGCGGAGG
GGGGGGCAGGTC
ACGTACTGGCAC
AGGTGGTCGTGC
GCATACCTTAAC
ACGGGACTTCAA
ACGTAGCTGGTA
TCGTCTATGGTC
GCATGCTGGATC
CCGTCAGCTCTC"""
        assert record[7].alphabet == "ACGT"
        assert len(record[7].alignment.sequences) == 20
        assert record[7].alignment.sequences[0] == "GAACCGAGGTCCGGTACGGGC"
        assert record[7].alignment.sequences[1] == "GCCCCCCGCATAGTAGGGGGA"
        assert record[7].alignment.sequences[2] == "GTCCCTGGGTAAGCTTGGGGC"
        assert record[7].alignment.sequences[3] == "ACTCCACGCTTCGACACGTGG"
        assert record[7].alignment.sequences[4] == "ATCCTCTGCGTCGCATGGCGG"
        assert record[7].alignment.sequences[5] == "GTTCAATGCTAAGCTCTGTGC"
        assert record[7].alignment.sequences[6] == "GCTCATAGGGACGTCGCGGAG"
        assert record[7].alignment.sequences[7] == "GTCCCGGGCCAATAGCGGCGC"
        assert record[7].alignment.sequences[8] == "GCACTTAGCAGCGTATCGTTA"
        assert record[7].alignment.sequences[9] == "GGCCCTCGGATCGCTTGGGAA"
        assert record[7].alignment.sequences[10] == "CTGCTGGACAACGGGCCGAGC"
        assert record[7].alignment.sequences[11] == "GGGCACTACATAGAGAGTTGC"
        assert record[7].alignment.sequences[12] == "AGCCTCCAGGTCGCATGGAGA"
        assert record[7].alignment.sequences[13] == "AATCGTAGATCAGAGGCGAGA"
        assert record[7].alignment.sequences[14] == "GAACTCCACTAAGACTTGAGA"
        assert record[7].alignment.sequences[15] == "GAGCAGCGATCAGCTTGTGGG"
        assert record[7].alignment.sequences[16] == "GCCAGGTACAAAGCGTCGTGC"
        assert record[7].alignment.sequences[17] == "AGTCAATGACACGCGCCTGGG"
        assert record[7].alignment.sequences[18] == "GGTCATGGAATCTTATGTAGC"
        assert record[7].alignment.sequences[19] == "GTAGATAACAGAGGTCGGGGG"
        assert record[7].mask == (1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1)
        assert record[7].score == pytest.approx(11.6098, abs=5e-8)
        assert str(record[7]) == """\
GAACCGAGGTCCGGTACGGGC
GCCCCCCGCATAGTAGGGGGA
GTCCCTGGGTAAGCTTGGGGC
ACTCCACGCTTCGACACGTGG
ATCCTCTGCGTCGCATGGCGG
GTTCAATGCTAAGCTCTGTGC
GCTCATAGGGACGTCGCGGAG
GTCCCGGGCCAATAGCGGCGC
GCACTTAGCAGCGTATCGTTA
GGCCCTCGGATCGCTTGGGAA
CTGCTGGACAACGGGCCGAGC
GGGCACTACATAGAGAGTTGC
AGCCTCCAGGTCGCATGGAGA
AATCGTAGATCAGAGGCGAGA
GAACTCCACTAAGACTTGAGA
GAGCAGCGATCAGCTTGTGGG
GCCAGGTACAAAGCGTCGTGC
AGTCAATGACACGCGCCTGGG
GGTCATGGAATCTTATGTAGC
GTAGATAACAGAGGTCGGGGG"""
        motif = record[7][2:-1]
        assert motif.alphabet == "ACGT"
        assert str(motif) == """\
ACCGAGGTCCGGTACGGG
CCCCCGCATAGTAGGGGG
CCCTGGGTAAGCTTGGGG
TCCACGCTTCGACACGTG
CCTCTGCGTCGCATGGCG
TCAATGCTAAGCTCTGTG
TCATAGGGACGTCGCGGA
CCCGGGCCAATAGCGGCG
ACTTAGCAGCGTATCGTT
CCCTCGGATCGCTTGGGA
GCTGGACAACGGGCCGAG
GCACTACATAGAGAGTTG
CCTCCAGGTCGCATGGAG
TCGTAGATCAGAGGCGAG
ACTCCACTAAGACTTGAG
GCAGCGATCAGCTTGTGG
CAGGTACAAAGCGTCGTG
TCAATGACACGCGCCTGG
TCATGGAATCTTATGTAG
AGATAACAGAGGTCGGGG"""
        assert record[8].alphabet == "ACGT"
        assert len(record[8].alignment.sequences) == 14
        assert record[8].alignment.sequences[0] == "CCGAGTAAAGGGCTG"
        assert record[8].alignment.sequences[1] == "GTGGTCATCGGGCAC"
        assert record[8].alignment.sequences[2] == "GATAACAGAGGTCGG"
        assert record[8].alignment.sequences[3] == "CGGCGCCGGAGTCTG"
        assert record[8].alignment.sequences[4] == "GCGCGTCCCGGGCCA"
        assert record[8].alignment.sequences[5] == "CTGGACAACGGGCCG"
        assert record[8].alignment.sequences[6] == "CGGATACTGGGGCAG"
        assert record[8].alignment.sequences[7] == "GGGAGCAGCGATCAG"
        assert record[8].alignment.sequences[8] == "CAGAACCGAGGTCCG"
        assert record[8].alignment.sequences[9] == "GGGTCCCTGGGTAAG"
        assert record[8].alignment.sequences[10] == "GTGCTCATAGGGACG"
        assert record[8].alignment.sequences[11] == "GAGATCCGGAGGAGG"
        assert record[8].alignment.sequences[12] == "GCGATCCGAGGGCCG"
        assert record[8].alignment.sequences[13] == "GAGTTCACATGGCTG"
        assert record[8].mask == (1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1)
        assert record[8].score == pytest.approx(11.2943, abs=5e-8)
        assert str(record[8]) == """\
CCGAGTAAAGGGCTG
GTGGTCATCGGGCAC
GATAACAGAGGTCGG
CGGCGCCGGAGTCTG
GCGCGTCCCGGGCCA
CTGGACAACGGGCCG
CGGATACTGGGGCAG
GGGAGCAGCGATCAG
CAGAACCGAGGTCCG
GGGTCCCTGGGTAAG
GTGCTCATAGGGACG
GAGATCCGGAGGAGG
GCGATCCGAGGGCCG
GAGTTCACATGGCTG"""
        motif = record[8][2:-1]
        assert motif.alphabet == "ACGT"
        assert str(motif) == """\
GAGTAAAGGGCT
GGTCATCGGGCA
TAACAGAGGTCG
GCGCCGGAGTCT
GCGTCCCGGGCC
GGACAACGGGCC
GATACTGGGGCA
GAGCAGCGATCA
GAACCGAGGTCC
GTCCCTGGGTAA
GCTCATAGGGAC
GATCCGGAGGAG
GATCCGAGGGCC
GTTCACATGGCT"""
        assert record[9].alphabet == "ACGT"
        assert len(record[9].alignment.sequences) == 18
        assert record[9].alignment.sequences[0] == "TAGAGGCGGTG"
        assert record[9].alignment.sequences[1] == "GCTAAGCTCTG"
        assert record[9].alignment.sequences[2] == "TGGAAGCAGTG"
        assert record[9].alignment.sequences[3] == "GCGAGGCTGTG"
        assert record[9].alignment.sequences[4] == "ACGACGCTTTG"
        assert record[9].alignment.sequences[5] == "GGGACGCGCAC"
        assert record[9].alignment.sequences[6] == "TCGAAGCGTGG"
        assert record[9].alignment.sequences[7] == "TGTATGCGGGG"
        assert record[9].alignment.sequences[8] == "GGTAAGCTTGG"
        assert record[9].alignment.sequences[9] == "TGTACGCTGGG"
        assert record[9].alignment.sequences[10] == "ACTATGCGGGG"
        assert record[9].alignment.sequences[11] == "GGTATGCGCTG"
        assert record[9].alignment.sequences[12] == "GGTACCCGGAG"
        assert record[9].alignment.sequences[13] == "GCGACGCAGAG"
        assert record[9].alignment.sequences[14] == "TGGCGGCGTGG"
        assert record[9].alignment.sequences[15] == "TCTAGGCGGGC"
        assert record[9].alignment.sequences[16] == "AGTATGCTTAG"
        assert record[9].alignment.sequences[17] == "TGGAGGCTTAG"
        assert record[9].mask == (1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1)
        assert record[9].score == pytest.approx(9.7924, abs=5e-8)
        assert str(record[9]) == """\
TAGAGGCGGTG
GCTAAGCTCTG
TGGAAGCAGTG
GCGAGGCTGTG
ACGACGCTTTG
GGGACGCGCAC
TCGAAGCGTGG
TGTATGCGGGG
GGTAAGCTTGG
TGTACGCTGGG
ACTATGCGGGG
GGTATGCGCTG
GGTACCCGGAG
GCGACGCAGAG
TGGCGGCGTGG
TCTAGGCGGGC
AGTATGCTTAG
TGGAGGCTTAG"""
        motif = record[9][2:-1]
        assert motif.alphabet == "ACGT"
        assert str(motif) == """\
GAGGCGGT
TAAGCTCT
GAAGCAGT
GAGGCTGT
GACGCTTT
GACGCGCA
GAAGCGTG
TATGCGGG
TAAGCTTG
TACGCTGG
TATGCGGG
TATGCGCT
TACCCGGA
GACGCAGA
GCGGCGTG
TAGGCGGG
TATGCTTA
GAGGCTTA"""
        assert record[10].alphabet == "ACGT"
        assert len(record[10].alignment.sequences) == 13
        assert record[10].alignment.sequences[0] == "GCACAGAGCTTAGCATTGAAC"
        assert record[10].alignment.sequences[1] == "GTCCGCGGATTCCCAACATGC"
        assert record[10].alignment.sequences[2] == "ATACACAGCCTCGCAAGCCAG"
        assert record[10].alignment.sequences[3] == "GGCCCGGGACGCGCACTAAGA"
        assert record[10].alignment.sequences[4] == "GCCCGTTGTCCAGCAGACGGC"
        assert record[10].alignment.sequences[5] == "GAGCAGCGATCAGCTTGTGGG"
        assert record[10].alignment.sequences[6] == "GAACCGAGGTCCGGTACGGGC"
        assert record[10].alignment.sequences[7] == "GTCCCTGGGTAAGCTTGGGGC"
        assert record[10].alignment.sequences[8] == "GACCTGCCCCCCGCATAGTAG"
        assert record[10].alignment.sequences[9] == "AACCAGCGCATACCTTAACAG"
        assert record[10].alignment.sequences[10] == "ATCCTCTGCGTCGCATGGCGG"
        assert record[10].alignment.sequences[11] == "GACCATAGACGAGCATCAAAG"
        assert record[10].alignment.sequences[12] == "GGCCCTCGGATCGCTTGGGAA"
        assert record[10].mask == (1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1)
        assert record[10].score == pytest.approx(9.01393, abs=5e-8)
        assert str(record[10]) == """\
GCACAGAGCTTAGCATTGAAC
GTCCGCGGATTCCCAACATGC
ATACACAGCCTCGCAAGCCAG
GGCCCGGGACGCGCACTAAGA
GCCCGTTGTCCAGCAGACGGC
GAGCAGCGATCAGCTTGTGGG
GAACCGAGGTCCGGTACGGGC
GTCCCTGGGTAAGCTTGGGGC
GACCTGCCCCCCGCATAGTAG
AACCAGCGCATACCTTAACAG
ATCCTCTGCGTCGCATGGCGG
GACCATAGACGAGCATCAAAG
GGCCCTCGGATCGCTTGGGAA"""
        motif = record[10][2:-1]
        assert motif.alphabet == "ACGT"
        assert str(motif) == """\
ACAGAGCTTAGCATTGAA
CCGCGGATTCCCAACATG
ACACAGCCTCGCAAGCCA
CCCGGGACGCGCACTAAG
CCGTTGTCCAGCAGACGG
GCAGCGATCAGCTTGTGG
ACCGAGGTCCGGTACGGG
CCCTGGGTAAGCTTGGGG
CCTGCCCCCCGCATAGTA
CCAGCGCATACCTTAACA
CCTCTGCGTCGCATGGCG
CCATAGACGAGCATCAAA
CCCTCGGATCGCTTGGGA"""
        assert record[11].alphabet == "ACGT"
        assert len(record[11].alignment.sequences) == 16
        assert record[11].alignment.sequences[0] == "GCCGTCCGTC"
        assert record[11].alignment.sequences[1] == "GGCGTGCGCG"
        assert record[11].alignment.sequences[2] == "GGCGCGTGTC"
        assert record[11].alignment.sequences[3] == "AGCGCGTGTG"
        assert record[11].alignment.sequences[4] == "GCGGTGCGTG"
        assert record[11].alignment.sequences[5] == "AGCGCGTGTC"
        assert record[11].alignment.sequences[6] == "AGCGTCCGCG"
        assert record[11].alignment.sequences[7] == "ACCGTCTGTG"
        assert record[11].alignment.sequences[8] == "GCCATGCGAC"
        assert record[11].alignment.sequences[9] == "ACCACCCGTC"
        assert record[11].alignment.sequences[10] == "GGCGCCGGAG"
        assert record[11].alignment.sequences[11] == "ACCACGTGTC"
        assert record[11].alignment.sequences[12] == "GGCTTGCGAG"
        assert record[11].alignment.sequences[13] == "GCGATCCGAG"
        assert record[11].alignment.sequences[14] == "AGTGCGCGTC"
        assert record[11].alignment.sequences[15] == "AGTGCCCGAG"
        assert record[11].mask == (1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
        assert record[11].score == pytest.approx(7.51121, abs=5e-8)
        assert str(record[11]) == """\
GCCGTCCGTC
GGCGTGCGCG
GGCGCGTGTC
AGCGCGTGTG
GCGGTGCGTG
AGCGCGTGTC
AGCGTCCGCG
ACCGTCTGTG
GCCATGCGAC
ACCACCCGTC
GGCGCCGGAG
ACCACGTGTC
GGCTTGCGAG
GCGATCCGAG
AGTGCGCGTC
AGTGCCCGAG"""
        motif = record[11][2:-1]
        assert motif.alphabet == "ACGT"
        assert str(motif) == """\
CGTCCGT
CGTGCGC
CGCGTGT
CGCGTGT
GGTGCGT
CGCGTGT
CGTCCGC
CGTCTGT
CATGCGA
CACCCGT
CGCCGGA
CACGTGT
CTTGCGA
GATCCGA
TGCGCGT
TGCCCGA"""
        assert record[12].alphabet == "ACGT"
        assert len(record[12].alignment.sequences) == 16
        assert record[12].alignment.sequences[0] == "GCCGACGGGTGGTCATCGGG"
        assert record[12].alignment.sequences[1] == "GCACGACGCTTTGTACCTGG"
        assert record[12].alignment.sequences[2] == "CCTGGGAGGGTTCAATAACG"
        assert record[12].alignment.sequences[3] == "GCGCGTCCCGGGCCAATAGC"
        assert record[12].alignment.sequences[4] == "GCCGTCTGCTGGACAACGGG"
        assert record[12].alignment.sequences[5] == "GTCCCTTCCGGTACATGAGG"
        assert record[12].alignment.sequences[6] == "GCTGCTCCCCGCATACAGCG"
        assert record[12].alignment.sequences[7] == "GCCCCAAGCTTACCCAGGGA"
        assert record[12].alignment.sequences[8] == "ACCGGCTGACGCTAATACGG"
        assert record[12].alignment.sequences[9] == "GCGGGGGGCAGGTCATTACA"
        assert record[12].alignment.sequences[10] == "GCTGGCAGCGTCTAAGAAGG"
        assert record[12].alignment.sequences[11] == "GCAGGTGGTCGTGCAATACG"
        assert record[12].alignment.sequences[12] == "GCTGGTTGAAGTCCCGTGCG"
        assert record[12].alignment.sequences[13] == "GCACGTAGCTGGTAAATAGG"
        assert record[12].alignment.sequences[14] == "GCGGCGTGGATTTCATACAG"
        assert record[12].alignment.sequences[15] == "CCTGGAGGCTTAGACTTGGG"
        assert record[12].mask == (1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1)
        assert record[12].score == pytest.approx(5.63667, abs=5e-8)
        assert str(record[12]) == """\
GCCGACGGGTGGTCATCGGG
GCACGACGCTTTGTACCTGG
CCTGGGAGGGTTCAATAACG
GCGCGTCCCGGGCCAATAGC
GCCGTCTGCTGGACAACGGG
GTCCCTTCCGGTACATGAGG
GCTGCTCCCCGCATACAGCG
GCCCCAAGCTTACCCAGGGA
ACCGGCTGACGCTAATACGG
GCGGGGGGCAGGTCATTACA
GCTGGCAGCGTCTAAGAAGG
GCAGGTGGTCGTGCAATACG
GCTGGTTGAAGTCCCGTGCG
GCACGTAGCTGGTAAATAGG
GCGGCGTGGATTTCATACAG
CCTGGAGGCTTAGACTTGGG"""
        motif = record[12][2:-1]
        assert motif.alphabet == "ACGT"
        assert str(motif) == """\
CGACGGGTGGTCATCGG
ACGACGCTTTGTACCTG
TGGGAGGGTTCAATAAC
GCGTCCCGGGCCAATAG
CGTCTGCTGGACAACGG
CCCTTCCGGTACATGAG
TGCTCCCCGCATACAGC
CCCAAGCTTACCCAGGG
CGGCTGACGCTAATACG
GGGGGGCAGGTCATTAC
TGGCAGCGTCTAAGAAG
AGGTGGTCGTGCAATAC
TGGTTGAAGTCCCGTGC
ACGTAGCTGGTAAATAG
GGCGTGGATTTCATACA
TGGAGGCTTAGACTTGG"""
        assert record[13].alphabet == "ACGT"
        assert len(record[13].alignment.sequences) == 15
        assert record[13].alignment.sequences[0] == "GCCGACGGGTGGTCATCGGG"
        assert record[13].alignment.sequences[1] == "ATCCGCGGACGCTTAGAGGG"
        assert record[13].alignment.sequences[2] == "ACGCTTTGTACCTGGCTTGC"
        assert record[13].alignment.sequences[3] == "ACGGACGGCACTTAGCAGCG"
        assert record[13].alignment.sequences[4] == "GCCGTCTGCTGGACAACGGG"
        assert record[13].alignment.sequences[5] == "ACACACAGACGGTTGAAAGG"
        assert record[13].alignment.sequences[6] == "GCCGATAGTGCTTAAGTTCG"
        assert record[13].alignment.sequences[7] == "CTTGCCCGTACCGGACCTCG"
        assert record[13].alignment.sequences[8] == "ACCGGCTGACGCTAATACGG"
        assert record[13].alignment.sequences[9] == "GCCCCCCGCATAGTAGGGGG"
        assert record[13].alignment.sequences[10] == "GCTGGCAGCGTCTAAGAAGG"
        assert record[13].alignment.sequences[11] == "GCAGGTGGTCGTGCAATACG"
        assert record[13].alignment.sequences[12] == "ACGCACGGGACTTCAACCAG"
        assert record[13].alignment.sequences[13] == "GCACGTAGCTGGTAAATAGG"
        assert record[13].alignment.sequences[14] == "ATCCTCTGCGTCGCATGGCG"
        assert record[13].mask == (1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1)
        assert record[13].score == pytest.approx(3.89842, abs=5e-8)
        assert str(record[13]) == """\
GCCGACGGGTGGTCATCGGG
ATCCGCGGACGCTTAGAGGG
ACGCTTTGTACCTGGCTTGC
ACGGACGGCACTTAGCAGCG
GCCGTCTGCTGGACAACGGG
ACACACAGACGGTTGAAAGG
GCCGATAGTGCTTAAGTTCG
CTTGCCCGTACCGGACCTCG
ACCGGCTGACGCTAATACGG
GCCCCCCGCATAGTAGGGGG
GCTGGCAGCGTCTAAGAAGG
GCAGGTGGTCGTGCAATACG
ACGCACGGGACTTCAACCAG
GCACGTAGCTGGTAAATAGG
ATCCTCTGCGTCGCATGGCG"""
        motif = record[13][2:-1]
        assert motif.alphabet == "ACGT"
        assert str(motif) == """\
CGACGGGTGGTCATCGG
CCGCGGACGCTTAGAGG
GCTTTGTACCTGGCTTG
GGACGGCACTTAGCAGC
CGTCTGCTGGACAACGG
ACACAGACGGTTGAAAG
CGATAGTGCTTAAGTTC
TGCCCGTACCGGACCTC
CGGCTGACGCTAATACG
CCCCCGCATAGTAGGGG
TGGCAGCGTCTAAGAAG
AGGTGGTCGTGCAATAC
GCACGGGACTTCAACCA
ACGTAGCTGGTAAATAG
CCTCTGCGTCGCATGGC"""
        assert record[14].alphabet == "ACGT"
        assert len(record[14].alignment.sequences) == 14
        assert record[14].alignment.sequences[0] == "GAGGCTGTGTAT"
        assert record[14].alignment.sequences[1] == "GAGGTCGGGGGT"
        assert record[14].alignment.sequences[2] == "GACGGACGGCAC"
        assert record[14].alignment.sequences[3] == "TTGGCCCGGGAC"
        assert record[14].alignment.sequences[4] == "GAGGCTCGGCCC"
        assert record[14].alignment.sequences[5] == "CACGCGCTGTAT"
        assert record[14].alignment.sequences[6] == "TAGGCCAGGTAT"
        assert record[14].alignment.sequences[7] == "GAGGTCCGGTAC"
        assert record[14].alignment.sequences[8] == "TACGCTGGGGAT"
        assert record[14].alignment.sequences[9] == "GTCGCGGAGGAT"
        assert record[14].alignment.sequences[10] == "TACGCACGGGAC"
        assert record[14].alignment.sequences[11] == "TACTCCGGGTAC"
        assert record[14].alignment.sequences[12] == "GACGCAGAGGAT"
        assert record[14].alignment.sequences[13] == "TAGGCGGGCCAT"
        assert record[14].mask == (1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1)
        assert record[14].score == pytest.approx(3.33444, abs=5e-8)
        assert str(record[14]) == """\
GAGGCTGTGTAT
GAGGTCGGGGGT
GACGGACGGCAC
TTGGCCCGGGAC
GAGGCTCGGCCC
CACGCGCTGTAT
TAGGCCAGGTAT
GAGGTCCGGTAC
TACGCTGGGGAT
GTCGCGGAGGAT
TACGCACGGGAC
TACTCCGGGTAC
GACGCAGAGGAT
TAGGCGGGCCAT"""
        motif = record[14][2:-1]
        assert motif.alphabet == "ACGT"
        assert str(motif) == """\
GGCTGTGTA
GGTCGGGGG
CGGACGGCA
GGCCCGGGA
GGCTCGGCC
CGCGCTGTA
GGCCAGGTA
GGTCCGGTA
CGCTGGGGA
CGCGGAGGA
CGCACGGGA
CTCCGGGTA
CGCAGAGGA
GGCGGGCCA"""
        assert record[15].alphabet == "ACGT"
        assert len(record[15].alignment.sequences) == 21
        assert record[15].alignment.sequences[0] == "CGGCTCAATCGTAGAGGC"
        assert record[15].alignment.sequences[1] == "CGACGGGTGGTCATCGGG"
        assert record[15].alignment.sequences[2] == "CGCTTAGAGGGCACAAGC"
        assert record[15].alignment.sequences[3] == "TGACACGCGCCTGGGAGG"
        assert record[15].alignment.sequences[4] == "CGATACGCTGCTAAGTGC"
        assert record[15].alignment.sequences[5] == "CGTCCCGGGCCAATAGCG"
        assert record[15].alignment.sequences[6] == "CCACGCTTCGACACGTGG"
        assert record[15].alignment.sequences[7] == "CGTCTGCTGGACAACGGG"
        assert record[15].alignment.sequences[8] == "ACACAGACGGTTGAAAGG"
        assert record[15].alignment.sequences[9] == "TGCTCCCCGCATACAGCG"
        assert record[15].alignment.sequences[10] == "TGAGGCTTGCCCGTACCG"
        assert record[15].alignment.sequences[11] == "TGCCCCAAGCTTACCCAG"
        assert record[15].alignment.sequences[12] == "CGGCTGACGCTAATACGG"
        assert record[15].alignment.sequences[13] == "CGCGACGTCCCTATGAGC"
        assert record[15].alignment.sequences[14] == "TGCCCCCCGCATAGTAGG"
        assert record[15].alignment.sequences[15] == "CGTTGCCTTCTTAGACGC"
        assert record[15].alignment.sequences[16] == "TGACTCAATCGTAGACCC"
        assert record[15].alignment.sequences[17] == "AGTCCCGTGCGTATGTGG"
        assert record[15].alignment.sequences[18] == "AGGCTCGCACGTAGCTGG"
        assert record[15].alignment.sequences[19] == "CCACGCCGCCATGCGACG"
        assert record[15].alignment.sequences[20] == "AGCCTCCAGGTCGCATGG"
        assert record[15].mask == (1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1)
        assert record[15].score == pytest.approx(1.0395, abs=5e-8)
        assert str(record[15]) == """\
CGGCTCAATCGTAGAGGC
CGACGGGTGGTCATCGGG
CGCTTAGAGGGCACAAGC
TGACACGCGCCTGGGAGG
CGATACGCTGCTAAGTGC
CGTCCCGGGCCAATAGCG
CCACGCTTCGACACGTGG
CGTCTGCTGGACAACGGG
ACACAGACGGTTGAAAGG
TGCTCCCCGCATACAGCG
TGAGGCTTGCCCGTACCG
TGCCCCAAGCTTACCCAG
CGGCTGACGCTAATACGG
CGCGACGTCCCTATGAGC
TGCCCCCCGCATAGTAGG
CGTTGCCTTCTTAGACGC
TGACTCAATCGTAGACCC
AGTCCCGTGCGTATGTGG
AGGCTCGCACGTAGCTGG
CCACGCCGCCATGCGACG
AGCCTCCAGGTCGCATGG"""
        motif = record[15][2:-1]
        assert motif.alphabet == "ACGT"
        assert str(motif) == """\
GCTCAATCGTAGAGG
ACGGGTGGTCATCGG
CTTAGAGGGCACAAG
ACACGCGCCTGGGAG
ATACGCTGCTAAGTG
TCCCGGGCCAATAGC
ACGCTTCGACACGTG
TCTGCTGGACAACGG
ACAGACGGTTGAAAG
CTCCCCGCATACAGC
AGGCTTGCCCGTACC
CCCCAAGCTTACCCA
GCTGACGCTAATACG
CGACGTCCCTATGAG
CCCCCCGCATAGTAG
TTGCCTTCTTAGACG
ACTCAATCGTAGACC
TCCCGTGCGTATGTG
GCTCGCACGTAGCTG
ACGCCGCCATGCGAC
CCTCCAGGTCGCATG"""


class TestClusterBuster(unittest.TestCase):
    """Testing parsing Cluster-Buster output files."""

    def test_clusterbuster_parsing_and_output(self):
        """Test if Bio.motifs can parse and output Cluster-Buster PFM files."""
        with open("motifs/clusterbuster.pfm") as stream:
            record = motifs.parse(stream, "clusterbuster")
            assert len(record) == 3
            motif = record[0]
            assert motif.name == "MA0004.1"
            assert motif.alphabet == "GATC"
            assert motif.consensus == "CACGTG"
            assert motif.degenerate_consensus == "CACGTG"
            assert np.allclose(
                    motif.relative_entropy,
                    np.array(
                        [1.278071905112638, 1.7136030428840439, 2.0, 2.0, 2.0, 2.0]
                    ),
                )
            assert motif[1:-2].consensus == "ACG"
            assert motif.length == 6
            assert motif.weight is None
            assert motif.gap is None
            assert motif.counts["G", 0] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["G", 1] == pytest.approx(1.0, abs=5e-8)
            assert motif.counts["G", 2] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["G", 3] == pytest.approx(20.0, abs=5e-8)
            assert motif.counts["G", 4] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["G", 5] == pytest.approx(20.0, abs=5e-8)
            assert motif.counts["A", 0] == pytest.approx(4.0, abs=5e-8)
            assert motif.counts["A", 1] == pytest.approx(19.0, abs=5e-8)
            assert motif.counts["A", 2] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["A", 3] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["A", 4] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["A", 5] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["T", 0] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["T", 1] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["T", 2] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["T", 3] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["T", 4] == pytest.approx(20.0, abs=5e-8)
            assert motif.counts["T", 5] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["C", 0] == pytest.approx(16.0, abs=5e-8)
            assert motif.counts["C", 1] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["C", 2] == pytest.approx(20.0, abs=5e-8)
            assert motif.counts["C", 3] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["C", 4] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["C", 5] == pytest.approx(0.0, abs=5e-8)
            motif = record[1]
            assert motif.name == "MA0006.1"
            assert motif.alphabet == "GATC"
            assert motif.consensus == "TGCGTG"
            assert motif.degenerate_consensus == "YGCGTG"
            assert np.allclose(
                    motif.relative_entropy,
                    np.array(
                        [
                            0.28206397041108283,
                            1.7501177071668148,
                            1.7501177071668148,
                            1.7501177071668148,
                            2.0,
                            2.0,
                        ]
                    ),
                )
            assert motif[1:-2].consensus == "GCG"
            assert motif.length == 6
            assert motif.weight is None
            assert motif.gap is None
            assert motif.counts["G", 0] == pytest.approx(2.0, abs=5e-8)
            assert motif.counts["G", 1] == pytest.approx(23.0, abs=5e-8)
            assert motif.counts["G", 2] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["G", 3] == pytest.approx(23.0, abs=5e-8)
            assert motif.counts["G", 4] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["G", 5] == pytest.approx(24.0, abs=5e-8)
            assert motif.counts["A", 0] == pytest.approx(3.0, abs=5e-8)
            assert motif.counts["A", 1] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["A", 2] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["A", 3] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["A", 4] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["A", 5] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["T", 0] == pytest.approx(11.0, abs=5e-8)
            assert motif.counts["T", 1] == pytest.approx(1.0, abs=5e-8)
            assert motif.counts["T", 2] == pytest.approx(1.0, abs=5e-8)
            assert motif.counts["T", 3] == pytest.approx(1.0, abs=5e-8)
            assert motif.counts["T", 4] == pytest.approx(24.0, abs=5e-8)
            assert motif.counts["T", 5] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["C", 0] == pytest.approx(8.0, abs=5e-8)
            assert motif.counts["C", 1] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["C", 2] == pytest.approx(23.0, abs=5e-8)
            assert motif.counts["C", 3] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["C", 4] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["C", 5] == pytest.approx(0.0, abs=5e-8)
            motif = record[2]
            assert motif.name == "MA0008.1"
            assert motif.alphabet == "GATC"
            assert motif.consensus == "CAATTATT"
            assert motif.degenerate_consensus == "CAATTATT"
            assert np.allclose(
                    motif.relative_entropy,
                    np.array(
                        [
                            0.2549535827226545,
                            1.2358859454459725,
                            2.0,
                            2.0,
                            1.278071905112638,
                            1.7577078109175852,
                            1.7577078109175852,
                            1.5978208097977271,
                        ]
                    ),
                )
            assert motif[1:-2].consensus == "AATTA"
            assert motif.length == 8
            assert motif.weight == 3.0
            assert motif.gap == 10.0
            assert motif.counts["G", 0] == pytest.approx(4.0, abs=5e-8)
            assert motif.counts["G", 1] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["G", 2] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["G", 3] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["G", 4] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["G", 5] == pytest.approx(1.0, abs=5e-8)
            assert motif.counts["G", 6] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["G", 7] == pytest.approx(2.0, abs=5e-8)
            assert motif.counts["A", 0] == pytest.approx(3.0, abs=5e-8)
            assert motif.counts["A", 1] == pytest.approx(21.0, abs=5e-8)
            assert motif.counts["A", 2] == pytest.approx(25.0, abs=5e-8)
            assert motif.counts["A", 3] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["A", 4] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["A", 5] == pytest.approx(24.0, abs=5e-8)
            assert motif.counts["A", 6] == pytest.approx(1.0, abs=5e-8)
            assert motif.counts["A", 7] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["T", 0] == pytest.approx(5.0, abs=5e-8)
            assert motif.counts["T", 1] == pytest.approx(3.0, abs=5e-8)
            assert motif.counts["T", 2] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["T", 3] == pytest.approx(25.0, abs=5e-8)
            assert motif.counts["T", 4] == pytest.approx(20.0, abs=5e-8)
            assert motif.counts["T", 5] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["T", 6] == pytest.approx(24.0, abs=5e-8)
            assert motif.counts["T", 7] == pytest.approx(23.0, abs=5e-8)
            assert motif.counts["C", 0] == pytest.approx(13.0, abs=5e-8)
            assert motif.counts["C", 1] == pytest.approx(1.0, abs=5e-8)
            assert motif.counts["C", 2] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["C", 3] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["C", 4] == pytest.approx(5.0, abs=5e-8)
            assert motif.counts["C", 5] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["C", 6] == pytest.approx(0.0, abs=5e-8)
            assert motif.counts["C", 7] == pytest.approx(0.0, abs=5e-8)
            stream.seek(0)
            assert motifs.write(record, "clusterbuster").split() == stream.read().split()
            stream.seek(0)
            assert motifs.write(record, "clusterbuster", precision=2).split("\n") == [
                    (
                        line
                        if (line.startswith(">") or line.startswith("#"))
                        else "\t".join([f"{x}.00" for x in line.split()])
                    )
                    for line in stream.read().split("\n")
                ]


class TestXMS(unittest.TestCase):
    """Testing parsing xms output files."""

    def test_xms_parsing(self):
        """Test if Bio.motifs can parse and output xms PFM files."""
        with open("motifs/abdb.xms") as stream:
            record = motifs.parse(stream, "xms")
        assert len(record) == 1
        motif = record[0]
        assert motif.name == "Abd-B"
        assert motif.length == 14
        assert motif.alphabet == "GATC"
        assert motif.counts["G", 0] == pytest.approx(0.333333333, abs=5e-8)
        assert motif.counts["G", 1] == pytest.approx(0.379310345, abs=5e-8)
        assert motif.counts["G", 2] == pytest.approx(0.264705882, abs=5e-8)
        assert motif.counts["G", 3] == pytest.approx(0.194444444, abs=5e-8)
        assert motif.counts["G", 4] == pytest.approx(0.102564103, abs=5e-8)
        assert motif.counts["G", 5] == pytest.approx(0.177777778, abs=5e-8)
        assert motif.counts["G", 6] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["G", 7] == pytest.approx(0.022222222, abs=5e-8)
        assert motif.counts["G", 8] == pytest.approx(0.697674419, abs=5e-8)
        assert motif.counts["G", 9] == pytest.approx(0.571428571, abs=5e-8)
        assert motif.counts["G", 10] == pytest.approx(0.150000000, abs=5e-8)
        assert motif.counts["G", 11] == pytest.approx(0.305555556, abs=5e-8)
        assert motif.counts["G", 12] == pytest.approx(0.258064516, abs=5e-8)
        assert motif.counts["G", 13] == pytest.approx(0.259259259, abs=5e-8)
        assert motif.counts["A", 0] == pytest.approx(0.333333333, abs=5e-8)
        assert motif.counts["A", 1] == pytest.approx(0.103448276, abs=5e-8)
        assert motif.counts["A", 2] == pytest.approx(0.264705882, abs=5e-8)
        assert motif.counts["A", 3] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["A", 4] == pytest.approx(0.102564103, abs=5e-8)
        assert motif.counts["A", 5] == pytest.approx(0.244444444, abs=5e-8)
        assert motif.counts["A", 6] == pytest.approx(0.800000000, abs=5e-8)
        assert motif.counts["A", 7] == pytest.approx(0.133333333, abs=5e-8)
        assert motif.counts["A", 8] == pytest.approx(0.046511628, abs=5e-8)
        assert motif.counts["A", 9] == pytest.approx(0.238095238, abs=5e-8)
        assert motif.counts["A", 10] == pytest.approx(0.025000000, abs=5e-8)
        assert motif.counts["A", 11] == pytest.approx(0.222222222, abs=5e-8)
        assert motif.counts["A", 12] == pytest.approx(0.354838710, abs=5e-8)
        assert motif.counts["A", 13] == pytest.approx(0.185185185, abs=5e-8)
        assert motif.counts["T", 0] == pytest.approx(0.125000000, abs=5e-8)
        assert motif.counts["T", 1] == pytest.approx(0.103448276, abs=5e-8)
        assert motif.counts["T", 2] == pytest.approx(0.205882353, abs=5e-8)
        assert motif.counts["T", 3] == pytest.approx(0.777777778, abs=5e-8)
        assert motif.counts["T", 4] == pytest.approx(0.743589744, abs=5e-8)
        assert motif.counts["T", 5] == pytest.approx(0.533333333, abs=5e-8)
        assert motif.counts["T", 6] == pytest.approx(0.155555556, abs=5e-8)
        assert motif.counts["T", 7] == pytest.approx(0.688888889, abs=5e-8)
        assert motif.counts["T", 8] == pytest.approx(0.209302326, abs=5e-8)
        assert motif.counts["T", 9] == pytest.approx(0.095238095, abs=5e-8)
        assert motif.counts["T", 10] == pytest.approx(0.025000000, abs=5e-8)
        assert motif.counts["T", 11] == pytest.approx(0.194444444, abs=5e-8)
        assert motif.counts["T", 12] == pytest.approx(0.129032258, abs=5e-8)
        assert motif.counts["T", 13] == pytest.approx(0.222222222, abs=5e-8)
        assert motif.counts["C", 0] == pytest.approx(0.208333333, abs=5e-8)
        assert motif.counts["C", 1] == pytest.approx(0.413793103, abs=5e-8)
        assert motif.counts["C", 2] == pytest.approx(0.264705882, abs=5e-8)
        assert motif.counts["C", 3] == pytest.approx(0.027777778, abs=5e-8)
        assert motif.counts["C", 4] == pytest.approx(0.051282051, abs=5e-8)
        assert motif.counts["C", 5] == pytest.approx(0.044444444, abs=5e-8)
        assert motif.counts["C", 6] == pytest.approx(0.044444444, abs=5e-8)
        assert motif.counts["C", 7] == pytest.approx(0.155555556, abs=5e-8)
        assert motif.counts["C", 8] == pytest.approx(0.046511628, abs=5e-8)
        assert motif.counts["C", 9] == pytest.approx(0.095238095, abs=5e-8)
        assert motif.counts["C", 10] == pytest.approx(0.800000000, abs=5e-8)
        assert motif.counts["C", 11] == pytest.approx(0.277777778, abs=5e-8)
        assert motif.counts["C", 12] == pytest.approx(0.258064516, abs=5e-8)
        assert motif.counts["C", 13] == pytest.approx(0.333333333, abs=5e-8)
        assert motif.consensus == "GCGTTTATGGCGAC"
        assert motif.degenerate_consensus == "NSNTTTATGGCNNN"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.09689283163718865,
                        0.26557323997556864,
                        0.007815379142180268,
                        1.1150033950025815,
                        0.78848108520697,
                        0.3768768552773923,
                        1.125231003810913,
                        0.7023990165752877,
                        0.7536432192801433,
                        0.3995487907017483,
                        1.0658162802208113,
                        0.022422587676774776,
                        0.07979555429087543,
                        0.03400971806422712,
                    ]
                ),
            )
        assert motif[3::2].consensus == "TTTGGC"
        assert motif[3::2].degenerate_consensus == "TTTGNN"
        assert np.allclose(
                motif[3::2].relative_entropy,
                np.array(
                    [
                        1.1150033950025815,
                        0.3768768552773923,
                        0.7023990165752877,
                        0.3995487907017483,
                        0.022422587676774776,
                        0.03400971806422712,
                    ]
                ),
            )


class TestJASPAR(unittest.TestCase):
    """Testing parsing JASPAR files."""

    def test_pfm_parsing(self):
        """Test if Bio.motifs can parse JASPAR-style pfm files."""
        with open("motifs/SRF.pfm") as stream:
            m = motifs.read(stream, "pfm")
        assert m.length == 12

    def test_pfm_four_columns_parsing(self):
        """Test if Bio.motifs.pfm can parse motifs in position frequency matrix format (4 columns)."""
        with open("motifs/fourcolumns.pfm") as stream:
            record = motifs.parse(stream, "pfm-four-columns")
        assert len(record) == 8
        motif = record[0]
        assert motif.name == ""
        assert motif.length == 8
        assert motif.alphabet == "GATC"
        assert motif.counts["G", 0] == pytest.approx(0.009615385, abs=5e-8)
        assert motif.counts["G", 1] == pytest.approx(0.009615385, abs=5e-8)
        assert motif.counts["G", 2] == pytest.approx(0.009615385, abs=5e-8)
        assert motif.counts["G", 3] == pytest.approx(0.009615385, abs=5e-8)
        assert motif.counts["G", 4] == pytest.approx(0.009615385, abs=5e-8)
        assert motif.counts["G", 5] == pytest.approx(0.009615385, abs=5e-8)
        assert motif.counts["G", 6] == pytest.approx(0.009615385, abs=5e-8)
        assert motif.counts["G", 7] == pytest.approx(0.009615385, abs=5e-8)
        assert motif.counts["A", 0] == pytest.approx(0.009615385, abs=5e-8)
        assert motif.counts["A", 1] == pytest.approx(0.009615385, abs=5e-8)
        assert motif.counts["A", 2] == pytest.approx(0.971153846, abs=5e-8)
        assert motif.counts["A", 3] == pytest.approx(0.009615385, abs=5e-8)
        assert motif.counts["A", 4] == pytest.approx(0.009615385, abs=5e-8)
        assert motif.counts["A", 5] == pytest.approx(0.971153846, abs=5e-8)
        assert motif.counts["A", 6] == pytest.approx(0.009615385, abs=5e-8)
        assert motif.counts["A", 7] == pytest.approx(0.009615385, abs=5e-8)
        assert motif.counts["T", 0] == pytest.approx(0.971153846, abs=5e-8)
        assert motif.counts["T", 1] == pytest.approx(0.971153846, abs=5e-8)
        assert motif.counts["T", 2] == pytest.approx(0.009615385, abs=5e-8)
        assert motif.counts["T", 3] == pytest.approx(0.971153846, abs=5e-8)
        assert motif.counts["T", 4] == pytest.approx(0.009615385, abs=5e-8)
        assert motif.counts["T", 5] == pytest.approx(0.009615385, abs=5e-8)
        assert motif.counts["T", 6] == pytest.approx(0.009615385, abs=5e-8)
        assert motif.counts["T", 7] == pytest.approx(0.971153846, abs=5e-8)
        assert motif.counts["C", 0] == pytest.approx(0.009615385, abs=5e-8)
        assert motif.counts["C", 1] == pytest.approx(0.009615385, abs=5e-8)
        assert motif.counts["C", 2] == pytest.approx(0.009615385, abs=5e-8)
        assert motif.counts["C", 3] == pytest.approx(0.009615385, abs=5e-8)
        assert motif.counts["C", 4] == pytest.approx(0.971153846, abs=5e-8)
        assert motif.counts["C", 5] == pytest.approx(0.009615385, abs=5e-8)
        assert motif.counts["C", 6] == pytest.approx(0.971153846, abs=5e-8)
        assert motif.counts["C", 7] == pytest.approx(0.009615385, abs=5e-8)
        assert motif.consensus == "TTATCACT"
        assert motif.degenerate_consensus == "TTATCACT"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        1.765707971839016,
                        1.765707971839016,
                        1.7657079718390165,
                        1.765707971839016,
                        1.7657079718390158,
                        1.7657079718390165,
                        1.7657079718390158,
                        1.765707971839016,
                    ]
                ),
            )
        assert motif[1:-2].consensus == "TATCA"
        motif = record[1]
        assert motif.name == "ENSG00000197372"
        assert motif.length == 20
        assert motif.alphabet == "GATC"
        assert motif.counts["G", 0] == pytest.approx(0.117054000, abs=5e-8)
        assert motif.counts["G", 1] == pytest.approx(0.364552000, abs=5e-8)
        assert motif.counts["G", 2] == pytest.approx(0.310520000, abs=5e-8)
        assert motif.counts["G", 3] == pytest.approx(0.131007000, abs=5e-8)
        assert motif.counts["G", 4] == pytest.approx(0.176504000, abs=5e-8)
        assert motif.counts["G", 5] == pytest.approx(0.197793000, abs=5e-8)
        assert motif.counts["G", 6] == pytest.approx(0.926202000, abs=5e-8)
        assert motif.counts["G", 7] == pytest.approx(0.983797000, abs=5e-8)
        assert motif.counts["G", 8] == pytest.approx(0.002387000, abs=5e-8)
        assert motif.counts["G", 9] == pytest.approx(0.002418000, abs=5e-8)
        assert motif.counts["G", 10] == pytest.approx(0.001991000, abs=5e-8)
        assert motif.counts["G", 11] == pytest.approx(0.002868000, abs=5e-8)
        assert motif.counts["G", 12] == pytest.approx(0.350783000, abs=5e-8)
        assert motif.counts["G", 13] == pytest.approx(1.000000000, abs=5e-8)
        assert motif.counts["G", 14] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["G", 15] == pytest.approx(1.000000000, abs=5e-8)
        assert motif.counts["G", 16] == pytest.approx(1.000000000, abs=5e-8)
        assert motif.counts["G", 17] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["G", 18] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["G", 19] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["A", 0] == pytest.approx(0.341303000, abs=5e-8)
        assert motif.counts["A", 1] == pytest.approx(0.283785000, abs=5e-8)
        assert motif.counts["A", 2] == pytest.approx(0.491055000, abs=5e-8)
        assert motif.counts["A", 3] == pytest.approx(0.492621000, abs=5e-8)
        assert motif.counts["A", 4] == pytest.approx(0.250645000, abs=5e-8)
        assert motif.counts["A", 5] == pytest.approx(0.276694000, abs=5e-8)
        assert motif.counts["A", 6] == pytest.approx(0.056317000, abs=5e-8)
        assert motif.counts["A", 7] == pytest.approx(0.004470000, abs=5e-8)
        assert motif.counts["A", 8] == pytest.approx(0.936213000, abs=5e-8)
        assert motif.counts["A", 9] == pytest.approx(0.004352000, abs=5e-8)
        assert motif.counts["A", 10] == pytest.approx(0.013277000, abs=5e-8)
        assert motif.counts["A", 11] == pytest.approx(0.968132000, abs=5e-8)
        assert motif.counts["A", 12] == pytest.approx(0.397623000, abs=5e-8)
        assert motif.counts["A", 13] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["A", 14] == pytest.approx(1.000000000, abs=5e-8)
        assert motif.counts["A", 15] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["A", 16] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["A", 17] == pytest.approx(1.000000000, abs=5e-8)
        assert motif.counts["A", 18] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["A", 19] == pytest.approx(1.000000000, abs=5e-8)
        assert motif.counts["T", 0] == pytest.approx(0.409215000, abs=5e-8)
        assert motif.counts["T", 1] == pytest.approx(0.274597000, abs=5e-8)
        assert motif.counts["T", 2] == pytest.approx(0.120217000, abs=5e-8)
        assert motif.counts["T", 3] == pytest.approx(0.300256000, abs=5e-8)
        assert motif.counts["T", 4] == pytest.approx(0.211387000, abs=5e-8)
        assert motif.counts["T", 5] == pytest.approx(0.027444000, abs=5e-8)
        assert motif.counts["T", 6] == pytest.approx(0.002850000, abs=5e-8)
        assert motif.counts["T", 7] == pytest.approx(0.003964000, abs=5e-8)
        assert motif.counts["T", 8] == pytest.approx(0.002613000, abs=5e-8)
        assert motif.counts["T", 9] == pytest.approx(0.989200000, abs=5e-8)
        assert motif.counts["T", 10] == pytest.approx(0.976567000, abs=5e-8)
        assert motif.counts["T", 11] == pytest.approx(0.026737000, abs=5e-8)
        assert motif.counts["T", 12] == pytest.approx(0.199577000, abs=5e-8)
        assert motif.counts["T", 13] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["T", 14] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["T", 15] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["T", 16] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["T", 17] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["T", 18] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["T", 19] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["C", 0] == pytest.approx(0.132427000, abs=5e-8)
        assert motif.counts["C", 1] == pytest.approx(0.077066000, abs=5e-8)
        assert motif.counts["C", 2] == pytest.approx(0.078208000, abs=5e-8)
        assert motif.counts["C", 3] == pytest.approx(0.076117000, abs=5e-8)
        assert motif.counts["C", 4] == pytest.approx(0.361464000, abs=5e-8)
        assert motif.counts["C", 5] == pytest.approx(0.498070000, abs=5e-8)
        assert motif.counts["C", 6] == pytest.approx(0.014631000, abs=5e-8)
        assert motif.counts["C", 7] == pytest.approx(0.007769000, abs=5e-8)
        assert motif.counts["C", 8] == pytest.approx(0.058787000, abs=5e-8)
        assert motif.counts["C", 9] == pytest.approx(0.004030000, abs=5e-8)
        assert motif.counts["C", 10] == pytest.approx(0.008165000, abs=5e-8)
        assert motif.counts["C", 11] == pytest.approx(0.002263000, abs=5e-8)
        assert motif.counts["C", 12] == pytest.approx(0.052017000, abs=5e-8)
        assert motif.counts["C", 13] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["C", 14] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["C", 15] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["C", 16] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["C", 17] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["C", 18] == pytest.approx(1.000000000, abs=5e-8)
        assert motif.counts["C", 19] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.consensus == "TGAACCGGATTAAGAGGACA"
        assert motif.degenerate_consensus == "WNRWNMGGATTANGAGGACA"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.1946677220077018,
                        0.1566211351816578,
                        0.31728135119311995,
                        0.3086747573287918,
                        0.053393542701508756,
                        0.381471417197324,
                        1.5505596169174871,
                        1.8558501430757017,
                        1.6274200195132635,
                        1.8972899364737197,
                        1.809312450637467,
                        1.7709547585539227,
                        0.2549373240046801,
                        2.0,
                        2.0,
                        2.0,
                        2.0,
                        2.0,
                        2.0,
                        2.0,
                    ]
                ),
            )
        assert motif[1:-2].consensus == "GAACCGGATTAAGAGGA"
        motif = record[2]
        assert motif.counts["G", 0] == pytest.approx(0.083333300, abs=5e-8)
        assert motif.counts["G", 1] == pytest.approx(0.083333300, abs=5e-8)
        assert motif.counts["G", 2] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["G", 3] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["G", 4] == pytest.approx(0.083333300, abs=5e-8)
        assert motif.counts["G", 5] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["G", 6] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["G", 7] == pytest.approx(0.333333000, abs=5e-8)
        assert motif.counts["G", 8] == pytest.approx(0.166667000, abs=5e-8)
        assert motif.counts["G", 9] == pytest.approx(0.166667000, abs=5e-8)
        assert motif.counts["G", 10] == pytest.approx(0.416667000, abs=5e-8)
        assert motif.counts["A", 0] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["A", 1] == pytest.approx(0.750000000, abs=5e-8)
        assert motif.counts["A", 2] == pytest.approx(0.833333000, abs=5e-8)
        assert motif.counts["A", 3] == pytest.approx(1.000000000, abs=5e-8)
        assert motif.counts["A", 4] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["A", 5] == pytest.approx(0.333333000, abs=5e-8)
        assert motif.counts["A", 6] == pytest.approx(0.833333000, abs=5e-8)
        assert motif.counts["A", 7] == pytest.approx(0.500000000, abs=5e-8)
        assert motif.counts["A", 8] == pytest.approx(0.500000000, abs=5e-8)
        assert motif.counts["A", 9] == pytest.approx(0.333333000, abs=5e-8)
        assert motif.counts["A", 10] == pytest.approx(0.166667000, abs=5e-8)
        assert motif.counts["T", 0] == pytest.approx(0.583333000, abs=5e-8)
        assert motif.counts["T", 1] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["T", 2] == pytest.approx(0.166667000, abs=5e-8)
        assert motif.counts["T", 3] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["T", 4] == pytest.approx(0.083333300, abs=5e-8)
        assert motif.counts["T", 5] == pytest.approx(0.666667000, abs=5e-8)
        assert motif.counts["T", 6] == pytest.approx(0.166667000, abs=5e-8)
        assert motif.counts["T", 7] == pytest.approx(0.166667000, abs=5e-8)
        assert motif.counts["T", 8] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["T", 9] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["T", 10] == pytest.approx(0.166667000, abs=5e-8)
        assert motif.counts["C", 0] == pytest.approx(0.083333300, abs=5e-8)
        assert motif.counts["C", 1] == pytest.approx(0.166667000, abs=5e-8)
        assert motif.counts["C", 2] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["C", 3] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["C", 4] == pytest.approx(0.833333000, abs=5e-8)
        assert motif.counts["C", 5] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["C", 6] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["C", 7] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["C", 8] == pytest.approx(0.083333300, abs=5e-8)
        assert motif.counts["C", 9] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["C", 10] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.name == "M1734_0.90"
        assert motif.length == 11
        assert motif.alphabet == "GATC"
        assert motif.consensus == "TAAACTAAAAG"
        assert motif.degenerate_consensus == "TAAACTARNNN"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.4489017067534855,
                        0.9591474871280075,
                        1.3499768043761913,
                        2.0,
                        1.1833109116849791,
                        1.0817044992792044,
                        1.3499768043761913,
                        0.5408517496401433,
                        0.2704258182036411,
                        0.04085174964014324,
                        0.1120812409282564,
                    ]
                ),
            )
        assert motif[1:-2].consensus == "AAACTAAA"
        motif = record[3]
        assert motif.name == "AbdA_Cell_FBgn0000014"
        assert motif.length == 7
        assert motif.alphabet == "GATC"
        assert motif.counts["G", 0] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["G", 1] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["G", 2] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["G", 3] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["G", 4] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["G", 5] == pytest.approx(6.000000000, abs=5e-8)
        assert motif.counts["G", 6] == pytest.approx(2.000000000, abs=5e-8)
        assert motif.counts["A", 0] == pytest.approx(1.000000000, abs=5e-8)
        assert motif.counts["A", 1] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["A", 2] == pytest.approx(16.000000000, abs=5e-8)
        assert motif.counts["A", 3] == pytest.approx(18.000000000, abs=5e-8)
        assert motif.counts["A", 4] == pytest.approx(1.000000000, abs=5e-8)
        assert motif.counts["A", 5] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["A", 6] == pytest.approx(15.000000000, abs=5e-8)
        assert motif.counts["T", 0] == pytest.approx(14.000000000, abs=5e-8)
        assert motif.counts["T", 1] == pytest.approx(18.000000000, abs=5e-8)
        assert motif.counts["T", 2] == pytest.approx(2.000000000, abs=5e-8)
        assert motif.counts["T", 3] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["T", 4] == pytest.approx(17.000000000, abs=5e-8)
        assert motif.counts["T", 5] == pytest.approx(12.000000000, abs=5e-8)
        assert motif.counts["T", 6] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["C", 0] == pytest.approx(3.000000000, abs=5e-8)
        assert motif.counts["C", 1] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["C", 2] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["C", 3] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["C", 4] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["C", 5] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["C", 6] == pytest.approx(1.000000000, abs=5e-8)
        assert motif.consensus == "TTAATTA"
        assert motif.degenerate_consensus == "TTAATKA"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        1.0555114658337947,
                        2.0,
                        1.4967416652243541,
                        2.0,
                        1.6904565708496748,
                        1.0817041659455104,
                        1.1969282726758976,
                    ]
                ),
            )
        assert motif[1:-2].consensus == "TAAT"
        motif = record[4]
        assert motif.name == "ATGACTCATC AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer    6.049537    -1.782996e+03   0   9805.3,5781.0,3085.1,2715.0,0.00e+00"
        assert motif.length == 10
        assert motif.alphabet == "GATC"
        assert motif.counts["G", 0] == pytest.approx(0.277000000, abs=5e-8)
        assert motif.counts["G", 1] == pytest.approx(0.001000000, abs=5e-8)
        assert motif.counts["G", 2] == pytest.approx(0.965000000, abs=5e-8)
        assert motif.counts["G", 3] == pytest.approx(0.001000000, abs=5e-8)
        assert motif.counts["G", 4] == pytest.approx(0.305000000, abs=5e-8)
        assert motif.counts["G", 5] == pytest.approx(0.001000000, abs=5e-8)
        assert motif.counts["G", 6] == pytest.approx(0.001000000, abs=5e-8)
        assert motif.counts["G", 7] == pytest.approx(0.001000000, abs=5e-8)
        assert motif.counts["G", 8] == pytest.approx(0.307000000, abs=5e-8)
        assert motif.counts["G", 9] == pytest.approx(0.211000000, abs=5e-8)
        assert motif.counts["A", 0] == pytest.approx(0.419000000, abs=5e-8)
        assert motif.counts["A", 1] == pytest.approx(0.001000000, abs=5e-8)
        assert motif.counts["A", 2] == pytest.approx(0.010000000, abs=5e-8)
        assert motif.counts["A", 3] == pytest.approx(0.984000000, abs=5e-8)
        assert motif.counts["A", 4] == pytest.approx(0.062000000, abs=5e-8)
        assert motif.counts["A", 5] == pytest.approx(0.026000000, abs=5e-8)
        assert motif.counts["A", 6] == pytest.approx(0.043000000, abs=5e-8)
        assert motif.counts["A", 7] == pytest.approx(0.980000000, abs=5e-8)
        assert motif.counts["A", 8] == pytest.approx(0.050000000, abs=5e-8)
        assert motif.counts["A", 9] == pytest.approx(0.149000000, abs=5e-8)
        assert motif.counts["T", 0] == pytest.approx(0.028000000, abs=5e-8)
        assert motif.counts["T", 1] == pytest.approx(0.997000000, abs=5e-8)
        assert motif.counts["T", 2] == pytest.approx(0.023000000, abs=5e-8)
        assert motif.counts["T", 3] == pytest.approx(0.012000000, abs=5e-8)
        assert motif.counts["T", 4] == pytest.approx(0.054000000, abs=5e-8)
        assert motif.counts["T", 5] == pytest.approx(0.972000000, abs=5e-8)
        assert motif.counts["T", 6] == pytest.approx(0.012000000, abs=5e-8)
        assert motif.counts["T", 7] == pytest.approx(0.014000000, abs=5e-8)
        assert motif.counts["T", 8] == pytest.approx(0.471000000, abs=5e-8)
        assert motif.counts["T", 9] == pytest.approx(0.195000000, abs=5e-8)
        assert motif.counts["C", 0] == pytest.approx(0.275000000, abs=5e-8)
        assert motif.counts["C", 1] == pytest.approx(0.001000000, abs=5e-8)
        assert motif.counts["C", 2] == pytest.approx(0.002000000, abs=5e-8)
        assert motif.counts["C", 3] == pytest.approx(0.003000000, abs=5e-8)
        assert motif.counts["C", 4] == pytest.approx(0.579000000, abs=5e-8)
        assert motif.counts["C", 5] == pytest.approx(0.001000000, abs=5e-8)
        assert motif.counts["C", 6] == pytest.approx(0.943000000, abs=5e-8)
        assert motif.counts["C", 7] == pytest.approx(0.005000000, abs=5e-8)
        assert motif.counts["C", 8] == pytest.approx(0.172000000, abs=5e-8)
        assert motif.counts["C", 9] == pytest.approx(0.444000000, abs=5e-8)
        assert motif.consensus == "ATGACTCATC"
        assert motif.degenerate_consensus == "NTGASTCAKN"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.30427230622817475,
                        1.9657810606529142,
                        1.7408585738061,
                        1.8654244261025423,
                        0.5449286810918202,
                        1.8033449015144003,
                        1.639502374827662,
                        1.8370335049436752,
                        0.3124728907316759,
                        0.13671828556764112,
                    ]
                ),
            )
        assert motif[1:-2].consensus == "TGACTCA"
        motif = record[5]
        assert motif.name == "AHR_si"
        assert motif.length == 9
        assert motif.alphabet == "GATC"
        assert motif.counts["G", 0] == pytest.approx(56.412537571, abs=5e-8)
        assert motif.counts["G", 1] == pytest.approx(34.663129823, abs=5e-8)
        assert motif.counts["G", 2] == pytest.approx(20.706746562, abs=5e-8)
        assert motif.counts["G", 3] == pytest.approx(145.863705132, abs=5e-8)
        assert motif.counts["G", 4] == pytest.approx(1.492783630, abs=5e-8)
        assert motif.counts["G", 5] == pytest.approx(149.376137203, abs=5e-8)
        assert motif.counts["G", 6] == pytest.approx(0.702486414, abs=5e-8)
        assert motif.counts["G", 7] == pytest.approx(153.958717377, abs=5e-8)
        assert motif.counts["G", 8] == pytest.approx(16.159862547, abs=5e-8)
        assert motif.counts["A", 0] == pytest.approx(40.513432405, abs=5e-8)
        assert motif.counts["A", 1] == pytest.approx(10.877470983, abs=5e-8)
        assert motif.counts["A", 2] == pytest.approx(21.716570782, abs=5e-8)
        assert motif.counts["A", 3] == pytest.approx(2.546513251, abs=5e-8)
        assert motif.counts["A", 4] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["A", 5] == pytest.approx(3.441039751, abs=5e-8)
        assert motif.counts["A", 6] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["A", 7] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["A", 8] == pytest.approx(43.079223333, abs=5e-8)
        assert motif.counts["T", 0] == pytest.approx(38.773634853, abs=5e-8)
        assert motif.counts["T", 1] == pytest.approx(96.547239851, abs=5e-8)
        assert motif.counts["T", 2] == pytest.approx(67.652320196, abs=5e-8)
        assert motif.counts["T", 3] == pytest.approx(4.231336967, abs=5e-8)
        assert motif.counts["T", 4] == pytest.approx(2.107459242, abs=5e-8)
        assert motif.counts["T", 5] == pytest.approx(0.351243207, abs=5e-8)
        assert motif.counts["T", 6] == pytest.approx(149.815191211, abs=5e-8)
        assert motif.counts["T", 7] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["T", 8] == pytest.approx(27.844049228, abs=5e-8)
        assert motif.counts["C", 0] == pytest.approx(18.259112548, abs=5e-8)
        assert motif.counts["C", 1] == pytest.approx(11.870876720, abs=5e-8)
        assert motif.counts["C", 2] == pytest.approx(43.883079838, abs=5e-8)
        assert motif.counts["C", 3] == pytest.approx(1.317162026, abs=5e-8)
        assert motif.counts["C", 4] == pytest.approx(150.358474505, abs=5e-8)
        assert motif.counts["C", 5] == pytest.approx(0.790297216, abs=5e-8)
        assert motif.counts["C", 6] == pytest.approx(3.441039751, abs=5e-8)
        assert motif.counts["C", 7] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["C", 8] == pytest.approx(66.875582269, abs=5e-8)
        assert motif.consensus == "GTTGCGTGC"
        assert motif.degenerate_consensus == "NTNGCGTGN"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.09662409645348236,
                        0.5383413903068038,
                        0.17471270188228985,
                        1.6270151623731723,
                        1.8170663607301638,
                        1.7760800937680195,
                        1.803660112630464,
                        2.0,
                        0.17577786614573548,
                    ]
                ),
            )
        assert motif[1:-2].consensus == "TTGCGT"
        motif = record[6]
        assert motif.name == ""
        assert motif.length == 8
        assert motif.alphabet == "GATC"
        assert motif.counts["G", 0] == pytest.approx(0.098612000, abs=5e-8)
        assert motif.counts["G", 1] == pytest.approx(0.025056000, abs=5e-8)
        assert motif.counts["G", 2] == pytest.approx(0.918728000, abs=5e-8)
        assert motif.counts["G", 3] == pytest.approx(0.029759000, abs=5e-8)
        assert motif.counts["G", 4] == pytest.approx(0.104968000, abs=5e-8)
        assert motif.counts["G", 5] == pytest.approx(0.006667000, abs=5e-8)
        assert motif.counts["G", 6] == pytest.approx(0.026928000, abs=5e-8)
        assert motif.counts["G", 7] == pytest.approx(0.005737000, abs=5e-8)
        assert motif.counts["A", 0] == pytest.approx(0.772949000, abs=5e-8)
        assert motif.counts["A", 1] == pytest.approx(0.026652000, abs=5e-8)
        assert motif.counts["A", 2] == pytest.approx(0.017663000, abs=5e-8)
        assert motif.counts["A", 3] == pytest.approx(0.919596000, abs=5e-8)
        assert motif.counts["A", 4] == pytest.approx(0.060312000, abs=5e-8)
        assert motif.counts["A", 5] == pytest.approx(0.037406000, abs=5e-8)
        assert motif.counts["A", 6] == pytest.approx(0.047316000, abs=5e-8)
        assert motif.counts["A", 7] == pytest.approx(0.948639000, abs=5e-8)
        assert motif.counts["T", 0] == pytest.approx(0.038860000, abs=5e-8)
        assert motif.counts["T", 1] == pytest.approx(0.943639000, abs=5e-8)
        assert motif.counts["T", 2] == pytest.approx(0.040264000, abs=5e-8)
        assert motif.counts["T", 3] == pytest.approx(0.025231000, abs=5e-8)
        assert motif.counts["T", 4] == pytest.approx(0.062462000, abs=5e-8)
        assert motif.counts["T", 5] == pytest.approx(0.935284000, abs=5e-8)
        assert motif.counts["T", 6] == pytest.approx(0.026732000, abs=5e-8)
        assert motif.counts["T", 7] == pytest.approx(0.026128000, abs=5e-8)
        assert motif.counts["C", 0] == pytest.approx(0.089579000, abs=5e-8)
        assert motif.counts["C", 1] == pytest.approx(0.004653000, abs=5e-8)
        assert motif.counts["C", 2] == pytest.approx(0.023344000, abs=5e-8)
        assert motif.counts["C", 3] == pytest.approx(0.025414000, abs=5e-8)
        assert motif.counts["C", 4] == pytest.approx(0.772259000, abs=5e-8)
        assert motif.counts["C", 5] == pytest.approx(0.020643000, abs=5e-8)
        assert motif.counts["C", 6] == pytest.approx(0.899024000, abs=5e-8)
        assert motif.counts["C", 7] == pytest.approx(0.019497000, abs=5e-8)
        assert motif.consensus == "ATGACTCA"
        assert motif.degenerate_consensus == "ATGACTCA"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.889358068874075,
                        1.6123293058245811,
                        1.471654165929799,
                        1.4693092198124151,
                        0.8764628815119266,
                        1.5686388858173408,
                        1.37357038822754,
                        1.6369796776980579,
                    ]
                ),
            )
        assert motif[1:-2].consensus == "TGACT"
        motif = record[7]
        assert motif.name == ""
        assert motif.length == 11
        assert motif.alphabet == "GATC"
        assert motif.counts["G", 0] == pytest.approx(28.0, abs=5e-8)
        assert motif.counts["G", 1] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["G", 2] == pytest.approx(14.0, abs=5e-8)
        assert motif.counts["G", 3] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["G", 4] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["G", 5] == pytest.approx(7.0, abs=5e-8)
        assert motif.counts["G", 6] == pytest.approx(11.0, abs=5e-8)
        assert motif.counts["G", 7] == pytest.approx(38.0, abs=5e-8)
        assert motif.counts["G", 8] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["G", 9] == pytest.approx(25.0, abs=5e-8)
        assert motif.counts["G", 10] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["A", 0] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["A", 1] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["A", 2] == pytest.approx(55.0, abs=5e-8)
        assert motif.counts["A", 3] == pytest.approx(99.0, abs=5e-8)
        assert motif.counts["A", 4] == pytest.approx(78.0, abs=5e-8)
        assert motif.counts["A", 5] == pytest.approx(52.0, abs=5e-8)
        assert motif.counts["A", 6] == pytest.approx(46.0, abs=5e-8)
        assert motif.counts["A", 7] == pytest.approx(60.0, abs=5e-8)
        assert motif.counts["A", 8] == pytest.approx(33.0, abs=5e-8)
        assert motif.counts["A", 9] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["A", 10] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["T", 0] == pytest.approx(30.0, abs=5e-8)
        assert motif.counts["T", 1] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["T", 2] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["T", 3] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["T", 4] == pytest.approx(20.0, abs=5e-8)
        assert motif.counts["T", 5] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["T", 6] == pytest.approx(19.0, abs=5e-8)
        assert motif.counts["T", 7] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["T", 8] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["T", 9] == pytest.approx(73.0, abs=5e-8)
        assert motif.counts["T", 10] == pytest.approx(99.0, abs=5e-8)
        assert motif.counts["C", 0] == pytest.approx(40.0, abs=5e-8)
        assert motif.counts["C", 1] == pytest.approx(99.0, abs=5e-8)
        assert motif.counts["C", 2] == pytest.approx(29.0, abs=5e-8)
        assert motif.counts["C", 3] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["C", 4] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["C", 5] == pytest.approx(39.0, abs=5e-8)
        assert motif.counts["C", 6] == pytest.approx(22.0, abs=5e-8)
        assert motif.counts["C", 7] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["C", 8] == pytest.approx(66.0, abs=5e-8)
        assert motif.counts["C", 9] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["C", 10] == pytest.approx(0.0, abs=5e-8)
        assert motif.consensus == "CCAAAAAACTT"
        assert motif.degenerate_consensus == "BCMAAMNRMTT"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.43314504855176084,
                        2.0,
                        0.6114044621231828,
                        2.0,
                        1.2699833698542062,
                        0.7139129756130338,
                        0.1909607288346033,
                        1.0366644543273158,
                        1.0817041659455104,
                        1.180735028768561,
                        2.0,
                    ]
                ),
            )
        assert motif[1:-2].consensus == "CAAAAAAC"

    def test_pfm_four_rows_parsing(self):
        """Test if Bio.motifs.pfm can parse motifs in position frequency matrix format (4 rows)."""
        with open("motifs/fourrows.pfm") as stream:
            record = motifs.parse(stream, "pfm-four-rows")
        assert len(record) == 9
        motif = record[0]
        assert motif.name == ""
        assert motif.length == 6
        assert motif.alphabet == "GATC"
        assert motif.counts["G", 0] == pytest.approx(5.0, abs=5e-8)
        assert motif.counts["G", 1] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["G", 2] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["G", 3] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["G", 4] == pytest.approx(3.0, abs=5e-8)
        assert motif.counts["G", 5] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["A", 0] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["A", 1] == pytest.approx(5.0, abs=5e-8)
        assert motif.counts["A", 2] == pytest.approx(6.0, abs=5e-8)
        assert motif.counts["A", 3] == pytest.approx(5.0, abs=5e-8)
        assert motif.counts["A", 4] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["A", 5] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["T", 0] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["T", 1] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["T", 2] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["T", 3] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["T", 4] == pytest.approx(2.0, abs=5e-8)
        assert motif.counts["T", 5] == pytest.approx(2.0, abs=5e-8)
        assert motif.counts["C", 0] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["C", 1] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["C", 2] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["C", 3] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["C", 4] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["C", 5] == pytest.approx(4.0, abs=5e-8)
        assert motif.consensus == "GAAAGC"
        assert motif.degenerate_consensus == "GAAAKY"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        1.349977578351646,
                        1.349977578351646,
                        2.0,
                        1.349977578351646,
                        0.5408520829727552,
                        1.0817041659455104,
                    ]
                ),
            )
        assert motif[:-2].consensus == "GAAA"
        motif = record[1]
        assert motif.name == ""
        assert motif.length == 15
        assert motif.alphabet == "GATC"
        assert motif.counts["G", 0] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["G", 1] == pytest.approx(1.000000000, abs=5e-8)
        assert motif.counts["G", 2] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["G", 3] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["G", 4] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["G", 5] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["G", 6] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["G", 7] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["G", 8] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["G", 9] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["G", 10] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["G", 11] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["G", 12] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["G", 13] == pytest.approx(1.000000000, abs=5e-8)
        assert motif.counts["G", 14] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["A", 0] == pytest.approx(0.500000000, abs=5e-8)
        assert motif.counts["A", 1] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["A", 2] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["A", 3] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["A", 4] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["A", 5] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["A", 6] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["A", 7] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["A", 8] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["A", 9] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["A", 10] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["A", 11] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["A", 12] == pytest.approx(0.500000000, abs=5e-8)
        assert motif.counts["A", 13] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["A", 14] == pytest.approx(0.083333333, abs=5e-8)
        assert motif.counts["T", 0] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["T", 1] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["T", 2] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["T", 3] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["T", 4] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["T", 5] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["T", 6] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["T", 7] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["T", 8] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["T", 9] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["T", 10] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["T", 11] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["T", 12] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["T", 13] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["T", 14] == pytest.approx(0.083333333, abs=5e-8)
        assert motif.counts["C", 0] == pytest.approx(0.500000000, abs=5e-8)
        assert motif.counts["C", 1] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["C", 2] == pytest.approx(1.000000000, abs=5e-8)
        assert motif.counts["C", 3] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["C", 4] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["C", 5] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["C", 6] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["C", 7] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["C", 8] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["C", 9] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["C", 10] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["C", 11] == pytest.approx(0.250000000, abs=5e-8)
        assert motif.counts["C", 12] == pytest.approx(0.500000000, abs=5e-8)
        assert motif.counts["C", 13] == pytest.approx(0.000000000, abs=5e-8)
        assert motif.counts["C", 14] == pytest.approx(0.583333333, abs=5e-8)
        assert motif.consensus == "AGCGGGGGGGGGAGC"
        assert motif.degenerate_consensus == "MGCNNNNNNNNNMGC"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        1.0,
                        2.0,
                        2.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        1.0,
                        2.0,
                        0.44890182844369547,
                    ]
                ),
            )
        assert motif[:-2].consensus == "AGCGGGGGGGGGA"
        motif = record[2]
        assert motif.name == ""
        assert motif.length == 15
        assert motif.alphabet == "GATC"
        assert motif.counts["G", 0] == pytest.approx(270.0, abs=5e-8)
        assert motif.counts["G", 1] == pytest.approx(398.0, abs=5e-8)
        assert motif.counts["G", 2] == pytest.approx(54.0, abs=5e-8)
        assert motif.counts["G", 3] == pytest.approx(164.0, abs=5e-8)
        assert motif.counts["G", 4] == pytest.approx(7.0, abs=5e-8)
        assert motif.counts["G", 5] == pytest.approx(659.0, abs=5e-8)
        assert motif.counts["G", 6] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["G", 7] == pytest.approx(750.0, abs=5e-8)
        assert motif.counts["G", 8] == pytest.approx(755.0, abs=5e-8)
        assert motif.counts["G", 9] == pytest.approx(65.0, abs=5e-8)
        assert motif.counts["G", 10] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["G", 11] == pytest.approx(41.0, abs=5e-8)
        assert motif.counts["G", 12] == pytest.approx(202.0, abs=5e-8)
        assert motif.counts["G", 13] == pytest.approx(234.0, abs=5e-8)
        assert motif.counts["G", 14] == pytest.approx(205.0, abs=5e-8)
        assert motif.counts["A", 0] == pytest.approx(92.0, abs=5e-8)
        assert motif.counts["A", 1] == pytest.approx(106.0, abs=5e-8)
        assert motif.counts["A", 2] == pytest.approx(231.0, abs=5e-8)
        assert motif.counts["A", 3] == pytest.approx(135.0, abs=5e-8)
        assert motif.counts["A", 4] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["A", 5] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["A", 6] == pytest.approx(780.0, abs=5e-8)
        assert motif.counts["A", 7] == pytest.approx(28.0, abs=5e-8)
        assert motif.counts["A", 8] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["A", 9] == pytest.approx(700.0, abs=5e-8)
        assert motif.counts["A", 10] == pytest.approx(739.0, abs=5e-8)
        assert motif.counts["A", 11] == pytest.approx(94.0, abs=5e-8)
        assert motif.counts["A", 12] == pytest.approx(60.0, abs=5e-8)
        assert motif.counts["A", 13] == pytest.approx(127.0, abs=5e-8)
        assert motif.counts["A", 14] == pytest.approx(130.0, abs=5e-8)
        assert motif.counts["T", 0] == pytest.approx(290.0, abs=5e-8)
        assert motif.counts["T", 1] == pytest.approx(204.0, abs=5e-8)
        assert motif.counts["T", 2] == pytest.approx(375.0, abs=5e-8)
        assert motif.counts["T", 3] == pytest.approx(411.0, abs=5e-8)
        assert motif.counts["T", 4] == pytest.approx(9.0, abs=5e-8)
        assert motif.counts["T", 5] == pytest.approx(127.0, abs=5e-8)
        assert motif.counts["T", 6] == pytest.approx(6.0, abs=5e-8)
        assert motif.counts["T", 7] == pytest.approx(11.0, abs=5e-8)
        assert motif.counts["T", 8] == pytest.approx(36.0, abs=5e-8)
        assert motif.counts["T", 9] == pytest.approx(20.0, abs=5e-8)
        assert motif.counts["T", 10] == pytest.approx(31.0, abs=5e-8)
        assert motif.counts["T", 11] == pytest.approx(605.0, abs=5e-8)
        assert motif.counts["T", 12] == pytest.approx(335.0, abs=5e-8)
        assert motif.counts["T", 13] == pytest.approx(307.0, abs=5e-8)
        assert motif.counts["T", 14] == pytest.approx(308.0, abs=5e-8)
        assert motif.counts["C", 0] == pytest.approx(138.0, abs=5e-8)
        assert motif.counts["C", 1] == pytest.approx(82.0, abs=5e-8)
        assert motif.counts["C", 2] == pytest.approx(129.0, abs=5e-8)
        assert motif.counts["C", 3] == pytest.approx(81.0, abs=5e-8)
        assert motif.counts["C", 4] == pytest.approx(774.0, abs=5e-8)
        assert motif.counts["C", 5] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["C", 6] == pytest.approx(3.0, abs=5e-8)
        assert motif.counts["C", 7] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["C", 8] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["C", 9] == pytest.approx(6.0, abs=5e-8)
        assert motif.counts["C", 10] == pytest.approx(17.0, abs=5e-8)
        assert motif.counts["C", 11] == pytest.approx(49.0, abs=5e-8)
        assert motif.counts["C", 12] == pytest.approx(193.0, abs=5e-8)
        assert motif.counts["C", 13] == pytest.approx(122.0, abs=5e-8)
        assert motif.counts["C", 14] == pytest.approx(148.0, abs=5e-8)
        assert motif.consensus == "TGTTCGAGGAATTTT"
        assert motif.degenerate_consensus == "NKWTCGAGGAATNNN"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.13892143832881046,
                        0.2692660952911542,
                        0.27915566353819243,
                        0.2665840150038887,
                        1.8371160692433293,
                        1.3354706334248059,
                        1.8856611660889357,
                        1.6600123906824402,
                        1.7329826640509962,
                        1.3601399752384014,
                        1.5978925123167893,
                        0.8698961051280728,
                        0.19290147849975406,
                        0.11003972948477392,
                        0.08469189143040626,
                    ]
                ),
            )
        assert motif[:-2].consensus == "TGTTCGAGGAATT"
        motif = record[3]
        assert motif.name == ""
        assert motif.length == 6
        assert motif.alphabet == "GATC"
        assert motif.counts["G", 0] == pytest.approx(2.0, abs=5e-8)
        assert motif.counts["G", 1] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["G", 2] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["G", 3] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["G", 4] == pytest.approx(97.0, abs=5e-8)
        assert motif.counts["G", 5] == pytest.approx(2.0, abs=5e-8)
        assert motif.counts["A", 0] == pytest.approx(9.0, abs=5e-8)
        assert motif.counts["A", 1] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["A", 2] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["A", 3] == pytest.approx(97.0, abs=5e-8)
        assert motif.counts["A", 4] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["A", 5] == pytest.approx(94.0, abs=5e-8)
        assert motif.counts["T", 0] == pytest.approx(80.0, abs=5e-8)
        assert motif.counts["T", 1] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["T", 2] == pytest.approx(97.0, abs=5e-8)
        assert motif.counts["T", 3] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["T", 4] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["T", 5] == pytest.approx(2.0, abs=5e-8)
        assert motif.counts["C", 0] == pytest.approx(9.0, abs=5e-8)
        assert motif.counts["C", 1] == pytest.approx(97.0, abs=5e-8)
        assert motif.counts["C", 2] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["C", 3] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["C", 4] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["C", 5] == pytest.approx(2.0, abs=5e-8)
        assert motif.consensus == "TCTAGA"
        assert motif.degenerate_consensus == "TCTAGA"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        1.0042727863947818,
                        1.758059267146789,
                        1.7580592671467892,
                        1.7580592671467892,
                        1.7580592671467892,
                        1.5774573308022544,
                    ]
                ),
            )
        assert motif[:-2].consensus == "TCTA"
        motif = record[4]
        assert motif.name == ""
        assert motif.length == 6
        assert motif.alphabet == "GATC"
        assert motif.counts["G", 0] == pytest.approx(0.02, abs=5e-8)
        assert motif.counts["G", 1] == pytest.approx(0.01, abs=5e-8)
        assert motif.counts["G", 2] == pytest.approx(0.01, abs=5e-8)
        assert motif.counts["G", 3] == pytest.approx(0.01, abs=5e-8)
        assert motif.counts["G", 4] == pytest.approx(0.97, abs=5e-8)
        assert motif.counts["G", 5] == pytest.approx(0.02, abs=5e-8)
        assert motif.counts["A", 0] == pytest.approx(0.09, abs=5e-8)
        assert motif.counts["A", 1] == pytest.approx(0.01, abs=5e-8)
        assert motif.counts["A", 2] == pytest.approx(0.01, abs=5e-8)
        assert motif.counts["A", 3] == pytest.approx(0.97, abs=5e-8)
        assert motif.counts["A", 4] == pytest.approx(0.01, abs=5e-8)
        assert motif.counts["A", 5] == pytest.approx(0.94, abs=5e-8)
        assert motif.counts["T", 0] == pytest.approx(0.80, abs=5e-8)
        assert motif.counts["T", 1] == pytest.approx(0.01, abs=5e-8)
        assert motif.counts["T", 2] == pytest.approx(0.97, abs=5e-8)
        assert motif.counts["T", 3] == pytest.approx(0.01, abs=5e-8)
        assert motif.counts["T", 4] == pytest.approx(0.01, abs=5e-8)
        assert motif.counts["T", 5] == pytest.approx(0.02, abs=5e-8)
        assert motif.counts["C", 0] == pytest.approx(0.09, abs=5e-8)
        assert motif.counts["C", 1] == pytest.approx(0.97, abs=5e-8)
        assert motif.counts["C", 2] == pytest.approx(0.01, abs=5e-8)
        assert motif.counts["C", 3] == pytest.approx(0.01, abs=5e-8)
        assert motif.counts["C", 4] == pytest.approx(0.01, abs=5e-8)
        assert motif.counts["C", 5] == pytest.approx(0.02, abs=5e-8)
        assert motif.consensus == "TCTAGA"
        assert motif.degenerate_consensus == "TCTAGA"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        1.0042727863947818,
                        1.758059267146789,
                        1.7580592671467892,
                        1.7580592671467892,
                        1.7580592671467892,
                        1.5774573308022544,
                    ]
                ),
            )
        assert motif[:-2].consensus == "TCTA"
        motif = record[5]
        assert motif.name == "abd-A"
        assert motif.length == 8
        assert motif.alphabet == "GATC"
        assert motif.counts["G", 0] == pytest.approx(0.455991516, abs=5e-8)
        assert motif.counts["G", 1] == pytest.approx(0.069194062, abs=5e-8)
        assert motif.counts["G", 2] == pytest.approx(0.010869565, abs=5e-8)
        assert motif.counts["G", 3] == pytest.approx(0.021739130, abs=5e-8)
        assert motif.counts["G", 4] == pytest.approx(0.028499470, abs=5e-8)
        assert motif.counts["G", 5] == pytest.approx(0.028499470, abs=5e-8)
        assert motif.counts["G", 6] == pytest.approx(0.016304348, abs=5e-8)
        assert motif.counts["G", 7] == pytest.approx(0.160127253, abs=5e-8)
        assert motif.counts["A", 0] == pytest.approx(0.218451750, abs=5e-8)
        assert motif.counts["A", 1] == pytest.approx(0.023064687, abs=5e-8)
        assert motif.counts["A", 2] == pytest.approx(0.656680806, abs=5e-8)
        assert motif.counts["A", 3] == pytest.approx(0.898197243, abs=5e-8)
        assert motif.counts["A", 4] == pytest.approx(0.040694592, abs=5e-8)
        assert motif.counts["A", 5] == pytest.approx(0.132953340, abs=5e-8)
        assert motif.counts["A", 6] == pytest.approx(0.749072110, abs=5e-8)
        assert motif.counts["A", 7] == pytest.approx(0.628313892, abs=5e-8)
        assert motif.counts["T", 0] == pytest.approx(0.235949099, abs=5e-8)
        assert motif.counts["T", 1] == pytest.approx(0.590402969, abs=5e-8)
        assert motif.counts["T", 2] == pytest.approx(0.010869565, abs=5e-8)
        assert motif.counts["T", 3] == pytest.approx(0.033934252, abs=5e-8)
        assert motif.counts["T", 4] == pytest.approx(0.880567338, abs=5e-8)
        assert motif.counts["T", 5] == pytest.approx(0.797852598, abs=5e-8)
        assert motif.counts["T", 6] == pytest.approx(0.206124072, abs=5e-8)
        assert motif.counts["T", 7] == pytest.approx(0.177624602, abs=5e-8)
        assert motif.counts["C", 0] == pytest.approx(0.089607635, abs=5e-8)
        assert motif.counts["C", 1] == pytest.approx(0.317338282, abs=5e-8)
        assert motif.counts["C", 2] == pytest.approx(0.321580064, abs=5e-8)
        assert motif.counts["C", 3] == pytest.approx(0.046129374, abs=5e-8)
        assert motif.counts["C", 4] == pytest.approx(0.050238600, abs=5e-8)
        assert motif.counts["C", 5] == pytest.approx(0.040694592, abs=5e-8)
        assert motif.counts["C", 6] == pytest.approx(0.028499470, abs=5e-8)
        assert motif.counts["C", 7] == pytest.approx(0.033934252, abs=5e-8)
        assert motif.consensus == "GTAATTAA"
        assert motif.degenerate_consensus == "NYAATTAA"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.2005361303021225,
                        0.6336277209668335,
                        0.933405467206956,
                        1.3704286046679186,
                        1.2873833086962072,
                        1.0187720746919493,
                        0.975022432438911,
                        0.547109562258496,
                    ]
                ),
            )
        assert motif[:-2].consensus == "GTAATT"
        motif = record[6]
        assert motif.name == "MA0001.1 AGL3"
        assert motif.length == 10
        assert motif.alphabet == "GATC"
        assert motif.counts["G", 0] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["G", 1] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["G", 2] == pytest.approx(3.0, abs=5e-8)
        assert motif.counts["G", 3] == pytest.approx(4.0, abs=5e-8)
        assert motif.counts["G", 4] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["G", 5] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["G", 6] == pytest.approx(5.0, abs=5e-8)
        assert motif.counts["G", 7] == pytest.approx(3.0, abs=5e-8)
        assert motif.counts["G", 8] == pytest.approx(28.0, abs=5e-8)
        assert motif.counts["G", 9] == pytest.approx(88.0, abs=5e-8)
        assert motif.counts["A", 0] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["A", 1] == pytest.approx(3.0, abs=5e-8)
        assert motif.counts["A", 2] == pytest.approx(79.0, abs=5e-8)
        assert motif.counts["A", 3] == pytest.approx(40.0, abs=5e-8)
        assert motif.counts["A", 4] == pytest.approx(66.0, abs=5e-8)
        assert motif.counts["A", 5] == pytest.approx(48.0, abs=5e-8)
        assert motif.counts["A", 6] == pytest.approx(65.0, abs=5e-8)
        assert motif.counts["A", 7] == pytest.approx(11.0, abs=5e-8)
        assert motif.counts["A", 8] == pytest.approx(65.0, abs=5e-8)
        assert motif.counts["A", 9] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["T", 0] == pytest.approx(2.0, abs=5e-8)
        assert motif.counts["T", 1] == pytest.approx(19.0, abs=5e-8)
        assert motif.counts["T", 2] == pytest.approx(11.0, abs=5e-8)
        assert motif.counts["T", 3] == pytest.approx(50.0, abs=5e-8)
        assert motif.counts["T", 4] == pytest.approx(29.0, abs=5e-8)
        assert motif.counts["T", 5] == pytest.approx(47.0, abs=5e-8)
        assert motif.counts["T", 6] == pytest.approx(22.0, abs=5e-8)
        assert motif.counts["T", 7] == pytest.approx(81.0, abs=5e-8)
        assert motif.counts["T", 8] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["T", 9] == pytest.approx(6.0, abs=5e-8)
        assert motif.counts["C", 0] == pytest.approx(94.0, abs=5e-8)
        assert motif.counts["C", 1] == pytest.approx(75.0, abs=5e-8)
        assert motif.counts["C", 2] == pytest.approx(4.0, abs=5e-8)
        assert motif.counts["C", 3] == pytest.approx(3.0, abs=5e-8)
        assert motif.counts["C", 4] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["C", 5] == pytest.approx(2.0, abs=5e-8)
        assert motif.counts["C", 6] == pytest.approx(5.0, abs=5e-8)
        assert motif.counts["C", 7] == pytest.approx(2.0, abs=5e-8)
        assert motif.counts["C", 8] == pytest.approx(3.0, abs=5e-8)
        assert motif.counts["C", 9] == pytest.approx(3.0, abs=5e-8)
        assert motif.consensus == "CCATAAATAG"
        assert motif.degenerate_consensus == "CCAWAWATAG"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        1.7725753233561499,
                        1.0972718180683638,
                        1.0578945228970464,
                        0.6353945886004412,
                        0.9651537633423314,
                        0.8757972203228152,
                        0.6864859661195083,
                        1.1561334005018244,
                        0.8724039945822116,
                        1.4691041160249607,
                    ]
                ),
            )
        assert motif[:-2].consensus == "CCATAAAT"
        motif = record[7]
        assert motif.name == "MA0001.1 AGL3"
        assert motif.length == 10
        assert motif.alphabet == "GATC"
        assert motif.counts["G", 0] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["G", 1] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["G", 2] == pytest.approx(3.0, abs=5e-8)
        assert motif.counts["G", 3] == pytest.approx(4.0, abs=5e-8)
        assert motif.counts["G", 4] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["G", 5] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["G", 6] == pytest.approx(5.0, abs=5e-8)
        assert motif.counts["G", 7] == pytest.approx(3.0, abs=5e-8)
        assert motif.counts["G", 8] == pytest.approx(28.0, abs=5e-8)
        assert motif.counts["G", 9] == pytest.approx(88.0, abs=5e-8)
        assert motif.counts["A", 0] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["A", 1] == pytest.approx(3.0, abs=5e-8)
        assert motif.counts["A", 2] == pytest.approx(79.0, abs=5e-8)
        assert motif.counts["A", 3] == pytest.approx(40.0, abs=5e-8)
        assert motif.counts["A", 4] == pytest.approx(66.0, abs=5e-8)
        assert motif.counts["A", 5] == pytest.approx(48.0, abs=5e-8)
        assert motif.counts["A", 6] == pytest.approx(65.0, abs=5e-8)
        assert motif.counts["A", 7] == pytest.approx(11.0, abs=5e-8)
        assert motif.counts["A", 8] == pytest.approx(65.0, abs=5e-8)
        assert motif.counts["A", 9] == pytest.approx(0.0, abs=5e-8)
        assert motif.counts["T", 0] == pytest.approx(2.0, abs=5e-8)
        assert motif.counts["T", 1] == pytest.approx(19.0, abs=5e-8)
        assert motif.counts["T", 2] == pytest.approx(11.0, abs=5e-8)
        assert motif.counts["T", 3] == pytest.approx(50.0, abs=5e-8)
        assert motif.counts["T", 4] == pytest.approx(29.0, abs=5e-8)
        assert motif.counts["T", 5] == pytest.approx(47.0, abs=5e-8)
        assert motif.counts["T", 6] == pytest.approx(22.0, abs=5e-8)
        assert motif.counts["T", 7] == pytest.approx(81.0, abs=5e-8)
        assert motif.counts["T", 8] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["T", 9] == pytest.approx(6.0, abs=5e-8)
        assert motif.counts["C", 0] == pytest.approx(94.0, abs=5e-8)
        assert motif.counts["C", 1] == pytest.approx(75.0, abs=5e-8)
        assert motif.counts["C", 2] == pytest.approx(4.0, abs=5e-8)
        assert motif.counts["C", 3] == pytest.approx(3.0, abs=5e-8)
        assert motif.counts["C", 4] == pytest.approx(1.0, abs=5e-8)
        assert motif.counts["C", 5] == pytest.approx(2.0, abs=5e-8)
        assert motif.counts["C", 6] == pytest.approx(5.0, abs=5e-8)
        assert motif.counts["C", 7] == pytest.approx(2.0, abs=5e-8)
        assert motif.counts["C", 8] == pytest.approx(3.0, abs=5e-8)
        assert motif.counts["C", 9] == pytest.approx(3.0, abs=5e-8)
        assert motif.consensus == "CCATAAATAG"
        assert motif.degenerate_consensus == "CCAWAWATAG"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        1.7725753233561499,
                        1.0972718180683638,
                        1.0578945228970464,
                        0.6353945886004412,
                        0.9651537633423314,
                        0.8757972203228152,
                        0.6864859661195083,
                        1.1561334005018244,
                        0.8724039945822116,
                        1.4691041160249607,
                    ]
                ),
            )
        assert motif[:-2].consensus == "CCATAAAT"
        motif = record[8]
        assert motif.name == ""
        assert motif.length == 9
        assert motif.alphabet == "GATC"
        assert motif.counts["G", 0] == pytest.approx(0.016, abs=5e-8)
        assert motif.counts["G", 1] == pytest.approx(0.020, abs=5e-8)
        assert motif.counts["G", 2] == pytest.approx(0.028, abs=5e-8)
        assert motif.counts["G", 3] == pytest.approx(0.016, abs=5e-8)
        assert motif.counts["G", 4] == pytest.approx(0.020, abs=5e-8)
        assert motif.counts["G", 5] == pytest.approx(0.028, abs=5e-8)
        assert motif.counts["G", 6] == pytest.approx(0.047, abs=5e-8)
        assert motif.counts["G", 7] == pytest.approx(0.045, abs=5e-8)
        assert motif.counts["G", 8] == pytest.approx(0.216, abs=5e-8)
        assert motif.counts["A", 0] == pytest.approx(0.116, abs=5e-8)
        assert motif.counts["A", 1] == pytest.approx(0.974, abs=5e-8)
        assert motif.counts["A", 2] == pytest.approx(0.444, abs=5e-8)
        assert motif.counts["A", 3] == pytest.approx(0.116, abs=5e-8)
        assert motif.counts["A", 4] == pytest.approx(0.974, abs=5e-8)
        assert motif.counts["A", 5] == pytest.approx(0.444, abs=5e-8)
        assert motif.counts["A", 6] == pytest.approx(0.667, abs=5e-8)
        assert motif.counts["A", 7] == pytest.approx(0.939, abs=5e-8)
        assert motif.counts["A", 8] == pytest.approx(0.068, abs=5e-8)
        assert motif.counts["T", 0] == pytest.approx(0.150, abs=5e-8)
        assert motif.counts["T", 1] == pytest.approx(0.001, abs=5e-8)
        assert motif.counts["T", 2] == pytest.approx(0.314, abs=5e-8)
        assert motif.counts["T", 3] == pytest.approx(0.150, abs=5e-8)
        assert motif.counts["T", 4] == pytest.approx(0.001, abs=5e-8)
        assert motif.counts["T", 5] == pytest.approx(0.314, abs=5e-8)
        assert motif.counts["T", 6] == pytest.approx(0.143, abs=5e-8)
        assert motif.counts["T", 7] == pytest.approx(0.009, abs=5e-8)
        assert motif.counts["T", 8] == pytest.approx(0.609, abs=5e-8)
        assert motif.counts["C", 0] == pytest.approx(0.718, abs=5e-8)
        assert motif.counts["C", 1] == pytest.approx(0.006, abs=5e-8)
        assert motif.counts["C", 2] == pytest.approx(0.214, abs=5e-8)
        assert motif.counts["C", 3] == pytest.approx(0.718, abs=5e-8)
        assert motif.counts["C", 4] == pytest.approx(0.006, abs=5e-8)
        assert motif.counts["C", 5] == pytest.approx(0.214, abs=5e-8)
        assert motif.counts["C", 6] == pytest.approx(0.143, abs=5e-8)
        assert motif.counts["C", 7] == pytest.approx(0.006, abs=5e-8)
        assert motif.counts["C", 8] == pytest.approx(0.107, abs=5e-8)
        assert motif.consensus == "CAACAAAAT"
        assert motif.degenerate_consensus == "CAWCAWAAT"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.79033346,
                        1.79461597,
                        0.33472715,
                        0.79033346,
                        1.79461597,
                        0.33472715,
                        0.60049374,
                        1.60901246,
                        0.47798759,
                    ]
                ),
            )
        assert motif[:-2].consensus == "CAACAAA"

    def test_sites_parsing(self):
        """Test if Bio.motifs can parse JASPAR-style sites files."""
        with open("motifs/Arnt.sites") as stream:
            m = motifs.read(stream, "sites")
        assert m.length == 6
        assert m.alignment.sequences[0] == "CACGTG"
        assert m.alignment.sequences[1] == "CACGTG"
        assert m.alignment.sequences[2] == "CACGTG"
        assert m.alignment.sequences[3] == "CACGTG"
        assert m.alignment.sequences[4] == "CACGTG"
        assert m.alignment.sequences[5] == "CACGTG"
        assert m.alignment.sequences[6] == "CACGTG"
        assert m.alignment.sequences[7] == "CACGTG"
        assert m.alignment.sequences[8] == "CACGTG"
        assert m.alignment.sequences[9] == "CACGTG"
        assert m.alignment.sequences[10] == "CACGTG"
        assert m.alignment.sequences[11] == "CACGTG"
        assert m.alignment.sequences[12] == "CACGTG"
        assert m.alignment.sequences[13] == "CACGTG"
        assert m.alignment.sequences[14] == "CACGTG"
        assert m.alignment.sequences[15] == "AACGTG"
        assert m.alignment.sequences[16] == "AACGTG"
        assert m.alignment.sequences[17] == "AACGTG"
        assert m.alignment.sequences[18] == "AACGTG"
        assert m.alignment.sequences[19] == "CGCGTG"
        assert m.counts["A", 0] == pytest.approx(4, abs=5e-8)
        assert m.counts["A", 1] == pytest.approx(19, abs=5e-8)
        assert m.counts["A", 2] == pytest.approx(0, abs=5e-8)
        assert m.counts["A", 3] == pytest.approx(0, abs=5e-8)
        assert m.counts["A", 4] == pytest.approx(0, abs=5e-8)
        assert m.counts["A", 5] == pytest.approx(0, abs=5e-8)
        assert m.counts["C", 0] == pytest.approx(16, abs=5e-8)
        assert m.counts["C", 1] == pytest.approx(0, abs=5e-8)
        assert m.counts["C", 2] == pytest.approx(20, abs=5e-8)
        assert m.counts["C", 3] == pytest.approx(0, abs=5e-8)
        assert m.counts["C", 4] == pytest.approx(0, abs=5e-8)
        assert m.counts["C", 5] == pytest.approx(0, abs=5e-8)
        assert m.counts["G", 0] == pytest.approx(0, abs=5e-8)
        assert m.counts["G", 1] == pytest.approx(1, abs=5e-8)
        assert m.counts["G", 2] == pytest.approx(0, abs=5e-8)
        assert m.counts["G", 3] == pytest.approx(20, abs=5e-8)
        assert m.counts["G", 4] == pytest.approx(0, abs=5e-8)
        assert m.counts["G", 5] == pytest.approx(20, abs=5e-8)
        assert m.counts["T", 0] == pytest.approx(0, abs=5e-8)
        assert m.counts["T", 1] == pytest.approx(0, abs=5e-8)
        assert m.counts["T", 2] == pytest.approx(0, abs=5e-8)
        assert m.counts["T", 3] == pytest.approx(0, abs=5e-8)
        assert m.counts["T", 4] == pytest.approx(20, abs=5e-8)
        assert m.counts["T", 5] == pytest.approx(0, abs=5e-8)
        assert m.consensus == "CACGTG"
        assert m.degenerate_consensus == "CACGTG"
        assert np.allclose(
                m.relative_entropy,
                np.array([1.278071905112638, 1.7136030428840439, 2.0, 2.0, 2.0, 2.0]),
            )
        assert m[::2].consensus == "CCT"


class TestMEME(unittest.TestCase):
    def test_meme_parser_1(self):
        """Parse motifs/meme.INO_up800.classic.oops.xml file."""
        with open("motifs/meme.INO_up800.classic.oops.xml") as stream:
            record = motifs.parse(stream, "meme")
        assert record.version == "5.0.1"
        assert record.datafile == "common/INO_up800.s"
        assert record.alphabet == "ACGT"
        assert len(record.sequences) == 7
        assert record.sequences[0] == "sequence_0"
        assert record.sequences[1] == "sequence_1"
        assert record.sequences[2] == "sequence_2"
        assert record.sequences[3] == "sequence_3"
        assert record.sequences[4] == "sequence_4"
        assert record.sequences[5] == "sequence_5"
        assert record.sequences[6] == "sequence_6"
        assert record.command == "meme common/INO_up800.s -oc results/meme10 -mod oops -dna -revcomp -bfile common/yeast.nc.6.freq -nmotifs 2 -objfun classic -minw 8 -nostatus "
        assert len(record) == 2
        motif = record[0]
        assert motif.name == "GSKGCATGTGAAA"
        assert record["GSKGCATGTGAAA"] == motif
        assert motif.num_occurrences == 7
        assert motif.evalue == pytest.approx(0.19, abs=5e-8)
        assert motif.alphabet == "ACGT"
        assert len(motif.alignment.sequences) == 7
        assert motif.alignment.sequences[0].pvalue == pytest.approx(1.21e-08, abs=5e-11)
        assert motif.alignment.sequences[1].pvalue == pytest.approx(1.87e-08, abs=5e-11)
        assert motif.alignment.sequences[2].pvalue == pytest.approx(6.62e-08, abs=5e-11)
        assert motif.alignment.sequences[3].pvalue == pytest.approx(1.05e-07, abs=5e-10)
        assert motif.alignment.sequences[4].pvalue == pytest.approx(1.69e-07, abs=5e-10)
        assert motif.alignment.sequences[5].pvalue == pytest.approx(5.62e-07, abs=5e-10)
        assert motif.alignment.sequences[6].pvalue == pytest.approx(1.08e-06, abs=5e-09)
        assert motif.alignment.sequences[0].sequence_name == "INO1"
        assert motif.alignment.sequences[1].sequence_name == "FAS1"
        assert motif.alignment.sequences[2].sequence_name == "ACC1"
        assert motif.alignment.sequences[3].sequence_name == "CHO2"
        assert motif.alignment.sequences[4].sequence_name == "CHO1"
        assert motif.alignment.sequences[5].sequence_name == "FAS2"
        assert motif.alignment.sequences[6].sequence_name == "OPI3"
        assert motif.alignment.sequences[0].sequence_id == "sequence_5"
        assert motif.alignment.sequences[1].sequence_id == "sequence_2"
        assert motif.alignment.sequences[2].sequence_id == "sequence_4"
        assert motif.alignment.sequences[3].sequence_id == "sequence_1"
        assert motif.alignment.sequences[4].sequence_id == "sequence_0"
        assert motif.alignment.sequences[5].sequence_id == "sequence_3"
        assert motif.alignment.sequences[6].sequence_id == "sequence_6"
        assert motif.alignment.sequences[0].strand == "+"
        assert motif.alignment.sequences[1].strand == "-"
        assert motif.alignment.sequences[2].strand == "-"
        assert motif.alignment.sequences[3].strand == "-"
        assert motif.alignment.sequences[4].strand == "-"
        assert motif.alignment.sequences[5].strand == "-"
        assert motif.alignment.sequences[6].strand == "+"
        assert motif.alignment.sequences[0].length == 13
        assert motif.alignment.sequences[1].length == 13
        assert motif.alignment.sequences[2].length == 13
        assert motif.alignment.sequences[3].length == 13
        assert motif.alignment.sequences[4].length == 13
        assert motif.alignment.sequences[5].length == 13
        assert motif.alignment.sequences[6].length == 13
        assert motif.alignment.sequences[0].start == 620
        assert motif.alignment.sequences[1].start == 94
        assert motif.alignment.sequences[2].start == 82
        assert motif.alignment.sequences[3].start == 353
        assert motif.alignment.sequences[4].start == 639
        assert motif.alignment.sequences[5].start == 566
        assert motif.alignment.sequences[6].start == 585
        assert motif.alignment.sequences[0] == "GCGGCATGTGAAA"
        assert motif.alignment.sequences[1] == "GCGGCATGTGAAG"
        assert motif.alignment.sequences[2] == "GGGCCATGTGAAG"
        assert motif.alignment.sequences[3] == "GCGGCATGAGAAA"
        assert motif.alignment.sequences[4] == "GGTCCATGTGAAA"
        assert motif.alignment.sequences[5] == "GTAGCATGTGAAA"
        assert motif.alignment.sequences[6] == "AGTGCATGTGGAA"
        assert motif.consensus == "GCGGCATGTGAAA"
        assert motif.degenerate_consensus == "GSKGCATGTGAAA"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        1.4083272214176723,
                        0.5511843642748154,
                        0.6212165065138244,
                        1.136879431433369,
                        2.0,
                        2.0,
                        2.0,
                        2.0,
                        1.4083272214176723,
                        2.0,
                        1.4083272214176723,
                        2.0,
                        1.136879431433369,
                    ]
                ),
            )
        assert motif[1::2].consensus == "CGAGGA"
        motif = record[1]
        assert motif.name == "TTGACWCYTGCYCWG"
        assert record["TTGACWCYTGCYCWG"] == motif
        assert motif.num_occurrences == 7
        assert motif.evalue == pytest.approx(54, abs=5e-8)
        assert motif.alphabet == "ACGT"
        assert len(motif.alignment.sequences) == 7
        assert motif.alignment.sequences[0].pvalue == pytest.approx(7.2e-10, abs=5e-12)
        assert motif.alignment.sequences[1].pvalue == pytest.approx(2.56e-08, abs=5e-11)
        assert motif.alignment.sequences[2].pvalue == pytest.approx(1.59e-07, abs=5e-10)
        assert motif.alignment.sequences[3].pvalue == pytest.approx(2.05e-07, abs=5e-10)
        assert motif.alignment.sequences[4].pvalue == pytest.approx(3.85e-07, abs=5e-10)
        assert motif.alignment.sequences[5].pvalue == pytest.approx(5.11e-07, abs=5e-10)
        assert motif.alignment.sequences[6].pvalue == pytest.approx(8.01e-07, abs=5e-10)
        assert motif.alignment.sequences[0].sequence_id == "sequence_1"
        assert motif.alignment.sequences[1].sequence_id == "sequence_6"
        assert motif.alignment.sequences[2].sequence_id == "sequence_4"
        assert motif.alignment.sequences[3].sequence_id == "sequence_0"
        assert motif.alignment.sequences[4].sequence_id == "sequence_2"
        assert motif.alignment.sequences[5].sequence_id == "sequence_3"
        assert motif.alignment.sequences[6].sequence_id == "sequence_5"
        assert motif.alignment.sequences[0].strand == "+"
        assert motif.alignment.sequences[1].strand == "-"
        assert motif.alignment.sequences[2].strand == "-"
        assert motif.alignment.sequences[3].strand == "+"
        assert motif.alignment.sequences[4].strand == "+"
        assert motif.alignment.sequences[5].strand == "-"
        assert motif.alignment.sequences[6].strand == "+"
        assert motif.alignment.sequences[0].length == 15
        assert motif.alignment.sequences[1].length == 15
        assert motif.alignment.sequences[2].length == 15
        assert motif.alignment.sequences[3].length == 15
        assert motif.alignment.sequences[4].length == 15
        assert motif.alignment.sequences[5].length == 15
        assert motif.alignment.sequences[6].length == 15
        assert motif.alignment.sequences[0].start == 104
        assert motif.alignment.sequences[1].start == 566
        assert motif.alignment.sequences[2].start == 585
        assert motif.alignment.sequences[3].start == 30
        assert motif.alignment.sequences[4].start == 54
        assert motif.alignment.sequences[5].start == 272
        assert motif.alignment.sequences[6].start == 214
        assert motif.alignment.sequences[0] == "TTGACACCTGCCCAG"
        assert motif.alignment.sequences[1] == "TTGACACCTACCCTG"
        assert motif.alignment.sequences[2] == "TTGTCTCTTGCTCTG"
        assert motif.alignment.sequences[3] == "TTGACACTTGATCAG"
        assert motif.alignment.sequences[4] == "TTCACTACTCCCCTG"
        assert motif.alignment.sequences[5] == "TTGACAACGGCTGGG"
        assert motif.alignment.sequences[6] == "TTCACGCTTGCTACG"
        assert motif.consensus == "TTGACACCTGCTCTG"
        assert motif.degenerate_consensus == "TTGACWCYTGCYCNG"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        2.0,
                        2.0,
                        1.136879431433369,
                        1.4083272214176723,
                        2.0,
                        0.6212165065138244,
                        1.136879431433369,
                        1.0147718639657484,
                        1.4083272214176723,
                        0.8511651457190834,
                        1.4083272214176723,
                        1.0147718639657484,
                        0.8511651457190834,
                        0.15762900682289133,
                        2.0,
                    ]
                ),
            )
        assert motif[1::2].consensus == "TAACGTT"

    def test_meme_parser_2(self):
        """Parsing motifs/meme.adh.classic.oops.xml file."""
        with open("motifs/meme.adh.classic.oops.xml") as stream:
            record = motifs.parse(stream, "meme")
        assert record.version == "5.0.1"
        assert record.datafile == "common/adh.s"
        assert record.alphabet == "ACDEFGHIKLMNPQRSTVWY"
        assert len(record.sequences) == 33
        assert record.sequences[0] == "sequence_0"
        assert record.sequences[1] == "sequence_1"
        assert record.sequences[2] == "sequence_2"
        assert record.sequences[3] == "sequence_3"
        assert record.sequences[4] == "sequence_4"
        assert record.sequences[5] == "sequence_5"
        assert record.sequences[6] == "sequence_6"
        assert record.sequences[7] == "sequence_7"
        assert record.sequences[8] == "sequence_8"
        assert record.sequences[9] == "sequence_9"
        assert record.sequences[10] == "sequence_10"
        assert record.sequences[11] == "sequence_11"
        assert record.sequences[12] == "sequence_12"
        assert record.sequences[13] == "sequence_13"
        assert record.sequences[14] == "sequence_14"
        assert record.sequences[15] == "sequence_15"
        assert record.sequences[16] == "sequence_16"
        assert record.sequences[17] == "sequence_17"
        assert record.sequences[18] == "sequence_18"
        assert record.sequences[19] == "sequence_19"
        assert record.sequences[20] == "sequence_20"
        assert record.sequences[21] == "sequence_21"
        assert record.sequences[22] == "sequence_22"
        assert record.sequences[23] == "sequence_23"
        assert record.sequences[24] == "sequence_24"
        assert record.sequences[25] == "sequence_25"
        assert record.sequences[26] == "sequence_26"
        assert record.sequences[27] == "sequence_27"
        assert record.sequences[28] == "sequence_28"
        assert record.sequences[29] == "sequence_29"
        assert record.sequences[30] == "sequence_30"
        assert record.sequences[31] == "sequence_31"
        assert record.sequences[32] == "sequence_32"
        assert record.command == "meme common/adh.s -oc results/meme4 -mod oops -protein -nmotifs 2 -objfun classic -minw 8 -nostatus "
        assert len(record) == 2
        motif = record[0]
        assert motif.id == "motif_1"
        assert motif.name == "GKVALVTGAASGJGKATAKAL"
        assert motif.alt_id == "MEME-1"
        assert record["GKVALVTGAASGJGKATAKAL"] == motif
        assert motif.num_occurrences == 33
        assert motif.evalue == pytest.approx(4.0e-129, abs=5e-131)
        assert motif.alphabet == "ACDEFGHIKLMNPQRSTVWY"
        assert len(motif.alignment.sequences) == 33
        assert motif.alignment.sequences[0].pvalue == pytest.approx(8.78e-18, abs=5e-21)
        assert motif.alignment.sequences[1].pvalue == pytest.approx(1.41e-17, abs=5e-20)
        assert motif.alignment.sequences[2].pvalue == pytest.approx(1.42e-16, abs=5e-19)
        assert motif.alignment.sequences[3].pvalue == pytest.approx(2.75e-16, abs=5e-19)
        assert motif.alignment.sequences[4].pvalue == pytest.approx(3.55e-16, abs=5e-19)
        assert motif.alignment.sequences[5].pvalue == pytest.approx(3.55e-16, abs=5e-19)
        assert motif.alignment.sequences[6].pvalue == pytest.approx(1.74e-15, abs=5e-18)
        assert motif.alignment.sequences[7].pvalue == pytest.approx(3.87e-15, abs=5e-18)
        assert motif.alignment.sequences[8].pvalue == pytest.approx(4.84e-15, abs=5e-18)
        assert motif.alignment.sequences[9].pvalue == pytest.approx(1.04e-14, abs=5e-17)
        assert motif.alignment.sequences[10].pvalue == pytest.approx(1.58e-14, abs=5e-17)
        assert motif.alignment.sequences[11].pvalue == pytest.approx(1.76e-14, abs=5e-17)
        assert motif.alignment.sequences[12].pvalue == pytest.approx(2.16e-14, abs=5e-17)
        assert motif.alignment.sequences[13].pvalue == pytest.approx(2.94e-14, abs=5e-17)
        assert motif.alignment.sequences[14].pvalue == pytest.approx(3.25e-14, abs=5e-17)
        assert motif.alignment.sequences[15].pvalue == pytest.approx(3.98e-14, abs=5e-17)
        assert motif.alignment.sequences[16].pvalue == pytest.approx(4.39e-14, abs=5e-17)
        assert motif.alignment.sequences[17].pvalue == pytest.approx(4.39e-14, abs=5e-17)
        assert motif.alignment.sequences[18].pvalue == pytest.approx(4.85e-14, abs=5e-17)
        assert motif.alignment.sequences[19].pvalue == pytest.approx(6.52e-14, abs=5e-17)
        assert motif.alignment.sequences[20].pvalue == pytest.approx(1.41e-13, abs=5e-16)
        assert motif.alignment.sequences[21].pvalue == pytest.approx(1.55e-13, abs=5e-16)
        assert motif.alignment.sequences[22].pvalue == pytest.approx(3.07e-12, abs=5e-15)
        assert motif.alignment.sequences[23].pvalue == pytest.approx(5.43e-12, abs=5e-15)
        assert motif.alignment.sequences[24].pvalue == pytest.approx(6.91e-12, abs=5e-15)
        assert motif.alignment.sequences[25].pvalue == pytest.approx(8.76e-12, abs=5e-15)
        assert motif.alignment.sequences[26].pvalue == pytest.approx(9.48e-12, abs=5e-15)
        assert motif.alignment.sequences[27].pvalue == pytest.approx(1.2e-11, abs=5e-13)
        assert motif.alignment.sequences[28].pvalue == pytest.approx(1.19e-09, abs=5e-12)
        assert motif.alignment.sequences[29].pvalue == pytest.approx(1.54e-09, abs=5e-12)
        assert motif.alignment.sequences[30].pvalue == pytest.approx(1.99e-09, abs=5e-12)
        assert motif.alignment.sequences[31].pvalue == pytest.approx(1.42e-06, abs=5e-09)
        assert motif.alignment.sequences[32].pvalue == pytest.approx(3.43e-06, abs=5e-09)
        assert motif.alignment.sequences[0].sequence_name == "BUDC_KLETE"
        assert motif.alignment.sequences[1].sequence_name == "YINL_LISMO"
        assert motif.alignment.sequences[2].sequence_name == "DHII_HUMAN"
        assert motif.alignment.sequences[3].sequence_name == "HDE_CANTR"
        assert motif.alignment.sequences[4].sequence_name == "YRTP_BACSU"
        assert motif.alignment.sequences[5].sequence_name == "ENTA_ECOLI"
        assert motif.alignment.sequences[6].sequence_name == "HDHA_ECOLI"
        assert motif.alignment.sequences[7].sequence_name == "RIDH_KLEAE"
        assert motif.alignment.sequences[8].sequence_name == "DHB2_HUMAN"
        assert motif.alignment.sequences[9].sequence_name == "FIXR_BRAJA"
        assert motif.alignment.sequences[10].sequence_name == "PCR_PEA"
        assert motif.alignment.sequences[11].sequence_name == "DHCA_HUMAN"
        assert motif.alignment.sequences[12].sequence_name == "BDH_HUMAN"
        assert motif.alignment.sequences[13].sequence_name == "3BHD_COMTE"
        assert motif.alignment.sequences[14].sequence_name == "DHGB_BACME"
        assert motif.alignment.sequences[15].sequence_name == "DHMA_FLAS1"
        assert motif.alignment.sequences[16].sequence_name == "FVT1_HUMAN"
        assert motif.alignment.sequences[17].sequence_name == "BA72_EUBSP"
        assert motif.alignment.sequences[18].sequence_name == "BPHB_PSEPS"
        assert motif.alignment.sequences[19].sequence_name == "DHB3_HUMAN"
        assert motif.alignment.sequences[20].sequence_name == "DHES_HUMAN"
        assert motif.alignment.sequences[21].sequence_name == "AP27_MOUSE"
        assert motif.alignment.sequences[22].sequence_name == "2BHD_STREX"
        assert motif.alignment.sequences[23].sequence_name == "NODG_RHIME"
        assert motif.alignment.sequences[24].sequence_name == "HMTR_LEIMA"
        assert motif.alignment.sequences[25].sequence_name == "LIGD_PSEPA"
        assert motif.alignment.sequences[26].sequence_name == "MAS1_AGRRA"
        assert motif.alignment.sequences[27].sequence_name == "RFBB_NEIGO"
        assert motif.alignment.sequences[28].sequence_name == "GUTD_ECOLI"
        assert motif.alignment.sequences[29].sequence_name == "ADH_DROME"
        assert motif.alignment.sequences[30].sequence_name == "FABI_ECOLI"
        assert motif.alignment.sequences[31].sequence_name == "CSGA_MYXXA"
        assert motif.alignment.sequences[32].sequence_name == "YURA_MYXXA"
        assert motif.alignment.sequences[0].sequence_id == "sequence_7"
        assert motif.alignment.sequences[1].sequence_id == "sequence_20"
        assert motif.alignment.sequences[2].sequence_id == "sequence_10"
        assert motif.alignment.sequences[3].sequence_id == "sequence_15"
        assert motif.alignment.sequences[4].sequence_id == "sequence_21"
        assert motif.alignment.sequences[5].sequence_id == "sequence_12"
        assert motif.alignment.sequences[6].sequence_id == "sequence_16"
        assert motif.alignment.sequences[7].sequence_id == "sequence_19"
        assert motif.alignment.sequences[8].sequence_id == "sequence_23"
        assert motif.alignment.sequences[9].sequence_id == "sequence_13"
        assert motif.alignment.sequences[10].sequence_id == "sequence_30"
        assert motif.alignment.sequences[11].sequence_id == "sequence_25"
        assert motif.alignment.sequences[12].sequence_id == "sequence_5"
        assert motif.alignment.sequences[13].sequence_id == "sequence_1"
        assert motif.alignment.sequences[14].sequence_id == "sequence_9"
        assert motif.alignment.sequences[15].sequence_id == "sequence_11"
        assert motif.alignment.sequences[16].sequence_id == "sequence_27"
        assert motif.alignment.sequences[17].sequence_id == "sequence_4"
        assert motif.alignment.sequences[18].sequence_id == "sequence_6"
        assert motif.alignment.sequences[19].sequence_id == "sequence_24"
        assert motif.alignment.sequences[20].sequence_id == "sequence_8"
        assert motif.alignment.sequences[21].sequence_id == "sequence_3"
        assert motif.alignment.sequences[22].sequence_id == "sequence_0"
        assert motif.alignment.sequences[23].sequence_id == "sequence_18"
        assert motif.alignment.sequences[24].sequence_id == "sequence_28"
        assert motif.alignment.sequences[25].sequence_id == "sequence_17"
        assert motif.alignment.sequences[26].sequence_id == "sequence_29"
        assert motif.alignment.sequences[27].sequence_id == "sequence_31"
        assert motif.alignment.sequences[28].sequence_id == "sequence_14"
        assert motif.alignment.sequences[29].sequence_id == "sequence_2"
        assert motif.alignment.sequences[30].sequence_id == "sequence_26"
        assert motif.alignment.sequences[31].sequence_id == "sequence_22"
        assert motif.alignment.sequences[32].sequence_id == "sequence_32"
        assert motif.alignment.sequences[0].strand == "+"
        assert motif.alignment.sequences[1].strand == "+"
        assert motif.alignment.sequences[2].strand == "+"
        assert motif.alignment.sequences[3].strand == "+"
        assert motif.alignment.sequences[4].strand == "+"
        assert motif.alignment.sequences[5].strand == "+"
        assert motif.alignment.sequences[6].strand == "+"
        assert motif.alignment.sequences[7].strand == "+"
        assert motif.alignment.sequences[8].strand == "+"
        assert motif.alignment.sequences[9].strand == "+"
        assert motif.alignment.sequences[10].strand == "+"
        assert motif.alignment.sequences[11].strand == "+"
        assert motif.alignment.sequences[12].strand == "+"
        assert motif.alignment.sequences[13].strand == "+"
        assert motif.alignment.sequences[14].strand == "+"
        assert motif.alignment.sequences[15].strand == "+"
        assert motif.alignment.sequences[16].strand == "+"
        assert motif.alignment.sequences[17].strand == "+"
        assert motif.alignment.sequences[18].strand == "+"
        assert motif.alignment.sequences[19].strand == "+"
        assert motif.alignment.sequences[20].strand == "+"
        assert motif.alignment.sequences[21].strand == "+"
        assert motif.alignment.sequences[22].strand == "+"
        assert motif.alignment.sequences[23].strand == "+"
        assert motif.alignment.sequences[24].strand == "+"
        assert motif.alignment.sequences[25].strand == "+"
        assert motif.alignment.sequences[26].strand == "+"
        assert motif.alignment.sequences[27].strand == "+"
        assert motif.alignment.sequences[28].strand == "+"
        assert motif.alignment.sequences[29].strand == "+"
        assert motif.alignment.sequences[30].strand == "+"
        assert motif.alignment.sequences[31].strand == "+"
        assert motif.alignment.sequences[32].strand == "+"
        assert motif.alignment.sequences[0].length == 21
        assert motif.alignment.sequences[1].length == 21
        assert motif.alignment.sequences[2].length == 21
        assert motif.alignment.sequences[3].length == 21
        assert motif.alignment.sequences[4].length == 21
        assert motif.alignment.sequences[5].length == 21
        assert motif.alignment.sequences[6].length == 21
        assert motif.alignment.sequences[7].length == 21
        assert motif.alignment.sequences[8].length == 21
        assert motif.alignment.sequences[9].length == 21
        assert motif.alignment.sequences[10].length == 21
        assert motif.alignment.sequences[11].length == 21
        assert motif.alignment.sequences[12].length == 21
        assert motif.alignment.sequences[13].length == 21
        assert motif.alignment.sequences[14].length == 21
        assert motif.alignment.sequences[15].length == 21
        assert motif.alignment.sequences[16].length == 21
        assert motif.alignment.sequences[17].length == 21
        assert motif.alignment.sequences[18].length == 21
        assert motif.alignment.sequences[19].length == 21
        assert motif.alignment.sequences[20].length == 21
        assert motif.alignment.sequences[21].length == 21
        assert motif.alignment.sequences[22].length == 21
        assert motif.alignment.sequences[23].length == 21
        assert motif.alignment.sequences[24].length == 21
        assert motif.alignment.sequences[25].length == 21
        assert motif.alignment.sequences[26].length == 21
        assert motif.alignment.sequences[27].length == 21
        assert motif.alignment.sequences[28].length == 21
        assert motif.alignment.sequences[29].length == 21
        assert motif.alignment.sequences[30].length == 21
        assert motif.alignment.sequences[31].length == 21
        assert motif.alignment.sequences[32].length == 21
        assert motif.alignment.sequences[0].start == 2
        assert motif.alignment.sequences[1].start == 5
        assert motif.alignment.sequences[2].start == 34
        assert motif.alignment.sequences[3].start == 322
        assert motif.alignment.sequences[4].start == 6
        assert motif.alignment.sequences[5].start == 5
        assert motif.alignment.sequences[6].start == 11
        assert motif.alignment.sequences[7].start == 14
        assert motif.alignment.sequences[8].start == 82
        assert motif.alignment.sequences[9].start == 36
        assert motif.alignment.sequences[10].start == 86
        assert motif.alignment.sequences[11].start == 4
        assert motif.alignment.sequences[12].start == 55
        assert motif.alignment.sequences[13].start == 6
        assert motif.alignment.sequences[14].start == 7
        assert motif.alignment.sequences[15].start == 14
        assert motif.alignment.sequences[16].start == 32
        assert motif.alignment.sequences[17].start == 6
        assert motif.alignment.sequences[18].start == 5
        assert motif.alignment.sequences[19].start == 48
        assert motif.alignment.sequences[20].start == 2
        assert motif.alignment.sequences[21].start == 7
        assert motif.alignment.sequences[22].start == 6
        assert motif.alignment.sequences[23].start == 6
        assert motif.alignment.sequences[24].start == 6
        assert motif.alignment.sequences[25].start == 6
        assert motif.alignment.sequences[26].start == 245
        assert motif.alignment.sequences[27].start == 6
        assert motif.alignment.sequences[28].start == 2
        assert motif.alignment.sequences[29].start == 6
        assert motif.alignment.sequences[30].start == 6
        assert motif.alignment.sequences[31].start == 13
        assert motif.alignment.sequences[32].start == 116
        assert motif.alignment.sequences[0] == "QKVALVTGAGQGIGKAIALRL"
        assert motif.alignment.sequences[1] == "NKVIIITGASSGIGKATALLL"
        assert motif.alignment.sequences[2] == "GKKVIVTGASKGIGREMAYHL"
        assert motif.alignment.sequences[3] == "DKVVLITGAGAGLGKEYAKWF"
        assert motif.alignment.sequences[4] == "HKTALITGGGRGIGRATALAL"
        assert motif.alignment.sequences[5] == "GKNVWVTGAGKGIGYATALAF"
        assert motif.alignment.sequences[6] == "GKCAIITGAGAGIGKEIAITF"
        assert motif.alignment.sequences[7] == "GKVAAITGAASGIGLECARTL"
        assert motif.alignment.sequences[8] == "QKAVLVTGGDCGLGHALCKYL"
        assert motif.alignment.sequences[9] == "PKVMLLTGASRGIGHATAKLF"
        assert motif.alignment.sequences[10] == "KGNVVITGASSGLGLATAKAL"
        assert motif.alignment.sequences[11] == "IHVALVTGGNKGIGLAIVRDL"
        assert motif.alignment.sequences[12] == "SKAVLVTGCDSGFGFSLAKHL"
        assert motif.alignment.sequences[13] == "GKVALVTGGASGVGLEVVKLL"
        assert motif.alignment.sequences[14] == "GKVVVITGSSTGLGKSMAIRF"
        assert motif.alignment.sequences[15] == "GKAAIVTGAAGGIGRATVEAY"
        assert motif.alignment.sequences[16] == "GAHVVVTGGSSGIGKCIAIEC"
        assert motif.alignment.sequences[17] == "DKVTIITGGTRGIGFAAAKIF"
        assert motif.alignment.sequences[18] == "GEAVLITGGASGLGRALVDRF"
        assert motif.alignment.sequences[19] == "GQWAVITGAGDGIGKAYSFEL"
        assert motif.alignment.sequences[20] == "RTVVLITGCSSGIGLHLAVRL"
        assert motif.alignment.sequences[21] == "GLRALVTGAGKGIGRDTVKAL"
        assert motif.alignment.sequences[22] == "GKTVIITGGARGLGAEAARQA"
        assert motif.alignment.sequences[23] == "GRKALVTGASGAIGGAIARVL"
        assert motif.alignment.sequences[24] == "VPVALVTGAAKRLGRSIAEGL"
        assert motif.alignment.sequences[25] == "DQVAFITGGASGAGFGQAKVF"
        assert motif.alignment.sequences[26] == "SPVILVSGSNRGVGKAIAEDL"
        assert motif.alignment.sequences[27] == "KKNILVTGGAGFIGSAVVRHI"
        assert motif.alignment.sequences[28] == "NQVAVVIGGGQTLGAFLCHGL"
        assert motif.alignment.sequences[29] == "NKNVIFVAGLGGIGLDTSKEL"
        assert motif.alignment.sequences[30] == "GKRILVTGVASKLSIAYGIAQ"
        assert motif.alignment.sequences[31] == "VDVLINNAGVSGLWCALGDVD"
        assert motif.alignment.sequences[32] == "IIDTNVTGAAATLSAVLPQMV"
        assert motif.consensus == "GKVALVTGAASGIGKATAKAL"
        assert motif[2:8].consensus == "VALVTG"
        motif = record[1]
        assert motif.name == "VGNPGASAYSASKAAVRGLTESLALELAP"
        assert motif.alt_id == "MEME-2"
        assert record["VGNPGASAYSASKAAVRGLTESLALELAP"] == motif
        assert motif.num_occurrences == 33
        assert motif.evalue == pytest.approx(3.1e-130, abs=5e-132)
        assert motif.alphabet == "ACDEFGHIKLMNPQRSTVWY"
        assert len(motif.alignment.sequences) == 33
        assert motif.alignment.sequences[0].pvalue == pytest.approx(2.09e-21, abs=5e-24)
        assert motif.alignment.sequences[1].pvalue == pytest.approx(7.63e-20, abs=5e-23)
        assert motif.alignment.sequences[2].pvalue == pytest.approx(6.49e-19, abs=5e-22)
        assert motif.alignment.sequences[3].pvalue == pytest.approx(1.92e-18, abs=5e-21)
        assert motif.alignment.sequences[4].pvalue == pytest.approx(5.46e-18, abs=5e-21)
        assert motif.alignment.sequences[5].pvalue == pytest.approx(6.21e-18, abs=5e-21)
        assert motif.alignment.sequences[6].pvalue == pytest.approx(4.52e-17, abs=5e-20)
        assert motif.alignment.sequences[7].pvalue == pytest.approx(4.52e-17, abs=5e-20)
        assert motif.alignment.sequences[8].pvalue == pytest.approx(9.21e-17, abs=5e-20)
        assert motif.alignment.sequences[9].pvalue == pytest.approx(1.65e-16, abs=5e-19)
        assert motif.alignment.sequences[10].pvalue == pytest.approx(2.07e-16, abs=5e-19)
        assert motif.alignment.sequences[11].pvalue == pytest.approx(3.65e-16, abs=5e-19)
        assert motif.alignment.sequences[12].pvalue == pytest.approx(5.7e-16, abs=5e-18)
        assert motif.alignment.sequences[13].pvalue == pytest.approx(5.7e-16, abs=5e-18)
        assert motif.alignment.sequences[14].pvalue == pytest.approx(7.93e-16, abs=5e-19)
        assert motif.alignment.sequences[15].pvalue == pytest.approx(8.85e-16, abs=5e-19)
        assert motif.alignment.sequences[16].pvalue == pytest.approx(1.1e-15, abs=5e-17)
        assert motif.alignment.sequences[17].pvalue == pytest.approx(1.69e-15, abs=5e-18)
        assert motif.alignment.sequences[18].pvalue == pytest.approx(3.54e-15, abs=5e-18)
        assert motif.alignment.sequences[19].pvalue == pytest.approx(4.83e-15, abs=5e-18)
        assert motif.alignment.sequences[20].pvalue == pytest.approx(7.27e-15, abs=5e-18)
        assert motif.alignment.sequences[21].pvalue == pytest.approx(9.85e-15, abs=5e-18)
        assert motif.alignment.sequences[22].pvalue == pytest.approx(2.41e-14, abs=5e-17)
        assert motif.alignment.sequences[23].pvalue == pytest.approx(2.66e-14, abs=5e-17)
        assert motif.alignment.sequences[24].pvalue == pytest.approx(1.22e-13, abs=5e-16)
        assert motif.alignment.sequences[25].pvalue == pytest.approx(5.18e-13, abs=5e-16)
        assert motif.alignment.sequences[26].pvalue == pytest.approx(1.24e-12, abs=5e-15)
        assert motif.alignment.sequences[27].pvalue == pytest.approx(1.35e-12, abs=5e-15)
        assert motif.alignment.sequences[28].pvalue == pytest.approx(5.59e-12, abs=5e-15)
        assert motif.alignment.sequences[29].pvalue == pytest.approx(1.44e-10, abs=5e-13)
        assert motif.alignment.sequences[30].pvalue == pytest.approx(1.61e-08, abs=5e-11)
        assert motif.alignment.sequences[31].pvalue == pytest.approx(4.26e-08, abs=5e-11)
        assert motif.alignment.sequences[32].pvalue == pytest.approx(1.16e-07, abs=5e-10)
        assert motif.alignment.sequences[0].sequence_name == "BUDC_KLETE"
        assert motif.alignment.sequences[1].sequence_name == "NODG_RHIME"
        assert motif.alignment.sequences[2].sequence_name == "FVT1_HUMAN"
        assert motif.alignment.sequences[3].sequence_name == "DHES_HUMAN"
        assert motif.alignment.sequences[4].sequence_name == "DHB3_HUMAN"
        assert motif.alignment.sequences[5].sequence_name == "YRTP_BACSU"
        assert motif.alignment.sequences[6].sequence_name == "HMTR_LEIMA"
        assert motif.alignment.sequences[7].sequence_name == "HDE_CANTR"
        assert motif.alignment.sequences[8].sequence_name == "DHGB_BACME"
        assert motif.alignment.sequences[9].sequence_name == "GUTD_ECOLI"
        assert motif.alignment.sequences[10].sequence_name == "3BHD_COMTE"
        assert motif.alignment.sequences[11].sequence_name == "DHII_HUMAN"
        assert motif.alignment.sequences[12].sequence_name == "BPHB_PSEPS"
        assert motif.alignment.sequences[13].sequence_name == "AP27_MOUSE"
        assert motif.alignment.sequences[14].sequence_name == "BDH_HUMAN"
        assert motif.alignment.sequences[15].sequence_name == "YINL_LISMO"
        assert motif.alignment.sequences[16].sequence_name == "FIXR_BRAJA"
        assert motif.alignment.sequences[17].sequence_name == "2BHD_STREX"
        assert motif.alignment.sequences[18].sequence_name == "RFBB_NEIGO"
        assert motif.alignment.sequences[19].sequence_name == "YURA_MYXXA"
        assert motif.alignment.sequences[20].sequence_name == "RIDH_KLEAE"
        assert motif.alignment.sequences[21].sequence_name == "DHMA_FLAS1"
        assert motif.alignment.sequences[22].sequence_name == "DHB2_HUMAN"
        assert motif.alignment.sequences[23].sequence_name == "HDHA_ECOLI"
        assert motif.alignment.sequences[24].sequence_name == "ENTA_ECOLI"
        assert motif.alignment.sequences[25].sequence_name == "LIGD_PSEPA"
        assert motif.alignment.sequences[26].sequence_name == "CSGA_MYXXA"
        assert motif.alignment.sequences[27].sequence_name == "BA72_EUBSP"
        assert motif.alignment.sequences[28].sequence_name == "ADH_DROME"
        assert motif.alignment.sequences[29].sequence_name == "MAS1_AGRRA"
        assert motif.alignment.sequences[30].sequence_name == "PCR_PEA"
        assert motif.alignment.sequences[31].sequence_name == "FABI_ECOLI"
        assert motif.alignment.sequences[32].sequence_name == "DHCA_HUMAN"
        assert motif.alignment.sequences[0].sequence_id == "sequence_7"
        assert motif.alignment.sequences[1].sequence_id == "sequence_18"
        assert motif.alignment.sequences[2].sequence_id == "sequence_27"
        assert motif.alignment.sequences[3].sequence_id == "sequence_8"
        assert motif.alignment.sequences[4].sequence_id == "sequence_24"
        assert motif.alignment.sequences[5].sequence_id == "sequence_21"
        assert motif.alignment.sequences[6].sequence_id == "sequence_28"
        assert motif.alignment.sequences[7].sequence_id == "sequence_15"
        assert motif.alignment.sequences[8].sequence_id == "sequence_9"
        assert motif.alignment.sequences[9].sequence_id == "sequence_14"
        assert motif.alignment.sequences[10].sequence_id == "sequence_1"
        assert motif.alignment.sequences[11].sequence_id == "sequence_10"
        assert motif.alignment.sequences[12].sequence_id == "sequence_6"
        assert motif.alignment.sequences[13].sequence_id == "sequence_3"
        assert motif.alignment.sequences[14].sequence_id == "sequence_5"
        assert motif.alignment.sequences[15].sequence_id == "sequence_20"
        assert motif.alignment.sequences[16].sequence_id == "sequence_13"
        assert motif.alignment.sequences[17].sequence_id == "sequence_0"
        assert motif.alignment.sequences[18].sequence_id == "sequence_31"
        assert motif.alignment.sequences[19].sequence_id == "sequence_32"
        assert motif.alignment.sequences[20].sequence_id == "sequence_19"
        assert motif.alignment.sequences[21].sequence_id == "sequence_11"
        assert motif.alignment.sequences[22].sequence_id == "sequence_23"
        assert motif.alignment.sequences[23].sequence_id == "sequence_16"
        assert motif.alignment.sequences[24].sequence_id == "sequence_12"
        assert motif.alignment.sequences[25].sequence_id == "sequence_17"
        assert motif.alignment.sequences[26].sequence_id == "sequence_22"
        assert motif.alignment.sequences[27].sequence_id == "sequence_4"
        assert motif.alignment.sequences[28].sequence_id == "sequence_2"
        assert motif.alignment.sequences[29].sequence_id == "sequence_29"
        assert motif.alignment.sequences[30].sequence_id == "sequence_30"
        assert motif.alignment.sequences[31].sequence_id == "sequence_26"
        assert motif.alignment.sequences[32].sequence_id == "sequence_25"
        assert motif.alignment.sequences[0].start == 144
        assert motif.alignment.sequences[1].start == 144
        assert motif.alignment.sequences[2].start == 178
        assert motif.alignment.sequences[3].start == 147
        assert motif.alignment.sequences[4].start == 190
        assert motif.alignment.sequences[5].start == 147
        assert motif.alignment.sequences[6].start == 185
        assert motif.alignment.sequences[7].start == 459
        assert motif.alignment.sequences[8].start == 152
        assert motif.alignment.sequences[9].start == 146
        assert motif.alignment.sequences[10].start == 143
        assert motif.alignment.sequences[11].start == 175
        assert motif.alignment.sequences[12].start == 145
        assert motif.alignment.sequences[13].start == 141
        assert motif.alignment.sequences[14].start == 200
        assert motif.alignment.sequences[15].start == 146
        assert motif.alignment.sequences[16].start == 181
        assert motif.alignment.sequences[17].start == 144
        assert motif.alignment.sequences[18].start == 157
        assert motif.alignment.sequences[19].start == 152
        assert motif.alignment.sequences[20].start == 152
        assert motif.alignment.sequences[21].start == 157
        assert motif.alignment.sequences[22].start == 224
        assert motif.alignment.sequences[23].start == 151
        assert motif.alignment.sequences[24].start == 136
        assert motif.alignment.sequences[25].start == 149
        assert motif.alignment.sequences[26].start == 80
        assert motif.alignment.sequences[27].start == 149
        assert motif.alignment.sequences[28].start == 144
        assert motif.alignment.sequences[29].start == 384
        assert motif.alignment.sequences[30].start == 18
        assert motif.alignment.sequences[31].start == 177
        assert motif.alignment.sequences[32].start == 144
        assert motif.alignment.sequences[0] == "VGNPELAVYSSSKFAVRGLTQTAARDLAP"
        assert motif.alignment.sequences[1] == "IGNPGQTNYCASKAGMIGFSKSLAQEIAT"
        assert motif.alignment.sequences[2] == "LGLFGFTAYSASKFAIRGLAEALQMEVKP"
        assert motif.alignment.sequences[3] == "MGLPFNDVYCASKFALEGLCESLAVLLLP"
        assert motif.alignment.sequences[4] == "FPWPLYSMYSASKAFVCAFSKALQEEYKA"
        assert motif.alignment.sequences[5] == "RGAAVTSAYSASKFAVLGLTESLMQEVRK"
        assert motif.alignment.sequences[6] == "QPLLGYTIYTMAKGALEGLTRSAALELAP"
        assert motif.alignment.sequences[7] == "YGNFGQANYSSSKAGILGLSKTMAIEGAK"
        assert motif.alignment.sequences[8] == "IPWPLFVHYAASKGGMKLMTETLALEYAP"
        assert motif.alignment.sequences[9] == "VGSKHNSGYSAAKFGGVGLTQSLALDLAE"
        assert motif.alignment.sequences[10] == "LPIEQYAGYSASKAAVSALTRAAALSCRK"
        assert motif.alignment.sequences[11] == "VAYPMVAAYSASKFALDGFFSSIRKEYSV"
        assert motif.alignment.sequences[12] == "YPNGGGPLYTAAKQAIVGLVRELAFELAP"
        assert motif.alignment.sequences[13] == "VTFPNLITYSSTKGAMTMLTKAMAMELGP"
        assert motif.alignment.sequences[14] == "MANPARSPYCITKFGVEAFSDCLRYEMYP"
        assert motif.alignment.sequences[15] == "KAYPGGAVYGATKWAVRDLMEVLRMESAQ"
        assert motif.alignment.sequences[16] == "VHPFAGSAYATSKAALASLTRELAHDYAP"
        assert motif.alignment.sequences[17] == "MGLALTSSYGASKWGVRGLSKLAAVELGT"
        assert motif.alignment.sequences[18] == "TPYAPSSPYSASKAAADHLVRAWQRTYRL"
        assert motif.alignment.sequences[19] == "FRGLPATRYSASKAFLSTFMESLRVDLRG"
        assert motif.alignment.sequences[20] == "VPVIWEPVYTASKFAVQAFVHTTRRQVAQ"
        assert motif.alignment.sequences[21] == "MAEPEAAAYVAAKGGVAMLTRAMAVDLAR"
        assert motif.alignment.sequences[22] == "APMERLASYGSSKAAVTMFSSVMRLELSK"
        assert motif.alignment.sequences[23] == "NKNINMTSYASSKAAASHLVRNMAFDLGE"
        assert motif.alignment.sequences[24] == "TPRIGMSAYGASKAALKSLALSVGLELAG"
        assert motif.alignment.sequences[25] == "MGSALAGPYSAAKAASINLMEGYRQGLEK"
        assert motif.alignment.sequences[26] == "NTDGGAYAYRMSKAALNMAVRSMSTDLRP"
        assert motif.alignment.sequences[27] == "FGSLSGVGYPASKASVIGLTHGLGREIIR"
        assert motif.alignment.sequences[28] == "NAIYQVPVYSGTKAAVVNFTSSLAKLAPI"
        assert motif.alignment.sequences[29] == "RVLNPLVGYNMTKHALGGLTKTTQHVGWD"
        assert motif.alignment.sequences[30] == "EGKIGASLKDSTLFGVSSLSDSLKGDFTS"
        assert motif.alignment.sequences[31] == "MGPEGVRVNAISAGPIRTLAASGIKDFRK"
        assert motif.alignment.sequences[32] == "RALKSCSPELQQKFRSETITEEELVGLMN"
        assert motif.consensus == "MGLPGASAYSASKAAVRGLTESLALELAP"
        assert motif[-8:-2].consensus == "SLALEL"

    def test_meme_parser_3(self):
        """Parse motifs/meme.farntrans5.classic.anr.xml file."""
        with open("motifs/meme.farntrans5.classic.anr.xml") as stream:
            record = motifs.parse(stream, "meme")
        assert record.version == "5.0.1"
        assert record.datafile == "common/farntrans5.s"
        assert record.alphabet == "ACDEFGHIKLMNPQRSTVWY"
        assert len(record.sequences) == 5
        assert record.sequences[0] == "sequence_0"
        assert record.sequences[1] == "sequence_1"
        assert record.sequences[2] == "sequence_2"
        assert record.sequences[3] == "sequence_3"
        assert record.sequences[4] == "sequence_4"
        assert record.command == "meme common/farntrans5.s -oc results/meme15 -mod anr -protein -nmotifs 2 -objfun classic -minw 8 -nostatus "
        assert len(record) == 2
        motif = record[0]
        assert motif.name == "GGFGGRPGKEVDLCYTYCALAALAJLGSLD"
        assert record["GGFGGRPGKEVDLCYTYCALAALAJLGSLD"] == motif
        assert motif.num_occurrences == 24
        assert motif.evalue == pytest.approx(2.2e-94, abs=5e-96)
        assert motif.alphabet == "ACDEFGHIKLMNPQRSTVWY"
        assert len(motif.alignment.sequences) == 24
        assert motif.alignment.sequences[0].pvalue == pytest.approx(6.98e-22, abs=5e-25)
        assert motif.alignment.sequences[1].pvalue == pytest.approx(4.67e-21, abs=5e-24)
        assert motif.alignment.sequences[2].pvalue == pytest.approx(1.25e-19, abs=5e-22)
        assert motif.alignment.sequences[3].pvalue == pytest.approx(1.56e-19, abs=5e-22)
        assert motif.alignment.sequences[4].pvalue == pytest.approx(2.44e-19, abs=5e-22)
        assert motif.alignment.sequences[5].pvalue == pytest.approx(6.47e-19, abs=5e-22)
        assert motif.alignment.sequences[6].pvalue == pytest.approx(8.9e-19, abs=5e-21)
        assert motif.alignment.sequences[7].pvalue == pytest.approx(2.53e-18, abs=5e-21)
        assert motif.alignment.sequences[8].pvalue == pytest.approx(1.27e-17, abs=5e-20)
        assert motif.alignment.sequences[9].pvalue == pytest.approx(2.77e-17, abs=5e-20)
        assert motif.alignment.sequences[10].pvalue == pytest.approx(4.93e-17, abs=5e-20)
        assert motif.alignment.sequences[11].pvalue == pytest.approx(7.19e-17, abs=5e-20)
        assert motif.alignment.sequences[12].pvalue == pytest.approx(8.68e-17, abs=5e-20)
        assert motif.alignment.sequences[13].pvalue == pytest.approx(2.62e-16, abs=5e-19)
        assert motif.alignment.sequences[14].pvalue == pytest.approx(2.87e-16, abs=5e-19)
        assert motif.alignment.sequences[15].pvalue == pytest.approx(7.66e-15, abs=5e-18)
        assert motif.alignment.sequences[16].pvalue == pytest.approx(2.21e-14, abs=5e-17)
        assert motif.alignment.sequences[17].pvalue == pytest.approx(3.29e-14, abs=5e-17)
        assert motif.alignment.sequences[18].pvalue == pytest.approx(7.21e-14, abs=5e-17)
        assert motif.alignment.sequences[19].pvalue == pytest.approx(1.14e-13, abs=5e-16)
        assert motif.alignment.sequences[20].pvalue == pytest.approx(1.67e-13, abs=5e-16)
        assert motif.alignment.sequences[21].pvalue == pytest.approx(4.42e-13, abs=5e-16)
        assert motif.alignment.sequences[22].pvalue == pytest.approx(5.11e-13, abs=5e-16)
        assert motif.alignment.sequences[23].pvalue == pytest.approx(2.82e-10, abs=5e-13)
        assert motif.alignment.sequences[0].sequence_name == "BET2_YEAST"
        assert motif.alignment.sequences[1].sequence_name == "RATRABGERB"
        assert motif.alignment.sequences[2].sequence_name == "CAL1_YEAST"
        assert motif.alignment.sequences[3].sequence_name == "PFTB_RAT"
        assert motif.alignment.sequences[4].sequence_name == "PFTB_RAT"
        assert motif.alignment.sequences[5].sequence_name == "RATRABGERB"
        assert motif.alignment.sequences[6].sequence_name == "RATRABGERB"
        assert motif.alignment.sequences[7].sequence_name == "BET2_YEAST"
        assert motif.alignment.sequences[8].sequence_name == "RATRABGERB"
        assert motif.alignment.sequences[9].sequence_name == "BET2_YEAST"
        assert motif.alignment.sequences[10].sequence_name == "RAM1_YEAST"
        assert motif.alignment.sequences[11].sequence_name == "BET2_YEAST"
        assert motif.alignment.sequences[12].sequence_name == "RAM1_YEAST"
        assert motif.alignment.sequences[13].sequence_name == "PFTB_RAT"
        assert motif.alignment.sequences[14].sequence_name == "RAM1_YEAST"
        assert motif.alignment.sequences[15].sequence_name == "PFTB_RAT"
        assert motif.alignment.sequences[16].sequence_name == "RATRABGERB"
        assert motif.alignment.sequences[17].sequence_name == "PFTB_RAT"
        assert motif.alignment.sequences[18].sequence_name == "BET2_YEAST"
        assert motif.alignment.sequences[19].sequence_name == "CAL1_YEAST"
        assert motif.alignment.sequences[20].sequence_name == "RAM1_YEAST"
        assert motif.alignment.sequences[21].sequence_name == "CAL1_YEAST"
        assert motif.alignment.sequences[22].sequence_name == "RAM1_YEAST"
        assert motif.alignment.sequences[23].sequence_name == "BET2_YEAST"
        assert motif.alignment.sequences[0].sequence_id == "sequence_2"
        assert motif.alignment.sequences[1].sequence_id == "sequence_3"
        assert motif.alignment.sequences[2].sequence_id == "sequence_4"
        assert motif.alignment.sequences[3].sequence_id == "sequence_1"
        assert motif.alignment.sequences[4].sequence_id == "sequence_1"
        assert motif.alignment.sequences[5].sequence_id == "sequence_3"
        assert motif.alignment.sequences[6].sequence_id == "sequence_3"
        assert motif.alignment.sequences[7].sequence_id == "sequence_2"
        assert motif.alignment.sequences[8].sequence_id == "sequence_3"
        assert motif.alignment.sequences[9].sequence_id == "sequence_2"
        assert motif.alignment.sequences[10].sequence_id == "sequence_0"
        assert motif.alignment.sequences[11].sequence_id == "sequence_2"
        assert motif.alignment.sequences[12].sequence_id == "sequence_0"
        assert motif.alignment.sequences[13].sequence_id == "sequence_1"
        assert motif.alignment.sequences[14].sequence_id == "sequence_0"
        assert motif.alignment.sequences[15].sequence_id == "sequence_1"
        assert motif.alignment.sequences[16].sequence_id == "sequence_3"
        assert motif.alignment.sequences[17].sequence_id == "sequence_1"
        assert motif.alignment.sequences[18].sequence_id == "sequence_2"
        assert motif.alignment.sequences[19].sequence_id == "sequence_4"
        assert motif.alignment.sequences[20].sequence_id == "sequence_0"
        assert motif.alignment.sequences[21].sequence_id == "sequence_4"
        assert motif.alignment.sequences[22].sequence_id == "sequence_0"
        assert motif.alignment.sequences[23].sequence_id == "sequence_2"
        assert motif.alignment.sequences[0].strand == "+"
        assert motif.alignment.sequences[1].strand == "+"
        assert motif.alignment.sequences[2].strand == "+"
        assert motif.alignment.sequences[3].strand == "+"
        assert motif.alignment.sequences[4].strand == "+"
        assert motif.alignment.sequences[5].strand == "+"
        assert motif.alignment.sequences[6].strand == "+"
        assert motif.alignment.sequences[7].strand == "+"
        assert motif.alignment.sequences[8].strand == "+"
        assert motif.alignment.sequences[9].strand == "+"
        assert motif.alignment.sequences[10].strand == "+"
        assert motif.alignment.sequences[11].strand == "+"
        assert motif.alignment.sequences[12].strand == "+"
        assert motif.alignment.sequences[13].strand == "+"
        assert motif.alignment.sequences[14].strand == "+"
        assert motif.alignment.sequences[15].strand == "+"
        assert motif.alignment.sequences[16].strand == "+"
        assert motif.alignment.sequences[17].strand == "+"
        assert motif.alignment.sequences[18].strand == "+"
        assert motif.alignment.sequences[19].strand == "+"
        assert motif.alignment.sequences[20].strand == "+"
        assert motif.alignment.sequences[21].strand == "+"
        assert motif.alignment.sequences[22].strand == "+"
        assert motif.alignment.sequences[23].strand == "+"
        assert motif.alignment.sequences[0].length == 30
        assert motif.alignment.sequences[1].length == 30
        assert motif.alignment.sequences[2].length == 30
        assert motif.alignment.sequences[3].length == 30
        assert motif.alignment.sequences[4].length == 30
        assert motif.alignment.sequences[5].length == 30
        assert motif.alignment.sequences[6].length == 30
        assert motif.alignment.sequences[7].length == 30
        assert motif.alignment.sequences[8].length == 30
        assert motif.alignment.sequences[9].length == 30
        assert motif.alignment.sequences[10].length == 30
        assert motif.alignment.sequences[11].length == 30
        assert motif.alignment.sequences[12].length == 30
        assert motif.alignment.sequences[13].length == 30
        assert motif.alignment.sequences[14].length == 30
        assert motif.alignment.sequences[15].length == 30
        assert motif.alignment.sequences[16].length == 30
        assert motif.alignment.sequences[17].length == 30
        assert motif.alignment.sequences[18].length == 30
        assert motif.alignment.sequences[19].length == 30
        assert motif.alignment.sequences[20].length == 30
        assert motif.alignment.sequences[21].length == 30
        assert motif.alignment.sequences[22].length == 30
        assert motif.alignment.sequences[23].length == 30
        assert motif.alignment.sequences[0].start == 223
        assert motif.alignment.sequences[1].start == 227
        assert motif.alignment.sequences[2].start == 275
        assert motif.alignment.sequences[3].start == 237
        assert motif.alignment.sequences[4].start == 138
        assert motif.alignment.sequences[5].start == 179
        assert motif.alignment.sequences[6].start == 131
        assert motif.alignment.sequences[7].start == 172
        assert motif.alignment.sequences[8].start == 276
        assert motif.alignment.sequences[9].start == 124
        assert motif.alignment.sequences[10].start == 247
        assert motif.alignment.sequences[11].start == 272
        assert motif.alignment.sequences[12].start == 145
        assert motif.alignment.sequences[13].start == 286
        assert motif.alignment.sequences[14].start == 296
        assert motif.alignment.sequences[15].start == 348
        assert motif.alignment.sequences[16].start == 83
        assert motif.alignment.sequences[17].start == 189
        assert motif.alignment.sequences[18].start == 73
        assert motif.alignment.sequences[19].start == 205
        assert motif.alignment.sequences[20].start == 198
        assert motif.alignment.sequences[21].start == 327
        assert motif.alignment.sequences[22].start == 349
        assert motif.alignment.sequences[23].start == 24
        assert motif.alignment.sequences[0] == "GGLNGRPSKLPDVCYSWWVLSSLAIIGRLD"
        assert motif.alignment.sequences[1] == "GGLNGRPEKLPDVCYSWWVLASLKIIGRLH"
        assert motif.alignment.sequences[2] == "GGFQGRENKFADTCYAFWCLNSLHLLTKDW"
        assert motif.alignment.sequences[3] == "GGIGGVPGMEAHGGYTFCGLAALVILKKER"
        assert motif.alignment.sequences[4] == "GGFGGGPGQYPHLAPTYAAVNALCIIGTEE"
        assert motif.alignment.sequences[5] == "GGFGCRPGSESHAGQIYCCTGFLAITSQLH"
        assert motif.alignment.sequences[6] == "GSFAGDIWGEIDTRFSFCAVATLALLGKLD"
        assert motif.alignment.sequences[7] == "GGFGLCPNAESHAAQAFTCLGALAIANKLD"
        assert motif.alignment.sequences[8] == "GGFADRPGDMVDPFHTLFGIAGLSLLGEEQ"
        assert motif.alignment.sequences[9] == "GSFQGDRFGEVDTRFVYTALSALSILGELT"
        assert motif.alignment.sequences[10] == "GFGSCPHVDEAHGGYTFCATASLAILRSMD"
        assert motif.alignment.sequences[11] == "GGISDRPENEVDVFHTVFGVAGLSLMGYDN"
        assert motif.alignment.sequences[12] == "GPFGGGPGQLSHLASTYAAINALSLCDNID"
        assert motif.alignment.sequences[13] == "GGFQGRCNKLVDGCYSFWQAGLLPLLHRAL"
        assert motif.alignment.sequences[14] == "RGFCGRSNKLVDGCYSFWVGGSAAILEAFG"
        assert motif.alignment.sequences[15] == "GGLLDKPGKSRDFYHTCYCLSGLSIAQHFG"
        assert motif.alignment.sequences[16] == "GGVSASIGHDPHLLYTLSAVQILTLYDSIH"
        assert motif.alignment.sequences[17] == "GSFLMHVGGEVDVRSAYCAASVASLTNIIT"
        assert motif.alignment.sequences[18] == "GAFAPFPRHDAHLLTTLSAVQILATYDALD"
        assert motif.alignment.sequences[19] == "YNGAFGAHNEPHSGYTSCALSTLALLSSLE"
        assert motif.alignment.sequences[20] == "GFKTCLEVGEVDTRGIYCALSIATLLNILT"
        assert motif.alignment.sequences[21] == "GGFSKNDEEDADLYHSCLGSAALALIEGKF"
        assert motif.alignment.sequences[22] == "PGLRDKPGAHSDFYHTNYCLLGLAVAESSY"
        assert motif.alignment.sequences[23] == "HNFEYWLTEHLRLNGIYWGLTALCVLDSPE"
        assert motif.consensus == "GGFGGRPGKEVDLCYTFCALAALALLGSLD"
        assert motif[3:-8].consensus == "GGRPGKEVDLCYTFCALAA"
        motif = record[1]
        assert motif.name == "JNKEKLLEYILSCQ"
        assert record["JNKEKLLEYILSCQ"] == motif
        assert motif.num_occurrences == 21
        assert motif.evalue == pytest.approx(6.1e-21, abs=5e-23)
        assert motif.alphabet == "ACDEFGHIKLMNPQRSTVWY"
        assert len(motif.alignment.sequences) == 21
        assert motif.alignment.sequences[0].pvalue == pytest.approx(2.71e-12, abs=5e-15)
        assert motif.alignment.sequences[1].pvalue == pytest.approx(5.7e-12, abs=5e-14)
        assert motif.alignment.sequences[2].pvalue == pytest.approx(6.43e-12, abs=5e-15)
        assert motif.alignment.sequences[3].pvalue == pytest.approx(2.61e-11, abs=5e-14)
        assert motif.alignment.sequences[4].pvalue == pytest.approx(6.3e-11, abs=5e-13)
        assert motif.alignment.sequences[5].pvalue == pytest.approx(2.7e-10, abs=5e-12)
        assert motif.alignment.sequences[6].pvalue == pytest.approx(4.03e-10, abs=5e-13)
        assert motif.alignment.sequences[7].pvalue == pytest.approx(1.27e-09, abs=5e-12)
        assert motif.alignment.sequences[8].pvalue == pytest.approx(3.17e-09, abs=5e-12)
        assert motif.alignment.sequences[9].pvalue == pytest.approx(6.39e-09, abs=5e-12)
        assert motif.alignment.sequences[10].pvalue == pytest.approx(6.96e-09, abs=5e-12)
        assert motif.alignment.sequences[11].pvalue == pytest.approx(1.06e-08, abs=5e-11)
        assert motif.alignment.sequences[12].pvalue == pytest.approx(1.26e-08, abs=5e-11)
        assert motif.alignment.sequences[13].pvalue == pytest.approx(1.37e-08, abs=5e-11)
        assert motif.alignment.sequences[14].pvalue == pytest.approx(2.07e-08, abs=5e-11)
        assert motif.alignment.sequences[15].pvalue == pytest.approx(4.96e-08, abs=5e-11)
        assert motif.alignment.sequences[16].pvalue == pytest.approx(1.15e-07, abs=5e-10)
        assert motif.alignment.sequences[17].pvalue == pytest.approx(1.44e-07, abs=5e-10)
        assert motif.alignment.sequences[18].pvalue == pytest.approx(1.55e-07, abs=5e-10)
        assert motif.alignment.sequences[19].pvalue == pytest.approx(1.93e-07, abs=5e-10)
        assert motif.alignment.sequences[20].pvalue == pytest.approx(5.2e-07, abs=5e-09)
        assert motif.alignment.sequences[0].sequence_name == "RATRABGERB"
        assert motif.alignment.sequences[1].sequence_name == "BET2_YEAST"
        assert motif.alignment.sequences[2].sequence_name == "RATRABGERB"
        assert motif.alignment.sequences[3].sequence_name == "RATRABGERB"
        assert motif.alignment.sequences[4].sequence_name == "CAL1_YEAST"
        assert motif.alignment.sequences[5].sequence_name == "RAM1_YEAST"
        assert motif.alignment.sequences[6].sequence_name == "PFTB_RAT"
        assert motif.alignment.sequences[7].sequence_name == "RATRABGERB"
        assert motif.alignment.sequences[8].sequence_name == "BET2_YEAST"
        assert motif.alignment.sequences[9].sequence_name == "PFTB_RAT"
        assert motif.alignment.sequences[10].sequence_name == "RAM1_YEAST"
        assert motif.alignment.sequences[11].sequence_name == "CAL1_YEAST"
        assert motif.alignment.sequences[12].sequence_name == "PFTB_RAT"
        assert motif.alignment.sequences[13].sequence_name == "BET2_YEAST"
        assert motif.alignment.sequences[14].sequence_name == "RAM1_YEAST"
        assert motif.alignment.sequences[15].sequence_name == "RAM1_YEAST"
        assert motif.alignment.sequences[16].sequence_name == "RATRABGERB"
        assert motif.alignment.sequences[17].sequence_name == "RAM1_YEAST"
        assert motif.alignment.sequences[18].sequence_name == "PFTB_RAT"
        assert motif.alignment.sequences[19].sequence_name == "BET2_YEAST"
        assert motif.alignment.sequences[20].sequence_name == "CAL1_YEAST"
        assert motif.alignment.sequences[0].sequence_id == "sequence_3"
        assert motif.alignment.sequences[1].sequence_id == "sequence_2"
        assert motif.alignment.sequences[2].sequence_id == "sequence_3"
        assert motif.alignment.sequences[3].sequence_id == "sequence_3"
        assert motif.alignment.sequences[4].sequence_id == "sequence_4"
        assert motif.alignment.sequences[5].sequence_id == "sequence_0"
        assert motif.alignment.sequences[6].sequence_id == "sequence_1"
        assert motif.alignment.sequences[7].sequence_id == "sequence_3"
        assert motif.alignment.sequences[8].sequence_id == "sequence_2"
        assert motif.alignment.sequences[9].sequence_id == "sequence_1"
        assert motif.alignment.sequences[10].sequence_id == "sequence_0"
        assert motif.alignment.sequences[11].sequence_id == "sequence_4"
        assert motif.alignment.sequences[12].sequence_id == "sequence_1"
        assert motif.alignment.sequences[13].sequence_id == "sequence_2"
        assert motif.alignment.sequences[14].sequence_id == "sequence_0"
        assert motif.alignment.sequences[15].sequence_id == "sequence_0"
        assert motif.alignment.sequences[16].sequence_id == "sequence_3"
        assert motif.alignment.sequences[17].sequence_id == "sequence_0"
        assert motif.alignment.sequences[18].sequence_id == "sequence_1"
        assert motif.alignment.sequences[19].sequence_id == "sequence_2"
        assert motif.alignment.sequences[20].sequence_id == "sequence_4"
        assert motif.alignment.sequences[0].strand == "+"
        assert motif.alignment.sequences[1].strand == "+"
        assert motif.alignment.sequences[2].strand == "+"
        assert motif.alignment.sequences[3].strand == "+"
        assert motif.alignment.sequences[4].strand == "+"
        assert motif.alignment.sequences[5].strand == "+"
        assert motif.alignment.sequences[6].strand == "+"
        assert motif.alignment.sequences[7].strand == "+"
        assert motif.alignment.sequences[8].strand == "+"
        assert motif.alignment.sequences[9].strand == "+"
        assert motif.alignment.sequences[10].strand == "+"
        assert motif.alignment.sequences[11].strand == "+"
        assert motif.alignment.sequences[12].strand == "+"
        assert motif.alignment.sequences[13].strand == "+"
        assert motif.alignment.sequences[14].strand == "+"
        assert motif.alignment.sequences[15].strand == "+"
        assert motif.alignment.sequences[16].strand == "+"
        assert motif.alignment.sequences[17].strand == "+"
        assert motif.alignment.sequences[18].strand == "+"
        assert motif.alignment.sequences[19].strand == "+"
        assert motif.alignment.sequences[20].strand == "+"
        assert motif.alignment.sequences[0].length == 14
        assert motif.alignment.sequences[1].length == 14
        assert motif.alignment.sequences[2].length == 14
        assert motif.alignment.sequences[3].length == 14
        assert motif.alignment.sequences[4].length == 14
        assert motif.alignment.sequences[5].length == 14
        assert motif.alignment.sequences[6].length == 14
        assert motif.alignment.sequences[7].length == 14
        assert motif.alignment.sequences[8].length == 14
        assert motif.alignment.sequences[9].length == 14
        assert motif.alignment.sequences[10].length == 14
        assert motif.alignment.sequences[11].length == 14
        assert motif.alignment.sequences[12].length == 14
        assert motif.alignment.sequences[13].length == 14
        assert motif.alignment.sequences[14].length == 14
        assert motif.alignment.sequences[15].length == 14
        assert motif.alignment.sequences[16].length == 14
        assert motif.alignment.sequences[17].length == 14
        assert motif.alignment.sequences[18].length == 14
        assert motif.alignment.sequences[19].length == 14
        assert motif.alignment.sequences[20].length == 14
        assert motif.alignment.sequences[0].start == 66
        assert motif.alignment.sequences[1].start == 254
        assert motif.alignment.sequences[2].start == 258
        assert motif.alignment.sequences[3].start == 162
        assert motif.alignment.sequences[4].start == 190
        assert motif.alignment.sequences[5].start == 278
        assert motif.alignment.sequences[6].start == 172
        assert motif.alignment.sequences[7].start == 114
        assert motif.alignment.sequences[8].start == 7
        assert motif.alignment.sequences[9].start == 268
        assert motif.alignment.sequences[10].start == 414
        assert motif.alignment.sequences[11].start == 126
        assert motif.alignment.sequences[12].start == 220
        assert motif.alignment.sequences[13].start == 55
        assert motif.alignment.sequences[14].start == 229
        assert motif.alignment.sequences[15].start == 330
        assert motif.alignment.sequences[16].start == 18
        assert motif.alignment.sequences[17].start == 180
        assert motif.alignment.sequences[18].start == 73
        assert motif.alignment.sequences[19].start == 107
        assert motif.alignment.sequences[20].start == 36
        assert motif.alignment.sequences[0] == "MNKEEILVFIKSCQ"
        assert motif.alignment.sequences[1] == "INYEKLTEFILKCQ"
        assert motif.alignment.sequences[2] == "IDREKLRSFILACQ"
        assert motif.alignment.sequences[3] == "INVEKAIEFVLSCM"
        assert motif.alignment.sequences[4] == "IDTEKLLGYIMSQQ"
        assert motif.alignment.sequences[5] == "INVEKLLEWSSARQ"
        assert motif.alignment.sequences[6] == "INREKLLQYLYSLK"
        assert motif.alignment.sequences[7] == "INVDKVVAYVQSLQ"
        assert motif.alignment.sequences[8] == "LLKEKHIRYIESLD"
        assert motif.alignment.sequences[9] == "LNLKSLLQWVTSRQ"
        assert motif.alignment.sequences[10] == "ENVRKIIHYFKSNL"
        assert motif.alignment.sequences[11] == "LDKRSLARFVSKCQ"
        assert motif.alignment.sequences[12] == "DLFEGTAEWIARCQ"
        assert motif.alignment.sequences[13] == "FVKEEVISFVLSCW"
        assert motif.alignment.sequences[14] == "ELTEGVLNYLKNCQ"
        assert motif.alignment.sequences[15] == "FNKHALRDYILYCC"
        assert motif.alignment.sequences[16] == "LLLEKHADYIASYG"
        assert motif.alignment.sequences[17] == "IDRKGIYQWLISLK"
        assert motif.alignment.sequences[18] == "LQREKHFHYLKRGL"
        assert motif.alignment.sequences[19] == "DRKVRLISFIRGNQ"
        assert motif.alignment.sequences[20] == "VNRMAIIFYSISGL"
        assert motif.consensus == "INKEKLIEYILSCQ"
        assert motif[3:-8].consensus == "EKL"

    def test_minimal_meme_parser(self):
        """Parse motifs/minimal_test.meme file."""
        with open("motifs/minimal_test.meme") as stream:
            record = motifs.parse(stream, "minimal")
        assert record.version == "4"
        assert record.alphabet == "ACGT"
        assert len(record.sequences) == 0
        assert record.command == ""
        assert len(record) == 3
        motif = record[0]
        assert motif.name == "KRP"
        assert record["KRP"] == motif
        assert motif.num_occurrences == 17
        assert motif.length == 19
        assert motif.background["A"] == pytest.approx(0.30269730269730266, abs=5e-8)
        assert motif.background["C"] == pytest.approx(0.1828171828171828, abs=5e-8)
        assert motif.background["G"] == pytest.approx(0.20879120879120877, abs=5e-8)
        assert motif.background["T"] == pytest.approx(0.30569430569430567, abs=5e-8)
        assert motif.evalue == pytest.approx(4.1e-09, abs=5e-11)
        assert motif.alphabet == "ACGT"
        assert motif.alignment is None
        assert motif.consensus == "TGTGATCGAGGTCACACTT"
        assert motif.degenerate_consensus == "TGTGANNNWGNTCACAYWW"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        1.1684297174927525,
                        0.9432809925744818,
                        1.4307101633876265,
                        1.1549413780465179,
                        0.9308256303218774,
                        0.009164393966550805,
                        0.20124190687894253,
                        0.17618542656995528,
                        0.36777933103380855,
                        0.6635834532368525,
                        0.07729943368061855,
                        0.9838293592717438,
                        1.72489868427398,
                        0.8397561713453014,
                        1.72489868427398,
                        0.8455332015343343,
                        0.3106481207768122,
                        0.7382733641762232,
                        0.537435993300495,
                    ]
                ),
            )
        assert motif[2:9].consensus == "TGATCGA"
        motif = record[1]
        assert motif.name == "IFXA"
        assert record["IFXA"] == motif
        assert motif.num_occurrences == 14
        assert motif.length == 18
        assert motif.background["A"] == pytest.approx(0.30269730269730266, abs=5e-8)
        assert motif.background["C"] == pytest.approx(0.1828171828171828, abs=5e-8)
        assert motif.background["G"] == pytest.approx(0.20879120879120877, abs=5e-8)
        assert motif.background["T"] == pytest.approx(0.30569430569430567, abs=5e-8)
        assert motif.evalue == pytest.approx(3.2e-35, abs=5e-37)
        assert motif.alphabet == "ACGT"
        assert motif.alignment is None
        assert motif.consensus == "TACTGTATATATATCCAG"
        assert motif.degenerate_consensus == "TACTGTATATAHAWMCAG"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.9632889858595118,
                        1.02677956765017,
                        2.451526420551951,
                        1.7098384161433415,
                        2.2598671267551107,
                        1.7098384161433415,
                        1.02677956765017,
                        1.391583804103081,
                        1.02677956765017,
                        1.1201961888781142,
                        0.27822438781180836,
                        0.36915366971717867,
                        1.7240522753630425,
                        0.3802185945622609,
                        0.790937683007783,
                        2.451526420551951,
                        1.7240522753630425,
                        1.3924085743645374,
                    ]
                ),
            )
        assert motif[2:9].consensus == "CTGTATA"
        with open("motifs/minimal_test.meme") as stream:
            record = motifs.parse(stream, "minimal")
        motif = record[2]
        assert motif.name == "IFXA_no_nsites_no_evalue"
        assert record["IFXA_no_nsites_no_evalue"] == motif
        assert motif.num_occurrences == 20
        assert motif.length == 18
        assert motif.background["A"] == pytest.approx(0.30269730269730266, abs=5e-8)
        assert motif.background["C"] == pytest.approx(0.1828171828171828, abs=5e-8)
        assert motif.background["G"] == pytest.approx(0.20879120879120877, abs=5e-8)
        assert motif.background["T"] == pytest.approx(0.30569430569430567, abs=5e-8)
        assert motif.evalue == pytest.approx(0.0, abs=5e-37)
        assert motif.alphabet == "ACGT"
        assert motif.alignment is None
        assert motif.consensus == "TACTGTATATATATCCAG"
        assert motif.degenerate_consensus == "TACTGTATATAHAWMCAG"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.99075309,
                        1.16078104,
                        2.45152642,
                        1.70983842,
                        2.25986713,
                        1.70983842,
                        1.16078104,
                        1.46052586,
                        1.16078104,
                        1.10213019,
                        0.29911041,
                        0.36915367,
                        1.72405228,
                        0.37696488,
                        0.85258086,
                        2.45152642,
                        1.72405228,
                        1.42793329,
                    ]
                ),
            )
        assert motif[2:9].consensus == "CTGTATA"

    def test_meme_parser_rna(self):
        """Test if Bio.motifs can parse MEME output files using RNA."""
        with open("motifs/minimal_test_rna.meme") as stream:
            record = motifs.parse(stream, "minimal")
        assert record.version == "4"
        assert record.alphabet == "ACGU"
        assert len(record.sequences) == 0
        assert record.command == ""
        assert len(record) == 3
        motif = record[0]
        assert motif.name == "KRP_fake_RNA"
        assert record["KRP_fake_RNA"] == motif
        assert motif.num_occurrences == 17
        assert motif.length == 19
        assert motif.background["A"] == pytest.approx(0.30269730269730266, abs=5e-8)
        assert motif.background["C"] == pytest.approx(0.1828171828171828, abs=5e-8)
        assert motif.background["G"] == pytest.approx(0.20879120879120877, abs=5e-8)
        assert motif.background["U"] == pytest.approx(0.30569430569430567, abs=5e-8)
        assert motif.evalue == pytest.approx(4.1e-09, abs=5e-11)
        assert motif.alphabet == "ACGU"
        assert motif.alignment is None
        assert motif.consensus == "UGUGAUCGAGGUCACACUU"
        assert motif.degenerate_consensus == "UGUGANNNWGNUCACAYWW"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        1.1684297174927525,
                        0.9432809925744818,
                        1.4307101633876265,
                        1.1549413780465179,
                        0.9308256303218774,
                        0.009164393966550805,
                        0.20124190687894253,
                        0.17618542656995528,
                        0.36777933103380855,
                        0.6635834532368525,
                        0.07729943368061855,
                        0.9838293592717438,
                        1.72489868427398,
                        0.8397561713453014,
                        1.72489868427398,
                        0.8455332015343343,
                        0.3106481207768122,
                        0.7382733641762232,
                        0.537435993300495,
                    ]
                ),
            )
        assert motif[2:9].consensus == "UGAUCGA"
        motif = record[1]
        assert motif.name == "IFXA_fake_RNA"
        assert record["IFXA_fake_RNA"] == motif
        assert motif.num_occurrences == 14
        assert motif.length == 18
        assert motif.background["A"] == pytest.approx(0.30269730269730266, abs=5e-8)
        assert motif.background["C"] == pytest.approx(0.1828171828171828, abs=5e-8)
        assert motif.background["G"] == pytest.approx(0.20879120879120877, abs=5e-8)
        assert motif.background["U"] == pytest.approx(0.30569430569430567, abs=5e-8)
        assert motif.evalue == pytest.approx(3.2e-35, abs=5e-37)
        assert motif.alphabet == "ACGU"
        assert motif.alignment is None
        assert motif.consensus == "UACUGUAUAUAUAUCCAG"
        assert motif.degenerate_consensus == "UACUGUAUAUAHAWMCAG"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.9632889858595118,
                        1.02677956765017,
                        2.451526420551951,
                        1.7098384161433415,
                        2.2598671267551107,
                        1.7098384161433415,
                        1.02677956765017,
                        1.391583804103081,
                        1.02677956765017,
                        1.1201961888781142,
                        0.27822438781180836,
                        0.36915366971717867,
                        1.7240522753630425,
                        0.3802185945622609,
                        0.790937683007783,
                        2.451526420551951,
                        1.7240522753630425,
                        1.3924085743645374,
                    ]
                ),
            )
        assert motif[2:9].consensus == "CUGUAUA"

        motif = record[2]
        assert motif.name == "IFXA_no_nsites_no_evalue_fake_RNA"
        assert record["IFXA_no_nsites_no_evalue_fake_RNA"] == motif
        assert motif.num_occurrences == 20
        assert motif.length == 18
        assert motif.background["A"] == pytest.approx(0.30269730269730266, abs=5e-8)
        assert motif.background["C"] == pytest.approx(0.1828171828171828, abs=5e-8)
        assert motif.background["G"] == pytest.approx(0.20879120879120877, abs=5e-8)
        assert motif.background["U"] == pytest.approx(0.30569430569430567, abs=5e-8)
        assert motif.evalue == pytest.approx(0.0, abs=5e-37)
        assert motif.alphabet == "ACGU"
        assert motif.alignment is None
        assert motif.consensus == "UACUGUAUAUAUAUCCAG"
        assert motif.degenerate_consensus == "UACUGUAUAUAHAWMCAG"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.99075309,
                        1.16078104,
                        2.45152642,
                        1.70983842,
                        2.25986713,
                        1.70983842,
                        1.16078104,
                        1.46052586,
                        1.16078104,
                        1.10213019,
                        0.29911041,
                        0.36915367,
                        1.72405228,
                        0.37696488,
                        0.85258086,
                        2.45152642,
                        1.72405228,
                        1.42793329,
                    ]
                ),
            )
        assert motif[2:9].consensus == "CUGUAUA"


class TestMAST(unittest.TestCase):
    """MAST format tests."""

    def test_mast_parser_1(self):
        """Parse motifs/mast.crp0.de.oops.txt.xml file."""
        with open("motifs/mast.crp0.de.oops.txt.xml") as stream:
            record = motifs.parse(stream, "MAST")
        assert record.version == "5.0.1"
        assert record.database == "common/crp0.s"
        assert record.alphabet == "DNA"
        assert len(record) == 2
        assert len(record.sequences) == 18
        assert record.sequences[0] == "lac"
        assert record.sequences[1] == "bglr1"
        assert record.sequences[2] == "tdc"
        assert record.sequences[3] == "deop2"
        assert record.sequences[4] == "pbr322"
        assert record.sequences[5] == "malk"
        assert record.sequences[6] == "tnaa"
        assert record.sequences[7] == "male"
        assert record.sequences[8] == "ara"
        assert record.sequences[9] == "cya"
        assert record.sequences[10] == "ompa"
        assert record.sequences[11] == "ilv"
        assert record.sequences[12] == "gale"
        assert record.sequences[13] == "malt"
        assert record.sequences[14] == "crp"
        assert record.sequences[15] == "ce1cg"
        assert record.sequences[16] == "trn9cat"
        assert record.sequences[17] == "uxu1"
        assert record.diagrams["lac"] == "[+1]-2-[-2]-79"
        assert record.diagrams["bglr1"] == "79-[+2]-14"
        assert record.diagrams["tdc"] == "30-[+1]-39-[+2]-12"
        assert record.diagrams["deop2"] == "19-[+1]-74"
        assert record.diagrams["pbr322"] == "58-[-2]-35"
        assert record.diagrams["malk"] == "32-[+2]-61"
        assert record.diagrams["tnaa"] == "105"
        assert record.diagrams["male"] == "105"
        assert record.diagrams["ara"] == "105"
        assert record.diagrams["cya"] == "105"
        assert record.diagrams["ompa"] == "105"
        assert record.diagrams["ilv"] == "105"
        assert record.diagrams["gale"] == "105"
        assert record.diagrams["malt"] == "105"
        assert record.diagrams["crp"] == "105"
        assert record.diagrams["ce1cg"] == "105"
        assert record.diagrams["trn9cat"] == "105"
        assert record.diagrams["uxu1"] == "105"
        motif = record[0]
        assert record["1"] is motif
        assert motif.name == "1"
        assert motif.alphabet == "DNA"
        assert motif.length == 12
        assert motif[1:-2].length == 9
        motif = record[1]
        assert record["2"] is motif
        assert motif.name == "2"
        assert motif.alphabet == "DNA"
        assert motif.length == 12
        assert motif[1:40:3].length == 4

    def test_mast_parser_2(self):
        """Parse motifs/mast.adh.de.oops.html.xml file."""
        with open("motifs/mast.adh.de.oops.html.xml") as stream:
            record = motifs.parse(stream, "MAST")
        assert record.version == "5.0.1"
        assert record.database == "common/adh.s"
        assert record.alphabet == "Protein"
        assert len(record.sequences) == 33
        assert record.sequences[0] == "ENTA_ECOLI"
        assert record.sequences[1] == "DHII_HUMAN"
        assert record.sequences[2] == "YINL_LISMO"
        assert record.sequences[3] == "FIXR_BRAJA"
        assert record.sequences[4] == "HDHA_ECOLI"
        assert record.sequences[5] == "BUDC_KLETE"
        assert record.sequences[6] == "AP27_MOUSE"
        assert record.sequences[7] == "FVT1_HUMAN"
        assert record.sequences[8] == "YRTP_BACSU"
        assert record.sequences[9] == "DHMA_FLAS1"
        assert record.sequences[10] == "HDE_CANTR"
        assert record.sequences[11] == "3BHD_COMTE"
        assert record.sequences[12] == "BDH_HUMAN"
        assert record.sequences[13] == "2BHD_STREX"
        assert record.sequences[14] == "BA72_EUBSP"
        assert record.sequences[15] == "RIDH_KLEAE"
        assert record.sequences[16] == "DHGB_BACME"
        assert record.sequences[17] == "PCR_PEA"
        assert record.sequences[18] == "RFBB_NEIGO"
        assert record.sequences[19] == "BPHB_PSEPS"
        assert record.sequences[20] == "DHB2_HUMAN"
        assert record.sequences[21] == "NODG_RHIME"
        assert record.sequences[22] == "MAS1_AGRRA"
        assert record.sequences[23] == "DHCA_HUMAN"
        assert record.sequences[24] == "DHES_HUMAN"
        assert record.sequences[25] == "DHB3_HUMAN"
        assert record.sequences[26] == "HMTR_LEIMA"
        assert record.sequences[27] == "ADH_DROME"
        assert record.sequences[28] == "YURA_MYXXA"
        assert record.sequences[29] == "LIGD_PSEPA"
        assert record.sequences[30] == "FABI_ECOLI"
        assert record.sequences[31] == "GUTD_ECOLI"
        assert record.sequences[32] == "CSGA_MYXXA"
        assert record.diagrams["ENTA_ECOLI"] == "[1]-[2]-224"
        assert record.diagrams["DHII_HUMAN"] == "29-[1]-[2]-239"
        assert record.diagrams["YINL_LISMO"] == "[1]-[2]-224"
        assert record.diagrams["FIXR_BRAJA"] == "43-[2]-149-[1]-62"
        assert record.diagrams["HDHA_ECOLI"] == "6-[1]-[2]-144-[1]-69"
        assert record.diagrams["BUDC_KLETE"] == "9-[2]-53-[1]-81-[1]-62"
        assert record.diagrams["AP27_MOUSE"] == "2-[1]-[2]-138-[1]-68"
        assert record.diagrams["FVT1_HUMAN"] == "39-[2]-150-[1]-119"
        assert record.diagrams["YRTP_BACSU"] == "1-[1]-[2]-145-[1]-56"
        assert record.diagrams["DHMA_FLAS1"] == "9-[1]-[2]-147-[1]-78"
        assert record.diagrams["HDE_CANTR"] == "3-[1]-[2]-290-[1]-[2]-565"
        assert record.diagrams["3BHD_COMTE"] == "1-[1]-[2]-50-[1]-166"
        assert record.diagrams["BDH_HUMAN"] == "50-[1]-[2]-269"
        assert record.diagrams["2BHD_STREX"] == "1-[1]-[2]-142-[1]-76"
        assert record.diagrams["BA72_EUBSP"] == "1-[1]-[2]-125-[2]-10-[1]-65"
        assert record.diagrams["RIDH_KLEAE"] == "9-[1]-[2]-216"
        assert record.diagrams["DHGB_BACME"] == "2-[1]-[2]-149-[1]-75"
        assert record.diagrams["PCR_PEA"] == "81-[1]-[2]-108-[1]-174"
        assert record.diagrams["RFBB_NEIGO"] == "1-[1]-[2]-321"
        assert record.diagrams["BPHB_PSEPS"] == "[1]-[2]-251"
        assert record.diagrams["DHB2_HUMAN"] == "77-[1]-[2]-286"
        assert record.diagrams["NODG_RHIME"] == "1-[1]-[2]-142-[1]-66"
        assert record.diagrams["MAS1_AGRRA"] == "252-[2]-36-[1]-164"
        assert record.diagrams["DHCA_HUMAN"] == "11-[2]-54-[1]-101-[1]-74"
        assert record.diagrams["DHES_HUMAN"] == "9-[2]-108-[1]-186"
        assert record.diagrams["DHB3_HUMAN"] == "55-[2]-146-[1]-85"
        assert record.diagrams["HMTR_LEIMA"] == "24-[2]-172-[1]-67"
        assert record.diagrams["ADH_DROME"] == "13-[2]-217-[1]-1"
        assert record.diagrams["YURA_MYXXA"] == "94-[2]-69-[1]-71"
        assert record.diagrams["LIGD_PSEPA"] == "1-[1]-[2]-280"
        assert record.diagrams["FABI_ECOLI"] == "1-[1]-161-[1]-76"
        assert record.diagrams["GUTD_ECOLI"] == "147-[2]-10-[1]-78"
        assert record.diagrams["CSGA_MYXXA"] == "12-[1]-53-[2]-77"
        assert len(record) == 2
        motif = record[0]
        assert record["1"] is motif
        assert motif.alphabet == "Protein"
        assert motif.length == 12
        assert motif.name == "1"
        assert motif[1:-2].length == 9
        motif = record[1]
        assert record["2"] is motif
        assert motif.alphabet == "Protein"
        assert motif.length == 12
        assert motif.name == "2"
        assert motif[-20:-2].length == 10

    def test_mast_parser_3(self):
        """Parse motifs/mast.Klf1-200.cd.oops.xml.xml file."""
        with open("motifs/mast.Klf1-200.cd.oops.xml.xml") as stream:
            record = motifs.parse(stream, "MAST")
        assert record.version == "5.0.1"
        assert record.database == "common/Klf1-200.fa"
        assert record.alphabet == "DNA"
        assert len(record.sequences) == 113
        assert record.sequences[0] == "chr3:104843905-104844405"
        assert record.sequences[1] == "chr12:114390660-114391160"
        assert record.sequences[2] == "chr12:27135944-27136444"
        assert record.sequences[3] == "chr10:59256089-59256589"
        assert record.sequences[4] == "chr4:135733850-135734350"
        assert record.sequences[5] == "chr1:137838164-137838664"
        assert record.sequences[6] == "chr17:47735006-47735506"
        assert record.sequences[7] == "chr6:72223026-72223526"
        assert record.sequences[8] == "chr13:3866266-3866766"
        assert record.sequences[9] == "chr1:133343883-133344383"
        assert record.sequences[10] == "chr11:117187372-117187872"
        assert record.sequences[11] == "chr13:76003199-76003699"
        assert record.sequences[12] == "chr5:65202593-65203093"
        assert record.sequences[13] == "chr14:79702844-79703344"
        assert record.sequences[14] == "chr12:112796794-112797294"
        assert record.sequences[15] == "chr13:112863645-112864145"
        assert record.sequences[16] == "chr7:111007530-111008030"
        assert record.sequences[17] == "chr1:43307690-43308190"
        assert record.sequences[18] == "chr14:47973722-47974222"
        assert record.sequences[19] == "chr9:120025371-120025871"
        assert record.sequences[20] == "chr7:105490727-105491227"
        assert record.sequences[21] == "chr5:37127175-37127675"
        assert record.sequences[22] == "chr5:45951565-45952065"
        assert record.sequences[23] == "chr7:91033422-91033922"
        assert record.sequences[24] == "chr4:154285745-154286245"
        assert record.sequences[25] == "chr13:100518008-100518508"
        assert record.sequences[26] == "chr1:36977019-36977519"
        assert record.sequences[27] == "chr7:151917814-151918314"
        assert record.sequences[28] == "chr7:110976195-110976695"
        assert record.sequences[29] == "chr15:58719281-58719781"
        assert record.sequences[30] == "chr11:57590460-57590960"
        assert record.sequences[31] == "chr8:83025150-83025650"
        assert record.sequences[32] == "chr13:54345922-54346422"
        assert record.sequences[33] == "chr12:82044358-82044858"
        assert record.sequences[34] == "chr11:105013714-105014214"
        assert record.sequences[35] == "chr10:93585404-93585904"
        assert record.sequences[36] == "chr7:19832207-19832707"
        assert record.sequences[37] == "chr8:97323995-97324495"
        assert record.sequences[38] == "chr10:126642277-126642777"
        assert record.sequences[39] == "chr1:156887119-156887619"
        assert record.sequences[40] == "chr15:81700367-81700867"
        assert record.sequences[41] == "chr6:121187425-121187925"
        assert record.sequences[42] == "chr4:43977111-43977611"
        assert record.sequences[43] == "chr11:102236405-102236905"
        assert record.sequences[44] == "chr17:5112057-5112557"
        assert record.sequences[45] == "chr10:110604369-110604869"
        assert record.sequences[46] == "chr1:169314208-169314708"
        assert record.sequences[47] == "chr9:57618594-57619094"
        assert record.sequences[48] == "chr10:128184604-128185104"
        assert record.sequences[49] == "chr4:109112541-109113041"
        assert record.sequences[50] == "chr3:97461668-97462168"
        assert record.sequences[51] == "chr9:102674395-102674895"
        assert record.sequences[52] == "chr17:24289205-24289705"
        assert record.sequences[53] == "chr17:28960252-28960752"
        assert record.sequences[54] == "chr2:73323093-73323593"
        assert record.sequences[55] == "chr11:32150818-32151318"
        assert record.sequences[56] == "chr7:103853792-103854292"
        assert record.sequences[57] == "chr16:49839621-49840121"
        assert record.sequences[58] == "chr6:135115628-135116128"
        assert record.sequences[59] == "chr3:88305500-88306000"
        assert record.sequences[60] == "chr18:57137388-57137888"
        assert record.sequences[61] == "chr5:97380648-97381148"
        assert record.sequences[62] == "chr15:91082416-91082916"
        assert record.sequences[63] == "chr14:61272713-61273213"
        assert record.sequences[64] == "chr5:33616214-33616714"
        assert record.sequences[65] == "chr18:23982470-23982970"
        assert record.sequences[66] == "chr9:24715045-24715545"
        assert record.sequences[67] == "chr10:116195445-116195945"
        assert record.sequences[68] == "chr11:77795184-77795684"
        assert record.sequences[69] == "chr16:32508975-32509475"
        assert record.sequences[70] == "chr18:80416880-80417380"
        assert record.sequences[71] == "chr10:57252236-57252736"
        assert record.sequences[72] == "chr5:34915767-34916267"
        assert record.sequences[73] == "chr9:98389943-98390443"
        assert record.sequences[74] == "chr19:5845899-5846399"
        assert record.sequences[75] == "chr3:151777796-151778296"
        assert record.sequences[76] == "chr4:76585120-76585620"
        assert record.sequences[77] == "chr7:104332488-104332988"
        assert record.sequences[78] == "chr5:138127197-138127697"
        assert record.sequences[79] == "chr11:60988820-60989320"
        assert record.sequences[80] == "chr8:19984030-19984530"
        assert record.sequences[81] == "chr11:31712262-31712762"
        assert record.sequences[82] == "chr15:41338514-41339014"
        assert record.sequences[83] == "chr9:21362671-21363171"
        assert record.sequences[84] == "chr18:58822702-58823202"
        assert record.sequences[85] == "chr1:173447614-173448114"
        assert record.sequences[86] == "chr6:81915769-81916269"
        assert record.sequences[87] == "chr1:169322898-169323398"
        assert record.sequences[88] == "chr12:70860461-70860961"
        assert record.sequences[89] == "chr9:59598186-59598686"
        assert record.sequences[90] == "chr3:19550495-19550995"
        assert record.sequences[91] == "chr7:36132953-36133453"
        assert record.sequences[92] == "chr7:38970375-38970875"
        assert record.sequences[93] == "chr15:78243390-78243890"
        assert record.sequences[94] == "chr7:87847381-87847881"
        assert record.sequences[95] == "chr1:33631214-33631714"
        assert record.sequences[96] == "chr4:135407873-135408373"
        assert record.sequences[97] == "chr7:101244829-101245329"
        assert record.sequences[98] == "chr10:60612190-60612690"
        assert record.sequences[99] == "chr19:56465963-56466463"
        assert record.sequences[100] == "chr4:41334759-41335259"
        assert record.sequences[101] == "chr8:92969521-92970021"
        assert record.sequences[102] == "chr6:145703215-145703715"
        assert record.sequences[103] == "chr13:57679178-57679678"
        assert record.sequences[104] == "chr19:45121628-45122128"
        assert record.sequences[105] == "chr15:79757891-79758391"
        assert record.sequences[106] == "chr1:134264178-134264678"
        assert record.sequences[107] == "chr13:81067500-81068000"
        assert record.sequences[108] == "chr11:69714224-69714724"
        assert record.sequences[109] == "chr2:103728071-103728571"
        assert record.sequences[110] == "chr5:105994747-105995247"
        assert record.sequences[111] == "chr17:84209565-84210065"
        assert record.sequences[112] == "chr7:16507689-16508189"
        assert record.diagrams["chr3:104843905-104844405"] == "115-[-1]-209-[-2]-126"
        assert record.diagrams["chr12:114390660-114391160"] == "3-[+2]-[+2]-3-[+1]-173-[+1]-3-[-2]-188"
        assert record.diagrams["chr12:27135944-27136444"] == "275-[-1]-89-[+2]-4-[+2]-52"
        assert record.diagrams["chr10:59256089-59256589"] == "247-[+2]-17-[-1]-186"
        assert record.diagrams["chr4:135733850-135734350"] == "183-[-1]-263-[+2]-4"
        assert record.diagrams["chr1:137838164-137838664"] == "192-[-2]-1-[+1]-44-[-1]-193"
        assert record.diagrams["chr17:47735006-47735506"] == "203-[+2]-15-[+1]-97-[-1]-115"
        assert record.diagrams["chr6:72223026-72223526"] == "52-[-2]-7-[+2]-162-[-1]-42-[-1]-137"
        assert record.diagrams["chr13:3866266-3866766"] == "241-[+1]-2-[-1]-217"
        assert record.diagrams["chr1:133343883-133344383"] == "190-[+2]-15-[+1]-245"
        assert record.diagrams["chr11:117187372-117187872"] == "242-[+1]-46-[-2]-71-[+1]-71"
        assert record.diagrams["chr13:76003199-76003699"] == "230-[+2]-15-[+2]-60-[-1]-115"
        assert record.diagrams["chr5:65202593-65203093"] == "24-[-2]-36-[+2]-193-[-1]-11-[+1]-10-[+1]-106"
        assert record.diagrams["chr14:79702844-79703344"] == "247-[-1]-46-[-2]-157"
        assert record.diagrams["chr12:112796794-112797294"] == "232-[+1]-41-[+1]-187"
        assert record.diagrams["chr13:112863645-112864145"] == "228-[+1]-20-[-1]-212"
        assert record.diagrams["chr7:111007530-111008030"] == "217-[+1]-83-[+2]-150"
        assert record.diagrams["chr1:43307690-43308190"] == "164-[-2]-52-[-2]-224"
        assert record.diagrams["chr14:47973722-47974222"] == "21-[+1]-181-[+1]-20-[-2]-208"
        assert record.diagrams["chr9:120025371-120025871"] == "110-[-2]-58-[+1]-282"
        assert record.diagrams["chr7:105490727-105491227"] == "100-[-2]-111-[-1]-239"
        assert record.diagrams["chr5:37127175-37127675"] == "234-[-2]-24-[+1]-192"
        assert record.diagrams["chr5:45951565-45952065"] == "261-[-1]-219"
        assert record.diagrams["chr7:91033422-91033922"] == "465-[-1]-15"
        assert record.diagrams["chr4:154285745-154286245"] == "235-[+1]-20-[-2]-195"
        assert record.diagrams["chr13:100518008-100518508"] == "226-[-2]-18-[-1]-206"
        assert record.diagrams["chr1:36977019-36977519"] == "88-[+1]-187-[+2]-60-[-1]-95"
        assert record.diagrams["chr7:151917814-151918314"] == "219-[+1]-80-[+2]-151"
        assert record.diagrams["chr7:110976195-110976695"] == "287-[+2]-12-[+1]-151"
        assert record.diagrams["chr15:58719281-58719781"] == "212-[-2]-258"
        assert record.diagrams["chr11:57590460-57590960"] == "56-[-1]-271-[-1]-75-[+2]-28"
        assert record.diagrams["chr8:83025150-83025650"] == "219-[+1]-87-[+2]-144"
        assert record.diagrams["chr13:54345922-54346422"] == "283-[-2]-161-[+1]-6"
        assert record.diagrams["chr12:82044358-82044858"] == "50-[+2]-160-[+1]-39-[+2]-171"
        assert record.diagrams["chr11:105013714-105014214"] == "115-[-2]-160-[+1]-26-[-1]-129"
        assert record.diagrams["chr10:93585404-93585904"] == "141-[+2]-48-[+1]-261"
        assert record.diagrams["chr7:19832207-19832707"] == "229-[-1]-251"
        assert record.diagrams["chr8:97323995-97324495"] == "177-[-1]-40-[-2]-139-[+1]-74"
        assert record.diagrams["chr10:126642277-126642777"] == "252-[-1]-92-[-2]-106"
        assert record.diagrams["chr1:156887119-156887619"] == "189-[-2]-78-[-1]-183"
        assert record.diagrams["chr15:81700367-81700867"] == "109-[-1]-99-[-1]-252"
        assert record.diagrams["chr6:121187425-121187925"] == "29-[+2]-313-[-1]-108"
        assert record.diagrams["chr4:43977111-43977611"] == "60-[+1]-148-[+1]-252"
        assert record.diagrams["chr11:102236405-102236905"] == "10-[+2]-145-[-1]-3-[-1]-6-[+2]-60-[+1]-156"
        assert record.diagrams["chr17:5112057-5112557"] == "249-[+1]-231"
        assert record.diagrams["chr10:110604369-110604869"] == "232-[+1]-248"
        assert record.diagrams["chr1:169314208-169314708"] == "192-[-1]-[-1]-11-[-2]-227"
        assert record.diagrams["chr9:57618594-57619094"] == "125-[+2]-151-[-1]-4-[-1]-150"
        assert record.diagrams["chr10:128184604-128185104"] == "30-[-2]-128-[+1]-292"
        assert record.diagrams["chr4:109112541-109113041"] == "21-[-1]-13-[+1]-94-[+2]-302"
        assert record.diagrams["chr3:97461668-97462168"] == "18-[+2]-256-[-1]-81-[+1]-21-[+1]-34"
        assert record.diagrams["chr9:102674395-102674895"] == "372-[+2]-98"
        assert record.diagrams["chr17:24289205-24289705"] == "262-[-1]-218"
        assert record.diagrams["chr17:28960252-28960752"] == "221-[+1]-81-[+1]-158"
        assert record.diagrams["chr2:73323093-73323593"] == "49-[-2]-421"
        assert record.diagrams["chr11:32150818-32151318"] == "151-[-1]-27-[-1]-118-[-2]-134"
        assert record.diagrams["chr7:103853792-103854292"] == "212-[-2]-42-[+1]-196"
        assert record.diagrams["chr16:49839621-49840121"] == "192-[+2]-47-[-1]-17-[+2]-164"
        assert record.diagrams["chr6:135115628-135116128"] == "231-[-1]-249"
        assert record.diagrams["chr3:88305500-88306000"] == "229-[+1]-251"
        assert record.diagrams["chr18:57137388-57137888"] == "296-[+2]-174"
        assert record.diagrams["chr5:97380648-97381148"] == "188-[-2]-282"
        assert record.diagrams["chr15:91082416-91082916"] == "239-[-1]-104-[-1]-73-[+2]-14"
        assert record.diagrams["chr14:61272713-61273213"] == "216-[+2]-104-[+1]-130"
        assert record.diagrams["chr5:33616214-33616714"] == "247-[-1]-233"
        assert record.diagrams["chr18:23982470-23982970"] == "285-[-1]-195"
        assert record.diagrams["chr9:24715045-24715545"] == "214-[-1]-153-[+1]-93"
        assert record.diagrams["chr10:116195445-116195945"] == "400-[+2]-70"
        assert record.diagrams["chr11:77795184-77795684"] == "247-[+1]-42-[-2]-67-[-2]-64"
        assert record.diagrams["chr16:32508975-32509475"] == "213-[+2]-29-[-1]-208"
        assert record.diagrams["chr18:80416880-80417380"] == "239-[-1]-241"
        assert record.diagrams["chr10:57252236-57252736"] == "155-[+1]-158-[+2]-137"
        assert record.diagrams["chr5:34915767-34916267"] == "179-[+2]-29-[-1]-242"
        assert record.diagrams["chr9:98389943-98390443"] == "252-[-1]-228"
        assert record.diagrams["chr19:5845899-5846399"] == "136-[+1]-193-[+1]-131"
        assert record.diagrams["chr3:151777796-151778296"] == "30-[-2]-58-[-1]-362"
        assert record.diagrams["chr4:76585120-76585620"] == "329-[+2]-141"
        assert record.diagrams["chr7:104332488-104332988"] == "164-[+2]-23-[-1]-222-[+1]-21"
        assert record.diagrams["chr5:138127197-138127697"] == "238-[+1]-242"
        assert record.diagrams["chr11:60988820-60989320"] == "115-[+1]-68-[+1]-47-[+1]-210"
        assert record.diagrams["chr8:19984030-19984530"] == "103-[-1]-81-[+2]-266"
        assert record.diagrams["chr11:31712262-31712762"] == "118-[+2]-53-[+2]-269"
        assert record.diagrams["chr15:41338514-41339014"] == "173-[+2]-75-[+2]-192"
        assert record.diagrams["chr9:21362671-21363171"] == "105-[+1]-131-[+1]-224"
        assert record.diagrams["chr18:58822702-58823202"] == "467-[-2]-3"
        assert record.diagrams["chr1:173447614-173448114"] == "369-[-1]-111"
        assert record.diagrams["chr6:81915769-81916269"] == "197-[+1]-283"
        assert record.diagrams["chr1:169322898-169323398"] == "253-[-1]-227"
        assert record.diagrams["chr12:70860461-70860961"] == "197-[+2]-22-[-1]-231"
        assert record.diagrams["chr9:59598186-59598686"] == "163-[-2]-10-[-1]-277"
        assert record.diagrams["chr3:19550495-19550995"] == "452-[-2]-18"
        assert record.diagrams["chr7:36132953-36133453"] == "157-[-1]-323"
        assert record.diagrams["chr7:38970375-38970875"] == "49-[+1]-114-[+1]-297"
        assert record.diagrams["chr15:78243390-78243890"] == "234-[+1]-246"
        assert record.diagrams["chr7:87847381-87847881"] == "99-[+2]-2-[-1]-230-[-1]-99"
        assert record.diagrams["chr1:33631214-33631714"] == "358-[-1]-122"
        assert record.diagrams["chr4:135407873-135408373"] == "116-[-1]-64-[+2]-270"
        assert record.diagrams["chr7:101244829-101245329"] == "311-[-2]-159"
        assert record.diagrams["chr10:60612190-60612690"] == "215-[+1]-265"
        assert record.diagrams["chr19:56465963-56466463"] == "306-[+1]-36-[+1]-18-[+1]-80"
        assert record.diagrams["chr4:41334759-41335259"] == "204-[+1]-276"
        assert record.diagrams["chr8:92969521-92970021"] == "453-[+2]-17"
        assert record.diagrams["chr6:145703215-145703715"] == "154-[-2]-58-[+2]-228"
        assert record.diagrams["chr13:57679178-57679678"] == "217-[-1]-263"
        assert record.diagrams["chr19:45121628-45122128"] == "35-[-2]-435"
        assert record.diagrams["chr15:79757891-79758391"] == "310-[+1]-170"
        assert record.diagrams["chr1:134264178-134264678"] == "23-[+2]-447"
        assert record.diagrams["chr13:81067500-81068000"] == "252-[+1]-228"
        assert record.diagrams["chr11:69714224-69714724"] == "145-[+2]-325"
        assert record.diagrams["chr2:103728071-103728571"] == "369-[+1]-111"
        assert record.diagrams["chr5:105994747-105995247"] == "93-[+2]-153-[-2]-194"
        assert record.diagrams["chr17:84209565-84210065"] == "64-[-2]-406"
        assert record.diagrams["chr7:16507689-16508189"] == "231-[+2]-239"
        assert len(record) == 2
        motif = record[0]
        assert record["1"] is motif
        assert motif.alphabet == "DNA"
        assert motif.length == 20
        assert motif.name == "1"
        assert motif[1:-2].length == 17
        motif = record[1]
        assert record["2"] is motif
        assert motif.alphabet == "DNA"
        assert motif.length == 30
        assert motif.name == "2"
        assert motif[10:20].length == 10


class TestTransfac(unittest.TestCase):
    """Transfac format tests."""

    def test_transfac_parser(self):
        """Parse motifs/transfac.dat file."""
        with open("motifs/transfac.dat") as stream:
            record = motifs.parse(stream, "TRANSFAC")
        motif = record[0]
        assert motif["ID"] == "motif1"
        assert len(motif.counts) == 4
        assert motif.counts.length == 12
        assert motif.counts["A", 0] == 1
        assert motif.counts["A", 1] == 2
        assert motif.counts["A", 2] == 3
        assert motif.counts["A", 3] == 0
        assert motif.counts["A", 4] == 5
        assert motif.counts["A", 5] == 0
        assert motif.counts["A", 6] == 0
        assert motif.counts["A", 7] == 0
        assert motif.counts["A", 8] == 0
        assert motif.counts["A", 9] == 0
        assert motif.counts["A", 10] == 0
        assert motif.counts["A", 11] == 1
        assert motif.counts["C", 0] == 2
        assert motif.counts["C", 1] == 1
        assert motif.counts["C", 2] == 0
        assert motif.counts["C", 3] == 5
        assert motif.counts["C", 4] == 0
        assert motif.counts["C", 5] == 0
        assert motif.counts["C", 6] == 1
        assert motif.counts["C", 7] == 0
        assert motif.counts["C", 8] == 0
        assert motif.counts["C", 9] == 1
        assert motif.counts["C", 10] == 2
        assert motif.counts["C", 11] == 0
        assert motif.counts["G", 0] == 2
        assert motif.counts["G", 1] == 2
        assert motif.counts["G", 2] == 1
        assert motif.counts["G", 3] == 0
        assert motif.counts["G", 4] == 0
        assert motif.counts["G", 5] == 4
        assert motif.counts["G", 6] == 4
        assert motif.counts["G", 7] == 0
        assert motif.counts["G", 8] == 5
        assert motif.counts["G", 9] == 2
        assert motif.counts["G", 10] == 0
        assert motif.counts["G", 11] == 3
        assert motif.counts["T", 0] == 0
        assert motif.counts["T", 1] == 0
        assert motif.counts["T", 2] == 1
        assert motif.counts["T", 3] == 0
        assert motif.counts["T", 4] == 0
        assert motif.counts["T", 5] == 1
        assert motif.counts["T", 6] == 0
        assert motif.counts["T", 7] == 5
        assert motif.counts["T", 8] == 0
        assert motif.counts["T", 9] == 2
        assert motif.counts["T", 10] == 3
        assert motif.counts["T", 11] == 1
        assert motif.degenerate_consensus == "SRACAGGTGKYG"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.4780719051126377,
                        0.4780719051126377,
                        0.6290494055453314,
                        2.0,
                        2.0,
                        1.278071905112638,
                        1.278071905112638,
                        2.0,
                        2.0,
                        0.4780719051126377,
                        1.0290494055453312,
                        0.6290494055453314,
                    ]
                ),
            )
        assert motif[1:-2].degenerate_consensus == "RACAGGTGK"
        assert np.allclose(
                motif[1:-2].relative_entropy,
                np.array(
                    [
                        0.4780719051126377,
                        0.6290494055453314,
                        2.0,
                        2.0,
                        1.278071905112638,
                        1.278071905112638,
                        2.0,
                        2.0,
                        0.4780719051126377,
                    ]
                ),
            )
        motif = record[1]
        assert motif["ID"] == "motif2"
        assert len(motif.counts) == 4
        assert motif.counts.length == 10
        assert motif.counts["A", 0] == 2
        assert motif.counts["A", 1] == 1
        assert motif.counts["A", 2] == 0
        assert motif.counts["A", 3] == 3
        assert motif.counts["A", 4] == 0
        assert motif.counts["A", 5] == 5
        assert motif.counts["A", 6] == 0
        assert motif.counts["A", 7] == 0
        assert motif.counts["A", 8] == 0
        assert motif.counts["A", 9] == 0
        assert motif.counts["C", 0] == 1
        assert motif.counts["C", 1] == 2
        assert motif.counts["C", 2] == 5
        assert motif.counts["C", 3] == 0
        assert motif.counts["C", 4] == 0
        assert motif.counts["C", 5] == 0
        assert motif.counts["C", 6] == 1
        assert motif.counts["C", 7] == 0
        assert motif.counts["C", 8] == 0
        assert motif.counts["C", 9] == 2
        assert motif.counts["G", 0] == 2
        assert motif.counts["G", 1] == 2
        assert motif.counts["G", 2] == 0
        assert motif.counts["G", 3] == 1
        assert motif.counts["G", 4] == 4
        assert motif.counts["G", 5] == 0
        assert motif.counts["G", 6] == 4
        assert motif.counts["G", 7] == 5
        assert motif.counts["G", 8] == 0
        assert motif.counts["G", 9] == 0
        assert motif.counts["T", 0] == 0
        assert motif.counts["T", 1] == 0
        assert motif.counts["T", 2] == 0
        assert motif.counts["T", 3] == 1
        assert motif.counts["T", 4] == 1
        assert motif.counts["T", 5] == 0
        assert motif.counts["T", 6] == 0
        assert motif.counts["T", 7] == 0
        assert motif.counts["T", 8] == 5
        assert motif.counts["T", 9] == 3
        assert motif.degenerate_consensus == "RSCAGAGGTY"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.4780719051126377,
                        0.4780719051126377,
                        2.0,
                        0.6290494055453314,
                        1.278071905112638,
                        2.0,
                        1.278071905112638,
                        2.0,
                        2.0,
                        1.0290494055453312,
                    ]
                ),
            )
        assert motif[::2].degenerate_consensus == "RCGGT"
        assert np.allclose(
                motif[::2].relative_entropy,
                np.array(
                    [0.4780719051126377, 2.0, 1.278071905112638, 1.278071905112638, 2.0]
                ),
            )

    def test_permissive_transfac_parser(self):
        """Parse the TRANSFAC-like file motifs/MA0056.1.transfac."""
        # The test file MA0056.1.transfac was obtained from the JASPAR database
        # in a TRANSFAC-like format.
        # Khan, A. et al. JASPAR 2018: update of the open-access database of
        # transcription factor binding profiles and its web framework.
        # Nucleic Acids Res. 2018; 46:D260-D266,
        path = "motifs/MA0056.1.transfac"
        with open(path) as stream:
            with pytest.raises(ValueError):
                motifs.parse(stream, "TRANSFAC")
        with open(path) as stream:
            records = motifs.parse(stream, "TRANSFAC", strict=False)
        motif = records[0]
        assert sorted(motif.keys()) == ["AC", "DE", "ID"]
        assert motif["AC"] == "MA0056.1"
        assert motif["DE"] == "MA0056.1 MZF1 ; From JASPAR 2018"
        assert motif["ID"] == "MZF1"
        assert motif.counts.length == 6
        assert len(motif.counts) == 4
        assert motif.counts["A", 0] == 3.0
        assert motif.counts["A", 1] == 0.0
        assert motif.counts["A", 2] == 2.0
        assert motif.counts["A", 3] == 0.0
        assert motif.counts["A", 4] == 0.0
        assert motif.counts["A", 5] == 18.0
        assert motif.counts["C", 0] == 5.0
        assert motif.counts["C", 1] == 0.0
        assert motif.counts["C", 2] == 0.0
        assert motif.counts["C", 3] == 0.0
        assert motif.counts["C", 4] == 0.0
        assert motif.counts["C", 5] == 0.0
        assert motif.counts["G", 0] == 4.0
        assert motif.counts["G", 1] == 19.0
        assert motif.counts["G", 2] == 18.0
        assert motif.counts["G", 3] == 19.0
        assert motif.counts["G", 4] == 20.0
        assert motif.counts["G", 5] == 2.0
        assert motif.counts["T", 0] == 8.0
        assert motif.counts["T", 1] == 1.0
        assert motif.counts["T", 2] == 0.0
        assert motif.counts["T", 3] == 1.0
        assert motif.counts["T", 4] == 0.0
        assert motif.counts["T", 5] == 0.0
        assert motif.consensus == "TGGGGA"
        assert motif.degenerate_consensus == "NGGGGA"
        assert np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.09629830394265171,
                        1.7136030428840439,
                        1.5310044064107189,
                        1.7136030428840439,
                        2.0,
                        1.5310044064107189,
                    ]
                ),
            )
        assert motif[1:-3].degenerate_consensus == "GG"
        assert np.allclose(
                motif[1:-3].relative_entropy,
                np.array([1.7136030428840439, 1.5310044064107189]),
            )

    def test_TFoutput(self):
        """Ensure that we can write proper TransFac output files."""
        m = motifs.create([Seq("ATATA")])
        with tempfile.TemporaryFile("w") as stream:
            stream.write(format(m, "transfac"))


class MotifTestPWM(unittest.TestCase):
    """PWM motif tests."""

    with open("motifs/SRF.pfm") as stream:
        m = motifs.read(stream, "pfm")

    s = Seq("ACGTGTGCGTAGTGCGT")

    def test_getitem(self):
        counts = self.m.counts
        python_integers = range(13)
        numpy_integers = np.array(python_integers)
        integers = {"python": python_integers, "numpy": numpy_integers}
        for int_type in ("python", "numpy"):
            i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12 = integers[int_type]
            msg = f"using {int_type} integers as indices"
            # slice, slice
            d = counts[i1::i2, i2:i12:i3]
            assert isinstance(d, dict), msg
            assert len(d) == 2, msg
            assert len(d["C"]) == 4, msg
            assert len(d["T"]) == 4, msg
            assert d["C"][i0] == pytest.approx(45.0, abs=5e-8), msg
            assert d["C"][i1] == pytest.approx(1.0, abs=5e-8), msg
            assert d["C"][i2] == pytest.approx(0.0, abs=5e-8), msg
            assert d["C"][i3] == pytest.approx(1.0, abs=5e-8), msg
            assert d["T"][i0] == pytest.approx(0.0, abs=5e-8), msg
            assert d["T"][i1] == pytest.approx(42.0, abs=5e-8), msg
            assert d["T"][i2] == pytest.approx(3.0, abs=5e-8), msg
            assert d["T"][i3] == pytest.approx(0.0, abs=5e-8), msg
            # slice, int
            d = counts[i1::i2, i4]
            assert isinstance(d, dict), msg
            assert len(d) == 2, msg
            assert d["C"] == pytest.approx(1.0, abs=5e-8), msg
            assert d["T"] == pytest.approx(13.0, abs=5e-8), msg
            # int, slice
            t = counts[i2, i3:i12:i2]
            assert isinstance(t, tuple), msg
            assert t[i0] == pytest.approx(0.0, abs=5e-8), msg
            assert t[i1] == pytest.approx(0.0, abs=5e-8), msg
            assert t[i2] == pytest.approx(0.0, abs=5e-8), msg
            assert t[i3] == pytest.approx(0.0, abs=5e-8), msg
            assert t[i4] == pytest.approx(43.0, abs=5e-8), msg
            # int, int
            v = counts[i1, i5]
            assert v == pytest.approx(1.0, abs=5e-8), msg
            # tuple, slice
            d = counts[(i0, i3), i3:i12:i2]
            assert isinstance(d, dict), msg
            assert len(d) == 2, msg
            assert len(d["A"]) == 5, msg
            assert len(d["T"]) == 5, msg
            assert d["A"][i0] == pytest.approx(1.0, abs=5e-8), msg
            assert d["A"][i1] == pytest.approx(3.0, abs=5e-8), msg
            assert d["A"][i2] == pytest.approx(1.0, abs=5e-8), msg
            assert d["A"][i3] == pytest.approx(15.0, abs=5e-8), msg
            assert d["A"][i4] == pytest.approx(2.0, abs=5e-8), msg
            assert d["T"][i0] == pytest.approx(0.0, abs=5e-8), msg
            assert d["T"][i1] == pytest.approx(42.0, abs=5e-8), msg
            assert d["T"][i2] == pytest.approx(45.0, abs=5e-8), msg
            assert d["T"][i3] == pytest.approx(30.0, abs=5e-8), msg
            assert d["T"][i4] == pytest.approx(0.0, abs=5e-8), msg
            # tuple, int
            d = counts[(i0, i3), i5]
            assert isinstance(d, dict), msg
            assert len(d) == 2, msg
            assert d["A"] == pytest.approx(3.0, abs=5e-8), msg
            assert d["T"] == pytest.approx(42.0, abs=5e-8), msg
            # str, slice
            t = counts["C", i2:i12:i4]
            assert isinstance(t, tuple), msg
            assert t[i0] == pytest.approx(45.0, abs=5e-8), msg
            assert t[i1] == pytest.approx(0.0, abs=5e-8), msg
            assert t[i2] == pytest.approx(0.0, abs=5e-8), msg
            # str, int
            assert counts["T", i4] == pytest.approx(13.0, abs=5e-8), msg

    def test_simple(self):
        """Test if Bio.motifs PWM scoring works."""
        counts = self.m.counts
        pwm = counts.normalize(pseudocounts=0.25)
        pssm = pwm.log_odds()
        result = pssm.calculate(self.s)
        assert 6 == len(result)
        # The fast C-code in Bio/motifs/_pwm.c stores all results as 32-bit
        # floats; the slower Python code in Bio/motifs/__init__.py uses 64-bit
        # doubles. The C-code and Python code results will therefore not be
        # exactly equal. Test the first 5 decimal places only to avoid either
        # the C-code or the Python code to inadvertently fail this test.
        assert result[0] == pytest.approx(-29.18363571, abs=5e-06)
        assert result[1] == pytest.approx(-38.3365097, abs=5e-06)
        assert result[2] == pytest.approx(-29.17756271, abs=5e-06)
        assert result[3] == pytest.approx(-38.04542542, abs=5e-06)
        assert result[4] == pytest.approx(-20.3014183, abs=5e-06)
        assert result[5] == pytest.approx(-25.18009186, abs=5e-06)

    def test_with_mixed_case(self):
        """Test if Bio.motifs PWM scoring works with mixed case."""
        counts = self.m.counts
        pwm = counts.normalize(pseudocounts=0.25)
        pssm = pwm.log_odds()
        result = pssm.calculate(Seq("AcGTgTGCGtaGTGCGT"))
        assert 6 == len(result)
        assert result[0] == pytest.approx(-29.18363571, abs=5e-06)
        assert result[1] == pytest.approx(-38.3365097, abs=5e-06)
        assert result[2] == pytest.approx(-29.17756271, abs=5e-06)
        assert result[3] == pytest.approx(-38.04542542, abs=5e-06)
        assert result[4] == pytest.approx(-20.3014183, abs=5e-06)
        assert result[5] == pytest.approx(-25.18009186, abs=5e-06)

    def test_with_bad_char(self):
        """Test if Bio.motifs PWM scoring works with unexpected letters like N."""
        counts = self.m.counts
        pwm = counts.normalize(pseudocounts=0.25)
        pssm = pwm.log_odds()
        result = pssm.calculate(Seq("ACGTGTGCGTAGTGCGTN"))
        assert 7 == len(result)
        assert result[0] == pytest.approx(-29.18363571, abs=5e-06)
        assert result[1] == pytest.approx(-38.3365097, abs=5e-06)
        assert result[2] == pytest.approx(-29.17756271, abs=5e-06)
        assert result[3] == pytest.approx(-38.04542542, abs=5e-06)
        assert result[4] == pytest.approx(-20.3014183, abs=5e-06)
        assert result[5] == pytest.approx(-25.18009186, abs=5e-06)
        assert math.isnan(result[6]), f"Expected nan, not {result[6]!r}"

    def test_calculate_pseudocounts(self):
        pseudocounts = motifs.jaspar.calculate_pseudocounts(self.m)
        assert pseudocounts["A"] == pytest.approx(1.695582495781317, abs=5e-06)
        assert pseudocounts["C"] == pytest.approx(1.695582495781317, abs=5e-06)
        assert pseudocounts["G"] == pytest.approx(1.695582495781317, abs=5e-06)
        assert pseudocounts["T"] == pytest.approx(1.695582495781317, abs=5e-06)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
