# Copyright 2022 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Align.bigbed module."""
import os
import sys
import tempfile
import unittest
import pytest
from io import StringIO

from Bio import Align
from Bio import SeqIO
from Bio.Align import bigbed
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

np = pytest.importorskip("numpy")
for i, argument in enumerate(sys.argv):
    if argument == "--big":
        big = True
        sys.argv.pop(i)
        break
else:
    big = False


class TestAlign_dna_rna(unittest.TestCase):
    # The bigBed file dna_rna.bb was generated using the commands
    # sort -k1,1 -k2,2n dna_rna.bed > dna_rna.sorted.bed
    # twoBitInfo hg38.2bit hg38.chrom.sizes
    # bedToBigBed dna_rna.sorted.bed hg38.chrom.sizes dna_rna.bb

    path = "Blat/dna_rna.bb"

    def setUp(self):
        data = {}
        records = SeqIO.parse("Blat/dna.fa", "fasta")
        for record in records:
            name, start_end = record.id.split(":")
            assert name == "chr3"
            start, end = start_end.split("-")
            start = int(start)
            end = int(end)
            sequence = str(record.seq)
            assert len(sequence) == end - start
            data[start] = sequence
        self.dna = Seq(data, length=198295559)  # hg38 chr3
        records = SeqIO.parse("Blat/rna.fa", "fasta")
        self.rna = {record.id: record.seq for record in records}

    def test_reading(self):
        """Test parsing dna_rna.bb."""
        path = "Blat/dna_rna.bb"
        alignments = Align.parse(path, "bigbed")
        self.check_alignments(alignments)
        alignments = iter(alignments)
        self.check_alignments(alignments)
        with Align.parse(path, "bigbed") as alignments:
            self.check_alignments(alignments)
        with pytest.raises(AttributeError):
            alignments._stream
        with Align.parse(path, "bigbed") as alignments:
            pass
        with pytest.raises(AttributeError):
            alignments._stream

    def test_writing(self):
        """Test writing dna_rna.bb."""
        alignments = Align.parse(self.path, "bigbed")
        with tempfile.TemporaryFile() as output:
            Align.write(alignments, output, "bigbed")
            output.flush()
            output.seek(0)
            alignments = Align.parse(output, "bigbed")
            self.check_alignments(alignments)

    def check_alignments(self, alignments):
        assert str(alignments.declaration) == """\
table bed
"Browser Extensible Data"
(
   string          chrom;          "Reference sequence chromosome or scaffold"
   uint            chromStart;     "Start position in chromosome"
   uint            chromEnd;       "End position in chromosome"
   string          name;           "Name of item."
   uint            score;          "Score (0-1000)"
   char[1]         strand;         "+ or - for strand"
   uint            thickStart;     "Start of where display should be thick (start codon)"
   uint            thickEnd;       "End of where display should be thick (stop codon)"
   uint            reserved;       "Used as itemRgb as of 2004-11-22"
   int             blockCount;     "Number of blocks"
   int[blockCount] blockSizes;     "Comma separated list of block sizes"
   int[blockCount] chromStarts;    "Start positions relative to chromStart"
)
"""
        assert len(alignments.targets) == 1
        assert alignments.targets[0].id == "chr3"
        assert len(alignments.targets[0]) == 198295559
        assert len(alignments) == 4
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.shape == (2, 1711)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr3"
        assert alignment.query.id == "NR_046654.1"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[42530895, 42530958, 42532020,
                           42532095, 42532563, 42532606],
                          [     181,      118,      118,
                                 43,       43,        0]])
                # fmt: on
            )
        alignment.target.seq = self.dna
        alignment.query.seq = self.rna[alignment.query.id]
        assert np.array_equal(
                alignment.substitutions,
                # fmt: off
            np.array([[36.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0., 40.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0., 57.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0., 42.,  0.,  0.,  0.,  0.],
                      [ 2.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  3.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                     ])
            )
        assert alignment.substitutions.alphabet == "ACGTacgt"
        # The modified RNAs have gaps in their sequence. As this information is
        # not stored in a BED file, we cannot calculate the substitution matrix.
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (181 aligned letters; 175 identities; 6 mismatches; 1530 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 181:
        identities = 175,
        mismatches = 6.
    gaps = 1530:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1530:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1530:
                open_internal_deletions = 2,
                extend_internal_deletions = 1528;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 1530
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 1530
        assert counts.insertions == 0
        assert counts.deletions == 1530
        assert counts.gaps == 1530
        assert counts.aligned == 181
        assert counts.identities == 175
        assert counts.mismatches == 6
        alignment = next(alignments)
        assert alignment.score == 978
        assert alignment.shape == (2, 1711)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr3"
        assert alignment.query.id == "NR_046654.1_modified"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[42530895, 42530922, 42530958, 42532020, 42532037,
                           42532039, 42532095, 42532563, 42532606],
                          [     179,      152,      116,      116,       99,
                                 99,       43,       43,        0]])
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (179 aligned letters; 0 identities; 0 mismatches; 1532 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 179:
        identities = 0,
        mismatches = 0.
    gaps = 1532:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1532:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1532:
                open_internal_deletions = 3,
                extend_internal_deletions = 1529;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 1532
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 1532
        assert counts.insertions == 0
        assert counts.deletions == 1532
        assert counts.gaps == 1532
        assert counts.aligned == 179
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.shape == (2, 5407)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr3"
        assert alignment.query.id == "NR_111921.1"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[48663767, 48663813, 48665640,
                           48665722, 48669098, 48669174],
                          [       0,       46,       46,
                                128,      128,      204]])
                # fmt: on
            )
        alignment.target.seq = self.dna
        alignment.query.seq = self.rna[alignment.query.id]
        assert np.array_equal(
                alignment.substitutions,
                # fmt: off
            np.array([[53.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0., 35.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0., 50.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0., 27.,  0.,  0.,  0.,  0.],
                      [ 9.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  7.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0., 16.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  7.,  0.,  0.,  0.,  0.],
                     ])
            )
        assert alignment.substitutions.alphabet == "ACGTacgt"
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (204 aligned letters; 165 identities; 39 mismatches; 5203 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 204:
        identities = 165,
        mismatches = 39.
    gaps = 5203:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 5203:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 5203:
                open_internal_deletions = 2,
                extend_internal_deletions = 5201;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 5203
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 5203
        assert counts.insertions == 0
        assert counts.deletions == 5203
        assert counts.gaps == 5203
        assert counts.aligned == 204
        assert counts.identities == 165
        assert counts.mismatches == 39
        alignment = next(alignments)
        assert alignment.score == 972
        assert alignment.shape == (2, 5407)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr3"
        assert alignment.query.id == "NR_111921.1_modified"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[48663767, 48663795, 48663796, 48663813, 48665640,
                           48665716, 48665722, 48669098, 48669174],
                          [       0,       28,       28,       45,       45,
                                121,      127,      127,      203]])
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (203 aligned letters; 0 identities; 0 mismatches; 5204 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 203:
        identities = 0,
        mismatches = 0.
    gaps = 5204:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 5204:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 5204:
                open_internal_deletions = 3,
                extend_internal_deletions = 5201;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 5204
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 5204
        assert counts.insertions == 0
        assert counts.deletions == 5204
        assert counts.gaps == 5204
        assert counts.aligned == 203
        with pytest.raises(StopIteration):
            next(alignments)


class TestAlign_dna(unittest.TestCase):
    # The bigBed file psl_34_001.bb was generated using the commands
    # sort -k1,1 -k2,2n psl_34_001.bed > psl_34_001.sorted.bed
    # twoBitInfo hg19.2bit hg19.chrom.sizes
    # bedToBigBed psl_34_001.sorted.bed hg19.chrom.sizes psl_34_001.bb

    # The bigBed file psl_34_003.bb was generated using the commands
    # sort -k1,1 -k2,2n psl_34_003.bed > psl_34_003.sorted.bed
    # twoBitInfo hg19.2bit hg19.chrom.sizes
    # bedToBigBed psl_34_003.sorted.bed hg19.chrom.sizes psl_34_003.bb

    # The bigBed file psl_34_004.bb was generated using the commands
    # sort -k1,1 -k2,2n psl_34_004.bed > psl_34_004.sorted.bed
    # twoBitInfo hg19.2bit hg19.chrom.sizes
    # bedToBigBed psl_34_004.sorted.bed hg19.chrom.sizes psl_34_004.bb

    # The bigBed file psl_34_005.bb was generated using the commands
    # sort -k1,1 -k2,2n psl_34_005.bed > psl_34_005.sorted.bed
    # twoBitInfo hg19.2bit hg19.chrom.sizes
    # bedToBigBed psl_34_005.sorted.bed hg19.chrom.sizes psl_34_005.bb

    def check_alignments_psl_34_001(self, alignments):
        assert str(alignments.declaration) == """\
table bed
"Browser Extensible Data"
(
   string          chrom;          "Reference sequence chromosome or scaffold"
   uint            chromStart;     "Start position in chromosome"
   uint            chromEnd;       "End position in chromosome"
   string          name;           "Name of item."
   uint            score;          "Score (0-1000)"
   char[1]         strand;         "+ or - for strand"
   uint            thickStart;     "Start of where display should be thick (start codon)"
   uint            thickEnd;       "End of where display should be thick (stop codon)"
   uint            reserved;       "Used as itemRgb as of 2004-11-22"
   int             blockCount;     "Number of blocks"
   int[blockCount] blockSizes;     "Comma separated list of block sizes"
   int[blockCount] chromStarts;    "Start positions relative to chromStart"
)
"""
        assert len(alignments.targets) == 10
        assert alignments.targets[0].id == "chr1"
        assert len(alignments.targets[0]) == 249250621
        assert alignments.targets[1].id == "chr10"
        assert len(alignments.targets[1]) == 135534747
        assert alignments.targets[2].id == "chr13"
        assert len(alignments.targets[2]) == 115169878
        assert alignments.targets[3].id == "chr18"
        assert len(alignments.targets[3]) == 78077248
        assert alignments.targets[4].id == "chr19"
        assert len(alignments.targets[4]) == 59128983
        assert alignments.targets[5].id == "chr2"
        assert len(alignments.targets[5]) == 243199373
        assert alignments.targets[6].id == "chr22"
        assert len(alignments.targets[6]) == 51304566
        assert alignments.targets[7].id == "chr4"
        assert len(alignments.targets[7]) == 191154276
        assert alignments.targets[8].id == "chr8"
        assert len(alignments.targets[8]) == 146364022
        assert alignments.targets[9].id == "chr9"
        assert len(alignments.targets[9]) == 141213431
        assert len(alignments) == 22
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.shape == (2, 50)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr1"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[1207056, 1207106],
                          [      0,      50]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (50 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 50:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 50
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.shape == (2, 33)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr1"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[10271783, 10271816],
                          [       0,       33]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (33 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 33:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 33
        alignment = next(alignments)
        assert alignment.score == 946
        assert alignment.shape == (2, 36)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr1"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[39368490, 39368526],
                          [      36,        0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 36
        alignment = next(alignments)
        assert alignment.score == 824
        assert alignment.shape == (2, 34)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr1"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[61700837, 61700871],
                          [       0,       34]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (34 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 34:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 34
        alignment = next(alignments)
        assert alignment.score == 942
        assert alignment.shape == (2, 34)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr1"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[220325687, 220325721],
                          [       34,         0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (34 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 34:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 34
        alignment = next(alignments)
        assert alignment.score == 834
        assert alignment.shape == (2, 36)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr10"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[99388555, 99388591],
                          [      36,        0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 36
        alignment = next(alignments)
        assert alignment.score == 920
        assert alignment.shape == (2, 25)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr10"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[112178171, 112178196],
                          [       25,         0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (25 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 25:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 25
        alignment = next(alignments)
        assert alignment.score == 912
        assert alignment.shape == (2, 51)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr13"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[52759147, 52759154, 52759160, 52759198],
                          [       0,        7,        7,       45]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (45 aligned letters; 0 identities; 0 mismatches; 6 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 45:
        identities = 0,
        mismatches = 0.
    gaps = 6:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 6:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 6:
                open_internal_deletions = 1,
                extend_internal_deletions = 5;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 6
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 6
        assert counts.insertions == 0
        assert counts.deletions == 6
        assert counts.gaps == 6
        assert counts.aligned == 45
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.shape == (2, 39)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr18"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[23891310, 23891349],
                          [       0,       39]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (39 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 39:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 39
        alignment = next(alignments)
        assert alignment.score == 930
        assert alignment.shape == (2, 28)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr18"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[43252217, 43252245],
                          [       0,       28]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (28 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 28:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 28
        alignment = next(alignments)
        assert alignment.score == 848
        assert alignment.shape == (2, 39)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr19"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[553742, 553781],
                          [    39,      0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (39 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 39:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 39
        alignment = next(alignments)
        assert alignment.score == 890
        assert alignment.shape == (2, 170)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr19"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[35483340, 35483365, 35483499, 35483510],
                          [       0,       25,       25,       36]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 134 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 134:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 134:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 134:
                open_internal_deletions = 1,
                extend_internal_deletions = 133;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 134
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 134
        assert counts.insertions == 0
        assert counts.deletions == 134
        assert counts.gaps == 134
        assert counts.aligned == 36
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.shape == (2, 39)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr19"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[54017130, 54017169],
                          [      39,        0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (39 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 39:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 39
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.shape == (2, 17)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr2"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[53575980, 53575997],
                          [      17,        0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (17 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 17:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 17
        alignment = next(alignments)
        assert alignment.score == 946
        assert alignment.shape == (2, 36)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr2"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[120641740, 120641776],
                          [       36,         0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 36
        alignment = next(alignments)
        assert alignment.score == 682
        assert alignment.shape == (2, 44)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr2"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[183925984, 183925990, 183926028],
                          [        0,         6,        44]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (44 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 44:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 44
        alignment = next(alignments)
        assert alignment.score == 834
        assert alignment.shape == (2, 36)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr22"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[42144400, 42144436],
                          [       0,       36]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 36
        alignment = next(alignments)
        assert alignment.score == 892
        assert alignment.shape == (2, 37)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr22"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[48997405, 48997442],
                          [      37,        0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (37 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 37:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 37
        alignment = next(alignments)
        assert alignment.score == 572
        assert alignment.shape == (2, 34)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr4"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[37558157, 37558167, 37558173, 37558191],
                          [      28,       18,       18,        0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (28 aligned letters; 0 identities; 0 mismatches; 6 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 28:
        identities = 0,
        mismatches = 0.
    gaps = 6:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 6:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 6:
                open_internal_deletions = 1,
                extend_internal_deletions = 5;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 6
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 6
        assert counts.insertions == 0
        assert counts.deletions == 6
        assert counts.gaps == 6
        assert counts.aligned == 28
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.shape == (2, 16)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr4"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[61646095, 61646111],
                          [       0,       16]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (16 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 16:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 16
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.shape == (2, 41)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr8"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[95160479, 95160520],
                          [       0,       41]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (41 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 41:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 41
        alignment = next(alignments)
        assert alignment.score == 854
        assert alignment.shape == (2, 41)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr9"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[85737865, 85737906],
                          [       0,       41]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (41 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 41:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 41
        with pytest.raises(StopIteration):
            next(alignments)

    def test_reading_psl_34_001(self):
        """Test reading psl_34_001.bb."""
        path = "Blat/psl_34_001.bb"
        alignments = Align.parse(path, "bigbed")
        self.check_alignments_psl_34_001(alignments)

    def test_writing_psl_34_001(self):
        """Test writing psl_34_001.bb."""
        path = "Blat/psl_34_001.bb"
        alignments = Align.parse(path, "bigbed")
        with tempfile.TemporaryFile() as output:
            Align.write(alignments, output, "bigbed")
            output.flush()
            output.seek(0)
            alignments = Align.parse(output, "bigbed")
            self.check_alignments_psl_34_001(alignments)

    def check_alignments_psl_34_003(self, alignments):
        assert str(alignments.declaration) == """\
table bed
"Browser Extensible Data"
(
   string          chrom;          "Reference sequence chromosome or scaffold"
   uint            chromStart;     "Start position in chromosome"
   uint            chromEnd;       "End position in chromosome"
   string          name;           "Name of item."
   uint            score;          "Score (0-1000)"
   char[1]         strand;         "+ or - for strand"
   uint            thickStart;     "Start of where display should be thick (start codon)"
   uint            thickEnd;       "End of where display should be thick (stop codon)"
   uint            reserved;       "Used as itemRgb as of 2004-11-22"
   int             blockCount;     "Number of blocks"
   int[blockCount] blockSizes;     "Comma separated list of block sizes"
   int[blockCount] chromStarts;    "Start positions relative to chromStart"
)
"""
        assert len(alignments.targets) == 3
        assert alignments.targets[0].id == "chr1"
        assert len(alignments.targets[0]) == 249250621
        assert alignments.targets[1].id == "chr2"
        assert len(alignments.targets[1]) == 243199373
        assert alignments.targets[2].id == "chr4"
        assert len(alignments.targets[2]) == 191154276
        assert len(alignments) == 3
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.shape == (2, 33)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr1"
        assert alignment.query.id == "hg18_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[10271783, 10271816],
                          [       0,       33]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (33 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 33:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 33
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.shape == (2, 17)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr2"
        assert alignment.query.id == "hg18_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[53575980, 53575997],
                          [      17,        0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (17 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 17:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 17
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.shape == (2, 16)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr4"
        assert alignment.query.id == "hg18_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[61646095, 61646111],
                          [       0,       16]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (16 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 16:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 16
        with pytest.raises(StopIteration):
            next(alignments)

    def test_reading_psl_34_003(self):
        """Test reading psl_34_003.bb."""
        path = "Blat/psl_34_003.bb"
        alignments = Align.parse(path, "bigbed")
        self.check_alignments_psl_34_003(alignments)

    def test_writing_psl_34_003(self):
        """Test writing psl_34_003.bb."""
        path = "Blat/psl_34_003.bb"
        alignments = Align.parse(path, "bigbed")
        with tempfile.TemporaryFile() as output:
            Align.write(alignments, output, "bigbed")
            output.flush()
            output.seek(0)
            alignments = Align.parse(output, "bigbed")
            self.check_alignments_psl_34_003(alignments)

    def check_alignments_psl_34_004(self, alignments):
        assert str(alignments.declaration) == """\
table bed
"Browser Extensible Data"
(
   string          chrom;          "Reference sequence chromosome or scaffold"
   uint            chromStart;     "Start position in chromosome"
   uint            chromEnd;       "End position in chromosome"
   string          name;           "Name of item."
   uint            score;          "Score (0-1000)"
   char[1]         strand;         "+ or - for strand"
   uint            thickStart;     "Start of where display should be thick (start codon)"
   uint            thickEnd;       "End of where display should be thick (stop codon)"
   uint            reserved;       "Used as itemRgb as of 2004-11-22"
   int             blockCount;     "Number of blocks"
   int[blockCount] blockSizes;     "Comma separated list of block sizes"
   int[blockCount] chromStarts;    "Start positions relative to chromStart"
)
"""
        assert len(alignments.targets) == 10
        assert alignments.targets[0].id == "chr1"
        assert len(alignments.targets[0]) == 249250621
        assert alignments.targets[1].id == "chr10"
        assert len(alignments.targets[1]) == 135534747
        assert alignments.targets[2].id == "chr13"
        assert len(alignments.targets[2]) == 115169878
        assert alignments.targets[3].id == "chr18"
        assert len(alignments.targets[3]) == 78077248
        assert alignments.targets[4].id == "chr19"
        assert len(alignments.targets[4]) == 59128983
        assert alignments.targets[5].id == "chr2"
        assert len(alignments.targets[5]) == 243199373
        assert alignments.targets[6].id == "chr22"
        assert len(alignments.targets[6]) == 51304566
        assert alignments.targets[7].id == "chr4"
        assert len(alignments.targets[7]) == 191154276
        assert alignments.targets[8].id == "chr8"
        assert len(alignments.targets[8]) == 146364022
        assert alignments.targets[9].id == "chr9"
        assert len(alignments.targets[9]) == 141213431
        assert len(alignments) == 19
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.shape == (2, 50)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr1"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[1207056, 1207106],
                          [      0,      50]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (50 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 50:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 50
        alignment = next(alignments)
        assert alignment.score == 946
        assert alignment.shape == (2, 36)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr1"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[39368490, 39368526],
                          [      36,        0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 36
        alignment = next(alignments)
        assert alignment.score == 824
        assert alignment.shape == (2, 34)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr1"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[61700837, 61700871],
                          [       0,       34]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (34 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 34:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 34
        alignment = next(alignments)
        assert alignment.score == 942
        assert alignment.shape == (2, 34)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr1"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[220325687, 220325721],
                          [       34,         0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (34 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 34:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 34
        alignment = next(alignments)
        assert alignment.score == 834
        assert alignment.shape == (2, 36)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr10"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[99388555, 99388591],
                          [      36,        0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 36
        alignment = next(alignments)
        assert alignment.score == 920
        assert alignment.shape == (2, 25)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr10"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[112178171, 112178196],
                          [       25,         0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (25 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 25:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 25
        alignment = next(alignments)
        assert alignment.score == 912
        assert alignment.shape == (2, 51)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr13"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[52759147, 52759154, 52759160, 52759198],
                          [       0,        7,        7,       45]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (45 aligned letters; 0 identities; 0 mismatches; 6 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 45:
        identities = 0,
        mismatches = 0.
    gaps = 6:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 6:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 6:
                open_internal_deletions = 1,
                extend_internal_deletions = 5;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 6
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 6
        assert counts.insertions == 0
        assert counts.deletions == 6
        assert counts.gaps == 6
        assert counts.aligned == 45
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.shape == (2, 39)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr18"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[23891310, 23891349],
                          [       0,       39]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (39 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 39:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 39
        alignment = next(alignments)
        assert alignment.score == 930
        assert alignment.shape == (2, 28)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr18"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[43252217, 43252245],
                          [       0,       28]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (28 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 28:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 28
        alignment = next(alignments)
        assert alignment.score == 848
        assert alignment.shape == (2, 39)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr19"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[553742, 553781],
                          [    39,      0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (39 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 39:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 39
        alignment = next(alignments)
        assert alignment.score == 890
        assert alignment.shape == (2, 170)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr19"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[35483340, 35483365, 35483499, 35483510],
                          [       0,       25,       25,       36]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 134 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 134:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 134:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 134:
                open_internal_deletions = 1,
                extend_internal_deletions = 133;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 134
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 134
        assert counts.insertions == 0
        assert counts.deletions == 134
        assert counts.gaps == 134
        assert counts.aligned == 36
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.shape == (2, 39)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr19"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[54017130, 54017169],
                          [      39,        0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (39 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 39:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 39
        alignment = next(alignments)
        assert alignment.score == 946
        assert alignment.shape == (2, 36)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr2"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[120641740, 120641776],
                          [       36,         0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 36
        alignment = next(alignments)
        assert alignment.score == 682
        assert alignment.shape == (2, 44)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr2"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[183925984, 183925990, 183926028],
                          [        0,         6,        44]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (44 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 44:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 44
        alignment = next(alignments)
        assert alignment.score == 834
        assert alignment.shape == (2, 36)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr22"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[42144400, 42144436],
                          [       0,       36]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 36
        alignment = next(alignments)
        assert alignment.score == 892
        assert alignment.shape == (2, 37)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr22"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[48997405, 48997442],
                          [      37,        0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (37 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 37:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 37
        alignment = next(alignments)
        assert alignment.score == 572
        assert alignment.shape == (2, 34)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr4"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[37558157, 37558167, 37558173, 37558191],
                          [      28,       18,       18,        0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (28 aligned letters; 0 identities; 0 mismatches; 6 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 28:
        identities = 0,
        mismatches = 0.
    gaps = 6:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 6:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 6:
                open_internal_deletions = 1,
                extend_internal_deletions = 5;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 6
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 6
        assert counts.insertions == 0
        assert counts.deletions == 6
        assert counts.gaps == 6
        assert counts.aligned == 28
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.shape == (2, 41)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr8"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[95160479, 95160520],
                          [       0,       41]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (41 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 41:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 41
        alignment = next(alignments)
        assert alignment.score == 854
        assert alignment.shape == (2, 41)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr9"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[85737865, 85737906],
                          [       0,       41]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (41 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 41:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 41
        with pytest.raises(StopIteration):
            next(alignments)

    def test_reading_psl_34_004(self):
        """Test reading psl_34_004.bb."""
        path = "Blat/psl_34_004.bb"
        alignments = Align.parse(path, "bigbed")
        self.check_alignments_psl_34_004(alignments)

    def test_writing_psl_34_004(self):
        """Test writing psl_34_004.bb."""
        path = "Blat/psl_34_004.bb"
        alignments = Align.parse(path, "bigbed")
        with tempfile.TemporaryFile() as output:
            Align.write(alignments, output, "bigbed")
            output.flush()
            output.seek(0)
            alignments = Align.parse(output, "bigbed")
            self.check_alignments_psl_34_004(alignments)

    def check_alignments_psl_34_005(self, alignments):
        assert str(alignments.declaration) == """\
table bed
"Browser Extensible Data"
(
   string          chrom;          "Reference sequence chromosome or scaffold"
   uint            chromStart;     "Start position in chromosome"
   uint            chromEnd;       "End position in chromosome"
   string          name;           "Name of item."
   uint            score;          "Score (0-1000)"
   char[1]         strand;         "+ or - for strand"
   uint            thickStart;     "Start of where display should be thick (start codon)"
   uint            thickEnd;       "End of where display should be thick (stop codon)"
   uint            reserved;       "Used as itemRgb as of 2004-11-22"
   int             blockCount;     "Number of blocks"
   int[blockCount] blockSizes;     "Comma separated list of block sizes"
   int[blockCount] chromStarts;    "Start positions relative to chromStart"
)
"""
        assert len(alignments.targets) == 10
        assert alignments.targets[0].id == "chr1"
        assert len(alignments.targets[0]) == 249250621
        assert alignments.targets[1].id == "chr10"
        assert len(alignments.targets[1]) == 135534747
        assert alignments.targets[2].id == "chr13"
        assert len(alignments.targets[2]) == 115169878
        assert alignments.targets[3].id == "chr18"
        assert len(alignments.targets[3]) == 78077248
        assert alignments.targets[4].id == "chr19"
        assert len(alignments.targets[4]) == 59128983
        assert alignments.targets[5].id == "chr2"
        assert len(alignments.targets[5]) == 243199373
        assert alignments.targets[6].id == "chr22"
        assert len(alignments.targets[6]) == 51304566
        assert alignments.targets[7].id == "chr4"
        assert len(alignments.targets[7]) == 191154276
        assert alignments.targets[8].id == "chr8"
        assert len(alignments.targets[8]) == 146364022
        assert alignments.targets[9].id == "chr9"
        assert len(alignments.targets[9]) == 141213431
        assert len(alignments) == 22
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.shape == (2, 50)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr1"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[1207056, 1207106],
                          [      0,      50]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (50 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 50:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 50
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.shape == (2, 33)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr1"
        assert alignment.query.id == "hg18_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[10271783, 10271816],
                          [       0,       33]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (33 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 33:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 33
        alignment = next(alignments)
        assert alignment.score == 946
        assert alignment.shape == (2, 36)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr1"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[39368490, 39368526],
                          [      36,        0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 36
        alignment = next(alignments)
        assert alignment.score == 824
        assert alignment.shape == (2, 34)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr1"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[61700837, 61700871],
                          [       0,       34]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (34 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 34:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 34
        alignment = next(alignments)
        assert alignment.score == 942
        assert alignment.shape == (2, 34)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr1"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[220325687, 220325721],
                          [       34,         0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (34 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 34:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 34
        alignment = next(alignments)
        assert alignment.score == 834
        assert alignment.shape == (2, 36)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr10"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[99388555, 99388591],
                          [      36,        0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 36
        alignment = next(alignments)
        assert alignment.score == 920
        assert alignment.shape == (2, 25)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr10"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[112178171, 112178196],
                          [       25,         0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (25 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 25:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 25
        alignment = next(alignments)
        assert alignment.score == 912
        assert alignment.shape == (2, 51)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr13"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[52759147, 52759154, 52759160, 52759198],
                          [       0,        7,        7,       45]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (45 aligned letters; 0 identities; 0 mismatches; 6 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 45:
        identities = 0,
        mismatches = 0.
    gaps = 6:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 6:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 6:
                open_internal_deletions = 1,
                extend_internal_deletions = 5;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 6
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 6
        assert counts.insertions == 0
        assert counts.deletions == 6
        assert counts.gaps == 6
        assert counts.aligned == 45
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.shape == (2, 39)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr18"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[23891310, 23891349],
                          [       0,       39]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (39 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 39:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 39
        alignment = next(alignments)
        assert alignment.score == 930
        assert alignment.shape == (2, 28)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr18"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[43252217, 43252245],
                          [       0,       28]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (28 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 28:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 28
        alignment = next(alignments)
        assert alignment.score == 848
        assert alignment.shape == (2, 39)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr19"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[553742, 553781],
                          [    39,      0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (39 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 39:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 39
        alignment = next(alignments)
        assert alignment.score == 890
        assert alignment.shape == (2, 170)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr19"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[35483340, 35483365, 35483499, 35483510],
                          [       0,       25,       25,       36]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 134 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 134:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 134:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 134:
                open_internal_deletions = 1,
                extend_internal_deletions = 133;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 134
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 134
        assert counts.insertions == 0
        assert counts.deletions == 134
        assert counts.gaps == 134
        assert counts.aligned == 36
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.shape == (2, 39)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr19"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[54017130, 54017169],
                          [      39,        0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (39 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 39:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 39
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.shape == (2, 17)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr2"
        assert alignment.query.id == "hg18_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[53575980, 53575997],
                          [      17,        0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (17 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 17:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 17
        alignment = next(alignments)
        assert alignment.score == 946
        assert alignment.shape == (2, 36)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr2"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[120641740, 120641776],
                          [       36,         0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 36
        alignment = next(alignments)
        assert alignment.score == 682
        assert alignment.shape == (2, 44)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr2"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[183925984, 183925990, 183926028],
                          [        0,         6,        44]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (44 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 44:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 44
        alignment = next(alignments)
        assert alignment.score == 834
        assert alignment.shape == (2, 36)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr22"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[42144400, 42144436],
                          [       0,       36]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 36
        alignment = next(alignments)
        assert alignment.score == 892
        assert alignment.shape == (2, 37)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr22"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[48997405, 48997442],
                          [      37,        0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (37 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 37:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 37
        alignment = next(alignments)
        assert alignment.score == 572
        assert alignment.shape == (2, 34)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr4"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[37558157, 37558167, 37558173, 37558191],
                          [      28,       18,       18,        0]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (28 aligned letters; 0 identities; 0 mismatches; 6 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 28:
        identities = 0,
        mismatches = 0.
    gaps = 6:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 6:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 6:
                open_internal_deletions = 1,
                extend_internal_deletions = 5;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 6
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 6
        assert counts.insertions == 0
        assert counts.deletions == 6
        assert counts.gaps == 6
        assert counts.aligned == 28
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.shape == (2, 16)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr4"
        assert alignment.query.id == "hg18_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[61646095, 61646111],
                          [       0,       16]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (16 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 16:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 16
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.shape == (2, 41)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr8"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[95160479, 95160520],
                          [       0,       41]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (41 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 41:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 41
        alignment = next(alignments)
        assert alignment.score == 854
        assert alignment.shape == (2, 41)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr9"
        assert alignment.query.id == "hg19_dna"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[85737865, 85737906],
                          [       0,       41]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (41 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 41:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 41
        with pytest.raises(StopIteration):
            next(alignments)

    def test_reading_psl_34_005(self):
        """Test reading psl_34_005.bb."""
        path = "Blat/psl_34_005.bb"
        alignments = Align.parse(path, "bigbed")
        self.check_alignments_psl_34_005(alignments)

    def test_writing_psl_34_005(self):
        """Test writing psl_34_005.bb."""
        path = "Blat/psl_34_005.bb"
        alignments = Align.parse(path, "bigbed")
        with tempfile.TemporaryFile() as output:
            Align.write(alignments, output, "bigbed")
            output.flush()
            output.seek(0)
            alignments = Align.parse(output, "bigbed")
            self.check_alignments_psl_34_005(alignments)


class TestAlign_dnax_prot(unittest.TestCase):
    # The bigBed file psl_35_001.bb was generated using the commands
    # sort -k1,1 -k2,2n psl_35_001.bed > psl_35_001.sorted.bed
    # twoBitInfo hg38.2bit hg38.chrom.sizes
    # bedToBigBed psl_35_001.sorted.bed hg38.chrom.sizes psl_35_001.bb

    path = "Blat/psl_35_001.bb"

    def check_alignments(self, alignments):
        assert str(alignments.declaration) == """\
table bed
"Browser Extensible Data"
(
   string          chrom;          "Reference sequence chromosome or scaffold"
   uint            chromStart;     "Start position in chromosome"
   uint            chromEnd;       "End position in chromosome"
   string          name;           "Name of item."
   uint            score;          "Score (0-1000)"
   char[1]         strand;         "+ or - for strand"
   uint            thickStart;     "Start of where display should be thick (start codon)"
   uint            thickEnd;       "End of where display should be thick (stop codon)"
   uint            reserved;       "Used as itemRgb as of 2004-11-22"
   int             blockCount;     "Number of blocks"
   int[blockCount] blockSizes;     "Comma separated list of block sizes"
   int[blockCount] chromStarts;    "Start positions relative to chromStart"
)
"""
        assert len(alignments.targets) == 2
        assert alignments.targets[0].id == "chr13"
        assert len(alignments.targets[0]) == 114364328
        assert alignments.targets[1].id == "chr4"
        assert len(alignments.targets[1]) == 190214555
        assert len(alignments) == 8
        alignment = next(alignments)
        assert alignment.score == 986
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr13"
        assert alignment.query.id == "CAG33136.1"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[75549820, 75549865, 75567225, 75567312],
                          [       0,       45,       45,      132]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (132 aligned letters; 0 identities; 0 mismatches; 17360 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 132:
        identities = 0,
        mismatches = 0.
    gaps = 17360:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 17360:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 17360:
                open_internal_deletions = 1,
                extend_internal_deletions = 17359;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 17360
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 17360
        assert counts.insertions == 0
        assert counts.deletions == 17360
        assert counts.gaps == 17360
        assert counts.aligned == 132
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr13"
        assert alignment.query.id == "CAG33136.1"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[75560749, 75560881],
                          [       0,      132]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (132 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 132:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 132
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr13"
        assert alignment.query.id == "CAG33136.1"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[75566694, 75566850],
                          [       0,      156]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (156 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 156:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 156
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr13"
        assert alignment.query.id == "CAG33136.1"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[75569459, 75569507],
                          [       0,       48]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (48 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 48:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 48
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr13"
        assert alignment.query.id == "CAG33136.1"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[75594914, 75594989],
                          [       0,       75]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (75 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 75:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 75
        alignment = next(alignments)
        assert alignment.score == 1000
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr13"
        assert alignment.query.id == "CAG33136.1"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[75604767, 75604827, 75605728, 75605809],
                          [       0,       60,       60,      141]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (141 aligned letters; 0 identities; 0 mismatches; 901 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 141:
        identities = 0,
        mismatches = 0.
    gaps = 901:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 901:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 901:
                open_internal_deletions = 1,
                extend_internal_deletions = 900;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 901
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 901
        assert counts.insertions == 0
        assert counts.deletions == 901
        assert counts.gaps == 901
        assert counts.aligned == 141
        alignment = next(alignments)
        assert alignment.score == 166
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr4"
        assert alignment.query.id == "CAG33136.1"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[41257605, 41257731, 41263227, 41263290],
                          [       0,      126,      126,      189]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (189 aligned letters; 0 identities; 0 mismatches; 5496 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 189:
        identities = 0,
        mismatches = 0.
    gaps = 5496:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 5496:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 5496:
                open_internal_deletions = 1,
                extend_internal_deletions = 5495;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 5496
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 5496
        assert counts.insertions == 0
        assert counts.deletions == 5496
        assert counts.gaps == 5496
        assert counts.aligned == 189
        alignment = next(alignments)
        assert alignment.score == 530
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr4"
        assert alignment.query.id == "CAG33136.1"
        assert np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[41260685, 41260787],
                          [       0,      102]]),
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (102 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 102:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 102
        with pytest.raises(StopIteration):
            next(alignments)

    def test_reading_psl_35_001(self):
        """Test parsing psl_35_001.bb."""
        alignments = Align.parse(self.path, "bigbed")
        self.check_alignments(alignments)

    def test_writing_psl_35_001(self):
        """Test writing psl_35_001.bb."""
        alignments = Align.parse(self.path, "bigbed")
        with tempfile.TemporaryFile() as output:
            Align.write(alignments, output, "bigbed")
            output.flush()
            output.seek(0)
            alignments = Align.parse(output, "bigbed")
            self.check_alignments(alignments)


class TestAlign_bed12(unittest.TestCase):
    # The bigBed files were generated using the commands
    # twoBitInfo hg19.2bit hg19.chrom.sizes
    # bedToBigBed bed3.bed hg19.chrom.sizes bed3.bb
    # bedToBigBed bed4.bed hg19.chrom.sizes bed4.bb
    # bedToBigBed bed5.bed hg19.chrom.sizes bed5.bb
    # bedToBigBed bed6.bed hg19.chrom.sizes bed6.bb
    # bedToBigBed bed7.bed hg19.chrom.sizes bed7.bb
    # bedToBigBed bed8.bed hg19.chrom.sizes bed8.bb
    # bedToBigBed bed9.bed hg19.chrom.sizes bed9.bb
    # bedToBigBed bed12.bed hg19.chrom.sizes bed12.bb

    def check_autosql(self, declaration, bedN, msg):
        if bedN == 3:
            assert str(declaration) == """\
table bed
"Browser Extensible Data"
(
   string chrom;         "Reference sequence chromosome or scaffold"
   uint   chromStart;    "Start position in chromosome"
   uint   chromEnd;      "End position in chromosome"
)
""", msg
        elif bedN == 4:
            assert str(declaration) == """\
table bed
"Browser Extensible Data"
(
   string chrom;         "Reference sequence chromosome or scaffold"
   uint   chromStart;    "Start position in chromosome"
   uint   chromEnd;      "End position in chromosome"
   string name;          "Name of item."
)
""", msg
        elif bedN == 5:
            assert str(declaration) == """\
table bed
"Browser Extensible Data"
(
   string chrom;         "Reference sequence chromosome or scaffold"
   uint   chromStart;    "Start position in chromosome"
   uint   chromEnd;      "End position in chromosome"
   string name;          "Name of item."
   uint   score;         "Score (0-1000)"
)
""", msg
        elif bedN == 6:
            assert str(declaration) == """\
table bed
"Browser Extensible Data"
(
   string  chrom;         "Reference sequence chromosome or scaffold"
   uint    chromStart;    "Start position in chromosome"
   uint    chromEnd;      "End position in chromosome"
   string  name;          "Name of item."
   uint    score;         "Score (0-1000)"
   char[1] strand;        "+ or - for strand"
)
""", msg
        elif bedN == 7:
            assert str(declaration) == """\
table bed
"Browser Extensible Data"
(
   string  chrom;         "Reference sequence chromosome or scaffold"
   uint    chromStart;    "Start position in chromosome"
   uint    chromEnd;      "End position in chromosome"
   string  name;          "Name of item."
   uint    score;         "Score (0-1000)"
   char[1] strand;        "+ or - for strand"
   uint    thickStart;    "Start of where display should be thick (start codon)"
)
""", msg
        elif bedN == 8:
            assert str(declaration) == """\
table bed
"Browser Extensible Data"
(
   string  chrom;         "Reference sequence chromosome or scaffold"
   uint    chromStart;    "Start position in chromosome"
   uint    chromEnd;      "End position in chromosome"
   string  name;          "Name of item."
   uint    score;         "Score (0-1000)"
   char[1] strand;        "+ or - for strand"
   uint    thickStart;    "Start of where display should be thick (start codon)"
   uint    thickEnd;      "End of where display should be thick (stop codon)"
)
""", msg
        elif bedN == 9:
            assert str(declaration) == """\
table bed
"Browser Extensible Data"
(
   string  chrom;         "Reference sequence chromosome or scaffold"
   uint    chromStart;    "Start position in chromosome"
   uint    chromEnd;      "End position in chromosome"
   string  name;          "Name of item."
   uint    score;         "Score (0-1000)"
   char[1] strand;        "+ or - for strand"
   uint    thickStart;    "Start of where display should be thick (start codon)"
   uint    thickEnd;      "End of where display should be thick (stop codon)"
   uint    reserved;      "Used as itemRgb as of 2004-11-22"
)
""", msg
        elif bedN == 12:
            assert str(declaration) == """\
table bed
"Browser Extensible Data"
(
   string          chrom;          "Reference sequence chromosome or scaffold"
   uint            chromStart;     "Start position in chromosome"
   uint            chromEnd;       "End position in chromosome"
   string          name;           "Name of item."
   uint            score;          "Score (0-1000)"
   char[1]         strand;         "+ or - for strand"
   uint            thickStart;     "Start of where display should be thick (start codon)"
   uint            thickEnd;       "End of where display should be thick (stop codon)"
   uint            reserved;       "Used as itemRgb as of 2004-11-22"
   int             blockCount;     "Number of blocks"
   int[blockCount] blockSizes;     "Comma separated list of block sizes"
   int[blockCount] chromStarts;    "Start positions relative to chromStart"
)
""", msg

    def check_alignments(self, alignments, bedN, msg):
        assert len(alignments) == 2
        alignment = next(alignments)
        if bedN >= 5:
            assert alignment.score == 960, msg
        assert alignment.shape == (2, 4000), msg
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1], msg
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1], msg
        assert len(alignment) == 2, msg
        assert alignment.sequences[0] is alignment.target, msg
        assert alignment.sequences[1] is alignment.query, msg
        assert alignment.target.id == "chr22", msg
        if bedN >= 4:
            assert alignment.query.id == "mRNA1", msg
        else:
            assert alignment.query.id is None, msg
        if bedN == 12:
            assert np.array_equal(
                    alignment.coordinates,
                    # fmt: off
                    np.array([[1000, 1567, 4512, 5000],
                              [   0,  567,  567, 1055]]),
                    # fmt: on
                ), msg
        else:
            assert np.array_equal(
                    alignment.coordinates,
                    np.array([[1000, 5000], [0, 4000]]),
                ), msg
        if bedN >= 7:
            assert alignment.thickStart == 1200, msg
        if bedN >= 8:
            assert alignment.thickEnd == 4900, msg
        if bedN >= 9:
            assert alignment.itemRgb == "255,0,0", msg
        alignment = next(alignments)
        if bedN >= 5:
            assert alignment.score == 900, msg
        assert alignment.shape == (2, 4000), msg
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1], msg
        if bedN >= 6:
            assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1], msg
        else:
            assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1], msg
        assert len(alignment) == 2, msg
        assert alignment.sequences[0] is alignment.target, msg
        assert alignment.sequences[1] is alignment.query, msg
        assert alignment.target.id == "chr22", msg
        if bedN >= 4:
            assert alignment.query.id == "mRNA2", msg
        else:
            assert alignment.query.id is None, msg
        if bedN == 12:
            assert np.array_equal(
                    alignment.coordinates,
                    # fmt: off
                    np.array([[2000, 2433, 5601, 6000],
                              [ 832,  399,  399,    0]])
                    # fmt: on
                ), msg
        elif bedN >= 6:
            assert np.array_equal(
                    alignment.coordinates,
                    np.array([[2000, 6000], [4000, 0]]),
                ), msg
        else:
            assert np.array_equal(
                    alignment.coordinates,
                    np.array([[2000, 6000], [0, 4000]]),
                ), msg
        if bedN >= 7:
            assert alignment.thickStart == 2300, msg
        if bedN >= 8:
            assert alignment.thickEnd == 5960, msg
        if bedN >= 9:
            assert alignment.itemRgb == "0,255,0", msg
        with pytest.raises(StopIteration) as cm:
            next(alignments)
            raise AssertionError("More than two alignments reported")

    def test_reading(self):
        """Test parsing alignments in file formats BED3 through BED12."""
        for bedN in (3, 4, 5, 6, 7, 8, 9, 12):
            filename = "bed%d.bb" % bedN
            path = os.path.join("Blat", filename)
            alignments = Align.parse(path, "bigbed")
            msg = "bed%d" % bedN
            self.check_autosql(alignments.declaration, bedN, msg)
            self.check_alignments(alignments, bedN, msg)

    def test_writing(self):
        """Test Writing alignments in file formats BED3 through BED12."""
        for bedN in (3, 4, 5, 6, 7, 8, 9, 12):
            filename = "bed%d.bb" % bedN
            path = os.path.join("Blat", filename)
            alignments = Align.parse(path, "bigbed")
            with tempfile.TemporaryFile() as output:
                Align.write(alignments, output, "bigbed", bedN=bedN)
                output.flush()
                output.seek(0)
                alignments = Align.parse(output, "bigbed")
                msg = "bed%d" % bedN
                self.check_autosql(alignments.declaration, bedN, msg)
                self.check_alignments(alignments, bedN, msg)


class TestAlign_extended_bed(unittest.TestCase):
    # The bigBed files bigbed_extended.littleendian.bb and
    # bigbed_extended.bigendian.bb are BED9+2 files, with nine predefined BED
    # fields and 2 extra (custom) fields. It was created by running
    #
    # bedToBigBed -as=bedExample2.as -type=bed9+2 -extraIndex=name,geneSymbol bedExample2.bed hg18.chrom.sizes bigbed_extended.bb
    #
    # where bedExample2.bed contains 10 lines selected from the example BED9+2
    # file bedExample2.bed downloaded from UCSC
    # (https://genome.ucsc.edu/goldenPath/help/examples/bedExample2.bed)
    # and bedExample2.as the associated AutoSQL file, also downloaded from UCSC
    # (https://genome.ucsc.edu/goldenPath/help/examples/bedExample2.as)
    # declaring the nine predefined BED fields and the two extra fields
    #
    # The bigBed file bigbed_extended.littleendian.bb was generated on a
    # little-endian machine; the bigBed file bigbed_extended.bigendian.bb was
    # generated on a big-endian machine.

    def test_reading(self):
        """Test parsing bigbed_extended.bb."""
        path = "Blat/bigbed_extended.littleendian.bb"
        alignments = Align.parse(path, "bigbed")
        assert alignments.byteorder == "<"
        self.check_alignments(alignments)
        path = "Blat/bigbed_extended.bigendian.bb"
        alignments = Align.parse(path, "bigbed")
        assert alignments.byteorder == ">"
        self.check_alignments(alignments)

    def check_alignments(self, alignments):
        assert str(alignments.declaration) == """\
table hg18KGchr7
"UCSC Genes for chr7 with color plus GeneSymbol and SwissProtID"
(
   string  chrom;         "Reference sequence chromosome or scaffold"
   uint    chromStart;    "Start position of feature on chromosome"
   uint    chromEnd;      "End position of feature on chromosome"
   string  name;          "Name of gene"
   uint    score;         "Score"
   char[1] strand;        "+ or - for strand"
   uint    thickStart;    "Coding region start"
   uint    thickEnd;      "Coding region end"
   uint    reserved;      "Green on + strand, Red on - strand"
   string  geneSymbol;    "Gene Symbol"
   string  spID;          "SWISS-PROT protein Accession number"
)
"""
        assert len(alignments.targets) == 1
        assert alignments.targets[0].id == "chr7"
        assert len(alignments.targets[0]) == 158821424
        assert len(alignments) == 10
        alignment = next(alignments)
        assert alignment.score == 0
        assert alignment.thickStart == 60328
        assert alignment.thickEnd == 60328
        assert alignment.itemRgb == "255,0,0"
        assert alignment.shape == (2, 1241)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr7"
        assert alignment.query.id == "uc010krx.1"
        assert np.array_equal(alignment.coordinates, np.array([[60328, 61569], [1241, 0]]))
        assert alignment.annotations["geneSymbol"] == "."
        assert alignment.annotations["spID"] == "PDGFA"
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (1241 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 1241:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 1241
        alignment = next(alignments)
        assert alignment.score == 0
        assert alignment.thickStart == 506606
        assert alignment.thickEnd == 525164
        assert alignment.itemRgb == "255,0,0"
        assert alignment.shape == (2, 22585)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr7"
        assert alignment.query.id == "uc003sir.1"
        assert np.array_equal(
                alignment.coordinates, np.array([[503422, 526007], [22585, 0]])
            )
        assert alignment.annotations["geneSymbol"] == "PDGFA"
        assert alignment.annotations["spID"] == "P04085"
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (22585 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 22585:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 22585
        alignment = next(alignments)
        assert alignment.score == 0
        assert alignment.thickStart == 504726
        assert alignment.thickEnd == 525164
        assert alignment.itemRgb == "255,0,0"
        assert alignment.shape == (2, 22585)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr7"
        assert alignment.query.id == "uc003sis.1"
        assert np.array_equal(
                alignment.coordinates, np.array([[503422, 526007], [22585, 0]])
            )
        assert alignment.annotations["geneSymbol"] == "PDGFA"
        assert alignment.annotations["spID"] == "P04085-2"
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (22585 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 22585:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 22585
        alignment = next(alignments)
        assert alignment.score == 0
        assert alignment.thickStart == 507195
        assert alignment.thickEnd == 518820
        assert alignment.itemRgb == "255,0,0"
        assert alignment.shape == (2, 12690)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr7"
        assert alignment.query.id == "uc003sit.1"
        assert np.array_equal(
                alignment.coordinates, np.array([[506940, 519630], [12690, 0]])
            )
        assert alignment.annotations["geneSymbol"] == "PDGFA"
        assert alignment.annotations["spID"] == "Q32M96"
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (12690 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 12690:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 12690
        alignment = next(alignments)
        assert alignment.score == 0
        assert alignment.thickStart == 556592
        assert alignment.thickEnd == 717668
        assert alignment.itemRgb == "255,0,0"
        assert alignment.shape == (2, 162747)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr7"
        assert alignment.query.id == "uc003siu.1"
        assert np.array_equal(
                alignment.coordinates, np.array([[555912, 718659], [162747, 0]])
            )
        assert alignment.annotations["geneSymbol"] == "PRKAR1B"
        assert alignment.annotations["spID"] == "Q8N422"
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (162747 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 162747:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 162747
        alignment = next(alignments)
        assert alignment.score == 0
        assert alignment.thickStart == 556592
        assert alignment.thickEnd == 717668
        assert alignment.itemRgb == "255,0,0"
        assert alignment.shape == (2, 163357)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr7"
        assert alignment.query.id == "uc003siv.1"
        assert np.array_equal(
                alignment.coordinates, np.array([[555912, 719269], [163357, 0]])
            )
        assert alignment.annotations["geneSymbol"] == "PRKAR1B"
        assert alignment.annotations["spID"] == "Q8N422"
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (163357 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 163357:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 163357
        alignment = next(alignments)
        assert alignment.score == 0
        assert alignment.thickStart == 556592
        assert alignment.thickEnd == 717668
        assert alignment.itemRgb == "255,0,0"
        assert alignment.shape == (2, 177901)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr7"
        assert alignment.query.id == "uc003siw.1"
        assert np.array_equal(
                alignment.coordinates, np.array([[555912, 733813], [177901, 0]])
            )
        assert alignment.annotations["geneSymbol"] == "PRKAR1B"
        assert alignment.annotations["spID"] == "Q8N422"
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (177901 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 177901:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 177901
        alignment = next(alignments)
        assert alignment.score == 0
        assert alignment.thickStart == 585418
        assert alignment.thickEnd == 585418
        assert alignment.itemRgb == "255,0,0"
        assert alignment.shape == (2, 22329)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] > alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr7"
        assert alignment.query.id == "uc003six.1"
        assert np.array_equal(
                alignment.coordinates, np.array([[585418, 607747], [22329, 0]])
            )
        assert alignment.annotations["geneSymbol"] == "."
        assert alignment.annotations["spID"] == "PRKAR1B"
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (22329 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 22329:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 22329
        alignment = next(alignments)
        assert alignment.score == 0
        assert alignment.thickStart == 733217
        assert alignment.thickEnd == 791816
        assert alignment.itemRgb == "0,255,0"
        assert alignment.shape == (2, 59779)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr7"
        assert alignment.query.id == "uc003siz.2"
        assert np.array_equal(
                alignment.coordinates, np.array([[732863, 792642], [0, 59779]])
            )
        assert alignment.annotations["geneSymbol"] == "."
        assert alignment.annotations["spID"] == "DKFZp762F1415"
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (59779 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 59779:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 59779
        alignment = next(alignments)
        assert alignment.score == 0
        assert alignment.thickStart == 732883
        assert alignment.thickEnd == 791816
        assert alignment.itemRgb == "0,255,0"
        assert alignment.shape == (2, 59779)
        assert alignment.coordinates[0, 0] < alignment.coordinates[0, -1]
        assert alignment.coordinates[1, 0] < alignment.coordinates[1, -1]
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr7"
        assert alignment.query.id == "uc010krz.1"
        assert np.array_equal(
                alignment.coordinates, np.array([[732863, 792642], [0, 59779]])
            )
        assert alignment.annotations["geneSymbol"] == "HEATR2"
        assert alignment.annotations["spID"] == "Q86Y56"
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (59779 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 59779:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 59779

    def test_writing(self):
        """Test writing bigbed_extended.bb."""
        byteorder = sys.byteorder  # "little" or "big"
        path = f"Blat/bigbed_extended.{byteorder}endian.bb"
        with open(path, "rb") as stream:
            correct = stream.read()
        alignments = Align.parse(path, "bigbed")
        with open("Blat/bedExample2.as") as stream:
            autosql_data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(autosql_data)
        with tempfile.TemporaryFile() as output:
            Align.write(
                alignments,
                output,
                "bigbed",
                bedN=9,
                declaration=declaration,
                extraIndex=["name", "geneSymbol"],
            )
            output.flush()
            output.seek(0)
            data = output.read()
        assert correct == data
        alignments = Align.parse(path, "bigbed")
        targets = alignments.targets
        with tempfile.TemporaryFile() as output:
            Align.write(
                alignments,
                output,
                "bigbed",
                bedN=9,
                declaration=declaration,
                targets=targets,
                extraIndex=["name", "geneSymbol"],
            )
            output.flush()
            output.seek(0)
            data = output.read()
        assert correct == data


class TestAlign_searching(unittest.TestCase):
    path = "Blat/bigbedtest.bb"

    # The bigBed file bigbedtest.bb contains the following data:
    # chr1     10     100     name1   1       +
    # chr1     29      39     name2   2       -
    # chr1    200     300     name3   3       +
    # chr2     50      50     name4   6       +
    # chr2    100     110     name5   4       +
    # chr2    200     210     name6   5       +
    # chr2    220     220     name7   6       +
    # chr3      0       0     name8   7       -

    def check_alignments(self, alignments):
        assert str(alignments.declaration) == """\
table bed
"Browser Extensible Data"
(
   string  chrom;         "Reference sequence chromosome or scaffold"
   uint    chromStart;    "Start position in chromosome"
   uint    chromEnd;      "End position in chromosome"
   string  name;          "Name of item."
   uint    score;         "Score (0-1000)"
   char[1] strand;        "+ or - for strand"
)
"""
        assert len(alignments.targets) == 3
        assert alignments.targets[0].id == "chr1"
        assert len(alignments.targets[0]) == 1000
        assert alignments.targets[1].id == "chr2"
        assert len(alignments.targets[1]) == 2000
        assert alignments.targets[2].id == "chr3"
        assert len(alignments.targets[2]) == 1000
        assert len(alignments) == 8
        alignment = next(alignments)
        assert alignment.score == 1
        assert alignment.shape == (2, 90)
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr1"
        assert alignment.query.id == "name1"
        assert np.array_equal(alignment.coordinates, np.array([[10, 100], [0, 90]]))
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (90 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 90:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 90
        alignment = next(alignments)
        assert alignment.score == 2
        assert alignment.shape == (2, 10)
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr1"
        assert alignment.query.id == "name2"
        assert np.array_equal(alignment.coordinates, np.array([[29, 39], [10, 0]]))
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (10 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 10:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 10
        alignment = next(alignments)
        assert alignment.score == 3
        assert alignment.shape == (2, 100)
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr1"
        assert alignment.query.id == "name3"
        assert np.array_equal(alignment.coordinates, np.array([[200, 300], [0, 100]]))
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (100 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 100:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 100
        alignment = next(alignments)
        assert alignment.score == 6
        assert alignment.shape == (2, 0)
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr2"
        assert alignment.query.id == "name4"
        assert np.array_equal(alignment.coordinates, np.array([[50, 50], [0, 0]]))
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (0 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 0:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 0
        assert counts.identities == 0
        assert counts.mismatches == 0
        alignment = next(alignments)
        assert alignment.score == 4
        assert alignment.shape == (2, 10)
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr2"
        assert alignment.query.id == "name5"
        assert np.array_equal(alignment.coordinates, np.array([[100, 110], [0, 10]]))
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (10 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 10:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 10
        alignment = next(alignments)
        assert alignment.score == 5
        assert alignment.shape == (2, 10)
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr2"
        assert alignment.query.id == "name6"
        assert np.array_equal(alignment.coordinates, np.array([[200, 210], [0, 10]]))
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (10 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 10:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 10
        alignment = next(alignments)
        assert alignment.score == 6
        assert alignment.shape == (2, 0)
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr2"
        assert alignment.query.id == "name7"
        assert np.array_equal(alignment.coordinates, np.array([[220, 220], [0, 0]]))
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (0 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 0:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 0
        assert counts.identities == 0
        assert counts.mismatches == 0
        alignment = next(alignments)
        assert alignment.score == 7
        assert alignment.shape == (2, 0)
        assert len(alignment) == 2
        assert alignment.sequences[0] is alignment.target
        assert alignment.sequences[1] is alignment.query
        assert alignment.target.id == "chr3"
        assert alignment.query.id == "name8"
        assert np.array_equal(alignment.coordinates, np.array([[0, 0], [0, 0]]))
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (0 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 0:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 0
        assert counts.left_deletions == 0
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 0
        assert counts.internal_deletions == 0
        assert counts.left_gaps == 0
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 0
        assert counts.insertions == 0
        assert counts.deletions == 0
        assert counts.gaps == 0
        assert counts.aligned == 0
        assert counts.identities == 0
        assert counts.mismatches == 0

    def test_reading(self):
        """Test reading bigbedtest.bb."""
        alignments = Align.parse(self.path, "bigbed")
        self.check_alignments(alignments)

    def test_writing(self):
        """Test writing bigbedtest.bb."""
        alignments = Align.parse(self.path, "bigbed")
        with tempfile.TemporaryFile() as output:
            Align.write(alignments, output, "bigbed", bedN=6)
            output.flush()
            output.seek(0)
            alignments = Align.parse(output, "bigbed")
            self.check_alignments(alignments)

    def test_search_chromosome(self):
        alignments = Align.parse(self.path, "bigbed")
        assert str(alignments.declaration) == """\
table bed
"Browser Extensible Data"
(
   string  chrom;         "Reference sequence chromosome or scaffold"
   uint    chromStart;    "Start position in chromosome"
   uint    chromEnd;      "End position in chromosome"
   string  name;          "Name of item."
   uint    score;         "Score (0-1000)"
   char[1] strand;        "+ or - for strand"
)
"""
        selected_alignments = alignments.search("chr2")
        names = [alignment.query.id for alignment in selected_alignments]
        assert names == ["name4", "name5", "name6", "name7"]

    def test_search_region(self):
        alignments = Align.parse(self.path, "bigbed")
        selected_alignments = alignments.search("chr2", 105, 1000)
        names = [alignment.query.id for alignment in selected_alignments]
        assert names == ["name5", "name6", "name7"]
        selected_alignments = alignments.search("chr2", 110, 1000)
        names = [alignment.query.id for alignment in selected_alignments]
        assert names == ["name6", "name7"]
        selected_alignments = alignments.search("chr2", 40, 50)
        names = [alignment.query.id for alignment in selected_alignments]
        assert names == ["name4"]
        selected_alignments = alignments.search("chr2", 50, 50)
        names = [alignment.query.id for alignment in selected_alignments]
        assert names == ["name4"]
        selected_alignments = alignments.search("chr2", 50, 200)
        names = [alignment.query.id for alignment in selected_alignments]
        assert names == ["name4", "name5"]
        selected_alignments = alignments.search("chr2", 200, 220)
        names = [alignment.query.id for alignment in selected_alignments]
        assert names == ["name6", "name7"]
        selected_alignments = alignments.search("chr2", 220, 220)
        names = [alignment.query.id for alignment in selected_alignments]
        assert names == ["name7"]

    def test_search_position(self):
        alignments = Align.parse(self.path, "bigbed")
        selected_alignments = alignments.search("chr1", 250)
        names = [alignment.query.id for alignment in selected_alignments]
        assert names == ["name3"]

    def test_three_iterators(self):
        """Create three iterators and use them concurrently."""
        alignments1 = Align.parse(self.path, "bigbed")
        alignments2 = alignments1.search("chr2")
        alignments3 = alignments1.search("chr2", 110, 1000)
        alignment1 = next(alignments1)
        assert alignment1.query.id == "name1"
        alignment1 = next(alignments1)
        assert alignment1.query.id == "name2"
        alignment2 = next(alignments2)
        assert alignment2.query.id == "name4"
        alignment2 = next(alignments2)
        assert alignment2.query.id == "name5"
        alignment2 = next(alignments2)
        assert alignment2.query.id == "name6"
        alignment3 = next(alignments3)
        assert alignment3.query.id == "name6"
        alignment3 = next(alignments3)
        assert alignment3.query.id == "name7"
        alignment1 = next(alignments1)
        assert alignment1.query.id == "name3"
        alignment1 = next(alignments1)
        assert alignment1.query.id == "name4"
        alignment1 = next(alignments1)
        assert alignment1.query.id == "name5"
        alignment2 = next(alignments2)
        assert alignment2.query.id == "name7"
        with pytest.raises(StopIteration):
            next(alignments2)
        alignment1 = next(alignments1)
        assert alignment1.query.id == "name6"
        alignment1 = next(alignments1)
        assert alignment1.query.id == "name7"
        with pytest.raises(StopIteration):
            next(alignments3)
        alignment1 = next(alignments1)
        assert alignment1.query.id == "name8"
        with pytest.raises(StopIteration):
            next(alignments1)


class BinaryTestBaseClass(unittest.TestCase):
    def assertBinaryEqual(self, file1, file2):
        blocksize = 1024
        n = 0
        while True:
            data1 = file1.read(blocksize)
            data2 = file2.read(blocksize)
            if data1 == b"" and data2 == b"":
                return
            n1 = len(data1)
            if data1 == data2:
                n += n1
                continue
            n2 = len(data2)
            if n1 < n2:
                return self.fail(f"unequal file sizes: {n1} bytes vs >= {n2} bytes")
            if n1 > n2:
                return self.fail(f"unequal file sizes: >= {n1} bytes vs {n2} bytes")
            for i, (c1, c2) in enumerate(zip(data1, data2)):
                if c1 != c2:
                    return self.fail(f"bytes at position {n + i} differ: {c1} vs {c2}")


@unittest.skipUnless(big is True, "big file; use --big to run")
class TestAlign_big(BinaryTestBaseClass):
    # BED files were downloaded from the UCSC table browser:
    #
    # ucsc.bed contains the GENCODE V43 Basic gene annotations for human genome
    # assembly hg38.
    #
    # anoGam3.bed contains the AUGUSTUS gene annotations for genome assembly
    # anoGam3 of Anopheles gambiae (African malaria mosquito).
    #
    # ailMel1.bed contains the NCBI RefSeq All gene annotations for genome
    # assembly ailMel1 of Ailuropoda melanoleuca (giant panda).
    #
    # bisBis1.bed contains the AUGUSTUS gene annotations for genome assembly
    # bisBis1 of Bison bison bison (American bison).
    #
    # bigBed files were generated using bedToBigBed v. 2.9 found in
    # jksrc.v445.zip of the 'kent' source tree provided by UCSC
    # (http://hgdownload.cse.ucsc.edu/admin/).

    def test_a_compressed(self):
        # grep -E -v 'fix|alt' Blat/ucsc.bed > ucsc.clean.bed
        # sort -k1,1 -k2,2n ucsc.clean.bed -o ucsc.clean.bed
        # bedToBigBed -as=Blat/bed12.as ucsc.clean.bed Align/hg38.chrom.sizes ucsc.bb
        with open("Blat/bed12.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_b_uncompressed(self):
        # grep -E -v 'fix|alt' Blat/ucsc.bed > ucsc.clean.bed
        # sort -k1,1 -k2,2n ucsc.clean.bed -o ucsc.clean.bed
        # bedToBigBed -as=Blat/bed12.as -unc ucsc.clean.bed Align/hg38.chrom.sizes ucsc.unc.bb
        with open("Blat/bed12.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        bigBedFileName = "ucsc.unc.bb"
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                declaration=declaration,
                compress=False,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_c_bed3(self):
        # grep -E -v 'fix|alt' Blat/ucsc.bed > ucsc.clean.bed
        # sort -k1,1 -k2,2n ucsc.clean.bed -o ucsc.clean.bed
        # cut -f 1-3 ucsc.clean.bed > ucsc.bed3.bed
        # bedToBigBed -as=Blat/bed3.as -type=bed3 ucsc.bed3.bed Align/hg38.chrom.sizes ucsc.bed3.bb
        with open("Blat/bed3.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        bigBedFileName = "ucsc.bed3.bb"
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                bedN=3,
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_d_bed4(self):
        # grep -E -v 'fix|alt' Blat/ucsc.bed > ucsc.clean.bed
        # sort -k1,1 -k2,2n ucsc.clean.bed -o ucsc.clean.bed
        # cut -f 1-4 ucsc.clean.bed > ucsc.bed4.bed
        # bedToBigBed -as=Blat/bed4.as -type=bed4 ucsc.bed4.bed Align/hg38.chrom.sizes ucsc.bed4.bb
        with open("Blat/bed4.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        bigBedFileName = "ucsc.bed4.bb"
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                bedN=4,
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_e_bed5(self):
        # grep -E -v 'fix|alt' Blat/ucsc.bed > ucsc.clean.bed
        # sort -k1,1 -k2,2n ucsc.clean.bed -o ucsc.clean.bed
        # cut -f 1-5 ucsc.clean.bed > ucsc.bed5.bed
        # bedToBigBed -as=Blat/bed5.as -type=bed5 ucsc.bed5.bed Align/hg38.chrom.sizes ucsc.bed5.bb
        with open("Blat/bed5.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        bigBedFileName = "ucsc.bed5.bb"
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                bedN=5,
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_f_bed6(self):
        # grep -E -v 'fix|alt' Blat/ucsc.bed > ucsc.clean.bed
        # sort -k1,1 -k2,2n ucsc.clean.bed -o ucsc.clean.bed
        # cut -f 1-6 ucsc.clean.bed > ucsc.bed6.bed
        # bedToBigBed -as=Blat/bed6.as -type=bed6 ucsc.bed6.bed Align/hg38.chrom.sizes ucsc.bed6.bb
        with open("Blat/bed6.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        bigBedFileName = "ucsc.bed6.bb"
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                bedN=6,
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_g_bed7(self):
        # grep -E -v 'fix|alt' Blat/ucsc.bed > ucsc.clean.bed
        # sort -k1,1 -k2,2n ucsc.clean.bed -o ucsc.clean.bed
        # cut -f 1-7 ucsc.clean.bed > ucsc.bed7.bed
        # bedToBigBed -as=Blat/bed7.as -type=bed7 ucsc.bed7.bed Align/hg38.chrom.sizes ucsc.bed7.bb
        with open("Blat/bed7.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        bigBedFileName = "ucsc.bed7.bb"
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                bedN=7,
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_h_bed8(self):
        # grep -E -v 'fix|alt' Blat/ucsc.bed > ucsc.clean.bed
        # sort -k1,1 -k2,2n ucsc.clean.bed -o ucsc.clean.bed
        # cut -f 1-8 ucsc.clean.bed > ucsc.bed8.bed
        # bedToBigBed -as=Blat/bed8.as -type=bed8 ucsc.bed8.bed Align/hg38.chrom.sizes ucsc.bed8.bb
        with open("Blat/bed8.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        bigBedFileName = "ucsc.bed8.bb"
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                bedN=8,
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_i_bed9(self):
        # grep -E -v 'fix|alt' Blat/ucsc.bed > ucsc.clean.bed
        # sort -k1,1 -k2,2n ucsc.clean.bed -o ucsc.clean.bed
        # cut -f 1-9 ucsc.clean.bed > ucsc.bed9.bed
        # bedToBigBed -as=Blat/bed9.as -type=bed9 ucsc.bed9.bed Align/hg38.chrom.sizes ucsc.bed9.bb
        with open("Blat/bed9.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        bigBedFileName = "ucsc.bed9.bb"
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                bedN=9,
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_j_extraindex(self):
        # If the names in the BED file are not unique, the binary contents of
        # the generated bigBed file will depend on the exact implementation of
        # the quicksort algorithm in the Standard Library of C (used by
        # bedToBigBed) and the quicksort algorithm in numpy (used by Biopython).
        # To prevent spurious errors when comparing bigBed files created by
        # bedToBigBed and by Biopython, we first remove lines with duplicated
        # names from the BED file.
        # grep -E -v 'fix|alt' Blat/ucsc.bed > ucsc.clean.bed
        # sort -u -k4,4 ucsc.clean.bed | sort -k1,1 -k2,2n > ucsc.unique.bed
        # bedToBigBed -as=Blat/bed12.as -extraIndex=name ucsc.unique.bed Align/hg38.chrom.sizes ucsc.indexed.bb
        with open("Blat/bed12.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.indexed.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                declaration=declaration,
                extraIndex=["name"],
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_k_anogam(self):
        # sort -k1,1 -k2,2n Blat/anoGam3.bed -o anoGam3.bed
        # bedToBigBed -as=Blat/bed12.as anoGam3.bed Blat/anoGam3.chrom.sizes anoGam3.bb
        with open("Blat/bed12.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "anoGam3.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_l_ailmel(self):
        # sort -k1,1 -k2,2n Blat/ailMel1.bed -o ailMel1.bed
        # bedToBigBed -as=Blat/bed12.as ailMel1.bed Blat/ailMel1.chrom.sizes ailMel1.bb
        with open("Blat/bed12.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ailMel1.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_m_bisbis(self):
        # sort -k1,1 -k2,2n bisBis1.bed -o bisBis1.bed
        # bedToBigBed -as=Blat/bed12.as bisBis1.bed Blat/bisBis1.chrom.sizes bisBis1.bb
        with open("Blat/bed12.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "bisBis1.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)


class TestDeclarations(unittest.TestCase):
    def test_declarations(self):
        for length in (3, 4, 5, 6, 7, 8, 9, 12):
            filename = "bed%d.as" % length
            path = os.path.join("Blat", filename)
            with open(path) as stream:
                data = stream.read()
            declaration = bigbed.AutoSQLTable.from_string(data)
            assert declaration.name == "bed", filename
            assert declaration.comment == "Browser Extensible Data", filename
            assert len(declaration) == length, filename
            field = declaration[0]
            assert field.as_type == "string", filename
            assert field.name == "chrom", filename
            assert field.comment == "Reference sequence chromosome or scaffold", filename
            field = declaration[1]
            assert field.as_type == "uint", filename
            assert field.name == "chromStart", filename
            assert field.comment == "Start position in chromosome", filename
            field = declaration[2]
            assert field.as_type == "uint", filename
            assert field.name == "chromEnd", filename
            assert field.comment == "End position in chromosome", filename
            if length == 3:
                return
            field = declaration[3]
            assert field.as_type == "string", filename
            assert field.name == "name", filename
            assert field.comment == "Name of item.", filename
            if length == 4:
                return
            field = declaration[4]
            assert field.as_type == "uint", filename
            assert field.name == "score", filename
            assert field.comment == "Score (0-1000)", filename
            if length == 5:
                return
            field = declaration[5]
            assert field.as_type == "char[1]", filename
            assert field.name == "strand", filename
            assert field.comment == "+ or - for strand", filename
            if length == 6:
                return
            field = declaration[6]
            assert field.as_type == "uint", filename
            assert field.name == "thickStart", filename
            assert field.comment == "Start of where display should be thick (start codon)", filename
            if length == 7:
                return
            field = declaration[7]
            assert field.as_type == "uint", filename
            assert field.name == "thickEnd", filename
            assert field.comment == "End of where display should be thick (stop codon)", filename
            if length == 8:
                return
            field = declaration[8]
            assert field.as_type == "uint", filename
            assert field.name == "reserved", filename
            assert field.comment == "Used as itemRgb as of 2004-11-22", filename
            if length == 9:
                return
            field = declaration[9]
            assert field.as_type == "int", filename
            assert field.name == "blockCount", filename
            assert field.comment == "Number of blocks", filename
            field = declaration[10]
            assert field.as_type == "int[blockCount]", filename
            assert field.name == "blockSizes", filename
            assert field.comment == "Comma separated list of block sizes", filename
            field = declaration[11]
            assert field.as_type == "int[blockCount]", filename
            assert field.name == "chromStarts", filename
            assert field.comment == "Start positions relative to chromStart", filename


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
