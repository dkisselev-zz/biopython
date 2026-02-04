# Copyright 2006-2014 by Peter Cock.  All rights reserved.
# Revisions copyright 2011 Brandon Invergo. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.Align.mauve module."""
import os
import unittest
import pytest
from io import StringIO

from Bio import Align
from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Seq import Seq

np = pytest.importorskip("numpy")
class TestCombinedFile(unittest.TestCase):
    # Generate the output file combined.xmfa by running
    # progressiveMauve combined.fa --output=combined.xmfa

    filename = "combined.fa"
    path = os.path.join("Mauve", filename)
    records = SeqIO.parse(path, "fasta")
    sequences = {str(index): record.seq for index, record in enumerate(records)}
    del filename
    del path
    del records

    def test_parse(self):
        path = os.path.join("Mauve", "combined.xmfa")
        with open(path) as stream:
            alignments = Align.parse(stream, "mauve")
            self.check_alignments(alignments)
            alignments = iter(alignments)
            self.check_alignments(alignments)
        with Align.parse(path, "mauve") as alignments:
            self.check_alignments(alignments)
        with pytest.raises(AttributeError):
            alignments._stream
        with Align.parse(path, "mauve") as alignments:
            pass
        with pytest.raises(AttributeError):
            alignments._stream

    def check_alignments(self, alignments):
        saved_alignments = []
        metadata = alignments.metadata
        assert len(metadata) == 3
        assert metadata["FormatVersion"] == "Mauve1"
        assert metadata["File"] == "combined.fa"
        assert metadata["BackboneFile"] == "combined.xmfa.bbcols"
        identifiers = alignments.identifiers
        alignment = next(alignments)
        saved_alignments.append(alignment)
        assert len(alignment) == 3
        assert len(alignment.sequences) == 3
        assert alignment.sequences[0].id == "0"
        assert repr(alignment.sequences[0].seq) == "Seq({1: 'AAAAGGAAAGTACGGCCCGGCCACTCCGGGTGTGTGCTAGGAGGGCTT'}, length=49)"
        sequence = self.sequences[alignment.sequences[0].id]
        start = len(sequence) - alignment.coordinates[0, 0]
        end = len(sequence) - alignment.coordinates[0, -1]
        assert start == 1
        assert end == 49
        assert alignment.sequences[0].seq[start:end] == sequence[start:end]
        assert alignment.sequences[1].id == "1"
        assert alignment.sequences[1].seq == ""
        start = alignment.coordinates[1, 0]
        end = alignment.coordinates[1, -1]
        assert start == 0
        assert end == 0
        assert alignment.sequences[2].id == "2"
        assert repr(alignment.sequences[2].seq) == "Seq({1: 'AAGCCCTGCGCGCTCAGCCGGAGTGTCCCGGGCCCTGCTTTCCTTTT'}, length=48)"
        start = alignment.coordinates[2, 0]
        end = alignment.coordinates[2, -1]
        assert start == 1
        assert end == 48
        sequence = self.sequences[alignment.sequences[2].id][start:end]
        assert alignment.sequences[2].seq[start:end] == sequence
        assert alignment[0] == "AAGCCCTCCTAGCACACACCCGGAGTGG-CCGGGCCGTACTTTCCTTTT"
        assert alignment[1] == "-------------------------------------------------"
        assert alignment[2] == "AAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGCTTTCCTTTT"
        assert np.array_equal(
                alignment.coordinates,
                np.array(
                    [
                        [49, 40, 38, 21, 21, 1],
                        [0, 0, 0, 0, 0, 0],
                        [1, 10, 10, 27, 28, 48],
                    ]
                ),
            )
        assert str(alignment) == """\
0                49 AAGCCCTCCTAGCACACACCCGGAGTGG-CCGGGCCGTACTTTCCTTTT  1
1                 0 -------------------------------------------------  0
2                 1 AAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGCTTTCCTTTT 48
"""
        assert alignment.format("mauve", metadata, identifiers) == """\
> 1:2-49 - combined.fa
AAGCCCTCCTAGCACACACCCGGAGTGG-CCGGGCCGTACTTTCCTTTT
> 2:0-0 + combined.fa
-------------------------------------------------
> 3:2-48 + combined.fa
AAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGCTTTCCTTTT
=
"""
        assert np.array_equal(
                np.array(alignment, "U"),
                # fmt: off
np.array([['A', 'A', 'G', 'C', 'C', 'C', 'T', 'C', 'C', 'T', 'A', 'G', 'C',
           'A', 'C', 'A', 'C', 'A', 'C', 'C', 'C', 'G', 'G', 'A', 'G', 'T',
           'G', 'G', '-', 'C', 'C', 'G', 'G', 'G', 'C', 'C', 'G', 'T', 'A',
           'C', 'T', 'T', 'T', 'C', 'C', 'T', 'T', 'T', 'T'],
          ['-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
           '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
           '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
           '-', '-', '-', '-', '-', '-', '-', '-', '-', '-'],
          ['A', 'A', 'G', 'C', 'C', 'C', 'T', 'G', 'C', '-', '-', 'G', 'C',
           'G', 'C', 'T', 'C', 'A', 'G', 'C', 'C', 'G', 'G', 'A', 'G', 'T',
           'G', 'T', 'C', 'C', 'C', 'G', 'G', 'G', 'C', 'C', 'C', 'T', 'G',
           'C', 'T', 'T', 'T', 'C', 'C', 'T', 'T', 'T', 'T']], dtype='U')
                # fmt: on
            )
        counts = alignment.counts()
        assert (repr(counts) == "<AlignmentCounts object (46 aligned letters; 39 identities; 7 mismatches; 98 gaps) at 0x%x>"
            % id(counts))
        assert str(counts) == """\
AlignmentCounts object with
    aligned = 46:
        identities = 39,
        mismatches = 7.
    gaps = 98:
        left_gaps = 95:
            left_insertions = 47:
                open_left_insertions = 1,
                extend_left_insertions = 46;
            left_deletions = 48:
                open_left_deletions = 1,
                extend_left_deletions = 47;
        internal_gaps = 3:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
        assert counts.left_insertions == 47
        assert counts.left_deletions == 48
        assert counts.right_insertions == 0
        assert counts.right_deletions == 0
        assert counts.internal_insertions == 1
        assert counts.internal_deletions == 2
        assert counts.left_gaps == 95
        assert counts.right_gaps == 0
        assert counts.internal_gaps == 3
        assert counts.insertions == 48
        assert counts.deletions == 50
        assert counts.gaps == 98
        assert counts.aligned == 46
        assert counts.identities == 39
        assert counts.mismatches == 7
        alignment = next(alignments)
        saved_alignments.append(alignment)
        assert len(alignment) == 1
        assert len(alignment.sequences) == 1
        assert alignment.sequences[0].id == "0"
        assert alignment.sequences[0].seq == "G"
        sequence = self.sequences[alignment.sequences[0].id]
        start = alignment.coordinates[0, 0]
        end = alignment.coordinates[0, -1]
        assert alignment.sequences[0].seq[start:end] == sequence[start:end]
        assert alignment[0] == "G"
        assert np.array_equal(alignment.coordinates, np.array([[0, 1]]))
        assert str(alignment) == """\
0                 0 G 1
"""
        assert alignment.format("mauve", metadata, identifiers) == """\
> 1:1-1 + combined.fa
G
=
"""
        assert np.array_equal(np.array(alignment, "U"), np.array([["G"]], dtype="U"))
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
        saved_alignments.append(alignment)
        assert len(alignment) == 1
        assert len(alignment.sequences) == 1
        assert alignment.sequences[0].id == "0"
        assert repr(alignment.sequences[0].seq) == "Seq({49: 'A'}, length=50)"
        sequence = self.sequences[alignment.sequences[0].id]
        start = alignment.coordinates[0, 0]
        end = alignment.coordinates[0, -1]
        assert alignment.sequences[0].seq[start:end] == sequence[start:end]
        assert alignment[0] == "A"
        assert np.array_equal(alignment.coordinates, np.array([[49, 50]]))
        assert str(alignment) == """\
0                49 A 50
"""
        assert alignment.format("mauve", metadata, identifiers) == """\
> 1:50-50 + combined.fa
A
=
"""
        assert np.array_equal(np.array(alignment, "U"), np.array([["A"]], dtype="U"))
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
        saved_alignments.append(alignment)
        assert len(alignment) == 1
        assert len(alignment.sequences) == 1
        assert alignment.sequences[0].id == "1"
        assert alignment.sequences[0].seq == "GAAGAGGAAAAGTAGATCCCTGGCGTCCGGAGCTGGGACGT"
        sequence = self.sequences[alignment.sequences[0].id]
        start = alignment.coordinates[0, 0]
        end = alignment.coordinates[0, -1]
        assert alignment.sequences[0].seq[start:end] == sequence[start:end]
        assert alignment[0] == "GAAGAGGAAAAGTAGATCCCTGGCGTCCGGAGCTGGGACGT"
        assert np.array_equal(alignment.coordinates, np.array([[0, 41]]))
        assert str(alignment) == """\
1                 0 GAAGAGGAAAAGTAGATCCCTGGCGTCCGGAGCTGGGACGT 41
"""
        assert alignment.format("mauve", metadata, identifiers) == """\
> 2:1-41 + combined.fa
GAAGAGGAAAAGTAGATCCCTGGCGTCCGGAGCTGGGACGT
=
"""
        assert np.array_equal(
                np.array(alignment, "U"),
                # fmt: off
np.array([['G', 'A', 'A', 'G', 'A', 'G', 'G', 'A', 'A', 'A', 'A', 'G', 'T',
           'A', 'G', 'A', 'T', 'C', 'C', 'C', 'T', 'G', 'G', 'C', 'G', 'T',
           'C', 'C', 'G', 'G', 'A', 'G', 'C', 'T', 'G', 'G', 'G', 'A', 'C',
           'G', 'T']], dtype='U')
                # fmt: on
            )
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
        saved_alignments.append(alignment)
        assert len(alignment) == 1
        assert len(alignment.sequences) == 1
        assert alignment.sequences[0].id == "2"
        assert alignment.sequences[0].seq == "C"
        sequence = self.sequences[alignment.sequences[0].id]
        start = alignment.coordinates[0, 0]
        end = alignment.coordinates[0, -1]
        assert alignment.sequences[0].seq[start:end] == sequence[start:end]
        assert alignment[0] == "C"
        assert np.array_equal(alignment.coordinates, np.array([[0, 1]]))
        assert str(alignment) == """\
2                 0 C 1
"""
        assert alignment.format("mauve", metadata, identifiers) == """\
> 3:1-1 + combined.fa
C
=
"""
        assert np.array_equal(np.array(alignment, "U"), np.array([["C"]], dtype="U"))
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
        saved_alignments.append(alignment)
        assert len(alignment) == 1
        assert len(alignment.sequences) == 1
        assert alignment.sequences[0].id == "2"
        assert repr(alignment.sequences[0].seq) == "Seq({48: 'C'}, length=49)"
        sequence = self.sequences[alignment.sequences[0].id]
        start = alignment.coordinates[0, 0]
        end = alignment.coordinates[0, -1]
        assert alignment.sequences[0].seq[start:end] == sequence[start:end]
        assert alignment[0] == "C"
        assert np.array_equal(alignment.coordinates, np.array([[48, 49]]))
        assert str(alignment) == """\
2                48 C 49
"""
        assert alignment.format("mauve", metadata, identifiers) == """\
> 3:49-49 + combined.fa
C
=
"""
        assert np.array_equal(np.array(alignment, "U"), np.array([["C"]], dtype="U"))
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
        with pytest.raises(StopIteration):
            next(alignments)
        # As each nucleotide in each sequence is stored exactly once in an XMFA
        # file, we can reconstitute the full sequences:
        assert len(saved_alignments) == 6
        maxindex = -1
        for alignment in saved_alignments:
            for record in alignment.sequences:
                index = int(record.id)
                if index > maxindex:
                    maxindex = index
        n = maxindex + 1
        assert n == 3
        lengths = [0] * n
        for alignment in saved_alignments:
            for record in alignment.sequences:
                index = int(record.id)
                length = len(record.seq)
                if length > lengths[index]:
                    lengths[index] = length
        assert lengths[0] == 50
        assert lengths[1] == 41
        assert lengths[2] == 49
        sequences = [None] * 3
        for index, length in enumerate(lengths):
            sequences[index] = MutableSeq("N" * length)
        # Now fill up the sequences:
        for alignment in saved_alignments:
            for row, record in zip(alignment.coordinates, alignment.sequences):
                index = int(record.id)
                start = row[0]
                end = row[-1]
                if start > end:
                    start, end = end, start
                sequences[index][start:end] = record.seq[start:end]
        # Confirm that the fully defined sequences agree with the Fasta file:
        for index, sequence in enumerate(sequences):
            sequences[index] = Seq(sequence)
            key = str(index)
            assert sequences[index] == self.sequences[key]
        # Make sure we can replace the partially defined sequences by these
        # fully defined sequences, and get the same alignment:
        alignment = saved_alignments[0]
        for record in alignment.sequences:
            index = int(record.id)
            record.seq = sequences[index]
            assert alignment[0] == "AAGCCCTCCTAGCACACACCCGGAGTGG-CCGGGCCGTACTTTCCTTTT"
            assert alignment[1] == "-------------------------------------------------"
            assert alignment[2] == "AAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGCTTTCCTTTT"
        alignment = saved_alignments[1]
        for record in alignment.sequences:
            index = int(record.id)
            record.seq = sequences[index]
            assert alignment[0] == "G"
        alignment = saved_alignments[2]
        for record in alignment.sequences:
            index = int(record.id)
            record.seq = sequences[index]
            assert alignment[0] == "A"
        alignment = saved_alignments[3]
        for record in alignment.sequences:
            index = int(record.id)
            record.seq = sequences[index]
            assert alignment[0] == "GAAGAGGAAAAGTAGATCCCTGGCGTCCGGAGCTGGGACGT"
        alignment = saved_alignments[4]
        for record in alignment.sequences:
            index = int(record.id)
            record.seq = sequences[index]
            assert alignment[0] == "C"
        alignment = saved_alignments[5]
        for record in alignment.sequences:
            index = int(record.id)
            record.seq = sequences[index]
            assert alignment[0] == "C"

    def test_write_read(self):
        path = os.path.join("Mauve", "combined.xmfa")
        with open(path) as stream:
            data = stream.read()

        stream = StringIO()
        stream.write(data)
        stream.seek(0)
        alignments = Align.parse(stream, "mauve")
        output = StringIO()
        n = Align.write(alignments, output, "mauve")
        assert n == 6
        output.seek(0)
        assert output.read() == data


class TestSeparateFiles(unittest.TestCase):
    # Generate the output file separate.xmfa by running
    # progressiveMauve --solid-seeds equCab1.fa canFam2.fa mm9.fa --output=separate.xmfa

    sequences = {}
    for species in ("equCab1", "canFam2", "mm9"):
        filename = f"{species}.fa"
        path = os.path.join("Mauve", filename)
        record = SeqIO.read(path, "fasta")
        sequences[filename] = record.seq
        del filename
        del path
        del record

    def test_parse(self):
        path = os.path.join("Mauve", "separate.xmfa")
        saved_alignments = []
        with open(path) as stream:
            alignments = Align.parse(stream, "mauve")
            metadata = alignments.metadata
            assert len(metadata) == 2
            assert metadata["FormatVersion"] == "Mauve1"
            assert metadata["BackboneFile"] == "separate.xmfa.bbcols"
            identifiers = alignments.identifiers
            alignment = next(alignments)
            saved_alignments.append(alignment)
            assert len(alignment) == 3
            assert len(alignment.sequences) == 3
            assert alignment.sequences[0].id == "equCab1.fa"
            assert alignment.sequences[0].seq == Seq("GAAAAGGAAAGTACGGCCCGGCCACTCCGGGTGTGTGCTAGGAGGGCTTA")
            start = alignment.coordinates[0, 0]
            end = alignment.coordinates[0, -1]
            assert start == 50
            assert end == 0
            assert alignment.sequences[1].id == "canFam2.fa"
            assert alignment.sequences[1].seq == Seq("CAAGCCCTGCGCGCTCAGCCGGAGTGTCCCGGGCCCTGCTTTCCTTTTC")
            start = alignment.coordinates[1, 0]
            end = alignment.coordinates[1, -1]
            sequence = self.sequences[alignment.sequences[1].id]
            assert start == 0
            assert end == 49
            assert alignment.sequences[1].seq[start:end] == sequence[start:end]
            assert alignment.sequences[2].id == "mm9.fa"

            sequence = alignment.sequences[2].seq
            start = len(sequence) - alignment.coordinates[2, 0]
            end = len(sequence) - alignment.coordinates[2, -1]
            assert start == 0
            assert end == 19
            sequence = self.sequences[alignment.sequences[2].id][start:end]
            assert alignment.sequences[2].seq[start:end] == sequence
            assert alignment[0] == "TAAGCCCTCCTAGCACACACCCGGAGTGGCC-GGGCCGTAC-TTTCCTTTTC"
            assert alignment[1] == "CAAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGC-TTTCCTTTTC"
            assert alignment[2] == "---------------------------------GGATCTACTTTTCCTCTTC"
            assert str(alignment) == """\
equCab1.f        50 TAAGCCCTCCTAGCACACACCCGGAGTGGCC-GGGCCGTAC-TTTCCTTTTC  0
canFam2.f         0 CAAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGC-TTTCCTTTTC 49
mm9.fa           19 ---------------------------------GGATCTACTTTTCCTCTTC  0
"""
            assert np.array_equal(
                    alignment.coordinates,
                    # fmt: off
                    np.array([[50, 40, 38, 19, 19, 18, 10, 10,  0],
                              [ 0, 10, 10, 29, 30, 31, 39, 39, 49],
                              [19, 19, 19, 19, 19, 19, 11, 10,  0]]),
                    # fmt: on
                )
            assert alignment.format("mauve", metadata, identifiers) == """\
> 1:1-50 - equCab1.fa
TAAGCCCTCCTAGCACACACCCGGAGTGGCC-GGGCCGTAC-TTTCCTTTTC
> 2:1-49 + canFam2.fa
CAAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGC-TTTCCTTTTC
> 3:1-19 - mm9.fa
---------------------------------GGATCTACTTTTCCTCTTC
=
"""
            assert np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['T', 'A', 'A', 'G', 'C', 'C', 'C', 'T', 'C', 'C', 'T', 'A', 'G',
           'C', 'A', 'C', 'A', 'C', 'A', 'C', 'C', 'C', 'G', 'G', 'A', 'G',
           'T', 'G', 'G', 'C', 'C', '-', 'G', 'G', 'G', 'C', 'C', 'G', 'T',
           'A', 'C', '-', 'T', 'T', 'T', 'C', 'C', 'T', 'T', 'T', 'T', 'C'],
          ['C', 'A', 'A', 'G', 'C', 'C', 'C', 'T', 'G', 'C', '-', '-', 'G',
           'C', 'G', 'C', 'T', 'C', 'A', 'G', 'C', 'C', 'G', 'G', 'A', 'G',
           'T', 'G', 'T', 'C', 'C', 'C', 'G', 'G', 'G', 'C', 'C', 'C', 'T',
           'G', 'C', '-', 'T', 'T', 'T', 'C', 'C', 'T', 'T', 'T', 'T', 'C'],
          ['-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
           '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
           '-', '-', '-', '-', '-', '-', '-', 'G', 'G', 'A', 'T', 'C', 'T',
           'A', 'C', 'T', 'T', 'T', 'T', 'C', 'C', 'T', 'C', 'T', 'T', 'C']],
         dtype='U')
                    # fmt: on
                )
            counts = alignment.counts()
            assert (repr(counts) == "<AlignmentCounts object (84 aligned letters; 68 identities; 16 mismatches; 68 gaps) at 0x%x>"
                % id(counts))
            assert str(counts) == """\
AlignmentCounts object with
    aligned = 84:
        identities = 68,
        mismatches = 16.
    gaps = 68:
        left_gaps = 63:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 63:
                open_left_deletions = 2,
                extend_left_deletions = 61;
        internal_gaps = 5:
            internal_insertions = 3:
                open_internal_insertions = 3,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
"""
            assert counts.left_insertions == 0
            assert counts.left_deletions == 63
            assert counts.right_insertions == 0
            assert counts.right_deletions == 0
            assert counts.internal_insertions == 3
            assert counts.internal_deletions == 2
            assert counts.left_gaps == 63
            assert counts.right_gaps == 0
            assert counts.internal_gaps == 5
            assert counts.insertions == 3
            assert counts.deletions == 65
            assert counts.gaps == 68
            assert counts.aligned == 84
            assert counts.identities == 68
            assert counts.mismatches == 16
            alignment = next(alignments)
            saved_alignments.append(alignment)
            assert len(alignment) == 1
            assert len(alignment.sequences) == 1
            assert alignment.sequences[0].id == "mm9.fa"
            assert repr(alignment.sequences[0].seq) == "Seq({19: 'CTGGCGTCCGGAGCTGGGACGT'}, length=41)"
            sequence = self.sequences[alignment.sequences[0].id]
            start = alignment.coordinates[0, 0]
            end = alignment.coordinates[0, -1]
            assert alignment.sequences[0].seq[start:end] == sequence[start:end]
            assert alignment[0] == "CTGGCGTCCGGAGCTGGGACGT"
            assert np.array_equal(alignment.coordinates, np.array([[19, 41]]))
            assert str(alignment) == """\
mm9.fa           19 CTGGCGTCCGGAGCTGGGACGT 41
"""
            assert alignment.format("mauve", metadata, identifiers) == """\
> 3:20-41 + mm9.fa
CTGGCGTCCGGAGCTGGGACGT
=
"""
            assert np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['C', 'T', 'G', 'G', 'C', 'G', 'T', 'C', 'C', 'G', 'G',
           'A', 'G', 'C', 'T', 'G', 'G', 'G', 'A', 'C', 'G', 'T']], dtype='U')
                    # fmt: on
                )
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
            with pytest.raises(StopIteration):
                next(alignments)
        # As each nucleotide in each sequence is stored exactly once in an XMFA
        # file, we can reconstitute the full sequences:
        assert len(saved_alignments) == 2
        filenames = []
        for alignment in saved_alignments:
            for record in alignment.sequences:
                filename = record.id
                filenames.append(filename)
        filenames = set(filenames)
        n = len(filenames)
        assert n == 3
        lengths = dict.fromkeys(filenames, 0)
        for alignment in saved_alignments:
            for record in alignment.sequences:
                filename = record.id
                length = len(record.seq)
                if length > lengths[filename]:
                    lengths[filename] = length
        assert lengths["equCab1.fa"] == 50
        assert lengths["canFam2.fa"] == 49
        assert lengths["mm9.fa"] == 41
        sequences = {}
        for filename, length in lengths.items():
            sequences[filename] = MutableSeq("N" * length)
        # Now fill up the sequences:
        for alignment in saved_alignments:
            for row, record in zip(alignment.coordinates, alignment.sequences):
                filename = record.id
                start = row[0]
                end = row[-1]
                if start > end:
                    start, end = end, start
                sequences[filename][start:end] = record.seq[start:end]
        # Confirm that the fully defined sequences agree with the Fasta file:
        for filename, sequence in sequences.items():
            sequences[filename] = Seq(sequence)
            assert sequences[filename] == self.sequences[filename]
        # Make sure we can replace the partially defined sequences by these
        # fully defined sequences, and get the same alignment:
        alignment = saved_alignments[0]
        for record in alignment.sequences:
            filename = record.id
            record.seq = sequences[filename]
            assert alignment[0] == "TAAGCCCTCCTAGCACACACCCGGAGTGGCC-GGGCCGTAC-TTTCCTTTTC"
            assert alignment[1] == "CAAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGC-TTTCCTTTTC"
            assert alignment[2] == "---------------------------------GGATCTACTTTTCCTCTTC"
        alignment = saved_alignments[1]
        for record in alignment.sequences:
            filename = record.id
            record.seq = sequences[filename]
            assert alignment[0] == "CTGGCGTCCGGAGCTGGGACGT"

    def test_write_read(self):
        path = os.path.join("Mauve", "separate.xmfa")
        with open(path) as stream:
            data = stream.read()

        stream = StringIO()
        stream.write(data)
        stream.seek(0)
        alignments = Align.parse(stream, "mauve")
        output = StringIO()
        n = Align.write(alignments, output, "mauve")
        assert n == 2
        output.seek(0)
        assert output.read() == data


class TestMauveBasic(unittest.TestCase):
    def test_empty(self):
        stream = StringIO()
        with pytest.raises(ValueError, match="Empty file."):
            Align.parse(stream, "mauve")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
