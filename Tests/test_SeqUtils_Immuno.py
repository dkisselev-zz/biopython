# Copyright 2024 by the Biopython contributors. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Tests for Bio.SeqUtils.Immuno module."""

import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.Immuno import cdr3_distance
from Bio.SeqUtils.Immuno import gene_usage
from Bio.SeqUtils.Immuno import group_clonotypes
from Bio.SeqUtils.Immuno import translate_junction


def _make_record(seq_id, v_call, j_call=None, junction_aa=None):
    """Helper to build a minimal SeqRecord with AIRR-style annotations."""
    annotations = {"v_call": v_call}
    if j_call is not None:
        annotations["j_call"] = j_call
    if junction_aa is not None:
        annotations["junction_aa"] = junction_aa
    return SeqRecord(
        Seq("ACGT"), id=seq_id, name=seq_id, description="",
        annotations=annotations,
    )


class TestTranslateJunction(unittest.TestCase):
    """Tests for translate_junction."""

    def test_valid_productive(self):
        self.assertEqual(translate_junction("TGTGCAAGAGATCGACTG"), "CARDRL")

    def test_valid_productive_carldy(self):
        self.assertEqual(translate_junction("TGTGCAAGACTGGACTAC"), "CARLDY")

    def test_stop_codon_productive(self):
        # TGT=C, TAA=stop -> None in productive mode
        self.assertIsNone(translate_junction("TGTTAA"))

    def test_non_mod3_productive(self):
        # 4 nt is not divisible by 3 -> None
        self.assertIsNone(translate_junction("TGTG"))

    def test_non_mod3_non_productive(self):
        # Truncates TGTG to TGT -> 'C'
        self.assertEqual(translate_junction("TGTG", productive=False), "C")

    def test_empty_string(self):
        self.assertIsNone(translate_junction(""))

    def test_empty_non_productive(self):
        self.assertIsNone(translate_junction("", productive=False))

    def test_valid_short(self):
        # TGTGCA = 6 nt = 2 codons, valid
        self.assertEqual(translate_junction("TGTGCA"), "CA")

    def test_stop_codon_non_productive(self):
        # Stop codon allowed when productive=False
        result = translate_junction("TGTTAA", productive=False)
        self.assertIn("*", result)


class TestCdr3Distance(unittest.TestCase):
    """Tests for cdr3_distance."""

    def test_hamming_identical(self):
        self.assertEqual(cdr3_distance("CARDRL", "CARDRL"), 0)

    def test_hamming_mismatch(self):
        self.assertEqual(cdr3_distance("CARDRL", "CARLDY"), 3)

    def test_hamming_length_mismatch_raises(self):
        self.assertRaises(ValueError, cdr3_distance, "CARDRL", "CARDRLL")

    def test_levenshtein_insertion(self):
        self.assertEqual(
            cdr3_distance("CARDRL", "CARDRLL", metric="levenshtein"), 1
        )

    def test_levenshtein_identical(self):
        self.assertEqual(
            cdr3_distance("CARDRL", "CARDRL", metric="levenshtein"), 0
        )

    def test_levenshtein_deletion(self):
        self.assertEqual(
            cdr3_distance("CARDRLL", "CARDRL", metric="levenshtein"), 1
        )

    def test_unknown_metric_raises(self):
        self.assertRaises(ValueError, cdr3_distance, "A", "B", metric="bogus")


class TestGroupClonotypes(unittest.TestCase):
    """Tests for group_clonotypes."""

    def test_basic_grouping(self):
        records = [
            _make_record("r1", "IGHV1-2*01", "IGHJ4*02", "CARDRL"),
            _make_record("r2", "IGHV3-30*01", "IGHJ6*01", "CARLDY"),
            _make_record("r3", "IGHV1-2*01", "IGHJ4*02", "CARDRL"),
        ]
        clones = group_clonotypes(records)
        self.assertEqual(len(clones), 2)
        key = ("IGHV1-2", "IGHJ4", "CARDRL")
        self.assertIn(key, clones)
        self.assertEqual(len(clones[key]), 2)

    def test_multi_allele_normalization(self):
        records = [
            _make_record("r1", "IGHV1-2*01,IGHV1-2*02", "IGHJ4*02", "CARDRL"),
        ]
        clones = group_clonotypes(records)
        self.assertIn(("IGHV1-2", "IGHJ4", "CARDRL"), clones)

    def test_allele_stripping(self):
        records = [
            _make_record("r1", "IGHV1-2*01", "IGHJ4*01", "CARDRL"),
            _make_record("r2", "IGHV1-2*02", "IGHJ4*02", "CARDRL"),
        ]
        clones = group_clonotypes(records)
        # Both should map to the same clonotype after allele stripping
        self.assertEqual(len(clones), 1)
        self.assertIn(("IGHV1-2", "IGHJ4", "CARDRL"), clones)


class TestGeneUsage(unittest.TestCase):
    """Tests for gene_usage."""

    def test_basic_usage(self):
        records = [
            _make_record("r1", "IGHV1-2*01"),
            _make_record("r2", "IGHV1-2*01"),
            _make_record("r3", "IGHV3-30*01"),
        ]
        counter = gene_usage(records)
        self.assertEqual(counter["IGHV1-2"], 2)
        self.assertEqual(counter["IGHV3-30"], 1)

    def test_no_strip_allele(self):
        records = [_make_record("r1", "IGHV1-2*01")]
        counter = gene_usage(records, strip_allele=False)
        self.assertEqual(counter["IGHV1-2*01"], 1)
        self.assertNotIn("IGHV1-2", counter)

    def test_multi_allele_splitting(self):
        records = [_make_record("r1", "IGHV1-2*01,IGHV1-2*02")]
        counter = gene_usage(records)
        # Each allele counted separately, both strip to IGHV1-2
        self.assertEqual(counter["IGHV1-2"], 2)

    def test_empty_call_skipped(self):
        records = [_make_record("r1", "")]
        counter = gene_usage(records)
        self.assertEqual(len(counter), 0)

    def test_none_call_skipped(self):
        r = SeqRecord(
            Seq("ACGT"), id="r1", name="r1", description="",
            annotations={"v_call": None},
        )
        counter = gene_usage([r])
        self.assertEqual(len(counter), 0)

    def test_j_call_field(self):
        records = [
            SeqRecord(
                Seq("ACGT"), id="r1", name="r1", description="",
                annotations={"v_call": "IGHV1-2*01", "j_call": "IGHJ4*02"},
            )
        ]
        counter = gene_usage(records, gene="j_call")
        self.assertEqual(counter["IGHJ4"], 1)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
