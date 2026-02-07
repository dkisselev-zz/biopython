# Copyright 2024 by Biopython Contributors. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Tests for Bio.SeqUtils.Immuno (immune repertoire utilities)."""

import unittest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.Immuno import (
    cdr3_distance,
    count_gene_usage,
    group_clonotypes,
    hamming_distance,
    levenshtein_distance,
    translate_cdr3,
)


class TestTranslateCDR3(unittest.TestCase):
    """Tests for translate_cdr3 function."""

    def test_translate_valid(self):
        """Test successful translation of valid CDR3."""
        # TGTGCGAGA = CAR (Cys-Ala-Arg)
        result = translate_cdr3("TGTGCGAGA")
        self.assertEqual(result, "CAR")

        # Longer sequence
        result = translate_cdr3("TGTGCGAGAGATAGGCTC")
        self.assertEqual(result, "CARDRL")

        # Full junction (must be divisible by 3)
        # 63 bases = 21 amino acids
        result = translate_cdr3(
            "TGTGCGAGAGATAGGCTCGGGACTAGCTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCA"
        )
        self.assertEqual(result, "CARDRLGTSYWGQGTLVTVSS")

    def test_translate_empty(self):
        """Test translation of empty sequence."""
        result = translate_cdr3("")
        self.assertEqual(result, "")

    def test_translate_invalid_length(self):
        """Test that frameshift sequences return None."""
        # Length not divisible by 3
        result = translate_cdr3("TGTGCGAGAT")  # 10 bases
        self.assertIsNone(result)

        result = translate_cdr3("TGTGCGAGATGG")  # 12 bases (should work, no stop)
        self.assertIsNotNone(result)

        result = translate_cdr3("TGTGCGAGATAGC")  # 13 bases
        self.assertIsNone(result)

    def test_translate_stop_codon(self):
        """Test that sequences with stop codons return None."""
        # TAA is stop codon
        result = translate_cdr3("TGTGCGAGATAATAA")  # 15 bases, contains TAA
        self.assertIsNone(result)

        # TAG is stop codon (in middle)
        result = translate_cdr3("TGTGCGAGATAG")  # 12 bases, TAG at end
        self.assertIsNone(result)

        # Sequence without stop codon should work
        result = translate_cdr3("TGTGCGAGATGG")  # 12 bases, no stop
        self.assertIsNotNone(result)

        # TGA is stop codon
        result = translate_cdr3("TGTGCGAGATGA")  # 12 bases, TGA at end
        self.assertIsNone(result)

    def test_translate_invalid_nucleotides(self):
        """Test that invalid nucleotides return None."""
        result = translate_cdr3("TGTGCGAGAXYZ")
        self.assertIsNone(result)

    def test_translate_lowercase(self):
        """Test that lowercase sequences work."""
        result = translate_cdr3("tgtgcgaga")
        self.assertEqual(result, "CAR")

    def test_translate_with_productive_false(self):
        """Test that productive=False returns None immediately."""
        # Even with valid sequence, should return None if productive=False
        result = translate_cdr3("TGTGCGAGAGATAGGCTC", productive=False)
        self.assertIsNone(result)

        # Empty sequence with productive=False
        result = translate_cdr3("", productive=False)
        self.assertIsNone(result)

    def test_translate_with_productive_true(self):
        """Test that productive=True proceeds with normal translation."""
        result = translate_cdr3("TGTGCGAGAGATAGGCTC", productive=True)
        self.assertEqual(result, "CARDRL")

        # Frameshifted sequence still returns None
        result = translate_cdr3("TGTGCGAGAT", productive=True)
        self.assertIsNone(result)

    def test_translate_with_seqrecord(self):
        """Test translation from SeqRecord with AIRR annotations."""
        # Create SeqRecord with AIRR annotations
        record = SeqRecord(Seq("ATGCGTACG"), id="SEQ001")
        record.annotations["airr"] = {
            "junction": "TGTGCGAGAGATAGGCTC",
            "productive": True,
        }

        result = translate_cdr3(record)
        self.assertEqual(result, "CARDRL")

    def test_translate_seqrecord_productive_false(self):
        """Test SeqRecord with productive=False in annotations."""
        record = SeqRecord(Seq("ATGCGTACG"), id="SEQ002")
        record.annotations["airr"] = {
            "junction": "TGTGCGAGAGATAGGCTC",
            "productive": False,  # Non-productive annotation
        }

        # Should return None immediately based on annotation
        result = translate_cdr3(record)
        self.assertIsNone(result)

    def test_translate_seqrecord_explicit_productive_override(self):
        """Test that explicit productive parameter overrides annotation."""
        record = SeqRecord(Seq("ATGCGTACG"), id="SEQ003")
        record.annotations["airr"] = {
            "junction": "TGTGCGAGAGATAGGCTC",
            "productive": True,  # Annotation says productive
        }

        # But explicit parameter says non-productive
        result = translate_cdr3(record, productive=False)
        self.assertIsNone(result)

    def test_translate_seqrecord_missing_junction(self):
        """Test SeqRecord without junction field."""
        record = SeqRecord(Seq("ATGCGTACG"), id="SEQ004")
        record.annotations["airr"] = {
            "v_call": "IGHV1-69*01",
            # No junction field
        }

        # Should warn and return None
        import warnings

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = translate_cdr3(record)

            self.assertIsNone(result)
            self.assertEqual(len(w), 1)
            self.assertIn("no 'junction' field", str(w[0].message))

    def test_translate_seqrecord_custom_junction_field(self):
        """Test SeqRecord with custom junction field name."""
        record = SeqRecord(Seq("ATGCGTACG"), id="SEQ005")
        record.annotations["airr"] = {
            "cdr3": "TGTGCGAGAGATAGGCTC",  # Using cdr3 instead of junction
            "productive": True,
        }

        # Should use cdr3 field when specified
        result = translate_cdr3(record, junction_field="cdr3")
        self.assertEqual(result, "CARDRL")

        # Should fail with default junction_field
        import warnings

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = translate_cdr3(record)  # Defaults to "junction"
            self.assertIsNone(result)

    def test_translate_seqrecord_from_airr_file(self):
        """Test translation using actual SeqRecord from AIRR file."""
        # Read first productive record from test file
        for record in SeqIO.parse("Airr/minimal.tsv", "airr"):
            airr = record.annotations["airr"]
            if airr.get("productive") is True:
                # Get junction from annotations
                junction = airr.get("junction")
                if junction:
                    # Translate using SeqRecord
                    result1 = translate_cdr3(record)

                    # Translate using string
                    result2 = translate_cdr3(junction)

                    # Should give same result
                    if result1 is not None and result2 is not None:
                        self.assertEqual(result1, result2)
                break


class TestHammingDistance(unittest.TestCase):
    """Tests for hamming_distance function."""

    def test_hamming_basic(self):
        """Test basic Hamming distance calculation."""
        # One mismatch
        self.assertEqual(hamming_distance("CARD", "CARE"), 1)

        # No mismatches
        self.assertEqual(hamming_distance("CARD", "CARD"), 0)

        # Two mismatches
        self.assertEqual(hamming_distance("CARD", "RARE"), 2)

        # All mismatches
        self.assertEqual(hamming_distance("AAAA", "TTTT"), 4)

    def test_hamming_case_insensitive(self):
        """Test that Hamming distance is case-insensitive."""
        self.assertEqual(hamming_distance("CARD", "card"), 0)
        self.assertEqual(hamming_distance("Card", "CaRe"), 1)

    def test_hamming_empty(self):
        """Test Hamming distance with empty sequences."""
        self.assertEqual(hamming_distance("", ""), 0)

    def test_hamming_different_lengths(self):
        """Test that different lengths raise ValueError."""
        with self.assertRaises(ValueError) as cm:
            hamming_distance("CARD", "CAR")
        self.assertIn("same length", str(cm.exception))

        with self.assertRaises(ValueError):
            hamming_distance("CARD", "CARDS")


class TestLevenshteinDistance(unittest.TestCase):
    """Tests for levenshtein_distance function."""

    def test_levenshtein_basic(self):
        """Test basic Levenshtein distance calculation."""
        # One deletion
        self.assertEqual(levenshtein_distance("CARD", "CAR"), 1)

        # One insertion
        self.assertEqual(levenshtein_distance("CAR", "CARD"), 1)

        # One substitution
        self.assertEqual(levenshtein_distance("CARD", "CARE"), 1)
        self.assertEqual(levenshtein_distance("CARD", "HARD"), 1)

        # No edits
        self.assertEqual(levenshtein_distance("CARD", "CARD"), 0)

    def test_levenshtein_multiple_edits(self):
        """Test multiple edits."""
        # Two substitutions
        self.assertEqual(levenshtein_distance("CARD", "HARE"), 2)

        # One deletion + one substitution
        self.assertEqual(levenshtein_distance("CARD", "HAR"), 2)

        # Complete replacement
        self.assertEqual(levenshtein_distance("AAAA", "TTTT"), 4)

    def test_levenshtein_empty(self):
        """Test Levenshtein distance with empty sequences."""
        self.assertEqual(levenshtein_distance("", ""), 0)
        self.assertEqual(levenshtein_distance("CARD", ""), 4)
        self.assertEqual(levenshtein_distance("", "CARD"), 4)

    def test_levenshtein_case_insensitive(self):
        """Test that Levenshtein distance is case-insensitive."""
        self.assertEqual(levenshtein_distance("CARD", "card"), 0)
        self.assertEqual(levenshtein_distance("Card", "CaRe"), 1)

    def test_levenshtein_longer_sequences(self):
        """Test with longer sequences."""
        seq1 = "CARDRLGTSWGQGTLVTVSS"
        seq2 = "CARDYYDSSGYYPRGATNYYGQGTLVTVSS"

        # Should be able to compute without error
        distance = levenshtein_distance(seq1, seq2)
        self.assertIsInstance(distance, int)
        self.assertGreater(distance, 0)


class TestCDR3Distance(unittest.TestCase):
    """Tests for unified cdr3_distance function."""

    def test_default_method(self):
        """Test default method (Levenshtein)."""
        # Default should be Levenshtein
        self.assertEqual(cdr3_distance("CARD", "CAR"), 1)
        self.assertEqual(cdr3_distance("CARD", "CARE"), 1)
        self.assertEqual(cdr3_distance("CARD", "CARD"), 0)

    def test_levenshtein_method(self):
        """Test explicit Levenshtein method."""
        self.assertEqual(cdr3_distance("CARD", "CAR", method="levenshtein"), 1)
        self.assertEqual(cdr3_distance("CARD", "CARE", method="levenshtein"), 1)
        self.assertEqual(
            cdr3_distance("CARDRLGTS", "CARDYYGTS", method="levenshtein"), 2
        )

    def test_hamming_method(self):
        """Test Hamming distance method."""
        self.assertEqual(cdr3_distance("CARD", "CARE", method="hamming"), 1)
        self.assertEqual(cdr3_distance("CARD", "CARD", method="hamming"), 0)
        self.assertEqual(cdr3_distance("CARD", "RARE", method="hamming"), 2)

    def test_case_insensitive_method_name(self):
        """Test that method names are case-insensitive."""
        self.assertEqual(cdr3_distance("CARD", "CARE", method="HAMMING"), 1)
        self.assertEqual(cdr3_distance("CARD", "CAR", method="Levenshtein"), 1)
        self.assertEqual(cdr3_distance("CARD", "CAR", method="LEVENSHTEIN"), 1)

    def test_case_insensitive_sequences(self):
        """Test that sequences are case-insensitive."""
        self.assertEqual(cdr3_distance("CARD", "card"), 0)
        self.assertEqual(cdr3_distance("Card", "CaRe", method="hamming"), 1)

    def test_invalid_method(self):
        """Test that invalid method raises ValueError."""
        with self.assertRaises(ValueError) as cm:
            cdr3_distance("CARD", "CARE", method="invalid")
        self.assertIn("Unknown distance method", str(cm.exception))
        self.assertIn("invalid", str(cm.exception))

    def test_hamming_different_lengths_raises(self):
        """Test that Hamming with different lengths raises ValueError."""
        with self.assertRaises(ValueError) as cm:
            cdr3_distance("CARD", "CAR", method="hamming")
        self.assertIn("same length", str(cm.exception))

    def test_cdr3_examples(self):
        """Test with realistic CDR3 sequences."""
        # Same clonotype (identical)
        cdr3_1 = "CARDRLGTS"
        cdr3_2 = "CARDRLGTS"
        self.assertEqual(cdr3_distance(cdr3_1, cdr3_2), 0)

        # Similar clonotypes (few differences)
        cdr3_3 = "CARDRLGTS"
        cdr3_4 = "CARDYYGTS"  # Different middle residues
        dist = cdr3_distance(cdr3_3, cdr3_4)
        self.assertGreater(dist, 0)
        self.assertLess(dist, 5)

        # Different length CDR3s (Levenshtein handles this)
        cdr3_5 = "CARDRLGTS"
        cdr3_6 = "CARDRLGTSWGQ"  # Longer
        dist = cdr3_distance(cdr3_5, cdr3_6)
        self.assertEqual(dist, 3)  # 3 extra characters (WGQ)


class TestGroupClonotypes(unittest.TestCase):
    """Tests for group_clonotypes function."""

    def test_group_basic(self):
        """Test basic clonotype grouping with nucleotide sequences (default)."""
        # Create test records with nucleotide junctions
        rec1 = SeqRecord(Seq("ATG"), id="SEQ001")
        rec1.annotations["airr"] = {
            "v_call": "IGHV1-69*01",
            "j_call": "IGHJ4*02",
            "junction": "TGTGCGAGAGATAGGCTC",  # CARDRL
            "junction_aa": "CARDRL",
        }

        rec2 = SeqRecord(Seq("ATG"), id="SEQ002")
        rec2.annotations["airr"] = {
            "v_call": "IGHV1-69*01",
            "j_call": "IGHJ4*02",
            "junction": "TGTGCGAGAGATAGGCTC",  # Same nucleotide sequence
            "junction_aa": "CARDRL",
        }

        rec3 = SeqRecord(Seq("ATG"), id="SEQ003")
        rec3.annotations["airr"] = {
            "v_call": "IGHV1-2*02",
            "j_call": "IGHJ6*02",
            "junction": "TGTGCGAGAGATTATTAT",  # CARDYY
            "junction_aa": "CARDYY",
        }

        records = [rec1, rec2, rec3]

        # Default: nucleotide-based clonotyping
        clonotypes = group_clonotypes(records)

        # Should have 2 clonotypes
        self.assertEqual(len(clonotypes), 2)

        # First clonotype should have 2 sequences (same nucleotide junction)
        key1 = ("IGHV1-69", "IGHJ4", "TGTGCGAGAGATAGGCTC")
        self.assertIn(key1, clonotypes)
        self.assertEqual(len(clonotypes[key1]), 2)

        # Second clonotype should have 1 sequence
        key2 = ("IGHV1-2", "IGHJ6", "TGTGCGAGAGATTATTAT")
        self.assertIn(key2, clonotypes)
        self.assertEqual(len(clonotypes[key2]), 1)

    def test_group_amino_acid_mode(self):
        """Test functional grouping with amino acid sequences."""
        # Create test records
        rec1 = SeqRecord(Seq("ATG"), id="SEQ001")
        rec1.annotations["airr"] = {
            "v_call": "IGHV1-69*01",
            "j_call": "IGHJ4*02",
            "junction": "TGTGCGAGAGATAGGCTC",  # CARDRL
            "junction_aa": "CARDRL",
        }

        rec2 = SeqRecord(Seq("ATG"), id="SEQ002")
        rec2.annotations["airr"] = {
            "v_call": "IGHV1-69*01",
            "j_call": "IGHJ4*02",
            "junction": "TGTGCGAGAGATAGGCTT",  # Different nucleotides, same AA
            "junction_aa": "CARDRL",
        }

        records = [rec1, rec2]

        # Nucleotide-based: should be 2 different clonotypes
        clonotypes_nt = group_clonotypes(records)
        self.assertEqual(len(clonotypes_nt), 2)

        # Amino acid-based: should be 1 functional group
        clonotypes_aa = group_clonotypes(records, cdr3_field="junction_aa")
        self.assertEqual(len(clonotypes_aa), 1)
        key_aa = ("IGHV1-69", "IGHJ4", "CARDRL")
        self.assertIn(key_aa, clonotypes_aa)
        self.assertEqual(len(clonotypes_aa[key_aa]), 2)

    def test_group_strip_alleles(self):
        """Test allele stripping in clonotype grouping."""
        rec1 = SeqRecord(Seq("ATG"), id="SEQ001")
        rec1.annotations["airr"] = {
            "v_call": "IGHV1-69*01",
            "j_call": "IGHJ4*02",
            "junction": "TGTGCGAGAGATAGGCTC",
            "junction_aa": "CARDRL",
        }

        rec2 = SeqRecord(Seq("ATG"), id="SEQ002")
        rec2.annotations["airr"] = {
            "v_call": "IGHV1-69*06",  # Different allele
            "j_call": "IGHJ4*01",  # Different allele
            "junction": "TGTGCGAGAGATAGGCTC",  # Same junction
            "junction_aa": "CARDRL",
        }

        records = [rec1, rec2]

        # With allele stripping (default)
        clonotypes = group_clonotypes(records, strip_alleles=True)
        self.assertEqual(len(clonotypes), 1)  # Same clonotype

        # Without allele stripping
        clonotypes = group_clonotypes(records, strip_alleles=False)
        self.assertEqual(len(clonotypes), 2)  # Different clonotypes

    def test_group_missing_data(self):
        """Test grouping with missing gene assignments."""
        rec1 = SeqRecord(Seq("ATG"), id="SEQ001")
        rec1.annotations["airr"] = {
            "v_call": "IGHV1-69*01",
            "j_call": None,  # Missing J gene
            "junction": "TGTGCGAGAGATAGGCTC",
            "junction_aa": "CARDRL",
        }

        rec2 = SeqRecord(Seq("ATG"), id="SEQ002")
        rec2.annotations["airr"] = {
            "v_call": "IGHV1-69*01",
            "j_call": None,  # Missing J gene
            "junction": "TGTGCGAGAGATAGGCTC",
            "junction_aa": "CARDRL",
        }

        records = [rec1, rec2]
        clonotypes = group_clonotypes(records)

        # Should group together with None J gene
        key = ("IGHV1-69", None, "TGTGCGAGAGATAGGCTC")
        self.assertIn(key, clonotypes)
        self.assertEqual(len(clonotypes[key]), 2)

    def test_group_codon_degeneracy(self):
        """Test that nucleotide mode distinguishes synonymous codons.

        This test demonstrates why nucleotide-based clonotyping is more
        precise than amino acid-based grouping.
        """
        # Create sequences with same amino acid but different nucleotides
        # ARG can be encoded as: CGT, CGC, CGA, CGG, AGA, AGG
        rec1 = SeqRecord(Seq("ATG"), id="SEQ001")
        rec1.annotations["airr"] = {
            "v_call": "IGHV1-69*01",
            "j_call": "IGHJ4*02",
            "junction": "TGTGCGCGT",  # TGT GCG CGT -> CAR (Cys-Ala-Arg)
            "junction_aa": "CAR",
        }

        rec2 = SeqRecord(Seq("ATG"), id="SEQ002")
        rec2.annotations["airr"] = {
            "v_call": "IGHV1-69*01",
            "j_call": "IGHJ4*02",
            "junction": "TGTGCGCGC",  # TGT GCG CGC -> CAR (same AA, different nt)
            "junction_aa": "CAR",
        }

        rec3 = SeqRecord(Seq("ATG"), id="SEQ003")
        rec3.annotations["airr"] = {
            "v_call": "IGHV1-69*01",
            "j_call": "IGHJ4*02",
            "junction": "TGTGCGAGA",  # TGT GCG AGA -> CAR (same AA, different nt)
            "junction_aa": "CAR",
        }

        records = [rec1, rec2, rec3]

        # Nucleotide-based (default): should be 3 separate clonotypes
        # This is more accurate for clonal lineage analysis
        clonotypes_nt = group_clonotypes(records)
        self.assertEqual(
            len(clonotypes_nt),
            3,
            "Nucleotide mode should distinguish synonymous codons",
        )

        # Amino acid-based: should be 1 functional group
        # Useful for identifying convergent recombination
        clonotypes_aa = group_clonotypes(records, cdr3_field="junction_aa")
        self.assertEqual(
            len(clonotypes_aa),
            1,
            "Amino acid mode should group synonymous variants",
        )

        # Verify the amino acid group contains all 3 sequences
        key_aa = ("IGHV1-69", "IGHJ4", "CAR")
        self.assertIn(key_aa, clonotypes_aa)
        self.assertEqual(len(clonotypes_aa[key_aa]), 3)

    def test_group_with_file(self):
        """Test grouping with actual AIRR file."""
        records = list(SeqIO.parse("Airr/standard.tsv", "airr"))
        clonotypes = group_clonotypes(records)

        # Should have multiple clonotypes
        self.assertGreater(len(clonotypes), 0)

        # Each clonotype should have at least one sequence
        for key, seqs in clonotypes.items():
            self.assertGreater(len(seqs), 0)
            # Check key structure
            self.assertEqual(len(key), 3)  # (V, J, CDR3)


class TestCountGeneUsage(unittest.TestCase):
    """Tests for count_gene_usage function."""

    def test_count_basic(self):
        """Test basic gene usage counting."""
        rec1 = SeqRecord(Seq("ATG"), id="SEQ001")
        rec1.annotations["airr"] = {"v_call": "IGHV1-69*01"}

        rec2 = SeqRecord(Seq("ATG"), id="SEQ002")
        rec2.annotations["airr"] = {"v_call": "IGHV1-69*01"}

        rec3 = SeqRecord(Seq("ATG"), id="SEQ003")
        rec3.annotations["airr"] = {"v_call": "IGHV1-2*02"}

        records = [rec1, rec2, rec3]
        v_usage = count_gene_usage(records, gene_field="v_call")

        # Should have 2 genes
        self.assertEqual(len(v_usage), 2)

        # Check counts
        self.assertEqual(v_usage["IGHV1-69"], 2)
        self.assertEqual(v_usage["IGHV1-2"], 1)

        # Should be sorted by count (descending)
        counts = list(v_usage.values())
        self.assertEqual(counts, sorted(counts, reverse=True))

    def test_count_strip_alleles(self):
        """Test allele stripping in gene counting."""
        rec1 = SeqRecord(Seq("ATG"), id="SEQ001")
        rec1.annotations["airr"] = {"v_call": "IGHV1-69*01"}

        rec2 = SeqRecord(Seq("ATG"), id="SEQ002")
        rec2.annotations["airr"] = {"v_call": "IGHV1-69*06"}

        records = [rec1, rec2]

        # With allele stripping (default)
        v_usage = count_gene_usage(records, strip_alleles=True)
        self.assertEqual(len(v_usage), 1)
        self.assertEqual(v_usage["IGHV1-69"], 2)

        # Without allele stripping
        v_usage = count_gene_usage(records, strip_alleles=False)
        self.assertEqual(len(v_usage), 2)
        self.assertEqual(v_usage["IGHV1-69*01"], 1)
        self.assertEqual(v_usage["IGHV1-69*06"], 1)

    def test_count_j_genes(self):
        """Test counting J genes."""
        rec1 = SeqRecord(Seq("ATG"), id="SEQ001")
        rec1.annotations["airr"] = {"j_call": "IGHJ4*02"}

        rec2 = SeqRecord(Seq("ATG"), id="SEQ002")
        rec2.annotations["airr"] = {"j_call": "IGHJ4*02"}

        rec3 = SeqRecord(Seq("ATG"), id="SEQ003")
        rec3.annotations["airr"] = {"j_call": "IGHJ6*02"}

        records = [rec1, rec2, rec3]
        j_usage = count_gene_usage(records, gene_field="j_call")

        self.assertEqual(len(j_usage), 2)
        self.assertEqual(j_usage["IGHJ4"], 2)
        self.assertEqual(j_usage["IGHJ6"], 1)

    def test_count_missing_genes(self):
        """Test counting with missing gene assignments."""
        rec1 = SeqRecord(Seq("ATG"), id="SEQ001")
        rec1.annotations["airr"] = {"v_call": "IGHV1-69*01"}

        rec2 = SeqRecord(Seq("ATG"), id="SEQ002")
        rec2.annotations["airr"] = {"v_call": None}

        rec3 = SeqRecord(Seq("ATG"), id="SEQ003")
        rec3.annotations["airr"] = {}  # Missing v_call

        records = [rec1, rec2, rec3]
        v_usage = count_gene_usage(records, gene_field="v_call")

        # Should count None separately
        self.assertIn("IGHV1-69", v_usage)
        self.assertIn(None, v_usage)
        self.assertEqual(v_usage[None], 2)

    def test_count_with_file(self):
        """Test counting with actual AIRR file."""
        records = list(SeqIO.parse("Airr/standard.tsv", "airr"))

        # Count V genes
        v_usage = count_gene_usage(records, gene_field="v_call")
        self.assertGreater(len(v_usage), 0)

        # Check sorting (descending by count)
        counts = list(v_usage.values())
        self.assertEqual(counts, sorted(counts, reverse=True))

        # Count J genes
        j_usage = count_gene_usage(records, gene_field="j_call")
        self.assertGreater(len(j_usage), 0)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
