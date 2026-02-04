"""Tests for Bio.SeqUtils.KmerCounter module."""

import unittest
from Bio import MissingExternalDependencyError
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.KmerCounter import KmerCounter, count_kmers


class KmerCounterBasicTests(unittest.TestCase):
    """Basic tests (work with both Rust and Python)."""

    def test_simple_count(self):
        """Count k-mers in simple sequence."""
        counter = KmerCounter(k=3)
        counts = counter.count("ACGTACGT")
        self.assertEqual(counts["ACG"], 2)
        self.assertEqual(counts["CGT"], 2)
        self.assertEqual(counts["GTA"], 1)
        self.assertEqual(counts["TAC"], 1)

    def test_canonical_kmers(self):
        """Test canonical k-mer counting."""
        counter = KmerCounter(k=3, canonical=True)
        counts = counter.count("ACGTACGT")
        # ACG and CGT are reverse complements
        # With canonical mode, they should be counted together
        # The canonical form will be ACG (lexicographically smaller)
        self.assertEqual(counts["ACG"], 4)

    def test_input_types(self):
        """Test str, Seq, SeqRecord inputs."""
        counter = KmerCounter(k=3)
        counts1 = counter.count("ACGTACGT")
        counts2 = counter.count(Seq("ACGTACGT"))
        rec = SeqRecord(Seq("ACGTACGT"), id="test")
        counts3 = counter.count(rec)
        self.assertEqual(counts1, counts2)
        self.assertEqual(counts2, counts3)

    def test_iterable_input(self):
        """Test iterable of sequences."""
        counter = KmerCounter(k=3)
        seqs = ["ACGT", "ACGT"]
        counts = counter.count(seqs)
        self.assertEqual(counts["ACG"], 2)

    def test_ambiguous_bases(self):
        """K-mers with N should be skipped."""
        counter = KmerCounter(k=3)
        counts = counter.count("ACNGT")
        self.assertEqual(len(counts), 0)

    def test_empty_sequence(self):
        """Empty sequence returns empty dict."""
        counter = KmerCounter(k=3)
        counts = counter.count("")
        self.assertEqual(counts, {})

    def test_sequence_shorter_than_k(self):
        """Sequence shorter than k returns empty dict."""
        counter = KmerCounter(k=10)
        counts = counter.count("ACGT")
        self.assertEqual(counts, {})

    def test_case_insensitive(self):
        """Test case insensitivity."""
        counter = KmerCounter(k=3)
        counts1 = counter.count("acgt")
        counts2 = counter.count("ACGT")
        self.assertEqual(counts1, counts2)

    def test_convenience_function(self):
        """Test count_kmers convenience function."""
        counts = count_kmers("ACGTACGT", k=3)
        self.assertEqual(counts["ACG"], 2)

    def test_invalid_k(self):
        """Test invalid k values."""
        with self.assertRaises(ValueError):
            KmerCounter(k=0)
        with self.assertRaises(ValueError):
            KmerCounter(k=-1)

    def test_multiple_sequences(self):
        """Test counting across multiple sequences."""
        counter = KmerCounter(k=3)
        seqs = ["ACGTACGT", "TGCATGCA"]
        counts = counter.count(seqs)
        # ACGTACGT: ACG=2, CGT=2, GTA=2, TAC=1
        # TGCATGCA: TGC=2, GCA=2, CAT=1, ATG=1
        self.assertEqual(counts["ACG"], 2)
        self.assertEqual(counts["TGC"], 2)

    def test_only_valid_bases(self):
        """Test that only ACGT bases are counted."""
        counter = KmerCounter(k=3)
        counts = counter.count("ACGTRYSWKM")
        # Should only count ACG and CGT (before first ambiguous base)
        self.assertEqual(counts.get("ACG", 0), 1)
        self.assertEqual(counts.get("CGT", 0), 1)
        # Ambiguous bases should stop k-mer extraction
        self.assertNotIn("GTR", counts)


class KmerCounterRustTests(unittest.TestCase):
    """Tests specific to Rust implementation."""

    def setUp(self):
        """Check if Rust extension available."""
        from Bio.SeqUtils.KmerCounter import RUST_AVAILABLE

        if not RUST_AVAILABLE:
            raise MissingExternalDependencyError(
                "Rust k-mer counter extension not available"
            )

    def test_rust_is_used(self):
        """Verify Rust implementation is being used."""
        counter = KmerCounter(k=3, use_rust=True)
        self.assertTrue(counter.using_rust)

    def test_rust_vs_python_consistency(self):
        """Rust and Python should give same results."""
        test_seq = "ACGTACGTNNACGTACGT" * 10

        counter_rust = KmerCounter(k=5, use_rust=True)
        counter_python = KmerCounter(k=5, use_rust=False)

        counts_rust = counter_rust.count(test_seq)
        counts_python = counter_python.count(test_seq)

        self.assertEqual(counts_rust, counts_python)

    def test_large_k_fallback(self):
        """k >= 32 should fall back to Python."""
        counter = KmerCounter(k=32, use_rust=True)
        self.assertFalse(counter.using_rust)

    def test_rust_canonical_vs_python_canonical(self):
        """Test canonical mode consistency between Rust and Python."""
        test_seq = "ACGTACGTTGCATGCA" * 5

        counter_rust = KmerCounter(k=7, canonical=True, use_rust=True)
        counter_python = KmerCounter(k=7, canonical=True, use_rust=False)

        counts_rust = counter_rust.count(test_seq)
        counts_python = counter_python.count(test_seq)

        self.assertEqual(counts_rust, counts_python)

    def test_rust_performance_smoke_test(self):
        """Smoke test that Rust implementation completes without error on large input."""
        # Generate a reasonably large sequence
        import random

        random.seed(42)
        seq = "".join(random.choice("ACGT") for _ in range(100000))

        counter = KmerCounter(k=21, use_rust=True)
        counts = counter.count(seq)

        # Just verify it completes and returns reasonable results
        self.assertGreater(len(counts), 0)
        self.assertLess(len(counts), 4**21)  # Can't have more than 4^k unique k-mers


class KmerCounterStreamingTests(unittest.TestCase):
    """Tests for streaming_count_kmers function."""

    def setUp(self):
        """Create temporary test FASTA file."""
        import tempfile
        from Bio import SeqIO

        self.temp_dir = tempfile.mkdtemp()
        self.fasta_path = f"{self.temp_dir}/test.fasta"

        # Create test FASTA with multiple records
        records = [
            SeqRecord(Seq("ACGTACGT"), id="seq1", description=""),
            SeqRecord(Seq("TGCATGCA"), id="seq2", description=""),
            SeqRecord(Seq("NNNNNNNN"), id="seq3", description=""),
        ]
        SeqIO.write(records, self.fasta_path, "fasta")

    def tearDown(self):
        """Clean up temporary files."""
        import shutil

        shutil.rmtree(self.temp_dir)

    def test_streaming_fasta(self):
        """Test streaming k-mer counting from FASTA file."""
        from Bio.SeqUtils.KmerCounter import streaming_count_kmers

        counts = streaming_count_kmers(self.fasta_path, k=3)

        # Verify expected k-mers from seq1 and seq2
        self.assertEqual(counts["ACG"], 2)
        self.assertEqual(counts["TGC"], 2)
        # seq3 should contribute nothing (all N's)

    def test_streaming_canonical(self):
        """Test streaming with canonical k-mers."""
        from Bio.SeqUtils.KmerCounter import streaming_count_kmers

        counts = streaming_count_kmers(self.fasta_path, k=3, canonical=True)

        # ACG and CGT should be counted together
        self.assertEqual(counts["ACG"], 4)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
