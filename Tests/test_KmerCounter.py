# Copyright 2025 by the Biopython Contributors.
# All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Tests for Bio.SeqUtils.KmerCounter."""

import os
import random
import shutil
import subprocess
import tempfile
import time
import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import KmerCounter
from Bio.SeqUtils._kmer_counter import _count_kmers_python

# ---------------------------------------------------------------------------
# Module-level feature flags — do NOT raise MissingPythonDependencyError;
# the core class always works via the Python fallback so we must never skip
# the entire module.
# ---------------------------------------------------------------------------
try:
    from Bio.SeqUtils._kmer_rust import count_kmers as _rust_count_kmers

    HAS_RUST = True
except ImportError:
    HAS_RUST = False

HAS_JELLYFISH = shutil.which("jellyfish") is not None


# ===========================================================================
# 1.  Input-type normalisation
# ===========================================================================
class TestKmerCounterInputTypes(unittest.TestCase):
    """All supported input types produce the same counts."""

    def setUp(self):
        self.seq_str = "ATCGATCG"
        self.expected = {"ATC": 2, "TCG": 2, "CGA": 1, "GAT": 1}
        self.counter = KmerCounter(k=3)

    def test_str(self):
        self.assertEqual(self.counter.count(self.seq_str), self.expected)

    def test_seq(self):
        self.assertEqual(self.counter.count(Seq(self.seq_str)), self.expected)

    def test_seqrecord(self):
        rec = SeqRecord(Seq(self.seq_str), id="x")
        self.assertEqual(self.counter.count(rec), self.expected)

    def test_list_of_str(self):
        # Two halves whose combined k-mers equal the single-string result.
        # "ATCG" yields ATC, TCG  and  "GATCG" yields GAT, ATC, TCG
        # Together: ATC:2, TCG:2, GAT:1 — but CGA only comes from the full
        # string.  So use the full string as a single-element list instead.
        self.assertEqual(self.counter.count([self.seq_str]), self.expected)

    def test_list_of_seqrecord(self):
        recs = [SeqRecord(Seq(self.seq_str), id="x")]
        self.assertEqual(self.counter.count(recs), self.expected)

    def test_mixed_iterable(self):
        items = [self.seq_str, Seq(self.seq_str)]
        # Each element contributes the full set; counts should double.
        doubled = {k: v * 2 for k, v in self.expected.items()}
        self.assertEqual(self.counter.count(items), doubled)

    def test_lowercase_input(self):
        self.assertEqual(self.counter.count(self.seq_str.lower()), self.expected)

    def test_mixed_case_input(self):
        self.assertEqual(self.counter.count("AtCgAtCg"), self.expected)


# ===========================================================================
# 2.  Hand-verified correctness
# ===========================================================================
class TestKmerCounterCorrectness(unittest.TestCase):
    """Counts match hand-computed values."""

    def test_k2(self):
        # ACGT -> AC CG GT
        self.assertEqual(KmerCounter(k=2).count("ACGT"), {"AC": 1, "CG": 1, "GT": 1})

    def test_k3(self):
        self.assertEqual(
            KmerCounter(k=3).count("ATCGATCG"),
            {"ATC": 2, "TCG": 2, "CGA": 1, "GAT": 1},
        )

    def test_k4(self):
        # ATCGA -> ATCG TCGA
        self.assertEqual(KmerCounter(k=4).count("ATCGA"), {"ATCG": 1, "TCGA": 1})

    def test_multi_sequence(self):
        # Two sequences; counts are summed.
        counts = KmerCounter(k=2).count(["AC", "AC"])
        self.assertEqual(counts, {"AC": 2})

    def test_repeated_kmer(self):
        # AAAA -> AAA AAA  (two overlapping 3-mers)
        self.assertEqual(KmerCounter(k=3).count("AAAA"), {"AAA": 2})


# ===========================================================================
# 3.  Edge cases
# ===========================================================================
class TestKmerCounterEdgeCases(unittest.TestCase):
    """Boundary and degenerate inputs."""

    def test_empty_string(self):
        self.assertEqual(KmerCounter(k=3).count(""), {})

    def test_empty_list(self):
        self.assertEqual(KmerCounter(k=3).count([]), {})

    def test_k_greater_than_length(self):
        self.assertEqual(KmerCounter(k=5).count("ATCG"), {})

    def test_k_equals_length(self):
        self.assertEqual(KmerCounter(k=4).count("ATCG"), {"ATCG": 1})

    def test_single_base(self):
        self.assertEqual(KmerCounter(k=1).count("A"), {"A": 1})

    def test_all_N(self):
        self.assertEqual(KmerCounter(k=3).count("NNNNN"), {})

    def test_N_in_middle(self):
        # "ATCNATCG" -> ATC (from first half), then N resets, then ATC TCG (second half)
        counts = KmerCounter(k=3).count("ATCNATCG")
        self.assertEqual(counts, {"ATC": 2, "TCG": 1})

    def test_N_at_start(self):
        counts = KmerCounter(k=2).count("NNATCG")
        self.assertEqual(counts, {"AT": 1, "TC": 1, "CG": 1})

    def test_N_at_end(self):
        counts = KmerCounter(k=2).count("ATCGNN")
        self.assertEqual(counts, {"AT": 1, "TC": 1, "CG": 1})

    def test_invalid_k_zero(self):
        with self.assertRaises(ValueError):
            KmerCounter(k=0)

    def test_invalid_k_33(self):
        with self.assertRaises(ValueError):
            KmerCounter(k=33)

    def test_k1_all_bases(self):
        self.assertEqual(
            KmerCounter(k=1).count("ACGT"), {"A": 1, "C": 1, "G": 1, "T": 1}
        )

    def test_even_k_palindrome(self):
        # "ATAT" is its own reverse complement (even-k palindrome).
        # Non-canonical: count it as-is.
        self.assertEqual(KmerCounter(k=4).count("ATAT"), {"ATAT": 1})

    def test_even_k_palindrome_catg(self):
        # "CATG" is its own reverse complement.
        self.assertEqual(KmerCounter(k=4).count("CATG"), {"CATG": 1})


# ===========================================================================
# 4.  Canonical mode
# ===========================================================================
class TestKmerCounterCanonical(unittest.TestCase):
    """Canonical mode collapses reverse-complement pairs correctly."""

    def test_simple_pair_collapse(self):
        # k=2: AC and GT are reverse complements.  canonical = min("AC","GT") = "AC".
        counts = KmerCounter(k=2, canonical=True).count("ACGT")
        # AC -> AC, CG -> CG (palindrome), GT -> AC
        self.assertEqual(counts, {"AC": 2, "CG": 1})

    def test_palindrome_counted_once(self):
        # CG is a k=2 palindrome (rc("CG") == "CG").
        counts = KmerCounter(k=2, canonical=True).count("CG")
        self.assertEqual(counts, {"CG": 1})

    def test_even_k_palindrome_canonical(self):
        # ATAT: rc = ATAT.  canonical("ATAT") == "ATAT".
        counts = KmerCounter(k=4, canonical=True).count("ATAT")
        self.assertEqual(counts, {"ATAT": 1})

    def test_canonical_sum_equals_noncanonical_sum(self):
        seq = "ATCGATCGATCG"
        k = 3
        non_canon = KmerCounter(k=k, canonical=False).count(seq)
        canon = KmerCounter(k=k, canonical=True).count(seq)
        # Total number of k-mers emitted must be the same.
        self.assertEqual(sum(non_canon.values()), sum(canon.values()))

    def test_forward_and_reverse_merged(self):
        # "AT" and its rc "AT" — wait, rc(AT)=AT (palindrome for k=2).
        # Use a non-palindromic example: "AAC" has rc "GTT".
        # min("AAC","GTT") = "AAC".
        # Sequence "AACGTT": AAC ACG CGT GTT
        # canonical: AAC -> AAC, ACG -> ACG (rc=CGT, min=ACG), CGT -> ACG, GTT -> AAC
        counts = KmerCounter(k=3, canonical=True).count("AACGTT")
        self.assertEqual(counts, {"AAC": 2, "ACG": 2})


# ===========================================================================
# 5.  Rust vs Python parity
# ===========================================================================
@unittest.skipUnless(HAS_RUST, "Rust extension not available")
class TestKmerCounterRustVsPython(unittest.TestCase):
    """Rust back-end produces identical output to the Python fallback."""

    def _python_count(self, sequences, k, canonical):
        """Count using only the Python fallback (bypass Rust)."""
        from Bio.SeqUtils._kmer_counter import _normalise_input

        byte_list = _normalise_input(sequences)
        merged = {}
        for seq_bytes in byte_list:
            for kmer, cnt in _count_kmers_python(seq_bytes, k, canonical).items():
                merged[kmer] = merged.get(kmer, 0) + cnt
        return merged

    def test_random_sequences(self):
        rng = random.Random(12345)
        bases = "ACGT"
        alphabet_with_n = "ACGTN"
        seqs = []
        for _ in range(100):
            length = rng.randint(1, 500)
            seq = "".join(rng.choices(alphabet_with_n, k=length))
            seqs.append(seq)
        for k in (3, 7, 15):
            for canonical in (False, True):
                rust_result = KmerCounter(k=k, canonical=canonical).count(seqs)
                py_result = self._python_count(seqs, k, canonical)
                self.assertEqual(
                    rust_result,
                    py_result,
                    f"Mismatch at k={k}, canonical={canonical}",
                )


# ===========================================================================
# 6.  Jellyfish accuracy + timing
# ===========================================================================
@unittest.skipUnless(HAS_JELLYFISH, "jellyfish not installed")
class TestKmerCounterJellyfish(unittest.TestCase):
    """Accuracy: output matches jellyfish bit-for-bit.  Timing: both are reported."""

    @classmethod
    def setUpClass(cls):
        # Write a deterministic FASTA to a temp file.
        rng = random.Random(42)
        bases = "ACGT"
        cls.tmpdir = tempfile.mkdtemp()
        cls.fasta_path = os.path.join(cls.tmpdir, "bench.fa")
        with open(cls.fasta_path, "w") as fh:
            for i in range(1000):
                seq = "".join(rng.choices(bases, k=150))
                fh.write(f">read_{i}\n{seq}\n")

    @classmethod
    def tearDownClass(cls):
        import shutil as _shutil

        _shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def _run_jellyfish(self, k):
        jf_path = os.path.join(self.tmpdir, "out.jf")
        subprocess.run(
            [
                "jellyfish",
                "count",
                "-m",
                str(k),
                "-s",
                "100M",
                "-t",
                "4",
                "-C",
                self.fasta_path,
                "-o",
                jf_path,
            ],
            check=True,
            capture_output=True,
        )
        dump = subprocess.run(
            ["jellyfish", "dump", "-c", jf_path], check=True, capture_output=True, text=True
        )
        counts = {}
        for line in dump.stdout.strip().splitlines():
            kmer, cnt = line.split()
            counts[kmer] = int(cnt)
        return counts

    def test_accuracy_k21(self):
        t0 = time.perf_counter()
        jf_counts = self._run_jellyfish(21)
        jf_time = time.perf_counter() - t0

        counter = KmerCounter(k=21, canonical=True)
        t0 = time.perf_counter()
        bp_counts = counter.count_file(self.fasta_path)
        bp_time = time.perf_counter() - t0

        print(
            f"\n  [Jellyfish k=21] jellyfish={jf_time:.3f}s  "
            f"KmerCounter={bp_time:.3f}s  ratio={jf_time / max(bp_time, 1e-9):.2f}x"
        )
        self.assertEqual(bp_counts, jf_counts)

    def test_accuracy_k5(self):
        jf_counts = self._run_jellyfish(5)
        bp_counts = KmerCounter(k=5, canonical=True).count_file(self.fasta_path)
        self.assertEqual(bp_counts, jf_counts)


# ===========================================================================
# 7.  Benchmark — Rust vs Python speedup & jellyfish throughput
# ===========================================================================
@unittest.skipUnless(HAS_RUST, "Rust extension not available")
class TestKmerCounterBenchmark(unittest.TestCase):
    """Timing comparison.  No minimum speedup asserted — ratios are env-dependent."""

    @classmethod
    def setUpClass(cls):
        # 10 000 reads x 150 bp, fixed seed for reproducibility.
        rng = random.Random(42)
        bases = "ACGT"
        cls.sequences = ["".join(rng.choices(bases, k=150)) for _ in range(10_000)]

    def _python_count(self, sequences, k, canonical):
        from Bio.SeqUtils._kmer_counter import _normalise_input

        byte_list = _normalise_input(sequences)
        merged = {}
        for seq_bytes in byte_list:
            for kmer, cnt in _count_kmers_python(seq_bytes, k, canonical).items():
                merged[kmer] = merged.get(kmer, 0) + cnt
        return merged

    def test_rust_vs_python_speedup(self):
        k, canonical = 21, True

        # Rust path
        counter = KmerCounter(k=k, canonical=canonical)
        t0 = time.perf_counter()
        rust_counts = counter.count(self.sequences)
        rust_time = time.perf_counter() - t0

        # Python fallback
        t0 = time.perf_counter()
        py_counts = self._python_count(self.sequences, k, canonical)
        py_time = time.perf_counter() - t0

        print(
            f"\n  [Benchmark k=21] Rust={rust_time:.3f}s  Python={py_time:.3f}s  "
            f"speedup={py_time / max(rust_time, 1e-9):.1f}x"
        )
        # Outputs must be identical — this IS asserted.
        self.assertEqual(rust_counts, py_counts)

    @unittest.skipUnless(HAS_JELLYFISH, "jellyfish not installed")
    def test_jellyfish_throughput(self):
        # Write sequences to a temp FASTA, run jellyfish, compare timing.
        tmpdir = tempfile.mkdtemp()
        fasta_path = os.path.join(tmpdir, "bench.fa")
        try:
            with open(fasta_path, "w") as fh:
                for i, seq in enumerate(self.sequences):
                    fh.write(f">read_{i}\n{seq}\n")

            # Jellyfish
            jf_path = os.path.join(tmpdir, "out.jf")
            t0 = time.perf_counter()
            subprocess.run(
                [
                    "jellyfish",
                    "count",
                    "-m",
                    "21",
                    "-s",
                    "100M",
                    "-t",
                    "4",
                    "-C",
                    fasta_path,
                    "-o",
                    jf_path,
                ],
                check=True,
                capture_output=True,
            )
            subprocess.run(
                ["jellyfish", "dump", "-c", jf_path],
                check=True,
                capture_output=True,
            )
            jf_time = time.perf_counter() - t0

            # KmerCounter on in-memory sequences (excludes file I/O for fairness)
            counter = KmerCounter(k=21, canonical=True)
            t0 = time.perf_counter()
            counter.count(self.sequences)
            bp_time = time.perf_counter() - t0

            print(
                f"\n  [Jellyfish throughput] jellyfish={jf_time:.3f}s  "
                f"KmerCounter={bp_time:.3f}s  ratio={jf_time / max(bp_time, 1e-9):.2f}x"
            )
        finally:
            import shutil as _shutil

            _shutil.rmtree(tmpdir, ignore_errors=True)


# ===========================================================================
# 8.  File streaming
# ===========================================================================
class TestKmerCounterFile(unittest.TestCase):
    """count_file() on a multi-record FASTA matches count() on the same records."""

    def test_multi_record_fasta(self):
        records = ["ATCGATCG", "GCTAGCTA", "ATCGATCG"]
        tmpdir = tempfile.mkdtemp()
        fasta_path = os.path.join(tmpdir, "test.fa")
        try:
            with open(fasta_path, "w") as fh:
                for i, seq in enumerate(records):
                    fh.write(f">seq{i}\n{seq}\n")

            counter = KmerCounter(k=3)
            file_counts = counter.count_file(fasta_path)
            mem_counts = counter.count(records)
            self.assertEqual(file_counts, mem_counts)
        finally:
            import shutil as _shutil

            _shutil.rmtree(tmpdir, ignore_errors=True)

    def test_fastq_format(self):
        tmpdir = tempfile.mkdtemp()
        fastq_path = os.path.join(tmpdir, "test.fq")
        try:
            with open(fastq_path, "w") as fh:
                fh.write("@read1\nATCGATCG\n+\nIIIIIIII\n")
                fh.write("@read2\nGCTAGCTA\n+\nIIIIIIII\n")

            counter = KmerCounter(k=3)
            file_counts = counter.count_file(fastq_path, format="fastq")
            mem_counts = counter.count(["ATCGATCG", "GCTAGCTA"])
            self.assertEqual(file_counts, mem_counts)
        finally:
            import shutil as _shutil

            _shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
