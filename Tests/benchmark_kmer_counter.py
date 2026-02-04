"""Benchmark k-mer counter performance and validate against Jellyfish.

This script performs both performance testing and accuracy verification.
Jellyfish validation is REQUIRED to ensure correctness.
"""

import subprocess
import sys
import tempfile
import time
from collections import Counter

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.KmerCounter import KmerCounter


def generate_test_sequence(length=1000000, seed=42):
    """Generate random DNA sequence with fixed seed for reproducibility."""
    import random

    random.seed(seed)
    return "".join(random.choice("ACGT") for _ in range(length))


def benchmark_python_counter(seq, k):
    """Benchmark collections.Counter."""
    start = time.perf_counter()
    kmers = [seq[i : i + k] for i in range(len(seq) - k + 1)]
    counts = Counter(kmers)
    elapsed = time.perf_counter() - start
    return elapsed, dict(counts)


def benchmark_kmer_counter(seq, k, use_rust):
    """Benchmark KmerCounter."""
    start = time.perf_counter()
    counter = KmerCounter(k=k, use_rust=use_rust)
    counts = counter.count(seq)
    elapsed = time.perf_counter() - start
    return elapsed, counts


def run_jellyfish(fasta_file, k, output_prefix):
    """Run Jellyfish and return counts.

    Returns
    -------
    dict or None
        Dictionary of k-mer counts, or None if Jellyfish fails
    """
    try:
        # Run jellyfish count
        count_cmd = [
            "jellyfish",
            "count",
            "-m",
            str(k),
            "-s",
            "100M",
            "-o",
            output_prefix + ".jf",
            fasta_file,
        ]
        result = subprocess.run(count_cmd, check=True, capture_output=True, text=True)

        # Dump counts
        dump_cmd = ["jellyfish", "dump", output_prefix + ".jf"]
        result = subprocess.run(dump_cmd, check=True, capture_output=True, text=True)

        # Parse output
        counts = {}
        lines = result.stdout.strip().split("\n")
        for i in range(0, len(lines), 2):
            if i + 1 < len(lines):
                count = int(lines[i].lstrip(">"))
                kmer = lines[i + 1]
                counts[kmer] = count

        return counts
    except subprocess.CalledProcessError as e:
        print(f"Jellyfish error: {e}")
        if e.stderr:
            print(f"  stderr: {e.stderr}")
        return None
    except FileNotFoundError:
        print(
            "Jellyfish not found. Install from: http://www.genome.umd.edu/jellyfish.html"
        )
        print("Or on macOS: brew install jellyfish")
        return None


def compare_counts(counts1, counts2, label1, label2, max_show=5):
    """Compare two count dictionaries and report differences.

    Returns
    -------
    bool
        True if counts match exactly, False otherwise
    """
    if counts1 == counts2:
        return True

    print(f"✗ Mismatch between {label1} and {label2}!")
    print(f"  {label1}: {len(counts1)} unique k-mers")
    print(f"  {label2}: {len(counts2)} unique k-mers")

    # Find differences
    all_kmers = set(counts1.keys()) | set(counts2.keys())
    mismatches = []
    for kmer in all_kmers:
        c1 = counts1.get(kmer, 0)
        c2 = counts2.get(kmer, 0)
        if c1 != c2:
            mismatches.append((kmer, c1, c2))

    if mismatches:
        print(f"  Total mismatches: {len(mismatches)}")
        for kmer, c1, c2 in mismatches[:max_show]:
            print(f"    {kmer}: {label1}={c1}, {label2}={c2}")
        if len(mismatches) > max_show:
            print(f"    ... and {len(mismatches) - max_show} more")

    return False


def main():
    """Run benchmarks and validation."""
    # Configuration
    seq_length = 1000000
    k = 21
    seed = 42

    print("K-mer Counter Benchmark and Validation")
    print("=" * 70)
    print(f"Sequence length: {seq_length:,} nucleotides")
    print(f"K-mer size: {k}")
    print(f"Random seed: {seed}")
    print()

    # Generate test sequence
    print("Generating test sequence...")
    seq = generate_test_sequence(seq_length, seed)
    print(f"Generated {len(seq):,} nucleotides")
    print()

    # Performance Benchmarks
    print("=" * 70)
    print("PERFORMANCE BENCHMARKS")
    print("=" * 70)

    # collections.Counter
    print(f"Running collections.Counter (k={k})...")
    time_counter, counts_counter = benchmark_python_counter(seq, k)
    print(
        f"  collections.Counter:   {time_counter:.3f}s ({len(counts_counter):,} unique k-mers)"
    )

    # KmerCounter (Python)
    print(f"Running KmerCounter (Python fallback, k={k})...")
    time_kc_python, counts_kc_python = benchmark_kmer_counter(seq, k, use_rust=False)
    print(
        f"  KmerCounter (Python):  {time_kc_python:.3f}s ({len(counts_kc_python):,} unique)"
    )

    # KmerCounter (Rust)
    print(f"Running KmerCounter (Rust, k={k})...")
    counts_rust = None
    try:
        time_kc_rust, counts_rust = benchmark_kmer_counter(seq, k, use_rust=True)
        print(
            f"  KmerCounter (Rust):    {time_kc_rust:.3f}s ({len(counts_rust):,} unique)"
        )

        # Calculate speedups
        speedup_counter = time_counter / time_kc_rust
        speedup_python = time_kc_python / time_kc_rust

        print()
        print(f"Speedup vs Counter:    {speedup_counter:.1f}x")
        print(f"Speedup vs Python KmerCounter: {speedup_python:.1f}x")
    except Exception as e:
        print(f"  Rust implementation not available: {e}")

    # Jellyfish Validation (REQUIRED)
    print()
    print("=" * 70)
    print("JELLYFISH VALIDATION (REQUIRED)")
    print("=" * 70)

    with tempfile.TemporaryDirectory() as tmpdir:
        # Write test FASTA
        fasta_path = f"{tmpdir}/test.fasta"
        rec = SeqRecord(Seq(seq), id="test", description="")
        SeqIO.write([rec], fasta_path, "fasta")
        print(f"Created test FASTA: {fasta_path}")

        # Run Jellyfish
        print(f"\nRunning Jellyfish (k={k})...")
        start_jf = time.perf_counter()
        counts_jf = run_jellyfish(fasta_path, k, f"{tmpdir}/test")
        time_jf = time.perf_counter() - start_jf

        if counts_jf is None:
            print()
            print("=" * 70)
            print("⚠ CRITICAL: Jellyfish validation FAILED")
            print("=" * 70)
            print("Jellyfish is REQUIRED for accuracy verification.")
            print("Without Jellyfish validation, we cannot guarantee correctness.")
            print()
            print("Installation:")
            print("  - macOS: brew install jellyfish")
            print("  - Linux: apt-get install jellyfish (or compile from source)")
            print("  - Source: http://www.genome.umd.edu/jellyfish.html")
            return 1

        print(f"  Jellyfish:             {time_jf:.3f}s ({len(counts_jf):,} unique k-mers)")

        # Accuracy Verification
        print()
        print("=" * 70)
        print("ACCURACY VERIFICATION")
        print("=" * 70)

        all_pass = True

        # Verify Rust vs Jellyfish
        if counts_rust:
            print("\n1. Rust implementation vs Jellyfish:")
            if compare_counts(counts_rust, counts_jf, "Rust", "Jellyfish"):
                print("   ✓ PASS: Rust implementation matches Jellyfish exactly")
            else:
                all_pass = False

        # Verify Python vs Jellyfish
        print("\n2. Python implementation vs Jellyfish:")
        if compare_counts(counts_kc_python, counts_jf, "Python", "Jellyfish"):
            print("   ✓ PASS: Python implementation matches Jellyfish exactly")
        else:
            all_pass = False

        # Verify Rust vs Python
        if counts_rust:
            print("\n3. Rust vs Python consistency:")
            if compare_counts(counts_rust, counts_kc_python, "Rust", "Python"):
                print("   ✓ PASS: Rust and Python implementations match")
            else:
                all_pass = False

        # Final verdict
        print()
        print("=" * 70)
        if all_pass:
            print("✓ ALL VALIDATIONS PASSED")
            print("=" * 70)
            print("Both implementations produce correct results verified by Jellyfish.")
            return 0
        else:
            print("✗ VALIDATION FAILED")
            print("=" * 70)
            print("Implementation does not match Jellyfish reference.")
            print("This indicates a bug that must be fixed.")
            return 1


if __name__ == "__main__":
    sys.exit(main())
