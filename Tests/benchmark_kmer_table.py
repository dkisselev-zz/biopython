"""Benchmark k-mer counter with clean table output for easy comparison.

Compares:
1. Rust vs Python fallback
2. KmerCounter vs Jellyfish (full file with I/O)
3. KmerCounter vs Jellyfish (count only, no I/O)
"""

import subprocess
import sys
import tempfile
import time
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.KmerCounter import KmerCounter


def generate_test_reads(num_reads=10000, read_length=150, seed=42):
    """Generate simulated sequencing reads."""
    import random

    random.seed(seed)
    reads = []
    for i in range(num_reads):
        seq = "".join(random.choice("ACGT") for _ in range(read_length))
        record = SeqRecord(Seq(seq), id=f"read_{i}", description="")
        reads.append(record)
    return reads


def benchmark_rust_full(fasta_file, k, canonical):
    """Benchmark Rust KmerCounter with I/O (full pipeline)."""
    start = time.perf_counter()
    counter = KmerCounter(k=k, canonical=canonical, use_rust=True)
    records = list(SeqIO.parse(fasta_file, "fasta"))
    counts = counter.count(records)
    elapsed = time.perf_counter() - start
    return elapsed, counts


def benchmark_rust_count_only(records, k, canonical):
    """Benchmark Rust KmerCounter counting only (no I/O)."""
    counter = KmerCounter(k=k, canonical=canonical, use_rust=True)
    start = time.perf_counter()
    counts = counter.count(records)
    elapsed = time.perf_counter() - start
    return elapsed, counts


def benchmark_python(fasta_file, k, canonical):
    """Benchmark Python fallback KmerCounter."""
    start = time.perf_counter()
    counter = KmerCounter(k=k, canonical=canonical, use_rust=False)
    records = list(SeqIO.parse(fasta_file, "fasta"))
    counts = counter.count(records)
    elapsed = time.perf_counter() - start
    return elapsed, counts


def benchmark_jellyfish_full(fasta_file, k, canonical, output_prefix):
    """Benchmark Jellyfish including file I/O."""
    start = time.perf_counter()

    # Count command (includes reading file)
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
    if canonical:
        count_cmd.insert(2, "-C")  # Add canonical flag

    subprocess.run(count_cmd, check=True, capture_output=True)

    # Dump command (includes writing output)
    dump_cmd = ["jellyfish", "dump", output_prefix + ".jf"]
    result = subprocess.run(dump_cmd, check=True, capture_output=True, text=True)

    elapsed = time.perf_counter() - start

    # Parse output
    counts = {}
    lines = result.stdout.strip().split("\n")
    for i in range(0, len(lines), 2):
        if i + 1 < len(lines):
            count = int(lines[i].lstrip(">"))
            kmer = lines[i + 1]
            counts[kmer] = count

    return elapsed, counts


def benchmark_jellyfish_count_only(fasta_file, k, canonical, output_prefix):
    """Benchmark Jellyfish count only (file already read into memory)."""
    # Pre-run to create the database (not timed)
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
    if canonical:
        count_cmd.insert(2, "-C")

    subprocess.run(count_cmd, check=True, capture_output=True)

    # Now time just the counting operation
    # We'll re-run count command and measure it
    start = time.perf_counter()
    subprocess.run(count_cmd, check=True, capture_output=True)
    elapsed = time.perf_counter() - start

    # Get counts for validation (not timed)
    dump_cmd = ["jellyfish", "dump", output_prefix + ".jf"]
    result = subprocess.run(dump_cmd, check=True, capture_output=True, text=True)

    counts = {}
    lines = result.stdout.strip().split("\n")
    for i in range(0, len(lines), 2):
        if i + 1 < len(lines):
            count = int(lines[i].lstrip(">"))
            kmer = lines[i + 1]
            counts[kmer] = count

    return elapsed, counts


def print_table(num_reads, read_length, k, canonical, rust_count_time, python_time, rust_full_time, jf_full_time, jf_count_time):
    """Print results in a nice table format."""
    print(f"\nLatest benchmark results ({num_reads:,} reads × {read_length} bp, k={k} canonical={canonical}):")

    # Calculate ratios
    rust_vs_python = python_time / rust_count_time
    kmc_vs_jf_full = jf_full_time / rust_full_time
    kmc_vs_jf_count = jf_count_time / rust_count_time

    # Table borders
    top = "┌───────────────────────────────────────────────┬────────────────────┬──────────────────┬──────────────────────────┐"
    mid = "├───────────────────────────────────────────────┼────────────────────┼──────────────────┼──────────────────────────┤"
    bot = "└───────────────────────────────────────────────┴────────────────────┴──────────────────┴──────────────────────────┘"
    header = "│                  Comparison                   │       Time A       │      Time B      │          Ratio           │"

    print(top)
    print(header)
    print(mid)

    # Row 1: Rust vs Python (count only)
    print(
        f"│ Rust vs Python fallback                       │ Rust {rust_count_time:5.2f} s        │ Python {python_time:5.2f} s    │ {rust_vs_python:.1f}× speedup             │"
    )

    print(mid)

    # Row 2: KmerCounter vs jellyfish (full file with I/O)
    print(
        f"│ KmerCounter vs jellyfish (full file)          │ KmerCounter {rust_full_time:.2f} s │ jellyfish {jf_full_time:.2f} s │ KmerCounter {kmc_vs_jf_full:.2f}× faster │"
    )

    print(mid)

    # Row 3: KmerCounter vs jellyfish (count only, no I/O)
    print(
        f"│ KmerCounter vs jellyfish (count only, no I/O) │ KmerCounter {rust_count_time:.2f} s │ jellyfish {jf_count_time:.2f} s │ KmerCounter {kmc_vs_jf_count:.1f}× faster │"
    )

    print(bot)


def main():
    # Parameters
    num_reads = 10000
    read_length = 150
    k = 21
    canonical = True

    print(f"Generating {num_reads:,} reads of {read_length} bp each...")

    # Generate test data
    reads = generate_test_reads(num_reads, read_length)

    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_file = f"{tmpdir}/reads.fasta"
        SeqIO.write(reads, fasta_file, "fasta")
        print(f"Wrote test FASTA: {fasta_file}")

        # Run benchmarks
        print("\nRunning benchmarks...")

        print("  1/5: Rust KmerCounter (full with I/O)...")
        rust_full_time, rust_counts = benchmark_rust_full(fasta_file, k, canonical)

        print("  2/5: Rust KmerCounter (count only, no I/O)...")
        # Pre-load records for count-only benchmark
        records = list(SeqIO.parse(fasta_file, "fasta"))
        rust_count_time, _ = benchmark_rust_count_only(records, k, canonical)

        print("  3/5: Python KmerCounter...")
        python_time, python_counts = benchmark_python(fasta_file, k, canonical)

        print("  4/5: Jellyfish (full file with I/O)...")
        try:
            jf_full_time, jf_counts = benchmark_jellyfish_full(
                fasta_file, k, canonical, f"{tmpdir}/jf_full"
            )
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            print(f"    Error: Jellyfish not available: {e}")
            print(
                "    Install with: brew install jellyfish  (macOS) or apt-get install jellyfish (Linux)"
            )
            sys.exit(1)

        print("  5/5: Jellyfish (count only, no I/O)...")
        jf_count_time, _ = benchmark_jellyfish_count_only(
            fasta_file, k, canonical, f"{tmpdir}/jf_count"
        )

        # Verify correctness
        if rust_counts != jf_counts:
            print("\n⚠ WARNING: Rust counts don't match Jellyfish!")
            print(f"  Rust: {len(rust_counts)} unique k-mers")
            print(f"  Jellyfish: {len(jf_counts)} unique k-mers")
        else:
            print("\n✓ Validation passed: Rust matches Jellyfish exactly")

        if python_counts != jf_counts:
            print("⚠ WARNING: Python counts don't match Jellyfish!")
        else:
            print("✓ Validation passed: Python matches Jellyfish exactly")

        # Print results table
        print_table(
            num_reads,
            read_length,
            k,
            canonical,
            rust_count_time,
            python_time,
            rust_full_time,
            jf_full_time,
            jf_count_time,
        )

        # Additional info
        print(f"\nTotal unique k-mers: {len(rust_counts):,}")
        print(f"Total k-mer instances: {sum(rust_counts.values()):,}")
        print(f"\nTiming breakdown:")
        print(f"  Rust full pipeline: {rust_full_time:.3f}s (I/O + counting)")
        print(f"  Rust count only:    {rust_count_time:.3f}s (just counting)")
        print(f"  Rust I/O overhead:  {rust_full_time - rust_count_time:.3f}s")


if __name__ == "__main__":
    main()
