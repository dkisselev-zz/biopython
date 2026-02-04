#!/usr/bin/env python
"""Examples of using Bio.SeqUtils.KmerCounter for k-mer counting.

This script demonstrates various use cases for the KmerCounter module,
including basic counting, canonical k-mers, streaming large files, and
working with different sequence formats.
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.KmerCounter import KmerCounter, count_kmers, streaming_count_kmers

print("K-mer Counting Examples")
print("=" * 70)

# Example 1: Basic k-mer counting
print("\n1. Basic k-mer counting")
print("-" * 70)

sequence = "ACGTACGTACGTACGT"
k = 3
counts = count_kmers(sequence, k=k)

print(f"Sequence: {sequence}")
print(f"K-mer size: {k}")
print(f"K-mer counts:")
for kmer in sorted(counts.keys()):
    print(f"  {kmer}: {counts[kmer]}")

# Example 2: Canonical k-mers
print("\n2. Canonical k-mer counting")
print("-" * 70)

sequence = "ACGTACGT"
k = 3

# Without canonical (k-mer and RC are different)
counter_normal = KmerCounter(k=k, canonical=False)
counts_normal = counter_normal.count(sequence)

# With canonical (k-mer and RC are same)
counter_canonical = KmerCounter(k=k, canonical=True)
counts_canonical = counter_canonical.count(sequence)

print(f"Sequence: {sequence}")
print(f"K-mer size: {k}")
print()
print("Without canonical mode:")
for kmer in sorted(counts_normal.keys()):
    print(f"  {kmer}: {counts_normal[kmer]}")

print()
print("With canonical mode:")
for kmer in sorted(counts_canonical.keys()):
    print(f"  {kmer}: {counts_canonical[kmer]}")
print()
print("Note: ACG and CGT are reverse complements, counted together in canonical mode")

# Example 3: Working with Bio.Seq objects
print("\n3. Working with Bio.Seq objects")
print("-" * 70)

seq_obj = Seq("ATCGATCGATCG")
counts = count_kmers(seq_obj, k=4)

print(f"Seq object: {seq_obj}")
print(f"K-mer counts (k=4):")
for kmer in sorted(counts.keys()):
    print(f"  {kmer}: {counts[kmer]}")

# Example 4: Multiple sequences
print("\n4. Counting across multiple sequences")
print("-" * 70)

sequences = [
    "ACGTACGT",
    "TGCATGCA",
    "AAAACCCC",
]

counter = KmerCounter(k=3)
counts = counter.count(sequences)

print("Sequences:")
for i, seq in enumerate(sequences, 1):
    print(f"  {i}. {seq}")

print(f"\nCombined k-mer counts (k=3):")
for kmer in sorted(counts.keys()):
    print(f"  {kmer}: {counts[kmer]}")

# Example 5: SeqRecord objects
print("\n5. Working with SeqRecord objects")
print("-" * 70)

records = [
    SeqRecord(Seq("ACGTACGT"), id="seq1", description="Test sequence 1"),
    SeqRecord(Seq("TGCATGCA"), id="seq2", description="Test sequence 2"),
]

counter = KmerCounter(k=4)
counts = counter.count(records)

print("SeqRecords:")
for record in records:
    print(f"  {record.id}: {record.seq}")

print(f"\nK-mer counts (k=4):")
for kmer in sorted(counts.keys()):
    print(f"  {kmer}: {counts[kmer]}")

# Example 6: Handling ambiguous bases
print("\n6. Handling ambiguous bases")
print("-" * 70)

# Sequence with N
sequence_with_n = "ACNGTACGT"
counts = count_kmers(sequence_with_n, k=3)

print(f"Sequence with N: {sequence_with_n}")
print(f"K-mer counts (k=3): {counts}")
print("Note: K-mers containing N are excluded")

# Sequence without N
sequence_clean = "ACGTACGT"
counts = count_kmers(sequence_clean, k=3)
print(f"\nSequence without N: {sequence_clean}")
print(f"K-mer counts (k=3): {len(counts)} unique k-mers")

# Example 7: Performance comparison
print("\n7. Performance comparison (Rust vs Python)")
print("-" * 70)

# Generate a longer sequence for performance test
import random

random.seed(42)
long_sequence = "".join(random.choice("ACGT") for _ in range(100000))

# Time Rust implementation
import time

counter_rust = KmerCounter(k=21, use_rust=True)
start = time.perf_counter()
counts_rust = counter_rust.count(long_sequence)
time_rust = time.perf_counter() - start

# Time Python implementation
counter_python = KmerCounter(k=21, use_rust=False)
start = time.perf_counter()
counts_python = counter_python.count(long_sequence)
time_python = time.perf_counter() - start

print(f"Sequence length: {len(long_sequence):,} nucleotides")
print(f"K-mer size: 21")
print(f"\nRust implementation:")
print(f"  Time: {time_rust:.3f}s")
print(f"  Unique k-mers: {len(counts_rust):,}")
print(f"  Using Rust: {counter_rust.using_rust}")

print(f"\nPython implementation:")
print(f"  Time: {time_python:.3f}s")
print(f"  Unique k-mers: {len(counts_python):,}")

if time_python > 0:
    speedup = time_python / time_rust if time_rust > 0 else 0
    print(f"\nSpeedup: {speedup:.1f}x")

# Verify results match
if counts_rust == counts_python:
    print("✓ Results match exactly")
else:
    print("✗ Results differ!")

# Example 8: Streaming large files (demonstration)
print("\n8. Streaming large FASTA files")
print("-" * 70)

# Create a temporary test file
import tempfile

with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
    temp_fasta = f.name
    # Write some test sequences
    f.write(">seq1\n")
    f.write("ACGTACGTACGTACGT\n")
    f.write(">seq2\n")
    f.write("TGCATGCATGCATGCA\n")
    f.write(">seq3\n")
    f.write("AAAACCCCGGGGTTTT\n")

# Count k-mers using streaming
counts = streaming_count_kmers(temp_fasta, k=5, format="fasta")

print(f"Counted k-mers from temporary FASTA file")
print(f"K-mer size: 5")
print(f"Unique k-mers: {len(counts)}")
print(f"\nTop 10 most common k-mers:")
for kmer, count in sorted(counts.items(), key=lambda x: x[1], reverse=True)[:10]:
    print(f"  {kmer}: {count}")

# Clean up
import os

os.unlink(temp_fasta)

print("\n" + "=" * 70)
print("Examples completed successfully!")
