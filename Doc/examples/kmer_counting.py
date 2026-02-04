# Copyright 2025 by the Biopython Contributors.
# All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Example: counting k-mers with Bio.SeqUtils.KmerCounter.

Demonstrates the three main entry points:
  1. Counting k-mers in a single string.
  2. Canonical mode — reverse-complement pairs are merged.
  3. Streaming a FASTA / FASTQ file without loading it fully into memory.
"""

from Bio.SeqUtils import KmerCounter

# --- 1. Count k-mers in a single string --------------------------------
counter = KmerCounter(k=5)
counts = counter.count("ATCGATCGATCGATCG")
print("=== Non-canonical k=5 ===")
for kmer, n in sorted(counts.items()):
    print(f"  {kmer}\t{n}")

# --- 2. Canonical mode (reverse-complement pairs merged) ----------------
counter_c = KmerCounter(k=5, canonical=True)
counts_c = counter_c.count("ATCGATCGATCGATCG")
print("\n=== Canonical k=5 ===")
for kmer, n in sorted(counts_c.items()):
    print(f"  {kmer}\t{n}")

# --- 3. Streaming a FASTA / FASTQ file ----------------------------------
# Uncomment and supply a real file path to try file streaming:
#
# counts_file = counter.count_file("reads.fastq", format="fastq")
# print("\n=== File counts ===")
# for kmer, n in sorted(counts_file.items()):
#     print(f"  {kmer}\t{n}")
