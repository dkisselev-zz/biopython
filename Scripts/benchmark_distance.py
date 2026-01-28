#!/usr/bin/env python
# Copyright 2024 The Biopython Contributors.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Benchmark distance matrix computation methods.

Generates synthetic alignments at configurable size profiles, runs each
registered method, and reports wall-clock time, peak memory, and RMSE
relative to the legacy pure-Python implementation.

Usage::

    python Scripts/benchmark_distance.py --profile amplicon --method all
    python Scripts/benchmark_distance.py --profile viral --method numpy --model identity
    python Scripts/benchmark_distance.py --profile metagenomic --method all --skip-ref-above 500 --output-format csv

Predefined profiles:
 - viral:        100 sequences x 30,000 columns  (DNA)
 - amplicon:     1,000 sequences x 1,500 columns (DNA)
 - metagenomic:  10,000 sequences x 300 columns  (DNA)

Output formats:
 - text (default): human-readable table
 - csv: comma-separated for CI pipelines and result archiving
"""

import argparse
import csv
import math
import random
import sys
import time
import tracemalloc
import warnings

from Bio.Align import MultipleSeqAlignment
from Bio import BiopythonExperimentalWarning
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

PROFILES = {
    "viral": (100, 30000),
    "amplicon": (1000, 1500),
    "metagenomic": (10000, 300),
}

DNA_ALPHABET = "ACGT"
PROTEIN_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"


def _generate_synthetic_msa(n_seqs, n_cols, alphabet=DNA_ALPHABET, seed=42):
    """Generate a synthetic MSA as a MultipleSeqAlignment.

    Arguments:
     - n_seqs - Number of sequences.
     - n_cols - Alignment length (columns).
     - alphabet - Character set for the consensus and mutations.
       Defaults to DNA; pass ``PROTEIN_ALPHABET`` for protein benchmarks.
     - seed - Random seed for reproducibility.

    A consensus sequence is generated first; each derived sequence
    accumulates point mutations at a 15 % rate.
    """
    rng = random.Random(seed)
    consensus = "".join(rng.choice(alphabet) for _ in range(n_cols))
    records = []
    for i in range(n_seqs):
        seq = list(consensus)
        for col in range(n_cols):
            if rng.random() < 0.15:
                seq[col] = rng.choice(alphabet)
        records.append(SeqRecord(Seq("".join(seq)), id=f"Seq{i:06d}"))
    return MultipleSeqAlignment(records)


def _run_benchmark(msa, method, model="identity"):
    """Run one benchmark iteration, returning (time_sec, peak_mem_mb, dm).

    Arguments:
     - msa - A MultipleSeqAlignment to compute distances for.
     - method - Method name (e.g. ``"python"``, ``"numpy"``).
     - model - Substitution model name.
    """
    tracemalloc.start()
    start = time.perf_counter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonExperimentalWarning)
        calc = DistanceCalculator(model, method=method, n_jobs=None)
    dm = calc.get_distance(msa)
    elapsed = time.perf_counter() - start
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return elapsed, peak / (1024 * 1024), dm


def _compute_rmse(dm_ref, dm_test):
    """Compute RMSE between two DistanceMatrix objects.

    Arguments:
     - dm_ref - Reference DistanceMatrix (legacy).
     - dm_test - Test DistanceMatrix (accelerated method).
    """
    n = len(dm_ref)
    total = 0.0
    count = 0
    for i in range(n):
        for j in range(i):
            diff = dm_ref[i, j] - dm_test[i, j]
            total += diff * diff
            count += 1
    return math.sqrt(total / count) if count else 0.0


def main():
    """Entry point for the benchmark CLI."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--profile",
        choices=PROFILES,
        default="amplicon",
        help="Synthetic alignment size profile.",
    )
    parser.add_argument(
        "--method",
        default="all",
        help="Method to benchmark, or 'all' for full comparison.",
    )
    parser.add_argument(
        "--model",
        default="identity",
        help="Substitution model (e.g. identity, blosum62).",
    )
    parser.add_argument(
        "--alphabet",
        choices=["dna", "protein"],
        default="dna",
        help="Alphabet for synthetic MSA generation.",
    )
    parser.add_argument(
        "--skip-ref-above",
        type=int,
        default=500,
        help="Skip the legacy reference run when sequence count exceeds this threshold.",
    )
    parser.add_argument(
        "--output-format",
        choices=["text", "csv"],
        default="text",
        help="Output format for results.",
    )
    args = parser.parse_args()

    n_seqs, n_cols = PROFILES[args.profile]
    alphabet = PROTEIN_ALPHABET if args.alphabet == "protein" else DNA_ALPHABET
    msa = _generate_synthetic_msa(n_seqs, n_cols, alphabet=alphabet)

    methods = (
        ["python", "numpy", "scipy", "onehot"]
        if args.method == "all"
        else [args.method]
    )

    # Determine whether to compute legacy reference
    compute_ref = (
        "python" in methods or args.method == "all"
    ) and n_seqs <= args.skip_ref_above
    ref_dm = None

    rows = []  # list of dicts for CSV output

    if compute_ref:
        t, mem, ref_dm = _run_benchmark(msa, "python", args.model)
        rows.append(
            {
                "method": "python",
                "time_sec": f"{t:.3f}",
                "mem_mb": f"{mem:.2f}",
                "rmse": "N/A (reference)",
            }
        )

    for method in methods:
        if method == "python" and compute_ref:
            continue  # already ran
        if method == "python" and not compute_ref:
            rows.append(
                {
                    "method": "python",
                    "time_sec": "SKIPPED",
                    "mem_mb": "SKIPPED",
                    "rmse": f"N > {args.skip_ref_above}",
                }
            )
            continue
        try:
            t, mem, dm = _run_benchmark(msa, method, args.model)
            rmse = f"{_compute_rmse(ref_dm, dm):.2e}" if ref_dm else "N/A (ref skipped)"
            rows.append(
                {
                    "method": method,
                    "time_sec": f"{t:.3f}",
                    "mem_mb": f"{mem:.2f}",
                    "rmse": rmse,
                }
            )
        except Exception as exc:
            rows.append(
                {
                    "method": method,
                    "time_sec": "FAILED",
                    "mem_mb": "FAILED",
                    "rmse": str(exc),
                }
            )

    # Output
    header = (
        f"Profile: {args.profile} ({n_seqs} seqs x {n_cols} cols), "
        f"model={args.model}, alphabet={args.alphabet}"
    )
    if args.output_format == "csv":
        writer = csv.DictWriter(
            sys.stdout, fieldnames=["method", "time_sec", "mem_mb", "rmse"]
        )
        writer.writeheader()
        writer.writerows(rows)
    else:
        print(header)
        print(f"  {'method':<10} {'time(s)':>9} {'mem(MB)':>9} {'rmse'}")
        print(f"  {'-' * 10} {'-' * 9} {'-' * 9} {'-' * 20}")
        for r in rows:
            print(
                f"  {r['method']:<10} {r['time_sec']:>9} {r['mem_mb']:>9} {r['rmse']}"
            )


if __name__ == "__main__":
    main()
