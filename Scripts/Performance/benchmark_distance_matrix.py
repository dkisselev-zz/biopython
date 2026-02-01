#!/usr/bin/env python
"""Benchmark Bio.Phylo distance matrix computation backends.

Two-phase runner
----------------
1. **Correctness** — every method is compared against the ``python`` (legacy)
   backend via RMSE on a validation subset.  The subset is capped at
   ``VALIDATION_N`` sequences so that the legacy baseline finishes quickly
   regardless of profile size.  RMSE is *always* reported; there is no way to
   skip it.

2. **Performance** — each method is timed on the *full* alignment for
   ``--repeats`` iterations using the public ``DistanceCalculator.get_distance``
   API (the same code path a user calls).  ``DistanceCalculator`` construction
   is performed once outside the timed loop so that one-time model-loading
   overhead is excluded.  Wall-clock ``time_min`` and ``time_mean`` plus peak
   memory are reported.  The ``python`` backend is only timed when
   N <= ``PYTHON_TIMING_MAX_N``; above that threshold its timing slot is marked
   ``N/A`` (correctness is still validated on the subset).

Predefined profiles
-------------------
=============  ======  =======  ================================================
Profile          N        L      Typical use-case
=============  ======  =======  ================================================
viral            100   30 000   Large viral genomes (e.g. influenza pan-genome)
amplicon       1 000    1 500   16S / amplicon sequencing
metagenomic   10 000      300   Short metagenomic reads
=============  ======  =======  ================================================

Custom ``NxL`` profiles (e.g. ``500x2000``) are also accepted.

Usage examples
--------------
::

    # Full validation + benchmark — all three baseline profiles, all methods
    python benchmark_distance_matrix.py

    # Viral profile only, blastn model, 5 repeats, CSV output
    python benchmark_distance_matrix.py --profiles viral --model blastn \\
        --repeats 5 --output viral_blastn.csv

    # Single custom profile
    python benchmark_distance_matrix.py --profiles 200x5000
"""

import argparse
import csv
import math
import random
import sys
import time
import tracemalloc

# ---------------------------------------------------------------------------
# Configuration constants
# ---------------------------------------------------------------------------

PROFILES = {
    "viral": {"n": 100, "l": 30000},
    "amplicon": {"n": 1000, "l": 1500},
    "metagenomic": {"n": 10000, "l": 300},
}

METHODS = ["python", "numpy", "scipy", "onehot"]

DNA_ALPHABET = "ACGT"
PROTEIN_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

# Correctness validation is performed on at most this many sequences so that
# the python (legacy) baseline finishes in seconds even for the largest profiles.
VALIDATION_N = 200

# The python backend is timed on the full alignment only when N is at or below
# this value.  Above it, python timing is skipped but correctness is still
# validated on the VALIDATION_N subset.
PYTHON_TIMING_MAX_N = 1000


# ---------------------------------------------------------------------------
# Synthetic alignment generator
# ---------------------------------------------------------------------------


def generate_alignment(n, l, alphabet="ACGT", mutation_rate=0.1, gap_rate=0.0, seed=42):
    """Generate a synthetic MSA of *n* sequences of length *l*.

    A random reference sequence is drawn from *alphabet* using the given
    *seed*.  Each subsequent sequence is an independent per-position mutation
    of that reference.  Gap insertion (``"-"``) is applied after mutation when
    *gap_rate* > 0.

    Parameters
    ----------
    n : int
        Number of sequences.
    l : int
        Alignment length (columns).
    alphabet : str
        Characters available for mutation.
    mutation_rate : float
        Per-position probability of a mutation (0–1).
    gap_rate : float
        Per-position probability of inserting a gap after mutation (0–1).
    seed : int
        Random seed for full reproducibility.

    Returns
    -------
    Bio.Align.MultipleSeqAlignment
        Aligned sequences, all of length *l*.
    """
    from Bio.Align import MultipleSeqAlignment
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    rng = random.Random(seed)
    ref = [rng.choice(alphabet) for _ in range(l)]
    records = [SeqRecord(Seq("".join(ref)), id="ref_0")]
    for idx in range(1, n):
        seq = list(ref)
        for i in range(l):
            if rng.random() < mutation_rate:
                seq[i] = rng.choice(alphabet)
            if gap_rate > 0 and rng.random() < gap_rate:
                seq[i] = "-"
        records.append(SeqRecord(Seq("".join(seq)), id=f"seq_{idx}"))
    return MultipleSeqAlignment(records)


# ---------------------------------------------------------------------------
# Measurement helpers
# ---------------------------------------------------------------------------


def _get_flat_distances(msa, model, method):
    """Return lower-triangular distances via the public DistanceCalculator API.

    No timing or memory instrumentation — used for correctness validation only.
    """
    from Bio.Phylo.TreeConstruction import DistanceCalculator

    dc = DistanceCalculator(model, method=method)
    dm = dc.get_distance(msa)
    names = dm.names
    flat = []
    for i in range(1, len(names)):
        for j in range(i):
            flat.append(dm[names[i], names[j]])
    return flat


def run_timed(dc, msa):
    """Time a single ``get_distance`` call on a pre-constructed calculator.

    Parameters
    ----------
    dc : DistanceCalculator
        Already-constructed calculator (model loading excluded from timing).
    msa : MultipleSeqAlignment
        The alignment to compute distances for.

    Returns
    -------
    time_s : float
        Wall-clock seconds for ``get_distance``.
    peak_mb : float
        Peak memory (tracemalloc) in megabytes.
    """
    tracemalloc.start()
    t0 = time.perf_counter()
    dc.get_distance(msa)
    t1 = time.perf_counter()
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return t1 - t0, peak / (1024 * 1024)


def compute_rmse(flat_a, flat_b):
    """Root-mean-square error between two lower-triangular distance vectors.

    Returns 0.0 for empty inputs.
    """
    n = len(flat_a)
    if n == 0:
        return 0.0
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(flat_a, flat_b)) / n)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _parse_profile(name):
    """Return (n, l) for a predefined or custom ``NxL`` profile.

    Returns ``None`` when the name cannot be interpreted.
    """
    if name in PROFILES:
        p = PROFILES[name]
        return p["n"], p["l"]
    parts = name.lower().split("x")
    if len(parts) == 2:
        try:
            return int(parts[0]), int(parts[1])
        except ValueError:
            pass
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark Bio.Phylo distance matrix backends"
    )
    parser.add_argument(
        "--profiles",
        nargs="+",
        default=list(PROFILES.keys()),
        help=(
            "Profiles to run: viral, amplicon, metagenomic, or NxL "
            "(default: all three baseline profiles)"
        ),
    )
    parser.add_argument(
        "--methods",
        nargs="+",
        default=METHODS,
        help="Backends to benchmark (default: all)",
    )
    parser.add_argument(
        "--model",
        default="identity",
        help="Scoring model — identity, blastn, blosum62, etc. (default: identity)",
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=3,
        help="Timing repetitions per method (default: 3)",
    )
    parser.add_argument(
        "--mutation-rate",
        type=float,
        default=0.1,
        help="Per-position mutation probability (default: 0.1)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for alignment generation (default: 42)",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Write results to this CSV file",
    )
    args = parser.parse_args()

    results = []

    for profile_name in args.profiles:
        parsed = _parse_profile(profile_name)
        if parsed is None:
            print(
                f"Unknown profile '{profile_name}'. Use a predefined name "
                f"({', '.join(PROFILES)}) or NxL format.",
                file=sys.stderr,
            )
            continue
        n, l = parsed

        print(f"\n{'=' * 62}")
        print(f"  Profile: {profile_name}  (N={n}, L={l})")
        print(f"{'=' * 62}")

        # Choose alphabet based on model
        alphabet = (
            DNA_ALPHABET if args.model in ("identity", "blastn") else PROTEIN_ALPHABET
        )
        msa = generate_alignment(
            n, l, alphabet=alphabet, mutation_rate=args.mutation_rate, seed=args.seed
        )

        # --------------------------------------------------------------------
        # Phase 1: Correctness validation (RMSE vs python legacy baseline)
        # --------------------------------------------------------------------
        # Cap the validation subset so that the python baseline is always fast.
        validate_n = min(n, VALIDATION_N)
        if validate_n < n:
            from Bio.Align import MultipleSeqAlignment

            val_msa = MultipleSeqAlignment(list(msa)[:validate_n])
        else:
            val_msa = msa

        print(f"\n  Correctness (RMSE vs python baseline, N={validate_n}):")

        baseline_flat = _get_flat_distances(val_msa, args.model, "python")

        # rmse_map collects the validated RMSE for each method so that Phase 2
        # can include it in the performance row without re-running validation.
        rmse_map = {}
        for method in args.methods:
            if method == "python":
                rmse_map["python"] = 0.0
                print(f"    [{'python':8s}] rmse=0.00e+00  (baseline)")
                continue
            try:
                val_flat = _get_flat_distances(val_msa, args.model, method)
                rmse = compute_rmse(val_flat, baseline_flat)
                rmse_map[method] = rmse
                status = "PASS" if rmse < 1e-10 else "WARN"
                print(f"    [{method:8s}] rmse={rmse:.2e}  {status}")
            except Exception as exc:
                rmse_map[method] = None
                print(f"    [{method:8s}] FAILED: {exc}")

        # --------------------------------------------------------------------
        # Phase 2: Performance benchmark (full alignment, timed repeats)
        # --------------------------------------------------------------------
        print(f"\n  Performance ({args.repeats} repeats, N={n}, L={l}):")

        from Bio.Phylo.TreeConstruction import DistanceCalculator

        for method in args.methods:
            rmse = rmse_map.get(method)

            # Method failed validation — skip performance entirely
            if rmse is None:
                continue

            # Skip python timing when N is impractical
            if method == "python" and n > PYTHON_TIMING_MAX_N:
                print(
                    f"    [{'python':8s}] timing skipped "
                    f"(N={n} > {PYTHON_TIMING_MAX_N})"
                )
                results.append(
                    {
                        "profile": profile_name,
                        "n": n,
                        "l": l,
                        "method": method,
                        "model": args.model,
                        "time_min": "N/A",
                        "time_mean": "N/A",
                        "mem_peak_mb": "N/A",
                        "rmse": f"{rmse:.2e}",
                    }
                )
                continue

            try:
                # Construct once — model loading excluded from timing
                dc = DistanceCalculator(args.model, method=method)

                times, mems = [], []
                for _ in range(args.repeats):
                    t, mem = run_timed(dc, msa)
                    times.append(t)
                    mems.append(mem)

                time_min = min(times)
                time_mean = sum(times) / len(times)
                mem_peak = max(mems)

                print(
                    f"    [{method:8s}] "
                    f"time_min={time_min:.3f}s  "
                    f"time_mean={time_mean:.3f}s  "
                    f"mem_peak={mem_peak:.1f} MB  "
                    f"rmse={rmse:.2e}"
                )
                results.append(
                    {
                        "profile": profile_name,
                        "n": n,
                        "l": l,
                        "method": method,
                        "model": args.model,
                        "time_min": f"{time_min:.4f}",
                        "time_mean": f"{time_mean:.4f}",
                        "mem_peak_mb": f"{mem_peak:.2f}",
                        "rmse": f"{rmse:.2e}",
                    }
                )
            except Exception as exc:
                print(f"    [{method:8s}] FAILED: {exc}")

    # Write CSV
    if args.output and results:
        with open(args.output, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        print(f"\nResults written to {args.output}")


if __name__ == "__main__":
    main()
