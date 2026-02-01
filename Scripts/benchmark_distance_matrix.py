#!/usr/bin/env python
# Copyright 2024 by Biopython contributors. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
r"""Benchmark distance matrix computation strategies.

This script measures wall-clock time, memory usage, and correctness for
different distance calculation strategies.

Usage
-----
Run single profile::

    python benchmark_distance_matrix.py --profile 16s

Run all profiles with all strategies::

    python benchmark_distance_matrix.py --all-profiles

Custom configuration::

    python benchmark_distance_matrix.py \
        --n-seqs 500 --length 2000 \
        --strategies numpy scipy \
        --n-jobs 1 2 4 \
        --output results.json

Quick correctness check::

    python benchmark_distance_matrix.py --profile viral --correctness-only

"""

import argparse
import json
import sys
import time
import tracemalloc
from typing import Any

import numpy as np

# Benchmark profiles
PROFILES = {
    "viral": {"n_seqs": 100, "length": 30000, "description": "Viral genomes"},
    "16s": {"n_seqs": 1000, "length": 1500, "description": "16S rRNA amplicons"},
    "metagenomic": {
        "n_seqs": 10000,
        "length": 300,
        "description": "Metagenomic reads",
    },
    "small": {"n_seqs": 50, "length": 500, "description": "Small test (quick)"},
}

STRATEGIES = ["python", "numpy", "scipy", "onehot"]
MODELS = ["identity", "blosum62"]

# Correctness thresholds
RMSE_THRESHOLD = 1e-10
MAX_ABS_ERROR_THRESHOLD = 1e-9


def generate_synthetic_alignment(
    n_seqs: int, length: int, gap_rate: float = 0.05, seed: int = 42
) -> tuple[np.ndarray, list[str]]:
    """Generate reproducible synthetic alignment for benchmarking.

    :param n_seqs: Number of sequences.
    :param length: Alignment length.
    :param gap_rate: Fraction of positions that are gaps.
    :param seed: Random seed for reproducibility.
    :returns: Tuple of (alignment array, sequence names).
    """
    rng = np.random.default_rng(seed)
    alphabet = b"ACGT"

    # Generate random sequences
    indices = rng.integers(0, len(alphabet), size=(n_seqs, length))
    alignment = np.array(
        [[alphabet[i : i + 1] for i in row] for row in indices], dtype="S1"
    )

    # Introduce gaps
    gap_mask = rng.random((n_seqs, length)) < gap_rate
    alignment[gap_mask] = b"-"

    # Generate unique sequence names
    names = [f"seq_{i:05d}" for i in range(n_seqs)]

    return alignment, names


def compute_rmse(distances_new: np.ndarray, distances_legacy: np.ndarray) -> float:
    """Compute Root Mean Square Error between two distance arrays."""
    diff = distances_new - distances_legacy
    return float(np.sqrt(np.mean(diff**2)))


def compute_max_abs_error(
    distances_new: np.ndarray, distances_legacy: np.ndarray
) -> float:
    """Compute maximum absolute error between two distance arrays."""
    return float(np.max(np.abs(distances_new - distances_legacy)))


def benchmark_strategy(
    strategy_name: str,
    alignment: np.ndarray,
    model: str = "identity",
    n_jobs: int = 1,
    n_runs: int = 3,
    baseline_distances: np.ndarray | None = None,
) -> dict[str, Any]:
    """Benchmark a single strategy.

    :param strategy_name: Name of the strategy to benchmark.
    :param alignment: Alignment array.
    :param model: Distance model name.
    :param n_jobs: Number of parallel workers.
    :param n_runs: Number of timing runs.
    :param baseline_distances: Optional baseline for correctness comparison.
    :returns: Dictionary with benchmark results.
    """
    from Bio.Phylo.DistanceMatrix._base import get_strategy
    from Bio.Align import substitution_matrices

    results: dict[str, Any] = {
        "strategy": strategy_name,
        "model": model,
        "n_jobs": n_jobs,
        "n_seqs": alignment.shape[0],
        "length": alignment.shape[1],
    }

    try:
        strategy = get_strategy(strategy_name)

        # Set up scoring matrix
        if model == "identity":
            scoring_matrix = None
            alphabet = None
            skip_letters = ()
        else:
            if model == "blastn":
                name = "NUC.4.4"
            else:
                name = model.upper()
            matrix = substitution_matrices.load(name)
            scoring_matrix = np.array(matrix, dtype=np.float64)
            alphabet = matrix.alphabet
            skip_letters = ("-", "*")

        # Warmup run (smaller subset)
        warmup_size = min(10, alignment.shape[0])
        _ = strategy.compute(
            alignment[:warmup_size],
            scoring_matrix=scoring_matrix,
            alphabet=alphabet,
            skip_letters=skip_letters,
            n_jobs=1,
        )

        # Timing runs
        times = []
        distances = None
        peak_memory = 0

        for run in range(n_runs):
            tracemalloc.start()
            start = time.perf_counter()

            distances = strategy.compute(
                alignment,
                scoring_matrix=scoring_matrix,
                alphabet=alphabet,
                skip_letters=skip_letters,
                n_jobs=n_jobs,
            )

            elapsed = time.perf_counter() - start
            current, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()

            times.append(elapsed)
            peak_memory = max(peak_memory, peak)

        results["time_mean_sec"] = float(np.mean(times))
        results["time_std_sec"] = float(np.std(times))
        results["time_min_sec"] = float(np.min(times))
        results["memory_peak_mb"] = peak_memory / 1024 / 1024
        results["success"] = True

        # Correctness check
        if baseline_distances is not None and distances is not None:
            rmse = compute_rmse(distances, baseline_distances)
            max_abs_error = compute_max_abs_error(distances, baseline_distances)
            results["rmse_vs_legacy"] = rmse
            results["max_abs_error"] = max_abs_error
            results["correctness_passed"] = (
                rmse < RMSE_THRESHOLD and max_abs_error < MAX_ABS_ERROR_THRESHOLD
            )

    except ImportError as e:
        results["success"] = False
        results["error"] = f"ImportError: {e}"
    except Exception as e:
        results["success"] = False
        results["error"] = str(e)

    return results


def run_benchmarks(
    profiles: list[tuple[str, dict]],
    strategies: list[str],
    models: list[str],
    n_jobs_list: list[int],
    n_runs: int,
    correctness_only: bool = False,
) -> list[dict[str, Any]]:
    """Run benchmarks for all configurations.

    :param profiles: List of (profile_name, profile_config) tuples.
    :param strategies: List of strategy names to benchmark.
    :param models: List of model names to benchmark.
    :param n_jobs_list: List of n_jobs values to test.
    :param n_runs: Number of timing runs.
    :param correctness_only: If True, only check correctness (skip timing).
    :returns: List of benchmark results.
    """
    results = []

    for profile_name, profile in profiles:
        print(f"\n{'=' * 60}")
        print(
            f"Profile: {profile_name} "
            f"(n_seqs={profile['n_seqs']}, length={profile['length']})"
        )
        print(f"{'=' * 60}")

        alignment, names = generate_synthetic_alignment(
            profile["n_seqs"], profile["length"]
        )

        for model in models:
            # Compute baseline (python strategy) for correctness comparison
            print(f"\n  Computing baseline ({model})...")
            baseline_result = benchmark_strategy(
                "python",
                alignment,
                model=model,
                n_jobs=1,
                n_runs=1 if correctness_only else n_runs,
            )

            if not baseline_result["success"]:
                print(f"    Baseline failed: {baseline_result.get('error', 'Unknown')}")
                continue

            # Get baseline distances for comparison
            from Bio.Phylo.DistanceMatrix._base import get_strategy
            from Bio.Align import substitution_matrices

            python_strategy = get_strategy("python")
            if model == "identity":
                scoring_matrix = None
                alphabet = None
                skip_letters = ()
            else:
                if model == "blastn":
                    name = "NUC.4.4"
                else:
                    name = model.upper()
                matrix = substitution_matrices.load(name)
                scoring_matrix = np.array(matrix, dtype=np.float64)
                alphabet = matrix.alphabet
                skip_letters = ("-", "*")

            baseline_distances = python_strategy.compute(
                alignment,
                scoring_matrix=scoring_matrix,
                alphabet=alphabet,
                skip_letters=skip_letters,
                n_jobs=1,
            )

            baseline_result["profile"] = profile_name
            results.append(baseline_result)
            baseline_time = baseline_result["time_mean_sec"]

            print(
                f"    python | {model:10s} | n_jobs=1  | "
                f"time={baseline_time:.3f}s | "
                f"mem={baseline_result['memory_peak_mb']:.1f}MB"
            )

            for strategy in strategies:
                if strategy == "python":
                    continue  # Already computed as baseline

                for n_jobs in n_jobs_list:
                    result = benchmark_strategy(
                        strategy,
                        alignment,
                        model=model,
                        n_jobs=n_jobs,
                        n_runs=1 if correctness_only else n_runs,
                        baseline_distances=baseline_distances,
                    )
                    result["profile"] = profile_name

                    if result["success"]:
                        speedup = baseline_time / result["time_mean_sec"]
                        result["speedup_vs_python"] = speedup

                        correctness_str = ""
                        if "correctness_passed" in result:
                            status = "PASS" if result["correctness_passed"] else "FAIL"
                            correctness_str = (
                                f" | rmse={result['rmse_vs_legacy']:.2e} [{status}]"
                            )

                        print(
                            f"    {strategy:7s} | {model:10s} | n_jobs={n_jobs:<2d} | "
                            f"time={result['time_mean_sec']:.3f}s | "
                            f"mem={result['memory_peak_mb']:.1f}MB | "
                            f"speedup={speedup:.1f}x{correctness_str}"
                        )
                    else:
                        print(
                            f"    {strategy:7s} | {model:10s} | n_jobs={n_jobs:<2d} | "
                            f"FAILED: {result.get('error', 'Unknown')}"
                        )

                    results.append(result)

    return results


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Benchmark distance matrix computation strategies"
    )
    parser.add_argument(
        "--profile",
        choices=list(PROFILES.keys()),
        help="Use predefined profile",
    )
    parser.add_argument(
        "--all-profiles",
        action="store_true",
        help="Run all profiles",
    )
    parser.add_argument(
        "--n-seqs",
        type=int,
        default=100,
        help="Number of sequences (default: 100)",
    )
    parser.add_argument(
        "--length",
        type=int,
        default=1000,
        help="Alignment length (default: 1000)",
    )
    parser.add_argument(
        "--strategies",
        nargs="+",
        default=STRATEGIES,
        help=f"Strategies to test (default: {STRATEGIES})",
    )
    parser.add_argument(
        "--models",
        nargs="+",
        default=["identity"],
        help="Distance models (default: identity)",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        nargs="+",
        default=[1],
        help="Parallelism levels (default: 1)",
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Output JSON file",
    )
    parser.add_argument(
        "--n-runs",
        type=int,
        default=3,
        help="Number of timing runs (default: 3)",
    )
    parser.add_argument(
        "--correctness-only",
        action="store_true",
        help="Only check correctness, skip timing",
    )

    args = parser.parse_args()

    # Determine profiles to run
    if args.all_profiles:
        profiles = list(PROFILES.items())
    elif args.profile:
        profiles = [(args.profile, PROFILES[args.profile])]
    else:
        profiles = [
            (
                "custom",
                {"n_seqs": args.n_seqs, "length": args.length, "description": "Custom"},
            )
        ]

    print("Distance Matrix Benchmark")
    print("=" * 60)
    print(f"Strategies: {args.strategies}")
    print(f"Models: {args.models}")
    print(f"n_jobs: {args.n_jobs}")
    print(f"n_runs: {args.n_runs}")

    results = run_benchmarks(
        profiles=profiles,
        strategies=args.strategies,
        models=args.models,
        n_jobs_list=args.n_jobs,
        n_runs=args.n_runs,
        correctness_only=args.correctness_only,
    )

    # Summary
    print(f"\n{'=' * 60}")
    print("Summary")
    print(f"{'=' * 60}")

    successful = [r for r in results if r.get("success", False)]
    failed = [r for r in results if not r.get("success", False)]
    correctness_failed = [r for r in successful if r.get("correctness_passed") is False]

    print(f"Total runs: {len(results)}")
    print(f"Successful: {len(successful)}")
    print(f"Failed: {len(failed)}")
    print(f"Correctness failures: {len(correctness_failed)}")

    if correctness_failed:
        print("\nCorrectness failures:")
        for r in correctness_failed:
            print(
                f"  - {r['strategy']} ({r['model']}): "
                f"rmse={r.get('rmse_vs_legacy', 'N/A')}, "
                f"max_abs={r.get('max_abs_error', 'N/A')}"
            )

    if args.output:
        with open(args.output, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nResults written to {args.output}")

    # Return non-zero if any correctness failures
    return 1 if correctness_failed else 0


if __name__ == "__main__":
    sys.exit(main())
