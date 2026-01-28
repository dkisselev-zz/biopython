#!/usr/bin/env python
"""Quick benchmark to verify NumPy speedup for distance matrix computation.

This script generates a synthetic alignment and compares the performance
of the python and numpy methods.
"""

import time
from Bio import Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator


def generate_synthetic_alignment(n_seqs, seq_len, divergence=0.1, seed=42):
    """Generate a synthetic alignment with controlled divergence.

    Arguments:
     - n_seqs - Number of sequences
     - seq_len - Length of sequences
     - divergence - Fraction of positions to mutate (0.0-1.0)
     - seed - Random seed for reproducibility

    Returns an Alignment object.
    """
    import random

    random.seed(seed)

    # Generate base sequence
    alphabet = "ACDEFGHIKLMNPQRSTVWY"  # 20 amino acids
    base_seq = "".join(random.choices(alphabet, k=seq_len))

    # Generate diverged sequences
    sequences = []
    for i in range(n_seqs):
        # Mutate random positions
        seq_list = list(base_seq)
        n_mutations = int(seq_len * divergence * random.random())
        positions = random.sample(range(seq_len), n_mutations)
        for pos in positions:
            seq_list[pos] = random.choice(alphabet)

        seq_str = "".join(seq_list)
        record = SeqRecord(Seq(seq_str), id=f"seq{i}")
        sequences.append(record)

    return Align.Alignment(sequences=sequences)


def benchmark_method(calculator, alignment, method_name, n_runs=3):
    """Benchmark a distance computation method.

    Arguments:
     - calculator - DistanceCalculator instance
     - alignment - Alignment object
     - method_name - Method to benchmark ("python" or "numpy")
     - n_runs - Number of runs for timing

    Returns tuple (mean_time, distance_matrix).
    """
    times = []
    dm = None

    for _ in range(n_runs):
        start = time.time()
        dm = calculator.get_distance(alignment, method=method_name)
        end = time.time()
        times.append(end - start)

    mean_time = sum(times) / len(times)
    return mean_time, dm


def compare_results(dm1, dm2, tolerance=1e-10):
    """Compare two distance matrices for numerical equivalence.

    Arguments:
     - dm1 - First DistanceMatrix
     - dm2 - Second DistanceMatrix
     - tolerance - Maximum allowed difference

    Returns tuple (max_diff, all_close).
    """
    names = dm1.names
    max_diff = 0.0

    for i, name1 in enumerate(names):
        for j in range(i):
            name2 = names[j]
            diff = abs(dm1[name1, name2] - dm2[name1, name2])
            max_diff = max(max_diff, diff)

    all_close = max_diff < tolerance
    return max_diff, all_close


def main():
    """Run quick benchmark."""
    print("Quick Benchmark: Distance Computation Methods")
    print("=" * 60)
    print()

    # Check if SciPy is available
    try:
        import scipy  # noqa: F401

        scipy_available = True
    except ImportError:
        scipy_available = False
        print("Note: SciPy not available, skipping scipy benchmarks")
        print()

    # Test scenarios
    scenarios = [
        {"name": "Small (N=20, L=100)", "n_seqs": 20, "seq_len": 100},
        {"name": "Medium (N=50, L=200)", "n_seqs": 50, "seq_len": 200},
        {"name": "Large (N=100, L=500)", "n_seqs": 100, "seq_len": 500},
    ]

    # Test with BLOSUM62 (all methods)
    print("=" * 60)
    print("BLOSUM62 Model (python vs numpy vs onehot)")
    print("=" * 60)
    print()

    calculator = DistanceCalculator("blosum62")

    for scenario in scenarios:
        print(f"Scenario: {scenario['name']}")
        print("-" * 60)

        # Generate alignment
        aln = generate_synthetic_alignment(
            scenario["n_seqs"], scenario["seq_len"], divergence=0.2
        )

        # Benchmark Python method
        python_time, dm_python = benchmark_method(calculator, aln, "python", n_runs=3)

        # Benchmark NumPy method
        numpy_time, dm_numpy = benchmark_method(calculator, aln, "numpy", n_runs=3)

        # Benchmark OneHot method
        onehot_time, dm_onehot = benchmark_method(calculator, aln, "onehot", n_runs=3)

        # Compare results
        max_diff_numpy, all_close_numpy = compare_results(
            dm_python, dm_numpy, tolerance=1e-10
        )
        max_diff_onehot, all_close_onehot = compare_results(
            dm_python, dm_onehot, tolerance=1e-5
        )  # Lower tolerance for float32

        # Calculate speedups
        speedup_numpy = python_time / numpy_time if numpy_time > 0 else 0
        speedup_onehot = python_time / onehot_time if onehot_time > 0 else 0

        # Display results
        print(f"  Python method: {python_time:.3f} seconds")
        print(
            f"  NumPy method:  {numpy_time:.3f} seconds (speedup: {speedup_numpy:.2f}x)"
        )
        print(
            f"  OneHot method: {onehot_time:.3f} seconds (speedup: {speedup_onehot:.2f}x)"
        )
        print(
            f"  NumPy max diff: {max_diff_numpy:.2e} ({'✓' if all_close_numpy else '✗'})"
        )
        print(
            f"  OneHot max diff: {max_diff_onehot:.2e} ({'✓' if all_close_onehot else '✗'})"
        )
        print()

    # Test with identity model (all methods)
    if scipy_available:
        print("=" * 60)
        print("Identity Model (python vs numpy vs scipy vs onehot)")
        print("=" * 60)
        print()

        calculator = DistanceCalculator("identity")

        for scenario in scenarios:
            print(f"Scenario: {scenario['name']}")
            print("-" * 60)

            # Generate alignment
            aln = generate_synthetic_alignment(
                scenario["n_seqs"], scenario["seq_len"], divergence=0.2
            )

            # Benchmark Python method
            python_time, dm_python = benchmark_method(
                calculator, aln, "python", n_runs=3
            )

            # Benchmark NumPy method
            numpy_time, dm_numpy = benchmark_method(calculator, aln, "numpy", n_runs=3)

            # Benchmark SciPy method
            scipy_time, dm_scipy = benchmark_method(calculator, aln, "scipy", n_runs=3)

            # Benchmark OneHot method
            onehot_time, dm_onehot = benchmark_method(
                calculator, aln, "onehot", n_runs=3
            )

            # Compare results
            max_diff_numpy, all_close_numpy = compare_results(
                dm_python, dm_numpy, tolerance=1e-10
            )
            max_diff_scipy, all_close_scipy = compare_results(
                dm_python, dm_scipy, tolerance=1e-10
            )
            max_diff_onehot, all_close_onehot = compare_results(
                dm_python, dm_onehot, tolerance=1e-5
            )  # Lower tolerance for float32

            # Calculate speedups
            speedup_numpy = python_time / numpy_time if numpy_time > 0 else 0
            speedup_scipy = python_time / scipy_time if scipy_time > 0 else 0
            speedup_onehot = python_time / onehot_time if onehot_time > 0 else 0

            # Display results
            print(f"  Python method: {python_time:.3f} seconds")
            print(
                f"  NumPy method:  {numpy_time:.3f} seconds (speedup: {speedup_numpy:.2f}x)"
            )
            print(
                f"  SciPy method:  {scipy_time:.3f} seconds (speedup: {speedup_scipy:.2f}x)"
            )
            print(
                f"  OneHot method: {onehot_time:.3f} seconds (speedup: {speedup_onehot:.2f}x)"
            )
            print(
                f"  NumPy max diff: {max_diff_numpy:.2e} ({'✓' if all_close_numpy else '✗'})"
            )
            print(
                f"  SciPy max diff: {max_diff_scipy:.2e} ({'✓' if all_close_scipy else '✗'})"
            )
            print(
                f"  OneHot max diff: {max_diff_onehot:.2e} ({'✓' if all_close_onehot else '✗'})"
            )
            print()


if __name__ == "__main__":
    main()
