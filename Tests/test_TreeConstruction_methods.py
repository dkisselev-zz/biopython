#!/usr/bin/env python
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for distance matrix computation methods in Bio.Phylo.TreeConstruction.

This module tests equivalence between different computation strategies
(python, numpy, scipy, onehot) to ensure all methods produce numerically
identical or near-identical results.
"""

import unittest

from Bio import Align
from Bio.Phylo.TreeConstruction import DistanceCalculator


class TestDistanceComputationMethods(unittest.TestCase):
    """Test equivalence of different distance computation methods."""

    def setUp(self):
        """Load test alignment data."""
        # Use existing test alignment from Biopython test suite
        self.aln = Align.read("TreeConstruction/msa.phy", "phylip")

    def test_identity_model_equivalence(self):
        """Test that numpy method matches python for identity model."""
        calculator = DistanceCalculator("identity")

        # Compute with both methods
        dm_python = calculator.get_distance(self.aln, method="python")
        dm_numpy = calculator.get_distance(self.aln, method="numpy")

        # Check all pairwise distances match
        names = dm_python.names
        for i, name1 in enumerate(names):
            for j in range(i):
                name2 = names[j]
                dist_python = dm_python[name1, name2]
                dist_numpy = dm_numpy[name1, name2]
                self.assertAlmostEqual(
                    dist_python,
                    dist_numpy,
                    places=12,
                    msg=f"Mismatch for {name1}-{name2}: "
                    f"python={dist_python}, numpy={dist_numpy}",
                )

    def test_scipy_identity_model_equivalence(self):
        """Test that scipy method matches python for identity model."""
        # Check if SciPy is available
        try:
            import scipy  # noqa: F401
        except ImportError:
            from Bio import MissingPythonDependencyError

            raise MissingPythonDependencyError(
                "SciPy required for this test. Install with: pip install scipy"
            )

        calculator = DistanceCalculator("identity")

        # Compute with both methods
        dm_python = calculator.get_distance(self.aln, method="python")
        dm_scipy = calculator.get_distance(self.aln, method="scipy")

        # Check all pairwise distances match
        names = dm_python.names
        for i, name1 in enumerate(names):
            for j in range(i):
                name2 = names[j]
                dist_python = dm_python[name1, name2]
                dist_scipy = dm_scipy[name1, name2]
                self.assertAlmostEqual(
                    dist_python,
                    dist_scipy,
                    places=12,
                    msg=f"Mismatch for {name1}-{name2}: "
                    f"python={dist_python}, scipy={dist_scipy}",
                )

    def test_blosum62_model_equivalence(self):
        """Test that numpy method matches python for BLOSUM62 model."""
        calculator = DistanceCalculator("blosum62")

        # Compute with both methods
        dm_python = calculator.get_distance(self.aln, method="python")
        dm_numpy = calculator.get_distance(self.aln, method="numpy")

        # Check all pairwise distances match
        names = dm_python.names
        for i, name1 in enumerate(names):
            for j in range(i):
                name2 = names[j]
                dist_python = dm_python[name1, name2]
                dist_numpy = dm_numpy[name1, name2]
                self.assertAlmostEqual(
                    dist_python,
                    dist_numpy,
                    places=12,
                    msg=f"Mismatch for {name1}-{name2}: "
                    f"python={dist_python}, numpy={dist_numpy}",
                )

    def test_pam250_model_equivalence(self):
        """Test that numpy method matches python for PAM250 model."""
        calculator = DistanceCalculator("pam250")

        # Compute with both methods
        dm_python = calculator.get_distance(self.aln, method="python")
        dm_numpy = calculator.get_distance(self.aln, method="numpy")

        # Check all pairwise distances match
        names = dm_python.names
        for i, name1 in enumerate(names):
            for j in range(i):
                name2 = names[j]
                dist_python = dm_python[name1, name2]
                dist_numpy = dm_numpy[name1, name2]
                self.assertAlmostEqual(
                    dist_python,
                    dist_numpy,
                    places=12,
                    msg=f"Mismatch for {name1}-{name2}: "
                    f"python={dist_python}, numpy={dist_numpy}",
                )

    def test_onehot_identity_model_equivalence(self):
        """Test that onehot method matches python for identity model."""
        calculator = DistanceCalculator("identity")

        # Compute with both methods
        dm_python = calculator.get_distance(self.aln, method="python")
        dm_onehot = calculator.get_distance(self.aln, method="onehot")

        # Check all pairwise distances match (within tolerance for float32)
        names = dm_python.names
        for i, name1 in enumerate(names):
            for j in range(i):
                name2 = names[j]
                dist_python = dm_python[name1, name2]
                dist_onehot = dm_onehot[name1, name2]
                self.assertAlmostEqual(
                    dist_python,
                    dist_onehot,
                    places=6,  # Lower precision due to float32 in onehot
                    msg=f"Mismatch for {name1}-{name2}: "
                    f"python={dist_python}, onehot={dist_onehot}",
                )

    def test_onehot_blosum62_model_equivalence(self):
        """Test that onehot method matches python for BLOSUM62 model."""
        calculator = DistanceCalculator("blosum62")

        # Compute with both methods
        dm_python = calculator.get_distance(self.aln, method="python")
        dm_onehot = calculator.get_distance(self.aln, method="onehot")

        # Check all pairwise distances match (within tolerance for float32)
        names = dm_python.names
        for i, name1 in enumerate(names):
            for j in range(i):
                name2 = names[j]
                dist_python = dm_python[name1, name2]
                dist_onehot = dm_onehot[name1, name2]
                self.assertAlmostEqual(
                    dist_python,
                    dist_onehot,
                    places=6,  # Lower precision due to float32 in onehot
                    msg=f"Mismatch for {name1}-{name2}: "
                    f"python={dist_python}, onehot={dist_onehot}",
                )

    def test_matrix_symmetry(self):
        """Test that distance matrices are symmetric."""
        calculator = DistanceCalculator("blosum62")

        for method in ["python", "numpy"]:
            dm = calculator.get_distance(self.aln, method=method)
            names = dm.names

            # Check symmetry: dm[i,j] == dm[j,i]
            for i, name1 in enumerate(names):
                for j in range(i):
                    name2 = names[j]
                    dist_ij = dm[name1, name2]
                    dist_ji = dm[name2, name1]
                    self.assertAlmostEqual(
                        dist_ij,
                        dist_ji,
                        places=12,
                        msg=f"Asymmetry in {method} method: "
                        f"dm[{name1},{name2}]={dist_ij} != "
                        f"dm[{name2},{name1}]={dist_ji}",
                    )

    def test_distance_bounds(self):
        """Test that distances are in valid range [0, 1]."""
        calculator = DistanceCalculator("blosum62")

        for method in ["python", "numpy"]:
            dm = calculator.get_distance(self.aln, method=method)
            names = dm.names

            for i, name1 in enumerate(names):
                for j in range(i):
                    name2 = names[j]
                    dist = dm[name1, name2]
                    self.assertGreaterEqual(
                        dist,
                        0.0,
                        msg=f"{method} method produced negative distance: "
                        f"dm[{name1},{name2}]={dist}",
                    )
                    self.assertLessEqual(
                        dist,
                        1.0,
                        msg=f"{method} method produced distance > 1: "
                        f"dm[{name1},{name2}]={dist}",
                    )

    def test_empty_alignment(self):
        """Test handling of empty alignment."""
        from Bio.Phylo.TreeConstruction import DistanceMatrix

        # Create empty alignment
        empty_aln = Align.Alignment(sequences=[])

        calculator = DistanceCalculator("identity")

        for method in ["python", "numpy"]:
            dm = calculator.get_distance(empty_aln, method=method)
            self.assertIsInstance(dm, DistanceMatrix)
            self.assertEqual(len(dm.names), 0)

    def test_single_sequence(self):
        """Test handling of single sequence alignment."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        # Create single sequence alignment
        seq = SeqRecord(Seq("ACGT"), id="seq1")
        single_aln = Align.Alignment(sequences=[seq])

        calculator = DistanceCalculator("identity")

        for method in ["python", "numpy"]:
            dm = calculator.get_distance(single_aln, method=method)
            self.assertEqual(len(dm.names), 1)
            self.assertEqual(dm.names[0], "seq1")


class TestMethodPerformanceCharacteristics(unittest.TestCase):
    """Test basic performance characteristics of methods.

    These tests verify that methods complete successfully and produce
    reasonable results, without strict timing requirements.
    """

    def setUp(self):
        """Load test alignment data."""
        self.aln = Align.read("TreeConstruction/msa.phy", "phylip")

    def test_numpy_completes(self):
        """Test that numpy method completes successfully."""
        calculator = DistanceCalculator("blosum62")
        dm = calculator.get_distance(self.aln, method="numpy")
        self.assertIsNotNone(dm)
        self.assertEqual(len(dm.names), 5)

    def test_python_completes(self):
        """Test that python method completes successfully."""
        calculator = DistanceCalculator("blosum62")
        dm = calculator.get_distance(self.aln, method="python")
        self.assertIsNotNone(dm)
        self.assertEqual(len(dm.names), 5)

    def test_scipy_only_identity_model(self):
        """Test that scipy method raises error with substitution matrices."""
        try:
            import scipy  # noqa: F401
        except ImportError:
            from Bio import MissingPythonDependencyError

            raise MissingPythonDependencyError(
                "SciPy required for this test. Install with: pip install scipy"
            )

        calculator = DistanceCalculator("blosum62")

        # SciPy method should raise ValueError for substitution matrices
        with self.assertRaises(ValueError) as context:
            calculator.get_distance(self.aln, method="scipy")

        self.assertIn("identity model", str(context.exception).lower())

    def test_onehot_completes(self):
        """Test that onehot method completes successfully."""
        calculator = DistanceCalculator("blosum62")
        dm = calculator.get_distance(self.aln, method="onehot")
        self.assertIsNotNone(dm)
        self.assertEqual(len(dm.names), 5)


class TestParallelExecution(unittest.TestCase):
    """Test parallel execution of distance computation."""

    def setUp(self):
        """Load test alignment data."""
        self.aln = Align.read("TreeConstruction/msa.phy", "phylip")

    def test_parallel_produces_same_results(self):
        """Test that parallel execution produces identical results."""
        calculator = DistanceCalculator("blosum62")

        # Compute with serial execution
        dm_serial = calculator.get_distance(self.aln, method="python", n_jobs=1)

        # Compute with parallel execution
        dm_parallel = calculator.get_distance(self.aln, method="python", n_jobs=2)

        # Check all pairwise distances match
        names = dm_serial.names
        for i, name1 in enumerate(names):
            for j in range(i):
                name2 = names[j]
                dist_serial = dm_serial[name1, name2]
                dist_parallel = dm_parallel[name1, name2]
                self.assertAlmostEqual(
                    dist_serial,
                    dist_parallel,
                    places=12,
                    msg=f"Mismatch for {name1}-{name2}: "
                    f"serial={dist_serial}, parallel={dist_parallel}",
                )

    def test_parallel_all_cpus(self):
        """Test that n_jobs=-1 uses all CPUs."""
        calculator = DistanceCalculator("identity")

        # Should not raise an error
        dm = calculator.get_distance(self.aln, method="python", n_jobs=-1)
        self.assertIsNotNone(dm)
        self.assertEqual(len(dm.names), 5)

    def test_parallel_invalid_n_jobs(self):
        """Test that invalid n_jobs raises error."""
        calculator = DistanceCalculator("identity")

        # n_jobs = 0 should raise ValueError
        with self.assertRaises(ValueError):
            calculator.get_distance(self.aln, method="python", n_jobs=0)

        # n_jobs < -1 should raise ValueError
        with self.assertRaises(ValueError):
            calculator.get_distance(self.aln, method="python", n_jobs=-2)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
