# Copyright 2024 by Biopython contributors. All rights reserved.
#
# This file is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.Phylo.DistanceMatrix module."""

import os
import unittest

import numpy as np

from Bio import Align
from Bio import AlignIO
from Bio.Phylo.DistanceMatrix import FastDistanceCalculator
from Bio.Phylo.DistanceMatrix import list_strategies
from Bio.Phylo.DistanceMatrix._base import condensed_to_squareform
from Bio.Phylo.DistanceMatrix._base import create_gap_mask
from Bio.Phylo.DistanceMatrix._parallel import get_n_jobs
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor


class TestParallelUtilities(unittest.TestCase):
    """Tests for parallelism utilities."""

    def test_get_n_jobs_default(self):
        """Test default n_jobs is 1."""
        self.assertEqual(get_n_jobs(), 1)

    def test_get_n_jobs_explicit(self):
        """Test explicit n_jobs argument."""
        self.assertEqual(get_n_jobs(4), 4)
        self.assertEqual(get_n_jobs(1), 1)

    def test_get_n_jobs_negative(self):
        """Test n_jobs=-1 returns CPU count."""
        n = get_n_jobs(-1)
        self.assertGreater(n, 0)

    def test_get_n_jobs_zero_raises(self):
        """Test n_jobs=0 raises ValueError."""
        self.assertRaises(ValueError, get_n_jobs, 0)

    def test_get_n_jobs_env_variable(self):
        """Test BIOPYTHON_NUM_THREADS environment variable."""
        old_env = os.environ.get("BIOPYTHON_NUM_THREADS")
        try:
            os.environ["BIOPYTHON_NUM_THREADS"] = "3"
            # Explicit argument takes precedence
            self.assertEqual(get_n_jobs(2), 2)
            # Env var used when no argument
            self.assertEqual(get_n_jobs(), 3)
        finally:
            if old_env is None:
                os.environ.pop("BIOPYTHON_NUM_THREADS", None)
            else:
                os.environ["BIOPYTHON_NUM_THREADS"] = old_env


class TestBaseUtilities(unittest.TestCase):
    """Tests for base utilities."""

    def test_create_gap_mask(self):
        """Test gap mask creation."""
        alignment = np.array([[b"A", b"C", b"-", b"G"], [b"-", b"T", b"A", b"-"]])
        mask = create_gap_mask(alignment, ["-"])
        expected = np.array([[False, False, True, False], [True, False, False, True]])
        np.testing.assert_array_equal(mask, expected)

    def test_condensed_to_squareform(self):
        """Test condensed to square matrix conversion."""
        condensed = np.array([0.1, 0.2, 0.3])  # 3 elements for n=3
        square = condensed_to_squareform(condensed, 3)
        expected = np.array([[0.0, 0.1, 0.2], [0.1, 0.0, 0.3], [0.2, 0.3, 0.0]])
        np.testing.assert_array_almost_equal(square, expected)

    def test_list_strategies(self):
        """Test listing available strategies."""
        strategies = list_strategies()
        self.assertIn("python", strategies)
        self.assertIn("numpy", strategies)
        self.assertIn("scipy", strategies)
        self.assertIn("onehot", strategies)


class DistanceStrategyTestMixin:
    """Mixin providing common tests for all strategies."""

    strategy_name = None

    def setUp(self):
        """Set up test data."""
        # Load test alignment
        self.alignment = Align.read("TreeConstruction/msa.phy", "phylip")
        self.msa = AlignIO.read("TreeConstruction/msa.phy", "phylip")

        # Known expected values from legacy implementation
        # Identity model - values from DistanceCalculator("identity")
        self.expected_identity = {
            ("Alpha", "Beta"): 0.23076923076923073,  # 1 - 10/13
            ("Alpha", "Gamma"): 0.3846153846153846,  # 1 - 8/13
            ("Alpha", "Delta"): 0.5384615384615384,  # 1 - 6/13
            ("Alpha", "Epsilon"): 0.6153846153846154,  # 1 - 5/13
        }
        # BLOSUM62 model
        self.expected_blosum62 = {
            ("Alpha", "Beta"): 1 - 53 / 84,
        }

    def test_identity_model(self):
        """Test identity distance matches known values."""
        if self.strategy_name is None:
            self.skipTest("Base mixin class")

        calc = FastDistanceCalculator("identity", method=self.strategy_name)
        dm = calc.get_distance(self.alignment)

        for (name1, name2), expected in self.expected_identity.items():
            with self.subTest(pair=(name1, name2)):
                self.assertAlmostEqual(dm[name1, name2], expected, places=10)

    def test_identity_model_msa(self):
        """Test identity distance with legacy MultipleSeqAlignment."""
        if self.strategy_name is None:
            self.skipTest("Base mixin class")

        calc = FastDistanceCalculator("identity", method=self.strategy_name)
        dm = calc.get_distance(self.msa)

        for (name1, name2), expected in self.expected_identity.items():
            with self.subTest(pair=(name1, name2)):
                self.assertAlmostEqual(dm[name1, name2], expected, places=10)

    def test_substitution_model(self):
        """Test substitution matrix distance matches known values."""
        if self.strategy_name is None:
            self.skipTest("Base mixin class")

        calc = FastDistanceCalculator("blosum62", method=self.strategy_name)
        dm = calc.get_distance(self.alignment)

        for (name1, name2), expected in self.expected_blosum62.items():
            with self.subTest(pair=(name1, name2)):
                self.assertAlmostEqual(dm[name1, name2], expected, places=10)

    def test_returns_distance_matrix(self):
        """Test that result is a DistanceMatrix."""
        if self.strategy_name is None:
            self.skipTest("Base mixin class")

        calc = FastDistanceCalculator("identity", method=self.strategy_name)
        dm = calc.get_distance(self.alignment)
        self.assertIsInstance(dm, DistanceMatrix)

    def test_diagonal_is_zero(self):
        """Test that diagonal elements are zero."""
        if self.strategy_name is None:
            self.skipTest("Base mixin class")

        calc = FastDistanceCalculator("identity", method=self.strategy_name)
        dm = calc.get_distance(self.alignment)
        for name in dm.names:
            self.assertEqual(dm[name, name], 0.0)

    def test_symmetry(self):
        """Test that distance matrix is symmetric."""
        if self.strategy_name is None:
            self.skipTest("Base mixin class")

        calc = FastDistanceCalculator("identity", method=self.strategy_name)
        dm = calc.get_distance(self.alignment)
        for i, name1 in enumerate(dm.names):
            for j, name2 in enumerate(dm.names):
                self.assertAlmostEqual(dm[name1, name2], dm[name2, name1], places=10)


class TestPythonStrategy(DistanceStrategyTestMixin, unittest.TestCase):
    """Tests for Python strategy."""

    strategy_name = "python"


class TestNumPyStrategy(DistanceStrategyTestMixin, unittest.TestCase):
    """Tests for NumPy strategy."""

    strategy_name = "numpy"


class TestSciPyStrategy(DistanceStrategyTestMixin, unittest.TestCase):
    """Tests for SciPy strategy."""

    strategy_name = "scipy"

    def setUp(self):
        """Check scipy is available."""
        try:
            import scipy  # noqa: F401
        except ImportError:
            self.skipTest("SciPy not available")
        super().setUp()


class TestOneHotStrategy(DistanceStrategyTestMixin, unittest.TestCase):
    """Tests for OneHot strategy."""

    strategy_name = "onehot"


class TestLegacyCompatibility(unittest.TestCase):
    """Tests ensuring new implementation matches legacy DistanceCalculator."""

    def setUp(self):
        """Set up test data."""
        self.alignment = Align.read("TreeConstruction/msa.phy", "phylip")

    def test_legacy_vs_python_identity(self):
        """Test Python strategy matches legacy for identity model."""
        legacy = DistanceCalculator("identity")
        dm_legacy = legacy.get_distance(self.alignment)

        new = FastDistanceCalculator("identity", method="python")
        dm_new = new.get_distance(self.alignment)

        for name1 in dm_legacy.names:
            for name2 in dm_legacy.names:
                self.assertAlmostEqual(
                    dm_legacy[name1, name2], dm_new[name1, name2], places=10
                )

    def test_legacy_vs_numpy_identity(self):
        """Test NumPy strategy matches legacy for identity model."""
        legacy = DistanceCalculator("identity")
        dm_legacy = legacy.get_distance(self.alignment)

        new = FastDistanceCalculator("identity", method="numpy")
        dm_new = new.get_distance(self.alignment)

        for name1 in dm_legacy.names:
            for name2 in dm_legacy.names:
                self.assertAlmostEqual(
                    dm_legacy[name1, name2], dm_new[name1, name2], places=10
                )

    def test_legacy_vs_python_blosum62(self):
        """Test Python strategy matches legacy for blosum62 model."""
        legacy = DistanceCalculator("blosum62")
        dm_legacy = legacy.get_distance(self.alignment)

        new = FastDistanceCalculator("blosum62", method="python")
        dm_new = new.get_distance(self.alignment)

        for name1 in dm_legacy.names:
            for name2 in dm_legacy.names:
                self.assertAlmostEqual(
                    dm_legacy[name1, name2], dm_new[name1, name2], places=10
                )

    def test_legacy_vs_numpy_blosum62(self):
        """Test NumPy strategy matches legacy for blosum62 model."""
        legacy = DistanceCalculator("blosum62")
        dm_legacy = legacy.get_distance(self.alignment)

        new = FastDistanceCalculator("blosum62", method="numpy")
        dm_new = new.get_distance(self.alignment)

        for name1 in dm_legacy.names:
            for name2 in dm_legacy.names:
                self.assertAlmostEqual(
                    dm_legacy[name1, name2], dm_new[name1, name2], places=10
                )

    def test_all_strategies_match(self):
        """Test all strategies produce same results."""
        strategies = ["python", "numpy", "onehot"]
        try:
            import scipy  # noqa: F401

            strategies.append("scipy")
        except ImportError:
            pass

        results = {}
        for strategy in strategies:
            calc = FastDistanceCalculator("identity", method=strategy)
            results[strategy] = calc.get_distance(self.alignment)

        # Compare all pairs
        baseline = results["python"]
        for strategy, dm in results.items():
            if strategy == "python":
                continue
            with self.subTest(strategy=strategy):
                for name1 in baseline.names:
                    for name2 in baseline.names:
                        self.assertAlmostEqual(
                            baseline[name1, name2],
                            dm[name1, name2],
                            places=10,
                            msg=f"Mismatch for {name1}, {name2}",
                        )


class TestIntegration(unittest.TestCase):
    """Integration tests with tree construction."""

    def setUp(self):
        """Set up test data."""
        self.alignment = Align.read("TreeConstruction/msa.phy", "phylip")

    def test_upgma_tree_construction(self):
        """Test that computed distance matrix works with UPGMA."""
        calc = FastDistanceCalculator("blosum62", method="numpy")
        dm = calc.get_distance(self.alignment)

        constructor = DistanceTreeConstructor()
        tree = constructor.upgma(dm)

        from Bio.Phylo import BaseTree

        self.assertIsInstance(tree, BaseTree.Tree)

    def test_nj_tree_construction(self):
        """Test that computed distance matrix works with Neighbor-Joining."""
        calc = FastDistanceCalculator("blosum62", method="numpy")
        dm = calc.get_distance(self.alignment)

        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)

        from Bio.Phylo import BaseTree

        self.assertIsInstance(tree, BaseTree.Tree)


class TestEdgeCases(unittest.TestCase):
    """Tests for edge cases."""

    def test_two_sequences(self):
        """Test with minimal two-sequence alignment."""
        from io import StringIO

        aln = Align.read(StringIO(">A\nACGT\n>B\nACGT"), "fasta")
        calc = FastDistanceCalculator("identity", method="numpy")
        dm = calc.get_distance(aln)

        self.assertEqual(dm["A", "B"], 0.0)  # Identical sequences

    def test_with_all_gaps(self):
        """Test handling of sequences that are all gaps at some positions."""
        from io import StringIO

        aln = Align.read(StringIO(">A\nA-G\n>B\n-C-"), "fasta")
        calc = FastDistanceCalculator("identity", method="numpy")
        dm = calc.get_distance(aln)

        # Should handle gracefully (no valid positions to compare)
        self.assertTrue(0.0 <= dm["A", "B"] <= 1.0)

    def test_invalid_model_raises(self):
        """Test that invalid model name raises ValueError."""
        self.assertRaises(
            ValueError, FastDistanceCalculator, "invalid_model", method="numpy"
        )

    def test_invalid_strategy_raises(self):
        """Test that invalid strategy name raises ValueError."""
        self.assertRaises(
            ValueError, FastDistanceCalculator, "identity", method="invalid"
        )


class TestNumericalAccuracy(unittest.TestCase):
    """Tests for numerical accuracy across strategies."""

    def setUp(self):
        """Generate larger test alignment for numerical accuracy testing."""
        np.random.seed(42)
        n_seqs = 20
        length = 100
        alphabet = b"ACGT"

        indices = np.random.randint(0, len(alphabet), size=(n_seqs, length))
        self.alignment = np.array(
            [[alphabet[i : i + 1] for i in row] for row in indices], dtype="S1"
        )

        # Add some gaps
        gap_mask = np.random.random((n_seqs, length)) < 0.05
        self.alignment[gap_mask] = b"-"

        self.names = [f"seq_{i}" for i in range(n_seqs)]

    def test_rmse_across_strategies(self):
        """Test RMSE between strategies is within tolerance."""
        from Bio.Phylo.DistanceMatrix._base import get_strategy

        strategies = ["python", "numpy", "onehot"]
        try:
            import scipy  # noqa: F401

            strategies.append("scipy")
        except ImportError:
            pass

        results = {}
        for name in strategies:
            strategy = get_strategy(name)
            results[name] = strategy.compute(
                self.alignment,
                scoring_matrix=None,
                alphabet=None,
                skip_letters=["-"],
                n_jobs=1,
            )

        # Compare all against python baseline
        baseline = results["python"]
        for name, distances in results.items():
            if name == "python":
                continue
            rmse = np.sqrt(np.mean((distances - baseline) ** 2))
            max_abs_error = np.max(np.abs(distances - baseline))
            with self.subTest(strategy=name):
                self.assertLess(rmse, 1e-10, f"RMSE too high for {name}")
                self.assertLess(
                    max_abs_error, 1e-9, f"Max abs error too high for {name}"
                )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
