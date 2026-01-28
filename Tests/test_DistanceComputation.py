# Copyright 2024 The Biopython Contributors.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Unit tests for Bio.Phylo.DistanceComputation strategies.

Validates that every accelerated method produces results numerically
identical to the legacy pure-Python loop, integrates correctly with
NJ and UPGMA tree constructors, and handles edge cases (gaps, custom
skip_letters, subset alignment) consistently.
"""

import os
import unittest
import warnings
from io import StringIO

from Bio import Align
from Bio import AlignIO
from Bio import BiopythonExperimentalWarning
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

try:
    import scipy  # noqa: F401

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

# Path relative to Tests/ working directory
MSA_PHY = "TreeConstruction/msa.phy"


class BackwardCompatibilityTest(unittest.TestCase):
    """Ensure the default API path is unchanged."""

    def test_default_method_is_python(self):
        """Default method leaves _strategy as None."""
        calc = DistanceCalculator("identity")
        self.assertIsNone(calc._strategy)

    @unittest.skipUnless(HAS_NUMPY, "numpy required")
    def test_isinstance_gate_passes_for_all_methods(self):
        """DistanceTreeConstructor accepts calculators with any method."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonExperimentalWarning)
            calc = DistanceCalculator("identity", method="numpy")
        self.assertIsInstance(calc, DistanceCalculator)
        # Must not raise
        DistanceTreeConstructor(distance_calculator=calc)

    def test_invalid_method_raises_valueerror(self):
        """Unrecognised method name raises ValueError."""
        with self.assertRaises(ValueError):
            DistanceCalculator("identity", method="nonexistent")

    def test_legacy_known_values_unchanged(self):
        """Legacy path still produces exact reference values."""
        msa = AlignIO.read(MSA_PHY, "phylip")
        dm = DistanceCalculator("identity", method="python").get_distance(msa)
        self.assertAlmostEqual(dm["Alpha", "Beta"], 1 - 10 / 13)

    @unittest.skipUnless(HAS_NUMPY, "numpy required")
    def test_experimental_warning_emitted(self):
        """Non-python methods emit BiopythonExperimentalWarning."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            DistanceCalculator("identity", method="numpy")
        self.assertTrue(
            any(issubclass(x.category, BiopythonExperimentalWarning) for x in w)
        )


def _build_method_test_class(method_name):
    """Factory: build a test class for a given strategy method.

    This generates NumpyStrategyTest, ScipyStrategyTest, and
    OneHotStrategyTest with identical logic, differing only in
    the method string passed to DistanceCalculator.
    """

    class _MethodTest(unittest.TestCase):
        """Tests for an accelerated distance computation method."""

        def setUp(self):
            """Load test alignment in both object types."""
            self.msa = AlignIO.read(MSA_PHY, "phylip")
            self.aln = Align.read(MSA_PHY, "phylip")

        def _calc(self, model="identity"):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", BiopythonExperimentalWarning)
                return DistanceCalculator(model, method=method_name)

        # --- Known-value assertions ---

        def test_identity_known_value(self):
            """Identity distance for Alpha/Beta matches hand-calculated value."""
            dm = self._calc("identity").get_distance(self.msa)
            self.assertAlmostEqual(dm["Alpha", "Beta"], 1 - 10 / 13)

        def test_blosum62_known_value(self):
            """BLOSUM62 distance for Alpha/Beta matches hand-calculated value."""
            dm = self._calc("blosum62").get_distance(self.msa)
            self.assertAlmostEqual(dm["Alpha", "Beta"], 1 - 53 / 84)

        def test_blastn_known_value(self):
            """BLASTN distance for Alpha/Beta matches hand-calculated value."""
            dm = self._calc("blastn").get_distance(self.msa)
            self.assertAlmostEqual(dm["Alpha", "Beta"], 1 - 38 / 65)

        # --- Legacy parity (all pairs, all models) ---

        def test_identity_matches_legacy_all_pairs(self):
            """All identity distances match the pure-Python loop exactly."""
            legacy = DistanceCalculator("identity", method="python").get_distance(
                self.msa
            )
            accel = self._calc("identity").get_distance(self.msa)
            for i in range(len(legacy)):
                for j in range(i):
                    self.assertAlmostEqual(
                        legacy[i, j],
                        accel[i, j],
                        msg=f"Mismatch at ({i},{j}) identity",
                    )

        def test_blosum62_matches_legacy_all_pairs(self):
            """All BLOSUM62 distances match the pure-Python loop exactly."""
            legacy = DistanceCalculator("blosum62", method="python").get_distance(
                self.msa
            )
            accel = self._calc("blosum62").get_distance(self.msa)
            for i in range(len(legacy)):
                for j in range(i):
                    self.assertAlmostEqual(
                        legacy[i, j],
                        accel[i, j],
                        msg=f"Mismatch at ({i},{j}) blosum62",
                    )

        # --- Both alignment object types ---

        def test_new_alignment_object(self):
            """Bio.Align.Alignment input produces correct distances."""
            dm = self._calc("identity").get_distance(self.aln)
            self.assertAlmostEqual(dm["Alpha", "Beta"], 1 - 10 / 13)

        # --- Gap / skip_letters handling ---

        def test_identity_gaps_compared_as_characters(self):
            """Identity model treats gaps as normal characters."""
            gapped = AlignIO.read(StringIO(">A\nA-A--\n>B\nAAA--\n"), "fasta")
            legacy = DistanceCalculator("identity", method="python").get_distance(
                gapped
            )
            accel = self._calc("identity").get_distance(gapped)
            self.assertAlmostEqual(legacy["A", "B"], accel["A", "B"])
            # 4 matches (A,A,-,-) out of 5 columns
            self.assertAlmostEqual(accel["A", "B"], 1 - 4 / 5)

        def test_scoring_model_skips_gaps(self):
            """Scoring models skip '-' and '*' characters."""
            gapped = AlignIO.read(StringIO(">Alpha\nA-A--\n>Gamma\n-Y-Y-\n"), "fasta")
            legacy = DistanceCalculator("blosum62", method="python").get_distance(
                gapped
            )
            accel = self._calc("blosum62").get_distance(gapped)
            self.assertAlmostEqual(legacy["Alpha", "Gamma"], accel["Alpha", "Gamma"])
            # No overlapping valid positions -> distance is 1.0
            self.assertEqual(accel["Alpha", "Gamma"], 1.0)

        # --- Subset-based validation ---

        def test_subset_equivalence(self):
            """Sub-matrix entries match those in the full matrix."""
            full_dm = self._calc("blosum62").get_distance(self.msa)
            sub_seqs = [self.msa[0], self.msa[1], self.msa[2]]
            sub_msa = AlignIO.read(
                StringIO(
                    f">{sub_seqs[0].id}\n{sub_seqs[0].seq}\n"
                    f">{sub_seqs[1].id}\n{sub_seqs[1].seq}\n"
                    f">{sub_seqs[2].id}\n{sub_seqs[2].seq}\n"
                ),
                "fasta",
            )
            sub_dm = self._calc("blosum62").get_distance(sub_msa)
            for a, b in [("Alpha", "Beta"), ("Alpha", "Gamma"), ("Beta", "Gamma")]:
                self.assertAlmostEqual(
                    full_dm[a, b],
                    sub_dm[a, b],
                    msg=f"Subset mismatch for ({a},{b})",
                )

        # --- Integration: NJ and UPGMA tree construction ---

        def test_nj_tree_topology(self):
            """Distance matrix feeds correctly into Neighbor-Joining."""
            dm = self._calc("identity").get_distance(self.msa)
            constructor = DistanceTreeConstructor()
            tree = constructor.nj(dm)
            terminals = tree.get_terminals()
            self.assertEqual(len(terminals), 5)
            terminal_names = sorted(t.name for t in terminals)
            self.assertEqual(
                terminal_names, ["Alpha", "Beta", "Delta", "Epsilon", "Gamma"]
            )

        def test_upgma_tree_topology(self):
            """Distance matrix feeds correctly into UPGMA."""
            dm = self._calc("identity").get_distance(self.msa)
            constructor = DistanceTreeConstructor()
            tree = constructor.upgma(dm)
            terminals = tree.get_terminals()
            self.assertEqual(len(terminals), 5)
            terminal_names = sorted(t.name for t in terminals)
            self.assertEqual(
                terminal_names, ["Alpha", "Beta", "Delta", "Epsilon", "Gamma"]
            )

        def test_nj_matches_legacy_tree(self):
            """NJ tree from accelerated DM matches tree from legacy DM."""
            legacy_dm = DistanceCalculator("identity", method="python").get_distance(
                self.msa
            )
            accel_dm = self._calc("identity").get_distance(self.msa)
            constructor = DistanceTreeConstructor()
            legacy_tree = constructor.nj(legacy_dm)
            accel_tree = constructor.nj(accel_dm)
            legacy_terms = sorted(legacy_tree.get_terminals(), key=lambda t: t.name)
            accel_terms = sorted(accel_tree.get_terminals(), key=lambda t: t.name)
            for lt, at in zip(legacy_terms, accel_terms):
                self.assertEqual(lt.name, at.name)
                self.assertAlmostEqual(
                    lt.branch_length,
                    at.branch_length,
                    msg=f"Branch length mismatch for {lt.name}",
                )

    _MethodTest.__name__ = f"{method_name.capitalize()}StrategyTest"
    _MethodTest.__qualname__ = _MethodTest.__name__
    return _MethodTest


# Conditionally define test classes for each method
if HAS_NUMPY:
    NumpyStrategyTest = _build_method_test_class("numpy")
    OneHotStrategyTest = _build_method_test_class("onehot")

if HAS_SCIPY:
    ScipyStrategyTest = _build_method_test_class("scipy")


class ParallelismTest(unittest.TestCase):
    """Test parallel execution paths and environment variable override."""

    @unittest.skipUnless(HAS_NUMPY, "numpy required")
    def test_n_jobs_1_sequential(self):
        """n_jobs=1 forces sequential execution without error."""
        msa = AlignIO.read(MSA_PHY, "phylip")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonExperimentalWarning)
            dm = DistanceCalculator("identity", method="numpy", n_jobs=1).get_distance(
                msa
            )
        self.assertAlmostEqual(dm["Alpha", "Beta"], 1 - 10 / 13)

    @unittest.skipUnless(HAS_NUMPY, "numpy required")
    def test_n_jobs_2_threaded(self):
        """n_jobs=2 runs with threading without error."""
        msa = AlignIO.read(MSA_PHY, "phylip")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonExperimentalWarning)
            dm = DistanceCalculator("identity", method="numpy", n_jobs=2).get_distance(
                msa
            )
        self.assertAlmostEqual(dm["Alpha", "Beta"], 1 - 10 / 13)

    @unittest.skipUnless(HAS_NUMPY, "numpy required")
    def test_env_var_overrides_n_jobs(self):
        """BIOPYTHON_DIST_JOBS environment variable takes precedence."""
        os.environ["BIOPYTHON_DIST_JOBS"] = "1"
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", BiopythonExperimentalWarning)
                calc = DistanceCalculator("identity", method="numpy", n_jobs=4)
            self.assertEqual(calc._strategy._n_jobs, 1)
        finally:
            del os.environ["BIOPYTHON_DIST_JOBS"]


if __name__ == "__main__":
    unittest.main()
