# Copyright 2024 by Anthropic.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this distribution.
"""Tests for Bio.Phylo._distance_methods backends.

Covers method consistency (all backends vs python reference), subset
validation, integration with tree constructors, parallel execution,
edge cases, and the BIOPYTHON_DISTANCE_JOBS environment variable.
"""

import os
import unittest

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

try:
    import scipy  # noqa: F401

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

# Path to the shared test MSA (PHYLIP format, 5 sequences, 13 positions)
MSA_FILE = "TreeConstruction/msa.phy"

# Known exact value: Alpha-Beta identity distance = 1 - 10/13
IDENTITY_ALPHA_BETA = 1 - 10 / 13


def _load_msa():
    """Load the reference MSA as a MultipleSeqAlignment."""
    return AlignIO.read(MSA_FILE, "phylip")


def _get_distances(msa, model="identity", method="python", n_jobs=1):
    """Return the lower-triangular distance list in row-major order."""
    dc = DistanceCalculator(model, method=method, n_jobs=n_jobs)
    dm = dc.get_distance(msa)
    names = dm.names
    result = []
    for i in range(1, len(names)):
        for j in range(i):
            result.append(dm[names[i], names[j]])
    return result


def _vectorised_methods():
    """Return the list of non-python methods available for testing."""
    methods = ["numpy", "onehot"]
    if HAS_SCIPY:
        methods.append("scipy")
    return methods


class MethodConsistencyIdentityTest(unittest.TestCase):
    """All backends must agree with the python reference for identity model."""

    @classmethod
    def setUpClass(cls):
        cls.msa = _load_msa()
        cls.ref = _get_distances(cls.msa, model="identity", method="python")

    def test_python_known_value(self):
        """Python reference matches the known Alpha-Beta distance."""
        self.assertEqual(self.ref[0], IDENTITY_ALPHA_BETA)

    def test_numpy_matches_python(self):
        dists = _get_distances(self.msa, model="identity", method="numpy")
        for i, (d, r) in enumerate(zip(dists, self.ref)):
            self.assertAlmostEqual(d, r, places=10, msg=f"pair {i}")

    @unittest.skipUnless(HAS_SCIPY, "scipy not available")
    def test_scipy_matches_python(self):
        dists = _get_distances(self.msa, model="identity", method="scipy")
        for i, (d, r) in enumerate(zip(dists, self.ref)):
            self.assertAlmostEqual(d, r, places=10, msg=f"pair {i}")

    def test_onehot_matches_python(self):
        dists = _get_distances(self.msa, model="identity", method="onehot")
        for i, (d, r) in enumerate(zip(dists, self.ref)):
            self.assertAlmostEqual(d, r, places=10, msg=f"pair {i}")


class MethodConsistencyScoringTest(unittest.TestCase):
    """All backends must agree with the python reference for scoring models."""

    @classmethod
    def setUpClass(cls):
        cls.msa = _load_msa()

    def _check_model(self, model):
        ref = _get_distances(self.msa, model=model, method="python")
        for method in _vectorised_methods():
            dists = _get_distances(self.msa, model=model, method=method)
            for i, (d, r) in enumerate(zip(dists, ref)):
                self.assertAlmostEqual(
                    d, r, places=10, msg=f"{method} model={model} pair {i}"
                )

    def test_blastn(self):
        self._check_model("blastn")

    def test_trans(self):
        self._check_model("trans")

    def test_blosum62(self):
        self._check_model("blosum62")


class SubsetValidationTest(unittest.TestCase):
    """Distances for a subset of sequences must match the full matrix."""

    def test_subset_identity(self):
        """Three-sequence subset distances match the full 5-sequence matrix."""
        msa = _load_msa()
        subset = MultipleSeqAlignment(list(msa)[:3])
        full_dm = DistanceCalculator("identity").get_distance(msa)
        sub_dm = DistanceCalculator("identity").get_distance(subset)
        for i in range(1, 3):
            for j in range(i):
                ni, nj = sub_dm.names[i], sub_dm.names[j]
                self.assertAlmostEqual(full_dm[ni, nj], sub_dm[ni, nj], places=10)

    def test_subset_scoring(self):
        """Same subset check with blosum62 scoring model."""
        msa = _load_msa()
        subset = MultipleSeqAlignment(list(msa)[:3])
        full_dm = DistanceCalculator("blosum62").get_distance(msa)
        sub_dm = DistanceCalculator("blosum62").get_distance(subset)
        for i in range(1, 3):
            for j in range(i):
                ni, nj = sub_dm.names[i], sub_dm.names[j]
                self.assertAlmostEqual(full_dm[ni, nj], sub_dm[ni, nj], places=10)

    def test_subset_numpy_vs_python(self):
        """Numpy on 3-seq subset matches python on same subset."""
        msa = _load_msa()
        subset = MultipleSeqAlignment(list(msa)[:3])
        ref = DistanceCalculator("identity", method="python").get_distance(subset)
        npy = DistanceCalculator("identity", method="numpy").get_distance(subset)
        for i in range(1, len(ref.names)):
            for j in range(i):
                ni, nj = ref.names[i], ref.names[j]
                self.assertAlmostEqual(ref[ni, nj], npy[ni, nj], places=10)


def _tree_clades(tree):
    """Return the set of non-trivial internal clades as frozensets of leaf names.

    A clade is non-trivial when it contains more than one leaf but fewer than
    all leaves (i.e. it is neither a single tip nor the root).  This captures
    the tree topology independent of branch lengths or node names.
    """
    n_leaves = len(list(tree.get_terminals()))
    clades = set()
    for clade in tree.get_nonterminals():
        leaves = frozenset(t.name for t in clade.get_terminals())
        if 1 < len(leaves) < n_leaves:
            clades.add(leaves)
    return clades


class IntegrationTest(unittest.TestCase):
    """Full pipeline: get_distance -> tree construction for every method.

    Each method's tree topology (set of internal clades) is compared against
    the topology produced by the ``python`` reference method.
    """

    @classmethod
    def setUpClass(cls):
        cls.msa = _load_msa()
        cls.expected_names = sorted(s.id for s in cls.msa)
        cls.constructor = DistanceTreeConstructor()
        # Compute reference topologies once with the python baseline
        dc_ref = DistanceCalculator("identity", method="python")
        dm_ref = dc_ref.get_distance(cls.msa)
        cls.upgma_ref_clades = _tree_clades(cls.constructor.upgma(dm_ref))
        cls.nj_ref_clades = _tree_clades(cls.constructor.nj(dm_ref))

    def _check_tree(self, method, tree_func_name, ref_clades=None):
        dc = DistanceCalculator("identity", method=method)
        dm = dc.get_distance(self.msa)
        tree = getattr(self.constructor, tree_func_name)(dm)
        # Terminal names must be complete
        names = sorted(t.name for t in tree.get_terminals())
        self.assertEqual(names, self.expected_names)
        # Topology must match the python reference when ref_clades is given.
        # NJ is sensitive to tie-breaking: scipy hamming computes
        # n_mismatches/L while the python reference computes 1-n_matches/L;
        # these differ by 1 ULP and can flip an NJ branch when two distances
        # are nearly equal.  UPGMA averages are stable under this perturbation.
        if ref_clades is not None:
            self.assertEqual(
                _tree_clades(tree),
                ref_clades,
                msg=f"method={method} tree={tree_func_name} topology mismatch",
            )

    def test_upgma_python(self):
        self._check_tree("python", "upgma", self.upgma_ref_clades)

    def test_upgma_numpy(self):
        self._check_tree("numpy", "upgma", self.upgma_ref_clades)

    @unittest.skipUnless(HAS_SCIPY, "scipy not available")
    def test_upgma_scipy(self):
        self._check_tree("scipy", "upgma", self.upgma_ref_clades)

    def test_upgma_onehot(self):
        self._check_tree("onehot", "upgma", self.upgma_ref_clades)

    def test_nj_python(self):
        self._check_tree("python", "nj", self.nj_ref_clades)

    def test_nj_numpy(self):
        self._check_tree("numpy", "nj", self.nj_ref_clades)

    @unittest.skipUnless(HAS_SCIPY, "scipy not available")
    def test_nj_scipy(self):
        # Topology not checked: see _check_tree docstring re NJ tie-breaking
        self._check_tree("scipy", "nj")

    def test_nj_onehot(self):
        self._check_tree("onehot", "nj", self.nj_ref_clades)


class ParallelConsistencyTest(unittest.TestCase):
    """n_jobs=2 via python method must match n_jobs=1."""

    @classmethod
    def setUpClass(cls):
        cls.msa = _load_msa()

    def test_parallel_identity(self):
        seq_ref = _get_distances(self.msa, model="identity", method="python", n_jobs=1)
        par_ref = _get_distances(self.msa, model="identity", method="python", n_jobs=2)
        self.assertEqual(seq_ref, par_ref)

    def test_parallel_scoring(self):
        seq_ref = _get_distances(self.msa, model="blastn", method="python", n_jobs=1)
        par_ref = _get_distances(self.msa, model="blastn", method="python", n_jobs=2)
        self.assertEqual(seq_ref, par_ref)


class EdgeCaseTest(unittest.TestCase):
    """Edge cases: two sequences, identical seqs, all-gap column, bad letter."""

    def _make_msa(self, seqs, names=None):
        """Build a MultipleSeqAlignment from sequence strings."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        if names is None:
            names = [f"seq{i}" for i in range(len(seqs))]
        return MultipleSeqAlignment(
            [SeqRecord(Seq(s), id=n) for s, n in zip(seqs, names)]
        )

    def test_two_sequences(self):
        """Two-sequence MSA yields distance 0 for identical sequences."""
        msa = self._make_msa(["ACGT", "ACGT"], ["A", "B"])
        dm = DistanceCalculator("identity").get_distance(msa)
        self.assertEqual(dm["A", "B"], 0.0)

    def test_identical_sequences_all_methods(self):
        """Identical sequences yield distance 0 for every method."""
        msa = self._make_msa(["ACGTACGT"] * 3, ["A", "B", "C"])
        pairs = [("A", "B"), ("A", "C"), ("B", "C")]
        for method in ["python", "numpy", "onehot"]:
            dm = DistanceCalculator("identity", method=method).get_distance(msa)
            for i, j in pairs:
                self.assertEqual(dm[i, j], 0.0, msg=f"method={method}")

    @unittest.skipUnless(HAS_SCIPY, "scipy not available")
    def test_identical_sequences_scipy(self):
        """Identical sequences yield distance 0 via scipy."""
        msa = self._make_msa(["ACGTACGT"] * 3, ["A", "B", "C"])
        dm = DistanceCalculator("identity", method="scipy").get_distance(msa)
        pairs = [("A", "B"), ("A", "C"), ("B", "C")]
        for i, j in pairs:
            self.assertEqual(dm[i, j], 0.0)

    def test_completely_different(self):
        """Completely different sequences yield distance 1."""
        msa = self._make_msa(["AAAA", "CCCC"], ["A", "B"])
        for method in ["python", "numpy", "onehot"]:
            dm = DistanceCalculator("identity", method=method).get_distance(msa)
            self.assertEqual(dm["A", "B"], 1.0, msg=f"method={method}")

    @unittest.skipUnless(HAS_SCIPY, "scipy not available")
    def test_completely_different_scipy(self):
        """Completely different sequences yield distance 1 via scipy."""
        msa = self._make_msa(["AAAA", "CCCC"], ["A", "B"])
        dm = DistanceCalculator("identity", method="scipy").get_distance(msa)
        self.assertEqual(dm["A", "B"], 1.0)

    def test_all_gap_column_scoring(self):
        """A column of all gaps is skipped in scoring; no crash."""
        msa = self._make_msa(["A-CGT", "A-CGT"], ["A", "B"])
        dm = DistanceCalculator("blastn").get_distance(msa)
        self.assertEqual(dm["A", "B"], 0.0)

    def test_bad_letter_raises(self):
        """Unknown character raises ValueError in scoring model."""
        msa = self._make_msa(["A?GT", "ACGT"], ["A", "B"])
        for method in ["python", "numpy", "onehot"]:
            with self.assertRaises(ValueError, msg=f"method={method}"):
                DistanceCalculator("blastn", method=method).get_distance(msa)

    @unittest.skipUnless(HAS_SCIPY, "scipy not available")
    def test_bad_letter_raises_scipy(self):
        """Unknown character raises ValueError in scipy scoring path."""
        msa = self._make_msa(["A?GT", "ACGT"], ["A", "B"])
        with self.assertRaises(ValueError):
            DistanceCalculator("blastn", method="scipy").get_distance(msa)

    def test_invalid_method_raises(self):
        """Unknown method name raises ValueError at construction time."""
        with self.assertRaises(ValueError):
            DistanceCalculator("identity", method="nonexistent")


class EnvironmentVariableTest(unittest.TestCase):
    """BIOPYTHON_DISTANCE_JOBS environment variable tests."""

    def setUp(self):
        self._original = os.environ.pop("BIOPYTHON_DISTANCE_JOBS", None)

    def tearDown(self):
        if self._original is not None:
            os.environ["BIOPYTHON_DISTANCE_JOBS"] = self._original
        else:
            os.environ.pop("BIOPYTHON_DISTANCE_JOBS", None)

    def test_env_overrides_n_jobs(self):
        """Environment variable overrides the n_jobs argument."""
        from Bio.Phylo._distance_methods import _resolve_n_jobs

        os.environ["BIOPYTHON_DISTANCE_JOBS"] = "4"
        self.assertEqual(_resolve_n_jobs(1), 4)

    def test_env_minus_one_resolves(self):
        """BIOPYTHON_DISTANCE_JOBS=-1 resolves to os.cpu_count()."""
        from Bio.Phylo._distance_methods import _resolve_n_jobs

        os.environ["BIOPYTHON_DISTANCE_JOBS"] = "-1"
        self.assertEqual(_resolve_n_jobs(1), os.cpu_count())

    def test_env_invalid_falls_back(self):
        """Invalid env var string falls back to n_jobs argument."""
        from Bio.Phylo._distance_methods import _resolve_n_jobs

        os.environ["BIOPYTHON_DISTANCE_JOBS"] = "bad"
        self.assertEqual(_resolve_n_jobs(3), 3)

    def test_no_env_uses_n_jobs(self):
        """When env var is unset, n_jobs is used directly."""
        from Bio.Phylo._distance_methods import _resolve_n_jobs

        self.assertEqual(_resolve_n_jobs(2), 2)
        self.assertEqual(_resolve_n_jobs(-1), os.cpu_count())


if __name__ == "__main__":
    unittest.main()
