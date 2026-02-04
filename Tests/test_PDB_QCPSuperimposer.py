# Copyright 2017 by Maximilian Greil.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for QCPSuperimposer module."""

import unittest
import pytest

np = pytest.importorskip("numpy")
from Bio.PDB import PDBParser
from Bio.PDB import Selection
from Bio.PDB.qcprot import QCPSuperimposer
from Bio.SVDSuperimposer import SVDSuperimposer


class QCPSuperimposerTest(unittest.TestCase):
    def setUp(self):
        self.x = np.array(
            [
                [51.65, -1.90, 50.07],
                [50.40, -1.23, 50.65],
                [50.68, -0.04, 51.54],
                [50.22, -0.02, 52.85],
            ]
        )

        self.y = np.array(
            [
                [51.30, -2.99, 46.54],
                [51.09, -1.88, 47.58],
                [52.36, -1.20, 48.03],
                [52.71, -1.18, 49.38],
            ]
        )

    # Public methods
    def test_set(self):
        """Test setting of initial parameters."""
        sup = QCPSuperimposer()
        sup.set(self.x, self.y)

        assert np.allclose(sup.reference_coords, self.x, atol=1e-6)
        assert np.allclose(sup.coords, self.y, atol=1e-6)
        assert sup.transformed_coords is None
        assert sup.rot is None
        assert sup.tran is None
        assert sup.rms is None
        assert sup.init_rms is None

    def test_run(self):
        """Test QCP on dummy data."""
        sup = QCPSuperimposer()
        sup.set(self.x, self.y)

        sup.run()
        assert np.allclose(sup.reference_coords, self.x, atol=1e-6)
        assert np.allclose(sup.coords, self.y, atol=1e-6)
        assert sup.transformed_coords is None

        calc_rot = [
            [0.683, 0.537, 0.495],
            [-0.523, 0.833, -0.181],
            [-0.510, -0.135, 0.849],
        ]
        calc_tran = [38.786, -20.655, -15.422]

        assert np.allclose(np.array(calc_rot), sup.rot, atol=1e-3)
        assert np.allclose(np.array(calc_tran), sup.tran, atol=1e-3)

        # We can reduce precision here since we do a similar calculation
        # for a full structure down below.
        assert sup.rms == pytest.approx(0.003, abs=0.0005)
        assert sup.init_rms is None

    def test_compare_to_svd(self):
        """Compare results of QCP to SVD."""
        sup = QCPSuperimposer()
        sup.set(self.x, self.y)
        sup.run()

        svd_sup = SVDSuperimposer()
        svd_sup.set(self.x, self.y)
        svd_sup.run()

        assert svd_sup.get_rms() == pytest.approx(sup.rms, abs=0.0005)
        assert np.allclose(svd_sup.rot, sup.rot, atol=1e-3)
        assert np.allclose(svd_sup.tran, sup.tran, atol=1e-3)
        assert np.allclose(svd_sup.get_transformed(), sup.get_transformed(), atol=1e-3)

    def test_compare_to_svd_lines(self):
        """Compare results of QCP to SVD using simple lines."""
        ref = np.linspace([0, 0, 0], [10, 10, 10], 5)
        mob = np.linspace([20, 10, 30], [20, 20, 30], 5)

        sup = QCPSuperimposer()
        sup.set(ref, mob)
        sup.run()

        svd_sup = SVDSuperimposer()
        svd_sup.set(ref, mob)
        svd_sup.run()

        assert svd_sup.get_rms() == pytest.approx(sup.rms, abs=0.0005)
        assert np.allclose(svd_sup.get_transformed(), sup.get_transformed(), atol=1e-3)

    def test_get_transformed(self):
        """Test transformation of coordinates after QCP."""
        sup = QCPSuperimposer()
        sup.set(self.x, self.y)

        sup.run()
        transformed_coords = [
            [51.652, -1.900, 50.071],
            [50.398, -1.229, 50.649],
            [50.680, -0.042, 51.537],
            [50.220, -0.019, 52.853],
        ]

        assert np.allclose(sup.get_transformed(), np.array(transformed_coords), atol=1e-3)

    def test_get_init_rms(self):
        """Test initial RMS calculation."""
        x = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        y = np.array([[2.0, 0.0, 0.0], [1.0, 1.0, 0.0], [1.0, 0.0, 1.0]])

        sup = QCPSuperimposer()
        sup.set(x, y)
        assert sup.init_rms is None

        expected_init_rms = 1.0
        assert sup.get_init_rms() == pytest.approx(expected_init_rms, abs=5e-07)

    def test_on_pdb(self):
        """Align a PDB to itself."""
        pdb1 = "PDB/1A8O.pdb"
        p = PDBParser()
        s1 = p.get_structure("FIXED", pdb1)
        fixed = Selection.unfold_entities(s1, "A")
        s2 = p.get_structure("MOVING", pdb1)
        moving = Selection.unfold_entities(s2, "A")

        rot = np.eye(3, dtype=np.float64)
        tran = np.array([1.0, 2.0, 3.0], dtype=np.float64)
        for atom in moving:
            atom.transform(rot, tran)

        sup = QCPSuperimposer()
        sup.set_atoms(fixed, moving)
        assert np.allclose(sup.rotran[0], rot, atol=1e-3)
        assert np.allclose(sup.rotran[1], -tran, atol=1e-3)
        assert sup.rms == pytest.approx(0.0, abs=5e-07)

    def test_compare_rmsd_to_transformed(self):
        """Compare RMSD from QCP algorithm to that from transformed"""
        ref = np.linspace([0, 0, 0], [10, 10, 10], 5)
        mob = np.linspace([20, 10, 30], [20, 20, 30], 5)

        sup = QCPSuperimposer()
        sup.set(ref, mob)
        sup.run()
        rms = sup.get_rms()
        mob_fitted = sup.get_transformed()
        rms_fitted = np.sqrt(((ref - mob_fitted) ** 2).sum() / ref.shape[0])
        assert rms == pytest.approx(rms_fitted, abs=5e-07)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
