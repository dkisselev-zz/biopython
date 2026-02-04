# Copyright 2017 by Maximilian Greil.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SVDSuperimposer module."""

import unittest
import pytest

pytest.importorskip("numpy")
from numpy import around
from numpy import array
from numpy import array_equal
from numpy import dot  # missing in old PyPy's micronumpy
from numpy.linalg import det  # Missing in PyPy 2.0 numpypy
from numpy.linalg import svd  # Missing in PyPy 2.0 numpypy
from Bio.SVDSuperimposer import SVDSuperimposer


class SVDSuperimposerTest(unittest.TestCase):
    def setUp(self):
        self.x = array(
            [
                [51.65, -1.90, 50.07],
                [50.40, -1.23, 50.65],
                [50.68, -0.04, 51.54],
                [50.22, -0.02, 52.85],
            ]
        )

        self.y = array(
            [
                [51.30, -2.99, 46.54],
                [51.09, -1.88, 47.58],
                [52.36, -1.20, 48.03],
                [52.71, -1.18, 49.38],
            ]
        )

        self.sup = SVDSuperimposer()
        self.sup.set(self.x, self.y)

    def test_get_init_rms(self):
        x = array([[1.19, 1.28, 1.37], [1.46, 1.55, 1.64], [1.73, 1.82, 1.91]])
        y = array([[1.91, 1.82, 1.73], [1.64, 1.55, 1.46], [1.37, 1.28, 1.19]])
        self.sup.set(x, y)
        assert self.sup.init_rms is None
        assert self.sup.get_init_rms() == pytest.approx(0.8049844719, abs=5e-8)

    def test_oldTest(self):
        assert array_equal(
                around(self.sup.reference_coords, decimals=3),
                around(self.x, decimals=3),
            )
        assert array_equal(around(self.sup.coords, decimals=3), around(self.y, decimals=3))
        assert self.sup.rot is None
        assert self.sup.tran is None
        assert self.sup.rms is None
        assert self.sup.init_rms is None

        self.sup.run()
        assert array_equal(
                around(self.sup.reference_coords, decimals=3),
                around(self.x, decimals=3),
            )
        assert array_equal(around(self.sup.coords, decimals=3), around(self.y, decimals=3))
        rot = array(
            [
                [0.68304983, 0.53664371, 0.49543563],
                [-0.52277295, 0.83293229, -0.18147242],
                [-0.51005037, -0.13504564, 0.84947707],
            ]
        )
        tran = array([38.78608157, -20.65451334, -15.42227366])
        assert array_equal(around(self.sup.rot, decimals=3), around(rot, decimals=3))
        assert array_equal(around(self.sup.tran, decimals=3), around(tran, decimals=3))
        assert self.sup.rms is None
        assert self.sup.init_rms is None

        assert self.sup.get_rms() == pytest.approx(0.00304266526014, abs=5e-8)

        rot_get, tran_get = self.sup.get_rotran()
        assert array_equal(around(rot_get, decimals=3), around(rot, decimals=3))
        assert array_equal(around(tran_get, decimals=3), around(tran, decimals=3))

        y_on_x1 = dot(self.y, rot) + tran
        y_x_solution = array(
            [
                [5.16518846e01, -1.90018270e00, 5.00708397e01],
                [5.03977138e01, -1.22877050e00, 5.06488200e01],
                [5.06801788e01, -4.16095666e-02, 5.15368866e01],
                [5.02202228e01, -1.94372374e-02, 5.28534537e01],
            ]
        )
        assert array_equal(around(y_on_x1, decimals=3), around(y_x_solution, decimals=3))

        y_on_x2 = self.sup.get_transformed()
        assert array_equal(around(y_on_x2, decimals=3), around(y_x_solution, decimals=3))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
