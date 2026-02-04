# Copyright 2014-2016 Marco Galardini.  All rights reserved.
# Adapted from test_Mymodule.py by Jeff Chang
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Tests for the Bio.phenotype module's fitting functionality."""

import pytest

try:
    import numpy as np

    del np
except ImportError:

    pytest.skip("Install NumPy if you want to use Bio.phenotype.", allow_module_level=True)
try:
    import scipy

    del scipy
    from scipy.optimize import OptimizeWarning
except ImportError:

    pytest.skip("Install SciPy if you want to use Bio.phenotype fit functionality.", allow_module_level=True)

import json
import unittest
import warnings

from Bio import BiopythonExperimentalWarning

with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonExperimentalWarning)
    from Bio import phenotype

# Example plate files
JSON_PLATE = "phenotype/Plate.json"


class TestPhenoMicro(unittest.TestCase):
    """Tests for phenotype module."""

    def test_WellRecord(self):
        """Test basic functionalities of WellRecord objects."""
        with open(JSON_PLATE) as handle:
            p = json.load(handle)

        times = p["measurements"]["Hour"]
        w = phenotype.phen_micro.WellRecord(
            "A10",
            signals={times[i]: p["measurements"]["A10"][i] for i in range(len(times))},
        )

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", OptimizeWarning)
            w.fit()
        assert w.area == pytest.approx(20879.5, abs=5e-8)
        assert w.model == "gompertz"
        assert w.lag == pytest.approx(6.0425868725090357, abs=5e-06)
        assert w.plateau == pytest.approx(188.51404344898586, abs=5e-05)
        assert w.slope == pytest.approx(48.190618284831132, abs=0.0005)
        assert w.v == pytest.approx(0.10000000000000001, abs=5e-06)
        assert w.y0 == pytest.approx(45.879770069807989, abs=5e-05)
        assert w.max == 313.0
        assert w.min == 29.0
        assert w.average_height == 217.82552083333334


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
