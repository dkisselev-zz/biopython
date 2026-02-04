# Copyright (C) 2011 by Brandon Invergo (b.invergo@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for Chi2 module."""

import unittest
import pytest

from Bio.Phylo.PAML import chi2


class ModTest(unittest.TestCase):
    def test_cdf_chi2(self):
        with pytest.raises(ValueError):
            chi2.cdf_chi2(df=0, stat=3.84)
        with pytest.raises(ValueError):
            chi2.cdf_chi2(df=1, stat=-3.84)
        with pytest.raises(TypeError):
            chi2.cdf_chi2(df="d", stat="stat")
        assert chi2.cdf_chi2(2, 3.84) == pytest.approx(0.1466070, abs=5e-06)

    def test_ln_gamma(self):
        with pytest.raises(ValueError):
            chi2._ln_gamma_function(-1)
        assert chi2._ln_gamma_function(10) == pytest.approx(12.80183, abs=5e-06)

    def test_incomplete_gamma(self):
        with pytest.raises(ValueError):
            chi2._incomplete_gamma(x=0.5, alpha=-1)
        assert chi2._incomplete_gamma(0.5, 0.5) == pytest.approx(0.6826895, abs=5e-06)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
