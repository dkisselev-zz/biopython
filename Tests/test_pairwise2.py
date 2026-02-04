# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for pairwise2 module using the default C functions.

This test file imports the TestCases from ``pairwise2_testCases.py``.
If you want to add more tests, do this over there.

"""

import unittest
import pytest

# Import all test classes deliberately:
from pairwise2_testCases import *  # noqa: F401, F403

# Implicitly using functions from C extension:
from Bio import pairwise2

if pairwise2.rint == pairwise2._python_rint:

    pytest.skip("Missing or non-compiled file: 'cpairwise2'", allow_module_level=True)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
