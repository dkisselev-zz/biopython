# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for pairwise2 module using the pure Python fallback functions.

This test file imports the TestCases from ``pairwise2_testCases.py``.
If you want to add more tests, do this over there.

"""

import unittest
import pytest

# Import all test classes deliberately:
from pairwise2_testCases import *  # noqa: F401, F403

from Bio import pairwise2

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
