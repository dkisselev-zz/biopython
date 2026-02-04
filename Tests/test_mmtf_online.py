# Copyright 2017 by Peter Cock.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Tests for mmtf online module."""

import unittest
import pytest
import warnings


from Bio.PDB.mmtf import MMTFParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

pytestmark = pytest.mark.online



class OnlineMMTF(unittest.TestCase):
    """Online tests for the MMTF code."""

    def test_from_url(self):
        """Check parser can fetch a record from its PDB ID."""
        parser = MMTFParser()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            struct = parser.get_structure_from_url("4ZHL")
        atoms = list(struct.get_atoms())
        assert len(atoms) == 2080


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
