# Copyright 2009-2011 by Eric Talevich.  All rights reserved.
# Revisions copyright 2009-2013 by Peter Cock.  All rights reserved.
# Revisions copyright 2013 Lenna X. Peterson. All rights reserved.
#
# Converted by Eric Talevich from an older unit test copyright 2002
# by Thomas Hamelryck.
#
# Merged related test files into one, by Joao Rodrigues (2020)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Unit tests for the Bio.PDB exposure classes."""

import unittest
import pytest
import warnings

np = pytest.importorskip("numpy")
from Bio.PDB import ExposureCN
from Bio.PDB import HSExposureCA
from Bio.PDB import HSExposureCB
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning


class Exposure(unittest.TestCase):
    """Testing Bio.PDB.HSExposure."""

    def setUp(self):
        pdb_filename = "PDB/a_structure.pdb"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            structure = PDBParser(PERMISSIVE=True).get_structure("X", pdb_filename)
        structure[1].detach_child("B")
        self.model = structure[1]
        # Look at first chain only
        a_residues = list(self.model["A"].child_list)
        assert 86 == len(a_residues)
        assert a_residues[0].get_resname() == "CYS"
        assert a_residues[1].get_resname() == "ARG"
        assert a_residues[2].get_resname() == "CYS"
        assert a_residues[3].get_resname() == "GLY"
        # ...
        assert a_residues[-3].get_resname() == "TYR"
        assert a_residues[-2].get_resname() == "ARG"
        assert a_residues[-1].get_resname() == "CYS"
        self.a_residues = a_residues
        self.radius = 13.0

    def test_HSExposureCA(self):
        """HSExposureCA."""
        _ = HSExposureCA(self.model, self.radius)
        residues = self.a_residues
        assert 0 == len(residues[0].xtra)
        assert 0 == len(residues[1].xtra)
        assert 3 == len(residues[2].xtra)
        assert 0.81250973133184456 == pytest.approx(residues[2].xtra["EXP_CB_PCB_ANGLE"], abs=5e-8)
        assert 14 == residues[2].xtra["EXP_HSE_A_D"]
        assert 14 == residues[2].xtra["EXP_HSE_A_U"]
        assert 3 == len(residues[3].xtra)
        assert 1.3383737 == pytest.approx(residues[3].xtra["EXP_CB_PCB_ANGLE"], abs=5e-8)
        assert 13 == residues[3].xtra["EXP_HSE_A_D"]
        assert 16 == residues[3].xtra["EXP_HSE_A_U"]
        # ...
        assert 3 == len(residues[-2].xtra)
        assert 0.77124014456278489 == pytest.approx(residues[-2].xtra["EXP_CB_PCB_ANGLE"], abs=5e-8)
        assert 24 == residues[-2].xtra["EXP_HSE_A_D"]
        assert 24 == residues[-2].xtra["EXP_HSE_A_U"]
        assert 0 == len(residues[-1].xtra)

    def test_HSExposureCB(self):
        """HSExposureCB."""
        _ = HSExposureCB(self.model, self.radius)
        residues = self.a_residues
        assert 0 == len(residues[0].xtra)
        assert 2 == len(residues[1].xtra)
        assert 20 == residues[1].xtra["EXP_HSE_B_D"]
        assert 5 == residues[1].xtra["EXP_HSE_B_U"]
        assert 2 == len(residues[2].xtra)
        assert 10 == residues[2].xtra["EXP_HSE_B_D"]
        assert 18 == residues[2].xtra["EXP_HSE_B_U"]
        assert 2 == len(residues[3].xtra)
        assert 7 == residues[3].xtra["EXP_HSE_B_D"]
        assert 22 == residues[3].xtra["EXP_HSE_B_U"]
        # ...
        assert 2 == len(residues[-2].xtra)
        assert 14 == residues[-2].xtra["EXP_HSE_B_D"]
        assert 34 == residues[-2].xtra["EXP_HSE_B_U"]
        assert 2 == len(residues[-1].xtra)
        assert 23 == residues[-1].xtra["EXP_HSE_B_D"]
        assert 15 == residues[-1].xtra["EXP_HSE_B_U"]

    def test_ExposureCN(self):
        """HSExposureCN."""
        _ = ExposureCN(self.model, self.radius)
        residues = self.a_residues
        assert 0 == len(residues[0].xtra)
        assert 1 == len(residues[1].xtra)
        assert 25 == residues[1].xtra["EXP_CN"]
        assert 1 == len(residues[2].xtra)
        assert 28 == residues[2].xtra["EXP_CN"]
        assert 1 == len(residues[3].xtra)
        assert 29 == residues[3].xtra["EXP_CN"]
        # ...
        assert 1 == len(residues[-2].xtra)
        assert 48 == residues[-2].xtra["EXP_CN"]
        assert 1 == len(residues[-1].xtra)
        assert 38 == residues[-1].xtra["EXP_CN"]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
