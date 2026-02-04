# Copyright 2019 by Sergio Valqui. All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Unit tests for the Bio.PDB.Dice Module."""

import os
import tempfile
import unittest
import pytest

from Bio.PDB import Dice
from Bio.PDB import PDBParser


class DiceTests(unittest.TestCase):
    """Tests for PDB.Dice module."""

    def test_dice(self):
        """Self test for PDB.Dice module."""
        # Load and dice 2BEG.pdb
        parser = PDBParser()
        pdb_file = "PDB/2BEG.pdb"
        structure = parser.get_structure("scr", pdb_file)
        file_number, file_pdb_diced = tempfile.mkstemp()
        os.close(file_number)
        # dice 2BEG.pdb, chain "B" from residue 18 to 20"
        Dice.extract(structure, "B", 18, 20, file_pdb_diced)
        try:
            # Load the diced file
            diced_str = parser.get_structure("2beg_diced", file_pdb_diced)
            l_diced_cha = list(diced_str.get_chains())
            l_diced_res = list(diced_str.get_residues())
            l_diced_ato = list(diced_str.get_atoms())

            # Chain checks
            assert len(l_diced_cha) == 1
            assert l_diced_cha[0].id == "B"

            # Residue Checks
            assert len(l_diced_res) == 3
            assert l_diced_res[0].id[1] == 18
            assert l_diced_res[2].id[1] == 20

            # Atom checks
            assert len(l_diced_ato) == 29
            assert l_diced_ato[0].name == "N"
            assert l_diced_ato[0].parent.resname == "VAL"
            assert l_diced_ato[0].parent.parent.id == "B"
            assert l_diced_ato[28].name == "CZ"
            assert l_diced_ato[28].parent.resname == "PHE"
            assert l_diced_ato[28].parent.parent.id == "B"

        finally:
            os.remove(file_pdb_diced)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
