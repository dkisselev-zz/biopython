# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for KeyWList."""

import os.path
import unittest
import pytest

from Bio.SwissProt import KeyWList


class KeyWListTest(unittest.TestCase):
    """Tests for KeyWList."""

    def test_parse(self):
        """Test parsing keywlist.txt works."""
        filename = os.path.join("SwissProt", "keywlist.txt")
        with open(filename) as handle:
            records = KeyWList.parse(handle)

            # Testing the first record
            record = next(records)
            assert record["ID"] == "2Fe-2S."
            assert record["AC"] == "KW-0001"
            assert record["DE"] == "Protein which contains at least one 2Fe-2S iron-sulfur cluster: 2 iron atoms complexed to 2 inorganic sulfides and 4 sulfur atoms of cysteines from the protein."
            assert record["SY"] == "Fe2S2; [2Fe-2S] cluster; [Fe2S2] cluster; Fe2/S2 (inorganic) cluster; Di-mu-sulfido-diiron; 2 iron, 2 sulfur cluster binding."
            assert len(record["GO"]) == 1
            assert record["GO"] == ["GO:0051537; 2 iron, 2 sulfur cluster binding"]
            assert len(record["HI"]) == 2
            assert record["HI"][0] == "Ligand: Iron; Iron-sulfur; 2Fe-2S."
            assert record["HI"][1] == "Ligand: Metal-binding; 2Fe-2S."
            assert record["CA"] == "Ligand."

            # Testing the second record
            record = next(records)
            assert record["IC"] == "Molecular function."
            assert record["AC"] == "KW-9992"
            assert record["DE"] == "Keywords assigned to proteins due to their particular molecular function."

            # Testing the third record
            record = next(records)
            assert record["ID"] == "Zymogen."
            assert record["AC"] == "KW-0865"
            assert record["DE"] == "The enzymatically inactive precursor of mostly proteolytic enzymes."
            assert record["SY"] == "Proenzyme."
            assert len(record["HI"]) == 1
            assert record["HI"][0] == "PTM: Zymogen."
            assert record["CA"] == "PTM."

    def test_parse2(self):
        """Parsing keywlist2.txt (without header and footer)."""
        filename = os.path.join("SwissProt", "keywlist2.txt")
        with open(filename) as handle:
            records = KeyWList.parse(handle)

            # Testing the first record
            record = next(records)
            assert record["ID"] == "2Fe-2S."
            assert record["AC"] == "KW-0001"
            assert record["DE"] == "Protein which contains at least one 2Fe-2S iron-sulfur cluster: 2 iron atoms complexed to 2 inorganic sulfides and 4 sulfur atoms of cysteines from the protein."
            assert record["SY"] == "Fe2S2; [2Fe-2S] cluster; [Fe2S2] cluster; Fe2/S2 (inorganic) cluster; Di-mu-sulfido-diiron; 2 iron, 2 sulfur cluster binding."
            assert len(record["GO"]) == 1
            assert record["GO"] == ["GO:0051537; 2 iron, 2 sulfur cluster binding"]
            assert len(record["HI"]) == 2
            assert record["HI"][0] == "Ligand: Iron; Iron-sulfur; 2Fe-2S."
            assert record["HI"][1] == "Ligand: Metal-binding; 2Fe-2S."
            assert record["CA"] == "Ligand."

            # Testing the second record
            record = next(records)
            assert record["ID"] == "3D-structure."
            assert record["AC"] == "KW-0002"
            assert record["DE"] == "Protein, or part of a protein, whose three-dimensional structure has been resolved experimentally (for example by X-ray crystallography or NMR spectroscopy) and whose coordinates are available in the PDB database. Can also be used for theoretical models."
            assert len(record["HI"]) == 1
            assert record["HI"][0] == "Technical term: 3D-structure."
            assert record["CA"] == "Technical term."

            # Testing the third record
            record = next(records)
            assert record["ID"] == "3Fe-4S."
            assert record["AC"] == "KW-0003"
            assert record["DE"] == "Protein which contains at least one 3Fe-4S iron-sulfur cluster: 3 iron atoms complexed to 4 inorganic sulfides and 3 sulfur atoms of cysteines from the protein. In a number of iron-sulfur proteins, the 4Fe-4S cluster can be reversibly converted by oxidation and loss of one iron ion to a 3Fe-4S cluster."
            assert record["SY"] == ""
            assert len(record["GO"]) == 1
            assert record["GO"] == ["GO:0051538; 3 iron, 4 sulfur cluster binding"]
            assert len(record["HI"]) == 2
            assert record["HI"][0] == "Ligand: Iron; Iron-sulfur; 3Fe-4S."
            assert record["HI"][1] == "Ligand: Metal-binding; 3Fe-4S."
            assert record["CA"] == "Ligand."


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
