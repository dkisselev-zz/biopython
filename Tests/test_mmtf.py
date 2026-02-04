# Copyright 2016 by Anthony Bradley.  All rights reserved.
# Revisions copyright 2017 by Peter Cock.  All rights reserved.
# Revisions copyright 2019 by Joe Greener.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Tests for mmtf module."""

import os
import tempfile
import unittest
import pytest
import warnings

mmtf = pytest.importorskip("mmtf")
from Bio.PDB import PDBParser
from Bio.PDB import Select
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.mmtf import MMTFIO
from Bio.PDB.mmtf import MMTFParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning


class ParseMMTF(unittest.TestCase):
    """Testing with real mmtf file(s)."""

    def check_atoms(self):
        """Check all atoms in self.mmtf_atoms and self.mmcif_atoms are equivalent."""
        assert len(self.mmcif_atoms) == len(self.mmtf_atoms)
        for i, e in enumerate(self.mmcif_atoms):
            mmtf_atom = self.mmtf_atoms[i]
            mmcif_atom = self.mmcif_atoms[i]
            assert mmtf_atom.name == mmcif_atom.name  # eg. CA, spaces are removed from atom name
            assert mmtf_atom.fullname == mmcif_atom.fullname  # e.g. " CA ", spaces included
            assert mmtf_atom.coord[0] == pytest.approx(mmcif_atom.coord[0], abs=0.0005)
            assert mmtf_atom.coord[1] == pytest.approx(mmcif_atom.coord[1], abs=0.0005)
            assert mmtf_atom.coord[2] == pytest.approx(mmcif_atom.coord[2], abs=0.0005)
            assert mmtf_atom.bfactor == mmcif_atom.bfactor
            assert mmtf_atom.occupancy == mmcif_atom.occupancy
            assert mmtf_atom.altloc == mmcif_atom.altloc
            assert mmtf_atom.full_id == mmcif_atom.full_id  # (structure id, model id, chain id, residue id, atom id)
            assert mmtf_atom.id == mmcif_atom.name  # id of atom is the atom name (e.g. "CA")
            # self.assertEqual(mmtf_atom.serial_number,mmcif_atom.serial_number) # mmCIF serial number is none
            assert mmtf_atom - mmtf_atom == 0
            assert mmtf_atom - mmcif_atom == 0

    def check_residues(self):
        """Check all residues in self.mmcif_res and self.mmtf_res are equivalent."""
        assert len(self.mmcif_res) == len(self.mmtf_res)
        for i, e in enumerate(self.mmcif_res):
            mmcif_r = self.mmcif_res[i]
            mmtf_r = self.mmtf_res[i]
            assert mmtf_r.level == mmcif_r.level
            assert mmtf_r.disordered == mmcif_r.disordered
            assert mmtf_r.resname == mmcif_r.resname
            assert mmtf_r.segid == mmcif_r.segid
            self.mmcif_atoms = list(mmcif_r.get_atoms())
            self.mmtf_atoms = list(mmtf_r.get_atoms())
            self.check_atoms()

    def check_mmtf_vs_cif(self, mmtf_filename, cif_filename):
        """Compare parsed structures for MMTF and CIF files."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            mmtf_struct = MMTFParser.get_structure(mmtf_filename)
        mmcif_parser = MMCIFParser()
        mmcif_struct = mmcif_parser.get_structure("4CUP", cif_filename)
        self.mmcif_atoms = list(mmcif_struct.get_atoms())
        self.mmtf_atoms = list(mmtf_struct.get_atoms())
        self.check_atoms()
        mmcif_chains = list(mmcif_struct.get_chains())
        mmtf_chains = list(mmtf_struct.get_chains())
        assert len(mmcif_chains) == len(mmtf_chains)
        for i, e in enumerate(mmcif_chains):
            self.mmcif_res = list(mmcif_chains[i].get_residues())
            self.mmtf_res = list(mmtf_chains[i].get_residues())
            self.check_residues()

        self.mmcif_res = list(mmcif_struct.get_residues())
        self.mmtf_res = list(mmtf_struct.get_residues())
        self.check_residues()
        assert sum(1 for _ in mmcif_struct.get_models()) == sum(1 for _ in mmtf_struct.get_models())

    def test_4CUP(self):
        """Compare parsing 4CUP.mmtf and 4CUP.cif."""
        self.check_mmtf_vs_cif("PDB/4CUP.mmtf", "PDB/4CUP.cif")


# TODO:
#    def test_1A8O(self):
#        """Compare parsing 1A8O.mmtf and 1A8O.cif"""
#        self.check_mmtf_vs_cif("PDB/1A8O.mmtf", "PDB/1A8O.cif")

# TODO:
#    def test_4ZHL(self):
#        """Compare parsing 4ZHL.mmtf and 4ZHL.cif"""
#        self.check_mmtf_vs_cif("PDB/4ZHL.mmtf", "PDB/4ZHL.cif")


class SimpleParseMMTF(unittest.TestCase):
    """Just parse some real mmtf files."""

    def test_4ZHL(self):
        """Parse 4ZHL.mmtf."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            structure = MMTFParser.get_structure("PDB/4ZHL.mmtf")

    def test_1A80(self):
        """Parse 1A8O.mmtf."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            structure = MMTFParser.get_structure("PDB/1A8O.mmtf")


class WriteMMTF(unittest.TestCase):
    """Write some MMTF files, read them back in and check them."""

    def test_write(self):
        """Test a simple structure object is written out correctly to MMTF."""
        parser = MMCIFParser()
        structure = parser.get_structure("1A8O", "PDB/1A8O.cif")
        io = MMTFIO()
        io.set_structure(structure)
        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        try:
            io.save(filename)
            struc_back = MMTFParser.get_structure(filename)
            dict_back = mmtf.parse(filename)
            assert dict_back.structure_id == "1A8O"
            assert dict_back.num_models == 1
            assert dict_back.num_chains == 2
            assert dict_back.num_groups == 158
            assert dict_back.num_atoms == 644
            assert len(dict_back.x_coord_list) == 644
            assert len(dict_back.y_coord_list) == 644
            assert len(dict_back.z_coord_list) == 644
            assert len(dict_back.b_factor_list) == 644
            assert len(dict_back.occupancy_list) == 644
            assert dict_back.x_coord_list[5] == 20.022
            assert set(dict_back.ins_code_list) == {"\x00"}
            assert set(dict_back.alt_loc_list) == {"\x00"}
            assert list(dict_back.atom_id_list) == list(range(1, 645))
            assert list(dict_back.sequence_index_list) == list(range(70)) + [-1] * 88
            assert dict_back.chain_id_list == ["A", "B"]
            assert dict_back.chain_name_list == ["A", "A"]
            assert dict_back.chains_per_model == [2]
            assert len(dict_back.group_list) == 21
            assert len(dict_back.group_id_list) == 158
            assert len(dict_back.group_type_list) == 158
            assert dict_back.groups_per_chain == [70, 88]
            assert len(dict_back.entity_list) == 2
            assert dict_back.entity_list[0]["type"] == "polymer"
            assert dict_back.entity_list[0]["chainIndexList"] == [0]
            assert dict_back.entity_list[0]["sequence"] == "MDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPGATLEEMMTACQG"
            assert dict_back.entity_list[1]["type"] == "water"
            assert dict_back.entity_list[1]["chainIndexList"] == [1]
            assert dict_back.entity_list[1]["sequence"] == ""
        finally:
            os.remove(filename)

    def test_multi_model_write(self):
        """Test multiple models are written out correctly to MMTF."""
        parser = PDBParser()
        structure = parser.get_structure("1SSU_mod", "PDB/1SSU_mod.pdb")
        io = MMTFIO()
        io.set_structure(structure)
        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        try:
            io.save(filename)
            struc_back = MMTFParser.get_structure(filename)
            dict_back = mmtf.parse(filename)
            assert dict_back.num_models == 2
            assert dict_back.num_chains == 4
            assert dict_back.num_groups == 4
            assert dict_back.num_atoms == 4
            assert list(dict_back.x_coord_list) == [-1.058, -0.025, 7.024, 6.259]
            assert dict_back.chain_id_list == ["A", "B", "A", "B"]
            assert dict_back.chain_name_list == ["A", "B", "A", "B"]
            assert dict_back.chains_per_model == [2, 2]
            assert len(dict_back.group_list) == 1
            assert len(dict_back.group_id_list) == 4
            assert len(dict_back.group_type_list) == 4
            assert dict_back.groups_per_chain == [1, 1, 1, 1]
            assert len(dict_back.entity_list) == 4
        finally:
            os.remove(filename)

    def test_selection_write(self):
        """Test the use of a Select subclass when writing MMTF files."""
        structure = MMTFParser.get_structure("PDB/4CUP.mmtf")
        io = MMTFIO()
        io.set_structure(structure)
        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)

        class CAonly(Select):
            """Accepts only CA residues."""

            def accept_atom(self, atom):
                if atom.name == "CA" and atom.element == "C":
                    return 1

        try:
            io.save(filename, CAonly())
            struc_back = MMTFParser.get_structure(filename)
            dict_back = mmtf.parse(filename)
            assert dict_back.num_atoms == 116
            assert len(dict_back.x_coord_list) == 116
            assert set(dict_back.alt_loc_list) == {"\x00", "A", "B"}
        finally:
            os.remove(filename)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
