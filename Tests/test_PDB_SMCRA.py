# Copyright 2009-2011 by Eric Talevich.  All rights reserved.
# Revisions copyright 2009-2013 by Peter Cock.  All rights reserved.
# Revisions copyright 2013 Lenna X. Peterson. All rights reserved.
# Revisions copyright 2020 Joao Rodrigues. All rights reserved.
#
# Converted by Eric Talevich from an older unit test copyright 2002
# by Thomas Hamelryck.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Generic unit tests for the SMCRA classes of the Bio.PDB module."""

import unittest
import pytest
import warnings
from copy import deepcopy

np = pytest.importorskip("numpy")
from numpy import dot  # Missing on old PyPy's micronumpy
del dot
from numpy.linalg import det  # Missing in PyPy 2.0 numpypy
from numpy.linalg import svd  # Missing in PyPy 2.0 numpypy
del svd, det
from Bio import BiopythonWarning
from Bio.PDB import Atom
from Bio.PDB import PDBParser
from Bio.PDB import rotmat
from Bio.PDB import Vector
from Bio.PDB.PDBExceptions import PDBConstructionWarning


class Atom_Element(unittest.TestCase):
    """induces Atom Element from Atom Name."""

    def test_atom_element_assignment(self):
        """Atom Element."""
        parser = PDBParser(PERMISSIVE=True, QUIET=True)
        structure = parser.get_structure("X", "PDB/a_structure.pdb")
        residue = structure[0]["A"][("H_PCA", 1, " ")]

        atoms = residue.child_list
        assert "N" == atoms[0].element  # N
        assert "C" == atoms[1].element  # Alpha Carbon
        assert "D" == atoms[4].element  # Deuterium
        assert "CA" == atoms[8].element  # Calcium

    def test_assign_unknown_element(self):
        """Unknown element is assigned 'X'."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            a = Atom.Atom(
                "XE1",
                None,
                None,
                None,
                None,
                " XE1",
                None,  # serial 5170 - 4CP4
            )
        assert a.element == "X"

    def test_ions(self):
        """Element for magnesium is assigned correctly."""
        parser = PDBParser(PERMISSIVE=True)
        structure = parser.get_structure("X", "PDB/ions.pdb")
        # check magnesium atom
        atoms = structure[0]["A"][("H_MG", 1, " ")].child_list
        assert "MG" == atoms[0].element

    def test_hydrogens(self):
        def quick_assign(fullname):
            return Atom.Atom(
                fullname.strip(), None, None, None, None, fullname, None
            ).element

        pdb_elements = {
            "H": (
                " H  ",
                " HA ",
                " HB ",
                " HD1",
                " HD2",
                " HE ",
                " HE1",
                " HE2",
                " HE3",
                " HG ",
                " HG1",
                " HH ",
                " HH2",
                " HZ ",
                " HZ2",
                " HZ3",
                "1H  ",
                "1HA ",
                "1HB ",
                "1HD ",
                "1HD1",
                "1HD2",
                "1HE ",
                "1HE2",
                "1HG ",
                "1HG1",
                "1HG2",
                "1HH1",
                "1HH2",
                "1HZ ",
                "2H  ",
                "2HA ",
                "2HB ",
                "2HD ",
                "2HD1",
                "2HD2",
                "2HE ",
                "2HE2",
                "2HG ",
                "2HG1",
                "2HG2",
                "2HH1",
                "2HH2",
                "2HZ ",
                "3H  ",
                "3HB ",
                "3HD1",
                "3HD2",
                "3HE ",
                "3HG1",
                "3HG2",
                "3HZ ",
                "HE21",
            ),
            "O": (" OH ",),  # noqa: E741
            "C": (" CH2",),
            "N": (" NH1", " NH2"),
        }

        for element, atom_names in pdb_elements.items():
            for fullname in atom_names:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", PDBConstructionWarning)
                    e = quick_assign(fullname)
                assert e == element


class SortingTests(unittest.TestCase):
    """Tests for sorting elements of the SMCRA representation."""

    def test_strict_equality(self):
        parser = PDBParser()
        structure = parser.get_structure("example", "PDB/1A8O.pdb")
        structure2 = parser.get_structure("example", "PDB/1A8O.pdb")

        assert structure.strictly_equals(structure2)
        assert structure2.strictly_equals(structure)  # Strict equality should be symmetric

        # Modify an atom
        structure2[0]["A"][(" ", 200, " ")]["CA"].name = "AC"

        assert not structure.strictly_equals(structure2)
        assert not structure2.strictly_equals(structure)  # Strict equality should be symmetric

        # Remove a chain from a model in the structure
        structure2[0].detach_child("A")

        assert not structure.strictly_equals(structure2)
        assert not structure2.strictly_equals(structure)  # Strict equality should be symmetric

        # Reset structure2
        structure2 = parser.get_structure("example", "PDB/1A8O.pdb")

        assert structure.strictly_equals(structure2)
        assert structure2.strictly_equals(structure)  # Strict equality should be symmetric

        # Change the coordinates of an atom in structure2
        structure2[0]["A"][(" ", 180, " ")]["C"].set_coord((0, 0, 0))

        assert structure.strictly_equals(structure2)
        assert structure2.strictly_equals(structure)  # Strict equality should be symmetric

        assert not structure.strictly_equals(structure2, compare_coordinates=True)
        assert not structure2.strictly_equals(structure, compare_coordinates=True)  # Strict equality should be symmetric

    def test_residue_sort(self):
        """Test atoms are sorted correctly in residues."""
        parser = PDBParser()
        structure = parser.get_structure("example", "PDB/1A8O.pdb")

        for residue in structure.get_residues():
            old = [a.name for a in residue]
            new = [a.name for a in sorted(residue)]
            special = []
            for a in ["N", "CA", "C", "O"]:
                if a in old:
                    special.append(a)
            special_len = len(special)
            assert new[0:special_len] == special, f"Sorted residue did not place N, CA, C, O first: {new}"
            assert new[special_len:] == sorted(new[special_len:]), f"After N, CA, C, O should be alphabet: {new}"

    # Tests for sorting methods
    def test_comparison_entities(self):
        """Test comparing and sorting the several SMCRA objects."""
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("example", "PDB/a_structure.pdb")

        # Test deepcopy of a structure with disordered atoms
        structure2 = deepcopy(structure)

        # Sorting (<, >, <=, <=)
        # Chains (same code as models)
        model = structure[1]
        chains = [c.id for c in sorted(model)]
        assert chains == ["A", "B", "C", " "]
        # Residues
        residues = [r.id[1] for r in sorted(structure[1]["C"])]
        assert residues == [1, 2, 3, 4, 0]
        # Atoms
        for residue in structure.get_residues():
            old = [a.name for a in residue]
            new = [a.name for a in sorted(residue)]

            special = [a for a in ("N", "CA", "C", "O") if a in old]
            len_special = len(special)
            # Placed N, CA, C, O first?
            assert new[:len_special] == special, f"Sorted residue did not place N, CA, C, O first: {new}"
            # Placed everyone else alphabetically?
            assert new[len_special:] == sorted(new[len_special:]), f"After N, CA, C, O order Should be alphabetical: {new}"
        # DisorderedResidue
        residues = [r.id[1] for r in sorted(structure[1]["A"])][79:81]
        assert residues == [80, 81]
        # Insertion code + hetflag + chain
        residues = list(structure[1]["B"]) + [structure[1]["A"][44]]
        assert [("{}" * 4).format(r.parent.id, *r.id) for r in sorted(residues)] == [
                "A 44 ",
                "B 44 ",
                "B 46 ",
                "B 47 ",
                "B 48 ",
                "B 49 ",
                "B 50 ",
                "B 51 ",
                "B 51A",
                "B 52 ",
                "BH_SEP45 ",
                "BW0 ",
            ]
        # DisorderedAtom
        atoms = [a.altloc for a in sorted(structure[1]["A"][74]["OD1"])]
        assert atoms == ["A", "B"]

        # Comparisons
        # Structure
        assert structure == structure2
        assert structure <= structure2
        assert structure >= structure2
        structure2.id = "new_id"
        assert structure != structure2
        assert structure < structure2
        assert structure <= structure2
        assert structure2 > structure
        assert structure2 >= structure

        # Model
        assert model == model  # __eq__ same type
        assert structure[0] != structure[1]

        assert structure[0] != []  # __eq__ diff. types
        assert structure != model

        # residues with same ID string should not be equal if the parent is not equal
        res1, res2, res3 = residues[0], residues[-1], structure2[1]["A"][44]
        assert res1.id == res2.id
        assert res2 == res3  # Equality of identical residues with different structure ID
        assert res1 != res2
        assert res1 > res2
        assert res1 >= res2
        assert res2 < res1
        assert res2 <= res1

        # atom should not be equal if the parent is not equal
        atom1, atom2, atom3 = res1["CA"], res2["CA"], res3["CA"]
        assert atom2 == atom3  # Equality of identical atoms with different structure ID
        assert atom1 > atom2
        assert atom1 >= atom2
        assert atom2 >= atom3
        assert atom1 != atom2
        assert atom2 < atom1
        assert atom2 <= atom1
        assert atom2 <= atom3


class IterationTests(unittest.TestCase):
    """Tests iterating over the SMCRA hierarchy."""

    @classmethod
    def setUpClass(cls):
        parser = PDBParser(PERMISSIVE=True, QUIET=True)
        cls.structure = parser.get_structure("X", "PDB/a_structure.pdb")

    def test_get_chains(self):
        """Yields chains from different models separately."""
        chains = [chain.id for chain in self.structure.get_chains()]
        assert chains == ["A", "A", "B", "C", " "]

    def test_get_residues(self):
        """Yields all residues from all models."""
        residues = [resi.id for resi in self.structure.get_residues()]
        assert len(residues) == 179

    def test_get_atoms(self):
        """Yields all atoms from the structure, excluding duplicates and ALTLOCs which are not parsed."""
        atoms = [
            "%12s" % str((atom.id, atom.altloc)) for atom in self.structure.get_atoms()
        ]
        assert len(atoms) == 835


class ChangingIdTests(unittest.TestCase):
    """Tests changing properties of SMCRA objects."""

    def setUp(self):
        parser = PDBParser(PERMISSIVE=True, QUIET=True)
        self.structure = parser.get_structure("X", "PDB/a_structure.pdb")

    def test_change_model_id(self):
        """Change the id of a model."""
        for model in self.structure:
            break  # Get first model in structure
        model.id = 2
        assert model.id == 2
        assert 2 in self.structure
        assert 0 not in self.structure

    def test_change_model_id_warns(self):
        """Warning when changing id to a value already in use by another child."""
        model = next(iter(self.structure))
        with pytest.warns(BiopythonWarning):
            model.id = 1
        # make sure children were not overwritten
        assert model.id == 1
        assert len(self.structure.child_list) == 2
        assert 1 in self.structure

    def test_change_chain_id(self):
        """Change the id of a model."""
        chain = next(iter(self.structure.get_chains()))
        chain.id = "R"
        assert chain.id == "R"
        model = next(iter(self.structure))
        assert "R" in model

    def test_change_id_to_self(self):
        """Changing the id to itself does nothing (does not raise)."""
        chain = next(iter(self.structure.get_chains()))
        chain_id = chain.id
        chain.id = chain_id
        assert chain.id == chain_id

    def test_change_residue_id(self):
        """Change the id of a residue."""
        chain = next(iter(self.structure.get_chains()))
        res = chain[("H_PCA", 1, " ")]
        res.id = (" ", 1, " ")

        assert res.id == (" ", 1, " ")
        assert (" ", 1, " ") in chain
        assert ("H_PCA", 1, " ") not in chain
        assert chain[(" ", 1, " ")] == res

    def test_full_id_is_updated_residue(self):
        """Invalidate cached full_ids if an id is changed."""
        atom = next(iter(self.structure.get_atoms()))

        # Generate the original full id.
        original_id = atom.get_full_id()
        assert original_id == ("X", 0, "A", ("H_PCA", 1, " "), ("N", " "))
        residue = next(iter(self.structure.get_residues()))

        # Make sure the full id was in fact cached,
        # so we need to invalidate it later.
        assert residue.full_id == ("X", 0, "A", ("H_PCA", 1, " "))

        # Changing the residue's id should lead to an updated full id.
        residue.id = (" ", 1, " ")
        new_id = atom.get_full_id()
        assert original_id != new_id
        assert new_id == ("X", 0, "A", (" ", 1, " "), ("N", " "))

    def test_full_id_is_updated_chain(self):
        """Invalidate cached full_ids if an id is changed."""
        atom = next(iter(self.structure.get_atoms()))

        # Generate the original full id.
        original_id = atom.get_full_id()
        assert original_id == ("X", 0, "A", ("H_PCA", 1, " "), ("N", " "))
        residue = next(iter(self.structure.get_residues()))

        # Make sure the full id was in fact cached,
        # so we need to invalidate it later.
        assert residue.full_id == ("X", 0, "A", ("H_PCA", 1, " "))
        chain = next(iter(self.structure.get_chains()))

        # Changing the chain's id should lead to an updated full id.
        chain.id = "Q"
        new_id = atom.get_full_id()
        assert original_id != new_id
        assert new_id == ("X", 0, "Q", ("H_PCA", 1, " "), ("N", " "))


class TransformTests(unittest.TestCase):
    """Tests transforming the coordinates of Atoms in a structure."""

    def setUp(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            self.s = PDBParser(PERMISSIVE=True).get_structure(
                "X", "PDB/a_structure.pdb"
            )
        self.m = self.s.get_list()[0]
        self.c = self.m.get_list()[0]
        self.r = self.c.get_list()[0]
        self.a = self.r.get_list()[0]

    def get_total_pos(self, o):
        """Sum of positions of atoms in an entity along with the number of atoms."""
        if hasattr(o, "get_coord"):
            return o.get_coord(), 1
        total_pos = np.array((0.0, 0.0, 0.0))
        total_count = 0
        for p in o.get_list():
            pos, count = self.get_total_pos(p)
            total_pos += pos
            total_count += count
        return total_pos, total_count

    def get_pos(self, o):
        """Average atom position in an entity."""
        pos, count = self.get_total_pos(o)
        return pos / count

    def test_transform(self):
        """Transform entities (rotation and translation)."""
        for o in (self.s, self.m, self.c, self.r, self.a):
            rotation = rotmat(Vector(1, 3, 5), Vector(1, 0, 0))
            translation = np.array((2.4, 0, 1), "f")
            oldpos = self.get_pos(o)
            o.transform(rotation, translation)
            newpos = self.get_pos(o)
            newpos_check = np.dot(oldpos, rotation) + translation
            for i in range(3):
                assert newpos[i] == pytest.approx(newpos_check[i], abs=5e-8)


class CopyTests(unittest.TestCase):
    """Tests copying SMCRA objects."""

    def setUp(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            self.s = PDBParser(PERMISSIVE=True).get_structure(
                "X", "PDB/a_structure.pdb"
            )
        self.m = self.s.get_list()[0]
        self.c = self.m.get_list()[0]
        self.r = self.c.get_list()[0]
        self.a = self.r.get_list()[0]

    def test_atom_copy(self):
        aa = self.a.copy()
        assert self.a is not aa
        assert self.a.get_coord() is not aa.get_coord()

    def test_entity_copy(self):
        """Make a copy of a residue."""
        for e in (self.s, self.m, self.c, self.r):
            ee = e.copy()
            assert e is not ee
            assert e.get_list()[0] is not ee.get_list()[0]


class IndexingTests(unittest.TestCase):
    """Tests for indexing SMCRA objects."""

    def setUp(self):
        parser = PDBParser(PERMISSIVE=True, QUIET=True)
        self.structure = parser.get_structure("a", "PDB/1LCD.pdb")

    def test_res_indexing(self):
        chain = self.structure[0]["A"]  # Get first chain in the first model
        all_residues = list(chain.get_residues())
        res_positions = [
            res.id[1] for res in all_residues if res.id[0] == " "
        ]  # get only standard residues, e.g. (' ', 1, ' ') vs ('W', 56, ' ')

        for res in res_positions:
            assert res in chain  # check if residue is in chain
            assert isinstance(res, int)  # check if residue is an int

            assert np.int64(res) in chain  # check if residue is recognised in chain as an `np.int64`

            with pytest.raises(KeyError):
                chain[oct(int(res))]  # the conversion should not work for octal strings


class CenterOfMassTests(unittest.TestCase):
    """Tests calculating centers of mass/geometry."""

    @classmethod
    def setUpClass(cls):
        cls.parser = parser = PDBParser()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            cls.structure = parser.get_structure("a", "PDB/1LCD.pdb")

    def test_structure_com(self):
        """Calculate Structure center of mass."""
        com = self.structure.center_of_mass()

        assert np.allclose(com, [19.870, 25.455, 28.753], atol=1e-3)

    def test_structure_cog(self):
        """Calculate Structure center of geometry."""
        cog = self.structure.center_of_mass(geometric=True)

        assert np.allclose(cog, [19.882, 25.842, 28.333], atol=1e-3)

    def test_chain_cog(self):
        """Calculate center of geometry of individual chains."""
        expected = {
            "A": [20.271, 30.191, 23.563],
            "B": [19.272, 21.163, 33.711],
            "C": [19.610, 20.599, 32.708],
        }

        for chain in self.structure[0].get_chains():  # one model only
            cog = chain.center_of_mass(geometric=True)
            assert np.allclose(cog, expected[chain.id], atol=1e-3)

    def test_com_empty_structure(self):
        """Center of mass of empty structure raises ValueError."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            s = self.parser.get_structure("b", "PDB/disordered.pdb")  # smaller

        for child in list(s):
            s.detach_child(child.id)

        with pytest.raises(ValueError):
            s.center_of_mass()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
