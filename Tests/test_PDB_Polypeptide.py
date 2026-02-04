# Copyright 2017 by Francesco Gastaldello. All rights reserved.
# Revisions copyright 2020 Joao Rodrigues. All rights reserved.
#
# Converted by Francesco Gastaldello from an older unit test copyright 2004
# by Thomas Hamelryck.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Unit tests for the Bio.PDB.Polypeptide module."""

import unittest
import pytest

from Bio.PDB import CaPPBuilder
from Bio.PDB import PDBParser
from Bio.PDB import PPBuilder
from Bio.Seq import Seq


class PolypeptideTests(unittest.TestCase):
    """Test Polypeptide module."""

    @classmethod
    def setUpClass(self):
        pdb1 = "PDB/1A8O.pdb"
        self.parser = PDBParser(PERMISSIVE=True)
        self.structure = self.parser.get_structure("scr", pdb1)

    def test_ppbuilder_real(self):
        """Test PPBuilder on real PDB file."""
        ppb = PPBuilder()
        pp = ppb.build_peptides(self.structure)

        assert len(pp) == 3

        # Check termini
        assert pp[0][0].get_id()[1] == 152
        assert pp[0][-1].get_id()[1] == 184
        assert pp[1][0].get_id()[1] == 186
        assert pp[1][-1].get_id()[1] == 213
        assert pp[2][0].get_id()[1] == 216
        assert pp[2][-1].get_id()[1] == 220

        # Now check sequences
        pp0_seq = pp[0].get_sequence()
        pp1_seq = pp[1].get_sequence()
        pp2_seq = pp[2].get_sequence()
        assert isinstance(pp0_seq, Seq)
        assert pp0_seq == "DIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNW"
        assert pp1_seq == "TETLLVQNANPDCKTILKALGPGATLEE"
        assert pp2_seq == "TACQG"

    def test_ppbuilder_real_nonstd(self):
        """Test PPBuilder on real PDB file allowing non-standard amino acids."""
        ppb = PPBuilder()
        pp = ppb.build_peptides(self.structure, False)

        assert len(pp) == 1

        # Check the start and end positions
        assert pp[0][0].get_id()[1] == 151
        assert pp[0][-1].get_id()[1] == 220

        # Check the sequence
        s = pp[0].get_sequence()
        assert isinstance(s, Seq)
        # Here non-standard MSE are shown as M
        assert "MDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPGATLEEMMTACQG" == s

    def test_ppbuilder_torsion(self):
        """Test phi/psi angles calculated with PPBuilder."""
        ppb = PPBuilder()
        pp = ppb.build_peptides(self.structure)

        phi_psi = pp[0].get_phi_psi_list()
        assert phi_psi[0][0] is None
        assert phi_psi[0][1] == pytest.approx(-0.46297171497725553, abs=0.0005)
        assert phi_psi[1][0] == pytest.approx(-1.0873937604007962, abs=0.0005)
        assert phi_psi[1][1] == pytest.approx(2.1337707832637109, abs=0.0005)
        assert phi_psi[2][0] == pytest.approx(-2.4052232743651878, abs=0.0005)
        assert phi_psi[2][1] == pytest.approx(2.3807316946081554, abs=0.0005)

        phi_psi = pp[1].get_phi_psi_list()
        assert phi_psi[0][0] is None
        assert phi_psi[0][1] == pytest.approx(-0.6810077089092923, abs=0.0005)
        assert phi_psi[1][0] == pytest.approx(-1.2654003477656888, abs=0.0005)
        assert phi_psi[1][1] == pytest.approx(-0.58689987042756309, abs=0.0005)
        assert phi_psi[2][0] == pytest.approx(-1.7467679151684763, abs=0.0005)
        assert phi_psi[2][1] == pytest.approx(-1.5655066256698336, abs=0.0005)

        phi_psi = pp[2].get_phi_psi_list()
        assert phi_psi[0][0] is None
        assert phi_psi[0][1] == pytest.approx(-0.73222884210889716, abs=0.0005)
        assert phi_psi[1][0] == pytest.approx(-1.1044740234566259, abs=0.0005)
        assert phi_psi[1][1] == pytest.approx(-0.69681334592782884, abs=0.0005)
        assert phi_psi[2][0] == pytest.approx(-1.8497413300164958, abs=0.0005)
        assert phi_psi[2][1] == pytest.approx(0.34762889834809058, abs=0.0005)

    def test_cappbuilder_real(self):
        """Test CaPPBuilder on real PDB file."""
        ppb = CaPPBuilder()
        pp = ppb.build_peptides(self.structure)

        pp0_seq = pp[0].get_sequence()
        pp1_seq = pp[1].get_sequence()
        pp2_seq = pp[2].get_sequence()
        assert pp0_seq == "DIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNW"
        assert pp1_seq == "TETLLVQNANPDCKTILKALGPGATLEE"
        assert pp2_seq == "TACQG"
        assert [ca.serial_number for ca in pp[0].get_ca_list()] == [
                10,
                18,
                26,
                37,
                46,
                50,
                57,
                66,
                75,
                82,
                93,
                104,
                112,
                124,
                131,
                139,
                150,
                161,
                173,
                182,
                189,
                197,
                208,
                213,
                222,
                231,
                236,
                242,
                251,
                260,
                267,
                276,
                284,
            ]

    def test_cappbuilder_real_nonstd(self):
        """Test CaPPBuilder on real PDB file allowing non-standard amino acids."""
        ppb = CaPPBuilder()
        pp = ppb.build_peptides(self.structure, False)

        assert len(pp) == 1

        # Check the start and end positions
        assert pp[0][0].get_id()[1] == 151
        assert pp[0][-1].get_id()[1] == 220

        # Check the sequence
        s = pp[0].get_sequence()
        assert isinstance(s, Seq)
        # Here non-standard MSE are shown as M
        assert "MDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPGATLEEMMTACQG" == s

    def test_cappbuilder_tau(self):
        """Test tau angles calculated with CaPPBuilder."""
        ppb = CaPPBuilder()
        pp = ppb.build_peptides(self.structure)

        taus = pp[1].get_tau_list()
        assert taus[0] == pytest.approx(0.3597907225123525, abs=0.0005)
        assert taus[1] == pytest.approx(0.43239284636769254, abs=0.0005)
        assert taus[2] == pytest.approx(0.99820157492712114, abs=0.0005)
        thetas = pp[2].get_theta_list()
        assert thetas[0] == pytest.approx(1.6610069445335354, abs=0.0005)
        assert thetas[1] == pytest.approx(1.7491703334817772, abs=0.0005)
        assert thetas[2] == pytest.approx(2.0702447422720143, abs=0.0005)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
