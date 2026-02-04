# Copyright 2001 by Brad Chapman.  All rights reserved.
# Revisions copyright 2008-2009 by Michiel de Hoon. All rights reserved.
# Revisions copyright 2010-2011 by Peter Cock. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for Primer-based programs in the Emboss suite."""

# standard library
import os
import unittest
import pytest

from Bio.Emboss import Primer3

# local stuff
from Bio.Emboss import PrimerSearch


class Primer3ParseTest(unittest.TestCase):
    def setUp(self):
        self.test_files = [
            os.path.join("Emboss", "bac_find.primer3"),
            os.path.join("Emboss", "cds_forward.primer3"),
            os.path.join("Emboss", "cds_reverse.primer3"),
            os.path.join("Emboss", "short.primer3"),
            os.path.join("Emboss", "internal_oligo.primer3"),
            os.path.join("Emboss", "no_oligo.primer3"),
        ]

    def test_simple_parse(self):
        """Make sure that we can use all single target primer3 files."""
        for file in self.test_files:
            # First using read...
            with open(file) as handle:
                Primer3.read(handle)
            # Now using parse...
            with open(file) as handle:
                assert 1 == len(list(Primer3.parse(handle)))

    def test_indepth_regular_parse(self):
        """Make sure we get the data from normal primer3 files okay."""
        regular_file = self.test_files[0]
        with open(regular_file) as handle:
            primer_info = Primer3.read(handle)

        assert len(primer_info.primers) == 5
        assert primer_info.comments == "# PRIMER3 RESULTS FOR AC074298\n"
        assert primer_info.primers[1].forward_seq == "CCGGTTTCTCTGGTTGAAAA"
        assert primer_info.primers[2].reverse_seq == "TCACATTCCCAAATGTAGATCG"
        assert primer_info.primers[0].size == 218
        assert len(primer_info.primers[0]) == 218
        assert primer_info.primers[3].forward_start == 112
        assert primer_info.primers[3].forward_length == 20
        assert primer_info.primers[3].forward_tm == 59.57
        assert primer_info.primers[3].forward_gc == 45.00
        assert primer_info.primers[4].reverse_start == 304
        assert primer_info.primers[4].reverse_length == 22
        assert primer_info.primers[4].reverse_tm == 59.61
        assert primer_info.primers[4].reverse_gc == 40.91

    def test_in_depth_single_parse(self):
        """Make sure we get info right from a single primer find."""
        file = self.test_files[1]
        with open(file) as handle:
            primer_info = Primer3.read(handle)

        assert len(primer_info.primers) == 5
        assert primer_info.comments == "# PRIMER3 RESULTS FOR 26964-28647#\n"
        assert primer_info.primers[1].reverse_seq == ""
        assert primer_info.primers[1].internal_seq == ""
        assert primer_info.primers[3].forward_seq == "TGTGATTGCTTGAGCTGGAC"
        assert primer_info.primers[3].internal_seq == ""
        assert primer_info.primers[3].forward_start == 253

    def test_internal_oligo_single_parse(self):
        """Make sure we can parse an internal oligo file correctly."""
        # these files are generated when designing hybridization probes.
        file = self.test_files[4]
        with open(file) as handle:
            primer_info = Primer3.read(handle)

        assert len(primer_info.primers) == 5
        assert primer_info.comments == "# EPRIMER3 RESULTS FOR YNL138W-A\n"
        assert primer_info.primers[0].internal_length == 22
        assert primer_info.primers[1].internal_seq == "TTGCGCTTTAGTTTGAATTGAA"
        assert primer_info.primers[2].internal_tm == 58.62
        assert primer_info.primers[3].internal_start == 16
        assert primer_info.primers[4].internal_gc == 35.00

    def test_multi_record_fwd(self):
        """Test parsing multiple primer sets (NirK forward)."""
        with open(os.path.join("Emboss", "NirK.primer3")) as handle:
            targets = list(Primer3.parse(handle))

        assert len(targets) == 16
        for target in targets:
            assert len(target.primers) == 5

        assert targets[0].primers[0].forward_seq == "GCAAACTGAAAAGCGGACTC"
        assert targets[0].primers[1].forward_seq == "GGGACGTACTTTCGCACAAT"
        assert targets[0].primers[2].forward_seq == "GTCTTATGCGTGGTGGAGGT"
        assert targets[0].primers[3].forward_seq == "GTACATCAACATCCGCAACG"
        assert targets[0].primers[4].forward_seq == "CGTACATCAACATCCGCAAC"

        assert targets[1].primers[0].forward_seq == "GGAAGTGCTTCTCGTTTTCG"
        assert targets[1].primers[1].forward_seq == "TACAGAGCGTCACGGATGAG"
        assert targets[1].primers[2].forward_seq == "TTGTCATCGTGCTCTTCGTC"
        assert targets[1].primers[3].forward_seq == "GACTCCAACCTCAGCTTTCG"
        assert targets[1].primers[4].forward_seq == "GGCACGAAGAAGGACAGAAG"

        assert targets[15].primers[0].forward_seq == "TGCTTGAAAATGACGCACTC"
        assert targets[15].primers[1].forward_seq == "CTCGCTGGCTAGGTCATAGG"
        assert targets[15].primers[2].forward_seq == "TATCGCACCAAACACGGTAA"
        assert targets[15].primers[3].forward_seq == "CGATTACCCTCACCGTCACT"
        assert targets[15].primers[4].forward_seq == "TATCGCAACCACTGAGCAAG"

    def test_multi_record_full(self):
        """Test parsing multiple primer sets (NirK full)."""
        with open(os.path.join("Emboss", "NirK_full.primer3")) as handle:
            targets = list(Primer3.parse(handle))

        assert len(targets) == 16
        for target in targets:
            assert len(target.primers) == 5

        assert targets[15].primers[0].forward_seq == "ACTCACTTCGGCTGAATGCT"
        assert targets[15].primers[1].forward_seq == "GGCGATTAGCGCTGTCTATC"
        assert targets[15].primers[2].forward_seq == "ACTCACTTCGGCTGAATGCT"
        assert targets[15].primers[3].forward_seq == "TAGGCGTATAGACCGGGTTG"
        assert targets[15].primers[4].forward_seq == "AGCAAGCTGACCACTGGTTT"

        assert targets[15].primers[0].reverse_seq == "CATTTAATCCGGATGCCAAC"
        assert targets[15].primers[1].reverse_seq == "TGGCCTTTCTCTCCTCTTCA"
        assert targets[15].primers[2].reverse_seq == "ATTTAATCCGGATGCCAACA"
        assert targets[15].primers[3].reverse_seq == "CACACATTATTGGCGGTCAC"
        assert targets[15].primers[4].reverse_seq == "TCTGAAACCACCAAGGAAGC"

        assert targets[15].primers[0].internal_seq == "CCCACCAATATTTGGCTAGC"
        assert targets[15].primers[1].internal_seq == "AATCTTCTGTGCACCTTGCC"
        assert targets[15].primers[2].internal_seq == "CCCACCAATATTTGGCTAGC"
        assert targets[15].primers[3].internal_seq == "TGAGCCTGTGTTCCACACAT"
        assert targets[15].primers[4].internal_seq == "CTATGCCCTTCTGCCACAAT"


class PrimersearchParseTest(unittest.TestCase):
    def setUp(self):
        self.test_files = [os.path.join("Emboss", "bac_find.psearch")]

    def test_simple_parse(self):
        """Make sure that we can parse all primersearch files."""
        for file in self.test_files:
            with open(file) as handle:
                PrimerSearch.read(handle)

    def test_in_depth_normal_parse(self):
        """Make sure the output from a simple primersearch file is correct."""
        file = self.test_files[0]
        with open(file) as handle:
            amp_info = PrimerSearch.read(handle)

        assert len(amp_info.amplifiers) == 1
        assert "Test" in amp_info.amplifiers
        assert len(amp_info.amplifiers["Test"]) == 1

        assert amp_info.amplifiers["Test"][0].length == 218
        assert (amp_info.amplifiers["Test"][0].hit_info == "AC074298 AC074298 \n"
            "\tTelomere associated sequence for Arabidopsis thaliana "
            "TEL1N from chromosome I, complete sequence.\n"
            "\tCCGGTTTCTCTGGTTGAAAA hits forward strand at 114 "
            "with 0 mismatches\n"
            "\tTCACATTCCCAAATGTAGATCG hits reverse strand at "
            "[114] with 0 mismatches")


class PrimerSearchInputTest(unittest.TestCase):
    """Test creating input files for primersearch."""

    def setUp(self):
        pass

    def test_primer_representation(self):
        """Make sure we can output primer information correctly."""
        p_info = PrimerSearch.InputRecord()
        p_info.add_primer_set("Test", "GATC", "CATG")
        p_info.add_primer_set("Test2", "AATA", "TTAT")

        output = str(p_info)
        assert output == "Test GATC CATG\nTest2 AATA TTAT\n"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
