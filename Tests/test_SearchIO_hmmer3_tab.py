# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO HmmerIO hmmer3-tab parser."""

import os
import unittest
import pytest

from Bio.SearchIO import parse

# test case files are in the Blast directory
TEST_DIR = "Hmmer"
FMT = "hmmer3-tab"


def get_file(filename):
    """Return the path of a test file."""
    return os.path.join(TEST_DIR, filename)


class HmmscanCases(unittest.TestCase):
    """Test parsing hmmscan output."""

    def test_31b1_hmmscan_001(self):
        """Test parsing hmmer3-tab, hmmscan 3.1b1, multiple queries (tab_31b1_hmmscan_001)."""
        tab_file = get_file("tab_31b1_hmmscan_001.out")
        qresults = list(parse(tab_file, FMT))
        assert 4 == len(qresults)

        # first qresult, first hit, first hsp
        qresult = qresults[0]
        assert 1 == len(qresult)
        assert "gi|4885477|ref|NP_005359.1|" == qresult.id
        assert "-" == qresult.accession
        hit = qresult[0]
        assert 1 == len(hit)
        assert "Globin" == hit.id
        assert "PF00042.17" == hit.accession
        assert 1e-22 == hit.evalue
        assert 80.5 == hit.bitscore
        assert 0.3 == hit.bias
        assert 1.3 == hit.domain_exp_num
        assert 1 == hit.region_num
        assert 0 == hit.cluster_num
        assert 0 == hit.overlap_num
        assert 1 == hit.env_num
        assert 1 == hit.domain_obs_num
        assert 1 == hit.domain_reported_num
        assert 1 == hit.domain_included_num
        assert "Globin" == hit.description
        hsp = hit.hsps[0]
        assert 1.6e-22 == hsp.evalue
        assert 79.8 == hsp.bitscore
        assert 0.3 == hsp.bias

        # last qresult, last hit, last hsp
        qresult = qresults[-1]
        assert 5 == len(qresult)
        assert "gi|125490392|ref|NP_038661.2|" == qresult.id
        assert "-" == qresult.accession
        hit = qresult[-1]
        assert 1 == len(hit)
        assert "DUF521" == hit.id
        assert "PF04412.8" == hit.accession
        assert 0.15 == hit.evalue
        assert 10.5 == hit.bitscore
        assert 0.1 == hit.bias
        assert 1.4 == hit.domain_exp_num
        assert 1 == hit.region_num
        assert 0 == hit.cluster_num
        assert 0 == hit.overlap_num
        assert 1 == hit.env_num
        assert 1 == hit.domain_obs_num
        assert 1 == hit.domain_reported_num
        assert 0 == hit.domain_included_num
        assert "Protein of unknown function (DUF521)" == hit.description
        hsp = hit.hsps[0]
        assert 0.28 == hsp.evalue
        assert 9.6 == hsp.bitscore
        assert 0.1 == hsp.bias

    def test_30_hmmscan_001(self):
        """Test parsing hmmer3-tab, hmmscan 3.0, multiple queries (tab_30_hmmscan_001)."""
        tab_file = get_file("tab_30_hmmscan_001.out")
        qresults = parse(tab_file, FMT)
        counter = 0

        # first qresult
        qresult = next(qresults)
        counter += 1
        assert 1 == len(qresult)
        assert "gi|4885477|ref|NP_005359.1|" == qresult.id
        assert "-" == qresult.accession
        hit = qresult[0]
        assert 1 == len(hit)
        assert "Globin" == hit.id
        assert "PF00042.17" == hit.accession
        assert 6e-21 == hit.evalue
        assert 74.6 == hit.bitscore
        assert 0.3 == hit.bias
        assert 1.3 == hit.domain_exp_num
        assert 1 == hit.region_num
        assert 0 == hit.cluster_num
        assert 0 == hit.overlap_num
        assert 1 == hit.env_num
        assert 1 == hit.domain_obs_num
        assert 1 == hit.domain_reported_num
        assert 1 == hit.domain_included_num
        assert "Globin" == hit.description
        hsp = hit.hsps[0]
        assert 9.2e-21 == hsp.evalue
        assert 74.0 == hsp.bitscore
        assert 0.2 == hsp.bias

        # second qresult
        qresult = next(qresults)
        counter += 1
        assert 2 == len(qresult)
        assert "gi|126362951:116-221" == qresult.id
        assert "-" == qresult.accession
        hit = qresult[0]
        assert 1 == len(hit)
        assert "Ig_3" == hit.id
        assert "PF13927.1" == hit.accession
        assert 1.4e-09 == hit.evalue
        assert 38.2 == hit.bitscore
        assert 0.4 == hit.bias
        assert 1.3 == hit.domain_exp_num
        assert 1 == hit.region_num
        assert 0 == hit.cluster_num
        assert 0 == hit.overlap_num
        assert 1 == hit.env_num
        assert 1 == hit.domain_obs_num
        assert 1 == hit.domain_reported_num
        assert 1 == hit.domain_included_num
        assert "Immunoglobulin domain" == hit.description
        # first hsp
        hsp = hit.hsps[0]
        assert 2.1e-09 == hsp.evalue
        assert 37.6 == hsp.bitscore
        assert 0.3 == hsp.bias
        hit = qresult[1]
        assert 1 == len(hit)
        assert "Ig_2" == hit.id
        assert "PF13895.1" == hit.accession
        assert 3.5e-05 == hit.evalue
        assert 23.7 == hit.bitscore
        assert 0.1 == hit.bias
        assert 1.1 == hit.domain_exp_num
        assert 1 == hit.region_num
        assert 0 == hit.cluster_num
        assert 0 == hit.overlap_num
        assert 1 == hit.env_num
        assert 1 == hit.domain_obs_num
        assert 1 == hit.domain_reported_num
        assert 1 == hit.domain_included_num
        assert "Immunoglobulin domain" == hit.description
        hsp = hit.hsps[0]
        assert 4.3e-05 == hsp.evalue
        assert 23.4 == hsp.bitscore
        assert 0.1 == hsp.bias

        # third qresult
        qresult = next(qresults)
        counter += 1
        assert 2 == len(qresult)
        assert "gi|22748937|ref|NP_065801.1|" == qresult.id
        assert "-" == qresult.accession
        hit = qresult[0]
        assert 1 == len(hit)
        assert "Xpo1" == hit.id
        assert "PF08389.7" == hit.accession
        assert 7.8e-34 == hit.evalue
        assert 116.6 == hit.bitscore
        assert 7.8 == hit.bias
        assert 2.8 == hit.domain_exp_num
        assert 2 == hit.region_num
        assert 0 == hit.cluster_num
        assert 0 == hit.overlap_num
        assert 2 == hit.env_num
        assert 2 == hit.domain_obs_num
        assert 2 == hit.domain_reported_num
        assert 1 == hit.domain_included_num
        assert "Exportin 1-like protein" == hit.description
        hsp = hit.hsps[0]
        assert 1.1e-33 == hsp.evalue
        assert 116.1 == hsp.bitscore
        assert 3.4 == hsp.bias
        hit = qresult[1]
        assert 1 == len(hit)
        assert "IBN_N" == hit.id
        assert "PF03810.14" == hit.accession
        assert 0.0039 == hit.evalue
        assert 16.9 == hit.bitscore
        assert 0.0 == hit.bias
        assert 2.7 == hit.domain_exp_num
        assert 2 == hit.region_num
        assert 0 == hit.cluster_num
        assert 0 == hit.overlap_num
        assert 2 == hit.env_num
        assert 2 == hit.domain_obs_num
        assert 2 == hit.domain_reported_num
        assert 1 == hit.domain_included_num
        assert "Importin-beta N-terminal domain" == hit.description
        hsp = hit.hsps[0]
        assert 0.033 == hsp.evalue
        assert 14.0 == hsp.bitscore
        assert 0.0 == hsp.bias

        # last qresult
        qresult = next(qresults)
        counter += 1
        assert 5 == len(qresult)
        assert "gi|125490392|ref|NP_038661.2|" == qresult.id
        assert "-" == qresult.accession
        # first hit
        hit = qresult[0]
        assert 1 == len(hit)
        assert "Pou" == hit.id
        assert "PF00157.12" == hit.accession
        assert 7e-37 == hit.evalue
        assert 124.8 == hit.bitscore
        assert 0.5 == hit.bias
        assert 1.5 == hit.domain_exp_num
        assert 1 == hit.region_num
        assert 0 == hit.cluster_num
        assert 0 == hit.overlap_num
        assert 1 == hit.env_num
        assert 1 == hit.domain_obs_num
        assert 1 == hit.domain_reported_num
        assert 1 == hit.domain_included_num
        assert "Pou domain - N-terminal to homeobox domain" == hit.description
        hsp = hit.hsps[0]
        assert 1.4e-36 == hsp.evalue
        assert 123.9 == hsp.bitscore
        assert 0.3 == hsp.bias
        # second hit
        hit = qresult[1]
        assert 1 == len(hit)
        assert "Homeobox" == hit.id
        assert "PF00046.24" == hit.accession
        assert 2.1e-18 == hit.evalue
        assert 65.5 == hit.bitscore
        assert 1.1 == hit.bias
        assert 1.5 == hit.domain_exp_num
        assert 1 == hit.region_num
        assert 0 == hit.cluster_num
        assert 0 == hit.overlap_num
        assert 1 == hit.env_num
        assert 1 == hit.domain_obs_num
        assert 1 == hit.domain_reported_num
        assert 1 == hit.domain_included_num
        assert "Homeobox domain" == hit.description
        hsp = hit.hsps[0]
        assert 4.1e-18 == hsp.evalue
        assert 64.6 == hsp.bitscore
        assert 0.7 == hsp.bias
        # third hit
        hit = qresult[2]
        assert 1 == len(hit)
        assert "HTH_31" == hit.id
        assert "PF13560.1" == hit.accession
        assert 0.012 == hit.evalue
        assert 15.6 == hit.bitscore
        assert 0.0 == hit.bias
        assert 2.2 == hit.domain_exp_num
        assert 2 == hit.region_num
        assert 0 == hit.cluster_num
        assert 0 == hit.overlap_num
        assert 2 == hit.env_num
        assert 2 == hit.domain_obs_num
        assert 2 == hit.domain_reported_num
        assert 0 == hit.domain_included_num
        assert "Helix-turn-helix domain" == hit.description
        hsp = hit.hsps[0]
        assert 0.16 == hsp.evalue
        assert 12.0 == hsp.bitscore
        assert 0.0 == hsp.bias
        # fourth hit
        hit = qresult[3]
        assert 1 == len(hit)
        assert "Homeobox_KN" == hit.id
        assert "PF05920.6" == hit.accession
        assert 0.039 == hit.evalue
        assert 13.5 == hit.bitscore
        assert 0.0 == hit.bias
        assert 1.6 == hit.domain_exp_num
        assert 1 == hit.region_num
        assert 0 == hit.cluster_num
        assert 0 == hit.overlap_num
        assert 1 == hit.env_num
        assert 1 == hit.domain_obs_num
        assert 1 == hit.domain_reported_num
        assert 0 == hit.domain_included_num
        assert "Homeobox KN domain" == hit.description
        hsp = hit.hsps[0]
        assert 0.095 == hsp.evalue
        assert 12.3 == hsp.bitscore
        assert 0.0 == hsp.bias
        # fifth hit
        hit = qresult[4]
        assert 1 == len(hit)
        assert "DUF521" == hit.id
        assert "PF04412.8" == hit.accession
        assert 0.14 == hit.evalue
        assert 10.5 == hit.bitscore
        assert 0.1 == hit.bias
        assert 1.4 == hit.domain_exp_num
        assert 1 == hit.region_num
        assert 0 == hit.cluster_num
        assert 0 == hit.overlap_num
        assert 1 == hit.env_num
        assert 1 == hit.domain_obs_num
        assert 1 == hit.domain_reported_num
        assert 0 == hit.domain_included_num
        assert "Protein of unknown function (DUF521)" == hit.description
        hsp = hit.hsps[0]
        assert 0.26 == hsp.evalue
        assert 9.6 == hsp.bitscore
        assert 0.1 == hsp.bias

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 4 == counter

    def test_30_hmmscan_002(self):
        """Test parsing hmmer3-tab, hmmscan 3.0, single query, no hits (tab_30_hmmscan_002)."""
        tab_file = get_file("tab_30_hmmscan_002.out")
        qresults = parse(tab_file, FMT)

        with pytest.raises(StopIteration):
            next(qresults)

    def test_30_hmmscan_003(self):
        """Test parsing hmmer3-tab, hmmscan 3.0, single query, single hit, single hsp (tab_30_hmmscan_003)."""
        tab_file = get_file("tab_30_hmmscan_003.out")
        qresults = parse(tab_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1
        assert 1 == len(qresult)
        assert "gi|4885477|ref|NP_005359.1|" == qresult.id
        assert "-" == qresult.accession
        hit = qresult[0]
        assert 1 == len(hit)
        assert "Globin" == hit.id
        assert "PF00042.17" == hit.accession
        assert 6e-21 == hit.evalue
        assert 74.6 == hit.bitscore
        assert 0.3 == hit.bias
        assert 1.3 == hit.domain_exp_num
        assert 1 == hit.region_num
        assert 0 == hit.cluster_num
        assert 0 == hit.overlap_num
        assert 1 == hit.env_num
        assert 1 == hit.domain_obs_num
        assert 1 == hit.domain_reported_num
        assert 1 == hit.domain_included_num
        assert "Globin" == hit.description
        hsp = hit.hsps[0]
        assert 9.2e-21 == hsp.evalue
        assert 74.0 == hsp.bitscore
        assert 0.2 == hsp.bias

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_30_hmmscan_004(self):
        """Test parsing hmmer3-tab, hmmscan 3.0, single query, multiple hits (tab_30_hmmscan_004)."""
        tab_file = get_file("tab_30_hmmscan_004.out")
        qresults = parse(tab_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1
        assert 2 == len(qresult)
        assert "gi|126362951:116-221" == qresult.id
        assert "-" == qresult.accession
        hit = qresult[0]
        assert 1 == len(hit)
        assert "Ig_3" == hit.id
        assert "PF13927.1" == hit.accession
        assert 1.4e-09 == hit.evalue
        assert 38.2 == hit.bitscore
        assert 0.4 == hit.bias
        assert 1.3 == hit.domain_exp_num
        assert 1 == hit.region_num
        assert 0 == hit.cluster_num
        assert 0 == hit.overlap_num
        assert 1 == hit.env_num
        assert 1 == hit.domain_obs_num
        assert 1 == hit.domain_reported_num
        assert 1 == hit.domain_included_num
        assert "Immunoglobulin domain" == hit.description
        hsp = hit.hsps[0]
        assert 2.1e-09 == hsp.evalue
        assert 37.6 == hsp.bitscore
        assert 0.3 == hsp.bias
        hit = qresult[1]
        assert 1 == len(hit)
        assert "Ig_2" == hit.id
        assert "PF13895.1" == hit.accession
        assert 3.5e-05 == hit.evalue
        assert 23.7 == hit.bitscore
        assert 0.1 == hit.bias
        assert 1.1 == hit.domain_exp_num
        assert 1 == hit.region_num
        assert 0 == hit.cluster_num
        assert 0 == hit.overlap_num
        assert 1 == hit.env_num
        assert 1 == hit.domain_obs_num
        assert 1 == hit.domain_reported_num
        assert 1 == hit.domain_included_num
        assert "Immunoglobulin domain" == hit.description
        hsp = hit.hsps[0]
        assert 4.3e-05 == hsp.evalue
        assert 23.4 == hsp.bitscore
        assert 0.1 == hsp.bias

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter


class HmmsearchCases(unittest.TestCase):
    """Tests for hmmsearch output."""

    def test_31b1_hmmsearch_001(self):
        """Test parsing hmmer3-tab, hmmsearch 3.1b1, multiple queries (tab_31b1_hmmscan_001)."""
        tab_file = get_file("tab_31b1_hmmsearch_001.out")
        qresults = list(parse(tab_file, FMT))
        assert 1 == len(qresults)

        # first qresult
        qresult = qresults[0]
        assert 4 == len(qresult)
        assert "Pkinase" == qresult.id
        assert "PF00069.17" == qresult.accession

        # first hit, first hsp
        hit = qresult[0]
        assert 1 == len(hit)
        assert "sp|Q9WUT3|KS6A2_MOUSE" == hit.id
        assert "-" == hit.accession
        assert 8.5e-147 == hit.evalue
        assert 492.3 == hit.bitscore
        assert 0.0 == hit.bias
        assert 2.1 == hit.domain_exp_num
        assert 2 == hit.region_num
        assert 0 == hit.cluster_num
        assert 0 == hit.overlap_num
        assert 2 == hit.env_num
        assert 2 == hit.domain_obs_num
        assert 2 == hit.domain_reported_num
        assert 2 == hit.domain_included_num
        assert "Ribosomal protein S6 kinase alpha-2 OS=Mus musculus GN=Rps6ka2 PE=1 SV=1" == hit.description
        hsp = hit.hsps[0]
        assert 1.2e-72 == hsp.evalue
        assert 249.3 == hsp.bitscore
        assert 0.0 == hsp.bias

        # last hit, last hsp
        hit = qresult[-1]
        assert 1 == len(hit)
        assert "sp|P18652|KS6AA_CHICK" == hit.id
        assert "-" == hit.accession
        assert 2.6e-145 == hit.evalue
        assert 487.5 == hit.bitscore
        assert 0.0 == hit.bias
        assert 2.1 == hit.domain_exp_num
        assert 2 == hit.region_num
        assert 0 == hit.cluster_num
        assert 0 == hit.overlap_num
        assert 2 == hit.env_num
        assert 2 == hit.domain_obs_num
        assert 2 == hit.domain_reported_num
        assert 2 == hit.domain_included_num
        assert "Ribosomal protein S6 kinase 2 alpha OS=Gallus gallus GN=RPS6KA PE=2 SV=1" == hit.description
        hsp = hit.hsps[-1]
        assert 7.6e-72 == hsp.evalue
        assert 246.7 == hsp.bitscore
        assert 0.0 == hsp.bias


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
