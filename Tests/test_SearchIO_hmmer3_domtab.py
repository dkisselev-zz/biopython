# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO HmmerIO hmmer3-domtab parsers."""

import os
import unittest
import pytest

from Bio.SearchIO import parse

# test case files are in the Blast directory
TEST_DIR = "Hmmer"


def get_file(filename):
    """Return the path of a test file."""
    return os.path.join(TEST_DIR, filename)


class HmmscanCases(unittest.TestCase):
    """Testing the hmmscan domtab output."""

    fmt = "hmmscan3-domtab"

    def test_domtab_31b1_hmmscan_001(self):
        """Parsing hmmscan-domtab, hmmscan 3.1b1, multiple queries (domtab_31b1_hmmscan_001)."""
        tab_file = get_file("domtab_31b1_hmmscan_001.out")
        qresults = list(parse(tab_file, self.fmt))
        assert 4 == len(qresults)

        # first qresult, first hit, first hsp
        qresult = qresults[0]
        assert 1 == len(qresult)
        assert "gi|4885477|ref|NP_005359.1|" == qresult.id
        assert "-" == qresult.accession
        assert 154 == qresult.seq_len
        hit = qresult[0]
        assert 1 == len(hit)
        assert "Globin" == hit.id
        assert "gi|4885477|ref|NP_005359.1|" == hit.query_id
        assert "PF00042.17" == hit.accession
        assert 110 == hit.seq_len
        assert 1e-22 == hit.evalue
        assert 80.5 == hit.bitscore
        assert 0.3 == hit.bias
        assert "Globin" == hit.description
        hsp = hit.hsps[0]
        assert "Globin" == hsp.hit_id
        assert "gi|4885477|ref|NP_005359.1|" == hsp.query_id
        assert 1 == hsp.domain_index
        assert 1.1e-26 == hsp.evalue_cond
        assert 1.6e-22 == hsp.evalue
        assert 79.8 == hsp.bitscore
        assert 0.3 == hsp.bias
        assert 0 == hsp.hit_start
        assert 109 == hsp.hit_end
        assert 6 == hsp.query_start
        assert 112 == hsp.query_end
        assert 6 == hsp.env_start
        assert 113 == hsp.env_end
        assert 0.97 == hsp.acc_avg

        # last qresult, last hit, last hsp
        qresult = qresults[-1]
        assert 5 == len(qresult)
        assert "gi|125490392|ref|NP_038661.2|" == qresult.id
        assert "-" == qresult.accession
        assert 352 == qresult.seq_len
        hit = qresult[-1]
        assert 1 == len(hit)
        assert "DUF521" == hit.id
        assert "gi|125490392|ref|NP_038661.2|" == hit.query_id
        assert "PF04412.8" == hit.accession
        assert 400 == hit.seq_len
        assert 0.15 == hit.evalue
        assert 10.5 == hit.bitscore
        assert 0.1 == hit.bias
        assert "Protein of unknown function (DUF521)" == hit.description
        hsp = hit.hsps[0]
        assert "DUF521" == hsp.hit_id
        assert "gi|125490392|ref|NP_038661.2|" == hsp.query_id
        assert 1 == hsp.domain_index
        assert 9.4e-05 == hsp.evalue_cond
        assert 0.28 == hsp.evalue
        assert 9.6 == hsp.bitscore
        assert 0.1 == hsp.bias
        assert 272 == hsp.hit_start
        assert 334 == hsp.hit_end
        assert 220 == hsp.query_start
        assert 280 == hsp.query_end
        assert 196 == hsp.env_start
        assert 294 == hsp.env_end
        assert 0.77 == hsp.acc_avg

    def test_domtab_30_hmmscan_001(self):
        """Parsing hmmscan-domtab, hmmscan 3.0, multiple queries (domtab_30_hmmscan_001)."""
        tab_file = get_file("domtab_30_hmmscan_001.out")
        qresults = parse(tab_file, self.fmt)
        counter = 0

        # first qresult
        qresult = next(qresults)
        counter += 1
        assert 1 == len(qresult)
        assert "gi|4885477|ref|NP_005359.1|" == qresult.id
        assert "-" == qresult.accession
        assert 154 == qresult.seq_len
        hit = qresult[0]
        assert 1 == len(hit)
        assert "Globin" == hit.id
        assert "gi|4885477|ref|NP_005359.1|" == hit.query_id
        assert "PF00042.17" == hit.accession
        assert 108 == hit.seq_len
        assert 6e-21 == hit.evalue
        assert 74.6 == hit.bitscore
        assert 0.3 == hit.bias
        assert "Globin" == hit.description
        hsp = hit.hsps[0]
        assert "Globin" == hsp.hit_id
        assert "gi|4885477|ref|NP_005359.1|" == hsp.query_id
        assert 1 == hsp.domain_index
        assert 6.7e-25 == hsp.evalue_cond
        assert 9.2e-21 == hsp.evalue
        assert 74.0 == hsp.bitscore
        assert 0.2 == hsp.bias
        assert 0 == hsp.hit_start
        assert 107 == hsp.hit_end
        assert 6 == hsp.query_start
        assert 112 == hsp.query_end
        assert 6 == hsp.env_start
        assert 113 == hsp.env_end
        assert 0.97 == hsp.acc_avg

        # second qresult
        qresult = next(qresults)
        counter += 1
        assert 2 == len(qresult)
        assert "gi|126362951:116-221" == qresult.id
        assert "-" == qresult.accession
        assert 106 == qresult.seq_len
        hit = qresult[0]
        assert 1 == len(hit)
        assert "Ig_3" == hit.id
        assert "gi|126362951:116-221" == hit.query_id
        assert "PF13927.1" == hit.accession
        assert 75 == hit.seq_len
        assert 1.4e-09 == hit.evalue
        assert 38.2 == hit.bitscore
        assert 0.4 == hit.bias
        assert "Immunoglobulin domain" == hit.description
        hsp = hit.hsps[0]
        assert "Ig_3" == hsp.hit_id
        assert "gi|126362951:116-221" == hsp.query_id
        assert 1 == hsp.domain_index
        assert 3e-13 == hsp.evalue_cond
        assert 2.1e-09 == hsp.evalue
        assert 37.6 == hsp.bitscore
        assert 0.3 == hsp.bias
        assert 0 == hsp.hit_start
        assert 73 == hsp.hit_end
        assert 8 == hsp.query_start
        assert 84 == hsp.query_end
        assert 8 == hsp.env_start
        assert 88 == hsp.env_end
        assert 0.94 == hsp.acc_avg
        hit = qresult[1]
        assert 1 == len(hit)
        assert "Ig_2" == hit.id
        assert "gi|126362951:116-221" == hit.query_id
        assert "PF13895.1" == hit.accession
        assert 80 == hit.seq_len
        assert 3.5e-05 == hit.evalue
        assert 23.7 == hit.bitscore
        assert 0.1 == hit.bias
        assert "Immunoglobulin domain" == hit.description
        hsp = hit.hsps[0]
        assert "Ig_2" == hsp.hit_id
        assert "gi|126362951:116-221" == hsp.query_id
        assert 1 == hsp.domain_index
        assert 6.2e-09 == hsp.evalue_cond
        assert 4.3e-05 == hsp.evalue
        assert 23.4 == hsp.bitscore
        assert 0.1 == hsp.bias
        assert 0 == hsp.hit_start
        assert 80 == hsp.hit_end
        assert 8 == hsp.query_start
        assert 104 == hsp.query_end
        assert 8 == hsp.env_start
        assert 104 == hsp.env_end
        assert 0.71 == hsp.acc_avg

        # third qresult
        qresult = next(qresults)
        counter += 1
        assert 2 == len(qresult)
        assert "gi|22748937|ref|NP_065801.1|" == qresult.id
        assert "-" == qresult.accession
        assert 1204 == qresult.seq_len
        hit = qresult[0]
        assert 2 == len(hit)
        assert "Xpo1" == hit.id
        assert "gi|22748937|ref|NP_065801.1|" == hit.query_id
        assert "PF08389.7" == hit.accession
        assert 148 == hit.seq_len
        assert 7.8e-34 == hit.evalue
        assert 116.6 == hit.bitscore
        assert 7.8 == hit.bias
        assert "Exportin 1-like protein" == hit.description
        hsp = hit.hsps[0]
        assert "Xpo1" == hsp.hit_id
        assert "gi|22748937|ref|NP_065801.1|" == hsp.query_id
        assert 1 == hsp.domain_index
        assert 1.6e-37 == hsp.evalue_cond
        assert 1.1e-33 == hsp.evalue
        assert 116.1 == hsp.bitscore
        assert 3.4 == hsp.bias
        assert 1 == hsp.hit_start
        assert 148 == hsp.hit_end
        assert 109 == hsp.query_start
        assert 271 == hsp.query_end
        assert 108 == hsp.env_start
        assert 271 == hsp.env_end
        assert 0.98 == hsp.acc_avg
        hsp = hit.hsps[1]
        assert "Xpo1" == hsp.hit_id
        assert "gi|22748937|ref|NP_065801.1|" == hsp.query_id
        assert 2 == hsp.domain_index
        assert 0.35 == hsp.evalue_cond
        assert 2.4e03 == hsp.evalue
        assert -1.8 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 111 == hsp.hit_start
        assert 139 == hsp.hit_end
        assert 498 == hsp.query_start
        assert 525 == hsp.query_end
        assert 495 == hsp.env_start
        assert 529 == hsp.env_end
        assert 0.86 == hsp.acc_avg
        # next hit in the third qresult
        hit = qresult[1]
        assert 2 == len(hit)
        assert "IBN_N" == hit.id
        assert "gi|22748937|ref|NP_065801.1|" == hit.query_id
        assert "PF03810.14" == hit.accession
        assert 77 == hit.seq_len
        assert 0.0039 == hit.evalue
        assert 16.9 == hit.bitscore
        assert 0.0 == hit.bias
        assert "Importin-beta N-terminal domain" == hit.description
        hsp = hit.hsps[0]
        assert "IBN_N" == hsp.hit_id
        assert "gi|22748937|ref|NP_065801.1|" == hsp.query_id
        assert 1 == hsp.domain_index
        assert 4.8e-06 == hsp.evalue_cond
        assert 0.033 == hsp.evalue
        assert 14.0 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 3 == hsp.hit_start
        assert 75 == hsp.hit_end
        assert 35 == hsp.query_start
        assert 98 == hsp.query_end
        assert 32 == hsp.env_start
        assert 100 == hsp.env_end
        assert 0.87 == hsp.acc_avg
        hsp = hit.hsps[1]
        assert "IBN_N" == hsp.hit_id
        assert "gi|22748937|ref|NP_065801.1|" == hsp.query_id
        assert 2 == hsp.domain_index
        assert 1.2 == hsp.evalue_cond
        assert 8e03 == hsp.evalue
        assert -3.3 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 56 == hsp.hit_start
        assert 75 == hsp.hit_end
        assert 167 == hsp.query_start
        assert 186 == hsp.query_end
        assert 164 == hsp.env_start
        assert 187 == hsp.env_end
        assert 0.85 == hsp.acc_avg

        # fourth qresult
        qresult = next(qresults)
        counter += 1
        assert 5 == len(qresult)
        assert "gi|125490392|ref|NP_038661.2|" == qresult.id
        assert "-" == qresult.accession
        assert 352 == qresult.seq_len
        # first hit
        hit = qresult[0]
        assert 1 == len(hit)
        assert "Pou" == hit.id
        assert "gi|125490392|ref|NP_038661.2|" == hit.query_id
        assert "PF00157.12" == hit.accession
        assert 75 == hit.seq_len
        assert 7e-37 == hit.evalue
        assert 124.8 == hit.bitscore
        assert 0.5 == hit.bias
        assert "Pou domain - N-terminal to homeobox domain" == hit.description
        hsp = hit.hsps[0]
        assert "Pou" == hsp.hit_id
        assert "gi|125490392|ref|NP_038661.2|" == hsp.query_id
        assert 1 == hsp.domain_index
        assert 5e-40 == hsp.evalue_cond
        assert 1.4e-36 == hsp.evalue
        assert 123.9 == hsp.bitscore
        assert 0.3 == hsp.bias
        assert 2 == hsp.hit_start
        assert 75 == hsp.hit_end
        assert 132 == hsp.query_start
        assert 205 == hsp.query_end
        assert 130 == hsp.env_start
        assert 205 == hsp.env_end
        assert 0.97 == hsp.acc_avg
        # second hit
        hit = qresult[1]
        assert 1 == len(hit)
        assert "Homeobox" == hit.id
        assert "gi|125490392|ref|NP_038661.2|" == hit.query_id
        assert "PF00046.24" == hit.accession
        assert 57 == hit.seq_len
        assert 2.1e-18 == hit.evalue
        assert 65.5 == hit.bitscore
        assert 1.1 == hit.bias
        assert "Homeobox domain" == hit.description
        hsp = hit.hsps[0]
        assert "Homeobox" == hsp.hit_id
        assert "gi|125490392|ref|NP_038661.2|" == hsp.query_id
        assert 1 == hsp.domain_index
        assert 1.5e-21 == hsp.evalue_cond
        assert 4.1e-18 == hsp.evalue
        assert 64.6 == hsp.bitscore
        assert 0.7 == hsp.bias
        assert 0 == hsp.hit_start
        assert 57 == hsp.hit_end
        assert 223 == hsp.query_start
        assert 280 == hsp.query_end
        assert 223 == hsp.env_start
        assert 280 == hsp.env_end
        assert 0.98 == hsp.acc_avg
        # third hit
        hit = qresult[2]
        assert 2 == len(hit)
        assert "HTH_31" == hit.id
        assert "gi|125490392|ref|NP_038661.2|" == hit.query_id
        assert "PF13560.1" == hit.accession
        assert 64 == hit.seq_len
        assert 0.012 == hit.evalue
        assert 15.6 == hit.bitscore
        assert 0.0 == hit.bias
        assert "Helix-turn-helix domain" == hit.description
        hsp = hit.hsps[0]
        assert "HTH_31" == hsp.hit_id
        assert "gi|125490392|ref|NP_038661.2|" == hsp.query_id
        assert 1 == hsp.domain_index
        assert 5.7e-05 == hsp.evalue_cond
        assert 0.16 == hsp.evalue
        assert 12.0 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 0 == hsp.hit_start
        assert 35 == hsp.hit_end
        assert 140 == hsp.query_start
        assert 181 == hsp.query_end
        assert 140 == hsp.env_start
        assert 184 == hsp.env_end
        assert 0.96 == hsp.acc_avg
        hsp = hit.hsps[1]
        assert "HTH_31" == hsp.hit_id
        assert "gi|125490392|ref|NP_038661.2|" == hsp.query_id
        assert 2 == hsp.domain_index
        assert 0.19 == hsp.evalue_cond
        assert 5.2e02 == hsp.evalue
        assert 0.8 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 38 == hsp.hit_start
        assert 62 == hsp.hit_end
        assert 244 == hsp.query_start
        assert 268 == hsp.query_end
        assert 242 == hsp.env_start
        assert 270 == hsp.env_end
        assert 0.86 == hsp.acc_avg
        # fourth hit
        hit = qresult[3]
        assert 1 == len(hit)
        assert "Homeobox_KN" == hit.id
        assert "gi|125490392|ref|NP_038661.2|" == hit.query_id
        assert "PF05920.6" == hit.accession
        assert 40 == hit.seq_len
        assert 0.039 == hit.evalue
        assert 13.5 == hit.bitscore
        assert 0.0 == hit.bias
        assert "Homeobox KN domain" == hit.description
        hsp = hit.hsps[0]
        assert "Homeobox_KN" == hsp.hit_id
        assert "gi|125490392|ref|NP_038661.2|" == hsp.query_id
        assert 1 == hsp.domain_index
        assert 3.5e-05 == hsp.evalue_cond
        assert 0.095 == hsp.evalue
        assert 12.3 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 6 == hsp.hit_start
        assert 39 == hsp.hit_end
        assert 243 == hsp.query_start
        assert 276 == hsp.query_end
        assert 240 == hsp.env_start
        assert 277 == hsp.env_end
        assert 0.91 == hsp.acc_avg
        # fifth hit
        hit = qresult[4]
        assert 1 == len(hit)
        assert "DUF521" == hit.id
        assert "gi|125490392|ref|NP_038661.2|" == hit.query_id
        assert "PF04412.8" == hit.accession
        assert 400 == hit.seq_len
        assert 0.14 == hit.evalue
        assert 10.5 == hit.bitscore
        assert 0.1 == hit.bias
        assert "Protein of unknown function (DUF521)" == hit.description
        hsp = hit.hsps[0]
        assert "DUF521" == hsp.hit_id
        assert "gi|125490392|ref|NP_038661.2|" == hsp.query_id
        assert 1 == hsp.domain_index
        assert 9.4e-05 == hsp.evalue_cond
        assert 0.26 == hsp.evalue
        assert 9.6 == hsp.bitscore
        assert 0.1 == hsp.bias
        assert 272 == hsp.hit_start
        assert 334 == hsp.hit_end
        assert 220 == hsp.query_start
        assert 280 == hsp.query_end
        assert 196 == hsp.env_start
        assert 294 == hsp.env_end
        assert 0.77 == hsp.acc_avg

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 4 == counter

    def test_domtab_30_hmmscan_002(self):
        """Parsing hmmscan-domtab, hmmscan 3.0, single query, no hits (domtab_30_hmmscan_002)."""
        tab_file = get_file("domtab_30_hmmscan_002.out")
        qresults = parse(tab_file, self.fmt)

        with pytest.raises(StopIteration):
            next(qresults)

    def test_domtab_30_hmmscan_003(self):
        """Parsing hmmscan-domtab, hmmscan 3.0, multiple queries (domtab_30_hmmscan_003)."""
        tab_file = get_file("domtab_30_hmmscan_003.out")
        qresults = parse(tab_file, self.fmt)
        counter = 0

        qresult = next(qresults)
        counter += 1
        assert 1 == len(qresult)
        assert "gi|4885477|ref|NP_005359.1|" == qresult.id
        assert "-" == qresult.accession
        assert 154 == qresult.seq_len
        hit = qresult[0]
        assert 1 == len(hit)
        assert "Globin" == hit.id
        assert "gi|4885477|ref|NP_005359.1|" == hit.query_id
        assert "PF00042.17" == hit.accession
        assert 108 == hit.seq_len
        assert 6e-21 == hit.evalue
        assert 74.6 == hit.bitscore
        assert 0.3 == hit.bias
        assert "Globin" == hit.description
        hsp = hit.hsps[0]
        assert "Globin" == hsp.hit_id
        assert "gi|4885477|ref|NP_005359.1|" == hsp.query_id
        assert 1 == hsp.domain_index
        assert 6.7e-25 == hsp.evalue_cond
        assert 9.2e-21 == hsp.evalue
        assert 74.0 == hsp.bitscore
        assert 0.2 == hsp.bias
        assert 0 == hsp.hit_start
        assert 107 == hsp.hit_end
        assert 6 == hsp.query_start
        assert 112 == hsp.query_end
        assert 6 == hsp.env_start
        assert 113 == hsp.env_end
        assert 0.97 == hsp.acc_avg

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_domtab_30_hmmscan_004(self):
        """Parsing hmmscan-domtab, hmmscan 3.0, multiple queries (domtab_30_hmmscan_004)."""
        tab_file = get_file("domtab_30_hmmscan_004.out")
        qresults = parse(tab_file, self.fmt)
        counter = 0

        qresult = next(qresults)
        counter += 1
        assert 2 == len(qresult)
        assert "gi|126362951:116-221" == qresult.id
        assert "-" == qresult.accession
        assert 106 == qresult.seq_len
        hit = qresult[0]
        assert 1 == len(hit)
        assert "Ig_3" == hit.id
        assert "gi|126362951:116-221" == hit.query_id
        assert "PF13927.1" == hit.accession
        assert 75 == hit.seq_len
        assert 1.4e-09 == hit.evalue
        assert 38.2 == hit.bitscore
        assert 0.4 == hit.bias
        assert "Immunoglobulin domain" == hit.description
        hsp = hit.hsps[0]
        assert "Ig_3" == hsp.hit_id
        assert "gi|126362951:116-221" == hsp.query_id
        assert 1 == hsp.domain_index
        assert 3e-13 == hsp.evalue_cond
        assert 2.1e-09 == hsp.evalue
        assert 37.6 == hsp.bitscore
        assert 0.3 == hsp.bias
        assert 0 == hsp.hit_start
        assert 73 == hsp.hit_end
        assert 8 == hsp.query_start
        assert 84 == hsp.query_end
        assert 8 == hsp.env_start
        assert 88 == hsp.env_end
        assert 0.94 == hsp.acc_avg
        hit = qresult[1]
        assert 1 == len(hit)
        assert "Ig_2" == hit.id
        assert "gi|126362951:116-221" == hit.query_id
        assert "PF13895.1" == hit.accession
        assert 80 == hit.seq_len
        assert 3.5e-05 == hit.evalue
        assert 23.7 == hit.bitscore
        assert 0.1 == hit.bias
        assert "Immunoglobulin domain" == hit.description
        hsp = hit.hsps[0]
        assert "Ig_2" == hsp.hit_id
        assert "gi|126362951:116-221" == hsp.query_id
        assert 1 == hsp.domain_index
        assert 6.2e-09 == hsp.evalue_cond
        assert 4.3e-05 == hsp.evalue
        assert 23.4 == hsp.bitscore
        assert 0.1 == hsp.bias
        assert 0 == hsp.hit_start
        assert 80 == hsp.hit_end
        assert 8 == hsp.query_start
        assert 104 == hsp.query_end
        assert 8 == hsp.env_start
        assert 104 == hsp.env_end
        assert 0.71 == hsp.acc_avg

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter


class HmmersearchCases(unittest.TestCase):
    """Test the hmmsearch domtab parser."""

    fmt = "hmmsearch3-domtab"

    def test_domtab_31b1_hmmsearch_001(self):
        """Parsing hmmsearch-domtab, hmmsearch 3.1b1, single query (domtab_31b1_hmmsearch_001)."""
        tab_file = get_file("domtab_31b1_hmmsearch_001.out")
        qresults = list(parse(tab_file, self.fmt))

        assert 1 == len(qresults)

        qresult = qresults[0]
        assert "Pkinase" == qresult.id
        assert "PF00069.17" == qresult.accession
        assert 260 == qresult.seq_len
        hit = qresult[0]
        assert 2 == len(hit)
        assert "sp|Q9WUT3|KS6A2_MOUSE" == hit.id
        assert "Pkinase" == hit.query_id
        assert "-" == hit.accession
        assert 733 == hit.seq_len
        assert 8.5e-147 == hit.evalue
        assert 492.3 == hit.bitscore
        assert 0.0 == hit.bias
        assert "Ribosomal protein S6 kinase alpha-2 OS=Mus musculus GN=Rps6ka2 PE=1 SV=1" == hit.description
        hsp = hit.hsps[0]
        assert "sp|Q9WUT3|KS6A2_MOUSE" == hsp.hit_id
        assert "Pkinase" == hsp.query_id
        assert 1 == hsp.domain_index
        assert 2.6e-75 == hsp.evalue_cond
        assert 3.6e-70 == hsp.evalue
        assert 241.2 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 58 == hsp.hit_start
        assert 318 == hsp.hit_end
        assert 0 == hsp.query_start
        assert 260 == hsp.query_end
        assert 58 == hsp.env_start
        assert 318 == hsp.env_end
        assert 0.95 == hsp.acc_avg

    def test_domtab_30_hmmsearch_001(self):
        """Parsing hmmsearch-domtab, hmmsearch 3.0, multiple queries (domtab_30_hmmsearch_001)."""
        tab_file = get_file("domtab_30_hmmsearch_001.out")
        qresults = parse(tab_file, self.fmt)

        # first qresult
        # we only want to check the coordinate switch actually
        # so checking the first hsp of the first hit of the qresult is enough
        qresult = next(qresults)
        assert 7 == len(qresult)
        assert "Pkinase" == qresult.id
        assert "PF00069.17" == qresult.accession
        assert 260 == qresult.seq_len
        hit = qresult[0]
        assert 2 == len(hit)
        assert "sp|Q9WUT3|KS6A2_MOUSE" == hit.id
        assert "Pkinase" == hit.query_id
        assert "-" == hit.accession
        assert 733 == hit.seq_len
        assert 8.4e-147 == hit.evalue
        assert 492.3 == hit.bitscore
        assert 0.0 == hit.bias
        assert "Ribosomal protein S6 kinase alpha-2 OS=Mus musculus GN=Rps6ka2 PE=2 SV=1" == hit.description
        hsp = hit.hsps[0]
        assert "sp|Q9WUT3|KS6A2_MOUSE" == hsp.hit_id
        assert "Pkinase" == hsp.query_id
        assert 1 == hsp.domain_index
        assert 4.6e-75 == hsp.evalue_cond
        assert 3.5e-70 == hsp.evalue
        assert 241.2 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 58 == hsp.hit_start
        assert 318 == hsp.hit_end
        assert 0 == hsp.query_start
        assert 260 == hsp.query_end
        assert 58 == hsp.env_start
        assert 318 == hsp.env_end
        assert 0.95 == hsp.acc_avg


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
