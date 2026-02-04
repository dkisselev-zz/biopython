# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO BlatIO parsers."""

import os
import unittest
import pytest

from Bio.SearchIO import parse

# test case files are in the Blast directory
TEST_DIR = "Blat"
FMT = "blat-psl"


def get_file(filename):
    """Return the path of a test file."""
    return os.path.join(TEST_DIR, filename)


class BlatPslCases(unittest.TestCase):
    def test_psl_34_001(self, testf="psl_34_001.psl", pslx=False):
        """Test parsing blat output (psl_34_001.psl)."""
        blat_file = get_file(testf)
        self.qresults = list(parse(blat_file, FMT, pslx=pslx))
        assert 2 == len(self.qresults)
        # check common attributes
        for qresult in self.qresults:
            for hit in qresult:
                assert qresult.id == hit.query_id
                for hsp in hit:
                    assert hit.id == hsp.hit_id
                    assert qresult.id == hsp.query_id

        # test first qresult
        qresult = self.qresults[0]
        assert "hg18_dna" == qresult.id
        assert "blat" == qresult.program
        assert 33 == qresult.seq_len
        assert 3 == len(qresult)
        # first qresult, first hit
        hit = qresult[0]
        assert "chr4" == hit.id
        assert 191154276 == hit.seq_len
        assert 1 == len(hit.hsps)
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 16 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 0 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert 1 == hsp[0].query_strand
        assert 11 == hsp.query_start
        assert 61646095 == hsp.hit_start
        assert 27 == hsp.query_end
        assert 61646111 == hsp.hit_end
        assert 1 == len(hsp)
        assert [16] == hsp.query_span_all
        assert [16] == hsp.hit_span_all
        assert [(11, 27)] == hsp.query_range_all
        assert [(61646095, 61646111)] == hsp.hit_range_all
        # first qresult, second hit
        hit = qresult[1]
        assert "chr1" == hit.id
        assert 249250621 == hit.seq_len
        assert 1 == len(hit.hsps)
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 33 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 0 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert 1 == hsp[0].query_strand
        assert 0 == hsp.query_start
        assert 10271783 == hsp.hit_start
        assert 33 == hsp.query_end
        assert 10271816 == hsp.hit_end
        assert 1 == len(hsp)
        assert [33] == hsp.query_span_all
        assert [33] == hsp.hit_span_all
        assert [(0, 33)] == hsp.query_range_all
        assert [(10271783, 10271816)] == hsp.hit_range_all
        # first qresult, third hit
        hit = qresult[2]
        assert "chr2" == hit.id
        assert 243199373 == hit.seq_len
        assert 1 == len(hit.hsps)
        # first qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        assert 17 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 0 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert -1 == hsp[0].query_strand
        assert 8 == hsp.query_start
        assert 53575980 == hsp.hit_start
        assert 25 == hsp.query_end
        assert 53575997 == hsp.hit_end
        assert 1 == len(hsp)
        assert [17] == hsp.query_span_all
        assert [17] == hsp.hit_span_all
        assert [(8, 25)] == hsp.query_range_all
        assert [(53575980, 53575997)] == hsp.hit_range_all

        # test second qresult
        qresult = self.qresults[1]
        assert "hg19_dna" == qresult.id
        assert "blat" == qresult.program
        assert 50 == qresult.seq_len
        assert 10 == len(qresult)
        # second qresult, first hit
        hit = qresult[0]
        assert "chr9" == hit.id
        assert 141213431 == hit.seq_len
        assert 1 == len(hit.hsps)
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 38 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 3 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert 1 == hsp[0].query_strand
        assert 9 == hsp.query_start
        assert 85737865 == hsp.hit_start
        assert 50 == hsp.query_end
        assert 85737906 == hsp.hit_end
        assert 1 == len(hsp)
        assert [41] == hsp.query_span_all
        assert [41] == hsp.hit_span_all
        assert [(9, 50)] == hsp.query_range_all
        assert [(85737865, 85737906)] == hsp.hit_range_all
        # second qresult, second hit
        hit = qresult[1]
        assert "chr8" == hit.id
        assert 146364022 == hit.seq_len
        assert 1 == len(hit.hsps)
        # second qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 41 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 0 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert 1 == hsp[0].query_strand
        assert 8 == hsp.query_start
        assert 95160479 == hsp.hit_start
        assert 49 == hsp.query_end
        assert 95160520 == hsp.hit_end
        assert 1 == len(hsp)
        assert [41] == hsp.query_span_all
        assert [41] == hsp.hit_span_all
        assert [(8, 49)] == hsp.query_range_all
        assert [(95160479, 95160520)] == hsp.hit_range_all
        # second qresult, third hit
        hit = qresult[2]
        assert "chr22" == hit.id
        assert 51304566 == hit.seq_len
        assert 2 == len(hit.hsps)
        # second qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        assert 33 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 3 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert 1 == hsp[0].query_strand
        assert 11 == hsp.query_start
        assert 42144400 == hsp.hit_start
        assert 47 == hsp.query_end
        assert 42144436 == hsp.hit_end
        assert 1 == len(hsp)
        assert [36] == hsp.query_span_all
        assert [36] == hsp.hit_span_all
        assert [(11, 47)] == hsp.query_range_all
        assert [(42144400, 42144436)] == hsp.hit_range_all
        # second qresult, third hit, second hsp
        hsp = qresult[2].hsps[1]
        assert 35 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 2 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert -1 == hsp[0].query_strand
        assert 12 == hsp.query_start
        assert 48997405 == hsp.hit_start
        assert 49 == hsp.query_end
        assert 48997442 == hsp.hit_end
        assert 1 == len(hsp)
        assert [37] == hsp.query_span_all
        assert [37] == hsp.hit_span_all
        assert [(12, 49)] == hsp.query_range_all
        assert [(48997405, 48997442)] == hsp.hit_range_all
        # second qresult, fourth hit
        hit = qresult[3]
        assert "chr2" == hit.id
        assert 243199373 == hit.seq_len
        assert 2 == len(hit.hsps)
        # second qresult, fourth hit, first hsp
        hsp = qresult[3].hsps[0]
        assert 43 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 1 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 1 == hsp.query_gapopen_num
        assert 4 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert 1 == hsp[0].query_strand
        assert 1 == hsp.query_start
        assert 183925984 == hsp.hit_start
        assert 49 == hsp.query_end
        assert 183926028 == hsp.hit_end
        assert 2 == len(hsp)
        assert [6, 38] == hsp.query_span_all
        assert [6, 38] == hsp.hit_span_all
        assert [(1, 7), (11, 49)] == hsp.query_range_all
        assert [1, 11] == hsp.query_start_all
        assert [7, 49] == hsp.query_end_all
        assert [(183925984, 183925990), (183925990, 183926028)] == hsp.hit_range_all
        assert [183925984, 183925990] == hsp.hit_start_all
        assert [183925990, 183926028] == hsp.hit_end_all
        # second qresult, fourth hit, second hsp
        hsp = qresult[3].hsps[1]
        assert 35 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 1 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert -1 == hsp[0].query_strand
        assert 13 == hsp.query_start
        assert 120641740 == hsp.hit_start
        assert 49 == hsp.query_end
        assert 120641776 == hsp.hit_end
        assert 1 == len(hsp)
        assert [36] == hsp.query_span_all
        assert [36] == hsp.hit_span_all
        assert [(13, 49)] == hsp.query_range_all
        assert [(120641740, 120641776)] == hsp.hit_range_all
        # second qresult, fifth hit
        hit = qresult[4]
        assert "chr19" == hit.id
        assert 59128983 == hit.seq_len
        assert 3 == len(hit.hsps)
        # second qresult, fifth hit, first hsp
        hsp = qresult[4].hsps[0]
        assert 34 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 2 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 1 == hsp.hit_gapopen_num
        assert 134 == hsp.hit_gap_num
        assert 1 == hsp[0].query_strand
        assert 10 == hsp.query_start
        assert 35483340 == hsp.hit_start
        assert 46 == hsp.query_end
        assert 35483510 == hsp.hit_end
        assert 2 == len(hsp)
        assert [25, 11] == hsp.query_span_all
        assert [25, 11] == hsp.hit_span_all
        assert [(10, 35), (35, 46)] == hsp.query_range_all
        assert [(35483340, 35483365), (35483499, 35483510)] == hsp.hit_range_all
        # second qresult, fifth hit, second hsp
        hsp = qresult[4].hsps[1]
        assert 39 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 0 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert -1 == hsp[0].query_strand
        assert 10 == hsp.query_start
        assert 54017130 == hsp.hit_start
        assert 49 == hsp.query_end
        assert 54017169 == hsp.hit_end
        assert 1 == len(hsp)
        assert [39] == hsp.query_span_all
        assert [39] == hsp.hit_span_all
        assert [(10, 49)] == hsp.query_range_all
        assert [(54017130, 54017169)] == hsp.hit_range_all
        # second qresult, fifth hit, third hsp
        hsp = qresult[4].hsps[2]
        assert 36 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 3 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert -1 == hsp[0].query_strand
        assert 10 == hsp.query_start
        assert 553742 == hsp.hit_start
        assert 49 == hsp.query_end
        assert 553781 == hsp.hit_end
        assert 1 == len(hsp)
        assert [39] == hsp.query_span_all
        assert [39] == hsp.hit_span_all
        assert [(10, 49)] == hsp.query_range_all
        assert [(553742, 553781)] == hsp.hit_range_all

    def test_psl_34_002(self, testf="psl_34_002.psl", pslx=False):
        """Test parsing blat output (psl_34_001.psl)."""
        blat_file = get_file(testf)
        self.qresults = list(parse(blat_file, FMT, pslx=pslx))
        assert 0 == len(self.qresults)

    def test_psl_34_003(self, testf="psl_34_003.psl", pslx=False):
        """Test parsing blat output (psl_34_003.psl)."""
        blat_file = get_file(testf)
        self.qresults = list(parse(blat_file, FMT, pslx=pslx))
        assert 1 == len(self.qresults)
        # check common attributes
        for qresult in self.qresults:
            for hit in qresult:
                assert qresult.id == hit.query_id
                for hsp in hit:
                    assert hit.id == hsp.hit_id
                    assert qresult.id == hsp.query_id

        # test first qresult
        qresult = self.qresults[0]
        assert "hg18_dna" == qresult.id
        assert "blat" == qresult.program
        assert 33 == qresult.seq_len
        assert 3 == len(qresult)
        # first qresult, first hit
        hit = qresult[0]
        assert "chr4" == hit.id
        assert 191154276 == hit.seq_len
        assert 1 == len(hit.hsps)
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 16 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 0 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert 1 == hsp[0].query_strand
        assert 11 == hsp.query_start
        assert 61646095 == hsp.hit_start
        assert 27 == hsp.query_end
        assert 61646111 == hsp.hit_end
        assert 1 == len(hsp)
        assert [16] == hsp.query_span_all
        assert [16] == hsp.hit_span_all
        assert [(11, 27)] == hsp.query_range_all
        assert [(61646095, 61646111)] == hsp.hit_range_all
        # first qresult, second hit
        hit = qresult[1]
        assert "chr1" == hit.id
        assert 249250621 == hit.seq_len
        assert 1 == len(hit.hsps)
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 33 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 0 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert 1 == hsp[0].query_strand
        assert 0 == hsp.query_start
        assert 10271783 == hsp.hit_start
        assert 33 == hsp.query_end
        assert 10271816 == hsp.hit_end
        assert 1 == len(hsp)
        assert [33] == hsp.query_span_all
        assert [33] == hsp.hit_span_all
        assert [(0, 33)] == hsp.query_range_all
        assert [(10271783, 10271816)] == hsp.hit_range_all
        # first qresult, third hit
        hit = qresult[2]
        assert "chr2" == hit.id
        assert 243199373 == hit.seq_len
        assert 1 == len(hit.hsps)
        # first qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        assert 17 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 0 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert -1 == hsp[0].query_strand
        assert 8 == hsp.query_start
        assert 53575980 == hsp.hit_start
        assert 25 == hsp.query_end
        assert 53575997 == hsp.hit_end
        assert 1 == len(hsp)
        assert [17] == hsp.query_span_all
        assert [17] == hsp.hit_span_all
        assert [(8, 25)] == hsp.query_range_all
        assert [(53575980, 53575997)] == hsp.hit_range_all

    def test_psl_34_004(self, testf="psl_34_004.psl", pslx=False):
        """Test parsing blat output (psl_34_004.psl)."""
        blat_file = get_file(testf)
        self.qresults = list(parse(blat_file, FMT, pslx=pslx))
        assert 1 == len(self.qresults)
        # check common attributes
        for qresult in self.qresults:
            for hit in qresult:
                assert qresult.id == hit.query_id
                for hsp in hit:
                    assert hit.id == hsp.hit_id
                    assert qresult.id == hsp.query_id

        # test first qresult
        qresult = self.qresults[0]
        assert "hg19_dna" == qresult.id
        assert "blat" == qresult.program
        assert 50 == qresult.seq_len
        assert 10 == len(qresult)
        # first qresult, first hit
        hit = qresult[0]
        assert "chr9" == hit.id
        assert 141213431 == hit.seq_len
        assert 1 == len(hit.hsps)
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 38 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 3 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert 1 == hsp[0].query_strand
        assert 9 == hsp.query_start
        assert 85737865 == hsp.hit_start
        assert 50 == hsp.query_end
        assert 85737906 == hsp.hit_end
        assert 1 == len(hsp)
        assert [41] == hsp.query_span_all
        assert [41] == hsp.hit_span_all
        assert [(9, 50)] == hsp.query_range_all
        assert [(85737865, 85737906)] == hsp.hit_range_all
        # first qresult, second hit
        hit = qresult[1]
        assert "chr8" == hit.id
        assert 146364022 == hit.seq_len
        assert 1 == len(hit.hsps)
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 41 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 0 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert 1 == hsp[0].query_strand
        assert 8 == hsp.query_start
        assert 95160479 == hsp.hit_start
        assert 49 == hsp.query_end
        assert 95160520 == hsp.hit_end
        assert 1 == len(hsp)
        assert [41] == hsp.query_span_all
        assert [41] == hsp.hit_span_all
        assert [(8, 49)] == hsp.query_range_all
        assert [(95160479, 95160520)] == hsp.hit_range_all
        # first qresult, third hit
        hit = qresult[2]
        assert "chr22" == hit.id
        assert 51304566 == hit.seq_len
        assert 2 == len(hit.hsps)
        # first qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        assert 33 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 3 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert 1 == hsp[0].query_strand
        assert 11 == hsp.query_start
        assert 42144400 == hsp.hit_start
        assert 47 == hsp.query_end
        assert 42144436 == hsp.hit_end
        assert 1 == len(hsp)
        assert [36] == hsp.query_span_all
        assert [36] == hsp.hit_span_all
        assert [(11, 47)] == hsp.query_range_all
        assert [(42144400, 42144436)] == hsp.hit_range_all
        # first qresult, third hit, second hsp
        hsp = qresult[2].hsps[1]
        assert 35 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 2 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert -1 == hsp[0].query_strand
        assert 12 == hsp.query_start
        assert 48997405 == hsp.hit_start
        assert 49 == hsp.query_end
        assert 48997442 == hsp.hit_end
        assert 1 == len(hsp)
        assert [37] == hsp.query_span_all
        assert [37] == hsp.hit_span_all
        assert [(12, 49)] == hsp.query_range_all
        assert [(48997405, 48997442)] == hsp.hit_range_all
        # first qresult, fourth hit
        hit = qresult[3]
        assert "chr2" == hit.id
        assert 243199373 == hit.seq_len
        assert 2 == len(hit.hsps)
        # first qresult, fourth hit, first hsp
        hsp = qresult[3].hsps[0]
        assert 43 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 1 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 1 == hsp.query_gapopen_num
        assert 4 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert 1 == hsp[0].query_strand
        assert 1 == hsp.query_start
        assert 183925984 == hsp.hit_start
        assert 49 == hsp.query_end
        assert 183926028 == hsp.hit_end
        assert 2 == len(hsp)
        assert [6, 38] == hsp.query_span_all
        assert [6, 38] == hsp.hit_span_all
        assert [(1, 7), (11, 49)] == hsp.query_range_all
        assert [(183925984, 183925990), (183925990, 183926028)] == hsp.hit_range_all
        # first qresult, fourth hit, second hsp
        hsp = qresult[3].hsps[1]
        assert 35 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 1 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert -1 == hsp[0].query_strand
        assert 13 == hsp.query_start
        assert 120641740 == hsp.hit_start
        assert 49 == hsp.query_end
        assert 120641776 == hsp.hit_end
        assert 1 == len(hsp)
        assert [36] == hsp.query_span_all
        assert [36] == hsp.hit_span_all
        assert [(13, 49)] == hsp.query_range_all
        assert [(120641740, 120641776)] == hsp.hit_range_all
        # first qresult, fifth hit
        hit = qresult[4]
        assert "chr19" == hit.id
        assert 59128983 == hit.seq_len
        assert 3 == len(hit.hsps)
        # first qresult, fifth hit, first hsp
        hsp = qresult[4].hsps[0]
        assert 34 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 2 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 1 == hsp.hit_gapopen_num
        assert 134 == hsp.hit_gap_num
        assert 1 == hsp[0].query_strand
        assert 10 == hsp.query_start
        assert 35483340 == hsp.hit_start
        assert 46 == hsp.query_end
        assert 35483510 == hsp.hit_end
        assert 2 == len(hsp)
        assert [25, 11] == hsp.query_span_all
        assert [25, 11] == hsp.hit_span_all
        assert [(10, 35), (35, 46)] == hsp.query_range_all
        assert [(35483340, 35483365), (35483499, 35483510)] == hsp.hit_range_all
        # first qresult, fifth hit, second hsp
        hsp = qresult[4].hsps[1]
        assert 39 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 0 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert -1 == hsp[0].query_strand
        assert 10 == hsp.query_start
        assert 54017130 == hsp.hit_start
        assert 49 == hsp.query_end
        assert 54017169 == hsp.hit_end
        assert 1 == len(hsp)
        assert [39] == hsp.query_span_all
        assert [39] == hsp.hit_span_all
        assert [(10, 49)] == hsp.query_range_all
        assert [(54017130, 54017169)] == hsp.hit_range_all
        # first qresult, fifth hit, third hsp
        hsp = qresult[4].hsps[2]
        assert 36 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 3 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert -1 == hsp[0].query_strand
        assert 10 == hsp.query_start
        assert 553742 == hsp.hit_start
        assert 49 == hsp.query_end
        assert 553781 == hsp.hit_end
        assert 1 == len(hsp)
        assert [39] == hsp.query_span_all
        assert [39] == hsp.hit_span_all
        assert [(10, 49)] == hsp.query_range_all
        assert [(553742, 553781)] == hsp.hit_range_all

    def test_psl_34_005(self, testf="psl_34_005.psl", pslx=False):
        """Test parsing blat output (psl_34_005.psl)."""
        blat_file = get_file(testf)
        self.qresults = list(parse(blat_file, FMT, pslx=pslx))
        assert 2 == len(self.qresults)
        # check common attributes
        for qresult in self.qresults:
            for hit in qresult:
                assert qresult.id == hit.query_id
                for hsp in hit:
                    assert hit.id == hsp.hit_id
                    assert qresult.id == hsp.query_id

        # test first qresult
        qresult = self.qresults[0]
        assert "hg18_dna" == qresult.id
        assert "blat" == qresult.program
        assert 33 == qresult.seq_len
        assert 3 == len(qresult)
        # first qresult, first hit
        hit = qresult[0]
        assert "chr4" == hit.id
        assert 191154276 == hit.seq_len
        assert 1 == len(hit.hsps)
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 16 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 0 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert 1 == hsp[0].query_strand
        assert 11 == hsp.query_start
        assert 61646095 == hsp.hit_start
        assert 27 == hsp.query_end
        assert 61646111 == hsp.hit_end
        assert 1 == len(hsp)
        assert [16] == hsp.query_span_all
        assert [16] == hsp.hit_span_all
        assert [(11, 27)] == hsp.query_range_all
        assert [(61646095, 61646111)] == hsp.hit_range_all
        # first qresult, second hit
        hit = qresult[1]
        assert "chr1" == hit.id
        assert 249250621 == hit.seq_len
        assert 1 == len(hit.hsps)
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 33 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 0 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert 1 == hsp[0].query_strand
        assert 0 == hsp.query_start
        assert 10271783 == hsp.hit_start
        assert 33 == hsp.query_end
        assert 10271816 == hsp.hit_end
        assert 1 == len(hsp)
        assert [33] == hsp.query_span_all
        assert [33] == hsp.hit_span_all
        assert [(0, 33)] == hsp.query_range_all
        assert [(10271783, 10271816)] == hsp.hit_range_all
        # first qresult, third hit
        hit = qresult[2]
        assert "chr2" == hit.id
        assert 243199373 == hit.seq_len
        assert 1 == len(hit.hsps)
        # first qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        assert 17 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 0 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert -1 == hsp[0].query_strand
        assert 8 == hsp.query_start
        assert 53575980 == hsp.hit_start
        assert 25 == hsp.query_end
        assert 53575997 == hsp.hit_end
        assert 1 == len(hsp)
        assert [17] == hsp.query_span_all
        assert [17] == hsp.hit_span_all
        assert [(8, 25)] == hsp.query_range_all
        assert [(53575980, 53575997)] == hsp.hit_range_all

        # test second qresult
        qresult = self.qresults[1]
        assert "hg19_dna" == qresult.id
        assert "blat" == qresult.program
        assert 50 == qresult.seq_len
        assert 10 == len(qresult)
        # second qresult, first hit
        hit = qresult[0]
        assert "chr9" == hit.id
        assert 141213431 == hit.seq_len
        assert 1 == len(hit.hsps)
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 38 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 3 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert 1 == hsp[0].query_strand
        assert 9 == hsp.query_start
        assert 85737865 == hsp.hit_start
        assert 50 == hsp.query_end
        assert 85737906 == hsp.hit_end
        assert 1 == len(hsp)
        assert [41] == hsp.query_span_all
        assert [41] == hsp.hit_span_all
        assert [(9, 50)] == hsp.query_range_all
        assert [(85737865, 85737906)] == hsp.hit_range_all
        # second qresult, second hit
        hit = qresult[1]
        assert "chr8" == hit.id
        assert 146364022 == hit.seq_len
        assert 1 == len(hit.hsps)
        # second qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 41 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 0 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert 1 == hsp[0].query_strand
        assert 8 == hsp.query_start
        assert 95160479 == hsp.hit_start
        assert 49 == hsp.query_end
        assert 95160520 == hsp.hit_end
        assert 1 == len(hsp)
        assert [41] == hsp.query_span_all
        assert [41] == hsp.hit_span_all
        assert [(8, 49)] == hsp.query_range_all
        assert [(95160479, 95160520)] == hsp.hit_range_all
        # second qresult, third hit
        hit = qresult[2]
        assert "chr22" == hit.id
        assert 51304566 == hit.seq_len
        assert 2 == len(hit.hsps)
        # second qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        assert 33 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 3 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert 1 == hsp[0].query_strand
        assert 11 == hsp.query_start
        assert 42144400 == hsp.hit_start
        assert 47 == hsp.query_end
        assert 42144436 == hsp.hit_end
        assert 1 == len(hsp)
        assert [36] == hsp.query_span_all
        assert [36] == hsp.hit_span_all
        assert [(11, 47)] == hsp.query_range_all
        assert [(42144400, 42144436)] == hsp.hit_range_all
        # second qresult, third hit, second hsp
        hsp = qresult[2].hsps[1]
        assert 35 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 2 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert -1 == hsp[0].query_strand
        assert 12 == hsp.query_start
        assert 48997405 == hsp.hit_start
        assert 49 == hsp.query_end
        assert 48997442 == hsp.hit_end
        assert 1 == len(hsp)
        assert [37] == hsp.query_span_all
        assert [37] == hsp.hit_span_all
        assert [(12, 49)] == hsp.query_range_all
        assert [(48997405, 48997442)] == hsp.hit_range_all
        # second qresult, fourth hit
        hit = qresult[3]
        assert "chr2" == hit.id
        assert 243199373 == hit.seq_len
        assert 2 == len(hit.hsps)
        # second qresult, fourth hit, first hsp
        hsp = qresult[3].hsps[0]
        assert 43 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 1 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 1 == hsp.query_gapopen_num
        assert 4 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert 1 == hsp[0].query_strand
        assert 1 == hsp.query_start
        assert 183925984 == hsp.hit_start
        assert 49 == hsp.query_end
        assert 183926028 == hsp.hit_end
        assert 2 == len(hsp)
        assert [6, 38] == hsp.query_span_all
        assert [6, 38] == hsp.hit_span_all
        assert [(1, 7), (11, 49)] == hsp.query_range_all
        assert [(183925984, 183925990), (183925990, 183926028)] == hsp.hit_range_all
        # second qresult, fourth hit, second hsp
        hsp = qresult[3].hsps[1]
        assert 35 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 1 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert -1 == hsp[0].query_strand
        assert 13 == hsp.query_start
        assert 120641740 == hsp.hit_start
        assert 49 == hsp.query_end
        assert 120641776 == hsp.hit_end
        assert 1 == len(hsp)
        assert [36] == hsp.query_span_all
        assert [36] == hsp.hit_span_all
        assert [(13, 49)] == hsp.query_range_all
        assert [(120641740, 120641776)] == hsp.hit_range_all
        # second qresult, fifth hit
        hit = qresult[4]
        assert "chr19" == hit.id
        assert 59128983 == hit.seq_len
        assert 3 == len(hit.hsps)
        # second qresult, fifth hit, first hsp
        hsp = qresult[4].hsps[0]
        assert 34 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 2 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 1 == hsp.hit_gapopen_num
        assert 134 == hsp.hit_gap_num
        assert 1 == hsp[0].query_strand
        assert 10 == hsp.query_start
        assert 35483340 == hsp.hit_start
        assert 46 == hsp.query_end
        assert 35483510 == hsp.hit_end
        assert 2 == len(hsp)
        assert [25, 11] == hsp.query_span_all
        assert [25, 11] == hsp.hit_span_all
        assert [(10, 35), (35, 46)] == hsp.query_range_all
        assert [(35483340, 35483365), (35483499, 35483510)] == hsp.hit_range_all
        # second qresult, fifth hit, second hsp
        hsp = qresult[4].hsps[1]
        assert 39 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 0 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert -1 == hsp[0].query_strand
        assert 10 == hsp.query_start
        assert 54017130 == hsp.hit_start
        assert 49 == hsp.query_end
        assert 54017169 == hsp.hit_end
        assert 1 == len(hsp)
        assert [39] == hsp.query_span_all
        assert [39] == hsp.hit_span_all
        assert [(10, 49)] == hsp.query_range_all
        assert [(54017130, 54017169)] == hsp.hit_range_all
        # second qresult, fifth hit, third hsp
        hsp = qresult[4].hsps[2]
        assert 36 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 3 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert -1 == hsp[0].query_strand
        assert 10 == hsp.query_start
        assert 553742 == hsp.hit_start
        assert 49 == hsp.query_end
        assert 553781 == hsp.hit_end
        assert 1 == len(hsp)
        assert [39] == hsp.query_span_all
        assert [39] == hsp.hit_span_all
        assert [(10, 49)] == hsp.query_range_all
        assert [(553742, 553781)] == hsp.hit_range_all

    def test_psl_35_001(self, testf="psl_35_001.psl", pslx=False):
        """Test parsing blat output (psl_35_001.psl)."""
        blat_file = get_file(testf)
        self.qresults = list(parse(blat_file, FMT, pslx=pslx))
        assert 1 == len(self.qresults)
        # check common attributes
        for qresult in self.qresults:
            for hit in qresult:
                assert qresult.id == hit.query_id
                for hsp in hit:
                    assert hit.id == hsp.hit_id
                    assert qresult.id == hsp.query_id

        # test first qresult
        qresult = self.qresults[0]
        assert "CAG33136.1" == qresult.id
        assert "blat" == qresult.program
        assert 230 == qresult.seq_len
        assert 2 == len(qresult)
        # first qresult, first hit
        hit = qresult[0]
        assert "chr13" == hit.id
        assert 114364328 == hit.seq_len
        assert 6 == len(hit.hsps)
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 52 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 0 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 0 == hsp.query_gapopen_num
        assert 0 == hsp.query_gap_num
        assert 0 == hsp.hit_gapopen_num
        assert 0 == hsp.hit_gap_num
        assert 0 == hsp[0].query_strand
        assert 61 == hsp.query_start
        assert 75566694 == hsp.hit_start
        assert 113 == hsp.query_end
        assert 75566850 == hsp.hit_end
        assert 1 == len(hsp)
        assert [52] == hsp.query_span_all
        assert [156] == hsp.hit_span_all
        assert [(61, 113)] == hsp.query_range_all
        assert [(75566694, 75566850)] == hsp.hit_range_all

    def test_psl_35_002(self, testf="psl_35_002.psl", pslx=False):
        """Test parsing blat output (psl_35_002.psl)."""
        blat_file = get_file(testf)
        self.qresults = list(parse(blat_file, FMT, pslx=pslx))
        assert 1 == len(self.qresults)
        # check common attributes
        for qresult in self.qresults:
            for hit in qresult:
                assert qresult.id == hit.query_id
                for hsp in hit:
                    assert hit.id == hsp.hit_id
                    assert qresult.id == hsp.query_id

        # test first qresult
        qresult = self.qresults[0]
        assert "CAG33136.1" == qresult.id
        assert "blat" == qresult.program
        assert 230 == qresult.seq_len
        assert 3 == len(qresult)
        # first qresult, last hit
        hit = qresult[-1]
        assert "KI537194" == hit.id
        assert 37111980 == hit.seq_len
        assert 1 == len(hit.hsps)
        # # first qresult, last hit, first hsp
        hsp = hit.hsps[-1]
        assert 204 == hsp.match_num
        assert 0 == hsp.match_rep_num
        assert 6 == hsp.mismatch_num
        assert 0 == hsp.n_num
        assert 1 == hsp.query_gapopen_num
        assert 20 == hsp.query_gap_num
        assert 1 == hsp.hit_gapopen_num
        assert 1 == hsp.hit_gap_num
        assert 0 == hsp[0].query_strand
        assert -1 == hsp[0].hit_strand
        assert 0 == hsp.query_start
        assert 20872390 == hsp.hit_start
        assert 230 == hsp.query_end
        assert 20873021 == hsp.hit_end
        assert 2 == len(hsp)
        assert [183, 27] == hsp.query_span_all
        assert [549, 81] == hsp.hit_span_all
        assert [(0, 183), (203, 230)] == hsp.query_range_all
        assert [(20872472, 20873021), (20872390, 20872471)] == hsp.hit_range_all


class BlatPslxCases(BlatPslCases):
    def test_pslx_34_001(self, testf="pslx_34_001.pslx"):
        """Test parsing blat output (pslx_34_001.pslx)."""
        BlatPslCases.test_psl_34_001(self, "pslx_34_001.pslx", pslx=True)

        # test first qresult
        qresult = self.qresults[0]
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert "aggtaaactgccttca" == hsp.query_all[0].seq
        assert "aggtaaactgccttca" == hsp.hit_all[0].seq
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert "atgagcttccaaggtaaactgccttcaagattc" == hsp.query_all[0].seq
        assert "atgagcttccaaggtaaactgccttcaagattc" == hsp.hit_all[0].seq
        # first qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        assert "aaggcagtttaccttgg" == hsp.query_all[0].seq
        assert "aaggcagtttaccttgg" == hsp.hit_all[0].seq

        # test second qresult
        qresult = self.qresults[1]
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert "acaaaggggctgggcgtggtggctcacacctgtaatcccaa" == hsp.query_all[0].seq
        assert "acaaaggggctgggcgcagtggctcacgcctgtaatcccaa" == hsp.hit_all[0].seq
        # second qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert "cacaaaggggctgggcgtggtggctcacacctgtaatccca" == hsp.query_all[0].seq
        assert "cacaaaggggctgggcgtggtggctcacacctgtaatccca" == hsp.hit_all[0].seq
        # second qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        assert "aaaggggctgggcgtggtggctcacacctgtaatcc" == hsp.query_all[0].seq
        assert "aaaggggctgggcgtggtagctcatgcctgtaatcc" == hsp.hit_all[0].seq
        # second qresult, third hit, second hsp
        hsp = qresult[2].hsps[1]
        assert "tgggattacaggtgtgagccaccacgcccagcccctt" == hsp.query_all[0].seq
        assert "tgggattacaggcgggagccaccacgcccagcccctt" == hsp.hit_all[0].seq
        # second qresult, fourth hit, first hsp
        hsp = qresult[3].hsps[0]
        assert "aaaaat" == hsp.query_all[0].seq
        assert "aaaaat" == hsp.hit_all[0].seq
        assert "aaaggggctgggcgtggtggctcacacctgtaatccca" == hsp.query_all[1].seq
        assert "aaaggggctgggcgtggtggctcacgcctgtaatccca" == hsp.hit_all[1].seq
        # second qresult, fourth hit, second hsp
        hsp = qresult[3].hsps[1]
        assert "tgggattacaggtgtgagccaccacgcccagcccct" == hsp.query_all[0].seq
        assert "tgggattacaggcgtgagccaccacgcccagcccct" == hsp.hit_all[0].seq
        # second qresult, fifth hit, first hsp
        hsp = qresult[4].hsps[0]
        assert "caaaggggctgggcgtggtggctca" == hsp.query_all[0].seq
        assert "caaaggggctgggcgtagtggctga" == hsp.hit_all[0].seq
        assert "cacctgtaatc" == hsp.query_all[1].seq
        assert "cacctgtaatc" == hsp.hit_all[1].seq
        # second qresult, fifth hit, second hsp
        hsp = qresult[4].hsps[1]
        assert "tgggattacaggtgtgagccaccacgcccagcccctttg" == hsp.query_all[0].seq
        assert "tgggattacaggtgtgagccaccacgcccagcccctttg" == hsp.hit_all[0].seq
        # second qresult, fifth hit, third hsp
        hsp = qresult[4].hsps[2]
        assert "tgggattacaggtgtgagccaccacgcccagcccctttg" == hsp.query_all[0].seq
        assert "tgggatgacaggggtgaggcaccacgcccagcccctttg" == hsp.hit_all[0].seq

    def test_pslx_34_002(self, testf="pslx_34_002.pslx"):
        """Test parsing blat output (pslx_34_002.pslx)."""
        BlatPslCases.test_psl_34_002(self, "pslx_34_002.pslx", pslx=True)

    def test_pslx_34_003(self, testf="pslx_34_003.pslx"):
        """Test parsing blat output (pslx_34_003.pslx)."""
        BlatPslCases.test_psl_34_003(self, "pslx_34_003.pslx", pslx=True)

        # test first qresult
        qresult = self.qresults[0]
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert "aggtaaactgccttca" == hsp.query_all[0].seq
        assert "aggtaaactgccttca" == hsp.hit_all[0].seq
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert "atgagcttccaaggtaaactgccttcaagattc" == hsp.query_all[0].seq
        assert "atgagcttccaaggtaaactgccttcaagattc" == hsp.hit_all[0].seq
        # first qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        assert "aaggcagtttaccttgg" == hsp.query_all[0].seq
        assert "aaggcagtttaccttgg" == hsp.hit_all[0].seq

    def test_pslx_34_004(self, testf="pslx_34_004.pslx"):
        """Test parsing blat output (pslx_34_004.pslx)."""
        BlatPslCases.test_psl_34_004(self, "pslx_34_004.pslx", pslx=True)

        # test first qresult
        qresult = self.qresults[0]
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert "acaaaggggctgggcgtggtggctcacacctgtaatcccaa" == hsp.query_all[0].seq
        assert "acaaaggggctgggcgcagtggctcacgcctgtaatcccaa" == hsp.hit_all[0].seq
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert "cacaaaggggctgggcgtggtggctcacacctgtaatccca" == hsp.query_all[0].seq
        assert "cacaaaggggctgggcgtggtggctcacacctgtaatccca" == hsp.hit_all[0].seq
        # first qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        assert "aaaggggctgggcgtggtggctcacacctgtaatcc" == hsp.query_all[0].seq
        assert "aaaggggctgggcgtggtagctcatgcctgtaatcc" == hsp.hit_all[0].seq
        # first qresult, third hit, second hsp
        hsp = qresult[2].hsps[1]
        assert "tgggattacaggtgtgagccaccacgcccagcccctt" == hsp.query_all[0].seq
        assert "tgggattacaggcgggagccaccacgcccagcccctt" == hsp.hit_all[0].seq
        # first qresult, fourth hit, first hsp
        hsp = qresult[3].hsps[0]
        assert "aaaaat" == hsp.query_all[0].seq
        assert "aaaaat" == hsp.hit_all[0].seq
        assert "aaaggggctgggcgtggtggctcacacctgtaatccca" == hsp.query_all[1].seq
        assert "aaaggggctgggcgtggtggctcacgcctgtaatccca" == hsp.hit_all[1].seq
        # first qresult, fourth hit, second hsp
        hsp = qresult[3].hsps[1]
        assert "tgggattacaggtgtgagccaccacgcccagcccct" == hsp.query_all[0].seq
        assert "tgggattacaggcgtgagccaccacgcccagcccct" == hsp.hit_all[0].seq
        # first qresult, fifth hit, first hsp
        hsp = qresult[4].hsps[0]
        assert "caaaggggctgggcgtggtggctca" == hsp.query_all[0].seq
        assert "caaaggggctgggcgtagtggctga" == hsp.hit_all[0].seq
        assert "cacctgtaatc" == hsp.query_all[1].seq
        assert "cacctgtaatc" == hsp.hit_all[1].seq
        # first qresult, fifth hit, second hsp
        hsp = qresult[4].hsps[1]
        assert "tgggattacaggtgtgagccaccacgcccagcccctttg" == hsp.query_all[0].seq
        assert "tgggattacaggtgtgagccaccacgcccagcccctttg" == hsp.hit_all[0].seq
        # first qresult, fifth hit, third hsp
        hsp = qresult[4].hsps[2]
        assert "tgggattacaggtgtgagccaccacgcccagcccctttg" == hsp.query_all[0].seq
        assert "tgggatgacaggggtgaggcaccacgcccagcccctttg" == hsp.hit_all[0].seq

    def test_pslx_34_005(self, testf="pslx_34_005.pslx"):
        """Test parsing blat output (pslx_34_005.pslx)."""
        BlatPslCases.test_psl_34_005(self, "pslx_34_005.pslx", pslx=True)

        # test first qresult
        qresult = self.qresults[0]
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert "aggtaaactgccttca" == hsp.query_all[0].seq
        assert "aggtaaactgccttca" == hsp.hit_all[0].seq
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert "atgagcttccaaggtaaactgccttcaagattc" == hsp.query_all[0].seq
        assert "atgagcttccaaggtaaactgccttcaagattc" == hsp.hit_all[0].seq
        # first qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        assert "aaggcagtttaccttgg" == hsp.query_all[0].seq
        assert "aaggcagtttaccttgg" == hsp.hit_all[0].seq

        # test second qresult
        qresult = self.qresults[1]
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert "acaaaggggctgggcgtggtggctcacacctgtaatcccaa" == hsp.query_all[0].seq
        assert "acaaaggggctgggcgcagtggctcacgcctgtaatcccaa" == hsp.hit_all[0].seq
        # second qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert "cacaaaggggctgggcgtggtggctcacacctgtaatccca" == hsp.query_all[0].seq
        assert "cacaaaggggctgggcgtggtggctcacacctgtaatccca" == hsp.hit_all[0].seq
        # second qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        assert "aaaggggctgggcgtggtggctcacacctgtaatcc" == hsp.query_all[0].seq
        assert "aaaggggctgggcgtggtagctcatgcctgtaatcc" == hsp.hit_all[0].seq
        # second qresult, third hit, second hsp
        hsp = qresult[2].hsps[1]
        assert "tgggattacaggtgtgagccaccacgcccagcccctt" == hsp.query_all[0].seq
        assert "tgggattacaggcgggagccaccacgcccagcccctt" == hsp.hit_all[0].seq
        # second qresult, fourth hit, first hsp
        hsp = qresult[3].hsps[0]
        assert "aaaaat" == hsp.query_all[0].seq
        assert "aaaaat" == hsp.hit_all[0].seq
        assert "aaaggggctgggcgtggtggctcacacctgtaatccca" == hsp.query_all[1].seq
        assert "aaaggggctgggcgtggtggctcacgcctgtaatccca" == hsp.hit_all[1].seq
        # second qresult, fourth hit, second hsp
        hsp = qresult[3].hsps[1]
        assert "tgggattacaggtgtgagccaccacgcccagcccct" == hsp.query_all[0].seq
        assert "tgggattacaggcgtgagccaccacgcccagcccct" == hsp.hit_all[0].seq
        # second qresult, fifth hit, first hsp
        hsp = qresult[4].hsps[0]
        assert "caaaggggctgggcgtggtggctca" == hsp.query_all[0].seq
        assert "caaaggggctgggcgtagtggctga" == hsp.hit_all[0].seq
        assert "cacctgtaatc" == hsp.query_all[1].seq
        assert "cacctgtaatc" == hsp.hit_all[1].seq
        # second qresult, fifth hit, second hsp
        hsp = qresult[4].hsps[1]
        assert "tgggattacaggtgtgagccaccacgcccagcccctttg" == hsp.query_all[0].seq
        assert "tgggattacaggtgtgagccaccacgcccagcccctttg" == hsp.hit_all[0].seq
        # second qresult, fifth hit, third hsp
        hsp = qresult[4].hsps[2]
        assert "tgggattacaggtgtgagccaccacgcccagcccctttg" == hsp.query_all[0].seq
        assert "tgggatgacaggggtgaggcaccacgcccagcccctttg" == hsp.hit_all[0].seq

    def test_pslx_35_002(self, testf="pslx_35_002.pslx"):
        """Test parsing blat output (pslx_35_002.pslx)."""
        BlatPslCases.test_psl_35_002(self, "pslx_35_002.pslx", pslx=True)

        # first qresult, last hit, first hsp
        qresult = self.qresults[0]
        hsp = qresult[-1].hsps[0]

        assert "MEGQRWLPLEANPEVTNQFLKQLGLHPNWQFVDVY" == hsp.query_all[0].seq[:35]
        assert "ETSAHEGQTEAPSIDEKVDLHFIALVHVDGHLYEL" == hsp.query_all[0].seq[-35:]
        assert "DAIEVCKKFMERDPDELRFNAIALSAA" == hsp.query_all[1].seq

        assert "MESQRWLPLEANPEVTNQFLKQLGLHPNWQCVDVY" == hsp.hit_all[0].seq[:35]
        assert "ETSAHEGQTEAPNIDEKVDLHFIALVHVDGHLYEL" == hsp.hit_all[0].seq[-35:]
        assert "DAIEVCKKFMERDPDELRFNAIALSAA" == hsp.hit_all[1].seq


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
