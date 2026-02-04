# Copyright 2024 by Samuel Prince. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.


"""Tests for SearchIO InfernalIO infernal-tab parser."""

import os
import unittest
import pytest
import itertools

from Bio.SearchIO import parse

# test case files are in the Infernal directory
TEST_DIR = "Infernal"
FMT = "infernal-tab"


def get_file(filename):
    """Return the path of a test file."""
    return os.path.join(TEST_DIR, filename)


def next_result(qresults, counter):
    """Iterate over the results and counter."""
    return next(qresults), next(counter)


class CmscanCases(unittest.TestCase):
    """Test parsing cmscan output."""

    def test_cmscan_mq_mm(self):
        """Test parsing infernal-tab, cmscan, multiple queries, multiple hit, one hsp, default format (IRES_5S_U2_Yeast)"""
        tab_file = get_file("cmscan_115_IRES_5S_U2_Yeast.tbl")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        # first qresult
        qresult, count = next_result(qresults, counter)
        assert len(qresult) == 1
        assert qresult.id == "ENA|BK006935|BK006935.2"
        assert qresult.accession == "-"
        hit = qresult[0]
        assert len(hit) == 2
        assert hit.id == "U2"
        assert hit.accession == "RF00004"
        assert hit.description == "U2 spliceosomal RNA"
        # first hsp
        hsp = hit[0]
        assert len(hsp) == 1
        assert hsp.evalue == 0.91
        assert hsp.bitscore == 13.5
        assert hsp.bias == 0.0
        assert hsp.gc == 0.44
        assert hsp.truncated == "no"
        assert hsp.model == "cm"
        assert hsp.pipeline_pass == 1
        assert not hsp.is_included
        frag = hsp[0]
        assert frag.query_start == 1
        assert frag.query_end == 193
        assert 52929 == frag.hit_start
        assert 53083 == frag.hit_end
        assert frag.hit_strand == 0
        # second hsp
        hsp = hit[1]
        assert len(hsp) == 1
        assert 1.3 == hsp.evalue
        assert 12.8 == hsp.bitscore
        assert 5.3 == hsp.bias
        assert hsp.gc == 0.33
        assert hsp.truncated == "no"
        assert hsp.pipeline_pass == 1
        assert not hsp.is_included
        frag = hsp[0]
        assert frag.query_start == 1
        assert frag.query_end == 193
        assert frag.hit_start == 196389
        assert frag.hit_end == 196571
        assert frag.hit_strand == -1

        # second qresult
        qresult, count = next_result(qresults, counter)
        assert len(qresult) == 1
        assert "ENA|BK006936|BK006936.2" == qresult.id
        assert qresult.accession == "-"
        hit = qresult[0]
        assert len(hit) == 1
        assert hit.id == "U2"
        assert hit.accession == "RF00004"
        assert hit.description == "U2 spliceosomal RNA"
        hsp = hit[0]
        assert len(hsp) == 1
        assert hsp.evalue == 1.2e-20
        assert hsp.bitscore == 98.7
        assert hsp.bias == 0.1
        assert hsp.gc == 0.33
        assert hsp.truncated == "no"
        assert hsp.model == "cm"
        assert hsp.pipeline_pass == 1
        assert hsp.is_included
        frag = hsp[0]
        assert frag.query_start == 1
        assert frag.query_end == 193
        assert frag.hit_start == 681747
        assert frag.hit_end == 681858
        assert frag.hit_strand == -1

        # third qresult
        qresult, count = next_result(qresults, counter)
        assert len(qresult) == 2
        assert qresult.id == "ENA|BK006937|BK006937.2"
        assert qresult.accession == "-"
        # first hit
        hit = qresult[0]
        assert len(hit) == 1
        assert hit.id == "5S_rRNA"
        assert hit.accession == "RF00001"
        assert hit.description == "5S ribosomal RNA"
        hsp = hit[0]
        assert len(hsp) == 1
        assert hsp.evalue == 2.4
        assert hsp.bitscore == 14.1
        assert hsp.bias == 0.3
        assert hsp.gc == 0.41
        assert hsp.truncated == "no"
        assert hsp.model == "cm"
        assert hsp.pipeline_pass == 1
        assert not hsp.is_included
        frag = hsp[0]
        assert frag.query_start == 1
        assert frag.query_end == 119
        assert frag.hit_start == 644
        assert frag.hit_end == 761
        assert frag.hit_strand == -1
        # second hit
        hit = qresult[1]
        assert len(hit) == 1
        assert hit.id == "U2"
        assert hit.accession == "RF00004"
        assert hit.description == "U2 spliceosomal RNA"
        hsp = hit[0]
        assert len(hsp) == 1
        assert hsp.evalue == 4.7
        assert hsp.bitscore == 11.1
        assert hsp.bias == 0.1
        assert hsp.gc == 0.32
        assert hsp.truncated == "no"
        assert hsp.pipeline_pass == 1
        assert not hsp.is_included
        frag = hsp[0]
        assert frag.query_start == 1
        assert frag.query_end == 193
        assert frag.hit_start == 229885
        assert frag.hit_end == 229986
        assert frag.hit_strand == -1

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert count == 3

    def test_cmscan_mq_mm_fmt2(self):
        """Test parsing infernal-tab, cmscan, multiple queries, multiple hit, one hsp, fmt 2 (IRES_5S_U2_Yeast_fmt_2)"""
        tab_file = get_file("cmscan_115_IRES_5S_U2_Yeast_fmt_2.tbl")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        # first qresult
        qresult, count = next_result(qresults, counter)
        assert len(qresult) == 1
        assert qresult.id == "ENA|BK006936|BK006936.2"
        assert qresult.accession == "-"
        assert qresult.clan == "-"
        assert qresult.seq_len == 813184
        hit = qresult[0]
        assert len(hit) == 1
        assert hit.id == "U2"
        assert hit.accession == "RF00004"
        assert hit.description == "U2 spliceosomal RNA"
        assert hit.seq_len == 193
        # first hsp
        hsp = hit[0]
        assert len(hsp) == 1
        assert hsp.evalue == 1.2e-20
        assert hsp.bitscore == 98.7
        assert hsp.bias == 0.1
        assert hsp.gc == 0.33
        assert hsp.truncated == "no"
        assert hsp.model == "cm"
        assert hsp.pipeline_pass == 1
        assert hsp.is_included
        assert "*" == hsp.olp
        assert None == hsp.anyidx
        assert None == hsp.afrct1
        assert None == hsp.afrct2
        assert None == hsp.winidx
        assert None == hsp.wfrct1
        assert None == hsp.wfrct2
        frag = hsp[0]
        assert frag.query_start == 1
        assert frag.query_end == 193
        assert frag.hit_start == 681747
        assert frag.hit_end == 681858
        assert frag.hit_strand == -1

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert count == 1

    def test_cmscan_mq_mm_fmt2_clan(self):
        """Test parsing infernal-tab, cmscan, mulitple queries, multiple hit, one hsp, fmt 2, multiple clan (SSU_clan_fmt_2)"""
        tab_file = get_file("cmscan_115_SSU_clan-fmt_2.tbl")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        # first qresult
        qresult, count = next_result(qresults, counter)
        assert len(qresult) == 2
        assert qresult.id == "ENA|BK006945|BK006945.2"
        assert qresult.accession == "-"
        assert qresult.clan == "CL00111"
        assert qresult.seq_len == 1078177
        hit = qresult[0]
        assert len(hit) == 1
        assert hit.id == "SSU_rRNA_eukarya"
        assert hit.accession == "RF01960"
        assert hit.description == "Eukaryotic small subunit ribosomal RNA"
        assert hit.seq_len == 1831
        hsp = hit[0]
        assert len(hsp) == 1
        assert hsp.evalue == 0
        assert hsp.bitscore == 1817.8
        assert hsp.bias == 9.5
        assert hsp.gc == 0.45
        assert hsp.truncated == "no"
        assert hsp.model == "cm"
        assert hsp.pipeline_pass == 1
        assert hsp.is_included
        assert "^" == hsp.olp
        assert None == hsp.anyidx
        assert None == hsp.afrct1
        assert None == hsp.afrct2
        assert None == hsp.winidx
        assert None == hsp.wfrct1
        assert None == hsp.wfrct2
        frag = hsp[0]
        assert frag.query_start == 1
        assert frag.query_end == 1831
        assert frag.hit_start == 455933
        assert frag.hit_end == 457732
        assert frag.hit_strand == -1

        # second hit (overlapping with the first hit)
        hit = qresult[1]
        assert len(hit) == 1
        assert hit.id == "SSU_rRNA_archaea"
        assert hit.accession == "RF01959"
        assert hit.description == "Archaeal small subunit ribosomal RNA"
        assert hit.seq_len == 1478
        # first hsp
        hsp = hit[0]
        assert len(hsp) == 1
        assert hsp.evalue == 5.8e-187
        assert hsp.bitscore == 633.9
        assert hsp.bias == 9.7
        assert hsp.gc == 0.45
        assert hsp.truncated == "no"
        assert hsp.model == "cm"
        assert hsp.pipeline_pass == 1
        assert hsp.is_included
        assert "=" == hsp.olp
        assert 1 == hsp.anyidx
        assert 0.998 == hsp.afrct1
        assert 1.000 == hsp.afrct2
        assert None == hsp.winidx
        assert None == hsp.wfrct1
        assert None == hsp.wfrct2
        frag = hsp[0]
        assert frag.query_start == 1
        assert frag.query_end == 1478
        assert frag.hit_start == 455930
        assert frag.hit_end == 457732
        assert frag.hit_strand == -1

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert count == 1

    def test_cmscan_mq_mm_fmt3(self):
        """Test parsing infernal-tab, cmscan, multiple queries, multiple hit, one hsp, fmt 3 (IRES_5S_U2_Yeast_fmt_3)"""
        tab_file = get_file("cmscan_115_IRES_5S_U2_Yeast_fmt_3.tbl")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        # first qresult
        qresult, count = next_result(qresults, counter)
        assert len(qresult) == 1
        assert "ENA|BK006936|BK006936.2" == qresult.id
        assert qresult.accession == "-"
        assert 813184 == qresult.seq_len
        hit = qresult[0]
        assert len(hit) == 1
        assert hit.id == "U2"
        assert hit.accession == "RF00004"
        assert hit.description == "U2 spliceosomal RNA"
        assert 193 == hit.seq_len
        # first hsp
        hsp = hit[0]
        assert len(hsp) == 1
        assert hsp.evalue == 1.2e-20
        assert hsp.bitscore == 98.7
        assert hsp.bias == 0.1
        assert hsp.gc == 0.33
        assert hsp.truncated == "no"
        assert hsp.model == "cm"
        assert hsp.pipeline_pass == 1
        assert hsp.is_included
        frag = hsp[0]
        assert frag.query_start == 1
        assert frag.query_end == 193
        assert frag.hit_start == 681747
        assert frag.hit_end == 681858
        assert frag.hit_strand == -1

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert count == 1


class CmsearchCases(unittest.TestCase):
    """Test parsing cmsearch output."""

    def test_1q_0m(self):
        """Test parsing infernal-tab, cmsearch, one query, no hits (IRES_Yeast)"""
        tab_file = get_file("cmsearch_114_IRES_Yeast.tbl")
        qresults = parse(tab_file, FMT)

        with pytest.raises(StopIteration):
            next(qresults)

    def test_cmsearch_1q_1m(self):
        """Test parsing infernal-tab, cmsearch, one queries, one hit, one hsp (U2_Yeast)"""
        tab_file = get_file("cmsearch_114_U2_Yeast.tbl")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, count = next_result(qresults, counter)
        assert len(qresult) == 1
        assert qresult.id == "U2"
        assert qresult.accession == "RF00004"
        hit = qresult[0]
        assert len(hit) == 1
        assert hit.id == "ENA|BK006936|BK006936.2"
        assert hit.accession == "-"
        assert hit.description == "TPA_inf: Saccharomyces cerevisiae S288C chromosome II, complete sequence."

        hsp = hit[0]
        assert len(hsp) == 1
        assert hsp.evalue == 5.9e-20
        assert hsp.bitscore == 98.7
        assert hsp.bias == 0.1
        assert hsp.is_included
        assert hsp.gc == 0.33
        assert hsp.truncated == "no"
        assert hsp.model == "cm"
        assert hsp.pipeline_pass == 1
        frag = hsp[0]
        assert frag.query_start == 1
        assert frag.query_end == 193
        assert frag.hit_start == 681747
        assert frag.hit_end == 681858
        assert frag.hit_strand == -1

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert count == 1

    def test_cmsearch_1q_mm(self):
        """Test parsing infernal-tab, cmsearch, one queries, multiple hit, one hsp (5S_Yeast)"""
        tab_file = get_file("cmsearch_114_5S_Yeast.tbl")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, count = next_result(qresults, counter)
        assert len(qresult) == 1
        assert qresult.id == "5S_rRNA"
        assert qresult.accession == "RF00001"
        # first hit
        hit = qresult[0]
        assert 6 == len(hit)
        assert hit.id == "ENA|BK006945|BK006945.2"
        assert hit.accession == "-"
        assert hit.description == "TPA_inf: Saccharomyces cerevisiae S288C chromosome XII, complete sequence."
        hsp = hit[0]
        assert len(hsp) == 1
        assert hsp.evalue == 1.6e-18
        assert hsp.bitscore == 88.8
        assert hsp.bias == 0.0
        assert hsp.gc == 0.52
        assert hsp.truncated == "no"
        assert hsp.model == "cm"
        assert hsp.pipeline_pass == 1
        assert hsp.is_included
        frag = hsp[0]
        assert frag.query_start == 1
        assert frag.query_end == 119
        assert frag.hit_start == 459676
        assert frag.hit_end == 459796
        assert frag.hit_strand == 0
        # last hit
        hsp = hit[-1]
        assert len(hsp) == 1
        assert hsp.evalue == 4.4e-17
        assert hsp.bitscore == 83.2
        assert hsp.bias == 0.0
        assert hsp.gc == 0.53
        assert hsp.truncated == "no"
        assert hsp.pipeline_pass == 1
        assert hsp.is_included
        frag = hsp[0]
        assert frag.query_start == 1
        assert frag.query_end == 119
        assert frag.hit_start == 485697
        assert frag.hit_end == 485817
        assert frag.hit_strand == 0

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert count == 1

    def test_cmsearch_1q_mm_shuf(self):
        """Test parsing infernal-tab, cmsearch, one queries, multiple non-consecutive hits, one hsp (U2_Yeast_full_shuffled)"""
        tab_file = get_file("cmsearch_114_U2_Yeast_full_shuffled.tbl")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, count = next_result(qresults, counter)
        assert 2 == len(qresult)
        assert qresult.id == "U2"
        assert qresult.accession == "RF00004"
        # first hit
        # first hit (3 hsps at rank 1,3 and 4)
        hit = qresult[0]
        assert 3 == len(hit)
        assert hit.id == "ENA|BK006936|BK006936.2"
        assert hit.description == "-"
        assert hit.query_id == "U2"
        # first hsp (rank 1)
        hsp = hit[0]
        assert hsp.query_start == 1
        assert hsp.query_end == 193
        assert hsp.hit_start == 681747
        assert hsp.hit_end == 681858
        # second hsp (rank 3)
        hsp = hit[1]
        assert hsp.query_start == 1
        assert hsp.query_end == 193
        assert hsp.hit_start == 1370418
        assert hsp.hit_end == 1370563
        # last hsp (rank 4)
        hsp = hit[2]
        assert hsp.query_start == 1
        assert hsp.query_end == 193
        assert hsp.hit_start == 1079243
        assert hsp.hit_end == 1079392
        # second hit
        hit = qresult[1]
        assert 3 == len(hit)
        assert hit.id == "ENA|BK006948|BK006948.2"
        assert hit.description == "-"
        assert hit.query_id == "U2"
        # first hsp (rank 2)
        hsp = hit[0]
        assert hsp.query_start == 1
        assert hsp.query_end == 193
        assert hsp.hit_start == 737324
        assert hsp.hit_end == 737498
        # second hsp (rank 5)
        hsp = hit[1]
        assert hsp.query_start == 1
        assert hsp.query_end == 193
        assert hsp.hit_start == 425490
        assert hsp.hit_end == 425693
        # last hsp (rank 6)
        hsp = hit[2]
        assert hsp.query_start == 1
        assert hsp.query_end == 193
        assert hsp.hit_start == 1073786
        assert hsp.hit_end == 1073950

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert count == 1


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
