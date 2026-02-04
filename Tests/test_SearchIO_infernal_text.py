# Copyright 2024 by Samuel Prince. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Tests for SearchIO InfernalIO infernal-text parser."""

import os
import unittest
import pytest
import itertools

from Bio.SearchIO import parse

# test case files are in the Infernal directory
TEST_DIR = "Infernal"
FMT = "infernal-text"


def get_file(filename):
    """Return the path of a test file."""
    return os.path.join(TEST_DIR, filename)


def next_result(qresults, counter):
    """Iterate over the results and counter."""
    return next(qresults), next(counter)


class CmscanCases(unittest.TestCase):
    """Test parsing cmsearch output."""

    def test_cmscan(self):
        """Test parsing infernal-text, cmscan, multiple queries"""
        tab_file = get_file("cmscan_115_IRES_5S_U2_Yeast.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        # first qresult
        qresult, count = next_result(qresults, counter)
        assert len(qresult) == 1
        assert qresult.id == "ENA|BK006935|BK006935.2"
        assert qresult.seq_len == 230218
        assert qresult.description == "<unknown description>"
        assert qresult.program == "cmscan"
        assert qresult.version == "1.1.5"
        assert qresult.target == "IRES_5S_U2.cm"
        # first hit
        hit = qresult[0]
        assert len(hit) == 2
        assert hit.id == "U2"
        assert hit.description == "U2 spliceosomal RNA"
        assert hit.query_id == "ENA|BK006935|BK006935.2"
        hsp = hit[0]
        assert len(hsp) == 2
        assert hsp.model == "cm"
        assert hsp.truncated == "no"
        assert hsp.gc == 0.44
        assert hsp.evalue == 0.91
        assert hsp.bitscore == 13.5
        assert hsp.bias == 0.0
        assert not hsp.is_included
        assert hsp.query_start == 1
        assert hsp.query_end == 193
        assert hsp.query_endtype == "[]"
        assert hsp.hit_start == 52929
        assert hsp.hit_end == 53083
        assert hsp.hit_endtype == ".."
        assert hsp.avg_acc == 0.80
        # first fragment
        frag = hsp[0]
        assert frag.query_start == 1
        assert frag.query_end == 46
        assert frag.hit_start == 52929
        assert frag.hit_end == 52974
        assert frag.hit_strand == 0
        # second fragment
        frag = hsp[1]
        assert frag.query_start == 84
        assert frag.query_end == 193
        assert frag.hit_start == 52977
        assert frag.hit_end == 53083
        assert frag.hit_strand == 0
        # second hit
        hsp = hit[1]
        assert hsp.query_start == 1
        assert hsp.query_end == 193
        assert hsp.hit_start == 196389
        assert hsp.hit_end == 196571
        # second qresult
        qresult, count = next_result(qresults, counter)
        assert len(qresult) == 1
        assert qresult.id == "ENA|BK006936|BK006936.2"
        # first hit
        hit = qresult[0]
        assert len(hit) == 1
        assert hit.id == "U2"
        assert hit.description == "U2 spliceosomal RNA"
        hsp = hit[0]
        assert hsp.query_start == 1
        assert hsp.query_end == 193
        assert hsp.hit_start == 681747
        assert hsp.hit_end == 681858
        # third (last) qresult
        qresult, count = next_result(qresults, counter)
        assert qresult.id == "ENA|BK006937|BK006937.2"

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 3 == count


class CmsearchCases(unittest.TestCase):
    """Test parsing cmsearch output."""

    def test_cmsearch_1q_0m(self):
        """Test parsing infernal-text, cmsearch, one query, no hits (IRES_Yeast)"""
        text_file = get_file("cmsearch_114_IRES_Yeast.txt")
        qresults = parse(text_file, FMT)
        counter = itertools.count(start=1)

        qresult, count = next_result(qresults, counter)
        assert len(qresult) == 0
        assert qresult.id == "IRES_HCV"
        assert qresult.seq_len == 352
        assert qresult.accession == "RF00061"
        assert qresult.description == "Hepatitis C virus internal ribosome entry site"
        assert qresult.program == "cmsearch"
        assert qresult.version == "1.1.4"
        assert qresult.target == "GCA_000146045.2.fasta"

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == count

    def test_cmsearch_1q_1m_1h(self):
        """Test parsing infernal-text, cmsearch, one queries, one hit, one hsp, multiple fragments (U2_Yeast)"""
        tab_file = get_file("cmsearch_114_U2_Yeast.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, count = next_result(qresults, counter)
        assert len(qresult) == 1
        assert qresult.id == "U2"
        assert qresult.seq_len == 193
        assert qresult.accession == "RF00004"
        assert qresult.description == "U2 spliceosomal RNA"
        assert qresult.program == "cmsearch"
        assert qresult.version == "1.1.4"
        assert qresult.target == "GCA_000146045.2.fasta"
        hit = qresult[0]
        assert len(hit) == 1
        assert hit.id == "ENA|BK006936|BK006936.2"
        assert hit.description == "TPA_inf: Saccharomyces cerevisiae S288C chromosome II, complete sequence."
        assert hit.query_id == "U2"
        hsp = hit[0]
        assert len(hsp) == 2
        assert hsp.model == "cm"
        assert hsp.truncated == "no"
        assert hsp.gc == 0.33
        assert hsp.evalue == 5.9e-20
        assert hsp.bitscore == 98.7
        assert hsp.bias == 0.1
        assert hsp.is_included
        assert hsp.query_start == 1
        assert hsp.query_end == 193
        assert hsp.query_endtype == "[]"
        assert hsp.hit_start == 681747
        assert hsp.hit_end == 681858
        assert hsp.hit_endtype == ".."
        assert hsp.avg_acc == 0.91
        # first fragment
        frag = hsp[0]
        assert frag.query_start == 1
        assert frag.query_end == 112
        assert frag.hit_start == 681763
        assert frag.hit_end == 681858
        assert frag.hit_strand == -1
        assert frag.query.seq == "AUacCUUCucgGCcUUUUgGCuaaGAUCAAGUGUAGUAUCUGUUCUUauCAGUuUAAuAuCUGauAuggcccccAuugggggccaau-uauaUUAaauuaAUUUUUggaacua"
        assert frag.hit.seq == "AUC---UCUUUGCCUUUUGGCUUAGAUCAAGUGUAGUAUCUGUUCUUUUCAGUGUAACAACUGAAAUGA-CCUCAAUGAGGCUCAUUaCCUUUUAAUUUG-------------"
        assert frag.aln_annotation["PP"] == "***...************************************************************999.****************86555555555443............."
        assert frag.aln_annotation["NC"] == "                                                                     v           v                               "
        assert frag.aln_annotation["CS"] == "::::::<<<-<<<<____>>>>->>>,,,,,,,,,,,,,,,,,,,,<<<<<<________>>>>>>,<<<<<<<___>>>>>>>,,,.,,,,,,,,,,,,,,,,,,,,,,,,,"
        assert frag.aln_annotation["similarity"] == "AU+   UCU+:GCCUUUUGGC:+AGAUCAAGUGUAGUAUCUGUUCUU:UCAGU+UAA+A+CUGA:AUG: CC:CA+UG:GG+:CA+U   U+UUAA+UU              "
        # second fragment
        frag = hsp[1]
        assert frag.query_start == 186
        assert frag.query_end == 193
        assert frag.hit_start == 681747
        assert frag.hit_end == 681754
        assert frag.hit_strand == -1
        assert frag.query.seq == "Acccuuu"
        assert frag.hit.seq == "ACAUUUU"
        assert frag.aln_annotation["PP"] == "*******"
        assert frag.aln_annotation["NC"] == "       "
        assert frag.aln_annotation["CS"] == ":::::::"
        assert frag.aln_annotation["similarity"] == "AC +UUU"

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == count

    def test_cmsearch_1q_1m_1h_notextw(self):
        """Test parsing infernal-text, cmsearch, one queries, one hit, one hsp, multiple fragments, notextw (U2_Yeast)"""
        tab_file = get_file("cmsearch_114_U2_Yeast.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, count = next_result(qresults, counter)
        assert len(qresult) == 1
        assert qresult.id == "U2"
        assert qresult.seq_len == 193
        assert qresult.accession == "RF00004"
        assert qresult.description == "U2 spliceosomal RNA"
        assert qresult.program == "cmsearch"
        assert qresult.version == "1.1.4"
        assert qresult.target == "GCA_000146045.2.fasta"
        hit = qresult[0]
        assert len(hit) == 1
        assert hit.id == "ENA|BK006936|BK006936.2"
        assert hit.description == "TPA_inf: Saccharomyces cerevisiae S288C chromosome II, complete sequence."
        assert hit.query_id == "U2"
        hsp = hit[0]
        assert len(hsp) == 2
        assert hsp.model == "cm"
        assert hsp.truncated == "no"
        assert hsp.gc == 0.33
        assert hsp.evalue == 5.9e-20
        assert hsp.bitscore == 98.7
        assert hsp.bias == 0.1
        assert hsp.is_included
        assert hsp.query_start == 1
        assert hsp.query_end == 193
        assert hsp.query_endtype == "[]"
        assert hsp.hit_start == 681747
        assert hsp.hit_end == 681858
        assert hsp.hit_endtype == ".."
        assert hsp.avg_acc == 0.91
        # first fragment
        frag = hsp[0]
        assert frag.query_start == 1
        assert frag.query_end == 112
        assert frag.hit_start == 681763
        assert frag.hit_end == 681858
        assert frag.hit_strand == -1
        assert frag.query.seq == "AUacCUUCucgGCcUUUUgGCuaaGAUCAAGUGUAGUAUCUGUUCUUauCAGUuUAAuAuCUGauAuggcccccAuugggggccaau-uauaUUAaauuaAUUUUUggaacua"
        assert frag.hit.seq == "AUC---UCUUUGCCUUUUGGCUUAGAUCAAGUGUAGUAUCUGUUCUUUUCAGUGUAACAACUGAAAUGA-CCUCAAUGAGGCUCAUUaCCUUUUAAUUUG-------------"
        assert frag.aln_annotation["PP"] == "***...************************************************************999.****************86555555555443............."
        assert frag.aln_annotation["NC"] == "                                                                     v           v                               "
        assert frag.aln_annotation["CS"] == "::::::<<<-<<<<____>>>>->>>,,,,,,,,,,,,,,,,,,,,<<<<<<________>>>>>>,<<<<<<<___>>>>>>>,,,.,,,,,,,,,,,,,,,,,,,,,,,,,"
        assert frag.aln_annotation["similarity"] == "AU+   UCU+:GCCUUUUGGC:+AGAUCAAGUGUAGUAUCUGUUCUU:UCAGU+UAA+A+CUGA:AUG: CC:CA+UG:GG+:CA+U   U+UUAA+UU              "
        # second fragment
        frag = hsp[1]
        assert frag.query_start == 186
        assert frag.query_end == 193
        assert frag.hit_start == 681747
        assert frag.hit_end == 681754
        assert frag.hit_strand == -1
        assert frag.query.seq == "Acccuuu"
        assert frag.hit.seq == "ACAUUUU"
        assert frag.aln_annotation["PP"] == "*******"
        assert frag.aln_annotation["NC"] == "       "
        assert frag.aln_annotation["CS"] == ":::::::"
        assert frag.aln_annotation["similarity"] == "AC +UUU"

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == count

    def test_cmsearch_1q_mm_1h(self):
        """Test parsing infernal-text, cmsearch, one queries, multiple hits, one hsp, multiple fragments (U2_Yeast_full)"""
        tab_file = get_file("cmsearch_114_U2_Yeast_full.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, count = next_result(qresults, counter)
        assert len(qresult) == 5
        assert qresult.id == "U2"
        assert qresult.seq_len == 193
        assert qresult.accession == "RF00004"
        assert qresult.description == "U2 spliceosomal RNA"
        assert qresult.program == "cmsearch"
        assert qresult.version == "1.1.4"
        assert qresult.target == "GCA_000146045.2.fasta"
        # skip first hit (equivalent to test_cmsearch_1q_1m)
        # second hit (3 hsp, reverse strand)
        hit = qresult[1]
        assert len(hit) == 1
        assert hit.id == "ENA|BK006948|BK006948.2"
        assert hit.description == "TPA_inf: Saccharomyces cerevisiae S288C chromosome XV, complete sequence."
        assert hit.query_id == "U2"
        hsp = hit[0]
        assert len(hsp) == 3
        assert hsp.model == "cm"
        assert hsp.truncated == "no"
        assert hsp.gc == 0.39
        assert hsp.evalue == 0.49
        assert hsp.bitscore == 19.8
        assert hsp.bias == 0.0
        assert not hsp.is_included
        assert hsp.query_start == 1
        assert hsp.query_end == 193
        assert hsp.query_endtype == "[]"
        assert hsp.hit_start == 737324
        assert hsp.hit_end == 737498
        assert hsp.hit_endtype == ".."
        assert hsp.avg_acc == 0.96
        # first fragment
        frag = hsp[0]
        assert frag.query_start == 1
        assert frag.query_end == 52
        assert frag.hit_start == 737448
        assert frag.hit_end == 737498
        assert frag.hit_strand == -1
        assert frag.query.seq == "AUacCUUCucgGCcUUUUgGCuaaGAUCAAGUGUAGUAUCUGUUCUUauCAG"
        assert frag.hit.seq == "AUCCCAUAUUUGCCAUC-GGCAUAUAUUAAGUAUAUUAGCAGUUCUAAUUAC"
        assert frag.aln_annotation["PP"] == "**************999.*******************************996"
        assert frag.aln_annotation["NC"] == "                                                    "
        assert frag.aln_annotation["CS"] == "::::::<<<-<<<<____>>>>->>>,,,,,,,,,,,,,,,,,,,,<<<<<<"
        assert frag.aln_annotation["similarity"] == "AU+CC U U+ GCC U  GGC +A AU AAGU UA UA C GUUCU:A::A "
        # second fragment
        frag = hsp[1]
        assert frag.query_start == 60
        assert frag.query_end == 112
        assert frag.hit_start == 737334
        assert frag.hit_end == 737360
        assert frag.hit_strand == -1
        assert frag.query.seq == "CUGauAuggcccccAuugggggccaauuauaUUAaauuaAUUUUUggaacua"
        assert frag.hit.seq == "GUAGUUGGAAGGAUACUAUCCUUUAU--------------------------"
        assert frag.aln_annotation["PP"] == "69999999999999999999999987.........................."
        assert frag.aln_annotation["NC"] == "                                                    "
        assert frag.aln_annotation["CS"] == ">>>>>>,<<<<<<<___>>>>>>>,,,,,,,,,,,,,,,,,,,,,,,,,,,,"
        assert frag.aln_annotation["similarity"] == " U::U:  ::::::A U:::::: A+                          "
        # third fragment
        frag = hsp[2]
        assert frag.query_start == 186
        assert frag.query_end == 193
        assert frag.hit_start == 737324
        assert frag.hit_end == 737331
        assert frag.hit_strand == -1
        assert frag.query.seq == "Acccuuu"
        assert frag.hit.seq == "AUCCCCU"
        assert frag.aln_annotation["PP"] == "*******"
        assert frag.aln_annotation["NC"] == "       "
        assert frag.aln_annotation["CS"] == ":::::::"
        assert frag.aln_annotation["similarity"] == "A CC++U"
        # third hit (2 hsp, forward strand)
        hit = qresult[2]
        assert len(hit) == 1
        assert hit.id == "ENA|BK006947|BK006947.3"
        assert hit.description == "TPA_inf: Saccharomyces cerevisiae S288C chromosome XIV, complete sequence."
        assert hit.query_id == "U2"
        hsp = hit[0]
        assert len(hsp) == 2
        assert hsp.model == "cm"
        assert hsp.truncated == "no"
        assert hsp.gc == 0.39
        assert hsp.evalue == 5.7
        assert hsp.bitscore == 15.3
        assert hsp.bias == 0.0
        assert not hsp.is_included
        assert hsp.query_start == 1
        assert hsp.query_end == 193
        assert hsp.query_endtype == "[]"
        assert hsp.hit_start == 266059
        assert hsp.hit_end == 266208
        assert hsp.hit_endtype == ".."
        assert hsp.avg_acc == 0.91
        # first fragment
        frag = hsp[0]
        assert frag.query_start == 1
        assert frag.query_end == 36
        assert frag.hit_start == 266059
        assert frag.hit_end == 266097
        assert frag.hit_strand == 0
        assert frag.query.seq == "AUacCU-UCu-cgGCcUUUUgGCuaaGAUCAA-GUGUAG"
        assert frag.hit.seq == "AUGUUGaUCUaUCGUCAAUUGACCCAGAUGAUaGUGUAG"
        assert frag.aln_annotation["PP"] == "*****9****999****************9988999987"
        assert frag.aln_annotation["NC"] == "            v          v               "
        assert frag.aln_annotation["CS"] == "::::::.<<<.-<<<<____>>>>->>>,,,,.,,,,,,"
        assert frag.aln_annotation["similarity"] == "AU     UCU + G C  UUG C  AGAU A  GUGUAG"
        # second fragment
        frag = hsp[1]
        assert frag.query_start == 84
        assert frag.query_end == 193
        assert frag.hit_start == 266098
        assert frag.hit_end == 266208
        assert frag.hit_strand == 0
        assert frag.query.seq == "aauuauaUUAaauuaAUUUUUggaacuaGugggggcauuu-uggGCUUGCccau--ugcccccaCacggguugaccuggcaUUGCAcUaccgccagguucagcccAcccuuu"
        assert frag.hit.seq == "-GUUAUAUAGUUUUGAUAUUUUGGCGAAAAGUUGAGAAUAuUGCGCUUGCGUAUauAUUCCAUUUGAGGUGGCACUAGAGCUCGCAUUAU-UACCAGUAGUGGCAGGAUUGC"
        assert frag.aln_annotation["PP"] == ".3377778888888888888888888888888888888888**************99999******************************.99*******************"
        assert frag.aln_annotation["NC"] == "                                v                           v      v  v     v v             v v      v  v       "
        assert frag.aln_annotation["CS"] == ",,,,,,,,,,,,,,,,,,,,,,,,,,,,<<<<<<<<----.<<<<<__>>>>>-..->>>>>>>>,,<<<<<<-<<<<<<___________>>>>>>-->>>>>>:::::::"
        assert frag.aln_annotation["similarity"] == "  UUAUAU   +UU AU UUU G   +A:: : G::A+U+ U :GCUUGC: AU  +::C : ::   G:  :AC: G   U GCA UA+   C :GU+:  :C    +U +"

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == count

    def test_cmsearch_1q_mm_1h_shuffled(self):
        """Test parsing infernal-text, cmsearch, one queries, multiple non-consecutive hits, one hsp, multiple fragments (U2_Yeast_full_shuffled)"""
        tab_file = get_file("cmsearch_114_U2_Yeast_full_shuffled.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, count = next_result(qresults, counter)
        assert len(qresult) == 2
        assert qresult.id == "U2"
        assert qresult.seq_len == 193
        assert qresult.accession == "RF00004"
        assert qresult.description == "U2 spliceosomal RNA"
        assert qresult.program == "cmsearch"
        assert qresult.version == "1.1.4"
        assert qresult.target == "BK006936_7-8.fasta"
        # first hit (3 hsps at rank 1,3 and 4)
        hit = qresult[0]
        assert len(hit) == 3
        assert hit.id == "ENA|BK006936|BK006936.2"
        assert hit.description == ""
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
        assert len(hit) == 3
        assert hit.id == "ENA|BK006948|BK006948.2"
        assert hit.description == ""
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
        assert 1 == count

    def test_cmsearch_1q_mm_mh(self):
        """Test parsing infernal-text, cmsearch, one queries, one hit, multiple hsp, one fragment (5S_Yeast)"""
        tab_file = get_file("cmsearch_114_5S_Yeast.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, count = next_result(qresults, counter)
        assert len(qresult) == 1
        assert qresult.id == "5S_rRNA"
        assert qresult.seq_len == 119
        assert qresult.accession == "RF00001"
        assert qresult.description == "5S ribosomal RNA"
        assert qresult.program == "cmsearch"
        assert qresult.version == "1.1.4"
        assert qresult.target == "GCA_000146045.2.fasta"
        hit = qresult[0]
        assert len(hit) == 6
        assert hit.id == "ENA|BK006945|BK006945.2"
        assert hit.description == "TPA_inf: Saccharomyces cerevisiae S288C chromosome XII, complete sequence."
        assert hit.query_id == "5S_rRNA"
        # first hit
        hsp = hit[0]
        assert len(hsp) == 1
        assert hsp.model == "cm"
        assert hsp.truncated == "no"
        assert hsp.gc == 0.52
        assert hsp.evalue == 1.6e-18
        assert hsp.bitscore == 88.8
        assert hsp.bias == 0.0
        assert hsp.is_included
        assert hsp.query_start == 1
        assert hsp.query_end == 119
        assert hsp.query_endtype == "[]"
        assert hsp.hit_start == 459676
        assert hsp.hit_end == 459796
        assert hsp.hit_endtype == ".."
        assert hsp.avg_acc == 0.99
        frag = hsp[0]
        assert frag.query_start == 1
        assert frag.query_end == 119
        assert frag.hit_start == 459676
        assert frag.hit_end == 459796
        assert frag.hit_strand == 0
        assert frag.query.seq == "gccuGcggcCAUAccagcgcgaAagcACcgGauCCCAUCcGaACuCc-gAAguUAAGcgcgcUugggCcagggUA-GUAcuagGaUGgGuGAcCuCcUGggAAgaccagGugccgCaggcc"
        assert frag.hit.seq == "GGUUGCGGCCAUAUCUACCAGAAAGCACCGUUUCCCGUCCGAUCAACuGUAGUUAAGCUGGUAAGAGCCUGACCGaGUAGUGUAGUGGGUGACCAUACGCGAAACUCAGGUGCUGCAAUCU"
        assert frag.aln_annotation["PP"] == "***********************************************99***********************8756***************************9*****************"
        assert frag.aln_annotation["NC"] == "                                                                               vv                  vv                    "
        assert frag.aln_annotation["CS"] == "(((((((((,,,,<<-<<<<<---<<--<<<<<<______>>-->>>.>-->>---->>>>>-->><<<-<<---.-<-<<-----<<____>>----->>->-->>->>>))))))))):"
        assert frag.aln_annotation["similarity"] == "G::UGC:GCCAUA:C :C::GAAAGCACCG :UCCC+UCCGA C: C G AGUUAAGC::G: +G:GCC G:    GUA  +  +UGGGUGACC+   G  AA  :CAGGUGC:GCA::C+"

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == count

    def test_cmsearch_1q_1m_1h_noali(self):
        """Test parsing infernal-text, cmsearch, one queries, one hit, one hsp, noali (U2_Yeast_noali)"""
        tab_file = get_file("cmsearch_114_U2_Yeast_noali.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, count = next_result(qresults, counter)
        assert len(qresult) == 1
        assert qresult.id == "U2"
        assert qresult.seq_len == 193
        assert qresult.accession == "RF00004"
        assert qresult.description == "U2 spliceosomal RNA"
        assert qresult.program == "cmsearch"
        assert qresult.version == "1.1.4"
        assert qresult.target == "GCA_000146045.2.fasta"
        hit = qresult[0]
        assert len(hit) == 1
        assert hit.id == "ENA|BK006936|BK006936.2"
        assert hit.description == "TPA_inf: Saccharomyces cerevisiae S288C chromosome II,"
        assert hit.query_id == "U2"
        hsp = hit[0]
        assert len(hsp) == 1
        assert hsp.model == "cm"
        assert hsp.truncated == "no"
        assert hsp.gc == 0.33
        assert hsp.evalue == 5.9e-20
        assert hsp.bitscore == 98.7
        assert hsp.bias == 0.1
        assert hsp.is_included
        assert hsp.hit_start == 681747
        assert hsp.hit_end == 681858
        # first fragment
        frag = hsp[0]
        assert frag.hit_start == 681747
        assert frag.hit_end == 681858
        assert frag.hit_strand == -1

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == count

    def test_cmsearch_1q_1m_mh_noali(self):
        """Test parsing infernal-text, cmsearch, one queries, one hit, multiple hsp, noali (5S_Yeast_noali)"""
        tab_file = get_file("cmsearch_114_5S_Yeast_noali.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, count = next_result(qresults, counter)
        assert len(qresult) == 1
        assert qresult.id == "5S_rRNA"
        assert qresult.seq_len == 119
        assert qresult.accession == "RF00001"
        assert qresult.description == "5S ribosomal RNA"
        assert qresult.program == "cmsearch"
        assert qresult.version == "1.1.4"
        assert qresult.target == "GCA_000146045.2.fasta"
        hit = qresult[0]
        assert len(hit) == 6
        assert hit.id == "ENA|BK006945|BK006945.2"
        assert hit.description == "TPA_inf: Saccharomyces cerevisiae S288C chromosome XII"
        assert hit.query_id == "5S_rRNA"
        # first hsp
        hsp = hit[0]
        assert len(hsp) == 1
        assert hsp.model == "cm"
        assert hsp.truncated == "no"
        assert hsp.gc == 0.52
        assert hsp.evalue == 1.6e-18
        assert hsp.bitscore == 88.8
        assert hsp.bias == 0.0
        assert hsp.is_included
        assert hsp.hit_start == 459676
        assert hsp.hit_end == 459796
        frag = hsp[0]
        assert frag.hit_start == 459676
        assert frag.hit_end == 459796
        assert frag.hit_strand == 0
        # last hsp
        hsp = hit[-1]
        assert len(hsp) == 1
        assert hsp.model == "cm"
        assert hsp.truncated == "no"
        assert hsp.gc == 0.53
        assert hsp.evalue == 4.4e-17
        assert hsp.bitscore == 83.2
        assert hsp.bias == 0.0
        assert hsp.is_included
        assert hsp.hit_start == 485697
        assert hsp.hit_end == 485817
        frag = hsp[0]
        assert frag.hit_start == 485697
        assert frag.hit_end == 485817
        assert frag.hit_strand == 0

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == count

    def test_cmsearch_1q_1m_mh_noali_inc(self):
        """Test parsing infernal-text, cmsearch, one queries, one hit, multiple hsp, noali, inclusion threshold (U2_Yeast_full_noali)"""
        tab_file = get_file("cmsearch_114_U2_Yeast_full_noali.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, count = next_result(qresults, counter)
        assert len(qresult) == 5
        assert qresult.id == "U2"
        assert qresult.seq_len == 193
        assert qresult.accession == "RF00004"
        assert qresult.description == "U2 spliceosomal RNA"
        assert qresult.program == "cmsearch"
        assert qresult.version == "1.1.4"
        assert qresult.target == "GCA_000146045.2.fasta"
        # first hit
        hit = qresult[0]
        assert len(hit) == 1
        assert hit.id == "ENA|BK006936|BK006936.2"
        assert hit.description == "TPA_inf: Saccharomyces cerevisiae S288C chromosome II,"
        assert hit.query_id == "U2"
        hsp = hit[0]
        assert len(hsp) == 1
        assert hsp.model == "cm"
        assert hsp.truncated == "no"
        assert hsp.gc == 0.33
        assert hsp.evalue == 5.9e-20
        assert hsp.bitscore == 98.7
        assert hsp.bias == 0.1
        assert hsp.is_included
        assert hsp.hit_start == 681747
        assert hsp.hit_end == 681858
        frag = hsp[0]
        assert frag.hit_start == 681747
        assert frag.hit_end == 681858
        assert frag.hit_strand == -1
        # second hit
        hit = qresult[1]
        assert len(hit) == 1
        assert hit.id == "ENA|BK006948|BK006948.2"
        assert hit.description == "TPA_inf: Saccharomyces cerevisiae S288C chromosome XV,"
        assert hit.query_id == "U2"
        hsp = hit[0]
        assert len(hsp) == 1
        assert hsp.model == "cm"
        assert hsp.truncated == "no"
        assert hsp.gc == 0.39
        assert hsp.evalue == 0.49
        assert hsp.bitscore == 19.8
        assert hsp.bias == 0.0
        assert not hsp.is_included
        assert hsp.hit_start == 737324
        assert hsp.hit_end == 737498
        frag = hsp[0]
        assert frag.hit_start == 737324
        assert frag.hit_end == 737498
        assert frag.hit_strand == -1
        # last hit
        hit = qresult[-1]
        assert len(hit) == 1
        assert hit.id == "ENA|BK006939|BK006939.2"
        assert hit.description == "TPA_inf: Saccharomyces cerevisiae S288C chromosome V,"
        assert hit.query_id == "U2"
        hsp = hit[0]
        assert len(hsp) == 1
        assert hsp.model == "cm"
        assert hsp.truncated == "no"
        assert hsp.gc == 0.41
        assert hsp.evalue == 7.1
        assert hsp.bitscore == 14.9
        assert hsp.bias == 0.0
        assert not hsp.is_included
        assert hsp.hit_start == 190882
        assert hsp.hit_end == 191043
        frag = hsp[0]
        assert frag.hit_start == 190882
        assert frag.hit_end == 191043
        assert frag.hit_strand == 0

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == count

    def test_cmsearch_1q_1m_1h_hmmonly(self):
        """Test parsing infernal-text, cmsearch, one queries, one hit, one hsp, one fragments, hmmonly (U2_Yeast_hmmonly)"""
        tab_file = get_file("cmsearch_114_U2_Yeast_hmmonly.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, count = next_result(qresults, counter)
        assert len(qresult) == 1
        assert qresult.id == "U2"
        assert qresult.seq_len == 193
        assert qresult.accession == "RF00004"
        assert qresult.description == "U2 spliceosomal RNA"
        assert qresult.program == "cmsearch"
        assert qresult.version == "1.1.4"
        assert qresult.target == "GCA_000146045.2.fasta"
        hit = qresult[0]
        assert len(hit) == 1
        assert hit.id == "ENA|BK006936|BK006936.2"
        assert hit.description == "TPA_inf: Saccharomyces cerevisiae S288C chromosome II, complete sequence."
        assert hit.query_id == "U2"
        hsp = hit[0]
        assert len(hsp) == 1
        assert hsp.model == "hmm"
        assert hsp.truncated == "-"
        assert hsp.gc == 0.35
        assert hsp.evalue == 1.5e-19
        assert hsp.bitscore == 73.1
        assert hsp.bias == 2.7
        assert hsp.is_included
        assert hsp.query_start == 7
        assert hsp.query_end == 100
        assert hsp.query_endtype == ".."
        assert hsp.hit_start == 681762
        assert hsp.hit_end == 681855
        assert hsp.hit_endtype == ".."
        assert hsp.avg_acc == 0.74
        # first fragment
        frag = hsp[0]
        assert frag.query_start == 7
        assert frag.query_end == 100
        assert frag.hit_start == 681762
        assert frag.hit_end == 681855
        assert frag.hit_strand == -1
        assert frag.query.seq == "UCucgGCcUUUUGGCUaaGAUCAAGUGUAGUAUCUGUUCUUuuCAGUuUAAuAuCUGauAuugucucuAuugggggccaau-uauaUUAaauuaA"
        assert frag.hit.seq == "UCUUUGCCUUUUGGCUUAGAUCAAGUGUAGUAUCUGUUCUUUUCAGUGUAACAACUGAAAUGACCUCAAU-GAGGCUCAUUaCCUUUUAAUUUGU"
        assert frag.aln_annotation["PP"] == "999*********************************************************7777666666.777776666544444444433221"
        assert frag.aln_annotation["CS"] == "<<<-<<<<____>>>>->>>,,,,,,,,,,,,,,,,,,,,<<<<<<________>>>>>>,<<<<<<<___>>>>>>>,,,.,,,,,,,,,,,,,"
        assert frag.aln_annotation["similarity"] == "UCu+ GCcUUUUGGCU+aGAUCAAGUGUAGUAUCUGUUCUUuuCAGU+UAA+A+CUGa+Au+ +cuc Au g+gg  ca+u   u+UUAa+uu  "

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == count

    def test_cmsearch_mq(self):
        """Test parsing infernal-text, cmsearch, multiple queries"""
        tab_file = get_file("cmsearch_114_IRES_5S_U2_Yeast.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        # First qresult (empty)
        qresult, count = next_result(qresults, counter)
        assert len(qresult) == 0
        # Second qresult (5S, multiple hits)
        qresult, count = next_result(qresults, counter)
        assert len(qresult) == 3
        assert qresult.id == "5S_rRNA"
        assert qresult.seq_len == 119
        assert qresult.accession == "RF00001"
        assert qresult.description == "5S ribosomal RNA"
        assert qresult.program == "cmsearch"
        assert qresult.version == "1.1.4"
        assert qresult.target == "GCA_000146045.2.fasta"
        # first hit
        hit = qresult[0]
        assert len(hit) == 6
        assert hit.id == "ENA|BK006945|BK006945.2"
        assert hit.description == "TPA_inf: Saccharomyces cerevisiae S288C chromosome XII, complete sequence."
        assert hit.query_id == "5S_rRNA"
        hsp = hit[0]
        assert len(hsp) == 1
        assert hsp.model == "cm"
        assert hsp.truncated == "no"
        assert hsp.gc == 0.52
        assert hsp.evalue == 1.6e-18
        assert hsp.bitscore == 88.8
        assert hsp.bias == 0.0
        assert hsp.is_included
        assert hsp.query_start == 1
        assert hsp.query_end == 119
        assert hsp.query_endtype == "[]"
        assert hsp.hit_start == 459676
        assert hsp.hit_end == 459796
        assert hsp.hit_endtype == ".."
        assert hsp.avg_acc == 0.99
        frag = hsp[0]
        assert frag.query_start == 1
        assert frag.query_end == 119
        assert frag.hit_start == 459676
        assert frag.hit_end == 459796
        assert frag.hit_strand == 0
        assert frag.query.seq == "gccuGcggcCAUAccagcgcgaAagcACcgGauCCCAUCcGaACuCc-gAAguUAAGcgcgcUugggCcagggUA-GUAcuagGaUGgGuGAcCuCcUGggAAgaccagGugccgCaggcc"
        assert frag.hit.seq == "GGUUGCGGCCAUAUCUACCAGAAAGCACCGUUUCCCGUCCGAUCAACuGUAGUUAAGCUGGUAAGAGCCUGACCGaGUAGUGUAGUGGGUGACCAUACGCGAAACUCAGGUGCUGCAAUCU"
        assert frag.aln_annotation["PP"] == "***********************************************99***********************8756***************************9*****************"
        assert frag.aln_annotation["NC"] == "                                                                               vv                  vv                    "
        assert frag.aln_annotation["CS"] == "(((((((((,,,,<<-<<<<<---<<--<<<<<<______>>-->>>.>-->>---->>>>>-->><<<-<<---.-<-<<-----<<____>>----->>->-->>->>>))))))))):"
        assert frag.aln_annotation["similarity"] == "G::UGC:GCCAUA:C :C::GAAAGCACCG :UCCC+UCCGA C: C G AGUUAAGC::G: +G:GCC G:    GUA  +  +UGGGUGACC+   G  AA  :CAGGUGC:GCA::C+"
        # last hit
        hsp = hit[-1]
        hit = qresult[-1]
        assert len(hit) == 1
        assert hit.id == "ENA|BK006947|BK006947.3"
        assert hit.description == "TPA_inf: Saccharomyces cerevisiae S288C chromosome XIV, complete sequence."
        assert hit.query_id == "5S_rRNA"
        hsp = hit[0]
        assert len(hsp) == 1
        assert hsp.model == "cm"
        assert hsp.truncated == "no"
        assert hsp.gc == 0.41
        assert hsp.evalue == 6.6
        assert hsp.bitscore == 16.7
        assert hsp.bias == 0.3
        assert not hsp.is_included
        assert hsp.query_start == 1
        assert hsp.query_end == 119
        assert hsp.query_endtype == "[]"
        assert hsp.hit_start == 6968
        assert hsp.hit_end == 7085
        assert hsp.hit_endtype == ".."
        assert hsp.avg_acc == 0.91
        frag = hsp[0]
        assert frag.query_start == 1
        assert frag.query_end == 119
        assert frag.hit_start == 6968
        assert frag.hit_end == 7085
        assert frag.hit_strand == -1
        assert frag.query.seq == "gccuGcggcCAUAccagc-gcg-aAagcACcgGa-uCCCAUCcGaACuCcgAAguUAAGcgcgcUugggCcagggUAGUAcuagGaUGgGuGAcCuCcUGggAAgaccagGu-gccgCaggcc"
        assert frag.hit.seq == "GAGAUGGUAUAUACUGUAgCAUcCGUGUACGUAUgACCGAUCAGA--AUACAAGUGAAGGUGAGUAUGGCAUGUG--GUAGUGGGAUUAGAG-UGGUAGGGUAAGUAUAUGUgUAUUAUUUAC"
        assert frag.aln_annotation["PP"] == "**************976325541459999****989999999999..89**********9999999*********..**********99988.689999************************"
        assert frag.aln_annotation["NC"] == "v             v  v                 v        v                 v   v                                                      v "
        assert frag.aln_annotation["CS"] == "(((((((((,,,,<<-<<.<<<.---<<--<<<<.<<______>>-->>>>-->>---->>>>>-->><<<-<<----<-<<-----<<____>>----->>->-->>->>>.))))))))):"
        assert frag.aln_annotation["similarity"] == " : :: ::: AUAC +   ::     G:AC::::  CC AUC+G   ::::AA:U AAG ::  U+ GGC:  :G  GUA U+GGAU :G G :     GG AAG+: A:GU ::: :: : C"
        # third qresult (U2, multiple hits)
        qresult, count = next_result(qresults, counter)
        assert len(qresult) == 5

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 3 == count


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
