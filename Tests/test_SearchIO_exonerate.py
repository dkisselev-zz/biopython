# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO ExonerateIO parsers."""

import os
import unittest
import pytest

from Bio.SearchIO import parse
from Bio.SearchIO import read

# test case files are in the Blast directory
TEST_DIR = "Exonerate"


def get_file(filename):
    """Return the path of a test file."""
    return os.path.join(TEST_DIR, filename)


class ExonerateSpcCases(unittest.TestCase):
    coord = ("start", "end")
    coords = ("inter_ranges",)
    stype = ("hit_", "query_")

    def check_vulgar_text(self, vulgar, text):
        """Compare coordinate parsing for vulgar and text formats."""
        vfile = get_file(vulgar)
        tfile = get_file(text)

        vqres = read(vfile, "exonerate-vulgar")
        tqres = read(tfile, "exonerate-text")

        # compare coordinates of vulgar and text formats
        # should be the same since the files are results of the same query
        # vs db search
        for vhit, thit in zip(vqres, tqres):
            for vhsp, thsp in zip(vhit.hsps, thit.hsps):
                assert vhsp.query_start == thsp.query_start
                assert vhsp.hit_start == thsp.hit_start
                assert vhsp.query_end == thsp.query_end
                assert vhsp.hit_end == thsp.hit_end
                assert vhsp.query_inter_ranges == thsp.query_inter_ranges
                assert vhsp.hit_inter_ranges == thsp.hit_inter_ranges
                assert vhsp.query_split_codons == thsp.query_split_codons
                assert vhsp.hit_split_codons == thsp.hit_split_codons
                assert vhsp.query_frame_all == thsp.query_frame_all
                assert vhsp.hit_frame_all == thsp.hit_frame_all

    def test_vulgar_text_similar_g2g(self):
        """Compares vulgar-text coordinate parsing for the genome2genome model."""
        self.check_vulgar_text("exn_22_o_vulgar.exn", "exn_22_m_genome2genome.exn")

    def test_vulgar_text_similar_c2c(self):
        """Compares vulgar-text coordinate parsing for the coding2coding model."""
        self.check_vulgar_text(
            "exn_22_o_vulgar_fshifts.exn", "exn_22_m_coding2coding_fshifts.exn"
        )

    def test_vulgar_text_similar_p2d(self):
        """Compares vulgar-text coordinate parsing for the protein2dna model."""
        self.check_vulgar_text(
            "exn_22_o_vulgar_fshifts2.exn", "exn_22_m_protein2dna_fshifts.exn"
        )


class ExonerateTextCases(unittest.TestCase):
    fmt = "exonerate-text"

    def test_exn_22_m_affine_local(self):
        """Test parsing exonerate output (exn_22_m_affine_local.exn)."""
        exn_file = get_file("exn_22_m_affine_local.exn")
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            assert qresult.id == hit.query_id
            for hsp in hit:
                assert hit.id == hsp.hit_id
                assert qresult.id == hsp.query_id

        assert "gi|296143771|ref|NM_001180731.1|" == qresult.id
        assert "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds" == qresult.description
        assert "exonerate" == qresult.program
        assert "affine:local:dna2dna" == qresult.model
        assert 3 == len(qresult)
        # first hit
        hit = qresult[0]
        assert "gi|330443520|ref|NC_001136.10|" == hit.id
        assert "Saccharomyces cerevisiae S288c chromosome IV, complete sequence" == hit.description
        assert 1 == len(hit.hsps)
        # first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 6150 == hsp.score
        assert 1 == hsp[0].query_strand
        assert -1 == hsp[0].hit_strand
        assert 0 == hsp.query_start
        assert 1318045 == hsp.hit_start
        assert 1230 == hsp.query_end
        assert 1319275 == hsp.hit_end
        assert [(0, 1230)] == hsp.query_range_all
        assert [(1318045, 1319275)] == hsp.hit_range_all
        assert [] == hsp.query_inter_ranges
        assert [] == hsp.hit_inter_ranges
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert 1 == len(hsp.query_all)
        assert 1 == len(hsp.hit_all)
        assert "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG" == hsp.query_all[0].seq[:40]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][:40]
        assert "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG" == hsp.hit_all[0].seq[:40]
        assert "TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG" == hsp.query_all[0].seq[-40:]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][-40:]
        assert "TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG" == hsp.hit_all[0].seq[-40:]

        # second hit
        hit = qresult[1]
        assert "gi|330443688|ref|NC_001145.3|" == hit.id
        assert "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence" == hit.description
        assert 1 == len(hit.hsps)
        # second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 359 == hsp.score
        assert 1 == hsp[0].query_strand
        assert 1 == hsp[0].hit_strand
        assert 83 == hsp.query_start
        assert 253990 == hsp.hit_start
        assert 552 == hsp.query_end
        assert 254474 == hsp.hit_end
        assert [(83, 552)] == hsp.query_range_all
        assert [(253990, 254474)] == hsp.hit_range_all
        assert [] == hsp.query_inter_ranges
        assert [] == hsp.hit_inter_ranges
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert 1 == len(hsp.query_all)
        assert 1 == len(hsp.hit_all)
        assert "ACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAA-AGGAT" == hsp.query_all[0].seq[:40]
        assert "||| |||| |||||   ||| | |  ||| |||  | |||" == hsp[0].aln_annotation["similarity"][:40]
        assert "ACCGAAGAAGAAGGGTAGCAAAACTAGCAAAAAGCAAGAT" == hsp.hit_all[0].seq[:40]
        assert "AAGTTATGTGGAACA--TAGGCTCATGGAACGCTCCCAGT" == hsp.query_all[0].seq[-40:]
        assert "| ||   | | ||||  ||   |||   || | ||| |||" == hsp[0].aln_annotation["similarity"][-40:]
        assert "ATGT--GGGGAAACAATTACCTTCACCAAATGATCCAAGT" == hsp.hit_all[0].seq[-40:]

        # third hit
        hit = qresult[2]
        assert "gi|330443715|ref|NC_001146.8|" == hit.id
        assert "Saccharomyces cerevisiae S288c chromosome XIV, complete sequence" == hit.description
        assert 1 == len(hit.hsps)
        # third hit, first hsp
        hsp = qresult[2].hsps[0]
        assert 219 == hsp.score
        assert 1 == hsp[0].query_strand
        assert 1 == hsp[0].hit_strand
        assert 60 == hsp.query_start
        assert 454073 == hsp.hit_start
        assert 517 == hsp.query_end
        assert 454531 == hsp.hit_end
        assert [(60, 517)] == hsp.query_range_all
        assert [(454073, 454531)] == hsp.hit_range_all
        assert [] == hsp.query_inter_ranges
        assert [] == hsp.hit_inter_ranges
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert 1 == len(hsp.query_all)
        assert 1 == len(hsp.hit_all)
        assert "ATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTG-" == hsp.query_all[0].seq[:40]
        assert "||| || || || | |||||   |   |||| ||    | " == hsp[0].aln_annotation["similarity"][:40]
        assert "ATGATGATATATTA-GATGGGG-ATG-AAGATGAGCCAGA" == hsp.hit_all[0].seq[:40]
        assert "G-TATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAG" == hsp.query_all[0].seq[-40:]
        assert "| |  |||||   | | | | |   | ||| | |||||||" == hsp[0].aln_annotation["similarity"][-40:]
        assert "GTTGAAGAAGCCAAACGGAAGAAAGACGAGGA-GAGAAAG" == hsp.hit_all[0].seq[-40:]

    def test_exn_22_m_cdna2genome(self):
        """Test parsing exonerate output (exn_22_m_cdna2genome.exn)."""
        exn_file = get_file("exn_22_m_cdna2genome.exn")
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            assert qresult.id == hit.query_id
            for hsp in hit:
                assert hit.id == hsp.hit_id
                assert qresult.id == hsp.query_id

        assert "gi|296143771|ref|NM_001180731.1|" == qresult.id
        assert "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds" == qresult.description
        assert "exonerate" == qresult.program
        assert "cdna2genome" == qresult.model
        assert 2 == len(qresult)
        # first hit
        hit = qresult[0]
        assert "gi|330443520|ref|NC_001136.10|" == hit.id
        assert "Saccharomyces cerevisiae S288c chromosome IV, complete sequence" == hit.description
        assert 2 == len(hit.hsps)
        # first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 6146 == hsp.score
        assert 1 == hsp[0].query_strand
        assert -1 == hsp[0].hit_strand
        assert 0 == hsp.query_start
        assert 1318045 == hsp.hit_start
        assert 1230 == hsp.query_end
        assert 1319275 == hsp.hit_end
        assert [(0, 1230)] == hsp.query_range_all
        assert [(1318045, 1319275)] == hsp.hit_range_all
        assert [] == hsp.query_inter_ranges
        assert [] == hsp.hit_inter_ranges
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert 1 == len(hsp.query_all)
        assert 1 == len(hsp.hit_all)
        assert "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG" == hsp.query_all[0].seq[:40]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][:40]
        assert "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG" == hsp.hit_all[0].seq[:40]
        assert "TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG" == hsp.query_all[0].seq[-40:]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][-40:]
        assert "TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG" == hsp.hit_all[0].seq[-40:]

        # first hit, second hsp
        hsp = qresult[0].hsps[1]
        assert 6146 == hsp.score
        assert -1 == hsp[0].query_strand
        assert 1 == hsp[0].hit_strand
        assert 0 == hsp.query_start
        assert 1318045 == hsp.hit_start
        assert 1230 == hsp.query_end
        assert 1319275 == hsp.hit_end
        assert [(0, 1230)] == hsp.query_range_all
        assert [(1318045, 1319275)] == hsp.hit_range_all
        assert [] == hsp.query_inter_ranges
        assert [] == hsp.hit_inter_ranges
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert 1 == len(hsp.query_all)
        assert 1 == len(hsp.hit_all)
        assert "CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGA" == hsp.query_all[0].seq[:40]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][:40]
        assert "CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGA" == hsp.hit_all[0].seq[:40]
        assert "CTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCAT" == hsp.query_all[0].seq[-40:]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][-40:]
        assert "CTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCAT" == hsp.hit_all[0].seq[-40:]

        # second hit
        hit = qresult[1]
        assert "gi|330443688|ref|NC_001145.3|" == hit.id
        assert "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence" == hit.description
        assert 1 == len(hit.hsps)
        # second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 518 == hsp.score
        assert 1 == hsp[0].query_strand
        assert 1 == hsp[0].hit_strand
        assert 0 == hsp.query_start
        assert 85010 == hsp.hit_start
        assert 516 == hsp.query_end
        assert 667216 == hsp.hit_end
        assert [(0, 65), (65, 225), (225, 320), (320, 346), (346, 516)] == hsp.query_range_all
        assert [
                (85010, 85066),
                (253974, 254135),
                (350959, 351052),
                (473170, 473201),
                (667040, 667216),
            ] == hsp.hit_range_all
        assert [(65, 65), (225, 225), (320, 320), (346, 346)] == hsp.query_inter_ranges
        assert [(85066, 253974), (254135, 350959), (351052, 473170), (473201, 667040)] == hsp.hit_inter_ranges
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert 5 == len(hsp.query_all)
        assert 5 == len(hsp.hit_all)
        # first block
        assert "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG" == hsp.query_all[0].seq[:40]
        assert "||||  ||  | ||||   | ||||||  |||| | | | " == hsp[0].aln_annotation["similarity"][:40]
        assert "ATGGTGAACCT-CTTCAAGACGGTCAG--AATA-A-TCAA" == hsp.hit_all[0].seq[:40]
        assert "AGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTT" == hsp.query_all[0].seq[-40:]
        assert "||  |||| | | | ||||    ||||||||||||| | |" == hsp[0].aln_annotation["similarity"][-40:]
        assert "AG--AATA-A-TCAACAGG----ATGAAGAAGCAAAAGAT" == hsp.hit_all[0].seq[-40:]
        # last block
        assert "TATTAGCCTTCC--TCGATGATCTGCA--A-GAACAACAG" == hsp.query_all[-1].seq[:40]
        assert "|  |||| || |  ||||| | || ||  | ||| | |  " == hsp[-1].aln_annotation["similarity"][:40]
        assert "TCATAGCGTTACGTTCGAT-ACCTTCACTACGAAGATCCA" == hsp.hit_all[-1].seq[:40]
        assert "AAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAA" == hsp.query_all[-1].seq[-40:]
        assert "   |||||||||||||     ||  ||| | ||  | |||" == hsp[-1].aln_annotation["similarity"][-40:]
        assert "TTCTATAGAAGTACAGTTATTCAAACAAAAAAAAAAAAAA" == hsp.hit_all[-1].seq[-40:]

    def test_exn_22_m_coding2coding(self):
        """Test parsing exonerate output (exn_22_m_coding2coding.exn)."""
        exn_file = get_file("exn_22_m_coding2coding.exn")
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            assert qresult.id == hit.query_id
            for hsp in hit:
                assert hit.id == hsp.hit_id
                assert qresult.id == hsp.query_id

        assert "gi|296143771|ref|NM_001180731.1|" == qresult.id
        assert "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds" == qresult.description
        assert "exonerate" == qresult.program
        assert "coding2coding" == qresult.model
        assert 2 == len(qresult)
        # first hit
        hit = qresult[0]
        assert "gi|330443520|ref|NC_001136.10|" == hit.id
        assert "Saccharomyces cerevisiae S288c chromosome IV, complete sequence" == hit.description
        assert 2 == len(hit.hsps)
        # first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 2151 == hsp.score
        assert -1 == hsp[0].query_strand
        assert 1 == hsp[0].hit_strand
        assert 1 == hsp.query_start
        assert 1318047 == hsp.hit_start
        assert 1228 == hsp.query_end
        assert 1319274 == hsp.hit_end
        assert [(1, 1228)] == hsp.query_range_all
        assert [(1318047, 1319274)] == hsp.hit_range_all
        assert [] == hsp.query_inter_ranges
        assert [] == hsp.hit_inter_ranges
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert 1 == len(hsp.query_all)
        assert 1 == len(hsp.hit_all)
        assert "ACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGC" == hsp.query_all[0].seq[:40]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][:40]
        assert "ACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGC" == hsp.hit_all[0].seq[:40]
        assert "GCTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCA" == hsp.query_all[0].seq[-40:]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][-40:]
        assert "GCTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCA" == hsp.hit_all[0].seq[-40:]

        # first hit, second hsp
        hsp = qresult[0].hsps[1]
        assert 2106 == hsp.score
        assert 1 == hsp[0].query_strand
        assert -1 == hsp[0].hit_strand
        assert 0 == hsp.query_start
        assert 1318045 == hsp.hit_start
        assert 1230 == hsp.query_end
        assert 1319275 == hsp.hit_end
        assert [(0, 1230)] == hsp.query_range_all
        assert [(1318045, 1319275)] == hsp.hit_range_all
        assert [] == hsp.query_inter_ranges
        assert [] == hsp.hit_inter_ranges
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert 1 == len(hsp.query_all)
        assert 1 == len(hsp.hit_all)
        assert "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG" == hsp.query_all[0].seq[:40]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][:40]
        assert "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG" == hsp.hit_all[0].seq[:40]
        assert "TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG" == hsp.query_all[0].seq[-40:]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][-40:]
        assert "TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG" == hsp.hit_all[0].seq[-40:]

        # second hit
        hit = qresult[1]
        assert "gi|330443688|ref|NC_001145.3|" == hit.id
        assert "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence" == hit.description
        assert 1 == len(hit.hsps)
        # second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 116 == hsp.score
        assert 1 == hsp[0].query_strand
        assert 1 == hsp[0].hit_strand
        assert 1065 == hsp.query_start
        assert 255638 == hsp.hit_start
        assert 1224 == hsp.query_end
        assert 255794 == hsp.hit_end
        assert [(1065, 1224)] == hsp.query_range_all
        assert [(255638, 255794)] == hsp.hit_range_all
        assert [] == hsp.query_inter_ranges
        assert [] == hsp.hit_inter_ranges
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert 1 == len(hsp.query_all)
        assert 1 == len(hsp.hit_all)
        assert "TGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAAT" == hsp.query_all[0].seq[:40]
        assert "||+!  .!.|||   !!:...||+:!::!:!  ||+||||" == hsp[0].aln_annotation["similarity"][:40]
        assert "TGTTCGGAAATTTGGGATAGAATAACAACACATCCGAAAT" == hsp.hit_all[0].seq[:40]
        assert "CAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTC" == hsp.query_all[0].seq[-40:]
        assert "!!!.||+...|||:!:||+   |||+||  !!::!!.:!:" == hsp[0].aln_annotation["similarity"][-40:]
        assert "CAATGCAGAAGACGTTCAATTAGCTTTGAATAAGCATATG" == hsp.hit_all[0].seq[-40:]

    def test_exn_22_m_coding2genome(self):
        """Test parsing exonerate output (exn_22_m_coding2genome.exn)."""
        exn_file = get_file("exn_22_m_coding2genome.exn")
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            assert qresult.id == hit.query_id
            for hsp in hit:
                assert hit.id == hsp.hit_id
                assert qresult.id == hsp.query_id

        assert "gi|296143771|ref|NM_001180731.1|" == qresult.id
        assert "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds" == qresult.description
        assert "exonerate" == qresult.program
        assert "coding2genome" == qresult.model
        assert 2 == len(qresult)
        # first hit
        hit = qresult[0]
        assert "gi|330443520|ref|NC_001136.10|" == hit.id
        assert "Saccharomyces cerevisiae S288c chromosome IV, complete sequence" == hit.description
        assert 2 == len(hit.hsps)
        # first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 2151 == hsp.score
        assert -1 == hsp[0].query_strand
        assert 1 == hsp[0].hit_strand
        assert 1 == hsp.query_start
        assert 1318047 == hsp.hit_start
        assert 1228 == hsp.query_end
        assert 1319274 == hsp.hit_end
        assert [(1, 1228)] == hsp.query_range_all
        assert [(1318047, 1319274)] == hsp.hit_range_all
        assert [] == hsp.query_inter_ranges
        assert [] == hsp.hit_inter_ranges
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert 1 == len(hsp.query_all)
        assert 1 == len(hsp.hit_all)
        assert "ACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGC" == hsp.query_all[0].seq[:40]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][:40]
        assert "ACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGC" == hsp.hit_all[0].seq[:40]
        assert "GCTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCA" == hsp.query_all[0].seq[-40:]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][-40:]
        assert "GCTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCA" == hsp.hit_all[0].seq[-40:]

        # first hit, second hsp
        hsp = qresult[0].hsps[1]
        assert 2106 == hsp.score
        assert 1 == hsp[0].query_strand
        assert -1 == hsp[0].hit_strand
        assert 0 == hsp.query_start
        assert 1318045 == hsp.hit_start
        assert 1230 == hsp.query_end
        assert 1319275 == hsp.hit_end
        assert [(0, 1230)] == hsp.query_range_all
        assert [(1318045, 1319275)] == hsp.hit_range_all
        assert [] == hsp.query_inter_ranges
        assert [] == hsp.hit_inter_ranges
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert 1 == len(hsp.query_all)
        assert 1 == len(hsp.hit_all)
        assert "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG" == hsp.query_all[0].seq[:40]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][:40]
        assert "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG" == hsp.hit_all[0].seq[:40]
        assert "TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG" == hsp.query_all[0].seq[-40:]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][-40:]
        assert "TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG" == hsp.hit_all[0].seq[-40:]

        # second hit
        hit = qresult[1]
        assert "gi|330443688|ref|NC_001145.3|" == hit.id
        assert "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence" == hit.description
        assert 1 == len(hit.hsps)
        # second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 116 == hsp.score
        assert 1 == hsp[0].query_strand
        assert 1 == hsp[0].hit_strand
        assert 1065 == hsp.query_start
        assert 255638 == hsp.hit_start
        assert 1224 == hsp.query_end
        assert 255794 == hsp.hit_end
        assert [(1065, 1224)] == hsp.query_range_all
        assert [(255638, 255794)] == hsp.hit_range_all
        assert [] == hsp.query_inter_ranges
        assert [] == hsp.hit_inter_ranges
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert 1 == len(hsp.query_all)
        assert 1 == len(hsp.hit_all)
        assert "TGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAAT" == hsp.query_all[0].seq[:40]
        assert "||+!  .!.|||   !!:...||+:!::!:!  ||+||||" == hsp[0].aln_annotation["similarity"][:40]
        assert "TGTTCGGAAATTTGGGATAGAATAACAACACATCCGAAAT" == hsp.hit_all[0].seq[:40]
        assert "CAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTC" == hsp.query_all[0].seq[-40:]
        assert "!!!.||+...|||:!:||+   |||+||  !!::!!.:!:" == hsp[0].aln_annotation["similarity"][-40:]
        assert "CAATGCAGAAGACGTTCAATTAGCTTTGAATAAGCATATG" == hsp.hit_all[0].seq[-40:]

    def test_exn_22_m_dna2protein(self):
        """Test parsing exonerate output (exn_22_m_dna2protein.exn)."""
        exn_file = get_file("exn_22_m_dna2protein.exn")
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            assert qresult.id == hit.query_id
            for hsp in hit:
                assert hit.id == hsp.hit_id
                assert qresult.id == hsp.query_id

        assert "dna" == qresult.id
        assert "" == qresult.description
        assert "exonerate" == qresult.program
        assert "ungapped:dna2protein" == qresult.model
        assert 1 == len(qresult)

        hit = qresult[0]
        assert "protein" == hit.id
        assert "" == hit.description
        assert 1 == len(hit.hsps)

        hsp = hit[0]
        assert 105 == hsp.score
        assert 1 == hsp[0].query_strand
        assert 0 == hsp[0].hit_strand
        assert 0 == hsp[0].query_start
        assert 93 == hsp[0].query_end
        assert 313 == hsp[0].hit_start
        assert 344 == hsp[0].hit_end
        assert 1 == len(hsp.query_all)
        assert 1 == len(hsp.hit_all)
        assert "NSPFXKGPLASVQNPVYHKQPLNPAPNAETH" == hsp[0].query.seq[:40]
        assert ["|||", "...", " !!", " !!", "! !"] == hsp[0].aln_annotation["similarity"][:5]
        assert "NQSVPKRPAGSVQNPVYHNQPLNPAPSRDPH" == hsp[0].hit.seq[:40]

    def test_exn_22_m_est2genome(self):
        """Test parsing exonerate output (exn_22_m_est2genome.exn)."""
        exn_file = get_file("exn_22_m_est2genome.exn")
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            assert qresult.id == hit.query_id
            for hsp in hit:
                assert hit.id == hsp.hit_id
                assert qresult.id == hsp.query_id

        assert "gi|296143771|ref|NM_001180731.1|" == qresult.id
        assert "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds" == qresult.description
        assert "exonerate" == qresult.program
        assert "est2genome" == qresult.model
        assert 2 == len(qresult)
        # first hit
        hit = qresult[0]
        assert "gi|330443520|ref|NC_001136.10|" == hit.id
        assert "Saccharomyces cerevisiae S288c chromosome IV, complete sequence" == hit.description
        assert 1 == len(hit.hsps)
        # first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 6150 == hsp.score
        assert 1 == hsp[0].query_strand
        assert -1 == hsp[0].hit_strand
        assert 0 == hsp.query_start
        assert 1318045 == hsp.hit_start
        assert 1230 == hsp.query_end
        assert 1319275 == hsp.hit_end
        assert [(0, 1230)] == hsp.query_range_all
        assert [(1318045, 1319275)] == hsp.hit_range_all
        assert [] == hsp.query_inter_ranges
        assert [] == hsp.hit_inter_ranges
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert 1 == len(hsp.query_all)
        assert 1 == len(hsp.hit_all)
        assert "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG" == hsp.query_all[0].seq[:40]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][:40]
        assert "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG" == hsp.hit_all[0].seq[:40]
        assert "TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG" == hsp.query_all[0].seq[-40:]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][-40:]
        assert "TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG" == hsp.hit_all[0].seq[-40:]

        # second hit
        hit = qresult[1]
        assert "gi|330443688|ref|NC_001145.3|" == hit.id
        assert "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence" == hit.description
        assert 2 == len(hit.hsps)
        # second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 439 == hsp.score
        assert 1 == hsp[0].query_strand
        assert 1 == hsp[0].hit_strand
        assert 0 == hsp.query_start
        assert 85010 == hsp.hit_start
        assert 346 == hsp.query_end
        assert 473201 == hsp.hit_end
        assert [(0, 65), (65, 225), (225, 320), (320, 346)] == hsp.query_range_all
        assert [(85010, 85066), (253974, 254135), (350959, 351052), (473170, 473201)] == hsp.hit_range_all
        assert [(65, 65), (225, 225), (320, 320)] == hsp.query_inter_ranges
        assert [(85066, 253974), (254135, 350959), (351052, 473170)] == hsp.hit_inter_ranges
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert 4 == len(hsp.query_all)
        assert 4 == len(hsp.hit_all)
        # first block
        assert "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG" == hsp.query_all[0].seq[:40]
        assert "||||  ||  | ||||   | ||||||  |||| | | | " == hsp[0].aln_annotation["similarity"][:40]
        assert "ATGGTGAACCT-CTTCAAGACGGTCAG--AATA-A-TCAA" == hsp.hit_all[0].seq[:40]
        assert "AGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTT" == hsp.query_all[0].seq[-40:]
        assert "||  |||| | | | ||||    ||||||||||||| | |" == hsp[0].aln_annotation["similarity"][-40:]
        assert "AG--AATA-A-TCAACAGG----ATGAAGAAGCAAAAGAT" == hsp.hit_all[0].seq[-40:]
        # last block
        assert "AGCTAAGAATTCTGATGATG-----AAAGAA" == hsp.query_all[-1].seq
        assert "|   |||||||||||| |||     ||||||" == hsp[-1].aln_annotation["similarity"]
        assert "ATGGAAGAATTCTGATAATGCTGTAAAAGAA" == hsp.hit_all[-1].seq

        # second hit, second hsp
        hsp = qresult[1].hsps[1]
        assert 263 == hsp.score
        assert 1 == hsp[0].query_strand
        assert -1 == hsp[0].hit_strand
        assert 25 == hsp.query_start
        assert 11338 == hsp.hit_start
        assert 406 == hsp.query_end
        assert 130198 == hsp.hit_end
        assert [(25, 183), (183, 252), (252, 406)] == hsp.query_range_all
        assert [(130038, 130198), (120612, 120681), (11338, 11487)] == hsp.hit_range_all
        assert [(183, 183), (252, 252)] == hsp.query_inter_ranges
        assert [(120681, 130038), (11487, 120612)] == hsp.hit_inter_ranges
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert 3 == len(hsp.query_all)
        assert 3 == len(hsp.hit_all)
        # first block
        assert "AGCAAATATATTTA-GCAGGTGACATGAAGAAGCAAATGT" == hsp.query_all[0].seq[:40]
        assert "| |||| |||   | ||||   | | || |||| | |  |" == hsp[0].aln_annotation["similarity"][:40]
        assert "ACCAAAGATAACAAGGCAG--AAAAAGAGGAAGAAGAAAT" == hsp.hit_all[0].seq[:40]
        assert "AG-GACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGA" == hsp.query_all[0].seq[-40:]
        assert "|| |||  ||||||  ||   |||  || ||   ||| ||" == hsp[0].aln_annotation["similarity"][-40:]
        assert "AGTGAC--CCCAGAGGAGCCAAGCAAAAAGA---TTCGGA" == hsp.hit_all[0].seq[-40:]
        # last block
        assert "AATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAA" == hsp.query_all[-1].seq[:40]
        assert "|||||||  | ||| |    |||| | |  | | ||    " == hsp[-1].aln_annotation["similarity"][:40]
        assert "AATAAGAGCAACACAG----TTTA-TCTTATATGTA----" == hsp.hit_all[-1].seq[:40]
        assert "CTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAA" == hsp.query_all[-1].seq[-40:]
        assert "|  | | || |  | || ||  ||||||||  ||  ||||" == hsp[-1].aln_annotation["similarity"][-40:]
        assert "CCACTAAAAAATTATAAGAGCCAAAACGAAGTAGATACAA" == hsp.hit_all[-1].seq[-40:]

    def test_exn_22_m_genome2genome(self):
        """Test parsing exonerate output (exn_22_m_genome2genome.exn)."""
        exn_file = get_file("exn_22_m_genome2genome.exn")
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            assert qresult.id == hit.query_id
            for hsp in hit:
                assert hit.id == hsp.hit_id
                assert qresult.id == hsp.query_id

        assert "sacCer3_dna" == qresult.id
        assert "range=chrIV:1319469-1319997 5'pad=0 3'pad=0 strand=+ repeatMasking=none" == qresult.description
        assert "exonerate" == qresult.program
        assert "genome2genome" == qresult.model
        assert 3 == len(qresult)
        # first hit
        hit = qresult[0]
        assert "gi|330443520|ref|NC_001136.10|" == hit.id
        assert "Saccharomyces cerevisiae S288c chromosome IV, complete sequence" == hit.description
        assert 2 == len(hit.hsps)
        # first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 2641 == hsp.score
        assert -1 == hsp[0].query_strand
        assert -1 == hsp[0].hit_strand
        assert 0 == hsp.query_start
        assert 1319468 == hsp.hit_start
        assert 529 == hsp.query_end
        assert 1319997 == hsp.hit_end
        assert [(0, 529)] == hsp.query_range_all
        assert [(1319468, 1319997)] == hsp.hit_range_all
        assert [] == hsp.query_inter_ranges
        assert [] == hsp.hit_inter_ranges
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert 1 == len(hsp.query_all)
        assert 1 == len(hsp.hit_all)
        assert "ATCCCTTATCTCTTTATCTTGTTGCCTGGTTCTCTTTTCC" == hsp.query_all[0].seq[:40]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][:40]
        assert "ATCCCTTATCTCTTTATCTTGTTGCCTGGTTCTCTTTTCC" == hsp.hit_all[0].seq[:40]
        assert "ACGGCAATACCTGGCATGTGATTGTCGGAAAGAACTTTGG" == hsp.query_all[0].seq[-40:]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][-40:]
        assert "ACGGCAATACCTGGCATGTGATTGTCGGAAAGAACTTTGG" == hsp.hit_all[0].seq[-40:]

        # first hit, second hsp
        hsp = qresult[0].hsps[1]
        assert 2641 == hsp.score
        assert 1 == hsp[0].query_strand
        assert 1 == hsp[0].hit_strand
        assert 0 == hsp.query_start
        assert 1319468 == hsp.hit_start
        assert 529 == hsp.query_end
        assert 1319997 == hsp.hit_end
        assert [(0, 529)] == hsp.query_range_all
        assert [(1319468, 1319997)] == hsp.hit_range_all
        assert [] == hsp.query_inter_ranges
        assert [] == hsp.hit_inter_ranges
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert 1 == len(hsp.query_all)
        assert 1 == len(hsp.hit_all)
        assert "CCAAAGTTCTTTCCGACAATCACATGCCAGGTATTGCCGT" == hsp.query_all[0].seq[:40]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][:40]
        assert "CCAAAGTTCTTTCCGACAATCACATGCCAGGTATTGCCGT" == hsp.hit_all[0].seq[:40]
        assert "GGAAAAGAGAACCAGGCAACAAGATAAAGAGATAAGGGAT" == hsp.query_all[0].seq[-40:]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][-40:]
        assert "GGAAAAGAGAACCAGGCAACAAGATAAAGAGATAAGGGAT" == hsp.hit_all[0].seq[-40:]

        # second hit
        hit = qresult[1]
        assert "gi|330443489|ref|NC_001135.5|" == hit.id
        assert "Saccharomyces cerevisiae S288c chromosome III, complete sequence" == hit.description
        assert 1 == len(hit.hsps)
        # second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 267 == hsp.score
        assert -1 == hsp[0].query_strand
        assert 1 == hsp[0].hit_strand
        assert 162 == hsp.query_start
        assert 23668 == hsp.hit_start
        assert 491 == hsp.query_end
        assert 115569 == hsp.hit_end
        assert [(462, 491), (413, 462), (378, 413), (302, 378), (162, 302)] == hsp.query_range_all
        assert [
                (23668, 23697),
                (32680, 32732),
                (42287, 42325),
                (97748, 97821),
                (115419, 115569),
            ] == hsp.hit_range_all
        assert [(462, 462), (413, 413), (378, 378), (302, 302)] == hsp.query_inter_ranges
        assert [(23697, 32680), (32732, 42287), (42325, 97748), (97821, 115419)] == hsp.hit_inter_ranges
        assert [(378, 379), (376, 378)] == hsp.query_split_codons
        assert [(42324, 42325), (97748, 97750)] == hsp.hit_split_codons
        assert 5 == len(hsp.query_all)
        assert 5 == len(hsp.hit_all)
        # first block
        assert "CCCTTTAAATGGAGATTACAAACTAGCGA" == hsp.query_all[0].seq
        assert "||  | ||| | |||  ||||| |  | |" == hsp[0].aln_annotation["similarity"]
        assert "CCGCTGAAAGGAAGAGAACAAAGTTACAA" == hsp.hit_all[0].seq
        # last block
        assert "TTTTCTTTACTAAC-TCGAGGAAGAGTGAGGTTTTCTTCC" == hsp.query_all[-1].seq[:40]
        assert "| ||    || | | |  |||||| |||| | | |  |||" == hsp[-1].aln_annotation["similarity"][:40]
        assert "TCTTGAAGACCAGCATGTAGGAAG-GTGATGATATGCTCC" == hsp.hit_all[-1].seq[:40]
        assert "TTTGTGTGTGTACATTTGAATATATATATTTAC-TAACAA" == hsp.query_all[-1].seq[-40:]
        assert " |||  ||| |   |||||||||||||   | | ||||||" == hsp[-1].aln_annotation["similarity"][-40:]
        assert "ATTGATTGTTTTGTTTTGAATATATATTGATGCTTAACAA" == hsp.hit_all[-1].seq[-40:]

        # third hit
        hit = qresult[2]
        assert "gi|330443667|ref|NC_001143.9|" == hit.id
        assert "Saccharomyces cerevisiae S288c chromosome XI, complete sequence" == hit.description
        assert 1 == len(hit.hsps)
        # third hit, first hsp
        hsp = qresult[2].hsps[0]
        assert 267 == hsp.score
        assert -1 == hsp[0].query_strand
        assert -1 == hsp[0].hit_strand
        assert 78 == hsp.query_start
        assert 71883 == hsp.hit_start
        assert 529 == hsp.query_end
        assert 641760 == hsp.hit_end
        assert [(449, 529), (319, 388), (198, 284), (161, 198), (78, 114)] == hsp.query_range_all
        assert [
                (641682, 641760),
                (487327, 487387),
                (386123, 386207),
                (208639, 208677),
                (71883, 71917),
            ] == hsp.hit_range_all
        assert [(388, 449), (284, 319), (198, 198), (114, 161)] == hsp.query_inter_ranges
        assert [(487387, 641682), (386207, 487327), (208677, 386123), (71917, 208639)] == hsp.hit_inter_ranges
        assert [(198, 200), (197, 198)] == hsp.query_split_codons
        assert [(386123, 386125), (208676, 208677)] == hsp.hit_split_codons
        assert 5 == len(hsp.query_all)
        assert 5 == len(hsp.hit_all)
        # first block
        assert "ATCCCTTATCTCTTTATCTTGTTGCCTGGTTCTCTTTTCC" == hsp.query_all[0].seq[:40]
        assert "||||||||||||||       |||  |||||   ||||  " == hsp[0].aln_annotation["similarity"][:40]
        assert "ATCCCTTATCTCTTCTAAAGATTGTGTGGTT---TTTT--" == hsp.hit_all[0].seq[:40]
        assert "AAATGGAGATTACAA---ACTAGCGAA-ACTGCAGAAAAG" == hsp.query_all[0].seq[-40:]
        assert "  ||     || |||    || || ||  || || | |||" == hsp[0].aln_annotation["similarity"][-40:]
        assert "GCATATTTTTTCCAACCTTCTTGCCAATTCTTCA-ACAAG" == hsp.hit_all[0].seq[-40:]
        # last block
        assert "TAAAGATGCTCTGGACAAGTACCAGTTGGAAAGAGA" == hsp.query_all[-1].seq
        assert " ||||||  |||  || | |  ||||||||||||||" == hsp[-1].aln_annotation["similarity"]
        assert "AAAAGATTTTCT--ACGACTTGCAGTTGGAAAGAGA" == hsp.hit_all[-1].seq

    def test_exn_22_m_ungapped(self):
        """Test parsing exonerate output (exn_22_m_ungapped.exn)."""
        exn_file = get_file("exn_22_m_ungapped.exn")
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            assert qresult.id == hit.query_id
            for hsp in hit:
                assert hit.id == hsp.hit_id
                assert qresult.id == hsp.query_id

        assert "gi|296143771|ref|NM_001180731.1|" == qresult.id
        assert "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds" == qresult.description
        assert "exonerate" == qresult.program
        assert "ungapped:dna2dna" == qresult.model
        assert 2 == len(qresult)
        # first hit
        hit = qresult[0]
        assert "gi|330443520|ref|NC_001136.10|" == hit.id
        assert "Saccharomyces cerevisiae S288c chromosome IV, complete sequence" == hit.description
        assert 1 == len(hit.hsps)
        # first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 6150 == hsp.score
        assert 1 == hsp[0].query_strand
        assert -1 == hsp[0].hit_strand
        assert 0 == hsp.query_start
        assert 1318045 == hsp.hit_start
        assert 1230 == hsp.query_end
        assert 1319275 == hsp.hit_end
        assert [(0, 1230)] == hsp.query_range_all
        assert [(1318045, 1319275)] == hsp.hit_range_all
        assert [] == hsp.query_inter_ranges
        assert [] == hsp.hit_inter_ranges
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert 1 == len(hsp.query_all)
        assert 1 == len(hsp.hit_all)
        assert "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG" == hsp.query_all[0].seq[:40]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][:40]
        assert "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG" == hsp.hit_all[0].seq[:40]
        assert "TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG" == hsp.query_all[0].seq[-40:]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][-40:]
        assert "TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG" == hsp.hit_all[0].seq[-40:]

        # second hit
        hit = qresult[1]
        assert "gi|330443688|ref|NC_001145.3|" == hit.id
        assert "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence" == hit.description
        assert 2 == len(hit.hsps)
        # second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 233 == hsp.score
        assert 1 == hsp[0].query_strand
        assert 1 == hsp[0].hit_strand
        assert 121 == hsp.query_start
        assert 254031 == hsp.hit_start
        assert 236 == hsp.query_end
        assert 254146 == hsp.hit_end
        assert [(121, 236)] == hsp.query_range_all
        assert [(254031, 254146)] == hsp.hit_range_all
        assert [] == hsp.query_inter_ranges
        assert [] == hsp.hit_inter_ranges
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert 1 == len(hsp.query_all)
        assert 1 == len(hsp.hit_all)
        assert "TTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAG" == hsp.query_all[0].seq[:40]
        assert "| ||  ||||| |||||   | |||||||||||| ||| |" == hsp[0].aln_annotation["similarity"][:40]
        assert "TGGATCCTGAAACTAAGCAGAAGAGGACTGCCCAAAATCG" == hsp.hit_all[0].seq[:40]
        assert "CCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGA" == hsp.query_all[0].seq[-40:]
        assert "  || ||||||   ||| |  ||| |||| |     ||||" == hsp[0].aln_annotation["similarity"][-40:]
        assert "GGAAGATGAAGGAATTGGAGAAGAAGGTACAAAGTTTAGA" == hsp.hit_all[0].seq[-40:]

        # second hit, second hsp
        hsp = qresult[1].hsps[1]
        assert 151 == hsp.score
        assert 1 == hsp[0].query_strand
        assert 1 == hsp[0].hit_strand
        assert 1098 == hsp.query_start
        assert 255671 == hsp.hit_start
        assert 1166 == hsp.query_end
        assert 255739 == hsp.hit_end
        assert [(1098, 1166)] == hsp.query_range_all
        assert [(255671, 255739)] == hsp.hit_range_all
        assert [] == hsp.query_inter_ranges
        assert [] == hsp.hit_inter_ranges
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert 1 == len(hsp.query_all)
        assert 1 == len(hsp.hit_all)
        assert "CCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCG" == hsp.query_all[0].seq[:40]
        assert "|| ||||| |||    | ||  | |||| ||||||   ||" == hsp[0].aln_annotation["similarity"][:40]
        assert "CCGAAATACTCAGATATTGATGTCGATGGTTTATGTTCCG" == hsp.hit_all[0].seq[:40]
        assert "ATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGA" == hsp.query_all[0].seq[-40:]
        assert " ||||||   |||  ||||    |||||||||||| ||||" == hsp[0].aln_annotation["similarity"][-40:]
        assert "GTTTATGTTCCGAGCTAATGGCAAAGGCAAAATGTTCAGA" == hsp.hit_all[0].seq[-40:]

    def test_exn_22_m_ungapped_trans(self):
        """Test parsing exonerate output (exn_22_m_ungapped_trans.exn)."""
        exn_file = get_file("exn_22_m_ungapped_trans.exn")
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            assert qresult.id == hit.query_id
            for hsp in hit:
                assert hit.id == hsp.hit_id
                assert qresult.id == hsp.query_id

        assert "gi|296143771|ref|NM_001180731.1|" == qresult.id
        assert "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds" == qresult.description
        assert "exonerate" == qresult.program
        assert "ungapped:codon" == qresult.model
        assert 1 == len(qresult)
        # first hit
        hit = qresult[0]
        assert "gi|330443520|ref|NC_001136.10|" == hit.id
        assert "Saccharomyces cerevisiae S288c chromosome IV, complete sequence" == hit.description
        assert 3 == len(hit.hsps)
        # first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 2151 == hsp.score
        assert -1 == hsp[0].query_strand
        assert 1 == hsp[0].hit_strand
        assert 1 == hsp.query_start
        assert 1318047 == hsp.hit_start
        assert 1228 == hsp.query_end
        assert 1319274 == hsp.hit_end
        assert [(1, 1228)] == hsp.query_range_all
        assert [(1318047, 1319274)] == hsp.hit_range_all
        assert [] == hsp.query_inter_ranges
        assert [] == hsp.hit_inter_ranges
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert 1 == len(hsp.query_all)
        assert 1 == len(hsp.hit_all)
        assert "ACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGC" == hsp.query_all[0].seq[:40]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][:40]
        assert "ACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGC" == hsp.hit_all[0].seq[:40]
        assert "GCTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCA" == hsp.query_all[0].seq[-40:]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][-40:]
        assert "GCTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCA" == hsp.hit_all[0].seq[-40:]

        # first hit, second hsp
        hsp = qresult[0].hsps[1]
        assert 2106 == hsp.score
        assert 1 == hsp[0].query_strand
        assert -1 == hsp[0].hit_strand
        assert 0 == hsp.query_start
        assert 1318045 == hsp.hit_start
        assert 1230 == hsp.query_end
        assert 1319275 == hsp.hit_end
        assert [(0, 1230)] == hsp.query_range_all
        assert [(1318045, 1319275)] == hsp.hit_range_all
        assert [] == hsp.query_inter_ranges
        assert [] == hsp.hit_inter_ranges
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert 1 == len(hsp.query_all)
        assert 1 == len(hsp.hit_all)
        assert "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG" == hsp.query_all[0].seq[:40]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][:40]
        assert "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG" == hsp.hit_all[0].seq[:40]
        assert "TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG" == hsp.query_all[0].seq[-40:]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][-40:]
        assert "TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG" == hsp.hit_all[0].seq[-40:]

        # first hit, third hsp
        hsp = qresult[0].hsps[2]
        assert 2072 == hsp.score
        assert -1 == hsp[0].query_strand
        assert 1 == hsp[0].hit_strand
        assert 0 == hsp.query_start
        assert 1318045 == hsp.hit_start
        assert 1230 == hsp.query_end
        assert 1319275 == hsp.hit_end
        assert [(0, 1230)] == hsp.query_range_all
        assert [(1318045, 1319275)] == hsp.hit_range_all
        assert [] == hsp.query_inter_ranges
        assert [] == hsp.hit_inter_ranges
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert 1 == len(hsp.query_all)
        assert 1 == len(hsp.hit_all)
        assert "CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGA" == hsp.query_all[0].seq[:40]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][:40]
        assert "CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGA" == hsp.hit_all[0].seq[:40]
        assert "CTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCAT" == hsp.query_all[0].seq[-40:]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][-40:]
        assert "CTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCAT" == hsp.hit_all[0].seq[-40:]

    def test_exn_22_m_ner(self):
        """Test parsing exonerate output (exn_22_m_ner.exn)."""
        exn_file = get_file("exn_22_m_ner.exn")
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            assert qresult.id == hit.query_id
            for hsp in hit:
                assert hit.id == hsp.hit_id
                assert qresult.id == hsp.query_id

        assert "gi|296143771|ref|NM_001180731.1|" == qresult.id
        assert "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds" == qresult.description
        assert "exonerate" == qresult.program
        assert "NER:affine:local:dna2dna" == qresult.model
        assert 2 == len(qresult)
        # first hit
        hit = qresult[0]
        assert "gi|330443520|ref|NC_001136.10|" == hit.id
        assert "Saccharomyces cerevisiae S288c chromosome IV, complete sequence" == hit.description
        assert 2 == len(hit.hsps)
        # first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 6150 == hsp.score
        assert 1 == hsp[0].query_strand
        assert -1 == hsp[0].hit_strand
        assert 0 == hsp.query_start
        assert 1318045 == hsp.hit_start
        assert 1230 == hsp.query_end
        assert 1319275 == hsp.hit_end
        assert [(0, 1230)] == hsp.query_range_all[:5]
        assert [(1318045, 1319275)] == hsp.hit_range_all[:5]
        assert [] == hsp.query_inter_ranges[:5]
        assert [] == hsp.hit_inter_ranges[:5]
        assert 1 == len(hsp.query_all)
        assert 1 == len(hsp.hit_all)
        assert "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG" == hsp.query_all[0].seq[:40]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][:40]
        assert "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG" == hsp.hit_all[0].seq[:40]
        assert "TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG" == hsp.query_all[0].seq[-40:]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][-40:]
        assert "TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG" == hsp.hit_all[0].seq[-40:]

        # first hit, second hsp
        hsp = qresult[0].hsps[1]
        assert 440 == hsp.score
        assert 1 == hsp[0].query_strand
        assert 1 == hsp[0].hit_strand
        assert 509 == hsp.query_start
        assert 183946 == hsp.hit_start
        assert 1192 == hsp.query_end
        assert 184603 == hsp.hit_end
        assert [(509, 514), (537, 547), (567, 595), (607, 617), (636, 650)] == hsp.query_range_all[:5]
        assert [
                (183946, 183951),
                (183977, 183987),
                (184002, 184030),
                (184044, 184054),
                (184066, 184080),
            ] == hsp.hit_range_all[:5]
        assert [(514, 537), (547, 567), (595, 607), (617, 636), (650, 667)] == hsp.query_inter_ranges[:5]
        assert [
                (183951, 183977),
                (183987, 184002),
                (184030, 184044),
                (184054, 184066),
                (184080, 184092),
            ] == hsp.hit_inter_ranges[:5]
        assert 24 == len(hsp.query_all)
        assert 24 == len(hsp.hit_all)
        # first block
        assert "TGAGA" == hsp.query_all[0].seq
        assert "|||||" == hsp[0].aln_annotation["similarity"]
        assert "TGAGA" == hsp.hit_all[0].seq
        # last block
        assert "GACTGCAAAATAGTAGTCAAAGCTC" == hsp.query_all[-1].seq
        assert "||| | ||||||||||||| | |||" == hsp[-1].aln_annotation["similarity"]
        assert "GACGGTAAAATAGTAGTCACACCTC" == hsp.hit_all[-1].seq

        # second hit
        hit = qresult[1]
        assert "gi|330443681|ref|NC_001144.5|" == hit.id
        assert "Saccharomyces cerevisiae S288c chromosome XII, complete sequence" == hit.description
        assert 1 == len(hit.hsps)
        # second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 502 == hsp.score
        assert 1 == hsp[0].query_strand
        assert 1 == hsp[0].hit_strand
        assert 110 == hsp.query_start
        assert 297910 == hsp.hit_start
        assert 1230 == hsp.query_end
        assert 318994 == hsp.hit_end
        assert [(110, 117), (148, 159), (169, 182), (184, 197), (227, 244)] == hsp.query_range_all[:5]
        assert [
                (297910, 297917),
                (297946, 297957),
                (297970, 297983),
                (297992, 298004),
                (298019, 298038),
            ] == hsp.hit_range_all[:5]
        assert [(117, 148), (159, 169), (182, 184), (197, 227), (244, 255)] == hsp.query_inter_ranges[:5]
        assert [
                (297917, 297946),
                (297957, 297970),
                (297983, 297992),
                (298004, 298019),
                (298038, 298049),
            ] == hsp.hit_inter_ranges[:5]
        assert 33 == len(hsp.query_all)
        assert 33 == len(hsp.hit_all)
        # first block
        assert "CAGAAAA" == hsp.query_all[0].seq
        assert "| |||||" == hsp[0].aln_annotation["similarity"]
        assert "CTGAAAA" == hsp.hit_all[0].seq
        # last block
        assert "TGGTTAGACAGCTCCTGTAG" == hsp.query_all[-1].seq
        assert "|| |  ||||||||||||||" == hsp[-1].aln_annotation["similarity"]
        assert "TGATAGGACAGCTCCTGTAG" == hsp.hit_all[-1].seq

    def test_exn_22_q_multiple(self):
        """Test parsing exonerate output (exn_22_q_multiple.exn)."""
        exn_file = get_file("exn_22_q_multiple.exn")
        qresults = list(parse(exn_file, self.fmt))
        assert 2 == len(qresults)
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                assert qresult.id == hit.query_id
                for hsp in hit:
                    assert hit.id == hsp.hit_id
                    assert qresult.id == hsp.query_id

        # test first qresult
        qresult = qresults[0]
        assert "gi|296142823|ref|NM_001178508.1|" == qresult.id
        assert "exonerate" == qresult.program
        assert "est2genome" == qresult.model
        assert 3 == len(qresult)
        # first qresult, first hit
        hit = qresult[0]
        assert "gi|330443482|ref|NC_001134.8|" == hit.id
        assert 1 == len(hit.hsps)
        # first qresult, first hit, first hsp
        # first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 4485 == hsp.score
        assert 1 == hsp[0].query_strand
        assert 1 == hsp[0].hit_strand
        assert 0 == hsp.query_start
        assert 560077 == hsp.hit_start
        assert 897 == hsp.query_end
        assert 560974 == hsp.hit_end
        assert [(0, 897)] == hsp.query_range_all
        assert [(560077, 560974)] == hsp.hit_range_all
        assert [] == hsp.query_inter_ranges
        assert [] == hsp.hit_inter_ranges
        assert 1 == len(hsp.query_all)
        assert 1 == len(hsp.hit_all)
        assert "ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAG" == hsp.query_all[0].seq[:40]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][:40]
        assert "ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAG" == hsp.hit_all[0].seq[:40]
        assert "CAGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATAA" == hsp.query_all[0].seq[-40:]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][-40:]
        assert "CAGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATAA" == hsp.hit_all[0].seq[-40:]

        # first qresult, second hit
        hit = qresult[1]
        assert "gi|330443753|ref|NC_001148.4|" == hit.id
        assert 1 == len(hit.hsps)
        # first qresult, second hit, first hsp
        # second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 941 == hsp.score
        assert 1 == hsp[0].query_strand
        assert -1 == hsp[0].hit_strand
        assert 2 == hsp.query_start
        assert 492033 == hsp.hit_start
        assert 896 == hsp.query_end
        assert 492933 == hsp.hit_end
        assert [(2, 896)] == hsp.query_range_all
        assert [(492033, 492933)] == hsp.hit_range_all
        assert [] == hsp.query_inter_ranges
        assert [] == hsp.hit_inter_ranges
        assert 1 == len(hsp.query_all)
        assert 1 == len(hsp.hit_all)
        assert "GAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTC" == hsp.query_all[0].seq[:40]
        assert "||||  |  || | || |||  |  ||  | || ||  | " == hsp[0].aln_annotation["similarity"][:40]
        assert "GAGC--TCTATGAACAGATTTAAGCAG-TTAGAAAAGCTT" == hsp.hit_all[0].seq[:40]
        assert "C-AGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATA" == hsp.query_all[0].seq[-40:]
        assert "| || ||      ||| ||||| |  ||   ||| |  ||" == hsp[0].aln_annotation["similarity"][-40:]
        assert "CAAGCAGGCTCTGCAT-CACCCTTGGTTTGCAGAGTACTA" == hsp.hit_all[0].seq[-40:]

        # first qresult, third hit
        hit = qresult[2]
        assert "gi|330443520|ref|NC_001136.10|" == hit.id
        assert 1 == len(hit.hsps)
        # first qresult, third hit, first hsp
        # third hit, first hsp
        hsp = qresult[2].hsps[0]
        assert 651 == hsp.score
        assert 1 == hsp[0].query_strand
        assert 1 == hsp[0].hit_strand
        assert 34 == hsp.query_start
        assert 267809 == hsp.hit_start
        assert 721 == hsp.query_end
        assert 300717 == hsp.hit_end
        assert [(34, 691), (691, 721)] == hsp.query_range_all
        assert [(267809, 268448), (300686, 300717)] == hsp.hit_range_all
        assert [(691, 691)] == hsp.query_inter_ranges
        assert [(268448, 300686)] == hsp.hit_inter_ranges
        assert 2 == len(hsp.query_all)
        assert 2 == len(hsp.hit_all)
        # first block
        assert "AGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGC" == hsp.query_all[0].seq[:40]
        assert "||||||| ||||| ||||| || |  ||||||||    | " == hsp[0].aln_annotation["similarity"][:40]
        assert "AGAAAGTTGGTGAGGGTACTTATGCGGTTGTTTA-CTTGG" == hsp.hit_all[0].seq[:40]
        assert "CGATCAGATTTTCAAG--ATATTCAGAGTATTGGGAACGC" == hsp.query_all[0].seq[-40:]
        assert "|||||| |  |  |||  |  ||||| |  || || || |" == hsp[0].aln_annotation["similarity"][-40:]
        assert "CGATCAAA--TGGAAGTAACGTTCAGGGCCTTAGGGACAC" == hsp.hit_all[0].seq[-40:]
        # last block
        assert "CGAATGAAGCTA-TATGGCCAGATATTGTCT" == hsp.query_all[-1].seq
        assert "| ||   || ||  ||||||||||||| |||" == hsp[-1].aln_annotation["similarity"]
        assert "CAAACCGAGATAGAATGGCCAGATATTCTCT" == hsp.hit_all[-1].seq

        # test second qresult
        qresult = qresults[1]
        assert "gi|296143771|ref|NM_001180731.1|" == qresult.id
        assert "exonerate" == qresult.program
        assert "est2genome" == qresult.model
        assert 2 == len(qresult)
        # second qresult, first hit
        hit = qresult[0]
        assert "gi|330443520|ref|NC_001136.10|" == hit.id
        assert 1 == len(hit.hsps)
        # second qresult, first hit, first hsp
        # first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 6150 == hsp.score
        assert 1 == hsp[0].query_strand
        assert -1 == hsp[0].hit_strand
        assert 0 == hsp.query_start
        assert 1318045 == hsp.hit_start
        assert 1230 == hsp.query_end
        assert 1319275 == hsp.hit_end
        assert [(0, 1230)] == hsp.query_range_all
        assert [(1318045, 1319275)] == hsp.hit_range_all
        assert [] == hsp.query_inter_ranges
        assert [] == hsp.hit_inter_ranges
        assert 1 == len(hsp.query_all)
        assert 1 == len(hsp.hit_all)
        assert "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG" == hsp.query_all[0].seq[:40]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][:40]
        assert "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG" == hsp.hit_all[0].seq[:40]
        assert "TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG" == hsp.query_all[0].seq[-40:]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][-40:]
        assert "TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG" == hsp.hit_all[0].seq[-40:]

        # second qresult, second hit
        hit = qresult[1]
        assert "gi|330443688|ref|NC_001145.3|" == hit.id
        assert 2 == len(hit.hsps)
        # second qresult, second hit, first hsp
        # second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 439 == hsp.score
        assert 1 == hsp[0].query_strand
        assert 1 == hsp[0].hit_strand
        assert 0 == hsp.query_start
        assert 85010 == hsp.hit_start
        assert 346 == hsp.query_end
        assert 473201 == hsp.hit_end
        assert [(0, 65), (65, 225), (225, 320), (320, 346)] == hsp.query_range_all
        assert [(85010, 85066), (253974, 254135), (350959, 351052), (473170, 473201)] == hsp.hit_range_all
        assert [(65, 65), (225, 225), (320, 320)] == hsp.query_inter_ranges
        assert [(85066, 253974), (254135, 350959), (351052, 473170)] == hsp.hit_inter_ranges
        assert 4 == len(hsp.query_all)
        assert 4 == len(hsp.hit_all)
        # first block
        assert "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG" == hsp.query_all[0].seq[:40]
        assert "||||  ||  | ||||   | ||||||  |||| | | | " == hsp[0].aln_annotation["similarity"][:40]
        assert "ATGGTGAACCT-CTTCAAGACGGTCAG--AATA-A-TCAA" == hsp.hit_all[0].seq[:40]
        assert "AGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTT" == hsp.query_all[0].seq[-40:]
        assert "||  |||| | | | ||||    ||||||||||||| | |" == hsp[0].aln_annotation["similarity"][-40:]
        assert "AG--AATA-A-TCAACAGG----ATGAAGAAGCAAAAGAT" == hsp.hit_all[0].seq[-40:]
        # last block
        assert "AGCTAAGAATTCTGATGATG-----AAAGAA" == hsp.query_all[-1].seq
        assert "|   |||||||||||| |||     ||||||" == hsp[-1].aln_annotation["similarity"]
        assert "ATGGAAGAATTCTGATAATGCTGTAAAAGAA" == hsp.hit_all[-1].seq

        # second qresult, second hit, second hsp
        # second hit, second hsp
        hsp = qresult[1].hsps[1]
        assert 263 == hsp.score
        assert 1 == hsp[0].query_strand
        assert -1 == hsp[0].hit_strand
        assert 25 == hsp.query_start
        assert 11338 == hsp.hit_start
        assert 406 == hsp.query_end
        assert 130198 == hsp.hit_end
        assert [(25, 183), (183, 252), (252, 406)] == hsp.query_range_all
        assert [(130038, 130198), (120612, 120681), (11338, 11487)] == hsp.hit_range_all
        assert [(183, 183), (252, 252)] == hsp.query_inter_ranges
        assert [(120681, 130038), (11487, 120612)] == hsp.hit_inter_ranges
        assert 3 == len(hsp.query_all)
        assert 3 == len(hsp.hit_all)
        # first block
        assert "AGCAAATATATTTA-GCAGGTGACATGAAGAAGCAAATGT" == hsp.query_all[0].seq[:40]
        assert "| |||| |||   | ||||   | | || |||| | |  |" == hsp[0].aln_annotation["similarity"][:40]
        assert "ACCAAAGATAACAAGGCAG--AAAAAGAGGAAGAAGAAAT" == hsp.hit_all[0].seq[:40]
        assert "AG-GACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGA" == hsp.query_all[0].seq[-40:]
        assert "|| |||  ||||||  ||   |||  || ||   ||| ||" == hsp[0].aln_annotation["similarity"][-40:]
        assert "AGTGAC--CCCAGAGGAGCCAAGCAAAAAGA---TTCGGA" == hsp.hit_all[0].seq[-40:]
        # last block
        assert "AATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAA" == hsp.query_all[-1].seq[:40]
        assert "|||||||  | ||| |    |||| | |  | | ||    " == hsp[-1].aln_annotation["similarity"][:40]
        assert "AATAAGAGCAACACAG----TTTA-TCTTATATGTA----" == hsp.hit_all[-1].seq[:40]
        assert "CTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAA" == hsp.query_all[-1].seq[-40:]
        assert "|  | | || |  | || ||  ||||||||  ||  ||||" == hsp[-1].aln_annotation["similarity"][-40:]
        assert "CCACTAAAAAATTATAAGAGCCAAAACGAAGTAGATACAA" == hsp.hit_all[-1].seq[-40:]

    def test_exn_22_m_coding2coding_fshifts(self):
        """Test parsing exonerate output (exn_22_m_coding2coding_fshifts.exn)."""
        exn_file = get_file("exn_22_m_coding2coding_fshifts.exn")
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            assert qresult.id == hit.query_id
            for hsp in hit:
                assert hit.id == hsp.hit_id
                assert qresult.id == hsp.query_id

        assert "gi|296143771|ref|NM_001180731.1|" == qresult.id
        assert "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds" == qresult.description
        assert "exonerate" == qresult.program
        assert "coding2coding" == qresult.model
        assert 1 == len(qresult)
        # first hit
        hit = qresult[0]
        assert "gi|296143771|ref|NM_001180731.1|" == hit.id
        assert "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds" == hit.description
        assert 2 == len(hit)
        # first hit, first hsp
        hsp = qresult[0][0]
        assert 213 == hsp.score
        assert [1, 1, 1, 1] == hsp.query_strand_all
        assert [1, 1, 1, 1] == hsp.hit_strand_all
        assert 0 == hsp.query_start
        assert 465 == hsp.hit_start
        assert 160 == hsp.query_end
        assert 630 == hsp.hit_end
        assert [(0, 93), (94, 127), (127, 139), (139, 160)] == hsp.query_range_all
        assert [(465, 558), (558, 591), (593, 605), (609, 630)] == hsp.hit_range_all
        assert [(93, 94), (127, 127), (139, 139)] == hsp.query_inter_ranges
        assert [(558, 558), (591, 593), (605, 609)] == hsp.hit_inter_ranges
        assert [1, 2, 2, 2] == hsp.query_frame_all
        assert [1, 1, 3, 1] == hsp.hit_frame_all
        assert 4 == len(hsp.query_all)
        assert 4 == len(hsp.hit_all)
        assert 4 == len(hsp.aln_annotation_all)
        # first block
        assert "ACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAG" == hsp[0].query.seq[:40]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][:40]
        assert "ACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAG" == hsp[0].hit.seq[:40]
        assert "TATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACC" == hsp[0].query.seq[-40:]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[0].aln_annotation["similarity"][-40:]
        assert "TATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACC" == hsp[0].hit.seq[-40:]
        # last block
        assert "GACGAAAGTATTAATGGTAGT" == hsp[-1].query.seq
        assert "|||||||||||||||||||||" == hsp[-1].aln_annotation["similarity"]
        assert "GACGAAAGTATTAATGGTAGT" == hsp[-1].hit.seq

        # first hit, second hsp
        hsp = qresult[0][1]
        assert 201 == hsp.score
        assert [-1, -1] == hsp.query_strand_all
        assert [-1, -1] == hsp.hit_strand_all
        assert 1 == hsp.query_start
        assert 466 == hsp.hit_start
        assert 158 == hsp.query_end
        assert 628 == hsp.hit_end
        assert [(95, 158), (1, 94)] == hsp.query_range_all
        assert [(559, 628), (466, 559)] == hsp.hit_range_all
        assert [(94, 95)] == hsp.query_inter_ranges
        assert [(559, 559)] == hsp.hit_inter_ranges
        assert [-3, -2] == hsp.query_frame_all
        assert [-2, -2] == hsp.hit_frame_all
        assert 2 == len(hsp.query_all)
        assert 2 == len(hsp.hit_all)
        assert 2 == len(hsp.aln_annotation_all)
        # first block
        assert "TACCATTAATACTTTCGTCATGGT<-><->AACGGCATGT" == hsp[0].query.seq[:40]
        assert "||||||||||||||||||||+ !       ...  !:!!|" == hsp[0].aln_annotation["similarity"][:40]
        assert "TACCATTAATACTTTCGTCACCGATGGTAACGGCACCTGT" == hsp[0].hit.seq[:40]
        # last block
        assert "TGGTTAAACTGGGAGCGTTCCATGAGCCTATGTTCCACAT" == hsp[-1].query.seq[:40]
        assert "||||||||||||||||||||||||||||||||||||||||" == hsp[-1].aln_annotation["similarity"][:40]
        assert "TGGTTAAACTGGGAGCGTTCCATGAGCCTATGTTCCACAT" == hsp[-1].hit.seq[:40]

    def test_exn_22_m_protein2dna_fshifts(self):
        """Test parsing exonerate output (exn_22_m_protein2dna_fshifts.exn)."""
        exn_file = get_file("exn_22_m_protein2dna_fshifts.exn")
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            assert qresult.id == hit.query_id
            for hsp in hit:
                assert hit.id == hsp.hit_id
                assert qresult.id == hsp.query_id

        assert "sp|P24813|YAP2_YEAST" == qresult.id
        assert ("AP-1-like transcription activator YAP2 OS=Saccharomyces cerevisiae (strain"
            " ATCC 204508 / S288c) GN=CAD1 PE=1 SV=2" == qresult.description)
        assert "exonerate" == qresult.program
        assert "protein2dna:local" == qresult.model
        assert 1 == len(qresult)
        # first hit
        hit = qresult[0]
        assert "gi|296143771|ref|NM_001180731.1|" == hit.id
        assert "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds" == hit.description
        assert 2 == len(hit)
        # first hit, first hsp
        hsp = qresult[0][0]
        assert 367 == hsp.score
        assert [0, 0] == hsp.query_strand_all
        assert [1, 1] == hsp.hit_strand_all
        assert 330 == hsp.query_start
        assert 216 == hsp.hit_start
        assert 409 == hsp.query_end
        assert 455 == hsp.hit_end
        assert [(330, 373), (373, 409)] == hsp.query_range_all
        assert [(216, 345), (347, 455)] == hsp.hit_range_all
        assert [(373, 373)] == hsp.query_inter_ranges
        assert [(345, 347)] == hsp.hit_inter_ranges
        assert [0, 0] == hsp.query_frame_all
        assert [1, 3] == hsp.hit_frame_all
        assert 2 == len(hsp.query_all)
        assert 2 == len(hsp.hit_all)
        assert 2 == len(hsp.aln_annotation_all)
        # first block
        assert "HTKTIRTQSEAIEHISSAISNGKASCYHILEEISSLPKYS" == hsp[0].query.seq[:40]
        assert "HTKTIRTQSEAIEHISSAISNGKASCYHILEEISSLPKYS" == hsp[0].hit.seq[:40]
        assert "TIRTQSEAIEHISSAISNGKASCYHILEEISSLPKYSSLD" == hsp[0].query.seq[-40:]
        assert "TIRTQSEAIEHISSAISNGKASCYHILEEISSLPKYSSLD" == hsp[0].hit.seq[-40:]
        # last block
        assert "IDDLCSELIIKAKCTDDCKIVVKARDLQSALVRQLL" == hsp[-1].query.seq
        assert "IDDLCSELIIKAKCTDDCKIVVKARDLQSALVRQLL" == hsp[-1].hit.seq

        # first hit, second hsp
        hsp = qresult[0][1]
        assert 322 == hsp.score
        assert [0] == hsp.query_strand_all
        assert [1] == hsp.hit_strand_all
        assert 6 == hsp.query_start
        assert 16 == hsp.hit_start
        assert 70 == hsp.query_end
        assert 208 == hsp.hit_end
        assert [(6, 70)] == hsp.query_range_all
        assert [(16, 208)] == hsp.hit_range_all
        assert [] == hsp.query_inter_ranges
        assert [] == hsp.hit_inter_ranges
        assert [0] == hsp.query_frame_all
        assert [2] == hsp.hit_frame_all
        assert 1 == len(hsp.query_all)
        assert 1 == len(hsp.hit_all)
        assert 1 == len(hsp.aln_annotation_all)
        assert "KGQQIYLAGDMKKQMLLNKDGTPKRKVGRPGRKRIDSEAK" == hsp[0].query.seq[:40]
        assert "KGQQIYLAGDMKKQMLLNKDGTPKRKVGRPGRKRIDSEAK" == hsp[0].hit.seq[:40]
        assert "RKVGRPGRKRIDSEAKSRRTAQNRAAQRAFRDRKEAKMKS" == hsp[0].query.seq[-40:]
        assert "RKVGRPGRKRIDSEAKSRRTAQNRAAQRAFRDRKEAKMKS" == hsp[0].hit.seq[-40:]

    def test_exn_22_m_protein2genome(self):
        """Test parsing exonerate output (exn_22_m_protein2genome.exn)."""
        exn_file = get_file("exn_22_m_protein2genome.exn")
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            assert qresult.id == hit.query_id
            for hsp in hit:
                assert hit.id == hsp.hit_id
                assert qresult.id == hsp.query_id

        assert "sp|P24813|YAP2_YEAST" == qresult.id
        assert ("AP-1-like transcription activator YAP2 OS=Saccharomyces cerevisiae (strain"
            " ATCC 204508 / S288c) GN=CAD1 PE=1 SV=2" == qresult.description)
        assert "exonerate" == qresult.program
        assert "protein2genome:local" == qresult.model
        assert 3 == len(qresult)
        # first hit
        hit = qresult[0]
        assert "gi|330443520|ref|NC_001136.10|" == hit.id
        assert "Saccharomyces cerevisiae S288c chromosome IV, complete sequence" == hit.description
        assert 1 == len(hit)
        # first hit, first hsp
        hsp = qresult[0][0]
        assert 2105 == hsp.score
        assert 0 == hsp.query_strand
        assert -1 == hsp.hit_strand
        assert 0 == hsp.query_start
        assert 1318048 == hsp.hit_start
        assert 409 == hsp.query_end
        assert 1319275 == hsp.hit_end
        assert "MGNILRKGQQIYLAGDMKKQMLLNKDGTPKRKVGRPGRKR" == hsp[0].query.seq[:40]
        assert "MGNILRKGQQIYLAGDMKKQMLLNKDGTPKRKVGRPGRKR" == hsp[0].hit.seq[:40]
        assert "SSLDIDDLCSELIIKAKCTDDCKIVVKARDLQSALVRQLL" == hsp[0].query.seq[-40:]
        assert "SSLDIDDLCSELIIKAKCTDDCKIVVKARDLQSALVRQLL" == hsp[0].hit.seq[-40:]

        # last hit
        hit = qresult[-1]
        assert "gi|330443590|ref|NC_001140.6|" == hit.id
        assert "Saccharomyces cerevisiae S288c chromosome VIII, complete sequence" == hit.description
        assert 1 == len(hit)
        # last hit, first hsp
        hsp = qresult[-1][0]
        assert 122 == hsp.score
        assert [0, 0] == hsp.query_strand_all
        assert [-1, -1] == hsp.hit_strand_all
        assert "RKRIDSEAKSRRTAQNRAAQRAFRDRKEAKMKSLQERX" == hsp[0].query.seq
        assert "NENVPDDSKAKKKAQNRAAQKAFRERKEARMKELQDKX" == hsp[0].hit.seq
        assert "!.!" == hsp.aln_annotation_all[0]["similarity"][0]
        assert ":!" == hsp.aln_annotation_all[0]["similarity"][-1]
        assert "AAT" == hsp.aln_annotation_all[0]["hit_annotation"][0]
        assert "TT" == hsp.aln_annotation_all[0]["hit_annotation"][-1]
        assert "XELLEQKDAQNKTTTDFLLCSLKSLLSEITKYRAKNSDDERILAFLDDLQE" == hsp[-1].query.seq
        assert "XNKILNRDPQFMSNSSFHQCVSLDSINTIEKDEEKNSDDDAGLQAATDARE" == hsp[-1].hit.seq
        assert "!" == hsp.aln_annotation_all[-1]["similarity"][0]
        assert "|||" == hsp.aln_annotation_all[-1]["similarity"][-1]
        assert "A" == hsp.aln_annotation_all[-1]["hit_annotation"][0]
        assert "GAA" == hsp.aln_annotation_all[-1]["hit_annotation"][-1]

        assert [(37, 74), (75, 125)] == hsp.query_range_all
        assert [(84533, 84646), (68450, 68601)] == hsp.hit_range_all
        assert [(74, 75)] == hsp.query_inter_ranges
        assert [(68601, 84533)] == hsp.hit_inter_ranges
        assert [0, 0] == hsp.query_frame_all
        assert [-3, -3] == hsp.hit_frame_all
        assert 2 == len(hsp.query_all)
        assert 2 == len(hsp.hit_all)
        assert 2 == len(hsp.aln_annotation_all)

    def test_exn_24_m_protein2genome_revcomp_fshifts(self):
        """Test parsing exonerate output (exn_24_m_protein2genome_revcomp_fshifts.exn)."""
        exn_file = get_file("exn_24_m_protein2genome_revcomp_fshifts.exn")
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            assert qresult.id == hit.query_id
            for hsp in hit:
                assert hit.id == hsp.hit_id
                assert qresult.id == hsp.query_id

        assert "Morus-gene026" == qresult.id
        assert "" == qresult.description
        assert "exonerate" == qresult.program
        assert "protein2genome:local" == qresult.model
        assert 1 == len(qresult)
        # first (only) hit
        hit = qresult[0]
        assert "NODE_2_length_1708_cov_48.590765" == hit.id
        assert "SPAdes NODE_2 contig" == hit.description
        assert 1 == len(hit)
        # first hit, first hsp
        hsp = qresult[0][0]
        assert 1308 == hsp.score
        assert [0, 0] == hsp.query_strand_all
        assert [-1, -1] == hsp.hit_strand_all
        assert 69 == hsp.query_start
        assert 331 == hsp.hit_start
        assert 441 == hsp.query_end
        assert 1416 == hsp.hit_end
        assert [(69, 402), (402, 441)] == hsp.query_range_all
        assert [(450, 1416), (331, 448)] == hsp.hit_range_all
        assert [(402, 402)] == hsp.query_inter_ranges
        assert [(448, 450)] == hsp.hit_inter_ranges
        assert [0, 0] == hsp.query_frame_all
        assert [-1, -2] == hsp.hit_frame_all
        assert 2 == len(hsp.query_all)
        assert 2 == len(hsp.hit_all)
        assert 2 == len(hsp.aln_annotation_all)
        # first block
        assert "PESPWTCSPLQT--PSPSLLYHCIASLHRHDGTIHSIAVS" == hsp[0].query.seq[:40]
        assert "PESPWTSSPLQTLSHSPSLLYHCIASLRRHDGTIYSIATS" == hsp[0].hit.seq[:40]
        assert "VEKMVFSGSEDTTIRIWRREEGSCLHECLAVLDGHRGPVK" == hsp[0].query.seq[-40:]
        assert "VEKMVFGGSEDTTIRIWRREEGGCFHKCLAVLDGHRXXXX" == hsp[0].hit.seq[-40:]
        # last block
        assert "CLAACLEVEKVVMMGFLVYSASLDQTFKVWRVKVLPDEE" == hsp[-1].query.seq
        assert "CLAAC*QVEKMVMMGFLIYSVSLDQTLKVWRVKILPDQE" == hsp[-1].hit.seq

    def test_exn_24_protein2genome_met_intron(self):
        """Test parsing exonerate output (exn_24_m_protein2genome_met_intron.exn)."""
        exn_file = get_file("exn_24_m_protein2genome_met_intron.exn")
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            assert qresult.id == hit.query_id
            for hsp in hit:
                assert hit.id == hsp.hit_id
                assert qresult.id == hsp.query_id

        assert "Morus-gene001" == qresult.id
        assert "" == qresult.description
        assert "exonerate" == qresult.program
        assert "protein2genome:local" == qresult.model
        assert 1 == len(qresult)
        # first hit
        hit = qresult[0]
        assert "NODE_1_length_2817_cov_100.387732" == hit.id
        assert "SPAdes contig NODE_1" == hit.description
        assert 1 == len(hit)
        # first hit, first hsp
        assert 1978 == hsp.score
        assert [0, 0, 0, 0, 0, 0] == hsp.query_strand_all
        assert [-1, -1, -1, -1, -1, -1] == hsp.hit_strand_all
        assert 48 == hsp.query_start
        assert 388 == hsp.hit_start
        assert 482 == hsp.query_end
        assert 2392 == hsp.hit_end
        assert [(48, 85), (85, 118), (118, 155), (155, 256), (257, 303), (303, 482)] == hsp.query_range_all
        assert [
                (2281, 2392),
                (2030, 2129),
                (1810, 1921),
                (1420, 1724),
                (1058, 1198),
                (388, 925),
            ] == hsp.hit_range_all
        assert [(2129, 2281), (1921, 2030), (1724, 1810), (1198, 1420), (925, 1058)] == hsp.hit_inter_ranges
        assert [(85, 85), (118, 118), (155, 155), (256, 257), (303, 303)] == hsp.query_inter_ranges
        assert "MVQTPLHVSAGNNRADIVKF" == hsp[0].query.seq[:20]
        assert "VKFLLEFPGPEKVELEAKNM" == hsp[0].query.seq[-20:]
        assert "|||" == hsp[0].aln_annotation["similarity"][0]
        assert "|||" == hsp[0].aln_annotation["similarity"][-1]
        assert "ATG" == hsp[0].aln_annotation["hit_annotation"][0]
        assert "ATG" == hsp[0].aln_annotation["hit_annotation"][-1]
        assert [0, 0, 0, 0, 0, 0] == hsp.query_frame_all
        assert [-2, -3, -2, -2, -3, -2] == hsp.hit_frame_all
        assert 6 == len(hsp.query_all)
        assert 6 == len(hsp.hit_all)
        assert 6 == len(hsp.aln_annotation_all)

    def test_exn_22_q_none(self):
        """Test parsing exonerate output (exn_22_q_none.exn)."""
        exn_file = get_file("exn_22_q_none.exn")
        qresults = parse(exn_file, "exonerate-text")
        with pytest.raises(StopIteration):
            next(qresults)


class ExonerateVulgarCases(unittest.TestCase):
    fmt = "exonerate-vulgar"

    def test_exn_22_o_vulgar(self):
        """Test parsing exonerate output (exn_22_o_vulgar.exn)."""
        exn_file = get_file("exn_22_o_vulgar.exn")
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            assert qresult.id == hit.query_id
            for hsp in hit:
                assert hit.id == hsp.hit_id
                assert qresult.id == hsp.query_id

        assert "sacCer3_dna" == qresult.id
        assert "<unknown description>" == qresult.description
        assert "exonerate" == qresult.program
        assert 3 == len(qresult)
        # first hit
        hit = qresult[0]
        assert "gi|330443520|ref|NC_001136.10|" == hit.id
        assert "<unknown description>" == hit.description
        assert 2 == len(hit.hsps)
        # first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 2641 == hsp.score
        assert -1 == hsp[0].query_strand
        assert -1 == hsp[0].hit_strand
        assert 0 == hsp.query_start
        assert 1319468 == hsp.hit_start
        assert 529 == hsp.query_end
        assert 1319997 == hsp.hit_end
        assert [(0, 529)] == hsp.query_range_all[:5]
        assert [(1319468, 1319997)] == hsp.hit_range_all[:5]
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert " M 26 26 C 3 3 M 500 500" == hsp.vulgar_comp
        # first hit, second hsp
        hsp = qresult[0].hsps[1]
        assert 2641 == hsp.score
        assert 1 == hsp[0].query_strand
        assert 1 == hsp[0].hit_strand
        assert 0 == hsp.query_start
        assert 1319468 == hsp.hit_start
        assert 529 == hsp.query_end
        assert 1319997 == hsp.hit_end
        assert [(0, 529)] == hsp.query_range_all[:5]
        assert [(1319468, 1319997)] == hsp.hit_range_all[:5]
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert " M 90 90 C 3 3 M 436 436" == hsp.vulgar_comp
        # second hit
        hit = qresult[1]
        assert "gi|330443489|ref|NC_001135.5|" == hit.id
        assert "<unknown description>" == hit.description
        assert 1 == len(hit.hsps)
        # second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 267 == hsp.score
        assert -1 == hsp[0].query_strand
        assert 1 == hsp[0].hit_strand
        assert 162 == hsp.query_start
        assert 23668 == hsp.hit_start
        assert 491 == hsp.query_end
        assert 115569 == hsp.hit_end
        assert [(462, 491), (413, 462), (378, 413), (302, 378), (162, 302)] == hsp.query_range_all[:5]
        assert [
                (23668, 23697),
                (32680, 32732),
                (42287, 42325),
                (97748, 97821),
                (115419, 115569),
            ] == hsp.hit_range_all[:5]
        assert [(378, 379), (376, 378)] == hsp.query_split_codons
        assert [(42324, 42325), (97748, 97750)] == hsp.hit_split_codons
        assert (" M 29 29 5 0 2 I 0 8979 3 0 2 M 32 32 G 0 2 M 2 2 G 0 1 M 15 15 5 0 2 I 0 "
            "9551 3 0 2 M 3 3 G 1 0 M 5 5 G 0 2 M 3 3 G 0 1 M 4 4 G 0 1 M 18 18 S 1 1 "
            "5 0 2 I 0 55419 3 0 2 S 2 2 C 3 3 M 22 22 G 3 0 M 46 46 5 0 2 I 0 17594 "
            "3 0 2 M 14 14 G 0 1 M 9 9 G 1 0 M 15 15 G 0 3 M 17 17 G 0 3 M 1 1 G 0 1 "
            "M 13 13 G 0 1 M 6 6 G 1 0 M 12 12 G 0 2 M 45 45 G 0 1 M 6 6" == hsp.vulgar_comp)
        # third hit
        hit = qresult[2]
        assert "gi|330443667|ref|NC_001143.9|" == hit.id
        assert "<unknown description>" == hit.description
        assert 1 == len(hit.hsps)
        # third hit, first hsp
        hsp = qresult[2].hsps[0]
        assert 267 == hsp.score
        assert -1 == hsp[0].query_strand
        assert -1 == hsp[0].hit_strand
        assert 78 == hsp.query_start
        assert 71883 == hsp.hit_start
        assert 529 == hsp.query_end
        assert 641760 == hsp.hit_end
        assert [(449, 529), (319, 388), (198, 284), (161, 198), (78, 114)] == hsp.query_range_all[:5]
        assert [
                (641682, 641760),
                (487327, 487387),
                (386123, 386207),
                (208639, 208677),
                (71883, 71917),
            ] == hsp.hit_range_all[:5]
        assert [(198, 200), (197, 198)] == hsp.query_split_codons
        assert [(386123, 386125), (208676, 208677)] == hsp.hit_split_codons
        assert (" M 31 31 G 3 0 M 4 4 G 2 0 M 19 19 G 0 3 M 9 9 G 0 1 M 6 6 G 1 0 M 5 5 5 "
            "2 2 I 0 154244 I 57 0 I 0 47 3 2 2 M 25 25 G 5 0 M 4 4 G 1 0 M 3 3 G 3 0 "
            "M 4 4 G 1 0 M 9 9 G 0 1 M 14 14 5 2 2 I 0 101116 I 31 0 3 2 2 M 23 23 G "
            "0 1 M 15 15 G 1 0 M 9 9 G 1 0 M 2 2 G 1 0 M 14 14 C 18 18 S 2 2 5 0 2 I "
            "0 177442 3 0 2 S 1 1 C 12 12 M 2 2 G 0 1 M 22 22 5 2 2 I 0 136697 I 7 0 "
            "I 0 6 I 1 0 I 0 1 I 1 0 I 0 1 I 1 0 I 0 1 I 1 0 I 0 1 I 2 0 I 0 1 I 1 0 "
            "I 0 1 I 1 0 I 0 1 I 3 0 I 0 1 I 2 0 I 0 1 I 1 0 I 0 1 I 1 0 I 0 1 I 2 0 "
            "I 0 2 I 2 0 I 0 2 I 17 0 3 2 2 M 12 12 G 2 0 M 22 22" == hsp.vulgar_comp)

    def test_exn_22_o_vulgar_fshifts(self):
        """Test parsing exonerate output (exn_22_o_vulgar_fshifts.exn)."""
        exn_file = get_file("exn_22_o_vulgar_fshifts.exn")
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            assert qresult.id == hit.query_id
            for hsp in hit:
                assert hit.id == hsp.hit_id
                assert qresult.id == hsp.query_id

        assert "gi|296143771|ref|NM_001180731.1|" == qresult.id
        assert "<unknown description>" == qresult.description
        assert "exonerate" == qresult.program
        assert 1 == len(qresult)
        # first hit
        hit = qresult[0]
        assert "gi|296143771|ref|NM_001180731.1|" == hit.id
        assert "<unknown description>" == hit.description
        assert 2 == len(hit)
        # first hit, first hsp
        hsp = qresult[0][0]
        assert 213 == hsp.score
        assert 1 == hsp[0].query_strand
        assert 1 == hsp[0].hit_strand
        assert 0 == hsp.query_start
        assert 465 == hsp.hit_start
        assert 160 == hsp.query_end
        assert 630 == hsp.hit_end
        assert [(0, 93), (94, 127), (127, 139), (139, 160)] == hsp.query_range_all[:5]
        assert [(465, 558), (558, 591), (593, 605), (609, 630)] == hsp.hit_range_all[:5]
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert " C 93 93 F 1 0 C 33 33 F 0 2 C 12 12 F 0 4 C 21 21" == hsp.vulgar_comp
        # first hit, second hsp
        hsp = qresult[0][1]
        assert 201 == hsp.score
        assert -1 == hsp[0].query_strand
        assert -1 == hsp[0].hit_strand
        assert 1 == hsp.query_start
        assert 466 == hsp.hit_start
        assert 158 == hsp.query_end
        assert 628 == hsp.hit_end
        assert [(95, 158), (1, 94)] == hsp.query_range_all[:5]
        assert [(559, 628), (466, 559)] == hsp.hit_range_all[:5]
        assert [] == hsp.query_split_codons
        assert [] == hsp.hit_split_codons
        assert " C 24 24 G 0 6 C 39 39 F 1 0 C 93 93" == hsp.vulgar_comp


class ExonerateCigarCases(unittest.TestCase):
    fmt = "exonerate-cigar"

    def test_exn_22_o_vulgar_cigar(self):
        """Test parsing exonerate output (exn_22_o_vulgar_cigar.exn)."""
        exn_file = get_file("exn_22_o_vulgar_cigar.exn")
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            assert qresult.id == hit.query_id
            for hsp in hit:
                assert hit.id == hsp.hit_id
                assert qresult.id == hsp.query_id

        assert "sacCer3_dna" == qresult.id
        assert "<unknown description>" == qresult.description
        assert "exonerate" == qresult.program
        assert 3 == len(qresult)
        # first hit
        hit = qresult[0]
        assert "gi|330443520|ref|NC_001136.10|" == hit.id
        assert "<unknown description>" == hit.description
        assert 2 == len(hit.hsps)
        # first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 2641 == hsp.score
        assert -1 == hsp[0].query_strand
        assert -1 == hsp[0].hit_strand
        assert 0 == hsp.query_start
        assert 1319468 == hsp.hit_start
        assert 529 == hsp.query_end
        assert 1319997 == hsp.hit_end
        assert "  M 26 M 3 M 500" == hsp.cigar_comp
        # first hit, second hsp
        hsp = qresult[0].hsps[1]
        assert 2641 == hsp.score
        assert 1 == hsp[0].query_strand
        assert 1 == hsp[0].hit_strand
        assert 0 == hsp.query_start
        assert 1319468 == hsp.hit_start
        assert 529 == hsp.query_end
        assert 1319997 == hsp.hit_end
        assert "  M 90 M 3 M 436" == hsp.cigar_comp
        # second hit
        hit = qresult[1]
        assert "gi|330443489|ref|NC_001135.5|" == hit.id
        assert "<unknown description>" == hit.description
        assert 1 == len(hit.hsps)
        # second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 267 == hsp.score
        assert -1 == hsp[0].query_strand
        assert 1 == hsp[0].hit_strand
        assert 162 == hsp.query_start
        assert 23668 == hsp.hit_start
        assert 491 == hsp.query_end
        assert 115569 == hsp.hit_end
        assert ("  M 29 D 8983 M 32 D 2 M 2 D 1 M 15 D 9555 M 3 I 1 M 5 D 2 M 3 D 1 M 4 D "
            "1 M 18 M 1 D 55423 M 5 M 22 I 3 M 46 D 17598 M 14 D 1 M 9 I 1 M 15 D 3 M "
            "17 D 3 M 1 D 1 M 13 D 1 M 6 I 1 M 12 D 2 M 45 D 1 M 6" == hsp.cigar_comp)
        # third hit
        hit = qresult[2]
        assert "gi|330443667|ref|NC_001143.9|" == hit.id
        assert "<unknown description>" == hit.description
        assert 1 == len(hit.hsps)
        # third hit, first hsp
        hsp = qresult[2].hsps[0]
        assert 267 == hsp.score
        assert -1 == hsp[0].query_strand
        assert -1 == hsp[0].hit_strand
        assert 78 == hsp.query_start
        assert 71883 == hsp.hit_start
        assert 529 == hsp.query_end
        assert 641760 == hsp.hit_end
        assert ("  M 31 I 3 M 4 I 2 M 19 D 3 M 9 D 1 M 6 I 1 M 7 D 154244 I 57 D 47 M 27 I "
            "5 M 4 I 1 M 3 I 3 M 4 I 1 M 9 D 1 M 16 D 101116 I 31 M 25 D 1 M 15 I 1 M "
            "9 I 1 M 2 I 1 M 14 M 20 D 177446 M 13 M 2 D 1 M 24 D 136697 I 7 D 6 I 1 D "
            "1 I 1 D 1 I 1 D 1 I 1 D 1 I 2 D 1 I 1 D 1 I 1 D 1 I 3 D 1 I 2 D 1 I 1 D 1 "
            "I 1 D 1 I 2 D 2 I 2 D 2 I 17 M 14 I 2 M 22" == hsp.cigar_comp)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
