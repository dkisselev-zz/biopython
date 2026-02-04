# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO BlastIO parsers."""

import os
import unittest
import pytest
import warnings

from Bio import BiopythonParserWarning
from Bio.SearchIO import parse

# test case files are in the Blast directory
TEST_DIR = "Blast"
FMT = "blast-xml"

REFERENCE = (
    "Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schäffer, "
    "Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "
    '"Gapped BLAST and PSI-BLAST: a new generation of protein database '
    'search programs", Nucleic Acids Res. 25:3389-3402.'
)


def get_file(filename):
    """Return the path of a test file."""
    return os.path.join(TEST_DIR, filename)


class BlastnXmlCases(unittest.TestCase):
    def test_xml_2212L_blastn_001(self):
        xml_file = get_file("xml_2212L_blastn_001.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = next(qresults)
        counter += 1

        assert "2.2.12" == qresult.version
        assert REFERENCE == qresult.reference
        assert 10.0 == qresult.param_evalue_threshold
        assert 1 == qresult.param_score_match
        assert -3 == qresult.param_score_mismatch
        assert 5 == qresult.param_gap_open
        assert 2 == qresult.param_gap_extend

        # test parsed values of qresult
        assert "gi|1348916|gb|G26684.1|G26684" == qresult.id
        assert "human STS STS_D11570, sequence tagged site" == qresult.description
        assert 285 == qresult.seq_len
        assert 371021 == qresult.stat_db_num
        assert 1233631384 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert 0.710603 == qresult.stat_kappa
        assert 1.37406 == qresult.stat_lambda
        assert 1.30725 == qresult.stat_entropy
        assert 2 == len(qresult)

        hit = qresult[0]
        assert "gi|9950606|gb|AE004854.1|" == hit.id
        assert "Pseudomonas aeruginosa PAO1, section 415 of 529 of the complete genome" == hit.description
        assert 11884 == hit.seq_len
        assert 1 == len(hit)
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        hsp = hit.hsps[0]
        assert 38.1576 == hsp.bitscore
        assert 19 == hsp.bitscore_raw
        assert 1.0598 == hsp.evalue
        assert 67 == hsp.query_start
        assert 86 == hsp.query_end
        assert 6011 == hsp.hit_start
        assert 6030 == hsp.hit_end
        assert 1 == hsp.query_frame
        assert 1 == hsp.hit_frame
        assert 19 == hsp.ident_num
        assert 19 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 19 == hsp.aln_span
        assert "CAGGCCAGCGACTTCTGGG" == hsp.query.seq
        assert "CAGGCCAGCGACTTCTGGG" == hsp.hit.seq
        assert "|||||||||||||||||||" == hsp.aln_annotation["similarity"]
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        # parse last hit
        hit = qresult[-1]
        assert "gi|15073988|emb|AL591786.1|SME591786" == hit.id
        assert "Sinorhizobium meliloti 1021 complete chromosome; segment 5/12" == hit.description
        assert 299350 == hit.seq_len
        assert 1 == len(hit)
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        hsp = hit.hsps[0]
        assert 36.1753 == hsp.bitscore
        assert 18 == hsp.bitscore_raw
        assert 4.18768 == hsp.evalue
        assert 203 == hsp.query_start
        assert 224 == hsp.query_end
        assert 83627 == hsp.hit_start
        assert 83648 == hsp.hit_end
        assert 1 == hsp.query_frame
        assert -1 == hsp.hit_frame
        assert 20 == hsp.ident_num
        assert 20 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 21 == hsp.aln_span
        assert "TGAAAGGAAATNAAAATGGAA" == hsp.query.seq
        assert "TGAAAGGAAATCAAAATGGAA" == hsp.hit.seq
        assert "||||||||||| |||||||||" == hsp.aln_annotation["similarity"]
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_xml_2226_blastn_001(self):
        xml_file = get_file("xml_2226_blastn_001.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = next(qresults)
        counter += 1

        # test meta variables, only for the first one
        assert "2.2.26+" == qresult.version
        assert ("Zheng Zhang, Scott Schwartz, Lukas Wagner, and "
            'Webb Miller (2000), "A greedy algorithm for '
            'aligning DNA sequences", J Comput Biol 2000; '
            "7(1-2):203-14." == qresult.reference)
        assert 1.0 == qresult.param_score_match
        assert -2.0 == qresult.param_score_mismatch
        assert 10.0 == qresult.param_evalue_threshold
        assert "L;m;" == qresult.param_filter
        assert 0.0 == qresult.param_gap_open
        assert 0.0 == qresult.param_gap_extend
        assert "blastn" == qresult.program
        assert "db/minirefseq_mrna" == qresult.target

        # test parsed values of the first qresult
        assert "random_s00" == qresult.id
        assert "" == qresult.description
        assert 128 == qresult.seq_len
        assert 23 == qresult.stat_db_num
        assert 67750 == qresult.stat_db_len
        assert 7616765.0 == qresult.stat_eff_space
        assert 0.46 == qresult.stat_kappa
        assert 1.28 == qresult.stat_lambda
        assert 0.85 == qresult.stat_entropy
        assert 0 == len(qresult)

        # test parsed values of the second qresult
        qresult = next(qresults)
        counter += 1
        assert "gi|356995852:1-490" == qresult.id
        assert "Mus musculus POU domain, class 5, transcription factor 1 (Pou5f1), transcript variant 1, mRNA" == qresult.description
        assert 490 == qresult.seq_len
        assert 23 == qresult.stat_db_num
        assert 67750 == qresult.stat_db_len
        assert 31860807.0 == qresult.stat_eff_space
        assert 0.46 == qresult.stat_kappa
        assert 1.28 == qresult.stat_lambda
        assert 0.85 == qresult.stat_entropy
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "gi|356995852|ref|NM_013633.3|" == hit.id
        assert ("Mus musculus POU "
            "domain, class 5, transcription factor 1 (Pou5f1), "
            "transcript variant 1, mRNA" == hit.description)
        assert 1353 == hit.seq_len
        assert 1 == len(hit)
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        hsp = hit.hsps[0]
        assert 905.979 == hsp.bitscore
        assert 490 == hsp.bitscore_raw
        assert 0 == hsp.evalue
        assert 0 == hsp.query_start
        assert 490 == hsp.query_end
        assert 0 == hsp.hit_start
        assert 490 == hsp.hit_end
        assert 1 == hsp.query_frame
        assert 1 == hsp.hit_frame
        assert 490 == hsp.ident_num
        assert 490 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 490 == hsp.aln_span
        assert "GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGGCCCCCGGCTCGGGGTGCCCACCTTCCCCATGGCTGGACACCTGGCTTCAGACTTCGCCTTCTCACCCCCACCAGGTGGGGGTGATGGGTCAGCAGGGCTGGAGCCGGGCTGGGTGGATCCTCGAACCTGGCTAAGCTTCCAAGGGCCTCCAGGTGGGCCTGGAATCGGACCAGGCTCAGAGGTATTGGGGATCTCCCCATGTCCGCCCGCATACGAGTTCTGCGGAGGGATGGCATACTGTGGACCTCAGGTTGGACTGGGCCTAGTCCCCCAAGTTGGCGTGGAGACTTTGCAGCCTGAGGGCCAGGCAGGAGCACGAGTGGAAAGCAACTCAGAGGGAACCTCCTCTGAGCCCTGTGCCGACCGCCCCAATGCCGTGAAGTTGGAGAAGGTGGAACCAACTCCCGAGGAGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA" == hsp.query.seq
        assert "GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGGCCCCCGGCTCGGGGTGCCCACCTTCCCCATGGCTGGACACCTGGCTTCAGACTTCGCCTTCTCACCCCCACCAGGTGGGGGTGATGGGTCAGCAGGGCTGGAGCCGGGCTGGGTGGATCCTCGAACCTGGCTAAGCTTCCAAGGGCCTCCAGGTGGGCCTGGAATCGGACCAGGCTCAGAGGTATTGGGGATCTCCCCATGTCCGCCCGCATACGAGTTCTGCGGAGGGATGGCATACTGTGGACCTCAGGTTGGACTGGGCCTAGTCCCCCAAGTTGGCGTGGAGACTTTGCAGCCTGAGGGCCAGGCAGGAGCACGAGTGGAAAGCAACTCAGAGGGAACCTCCTCTGAGCCCTGTGCCGACCGCCCCAATGCCGTGAAGTTGGAGAAGGTGGAACCAACTCCCGAGGAGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA" == hsp.hit.seq
        assert "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" == hsp.aln_annotation["similarity"]
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        # test parsed values of the third qresult
        qresult = next(qresults)
        counter += 1
        assert "hg19_dna" == qresult.id
        assert "range=chr1:1207307-1207372 5'pad=0 3'pad=0 strand=+ repeatMasking=none" == qresult.description
        assert 66 == qresult.seq_len
        assert 23 == qresult.stat_db_num
        assert 67750 == qresult.stat_db_len
        assert 3506256.0 == qresult.stat_eff_space
        assert 0.46 == qresult.stat_kappa
        assert 1.28 == qresult.stat_lambda
        assert 0.85 == qresult.stat_entropy
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "gi|94721341|ref|NM_001040441.1|" == hit.id
        assert "Homo sapiens zinc finger and BTB domain containing 8A (ZBTB8A), mRNA" == hit.description
        assert 7333 == hit.seq_len
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 115.613 == hsp.bitscore
        assert 62 == hsp.bitscore_raw
        assert 5.52066e-29 == hsp.evalue
        assert 4 == hsp.query_start
        assert 66 == hsp.query_end
        assert 3676 == hsp.hit_start
        assert 3738 == hsp.hit_end
        assert 1 == hsp.query_frame
        assert 1 == hsp.hit_frame
        assert 62 == hsp.ident_num
        assert 62 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 62 == hsp.aln_span
        assert "GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA" == hsp.query.seq
        assert "GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA" == hsp.hit.seq
        assert "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[-1]
        assert 98.9927 == hsp.bitscore
        assert 53 == hsp.bitscore_raw
        assert 5.55986e-24 == hsp.evalue
        assert 5 == hsp.query_start
        assert 58 == hsp.query_end
        assert 2823 == hsp.hit_start
        assert 2876 == hsp.hit_end
        assert 1 == hsp.query_frame
        assert 1 == hsp.hit_frame
        assert 53 == hsp.ident_num
        assert 53 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 53 == hsp.aln_span
        assert "CCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAA" == hsp.query.seq
        assert "CCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAA" == hsp.hit.seq
        assert "|||||||||||||||||||||||||||||||||||||||||||||||||||||" == hsp.aln_annotation["similarity"]

        hit = qresult[-1]
        assert "gi|332865372|ref|XM_003318468.1|" == hit.id
        assert ("PREDICTED: Pan "
            "troglodytes zinc finger protein 273, transcript "
            "variant 1 (ZNF273), mRNA" == hit.description)
        assert 4430 == hit.seq_len
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 111.919 == hsp.bitscore
        assert 60 == hsp.bitscore_raw
        assert 7.14143e-28 == hsp.evalue
        assert 0 == hsp.query_start
        assert 66 == hsp.query_end
        assert 2734 == hsp.hit_start
        assert 2800 == hsp.hit_end
        assert 1 == hsp.query_frame
        assert -1 == hsp.hit_frame
        assert 64 == hsp.ident_num
        assert 64 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 66 == hsp.aln_span
        assert "TCAAGCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA" == hsp.query.seq
        assert "TCACGCCATTGCACTCCAGCCTGGGCAACAAGAGTGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA" == hsp.hit.seq
        assert "||| |||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||" == hsp.aln_annotation["similarity"]

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 3 == counter

    def test_xml_2226_blastn_002(self):
        xml_file = get_file("xml_2226_blastn_002.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = next(qresults)
        counter += 1

        assert "2.2.26+" == qresult.version
        assert ("Zheng Zhang, Scott Schwartz, Lukas Wagner, and "
            'Webb Miller (2000), "A greedy algorithm for '
            'aligning DNA sequences", J Comput Biol 2000; '
            "7(1-2):203-14." == qresult.reference)
        assert 1.0 == qresult.param_score_match
        assert -2.0 == qresult.param_score_mismatch
        assert 10.0 == qresult.param_evalue_threshold
        assert "L;m;" == qresult.param_filter
        assert 0.0 == qresult.param_gap_open
        assert 0.0 == qresult.param_gap_extend
        assert "blastn" == qresult.program
        assert "db/minirefseq_mrna" == qresult.target
        assert "random_s00" == qresult.id
        assert "" == qresult.description
        assert 128 == qresult.seq_len
        assert 23 == qresult.stat_db_num
        assert 67750 == qresult.stat_db_len
        assert 7616765.0 == qresult.stat_eff_space
        assert 0.46 == qresult.stat_kappa
        assert 1.28 == qresult.stat_lambda
        assert 0.85 == qresult.stat_entropy
        assert 0 == len(qresult)
        assert [] == list(qresult.hits)

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_xml_2226_blastn_003(self):
        xml_file = get_file("xml_2226_blastn_003.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = next(qresults)
        counter += 1

        assert "2.2.26+" == qresult.version
        assert ("Zheng Zhang, Scott Schwartz, Lukas Wagner, and "
            'Webb Miller (2000), "A greedy algorithm for '
            'aligning DNA sequences", J Comput Biol 2000; '
            "7(1-2):203-14." == qresult.reference)
        assert 1.0 == qresult.param_score_match
        assert -2.0 == qresult.param_score_mismatch
        assert 10.0 == qresult.param_evalue_threshold
        assert "L;m;" == qresult.param_filter
        assert 0.0 == qresult.param_gap_open
        assert 0.0 == qresult.param_gap_extend
        assert "blastn" == qresult.program
        assert "db/minirefseq_mrna" == qresult.target

        assert "gi|356995852:1-490" == qresult.id
        assert ("Mus musculus POU domain, class 5, transcription "
            "factor 1 (Pou5f1), transcript variant 1, mRNA" == qresult.description)
        assert 490 == qresult.seq_len
        assert 23 == qresult.stat_db_num
        assert 67750 == qresult.stat_db_len
        assert 31860807.0 == qresult.stat_eff_space
        assert 0.46 == qresult.stat_kappa
        assert 1.28 == qresult.stat_lambda
        assert 0.85 == qresult.stat_entropy
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "gi|356995852|ref|NM_013633.3|" == hit.id
        assert ("Mus musculus POU "
            "domain, class 5, transcription factor 1 (Pou5f1), "
            "transcript variant 1, mRNA" == hit.description)
        assert 1353 == hit.seq_len
        assert 1 == len(hit)
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        hsp = hit.hsps[0]
        assert 905.979 == hsp.bitscore
        assert 490 == hsp.bitscore_raw
        assert 0 == hsp.evalue
        assert 0 == hsp.query_start
        assert 490 == hsp.query_end
        assert 0 == hsp.hit_start
        assert 490 == hsp.hit_end
        assert 1 == hsp.query_frame
        assert 1 == hsp.hit_frame
        assert 490 == hsp.ident_num
        assert 490 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 490 == hsp.aln_span
        assert "GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGGCCCCCGGCTCGGGGTGCCCACCTTCCCCATGGCTGGACACCTGGCTTCAGACTTCGCCTTCTCACCCCCACCAGGTGGGGGTGATGGGTCAGCAGGGCTGGAGCCGGGCTGGGTGGATCCTCGAACCTGGCTAAGCTTCCAAGGGCCTCCAGGTGGGCCTGGAATCGGACCAGGCTCAGAGGTATTGGGGATCTCCCCATGTCCGCCCGCATACGAGTTCTGCGGAGGGATGGCATACTGTGGACCTCAGGTTGGACTGGGCCTAGTCCCCCAAGTTGGCGTGGAGACTTTGCAGCCTGAGGGCCAGGCAGGAGCACGAGTGGAAAGCAACTCAGAGGGAACCTCCTCTGAGCCCTGTGCCGACCGCCCCAATGCCGTGAAGTTGGAGAAGGTGGAACCAACTCCCGAGGAGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA" == hsp.query.seq
        assert "GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGGCCCCCGGCTCGGGGTGCCCACCTTCCCCATGGCTGGACACCTGGCTTCAGACTTCGCCTTCTCACCCCCACCAGGTGGGGGTGATGGGTCAGCAGGGCTGGAGCCGGGCTGGGTGGATCCTCGAACCTGGCTAAGCTTCCAAGGGCCTCCAGGTGGGCCTGGAATCGGACCAGGCTCAGAGGTATTGGGGATCTCCCCATGTCCGCCCGCATACGAGTTCTGCGGAGGGATGGCATACTGTGGACCTCAGGTTGGACTGGGCCTAGTCCCCCAAGTTGGCGTGGAGACTTTGCAGCCTGAGGGCCAGGCAGGAGCACGAGTGGAAAGCAACTCAGAGGGAACCTCCTCTGAGCCCTGTGCCGACCGCCCCAATGCCGTGAAGTTGGAGAAGGTGGAACCAACTCCCGAGGAGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA" == hsp.hit.seq
        assert "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" == hsp.aln_annotation["similarity"]
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_xml_2226_blastn_004(self):
        xml_file = get_file("xml_2226_blastn_004.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1
        assert "hg19_dna" == qresult.id
        assert "range=chr1:1207307-1207372 5'pad=0 3'pad=0 strand=+ repeatMasking=none" == qresult.description
        assert 66 == qresult.seq_len
        assert 23 == qresult.stat_db_num
        assert 67750 == qresult.stat_db_len
        assert 3506256.0 == qresult.stat_eff_space
        assert 0.46 == qresult.stat_kappa
        assert 1.28 == qresult.stat_lambda
        assert 0.85 == qresult.stat_entropy
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "gi|94721341|ref|NM_001040441.1|" == hit.id
        assert "Homo sapiens zinc finger and BTB domain containing 8A (ZBTB8A), mRNA" == hit.description
        assert 7333 == hit.seq_len
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 115.613 == hsp.bitscore
        assert 62 == hsp.bitscore_raw
        assert 5.52066e-29 == hsp.evalue
        assert 4 == hsp.query_start
        assert 66 == hsp.query_end
        assert 3676 == hsp.hit_start
        assert 3738 == hsp.hit_end
        assert 1 == hsp.query_frame
        assert 1 == hsp.hit_frame
        assert 62 == hsp.ident_num
        assert 62 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 62 == hsp.aln_span
        assert "GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA" == hsp.query.seq
        assert "GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA" == hsp.hit.seq
        assert "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[-1]
        assert 98.9927 == hsp.bitscore
        assert 53 == hsp.bitscore_raw
        assert 5.55986e-24 == hsp.evalue
        assert 5 == hsp.query_start
        assert 58 == hsp.query_end
        assert 2823 == hsp.hit_start
        assert 2876 == hsp.hit_end
        assert 1 == hsp.query_frame
        assert 1 == hsp.hit_frame
        assert 53 == hsp.ident_num
        assert 53 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 53 == hsp.aln_span
        assert "CCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAA" == hsp.query.seq
        assert "CCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAA" == hsp.hit.seq
        assert "|||||||||||||||||||||||||||||||||||||||||||||||||||||" == hsp.aln_annotation["similarity"]

        hit = qresult[-1]
        assert "gi|332865372|ref|XM_003318468.1|" == hit.id
        assert ("PREDICTED: Pan "
            "troglodytes zinc finger protein 273, transcript "
            "variant 1 (ZNF273), mRNA" == hit.description)
        assert 4430 == hit.seq_len
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 111.919 == hsp.bitscore
        assert 60 == hsp.bitscore_raw
        assert 7.14143e-28 == hsp.evalue
        assert 0 == hsp.query_start
        assert 66 == hsp.query_end
        assert 2734 == hsp.hit_start
        assert 2800 == hsp.hit_end
        assert 1 == hsp.query_frame
        assert -1 == hsp.hit_frame
        assert 64 == hsp.ident_num
        assert 64 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 66 == hsp.aln_span
        assert "TCAAGCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA" == hsp.query.seq
        assert "TCACGCCATTGCACTCCAGCCTGGGCAACAAGAGTGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA" == hsp.hit.seq
        assert "||| |||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||" == hsp.aln_annotation["similarity"]

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_xml_2226_blastn_005(self):
        xml_file = get_file("xml_2226_blastn_005.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = next(qresults)
        counter += 1

        # test meta variables, only for the first one
        assert "2.2.26+" == qresult.version
        assert ("Zheng Zhang, Scott Schwartz, Lukas Wagner, and "
            'Webb Miller (2000), "A greedy algorithm for '
            'aligning DNA sequences", J Comput Biol 2000; '
            "7(1-2):203-14." == qresult.reference)
        assert 1.0 == qresult.param_score_match
        assert -2.0 == qresult.param_score_mismatch
        assert 10.0 == qresult.param_evalue_threshold
        assert "L;m;" == qresult.param_filter
        assert 0.0 == qresult.param_gap_open
        assert 0.0 == qresult.param_gap_extend
        assert "blastn" == qresult.program
        assert "refseq_rna" == qresult.target

        # test parsed values of the first qresult
        assert "random_s00" == qresult.id
        assert "" == qresult.description
        assert 128 == qresult.seq_len
        assert 2933984 == qresult.stat_db_num
        assert 4726730735 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert 0.46 == qresult.stat_kappa
        assert 1.28 == qresult.stat_lambda
        assert 0.85 == qresult.stat_entropy
        assert 0 == len(qresult)

        # test parsed values of the second qresult
        qresult = next(qresults)
        counter += 1
        assert "gi|356995852:1-490" == qresult.id
        assert ("Mus musculus POU domain, class 5, transcription "
            "factor 1 (Pou5f1), transcript variant 1, mRNA" == qresult.description)
        assert 490 == qresult.seq_len
        assert 2933984 == qresult.stat_db_num
        assert 4726730735 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert 0.46 == qresult.stat_kappa
        assert 1.28 == qresult.stat_lambda
        assert 0.85 == qresult.stat_entropy
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "gi|356995852|ref|NM_013633.3|" == hit.id
        assert ("Mus musculus POU "
            "domain, class 5, transcription factor 1 (Pou5f1), "
            "transcript variant 1, mRNA" == hit.description)
        assert 1353 == hit.seq_len
        assert 1 == len(hit)
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        hsp = hit.hsps[0]
        assert 905.979 == hsp.bitscore
        assert 490 == hsp.bitscore_raw
        assert 0 == hsp.evalue
        assert 0 == hsp.query_start
        assert 490 == hsp.query_end
        assert 0 == hsp.hit_start
        assert 490 == hsp.hit_end
        assert 1 == hsp.query_frame
        assert 1 == hsp.hit_frame
        assert 490 == hsp.ident_num
        assert 490 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 490 == hsp.aln_span
        assert "GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGGCCCCCGGCTCGGGGTGCCCACCTTCCCCATGGCTGGACACCTGGCTTCAGACTTCGCCTTCTCACCCCCACCAGGTGGGGGTGATGGGTCAGCAGGGCTGGAGCCGGGCTGGGTGGATCCTCGAACCTGGCTAAGCTTCCAAGGGCCTCCAGGTGGGCCTGGAATCGGACCAGGCTCAGAGGTATTGGGGATCTCCCCATGTCCGCCCGCATACGAGTTCTGCGGAGGGATGGCATACTGTGGACCTCAGGTTGGACTGGGCCTAGTCCCCCAAGTTGGCGTGGAGACTTTGCAGCCTGAGGGCCAGGCAGGAGCACGAGTGGAAAGCAACTCAGAGGGAACCTCCTCTGAGCCCTGTGCCGACCGCCCCAATGCCGTGAAGTTGGAGAAGGTGGAACCAACTCCCGAGGAGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA" == hsp.query.seq
        assert "GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGGCCCCCGGCTCGGGGTGCCCACCTTCCCCATGGCTGGACACCTGGCTTCAGACTTCGCCTTCTCACCCCCACCAGGTGGGGGTGATGGGTCAGCAGGGCTGGAGCCGGGCTGGGTGGATCCTCGAACCTGGCTAAGCTTCCAAGGGCCTCCAGGTGGGCCTGGAATCGGACCAGGCTCAGAGGTATTGGGGATCTCCCCATGTCCGCCCGCATACGAGTTCTGCGGAGGGATGGCATACTGTGGACCTCAGGTTGGACTGGGCCTAGTCCCCCAAGTTGGCGTGGAGACTTTGCAGCCTGAGGGCCAGGCAGGAGCACGAGTGGAAAGCAACTCAGAGGGAACCTCCTCTGAGCCCTGTGCCGACCGCCCCAATGCCGTGAAGTTGGAGAAGGTGGAACCAACTCCCGAGGAGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA" == hsp.hit.seq
        assert "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" == hsp.aln_annotation["similarity"]
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        # test parsed values of the third qresult
        qresult = next(qresults)
        counter += 1
        assert "hg19_dna" == qresult.id
        assert "range=chr1:1207307-1207372 5'pad=0 3'pad=0 strand=+ repeatMasking=none" == qresult.description
        assert 66 == qresult.seq_len
        assert 2933984 == qresult.stat_db_num
        assert 4726730735 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert 0.46 == qresult.stat_kappa
        assert 1.28 == qresult.stat_lambda
        assert 0.85 == qresult.stat_entropy
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "gi|332237160|ref|XM_003267724.1|" == hit.id
        assert ("PREDICTED: Nomascus leucogenys ATG14 autophagy "
            "related 14 homolog (S. cerevisiae) (ATG14), mRNA" == hit.description)
        assert 4771 == hit.seq_len
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 115.613 == hsp.bitscore
        assert 62 == hsp.bitscore_raw
        assert 3.35972e-23 == hsp.evalue
        assert 4 == hsp.query_start
        assert 66 == hsp.query_end
        assert 2864 == hsp.hit_start
        assert 2926 == hsp.hit_end
        assert 1 == hsp.query_frame
        assert 1 == hsp.hit_frame
        assert 62 == hsp.ident_num
        assert 62 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 62 == hsp.aln_span
        assert "GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA" == hsp.query.seq
        assert "GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA" == hsp.hit.seq
        assert "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" == hsp.aln_annotation["similarity"]

        hit = qresult[-1]
        assert "gi|332254616|ref|XM_003276378.1|" == hit.id
        assert ("PREDICTED: Nomascus leucogenys S100P binding "
            "protein, transcript variant 2 (S100PBP), mRNA" == hit.description)
        assert 4345 == hit.seq_len
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 111.919 == hsp.bitscore
        assert 60 == hsp.bitscore_raw
        assert 4.34607e-22 == hsp.evalue
        assert 0 == hsp.query_start
        assert 66 == hsp.query_end
        assert 2791 == hsp.hit_start
        assert 2857 == hsp.hit_end
        assert 1 == hsp.query_frame
        assert -1 == hsp.hit_frame
        assert 64 == hsp.ident_num
        assert 64 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 66 == hsp.aln_span
        assert "TCAAGCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA" == hsp.query.seq
        assert "TCATGCCACTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA" == hsp.hit.seq
        assert "||| |||| |||||||||||||||||||||||||||||||||||||||||||||||||||||||||" == hsp.aln_annotation["similarity"]

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 3 == counter


class BlastpXmlCases(unittest.TestCase):
    def test_xml_2212L_blastp_001(self):
        xml_file = get_file("xml_2212L_blastp_001.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1

        assert "2.2.12" == qresult.version
        assert REFERENCE == qresult.reference
        assert "BLOSUM62" == qresult.param_matrix
        assert 10.0 == qresult.param_evalue_threshold
        assert "L;" == qresult.param_filter
        assert 11 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "blastp" == qresult.program
        assert "nr" == qresult.target

        assert "gi|49176427|ref|NP_418280.3|" == qresult.id
        assert "component of Sec-independent translocase [Escherichia coli K12]" == qresult.description
        assert 103 == qresult.seq_len
        assert 2934173 == qresult.stat_db_num
        assert 1011751523 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 212 == len(qresult)
        # check for alternative ID results
        assert qresult["gi|49176427|ref|NP_418280.3|"] == qresult["gi|26250604|ref|NP_756644.1|"]

        assert qresult["gi|49176427|ref|NP_418280.3|"] == qresult["gi|30064867|ref|NP_839038.1|"]

        assert qresult["gi|49176427|ref|NP_418280.3|"] == qresult["gi|24115132|ref|NP_709642.1|"]

        assert qresult["gi|49176427|ref|NP_418280.3|"] == qresult["gi|24054404|gb|AAN45349.1|"]

        assert qresult["gi|49176427|ref|NP_418280.3|"] == qresult["gi|2367310|gb|AAC76839.1|"]

        assert qresult["gi|49176427|ref|NP_418280.3|"] == qresult["gi|30043127|gb|AAP18849.1|"]

        assert qresult["gi|49176427|ref|NP_418280.3|"] == qresult["gi|26111035|gb|AAN83218.1|"]

        assert qresult["gi|49176427|ref|NP_418280.3|"] == qresult["gi|3193217|gb|AAC19240.1|"]

        assert qresult["gi|49176427|ref|NP_418280.3|"] == qresult["gi|7444818|pir||E65188"]

        hit = qresult[0]
        assert "gi|49176427|ref|NP_418280.3|" == hit.id
        assert "component of Sec-independent translocase [Escherichia coli K12]" == hit.description
        assert 10 == len(hit.id_all)
        assert 10 == len(hit.description_all)
        assert 103 == hit.seq_len
        assert 1 == len(hit)
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        hsp = hit.hsps[0]
        assert 185.267 == hsp.bitscore
        assert 469 == hsp.bitscore_raw
        assert 4.20576e-46 == hsp.evalue
        assert 0 == hsp.query_start
        assert 103 == hsp.query_end
        assert 0 == hsp.hit_start
        assert 103 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 103 == hsp.ident_num
        assert 103 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 103 == hsp.aln_span
        assert "MRLCLIIIYHRGTCMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV" == hsp.query.seq
        assert "MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV" == hsp.hit.seq
        assert "MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV" == hsp.aln_annotation["similarity"]
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        # parse last hit
        hit = qresult[-1]
        assert "gi|39593039|emb|CAE64508.1|" == hit.id
        assert "Hypothetical protein CBG09238 [Caenorhabditis briggsae]" == hit.description
        assert 960 == hit.seq_len
        assert 1 == len(hit)
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        hsp = hit.hsps[0]
        assert 31.5722 == hsp.bitscore
        assert 70 == hsp.bitscore_raw
        assert 7.7721 == hsp.evalue
        assert 54 == hsp.query_start
        assert 102 == hsp.query_end
        assert 409 == hsp.hit_start
        assert 459 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 19 == hsp.ident_num
        assert 33 == hsp.pos_num
        assert 4 == hsp.gap_num
        assert 51 == hsp.aln_span
        assert "KAMSDDEPKQD---KTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ" == hsp.query.seq
        assert "KKEADDKAKKDLEAKTKKEADEKAKKEADEKA-KKEAEAKTKEAEAKTKKE" == hsp.hit.seq
        assert "K  +DD+ K+D   KT ++AD  AK  AD++A   + +AKT++A+   K++" == hsp.aln_annotation["similarity"]
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_xml_2218_blastp_001(self):
        xml_file = get_file("xml_2218_blastp_001.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1

        assert "2.2.18+" == qresult.version
        assert REFERENCE == qresult.reference
        assert "BLOSUM62" == qresult.param_matrix
        assert 10.0 == qresult.param_evalue_threshold
        assert 11 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "blastp" == qresult.program
        assert "nr" == qresult.target

        assert "31493" == qresult.id
        assert "unnamed protein product" == qresult.description
        assert 70 == qresult.seq_len
        assert 15287 == qresult.stat_db_num
        assert 7033566 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 10 == len(qresult)

        hit = qresult[0]
        assert "gi|151942244|gb|EDN60600.1|" == hit.id
        assert "cytosolic iron-sulfur protein assembly protein [Saccharomyces cerevisiae YJM789]" == hit.description
        assert 330 == hit.seq_len
        assert 1 == len(hit)
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        hsp = hit.hsps[0]
        assert 33.113 == hsp.bitscore
        assert 74 == hsp.bitscore_raw
        assert 0.0185319 == hsp.evalue
        assert 14 == hsp.query_start
        assert 62 == hsp.query_end
        assert 113 == hsp.hit_start
        assert 163 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 16 == hsp.ident_num
        assert 27 == hsp.pos_num
        assert 2 == hsp.gap_num
        assert 50 == hsp.aln_span
        assert "AWNKDRTQIAICPNNHEVHIYE--KSGAKWNKVHELKEHNGQVTGIDWAP" == hsp.query.seq
        assert "AWSNDGYYLATCSRDKSVWIWETDESGEEYECISVLQEHSQDVKHVIWHP" == hsp.hit.seq
        assert "AW+ D   +A C  +  V I+E  +SG ++  +  L+EH+  V  + W P" == hsp.aln_annotation["similarity"]
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        # parse last hit
        hit = qresult[-1]
        assert "gi|151567870|pdb|2PM9|B" == hit.id
        assert "Chain B, Crystal Structure Of Yeast Sec1331 VERTEX ELEMENT OF THE Copii Vesicular Coat" == hit.description
        assert 297 == hit.seq_len
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 30.8018 == hsp.bitscore
        assert 68 == hsp.bitscore_raw
        assert 0.0919731 == hsp.evalue
        assert 20 == hsp.query_start
        assert 62 == hsp.query_end
        assert 67 == hsp.hit_start
        assert 109 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 11 == hsp.ident_num
        assert 23 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 42 == hsp.aln_span
        assert "TQIAICPNNHEVHIYEKSGAKWNKVHELKEHNGQVTGIDWAP" == hsp.query.seq
        assert "TILASCSYDGKVMIWKEENGRWSQIAVHAVHSASVNSVQWAP" == hsp.hit.seq
        assert "T +A C  + +V I+++   +W+++     H+  V  + WAP" == hsp.aln_annotation["similarity"]

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_xml_2218_blastp_002(self):
        xml_file = get_file("xml_2218_blastp_002.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1

        # test meta variables, only for the first one
        assert "2.2.18+" == qresult.version
        assert REFERENCE == qresult.reference
        assert "BLOSUM62" == qresult.param_matrix
        assert 0.01 == qresult.param_evalue_threshold
        assert 11 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "m L; R -d repeat/repeat_9606;" == qresult.param_filter
        assert "blastp" == qresult.program
        assert "gpipe/9606/Previous/protein" == qresult.target

        # test parsed values of the first qresult
        assert "gi|585505|sp|Q08386|MOPB_RHOCA" == qresult.id
        assert "Molybdenum-pterin-binding protein mopB >gi|310278|gb|AAA71913.1| molybdenum-pterin-binding protein" == qresult.description
        assert 270 == qresult.seq_len
        assert 27252 == qresult.stat_db_num
        assert 13958303 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 0 == len(qresult)

        # second qresult
        qresult = next(qresults)
        counter += 1

        assert "gi|129628|sp|P07175.1|PARA_AGRTU" == qresult.id
        assert "Protein parA" == qresult.description
        assert 222 == qresult.seq_len
        assert 27252 == qresult.stat_db_num
        assert 13958303 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 0 == len(qresult)

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 2 == counter

    def test_xml_2218L_blastp_001(self):
        xml_file = get_file("xml_2218L_blastp_001.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = next(qresults)
        counter += 1

        # test meta variables, only for the first one
        assert "2.2.18" == qresult.version
        assert ("~Reference: Altschul, Stephen F., "
            "Thomas L. Madden, Alejandro A. Schaffer, "
            "~Jinghui Zhang, Zheng Zhang, Webb Miller, "
            'and David J. Lipman (1997), ~"Gapped BLAST '
            "and PSI-BLAST: a new generation of protein "
            'database search~programs",  '
            "Nucleic Acids Res. 25:3389-3402." == qresult.reference)
        assert "BLOSUM62" == qresult.param_matrix
        assert 1e-05 == qresult.param_evalue_threshold
        assert "F" == qresult.param_filter
        assert 11 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "blastp" == qresult.program
        assert "/Users/pjcock/Downloads/Software/blast-2.2.18/data/nr" == qresult.target

        # test parsed values of the first qresult
        assert "Fake" == qresult.id
        assert "" == qresult.description
        assert 9 == qresult.seq_len
        assert 6589360 == qresult.stat_db_num
        assert 2253133281 == qresult.stat_db_len
        assert 2.02782e10 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 0 == len(qresult)

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_xml_2222_blastp_001(self):
        xml_file = get_file("xml_2222_blastp_001.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1

        assert "2.2.22+" == qresult.version
        assert 'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.' == qresult.reference
        assert "BLOSUM62" == qresult.param_matrix
        assert 1e-06 == qresult.param_evalue_threshold
        assert "F" == qresult.param_filter
        assert 11 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "blastp" == qresult.program
        assert "nr" == qresult.target

        assert "1" == qresult.id
        assert "gi|3298468|dbj|BAA31520.1| SAMIPF" == qresult.description
        assert 107 == qresult.seq_len
        assert 8994603 == qresult.stat_db_num
        assert -1216159329 == qresult.stat_db_len
        assert 76934807744 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 10 == len(qresult)

        hit = qresult[0]
        assert "gi|3298468|dbj|BAA31520.1|" == hit.id
        assert "SAMIPF [Aster tripolium]" == hit.description
        assert 107 == hit.seq_len
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 204.912011757068 == hsp.bitscore
        assert 520 == hsp.bitscore_raw
        assert 1.77242652875017e-51 == hsp.evalue
        assert 0 == hsp.query_start
        assert 107 == hsp.query_end
        assert 0 == hsp.hit_start
        assert 107 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 107 == hsp.ident_num
        assert 107 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 107 == hsp.aln_span
        assert "GGHVNPAVTFGAFVGGNITLLRGIVYIIAQLLGSTVACLLLKFVTNDMAVGVFSLSAGVGVTNALVFEIVMTFGLVYTVYATAIDPKKGSLGTIAPIAIGFIVGANI" == hsp.query.seq
        assert "GGHVNPAVTFGAFVGGNITLLRGIVYIIAQLLGSTVACLLLKFVTNDMAVGVFSLSAGVGVTNALVFEIVMTFGLVYTVYATAIDPKKGSLGTIAPIAIGFIVGANI" == hsp.hit.seq
        assert "GGHVNPAVTFGAFVGGNITLLRGIVYIIAQLLGSTVACLLLKFVTNDMAVGVFSLSAGVGVTNALVFEIVMTFGLVYTVYATAIDPKKGSLGTIAPIAIGFIVGANI" == hsp.aln_annotation["similarity"]
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        # parse last hit
        hit = qresult[-1]
        assert "gi|162809290|dbj|BAF95576.1|" == hit.id
        assert "tonoplast intrinsic protein [Nicotiana tabacum]" == hit.description
        assert 251 == hit.seq_len
        assert 1 == len(hit)
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        hsp = hit.hsps[0]
        assert 177.948041442853 == hsp.bitscore
        assert 450 == hsp.bitscore_raw
        assert 2.0302699895292e-43 == hsp.evalue
        assert 0 == hsp.query_start
        assert 107 == hsp.query_end
        assert 80 == hsp.hit_start
        assert 187 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 91 == hsp.ident_num
        assert 95 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 107 == hsp.aln_span
        assert "GGHVNPAVTFGAFVGGNITLLRGIVYIIAQLLGSTVACLLLKFVTNDMAVGVFSLSAGVGVTNALVFEIVMTFGLVYTVYATAIDPKKGSLGTIAPIAIGFIVGANI" == hsp.query.seq
        assert "GGHVNPAVTFGAFVGGNITLFRGILYIIAQLLGSTVACFLLEFATGGMSTGAFALSAGVSVWNAFVFEIVMTFGLVYTVYATAIDPKKGDLGVIAPIAIGFIVGANI" == hsp.hit.seq
        assert "GGHVNPAVTFGAFVGGNITL RGI+YIIAQLLGSTVAC LL+F T  M+ G F+LSAGV V NA VFEIVMTFGLVYTVYATAIDPKKG LG IAPIAIGFIVGANI" == hsp.aln_annotation["similarity"]
        with pytest.raises(IndexError):
            hit.__getitem__(1)

    # def test_xml_2218L_rpsblast_001(self):
    # PSI-blast, handle later

    def test_xml_2226_blastp_001(self):
        xml_file = get_file("xml_2226_blastp_001.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = next(qresults)
        counter += 1

        # test meta variables, only for the first one
        assert "2.2.26+" == qresult.version
        assert ("Stephen F. Altschul, Thomas L. Madden, Alejandro "
            "A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb "
            'Miller, and David J. Lipman (1997), "Gapped BLAST '
            "and PSI-BLAST: a new generation of protein "
            'database search programs", '
            "Nucleic Acids Res. 25:3389-3402." == qresult.reference)
        assert "BLOSUM62" == qresult.param_matrix
        assert 10.0 == qresult.param_evalue_threshold
        assert "F" == qresult.param_filter
        assert 11 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "blastp" == qresult.program
        assert "db/minirefseq_prot" == qresult.target

        # test parsed values of the first qresult
        assert "random_s00" == qresult.id
        assert "" == qresult.description
        assert 32 == qresult.seq_len
        assert 20 == qresult.stat_db_num
        assert 6406 == qresult.stat_db_len
        assert 156650.0 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 0 == len(qresult)

        # test parsed values of the second qresult
        qresult = next(qresults)
        counter += 1
        assert "gi|16080617|ref|NP_391444.1|" == qresult.id
        assert "membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. 168]" == qresult.description
        assert 102 == qresult.seq_len
        assert 20 == qresult.stat_db_num
        assert 6406 == qresult.stat_db_len
        assert 361344 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "gi|308175296|ref|YP_003922001.1|" == hit.id
        assert "membrane bound lipoprotein [Bacillus amyloliquefaciens DSM 7]" == hit.description
        assert 100 == hit.seq_len
        assert 1 == len(hit)
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        hsp = hit.hsps[0]
        assert 139.428 == hsp.bitscore
        assert 350 == hsp.bitscore_raw
        assert 1.99275e-46 == hsp.evalue
        assert 0 == hsp.query_start
        assert 102 == hsp.query_end
        assert 0 == hsp.hit_start
        assert 100 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 69 == hsp.ident_num
        assert 81 == hsp.pos_num
        assert 2 == hsp.gap_num
        assert 102 == hsp.aln_span
        assert "MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERAN" == hsp.query.seq
        assert "MKKIFGCLFFILLLAGCGVTNEKSQGEDAG--EKLVTKEGTYVGLADTHTIEVTVDHEPVSFDITEESADDVKNLNNGEKVTVKYQKNSKGQLVLKDIEPAN" == hsp.hit.seq
        assert "MKK    LFFILLL+GCGV ++KSQGED      + TKEGTYVGLADTHTIEVTVD+EPVS DITEES  D+   N+G+KVT+ Y+KN +GQL+LKDIE AN" == hsp.aln_annotation["similarity"]
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        # test parsed values of the third qresult
        qresult = next(qresults)
        counter += 1
        assert "gi|11464971:4-101" == qresult.id
        assert "pleckstrin [Mus musculus]" == qresult.description
        assert 98 == qresult.seq_len
        assert 20 == qresult.stat_db_num
        assert 6406 == qresult.stat_db_len
        assert 345626 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "gi|11464971|ref|NP_062422.1|" == hit.id
        assert "pleckstrin [Mus musculus]" == hit.description
        assert 350 == hit.seq_len
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 205.682 == hsp.bitscore
        assert 522 == hsp.bitscore_raw
        assert 2.24956e-69 == hsp.evalue
        assert 0 == hsp.query_start
        assert 98 == hsp.query_end
        assert 3 == hsp.hit_start
        assert 101 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 98 == hsp.ident_num
        assert 98 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 98 == hsp.aln_span
        assert "KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK" == hsp.hit.seq
        assert "KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK" == hsp.query.seq
        assert "KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[-1]
        assert 43.5134 == hsp.bitscore
        assert 101 == hsp.bitscore_raw
        assert 2.90061e-09 == hsp.evalue
        assert 2 == hsp.query_start
        assert 96 == hsp.query_end
        assert 245 == hsp.hit_start
        assert 345 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 29 == hsp.ident_num
        assert 48 == hsp.pos_num
        assert 6 == hsp.gap_num
        assert 100 == hsp.aln_span
        assert "IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNSPKGMIPLKGSTLTS--PCQDFGK--RMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA" == hsp.query.seq
        assert "IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGEDPLGAVHLRGCVVTSVESSHDVKKSDEENLFEIITADEVHYYLQAATSKERTEWIKAIQVA" == hsp.hit.seq
        assert "I++G L+K+G     WK    +L ED   + +Y       P G + L+G  +TS     D  K     + +I T  +  ++ QAA  +ER  W++ I+ A" == hsp.aln_annotation["similarity"]

        hit = qresult[-1]
        assert "gi|350596020|ref|XP_003360649.2|" == hit.id
        assert "PREDICTED: pleckstrin-like [Sus scrofa]" == hit.description
        assert 228 == hit.seq_len
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 199.519 == hsp.bitscore
        assert 506 == hsp.bitscore_raw
        assert 1.97058e-68 == hsp.evalue
        assert 0 == hsp.query_start
        assert 98 == hsp.query_end
        assert 3 == hsp.hit_start
        assert 101 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 94 == hsp.ident_num
        assert 96 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 98 == hsp.aln_span
        assert "KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK" == hsp.query.seq
        assert "KRIREGYLVKKGSMFNTWKPMWVILLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDGWVRDIKKAIK" == hsp.hit.seq
        assert "KRIREGYLVKKGS+FNTWKPMWV+LLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERD WVRDIKKAIK" == hsp.aln_annotation["similarity"]

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 3 == counter

    def test_xml_2226_blastp_002(self):
        xml_file = get_file("xml_2226_blastp_002.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = next(qresults)
        counter += 1

        # test meta variables, only for the first one
        assert "2.2.26+" == qresult.version
        assert ("Stephen F. Altschul, Thomas L. Madden, Alejandro "
            "A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb "
            'Miller, and David J. Lipman (1997), "Gapped BLAST '
            "and PSI-BLAST: a new generation of protein "
            'database search programs", '
            "Nucleic Acids Res. 25:3389-3402." == qresult.reference)
        assert "BLOSUM62" == qresult.param_matrix
        assert 10.0 == qresult.param_evalue_threshold
        assert "F" == qresult.param_filter
        assert 11 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "blastp" == qresult.program
        assert "db/minirefseq_prot" == qresult.target

        # test parsed values of the first qresult
        assert "random_s00" == qresult.id
        assert "" == qresult.description
        assert 32 == qresult.seq_len
        assert 20 == qresult.stat_db_num
        assert 6406 == qresult.stat_db_len
        assert 156650.0 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 0 == len(qresult)

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_xml_2226_blastp_003(self):
        xml_file = get_file("xml_2226_blastp_003.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = next(qresults)
        counter += 1

        assert "2.2.26+" == qresult.version
        assert ("Stephen F. Altschul, Thomas L. Madden, Alejandro "
            "A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb "
            'Miller, and David J. Lipman (1997), "Gapped BLAST '
            "and PSI-BLAST: a new generation of protein "
            'database search programs", '
            "Nucleic Acids Res. 25:3389-3402." == qresult.reference)
        assert "BLOSUM62" == qresult.param_matrix
        assert 10.0 == qresult.param_evalue_threshold
        assert "F" == qresult.param_filter
        assert 11 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "blastp" == qresult.program
        assert "db/minirefseq_prot" == qresult.target

        assert "gi|16080617|ref|NP_391444.1|" == qresult.id
        assert "membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. 168]" == qresult.description
        assert 102 == qresult.seq_len
        assert 20 == qresult.stat_db_num
        assert 6406 == qresult.stat_db_len
        assert 361344 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "gi|308175296|ref|YP_003922001.1|" == hit.id
        assert "membrane bound lipoprotein [Bacillus amyloliquefaciens DSM 7]" == hit.description
        assert 100 == hit.seq_len
        assert 1 == len(hit)
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        hsp = hit.hsps[0]
        assert 139.428 == hsp.bitscore
        assert 350 == hsp.bitscore_raw
        assert 1.99275e-46 == hsp.evalue
        assert 0 == hsp.query_start
        assert 102 == hsp.query_end
        assert 0 == hsp.hit_start
        assert 100 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 69 == hsp.ident_num
        assert 81 == hsp.pos_num
        assert 2 == hsp.gap_num
        assert 102 == hsp.aln_span
        assert "MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERAN" == hsp.query.seq
        assert "MKKIFGCLFFILLLAGCGVTNEKSQGEDAG--EKLVTKEGTYVGLADTHTIEVTVDHEPVSFDITEESADDVKNLNNGEKVTVKYQKNSKGQLVLKDIEPAN" == hsp.hit.seq
        assert "MKK    LFFILLL+GCGV ++KSQGED      + TKEGTYVGLADTHTIEVTVD+EPVS DITEES  D+   N+G+KVT+ Y+KN +GQL+LKDIE AN" == hsp.aln_annotation["similarity"]
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_xml_2226_blastp_004(self):
        xml_file = get_file("xml_2226_blastp_004.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = next(qresults)
        counter += 1

        assert "2.2.26+" == qresult.version
        assert ("Stephen F. Altschul, Thomas L. Madden, Alejandro "
            "A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb "
            'Miller, and David J. Lipman (1997), "Gapped BLAST '
            "and PSI-BLAST: a new generation of protein "
            'database search programs", '
            "Nucleic Acids Res. 25:3389-3402." == qresult.reference)
        assert "BLOSUM62" == qresult.param_matrix
        assert 10.0 == qresult.param_evalue_threshold
        assert "F" == qresult.param_filter
        assert 11 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "blastp" == qresult.program
        assert "db/minirefseq_prot" == qresult.target

        assert "gi|11464971:4-101" == qresult.id
        assert "pleckstrin [Mus musculus]" == qresult.description
        assert 98 == qresult.seq_len
        assert 20 == qresult.stat_db_num
        assert 6406 == qresult.stat_db_len
        assert 345626 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "gi|11464971|ref|NP_062422.1|" == hit.id
        assert "pleckstrin [Mus musculus]" == hit.description
        assert 350 == hit.seq_len
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 205.682 == hsp.bitscore
        assert 522 == hsp.bitscore_raw
        assert 2.24956e-69 == hsp.evalue
        assert 0 == hsp.query_start
        assert 98 == hsp.query_end
        assert 3 == hsp.hit_start
        assert 101 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 98 == hsp.ident_num
        assert 98 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 98 == hsp.aln_span
        assert "KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK" == hsp.hit.seq
        assert "KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK" == hsp.query.seq
        assert "KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[-1]
        assert 43.5134 == hsp.bitscore
        assert 101 == hsp.bitscore_raw
        assert 2.90061e-09 == hsp.evalue
        assert 2 == hsp.query_start
        assert 96 == hsp.query_end
        assert 245 == hsp.hit_start
        assert 345 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 29 == hsp.ident_num
        assert 48 == hsp.pos_num
        assert 6 == hsp.gap_num
        assert 100 == hsp.aln_span
        assert "IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNSPKGMIPLKGSTLTS--PCQDFGK--RMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA" == hsp.query.seq
        assert "IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGEDPLGAVHLRGCVVTSVESSHDVKKSDEENLFEIITADEVHYYLQAATSKERTEWIKAIQVA" == hsp.hit.seq
        assert "I++G L+K+G     WK    +L ED   + +Y       P G + L+G  +TS     D  K     + +I T  +  ++ QAA  +ER  W++ I+ A" == hsp.aln_annotation["similarity"]

        hit = qresult[-1]
        assert "gi|350596020|ref|XP_003360649.2|" == hit.id
        assert "PREDICTED: pleckstrin-like [Sus scrofa]" == hit.description
        assert 228 == hit.seq_len
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 199.519 == hsp.bitscore
        assert 506 == hsp.bitscore_raw
        assert 1.97058e-68 == hsp.evalue
        assert 0 == hsp.query_start
        assert 98 == hsp.query_end
        assert 3 == hsp.hit_start
        assert 101 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 94 == hsp.ident_num
        assert 96 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 98 == hsp.aln_span
        assert "KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK" == hsp.query.seq
        assert "KRIREGYLVKKGSMFNTWKPMWVILLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDGWVRDIKKAIK" == hsp.hit.seq
        assert "KRIREGYLVKKGS+FNTWKPMWV+LLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERD WVRDIKKAIK" == hsp.aln_annotation["similarity"]

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_xml_2226_blastp_005(self):
        xml_file = get_file("xml_2226_blastp_005.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = next(qresults)
        counter += 1

        # test meta variables, only for the first one
        assert "2.2.26+" == qresult.version
        assert ("Stephen F. Altschul, Thomas L. Madden, Alejandro "
            "A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb "
            'Miller, and David J. Lipman (1997), "Gapped BLAST '
            "and PSI-BLAST: a new generation of protein "
            'database search programs", '
            "Nucleic Acids Res. 25:3389-3402." == qresult.reference)
        assert "BLOSUM62" == qresult.param_matrix
        assert 10.0 == qresult.param_evalue_threshold
        assert "F" == qresult.param_filter
        assert 11 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "blastp" == qresult.program
        assert "refseq_protein" == qresult.target

        # test parsed values of the first qresult
        assert "random_s00" == qresult.id
        assert "" == qresult.description
        assert 32 == qresult.seq_len
        assert 12646943 == qresult.stat_db_num
        assert 4397139428 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 0 == len(qresult)

        # test parsed values of the second qresult
        qresult = next(qresults)
        counter += 1
        assert "gi|16080617|ref|NP_391444.1|" == qresult.id
        assert "membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. 168]" == qresult.description
        assert 102 == qresult.seq_len
        assert 12646943 == qresult.stat_db_num
        assert 4397139428 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 5 == len(qresult)
        # check for alternative ID results
        assert qresult["gi|16080617|ref|NP_391444.1|"] == qresult["gi|221311516|ref|ZP_03593363.1|"]

        assert qresult["gi|16080617|ref|NP_391444.1|"] == qresult["gi|221315843|ref|ZP_03597648.1|"]

        assert qresult["gi|16080617|ref|NP_391444.1|"] == qresult["gi|221320757|ref|ZP_03602051.1|"]

        assert qresult["gi|16080617|ref|NP_391444.1|"] == qresult["gi|221325043|ref|ZP_03606337.1|"]

        assert qresult["gi|16080617|ref|NP_391444.1|"] == qresult["gi|321313111|ref|YP_004205398.1|"]

        hit = qresult[0]
        assert "gi|16080617|ref|NP_391444.1|" == hit.id
        assert "membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. 168]" == hit.description
        assert 6 == len(hit.id_all)
        assert 6 == len(hit.description_all)
        assert 102 == hit.seq_len
        assert 1 == len(hit)
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        hsp = hit.hsps[0]
        assert 205.297 == hsp.bitscore
        assert 521 == hsp.bitscore_raw
        assert 1.45285e-66 == hsp.evalue
        assert 0 == hsp.query_start
        assert 102 == hsp.query_end
        assert 0 == hsp.hit_start
        assert 102 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 102 == hsp.ident_num
        assert 102 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 102 == hsp.aln_span
        assert "MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERAN" == hsp.query.seq
        assert "MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERAN" == hsp.hit.seq
        assert "MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERAN" == hsp.aln_annotation["similarity"]
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        # test parsed values of the third qresult
        qresult = next(qresults)
        counter += 1
        assert "gi|11464971:4-101" == qresult.id
        assert "pleckstrin [Mus musculus]" == qresult.description
        assert 98 == qresult.seq_len
        assert 12646943 == qresult.stat_db_num
        assert 4397139428 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "gi|11464971|ref|NP_062422.1|" == hit.id
        assert "pleckstrin [Mus musculus]" == hit.description
        assert 350 == hit.seq_len
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 205.682 == hsp.bitscore
        assert 522 == hsp.bitscore_raw
        assert 1.54412e-63 == hsp.evalue
        assert 0 == hsp.query_start
        assert 98 == hsp.query_end
        assert 3 == hsp.hit_start
        assert 101 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 98 == hsp.ident_num
        assert 98 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 98 == hsp.aln_span
        assert "KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK" == hsp.hit.seq
        assert "KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK" == hsp.query.seq
        assert "KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[-1]
        assert 43.5134 == hsp.bitscore
        assert 101 == hsp.bitscore_raw
        assert 0.00199101 == hsp.evalue
        assert 2 == hsp.query_start
        assert 96 == hsp.query_end
        assert 245 == hsp.hit_start
        assert 345 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 29 == hsp.ident_num
        assert 48 == hsp.pos_num
        assert 6 == hsp.gap_num
        assert 100 == hsp.aln_span
        assert "IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNSPKGMIPLKGSTLTS--PCQDFGK--RMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA" == hsp.query.seq
        assert "IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGEDPLGAVHLRGCVVTSVESSHDVKKSDEENLFEIITADEVHYYLQAATSKERTEWIKAIQVA" == hsp.hit.seq
        assert "I++G L+K+G     WK    +L ED   + +Y       P G + L+G  +TS     D  K     + +I T  +  ++ QAA  +ER  W++ I+ A" == hsp.aln_annotation["similarity"]

        hit = qresult[-1]
        assert "gi|350596020|ref|XP_003360649.2|" == hit.id
        assert "PREDICTED: pleckstrin-like [Sus scrofa]" == hit.description
        assert 228 == hit.seq_len
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 199.519 == hsp.bitscore
        assert 506 == hsp.bitscore_raw
        assert 1.35263e-62 == hsp.evalue
        assert 0 == hsp.query_start
        assert 98 == hsp.query_end
        assert 3 == hsp.hit_start
        assert 101 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 94 == hsp.ident_num
        assert 96 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 98 == hsp.aln_span
        assert "KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK" == hsp.query.seq
        assert "KRIREGYLVKKGSMFNTWKPMWVILLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDGWVRDIKKAIK" == hsp.hit.seq
        assert "KRIREGYLVKKGS+FNTWKPMWV+LLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERD WVRDIKKAIK" == hsp.aln_annotation["similarity"]

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 3 == counter


class BlastxXmlCases(unittest.TestCase):
    def test_xml_2212L_blastx_001(self):
        xml_file = get_file("xml_2212L_blastx_001.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1

        assert "2.2.12" == qresult.version
        assert REFERENCE == qresult.reference
        assert "BLOSUM62" == qresult.param_matrix
        assert 10.0 == qresult.param_evalue_threshold
        assert "L;" == qresult.param_filter
        assert 11 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "blastx" == qresult.program
        assert "nr" == qresult.target

        # test parsed values of the first qresult
        assert "gi|1347369|gb|G25137.1|G25137" == qresult.id
        assert "human STS EST48004, sequence tagged site" == qresult.description
        assert 556 == qresult.seq_len
        assert 2934173 == qresult.stat_db_num
        assert 1011751523 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy

        # test parsed values of the first hit
        hit = qresult[0]
        assert "gi|12654095|gb|AAH00859.1|" == hit.id
        assert "Unknown (protein for IMAGE:3459481) [Homo sapiens]" == hit.description
        assert 319 == hit.seq_len
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 247.284 == hsp.bitscore
        assert 630 == hsp.bitscore_raw
        assert 1.69599e-64 == hsp.evalue
        assert 0 == hsp.query_start
        assert 399 == hsp.query_end
        assert 155 == hsp.hit_start
        assert 288 == hsp.hit_end
        assert 1 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 122 == hsp.ident_num
        assert 123 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 133 == hsp.aln_span
        assert "DLQLLIKAVNLFPAGTNSRWEVIANYMNIHSSSGVKRTAKDVIGKAKSLQKLDPHQKDDINKKAFDKFKKEHGVVPQADNATPSERFXGPYTDFTPXTTEXQKLXEQALNTYPVNTXERWXXIAVAVPGRXKE" == hsp.query.seq
        assert "DLQLLIKAVNLFPAGTNSRWEVIANYMNIHSSSGVKRTAKDVIGKAKSLQKLDPHQKDDINKKAFDKFKKEHGVVPQADNATPSERFEGPYTDFTPWTTEEQKLLEQALKTYPVNTPERWEKIAEAVPGRTKK" == hsp.hit.seq
        assert "DLQLLIKAVNLFPAGTNSRWEVIANYMNIHSSSGVKRTAKDVIGKAKSLQKLDPHQKDDINKKAFDKFKKEHGVVPQADNATPSERF GPYTDFTP TTE QKL EQAL TYPVNT ERW  IA AVPGR K+" == hsp.aln_annotation["similarity"]

        # test parsed values of last hit
        hit = qresult[-1]
        assert "gi|72081091|ref|XP_800619.1|" == hit.id
        assert "PREDICTED: hypothetical protein XP_795526 [Strongylocentrotus purpuratus]" == hit.description
        assert 337 == hit.seq_len

        hsp = hit.hsps[0]
        assert 32.3426 == hsp.bitscore
        assert 72 == hsp.bitscore_raw
        assert 8.57476 == hsp.evalue
        assert 39 == hsp.query_start
        assert 231 == hsp.query_end
        assert 105 == hsp.hit_start
        assert 172 == hsp.hit_end
        assert 1 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 21 == hsp.ident_num
        assert 37 == hsp.pos_num
        assert 3 == hsp.gap_num
        assert "AGTNSRWEVIANYMNI--HSSSGVKRT-AKDVIGKAKSLQKLDPHQKDDINKKAFDKFKKEHGVVPQ" == hsp.query.seq
        assert "SSSNSSSKASASSSNVGASSSSGTKKSDSKSSNESSKSKRDKEDHKEGSINRSKDEKVSKEHRVVKE" == hsp.hit.seq
        assert "+ +NS  +  A+  N+   SSSG K++ +K     +KS +  + H++  IN+   +K  KEH VV +" == hsp.aln_annotation["similarity"]

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_xml_2222_blastx_001(self):
        xml_file = get_file("xml_2222_blastx_001.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1

        assert "2.2.22+" == qresult.version
        assert 'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.' == qresult.reference
        assert "BLOSUM62" == qresult.param_matrix
        assert 0.0001 == qresult.param_evalue_threshold
        assert "L;" == qresult.param_filter
        assert 11 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "blastx" == qresult.program
        assert "nr" == qresult.target

        # test parsed values of the first qresult
        assert "1" == qresult.id
        assert "gi|4104054|gb|AH007193.1|SEG_CVIGS Centaurea vallesiaca 18S ribosomal RNA gene, partial sequence" == qresult.description
        assert 1002 == qresult.seq_len
        assert 8994603 == qresult.stat_db_num
        assert -1216159329 == qresult.stat_db_len
        assert 367397307882 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy

        # test parsed values of the first hit
        hit = qresult[0]
        assert "gi|149390769|gb|ABR25402.1|" == hit.id
        assert "unknown [Oryza sativa (indica cultivar-group)]" == hit.description
        assert 26 == hit.seq_len
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 54.2989775733826 == hsp.bitscore
        assert 129 == hsp.bitscore_raw
        assert 1.83262460293058e-05 == hsp.evalue
        assert 910 == hsp.query_start
        assert 988 == hsp.query_end
        assert 0 == hsp.hit_start
        assert 26 == hsp.hit_end
        assert 2 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 24 == hsp.ident_num
        assert 25 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 26 == hsp.aln_span
        assert "HMLVSKIKPCMCKYEQIQTVKLRMAH" == hsp.query.seq
        assert "HMLVSKIKPCMCKYELIRTVKLRMAH" == hsp.hit.seq
        assert "HMLVSKIKPCMCKYE I+TVKLRMAH" == hsp.aln_annotation["similarity"]

    def test_xml_2226_blastx_001(self):
        xml_file = get_file("xml_2226_blastx_001.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = next(qresults)
        counter += 1

        # test meta variables, only for the first one
        assert "2.2.26+" == qresult.version
        assert ("Stephen F. Altschul, Thomas L. Madden, Alejandro "
            "A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb "
            'Miller, and David J. Lipman (1997), "Gapped BLAST '
            "and PSI-BLAST: a new generation of protein "
            'database search programs", '
            "Nucleic Acids Res. 25:3389-3402." == qresult.reference)
        assert "BLOSUM62" == qresult.param_matrix
        assert 10.0 == qresult.param_evalue_threshold
        assert "L;" == qresult.param_filter
        assert 11 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "blastx" == qresult.program
        assert "db/minirefseq_prot" == qresult.target

        # test parsed values of the first qresult
        assert "random_s00" == qresult.id
        assert "" == qresult.description
        assert 128 == qresult.seq_len
        assert 20 == qresult.stat_db_num
        assert 6406 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert -1 == qresult.stat_kappa
        assert -1 == qresult.stat_lambda
        assert -1 == qresult.stat_entropy
        assert 0 == len(qresult)

        # test parsed values of the third qresult
        qresult = next(qresults)
        counter += 1
        assert "hg19_dna" == qresult.id
        assert "range=chr1:1207057-1207541 5'pad=0 3'pad=0 strand=+ repeatMasking=none" == qresult.description
        assert 485 == qresult.seq_len
        assert 20 == qresult.stat_db_num
        assert 6406 == qresult.stat_db_len
        assert 662354 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "gi|332258565|ref|XP_003278367.1|" == hit.id
        assert "PREDICTED: UPF0764 protein C16orf89-like [Nomascus leucogenys]" == hit.description
        assert 132 == hit.seq_len
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 121.709 == hsp.bitscore
        assert 304 == hsp.bitscore_raw
        assert 2.9522e-38 == hsp.evalue
        assert 15 == hsp.query_start
        assert 300 == hsp.query_end
        assert 24 == hsp.hit_start
        assert 119 == hsp.hit_end
        assert -3 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 69 == hsp.ident_num
        assert 74 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 95 == hsp.aln_span
        assert "LRRSFALVAQAGVQWLDLGXXXXXXPGFK*FSCLSHPSSWDYRHMPPCLINFVFLVETGFYHVGQAGLEPPISGNLPAWASQSVGITGVSHHAQP" == hsp.query.seq
        assert "LRRSFALVAQTRVQWYNLGSPQPPPPGFKRFSCLSLLSSWEYRHVPPHLANFLFLVEMGFLHVGQAGLELVTSGDPPTLTSQSAGIIGVSHCAQP" == hsp.hit.seq
        assert "LRRSFALVAQ  VQW +LG PQPPPPGFK FSCLS  SSW+YRH+PP L NF+FLVE GF HVGQAGLE   SG+ P   SQS GI GVSH AQP" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[-1]
        assert 51.6026 == hsp.bitscore
        assert 122 == hsp.bitscore_raw
        assert 2.73605e-12 == hsp.evalue
        assert 243 == hsp.query_start
        assert 459 == hsp.query_end
        assert 31 == hsp.hit_start
        assert 98 == hsp.hit_end
        assert -3 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 34 == hsp.ident_num
        assert 41 == hsp.pos_num
        assert 5 == hsp.gap_num
        assert 72 == hsp.aln_span
        assert "VGPARVQ*HDLSSLQPPAPEFK*FSHLSLQSSWDCRCPPPHPANXXXXXXXXFLRRSFALVAQAGVQWLDLG" == hsp.query.seq
        assert "VAQTRVQWYNLGSPQPPPPGFKRFSCLSLLSSWEYRHVPPHLAN-----FLFLVEMGFLHVGQAGLELVTSG" == hsp.hit.seq
        assert "V   RVQ ++L S QPP P FK FS LSL SSW+ R  PPH AN     F F +   F  V QAG++ +  G" == hsp.aln_annotation["similarity"]

        hit = qresult[-1]
        assert "gi|33188429|ref|NP_872601.1|" == hit.id
        assert "histone demethylase UTY isoform 1 [Homo sapiens]" == hit.description
        assert 1079 == hit.seq_len
        assert 6 == len(hit)

        hsp = hit.hsps[0]
        assert 104.375 == hsp.bitscore
        assert 259 == hsp.bitscore_raw
        assert 6.31914e-29 == hsp.evalue
        assert 18 == hsp.query_start
        assert 291 == hsp.query_end
        assert 988 == hsp.hit_start
        assert 1079 == hsp.hit_end
        assert -3 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 59 == hsp.ident_num
        assert 66 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 91 == hsp.aln_span
        assert "SFALVAQAGVQWLDLGXXXXXXPGFK*FSCLSHPSSWDYRHMPPCLINFVFLVETGFYHVGQAGLEPPISGNLPAWASQSVGITGVSHHAQ" == hsp.query.seq
        assert "SFQESLRAGMQWCDLSSLQPPPPGFKRFSHLSLPNSWNYRHLPSCPTNFCIFVETGFHHVGQACLELLTSGGLLASASQSAGITGVSHHAR" == hsp.hit.seq
        assert "SF    +AG+QW DL   QPPPPGFK FS LS P+SW+YRH+P C  NF   VETGF+HVGQA LE   SG L A ASQS GITGVSHHA+" == hsp.aln_annotation["similarity"]

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 2 == counter

    def test_xml_2226_blastx_002(self):
        xml_file = get_file("xml_2226_blastx_002.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1

        assert "2.2.26+" == qresult.version
        assert ("Stephen F. Altschul, Thomas L. Madden, Alejandro "
            "A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb "
            'Miller, and David J. Lipman (1997), "Gapped BLAST '
            "and PSI-BLAST: a new generation of protein database "
            'search programs", Nucleic Acids Res. 25:3389-3402.' == qresult.reference)
        assert "BLOSUM62" == qresult.param_matrix
        assert 10.0 == qresult.param_evalue_threshold
        assert "L;" == qresult.param_filter
        assert 11 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "blastx" == qresult.program
        assert "db/minirefseq_prot" == qresult.target

        assert "random_s00" == qresult.id
        assert "" == qresult.description
        assert 128 == qresult.seq_len
        assert 20 == qresult.stat_db_num
        assert 6406 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert -1 == qresult.stat_kappa
        assert -1 == qresult.stat_lambda
        assert -1 == qresult.stat_entropy
        assert 0 == len(qresult)

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_xml_2226_blastx_003(self):
        xml_file = get_file("xml_2226_blastx_003.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1

        assert "2.2.26+" == qresult.version
        assert ("Stephen F. Altschul, Thomas L. Madden, Alejandro "
            "A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb "
            'Miller, and David J. Lipman (1997), "Gapped BLAST '
            "and PSI-BLAST: a new generation of protein database "
            'search programs", Nucleic Acids Res. 25:3389-3402.' == qresult.reference)
        assert "BLOSUM62" == qresult.param_matrix
        assert 10.0 == qresult.param_evalue_threshold
        assert "L;" == qresult.param_filter
        assert 11 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "blastx" == qresult.program
        assert "db/minirefseq_prot" == qresult.target

        assert "hg19_dna" == qresult.id
        assert "range=chr1:1207057-1207541 5'pad=0 3'pad=0 strand=+ repeatMasking=none" == qresult.description
        assert 485 == qresult.seq_len
        assert 20 == qresult.stat_db_num
        assert 6406 == qresult.stat_db_len
        assert 662354 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "gi|332258565|ref|XP_003278367.1|" == hit.id
        assert "PREDICTED: UPF0764 protein C16orf89-like [Nomascus leucogenys]" == hit.description
        assert 132 == hit.seq_len
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 121.709 == hsp.bitscore
        assert 304 == hsp.bitscore_raw
        assert 2.9522e-38 == hsp.evalue
        assert 15 == hsp.query_start
        assert 300 == hsp.query_end
        assert 24 == hsp.hit_start
        assert 119 == hsp.hit_end
        assert -3 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 69 == hsp.ident_num
        assert 74 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 95 == hsp.aln_span
        assert "LRRSFALVAQAGVQWLDLGXXXXXXPGFK*FSCLSHPSSWDYRHMPPCLINFVFLVETGFYHVGQAGLEPPISGNLPAWASQSVGITGVSHHAQP" == hsp.query.seq
        assert "LRRSFALVAQTRVQWYNLGSPQPPPPGFKRFSCLSLLSSWEYRHVPPHLANFLFLVEMGFLHVGQAGLELVTSGDPPTLTSQSAGIIGVSHCAQP" == hsp.hit.seq
        assert "LRRSFALVAQ  VQW +LG PQPPPPGFK FSCLS  SSW+YRH+PP L NF+FLVE GF HVGQAGLE   SG+ P   SQS GI GVSH AQP" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[-1]
        assert 51.6026 == hsp.bitscore
        assert 122 == hsp.bitscore_raw
        assert 2.73605e-12 == hsp.evalue
        assert 243 == hsp.query_start
        assert 459 == hsp.query_end
        assert 31 == hsp.hit_start
        assert 98 == hsp.hit_end
        assert -3 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 34 == hsp.ident_num
        assert 41 == hsp.pos_num
        assert 5 == hsp.gap_num
        assert 72 == hsp.aln_span
        assert "VGPARVQ*HDLSSLQPPAPEFK*FSHLSLQSSWDCRCPPPHPANXXXXXXXXFLRRSFALVAQAGVQWLDLG" == hsp.query.seq
        assert "VAQTRVQWYNLGSPQPPPPGFKRFSCLSLLSSWEYRHVPPHLAN-----FLFLVEMGFLHVGQAGLELVTSG" == hsp.hit.seq
        assert "V   RVQ ++L S QPP P FK FS LSL SSW+ R  PPH AN     F F +   F  V QAG++ +  G" == hsp.aln_annotation["similarity"]

        hit = qresult[-1]
        assert "gi|33188429|ref|NP_872601.1|" == hit.id
        assert "histone demethylase UTY isoform 1 [Homo sapiens]" == hit.description
        assert 1079 == hit.seq_len
        assert 6 == len(hit)

        hsp = hit.hsps[0]
        assert 104.375 == hsp.bitscore
        assert 259 == hsp.bitscore_raw
        assert 6.31914e-29 == hsp.evalue
        assert 18 == hsp.query_start
        assert 291 == hsp.query_end
        assert 988 == hsp.hit_start
        assert 1079 == hsp.hit_end
        assert -3 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 59 == hsp.ident_num
        assert 66 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 91 == hsp.aln_span
        assert "SFALVAQAGVQWLDLGXXXXXXPGFK*FSCLSHPSSWDYRHMPPCLINFVFLVETGFYHVGQAGLEPPISGNLPAWASQSVGITGVSHHAQ" == hsp.query.seq
        assert "SFQESLRAGMQWCDLSSLQPPPPGFKRFSHLSLPNSWNYRHLPSCPTNFCIFVETGFHHVGQACLELLTSGGLLASASQSAGITGVSHHAR" == hsp.hit.seq
        assert "SF    +AG+QW DL   QPPPPGFK FS LS P+SW+YRH+P C  NF   VETGF+HVGQA LE   SG L A ASQS GITGVSHHA+" == hsp.aln_annotation["similarity"]

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_xml_2226_blastx_004(self):
        xml_file = get_file("xml_2226_blastx_004.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = next(qresults)
        counter += 1

        # test meta variables, only for the first one
        assert "2.2.26+" == qresult.version
        assert ("Stephen F. Altschul, Thomas L. Madden, Alejandro "
            "A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb "
            'Miller, and David J. Lipman (1997), "Gapped BLAST '
            "and PSI-BLAST: a new generation of protein database "
            'search programs", Nucleic Acids Res. 25:3389-3402.' == qresult.reference)
        assert "BLOSUM62" == qresult.param_matrix
        assert 10.0 == qresult.param_evalue_threshold
        assert "L;" == qresult.param_filter
        assert 11 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "blastx" == qresult.program
        assert "refseq_protein" == qresult.target

        # test parsed values of the first qresult
        assert "random_s00" == qresult.id
        assert "" == qresult.description
        assert 128 == qresult.seq_len
        assert 12646943 == qresult.stat_db_num
        assert 4397139428 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 0 == len(qresult)

        # test parsed values of the second qresult
        qresult = next(qresults)
        counter += 1
        assert "hg19_dna" == qresult.id
        assert "range=chr1:1207057-1207541 5'pad=0 3'pad=0 strand=+ repeatMasking=none" == qresult.description
        assert 485 == qresult.seq_len
        assert 12646943 == qresult.stat_db_num
        assert 4397139428 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "gi|332258565|ref|XP_003278367.1|" == hit.id
        assert "PREDICTED: UPF0764 protein C16orf89-like [Nomascus leucogenys]" == hit.description
        assert 132 == hit.seq_len
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 121.709 == hsp.bitscore
        assert 304 == hsp.bitscore_raw
        assert 2.02642e-32 == hsp.evalue
        assert 15 == hsp.query_start
        assert 300 == hsp.query_end
        assert 24 == hsp.hit_start
        assert 119 == hsp.hit_end
        assert -3 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 69 == hsp.ident_num
        assert 74 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 95 == hsp.aln_span
        assert "LRRSFALVAQAGVQWLDLGXXXXXXPGFK*FSCLSHPSSWDYRHMPPCLINFVFLVETGFYHVGQAGLEPPISGNLPAWASQSVGITGVSHHAQP" == hsp.query.seq
        assert "LRRSFALVAQTRVQWYNLGSPQPPPPGFKRFSCLSLLSSWEYRHVPPHLANFLFLVEMGFLHVGQAGLELVTSGDPPTLTSQSAGIIGVSHCAQP" == hsp.hit.seq
        assert "LRRSFALVAQ  VQW +LG PQPPPPGFK FSCLS  SSW+YRH+PP L NF+FLVE GF HVGQAGLE   SG+ P   SQS GI GVSH AQP" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[-1]
        assert 51.6026 == hsp.bitscore
        assert 122 == hsp.bitscore_raw
        assert 1.87805e-06 == hsp.evalue
        assert 243 == hsp.query_start
        assert 459 == hsp.query_end
        assert 31 == hsp.hit_start
        assert 98 == hsp.hit_end
        assert -3 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 34 == hsp.ident_num
        assert 41 == hsp.pos_num
        assert 5 == hsp.gap_num
        assert 72 == hsp.aln_span
        assert "VGPARVQ*HDLSSLQPPAPEFK*FSHLSLQSSWDCRCPPPHPANXXXXXXXXFLRRSFALVAQAGVQWLDLG" == hsp.query.seq
        assert "VAQTRVQWYNLGSPQPPPPGFKRFSCLSLLSSWEYRHVPPHLAN-----FLFLVEMGFLHVGQAGLELVTSG" == hsp.hit.seq
        assert "V   RVQ ++L S QPP P FK FS LSL SSW+ R  PPH AN     F F +   F  V QAG++ +  G" == hsp.aln_annotation["similarity"]

        hit = qresult[-1]
        assert "gi|332815399|ref|XP_003309509.1|" == hit.id
        assert "PREDICTED: histone demethylase UTY-like [Pan troglodytes]" == hit.description
        assert 101 == hit.seq_len
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 97.0561 == hsp.bitscore
        assert 240 == hsp.bitscore_raw
        assert 2.76414e-23 == hsp.evalue
        assert 6 == hsp.query_start
        assert 279 == hsp.query_end
        assert 9 == hsp.hit_start
        assert 100 == hsp.hit_end
        assert -3 == hsp.query_frame
        assert 0 == hsp.hit_frame
        assert 56 == hsp.ident_num
        assert 62 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 91 == hsp.aln_span
        assert "VAQAGVQWLDLGXXXXXXPGFK*FSCLSHPSSWDYRHMPPCLINFVFLVETGFYHVGQAGLEPPISGNLPAWASQSVGITGVSHHAQPLCE" == hsp.query.seq
        assert "VPHAGVQWHNLSSLQPPPSRFKPFSYLSLLSSWDQRRPPPCLVTFVFLIETGFRHVGQAGLKLLTSGDPSASASQSAGIRGVSHCTWPECQ" == hsp.hit.seq
        assert "V  AGVQW +L   QPPP  FK FS LS  SSWD R  PPCL+ FVFL+ETGF HVGQAGL+   SG+  A ASQS GI GVSH   P C+" == hsp.aln_annotation["similarity"]

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 2 == counter


class TblastnXmlCases(unittest.TestCase):
    def test_xml_2212L_tblastn_001(self):
        xml_file = get_file("xml_2212L_tblastn_001.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1

        assert "2.2.12" == qresult.version
        assert REFERENCE == qresult.reference
        assert "BLOSUM62" == qresult.param_matrix
        assert 0.001 == qresult.param_evalue_threshold
        assert 11 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "L;" == qresult.param_filter
        assert "tblastn" == qresult.program
        assert "nr" == qresult.target

        # test parsed values of qresult
        assert "gi|729325|sp|P39483|DHG2_BACME" == qresult.id
        assert "Glucose 1-dehydrogenase II (GLCDH-II)" == qresult.description
        assert 261 == qresult.seq_len
        assert 251887 == qresult.stat_db_num
        assert 438542399 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 100 == len(qresult)

        hit = qresult[0]
        assert "gi|58264321|ref|XM_569317.1|" == hit.id
        assert "Filobasidiella neoformans glucose 1-dehydrogenase, putative (CNB05760) mRNA, complete cds" == hit.description
        assert 904 == hit.seq_len
        assert 1 == len(hit)
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        hsp = hit.hsps[0]
        assert 148.288 == hsp.bitscore
        assert 373 == hsp.bitscore_raw
        assert 1.46834e-35 == hsp.evalue
        assert 4 == hsp.query_start
        assert 250 == hsp.query_end
        assert 15 == hsp.hit_start
        assert 762 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 1 == hsp.hit_frame
        assert 84 == hsp.ident_num
        assert 143 == hsp.pos_num
        assert 9 == hsp.gap_num
        assert 252 == hsp.aln_span
        assert "LKDKVVVVTGGSKGLGRAMAVRFGQEQSKVVVNYRSNXXXXXXXXXXXXXXXXGGQAIIVRGDVTKEEDVVNLVETAVKEFGSLDVMINNAGVENPVPSH---ELSLENWNQVIDTNLTGAFLGSREAIKYFVENDIKG-NVINMSSVHEMIPWPLFVHYAASKGGMKLMTETLALEYAPKGIRVNNIGPGAIDTPINAEKFADPEQRADVESMIPMGYIGKPEEIASVAAFLASSQASYVTGITLFADGGM" == hsp.query.seq
        assert "LQGKVVAITGCSTGIGRAIAIGAAKNGANVVLHHLGDSTASDIAQVQEECKQAGAKTVVVPGDIAEAKTANEIVSAAVSSFSRIDVLISNAGI---CPFHSFLDLPHPLWKRVQDVNLNGSFYVVQAVANQMAKQEPKGGSIVAVSSISALMGGGEQCHYTPTKAGIKSLMESCAIALGPMGIRCNSVLPGTIETNINKEDLSNPEKRADQIRRVPLGRLGKPEDLVGPTLFFASDLSNYCTGASVLVDGGM" == hsp.hit.seq
        assert "L+ KVV +TG S G+GRA+A+   +  + VV+++  +   +   + + E  +AG + ++V GD+ + +    +V  AV  F  +DV+I+NAG+    P H   +L    W +V D NL G+F   +       + + KG +++ +SS+  ++      HY  +K G+K + E+ A+   P GIR N++ PG I+T IN E  ++PE+RAD    +P+G +GKPE++     F AS  ++Y TG ++  DGGM" == hsp.aln_annotation["similarity"]

        # parse last hit
        hit = qresult[-1]
        assert "gi|450259|gb|L27825.1|EMEVERA1AA" == hit.id
        assert "Emericella nidulans (verA) gene, complete cds, ORF 1 gene, complete cds, and ORF 2 gene, 5' end" == hit.description
        assert 4310 == hit.seq_len
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 91.2781 == hsp.bitscore
        assert 225 == hsp.bitscore_raw
        assert 1.31998e-20 == hsp.evalue
        assert 4 == hsp.query_start
        assert 204 == hsp.query_end
        assert 578 == hsp.hit_start
        assert 1253 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 3 == hsp.hit_frame
        assert 76 == hsp.ident_num
        assert 113 == hsp.pos_num
        assert 31 == hsp.gap_num
        assert 228 == hsp.aln_span
        assert "LKDKVVVVTGGSKGLGRAMAVRFGQEQSKVVVNYRSNXXXXXXXXXXXXXXGGQAIIVRGDVTKEEDVVNLVETAVKEFGSLDVMINNAGV-----------ENPVPS-HELSLE-----NWNQVIDTNLTGAFLGSREAIKYFVENDIKGNVINMSSVHEMIPW-PLFVHYAASKGGMKLMTETLALEYAPKGIRVNNIGPGAIDTPI----------NAEKFADPE" == hsp.query.seq
        assert "LDGKVALVTGAGRGIGAAIAVALGQPGAKVVVNYANSREAAEKVVDEIKSNAQSAISIQADVGDPDAVTKLMDQAVEHFGYLDIVSSNAGIVSFGHVKDVTPDVCVPSPYESPVEL*PQQEFDRVFRVNTRGQFFVAREAYRHLREG---GRIILTSSNTASVKGVPRHAVYSGSKGAIDTFVRCLAIDCGDKKITVNAVAPGAIKTDMFLSVSREYIPNGETFTDEQ" == hsp.hit.seq
        assert "L  KV +VTG  +G+G A+AV  GQ  +KVVVNY ++ E A +V  EI+     AI ++ DV   + V  L++ AV+ FG LD++ +NAG+           +  VPS +E  +E      +++V   N  G F  +REA ++  E    G +I  SS    +   P    Y+ SKG +      LA++   K I VN + PGAI T +          N E F D +" == hsp.aln_annotation["similarity"]

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_xml_2226_tblastn_001(self):
        xml_file = get_file("xml_2226_tblastn_001.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = next(qresults)
        counter += 1

        # test meta variables, only for the first one
        assert "2.2.26+" == qresult.version
        assert ("Stephen F. Altschul, Thomas L. Madden, Alejandro "
            "A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, "
            "Webb Miller, and David J. Lipman (1997), "
            '"Gapped BLAST and PSI-BLAST: a new generation of '
            'protein database search programs", '
            "Nucleic Acids Res. 25:3389-3402." == qresult.reference)
        assert "BLOSUM62" == qresult.param_matrix
        assert 10.0 == qresult.param_evalue_threshold
        assert "L;" == qresult.param_filter
        assert 11 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "tblastn" == qresult.program
        assert "db/minirefseq_mrna" == qresult.target

        # test parsed values of the first qresult
        assert "random_s00" == qresult.id
        assert "" == qresult.description
        assert 32 == qresult.seq_len
        assert 23 == qresult.stat_db_num
        assert 67750 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert -1 == qresult.stat_kappa
        assert -1 == qresult.stat_lambda
        assert -1 == qresult.stat_entropy
        assert 0 == len(qresult)

        # test parsed values of the second qresult
        qresult = next(qresults)
        counter += 1
        assert "gi|16080617|ref|NP_391444.1|" == qresult.id
        assert "membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. 168]" == qresult.description
        assert 102 == qresult.seq_len
        assert 23 == qresult.stat_db_num
        assert 67750 == qresult.stat_db_len
        assert 1205400.0 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 3 == len(qresult)

        hit = qresult[0]
        assert "gi|145479850|ref|XM_001425911.1|" == hit.id
        assert "Paramecium tetraurelia hypothetical protein (GSPATT00004923001) partial mRNA" == hit.description
        assert 4632 == hit.seq_len
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 34.6538 == hsp.bitscore
        assert 78 == hsp.bitscore_raw
        assert 1.08241e-05 == hsp.evalue
        assert 30 == hsp.query_start
        assert 73 == hsp.query_end
        assert 1743 == hsp.hit_start
        assert 1872 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 1 == hsp.hit_frame
        assert 15 == hsp.ident_num
        assert 26 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 43 == hsp.aln_span
        assert "PDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLD" == hsp.query.seq
        assert "PKTATGTKKGTIIGLLSIHTILFILTSHALSLEVKEQT*KDID" == hsp.hit.seq
        assert "P +   TK+GT +GL   HTI   + +  +SL++ E++  D+D" == hsp.aln_annotation["similarity"]
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        # test parsed values of the third qresult
        qresult = next(qresults)
        counter += 1
        assert "gi|11464971:4-101" == qresult.id
        assert "pleckstrin [Mus musculus]" == qresult.description
        assert 98 == qresult.seq_len
        assert 23 == qresult.stat_db_num
        assert 67750 == qresult.stat_db_len
        assert 1119300.0 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "gi|350596019|ref|XM_003360601.2|" == hit.id
        assert "PREDICTED: Sus scrofa pleckstrin-like (LOC100626968), mRNA" == hit.description
        assert 772 == hit.seq_len
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 199.519 == hsp.bitscore
        assert 506 == hsp.bitscore_raw
        assert 1.57249e-67 == hsp.evalue
        assert 0 == hsp.query_start
        assert 98 == hsp.query_end
        assert 94 == hsp.hit_start
        assert 388 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 2 == hsp.hit_frame
        assert 94 == hsp.ident_num
        assert 96 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 98 == hsp.aln_span
        assert "KRIREGYLVKKGSMFNTWKPMWVILLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDGWVRDIKKAIK" == hsp.hit.seq
        assert "KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK" == hsp.query.seq
        assert "KRIREGYLVKKGS+FNTWKPMWV+LLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERD WVRDIKKAIK" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[-1]
        assert 32.7278 == hsp.bitscore
        assert 73 == hsp.bitscore_raw
        assert 4.07518e-05 == hsp.evalue
        assert 29 == hsp.query_start
        assert 96 == hsp.query_end
        assert 541 == hsp.hit_start
        assert 754 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 2 == hsp.hit_frame
        assert 21 == hsp.ident_num
        assert 33 == hsp.pos_num
        assert 4 == hsp.gap_num
        assert 71 == hsp.aln_span
        assert "IEFYKKKSDNSPKGMIPLKGSTLTS-PCQDFGKRMFVLK---ITTTKQQDHFFQAAFLEERDAWVRDIKKA" == hsp.query.seq
        assert "LHYYDPAGGEDPLGAIHLRGCVVTSVESNTDGKNGFLWERAXXITADEVHYFLQAANPKERTEWIKAIQVA" == hsp.hit.seq
        assert "+ +Y       P G I L+G  +TS      GK  F+ +     T  +  +F QAA  +ER  W++ I+ A" == hsp.aln_annotation["similarity"]

        hit = qresult[-1]
        assert "gi|365982352|ref|XM_003667962.1|" == hit.id
        assert "Naumovozyma dairenensis CBS 421 hypothetical protein (NDAI0A06120), mRNA" == hit.description
        assert 4932 == hit.seq_len
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 19.631 == hsp.bitscore
        assert 39 == hsp.bitscore_raw
        assert 1.65923 == hsp.evalue
        assert 11 == hsp.query_start
        assert 54 == hsp.query_end
        assert 3180 == hsp.hit_start
        assert 3336 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 1 == hsp.hit_frame
        assert 16 == hsp.ident_num
        assert 23 == hsp.pos_num
        assert 9 == hsp.gap_num
        assert 52 == hsp.aln_span
        assert "GSVFNTWKPMWVVLL---------EDGIEFYKKKSDNSPKGMIPLKGSTLTS" == hsp.query.seq
        assert "GSCFPTWDLIFIEVLNPFLKEKLWEADNEEISKFVDLTLKGLVDLYPSHFTS" == hsp.hit.seq
        assert "GS F TW  +++ +L         E   E   K  D + KG++ L  S  TS" == hsp.aln_annotation["similarity"]

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 3 == counter

    def test_xml_2226_tblastn_002(self):
        xml_file = get_file("xml_2226_tblastn_002.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1

        assert "2.2.26+" == qresult.version
        assert ("Stephen F. Altschul, Thomas L. Madden, Alejandro "
            "A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, "
            "Webb Miller, and David J. Lipman (1997), "
            '"Gapped BLAST and PSI-BLAST: a new generation of '
            'protein database search programs", '
            "Nucleic Acids Res. 25:3389-3402." == qresult.reference)
        assert "BLOSUM62" == qresult.param_matrix
        assert 10.0 == qresult.param_evalue_threshold
        assert "L;" == qresult.param_filter
        assert 11 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "tblastn" == qresult.program
        assert "db/minirefseq_mrna" == qresult.target

        assert "random_s00" == qresult.id
        assert "" == qresult.description
        assert 32 == qresult.seq_len
        assert 23 == qresult.stat_db_num
        assert 67750 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert -1 == qresult.stat_kappa
        assert -1 == qresult.stat_lambda
        assert -1 == qresult.stat_entropy
        assert 0 == len(qresult)

        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_xml_2226_tblastn_003(self):
        xml_file = get_file("xml_2226_tblastn_003.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1

        assert "2.2.26+" == qresult.version
        assert ("Stephen F. Altschul, Thomas L. Madden, Alejandro "
            "A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, "
            "Webb Miller, and David J. Lipman (1997), "
            '"Gapped BLAST and PSI-BLAST: a new generation of '
            'protein database search programs", '
            "Nucleic Acids Res. 25:3389-3402." == qresult.reference)
        assert "BLOSUM62" == qresult.param_matrix
        assert 10.0 == qresult.param_evalue_threshold
        assert "L;" == qresult.param_filter
        assert 11 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "tblastn" == qresult.program
        assert "db/minirefseq_mrna" == qresult.target

        assert "gi|16080617|ref|NP_391444.1|" == qresult.id
        assert "membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. 168]" == qresult.description
        assert 102 == qresult.seq_len
        assert 23 == qresult.stat_db_num
        assert 67750 == qresult.stat_db_len
        assert 1205400.0 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 3 == len(qresult)

        hit = qresult[0]
        assert "gi|145479850|ref|XM_001425911.1|" == hit.id
        assert "Paramecium tetraurelia hypothetical protein (GSPATT00004923001) partial mRNA" == hit.description
        assert 4632 == hit.seq_len
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 34.6538 == hsp.bitscore
        assert 78 == hsp.bitscore_raw
        assert 1.08241e-05 == hsp.evalue
        assert 30 == hsp.query_start
        assert 73 == hsp.query_end
        assert 1743 == hsp.hit_start
        assert 1872 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 1 == hsp.hit_frame
        assert 15 == hsp.ident_num
        assert 26 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 43 == hsp.aln_span
        assert "PDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLD" == hsp.query.seq
        assert "PKTATGTKKGTIIGLLSIHTILFILTSHALSLEVKEQT*KDID" == hsp.hit.seq
        assert "P +   TK+GT +GL   HTI   + +  +SL++ E++  D+D" == hsp.aln_annotation["similarity"]

        with pytest.raises(IndexError):
            hit.__getitem__(1)

        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_xml_2226_tblastn_004(self):
        xml_file = get_file("xml_2226_tblastn_004.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1

        assert "2.2.26+" == qresult.version
        assert ("Stephen F. Altschul, Thomas L. Madden, Alejandro "
            "A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, "
            "Webb Miller, and David J. Lipman (1997), "
            '"Gapped BLAST and PSI-BLAST: a new generation of '
            'protein database search programs", '
            "Nucleic Acids Res. 25:3389-3402." == qresult.reference)
        assert "BLOSUM62" == qresult.param_matrix
        assert 10.0 == qresult.param_evalue_threshold
        assert "L;" == qresult.param_filter
        assert 11 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "tblastn" == qresult.program
        assert "db/minirefseq_mrna" == qresult.target

        assert "gi|11464971:4-101" == qresult.id
        assert "pleckstrin [Mus musculus]" == qresult.description
        assert 98 == qresult.seq_len
        assert 23 == qresult.stat_db_num
        assert 67750 == qresult.stat_db_len
        assert 1119300.0 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "gi|350596019|ref|XM_003360601.2|" == hit.id
        assert "PREDICTED: Sus scrofa pleckstrin-like (LOC100626968), mRNA" == hit.description
        assert 772 == hit.seq_len
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 199.519 == hsp.bitscore
        assert 506 == hsp.bitscore_raw
        assert 1.57249e-67 == hsp.evalue
        assert 0 == hsp.query_start
        assert 98 == hsp.query_end
        assert 94 == hsp.hit_start
        assert 388 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 2 == hsp.hit_frame
        assert 94 == hsp.ident_num
        assert 96 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 98 == hsp.aln_span
        assert "KRIREGYLVKKGSMFNTWKPMWVILLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDGWVRDIKKAIK" == hsp.hit.seq
        assert "KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK" == hsp.query.seq
        assert "KRIREGYLVKKGS+FNTWKPMWV+LLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERD WVRDIKKAIK" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[-1]
        assert 32.7278 == hsp.bitscore
        assert 73 == hsp.bitscore_raw
        assert 4.07518e-05 == hsp.evalue
        assert 29 == hsp.query_start
        assert 96 == hsp.query_end
        assert 541 == hsp.hit_start
        assert 754 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 2 == hsp.hit_frame
        assert 21 == hsp.ident_num
        assert 33 == hsp.pos_num
        assert 4 == hsp.gap_num
        assert 71 == hsp.aln_span
        assert "IEFYKKKSDNSPKGMIPLKGSTLTS-PCQDFGKRMFVLK---ITTTKQQDHFFQAAFLEERDAWVRDIKKA" == hsp.query.seq
        assert "LHYYDPAGGEDPLGAIHLRGCVVTSVESNTDGKNGFLWERAXXITADEVHYFLQAANPKERTEWIKAIQVA" == hsp.hit.seq
        assert "+ +Y       P G I L+G  +TS      GK  F+ +     T  +  +F QAA  +ER  W++ I+ A" == hsp.aln_annotation["similarity"]

        hit = qresult[-1]
        assert "gi|365982352|ref|XM_003667962.1|" == hit.id
        assert "Naumovozyma dairenensis CBS 421 hypothetical protein (NDAI0A06120), mRNA" == hit.description
        assert 4932 == hit.seq_len
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 19.631 == hsp.bitscore
        assert 39 == hsp.bitscore_raw
        assert 1.65923 == hsp.evalue
        assert 11 == hsp.query_start
        assert 54 == hsp.query_end
        assert 3180 == hsp.hit_start
        assert 3336 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 1 == hsp.hit_frame
        assert 16 == hsp.ident_num
        assert 23 == hsp.pos_num
        assert 9 == hsp.gap_num
        assert 52 == hsp.aln_span
        assert "GSVFNTWKPMWVVLL---------EDGIEFYKKKSDNSPKGMIPLKGSTLTS" == hsp.query.seq
        assert "GSCFPTWDLIFIEVLNPFLKEKLWEADNEEISKFVDLTLKGLVDLYPSHFTS" == hsp.hit.seq
        assert "GS F TW  +++ +L         E   E   K  D + KG++ L  S  TS" == hsp.aln_annotation["similarity"]

        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_xml_2226_tblastn_005(self):
        xml_file = get_file("xml_2226_tblastn_005.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = next(qresults)
        counter += 1

        # test meta variables, only for the first one
        assert "2.2.26+" == qresult.version
        assert ("Stephen F. Altschul, Thomas L. Madden, Alejandro "
            "A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, "
            "Webb Miller, and David J. Lipman (1997), "
            '"Gapped BLAST and PSI-BLAST: a new generation of '
            'protein database search programs", '
            "Nucleic Acids Res. 25:3389-3402." == qresult.reference)
        assert "BLOSUM62" == qresult.param_matrix
        assert 10.0 == qresult.param_evalue_threshold
        assert "L;" == qresult.param_filter
        assert 11 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "tblastn" == qresult.program
        assert "refseq_rna" == qresult.target

        # test parsed values of the first qresult
        assert "random_s00" == qresult.id
        assert "" == qresult.description
        assert 32 == qresult.seq_len
        assert 2933984 == qresult.stat_db_num
        assert 4726730735 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 0 == len(qresult)

        # test parsed values of the second qresult
        qresult = next(qresults)
        counter += 1
        assert "gi|16080617|ref|NP_391444.1|" == qresult.id
        assert "membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. 168]" == qresult.description
        assert 102 == qresult.seq_len
        assert 2933984 == qresult.stat_db_num
        assert 4726730735 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 3 == len(qresult)

        hit = qresult[0]
        assert "gi|145479850|ref|XM_001425911.1|" == hit.id
        assert "Paramecium tetraurelia hypothetical protein (GSPATT00004923001) partial mRNA" == hit.description
        assert 4632 == hit.seq_len
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 34.6538 == hsp.bitscore
        assert 78 == hsp.bitscore_raw
        assert 0.755176 == hsp.evalue
        assert 30 == hsp.query_start
        assert 73 == hsp.query_end
        assert 1743 == hsp.hit_start
        assert 1872 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 1 == hsp.hit_frame
        assert 15 == hsp.ident_num
        assert 26 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 43 == hsp.aln_span
        assert "PDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLD" == hsp.query.seq
        assert "PKTATGTKKGTIIGLLSIHTILFILTSHALSLEVKEQT*KDID" == hsp.hit.seq
        assert "P +   TK+GT +GL   HTI   + +  +SL++ E++  D+D" == hsp.aln_annotation["similarity"]
        with pytest.raises(IndexError):
            hit.__getitem__(1)

        # test parsed values of the third qresult
        qresult = next(qresults)
        counter += 1
        assert "gi|11464971:4-101" == qresult.id
        assert "pleckstrin [Mus musculus]" == qresult.description
        assert 98 == qresult.seq_len
        assert 2933984 == qresult.stat_db_num
        assert 4726730735 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert 0.041 == qresult.stat_kappa
        assert 0.267 == qresult.stat_lambda
        assert 0.14 == qresult.stat_entropy
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "gi|354480463|ref|XM_003502378.1|" == hit.id
        assert "PREDICTED: Cricetulus griseus pleckstrin-like (LOC100773128), mRNA" == hit.description
        assert 1119 == hit.seq_len
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 205.297 == hsp.bitscore
        assert 521 == hsp.bitscore_raw
        assert 1.46172e-63 == hsp.evalue
        assert 0 == hsp.query_start
        assert 98 == hsp.query_end
        assert 75 == hsp.hit_start
        assert 369 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 1 == hsp.hit_frame
        assert 98 == hsp.ident_num
        assert 98 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 98 == hsp.aln_span
        assert "KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK" == hsp.hit.seq
        assert "KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK" == hsp.query.seq
        assert "KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[-1]
        assert 43.8986 == hsp.bitscore
        assert 102 == hsp.bitscore_raw
        assert 0.00054161 == hsp.evalue
        assert 2 == hsp.query_start
        assert 96 == hsp.query_end
        assert 801 == hsp.hit_start
        assert 1101 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 1 == hsp.hit_frame
        assert 30 == hsp.ident_num
        assert 50 == hsp.pos_num
        assert 6 == hsp.gap_num
        assert 100 == hsp.aln_span
        assert "IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNSPKGMIPLKGSTLTSPCQDF-GKRM---FVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA" == hsp.query.seq
        assert "IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGEDPLGAIHLRGCVVTSVESNHDGKKSDDENLFEIITADEVHYYLQAAAPKERTEWIKAIQVA" == hsp.hit.seq
        assert "I++G L+K+G     WK    +L ED   + +Y       P G I L+G  +TS   +  GK+     + +I T  +  ++ QAA  +ER  W++ I+ A" == hsp.aln_annotation["similarity"]

        hit = qresult[-1]
        assert "gi|390474391|ref|XM_002757683.2|" == hit.id
        assert "PREDICTED: Callithrix jacchus pleckstrin (PLEK), mRNA" == hit.description
        assert 1402 == hit.seq_len
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 202.986 == hsp.bitscore
        assert 515 == hsp.bitscore_raw
        assert 1.27031e-61 == hsp.evalue
        assert 0 == hsp.query_start
        assert 98 == hsp.query_end
        assert 160 == hsp.hit_start
        assert 454 == hsp.hit_end
        assert 0 == hsp.query_frame
        assert 2 == hsp.hit_frame
        assert 96 == hsp.ident_num
        assert 97 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 98 == hsp.aln_span
        assert "KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK" == hsp.query.seq
        assert "KRIREGYLVKKGSMFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK" == hsp.hit.seq
        assert "KRIREGYLVKKGS+FNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERDAWVRDIKKAIK" == hsp.aln_annotation["similarity"]

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 3 == counter


class TblastxXmlCases(unittest.TestCase):
    def test_xml_2212L_tblastx_001(self):
        xml_file = get_file("xml_2212L_tblastx_001.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        # test the first qresult
        qresult = next(qresults)
        counter += 1

        # test meta variables, only for the first one
        assert "2.2.12" == qresult.version
        assert REFERENCE == qresult.reference
        assert "BLOSUM80" == qresult.param_matrix
        assert 1 == qresult.param_evalue_threshold
        assert "L;" == qresult.param_filter
        assert 10 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "tblastx" == qresult.program
        assert "nr" == qresult.target

        # test parsed values of the first qresult
        assert "gi|1348853|gb|G26621.1|G26621" == qresult.id
        assert "human STS STS_D12006, sequence tagged site" == qresult.description
        assert 615 == qresult.seq_len
        assert 3533718 == qresult.stat_db_num
        # why is the value negative? is this a blast bug?
        assert -1496331888 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert 0.177051 == qresult.stat_kappa
        assert 0.342969 == qresult.stat_lambda
        assert 0.656794 == qresult.stat_entropy
        assert 10 == len(qresult)

        hit = qresult[0]
        assert "gi|18072170|gb|AC010333.7|" == hit.id
        assert "Homo sapiens chromosome 16 clone CTD-3037G24, complete sequence" == hit.description
        assert 159870 == hit.seq_len
        assert 13 == len(hit)

        hsp = hit.hsps[0]
        assert 329.561 == hsp.bitscore
        assert 661.0 == hsp.bitscore_raw
        assert 5.29552e-90 == hsp.evalue
        assert 1 == hsp.query_start
        assert 355 == hsp.query_end
        assert 44323 == hsp.hit_start
        assert 44677 == hsp.hit_end
        assert -3 == hsp.query_frame
        assert -3 == hsp.hit_frame
        assert 117 == hsp.ident_num
        assert 117 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 118 == hsp.aln_span
        assert "ECXFIMLYIFPARIWST*VICPPEQWL*RRKLSSGQKLLRRCGKTGYIKNNAGLK*PMRFCQGILH*CI*SKSGPIQRMWLKAPFWPFLFLLRALHTFPLFLSKWTK*RVS*VEGHSD" == hsp.query.seq
        assert "ECCFIMLYIFPARIWST*VICPPEQWL*RRKLSSGQKLLRRCGKTGYIKNNAGLK*PMRFCQGILH*CI*SKSGPIQRMWLKAPFWPFLFLLRALHTFPLFLSKWTK*RVS*VEGHSD" == hsp.hit.seq
        assert "EC FIMLYIFPARIWST*VICPPEQWL*RRKLSSGQKLLRRCGKTGYIKNNAGLK*PMRFCQGILH*CI*SKSGPIQRMWLKAPFWPFLFLLRALHTFPLFLSKWTK*RVS*VEGHSD" == hsp.aln_annotation["similarity"]

        hit = qresult[-1]
        assert "gi|4309961|gb|AC005993.2|AC005993" == hit.id
        assert "Homo sapiens PAC clone RP6-114E22 from 14, complete sequence" == hit.description
        assert 143943 == hit.seq_len
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 41.0922 == hsp.bitscore
        assert 78 == hsp.bitscore_raw
        assert 0.716571 == hsp.evalue
        assert 166 == hsp.query_start
        assert 250 == hsp.query_end
        assert 43679 == hsp.hit_start
        assert 43763 == hsp.hit_end
        assert 2 == hsp.query_frame
        assert 3 == hsp.hit_frame
        assert 13 == hsp.ident_num
        assert 19 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 28 == hsp.aln_span
        assert "PLTKAHRLFQTSIVFYVTCFTASSQQLL" == hsp.query.seq
        assert "PLNKYHTIFQISLCFYLFCYNMAQKQLL" == hsp.hit.seq
        assert "PL K H +FQ S+ FY+ C+  + +QLL" == hsp.aln_annotation["similarity"]

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_xml_2226_tblastx_001(self):
        xml_file = get_file("xml_2226_tblastx_001.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = next(qresults)
        counter += 1

        # test meta variables, only for the first one
        assert "2.2.26+" == qresult.version
        assert ("Stephen F. Altschul, Thomas L. Madden, Alejandro "
            "A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, "
            "Webb Miller, and David J. Lipman (1997), "
            '"Gapped BLAST and PSI-BLAST: a new generation of '
            'protein database search programs", '
            "Nucleic Acids Res. 25:3389-3402." == qresult.reference)
        assert "BLOSUM62" == qresult.param_matrix
        assert 10.0 == qresult.param_evalue_threshold
        assert "L;" == qresult.param_filter
        assert 11 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "tblastx" == qresult.program
        assert "db/minirefseq_mrna" == qresult.target

        # test parsed values of the first qresult
        assert "random_s00" == qresult.id
        assert "" == qresult.description
        assert 128 == qresult.seq_len
        assert 23 == qresult.stat_db_num
        assert 67750 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert -1 == qresult.stat_kappa
        assert -1 == qresult.stat_lambda
        assert -1 == qresult.stat_entropy
        assert 0 == len(qresult)

        # test parsed values of the second qresult
        qresult = next(qresults)
        counter += 1
        assert "gi|296147483:1-350" == qresult.id
        assert "Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA, complete cds" == qresult.description
        assert 350 == qresult.seq_len
        assert 23 == qresult.stat_db_num
        assert 67750 == qresult.stat_db_len
        assert 1954618.0 == qresult.stat_eff_space
        assert 0.133956144488482 == qresult.stat_kappa
        assert 0.317605957635731 == qresult.stat_lambda
        assert 0.401214524497119 == qresult.stat_entropy
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "gi|296147483|ref|NM_001183135.1|" == hit.id
        assert "Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA, complete cds" == hit.description
        assert 4911 == hit.seq_len
        assert 8 == len(hit)

        hsp = hit.hsps[0]
        assert 289.739 == hsp.bitscore
        assert 626 == hsp.bitscore_raw
        assert 2.35531e-81 == hsp.evalue
        assert 1 == hsp.query_start
        assert 349 == hsp.query_end
        assert 1 == hsp.hit_start
        assert 349 == hsp.hit_end
        assert 2 == hsp.query_frame
        assert 2 == hsp.hit_frame
        assert 116 == hsp.ident_num
        assert 116 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 116 == hsp.aln_span
        assert "WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K" == hsp.hit.seq
        assert "WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K" == hsp.query.seq
        assert "WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[-1]
        assert 18.9375 == hsp.bitscore
        assert 35 == hsp.bitscore_raw
        assert 7.78658 == hsp.evalue
        assert 292 == hsp.query_start
        assert 325 == hsp.query_end
        assert 340 == hsp.hit_start
        assert 373 == hsp.hit_end
        assert 2 == hsp.query_frame
        assert -3 == hsp.hit_frame
        assert 6 == hsp.ident_num
        assert 9 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 11 == hsp.aln_span
        assert "KFWMPSLRLLI" == hsp.query.seq
        assert "KKWVPPVKLLI" == hsp.hit.seq
        assert "K W+P ++LLI" == hsp.aln_annotation["similarity"]

        hit = qresult[-1]
        assert "gi|254579534|ref|XM_002495708.1|" == hit.id
        assert "Zygosaccharomyces rouxii hypothetical protein (ZYRO0C02266g) mRNA, complete cds" == hit.description
        assert 4866 == hit.seq_len
        assert 6 == len(hit)

        hsp = hit.hsps[0]
        assert 141.279 == hsp.bitscore
        assert 302 == hsp.bitscore_raw
        assert 1.15566e-36 == hsp.evalue
        assert 96 == hsp.query_start
        assert 348 == hsp.query_end
        assert 96 == hsp.hit_start
        assert 348 == hsp.hit_end
        assert 1 == hsp.query_frame
        assert 1 == hsp.hit_frame
        assert 57 == hsp.ident_num
        assert 72 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 84 == hsp.aln_span
        assert "IRHASDKSIEILKRVHSFEELERHPDFALPFVLACQSRNAKMTTLAMQCLQGLSTVPSIPRSRLSEILDAFIEATHLAMEIQLK" == hsp.query.seq
        assert "IRNASDKSIEILKVVHSYEELSRHPDFIVPLVMSCASKNAKLTTISMQCFQKLATVPCIPVDKLSDVLDAFIEANQLAMDIKLK" == hsp.hit.seq
        assert "IR+ASDKSIEILK VHS+EEL RHPDF +P V++C S+NAK+TT++MQC Q L+TVP IP  +LS++LDAFIEA  LAM+I+LK" == hsp.aln_annotation["similarity"]

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 2 == counter

    def test_xml_2226_tblastx_002(self):
        xml_file = get_file("xml_2226_tblastx_002.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1

        assert "2.2.26+" == qresult.version
        assert ("Stephen F. Altschul, Thomas L. Madden, Alejandro "
            "A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, "
            "Webb Miller, and David J. Lipman (1997), "
            '"Gapped BLAST and PSI-BLAST: a new generation of '
            'protein database search programs", '
            "Nucleic Acids Res. 25:3389-3402." == qresult.reference)
        assert "BLOSUM62" == qresult.param_matrix
        assert 10.0 == qresult.param_evalue_threshold
        assert "L;" == qresult.param_filter
        assert 11 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "tblastx" == qresult.program
        assert "db/minirefseq_mrna" == qresult.target

        assert "random_s00" == qresult.id
        assert "" == qresult.description
        assert 128 == qresult.seq_len
        assert 23 == qresult.stat_db_num
        assert 67750 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert -1 == qresult.stat_kappa
        assert -1 == qresult.stat_lambda
        assert -1 == qresult.stat_entropy
        assert 0 == len(qresult)

        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_xml_2226_tblastx_003(self):
        xml_file = get_file("xml_2226_tblastx_003.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1

        assert "2.2.26+" == qresult.version
        assert ("Stephen F. Altschul, Thomas L. Madden, Alejandro "
            "A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, "
            "Webb Miller, and David J. Lipman (1997), "
            '"Gapped BLAST and PSI-BLAST: a new generation of '
            'protein database search programs", '
            "Nucleic Acids Res. 25:3389-3402." == qresult.reference)
        assert "BLOSUM62" == qresult.param_matrix
        assert 10.0 == qresult.param_evalue_threshold
        assert "L;" == qresult.param_filter
        assert 11 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "tblastx" == qresult.program
        assert "db/minirefseq_mrna" == qresult.target

        assert "gi|296147483:1-350" == qresult.id
        assert "Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA, complete cds" == qresult.description
        assert 350 == qresult.seq_len
        assert 23 == qresult.stat_db_num
        assert 67750 == qresult.stat_db_len
        assert 1954618.0 == qresult.stat_eff_space
        assert 0.133956144488482 == qresult.stat_kappa
        assert 0.317605957635731 == qresult.stat_lambda
        assert 0.401214524497119 == qresult.stat_entropy
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "gi|296147483|ref|NM_001183135.1|" == hit.id
        assert "Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA, complete cds" == hit.description
        assert 4911 == hit.seq_len
        assert 8 == len(hit)

        hsp = hit.hsps[0]
        assert 289.739 == hsp.bitscore
        assert 626 == hsp.bitscore_raw
        assert 2.35531e-81 == hsp.evalue
        assert 1 == hsp.query_start
        assert 349 == hsp.query_end
        assert 1 == hsp.hit_start
        assert 349 == hsp.hit_end
        assert 2 == hsp.query_frame
        assert 2 == hsp.hit_frame
        assert 116 == hsp.ident_num
        assert 116 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 116 == hsp.aln_span
        assert "WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K" == hsp.hit.seq
        assert "WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K" == hsp.query.seq
        assert "WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[-1]
        assert 18.9375 == hsp.bitscore
        assert 35 == hsp.bitscore_raw
        assert 7.78658 == hsp.evalue
        assert 292 == hsp.query_start
        assert 325 == hsp.query_end
        assert 340 == hsp.hit_start
        assert 373 == hsp.hit_end
        assert 2 == hsp.query_frame
        assert -3 == hsp.hit_frame
        assert 6 == hsp.ident_num
        assert 9 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 11 == hsp.aln_span
        assert "KFWMPSLRLLI" == hsp.query.seq
        assert "KKWVPPVKLLI" == hsp.hit.seq
        assert "K W+P ++LLI" == hsp.aln_annotation["similarity"]

        hit = qresult[-1]
        assert "gi|254579534|ref|XM_002495708.1|" == hit.id
        assert "Zygosaccharomyces rouxii hypothetical protein (ZYRO0C02266g) mRNA, complete cds" == hit.description
        assert 4866 == hit.seq_len
        assert 6 == len(hit)

        hsp = hit.hsps[0]
        assert 141.279 == hsp.bitscore
        assert 302 == hsp.bitscore_raw
        assert 1.15566e-36 == hsp.evalue
        assert 96 == hsp.query_start
        assert 348 == hsp.query_end
        assert 96 == hsp.hit_start
        assert 348 == hsp.hit_end
        assert 1 == hsp.query_frame
        assert 1 == hsp.hit_frame
        assert 57 == hsp.ident_num
        assert 72 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 84 == hsp.aln_span
        assert "IRHASDKSIEILKRVHSFEELERHPDFALPFVLACQSRNAKMTTLAMQCLQGLSTVPSIPRSRLSEILDAFIEATHLAMEIQLK" == hsp.query.seq
        assert "IRNASDKSIEILKVVHSYEELSRHPDFIVPLVMSCASKNAKLTTISMQCFQKLATVPCIPVDKLSDVLDAFIEANQLAMDIKLK" == hsp.hit.seq
        assert "IR+ASDKSIEILK VHS+EEL RHPDF +P V++C S+NAK+TT++MQC Q L+TVP IP  +LS++LDAFIEA  LAM+I+LK" == hsp.aln_annotation["similarity"]

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_xml_2226_tblastx_004(self):
        xml_file = get_file("xml_2226_tblastx_004.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = next(qresults)
        counter += 1

        # test meta variables, only for the first one
        assert "2.2.26+" == qresult.version
        assert ("Stephen F. Altschul, Thomas L. Madden, Alejandro "
            "A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, "
            "Webb Miller, and David J. Lipman (1997), "
            '"Gapped BLAST and PSI-BLAST: a new generation of '
            'protein database search programs", '
            "Nucleic Acids Res. 25:3389-3402." == qresult.reference)
        assert "BLOSUM62" == qresult.param_matrix
        assert 10.0 == qresult.param_evalue_threshold
        assert "L;" == qresult.param_filter
        assert 11 == qresult.param_gap_open
        assert 1 == qresult.param_gap_extend
        assert "tblastx" == qresult.program
        assert "refseq_rna" == qresult.target

        # test parsed values of the first qresult
        assert "random_s00" == qresult.id
        assert "" == qresult.description
        assert 128 == qresult.seq_len
        assert 2933984 == qresult.stat_db_num
        assert 4726730735 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert 0 == qresult.stat_kappa
        assert 0 == qresult.stat_lambda
        assert 0 == qresult.stat_entropy
        assert 0 == len(qresult)

        # test parsed values of the second qresult
        qresult = next(qresults)
        counter += 1
        assert "gi|296147483:1-350" == qresult.id
        assert "Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA, complete cds" == qresult.description
        assert 350 == qresult.seq_len
        assert 2933984 == qresult.stat_db_num
        assert 4726730735 == qresult.stat_db_len
        assert 0 == qresult.stat_eff_space
        assert 0 == qresult.stat_kappa
        assert 0 == qresult.stat_lambda
        assert 0 == qresult.stat_entropy
        assert 5 == len(qresult)
        # check for alternative ID results
        assert qresult["gi|296147483|ref|NM_001183135.1|"] == qresult["gi|116616412|gb|EF059095.1|"]

        hit = qresult[0]
        assert "gi|296147483|ref|NM_001183135.1|" == hit.id
        assert "Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA, complete cds" == hit.description
        assert "gi|116616412|gb|EF059095.1|" == hit.id_all[1]
        assert ("Synthetic construct Saccharomyces cerevisiae "
            "clone FLH203015.01X MON2, complete sequence" == hit.description_all[1])
        assert 4911 == hit.seq_len
        assert 7 == len(hit)

        hsp = hit.hsps[0]
        assert 289.739 == hsp.bitscore
        assert 626 == hsp.bitscore_raw
        assert 1.05874e-76 == hsp.evalue
        assert 1 == hsp.query_start
        assert 349 == hsp.query_end
        assert 1 == hsp.hit_start
        assert 349 == hsp.hit_end
        assert 2 == hsp.query_frame
        assert 2 == hsp.hit_frame
        assert 116 == hsp.ident_num
        assert 116 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 116 == hsp.aln_span
        assert "WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K" == hsp.hit.seq
        assert "WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K" == hsp.query.seq
        assert "WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[-1]
        assert 36.3494 == hsp.bitscore
        assert 73 == hsp.bitscore_raw
        assert 9.12288e-54 == hsp.evalue
        assert 0 == hsp.query_start
        assert 42 == hsp.query_end
        assert 0 == hsp.hit_start
        assert 42 == hsp.hit_end
        assert 1 == hsp.query_frame
        assert 1 == hsp.hit_frame
        assert 14 == hsp.ident_num
        assert 14 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 14 == hsp.aln_span
        assert "MAMNTGGFDSMQRQ" == hsp.query.seq
        assert "MAMNTGGFDSMQRQ" == hsp.hit.seq
        assert "MAMNTGGFDSMQRQ" == hsp.aln_annotation["similarity"]

        hit = qresult[-1]
        assert "gi|254579534|ref|XM_002495708.1|" == hit.id
        assert ("Zygosaccharomyces rouxii hypothetical protein "
            "(ZYRO0C02266g) mRNA, complete cds" == hit.description)
        assert 4866 == hit.seq_len
        assert 4 == len(hit)

        hsp = hit.hsps[0]
        assert 141.279 == hsp.bitscore
        assert 302 == hsp.bitscore_raw
        assert 5.19486e-32 == hsp.evalue
        assert 96 == hsp.query_start
        assert 348 == hsp.query_end
        assert 96 == hsp.hit_start
        assert 348 == hsp.hit_end
        assert 1 == hsp.query_frame
        assert 1 == hsp.hit_frame
        assert 57 == hsp.ident_num
        assert 72 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 84 == hsp.aln_span
        assert "IRHASDKSIEILKRVHSFEELERHPDFALPFVLACQSRNAKMTTLAMQCLQGLSTVPSIPRSRLSEILDAFIEATHLAMEIQLK" == hsp.query.seq
        assert "IRNASDKSIEILKVVHSYEELSRHPDFIVPLVMSCASKNAKLTTISMQCFQKLATVPCIPVDKLSDVLDAFIEANQLAMDIKLK" == hsp.hit.seq
        assert "IR+ASDKSIEILK VHS+EEL RHPDF +P V++C S+NAK+TT++MQC Q L+TVP IP  +LS++LDAFIEA  LAM+I+LK" == hsp.aln_annotation["similarity"]

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 2 == counter


class BlastXmlSpecialCases(unittest.TestCase):
    def test_xml_2226_blastn_006(self):
        xml_file = get_file("xml_2226_blastn_006.xml")
        qresults = parse(xml_file, FMT)

        exp_warning = 1
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", BiopythonParserWarning)
            qresult = next(qresults)
            assert exp_warning == len(w), f"Expected {exp_warning} warning(s), got {len(w)}"

        assert qresult.blast_id == "Query_1"
        hit1 = qresult[0]
        hit2 = qresult[1]
        assert "gnl|BL_ORD_ID|18" == hit1.blast_id
        assert "gi|347972582|ref|XM_309352.4|" == hit1.id
        assert "Anopheles gambiae str. PEST AGAP011294-PA (DEFI_ANOGA) mRNA, complete cds" == hit1.description
        assert "gnl|BL_ORD_ID|17" == hit2.blast_id
        assert "gnl|BL_ORD_ID|17" == hit2.id
        assert "gi|347972582|ref|XM_309352.4| Anopheles gambiae str. PEST AGAP011294-PA (DEFI_ANOGA) mRNA, complete cds" == hit2.description

    def test_xml_2226_blastn_006_use_raw_hit_ids(self):
        xml_file = get_file("xml_2226_blastn_006.xml")
        qresults = parse(xml_file, FMT, use_raw_hit_ids=True)

        exp_warning = 0
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", BiopythonParserWarning)
            qresult = next(qresults)
            assert exp_warning == len(w), f"Expected {exp_warning} warning(s), got {len(w)}"

        assert qresult.blast_id == "Query_1"
        hit1 = qresult[0]
        hit2 = qresult[1]
        assert "gnl|BL_ORD_ID|18" == hit1.blast_id
        assert "gnl|BL_ORD_ID|18" == hit1.id
        assert "gi|347972582|ref|XM_309352.4| Anopheles gambiae str. PEST AGAP011294-PA (DEFI_ANOGA) mRNA, complete cds" == hit1.description
        assert "gnl|BL_ORD_ID|17" == hit2.blast_id
        assert "gnl|BL_ORD_ID|17" == hit2.id
        assert "gi|347972582|ref|XM_309352.4| Anopheles gambiae str. PEST AGAP011294-PA (DEFI_ANOGA) mRNA, complete cds" == hit2.description

    def test_xml_2226_blastn_006_use_raw_query_ids(self):
        xml_file = get_file("xml_2226_blastn_006.xml")
        qresults = parse(xml_file, FMT, use_raw_query_ids=True)

        exp_warning = 1
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", BiopythonParserWarning)
            qresult = next(qresults)
            assert exp_warning == len(w), f"Expected {exp_warning} warning(s), got {len(w)}"

        assert qresult.id == "Query_1"
        assert qresult.description == "gi|347972582|ref|XM_309352.4| Anopheles gambiae str. PEST AGAP011294-PA (DEFI_ANOGA) mRNA, complete cds"
        assert qresult.blast_id == "Query_1"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
