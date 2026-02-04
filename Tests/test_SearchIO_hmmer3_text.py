# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO HmmerIO parsers."""

import os
import unittest
import pytest

from Bio.SearchIO import parse

# test case files are in the Blast directory
TEST_DIR = "Hmmer"
FMT = "hmmer3-text"


def get_file(filename):
    """Return the path of a test file."""
    return os.path.join(TEST_DIR, filename)


class HmmscanCases(unittest.TestCase):
    """Testing hmmscan output."""

    def test_31b1_hmmscan_001(self):
        """Parsing hmmscan 3.1b1 (text_31b1_hmmscan_001)."""
        txt_file = get_file("text_31b1_hmmscan_001.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        # test first qresult
        qresult = next(qresults)
        counter += 1

        assert "hmmscan" == qresult.program
        assert "/home/bow/db/hmmer/protdb/Pfam-A.hmm" == qresult.target
        assert "3.1b1" == qresult.version
        assert "random_s00" == qresult.id
        assert 22 == qresult.seq_len
        assert 0 == len(qresult)

        # test fifth result
        for qresult in qresults:
            counter += 1

        assert "hmmscan" == qresult.program
        assert "/home/bow/db/hmmer/protdb/Pfam-A.hmm" == qresult.target
        assert "3.1b1" == qresult.version
        assert "gi|125490392|ref|NP_038661.2|" == qresult.id
        assert "POU domain, class 5, transcription factor 1 isoform 1 [Mus musculus]" == qresult.description
        assert 352 == qresult.seq_len
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "Pou" == hit.id
        assert "Pou domain - N-terminal to homeobox domain" == hit.description
        assert hit.is_included
        assert 7.6e-37 == hit.evalue
        assert 124.8 == hit.bitscore
        assert 0.5 == hit.bias
        assert 1.5 == hit.domain_exp_num
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 123.9 == hsp.bitscore
        assert 0.5 == hsp.bias
        assert 5e-40 == hsp.evalue_cond
        assert 1.5e-36 == hsp.evalue
        assert 2 == hsp.hit_start
        assert 75 == hsp.hit_end
        assert ".]" == hsp.hit_endtype
        assert 132 == hsp.query_start
        assert 205 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 130 == hsp.env_start
        assert 205 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.97 == hsp.acc_avg
        assert "eldleeleefakefkqrrikLgltqadvgsalgalyGkefsqttIcrFEalqLslknmckLkpllekWLeeae" == hsp.hit.seq
        assert "KALQKELEQFAKLLKQKRITLGYTQADVGLTLGVLFGKVFSQTTICRFEALQLSLKNMCKLRPLLEKWVEEAD" == hsp.query.seq
        assert "67899******************************************************************96" == hsp.aln_annotation["PP"]
        assert "++ ++ele+fak +kq+ri+Lg+tqadvg +lg+l+Gk+fsqttIcrFEalqLslknmckL+pllekW+eea+" == hsp.aln_annotation["similarity"]

        # last hit
        hit = qresult[4]
        assert "DUF521" == hit.id
        assert "Protein of unknown function (DUF521)" == hit.description
        assert not hit.is_included
        assert 0.15 == hit.evalue
        assert 10.5 == hit.bitscore
        assert 0.1 == hit.bias
        assert 1.4 == hit.domain_exp_num
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert not hsp.is_included
        assert 9.6 == hsp.bitscore
        assert 0.1 == hsp.bias
        assert 9.4e-05 == hsp.evalue_cond
        assert 0.28 == hsp.evalue
        assert 272 == hsp.hit_start
        assert 334 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 220 == hsp.query_start
        assert 280 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 196 == hsp.env_start
        assert 294 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.77 == hsp.acc_avg
        assert "adlaavleelnkakkeevdlvvlGcPhlsleeleelaellkgrkkkvsvelvvttsravlsk" == hsp.hit.seq
        assert "QARKRKRTSIENRVRWSLETMFLKCPKPSLQQITHIANQLGLEK--DVVRVWFCNRRQKGKR" == hsp.query.seq
        assert "345666667778888899************************99..9999999988876554" == hsp.aln_annotation["PP"]
        assert "+ +++ + +++++   +++ ++l cP  sl++++++a++l  +k    v+++ +  r+  ++" == hsp.aln_annotation["similarity"]

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 5 == counter

    def test_30_hmmscan_001(self):
        """Parsing hmmscan 3.0 (text_30_hmmscan_001)."""
        txt_file = get_file("text_30_hmmscan_001.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        # test first qresult
        qresult = next(qresults)
        counter += 1

        assert "hmmscan" == qresult.program
        assert "/home/bow/db/hmmer/Pfam-A.hmm" == qresult.target
        assert "3.0" == qresult.version
        assert "random_s00" == qresult.id
        assert 32 == qresult.seq_len
        assert 0 == len(qresult)

        # test second result
        qresult = next(qresults)
        counter += 1

        assert "hmmscan" == qresult.program
        assert "/home/bow/db/hmmer/Pfam-A.hmm" == qresult.target
        assert "3.0" == qresult.version
        assert "gi|4885477|ref|NP_005359.1|" == qresult.id
        assert "myoglobin [Homo sapiens]" == qresult.description
        assert 154 == qresult.seq_len
        assert 1 == len(qresult)

        hit = qresult[0]
        assert "Globin" == hit.id
        assert "Globin" == hit.description
        assert hit.is_included
        assert 6e-21 == hit.evalue
        assert 74.6 == hit.bitscore
        assert 0.3 == hit.bias
        assert 1.3 == hit.domain_exp_num
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 74.0 == hsp.bitscore
        assert 0.2 == hsp.bias
        assert 6.7e-25 == hsp.evalue_cond
        assert 9.2e-21 == hsp.evalue
        assert 0 == hsp.hit_start
        assert 107 == hsp.hit_end
        assert "[." == hsp.hit_endtype
        assert 6 == hsp.query_start
        assert 112 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 6 == hsp.env_start
        assert 113 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.97 == hsp.acc_avg
        assert "HHHHHHHHHHHHCHHHHHHHHHHHHHHHHHSGGGGGGGCCCTTTT.HHHHHTSCHHHHHHHHHHHHHHHHHHCTTSHHHHHHHHHHHHHHHHTT-.--HHHHCCHHHHH" == hsp.aln_annotation["CS"]
        assert "qkalvkaswekvkanaeeigaeilkrlfkaypdtkklFkkfgdls.aedlksspkfkahakkvlaaldeavknldnddnlkaalkklgarHakrg.vdpanfklfgeal" == hsp.hit.seq
        assert "EWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKsEDEMKASEDLKKHGATVLTALGGILKK---KGHHEAEIKPLAQSHATKHkIPVKYLEFISECI" == hsp.query.seq
        assert "5789*********************************************************************...6899***********************999998" == hsp.aln_annotation["PP"]
        assert "+++lv   w+kv+a+++ +g+e+l rlfk +p+t ++F kf+ l+  +++k s+++k+h+++vl al+ ++k+   ++ ++a++k l+++Ha+++ ++ ++ + ++e++" == hsp.aln_annotation["similarity"]

        # test third result
        qresult = next(qresults)
        counter += 1

        assert "hmmscan" == qresult.program
        assert "/home/bow/db/hmmer/Pfam-A.hmm" == qresult.target
        assert "3.0" == qresult.version
        assert "gi|126362951:116-221" == qresult.id
        assert "leukocyte immunoglobulin-like receptor subfamily B member 1 isoform 2 precursor [Homo sapiens]" == qresult.description
        assert 106 == qresult.seq_len
        assert 2 == len(qresult)

        hit = qresult[0]
        assert "Ig_3" == hit.id
        assert "Immunoglobulin domain" == hit.description
        assert hit.is_included
        assert 1.4e-09 == hit.evalue
        assert 38.2 == hit.bitscore
        assert 0.4 == hit.bias
        assert 1.3 == hit.domain_exp_num
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 37.6 == hsp.bitscore
        assert 0.3 == hsp.bias
        assert 3e-13 == hsp.evalue_cond
        assert 2.1e-09 == hsp.evalue
        assert 0 == hsp.hit_start
        assert 73 == hsp.hit_end
        assert "[." == hsp.hit_endtype
        assert 8 == hsp.query_start
        assert 84 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 8 == hsp.env_start
        assert 88 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.94 == hsp.acc_avg
        assert "kPvisvspsptvtsggnvtLtCsaeggpppptisWy.....ietppelqgsegssssestLtissvtsedsgtYtCva" == hsp.hit.seq
        assert "KPTLSAQPSPVVNSGGNVILQCDSQVA--FDGFSLCkegedEHPQCLNSQPHARGSSRAIFSVGPVSPSRRWWYRCYA" == hsp.query.seq
        assert "8************************99..78888888****************************************9" == hsp.aln_annotation["PP"]
        assert "kP++s++psp+v+sggnv L+C ++     + +s +     +++   +++++ ++ss +++++ +v+++ +  Y+C+a" == hsp.aln_annotation["similarity"]

        hit = qresult[-1]
        assert "Ig_2" == hit.id
        assert "Immunoglobulin domain" == hit.description
        assert 3.5e-05 == hit.evalue
        assert 23.7 == hit.bitscore
        assert 0.1 == hit.bias
        assert 1.1 == hit.domain_exp_num
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 23.4 == hsp.bitscore
        assert 0.1 == hsp.bias
        assert 6.2e-09 == hsp.evalue_cond
        assert 4.3e-05 == hsp.evalue
        assert 0 == hsp.hit_start
        assert 80 == hsp.hit_end
        assert "[]" == hsp.hit_endtype
        assert 8 == hsp.query_start
        assert 104 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 8 == hsp.env_start
        assert 104 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.71 == hsp.acc_avg
        assert "kpvlvapp.svvtegenvtLtCsapgnptprvqwykdg.vels......qsqnq........lfipnvsaedsgtYtCra....rnseggktstsveltv" == hsp.hit.seq
        assert "KPTLSAQPsPVVNSGGNVILQCDSQVA-FDGFSLCKEGeDEHPqclnsqP---HargssraiFSVGPVSPSRRWWYRCYAydsnSPYEWSLPSDLLELLV" == hsp.query.seq
        assert "799998885779*************85.899***9988655554443320...134455543444669************88443344588888888766" == hsp.aln_annotation["PP"]
        assert "kp+l+a+p +vv++g nv L+C ++    + +++ k+g +e +      +            + +  vs++    Y+C+a    + +e++ +s+ +el v" == hsp.aln_annotation["similarity"]

        # test fourth result
        qresult = next(qresults)
        counter += 1

        assert "hmmscan" == qresult.program
        assert "/home/bow/db/hmmer/Pfam-A.hmm" == qresult.target
        assert "3.0" == qresult.version
        assert "gi|22748937|ref|NP_065801.1|" == qresult.id
        assert "exportin-5 [Homo sapiens]" == qresult.description
        assert 1204 == qresult.seq_len
        assert 2 == len(qresult)

        hit = qresult[0]
        assert "Xpo1" == hit.id
        assert "Exportin 1-like protein" == hit.description
        assert hit.is_included
        assert 7.8e-34 == hit.evalue
        assert 116.6 == hit.bitscore
        assert 7.8 == hit.bias
        assert 2.8 == hit.domain_exp_num
        assert 2 == hit.domain_obs_num
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 116.1 == hsp.bitscore
        assert 3.4 == hsp.bias
        assert 1.6e-37 == hsp.evalue_cond
        assert 1.1e-33 == hsp.evalue
        assert 1 == hsp.hit_start
        assert 148 == hsp.hit_end
        assert ".]" == hsp.hit_endtype
        assert 109 == hsp.query_start
        assert 271 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 108 == hsp.env_start
        assert 271 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.98 == hsp.acc_avg
        assert "HHHHHHHHHHHHHHHHHHTTTTSTTHHHHHHHHHHG-HHHHHHHHHHHHHHHHHHCCS-TTTS-CCCHHHHHHHCHHHHHHHHHHHHHHHC-TT-..................HHHHHHHHHHHHHHCTTS-CHHCHCS...HHHHHCHHCCSCCCHHHHHHHH" == hsp.aln_annotation["CS"]
        assert "kflrnklaealaelflqeypnqWpsffddllsllssspsglelllriLkvlpeEiadfsrskleqerrnelkdllrsqvqkilelllqileqsvskk...............sselveatLkclsswvswidiglivnsp..llsllfqlLndpelreaAvecL" == hsp.hit.seq
        assert "NHIKDALSRIVVEMIKREWPQHWPDMLIELDTLSKQGETQTELVMFILLRLAEDVVTF--QTLPPQRRRDIQQTLTQNMERIFSFLLNTLQENVNKYqqvktdtsqeskaqaNCRVGVAALNTLAGYIDWVSMSHITAENckLLEILCLLLNEQELQLGAAECL" == hsp.query.seq
        assert "89******************************************************99..79*********************************99*****************************************8889*********************8" == hsp.aln_annotation["PP"]
        assert "+++++ l+++++e++++e+p++Wp+++ +l  l++++++++el++ iL++l+e++++f  ++l  +rr++++++l++++++i+++ll+ l+++v+k+               ++++  a+L++l+ +++w+++++i +++  ll++l+ lLn++el+  A+ecL" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[-1]
        assert 2 == hsp.domain_index
        assert not hsp.is_included
        assert -1.8 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 0.35 == hsp.evalue_cond
        assert 2.4e03 == hsp.evalue
        assert 111 == hsp.hit_start
        assert 139 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 498 == hsp.query_start
        assert 525 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 495 == hsp.env_start
        assert 529 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.86 == hsp.acc_avg
        assert "HHCTTS-CHHCHCS.HHHHHCHHCCSCC" == hsp.aln_annotation["CS"]
        assert "swvswidiglivnspllsllfqlLndpe" == hsp.hit.seq
        assert "SFVQWEAMTLFLES-VITQMFRTLNREE" == hsp.query.seq
        assert "899*********98.8888899998776" == hsp.aln_annotation["PP"]
        assert "s+v+w  ++l+++s +++ +f+ Ln++e" == hsp.aln_annotation["similarity"]

        hit = qresult[-1]
        assert "IBN_N" == hit.id
        assert "Importin-beta N-terminal domain" == hit.description
        assert hit.is_included
        assert 0.0039 == hit.evalue
        assert 16.9 == hit.bitscore
        assert 0.0 == hit.bias
        assert 2.7 == hit.domain_exp_num
        assert 2 == hit.domain_obs_num
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 14.0 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 4.8e-06 == hsp.evalue_cond
        assert 0.033 == hsp.evalue
        assert 3 == hsp.hit_start
        assert 75 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 35 == hsp.query_start
        assert 98 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 32 == hsp.env_start
        assert 100 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.87 == hsp.acc_avg
        assert "HHHHHHHSCTHHHHHHHHHHHTTTSTHHHHHHHHHHHHHHHHHSCCHHHHHHHHCS-HHHHHHHHHHHHHHH" == hsp.aln_annotation["CS"]
        assert "qLnqlekqkPgflsallqilanksldlevRqlAalyLknlItkhWkseeaqrqqqlpeeekelIrnnllnll" == hsp.hit.seq
        assert "FCEEFKEKCPICVPCGLRLA-EKTQVAIVRHFGLQILEHVVKFRWN--------GMSRLEKVYLKNSVMELI" == hsp.query.seq
        assert "56788886699*********.6555899******************........999999****99999887" == hsp.aln_annotation["PP"]
        assert " ++++++  P ++++ l+++ +k+    vR++++++L++ ++ +W+         ++  ek +++n++ +l+" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[-1]
        assert 2 == hsp.domain_index
        assert not hsp.is_included
        assert -3.3 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 1.2 == hsp.evalue_cond
        assert 8e03 == hsp.evalue
        assert 56 == hsp.hit_start
        assert 75 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 167 == hsp.query_start
        assert 186 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 164 == hsp.env_start
        assert 187 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.85 == hsp.acc_avg
        assert "HCS-HHHHHHHHHHHHHHH" == hsp.aln_annotation["CS"]
        assert "qqlpeeekelIrnnllnll" == hsp.hit.seq
        assert "QTLPPQRRRDIQQTLTQNM" == hsp.query.seq
        assert "6899*******99998865" == hsp.aln_annotation["PP"]
        assert "q+lp++ + +I ++l + +" == hsp.aln_annotation["similarity"]

        # test fifth result
        qresult = next(qresults)
        counter += 1

        assert "hmmscan" == qresult.program
        assert "/home/bow/db/hmmer/Pfam-A.hmm" == qresult.target
        assert "3.0" == qresult.version
        assert "gi|125490392|ref|NP_038661.2|" == qresult.id
        assert "POU domain, class 5, transcription factor 1 isoform 1 [Mus musculus]" == qresult.description
        assert 352 == qresult.seq_len
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "Pou" == hit.id
        assert "Pou domain - N-terminal to homeobox domain" == hit.description
        assert hit.is_included
        assert 7e-37 == hit.evalue
        assert 124.8 == hit.bitscore
        assert 0.5 == hit.bias
        assert 1.5 == hit.domain_exp_num
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 123.9 == hsp.bitscore
        assert 0.3 == hsp.bias
        assert 5e-40 == hsp.evalue_cond
        assert 1.4e-36 == hsp.evalue
        assert 2 == hsp.hit_start
        assert 75 == hsp.hit_end
        assert ".]" == hsp.hit_endtype
        assert 132 == hsp.query_start
        assert 205 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 130 == hsp.env_start
        assert 205 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.97 == hsp.acc_avg
        assert "eldleeleefakefkqrrikLgltqadvgsalgalyGkefsqttIcrFEalqLslknmckLkpllekWLeeae" == hsp.hit.seq
        assert "KALQKELEQFAKLLKQKRITLGYTQADVGLTLGVLFGKVFSQTTICRFEALQLSLKNMCKLRPLLEKWVEEAD" == hsp.query.seq
        assert "67899******************************************************************96" == hsp.aln_annotation["PP"]
        assert "++ ++ele+fak +kq+ri+Lg+tqadvg +lg+l+Gk+fsqttIcrFEalqLslknmckL+pllekW+eea+" == hsp.aln_annotation["similarity"]

        hit = qresult[1]
        assert "Homeobox" == hit.id
        assert "Homeobox domain" == hit.description
        assert hit.is_included
        assert 2.1e-18 == hit.evalue
        assert 65.5 == hit.bitscore
        assert 1.1 == hit.bias
        assert 1.5 == hit.domain_exp_num
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 64.6 == hsp.bitscore
        assert 0.7 == hsp.bias
        assert 1.5e-21 == hsp.evalue_cond
        assert 4.1e-18 == hsp.evalue
        assert 0 == hsp.hit_start
        assert 57 == hsp.hit_end
        assert "[]" == hsp.hit_endtype
        assert 223 == hsp.query_start
        assert 280 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 223 == hsp.env_start
        assert 280 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.98 == hsp.acc_avg
        assert "SS--SS--HHHHHHHHHHCCTSSS--HHHHHHHHHH----HHHHHHHHHHHHHHHHH" == hsp.aln_annotation["CS"]
        assert "rrkRttftkeqleeLeelFeknrypsaeereeLAkklgLterqVkvWFqNrRakekk" == hsp.hit.seq
        assert "KRKRTSIENRVRWSLETMFLKCPKPSLQQITHIANQLGLEKDVVRVWFCNRRQKGKR" == hsp.query.seq
        assert "79****************************************************997" == hsp.aln_annotation["PP"]
        assert "+rkRt++++     Le +F k+++ps ++++++A++lgL++++V+vWF+NrR+k k+" == hsp.aln_annotation["similarity"]

        hit = qresult[2]
        assert "HTH_31" == hit.id
        assert "Helix-turn-helix domain" == hit.description
        assert not hit.is_included
        assert 0.012 == hit.evalue
        assert 15.6 == hit.bitscore
        assert 0.0 == hit.bias
        assert 2.2 == hit.domain_exp_num
        assert 2 == hit.domain_obs_num
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert not hsp.is_included
        assert 12.0 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 5.7e-05 == hsp.evalue_cond
        assert 0.16 == hsp.evalue
        assert 0 == hsp.hit_start
        assert 35 == hsp.hit_end
        assert "[." == hsp.hit_endtype
        assert 140 == hsp.query_start
        assert 181 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 140 == hsp.env_start
        assert 184 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.96 == hsp.acc_avg
        assert "aLGarLralReraGLtqeevAerlg......vSastlsrlE" == hsp.hit.seq
        assert "QFAKLLKQKRITLGYTQADVGLTLGvlfgkvFSQTTICRFE" == hsp.query.seq
        assert "6999***********************************99" == hsp.aln_annotation["PP"]
        assert "+++ +L++ R + G tq++v+  lg      +S++t++r E" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[1]
        assert 2 == hsp.domain_index
        assert not hsp.is_included
        assert 0.8 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 0.19 == hsp.evalue_cond
        assert 5.2e02 == hsp.evalue
        assert 38 == hsp.hit_start
        assert 62 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 244 == hsp.query_start
        assert 268 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 242 == hsp.env_start
        assert 270 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.86 == hsp.acc_avg
        assert "rgrpsaavlaalaralgldpaera" == hsp.hit.seq
        assert "CPKPSLQQITHIANQLGLEKDVVR" == hsp.query.seq
        assert "678**************9988765" == hsp.aln_annotation["PP"]
        assert "++ ps+++++ +a+ lgl+ + ++" == hsp.aln_annotation["similarity"]

        hit = qresult[3]
        assert "Homeobox_KN" == hit.id
        assert "Homeobox KN domain" == hit.description
        assert not hit.is_included
        assert 0.039 == hit.evalue
        assert 13.5 == hit.bitscore
        assert 0.0 == hit.bias
        assert 1.6 == hit.domain_exp_num
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert not hsp.is_included
        assert 12.3 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 3.5e-05 == hsp.evalue_cond
        assert 0.095 == hsp.evalue
        assert 6 == hsp.hit_start
        assert 39 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 243 == hsp.query_start
        assert 276 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 240 == hsp.env_start
        assert 277 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.91 == hsp.acc_avg
        assert "hnPYPskevkeelakqTglsrkqidnWFiNaRr" == hsp.hit.seq
        assert "KCPKPSLQQITHIANQLGLEKDVVRVWFCNRRQ" == hsp.query.seq
        assert "56779*************************996" == hsp.aln_annotation["PP"]
        assert "+ P Ps +++  +a+q gl  + +  WF N R " == hsp.aln_annotation["similarity"]

        hit = qresult[4]
        assert "DUF521" == hit.id
        assert "Protein of unknown function (DUF521)" == hit.description
        assert not hit.is_included
        assert 0.14 == hit.evalue
        assert 10.5 == hit.bitscore
        assert 0.1 == hit.bias
        assert 1.4 == hit.domain_exp_num
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert not hsp.is_included
        assert 9.6 == hsp.bitscore
        assert 0.1 == hsp.bias
        assert 9.4e-05 == hsp.evalue_cond
        assert 0.26 == hsp.evalue
        assert 272 == hsp.hit_start
        assert 334 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 220 == hsp.query_start
        assert 280 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 196 == hsp.env_start
        assert 294 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.77 == hsp.acc_avg
        assert "adlaavleelnkakkeevdlvvlGcPhlsleeleelaellkgrkkkvsvelvvttsravlsk" == hsp.hit.seq
        assert "QARKRKRTSIENRVRWSLETMFLKCPKPSLQQITHIANQLGLEK--DVVRVWFCNRRQKGKR" == hsp.query.seq
        assert "345666667778888899************************99..9999999988876554" == hsp.aln_annotation["PP"]
        assert "+ +++ + +++++   +++ ++l cP  sl++++++a++l  +k    v+++ +  r+  ++" == hsp.aln_annotation["similarity"]

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 5 == counter

    def test_30_hmmscan_002(self):
        """Parsing hmmscan 3.0 (text_30_hmmscan_002)."""
        txt_file = get_file("text_30_hmmscan_002.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        # test first qresult
        qresult = next(qresults)
        counter += 1

        assert "hmmscan" == qresult.program
        assert "/home/bow/db/hmmer/Pfam-A.hmm" == qresult.target
        assert "3.0" == qresult.version
        assert "random_s00" == qresult.id
        assert 32 == qresult.seq_len
        assert 0 == len(qresult)

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_30_hmmscan_003(self):
        """Parsing hmmscan 3.0 (text_30_hmmscan_003)."""
        txt_file = get_file("text_30_hmmscan_003.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        # test first qresult
        qresult = next(qresults)
        counter += 1

        assert "hmmscan" == qresult.program
        assert "/home/bow/db/hmmer/Pfam-A.hmm" == qresult.target
        assert "3.0" == qresult.version
        assert "gi|4885477|ref|NP_005359.1|" == qresult.id
        assert "myoglobin [Homo sapiens]" == qresult.description
        assert 154 == qresult.seq_len
        assert 1 == len(qresult)

        hit = qresult[0]
        assert "Globin" == hit.id
        assert "Globin" == hit.description
        assert hit.is_included
        assert 6e-21 == hit.evalue
        assert 74.6 == hit.bitscore
        assert 0.3 == hit.bias
        assert 1.3 == hit.domain_exp_num
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 74.0 == hsp.bitscore
        assert 0.2 == hsp.bias
        assert 6.7e-25 == hsp.evalue_cond
        assert 9.2e-21 == hsp.evalue
        assert 0 == hsp.hit_start
        assert 107 == hsp.hit_end
        assert "[." == hsp.hit_endtype
        assert 6 == hsp.query_start
        assert 112 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 6 == hsp.env_start
        assert 113 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.97 == hsp.acc_avg
        assert "HHHHHHHHHHHHCHHHHHHHHHHHHHHHHHSGGGGGGGCCCTTTT.HHHHHTSCHHHHHHHHHHHHHHHHHHCTTSHHHHHHHHHHHHHHHHTT-.--HHHHCCHHHHH" == hsp.aln_annotation["CS"]
        assert "qkalvkaswekvkanaeeigaeilkrlfkaypdtkklFkkfgdls.aedlksspkfkahakkvlaaldeavknldnddnlkaalkklgarHakrg.vdpanfklfgeal" == hsp.hit.seq
        assert "EWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKsEDEMKASEDLKKHGATVLTALGGILKK---KGHHEAEIKPLAQSHATKHkIPVKYLEFISECI" == hsp.query.seq
        assert "5789*********************************************************************...6899***********************999998" == hsp.aln_annotation["PP"]
        assert "+++lv   w+kv+a+++ +g+e+l rlfk +p+t ++F kf+ l+  +++k s+++k+h+++vl al+ ++k+   ++ ++a++k l+++Ha+++ ++ ++ + ++e++" == hsp.aln_annotation["similarity"]

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_30_hmmscan_004(self):
        """Parsing hmmscan 3.0 (text_30_hmmscan_004)."""
        txt_file = get_file("text_30_hmmscan_004.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        # test first qresult
        qresult = next(qresults)
        counter += 1

        assert "hmmscan" == qresult.program
        assert "/home/bow/db/hmmer/Pfam-A.hmm" == qresult.target
        assert "3.0" == qresult.version
        assert "gi|126362951:116-221" == qresult.id
        assert "leukocyte immunoglobulin-like receptor subfamily B member 1 isoform 2 precursor [Homo sapiens]" == qresult.description
        assert 106 == qresult.seq_len
        assert 2 == len(qresult)

        hit = qresult[0]
        assert "Ig_3" == hit.id
        assert "Immunoglobulin domain" == hit.description
        assert hit.is_included
        assert 1.4e-09 == hit.evalue
        assert 38.2 == hit.bitscore
        assert 0.4 == hit.bias
        assert 1.3 == hit.domain_exp_num
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 37.6 == hsp.bitscore
        assert 0.3 == hsp.bias
        assert 3e-13 == hsp.evalue_cond
        assert 2.1e-09 == hsp.evalue
        assert 0 == hsp.hit_start
        assert 73 == hsp.hit_end
        assert "[." == hsp.hit_endtype
        assert 8 == hsp.query_start
        assert 84 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 8 == hsp.env_start
        assert 88 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.94 == hsp.acc_avg
        assert "kPvisvspsptvtsggnvtLtCsaeggpppptisWy.....ietppelqgsegssssestLtissvtsedsgtYtCva" == hsp.hit.seq
        assert "KPTLSAQPSPVVNSGGNVILQCDSQVA--FDGFSLCkegedEHPQCLNSQPHARGSSRAIFSVGPVSPSRRWWYRCYA" == hsp.query.seq
        assert "8************************99..78888888****************************************9" == hsp.aln_annotation["PP"]
        assert "kP++s++psp+v+sggnv L+C ++     + +s +     +++   +++++ ++ss +++++ +v+++ +  Y+C+a" == hsp.aln_annotation["similarity"]

        hit = qresult[-1]
        assert "Ig_2" == hit.id
        assert "Immunoglobulin domain" == hit.description
        assert 3.5e-05 == hit.evalue
        assert 23.7 == hit.bitscore
        assert 0.1 == hit.bias
        assert 1.1 == hit.domain_exp_num
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 23.4 == hsp.bitscore
        assert 0.1 == hsp.bias
        assert 6.2e-09 == hsp.evalue_cond
        assert 4.3e-05 == hsp.evalue
        assert 0 == hsp.hit_start
        assert 80 == hsp.hit_end
        assert "[]" == hsp.hit_endtype
        assert 8 == hsp.query_start
        assert 104 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 8 == hsp.env_start
        assert 104 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.71 == hsp.acc_avg
        assert "kpvlvapp.svvtegenvtLtCsapgnptprvqwykdg.vels......qsqnq........lfipnvsaedsgtYtCra....rnseggktstsveltv" == hsp.hit.seq
        assert "KPTLSAQPsPVVNSGGNVILQCDSQVA-FDGFSLCKEGeDEHPqclnsqP---HargssraiFSVGPVSPSRRWWYRCYAydsnSPYEWSLPSDLLELLV" == hsp.query.seq
        assert "799998885779*************85.899***9988655554443320...134455543444669************88443344588888888766" == hsp.aln_annotation["PP"]
        assert "kp+l+a+p +vv++g nv L+C ++    + +++ k+g +e +      +            + +  vs++    Y+C+a    + +e++ +s+ +el v" == hsp.aln_annotation["similarity"]

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_30_hmmscan_005(self):
        """Parsing hmmscan 3.0 (text_30_hmmscan_005)."""
        txt_file = get_file("text_30_hmmscan_005.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        # test first result
        qresult = next(qresults)
        counter += 1

        assert "hmmscan" == qresult.program
        assert "/home/bow/db/hmmer/Pfam-A.hmm" == qresult.target
        assert "3.0" == qresult.version
        assert "gi|22748937|ref|NP_065801.1|" == qresult.id
        assert "exportin-5 [Homo sapiens]" == qresult.description
        assert 1204 == qresult.seq_len
        assert 2 == len(qresult)

        hit = qresult[0]
        assert "Xpo1" == hit.id
        assert "Exportin 1-like protein" == hit.description
        assert hit.is_included
        assert 7.8e-34 == hit.evalue
        assert 116.6 == hit.bitscore
        assert 7.8 == hit.bias
        assert 2.8 == hit.domain_exp_num
        assert 2 == hit.domain_obs_num
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 116.1 == hsp.bitscore
        assert 3.4 == hsp.bias
        assert 1.6e-37 == hsp.evalue_cond
        assert 1.1e-33 == hsp.evalue
        assert 1 == hsp.hit_start
        assert 148 == hsp.hit_end
        assert ".]" == hsp.hit_endtype
        assert 109 == hsp.query_start
        assert 271 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 108 == hsp.env_start
        assert 271 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.98 == hsp.acc_avg
        assert "HHHHHHHHHHHHHHHHHHTTTTSTTHHHHHHHHHHG-HHHHHHHHHHHHHHHHHHCCS-TTTS-CCCHHHHHHHCHHHHHHHHHHHHHHHC-TT-..................HHHHHHHHHHHHHHCTTS-CHHCHCS...HHHHHCHHCCSCCCHHHHHHHH" == hsp.aln_annotation["CS"]
        assert "kflrnklaealaelflqeypnqWpsffddllsllssspsglelllriLkvlpeEiadfsrskleqerrnelkdllrsqvqkilelllqileqsvskk...............sselveatLkclsswvswidiglivnsp..llsllfqlLndpelreaAvecL" == hsp.hit.seq
        assert "NHIKDALSRIVVEMIKREWPQHWPDMLIELDTLSKQGETQTELVMFILLRLAEDVVTF--QTLPPQRRRDIQQTLTQNMERIFSFLLNTLQENVNKYqqvktdtsqeskaqaNCRVGVAALNTLAGYIDWVSMSHITAENckLLEILCLLLNEQELQLGAAECL" == hsp.query.seq
        assert "89******************************************************99..79*********************************99*****************************************8889*********************8" == hsp.aln_annotation["PP"]
        assert "+++++ l+++++e++++e+p++Wp+++ +l  l++++++++el++ iL++l+e++++f  ++l  +rr++++++l++++++i+++ll+ l+++v+k+               ++++  a+L++l+ +++w+++++i +++  ll++l+ lLn++el+  A+ecL" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[-1]
        assert 2 == hsp.domain_index
        assert not hsp.is_included
        assert -1.8 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 0.35 == hsp.evalue_cond
        assert 2.4e03 == hsp.evalue
        assert 111 == hsp.hit_start
        assert 139 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 498 == hsp.query_start
        assert 525 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 495 == hsp.env_start
        assert 529 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.86 == hsp.acc_avg
        assert "HHCTTS-CHHCHCS.HHHHHCHHCCSCC" == hsp.aln_annotation["CS"]
        assert "swvswidiglivnspllsllfqlLndpe" == hsp.hit.seq
        assert "SFVQWEAMTLFLES-VITQMFRTLNREE" == hsp.query.seq
        assert "899*********98.8888899998776" == hsp.aln_annotation["PP"]
        assert "s+v+w  ++l+++s +++ +f+ Ln++e" == hsp.aln_annotation["similarity"]

        hit = qresult[-1]
        assert "IBN_N" == hit.id
        assert "Importin-beta N-terminal domain" == hit.description
        assert hit.is_included
        assert 0.0039 == hit.evalue
        assert 16.9 == hit.bitscore
        assert 0.0 == hit.bias
        assert 2.7 == hit.domain_exp_num
        assert 2 == hit.domain_obs_num
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 14.0 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 4.8e-06 == hsp.evalue_cond
        assert 0.033 == hsp.evalue
        assert 3 == hsp.hit_start
        assert 75 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 35 == hsp.query_start
        assert 98 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 32 == hsp.env_start
        assert 100 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.87 == hsp.acc_avg
        assert "HHHHHHHSCTHHHHHHHHHHHTTTSTHHHHHHHHHHHHHHHHHSCCHHHHHHHHCS-HHHHHHHHHHHHHHH" == hsp.aln_annotation["CS"]
        assert "qLnqlekqkPgflsallqilanksldlevRqlAalyLknlItkhWkseeaqrqqqlpeeekelIrnnllnll" == hsp.hit.seq
        assert "FCEEFKEKCPICVPCGLRLA-EKTQVAIVRHFGLQILEHVVKFRWN--------GMSRLEKVYLKNSVMELI" == hsp.query.seq
        assert "56788886699*********.6555899******************........999999****99999887" == hsp.aln_annotation["PP"]
        assert " ++++++  P ++++ l+++ +k+    vR++++++L++ ++ +W+         ++  ek +++n++ +l+" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[-1]
        assert 2 == hsp.domain_index
        assert not hsp.is_included
        assert -3.3 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 1.2 == hsp.evalue_cond
        assert 8e03 == hsp.evalue
        assert 56 == hsp.hit_start
        assert 75 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 167 == hsp.query_start
        assert 186 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 164 == hsp.env_start
        assert 187 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.85 == hsp.acc_avg
        assert "HCS-HHHHHHHHHHHHHHH" == hsp.aln_annotation["CS"]
        assert "qqlpeeekelIrnnllnll" == hsp.hit.seq
        assert "QTLPPQRRRDIQQTLTQNM" == hsp.query.seq
        assert "6899*******99998865" == hsp.aln_annotation["PP"]
        assert "q+lp++ + +I ++l + +" == hsp.aln_annotation["similarity"]

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_30_hmmscan_006(self):
        """Parsing hmmscan 3.0 (text_30_hmmscan_006)."""
        txt_file = get_file("text_30_hmmscan_006.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        # test first qresult
        qresult = next(qresults)
        counter += 1

        assert "hmmscan" == qresult.program
        assert "/home/bow/db/hmmer/Pfam-A.hmm" == qresult.target
        assert "3.0" == qresult.version
        assert "gi|125490392|ref|NP_038661.2|" == qresult.id
        assert "POU domain, class 5, transcription factor 1 isoform 1 [Mus musculus]" == qresult.description
        assert 352 == qresult.seq_len
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "Pou" == hit.id
        assert "Pou domain - N-terminal to homeobox domain" == hit.description
        assert hit.is_included
        assert 7e-37 == hit.evalue
        assert 124.8 == hit.bitscore
        assert 0.5 == hit.bias
        assert 1.5 == hit.domain_exp_num
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 123.9 == hsp.bitscore
        assert 0.3 == hsp.bias
        assert 5e-40 == hsp.evalue_cond
        assert 1.4e-36 == hsp.evalue
        assert 2 == hsp.hit_start
        assert 75 == hsp.hit_end
        assert ".]" == hsp.hit_endtype
        assert 132 == hsp.query_start
        assert 205 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 130 == hsp.env_start
        assert 205 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.97 == hsp.acc_avg
        assert "eldleeleefakefkqrrikLgltqadvgsalgalyGkefsqttIcrFEalqLslknmckLkpllekWLeeae" == hsp.hit.seq
        assert "KALQKELEQFAKLLKQKRITLGYTQADVGLTLGVLFGKVFSQTTICRFEALQLSLKNMCKLRPLLEKWVEEAD" == hsp.query.seq
        assert "67899******************************************************************96" == hsp.aln_annotation["PP"]
        assert "++ ++ele+fak +kq+ri+Lg+tqadvg +lg+l+Gk+fsqttIcrFEalqLslknmckL+pllekW+eea+" == hsp.aln_annotation["similarity"]

        hit = qresult[1]
        assert "Homeobox" == hit.id
        assert "Homeobox domain" == hit.description
        assert hit.is_included
        assert 2.1e-18 == hit.evalue
        assert 65.5 == hit.bitscore
        assert 1.1 == hit.bias
        assert 1.5 == hit.domain_exp_num
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 64.6 == hsp.bitscore
        assert 0.7 == hsp.bias
        assert 1.5e-21 == hsp.evalue_cond
        assert 4.1e-18 == hsp.evalue
        assert 0 == hsp.hit_start
        assert 57 == hsp.hit_end
        assert "[]" == hsp.hit_endtype
        assert 223 == hsp.query_start
        assert 280 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 223 == hsp.env_start
        assert 280 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.98 == hsp.acc_avg
        assert "SS--SS--HHHHHHHHHHCCTSSS--HHHHHHHHHH----HHHHHHHHHHHHHHHHH" == hsp.aln_annotation["CS"]
        assert "rrkRttftkeqleeLeelFeknrypsaeereeLAkklgLterqVkvWFqNrRakekk" == hsp.hit.seq
        assert "KRKRTSIENRVRWSLETMFLKCPKPSLQQITHIANQLGLEKDVVRVWFCNRRQKGKR" == hsp.query.seq
        assert "79****************************************************997" == hsp.aln_annotation["PP"]
        assert "+rkRt++++     Le +F k+++ps ++++++A++lgL++++V+vWF+NrR+k k+" == hsp.aln_annotation["similarity"]

        hit = qresult[2]
        assert "HTH_31" == hit.id
        assert "Helix-turn-helix domain" == hit.description
        assert not hit.is_included
        assert 0.012 == hit.evalue
        assert 15.6 == hit.bitscore
        assert 0.0 == hit.bias
        assert 2.2 == hit.domain_exp_num
        assert 2 == hit.domain_obs_num
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert not hsp.is_included
        assert 12.0 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 5.7e-05 == hsp.evalue_cond
        assert 0.16 == hsp.evalue
        assert 0 == hsp.hit_start
        assert 35 == hsp.hit_end
        assert "[." == hsp.hit_endtype
        assert 140 == hsp.query_start
        assert 181 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 140 == hsp.env_start
        assert 184 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.96 == hsp.acc_avg
        assert "aLGarLralReraGLtqeevAerlg......vSastlsrlE" == hsp.hit.seq
        assert "QFAKLLKQKRITLGYTQADVGLTLGvlfgkvFSQTTICRFE" == hsp.query.seq
        assert "6999***********************************99" == hsp.aln_annotation["PP"]
        assert "+++ +L++ R + G tq++v+  lg      +S++t++r E" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[1]
        assert 2 == hsp.domain_index
        assert not hsp.is_included
        assert 0.8 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 0.19 == hsp.evalue_cond
        assert 5.2e02 == hsp.evalue
        assert 38 == hsp.hit_start
        assert 62 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 244 == hsp.query_start
        assert 268 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 242 == hsp.env_start
        assert 270 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.86 == hsp.acc_avg
        assert "rgrpsaavlaalaralgldpaera" == hsp.hit.seq
        assert "CPKPSLQQITHIANQLGLEKDVVR" == hsp.query.seq
        assert "678**************9988765" == hsp.aln_annotation["PP"]
        assert "++ ps+++++ +a+ lgl+ + ++" == hsp.aln_annotation["similarity"]

        hit = qresult[3]
        assert "Homeobox_KN" == hit.id
        assert "Homeobox KN domain" == hit.description
        assert not hit.is_included
        assert 0.039 == hit.evalue
        assert 13.5 == hit.bitscore
        assert 0.0 == hit.bias
        assert 1.6 == hit.domain_exp_num
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert not hsp.is_included
        assert 12.3 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 3.5e-05 == hsp.evalue_cond
        assert 0.095 == hsp.evalue
        assert 6 == hsp.hit_start
        assert 39 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 243 == hsp.query_start
        assert 276 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 240 == hsp.env_start
        assert 277 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.91 == hsp.acc_avg
        assert "hnPYPskevkeelakqTglsrkqidnWFiNaRr" == hsp.hit.seq
        assert "KCPKPSLQQITHIANQLGLEKDVVRVWFCNRRQ" == hsp.query.seq
        assert "56779*************************996" == hsp.aln_annotation["PP"]
        assert "+ P Ps +++  +a+q gl  + +  WF N R " == hsp.aln_annotation["similarity"]

        hit = qresult[4]
        assert "DUF521" == hit.id
        assert "Protein of unknown function (DUF521)" == hit.description
        assert not hit.is_included
        assert 0.14 == hit.evalue
        assert 10.5 == hit.bitscore
        assert 0.1 == hit.bias
        assert 1.4 == hit.domain_exp_num
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert not hsp.is_included
        assert 9.6 == hsp.bitscore
        assert 0.1 == hsp.bias
        assert 9.4e-05 == hsp.evalue_cond
        assert 0.26 == hsp.evalue
        assert 272 == hsp.hit_start
        assert 334 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 220 == hsp.query_start
        assert 280 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 196 == hsp.env_start
        assert 294 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.77 == hsp.acc_avg
        assert "adlaavleelnkakkeevdlvvlGcPhlsleeleelaellkgrkkkvsvelvvttsravlsk" == hsp.hit.seq
        assert "QARKRKRTSIENRVRWSLETMFLKCPKPSLQQITHIANQLGLEK--DVVRVWFCNRRQKGKR" == hsp.query.seq
        assert "345666667778888899************************99..9999999988876554" == hsp.aln_annotation["PP"]
        assert "+ +++ + +++++   +++ ++l cP  sl++++++a++l  +k    v+++ +  r+  ++" == hsp.aln_annotation["similarity"]

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_30_hmmscan_007(self):
        """Parsing hmmscan 3.0 (text_30_hmmscan_007)."""
        txt_file = get_file("text_30_hmmscan_007.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        # test first qresult
        qresult = next(qresults)
        counter += 1

        assert "hmmscan" == qresult.program
        assert "/home/bow/db/hmmer/Pfam-A.hmm" == qresult.target
        assert "3.0" == qresult.version
        assert "gi|125490392|ref|NP_038661.2|" == qresult.id
        assert "POU domain, class 5, transcription factor 1 isoform 1 [Mus musculus]" == qresult.description
        assert 352 == qresult.seq_len
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "Pou" == hit.id
        assert "Pou domain - N-terminal to homeobox domain" == hit.description
        assert hit.is_included
        assert 7e-37 == hit.evalue
        assert 124.8 == hit.bitscore
        assert 0.5 == hit.bias
        assert 1.5 == hit.domain_exp_num
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 123.9 == hsp.bitscore
        assert 0.3 == hsp.bias
        assert 5e-40 == hsp.evalue_cond
        assert 1.4e-36 == hsp.evalue
        assert 2 == hsp.hit_start
        assert 75 == hsp.hit_end
        assert ".]" == hsp.hit_endtype
        assert 132 == hsp.query_start
        assert 205 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 130 == hsp.env_start
        assert 205 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.97 == hsp.acc_avg

        hit = qresult[1]
        assert "Homeobox" == hit.id
        assert "Homeobox domain" == hit.description
        assert hit.is_included
        assert 2.1e-18 == hit.evalue
        assert 65.5 == hit.bitscore
        assert 1.1 == hit.bias
        assert 1.5 == hit.domain_exp_num
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 64.6 == hsp.bitscore
        assert 0.7 == hsp.bias
        assert 1.5e-21 == hsp.evalue_cond
        assert 4.1e-18 == hsp.evalue
        assert 0 == hsp.hit_start
        assert 57 == hsp.hit_end
        assert "[]" == hsp.hit_endtype
        assert 223 == hsp.query_start
        assert 280 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 223 == hsp.env_start
        assert 280 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.98 == hsp.acc_avg

        hit = qresult[2]
        assert "HTH_31" == hit.id
        assert "Helix-turn-helix domain" == hit.description
        assert not hit.is_included
        assert 0.012 == hit.evalue
        assert 15.6 == hit.bitscore
        assert 0.0 == hit.bias
        assert 2.2 == hit.domain_exp_num
        assert 2 == hit.domain_obs_num
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert not hsp.is_included
        assert 12.0 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 5.7e-05 == hsp.evalue_cond
        assert 0.16 == hsp.evalue
        assert 0 == hsp.hit_start
        assert 35 == hsp.hit_end
        assert "[." == hsp.hit_endtype
        assert 140 == hsp.query_start
        assert 181 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 140 == hsp.env_start
        assert 184 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.96 == hsp.acc_avg

        hsp = hit.hsps[1]
        assert 2 == hsp.domain_index
        assert not hsp.is_included
        assert 0.8 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 0.19 == hsp.evalue_cond
        assert 5.2e02 == hsp.evalue
        assert 38 == hsp.hit_start
        assert 62 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 244 == hsp.query_start
        assert 268 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 242 == hsp.env_start
        assert 270 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.86 == hsp.acc_avg

        hit = qresult[3]
        assert "Homeobox_KN" == hit.id
        assert "Homeobox KN domain" == hit.description
        assert not hit.is_included
        assert 0.039 == hit.evalue
        assert 13.5 == hit.bitscore
        assert 0.0 == hit.bias
        assert 1.6 == hit.domain_exp_num
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert not hsp.is_included
        assert 12.3 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 3.5e-05 == hsp.evalue_cond
        assert 0.095 == hsp.evalue
        assert 6 == hsp.hit_start
        assert 39 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 243 == hsp.query_start
        assert 276 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 240 == hsp.env_start
        assert 277 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.91 == hsp.acc_avg

        hit = qresult[4]
        assert "DUF521" == hit.id
        assert "Protein of unknown function (DUF521)" == hit.description
        assert not hit.is_included
        assert 0.14 == hit.evalue
        assert 10.5 == hit.bitscore
        assert 0.1 == hit.bias
        assert 1.4 == hit.domain_exp_num
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert not hsp.is_included
        assert 9.6 == hsp.bitscore
        assert 0.1 == hsp.bias
        assert 9.4e-05 == hsp.evalue_cond
        assert 0.26 == hsp.evalue
        assert 272 == hsp.hit_start
        assert 334 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 220 == hsp.query_start
        assert 280 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 196 == hsp.env_start
        assert 294 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.77 == hsp.acc_avg

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_30_hmmscan_008(self):
        """Parsing hmmscan 3.0 (text_30_hmmscan_008)."""
        txt_file = get_file("text_30_hmmscan_008.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        # test first result
        qresult = next(qresults)
        counter += 1

        assert "hmmscan" == qresult.program
        assert "/home/bow/db/hmmer/Pfam-A.hmm" == qresult.target
        assert "3.0" == qresult.version
        assert "gi|22748937|ref|NP_065801.1|" == qresult.id
        assert "exportin-5 [Homo sapiens]" == qresult.description
        assert 1204 == qresult.seq_len
        assert 2 == len(qresult)

        hit = qresult[0]
        assert "Xpo1" == hit.id
        assert "Exportin 1-like protein" == hit.description
        assert hit.is_included
        assert 7.8e-34 == hit.evalue
        assert 116.6 == hit.bitscore
        assert 7.8 == hit.bias
        assert 2.8 == hit.domain_exp_num
        assert 2 == hit.domain_obs_num
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 116.1 == hsp.bitscore
        assert 3.4 == hsp.bias
        assert 1.6e-37 == hsp.evalue_cond
        assert 1.1e-33 == hsp.evalue
        assert 1 == hsp.hit_start
        assert 148 == hsp.hit_end
        assert ".]" == hsp.hit_endtype
        assert 109 == hsp.query_start
        assert 271 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 108 == hsp.env_start
        assert 271 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.98 == hsp.acc_avg
        assert "HHHHHHHHHHHHHHHHHHTTTTSTTHHHHHHHHHHG-HHHHHHHHHHHHHHHHHHCCS-TTTS-CCCHHHHHHHCHHHHHHHHHHHHHHHC-TT-..................HHHHHHHHHHHHHHCTTS-CHHCHCS...HHHHHCHHCCSCCCHHHHHHHH" == hsp.aln_annotation["CS"]
        assert "kflrnklaealaelflqeypnqWpsffddllsllssspsglelllriLkvlpeEiadfsrskleqerrnelkdllrsqvqkilelllqileqsvskk...............sselveatLkclsswvswidiglivnsp..llsllfqlLndpelreaAvecL" == hsp.hit.seq
        assert "NHIKDALSRIVVEMIKREWPQHWPDMLIELDTLSKQGETQTELVMFILLRLAEDVVTF--QTLPPQRRRDIQQTLTQNMERIFSFLLNTLQENVNKYqqvktdtsqeskaqaNCRVGVAALNTLAGYIDWVSMSHITAENckLLEILCLLLNEQELQLGAAECL" == hsp.query.seq
        assert "89******************************************************99..79*********************************99*****************************************8889*********************8" == hsp.aln_annotation["PP"]
        assert "+++++ l+++++e++++e+p++Wp+++ +l  l++++++++el++ iL++l+e++++f  ++l  +rr++++++l++++++i+++ll+ l+++v+k+               ++++  a+L++l+ +++w+++++i +++  ll++l+ lLn++el+  A+ecL" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[-1]
        assert 2 == hsp.domain_index
        assert not hsp.is_included
        assert -1.8 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 0.35 == hsp.evalue_cond
        assert 2.4e03 == hsp.evalue
        assert 111 == hsp.hit_start
        assert 139 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 498 == hsp.query_start
        assert 525 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 495 == hsp.env_start
        assert 529 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.86 == hsp.acc_avg
        assert "HHCTTS-CHHCHCS.HHHHHCHHCCSCC" == hsp.aln_annotation["CS"]
        assert "swvswidiglivnspllsllfqlLndpe" == hsp.hit.seq
        assert "SFVQWEAMTLFLES-VITQMFRTLNREE" == hsp.query.seq
        assert "899*********98.8888899998776" == hsp.aln_annotation["PP"]
        assert "s+v+w  ++l+++s +++ +f+ Ln++e" == hsp.aln_annotation["similarity"]

        hit = qresult[-1]
        assert "IBN_N" == hit.id
        assert "Importin-beta N-terminal domain" == hit.description
        assert hit.is_included
        assert 0.0039 == hit.evalue
        assert 16.9 == hit.bitscore
        assert 0.0 == hit.bias
        assert 2.7 == hit.domain_exp_num
        assert 2 == hit.domain_obs_num
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 14.0 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 4.8e-06 == hsp.evalue_cond
        assert 0.033 == hsp.evalue
        assert 3 == hsp.hit_start
        assert 75 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 35 == hsp.query_start
        assert 98 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 32 == hsp.env_start
        assert 100 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.87 == hsp.acc_avg
        assert "HHHHHHHSCTHHHHHHHHHHHTTTSTHHHHHHHHHHHHHHHHHSCCHHHHHHHHCS-HHHHHHHHHHHHHHH" == hsp.aln_annotation["CS"]
        assert "qLnqlekqkPgflsallqilanksldlevRqlAalyLknlItkhWkseeaqrqqqlpeeekelIrnnllnll" == hsp.hit.seq
        assert "FCEEFKEKCPICVPCGLRLA-EKTQVAIVRHFGLQILEHVVKFRWN--------GMSRLEKVYLKNSVMELI" == hsp.query.seq
        assert "56788886699*********.6555899******************........999999****99999887" == hsp.aln_annotation["PP"]
        assert " ++++++  P ++++ l+++ +k+    vR++++++L++ ++ +W+         ++  ek +++n++ +l+" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[-1]
        assert 2 == hsp.domain_index
        assert not hsp.is_included
        assert -3.3 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 1.2 == hsp.evalue_cond
        assert 8e03 == hsp.evalue
        assert 56 == hsp.hit_start
        assert 75 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 167 == hsp.query_start
        assert 186 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 164 == hsp.env_start
        assert 187 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.85 == hsp.acc_avg
        assert "HCS-HHHHHHHHHHHHHHH" == hsp.aln_annotation["CS"]
        assert "qqlpeeekelIrnnllnll" == hsp.hit.seq
        assert "QTLPPQRRRDIQQTLTQNM" == hsp.query.seq
        assert "6899*******99998865" == hsp.aln_annotation["PP"]
        assert "q+lp++ + +I ++l + +" == hsp.aln_annotation["similarity"]

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_30_hmmscan_009(self):
        """Parsing hmmscan 3.0 (text_30_hmmscan_009)."""
        hmmer_file = get_file("text_30_hmmscan_009.out")
        qresults = parse(hmmer_file, FMT)
        counter = 0

        # test first qresult
        qresult = next(qresults)
        counter += 1
        assert "SCO3574" == qresult.id
        assert 5 == len(qresult.hits)
        assert "Esterase" == qresult.hits[3].id

    def test_30_hmmscan_010(self):
        """Parsing hmmscan 3.0 (text_30_hmmscan_010)."""
        hmmer_file = get_file("text_30_hmmscan_010.out")
        qresults = list(parse(hmmer_file, FMT))

        # test the Hit object without HSPs
        hit = qresults[0][-1]
        assert not hit
        assert "NRPS-COM_Cterm" == hit.id
        assert "" == hit.description
        assert "bpsA" == hit.query_id
        assert "<unknown description>" == hit.query_description
        assert 4.4e-11 == hit.evalue
        assert 33.6 == hit.bitscore
        assert 10.2 == hit.bias
        assert 2.9 == hit.domain_exp_num
        assert 0 == hit.domain_obs_num
        assert 0 == len(hit)

    def test_30_hmmscan_011(self):
        """Parsing hmmscan 3.0 (text_30_hmmscan_011)."""
        hmmer_file = get_file("text_30_hmmscan_011.out")
        qresults = list(parse(hmmer_file, FMT))
        assert len(qresults) == 2

        result = qresults[0]
        hit = result[-1]
        assert len(result) == 29
        assert hit.hsps[0].bitscore == 22.4

        result = qresults[1]
        hit = result[-1]
        assert len(result) == 29
        assert hit.hsps[0].bitscore == 20.8

        # Test getting program name when preamble contains full path and whitespace
        expected_db_target = (
            "C:\\msys64\\scr\\byerly\\builds\\NB\\20231025 _\\internal\\lib\\"
            "site-packages\\anarci\\dat\\HMMs\\ALL.hmm"
        )
        for idx, qresult in enumerate(qresults):
            assert qresult.program == "hmmscan", (
                "Expected program 'hmmscan' for item at index "
                f"{idx}, found {qresult.program}"
            )
            # hmm db path also contains space
            assert qresult.target == expected_db_target


class HmmersearchCases(unittest.TestCase):
    """Test hmmsearch output."""

    def test_31b2_hmmsearch_001(self):
        """Test parsing hmmsearch 3.1b2 (text_31b2_hmmsearch_001)."""
        txt_file = get_file("text_31b2_hmmsearch_001.out")
        qresults = list(parse(txt_file, FMT))

        # Assert we deal with only 1 query
        assert len(qresults) == 1

        # Test whether proper query id is read
        assert qresults[0].id == "infile_sto"

        # Test if proper number of hits is read
        assert len(qresults[0].hits) == 10

    def test_31b2_hmmsearch_002(self):
        """Test parsing hmmsearch 3.1b2 (text_31b2_hmmsearch_002)."""
        txt_file = get_file("text_31b2_hmmsearch_002.out")
        qresults = list(parse(txt_file, FMT))

        # Assert we deal with 2 queries
        assert len(qresults) == 2

        # Test whether proper query id is read
        assert [x.id for x in qresults] == ["infile_sto", "infile_sto2"]

        # Test if proper number of hits is read
        assert [len(x.hits) for x in qresults] == [10, 10]

    def test_31b2_hmmscan_001(self):
        """Test parsing hmmscan 3.1b2 (text_31b2_hmmscan_001)."""
        txt_file = get_file("text_31b2_hmmscan_001.out")
        qresults = list(parse(txt_file, FMT))

        # expects only 1 query result
        assert 1 == len(qresults)

        # in that query result, there are 108 hits
        qresult = qresults[0]
        assert 108 == len(qresult)

        # making sure all query IDs are equal
        assert ("Protein-arginine kinase activator protein "
            "OS=Bacillus subtilis (strain 168) GN=mcsA PE=1 SV=1" == qresult.description)

        # and all hits have no descriptions
        for hit in qresult:
            assert "" == hit.description

    def test_31b1_hmmsearch_001(self):
        """Test parsing hmmsearch 3.1b1 (text_31b1_hmmsearch_001)."""
        txt_file = get_file("text_31b1_hmmsearch_001.out")
        qresults = list(parse(txt_file, FMT))

        assert 2 == len(qresults)

        # first qresult is empty
        qresult = qresults[0]
        assert "Globins" == qresult.id
        assert 149 == qresult.seq_len
        assert 0 == len(qresult)

        # second qresult
        qresult = qresults[1]
        assert "Pkinase" == qresult.id
        assert 260 == qresult.seq_len
        assert "PF00069.17" == qresult.accession
        assert "Protein kinase domain" == qresult.description
        assert 4 == len(qresult)

        # first hit, first hsp
        hit = qresult[0]
        assert "sp|Q9WUT3|KS6A2_MOUSE" == hit.id
        assert "Ribosomal protein S6 kinase alpha-2 OS=Mus musculus GN=Rps6ka2 PE=1 SV=1" == hit.description
        assert hit.is_included
        assert 8.5e-147 == hit.evalue
        assert 492.3 == hit.bitscore
        assert 0.0 == hit.bias
        assert 2.1 == hit.domain_exp_num
        assert 2 == hit.domain_obs_num
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 241.2 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 2.6e-75 == hsp.evalue_cond
        assert 3.6e-70 == hsp.evalue
        assert 0 == hsp.query_start
        assert 260 == hsp.query_end
        assert "[]" == hsp.query_endtype
        assert 58 == hsp.hit_start
        assert 318 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 58 == hsp.env_start
        assert 318 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.95 == hsp.acc_avg
        assert "EEEEEEEEEETTEEEEEEEE...TTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHHSSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEEEE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTT.....HHHHHTSGGG" == hsp.aln_annotation["CS"]
        assert "yelleklGsGsfGkVykakk...kktgkkvAvKilkkeeekskkektavrElkilkklsHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgevkiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRlt.....aeeilkhpwl" == hsp.query.seq
        assert "FELLKVLGQGSYGKVFLVRKvtgSDAGQLYAMKVLKKATLKVRDRVRSKMERDILAEVNHPFIVKLHYAFQTEGKLYLILDFLRGGDLFTRLSKEVMFTEEDVKFYLAELALALDHLHGLGIIYRDLKPENILLDEEGHIKITDFGLSKEATDHDKRAYSFCGTIEYMAPEVVN-RRGHTQSADWWSFGVLMFEMLTGSLPFQGK---DRKETMALILKAKLGMPQFLS----AEAQSLLRALFKRNPCNRLGagvdgVEEIKRHPFF" == hsp.hit.seq
        assert "67899**********7666611155667*****************99999****************************************************************************************************************************.******************************...999999999999999998866....99******************9999999*****997" == hsp.aln_annotation["PP"]
        assert "+ell+ lG+Gs+GkV+ ++k   ++ g+ +A+K+lkk + k + + +++ E  il++++Hp+ivkl+ +f+t+ +lyl+l++++ggdlf+ l+ke  ++ee++k+++ +++ +l++lH  gii+rDLKpeNiLld++g++ki+DFGl+k+++ +++++++++gt eYmAPEv++ ++++t+++D+Ws+Gv+++e+ltg+lpf+g+   d++e+++ ilk kl  ++  s    +e+++l++ l++++p +Rl      +eei++hp++" == hsp.aln_annotation["similarity"]

    def test_30_hmmsearch_001(self):
        """Parsing hmmersearch 3.0 (text_30_hmmsearch_001)."""
        txt_file = get_file("text_30_hmmsearch_001.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1

        assert "hmmsearch" == qresult.program
        assert "/home/bow/db/hmmer/uniprot_sprot.fasta" == qresult.target
        assert "3.0" == qresult.version
        assert "globins4" == qresult.id
        assert 149 == qresult.seq_len
        assert 0 == len(qresult)

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_30_hmmsearch_002(self):
        """Parsing hmmersearch 3.0 (text_30_hmmsearch_002)."""
        txt_file = get_file("text_30_hmmsearch_002.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1

        assert "hmmsearch" == qresult.program
        assert "/home/bow/db/hmmer/uniprot_sprot.fasta" == qresult.target
        assert "3.0" == qresult.version
        assert "Pkinase" == qresult.id
        assert "PF00069.17" == qresult.accession
        assert "Protein kinase domain" == qresult.description
        assert 260 == qresult.seq_len
        assert 7 == len(qresult)

        hit = qresult[0]
        assert "sp|Q9WUT3|KS6A2_MOUSE" == hit.id
        assert "Ribosomal protein S6 kinase alpha-2 OS=Mus musculus GN=Rps6ka2 PE=2 SV=1" == hit.description
        assert hit.is_included
        assert 8.4e-147 == hit.evalue
        assert 492.3 == hit.bitscore
        assert 0.0 == hit.bias
        assert 2.1 == hit.domain_exp_num
        assert 2 == hit.domain_obs_num
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 241.2 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 4.6e-75 == hsp.evalue_cond
        assert 3.5e-70 == hsp.evalue
        assert 0 == hsp.query_start
        assert 260 == hsp.query_end
        assert "[]" == hsp.query_endtype
        assert 58 == hsp.hit_start
        assert 318 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 58 == hsp.env_start
        assert 318 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.95 == hsp.acc_avg
        assert "EEEEEEEEEETTEEEEEEEE...TTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHHSSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEEEE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTT.....HHHHHTSGGG" == hsp.aln_annotation["CS"]
        assert "yelleklGsGsfGkVykakk...kktgkkvAvKilkkeeekskkektavrElkilkklsHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgevkiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRlt.....aeeilkhpwl" == hsp.query.seq
        assert "FELLKVLGQGSYGKVFLVRKvtgSDAGQLYAMKVLKKATLKVRDRVRSKMERDILAEVNHPFIVKLHYAFQTEGKLYLILDFLRGGDLFTRLSKEVMFTEEDVKFYLAELALALDHLHGLGIIYRDLKPENILLDEEGHIKITDFGLSKEATDHDKRAYSFCGTIEYMAPEVVN-RRGHTQSADWWSFGVLMFEMLTGSLPFQGK---DRKETMALILKAKLGMPQFLS----AEAQSLLRALFKRNPCNRLGagvdgVEEIKRHPFF" == hsp.hit.seq
        assert "67899**********7666611155667*****************99999****************************************************************************************************************************.******************************...999999999999999998866....99******************9999999*****997" == hsp.aln_annotation["PP"]
        assert "+ell+ lG+Gs+GkV+ ++k   ++ g+ +A+K+lkk + k + + +++ E  il++++Hp+ivkl+ +f+t+ +lyl+l++++ggdlf+ l+ke  ++ee++k+++ +++ +l++lH  gii+rDLKpeNiLld++g++ki+DFGl+k+++ +++++++++gt eYmAPEv++ ++++t+++D+Ws+Gv+++e+ltg+lpf+g+   d++e+++ ilk kl  ++  s    +e+++l++ l++++p +Rl      +eei++hp++" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[1]
        assert 2 == hsp.domain_index
        assert hsp.is_included
        assert 249.3 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 1.5e-77 == hsp.evalue_cond
        assert 1.1e-72 == hsp.evalue
        assert 0 == hsp.query_start
        assert 260 == hsp.query_end
        assert "[]" == hsp.query_endtype
        assert 414 == hsp.hit_start
        assert 672 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 414 == hsp.env_start
        assert 672 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.97 == hsp.acc_avg
        assert "EEEEEEEEEETTEEEEEEEETTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHH.SSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEE....EE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTTHHHHHTSGGG" == hsp.aln_annotation["CS"]
        assert "yelleklGsGsfGkVykakkkktgkkvAvKilkkeeekskkektavrElkilkkl.sHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgev....kiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRltaeeilkhpwl" == hsp.query.seq
        assert "YEIKEDIGVGSYSVCKRCVHKATDAEYAVKIIDKSKRDPSE------EIEILLRYgQHPNIITLKDVYDDGKYVYLVMELMRGGELLDRILRQRCFSEREASDVLYTIARTMDYLHSQGVVHRDLKPSNILYMDESGNpesiRICDFGFAKQLRAENGLLMTPCYTANFVAPEVLK-RQGYDAACDVWSLGILLYTMLAGFTPFANGPDDTPEEILARIGSGKYALSGGNWDSISDAAKDVVSKMLHVDPQQRLTAVQVLKHPWI" == hsp.hit.seq
        assert "7899***********************************98......9*******99**************************************************************************98544444888**********************************.***************************************************************************************7" == hsp.aln_annotation["PP"]
        assert "ye++e +G Gs+++++++++k t  ++AvKi++k++++ ++      E++il +  +Hpni++l +v+ + +++ylv+e+++gg+l d + +++++se+e+ +++++i++ ++ylHs+g++HrDLKp+NiL  +++      +i+DFG+ak+l  +++ l t + t +++APEvl+ +++y++++DvWslG++ly++l g +pf +  +++ +e++++i ++k+  +  +++s s+ +kd+++k+l+ dp++Rlta ++lkhpw+" == hsp.aln_annotation["similarity"]

        hit = qresult[-1]
        assert "sp|P18654|KS6A3_MOUSE" == hit.id
        assert "Ribosomal protein S6 kinase alpha-3 OS=Mus musculus GN=Rps6ka3 PE=1 SV=2" == hit.description
        assert hit.is_included
        assert 5e-144 == hit.evalue
        assert 483.2 == hit.bitscore
        assert 0.0 == hit.bias
        assert 2.2 == hit.domain_exp_num
        assert 2 == hit.domain_obs_num
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 240.3 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 8.7e-75 == hsp.evalue_cond
        assert 6.6e-70 == hsp.evalue
        assert 0 == hsp.query_start
        assert 260 == hsp.query_end
        assert "[]" == hsp.query_endtype
        assert 67 == hsp.hit_start
        assert 327 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 67 == hsp.env_start
        assert 327 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.95 == hsp.acc_avg
        assert "EEEEEEEEEETTEEEEEEEETTTTE...EEEEEEEEHHHCCCCCCHHHHHHHHHHHHHSSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEEEE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTT.....HHHHHTSGGG" == hsp.aln_annotation["CS"]
        assert "yelleklGsGsfGkVykakkkktgk...kvAvKilkkeeekskkektavrElkilkklsHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgevkiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRlt.....aeeilkhpwl" == hsp.query.seq
        assert "FELLKVLGQGSFGKVFLVKKISGSDarqLYAMKVLKKATLKVRDRVRTKMERDILVEVNHPFIVKLHYAFQTEGKLYLILDFLRGGDLFTRLSKEVMFTEEDVKFYLAELALALDHLHSLGIIYRDLKPENILLDEEGHIKLTDFGLSKESIDHEKKAYSFCGTVEYMAPEVVN-RRGHTQSADWWSFGVLMFEMLTGTLPFQGK---DRKETMTMILKAKLGMPQFLS----PEAQSLLRMLFKRNPANRLGagpdgVEEIKRHSFF" == hsp.hit.seq
        assert "67899**********8888876655455****************999999****************************************************************************************************9***********************.******************************...999999999999988888866....9******************9888888999999886" == hsp.aln_annotation["PP"]
        assert "+ell+ lG+GsfGkV+ +kk + ++    +A+K+lkk + k + + +++ E  il +++Hp+ivkl+ +f+t+ +lyl+l++++ggdlf+ l+ke  ++ee++k+++ +++ +l++lHs gii+rDLKpeNiLld++g++k++DFGl+k+   +++k+++++gt eYmAPEv++ ++++t+++D+Ws+Gv+++e+ltg+lpf+g+   d++e++ +ilk kl  ++  s    +e+++l++ l++++pa+Rl      +eei++h+++" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[1]
        assert 2 == hsp.domain_index
        assert hsp.is_included
        assert 241.0 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 5.1e-75 == hsp.evalue_cond
        assert 3.9e-70 == hsp.evalue
        assert 0 == hsp.query_start
        assert 260 == hsp.query_end
        assert "[]" == hsp.query_endtype
        assert 421 == hsp.hit_start
        assert 679 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 421 == hsp.env_start
        assert 679 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.97 == hsp.acc_avg
        assert "EEEEEEEEEETTEEEEEEEETTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHH.SSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEE....EE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTTHHHHHTSGGG" == hsp.aln_annotation["CS"]
        assert "yelleklGsGsfGkVykakkkktgkkvAvKilkkeeekskkektavrElkilkkl.sHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgev....kiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRltaeeilkhpwl" == hsp.query.seq
        assert "YEVKEDIGVGSYSVCKRCIHKATNMEFAVKIIDKSKRDPTE------EIEILLRYgQHPNIITLKDVYDDGKYVYVVTELMKGGELLDKILRQKFFSEREASAVLFTITKTVEYLHAQGVVHRDLKPSNILYVDESGNpesiRICDFGFAKQLRAENGLLMTPCYTANFVAPEVLK-RQGYDAACDIWSLGVLLYTMLTGYTPFANGPDDTPEEILARIGSGKFSLSGGYWNSVSDTAKDLVSKMLHVDPHQRLTAALVLRHPWI" == hsp.hit.seq
        assert "7899***********************************88......9*******99**********************************9***************************************98554444888**********************************.***************************************************************************************7" == hsp.aln_annotation["PP"]
        assert "ye++e +G Gs+++++++++k t+ ++AvKi++k++++ ++      E++il +  +Hpni++l +v+ + +++y+v+e+++gg+l d + +++ +se+e+  ++ +i + +eylH +g++HrDLKp+NiL  +++      +i+DFG+ak+l  +++ l t + t +++APEvl+ +++y++++D+WslGv+ly++ltg +pf +  +++ +e++++i ++k + +   ++s s+++kdl++k+l+ dp++Rlta+ +l+hpw+" == hsp.aln_annotation["similarity"]

        # test if we've properly finished iteration
        # self.assertRaises(StopIteration, next, qresults)
        # self.assertEqual(1, counter)

    def test_30_hmmsearch_003(self):
        """Parsing hmmersearch 3.0 (text_30_hmmsearch_003)."""
        txt_file = get_file("text_30_hmmsearch_003.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1

        assert "hmmsearch" == qresult.program
        assert "/home/bow/db/hmmer/uniprot_sprot.fasta" == qresult.target
        assert "3.0" == qresult.version
        assert "Pkinase" == qresult.id
        assert "PF00069.17" == qresult.accession
        assert "Protein kinase domain" == qresult.description
        assert 260 == qresult.seq_len
        assert 7 == len(qresult)

        hit = qresult[0]
        assert "sp|Q9WUT3|KS6A2_MOUSE" == hit.id
        assert "Ribosomal protein S6 kinase alpha-2 OS=Mus musculus GN=Rps6ka2 PE=2 SV=1" == hit.description
        assert hit.is_included
        assert 8.4e-147 == hit.evalue
        assert 492.3 == hit.bitscore
        assert 0.0 == hit.bias
        assert 2.1 == hit.domain_exp_num
        assert 2 == hit.domain_obs_num
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 241.2 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 4.6e-75 == hsp.evalue_cond
        assert 3.5e-70 == hsp.evalue
        assert 0 == hsp.query_start
        assert 260 == hsp.query_end
        assert "[]" == hsp.query_endtype
        assert 58 == hsp.hit_start
        assert 318 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 58 == hsp.env_start
        assert 318 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.95 == hsp.acc_avg

        hsp = hit.hsps[1]
        assert 2 == hsp.domain_index
        assert hsp.is_included
        assert 249.3 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 1.5e-77 == hsp.evalue_cond
        assert 1.1e-72 == hsp.evalue
        assert 0 == hsp.query_start
        assert 260 == hsp.query_end
        assert "[]" == hsp.query_endtype
        assert 414 == hsp.hit_start
        assert 672 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 414 == hsp.env_start
        assert 672 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.97 == hsp.acc_avg

        hit = qresult[-1]
        assert "sp|P18654|KS6A3_MOUSE" == hit.id
        assert "Ribosomal protein S6 kinase alpha-3 OS=Mus musculus GN=Rps6ka3 PE=1 SV=2" == hit.description
        assert hit.is_included
        assert 5e-144 == hit.evalue
        assert 483.2 == hit.bitscore
        assert 0.0 == hit.bias
        assert 2.2 == hit.domain_exp_num
        assert 2 == hit.domain_obs_num
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 240.3 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 8.7e-75 == hsp.evalue_cond
        assert 6.6e-70 == hsp.evalue
        assert 0 == hsp.query_start
        assert 260 == hsp.query_end
        assert "[]" == hsp.query_endtype
        assert 67 == hsp.hit_start
        assert 327 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 67 == hsp.env_start
        assert 327 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.95 == hsp.acc_avg

        hsp = hit.hsps[1]
        assert 2 == hsp.domain_index
        assert hsp.is_included
        assert 241.0 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 5.1e-75 == hsp.evalue_cond
        assert 3.9e-70 == hsp.evalue
        assert 0 == hsp.query_start
        assert 260 == hsp.query_end
        assert "[]" == hsp.query_endtype
        assert 421 == hsp.hit_start
        assert 679 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 421 == hsp.env_start
        assert 679 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.97 == hsp.acc_avg

        # test if we've properly finished iteration
        # self.assertRaises(StopIteration, next, qresults)
        # self.assertEqual(1, counter)

    def test_30_hmmsearch_004(self):
        """Parsing hmmersearch 3.0 (text_30_hmmsearch_004)."""
        txt_file = get_file("text_30_hmmsearch_004.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1

        assert "hmmsearch" == qresult.program
        assert "/home/bow/db/hmmer/uniprot_sprot.fasta" == qresult.target
        assert "3.0" == qresult.version
        assert "Pkinase" == qresult.id
        assert "PF00069.17" == qresult.accession
        assert "Protein kinase domain" == qresult.description
        assert 260 == qresult.seq_len
        assert 7 == len(qresult)

        hit = qresult[0]
        assert "sp|Q9WUT3|KS6A2_MOUSE" == hit.id
        assert "Ribosomal protein S6 kinase alpha-2 OS=Mus musculus GN=Rps6ka2 PE=2 SV=1" == hit.description
        assert hit.is_included
        assert 8.4e-147 == hit.evalue
        assert 492.3 == hit.bitscore
        assert 0.0 == hit.bias
        assert 2.1 == hit.domain_exp_num
        assert 2 == hit.domain_obs_num
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 241.2 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 4.6e-75 == hsp.evalue_cond
        assert 3.5e-70 == hsp.evalue
        assert 0 == hsp.query_start
        assert 260 == hsp.query_end
        assert "[]" == hsp.query_endtype
        assert 58 == hsp.hit_start
        assert 318 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 58 == hsp.env_start
        assert 318 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.95 == hsp.acc_avg
        assert "EEEEEEEEEETTEEEEEEEE...TTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHHSSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEEEE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTT.....HHHHHTSGGG" == hsp.aln_annotation["CS"]
        assert "yelleklGsGsfGkVykakk...kktgkkvAvKilkkeeekskkektavrElkilkklsHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgevkiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRlt.....aeeilkhpwl" == hsp.query.seq
        assert "FELLKVLGQGSYGKVFLVRKvtgSDAGQLYAMKVLKKATLKVRDRVRSKMERDILAEVNHPFIVKLHYAFQTEGKLYLILDFLRGGDLFTRLSKEVMFTEEDVKFYLAELALALDHLHGLGIIYRDLKPENILLDEEGHIKITDFGLSKEATDHDKRAYSFCGTIEYMAPEVVN-RRGHTQSADWWSFGVLMFEMLTGSLPFQGK---DRKETMALILKAKLGMPQFLS----AEAQSLLRALFKRNPCNRLGagvdgVEEIKRHPFF" == hsp.hit.seq
        assert "67899**********7666611155667*****************99999****************************************************************************************************************************.******************************...999999999999999998866....99******************9999999*****997" == hsp.aln_annotation["PP"]
        assert "+ell+ lG+Gs+GkV+ ++k   ++ g+ +A+K+lkk + k + + +++ E  il++++Hp+ivkl+ +f+t+ +lyl+l++++ggdlf+ l+ke  ++ee++k+++ +++ +l++lH  gii+rDLKpeNiLld++g++ki+DFGl+k+++ +++++++++gt eYmAPEv++ ++++t+++D+Ws+Gv+++e+ltg+lpf+g+   d++e+++ ilk kl  ++  s    +e+++l++ l++++p +Rl      +eei++hp++" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[1]
        assert 2 == hsp.domain_index
        assert hsp.is_included
        assert 249.3 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 1.5e-77 == hsp.evalue_cond
        assert 1.1e-72 == hsp.evalue
        assert 0 == hsp.query_start
        assert 260 == hsp.query_end
        assert "[]" == hsp.query_endtype
        assert 414 == hsp.hit_start
        assert 672 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 414 == hsp.env_start
        assert 672 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.97 == hsp.acc_avg
        assert "EEEEEEEEEETTEEEEEEEETTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHH.SSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEE....EE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTTHHHHHTSGGG" == hsp.aln_annotation["CS"]
        assert "yelleklGsGsfGkVykakkkktgkkvAvKilkkeeekskkektavrElkilkkl.sHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgev....kiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRltaeeilkhpwl" == hsp.query.seq
        assert "YEIKEDIGVGSYSVCKRCVHKATDAEYAVKIIDKSKRDPSE------EIEILLRYgQHPNIITLKDVYDDGKYVYLVMELMRGGELLDRILRQRCFSEREASDVLYTIARTMDYLHSQGVVHRDLKPSNILYMDESGNpesiRICDFGFAKQLRAENGLLMTPCYTANFVAPEVLK-RQGYDAACDVWSLGILLYTMLAGFTPFANGPDDTPEEILARIGSGKYALSGGNWDSISDAAKDVVSKMLHVDPQQRLTAVQVLKHPWI" == hsp.hit.seq
        assert "7899***********************************98......9*******99**************************************************************************98544444888**********************************.***************************************************************************************7" == hsp.aln_annotation["PP"]
        assert "ye++e +G Gs+++++++++k t  ++AvKi++k++++ ++      E++il +  +Hpni++l +v+ + +++ylv+e+++gg+l d + +++++se+e+ +++++i++ ++ylHs+g++HrDLKp+NiL  +++      +i+DFG+ak+l  +++ l t + t +++APEvl+ +++y++++DvWslG++ly++l g +pf +  +++ +e++++i ++k+  +  +++s s+ +kd+++k+l+ dp++Rlta ++lkhpw+" == hsp.aln_annotation["similarity"]

        hit = qresult[-1]
        assert "sp|P18654|KS6A3_MOUSE" == hit.id
        assert "Ribosomal protein S6 kinase alpha-3 OS=Mus musculus GN=Rps6ka3 PE=1 SV=2" == hit.description
        assert hit.is_included
        assert 5e-144 == hit.evalue
        assert 483.2 == hit.bitscore
        assert 0.0 == hit.bias
        assert 2.2 == hit.domain_exp_num
        assert 2 == hit.domain_obs_num
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 240.3 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 8.7e-75 == hsp.evalue_cond
        assert 6.6e-70 == hsp.evalue
        assert 0 == hsp.query_start
        assert 260 == hsp.query_end
        assert "[]" == hsp.query_endtype
        assert 67 == hsp.hit_start
        assert 327 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 67 == hsp.env_start
        assert 327 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.95 == hsp.acc_avg
        assert "EEEEEEEEEETTEEEEEEEETTTTE...EEEEEEEEHHHCCCCCCHHHHHHHHHHHHHSSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEEEE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTT.....HHHHHTSGGG" == hsp.aln_annotation["CS"]
        assert "yelleklGsGsfGkVykakkkktgk...kvAvKilkkeeekskkektavrElkilkklsHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgevkiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRlt.....aeeilkhpwl" == hsp.query.seq
        assert "FELLKVLGQGSFGKVFLVKKISGSDarqLYAMKVLKKATLKVRDRVRTKMERDILVEVNHPFIVKLHYAFQTEGKLYLILDFLRGGDLFTRLSKEVMFTEEDVKFYLAELALALDHLHSLGIIYRDLKPENILLDEEGHIKLTDFGLSKESIDHEKKAYSFCGTVEYMAPEVVN-RRGHTQSADWWSFGVLMFEMLTGTLPFQGK---DRKETMTMILKAKLGMPQFLS----PEAQSLLRMLFKRNPANRLGagpdgVEEIKRHSFF" == hsp.hit.seq
        assert "67899**********8888876655455****************999999****************************************************************************************************9***********************.******************************...999999999999988888866....9******************9888888999999886" == hsp.aln_annotation["PP"]
        assert "+ell+ lG+GsfGkV+ +kk + ++    +A+K+lkk + k + + +++ E  il +++Hp+ivkl+ +f+t+ +lyl+l++++ggdlf+ l+ke  ++ee++k+++ +++ +l++lHs gii+rDLKpeNiLld++g++k++DFGl+k+   +++k+++++gt eYmAPEv++ ++++t+++D+Ws+Gv+++e+ltg+lpf+g+   d++e++ +ilk kl  ++  s    +e+++l++ l++++pa+Rl      +eei++h+++" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[1]
        assert 2 == hsp.domain_index
        assert hsp.is_included
        assert 241.0 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 5.1e-75 == hsp.evalue_cond
        assert 3.9e-70 == hsp.evalue
        assert 0 == hsp.query_start
        assert 260 == hsp.query_end
        assert "[]" == hsp.query_endtype
        assert 421 == hsp.hit_start
        assert 679 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 421 == hsp.env_start
        assert 679 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.97 == hsp.acc_avg
        assert "EEEEEEEEEETTEEEEEEEETTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHH.SSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEE....EE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTTHHHHHTSGGG" == hsp.aln_annotation["CS"]
        assert "yelleklGsGsfGkVykakkkktgkkvAvKilkkeeekskkektavrElkilkkl.sHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgev....kiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRltaeeilkhpwl" == hsp.query.seq
        assert "YEVKEDIGVGSYSVCKRCIHKATNMEFAVKIIDKSKRDPTE------EIEILLRYgQHPNIITLKDVYDDGKYVYVVTELMKGGELLDKILRQKFFSEREASAVLFTITKTVEYLHAQGVVHRDLKPSNILYVDESGNpesiRICDFGFAKQLRAENGLLMTPCYTANFVAPEVLK-RQGYDAACDIWSLGVLLYTMLTGYTPFANGPDDTPEEILARIGSGKFSLSGGYWNSVSDTAKDLVSKMLHVDPHQRLTAALVLRHPWI" == hsp.hit.seq
        assert "7899***********************************88......9*******99**********************************9***************************************98554444888**********************************.***************************************************************************************7" == hsp.aln_annotation["PP"]
        assert "ye++e +G Gs+++++++++k t+ ++AvKi++k++++ ++      E++il +  +Hpni++l +v+ + +++y+v+e+++gg+l d + +++ +se+e+  ++ +i + +eylH +g++HrDLKp+NiL  +++      +i+DFG+ak+l  +++ l t + t +++APEvl+ +++y++++D+WslGv+ly++ltg +pf +  +++ +e++++i ++k + +   ++s s+++kdl++k+l+ dp++Rlta+ +l+hpw+" == hsp.aln_annotation["similarity"]

        # test if we've properly finished iteration
        # self.assertRaises(StopIteration, next, qresults)
        # self.assertEqual(1, counter)

    def test_30_hmmsearch_005(self):
        """Parsing hmmersearch 3.0 (text_30_hmmsearch_005)."""
        txt_file = get_file("text_30_hmmsearch_005.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        # test first qresult
        qresult = next(qresults)
        counter += 1

        assert "hmmsearch" == qresult.program
        assert "/home/bow/db/hmmer/uniprot_sprot.fasta" == qresult.target
        assert "3.0" == qresult.version
        assert "globins4" == qresult.id
        assert 149 == qresult.seq_len
        assert 0 == len(qresult)

        # test second qresult
        qresult = next(qresults)
        counter += 1

        assert "hmmsearch" == qresult.program
        assert "/home/bow/db/hmmer/uniprot_sprot.fasta" == qresult.target
        assert "3.0" == qresult.version
        assert "Pkinase" == qresult.id
        assert "PF00069.17" == qresult.accession
        assert "Protein kinase domain" == qresult.description
        assert 260 == qresult.seq_len
        assert 7 == len(qresult)

        hit = qresult[0]
        assert "sp|Q9WUT3|KS6A2_MOUSE" == hit.id
        assert "Ribosomal protein S6 kinase alpha-2 OS=Mus musculus GN=Rps6ka2 PE=2 SV=1" == hit.description
        assert hit.is_included
        assert 8.4e-147 == hit.evalue
        assert 492.3 == hit.bitscore
        assert 0.0 == hit.bias
        assert 2.1 == hit.domain_exp_num
        assert 2 == hit.domain_obs_num
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 241.2 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 4.6e-75 == hsp.evalue_cond
        assert 3.5e-70 == hsp.evalue
        assert 0 == hsp.query_start
        assert 260 == hsp.query_end
        assert "[]" == hsp.query_endtype
        assert 58 == hsp.hit_start
        assert 318 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 58 == hsp.env_start
        assert 318 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.95 == hsp.acc_avg
        assert "EEEEEEEEEETTEEEEEEEE...TTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHHSSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEEEE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTT.....HHHHHTSGGG" == hsp.aln_annotation["CS"]
        assert "yelleklGsGsfGkVykakk...kktgkkvAvKilkkeeekskkektavrElkilkklsHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgevkiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRlt.....aeeilkhpwl" == hsp.query.seq
        assert "FELLKVLGQGSYGKVFLVRKvtgSDAGQLYAMKVLKKATLKVRDRVRSKMERDILAEVNHPFIVKLHYAFQTEGKLYLILDFLRGGDLFTRLSKEVMFTEEDVKFYLAELALALDHLHGLGIIYRDLKPENILLDEEGHIKITDFGLSKEATDHDKRAYSFCGTIEYMAPEVVN-RRGHTQSADWWSFGVLMFEMLTGSLPFQGK---DRKETMALILKAKLGMPQFLS----AEAQSLLRALFKRNPCNRLGagvdgVEEIKRHPFF" == hsp.hit.seq
        assert "67899**********7666611155667*****************99999****************************************************************************************************************************.******************************...999999999999999998866....99******************9999999*****997" == hsp.aln_annotation["PP"]
        assert "+ell+ lG+Gs+GkV+ ++k   ++ g+ +A+K+lkk + k + + +++ E  il++++Hp+ivkl+ +f+t+ +lyl+l++++ggdlf+ l+ke  ++ee++k+++ +++ +l++lH  gii+rDLKpeNiLld++g++ki+DFGl+k+++ +++++++++gt eYmAPEv++ ++++t+++D+Ws+Gv+++e+ltg+lpf+g+   d++e+++ ilk kl  ++  s    +e+++l++ l++++p +Rl      +eei++hp++" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[1]
        assert 2 == hsp.domain_index
        assert hsp.is_included
        assert 249.3 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 1.5e-77 == hsp.evalue_cond
        assert 1.1e-72 == hsp.evalue
        assert 0 == hsp.query_start
        assert 260 == hsp.query_end
        assert "[]" == hsp.query_endtype
        assert 414 == hsp.hit_start
        assert 672 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 414 == hsp.env_start
        assert 672 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.97 == hsp.acc_avg
        assert "EEEEEEEEEETTEEEEEEEETTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHH.SSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEE....EE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTTHHHHHTSGGG" == hsp.aln_annotation["CS"]
        assert "yelleklGsGsfGkVykakkkktgkkvAvKilkkeeekskkektavrElkilkkl.sHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgev....kiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRltaeeilkhpwl" == hsp.query.seq
        assert "YEIKEDIGVGSYSVCKRCVHKATDAEYAVKIIDKSKRDPSE------EIEILLRYgQHPNIITLKDVYDDGKYVYLVMELMRGGELLDRILRQRCFSEREASDVLYTIARTMDYLHSQGVVHRDLKPSNILYMDESGNpesiRICDFGFAKQLRAENGLLMTPCYTANFVAPEVLK-RQGYDAACDVWSLGILLYTMLAGFTPFANGPDDTPEEILARIGSGKYALSGGNWDSISDAAKDVVSKMLHVDPQQRLTAVQVLKHPWI" == hsp.hit.seq
        assert "7899***********************************98......9*******99**************************************************************************98544444888**********************************.***************************************************************************************7" == hsp.aln_annotation["PP"]
        assert "ye++e +G Gs+++++++++k t  ++AvKi++k++++ ++      E++il +  +Hpni++l +v+ + +++ylv+e+++gg+l d + +++++se+e+ +++++i++ ++ylHs+g++HrDLKp+NiL  +++      +i+DFG+ak+l  +++ l t + t +++APEvl+ +++y++++DvWslG++ly++l g +pf +  +++ +e++++i ++k+  +  +++s s+ +kd+++k+l+ dp++Rlta ++lkhpw+" == hsp.aln_annotation["similarity"]

        hit = qresult[-1]
        assert "sp|P18654|KS6A3_MOUSE" == hit.id
        assert "Ribosomal protein S6 kinase alpha-3 OS=Mus musculus GN=Rps6ka3 PE=1 SV=2" == hit.description
        assert hit.is_included
        assert 5e-144 == hit.evalue
        assert 483.2 == hit.bitscore
        assert 0.0 == hit.bias
        assert 2.2 == hit.domain_exp_num
        assert 2 == hit.domain_obs_num
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 240.3 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 8.7e-75 == hsp.evalue_cond
        assert 6.6e-70 == hsp.evalue
        assert 0 == hsp.query_start
        assert 260 == hsp.query_end
        assert "[]" == hsp.query_endtype
        assert 67 == hsp.hit_start
        assert 327 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 67 == hsp.env_start
        assert 327 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.95 == hsp.acc_avg
        assert "EEEEEEEEEETTEEEEEEEETTTTE...EEEEEEEEHHHCCCCCCHHHHHHHHHHHHHSSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEEEE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTT.....HHHHHTSGGG" == hsp.aln_annotation["CS"]
        assert "yelleklGsGsfGkVykakkkktgk...kvAvKilkkeeekskkektavrElkilkklsHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgevkiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRlt.....aeeilkhpwl" == hsp.query.seq
        assert "FELLKVLGQGSFGKVFLVKKISGSDarqLYAMKVLKKATLKVRDRVRTKMERDILVEVNHPFIVKLHYAFQTEGKLYLILDFLRGGDLFTRLSKEVMFTEEDVKFYLAELALALDHLHSLGIIYRDLKPENILLDEEGHIKLTDFGLSKESIDHEKKAYSFCGTVEYMAPEVVN-RRGHTQSADWWSFGVLMFEMLTGTLPFQGK---DRKETMTMILKAKLGMPQFLS----PEAQSLLRMLFKRNPANRLGagpdgVEEIKRHSFF" == hsp.hit.seq
        assert "67899**********8888876655455****************999999****************************************************************************************************9***********************.******************************...999999999999988888866....9******************9888888999999886" == hsp.aln_annotation["PP"]
        assert "+ell+ lG+GsfGkV+ +kk + ++    +A+K+lkk + k + + +++ E  il +++Hp+ivkl+ +f+t+ +lyl+l++++ggdlf+ l+ke  ++ee++k+++ +++ +l++lHs gii+rDLKpeNiLld++g++k++DFGl+k+   +++k+++++gt eYmAPEv++ ++++t+++D+Ws+Gv+++e+ltg+lpf+g+   d++e++ +ilk kl  ++  s    +e+++l++ l++++pa+Rl      +eei++h+++" == hsp.aln_annotation["similarity"]

        hsp = hit.hsps[1]
        assert 2 == hsp.domain_index
        assert hsp.is_included
        assert 241.0 == hsp.bitscore
        assert 0.0 == hsp.bias
        assert 5.1e-75 == hsp.evalue_cond
        assert 3.9e-70 == hsp.evalue
        assert 0 == hsp.query_start
        assert 260 == hsp.query_end
        assert "[]" == hsp.query_endtype
        assert 421 == hsp.hit_start
        assert 679 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 421 == hsp.env_start
        assert 679 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.97 == hsp.acc_avg
        assert "EEEEEEEEEETTEEEEEEEETTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHH.SSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEE....EE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTTHHHHHTSGGG" == hsp.aln_annotation["CS"]
        assert "yelleklGsGsfGkVykakkkktgkkvAvKilkkeeekskkektavrElkilkkl.sHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgev....kiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRltaeeilkhpwl" == hsp.query.seq
        assert "YEVKEDIGVGSYSVCKRCIHKATNMEFAVKIIDKSKRDPTE------EIEILLRYgQHPNIITLKDVYDDGKYVYVVTELMKGGELLDKILRQKFFSEREASAVLFTITKTVEYLHAQGVVHRDLKPSNILYVDESGNpesiRICDFGFAKQLRAENGLLMTPCYTANFVAPEVLK-RQGYDAACDIWSLGVLLYTMLTGYTPFANGPDDTPEEILARIGSGKFSLSGGYWNSVSDTAKDLVSKMLHVDPHQRLTAALVLRHPWI" == hsp.hit.seq
        assert "7899***********************************88......9*******99**********************************9***************************************98554444888**********************************.***************************************************************************************7" == hsp.aln_annotation["PP"]
        assert "ye++e +G Gs+++++++++k t+ ++AvKi++k++++ ++      E++il +  +Hpni++l +v+ + +++y+v+e+++gg+l d + +++ +se+e+  ++ +i + +eylH +g++HrDLKp+NiL  +++      +i+DFG+ak+l  +++ l t + t +++APEvl+ +++y++++D+WslGv+ly++ltg +pf +  +++ +e++++i ++k + +   ++s s+++kdl++k+l+ dp++Rlta+ +l+hpw+" == hsp.aln_annotation["similarity"]

        # test if we've properly finished iteration
        with pytest.raises(StopIteration):
            next(qresults)
        assert 2 == counter


class PhmmerCases(unittest.TestCase):
    """Testing phmmer output."""

    def test_31b2_phmmer_001(self):
        """Parsing phmmer 3.1b2 (text_31b2_phmmer_001)."""
        txt_file = get_file("text_31b2_phmmer_001.out")
        qresults = list(parse(txt_file, FMT))

        # first qresult
        qresult = qresults[0]

        assert "phmmer" == qresult.program
        assert "/home/bow/devel/sandbox/biopython-sandbox/db/hmmer/protdb/uniprot_sprot.fasta" == qresult.target
        assert "3.1b2" == qresult.version
        assert "sp|Q6GZX4|001R_FRG3G" == qresult.id
        assert 256 == qresult.seq_len
        assert 13 == len(qresult)

        hit = qresult[0]
        assert "sp|Q6GZX4|001R_FRG3G" == hit.id
        assert "Putative transcription factor 001R OS=Frog virus 3 (isolate Goorha) GN=FV3-001R PE=4 SV=1" == hit.description
        assert hit.is_included
        assert 1.1e-176 == hit.evalue
        assert 590.0 == hit.bitscore
        assert 1.4 == hit.bias
        assert 1.0 == hit.domain_exp_num
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert 1 == hsp.domain_index
        assert hsp.is_included
        assert 589.9 == hsp.bitscore
        assert 1.4 == hsp.bias
        assert 3e-181 == hsp.evalue_cond
        assert 1.3e-176 == hsp.evalue
        assert 0 == hsp.hit_start
        assert 256 == hsp.hit_end
        assert "[]" == hsp.hit_endtype
        assert 0 == hsp.query_start
        assert 256 == hsp.query_end
        assert "[]" == hsp.query_endtype
        assert 0 == hsp.env_start
        assert 256 == hsp.env_end
        assert "[]" == hsp.env_endtype
        assert 1.00 == hsp.acc_avg
        assert "mafsaedvlkeydrrrrmealllslyypndrklldykewspprvqvecpkapvewnnppsekglivghfsgikykgekaqasevdvnkm" == hsp.query.seq[:89]
        assert "MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPKAPVEWNNPPSEKGLIVGHFSGIKYKGEKAQASEVDVNKM" == hsp.hit.seq[:89]
        assert "89***************************************************************************************" == hsp.aln_annotation["PP"][:89]
        assert "mafsaedvlkeydrrrrmealllslyypndrklldykewspprvqvecpkapvewnnppsekglivghfsgikykgekaqasevdvnkm" == hsp.aln_annotation["similarity"][:89]

        # last query, last hit
        qresult = qresults[-1]

        assert "phmmer" == qresult.program
        assert "/home/bow/devel/sandbox/biopython-sandbox/db/hmmer/protdb/uniprot_sprot.fasta" == qresult.target
        assert "3.1b2" == qresult.version
        assert "sp|Q197F7|003L_IIV3" == qresult.id
        assert 156 == qresult.seq_len
        assert 6 == len(qresult)

        hit = qresult[-1]
        assert "sp|P04060|RNAS1_HYSCR" == hit.id
        assert "Ribonuclease pancreatic OS=Hystrix cristata GN=RNASE1 PE=1 SV=1" == hit.description
        assert not hit.is_included
        assert 2.1 == hit.evalue
        assert 13.1 == hit.bitscore
        assert 0.7 == hit.bias
        assert 1.8 == hit.domain_exp_num
        assert 2 == hit.domain_obs_num
        assert 2 == len(hit)

        hsp = hit.hsps[-1]
        assert 2 == hsp.domain_index
        assert not hsp.is_included
        assert 12.1 == hsp.bitscore
        assert 0.1 == hsp.bias
        assert 4.7e-05 == hsp.evalue_cond
        assert 4.2 == hsp.evalue
        assert 91 == hsp.hit_start
        assert 121 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 7 == hsp.query_start
        assert 37 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 84 == hsp.env_start
        assert 127 == hsp.env_end
        assert ".." == hsp.env_endtype
        assert 0.84 == hsp.acc_avg
        assert "cpqswygspqlereivckmsgaphypnyyp" == hsp.query.seq
        assert "YPDCSYGMSQLERSIVVACEGSPYVPVHFD" == hsp.hit.seq
        assert "68889******************9887764" == hsp.aln_annotation["PP"]
        assert " p   yg  qler iv    g+p+ p ++ " == hsp.aln_annotation["similarity"]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
