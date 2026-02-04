# Copyright 2012 by Kai Blin.
# Revisions copyright 2012 by Wibowo Arindrarto
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO hmmer2 text module."""

import unittest
import pytest
from os import path

from Bio.SearchIO import parse
from Bio.SearchIO import read


class HmmpfamTests(unittest.TestCase):
    fmt = "hmmer2-text"

    def test_hmmpfam_21(self):
        """Test parsing hmmpfam 2.1 file (text_21_hmmpfam_001.out)."""
        results = parse(path.join("Hmmer", "text_21_hmmpfam_001.out"), self.fmt)
        res = next(results)
        assert "roa1_drome" == res.id
        assert "<unknown description>" == res.description
        assert "hmmpfam" == res.program
        assert "2.1.1" == res.version
        assert "pfam" == res.target
        assert 1 == len(res)

        hit = res[0]
        assert "SEED" == hit.id
        assert "<unknown description>" == hit.description
        assert 146.1 == pytest.approx(hit.bitscore, abs=5e-8)
        assert 6.3e-40 == pytest.approx(hit.evalue, abs=5e-42)
        assert 2 == hit.domain_obs_num
        assert 2 == len(hit)

        hsp = hit[0]
        assert 1 == hsp.domain_index
        assert 0 == hsp.hit_start
        assert 77 == hsp.hit_end
        assert "[]" == hsp.hit_endtype
        assert 32 == hsp.query_start
        assert 103 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 71.2 == pytest.approx(hsp.bitscore, abs=5e-8)
        assert 2.2e-17 == pytest.approx(hsp.evalue, abs=5e-19)
        assert "lfVgNLppdvteedLkdlFskfGpivsikivrDiiekpketgkskGfaFVeFeseedAekAlealnG.kelggrklrv" == hsp.hit.seq
        assert "lf+g+L + +t+e Lk++F+k G iv++ +++D     + t++s+Gf+F+++  ++  + A +    +++++gr+++ " == str(hsp.aln_annotation["similarity"])
        assert "LFIGGLDYRTTDENLKAHFEKWGNIVDVVVMKD-----PRTKRSRGFGFITYSHSSMIDEAQK--SRpHKIDGRVVEP" == hsp.query.seq

        hsp = hit[1]
        assert 2 == hsp.domain_index
        assert 0 == hsp.hit_start
        assert 77 == hsp.hit_end
        assert "[]" == hsp.hit_endtype
        assert 123 == hsp.query_start
        assert 194 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 75.5 == pytest.approx(hsp.bitscore, abs=5e-8)
        assert 1.1e-18 == pytest.approx(hsp.evalue, abs=5e-20)
        assert "lfVgNLppdvteedLkdlFskfGpivsikivrDiiekpketgkskGfaFVeFeseedAekAlealnGkelggrklrv" == hsp.hit.seq
        assert "lfVg L  d +e+ ++d+F++fG iv+i+iv+D     ketgk +GfaFVeF++++ ++k +     ++l+g+ + v" == str(hsp.aln_annotation["similarity"])
        assert "LFVGALKDDHDEQSIRDYFQHFGNIVDINIVID-----KETGKKRGFAFVEFDDYDPVDKVVL-QKQHQLNGKMVDV" == hsp.query.seq

    def test_hmmpfam_22(self):
        """Test parsing hmmpfam 2.2 file (text_22_hmmpfam_001.out)."""
        results = parse(path.join("Hmmer", "text_22_hmmpfam_001.out"), self.fmt)
        res = next(results)
        assert "gi|1522636|gb|AAC37060.1|" == res.id
        assert "M. jannaschii predicted coding region MJECS02 [Methanococcus jannaschii]" == res.description
        assert "[none]" == res.accession
        assert "hmmpfam" == res.program
        assert "2.2g" == res.version
        assert "Pfam" == res.target
        assert 1 == len(res)

        hit = res[0]
        assert "Methylase_M" == hit.id
        assert "Type I restriction modification system, M" == hit.description
        assert -105.2 == pytest.approx(hit.bitscore, abs=5e-8)
        assert 0.0022 == pytest.approx(hit.evalue, abs=5e-8)
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        hsp = hit[0]
        assert 1 == hsp.domain_index
        assert 0 == hsp.hit_start
        assert 279 == hsp.hit_end
        assert "[]" == hsp.hit_endtype
        assert 279 == hsp.query_start
        assert 481 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert -105.2 == pytest.approx(hsp.bitscore, abs=5e-8)
        assert 0.0022 == pytest.approx(hsp.evalue, abs=5e-8)
        assert ("lrnELentLWavADkLRGsmDaseYKdyVLGLlFlKYiSdkFlerrieieerktdtesepsldyakledqyeql"
            "ededlekedfyqkkGvFilPsqlFwdfikeaeknkldedigtdldkifseledqialgypaSeedfkGlfpdld"
            "fnsnkLgskaqarnetLtelidlfselelgtPmHNG.dfeelgikDlfGDaYEYLLgkFAeneGKsGGeFYTPq"
            "eVSkLiaeiLtigqpsegdfsIYDPAcGSGSLllqaskflgehdgkrnaisyYGQEsn" == hsp.hit.seq)
        assert (" ++EL+++  av+   R              L+F K++ dk      +i+         p +   + +++y   "
            "++   ++ ++y ++      + lF++++   e ++  ++++ + +    ++      + +       Glf +++"
            "+  ++ +s+   +ne ++e+i+ +++ +++     G++ +el   D++G +YE L+   Ae   K+ G +YTP "
            "e++  ia+ + i+  ++                  +++ ++    k+n+i +    s+" == str(hsp.aln_annotation["similarity"]))
        assert ("NTSELDKKKFAVLLMNR--------------LIFIKFLEDK------GIV---------PRDLLRRTYEDY---"
            "KKSNVLI-NYYDAY-L----KPLFYEVLNTPEDER--KENIRT-NPYYKDIPYL---N-G-------GLFRSNN"
            "V--PNELSFTIKDNEIIGEVINFLERYKFTLSTSEGsEEVELNP-DILGYVYEKLINILAEKGQKGLGAYYTPD"
            "EITSYIAKNT-IEPIVVE----------------RFKEIIK--NWKINDINF----ST" == hsp.query.seq)

    def test_hmmpfam_23(self):
        """Test parsing hmmpfam 2.3 file (text_23_hmmpfam_001.out)."""
        results = parse(path.join("Hmmer", "text_23_hmmpfam_001.out"), self.fmt)
        res = next(results)
        assert "gi|90819130|dbj|BAE92499.1|" == res.id
        assert "glutamate synthase [Porphyra yezoensis]" == res.description
        assert "[none]" == res.accession
        assert "hmmpfam" == res.program
        assert "2.3.2" == res.version
        assert "../Shared/Pfam_fs" == res.target
        assert 54 == len(res)

        hit = res[0]
        assert "Glu_synthase" == hit.id
        assert "Conserved region in glutamate synthas" == hit.description
        assert 858.6 == pytest.approx(hit.bitscore, abs=5e-8)
        assert 3.6e-255 == pytest.approx(hit.evalue, abs=5e-257)
        assert 2 == hit.domain_obs_num
        assert 2 == len(hit)

        hsp = hit[0]
        assert 1 == hsp.domain_index
        assert 296 == hsp.hit_start
        assert 323 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert 649 == hsp.query_start
        assert 676 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 1.3 == pytest.approx(hsp.bitscore, abs=5e-8)
        assert 3 == pytest.approx(hsp.evalue, abs=5e-8)
        assert "lPwelgLaevhqtLvengLRdrVsLia" == hsp.hit.seq
        assert "+P  l++ +vh  L++ gLR + s+ +" == str(hsp.aln_annotation["similarity"])
        assert "IPPLLAVGAVHHHLINKGLRQEASILV" == hsp.query.seq

        hsp = hit[1]
        assert 2 == hsp.domain_index
        assert 0 == hsp.hit_start
        assert 412 == hsp.hit_end
        assert "[]" == hsp.hit_endtype
        assert 829 == hsp.query_start
        assert 1216 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert 857.3 == pytest.approx(hsp.bitscore, abs=5e-8)
        assert 9e-255 == pytest.approx(hsp.evalue, abs=5e-256)

    def test_hmmpfam_23_no_match(self):
        """Test parsing hmmpfam 2.3 file (text_23_hmmpfam_002.out)."""
        results = parse(path.join("Hmmer", "text_23_hmmpfam_002.out"), self.fmt)
        res = next(results)

        assert "SEQ0001" == res.id
        assert 0 == len(res.hits)

        res = next(results)

        assert "SEQ0002" == res.id
        assert 0 == len(res.hits)

    def test_hmmpfam_23_missing_consensus(self):
        """Test parsing hmmpfam 2.3 file (text_23_hmmpfam_003.out)."""
        results = parse(path.join("Hmmer", "text_23_hmmpfam_003.out"), self.fmt)
        res = next(results)

        assert "small_input" == res.id
        assert "[none]" == res.description
        assert "[none]" == res.accession
        assert "hmmpfam" == res.program
        assert "2.3.2" == res.version
        assert "antismash/specific_modules/lantipeptides/ClassIVLanti.hmm" == res.target
        assert 1 == len(res)

        hit = res[0]
        assert "ClassIVLanti" == hit.id
        assert "Class-IV" == hit.description
        assert -79.3 == pytest.approx(hit.bitscore, abs=5e-8)
        assert 1 == pytest.approx(hit.evalue, abs=5e-8)
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        hsp = hit[0]
        assert 1 == hsp.domain_index
        assert 0 == hsp.hit_start
        assert 66 == hsp.hit_end
        assert "[]" == hsp.hit_endtype
        assert 5 == hsp.query_start
        assert 20 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert -79.3 == pytest.approx(hsp.bitscore, abs=5e-8)
        assert 1 == pytest.approx(hsp.evalue, abs=5e-8)
        assert len(hsp.query.seq) == len(hsp.hit.seq)
        assert len(hsp.query.seq) == len(hsp.aln_annotation["similarity"])
        assert "msEEqLKAFiAKvqaDtsLqEqLKaEGADvvaiAKAaGFtitteDLnahiqakeLsdeeLEgvaGg" == hsp.hit.seq
        assert "        F+                           G  +t   Ln                   " == str(hsp.aln_annotation["similarity"])
        assert "-------CFL---------------------------GCLVTNWVLNRS-----------------" == hsp.query.seq

    def test_hmmpfam_23_break_in_end_of_seq(self):
        """Test parsing hmmpfam 2.3 file with a line break in the end of seq marker.

        file (text_23_hmmpfam_004.out)
        """
        results = parse(path.join("Hmmer", "text_23_hmmpfam_004.out"), self.fmt)
        res = next(results)
        assert "PKSI-KS" == res[0].id
        assert "PKSI-FK" == res[1].id

    def test_hmmpfam_24(self):
        """Test parsing hmmpfam 2.4 file (text_24_hmmpfam_001.out)."""
        results = list(parse(path.join("Hmmer", "text_24_hmmpfam_001.out"), self.fmt))
        assert 5 == len(results)

        # first qresult
        res = results[0]
        assert "random_s00" == res.id
        assert "[none]" == res.accession
        assert "[none]" == res.description
        assert "hmmpfam" == res.program
        assert "2.4i" == res.version
        assert "/home/bow/db/hmmer/Pfam_fs" == res.target
        assert 0 == len(res)

        # fourth qresult
        res = results[3]
        assert "gi|22748937|ref|NP_065801.1|" == res.id
        assert "[none]" == res.accession
        assert "exportin-5 [Homo sapiens]" == res.description
        assert "hmmpfam" == res.program
        assert "2.4i" == res.version
        assert "/home/bow/db/hmmer/Pfam_fs" == res.target
        assert 33 == len(res)

        # fourth qresult, first hit
        hit = res[0]
        assert "Xpo1" == hit.id
        assert "Exportin 1-like protein" == hit.description
        assert 170.1 == pytest.approx(hit.bitscore, abs=5e-8)
        assert 5.1e-48 == pytest.approx(hit.evalue, abs=5e-50)
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        # fourth qresult, first hit, first hsp
        hsp = hit[0]
        assert 1 == hsp.domain_index
        assert 170.1 == pytest.approx(hsp.bitscore, abs=5e-8)
        assert 5.1e-48 == pytest.approx(hsp.evalue, abs=5e-50)
        assert 108 == hsp.query_start
        assert 271 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert "ENHIKDALSRIVVEMIKREWPQHWPDMLIELDTLSKQG--" == hsp.query.seq[:40]
        assert "+++++  L+++++e++k+ewP++Wp+ + +l  l++++  " == str(hsp.aln_annotation["similarity"])[:40]
        assert "WVSMSHITA-ENCkLLEILCLLL----NEQELQLGAAECL" == hsp.query.seq[-40:]
        assert 0 == hsp.hit_start
        assert 178 == hsp.hit_end
        assert "[]" == hsp.hit_endtype
        assert "pkflrnKLalalaelakqewPsnWpsffpdlvsllsssss" == hsp.hit.seq[:40]
        assert "W+++++i + ++++ll++l+ lL    +  +l++ A+eCL" == str(hsp.aln_annotation["similarity"])[-40:]
        assert "Wipiglianvnpi.llnllfslLsgpesdpdlreaAveCL" == hsp.hit.seq[-40:]

        # fourth qresult, second from last hit
        hit = res[-2]
        assert "Rad50_zn_hook" == hit.id
        assert "Rad50 zinc hook motif" == hit.description
        assert 2.2 == pytest.approx(hit.bitscore, abs=5e-8)
        assert 9.2 == pytest.approx(hit.evalue, abs=5e-8)
        assert 2 == hit.domain_obs_num
        assert 2 == len(hit)

        # fourth qresult, second from last hit, first hsp
        hsp = hit[0]
        assert 1 == hsp.domain_index
        assert 0.8 == pytest.approx(hsp.bitscore, abs=5e-8)
        assert 22 == pytest.approx(hsp.evalue, abs=5e-8)
        assert 20 == hsp.query_start
        assert 47 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert "MDPNSTQRYRLEALKFCEEFKE-KCPIC" == hsp.query.seq
        assert 0 == hsp.hit_start
        assert 28 == hsp.hit_end
        assert "[." == hsp.hit_endtype
        assert "galesekaelkkaieeleeeesscCPvC" == hsp.hit.seq

        # fourth qresult, second from last hit, last hsp
        hsp = hit[-1]
        assert 2 == hsp.domain_index
        assert 1.3 == pytest.approx(hsp.bitscore, abs=5e-8)
        assert 16 == pytest.approx(hsp.evalue, abs=5e-8)
        assert 789 == hsp.query_start
        assert 811 == hsp.query_end
        assert ".." == hsp.query_endtype
        assert "EMLAKMAEPFTKALDMLDAEKS" == hsp.query.seq
        assert 0 == hsp.hit_start
        assert 22 == hsp.hit_end
        assert "[." == hsp.hit_endtype
        assert "galesekaelkkaieeleeees" == hsp.hit.seq


class HmmsearchTests(unittest.TestCase):
    fmt = "hmmer2-text"

    def test_hmmsearch_20(self):
        """Test parsing hmmsearch 2.0 file (text_20_hmmsearch_001.out)."""
        res = read(path.join("Hmmer", "text_20_hmmsearch_001.out"), self.fmt)

        # first query
        assert "SEED" == res.id
        assert "<unknown description>" == res.description
        assert "hmmsearch" == res.program
        assert "2.0" == res.version
        assert "HMM.dbtemp.29591" == res.target
        assert 751 == len(res)

        # first hit
        hit = res[0]
        assert "PAB2_ARATH" == hit.id
        assert "P42731 POLYADENYLATE-BINDING PROTEIN 2 (PO" == hit.description
        assert 393.8 == pytest.approx(hit.bitscore, abs=5e-8)
        assert 6.1e-114 == pytest.approx(hit.evalue, abs=5e-146)
        assert 4 == hit.domain_obs_num
        assert 4 == len(hit)

        # first hit, first hsp
        hsp = hit[0]
        assert "SEED" == hsp.query_id
        assert "<unknown description>" == hsp.query_description
        assert "PAB2_ARATH" == hsp.hit_id
        assert "P42731 POLYADENYLATE-BINDING PROTEIN 2 (PO" == hsp.hit_description
        assert 3 == hsp.domain_index
        assert 109.1 == pytest.approx(hsp.bitscore, abs=5e-8)
        assert 3e-28 == pytest.approx(hsp.evalue, abs=5e-29)
        assert 0 == hsp.query_start
        assert 77 == hsp.query_end
        assert "[]" == hsp.query_endtype
        assert 216 == hsp.hit_start
        assert 287 == hsp.hit_end
        assert ".." == hsp.hit_endtype

        # first hit, last hsp
        hsp = hit[-1]
        assert "SEED" == hsp.query_id
        assert "<unknown description>" == hsp.query_description
        assert "PAB2_ARATH" == hsp.hit_id
        assert "P42731 POLYADENYLATE-BINDING PROTEIN 2 (PO" == hsp.hit_description
        assert 1 == hsp.domain_index
        assert 92.1 == pytest.approx(hsp.bitscore, abs=5e-8)
        assert 3.9e-23 == pytest.approx(hsp.evalue, abs=5e-25)
        assert 0 == hsp.query_start
        assert 77 == hsp.query_end
        assert "[]" == hsp.query_endtype
        assert 37 == hsp.hit_start
        assert 109 == hsp.hit_end
        assert ".." == hsp.hit_endtype

        # last hit
        hit = res[-1]
        assert "O00369" == hit.id
        assert "O00369 L1 ELEMENT L1.20 P40 AND PUTATIVE P" == hit.description
        assert -23.8 == pytest.approx(hit.bitscore, abs=5e-8)
        assert 9.9e02 == pytest.approx(hit.evalue, abs=5e-8)
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        # last hit, first hsp
        hsp = hit[0]
        assert "SEED" == hsp.query_id
        assert "<unknown description>" == hsp.query_description
        assert "O00369" == hsp.hit_id
        assert "O00369 L1 ELEMENT L1.20 P40 AND PUTATIVE P" == hsp.hit_description
        assert 1 == hsp.domain_index
        assert -23.8 == pytest.approx(hsp.bitscore, abs=5e-8)
        assert 9.9e02 == pytest.approx(hsp.evalue, abs=5e-8)
        assert 0 == hsp.query_start
        assert 77 == hsp.query_end
        assert "[]" == hsp.query_endtype
        assert 180 == hsp.hit_start
        assert 249 == hsp.hit_end
        assert ".." == hsp.hit_endtype

    def test_hmmsearch_22(self):
        """Test parsing hmmsearch 2.2 file (text_22_hmmsearch_001.out)."""
        res = read(path.join("Hmmer", "text_22_hmmsearch_001.out"), self.fmt)

        # first query
        assert "Peptidase_C1" == res.id
        assert "Papain family cysteine protease" == res.description
        assert "hmmsearch" == res.program
        assert "2.2g" == res.version
        assert "cysprot1b.fa" == res.target
        assert 4 == len(res)

        # first hit
        hit = res[0]
        assert "CATL_RAT" == hit.id
        assert "<unknown description>" == hit.description
        assert 449.4 == pytest.approx(hit.bitscore, abs=5e-8)
        assert 2e-135 == pytest.approx(hit.evalue, abs=5e-136)
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        # first hit, first hsp
        hsp = hit[0]
        assert "Peptidase_C1" == hsp.query_id
        assert "Papain family cysteine protease" == hsp.query_description
        assert "CATL_RAT" == hit.id
        assert "<unknown description>" == hit.description
        assert 1 == hsp.domain_index
        assert 449.4 == pytest.approx(hsp.bitscore, abs=5e-8)
        assert 2e-135 == pytest.approx(hsp.evalue, abs=5e-136)
        assert 0 == hsp.query_start
        assert 337 == hsp.query_end
        assert "[]" == hsp.query_endtype
        assert "lPesfDWReWkggaVtpVKdQGiqCGSCWAFSavgalEgr" == hsp.query.seq[:40]
        assert "IVKNSWGtdWGEnGYfriaRgknksgkneCGIaseasypi" == hsp.query.seq[-40:]
        assert 337 == len(hsp.query.seq)
        assert "+P+++DWRe kg  VtpVK+QG qCGSCWAFSa g lEg+" == str(hsp.aln_annotation["similarity"])[:40]
        assert "+VKNSWG++WG++GY++ia+++n    n+CG+a+ asypi" == str(hsp.aln_annotation["similarity"])[-40:]
        assert 113 == hsp.hit_start
        assert 332 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert "IPKTVDWRE-KG-CVTPVKNQG-QCGSCWAFSASGCLEGQ" == hsp.hit.seq[:40]
        assert "LVKNSWGKEWGMDGYIKIAKDRN----NHCGLATAASYPI" == hsp.hit.seq[-40:]
        assert 337 == len(hsp.hit.seq)

        # last hit
        hit = res[-1]
        assert "PAPA_CARPA" == hit.id
        assert "<unknown description>" == hit.description
        assert 337.7 == pytest.approx(hit.bitscore, abs=5e-8)
        assert 9e-102 == pytest.approx(hit.evalue, abs=5e-103)
        assert 1 == hit.domain_obs_num
        assert 1 == len(hit)

        # last hit, last hsp
        hsp = hit[-1]
        assert "Peptidase_C1" == hsp.query_id
        assert "Papain family cysteine protease" == hsp.query_description
        assert "PAPA_CARPA" == hit.id
        assert "<unknown description>" == hit.description
        assert 1 == hsp.domain_index
        assert 337.7 == pytest.approx(hsp.bitscore, abs=5e-8)
        assert 9e-102 == pytest.approx(hsp.evalue, abs=5e-103)
        assert 0 == hsp.query_start
        assert 337 == hsp.query_end
        assert "[]" == hsp.query_endtype
        assert "lPesfDWReWkggaVtpVKdQGiqCGSCWAFSavgalEgr" == hsp.query.seq[:40]
        assert "IVKNSWGtdWGEnGYfriaRgknksgkneCGIaseasypi" == hsp.query.seq[-40:]
        assert 337 == len(hsp.query.seq)
        assert "+Pe +DWR+ kg aVtpVK+QG +CGSCWAFSav ++Eg+" == str(hsp.aln_annotation["similarity"])[:40]
        assert "++KNSWGt WGEnGY+ri+Rg+++s ++ CG+ ++  yp+" == str(hsp.aln_annotation["similarity"])[-40:]
        assert 133 == hsp.hit_start
        assert 343 == hsp.hit_end
        assert ".." == hsp.hit_endtype
        assert "IPEYVDWRQ-KG-AVTPVKNQG-SCGSCWAFSAVVTIEGI" == hsp.hit.seq[:40]
        assert "LIKNSWGTGWGENGYIRIKRGTGNS-YGVCGLYTSSFYPV" == hsp.hit.seq[-40:]
        assert 337 == len(hsp.hit.seq)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
