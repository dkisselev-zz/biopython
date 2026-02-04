# Copyright 2019 by Jens Thomas.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO HhsuiteIO parsers."""

import os
import unittest
import pytest

from Bio.SearchIO import parse

# test case files are in the Blast directory
TEST_DIR = "HHsuite"
FMT = "hhsuite2-text"


def get_file(filename):
    """Return the path of a test file."""
    return os.path.join(TEST_DIR, filename)


class HhsuiteCases(unittest.TestCase):
    """Test hhsuite2 output."""

    def test_2uvo(self):
        """Parsing 2uvo."""
        txt_file = get_file("2uvo_hhblits.hhr")
        qresults = parse(txt_file, FMT)

        # test first and only qresult
        qresult = next(qresults)

        num_hits = 16
        assert "HHSUITE" == qresult.program
        assert "2UVO:A|PDBID|CHAIN|SEQUENCE" == qresult.id
        assert 171 == qresult.seq_len
        assert num_hits == len(qresult)

        hit = qresult[0]
        assert "2uvo_A" == hit.id
        assert ("Agglutinin isolectin 1; carbohydrate-binding protein, hevein domain, chitin-binding,"
            " GERM agglutinin, chitin-binding protein; HET: NDG NAG GOL; 1.40A {Triticum aestivum}"
            " PDB: 1wgc_A* 2cwg_A* 2x3t_A* 4aml_A* 7wga_A 9wga_A 2wgc_A 1wgt_A 1k7t_A* 1k7v_A* 1k7u_A"
            " 2x52_A* 1t0w_A*" == hit.description)
        assert hit.is_included
        assert 3.7e-34 == hit.evalue
        assert 210.31 == hit.score
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert hsp.is_included
        assert 0 == hsp.output_index
        assert 99.95 == hsp.prob
        assert 210.31 == hsp.score
        assert 3.7e-34 == hsp.evalue
        assert 0 == hsp.hit_start
        assert 171 == hsp.hit_end
        assert 0 == hsp.query_start
        assert 171 == hsp.query_end
        assert ("ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYC"
            "GAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCS"
            "KWGSCGIGPGYCGAGCQSGGCDG" == hsp.hit.seq)
        assert ("ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYC"
            "GAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCS"
            "KWGSCGIGPGYCGAGCQSGGCDG" == hsp.query.seq)

        # Check last hit
        hit = qresult[num_hits - 1]
        assert "4z8i_A" == hit.id
        assert ("BBTPGRP3, peptidoglycan recognition protein 3; chitin-binding domain, "
            "AM hydrolase; 2.70A {Branchiostoma belcheri tsingtauense}" == hit.description)
        assert hit.is_included
        assert 0.11 == hit.evalue
        assert 36.29 == hit.score
        assert 2 == len(hit)

        # Check we can get the original last HSP from the file.
        num_hsps = 32
        assert num_hsps == len(qresult.hsps)

        hsp = qresult.hsps[-1]
        assert hsp.is_included
        assert num_hsps - 1 == hsp.output_index
        assert 2.6 == hsp.evalue
        assert 25.90 == hsp.score
        assert 40.43 == hsp.prob
        assert 10 == hsp.hit_start
        assert 116 == hsp.hit_end
        assert 53 == hsp.query_start
        assert 163 == hsp.query_end
        assert ("XCXXXXCCXXXXXCXXXXXXCXXXCXXXXCXXXXXCXXX--XXXCXXXXCCXXXXXCXXXXXXCXXXCXXXXCXXXXXCX"
            "XX--XXXCXXXXCCXXXXXCXXXXXXCXXX" == hsp.hit.seq)
        assert ("TCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCG"
            "KDAGGRVCTNNYCCSKWGSCGIGPGYCGAG" == hsp.query.seq)

    def test_2uvo_onlyheader(self):
        """Parsing 4uvo with only header present."""
        txt_file = get_file("2uvo_hhblits_onlyheader.hhr")
        qresults = parse(txt_file, FMT)

        with pytest.raises(RuntimeError):
            next(qresults)

    def test_2uvo_emptytable(self):
        """Parsing 4uvo with empty results table."""
        txt_file = get_file("2uvo_hhblits_emptytable.hhr")
        qresults = parse(txt_file, FMT)

        with pytest.raises(RuntimeError):
            next(qresults)

    def test_allx(self):
        """Parsing allx.hhr file."""
        txt_file = get_file("allx.hhr")
        qresults = parse(txt_file, FMT)

        # test first and only qresult
        qresult = next(qresults)

        num_hits = 10
        assert "HHSUITE" == qresult.program
        assert "Only X amino acids" == qresult.id
        assert 39 == qresult.seq_len
        assert num_hits == len(qresult)

        hit = qresult[0]
        assert "1klr_A" == hit.id
        assert ("Zinc finger Y-chromosomal protein; transcription; NMR {Synthetic} SCOP: g.37.1.1 PDB: "
            "5znf_A 1kls_A 1xrz_A* 7znf_A" == hit.description)
        assert hit.is_included
        assert 3.4e04 == hit.evalue
        assert -0.01 == hit.score
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert hsp.is_included
        assert 0 == hsp.output_index
        assert 3.4e04 == hsp.evalue
        assert -0.01 == hsp.score
        assert 0.04 == hsp.prob
        assert 23 == hsp.hit_start
        assert 24 == hsp.hit_end
        assert 38 == hsp.query_start
        assert 39 == hsp.query_end
        assert "T" == hsp.hit.seq
        assert "X" == hsp.query.seq

        # Check last hit
        hit = qresult[num_hits - 1]
        assert "1zfd_A" == hit.id
        assert ("SWI5; DNA binding motif, zinc finger DNA binding domain; NMR {Saccharomyces cerevisiae}"
            " SCOP: g.37.1.1" == hit.description)
        assert hit.is_included
        assert 3.6e04 == hit.evalue
        assert 0.03 == hit.score
        assert 1 == len(hit)

        # Check we can get the original last HSP from the file.
        num_hsps = num_hits
        assert num_hsps == len(qresult.hsps)
        hsp = qresult.hsps[-1]

        assert hsp.is_included
        assert num_hsps - 1 == hsp.output_index
        assert 3.6e04 == hsp.evalue
        assert 0.03 == hsp.score
        assert 0.03 == hsp.prob
        assert 0 == hsp.hit_start
        assert 1 == hsp.hit_end
        assert 3 == hsp.query_start
        assert 4 == hsp.query_end
        assert "D" == hsp.hit.seq
        assert "X" == hsp.query.seq

    def test_4y9h_nossm(self):
        """Parsing 4y9h_hhsearch_server_NOssm.hhr file."""
        txt_file = get_file("4y9h_hhsearch_server_NOssm.hhr")
        qresults = parse(txt_file, FMT)

        # test first and only qresult
        qresult = next(qresults)

        num_hits = 29
        assert "HHSUITE" == qresult.program
        assert "4Y9H:A|PDBID|CHAIN|SEQUENCE" == qresult.id
        assert 226 == qresult.seq_len
        assert num_hits == len(qresult)

        hit = qresult[0]
        assert "5ZIM_A" == hit.id
        assert ("Bacteriorhodopsin; proton pump, membrane protein, PROTON; HET: L2P, RET; 1.25A {Halobacterium"
            " salinarum}; Related PDB entries: 1R84_A 1KG8_A 1KME_B 1KGB_A 1KG9_A 1KME_A 4X31_A 5ZIL_A 1E0P_A "
            "4X32_A 5ZIN_A 1S53_B 1S51_B 1S53_A 1S54_A 1F50_A 1S54_B 1S51_A 1F4Z_A 5J7A_A 1S52_B 1S52_A 4Y9H_A "
            "3T45_A 3T45_C 3T45_B 1C3W_A 1L0M_A" == hit.description)
        assert hit.is_included
        assert 2.1e-48 == hit.evalue
        assert 320.44 == hit.score
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert hsp.is_included
        assert 0 == hsp.output_index
        assert 2.1e-48 == hsp.evalue
        assert 320.44 == hsp.score
        assert 100.00 == hsp.prob
        assert 1 == hsp.hit_start
        assert 227 == hsp.hit_end
        assert 0 == hsp.query_start
        assert 226 == hsp.query_end
        assert ("GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYAD"
            "WLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRP"
            "EVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG" == hsp.hit.seq)
        assert ("GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYAD"
            "WLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRP"
            "EVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG" == hsp.query.seq)

        # Check last hit
        hit = qresult[num_hits - 1]
        assert "5ABB_Z" == hit.id
        assert ("PROTEIN TRANSLOCASE SUBUNIT SECY, PROTEIN; TRANSLATION, RIBOSOME, MEMBRANE PROTEIN, "
            "TRANSLOCON; 8.0A {ESCHERICHIA COLI}" == hit.description)
        assert hit.is_included
        assert 3.3e-05 == hit.evalue
        assert 51.24 == hit.score
        assert 1 == len(hit)

        # Check we can get the original last HSP from the file.
        num_hsps = num_hits
        assert num_hsps == len(qresult.hsps)
        hsp = qresult.hsps[-1]

        assert hsp.is_included
        assert num_hsps - 1 == hsp.output_index
        assert 3.3e-05 == hsp.evalue
        assert 51.24 == hsp.score
        assert 96.55 == hsp.prob
        assert 14 == hsp.hit_start
        assert 65 == hsp.hit_end
        assert 7 == hsp.query_start
        assert 59 == hsp.query_end
        assert "FWLVTAALLASTVFFFVERDRVS-AKWKTSLTVSGLVTGIAFWHYMYMRGVW" == hsp.hit.seq
        assert "LALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGY" == hsp.query.seq

    def test_q9bsu1(self):
        """Parsing hhsearch_q9bsu1_uniclust_w_ss_pfamA_30.hhr file."""
        txt_file = get_file("hhsearch_q9bsu1_uniclust_w_ss_pfamA_30.hhr")
        qresults = parse(txt_file, FMT)

        # test first and only qresult
        qresult = next(qresults)

        num_hits = 12
        assert "HHSUITE" == qresult.program
        assert ("sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70"
            " PE=1 SV=1" == qresult.id)
        assert 422 == qresult.seq_len
        assert num_hits == len(qresult)

        hit = qresult[0]
        assert "PF03676.13" == hit.id
        assert "UPF0183 ; Uncharacterised protein family (UPF0183)" == hit.description
        assert hit.is_included
        assert 2e-106 == hit.evalue
        assert 822.75 == hit.score
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert hsp.is_included
        assert 0 == hsp.output_index
        assert 2e-106 == hsp.evalue
        assert 822.75 == hsp.score
        assert 100.00 == hsp.prob
        assert 0 == hsp.hit_start
        assert 395 == hsp.hit_end
        assert 10 == hsp.query_start
        assert 407 == hsp.query_end
        assert ("SLGNEQWEFTLGMPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCDLTKVKLK"
            "YCGVHFNSQAIAPTIEQIDQSFGATHPGVYNSAEQLFHLNFRGLSFSFQLDSWTEAPKYEPNFAHGLASLQIPHGATVKR"
            "MYIYSGNSLQDTKAPMMPLSCFLGNVYAESVDVLRDGTGPAGLRLRLLAAGCGPGLLADAKMRVFERSVYFGDSCQDVLS"
            "MLGSPHKVFYKSEDKMKIHSPSPHKQVPSKCNDYFFNYFTLGVDILFDANTHKVKKFVLHTNYPGHYNFNIYHRCEFKIP"
            "LAIKKENADGQTE--TCTTYSKWDNIQELLGHPVEKPVVLHRSSSPNNTNPFGSTFCFGLQRMIFEVMQNNHIASVTLY" == hsp.query.seq)
        assert ("EQWE----FALGMPLAQAISILQKHCRIIKNVQVLYSEQMPLSHDLILNLTQDGIKLLFDACNQRLKVIEVYDLTKVKLK"
            "YCGVHFNSQAIAPTIEQIDQSFGATHPGVYNAAEQLFHLNFRGLSFSFQLDSWSEAPKYEPNFAHGLASLQIPHGATVKR"
            "MYIYSGNNLQETKAPAMPLACFLGNVYAECVEVLRDGAGPLGLKLRLLTAGCGPGVLADTKVRAVERSIYFGDSCQDVLS"
            "ALGSPHKVFYKSEDKMKIHSPSPHKQVPSKCNDYFFNYYILGVDILFDSTTHLVKKFVLHTNFPGHYNFNIYHRCDFKIP"
            "LIIKKDGADAHSEDCILTTYSKWDQIQELLGHPMEKPVVLHRSSSANNTNPFGSTFCFGLQRMIFEVMQNNHIASVTLY" == hsp.hit.seq)

        # Check last hit
        hit = qresult[num_hits - 1]
        assert "PF10049.8" == hit.id
        assert "DUF2283 ; Protein of unknown function (DUF2283)" == hit.description
        assert hit.is_included
        assert 78 == hit.evalue
        assert 19.81 == hit.score
        assert 1 == len(hit)

        # Check we can get the original last HSP from the file.
        num_hsps = 16
        assert num_hsps == len(qresult.hsps)

        hsp = qresult.hsps[-1]
        assert hsp.is_included
        assert num_hsps - 1 == hsp.output_index
        assert 78 == hsp.evalue
        assert 19.81 == hsp.score
        assert 20.88 == hsp.prob
        assert 25 == hsp.hit_start
        assert 48 == hsp.hit_end
        assert 61 == hsp.query_start
        assert 85 == hsp.query_end
        assert "APNVIFDYDA-EGRIVGIELLDAR" == hsp.hit.seq
        assert "QDGIKLMFDAFNQRLKVIEVCDLT" == hsp.query.seq

    def test_4p79(self):
        """Parsing 4p79_hhsearch_server_NOssm.hhr file."""
        txt_file = get_file("4p79_hhsearch_server_NOssm.hhr")
        qresults = parse(txt_file, FMT)

        # test first and only qresult
        qresult = next(qresults)

        num_hits = 8
        assert "HHSUITE" == qresult.program
        assert "4P79:A|PDBID|CHAIN|SEQUENCE" == qresult.id
        assert 198 == qresult.seq_len
        assert num_hits == len(qresult)

        hit = qresult[0]
        assert "4P79_A" == hit.id
        assert ("cell adhesion protein; cell adhesion, tight junction, membrane; HET: OLC"
            ", MSE; 2.4A {Mus musculus}" == hit.description)
        assert hit.is_included
        assert 6.8e-32 == hit.evalue
        assert 194.63 == hit.score
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert hsp.is_included
        assert 0 == hsp.output_index
        assert 6.8e-32 == hsp.evalue
        assert 194.63 == hsp.score
        assert 99.94 == hsp.prob
        assert 0 == hsp.hit_start
        assert 198 == hsp.hit_end
        assert 0 == hsp.query_start
        assert 198 == hsp.query_end
        assert ("GSEFMSVAVETFGFFMSALGLLMLGLTLSNSYWRVSTVHGNVITTNTIFENLWYSCATDSLGVSNCWDFPSMLALSGYVQ"
            "GCRALMITAILLGFLGLFLGMVGLRATNVGNMDLSKKAKLLAIAGTLHILAGACGMVAISWYAVNITTDFFNPLYAGTKY"
            "ELGPALYLGWSASLLSILGGICVFSTAAASSKEEPATR" == hsp.query.seq)
        assert ("GSEFMSVAVETFGFFMSALGLLMLGLTLSNSYWRVSTVHGNVITTNTIFENLWYSCATDSLGVSNCWDFPSMLALSGYVQ"
            "GCRALMITAILLGFLGLFLGMVGLRATNVGNMDLSKKAKLLAIAGTLHILAGACGMVAISWYAVNITTDFFNPLYAGTKY"
            "ELGPALYLGWSASLLSILGGICVFSTAAASSKEEPATR" == hsp.hit.seq)

        # Check last hit
        hit = qresult[num_hits - 1]
        assert "5YQ7_F" == hit.id
        assert ("Beta subunit of light-harvesting 1; Photosynthetic core complex, PHOTOSYNTHESIS; "
            "HET: MQE, BCL, HEM, KGD, BPH;{Roseiflexus castenholzii}; Related PDB entries: 5YQ7_V"
            " 5YQ7_3 5YQ7_T 5YQ7_J 5YQ7_9 5YQ7_N 5YQ7_A 5YQ7_P 5YQ7_H 5YQ7_D 5YQ7_5 5YQ7_7 5YQ7_1 "
            "5YQ7_R" == hit.description)
        assert hit.is_included
        assert 6.7 == hit.evalue
        assert 20.51 == hit.score
        assert 1 == len(hit)

        # Check we can get the original last HSP from the file.
        num_hsps = num_hits
        assert num_hsps == len(qresult.hsps)

        hsp = qresult.hsps[-1]
        assert hsp.is_included
        assert num_hsps - 1 == hsp.output_index
        assert 6.7 == hsp.evalue
        assert 20.51 == hsp.score
        assert 52.07 == hsp.prob
        assert 8 == hsp.hit_start
        assert 42 == hsp.hit_end
        assert 5 == hsp.query_start
        assert 37 == hsp.query_end
        assert "RTSVVVSTLLGLVMALLIHFVVLSSGAFNWLRAP" == hsp.hit.seq
        assert "SVAVETFGFFMSALGLLMLGLTLSNS--YWRVST" == hsp.query.seq

    def test_9590198(self):
        """Parsing hhpred_9590198.hhr file."""
        txt_file = get_file("hhpred_9590198.hhr")
        qresults = parse(txt_file, FMT)

        # test first and only qresult
        qresult = next(qresults)

        num_hits = 22
        assert "HHSUITE" == qresult.program
        assert ("sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70"
            " PE=1 SV=1" == qresult.id)
        assert 422 == qresult.seq_len
        assert num_hits == len(qresult)

        hit = qresult[0]
        assert "PF03676.14" == hit.id
        assert "UPF0183 ; Uncharacterised protein family (UPF0183)" == hit.description
        assert hit.is_included
        assert 9.9e-102 == hit.evalue
        assert 792.76 == hit.score
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert hsp.is_included
        assert 0 == hsp.output_index
        assert 9.9e-102 == hsp.evalue
        assert 792.76 == hsp.score
        assert 100.00 == hsp.prob
        assert 0 == hsp.hit_start
        assert 394 == hsp.hit_end
        assert 21 == hsp.query_start
        assert 407 == hsp.query_end
        assert ("GMHFSQSVAIIQSQVGTIRGVQVLYSDQNPLSVDLVINMPQDGMRLIFDPVAQRLKIIEIYNMKLVKLRYSGMCFNSPEI"
            "TPSIEQVEHCFGATHPGLYDSQRHLFALNFRGLSFYFPVDS-----KFEPGYAHGLGSLQFPNGGSPVVSRTTIYYGSQH"
            "QLSSNTSSRVSGVPLPDLPLSCYRQQLHLRRCDVLRNTTSTMGLRLHMFTEGT--SRALEPSQVALVRVVRFGDSCQGVA"
            "RALGAPARLYYKADDKMRIHRPTARRR-PPPASDYLFNYFTLGLDVLFDARTNQVKKFVLHTNYPGHYNFNMYHRCEFEL"
            "TVQPD-KSEAHSLVESGGGVAVTAYSKWEVVSRAL-RVCERPVVLNRASSTNTTNPFGSTFCYGYQDIIFEVMSNNYIAS"
            "ITLY" == hsp.hit.seq)
        assert ("GMPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCDLTKVKLKYCGVHFNSQAI"
            "APTIEQIDQSFGATHPGVYNSAEQLFHLNFRGLSFSFQLDSWTEAPKYEPNFAHGLASLQIPHGA--TVKRMYIYSGNSL"
            "Q---------DTKA-PMMPLSCFLGNVYAESVDVLRDGTGPAGLRLRLLAAGCGPGLLADAKMRVFERSVYFGDSCQDVL"
            "SMLGSPHKVFYKSEDKMKIHSPSPHKQVPSKCNDYFFNYFTLGVDILFDANTHKVKKFVLHTNYPGHYNFNIYHRCEFKI"
            "PLAIKKENADG------QTETCTTYSKWDNIQELLGHPVEKPVVLHRSSSPNNTNPFGSTFCFGLQRMIFEVMQNNHIAS"
            "VTLY" == hsp.query.seq)

        # Check last hit
        hit = qresult[num_hits - 1]
        assert "4IL7_A" == hit.id
        assert ("Putative uncharacterized protein; partial jelly roll fold, hypothetical; 1.4A "
            "{Sulfolobus turreted icosahedral virus}" == hit.description)
        assert hit.is_included
        assert 6.8e02 == hit.evalue
        assert 22.72 == hit.score
        assert 1 == len(hit)

        # Check we can get the original last HSP from the file.
        num_hsps = 34
        assert num_hsps == len(qresult.hsps)

        hsp = qresult.hsps[-1]
        assert hsp.is_included
        assert num_hsps - 1 == hsp.output_index
        assert 3.9e02 == hsp.evalue
        assert 22.84 == hsp.score
        assert 21.56 == hsp.prob
        assert 7 == hsp.hit_start
        assert 96 == hsp.hit_end
        assert 18 == hsp.query_start
        assert 114 == hsp.query_end
        assert ("FTLGMPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCDLTKVKLKYCGVH-FN"
            "SQAIAPTIEQIDQSFGA" == hsp.query.seq)
        assert ("IQFGMDRTLVWQLAGADQSCSDQVERIICYNNPDH-------YGPQGHFFFNA-ADKLIHKRQMELFPAPKPTMRLATYN"
            "KTQTGMTEAQFWAAVPS" == hsp.hit.seq)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
