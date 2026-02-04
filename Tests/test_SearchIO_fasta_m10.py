# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO FastaIO parsers."""

import os
import unittest
import pytest

from Bio.SearchIO import parse

# test case files are in the Blast directory
TEST_DIR = "Fasta"
FMT = "fasta-m10"


def get_file(filename):
    """Return the path of a test file."""
    return os.path.join(TEST_DIR, filename)


class Fasta34Cases(unittest.TestCase):
    def test_output002(self):
        """Test parsing fasta34 output (output002.m10)."""
        m10_file = get_file("output002.m10")
        qresults = list(parse(m10_file, FMT))
        assert 3 == len(qresults)
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                assert qresult.id == hit.query_id
                for hsp in hit:
                    assert hit.id == hsp.hit_id
                    assert qresult.id == hsp.query_id

        # test first qresult
        qresult = qresults[0]
        assert "gi|10955263|ref|NP_052604.1|" == qresult.id
        assert "fasta" == qresult.program
        assert "34.26" == qresult.version
        assert "NC_002695.faa" == qresult.target
        assert 107 == qresult.seq_len
        assert "plasmid mobilization [Escherichia coli O157:H7 s 107 aa" == qresult.description
        assert 2 == len(qresult)
        # first qresult, first hit
        hit = qresult[0]
        assert "gi|162139799|ref|NP_309634.2|" == hit.id
        assert "23S rRNA pseudouridine synthase E [Escherichia coli O157:H7 str. Sakai]" == hit.description
        assert 207 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 55 == hsp.initn_score
        assert 55 == hsp.init1_score
        assert 77 == hsp.opt_score
        assert 110.8 == hsp.z_score
        assert 26.5 == hsp.bitscore
        assert 1.2 == hsp.evalue
        assert 77 == hsp.sw_score
        assert 28.4 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 54.5 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 88 == hsp.aln_span
        assert 4 == hsp.query_start
        assert 89 == hsp.query_end
        assert "SGSNTRRRAISRPVR--LTAEEDQEIRKRAAECG-KTVSGFLRAAALGKKVNSLTDDRVLKEVMRLGALQKKLFIDGKRVGDREYAEV" == hsp.query.seq
        assert 15 == hsp.hit_start
        assert 103 == hsp.hit_end
        assert "SQRSTRRKPENQPTRVILFNKPYDVLPQFTDEAGRKTLKEFIPVQGVYAAGRLDRDSEGLLVLTNNGALQARLTQPGKRTGKIYYVQV" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {} == hsp.aln_annotation
        # first qresult, second hit
        hit = qresult[1]
        assert "gi|15831859|ref|NP_310632.1|" == hit.id
        assert "trehalose-6-phosphate phosphatase [Escherichia coli O157:H7 str. Sakai]" == hit.description
        assert 266 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 43 == hsp.initn_score
        assert 43 == hsp.init1_score
        assert 69 == hsp.opt_score
        assert 98.6 == hsp.z_score
        assert 24.6 == hsp.bitscore
        assert 5.8 == hsp.evalue
        assert 69 == hsp.sw_score
        assert 28.3 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 66.0 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 53 == hsp.aln_span
        assert 26 == hsp.query_start
        assert 74 == hsp.query_end
        assert "EIRKRAAECGKTVSGFLRAAA-LGKKV----NSLTDDRVLKEVMRLGALQKKL" == hsp.query.seq
        assert 166 == hsp.hit_start
        assert 219 == hsp.hit_end
        assert "EIKPRGTSKGEAIAAFMQEAPFIGRTPVFLGDDLTDESGFAVVNRLGGMSVKI" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {} == hsp.aln_annotation

        # test second qresult
        qresult = qresults[1]
        assert "gi|10955264|ref|NP_052605.1|" == qresult.id
        assert "fasta" == qresult.program
        assert "34.26" == qresult.version
        assert "NC_002695.faa" == qresult.target
        assert 126 == qresult.seq_len
        assert "hypothetical protein pOSAK1_02 [Escherichia coli O157:H7 s 126 aa" == qresult.description
        assert 2 == len(qresult)
        # second qresult, first hit
        hit = qresult[0]
        assert "gi|15829419|ref|NP_308192.1|" == hit.id
        assert "serine endoprotease [Escherichia coli O157:H7 str. Sakai]" == hit.description
        assert 474 == hit.seq_len
        assert 1 == len(hit)
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 64 == hsp.initn_score
        assert 40 == hsp.init1_score
        assert 77 == hsp.opt_score
        assert 105.8 == hsp.z_score
        assert 27.0 == hsp.bitscore
        assert 2.3 == hsp.evalue
        assert 77 == hsp.sw_score
        assert 25.0 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 62.0 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 92 == hsp.aln_span
        assert 30 == hsp.query_start
        assert 117 == hsp.query_end
        assert "SEFFSKIESDLKKKKSKGDVFFDLIIPNG-----GKKDRYVYTSFNGEKFSSYTLNKVTKTDEYNDLSELSASFFKKNFDKINVNLLSKATS" == hsp.query.seq
        assert 295 == hsp.hit_start
        assert 384 == hsp.hit_end
        assert "TELNSELAKAMKVDAQRG-AFVSQVLPNSSAAKAGIKAGDVITSLNGKPISSFAALRA-QVGTMPVGSKLTLGLLRDG-KQVNVNLELQQSS" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {} == hsp.aln_annotation
        # second qresult, second hit
        hit = qresult[1]
        assert "gi|15832592|ref|NP_311365.1|" == hit.id
        assert "phosphoribosylaminoimidazole-succinocarboxamide synthase [Escherichia coli O157:H7 str. Sakai]" == hit.description
        assert 237 == hit.seq_len
        assert 1 == len(hit)
        # second qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 73 == hsp.initn_score
        assert 45 == hsp.init1_score
        assert 74 == hsp.opt_score
        assert 105.5 == hsp.z_score
        assert 26.0 == hsp.bitscore
        assert 2.4 == hsp.evalue
        assert 74 == hsp.sw_score
        assert 27.4 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 58.9 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 73 == hsp.aln_span
        assert 50 == hsp.query_start
        assert 123 == hsp.query_end
        assert "FFDLIIPNGGKKDRYVYTSFNGEKFSSYTLNKVTKTDEYNDLSELSASFFKKNFDKINVNLLSKATSFALKKG" == hsp.query.seq
        assert 116 == hsp.hit_start
        assert 185 == hsp.hit_end
        assert "LFDLFLKNDAMHDPMVNESYC-ETFGWVSKENLARMKE---LTYKANDVLKKLFDDAGLILVDFKLEFGLYKG" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {} == hsp.aln_annotation

        # test third qresult
        qresult = qresults[2]
        assert "gi|10955265|ref|NP_052606.1|" == qresult.id
        assert "fasta" == qresult.program
        assert "34.26" == qresult.version
        assert "NC_002695.faa" == qresult.target
        assert 346 == qresult.seq_len
        assert "hypothetical protein pOSAK1_03 [Escherichia coli O157:H7 s 346 aa" == qresult.description
        assert 2 == len(qresult)
        # third qresult, first hit
        hit = qresult[0]
        assert "gi|38704138|ref|NP_311957.2|" == hit.id
        assert "hypothetical protein ECs3930 [Escherichia coli O157:H7 str. Sakai]" == hit.description
        assert 111 == hit.seq_len
        assert 1 == len(hit)
        # third qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 50 == hsp.initn_score
        assert 50 == hsp.init1_score
        assert 86 == hsp.opt_score
        assert 117.5 == hsp.z_score
        assert 28.6 == hsp.bitscore
        assert 0.51 == hsp.evalue
        assert 86 == hsp.sw_score
        assert 30.2 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 63.5 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 63 == hsp.aln_span
        assert 187 == hsp.query_start
        assert 246 == hsp.query_end
        assert "VDIKK-ETIESELHSKLPKSIDKIHEDIKKQLSCSLI--MKKID-VEMEDYSTYCFSALRAIE" == hsp.query.seq
        assert 13 == hsp.hit_start
        assert 76 == hsp.hit_end
        assert "IDPKKIEQIARQVHESMPKGIREFGEDVEKKIRQTLQAQLTRLDLVSREEFDVQTQVLLRTRE" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {} == hsp.aln_annotation
        # third qresult, second hit
        hit = qresult[1]
        assert "gi|15833861|ref|NP_312634.1|" == hit.id
        assert "hypothetical protein ECs4607 [Escherichia coli O157:H7 str. Sakai]" == hit.description
        assert 330 == hit.seq_len
        assert 1 == len(hit)
        # third qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 32 == hsp.initn_score
        assert 32 == hsp.init1_score
        assert 87 == hsp.opt_score
        assert 112.7 == hsp.z_score
        assert 29.2 == hsp.bitscore
        assert 0.95 == hsp.evalue
        assert 87 == hsp.sw_score
        assert 21.0 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 58.0 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 157 == hsp.aln_span
        assert 130 == hsp.query_start
        assert 281 == hsp.query_end
        assert "QYIMTTSNGDRVRAKIYKRGSIQFQGKYLQIASLINDFMCSILNMKEIVEQKNKEFNVDI---KKETI-ESELHSKLPKSIDKIHEDIKKQLSCSLIMKKIDV-EMEDYSTYCFSALRA-IEGFIYQILNDVCNPSSSKNLGEYFTENKPKYIIREI" == hsp.query.seq
        assert 9 == hsp.hit_start
        assert 155 == hsp.hit_end
        assert "EFIRLLSDHDQFEKDQISELTVAANALKLEVAK--NNY-----NMKYSFDTQTERRMIELIREQKDLIPEKYLHQSGIKKL-KLHED---EFSSLLVDAERQVLEGSSFVLCCGEKINSTISELLSKKITDLTHPTESFTLSEYFSYDVYEEIFKKV" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {} == hsp.aln_annotation

    def test_output003(self):
        """Test parsing fasta34 output (output003.m10)."""
        m10_file = get_file("output003.m10")
        qresults = list(parse(m10_file, FMT))
        assert 5 == len(qresults)
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                assert qresult.id == hit.query_id
                for hsp in hit:
                    assert hit.id == hsp.hit_id
                    assert qresult.id == hsp.query_id

        # test first qresult
        qresult = qresults[0]
        assert "gi|152973837|ref|YP_001338874.1|" == qresult.id
        assert "fasta" == qresult.program
        assert "34.26" == qresult.version
        assert "NC_002127.faa" == qresult.target
        assert 183 == qresult.seq_len
        assert "hypothetical protein KPN_pKPN7p10262 [Klebsiella pneumoniae subsp. pneumonia 183 aa" == qresult.description
        assert 1 == len(qresult)
        # first qresult, first hit
        hit = qresult[0]
        assert "gi|10955263|ref|NP_052604.1|" == hit.id
        assert "plasmid mobilization [Escherichia coli O157:H7 str. Sakai]" == hit.description
        assert 107 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 43 == hsp.initn_score
        assert 43 == hsp.init1_score
        assert 45 == hsp.opt_score
        assert 64.1 == hsp.z_score
        assert 17.7 == hsp.bitscore
        assert 0.26 == hsp.evalue
        assert 59 == hsp.sw_score
        assert 25.5 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 67.3 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 55 == hsp.aln_span
        assert 86 == hsp.query_start
        assert 141 == hsp.query_end
        assert "ISISNNKDQYEELQKEQGERDLKTVDQLVRIAAAGGGLRLSASTKTVDQLVRIAA" == hsp.query.seq
        assert 17 == hsp.hit_start
        assert 69 == hsp.hit_end
        assert "VRLTAEEDQ--EIRKRAAECG-KTVSGFLRAAALGKKVNSLTDDRVLKEVMRLGA" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {} == hsp.aln_annotation

        # test second qresult
        qresult = qresults[1]
        assert "gi|152973838|ref|YP_001338875.1|" == qresult.id
        assert "fasta" == qresult.program
        assert "34.26" == qresult.version
        assert "NC_002127.faa" == qresult.target
        assert 76 == qresult.seq_len
        assert "hypothetical protein KPN_pKPN7p10263 [Klebsiella pneumoniae subsp. pneumonia 76 aa" == qresult.description
        assert 0 == len(qresult)

        # test third qresult
        qresult = qresults[2]
        assert "gi|152973839|ref|YP_001338876.1|" == qresult.id
        assert "fasta" == qresult.program
        assert "34.26" == qresult.version
        assert "NC_002127.faa" == qresult.target
        assert 112 == qresult.seq_len
        assert "hypothetical protein KPN_pKPN7p10264 [Klebsiella pneumoniae subsp. pneumonia 112 aa" == qresult.description
        assert 0 == len(qresult)

        # test fourth qresult
        qresult = qresults[3]
        assert "gi|152973840|ref|YP_001338877.1|" == qresult.id
        assert "fasta" == qresult.program
        assert "34.26" == qresult.version
        assert "NC_002127.faa" == qresult.target
        assert 63 == qresult.seq_len
        assert "RNA one modulator-like protein [Klebsiella pneumoniae subsp. pneumoniae  63 aa" == qresult.description
        assert 1 == len(qresult)
        # fourth qresult, first hit
        hit = qresult[0]
        assert "gi|10955265|ref|NP_052606.1|" == hit.id
        assert "hypothetical protein pOSAK1_03 [Escherichia coli O157:H7 str. Sakai]" == hit.description
        assert 346 == hit.seq_len
        assert 1 == len(hit)
        # fourth qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 35 == hsp.initn_score
        assert 35 == hsp.init1_score
        assert 38 == hsp.opt_score
        assert 71.3 == hsp.z_score
        assert 19.2 == hsp.bitscore
        assert 0.11 == hsp.evalue
        assert 38 == hsp.sw_score
        assert 36.4 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 63.6 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 22 == hsp.aln_span
        assert 42 == hsp.query_start
        assert 63 == hsp.query_end
        assert "DDAEHLFRTLSSR-LDALQDGN" == hsp.query.seq
        assert 101 == hsp.hit_start
        assert 123 == hsp.hit_end
        assert "DDRANLFEFLSEEGITITEDNN" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {} == hsp.aln_annotation

        # test fifth qresult
        qresult = qresults[4]
        assert "gi|152973841|ref|YP_001338878.1|" == qresult.id
        assert "fasta" == qresult.program
        assert "34.26" == qresult.version
        assert "NC_002127.faa" == qresult.target
        assert 133 == qresult.seq_len
        assert "Excl1 protein [Klebsiella pneumoniae subsp. pneumoniae  133 aa" == qresult.description
        assert 1 == len(qresult)
        # fifth qresult, first hit
        hit = qresult[0]
        assert "gi|10955264|ref|NP_052605.1|" == hit.id
        assert "hypothetical protein pOSAK1_02 [Escherichia coli O157:H7 str. Sakai]" == hit.description
        assert 126 == hit.seq_len
        assert 1 == len(hit)
        # fifth qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 37 == hsp.initn_score
        assert 37 == hsp.init1_score
        assert 57 == hsp.opt_score
        assert 80.0 == hsp.z_score
        assert 20.4 == hsp.bitscore
        assert 0.036 == hsp.evalue
        assert 57 == hsp.sw_score
        assert 25.4 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 65.1 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 63 == hsp.aln_span
        assert 48 == hsp.query_start
        assert 109 == hsp.query_end
        assert "VFGSFEQPKGEHLSGQVSEQ--RDTAFADQNEQVIRHLKQEIEHLNTLLLSKDSHIDSLKQAM" == hsp.query.seq
        assert 65 == hsp.hit_start
        assert 124 == hsp.hit_end
        assert "VYTSFN---GEKFSSYTLNKVTKTDEYNDLSELSASFFKKNFDKINVNLLSKATSF-ALKKGI" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {} == hsp.aln_annotation


class Fasta35Cases(unittest.TestCase):
    def test_output001(self):
        """Test parsing fasta35 output (output001.m10)."""
        m10_file = get_file("output001.m10")
        qresults = list(parse(m10_file, FMT))
        assert 3 == len(qresults)
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                assert qresult.id == hit.query_id
                for hsp in hit:
                    assert hit.id == hsp.hit_id
                    assert qresult.id == hsp.query_id

        # test first qresult
        qresult = qresults[0]
        assert "gi|10955263|ref|NP_052604.1|" == qresult.id
        assert "fasta" == qresult.program
        assert "35.03" == qresult.version
        assert "NC_009649.faa" == qresult.target
        assert 107 == qresult.seq_len
        assert "plasmid mobilization [Escherichia coli O157:H7 s 107 aa" == qresult.description
        assert 2 == len(qresult)
        # first qresult, first hit
        hit = qresult[0]
        assert "gi|152973457|ref|YP_001338508.1|" == hit.id
        assert "ATPase with chaperone activity, ATP-binding subunit [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]" == hit.description
        assert 931 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 65 == hsp.initn_score
        assert 43 == hsp.init1_score
        assert 71 == hsp.opt_score
        assert 92.7 == hsp.z_score
        assert 25.3 == hsp.bitscore
        assert 0.43 == hsp.evalue
        assert 71 == hsp.sw_score
        assert 25.0 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 57.4 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 108 == hsp.aln_span
        assert 4 == hsp.query_start
        assert 103 == hsp.query_end
        assert "SGSNT-RRRAISRPVRLTAEED---QEIRKRAAECGKTVSGFLRAAALGKKVNSLTDDRVLKEVM-----RLGALQKKLFIDGKRVGDREYAEVLIAITEYHRALLSR" == hsp.query.seq
        assert 95 == hsp.hit_start
        assert 195 == hsp.hit_end
        assert "AGSGAPRRRGSGLASRISEQSEALLQEAAKHAAEFGRS------EVDTEHLLLALADSDVVKTILGQFKIKVDDLKRQIESEAKR-GDKPF-EGEIGVSPRVKDALSR" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {} == hsp.aln_annotation
        # first qresult, second hit
        hit = qresult[1]
        assert "gi|152973588|ref|YP_001338639.1|" == hit.id
        assert "F pilus assembly protein [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]" == hit.description
        assert 459 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 33 == hsp.initn_score
        assert 33 == hsp.init1_score
        assert 63 == hsp.opt_score
        assert 87.7 == hsp.z_score
        assert 23.4 == hsp.bitscore
        assert 0.81 == hsp.evalue
        assert 63 == hsp.sw_score
        assert 26.6 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 65.6 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 64 == hsp.aln_span
        assert 31 == hsp.query_start
        assert 94 == hsp.query_end
        assert "AAECGKTVSGFLRAAALGKKVNSLTDDRVLKEV-MRLGALQKKLFIDGKRVGDREYAEVLIAIT" == hsp.query.seq
        assert 190 == hsp.hit_start
        assert 248 == hsp.hit_end
        assert "ASRQGCTVGG--KMDSVQDKASDKDKERVMKNINIMWNALSKNRLFDG----NKELKEFIMTLT" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {} == hsp.aln_annotation

        # test second qresult
        qresult = qresults[1]
        assert "gi|10955264|ref|NP_052605.1|" == qresult.id
        assert "fasta" == qresult.program
        assert "35.03" == qresult.version
        assert "NC_009649.faa" == qresult.target
        assert 126 == qresult.seq_len
        assert "hypothetical protein pOSAK1_02 [Escherichia coli O157:H7 s 126 aa" == qresult.description
        assert 1 == len(qresult)
        # second qresult, first hit
        hit = qresult[0]
        assert "gi|152973462|ref|YP_001338513.1|" == hit.id
        assert "hypothetical protein KPN_pKPN3p05904 [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]" == hit.description
        assert 101 == hit.seq_len
        assert 1 == len(hit)
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 50 == hsp.initn_score
        assert 50 == hsp.init1_score
        assert 58 == hsp.opt_score
        assert 91.6 == hsp.z_score
        assert 22.2 == hsp.bitscore
        assert 0.49 == hsp.evalue
        assert 58 == hsp.sw_score
        assert 28.9 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 63.2 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 38 == hsp.aln_span
        assert 0 == hsp.query_start
        assert 38 == hsp.query_end
        assert "MKKDKKYQIEAIKNKDKTLFIVYATDIYSPSEFFSKIE" == hsp.query.seq
        assert 43 == hsp.hit_start
        assert 81 == hsp.hit_end
        assert "IKKDLGVSFLKLKNREKTLIVDALKKKYPVAELLSVLQ" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {} == hsp.aln_annotation

        # test third qresult
        qresult = qresults[2]
        assert "gi|10955265|ref|NP_052606.1|" == qresult.id
        assert "fasta" == qresult.program
        assert "35.03" == qresult.version
        assert "NC_009649.faa" == qresult.target
        assert 346 == qresult.seq_len
        assert "hypothetical protein pOSAK1_03 [Escherichia coli O157:H7 s 346 aa" == qresult.description
        assert 1 == len(qresult)
        # third qresult, first hit
        hit = qresult[0]
        assert "gi|152973545|ref|YP_001338596.1|" == hit.id
        assert "putative plasmid SOS inhibition protein A [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]" == hit.description
        assert 242 == hit.seq_len
        assert 1 == len(hit)
        # third qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 52 == hsp.initn_score
        assert 52 == hsp.init1_score
        assert 70 == hsp.opt_score
        assert 94.0 == hsp.z_score
        assert 25.3 == hsp.bitscore
        assert 0.36 == hsp.evalue
        assert 70 == hsp.sw_score
        assert 27.9 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 65.1 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 43 == hsp.aln_span
        assert 196 == hsp.query_start
        assert 238 == hsp.query_end
        assert "SELHSKLPKSIDKIHEDIKKQLSC-SLIMKKIDVEMEDYSTYC" == hsp.query.seq
        assert 51 == hsp.hit_start
        assert 94 == hsp.hit_end
        assert "SRINSDVARRIPGIHRDPKDRLSSLKQVEEALDMLISSHGEYC" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {} == hsp.aln_annotation

    def test_output004(self):
        """Test parsing fasta35 output (output004.m10)."""
        m10_file = get_file("output004.m10")
        qresults = list(parse(m10_file, FMT))
        assert 3 == len(qresults)
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                assert qresult.id == hit.query_id
                for hsp in hit:
                    assert hit.id == hsp.hit_id
                    assert qresult.id == hsp.query_id

        # test first qresult
        qresult = qresults[0]
        assert "ref|NC_002127.1|:413-736" == qresult.id
        assert "fasta" == qresult.program
        assert "35.04" == qresult.version
        assert "NC_002695.ffn" == qresult.target
        assert 324 == qresult.seq_len
        assert "" == qresult.description
        assert 0 == len(qresult)

        # test second qresult
        qresult = qresults[1]
        assert "ref|NC_002127.1|:c1351-971" == qresult.id
        assert "fasta" == qresult.program
        assert "35.04" == qresult.version
        assert "NC_002695.ffn" == qresult.target
        assert 381 == qresult.seq_len
        assert "" == qresult.description
        assert 1 == len(qresult)
        # second qresult, first hit
        hit = qresult[0]
        assert "ref|NC_002695.1|:1970775-1971404" == hit.id
        assert "" == hit.description
        assert 630 == hit.seq_len
        assert 1 == len(hit)
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 54 == hsp.initn_score
        assert 54 == hsp.init1_score
        assert 91 == hsp.opt_score
        assert 139.3 == hsp.z_score
        assert 35.2 == hsp.bitscore
        assert 0.045 == hsp.evalue
        assert 57.8 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 57.8 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 102 == hsp.aln_span
        assert 3 == hsp.query_start
        assert 105 == hsp.query_end
        assert "AAAAAAGATAAAAAATATCAAATAGAAGCAATAAAAAATAAAGATAAAACTTTATTTATTGTCTATGCTACTGATATTTATAGCCCGAGCGAATTTTTCTCA" == hsp.query.seq
        assert 312 == hsp.hit_start
        assert 414 == hsp.hit_end
        assert "AGAGAAAATAAAACAAGTAATAAAATATTAATGGAAAAAATAAATTCTTGTTTATTTAGACCTGATTCTAATCACTTTTCTTGCCCGGAGTCATTTTTGACA" == hsp.hit.seq
        assert 1 == hsp.query_strand
        assert {
                "similarity": ": : :: :::::: :  : : : :  :  :::  :::: : : ::     ::::::::      :: ::: : :  ::: : :::::     ::::::  ::"
            } == hsp.aln_annotation

        # test third qresult
        qresult = qresults[2]
        assert "ref|NC_002127.1|:c2388-1348" == qresult.id
        assert "fasta" == qresult.program
        assert "35.04" == qresult.version
        assert "NC_002695.ffn" == qresult.target
        assert 1041 == qresult.seq_len
        assert "" == qresult.description
        assert 0 == len(qresult)

    def test_output005(self):
        """Test parsing ssearch35 output (output005.m10)."""
        m10_file = get_file("output005.m10")
        qresults = list(parse(m10_file, FMT))
        assert 3 == len(qresults)
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                assert qresult.id == hit.query_id
                for hsp in hit:
                    assert hit.id == hsp.hit_id
                    assert qresult.id == hsp.query_id

        # test first qresult
        qresult = qresults[0]
        assert "gi|10955263|ref|NP_052604.1|" == qresult.id
        assert "ssearch" == qresult.program
        assert "35.03" == qresult.version
        assert "NC_002128.faa" == qresult.target
        assert 107 == qresult.seq_len
        assert "plasmid mobilization [Escherichia coli O157:H7 s 107 aa" == qresult.description
        assert 0 == len(qresult)

        # test second qresult
        qresult = qresults[1]
        assert "gi|10955264|ref|NP_052605.1|" == qresult.id
        assert "ssearch" == qresult.program
        assert "35.03" == qresult.version
        assert "NC_002128.faa" == qresult.target
        assert 126 == qresult.seq_len
        assert "hypothetical protein pOSAK1_02 [Escherichia coli O157:H7 s 126 aa" == qresult.description
        assert 1 == len(qresult)
        # second qresult, first hit
        hit = qresult[0]
        assert "gi|10955282|ref|NP_052623.1|" == hit.id
        assert "hemolysin C [Escherichia coli O157:H7 str. Sakai]" == hit.description
        assert 163 == hit.seq_len
        assert 1 == len(hit)
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 69 == hsp.opt_score
        assert 108.8 == hsp.z_score
        assert 26.0 == hsp.bitscore
        assert 0.025 == hsp.evalue
        assert 69 == hsp.sw_score
        assert 20.9 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 55.5 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 110 == hsp.aln_span
        assert 11 == hsp.query_start
        assert 114 == hsp.query_end
        assert "IKNKDKTLFIVYAT-DIYSPSEFFSKIESDLKKKKSKGDV--FFDLIIPNGGKKD--RYVYTSFNGEKFSSYTLNKVTKTDEYNDL--SELSASFFKKNFDKINVNLLSK" == hsp.query.seq
        assert 38 == hsp.hit_start
        assert 148 == hsp.hit_end
        assert "IKDELPVAFCSWASLDLECEVKYINDVTSLYAKDWMSGERKWFIDWIAPFGHNMELYKYMRKKYPYELFRAIRLDESSKTGKIAEFHGGGIDKKLASKIFRQYHHELMSE" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {} == hsp.aln_annotation

        # test third qresult
        qresult = qresults[2]
        assert "gi|10955265|ref|NP_052606.1|" == qresult.id
        assert "ssearch" == qresult.program
        assert "35.03" == qresult.version
        assert "NC_002128.faa" == qresult.target
        assert 346 == qresult.seq_len
        assert "hypothetical protein pOSAK1_03 [Escherichia coli O157:H7 s 346 aa" == qresult.description
        assert 0 == len(qresult)

    def test_output006(self):
        """Test parsing fasta35 output (output006.m10)."""
        m10_file = get_file("output006.m10")
        qresults = list(parse(m10_file, FMT))
        assert 1 == len(qresults)
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                assert qresult.id == hit.query_id
                for hsp in hit:
                    assert hit.id == hsp.hit_id
                    assert qresult.id == hsp.query_id

        # test first qresult
        qresult = qresults[0]
        assert "query" == qresult.id
        assert "fasta" == qresult.program
        assert "35.04" == qresult.version
        assert "orchid_cds.txt" == qresult.target
        assert 131 == qresult.seq_len
        assert "" == qresult.description
        assert 1 == len(qresult)
        # first qresult, first hit
        hit = qresult[0]
        assert "gi|116660610|gb|EG558221.1|EG558221" == hit.id
        assert "CR03001A07 Root CR03 cDNA library Catharanthus roseus cDNA clone CR03001A07 5', mRNA sequence" == hit.description
        assert 573 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 646 == hsp.initn_score
        assert 646 == hsp.init1_score
        assert 646 == hsp.opt_score
        assert 712.3 == hsp.z_score
        assert 139.6 == hsp.bitscore
        assert 7.2e-38 == hsp.evalue
        assert 99.2 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 99.2 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 131 == hsp.aln_span
        assert 0 == hsp.query_start
        assert 131 == hsp.query_end
        assert "GCAACGCTTCAAGAACTGGAATTAGGAACCGTGACAACGATTAATGAGGAGATTTATGAAGAGGGTTCTTCGATTTTAGGCCAATCGGAAGGAATTATGTAGCAAGTCCATCAGAAAATGGAAGAAGTCAT" == hsp.query.seq
        assert 359 == hsp.hit_start
        assert 490 == hsp.hit_end
        assert "GCAACGCTTCAAGAACTGGAATTAGGAACCGTGACAACGATTAATGAGGAGATTTATGAAGAGGGTTCTTCGATTTTAGGCCAATCGGAAGGAATTATGTAGCAAGTCCATCAGAAAATGGAAGTAGTCAT" == hsp.hit.seq
        assert -1 == hsp.query_strand
        assert {
                "similarity": ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ::::::"
            } == hsp.aln_annotation


class Fasta36Cases(unittest.TestCase):
    def test_output007(self):
        """Test parsing fasta36 output (output007.m10)."""
        m10_file = get_file("output007.m10")
        qresults = list(parse(m10_file, FMT))
        assert 3 == len(qresults)
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                assert qresult.id == hit.query_id
                for hsp in hit:
                    assert hit.id == hsp.hit_id
                    assert qresult.id == hsp.query_id

        # test first qresult
        qresult = qresults[0]
        assert "gi|10955263|ref|NP_052604.1|" == qresult.id
        assert "fasta" == qresult.program
        assert "36.3.4" == qresult.version
        assert "NC_009649.faa" == qresult.target
        assert 107 == qresult.seq_len
        assert "plasmid mobilization [Escherichia coli O157:H7 s" == qresult.description
        assert 3 == len(qresult)
        # first qresult, first hit
        hit = qresult[0]
        assert "gi|152973457|ref|YP_001338508.1|" == hit.id
        assert "ATPase with chaperone activity, ATP-binding subunit [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]" == hit.description
        assert 931 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 97 == hsp.initn_score
        assert 43 == hsp.init1_score
        assert 71 == hsp.opt_score
        assert 109.6 == hsp.z_score
        assert 28.5 == hsp.bitscore
        assert 0.048 == hsp.evalue
        assert 71 == hsp.sw_score
        assert 25.0 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 57.4 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 108 == hsp.aln_span
        assert 4 == hsp.query_start
        assert 103 == hsp.query_end
        assert "SGSNT-RRRAISRPVRLTAEED---QEIRKRAAECGKTVSGFLRAAALGKKVNSLTDDRVLKEVM-----RLGALQKKLFIDGKRVGDREYAEVLIAITEYHRALLSR" == hsp.query.seq
        assert 95 == hsp.hit_start
        assert 195 == hsp.hit_end
        assert "AGSGAPRRRGSGLASRISEQSEALLQEAAKHAAEFGRS------EVDTEHLLLALADSDVVKTILGQFKIKVDDLKRQIESEAKR-GDKPF-EGEIGVSPRVKDALSR" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {
                "similarity": ".::..-:::. .   :.. . .---::  :.::: :..------ .   . . .:.:. :.: ..-----..  :....  ..::-::. .-:  :...   .  :::"
            } == hsp.aln_annotation
        # first qresult, second hit
        hit = qresult[1]
        assert "gi|152973588|ref|YP_001338639.1|" == hit.id
        assert "F pilus assembly protein [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]" == hit.description
        assert 459 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 66 == hsp.initn_score
        assert 33 == hsp.init1_score
        assert 63 == hsp.opt_score
        assert 101.4 == hsp.z_score
        assert 25.9 == hsp.bitscore
        assert 0.14 == hsp.evalue
        assert 63 == hsp.sw_score
        assert 26.6 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 65.6 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 64 == hsp.aln_span
        assert 31 == hsp.query_start
        assert 94 == hsp.query_end
        assert "AAECGKTVSGFLRAAALGKKVNSLTDDRVLKEV-MRLGALQKKLFIDGKRVGDREYAEVLIAIT" == hsp.query.seq
        assert 190 == hsp.hit_start
        assert 248 == hsp.hit_end
        assert "ASRQGCTVGG--KMDSVQDKASDKDKERVMKNINIMWNALSKNRLFDG----NKELKEFIMTLT" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {
                "similarity": ":.. : ::.:--.  ..  :...   .::.:..-.  .::.:. ..::----..:  : ....:"
            } == hsp.aln_annotation
        # first qresult, third hit
        hit = qresult[2]
        assert "gi|152973480|ref|YP_001338531.1|" == hit.id
        assert "Arsenate reductase (Arsenical pump modifier) [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]" == hit.description
        assert 141 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        assert 45 == hsp.initn_score
        assert 37 == hsp.init1_score
        assert 51 == hsp.opt_score
        assert 89.6 == hsp.z_score
        assert 22.0 == hsp.bitscore
        assert 0.63 == hsp.evalue
        assert 51 == hsp.sw_score
        assert 26.7 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 62.2 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 45 == hsp.aln_span
        assert 26 == hsp.query_start
        assert 66 == hsp.query_end
        assert "EIRKRAAECGKTVSGFLRAAA-----LGKKVNSLTDDRVLKEVMR" == hsp.query.seq
        assert 42 == hsp.hit_start
        assert 87 == hsp.hit_end
        assert "ELVKLIADMGISVRALLRKNVEPYEELGLEEDKFTDDQLIDFMLQ" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {"similarity": ":. :  :. : .: ..::  .-----:: . ...:::...  ..."} == hsp.aln_annotation

        # test second qresult
        qresult = qresults[1]
        assert "gi|10955264|ref|NP_052605.1|" == qresult.id
        assert "fasta" == qresult.program
        assert "36.3.4" == qresult.version
        assert "NC_009649.faa" == qresult.target
        assert 126 == qresult.seq_len
        assert "hypothetical protein pOSAK1_02 [Escherichia coli O157:H7 s" == qresult.description
        assert 4 == len(qresult)
        # second qresult, first hit
        hit = qresult[0]
        assert "gi|152973462|ref|YP_001338513.1|" == hit.id
        assert "hypothetical protein KPN_pKPN3p05904 [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]" == hit.description
        assert 101 == hit.seq_len
        assert 1 == len(hit)
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 78 == hsp.initn_score
        assert 50 == hsp.init1_score
        assert 58 == hsp.opt_score
        assert 100.8 == hsp.z_score
        assert 23.9 == hsp.bitscore
        assert 0.15 == hsp.evalue
        assert 58 == hsp.sw_score
        assert 28.9 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 63.2 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 38 == hsp.aln_span
        assert 0 == hsp.query_start
        assert 38 == hsp.query_end
        assert "MKKDKKYQIEAIKNKDKTLFIVYATDIYSPSEFFSKIE" == hsp.query.seq
        assert 43 == hsp.hit_start
        assert 81 == hsp.hit_end
        assert "IKKDLGVSFLKLKNREKTLIVDALKKKYPVAELLSVLQ" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {"similarity": ".:::   ..  .::..:::..      :  .:..: .."} == hsp.aln_annotation
        # second qresult, second hit
        hit = qresult[1]
        assert "gi|152973509|ref|YP_001338560.1|" == hit.id
        assert "probable sensor kinase (silver resistance) [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]" == hit.description
        assert 448 == hit.seq_len
        assert 1 == len(hit)
        # second qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 73 == hsp.initn_score
        assert 56 == hsp.init1_score
        assert 56 == hsp.opt_score
        assert 89.9 == hsp.z_score
        assert 24.0 == hsp.bitscore
        assert 0.6 == hsp.evalue
        assert 56 == hsp.sw_score
        assert 72.7 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 81.8 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 11 == hsp.aln_span
        assert 50 == hsp.query_start
        assert 61 == hsp.query_end
        assert "FFDLIIPNGGK" == hsp.query.seq
        assert 407 == hsp.hit_start
        assert 418 == hsp.hit_end
        assert "FFDLVIENPGK" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {"similarity": "::::.: : ::"} == hsp.aln_annotation
        # second qresult, third hit
        hit = qresult[2]
        assert "gi|152973581|ref|YP_001338632.1|" == hit.id
        assert "inner membrane protein [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]" == hit.description
        assert 84 == hit.seq_len
        assert 1 == len(hit)
        # second qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        assert 61 == hsp.initn_score
        assert 46 == hsp.init1_score
        assert 48 == hsp.opt_score
        assert 88.5 == hsp.z_score
        assert 21.3 == hsp.bitscore
        assert 0.72 == hsp.evalue
        assert 48 == hsp.sw_score
        assert 30.0 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 67.5 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 40 == hsp.aln_span
        assert 15 == hsp.query_start
        assert 53 == hsp.query_end
        assert "DKTLFIVYATDIYSPSE-FFSKIESDLKKKKSKGD-VFFD" == hsp.query.seq
        assert 44 == hsp.hit_start
        assert 84 == hsp.hit_end
        assert "ESVVFILMAGFAMSVCYLFFSVLEKVINARKSKDESIYHD" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {"similarity": "....::..:   .:   -::: .:. .. .::: .-.. :"} == hsp.aln_annotation
        # second qresult, fourth hit
        hit = qresult[3]
        assert "gi|152973536|ref|YP_001338587.1|" == hit.id
        assert "putative inner membrane protein [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]" == hit.description
        assert 84 == hit.seq_len
        assert 1 == len(hit)
        # second qresult, fourth hit, first hsp
        hsp = qresult[3].hsps[0]
        assert 63 == hsp.initn_score
        assert 42 == hsp.init1_score
        assert 48 == hsp.opt_score
        assert 88.5 == hsp.z_score
        assert 21.3 == hsp.bitscore
        assert 0.72 == hsp.evalue
        assert 48 == hsp.sw_score
        assert 26.7 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 66.7 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 30 == hsp.aln_span
        assert 96 == hsp.query_start
        assert 126 == hsp.query_end
        assert "ASFFKKNFDKINVNLLSKATSFALKKGIPI" == hsp.query.seq
        assert 6 == hsp.hit_start
        assert 36 == hsp.hit_end
        assert "ASFSKEEQDKVAVDKVAADVAWQERMNKPV" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {"similarity": "::: :.. ::. :. ..  ...  . . :."} == hsp.aln_annotation

        # test third qresult
        qresult = qresults[2]
        assert "gi|10955265|ref|NP_052606.1|" == qresult.id
        assert "fasta" == qresult.program
        assert "36.3.4" == qresult.version
        assert "NC_009649.faa" == qresult.target
        assert 346 == qresult.seq_len
        assert "hypothetical protein pOSAK1_03 [Escherichia coli O157:H7 s" == qresult.description
        assert 2 == len(qresult)
        # third qresult, first hit
        hit = qresult[0]
        assert "gi|152973545|ref|YP_001338596.1|" == hit.id
        assert "putative plasmid SOS inhibition protein A [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]" == hit.description
        assert 242 == hit.seq_len
        assert 1 == len(hit)
        # third qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 72 == hsp.initn_score
        assert 52 == hsp.init1_score
        assert 70 == hsp.opt_score
        assert 110.9 == hsp.z_score
        assert 28.4 == hsp.bitscore
        assert 0.041 == hsp.evalue
        assert 70 == hsp.sw_score
        assert 27.9 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 65.1 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 43 == hsp.aln_span
        assert 196 == hsp.query_start
        assert 238 == hsp.query_end
        assert "SELHSKLPKSIDKIHEDIKKQLSC-SLIMKKIDVEMEDYSTYC" == hsp.query.seq
        assert 51 == hsp.hit_start
        assert 94 == hsp.hit_end
        assert "SRINSDVARRIPGIHRDPKDRLSSLKQVEEALDMLISSHGEYC" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {"similarity": ":...: . . :  ::.: : .:: -. . . .:. . ... ::"} == hsp.aln_annotation
        # third qresult, second hit
        hit = qresult[1]
        assert "gi|152973505|ref|YP_001338556.1|" == hit.id
        assert "putative membrane fusion protein SilB [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]" == hit.description
        assert 430 == hit.seq_len
        assert 1 == len(hit)
        # third qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 95 == hsp.initn_score
        assert 52 == hsp.init1_score
        assert 57 == hsp.opt_score
        assert 90.1 == hsp.z_score
        assert 25.4 == hsp.bitscore
        assert 0.59 == hsp.evalue
        assert 57 == hsp.sw_score
        assert 23.4 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 60.9 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 64 == hsp.aln_span
        assert 39 == hsp.query_start
        assert 101 == hsp.query_end
        assert "ISGTYKGIDFLIKLMPSGGNTTIGRASGQNNTYFDEIALIIKENCLY--SDTKNFEYTIPKFSD" == hsp.query.seq
        assert 221 == hsp.hit_start
        assert 281 == hsp.hit_end
        assert "IDGVITAFD-LRTGMNISKDKVVAQIQGMDPVW---ISAAVPESIAYLLKDTSQFEISVPAYPD" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {
                "similarity": ":.:.  ..:-:   :  . . .... .:.. ..---:.  . :.  :--.::..:: ..: . :"
            } == hsp.aln_annotation

    def test_output008(self):
        """Test parsing tfastx36 output (output008.m10)."""
        m10_file = get_file("output008.m10")
        qresults = list(parse(m10_file, FMT))
        assert 4 == len(qresults)
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                assert qresult.id == hit.query_id
                for hsp in hit:
                    assert hit.id == hsp.hit_id
                    assert qresult.id == hsp.query_id

        # test first qresult
        qresult = qresults[0]
        assert "sp|Q9BS26|ERP44_HUMAN" == qresult.id
        assert "tfastx" == qresult.program
        assert "36.3.4" == qresult.version
        assert "rhodopsin_nucs.fasta" == qresult.target
        assert 406 == qresult.seq_len
        assert "Endoplasmic reticulum resident protein 44 OS=Homo sapiens GN=ERP44" == qresult.description
        assert 0 == len(qresult)

        # test second qresult
        qresult = qresults[1]
        assert "sp|Q9NSY1|BMP2K_HUMAN" == qresult.id
        assert "tfastx" == qresult.program
        assert "36.3.4" == qresult.version
        assert "rhodopsin_nucs.fasta" == qresult.target
        assert 1161 == qresult.seq_len
        assert "BMP-2-inducible protein kinase OS=Homo sapiens GN=BMP2K" == qresult.description
        assert 2 == len(qresult)
        # second qresult, first hit
        hit = qresult[0]
        assert "gi|283855822|gb|GQ290312.1|" == hit.id
        assert "Myotis ricketti voucher GQX10 rhodopsin (RHO) mRNA, partial cds" == hit.description
        assert 983 == hit.seq_len
        assert 1 == len(hit)
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 106 == hsp.initn_score
        assert 53 == hsp.init1_score
        assert 70 == hsp.opt_score
        assert 88.0 == hsp.z_score
        assert 28.0 == hsp.bitscore
        assert 0.026 == hsp.evalue
        assert 70 == hsp.sw_score
        assert 23.1 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 46.2 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 65 == hsp.aln_span
        assert 452 == hsp.query_start
        assert 514 == hsp.query_end
        assert "LQHRHPHQQQQQQQQQQQQQQQQQQQQQQQQQQQH---HHHHHHHLLQDAYMQQYQHATQQQQML" == hsp.query.seq
        assert 122 == hsp.hit_start
        assert 317 == hsp.hit_end
        assert "IPHQLPHALRHRPAQEAAHASQLHPAQPGCGQPLHGLWRLHHHPVYLYAWILRLRGHGMQSGGLL" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {
                "similarity": ". :. ::  ...  :.  . .: .  :    :  :---. :::   :    ..   :. :.  .:"
            } == hsp.aln_annotation
        # second qresult, second hit
        hit = qresult[1]
        assert "gi|57163782|ref|NM_001009242.1|" == hit.id
        assert "Felis catus rhodopsin (RHO), mRNA" == hit.description
        assert 1047 == hit.seq_len
        assert 1 == len(hit)
        # second qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 105 == hsp.initn_score
        assert 57 == hsp.init1_score
        assert 68 == hsp.opt_score
        assert 85.8 == hsp.z_score
        assert 27.7 == hsp.bitscore
        assert 0.034 == hsp.evalue
        assert 72 == hsp.sw_score
        assert 23.5 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 45.0 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 201 == hsp.aln_span
        assert 417 == hsp.query_start
        assert 599 == hsp.query_end
        assert "GPEIL---LGQ-GPPQQPPQQHRVLQQLQQGDWRLQQLH-------LQHRHPHQQQQQQQQQQQQQQQQQQQQQQQQQQQH-----HHHHHH-HLLQDAYMQQYQHATQQQQMLQQQF-LMHSVYQPQPSASQYPTMMPQYQQAFFQQQMLAQHQPSQQQASPEYLTSPQEFSPALVSYTSSLPA-QVGTIMDSSYSANRS" == hsp.query.seq
        assert 14 == hsp.hit_start
        assert 595 == hsp.hit_end
        assert "GPELLRALLQQNGCGTQPLRVPTVLPG*AMAVLHAGRLHVPAHRAWLPHQLPHALRHGPAQEAAHASQLHPAQPGRG*PLHGLRWLHHHPLH/PLCMDTLSLGPQDAIWRASLPHWAVKLPCGLWWSWPLSGTWWCVSP*ATSA------LGRTMP*WASLSPGSWHWPALHPPSLVGPGTSLKACSVHAGSTTTHSSQKS" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {
                "similarity": ":::.:---: :-:   :: .   ::    ..  .  .::-------: :. ::  ..   :.  . .: .  :  .    :-----:::  :- : .:.     : :  . .. .   -:  ...   : .. .  . :   .:------:.. .:   . ::     :    :.::.  .:: :-.: .   ...:...:"
            } == hsp.aln_annotation

        # test third qresult
        qresult = qresults[2]
        assert "sp|P06213|INSR_HUMAN" == qresult.id
        assert "tfastx" == qresult.program
        assert "36.3.4" == qresult.version
        assert "rhodopsin_nucs.fasta" == qresult.target
        assert 1382 == qresult.seq_len
        assert "Insulin receptor OS=Homo sapiens GN=INSR" == qresult.description
        assert 0 == len(qresult)

        # test fourth qresult
        qresult = qresults[3]
        assert "sp|P08100|OPSD_HUMAN" == qresult.id
        assert "tfastx" == qresult.program
        assert "36.3.4" == qresult.version
        assert "rhodopsin_nucs.fasta" == qresult.target
        assert 348 == qresult.seq_len
        assert "Rhodopsin OS=Homo sapiens GN=RHO" == qresult.description
        assert 6 == len(qresult)
        # fourth qresult, first hit
        hit = qresult[0]
        assert "gi|57163782|ref|NM_001009242.1|" == hit.id
        assert "Felis catus rhodopsin (RHO), mRNA" == hit.description
        assert 1047 == hit.seq_len
        assert 1 == len(hit)
        # fourth qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 2298 == hsp.initn_score
        assert 2298 == hsp.init1_score
        assert 2298 == hsp.opt_score
        assert 3150.5 == hsp.z_score
        assert 593.0 == hsp.bitscore
        assert 6.7e-173 == hsp.evalue
        assert 2298 == hsp.sw_score
        assert 96.6 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 99.4 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 348 == hsp.aln_span
        assert 0 == hsp.query_start
        assert 348 == hsp.query_end
        assert "MNGTEGPNFYVPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVLGGFTSTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSRYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAIYNPVIYIMMNKQFRNCMLTTICCGKNPLGDDEASATVSKTETSQVAPA" == hsp.query.seq
        assert 0 == hsp.hit_start
        assert 1044 == hsp.hit_end
        assert "MNGTEGPNFYVPFSNKTGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTLPAFFAKSSSIYNPVIYIMMNKQFRNCMLTTLCCGKNPLGDDEASTTGSKTETSQVAPA" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {
                "similarity": "::::::::::::::: :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::.::::.:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::.:::::::::.::::::::::::::::::::::::::::::::::.:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::.:::::::..:::::::::::::::::::::.:::::::::::::.: :::::::::::"
            } == hsp.aln_annotation
        # fourth qresult, second hit
        hit = qresult[1]
        assert "gi|18148870|dbj|AB062417.1|" == hit.id
        assert "Synthetic construct Bos taurus gene for rhodopsin, complete cds" == hit.description
        assert 1047 == hit.seq_len
        assert 1 == len(hit)
        # fourth qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 2237 == hsp.initn_score
        assert 2237 == hsp.init1_score
        assert 2237 == hsp.opt_score
        assert 3067.2 == hsp.z_score
        assert 577.6 == hsp.bitscore
        assert 2.9e-168 == hsp.evalue
        assert 2237 == hsp.sw_score
        assert 93.4 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 98.6 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 348 == hsp.aln_span
        assert 0 == hsp.query_start
        assert 348 == hsp.query_end
        assert "MNGTEGPNFYVPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVLGGFTSTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSRYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAIYNPVIYIMMNKQFRNCMLTTICCGKNPLGDDEASATVSKTETSQVAPA" == hsp.query.seq
        assert 0 == hsp.hit_start
        assert 1044 == hsp.hit_end
        assert "MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # fourth qresult, third hit
        hit = qresult[2]
        assert "gi|283855822|gb|GQ290312.1|" == hit.id
        assert "Myotis ricketti voucher GQX10 rhodopsin (RHO) mRNA, partial cds" == hit.description
        assert 983 == hit.seq_len
        assert 2 == len(hit)
        assert {
                "similarity": "::::::::::::::: ::::::::: ::::::::::::::::::::::.::::::::::::::::::::::::::::::::::::::.::::.:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::.:::::::::.:::::::::: . :.:::::::::::::: ::.:.:::::::::::::::::::::::::::::::::::::::::::::::.:::.:::::::::::.::::::::::::::..:.:::::::::::::::::.::.:::::::::::::.:::::::::::::"
            } == hsp.aln_annotation
        # fourth qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        assert 2138 == hsp.initn_score
        assert 2138 == hsp.init1_score
        assert 2143 == hsp.opt_score
        assert 2939.0 == hsp.z_score
        assert 553.8 == hsp.bitscore
        assert 4.1e-161 == hsp.evalue
        assert 2143 == hsp.sw_score
        assert 95.1 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 99.4 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 326 == hsp.aln_span
        assert 10 == hsp.query_start
        assert 336 == hsp.query_end
        assert "VPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVLGGFTSTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSRYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAIYNPVIYIMMNKQFRNCMLTTICCGKNPLGDDEASAT" == hsp.query.seq
        assert 0 == hsp.hit_start
        assert 978 == hsp.hit_end
        assert "VPFSNKTGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVANLFMVFGGFTTTLYTSMHGYFVFGATGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGLAFTWVMALACAAPPLAGWSRYIPEGMQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVVAFLICWLPYASVAFYIFTHQGSNFGPVFMTIPAFFAKSSSIYNPVIYIMMNKQFRNCMLTTLCCGKNPLGDDEASTT" == hsp.hit.seq
        assert 0 == hsp.query_strand
        assert {
                "similarity": "::::: ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::.::::.::::.:::::.::::::: :::::::::::::::::::::::::::::::::::::::::::::::::.:::::::::::::::::::::::::.::::::::::::::::::::::::::::::::::.::::::::::::::::::::::::::::::::::::::::.::::::.:::::::::::::::::::.:::::::::::..:::::::::::::::::::::.:::::::::::::.:"
            } == hsp.aln_annotation
        # fourth qresult, third hit, second hsp
        hsp = qresult[2].hsps[1]
        assert 74 == hsp.initn_score
        assert 58 == hsp.init1_score
        assert 59 == hsp.opt_score
        assert 91.3 == hsp.z_score
        assert 26.9 == hsp.bitscore
        assert 0.017 == hsp.evalue
        assert 59 == hsp.sw_score
        assert 35.5 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 61.3 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 31 == hsp.aln_span
        assert 234 == hsp.query_start
        assert 265 == hsp.query_end
        assert "AQQQESATTQKAEKEVTRMVIIMVIAFLICW" == hsp.query.seq
        assert 674 == hsp.hit_start
        assert 767 == hsp.hit_end
        assert "SQQIRNATTMMMTMRVTSFSAFWVVADSCCW" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # fourth qresult, fourth hit
        hit = qresult[3]
        assert "gi|2734705|gb|U59921.1|BBU59921" == hit.id
        assert "Bufo bufo rhodopsin mRNA, complete cds" == hit.description
        assert 1574 == hit.seq_len
        assert 1 == len(hit)
        # fourth qresult, fourth hit, first hsp
        hsp = qresult[3].hsps[0]
        assert 2080 == hsp.initn_score
        assert 2031 == hsp.init1_score
        assert 2057 == hsp.opt_score
        assert 2819.6 == hsp.z_score
        assert 532.4 == hsp.bitscore
        assert 1.8e-154 == hsp.evalue
        assert 2057 == hsp.sw_score
        assert 83.3 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 95.2 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 354 == hsp.aln_span
        assert 0 == hsp.query_start
        assert 348 == hsp.query_end
        assert "MNGTEGPNFYVPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVLGGFTSTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSRYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAIYNPVIYIMMNKQFRNCMLTTICCGKNPLGDDEAS-ATVSKTE-----TSQVAPA" == hsp.query.seq
        assert 41 == hsp.hit_start
        assert 1103 == hsp.hit_end
        assert "MNGTEGPNFYIPMSNKTGVVRSPFEYPQYYLAEPWQYSILCAYMFLLILLGFPINFMTLYVTIQHKKLRTPLNYILLNLAFANHFMVLCGFTVTMYSSMNGYFILGATGCYVEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFSENHAVMGVAFTWIMALSCAVPPLLGWSRYIPEGMQCSCGVDYYTLKPEVNNESFVIYMFVVHFTIPLIIIFFCYGRLVCTVKEAAAQQQESATTQKAEKEVTRMVIIMVVFFLICWVPYASVAFFIFSNQGSEFGPIFMTVPAFFAKSSSIYNPVIYIMLNKQFRNCMITTLCCGKNPFGEDDASSAATSKTEASSVSSSQVSPA" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # fourth qresult, fifth hit
        hit = qresult[4]
        assert "gi|12583664|dbj|AB043817.1|" == hit.id
        assert "Conger myriaster conf gene for fresh water form rod opsin, complete cds" == hit.description
        assert 1344 == hit.seq_len
        assert 1 == len(hit)
        # fourth qresult, fifth hit, first hsp
        hsp = qresult[4].hsps[0]
        assert 1975 == hsp.initn_score
        assert 1951 == hsp.init1_score
        assert 1993 == hsp.opt_score
        assert 2732.8 == hsp.z_score
        assert 516.1 == hsp.bitscore
        assert 1.2e-149 == hsp.evalue
        assert 1993 == hsp.sw_score
        assert 81.6 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 93.9 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 347 == hsp.aln_span
        assert 0 == hsp.query_start
        assert 346 == hsp.query_end
        assert "MNGTEGPNFYVPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVLGGFTSTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSRYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAIYNPVIYIMMNKQFRNCMLTTICCGKNPLGD-DEASATVSKTETSQVA" == hsp.query.seq
        assert 22 == hsp.hit_start
        assert 1063 == hsp.hit_end
        assert "MNGTEGPNFYIPMSNATGVVRSPFEYPQYYLAEPWAFSALSAYMFFLIIAGFPINFLTLYVTIEHKKLRTPLNYILLNLAVADLFMVFGGFTTTMYTSMHGYFVFGPTGCNIEGFFATLGGEIALWCLVVLAIERWMVVCKPVTNFRFGESHAIMGVMVTWTMALACALPPLFGWSRYIPEGLQCSCGIDYYTRAPGINNESFVIYMFTCHFSIPLAVISFCYGRLVCTVKEAAAQQQESETTQRAEREVTRMVVIMVISFLVCWVPYASVAWYIFTHQGSTFGPIFMTIPSFFAKSSALYNPMIYICMNKQFRHCMITTLCCGKNPFEEEDGASATSSKTEASSVS" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # fourth qresult, sixth hit
        hit = qresult[5]
        assert "gi|283855845|gb|GQ290303.1|" == hit.id
        assert "Cynopterus brachyotis voucher 20020434 rhodopsin (RHO) gene, exons 1 through 5 and partial cds" == hit.description
        assert 4301 == hit.seq_len
        assert 4 == len(hit)
        # fourth qresult, sixth hit, first hsp
        hsp = qresult[5].hsps[0]
        assert 2094 == hsp.initn_score
        assert 723 == hsp.init1_score
        assert 723 == hsp.opt_score
        assert 992.9 == hsp.z_score
        assert 195.8 == hsp.bitscore
        assert 1e-52 == hsp.evalue
        assert 723 == hsp.sw_score
        assert 96.4 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 99.1 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 111 == hsp.aln_span
        assert 10 == hsp.query_start
        assert 121 == hsp.query_end
        assert "VPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVLGGFTSTLYTSLHGYFVFGPTGCNLEGFFATLGG" == hsp.query.seq
        assert 0 == hsp.hit_start
        assert 333 == hsp.hit_end
        assert "VPFSNKTGVVRSPFEHPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGG" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # fourth qresult, sixth hit, second hsp
        hsp = qresult[5].hsps[1]
        assert 1411 == hsp.initn_score
        assert 499 == hsp.init1_score
        assert 501 == hsp.opt_score
        assert 992.9 == hsp.z_score
        assert 195.8 == hsp.bitscore
        assert 1e-52 == hsp.evalue
        assert 783 == hsp.sw_score
        assert 75.4 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 79.5 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 172 == hsp.aln_span
        assert 176 == hsp.query_start
        assert 312 == hsp.query_end
        assert "RYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKE------------------------------------AAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAIYNPVIYIMMNKQ" == hsp.query.seq
        assert 2854 == hsp.hit_start
        assert 3368 == hsp.hit_end
        assert r"RYIPEGMQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIVIFFCYGQLVFTVKEVRSCVGHWGHAH*VNGAQLHSQSCHSLDT*PCVPA\AAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSNFGPIFMTLPAFFAKSSSIYNPVIYIMMNKQ" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # fourth qresult, sixth hit, third hsp
        hsp = qresult[5].hsps[2]
        assert 431 == hsp.initn_score
        assert 379 == hsp.init1_score
        assert 388 == hsp.opt_score
        assert 992.9 == hsp.z_score
        assert 195.8 == hsp.bitscore
        assert 1e-52 == hsp.evalue
        assert 388 == hsp.sw_score
        assert 80.8 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 90.4 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 73 == hsp.aln_span
        assert 118 == hsp.query_start
        assert 189 == hsp.query_end
        assert "LGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSR--YIPEGLQCSCGI" == hsp.query.seq
        assert 1403 == hsp.hit_start
        assert 1619 == hsp.hit_end
        assert "LAGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGLALTWVMALACAAPPLVGWSR*WH*TEG-KCL*GL" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # fourth qresult, sixth hit, fourth hsp
        hsp = qresult[5].hsps[3]
        assert 213 == hsp.initn_score
        assert 171 == hsp.init1_score
        assert 176 == hsp.opt_score
        assert 992.9 == hsp.z_score
        assert 195.8 == hsp.bitscore
        assert 1e-52 == hsp.evalue
        assert 176 == hsp.sw_score
        assert 76.7 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 93.3 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 30 == hsp.aln_span
        assert 306 == hsp.query_start
        assert 336 == hsp.query_end
        assert "IMMNKQFRNCMLTTICCGKNPLGDDEASAT" == hsp.query.seq
        assert 4206 == hsp.hit_start
        assert 4296 == hsp.hit_end
        assert "MLLAFQFRNCMLTTLCCGKNPLGDDEASTT" == hsp.hit.seq
        assert 0 == hsp.query_strand

    def test_output009(self):
        """Test parsing fasta36 output (output009.m10)."""
        m10_file = get_file("output009.m10")
        qresults = list(parse(m10_file, FMT))
        assert 3 == len(qresults)
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                assert qresult.id == hit.query_id
                for hsp in hit:
                    assert hit.id == hsp.hit_id
                    assert qresult.id == hsp.query_id

        # test first qresult
        qresult = qresults[0]
        assert "random_s00" == qresult.id
        assert "fasta" == qresult.program
        assert "36.3.5c" == qresult.version
        assert "mrnalib.fasta" == qresult.target
        assert 15 == qresult.seq_len
        assert "" == qresult.description
        assert 0 == len(qresult)

        # test second qresult
        qresult = qresults[1]
        assert "gi|255708421:1-99" == qresult.id
        assert "fasta" == qresult.program
        assert "36.3.5c" == qresult.version
        assert "mrnalib.fasta" == qresult.target
        assert 99 == qresult.seq_len
        assert "Mus musculus myoglobin (Mb), transcript varia" == qresult.description
        assert 1 == len(qresult)
        # second qresult, first hit
        hit = qresult[0]
        assert "gi|23308614|ref|NM_152952.1|" == hit.id
        assert "Danio rerio cytoglobin 1 (cygb1), mRNA" == hit.description
        assert 5188 == hit.seq_len
        assert 1 == len(hit)
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 124 == hsp.initn_score
        assert 74 == hsp.init1_score
        assert 74 == hsp.opt_score
        assert 70.3 == hsp.z_score
        assert 23.6 == hsp.bitscore
        assert 1.4 == hsp.evalue
        assert 81.8 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 81.8 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 22 == hsp.aln_span
        assert 7 == hsp.query_start
        assert 29 == hsp.query_end
        assert "TGATGTTCTGTTTCTAAAACAG" == hsp.query.seq
        assert 3483 == hsp.hit_start
        assert 3505 == hsp.hit_end
        assert "TGATTTTTTTTGTCTAAAACAG" == hsp.hit.seq
        assert -1 == hsp.query_strand

        # test third qresult
        qresult = qresults[2]
        assert "gi|156718121:2361-2376" == qresult.id
        assert "fasta" == qresult.program
        assert "36.3.5c" == qresult.version
        assert "mrnalib.fasta" == qresult.target
        assert 16 == qresult.seq_len
        assert "Bos taurus nucleoporin 43kDa (NU" == qresult.description
        assert 5 == len(qresult)
        # third qresult, first hit
        hit = qresult[0]
        assert "gi|47271416|ref|NM_131257.2|" == hit.id
        assert "Danio rerio hemoglobin alpha adult-1 (hbaa1), mRNA" == hit.description
        assert 597 == hit.seq_len
        assert 1 == len(hit)
        # third qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 52 == hsp.initn_score
        assert 52 == hsp.init1_score
        assert 52 == hsp.opt_score
        assert 87.3 == hsp.z_score
        assert 21.0 == hsp.bitscore
        assert 1.4 == hsp.evalue
        assert 85.7 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 85.7 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 14 == hsp.aln_span
        assert 2 == hsp.query_start
        assert 16 == hsp.query_end
        assert "AGAAGGAAAAAAAA" == hsp.query.seq
        assert 572 == hsp.hit_start
        assert 586 == hsp.hit_end
        assert "AGAACTAAAAAAAA" == hsp.hit.seq
        assert 1 == hsp.query_strand
        # third qresult, second hit
        hit = qresult[1]
        assert "gi|332859474|ref|XM_001156938.2|" == hit.id
        assert "PREDICTED: Pan troglodytes myoglobin, transcript variant 11 (MB), mRNA" == hit.description
        assert 762 == hit.seq_len
        assert 1 == len(hit)
        # third qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 52 == hsp.initn_score
        assert 52 == hsp.init1_score
        assert 52 == hsp.opt_score
        assert 85.3 == hsp.z_score
        assert 20.9 == hsp.bitscore
        assert 1.4 == hsp.evalue
        assert 85.7 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 85.7 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 14 == hsp.aln_span
        assert 2 == hsp.query_start
        assert 16 == hsp.query_end
        assert "AGAAGGAAAAAAAA" == hsp.query.seq
        assert 84 == hsp.hit_start
        assert 98 == hsp.hit_end
        assert "AGAAGGTATAAAAA" == hsp.hit.seq
        assert 1 == hsp.query_strand
        # third qresult, third hit
        hit = qresult[2]
        assert "gi|332211534|ref|XM_003254825.1|" == hit.id
        assert "PREDICTED: Nomascus leucogenys hemoglobin subunit gamma-2-like (LOC100581638), mRNA" == hit.description
        assert 805 == hit.seq_len
        assert 1 == len(hit)
        # third qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        assert 52 == hsp.initn_score
        assert 52 == hsp.init1_score
        assert 52 == hsp.opt_score
        assert 84.9 == hsp.z_score
        assert 20.9 == hsp.bitscore
        assert 1.4 == hsp.evalue
        assert 85.7 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 85.7 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 14 == hsp.aln_span
        assert 1 == hsp.query_start
        assert 15 == hsp.query_end
        assert "TTTTTTTCCTTCTT" == hsp.query.seq
        assert 633 == hsp.hit_start
        assert 647 == hsp.hit_end
        assert "TTTTTTTACATCTT" == hsp.hit.seq
        assert -1 == hsp.query_strand
        # third qresult, fourth hit
        hit = qresult[3]
        assert "gi|23308614|ref|NM_152952.1|" == hit.id
        assert "Danio rerio cytoglobin 1 (cygb1), mRNA" == hit.description
        assert 5188 == hit.seq_len
        assert 1 == len(hit)
        # third qresult, fourth hit, first hsp
        hsp = qresult[3].hsps[0]
        assert 52 == hsp.initn_score
        assert 52 == hsp.init1_score
        assert 52 == hsp.opt_score
        assert 70.2 == hsp.z_score
        assert 20.9 == hsp.bitscore
        assert 1.4 == hsp.evalue
        assert 85.7 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 85.7 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 14 == hsp.aln_span
        assert 1 == hsp.query_start
        assert 15 == hsp.query_end
        assert "AAGAAGGAAAAAAA" == hsp.query.seq
        assert 3547 == hsp.hit_start
        assert 3561 == hsp.hit_end
        assert "AATAAGTAAAAAAA" == hsp.hit.seq
        assert 1 == hsp.query_strand
        # third qresult, fifth hit
        hit = qresult[4]
        assert "gi|297689475|ref|XM_002822130.1|" == hit.id
        assert "PREDICTED: Pongo abelii hemoglobin subunit gamma-like (LOC100439631), mRNA" == hit.description
        assert 1158 == hit.seq_len
        assert 2 == len(hit)
        # third qresult, fifth hit, first hsp
        hsp = qresult[4].hsps[0]
        assert 52 == hsp.initn_score
        assert 52 == hsp.init1_score
        assert 52 == hsp.opt_score
        assert 82.0 == hsp.z_score
        assert 20.9 == hsp.bitscore
        assert 1.4 == hsp.evalue
        assert 85.7 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 85.7 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 14 == hsp.aln_span
        assert 1 == hsp.query_start
        assert 15 == hsp.query_end
        assert "TTTTTTTCCTTCTT" == hsp.query.seq
        assert 983 == hsp.hit_start
        assert 997 == hsp.hit_end
        assert "TTTTTTTACATCTT" == hsp.hit.seq
        assert -1 == hsp.query_strand
        # third qresult, fifth hit, second hsp
        hsp = qresult[4].hsps[1]
        assert 52 == hsp.initn_score
        assert 52 == hsp.init1_score
        assert 52 == hsp.opt_score
        assert 82.0 == hsp.z_score
        assert 20.9 == hsp.bitscore
        assert 1.4 == hsp.evalue
        assert 85.7 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 85.7 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 14 == hsp.aln_span
        assert 1 == hsp.query_start
        assert 15 == hsp.query_end
        assert "AAGAAGGAAAAAAA" == hsp.query.seq
        assert 19 == hsp.hit_start
        assert 33 == hsp.hit_end
        assert "AAGAAGGTAAAAGA" == hsp.hit.seq
        assert 1 == hsp.query_strand

    def test_output010(self):
        """Test parsing fasta36 output (output010.m10)."""
        m10_file = get_file("output010.m10")
        qresults = list(parse(m10_file, FMT))
        assert 1 == len(qresults)
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                assert qresult.id == hit.query_id
                for hsp in hit:
                    assert hit.id == hsp.hit_id
                    assert qresult.id == hsp.query_id

        # test first qresult
        qresult = qresults[0]
        assert "random_s00" == qresult.id
        assert "fasta" == qresult.program
        assert "36.3.5c" == qresult.version
        assert "mrnalib.fasta" == qresult.target
        assert 15 == qresult.seq_len
        assert "" == qresult.description
        assert 0 == len(qresult)

    def test_output011(self):
        """Test parsing fasta36 output (output011.m10)."""
        m10_file = get_file("output011.m10")
        qresults = list(parse(m10_file, FMT))
        assert 1 == len(qresults)
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                assert qresult.id == hit.query_id
                for hsp in hit:
                    assert hit.id == hsp.hit_id
                    assert qresult.id == hsp.query_id

        # test first qresult
        qresult = qresults[0]
        assert "gi|255708421:1-99" == qresult.id
        assert "fasta" == qresult.program
        assert "36.3.5c" == qresult.version
        assert "mrnalib.fasta" == qresult.target
        assert 99 == qresult.seq_len
        assert "Mus musculus myoglobin (Mb), transcript varia" == qresult.description
        assert 5 == len(qresult)
        # first qresult, first hit
        hit = qresult[0]
        assert "gi|284005422|ref|NM_001171502.1|" == hit.id
        assert "Oryctolagus cuniculus hemoglobin, zeta (HBZ_2), mRNA" == hit.description
        assert 429 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 79 == hsp.initn_score
        assert 66 == hsp.init1_score
        assert 99 == hsp.opt_score
        assert 90.0 == hsp.z_score
        assert 23.6 == hsp.bitscore
        assert 1.4 == hsp.evalue
        assert 73.8 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 73.8 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 42 == hsp.aln_span
        assert 33 == hsp.query_start
        assert 75 == hsp.query_end
        assert "CAACATCCAGAGGACTGTCATCCTTGTCCCTGTGGGTGAGGG" == hsp.query.seq
        assert 11 == hsp.hit_start
        assert 52 == hsp.hit_end
        assert "CAAGAGCGAGAGGACCATCAT-CATGTCCCTCTGGGACAAGG" == hsp.hit.seq
        assert 1 == hsp.query_strand
        # first qresult, second hit
        hit = qresult[1]
        assert "gi|284005386|ref|NM_001171415.1|" == hit.id
        assert "Oryctolagus cuniculus zeta globin (HBZ0), mRNA" == hit.description
        assert 429 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 73 == hsp.initn_score
        assert 66 == hsp.init1_score
        assert 99 == hsp.opt_score
        assert 90.0 == hsp.z_score
        assert 23.6 == hsp.bitscore
        assert 1.4 == hsp.evalue
        assert 73.8 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 73.8 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 42 == hsp.aln_span
        assert 33 == hsp.query_start
        assert 75 == hsp.query_end
        assert "CAACATCCAGAGGACTGTCATCCTTGTCCCTGTGGGTGAGGG" == hsp.query.seq
        assert 11 == hsp.hit_start
        assert 52 == hsp.hit_end
        assert "CAAGAGCGAGAGGACCATCAT-CATGTCCCTCTGGGACAAGG" == hsp.hit.seq
        assert 1 == hsp.query_strand
        # first qresult, third hit
        hit = qresult[2]
        assert "gi|284005381|ref|NM_001171414.1|" == hit.id
        assert "Oryctolagus cuniculus hemoglobin subunit zeta (RA_M008_JSM295ECF), mRNA" == hit.description
        assert 429 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        assert 90 == hsp.initn_score
        assert 66 == hsp.init1_score
        assert 99 == hsp.opt_score
        assert 90.0 == hsp.z_score
        assert 23.6 == hsp.bitscore
        assert 1.4 == hsp.evalue
        assert 73.8 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 73.8 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 42 == hsp.aln_span
        assert 33 == hsp.query_start
        assert 75 == hsp.query_end
        assert "CAACATCCAGAGGACTGTCATCCTTGTCCCTGTGGGTGAGGG" == hsp.query.seq
        assert 11 == hsp.hit_start
        assert 52 == hsp.hit_end
        assert "CAAGAGCGAGAGGACCATCAT-CATGTCCCTCTGGGACAAGG" == hsp.hit.seq
        assert 1 == hsp.query_strand
        # first qresult, fourth hit
        hit = qresult[3]
        assert "gi|23308614|ref|NM_152952.1|" == hit.id
        assert "Danio rerio cytoglobin 1 (cygb1), mRNA" == hit.description
        assert 5188 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, fourth hit, first hsp
        hsp = qresult[3].hsps[0]
        assert 115 == hsp.initn_score
        assert 78 == hsp.init1_score
        assert 80 == hsp.opt_score
        assert 70.4 == hsp.z_score
        assert 23.6 == hsp.bitscore
        assert 1.4 == hsp.evalue
        assert 80.0 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 80.0 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 25 == hsp.aln_span
        assert 12 == hsp.query_start
        assert 37 == hsp.query_end
        assert "TTAGAAACAGAACATCATCTTCAAC" == hsp.query.seq
        assert 415 == hsp.hit_start
        assert 440 == hsp.hit_end
        assert "TGACAAACTCAACACCATCTTCAAC" == hsp.hit.seq
        assert 1 == hsp.query_strand
        # first qresult, fifth hit
        hit = qresult[4]
        assert "gi|291415427|ref|XM_002723908.1|" == hit.id
        assert "PREDICTED: Oryctolagus cuniculus zeta globin-like (LOC100357627), mRNA" == hit.description
        assert 423 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, fifth hit, first hsp
        hsp = qresult[4].hsps[0]
        assert 91 == hsp.initn_score
        assert 66 == hsp.init1_score
        assert 99 == hsp.opt_score
        assert 90.0 == hsp.z_score
        assert 23.6 == hsp.bitscore
        assert 1.4 == hsp.evalue
        assert 73.8 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 73.8 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 42 == hsp.aln_span
        assert 33 == hsp.query_start
        assert 75 == hsp.query_end
        assert "CAACATCCAGAGGACTGTCATCCTTGTCCCTGTGGGTGAGGG" == hsp.query.seq
        assert 11 == hsp.hit_start
        assert 52 == hsp.hit_end
        assert "CAAGAGCGAGAGGACCATCAT-CATGTCCCTCTGGGACAAGG" == hsp.hit.seq
        assert 1 == hsp.query_strand

    def test_output012(self):
        """Test parsing fasta36 output (output012.m10)."""
        m10_file = get_file("output012.m10")
        qresults = list(parse(m10_file, FMT))
        assert 1 == len(qresults)
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                assert qresult.id == hit.query_id
                for hsp in hit:
                    assert hit.id == hsp.hit_id
                    assert qresult.id == hsp.query_id

        # test first qresult
        qresult = qresults[0]
        assert "gi|156718121:2361-2376" == qresult.id
        assert "fasta" == qresult.program
        assert "36.3.5c" == qresult.version
        assert "mrnalib.fasta" == qresult.target
        assert 16 == qresult.seq_len
        assert "Bos taurus nucleoporin 43kDa (NU" == qresult.description
        assert 8 == len(qresult)
        # first qresult, first hit
        hit = qresult[0]
        assert "gi|77681427|ref|NM_001016495.2|" == hit.id
        assert "Xenopus (Silurana) tropicalis hemoglobin, epsilon 1 (hbe1), mRNA" == hit.description
        assert 554 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 57 == hsp.initn_score
        assert 57 == hsp.init1_score
        assert 57 == hsp.opt_score
        assert 88.1 == hsp.z_score
        assert 21.0 == hsp.bitscore
        assert 1.4 == hsp.evalue
        assert 86.7 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 86.7 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 15 == hsp.aln_span
        assert 1 == hsp.query_start
        assert 16 == hsp.query_end
        assert "AAGAAGGAAAAAAAA" == hsp.query.seq
        assert 529 == hsp.hit_start
        assert 544 == hsp.hit_end
        assert "AAGAATAAAAAAAAA" == hsp.hit.seq
        assert 1 == hsp.query_strand
        # first qresult, second hit
        hit = qresult[1]
        assert "gi|53749657|ref|NM_182940.2|" == hit.id
        assert "Danio rerio hemoglobin alpha embryonic-1 (hbae1), mRNA" == hit.description
        assert 564 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 57 == hsp.initn_score
        assert 57 == hsp.init1_score
        assert 57 == hsp.opt_score
        assert 87.9 == hsp.z_score
        assert 21.0 == hsp.bitscore
        assert 1.4 == hsp.evalue
        assert 86.7 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 86.7 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 15 == hsp.aln_span
        assert 1 == hsp.query_start
        assert 16 == hsp.query_end
        assert "AAGAAGGAAAAAAAA" == hsp.query.seq
        assert 531 == hsp.hit_start
        assert 546 == hsp.hit_end
        assert "AAGAAAAAAAAAAAA" == hsp.hit.seq
        assert 1 == hsp.query_strand
        # first qresult, third hit
        hit = qresult[2]
        assert "gi|158508631|ref|NM_030206.4|" == hit.id
        assert "Mus musculus cytoglobin (Cygb), mRNA" == hit.description
        assert 2331 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        assert 62 == hsp.initn_score
        assert 62 == hsp.init1_score
        assert 62 == hsp.opt_score
        assert 76.3 == hsp.z_score
        assert 20.9 == hsp.bitscore
        assert 1.5 == hsp.evalue
        assert 87.5 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 87.5 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 16 == hsp.aln_span
        assert 0 == hsp.query_start
        assert 16 == hsp.query_end
        assert "GAAGAAGGAAAAAAAA" == hsp.query.seq
        assert 1781 == hsp.hit_start
        assert 1797 == hsp.hit_end
        assert "GAAGAAAAAAAAAAAA" == hsp.hit.seq
        assert 1 == hsp.query_strand
        # first qresult, fourth hit
        hit = qresult[3]
        assert "gi|296217287|ref|XM_002754912.1|" == hit.id
        assert "PREDICTED: Callithrix jacchus hemoglobin subunit gamma-like (LOC100389093), mRNA" == hit.description
        assert 627 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, fourth hit, first hsp
        hsp = qresult[3].hsps[0]
        assert 57 == hsp.initn_score
        assert 57 == hsp.init1_score
        assert 57 == hsp.opt_score
        assert 86.3 == hsp.z_score
        assert 20.8 == hsp.bitscore
        assert 1.5 == hsp.evalue
        assert 86.7 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 86.7 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 15 == hsp.aln_span
        assert 0 == hsp.query_start
        assert 15 == hsp.query_end
        assert "GAAGAAGGAAAAAAA" == hsp.query.seq
        assert 131 == hsp.hit_start
        assert 146 == hsp.hit_end
        assert "GAAGAAGGTAAAACA" == hsp.hit.seq
        assert 1 == hsp.query_strand
        # first qresult, fifth hit
        hit = qresult[4]
        assert "gi|332211540|ref|XM_003254828.1|" == hit.id
        assert "PREDICTED: Nomascus leucogenys hemoglobin subunit gamma-1-like, transcript variant 2 (LOC100582529), mRNA" == hit.description
        assert 640 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, fifth hit, first hsp
        hsp = qresult[4].hsps[0]
        assert 57 == hsp.initn_score
        assert 57 == hsp.init1_score
        assert 57 == hsp.opt_score
        assert 86.0 == hsp.z_score
        assert 20.8 == hsp.bitscore
        assert 1.5 == hsp.evalue
        assert 86.7 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 86.7 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 15 == hsp.aln_span
        assert 1 == hsp.query_start
        assert 16 == hsp.query_end
        assert "TTTTTTTTCCTTCTT" == hsp.query.seq
        assert 589 == hsp.hit_start
        assert 604 == hsp.hit_end
        assert "TTTTTTTTACATCTT" == hsp.hit.seq
        assert -1 == hsp.query_strand
        # first qresult, sixth hit
        hit = qresult[5]
        assert "gi|147903656|ref|NM_001086277.1|" == hit.id
        assert "Xenopus laevis hemoglobin, alpha 2 (hba2), mRNA" == hit.description
        assert 677 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, sixth hit, first hsp
        hsp = qresult[5].hsps[0]
        assert 57 == hsp.initn_score
        assert 57 == hsp.init1_score
        assert 57 == hsp.opt_score
        assert 85.2 == hsp.z_score
        assert 20.8 == hsp.bitscore
        assert 1.6 == hsp.evalue
        assert 86.7 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 86.7 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 15 == hsp.aln_span
        assert 1 == hsp.query_start
        assert 16 == hsp.query_end
        assert "AAGAAGGAAAAAAAA" == hsp.query.seq
        assert 642 == hsp.hit_start
        assert 657 == hsp.hit_end
        assert "AAAAAAGAAAAAAAA" == hsp.hit.seq
        assert 1 == hsp.query_strand
        # first qresult, seventh hit
        hit = qresult[6]
        assert "gi|380013536|ref|XM_003690762.1|" == hit.id
        assert "PREDICTED: Apis florea globin-like (LOC100870092), mRNA" == hit.description
        assert 707 == hit.seq_len
        assert 3 == len(hit)
        # first qresult, seventh hit, first hsp
        hsp = qresult[6].hsps[0]
        assert 57 == hsp.initn_score
        assert 57 == hsp.init1_score
        assert 57 == hsp.opt_score
        assert 84.6 == hsp.z_score
        assert 20.7 == hsp.bitscore
        assert 1.7 == hsp.evalue
        assert 86.7 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 86.7 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 15 == hsp.aln_span
        assert 1 == hsp.query_start
        assert 16 == hsp.query_end
        assert "TTTTTTTTCCTTCTT" == hsp.query.seq
        assert 6 == hsp.hit_start
        assert 21 == hsp.hit_end
        assert "TTTTTTTTTTTTCTT" == hsp.hit.seq
        assert -1 == hsp.query_strand
        # first qresult, seventh hit, second hsp
        hsp = qresult[6].hsps[1]
        assert 57 == hsp.initn_score
        assert 57 == hsp.init1_score
        assert 57 == hsp.opt_score
        assert 84.6 == hsp.z_score
        assert 20.7 == hsp.bitscore
        assert 1.7 == hsp.evalue
        assert 86.7 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 86.7 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 15 == hsp.aln_span
        assert 1 == hsp.query_start
        assert 16 == hsp.query_end
        assert "TTTTTTTTCCTTCTT" == hsp.query.seq
        assert 45 == hsp.hit_start
        assert 60 == hsp.hit_end
        assert "TTTTTTTTTTTTCTT" == hsp.hit.seq
        assert -1 == hsp.query_strand
        # first qresult, seventh hit, third hsp
        hsp = qresult[6].hsps[2]
        assert 57 == hsp.initn_score
        assert 57 == hsp.init1_score
        assert 57 == hsp.opt_score
        assert 84.6 == hsp.z_score
        assert 20.7 == hsp.bitscore
        assert 1.7 == hsp.evalue
        assert 86.7 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 86.7 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 15 == hsp.aln_span
        assert 0 == hsp.query_start
        assert 15 == hsp.query_end
        assert "TTTTTTTCCTTCTTC" == hsp.query.seq
        assert 644 == hsp.hit_start
        assert 659 == hsp.hit_end
        assert "TTTTTTTCCTCCTCC" == hsp.hit.seq
        assert -1 == hsp.query_strand
        # first qresult, eighth hit
        hit = qresult[7]
        assert "gi|332211538|ref|XM_003254827.1|" == hit.id
        assert "PREDICTED: Nomascus leucogenys hemoglobin subunit gamma-1-like, transcript variant 1 (LOC100582529), mRNA" == hit.description
        assert 713 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, eighth hit, first hsp
        hsp = qresult[7].hsps[0]
        assert 57 == hsp.initn_score
        assert 57 == hsp.init1_score
        assert 57 == hsp.opt_score
        assert 84.5 == hsp.z_score
        assert 20.7 == hsp.bitscore
        assert 1.7 == hsp.evalue
        assert 86.7 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 86.7 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 15 == hsp.aln_span
        assert 1 == hsp.query_start
        assert 16 == hsp.query_end
        assert "TTTTTTTTCCTTCTT" == hsp.query.seq
        assert 662 == hsp.hit_start
        assert 677 == hsp.hit_end
        assert "TTTTTTTTACATCTT" == hsp.hit.seq
        assert -1 == hsp.query_strand

    def test_output013(self):
        """Test parsing fasta36 output (output013.m10)."""
        m10_file = get_file("output013.m10")
        qresults = list(parse(m10_file, FMT))
        assert 3 == len(qresults)
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                assert qresult.id == hit.query_id
                for hsp in hit:
                    assert hit.id == hsp.hit_id
                    assert qresult.id == hsp.query_id

        # test first qresult
        qresult = qresults[0]
        assert "random_s00" == qresult.id
        assert "fasta" == qresult.program
        assert "36.3.5c" == qresult.version
        assert "protlib.fasta" == qresult.target
        assert 16 == qresult.seq_len
        assert "" == qresult.description
        assert 0 == len(qresult)

        # test second qresult
        qresult = qresults[1]
        assert "sp|Q9Y2H6|68-133" == qresult.id
        assert "fasta" == qresult.program
        assert "36.3.5c" == qresult.version
        assert "protlib.fasta" == qresult.target
        assert 66 == qresult.seq_len
        assert "" == qresult.description
        assert 1 == len(qresult)
        # second qresult, first hit
        hit = qresult[0]
        assert "gi|291391832|ref|XP_002712264.1|" == hit.id
        assert "PREDICTED: titin [Oryctolagus cuniculus]" == hit.description
        assert 33406 == hit.seq_len
        assert 1 == len(hit)
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 98 == hsp.initn_score
        assert 98 == hsp.init1_score
        assert 109 == hsp.opt_score
        assert 95.1 == hsp.z_score
        assert 30.2 == hsp.bitscore
        assert 0.43 == hsp.evalue
        assert 109 == hsp.sw_score
        assert 26.8 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 54.9 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 71 == hsp.aln_span
        assert 0 == hsp.query_start
        assert 66 == hsp.query_end
        assert "PNGSVPPIY-----VPPGYAPQVIEDNGVRRVVVVPQAPEFHPGSHTVLHRSPHPPLPGFIPVPTMMPPPP" == hsp.query.seq
        assert 10704 == hsp.hit_start
        assert 10775 == hsp.hit_end
        assert "PEKKVPPAVPKKPEAPPAKVPEAPKEVVPEKKIAVPKKPEVPPAKVPEVPKKPVIEEKPVIPVPKKVESPP" == hsp.hit.seq
        assert 0 == hsp.query_strand

        # test third qresult
        qresult = qresults[2]
        assert "sp|Q9Y2H6|265-345" == qresult.id
        assert "fasta" == qresult.program
        assert "36.3.5c" == qresult.version
        assert "protlib.fasta" == qresult.target
        assert 81 == qresult.seq_len
        assert "" == qresult.description
        assert 4 == len(qresult)
        # third qresult, first hit
        hit = qresult[0]
        assert "gi|260806189|ref|XP_002597967.1|" == hit.id
        assert "hypothetical protein BRAFLDRAFT_79792 [Branchiostoma floridae]" == hit.description
        assert 23830 == hit.seq_len
        assert 1 == len(hit)
        # third qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 220 == hsp.initn_score
        assert 62 == hsp.init1_score
        assert 92 == hsp.opt_score
        assert 97.4 == hsp.z_score
        assert 30.5 == hsp.bitscore
        assert 0.32 == hsp.evalue
        assert 92 == hsp.sw_score
        assert 31.6 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 60.8 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 79 == hsp.aln_span
        assert 1 == hsp.query_start
        assert 79 == hsp.query_end
        assert "LSNIVKPVASDIQARTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYVG-EETNITLNDLKPAMDY" == hsp.query.seq
        assert 22430 == hsp.hit_start
        assert 22499 == hsp.hit_end
        assert "VSNI-RPAASDISPHTLTLTWDTP------EDDGGSLITSYVVEMFDVS---DGKWQTLTTTCRRPPYPVKGLNPSATY" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # third qresult, second hit
        hit = qresult[1]
        assert "gi|348553521|ref|XP_003462575.1|" == hit.id
        assert "PREDICTED: receptor-type tyrosine-protein phosphatase F isoform 1 [Cavia porcellus]" == hit.description
        assert 1899 == hit.seq_len
        assert 1 == len(hit)
        # third qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 104 == hsp.initn_score
        assert 75 == hsp.init1_score
        assert 75 == hsp.opt_score
        assert 96.6 == hsp.z_score
        assert 26.7 == hsp.bitscore
        assert 0.36 == hsp.evalue
        assert 75 == hsp.sw_score
        assert 32.4 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 64.9 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 37 == hsp.aln_span
        assert 43 == hsp.query_start
        assert 80 == hsp.query_end
        assert "YEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDYH" == hsp.query.seq
        assert 542 == hsp.hit_start
        assert 579 == hsp.hit_end
        assert "YELVYWAAEEEGQQRKVTFDPTSSYTLEDLKPDTLYH" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # third qresult, third hit
        hit = qresult[2]
        assert "gi|348553523|ref|XP_003462576.1|" == hit.id
        assert "PREDICTED: receptor-type tyrosine-protein phosphatase F isoform 2 [Cavia porcellus]" == hit.description
        assert 1908 == hit.seq_len
        assert 1 == len(hit)
        # third qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        assert 104 == hsp.initn_score
        assert 75 == hsp.init1_score
        assert 75 == hsp.opt_score
        assert 96.6 == hsp.z_score
        assert 26.7 == hsp.bitscore
        assert 0.36 == hsp.evalue
        assert 75 == hsp.sw_score
        assert 32.4 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 64.9 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 37 == hsp.aln_span
        assert 43 == hsp.query_start
        assert 80 == hsp.query_end
        assert "YEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDYH" == hsp.query.seq
        assert 542 == hsp.hit_start
        assert 579 == hsp.hit_end
        assert "YELVYWAAEEEGQQRKVTFDPTSSYTLEDLKPDTLYH" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # third qresult, fourth hit
        hit = qresult[3]
        assert "gi|221124183|ref|XP_002154464.1|" == hit.id
        assert "PREDICTED: similar to FAD104 [Hydra magnipapillata]" == hit.description
        assert 860 == hit.seq_len
        assert 1 == len(hit)
        # third qresult, fourth hit, first hsp
        hsp = qresult[3].hsps[0]
        assert 85 == hsp.initn_score
        assert 66 == hsp.init1_score
        assert 70 == hsp.opt_score
        assert 95.1 == hsp.z_score
        assert 25.3 == hsp.bitscore
        assert 0.43 == hsp.evalue
        assert 70 == hsp.sw_score
        assert 27.1 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 58.6 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 70 == hsp.aln_span
        assert 9 == hsp.query_start
        assert 79 == hsp.query_end
        assert "ASDIQARTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDY" == hsp.query.seq
        assert 615 == hsp.hit_start
        assert 673 == hsp.hit_end
        assert "ASSISYHSIKLKWGHQSS-------KKSI-----LNHTLQMQNKSGSFNTVYSGMDTSFTLSKLKELTPY" == hsp.hit.seq
        assert 0 == hsp.query_strand

    def test_output014(self):
        """Test parsing fasta36 output (output014.m10)."""
        m10_file = get_file("output014.m10")
        qresults = list(parse(m10_file, FMT))
        assert 1 == len(qresults)
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                assert qresult.id == hit.query_id
                for hsp in hit:
                    assert hit.id == hsp.hit_id
                    assert qresult.id == hsp.query_id

        # test first qresult
        qresult = qresults[0]
        assert "random_s00" == qresult.id
        assert "fasta" == qresult.program
        assert "36.3.5c" == qresult.version
        assert "protlib.fasta" == qresult.target
        assert 16 == qresult.seq_len
        assert "" == qresult.description
        assert 0 == len(qresult)

    def test_output015(self):
        """Test parsing fasta36 output (output015.m10)."""
        m10_file = get_file("output015.m10")
        qresults = list(parse(m10_file, FMT))
        assert 1 == len(qresults)
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                assert qresult.id == hit.query_id
                for hsp in hit:
                    assert hit.id == hsp.hit_id
                    assert qresult.id == hsp.query_id

        # test first qresult
        qresult = qresults[0]
        assert "sp|Q9Y2H6|68-133" == qresult.id
        assert "fasta" == qresult.program
        assert "36.3.5c" == qresult.version
        assert "protlib.fasta" == qresult.target
        assert 66 == qresult.seq_len
        assert "" == qresult.description
        assert 2 == len(qresult)
        # first qresult, first hit
        hit = qresult[0]
        assert "gi|194762369|ref|XP_001963317.1|" == hit.id
        assert "GF14002 [Drosophila ananassae]" == hit.description
        assert 1761 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 88 == hsp.initn_score
        assert 68 == hsp.init1_score
        assert 85 == hsp.opt_score
        assert 95.3 == hsp.z_score
        assert 26.0 == hsp.bitscore
        assert 0.42 == hsp.evalue
        assert 85 == hsp.sw_score
        assert 31.0 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 49.3 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 71 == hsp.aln_span
        assert 5 == hsp.query_start
        assert 66 == hsp.query_end
        assert "PPIY----VPPGYA---PQVIEDNGVRRVVVVPQAPEFH---PGSHTVLHRSPHPPLPGFIPVPTMMPPPP" == hsp.query.seq
        assert 128 == hsp.hit_start
        assert 195 == hsp.hit_end
        assert "PPLLQQTATPPQGAQIVPPVCALHHPQQQLALMAAMQHHHPLPPPHA-LHHAPLPPPP---PLPLNPGPPP" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # first qresult, second hit
        hit = qresult[1]
        assert "gi|77812697|ref|NP_035782.3|" == hit.id
        assert "titin isoform N2-A [Mus musculus]" == hit.description
        assert 33467 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 104 == hsp.initn_score
        assert 92 == hsp.init1_score
        assert 106 == hsp.opt_score
        assert 94.9 == hsp.z_score
        assert 30.2 == hsp.bitscore
        assert 0.45 == hsp.evalue
        assert 106 == hsp.sw_score
        assert 29.4 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 57.4 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 68 == hsp.aln_span
        assert 0 == hsp.query_start
        assert 66 == hsp.query_end
        assert "PNGSVPPIY--VPPGYAPQVIEDNGVRRVVVVPQAPEFHPGSHTVLHRSPHPPLPGFIPVPTMMPPPP" == hsp.query.seq
        assert 10780 == hsp.hit_start
        assert 10848 == hsp.hit_end
        assert "PEKKVPPKKPEAPPAKVPEVPKEVVTEKKVAVPKKPEVPPAKVPEVPKKPVIEEKPAIPVVEKVASPP" == hsp.hit.seq
        assert 0 == hsp.query_strand

    def test_output016(self):
        """Test parsing fasta36 output (output016.m10)."""
        m10_file = get_file("output016.m10")
        # only check 8, 10 is too many
        qresults = list(parse(m10_file, FMT))[:8]
        assert 1 == len(qresults)
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                assert qresult.id == hit.query_id
                for hsp in hit:
                    assert hit.id == hsp.hit_id
                    assert qresult.id == hsp.query_id

        # test first qresult
        qresult = qresults[0]
        assert "sp|Q9Y2H6|265-345" == qresult.id
        assert "fasta" == qresult.program
        assert "36.3.5c" == qresult.version
        assert "protlib.fasta" == qresult.target
        assert 81 == qresult.seq_len
        assert "" == qresult.description
        assert 17 == len(qresult)
        # first qresult, first hit
        hit = qresult[0]
        assert "gi|167518632|ref|XP_001743656.1|" == hit.id
        assert "hypothetical protein [Monosiga brevicollis MX1]" == hit.description
        assert 1145 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        assert 88 == hsp.initn_score
        assert 68 == hsp.init1_score
        assert 68 == hsp.opt_score
        assert 97.5 == hsp.z_score
        assert 26.1 == hsp.bitscore
        assert 0.32 == hsp.evalue
        assert 68 == hsp.sw_score
        assert 43.5 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 65.2 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 23 == hsp.aln_span
        assert 56 == hsp.query_start
        assert 79 == hsp.query_end
        assert "YKSVYVGEETNITLNDLKPAMDY" == hsp.query.seq
        assert 424 == hsp.hit_start
        assert 447 == hsp.hit_end
        assert "FRPVYTGIDTNYKVVDLTPNCDY" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # first qresult, second hit
        hit = qresult[1]
        assert "gi|9507013|ref|NP_062122.1|" == hit.id
        assert "receptor-type tyrosine-protein phosphatase F precursor [Rattus norvegicus]" == hit.description
        assert 1898 == hit.seq_len
        assert 2 == len(hit)
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        assert 83 == hsp.initn_score
        assert 43 == hsp.init1_score
        assert 72 == hsp.opt_score
        assert 96.2 == hsp.z_score
        assert 26.6 == hsp.bitscore
        assert 0.37 == hsp.evalue
        assert 72 == hsp.sw_score
        assert 26.8 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 54.9 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 71 == hsp.aln_span
        assert 8 == hsp.query_start
        assert 79 == hsp.query_end
        assert "VASDIQARTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDY" == hsp.query.seq
        assert 325 == hsp.hit_start
        assert 385 == hsp.hit_end
        assert "VVTETTATSVTLTWD------SGNTEPVS---FYG--IQYRAAGTDGPFQEVDGVASTRYSIGGLSPFSEY" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # first qresult, second hit, second hsp
        hsp = qresult[1].hsps[1]
        assert 98 == hsp.initn_score
        assert 70 == hsp.init1_score
        assert 70 == hsp.opt_score
        assert 96.2 == hsp.z_score
        assert 26.6 == hsp.bitscore
        assert 0.37 == hsp.evalue
        assert 70 == hsp.sw_score
        assert 32.4 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 62.2 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 37 == hsp.aln_span
        assert 43 == hsp.query_start
        assert 80 == hsp.query_end
        assert "YEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDYH" == hsp.query.seq
        assert 542 == hsp.hit_start
        assert 579 == hsp.hit_end
        assert "YELVYWAAEDEGQQHKVTFDPTSSYTLEDLKPDTLYH" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # first qresult, third hit
        hit = qresult[2]
        assert "gi|115648048|ref|NP_035343.2|" == hit.id
        assert "receptor-type tyrosine-protein phosphatase F precursor [Mus musculus]" == hit.description
        assert 1898 == hit.seq_len
        assert 2 == len(hit)
        # first qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        assert 98 == hsp.initn_score
        assert 70 == hsp.init1_score
        assert 73 == hsp.opt_score
        assert 96.2 == hsp.z_score
        assert 26.6 == hsp.bitscore
        assert 0.37 == hsp.evalue
        assert 73 == hsp.sw_score
        assert 25.6 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 54.9 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 82 == hsp.aln_span
        assert 7 == hsp.query_start
        assert 80 == hsp.query_end
        assert "PVASDIQARTVVLTWSPPSSL-INGETDESS-----VP---ELYGYEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDYH" == hsp.query.seq
        assert 497 == hsp.hit_start
        assert 579 == hsp.hit_end
        assert "PPSPTIQVKTQQGVPAQPADFQANAESDTRIQLSWLLPPQERIVKYELVYWAAEDEGQQHKVTFDPTSSYTLEDLKPDTLYH" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # first qresult, third hit, second hsp
        hsp = qresult[2].hsps[1]
        assert 76 == hsp.initn_score
        assert 43 == hsp.init1_score
        assert 72 == hsp.opt_score
        assert 96.2 == hsp.z_score
        assert 26.6 == hsp.bitscore
        assert 0.37 == hsp.evalue
        assert 72 == hsp.sw_score
        assert 26.8 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 54.9 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 71 == hsp.aln_span
        assert 8 == hsp.query_start
        assert 79 == hsp.query_end
        assert "VASDIQARTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDY" == hsp.query.seq
        assert 325 == hsp.hit_start
        assert 385 == hsp.hit_end
        assert "VVTETTATSVTLTWD------SGNTEPVS---FYG--IQYRAAGTDGPFQEVDGVASTRYSIGGLSPFSEY" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # first qresult, fourth hit
        hit = qresult[3]
        assert "gi|354481005|ref|XP_003502693.1|" == hit.id
        assert "PREDICTED: LOW QUALITY PROTEIN: receptor-type tyrosine-protein phosphatase F-like [Cricetulus griseus]" == hit.description
        assert 1898 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, fourth hit, first hsp
        hsp = qresult[3].hsps[0]
        assert 98 == hsp.initn_score
        assert 70 == hsp.init1_score
        assert 70 == hsp.opt_score
        assert 96.2 == hsp.z_score
        assert 26.6 == hsp.bitscore
        assert 0.37 == hsp.evalue
        assert 70 == hsp.sw_score
        assert 32.4 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 62.2 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 37 == hsp.aln_span
        assert 43 == hsp.query_start
        assert 80 == hsp.query_end
        assert "YEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDYH" == hsp.query.seq
        assert 542 == hsp.hit_start
        assert 579 == hsp.hit_end
        assert "YELVYWAAEDEGQQHKVTFDPTSSYTLEDLKPDTVYH" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # first qresult, fifth hit
        hit = qresult[4]
        assert "gi|328789682|ref|XP_003251305.1|" == hit.id
        assert "PREDICTED: LOW QUALITY PROTEIN: twitchin [Apis mellifera]" == hit.description
        assert 8619 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, fifth hit, first hsp
        hsp = qresult[4].hsps[0]
        assert 70 == hsp.initn_score
        assert 70 == hsp.init1_score
        assert 78 == hsp.opt_score
        assert 95.2 == hsp.z_score
        assert 28.6 == hsp.bitscore
        assert 0.43 == hsp.evalue
        assert 78 == hsp.sw_score
        assert 28.6 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 54.3 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 70 == hsp.aln_span
        assert 9 == hsp.query_start
        assert 79 == hsp.query_end
        assert "ASDIQARTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDY" == hsp.query.seq
        assert 4760 == hsp.hit_start
        assert 4823 == hsp.hit_end
        assert "ASDVHAEGCTLTWKPP------EDDGGQPIDKYVVEKMDEATGRWVPAGETD-GPQTSLQVEGLTPGHKY" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # first qresult, sixth hit
        hit = qresult[5]
        assert "gi|260828627|ref|XP_002609264.1|" == hit.id
        assert "hypothetical protein BRAFLDRAFT_124749 [Branchiostoma floridae]" == hit.description
        assert 4389 == hit.seq_len
        assert 7 == len(hit)
        # first qresult, sixth hit, first hsp
        hsp = qresult[5].hsps[0]
        assert 81 == hsp.initn_score
        assert 73 == hsp.init1_score
        assert 97 == hsp.opt_score
        assert 95.1 == hsp.z_score
        assert 27.6 == hsp.bitscore
        assert 0.43 == hsp.evalue
        assert 97 == hsp.sw_score
        assert 21.4 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 67.1 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 70 == hsp.aln_span
        assert 9 == hsp.query_start
        assert 79 == hsp.query_end
        assert "ASDIQARTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDY" == hsp.query.seq
        assert 2241 == hsp.hit_start
        assert 2302 == hsp.hit_end
        assert "ANAVDSQSIRINWQPPTE-PNGN--------VLGYNIFYTTEGESGNNQQTVGPDDTTYVIEGLRPATQY" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # first qresult, sixth hit, second hsp
        hsp = qresult[5].hsps[1]
        assert 177 == hsp.initn_score
        assert 55 == hsp.init1_score
        assert 90 == hsp.opt_score
        assert 95.1 == hsp.z_score
        assert 27.6 == hsp.bitscore
        assert 0.43 == hsp.evalue
        assert 90 == hsp.sw_score
        assert 30.6 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 56.9 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 72 == hsp.aln_span
        assert 8 == hsp.query_start
        assert 79 == hsp.query_end
        assert "VASDIQA-RTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDY" == hsp.query.seq
        assert 2818 == hsp.hit_start
        assert 2881 == hsp.hit_end
        assert "VTADGQAPDTVVVTWQSPAET-NGD--------LLGYYIYYQVVGSTETSQAETGPDETTYSISGLRPATEY" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # first qresult, sixth hit, third hsp
        hsp = qresult[5].hsps[2]
        assert 196 == hsp.initn_score
        assert 61 == hsp.init1_score
        assert 84 == hsp.opt_score
        assert 95.1 == hsp.z_score
        assert 27.6 == hsp.bitscore
        assert 0.43 == hsp.evalue
        assert 84 == hsp.sw_score
        assert 27.8 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 56.9 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 72 == hsp.aln_span
        assert 8 == hsp.query_start
        assert 79 == hsp.query_end
        assert "VASDIQA-RTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDY" == hsp.query.seq
        assert 3300 == hsp.hit_start
        assert 3363 == hsp.hit_end
        assert "VTAEGQAPDTITVTWQSPAET-NGD--------LLGYYIYYQVVGSTEDVRAEAGPEETTYSISGLRPATEY" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # first qresult, sixth hit, fourth hsp
        hsp = qresult[5].hsps[3]
        assert 79 == hsp.initn_score
        assert 49 == hsp.init1_score
        assert 83 == hsp.opt_score
        assert 95.1 == hsp.z_score
        assert 27.6 == hsp.bitscore
        assert 0.43 == hsp.evalue
        assert 83 == hsp.sw_score
        assert 27.9 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 57.4 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 68 == hsp.aln_span
        assert 12 == hsp.query_start
        assert 80 == hsp.query_end
        assert "IQARTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDYH" == hsp.query.seq
        assert 3686 == hsp.hit_start
        assert 3747 == hsp.hit_end
        assert "IDSTTIELQWMPPSP------DEQN-GVIKGYKILYKKVGEEGENEEDAGLLDLMYTLSDLEKWTEYN" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # first qresult, sixth hit, fifth hsp
        hsp = qresult[5].hsps[4]
        assert 100 == hsp.initn_score
        assert 50 == hsp.init1_score
        assert 81 == hsp.opt_score
        assert 95.1 == hsp.z_score
        assert 27.6 == hsp.bitscore
        assert 0.43 == hsp.evalue
        assert 81 == hsp.sw_score
        assert 25.7 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 57.1 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 70 == hsp.aln_span
        assert 9 == hsp.query_start
        assert 79 == hsp.query_end
        assert "ASDIQARTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDY" == hsp.query.seq
        assert 3398 == hsp.hit_start
        assert 3459 == hsp.hit_end
        assert "ASSLGSEAIEVSWQPPPQS-NGE--------ILGYRLHYQIVGEESASTQEVEGYETFYLLRGLRPVTEY" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # first qresult, sixth hit, sixth hsp
        hsp = qresult[5].hsps[5]
        assert 178 == hsp.initn_score
        assert 58 == hsp.init1_score
        assert 81 == hsp.opt_score
        assert 95.1 == hsp.z_score
        assert 27.6 == hsp.bitscore
        assert 0.43 == hsp.evalue
        assert 81 == hsp.sw_score
        assert 27.1 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 55.7 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 70 == hsp.aln_span
        assert 9 == hsp.query_start
        assert 79 == hsp.query_end
        assert "ASDIQARTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDY" == hsp.query.seq
        assert 2145 == hsp.hit_start
        assert 2206 == hsp.hit_end
        assert "ATPVDPRTVRVEWQPPQQ-PNGE--------IQGYNIYYRTTESDEDALQQAGAQDIFLTLTGLSPFTEY" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # first qresult, sixth hit, seventh hsp
        hsp = qresult[5].hsps[6]
        assert 102 == hsp.initn_score
        assert 48 == hsp.init1_score
        assert 79 == hsp.opt_score
        assert 95.1 == hsp.z_score
        assert 27.6 == hsp.bitscore
        assert 0.43 == hsp.evalue
        assert 79 == hsp.sw_score
        assert 29.4 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 54.4 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 68 == hsp.aln_span
        assert 12 == hsp.query_start
        assert 79 == hsp.query_end
        assert "IQARTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYVGE-ETNITLNDLKPAMDY" == hsp.query.seq
        assert 3497 == hsp.hit_start
        assert 3555 == hsp.hit_end
        assert "VEPTTITVDWQPPLE-INGV--------LLGYKVIYMPENA-AEFSTVELGPAELSTMLLDLEPATTY" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # first qresult, seventh hit
        hit = qresult[6]
        assert "gi|119220552|ref|NP_689957.3|" == hit.id
        assert "protein sidekick-1 isoform 1 [Homo sapiens]" == hit.description
        assert 2213 == hit.seq_len
        assert 1 == len(hit)
        # first qresult, seventh hit, first hsp
        hsp = qresult[6].hsps[0]
        assert 87 == hsp.initn_score
        assert 51 == hsp.init1_score
        assert 70 == hsp.opt_score
        assert 95.0 == hsp.z_score
        assert 26.6 == hsp.bitscore
        assert 0.43 == hsp.evalue
        assert 70 == hsp.sw_score
        assert 29.9 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 58.2 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 67 == hsp.aln_span
        assert 8 == hsp.query_start
        assert 73 == hsp.query_end
        assert "VASDIQARTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYV-GEETNITL-NDL" == hsp.query.seq
        assert 775 == hsp.hit_start
        assert 835 == hsp.hit_end
        assert "VASGRTNQSIMVQWQPPP-----ETEHNGV--LRGYILRYRLAGLPGEYQQRNITSPEVNYCLVTDL" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # first qresult, eighth hit
        hit = qresult[7]
        assert "gi|332864595|ref|XP_518946.3|" == hit.id
        assert "PREDICTED: protein sidekick-1 [Pan troglodytes]" == hit.description
        assert 2213 == hit.seq_len
        assert 2 == len(hit)
        # first qresult, eighth hit, first hsp
        hsp = qresult[7].hsps[0]
        assert 68 == hsp.initn_score
        assert 45 == hsp.init1_score
        assert 76 == hsp.opt_score
        assert 95.0 == hsp.z_score
        assert 26.6 == hsp.bitscore
        assert 0.43 == hsp.evalue
        assert 76 == hsp.sw_score
        assert 32.5 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 61.0 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 77 == hsp.aln_span
        assert 4 == hsp.query_start
        assert 80 == hsp.query_end
        assert "IVKPVASDIQARTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYVGEE-TNITLNDLKPAMDYH" == hsp.query.seq
        assert 674 == hsp.hit_start
        assert 740 == hsp.hit_end
        assert "LASPNSS--HSHAVVLSWVRP---FDGNS-----PILY-YIVELSENNSPWKVHLSNVGPEMTGITVSGLTPARTYQ" == hsp.hit.seq
        assert 0 == hsp.query_strand
        # first qresult, eighth hit, second hsp
        hsp = qresult[7].hsps[1]
        assert 87 == hsp.initn_score
        assert 51 == hsp.init1_score
        assert 70 == hsp.opt_score
        assert 95.0 == hsp.z_score
        assert 26.6 == hsp.bitscore
        assert 0.43 == hsp.evalue
        assert 70 == hsp.sw_score
        assert 29.9 == pytest.approx(hsp.ident_pct, abs=5e-8)
        assert 58.2 == pytest.approx(hsp.pos_pct, abs=5e-8)
        assert 67 == hsp.aln_span
        assert 8 == hsp.query_start
        assert 73 == hsp.query_end
        assert "VASDIQARTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYV-GEETNITL-NDL" == hsp.query.seq
        assert 775 == hsp.hit_start
        assert 835 == hsp.hit_end
        assert "VASGRTNQSIMVQWQPPP-----ETEHNGV--LRGYILRYRLAGLPGEYQQRNITSPEVNYCLVTDL" == hsp.hit.seq
        assert 0 == hsp.query_strand


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
