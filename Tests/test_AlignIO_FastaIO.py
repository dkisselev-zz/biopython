# Copyright 2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Unit tests for the Bio.AlignIO.FastaIO module."""

import unittest
import pytest

from Bio.AlignIO import FastaIO


class FastaIOTests(unittest.TestCase):
    """Test FastaIO module."""

    def test_output001(self):
        """Check output001.m10 file."""
        fasta_file = "Fasta/output001.m10"
        with open(fasta_file) as handle:
            alignments = list(FastaIO.FastaM10Iterator(handle))
            assert len(alignments) == 4
            assert len(alignments[0]) == 2
            assert alignments[0].get_alignment_length() == 108
            assert (alignments[0][0].seq == "SGSNT-RRRAISRPVRLTAEED---QE"
                "IRKRAAECGKTVSGFLRAAALGKKVNS"
                "LTDDRVLKEVM-----RLGALQKKLFI"
                "DGKRVGDREYAEVLIAITEYHRALLSR")
            assert alignments[0][0].id == "gi|10955263|ref|NP_052604.1|"
            assert alignments[0][0].annotations["original_length"] == 107
            assert (alignments[0][1].seq == "AGSGAPRRRGSGLASRISEQSEALLQE"
                "AAKHAAEFGRS------EVDTEHLLLA"
                "LADSDVVKTILGQFKIKVDDLKRQIES"
                "EAKR-GDKPF-EGEIGVSPRVKDALSR")
            assert alignments[0][1].id == "gi|152973457|ref|YP_001338508.1|"
            assert alignments[0][1].annotations["original_length"] == 931
            assert len(alignments[1]) == 2
            assert alignments[1].get_alignment_length() == 64
            assert alignments[1][0].seq == "AAECGKTVSGFLRAAALGKKVNSLTDDRVLKEV-MRLGALQKKLFIDGKRVGDREYAEVLIAIT"
            assert alignments[1][0].id == "gi|10955263|ref|NP_052604.1|"
            assert alignments[1][0].annotations["original_length"] == 107
            assert alignments[1][1].seq == "ASRQGCTVGG--KMDSVQDKASDKDKERVMKNINIMWNALSKNRLFDG----NKELKEFIMTLT"
            assert alignments[1][1].id == "gi|152973588|ref|YP_001338639.1|"
            assert alignments[1][1].annotations["original_length"] == 459
            assert len(alignments[2]) == 2
            assert alignments[2].get_alignment_length() == 38
            assert alignments[2][0].seq == "MKKDKKYQIEAIKNKDKTLFIVYATDIYSPSEFFSKIE"
            assert alignments[2][0].id == "gi|10955264|ref|NP_052605.1|"
            assert alignments[2][0].annotations["original_length"] == 126
            assert alignments[2][1].seq == "IKKDLGVSFLKLKNREKTLIVDALKKKYPVAELLSVLQ"
            assert alignments[2][1].id == "gi|152973462|ref|YP_001338513.1|"
            assert alignments[2][1].annotations["original_length"] == 101
            assert len(alignments[3]) == 2
            assert alignments[3].get_alignment_length() == 43
            assert alignments[3][0].seq == "SELHSKLPKSIDKIHEDIKKQLSC-SLIMKKIDVEMEDYSTYC"
            assert alignments[3][0].id == "gi|10955265|ref|NP_052606.1|"
            assert alignments[3][0].annotations["original_length"] == 346
            assert alignments[3][1].seq == "SRINSDVARRIPGIHRDPKDRLSSLKQVEEALDMLISSHGEYC"
            assert alignments[3][1].id == "gi|152973545|ref|YP_001338596.1|"
            assert alignments[3][1].annotations["original_length"] == 242

    def test_output002(self):
        """Check output002.m10 file."""
        fasta_file = "Fasta/output002.m10"
        with open(fasta_file) as handle:
            alignments = list(FastaIO.FastaM10Iterator(handle))
            assert len(alignments) == 6
            assert len(alignments[0]) == 2
            assert alignments[0].get_alignment_length() == 88
            assert (alignments[0][0].seq == "SGSNTRRRAISRPVR--LTAEED"
                "QEIRKRAAECG-KTVSGFLRAAA"
                "LGKKVNSLTDDRVLKEVMRLGAL"
                "QKKLFIDGKRVGDREYAEV")
            assert alignments[0][0].id == "gi|10955263|ref|NP_052604.1|"
            assert alignments[0][0].annotations["original_length"] == 107
            assert (alignments[0][1].seq == "SQRSTRRKPENQPTRVILFNKPY"
                "DVLPQFTDEAGRKTLKEFIPVQG"
                "VYAAGRLDRDSEGLLVLTNNGAL"
                "QARLTQPGKRTGKIYYVQV")
            assert alignments[0][1].id == "gi|162139799|ref|NP_309634.2|"
            assert alignments[0][1].annotations["original_length"] == 207
            assert len(alignments[1]) == 2
            assert alignments[1].get_alignment_length() == 53
            assert alignments[1][0].seq == "EIRKRAAECGKTVSGFLRAAA-LGKKV----NSLTDDRVLKEVMRLGALQKKL"
            assert alignments[1][0].id == "gi|10955263|ref|NP_052604.1|"
            assert alignments[1][0].annotations["original_length"] == 107
            assert alignments[1][1].seq == "EIKPRGTSKGEAIAAFMQEAPFIGRTPVFLGDDLTDESGFAVVNRLGGMSVKI"
            assert alignments[1][1].id == "gi|15831859|ref|NP_310632.1|"
            assert alignments[1][1].annotations["original_length"] == 266
            assert len(alignments[2]) == 2
            assert alignments[2].get_alignment_length() == 92
            assert (alignments[2][0].seq == "SEFFSKIESDLKKKKSKGDVFFD"
                "LIIPNG-----GKKDRYVYTSFN"
                "GEKFSSYTLNKVTKTDEYNDLSE"
                "LSASFFKKNFDKINVNLLSKATS")
            assert alignments[2][0].id == "gi|10955264|ref|NP_052605.1|"
            assert alignments[2][0].annotations["original_length"] == 126
            assert (alignments[2][1].seq == "TELNSELAKAMKVDAQRG-AFVS"
                "QVLPNSSAAKAGIKAGDVITSL"
                "NGKPISSFAALRA-QVGTMPVG"
                "SKLTLGLLRDG-KQVNVNLELQ"
                "QSS")
            assert alignments[2][1].id == "gi|15829419|ref|NP_308192.1|"
            assert alignments[2][1].annotations["original_length"] == 474
            assert len(alignments[3]) == 2
            assert alignments[3].get_alignment_length() == 73
            assert (alignments[3][0].seq == "FFDLIIPNGGKKDRYVYTSFNGE"
                "KFSSYTLNKVTKTDEYNDLSELS"
                "ASFFKKNFDKINVNLLSKATSFA"
                "LKKG")
            assert alignments[3][0].id == "gi|10955264|ref|NP_052605.1|"
            assert alignments[3][0].annotations["original_length"] == 126
            assert (alignments[3][1].seq == "LFDLFLKNDAMHDPMVNESYC-E"
                "TFGWVSKENLARMKE---LTYKA"
                "NDVLKKLFDDAGLILVDFKLEFG"
                "LYKG")
            assert alignments[3][1].id == "gi|15832592|ref|NP_311365.1|"
            assert alignments[3][1].annotations["original_length"] == 237
            assert len(alignments[4]) == 2
            assert alignments[4].get_alignment_length() == 63
            assert alignments[4][0].seq == "VDIKK-ETIESELHSKLPKSIDKIHEDIKKQLSCSLI--MKKID-VEMEDYSTYCFSALRAIE"
            assert alignments[4][0].id == "gi|10955265|ref|NP_052606.1|"
            assert alignments[4][0].annotations["original_length"] == 346
            assert alignments[4][1].seq == "IDPKKIEQIARQVHESMPKGIREFGEDVEKKIRQTLQAQLTRLDLVSREEFDVQTQVLLRTRE"
            assert alignments[4][1].id == "gi|38704138|ref|NP_311957.2|"
            assert alignments[4][1].annotations["original_length"] == 111
            assert len(alignments[5]) == 2
            assert alignments[5].get_alignment_length() == 157
            assert (alignments[5][0].seq == "QYIMTTSNGDRVRAKIYKRGSIQ"
                "FQGKYLQIASLINDFMCSILNMK"
                "EIVEQKNKEFNVDI---KKETI-"
                "ESELHSKLPKSIDKIHEDIKKQL"
                "SCSLIMKKIDV-EMEDYSTYCFS"
                "ALRA-IEGFIYQILNDVCNPSSS"
                "KNLGEYFTENKPKYIIREI")
            assert alignments[5][0].id == "gi|10955265|ref|NP_052606.1|"
            assert alignments[5][0].annotations["original_length"] == 346
            assert (alignments[5][1].seq == "EFIRLLSDHDQFEKDQISELTVA"
                "ANALKLEVAK--NNY-----NMK"
                "YSFDTQTERRMIELIREQKDLIP"
                "EKYLHQSGIKKL-KLHED---EF"
                "SSLLVDAERQVLEGSSFVLCCGE"
                "KINSTISELLSKKITDLTHPTES"
                "FTLSEYFSYDVYEEIFKKV")
            assert alignments[5][1].id == "gi|15833861|ref|NP_312634.1|"
            assert alignments[5][1].annotations["original_length"] == 330

    def test_output003(self):
        """Check output003.m10 file."""
        fasta_file = "Fasta/output003.m10"
        with open(fasta_file) as handle:
            alignments = list(FastaIO.FastaM10Iterator(handle))
            assert len(alignments) == 3
            assert len(alignments[0]) == 2
            assert alignments[0].get_alignment_length() == 55
            assert alignments[0][0].seq == "ISISNNKDQYEELQKEQGERDLKTVDQLVRIAAAGGGLRLSASTKTVDQLVRIAA"
            assert alignments[0][0].id == "gi|152973837|ref|YP_001338874.1|"
            assert alignments[0][0].annotations["original_length"] == 183
            assert alignments[0][1].seq == "VRLTAEEDQ--EIRKRAAECG-KTVSGFLRAAALGKKVNSLTDDRVLKEVMRLGA"
            assert alignments[0][1].id == "gi|10955263|ref|NP_052604.1|"
            assert alignments[0][1].annotations["original_length"] == 107
            assert len(alignments[1]) == 2
            assert alignments[1].get_alignment_length() == 22
            assert alignments[1][0].seq == "DDAEHLFRTLSSR-LDALQDGN"
            assert alignments[1][0].id == "gi|152973840|ref|YP_001338877.1|"
            assert alignments[1][0].annotations["original_length"] == 63
            assert alignments[1][1].seq == "DDRANLFEFLSEEGITITEDNN"
            assert alignments[1][1].id == "gi|10955265|ref|NP_052606.1|"
            assert alignments[1][1].annotations["original_length"] == 346
            assert len(alignments[2]) == 2
            assert alignments[2].get_alignment_length() == 63
            assert alignments[2][0].seq == "VFGSFEQPKGEHLSGQVSEQ--RDTAFADQNEQVIRHLKQEIEHLNTLLLSKDSHIDSLKQAM"
            assert alignments[2][0].id == "gi|152973841|ref|YP_001338878.1|"
            assert alignments[2][0].annotations["original_length"] == 133
            assert alignments[2][1].seq == "VYTSFN---GEKFSSYTLNKVTKTDEYNDLSELSASFFKKNFDKINVNLLSKATSF-ALKKGI"
            assert alignments[2][1].id == "gi|10955264|ref|NP_052605.1|"
            assert alignments[2][1].annotations["original_length"] == 126

    def test_output004(self):
        """Check output004.m10 file."""
        fasta_file = "Fasta/output004.m10"
        with open(fasta_file) as handle:
            alignments = list(FastaIO.FastaM10Iterator(handle))
            assert len(alignments) == 1
            assert len(alignments[0]) == 2
            assert alignments[0].get_alignment_length() == 102
            assert (alignments[0][0].seq == "AAAAAAGATAAAAAATATCAAAT"
                "AGAAGCAATAAAAAATAAAGATA"
                "AAACTTTATTTATTGTCTATGCT"
                "ACTGATATTTATAGCCCGAGCGA"
                "ATTTTTCTCA")
            assert alignments[0][0].id == "ref|NC_002127.1|:c1351-971"
            assert alignments[0][0].annotations["original_length"] == 381
            assert (alignments[0][1].seq == "AGAGAAAATAAAACAAGTAATAA"
                "AATATTAATGGAAAAAATAAATT"
                "CTTGTTTATTTAGACCTGATTCT"
                "AATCACTTTTCTTGCCCGGAGTC"
                "ATTTTTGACA")
            assert alignments[0][1].id == "ref|NC_002695.1|:1970775-1971404"
            assert alignments[0][1].annotations["original_length"] == 630

    def test_output005(self):
        """Check output005.m10 file."""
        fasta_file = "Fasta/output005.m10"
        with open(fasta_file) as handle:
            alignments = list(FastaIO.FastaM10Iterator(handle))
            assert len(alignments) == 1
            assert len(alignments[0]) == 2
            assert alignments[0].get_alignment_length() == 110
            assert (alignments[0][0].seq == "IKNKDKTLFIVYAT-DIYSPSEF"
                "FSKIESDLKKKKSKGDV--FFDL"
                "IIPNGGKKD--RYVYTSFNGEKF"
                "SSYTLNKVTKTDEYNDL--SELS"
                "ASFFKKNFDKINVNLLSK")
            assert alignments[0][0].id == "gi|10955264|ref|NP_052605.1|"
            assert alignments[0][0].annotations["original_length"] == 126
            assert (alignments[0][1].seq == "IKDELPVAFCSWASLDLECEVKY"
                "INDVTSLYAKDWMSGERKWFIDW"
                "IAPFGHNMELYKYMRKKYPYELF"
                "RAIRLDESSKTGKIAEFHGGGID"
                "KKLASKIFRQYHHELMSE")
            assert alignments[0][1].id == "gi|10955282|ref|NP_052623.1|"
            assert alignments[0][1].annotations["original_length"] == 163

    def test_output006(self):
        """Check output006.m10 file."""
        fasta_file = "Fasta/output006.m10"
        with open(fasta_file) as handle:
            alignments = list(FastaIO.FastaM10Iterator(handle))
            assert len(alignments) == 1
            assert len(alignments[0]) == 2
            assert alignments[0].get_alignment_length() == 131
            assert (alignments[0][0].seq == "GCAACGCTTCAAGAACTGGAATT"
                "AGGAACCGTGACAACGATTAATG"
                "AGGAGATTTATGAAGAGGGTTCT"
                "TCGATTTTAGGCCAATCGGAAGG"
                "AATTATGTAGCAAGTCCATCAGA"
                "AAATGGAAGAAGTCAT")
            assert alignments[0][0].id == "query"
            assert alignments[0][0].annotations["original_length"] == 131
            assert (alignments[0][1].seq == "GCAACGCTTCAAGAACTGGAATT"
                "AGGAACCGTGACAACGATTAATG"
                "AGGAGATTTATGAAGAGGGTTCT"
                "TCGATTTTAGGCCAATCGGAAGG"
                "AATTATGTAGCAAGTCCATCAGA"
                "AAATGGAAGTAGTCAT")
            assert alignments[0][1].id == "gi|116660610|gb|EG558221.1|EG558221"
            assert alignments[0][1].annotations["original_length"] == 573

    def test_output007(self):
        """Check output007.m10 file."""
        fasta_file = "Fasta/output007.m10"
        with open(fasta_file) as handle:
            alignments = list(FastaIO.FastaM10Iterator(handle))
            assert len(alignments) == 9
            assert len(alignments[0]) == 2
            assert alignments[0].get_alignment_length() == 108
            assert (alignments[0][0].seq == "SGSNT-RRRAISRPVRLTAEED-"
                "--QEIRKRAAECGKTVSGFLRAA"
                "ALGKKVNSLTDDRVLKEVM----"
                "-RLGALQKKLFIDGKRVGDREYA"
                "EVLIAITEYHRALLSR")
            assert alignments[0][0].id == "gi|10955263|ref|NP_052604.1|"
            assert alignments[0][0].annotations["original_length"] == 107
            assert (alignments[0][1].seq == "AGSGAPRRRGSGLASRISEQSEA"
                "LLQEAAKHAAEFGRS------EV"
                "DTEHLLLALADSDVVKTILGQFK"
                "IKVDDLKRQIESEAKR-GDKPF-"
                "EGEIGVSPRVKDALSR")
            assert alignments[0][1].id == "gi|152973457|ref|YP_001338508.1|"
            assert alignments[0][1].annotations["original_length"] == 931
            assert len(alignments[1]) == 2
            assert alignments[1].get_alignment_length() == 64
            assert alignments[1][0].seq == "AAECGKTVSGFLRAAALGKKVNSLTDDRVLKEV-MRLGALQKKLFIDGKRVGDREYAEVLIAIT"
            assert alignments[1][0].id == "gi|10955263|ref|NP_052604.1|"
            assert alignments[1][0].annotations["original_length"] == 107
            assert alignments[1][1].seq == "ASRQGCTVGG--KMDSVQDKASDKDKERVMKNINIMWNALSKNRLFDG----NKELKEFIMTLT"
            assert alignments[1][1].id == "gi|152973588|ref|YP_001338639.1|"
            assert alignments[1][1].annotations["original_length"] == 459
            assert len(alignments[2]) == 2
            assert alignments[2].get_alignment_length() == 45
            assert alignments[2][0].seq == "EIRKRAAECGKTVSGFLRAAA-----LGKKVNSLTDDRVLKEVMR"
            assert alignments[2][0].id == "gi|10955263|ref|NP_052604.1|"
            assert alignments[2][0].annotations["original_length"] == 107
            assert alignments[2][1].seq == "ELVKLIADMGISVRALLRKNVEPYEELGLEEDKFTDDQLIDFMLQ"
            assert alignments[2][1].id == "gi|152973480|ref|YP_001338531.1|"
            assert alignments[2][1].annotations["original_length"] == 141
            assert len(alignments[3]) == 2
            assert alignments[3].get_alignment_length() == 38
            assert alignments[3][0].seq == "MKKDKKYQIEAIKNKDKTLFIVYATDIYSPSEFFSKIE"
            assert alignments[3][0].id == "gi|10955264|ref|NP_052605.1|"
            assert alignments[3][0].annotations["original_length"] == 126
            assert alignments[3][1].seq == "IKKDLGVSFLKLKNREKTLIVDALKKKYPVAELLSVLQ"
            assert alignments[3][1].id == "gi|152973462|ref|YP_001338513.1|"
            assert alignments[3][1].annotations["original_length"] == 101
            assert len(alignments[4]) == 2
            assert alignments[4].get_alignment_length() == 11
            assert alignments[4][0].seq == "FFDLIIPNGGK"
            assert alignments[4][0].id == "gi|10955264|ref|NP_052605.1|"
            assert alignments[4][0].annotations["original_length"] == 126
            assert alignments[4][1].seq == "FFDLVIENPGK"
            assert alignments[4][1].id == "gi|152973509|ref|YP_001338560.1|"
            assert alignments[4][1].annotations["original_length"] == 448
            assert len(alignments[5]) == 2
            assert alignments[5].get_alignment_length() == 40
            assert alignments[5][0].seq == "DKTLFIVYATDIYSPSE-FFSKIESDLKKKKSKGD-VFFD"
            assert alignments[5][0].id == "gi|10955264|ref|NP_052605.1|"
            assert alignments[5][0].annotations["original_length"] == 126
            assert alignments[5][1].seq == "ESVVFILMAGFAMSVCYLFFSVLEKVINARKSKDESIYHD"
            assert alignments[5][1].id == "gi|152973581|ref|YP_001338632.1|"
            assert alignments[5][1].annotations["original_length"] == 84
            assert len(alignments[6]) == 2
            assert alignments[6].get_alignment_length() == 30
            assert alignments[6][0].seq == "ASFFKKNFDKINVNLLSKATSFALKKGIPI"
            assert alignments[6][0].id == "gi|10955264|ref|NP_052605.1|"
            assert alignments[6][0].annotations["original_length"] == 126
            assert alignments[6][1].seq == "ASFSKEEQDKVAVDKVAADVAWQERMNKPV"
            assert alignments[6][1].id == "gi|152973536|ref|YP_001338587.1|"
            assert alignments[6][1].annotations["original_length"] == 84
            assert len(alignments[7]) == 2
            assert alignments[7].get_alignment_length() == 43
            assert alignments[7][0].seq == "SELHSKLPKSIDKIHEDIKKQLSC-SLIMKKIDVEMEDYSTYC"
            assert alignments[7][0].id == "gi|10955265|ref|NP_052606.1|"
            assert alignments[7][0].annotations["original_length"] == 346
            assert alignments[7][1].seq == "SRINSDVARRIPGIHRDPKDRLSSLKQVEEALDMLISSHGEYC"
            assert alignments[7][1].id == "gi|152973545|ref|YP_001338596.1|"
            assert alignments[7][1].annotations["original_length"] == 242
            assert len(alignments[8]) == 2
            assert alignments[8].get_alignment_length() == 64
            assert alignments[8][0].seq == "ISGTYKGIDFLIKLMPSGGNTTIGRASGQNNTYFDEIALIIKENCLY--SDTKNFEYTIPKFSD"
            assert alignments[8][0].id == "gi|10955265|ref|NP_052606.1|"
            assert alignments[8][0].annotations["original_length"] == 346
            assert alignments[8][1].seq == "IDGVITAFD-LRTGMNISKDKVVAQIQGMDPVW---ISAAVPESIAYLLKDTSQFEISVPAYPD"
            assert alignments[8][1].id == "gi|152973505|ref|YP_001338556.1|"
            assert alignments[8][1].annotations["original_length"] == 430

    def test_output008(self):
        """Check output008.m10 file."""
        fasta_file = "Fasta/output008.m10"
        with open(fasta_file) as handle:
            alignments = list(FastaIO.FastaM10Iterator(handle))
            assert len(alignments) == 12
            assert len(alignments[0]) == 2
            assert alignments[0].get_alignment_length() == 65
            assert alignments[0][0].seq == "LQHRHPHQQQQQQQQQQQQQQQQQQQQQQQQQQQH---HHHHHHHLLQDAYMQQYQHATQQQQML"
            assert alignments[0][0].id == "sp|Q9NSY1|BMP2K_HUMAN"
            assert alignments[0][0].annotations["original_length"] == 1161
            assert alignments[0][1].seq == "IPHQLPHALRHRPAQEAAHASQLHPAQPGCGQPLHGLWRLHHHPVYLYAWILRLRGHGMQSGGLL"
            assert alignments[0][1].id == "gi|283855822|gb|GQ290312.1|"
            assert alignments[0][1].annotations["original_length"] == 983
            assert len(alignments[1]) == 2
            assert alignments[1].get_alignment_length() == 201
            assert (alignments[1][0].seq == "GPEIL---LGQ-GPPQQPPQQHR"
                "VLQQLQQGDWRLQQLH-------"
                "LQHRHPHQQQQQQQQQQQQQQQQ"
                "QQQQQQQQQQQH-----HHHHHH"
                "-HLLQDAYMQQYQHATQQQQMLQ"
                "QQF-LMHSVYQPQPSASQYPTMM"
                "PQYQQAFFQQQMLAQHQPSQQQA"
                "SPEYLTSPQEFSPALVSYTSSLP"
                "A-QVGTIMDSSYSANRS")
            assert alignments[1][0].id == "sp|Q9NSY1|BMP2K_HUMAN"
            assert alignments[1][0].annotations["original_length"] == 1161
            assert (alignments[1][1].seq == "GPELLRALLQQNGCGTQPLRVPT"
                "VLPG*AMAVLHAGRLHVPAHRAW"
                "LPHQLPHALRHGPAQEAAHASQL"
                "HPAQPGRG*PLHGLRWLHHHPLH"
                "/PLCMDTLSLGPQDAIWRASLPH"
                "WAVKLPCGLWWSWPLSGTWWCVS"
                "P*ATSA------LGRTMP*WASL"
                "SPGSWHWPALHPPSLVGPGTSLK"
                "ACSVHAGSTTTHSSQKS")
            assert alignments[1][1].id == "gi|57163782|ref|NM_001009242.1|"
            assert alignments[1][1].annotations["original_length"] == 1047
            assert len(alignments[2]) == 2
            assert alignments[2].get_alignment_length() == 348
            assert (alignments[2][0].seq == "MNGTEGPNFYVPFSNATGVVRSP"
                "FEYPQYYLAEPWQFSMLAAYMFL"
                "LIVLGFPINFLTLYVTVQHKKLR"
                "TPLNYILLNLAVADLFMVLGGFT"
                "STLYTSLHGYFVFGPTGCNLEGF"
                "FATLGGEIALWSLVVLAIERYVV"
                "VCKPMSNFRFGENHAIMGVAFTW"
                "VMALACAAPPLAGWSRYIPEGLQ"
                "CSCGIDYYTLKPEVNNESFVIYM"
                "FVVHFTIPMIIIFFCYGQLVFTV"
                "KEAAAQQQESATTQKAEKEVTRM"
                "VIIMVIAFLICWVPYASVAFYIF"
                "THQGSNFGPIFMTIPAFFAKSAA"
                "IYNPVIYIMMNKQFRNCMLTTIC"
                "CGKNPLGDDEASATVSKTETSQV"
                "APA")
            assert alignments[2][0].id == "sp|P08100|OPSD_HUMAN"
            assert alignments[2][0].annotations["original_length"] == 348
            assert (alignments[2][1].seq == "MNGTEGPNFYVPFSNKTGVVRSP"
                "FEYPQYYLAEPWQFSMLAAYMFL"
                "LIVLGFPINFLTLYVTVQHKKLR"
                "TPLNYILLNLAVADLFMVFGGFT"
                "TTLYTSLHGYFVFGPTGCNLEGF"
                "FATLGGEIALWSLVVLAIERYVV"
                "VCKPMSNFRFGENHAIMGVAFTW"
                "VMALACAAPPLVGWSRYIPEGMQ"
                "CSCGIDYYTLKPEVNNESFVIYM"
                "FVVHFTIPMIVIFFCYGQLVFTV"
                "KEAAAQQQESATTQKAEKEVTRM"
                "VIIMVIAFLICWVPYASVAFYIF"
                "THQGSNFGPIFMTLPAFFAKSSS"
                "IYNPVIYIMMNKQFRNCMLTTLC"
                "CGKNPLGDDEASTTGSKTETSQV"
                "APA")
            assert alignments[2][1].id == "gi|57163782|ref|NM_001009242.1|"
            assert alignments[2][1].annotations["original_length"] == 1047
            assert len(alignments[3]) == 2
            assert alignments[3].get_alignment_length() == 348
            assert (alignments[3][0].seq == "MNGTEGPNFYVPFSNATGVVRSP"
                "FEYPQYYLAEPWQFSMLAAYMFL"
                "LIVLGFPINFLTLYVTVQHKKLR"
                "TPLNYILLNLAVADLFMVLGGFT"
                "STLYTSLHGYFVFGPTGCNLEGF"
                "FATLGGEIALWSLVVLAIERYVV"
                "VCKPMSNFRFGENHAIMGVAFTW"
                "VMALACAAPPLAGWSRYIPEGLQ"
                "CSCGIDYYTLKPEVNNESFVIYM"
                "FVVHFTIPMIIIFFCYGQLVFTV"
                "KEAAAQQQESATTQKAEKEVTRM"
                "VIIMVIAFLICWVPYASVAFYIF"
                "THQGSNFGPIFMTIPAFFAKSAA"
                "IYNPVIYIMMNKQFRNCMLTTIC"
                "CGKNPLGDDEASATVSKTETSQV"
                "APA")
            assert alignments[3][0].id == "sp|P08100|OPSD_HUMAN"
            assert alignments[3][0].annotations["original_length"] == 348
            assert (alignments[3][1].seq == "MNGTEGPNFYVPFSNKTGVVRSP"
                "FEAPQYYLAEPWQFSMLAAYMFL"
                "LIMLGFPINFLTLYVTVQHKKLR"
                "TPLNYILLNLAVADLFMVFGGFT"
                "TTLYTSLHGYFVFGPTGCNLEGF"
                "FATLGGEIALWSLVVLAIERYVV"
                "VCKPMSNFRFGENHAIMGVAFTW"
                "VMALACAAPPLVGWSRYIPEGMQ"
                "CSCGIDYYTPHEETNNESFVIYM"
                "FVVHFIIPLIVIFFCYGQLVFTV"
                "KEAAAQQQESATTQKAEKEVTRM"
                "VIIMVIAFLICWLPYAGVAFYIF"
                "THQGSDFGPIFMTIPAFFAKTSA"
                "VYNPVIYIMMNKQFRNCMVTTLC"
                "CGKNPLGDDEASTTVSKTETSQV"
                "APA")
            assert alignments[3][1].id == "gi|18148870|dbj|AB062417.1|"
            assert alignments[3][1].annotations["original_length"] == 1047
            assert len(alignments[4]) == 2
            assert alignments[4].get_alignment_length() == 326
            assert (alignments[4][0].seq == "VPFSNATGVVRSPFEYPQYYLAE"
                "PWQFSMLAAYMFLLIVLGFPINF"
                "LTLYVTVQHKKLRTPLNYILLNL"
                "AVADLFMVLGGFTSTLYTSLHGY"
                "FVFGPTGCNLEGFFATLGGEIAL"
                "WSLVVLAIERYVVVCKPMSNFRF"
                "GENHAIMGVAFTWVMALACAAPP"
                "LAGWSRYIPEGLQCSCGIDYYTL"
                "KPEVNNESFVIYMFVVHFTIPMI"
                "IIFFCYGQLVFTVKEAAAQQQES"
                "ATTQKAEKEVTRMVIIMVIAFLI"
                "CWVPYASVAFYIFTHQGSNFGPI"
                "FMTIPAFFAKSAAIYNPVIYIMM"
                "NKQFRNCMLTTICCGKNPLGDDE"
                "ASAT")
            assert alignments[4][0].id == "sp|P08100|OPSD_HUMAN"
            assert alignments[4][0].annotations["original_length"] == 348
            assert (alignments[4][1].seq == "VPFSNKTGVVRSPFEYPQYYLAE"
                "PWQFSMLAAYMFLLIVLGFPINF"
                "LTLYVTVQHKKLRTPLNYILLNL"
                "AVANLFMVFGGFTTTLYTSMHGY"
                "FVFGATGCNLEGFFATLGGEIAL"
                "WSLVVLAIERYVVVCKPMSNFRF"
                "GENHAIMGLAFTWVMALACAAPP"
                "LAGWSRYIPEGMQCSCGIDYYTL"
                "KPEVNNESFVIYMFVVHFTIPMI"
                "VIFFCYGQLVFTVKEAAAQQQES"
                "ATTQKAEKEVTRMVIIMVVAFLI"
                "CWLPYASVAFYIFTHQGSNFGPV"
                "FMTIPAFFAKSSSIYNPVIYIMM"
                "NKQFRNCMLTTLCCGKNPLGDDE"
                "ASTT")
            assert alignments[4][1].id == "gi|283855822|gb|GQ290312.1|"
            assert alignments[4][1].annotations["original_length"] == 983
            assert len(alignments[5]) == 2
            assert alignments[5].get_alignment_length() == 354
            assert (alignments[5][0].seq == "MNGTEGPNFYVPFSNATGVVRSP"
                "FEYPQYYLAEPWQFSMLAAYMFL"
                "LIVLGFPINFLTLYVTVQHKKLR"
                "TPLNYILLNLAVADLFMVLGGFT"
                "STLYTSLHGYFVFGPTGCNLEGF"
                "FATLGGEIALWSLVVLAIERYVV"
                "VCKPMSNFRFGENHAIMGVAFTW"
                "VMALACAAPPLAGWSRYIPEGLQ"
                "CSCGIDYYTLKPEVNNESFVIYM"
                "FVVHFTIPMIIIFFCYGQLVFTV"
                "KEAAAQQQESATTQKAEKEVTRM"
                "VIIMVIAFLICWVPYASVAFYIF"
                "THQGSNFGPIFMTIPAFFAKSAA"
                "IYNPVIYIMMNKQFRNCMLTTIC"
                "CGKNPLGDDEAS-ATVSKTE---"
                "--TSQVAPA")
            assert alignments[5][0].id == "sp|P08100|OPSD_HUMAN"
            assert alignments[5][0].annotations["original_length"] == 348
            assert (alignments[5][1].seq == "MNGTEGPNFYIPMSNKTGVVRSP"
                "FEYPQYYLAEPWQYSILCAYMFL"
                "LILLGFPINFMTLYVTIQHKKLR"
                "TPLNYILLNLAFANHFMVLCGFT"
                "VTMYSSMNGYFILGATGCYVEGF"
                "FATLGGEIALWSLVVLAIERYVV"
                "VCKPMSNFRFSENHAVMGVAFTW"
                "IMALSCAVPPLLGWSRYIPEGMQ"
                "CSCGVDYYTLKPEVNNESFVIYM"
                "FVVHFTIPLIIIFFCYGRLVCTV"
                "KEAAAQQQESATTQKAEKEVTRM"
                "VIIMVVFFLICWVPYASVAFFIF"
                "SNQGSEFGPIFMTVPAFFAKSSS"
                "IYNPVIYIMLNKQFRNCMITTLC"
                "CGKNPFGEDDASSAATSKTEASS"
                "VSSSQVSPA")
            assert alignments[5][1].id == "gi|2734705|gb|U59921.1|BBU59921"
            assert alignments[5][1].annotations["original_length"] == 1574
            assert len(alignments[6]) == 2
            assert alignments[6].get_alignment_length() == 347
            assert (alignments[6][0].seq == "MNGTEGPNFYVPFSNATGVVRSP"
                "FEYPQYYLAEPWQFSMLAAYMFL"
                "LIVLGFPINFLTLYVTVQHKKLR"
                "TPLNYILLNLAVADLFMVLGGFT"
                "STLYTSLHGYFVFGPTGCNLEGF"
                "FATLGGEIALWSLVVLAIERYVV"
                "VCKPMSNFRFGENHAIMGVAFTW"
                "VMALACAAPPLAGWSRYIPEGLQ"
                "CSCGIDYYTLKPEVNNESFVIYM"
                "FVVHFTIPMIIIFFCYGQLVFTV"
                "KEAAAQQQESATTQKAEKEVTRM"
                "VIIMVIAFLICWVPYASVAFYIF"
                "THQGSNFGPIFMTIPAFFAKSAA"
                "IYNPVIYIMMNKQFRNCMLTTIC"
                "CGKNPLGD-DEASATVSKTETSQ"
                "VA")
            assert alignments[6][0].id == "sp|P08100|OPSD_HUMAN"
            assert alignments[6][0].annotations["original_length"] == 348
            assert (alignments[6][1].seq == "MNGTEGPNFYIPMSNATGVVRSP"
                "FEYPQYYLAEPWAFSALSAYMFF"
                "LIIAGFPINFLTLYVTIEHKKLR"
                "TPLNYILLNLAVADLFMVFGGFT"
                "TTMYTSMHGYFVFGPTGCNIEGF"
                "FATLGGEIALWCLVVLAIERWMV"
                "VCKPVTNFRFGESHAIMGVMVTW"
                "TMALACALPPLFGWSRYIPEGLQ"
                "CSCGIDYYTRAPGINNESFVIYM"
                "FTCHFSIPLAVISFCYGRLVCTV"
                "KEAAAQQQESETTQRAEREVTRM"
                "VVIMVISFLVCWVPYASVAWYIF"
                "THQGSTFGPIFMTIPSFFAKSSA"
                "LYNPMIYICMNKQFRHCMITTLC"
                "CGKNPFEEEDGASATSSKTEASS"
                "VS")
            assert alignments[6][1].id == "gi|12583664|dbj|AB043817.1|"
            assert alignments[6][1].annotations["original_length"] == 1344
            assert len(alignments[7]) == 2
            assert alignments[7].get_alignment_length() == 111
            assert (alignments[7][0].seq == "VPFSNATGVVRSPFEYPQYYLAE"
                "PWQFSMLAAYMFLLIVLGFPINF"
                "LTLYVTVQHKKLRTPLNYILLNL"
                "AVADLFMVLGGFTSTLYTSLHGY"
                "FVFGPTGCNLEGFFATLGG")
            assert alignments[7][0].id == "sp|P08100|OPSD_HUMAN"
            assert alignments[7][0].annotations["original_length"] == 348
            assert (alignments[7][1].seq == "VPFSNKTGVVRSPFEHPQYYLAE"
                "PWQFSMLAAYMFLLIVLGFPINF"
                "LTLYVTVQHKKLRTPLNYILLNL"
                "AVADLFMVFGGFTTTLYTSLHGY"
                "FVFGPTGCNLEGFFATLGG")
            assert alignments[7][1].id == "gi|283855845|gb|GQ290303.1|"
            assert alignments[7][1].annotations["original_length"] == 4301
            assert len(alignments[8]) == 2
            assert alignments[8].get_alignment_length() == 172
            assert (alignments[8][0].seq == "RYIPEGLQCSCGIDYYTLKPEVN"
                "NESFVIYMFVVHFTIPMIIIFFC"
                "YGQLVFTVKE-------------"
                "-----------------------"
                "AAAQQQESATTQKAEKEVTRMVI"
                "IMVIAFLICWVPYASVAFYIFTH"
                "QGSNFGPIFMTIPAFFAKSAAIY"
                "NPVIYIMMNKQ")
            assert alignments[8][0].id == "sp|P08100|OPSD_HUMAN"
            assert alignments[8][0].annotations["original_length"] == 348
            assert (alignments[8][1].seq == "RYIPEGMQCSCGIDYYTLKPEVN"
                "NESFVIYMFVVHFTIPMIVIFFC"
                "YGQLVFTVKEVRSCVGHWGHAH*"
                "VNGAQLHSQSCHSLDT*PCVPA\\"
                "AAAQQQESATTQKAEKEVTRMVI"
                "IMVIAFLICWLPYAGVAFYIFTH"
                "QGSNFGPIFMTLPAFFAKSSSIY"
                "NPVIYIMMNKQ")
            assert alignments[8][1].id == "gi|283855845|gb|GQ290303.1|"
            assert alignments[8][1].annotations["original_length"] == 4301
            assert len(alignments[9]) == 2
            assert alignments[9].get_alignment_length() == 73
            assert (alignments[9][0].seq == "LGGEIALWSLVVLAIERYVVVCK"
                "PMSNFRFGENHAIMGVAFTWVMA"
                "LACAAPPLAGWSR--YIPEGLQC"
                "SCGI")
            assert alignments[9][0].id == "sp|P08100|OPSD_HUMAN"
            assert alignments[9][0].annotations["original_length"] == 348
            assert (alignments[9][1].seq == "LAGEIALWSLVVLAIERYVVVCK"
                "PMSNFRFGENHAIMGLALTWVMA"
                "LACAAPPLVGWSR*WH*TEG-KC"
                "L*GL")
            assert alignments[9][1].id == "gi|283855845|gb|GQ290303.1|"
            assert alignments[9][1].annotations["original_length"] == 4301
            assert len(alignments[10]) == 2
            assert alignments[10].get_alignment_length() == 30
            assert alignments[10][0].seq == "IMMNKQFRNCMLTTICCGKNPLGDDEASAT"
            assert alignments[10][0].id == "sp|P08100|OPSD_HUMAN"
            assert alignments[10][0].annotations["original_length"] == 348
            assert alignments[10][1].seq == "MLLAFQFRNCMLTTLCCGKNPLGDDEASTT"
            assert alignments[10][1].id == "gi|283855845|gb|GQ290303.1|"
            assert alignments[10][1].annotations["original_length"] == 4301
            assert len(alignments[11]) == 2
            assert alignments[11].get_alignment_length() == 31
            assert alignments[11][0].seq == "AQQQESATTQKAEKEVTRMVIIMVIAFLICW"
            assert alignments[11][0].id == "sp|P08100|OPSD_HUMAN"
            assert alignments[11][0].annotations["original_length"] == 348
            assert alignments[11][1].seq == "SQQIRNATTMMMTMRVTSFSAFWVVADSCCW"
            assert alignments[11][1].id == "gi|283855822|gb|GQ290312.1|"
            assert alignments[11][1].annotations["original_length"] == 983

    def test_output009(self):
        """Check output009.m10 file."""
        fasta_file = "Fasta/output009.m10"
        with open(fasta_file) as handle:
            alignments = list(FastaIO.FastaM10Iterator(handle))
            assert len(alignments) == 7
            assert len(alignments[0]) == 2
            assert alignments[0].get_alignment_length() == 22
            assert alignments[0][0].seq == "TGATGTTCTGTTTCTAAAACAG"
            assert alignments[0][0].id == "gi|255708421:1-99"
            assert alignments[0][0].annotations["original_length"] == 99
            assert alignments[0][1].seq == "TGATTTTTTTTGTCTAAAACAG"
            assert alignments[0][1].id == "gi|23308614|ref|NM_152952.1|"
            assert alignments[0][1].annotations["original_length"] == 5188
            assert len(alignments[1]) == 2
            assert alignments[1].get_alignment_length() == 14
            assert alignments[1][0].seq == "AGAAGGAAAAAAAA"
            assert alignments[1][0].id == "gi|156718121:2361-2376"
            assert alignments[1][0].annotations["original_length"] == 16
            assert alignments[1][1].seq == "AGAACTAAAAAAAA"
            assert alignments[1][1].id == "gi|47271416|ref|NM_131257.2|"
            assert alignments[1][1].annotations["original_length"] == 597
            assert len(alignments[2]) == 2
            assert alignments[2].get_alignment_length() == 14
            assert alignments[2][0].seq == "AGAAGGAAAAAAAA"
            assert alignments[2][0].id == "gi|156718121:2361-2376"
            assert alignments[2][0].annotations["original_length"] == 16
            assert alignments[2][1].seq == "AGAAGGTATAAAAA"
            assert alignments[2][1].id == "gi|332859474|ref|XM_001156938.2|"
            assert alignments[2][1].annotations["original_length"] == 762
            assert len(alignments[3]) == 2
            assert alignments[3].get_alignment_length() == 14
            assert alignments[3][0].seq == "TTTTTTTCCTTCTT"
            assert alignments[3][0].id == "gi|156718121:2361-2376"
            assert alignments[3][0].annotations["original_length"] == 16
            assert alignments[3][1].seq == "TTTTTTTACATCTT"
            assert alignments[3][1].id == "gi|332211534|ref|XM_003254825.1|"
            assert alignments[3][1].annotations["original_length"] == 805
            assert len(alignments[4]) == 2
            assert alignments[4].get_alignment_length() == 14
            assert alignments[4][0].seq == "AAGAAGGAAAAAAA"
            assert alignments[4][0].id == "gi|156718121:2361-2376"
            assert alignments[4][0].annotations["original_length"] == 16
            assert alignments[4][1].seq == "AATAAGTAAAAAAA"
            assert alignments[4][1].id == "gi|23308614|ref|NM_152952.1|"
            assert alignments[4][1].annotations["original_length"] == 5188
            assert len(alignments[5]) == 2
            assert alignments[5].get_alignment_length() == 14
            assert alignments[5][0].seq == "TTTTTTTCCTTCTT"
            assert alignments[5][0].id == "gi|156718121:2361-2376"
            assert alignments[5][0].annotations["original_length"] == 16
            assert alignments[5][1].seq == "TTTTTTTACATCTT"
            assert alignments[5][1].id == "gi|297689475|ref|XM_002822130.1|"
            assert alignments[5][1].annotations["original_length"] == 1158
            assert len(alignments[6]) == 2
            assert alignments[6].get_alignment_length() == 14
            assert alignments[6][0].seq == "AAGAAGGAAAAAAA"
            assert alignments[6][0].id == "gi|156718121:2361-2376"
            assert alignments[6][0].annotations["original_length"] == 16
            assert alignments[6][1].seq == "AAGAAGGTAAAAGA"
            assert alignments[6][1].id == "gi|297689475|ref|XM_002822130.1|"
            assert alignments[6][1].annotations["original_length"] == 1158


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
