# Copyright 2005 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Test the NCBI XML parser."""

import io
import os
import unittest
import pytest
import warnings

from Bio import BiopythonParserWarning
from Bio.Blast import NCBIXML


class TestNCBIXML(unittest.TestCase):
    """Tests for the NCBI XML parser."""

    def test_xml_2212L_blastp_001(self):
        """Parsing BLASTP 2.2.12, gi|49176427|ref|NP_418280.3| (xml_2212L_blastp_001)."""
        filename = "xml_2212L_blastp_001.xml"
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)
            record = next(records)
            with pytest.raises(StopIteration):
                next(records)
        self.check_xml_2212L_blastp_001(record)

        with open(datafile, "rb") as handle:
            record = NCBIXML.read(handle)
            self.check_xml_2212L_blastp_001(record)

    def check_xml_2212L_blastp_001(self, record):
        assert record.reference == 'Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schäffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.'
        assert record.date == "Aug-07-2005"
        assert record.version == "2.2.12"
        assert record.database == "nr"
        assert record.application == "BLASTP"
        alignments = record.alignments
        assert len(alignments) == 212
        assert record.query_id == "gi|49176427|ref|NP_418280.3|"
        assert sum(len(a.hsps) for a in alignments) == 212

        alignment = alignments[0]
        assert alignment.title[:50] == "gi|49176427|ref|NP_418280.3| component of Sec-inde"
        assert alignment.length == 103
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.expect, 4.20576e-46
        assert hsp.query[:75] == "MRLCLIIIYHRGTCMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTA"
        assert hsp.match[:75] == "MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTA"
        assert hsp.sbjct[:75] == "MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTA"

        alignment = alignments[1]
        assert alignment.title[:50] == "gi|15804428|ref|NP_290468.1| twin arginine translo"
        assert alignment.length == 103
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.expect, 2.72609e-45
        assert hsp.query[:75] == "MRLCLIIIYHRGTCMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTA"
        assert hsp.match[:75] == "MRLCLIIIYHR TCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTA"
        assert hsp.sbjct[:75] == "MRLCLIIIYHRXTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTA"

        alignment = alignments[2]
        assert alignment.title[:50] == "gi|74314349|ref|YP_312768.1| hypothetical protein "
        assert alignment.length == 103
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.expect, 2.72609e-45
        assert hsp.query[:75] == "MRLCLIIIYHRGTCMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTA"
        assert hsp.match[:75] == "MR CLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTA"
        assert hsp.sbjct[:75] == "MRPCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTA"

        alignment = alignments[3]
        assert alignment.title[:50] == "gi|75256240|ref|ZP_00727918.1| COG1826: Sec-indepe"
        assert alignment.length == 89
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.expect, 6.0872e-37
        assert hsp.query[:75] == "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ"
        assert hsp.match[:75] == "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ"
        assert hsp.sbjct[:75] == "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ"

        alignment = alignments[4]
        assert alignment.title[:50] == "gi|148236|gb|AAA67633.1| o261 [Escherichia coli]"
        assert alignment.length == 261
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.expect, 6.74582e-28
        assert hsp.query[:75] == "FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDA"
        assert hsp.match[:75] == "FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDA"
        assert hsp.sbjct[:75] == "FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDA"

        alignment = alignments[5]
        assert alignment.title[:50] == "gi|29143650|ref|NP_806992.1| sec-independent prote"
        assert alignment.length == 84
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.expect, 4.37251e-27
        assert hsp.query[:75] == "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ"
        assert hsp.match[:75] == "MGGISIWQLLI+AVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDD+ KQDKTSQDADFTAK+IADKQ      +"
        assert hsp.sbjct[:75] == "MGGISIWQLLIVAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDDAKQDKTSQDADFTAKSIADKQG-----E"

        alignment = alignments[6]
        assert alignment.title[:50] == "gi|49609685|emb|CAG73118.1| sec-independent protei"
        assert alignment.length == 86
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.expect, 4.10205e-17
        assert hsp.query[:75] == "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEP--KQDKTSQDADFTAKTIADKQADTNQ"
        assert hsp.match[:75] == "MGGIS+W LLIIAVIV+LLFGT KL ++GSDLGASIKGFKKAM DD+P    DK   DADF+ K+IAD Q+D   "
        assert hsp.sbjct[:75] == "MGGISLWNLLIIAVIVILLFGTNKLRTLGSDLGASIKGFKKAMGDDQPSTNADKAQPDADFSTKSIADNQSD---"

        alignment = alignments[7]
        assert alignment.title[:50] == "gi|37528238|ref|NP_931583.1| Sec-independent prote"
        assert alignment.length == 86
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.expect, 2.25087e-15
        assert hsp.query[:75] == "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDD-EPKQ-DKTSQDADFTAKTIADKQADTNQ"
        assert hsp.match[:75] == "MGGISIWQLLIIAVIVVLLFGT KL ++GSDLGASIKGFKKA+ DD +P+Q  KTS DADF  K I +KQ+    "
        assert hsp.sbjct[:75] == "MGGISIWQLLIIAVIVVLLFGTNKLRTLGSDLGASIKGFKKAIGDDNQPQQAQKTSSDADFETKNITEKQS----"

        alignment = alignments[8]
        assert alignment.title[:50] == "gi|59710656|ref|YP_203432.1| Sec-independent prote"
        assert alignment.length == 82
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.expect, 5.01441e-15
        assert hsp.query[:75] == "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ"
        assert hsp.match[:75] == "MGGISIWQLLIIAVI+VLLFGTKKL  +GSDLG+++KGFKKA+S+DEP ++   +DADF  + +  K+A+T ++Q"
        assert hsp.sbjct[:75] == "MGGISIWQLLIIAVIIVLLFGTKKLRGVGSDLGSAVKGFKKAISEDEPAKE-AKKDADFVPQNLEKKEAETVEKQ"

        alignment = alignments[9]
        assert alignment.title[:50] == "gi|54307340|ref|YP_128360.1| putative TatA protein"
        assert alignment.length == 87
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.expect, 7.2408e-14
        assert hsp.query[:75] == "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE--PKQDKTSQDADFTAKTIADKQADTNQ"
        assert hsp.match[:75] == "MGGISIWQLLIIA+I+VLLFGTKKL S+G DLG+++KGFKKA+ D+E   K+D T  DADF  KT++ ++  +  "
        assert hsp.sbjct[:75] == "MGGISIWQLLIIALIIVLLFGTKKLRSLGGDLGSAVKGFKKAIGDEELTVKKDNTEADADFEQKTLSKEEQQSED"

        alignment = alignments[10]
        assert alignment.title[:50] == "gi|45437890|gb|AAS63439.1| Sec-independent protein"
        assert alignment.length == 88
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.expect, 1.61308e-13
        assert hsp.query[:75] == "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE----PKQDKTSQDADFTAKTIADKQADT"
        assert hsp.match[:75] == "MG I   QLLIIAVIVVLLFGT KL ++GSDLGASIKGFKKAM DD        DKTS DADF AK+I +KQ   "
        assert hsp.sbjct[:75] == "MGSIGWAQLLIIAVIVVLLFGTNKLRTLGSDLGASIKGFKKAMGDDSQTPPTNVDKTSNDADF-AKSITEKQ---"

        alignment = alignments[11]
        assert alignment.title[:50] == "gi|75856473|ref|ZP_00764101.1| COG1826: Sec-indepe"
        assert alignment.length == 81
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.expect, 2.10675e-13
        assert hsp.query[:75] == "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ"
        assert hsp.match[:75] == "MGGIS+WQLLIIAVIVVLLFGTKKL  IG DLG ++KGFKKAMS+DEP   K  +DADF  K++ ++Q    +++"
        assert hsp.sbjct[:75] == "MGGISVWQLLIIAVIVVLLFGTKKLRGIGGDLGGAVKGFKKAMSEDEPA--KNDKDADFEPKSLEEQQ----KKE"

        alignment = alignments[12]
        assert alignment.title[:50] == "gi|75829371|ref|ZP_00758676.1| COG1826: Sec-indepe"
        assert alignment.length == 82
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.expect, 2.7515e-13
        assert hsp.query[:75] == "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ"
        assert hsp.match[:75] == "MGGISIWQLLIIAVIVVLLFGTKKL  IGSDLG+++KGFKKAMS++E       +DADF  K         N EQ"
        assert hsp.sbjct[:75] == "MGGISIWQLLIIAVIVVLLFGTKKLRGIGSDLGSAVKGFKKAMSEEESNSAANQKDADFETK---------NLEQ"

        alignment = alignments[13]
        assert alignment.title[:50] == "gi|75820019|ref|ZP_00750077.1| COG1826: Sec-indepe"
        assert alignment.length == 82
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.expect, 3.59357e-13
        assert hsp.query[:75] == "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ"
        assert hsp.match[:75] == "MGGISIWQLLIIAVIVVLLFGTKKL  IGSDLG+++KGFKKAMS++E       +DADF  K         N EQ"
        assert hsp.sbjct[:75] == "MGGISIWQLLIIAVIVVLLFGTKKLRGIGSDLGSAVKGFKKAMSEEESNSAANQKDADFETK---------NLEQ"

        alignment = alignments[14]
        assert alignment.title[:50] == "gi|28896872|ref|NP_796477.1| TatA protein [Vibrio "
        assert alignment.length == 81
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.expect, 2.32928e-12
        assert hsp.query[:75] == "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ"
        assert hsp.match[:75] == "MGGIS+WQLLIIAVIVVLLFGTKKL  IG DLG+++KGFKKAMSD++    K  +DADF  K++  +Q    Q++"
        assert hsp.sbjct[:75] == "MGGISVWQLLIIAVIVVLLFGTKKLRGIGGDLGSAVKGFKKAMSDED--SAKNEKDADFEPKSLEKQQ----QKE"

        alignment = alignments[15]
        assert alignment.title[:50] == "gi|27364353|ref|NP_759881.1| Sec-independent prote"
        assert alignment.length == 78
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.expect, 3.97316e-12
        assert hsp.query[:75] == "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ"
        assert hsp.match[:75] == "MGGISIWQLLIIAVIVVLLFGTKKL  IGSDLG +IKGFKKAM+++E ++    +DADF  K++     +   +Q"
        assert hsp.sbjct[:75] == "MGGISIWQLLIIAVIVVLLFGTKKLRGIGSDLGGAIKGFKKAMNEEESEK----KDADFEPKSL-----EQQSKQ"

        alignment = alignments[16]
        assert alignment.title[:50] == "gi|37678364|ref|NP_932973.1| Sec-independent prote"
        assert alignment.length == 78
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.expect, 3.97316e-12
        assert hsp.query[:75] == "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ"
        assert hsp.match[:75] == "MGGISIWQLLIIAVIVVLLFGTKKL  IGSDLG +IKGFKKAM+++E ++    +DADF  K++     +   +Q"
        assert hsp.sbjct[:75] == "MGGISIWQLLIIAVIVVLLFGTKKLRGIGSDLGGAIKGFKKAMNEEESEK----KDADFEPKSL-----EQQNKQ"

        alignment = alignments[17]
        assert alignment.title[:50] == "gi|71277787|ref|YP_266931.1| Sec-independent prote"
        assert alignment.length == 85
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.expect, 3.97316e-12
        assert hsp.query[:75] == "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ"
        assert hsp.match[:75] == "MGGI IWQL+I+AVIVVLLFGTKKL +IG DLG++IKGFK A+ +D  K+ K S  A+ T+ T+AD    T +E "
        assert hsp.sbjct[:75] == "MGGIGIWQLVIVAVIVVLLFGTKKLRNIGGDLGSAIKGFKSAIGED--KEQKNS--AEKTSDTLADSSKSTTEEV"

        alignment = alignments[18]
        assert alignment.title[:50] == "gi|68541995|ref|ZP_00581733.1| Twin-arginine trans"
        assert alignment.length == 79
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.expect, 2.57533e-11
        assert hsp.query[:75] == "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKT---IADKQADTN"
        assert hsp.match[:75] == "MGGISIWQLLI+A+IVVLLFGTKKL S+G DLG ++KGFK AMS +E K+     +A  TA+T     +K+ ++N"
        assert hsp.sbjct[:75] == "MGGISIWQLLIVALIVVLLFGTKKLRSLGGDLGGAVKGFKNAMSSEEDKKALEDTEAAKTAQTTQQATEKKPESN"

        alignment = alignments[19]
        assert alignment.title[:50] == "gi|77813363|ref|ZP_00812641.1| Twin-arginine trans"
        assert alignment.length == 79
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.expect, 2.57533e-11
        assert hsp.query[:75] == "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKT---IADKQADTN"
        assert hsp.match[:75] == "MGGISIWQLLIIA+IVVLLFGTKKL S+G DLG ++KGFK AMS +E K+     +A  TA+T     +K+ ++N"
        assert hsp.sbjct[:75] == "MGGISIWQLLIIALIVVLLFGTKKLRSLGGDLGGAVKGFKNAMSSEEDKKALEDTEAAKTAQTTQQATEKKPESN"

        alignment = alignments[20]
        assert alignment.title[:50] == "gi|52306607|gb|AAU37107.1| TatA protein [Mannheimi"
        assert alignment.length == 75
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.expect, 3.36348e-11
        assert hsp.query[:75] == "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ"
        assert hsp.match[:75] == "MGGISIWQLLII  I+VLLFGTKKL ++G+DLG S+KGFKKAM++DEPK      DA+F +    D+ A    E+"
        assert hsp.sbjct[:75] == "MGGISIWQLLIIVAIIVLLFGTKKLRTLGTDLGESVKGFKKAMNEDEPK------DAEFKSLN-KDESATAGSEK"

        alignment = alignments[21]
        assert alignment.title[:50] == "gi|75429751|ref|ZP_00732413.1| sec-independent pro"
        assert alignment.length == 74
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.expect, 3.36348e-11
        assert hsp.query[:75] == "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ"
        assert hsp.match[:75] == "MGGISIWQLLII  IVVLLFGTKKL ++GSDLG S+KGFKKAM+ +EPK      DA+F +   A+  A T +E+"
        assert hsp.sbjct[:75] == "MGGISIWQLLIIVAIVVLLFGTKKLRTLGSDLGESVKGFKKAMA-EEPK------DAEFKSLDKAENTAQTKKEE"

        alignment = alignments[22]
        assert alignment.title[:50] == "gi|32033565|ref|ZP_00133892.1| COG1826: Sec-indepe"
        assert alignment.length == 76
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.expect, 7.49305e-11
        assert hsp.query[:75] == "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ"
        assert hsp.match[:75] == "MGGISIWQLLII  I+VLLFGTKKL ++G+DLG S+KGFKKAM+DD+      SQ  D + + +  K+A + +++"
        assert hsp.sbjct[:75] == "MGGISIWQLLIIVAIIVLLFGTKKLRTLGTDLGESVKGFKKAMADDK------SQPQDASFEKVEAKEAASTEQK"

        alignment = alignments[23]
        assert alignment.title[:50] == "gi|12722097|gb|AAK03773.1| unknown [Pasteurella mu"
        assert alignment.length == 76
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(1.66928e-10, abs=5e-16)

        alignment = alignments[24]
        assert alignment.title[:50] == "gi|68546478|ref|ZP_00586025.1| Twin-arginine trans"
        assert alignment.length == 77
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(2.18014e-10, abs=5e-16)

        alignment = alignments[25]
        assert alignment.title[:50] == "gi|33151888|ref|NP_873241.1| sec-independent prote"
        assert alignment.length == 74
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(3.71876e-10, abs=5e-16)

        alignment = alignments[26]
        assert alignment.title[:50] == "gi|24375687|ref|NP_719730.1| Sec-independent prote"
        assert alignment.length == 88
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(3.71876e-10, abs=5e-16)

        alignment = alignments[27]
        assert alignment.title[:50] == "gi|71278553|ref|YP_269744.1| Sec-independent prote"
        assert alignment.length == 79
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(4.85685e-10, abs=5e-16)

        alignment = alignments[28]
        assert alignment.title[:50] == "gi|69159855|gb|EAN71956.1| Twin-arginine transloca"
        assert alignment.length == 79
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(6.34325e-10, abs=5e-16)

        alignment = alignments[29]
        assert alignment.title[:50] == "gi|69949858|ref|ZP_00637822.1| Twin-arginine trans"
        assert alignment.length == 81
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(1.08199e-09, abs=5e-15)

        alignment = alignments[30]
        assert alignment.title[:50] == "gi|48863844|ref|ZP_00317737.1| hypothetical protei"
        assert alignment.length == 83
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(1.41313e-09, abs=5e-15)

        alignment = alignments[31]
        assert alignment.title[:50] == "gi|77361831|ref|YP_341406.1| twin-arginine translo"
        assert alignment.length == 82
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(1.5624e-08, abs=5e-13)

        alignment = alignments[32]
        assert alignment.title[:50] == "gi|67676224|ref|ZP_00472975.1| Twin-arginine trans"
        assert alignment.length == 90
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(2.04055e-08, abs=5e-14)

        alignment = alignments[33]
        assert alignment.title[:50] == "gi|74317722|ref|YP_315462.1| twin-arginine translo"
        assert alignment.length == 70
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(3.48066e-08, abs=5e-14)

        alignment = alignments[34]
        assert alignment.title[:50] == "gi|77166504|ref|YP_345029.1| Twin-arginine translo"
        assert alignment.length == 90
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(5.9371e-08, abs=5e-13)

        alignment = alignments[35]
        assert alignment.title[:50] == "gi|16128610|ref|NP_415160.1| component of Sec-inde"
        assert alignment.length == 67
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(1.32265e-07, abs=5e-13)

        alignment = alignments[36]
        assert alignment.title[:50] == "gi|12831974|emb|CAC29147.1| TatA protein [Pseudomo"
        assert alignment.length == 76
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(1.72743e-07, abs=5e-13)

        alignment = alignments[37]
        assert alignment.title[:50] == "gi|32029972|ref|ZP_00132908.1| COG1826: Sec-indepe"
        assert alignment.length == 73
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(2.2561e-07, abs=5e-12)

        alignment = alignments[38]
        assert alignment.title[:50] == "gi|455172|gb|AAA24073.1| ORF; putative"
        assert alignment.length == 67
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(2.94655e-07, abs=5e-13)

        alignment = alignments[39]
        assert alignment.title[:50] == "gi|1224007|gb|AAA92108.1| ORF4"
        assert alignment.length == 192
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(3.84832e-07, abs=5e-13)

        alignment = alignments[40]
        assert alignment.title[:50] == "gi|68056990|gb|AAX87243.1| Sec-independent protein"
        assert alignment.length == 95
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(6.56423e-07, abs=5e-13)

        alignment = alignments[41]
        assert alignment.title[:50] == "gi|56461470|ref|YP_156751.1| Sec-independent prote"
        assert alignment.length == 73
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(1.11969e-06, abs=5e-12)

        alignment = alignments[42]
        assert alignment.title[:50] == "gi|76793313|ref|ZP_00775802.1| Twin-arginine trans"
        assert alignment.length == 84
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(1.46236e-06, abs=5e-12)

        alignment = alignments[43]
        assert alignment.title[:50] == "gi|42630489|ref|ZP_00156028.1| COG1826: Sec-indepe"
        assert alignment.length == 75
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(1.9099e-06, abs=5e-11)

        alignment = alignments[44]
        assert alignment.title[:50] == "gi|1074302|pir||B64145 hypothetical protein HI0187"
        assert alignment.length == 109
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(2.49441e-06, abs=5e-12)

        alignment = alignments[45]
        assert alignment.title[:50] == "gi|67641583|ref|ZP_00440359.1| COG1826: Sec-indepe"
        assert alignment.length == 77
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(3.25779e-06, abs=5e-12)

        alignment = alignments[46]
        assert alignment.title[:50] == "gi|67545726|ref|ZP_00423646.1| Twin-arginine trans"
        assert alignment.length == 76
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(3.25779e-06, abs=5e-12)

        alignment = alignments[47]
        assert alignment.title[:50] == "gi|45435806|gb|AAS61363.1| sec-independent protein"
        assert alignment.length == 85
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(3.25779e-06, abs=5e-12)

        alignment = alignments[48]
        assert alignment.title[:50] == "gi|49610761|emb|CAG74206.1| Sec-independent protei"
        assert alignment.length == 65
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(3.25779e-06, abs=5e-12)

        alignment = alignments[49]
        assert alignment.title[:50] == "gi|67663266|ref|ZP_00460549.1| Twin-arginine trans"
        assert alignment.length == 76
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(4.25481e-06, abs=5e-12)

        alignment = alignments[50]
        assert alignment.title[:50] == "gi|33594634|ref|NP_882278.1| Sec-independent prote"
        assert alignment.length == 75
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(4.25481e-06, abs=5e-12)

        alignment = alignments[51]
        assert alignment.title[:50] == "gi|46310681|ref|ZP_00211309.1| COG1826: Sec-indepe"
        assert alignment.length == 76
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(4.25481e-06, abs=5e-12)

        alignment = alignments[52]
        assert alignment.title[:50] == "gi|58584031|ref|YP_203047.1| sec-independent prote"
        assert alignment.length == 75
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(5.55696e-06, abs=5e-12)

        alignment = alignments[53]
        assert alignment.title[:50] == "gi|17429965|emb|CAD16649.1| PROBABLE SIGNAL PEPTID"
        assert alignment.length == 85
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(5.55696e-06, abs=5e-12)

        alignment = alignments[54]
        assert alignment.title[:50] == "gi|47573371|ref|ZP_00243410.1| COG1826: Sec-indepe"
        assert alignment.length == 76
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(5.55696e-06, abs=5e-12)

        alignment = alignments[55]
        assert alignment.title[:50] == "gi|16273687|ref|NP_438355.1| Sec-independent prote"
        assert alignment.length == 89
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(7.25761e-06, abs=5e-12)

        alignment = alignments[56]
        assert alignment.title[:50] == "gi|73542784|ref|YP_297304.1| Twin-arginine translo"
        assert alignment.length == 73
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(9.47873e-06, abs=5e-12)

        alignment = alignments[57]
        assert alignment.title[:50] == "gi|26987777|ref|NP_743202.1| Sec-independent prote"
        assert alignment.length == 77
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(9.47873e-06, abs=5e-12)

        alignment = alignments[58]
        assert alignment.title[:50] == "gi|29142636|ref|NP_805978.1| sec-independent prote"
        assert alignment.length == 67
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(1.61683e-05, abs=5e-11)

        alignment = alignments[59]
        assert alignment.title[:50] == "gi|18389921|gb|AAL68797.1| TatA [Ralstonia eutroph"
        assert alignment.length == 77
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(1.61683e-05, abs=5e-11)

        alignment = alignments[60]
        assert alignment.title[:50] == "gi|48781637|ref|ZP_00278228.1| COG1826: Sec-indepe"
        assert alignment.length == 79
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(2.75789e-05, abs=5e-11)

        alignment = alignments[61]
        assert alignment.title[:50] == "gi|77456610|ref|YP_346115.1| Twin-arginine translo"
        assert alignment.length == 92
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(2.75789e-05, abs=5e-11)

        alignment = alignments[62]
        assert alignment.title[:50] == "gi|1684735|emb|CAA98158.1| ORF57 protein [Pseudomo"
        assert alignment.length == 57
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(2.75789e-05, abs=5e-11)

        alignment = alignments[63]
        assert alignment.title[:50] == "gi|56476124|ref|YP_157713.1| Sec-independent prote"
        assert alignment.length == 75
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(3.60191e-05, abs=5e-11)

        alignment = alignments[64]
        assert alignment.title[:50] == "gi|34496078|ref|NP_900293.1| Sec-independent prote"
        assert alignment.length == 68
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(3.60191e-05, abs=5e-11)

        alignment = alignments[65]
        assert alignment.title[:50] == "gi|67848115|ref|ZP_00503233.1| Twin-arginine trans"
        assert alignment.length == 83
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(4.70425e-05, abs=5e-11)

        alignment = alignments[66]
        assert alignment.title[:50] == "gi|26991692|ref|NP_747117.1| Sec-independent prote"
        assert alignment.length == 90
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(4.70425e-05, abs=5e-11)

        alignment = alignments[67]
        assert alignment.title[:50] == "gi|15601293|ref|NP_232924.1| tatA protein [Vibrio "
        assert alignment.length == 78
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(6.14393e-05, abs=5e-11)

        alignment = alignments[68]
        assert alignment.title[:50] == "gi|66770480|ref|YP_245242.1| sec-independent prote"
        assert alignment.length == 75
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(8.02423e-05, abs=5e-11)

        alignment = alignments[69]
        assert alignment.title[:50] == "gi|53804435|ref|YP_113945.1| Sec-independent prote"
        assert alignment.length == 70
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(8.02423e-05, abs=5e-11)

        alignment = alignments[70]
        assert alignment.title[:50] == "gi|75825357|ref|ZP_00754793.1| COG1826: Sec-indepe"
        assert alignment.length == 80
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(8.02423e-05, abs=5e-11)

        alignment = alignments[71]
        assert alignment.title[:50] == "gi|71908987|ref|YP_286574.1| Twin-arginine translo"
        assert alignment.length == 76
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(8.02423e-05, abs=5e-11)

        alignment = alignments[72]
        assert alignment.title[:50] == "gi|68526571|gb|EAN49542.1| Twin-arginine transloca"
        assert alignment.length == 77
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.0001048, abs=5e-8)

        alignment = alignments[73]
        assert alignment.title[:50] == "gi|71736448|ref|YP_272670.1| sec-independent prote"
        assert alignment.length == 91
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.0001048, abs=5e-8)

        alignment = alignments[74]
        assert alignment.title[:50] == "gi|56460344|ref|YP_155625.1| Sec-independent prote"
        assert alignment.length == 72
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.0001048, abs=5e-8)

        alignment = alignments[75]
        assert alignment.title[:50] == "gi|68214708|ref|ZP_00566522.1| Twin-arginine trans"
        assert alignment.length == 72
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.000136873, abs=5e-8)

        alignment = alignments[76]
        assert alignment.title[:50] == "gi|30248650|ref|NP_840720.1| mttA/Hcf106 family [N"
        assert alignment.length == 76
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.000136873, abs=5e-8)

        alignment = alignments[77]
        assert alignment.title[:50] == "gi|75822907|ref|ZP_00752458.1| COG1826: Sec-indepe"
        assert alignment.length == 78
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.000136873, abs=5e-8)

        alignment = alignments[78]
        assert alignment.title[:50] == "gi|70733926|ref|YP_257566.1| sec-independent prote"
        assert alignment.length == 93
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.000136873, abs=5e-8)

        alignment = alignments[79]
        assert alignment.title[:50] == "gi|63254358|gb|AAY35454.1| Twin-arginine transloca"
        assert alignment.length == 91
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.000178761, abs=5e-8)

        alignment = alignments[80]
        assert alignment.title[:50] == "gi|73354814|gb|AAZ75668.1| TatA [Pseudomonas syrin"
        assert alignment.length == 91
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.000178761, abs=5e-8)

        alignment = alignments[81]
        assert alignment.title[:50] == "gi|50083761|ref|YP_045271.1| Sec-independent prote"
        assert alignment.length == 78
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.000233469, abs=5e-8)

        alignment = alignments[82]
        assert alignment.title[:50] == "gi|71548504|ref|ZP_00668728.1| Twin-arginine trans"
        assert alignment.length == 77
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.00030492, abs=5e-8)

        alignment = alignments[83]
        assert alignment.title[:50] == "gi|55247002|gb|EAL42253.1| ENSANGP00000028218 [Ano"
        assert alignment.length == 53
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.000398238, abs=5e-8)

        alignment = alignments[84]
        assert alignment.title[:50] == "gi|50084688|ref|YP_046198.1| Sec-independent prote"
        assert alignment.length == 71
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.000520115, abs=5e-8)

        alignment = alignments[85]
        assert alignment.title[:50] == "gi|28872267|ref|NP_794886.1| sec-independent prote"
        assert alignment.length == 91
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.000679292, abs=5e-8)

        alignment = alignments[86]
        assert alignment.title[:50] == "gi|49082486|gb|AAT50643.1| PA5068 [synthetic const"
        assert alignment.length == 83
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.000679292, abs=5e-8)

        alignment = alignments[87]
        assert alignment.title[:50] == "gi|53726598|ref|ZP_00141543.2| COG1826: Sec-indepe"
        assert alignment.length == 82
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.000679292, abs=5e-8)

        alignment = alignments[88]
        assert alignment.title[:50] == "gi|68213616|ref|ZP_00565447.1| Twin-arginine trans"
        assert alignment.length == 54
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.000887182, abs=5e-8)

        alignment = alignments[89]
        assert alignment.title[:50] == "gi|74023810|ref|ZP_00694377.1| Twin-arginine trans"
        assert alignment.length == 79
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.0011587, abs=5e-8)

        alignment = alignments[90]
        assert alignment.title[:50] == "gi|71066554|ref|YP_265281.1| twin-arginine translo"
        assert alignment.length == 87
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.0015133, abs=5e-8)

        alignment = alignments[91]
        assert alignment.title[:50] == "gi|15611372|ref|NP_223023.1| hypothetical protein "
        assert alignment.length == 79
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.00197644, abs=5e-8)

        alignment = alignments[92]
        assert alignment.title[:50] == "gi|13471183|ref|NP_102752.1| sec-independent prote"
        assert alignment.length == 73
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.00337129, abs=5e-8)

        alignment = alignments[93]
        assert alignment.title[:50] == "gi|42523995|ref|NP_969375.1| twin-arginine-depende"
        assert alignment.length == 81
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.00337129, abs=5e-8)

        alignment = alignments[94]
        assert alignment.title[:50] == "gi|67158086|ref|ZP_00419176.1| Twin-arginine trans"
        assert alignment.length == 85
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.00440304, abs=5e-8)

        alignment = alignments[95]
        assert alignment.title[:50] == "gi|15644948|ref|NP_207118.1| conserved hypothetica"
        assert alignment.length == 79
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.00440304, abs=5e-8)

        alignment = alignments[96]
        assert alignment.title[:50] == "gi|13277311|emb|CAC34414.1| putative TatA protein "
        assert alignment.length == 61
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.00751045, abs=5e-8)

        alignment = alignments[97]
        assert alignment.title[:50] == "gi|54298906|ref|YP_125275.1| Putative TatA protein"
        assert alignment.length == 61
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.00751045, abs=5e-8)

        alignment = alignments[98]
        assert alignment.title[:50] == "gi|71363513|ref|ZP_00654157.1| Twin-arginine trans"
        assert alignment.length == 94
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.00751045, abs=5e-8)

        alignment = alignments[99]
        assert alignment.title[:50] == "gi|71362217|ref|ZP_00653377.1| Twin-arginine trans"
        assert alignment.length == 80
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.00751045, abs=5e-8)

        alignment = alignments[100]
        assert alignment.title[:50] == "gi|27379862|ref|NP_771391.1| hypothetical protein "
        assert alignment.length == 77
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.00980895, abs=5e-8)

        alignment = alignments[101]
        assert alignment.title[:50] == "gi|39935914|ref|NP_948190.1| putative sec-independ"
        assert alignment.length == 78
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.0128109, abs=5e-8)

        alignment = alignments[102]
        assert alignment.title[:50] == "gi|17935600|ref|NP_532390.1| SEC-independent prote"
        assert alignment.length == 70
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.0128109, abs=5e-8)

        alignment = alignments[103]
        assert alignment.title[:50] == "gi|62289827|ref|YP_221620.1| Sec-independent prote"
        assert alignment.length == 72
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.0167315, abs=5e-8)

        alignment = alignments[104]
        assert alignment.title[:50] == "gi|23347697|gb|AAN29810.1| Sec-independent protein"
        assert alignment.length == 80
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.0167315, abs=5e-8)

        alignment = alignments[105]
        assert alignment.title[:50] == "gi|75675971|ref|YP_318392.1| twin-arginine translo"
        assert alignment.length == 78
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.0167315, abs=5e-8)

        alignment = alignments[106]
        assert alignment.title[:50] == "gi|69928230|ref|ZP_00625391.1| Twin-arginine trans"
        assert alignment.length == 79
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.0167315, abs=5e-8)

        alignment = alignments[107]
        assert alignment.title[:50] == "gi|77689454|ref|ZP_00804635.1| Twin-arginine trans"
        assert alignment.length == 79
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.0218521, abs=5e-8)

        alignment = alignments[108]
        assert alignment.title[:50] == "gi|77743614|ref|ZP_00812071.1| Twin-arginine trans"
        assert alignment.length == 78
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.0285397, abs=5e-8)

        alignment = alignments[109]
        assert alignment.title[:50] == "gi|71066141|ref|YP_264868.1| twin-arginine translo"
        assert alignment.length == 89
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.037274, abs=5e-8)

        alignment = alignments[110]
        assert alignment.title[:50] == "gi|28199457|ref|NP_779771.1| SEC-independent prote"
        assert alignment.length == 71
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.0486813, abs=5e-8)

        alignment = alignments[111]
        assert alignment.title[:50] == "gi|15837166|ref|NP_297854.1| hypothetical protein "
        assert alignment.length == 71
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.0486813, abs=5e-8)

        alignment = alignments[112]
        assert alignment.title[:50] == "gi|15074462|emb|CAC46108.1| HYPOTHETICAL TRANSMEMB"
        assert alignment.length == 68
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.0486813, abs=5e-8)

        alignment = alignments[113]
        assert alignment.title[:50] == "gi|27462871|gb|AAO15625.1| Sec-independent protein"
        assert alignment.length == 63
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.0830377, abs=5e-8)

        alignment = alignments[114]
        assert alignment.title[:50] == "gi|35211273|dbj|BAC88652.1| gsl0711 [Gloeobacter v"
        assert alignment.length == 72
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.0830377, abs=5e-8)

        alignment = alignments[115]
        assert alignment.title[:50] == "gi|34482347|emb|CAE09348.1| hypothetical protein ["
        assert alignment.length == 80
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.108451, abs=5e-8)

        alignment = alignments[116]
        assert alignment.title[:50] == "gi|32262257|gb|AAP77305.1| component of Sec-indepe"
        assert alignment.length == 82
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.108451, abs=5e-8)

        alignment = alignments[117]
        assert alignment.title[:50] == "gi|76261408|ref|ZP_00769019.1| Twin-arginine trans"
        assert alignment.length == 62
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.141641, abs=5e-8)

        alignment = alignments[118]
        assert alignment.title[:50] == "gi|69933726|ref|ZP_00628928.1| sec-independent tra"
        assert alignment.length == 159
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.141641, abs=5e-8)

        alignment = alignments[119]
        assert alignment.title[:50] == "gi|15605662|ref|NP_213037.1| hypothetical protein "
        assert alignment.length == 59
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.184989, abs=5e-8)

        alignment = alignments[120]
        assert alignment.title[:50] == "gi|68538777|ref|ZP_00578553.1| Twin-arginine trans"
        assert alignment.length == 77
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.241603, abs=5e-8)

        alignment = alignments[121]
        assert alignment.title[:50] == "gi|68136098|ref|ZP_00544086.1| Twin-arginine trans"
        assert alignment.length == 130
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.315543, abs=5e-8)

        alignment = alignments[122]
        assert alignment.title[:50] == "gi|20259265|gb|AAM14368.1| putative Tha4 protein ["
        assert alignment.length == 147
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.412112, abs=5e-8)

        alignment = alignments[123]
        assert alignment.title[:50] == "gi|75910646|ref|YP_324942.1| Twin-arginine translo"
        assert alignment.length == 90
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.412112, abs=5e-8)

        alignment = alignments[124]
        assert alignment.title[:50] == "gi|39982657|gb|AAR34117.1| twin-arginine transloca"
        assert alignment.length == 57
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.412112, abs=5e-8)

        alignment = alignments[125]
        assert alignment.title[:50] == "gi|33635687|emb|CAE22011.1| mttA/Hcf106 family [Pr"
        assert alignment.length == 91
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.538235, abs=5e-8)

        alignment = alignments[126]
        assert alignment.title[:50] == "gi|76791934|ref|ZP_00774438.1| Twin-arginine trans"
        assert alignment.length == 68
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.538235, abs=5e-8)

        alignment = alignments[127]
        assert alignment.title[:50] == "gi|23129516|ref|ZP_00111343.1| COG1826: Sec-indepe"
        assert alignment.length == 91
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.538235, abs=5e-8)

        alignment = alignments[128]
        assert alignment.title[:50] == "gi|48764199|ref|ZP_00268751.1| COG1826: Sec-indepe"
        assert alignment.length == 96
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.702957, abs=5e-8)

        alignment = alignments[129]
        assert alignment.title[:50] == "gi|15677995|ref|NP_273645.1| hypothetical protein "
        assert alignment.length == 67
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.702957, abs=5e-8)

        alignment = alignments[130]
        assert alignment.title[:50] == "gi|50917153|ref|XP_468973.1| putative sec-independ"
        assert alignment.length == 170
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.702957, abs=5e-8)

        alignment = alignments[131]
        assert alignment.title[:50] == "gi|16329622|ref|NP_440350.1| hypothetical protein "
        assert alignment.length == 126
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.702957, abs=5e-8)

        alignment = alignments[132]
        assert alignment.title[:50] == "gi|71083667|ref|YP_266387.1| Twin-arginine translo"
        assert alignment.length == 66
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.702957, abs=5e-8)

        alignment = alignments[133]
        assert alignment.title[:50] == "gi|17130190|dbj|BAB72802.1| asl0845 [Nostoc sp. PC"
        assert alignment.length == 90
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.702957, abs=5e-8)

        alignment = alignments[134]
        assert alignment.title[:50] == "gi|68246031|gb|EAN28138.1| Twin-arginine transloca"
        assert alignment.length == 69
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.91809, abs=5e-8)

        alignment = alignments[135]
        assert alignment.title[:50] == "gi|15604583|ref|NP_221101.1| hypothetical protein "
        assert alignment.length == 54
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.91809, abs=5e-8)

        alignment = alignments[136]
        assert alignment.title[:50] == "gi|77685166|ref|ZP_00800574.1| Twin-arginine trans"
        assert alignment.length == 69
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(0.91809, abs=5e-8)

        alignment = alignments[137]
        assert alignment.title[:50] == "gi|39985226|gb|AAR36581.1| twin-arginine transloca"
        assert alignment.length == 78
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(1.19906, abs=5e-8)

        alignment = alignments[138]
        assert alignment.title[:50] == "gi|1825636|gb|AAB42258.1| Hypothetical protein ZK3"
        assert alignment.length == 312
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(1.56602, abs=5e-8)

        alignment = alignments[139]
        assert alignment.title[:50] == "gi|65321915|ref|ZP_00394874.1| COG5386: Cell surfa"
        assert alignment.length == 237
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(1.56602, abs=5e-8)

        alignment = alignments[140]
        assert alignment.title[:50] == "gi|30022625|ref|NP_834256.1| Cell surface protein "
        assert alignment.length == 237
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(1.56602, abs=5e-8)

        alignment = alignments[141]
        assert alignment.title[:50] == "gi|55623442|ref|XP_517520.1| PREDICTED: similar to"
        assert alignment.length == 234
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(1.56602, abs=5e-8)

        alignment = alignments[142]
        assert alignment.title[:50] == "gi|75762866|ref|ZP_00742681.1| Cell surface protei"
        assert alignment.length == 237
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(1.56602, abs=5e-8)

        alignment = alignments[143]
        assert alignment.title[:50] == "gi|22945598|gb|AAN10511.1| CG18497-PC, isoform C ["
        assert alignment.length == 5476
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(1.56602, abs=5e-8)

        alignment = alignments[144]
        assert alignment.title[:50] == "gi|10727420|gb|AAF51534.2| CG18497-PB, isoform B ["
        assert alignment.length == 5533
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(1.56602, abs=5e-8)

        alignment = alignments[145]
        assert alignment.title[:50] == "gi|10727421|gb|AAF51535.2| CG18497-PA, isoform A ["
        assert alignment.length == 5560
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(1.56602, abs=5e-8)

        alignment = alignments[146]
        assert alignment.title[:50] == "gi|71481981|ref|ZP_00661682.1| Twin-arginine trans"
        assert alignment.length == 69
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(1.56602, abs=5e-8)

        alignment = alignments[147]
        assert alignment.title[:50] == "gi|71150623|ref|ZP_00649545.1| Twin-arginine trans"
        assert alignment.length == 81
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(1.56602, abs=5e-8)

        alignment = alignments[148]
        assert alignment.title[:50] == "gi|20151563|gb|AAM11141.1| LD15253p [Drosophila me"
        assert alignment.length == 1521
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(1.56602, abs=5e-8)

        alignment = alignments[149]
        assert alignment.title[:50] == "gi|6979936|gb|AAF34661.1| split ends long isoform "
        assert alignment.length == 5554
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(1.56602, abs=5e-8)

        alignment = alignments[150]
        assert alignment.title[:50] == "gi|6467825|gb|AAF13218.1| Spen RNP motif protein l"
        assert alignment.length == 5533
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(1.56602, abs=5e-8)

        alignment = alignments[151]
        assert alignment.title[:50] == "gi|61102013|ref|ZP_00377467.1| hypothetical protei"
        assert alignment.length == 80
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(1.56602, abs=5e-8)

        alignment = alignments[152]
        assert alignment.title[:50] == "gi|68056232|ref|ZP_00540361.1| Twin-arginine trans"
        assert alignment.length == 68
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(2.04529, abs=5e-8)

        alignment = alignments[153]
        assert alignment.title[:50] == "gi|68190120|gb|EAN04781.1| Twin-arginine transloca"
        assert alignment.length == 71
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(2.04529, abs=5e-8)

        alignment = alignments[154]
        assert alignment.title[:50] == "gi|15605663|ref|NP_213038.1| hypothetical protein "
        assert alignment.length == 77
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(2.04529, abs=5e-8)

        alignment = alignments[155]
        assert alignment.title[:50] == "gi|60493413|emb|CAH08199.1| aerotolerance-related "
        assert alignment.length == 238
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(2.04529, abs=5e-8)

        alignment = alignments[156]
        assert alignment.title[:50] == "gi|50877510|emb|CAG37350.1| related to Sec-indepen"
        assert alignment.length == 84
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(2.04529, abs=5e-8)

        alignment = alignments[157]
        assert alignment.title[:50] == "gi|42739647|gb|AAS43573.1| conserved domain protei"
        assert alignment.length == 236
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(2.04529, abs=5e-8)

        alignment = alignments[158]
        assert alignment.title[:50] == "gi|53713708|ref|YP_099700.1| conserved hypothetica"
        assert alignment.length == 238
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(2.04529, abs=5e-8)

        alignment = alignments[159]
        assert alignment.title[:50] == "gi|33860901|ref|NP_892462.1| mttA/Hcf106 family [P"
        assert alignment.length == 96
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(2.67123, abs=5e-8)

        alignment = alignments[160]
        assert alignment.title[:50] == "gi|48851224|ref|ZP_00305466.1| COG1826: Sec-indepe"
        assert alignment.length == 83
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(2.67123, abs=5e-8)

        alignment = alignments[161]
        assert alignment.title[:50] == "gi|67938449|ref|ZP_00530974.1| Twin-arginine trans"
        assert alignment.length == 69
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(2.67123, abs=5e-8)

        alignment = alignments[162]
        assert alignment.title[:50] == "gi|45657833|ref|YP_001919.1| sec-independent prote"
        assert alignment.length == 90
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(2.67123, abs=5e-8)

        alignment = alignments[163]
        assert alignment.title[:50] == "gi|57238048|ref|YP_179297.1| twin-arginine translo"
        assert alignment.length == 79
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(2.67123, abs=5e-8)

        alignment = alignments[164]
        assert alignment.title[:50] == "gi|56962648|ref|YP_174374.1| sec-independent prote"
        assert alignment.length == 63
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(2.67123, abs=5e-8)

        alignment = alignments[165]
        assert alignment.title[:50] == "gi|33239734|ref|NP_874676.1| Sec-independent prote"
        assert alignment.length == 84
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(2.67123, abs=5e-8)

        alignment = alignments[166]
        assert alignment.title[:50] == "gi|21674434|ref|NP_662499.1| Sec-independent prote"
        assert alignment.length == 67
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(2.67123, abs=5e-8)

        alignment = alignments[167]
        assert alignment.title[:50] == "gi|39968009|ref|XP_365395.1| hypothetical protein "
        assert alignment.length == 823
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(3.48874, abs=5e-8)

        alignment = alignments[168]
        assert alignment.title[:50] == "gi|4877986|gb|AAD31523.1| THA9 [Zea mays]"
        assert alignment.length == 169
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(3.48874, abs=5e-8)

        alignment = alignments[169]
        assert alignment.title[:50] == "gi|67934419|ref|ZP_00527476.1| Twin-arginine trans"
        assert alignment.length == 56
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(3.48874, abs=5e-8)

        alignment = alignments[170]
        assert alignment.title[:50] == "gi|42523658|ref|NP_969038.1| twin-argine protein t"
        assert alignment.length == 79
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(3.48874, abs=5e-8)

        alignment = alignments[171]
        assert alignment.title[:50] == "gi|71546080|ref|ZP_00666945.1| Twin-arginine trans"
        assert alignment.length == 73
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(3.48874, abs=5e-8)

        alignment = alignments[172]
        assert alignment.title[:50] == "gi|68002197|ref|ZP_00534828.1| Twin-arginine trans"
        assert alignment.length == 60
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(4.55643, abs=5e-8)

        alignment = alignments[173]
        assert alignment.title[:50] == "gi|67481641|ref|XP_656170.1| hypothetical protein "
        assert alignment.length == 434
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(4.55643, abs=5e-8)

        alignment = alignments[174]
        assert alignment.title[:50] == "gi|50935447|ref|XP_477251.1| putative Calreticulin"
        assert alignment.length == 424
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(4.55643, abs=5e-8)

        alignment = alignments[175]
        assert alignment.title[:50] == "gi|50978634|ref|NP_001003013.1| acidic (leucine-ri"
        assert alignment.length == 249
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(4.55643, abs=5e-8)

        alignment = alignments[176]
        assert alignment.title[:50] == "gi|70936814|ref|XP_739300.1| 40S ribosomal subunit"
        assert alignment.length == 184
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(4.55643, abs=5e-8)

        alignment = alignments[177]
        assert alignment.title[:50] == "gi|68075857|ref|XP_679848.1| hypothetical protein "
        assert alignment.length == 340
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(4.55643, abs=5e-8)

        alignment = alignments[178]
        assert alignment.title[:50] == "gi|39594005|emb|CAE70115.1| Hypothetical protein C"
        assert alignment.length == 192
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(4.55643, abs=5e-8)

        alignment = alignments[179]
        assert alignment.title[:50] == "gi|66809957|ref|XP_638702.1| hypothetical protein "
        assert alignment.length == 721
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(4.55643, abs=5e-8)

        alignment = alignments[180]
        assert alignment.title[:50] == "gi|68550463|ref|ZP_00589911.1| Twin-arginine trans"
        assert alignment.length == 69
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(4.55643, abs=5e-8)

        alignment = alignments[181]
        assert alignment.title[:50] == "gi|51473916|ref|YP_067673.1| TatA/E-like Sec-indep"
        assert alignment.length == 53
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(4.55643, abs=5e-8)

        alignment = alignments[182]
        assert alignment.title[:50] == "gi|61857708|ref|XP_612559.1| PREDICTED: similar to"
        assert alignment.length == 236
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(4.55643, abs=5e-8)

        alignment = alignments[183]
        assert alignment.title[:50] == "gi|39982651|gb|AAR34111.1| twin-arginine transloca"
        assert alignment.length == 59
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(4.55643, abs=5e-8)

        alignment = alignments[184]
        assert alignment.title[:50] == "gi|50877509|emb|CAG37349.1| related to Sec-indepen"
        assert alignment.length == 66
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(4.55643, abs=5e-8)

        alignment = alignments[185]
        assert alignment.title[:50] == "gi|52699323|ref|ZP_00340731.1| COG1826: Sec-indepe"
        assert alignment.length == 53
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(4.55643, abs=5e-8)

        alignment = alignments[186]
        assert alignment.title[:50] == "gi|62426215|ref|ZP_00381343.1| COG1826: Sec-indepe"
        assert alignment.length == 93
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(4.55643, abs=5e-8)

        alignment = alignments[187]
        assert alignment.title[:50] == "gi|11131838|sp|Q9SLY8|CRTC_ORYSA Calreticulin prec"
        assert alignment.length == 424
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(4.55643, abs=5e-8)

        alignment = alignments[188]
        assert alignment.title[:50] == "gi|56543690|gb|AAV89844.1| Sec-independent protein"
        assert alignment.length == 87
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(4.55643, abs=5e-8)

        alignment = alignments[189]
        assert alignment.title[:50] == "gi|67923730|ref|ZP_00517196.1| Twin-arginine trans"
        assert alignment.length == 95
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(5.95088, abs=5e-8)

        alignment = alignments[190]
        assert alignment.title[:50] == "gi|67462585|ref|XP_647954.1| hypothetical protein "
        assert alignment.length == 140
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(5.95088, abs=5e-8)

        alignment = alignments[191]
        assert alignment.title[:50] == "gi|51970620|dbj|BAD44002.1| unknown protein [Arabi"
        assert alignment.length == 784
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(5.95088, abs=5e-8)

        alignment = alignments[192]
        assert alignment.title[:50] == "gi|34581241|ref|ZP_00142721.1| hypothetical protei"
        assert alignment.length == 53
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(5.95088, abs=5e-8)

        alignment = alignments[193]
        assert alignment.title[:50] == "gi|4877984|gb|AAD31522.1| THA4 [Zea mays]"
        assert alignment.length == 170
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(5.95088, abs=5e-8)

        alignment = alignments[194]
        assert alignment.title[:50] == "gi|9757886|dbj|BAB08393.1| unnamed protein product"
        assert alignment.length == 707
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(5.95088, abs=5e-8)

        alignment = alignments[195]
        assert alignment.title[:50] == "gi|32422107|ref|XP_331497.1| predicted protein [Ne"
        assert alignment.length == 216
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(5.95088, abs=5e-8)

        alignment = alignments[196]
        assert alignment.title[:50] == "gi|68552035|ref|ZP_00591428.1| Twin-arginine trans"
        assert alignment.length == 70
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(5.95088, abs=5e-8)

        alignment = alignments[197]
        assert alignment.title[:50] == "gi|68177649|ref|ZP_00550794.1| Twin-arginine trans"
        assert alignment.length == 58
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(5.95088, abs=5e-8)

        alignment = alignments[198]
        assert alignment.title[:50] == "gi|67934756|ref|ZP_00527782.1| Twin-arginine trans"
        assert alignment.length == 65
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(5.95088, abs=5e-8)

        alignment = alignments[199]
        assert alignment.title[:50] == "gi|42550455|gb|EAA73298.1| hypothetical protein FG"
        assert alignment.length == 297
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(5.95088, abs=5e-8)

        alignment = alignments[200]
        assert alignment.title[:50] == "gi|15893083|ref|NP_360797.1| hypothetical protein "
        assert alignment.length == 53
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(5.95088, abs=5e-8)

        alignment = alignments[201]
        assert alignment.title[:50] == "gi|57233621|ref|YP_182297.1| twin-arginine translo"
        assert alignment.length == 65
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(5.95088, abs=5e-8)

        alignment = alignments[202]
        assert alignment.title[:50] == "gi|75908036|ref|YP_322332.1| Twin-arginine translo"
        assert alignment.length == 56
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(5.95088, abs=5e-8)

        alignment = alignments[203]
        assert alignment.title[:50] == "gi|72383453|ref|YP_292808.1| Twin-arginine translo"
        assert alignment.length == 71
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(5.95088, abs=5e-8)

        alignment = alignments[204]
        assert alignment.title[:50] == "gi|1666185|emb|CAB04766.1| ORF13(1) [Rhodococcus e"
        assert alignment.length == 98
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(5.95088, abs=5e-8)

        alignment = alignments[205]
        assert alignment.title[:50] == "gi|72138252|ref|XP_800288.1| PREDICTED: hypothetic"
        assert alignment.length == 946
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(5.95088, abs=5e-8)

        alignment = alignments[206]
        assert alignment.title[:50] == "gi|67923190|ref|ZP_00516678.1| Twin-arginine trans"
        assert alignment.length == 50
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(7.7721, abs=5e-8)

        alignment = alignments[207]
        assert alignment.title[:50] == "gi|3329623|gb|AAC26930.1| Hypothetical protein F36"
        assert alignment.length == 335
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(7.7721, abs=5e-8)

        alignment = alignments[208]
        assert alignment.title[:50] == "gi|39597929|emb|CAE68621.1| Hypothetical protein C"
        assert alignment.length == 2691
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(7.7721, abs=5e-8)

        alignment = alignments[209]
        assert alignment.title[:50] == "gi|68182025|ref|ZP_00555006.1| hypothetical protei"
        assert alignment.length == 438
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(7.7721, abs=5e-8)

        alignment = alignments[210]
        assert alignment.title[:50] == "gi|21204492|dbj|BAB95189.1| ebh [Staphylococcus au"
        assert alignment.length == 9904
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(7.7721, abs=5e-8)

        alignment = alignments[211]
        assert alignment.title[:50] == "gi|39593039|emb|CAE64508.1| Hypothetical protein C"
        assert alignment.length == 960
        assert len(alignment.hsps) == 1
        assert alignment.hsps[0].expect == pytest.approx(7.7721, abs=5e-8)

    def test_xml_2212L_blastn_001(self):
        """Parsing BLASTN 2.2.12, gi|1348916|gb|G26684.1|G26684 (xml_2212L_blastn_001)."""
        filename = "xml_2212L_blastn_001.xml"
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)
            record = next(records)
            with pytest.raises(StopIteration):
                next(records)
        self.check_xml_2212L_blastn_001(record)

        with open(datafile, "rb") as handle:
            record = NCBIXML.read(handle)
        self.check_xml_2212L_blastn_001(record)

    def check_xml_2212L_blastn_001(self, record):
        assert record.reference == 'Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schäffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.'
        assert record.date == "Aug-07-2005"
        assert record.version == "2.2.12"
        assert record.database == "nr"
        assert record.application == "BLASTN"
        alignments = record.alignments

        assert record.query_id == "gi|1348916|gb|G26684.1|G26684"
        assert len(alignments) == 2
        assert sum(len(a.hsps) for a in alignments) == 2
        assert alignments[0].title[:50] == "gi|9950606|gb|AE004854.1| Pseudomonas aeruginosa P"
        assert alignments[0].length == 11884
        assert len(alignments[0].hsps) == 1
        assert alignments[0].hsps[0].expect == pytest.approx(1.0598, abs=5e-8)
        assert alignments[1].title[:50] == "gi|15073988|emb|AL591786.1|SME591786 Sinorhizobium"
        assert alignments[1].length == 299350
        assert len(alignments[1].hsps) == 1
        assert alignments[1].hsps[0].expect == pytest.approx(4.18768, abs=5e-8)

    def test_xml_2212L_blastx_001(self):
        """Parsing BLASTX 2.2.12, gi|1347369|gb|G25137.1|G25137 (xml_2212L_blastx_001)."""
        filename = "xml_2212L_blastx_001.xml"
        datafile = os.path.join("Blast", filename)

        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)
            record = next(records)
            with pytest.raises(StopIteration):
                next(records)
        self.check_xml_2212L_blastx_001(record)

        with open(datafile, "rb") as handle:
            record = NCBIXML.read(handle)
        self.check_xml_2212L_blastx_001(record)

    def check_xml_2212L_blastx_001(self, record):
        assert record.reference == 'Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schäffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.'
        assert record.date == "Aug-07-2005"
        assert record.version == "2.2.12"
        assert record.database == "nr"
        assert record.application == "BLASTX"
        alignments = record.alignments
        assert record.query_id == "gi|1347369|gb|G25137.1|G25137"
        assert len(alignments) == 78
        assert sum(len(a.hsps) for a in alignments) == 84

        hsp = record.alignments[0].hsps[0]
        assert hsp.score == 630.0
        assert hsp.bits == 247.284

    def test_xml_2212L_tblastn_001(self):
        """Parsing TBLASTN 2.2.12, gi|729325|sp|P39483|DHG2_BACME (xml_2212L_tblastn_001)."""
        filename = "xml_2212L_tblastn_001.xml"
        datafile = os.path.join("Blast", filename)

        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)
            record = next(records)
            with pytest.raises(StopIteration):
                next(records)
        self.check_xml_2212L_tblastn_001(record)

        with open(datafile, "rb") as handle:
            record = NCBIXML.read(handle)
        self.check_xml_2212L_tblastn_001(record)

    def check_xml_2212L_tblastn_001(self, record):
        assert record.reference == 'Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schäffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.'
        assert record.date == "Aug-07-2005"
        assert record.version == "2.2.12"
        assert record.database == "nr"
        assert record.application == "TBLASTN"
        alignments = record.alignments
        assert record.query_id == "gi|729325|sp|P39483|DHG2_BACME"
        assert len(alignments) == 100
        assert sum(len(a.hsps) for a in alignments) == 127

    def test_xml_2212L_tblastx_001(self):
        """Parsing TBLASTX 2.2.12, gi|1348853|gb|G26621.1|G26621, BLOSUM80 (xml_2212L_tblastx_001)."""
        filename = "xml_2212L_tblastx_001.xml"
        datafile = os.path.join("Blast", filename)

        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)
            record = next(records)
            with pytest.raises(StopIteration):
                next(records)
        self.check_xml_2212L_tblastx_001(record)

        with open(datafile, "rb") as handle:
            record = NCBIXML.read(handle)
        self.check_xml_2212L_tblastx_001(record)

    def check_xml_2212L_tblastx_001(self, record):
        assert record.reference == 'Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schäffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.'
        assert record.date == "Aug-07-2005"
        assert record.version == "2.2.12"
        assert record.database == "nr"
        assert record.application == "TBLASTX"
        alignments = record.alignments
        assert record.query_id == "gi|1348853|gb|G26621.1|G26621"
        assert len(alignments) == 10
        assert sum(len(a.hsps) for a in alignments) == 102

    def test_xml_2218_blastp_001(self):
        """Parsing BLASTP 2.2.18+, gi|160837788|ref|NP_075631.2| (xml_2218_blastp_001)."""
        # NOTE - no date in version field, downloaded 2008/05/08

        filename = "xml_2218_blastp_001.xml"
        datafile = os.path.join("Blast", filename)

        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)
            record = next(records)
            with pytest.raises(StopIteration):
                next(records)
        self.check_xml_2218_blastp_001(record)

        with open(datafile, "rb") as handle:
            record = NCBIXML.read(handle)
            handle.close()
        self.check_xml_2218_blastp_001(record)

    def check_xml_2218_blastp_001(self, record):
        assert record.reference == 'Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schäffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.'
        assert record.date == ""
        assert record.version == "2.2.18+"
        assert record.database == "nr"
        assert record.application == "BLASTP"
        alignments = record.alignments
        assert record.query_id == "31493"
        assert len(alignments) == 10
        assert sum(len(a.hsps) for a in alignments) == 14
        assert alignments[0].title[:50] == "gi|151942244|gb|EDN60600.1| cytosolic iron-sulfur "
        assert alignments[0].length == 330
        assert len(alignments[0].hsps) == 1
        assert alignments[0].hsps[0].expect == pytest.approx(0.0185319, abs=5e-8)
        assert alignments[1].title[:50] == "gi|476059|emb|CAA55606.1| YBR0832 [Saccharomyces c"
        assert alignments[1].length == 535
        assert len(alignments[1].hsps) == 1
        assert alignments[1].hsps[0].expect == pytest.approx(0.0185319, abs=5e-8)
        assert alignments[2].title[:50] == "gi|6320473|ref|NP_010553.1| Essential protein invo"
        assert alignments[2].length == 330
        assert len(alignments[2].hsps) == 1
        assert alignments[2].hsps[0].expect == pytest.approx(0.0185319, abs=5e-8)
        assert alignments[3].title[:50] == "gi|61679798|pdb|1R5M|A Chain A, Crystal Structure "
        assert alignments[3].length == 425
        assert len(alignments[3].hsps) == 1
        assert alignments[3].hsps[0].expect == pytest.approx(0.0185319, abs=5e-8)
        assert alignments[4].title[:50] == "gi|6319579|ref|NP_009661.1| WD40 repeat-containing"
        assert alignments[4].length == 535
        assert len(alignments[4].hsps) == 1
        assert alignments[4].hsps[0].expect == pytest.approx(0.0185319, abs=5e-8)
        assert alignments[5].title[:50] == "gi|151946495|gb|EDN64717.1| Sir4p-interacting fact"
        assert alignments[5].length == 535
        assert len(alignments[5].hsps) == 1
        assert alignments[5].hsps[0].expect == pytest.approx(0.0412849, abs=5e-8)
        assert alignments[6].title[:50] == "gi|151943708|gb|EDN62018.1| nuclear pore complex s"
        assert alignments[6].length == 349
        assert len(alignments[6].hsps) == 2
        assert alignments[6].hsps[0].expect == pytest.approx(0.0704213, abs=5e-8)
        assert alignments[6].hsps[1].expect == pytest.approx(2.26538, abs=5e-8)
        assert alignments[7].title[:50] == "gi|151567866|pdb|2PM7|B Chain B, Crystal Structure"
        assert alignments[7].length == 297
        assert len(alignments[7].hsps) == 2
        assert alignments[7].hsps[0].expect == pytest.approx(0.0704213, abs=5e-8)
        assert alignments[7].hsps[1].expect == pytest.approx(0.456458, abs=5e-8)
        assert alignments[8].title[:50] == "gi|6321338|ref|NP_011415.1| Nuclear pore protein t"
        assert alignments[8].length == 349
        assert len(alignments[8].hsps) == 2
        assert alignments[8].hsps[0].expect == pytest.approx(0.0704213, abs=5e-8)
        assert alignments[8].hsps[1].expect == pytest.approx(2.26538, abs=5e-8)
        assert alignments[9].title[:50] == "gi|151567870|pdb|2PM9|B Chain B, Crystal Structure"
        assert alignments[9].length == 297
        assert len(alignments[9].hsps) == 2
        assert alignments[9].hsps[0].expect == pytest.approx(0.0919731, abs=5e-8)
        assert alignments[9].hsps[1].expect == pytest.approx(0.267601, abs=5e-8)

    def test_xml_2218_blastp_002(self):
        """Parsing BLASTP 2.2.18+, SwissProt Q08386 and P07175, no hits (xml_2218_blastp_002)."""
        filename = "xml_2218_blastp_002.xml"
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)
            record = next(records)
            assert record.reference == 'Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schäffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.'
            assert record.date == ""
            assert record.version == "2.2.18+"
            assert record.database == "gpipe/9606/Previous/protein"
            assert record.application == "BLASTP"
            assert record.query_id == "gi|585505|sp|Q08386|MOPB_RHOCA"
            assert len(record.alignments) == 0
            record = next(records)
            assert record.reference == 'Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schäffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.'
            assert record.date == ""
            assert record.version == "2.2.18+"
            assert record.database == "gpipe/9606/Previous/protein"
            assert record.application == "BLASTP"
            assert record.query_id == "gi|129628|sp|P07175.1|PARA_AGRTU"
            assert len(record.alignments) == 0
            with pytest.raises(StopIteration):
                next(records)

    def test_xml_2218L_blastp_001(self):
        """Parsing BLASTP 2.2.18, Fake query (xml_2218L_blastp_001)."""
        filename = "xml_2218L_blastp_001.xml"
        datafile = os.path.join("Blast", filename)

        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)
            record = next(records)
            with pytest.raises(StopIteration):
                next(records)
        self.check_xml_2218L_blastp_001(record)

        with open(datafile, "rb") as handle:
            record = NCBIXML.read(handle)
        self.check_xml_2218L_blastp_001(record)

    def check_xml_2218L_blastp_001(self, record):
        assert record.reference == '~Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaffer, ~Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), ~"Gapped BLAST and PSI-BLAST: a new generation of protein database search~programs",  Nucleic Acids Res. 25:3389-3402.'
        assert record.date == "Mar-02-2008"
        assert record.version == "2.2.18"
        assert record.database == "/Users/pjcock/Downloads/Software/blast-2.2.18/data/nr"
        assert record.application == "BLASTP"
        assert record.query_id == "lcl|1_0"
        alignments = record.alignments
        assert len(alignments) == 0

    def test_xml_2222_blastx_001(self):
        """Parsing BLASTX 2.2.22+, multiple queries against NR (xml_2222_blastx_001)."""
        # See also plain text file bt081.txt (matching output from blastx tool)

        filename = "xml_2222_blastx_001.xml"
        datafile = os.path.join("Blast", filename)

        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)

            record = next(records)
            assert record.reference == 'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.'
            assert record.database == "nr"
            assert record.application == "BLASTX"
            assert record.version == "2.2.22+"
            assert record.date == ""
            assert record.query == "gi|4104054|gb|AH007193.1|SEG_CVIGS Centaurea vallesiaca 18S ribosomal RNA gene, partial sequence"
            assert record.query_letters == 1002
            assert record.num_sequences_in_database == 8994603
            assert record.database_sequences == 8994603
            # self.assertEqual(record.database_length, 3078807967)
            assert record.database_length == -1216159329  # NCBI bug!
            assert len(record.descriptions) == 1
            assert len(record.alignments) == 1
            assert len(record.alignments[0].hsps) == 1

            record = next(records)
            assert record.reference == 'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.'
            assert record.database == "nr"
            assert record.application == "BLASTX"
            assert record.version == "2.2.22+"
            assert record.date == ""
            assert record.query == "gi|4218935|gb|AF074388.1|AF074388 Sambucus nigra hevein-like protein HLPf gene, partial cds"
            assert record.query_letters == 2050
            assert record.num_sequences_in_database == 8994603
            assert record.database_sequences == 8994603
            # self.assertEqual(record.database_length, 3078807967)
            assert record.database_length == -1216159329  # NCBI bug!
            # I used -num_descriptions 10 and -num_alignments 1
            assert len(record.descriptions) == 10
            assert len(record.alignments) == 10
            assert len(record.alignments[0].hsps) == 2
            assert len(record.alignments[1].hsps) == 2
            assert len(record.alignments[9].hsps) == 2

            record = next(records)
            assert record.reference == 'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.'
            assert record.database == "nr"
            assert record.application == "BLASTX"
            assert record.version == "2.2.22+"
            assert record.date == ""
            assert record.query == "gi|5690369|gb|AF158246.1|AF158246 Cricetulus griseus glucose phosphate isomerase (GPI) gene, partial intron sequence"
            assert record.query_letters == 550
            assert record.num_sequences_in_database == 8994603
            assert record.database_sequences == 8994603
            # self.assertEqual(record.database_length, 3078807967)
            assert record.database_length == -1216159329  # NCBI bug!
            assert len(record.descriptions) == 0
            assert len(record.alignments) == 0

            record = next(records)
            assert record.reference == 'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.'
            assert record.database == "nr"
            assert record.application == "BLASTX"
            assert record.version == "2.2.22+"
            assert record.date == ""
            assert record.query == "gi|5049839|gb|AI730987.1|AI730987 BNLGHi8354 Six-day Cotton fiber Gossypium hirsutum cDNA 5' similar to TUBULIN BETA-1 CHAIN gi|486734|pir|S35142 tubulin beta chain - white lupine gi|402636 (X70184) Beta tubulin 1 [Lupinus albus], mRNA sequence"
            assert record.query_letters == 655
            assert record.num_sequences_in_database == 8994603
            assert record.database_sequences == 8994603
            # self.assertEqual(record.database_length, 3078807967)
            assert record.database_length == -1216159329  # NCBI bug!
            assert len(record.descriptions) == 10
            assert len(record.alignments) == 10
            assert len(record.alignments[0].hsps) == 1
            assert len(record.alignments[9].hsps) == 1

            record = next(records)
            assert record.reference == 'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.'
            assert record.database == "nr"
            assert record.application == "BLASTX"
            assert record.version == "2.2.22+"
            assert record.date == ""
            assert record.query == "gi|5052071|gb|AF067555.1|AF067555 Phlox stansburyi internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence"
            assert record.query_letters == 623
            assert record.num_sequences_in_database == 8994603
            assert record.database_sequences == 8994603
            # self.assertEqual(record.database_length, 3078807967)
            assert record.database_length == -1216159329  # NCBI bug!
            assert len(record.descriptions) == 10
            assert len(record.alignments) == 10
            assert len(record.alignments[0].hsps) == 2
            assert len(record.alignments[9].hsps) == 1

            record = next(records)
            assert record.reference == 'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.'
            assert record.database == "nr"
            assert record.application == "BLASTX"
            assert record.version == "2.2.22+"
            assert record.date == ""
            assert record.query == "gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA) gene, partial cds"
            assert record.query_letters == 309
            assert record.num_sequences_in_database == 8994603
            assert record.database_sequences == 8994603
            # self.assertEqual(record.database_length, 3078807967)
            assert record.database_length == -1216159329  # NCBI bug!
            assert len(record.descriptions) == 10
            assert len(record.alignments) == 10
            assert len(record.alignments[0].hsps) == 1
            assert len(record.alignments[9].hsps) == 1

            record = next(records)
            assert record.reference == 'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.'
            assert record.database == "nr"
            assert record.application == "BLASTX"
            assert record.version == "2.2.22+"
            assert record.date == ""
            assert record.query == "gi|5817701|gb|AF142731.1|AF142731 Wisteria frutescens maturase-like protein (matK) gene, complete cds; chloroplast gene for chloroplast product"
            assert record.query_letters == 2551
            assert record.num_sequences_in_database == 8994603
            assert record.database_sequences == 8994603
            # self.assertEqual(record.database_length, 3078807967)
            assert record.database_length == -1216159329  # NCBI bug!
            assert len(record.descriptions) == 10
            assert len(record.alignments) == 10
            assert len(record.alignments[0].hsps) == 1
            assert len(record.alignments[9].hsps) == 1

            with pytest.raises(StopIteration):
                next(records)

    def test_xml_2222_blastp_001(self):
        """Parsing BLASTP 2.2.22+, multiple queries against NR (xml_2222_blastp_001)."""
        # This is from blastp NOT blastall

        filename = "xml_2222_blastp_001.xml"
        datafile = os.path.join("Blast", filename)

        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)

            record = next(records)
            assert record.application == "BLASTP"
            assert record.version == "2.2.22+"
            assert record.date == ""
            assert record.reference == 'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.'
            assert record.database == "nr"
            assert record.query == "gi|3298468|dbj|BAA31520.1| SAMIPF"
            assert record.query_letters == 107
            assert record.num_sequences_in_database == 8994603
            assert record.database_sequences == 8994603
            # self.assertEqual(record.database_length, 3078807967)
            assert record.database_length == -1216159329  # NCBI bug!
            assert len(record.descriptions) == 10
            assert len(record.alignments) == 10
            assert len(record.alignments[0].hsps) == 1

            record = next(records)
            assert record.application == "BLASTP"
            assert record.version == "2.2.22+"
            assert record.date == ""
            assert record.reference == 'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.'
            assert record.database == "nr"
            assert record.query == "gi|2781234|pdb|1JLY|B Chain B, Crystal Structure Of Amaranthus Caudatus Agglutinin"
            assert record.query_letters == 304

            record = next(records)
            assert record.application == "BLASTP"
            assert record.version == "2.2.22+"
            assert record.date == ""
            assert record.reference == 'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.'
            assert record.database == "nr"
            assert record.query == "gi|4959044|gb|AAD34209.1|AF069992_1 LIM domain interacting RING finger protein"
            assert record.query_letters == 600

            record = next(records)
            assert record.application == "BLASTP"
            assert record.version == "2.2.22+"
            assert record.date == ""
            assert record.reference == 'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.'
            assert record.database == "nr"
            assert record.query == "gi|671626|emb|CAA85685.1| rubisco large subunit"
            assert record.query_letters == 473

            with pytest.raises(StopIteration):
                next(records)

    def test_xml_2218L_rpsblast_001(self):
        """Parsing PSI-BLASTP 2.2.18, single query which converges in 3 iterations (xml_2218L_rpsblast_001)."""
        # This is from old pgpblast command line tool, NOT new psiblast
        # NOTE - The parser currently returns three BLAST record objects.
        # The old text parser would return a single PSI BLAST record object with three rounds.
        # This may change... although it may require a PSI BLAST specific XML parser.

        filename = "xml_2218L_rpsblast_001.xml"
        datafile = os.path.join("Blast", filename)

        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)

            record = next(records)
            assert record.reference == '~Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaffer, ~Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), ~"Gapped BLAST and PSI-BLAST: a new generation of protein database search~programs",  Nucleic Acids Res. 25:3389-3402.'
            assert record.database == "/opt/BlastDBs/nr"
            assert record.application == "BLASTP"
            assert record.version == "2.2.18"
            assert record.date == "Mar-02-2008"
            assert record.query == "tr|Q3V4Q3|Q3V4Q3_9VIRU"
            assert record.query_letters == 131
            assert record.num_sequences_in_database == 2563094
            assert record.database_sequences == 2563094
            assert record.database_length == 864488805
            assert len(record.descriptions) == 11
            assert len(record.alignments) == 11
            assert len(record.alignments[0].hsps) == 1
            hsp = record.alignments[0].hsps[0]
            assert hsp.align_length == 131
            assert hsp.identities == 131
            assert hsp.positives == 131
            assert hsp.query == "MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQAMAS"
            assert hsp.sbjct == "MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQAMAS"
            assert hsp.match == "MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQAMAS"
            assert hsp.score == 680
            assert hsp.expect == 4.72196e-70
            assert hsp.query_start == 1
            assert hsp.query_end == 131
            assert hsp.sbjct_start == 1
            assert hsp.sbjct_end == 131
            assert len(record.alignments[1].hsps) == 1
            hsp = record.alignments[1].hsps[0]
            assert hsp.align_length == 77
            assert hsp.identities == 36
            assert hsp.positives == 49
            assert hsp.query == "MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQ"
            assert hsp.sbjct == "MAREEPYKGDYVGGVAKILQGYFANYYGFPNVSLRLAGEEANLSKTGHANAKAIVHEMIKVIKEASKPLR-RGKGFK"
            assert hsp.match == "MA+ EP KGDY GG  KIL  +     G+P V+L+LAGEEAN  + G    K  +H ++K+I +A KP R +G GF+"
            assert hsp.score == 181
            assert hsp.expect == 3.03476e-12
            assert hsp.query_start == 1
            assert hsp.query_end == 77
            assert hsp.sbjct_start == 1
            assert hsp.sbjct_end == 76

            record = next(records)
            assert record.reference == '~Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaffer, ~Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), ~"Gapped BLAST and PSI-BLAST: a new generation of protein database search~programs",  Nucleic Acids Res. 25:3389-3402.'
            assert record.database == "/opt/BlastDBs/nr"
            assert record.application == "BLASTP"
            assert record.version == "2.2.18"
            assert record.date == "Mar-02-2008"
            assert record.query == "tr|Q3V4Q3|Q3V4Q3_9VIRU"
            assert record.query_letters == 131
            assert record.num_sequences_in_database == 2563094
            assert record.database_sequences == 2563094
            assert record.database_length == 864488805
            assert len(record.descriptions) == 19
            assert len(record.alignments) == 19
            assert len(record.alignments[0].hsps) == 1
            hsp = record.alignments[0].hsps[0]
            assert hsp.align_length == 131
            assert hsp.identities == 131
            assert hsp.positives == 131
            assert hsp.query == "MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQAMAS"
            assert hsp.sbjct == "MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQAMAS"
            assert hsp.match == "MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQAMAS"
            assert hsp.score == 590
            assert hsp.expect == 1.28615e-59
            assert hsp.query_start == 1
            assert hsp.query_end == 131
            assert hsp.sbjct_start == 1
            assert hsp.sbjct_end == 131

            record = next(records)
            assert record.reference == '~Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaffer, ~Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), ~"Gapped BLAST and PSI-BLAST: a new generation of protein database search~programs",  Nucleic Acids Res. 25:3389-3402.'
            assert record.database == "/opt/BlastDBs/nr"
            assert record.application == "BLASTP"
            assert record.version == "2.2.18"
            assert record.date == "Mar-02-2008"
            assert record.query == "tr|Q3V4Q3|Q3V4Q3_9VIRU"
            assert record.query_letters == 131
            assert record.num_sequences_in_database == 2563094
            assert record.database_sequences == 2563094
            assert record.database_length == 864488805
            assert len(record.descriptions) == 9
            assert len(record.alignments) == 9
            assert len(record.alignments[0].hsps) == 1
            hsp = record.alignments[0].hsps[0]
            assert hsp.align_length == 131
            assert hsp.identities == 131
            assert hsp.positives == 131
            assert hsp.query == "MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQAMAS"
            assert hsp.sbjct == "MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQAMAS"
            assert hsp.match == "MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQAMAS"
            assert hsp.score == 535
            assert hsp.expect == 3.43623e-53
            assert hsp.query_start == 1
            assert hsp.query_end == 131
            assert hsp.sbjct_start == 1
            assert hsp.sbjct_end == 131

            # TODO - Can we detect the convergence status:
            # <Iteration_message>CONVERGED</Iteration_message>
            with pytest.raises(StopIteration):
                next(records)

    def test_xml_2900_blastp_001_v1(self):
        record = self._test_xml_2900_blastp_001("xml_2900_blastp_001.xml")

        description = record.descriptions[0]
        assert len(description.title) == 4706
        assert description.title[:300] == "gi|447157535|ref|WP_001234791.1| MULTISPECIES: Sec-independent protein translocase subunit TatA [Shigella] >gi|24115132|ref|NP_709642.1| twin-arginine translocation protein TatA [Shigella flexneri 2a str. 301] >gi|82778983|ref|YP_405332.1| twin-arginine translocation protein TatA [Shigella dysenteri"
        description = record.descriptions[1]
        assert len(description.title) == 106
        assert description.title == "gi|91074959|gb|ABE09840.1| sec-independent twin-arginine translocase subunit tatA [Escherichia coli UTI89]"
        description = record.descriptions[2]
        assert len(description.title) == 979
        assert description.title[:300] == "gi|73857826|gb|AAZ90533.1| conserved hypothetical protein [Shigella sonnei Ss046] >gi|331067587|gb|EGI38991.1| Sec-independent protein translocase protein TatA [Escherichia coli TA280] >gi|412965214|emb|CCK49144.1| sec-independent protein translocase protein tata/e homolog 2 [Escherichia coli chi712"
        description = record.descriptions[3]
        assert len(description.title) == 260
        assert description.title == "gi|25302684|pir||D86071 hypothetical protein tatA [imported] - Escherichia coli (strain O157:H7, substrain EDL933) >gi|12518713|gb|AAG59032.1|AE005614_12 twin arginine translocation protein; sec-independent protein export [Escherichia coli O157:H7 str. EDL933]"
        description = record.descriptions[4]
        assert len(description.title) == 100
        assert description.title == "gi|331072495|gb|EGI43827.1| Sec-independent protein translocase protein TatA [Escherichia coli H591]"
        description = record.descriptions[5]
        assert len(description.title) == 120
        assert description.title == "gi|808042844|pdb|2MN7|A Chain A, Solution structure of monomeric TatA of twin-arginine translocation system from E. coli"
        description = record.descriptions[6]
        assert len(description.title) == 774516
        assert description.title[:300] == "gi|481023661|ref|WP_001295260.1| MULTISPECIES: Sec-independent protein translocase subunit TatA [Proteobacteria] >gi|15834020|ref|NP_312793.1| TatABCE protein translocation system subunit TatA [Escherichia coli O157:H7 str. Sakai] >gi|90111653|ref|NP_418280.4| twin arginine protein translocation sys"
        description = record.descriptions[7]
        assert len(description.title) == 238
        assert description.title == "gi|808042842|pdb|2MN6|B Chain B, Solution structure of dimeric TatA of twin-arginine translocation system from E. coli >gi|808042843|pdb|2MN6|A Chain A, Solution structure of dimeric TatA of twin-arginine translocation system from E. coli"
        description = record.descriptions[8]
        assert len(description.title) == 2111
        assert description.title[:300] == "gi|491167042|ref|WP_005025412.1| MULTISPECIES: twin-arginine translocase subunit TatA [Enterobacteriaceae] >gi|320176770|gb|EFW51804.1| Twin-arginine translocation protein TatA [Shigella dysenteriae CDC 74-1112] >gi|391297993|gb|EIQ56018.1| twin arginine-targeting translocase, TatA/E family protein "
        description = record.descriptions[9]
        assert len(description.title) == 384
        assert description.title[:300] == "gi|332996954|gb|EGK16572.1| sec-independent translocase protein tatA [Shigella flexneri VA-6] >gi|391245379|gb|EIQ04650.1| twin arginine-targeting translocase, TatA/E family protein [Shigella flexneri K-1770] >gi|1411457050|emb|SRN34259.1| twin arginine translocase protein A [Shigella flexneri] >gi|"

    def test_xml_2900_blastp_001_v2(self):
        record = self._test_xml_2900_blastp_001("xml_2900_blastp_001_v2.xml")

        alignment = record.alignments[0]
        assert alignment.title == "gi|447157535|ref|WP_001234791.1| MULTISPECIES: Sec-independent protein translocase subunit TatA [Shigella]"

        description = record.descriptions[0]
        assert description.title == "gi|447157535|ref|WP_001234791.1| MULTISPECIES: Sec-independent protein translocase subunit TatA [Shigella]"
        assert len(description.items) == 48
        description_item = description.items[0]
        assert description_item.id == "gi|447157535|ref|WP_001234791.1|"
        assert description_item.accession == "WP_001234791"
        assert description_item.title == "MULTISPECIES: Sec-independent protein translocase subunit TatA [Shigella]"
        assert description_item.taxid == 620
        assert description_item.sciname == "Shigella"

    def _test_xml_2900_blastp_001(self, filename):
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)
            record = next(records)
            with pytest.raises(StopIteration):
                next(records)

        assert record.application == "BLASTP"
        assert record.version == "2.9.0+"
        assert record.reference == 'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.'
        assert record.database == "nr"
        assert record.date == ""
        assert record.query == "twin argininte translocase protein A [Escherichia coli K12]"
        assert record.query_letters == 103
        assert record.num_sequences_in_database == 194611632
        assert record.database_sequences == 194611632
        assert record.database_length == 2104817704

        alignment = record.alignments[0]
        assert alignment.hit_id == "gi|447157535|ref|WP_001234791.1|"
        assert alignment.accession == "WP_001234791"
        assert alignment.length == 103
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 103
        assert hsp.identities == 103
        assert hsp.positives == 103
        assert hsp.query == "MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV"
        assert hsp.sbjct == "MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV"
        assert hsp.match == "MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV"
        assert hsp.score == 532.0
        assert hsp.expect == 1.71849e-68
        assert hsp.query_start == 1
        assert hsp.query_end == 103
        assert hsp.sbjct_start == 1
        assert hsp.sbjct_end == 103

        alignment = record.alignments[1]
        assert alignment.hit_id == "gi|91074959|gb|ABE09840.1|"
        assert alignment.accession == "ABE09840"
        assert alignment.title == "gi|91074959|gb|ABE09840.1| sec-independent twin-arginine translocase subunit tatA [Escherichia coli UTI89]"
        assert alignment.length == 103
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 103
        assert hsp.identities == 102
        assert hsp.positives == 102
        assert hsp.query == "MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV"
        assert hsp.sbjct == "MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKIEDAKRHDKEQV"
        assert hsp.match == "MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAK EDAKRHDKEQV"
        assert hsp.score == 526.0
        assert hsp.expect == 1.44614e-67
        assert hsp.query_start == 1
        assert hsp.query_end == 103
        assert hsp.sbjct_start == 1
        assert hsp.sbjct_end == 103

        alignment = record.alignments[2]
        assert alignment.hit_id == "gi|73857826|gb|AAZ90533.1|"
        assert alignment.accession == "AAZ90533"
        assert alignment.length == 103
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 103
        assert hsp.identities == 102
        assert hsp.positives == 102
        assert hsp.query == "MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV"
        assert hsp.sbjct == "MRPCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV"
        assert hsp.match == "MR CLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV"
        assert hsp.score == 525.0
        assert hsp.expect == 1.80126e-67
        assert hsp.query_start == 1
        assert hsp.query_end == 103
        assert hsp.sbjct_start == 1
        assert hsp.sbjct_end == 103

        alignment = record.alignments[3]
        assert alignment.hit_id == "gi|25302684|pir||D86071"
        assert alignment.accession == "D86071"
        assert alignment.length == 103
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 103
        assert hsp.identities == 102
        assert hsp.positives == 102
        assert hsp.query == "MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV"
        assert hsp.sbjct == "MRLCLIIIYHRXTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV"
        assert hsp.match == "MRLCLIIIYHR TCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV"
        assert hsp.score == 525.0
        assert hsp.expect == 2.10054e-67
        assert hsp.query_start == 1
        assert hsp.query_end == 103
        assert hsp.sbjct_start == 1
        assert hsp.sbjct_end == 103

        alignment = record.alignments[4]
        assert alignment.hit_id == "gi|331072495|gb|EGI43827.1|"
        assert alignment.accession == "EGI43827"
        assert alignment.length == 103
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 103
        assert hsp.identities == 101
        assert hsp.positives == 102
        assert hsp.query == "MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV"
        assert hsp.sbjct == "MRPCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQANTNQEQAKTEDAKRHDKEQV"
        assert hsp.match == "MR CLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQA+TNQEQAKTEDAKRHDKEQV"
        assert hsp.score == 522.0
        assert hsp.expect == 5.11164e-67
        assert hsp.query_start == 1
        assert hsp.query_end == 103
        assert hsp.sbjct_start == 1
        assert hsp.sbjct_end == 103

        alignment = record.alignments[5]
        assert alignment.hit_id == "gi|808042844|pdb|2MN7|A"
        assert alignment.accession == "2MN7_A"
        assert alignment.length == 97
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 89
        assert hsp.identities == 89
        assert hsp.positives == 89
        assert hsp.query == "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV"
        assert hsp.sbjct == "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV"
        assert hsp.match == "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV"
        assert hsp.score == 449.0
        assert hsp.expect == 5.99855e-56
        assert hsp.query_start == 15
        assert hsp.query_end == 103
        assert hsp.sbjct_start == 1
        assert hsp.sbjct_end == 89

        alignment = record.alignments[6]
        assert alignment.hit_id == "gi|481023661|ref|WP_001295260.1|"
        assert alignment.accession == "WP_001295260"
        assert alignment.length == 89
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 89
        assert hsp.identities == 89
        assert hsp.positives == 89
        assert hsp.query == "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV"
        assert hsp.sbjct == "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV"
        assert hsp.match == "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV"
        assert hsp.score == 448.0
        assert hsp.expect == 8.86198e-56
        assert hsp.query_start == 15
        assert hsp.query_end == 103
        assert hsp.sbjct_start == 1
        assert hsp.sbjct_end == 89

        alignment = record.alignments[7]
        assert alignment.hit_id == "gi|808042842|pdb|2MN6|B"
        assert alignment.accession == "2MN6_B"
        assert alignment.length == 100
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 89
        assert hsp.identities == 89
        assert hsp.positives == 89
        assert hsp.query == "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV"
        assert hsp.sbjct == "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV"
        assert hsp.match == "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV"
        assert hsp.score == 447.0
        assert hsp.expect == 1.55019e-55
        assert hsp.query_start == 15
        assert hsp.query_end == 103
        assert hsp.sbjct_start == 4
        assert hsp.sbjct_end == 92

        alignment = record.alignments[8]
        assert alignment.hit_id == "gi|491167042|ref|WP_005025412.1|"
        assert alignment.accession == "WP_005025412"
        assert alignment.length == 89
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 89
        assert hsp.identities == 88
        assert hsp.positives == 89
        assert hsp.query == "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV"
        assert hsp.sbjct == "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDSDFTAKTIADKQADTNQEQAKTEDAKRHDKEQV"
        assert hsp.match == "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQD+DFTAKTIADKQADTNQEQAKTEDAKRHDKEQV"
        assert hsp.score == 445.0
        assert hsp.expect == 2.25306e-55
        assert hsp.query_start == 15
        assert hsp.query_end == 103
        assert hsp.sbjct_start == 1
        assert hsp.sbjct_end == 89

        alignment = record.alignments[9]
        assert alignment.hit_id == "gi|332996954|gb|EGK16572.1|"
        assert alignment.accession == "EGK16572"
        assert alignment.length == 89
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 89
        assert hsp.identities == 88
        assert hsp.positives == 89
        assert hsp.query == "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV"
        assert hsp.sbjct == "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTISDKQADTNQEQAKTEDAKRHDKEQV"
        assert hsp.match == "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTI+DKQADTNQEQAKTEDAKRHDKEQV"
        assert hsp.score == 445.0
        assert hsp.expect == 2.25306e-55
        assert hsp.query_start == 15
        assert hsp.query_end == 103
        assert hsp.sbjct_start == 1
        assert hsp.sbjct_end == 89

        assert len(record.descriptions) == 10

        description = record.descriptions[0]
        assert description.score == 532.0
        assert description.e == 1.71849e-68
        assert description.num_alignments == 1
        description = record.descriptions[1]
        assert description.score == 526.0
        assert description.e == 1.44614e-67
        assert description.num_alignments == 1
        description = record.descriptions[2]
        assert description.score == 525.0
        assert description.e == 1.80126e-67
        assert description.num_alignments == 1
        description = record.descriptions[3]
        assert description.score == 525.0
        assert description.e == 2.10054e-67
        assert description.num_alignments == 1
        description = record.descriptions[4]
        assert description.score == 522.0
        assert description.e == 5.11164e-67
        assert description.num_alignments == 1
        description = record.descriptions[5]
        assert description.score == 449.0
        assert description.e == 5.99855e-56
        assert description.num_alignments == 1
        description = record.descriptions[6]
        assert description.score == 448.0
        assert description.e == 8.86198e-56
        assert description.num_alignments == 1
        description = record.descriptions[7]
        assert description.score == 447.0
        assert description.e == 1.55019e-55
        assert description.num_alignments == 1
        description = record.descriptions[8]
        assert description.score == 445.0
        assert description.e == 2.25306e-55
        assert description.num_alignments == 1
        description = record.descriptions[9]
        assert description.score == 445.0
        assert description.e == 2.25306e-55
        assert description.num_alignments == 1

        return record

    def test_xml_2900_blastn_001_v1(self):
        self._test_xml_2900_blastn_001("xml_2900_blastn_001.xml")

    def test_xml_2900_blastn_001_v2(self):
        self._test_xml_2900_blastn_001("xml_2900_blastn_001_v2.xml")

    def _test_xml_2900_blastn_001(self, filename):
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)
            record = next(records)
            with pytest.raises(StopIteration):
                next(records)

        assert record.application == "BLASTN"
        assert record.version == "2.9.0+"
        assert record.reference == 'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.'
        assert record.database == "GPIPE/10090/current/all_top_level GPIPE/10090/current/rna"
        assert record.date == ""
        assert record.query == "human STS STS_D11570, sequence tagged site"
        assert record.query_letters == 285
        assert record.num_sequences_in_database == 107382
        assert record.database_sequences == 107382
        assert record.database_length == 3164670549

        alignment = record.alignments[0]
        assert alignment.hit_id == "gi|372099107|ref|NC_000069.6|"
        assert alignment.accession == "NC_000069"
        assert alignment.length == 160039680
        assert alignment.title == "gi|372099107|ref|NC_000069.6| Mus musculus strain C57BL/6J chromosome 3, GRCm38.p4 C57BL/6J"
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 34
        assert hsp.identities == 30
        assert hsp.positives == 30
        assert hsp.query == "GAATCCTAGAGGCTTGATTGGCCCAGG-CTGCTG"
        assert hsp.sbjct == "GAATCCTAGAGGCTGGACTGGCCCTGGCCTGCTG"
        assert hsp.match == "|||||||||||||| || |||||| || ||||||"
        assert hsp.score == 44.0
        assert hsp.expect == 0.375311
        assert hsp.query_start == 134
        assert hsp.query_end == 166
        assert hsp.sbjct_start == 101449177
        assert hsp.sbjct_end == 101449144
        assert hsp.frame == (1, -1)
        assert hsp.strand == ("Plus", "Minus")

        alignment = record.alignments[1]
        assert alignment.hit_id == "gi|372099103|ref|NC_000073.6|"
        assert alignment.accession == "NC_000073"
        assert alignment.length == 145441459
        assert alignment.title == "gi|372099103|ref|NC_000073.6| Mus musculus strain C57BL/6J chromosome 7, GRCm38.p4 C57BL/6J"
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 29
        assert hsp.identities == 26
        assert hsp.positives == 26
        assert hsp.query == "GAAAGGAAATNAAAATGGAAAGTTCTTGT"
        assert hsp.sbjct == "GAAAGGAAAAAAAAATGGAAAGTTCTGGT"
        assert hsp.match == "|||||||||  ||||||||||||||| ||"
        assert hsp.score == 44.0
        assert hsp.expect == 0.375311
        assert hsp.query_start == 205
        assert hsp.query_end == 233
        assert hsp.sbjct_start == 131772185
        assert hsp.sbjct_end == 131772157
        assert hsp.frame == (1, -1)
        assert hsp.strand == ("Plus", "Minus")

        alignment = record.alignments[2]
        assert alignment.hit_id == "gi|372099106|ref|NC_000070.6|"
        assert alignment.accession == "NC_000070"
        assert alignment.length == 156508116
        assert alignment.title == "gi|372099106|ref|NC_000070.6| Mus musculus strain C57BL/6J chromosome 4, GRCm38.p4 C57BL/6J"
        assert len(alignment.hsps) == 2
        hsp = alignment.hsps[0]
        assert hsp.align_length == 24
        assert hsp.identities == 23
        assert hsp.positives == 23
        assert hsp.query == "CCAACACAGGCCAGCGACTTCTGG"
        assert hsp.sbjct == "CCAACACAGGCCAGCGGCTTCTGG"
        assert hsp.match == "|||||||||||||||| |||||||"
        assert hsp.score == 43.0
        assert hsp.expect == 1.30996
        assert hsp.query_start == 62
        assert hsp.query_end == 85
        assert hsp.sbjct_start == 9607562
        assert hsp.sbjct_end == 9607539
        assert hsp.frame == (1, -1)
        assert hsp.strand == ("Plus", "Minus")
        hsp = alignment.hsps[1]
        assert hsp.align_length == 32
        assert hsp.identities == 28
        assert hsp.positives == 28
        assert hsp.query == "GCCTGACATGG-GTAGCTGCTCAATAAATGCT"
        assert hsp.sbjct == "GCCTGGCATGAAGTAACTGCTCAATAAATGCT"
        assert hsp.match == "||||| ||||  ||| ||||||||||||||||"
        assert hsp.score == 40.0
        assert hsp.expect == 4.57222
        assert hsp.query_start == 242
        assert hsp.query_end == 272
        assert hsp.sbjct_start == 142902532
        assert hsp.sbjct_end == 142902563
        assert hsp.frame == (1, 1)
        assert hsp.strand == ("Plus", "Plus")

        alignment = record.alignments[3]
        assert alignment.hit_id == "gi|372099108|ref|NC_000068.7|"
        assert alignment.accession == "NC_000068"
        assert alignment.length == 182113224
        assert alignment.title == "gi|372099108|ref|NC_000068.7| Mus musculus strain C57BL/6J chromosome 2, GRCm38.p4 C57BL/6J"
        assert len(alignment.hsps) == 2
        hsp = alignment.hsps[0]
        assert hsp.align_length == 31
        assert hsp.identities == 27
        assert hsp.positives == 27
        assert hsp.query == "AAGGCCTGACATGGGTAGCTGCTCAATAAAT"
        assert hsp.sbjct == "AAGTCCTGGCATGAGTAGTTGCTCAATAAAT"
        assert hsp.match == "||| |||| |||| |||| ||||||||||||"
        assert hsp.score == 42.0
        assert hsp.expect == 1.30996
        assert hsp.query_start == 239
        assert hsp.query_end == 269
        assert hsp.sbjct_start == 3799647
        assert hsp.sbjct_end == 3799677
        assert hsp.frame == (1, 1)
        assert hsp.strand == ("Plus", "Plus")
        hsp = alignment.hsps[1]
        assert hsp.align_length == 25
        assert hsp.identities == 23
        assert hsp.positives == 23
        assert hsp.query == "AAATNAAAATGGAAAGTTCTTGTAG"
        assert hsp.sbjct == "AAATGAAAATGGAAAGTTCTTATAG"
        assert hsp.match == "|||| |||||||||||||||| |||"
        assert hsp.score == 41.0
        assert hsp.expect == 4.57222
        assert hsp.query_start == 211
        assert hsp.query_end == 235
        assert hsp.sbjct_start == 70278960
        assert hsp.sbjct_end == 70278984
        assert hsp.frame == (1, 1)
        assert hsp.strand == ("Plus", "Plus")

        alignment = record.alignments[4]
        assert alignment.hit_id == "gi|372099097|ref|NC_000079.6|"
        assert alignment.accession == "NC_000079"
        assert alignment.length == 120421639
        assert alignment.title == "gi|372099097|ref|NC_000079.6| Mus musculus strain C57BL/6J chromosome 13, GRCm38.p4 C57BL/6J"
        assert len(alignment.hsps) == 2
        hsp = alignment.hsps[0]
        assert hsp.align_length == 28
        assert hsp.identities == 25
        assert hsp.positives == 25
        assert hsp.query == "AAGGAAATNAAAATGGAAAGTTCTTGTA"
        assert hsp.sbjct == "AAGGACATCAAAATGGAAAGTTCTTCTA"
        assert hsp.match == "||||| || |||||||||||||||| ||"
        assert hsp.score == 42.0
        assert hsp.expect == 1.30996
        assert hsp.query_start == 207
        assert hsp.query_end == 234
        assert hsp.sbjct_start == 26806584
        assert hsp.sbjct_end == 26806557
        assert hsp.frame == (1, -1)
        assert hsp.strand == ("Plus", "Minus")
        hsp = alignment.hsps[1]
        assert hsp.align_length == 40
        assert hsp.identities == 32
        assert hsp.positives == 32
        assert hsp.query == "AGCGCAAGGCCTGACATGGGTAGCTGCTCAATAAATGCTA"
        assert hsp.sbjct == "AGCGCAAGGCCTGACATAGGAAAATGTTCAGTGAATACTA"
        assert hsp.match == "||||||||||||||||| || |  || ||| | ||| |||"
        assert hsp.score == 40.0
        assert hsp.expect == 4.57222
        assert hsp.query_start == 234
        assert hsp.query_end == 273
        assert hsp.sbjct_start == 56840340
        assert hsp.sbjct_end == 56840301
        assert hsp.frame == (1, -1)
        assert hsp.strand == ("Plus", "Minus")

        alignment = record.alignments[5]
        assert alignment.hit_id == "gi|372099098|ref|NC_000078.6|"
        assert alignment.accession == "NC_000078"
        assert alignment.length == 120129022
        assert alignment.title == "gi|372099098|ref|NC_000078.6| Mus musculus strain C57BL/6J chromosome 12, GRCm38.p4 C57BL/6J"
        assert len(alignment.hsps) == 2
        hsp = alignment.hsps[0]
        assert hsp.align_length == 23
        assert hsp.identities == 22
        assert hsp.positives == 22
        assert hsp.query == "CATCCATTCACACCCAACACAGG"
        assert hsp.sbjct == "CATCCATTCACACCCAGCACAGG"
        assert hsp.match == "|||||||||||||||| ||||||"
        assert hsp.score == 41.0
        assert hsp.expect == 4.57222
        assert hsp.query_start == 49
        assert hsp.query_end == 71
        assert hsp.sbjct_start == 113030663
        assert hsp.sbjct_end == 113030685
        assert hsp.frame == (1, 1)
        assert hsp.strand == ("Plus", "Plus")
        hsp = alignment.hsps[1]
        assert hsp.align_length == 32
        assert hsp.identities == 28
        assert hsp.positives == 28
        assert hsp.query == "TGTAGCGCAAGGCCTGACATGGGTAGCTGCTC"
        assert hsp.sbjct == "TGTAGCTCTAGGCCTGACATGGGT-GCTGGTC"
        assert hsp.match == "|||||| | ||||||||||||||| |||| ||"
        assert hsp.score == 40.0
        assert hsp.expect == 4.57222
        assert hsp.query_start == 231
        assert hsp.query_end == 262
        assert hsp.sbjct_start == 108990272
        assert hsp.sbjct_end == 108990242
        assert hsp.frame == (1, -1)
        assert hsp.strand == ("Plus", "Minus")

        alignment = record.alignments[6]
        assert alignment.hit_id == "gi|372099109|ref|NC_000067.6|"
        assert alignment.accession == "NC_000067"
        assert alignment.length == 195471971
        assert alignment.title == "gi|372099109|ref|NC_000067.6| Mus musculus strain C57BL/6J chromosome 1, GRCm38.p4 C57BL/6J"
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 43
        assert hsp.identities == 35
        assert hsp.positives == 35
        assert hsp.query == "GCTCAGCCACAGACATGGTTTGTNACTNTTGAGCTTCTGTTCC"
        assert hsp.sbjct == "GCTCAGCCACATACATGGTTT-TAAGTGTTGAGGCTCT-TTCC"
        assert hsp.match == "||||||||||| ||||||||| | | | |||||  ||| ||||"
        assert hsp.score == 40.0
        assert hsp.expect == 4.57222
        assert hsp.query_start == 87
        assert hsp.query_end == 129
        assert hsp.sbjct_start == 65190108
        assert hsp.sbjct_end == 65190148
        assert hsp.frame == (1, 1)
        assert hsp.strand == ("Plus", "Plus")

        alignment = record.alignments[7]
        assert alignment.hit_id == "gi|372099101|ref|NC_000075.6|"
        assert alignment.accession == "NC_000075"
        assert alignment.length == 124595110
        assert alignment.title == "gi|372099101|ref|NC_000075.6| Mus musculus strain C57BL/6J chromosome 9, GRCm38.p4 C57BL/6J"
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 47
        assert hsp.identities == 36
        assert hsp.positives == 36
        assert hsp.query == "CAAGGCCTGACATGGGTAGCTGCTCAATAAATGCTAGTNTGTTATTT"
        assert hsp.sbjct == "CAAAGCCTGACAGGTATGACTGCTCAATAAATACTATTTTTTTTTTT"
        assert hsp.match == "||| |||||||| |  |  ||||||||||||| ||| | | || |||"
        assert hsp.score == 40.0
        assert hsp.expect == 4.57222
        assert hsp.query_start == 238
        assert hsp.query_end == 284
        assert hsp.sbjct_start == 58227241
        assert hsp.sbjct_end == 58227195
        assert hsp.frame == (1, -1)
        assert hsp.strand == ("Plus", "Minus")

        alignment = record.alignments[8]
        assert alignment.hit_id == "gi|372099100|ref|NC_000076.6|"
        assert alignment.accession == "NC_000076"
        assert alignment.length == 130694993
        assert alignment.title == "gi|372099100|ref|NC_000076.6| Mus musculus strain C57BL/6J chromosome 10, GRCm38.p4 C57BL/6J"
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 20
        assert hsp.identities == 20
        assert hsp.positives == 20
        assert hsp.query == "AGCTGCTCAATAAATGCTAG"
        assert hsp.sbjct == "AGCTGCTCAATAAATGCTAG"
        assert hsp.match == "||||||||||||||||||||"
        assert hsp.score == 40.0
        assert hsp.expect == 4.57222
        assert hsp.query_start == 255
        assert hsp.query_end == 274
        assert hsp.sbjct_start == 119337186
        assert hsp.sbjct_end == 119337205
        assert hsp.frame == (1, 1)
        assert hsp.strand == ("Plus", "Plus")

        alignment = record.alignments[9]
        assert alignment.hit_id == "gi|372099094|ref|NC_000082.6|"
        assert alignment.accession == "NC_000082"
        assert alignment.length == 98207768
        assert alignment.title == "gi|372099094|ref|NC_000082.6| Mus musculus strain C57BL/6J chromosome 16, GRCm38.p4 C57BL/6J"
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 56
        assert hsp.identities == 43
        assert hsp.positives == 43
        assert hsp.query == "GGAGGCAAAGAATCCCTACCTCCT-AGGGGTGA-AAGGAAATNAAAATGGAAAGTT"
        assert hsp.sbjct == "GGAGGCAAAGAATCCCTACATTGTGACAGCTGATAAAGAAGGTAAAATGGAAAATT"
        assert hsp.match == "||||||||||||||||||| |  | |  | ||| || |||   |||||||||| ||"
        assert hsp.score == 40.0
        assert hsp.expect == 4.57222
        assert hsp.query_start == 175
        assert hsp.query_end == 228
        assert hsp.sbjct_start == 18854780
        assert hsp.sbjct_end == 18854835
        assert hsp.frame == (1, 1)
        assert hsp.strand == ("Plus", "Plus")
        assert len(record.descriptions) == 10
        description = record.descriptions[0]
        assert description.score == 44.0
        assert description.e == 0.375311
        assert description.num_alignments == 1
        description = record.descriptions[1]
        assert description.score == 44.0
        assert description.e == 0.375311
        assert description.num_alignments == 1
        description = record.descriptions[2]
        assert description.score == 43.0
        assert description.e == 1.30996
        assert description.num_alignments == 2
        description = record.descriptions[3]
        assert description.score == 42.0
        assert description.e == 1.30996
        assert description.num_alignments == 2
        description = record.descriptions[4]
        assert description.score == 42.0
        assert description.e == 1.30996
        assert description.num_alignments == 2
        description = record.descriptions[5]
        assert description.score == 41.0
        assert description.e == 4.57222
        assert description.num_alignments == 2
        description = record.descriptions[6]
        assert description.score == 40.0
        assert description.e == 4.57222
        assert description.num_alignments == 1
        description = record.descriptions[7]
        assert description.score == 40.0
        assert description.e == 4.57222
        assert description.num_alignments == 1
        description = record.descriptions[8]
        assert description.score == 40.0
        assert description.e == 4.57222
        assert description.num_alignments == 1
        description = record.descriptions[9]
        assert description.score == 40.0
        assert description.e == 4.57222
        assert description.num_alignments == 1

    def test_xml_2900_blastx_001_v1(self):
        self._test_xml_2900_blastx_001("xml_2900_blastx_001.xml")

    def test_xml_2900_blastx_001_v2(self):
        self._test_xml_2900_blastx_001("xml_2900_blastx_001_v2.xml")

    def _test_xml_2900_blastx_001(self, filename):
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)
            record = next(records)

        assert record.application == "BLASTX"
        assert record.version == "2.9.0+"
        assert record.reference == 'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.'
        assert record.database == "nr"
        assert record.date == ""
        assert record.query == "MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA, complete cds, mRNA sequence"
        assert record.query_letters == 365

        alignment = record.alignments[0]
        assert alignment.hit_id == "gi|1530504495|emb|VDM03167.1|"
        assert alignment.accession == "VDM03167"
        assert alignment.length == 132
        assert alignment.title == "gi|1530504495|emb|VDM03167.1| unnamed protein product, partial [Schistocephalus solidus]"
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 108
        assert hsp.identities == 81
        assert hsp.positives == 83
        assert hsp.query == "MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.sbjct == "MADEEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.match == "MADEEVQALVVDNGSGMCKAG       ++  P               G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.score == 408.0
        assert hsp.bits == 161.77
        assert hsp.expect == 8.11609e-49
        assert hsp.query_start == 20
        assert hsp.query_end == 343
        assert hsp.sbjct_start == 1
        assert hsp.sbjct_end == 108
        assert hsp.frame == (2, 0)

        alignment = record.alignments[1]
        assert alignment.hit_id == "gi|510859078|gb|EPB74633.1|"
        assert alignment.accession == "EPB74633"
        assert alignment.length == 119
        assert alignment.title == "gi|510859078|gb|EPB74633.1| hypothetical protein ANCCEY_06263 [Ancylostoma ceylanicum]"
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 115
        assert hsp.identities == 81
        assert hsp.positives == 85
        assert hsp.query == "MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTELHCIRKP"
        assert hsp.sbjct == "MCDDDVAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTEAHSILKP"
        assert hsp.match == "M D++V ALVVDNGSGMCKAG       ++  P               G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE H I KP"
        assert hsp.score == 405.0
        assert hsp.expect == 1.40046e-48
        assert hsp.query_start == 20
        assert hsp.query_end == 364
        assert hsp.sbjct_start == 1
        assert hsp.sbjct_end == 115
        assert hsp.frame == (2, 0)

        alignment = record.alignments[2]
        assert alignment.hit_id == "gi|684409690|ref|XP_009175831.1|"
        assert alignment.accession == "XP_009175831"
        assert alignment.length == 246
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 108
        assert hsp.identities == 81
        assert hsp.positives == 83
        assert hsp.query == "MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.sbjct == "MADEEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.match == "MADEEVQALVVDNGSGMCKAG       ++  P               G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.score == 413.0
        assert hsp.expect == 4.40404e-48
        assert hsp.query_start == 20
        assert hsp.query_end == 343
        assert hsp.sbjct_start == 1
        assert hsp.sbjct_end == 108
        assert hsp.frame == (2, 0)

        alignment = record.alignments[3]
        assert alignment.hit_id == "gi|449710331|gb|EMD49430.1|"
        assert alignment.accession == "EMD49430"
        assert alignment.length == 124
        assert alignment.title == "gi|449710331|gb|EMD49430.1| actin, putative, partial [Entamoeba histolytica KU27]"
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 108
        assert hsp.identities == 78
        assert hsp.positives == 81
        assert hsp.query == "MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.sbjct == "MGDEEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHVSVMAGMGQKDAYVGDEAQSKRGILTLKYPIEHGIVNNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.match == "M DEEVQALVVDNGSGMCKAG       ++  P               G KD+YVGDEAQSKRGILTLKYPIEHGIV NWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.score == 401.0
        assert hsp.expect == 9.0486e-48
        assert hsp.query_start == 20
        assert hsp.query_end == 343
        assert hsp.sbjct_start == 1
        assert hsp.sbjct_end == 108
        assert hsp.frame == (2, 0)

        alignment = record.alignments[4]
        assert alignment.hit_id == "gi|257215766|emb|CAX83035.1|"
        assert alignment.accession == "CAX83035"
        assert alignment.length == 252
        assert alignment.title == "gi|257215766|emb|CAX83035.1| Actin-2, partial [Schistosoma japonicum]"
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 108
        assert hsp.identities == 81
        assert hsp.positives == 83
        assert hsp.query == "MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.sbjct == "MADEEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.match == "MADEEVQALVVDNGSGMCKAG       ++  P               G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.score == 411.0
        assert hsp.expect == 1.00219e-47
        assert hsp.query_start == 20
        assert hsp.query_end == 343
        assert hsp.sbjct_start == 1
        assert hsp.sbjct_end == 108
        assert hsp.frame == (2, 0)

        alignment = record.alignments[5]
        assert alignment.hit_id == "gi|1535393712|emb|VDP83060.1|"
        assert alignment.accession == "VDP83060"
        assert alignment.length == 209
        assert alignment.title == "gi|1535393712|emb|VDP83060.1| unnamed protein product, partial [Echinostoma caproni]"
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 108
        assert hsp.identities == 80
        assert hsp.positives == 83
        assert hsp.query == "MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.sbjct == "MADDEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.match == "MAD+EVQALVVDNGSGMCKAG       ++  P               G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.score == 407.0
        assert hsp.expect == 1.16397e-47
        assert hsp.query_start == 20
        assert hsp.query_end == 343
        assert hsp.sbjct_start == 1
        assert hsp.sbjct_end == 108
        assert hsp.frame == (2, 0)

        alignment = record.alignments[6]
        assert alignment.hit_id == "gi|312773|emb|CAA50205.1|"
        assert alignment.accession == "CAA50205"
        assert alignment.length == 137
        assert alignment.title == "gi|312773|emb|CAA50205.1| actin, partial [Entamoeba histolytica]"
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 108
        assert hsp.identities == 78
        assert hsp.positives == 81
        assert hsp.query == "MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.sbjct == "MGDEEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHVSVMAGMGQKDAYVGDEAQSKRGILTLKYPIEHGIVNNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.match == "M DEEVQALVVDNGSGMCKAG       ++  P               G KD+YVGDEAQSKRGILTLKYPIEHGIV NWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.score == 401.0
        assert hsp.expect == 1.25869e-47
        assert hsp.query_start == 20
        assert hsp.query_end == 343
        assert hsp.sbjct_start == 1
        assert hsp.sbjct_end == 108
        assert hsp.frame == (2, 0)

        alignment = record.alignments[7]
        assert alignment.hit_id == "gi|1530341495|emb|VDN44756.1|"
        assert alignment.accession == "VDN44756"
        assert alignment.length == 145
        assert alignment.title == "gi|1530341495|emb|VDN44756.1| unnamed protein product, partial [Dibothriocephalus latus]"
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 108
        assert hsp.identities == 78
        assert hsp.positives == 82
        assert hsp.query == "MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.sbjct == "MGDEDVQALVIDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.match == "M DE+VQALV+DNGSGMCKAG       ++  P               G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.score == 400.0
        assert hsp.expect == 1.78336e-47
        assert hsp.query_start == 20
        assert hsp.query_end == 343
        assert hsp.sbjct_start == 1
        assert hsp.sbjct_end == 108
        assert hsp.frame == (2, 0)

        alignment = record.alignments[8]
        assert alignment.hit_id == "gi|1524877828|ref|XP_027046469.1|"
        assert alignment.accession == "XP_027046469"
        assert alignment.length == 122
        assert alignment.title == "gi|1524877828|ref|XP_027046469.1| actin-1, partial [Pocillopora damicornis]"
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 108
        assert hsp.identities == 78
        assert hsp.positives == 82
        assert hsp.query == "MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.sbjct == "MADEEVAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRIAPEEHPILLTE"
        assert hsp.match == "MADEEV ALVVDNGSGMCKAG       ++  P               G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELR+APEEHP+LLTE"
        assert hsp.score == 398.0
        assert hsp.expect == 1.93331e-47
        assert hsp.query_start == 20
        assert hsp.query_end == 343
        assert hsp.sbjct_start == 1
        assert hsp.sbjct_end == 108
        assert hsp.frame == (2, 0)

        alignment = record.alignments[9]
        assert alignment.hit_id == "gi|1524877860|ref|XP_027046487.1|"
        assert alignment.accession == "XP_027046487"
        assert alignment.length == 134
        assert alignment.title == "gi|1524877860|ref|XP_027046487.1| actin-1-like [Pocillopora damicornis]"
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 108
        assert hsp.identities == 79
        assert hsp.positives == 82
        assert hsp.query == "MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.sbjct == "MADEDVAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.match == "MADE+V ALVVDNGSGMCKAG       ++  P               G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.score == 399.0
        assert hsp.expect == 2.36088e-47
        assert hsp.query_start == 20
        assert hsp.query_end == 343
        assert hsp.sbjct_start == 1
        assert hsp.sbjct_end == 108
        assert hsp.frame == (2, 0)

        assert len(record.descriptions) == 10
        description = record.descriptions[0]
        assert description.title == "gi|1530504495|emb|VDM03167.1| unnamed protein product, partial [Schistocephalus solidus]"
        assert description.score == 408.0
        assert description.e == 8.11609e-49
        assert description.num_alignments == 1
        description = record.descriptions[1]
        assert description.title == "gi|510859078|gb|EPB74633.1| hypothetical protein ANCCEY_06263 [Ancylostoma ceylanicum]"
        assert description.score == 405.0
        assert description.e == 1.40046e-48
        assert description.num_alignments == 1
        description = record.descriptions[2]
        assert description.score == 413.0
        assert description.e == 4.40404e-48
        assert description.num_alignments == 1
        description = record.descriptions[3]
        assert description.title == "gi|449710331|gb|EMD49430.1| actin, putative, partial [Entamoeba histolytica KU27]"
        assert description.score == 401.0
        assert description.e == 9.0486e-48
        assert description.num_alignments == 1
        description = record.descriptions[4]
        assert description.title == "gi|257215766|emb|CAX83035.1| Actin-2, partial [Schistosoma japonicum]"
        assert description.score == 411.0
        assert description.e == 1.00219e-47
        assert description.num_alignments == 1
        description = record.descriptions[5]
        assert description.title == "gi|1535393712|emb|VDP83060.1| unnamed protein product, partial [Echinostoma caproni]"
        assert description.score == 407.0
        assert description.e == 1.16397e-47
        assert description.num_alignments == 1
        description = record.descriptions[6]
        assert description.title == "gi|312773|emb|CAA50205.1| actin, partial [Entamoeba histolytica]"
        assert description.score == 401.0
        assert description.e == 1.25869e-47
        assert description.num_alignments == 1
        description = record.descriptions[7]
        assert description.title == "gi|1530341495|emb|VDN44756.1| unnamed protein product, partial [Dibothriocephalus latus]"
        assert description.score == 400.0
        assert description.e == 1.78336e-47
        assert description.num_alignments == 1
        description = record.descriptions[8]
        assert description.title == "gi|1524877828|ref|XP_027046469.1| actin-1, partial [Pocillopora damicornis]"
        assert description.score == 398.0
        assert description.e == 1.93331e-47
        assert description.num_alignments == 1
        description = record.descriptions[9]
        assert description.title == "gi|1524877860|ref|XP_027046487.1| actin-1-like [Pocillopora damicornis]"
        assert description.score == 399.0
        assert description.e == 2.36088e-47
        assert description.num_alignments == 1

    def test_xml_2900_tblastn_001_v1(self):
        self._test_xml_2900_tblastn_001("xml_2900_tblastn_001.xml")

    def test_xml_2900_tblastn_001_v2(self):
        self._test_xml_2900_tblastn_001("xml_2900_tblastn_001_v2.xml")

    def _test_xml_2900_tblastn_001(self, filename):
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)
            record = next(records)
            with pytest.raises(StopIteration):
                next(records)

        assert record.application == "TBLASTN"
        assert record.version == "2.9.0+"
        assert record.reference == 'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.'
        assert record.database == "nr"
        assert record.date == ""
        assert record.query == "tim [Helicobacter acinonychis str. Sheeba]"
        assert record.query_letters == 234

        alignment = record.alignments[4]
        assert alignment.hit_id == "gi|1143706535|gb|CP018823.1|"
        assert alignment.accession == "CP018823"
        assert alignment.length == 1618480
        assert alignment.title == "gi|1143706535|gb|CP018823.1| Helicobacter pylori strain PMSS1 complete genome"
        assert len(alignment.hsps) == 1
        hsp = alignment.hsps[0]
        assert hsp.align_length == 234
        assert hsp.identities == 218
        assert hsp.positives == 223
        assert hsp.query == "MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQNAYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNFLKEKFDFFKDKKFKIVYCIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHGFLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL"
        assert hsp.sbjct == "MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHFDRVFVFPDFLGLLPNSFLHFTLGVQNAYPRDCGAFTGEITSKHLEELKIHTLLIGHSERRVLLKESPSFLKEKFDFFKDKNFKIVYCIGEDLTTREKGFKAVKEFLNEQLENIDLNYSNLIVAYEPIWAIGTKKSASLEDIYLTHGFLKQILNQKTPLLYGGSVNTQNAKEILGIDSVDGLLIGSASWELENFKTIISFL"
        assert hsp.match == "MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQH DRVFVFPDFLGLLPN+FLHFTLGVQNAYP+DCGAFTGEITSKHLEELKI+TLLIGHSERRVLLKESP+FLKEKFDFFKDK FKIVYCIGEDL TREKG  AVKEFLNEQLENIDL+Y NLIVAYEPIWAIGT KSASLEDIYLTHGFLKQ LNQK PLLYGGSVNTQNAKEILGIDSVDGLLIGS S ELENFKTIISFL"
        assert hsp.score == 1136.0
        assert hsp.bits == 442.195
        assert hsp.expect == 2.08707e-139
        assert hsp.query_start == 1
        assert hsp.query_end == 234
        assert hsp.sbjct_start == 190464
        assert hsp.sbjct_end == 191165
        assert hsp.frame == (0, 3)

        assert len(record.descriptions) == 10
        description = record.descriptions[4]
        assert description.title == "gi|1143706535|gb|CP018823.1| Helicobacter pylori strain PMSS1 complete genome"
        assert description.score == 1136.0
        assert description.e == 2.08707e-139
        assert description.num_alignments == 1

    def test_xml_2900_tblastx_001_v1(self):
        self._test_xml_2900_tblastx_001("xml_2900_tblastx_001.xml")

    def test_xml_2900_tblastx_001_v2(self):
        self._test_xml_2900_tblastx_001("xml_2900_tblastx_001_v2.xml")

    def _test_xml_2900_tblastx_001(self, filename):
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)
            record = next(records)

        assert record.application == "TBLASTX"
        assert record.version == "2.9.0+"
        assert record.reference == 'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.'
        assert record.database == "nr"
        assert record.date == ""
        assert record.query == "MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA, complete cds, mRNA sequence"
        assert record.query_letters == 365

        alignment = record.alignments[4]
        assert alignment.hit_id == "gi|1590279025|ref|XM_028307351.1|"
        assert alignment.accession == "XM_028307351"
        assert alignment.length == 1599
        assert alignment.title == "gi|1590279025|ref|XM_028307351.1| PREDICTED: Ostrinia furnacalis actin, cytoplasmic A3a (LOC114354791), mRNA"
        assert len(alignment.hsps) == 12
        hsp = alignment.hsps[0]
        assert hsp.align_length == 60
        assert hsp.identities == 59
        assert hsp.positives == 59
        assert hsp.query == "GSKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.sbjct == "GQKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.match == "G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        assert hsp.score == 325.0
        assert hsp.expect == 7.83523e-61
        assert hsp.query_start == 164
        assert hsp.query_end == 343
        assert hsp.sbjct_start == 389
        assert hsp.sbjct_end == 568
        assert hsp.frame == (2, 2)
        hsp = alignment.hsps[1]
        assert hsp.align_length == 36
        assert hsp.identities == 33
        assert hsp.positives == 33
        assert hsp.query == "GCAKLGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQK"
        assert hsp.sbjct == "GMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQK"
        assert hsp.match == "G  K GFAGDDAPRAVFPSIVGRPRHQGVMVGMGQK"
        assert hsp.score == 174.0
        assert hsp.expect == 7.83523e-61
        assert hsp.query_start == 66
        assert hsp.query_end == 173
        assert hsp.sbjct_start == 290
        assert hsp.sbjct_end == 397
        assert hsp.frame == (3, 2)
        hsp = alignment.hsps[2]
        assert hsp.align_length == 21
        assert hsp.identities == 19
        assert hsp.positives == 19
        assert hsp.query == "MADEEVQALVVDNGSGMCKAG"
        assert hsp.sbjct == "MCDEEVAALVVDNGSGMCKAG"
        assert hsp.match == "M DEEV ALVVDNGSGMCKAG"
        assert hsp.score == 97.0
        assert hsp.expect == 7.83523e-61
        assert hsp.query_start == 20
        assert hsp.query_end == 82
        assert hsp.sbjct_start == 245
        assert hsp.sbjct_end == 307
        assert hsp.frame == (2, 2)
        hsp = alignment.hsps[3]
        assert hsp.align_length == 61
        assert hsp.identities == 48
        assert hsp.positives == 52
        assert hsp.query == "SSVNRTGCSSGATRNSL*NV*CQIFSMSSQFVTIPCSIGYFSVRIPRFDCASSPT*LSFDP"
        assert hsp.sbjct == "ASVRRTGCSSGATRSSL*KV*CQIFSMSSQFVTIPCSMGYLSVRMPLLLWASSPT*ESFCP"
        assert hsp.match == "+SV RTGCSSGATR+SL* V*CQIFSMSSQFVTIPCS+GY SVR+P    ASSPT* SF P"
        assert hsp.score == 225.0
        assert hsp.expect == 3.26469e-39
        assert hsp.query_start == 163
        assert hsp.query_end == 345
        assert hsp.sbjct_start == 388
        assert hsp.sbjct_end == 570
        assert hsp.frame == (-3, -1)
        hsp = alignment.hsps[4]
        assert hsp.align_length == 36
        assert hsp.identities == 24
        assert hsp.positives == 28
        assert hsp.query == "F*PIPTITP*CRGRPTMEGNTALGASSPANPSFAHP"
        assert hsp.sbjct == "FCPMPTITPWWRGRPTIDGNTARGASSPAKPALHIP"
        assert hsp.match == "F P+PTITP  RGRPT++GNTA GASSPA P+   P"
        assert hsp.score == 121.0
        assert hsp.expect == 3.26469e-39
        assert hsp.query_start == 65
        assert hsp.query_end == 172
        assert hsp.sbjct_start == 289
        assert hsp.sbjct_end == 396
        assert hsp.frame == (-2, -1)
        hsp = alignment.hsps[5]
        assert hsp.align_length == 26
        assert hsp.identities == 18
        assert hsp.positives == 22
        assert hsp.query == "PALHIPDPLSTTRA*TSSSAMIIFQL"
        assert hsp.sbjct == "PALHIPDPLSTTNAATSSSHILVYLL"
        assert hsp.match == "PALHIPDPLSTT A TSSS ++++ L"
        assert hsp.score == 91.0
        assert hsp.expect == 3.26469e-39
        assert hsp.query_start == 4
        assert hsp.query_end == 81
        assert hsp.sbjct_start == 229
        assert hsp.sbjct_end == 306
        assert hsp.frame == (-3, -1)
        hsp = alignment.hsps[6]
        assert hsp.align_length == 57
        assert hsp.identities == 34
        assert hsp.positives == 43
        assert hsp.query == "GQQDRVFFWSHTQFIVECVMPDLLHVIPVRHNTVFDWVFQCENTTFRLCFITDVAVF"
        assert hsp.sbjct == "GKKDWVFLGSDTEFIVEGVMPDLLHVIPVCDDSVFDGVFECEDASLALGLISDVGVF"
        assert hsp.match == "G++D VF  S T+FIVE VMPDLLHVIPV  ++VFD VF+CE+ +  L  I+DV VF"
        assert hsp.score == 167.0
        assert hsp.expect == 2.27039e-16
        assert hsp.query_start == 170
        assert hsp.query_end == 340
        assert hsp.sbjct_start == 395
        assert hsp.sbjct_end == 565
        assert hsp.frame == (-2, -3)
        hsp = alignment.hsps[7]
        assert hsp.align_length == 32
        assert hsp.identities == 13
        assert hsp.positives == 23
        assert hsp.query == "LTHTNHHTLMSRSSNDGREYCSWCIITSESQL"
        assert hsp.sbjct == "LSHANHHTLVARAAHDRWEHGAGRVIARETGL"
        assert hsp.match == "L+H NHHTL++R+++D  E+ +  +I  E+ L"
        assert hsp.score == 69.0
        assert hsp.expect == 2.27039e-16
        assert hsp.query_start == 75
        assert hsp.query_end == 170
        assert hsp.sbjct_start == 299
        assert hsp.sbjct_end == 394
        assert hsp.frame == (-1, -3)
        hsp = alignment.hsps[8]
        assert hsp.align_length == 64
        assert hsp.identities == 36
        assert hsp.positives == 44
        assert hsp.query == "KRQLRR**STIETWYSHTEIPNRTRYCDELG*HGEDLASHILQ*IACGSRRTPCPVDRAPLYPK"
        assert hsp.sbjct == "KRLLRRR*GPEQERHPHTQIPHRTRNRHKLG*HGEDLASHLLQ*TPCRSRGTPSPSYRSPPEPQ"
        assert hsp.match == "KR LRR *   +  + HT+IP+RTR   +LG*HGEDLASH+LQ*  C SR TP P  R+P  P+"
        assert hsp.score == 160.0
        assert hsp.expect == 2.96635e-13
        assert hsp.query_start == 169
        assert hsp.query_end == 360
        assert hsp.sbjct_start == 394
        assert hsp.sbjct_end == 585
        assert hsp.frame == (1, 1)
        hsp = alignment.hsps[9]
        assert hsp.align_length == 20
        assert hsp.identities == 8
        assert hsp.positives == 12
        assert hsp.query == "IPFHRWTTSTSRCDGWYGSK"
        assert hsp.sbjct == "VPIDRGPPAPPGCDGWHGTK"
        assert hsp.match == "+P  R   +   CDGW+G+K"
        assert hsp.score == 53.0
        assert hsp.expect == 2.96635e-13
        assert hsp.query_start == 113
        assert hsp.query_end == 172
        assert hsp.sbjct_start == 337
        assert hsp.sbjct_end == 396
        assert hsp.frame == (2, 1)
        hsp = alignment.hsps[10]
        assert hsp.align_length == 48
        assert hsp.identities == 29
        assert hsp.positives == 33
        assert hsp.query == "FGYNGARSTGQGVLLEPHAIHCRMCDARSSPCHPSSSQYRVRLGISV*"
        assert hsp.sbjct == "WGSGGLR*EGLGVPRERHGVHCRRCDARSSPCHPSL*RFRVRWGI*V*"
        assert hsp.match == "+G  G R  G GV  E H +HCR CDARSSPCHPS  ++RVR GI V*"
        assert hsp.score == 142.0
        assert hsp.expect == 2.61074e-08
        assert hsp.query_start == 216
        assert hsp.query_end == 359
        assert hsp.sbjct_start == 441
        assert hsp.sbjct_end == 584
        assert hsp.frame == (-1, -2)
        hsp = alignment.hsps[11]
        assert hsp.align_length == 58
        assert hsp.identities == 33
        assert hsp.positives == 37
        assert hsp.query == "QKTATSVMKHNRNVVFSH*NTQSNTVL*RTGMTWRRSGITHSTMNCVWLQKNTLSC*P"
        assert hsp.sbjct == "KKTPTSEMRPRAREASSHSNTPSNTESSQTGMTWRRSGITPSTMNSVSLPRNTQSFLP"
        assert hsp.match == "+KT TS M+       SH NT SNT   +TGMTWRRSGIT STMN V L +NT S  P"
        assert hsp.score == 141.0
        assert hsp.expect == 3.58672e-08
        assert hsp.query_start == 168
        assert hsp.query_end == 341
        assert hsp.sbjct_start == 393
        assert hsp.sbjct_end == 566
        assert hsp.frame == (3, 3)
        assert len(record.descriptions) == 10
        description = record.descriptions[4]
        assert (description.title == "gi|1590279025|ref|XM_028307351.1| PREDICTED: Ostrinia furnacalis actin, cytoplasmic A3a "
            "(LOC114354791), mRNA")
        assert description.score == 325.0
        assert description.e == 7.83523e-61
        assert description.num_alignments == 12

    def test_create_view_warning(self):
        """Test that CREATE_VIEW between Hit tags raises a BiopythonParserWarning."""
        # Minimal BLAST XML with CREATE_VIEW text between Hit tags
        xml_data = b"""\
<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastp</BlastOutput_program>
  <BlastOutput_version>BLASTP 2.2.12 [Aug-07-2005]</BlastOutput_version>
  <BlastOutput_reference>Altschul et al.</BlastOutput_reference>
  <BlastOutput_db>nr</BlastOutput_db>
  <BlastOutput_query-ID>test_query</BlastOutput_query-ID>
  <BlastOutput_query-def>test query</BlastOutput_query-def>
  <BlastOutput_query-len>100</BlastOutput_query-len>
  <BlastOutput_param>
    <Parameters>
      <Parameters_matrix>BLOSUM62</Parameters_matrix>
      <Parameters_expect>10</Parameters_expect>
      <Parameters_gap-open>11</Parameters_gap-open>
      <Parameters_gap-extend>1</Parameters_gap-extend>
      <Parameters_filter>L;</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_hits>
        CREATE_VIEW
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>test_hit</Hit_id>
          <Hit_def>test hit definition</Hit_def>
          <Hit_accession>TEST123</Hit_accession>
          <Hit_len>100</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>50.0</Hsp_bit-score>
              <Hsp_score>100</Hsp_score>
              <Hsp_evalue>1e-10</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>100</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>100</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>100</Hsp_identity>
              <Hsp_positive>100</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>100</Hsp_align-len>
              <Hsp_qseq>TESTSEQUENCE</Hsp_qseq>
              <Hsp_hseq>TESTSEQUENCE</Hsp_hseq>
              <Hsp_midline>TESTSEQUENCE</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
      <Iteration_stat>
        <Statistics>
          <Statistics_db-num>1000</Statistics_db-num>
          <Statistics_db-len>100000</Statistics_db-len>
          <Statistics_hsp-len>10</Statistics_hsp-len>
          <Statistics_eff-space>1000000</Statistics_eff-space>
          <Statistics_kappa>0.041</Statistics_kappa>
          <Statistics_lambda>0.267</Statistics_lambda>
          <Statistics_entropy>0.14</Statistics_entropy>
        </Statistics>
      </Iteration_stat>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", BiopythonParserWarning)
            handle = io.BytesIO(xml_data)
            records = list(NCBIXML.parse(handle))
            assert len(records) == 1
            # Check that the CREATE_VIEW warning was raised
            create_view_warnings = [
                w
                for w in caught
                if issubclass(w.category, BiopythonParserWarning)
                and "CREATE_VIEW" in str(w.message)
            ]
            assert len(create_view_warnings) == 1
            assert "NCBIXML: Ignored:" in str(create_view_warnings[0].message)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
