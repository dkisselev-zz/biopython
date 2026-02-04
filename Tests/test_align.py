# Copyright 2000-2001 by Brad Chapman.  All rights reserved.
# Revisions copyright 2007-2003 by Peter Cock. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""Test alignment stuff.

Right now we've got tests for:

- Reading and Writing clustal format
- Reading and Writing fasta format
- Converting between formats

"""

# standard library
import os
import unittest
import pytest
import warnings
from io import StringIO

from Bio import Align
from Bio import AlignIO

# biopython
from Bio import BiopythonDeprecationWarning
from Bio import motifs
from Bio.Align import AlignInfo
from Bio.Align import Alignment
from Bio.Align import MultipleSeqAlignment
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np


class TestBasics(unittest.TestCase):
    def test_empty_alignment(self):
        """Very simple tests on an empty alignment."""
        alignment = MultipleSeqAlignment([])
        assert alignment.get_alignment_length() == 0
        assert len(alignment) == 0
        alignment = alignment.alignment  # new-style Alignment object
        assert alignment.length == 0
        assert len(alignment) == 0

    def test_basic_alignment(self):
        """Basic tests on a simple alignment of three sequences."""
        msa = MultipleSeqAlignment([])
        letters = "AbcDefGhiJklMnoPqrStuVwxYz"
        msa.append(SeqRecord(Seq(letters), id="mixed"))
        msa.append(SeqRecord(Seq(letters.lower()), id="lower"))
        msa.append(SeqRecord(Seq(letters.upper()), id="upper"))
        msa.append(SeqRecord(Seq(letters), id="duplicate"))
        del msa[3]
        assert msa.get_alignment_length() == 26
        assert len(msa) == 3
        assert msa[0].seq == letters
        assert msa[1].seq == letters.lower()
        assert msa[2].seq == letters.upper()
        assert msa[0].id == "mixed"
        assert msa[1].id == "lower"
        assert msa[2].id == "upper"
        for col, letter in enumerate(letters):
            assert msa[:, col] == letter + letter.lower() + letter.upper()
        # Check row extractions:
        assert msa[0].id == "mixed"
        assert msa[-1].id == "upper"
        # Check sub-alignment extraction by row slicing:
        assert isinstance(msa[::-1], MultipleSeqAlignment)
        assert msa[::-1][0].id == "upper"
        assert msa[::-1][2].id == "mixed"
        # create a new-style Alignment object
        alignment = msa.alignment
        assert alignment.shape == (3, 26)
        assert len(alignment) == 3
        assert alignment.sequences[0].seq == letters
        assert alignment.sequences[1].seq == letters.lower()
        assert alignment.sequences[2].seq == letters.upper()
        assert alignment.sequences[0].id == "mixed"
        assert alignment.sequences[1].id == "lower"
        assert alignment.sequences[2].id == "upper"
        for col, letter in enumerate(letters):
            assert alignment[:, col] == letter + letter.lower() + letter.upper()
        # Check row extractions:
        assert alignment[0] == letters
        assert alignment[-1] == letters.upper()
        # Check sub-alignment extraction by row slicing:
        assert isinstance(alignment[::-1], Alignment)
        assert alignment[::-1].sequences[0].id == "upper"
        assert alignment[::-1].sequences[2].id == "mixed"


class TestReading(unittest.TestCase):
    def test_read_clustal1(self):
        """Parse an alignment file and get an alignment object."""
        opuntia_clustal_header = """\
CLUSTAL X (1.81) multiple sequence alignment


"""
        opuntia_clustal_body = """\
gi|6273285|gb|AF191659.1|AF191      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
gi|6273284|gb|AF191658.1|AF191      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
gi|6273287|gb|AF191661.1|AF191      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
gi|6273286|gb|AF191660.1|AF191      TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
gi|6273290|gb|AF191664.1|AF191      TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
gi|6273289|gb|AF191663.1|AF191      TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
gi|6273291|gb|AF191665.1|AF191      TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
                                    ******* **** *************************************

gi|6273285|gb|AF191659.1|AF191      TATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATA
gi|6273284|gb|AF191658.1|AF191      TATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATA
gi|6273287|gb|AF191661.1|AF191      TATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATA
gi|6273286|gb|AF191660.1|AF191      TATATA----------ATATATTTATAATTTCCTTATATATCCAAATATA
gi|6273290|gb|AF191664.1|AF191      TATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
gi|6273289|gb|AF191663.1|AF191      TATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
gi|6273291|gb|AF191665.1|AF191      TATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATA
                                    ******          ********  **** ********* *********

gi|6273285|gb|AF191659.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGT
gi|6273284|gb|AF191658.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGT
gi|6273287|gb|AF191661.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGT
gi|6273286|gb|AF191660.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGT
gi|6273290|gb|AF191664.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGT
gi|6273289|gb|AF191663.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTAT
gi|6273291|gb|AF191665.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGT
                                    ************************************ *********** *

gi|6273285|gb|AF191659.1|AF191      ACCAGA
gi|6273284|gb|AF191658.1|AF191      ACCAGA
gi|6273287|gb|AF191661.1|AF191      ACCAGA
gi|6273286|gb|AF191660.1|AF191      ACCAGA
gi|6273290|gb|AF191664.1|AF191      ACCAGA
gi|6273289|gb|AF191663.1|AF191      ACCAGA
gi|6273291|gb|AF191665.1|AF191      ACCAGA
                                    ******


"""  # noqa : W291
        opuntia_fasta = """\
>gi|6273285|gb|AF191659.1|AF191
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGA
TGAATATCAAAGAATCCATTGATTTAGTGTACCAGA
>gi|6273284|gb|AF191658.1|AF191
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--
------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGA
TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
>gi|6273287|gb|AF191661.1|AF191
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
------ATATATTTCAAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
>gi|6273286|gb|AF191660.1|AF191
TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
------ATATATTTATAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
>gi|6273290|gb|AF191664.1|AF191
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA
------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
>gi|6273289|gb|AF191663.1|AF191
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA
------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
TGAATATCAAAGAATCTATTGATTTAGTATACCAGA
>gi|6273291|gb|AF191665.1|AF191
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA
TATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
"""
        opuntia_fasta_oneline = """\
>gi|6273285|gb|AF191659.1|AF191
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCAGA
>gi|6273284|gb|AF191658.1|AF191
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
>gi|6273287|gb|AF191661.1|AF191
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
>gi|6273286|gb|AF191660.1|AF191
TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTATAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
>gi|6273290|gb|AF191664.1|AF191
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
>gi|6273289|gb|AF191663.1|AF191
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTATACCAGA
>gi|6273291|gb|AF191665.1|AF191
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
"""
        opuntia_fasta_oneline_with_description = """\
>gi|6273285|gb|AF191659.1|AF191 gi|6273285|gb|AF191659.1|AF191
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCAGA
>gi|6273284|gb|AF191658.1|AF191 gi|6273284|gb|AF191658.1|AF191
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
>gi|6273287|gb|AF191661.1|AF191 gi|6273287|gb|AF191661.1|AF191
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
>gi|6273286|gb|AF191660.1|AF191 gi|6273286|gb|AF191660.1|AF191
TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTATAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
>gi|6273290|gb|AF191664.1|AF191 gi|6273290|gb|AF191664.1|AF191
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
>gi|6273289|gb|AF191663.1|AF191 gi|6273289|gb|AF191663.1|AF191
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTATACCAGA
>gi|6273291|gb|AF191665.1|AF191 gi|6273291|gb|AF191665.1|AF191
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
"""
        path = os.path.join(os.getcwd(), "Clustalw", "opuntia.aln")
        msa = AlignIO.read(path, "clustal")
        opuntia_clustal = opuntia_clustal_header + opuntia_clustal_body
        assert format(msa, "clustal") == opuntia_clustal
        assert format(msa, "fasta") == opuntia_fasta
        # create a new-style Alignment object
        alignment = msa.alignment
        assert format(alignment, "clustal") == opuntia_clustal_body
        # New-style Alignment objects generate FASTA format with the sequence
        # on one line. Also, the clustal parser in Bio.AlignIO generates
        # SeqRecords with an (identical) ID and a description; the clustal
        # parser in Bio.Align generates SeqRecords with an ID only.
        assert format(alignment, "fasta") == opuntia_fasta_oneline_with_description
        alignment = Align.read(path, "clustal")
        assert format(alignment, "fasta") == opuntia_fasta_oneline

    def test_read_clustal2(self):
        """Parse an alignment file and get an alignment object."""
        clustalw_clustal_header = """\
CLUSTAL X (1.81) multiple sequence alignment


"""
        clustalw_clustal_body = """\
gi|4959044|gb|AAD34209.1|AF069      MENSDSNDKGSDQSAAQRRSQMDRLDREEAFYQFVNNLSEEDYRLMRDNN
gi|671626|emb|CAA85685.1|           ---------MSPQTETKASVGFKAGVKEYKLTYYTPEYETKDTDILAAFR
                                              * *: ::    :.   :*  :  :. : . :*  ::   .

gi|4959044|gb|AAD34209.1|AF069      LLGTPGESTEEELLRRLQQIKEGPPPQSPDENRAGESSDDVTNSDSIIDW
gi|671626|emb|CAA85685.1|           VTPQPG-----------------VPPEEAGAAVAAESSTGT---------
                                    :   **                  **:...   *.*** ..         

gi|4959044|gb|AAD34209.1|AF069      LNSVRQTGNTTRSRQRGNQSWRAVSRTNPNSGDFRFSLEINVNRNNGSQT
gi|671626|emb|CAA85685.1|           WTTVWTDGLTSLDRYKG-----RCYHIEPVPG------------------
                                     .:*   * *: .* :*        : :* .*                  

gi|4959044|gb|AAD34209.1|AF069      SENESEPSTRRLSVENMESSSQRQMENSASESASARPSRAERNSTEAVTE
gi|671626|emb|CAA85685.1|           -EKDQCICYVAYPLDLFEEGSVTNMFTSIVGNVFGFKALRALRLEDLRIP
                                     *::.  .    .:: :*..*  :* .*   .. .  :    .  :    

gi|4959044|gb|AAD34209.1|AF069      VPTTRAQRRARSRSPEHRRTRARAERSMSPLQPTSEIPRRAPTLEQSSEN
gi|671626|emb|CAA85685.1|           VAYVKTFQGPPHGIQVERDKLNKYGRPLLGCTIKPKLGLSAKNYGRAVYE
                                    *. .:: : .      .* .  :  *.:     ..::   * .  ::  :

gi|4959044|gb|AAD34209.1|AF069      EPEGSSRTRHHVTLRQQISGPELLGRGLFAASGSRNPSQGTSSSDTGSNS
gi|671626|emb|CAA85685.1|           CLRGGLDFTKDDENVNSQPFMRWRDRFLFCAEAIYKAQAETGEIKGHYLN
                                      .*.    :.    :. .  .  .* **.*..  :..  *.. .    .

gi|4959044|gb|AAD34209.1|AF069      ESSGSGQRPPTIVLDLQVRRVRPGEYRQRDSIASRTRSRSQAPNNTVTYE
gi|671626|emb|CAA85685.1|           ATAG-----------------------TCEEMIKRAIFARELGVPIVMHD
                                     ::*                         :.: .*:    :     * ::

gi|4959044|gb|AAD34209.1|AF069      SERGGFRRTFSRSERAGVRTYVSTIRIPIRRILNTGLSETTSVAIQTMLR
gi|671626|emb|CAA85685.1|           YLTGGFTANTSLAHYCRDNGLLLHIHRAMHAVIDRQKNHGMHFRVLAKAL
                                       ***  . * :. .  .  :  *: .:: :::   ..   . : :   

gi|4959044|gb|AAD34209.1|AF069      QIMTGFGELSYFMYSDSDSEPSASVSSRNVERVESRNGRGSSGGGNSSGS
gi|671626|emb|CAA85685.1|           RLSGGDHIHSGTVVGKLEGERDITLGFVDLLRDDFIEKDRSRGIYFTQDW
                                    ::  *    *  : .. :.* . ::.  :: * :  :   * *   :.. 

gi|4959044|gb|AAD34209.1|AF069      SSSSSPSPSSSGESSESSSKMFEGSSEGGSSGPSRKDGRHRAPVTFDESG
gi|671626|emb|CAA85685.1|           VSLPGVIPVASG-----------------------------GIHVWHMPA
                                     * ..  * :**                             .  .:. ..

gi|4959044|gb|AAD34209.1|AF069      SLPFFSLAQFFLLNEDDEDQPRGLTKEQIDNLAMRSFGENDALKTCSVCI
gi|671626|emb|CAA85685.1|           LTEIFGDDSVLQFGGGTLGHPWGNAPGAVANRVA-----------VEACV
                                       :*.  ..: :. .  .:* * :   : * .             ..*:

gi|4959044|gb|AAD34209.1|AF069      TEYTEGDKLRKLPCSHEFHVHCIDRWLSE-NSTCPICRRAVLSSGNRESV
gi|671626|emb|CAA85685.1|           KARNEG---RDLAAEGNAIIREACKWSPELAAACEVWKEIKFEFPAMD--
                                    .  .**   *.*... :  ::   :* .*  ::* : :.  :.    :  

gi|4959044|gb|AAD34209.1|AF069      V
gi|671626|emb|CAA85685.1|           -
                                     


"""  # noqa : W291

        path = os.path.join(os.curdir, "Clustalw", "clustalw.aln")
        msa = AlignIO.read(path, "clustal")
        clustalw_clustal = clustalw_clustal_header + clustalw_clustal_body
        assert format(msa, "clustal") == clustalw_clustal
        # create a new-style Alignment object
        alignment = msa.alignment
        assert format(alignment, "clustal") == clustalw_clustal_body

    def test_read_write_clustal(self):
        """Test the base alignment stuff."""
        path = os.path.join(os.getcwd(), "Clustalw", "opuntia.aln")
        msa = AlignIO.read(path, "clustal")
        assert len(msa) == 7
        seq_record = msa[0]
        assert seq_record.description == "gi|6273285|gb|AF191659.1|AF191"
        assert seq_record.seq == Seq(
                "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCAGA"
            )
        seq_record = msa[1]
        assert seq_record.description == "gi|6273284|gb|AF191658.1|AF191"
        assert seq_record.seq == "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA"
        seq_record = msa[2]
        assert seq_record.description == "gi|6273287|gb|AF191661.1|AF191"
        assert seq_record.seq == "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA"
        seq_record = msa[3]
        assert seq_record.description == "gi|6273286|gb|AF191660.1|AF191"
        assert seq_record.seq == "TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTATAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA"
        seq_record = msa[4]
        assert seq_record.description == "gi|6273290|gb|AF191664.1|AF191"
        assert seq_record.seq == "TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA"
        seq_record = msa[5]
        assert seq_record.description == "gi|6273289|gb|AF191663.1|AF191"
        assert seq_record.seq == "TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTATACCAGA"
        seq_record = msa[6]
        assert seq_record.description == "gi|6273291|gb|AF191665.1|AF191"
        assert seq_record.seq == "TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA"
        assert msa.get_alignment_length() == 156
        alignment = msa.alignment
        dictionary = alignment.substitutions
        assert len(dictionary) == 4
        assert dictionary.shape == (4, 4)
        assert len(dictionary.keys()) == 16
        assert dictionary[("A", "A")] == pytest.approx(1395, abs=5e-8)
        assert dictionary[("A", "C")] == pytest.approx(3, abs=5e-8)
        assert dictionary[("A", "G")] == pytest.approx(13, abs=5e-8)
        assert dictionary[("A", "T")] == pytest.approx(6, abs=5e-8)
        assert dictionary[("C", "A")] == pytest.approx(3, abs=5e-8)
        assert dictionary[("C", "C")] == pytest.approx(271, abs=5e-8)
        assert dictionary[("C", "G")] == pytest.approx(0, abs=5e-8)
        assert dictionary[("C", "T")] == pytest.approx(16, abs=5e-8)
        assert dictionary[("G", "A")] == pytest.approx(5, abs=5e-8)
        assert dictionary[("G", "C")] == pytest.approx(0, abs=5e-8)
        assert dictionary[("G", "G")] == pytest.approx(480, abs=5e-8)
        assert dictionary[("G", "T")] == pytest.approx(0, abs=5e-8)
        assert dictionary[("T", "A")] == pytest.approx(6, abs=5e-8)
        assert dictionary[("T", "C")] == pytest.approx(12, abs=5e-8)
        assert dictionary[("T", "G")] == pytest.approx(0, abs=5e-8)
        assert dictionary[("T", "T")] == pytest.approx(874, abs=5e-8)

        motif = motifs.Motif("ACGT", alignment)
        counts = motif.counts
        assert counts.calculate_consensus(identity=0.7) == "TATACATTAAAGNAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATATATATAATATATTTCAAATTNCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA"
        assert str(counts) == """\
        0      1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16     17     18     19     20     21     22     23     24     25     26     27     28     29     30     31     32     33     34     35     36     37     38     39     40     41     42     43     44     45     46     47     48     49     50     51     52     53     54     55     56     57     58     59     60     61     62     63     64     65     66     67     68     69     70     71     72     73     74     75     76     77     78     79     80     81     82     83     84     85     86     87     88     89     90     91     92     93     94     95     96     97     98     99    100    101    102    103    104    105    106    107    108    109    110    111    112    113    114    115    116    117    118    119    120    121    122    123    124    125    126    127    128    129    130    131    132    133    134    135    136    137    138    139    140    141    142    143    144    145    146    147    148    149    150    151    152    153    154    155
A:   0.00   7.00   0.00   7.00   0.00   7.00   0.00   1.00   7.00   7.00   7.00   0.00   4.00   7.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   7.00   7.00   7.00   0.00   0.00   0.00   7.00   7.00   7.00   0.00   0.00   0.00   0.00   7.00   7.00   7.00   0.00   7.00   7.00   7.00   0.00   7.00   7.00   0.00   7.00   0.00   7.00   0.00   7.00   0.00   4.00   0.00   3.00   0.00   1.00   0.00   1.00   0.00   1.00   7.00   0.00   7.00   0.00   7.00   0.00   0.00   0.00   1.00   6.00   7.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   7.00   0.00   7.00   0.00   0.00   0.00   7.00   7.00   7.00   0.00   7.00   0.00   7.00   7.00   7.00   7.00   7.00   0.00   7.00   0.00   0.00   0.00   7.00   7.00   0.00   7.00   7.00   7.00   0.00   0.00   7.00   0.00   7.00   0.00   0.00   7.00   7.00   0.00   7.00   0.00   0.00   7.00   7.00   7.00   0.00   7.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   1.00   0.00   7.00   0.00   0.00   7.00   0.00   7.00
C:   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   6.00   0.00   0.00   0.00   0.00   0.00   3.00   7.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   2.00   7.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   1.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   7.00   0.00   0.00   0.00
G:   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   3.00   0.00   7.00   7.00   7.00   7.00   7.00   0.00   0.00   7.00   0.00   7.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   7.00   0.00   0.00   0.00   7.00   7.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   6.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00
T:   7.00   0.00   7.00   0.00   0.00   0.00   7.00   6.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   7.00   0.00   7.00   0.00   4.00   0.00   3.00   0.00   1.00   0.00   1.00   0.00   1.00   0.00   0.00   7.00   0.00   7.00   0.00   7.00   7.00   7.00   0.00   1.00   0.00   0.00   7.00   7.00   4.00   0.00   0.00   7.00   7.00   0.00   7.00   0.00   7.00   0.00   5.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   7.00   0.00   7.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   6.00   0.00   7.00   7.00   0.00   0.00   7.00   7.00   7.00   0.00   0.00   7.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00
"""

        alignment = msa.alignment
        motif = motifs.Motif("ACGT", alignment)
        counts = motif.counts
        assert str(counts) == """\
        0      1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16     17     18     19     20     21     22     23     24     25     26     27     28     29     30     31     32     33     34     35     36     37     38     39     40     41     42     43     44     45     46     47     48     49     50     51     52     53     54     55     56     57     58     59     60     61     62     63     64     65     66     67     68     69     70     71     72     73     74     75     76     77     78     79     80     81     82     83     84     85     86     87     88     89     90     91     92     93     94     95     96     97     98     99    100    101    102    103    104    105    106    107    108    109    110    111    112    113    114    115    116    117    118    119    120    121    122    123    124    125    126    127    128    129    130    131    132    133    134    135    136    137    138    139    140    141    142    143    144    145    146    147    148    149    150    151    152    153    154    155
A:   0.00   7.00   0.00   7.00   0.00   7.00   0.00   1.00   7.00   7.00   7.00   0.00   4.00   7.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   7.00   7.00   7.00   0.00   0.00   0.00   7.00   7.00   7.00   0.00   0.00   0.00   0.00   7.00   7.00   7.00   0.00   7.00   7.00   7.00   0.00   7.00   7.00   0.00   7.00   0.00   7.00   0.00   7.00   0.00   4.00   0.00   3.00   0.00   1.00   0.00   1.00   0.00   1.00   7.00   0.00   7.00   0.00   7.00   0.00   0.00   0.00   1.00   6.00   7.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   7.00   0.00   7.00   0.00   0.00   0.00   7.00   7.00   7.00   0.00   7.00   0.00   7.00   7.00   7.00   7.00   7.00   0.00   7.00   0.00   0.00   0.00   7.00   7.00   0.00   7.00   7.00   7.00   0.00   0.00   7.00   0.00   7.00   0.00   0.00   7.00   7.00   0.00   7.00   0.00   0.00   7.00   7.00   7.00   0.00   7.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   1.00   0.00   7.00   0.00   0.00   7.00   0.00   7.00
C:   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   6.00   0.00   0.00   0.00   0.00   0.00   3.00   7.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   2.00   7.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   1.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   7.00   0.00   0.00   0.00
G:   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   3.00   0.00   7.00   7.00   7.00   7.00   7.00   0.00   0.00   7.00   0.00   7.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   7.00   0.00   0.00   0.00   7.00   7.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   6.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00
T:   7.00   0.00   7.00   0.00   0.00   0.00   7.00   6.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   7.00   0.00   7.00   0.00   4.00   0.00   3.00   0.00   1.00   0.00   1.00   0.00   1.00   0.00   0.00   7.00   0.00   7.00   0.00   7.00   7.00   7.00   0.00   1.00   0.00   0.00   7.00   7.00   4.00   0.00   0.00   7.00   7.00   0.00   7.00   0.00   7.00   0.00   5.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   7.00   0.00   7.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   6.00   0.00   7.00   7.00   0.00   0.00   7.00   7.00   7.00   0.00   0.00   7.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00
"""

        with pytest.warns(BiopythonDeprecationWarning):
            align_info = AlignInfo.SummaryInfo(msa)
        assert align_info.get_column(1) == "AAAAAAA"
        assert align_info.get_column(7) == "TTTATTT"

        alignment = msa.alignment
        motif = motifs.Motif("ACGT", alignment)
        counts = motif.counts
        assert str(counts) == """\
        0      1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16     17     18     19     20     21     22     23     24     25     26     27     28     29     30     31     32     33     34     35     36     37     38     39     40     41     42     43     44     45     46     47     48     49     50     51     52     53     54     55     56     57     58     59     60     61     62     63     64     65     66     67     68     69     70     71     72     73     74     75     76     77     78     79     80     81     82     83     84     85     86     87     88     89     90     91     92     93     94     95     96     97     98     99    100    101    102    103    104    105    106    107    108    109    110    111    112    113    114    115    116    117    118    119    120    121    122    123    124    125    126    127    128    129    130    131    132    133    134    135    136    137    138    139    140    141    142    143    144    145    146    147    148    149    150    151    152    153    154    155
A:   0.00   7.00   0.00   7.00   0.00   7.00   0.00   1.00   7.00   7.00   7.00   0.00   4.00   7.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   7.00   7.00   7.00   0.00   0.00   0.00   7.00   7.00   7.00   0.00   0.00   0.00   0.00   7.00   7.00   7.00   0.00   7.00   7.00   7.00   0.00   7.00   7.00   0.00   7.00   0.00   7.00   0.00   7.00   0.00   4.00   0.00   3.00   0.00   1.00   0.00   1.00   0.00   1.00   7.00   0.00   7.00   0.00   7.00   0.00   0.00   0.00   1.00   6.00   7.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   7.00   0.00   7.00   0.00   0.00   0.00   7.00   7.00   7.00   0.00   7.00   0.00   7.00   7.00   7.00   7.00   7.00   0.00   7.00   0.00   0.00   0.00   7.00   7.00   0.00   7.00   7.00   7.00   0.00   0.00   7.00   0.00   7.00   0.00   0.00   7.00   7.00   0.00   7.00   0.00   0.00   7.00   7.00   7.00   0.00   7.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   1.00   0.00   7.00   0.00   0.00   7.00   0.00   7.00
C:   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   6.00   0.00   0.00   0.00   0.00   0.00   3.00   7.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   2.00   7.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   1.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   7.00   0.00   0.00   0.00
G:   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   3.00   0.00   7.00   7.00   7.00   7.00   7.00   0.00   0.00   7.00   0.00   7.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   7.00   0.00   0.00   0.00   7.00   7.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   6.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00
T:   7.00   0.00   7.00   0.00   0.00   0.00   7.00   6.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   7.00   0.00   7.00   0.00   4.00   0.00   3.00   0.00   1.00   0.00   1.00   0.00   1.00   0.00   0.00   7.00   0.00   7.00   0.00   7.00   7.00   7.00   0.00   1.00   0.00   0.00   7.00   7.00   4.00   0.00   0.00   7.00   7.00   0.00   7.00   0.00   7.00   0.00   5.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   7.00   0.00   7.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   6.00   0.00   7.00   7.00   0.00   0.00   7.00   7.00   7.00   0.00   0.00   7.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00
"""
        value = sum(motif[5:50].relative_entropy)
        assert value == pytest.approx(88.42309908538343, abs=5e-8)  # Alignment
        relative_entropy = motif.relative_entropy
        value = sum(relative_entropy)
        assert value == pytest.approx(306.2080592664532, abs=5e-8)  # Alignment
        assert relative_entropy[0] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[1] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[2] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[3] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[4] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[5] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[6] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[7] == pytest.approx(1.4083272214176723, abs=5e-8)
        assert relative_entropy[8] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[9] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[10] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[11] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[12] == pytest.approx(1.0147718639657484, abs=5e-8)
        assert relative_entropy[13] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[14] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[15] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[16] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[17] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[18] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[19] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[20] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[21] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[22] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[23] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[24] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[25] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[26] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[27] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[28] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[29] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[30] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[31] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[32] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[33] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[34] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[35] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[36] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[37] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[38] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[39] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[40] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[41] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[42] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[43] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[44] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[45] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[46] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[47] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[48] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[49] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[50] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[51] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[52] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[53] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[54] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[55] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[56] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[57] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[58] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[59] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[60] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[61] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[62] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[63] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[64] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[65] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[66] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[67] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[68] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[69] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[70] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[71] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[72] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[73] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[74] == pytest.approx(1.4083272214176723, abs=5e-8)
        assert relative_entropy[75] == pytest.approx(1.4083272214176723, abs=5e-8)
        assert relative_entropy[76] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[77] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[78] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[79] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[80] == pytest.approx(1.0147718639657484, abs=5e-8)
        assert relative_entropy[81] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[82] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[83] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[84] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[85] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[86] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[87] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[88] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[89] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[90] == pytest.approx(1.136879431433369, abs=5e-8)
        assert relative_entropy[91] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[92] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[93] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[94] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[95] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[96] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[97] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[98] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[99] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[100] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[101] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[102] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[103] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[104] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[105] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[106] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[107] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[108] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[109] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[110] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[111] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[112] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[113] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[114] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[115] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[116] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[117] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[118] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[119] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[120] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[121] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[122] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[123] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[124] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[125] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[126] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[127] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[128] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[129] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[130] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[131] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[132] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[133] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[134] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[135] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[136] == pytest.approx(1.4083272214176723, abs=5e-8)
        assert relative_entropy[137] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[138] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[139] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[140] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[141] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[142] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[143] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[144] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[145] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[146] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[147] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[148] == pytest.approx(1.4083272214176723, abs=5e-8)
        assert relative_entropy[149] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[150] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[151] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[152] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[153] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[154] == pytest.approx(2.0, abs=5e-8)
        assert relative_entropy[155] == pytest.approx(2.0, abs=5e-8)
        # create a new-style Alignment object
        del seq_record
        del align_info
        del dictionary
        del value
        alignment = msa.alignment
        assert len(alignment) == 7
        seq_record = alignment.sequences[0]
        assert seq_record.description == "gi|6273285|gb|AF191659.1|AF191"
        assert alignment[0] == "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCAGA"
        assert seq_record.seq == Seq(
                "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATAATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCAGA"
            )
        seq_record = alignment.sequences[1]
        assert seq_record.description == "gi|6273284|gb|AF191658.1|AF191"
        assert alignment[1] == "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA"
        assert seq_record.seq == "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATAATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA"
        seq_record = alignment.sequences[2]
        assert seq_record.description == "gi|6273287|gb|AF191661.1|AF191"
        assert alignment[2] == "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA"
        assert seq_record.seq == "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATAATATATTTCAAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA"
        seq_record = alignment.sequences[3]
        assert seq_record.description == "gi|6273286|gb|AF191660.1|AF191"
        assert alignment[3] == "TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTATAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA"
        assert seq_record.seq == "TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATAATATATTTATAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA"
        seq_record = alignment.sequences[4]
        assert seq_record.description == "gi|6273290|gb|AF191664.1|AF191"
        assert alignment[4] == "TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA"
        assert seq_record.seq == "TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA"
        seq_record = alignment.sequences[5]
        assert seq_record.description == "gi|6273289|gb|AF191663.1|AF191"
        assert alignment[5] == "TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTATACCAGA"
        assert seq_record.seq == "TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTATACCAGA"
        seq_record = alignment.sequences[6]
        assert seq_record.description == "gi|6273291|gb|AF191665.1|AF191"
        assert alignment[6] == "TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA"
        assert seq_record.seq == "TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA"
        assert alignment.shape == (7, 156)
        substitutions = alignment.substitutions
        assert len(substitutions) == 4
        assert substitutions.shape == (4, 4)
        assert substitutions[("A", "A")] == pytest.approx(1395, abs=5e-8)
        assert substitutions[("A", "C")] == pytest.approx(3, abs=5e-8)
        assert substitutions[("A", "G")] == pytest.approx(13, abs=5e-8)
        assert substitutions[("A", "T")] == pytest.approx(6, abs=5e-8)
        assert substitutions[("C", "A")] == pytest.approx(3, abs=5e-8)
        assert substitutions[("C", "C")] == pytest.approx(271, abs=5e-8)
        assert substitutions[("C", "G")] == pytest.approx(0, abs=5e-8)
        assert substitutions[("C", "T")] == pytest.approx(16, abs=5e-8)
        assert substitutions[("G", "A")] == pytest.approx(5, abs=5e-8)
        assert substitutions[("G", "C")] == pytest.approx(0, abs=5e-8)
        assert substitutions[("G", "G")] == pytest.approx(480, abs=5e-8)
        assert substitutions[("G", "T")] == pytest.approx(0, abs=5e-8)
        assert substitutions[("T", "A")] == pytest.approx(6, abs=5e-8)
        assert substitutions[("T", "C")] == pytest.approx(12, abs=5e-8)
        assert substitutions[("T", "G")] == pytest.approx(0, abs=5e-8)
        assert substitutions[("T", "T")] == pytest.approx(874, abs=5e-8)
        motif = motifs.Motif(alphabet="ACGT", alignment=alignment)
        assert motif.consensus == "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATATATATAATATATTTCAAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA"
        assert motif.degenerate_consensus == "TATACATTAAAGRAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATATATATAATATATTTCAAATTYCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA"
        matrix = motif.counts
        assert str(matrix) == """\
        0      1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16     17     18     19     20     21     22     23     24     25     26     27     28     29     30     31     32     33     34     35     36     37     38     39     40     41     42     43     44     45     46     47     48     49     50     51     52     53     54     55     56     57     58     59     60     61     62     63     64     65     66     67     68     69     70     71     72     73     74     75     76     77     78     79     80     81     82     83     84     85     86     87     88     89     90     91     92     93     94     95     96     97     98     99    100    101    102    103    104    105    106    107    108    109    110    111    112    113    114    115    116    117    118    119    120    121    122    123    124    125    126    127    128    129    130    131    132    133    134    135    136    137    138    139    140    141    142    143    144    145    146    147    148    149    150    151    152    153    154    155
A:   0.00   7.00   0.00   7.00   0.00   7.00   0.00   1.00   7.00   7.00   7.00   0.00   4.00   7.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   7.00   7.00   7.00   0.00   0.00   0.00   7.00   7.00   7.00   0.00   0.00   0.00   0.00   7.00   7.00   7.00   0.00   7.00   7.00   7.00   0.00   7.00   7.00   0.00   7.00   0.00   7.00   0.00   7.00   0.00   4.00   0.00   3.00   0.00   1.00   0.00   1.00   0.00   1.00   7.00   0.00   7.00   0.00   7.00   0.00   0.00   0.00   1.00   6.00   7.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   7.00   0.00   7.00   0.00   0.00   0.00   7.00   7.00   7.00   0.00   7.00   0.00   7.00   7.00   7.00   7.00   7.00   0.00   7.00   0.00   0.00   0.00   7.00   7.00   0.00   7.00   7.00   7.00   0.00   0.00   7.00   0.00   7.00   0.00   0.00   7.00   7.00   0.00   7.00   0.00   0.00   7.00   7.00   7.00   0.00   7.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   1.00   0.00   7.00   0.00   0.00   7.00   0.00   7.00
C:   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   6.00   0.00   0.00   0.00   0.00   0.00   3.00   7.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   2.00   7.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   1.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   7.00   0.00   0.00   0.00
G:   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   3.00   0.00   7.00   7.00   7.00   7.00   7.00   0.00   0.00   7.00   0.00   7.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   7.00   0.00   0.00   0.00   7.00   7.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   6.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00
T:   7.00   0.00   7.00   0.00   0.00   0.00   7.00   6.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   7.00   0.00   7.00   0.00   4.00   0.00   3.00   0.00   1.00   0.00   1.00   0.00   1.00   0.00   0.00   7.00   0.00   7.00   0.00   7.00   7.00   7.00   0.00   1.00   0.00   0.00   7.00   7.00   4.00   0.00   0.00   7.00   7.00   0.00   7.00   0.00   7.00   0.00   5.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   7.00   0.00   7.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   7.00   0.00   0.00   0.00   7.00   0.00   0.00   0.00   7.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   7.00   0.00   6.00   0.00   7.00   7.00   0.00   0.00   7.00   7.00   7.00   0.00   0.00   7.00   0.00   7.00   0.00   0.00   0.00   0.00   0.00   0.00
"""
        assert format(motif, "transfac") == """\
P0      A      C      G      T
01      0      0      0      7      T
02      7      0      0      0      A
03      0      0      0      7      T
04      7      0      0      0      A
05      0      7      0      0      C
06      7      0      0      0      A
07      0      0      0      7      T
08      1      0      0      6      T
09      7      0      0      0      A
10      7      0      0      0      A
11      7      0      0      0      A
12      0      0      7      0      G
13      4      0      3      0      R
14      7      0      0      0      A
15      0      0      7      0      G
16      0      0      7      0      G
17      0      0      7      0      G
18      0      0      7      0      G
19      0      0      7      0      G
20      7      0      0      0      A
21      0      0      0      7      T
22      0      0      7      0      G
23      0      7      0      0      C
24      0      0      7      0      G
25      0      0      7      0      G
26      7      0      0      0      A
27      0      0      0      7      T
28      7      0      0      0      A
29      7      0      0      0      A
30      7      0      0      0      A
31      0      0      0      7      T
32      0      0      7      0      G
33      0      0      7      0      G
34      7      0      0      0      A
35      7      0      0      0      A
36      7      0      0      0      A
37      0      0      7      0      G
38      0      0      7      0      G
39      0      7      0      0      C
40      0      0      7      0      G
41      7      0      0      0      A
42      7      0      0      0      A
43      7      0      0      0      A
44      0      0      7      0      G
45      7      0      0      0      A
46      7      0      0      0      A
47      7      0      0      0      A
48      0      0      7      0      G
49      7      0      0      0      A
50      7      0      0      0      A
51      0      0      0      7      T
52      7      0      0      0      A
53      0      0      0      7      T
54      7      0      0      0      A
55      0      0      0      7      T
56      7      0      0      0      A
57      0      0      0      4      T
58      4      0      0      0      A
59      0      0      0      3      T
60      3      0      0      0      A
61      0      0      0      1      T
62      1      0      0      0      A
63      0      0      0      1      T
64      1      0      0      0      A
65      0      0      0      1      T
66      1      0      0      0      A
67      7      0      0      0      A
68      0      0      0      7      T
69      7      0      0      0      A
70      0      0      0      7      T
71      7      0      0      0      A
72      0      0      0      7      T
73      0      0      0      7      T
74      0      0      0      7      T
75      1      6      0      0      C
76      6      0      0      1      A
77      7      0      0      0      A
78      7      0      0      0      A
79      0      0      0      7      T
80      0      0      0      7      T
81      0      3      0      4      Y
82      0      7      0      0      C
83      0      7      0      0      C
84      0      0      0      7      T
85      0      0      0      7      T
86      7      0      0      0      A
87      0      0      0      7      T
88      7      0      0      0      A
89      0      0      0      7      T
90      7      0      0      0      A
91      0      2      0      5      T
92      0      7      0      0      C
93      0      7      0      0      C
94      7      0      0      0      A
95      7      0      0      0      A
96      7      0      0      0      A
97      0      0      0      7      T
98      7      0      0      0      A
99      0      0      0      7      T
100      7      0      0      0      A
101      7      0      0      0      A
102      7      0      0      0      A
103      7      0      0      0      A
104      7      0      0      0      A
105      0      0      0      7      T
106      7      0      0      0      A
107      0      0      0      7      T
108      0      7      0      0      C
109      0      0      0      7      T
110      7      0      0      0      A
111      7      0      0      0      A
112      0      0      0      7      T
113      7      0      0      0      A
114      7      0      0      0      A
115      7      0      0      0      A
116      0      0      0      7      T
117      0      0      0      7      T
118      7      0      0      0      A
119      0      0      7      0      G
120      7      0      0      0      A
121      0      0      0      7      T
122      0      0      7      0      G
123      7      0      0      0      A
124      7      0      0      0      A
125      0      0      0      7      T
126      7      0      0      0      A
127      0      0      0      7      T
128      0      7      0      0      C
129      7      0      0      0      A
130      7      0      0      0      A
131      7      0      0      0      A
132      0      0      7      0      G
133      7      0      0      0      A
134      7      0      0      0      A
135      0      0      0      7      T
136      0      7      0      0      C
137      0      1      0      6      T
138      7      0      0      0      A
139      0      0      0      7      T
140      0      0      0      7      T
141      0      0      7      0      G
142      7      0      0      0      A
143      0      0      0      7      T
144      0      0      0      7      T
145      0      0      0      7      T
146      7      0      0      0      A
147      0      0      7      0      G
148      0      0      0      7      T
149      1      0      6      0      G
150      0      0      0      7      T
151      7      0      0      0      A
152      0      7      0      0      C
153      0      7      0      0      C
154      7      0      0      0      A
155      0      0      7      0      G
156      7      0      0      0      A
XX
//
"""
        assert sum(motif[5:50].relative_entropy) == pytest.approx(88.42309908538343, abs=5e-8)
        relative_entropy = motif.relative_entropy
        assert sum(relative_entropy[5:50]) == pytest.approx(88.42309908538343, abs=5e-8)
        assert sum(relative_entropy) == pytest.approx(306.20805926645323, abs=5e-8)
        assert alignment[:, 1] == "AAAAAAA"
        assert motif.relative_entropy[1] == pytest.approx(2.0, abs=5e-8)
        assert alignment[:, 7] == "TTTATTT"
        assert relative_entropy[7] == pytest.approx(1.4083272214176723, abs=5e-8)

    def test_read_fasta(self):
        path = os.path.join(os.curdir, "Quality", "example.fasta")
        msa = AlignIO.read(path, "fasta")
        assert len(msa) == 3
        seq_record = msa[0]
        assert seq_record.description == "EAS54_6_R1_2_1_413_324"
        assert seq_record.seq == "CCCTTCTTGTCTTCAGCGTTTCTCC"
        seq_record = msa[1]
        assert seq_record.description == "EAS54_6_R1_2_1_540_792"
        assert seq_record.seq == "TTGGCAGGCCAAGGCCGATGGATCA"
        seq_record = msa[2]
        assert seq_record.description == "EAS54_6_R1_2_1_443_348"
        assert seq_record.seq == "GTTGCTTCTGGCGTGGGTGGGGGGG"
        assert msa.get_alignment_length() == 25
        assert str(msa) == """\
Alignment with 3 rows and 25 columns
CCCTTCTTGTCTTCAGCGTTTCTCC EAS54_6_R1_2_1_413_324
TTGGCAGGCCAAGGCCGATGGATCA EAS54_6_R1_2_1_540_792
GTTGCTTCTGGCGTGGGTGGGGGGG EAS54_6_R1_2_1_443_348"""
        alignment = msa.alignment
        assert len(alignment) == 3
        seq_record = alignment.sequences[0]
        assert seq_record.description == "EAS54_6_R1_2_1_413_324"
        assert seq_record.seq == "CCCTTCTTGTCTTCAGCGTTTCTCC"
        seq_record = alignment.sequences[1]
        assert seq_record.description == "EAS54_6_R1_2_1_540_792"
        assert seq_record.seq == "TTGGCAGGCCAAGGCCGATGGATCA"
        seq_record = alignment.sequences[2]
        assert seq_record.description == "EAS54_6_R1_2_1_443_348"
        assert seq_record.seq == "GTTGCTTCTGGCGTGGGTGGGGGGG"
        assert alignment.length == 25
        motif = motifs.Motif(alphabet="ACGT", alignment=alignment)
        assert motif.consensus == "CTCGCATCCCAAGCAGGATGGATCA"
        assert motif.degenerate_consensus == "BYBKYHKBBBVHKBVSSDKKKVKSV"
        assert str(alignment) == """\
EAS54_6_R         0 CCCTTCTTGTCTTCAGCGTTTCTCC 25
EAS54_6_R         0 TTGGCAGGCCAAGGCCGATGGATCA 25
EAS54_6_R         0 GTTGCTTCTGGCGTGGGTGGGGGGG 25
"""
        assert motif.counts.calculate_consensus(identity=0.6) == "NTNGCNTNNNNNGNNGGNTGGNTCN"


class TestFromAlignmentsWithSameReference(unittest.TestCase):

    def test_pwas_built_with_strings(self):
        """Test that from_alignments_with_same_reference creates the correct Alignment when PWA inputs were built with strings"""
        aligner = PairwiseAligner()

        # Input sequences
        reference_str = "ACGT"
        seq1_str = "ACGGT"
        seq2_str = "AT"

        # Generate pairwise alignments
        pwa1 = next(aligner.align(reference_str, seq1_str))
        pwa2 = next(aligner.align(reference_str, seq2_str))

        # Use the method being tested
        msa = Alignment.from_alignments_with_same_reference([pwa1, pwa2])

        # Check that the output is of correct type
        assert isinstance(msa, Alignment)

        # Check the number of sequences in the MSA
        assert len(msa) == 3

        # Validate that from_alignments_with_same_reference gives the right msa
        assert str(msa[0]) == "ACG-T"
        assert str(msa[1]) == "ACGGT"
        assert str(msa[2]) == "A---T"

    def test_pwas_built_with_seqs(self):
        """Test that from_alignments_with_same_reference works with PWAs built with Seq objects."""
        aligner = PairwiseAligner()

        # Input sequences
        reference_seq = Seq("ACGT")
        seq1_seq = Seq("ACGGT")
        seq2_seq = Seq("AT")

        # Generate pairwise alignments
        pwa1 = next(aligner.align(reference_seq, seq1_seq))
        pwa2 = next(aligner.align(reference_seq, seq2_seq))

        # Use the method being tested
        msa = Alignment.from_alignments_with_same_reference([pwa1, pwa2])

        # Check that the output is of correct type
        assert isinstance(msa, Alignment)

        # Validate that from_alignments_with_same_reference gives the right msa
        assert str(msa[0]) == "ACG-T"
        assert str(msa[1]) == "ACGGT"
        assert str(msa[2]) == "A---T"

        # Works with undefined sequences
        pwa1_undefined = Alignment([Seq(None, 4), Seq(None, 5)], pwa1.coordinates)
        pwa2_undefined = Alignment([Seq(None, 4), Seq(None, 2)], pwa2.coordinates)
        msa_undefined = Alignment.from_alignments_with_same_reference(
            [pwa1_undefined, pwa2_undefined]
        )
        assert (msa_undefined.coordinates == msa.coordinates).all()

    def test_metadata_preservation(self):
        """Test that sequence metadata (IDs and descriptions) are preserved in the Alignment.

        Sequence metadata may exist if pairwise alignments are built with SeqRecords.
        This also tests that from_alignments_with_same_reference works with pairwise alignments built with SeqRecord objects.
        """
        aligner = PairwiseAligner()

        # Input sequences
        reference_seqr = SeqRecord(Seq("ACGT"), id="reference", description="desc 1")
        seq1_seqr = SeqRecord(Seq("ACGGT"), id="seq1", description="desc 2")
        seq2_seqr = SeqRecord(Seq("AT"), id="seq2", description="desc 3")

        # Generate pairwise alignments
        pwa1 = next(aligner.align(reference_seqr, seq1_seqr))
        pwa2 = next(aligner.align(reference_seqr, seq2_seqr))

        # Use the method being tested
        msa = Alignment.from_alignments_with_same_reference([pwa1, pwa2])

        # Validate that from_alignments_with_same_reference gives the right msa
        assert str(msa[0]) == "ACG-T"
        assert str(msa[1]) == "ACGGT"
        assert str(msa[2]) == "A---T"

        # Check that id and description are retained
        assert msa.sequences[0].id == "reference"
        assert msa.sequences[1].id == "seq1"
        assert msa.sequences[2].id == "seq2"
        assert msa.sequences[0].description == "desc 1"
        assert msa.sequences[1].description == "desc 2"

    def test_mismatched_references(self):
        """Test that mismatched reference sequences raise a ValueError."""
        aligner = PairwiseAligner()
        pwa1 = next(aligner.align("ACGAAAT", "ACGGT"))
        pwa2 = next(aligner.align("ACCT", "AT"))

        with pytest.raises(ValueError) as context:
            Alignment.from_alignments_with_same_reference([pwa1, pwa2])
        assert "All reference sequences must" in str(context.value)

        # Works with undefined sequences
        pwa1_undefined = Alignment([Seq(None, 7), Seq(None, 5)], pwa1.coordinates)
        pwa2_undefined = Alignment([Seq(None, 4), Seq(None, 2)], pwa2.coordinates)
        with pytest.raises(ValueError) as context:
            Alignment.from_alignments_with_same_reference(
                [pwa1_undefined, pwa2_undefined]
            )
        assert "All reference sequences must" in str(context.value)

        # Works with mix
        pwa1_undefined = Alignment([Seq(None, 7), Seq(None, 5)], pwa1.coordinates)
        with pytest.raises(ValueError) as context:
            Alignment.from_alignments_with_same_reference([pwa1_undefined, pwa2])
        assert "All reference sequences must" in str(context.value)

    def test_empty_input(self):
        """Test that empty input raises a ValueError."""
        with pytest.raises(ValueError):
            Alignment.from_alignments_with_same_reference([])

    def test_input_with_three_sequences(self):
        """Test that the method works with input Alignments of three sequences"""
        # Input sequences
        reference_str = "ACGT"
        seq1_str = "ACGGT"
        seq2_str = "AT"
        seq3_str = "AGT"
        seq4_str = "ACT"

        # Coordinates for each alignment
        coords1 = np.array([[0, 1, 2, 3, 3, 4], [0, 1, 2, 3, 4, 5], [0, 1, 1, 1, 1, 2]])

        coords2 = np.array([[0, 1, 2, 3, 4], [0, 1, 1, 2, 3], [0, 1, 2, 2, 3]])

        # Generate input alignments
        alignment1 = Alignment([reference_str, seq1_str, seq2_str], coords1)
        alignment2 = Alignment([reference_str, seq3_str, seq4_str], coords2)

        # Use the method being tested
        msa = Alignment.from_alignments_with_same_reference([alignment1, alignment2])

        # Check that the output is of correct type
        assert isinstance(msa, Alignment)

        # Check the number of sequences in the MSA
        assert len(msa) == 5

        # Validate that from_alignments_with_same_reference gives the right msa
        assert str(msa[0]) == "ACG-T"
        assert str(msa[1]) == "ACGGT"
        assert str(msa[2]) == "A---T"
        assert str(msa[3]) == "A-G-T"
        assert str(msa[4]) == "AC--T"

    def test_mixed_input(self):
        """Test that the method works with input Alignments of mixed number of sequences"""
        aligner = PairwiseAligner()

        # Input sequences
        reference_str = "ACGT"
        seq1_str = "ACT"
        seq2_str = "ACGGT"
        seq3_str = "AT"

        # Coordinates for an alignment of three sequences
        coords = np.array([[0, 1, 2, 3, 3, 4], [0, 1, 2, 3, 4, 5], [0, 1, 1, 1, 1, 2]])

        # Generate input alignments
        pwa = next(aligner.align(reference_str, seq1_str))
        not_pwa = Alignment([reference_str, seq2_str, seq3_str], coords)

        # Use the method being tested
        msa = Alignment.from_alignments_with_same_reference([pwa, not_pwa])

        # Check that the output is of correct type
        assert isinstance(msa, Alignment)

        # Check the number of sequences in the MSA
        assert len(msa) == 4

        # Validate that from_alignments_with_same_reference gives the right msa
        assert str(msa[0]) == "ACG-T"
        assert str(msa[1]) == "AC--T"
        assert str(msa[2]) == "ACGGT"
        assert str(msa[3]) == "A---T"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
