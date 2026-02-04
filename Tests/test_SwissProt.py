# Copyright 2009 by Michiel de Hoon.  All rights reserved.
# Revisions copyright 2010 by Peter Cock. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Test for the SwissProt parser on SwissProt files."""

import os
import unittest
import pytest

from Bio import SeqIO
from Bio import SwissProt
from Bio.SeqRecord import SeqRecord


class TestSwissProt(unittest.TestCase):
    def test_Q13454(self):
        """Parsing SwissProt file Q13454.txt."""
        filename = "Q13454.txt"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        assert isinstance(seq_record, SeqRecord)

        assert seq_record.id == "Q13454"
        assert seq_record.name == "TUSC3_HUMAN"
        assert seq_record.description == "RecName: Full=Tumor suppressor candidate 3; AltName: Full=Dolichyl-diphosphooligosaccharide--protein glycosyltransferase subunit TUSC3; Short=Oligosaccharyl transferase subunit TUSC3; AltName: Full=Magnesium uptake/transporter TUSC3; AltName: Full=Protein N33; Flags: Precursor;"
        assert repr(seq_record.seq) == "Seq('MGARGAPSRRRQAGRRLRYLPTGSFPFLLLLLLLCIQLGGGQKKKENLLAEKVE...DFE')"

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        assert record.entry_name == "TUSC3_HUMAN"
        assert record.accessions == ["Q13454", "A8MSM0", "D3DSP2", "Q14911", "Q14912", "Q96FW0"]
        assert record.gene_name == [{"Name": "TUSC3", "Synonyms": ["N33"]}]
        assert record.organism_classification == [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Euteleostomi",
                "Mammalia",
                "Eutheria",
                "Euarchontoglires",
                "Primates",
                "Haplorrhini",
                "Catarrhini",
                "Hominidae",
                "Homo",
            ]
        assert record.seqinfo == (348, 39676, "16D97CB1E00C5190")

        assert len(record.features) == 32
        feature = record.features[0]
        assert feature.type == "SIGNAL"
        assert feature.location.start == 0
        assert feature.location.end == 41
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        assert feature.id is None
        feature = record.features[1]
        assert feature.type == "CHAIN"
        assert feature.location.start == 41
        assert feature.location.end == 348
        assert feature.qualifiers["note"] == "Tumor suppressor candidate 3"
        assert feature.id == "PRO_0000215300"
        feature = record.features[2]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 41
        assert feature.location.end == 196
        assert feature.qualifiers["note"] == "Lumenal"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        assert feature.id is None
        feature = record.features[3]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 196
        assert feature.location.end == 217
        assert feature.qualifiers["note"] == "Helical"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        assert feature.id is None
        feature = record.features[4]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 217
        assert feature.location.end == 221
        assert feature.qualifiers["note"] == "Cytoplasmic"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        assert feature.id is None
        feature = record.features[5]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 221
        assert feature.location.end == 242
        assert feature.qualifiers["note"] == "Helical"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        assert feature.id is None
        feature = record.features[6]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 242
        assert feature.location.end == 276
        assert feature.qualifiers["note"] == "Lumenal"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        assert feature.id is None
        feature = record.features[7]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 276
        assert feature.location.end == 297
        assert feature.qualifiers["note"] == "Helical"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        assert feature.id is None
        feature = record.features[8]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 297
        assert feature.location.end == 312
        assert feature.qualifiers["note"] == "Cytoplasmic"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        assert feature.id is None
        feature = record.features[9]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 312
        assert feature.location.end == 333
        assert feature.qualifiers["note"] == "Helical"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        assert feature.id is None
        feature = record.features[10]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 333
        assert feature.location.end == 348
        assert feature.qualifiers["note"] == "Lumenal"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        assert feature.id is None
        feature = record.features[11]
        assert feature.type == "DOMAIN"
        assert feature.location.start == 58
        assert feature.location.end == 187
        assert feature.qualifiers["note"] == "Thioredoxin"
        assert feature.id is None
        feature = record.features[12]
        assert feature.type == "DISULFID"
        assert feature.location.start == 98
        assert feature.location.end == 102
        assert feature.qualifiers["note"] == "Redox-active"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:24685145,ECO:0007744|PDB:4M8G, ECO:0007744|PDB:4M90"
        assert feature.id is None
        feature = record.features[13]
        assert feature.type == "VAR_SEQ"
        assert feature.location.start == 343
        assert feature.location.end == 348
        assert feature.qualifiers["note"] == "DLDFE -> FLIK (in isoform 2)"
        assert feature.qualifiers["evidence"] == "ECO:0000303|PubMed:15489334,ECO:0000303|PubMed:8661104"
        assert feature.id == "VSP_003776"
        feature = record.features[14]
        assert feature.type == "VARIANT"
        assert feature.location.start == 64
        assert feature.location.end == 65
        assert feature.qualifiers["note"] == "I -> V (in dbSNP:rs11545035)"
        assert feature.id == "VAR_045836"
        feature = record.features[15]
        assert feature.type == "VARIANT"
        assert feature.location.start == 246
        assert feature.location.end == 247
        assert feature.qualifiers["note"] == "M -> V"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:23033978"
        assert feature.id == "VAR_069369"
        feature = record.features[16]
        assert feature.type == "MUTAGEN"
        assert feature.location.start == 98
        assert feature.location.end == 99
        assert feature.qualifiers["note"] == "C->S: Reduces N-glycosylation of cysteine-proximal acceptor sites; when associated with S-102."
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:25135935"
        assert feature.id is None
        feature = record.features[17]
        assert feature.type == "MUTAGEN"
        assert feature.location.start == 101
        assert feature.location.end == 102
        assert feature.qualifiers["note"] == "C->S: Reduces N-glycosylation of cysteine-proximal acceptor sites; when associated with S-99."
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:25135935"
        assert feature.id is None
        feature = record.features[18]
        assert feature.type == "HELIX"
        assert feature.location.start == 43
        assert feature.location.end == 62
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:4M91"
        assert feature.id is None
        feature = record.features[19]
        assert feature.type == "STRAND"
        assert feature.location.start == 63
        assert feature.location.end == 67
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:4M91"
        assert feature.id is None
        feature = record.features[20]
        assert feature.type == "HELIX"
        assert feature.location.start == 69
        assert feature.location.end == 76
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:4M91"
        assert feature.id is None
        feature = record.features[21]
        assert feature.type == "STRAND"
        assert feature.location.start == 83
        assert feature.location.end == 91
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:4M91"
        assert feature.id is None
        feature = record.features[22]
        assert feature.type == "HELIX"
        assert feature.location.start == 94
        assert feature.location.end == 97
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:4M91"
        assert feature.id is None
        feature = record.features[23]
        assert feature.type == "HELIX"
        assert feature.location.start == 99
        assert feature.location.end == 118
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:4M91"
        assert feature.id is None
        feature = record.features[24]
        assert feature.type == "STRAND"
        assert feature.location.start == 125
        assert feature.location.end == 132
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:4M91"
        assert feature.id is None
        feature = record.features[25]
        assert feature.type == "TURN"
        assert feature.location.start == 132
        assert feature.location.end == 135
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:4M91"
        assert feature.id is None
        feature = record.features[26]
        assert feature.type == "HELIX"
        assert feature.location.start == 136
        assert feature.location.end == 142
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:4M91"
        assert feature.id is None
        feature = record.features[27]
        assert feature.type == "STRAND"
        assert feature.location.start == 149
        assert feature.location.end == 154
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:4M91"
        assert feature.id is None
        feature = record.features[28]
        assert feature.type == "STRAND"
        assert feature.location.start == 155
        assert feature.location.end == 158
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:4M91"
        assert feature.id is None
        feature = record.features[29]
        assert feature.type == "HELIX"
        assert feature.location.start == 161
        assert feature.location.end == 164
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:4M91"
        assert feature.id is None
        feature = record.features[30]
        assert feature.type == "HELIX"
        assert feature.location.start == 167
        assert feature.location.end == 171
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:4M91"
        assert feature.id is None
        feature = record.features[31]
        assert feature.type == "HELIX"
        assert feature.location.start == 174
        assert feature.location.end == 186
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:4M91"
        assert feature.id is None
        assert len(record.references) == 12
        assert record.references[0].authors == "Macgrogan D., Levy A., Bova G.S., Isaacs W.B., Bookstein R."
        assert record.references[0].title == "Structure and methylation-associated silencing of a gene within a homozygously deleted region of human chromosome band 8p22."
        assert record.references[0].positions == ["NUCLEOTIDE SEQUENCE [GENOMIC DNA / MRNA] (ISOFORMS 1 AND 2)."]
        assert len(record.references[0].references) == 2
        assert record.references[0].references[0] == ("PubMed", "8661104")
        assert record.references[0].references[1] == ("DOI", "10.1006/geno.1996.0322")
        assert record.references[1].authors == "Kalnine N., Chen X., Rolfs A., Halleck A., Hines L., Eisenstein S., Koundinya M., Raphael J., Moreira D., Kelley T., LaBaer J., Lin Y., Phelan M., Farmer A."
        assert record.references[1].title == "Cloning of human full-length CDSs in BD Creator(TM) system donor vector."
        assert record.references[1].positions == ["NUCLEOTIDE SEQUENCE [LARGE SCALE MRNA] (ISOFORM 1)."]
        assert len(record.references[1].references) == 0
        assert record.references[2].authors == "Nusbaum C., Mikkelsen T.S., Zody M.C., Asakawa S., Taudien S., Garber M., Kodira C.D., Schueler M.G., Shimizu A., Whittaker C.A., Chang J.L., Cuomo C.A., Dewar K., FitzGerald M.G., Yang X., Allen N.R., Anderson S., Asakawa T., Blechschmidt K., Bloom T., Borowsky M.L., Butler J., Cook A., Corum B., DeArellano K., DeCaprio D., Dooley K.T., Dorris L. III, Engels R., Gloeckner G., Hafez N., Hagopian D.S., Hall J.L., Ishikawa S.K., Jaffe D.B., Kamat A., Kudoh J., Lehmann R., Lokitsang T., Macdonald P., Major J.E., Matthews C.D., Mauceli E., Menzel U., Mihalev A.H., Minoshima S., Murayama Y., Naylor J.W., Nicol R., Nguyen C., O'Leary S.B., O'Neill K., Parker S.C.J., Polley A., Raymond C.K., Reichwald K., Rodriguez J., Sasaki T., Schilhabel M., Siddiqui R., Smith C.L., Sneddon T.P., Talamas J.A., Tenzin P., Topham K., Venkataraman V., Wen G., Yamazaki S., Young S.K., Zeng Q., Zimmer A.R., Rosenthal A., Birren B.W., Platzer M., Shimizu N., Lander E.S."
        assert record.references[2].title == "DNA sequence and analysis of human chromosome 8."
        assert record.references[2].positions == ["NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA]."]
        assert len(record.references[2].references) == 2
        assert record.references[2].references[0] == ("PubMed", "16421571")
        assert record.references[2].references[1] == ("DOI", "10.1038/nature04406")
        assert record.references[3].authors == "Mural R.J., Istrail S., Sutton G.G., Florea L., Halpern A.L., Mobarry C.M., Lippert R., Walenz B., Shatkay H., Dew I., Miller J.R., Flanigan M.J., Edwards N.J., Bolanos R., Fasulo D., Halldorsson B.V., Hannenhalli S., Turner R., Yooseph S., Lu F., Nusskern D.R., Shue B.C., Zheng X.H., Zhong F., Delcher A.L., Huson D.H., Kravitz S.A., Mouchard L., Reinert K., Remington K.A., Clark A.G., Waterman M.S., Eichler E.E., Adams M.D., Hunkapiller M.W., Myers E.W., Venter J.C."
        assert record.references[3].title == ""
        assert record.references[3].positions == ["NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA]."]
        assert len(record.references[3].references) == 0
        assert record.references[4].authors == "The MGC Project Team"
        assert record.references[4].title == "The status, quality, and expansion of the NIH full-length cDNA project: the Mammalian Gene Collection (MGC)."
        assert record.references[4].positions == ["NUCLEOTIDE SEQUENCE [LARGE SCALE MRNA] (ISOFORM 2)."]
        assert len(record.references[4].references) == 2
        assert record.references[4].references[0] == ("PubMed", "15489334")
        assert record.references[4].references[1] == ("DOI", "10.1101/gr.2596504")
        assert record.references[5].authors == "Kelleher D.J., Karaoglu D., Mandon E.C., Gilmore R."
        assert record.references[5].title == "Oligosaccharyltransferase isoforms that contain different catalytic STT3 subunits have distinct enzymatic properties."
        assert record.references[5].positions == [
                "IDENTIFICATION IN THE OLIGOSACCHARYLTRANSFERASE (OST) COMPLEX, AND TISSUE",
                "SPECIFICITY.",
            ]
        assert len(record.references[5].references) == 2
        assert record.references[5].references[0] == ("PubMed", "12887896")
        assert record.references[5].references[1] == ("DOI", "10.1016/s1097-2765(03)00243-0")
        assert record.references[6].authors == "Molinari F., Foulquier F., Tarpey P.S., Morelle W., Boissel S., Teague J., Edkins S., Futreal P.A., Stratton M.R., Turner G., Matthijs G., Gecz J., Munnich A., Colleaux L."
        assert record.references[6].title == "Oligosaccharyltransferase-subunit mutations in nonsyndromic mental retardation."
        assert record.references[6].positions == ["INVOLVEMENT IN MRT7."]
        assert len(record.references[6].references) == 2
        assert record.references[6].references[0] == ("PubMed", "18455129")
        assert record.references[6].references[1] == ("DOI", "10.1016/j.ajhg.2008.03.021")
        assert record.references[7].authors == "Garshasbi M., Hadavi V., Habibi H., Kahrizi K., Kariminejad R., Behjati F., Tzschach A., Najmabadi H., Ropers H.H., Kuss A.W."
        assert record.references[7].title == "A defect in the TUSC3 gene is associated with autosomal recessive mental retardation."
        assert record.references[7].positions == ["INVOLVEMENT IN MRT7."]
        assert len(record.references[7].references) == 2
        assert record.references[7].references[0] == ("PubMed", "18452889")
        assert record.references[7].references[1] == ("DOI", "10.1016/j.ajhg.2008.03.018")
        assert record.references[8].authors == "Zhou H., Clapham D.E."
        assert record.references[8].title == "Mammalian MagT1 and TUSC3 are required for cellular magnesium uptake and vertebrate embryonic development."
        assert record.references[8].positions == ["FUNCTION IN MAGNESIUM UPTAKE."]
        assert len(record.references[8].references) == 2
        assert record.references[8].references[0] == ("PubMed", "19717468")
        assert record.references[8].references[1] == ("DOI", "10.1073/pnas.0908332106")
        assert record.references[9].authors == "Cherepanova N.A., Shrimal S., Gilmore R."
        assert record.references[9].title == "Oxidoreductase activity is necessary for N-glycosylation of cysteine-proximal acceptor sites in glycoproteins."
        assert record.references[9].positions == ["FUNCTION, AND MUTAGENESIS OF CYS-99 AND CYS-102."]
        assert len(record.references[9].references) == 2
        assert record.references[9].references[0] == ("PubMed", "25135935")
        assert record.references[9].references[1] == ("DOI", "10.1083/jcb.201404083")
        assert record.references[10].authors == "Mohorko E., Owen R.L., Malojcic G., Brozzo M.S., Aebi M., Glockshuber R."
        assert record.references[10].title == "Structural basis of substrate specificity of human oligosaccharyl transferase subunit N33/Tusc3 and its role in regulating protein N-glycosylation."
        assert record.references[10].positions == [
                "X-RAY CRYSTALLOGRAPHY (1.10 ANGSTROMS) OF 44-194, DISULFIDE BOND, PROPOSED",
                "FUNCTION, AND SUBUNIT.",
            ]
        assert len(record.references[10].references) == 2
        assert record.references[10].references[0] == ("PubMed", "24685145")
        assert record.references[10].references[1] == ("DOI", "10.1016/j.str.2014.02.013")
        assert record.references[11].authors == "de Ligt J., Willemsen M.H., van Bon B.W., Kleefstra T., Yntema H.G., Kroes T., Vulto-van Silfhout A.T., Koolen D.A., de Vries P., Gilissen C., del Rosario M., Hoischen A., Scheffer H., de Vries B.B., Brunner H.G., Veltman J.A., Vissers L.E."
        assert record.references[11].title == "Diagnostic exome sequencing in persons with severe intellectual disability."
        assert record.references[11].positions == ["VARIANT VAL-247."]
        assert len(record.references[11].references) == 2
        assert record.references[11].references[0] == ("PubMed", "23033978")
        assert record.references[11].references[1] == ("DOI", "10.1056/nejmoa1206524")

        # Check that the two parsers agree on the essentials
        assert seq_record.seq == record.sequence
        assert seq_record.description == record.description
        assert seq_record.name == record.entry_name
        assert seq_record.id in record.accessions

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        with open(datafile) as test_handle:
            records = list(SeqIO.parse(test_handle, "swiss"))

        assert len(records) == 1
        assert isinstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        assert records[0].seq == seq_record.seq
        assert records[0].description == seq_record.description
        assert records[0].name == seq_record.name
        assert records[0].id == seq_record.id

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        assert len(records) == 1
        assert isinstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        assert records[0].sequence == record.sequence
        assert records[0].description == record.description
        assert records[0].entry_name == record.entry_name
        assert records[0].accessions == record.accessions

    def test_P60904(self):
        """Parsing SwissProt file P60904.txt."""
        filename = "P60904.txt"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        assert isinstance(seq_record, SeqRecord)

        assert seq_record.id == "P60904"
        assert seq_record.name == "DNJC5_MOUSE"
        assert seq_record.description == "RecName: Full=DnaJ homolog subfamily C member 5 {ECO:0000305}; AltName: Full=Cysteine string protein {ECO:0000250|UniProtKB:Q9H3Z4}; Short=CSP {ECO:0000250|UniProtKB:Q9H3Z4}; AltName: Full=Cysteine-string protein isoform alpha {ECO:0000303|PubMed:20847230};"
        assert repr(seq_record.seq) == "Seq('MADQRQRSLSTSGESLYHVLGLDKNATSDDIKKSYRKLALKYHPDKNPDNPEAA...GFN')"

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        assert record.entry_name == "DNJC5_MOUSE"
        assert record.accessions == ["P60904", "P54101"]
        assert record.organism_classification == [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Euteleostomi",
                "Mammalia",
                "Eutheria",
                "Euarchontoglires",
                "Glires",
                "Rodentia",
                "Myomorpha",
                "Muroidea",
                "Muridae",
                "Murinae",
                "Mus",
                "Mus",
            ]
        assert record.seqinfo == (198, 22101, "52F98261FBAD978F")

        assert len(record.features) == 23
        feature = record.features[0]
        assert feature.type == "CHAIN"
        assert feature.location.start == 0
        assert feature.location.end == 198
        assert feature.qualifiers["note"] == "DnaJ homolog subfamily C member 5"
        feature = record.features[1]
        assert feature.type == "DOMAIN"
        assert feature.location.start == 12
        assert feature.location.end == 82
        assert feature.qualifiers["note"] == "J"
        assert feature.qualifiers["evidence"] == "ECO:0000255|PROSITE-ProRule:PRU00286"
        feature = record.features[2]
        assert feature.type == "MOD_RES"
        assert feature.location.start == 7
        assert feature.location.end == 8
        assert feature.qualifiers["note"] == "Phosphoserine"
        assert feature.qualifiers["evidence"] == "ECO:0000250|UniProtKB:Q9H3Z4"
        feature = record.features[3]
        assert feature.type == "MOD_RES"
        assert feature.location.start == 9
        assert feature.location.end == 10
        assert feature.qualifiers["note"] == "Phosphoserine"
        assert feature.qualifiers["evidence"] == "ECO:0007744|PubMed:19131326"
        feature = record.features[4]
        assert feature.type == "MOD_RES"
        assert feature.location.start == 11
        assert feature.location.end == 12
        assert feature.qualifiers["note"] == "Phosphoserine"
        assert feature.qualifiers["evidence"] == "ECO:0000250|UniProtKB:Q9H3Z4"
        feature = record.features[5]
        assert feature.type == "MOD_RES"
        assert feature.location.start == 14
        assert feature.location.end == 15
        assert feature.qualifiers["note"] == "Phosphoserine"
        assert feature.qualifiers["evidence"] == "ECO:0007744|PubMed:21183079"
        feature = record.features[6]
        assert feature.type == "MOD_RES"
        assert feature.location.start == 16
        assert feature.location.end == 17
        assert feature.qualifiers["note"] == "Phosphotyrosine"
        assert feature.qualifiers["evidence"] == "ECO:0000250|UniProtKB:Q9H3Z4"
        feature = record.features[7]
        assert feature.type == "MOD_RES"
        assert feature.location.start == 55
        assert feature.location.end == 56
        assert feature.qualifiers["note"] == "N6-acetyllysine"
        assert feature.qualifiers["evidence"] == "ECO:0000250|UniProtKB:Q9H3Z4"
        feature = record.features[8]
        assert feature.type == "MOD_RES"
        assert feature.location.start == 150
        assert feature.location.end == 151
        assert feature.qualifiers["note"] == "Phosphoserine"
        assert feature.qualifiers["evidence"] == "ECO:0007744|PubMed:21183079"
        feature = record.features[9]
        assert feature.type == "MUTAGEN"
        assert feature.location.start == 9
        assert feature.location.end == 10
        assert feature.qualifiers["note"] == "S->D: Reduced interaction with SYT9, but no effect on the interaction with HSC70."
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:20847230"
        feature = record.features[10]
        assert feature.type == "MUTAGEN"
        assert feature.location.start == 92
        assert feature.location.end == 93
        assert feature.qualifiers["note"] == "E->V: Reduced interaction with SYT9."
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:20847230"
        feature = record.features[11]
        assert feature.type == "MUTAGEN"
        assert feature.location.start == 112
        assert feature.location.end == 113
        assert feature.qualifiers["note"] == "C->V: No effect on palmitoylation. No change in subcellular location; when associated with G-118 and F-121."
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:17034881"
        feature = record.features[12]
        assert feature.type == "MUTAGEN"
        assert feature.location.start == 117
        assert feature.location.end == 118
        assert feature.qualifiers["note"] == "C->G: No effect on palmitoylation. No change in subcellular location; when associated with V-113 and F-121."
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:17034881"
        feature = record.features[13]
        assert feature.type == "MUTAGEN"
        assert feature.location.start == 120
        assert feature.location.end == 121
        assert feature.qualifiers["note"] == "C->F: No effect on palmitoylation. No change in subcellular location; when associated with V-113 and G-118."
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:17034881"
        feature = record.features[14]
        assert feature.type == "MUTAGEN"
        assert feature.location.start == 128
        assert feature.location.end == 129
        assert feature.qualifiers["note"] == "F->C: No effect on palmitoylation. No change in subcellular location; when associated with H-135."
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:17034881"
        feature = record.features[15]
        assert feature.type == "MUTAGEN"
        assert feature.location.start == 134
        assert feature.location.end == 135
        assert feature.qualifiers["note"] == "K->H: No effect on palmitoylation. No change in subcellular location; when associated with C-129."
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:17034881"
        feature = record.features[16]
        assert feature.type == "HELIX"
        assert feature.location.start == 15
        assert feature.location.end == 20
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:2CTW"
        feature = record.features[17]
        assert feature.type == "HELIX"
        assert feature.location.start == 27
        assert feature.location.end == 41
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:2CTW"
        feature = record.features[18]
        assert feature.type == "TURN"
        assert feature.location.start == 43
        assert feature.location.end == 46
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:2CTW"
        feature = record.features[19]
        assert feature.type == "HELIX"
        assert feature.location.start == 50
        assert feature.location.end == 67
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:2CTW"
        feature = record.features[20]
        assert feature.type == "HELIX"
        assert feature.location.start == 69
        assert feature.location.end == 78
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:2CTW"
        feature = record.features[21]
        assert feature.type == "HELIX"
        assert feature.location.start == 80
        assert feature.location.end == 89
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:2CTW"
        feature = record.features[22]
        assert feature.type == "HELIX"
        assert feature.location.start == 93
        assert feature.location.end == 100
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:2CTW"
        assert feature.id is None
        assert len(record.references) == 16
        reference = record.references[0]
        assert reference.authors == "Qin N., Lin T., Birnbaumer L."
        assert reference.title == ""
        assert len(reference.references) == 0
        reference = record.references[1]
        assert reference.authors == "Carninci P., Kasukawa T., Katayama S., Gough J., Frith M.C., Maeda N., Oyama R., Ravasi T., Lenhard B., Wells C., Kodzius R., Shimokawa K., Bajic V.B., Brenner S.E., Batalov S., Forrest A.R., Zavolan M., Davis M.J., Wilming L.G., Aidinis V., Allen J.E., Ambesi-Impiombato A., Apweiler R., Aturaliya R.N., Bailey T.L., Bansal M., Baxter L., Beisel K.W., Bersano T., Bono H., Chalk A.M., Chiu K.P., Choudhary V., Christoffels A., Clutterbuck D.R., Crowe M.L., Dalla E., Dalrymple B.P., de Bono B., Della Gatta G., di Bernardo D., Down T., Engstrom P., Fagiolini M., Faulkner G., Fletcher C.F., Fukushima T., Furuno M., Futaki S., Gariboldi M., Georgii-Hemming P., Gingeras T.R., Gojobori T., Green R.E., Gustincich S., Harbers M., Hayashi Y., Hensch T.K., Hirokawa N., Hill D., Huminiecki L., Iacono M., Ikeo K., Iwama A., Ishikawa T., Jakt M., Kanapin A., Katoh M., Kawasawa Y., Kelso J., Kitamura H., Kitano H., Kollias G., Krishnan S.P., Kruger A., Kummerfeld S.K., Kurochkin I.V., Lareau L.F., Lazarevic D., Lipovich L., Liu J., Liuni S., McWilliam S., Madan Babu M., Madera M., Marchionni L., Matsuda H., Matsuzawa S., Miki H., Mignone F., Miyake S., Morris K., Mottagui-Tabar S., Mulder N., Nakano N., Nakauchi H., Ng P., Nilsson R., Nishiguchi S., Nishikawa S., Nori F., Ohara O., Okazaki Y., Orlando V., Pang K.C., Pavan W.J., Pavesi G., Pesole G., Petrovsky N., Piazza S., Reed J., Reid J.F., Ring B.Z., Ringwald M., Rost B., Ruan Y., Salzberg S.L., Sandelin A., Schneider C., Schoenbach C., Sekiguchi K., Semple C.A., Seno S., Sessa L., Sheng Y., Shibata Y., Shimada H., Shimada K., Silva D., Sinclair B., Sperling S., Stupka E., Sugiura K., Sultana R., Takenaka Y., Taki K., Tammoja K., Tan S.L., Tang S., Taylor M.S., Tegner J., Teichmann S.A., Ueda H.R., van Nimwegen E., Verardo R., Wei C.L., Yagi K., Yamanishi H., Zabarovsky E., Zhu S., Zimmer A., Hide W., Bult C., Grimmond S.M., Teasdale R.D., Liu E.T., Brusic V., Quackenbush J., Wahlestedt C., Mattick J.S., Hume D.A., Kai C., Sasaki D., Tomaru Y., Fukuda S., Kanamori-Katayama M., Suzuki M., Aoki J., Arakawa T., Iida J., Imamura K., Itoh M., Kato T., Kawaji H., Kawagashira N., Kawashima T., Kojima M., Kondo S., Konno H., Nakano K., Ninomiya N., Nishio T., Okada M., Plessy C., Shibata K., Shiraki T., Suzuki S., Tagami M., Waki K., Watahiki A., Okamura-Oho Y., Suzuki H., Kawai J., Hayashizaki Y."
        assert reference.title == "The transcriptional landscape of the mammalian genome."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "16141072")
        assert reference.references[1] == ("DOI", "10.1126/science.1112014")
        reference = record.references[2]
        assert reference.authors == "Lubec G., Kang S.U."
        assert reference.title == ""
        assert len(reference.references) == 0
        reference = record.references[3]
        assert reference.authors == "Boal F., Le Pevelen S., Cziepluch C., Scotti P., Lang J."
        assert reference.title == "Cysteine-string protein isoform beta (Cspbeta) is targeted to the trans-Golgi network as a non-palmitoylated CSP in clonal beta-cells."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "17034881")
        assert reference.references[1] == ("DOI", "10.1016/j.bbamcr.2006.08.054")
        reference = record.references[4]
        assert reference.authors == "Lee J., Xu Y., Chen Y., Sprung R., Kim S.C., Xie S., Zhao Y."
        assert reference.title == "Mitochondrial phosphoproteome revealed by an improved IMAC method and MS/MS/MS."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "17208939")
        assert reference.references[1] == ("DOI", "10.1074/mcp.m600218-mcp200")
        reference = record.references[5]
        assert reference.authors == "Villen J., Beausoleil S.A., Gerber S.A., Gygi S.P."
        assert reference.title == "Large-scale phosphorylation analysis of mouse liver."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "17242355")
        assert reference.references[1] == ("DOI", "10.1073/pnas.0609836104")
        reference = record.references[6]
        assert reference.authors == "Zhou H., Ye M., Dong J., Han G., Jiang X., Wu R., Zou H."
        assert reference.title == "Specific phosphopeptide enrichment with immobilized titanium ion affinity chromatography adsorbent for phosphoproteome analysis."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "18630941")
        assert reference.references[1] == ("DOI", "10.1021/pr800223m")
        reference = record.references[7]
        assert reference.authors == "Trost M., English L., Lemieux S., Courcelles M., Desjardins M., Thibault P."
        assert reference.title == "The phagosomal proteome in interferon-gamma-activated macrophages."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "19144319")
        assert reference.references[1] == ("DOI", "10.1016/j.immuni.2008.11.006")
        reference = record.references[8]
        assert reference.authors == "Sweet S.M., Bailey C.M., Cunningham D.L., Heath J.K., Cooper H.J."
        assert reference.title == "Large scale localization of protein phosphorylation by use of electron capture dissociation mass spectrometry."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "19131326")
        assert reference.references[1] == ("DOI", "10.1074/mcp.m800451-mcp200")
        reference = record.references[9]
        assert reference.authors == "Huttlin E.L., Jedrychowski M.P., Elias J.E., Goswami T., Rad R., Beausoleil S.A., Villen J., Haas W., Sowa M.E., Gygi S.P."
        assert reference.title == "A tissue-specific atlas of mouse protein phosphorylation and expression."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "21183079")
        assert reference.references[1] == ("DOI", "10.1016/j.cell.2010.12.001")
        reference = record.references[10]
        assert reference.authors == "Boal F., Laguerre M., Milochau A., Lang J., Scotti P.A."
        assert reference.title == "A charged prominence in the linker domain of the cysteine-string protein Cspalpha mediates its regulated interaction with the calcium sensor synaptotagmin 9 during exocytosis."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "20847230")
        assert reference.references[1] == ("DOI", "10.1096/fj.09-152033")
        reference = record.references[11]
        assert reference.authors == "Sharma M., Burre J., Bronk P., Zhang Y., Xu W., Suedhof T.C."
        assert reference.title == "CSPalpha knockout causes neurodegeneration by impairing SNAP-25 function."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "22187053")
        assert reference.references[1] == ("DOI", "10.1038/emboj.2011.467")
        reference = record.references[12]
        assert reference.authors == "Lemonidis K., Gorleku O.A., Sanchez-Perez M.C., Grefen C., Chamberlain L.H."
        assert reference.title == "The Golgi S-acylation machinery comprises zDHHC enzymes with major differences in substrate affinity and S-acylation activity."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "25253725")
        assert reference.references[1] == ("DOI", "10.1091/mbc.e14-06-1169")
        reference = record.references[13]
        assert reference.authors == "Lemonidis K., Sanchez-Perez M.C., Chamberlain L.H."
        assert reference.title == "Identification of a novel sequence motif recognized by the ankyrin repeat domain of zDHHC17/13 S-acyltransferases."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "26198635")
        assert reference.references[1] == ("DOI", "10.1074/jbc.m115.657668")
        reference = record.references[14]
        assert reference.authors == "Burgoyne R.D., Morgan A."
        assert reference.title == "Cysteine string protein (CSP) and its role in preventing neurodegeneration."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "25800794")
        assert reference.references[1] == ("DOI", "10.1016/j.semcdb.2015.03.008")
        reference = record.references[15]
        assert reference.authors == "RIKEN structural genomics initiative (RSGI)"
        assert reference.title == "Solution structure of J-domain from mouse DnaJ subfamily C member 5."
        assert len(reference.references) == 0

        # Check the two parsers agree on the essentials
        assert seq_record.seq == record.sequence
        assert seq_record.description == record.description
        assert seq_record.name == record.entry_name
        assert seq_record.id in record.accessions

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        with open(datafile) as test_handle:
            records = list(SeqIO.parse(test_handle, "swiss"))

        assert len(records) == 1
        assert isinstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        assert records[0].seq == seq_record.seq
        assert records[0].description == seq_record.description
        assert records[0].name == seq_record.name
        assert records[0].id == seq_record.id

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        assert len(records) == 1
        assert isinstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        assert records[0].sequence == record.sequence
        assert records[0].description == record.description
        assert records[0].entry_name == record.entry_name
        assert records[0].accessions == record.accessions

    def test_P62258(self):
        """Parsing SwissProt file P62258."""
        filename = "P62258.txt"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        assert isinstance(seq_record, SeqRecord)

        assert seq_record.id == "P62258"
        assert seq_record.name == "1433E_HUMAN"
        assert seq_record.description == "RecName: Full=14-3-3 protein epsilon; Short=14-3-3E;"
        assert repr(seq_record.seq) == "Seq('MDDREDLVYQAKLAEQAERYDEMVESMKKVAGMDVELTVEERNLLSVAYKNVIG...ENQ')"

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        assert record.entry_name == "1433E_HUMAN"
        assert record.accessions == [
                "P62258",
                "B3KY71",
                "D3DTH5",
                "P29360",
                "P42655",
                "Q4VJB6",
                "Q53XZ5",
                "Q63631",
                "Q7M4R4",
            ]
        assert record.organism_classification == [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Euteleostomi",
                "Mammalia",
                "Eutheria",
                "Euarchontoglires",
                "Primates",
                "Haplorrhini",
                "Catarrhini",
                "Hominidae",
                "Homo",
            ]
        assert record.seqinfo == (255, 29174, "07817CCBD1F75B26")

        assert len(record.features) == 32
        assert len(record.references) == 46
        reference = record.references[0]
        assert reference.authors == "Conklin D.S., Galaktionov K., Beach D."
        assert reference.title == "14-3-3 proteins associate with cdc25 phosphatases."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "7644510")
        assert reference.references[1] == ("DOI", "10.1073/pnas.92.17.7892")
        reference = record.references[1]
        assert reference.authors == "Chong S.S., Tanigami A., Roschke A.V., Ledbetter D.H."
        assert reference.title == "14-3-3 epsilon has no homology to LIS1 and lies telomeric to it on chromosome 17p13.3 outside the Miller-Dieker syndrome chromosome region."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "8858348")
        assert reference.references[1] == ("DOI", "10.1101/gr.6.8.735")
        reference = record.references[2]
        assert reference.authors == "Jin D.-Y., Lyu M.S., Kozak C.A., Jeang K.-T."
        assert reference.title == "Function of 14-3-3 proteins."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "8684458")
        assert reference.references[1] == ("DOI", "10.1038/382308a0")
        reference = record.references[3]
        assert reference.authors == "Han D., Ye G., Liu T., Chen C., Yang X., Wan B., Pan Y., Yu L."
        assert reference.title == "Functional identification of a novel 14-3-3 epsilon splicing variant suggests dimerization is not necessary for 14-3-3 epsilon to inhibit UV-induced apoptosis."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "20417184")
        assert reference.references[1] == ("DOI", "10.1016/j.bbrc.2010.04.104")
        reference = record.references[4]
        assert reference.authors == "Luk S.C.W., Lee C.Y., Waye M.M.Y."
        assert reference.title == "Sequence determination of human epsilon 14-3-3 protein."
        assert len(reference.references) == 0
        reference = record.references[5]
        assert reference.authors == "Tanigami A., Chong S.S., Ledbetter D.H."
        assert reference.title == "14-3-3 epsilon genomic sequence."
        assert len(reference.references) == 0
        reference = record.references[6]
        assert reference.authors == "Ota T., Suzuki Y., Nishikawa T., Otsuki T., Sugiyama T., Irie R., Wakamatsu A., Hayashi K., Sato H., Nagai K., Kimura K., Makita H., Sekine M., Obayashi M., Nishi T., Shibahara T., Tanaka T., Ishii S., Yamamoto J., Saito K., Kawai Y., Isono Y., Nakamura Y., Nagahari K., Murakami K., Yasuda T., Iwayanagi T., Wagatsuma M., Shiratori A., Sudo H., Hosoiri T., Kaku Y., Kodaira H., Kondo H., Sugawara M., Takahashi M., Kanda K., Yokoi T., Furuya T., Kikkawa E., Omura Y., Abe K., Kamihara K., Katsuta N., Sato K., Tanikawa M., Yamazaki M., Ninomiya K., Ishibashi T., Yamashita H., Murakawa K., Fujimori K., Tanai H., Kimata M., Watanabe M., Hiraoka S., Chiba Y., Ishida S., Ono Y., Takiguchi S., Watanabe S., Yosida M., Hotuta T., Kusano J., Kanehori K., Takahashi-Fujii A., Hara H., Tanase T.-O., Nomura Y., Togiya S., Komai F., Hara R., Takeuchi K., Arita M., Imose N., Musashino K., Yuuki H., Oshima A., Sasaki N., Aotsuka S., Yoshikawa Y., Matsunawa H., Ichihara T., Shiohata N., Sano S., Moriya S., Momiyama H., Satoh N., Takami S., Terashima Y., Suzuki O., Nakagawa S., Senoh A., Mizoguchi H., Goto Y., Shimizu F., Wakebe H., Hishigaki H., Watanabe T., Sugiyama A., Takemoto M., Kawakami B., Yamazaki M., Watanabe K., Kumagai A., Itakura S., Fukuzumi Y., Fujimori Y., Komiyama M., Tashiro H., Tanigami A., Fujiwara T., Ono T., Yamada K., Fujii Y., Ozaki K., Hirao M., Ohmori Y., Kawabata A., Hikiji T., Kobatake N., Inagaki H., Ikema Y., Okamoto S., Okitani R., Kawakami T., Noguchi S., Itoh T., Shigeta K., Senba T., Matsumura K., Nakajima Y., Mizuno T., Morinaga M., Sasaki M., Togashi T., Oyama M., Hata H., Watanabe M., Komatsu T., Mizushima-Sugano J., Satoh T., Shirai Y., Takahashi Y., Nakagawa K., Okumura K., Nagase T., Nomura N., Kikuchi H., Masuho Y., Yamashita R., Nakai K., Yada T., Nakamura Y., Ohara O., Isogai T., Sugano S."
        assert reference.title == "Complete sequencing and characterization of 21,243 full-length human cDNAs."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "14702039")
        assert reference.references[1] == ("DOI", "10.1038/ng1285")
        reference = record.references[7]
        assert reference.authors == "Kalnine N., Chen X., Rolfs A., Halleck A., Hines L., Eisenstein S., Koundinya M., Raphael J., Moreira D., Kelley T., LaBaer J., Lin Y., Phelan M., Farmer A."
        assert reference.title == "Cloning of human full-length CDSs in BD Creator(TM) system donor vector."
        assert len(reference.references) == 0
        reference = record.references[8]
        assert reference.authors == "Mural R.J., Istrail S., Sutton G.G., Florea L., Halpern A.L., Mobarry C.M., Lippert R., Walenz B., Shatkay H., Dew I., Miller J.R., Flanigan M.J., Edwards N.J., Bolanos R., Fasulo D., Halldorsson B.V., Hannenhalli S., Turner R., Yooseph S., Lu F., Nusskern D.R., Shue B.C., Zheng X.H., Zhong F., Delcher A.L., Huson D.H., Kravitz S.A., Mouchard L., Reinert K., Remington K.A., Clark A.G., Waterman M.S., Eichler E.E., Adams M.D., Hunkapiller M.W., Myers E.W., Venter J.C."
        assert reference.title == ""
        assert len(reference.references) == 0
        reference = record.references[9]
        assert reference.authors == "The MGC Project Team"
        assert reference.title == "The status, quality, and expansion of the NIH full-length cDNA project: the Mammalian Gene Collection (MGC)."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "15489334")
        assert reference.references[1] == ("DOI", "10.1101/gr.2596504")
        reference = record.references[10]
        assert reference.authors == "Gevaert K., Goethals M., Martens L., Van Damme J., Staes A., Thomas G.R., Vandekerckhove J."
        assert reference.title == "Exploring proteomes and analyzing protein processing by mass spectrometric identification of sorted N-terminal peptides."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "12665801")
        assert reference.references[1] == ("DOI", "10.1038/nbt810")
        reference = record.references[11]
        assert reference.authors == "Greninger A.L., Knudsen G.M., Betegon M., Burlingame A.L., DeRisi J.L."
        assert reference.title == "ACBD3 interaction with TBC1 domain 22 protein is differentially affected by enteroviral and kobuviral 3A protein binding."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "23572552")
        assert reference.references[1] == ("DOI", "10.1128/mbio.00098-13")
        reference = record.references[12]
        assert reference.authors == "Bienvenut W.V."
        assert reference.title == ""
        assert len(reference.references) == 0
        reference = record.references[13]
        assert reference.authors == "Stewart S., Sundaram M., Zhang Y., Lee J., Han M., Guan K.L."
        assert reference.title == "Kinase suppressor of Ras forms a multiprotein signaling complex and modulates MEK localization."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "10409742")
        assert reference.references[1] == ("DOI", "10.1128/mcb.19.8.5523")
        reference = record.references[14]
        assert reference.authors == "Demeter J., Medzihradszky D., Kha H., Goetzl E.J., Turck C.W."
        assert reference.title == "Isolation and partial characterization of the structures of fibroblast activating factor-related proteins from U937 cells."
        assert len(reference.references) == 1
        assert reference.references[0] == ("PubMed", "2026444")
        reference = record.references[15]
        assert reference.authors == "Lubec G., Afjehi-Sadat L."
        assert reference.title == ""
        assert len(reference.references) == 0
        reference = record.references[16]
        assert reference.authors == "Aoki H., Hayashi J., Moriyama M., Arakawa Y., Hino O."
        assert reference.title == "Hepatitis C virus core protein interacts with 14-3-3 protein and activates the kinase Raf-1."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "10644344")
        assert reference.references[1] == ("DOI", "10.1128/jvi.74.4.1736-1741.2000")
        reference = record.references[17]
        assert reference.authors == "Ganguly S., Gastel J.A., Weller J.L., Schwartz C., Jaffe H., Namboodiri M.A., Coon S.L., Hickman A.B., Rollag M., Obsil T., Beauverger P., Ferry G., Boutin J.A., Klein D.C."
        assert reference.title == "Role of a pineal cAMP-operated arylalkylamine N-acetyltransferase/14-3-3-binding switch in melatonin synthesis."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "11427721")
        assert reference.references[1] == ("DOI", "10.1073/pnas.141118798")
        reference = record.references[18]
        assert reference.authors == "Fujita N., Sato S., Katayama K., Tsuruo T."
        assert reference.title == "Akt-dependent phosphorylation of p27Kip1 promotes binding to 14-3-3 and cytoplasmic localization."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "12042314")
        assert reference.references[1] == ("DOI", "10.1074/jbc.m203668200")
        reference = record.references[19]
        assert reference.authors == "Wang X., Grammatikakis N., Siganou A., Calderwood S.K."
        assert reference.title == "Regulation of molecular chaperone gene transcription involves the serine phosphorylation, 14-3-3 epsilon binding, and cytoplasmic sequestration of heat shock factor 1."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "12917326")
        assert reference.references[1] == ("DOI", "10.1128/mcb.23.17.6013-6026.2003")
        reference = record.references[20]
        assert reference.authors == "Andersen J.S., Wilkinson C.J., Mayor T., Mortensen P., Nigg E.A., Mann M."
        assert reference.title == "Proteomic characterization of the human centrosome by protein correlation profiling."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "14654843")
        assert reference.references[1] == ("DOI", "10.1038/nature02166")
        reference = record.references[21]
        assert reference.authors == "Urschel S., Bassermann F., Bai R.Y., Munch S., Peschel C., Duyster J."
        assert reference.title == "Phosphorylation of grb10 regulates its interaction with 14-3-3."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "15722337")
        assert reference.references[1] == ("DOI", "10.1074/jbc.m501477200")
        reference = record.references[22]
        assert reference.authors == "Yoshida K., Yamaguchi T., Natsume T., Kufe D., Miki Y."
        assert reference.title == "JNK phosphorylation of 14-3-3 proteins regulates nuclear targeting of c-Abl in the apoptotic response to DNA damage."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "15696159")
        assert reference.references[1] == ("DOI", "10.1038/ncb1228")
        reference = record.references[23]
        assert reference.authors == "Gu Y.-M., Jin Y.-H., Choi J.-K., Baek K.-H., Yeo C.-Y., Lee K.-Y."
        assert reference.title == "Protein kinase A phosphorylates and regulates dimerization of 14-3-3 epsilon."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "16376338")
        assert reference.references[1] == ("DOI", "10.1016/j.febslet.2005.12.024")
        reference = record.references[24]
        assert reference.authors == "Chi A., Valencia J.C., Hu Z.-Z., Watabe H., Yamaguchi H., Mangini N.J., Huang H., Canfield V.A., Cheng K.C., Yang F., Abe R., Yamagishi S., Shabanowitz J., Hearing V.J., Wu C., Appella E., Hunt D.F."
        assert reference.title == "Proteomic and bioinformatic characterization of the biogenesis and function of melanosomes."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "17081065")
        assert reference.references[1] == ("DOI", "10.1021/pr060363j")
        reference = record.references[25]
        assert reference.authors == "Linde C.I., Di Leva F., Domi T., Tosatto S.C., Brini M., Carafoli E."
        assert reference.title == "Inhibitory interaction of the 14-3-3 proteins with ubiquitous (PMCA1) and tissue-specific (PMCA3) isoforms of the plasma membrane Ca2+ pump."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "18029012")
        assert reference.references[1] == ("DOI", "10.1016/j.ceca.2007.09.003")
        reference = record.references[26]
        assert reference.authors == "Brummer T., Larance M., Herrera Abreu M.T., Lyons R.J., Timpson P., Emmerich C.H., Fleuren E.D.G., Lehrbach G.M., Schramek D., Guilhaus M., James D.E., Daly R.J."
        assert reference.title == "Phosphorylation-dependent binding of 14-3-3 terminates signalling by the Gab2 docking protein."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "19172738")
        assert reference.references[1] == ("DOI", "10.1038/emboj.2008.159")
        reference = record.references[27]
        assert reference.authors == "Han G., Ye M., Zhou H., Jiang X., Feng S., Jiang X., Tian R., Wan D., Zou H., Gu J."
        assert reference.title == "Large-scale phosphoproteome analysis of human liver tissue by enrichment and fractionation of phosphopeptides with strong anion exchange chromatography."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "18318008")
        assert reference.references[1] == ("DOI", "10.1002/pmic.200700884")
        reference = record.references[28]
        assert reference.authors == "Gauci S., Helbig A.O., Slijper M., Krijgsveld J., Heck A.J., Mohammed S."
        assert reference.title == "Lys-N and trypsin cover complementary parts of the phosphoproteome in a refined SCX-based approach."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "19413330")
        assert reference.references[1] == ("DOI", "10.1021/ac9004309")
        reference = record.references[29]
        assert reference.authors == "Kajiwara Y., Buxbaum J.D., Grice D.E."
        assert reference.title == "SLITRK1 binds 14-3-3 and regulates neurite outgrowth in a phosphorylation-dependent manner."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "19640509")
        assert reference.references[1] == ("DOI", "10.1016/j.biopsych.2009.05.033")
        reference = record.references[30]
        assert reference.authors == "Jang S.W., Liu X., Fu H., Rees H., Yepes M., Levey A., Ye K."
        assert reference.title == "Interaction of Akt-phosphorylated SRPK2 with 14-3-3 mediates cell cycle and cell death in neurons."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "19592491")
        assert reference.references[1] == ("DOI", "10.1074/jbc.m109.026237")
        reference = record.references[31]
        assert reference.authors == "Choudhary C., Kumar C., Gnad F., Nielsen M.L., Rehman M., Walther T.C., Olsen J.V., Mann M."
        assert reference.title == "Lysine acetylation targets protein complexes and co-regulates major cellular functions."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "19608861")
        assert reference.references[1] == ("DOI", "10.1126/science.1175371")
        reference = record.references[32]
        assert reference.authors == "Olsen J.V., Vermeulen M., Santamaria A., Kumar C., Miller M.L., Jensen L.J., Gnad F., Cox J., Jensen T.S., Nigg E.A., Brunak S., Mann M."
        assert reference.title == "Quantitative phosphoproteomics reveals widespread full phosphorylation site occupancy during mitosis."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "20068231")
        assert reference.references[1] == ("DOI", "10.1126/scisignal.2000475")
        reference = record.references[33]
        assert reference.authors == "Burkard T.R., Planyavsky M., Kaupe I., Breitwieser F.P., Buerckstuemmer T., Bennett K.L., Superti-Furga G., Colinge J."
        assert reference.title == "Initial characterization of the human central proteome."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "21269460")
        assert reference.references[1] == ("DOI", "10.1186/1752-0509-5-17")
        reference = record.references[34]
        assert reference.authors == "Rigbolt K.T., Prokhorova T.A., Akimov V., Henningsen J., Johansen P.T., Kratchmarova I., Kassem M., Mann M., Olsen J.V., Blagoev B."
        assert reference.title == "System-wide temporal characterization of the proteome and phosphoproteome of human embryonic stem cell differentiation."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "21406692")
        assert reference.references[1] == ("DOI", "10.1126/scisignal.2001570")
        reference = record.references[35]
        assert reference.authors == "Zhou H., Di Palma S., Preisinger C., Peng M., Polat A.N., Heck A.J., Mohammed S."
        assert reference.title == "Toward a comprehensive characterization of a human cancer cell phosphoproteome."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "23186163")
        assert reference.references[1] == ("DOI", "10.1021/pr300630k")
        reference = record.references[36]
        assert reference.authors == "Bian Y., Song C., Cheng K., Dong M., Wang F., Huang J., Sun D., Wang L., Ye M., Zou H."
        assert reference.title == "An enzyme assisted RP-RPLC approach for in-depth analysis of human liver phosphoproteome."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "24275569")
        assert reference.references[1] == ("DOI", "10.1016/j.jprot.2013.11.014")
        reference = record.references[37]
        assert reference.authors == "Yuasa K., Ota R., Matsuda S., Isshiki K., Inoue M., Tsuji A."
        assert reference.title == "Suppression of death-associated protein kinase 2 by interaction with 14-3-3 proteins."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "26047703")
        assert reference.references[1] == ("DOI", "10.1016/j.bbrc.2015.05.105")
        reference = record.references[38]
        assert reference.authors == "Kulasekaran G., Nossova N., Marat A.L., Lund I., Cremer C., Ioannou M.S., McPherson P.S."
        assert reference.title == "Phosphorylation-dependent regulation of Connecdenn/DENND1 guanine nucleotide exchange factors."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "26055712")
        assert reference.references[1] == ("DOI", "10.1074/jbc.m115.636712")
        reference = record.references[39]
        assert reference.authors == "Gao K., Tang W., Li Y., Zhang P., Wang D., Yu L., Wang C., Wu D."
        assert reference.title == "Front-signal-dependent accumulation of the RHOA inhibitor FAM65B at leading edges polarizes neutrophils."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "25588844")
        assert reference.references[1] == ("DOI", "10.1242/jcs.161497")
        reference = record.references[40]
        assert reference.authors == "Vaca Jacome A.S., Rabilloud T., Schaeffer-Reiss C., Rompais M., Ayoub D., Lane L., Bairoch A., Van Dorsselaer A., Carapito C."
        assert reference.title == "N-terminome analysis of the human mitochondrial proteome."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "25944712")
        assert reference.references[1] == ("DOI", "10.1002/pmic.201400617")
        reference = record.references[41]
        assert reference.authors == "Masters S.L., Lagou V., Jeru I., Baker P.J., Van Eyck L., Parry D.A., Lawless D., De Nardo D., Garcia-Perez J.E., Dagley L.F., Holley C.L., Dooley J., Moghaddas F., Pasciuto E., Jeandel P.Y., Sciot R., Lyras D., Webb A.I., Nicholson S.E., De Somer L., van Nieuwenhove E., Ruuth-Praz J., Copin B., Cochet E., Medlej-Hashim M., Megarbane A., Schroder K., Savic S., Goris A., Amselem S., Wouters C., Liston A."
        assert reference.title == "Familial autoinflammation with neutrophilic dermatosis reveals a regulatory mechanism of pyrin activation."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "27030597")
        assert reference.references[1] == ("DOI", "10.1126/scitranslmed.aaf1471")
        reference = record.references[42]
        assert reference.authors == "Hendriks I.A., Lyon D., Young C., Jensen L.J., Vertegaal A.C., Nielsen M.L."
        assert reference.title == "Site-specific mapping of the human SUMO proteome reveals co-modification with phosphorylation."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "28112733")
        assert reference.references[1] == ("DOI", "10.1038/nsmb.3366")
        reference = record.references[43]
        assert reference.authors == "Sonntag T., Ostojic J., Vaughan J.M., Moresco J.J., Yoon Y.S., Yates J.R. III, Montminy M."
        assert reference.title == "Mitogenic Signals Stimulate the CREB Coactivator CRTC3 through PP2A Recruitment."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "30611118")
        assert reference.references[1] == ("DOI", "10.1016/j.isci.2018.12.012")
        reference = record.references[44]
        assert reference.authors == "Chen J., Ou Y., Yang Y., Li W., Xu Y., Xie Y., Liu Y."
        assert reference.title == "KLHL22 activates amino-acid-dependent mTORC1 signalling to promote tumorigenesis and ageing."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "29769719")
        assert reference.references[1] == ("DOI", "10.1038/s41586-018-0128-9")
        reference = record.references[45]
        assert reference.authors == "Yang X., Lee W.H., Sobott F., Papagrigoriou E., Robinson C.V., Grossmann J.G., Sundstroem M., Doyle D.A., Elkins J.M."
        assert reference.title == "Structural basis for protein-protein interactions in the 14-3-3 protein family."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "17085597")
        assert reference.references[1] == ("DOI", "10.1073/pnas.0605779103")

        # Check the two parsers agree on the essentials
        assert seq_record.seq == record.sequence
        assert seq_record.description == record.description
        assert seq_record.name == record.entry_name
        assert seq_record.id in record.accessions

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        with open(datafile) as test_handle:
            records = list(SeqIO.parse(test_handle, "swiss"))

        assert len(records) == 1
        assert isinstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        assert records[0].seq == seq_record.seq
        assert records[0].description == seq_record.description
        assert records[0].name == seq_record.name
        assert records[0].id == seq_record.id

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        assert len(records) == 1
        assert isinstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        assert records[0].sequence == record.sequence
        assert records[0].description == record.description
        assert records[0].entry_name == record.entry_name
        assert records[0].accessions == record.accessions

    def test_P0A186(self):
        """Parsing SwissProt file P0A186.txt."""
        filename = "P0A186.txt"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        assert isinstance(seq_record, SeqRecord)

        assert seq_record.id == "P0A186"
        assert seq_record.name == "NDOA_PSEU8"
        assert seq_record.description == "RecName: Full=Naphthalene 1,2-dioxygenase system, ferredoxin component {ECO:0000303|PubMed:8226631};"
        assert repr(seq_record.seq) == "Seq('MTVKWIEAVALSDILEGDVLGVTVEGKELALYEVEGEIYATDNLCTHGSARMSD...DLS')"

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        assert record.entry_name == "NDOA_PSEU8"
        assert record.accessions == ["P0A186", "O07829", "P23082", "Q52123"]
        assert record.organism_classification == ["Bacteria", "Proteobacteria"]
        assert record.seqinfo == (104, 11446, "475625DCC3EDCD41")

        assert len(record.features) == 7
        feature = record.features[0]
        assert feature.type == "INIT_MET"
        assert feature.location.start == 0
        assert feature.location.end == 1
        assert feature.qualifiers["note"] == "Removed"
        assert feature.qualifiers["evidence"] == "ECO:0000250|UniProtKB:P0A185"
        feature = record.features[1]
        assert feature.id == "PRO_0000201694"
        assert feature.type == "CHAIN"
        assert feature.location.start == 1
        assert feature.location.end == 104
        assert feature.qualifiers["note"] == "Naphthalene 1,2-dioxygenase system, ferredoxin component"
        feature = record.features[2]
        assert feature.type == "DOMAIN"
        assert feature.location.start == 5
        assert feature.location.end == 101
        assert feature.qualifiers["note"] == "Rieske"
        assert feature.qualifiers["evidence"] == "ECO:0000255|PROSITE-ProRule:PRU00628"
        feature = record.features[3]
        assert feature.type == "METAL"
        assert feature.location.start == 44
        assert feature.location.end == 45
        assert feature.qualifiers["note"] == "Iron-sulfur (2Fe-2S)"
        assert feature.qualifiers["evidence"] == "ECO:0000250|UniProtKB:P0A185,ECO:0000255|PROSITE-ProRule:PRU00628"
        assert feature.id is None
        feature = record.features[4]
        assert feature.type == "METAL"
        assert feature.location.start == 46
        assert feature.location.end == 47
        assert feature.qualifiers["note"] == "Iron-sulfur (2Fe-2S); via pros nitrogen"
        assert feature.qualifiers["evidence"] == "ECO:0000250|UniProtKB:P0A185,ECO:0000255|PROSITE-ProRule:PRU00628"
        assert feature.id is None
        feature = record.features[5]
        assert feature.type == "METAL"
        assert feature.location.start == 63
        assert feature.location.end == 64
        assert feature.qualifiers["note"] == "Iron-sulfur (2Fe-2S)"
        assert feature.qualifiers["evidence"] == "ECO:0000250|UniProtKB:P0A185,ECO:0000255|PROSITE-ProRule:PRU00628"
        assert feature.id is None
        feature = record.features[6]
        assert feature.type == "METAL"
        assert feature.location.start == 66
        assert feature.location.end == 67
        assert feature.qualifiers["note"] == "Iron-sulfur (2Fe-2S); via pros nitrogen"
        assert feature.id is None

        assert len(record.references) == 1
        assert record.references[0].authors == "Denome S.A., Stanley D.C., Olson E.S., Young K.D."
        assert record.references[0].title == "Metabolism of dibenzothiophene and naphthalene in Pseudomonas strains: complete DNA sequence of an upper naphthalene catabolic pathway."
        assert len(record.references[0].references) == 2
        assert record.references[0].references[0] == ("PubMed", "8226631")
        assert record.references[0].references[1] == ("DOI", "10.1128/jb.175.21.6890-6901.1993")

        # Check the two parsers agree on the essentials
        assert seq_record.seq == record.sequence
        assert seq_record.description == record.description
        assert seq_record.name == record.entry_name
        assert seq_record.id in record.accessions

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        with open(datafile) as test_handle:
            records = list(SeqIO.parse(test_handle, "swiss"))

        assert len(records) == 1
        assert isinstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        assert records[0].seq == seq_record.seq
        assert records[0].description == seq_record.description
        assert records[0].name == seq_record.name
        assert records[0].id == seq_record.id

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        assert len(records) == 1
        assert isinstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        assert records[0].sequence == record.sequence
        assert records[0].description == record.description
        assert records[0].entry_name == record.entry_name
        assert records[0].accessions == record.accessions

    def test_P68308(self):
        """Parsing SwissProt file P68308.txt."""
        filename = "P68308.txt"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        assert isinstance(seq_record, SeqRecord)

        assert seq_record.id == "P68308"
        assert seq_record.name == "NU3M_BALPH"
        assert seq_record.description == "RecName: Full=NADH-ubiquinone oxidoreductase chain 3 {ECO:0000250|UniProtKB:P03897}; EC=7.1.1.2 {ECO:0000250|UniProtKB:P03897}; AltName: Full=NADH dehydrogenase subunit 3;"
        assert repr(seq_record.seq) == "Seq('MNLLLTLLTNTTLALLLVFIAFWLPQLNVYAEKTSPYECGFDPMGSARLPFSMK...WAE')"

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        assert record.entry_name == "NU3M_BALPH"
        assert record.accessions == ["P68308", "P24973"]
        assert record.organism_classification == [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Euteleostomi",
                "Mammalia",
                "Eutheria",
                "Laurasiatheria",
                "Artiodactyla",
                "Whippomorpha",
                "Cetacea",
                "Mysticeti",
                "Balaenopteridae",
                "Balaenoptera",
            ]
        assert record.seqinfo == (115, 13022, "405197D2F5D0AC4B")

        assert len(record.features) == 4
        feature = record.features[0]
        assert feature.type == "CHAIN"
        assert feature.location.start == 0
        assert feature.location.end == 115
        assert feature.qualifiers["note"] == "NADH-ubiquinone oxidoreductase chain 3"
        feature = record.features[1]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 2
        assert feature.location.end == 23
        assert feature.qualifiers["note"] == "Helical"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[2]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 54
        assert feature.location.end == 75
        assert feature.qualifiers["note"] == "Helical"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[3]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 83
        assert feature.location.end == 104
        assert feature.qualifiers["note"] == "Helical"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        assert feature.id is None

        assert len(record.references) == 1
        reference = record.references[0]
        assert reference.authors == "Arnason U., Gullberg A., Widegren B."
        assert reference.title == "The complete nucleotide sequence of the mitochondrial DNA of the fin whale, Balaenoptera physalus."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "1779436")
        assert reference.references[1] == ("DOI", "10.1007/bf02102808")

        # Check the two parsers agree on the essentials
        assert seq_record.seq == record.sequence
        assert seq_record.description == record.description
        assert seq_record.name == record.entry_name
        assert seq_record.id in record.accessions

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        with open(datafile) as test_handle:
            records = list(SeqIO.parse(test_handle, "swiss"))

        assert len(records) == 1
        assert isinstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        assert records[0].seq == seq_record.seq
        assert records[0].description == seq_record.description
        assert records[0].name == seq_record.name
        assert records[0].id == seq_record.id

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        assert len(records) == 1
        assert isinstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        assert records[0].sequence == record.sequence
        assert records[0].description == record.description
        assert records[0].entry_name == record.entry_name
        assert records[0].accessions == record.accessions

    def test_P39896(self):
        """Parsing SwissProt file P39896.txt."""
        filename = "P39896.txt"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        assert isinstance(seq_record, SeqRecord)

        assert seq_record.id == "P39896"
        assert seq_record.name == "TCMO_STRGA"
        assert seq_record.description == "RecName: Full=Tetracenomycin polyketide synthesis 8-O-methyl transferase TcmO; EC=2.1.1.-;"
        assert repr(seq_record.seq) == "Seq('MTPHTHVRGPGDILQLTMAFYGSRALISAVELDLFTLLAGKPLPLGELCERAGI...KPR')"

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        assert record.entry_name == "TCMO_STRGA"
        assert record.accessions == ["P39896"]
        assert record.organism_classification == [
                "Bacteria",
                "Actinobacteria",
                "Streptomycetales",
                "Streptomycetaceae",
                "Streptomyces",
            ]
        assert record.seqinfo == (339, 37035, "B228B66B24217F80")

        assert len(record.features) == 4
        feature = record.features[0]
        assert feature.type == "CHAIN"
        assert feature.location.start == 0
        assert feature.location.end == 339
        assert feature.qualifiers["note"] == "Tetracenomycin polyketide synthesis 8-O-methyl transferase TcmO"
        feature = record.features[1]
        assert feature.type == "ACT_SITE"
        assert feature.location.start == 245
        assert feature.location.end == 246
        assert feature.qualifiers["note"] == "Proton acceptor"
        assert feature.qualifiers["evidence"] == "ECO:0000255|PROSITE-ProRule:PRU01020"
        feature = record.features[2]
        assert feature.type == "BINDING"
        assert feature.location.start == 199
        assert feature.location.end == 200
        assert feature.qualifiers["ligand"] == "S-adenosyl-L-methionine"
        assert feature.qualifiers["ligand_id"] == "ChEBI:CHEBI:59789"
        assert feature.qualifiers["evidence"] == "ECO:0000255|PROSITE-ProRule:PRU01020"
        feature = record.features[3]
        assert feature.type == "BINDING"
        assert feature.location.start == 225
        assert feature.location.end == 228
        assert feature.qualifiers["ligand"] == "S-adenosyl-L-methionine"
        assert feature.qualifiers["ligand_id"] == "ChEBI:CHEBI:59789"
        assert feature.qualifiers["evidence"] == "ECO:0000255|PROSITE-ProRule:PRU01020"
        assert feature.id is None

        assert len(record.references) == 1
        reference = record.references[0]
        assert reference.authors == "Summers R.G., Wendt-Pienkowski E., Motamedi H., Hutchinson C.R."
        assert reference.title == "Nucleotide sequence of the tcmII-tcmIV region of the tetracenomycin C biosynthetic gene cluster of Streptomyces glaucescens and evidence that the tcmN gene encodes a multifunctional cyclase-dehydratase-O-methyl transferase."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "1548230")
        assert reference.references[1] == ("DOI", "10.1128/jb.174.6.1810-1820.1992")

        # Check the two parsers agree on the essentials
        assert seq_record.seq == record.sequence
        assert seq_record.description == record.description
        assert seq_record.name == record.entry_name
        assert seq_record.id in record.accessions

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        with open(datafile) as test_handle:
            records = list(SeqIO.parse(test_handle, "swiss"))

        assert len(records) == 1
        assert isinstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        assert records[0].seq == seq_record.seq
        assert records[0].description == seq_record.description
        assert records[0].name == seq_record.name
        assert records[0].id == seq_record.id

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        assert len(records) == 1
        assert isinstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        assert records[0].sequence == record.sequence
        assert records[0].description == record.description
        assert records[0].entry_name == record.entry_name
        assert records[0].accessions == record.accessions

    def test_O95832(self):
        """Parsing SwissProt file O95832.txt."""
        filename = "O95832.txt"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        assert isinstance(seq_record, SeqRecord)

        assert seq_record.id == "O95832"
        assert seq_record.name == "CLD1_HUMAN"
        assert seq_record.description == "RecName: Full=Claudin-1; AltName: Full=Senescence-associated epithelial membrane protein;"
        assert repr(seq_record.seq) == "Seq('MANAGLQLLGFILAFLGWIGAIVSTALPQWRIYSYAGDNIVTAQAMYEGLWMSC...DYV')"

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        assert record.entry_name == "CLD1_HUMAN"
        assert record.accessions == ["O95832"]
        assert record.organism_classification == [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Euteleostomi",
                "Mammalia",
                "Eutheria",
                "Euarchontoglires",
                "Primates",
                "Haplorrhini",
                "Catarrhini",
                "Hominidae",
                "Homo",
            ]
        assert record.seqinfo == (211, 22744, "07269000E6C214F0")

        assert len(record.features) == 17
        feature = record.features[0]
        assert feature.type == "CHAIN"
        assert feature.location.start == 0
        assert feature.location.end == 211
        assert feature.qualifiers["note"] == "Claudin-1"
        feature = record.features[1]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 0
        assert feature.location.end == 7
        assert feature.qualifiers["note"] == "Cytoplasmic"
        assert feature.qualifiers["evidence"] == "ECO:0000250"
        feature = record.features[2]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 7
        assert feature.location.end == 28
        assert feature.qualifiers["note"] == "Helical"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[3]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 28
        assert feature.location.end == 81
        assert feature.qualifiers["note"] == "Extracellular"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[4]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 81
        assert feature.location.end == 102
        assert feature.qualifiers["note"] == "Helical"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[5]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 102
        assert feature.location.end == 115
        assert feature.qualifiers["note"] == "Cytoplasmic"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[6]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 115
        assert feature.location.end == 136
        assert feature.qualifiers["note"] == "Helical"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[7]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 136
        assert feature.location.end == 163
        assert feature.qualifiers["note"] == "Extracellular"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[8]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 163
        assert feature.location.end == 184
        assert feature.qualifiers["note"] == "Helical"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[9]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 184
        assert feature.location.end == 211
        assert feature.qualifiers["note"] == "Cytoplasmic"
        assert feature.qualifiers["evidence"] == "ECO:0000250"
        feature = record.features[10]
        assert feature.type == "REGION"
        assert feature.location.start == 191
        assert feature.location.end == 211
        assert feature.qualifiers["note"] == "Disordered"
        assert feature.qualifiers["evidence"] == "ECO:0000256|SAM:MobiDB-lite"
        feature = record.features[11]
        assert feature.type == "REGION"
        assert feature.location.start == 209
        assert feature.location.end == 211
        assert feature.qualifiers["note"] == "Interactions with TJP1, TJP2, TJP3 and PATJ"
        assert feature.qualifiers["evidence"] == "ECO:0000250"
        feature = record.features[12]
        assert feature.type == "DISULFID"
        assert feature.location.start == 53
        assert feature.location.end == 64
        assert feature.qualifiers["evidence"] == "ECO:0000250"
        feature = record.features[13]
        assert feature.type == "MUTAGEN"
        assert feature.location.start == 31
        assert feature.location.end == 32
        assert feature.qualifiers["note"] == "I->M: Loss of HCV receptor activity. Significant loss of interaction with CD81. Reduced interaction with OCLN."
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:17325668,ECO:0000269|PubMed:20375010"
        feature = record.features[14]
        assert feature.type == "MUTAGEN"
        assert feature.location.start == 47
        assert feature.location.end == 48
        assert feature.qualifiers["note"] == "E->K: Loss of HCV receptor activity. Significant loss of interaction with CD81. Reduced interaction with OCLN. According to PubMed:17325668 no effect observed on HCV infection susceptibility in cell culture."
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:17325668,ECO:0000269|PubMed:20375010"
        feature = record.features[15]
        assert feature.type == "CONFLICT"
        assert feature.location.start == 61
        assert feature.location.end == 62
        assert feature.qualifiers["note"] == "I -> V (in Ref. 2; AAD22962)"
        assert feature.qualifiers["evidence"] == "ECO:0000305"
        feature = record.features[16]
        assert feature.type == "CONFLICT"
        assert feature.location.start == 134
        assert feature.location.end == 135
        assert feature.qualifiers["note"] == "V -> A (in Ref. 2; AAD22962)"
        assert feature.qualifiers["evidence"] == "ECO:0000305"
        assert feature.id is None

        assert len(record.references) == 17
        reference = record.references[0]
        assert reference.authors == "Swisshelm K.L., Machl A., Planitzer S., Robertson R., Kubbies M., Hosier S."
        assert reference.title == "SEMP1, a senescence-associated cDNA isolated from human mammary epithelial cells, is a member of an epithelial membrane protein superfamily."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "9931503")
        assert reference.references[1] == ("DOI", "10.1016/s0378-1119(98)00553-8")
        reference = record.references[1]
        assert reference.authors == "Mitic L.M., Anderson J.M."
        assert reference.title == "Human claudin-1 isolated from Caco-2 mRNA."
        assert len(reference.references) == 0
        reference = record.references[2]
        assert reference.authors == "Halford S., Spencer P., Greenwood J., Winton H., Hunt D.M., Adamson P."
        assert reference.title == "Assignment of claudin-1 (CLDN1) to human chromosome 3q28-->q29 with somatic cell hybrids."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "10828592")
        assert reference.references[1] == ("DOI", "10.1159/000015553")
        reference = record.references[3]
        assert reference.authors == "Kraemer F., White K., Kubbies M., Swisshelm K.L., Weber B.H.F."
        assert reference.title == "Genomic organization of claudin-1 and its assessment in hereditary and sporadic breast cancer."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "11071387")
        assert reference.references[1] == ("DOI", "10.1007/s004390000375")
        reference = record.references[4]
        assert reference.authors == "Clark H.F., Gurney A.L., Abaya E., Baker K., Baldwin D.T., Brush J., Chen J., Chow B., Chui C., Crowley C., Currell B., Deuel B., Dowd P., Eaton D., Foster J.S., Grimaldi C., Gu Q., Hass P.E., Heldens S., Huang A., Kim H.S., Klimowski L., Jin Y., Johnson S., Lee J., Lewis L., Liao D., Mark M.R., Robbie E., Sanchez C., Schoenfeld J., Seshagiri S., Simmons L., Singh J., Smith V., Stinson J., Vagts A., Vandlen R.L., Watanabe C., Wieand D., Woods K., Xie M.-H., Yansura D.G., Yi S., Yu G., Yuan J., Zhang M., Zhang Z., Goddard A.D., Wood W.I., Godowski P.J., Gray A.M."
        assert reference.title == "The secreted protein discovery initiative (SPDI), a large-scale effort to identify novel human secreted and transmembrane proteins: a bioinformatics assessment."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "12975309")
        assert reference.references[1] == ("DOI", "10.1101/gr.1293003")
        reference = record.references[5]
        assert reference.authors == "Ebert L., Schick M., Neubert P., Schatten R., Henze S., Korn B."
        assert reference.title == "Cloning of human full open reading frames in Gateway(TM) system entry vector (pDONR201)."
        assert len(reference.references) == 0
        reference = record.references[6]
        assert reference.authors == "The MGC Project Team"
        assert reference.title == "The status, quality, and expansion of the NIH full-length cDNA project: the Mammalian Gene Collection (MGC)."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "15489334")
        assert reference.references[1] == ("DOI", "10.1101/gr.2596504")
        reference = record.references[7]
        assert reference.authors == "Hadj-Rabia S., Baala L., Vabres P., Hamel-Teillac D., Jacquemin E., Fabre M., Lyonnet S., De Prost Y., Munnich A., Hadchouel M., Smahi A."
        assert reference.title == "Claudin-1 gene mutations in neonatal sclerosing cholangitis associated with ichthyosis: a tight junction disease."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "15521008")
        assert reference.references[1] == ("DOI", "10.1053/j.gastro.2004.07.022")
        reference = record.references[8]
        assert reference.authors == "Feldmeyer L., Huber M., Fellmann F., Beckmann J.S., Frenk E., Hohl D."
        assert reference.title == "Confirmation of the origin of NISCH syndrome."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "16619213")
        assert reference.references[1] == ("DOI", "10.1002/humu.20333")
        reference = record.references[9]
        assert reference.authors == "Evans M.J., von Hahn T., Tscherne D.M., Syder A.J., Panis M., Wolk B., Hatziioannou T., McKeating J.A., Bieniasz P.D., Rice C.M."
        assert reference.title == "Claudin-1 is a hepatitis C virus co-receptor required for a late step in entry."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "17325668")
        assert reference.references[1] == ("DOI", "10.1038/nature05654")
        reference = record.references[10]
        assert reference.authors == "Harris H.J., Davis C., Mullins J.G., Hu K., Goodall M., Farquhar M.J., Mee C.J., McCaffrey K., Young S., Drummer H., Balfe P., McKeating J.A."
        assert reference.title == "Claudin association with CD81 defines hepatitis C virus entry."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "20375010")
        assert reference.references[1] == ("DOI", "10.1074/jbc.m110.104836")
        reference = record.references[11]
        assert reference.authors == "Burkard T.R., Planyavsky M., Kaupe I., Breitwieser F.P., Buerckstuemmer T., Bennett K.L., Superti-Furga G., Colinge J."
        assert reference.title == "Initial characterization of the human central proteome."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "21269460")
        assert reference.references[1] == ("DOI", "10.1186/1752-0509-5-17")
        reference = record.references[12]
        assert reference.authors == "Lupberger J., Zeisel M.B., Xiao F., Thumann C., Fofana I., Zona L., Davis C., Mee C.J., Turek M., Gorke S., Royer C., Fischer B., Zahid M.N., Lavillette D., Fresquet J., Cosset F.L., Rothenberg S.M., Pietschmann T., Patel A.H., Pessaux P., Doffoel M., Raffelsberger W., Poch O., McKeating J.A., Brino L., Baumert T.F."
        assert reference.title == "EGFR and EphA2 are host factors for hepatitis C virus entry and possible targets for antiviral therapy."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "21516087")
        assert reference.references[1] == ("DOI", "10.1038/nm.2341")
        reference = record.references[13]
        assert reference.authors == "Kirschner N., Rosenthal R., Furuse M., Moll I., Fromm M., Brandner J.M."
        assert reference.title == "Contribution of tight junction proteins to ion, macromolecule, and water barrier in keratinocytes."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "23407391")
        assert reference.references[1] == ("DOI", "10.1038/jid.2012.507")
        reference = record.references[14]
        assert reference.authors == "Che P., Tang H., Li Q."
        assert reference.title == "The interaction between claudin-1 and dengue viral prM/M protein for its entry."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "24074594")
        assert reference.references[1] == ("DOI", "10.1016/j.virol.2013.08.009")
        reference = record.references[15]
        assert reference.authors == "Bonander N., Jamshad M., Oberthuer D., Clare M., Barwell J., Hu K., Farquhar M.J., Stamataki Z., Harris H.J., Dierks K., Dafforn T.R., Betzel C., McKeating J.A., Bill R.M."
        assert reference.title == "Production, purification and characterization of recombinant, full-length human claudin-1."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "23704991")
        assert reference.references[1] == ("DOI", "10.1371/journal.pone.0064517")
        reference = record.references[16]
        assert reference.authors == "Douam F., Dao Thi V.L., Maurin G., Fresquet J., Mompelat D., Zeisel M.B., Baumert T.F., Cosset F.L., Lavillette D."
        assert reference.title == "Critical interaction between E1 and E2 glycoproteins determines binding and fusion properties of hepatitis C virus during cell entry."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "24038151")
        assert reference.references[1] == ("DOI", "10.1002/hep.26733")

        # Check the two parsers agree on the essentials
        assert seq_record.seq == record.sequence
        assert seq_record.description == record.description
        assert seq_record.name == record.entry_name
        assert seq_record.id in record.accessions

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        with open(datafile) as test_handle:
            records = list(SeqIO.parse(test_handle, "swiss"))

        assert len(records) == 1
        assert isinstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        assert records[0].seq == seq_record.seq
        assert records[0].description == seq_record.description
        assert records[0].name == seq_record.name
        assert records[0].id == seq_record.id

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        assert len(records) == 1
        assert isinstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        assert records[0].sequence == record.sequence
        assert records[0].description == record.description
        assert records[0].entry_name == record.entry_name
        assert records[0].accessions == record.accessions

    def test_P04439(self):
        """Parsing SwissProt file P04439.txt."""
        filename = "P04439.txt"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        assert isinstance(seq_record, SeqRecord)

        assert seq_record.id == "P04439"
        assert seq_record.name == "HLAA_HUMAN"
        assert seq_record.description == "RecName: Full=HLA class I histocompatibility antigen, A alpha chain; AltName: Full=Human leukocyte antigen A; Short=HLA-A; Flags: Precursor;"
        assert repr(seq_record.seq) == "Seq('MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDD...CKV')"

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        assert record.entry_name == "HLAA_HUMAN"
        assert record.accessions == [
                "P04439",
                "B1PKZ3",
                "O02939",
                "O02954",
                "O02955",
                "O02963",
                "O19509",
                "O19546",
                "O19598",
                "O19605",
                "O19606",
                "O19619",
                "O19647",
                "O19673",
                "O19687",
                "O19695",
                "O19756",
                "O19794",
                "O19795",
                "O43906",
                "O43907",
                "O46874",
                "O62921",
                "O62924",
                "O77937",
                "O77938",
                "O77964",
                "O78073",
                "O78171",
                "O98009",
                "O98010",
                "O98011",
                "O98137",
                "P01891",
                "P01892",
                "P05534",
                "P06338",
                "P10313",
                "P10314",
                "P10315",
                "P10316",
                "P13746",
                "P16188",
                "P16189",
                "P16190",
                "P18462",
                "P30443",
                "P30444",
                "P30445",
                "P30446",
                "P30447",
                "P30448",
                "P30449",
                "P30450",
                "P30451",
                "P30452",
                "P30453",
                "P30454",
                "P30455",
                "P30456",
                "P30457",
                "P30458",
                "P30459",
                "P30512",
                "P30514",
                "P79505",
                "P79562",
                "P79563",
                "Q09160",
                "Q29680",
                "Q29747",
                "Q29835",
                "Q29837",
                "Q29838",
                "Q29899",
                "Q29908",
                "Q29909",
                "Q29910",
                "Q30208",
                "Q31623",
                "Q5S3G1",
                "Q65A82",
                "Q8MHM1",
                "Q8MHN9",
                "Q95352",
                "Q95355",
                "Q95362",
                "Q95377",
                "Q95380",
                "Q95IZ5",
                "Q9BCN0",
                "Q9BD15",
                "Q9BD19",
                "Q9GJE6",
                "Q9GJE7",
                "Q9GJE8",
                "Q9MW42",
                "Q9MY89",
                "Q9MYA3",
                "Q9MYA5",
                "Q9MYC4",
                "Q9MYE6",
                "Q9MYE9",
                "Q9MYG4",
                "Q9MYG5",
                "Q9MYI5",
                "Q9TP25",
                "Q9TPQ3",
                "Q9TPR8",
                "Q9TPX8",
                "Q9TPX9",
                "Q9TPY0",
                "Q9TQ24",
                "Q9TQE8",
                "Q9TQE9",
                "Q9TQF1",
                "Q9TQF5",
                "Q9TQF8",
                "Q9TQF9",
                "Q9TQG0",
                "Q9TQG5",
                "Q9TQG7",
                "Q9TQH5",
                "Q9TQI3",
                "Q9TQK5",
                "Q9TQM6",
                "Q9TQN5",
                "Q9TQP5",
                "Q9TQP6",
                "Q9TQP7",
                "Q9UIN1",
                "Q9UIN2",
                "Q9UIP7",
                "Q9UQU3",
                "Q9UQU6",
                "Q9UQU7",
            ]
        assert record.organism_classification == [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Euteleostomi",
                "Mammalia",
                "Eutheria",
                "Euarchontoglires",
                "Primates",
                "Haplorrhini",
                "Catarrhini",
                "Hominidae",
                "Homo",
            ]
        assert record.seqinfo == (365, 40841, "DEDFCEC4450E0580")

        assert len(record.features) == 161
        feature = record.features[0]
        assert feature.type == "SIGNAL"
        assert feature.location.start == 0
        assert feature.location.end == 24
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:92029"
        feature = record.features[1]
        assert feature.type == "CHAIN"
        assert feature.location.start == 24
        assert feature.location.end == 365
        assert feature.qualifiers["note"] == "HLA class I histocompatibility antigen, A alpha chain"
        feature = record.features[2]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 24
        assert feature.location.end == 308
        assert feature.qualifiers["note"] == "Extracellular"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[3]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 308
        assert feature.location.end == 332
        assert feature.qualifiers["note"] == "Helical"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[4]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 332
        assert feature.location.end == 365
        assert feature.qualifiers["note"] == "Cytoplasmic"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[5]
        assert feature.type == "DOMAIN"
        assert feature.location.start == 208
        assert feature.location.end == 295
        assert feature.qualifiers["note"] == "Ig-like C1-type"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[6]
        assert feature.type == "REGION"
        assert feature.location.start == 24
        assert feature.location.end == 114
        assert feature.qualifiers["note"] == "Alpha-1"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[7]
        assert feature.type == "REGION"
        assert feature.location.start == 114
        assert feature.location.end == 206
        assert feature.qualifiers["note"] == "Alpha-2"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[8]
        assert feature.type == "REGION"
        assert feature.location.start == 206
        assert feature.location.end == 298
        assert feature.qualifiers["note"] == "Alpha-3"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[9]
        assert feature.type == "REGION"
        assert feature.location.start == 298
        assert feature.location.end == 308
        assert feature.qualifiers["note"] == "Connecting peptide"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[10]
        assert feature.type == "REGION"
        assert feature.location.start == 338
        assert feature.location.end == 365
        assert feature.qualifiers["note"] == "Disordered"
        assert feature.qualifiers["evidence"] == "ECO:0000256|SAM:MobiDB-lite"
        feature = record.features[11]
        assert feature.type == "COMPBIAS"
        assert feature.location.start == 340
        assert feature.location.end == 359
        assert feature.qualifiers["note"] == "Polar residues"
        assert feature.qualifiers["evidence"] == "ECO:0000256|SAM:MobiDB-lite"
        feature = record.features[12]
        assert feature.type == "BINDING"
        assert feature.location.start == 30
        assert feature.location.end == 31
        assert feature.qualifiers["ligand"] == "a peptide antigen"
        assert feature.qualifiers["ligand_id"] == "ChEBI:CHEBI:166823"
        assert feature.qualifiers["ligand_label"] == "1"
        assert feature.qualifiers["ligand_note"] == "pathogen-derived peptide antigen"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:21943705"
        feature = record.features[13]
        assert feature.type == "BINDING"
        assert feature.location.start == 96
        assert feature.location.end == 97
        assert feature.qualifiers["ligand"] == "a peptide antigen"
        assert feature.qualifiers["ligand_id"] == "ChEBI:CHEBI:166823"
        assert feature.qualifiers["ligand_label"] == "1"
        assert feature.qualifiers["ligand_note"] == "pathogen-derived peptide antigen"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:21943705"
        feature = record.features[14]
        assert feature.type == "BINDING"
        assert feature.location.start == 107
        assert feature.location.end == 108
        assert feature.qualifiers["ligand"] == "a peptide antigen"
        assert feature.qualifiers["ligand_id"] == "ChEBI:CHEBI:166823"
        assert feature.qualifiers["ligand_label"] == "1"
        assert feature.qualifiers["ligand_note"] == "pathogen-derived peptide antigen"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:21943705"
        feature = record.features[15]
        assert feature.type == "BINDING"
        assert feature.location.start == 139
        assert feature.location.end == 140
        assert feature.qualifiers["ligand"] == "a peptide antigen"
        assert feature.qualifiers["ligand_id"] == "ChEBI:CHEBI:166823"
        assert feature.qualifiers["ligand_label"] == "2"
        assert feature.qualifiers["ligand_note"] == "self-peptide antigen"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:21543847"
        feature = record.features[16]
        assert feature.type == "BINDING"
        assert feature.location.start == 166
        assert feature.location.end == 167
        assert feature.qualifiers["ligand"] == "a peptide antigen"
        assert feature.qualifiers["ligand_id"] == "ChEBI:CHEBI:166823"
        assert feature.qualifiers["ligand_label"] == "1"
        assert feature.qualifiers["ligand_note"] == "pathogen-derived peptide antigen"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:21943705"
        feature = record.features[17]
        assert feature.type == "BINDING"
        assert feature.location.start == 169
        assert feature.location.end == 170
        assert feature.qualifiers["ligand"] == "a peptide antigen"
        assert feature.qualifiers["ligand_id"] == "ChEBI:CHEBI:166823"
        assert feature.qualifiers["ligand_label"] == "1"
        assert feature.qualifiers["ligand_note"] == "pathogen-derived peptide antigen"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:21943705"
        feature = record.features[18]
        assert feature.type == "BINDING"
        assert feature.location.start == 182
        assert feature.location.end == 183
        assert feature.qualifiers["ligand"] == "a peptide antigen"
        assert feature.qualifiers["ligand_id"] == "ChEBI:CHEBI:166823"
        assert feature.qualifiers["ligand_label"] == "1"
        assert feature.qualifiers["ligand_note"] == "pathogen-derived peptide antigen"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:21943705"
        feature = record.features[19]
        assert feature.type == "BINDING"
        assert feature.location.start == 182
        assert feature.location.end == 183
        assert feature.qualifiers["ligand"] == "a peptide antigen"
        assert feature.qualifiers["ligand_id"] == "ChEBI:CHEBI:166823"
        assert feature.qualifiers["ligand_label"] == "2"
        assert feature.qualifiers["ligand_note"] == "self-peptide antigen"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:21543847"
        feature = record.features[20]
        assert feature.type == "BINDING"
        assert feature.location.start == 194
        assert feature.location.end == 195
        assert feature.qualifiers["ligand"] == "a peptide antigen"
        assert feature.qualifiers["ligand_id"] == "ChEBI:CHEBI:166823"
        assert feature.qualifiers["ligand_label"] == "1"
        assert feature.qualifiers["ligand_note"] == "pathogen-derived peptide antigen"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:21943705"
        feature = record.features[21]
        assert feature.type == "MOD_RES"
        assert feature.location.start == 82
        assert feature.location.end == 83
        assert feature.qualifiers["note"] == "Sulfotyrosine"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[22]
        assert feature.type == "MOD_RES"
        assert feature.location.start == 342
        assert feature.location.end == 343
        assert feature.qualifiers["note"] == "Phosphoserine"
        assert feature.qualifiers["evidence"] == "ECO:0007744|PubMed:24275569"
        feature = record.features[23]
        assert feature.type == "MOD_RES"
        assert feature.location.start == 343
        assert feature.location.end == 344
        assert feature.qualifiers["note"] == "Phosphotyrosine"
        assert feature.qualifiers["evidence"] == "ECO:0007744|PubMed:24275569"
        feature = record.features[24]
        assert feature.type == "MOD_RES"
        assert feature.location.start == 348
        assert feature.location.end == 349
        assert feature.qualifiers["note"] == "Phosphoserine"
        assert feature.qualifiers["evidence"] == "ECO:0007744|PubMed:24275569"
        feature = record.features[25]
        assert feature.type == "MOD_RES"
        assert feature.location.start == 349
        assert feature.location.end == 350
        assert feature.qualifiers["note"] == "Phosphoserine"
        assert feature.qualifiers["evidence"] == "ECO:0007744|PubMed:24275569"
        feature = record.features[26]
        assert feature.type == "MOD_RES"
        assert feature.location.start == 351
        assert feature.location.end == 352
        assert feature.qualifiers["note"] == "Phosphoserine"
        assert feature.qualifiers["evidence"] == "ECO:0007744|PubMed:24275569"
        feature = record.features[27]
        assert feature.type == "MOD_RES"
        assert feature.location.start == 355
        assert feature.location.end == 356
        assert feature.qualifiers["note"] == "Phosphoserine"
        assert feature.qualifiers["evidence"] == "ECO:0007744|PubMed:23186163,ECO:0007744|PubMed:24275569"
        feature = record.features[28]
        assert feature.type == "MOD_RES"
        assert feature.location.start == 358
        assert feature.location.end == 359
        assert feature.qualifiers["note"] == "Phosphoserine"
        assert feature.qualifiers["evidence"] == "ECO:0007744|PubMed:23186163,ECO:0007744|PubMed:24275569"
        feature = record.features[29]
        assert feature.type == "CARBOHYD"
        assert feature.location.start == 109
        assert feature.location.end == 110
        assert feature.qualifiers["note"] == "N-linked (GlcNAc...) asparagine"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:19159218"
        feature = record.features[30]
        assert feature.type == "DISULFID"
        assert feature.location.start == 124
        assert feature.location.end == 188
        assert feature.qualifiers["evidence"] == "ECO:0000255|PROSITE-ProRule:PRU00114,ECO:0000269|PubMed:16041067, ECO:0000269|PubMed:19177349,ECO:0000269|PubMed:20844028, ECO:0000269|PubMed:21543847,ECO:0000269|PubMed:21943705, ECO:0000269|PubMed:26758806,ECO:0000269|PubMed:28250417, ECO:0000269|PubMed:7694806"
        feature = record.features[31]
        assert feature.type == "DISULFID"
        assert feature.location.start == 226
        assert feature.location.end == 283
        assert feature.qualifiers["evidence"] == "ECO:0000255|PROSITE-ProRule:PRU00114,ECO:0000269|PubMed:16041067, ECO:0000269|PubMed:19177349,ECO:0000269|PubMed:20844028, ECO:0000269|PubMed:21543847,ECO:0000269|PubMed:21943705, ECO:0000269|PubMed:26758806,ECO:0000269|PubMed:28250417, ECO:0000269|PubMed:7694806"
        feature = record.features[32]
        assert feature.type == "VAR_SEQ"
        assert feature.location.start == 175
        assert feature.location.end == 187
        assert feature.qualifiers["note"] == "EAEQLRAYLDGT -> AAEQQRAYLEGR (in isoform 2)"
        feature = record.features[33]
        assert feature.type == "VAR_SEQ"
        assert feature.location.start == 336
        assert feature.location.end == 337
        assert feature.qualifiers["note"] == "S -> SGGEGVK (in isoform 2)"
        feature = record.features[34]
        assert feature.type == "VARIANT"
        assert feature.location.start == 2
        assert feature.location.end == 3
        assert feature.qualifiers["note"] == "V -> I (in allele A*34:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115"
        feature = record.features[35]
        assert feature.type == "VARIANT"
        assert feature.location.start == 4
        assert feature.location.end == 5
        assert feature.qualifiers["note"] == "A -> P (in allele A*80:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791"
        feature = record.features[36]
        assert feature.type == "VARIANT"
        assert feature.location.start == 9
        assert feature.location.end == 10
        assert feature.qualifiers["note"] == "L -> V (in allele A*02:01, allele A*02:05, allele A*23:01, allele A*24:02, allele A*25:01, allele A*26:01, allele A*34:01, allele A*43:01, allele A*66:01, allele A*68:01 and allele A*69:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:1729171, ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2982951, ECO:0000269|PubMed:3496393,ECO:0000269|PubMed:3877632, ECO:0000269|PubMed:7836067,ECO:0000269|PubMed:8026990, ECO:0000269|PubMed:8475492,ECO:0000269|PubMed:9349616"
        feature = record.features[37]
        assert feature.type == "VARIANT"
        assert feature.location.start == 13
        assert feature.location.end == 14
        assert feature.qualifiers["note"] == "S -> L (in allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01 and allele A*74:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1782566,ECO:0000269|PubMed:2431040, ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:8795145, ECO:0000269|Ref.29"
        feature = record.features[38]
        assert feature.type == "VARIANT"
        assert feature.location.start == 22
        assert feature.location.end == 23
        assert feature.qualifiers["note"] == "W -> R (in allele A*74:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115"
        feature = record.features[39]
        assert feature.type == "VARIANT"
        assert feature.location.start == 32
        assert feature.location.end == 33
        assert feature.qualifiers["note"] == "F -> S (in allele A*23:01, allele A*24:02 and allele A*30:01; dbSNP:rs2075684)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:7871528,ECO:0000269|PubMed:9349616"
        feature = record.features[40]
        assert feature.type == "VARIANT"
        assert feature.location.start == 32
        assert feature.location.end == 33
        assert feature.qualifiers["note"] == "F -> T (in allele A*29:02, allele A*31:01 and allele A*33:01; requires 2 nucleotide substitutions)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1782566, ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:8795145"
        feature = record.features[41]
        assert feature.type == "VARIANT"
        assert feature.location.start == 32
        assert feature.location.end == 33
        assert feature.qualifiers["note"] == "F -> Y (in allele A*02:05, allele A*11:01, allele A*25:01, allele A*26:01, allele A*34:01, allele A*43:01, allele A*66:01, allele A*68:01 and allele A*69:01; dbSNP:rs2075684)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2437024,ECO:0000269|PubMed:2460344, ECO:0000269|PubMed:3496393,ECO:0000269|PubMed:3877632, ECO:0000269|PubMed:8016845,ECO:0000269|PubMed:8026990, ECO:0000269|PubMed:8475492"
        feature = record.features[42]
        assert feature.type == "VARIANT"
        assert feature.location.start == 40
        assert feature.location.end == 41
        assert feature.qualifiers["note"] == "R -> S (in allele A*30:01; dbSNP:rs1059423)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:7871528"
        feature = record.features[43]
        assert feature.type == "VARIANT"
        assert feature.location.start == 54
        assert feature.location.end == 55
        assert feature.qualifiers["note"] == "T -> S (in allele A*80:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791"
        feature = record.features[44]
        assert feature.type == "VARIANT"
        assert feature.location.start == 58
        assert feature.location.end == 59
        assert feature.qualifiers["note"] == "R -> Q (in allele A*80:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791"
        feature = record.features[45]
        assert feature.type == "VARIANT"
        assert feature.location.start == 66
        assert feature.location.end == 67
        assert feature.qualifiers["note"] == "Q -> R (in allele A*02:05; dbSNP:rs41559117)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:3496393"
        feature = record.features[46]
        assert feature.type == "VARIANT"
        assert feature.location.start == 67
        assert feature.location.end == 68
        assert feature.qualifiers["note"] == "R -> K (in alleles A*01:01 and allele A*36:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2251137, ECO:0000269|PubMed:2715640,ECO:0000269|PubMed:9349617"
        feature = record.features[47]
        assert feature.type == "VARIANT"
        assert feature.location.start == 79
        assert feature.location.end == 80
        assert feature.qualifiers["note"] == "G -> E (in allele A*80:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791"
        feature = record.features[48]
        assert feature.type == "VARIANT"
        assert feature.location.start == 79
        assert feature.location.end == 80
        assert feature.qualifiers["note"] == "G -> R (in allele A*30:01 and allele A*31:01; dbSNP:rs1059449)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:7871528,ECO:0000269|PubMed:8795145"
        feature = record.features[49]
        assert feature.type == "VARIANT"
        assert feature.location.start == 85
        assert feature.location.end == 86
        assert feature.qualifiers["note"] == "Q -> E (in allele A*23:01, allele 24:02 and allele A*80:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:8188325, ECO:0000269|PubMed:8284791,ECO:0000269|PubMed:9349616"
        feature = record.features[50]
        assert feature.type == "VARIANT"
        assert feature.location.start == 85
        assert feature.location.end == 86
        assert feature.qualifiers["note"] == "Q -> G (in allele A*02:01 and allele A*02:05; requires 2 nucleotide substitutions)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2982951, ECO:0000269|PubMed:3496393,ECO:0000269|PubMed:7836067"
        feature = record.features[51]
        assert feature.type == "VARIANT"
        assert feature.location.start == 85
        assert feature.location.end == 86
        assert feature.qualifiers["note"] == "Q -> L (in alleles A*29:02 and allele A*43:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:1782566"
        feature = record.features[52]
        assert feature.type == "VARIANT"
        assert feature.location.start == 85
        assert feature.location.end == 86
        assert feature.qualifiers["note"] == "Q -> R (in allele A*25:01, allele A*26:01, allele A*33:01, allele A*34:01, allele A*66:01, allele A*68:01 and allele A*69:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:3877632, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8475492"
        feature = record.features[53]
        assert feature.type == "VARIANT"
        assert feature.location.start == 86
        assert feature.location.end == 87
        assert feature.qualifiers["note"] == "E -> N (in alleles A*25:01, allele A*26:01, allele A*33:01, allele A*34:01, allele A*66:01, allele A*68:01 and allele A*69:01; requires 2 nucleotide substitutions)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:3877632, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8475492"
        feature = record.features[54]
        assert feature.type == "VARIANT"
        assert feature.location.start == 86
        assert feature.location.end == 87
        assert feature.qualifiers["note"] == "E -> Q (in allele A*29:02 and allele A*43:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:1782566"
        feature = record.features[55]
        assert feature.type == "VARIANT"
        assert feature.location.start == 88
        assert feature.location.end == 89
        assert feature.qualifiers["note"] == "R -> G (in allele A*23:01 and allele 24:02; dbSNP:rs199474430)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:9349616"
        feature = record.features[56]
        assert feature.type == "VARIANT"
        assert feature.location.start == 89
        assert feature.location.end == 90
        assert feature.qualifiers["note"] == "N -> K (in allele A*02:01, allele A*02:05, allele A*23:01, allele 24:02 and allele A*34:01; dbSNP:rs199474436)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:1729171, ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2982951, ECO:0000269|PubMed:3496393,ECO:0000269|PubMed:7836067, ECO:0000269|PubMed:9349616"
        feature = record.features[57]
        assert feature.type == "VARIANT"
        assert feature.location.start == 90
        assert feature.location.end == 91
        assert feature.qualifiers["note"] == "V -> M (in allele A*01:01 and allele A*36:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2251137, ECO:0000269|PubMed:2715640,ECO:0000269|PubMed:9349617"
        feature = record.features[58]
        assert feature.type == "VARIANT"
        assert feature.location.start == 93
        assert feature.location.end == 94
        assert feature.qualifiers["note"] == "Q -> H (in allele A*01:01, allele A*02:01, allele A*02:05, allele A*23:01, allele 24:02, allele A*25:01, allele A*26:01, allele A*31:01, allele A*32:01, allele A*33:01, allele A*36:01, allele A*43:01, allele A*74:01 and allele A*80:01; dbSNP:rs78306866)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:2251137, ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2431040, ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:2715640, ECO:0000269|PubMed:2982951,ECO:0000269|PubMed:3496393, ECO:0000269|PubMed:7836067,ECO:0000269|PubMed:8026990, ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791, ECO:0000269|PubMed:8475492,ECO:0000269|PubMed:8795145, ECO:0000269|PubMed:9349616,ECO:0000269|PubMed:9349617, ECO:0000269|Ref.29"
        feature = record.features[59]
        assert feature.type == "VARIANT"
        assert feature.location.start == 96
        assert feature.location.end == 97
        assert feature.qualifiers["note"] == "T -> I (in allele A*31:01 and allele A*33:01; dbSNP:rs199474457)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:8795145"
        feature = record.features[60]
        assert feature.type == "VARIANT"
        assert feature.location.start == 97
        assert feature.location.end == 98
        assert feature.qualifiers["note"] == "D -> H (in allele A*02:01 and allele A*02:05)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2982951, ECO:0000269|PubMed:3496393,ECO:0000269|PubMed:7836067"
        feature = record.features[61]
        assert feature.type == "VARIANT"
        assert feature.location.start == 97
        assert feature.location.end == 98
        assert feature.qualifiers["note"] == "D -> N (in allele A*80:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791"
        feature = record.features[62]
        assert feature.type == "VARIANT"
        assert feature.location.start == 99
        assert feature.location.end == 100
        assert feature.qualifiers["note"] == "V -> A (in allele A*01:01, allele A*26:01, allele A*29:02, allele A*36:01, allele A*43:01 and allele A*80:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:1782566, ECO:0000269|PubMed:2251137,ECO:0000269|PubMed:2715640, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8188325, ECO:0000269|PubMed:8284791,ECO:0000269|PubMed:8475492, ECO:0000269|PubMed:9349617"
        feature = record.features[63]
        assert feature.type == "VARIANT"
        assert feature.location.start == 99
        assert feature.location.end == 100
        assert feature.qualifiers["note"] == "V -> E (in allele A*23:01, allele A*24:02, allele A*25:01 and allele A*32:01; dbSNP:rs1071742)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2431040,ECO:0000269|PubMed:9349616, ECO:0000269|Ref.29"
        feature = record.features[64]
        assert feature.type == "VARIANT"
        assert feature.location.start == 100
        assert feature.location.end == 101
        assert feature.qualifiers["note"] == "D -> N (allele A*01:01, allele A*23:01, allele A*24:02, allele A*26:01, allele A*29:02, allele A*36:01, allele A*43:01 and allele A*80:01; dbSNP:rs1136688)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:1729171, ECO:0000269|PubMed:1782566,ECO:0000269|PubMed:2251137, ECO:0000269|PubMed:2715640,ECO:0000269|PubMed:8026990, ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791, ECO:0000269|PubMed:8475492,ECO:0000269|PubMed:9349616, ECO:0000269|PubMed:9349617"
        feature = record.features[65]
        assert feature.type == "VARIANT"
        assert feature.location.start == 100
        assert feature.location.end == 101
        assert feature.qualifiers["note"] == "D -> S (in allele A*25:01 and allele A*32:01; requires 2 nucleotide substitutions)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2431040, ECO:0000269|Ref.29"
        feature = record.features[66]
        assert feature.type == "VARIANT"
        assert feature.location.start == 102
        assert feature.location.end == 107
        assert feature.qualifiers["note"] == "GTLRG -> RIALR (in allele A*23:01, allele A*24:02, allele A*25:01 and allele A*32:01; Bw4 motif RIALR is involved in the recognition of NK cell inhibitory receptor KIR3DL1)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:17182537,ECO:0000269|PubMed:1729171, ECO:0000269|PubMed:18502829,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2431040,ECO:0000269|PubMed:9349616, ECO:0000269|Ref.29"
        feature = record.features[67]
        assert feature.type == "VARIANT"
        assert feature.location.start == 113
        assert feature.location.end == 114
        assert feature.qualifiers["note"] == "A -> D (in allele A*01:01, allele A*11:01, allele A*25:01, allele A*26:01, allele A*34:01, allele A*36:01, allele A*43:01, allele A*66:01 and allele A*80:01; dbSNP:rs1136692)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2251137, ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2437024, ECO:0000269|PubMed:2460344,ECO:0000269|PubMed:2715640, ECO:0000269|PubMed:8016845,ECO:0000269|PubMed:8026990, ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791, ECO:0000269|PubMed:8475492,ECO:0000269|PubMed:9349617"
        feature = record.features[68]
        assert feature.type == "VARIANT"
        assert feature.location.start == 118
        assert feature.location.end == 119
        assert feature.qualifiers["note"] == "I -> L (in allele A*02:05, allele A*23:01 and allele 24:02; dbSNP:rs1071743)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:3496393, ECO:0000269|PubMed:9349616"
        feature = record.features[69]
        assert feature.type == "VARIANT"
        assert feature.location.start == 118
        assert feature.location.end == 119
        assert feature.qualifiers["note"] == "I -> V (in allele A*02:01 and allele A*69:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2982951, ECO:0000269|PubMed:3877632,ECO:0000269|PubMed:7836067"
        feature = record.features[70]
        assert feature.type == "VARIANT"
        assert feature.location.start == 120
        assert feature.location.end == 121
        assert feature.qualifiers["note"] == "I -> M (in allele A*23:01, allele 24:02, allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01, allele A*68:01 and allele A*74:01; dbSNP:rs1136695)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:1782566, ECO:0000269|PubMed:2431040,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:3877632,ECO:0000269|PubMed:8795145, ECO:0000269|PubMed:9349616,ECO:0000269|Ref.29"
        feature = record.features[71]
        assert feature.type == "VARIANT"
        assert feature.location.start == 120
        assert feature.location.end == 121
        assert feature.qualifiers["note"] == "I -> R (in allele A*02:01, allele A*02:05, allele A*25:01, allele A*26:01, allele A*34:01, allele A*43:01, allele A*66:01 and allele A*69:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2982951,ECO:0000269|PubMed:3496393, ECO:0000269|PubMed:3877632,ECO:0000269|PubMed:7836067, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8475492"
        feature = record.features[72]
        assert feature.type == "VARIANT"
        assert feature.location.start == 122
        assert feature.location.end == 123
        assert feature.qualifiers["note"] == "Y -> F (in allele A*23:01, allele 24:02; dbSNP:rs1136697)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:9349616"
        feature = record.features[73]
        assert feature.type == "VARIANT"
        assert feature.location.start == 128
        assert feature.location.end == 129
        assert feature.qualifiers["note"] == "S -> P (in allele A*01:01, allele A*11:01, allele A*25:01, allele A*26:01, allele A*32:01, allele A*34:01, allele A*36:01, allele A*43:01, allele A*66:01 and allele A*74:01; dbSNP:rs1136700)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2251137, ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2431040, ECO:0000269|PubMed:2437024,ECO:0000269|PubMed:2460344, ECO:0000269|PubMed:2715640,ECO:0000269|PubMed:8016845, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8475492, ECO:0000269|PubMed:9349617,ECO:0000269|Ref.29"
        feature = record.features[74]
        assert feature.type == "VARIANT"
        assert feature.location.start == 130
        assert feature.location.end == 131
        assert feature.qualifiers["note"] == "G -> W (in allele A*02:01, allele A*02:05 and allele A*69:01; dbSNP:rs1136702)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2982951, ECO:0000269|PubMed:3496393,ECO:0000269|PubMed:3877632, ECO:0000269|PubMed:7836067"
        feature = record.features[75]
        assert feature.type == "VARIANT"
        assert feature.location.start == 132
        assert feature.location.end == 133
        assert feature.qualifiers["note"] == "F -> L (in allele A*32:01 and allele A*74:01; dbSNP:rs1059488)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2431040, ECO:0000269|Ref.29"
        feature = record.features[76]
        assert feature.type == "VARIANT"
        assert feature.location.start == 137
        assert feature.location.end == 138
        assert feature.qualifiers["note"] == "R -> E (in allele A*30:01; requires 2 nucleotide substitutions)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:7871528"
        feature = record.features[77]
        assert feature.type == "VARIANT"
        assert feature.location.start == 137
        assert feature.location.end == 138
        assert feature.qualifiers["note"] == "R -> H (in allele A*02:01, allele A*02:05, allele A*23:01, allele A*24:02, allele A*69:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2982951,ECO:0000269|PubMed:3496393, ECO:0000269|PubMed:3877632,ECO:0000269|PubMed:7836067, ECO:0000269|PubMed:9349616"
        feature = record.features[78]
        assert feature.type == "VARIANT"
        assert feature.location.start == 137
        assert feature.location.end == 138
        assert feature.qualifiers["note"] == "R -> Q (in allele A*25:01, allele A*26:01, allele A*31:01, allele A*32:01, allele A*33:01, allele A*34:01, allele A*43:01, allele A*66:01, allele A*74:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2431040, ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:8026990, ECO:0000269|PubMed:8475492,ECO:0000269|PubMed:8795145, ECO:0000269|Ref.29"
        feature = record.features[79]
        assert feature.type == "VARIANT"
        assert feature.location.start == 139
        assert feature.location.end == 140
        assert feature.qualifiers["note"] == "D -> H (in allele A*30:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:7871528"
        feature = record.features[80]
        assert feature.type == "VARIANT"
        assert feature.location.start == 139
        assert feature.location.end == 140
        assert feature.qualifiers["note"] == "D -> Y (in allele A*02:01, allele A*02:05, allele A*23:01, allele A*24:02 and allele A*69:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2982951,ECO:0000269|PubMed:3496393, ECO:0000269|PubMed:3877632,ECO:0000269|PubMed:7836067, ECO:0000269|PubMed:9349616"
        feature = record.features[81]
        assert feature.type == "VARIANT"
        assert feature.location.start == 150
        assert feature.location.end == 151
        assert feature.qualifiers["note"] == "N -> K (in allele A*02:01, allele A*02:05, allele A*23:01, allele A*24:02, allele A*68:01 and allele A*69:01; dbSNP:rs1059509)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2982951,ECO:0000269|PubMed:3496393, ECO:0000269|PubMed:3877632,ECO:0000269|PubMed:7836067, ECO:0000269|PubMed:9349616"
        feature = record.features[82]
        assert feature.type == "VARIANT"
        assert feature.location.start == 165
        assert feature.location.end == 166
        assert feature.qualifiers["note"] == "I -> T (in allele A*02:01, allele A*02:05, allele A*68:01 and allele A*69:01; dbSNP:rs1059516)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2982951, ECO:0000269|PubMed:3496393,ECO:0000269|PubMed:3877632, ECO:0000269|PubMed:7836067"
        feature = record.features[83]
        assert feature.type == "VARIANT"
        assert feature.location.start == 167
        assert feature.location.end == 168
        assert feature.qualifiers["note"] == "K -> Q (in allele A*23:01, allele A*25:01, allele A*26:01, allele A*29:02, allele A*30:01, allele A*31:01, allele A*32:01, allele A*33:01, allele A*34:01, allele A*43:01, allele A*66:01 and allele A*74:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:1782566, ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2431040, ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:7871528, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8475492, ECO:0000269|PubMed:8795145,ECO:0000269|Ref.29"
        feature = record.features[84]
        assert feature.type == "VARIANT"
        assert feature.location.start == 168
        assert feature.location.end == 169
        assert feature.qualifiers["note"] == "R -> H (in allele A*02:01, allele A*02:05, allele A*68:01 and allele A*69:01; dbSNP:rs1059520)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2982951, ECO:0000269|PubMed:3496393,ECO:0000269|PubMed:3877632, ECO:0000269|PubMed:7836067"
        feature = record.features[85]
        assert feature.type == "VARIANT"
        assert feature.location.start == 172
        assert feature.location.end == 173
        assert feature.qualifiers["note"] == "A -> T (in allele A*25:01, allele A*26:01, allele A*34:01, allele A*43:01 and allele A*66:01; dbSNP:rs1059526)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8475492"
        feature = record.features[86]
        assert feature.type == "VARIANT"
        assert feature.location.start == 173
        assert feature.location.end == 174
        assert feature.qualifiers["note"] == "A -> V (in allele A*01:01 and allele A*36:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2251137, ECO:0000269|PubMed:2715640,ECO:0000269|PubMed:9349617"
        feature = record.features[87]
        assert feature.type == "VARIANT"
        assert feature.location.start == 174
        assert feature.location.end == 175
        assert feature.qualifiers["note"] == "H -> R (in allele A*23:01, allele A*29:02, allele A*30:01, allele A*31:01, allele A*32:01, allele A*33:01, allele A*74:01 and allele A*80:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:1782566, ECO:0000269|PubMed:2431040,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:7871528,ECO:0000269|PubMed:8188325, ECO:0000269|PubMed:8284791,ECO:0000269|PubMed:8795145, ECO:0000269|Ref.29"
        feature = record.features[88]
        assert feature.type == "VARIANT"
        assert feature.location.start == 175
        assert feature.location.end == 176
        assert feature.qualifiers["note"] == "E -> A (in allele A*01:01, allele A*11:01 and allele A*36:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2251137, ECO:0000269|PubMed:2437024,ECO:0000269|PubMed:2460344, ECO:0000269|PubMed:2715640,ECO:0000269|PubMed:8016845, ECO:0000269|PubMed:9349617"
        feature = record.features[89]
        assert feature.type == "VARIANT"
        assert feature.location.start == 175
        assert feature.location.end == 176
        assert feature.qualifiers["note"] == "E -> R (in allele A*80:01; requires 2 nucleotide substitutions)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791"
        feature = record.features[90]
        assert feature.type == "VARIANT"
        assert feature.location.start == 175
        assert feature.location.end == 176
        assert feature.qualifiers["note"] == "E -> V (in allele A*02:01, allele A*02:05, allele A*23:01, allele A*24:02, allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01, allele A*68:01, allele A*69:01 and allele A*74:01; results in inefficient T cell recognition of epitopes derived from influenza A virus.; dbSNP:rs9256983)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:1782566, ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2431040, ECO:0000269|PubMed:2456340,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:2982951,ECO:0000269|PubMed:3496393, ECO:0000269|PubMed:3877632,ECO:0000269|PubMed:7836067, ECO:0000269|PubMed:8795145,ECO:0000269|PubMed:9349616, ECO:0000269|Ref.29"
        feature = record.features[91]
        assert feature.type == "VARIANT"
        assert feature.location.start == 175
        assert feature.location.end == 176
        assert feature.qualifiers["note"] == "E -> W (in allele A*30:01; requires 2 nucleotide substitutions)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:7871528"
        feature = record.features[92]
        assert feature.type == "VARIANT"
        assert feature.location.start == 179
        assert feature.location.end == 180
        assert feature.qualifiers["note"] == "L -> Q (in allele A*11:01 and allele A*24:02)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:2437024, ECO:0000269|PubMed:2460344,ECO:0000269|PubMed:8016845, ECO:0000269|PubMed:9349616"
        feature = record.features[93]
        assert feature.type == "VARIANT"
        assert feature.location.start == 179
        assert feature.location.end == 180
        assert feature.qualifiers["note"] == "L -> R (in allele A*01:01 and allele A*36:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2251137, ECO:0000269|PubMed:2715640,ECO:0000269|PubMed:9349617"
        feature = record.features[94]
        assert feature.type == "VARIANT"
        assert feature.location.start == 179
        assert feature.location.end == 180
        assert feature.qualifiers["note"] == "L -> W (in allele A*02:05, allele A*25:01, allele A*26:01, allele A*34:01, allele A*43:01, allele A*66:01, allele A*68:01; dbSNP:rs9260156)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:3496393,ECO:0000269|PubMed:3877632, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8475492"
        feature = record.features[95]
        assert feature.type == "VARIANT"
        assert feature.location.start == 181
        assert feature.location.end == 182
        assert feature.qualifiers["note"] == "A -> V (in allele A*01:01 and allele A*36:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2251137, ECO:0000269|PubMed:2715640,ECO:0000269|PubMed:9349617"
        feature = record.features[96]
        assert feature.type == "VARIANT"
        assert feature.location.start == 184
        assert feature.location.end == 185
        assert feature.qualifiers["note"] == "D -> E (in allele A*01:01, allele A*02:01, allele A*02:05, allele A*11:01, allele A*23:01, allele A*24:02, allele A*25:01, allele A*26:01, allele A*29:02, allele A*30:01, allele A*31:01, allele A*32:01, allele A*33:01, allele A*34:01, allele A*36:01, allele A*43:01, allele A*66:01, allele A*68:01, allele A*69:01, allele A*74:01 and allele A*80:01; dbSNP:rs1059542)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:1782566, ECO:0000269|PubMed:2251137,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2431040,ECO:0000269|PubMed:2437024, ECO:0000269|PubMed:2460344,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:2715640,ECO:0000269|PubMed:2982951, ECO:0000269|PubMed:3496393,ECO:0000269|PubMed:3877632, ECO:0000269|PubMed:7836067,ECO:0000269|PubMed:7871528, ECO:0000269|PubMed:8016845,ECO:0000269|PubMed:8026990, ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791, ECO:0000269|PubMed:8475492,ECO:0000269|PubMed:8795145, ECO:0000269|PubMed:9349616,ECO:0000269|PubMed:9349617, ECO:0000269|Ref.29"
        feature = record.features[97]
        assert feature.type == "VARIANT"
        assert feature.location.start == 186
        assert feature.location.end == 187
        assert feature.qualifiers["note"] == "T -> E (in allele A*80:01; requires 2 nucleotide substitutions)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791"
        feature = record.features[98]
        assert feature.type == "VARIANT"
        assert feature.location.start == 186
        assert feature.location.end == 187
        assert feature.qualifiers["note"] == "T -> R (in allele A*01:01, allele A*11:01, allele A*25:01, allele A*26:01, allele A*43:01 and allele A*66:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2251137, ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2437024, ECO:0000269|PubMed:2460344,ECO:0000269|PubMed:2715640, ECO:0000269|PubMed:8016845,ECO:0000269|PubMed:8026990, ECO:0000269|PubMed:8475492,ECO:0000269|PubMed:9349617"
        feature = record.features[99]
        assert feature.type == "VARIANT"
        assert feature.location.start == 189
        assert feature.location.end == 190
        assert feature.qualifiers["note"] == "E -> D (in allele A*01:01, allele A*23:01, allele A*24:02 and allele A*80:01; dbSNP:rs879577815)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:2251137, ECO:0000269|PubMed:2715640,ECO:0000269|PubMed:8188325, ECO:0000269|PubMed:8284791,ECO:0000269|PubMed:9349616, ECO:0000269|PubMed:9349617"
        feature = record.features[100]
        assert feature.type == "VARIANT"
        assert feature.location.start == 190
        assert feature.location.end == 191
        assert feature.qualifiers["note"] == "W -> G (in allele A*01:01, allele A*23:01, allele A*24:02 and allele A*80:01; dbSNP:rs3098019)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:2251137, ECO:0000269|PubMed:2715640,ECO:0000269|PubMed:8188325, ECO:0000269|PubMed:8284791,ECO:0000269|PubMed:9349616, ECO:0000269|PubMed:9349617"
        feature = record.features[101]
        assert feature.type == "VARIANT"
        assert feature.location.start == 194
        assert feature.location.end == 195
        assert feature.qualifiers["note"] == "Y -> H (in allele A*33:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:2478623"
        feature = record.features[102]
        assert feature.type == "VARIANT"
        assert feature.location.start == 207
        assert feature.location.end == 208
        assert feature.qualifiers["note"] == "P -> A (in allele A*02:01, allele A*02:05, allele A*25:01, allele A*26:01, allele A*29:02, allele A*32:01, allele A*34:01, allele A*43:01, allele A*66:01, allele A*68:01, allele A*69:01 and allele A*74:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:1782566, ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2431040, ECO:0000269|PubMed:2982951,ECO:0000269|PubMed:3496393, ECO:0000269|PubMed:3877632,ECO:0000269|PubMed:7836067, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8475492, ECO:0000269|Ref.29"
        feature = record.features[103]
        assert feature.type == "VARIANT"
        assert feature.location.start == 209
        assert feature.location.end == 210
        assert feature.qualifiers["note"] == "K -> R (in allele A*33:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:2478623"
        feature = record.features[104]
        assert feature.type == "VARIANT"
        assert feature.location.start == 216
        assert feature.location.end == 217
        assert feature.qualifiers["note"] == "P -> A (in allele A*02:01, allele A*02:05, allele A*25:01, allele A*26:01, allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01, allele A*34:01, allele A*43:01, allele A*66:01, allele A*68:01, allele A*69:01, allele A*74:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1782566,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2431040,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:2982951,ECO:0000269|PubMed:3496393, ECO:0000269|PubMed:3877632,ECO:0000269|PubMed:7836067, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8475492, ECO:0000269|PubMed:8795145,ECO:0000269|Ref.29"
        feature = record.features[105]
        assert feature.type == "VARIANT"
        assert feature.location.start == 217
        assert feature.location.end == 218
        assert feature.qualifiers["note"] == "I -> V (in allele A*02:01, allele A*02:05, allele A*25:01, allele A*26:01, allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01, allele A*34:01, allele A*43:01, allele A*66:01, allele A*68:01, allele A*69:01 and allele A*74:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1782566,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2431040,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:2982951,ECO:0000269|PubMed:3496393, ECO:0000269|PubMed:3877632,ECO:0000269|PubMed:7836067, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8475492, ECO:0000269|PubMed:8795145,ECO:0000269|Ref.29"
        feature = record.features[106]
        assert feature.type == "VARIANT"
        assert feature.location.start == 230
        assert feature.location.end == 231
        assert feature.qualifiers["note"] == "G -> S (in allele A*02:01, allele A*02:05, allele A*25:01, allele A*26:01, allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01, allele A*34:01, allele A*43:01, allele A*66:01, allele A*68:01, allele A*69:01, allele A*74:01 and allele A*80:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1782566,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2431040,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:2982951,ECO:0000269|PubMed:3496393, ECO:0000269|PubMed:3877632,ECO:0000269|PubMed:7836067, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8188325, ECO:0000269|PubMed:8284791,ECO:0000269|PubMed:8475492, ECO:0000269|PubMed:8795145,ECO:0000269|Ref.29"
        feature = record.features[107]
        assert feature.type == "VARIANT"
        assert feature.location.start == 268
        assert feature.location.end == 269
        assert feature.qualifiers["note"] == "A -> V (in allele A*68:01; impairs binding to CD8A and reduces recognition by antigen-specific CD8-positive T cells)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:2784196,ECO:0000269|PubMed:3877632"
        feature = record.features[108]
        assert feature.type == "VARIANT"
        assert feature.location.start == 269
        assert feature.location.end == 270
        assert feature.qualifiers["note"] == "A -> S (in allele A*25:01, allele A*26:01, allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01, allele A*34:01, allele A*43:01, allele A*66:01 and allele A*74:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:1782566, ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2431040, ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:8026990, ECO:0000269|PubMed:8475492,ECO:0000269|Ref.29"
        feature = record.features[109]
        assert feature.type == "VARIANT"
        assert feature.location.start == 276
        assert feature.location.end == 277
        assert feature.qualifiers["note"] == "E -> K (in allele A*80:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791"
        feature = record.features[110]
        assert feature.type == "VARIANT"
        assert feature.location.start == 276
        assert feature.location.end == 277
        assert feature.qualifiers["note"] == "E -> Q (in allele A*02:01, allele A*02:05, allele A*25:01, allele A*26:01, allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01, allele A*34:01, allele A*43:01, allele A*66:01, allele A*68:01, allele A*69:01 and allele A*74:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1782566,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2431040,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:2982951,ECO:0000269|PubMed:3496393, ECO:0000269|PubMed:3877632,ECO:0000269|PubMed:7836067, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8475492, ECO:0000269|PubMed:8795145,ECO:0000269|Ref.29"
        feature = record.features[111]
        assert feature.type == "VARIANT"
        assert feature.location.start == 278
        assert feature.location.end == 279
        assert feature.qualifiers["note"] == "Q -> K (in allele A*80:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791"
        feature = record.features[112]
        assert feature.type == "VARIANT"
        assert feature.location.start == 291
        assert feature.location.end == 292
        assert feature.qualifiers["note"] == "K -> E (in allele A*80:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791"
        feature = record.features[113]
        assert feature.type == "VARIANT"
        assert feature.location.start == 299
        assert feature.location.end == 300
        assert feature.qualifiers["note"] == "L -> P (in allele A*02:01, allele A*02:05, allele A*23:01, allele A*24:02, allele A*25:01, allele A*26:01, allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01, allele A*34:01, allele A*43:01, allele A*66:01, allele A*68:01, allele A*69:01, allele A*74:01 and allele A*80:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:1782566, ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2431040, ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:2982951, ECO:0000269|PubMed:3496393,ECO:0000269|PubMed:3877632, ECO:0000269|PubMed:7836067,ECO:0000269|PubMed:8026990, ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791, ECO:0000269|PubMed:8475492,ECO:0000269|PubMed:8795145, ECO:0000269|PubMed:9349616,ECO:0000269|Ref.29"
        feature = record.features[114]
        assert feature.type == "VARIANT"
        assert feature.location.start == 305
        assert feature.location.end == 306
        assert feature.qualifiers["note"] == "I -> V (in allele A*23:01 and allele A*24:02)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:9349616"
        feature = record.features[115]
        assert feature.type == "VARIANT"
        assert feature.location.start == 306
        assert feature.location.end == 307
        assert feature.qualifiers["note"] == "P -> H (in allele A*23:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1729171"
        feature = record.features[116]
        assert feature.type == "VARIANT"
        assert feature.location.start == 311
        assert feature.location.end == 312
        assert feature.qualifiers["note"] == "I -> L (in allele A*34:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115"
        feature = record.features[117]
        assert feature.type == "VARIANT"
        assert feature.location.start == 317
        assert feature.location.end == 318
        assert feature.qualifiers["note"] == "L -> F (in allele A*02:01, allele A*02:05, allele A*25:01, allele A*26:01, allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01, allele A*34:01, allele A*43:01, allele A*66:01, allele A*68:01, allele A*69:01 and allele A*74:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1782566,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2431040,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:2982951,ECO:0000269|PubMed:3496393, ECO:0000269|PubMed:3877632,ECO:0000269|PubMed:7836067, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8475492, ECO:0000269|PubMed:8795145,ECO:0000269|Ref.29"
        feature = record.features[118]
        assert feature.type == "VARIANT"
        assert feature.location.start == 320
        assert feature.location.end == 321
        assert feature.qualifiers["note"] == "V -> M (in allele A*32:01 and allele A*74:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2431040, ECO:0000269|Ref.29"
        feature = record.features[119]
        assert feature.type == "VARIANT"
        assert feature.location.start == 321
        assert feature.location.end == 322
        assert feature.qualifiers["note"] == "I -> F (in allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01 and allele A*74:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1782566,ECO:0000269|PubMed:2431040, ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:8795145, ECO:0000269|Ref.29"
        feature = record.features[120]
        assert feature.type == "VARIANT"
        assert feature.location.start == 322
        assert feature.location.end == 323
        assert feature.qualifiers["note"] == "T -> A (in allele A*25:01, allele A*26:01, allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01, allele A*34:01, allele A*43:01, allele A*66:01, allele A*74:01 and allele A*80:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1782566,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2431040,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8188325, ECO:0000269|PubMed:8284791,ECO:0000269|PubMed:8475492, ECO:0000269|PubMed:8795145,ECO:0000269|Ref.29"
        feature = record.features[121]
        assert feature.type == "VARIANT"
        assert feature.location.start == 330
        assert feature.location.end == 331
        assert feature.qualifiers["note"] == "M -> R (in allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01 and allele A*74:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1782566,ECO:0000269|PubMed:2431040, ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:8795145, ECO:0000269|Ref.29"
        feature = record.features[122]
        assert feature.type == "VARIANT"
        assert feature.location.start == 333
        assert feature.location.end == 334
        assert feature.qualifiers["note"] == "R -> K (allele A*80:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791"
        feature = record.features[123]
        assert feature.type == "VARIANT"
        assert feature.location.start == 334
        assert feature.location.end == 335
        assert feature.qualifiers["note"] == "K -> N (in allele A*23:01 and allele A*24:02)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:9349616"
        feature = record.features[124]
        assert feature.type == "VARIANT"
        assert feature.location.start == 337
        assert feature.location.end == 338
        assert feature.qualifiers["note"] == "D -> V (allele A*80:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791"
        feature = record.features[125]
        assert feature.type == "VARIANT"
        assert feature.location.start == 344
        assert feature.location.end == 345
        assert feature.qualifiers["note"] == "T -> S (in allele A*02:01, allele A*02:05, allele A*23:01, allele A*24:02, allele A*25:01, allele A*26:01, allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01, allele A*34:01, allele A*43:01, allele A*66:01, allele A*68:01 allele A*69:01, allele A*74:01 and allele A*80:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:1782566, ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2431040, ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:2982951, ECO:0000269|PubMed:3496393,ECO:0000269|PubMed:3877632, ECO:0000269|PubMed:7836067,ECO:0000269|PubMed:8026990, ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791, ECO:0000269|PubMed:8475492,ECO:0000269|PubMed:8795145, ECO:0000269|PubMed:9349616,ECO:0000269|Ref.29"
        feature = record.features[126]
        assert feature.type == "VARIANT"
        assert feature.location.start == 357
        assert feature.location.end == 358
        assert feature.qualifiers["note"] == "V -> M (in allele A*25:01, allele A*26:01, allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01, allele A*34:01, allele A*43:01, allele A*66:01 and allele A*74:01)"
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1782566,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2431040,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8475492, ECO:0000269|PubMed:8795145,ECO:0000269|Ref.29"
        feature = record.features[127]
        assert feature.type == "MUTAGEN"
        assert feature.location.start == 109
        assert feature.location.end == 110
        assert feature.qualifiers["note"] == "N->Q: Impairs the recruitment of HLA-A*02 in the peptide-loading complex."
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:21263072"
        feature = record.features[128]
        assert feature.type == "MUTAGEN"
        assert feature.location.start == 155
        assert feature.location.end == 156
        assert feature.qualifiers["note"] == "S->C: Impairs the maturation of a peptide-receptive HLA-A*02-B2M complex."
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:8805302"
        feature = record.features[129]
        assert feature.type == "MUTAGEN"
        assert feature.location.start == 157
        assert feature.location.end == 158
        assert feature.qualifiers["note"] == "T->K: Impairs binding to TAP1-TAP2 transporter, resulting in impaired presentation of intracellular peptides."
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:8630735,ECO:0000269|PubMed:8805302"
        feature = record.features[130]
        assert feature.type == "STRAND"
        assert feature.location.start == 26
        assert feature.location.end == 36
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:3MRE"
        feature = record.features[131]
        assert feature.type == "STRAND"
        assert feature.location.start == 40
        assert feature.location.end == 43
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:3MRE"
        feature = record.features[132]
        assert feature.type == "STRAND"
        assert feature.location.start == 44
        assert feature.location.end == 52
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:3MRE"
        feature = record.features[133]
        assert feature.type == "STRAND"
        assert feature.location.start == 54
        assert feature.location.end == 61
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:3MRE"
        feature = record.features[134]
        assert feature.type == "STRAND"
        assert feature.location.start == 63
        assert feature.location.end == 66
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:3MRE"
        feature = record.features[135]
        assert feature.type == "STRAND"
        assert feature.location.start == 69
        assert feature.location.end == 73
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:4F7T"
        feature = record.features[136]
        assert feature.type == "HELIX"
        assert feature.location.start == 73
        assert feature.location.end == 78
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:3MRE"
        feature = record.features[137]
        assert feature.type == "HELIX"
        assert feature.location.start == 80
        assert feature.location.end == 108
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:3MRE"
        feature = record.features[138]
        assert feature.type == "STRAND"
        assert feature.location.start == 112
        assert feature.location.end == 115
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:3D25"
        feature = record.features[139]
        assert feature.type == "STRAND"
        assert feature.location.start == 117
        assert feature.location.end == 127
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:3MRE"
        feature = record.features[140]
        assert feature.type == "STRAND"
        assert feature.location.start == 130
        assert feature.location.end == 142
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:3MRE"
        feature = record.features[141]
        assert feature.type == "STRAND"
        assert feature.location.start == 144
        assert feature.location.end == 150
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:3MRE"
        feature = record.features[142]
        assert feature.type == "STRAND"
        assert feature.location.start == 151
        assert feature.location.end == 155
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:6EWA"
        feature = record.features[143]
        assert feature.type == "STRAND"
        assert feature.location.start == 156
        assert feature.location.end == 159
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:3MRE"
        feature = record.features[144]
        assert feature.type == "HELIX"
        assert feature.location.start == 161
        assert feature.location.end == 173
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:3MRE"
        feature = record.features[145]
        assert feature.type == "HELIX"
        assert feature.location.start == 175
        assert feature.location.end == 184
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:3MRE"
        feature = record.features[146]
        assert feature.type == "HELIX"
        assert feature.location.start == 186
        assert feature.location.end == 198
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:3MRE"
        feature = record.features[147]
        assert feature.type == "HELIX"
        assert feature.location.start == 199
        assert feature.location.end == 203
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:3MRE"
        feature = record.features[148]
        assert feature.type == "STRAND"
        assert feature.location.start == 209
        assert feature.location.end == 235
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:3MRE"
        feature = record.features[149]
        assert feature.type == "STRAND"
        assert feature.location.start == 237
        assert feature.location.end == 243
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:3MRE"
        feature = record.features[150]
        assert feature.type == "STRAND"
        assert feature.location.start == 245
        assert feature.location.end == 248
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:3MRE"
        feature = record.features[151]
        assert feature.type == "HELIX"
        assert feature.location.start == 248
        assert feature.location.end == 251
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:2GTW"
        feature = record.features[152]
        assert feature.type == "STRAND"
        assert feature.location.start == 251
        assert feature.location.end == 254
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:3MRE"
        feature = record.features[153]
        assert feature.type == "STRAND"
        assert feature.location.start == 260
        assert feature.location.end == 263
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:3MRE"
        feature = record.features[154]
        assert feature.type == "STRAND"
        assert feature.location.start == 264
        assert feature.location.end == 274
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:3MRE"
        feature = record.features[155]
        assert feature.type == "STRAND"
        assert feature.location.start == 274
        assert feature.location.end == 277
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:4JFD"
        feature = record.features[156]
        assert feature.type == "HELIX"
        assert feature.location.start == 277
        assert feature.location.end == 280
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:3MRE"
        feature = record.features[157]
        assert feature.type == "STRAND"
        assert feature.location.start == 280
        assert feature.location.end == 286
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:3MRE"
        feature = record.features[158]
        assert feature.type == "STRAND"
        assert feature.location.start == 289
        assert feature.location.end == 292
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:2V2X"
        feature = record.features[159]
        assert feature.type == "STRAND"
        assert feature.location.start == 293
        assert feature.location.end == 297
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:3MRE"
        feature = record.features[160]
        assert feature.type == "STRAND"
        assert feature.location.start == 347
        assert feature.location.end == 350
        assert feature.qualifiers["evidence"] == "ECO:0007829|PDB:4EN2"
        assert feature.id is None

        assert len(record.references) == 103
        reference = record.references[0]
        assert reference.authors == "Wan A.M., Ennis P., Parham P., Holmes N."
        assert reference.title == "The primary structure of HLA-A32 suggests a region involved in formation of the Bw4/Bw6 epitopes."
        assert len(reference.references) == 1
        assert reference.references[0] == ("PubMed", "2431040")
        reference = record.references[1]
        assert reference.authors == "Holmes N., Ennis P., Wan A.M., Denney D.W., Parham P."
        assert reference.title == "Multiple genetic mechanisms have contributed to the generation of the HLA-A2/A28 family of class I MHC molecules."
        assert len(reference.references) == 1
        assert reference.references[0] == ("PubMed", "3496393")
        reference = record.references[2]
        assert reference.authors == "Mayer W.E., Jonker M., Klein D., Ivanyi P., van Seventer G., Klein J."
        assert reference.title == "Nucleotide sequences of chimpanzee MHC class I alleles: evidence for trans-species mode of evolution."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "2460344")
        assert reference.references[1] == ("DOI", "10.1002/j.1460-2075.1988.tb03131.x")
        reference = record.references[3]
        assert reference.authors == "Trapani J.A., Mizuno S., Kang S.H., Yang S.Y., Dupont B."
        assert reference.title == "Molecular mapping of a new public HLA class I epitope shared by all HLA-B and HLA-C antigens and defined by a monoclonal antibody."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "2461903")
        assert reference.references[1] == ("DOI", "10.1007/bf02341610")
        reference = record.references[4]
        assert reference.authors == "Kato K., Trapani J.A., Allopenna J., Dupont B., Yang S.Y."
        assert reference.title == "Molecular analysis of the serologically defined HLA-Aw19 antigens. A genetically distinct family of HLA-A antigens comprising A29, A31, A32, and Aw33, but probably not A30."
        assert len(reference.references) == 1
        assert reference.references[0] == ("PubMed", "2478623")
        reference = record.references[5]
        assert reference.authors == "Parham P., Lawlor D.A., Lomen C.E., Ennis P.D."
        assert reference.title == "Diversity and diversification of HLA-A,B,C alleles."
        assert len(reference.references) == 1
        assert reference.references[0] == ("PubMed", "2715640")
        reference = record.references[6]
        assert reference.authors == "Ennis P.D., Zemmour J., Salter R.D., Parham P."
        assert reference.title == "Rapid cloning of HLA-A,B cDNA by using the polymerase chain reaction: frequency and nature of errors produced in amplification."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "2320591")
        assert reference.references[1] == ("DOI", "10.1073/pnas.87.7.2833")
        reference = record.references[7]
        assert reference.authors == "Tabary T., Prochnicka-Chalufour A., Cornillet P., Lehoang P., Betuel H., Cohen H.M."
        assert reference.title == "HLA-A29 sub-types and 'Birdshot' choroido-retinopathy susceptibility: a possible 'resistance motif' in the HLA-A29.1 molecule."
        assert len(reference.references) == 1
        assert reference.references[0] == ("PubMed", "1782566")
        reference = record.references[8]
        assert reference.authors == "Little A.-M., Madrigal J.A., Parham P."
        assert reference.title == "Molecular definition of an elusive third HLA-A9 molecule: HLA-A9.3."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "1729171")
        assert reference.references[1] == ("DOI", "10.1007/bf00216625")
        reference = record.references[9]
        assert reference.authors == "Madrigal J.A., Belich M.P., Hildebrand W.H., Benjamin R.J., Little A.-M., Zemmour J., Ennis P.D., Ward F.E., Petzl-Erler M.L., Martell R.W., du Toit E.D., Parham P."
        assert reference.title == "Distinctive HLA-A,B antigens of black populations formed by interallelic conversion."
        assert len(reference.references) == 1
        assert reference.references[0] == ("PubMed", "1431115")
        reference = record.references[10]
        assert reference.authors == "Belich M.P., Madrigal J.A., Hildebrand W.H., Zemmour J., Williams R.C., Luz R., Petzl-Erler M.L., Parham P."
        assert reference.title == "Unusual HLA-B alleles in two tribes of Brazilian Indians."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "1317015")
        assert reference.references[1] == ("DOI", "10.1038/357326a0")
        reference = record.references[11]
        assert reference.authors == "Madrigal J.A., Hildebrand W.H., Belich M.P., Benjamin R.J., Little A.-M., Zemmour J., Ennis P.D., Ward F.E., Petzl-Erler M.L., du Toit E.D., Parham P."
        assert reference.title == "Structural diversity in the HLA-A10 family of alleles: correlations with serology."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "8475492")
        assert reference.references[1] == ("DOI", "10.1111/j.1399-0039.1993.tb01982.x")
        reference = record.references[12]
        assert reference.authors == "Domena J.D., Hildebrand W.H., Bias W.B., Parham P."
        assert reference.title == "A sixth family of HLA-A alleles defined by HLA-A*8001."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "8284791")
        assert reference.references[1] == ("DOI", "10.1111/j.1399-0039.1993.tb02186.x")
        reference = record.references[13]
        assert reference.authors == "Ishikawa Y., Tokunaga K., Lin L., Imanishi T., Saitou S., Kashiwase K., Akaza T., Tadokoro K., Juji T."
        assert reference.title == "Sequences of four splits of HLA-A10 group. Implications for serologic cross-reactivities and their evolution."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "8026990")
        assert reference.references[1] == ("DOI", "10.1016/0198-8859(94)90263-1")
        reference = record.references[14]
        assert reference.authors == "Balas A., Garcia-Sanchez F., Gomez-Reino F., Vicario J.L."
        assert reference.title == "Characterization of a new and highly distinguishable HLA-A allele in a Spanish family."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "8188325")
        assert reference.references[1] == ("DOI", "10.1007/bf00176169")
        reference = record.references[15]
        assert reference.authors == "Lin L., Tokunaga K., Ishikawa Y., Bannai M., Kashiwase K., Kuwata S., Akaza T., Tadokoro K., Shibata Y., Juji T."
        assert reference.title == "Sequence analysis of serological HLA-A11 split antigens, A11.1 and A11.2."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "8016845")
        assert reference.references[1] == ("DOI", "10.1111/j.1399-0039.1994.tb02304.x")
        reference = record.references[16]
        assert reference.authors == "Olerup O., Daniels T.J., Baxter-Lowe L."
        assert reference.title == "Correct sequence of the A*3001 allele obtained by PCR-SSP typing and automated nucleotide sequencing."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "7871528")
        assert reference.references[1] == ("DOI", "10.1111/j.1399-0039.1994.tb02393.x")
        reference = record.references[17]
        assert reference.authors == "Sun Y., Liu S., Luo Y., Liang F., Xi Y."
        assert reference.title == "Identification and frequency of a novel HLA-A allele, A*110104."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "17092262")
        assert reference.references[1] == ("DOI", "10.1111/j.1399-0039.2006.00687.x")
        reference = record.references[18]
        assert reference.authors == "Strachan T., Sodoyer R., Damotte M., Jordan B.R."
        assert reference.title == "Complete nucleotide sequence of a functional class I HLA gene, HLA-A3: implications for the evolution of HLA genes."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "6609814")
        assert reference.references[1] == ("DOI", "10.1002/j.1460-2075.1984.tb01901.x")
        reference = record.references[19]
        assert reference.authors == "Holmes N., Parham P."
        assert reference.title == "Exon shuffling in vivo can generate novel HLA class I molecules."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "3877632")
        assert reference.references[1] == ("DOI", "10.1002/j.1460-2075.1985.tb04013.x")
        reference = record.references[20]
        assert reference.authors == "Koller B.H., Orr H.T."
        assert reference.title == "Cloning and complete sequence of an HLA-A2 gene: analysis of two HLA-A alleles at the nucleotide level."
        assert len(reference.references) == 1
        assert reference.references[0] == ("PubMed", "2982951")
        reference = record.references[21]
        assert reference.authors == "Cowan E.P., Jelachich M.L., Biddison W.E., Coligan J.E."
        assert reference.title == "DNA sequence of HLA-A11: remarkable homology with HLA-A3 allows identification of residues involved in epitopes recognized by antibodies and T cells."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "2437024")
        assert reference.references[1] == ("DOI", "10.1007/bf00404694")
        reference = record.references[22]
        assert reference.authors == "Girdlestone J."
        assert reference.title == "Nucleotide sequence of an HLA-A1 gene."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "2251137")
        assert reference.references[1] == ("DOI", "10.1093/nar/18.22.6701")
        reference = record.references[23]
        assert reference.authors == "Balas A., Garcia-Sanchez F., Gomez-Reino F., Vicario J.L."
        assert reference.title == "HLA class I allele (HLA-A2) expression defect associated with a mutation in its enhancer B inverted CAT box in two families."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "7836067")
        assert reference.references[1] == ("DOI", "10.1016/0198-8859(94)90087-6")
        reference = record.references[24]
        assert reference.authors == "Arnett K.L., Adams E.J., Parham P."
        assert reference.title == "On the sequence of A*3101."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "8795145")
        assert reference.references[1] == ("DOI", "10.1111/j.1399-0039.1996.tb02580.x")
        reference = record.references[25]
        assert reference.authors == "Laforet M., Froelich N., Parissiadis A., Bausinger H., Pfeiffer B., Tongio M.M."
        assert reference.title == "An intronic mutation responsible for a low level of expression of an HLA-A*24 allele."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "9349616")
        assert reference.references[1] == ("DOI", "10.1111/j.1399-0039.1997.tb02884.x")
        reference = record.references[26]
        assert reference.authors == "Laforet M., Froelich N., Parissiadis A., Pfeiffer B., Schell A., Faller B., Woehl-Jaegle M.L., Cazenave J.-P., Tongio M.M."
        assert reference.title == "A nucleotide insertion in exon 4 is responsible for the absence of expression of an HLA-A*01 allele."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "9349617")
        assert reference.references[1] == ("DOI", "10.1111/j.1399-0039.1997.tb02885.x")
        reference = record.references[27]
        assert reference.authors == "Zhu F., He Y., Zhang W., He J., He J., Xu X., Yan L."
        assert reference.title == "Analysis of the complete genomic sequence of HLA-A alleles in the Chinese Han population."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "19735485")
        assert reference.references[1] == ("DOI", "10.1111/j.1744-313x.2009.00874.x")
        reference = record.references[28]
        assert reference.authors == "Domena J.D."
        assert reference.title == ""
        assert len(reference.references) == 0
        reference = record.references[29]
        assert reference.authors == "Hurley C.K."
        assert reference.title == ""
        assert len(reference.references) == 0
        reference = record.references[30]
        assert reference.authors == "Ellexson M.E., Hildebrand W.H."
        assert reference.title == ""
        assert len(reference.references) == 0
        reference = record.references[31]
        assert reference.authors == "Mayor N.P."
        assert reference.title == "Full length sequence of an HLA-A*0301 intron 2 variant."
        assert len(reference.references) == 0
        reference = record.references[32]
        assert reference.authors == "The MGC Project Team"
        assert reference.title == "The status, quality, and expansion of the NIH full-length cDNA project: the Mammalian Gene Collection (MGC)."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "15489334")
        assert reference.references[1] == ("DOI", "10.1101/gr.2596504")
        reference = record.references[33]
        assert reference.authors == "Mungall A.J., Palmer S.A., Sims S.K., Edwards C.A., Ashurst J.L., Wilming L., Jones M.C., Horton R., Hunt S.E., Scott C.E., Gilbert J.G.R., Clamp M.E., Bethel G., Milne S., Ainscough R., Almeida J.P., Ambrose K.D., Andrews T.D., Ashwell R.I.S., Babbage A.K., Bagguley C.L., Bailey J., Banerjee R., Barker D.J., Barlow K.F., Bates K., Beare D.M., Beasley H., Beasley O., Bird C.P., Blakey S.E., Bray-Allen S., Brook J., Brown A.J., Brown J.Y., Burford D.C., Burrill W., Burton J., Carder C., Carter N.P., Chapman J.C., Clark S.Y., Clark G., Clee C.M., Clegg S., Cobley V., Collier R.E., Collins J.E., Colman L.K., Corby N.R., Coville G.J., Culley K.M., Dhami P., Davies J., Dunn M., Earthrowl M.E., Ellington A.E., Evans K.A., Faulkner L., Francis M.D., Frankish A., Frankland J., French L., Garner P., Garnett J., Ghori M.J., Gilby L.M., Gillson C.J., Glithero R.J., Grafham D.V., Grant M., Gribble S., Griffiths C., Griffiths M.N.D., Hall R., Halls K.S., Hammond S., Harley J.L., Hart E.A., Heath P.D., Heathcott R., Holmes S.J., Howden P.J., Howe K.L., Howell G.R., Huckle E., Humphray S.J., Humphries M.D., Hunt A.R., Johnson C.M., Joy A.A., Kay M., Keenan S.J., Kimberley A.M., King A., Laird G.K., Langford C., Lawlor S., Leongamornlert D.A., Leversha M., Lloyd C.R., Lloyd D.M., Loveland J.E., Lovell J., Martin S., Mashreghi-Mohammadi M., Maslen G.L., Matthews L., McCann O.T., McLaren S.J., McLay K., McMurray A., Moore M.J.F., Mullikin J.C., Niblett D., Nickerson T., Novik K.L., Oliver K., Overton-Larty E.K., Parker A., Patel R., Pearce A.V., Peck A.I., Phillimore B.J.C.T., Phillips S., Plumb R.W., Porter K.M., Ramsey Y., Ranby S.A., Rice C.M., Ross M.T., Searle S.M., Sehra H.K., Sheridan E., Skuce C.D., Smith S., Smith M., Spraggon L., Squares S.L., Steward C.A., Sycamore N., Tamlyn-Hall G., Tester J., Theaker A.J., Thomas D.W., Thorpe A., Tracey A., Tromans A., Tubby B., Wall M., Wallis J.M., West A.P., White S.S., Whitehead S.L., Whittaker H., Wild A., Willey D.J., Wilmer T.E., Wood J.M., Wray P.W., Wyatt J.C., Young L., Younger R.M., Bentley D.R., Coulson A., Durbin R.M., Hubbard T., Sulston J.E., Dunham I., Rogers J., Beck S."
        assert reference.title == "The DNA sequence and analysis of human chromosome 6."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "14574404")
        assert reference.references[1] == ("DOI", "10.1038/nature02055")
        reference = record.references[34]
        assert reference.authors == "Orr H.T., Lopez de Castro J.A., Parham P., Ploegh H.L., Strominger J.L."
        assert reference.title == "Comparison of amino acid sequences of two human histocompatibility antigens, HLA-A2 and HLA-B7: location of putative alloantigenic sites."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "92029")
        assert reference.references[1] == ("DOI", "10.1073/pnas.76.9.4395")
        reference = record.references[35]
        assert reference.authors == "Lopez de Castro J.A., Strominger J.L., Strong D.M., Orr H.T."
        assert reference.title == "Structure of crossreactive human histocompatibility antigens HLA-A28 and HLA-A2: possible implications for the generation of HLA polymorphism."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "6179086")
        assert reference.references[1] == ("DOI", "10.1073/pnas.79.12.3813")
        reference = record.references[36]
        assert reference.authors == "Jelachich M.L., Cowan E.P., Turner R.V., Coligan J.E., Biddison W.E."
        assert reference.title == "Analysis of the molecular basis of HLA-A3 recognition by cytotoxic T cells using defined mutants of the HLA-A3 molecule."
        assert len(reference.references) == 1
        assert reference.references[0] == ("PubMed", "2456340")
        reference = record.references[37]
        assert reference.authors == "Salter R.D., Norment A.M., Chen B.P., Clayberger C., Krensky A.M., Littman D.R., Parham P."
        assert reference.title == "Polymorphism in the alpha 3 domain of HLA-A molecules affects binding to CD8."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "2784196")
        assert reference.references[1] == ("DOI", "10.1038/338345a0")
        reference = record.references[38]
        assert reference.authors == "Traversari C., van der Bruggen P., Luescher I.F., Lurquin C., Chomez P., Van Pel A., De Plaen E., Amar-Costesec A., Boon T."
        assert reference.title == "A nonapeptide encoded by human gene MAGE-1 is recognized on HLA-A1 by cytolytic T lymphocytes directed against tumor antigen MZ2-E."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "1402688")
        assert reference.references[1] == ("DOI", "10.1084/jem.176.5.1453")
        reference = record.references[39]
        assert reference.authors == "DiBrino M., Tsuchida T., Turner R.V., Parker K.C., Coligan J.E., Biddison W.E."
        assert reference.title == "HLA-A1 and HLA-A3 T cell epitopes derived from influenza virus proteins predicted from peptide binding motifs."
        assert len(reference.references) == 1
        assert reference.references[0] == ("PubMed", "7504010")
        reference = record.references[40]
        assert reference.authors == "DiBrino M., Parker K.C., Shiloach J., Knierman M., Lukszo J., Turner R.V., Biddison W.E., Coligan J.E."
        assert reference.title == "Endogenous peptides bound to HLA-A3 possess a specific combination of anchor residues that permit identification of potential antigenic peptides."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "7679507")
        assert reference.references[1] == ("DOI", "10.1073/pnas.90.4.1508")
        reference = record.references[41]
        assert reference.authors == "DiBrino M., Parker K.C., Shiloach J., Turner R.V., Tsuchida T., Garfield M., Biddison W.E., Coligan J.E."
        assert reference.title == "Endogenous peptides with distinct amino acid anchor residue motifs bind to HLA-A1 and HLA-B8."
        assert len(reference.references) == 1
        assert reference.references[0] == ("PubMed", "7506728")
        reference = record.references[42]
        assert reference.authors == "Lewis J.W., Neisig A., Neefjes J., Elliott T."
        assert reference.title == "Point mutations in the alpha 2 domain of HLA-A2.1 define a functionally relevant interaction with TAP."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "8805302")
        assert reference.references[1] == ("DOI", "10.1016/s0960-9822(02)00611-5")
        reference = record.references[43]
        assert reference.authors == "Peace-Brewer A.L., Tussey L.G., Matsui M., Li G., Quinn D.G., Frelinger J.A."
        assert reference.title == "A point mutation in HLA-A*0201 results in failure to bind the TAP complex and to present virus-derived peptides to CTL."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "8630735")
        assert reference.references[1] == ("DOI", "10.1016/s1074-7613(00)80416-1")
        reference = record.references[44]
        assert reference.authors == "Boisgerault F., Khalil I., Tieng V., Connan F., Tabary T., Cohen J.H., Choppin J., Charron D., Toubert A."
        assert reference.title == "Definition of the HLA-A29 peptide ligand motif allows prediction of potential T-cell epitopes from the retinal soluble antigen, a candidate autoantigen in birdshot retinopathy."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "8622959")
        assert reference.references[1] == ("DOI", "10.1073/pnas.93.8.3466")
        reference = record.references[45]
        assert reference.authors == "Ikeda H., Lethe B.G., Lehmann F., van Baren N., Baurain J.-F., de Smet C., Chambost H., Vitale M., Moretta A., Boon T., Coulie P.G."
        assert reference.title == "Characterization of an antigen that is recognized on a melanoma showing partial HLA loss by CTL expressing an NK inhibitory receptor."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "9047241")
        assert reference.references[1] == ("DOI", "10.1016/s1074-7613(00)80426-4")
        reference = record.references[46]
        assert reference.authors == "Kawakami Y., Robbins P.F., Wang X., Tupesis J.P., Parkhurst M.R., Kang X., Sakaguchi K., Appella E., Rosenberg S.A."
        assert reference.title == "Identification of new melanoma epitopes on melanosomal proteins recognized by tumor infiltrating T lymphocytes restricted by HLA-A1, -A2, and -A3 alleles."
        assert len(reference.references) == 1
        assert reference.references[0] == ("PubMed", "9862734")
        reference = record.references[47]
        assert reference.authors == "Fukada K., Chujoh Y., Tomiyama H., Miwa K., Kaneko Y., Oka S., Takiguchi M."
        assert reference.title == "HLA-A*1101-restricted cytotoxic T lymphocyte recognition of HIV-1 Pol protein."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "10449296")
        assert reference.references[1] == ("DOI", "10.1097/00002030-199907300-00021")
        reference = record.references[48]
        assert reference.authors == "Johnson J.M., Nicot C., Fullen J., Ciminale V., Casareto L., Mulloy J.C., Jacobson S., Franchini G."
        assert reference.title == "Free major histocompatibility complex class I heavy chain is preferentially targeted for degradation by human T-cell leukemia/lymphotropic virus type 1 p12(I) protein."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "11390610")
        assert reference.references[1] == ("DOI", "10.1128/jvi.75.13.6086-6094.2001")
        reference = record.references[49]
        assert reference.authors == "Hewitt E.W., Duncan L., Mufti D., Baker J., Stevenson P.G., Lehner P.J."
        assert reference.title == "Ubiquitylation of MHC class I by the K3 viral protein signals internalization and TSG101-dependent degradation."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "12006494")
        assert reference.references[1] == ("DOI", "10.1093/emboj/21.10.2418")
        reference = record.references[50]
        assert reference.authors == "Nagata Y., Ono S., Matsuo M., Gnjatic S., Valmori D., Ritter G., Garrett W., Old L.J., Mellman I."
        assert reference.title == "Differential presentation of a soluble exogenous tumor antigen, NY-ESO-1, by distinct human dendritic cell populations."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "12138174")
        assert reference.references[1] == ("DOI", "10.1073/pnas.112331099")
        reference = record.references[51]
        assert reference.authors == "Kuzushima K., Hayashi N., Kudoh A., Akatsuka Y., Tsujimura K., Morishima Y., Tsurumi T."
        assert reference.title == "Tetramer-assisted identification and characterization of epitopes recognized by HLA A*2402-restricted Epstein-Barr virus-specific CD8+ T cells."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "12393434")
        assert reference.references[1] == ("DOI", "10.1182/blood-2002-04-1240")
        reference = record.references[52]
        assert reference.authors == "Satoh M., Takamiya Y., Oka S., Tokunaga K., Takiguchi M."
        assert reference.title == "Identification and characterization of HIV-1-specific CD8+ T cell epitopes presented by HLA-A*2601."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "15893615")
        assert reference.references[1] == ("DOI", "10.1016/j.vaccine.2005.02.022")
        reference = record.references[53]
        assert reference.authors == "Asemissen A.M., Keilholz U., Tenzer S., Mueller M., Walter S., Stevanovic S., Schild H., Letsch A., Thiel E., Rammensee H.G., Scheibenbogen C."
        assert reference.title == "Identification of a highly immunogenic HLA-A*01-binding T cell epitope of WT1."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "17189421")
        assert reference.references[1] == ("DOI", "10.1158/1078-0432.ccr-06-1337")
        reference = record.references[54]
        assert reference.authors == "Thananchai H., Gillespie G., Martin M.P., Bashirova A., Yawata N., Yawata M., Easterbrook P., McVicar D.W., Maenaka K., Parham P., Carrington M., Dong T., Rowland-Jones S."
        assert reference.title == "Cutting Edge: Allele-specific and peptide-dependent interactions between KIR3DL1 and HLA-A and HLA-B."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "17182537")
        assert reference.references[1] == ("DOI", "10.4049/jimmunol.178.1.33")
        reference = record.references[55]
        assert reference.authors == "Robek M.D., Garcia M.L., Boyd B.S., Chisari F.V."
        assert reference.title == "Role of immunoproteasome catalytic subunits in the immune response to hepatitis B virus."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "17079320")
        assert reference.references[1] == ("DOI", "10.1128/jvi.01779-06")
        reference = record.references[56]
        assert reference.authors == "Stern M., Ruggeri L., Capanni M., Mancusi A., Velardi A."
        assert reference.title == "Human leukocyte antigens A23, A24, and A32 but not A25 are ligands for KIR3DL1."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "18502829")
        assert reference.references[1] == ("DOI", "10.1182/blood-2008-02-137521")
        reference = record.references[57]
        assert reference.authors == "Brennan R.M., Burrows S.R."
        assert reference.title == "A mechanism for the HLA-A*01-associated risk for EBV+ Hodgkin lymphoma and infectious mononucleosis."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "18779413")
        assert reference.references[1] == ("DOI", "10.1182/blood-2008-06-162883")
        reference = record.references[58]
        assert reference.authors == "Chen R., Jiang X., Sun D., Han G., Wang F., Ye M., Wang L., Zou H."
        assert reference.title == "Glycoproteomics analysis of human liver tissue by combination of multiple enzyme digestion and hydrazide chemistry."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "19159218")
        assert reference.references[1] == ("DOI", "10.1021/pr8008012")
        reference = record.references[59]
        assert reference.authors == "Hadrup S.R., Bakker A.H., Shu C.J., Andersen R.S., van Veluw J., Hombrink P., Castermans E., Thor Straten P., Blank C., Haanen J.B., Heemskerk M.H., Schumacher T.N."
        assert reference.title == "Parallel detection of antigen-specific T-cell responses by multidimensional encoding of MHC multimers."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "19543285")
        assert reference.references[1] == ("DOI", "10.1038/nmeth.1345")
        reference = record.references[60]
        assert reference.authors == "Parmentier N., Stroobant V., Colau D., de Diesbach P., Morel S., Chapiro J., van Endert P., Van den Eynde B.J."
        assert reference.title == "Production of an antigenic peptide by insulin-degrading enzyme."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "20364150")
        assert reference.references[1] == ("DOI", "10.1038/ni.1862")
        reference = record.references[61]
        assert reference.authors == "Marsh S.G., Albert E.D., Bodmer W.F., Bontrop R.E., Dupont B., Erlich H.A., Fernandez-Vina M., Geraghty D.E., Holdsworth R., Hurley C.K., Lau M., Lee K.W., Mach B., Maiers M., Mayr W.R., Mueller C.R., Parham P., Petersdorf E.W., Sasazuki T., Strominger J.L., Svejgaard A., Terasaki P.I., Tiercy J.M., Trowsdale J."
        assert reference.title == "Nomenclature for factors of the HLA system, 2010."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "20356336")
        assert reference.references[1] == ("DOI", "10.1111/j.1399-0039.2010.01466.x")
        reference = record.references[62]
        assert reference.authors == "Rizvi S.M., Del Cid N., Lybarger L., Raghavan M."
        assert reference.title == "Distinct functions for the glycans of tapasin and heavy chains in the assembly of MHC class I molecules."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "21263072")
        assert reference.references[1] == ("DOI", "10.4049/jimmunol.1002959")
        reference = record.references[63]
        assert reference.authors == "Matthews P.C., Adland E., Listgarten J., Leslie A., Mkhwanazi N., Carlson J.M., Harndahl M., Stryhn A., Payne R.P., Ogwu A., Huang K.H., Frater J., Paioni P., Kloverpris H., Jooste P., Goedhals D., van Vuuren C., Steyn D., Riddell L., Chen F., Luzzi G., Balachandran T., Ndung'u T., Buus S., Carrington M., Shapiro R., Heckerman D., Goulder P.J."
        assert reference.title == "HLA-A*7401-mediated control of HIV viremia is independent of its linkage disequilibrium with HLA-B*5703."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "21498667")
        assert reference.references[1] == ("DOI", "10.4049/jimmunol.1003711")
        reference = record.references[64]
        assert reference.authors == "Zhou H., Di Palma S., Preisinger C., Peng M., Polat A.N., Heck A.J., Mohammed S."
        assert reference.title == "Toward a comprehensive characterization of a human cancer cell phosphoproteome."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "23186163")
        assert reference.references[1] == ("DOI", "10.1021/pr300630k")
        reference = record.references[65]
        assert reference.authors == "Shimizu A., Kawana-Tachikawa A., Yamagata A., Han C., Zhu D., Sato Y., Nakamura H., Koibuchi T., Carlson J., Martin E., Brumme C.J., Shi Y., Gao G.F., Brumme Z.L., Fukai S., Iwamoto A."
        assert reference.title == "Structure of TCR and antigen complexes at an immunodominant CTL epitope in HIV-1 infection."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "24192765")
        assert reference.references[1] == ("DOI", "10.1038/srep03097")
        reference = record.references[66]
        assert reference.authors == "Bian Y., Song C., Cheng K., Dong M., Wang F., Huang J., Sun D., Wang L., Ye M., Zou H."
        assert reference.title == "An enzyme assisted RP-RPLC approach for in-depth analysis of human liver phosphoproteome."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "24275569")
        assert reference.references[1] == ("DOI", "10.1016/j.jprot.2013.11.014")
        reference = record.references[67]
        assert reference.authors == "Vaca Jacome A.S., Rabilloud T., Schaeffer-Reiss C., Rompais M., Ayoub D., Lane L., Bairoch A., Van Dorsselaer A., Carapito C."
        assert reference.title == "N-terminome analysis of the human mitochondrial proteome."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "25944712")
        assert reference.references[1] == ("DOI", "10.1002/pmic.201400617")
        reference = record.references[68]
        assert reference.authors == "Giam K., Ayala-Perez R., Illing P.T., Schittenhelm R.B., Croft N.P., Purcell A.W., Dudek N.L."
        assert reference.title == "A comprehensive analysis of peptides presented by HLA-A1."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "25880248")
        assert reference.references[1] == ("DOI", "10.1111/tan.12565")
        reference = record.references[69]
        assert reference.authors == "Morozov G.I., Zhao H., Mage M.G., Boyd L.F., Jiang J., Dolan M.A., Venna R., Norcross M.A., McMurtrey C.P., Hildebrand W., Schuck P., Natarajan K., Margulies D.H."
        assert reference.title == "Interaction of TAPBPR, a tapasin homolog, with MHC-I molecules promotes peptide editing."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "26869717")
        assert reference.references[1] == ("DOI", "10.1073/pnas.1519894113")
        reference = record.references[70]
        assert reference.authors == "Tripathi S.C., Peters H.L., Taguchi A., Katayama H., Wang H., Momin A., Jolly M.K., Celiktas M., Rodriguez-Canales J., Liu H., Behrens C., Wistuba I.I., Ben-Jacob E., Levine H., Molldrem J.J., Hanash S.M., Ostrin E.J."
        assert reference.title == "Immunoproteasome deficiency is a feature of non-small cell lung cancer with a mesenchymal phenotype and is associated with a poor outcome."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "26929325")
        assert reference.references[1] == ("DOI", "10.1073/pnas.1521812113")
        reference = record.references[71]
        assert reference.authors == "Ebstein F., Textoris-Taube K., Keller C., Golnik R., Vigneron N., Van den Eynde B.J., Schuler-Thurner B., Schadendorf D., Lorenz F.K., Uckert W., Urban S., Lehmann A., Albrecht-Koepke N., Janek K., Henklein P., Niewienda A., Kloetzel P.M., Mishto M."
        assert reference.title == "Proteasomes generate spliced epitopes by two different mechanisms and as efficiently as non-spliced epitopes."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "27049119")
        assert reference.references[1] == ("DOI", "10.1038/srep24032")
        reference = record.references[72]
        assert reference.authors == "Keib A., Mei Y.F., Cicin-Sain L., Busch D.H., Dennehy K.M."
        assert reference.title == "Measuring Antiviral Capacity of T Cell Responses to Adenovirus."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "30530481")
        assert reference.references[1] == ("DOI", "10.4049/jimmunol.1801003")
        reference = record.references[73]
        assert reference.authors == "Oxford Immunology Network Covid-19 Response T cell Consortium; ISARIC4C Investigators; Peng Y., Mentzer A.J., Liu G., Yao X., Yin Z., Dong D., Dejnirattisai W., Rostron T., Supasa P., Liu C., Lopez-Camacho C., Slon-Campos J., Zhao Y., Stuart D.I., Paesen G.C., Grimes J.M., Antson A.A., Bayfield O.W., Hawkins D.E.D.P., Ker D.S., Wang B., Turtle L., Subramaniam K., Thomson P., Zhang P., Dold C., Ratcliff J., Simmonds P., de Silva T., Sopp P., Wellington D., Rajapaksa U., Chen Y.L., Salio M., Napolitani G., Paes W., Borrow P., Kessler B.M., Fry J.W., Schwabe N.F., Semple M.G., Baillie J.K., Moore S.C., Openshaw P.J.M., Ansari M.A., Dunachie S., Barnes E., Frater J., Kerr G., Goulder P., Lockett T., Levin R., Zhang Y., Jing R., Ho L.P., Cornall R.J., Conlon C.P., Klenerman P., Screaton G.R., Mongkolsapaya J., McMichael A., Knight J.C., Ogg G., Dong T."
        assert reference.title == "Broad and strong memory CD4+ and CD8+ T cells induced by SARS-CoV-2 in UK convalescent individuals following COVID-19."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "32887977")
        assert reference.references[1] == ("DOI", "10.1038/s41590-020-0782-6")
        reference = record.references[74]
        assert reference.authors == "Guo H.-C., Jardetzky T.S., Garrett T.P.J., Lane W.S., Strominger J.L., Wiley D.C."
        assert reference.title == "Different length peptides bind to HLA-Aw68 similarly at their ends but bulge out in the middle."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "1448153")
        assert reference.references[1] == ("DOI", "10.1038/360364a0")
        reference = record.references[75]
        assert reference.authors == "Silver M.L., Guo H.-C., Strominger J.L., Wiley D.C."
        assert reference.title == "Atomic structure of a human MHC molecule presenting an influenza virus peptide."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "1448154")
        assert reference.references[1] == ("DOI", "10.1038/360367a0")
        reference = record.references[76]
        assert reference.authors == "Madden D.R., Garboczi D.N., Wiley D.C."
        assert reference.title == "The antigenic identity of peptide-MHC complexes: a comparison of the conformations of five viral peptides presented by HLA-A2."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "7694806")
        assert reference.references[1] == ("DOI", "10.1016/0092-8674(93)90490-h")
        reference = record.references[77]
        assert reference.authors == "Collins E.J., Garboczi D.N., Wiley D.C."
        assert reference.title == "Three-dimensional structure of a peptide extending from one end of a class I MHC binding site."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "7935798")
        assert reference.references[1] == ("DOI", "10.1038/371626a0")
        reference = record.references[78]
        assert reference.authors == "Garboczi D.N., Ghosh P., Utz U., Fan Q.R., Biddison W.E., Wiley D.C."
        assert reference.title == "Structure of the complex between human T-cell receptor, viral peptide and HLA-A2."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "8906788")
        assert reference.references[1] == ("DOI", "10.1038/384134a0")
        reference = record.references[79]
        assert reference.authors == "Gao G.F., Tormo J., Gerth U.C., Wyer J.R., McMichael A.J., Stuart D.I., Bell J.I., Jones E.Y., Jakobsen B.K."
        assert reference.title == "Crystal structure of the complex between human CD8alpha(alpha) and HLA-A2."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "9177355")
        assert reference.references[1] == ("DOI", "10.1038/42523")
        reference = record.references[80]
        assert reference.authors == "Hillig R.C., Coulie P.G., Stroobant V., Saenger W., Ziegler A., Hulsmeyer M."
        assert reference.title == "High-resolution structure of HLA-A*0201 in complex with a tumour-specific antigenic peptide encoded by the MAGE-A4 gene."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "11502003")
        assert reference.references[1] == ("DOI", "10.1006/jmbi.2001.4816")
        reference = record.references[81]
        assert reference.authors == "Stewart-Jones G.B.E., McMichael A.J., Bell J.I., Stuart D.I., Jones E.Y."
        assert reference.title == "A structural basis for immunodominant human T cell receptor recognition."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "12796775")
        assert reference.references[1] == ("DOI", "10.1038/ni942")
        reference = record.references[82]
        assert reference.authors == "Blicher T., Kastrup J.S., Buus S., Gajhede M."
        assert reference.title == "High-resolution structure of HLA-A*1101 in complex with SARS nucleocapsid peptide."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "16041067")
        assert reference.references[1] == ("DOI", "10.1107/s0907444905013090")
        reference = record.references[83]
        assert reference.authors == "Ishizuka J., Stewart-Jones G.B., van der Merwe A., Bell J.I., McMichael A.J., Jones E.Y."
        assert reference.title == "The structural dynamics and energetics of an immunodominant T cell receptor are programmed by its Vbeta domain."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "18275829")
        assert reference.references[1] == ("DOI", "10.1016/j.immuni.2007.12.018")
        reference = record.references[84]
        assert reference.authors == "Gras S., Saulquin X., Reiser J.B., Debeaupuis E., Echasserieau K., Kissenpfennig A., Legoux F., Chouquet A., Le Gorrec M., Machillot P., Neveu B., Thielens N., Malissen B., Bonneville M., Housset D."
        assert reference.title == "Structural bases for the affinity-driven selection of a public TCR against a dominant human cytomegalovirus epitope."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "19542454")
        assert reference.references[1] == ("DOI", "10.4049/jimmunol.0900556")
        reference = record.references[85]
        assert reference.authors == "Kumar P., Vahedi-Faridi A., Saenger W., Ziegler A., Uchanska-Ziegler B."
        assert reference.title == "Conformational changes within the HLA-A1:MAGE-A1 complex induced by binding of a recombinant antibody fragment with TCR-like specificity."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "19177349")
        assert reference.references[1] == ("DOI", "10.1002/pro.4")
        reference = record.references[86]
        assert reference.authors == "Liu J., Wu P., Gao F., Qi J., Kawana-Tachikawa A., Xie J., Vavricka C.J., Iwamoto A., Li T., Gao G.F."
        assert reference.title == "Novel immunodominant peptide presentation strategy: a featured HLA-A*2402-restricted cytotoxic T-lymphocyte epitope stabilized by intrachain hydrogen bonds from severe acute respiratory syndrome coronavirus nucleocapsid protein."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "20844028")
        assert reference.references[1] == ("DOI", "10.1128/jvi.01464-10")
        reference = record.references[87]
        assert reference.authors == "Borbulevych O.Y., Do P., Baker B.M."
        assert reference.title == "Structures of native and affinity-enhanced WT1 epitopes bound to HLA-A*0201: implications for WT1-based cancer therapeutics."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "20619457")
        assert reference.references[1] == ("DOI", "10.1016/j.molimm.2010.06.005")
        reference = record.references[88]
        assert reference.authors == "McMahon R.M., Friis L., Siebold C., Friese M.A., Fugger L., Jones E.Y."
        assert reference.title == "Structure of HLA-A*0301 in complex with a peptide of proteolipid protein: insights into the role of HLA-A alleles in susceptibility to multiple sclerosis."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "21543847")
        assert reference.references[1] == ("DOI", "10.1107/s0907444911007888")
        reference = record.references[89]
        assert reference.authors == "Zhang S., Liu J., Cheng H., Tan S., Qi J., Yan J., Gao G.F."
        assert reference.title == "Structural basis of cross-allele presentation by HLA-A*0301 and HLA-A*1101 revealed by two HIV-derived peptide complexes."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "21943705")
        assert reference.references[1] == ("DOI", "10.1016/j.molimm.2011.08.015")
        reference = record.references[90]
        assert reference.authors == "Bulek A.M., Cole D.K., Skowera A., Dolton G., Gras S., Madura F., Fuller A., Miles J.J., Gostick E., Price D.A., Drijfhout J.W., Knight R.R., Huang G.C., Lissin N., Molloy P.E., Wooldridge L., Jakobsen B.K., Rossjohn J., Peakman M., Rizkallah P.J., Sewell A.K."
        assert reference.title == "Structural basis for the killing of human beta cells by CD8(+) T cells in type 1 diabetes."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "22245737")
        assert reference.references[1] == ("DOI", "10.1038/ni.2206")
        reference = record.references[91]
        assert reference.authors == "Quinones-Parra S., Grant E., Loh L., Nguyen T.H., Campbell K.A., Tong S.Y., Miller A., Doherty P.C., Vijaykrishna D., Rossjohn J., Gras S., Kedzierska K."
        assert reference.title == "Preexisting CD8+ T-cell immunity to the H7N9 influenza A virus varies across ethnicities."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "24395804")
        assert reference.references[1] == ("DOI", "10.1073/pnas.1322229111")
        reference = record.references[92]
        assert reference.authors == "Raman M.C., Rizkallah P.J., Simmons R., Donnellan Z., Dukes J., Bossi G., Le Provost G.S., Todorov P., Baston E., Hickman E., Mahon T., Hassan N., Vuidepot A., Sami M., Cole D.K., Jakobsen B.K."
        assert reference.title == "Direct molecular mimicry enables off-target cardiovascular toxicity by an enhanced affinity TCR designed for cancer immunotherapy."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "26758806")
        assert reference.references[1] == ("DOI", "10.1038/srep18851")
        reference = record.references[93]
        assert reference.authors == "Song I., Gil A., Mishra R., Ghersi D., Selin L.K., Stern L.J."
        assert reference.title == "Broad TCR repertoire and diverse structural solutions for recognition of an immunodominant CD8+ T cell epitope."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "28250417")
        assert reference.references[1] == ("DOI", "10.1038/nsmb.3383")
        reference = record.references[94]
        assert reference.authors == "LeHoang P., Ozdemir N., Benhamou A., Tabary T., Edelson C., Betuel H., Semiglia R., Cohen J.H."
        assert reference.title == "HLA-A29.2 subtype associated with birdshot retinochoroidopathy."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "1728143")
        assert reference.references[1] == ("DOI", "10.1016/s0002-9394(14)75749-6")
        reference = record.references[95]
        assert reference.authors == "Fogdell-Hahn A., Ligers A., Groenning M., Hillert J., Olerup O."
        assert reference.title == "Multiple sclerosis: a modifying influence of HLA class I genes in an HLA class II associated autoimmune disease."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "10746785")
        assert reference.references[1] == ("DOI", "10.1034/j.1399-0039.2000.550205.x")
        reference = record.references[96]
        assert reference.authors == "Nakanishi K., Inoko H."
        assert reference.title == "Combination of HLA-A24, -DQA1*03, and -DR9 contributes to acute-onset and early complete beta-cell destruction in type 1 diabetes: longitudinal study of residual beta-cell function."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "16731854")
        assert reference.references[1] == ("DOI", "10.2337/db05-1049")
        reference = record.references[97]
        assert reference.authors == "Skowera A., Ellis R.J., Varela-Calvino R., Arif S., Huang G.C., Van-Krinks C., Zaremba A., Rackham C., Allen J.S., Tree T.I., Zhao M., Dayan C.M., Sewell A.K., Unger W.W., Unger W., Drijfhout J.W., Ossendorp F., Roep B.O., Peakman M."
        assert reference.title == "CTLs are targeted to kill beta cells in patients with type 1 diabetes through recognition of a glucose-regulated preproinsulin epitope."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "18802479")
        assert reference.references[1] == ("DOI", "10.1172/jci35449")
        reference = record.references[98]
        assert reference.authors == "Friese M.A., Jakobsen K.B., Friis L., Etzensperger R., Craner M.J., McMahon R.M., Jensen L.T., Huygelen V., Jones E.Y., Bell J.I., Fugger L."
        assert reference.title == "Opposing effects of HLA class I molecules in tuning autoreactive CD8+ T cells in multiple sclerosis."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "18953350")
        assert reference.references[1] == ("DOI", "10.1038/nm.1881")
        reference = record.references[99]
        assert reference.authors == "Kronenberg D., Knight R.R., Estorninho M., Ellis R.J., Kester M.G., de Ru A., Eichmann M., Huang G.C., Powrie J., Dayan C.M., Skowera A., van Veelen P.A., Peakman M."
        assert reference.title == "Circulating preproinsulin signal peptide-specific CD8 T cells restricted by the susceptibility molecule HLA-A24 are expanded at onset of type 1 diabetes and kill beta-cells."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "22522618")
        assert reference.references[1] == ("DOI", "10.2337/db11-1520")
        reference = record.references[100]
        assert reference.authors == "McCormack M., Alfirevic A., Bourgeois S., Farrell J.J., Kasperaviciute D., Carrington M., Sills G.J., Marson T., Jia X., de Bakker P.I., Chinthapalli K., Molokhia M., Johnson M.R., O'Connor G.D., Chaila E., Alhusaini S., Shianna K.V., Radtke R.A., Heinzen E.L., Walley N., Pandolfo M., Pichler W., Park B.K., Depondt C., Sisodiya S.M., Goldstein D.B., Deloukas P., Delanty N., Cavalleri G.L., Pirmohamed M."
        assert reference.title == "HLA-A*3101 and carbamazepine-induced hypersensitivity reactions in Europeans."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "21428769")
        assert reference.references[1] == ("DOI", "10.1056/nejmoa1013297")
        reference = record.references[101]
        assert reference.authors == "Nakamura J., Meguro A., Ishii G., Mihara T., Takeuchi M., Mizuki Y., Yuda K., Yamane T., Kawagoe T., Ota M., Mizuki N."
        assert reference.title == "The association analysis between HLA-A*26 and Behcet's disease."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "30872678")
        assert reference.references[1] == ("DOI", "10.1038/s41598-019-40824-y")
        reference = record.references[102]
        assert reference.authors == "Robinson J., Guethlein L.A., Cereb N., Yang S.Y., Norman P.J., Marsh S.G.E., Parham P."
        assert reference.title == "Distinguishing functional polymorphism from random variation in the sequences of >10,000 HLA-A, -B and -C alleles."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "28650991")
        assert reference.references[1] == ("DOI", "10.1371/journal.pgen.1006862")

        # Check the two parsers agree on the essentials
        assert seq_record.seq == record.sequence
        assert seq_record.description == record.description
        assert seq_record.name == record.entry_name
        assert seq_record.id in record.accessions

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        with open(datafile) as test_handle:
            records = list(SeqIO.parse(test_handle, "swiss"))

        assert len(records) == 1
        assert isinstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        assert records[0].seq == seq_record.seq
        assert records[0].description == seq_record.description
        assert records[0].name == seq_record.name
        assert records[0].id == seq_record.id

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        assert len(records) == 1
        assert isinstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        assert records[0].sequence == record.sequence
        assert records[0].description == record.description
        assert records[0].entry_name == record.entry_name
        assert records[0].accessions == record.accessions

    def test_O23729(self):
        """Parsing SwissProt file O23729.txt."""
        filename = "O23729.txt"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        assert isinstance(seq_record, SeqRecord)

        assert seq_record.id == "O23729"
        assert seq_record.name == "CHS3_BROFI"
        assert seq_record.description == "RecName: Full=Chalcone synthase 3; EC=2.3.1.74; AltName: Full=Naringenin-chalcone synthase 3;"
        assert repr(seq_record.seq) == "Seq('MAPAMEEIRQAQRAEGPAAVLAIGTSTPPNALYQADYPDYYFRITKSEHLTELK...GAE')"

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        assert record.entry_name == "CHS3_BROFI"
        assert record.accessions == ["O23729"]
        assert record.organism_classification == [
                "Eukaryota",
                "Viridiplantae",
                "Streptophyta",
                "Embryophyta",
                "Tracheophyta",
                "Spermatophyta",
                "Magnoliopsida",
                "Liliopsida",
                "Asparagales",
                "Orchidaceae",
                "Epidendroideae",
                "Vandeae",
                "Adrorhizinae",
                "Bromheadia",
            ]
        assert record.seqinfo == (394, 42942, "2F8D14AF4870BBB2")

        assert len(record.features) == 2
        feature = record.features[0]
        assert feature.type == "CHAIN"
        assert feature.location.start == 0
        assert feature.location.end == 394
        assert feature.qualifiers["note"] == "Chalcone synthase 3"
        assert feature.id == "PRO_0000215956"
        feature = record.features[1]
        assert feature.type == "ACT_SITE"
        assert feature.location.start == 164
        assert feature.location.end == 165
        assert feature.qualifiers["evidence"] == "ECO:0000255|PROSITE-ProRule:PRU10023"
        assert feature.id is None

        assert len(record.references) == 1
        assert record.references[0].authors == "Liew C.F., Lim S.H., Loh C.S., Goh C.J."
        assert record.references[0].title == "Molecular cloning and sequence analysis of chalcone synthase cDNAs of Bromheadia finlaysoniana."
        assert len(record.references[0].references) == 0

        # Check the two parsers agree on the essentials
        assert seq_record.seq == record.sequence
        assert seq_record.description == record.description
        assert seq_record.name == record.entry_name
        assert seq_record.id in record.accessions

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        with open(datafile) as test_handle:
            records = list(SeqIO.parse(test_handle, "swiss"))

        assert len(records) == 1
        assert isinstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        assert records[0].seq == seq_record.seq
        assert records[0].description == seq_record.description
        assert records[0].name == seq_record.name
        assert records[0].id == seq_record.id

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        assert len(records) == 1
        assert isinstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        assert records[0].sequence == record.sequence
        assert records[0].description == record.description
        assert records[0].entry_name == record.entry_name
        assert records[0].accessions == record.accessions

    def test_Q13639(self):
        """Parsing SwissProt file Q13639."""
        filename = "Q13639.txt"

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        assert isinstance(seq_record, SeqRecord)

        assert seq_record.id == "Q13639"
        assert seq_record.name == "5HT4R_HUMAN"
        assert seq_record.description == "RecName: Full=5-hydroxytryptamine receptor 4; Short=5-HT-4; Short=5-HT4; AltName: Full=Serotonin receptor 4;"
        assert repr(seq_record.seq) == "Seq('MDKLDANVSSEEGFGSVEKVVLLTFLSTVILMAILGNLLVMVAVCWDRQLRKIK...SDT')"

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        assert record.entry_name == "5HT4R_HUMAN"
        assert record.accessions == [
                "Q13639",
                "Q96KH9",
                "Q96KI0",
                "Q9H199",
                "Q9NY73",
                "Q9UBM6",
                "Q9UBT4",
                "Q9UE22",
                "Q9UE23",
                "Q9UQR6",
            ]
        assert record.organism_classification == [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Euteleostomi",
                "Mammalia",
                "Eutheria",
                "Euarchontoglires",
                "Primates",
                "Haplorrhini",
                "Catarrhini",
                "Hominidae",
                "Homo",
            ]
        assert record.seqinfo == (388, 43761, "7FCFEC60E7BDF560")

        assert len(record.features) == 26
        feature = record.features[0]
        assert feature.type == "CHAIN"
        assert feature.location.start == 0
        assert feature.location.end == 388
        assert feature.qualifiers["description"] == "5-hydroxytryptamine receptor 4."
        assert feature.id == "PRO_0000068965"
        feature = record.features[1]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 0
        assert feature.location.end == 19
        assert feature.qualifiers["description"] == "Extracellular (By similarity)."
        assert feature.id is None
        feature = record.features[2]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 19
        assert feature.location.end == 40
        assert feature.qualifiers["description"] == "1 (By similarity)."
        assert feature.id is None
        feature = record.features[3]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 40
        assert feature.location.end == 58
        assert feature.qualifiers["description"] == "Cytoplasmic (By similarity)."
        assert feature.id is None
        feature = record.features[4]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 58
        assert feature.location.end == 79
        assert feature.qualifiers["description"] == "2 (By similarity)."
        assert feature.id is None
        feature = record.features[5]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 79
        assert feature.location.end == 93
        assert feature.qualifiers["description"] == "Extracellular (By similarity)."
        assert feature.id is None
        feature = record.features[6]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 93
        assert feature.location.end == 116
        assert feature.qualifiers["description"] == "3 (By similarity)."
        assert feature.id is None
        feature = record.features[7]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 116
        assert feature.location.end == 137
        assert feature.qualifiers["description"] == "Cytoplasmic (By similarity)."
        assert feature.id is None
        feature = record.features[8]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 137
        assert feature.location.end == 158
        assert feature.qualifiers["description"] == "4 (By similarity)."
        assert feature.id is None
        feature = record.features[9]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 158
        assert feature.location.end == 192
        assert feature.qualifiers["description"] == "Extracellular (By similarity)."
        assert feature.id is None
        feature = record.features[10]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 192
        assert feature.location.end == 213
        assert feature.qualifiers["description"] == "5 (By similarity)."
        assert feature.id is None
        feature = record.features[11]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 213
        assert feature.location.end == 260
        assert feature.qualifiers["description"] == "Cytoplasmic (By similarity)."
        assert feature.id is None
        feature = record.features[12]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 260
        assert feature.location.end == 281
        assert feature.qualifiers["description"] == "6 (By similarity)."
        assert feature.id is None
        feature = record.features[13]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 281
        assert feature.location.end == 294
        assert feature.qualifiers["description"] == "Extracellular (By similarity)."
        assert feature.id is None
        feature = record.features[14]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 294
        assert feature.location.end == 315
        assert feature.qualifiers["description"] == "7 (By similarity)."
        assert feature.id is None
        feature = record.features[15]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 315
        assert feature.location.end == 388
        assert feature.qualifiers["description"] == "Cytoplasmic (By similarity)."
        assert feature.id is None
        feature = record.features[16]
        assert feature.type == "LIPID"
        assert feature.location.start == 328
        assert feature.location.end == 329
        assert feature.qualifiers["description"] == "S-palmitoyl cysteine (By similarity)."
        assert feature.id is None
        feature = record.features[17]
        assert feature.type == "CARBOHYD"
        assert feature.location.start == 6
        assert feature.location.end == 7
        assert feature.qualifiers["description"] == "N-linked (GlcNAc...) (Potential)."
        assert feature.id is None
        feature = record.features[18]
        assert feature.type == "DISULFID"
        assert feature.location.start == 92
        assert feature.location.end == 184
        assert feature.qualifiers["description"] == "By similarity."
        assert feature.id is None
        feature = record.features[19]
        assert feature.type == "VAR_SEQ"
        assert feature.location.start == 168
        assert feature.location.end == 169
        assert feature.qualifiers["description"] == "L -> LERSLNQGLGQDFHA (in isoform 5-HT4(F))."
        assert feature.id == "VSP_001845"
        feature = record.features[20]
        assert feature.type == "VAR_SEQ"
        assert feature.location.start == 358
        assert feature.location.end == 388
        assert feature.qualifiers["description"] == "RDAVECGGQWESQCHPPATSPLVAAQPSDT -> SSGTETDRRNFGIRKRRLTKPS (in isoform 5-HT4(D))."
        assert feature.id == "VSP_001847"
        feature = record.features[21]
        assert feature.type == "VAR_SEQ"
        assert feature.location.start == 358
        assert feature.location.end == 388
        assert feature.qualifiers["description"] == "RDAVECGGQWESQCHPPATSPLVAAQPSDT -> SGCSPVSSFLLLFCNRPVPV (in isoform 5-HT4(E))."
        assert feature.id == "VSP_001846"
        feature = record.features[22]
        assert feature.type == "VAR_SEQ"
        assert feature.location.start == 359
        assert feature.location.end == 388
        assert feature.qualifiers["description"] == "DAVECGGQWESQCHPPATSPLVAAQPSDT -> YTVLHRGHHQELEKLPIHNDPESLESCF (in isoform 5-HT4(A))."
        assert feature.id == "VSP_001849"
        feature = record.features[23]
        assert feature.type == "VAR_SEQ"
        assert feature.location.start == 359
        assert feature.location.end == 388
        assert feature.qualifiers["description"] == "DAVECGGQWESQCHPPATSPLVAAQPSDT -> F (in isoform 5-HT4(C))."
        assert feature.id == "VSP_001848"
        feature = record.features[24]
        assert feature.type == "VAR_SEQ"
        assert feature.location.start == 359
        assert feature.location.end == 388
        assert feature.qualifiers["description"] == "Missing (in isoform 5-HT4(G))."
        assert feature.id == "VSP_001850"
        feature = record.features[25]
        assert feature.type == "VARIANT"
        assert feature.location.start == 371
        assert feature.location.end == 372
        assert feature.qualifiers["description"] == "C -> Y (in dbSNP:rs34826744)."
        assert feature.id == "VAR_049364"
        assert len(record.references) == 8

        assert record.references[0].authors == "Blondel O., Gastineau M., Dahmoune Y., Langlois M., Fischmeister R."
        assert record.references[0].title == "Cloning, expression, and pharmacology of four human 5-hydroxytryptamine 4 receptor isoforms produced by alternative splicing in the carboxyl terminus."
        assert len(record.references[0].references) == 2
        assert record.references[0].references[0] == ("MEDLINE", "98264328")
        assert record.references[0].references[1] == ("PubMed", "9603189")
        assert record.references[1].authors == "Van den Wyngaert I., Gommeren W., Jurzak M., Verhasselt P., Gordon R., Leysen J., Luyten W., Bender E."
        assert record.references[1].title == "Cloning and expression of 5-HT4 receptor species and splice variants."
        assert len(record.references[1].references) == 0
        assert record.references[2].authors == "Claeysen S., Faye P., Sebben M., Lemaire S., Bockaert J., Dumuis A."
        assert record.references[2].title == "Cloning and expression of human 5-HT4S receptors. Effect of receptor density on their coupling to adenylyl cyclase."
        assert len(record.references[2].references) == 2
        assert record.references[2].references[0] == ("MEDLINE", "98012006")
        assert record.references[2].references[1] == ("PubMed", "9351641")
        assert record.references[3].authors == "Claeysen S., Sebben M., Becamel C., Bockaert J., Dumuis A."
        assert record.references[3].title == "Novel brain-specific 5-HT4 receptor splice variants show marked constitutive activity: role of the C-terminal intracellular domain."
        assert len(record.references[3].references) == 2
        assert record.references[3].references[0] == ("MEDLINE", "99238795")
        assert record.references[3].references[1] == ("PubMed", "10220570")
        assert record.references[4].authors == "Vilaro M.T., Domenech T., Palacios J.M., Mengod G."
        assert record.references[4].title == "Cloning and characterization of multiple human 5-HT4 receptor variants including a novel variant that lacks the alternatively spliced C-terminal exon."
        assert record.references[4].location == "Submitted (SEP-2000) to the EMBL/GenBank/DDBJ databases."
        assert len(record.references[4].comments) == 1
        assert record.references[4].comments[0] == ("TISSUE", "Hippocampus")
        assert len(record.references[4].positions) == 1
        assert record.references[4].positions[0] == "NUCLEOTIDE SEQUENCE [MRNA] (ISOFORMS 5-HT4(A); 5-HT4(E) AND 5-HT4(G))."
        assert len(record.references[4].references) == 0
        assert len(record.references[5].positions) == 1
        assert record.references[5].positions[0] == "NUCLEOTIDE SEQUENCE [LARGE SCALE MRNA] (ISOFORM 5-HT4(B))."
        assert len(record.references[5].references) == 2
        assert record.references[5].references[0] == ("PubMed", "15489334")
        assert record.references[5].references[1] == ("DOI", "10.1101/gr.2596504")
        assert record.references[5].authors == "The MGC Project Team"
        assert record.references[5].title == "The status, quality, and expansion of the NIH full-length cDNA project: the Mammalian Gene Collection (MGC)."
        assert record.references[5].location == "Genome Res. 14:2121-2127(2004)."
        assert record.references[6].authors == "Bender E., Pindon A., van Oers I., Zhang Y.B., Gommeren W., Verhasselt P., Jurzak M., Leysen J., Luyten W."
        assert record.references[6].title == "Structure of the human serotonin 5-HT4 receptor gene and cloning of a novel 5-HT4 splice variant."
        assert len(record.references[6].references) == 3
        assert record.references[6].references[0] == ("MEDLINE", "20110418")
        assert record.references[6].references[1] == ("PubMed", "10646498")
        assert record.references[6].references[2] == ("DOI", "10.1046/j.1471-4159.2000.740478.x")
        assert record.references[7].authors == "Ullmer C., Schmuck K., Kalkman H.O., Luebbert H."
        assert record.references[7].title == "Expression of serotonin receptor mRNAs in blood vessels."
        assert len(record.references[7].references) == 3
        assert record.references[7].references[0] == ("MEDLINE", "95385798")
        assert record.references[7].references[1] == ("PubMed", "7656980")
        assert record.references[7].references[2] == ("DOI", "10.1016/0014-5793(95)00828-W")

        # Check the two parsers agree on the essentials
        assert seq_record.seq == record.sequence
        assert seq_record.description == record.description
        assert seq_record.name == record.entry_name
        assert seq_record.id in record.accessions

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        records = list(SeqIO.parse(datafile, "swiss"))

        assert len(records) == 1
        assert isinstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        assert records[0].seq == seq_record.seq
        assert records[0].description == seq_record.description
        assert records[0].name == seq_record.name
        assert records[0].id == seq_record.id

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        assert len(records) == 1
        assert isinstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        assert records[0].sequence == record.sequence
        assert records[0].description == record.description
        assert records[0].entry_name == record.entry_name
        assert records[0].accessions == record.accessions

    def test_P16235(self):
        """Parsing SwissProt file P16235.txt."""
        filename = "P16235.txt"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)
        seq_record = SeqIO.read(datafile, "swiss")

        assert isinstance(seq_record, SeqRecord)

        assert seq_record.id == "P16235"
        assert seq_record.name == "LSHR_RAT"
        assert seq_record.description == "RecName: Full=Lutropin-choriogonadotropic hormone receptor; Short=LH/CG-R; AltName: Full=Luteinizing hormone receptor; Short=LSH-R; Flags: Precursor;"
        assert repr(seq_record.seq) == "Seq('MGRRVPALRQLLVLAVLLLKPSQLQSRELSGSRCPEPCDCAPDGALRCPGPRAG...LTH')"

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        assert record.entry_name == "LSHR_RAT"
        assert record.accessions == ["P16235", "P70646", "Q63807", "Q63808", "Q63809", "Q6LDI7"]
        assert record.organism_classification == [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Euteleostomi",
                "Mammalia",
                "Eutheria",
                "Euarchontoglires",
                "Glires",
                "Rodentia",
                "Myomorpha",
                "Muroidea",
                "Muridae",
                "Murinae",
                "Rattus",
            ]
        assert record.seqinfo == (700, 78036, "31807E73BAC94F1F")

        assert len(record.features) == 56
        feature = record.features[0]
        assert feature.type == "SIGNAL"
        assert feature.location.start == 0
        assert feature.location.end == 26
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:2601325,ECO:0000269|PubMed:2925659"
        feature = record.features[1]
        assert feature.type == "CHAIN"
        assert feature.location.start == 26
        assert feature.location.end == 700
        assert feature.qualifiers["note"] == "Lutropin-choriogonadotropic hormone receptor"
        feature = record.features[2]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 26
        assert feature.location.end == 362
        assert feature.qualifiers["note"] == "Extracellular"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[3]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 362
        assert feature.location.end == 390
        assert feature.qualifiers["note"] == "Helical; Name=1"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[4]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 390
        assert feature.location.end == 399
        assert feature.qualifiers["note"] == "Cytoplasmic"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[5]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 399
        assert feature.location.end == 422
        assert feature.qualifiers["note"] == "Helical; Name=2"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[6]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 422
        assert feature.location.end == 443
        assert feature.qualifiers["note"] == "Extracellular"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[7]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 443
        assert feature.location.end == 466
        assert feature.qualifiers["note"] == "Helical; Name=3"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[8]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 466
        assert feature.location.end == 486
        assert feature.qualifiers["note"] == "Cytoplasmic"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[9]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 486
        assert feature.location.end == 509
        assert feature.qualifiers["note"] == "Helical; Name=4"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[10]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 509
        assert feature.location.end == 529
        assert feature.qualifiers["note"] == "Extracellular"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[11]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 529
        assert feature.location.end == 551
        assert feature.qualifiers["note"] == "Helical; Name=5"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[12]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 551
        assert feature.location.end == 574
        assert feature.qualifiers["note"] == "Cytoplasmic"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[13]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 574
        assert feature.location.end == 598
        assert feature.qualifiers["note"] == "Helical; Name=6"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[14]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 598
        assert feature.location.end == 609
        assert feature.qualifiers["note"] == "Extracellular"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[15]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 609
        assert feature.location.end == 631
        assert feature.qualifiers["note"] == "Helical; Name=7"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[16]
        assert feature.type == "TOPO_DOM"
        assert feature.location.start == 631
        assert feature.location.end == 700
        assert feature.qualifiers["note"] == "Cytoplasmic"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[17]
        assert feature.type == "REPEAT"
        assert feature.location.start == 125
        assert feature.location.end == 150
        assert feature.qualifiers["note"] == "LRR 1"
        feature = record.features[18]
        assert feature.type == "REPEAT"
        assert feature.location.start == 151
        assert feature.location.end == 175
        assert feature.qualifiers["note"] == "LRR 2"
        feature = record.features[19]
        assert feature.type == "REPEAT"
        assert feature.location.start == 175
        assert feature.location.end == 200
        assert feature.qualifiers["note"] == "LRR 3"
        feature = record.features[20]
        assert feature.type == "REPEAT"
        assert feature.location.start == 201
        assert feature.location.end == 224
        assert feature.qualifiers["note"] == "LRR 4"
        feature = record.features[21]
        assert feature.type == "REPEAT"
        assert feature.location.start == 224
        assert feature.location.end == 248
        assert feature.qualifiers["note"] == "LRR 5"
        feature = record.features[22]
        assert feature.type == "REPEAT"
        assert feature.location.start == 249
        assert feature.location.end == 271
        assert feature.qualifiers["note"] == "LRR 6"
        feature = record.features[23]
        assert feature.type == "MOD_RES"
        assert feature.location.start == 334
        assert feature.location.end == 335
        assert feature.qualifiers["note"] == "Sulfotyrosine"
        assert feature.qualifiers["evidence"] == "ECO:0000250|UniProtKB:P22888"
        feature = record.features[24]
        assert feature.type == "LIPID"
        assert feature.location.start == 646
        assert feature.location.end == 647
        assert feature.qualifiers["note"] == "S-palmitoyl cysteine"
        assert feature.qualifiers["evidence"] == "ECO:0000305|PubMed:7776964"
        feature = record.features[25]
        assert feature.type == "LIPID"
        assert feature.location.start == 647
        assert feature.location.end == 648
        assert feature.qualifiers["note"] == "S-palmitoyl cysteine"
        assert feature.qualifiers["evidence"] == "ECO:0000305|PubMed:7776964"
        feature = record.features[26]
        assert feature.type == "CARBOHYD"
        assert feature.location.start == 102
        assert feature.location.end == 103
        assert feature.qualifiers["note"] == "N-linked (GlcNAc...) asparagine"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[27]
        assert feature.type == "CARBOHYD"
        assert feature.location.start == 177
        assert feature.location.end == 178
        assert feature.qualifiers["note"] == "N-linked (GlcNAc...) asparagine"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[28]
        assert feature.type == "CARBOHYD"
        assert feature.location.start == 198
        assert feature.location.end == 199
        assert feature.qualifiers["note"] == "N-linked (GlcNAc...) asparagine"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[29]
        assert feature.type == "CARBOHYD"
        assert feature.location.start == 294
        assert feature.location.end == 295
        assert feature.qualifiers["note"] == "N-linked (GlcNAc...) asparagine"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[30]
        assert feature.type == "CARBOHYD"
        assert feature.location.start == 302
        assert feature.location.end == 303
        assert feature.qualifiers["note"] == "N-linked (GlcNAc...) asparagine"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[31]
        assert feature.type == "CARBOHYD"
        assert feature.location.start == 316
        assert feature.location.end == 317
        assert feature.qualifiers["note"] == "N-linked (GlcNAc...) asparagine"
        assert feature.qualifiers["evidence"] == "ECO:0000255"
        feature = record.features[32]
        assert feature.type == "DISULFID"
        assert feature.location.start == 442
        assert feature.location.end == 518
        assert feature.qualifiers["evidence"] == "ECO:0000255|PROSITE-ProRule:PRU00521"
        feature = record.features[33]
        assert feature.type == "VAR_SEQ"
        assert feature.location.start == 82
        assert feature.location.end == 132
        assert feature.qualifiers["note"] == "Missing (in isoform 1950)"
        assert feature.qualifiers["evidence"] == "ECO:0000305"
        feature = record.features[34]
        assert feature.type == "VAR_SEQ"
        assert feature.location.start == 132
        assert feature.location.end == 157
        assert feature.qualifiers["note"] == "Missing (in isoform 1759)"
        assert feature.qualifiers["evidence"] == "ECO:0000305"
        feature = record.features[35]
        assert feature.type == "VAR_SEQ"
        assert feature.location.start == 183
        assert feature.location.end == 700
        assert feature.qualifiers["note"] == "Missing (in isoform C2)"
        assert feature.qualifiers["evidence"] == "ECO:0000305"
        feature = record.features[36]
        assert feature.type == "VAR_SEQ"
        assert feature.location.start == 231
        assert feature.location.end == 293
        assert feature.qualifiers["note"] == "Missing (in isoform EA2, isoform EB and isoform B1)"
        assert feature.qualifiers["evidence"] == "ECO:0000305"
        feature = record.features[37]
        assert feature.type == "VAR_SEQ"
        assert feature.location.start == 231
        assert feature.location.end == 251
        assert feature.qualifiers["note"] == "DISSTKLQALPSHGLESIQT -> PCRATGWSPFRRSSPCLPTH (in isoform 2075)"
        assert feature.qualifiers["evidence"] == "ECO:0000305"
        feature = record.features[38]
        assert feature.type == "VAR_SEQ"
        assert feature.location.start == 251
        assert feature.location.end == 700
        assert feature.qualifiers["note"] == "Missing (in isoform 2075)"
        assert feature.qualifiers["evidence"] == "ECO:0000305"
        feature = record.features[39]
        assert feature.type == "VAR_SEQ"
        assert feature.location.start == 293
        assert feature.location.end == 367
        assert feature.qualifiers["note"] == "QNFSFSIFENFSKQCESTVRKADNETLYSAIFEENELSGWDYDYGFCSPKTLQCAPEPDAFNPCEDIMGYAFLR -> IFHFPFLKTSPNNAKAQLEKQITRRFIPPSLRRMNSVAGIMIMASVHPRHSNVLQNQMLSTPVKILWAMPSLGS (in isoform B1 and isoform B3)"
        assert feature.qualifiers["evidence"] == "ECO:0000305"
        feature = record.features[40]
        assert feature.type == "VAR_SEQ"
        assert feature.location.start == 293
        assert feature.location.end == 294
        assert feature.qualifiers["note"] == "Q -> P (in isoform C1)"
        assert feature.qualifiers["evidence"] == "ECO:0000305"
        feature = record.features[41]
        assert feature.type == "VAR_SEQ"
        assert feature.location.start == 294
        assert feature.location.end == 700
        assert feature.qualifiers["note"] == "Missing (in isoform C1)"
        assert feature.qualifiers["evidence"] == "ECO:0000305"
        feature = record.features[42]
        assert feature.type == "VAR_SEQ"
        assert feature.location.start == 320
        assert feature.location.end == 342
        assert feature.qualifiers["note"] == "YSAIFEENELSGWDYDYGFCSP -> LHGALPAAHCLRGLPNKRPVL (in isoform 1834, isoform 1759 and isoform EB)"
        assert feature.qualifiers["evidence"] == "ECO:0000305"
        feature = record.features[43]
        assert feature.type == "VAR_SEQ"
        assert feature.location.start == 342
        assert feature.location.end == 700
        assert feature.qualifiers["note"] == "Missing (in isoform 1834, isoform 1759 and isoform EB)"
        assert feature.qualifiers["evidence"] == "ECO:0000305"
        feature = record.features[44]
        assert feature.type == "VAR_SEQ"
        assert feature.location.start == 367
        assert feature.location.end == 700
        assert feature.qualifiers["note"] == "Missing (in isoform B1 and isoform B3)"
        assert feature.qualifiers["evidence"] == "ECO:0000305"
        feature = record.features[45]
        assert feature.type == "VARIANT"
        assert feature.location.start == 81
        assert feature.location.end == 82
        assert feature.qualifiers["note"] == "I -> M (in isoform 1950)"
        feature = record.features[46]
        assert feature.type == "VARIANT"
        assert feature.location.start == 178
        assert feature.location.end == 179
        assert feature.qualifiers["note"] == "E -> G (in isoform 1759)"
        feature = record.features[47]
        assert feature.type == "VARIANT"
        assert feature.location.start == 232
        assert feature.location.end == 233
        assert feature.qualifiers["note"] == "I -> T (in isoform 1950)"
        feature = record.features[48]
        assert feature.type == "VARIANT"
        assert feature.location.start == 645
        assert feature.location.end == 646
        assert feature.qualifiers["note"] == "G -> S (in isoform 1950)"
        feature = record.features[49]
        assert feature.type == "MUTAGEN"
        assert feature.location.start == 408
        assert feature.location.end == 409
        assert feature.qualifiers["note"] == "D->N: Significant reduction of binding."
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1714448"
        feature = record.features[50]
        assert feature.type == "MUTAGEN"
        assert feature.location.start == 435
        assert feature.location.end == 436
        assert feature.qualifiers["note"] == "D->N: No change in binding or cAMP production."
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1714448"
        feature = record.features[51]
        assert feature.type == "MUTAGEN"
        assert feature.location.start == 454
        assert feature.location.end == 455
        assert feature.qualifiers["note"] == "E->Q: No change in binding or cAMP production."
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1714448"
        feature = record.features[52]
        assert feature.type == "MUTAGEN"
        assert feature.location.start == 581
        assert feature.location.end == 582
        assert feature.qualifiers["note"] == "D->N: No change in binding or cAMP production."
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:1714448"
        feature = record.features[53]
        assert feature.type == "MUTAGEN"
        assert feature.location.start == 646
        assert feature.location.end == 647
        assert feature.qualifiers["note"] == "C->A: Trapped intracellularly and does not appear to become mature; when associated with A-648."
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:7776964"
        feature = record.features[54]
        assert feature.type == "MUTAGEN"
        assert feature.location.start == 647
        assert feature.location.end == 648
        assert feature.qualifiers["note"] == "C->A: Trapped intracellularly and does not appear to become mature; when associated with A-647."
        assert feature.qualifiers["evidence"] == "ECO:0000269|PubMed:7776964"
        feature = record.features[55]
        assert feature.type == "CONFLICT"
        assert feature.location.start == 32
        assert feature.location.end == 33
        assert feature.qualifiers["note"] == "R -> L (in Ref. 9; AA sequence)"
        assert feature.qualifiers["evidence"] == "ECO:0000305"
        assert feature.id is None
        assert len(record.references) == 11
        reference = record.references[0]
        assert reference.authors == "McFarland K.C., Sprengel R., Phillips H.S., Koehler M., Rosemblit N., Nikolics K., Segaloff D.L., Seeburg P.H."
        assert reference.title == "Lutropin-choriogonadotropin receptor: an unusual member of the G protein-coupled receptor family."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "2502842")
        assert reference.references[1] == ("DOI", "10.1126/science.2502842")
        reference = record.references[1]
        assert reference.authors == "Aatsinki J.T., Pietila E.M., Lakkakorpi J.T., Rajaniemi H.J."
        assert reference.title == "Expression of the LH/CG receptor gene in rat ovarian tissue is regulated by an extensive alternative splicing of the primary transcript."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "1353463")
        assert reference.references[1] == ("DOI", "10.1016/0303-7207(92)90079-l")
        reference = record.references[2]
        assert reference.authors == "Koo Y.B., Slaughter R.G., Ji T.H."
        assert reference.title == "Structure of the luteinizing hormone receptor gene and multiple exons of the coding sequence."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "2019252")
        assert reference.references[1] == ("DOI", "10.1210/endo-128-5-2297")
        reference = record.references[3]
        assert reference.authors == "Bernard M.P., Myers R.V., Moyle W.R."
        assert reference.title == "Cloning of rat lutropin (LH) receptor analogs lacking the soybean lectin domain."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "1976554")
        assert reference.references[1] == ("DOI", "10.1016/0303-7207(90)90034-6")
        reference = record.references[4]
        assert reference.authors == "Segaloff D.L., Sprengel R., Nikolics K., Ascoli M."
        assert reference.title == "Structure of the lutropin/choriogonadotropin receptor."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "2281186")
        assert reference.references[1] == ("DOI", "10.1016/b978-0-12-571146-3.50014-6")
        reference = record.references[5]
        assert reference.authors == "Tsai-Morris C.H., Buczko E., Wang W., Xie X.-Z., Dufau M.L."
        assert reference.title == "Structural organization of the rat luteinizing hormone (LH) receptor gene."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "2040640")
        assert reference.references[1] == ("DOI", "10.1016/s0021-9258(18)99170-2")
        reference = record.references[6]
        assert reference.authors == "Tsai-Morris C.H., Buczko E., Wang W., Dufau M.L."
        assert reference.title == "Intronic nature of the rat luteinizing hormone receptor gene defines a soluble receptor subspecies with hormone binding activity."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "2174034")
        assert reference.references[1] == ("DOI", "10.1016/s0021-9258(17)45380-4")
        reference = record.references[7]
        assert reference.authors == "Dufau M.L., Minegishi T., Buczko E.S., Delgado C.J., Zhang R."
        assert reference.title == "Characterization and structure of ovarian and testicular LH/hCG receptors."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "2601325")
        assert reference.references[1] == ("DOI", "10.1016/0022-4731(89)90482-2")
        reference = record.references[8]
        assert reference.authors == "Roche P.C., Ryan R.J."
        assert reference.title == "Purification, characterization, and amino-terminal sequence of rat ovarian receptor for luteinizing hormone/human choriogonadotropin."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "2925659")
        assert reference.references[1] == ("DOI", "10.1016/s0021-9258(18)83790-5")
        reference = record.references[9]
        assert reference.authors == "Ji I., Ji T.H."
        assert reference.title == "Asp383 in the second transmembrane domain of the lutropin receptor is important for high affinity hormone binding and cAMP production."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "1714448")
        assert reference.references[1] == ("DOI", "10.1016/s0021-9258(18)98570-4")
        reference = record.references[10]
        assert reference.authors == "Zhu H., Wang H., Ascoli M."
        assert reference.title == "The lutropin/choriogonadotropin receptor is palmitoylated at intracellular cysteine residues."
        assert len(reference.references) == 2
        assert reference.references[0] == ("PubMed", "7776964")
        assert reference.references[1] == ("DOI", "10.1210/mend.9.2.7776964")

        # Check the two parsers agree on the essentials
        assert seq_record.seq == record.sequence
        assert seq_record.description == record.description
        assert seq_record.name == record.entry_name
        assert seq_record.id in record.accessions

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        records = list(SeqIO.parse(datafile, "swiss"))

        assert len(records) == 1
        assert isinstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        assert records[0].seq == seq_record.seq
        assert records[0].description == seq_record.description
        assert records[0].name == seq_record.name
        assert records[0].id == seq_record.id

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        assert len(records) == 1
        assert isinstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        assert records[0].sequence == record.sequence
        assert records[0].description == record.description
        assert records[0].entry_name == record.entry_name
        assert records[0].accessions == record.accessions

    def test_sp012(self):
        """Parsing SwissProt file sp012."""
        filename = "sp012"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)
        seq_record = SeqIO.read(datafile, "swiss")

        assert isinstance(seq_record, SeqRecord)

        assert seq_record.id == "Q9Y736"
        assert seq_record.name == "Q9Y736"
        assert seq_record.description == "UBIQUITIN."
        assert repr(seq_record.seq) == "Seq('MQIFVKTLTGKTITLEVESSDTIDNVKTKIQDKEGIPPDQQRLIFAGKQLEDGR...GGN')"

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        assert record.entry_name == "Q9Y736"
        assert record.accessions == ["Q9Y736"]
        assert record.organism_classification == [
                "Eukaryota",
                "Fungi",
                "Ascomycota",
                "Pezizomycotina",
                "Eurotiomycetes",
                "Onygenales",
                "Arthrodermataceae",
                "mitosporic Arthrodermataceae",
                "Trichophyton",
            ]
        assert record.seqinfo == (153, 17238, "01153CF30C2DEDFF")

        assert len(record.features) == 0

        assert len(record.references) == 2
        assert record.references[0].authors == "Kano R., Nakamura Y., Watanabe S., Hasegawa A."
        assert record.references[0].title == "Trichophyton mentagrophytes mRNA for ubiquitin."
        assert len(record.references[0].references) == 0
        assert record.references[1].authors == "Kano R."
        assert record.references[1].title == "Microsporum canis mRNA for ubiquitin, complete cds."
        assert len(record.references[1].references) == 0

        # Check the two parsers agree on the essentials
        assert seq_record.seq == record.sequence
        assert seq_record.description == record.description
        assert seq_record.name == record.entry_name
        assert seq_record.id in record.accessions

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        records = list(SeqIO.parse(datafile, "swiss"))

        assert len(records) == 1
        assert isinstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        assert records[0].seq == seq_record.seq
        assert records[0].description == seq_record.description
        assert records[0].name == seq_record.name
        assert records[0].id == seq_record.id

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        assert len(records) == 1
        assert isinstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        assert records[0].sequence == record.sequence
        assert records[0].description == record.description
        assert records[0].entry_name == record.entry_name
        assert records[0].accessions == record.accessions

    def test_sp013(self):
        """Parsing SwissProt file sp013."""
        filename = "sp013"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)
        seq_record = SeqIO.read(datafile, "swiss")

        assert isinstance(seq_record, SeqRecord)

        assert seq_record.id == "P82909"
        assert seq_record.name == "P82909"
        assert seq_record.description == "MITOCHONDRIAL 28S RIBOSOMAL PROTEIN S36 (MRP-S36)."
        assert repr(seq_record.seq) == "Seq('MGSKMASASRVVQVVKPHTPLIRFPDRRDNPKPNVSEALRSAGLPSHSSVISQH...GPE')"

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        assert record.entry_name == "P82909"
        assert record.accessions == ["P82909"]
        assert record.organism_classification == [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Euteleostomi",
                "Mammalia",
                "Eutheria",
                "Primates",
                "Catarrhini",
                "Hominidae",
                "Homo",
            ]
        assert record.seqinfo == (102, 11335, "83EF107B42E2FCFD")

        assert len(record.features) == 0

        assert len(record.references) == 2
        assert record.references[0].authors == "Strausberg R."
        assert record.references[0].title == ""
        assert len(record.references[0].references) == 0
        assert record.references[1].authors == "Koc E.C., Burkhart W., Blackburn K., Moseley A., Spremulli L.L."
        assert record.references[1].title == "The small subunit of the mammalian mitochondrial ribosome. Identification of the full complement ribosomal proteins present."
        assert len(record.references[1].references) == 0

        # Check the two parsers agree on the essentials
        assert seq_record.seq == record.sequence
        assert seq_record.description == record.description
        assert seq_record.name == record.entry_name
        assert seq_record.id in record.accessions

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        records = list(SeqIO.parse(datafile, "swiss"))

        assert len(records) == 1
        assert isinstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        assert records[0].seq == seq_record.seq
        assert records[0].description == seq_record.description
        assert records[0].name == seq_record.name
        assert records[0].id == seq_record.id

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        assert len(records) == 1
        assert isinstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        assert records[0].sequence == record.sequence
        assert records[0].description == record.description
        assert records[0].entry_name == record.entry_name
        assert records[0].accessions == record.accessions

    def test_P60137(self):
        """Parsing SwissProt file P60137.txt."""
        filename = "P60137.txt"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)
        seq_record = SeqIO.read(datafile, "swiss")

        assert isinstance(seq_record, SeqRecord)

        assert seq_record.id == "P60137"
        assert seq_record.name == "PSBL_ORYSJ"
        assert seq_record.description == "RecName: Full=Photosystem II reaction center protein L {ECO:0000255|HAMAP-Rule:MF_01317}; Short=PSII-L {ECO:0000255|HAMAP-Rule:MF_01317};"
        assert repr(seq_record.seq) == "Seq('MTQSNPNEQNVELNRTSLYWGLLLIFVLAVLFSNYFFN')"

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        assert record.entry_name == "PSBL_ORYSJ"
        assert record.accessions == ["P60137", "O47030", "P12166", "P12167", "Q34007"]
        assert record.organism_classification == [
                "Eukaryota",
                "Viridiplantae",
                "Streptophyta",
                "Embryophyta",
                "Tracheophyta",
                "Spermatophyta",
                "Magnoliopsida",
                "Liliopsida",
                "Poales",
                "Poaceae",
                "BOP clade",
                "Oryzoideae",
                "Oryzeae",
                "Oryzinae",
                "Oryza",
                "Oryza sativa",
            ]
        assert record.seqinfo == (38, 4497, "55537AEC50D25E8D")

        assert len(record.features) == 2
        feature = record.features[0]
        assert feature.type == "CHAIN"
        assert feature.location.start == 0
        assert feature.location.end == 38
        assert feature.qualifiers["note"] == "Photosystem II reaction center protein L"
        assert feature.id == "PRO_0000219754"
        feature = record.features[1]
        assert feature.type == "TRANSMEM"
        assert feature.location.start == 16
        assert feature.location.end == 37
        assert feature.qualifiers["note"] == "Helical"
        assert feature.qualifiers["evidence"] == "ECO:0000255|HAMAP-Rule:MF_01317"

        assert len(record.references) == 2
        assert record.references[0].authors == "Hiratsuka J., Shimada H., Whittier R., Ishibashi T., Sakamoto M., Mori M., Kondo C., Honji Y., Sun C.-R., Meng B.-Y., Li Y.-Q., Kanno A., Nishizawa Y., Hirai A., Shinozaki K., Sugiura M."
        assert record.references[0].title == "The complete sequence of the rice (Oryza sativa) chloroplast genome: intermolecular recombination between distinct tRNA genes accounts for a major plastid DNA inversion during the evolution of the cereals."
        assert len(record.references[0].references) == 2
        assert record.references[0].references[0] == ("PubMed", "2770692")
        assert record.references[0].references[1] == ("DOI", "10.1007/bf02464880")
        assert record.references[1].authors == "Tang J., Xia H., Cao M., Zhang X., Zeng W., Hu S., Tong W., Wang J., Wang J., Yu J., Yang H., Zhu L."
        assert record.references[1].title == "A comparison of rice chloroplast genomes."
        assert len(record.references[1].references) == 2
        assert record.references[1].references[0] == ("PubMed", "15122023")
        assert record.references[1].references[1] == ("DOI", "10.1104/pp.103.031245")

        # Check the two parsers agree on the essentials
        assert seq_record.seq == record.sequence
        assert seq_record.description == record.description
        assert seq_record.name == record.entry_name
        assert seq_record.id in record.accessions

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        records = list(SeqIO.parse(datafile, "swiss"))

        assert len(records) == 1
        assert isinstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        assert records[0].seq == seq_record.seq
        assert records[0].description == seq_record.description
        assert records[0].name == seq_record.name
        assert records[0].id == seq_record.id

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        assert len(records) == 1
        assert isinstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        assert records[0].sequence == record.sequence
        assert records[0].description == record.description
        assert records[0].entry_name == record.entry_name
        assert records[0].accessions == record.accessions

    def test_sp015(self):
        """Parsing SwissProt file sp015."""
        filename = "sp015"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        assert isinstance(seq_record, SeqRecord)

        assert seq_record.id == "IPI00383150"
        assert seq_record.name == "IPI00383150.2"
        assert seq_record.description == ""
        assert repr(seq_record.seq) == "Seq('MSFQAPRRLLELAGQSLLRDQALAISVLDELPRELFPRLFVEAFTSRRCEVLKV...TPC')"

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        assert record.entry_name == "IPI00383150.2"
        assert record.accessions == ["IPI00383150"]
        assert record.organism_classification == [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Euteleostomi",
                "Mammalia",
                "Eutheria",
                "Primates",
                "Catarrhini",
                "Hominidae",
                "Homo",
            ]
        assert record.seqinfo == (457, 52856, "5C3151AAADBDE232")

        assert len(record.features) == 0
        assert len(record.references) == 0

        # Check the two parsers agree on the essentials
        assert seq_record.seq == record.sequence
        assert seq_record.description == record.description
        assert seq_record.name == record.entry_name
        assert seq_record.id in record.accessions

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        with open(datafile) as test_handle:
            records = list(SeqIO.parse(test_handle, "swiss"))

        assert len(records) == 1
        assert isinstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        assert records[0].seq == seq_record.seq
        assert records[0].description == seq_record.description
        assert records[0].name == seq_record.name
        assert records[0].id == seq_record.id

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        assert len(records) == 1
        assert isinstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        assert records[0].sequence == record.sequence
        assert records[0].description == record.description
        assert records[0].entry_name == record.entry_name
        assert records[0].accessions == record.accessions

    def test_P0CK95(self):
        """Parsing SwissProt file P0CK95.txt."""
        filename = "P0CK95.txt"
        datafile = os.path.join("SwissProt", filename)
        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)
        # Check the simple variant
        assert record.features[5].qualifiers["note"] == "N -> G (in strain: O15:H- / 83/39 /ETEC)"
        # Check a FT where the 2nd line starts with /
        assert record.features[6].qualifiers["note"] == "DGTPLPEFYSE -> EGELPKFFSD (in strain: O15:H- / 83/39 / ETEC)"

    def test_Q7Z739(self):
        """Parsing SwissProt file Q7Z739.txt, which has new qualifiers for ligands from Uniprot version 2022_03."""
        filename = "Q7Z739.txt"
        datafile = os.path.join("SwissProt", filename)
        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)
        assert record.gene_name == [
                {
                    "Name": "YTHDF3 {ECO:0000303|PubMed:28106072, ECO:0000312|HGNC:HGNC:26465}"
                }
            ]
        # Check the new ligand feature
        assert record.features[10].qualifiers["ligand"] == "RNA"
        assert record.features[10].qualifiers["ligand_id"] == "ChEBI:CHEBI:33697"
        assert record.features[10].qualifiers["ligand_part"] == "N(6)-methyladenosine 5'-phosphate residue"
        assert record.features[10].qualifiers["ligand_part_id"] == "ChEBI:CHEBI:74449"

    def test_ft_line(self):
        """Parsing SwissProt file O23729, which has a new-style FT line."""
        filename = "O23729.txt"
        datafile = os.path.join("SwissProt", filename)
        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
