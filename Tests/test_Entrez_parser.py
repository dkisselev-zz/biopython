# Copyright 2008-2010 by Michiel de Hoon.  All rights reserved.
# Revisions copyright 2009-2016 by Peter Cock. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Testing code for Bio.Entrez parsers."""

import os
import pickle
import unittest
import pytest
from io import BytesIO

from Bio import Entrez
from Bio import StreamModeError
from Bio.Entrez.Parser import NoneElement


class GeneralTests(unittest.TestCase):
    """General tests for Bio.Entrez."""

    def test_closed_file(self):
        """Test parsing closed file fails gracefully."""
        stream = open("Entrez/taxonomy.xml", "rb")
        stream.close()
        with pytest.raises(ValueError):
            Entrez.read(stream)

    def test_read_bytes_stream(self):
        """Test reading a file opened in binary mode."""
        with open("Entrez/taxonomy.xml", "rb") as stream:
            records = Entrez.read(stream)
        assert len(records) == 2
        for record in records:
            assert "TaxId" in record
            assert "ScientificName" in record

    def test_parse_bytes_stream(self):
        """Test parsing a file opened in binary mode."""
        with open("Entrez/taxonomy.xml", "rb") as stream:
            records = Entrez.parse(stream)
            n = 0
            for record in records:
                assert "TaxId" in record
                assert "ScientificName" in record
                n += 1
        assert n == 2

    def test_read_text_file(self):
        """Test reading a file opened in text mode."""
        message = "^the XML file must be opened in binary mode.$"
        with open("Entrez/taxonomy.xml") as stream:
            with pytest.raises(StreamModeError, match=message):
                Entrez.read(stream)

    def test_parse_text_file(self):
        """Test parsing a file opened in text mode."""
        message = "^the XML file must be opened in binary mode.$"
        with open("Entrez/taxonomy.xml") as stream:
            records = Entrez.parse(stream)
            with pytest.raises(StreamModeError, match=message):
                next(records)

    def test_BytesIO(self):
        """Test parsing a BytesIO stream (bytes not string)."""
        with open("Entrez/taxonomy.xml", "rb") as stream:
            data = stream.read()
        stream = BytesIO(data)
        records = Entrez.read(stream)
        assert len(records) == 2
        assert records[0]["ScientificName"] == "Canis lupus familiaris"
        assert records[1]["ScientificName"] == "Felis catus"
        stream.close()

    def test_pickle(self):
        """Test if records created by the parser can be pickled."""
        directory = "Entrez"
        filenames = os.listdir(directory)
        for filename in sorted(filenames):
            basename, extension = os.path.splitext(filename)
            if extension != ".xml":
                continue
            if filename in (
                "biosample.xml",  # DTD not specified in XML file
                "einfo4.xml",  # XML corrupted
                "journals.xml",  # Missing XML declaration
            ):
                continue
            path = os.path.join(directory, filename)
            with open(path, "rb") as stream:
                if filename in ("epost2.xml", "epost3.xml", "esummary8.xml"):
                    # these include an ErrorElement
                    record = Entrez.read(stream, ignore_errors=True)
                else:
                    record = Entrez.read(stream)
            with BytesIO() as stream:
                pickle.dump(record, stream)
                stream.seek(0)
                pickled_record = pickle.load(stream)
            assert record == pickled_record, filename


class EInfoTest(unittest.TestCase):
    """Tests for parsing XML output returned by EInfo."""

    def test_list(self):
        """Test parsing database list returned by EInfo."""
        # To create the XML file, use
        # >>> Bio.Entrez.einfo()
        with open("Entrez/einfo1.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert record["DbList"] == [
                "pubmed",
                "protein",
                "nuccore",
                "ipg",
                "nucleotide",
                "structure",
                "genome",
                "annotinfo",
                "assembly",
                "bioproject",
                "biosample",
                "blastdbinfo",
                "books",
                "cdd",
                "clinvar",
                "gap",
                "gapplus",
                "grasp",
                "dbvar",
                "gene",
                "gds",
                "geoprofiles",
                "medgen",
                "mesh",
                "nlmcatalog",
                "omim",
                "orgtrack",
                "pmc",
                "proteinclusters",
                "pcassay",
                "protfam",
                "pccompound",
                "pcsubstance",
                "seqannot",
                "snp",
                "sra",
                "taxonomy",
                "biocollections",
                "gtr",
            ]

    def test_pubmed1(self):
        """Test parsing database info returned by EInfo."""
        # To create the XML file, use
        # >>> Bio.Entrez.einfo(db="pubmed")
        with open("Entrez/einfo2.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert len(record["DbInfo"]) == 1
        assert record["DbInfo"][0]["DbName"] == "pubmed"
        assert record["DbInfo"][0]["MenuName"] == "PubMed"
        assert record["DbInfo"][0]["Description"] == "PubMed bibliographic record"
        assert record["DbInfo"][0]["Count"] == "39730388"
        assert record["DbInfo"][0]["LastUpdate"] == "2025/11/27 06:33"

        assert len(record["DbInfo"][0]["FieldList"]) == 50

        assert record["DbInfo"][0]["FieldList"][0]["Name"] == "ALL"
        assert record["DbInfo"][0]["FieldList"][0]["FullName"] == "All Fields"
        assert record["DbInfo"][0]["FieldList"][0]["Description"] == "All terms from all searchable fields"

        assert record["DbInfo"][0]["FieldList"][0]["IsNumerical"] == "N"
        assert record["DbInfo"][0]["FieldList"][0]["SingleToken"] == "N"
        assert record["DbInfo"][0]["FieldList"][0]["Hierarchy"] == "N"
        assert record["DbInfo"][0]["FieldList"][0]["IsHidden"] == "N"

        assert len(record["DbInfo"][0]["LinkList"]) == 57

        assert record["DbInfo"][0]["LinkList"][0]["Name"] == "pubmed_assembly"
        assert record["DbInfo"][0]["LinkList"][0]["Menu"] == "Assembly"
        assert record["DbInfo"][0]["LinkList"][0]["Description"] == "Assembly"
        assert record["DbInfo"][0]["LinkList"][0]["DbTo"] == "assembly"

        assert record["DbInfo"][0]["LinkList"][56]["Name"] == "pubmed_taxonomy_entrez"
        assert record["DbInfo"][0]["LinkList"][56]["Menu"] == "Taxonomy via GenBank"
        assert record["DbInfo"][0]["LinkList"][56]["Description"] == "Related Taxonomy entry computed using other Entrez links"
        assert record["DbInfo"][0]["LinkList"][56]["DbTo"] == "taxonomy"

    def test_corrupted(self):
        """Test if corrupted XML is handled correctly."""
        # To create the XML file, use
        # >>> Bio.Entrez.einfo()
        # and manually delete the last couple of lines
        from Bio.Entrez import Parser

        with open("Entrez/einfo4.xml", "rb") as stream:
            with pytest.raises(Parser.CorruptedXMLError):
                Entrez.read(stream)


class ESearchTest(unittest.TestCase):
    """Tests for parsing XML output returned by ESearch."""

    def test_pubmed1(self):
        """Test parsing XML returned by ESearch from PubMed (first test)."""
        # To create the XML file, use
        # >>> Bio.Entrez.esearch(db="pubmed", term="biopython")
        with open("Entrez/esearch1.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert record["Count"] == "63"
        assert record["RetMax"] == "20"
        assert record["RetStart"] == "0"
        assert len(record["IdList"]) == 20
        assert record["IdList"][0] == "41282813"
        assert record["IdList"][1] == "41148224"
        assert record["IdList"][2] == "41011574"
        assert record["IdList"][3] == "40959146"
        assert record["IdList"][4] == "40937394"
        assert record["IdList"][5] == "40657423"
        assert record["IdList"][6] == "40651330"
        assert record["IdList"][7] == "40572159"
        assert record["IdList"][8] == "40160861"
        assert record["IdList"][9] == "39883659"
        assert record["IdList"][10] == "39882099"
        assert record["IdList"][11] == "39717221"
        assert record["IdList"][12] == "39546778"
        assert record["IdList"][13] == "39507944"
        assert record["IdList"][14] == "39445816"
        assert record["IdList"][15] == "38808697"
        assert record["IdList"][16] == "38650605"
        assert record["IdList"][17] == "38365590"
        assert record["IdList"][18] == "38235175"
        assert record["IdList"][19] == "37810457"
        assert len(record["TranslationSet"]) == 0
        assert record["QueryTranslation"] == '"biopython"[All Fields]'

    def test_pubmed2(self):
        """Test parsing XML returned by ESearch from PubMed (second test)."""
        # Search in PubMed for the term cancer for the entrez date from
        # the last 60 days and retrieve the first 100 IDs and translations
        # using the history parameter.
        # To create the XML file, use
        # >>> Bio.Entrez.esearch(db="pubmed", term="cancer", reldate=60,
        #                        datetype="edat", retmax=100, usehistory="y")
        with open("Entrez/esearch2.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert record["Count"] == "42249"
        assert record["RetMax"] == "100"
        assert record["RetStart"] == "0"
        assert record["QueryKey"] == "1"
        assert record["WebEnv"] == "MCID_6927d6e7fee3e90f880ec190"
        assert len(record["IdList"]) == 100
        assert record["IdList"][0] == "41297076"
        assert record["IdList"][1] == "41297074"
        assert record["IdList"][2] == "41297068"
        assert record["IdList"][3] == "41297055"
        assert record["IdList"][4] == "41297046"
        assert record["IdList"][5] == "41297043"
        assert record["IdList"][6] == "41297032"
        assert record["IdList"][7] == "41297028"
        assert record["IdList"][8] == "41297023"
        assert record["IdList"][9] == "41297015"
        assert record["IdList"][10] == "41297014"
        assert record["IdList"][11] == "41296995"
        assert record["IdList"][12] == "41296994"
        assert record["IdList"][13] == "41296991"
        assert record["IdList"][14] == "41296990"
        assert record["IdList"][15] == "41296985"
        assert record["IdList"][16] == "41296927"
        assert record["IdList"][17] == "41296920"
        assert record["IdList"][18] == "41296918"
        assert record["IdList"][19] == "41296886"
        assert record["IdList"][20] == "41296881"
        assert record["IdList"][21] == "41296866"
        assert record["IdList"][22] == "41296862"
        assert record["IdList"][23] == "41296854"
        assert record["IdList"][24] == "41296842"
        assert record["IdList"][25] == "41296829"
        assert record["IdList"][26] == "41296826"
        assert record["IdList"][27] == "41296824"
        assert record["IdList"][28] == "41296822"
        assert record["IdList"][29] == "41296821"
        assert record["IdList"][30] == "41296811"
        assert record["IdList"][31] == "41296783"
        assert record["IdList"][32] == "41296773"
        assert record["IdList"][33] == "41296763"
        assert record["IdList"][34] == "41296753"
        assert record["IdList"][35] == "41296747"
        assert record["IdList"][36] == "41296732"
        assert record["IdList"][37] == "41296731"
        assert record["IdList"][38] == "41296714"
        assert record["IdList"][39] == "41296711"
        assert record["IdList"][40] == "41296710"
        assert record["IdList"][41] == "41296689"
        assert record["IdList"][42] == "41296685"
        assert record["IdList"][43] == "41296683"
        assert record["IdList"][44] == "41296675"
        assert record["IdList"][45] == "41296667"
        assert record["IdList"][46] == "41296665"
        assert record["IdList"][47] == "41296663"
        assert record["IdList"][48] == "41296658"
        assert record["IdList"][49] == "41296656"
        assert record["IdList"][50] == "41296616"
        assert record["IdList"][51] == "41296609"
        assert record["IdList"][52] == "41296608"
        assert record["IdList"][53] == "41296606"
        assert record["IdList"][54] == "41296594"
        assert record["IdList"][55] == "41296590"
        assert record["IdList"][56] == "41296586"
        assert record["IdList"][57] == "41296583"
        assert record["IdList"][58] == "41296580"
        assert record["IdList"][59] == "41296569"
        assert record["IdList"][60] == "41296567"
        assert record["IdList"][61] == "41296563"
        assert record["IdList"][62] == "41296561"
        assert record["IdList"][63] == "41296555"
        assert record["IdList"][64] == "41296544"
        assert record["IdList"][65] == "41296517"
        assert record["IdList"][66] == "41296514"
        assert record["IdList"][67] == "41296508"
        assert record["IdList"][68] == "41296505"
        assert record["IdList"][69] == "41296501"
        assert record["IdList"][70] == "41296494"
        assert record["IdList"][71] == "41296485"
        assert record["IdList"][72] == "41296462"
        assert record["IdList"][73] == "41296458"
        assert record["IdList"][74] == "41296454"
        assert record["IdList"][75] == "41296452"
        assert record["IdList"][76] == "41296443"
        assert record["IdList"][77] == "41296441"
        assert record["IdList"][78] == "41296440"
        assert record["IdList"][79] == "41296437"
        assert record["IdList"][80] == "41296436"
        assert record["IdList"][81] == "41296434"
        assert record["IdList"][82] == "41296433"
        assert record["IdList"][83] == "41296432"
        assert record["IdList"][84] == "41296430"
        assert record["IdList"][85] == "41296429"
        assert record["IdList"][86] == "41296425"
        assert record["IdList"][87] == "41296417"
        assert record["IdList"][88] == "41296413"
        assert record["IdList"][89] == "41296412"
        assert record["IdList"][90] == "41296408"
        assert record["IdList"][91] == "41296406"
        assert record["IdList"][92] == "41296402"
        assert record["IdList"][93] == "41296396"
        assert record["IdList"][94] == "41296387"
        assert record["IdList"][95] == "41296385"
        assert record["IdList"][96] == "41296382"
        assert record["IdList"][97] == "41296379"
        assert record["IdList"][98] == "41296369"
        assert record["IdList"][99] == "41296368"
        assert len(record["TranslationSet"]) == 1
        assert record["TranslationSet"][0]["From"] == "cancer"
        assert record["TranslationSet"][0]["To"] == '"cancer\'s"[All Fields] OR "cancerated"[All Fields] OR "canceration"[All Fields] OR "cancerization"[All Fields] OR "cancerized"[All Fields] OR "cancerous"[All Fields] OR "neoplasms"[MeSH Terms] OR "neoplasms"[All Fields] OR "cancer"[All Fields] OR "cancers"[All Fields]'
        assert record["QueryTranslation"] == '("cancer s"[All Fields] OR "cancerated"[All Fields] OR "canceration"[All Fields] OR "cancerization"[All Fields] OR "cancerized"[All Fields] OR "cancerous"[All Fields] OR "neoplasms"[MeSH Terms] OR "neoplasms"[All Fields] OR "cancer"[All Fields] OR "cancers"[All Fields]) AND 2025/09/27:2025/11/26[Date - Entry]'

    def test_pubmed3(self):
        """Test parsing XML returned by ESearch from PubMed (third test)."""
        # Search in PubMed for the journal PNAS Volume 97, and retrieve
        # 6 IDs starting at ID 7.
        # To create the XML file, use
        # >>> Bio.Entrez.esearch(db="pubmed", term="PNAS[ta] AND 97[vi]",
        #                        retstart=6, retmax=6)
        with open("Entrez/esearch3.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert record["Count"] == "2651"
        assert record["RetMax"] == "6"
        assert record["RetStart"] == "6"
        assert len(record["IdList"]) == 6
        assert record["IdList"][0] == "11121077"
        assert record["IdList"][1] == "11121076"
        assert record["IdList"][2] == "11121075"
        assert record["IdList"][3] == "11121074"
        assert record["IdList"][4] == "11121073"
        assert record["IdList"][5] == "11121072"
        assert len(record["TranslationSet"]) == 1
        assert record["TranslationSet"][0]["From"] == "PNAS[ta]"
        assert record["TranslationSet"][0]["To"] == '"Proc Natl Acad Sci U S A"[Journal:__jid7505876]'
        assert record["QueryTranslation"] == '"proc natl acad sci u s a"[Journal] AND "97"[Volume]'

    def test_pmc(self):
        """Test parsing XML returned by ESearch from PubMed Central."""
        # Search in PubMed Central for stem cells in articles with an abstract.
        # To create the XML file, use
        # >>> Bio.Entrez.esearch(db="pmc", term="stem cells AND hasabstract")
        with open("Entrez/esearch5.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert record["Count"] == "5"
        assert record["RetMax"] == "5"
        assert record["RetStart"] == "0"
        assert len(record["IdList"]) == 5
        assert record["IdList"][0] == "6540636"
        assert record["IdList"][1] == "5830998"
        assert record["IdList"][2] == "4460425"
        assert record["IdList"][3] == "4466840"
        assert record["IdList"][4] == "3443459"
        assert len(record["TranslationSet"]) == 1
        assert record["TranslationSet"][0]["From"] == "stem cells"
        assert record["TranslationSet"][0]["To"] == '"stem cells"[MeSH Terms] OR ("stem"[All Fields] AND "cells"[All Fields]) OR "stem cells"[All Fields]'
        assert len(record["TranslationStack"]) == 11
        assert record["TranslationStack"][0]["Term"] == '"stem cells"[MeSH Terms]'
        assert record["TranslationStack"][0]["Field"] == "MeSH Terms"
        assert record["TranslationStack"][0]["Count"] == "109831"
        assert record["TranslationStack"][0]["Explode"] == "Y"
        assert record["TranslationStack"][0].tag == "TermSet"
        assert record["TranslationStack"][1]["Term"] == '"stem"[All Fields]'
        assert record["TranslationStack"][1]["Field"] == "All Fields"
        assert record["TranslationStack"][1]["Count"] == "1528859"
        assert record["TranslationStack"][1]["Explode"] == "N"
        assert record["TranslationStack"][1].tag == "TermSet"
        assert record["TranslationStack"][2] == {
                "Term": '"cells"[All Fields]',
                "Field": "All Fields",
                "Count": "5033502",
                "Explode": "N",
            }
        assert record["TranslationStack"][2].tag == "TermSet"
        assert record["TranslationStack"][3] == "AND"
        assert record["TranslationStack"][3].tag == "OP"
        assert record["TranslationStack"][4] == "GROUP"
        assert record["TranslationStack"][4].tag == "OP"
        assert record["TranslationStack"][5] == "OR"
        assert record["TranslationStack"][5].tag == "OP"
        assert record["TranslationStack"][6] == {
                "Term": '"stem cells"[All Fields]',
                "Field": "All Fields",
                "Count": "715386",
                "Explode": "N",
            }
        assert record["TranslationStack"][6].tag == "TermSet"
        assert record["TranslationStack"][7] == "OR"
        assert record["TranslationStack"][7].tag == "OP"
        assert record["TranslationStack"][8] == "GROUP"
        assert record["TranslationStack"][8].tag == "OP"
        assert record["TranslationStack"][9]["Term"] == "hasabstract[All Fields]"
        assert record["TranslationStack"][9]["Field"] == "All Fields"
        assert record["TranslationStack"][9]["Count"] == "83"
        assert record["TranslationStack"][9]["Explode"] == "N"
        assert record["TranslationStack"][9].tag == "TermSet"
        assert record["TranslationStack"][10] == "AND"
        assert record["TranslationStack"][10].tag == "OP"
        assert record["QueryTranslation"] == '("stem cells"[MeSH Terms] OR ("stem"[All Fields] AND "cells"[All Fields]) OR "stem cells"[All Fields]) AND hasabstract[All Fields]'

    def test_nucleotide(self):
        """Test parsing XML returned by ESearch from the Nucleotide database."""
        # Search in Nucleotide for a property of the sequence,
        # To create the XML file, use
        # >>> Bio.Entrez.esearch(db="nucleotide", term="biomol trna[prop]")
        with open("Entrez/esearch6.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert record["Count"] == "1018"
        assert record["RetMax"] == "20"
        assert record["RetStart"] == "0"
        assert len(record["IdList"]) == 20
        assert record["IdList"][0] == "3069587804"
        assert record["IdList"][1] == "3069587774"
        assert record["IdList"][2] == "3069587771"
        assert record["IdList"][3] == "3069587741"
        assert record["IdList"][4] == "3069402178"
        assert record["IdList"][5] == "3069218576"
        assert record["IdList"][6] == "2737963026"
        assert record["IdList"][7] == "2586967820"
        assert record["IdList"][8] == "2274792564"
        assert record["IdList"][9] == "2274792563"
        assert record["IdList"][10] == "2274792562"
        assert record["IdList"][11] == "2274792561"
        assert record["IdList"][12] == "2274792560"
        assert record["IdList"][13] == "2274792559"
        assert record["IdList"][14] == "2274792558"
        assert record["IdList"][15] == "2274792557"
        assert record["IdList"][16] == "2274792556"
        assert record["IdList"][17] == "2274792555"
        assert record["IdList"][18] == "2274792554"
        assert record["IdList"][19] == "2274792553"
        assert len(record["TranslationSet"]) == 0
        assert record["QueryTranslation"] == "biomol trna[prop]"

    def test_protein(self):
        """Test parsing XML returned by ESearch from the Protein database."""
        # Search in Protein for a molecular weight
        # To create the XML file, use
        # >>> Bio.Entrez.esearch(db="protein", term="200020[molecular weight]")
        with open("Entrez/esearch7.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert record["Count"] == "330"
        assert record["RetMax"] == "20"
        assert record["RetStart"] == "0"
        assert len(record["IdList"]) == 20
        assert record["IdList"][0] == "3108048616"
        assert record["IdList"][1] == "3076096315"
        assert record["IdList"][2] == "3070972985"
        assert record["IdList"][3] == "3070963398"
        assert record["IdList"][4] == "3042783758"
        assert record["IdList"][5] == "3042778745"
        assert record["IdList"][6] == "3042773739"
        assert record["IdList"][7] == "3042768843"
        assert record["IdList"][8] == "3042764246"
        assert record["IdList"][9] == "3042754755"
        assert record["IdList"][10] == "3042750053"
        assert record["IdList"][11] == "3037813595"
        assert record["IdList"][12] == "3032309824"
        assert record["IdList"][13] == "1468867753"
        assert record["IdList"][14] == "2994907531"
        assert record["IdList"][15] == "2994891115"
        assert record["IdList"][16] == "2993271659"
        assert record["IdList"][17] == "2993266811"
        assert record["IdList"][18] == "2993262038"
        assert record["IdList"][19] == "2174081062"
        assert len(record["TranslationSet"]) == 0
        assert len(record["TranslationStack"]) == 2
        assert record["TranslationStack"][0]["Term"] == "000200020[molecular weight]"
        assert record["TranslationStack"][0]["Field"] == "molecular weight"
        assert record["TranslationStack"][0]["Count"] == "330"
        assert record["TranslationStack"][0]["Explode"] == "N"
        assert record["TranslationStack"][0].tag == "TermSet"
        assert record["TranslationStack"][1] == "GROUP"
        assert record["TranslationStack"][1].tag == "OP"
        assert record["QueryTranslation"] == "000200020[molecular weight]"

    def test_notfound(self):
        """Test parsing XML returned by ESearch when no items were found."""
        # To create the XML file, use
        # >>> Bio.Entrez.esearch(db="protein", term="abcXYZ")
        with open("Entrez/esearch8.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert record["Count"] == "0"
        assert record["RetMax"] == "0"
        assert record["RetStart"] == "0"
        assert len(record["IdList"]) == 0
        assert len(record["TranslationSet"]) == 0
        assert record["QueryTranslation"] == "(abcXYZ[All Fields])"
        assert len(record["ErrorList"]) == 2
        assert "PhraseNotFound" in record["ErrorList"]
        assert "FieldNotFound" in record["ErrorList"]
        assert len(record["ErrorList"]["PhraseNotFound"]) == 1
        assert len(record["ErrorList"]["FieldNotFound"]) == 0
        assert record["ErrorList"]["PhraseNotFound"][0] == "abcXYZ"
        assert len(record["WarningList"]) == 3
        assert "PhraseIgnored" in record["WarningList"]
        assert "QuotedPhraseNotFound" in record["WarningList"]
        assert "OutputMessage" in record["WarningList"]
        assert len(record["WarningList"]["PhraseIgnored"]) == 0
        assert len(record["WarningList"]["QuotedPhraseNotFound"]) == 0
        assert len(record["WarningList"]["OutputMessage"]) == 1
        assert record["WarningList"]["OutputMessage"][0] == "No items found."


class EPostTest(unittest.TestCase):
    """Tests for parsing XML output returned by EPost."""

    # Don't know how to get an InvalidIdList in the XML returned by EPost;
    # unable to test if we are parsing it correctly.
    def test_epost(self):
        """Test parsing XML returned by EPost."""
        # To create the XML file, use
        # >>> Bio.Entrez.epost(db="pubmed", id="11237011")
        with open("Entrez/epost1.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert record["QueryKey"] == "1"
        assert record["WebEnv"] == "MCID_692851c130aec64bed0a8c2a"

    def test_wrong(self):
        """Test parsing XML returned by EPost with incorrect arguments."""
        # To create the XML file, use
        # >>> Bio.Entrez.epost(db="nothing")
        with open("Entrez/epost2.xml", "rb") as stream:
            with pytest.raises(RuntimeError):
                Entrez.read(stream)
        with open("Entrez/epost2.xml", "rb") as stream:
            record = Entrez.read(stream, ignore_errors=True)
        assert len(record) == 1
        assert len(record.attributes) == 0
        assert record["ERROR"] == "Invalid db name specified: nothing"
        assert record["ERROR"].tag == "ERROR"

    def test_invalid(self):
        """Test parsing XML returned by EPost with invalid id (overflow tag)."""
        # To create the XML file, use
        # >>> Bio.Entrez.epost(db="pubmed", id=99999999999999999999999999999999)
        with open("Entrez/epost3.xml", "rb") as stream:
            with pytest.raises(RuntimeError):
                Entrez.read(stream)
        with open("Entrez/epost3.xml", "rb") as stream:
            record = Entrez.read(stream, ignore_errors=True)
        assert len(record) == 1
        assert len(record.attributes) == 0
        assert record["ERROR"] == "Some IDs have invalid value and were omitted. Maximum ID value 18446744073709551615"
        assert record["ERROR"].tag == "ERROR"
        # Note that the first ERROR element is lost. Strictly speaking, the XML
        # is not consistent with the DTD, which allows only one ERROR element.


class ESummaryTest(unittest.TestCase):
    """Tests for parsing XML output returned by ESummary."""

    # Items have a type, which can be
    # (Integer|Date|String|Structure|List|Flags|Qualifier|Enumerator|Unknown)
    # I don't have an XML file where the type "Flags", "Qualifier",
    # "Enumerator", or "Unknown" is used, so they are not tested here.
    def test_pubmed(self):
        """Test parsing XML returned by ESummary from PubMed."""
        # In PubMed display records for PMIDs 11850928 and 11482001 in
        # xml retrieval mode
        # To create the XML file, use
        # >>> Bio.Entrez.esummary(db="pubmed", id=["11850928","11482001"],
        #                         retmode="xml")
        with open("Entrez/esummary1.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert record[0]["Id"] == "11850928"
        assert record[0]["PubDate"] == "1965 Aug"
        assert record[0]["EPubDate"] == ""
        assert record[0]["Source"] == "Arch Dermatol"
        assert len(record[0]["AuthorList"]) == 2
        assert record[0]["AuthorList"][0] == "LoPresti PJ"
        assert record[0]["AuthorList"][1] == "Hambrick GW Jr"
        assert record[0]["LastAuthor"] == "Hambrick GW Jr"
        assert record[0]["Title"] == "Zirconium granuloma following treatment of rhus dermatitis."
        assert record[0]["Volume"] == "92"
        assert record[0]["Issue"] == "2"
        assert record[0]["Pages"] == "188-91"
        assert record[0]["LangList"] == ["English"]
        assert record[0]["NlmUniqueID"] == "0372433"
        assert record[0]["ISSN"] == "0003-987X"
        assert record[0]["ESSN"] == ""
        assert len(record[0]["PubTypeList"]) == 1
        assert record[0]["PubTypeList"][0] == "Journal Article"
        assert record[0]["RecordStatus"] == "PubMed - indexed for MEDLINE"
        assert record[0]["PubStatus"] == "ppublish"
        assert len(record[0]["ArticleIds"]) == 2
        assert record[0]["ArticleIds"]["pubmed"] == ["11850928"]
        assert record[0]["ArticleIds"]["medline"] == []
        assert len(record[0]["History"]) == 3
        assert record[0]["History"]["pubmed"] == ["1965/08/01 00:00"]
        assert record[0]["History"]["medline"] == ["2002/03/09 10:01"]
        assert record[0]["History"]["entrez"] == "1965/08/01 00:00"
        assert len(record[0]["References"]) == 0
        assert record[0]["HasAbstract"] == 1
        assert record[0]["PmcRefCount"] == 0
        assert record[0]["FullJournalName"] == "Archives of dermatology"
        assert record[0]["ELocationID"] == ""
        assert record[0]["SO"] == "1965 Aug;92(2):188-91"

        assert record[1]["Id"] == "11482001"
        assert record[1]["PubDate"] == "2001 Jun"
        assert record[1]["EPubDate"] == ""
        assert record[1]["Source"] == "Adverse Drug React Toxicol Rev"
        assert len(record[1]["AuthorList"]) == 3
        assert record[1]["AuthorList"][0] == "Mantle D"
        assert record[1]["AuthorList"][1] == "Gok MA"
        assert record[1]["AuthorList"][2] == "Lennard TW"
        assert record[1]["LastAuthor"] == "Lennard TW"
        assert record[1]["Title"] == "Adverse and beneficial effects of plant extracts on skin and skin disorders."
        assert record[1]["Volume"] == "20"
        assert record[1]["Issue"] == "2"
        assert record[1]["Pages"] == "89-103"
        assert len(record[1]["LangList"]) == 1
        assert record[1]["LangList"][0] == "English"
        assert record[1]["NlmUniqueID"] == "9109474"
        assert record[1]["ISSN"] == "0964-198X"
        assert record[1]["ESSN"] == ""
        assert len(record[1]["PubTypeList"]) == 2
        assert record[1]["PubTypeList"][0] == "Journal Article"
        assert record[1]["PubTypeList"][1] == "Review"
        assert record[1]["RecordStatus"] == "PubMed - indexed for MEDLINE"
        assert record[1]["PubStatus"] == "ppublish"
        assert len(record[1]["ArticleIds"]) == 2
        assert record[1]["ArticleIds"]["pubmed"] == ["11482001"]
        assert record[1]["ArticleIds"]["medline"] == []
        assert len(record[1]["History"]) == 3
        assert record[1]["History"]["pubmed"] == ["2001/08/03 10:00"]
        assert record[1]["History"]["medline"] == ["2002/01/23 10:01"]
        assert record[1]["History"]["entrez"] == "2001/08/03 10:00"
        assert len(record[1]["References"]) == 0
        assert record[1]["HasAbstract"] == 1
        assert record[1]["PmcRefCount"] == 0
        assert record[1]["FullJournalName"] == "Adverse drug reactions and toxicological reviews"
        assert record[1]["ELocationID"] == ""
        assert record[1]["SO"] == "2001 Jun;20(2):89-103"

    def test_protein(self):
        """Test parsing XML returned by ESummary from the Protein database."""
        # In Protein display records for GIs 28800982 and 28628843 in xml retrieval mode
        # To create the XML file, use
        # >>> Bio.Entrez.esummary(db="protein", id="28800982,28628843", retmode="xml")
        with open("Entrez/esummary3.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert record[0]["Id"] == "28800982"
        assert record[0]["Caption"] == "AAO47091"
        assert record[0]["Title"] == "hemochromatosis, partial [Homo sapiens]"
        assert record[0]["Extra"] == "gi|28800982|gb|AAO47091.1|[28800982]"
        assert record[0]["Gi"] == 28800982
        assert record[0]["CreateDate"] == "2003/03/03"
        assert record[0]["UpdateDate"] == "2016/07/25"
        assert record[0]["Flags"] == NoneElement(None, None, None)
        assert record[0]["TaxId"] == 9606
        assert record[0]["Length"] == 268
        assert record[0]["Status"] == "live"
        assert record[0]["ReplacedBy"] == ""
        assert record[0]["Comment"] == "  "

        assert record[1]["Id"] == "28628843"
        assert record[1]["Caption"] == "AAO49381"
        assert record[1]["Title"] == "erythroid associated factor [Homo sapiens]"
        assert record[1]["Extra"] == "gi|28628843|gb|AAO49381.1|AF485325_1[28628843]"
        assert record[1]["Gi"] == 28628843
        assert record[1]["CreateDate"] == "2003/03/02"
        assert record[1]["UpdateDate"] == "2003/03/02"
        assert record[1]["Flags"] == NoneElement(None, None, None)
        assert record[1]["TaxId"] == 9606
        assert record[1]["Length"] == 102
        assert record[1]["Status"] == "live"
        assert record[1]["ReplacedBy"] == ""
        assert record[1]["Comment"] == "  "

    def test_nucleotide(self):
        """Test parsing XML returned by ESummary from the Nucleotide database."""
        # In Nucleotide display records for GIs 28864546 and 28800981
        # in xml retrieval mode
        # To create the XML file, use
        # >>> Bio.Entrez.esummary(db="nucleotide", id="28864546,28800981",
        #                         retmode="xml")
        with open("Entrez/esummary4.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert record[0]["Id"] == "28864546"
        assert record[0]["Caption"] == "AY207443"
        assert record[0]["Title"] == "Homo sapiens alpha hemoglobin (HBZP) pseudogene 3' UTR/AluJo repeat breakpoint junction"
        assert record[0]["Extra"] == "gi|28864546|gb|AY207443.1|[28864546]"
        assert record[0]["Gi"] == 28864546
        assert record[0]["CreateDate"] == "2003/03/05"
        assert record[0]["UpdateDate"] == "2003/03/05"
        assert record[0]["Flags"] == NoneElement(None, None, None)
        assert record[0]["TaxId"] == 9606
        assert record[0]["Length"] == 491
        assert record[0]["Status"] == "live"
        assert record[0]["ReplacedBy"] == ""
        assert record[0]["Comment"] == "  "

        assert record[1]["Id"] == "28800981"
        assert record[1]["Caption"] == "AY205604"
        assert record[1]["Title"] == "Homo sapiens hemochromatosis (HFE) mRNA, partial cds"
        assert record[1]["Extra"] == "gi|28800981|gb|AY205604.1|[28800981]"
        assert record[1]["Gi"] == 28800981
        assert record[1]["CreateDate"] == "2003/03/03"
        assert record[1]["UpdateDate"] == "2016/07/25"
        assert record[1]["Flags"] == NoneElement(None, None, None)
        assert record[1]["TaxId"] == 9606
        assert record[1]["Length"] == 860
        assert record[1]["Status"] == "live"
        assert record[1]["ReplacedBy"] == ""
        assert record[1]["Comment"] == "  "

    def test_structure(self):
        """Test parsing XML returned by ESummary from the Structure database."""
        # In Nucleotide display records for GIs 28864546 and 28800981
        # in xml retrieval mode
        # To create the XML file, use
        # >>> Bio.Entrez.esummary(db="structure", id=["19923","12120"],
        #                         retmode="xml")
        with open("Entrez/esummary5.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert record[0]["Id"] == "19923"
        assert record[0]["PdbAcc"] == "1L5J"
        assert record[0]["PdbDescr"] == "CRYSTAL STRUCTURE OF E. COLI ACONITASE B"
        assert record[0]["EC"] == ""
        assert record[0]["Resolution"] == "2.4"
        assert record[0]["ExpMethod"] == "X-ray Diffraction"
        assert record[0]["PdbClass"] == "LYASE"
        assert record[0]["PdbDepositDate"] == "2002/03/07 00:00"
        assert record[0]["MMDBEntryDate"] == "2002/07/11 00:00"
        assert record[0]["OrganismList"] == ["Escherichia coli"]
        assert record[0]["LigCode"] == "F3S|TRA"
        assert record[0]["LigCount"] == "2"
        assert record[0]["ModProteinResCount"] == "0"
        assert record[0]["ModDNAResCount"] == "0"
        assert record[0]["ModRNAResCount"] == "0"
        assert record[0]["ProteinChainCount"] == ""
        assert record[0]["DNAChainCount"] == ""
        assert record[0]["RNAChainCount"] == ""

        assert record[1]["Id"] == "12120"
        assert record[1]["PdbAcc"] == "1B0K"
        assert record[1]["PdbDescr"] == "S642A:FLUOROCITRATE COMPLEX OF ACONITASE"
        assert record[1]["EC"] == ""
        assert record[1]["Resolution"] == "2.5"
        assert record[1]["ExpMethod"] == "X-ray Diffraction"
        assert record[1]["PdbClass"] == "LYASE"
        assert record[1]["PdbDepositDate"] == "1998/11/11 00:00"
        assert record[1]["MMDBEntryDate"] == "2000/01/24 00:00"
        assert record[1]["OrganismList"] == ["Sus scrofa"]
        assert record[1]["LigCode"] == "FLC|O|SF4"
        assert record[1]["LigCount"] == "3"
        assert record[1]["ModProteinResCount"] == "0"
        assert record[1]["ModDNAResCount"] == "0"
        assert record[1]["ModRNAResCount"] == "0"
        assert record[1]["ProteinChainCount"] == ""
        assert record[1]["DNAChainCount"] == ""
        assert record[1]["RNAChainCount"] == ""

    def test_taxonomy(self):
        """Test parsing XML returned by ESummary from the Taxonomy database."""
        # In Taxonomy display records for TAXIDs 9913 and 30521 in
        # xml retrieval mode
        # To create the XML file, use
        # >>> Bio.Entrez.esummary(db="taxonomy", id=["9913","30521"],
        #                         retmode="xml")
        with open("Entrez/esummary6.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert record[0]["Id"] == "9913"
        assert record[0]["Status"] == "active"
        assert record[0]["Rank"] == "species"
        assert record[0]["Division"] == "even-toed ungulates & whales"
        assert record[0]["ScientificName"] == "Bos taurus"
        assert record[0]["CommonName"] == "domestic cattle"
        assert record[0]["TaxId"] == 9913
        assert record[0]["AkaTaxId"] == 0
        assert record[0]["Genus"] == ""
        assert record[0]["Species"] == ""
        assert record[0]["Subsp"] == ""
        assert record[0]["ModificationDate"] == "2024/08/09 00:00"

        assert record[1]["Id"] == "30521"
        assert record[1]["Status"] == "active"
        assert record[1]["Rank"] == "species"
        assert record[1]["Division"] == "even-toed ungulates & whales"
        assert record[1]["ScientificName"] == "Bos grunniens"
        assert record[1]["CommonName"] == "domestic yak"
        assert record[1]["TaxId"] == 30521
        assert record[1]["AkaTaxId"] == 0
        assert record[1]["Genus"] == ""
        assert record[1]["Species"] == ""
        assert record[1]["Subsp"] == ""
        assert record[1]["ModificationDate"] == "2020/04/29 00:00"

    def test_wrong(self):
        """Test parsing XML returned by ESummary with incorrect arguments."""
        # To create the XML file, use
        # >>> Bio.Entrez.esummary()
        with open("Entrez/esummary8.xml", "rb") as stream:
            with pytest.raises(RuntimeError):
                Entrez.read(stream)
        with open("Entrez/esummary8.xml", "rb") as stream:
            record = Entrez.read(stream, ignore_errors=True)
        assert len(record) == 1
        assert len(record.attributes) == 0
        assert record[0] == "Neither query_key nor id specified"
        assert record[0].tag == "ERROR"

    def test_integer_none(self):
        """Test parsing ESummary XML where an Integer is not defined."""
        # To create the XML file, use
        # >>> Entrez.esummary(db='pccompound', id='7488')
        with open("Entrez/esummary9.xml", "rb") as stream:
            records = Entrez.read(stream)
        assert len(records) == 1
        record = records[0]
        assert record["Id"] == "7488"
        assert record["CID"] == 7488
        assert record["SourceNameList"] == []
        assert record["SourceIDList"] == []
        assert len(record["SourceCategoryList"]) == 8
        assert record["SourceCategoryList"][0] == "Chemical Vendors"
        assert record["SourceCategoryList"][1] == "Research and Development"
        assert record["SourceCategoryList"][2] == "Curation Efforts"
        assert record["SourceCategoryList"][3] == "Governmental Organizations"
        assert record["SourceCategoryList"][4] == "Legacy Depositors"
        assert record["SourceCategoryList"][5] == "Subscription Services"
        assert record["SourceCategoryList"][6] == "Journal Publishers"
        assert record["SourceCategoryList"][7] == "NIH Initiatives"
        assert record["CreateDate"] == "2005/03/26 00:00"
        assert len(record["SynonymList"]) == 77
        assert record["SynonymList"][0] == "Terephthaloyl chloride"
        assert record["SynonymList"][1] == "100-20-9"
        assert record["SynonymList"][2] == "Terephthaloyl dichloride"
        assert record["SynonymList"][3] == "1,4-BENZENEDICARBONYL DICHLORIDE"
        assert record["SynonymList"][4] == "Terephthalic acid dichloride"
        assert record["SynonymList"][5] == "Terephthalic dichloride"
        assert record["SynonymList"][6] == "p-Phthaloyl chloride"
        assert record["SynonymList"][7] == "Terephthalic acid chloride"
        assert record["SynonymList"][8] == "p-Phthalyl dichloride"
        assert record["SynonymList"][9] == "p-Phthaloyl dichloride"
        assert record["SynonymList"][10] == "Terephthalyl dichloride"
        assert record["SynonymList"][11] == "1,4-Benzenedicarbonyl chloride"
        assert record["SynonymList"][12] == "p-Phenylenedicarbonyl dichloride"
        assert record["SynonymList"][13] == "benzene-1,4-dicarbonyl chloride"
        assert record["SynonymList"][14] == "NSC 41885"
        assert record["SynonymList"][15] == "terephthaloylchloride"
        assert record["SynonymList"][16] == "UNII-G247CO9608"
        assert record["SynonymList"][17] == "HSDB 5332"
        assert record["SynonymList"][18] == "EINECS 202-829-5"
        assert record["SynonymList"][19] == "BRN 0607796"
        assert record["SynonymList"][20] == "LXEJRKJRKIFVNY-UHFFFAOYSA-N"
        assert record["SynonymList"][21] == "MFCD00000693"
        assert record["SynonymList"][22] == "G247CO9608"
        assert record["SynonymList"][23] == "DSSTox_CID_6653"
        assert record["SynonymList"][24] == "DSSTox_RID_78175"
        assert record["SynonymList"][25] == "DSSTox_GSID_26653"
        assert record["SynonymList"][26] == "Q-201791"
        assert record["SynonymList"][27] == "Terephthaloyl chloride, 99+%"
        assert record["SynonymList"][28] == "CAS-100-20-9"
        assert record["SynonymList"][29] == "CCRIS 8626"
        assert record["SynonymList"][30] == "p-Phthalyl chloride"
        assert record["SynonymList"][31] == "terephthalic chloride"
        assert record["SynonymList"][32] == "tere-phthaloyl chloride"
        assert record["SynonymList"][33] == "AC1L1OVG"
        assert record["SynonymList"][34] == "ACMC-2097nf"
        assert record["SynonymList"][35] == "EC 202-829-5"
        assert record["SynonymList"][36] == "1,4-Dichloroformyl benzene"
        assert record["SynonymList"][37] == "SCHEMBL68148"
        assert record["SynonymList"][38] == "4-09-00-03318 (Beilstein Handbook Reference)"
        assert record["SynonymList"][39] == "KSC174E9T"
        assert record["SynonymList"][40] == "CHEMBL1893301"
        assert record["SynonymList"][41] == "DTXSID7026653"
        assert record["SynonymList"][42] == "KS-00000VAD"
        assert record["SynonymList"][43] == "benzene-1,4-dicarbonyl dichloride"
        assert record["SynonymList"][44] == "MolPort-003-926-079"
        assert record["SynonymList"][45] == "BCP27385"
        assert record["SynonymList"][46] == "NSC41885"
        assert record["SynonymList"][47] == "Tox21_201899"
        assert record["SynonymList"][48] == "Tox21_303166"
        assert record["SynonymList"][49] == "ANW-14185"
        assert record["SynonymList"][50] == "NSC-41885"
        assert record["SynonymList"][51] == "ZINC38141445"
        assert record["SynonymList"][52] == "AKOS015890038"
        assert record["SynonymList"][53] == "FCH1319904"
        assert record["SynonymList"][54] == "MCULE-9481285116"
        assert record["SynonymList"][55] == "RP25985"
        assert record["SynonymList"][56] == "Terephthaloyl chloride, >=99%, flakes"
        assert record["SynonymList"][57] == "NCGC00164045-01"
        assert record["SynonymList"][58] == "NCGC00164045-02"
        assert record["SynonymList"][59] == "NCGC00257127-01"
        assert record["SynonymList"][60] == "NCGC00259448-01"
        assert record["SynonymList"][61] == "AN-24545"
        assert record["SynonymList"][62] == "I764"
        assert record["SynonymList"][63] == "KB-10499"
        assert record["SynonymList"][64] == "OR315758"
        assert record["SynonymList"][65] == "SC-19185"
        assert record["SynonymList"][66] == "LS-148753"
        assert record["SynonymList"][67] == "RT-000669"
        assert record["SynonymList"][68] == "ST51037908"
        assert record["SynonymList"][69] == "6804-EP1441224A2"
        assert record["SynonymList"][70] == "I01-5090"
        assert record["SynonymList"][71] == "InChI=1/C8H4Cl2O2/c9-7(11)5-1-2-6(4-3-5)8(10)12/h1-4"
        assert record["SynonymList"][72] == "106158-15-0"
        assert record["SynonymList"][73] == "108454-76-8"
        assert record["SynonymList"][74] == "1640987-72-9"
        assert record["SynonymList"][75] == "188665-55-6"
        assert record["SynonymList"][76] == "1927884-58-9"
        assert len(record["MeSHHeadingList"]) == 1
        assert record["MeSHHeadingList"][0] == "terephthaloyl chloride"
        assert len(record["MeSHTermList"]) == 5
        assert record["MeSHTermList"][0] == "p-phthaloyl dichloride"
        assert record["MeSHTermList"][1] == "terephthaloyl dichloride"
        assert record["MeSHTermList"][2] == "1,4-benzenedicarbonyl dichloride"
        assert record["MeSHTermList"][3] == "1,4-phthaloyl dichloride"
        assert record["MeSHTermList"][4] == "terephthaloyl chloride"
        assert len(record["PharmActionList"]) == 0
        assert record["CommentList"] == []
        assert record["IUPACName"] == "benzene-1,4-dicarbonyl chloride"
        assert record["CanonicalSmiles"] == "C1=CC(=CC=C1C(=O)Cl)C(=O)Cl"
        assert record["IsomericSmiles"] == "C1=CC(=CC=C1C(=O)Cl)C(=O)Cl"
        assert record["RotatableBondCount"] == 2
        assert record["MolecularFormula"] == "C8H4Cl2O2"
        assert record["MolecularWeight"] == "203.018"
        assert record["TotalFormalCharge"] == 0
        assert record["XLogP"] == "4"
        assert record["HydrogenBondDonorCount"] == 0
        assert record["HydrogenBondAcceptorCount"] == 2
        assert record["Complexity"] == "173.000"
        assert record["HeavyAtomCount"] == 12
        assert record["AtomChiralCount"] == 0
        assert record["AtomChiralDefCount"] == 0
        assert record["AtomChiralUndefCount"] == 0
        assert record["BondChiralCount"] == 0
        assert record["BondChiralDefCount"] == 0
        assert record["BondChiralUndefCount"] == 0
        assert record["IsotopeAtomCount"] == 0
        assert record["CovalentUnitCount"] == 1
        assert record["TautomerCount"] == NoneElement(None, None, None)
        assert record["SubstanceIDList"] == []
        assert record["TPSA"] == "34.1"
        assert record["AssaySourceNameList"] == []
        assert record["MinAC"] == ""
        assert record["MaxAC"] == ""
        assert record["MinTC"] == ""
        assert record["MaxTC"] == ""
        assert record["ActiveAidCount"] == 1
        assert record["InactiveAidCount"] == NoneElement(None, None, None)
        assert record["TotalAidCount"] == 243
        assert record["InChIKey"] == "LXEJRKJRKIFVNY-UHFFFAOYSA-N"
        assert record["InChI"] == "InChI=1S/C8H4Cl2O2/c9-7(11)5-1-2-6(4-3-5)8(10)12/h1-4H"


class ELinkTest(unittest.TestCase):
    """Tests for parsing XML output returned by ELink."""

    def test_pubmed1(self):
        """Test parsing pubmed links returned by ELink (first test)."""
        # Retrieve IDs from PubMed for PMID 9298984 to the PubMed database
        # To create the XML file, use
        # >>> Bio.Entrez.elink(dbfrom="pubmed", id="9298984", cmd="neighbor")
        with open("Entrez/elink1.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert len(record) == 1
        assert len(record[0]) == 5
        assert record[0]["DbFrom"] == "pubmed"
        assert record[0]["IdList"] == ["9298984"]
        assert len(record[0]["LinkSetDb"]) == 7
        assert record[0]["LinkSetDb"][0]["DbTo"] == "pubmed"
        assert record[0]["LinkSetDb"][0]["LinkName"] == "pubmed_pubmed"
        assert len(record[0]["LinkSetDb"][0]["Link"]) == 101
        assert record[0]["LinkSetDb"][0]["Link"][0]["Id"] == "9298984"
        assert record[0]["LinkSetDb"][0]["Link"][1]["Id"] == "8794856"
        assert record[0]["LinkSetDb"][0]["Link"][2]["Id"] == "9700164"
        assert record[0]["LinkSetDb"][0]["Link"][3]["Id"] == "7914521"
        assert record[0]["LinkSetDb"][0]["Link"][4]["Id"] == "9914369"
        assert record[0]["LinkSetDb"][0]["Link"][5]["Id"] == "1339459"
        assert record[0]["LinkSetDb"][0]["Link"][6]["Id"] == "11590237"
        assert record[0]["LinkSetDb"][0]["Link"][7]["Id"] == "2211822"
        assert record[0]["LinkSetDb"][0]["Link"][8]["Id"] == "12686595"
        assert record[0]["LinkSetDb"][0]["Link"][9]["Id"] == "20980244"
        assert record[0]["LinkSetDb"][0]["Link"][10]["Id"] == "11146659"
        assert record[0]["LinkSetDb"][0]["Link"][11]["Id"] == "8978614"
        assert record[0]["LinkSetDb"][0]["Link"][12]["Id"] == "10893249"
        assert record[0]["LinkSetDb"][0]["Link"][13]["Id"] == "10402457"
        assert record[0]["LinkSetDb"][0]["Link"][14]["Id"] == "15371539"
        assert record[0]["LinkSetDb"][0]["Link"][15]["Id"] == "9074495"
        assert record[0]["LinkSetDb"][0]["Link"][16]["Id"] == "10806105"
        assert record[0]["LinkSetDb"][0]["Link"][17]["Id"] == "9490715"
        assert record[0]["LinkSetDb"][0]["Link"][18]["Id"] == "15915585"
        assert record[0]["LinkSetDb"][0]["Link"][19]["Id"] == "10545493"
        assert record[0]["LinkSetDb"][0]["Link"][20]["Id"] == "10523511"
        assert record[0]["LinkSetDb"][0]["Link"][21]["Id"] == "11483958"
        assert record[0]["LinkSetDb"][0]["Link"][22]["Id"] == "9869638"
        assert record[0]["LinkSetDb"][0]["Link"][23]["Id"] == "7690762"
        assert record[0]["LinkSetDb"][0]["Link"][24]["Id"] == "9425896"
        assert record[0]["LinkSetDb"][0]["Link"][25]["Id"] == "26892014"
        assert record[0]["LinkSetDb"][0]["Link"][26]["Id"] == "9378750"
        assert record[0]["LinkSetDb"][0]["Link"][27]["Id"] == "12515822"
        assert record[0]["LinkSetDb"][0]["Link"][28]["Id"] == "25919583"
        assert record[0]["LinkSetDb"][0]["Link"][29]["Id"] == "38830800"
        assert record[0]["LinkSetDb"][0]["Link"][30]["Id"] == "25081981"
        assert record[0]["LinkSetDb"][0]["Link"][31]["Id"] == "1691829"
        assert record[0]["LinkSetDb"][0]["Link"][32]["Id"] == "11146661"
        assert record[0]["LinkSetDb"][0]["Link"][33]["Id"] == "11685532"
        assert record[0]["LinkSetDb"][0]["Link"][34]["Id"] == "12080088"
        assert record[0]["LinkSetDb"][0]["Link"][35]["Id"] == "12034769"
        assert record[0]["LinkSetDb"][0]["Link"][36]["Id"] == "9852156"
        assert record[0]["LinkSetDb"][0]["Link"][37]["Id"] == "22733107"
        assert record[0]["LinkSetDb"][0]["Link"][38]["Id"] == "25194162"
        assert record[0]["LinkSetDb"][0]["Link"][39]["Id"] == "8923204"
        assert record[0]["LinkSetDb"][0]["Link"][40]["Id"] == "2022189"
        assert record[0]["LinkSetDb"][0]["Link"][41]["Id"] == "10985388"
        assert record[0]["LinkSetDb"][0]["Link"][42]["Id"] == "38402459"
        assert record[0]["LinkSetDb"][0]["Link"][43]["Id"] == "35609608"
        assert record[0]["LinkSetDb"][0]["Link"][44]["Id"] == "17222555"
        assert record[0]["LinkSetDb"][0]["Link"][45]["Id"] == "16741559"
        assert record[0]["LinkSetDb"][0]["Link"][46]["Id"] == "18936247"
        assert record[0]["LinkSetDb"][0]["Link"][47]["Id"] == "10749938"
        assert record[0]["LinkSetDb"][0]["Link"][48]["Id"] == "31150390"
        assert record[0]["LinkSetDb"][0]["Link"][49]["Id"] == "32794572"
        assert record[0]["LinkSetDb"][0]["Link"][50]["Id"] == "17895365"
        assert record[0]["LinkSetDb"][0]["Link"][51]["Id"] == "23300382"
        assert record[0]["LinkSetDb"][0]["Link"][52]["Id"] == "11914278"
        assert record[0]["LinkSetDb"][0]["Link"][53]["Id"] == "22563370"
        assert record[0]["LinkSetDb"][0]["Link"][54]["Id"] == "1541637"
        assert record[0]["LinkSetDb"][0]["Link"][55]["Id"] == "29706521"
        assert record[0]["LinkSetDb"][0]["Link"][56]["Id"] == "21118145"
        assert record[0]["LinkSetDb"][0]["Link"][57]["Id"] == "16732327"
        assert record[0]["LinkSetDb"][0]["Link"][58]["Id"] == "12388768"
        assert record[0]["LinkSetDb"][0]["Link"][59]["Id"] == "18202360"
        assert record[0]["LinkSetDb"][0]["Link"][60]["Id"] == "7585942"
        assert record[0]["LinkSetDb"][0]["Link"][61]["Id"] == "11179694"
        assert record[0]["LinkSetDb"][0]["Link"][62]["Id"] == "29158400"
        assert record[0]["LinkSetDb"][0]["Link"][63]["Id"] == "11352945"
        assert record[0]["LinkSetDb"][0]["Link"][64]["Id"] == "8056842"
        assert record[0]["LinkSetDb"][0]["Link"][65]["Id"] == "29046402"
        assert record[0]["LinkSetDb"][0]["Link"][66]["Id"] == "10398680"
        assert record[0]["LinkSetDb"][0]["Link"][67]["Id"] == "11267866"
        assert record[0]["LinkSetDb"][0]["Link"][68]["Id"] == "16516834"
        assert record[0]["LinkSetDb"][0]["Link"][69]["Id"] == "16839185"
        assert record[0]["LinkSetDb"][0]["Link"][70]["Id"] == "15616189"
        assert record[0]["LinkSetDb"][0]["Link"][71]["Id"] == "11266459"
        assert record[0]["LinkSetDb"][0]["Link"][72]["Id"] == "19641019"
        assert record[0]["LinkSetDb"][0]["Link"][73]["Id"] == "25976696"
        assert record[0]["LinkSetDb"][0]["Link"][74]["Id"] == "31250100"
        assert record[0]["LinkSetDb"][0]["Link"][75]["Id"] == "30784092"
        assert record[0]["LinkSetDb"][0]["Link"][76]["Id"] == "38019881"
        assert record[0]["LinkSetDb"][0]["Link"][77]["Id"] == "17182852"
        assert record[0]["LinkSetDb"][0]["Link"][78]["Id"] == "2211824"
        assert record[0]["LinkSetDb"][0]["Link"][79]["Id"] == "14522947"
        assert record[0]["LinkSetDb"][0]["Link"][80]["Id"] == "15268859"
        assert record[0]["LinkSetDb"][0]["Link"][81]["Id"] == "11252055"
        assert record[0]["LinkSetDb"][0]["Link"][82]["Id"] == "8175879"
        assert record[0]["LinkSetDb"][0]["Link"][83]["Id"] == "11102811"
        assert record[0]["LinkSetDb"][0]["Link"][84]["Id"] == "7904902"
        assert record[0]["LinkSetDb"][0]["Link"][85]["Id"] == "9606208"
        assert record[0]["LinkSetDb"][0]["Link"][86]["Id"] == "18460473"
        assert record[0]["LinkSetDb"][0]["Link"][87]["Id"] == "36920098"
        assert record[0]["LinkSetDb"][0]["Link"][88]["Id"] == "16510521"
        assert record[0]["LinkSetDb"][0]["Link"][89]["Id"] == "11092768"
        assert record[0]["LinkSetDb"][0]["Link"][90]["Id"] == "15824131"
        assert record[0]["LinkSetDb"][0]["Link"][91]["Id"] == "12235289"
        assert record[0]["LinkSetDb"][0]["Link"][92]["Id"] == "11266451"
        assert record[0]["LinkSetDb"][0]["Link"][93]["Id"] == "15485811"
        assert record[0]["LinkSetDb"][0]["Link"][94]["Id"] == "10898791"
        assert record[0]["LinkSetDb"][0]["Link"][95]["Id"] == "20729837"
        assert record[0]["LinkSetDb"][0]["Link"][96]["Id"] == "8548823"
        assert record[0]["LinkSetDb"][0]["Link"][97]["Id"] == "6807996"
        assert record[0]["LinkSetDb"][0]["Link"][98]["Id"] == "6791901"
        assert record[0]["LinkSetDb"][0]["Link"][99]["Id"] == "11715021"
        assert record[0]["LinkSetDb"][0]["Link"][100]["Id"] == "17333235"
        assert len(record[0]["LinkSetDb"][1]["Link"]) == 39
        assert record[0]["LinkSetDb"][1]["Link"][0]["Id"] == "38830800"
        assert record[0]["LinkSetDb"][1]["Link"][1]["Id"] == "38188366"
        assert record[0]["LinkSetDb"][1]["Link"][2]["Id"] == "37424454"
        assert record[0]["LinkSetDb"][1]["Link"][3]["Id"] == "34205694"
        assert record[0]["LinkSetDb"][1]["Link"][4]["Id"] == "32052088"
        assert record[0]["LinkSetDb"][1]["Link"][5]["Id"] == "29915359"
        assert record[0]["LinkSetDb"][1]["Link"][6]["Id"] == "29475948"
        assert record[0]["LinkSetDb"][1]["Link"][7]["Id"] == "29423089"
        assert record[0]["LinkSetDb"][1]["Link"][8]["Id"] == "29192061"
        assert record[0]["LinkSetDb"][1]["Link"][9]["Id"] == "28320824"
        assert record[0]["LinkSetDb"][1]["Link"][10]["Id"] == "28125061"
        assert record[0]["LinkSetDb"][1]["Link"][11]["Id"] == "20439434"
        assert record[0]["LinkSetDb"][1]["Link"][12]["Id"] == "19273145"
        assert record[0]["LinkSetDb"][1]["Link"][13]["Id"] == "19177000"
        assert record[0]["LinkSetDb"][1]["Link"][14]["Id"] == "18936247"
        assert record[0]["LinkSetDb"][1]["Link"][15]["Id"] == "18268100"
        assert record[0]["LinkSetDb"][1]["Link"][16]["Id"] == "17699596"
        assert record[0]["LinkSetDb"][1]["Link"][17]["Id"] == "16563186"
        assert record[0]["LinkSetDb"][1]["Link"][18]["Id"] == "16505164"
        assert record[0]["LinkSetDb"][1]["Link"][19]["Id"] == "16107559"
        assert record[0]["LinkSetDb"][1]["Link"][20]["Id"] == "15824131"
        assert record[0]["LinkSetDb"][1]["Link"][21]["Id"] == "15289669"
        assert record[0]["LinkSetDb"][1]["Link"][22]["Id"] == "15029241"
        assert record[0]["LinkSetDb"][1]["Link"][23]["Id"] == "12906131"
        assert record[0]["LinkSetDb"][1]["Link"][24]["Id"] == "12686595"
        assert record[0]["LinkSetDb"][1]["Link"][25]["Id"] == "12498345"
        assert record[0]["LinkSetDb"][1]["Link"][26]["Id"] == "11756470"
        assert record[0]["LinkSetDb"][1]["Link"][27]["Id"] == "11553716"
        assert record[0]["LinkSetDb"][1]["Link"][28]["Id"] == "11500386"
        assert record[0]["LinkSetDb"][1]["Link"][29]["Id"] == "11402076"
        assert record[0]["LinkSetDb"][1]["Link"][30]["Id"] == "11331754"
        assert record[0]["LinkSetDb"][1]["Link"][31]["Id"] == "10780705"
        assert record[0]["LinkSetDb"][1]["Link"][32]["Id"] == "10545493"
        assert record[0]["LinkSetDb"][1]["Link"][33]["Id"] == "10402457"
        assert record[0]["LinkSetDb"][1]["Link"][34]["Id"] == "10402425"
        assert record[0]["LinkSetDb"][1]["Link"][35]["Id"] == "9914368"
        assert record[0]["LinkSetDb"][1]["Link"][36]["Id"] == "9763420"
        assert record[0]["LinkSetDb"][1]["Link"][37]["Id"] == "9700166"
        assert record[0]["LinkSetDb"][1]["Link"][38]["Id"] == "9700164"
        assert len(record[0]["LinkSetDb"][2]["Link"]) == 5
        assert record[0]["LinkSetDb"][2]["Link"][0]["Id"] == "9298984"
        assert record[0]["LinkSetDb"][2]["Link"][1]["Id"] == "8794856"
        assert record[0]["LinkSetDb"][2]["Link"][2]["Id"] == "9700164"
        assert record[0]["LinkSetDb"][2]["Link"][3]["Id"] == "7914521"
        assert record[0]["LinkSetDb"][2]["Link"][4]["Id"] == "38830800"
        assert len(record[0]["LinkSetDb"][3]["Link"]) == 5
        assert record[0]["LinkSetDb"][3]["Link"][0]["Id"] == "9298984"
        assert record[0]["LinkSetDb"][3]["Link"][1]["Id"] == "8794856"
        assert record[0]["LinkSetDb"][3]["Link"][2]["Id"] == "9700164"
        assert record[0]["LinkSetDb"][3]["Link"][3]["Id"] == "7914521"
        assert record[0]["LinkSetDb"][3]["Link"][4]["Id"] == "9914369"
        assert len(record[0]["LinkSetDb"][4]["Link"]) == 56
        assert record[0]["LinkSetDb"][4]["Link"][0]["Id"] == "14732139"
        assert record[0]["LinkSetDb"][4]["Link"][1]["Id"] == "8909532"
        assert record[0]["LinkSetDb"][4]["Link"][2]["Id"] == "8898221"
        assert record[0]["LinkSetDb"][4]["Link"][3]["Id"] == "8824189"
        assert record[0]["LinkSetDb"][4]["Link"][4]["Id"] == "8824188"
        assert record[0]["LinkSetDb"][4]["Link"][5]["Id"] == "8794856"
        assert record[0]["LinkSetDb"][4]["Link"][6]["Id"] == "8763498"
        assert record[0]["LinkSetDb"][4]["Link"][7]["Id"] == "8706132"
        assert record[0]["LinkSetDb"][4]["Link"][8]["Id"] == "8706131"
        assert record[0]["LinkSetDb"][4]["Link"][9]["Id"] == "8647893"
        assert record[0]["LinkSetDb"][4]["Link"][10]["Id"] == "8617505"
        assert record[0]["LinkSetDb"][4]["Link"][11]["Id"] == "8560259"
        assert record[0]["LinkSetDb"][4]["Link"][12]["Id"] == "8521491"
        assert record[0]["LinkSetDb"][4]["Link"][13]["Id"] == "8505381"
        assert record[0]["LinkSetDb"][4]["Link"][14]["Id"] == "8485583"
        assert record[0]["LinkSetDb"][4]["Link"][15]["Id"] == "8416984"
        assert record[0]["LinkSetDb"][4]["Link"][16]["Id"] == "8267981"
        assert record[0]["LinkSetDb"][4]["Link"][17]["Id"] == "8143084"
        assert record[0]["LinkSetDb"][4]["Link"][18]["Id"] == "8023161"
        assert record[0]["LinkSetDb"][4]["Link"][19]["Id"] == "8005447"
        assert record[0]["LinkSetDb"][4]["Link"][20]["Id"] == "7914521"
        assert record[0]["LinkSetDb"][4]["Link"][21]["Id"] == "7906398"
        assert record[0]["LinkSetDb"][4]["Link"][22]["Id"] == "7860624"
        assert record[0]["LinkSetDb"][4]["Link"][23]["Id"] == "7854443"
        assert record[0]["LinkSetDb"][4]["Link"][24]["Id"] == "7854422"
        assert record[0]["LinkSetDb"][4]["Link"][25]["Id"] == "7846151"
        assert record[0]["LinkSetDb"][4]["Link"][26]["Id"] == "7821090"
        assert record[0]["LinkSetDb"][4]["Link"][27]["Id"] == "7758115"
        assert record[0]["LinkSetDb"][4]["Link"][28]["Id"] == "7739381"
        assert record[0]["LinkSetDb"][4]["Link"][29]["Id"] == "7704412"
        assert record[0]["LinkSetDb"][4]["Link"][30]["Id"] == "7698647"
        assert record[0]["LinkSetDb"][4]["Link"][31]["Id"] == "7664339"
        assert record[0]["LinkSetDb"][4]["Link"][32]["Id"] == "7642709"
        assert record[0]["LinkSetDb"][4]["Link"][33]["Id"] == "7642708"
        assert record[0]["LinkSetDb"][4]["Link"][34]["Id"] == "7579695"
        assert record[0]["LinkSetDb"][4]["Link"][35]["Id"] == "7542657"
        assert record[0]["LinkSetDb"][4]["Link"][36]["Id"] == "7502067"
        assert record[0]["LinkSetDb"][4]["Link"][37]["Id"] == "7172865"
        assert record[0]["LinkSetDb"][4]["Link"][38]["Id"] == "6966403"
        assert record[0]["LinkSetDb"][4]["Link"][39]["Id"] == "6793236"
        assert record[0]["LinkSetDb"][4]["Link"][40]["Id"] == "6684600"
        assert record[0]["LinkSetDb"][4]["Link"][41]["Id"] == "3928429"
        assert record[0]["LinkSetDb"][4]["Link"][42]["Id"] == "3670292"
        assert record[0]["LinkSetDb"][4]["Link"][43]["Id"] == "2686123"
        assert record[0]["LinkSetDb"][4]["Link"][44]["Id"] == "2683077"
        assert record[0]["LinkSetDb"][4]["Link"][45]["Id"] == "2512302"
        assert record[0]["LinkSetDb"][4]["Link"][46]["Id"] == "2498337"
        assert record[0]["LinkSetDb"][4]["Link"][47]["Id"] == "2195725"
        assert record[0]["LinkSetDb"][4]["Link"][48]["Id"] == "2185478"
        assert record[0]["LinkSetDb"][4]["Link"][49]["Id"] == "2139718"
        assert record[0]["LinkSetDb"][4]["Link"][50]["Id"] == "2139717"
        assert record[0]["LinkSetDb"][4]["Link"][51]["Id"] == "2022189"
        assert record[0]["LinkSetDb"][4]["Link"][52]["Id"] == "1999466"
        assert record[0]["LinkSetDb"][4]["Link"][53]["Id"] == "1684022"
        assert record[0]["LinkSetDb"][4]["Link"][54]["Id"] == "1406971"
        assert record[0]["LinkSetDb"][4]["Link"][55]["Id"] == "1339459"
        assert len(record[0]["LinkSetDb"][5]["Link"]) == 1
        assert record[0]["LinkSetDb"][5]["Link"][0]["Id"] == "9298984"
        assert len(record[0]["LinkSetDb"][6]["Link"]) == 5
        assert record[0]["LinkSetDb"][6]["Link"][0]["Id"] == "9298984"
        assert record[0]["LinkSetDb"][6]["Link"][1]["Id"] == "8794856"
        assert record[0]["LinkSetDb"][6]["Link"][2]["Id"] == "9700164"
        assert record[0]["LinkSetDb"][6]["Link"][3]["Id"] == "7914521"
        assert record[0]["LinkSetDb"][6]["Link"][4]["Id"] == "9914369"

    def test_nucleotide(self):
        """Test parsing Nucleotide to Protein links returned by ELink."""
        # Retrieve IDs from Nucleotide for GI  48819, 7140345 to Protein
        # To create the XML file, use
        # >>> Bio.Entrez.elink(dbfrom="nucleotide", db="protein",
        #                      id="48819,7140345")
        with open("Entrez/elink2.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert len(record) == 1
        assert len(record[0]) == 5
        assert record[0]["DbFrom"] == "nuccore"
        assert record[0]["IdList"] == ["48819", "7140345"]
        assert len(record[0]["LinkSetDb"]) == 1
        assert len(record[0]["LinkSetDb"][0]) == 3
        assert record[0]["LinkSetDb"][0]["DbTo"] == "protein"
        assert record[0]["LinkSetDb"][0]["LinkName"] == "nuccore_protein"
        assert len(record[0]["LinkSetDb"][0]["Link"]) == 1
        assert record[0]["LinkSetDb"][0]["Link"][0]["Id"] == "48820"

    def test_pubmed2(self):
        """Test parsing pubmed links returned by ELink (second test)."""
        # Retrieve PubMed related articles for PMIDs 11812492 11774222
        # with a publication date from 1995 to the present
        # To create the XML file, use
        # >>> Bio.Entrez.elink(dbfrom="pubmed", id="11812492,11774222",
        #                      db="pubmed", mindate="1995", datetype="pdat")
        with open("Entrez/elink3.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert len(record) == 1
        assert record[0]["DbFrom"] == "pubmed"
        assert len(record[0]["IdList"]) == 2
        assert record[0]["IdList"][0] == "11812492"
        assert record[0]["IdList"][1] == "11774222"
        assert record[0]["LinkSetDb"][0]["DbTo"] == "pubmed"
        assert record[0]["LinkSetDb"][0]["LinkName"] == "pubmed_pubmed"
        assert len(record[0]["LinkSetDb"]) == 6
        assert len(record[0]["LinkSetDb"][0]["Link"]) == 284
        assert record[0]["LinkSetDb"][0]["Link"][0]["Id"] == "39386366"
        assert record[0]["LinkSetDb"][0]["Link"][1]["Id"] == "39338906"
        assert record[0]["LinkSetDb"][0]["Link"][2]["Id"] == "39106844"
        assert record[0]["LinkSetDb"][0]["Link"][281]["Id"] == "10092480"
        assert record[0]["LinkSetDb"][0]["Link"][282]["Id"] == "9830540"
        assert record[0]["LinkSetDb"][0]["Link"][283]["Id"] == "9274032"
        assert len(record[0]["LinkSetDb"][1]["Link"]) == 3
        assert record[0]["LinkSetDb"][1]["Link"][0]["Id"] == "18508935"
        assert record[0]["LinkSetDb"][1]["Link"][1]["Id"] == "15780005"
        assert record[0]["LinkSetDb"][1]["Link"][2]["Id"] == "15024419"
        assert len(record[0]["LinkSetDb"][2]["Link"]) == 10
        assert record[0]["LinkSetDb"][2]["Link"][0]["Id"] == "16005284"
        assert record[0]["LinkSetDb"][2]["Link"][1]["Id"] == "15780005"
        assert record[0]["LinkSetDb"][2]["Link"][2]["Id"] == "15111095"
        assert record[0]["LinkSetDb"][2]["Link"][7]["Id"] == "11668631"
        assert record[0]["LinkSetDb"][2]["Link"][8]["Id"] == "10731564"
        assert record[0]["LinkSetDb"][2]["Link"][9]["Id"] == "10612825"
        assert len(record[0]["LinkSetDb"][3]["Link"]) == 10
        assert record[0]["LinkSetDb"][3]["Link"][0]["Id"] == "28358880"
        assert record[0]["LinkSetDb"][3]["Link"][1]["Id"] == "24053607"
        assert record[0]["LinkSetDb"][3]["Link"][2]["Id"] == "15780005"
        assert record[0]["LinkSetDb"][3]["Link"][7]["Id"] == "11668631"
        assert record[0]["LinkSetDb"][3]["Link"][8]["Id"] == "10731564"
        assert record[0]["LinkSetDb"][3]["Link"][9]["Id"] == "10612825"
        assert len(record[0]["LinkSetDb"][4]["Link"]) == 282
        assert record[0]["LinkSetDb"][4]["Link"][0]["Id"] == "40503031"
        assert record[0]["LinkSetDb"][4]["Link"][1]["Id"] == "39530225"
        assert record[0]["LinkSetDb"][4]["Link"][2]["Id"] == "39386366"
        assert record[0]["LinkSetDb"][4]["Link"][279]["Id"] == "10092480"
        assert record[0]["LinkSetDb"][4]["Link"][280]["Id"] == "9830540"
        assert record[0]["LinkSetDb"][4]["Link"][281]["Id"] == "9274032"
        assert len(record[0]["LinkSetDb"][5]["Link"]) == 10
        assert record[0]["LinkSetDb"][5]["Link"][0]["Id"] == "28358880"
        assert record[0]["LinkSetDb"][5]["Link"][1]["Id"] == "24053607"
        assert record[0]["LinkSetDb"][5]["Link"][2]["Id"] == "15780005"
        assert record[0]["LinkSetDb"][5]["Link"][7]["Id"] == "11668631"
        assert record[0]["LinkSetDb"][5]["Link"][8]["Id"] == "10731564"
        assert record[0]["LinkSetDb"][5]["Link"][9]["Id"] == "10612825"

    def test_medline(self):
        """Test parsing medline indexed articles returned by ELink."""
        # Retrieve MEDLINE indexed only related articles for PMID 12242737
        # To create the XML file, use
        # >>> Bio.Entrez.elink(dbfrom="pubmed", id="12242737", db="pubmed",
        #                      term="medline[sb]")
        with open("Entrez/elink4.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert len(record) == 1
        assert record[0]["DbFrom"] == "pubmed"
        assert record[0]["IdList"] == ["12242737"]
        assert len(record[0]["LinkSetDb"]) == 6
        assert record[0]["LinkSetDb"][0]["DbTo"] == "pubmed"
        assert record[0]["LinkSetDb"][0]["Link"][0]["Id"] == "38997184"
        assert record[0]["LinkSetDb"][0]["Link"][1]["Id"] == "37474462"
        assert record[0]["LinkSetDb"][0]["Link"][2]["Id"] == "36642882"
        assert record[0]["LinkSetDb"][0]["Link"][3]["Id"] == "33097552"
        assert record[0]["LinkSetDb"][0]["Link"][4]["Id"] == "32874413"
        assert record[0]["LinkSetDb"][0]["Link"][5]["Id"] == "32345945"
        assert record[0]["LinkSetDb"][0]["Link"][6]["Id"] == "31558954"
        assert record[0]["LinkSetDb"][0]["Link"][7]["Id"] == "31415410"
        assert record[0]["LinkSetDb"][0]["Link"][8]["Id"] == "31011388"
        assert record[0]["LinkSetDb"][0]["Link"][9]["Id"] == "30255446"
        assert record[0]["LinkSetDb"][0]["Link"][10]["Id"] == "30255443"
        assert record[0]["LinkSetDb"][0]["Link"][11]["Id"] == "30115443"
        assert record[0]["LinkSetDb"][0]["Link"][12]["Id"] == "29340557"
        assert record[0]["LinkSetDb"][0]["Link"][13]["Id"] == "29331013"
        assert record[0]["LinkSetDb"][0]["Link"][14]["Id"] == "29233545"
        assert record[0]["LinkSetDb"][0]["Link"][15]["Id"] == "29115651"
        assert record[0]["LinkSetDb"][0]["Link"][16]["Id"] == "29022115"
        assert record[0]["LinkSetDb"][0]["Link"][17]["Id"] == "27597304"
        assert record[0]["LinkSetDb"][0]["Link"][18]["Id"] == "27315096"
        assert record[0]["LinkSetDb"][0]["Link"][19]["Id"] == "27304929"
        assert record[0]["LinkSetDb"][0]["Link"][20]["Id"] == "27142382"
        assert record[0]["LinkSetDb"][0]["Link"][21]["Id"] == "26965844"
        assert record[0]["LinkSetDb"][0]["Link"][22]["Id"] == "26481976"
        assert record[0]["LinkSetDb"][0]["Link"][23]["Id"] == "26427946"
        assert record[0]["LinkSetDb"][0]["Link"][24]["Id"] == "26331169"
        assert record[0]["LinkSetDb"][0]["Link"][25]["Id"] == "25690945"
        assert record[0]["LinkSetDb"][0]["Link"][26]["Id"] == "25669229"
        assert record[0]["LinkSetDb"][0]["Link"][27]["Id"] == "25646204"
        assert record[0]["LinkSetDb"][0]["Link"][28]["Id"] == "25584961"
        assert record[0]["LinkSetDb"][0]["Link"][29]["Id"] == "25555004"
        assert record[0]["LinkSetDb"][0]["Link"][30]["Id"] == "25470877"
        assert record[0]["LinkSetDb"][0]["Link"][31]["Id"] == "25223134"
        assert record[0]["LinkSetDb"][0]["Link"][32]["Id"] == "25167349"
        assert record[0]["LinkSetDb"][0]["Link"][33]["Id"] == "24934824"
        assert record[0]["LinkSetDb"][0]["Link"][34]["Id"] == "24879722"
        assert record[0]["LinkSetDb"][0]["Link"][35]["Id"] == "24836494"
        assert record[0]["LinkSetDb"][0]["Link"][36]["Id"] == "24309417"
        assert record[0]["LinkSetDb"][0]["Link"][37]["Id"] == "23978699"
        assert record[0]["LinkSetDb"][0]["Link"][38]["Id"] == "23759294"
        assert record[0]["LinkSetDb"][0]["Link"][39]["Id"] == "23570763"
        assert record[0]["LinkSetDb"][0]["Link"][40]["Id"] == "23255877"
        assert record[0]["LinkSetDb"][0]["Link"][41]["Id"] == "22688104"
        assert record[0]["LinkSetDb"][0]["Link"][42]["Id"] == "22661362"
        assert record[0]["LinkSetDb"][0]["Link"][43]["Id"] == "22648258"
        assert record[0]["LinkSetDb"][0]["Link"][44]["Id"] == "22521021"
        assert record[0]["LinkSetDb"][0]["Link"][45]["Id"] == "22424988"
        assert record[0]["LinkSetDb"][0]["Link"][46]["Id"] == "22369817"
        assert record[0]["LinkSetDb"][0]["Link"][47]["Id"] == "22368911"
        assert record[0]["LinkSetDb"][0]["Link"][48]["Id"] == "22194507"
        assert record[0]["LinkSetDb"][0]["Link"][49]["Id"] == "22156652"
        assert record[0]["LinkSetDb"][0]["Link"][50]["Id"] == "22109321"
        assert record[0]["LinkSetDb"][0]["Link"][51]["Id"] == "21984464"
        assert record[0]["LinkSetDb"][0]["Link"][52]["Id"] == "21944608"
        assert record[0]["LinkSetDb"][0]["Link"][53]["Id"] == "21908142"
        assert record[0]["LinkSetDb"][0]["Link"][54]["Id"] == "21715237"
        assert record[0]["LinkSetDb"][0]["Link"][55]["Id"] == "21153952"
        assert record[0]["LinkSetDb"][0]["Link"][56]["Id"] == "20860230"
        assert record[0]["LinkSetDb"][0]["Link"][57]["Id"] == "20718377"
        assert record[0]["LinkSetDb"][0]["Link"][58]["Id"] == "20674629"
        assert record[0]["LinkSetDb"][0]["Link"][59]["Id"] == "20558858"
        assert record[0]["LinkSetDb"][0]["Link"][60]["Id"] == "20533237"
        assert record[0]["LinkSetDb"][0]["Link"][61]["Id"] == "20016426"
        assert record[0]["LinkSetDb"][0]["Link"][62]["Id"] == "19843737"
        assert record[0]["LinkSetDb"][0]["Link"][63]["Id"] == "19616724"
        assert record[0]["LinkSetDb"][0]["Link"][64]["Id"] == "19520357"
        assert record[0]["LinkSetDb"][0]["Link"][65]["Id"] == "18783095"
        assert record[0]["LinkSetDb"][0]["Link"][66]["Id"] == "18582671"
        assert record[0]["LinkSetDb"][0]["Link"][67]["Id"] == "18554854"
        assert record[0]["LinkSetDb"][0]["Link"][68]["Id"] == "18053822"
        assert record[0]["LinkSetDb"][0]["Link"][69]["Id"] == "18021675"
        assert record[0]["LinkSetDb"][0]["Link"][70]["Id"] == "17875143"
        assert record[0]["LinkSetDb"][0]["Link"][71]["Id"] == "17875142"
        assert record[0]["LinkSetDb"][0]["Link"][72]["Id"] == "17879696"
        assert record[0]["LinkSetDb"][0]["Link"][73]["Id"] == "17602359"
        assert record[0]["LinkSetDb"][0]["Link"][74]["Id"] == "17601500"
        assert record[0]["LinkSetDb"][0]["Link"][75]["Id"] == "17376366"
        assert record[0]["LinkSetDb"][0]["Link"][76]["Id"] == "17354190"
        assert record[0]["LinkSetDb"][0]["Link"][77]["Id"] == "17325998"
        assert record[0]["LinkSetDb"][0]["Link"][78]["Id"] == "17243036"
        assert record[0]["LinkSetDb"][0]["Link"][79]["Id"] == "17205643"
        assert record[0]["LinkSetDb"][0]["Link"][80]["Id"] == "17193860"
        assert record[0]["LinkSetDb"][0]["Link"][81]["Id"] == "17174054"
        assert record[0]["LinkSetDb"][0]["Link"][82]["Id"] == "17040637"
        assert record[0]["LinkSetDb"][0]["Link"][83]["Id"] == "16999328"
        assert record[0]["LinkSetDb"][0]["Link"][84]["Id"] == "16988291"
        assert record[0]["LinkSetDb"][0]["Link"][85]["Id"] == "16580806"
        assert record[0]["LinkSetDb"][0]["Link"][86]["Id"] == "16566645"
        assert record[0]["LinkSetDb"][0]["Link"][87]["Id"] == "16552382"
        assert record[0]["LinkSetDb"][0]["Link"][88]["Id"] == "16357381"
        assert record[0]["LinkSetDb"][0]["Link"][89]["Id"] == "16284132"
        assert record[0]["LinkSetDb"][0]["Link"][90]["Id"] == "16133609"
        assert record[0]["LinkSetDb"][0]["Link"][91]["Id"] == "16096604"
        assert record[0]["LinkSetDb"][0]["Link"][92]["Id"] == "16046437"
        assert record[0]["LinkSetDb"][0]["Link"][93]["Id"] == "15835031"
        assert record[0]["LinkSetDb"][0]["Link"][94]["Id"] == "15788585"
        assert record[0]["LinkSetDb"][0]["Link"][95]["Id"] == "15788584"
        assert record[0]["LinkSetDb"][0]["Link"][96]["Id"] == "15505294"
        assert record[0]["LinkSetDb"][0]["Link"][97]["Id"] == "15278705"
        assert record[0]["LinkSetDb"][0]["Link"][98]["Id"] == "15236131"
        assert record[0]["LinkSetDb"][0]["Link"][99]["Id"] == "15143223"
        assert record[0]["LinkSetDb"][0]["Link"][100]["Id"] == "15141648"
        assert record[0]["LinkSetDb"][0]["Link"][101]["Id"] == "15136027"
        assert record[0]["LinkSetDb"][0]["Link"][102]["Id"] == "15094630"
        assert record[0]["LinkSetDb"][0]["Link"][103]["Id"] == "15022983"
        assert record[0]["LinkSetDb"][0]["Link"][104]["Id"] == "14661668"
        assert record[0]["LinkSetDb"][0]["Link"][105]["Id"] == "14661661"
        assert record[0]["LinkSetDb"][0]["Link"][106]["Id"] == "14650118"
        assert record[0]["LinkSetDb"][0]["Link"][107]["Id"] == "12878072"
        assert record[0]["LinkSetDb"][0]["Link"][108]["Id"] == "12846253"
        assert record[0]["LinkSetDb"][0]["Link"][109]["Id"] == "12822521"
        assert record[0]["LinkSetDb"][0]["Link"][110]["Id"] == "12733684"
        assert record[0]["LinkSetDb"][0]["Link"][111]["Id"] == "12719915"
        assert record[0]["LinkSetDb"][0]["Link"][112]["Id"] == "12563154"
        assert record[0]["LinkSetDb"][0]["Link"][113]["Id"] == "12242737"
        assert record[0]["LinkSetDb"][0]["Link"][114]["Id"] == "12226761"
        assert record[0]["LinkSetDb"][0]["Link"][115]["Id"] == "12164574"
        assert record[0]["LinkSetDb"][0]["Link"][116]["Id"] == "12069469"
        assert record[0]["LinkSetDb"][0]["Link"][117]["Id"] == "11973040"
        assert record[0]["LinkSetDb"][0]["Link"][118]["Id"] == "11895298"
        assert record[0]["LinkSetDb"][0]["Link"][119]["Id"] == "11781922"
        assert record[0]["LinkSetDb"][0]["Link"][120]["Id"] == "11775722"
        assert record[0]["LinkSetDb"][0]["Link"][121]["Id"] == "11762248"
        assert record[0]["LinkSetDb"][0]["Link"][122]["Id"] == "11702119"
        assert record[0]["LinkSetDb"][0]["Link"][123]["Id"] == "11368937"
        assert record[0]["LinkSetDb"][0]["Link"][124]["Id"] == "11329656"
        assert record[0]["LinkSetDb"][0]["Link"][125]["Id"] == "11329655"
        assert record[0]["LinkSetDb"][0]["Link"][126]["Id"] == "11329162"
        assert record[0]["LinkSetDb"][0]["Link"][127]["Id"] == "11274884"
        assert record[0]["LinkSetDb"][0]["Link"][128]["Id"] == "11218011"
        assert record[0]["LinkSetDb"][0]["Link"][129]["Id"] == "11125632"
        assert record[0]["LinkSetDb"][0]["Link"][130]["Id"] == "11016058"
        assert record[0]["LinkSetDb"][0]["Link"][131]["Id"] == "10688063"
        assert record[0]["LinkSetDb"][0]["Link"][132]["Id"] == "10499696"
        assert record[0]["LinkSetDb"][0]["Link"][133]["Id"] == "10222515"
        assert record[0]["LinkSetDb"][0]["Link"][134]["Id"] == "10222514"
        assert record[0]["LinkSetDb"][0]["Link"][135]["Id"] == "10024396"
        assert record[0]["LinkSetDb"][0]["Link"][136]["Id"] == "9793138"
        assert record[0]["LinkSetDb"][0]["Link"][137]["Id"] == "9757294"
        assert record[0]["LinkSetDb"][0]["Link"][138]["Id"] == "9575723"
        assert record[0]["LinkSetDb"][0]["Link"][139]["Id"] == "9510579"
        assert record[0]["LinkSetDb"][0]["Link"][140]["Id"] == "9456947"
        assert record[0]["LinkSetDb"][0]["Link"][141]["Id"] == "9314960"
        assert record[0]["LinkSetDb"][0]["Link"][142]["Id"] == "9314959"
        assert record[0]["LinkSetDb"][0]["Link"][143]["Id"] == "9269670"
        assert record[0]["LinkSetDb"][0]["Link"][144]["Id"] == "9193407"
        assert record[0]["LinkSetDb"][0]["Link"][145]["Id"] == "8872409"
        assert record[0]["LinkSetDb"][0]["Link"][146]["Id"] == "8756148"
        assert record[0]["LinkSetDb"][0]["Link"][147]["Id"] == "8903064"
        assert record[0]["LinkSetDb"][0]["Link"][148]["Id"] == "8599783"
        assert record[0]["LinkSetDb"][0]["Link"][149]["Id"] == "8153333"
        assert record[0]["LinkSetDb"][0]["Link"][150]["Id"] == "1343378"
        assert record[0]["LinkSetDb"][0]["Link"][151]["Id"] == "1535863"
        assert record[0]["LinkSetDb"][0]["Link"][152]["Id"] == "4818442"
        assert record[0]["LinkSetDb"][0]["Link"][153]["Id"] == "4808999"
        assert record[0]["LinkSetDb"][0]["Link"][154]["Id"] == "13969511"
        assert record[0]["LinkSetDb"][0]["Link"][155]["Id"] == "13808134"
        assert record[0]["LinkSetDb"][0]["LinkName"] == "pubmed_pubmed"
        assert record[0]["LinkSetDb"][1]["DbTo"] == "pubmed"
        assert record[0]["LinkSetDb"][1]["Link"][0]["Id"] == "40900773"
        assert record[0]["LinkSetDb"][1]["Link"][1]["Id"] == "28779191"
        assert record[0]["LinkSetDb"][1]["Link"][2]["Id"] == "21102533"
        assert record[0]["LinkSetDb"][1]["Link"][3]["Id"] == "19517148"
        assert record[0]["LinkSetDb"][1]["Link"][4]["Id"] == "19132488"
        assert record[0]["LinkSetDb"][1]["Link"][5]["Id"] == "15278705"
        assert record[0]["LinkSetDb"][1]["LinkName"] == "pubmed_pubmed_citedin"
        assert record[0]["LinkSetDb"][2]["DbTo"] == "pubmed"
        assert record[0]["LinkSetDb"][2]["Link"][0]["Id"] == "20718377"
        assert record[0]["LinkSetDb"][2]["Link"][1]["Id"] == "12242737"
        assert record[0]["LinkSetDb"][2]["Link"][2]["Id"] == "11329656"
        assert record[0]["LinkSetDb"][2]["Link"][3]["Id"] == "11218011"
        assert record[0]["LinkSetDb"][2]["Link"][4]["Id"] == "9757294"
        assert record[0]["LinkSetDb"][2]["LinkName"] == "pubmed_pubmed_combined"
        assert record[0]["LinkSetDb"][3]["DbTo"] == "pubmed"
        assert record[0]["LinkSetDb"][3]["Link"][0]["Id"] == "20718377"
        assert record[0]["LinkSetDb"][3]["Link"][1]["Id"] == "17193860"
        assert record[0]["LinkSetDb"][3]["Link"][2]["Id"] == "12242737"
        assert record[0]["LinkSetDb"][3]["Link"][3]["Id"] == "11218011"
        assert record[0]["LinkSetDb"][3]["Link"][4]["Id"] == "9757294"
        assert record[0]["LinkSetDb"][3]["LinkName"] == "pubmed_pubmed_five"
        assert record[0]["LinkSetDb"][4]["DbTo"] == "pubmed"
        assert record[0]["LinkSetDb"][4]["Link"][0]["Id"] == "12242737"
        assert record[0]["LinkSetDb"][4]["LinkName"] == "pubmed_pubmed_reviews"
        assert record[0]["LinkSetDb"][5]["DbTo"] == "pubmed"
        assert record[0]["LinkSetDb"][5]["Link"][0]["Id"] == "20718377"
        assert record[0]["LinkSetDb"][5]["Link"][1]["Id"] == "17193860"
        assert record[0]["LinkSetDb"][5]["Link"][2]["Id"] == "12242737"
        assert record[0]["LinkSetDb"][5]["Link"][3]["Id"] == "11218011"
        assert record[0]["LinkSetDb"][5]["Link"][4]["Id"] == "9757294"
        assert record[0]["LinkSetDb"][5]["LinkName"] == "pubmed_pubmed_reviews_five"

    def test_pubmed3(self):
        """Test parsing pubmed link returned by ELink (third test)."""
        # Create a hyperlink to the first link available for PMID 10611131
        # in PubMed
        # To create the XML file, use
        # >>> Bio.Entrez.elink(dbfrom="pubmed", id="10611131", cmd="prlinks")

        with open("Entrez/elink5.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert len(record) == 1
        assert len(record[0]) == 5
        assert record[0]["DbFrom"] == "pubmed"
        assert len(record[0]["LinkSetDb"]) == 0
        assert len(record[0]["LinkSetDbHistory"]) == 0
        assert len(record[0]["ERROR"]) == 0
        assert len(record[0]["IdUrlList"]) == 2
        assert len(record[0]["IdUrlList"]["FirstChars"]) == 0
        assert len(record[0]["IdUrlList"]["IdUrlSet"]) == 1

        assert record[0]["IdUrlList"]["IdUrlSet"][0]["Id"] == "10611131"
        assert len(record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"]) == 1
        assert record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["Url"] == "https://academic.oup.com/brain/article-lookup/doi/10.1093/brain/123.1.171"
        assert record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["Url"].attributes == {"LNG": "EN"}
        assert record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["IconUrl"] == "/corehtml/query/egifs/https:--academic.oup.com-images-oup_pubmed.png"
        assert record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["IconUrl"].attributes == {"LNG": "EN"}
        assert record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["SubjectType"] == []
        assert len(record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["Attribute"]) == 1
        assert record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["Attribute"][0] == "subscription/membership/fee required"
        assert record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["Provider"]["Name"] == "Silverchair Information Systems"
        assert record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["Provider"]["NameAbbr"] == "silverchair"
        assert record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["Provider"]["Id"] == "7898"

    def test_pubmed4(self):
        """Test parsing pubmed links returned by ELink (fourth test)."""
        # List all available links in PubMed, except for libraries, for
        # PMIDs 12085856 and 12085853
        # To create the XML file, use
        # >>> Bio.Entrez.elink(dbfrom="pubmed", id="12085856,12085853", cmd="llinks")
        with open("Entrez/elink6.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert record[0]["DbFrom"] == "pubmed"
        assert len(record[0]["IdUrlList"]) == 2
        assert record[0]["IdUrlList"]["IdUrlSet"][0]["Id"] == "12085856"
        assert len(record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"]) == 1
        assert record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["Category"] == ["Medical"]
        assert record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["Url"] == "https://medlineplus.gov/coronaryarterybypasssurgery.html"
        assert record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["Attribute"] == ["free resource"]
        assert record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["SubjectType"] == []
        assert record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["IconUrl"] == "/corehtml/query/egifs/https:--medlineplus.gov-images-linkout_sm.gif"
        assert record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["Provider"]["Name"] == "MedlinePlus Health Information"
        assert record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["Provider"]["NameAbbr"] == "MEDPLUS"
        assert record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["Provider"]["Id"] == "3162"
        assert record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["LinkName"] == "Coronary Artery Bypass Surgery"
        assert len(record[0]["IdUrlList"]["IdUrlSet"][1]) == 2
        assert record[0]["IdUrlList"]["IdUrlSet"][1]["Id"] == "12085853"
        assert len(record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"]) == 2
        assert record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][0]["Category"] == ["Medical"]
        assert record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][0]["Url"] == "https://medlineplus.gov/exerciseandphysicalfitness.html"
        assert record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][0]["IconUrl"] == "/corehtml/query/egifs/https:--medlineplus.gov-images-linkout_sm.gif"
        assert record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][0]["Attribute"] == ["free resource"]
        assert record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][0]["SubjectType"] == []
        assert record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][0]["LinkName"] == "Exercise and Physical Fitness"
        assert record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][0]["Provider"]["Name"] == "MedlinePlus Health Information"
        assert record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][0]["Provider"]["NameAbbr"] == "MEDPLUS"
        assert record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][0]["Provider"]["Id"] == "3162"
        assert record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][1]["Category"] == ["Medical"]
        assert record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][1]["Attribute"] == ["free resource"]
        assert record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][1]["Url"] == "https://medlineplus.gov/exerciseandphysicalfitness.html"
        assert record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][1]["IconUrl"] == "/corehtml/query/egifs/https:--medlineplus.gov-images-linkout_sm.gif"
        assert record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][1]["LinkName"] == "Exercise and Physical Fitness"
        assert record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][1]["SubjectType"] == []
        assert record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][1]["Provider"]["Name"] == "MedlinePlus Consumer Health Information"
        assert record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][1]["Provider"]["NameAbbr"] == "medlineplus2"
        assert record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][1]["Provider"]["Id"] == "10405"

    def test_pubmed5(self):
        """Test parsing pubmed links returned by ELink (fifth test)."""
        # List Entrez database links for PubMed PMIDs 12169658 and 11748140
        # To create the XML file, use
        # >>> Bio.Entrez.elink(dbfrom="pubmed", id="12169658,11748140",
        #                      cmd="acheck")
        with open("Entrez/elink7.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert len(record) == 1
        assert record[0]["DbFrom"] == "pubmed"
        assert len(record[0]["IdCheckList"]) == 2
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["Id"] == "12169658"
        assert len(record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"]) == 16
        assert len(record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][0]) == 5
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][0]["DbTo"] == "books"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][0]["LinkName"] == "pubmed_books_refs"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][0]["MenuTag"] == "Cited in Books"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][0]["HtmlTag"] == "Cited in Books"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][0]["Priority"] == "185"
        assert len(record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][1]) == 5
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][1]["DbTo"] == "cdd"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][1]["LinkName"] == "pubmed_cdd"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][1]["MenuTag"] == "Conserved Domain Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][1]["HtmlTag"] == "Conserved Domains"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][1]["Priority"] == "130"
        assert len(record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][2]) == 5
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][2]["DbTo"] == "gene"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][2]["LinkName"] == "pubmed_gene"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][2]["MenuTag"] == "Gene Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][2]["HtmlTag"] == "Gene"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][2]["Priority"] == "128"
        assert len(record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][3]) == 5
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][3]["DbTo"] == "geoprofiles"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][3]["LinkName"] == "pubmed_geoprofiles"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][3]["MenuTag"] == "GEO Profile Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][3]["HtmlTag"] == "GEO Profiles"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][3]["Priority"] == "170"
        assert len(record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][4]) == 5
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][4]["DbTo"] == "nuccore"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][4]["LinkName"] == "pubmed_nuccore"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][4]["MenuTag"] == "Nucleotide Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][4]["HtmlTag"] == "Nucleotide"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][4]["Priority"] == "128"
        assert len(record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"]) == 15
        assert len(record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][0]) == 5
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][0]["DbTo"] == "books"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][0]["LinkName"] == "pubmed_books_refs"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][0]["MenuTag"] == "Cited in Books"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][0]["HtmlTag"] == "Cited in Books"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][0]["Priority"] == "185"
        assert len(record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][1]) == 5
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][1]["DbTo"] == "gene"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][1]["LinkName"] == "pubmed_gene"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][1]["MenuTag"] == "Gene Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][1]["HtmlTag"] == "Gene"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][1]["Priority"] == "128"
        assert len(record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][2]) == 5
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][2]["DbTo"] == "geoprofiles"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][2]["LinkName"] == "pubmed_geoprofiles"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][2]["MenuTag"] == "GEO Profile Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][2]["HtmlTag"] == "GEO Profiles"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][2]["Priority"] == "170"
        assert len(record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][3]) == 5
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][3]["DbTo"] == "nuccore"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][3]["LinkName"] == "pubmed_nuccore"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][3]["MenuTag"] == "Nucleotide Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][3]["HtmlTag"] == "Nucleotide"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][3]["Priority"] == "128"
        assert len(record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][4]) == 5
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][4]["DbTo"] == "nuccore"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][4]["LinkName"] == "pubmed_nuccore_refseq"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][4]["MenuTag"] == "Nucleotide (RefSeq) Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][4]["HtmlTag"] == "Nucleotide (RefSeq)"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][4]["Priority"] == "128"

    def test_pubmed6(self):
        """Test parsing pubmed links returned by ELink (sixth test)."""
        # Check for the existence of a Related Articles link for PMID
        # 12068369.
        # To create the XML file, use
        # >>> Bio.Entrez.elink(dbfrom="pubmed", id="12068369", cmd="ncheck")

        with open("Entrez/elink8.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert len(record) == 1
        assert record[0]["DbFrom"] == "pubmed"
        assert len(record[0]["IdCheckList"]) == 2
        assert len(record[0]["IdCheckList"]["Id"]) == 1
        assert record[0]["IdCheckList"]["Id"][0] == "12068369"
        assert len(record[0]["IdCheckList"]["Id"][0].attributes) == 1
        assert record[0]["IdCheckList"]["Id"][0].attributes["HasNeighbor"] == "Y"
        assert len(record[0]["IdCheckList"]["IdLinkSet"]) == 0


class ESpellTest(unittest.TestCase):
    """Tests for parsing XML output returned by ESpell."""

    def test_espell(self):
        """Test parsing XML output returned by ESpell."""
        # Request suggestions for the PubMed search biopythooon
        # To create the XML file, use
        # >>> Bio.Entrez.espell(db="pubmed", term="biopythooon")
        with open("Entrez/espell.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert record["Database"] == "pubmed"
        assert record["Query"] == "biopythooon"
        assert record["CorrectedQuery"] == "biopython"
        assert len(record["SpelledQuery"]) == 2
        assert record["SpelledQuery"][0] == ""
        assert record["SpelledQuery"][0].tag == "Original"
        assert record["SpelledQuery"][1] == "biopython"
        assert record["SpelledQuery"][1].tag == "Replaced"


class EFetchTest(unittest.TestCase):
    """Tests for parsing XML output returned by EFetch."""

    def test_pubmed1(self):
        """Test parsing XML returned by EFetch, PubMed database (first test)."""
        # In PubMed display PMIDs 12091962 and 9997 in xml retrieval mode
        # and abstract retrieval type.
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db='pubmed', id='12091962,9997',
        #                       retmode='xml', rettype='abstract')
        with open("Entrez/pubmed1.xml", "rb") as stream:
            record = Entrez.read(stream)
        # fmt: off
        assert record["PubmedBookArticle"] == []
        record = record["PubmedArticle"]
        assert record[0]["MedlineCitation"].attributes["Owner"] == "KIE"
        assert record[0]["MedlineCitation"].attributes["Status"] == "MEDLINE"
        assert record[0]["MedlineCitation"]["PMID"] == "12091962"
        assert record[0]["MedlineCitation"]["DateCompleted"]["Year"] == "1991"
        assert record[0]["MedlineCitation"]["DateCompleted"]["Month"] == "01"
        assert record[0]["MedlineCitation"]["DateCompleted"]["Day"] == "22"
        assert record[0]["MedlineCitation"]["DateRevised"]["Year"] == "2007"
        assert record[0]["MedlineCitation"]["DateRevised"]["Month"] == "11"
        assert record[0]["MedlineCitation"]["DateRevised"]["Day"] == "15"
        assert record[0]["MedlineCitation"]["Article"].attributes["PubModel"] == "Print"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["ISSN"] == "1043-1578"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["ISSN"].attributes[
                "IssnType"
            ] == "Print"
        assert record[0]["MedlineCitation"]["Article"]["Journal"][
                "JournalIssue"
            ].attributes["CitedMedium"] == "Print"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                "Volume"
            ] == "17"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["Issue"] == "1"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                "PubDate"
            ]["Year"] == "1990"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                "PubDate"
            ]["Season"] == "Spring"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["Title"] == "Social justice (San Francisco, Calif.)"
        assert record[0]["MedlineCitation"]["Article"]["ArticleTitle"] == "The treatment of AIDS behind the walls of correctional facilities."
        assert record[0]["MedlineCitation"]["Article"]["Pagination"]["MedlinePgn"] == "113-25"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"].attributes[
                "CompleteYN"
            ] == "Y"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][0].attributes[
                "ValidYN"
            ] == "Y"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][0]["LastName"] == "Olivero"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][0]["ForeName"] == "J Michael"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][0]["Initials"] == "JM"
        assert record[0]["MedlineCitation"]["Article"]["Language"] == ["eng"]
        assert record[0]["MedlineCitation"]["Article"]["PublicationTypeList"] == ["Journal Article", "Review"]
        assert record[0]["MedlineCitation"]["MedlineJournalInfo"]["Country"] == "United States"
        assert record[0]["MedlineCitation"]["MedlineJournalInfo"]["MedlineTA"] == "Soc Justice"
        assert record[0]["MedlineCitation"]["MedlineJournalInfo"]["NlmUniqueID"] == "9891830"
        assert record[0]["MedlineCitation"]["CitationSubset"] == []
        assert record[0]["MedlineCitation"]["MeshHeadingList"][0]["DescriptorName"] == "AIDS Serodiagnosis"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][0][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][1]["DescriptorName"] == "Acquired Immunodeficiency Syndrome"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][1][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "Y"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][2]["DescriptorName"] == "Civil Rights"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][2][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][3]["DescriptorName"] == "HIV Seropositivity"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][3][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "Y"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][4]["DescriptorName"] == "Humans"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][4][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][5]["DescriptorName"] == "Jurisprudence"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][5][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "Y"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][6]["DescriptorName"] == "Law Enforcement"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][6][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][7]["DescriptorName"] == "Mass Screening"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][7][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][8]["DescriptorName"] == "Minority Groups"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][8][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][9]["DescriptorName"] == "Organizational Policy"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][9][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][10]["DescriptorName"] == "Patient Care"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][10][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][11]["DescriptorName"] == "Prejudice"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][11][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][12]["DescriptorName"] == "Prisoners"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][12][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "Y"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][13]["DescriptorName"] == "Public Policy"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][13][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "Y"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][14]["DescriptorName"] == "Quarantine"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][14][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][15]["DescriptorName"] == "Social Control, Formal"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][15][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][16]["DescriptorName"] == "Statistics as Topic"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][16][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][17]["DescriptorName"] == "Stereotyping"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][17][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][18]["DescriptorName"] == "United States"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][18][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["NumberOfReferences"] == "63"
        assert record[0]["MedlineCitation"]["OtherID"][0] == "31840"
        assert record[0]["MedlineCitation"]["OtherID"][0].attributes["Source"] == "KIE"
        assert record[0]["MedlineCitation"]["KeywordList"][0].attributes["Owner"] == "KIE"
        assert record[0]["MedlineCitation"]["KeywordList"][0][0] == "Health Care and Public Health"
        assert record[0]["MedlineCitation"]["KeywordList"][0][0].attributes[
                "MajorTopicYN"
            ] == "N"
        assert record[0]["MedlineCitation"]["KeywordList"][0][1] == "Legal Approach"
        assert record[0]["MedlineCitation"]["KeywordList"][0][1].attributes[
                "MajorTopicYN"
            ] == "N"
        assert record[0]["MedlineCitation"]["GeneralNote"][0] == "14 fn."
        assert record[0]["MedlineCitation"]["GeneralNote"][0].attributes["Owner"] == "KIE"
        assert record[0]["MedlineCitation"]["GeneralNote"][1] == "KIE BoB Subject Heading: AIDS"
        assert record[0]["MedlineCitation"]["GeneralNote"][1].attributes["Owner"] == "KIE"
        assert record[0]["MedlineCitation"]["GeneralNote"][2] == "63 refs."
        assert record[0]["MedlineCitation"]["GeneralNote"][2].attributes["Owner"] == "KIE"
        assert record[0]["PubmedData"]["History"][0].attributes["PubStatus"] == "pubmed"
        assert record[0]["PubmedData"]["History"][0]["Year"] == "1990"
        assert record[0]["PubmedData"]["History"][0]["Month"] == "4"
        assert record[0]["PubmedData"]["History"][0]["Day"] == "1"
        assert record[0]["PubmedData"]["History"][0]["Hour"] == "0"
        assert record[0]["PubmedData"]["History"][0]["Minute"] == "0"
        assert record[0]["PubmedData"]["History"][1].attributes["PubStatus"] == "medline"
        assert record[0]["PubmedData"]["History"][1]["Year"] == "2002"
        assert record[0]["PubmedData"]["History"][1]["Month"] == "7"
        assert record[0]["PubmedData"]["History"][1]["Day"] == "16"
        assert record[0]["PubmedData"]["History"][1]["Hour"] == "10"
        assert record[0]["PubmedData"]["History"][1]["Minute"] == "1"
        assert record[0]["PubmedData"]["PublicationStatus"] == "ppublish"
        assert len(record[0]["PubmedData"]["ArticleIdList"]) == 1
        assert record[0]["PubmedData"]["ArticleIdList"][0] == "12091962"
        assert record[0]["PubmedData"]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert record[1]["MedlineCitation"].attributes["Owner"] == "NLM"
        assert record[1]["MedlineCitation"].attributes["Status"] == "MEDLINE"
        assert record[1]["MedlineCitation"]["PMID"] == "9997"
        assert record[1]["MedlineCitation"]["DateCompleted"]["Year"] == "1976"
        assert record[1]["MedlineCitation"]["DateCompleted"]["Month"] == "12"
        assert record[1]["MedlineCitation"]["DateCompleted"]["Day"] == "30"
        assert record[1]["MedlineCitation"]["DateRevised"]["Year"] == "2019"
        assert record[1]["MedlineCitation"]["DateRevised"]["Month"] == "06"
        assert record[1]["MedlineCitation"]["DateRevised"]["Day"] == "09"
        assert record[1]["MedlineCitation"]["Article"].attributes["PubModel"] == "Print"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["ISSN"] == "0006-3002"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["ISSN"].attributes["IssnType"] == "Print"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"].attributes["CitedMedium"] == "Print"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["Volume"] == "446"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["Issue"] == "1"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["Year"] == "1976"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["Month"] == "Sep"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["Day"] == "28"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["Title"] == "Biochimica et biophysica acta"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["ISOAbbreviation"] == "Biochim Biophys Acta"
        assert record[1]["MedlineCitation"]["Article"]["ArticleTitle"] == "Magnetic studies of Chromatium flavocytochrome C552. A mechanism for heme-flavin interaction."
        assert record[1]["MedlineCitation"]["Article"]["Pagination"]["MedlinePgn"] == "179-91"
        assert record[1]["MedlineCitation"]["Article"]["Abstract"]["AbstractText"] == ["Electron paramagnetic resonance and magnetic susceptibility studies of Chromatium flavocytochrome C552 and its diheme flavin-free subunit at temperatures below 45 degrees K are reported. The results show that in the intact protein and the subunit the two low-spin (S = 1/2) heme irons are distinguishable, giving rise to separate EPR signals. In the intact protein only, one of the heme irons exists in two different low spin environments in the pH range 5.5 to 10.5, while the other remains in a constant environment. Factors influencing the variable heme iron environment also influence flavin reactivity, indicating the existence of a mechanism for heme-flavin interaction."]
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"].attributes["CompleteYN"] == "Y"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][0].attributes["ValidYN"] == "Y"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][0]["LastName"] == "Strekas"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][0]["ForeName"] == "T C"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][0]["Initials"] == "TC"
        assert record[1]["MedlineCitation"]["Article"]["Language"] == ["eng"]
        assert record[1]["MedlineCitation"]["Article"]["PublicationTypeList"] == ["Journal Article"]
        assert record[1]["MedlineCitation"]["MedlineJournalInfo"]["Country"] == "Netherlands"
        assert record[1]["MedlineCitation"]["MedlineJournalInfo"]["MedlineTA"] == "Biochim Biophys Acta"
        assert record[1]["MedlineCitation"]["MedlineJournalInfo"]["NlmUniqueID"] == "0217513"
        assert record[1]["MedlineCitation"]["ChemicalList"][0]["RegistryNumber"] == "0"
        assert record[1]["MedlineCitation"]["ChemicalList"][0]["NameOfSubstance"] == "Cytochrome c Group"
        assert record[1]["MedlineCitation"]["ChemicalList"][1]["RegistryNumber"] == "0"
        assert record[1]["MedlineCitation"]["ChemicalList"][1]["NameOfSubstance"] == "Flavins"
        assert record[1]["MedlineCitation"]["ChemicalList"][2]["RegistryNumber"] == "42VZT0U6YR"
        assert record[1]["MedlineCitation"]["ChemicalList"][2]["NameOfSubstance"] == "Heme"
        assert record[1]["MedlineCitation"]["ChemicalList"][3]["RegistryNumber"] == "E1UOL152H7"
        assert record[1]["MedlineCitation"]["ChemicalList"][3]["NameOfSubstance"] == "Iron"
        assert record[1]["MedlineCitation"]["CitationSubset"] == ["IM"]
        assert record[1]["MedlineCitation"]["MeshHeadingList"][0]["DescriptorName"] == "Binding Sites"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][0]["DescriptorName"].attributes["MajorTopicYN"] == "N"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][1]["DescriptorName"] == "Chromatium"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][1]["DescriptorName"].attributes["MajorTopicYN"] == "N"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][1]["QualifierName"][0] == "enzymology"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][1]["QualifierName"][0].attributes["MajorTopicYN"] == "Y"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][2]["DescriptorName"] == "Cytochrome c Group"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][2]["DescriptorName"].attributes["MajorTopicYN"] == "Y"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][3]["DescriptorName"] == "Electron Spin Resonance Spectroscopy"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][3]["DescriptorName"].attributes["MajorTopicYN"] == "N"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][4]["DescriptorName"] == "Flavins"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][4]["DescriptorName"].attributes["MajorTopicYN"] == "N"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][5]["DescriptorName"] == "Heme"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][5]["DescriptorName"].attributes["MajorTopicYN"] == "N"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][6]["DescriptorName"] == "Hydrogen-Ion Concentration"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][6]["DescriptorName"].attributes["MajorTopicYN"] == "N"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][7]["DescriptorName"] == "Iron"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][7]["DescriptorName"].attributes["MajorTopicYN"] == "N"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][7]["QualifierName"][0] == "analysis"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][7]["QualifierName"][0].attributes["MajorTopicYN"] == "N"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][8]["DescriptorName"] == "Magnetics"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][8]["DescriptorName"].attributes["MajorTopicYN"] == "N"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][9]["DescriptorName"] == "Oxidation-Reduction"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][9]["DescriptorName"].attributes["MajorTopicYN"] == "N"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][10]["DescriptorName"] == "Protein Binding"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][10]["DescriptorName"].attributes["MajorTopicYN"] == "N"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][11]["DescriptorName"] == "Protein Conformation"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][11]["DescriptorName"].attributes["MajorTopicYN"] == "N"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][12]["DescriptorName"] == "Temperature"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][12]["DescriptorName"].attributes["MajorTopicYN"] == "N"
        assert record[1]["PubmedData"]["History"][0].attributes["PubStatus"] == "pubmed"
        assert record[1]["PubmedData"]["History"][0]["Year"] == "1976"
        assert record[1]["PubmedData"]["History"][0]["Month"] == "9"
        assert record[1]["PubmedData"]["History"][0]["Day"] == "28"
        assert record[1]["PubmedData"]["History"][1].attributes["PubStatus"] == "medline"
        assert record[1]["PubmedData"]["History"][1]["Year"] == "1976"
        assert record[1]["PubmedData"]["History"][1]["Month"] == "9"
        assert record[1]["PubmedData"]["History"][1]["Day"] == "28"
        assert record[1]["PubmedData"]["History"][1]["Hour"] == "0"
        assert record[1]["PubmedData"]["History"][1]["Minute"] == "1"
        assert record[1]["PubmedData"]["PublicationStatus"] == "ppublish"
        assert len(record[1]["PubmedData"]["ArticleIdList"]) == 3
        assert record[1]["PubmedData"]["ArticleIdList"][0] == "9997"
        assert record[1]["PubmedData"]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert record[1]["PubmedData"]["ArticleIdList"][1] == "10.1016/0005-2795(76)90109-4"
        assert record[1]["PubmedData"]["ArticleIdList"][1].attributes["IdType"] == "doi"
        assert record[1]["PubmedData"]["ArticleIdList"][2] == "0005-2795(76)90109-4"
        assert record[1]["PubmedData"]["ArticleIdList"][2].attributes["IdType"] == "pii"
        # fmt: on

    def test_pubmed2(self):
        """Test parsing XML returned by EFetch, PubMed database (second test)."""
        # In PubMed display PMIDs in xml retrieval mode.
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db='pubmed', id="11748933,11700088",
        #                       retmode="xml")
        with open("Entrez/pubmed2.xml", "rb") as stream:
            record = Entrez.read(stream)
        # fmt: off
        assert record["PubmedBookArticle"] == []
        record = record["PubmedArticle"]
        assert record[0]["MedlineCitation"].attributes["Owner"] == "NLM"
        assert record[0]["MedlineCitation"].attributes["Status"] == "MEDLINE"
        assert record[0]["MedlineCitation"]["PMID"] == "11748933"
        assert record[0]["MedlineCitation"]["DateCompleted"]["Year"] == "2002"
        assert record[0]["MedlineCitation"]["DateCompleted"]["Month"] == "03"
        assert record[0]["MedlineCitation"]["DateCompleted"]["Day"] == "04"
        assert record[0]["MedlineCitation"]["DateRevised"]["Year"] == "2006"
        assert record[0]["MedlineCitation"]["DateRevised"]["Month"] == "11"
        assert record[0]["MedlineCitation"]["DateRevised"]["Day"] == "15"
        assert record[0]["MedlineCitation"]["Article"].attributes["PubModel"] == "Print"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["ISSN"] == "0011-2240"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["ISSN"].attributes[
                "IssnType"
            ] == "Print"
        assert record[0]["MedlineCitation"]["Article"]["Journal"][
                "JournalIssue"
            ].attributes["CitedMedium"] == "Print"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                "Volume"
            ] == "42"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["Issue"] == "4"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                "PubDate"
            ]["Year"] == "2001"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                "PubDate"
            ]["Month"] == "Jun"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["Title"] == "Cryobiology"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["ISOAbbreviation"] == "Cryobiology"
        assert record[0]["MedlineCitation"]["Article"]["ArticleTitle"] == "Is cryopreservation a homogeneous process? Ultrastructure and motility of untreated, prefreezing, and postthawed spermatozoa of Diplodus puntazzo (Cetti)."
        assert record[0]["MedlineCitation"]["Article"]["Pagination"]["MedlinePgn"] == "244-55"
        assert record[0]["MedlineCitation"]["Article"]["Abstract"]["AbstractText"] == ["This study subdivides the cryopreservation procedure for Diplodus puntazzo spermatozoa into three key phases, fresh, prefreezing (samples equilibrated in cryosolutions), and postthawed stages, and examines the ultrastructural anomalies and motility profiles of spermatozoa in each stage, with different cryodiluents. Two simple cryosolutions were evaluated: 0.17 M sodium chloride containing a final concentration of 15% dimethyl sulfoxide (Me(2)SO) (cryosolution A) and 0.1 M sodium citrate containing a final concentration of 10% Me(2)SO (cryosolution B). Ultrastructural anomalies of the plasmatic and nuclear membranes of the sperm head were common and the severity of the cryoinjury differed significantly between the pre- and the postfreezing phases and between the two cryosolutions. In spermatozoa diluted with cryosolution A, during the prefreezing phase, the plasmalemma of 61% of the cells was absent or damaged compared with 24% in the fresh sample (P < 0.001). In spermatozoa diluted with cryosolution B, there was a pronounced increase in the number of cells lacking the head plasmatic membrane from the prefreezing to the postthawed stages (from 32 to 52%, P < 0.01). In both cryosolutions, damages to nuclear membrane were significantly higher after freezing (cryosolution A: 8 to 23%, P < 0.01; cryosolution B: 5 to 38%, P < 0.001). With cryosolution A, the after-activation motility profile confirmed a consistent drop from fresh at the prefreezing stage, whereas freezing and thawing did not affect the motility much further and 50% of the cells were immotile by 60-90 s after activation. With cryosolution B, only the postthawing stage showed a sharp drop of motility profile. This study suggests that the different phases of the cryoprocess should be investigated to better understand the process of sperm damage.",]
        assert record[0]["MedlineCitation"]["Article"]["Abstract"]["CopyrightInformation"] == "Copyright 2001 Elsevier Science."
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"].attributes["CompleteYN"] == "Y"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][0].attributes["ValidYN"] == "Y"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][0]["LastName"] == "Taddei"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][0]["ForeName"] == "A R"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][0]["Initials"] == "AR"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][0]["AffiliationInfo"][0]["Affiliation"] == "Dipartimento di Scienze Ambientali, Universit\xe0 degli Studi della Tuscia, 01100 Viterbo, Italy."
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][0]["AffiliationInfo"][0]["Identifier"] == []
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][1].attributes["ValidYN"] == "Y"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][1]["LastName"] == "Barbato"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][1]["ForeName"] == "F"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][1]["Initials"] == "F"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][2].attributes["ValidYN"] == "Y"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][2]["LastName"] == "Abelli"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][2]["ForeName"] == "L"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][2]["Initials"] == "L"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][3].attributes["ValidYN"] == "Y"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][3]["LastName"] == "Canese"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][3]["ForeName"] == "S"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][3]["Initials"] == "S"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][4].attributes["ValidYN"] == "Y"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][4]["LastName"] == "Moretti"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][4]["ForeName"] == "F"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][4]["Initials"] == "F"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][5].attributes["ValidYN"] == "Y"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][5]["LastName"] == "Rana"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][5]["ForeName"] == "K J"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][5]["Initials"] == "KJ"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][6].attributes[
                "ValidYN"
            ] == "Y"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][6]["LastName"] == "Fausto"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][6]["ForeName"] == "A M"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][6]["Initials"] == "AM"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][7].attributes[
                "ValidYN"
            ] == "Y"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][7]["LastName"] == "Mazzini"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][7]["ForeName"] == "M"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][7]["Initials"] == "M"
        assert record[0]["MedlineCitation"]["Article"]["Language"] == ["eng"]
        assert record[0]["MedlineCitation"]["Article"]["PublicationTypeList"][0] == "Journal Article"
        assert record[0]["MedlineCitation"]["Article"]["PublicationTypeList"][1] == "Research Support, Non-U.S. Gov't"
        assert record[0]["MedlineCitation"]["MedlineJournalInfo"]["Country"] == "Netherlands"
        assert record[0]["MedlineCitation"]["MedlineJournalInfo"]["MedlineTA"] == "Cryobiology"
        assert record[0]["MedlineCitation"]["MedlineJournalInfo"]["NlmUniqueID"] == "0006252"
        assert record[0]["MedlineCitation"]["CitationSubset"] == ["IM"]
        assert record[0]["MedlineCitation"]["MeshHeadingList"][0]["DescriptorName"] == "Animals"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][0][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][1]["DescriptorName"] == "Cell Membrane"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][1][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][1]["QualifierName"][0] == "ultrastructure"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][1]["QualifierName"][
                0
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][2]["DescriptorName"] == "Cryopreservation"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][2][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][2]["QualifierName"][0] == "methods"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][2]["QualifierName"][
                0
            ].attributes["MajorTopicYN"] == "Y"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][3]["DescriptorName"] == "Male"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][3][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][4]["DescriptorName"] == "Microscopy, Electron"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][4][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][5]["DescriptorName"] == "Microscopy, Electron, Scanning"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][5][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][6]["DescriptorName"] == "Nuclear Envelope"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][6][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][6]["QualifierName"][0] == "ultrastructure"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][6]["QualifierName"][
                0
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][7]["DescriptorName"] == "Sea Bream"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][7][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][7]["QualifierName"][0] == "anatomy & histology"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][7]["QualifierName"][
                0
            ].attributes["MajorTopicYN"] == "Y"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][7]["QualifierName"][1] == "physiology"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][7]["QualifierName"][
                1
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][8]["DescriptorName"] == "Semen Preservation"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][8][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][8]["QualifierName"][0] == "adverse effects"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][8]["QualifierName"][
                0
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][8]["QualifierName"][1] == "methods"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][8]["QualifierName"][
                1
            ].attributes["MajorTopicYN"] == "Y"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][9]["DescriptorName"] == "Sperm Motility"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][9][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "Y"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][10]["DescriptorName"] == "Spermatozoa"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][10][
                "DescriptorName"
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][10]["QualifierName"][0] == "physiology"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][10]["QualifierName"][
                0
            ].attributes["MajorTopicYN"] == "N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][10]["QualifierName"][1] == "ultrastructure"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][10]["QualifierName"][
                1
            ].attributes["MajorTopicYN"] == "Y"
        assert record[0]["PubmedData"]["History"][0].attributes["PubStatus"] == "pubmed"
        assert record[0]["PubmedData"]["History"][0]["Year"] == "2001"
        assert record[0]["PubmedData"]["History"][0]["Month"] == "12"
        assert record[0]["PubmedData"]["History"][0]["Day"] == "26"
        assert record[0]["PubmedData"]["History"][0]["Hour"] == "10"
        assert record[0]["PubmedData"]["History"][0]["Minute"] == "0"
        assert record[0]["PubmedData"]["History"][1].attributes["PubStatus"] == "medline"
        assert record[0]["PubmedData"]["History"][1]["Year"] == "2002"
        assert record[0]["PubmedData"]["History"][1]["Month"] == "3"
        assert record[0]["PubmedData"]["History"][1]["Day"] == "5"
        assert record[0]["PubmedData"]["History"][1]["Hour"] == "10"
        assert record[0]["PubmedData"]["History"][1]["Minute"] == "1"
        assert record[0]["PubmedData"]["PublicationStatus"] == "ppublish"
        assert record[0]["PubmedData"]["ArticleIdList"][0] == "11748933"
        assert record[0]["PubmedData"]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert record[0]["PubmedData"]["ArticleIdList"][1] == "10.1006/cryo.2001.2328"
        assert record[0]["PubmedData"]["ArticleIdList"][1].attributes["IdType"] == "doi"
        assert record[0]["PubmedData"]["ArticleIdList"][2] == "S0011-2240(01)92328-4"
        assert record[0]["PubmedData"]["ArticleIdList"][2].attributes["IdType"] == "pii"

        assert record[1]["MedlineCitation"].attributes["Owner"] == "NLM"
        assert record[1]["MedlineCitation"].attributes["Status"] == "PubMed-not-MEDLINE"
        assert record[1]["MedlineCitation"]["PMID"] == "11700088"
        assert record[1]["MedlineCitation"]["DateCompleted"]["Year"] == "2001"
        assert record[1]["MedlineCitation"]["DateCompleted"]["Month"] == "12"
        assert record[1]["MedlineCitation"]["DateCompleted"]["Day"] == "20"
        assert record[1]["MedlineCitation"]["DateRevised"]["Year"] == "2003"
        assert record[1]["MedlineCitation"]["DateRevised"]["Month"] == "10"
        assert record[1]["MedlineCitation"]["DateRevised"]["Day"] == "31"
        assert record[1]["MedlineCitation"]["Article"].attributes["PubModel"] == "Print"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["ISSN"] == "1090-7807"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["ISSN"].attributes[
                "IssnType"
            ] == "Print"
        assert record[1]["MedlineCitation"]["Article"]["Journal"][
                "JournalIssue"
            ].attributes["CitedMedium"] == "Print"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                "Volume"
            ] == "153"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["Issue"] == "1"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                "PubDate"
            ]["Year"] == "2001"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                "PubDate"
            ]["Month"] == "Nov"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["Title"] == "Journal of magnetic resonance (San Diego, Calif. : 1997)"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["ISOAbbreviation"] == "J Magn Reson"
        assert record[1]["MedlineCitation"]["Article"]["ArticleTitle"] == "Proton MRI of (13)C distribution by J and chemical shift editing."
        assert record[1]["MedlineCitation"]["Article"]["Pagination"]["MedlinePgn"] == "117-23"
        assert record[1]["MedlineCitation"]["Article"]["Abstract"]["AbstractText"] == ["The sensitivity of (13)C NMR imaging can be considerably favored by detecting the (1)H nuclei bound to (13)C nuclei via scalar J-interaction (X-filter). However, the J-editing approaches have difficulty in discriminating between compounds with similar J-constant as, for example, different glucose metabolites. In such cases, it is almost impossible to get J-edited images of a single-compound distribution, since the various molecules are distinguishable only via their chemical shift. In a recent application of J-editing to high-resolution spectroscopy, it has been shown that a more efficient chemical selectivity could be obtained by utilizing the larger chemical shift range of (13)C. This has been made by introducing frequency-selective (13)C pulses that allow a great capability of indirect chemical separation. Here a double-resonance imaging approach is proposed, based on both J-editing and (13)C chemical shift editing, which achieves a powerful chemical selectivity and is able to produce full maps of specific chemical compounds. Results are presented on a multicompartments sample containing solutions of glucose and lactic and glutamic acid in water."]
        assert record[1]["MedlineCitation"]["Article"]["Abstract"]["CopyrightInformation"] == "Copyright 2001 Academic Press."
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"].attributes["CompleteYN"] == "Y"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][0].attributes["ValidYN"] == "Y"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][0]["LastName"] == "Casieri"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][0]["ForeName"] == "C"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][0]["Initials"] == "C"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][0]["Identifier"] == []
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][0]["AffiliationInfo"][0]["Affiliation"] == "INFM and Department of Physics, University of L'Aquila, I-67100 L'Aquila, Italy."
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][1].attributes["ValidYN"] == "Y"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][1]["LastName"] == "Testa"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][1]["ForeName"] == "C"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][1]["Initials"] == "C"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][2].attributes["ValidYN"] == "Y"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][2]["LastName"] == "Carpinelli"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][2]["ForeName"] == "G"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][2]["Initials"] == "G"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][3].attributes["ValidYN"] == "Y"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][3]["LastName"] == "Canese"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][3]["ForeName"] == "R"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][3]["Initials"] == "R"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][4].attributes["ValidYN"] == "Y"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][4]["LastName"] == "Podo"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][4]["ForeName"] == "F"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][4]["Initials"] == "F"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][5].attributes["ValidYN"] == "Y"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][5]["LastName"] == "De Luca"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][5]["ForeName"] == "F"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][5]["Initials"] == "F"
        assert record[1]["MedlineCitation"]["Article"]["Language"] == ["eng"]
        assert record[1]["MedlineCitation"]["Article"]["PublicationTypeList"][0] == "Journal Article"
        assert record[1]["MedlineCitation"]["MedlineJournalInfo"]["Country"] == "United States"
        assert record[1]["MedlineCitation"]["MedlineJournalInfo"]["MedlineTA"] == "J Magn Reson"
        assert record[1]["MedlineCitation"]["MedlineJournalInfo"]["NlmUniqueID"] == "9707935"
        assert record[1]["PubmedData"]["History"][0].attributes["PubStatus"] == "pubmed"
        assert record[1]["PubmedData"]["History"][0]["Year"] == "2001"
        assert record[1]["PubmedData"]["History"][0]["Month"] == "11"
        assert record[1]["PubmedData"]["History"][0]["Day"] == "9"
        assert record[1]["PubmedData"]["History"][0]["Hour"] == "10"
        assert record[1]["PubmedData"]["History"][0]["Minute"] == "0"
        assert record[1]["PubmedData"]["History"][1].attributes["PubStatus"] == "medline"
        assert record[1]["PubmedData"]["History"][1]["Year"] == "2001"
        assert record[1]["PubmedData"]["History"][1]["Month"] == "11"
        assert record[1]["PubmedData"]["History"][1]["Day"] == "9"
        assert record[1]["PubmedData"]["History"][1]["Hour"] == "10"
        assert record[1]["PubmedData"]["History"][1]["Minute"] == "1"
        assert record[1]["PubmedData"]["PublicationStatus"] == "ppublish"
        assert record[1]["PubmedData"]["ArticleIdList"][0] == "11700088"
        assert record[1]["PubmedData"]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert record[1]["PubmedData"]["ArticleIdList"][1] == "10.1006/jmre.2001.2429"
        assert record[1]["PubmedData"]["ArticleIdList"][1].attributes["IdType"] == "doi"
        assert record[1]["PubmedData"]["ArticleIdList"][2] == "S1090-7807(01)92429-2"
        assert record[1]["PubmedData"]["ArticleIdList"][2].attributes["IdType"] == "pii"
        # fmt: on

    def test_pubmed_html_tags(self):
        """Test parsing XML returned by EFetch, PubMed database with HTML tags."""
        # In PubMed display PMIDs in xml retrieval mode.
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db='pubmed', retmode='xml', id='29106400')
        with open("Entrez/pubmed4.xml", "rb") as stream:
            records = Entrez.read(stream)
        # fmt: off
        assert len(records) == 2
        assert len(records["PubmedBookArticle"]) == 0
        assert len(records["PubmedArticle"]) == 1
        assert records["PubmedArticle"][0]["MedlineCitation"].attributes["Status"] == "MEDLINE"
        assert records["PubmedArticle"][0]["MedlineCitation"].attributes["Owner"] == "NLM"
        assert records["PubmedArticle"][0]["MedlineCitation"]["PMID"] == "27797938"
        assert records["PubmedArticle"][0]["MedlineCitation"]["PMID"].attributes[
                "Version"
            ] == "1"
        assert records["PubmedArticle"][0]["MedlineCitation"]["DateCompleted"]["Year"] == "2017"
        assert records["PubmedArticle"][0]["MedlineCitation"]["DateCompleted"]["Month"] == "08"
        assert records["PubmedArticle"][0]["MedlineCitation"]["DateCompleted"]["Day"] == "03"
        assert records["PubmedArticle"][0]["MedlineCitation"]["DateRevised"]["Year"] == "2018"
        assert records["PubmedArticle"][0]["MedlineCitation"]["DateRevised"]["Month"] == "04"
        assert records["PubmedArticle"][0]["MedlineCitation"]["DateRevised"]["Day"] == "17"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"].attributes[
                "PubModel"
            ] == "Print-Electronic"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Journal"][
                "ISSN"
            ] == "1468-3288"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Journal"][
                "ISSN"
            ].attributes["IssnType"] == "Electronic"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Journal"][
                "JournalIssue"
            ].attributes["CitedMedium"] == "Internet"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Journal"][
                "JournalIssue"
            ]["Volume"] == "66"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Journal"][
                "JournalIssue"
            ]["Issue"] == "6"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Journal"][
                "JournalIssue"
            ]["PubDate"]["Year"] == "2017"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Journal"][
                "JournalIssue"
            ]["PubDate"]["Month"] == "06"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Journal"][
                "Title"
            ] == "Gut"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Journal"][
                "ISOAbbreviation"
            ] == "Gut"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleTitle"] == "Leucocyte telomere length, genetic variants at the <i>TERT</i> gene region and risk of pancreatic cancer."
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Pagination"][
                "MedlinePgn"
            ] == "1116-1122"
        assert len(
                records["PubmedArticle"][0]["MedlineCitation"]["Article"]["ELocationID"]
            ) == 1
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"]["ELocationID"][0] == "10.1136/gutjnl-2016-312510"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"]["ELocationID"][
                0
            ].attributes["EIdType"] == "doi"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"]["ELocationID"][
                0
            ].attributes["ValidYN"] == "Y"
        assert len(records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"]) == 2
        assert len(
                records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"][
                    "AbstractText"
                ]
            ) == 4
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"][
                "AbstractText"
            ][0] == "Telomere shortening occurs as an early event in pancreatic tumorigenesis, and genetic variants at the telomerase reverse transcriptase (<i>TERT</i>) gene region have been associated with pancreatic cancer risk. However, it is unknown whether prediagnostic leucocyte telomere length is associated with subsequent risk of pancreatic cancer."
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"][
                "AbstractText"
            ][0].attributes["Label"] == "OBJECTIVE"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"][
                "AbstractText"
            ][1] == "We measured prediagnostic leucocyte telomere length in 386 pancreatic cancer cases and 896 matched controls from five prospective US cohorts. ORs and 95% CIs were calculated using conditional logistic regression. Matching factors included year of birth, cohort (which also matches on sex), smoking status, fasting status and month/year of blood collection. We additionally examined single-nucleotide polymorphisms (SNPs) at the <i>TERT</i> region in relation to pancreatic cancer risk and leucocyte telomere length using logistic and linear regression, respectively."
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"][
                "AbstractText"
            ][1].attributes["Label"] == "DESIGN"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"][
                "AbstractText"
            ][2] == "Shorter prediagnostic leucocyte telomere length was associated with higher risk of pancreatic cancer (comparing extreme quintiles of telomere length, OR 1.72; 95% CI 1.07 to 2.78; p<sub>trend</sub>=0.048). Results remained unchanged after adjustment for diabetes, body mass index and physical activity. Three SNPs at <i>TERT</i> (linkage disequilibrium r<sup>2</sup><0.25) were associated with pancreatic cancer risk, including rs401681 (per minor allele OR 1.33; 95% CI 1.12 to 1.59; p=0.002), rs2736100 (per minor allele OR 1.36; 95% CI 1.13 to 1.63; p=0.001) and rs2736098 (per minor allele OR 0.75; 95% CI 0.63 to 0.90; p=0.002). The minor allele for rs401681 was associated with shorter telomere length (p=0.023)."
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"][
                "AbstractText"
            ][2].attributes["Label"] == "RESULTS"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"][
                "AbstractText"
            ][3] == "Prediagnostic leucocyte telomere length and genetic variants at the <i>TERT</i> gene region were associated with risk of pancreatic cancer."
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"][
                "AbstractText"
            ][3].attributes["Label"] == "CONCLUSIONS"
        assert len(
                records["PubmedArticle"][0]["MedlineCitation"]["Article"]["AuthorList"]
            ) == 22
        assert len(records["PubmedArticle"][0]["MedlineCitation"]["Article"]["GrantList"]) == 35
        assert len(
                records["PubmedArticle"][0]["MedlineCitation"]["Article"][
                    "PublicationTypeList"
                ]
            ) == 5
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"][
                "PublicationTypeList"
            ][0] == "Journal Article"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"][
                "PublicationTypeList"
            ][0].attributes["UI"] == "D016428"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"][
                "PublicationTypeList"
            ][1] == "Observational Study"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"][
                "PublicationTypeList"
            ][1].attributes["UI"] == "D064888"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"][
                "PublicationTypeList"
            ][2] == "Research Support, N.I.H., Extramural"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"][
                "PublicationTypeList"
            ][2].attributes["UI"] == "D052061"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"][
                "PublicationTypeList"
            ][3] == "Research Support, U.S. Gov't, Non-P.H.S."
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"][
                "PublicationTypeList"
            ][3].attributes["UI"] == "D013486"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"][
                "PublicationTypeList"
            ][4] == "Research Support, Non-U.S. Gov't"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"][
                "PublicationTypeList"
            ][4].attributes["UI"] == "D013485"
        assert len(
                records["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleDate"]
            ) == 1
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleDate"][
                0
            ].attributes["DateType"] == "Electronic"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleDate"][0][
                "Year"
            ] == "2016"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleDate"][0][
                "Month"
            ] == "10"
        assert records["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleDate"][0][
                "Day"
            ] == "21"
        assert records["PubmedArticle"][0]["MedlineCitation"]["MedlineJournalInfo"][
                "Country"
            ] == "England"
        assert records["PubmedArticle"][0]["MedlineCitation"]["MedlineJournalInfo"][
                "MedlineTA"
            ] == "Gut"
        assert records["PubmedArticle"][0]["MedlineCitation"]["MedlineJournalInfo"][
                "NlmUniqueID"
            ] == "2985108R"
        assert records["PubmedArticle"][0]["MedlineCitation"]["MedlineJournalInfo"][
                "ISSNLinking"
            ] == "0017-5749"
        assert len(records["PubmedArticle"][0]["MedlineCitation"]["ChemicalList"]) == 2
        # fmt: on

    def test_pubmed_html_escaping(self):
        """Test parsing XML returned by EFetch, PubMed database with HTML tags and HTML escape characters."""
        # In PubMed display PMIDs in xml retrieval mode.
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db='pubmed', retmode='xml', id='28775130')
        with open("Entrez/pubmed5.xml", "rb") as stream:
            record = Entrez.read(stream, escape=True)
        assert len(record) == 2
        assert len(record["PubmedArticle"]) == 1
        assert len(record["PubmedBookArticle"]) == 0
        article = record["PubmedArticle"][0]
        assert len(article) == 2
        assert len(article["PubmedData"]) == 3
        assert len(article["PubmedData"]["ArticleIdList"]) == 5
        assert article["PubmedData"]["ArticleIdList"][0] == "28775130"
        assert article["PubmedData"]["ArticleIdList"][0].attributes == {"IdType": "pubmed"}
        assert article["PubmedData"]["ArticleIdList"][1] == "oemed-2017-104431"
        assert article["PubmedData"]["ArticleIdList"][1].attributes == {"IdType": "pii"}
        assert article["PubmedData"]["ArticleIdList"][2] == "10.1136/oemed-2017-104431"
        assert article["PubmedData"]["ArticleIdList"][2].attributes == {"IdType": "doi"}
        assert article["PubmedData"]["ArticleIdList"][3] == "PMC5771820"
        assert article["PubmedData"]["ArticleIdList"][3].attributes == {"IdType": "pmc"}
        assert article["PubmedData"]["ArticleIdList"][4] == "NIHMS932407"
        assert article["PubmedData"]["ArticleIdList"][4].attributes == {"IdType": "mid"}
        assert article["PubmedData"]["PublicationStatus"] == "ppublish"
        assert len(article["PubmedData"]["History"]) == 7
        assert len(article["PubmedData"]["History"][0]) == 3
        assert article["PubmedData"]["History"][0]["Year"] == "2017"
        assert article["PubmedData"]["History"][0]["Month"] == "03"
        assert article["PubmedData"]["History"][0]["Day"] == "10"
        assert article["PubmedData"]["History"][0].attributes == {"PubStatus": "received"}
        assert len(article["PubmedData"]["History"][1]) == 3
        assert article["PubmedData"]["History"][1]["Year"] == "2017"
        assert article["PubmedData"]["History"][1]["Month"] == "06"
        assert article["PubmedData"]["History"][1]["Day"] == "13"
        assert article["PubmedData"]["History"][1].attributes == {"PubStatus": "revised"}
        assert len(article["PubmedData"]["History"][2]) == 3
        assert article["PubmedData"]["History"][2]["Year"] == "2017"
        assert article["PubmedData"]["History"][2]["Month"] == "06"
        assert article["PubmedData"]["History"][2]["Day"] == "22"
        assert article["PubmedData"]["History"][2].attributes == {"PubStatus": "accepted"}
        assert len(article["PubmedData"]["History"][3]) == 3
        assert article["PubmedData"]["History"][3]["Year"] == "2019"
        assert article["PubmedData"]["History"][3]["Month"] == "02"
        assert article["PubmedData"]["History"][3]["Day"] == "01"
        assert article["PubmedData"]["History"][3].attributes == {"PubStatus": "pmc-release"}
        assert len(article["PubmedData"]["History"][4]) == 5
        assert article["PubmedData"]["History"][4]["Year"] == "2017"
        assert article["PubmedData"]["History"][4]["Month"] == "8"
        assert article["PubmedData"]["History"][4]["Day"] == "5"
        assert article["PubmedData"]["History"][4]["Hour"] == "6"
        assert article["PubmedData"]["History"][4]["Minute"] == "0"
        assert article["PubmedData"]["History"][4].attributes == {"PubStatus": "pubmed"}
        assert len(article["PubmedData"]["History"][5]) == 5
        assert article["PubmedData"]["History"][5]["Year"] == "2017"
        assert article["PubmedData"]["History"][5]["Month"] == "8"
        assert article["PubmedData"]["History"][5]["Day"] == "5"
        assert article["PubmedData"]["History"][5]["Hour"] == "6"
        assert article["PubmedData"]["History"][5]["Minute"] == "0"
        assert article["PubmedData"]["History"][5].attributes == {"PubStatus": "medline"}
        assert len(article["PubmedData"]["History"][6]) == 5
        assert article["PubmedData"]["History"][6]["Year"] == "2017"
        assert article["PubmedData"]["History"][6]["Month"] == "8"
        assert article["PubmedData"]["History"][6]["Day"] == "5"
        assert article["PubmedData"]["History"][6]["Hour"] == "6"
        assert article["PubmedData"]["History"][6]["Minute"] == "0"
        assert article["PubmedData"]["History"][6].attributes == {"PubStatus": "entrez"}
        assert len(article["MedlineCitation"]) == 12
        assert len(article["MedlineCitation"]["CitationSubset"]) == 0
        assert article["MedlineCitation"]["CoiStatement"] == "Competing interests: None declared."
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"]) == 40
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][0]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][0]["RefSource"] == "J Toxicol Environ Health A. 2003 Jun 13;66(11):965-86"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][0]["PMID"] == "12775511"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][0].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][1]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][1]["RefSource"] == "Ann Intern Med. 2015 May 5;162(9):641-50"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][1]["PMID"] == "25798805"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][1].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][2]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][2]["RefSource"] == "Cancer Causes Control. 1999 Dec;10(6):583-95"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][2]["PMID"] == "10616827"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][2].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][3]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][3]["RefSource"] == "Thyroid. 2010 Jul;20(7):755-61"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][3]["PMID"] == "20578899"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][3].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][4]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][4]["RefSource"] == "Environ Health Perspect. 1999 Mar;107(3):205-11"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][4]["PMID"] == "10064550"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][4].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][5]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][5]["RefSource"] == "J Clin Endocrinol Metab. 2006 Nov;91(11):4295-301"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][5]["PMID"] == "16868053"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][5].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][6]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][6]["RefSource"] == "Endocrinology. 1998 Oct;139(10):4252-63"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][6]["PMID"] == "9751507"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][6].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][7]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][7]["RefSource"] == "Eur J Endocrinol. 2016 Apr;174(4):409-14"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][7]["PMID"] == "26863886"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][7].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][8]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][8]["RefSource"] == "Eur J Endocrinol. 2000 Nov;143(5):639-47"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][8]["PMID"] == "11078988"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][8].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][9]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][9]["RefSource"] == "Environ Res. 2016 Nov;151:389-398"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][9]["PMID"] == "27540871"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][9].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][10]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][10]["RefSource"] == "Am J Epidemiol. 2010 Jan 15;171(2):242-52"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][10]["PMID"] == "19951937"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][10].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][11]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][11]["RefSource"] == "Thyroid. 1998 Sep;8(9):827-56"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][11]["PMID"] == "9777756"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][11].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][12]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][12]["RefSource"] == "Curr Opin Pharmacol. 2001 Dec;1(6):626-31"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][12]["PMID"] == "11757819"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][12].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][13]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][13]["RefSource"] == "Breast Cancer Res Treat. 2012 Jun;133(3):1169-77"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][13]["PMID"] == "22434524"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][13].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][14]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][14]["RefSource"] == "Int J Environ Res Public Health. 2011 Dec;8(12 ):4608-22"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][14]["PMID"] == "22408592"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][14].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][15]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][15]["RefSource"] == "Ann Oncol. 2014 Oct;25(10):2025-30"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][15]["PMID"] == "25081899"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][15].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][16]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][16]["RefSource"] == "Environ Health. 2006 Dec 06;5:32"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][16]["PMID"] == "17147831"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][16].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][17]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][17]["RefSource"] == "Environ Health Perspect. 1998 Aug;106(8):437-45"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][17]["PMID"] == "9681970"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][17].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][18]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][18]["RefSource"] == "Arch Intern Med. 2000 Feb 28;160(4):526-34"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][18]["PMID"] == "10695693"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][18].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][19]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][19]["RefSource"] == "Endocrine. 2011 Jun;39(3):259-65"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][19]["PMID"] == "21161440"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][19].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][20]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][20]["RefSource"] == "Cancer Epidemiol Biomarkers Prev. 2008 Aug;17(8):1880-3"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][20]["PMID"] == "18708375"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][20].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][21]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][21]["RefSource"] == "Am J Epidemiol. 2010 Feb 15;171(4):455-64"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][21]["PMID"] == "20061368"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][21].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][22]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][22]["RefSource"] == "J Clin Endocrinol Metab. 2002 Feb;87(2):489-99"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][22]["PMID"] == "11836274"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][22].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][23]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][23]["RefSource"] == "J Toxicol Environ Health A. 2015 ;78(21-22):1338-47"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][23]["PMID"] == "26555155"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][23].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][24]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][24]["RefSource"] == "Toxicol Sci. 2002 Jun;67(2):207-18"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][24]["PMID"] == "12011480"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][24].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][25]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][25]["RefSource"] == "Natl Cancer Inst Carcinog Tech Rep Ser. 1978;21:1-184"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][25]["PMID"] == "12844187"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][25].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][26]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][26]["RefSource"] == "Environ Res. 2013 Nov;127:7-15"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][26]["PMID"] == "24183346"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][26].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][27]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][27]["RefSource"] == "JAMA. 2004 Jan 14;291(2):228-38"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][27]["PMID"] == "14722150"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][27].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][28]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][28]["RefSource"] == "J Expo Sci Environ Epidemiol. 2010 Sep;20(6):559-69"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][28]["PMID"] == "19888312"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][28].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][29]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][29]["RefSource"] == "Environ Health Perspect. 1996 Apr;104(4):362-9"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][29]["PMID"] == "8732939"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][29].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][30]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][30]["RefSource"] == "Lancet. 2012 Mar 24;379(9821):1142-54"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][30]["PMID"] == "22273398"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][30].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][31]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][31]["RefSource"] == "JAMA. 1995 Mar 8;273(10):808-12"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][31]["PMID"] == "7532241"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][31].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][32]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][32]["RefSource"] == "Sci Total Environ. 2002 Aug 5;295(1-3):207-15"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][32]["PMID"] == "12186288"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][32].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][33]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][33]["RefSource"] == "Eur J Endocrinol. 2006 May;154(5):599-611"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][33]["PMID"] == "16645005"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][33].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][34]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][34]["RefSource"] == "J Occup Environ Med. 2013 Oct;55(10):1171-8"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][34]["PMID"] == "24064777"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][34].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][35]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][35]["RefSource"] == "Thyroid. 2007 Sep;17(9):811-7"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][35]["PMID"] == "17956155"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][35].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][36]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][36]["RefSource"] == "Rev Environ Contam Toxicol. 1991;120:1-82"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][36]["PMID"] == "1899728"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][36].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][37]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][37]["RefSource"] == "Environ Health Perspect. 1997 Oct;105(10):1126-30"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][37]["PMID"] == "9349837"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][37].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][38]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][38]["RefSource"] == "J Biochem Mol Toxicol. 2005;19(3):175"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][38]["PMID"] == "15977190"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][38].attributes == {"RefType": "Cites"}
        assert len(article["MedlineCitation"]["CommentsCorrectionsList"][39]) == 2
        assert article["MedlineCitation"]["CommentsCorrectionsList"][39]["RefSource"] == "Immunogenetics. 2002 Jun;54(3):141-57"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][39]["PMID"] == "12073143"
        assert article["MedlineCitation"]["CommentsCorrectionsList"][39].attributes == {"RefType": "Cites"}
        assert article["MedlineCitation"]["DateRevised"]["Year"] == "2018"
        assert article["MedlineCitation"]["DateRevised"]["Month"] == "04"
        assert article["MedlineCitation"]["DateRevised"]["Day"] == "25"
        assert len(article["MedlineCitation"]["DateRevised"].attributes) == 0
        assert len(article["MedlineCitation"]["GeneralNote"]) == 0
        assert len(article["MedlineCitation"]["KeywordList"]) == 1
        assert len(article["MedlineCitation"]["KeywordList"][0]) == 5
        assert article["MedlineCitation"]["KeywordList"][0][0] == "agriculture"
        assert article["MedlineCitation"]["KeywordList"][0][0].attributes == {"MajorTopicYN": "N"}
        assert article["MedlineCitation"]["KeywordList"][0][1] == "hypothyroidism"
        assert article["MedlineCitation"]["KeywordList"][0][1].attributes == {"MajorTopicYN": "N"}
        assert article["MedlineCitation"]["KeywordList"][0][2] == "pesticides"
        assert article["MedlineCitation"]["KeywordList"][0][2].attributes == {"MajorTopicYN": "N"}
        assert article["MedlineCitation"]["KeywordList"][0][3] == "thyroid disease"
        assert article["MedlineCitation"]["KeywordList"][0][3].attributes == {"MajorTopicYN": "N"}
        assert article["MedlineCitation"]["KeywordList"][0][4] == "thyroid stimulating hormone"
        assert article["MedlineCitation"]["KeywordList"][0][4].attributes == {"MajorTopicYN": "N"}
        assert len(article["MedlineCitation"]["MedlineJournalInfo"]) == 4
        assert article["MedlineCitation"]["MedlineJournalInfo"]["MedlineTA"] == "Occup Environ Med"
        assert article["MedlineCitation"]["MedlineJournalInfo"]["Country"] == "England"
        assert article["MedlineCitation"]["MedlineJournalInfo"]["NlmUniqueID"] == "9422759"
        assert article["MedlineCitation"]["MedlineJournalInfo"]["ISSNLinking"] == "1351-0711"
        assert len(article["MedlineCitation"]["MedlineJournalInfo"].attributes) == 0
        assert len(article["MedlineCitation"]["OtherAbstract"]) == 0
        assert len(article["MedlineCitation"]["OtherID"]) == 0
        assert article["MedlineCitation"]["PMID"] == "28775130"
        assert len(article["MedlineCitation"]["SpaceFlightMission"]) == 0
        assert len(article["MedlineCitation"]["Article"]["ArticleDate"]) == 1
        assert len(article["MedlineCitation"]["Article"]["ArticleDate"][0]) == 3
        assert article["MedlineCitation"]["Article"]["ArticleDate"][0]["Month"] == "08"
        assert article["MedlineCitation"]["Article"]["ArticleDate"][0]["Day"] == "03"
        assert article["MedlineCitation"]["Article"]["ArticleDate"][0]["Year"] == "2017"
        assert article["MedlineCitation"]["Article"]["ArticleDate"][0].attributes == {"DateType": "Electronic"}
        assert len(article["MedlineCitation"]["Article"]["Pagination"]) == 1
        assert article["MedlineCitation"]["Article"]["Pagination"]["MedlinePgn"] == "79-89"
        assert len(article["MedlineCitation"]["Article"]["Pagination"].attributes) == 0
        assert len(article["MedlineCitation"]["Article"]["AuthorList"]) == 12
        assert article["MedlineCitation"]["Article"]["AuthorList"][0]["LastName"] == "Lerro"
        assert article["MedlineCitation"]["Article"]["AuthorList"][0]["Initials"] == "CC"
        assert article["MedlineCitation"]["Article"]["AuthorList"][0]["Identifier"] == []
        assert len(
                article["MedlineCitation"]["Article"]["AuthorList"][0][
                    "AffiliationInfo"
                ]
            ) == 1
        assert article["MedlineCitation"]["Article"]["AuthorList"][0]["AffiliationInfo"][
                0
            ]["Affiliation"] == "Division of Cancer Epidemiology and Genetics, National Cancer Institute, Rockville, Maryland, USA."
        assert article["MedlineCitation"]["Article"]["AuthorList"][0]["AffiliationInfo"][
                0
            ]["Identifier"] == []
        assert len(
                article["MedlineCitation"]["Article"]["AuthorList"][0][
                    "AffiliationInfo"
                ][0].attributes
            ) == 0
        assert article["MedlineCitation"]["Article"]["AuthorList"][0]["ForeName"] == "Catherine C"
        assert article["MedlineCitation"]["Article"]["AuthorList"][0].attributes == {"ValidYN": "Y"}
        assert article["MedlineCitation"]["Article"]["AuthorList"][1]["LastName"] == "Beane Freeman"
        assert article["MedlineCitation"]["Article"]["AuthorList"][1]["Initials"] == "LE"
        assert article["MedlineCitation"]["Article"]["AuthorList"][1]["Identifier"] == []
        assert len(
                article["MedlineCitation"]["Article"]["AuthorList"][1][
                    "AffiliationInfo"
                ]
            ) == 1
        assert article["MedlineCitation"]["Article"]["AuthorList"][1]["AffiliationInfo"][
                0
            ]["Affiliation"] == "Division of Cancer Epidemiology and Genetics, National Cancer Institute, Rockville, Maryland, USA."
        assert article["MedlineCitation"]["Article"]["AuthorList"][1]["AffiliationInfo"][
                0
            ]["Identifier"] == []
        assert len(
                article["MedlineCitation"]["Article"]["AuthorList"][1][
                    "AffiliationInfo"
                ][0].attributes
            ) == 0
        assert article["MedlineCitation"]["Article"]["AuthorList"][1]["ForeName"] == "Laura E"
        assert article["MedlineCitation"]["Article"]["AuthorList"][1].attributes == {"ValidYN": "Y"}
        assert article["MedlineCitation"]["Article"]["AuthorList"][2]["LastName"] == "DellaValle"
        assert article["MedlineCitation"]["Article"]["AuthorList"][2]["Initials"] == "CT"
        assert article["MedlineCitation"]["Article"]["AuthorList"][2]["Identifier"] == []
        assert len(
                article["MedlineCitation"]["Article"]["AuthorList"][2][
                    "AffiliationInfo"
                ]
            ) == 2
        assert article["MedlineCitation"]["Article"]["AuthorList"][2]["AffiliationInfo"][
                0
            ]["Affiliation"] == "Division of Cancer Epidemiology and Genetics, National Cancer Institute, Rockville, Maryland, USA."
        assert article["MedlineCitation"]["Article"]["AuthorList"][2]["AffiliationInfo"][
                0
            ]["Identifier"] == []
        assert len(
                article["MedlineCitation"]["Article"]["AuthorList"][2][
                    "AffiliationInfo"
                ][0].attributes
            ) == 0
        assert article["MedlineCitation"]["Article"]["AuthorList"][2]["AffiliationInfo"][
                1
            ]["Affiliation"] == "Environmental Working Group, Washington, DC, USA."
        assert article["MedlineCitation"]["Article"]["AuthorList"][2]["AffiliationInfo"][
                1
            ]["Identifier"] == []
        assert len(
                article["MedlineCitation"]["Article"]["AuthorList"][2][
                    "AffiliationInfo"
                ][1].attributes
            ) == 0
        assert article["MedlineCitation"]["Article"]["AuthorList"][2]["ForeName"] == "Curt T"
        assert article["MedlineCitation"]["Article"]["AuthorList"][2].attributes == {"ValidYN": "Y"}
        assert article["MedlineCitation"]["Article"]["AuthorList"][3]["LastName"] == "Kibriya"
        assert article["MedlineCitation"]["Article"]["AuthorList"][3]["Initials"] == "MG"
        assert article["MedlineCitation"]["Article"]["AuthorList"][3]["Identifier"] == []
        assert len(
                article["MedlineCitation"]["Article"]["AuthorList"][3][
                    "AffiliationInfo"
                ]
            ) == 1
        assert article["MedlineCitation"]["Article"]["AuthorList"][3]["AffiliationInfo"][
                0
            ]["Affiliation"] == "Department of Public Health Sciences, The University of Chicago, Chicago, Illinois, USA."
        assert article["MedlineCitation"]["Article"]["AuthorList"][3]["AffiliationInfo"][
                0
            ]["Identifier"] == []
        assert len(
                article["MedlineCitation"]["Article"]["AuthorList"][3][
                    "AffiliationInfo"
                ][0].attributes
            ) == 0
        assert article["MedlineCitation"]["Article"]["AuthorList"][3]["ForeName"] == "Muhammad G"
        assert article["MedlineCitation"]["Article"]["AuthorList"][3].attributes == {"ValidYN": "Y"}
        assert article["MedlineCitation"]["Article"]["AuthorList"][4]["LastName"] == "Aschebrook-Kilfoy"
        assert article["MedlineCitation"]["Article"]["AuthorList"][4]["Initials"] == "B"
        assert article["MedlineCitation"]["Article"]["AuthorList"][4]["Identifier"] == []
        assert len(
                article["MedlineCitation"]["Article"]["AuthorList"][4][
                    "AffiliationInfo"
                ]
            ) == 1
        assert article["MedlineCitation"]["Article"]["AuthorList"][4]["AffiliationInfo"][
                0
            ]["Affiliation"] == "Department of Public Health Sciences, The University of Chicago, Chicago, Illinois, USA."
        assert article["MedlineCitation"]["Article"]["AuthorList"][4]["AffiliationInfo"][
                0
            ]["Identifier"] == []
        assert len(
                article["MedlineCitation"]["Article"]["AuthorList"][4][
                    "AffiliationInfo"
                ][0].attributes
            ) == 0
        assert article["MedlineCitation"]["Article"]["AuthorList"][4]["ForeName"] == "Briseis"
        assert article["MedlineCitation"]["Article"]["AuthorList"][4].attributes == {"ValidYN": "Y"}
        assert article["MedlineCitation"]["Article"]["AuthorList"][5]["LastName"] == "Jasmine"
        assert article["MedlineCitation"]["Article"]["AuthorList"][5]["Initials"] == "F"
        assert article["MedlineCitation"]["Article"]["AuthorList"][5]["Identifier"] == []
        assert len(
                article["MedlineCitation"]["Article"]["AuthorList"][5][
                    "AffiliationInfo"
                ]
            ) == 1
        assert article["MedlineCitation"]["Article"]["AuthorList"][5]["AffiliationInfo"][
                0
            ]["Affiliation"] == "Department of Public Health Sciences, The University of Chicago, Chicago, Illinois, USA."
        assert article["MedlineCitation"]["Article"]["AuthorList"][5]["AffiliationInfo"][
                0
            ]["Identifier"] == []
        assert len(
                article["MedlineCitation"]["Article"]["AuthorList"][5][
                    "AffiliationInfo"
                ][0].attributes
            ) == 0
        assert article["MedlineCitation"]["Article"]["AuthorList"][5]["ForeName"] == "Farzana"
        assert article["MedlineCitation"]["Article"]["AuthorList"][5].attributes == {"ValidYN": "Y"}
        assert article["MedlineCitation"]["Article"]["AuthorList"][6]["LastName"] == "Koutros"
        assert article["MedlineCitation"]["Article"]["AuthorList"][6]["Initials"] == "S"
        assert article["MedlineCitation"]["Article"]["AuthorList"][6]["Identifier"] == []
        assert len(
                article["MedlineCitation"]["Article"]["AuthorList"][6][
                    "AffiliationInfo"
                ]
            ) == 1
        assert article["MedlineCitation"]["Article"]["AuthorList"][6]["AffiliationInfo"][
                0
            ]["Affiliation"] == "Division of Cancer Epidemiology and Genetics, National Cancer Institute, Rockville, Maryland, USA."
        assert article["MedlineCitation"]["Article"]["AuthorList"][6]["AffiliationInfo"][
                0
            ]["Identifier"] == []
        assert len(
                article["MedlineCitation"]["Article"]["AuthorList"][6][
                    "AffiliationInfo"
                ][0].attributes
            ) == 0
        assert article["MedlineCitation"]["Article"]["AuthorList"][6]["ForeName"] == "Stella"
        assert article["MedlineCitation"]["Article"]["AuthorList"][6].attributes == {"ValidYN": "Y"}
        assert article["MedlineCitation"]["Article"]["AuthorList"][7]["LastName"] == "Parks"
        assert article["MedlineCitation"]["Article"]["AuthorList"][7]["Initials"] == "CG"
        assert article["MedlineCitation"]["Article"]["AuthorList"][7]["Identifier"] == []
        assert len(
                article["MedlineCitation"]["Article"]["AuthorList"][7][
                    "AffiliationInfo"
                ]
            ) == 1
        assert article["MedlineCitation"]["Article"]["AuthorList"][7]["AffiliationInfo"][
                0
            ]["Affiliation"] == "National Institute of Environmental Health Sciences, Research Triangle Park, North Carolina, USA."
        assert article["MedlineCitation"]["Article"]["AuthorList"][7]["AffiliationInfo"][
                0
            ]["Identifier"] == []
        assert len(
                article["MedlineCitation"]["Article"]["AuthorList"][7][
                    "AffiliationInfo"
                ][0].attributes
            ) == 0
        assert article["MedlineCitation"]["Article"]["AuthorList"][7]["ForeName"] == "Christine G"
        assert article["MedlineCitation"]["Article"]["AuthorList"][7].attributes == {"ValidYN": "Y"}
        assert article["MedlineCitation"]["Article"]["AuthorList"][8]["LastName"] == "Sandler"
        assert article["MedlineCitation"]["Article"]["AuthorList"][8]["Initials"] == "DP"
        assert article["MedlineCitation"]["Article"]["AuthorList"][8]["Identifier"] == []
        assert len(
                article["MedlineCitation"]["Article"]["AuthorList"][8][
                    "AffiliationInfo"
                ]
            ) == 1
        assert article["MedlineCitation"]["Article"]["AuthorList"][8]["AffiliationInfo"][
                0
            ]["Affiliation"] == "National Institute of Environmental Health Sciences, Research Triangle Park, North Carolina, USA."
        assert article["MedlineCitation"]["Article"]["AuthorList"][8]["AffiliationInfo"][
                0
            ]["Identifier"] == []
        assert len(
                article["MedlineCitation"]["Article"]["AuthorList"][8][
                    "AffiliationInfo"
                ][0].attributes
            ) == 0
        assert article["MedlineCitation"]["Article"]["AuthorList"][8]["ForeName"] == "Dale P"
        assert article["MedlineCitation"]["Article"]["AuthorList"][8].attributes == {"ValidYN": "Y"}
        assert article["MedlineCitation"]["Article"]["AuthorList"][9]["LastName"] == "Alavanja"
        assert article["MedlineCitation"]["Article"]["AuthorList"][9]["Initials"] == "MCR"
        assert article["MedlineCitation"]["Article"]["AuthorList"][9]["Identifier"] == []
        assert len(
                article["MedlineCitation"]["Article"]["AuthorList"][9][
                    "AffiliationInfo"
                ]
            ) == 2
        assert article["MedlineCitation"]["Article"]["AuthorList"][9]["AffiliationInfo"][
                0
            ]["Affiliation"] == "Division of Cancer Epidemiology and Genetics, National Cancer Institute, Rockville, Maryland, USA."
        assert article["MedlineCitation"]["Article"]["AuthorList"][9]["AffiliationInfo"][
                0
            ]["Identifier"] == []
        assert len(
                article["MedlineCitation"]["Article"]["AuthorList"][9][
                    "AffiliationInfo"
                ][0].attributes
            ) == 0
        assert article["MedlineCitation"]["Article"]["AuthorList"][9]["AffiliationInfo"][
                1
            ]["Affiliation"] == "Department of Biology, Hood College, Frederick, Maryland, USA."
        assert article["MedlineCitation"]["Article"]["AuthorList"][9]["AffiliationInfo"][
                1
            ]["Identifier"] == []
        assert len(
                article["MedlineCitation"]["Article"]["AuthorList"][9][
                    "AffiliationInfo"
                ][1].attributes
            ) == 0
        assert article["MedlineCitation"]["Article"]["AuthorList"][9]["ForeName"] == "Michael C R"
        assert article["MedlineCitation"]["Article"]["AuthorList"][9].attributes == {"ValidYN": "Y"}
        assert article["MedlineCitation"]["Article"]["AuthorList"][10]["LastName"] == "Hofmann"
        assert article["MedlineCitation"]["Article"]["AuthorList"][10]["Initials"] == "JN"
        assert article["MedlineCitation"]["Article"]["AuthorList"][10]["Identifier"] == []
        assert len(
                article["MedlineCitation"]["Article"]["AuthorList"][10][
                    "AffiliationInfo"
                ]
            ) == 1
        assert article["MedlineCitation"]["Article"]["AuthorList"][10]["AffiliationInfo"][
                0
            ]["Affiliation"] == "Division of Cancer Epidemiology and Genetics, National Cancer Institute, Rockville, Maryland, USA."
        assert article["MedlineCitation"]["Article"]["AuthorList"][10]["AffiliationInfo"][
                0
            ]["Identifier"] == []
        assert len(
                article["MedlineCitation"]["Article"]["AuthorList"][10][
                    "AffiliationInfo"
                ][0].attributes
            ) == 0
        assert article["MedlineCitation"]["Article"]["AuthorList"][10]["ForeName"] == "Jonathan N"
        assert article["MedlineCitation"]["Article"]["AuthorList"][10].attributes == {"ValidYN": "Y"}
        assert article["MedlineCitation"]["Article"]["AuthorList"][11]["LastName"] == "Ward"
        assert article["MedlineCitation"]["Article"]["AuthorList"][11]["Initials"] == "MH"
        assert article["MedlineCitation"]["Article"]["AuthorList"][11]["Identifier"] == []
        assert len(
                article["MedlineCitation"]["Article"]["AuthorList"][11][
                    "AffiliationInfo"
                ]
            ) == 1
        assert article["MedlineCitation"]["Article"]["AuthorList"][11]["AffiliationInfo"][
                0
            ]["Affiliation"] == "Division of Cancer Epidemiology and Genetics, National Cancer Institute, Rockville, Maryland, USA."
        assert article["MedlineCitation"]["Article"]["AuthorList"][11]["AffiliationInfo"][
                0
            ]["Identifier"] == []
        assert len(
                article["MedlineCitation"]["Article"]["AuthorList"][11][
                    "AffiliationInfo"
                ][0].attributes
            ) == 0
        assert article["MedlineCitation"]["Article"]["AuthorList"][11]["ForeName"] == "Mary H"
        assert article["MedlineCitation"]["Article"]["AuthorList"][11].attributes == {"ValidYN": "Y"}
        assert len(article["MedlineCitation"]["Article"]["Language"]) == 1
        assert article["MedlineCitation"]["Article"]["Language"][0] == "eng"
        assert len(article["MedlineCitation"]["Article"]["PublicationTypeList"]) == 1
        assert article["MedlineCitation"]["Article"]["PublicationTypeList"][0] == "Journal Article"
        assert article["MedlineCitation"]["Article"]["PublicationTypeList"][0].attributes == {"UI": "D016428"}
        assert len(article["MedlineCitation"]["Article"]["Journal"]) == 4
        assert article["MedlineCitation"]["Article"]["Journal"]["ISSN"] == "1470-7926"
        assert article["MedlineCitation"]["Article"]["Journal"]["ISSN"].attributes == {"IssnType": "Electronic"}
        assert article["MedlineCitation"]["Article"]["Journal"]["ISOAbbreviation"] == "Occup Environ Med"
        assert len(article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]) == 3
        assert article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["Volume"] == "75"
        assert article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["Issue"] == "2"
        assert len(
                article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                    "PubDate"
                ]
            ) == 2
        assert article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"][
                "Month"
            ] == "Feb"
        assert article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"][
                "Year"
            ] == "2018"
        assert len(
                article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                    "PubDate"
                ].attributes
            ) == 0
        assert article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"].attributes == {"CitedMedium": "Internet"}
        assert article["MedlineCitation"]["Article"]["Journal"]["Title"] == "Occupational and environmental medicine"
        assert article["MedlineCitation"]["Article"]["Journal"]["ISSN"].attributes == {"IssnType": "Electronic"}
        assert article["MedlineCitation"]["Article"]["ArticleTitle"] == "Occupational pesticide exposure and subclinical hypothyroidism among male pesticide applicators."
        assert len(article["MedlineCitation"]["Article"]["ELocationID"]) == 1
        assert article["MedlineCitation"]["Article"]["ELocationID"][0] == "10.1136/oemed-2017-104431"
        assert len(article["MedlineCitation"]["Article"]["ELocationID"][0].attributes) == 2
        assert article["MedlineCitation"]["Article"]["ELocationID"][0].attributes[
                "ValidYN"
            ] == "Y"
        assert article["MedlineCitation"]["Article"]["ELocationID"][0].attributes[
                "EIdType"
            ] == "doi"
        assert len(article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]) == 4
        assert article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][0] == "Animal studies suggest that exposure to pesticides may alter thyroid function; however, few epidemiologic studies have examined this association. We evaluated the relationship between individual pesticides and thyroid function in 679 men enrolled in a substudy of the Agricultural Health Study, a cohort of licensed pesticide applicators."
        assert article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][
                0
            ].attributes == {"NlmCategory": "OBJECTIVE", "Label": "OBJECTIVES"}
        assert article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][1] == "Self-reported lifetime pesticide use was obtained at cohort enrolment (1993-1997). Intensity-weighted lifetime days were computed for 33 pesticides, which adjusts cumulative days of pesticide use for factors that modify exposure (eg, use of personal protective equipment). Thyroid-stimulating hormone (TSH), thyroxine (T4), triiodothyronine (T3) and antithyroid peroxidase (anti-TPO) autoantibodies were measured in serum collected in 2010-2013. We used multivariate logistic regression to estimate ORs and 95% CIs for subclinical hypothyroidism (TSH &gt;4.5 mIU/L) compared with normal TSH (0.4-<u>&lt;</u>4.5 mIU/L) and for anti-TPO positivity. We also examined pesticide associations with TSH, T4 and T3 in multivariate linear regression models."
        assert article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][
                1
            ].attributes == {"NlmCategory": "METHODS", "Label": "METHODS"}
        assert article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][2] == "Higher exposure to the insecticide aldrin (third and fourth quartiles of intensity-weighted days vs no exposure) was positively associated with subclinical hypothyroidism (OR<sub>Q3</sub>=4.15, 95% CI 1.56 to 11.01, OR<sub>Q4</sub>=4.76, 95% CI 1.53 to 14.82, p<sub>trend</sub> &lt;0.01), higher TSH (p<sub>trend</sub>=0.01) and lower T4 (p<sub>trend</sub>=0.04). Higher exposure to the herbicide pendimethalin was associated with subclinical hypothyroidism (fourth quartile vs no exposure: OR<sub>Q4</sub>=2.78, 95% CI 1.30 to 5.95, p<sub>trend</sub>=0.02), higher TSH (p<sub>trend</sub>=0.04) and anti-TPO positivity (p<sub>trend</sub>=0.01). The fumigant methyl bromide was inversely associated with TSH (p<sub>trend</sub>=0.02) and positively associated with T4 (p<sub>trend</sub>=0.01)."
        assert article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][
                2
            ].attributes == {"NlmCategory": "RESULTS", "Label": "RESULTS"}
        assert article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][3] == "Our results suggest that long-term exposure to aldrin, pendimethalin and methyl bromide may alter thyroid function among male pesticide applicators."
        assert article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][
                3
            ].attributes == {"NlmCategory": "CONCLUSIONS", "Label": "CONCLUSIONS"}
        assert article["MedlineCitation"]["Article"]["Abstract"]["CopyrightInformation"] == "\xa9 Article author(s) (or their employer(s) unless otherwise stated in the text of the article) 2018. All rights reserved. No commercial use is permitted unless otherwise expressly granted."
        assert len(article["MedlineCitation"]["Article"]["GrantList"]) == 3
        assert len(article["MedlineCitation"]["Article"]["GrantList"][0]) == 4
        assert article["MedlineCitation"]["Article"]["GrantList"][0]["Acronym"] == "CP"
        assert article["MedlineCitation"]["Article"]["GrantList"][0]["Country"] == "United States"
        assert article["MedlineCitation"]["Article"]["GrantList"][0]["Agency"] == "NCI NIH HHS"
        assert article["MedlineCitation"]["Article"]["GrantList"][0]["GrantID"] == "Z01 CP010119"
        assert len(article["MedlineCitation"]["Article"]["GrantList"][0].attributes) == 0
        assert len(article["MedlineCitation"]["Article"]["GrantList"][1]) == 4
        assert article["MedlineCitation"]["Article"]["GrantList"][1]["Acronym"] == "ES"
        assert article["MedlineCitation"]["Article"]["GrantList"][1]["Country"] == "United States"
        assert article["MedlineCitation"]["Article"]["GrantList"][1]["Agency"] == "NIEHS NIH HHS"
        assert article["MedlineCitation"]["Article"]["GrantList"][1]["GrantID"] == "Z01 ES049030"
        assert len(article["MedlineCitation"]["Article"]["GrantList"][1].attributes) == 0
        assert len(article["MedlineCitation"]["Article"]["GrantList"][2]) == 4
        assert article["MedlineCitation"]["Article"]["GrantList"][2]["Acronym"] == "NULL"
        assert article["MedlineCitation"]["Article"]["GrantList"][2]["Country"] == "United States"
        assert article["MedlineCitation"]["Article"]["GrantList"][2]["Agency"] == "Intramural NIH HHS"
        assert article["MedlineCitation"]["Article"]["GrantList"][2]["GrantID"] == "Z99 CA999999"
        assert len(article["MedlineCitation"]["Article"]["GrantList"][2].attributes) == 0
        assert article["MedlineCitation"]["Article"]["GrantList"].attributes == {"CompleteYN": "Y"}

    def test_pubmed_html_mathml_tags(self):
        """Test parsing XML returned by EFetch, PubMed database, with both HTML and MathML tags."""
        # In PubMed display PMID 30108519 in xml retrieval mode, containing
        # both HTML and MathML tags in the abstract text.
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db="pubmed", id='30108519', rettype="null",
        #                       retmode="xml", parsed=True)
        with open("Entrez/pubmed6.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert len(record) == 2
        assert record["PubmedBookArticle"] == []
        assert len(record["PubmedArticle"]) == 1
        pubmed_article = record["PubmedArticle"][0]
        assert len(pubmed_article) == 2
        assert len(pubmed_article["PubmedData"]) == 3
        assert len(pubmed_article["PubmedData"]["ArticleIdList"]) == 3
        assert pubmed_article["PubmedData"]["ArticleIdList"][0] == "30108519"
        assert pubmed_article["PubmedData"]["ArticleIdList"][0].attributes == {"IdType": "pubmed"}
        assert pubmed_article["PubmedData"]["ArticleIdList"][1] == "10.3389/fphys.2018.01034"
        assert pubmed_article["PubmedData"]["ArticleIdList"][1].attributes == {"IdType": "doi"}
        assert pubmed_article["PubmedData"]["ArticleIdList"][2] == "PMC6079548"
        assert pubmed_article["PubmedData"]["ArticleIdList"][2].attributes == {"IdType": "pmc"}
        assert pubmed_article["PubmedData"]["PublicationStatus"] == "epublish"
        assert len(pubmed_article["PubmedData"]["History"]) == 5
        assert len(pubmed_article["PubmedData"]["History"][0]) == 3
        assert pubmed_article["PubmedData"]["History"][0]["Year"] == "2018"
        assert pubmed_article["PubmedData"]["History"][0]["Month"] == "05"
        assert pubmed_article["PubmedData"]["History"][0]["Day"] == "22"
        assert pubmed_article["PubmedData"]["History"][0].attributes == {"PubStatus": "received"}
        assert len(pubmed_article["PubmedData"]["History"][1]) == 3
        assert pubmed_article["PubmedData"]["History"][1]["Year"] == "2018"
        assert pubmed_article["PubmedData"]["History"][1]["Month"] == "07"
        assert pubmed_article["PubmedData"]["History"][1]["Day"] == "11"
        assert pubmed_article["PubmedData"]["History"][1].attributes == {"PubStatus": "accepted"}
        assert len(pubmed_article["PubmedData"]["History"][2]) == 5
        assert pubmed_article["PubmedData"]["History"][2]["Year"] == "2018"
        assert pubmed_article["PubmedData"]["History"][2]["Month"] == "8"
        assert pubmed_article["PubmedData"]["History"][2]["Day"] == "16"
        assert pubmed_article["PubmedData"]["History"][2]["Hour"] == "6"
        assert pubmed_article["PubmedData"]["History"][2]["Minute"] == "0"
        assert pubmed_article["PubmedData"]["History"][2].attributes == {"PubStatus": "entrez"}
        assert len(pubmed_article["PubmedData"]["History"][3]) == 5
        assert pubmed_article["PubmedData"]["History"][3]["Year"] == "2018"
        assert pubmed_article["PubmedData"]["History"][3]["Month"] == "8"
        assert pubmed_article["PubmedData"]["History"][3]["Day"] == "16"
        assert pubmed_article["PubmedData"]["History"][3]["Hour"] == "6"
        assert pubmed_article["PubmedData"]["History"][3]["Minute"] == "0"
        assert pubmed_article["PubmedData"]["History"][3].attributes == {"PubStatus": "pubmed"}
        assert len(pubmed_article["PubmedData"]["History"][4]) == 5
        assert pubmed_article["PubmedData"]["History"][4]["Year"] == "2018"
        assert pubmed_article["PubmedData"]["History"][4]["Month"] == "8"
        assert pubmed_article["PubmedData"]["History"][4]["Day"] == "16"
        assert pubmed_article["PubmedData"]["History"][4]["Hour"] == "6"
        assert pubmed_article["PubmedData"]["History"][4]["Minute"] == "1"
        assert pubmed_article["PubmedData"]["History"][4].attributes == {"PubStatus": "medline"}
        medline_citation = pubmed_article["MedlineCitation"]
        assert len(medline_citation) == 11
        assert medline_citation["GeneralNote"] == []
        assert len(medline_citation["KeywordList"]) == 1
        assert len(medline_citation["KeywordList"][0]) == 8
        assert medline_citation["KeywordList"][0][0] == "Owles' point"
        assert medline_citation["KeywordList"][0][0].attributes == {"MajorTopicYN": "N"}
        assert medline_citation["KeywordList"][0][1] == "aerobic capacity"
        assert medline_citation["KeywordList"][0][1].attributes == {"MajorTopicYN": "N"}
        assert medline_citation["KeywordList"][0][2] == "aerobic threshold"
        assert medline_citation["KeywordList"][0][2].attributes == {"MajorTopicYN": "N"}
        assert medline_citation["KeywordList"][0][3] == "anaerobic threshold"
        assert medline_citation["KeywordList"][0][3].attributes == {"MajorTopicYN": "N"}
        assert medline_citation["KeywordList"][0][4] == "endurance assessment"
        assert medline_citation["KeywordList"][0][4].attributes == {"MajorTopicYN": "N"}
        assert medline_citation["KeywordList"][0][5] == "lactate threshold"
        assert medline_citation["KeywordList"][0][5].attributes == {"MajorTopicYN": "N"}
        assert medline_citation["KeywordList"][0][6] == "oxygen endurance performance limit"
        assert medline_citation["KeywordList"][0][6].attributes == {"MajorTopicYN": "N"}
        assert medline_citation["KeywordList"][0][7] == "submaximal exercise testing"
        assert medline_citation["KeywordList"][0][7].attributes == {"MajorTopicYN": "N"}
        assert medline_citation["CitationSubset"] == []
        assert medline_citation["OtherAbstract"] == []
        assert medline_citation["OtherID"] == []
        assert medline_citation["SpaceFlightMission"] == []
        assert medline_citation["PMID"] == "30108519"
        assert medline_citation["PMID"].attributes == {"Version": "1"}
        assert len(medline_citation["DateRevised"]) == 3
        assert medline_citation["DateRevised"]["Year"] == "2018"
        assert medline_citation["DateRevised"]["Month"] == "08"
        assert medline_citation["DateRevised"]["Day"] == "17"
        assert medline_citation["DateRevised"].attributes == {}
        assert len(medline_citation["MedlineJournalInfo"]) == 4
        assert medline_citation["MedlineJournalInfo"]["Country"] == "Switzerland"
        assert medline_citation["MedlineJournalInfo"]["MedlineTA"] == "Front Physiol"
        assert medline_citation["MedlineJournalInfo"]["NlmUniqueID"] == "101549006"
        assert medline_citation["MedlineJournalInfo"]["ISSNLinking"] == "1664-042X"
        assert medline_citation["MedlineJournalInfo"].attributes == {}
        assert len(medline_citation["CommentsCorrectionsList"]) == 53
        assert len(medline_citation["CommentsCorrectionsList"][0]) == 2
        assert medline_citation["CommentsCorrectionsList"][0]["RefSource"] == "Stat Med. 2008 Feb 28;27(5):778-80"
        assert medline_citation["CommentsCorrectionsList"][0]["PMID"] == "17907247"
        assert medline_citation["CommentsCorrectionsList"][0]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][0].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][1]) == 2
        assert medline_citation["CommentsCorrectionsList"][1]["RefSource"] == "Int J Sports Med. 2009 Jan;30(1):40-45"
        assert medline_citation["CommentsCorrectionsList"][1]["PMID"] == "19202577"
        assert medline_citation["CommentsCorrectionsList"][1]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][1].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][2]) == 2
        assert medline_citation["CommentsCorrectionsList"][2]["RefSource"] == "Med Sci Sports Exerc. 1995 Jun;27(6):863-7"
        assert medline_citation["CommentsCorrectionsList"][2]["PMID"] == "7658947"
        assert medline_citation["CommentsCorrectionsList"][2]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][2].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][3]) == 2
        assert medline_citation["CommentsCorrectionsList"][3]["RefSource"] == "Eur J Appl Physiol. 2010 Apr;108(6):1153-67"
        assert medline_citation["CommentsCorrectionsList"][3]["PMID"] == "20033207"
        assert medline_citation["CommentsCorrectionsList"][3]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][3].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][4]) == 2
        assert medline_citation["CommentsCorrectionsList"][4]["RefSource"] == "Med Sci Sports Exerc. 1999 Apr;31(4):578-82"
        assert medline_citation["CommentsCorrectionsList"][4]["PMID"] == "10211855"
        assert medline_citation["CommentsCorrectionsList"][4]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][4].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][5]) == 2
        assert medline_citation["CommentsCorrectionsList"][5]["RefSource"] == "Br J Sports Med. 1988 Jun;22(2):51-4"
        assert medline_citation["CommentsCorrectionsList"][5]["PMID"] == "3167501"
        assert medline_citation["CommentsCorrectionsList"][5]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][5].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][6]) == 2
        assert medline_citation["CommentsCorrectionsList"][6]["RefSource"] == "Front Physiol. 2017 Jun 08;8:389"
        assert medline_citation["CommentsCorrectionsList"][6]["PMID"] == "28642717"
        assert medline_citation["CommentsCorrectionsList"][6]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][6].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][7]) == 2
        assert medline_citation["CommentsCorrectionsList"][7]["RefSource"] == "Med Sci Sports Exerc. 1999 Sep;31(9):1342-5"
        assert medline_citation["CommentsCorrectionsList"][7]["PMID"] == "10487378"
        assert medline_citation["CommentsCorrectionsList"][7]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][7].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][8]) == 2
        assert medline_citation["CommentsCorrectionsList"][8]["RefSource"] == "Med Sci Sports Exerc. 1998 Aug;30(8):1304-13"
        assert medline_citation["CommentsCorrectionsList"][8]["PMID"] == "9710874"
        assert medline_citation["CommentsCorrectionsList"][8]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][8].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][9]) == 2
        assert medline_citation["CommentsCorrectionsList"][9]["RefSource"] == "Med Sci Sports. 1979 Winter;11(4):338-44"
        assert medline_citation["CommentsCorrectionsList"][9]["PMID"] == "530025"
        assert medline_citation["CommentsCorrectionsList"][9]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][9].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][10]) == 2
        assert medline_citation["CommentsCorrectionsList"][10]["RefSource"] == "J Strength Cond Res. 2005 May;19(2):364-8"
        assert medline_citation["CommentsCorrectionsList"][10]["PMID"] == "15903376"
        assert medline_citation["CommentsCorrectionsList"][10]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][10].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][11]) == 2
        assert medline_citation["CommentsCorrectionsList"][11]["RefSource"] == "Eur J Appl Physiol Occup Physiol. 1984;53(3):196-9"
        assert medline_citation["CommentsCorrectionsList"][11]["PMID"] == "6542852"
        assert medline_citation["CommentsCorrectionsList"][11]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][11].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][12]) == 2
        assert medline_citation["CommentsCorrectionsList"][12]["RefSource"] == "Eur J Appl Physiol Occup Physiol. 1978 Oct 20;39(4):219-27"
        assert medline_citation["CommentsCorrectionsList"][12]["PMID"] == "710387"
        assert medline_citation["CommentsCorrectionsList"][12]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][12].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][13]) == 2
        assert medline_citation["CommentsCorrectionsList"][13]["RefSource"] == "J Appl Physiol Respir Environ Exerc Physiol. 1980 Mar;48(3):523-7"
        assert medline_citation["CommentsCorrectionsList"][13]["PMID"] == "7372524"
        assert medline_citation["CommentsCorrectionsList"][13]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][13].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][14]) == 2
        assert medline_citation["CommentsCorrectionsList"][14]["RefSource"] == "Int J Sports Med. 2015 Dec;36(14):1142-8"
        assert medline_citation["CommentsCorrectionsList"][14]["PMID"] == "26332904"
        assert medline_citation["CommentsCorrectionsList"][14]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][14].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][15]) == 2
        assert medline_citation["CommentsCorrectionsList"][15]["RefSource"] == "J Physiol. 1930 Apr 14;69(2):214-37"
        assert medline_citation["CommentsCorrectionsList"][15]["PMID"] == "16994099"
        assert medline_citation["CommentsCorrectionsList"][15]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][15].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][16]) == 2
        assert medline_citation["CommentsCorrectionsList"][16]["RefSource"] == "J Strength Cond Res. 2015 Oct;29(10):2794-801"
        assert medline_citation["CommentsCorrectionsList"][16]["PMID"] == "25844867"
        assert medline_citation["CommentsCorrectionsList"][16]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][16].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][17]) == 2
        assert medline_citation["CommentsCorrectionsList"][17]["RefSource"] == "PLoS One. 2018 Mar 13;13(3):e0194313"
        assert medline_citation["CommentsCorrectionsList"][17]["PMID"] == "29534108"
        assert medline_citation["CommentsCorrectionsList"][17]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][17].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][18]) == 2
        assert medline_citation["CommentsCorrectionsList"][18]["RefSource"] == "J Cardiopulm Rehabil Prev. 2012 Nov-Dec;32(6):327-50"
        assert medline_citation["CommentsCorrectionsList"][18]["PMID"] == "23103476"
        assert medline_citation["CommentsCorrectionsList"][18]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][18].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][19]) == 2
        assert medline_citation["CommentsCorrectionsList"][19]["RefSource"] == "Exerc Sport Sci Rev. 1982;10:49-83"
        assert medline_citation["CommentsCorrectionsList"][19]["PMID"] == "6811284"
        assert medline_citation["CommentsCorrectionsList"][19]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][19].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][20]) == 2
        assert medline_citation["CommentsCorrectionsList"][20]["RefSource"] == "Int J Sports Physiol Perform. 2010 Sep;5(3):276-91"
        assert medline_citation["CommentsCorrectionsList"][20]["PMID"] == "20861519"
        assert medline_citation["CommentsCorrectionsList"][20]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][20].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][21]) == 2
        assert medline_citation["CommentsCorrectionsList"][21]["RefSource"] == "Eur J Appl Physiol Occup Physiol. 1990;60(4):249-53"
        assert medline_citation["CommentsCorrectionsList"][21]["PMID"] == "2357979"
        assert medline_citation["CommentsCorrectionsList"][21]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][21].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][22]) == 2
        assert medline_citation["CommentsCorrectionsList"][22]["RefSource"] == "Med Sci Sports Exerc. 2004 Oct;36(10):1737-42"
        assert medline_citation["CommentsCorrectionsList"][22]["PMID"] == "15595295"
        assert medline_citation["CommentsCorrectionsList"][22]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][22].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][23]) == 2
        assert medline_citation["CommentsCorrectionsList"][23]["RefSource"] == "Int J Sports Med. 2016 Jun;37(7):539-46"
        assert medline_citation["CommentsCorrectionsList"][23]["PMID"] == "27116348"
        assert medline_citation["CommentsCorrectionsList"][23]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][23].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][24]) == 2
        assert medline_citation["CommentsCorrectionsList"][24]["RefSource"] == "Scand J Med Sci Sports. 2017 May;27(5):462-473"
        assert medline_citation["CommentsCorrectionsList"][24]["PMID"] == "28181710"
        assert medline_citation["CommentsCorrectionsList"][24]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][24].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][25]) == 2
        assert medline_citation["CommentsCorrectionsList"][25]["RefSource"] == "Int J Sports Med. 1983 Nov;4(4):226-30"
        assert medline_citation["CommentsCorrectionsList"][25]["PMID"] == "6654546"
        assert medline_citation["CommentsCorrectionsList"][25]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][25].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][26]) == 2
        assert medline_citation["CommentsCorrectionsList"][26]["RefSource"] == "J Appl Physiol (1985). 1988 Jun;64(6):2622-30"
        assert medline_citation["CommentsCorrectionsList"][26]["PMID"] == "3403447"
        assert medline_citation["CommentsCorrectionsList"][26]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][26].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][27]) == 2
        assert medline_citation["CommentsCorrectionsList"][27]["RefSource"] == "Med Sci Sports Exerc. 2009 Jan;41(1):3-13"
        assert medline_citation["CommentsCorrectionsList"][27]["PMID"] == "19092709"
        assert medline_citation["CommentsCorrectionsList"][27]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][27].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][28]) == 2
        assert medline_citation["CommentsCorrectionsList"][28]["RefSource"] == "Int J Sports Med. 2009 Sep;30(9):643-6"
        assert medline_citation["CommentsCorrectionsList"][28]["PMID"] == "19569005"
        assert medline_citation["CommentsCorrectionsList"][28]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][28].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][29]) == 2
        assert medline_citation["CommentsCorrectionsList"][29]["RefSource"] == "Eur J Appl Physiol Occup Physiol. 1988;57(4):420-4"
        assert medline_citation["CommentsCorrectionsList"][29]["PMID"] == "3396556"
        assert medline_citation["CommentsCorrectionsList"][29]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][29].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][30]) == 2
        assert medline_citation["CommentsCorrectionsList"][30]["RefSource"] == "J Physiol. 2004 Jul 1;558(Pt 1):5-30"
        assert medline_citation["CommentsCorrectionsList"][30]["PMID"] == "15131240"
        assert medline_citation["CommentsCorrectionsList"][30]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][30].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][31]) == 2
        assert medline_citation["CommentsCorrectionsList"][31]["RefSource"] == "Int J Sports Med. 1990 Feb;11(1):26-32"
        assert medline_citation["CommentsCorrectionsList"][31]["PMID"] == "2318561"
        assert medline_citation["CommentsCorrectionsList"][31]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][31].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][32]) == 2
        assert medline_citation["CommentsCorrectionsList"][32]["RefSource"] == "J Appl Physiol. 1973 Aug;35(2):236-43"
        assert medline_citation["CommentsCorrectionsList"][32]["PMID"] == "4723033"
        assert medline_citation["CommentsCorrectionsList"][32]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][32].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][33]) == 2
        assert medline_citation["CommentsCorrectionsList"][33]["RefSource"] == "Int J Sports Med. 1987 Dec;8(6):401-6"
        assert medline_citation["CommentsCorrectionsList"][33]["PMID"] == "3429086"
        assert medline_citation["CommentsCorrectionsList"][33]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][33].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][34]) == 2
        assert medline_citation["CommentsCorrectionsList"][34]["RefSource"] == "J Sci Med Sport. 2008 Jun;11(3):280-6"
        assert medline_citation["CommentsCorrectionsList"][34]["PMID"] == "17553745"
        assert medline_citation["CommentsCorrectionsList"][34]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][34].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][35]) == 2
        assert medline_citation["CommentsCorrectionsList"][35]["RefSource"] == "J Appl Physiol Respir Environ Exerc Physiol. 1984 May;56(5):1260-4"
        assert medline_citation["CommentsCorrectionsList"][35]["PMID"] == "6725086"
        assert medline_citation["CommentsCorrectionsList"][35]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][35].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][36]) == 2
        assert medline_citation["CommentsCorrectionsList"][36]["RefSource"] == "Int J Sports Med. 2008 Jun;29(6):475-9"
        assert medline_citation["CommentsCorrectionsList"][36]["PMID"] == "18302077"
        assert medline_citation["CommentsCorrectionsList"][36]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][36].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][37]) == 2
        assert medline_citation["CommentsCorrectionsList"][37]["RefSource"] == "Med Sci Sports Exerc. 1985 Feb;17(1):22-34"
        assert medline_citation["CommentsCorrectionsList"][37]["PMID"] == "3884959"
        assert medline_citation["CommentsCorrectionsList"][37]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][37].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][38]) == 2
        assert medline_citation["CommentsCorrectionsList"][38]["RefSource"] == "Sports Med. 2009;39(6):469-90"
        assert medline_citation["CommentsCorrectionsList"][38]["PMID"] == "19453206"
        assert medline_citation["CommentsCorrectionsList"][38]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][38].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][39]) == 2
        assert medline_citation["CommentsCorrectionsList"][39]["RefSource"] == "Int J Sports Med. 2004 Aug;25(6):403-8"
        assert medline_citation["CommentsCorrectionsList"][39]["PMID"] == "15346226"
        assert medline_citation["CommentsCorrectionsList"][39]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][39].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][40]) == 2
        assert medline_citation["CommentsCorrectionsList"][40]["RefSource"] == "J Sports Med Phys Fitness. 2004 Jun;44(2):132-40"
        assert medline_citation["CommentsCorrectionsList"][40]["PMID"] == "15470310"
        assert medline_citation["CommentsCorrectionsList"][40]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][40].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][41]) == 2
        assert medline_citation["CommentsCorrectionsList"][41]["RefSource"] == "Int J Sports Med. 1985 Jun;6(3):117-30"
        assert medline_citation["CommentsCorrectionsList"][41]["PMID"] == "4030186"
        assert medline_citation["CommentsCorrectionsList"][41]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][41].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][42]) == 2
        assert medline_citation["CommentsCorrectionsList"][42]["RefSource"] == "Int J Sports Med. 1999 Feb;20(2):122-7"
        assert medline_citation["CommentsCorrectionsList"][42]["PMID"] == "10190774"
        assert medline_citation["CommentsCorrectionsList"][42]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][42].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][43]) == 2
        assert medline_citation["CommentsCorrectionsList"][43]["RefSource"] == "Int J Sports Med. 2006 May;27(5):368-72"
        assert medline_citation["CommentsCorrectionsList"][43]["PMID"] == "16729378"
        assert medline_citation["CommentsCorrectionsList"][43]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][43].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][44]) == 2
        assert medline_citation["CommentsCorrectionsList"][44]["RefSource"] == "Int J Sports Med. 1985 Jun;6(3):109-16"
        assert medline_citation["CommentsCorrectionsList"][44]["PMID"] == "3897079"
        assert medline_citation["CommentsCorrectionsList"][44]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][44].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][45]) == 2
        assert medline_citation["CommentsCorrectionsList"][45]["RefSource"] == "Pneumologie. 1990 Jan;44(1):2-13"
        assert medline_citation["CommentsCorrectionsList"][45]["PMID"] == "2408033"
        assert medline_citation["CommentsCorrectionsList"][45]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][45].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][46]) == 2
        assert medline_citation["CommentsCorrectionsList"][46]["RefSource"] == "Eur J Appl Physiol. 2018 Apr;118(4):691-728"
        assert medline_citation["CommentsCorrectionsList"][46]["PMID"] == "29322250"
        assert medline_citation["CommentsCorrectionsList"][46]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][46].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][47]) == 2
        assert medline_citation["CommentsCorrectionsList"][47]["RefSource"] == "J Appl Physiol Respir Environ Exerc Physiol. 1983 Oct;55(4):1178-86"
        assert medline_citation["CommentsCorrectionsList"][47]["PMID"] == "6629951"
        assert medline_citation["CommentsCorrectionsList"][47]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][47].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][48]) == 2
        assert medline_citation["CommentsCorrectionsList"][48]["RefSource"] == "Sports Med. 2014 Nov;44 Suppl 2:S139-47"
        assert medline_citation["CommentsCorrectionsList"][48]["PMID"] == "25200666"
        assert medline_citation["CommentsCorrectionsList"][48]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][48].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][49]) == 2
        assert medline_citation["CommentsCorrectionsList"][49]["RefSource"] == "Front Physiol. 2015 Oct 30;6:308"
        assert medline_citation["CommentsCorrectionsList"][49]["PMID"] == "26578980"
        assert medline_citation["CommentsCorrectionsList"][49]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][49].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][50]) == 2
        assert medline_citation["CommentsCorrectionsList"][50]["RefSource"] == "Int J Sports Med. 2013 Mar;34(3):196-9"
        assert medline_citation["CommentsCorrectionsList"][50]["PMID"] == "22972242"
        assert medline_citation["CommentsCorrectionsList"][50]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][50].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][51]) == 2
        assert medline_citation["CommentsCorrectionsList"][51]["RefSource"] == "Int J Sports Med. 1992 Oct;13(7):518-22"
        assert medline_citation["CommentsCorrectionsList"][51]["PMID"] == "1459746"
        assert medline_citation["CommentsCorrectionsList"][51]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][51].attributes == {"RefType": "Cites"}
        assert len(medline_citation["CommentsCorrectionsList"][52]) == 2
        assert medline_citation["CommentsCorrectionsList"][52]["RefSource"] == "Med Sci Sports Exerc. 1993 May;25(5):620-7"
        assert medline_citation["CommentsCorrectionsList"][52]["PMID"] == "8492691"
        assert medline_citation["CommentsCorrectionsList"][52]["PMID"].attributes == {"Version": "1"}
        assert medline_citation["CommentsCorrectionsList"][52].attributes == {"RefType": "Cites"}
        article = medline_citation["Article"]
        assert len(article["ELocationID"]) == 1
        assert article["ELocationID"][0] == "10.3389/fphys.2018.01034"
        assert article["ELocationID"][0].attributes == {"EIdType": "doi", "ValidYN": "Y"}
        assert len(article["ArticleDate"]) == 1
        assert len(article["ArticleDate"][0]) == 3
        assert article["ArticleDate"][0]["Year"] == "2018"
        assert article["ArticleDate"][0]["Month"] == "07"
        assert article["ArticleDate"][0]["Day"] == "31"
        assert article["ArticleDate"][0].attributes == {"DateType": "Electronic"}
        assert article["Language"] == ["eng"]
        assert len(article["Journal"]) == 4
        assert article["Journal"]["ISSN"] == "1664-042X"
        assert article["Journal"]["ISSN"].attributes == {"IssnType": "Print"}
        assert article["Journal"]["JournalIssue"]["Volume"] == "9"
        assert article["Journal"]["JournalIssue"]["PubDate"]["Year"] == "2018"
        assert article["Journal"]["JournalIssue"]["PubDate"].attributes == {}
        assert article["Journal"]["JournalIssue"].attributes == {"CitedMedium": "Print"}
        assert article["Journal"]["Title"] == "Frontiers in physiology"
        assert article["Journal"]["ISOAbbreviation"] == "Front Physiol"
        assert article["Journal"].attributes == {}
        assert len(article["PublicationTypeList"]) == 1
        assert article["PublicationTypeList"][0] == "Journal Article"
        assert article["PublicationTypeList"][0].attributes == {"UI": "D016428"}
        assert article["ArticleTitle"] == 'A "<i>Blood Relationship"</i> Between the Overlooked Minimum Lactate Equivalent and Maximal Lactate Steady State in Trained Runners. Back to the Old Days?'
        assert len(article["Pagination"]) == 1
        assert article["Pagination"]["MedlinePgn"] == "1034"
        assert article["Pagination"].attributes == {}
        assert len(article["AuthorList"]) == 2
        assert len(article["AuthorList"][0]) == 5
        assert article["AuthorList"][0]["Identifier"] == []
        assert len(article["AuthorList"][0]["AffiliationInfo"]) == 1
        assert len(article["AuthorList"][0]["AffiliationInfo"][0]) == 2
        assert article["AuthorList"][0]["AffiliationInfo"][0]["Identifier"] == []
        assert article["AuthorList"][0]["AffiliationInfo"][0]["Affiliation"] == "Studies, Research and Sports Medicine Center, Government of Navarre, Pamplona, Spain."
        assert article["AuthorList"][0]["AffiliationInfo"][0].attributes == {}
        assert article["AuthorList"][0]["LastName"] == "Garcia-Tabar"
        assert article["AuthorList"][0]["ForeName"] == "Ibai"
        assert article["AuthorList"][0]["Initials"] == "I"
        assert article["AuthorList"][0].attributes == {"ValidYN": "Y"}
        assert len(article["AuthorList"][1]) == 5
        assert article["AuthorList"][1]["Identifier"] == []
        assert len(article["AuthorList"][1]["AffiliationInfo"]) == 1
        assert len(article["AuthorList"][1]["AffiliationInfo"][0]) == 2
        assert article["AuthorList"][1]["AffiliationInfo"][0]["Identifier"] == []
        assert article["AuthorList"][1]["AffiliationInfo"][0]["Affiliation"] == "Studies, Research and Sports Medicine Center, Government of Navarre, Pamplona, Spain."
        assert article["AuthorList"][1]["AffiliationInfo"][0].attributes == {}
        assert article["AuthorList"][1]["LastName"] == "Gorostiaga"
        assert article["AuthorList"][1]["ForeName"] == "Esteban M"
        assert article["AuthorList"][1]["Initials"] == "EM"
        assert article["AuthorList"][1].attributes == {"ValidYN": "Y"}
        assert len(article["Abstract"]) == 1
        assert article["Abstract"]["AbstractText"][0] == """\
Maximal Lactate Steady State (MLSS) and Lactate Threshold (LT) are physiologically-related and fundamental concepts within the sports and exercise sciences. Literature supporting their relationship, however, is scarce. Among the recognized LTs, we were particularly interested in the disused "Minimum Lactate Equivalent" (LE<sub>min</sub>), first described in the early 1980s. We hypothesized that velocity at LT, conceptually comprehended as in the old days (LE<sub>min</sub>), could predict velocity at MLSS (<sub>V</sub>MLSS) more accurate than some other blood lactate-related thresholds (BL<sub>R</sub>Ts) routinely used nowadays by many sport science practitioners. Thirteen male endurance-trained [<sub>V</sub>MLSS 15.0 ± 1.1 km·h<sup>-1</sup>; maximal oxygen uptake ( <math xmlns="http://www.w3.org/1998/Math/MathML">
                        <msub>
                            <mrow>
                                <mover>
                                    <mrow>
                                        <mi>V</mi>
                                    </mrow>
                                    <mo>.</mo>
                                </mover>
                                <mi>O</mi>
                            </mrow>
                            <mrow>
                                <mn>2</mn>
                                <mi>m</mi>
                                <mi>a</mi>
                                <mi>x</mi>
                            </mrow>
                        </msub>
                    </math> ) 67.6 ± 4.1 ml·kg<sup>-1</sup>·min<sup>-1</sup>] homogeneous (coefficient of variation: ≈7%) runners conducted 1) a submaximal discontinuous incremental running test to determine several BL<sub>R</sub>Ts followed by a maximal ramp incremental running test for <math xmlns="http://www.w3.org/1998/Math/MathML">
                        <msub>
                            <mrow>
                                <mover>
                                    <mrow>
                                        <mi>V</mi>
                                    </mrow>
                                    <mo>.</mo>
                                </mover>
                                <mi>O</mi>
                            </mrow>
                            <mrow>
                                <mn>2</mn>
                                <mi>m</mi>
                                <mi>a</mi>
                                <mi>x</mi>
                            </mrow>
                        </msub>
                        <mtext> </mtext>
                    </math> determination, and 2) several (4-5) constant velocity running tests to determine <sub>V</sub>MLSS with a precision of 0.20 km·h<sup>-1</sup>. Determined BL<sub>R</sub>Ts include LE<sub>min</sub> and LE<sub>min</sub>-related LE<sub>min</sub> plus 1 (LE<sub>min+1mM</sub>) and 1.5 mmol·L<sup>-1</sup> (LE<sub>min+1.5mM</sub>), along with well-established BL<sub>R</sub>Ts such as conventionally-calculated LT, D<sub>max</sub> and fixed blood lactate concentration thresholds. LE<sub>min</sub> did not differ from LT (<i>P</i> = 0.71; ES: 0.08) and was 27% lower than MLSS (<i>P</i> < 0.001; ES: 3.54). LE<sub>min+1mM</sub> was not different from MLSS (<i>P</i> = 0.47; ES: 0.09). LE<sub>min</sub> was the best predictor of <sub>V</sub>MLSS (<i>r</i> = 0.91; <i>P</i> < 0.001; SEE = 0.47 km·h<sup>-1</sup>), followed by LE<sub>min+1mM</sub> (<i>r</i> = 0.86; <i>P</i> < 0.001; SEE = 0.58 km·h<sup>-1</sup>) and LE<sub>min+1.5mM</sub> (<i>r</i> = 0.84; <i>P</i> < 0.001; SEE = 0.86 km·h<sup>-1</sup>). There was no statistical difference between MLSS and estimated MLSS using LE<sub>min</sub> prediction formula (<i>P</i> = 0.99; ES: 0.001). Mean bias and limits of agreement were 0.00 ± 0.45 km·h<sup>-1</sup> and ±0.89 km·h<sup>-1</sup>. Additionally, LE<sub>min</sub>, LE<sub>min+1mM</sub> and LE<sub>min+1.5mM</sub> were the best predictors of <math xmlns="http://www.w3.org/1998/Math/MathML">
                        <msub>
                            <mrow>
                                <mover>
                                    <mrow>
                                        <mi>V</mi>
                                    </mrow>
                                    <mo>.</mo>
                                </mover>
                                <mi>O</mi>
                            </mrow>
                            <mrow>
                                <mn>2</mn>
                                <mi>m</mi>
                                <mi>a</mi>
                                <mi>x</mi>
                            </mrow>
                        </msub>
                    </math> (<i>r</i> = 0.72-0.79; <i>P</i> < 0.001). These results support LE<sub>min</sub>, an objective submaximal overlooked and underused BL<sub>R</sub>T, to be one of the best single MLSS predictors in endurance trained runners. Our study advocates factors controlling LE<sub>min</sub> to be shared, at least partly, with those controlling MLSS."""

    def test_pubmed_mathml_tags(self):
        """Test parsing XML returned by EFetch, PubMed database, with extensive MathML tags."""
        # In PubMed display PMID 29963580 in xml retrieval mode, containing
        # extensive MathML tags in the abstract text.
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db="pubmed", id="29963580", retmode="xml")
        with open("Entrez/pubmed7.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert len(record) == 2
        assert record["PubmedBookArticle"] == []
        assert len(record["PubmedArticle"]) == 1
        pubmed_article = record["PubmedArticle"][0]
        assert len(pubmed_article) == 2
        assert len(pubmed_article["MedlineCitation"].attributes) == 2
        assert pubmed_article["MedlineCitation"].attributes["Status"] == "PubMed-not-MEDLINE"
        assert pubmed_article["MedlineCitation"].attributes["Owner"] == "NLM"
        assert pubmed_article["MedlineCitation"]["PMID"] == "29963580"
        assert pubmed_article["MedlineCitation"]["PMID"].attributes["Version"] == "1"
        assert pubmed_article["MedlineCitation"]["DateRevised"].attributes == {}
        assert pubmed_article["MedlineCitation"]["DateRevised"]["Year"] == "2018"
        assert pubmed_article["MedlineCitation"]["DateRevised"]["Month"] == "11"
        assert pubmed_article["MedlineCitation"]["DateRevised"]["Day"] == "14"
        assert len(pubmed_article["MedlineCitation"]["Article"].attributes) == 1
        assert pubmed_article["MedlineCitation"]["Article"].attributes["PubModel"] == "Print-Electronic"
        assert pubmed_article["MedlineCitation"]["Article"]["Journal"].attributes == {}
        assert len(
                pubmed_article["MedlineCitation"]["Article"]["Journal"][
                    "ISSN"
                ].attributes
            ) == 1
        assert pubmed_article["MedlineCitation"]["Article"]["Journal"]["ISSN"].attributes[
                "IssnType"
            ] == "Print"
        assert pubmed_article["MedlineCitation"]["Article"]["Journal"]["ISSN"] == "2329-4302"
        assert len(
                pubmed_article["MedlineCitation"]["Article"]["Journal"][
                    "JournalIssue"
                ].attributes
            ) == 1
        assert pubmed_article["MedlineCitation"]["Article"]["Journal"][
                "JournalIssue"
            ].attributes["CitedMedium"] == "Print"
        assert pubmed_article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                "Volume"
            ] == "5"
        assert pubmed_article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                "Issue"
            ] == "2"
        assert pubmed_article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                "PubDate"
            ]["Year"] == "2018"
        assert pubmed_article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                "PubDate"
            ]["Month"] == "Apr"
        assert pubmed_article["MedlineCitation"]["Article"]["Journal"]["Title"] == "Journal of medical imaging (Bellingham, Wash.)"
        assert pubmed_article["MedlineCitation"]["Article"]["Journal"]["ISOAbbreviation"] == "J Med Imaging (Bellingham)"
        assert pubmed_article["MedlineCitation"]["Article"]["ArticleTitle"] == "Development of a pulmonary imaging biomarker pipeline for phenotyping of chronic lung disease."
        assert len(pubmed_article["MedlineCitation"]["Article"]["Pagination"]) == 1
        assert pubmed_article["MedlineCitation"]["Article"]["Pagination"]["MedlinePgn"] == "026002"
        assert len(pubmed_article["MedlineCitation"]["Article"]["ELocationID"]) == 1
        assert pubmed_article["MedlineCitation"]["Article"]["ELocationID"][0].attributes[
                "EIdType"
            ] == "doi"
        assert pubmed_article["MedlineCitation"]["Article"]["ELocationID"][0].attributes[
                "ValidYN"
            ] == "Y"
        assert pubmed_article["MedlineCitation"]["Article"]["ELocationID"][0] == "10.1117/1.JMI.5.2.026002"
        assert len(pubmed_article["MedlineCitation"]["Article"]["Abstract"]) == 1
        assert len(
                pubmed_article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]
            ) == 1
        assert pubmed_article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][0] == """\
We designed and generated pulmonary imaging biomarker pipelines to facilitate high-throughput research and point-of-care use in patients with chronic lung disease. Image processing modules and algorithm pipelines were embedded within a graphical user interface (based on the .NET framework) for pulmonary magnetic resonance imaging (MRI) and x-ray computed-tomography (CT) datasets. The software pipelines were generated using C++ and included: (1) inhaled <math xmlns="http://www.w3.org/1998/Math/MathML">
                        <mrow>
                            <mmultiscripts>
                                <mrow>
                                    <mi>He</mi>
                                </mrow>
                                <mprescripts></mprescripts>
                                <none></none>
                                <mrow>
                                    <mn>3</mn>
                                </mrow>
                            </mmultiscripts>
                            <mo>/</mo>
                            <mmultiscripts>
                                <mrow>
                                    <mi>Xe</mi>
                                </mrow>
                                <mprescripts></mprescripts>
                                <none></none>
                                <mrow>
                                    <mn>129</mn>
                                </mrow>
                            </mmultiscripts>
                            <mtext> </mtext>
                            <mi>MRI</mi>
                        </mrow>
                    </math> ventilation and apparent diffusion coefficients, (2) CT-MRI coregistration for lobar and segmental ventilation and perfusion measurements, (3) ultrashort echo-time <math xmlns="http://www.w3.org/1998/Math/MathML">
                        <mrow>
                            <mmultiscripts>
                                <mrow>
                                    <mi>H</mi>
                                </mrow>
                                <mprescripts></mprescripts>
                                <none></none>
                                <mrow>
                                    <mn>1</mn>
                                </mrow>
                            </mmultiscripts>
                            <mtext> </mtext>
                            <mi>MRI</mi>
                        </mrow>
                    </math> proton density measurements, (4) free-breathing Fourier-decomposition <math xmlns="http://www.w3.org/1998/Math/MathML">
                        <mrow>
                            <mmultiscripts>
                                <mrow>
                                    <mi>H</mi>
                                </mrow>
                                <mprescripts></mprescripts>
                                <none></none>
                                <mrow>
                                    <mn>1</mn>
                                </mrow>
                            </mmultiscripts>
                            <mtext> </mtext>
                            <mi>MRI</mi>
                        </mrow>
                    </math> ventilation/perfusion and free-breathing <math xmlns="http://www.w3.org/1998/Math/MathML">
                        <mrow>
                            <mmultiscripts>
                                <mrow>
                                    <mi>H</mi>
                                </mrow>
                                <mprescripts></mprescripts>
                                <none></none>
                                <mrow>
                                    <mn>1</mn>
                                </mrow>
                            </mmultiscripts>
                            <mtext> </mtext>
                            <mi>MRI</mi>
                        </mrow>
                    </math> specific ventilation, (5)\u00a0multivolume CT and MRI parametric response maps, and (6)\u00a0MRI and CT texture analysis and radiomics. The image analysis framework was implemented on a desktop workstation/tablet to generate biomarkers of regional lung structure and function related to ventilation, perfusion, lung tissue texture, and integrity as well as multiparametric measures of gas trapping and airspace enlargement. All biomarkers were generated within 10 min with measurement reproducibility consistent with clinical and research requirements. The resultant pulmonary imaging biomarker pipeline provides real-time and automated lung imaging measurements for point-of-care and high-throughput research."""
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"].attributes[
                "CompleteYN"
            ] == "Y"
        assert len(pubmed_article["MedlineCitation"]["Article"]["AuthorList"]) == 9
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][0].attributes[
                "ValidYN"
            ] == "Y"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][0]["LastName"] == "Guo"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][0]["ForeName"] == "Fumin"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][0]["Initials"] == "F"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][0][
                "AffiliationInfo"
            ][0]["Affiliation"] == "University of Western Ontario, Robarts Research Institute, London, Ontario, Canada."
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][0][
                "AffiliationInfo"
            ][1]["Affiliation"] == "University of Western Ontario, Graduate Program in Biomedical Engineering, London, Ontario, Canada."
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][0][
                "AffiliationInfo"
            ][2]["Affiliation"] == "University of Toronto, Sunnybrook Research Institute, Toronto, Canada."
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][1].attributes[
                "ValidYN"
            ] == "Y"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][1]["LastName"] == "Capaldi"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][1]["ForeName"] == "Dante"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][1]["Initials"] == "D"
        assert len(
                pubmed_article["MedlineCitation"]["Article"]["AuthorList"][1][
                    "Identifier"
                ]
            ) == 1
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][1]["Identifier"][
                0
            ].attributes["Source"] == "ORCID"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][1]["Identifier"][
                0
            ] == "https://orcid.org/0000-0002-4590-7461"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][1][
                "AffiliationInfo"
            ][0]["Affiliation"] == "University of Western Ontario, Robarts Research Institute, London, Ontario, Canada."
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][1][
                "AffiliationInfo"
            ][1]["Affiliation"] == "University of Western Ontario, Department of Medical Biophysics, London, Ontario, Canada."
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][2].attributes[
                "ValidYN"
            ] == "Y"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][2]["LastName"] == "Kirby"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][2]["ForeName"] == "Miranda"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][2]["Initials"] == "M"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][2][
                "AffiliationInfo"
            ][0]["Affiliation"] == "University of British Columbia, St. Paul's Hospital, Centre for Heart Lung Innovation, Vancouver, Canada."
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][3].attributes[
                "ValidYN"
            ] == "Y"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][3]["LastName"] == "Sheikh"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][3]["ForeName"] == "Khadija"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][3]["Initials"] == "K"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][3][
                "AffiliationInfo"
            ][0]["Affiliation"] == "University of Western Ontario, Robarts Research Institute, London, Ontario, Canada."
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][4].attributes[
                "ValidYN"
            ] == "Y"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][4]["LastName"] == "Svenningsen"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][4]["ForeName"] == "Sarah"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][4]["Initials"] == "S"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][4][
                "AffiliationInfo"
            ][0]["Affiliation"] == "University of Western Ontario, Robarts Research Institute, London, Ontario, Canada."
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][5].attributes[
                "ValidYN"
            ] == "Y"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][5]["LastName"] == "McCormack"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][5]["ForeName"] == "David G"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][5]["Initials"] == "DG"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][5][
                "AffiliationInfo"
            ][0]["Affiliation"] == "University of Western Ontario, Division of Respirology, Department of Medicine, London, Ontario, Canada."
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][6].attributes[
                "ValidYN"
            ] == "Y"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][6]["LastName"] == "Fenster"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][6]["ForeName"] == "Aaron"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][6]["Initials"] == "A"
        assert len(
                pubmed_article["MedlineCitation"]["Article"]["AuthorList"][6][
                    "Identifier"
                ]
            ) == 1
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][6]["Identifier"][
                0
            ].attributes["Source"] == "ORCID"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][6]["Identifier"][
                0
            ] == "https://orcid.org/0000-0003-3525-2788"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][6][
                "AffiliationInfo"
            ][0]["Affiliation"] == "University of Western Ontario, Robarts Research Institute, London, Ontario, Canada."
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][6][
                "AffiliationInfo"
            ][1]["Affiliation"] == "University of Western Ontario, Graduate Program in Biomedical Engineering, London, Ontario, Canada."
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][6][
                "AffiliationInfo"
            ][2]["Affiliation"] == "University of Western Ontario, Department of Medical Biophysics, London, Ontario, Canada."
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][7].attributes[
                "ValidYN"
            ] == "Y"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][7]["LastName"] == "Parraga"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][7]["ForeName"] == "Grace"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][7]["Initials"] == "G"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][7][
                "AffiliationInfo"
            ][0]["Affiliation"] == "University of Western Ontario, Robarts Research Institute, London, Ontario, Canada."
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][7][
                "AffiliationInfo"
            ][1]["Affiliation"] == "University of Western Ontario, Graduate Program in Biomedical Engineering, London, Ontario, Canada."
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][7][
                "AffiliationInfo"
            ][2]["Affiliation"] == "University of Western Ontario, Department of Medical Biophysics, London, Ontario, Canada."
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][8].attributes[
                "ValidYN"
            ] == "Y"
        assert pubmed_article["MedlineCitation"]["Article"]["AuthorList"][8][
                "CollectiveName"
            ] == "Canadian Respiratory Research Network"
        assert len(pubmed_article["MedlineCitation"]["Article"]["Language"]) == 1
        assert pubmed_article["MedlineCitation"]["Article"]["Language"][0] == "eng"
        assert len(pubmed_article["MedlineCitation"]["Article"]["PublicationTypeList"]) == 1
        assert pubmed_article["MedlineCitation"]["Article"]["PublicationTypeList"][
                0
            ].attributes["UI"] == "D016428"
        assert pubmed_article["MedlineCitation"]["Article"]["PublicationTypeList"][0] == "Journal Article"
        assert len(pubmed_article["MedlineCitation"]["Article"]["ArticleDate"]) == 1
        assert pubmed_article["MedlineCitation"]["Article"]["ArticleDate"][0].attributes[
                "DateType"
            ] == "Electronic"
        assert pubmed_article["MedlineCitation"]["Article"]["ArticleDate"][0]["Year"] == "2018"
        assert pubmed_article["MedlineCitation"]["Article"]["ArticleDate"][0]["Month"] == "06"
        assert pubmed_article["MedlineCitation"]["Article"]["ArticleDate"][0]["Day"] == "28"
        assert len(pubmed_article["MedlineCitation"]["MedlineJournalInfo"]) == 4
        assert pubmed_article["MedlineCitation"]["MedlineJournalInfo"]["Country"] == "United States"
        assert pubmed_article["MedlineCitation"]["MedlineJournalInfo"]["MedlineTA"] == "J Med Imaging (Bellingham)"
        assert pubmed_article["MedlineCitation"]["MedlineJournalInfo"]["NlmUniqueID"] == "101643461"
        assert pubmed_article["MedlineCitation"]["MedlineJournalInfo"]["ISSNLinking"] == "2329-4302"
        assert len(pubmed_article["MedlineCitation"]["KeywordList"]) == 1
        assert pubmed_article["MedlineCitation"]["KeywordList"][0].attributes["Owner"] == "NOTNLM"
        assert len(pubmed_article["MedlineCitation"]["KeywordList"][0]) == 5
        assert pubmed_article["MedlineCitation"]["KeywordList"][0][0].attributes[
                "MajorTopicYN"
            ] == "N"
        assert pubmed_article["MedlineCitation"]["KeywordList"][0][0] == "asthma"
        assert pubmed_article["MedlineCitation"]["KeywordList"][0][1].attributes[
                "MajorTopicYN"
            ] == "N"
        assert pubmed_article["MedlineCitation"]["KeywordList"][0][1] == "chronic obstructive lung disease"
        assert pubmed_article["MedlineCitation"]["KeywordList"][0][2].attributes[
                "MajorTopicYN"
            ] == "N"
        assert pubmed_article["MedlineCitation"]["KeywordList"][0][2] == "image processing, biomarkers"
        assert pubmed_article["MedlineCitation"]["KeywordList"][0][3].attributes[
                "MajorTopicYN"
            ] == "N"
        assert pubmed_article["MedlineCitation"]["KeywordList"][0][3] == "magnetic resonance imaging"
        assert pubmed_article["MedlineCitation"]["KeywordList"][0][4].attributes[
                "MajorTopicYN"
            ] == "N"
        assert pubmed_article["MedlineCitation"]["KeywordList"][0][4] == "thoracic computed tomography"
        assert pubmed_article["PubmedData"]["History"][0].attributes["PubStatus"] == "received"
        assert pubmed_article["PubmedData"]["History"][0]["Year"] == "2017"
        assert pubmed_article["PubmedData"]["History"][0]["Month"] == "12"
        assert pubmed_article["PubmedData"]["History"][0]["Day"] == "12"
        assert pubmed_article["PubmedData"]["History"][1].attributes["PubStatus"] == "accepted"
        assert pubmed_article["PubmedData"]["History"][1]["Year"] == "2018"
        assert pubmed_article["PubmedData"]["History"][1]["Month"] == "06"
        assert pubmed_article["PubmedData"]["History"][1]["Day"] == "14"
        assert pubmed_article["PubmedData"]["History"][2].attributes["PubStatus"] == "pmc-release"
        assert pubmed_article["PubmedData"]["History"][2]["Year"] == "2019"
        assert pubmed_article["PubmedData"]["History"][2]["Month"] == "06"
        assert pubmed_article["PubmedData"]["History"][2]["Day"] == "28"
        assert pubmed_article["PubmedData"]["History"][3].attributes["PubStatus"] == "entrez"
        assert pubmed_article["PubmedData"]["History"][3]["Year"] == "2018"
        assert pubmed_article["PubmedData"]["History"][3]["Month"] == "7"
        assert pubmed_article["PubmedData"]["History"][3]["Day"] == "3"
        assert pubmed_article["PubmedData"]["History"][3]["Hour"] == "6"
        assert pubmed_article["PubmedData"]["History"][3]["Minute"] == "0"
        assert pubmed_article["PubmedData"]["History"][4].attributes["PubStatus"] == "pubmed"
        assert pubmed_article["PubmedData"]["History"][4]["Year"] == "2018"
        assert pubmed_article["PubmedData"]["History"][4]["Month"] == "7"
        assert pubmed_article["PubmedData"]["History"][4]["Day"] == "3"
        assert pubmed_article["PubmedData"]["History"][4]["Hour"] == "6"
        assert pubmed_article["PubmedData"]["History"][4]["Minute"] == "0"
        assert pubmed_article["PubmedData"]["History"][5].attributes["PubStatus"] == "medline"
        assert pubmed_article["PubmedData"]["History"][5]["Year"] == "2018"
        assert pubmed_article["PubmedData"]["History"][5]["Month"] == "7"
        assert pubmed_article["PubmedData"]["History"][5]["Day"] == "3"
        assert pubmed_article["PubmedData"]["History"][5]["Hour"] == "6"
        assert pubmed_article["PubmedData"]["History"][5]["Minute"] == "1"
        assert pubmed_article["PubmedData"]["PublicationStatus"] == "ppublish"
        assert len(pubmed_article["PubmedData"]["ArticleIdList"]) == 4
        assert pubmed_article["PubmedData"]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert pubmed_article["PubmedData"]["ArticleIdList"][0] == "29963580"
        assert pubmed_article["PubmedData"]["ArticleIdList"][1].attributes["IdType"] == "doi"
        assert pubmed_article["PubmedData"]["ArticleIdList"][1] == "10.1117/1.JMI.5.2.026002"
        assert pubmed_article["PubmedData"]["ArticleIdList"][2].attributes["IdType"] == "pii"
        assert pubmed_article["PubmedData"]["ArticleIdList"][2] == "17360RR"
        assert pubmed_article["PubmedData"]["ArticleIdList"][3].attributes["IdType"] == "pmc"
        assert pubmed_article["PubmedData"]["ArticleIdList"][3] == "PMC6022861"
        assert len(pubmed_article["PubmedData"]["ReferenceList"]) == 1
        assert len(pubmed_article["PubmedData"]["ReferenceList"][0]) == 2
        assert len(pubmed_article["PubmedData"]["ReferenceList"][0]["ReferenceList"]) == 0
        references = pubmed_article["PubmedData"]["ReferenceList"][0]["Reference"]
        assert len(references) == 49
        assert references[0]["Citation"] == "Radiology. 2015 Jan;274(1):250-9"
        assert len(references[0]["ArticleIdList"]) == 1
        assert references[0]["ArticleIdList"][0] == "25144646"
        assert references[0]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[1]["Citation"] == "Nature. 1994 Jul 21;370(6486):199-201"
        assert len(references[1]["ArticleIdList"]) == 1
        assert references[1]["ArticleIdList"][0] == "8028666"
        assert references[1]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[2]["Citation"] == "Magn Reson Med. 2009 Sep;62(3):656-64"
        assert len(references[2]["ArticleIdList"]) == 1
        assert references[2]["ArticleIdList"][0] == "19585597"
        assert references[2]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[3]["Citation"] == "Radiology. 2016 Feb;278(2):563-77"
        assert len(references[3]["ArticleIdList"]) == 1
        assert references[3]["ArticleIdList"][0] == "26579733"
        assert references[3]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[4]["Citation"] == "Radiology. 1991 Jun;179(3):777-81"
        assert len(references[4]["ArticleIdList"]) == 1
        assert references[4]["ArticleIdList"][0] == "2027991"
        assert references[4]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[5]["Citation"] == "Radiology. 2010 Jul;256(1):280-9"
        assert len(references[5]["ArticleIdList"]) == 1
        assert references[5]["ArticleIdList"][0] == "20574101"
        assert references[5]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[6]["Citation"] == "Phys Med Biol. 2001 May;46(5):R67-99"
        assert len(references[6]["ArticleIdList"]) == 1
        assert references[6]["ArticleIdList"][0] == "11384074"
        assert references[6]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[7]["Citation"] == "IEEE Trans Med Imaging. 2011 Nov;30(11):1901-20"
        assert len(references[7]["ArticleIdList"]) == 1
        assert references[7]["ArticleIdList"][0] == "21632295"
        assert references[7]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[8]["Citation"] == "Eur J Cancer. 2012 Mar;48(4):441-6"
        assert len(references[8]["ArticleIdList"]) == 1
        assert references[8]["ArticleIdList"][0] == "22257792"
        assert references[8]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[9]["Citation"] == "Med Phys. 2017 May;44(5):1718-1733"
        assert len(references[9]["ArticleIdList"]) == 1
        assert references[9]["ArticleIdList"][0] == "28206676"
        assert references[9]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[10]["Citation"] == "COPD. 2014 Apr;11(2):125-32"
        assert len(references[10]["ArticleIdList"]) == 1
        assert references[10]["ArticleIdList"][0] == "22433011"
        assert references[10]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[11]["Citation"] == "Am J Respir Crit Care Med. 2015 Nov 15;192(10):1215-22"
        assert len(references[11]["ArticleIdList"]) == 1
        assert references[11]["ArticleIdList"][0] == "26186608"
        assert references[11]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[12]["Citation"] == "N Engl J Med. 2016 May 12;374(19):1811-21"
        assert len(references[12]["ArticleIdList"]) == 1
        assert references[12]["ArticleIdList"][0] == "27168432"
        assert references[12]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[13]["Citation"] == "Am J Epidemiol. 2002 Nov 1;156(9):871-81"
        assert len(references[13]["ArticleIdList"]) == 1
        assert references[13]["ArticleIdList"][0] == "12397006"
        assert references[13]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[14]["Citation"] == "BMC Cancer. 2014 Dec 11;14:934"
        assert len(references[14]["ArticleIdList"]) == 1
        assert references[14]["ArticleIdList"][0] == "25496482"
        assert references[14]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[15]["Citation"] == "Acad Radiol. 2015 Mar;22(3):320-9"
        assert len(references[15]["ArticleIdList"]) == 1
        assert references[15]["ArticleIdList"][0] == "25491735"
        assert references[15]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[16]["Citation"] == "Chest. 1999 Dec;116(6):1750-61"
        assert len(references[16]["ArticleIdList"]) == 1
        assert references[16]["ArticleIdList"][0] == "10593802"
        assert references[16]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[17]["Citation"] == "Acad Radiol. 2012 Feb;19(2):141-52"
        assert len(references[17]["ArticleIdList"]) == 1
        assert references[17]["ArticleIdList"][0] == "22104288"
        assert references[17]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[18]["Citation"] == "Med Phys. 2014 Mar;41(3):033502"
        assert len(references[18]["ArticleIdList"]) == 1
        assert references[18]["ArticleIdList"][0] == "24593744"
        assert references[18]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[19]["Citation"] == "Med Phys. 2008 Oct;35(10):4695-707"
        assert len(references[19]["ArticleIdList"]) == 1
        assert references[19]["ArticleIdList"][0] == "18975715"
        assert references[19]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[20]["Citation"] == "Thorax. 2017 May;72(5):475-477"
        assert len(references[20]["ArticleIdList"]) == 1
        assert references[20]["ArticleIdList"][0] == "28258250"
        assert references[20]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[21]["Citation"] == "Nat Med. 1996 Nov;2(11):1236-9"
        assert len(references[21]["ArticleIdList"]) == 1
        assert references[21]["ArticleIdList"][0] == "8898751"
        assert references[21]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[22]["Citation"] == "J Magn Reson Imaging. 2015 May;41(5):1465-74"
        assert len(references[22]["ArticleIdList"]) == 1
        assert references[22]["ArticleIdList"][0] == "24965907"
        assert references[22]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[23]["Citation"] == "Acad Radiol. 2008 Jun;15(6):776-85"
        assert len(references[23]["ArticleIdList"]) == 1
        assert references[23]["ArticleIdList"][0] == "18486013"
        assert references[23]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[24]["Citation"] == "Magn Reson Med. 2000 Aug;44(2):174-9"
        assert len(references[24]["ArticleIdList"]) == 1
        assert references[24]["ArticleIdList"][0] == "10918314"
        assert references[24]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[25]["Citation"] == "Am J Respir Crit Care Med. 2014 Jul 15;190(2):135-44"
        assert len(references[25]["ArticleIdList"]) == 1
        assert references[25]["ArticleIdList"][0] == "24873985"
        assert references[25]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[26]["Citation"] == "Med Phys. 2016 Jun;43(6):2911-2926"
        assert len(references[26]["ArticleIdList"]) == 1
        assert references[26]["ArticleIdList"][0] == "27277040"
        assert references[26]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[27]["Citation"] == "Nat Med. 2009 May;15(5):572-6"
        assert len(references[27]["ArticleIdList"]) == 1
        assert references[27]["ArticleIdList"][0] == "19377487"
        assert references[27]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[28]["Citation"] == "Eur J Radiol. 2014 Nov;83(11):2093-101"
        assert len(references[28]["ArticleIdList"]) == 1
        assert references[28]["ArticleIdList"][0] == "25176287"
        assert references[28]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[29]["Citation"] == "Radiology. 2004 Sep;232(3):739-48"
        assert len(references[29]["ArticleIdList"]) == 1
        assert references[29]["ArticleIdList"][0] == "15333795"
        assert references[29]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[30]["Citation"] == "Med Image Anal. 2015 Jul;23(1):43-55"
        assert len(references[30]["ArticleIdList"]) == 1
        assert references[30]["ArticleIdList"][0] == "25958028"
        assert references[30]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[31]["Citation"] == "Radiology. 2015 Oct;277(1):192-205"
        assert len(references[31]["ArticleIdList"]) == 1
        assert references[31]["ArticleIdList"][0] == "25961632"
        assert references[31]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[32]["Citation"] == "Med Image Anal. 2012 Oct;16(7):1423-35"
        assert len(references[32]["ArticleIdList"]) == 1
        assert references[32]["ArticleIdList"][0] == "22722056"
        assert references[32]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[33]["Citation"] == "Radiology. 2016 May;279(2):597-608"
        assert len(references[33]["ArticleIdList"]) == 1
        assert references[33]["ArticleIdList"][0] == "26744928"
        assert references[33]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[34]["Citation"] == "J Allergy Clin Immunol. 2003 Jun;111(6):1205-11"
        assert len(references[34]["ArticleIdList"]) == 1
        assert references[34]["ArticleIdList"][0] == "12789218"
        assert references[34]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[35]["Citation"] == "J Magn Reson Imaging. 2016 Mar;43(3):544-57"
        assert len(references[35]["ArticleIdList"]) == 1
        assert references[35]["ArticleIdList"][0] == "26199216"
        assert references[35]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[36]["Citation"] == "Am J Respir Crit Care Med. 2016 Oct 1;194(7):794-806"
        assert len(references[36]["ArticleIdList"]) == 1
        assert references[36]["ArticleIdList"][0] == "27482984"
        assert references[36]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[37]["Citation"] == "Radiology. 1996 Nov;201(2):564-8"
        assert len(references[37]["ArticleIdList"]) == 1
        assert references[37]["ArticleIdList"][0] == "8888259"
        assert references[37]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[38]["Citation"] == "Thorax. 2014 May;69(5):491-4"
        assert len(references[38]["ArticleIdList"]) == 1
        assert references[38]["ArticleIdList"][0] == "24029743"
        assert references[38]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[39]["Citation"] == "J Magn Reson Imaging. 2017 Apr;45(4):1204-1215"
        assert len(references[39]["ArticleIdList"]) == 1
        assert references[39]["ArticleIdList"][0] == "27731948"
        assert references[39]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[40]["Citation"] == "J Appl Physiol (1985). 2009 Oct;107(4):1258-65"
        assert len(references[40]["ArticleIdList"]) == 1
        assert references[40]["ArticleIdList"][0] == "19661452"
        assert references[40]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[41]["Citation"] == "Acad Radiol. 2016 Feb;23(2):176-85"
        assert len(references[41]["ArticleIdList"]) == 1
        assert references[41]["ArticleIdList"][0] == "26601971"
        assert references[41]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[42]["Citation"] == "Radiology. 2018 May;287(2):693-704"
        assert len(references[42]["ArticleIdList"]) == 1
        assert references[42]["ArticleIdList"][0] == "29470939"
        assert references[42]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[43]["Citation"] == "Eur Respir J. 2016 Aug;48(2):370-9"
        assert len(references[43]["ArticleIdList"]) == 1
        assert references[43]["ArticleIdList"][0] == "27174885"
        assert references[43]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[44]["Citation"] == "Radiology. 2011 Oct;261(1):283-92"
        assert len(references[44]["ArticleIdList"]) == 1
        assert references[44]["ArticleIdList"][0] == "21813741"
        assert references[44]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[45]["Citation"] == "Am J Respir Crit Care Med. 2014 Mar 15;189(6):650-7"
        assert len(references[45]["ArticleIdList"]) == 1
        assert references[45]["ArticleIdList"][0] == "24401150"
        assert references[45]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[46]["Citation"] == "Am J Respir Crit Care Med. 2012 Feb 15;185(4):356-62"
        assert len(references[46]["ArticleIdList"]) == 1
        assert references[46]["ArticleIdList"][0] == "22095547"
        assert references[46]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[47]["Citation"] == "COPD. 2010 Feb;7(1):32-43"
        assert len(references[47]["ArticleIdList"]) == 1
        assert references[47]["ArticleIdList"][0] == "20214461"
        assert references[47]["ArticleIdList"][0].attributes["IdType"] == "pubmed"
        assert references[48]["Citation"] == "Eur Respir J. 2008 Apr;31(4):869-73"
        assert len(references[48]["ArticleIdList"]) == 1
        assert references[48]["ArticleIdList"][0] == "18216052"
        assert references[48]["ArticleIdList"][0].attributes["IdType"] == "pubmed"

    def test_pmc(self):
        """Test parsing XML returned by EFetch from PubMed Central."""
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db='pmc', id="2682512,3468381")
        with open("Entrez/efetch_pmc.xml", "rb") as stream:
            records = Entrez.parse(stream)
            records = list(records)
        assert len(records) == 1
        record = records[0]
        assert len(record) == 2
        assert len(record["front"]) == 9
        assert len(record["front"]["journal-meta"]) == 10
        assert len(record["front"]["journal-meta"]["journal-id"]) == 4
        assert record["front"]["journal-meta"]["journal-id"][0] == "ERJ Open Res"
        assert record["front"]["journal-meta"]["journal-id"][0].attributes == {"journal-id-type": "nlm-ta"}
        assert record["front"]["journal-meta"]["journal-id"][1] == "ERJ Open Res"
        assert record["front"]["journal-meta"]["journal-id"][1].attributes == {"journal-id-type": "iso-abbrev"}
        assert record["front"]["journal-meta"]["journal-id"][2] == "ERJOR"
        assert record["front"]["journal-meta"]["journal-id"][2].attributes == {"journal-id-type": "publisher-id"}
        assert record["front"]["journal-meta"]["journal-id"][3] == "erjor"
        assert record["front"]["journal-meta"]["journal-id"][3].attributes == {"journal-id-type": "hwp"}
        assert len(record["front"]["journal-meta"]["journal-title-group"]) == 1
        journal_title_group = record["front"]["journal-meta"]["journal-title-group"][0]
        assert len(journal_title_group) == 4
        assert journal_title_group["journal-title"] == ["ERJ Open Research"]
        assert journal_title_group["journal-subtitle"] == []
        assert journal_title_group["abbrev-journal-title"] == []
        assert journal_title_group["trans-title-group"] == []
        assert len(record["front"]["journal-meta"]["issn"]) == 1
        assert record["front"]["journal-meta"]["issn"][0] == "2312-0541"
        assert record["front"]["journal-meta"]["issn"][0].attributes == {"pub-type": "epub"}
        assert len(record["front"]["journal-meta"]["publisher"]) == 1
        assert len(record["front"]["journal-meta"]["publisher"][0]) == 1
        assert record["front"]["journal-meta"]["publisher"][0][0] == "European Respiratory Society"
        assert record["front"]["journal-meta"]["publisher"][0][0].tag == "publisher-name"
        assert record["front"]["journal-meta"]["contrib-group"] == []
        assert record["front"]["journal-meta"]["notes"] == []
        assert record["front"]["journal-meta"]["aff"] == []
        assert record["front"]["journal-meta"]["aff-alternatives"] == []
        assert record["front"]["journal-meta"]["self-uri"] == []
        assert record["front"]["journal-meta"]["isbn"] == []
        assert len(record["front"]["article-meta"]) == 34
        assert record["front"]["article-meta"]["abstract"] == []
        assert record["front"]["article-meta"]["funding-group"] == []
        assert record["front"]["article-meta"]["aff"] == []
        assert record["front"]["article-meta"]["issue-title"] == []
        assert len(record["front"]["article-meta"]["pub-date"]) == 3
        assert record["front"]["article-meta"]["pub-date"][0] == ["7", "2021"]
        assert record["front"]["article-meta"]["pub-date"][0].attributes == {"pub-type": "collection"}
        assert record["front"]["article-meta"]["pub-date"][1] == ["13", "9", "2021"]
        assert record["front"]["article-meta"]["pub-date"][1].attributes == {"pub-type": "epub"}
        assert record["front"]["article-meta"]["pub-date"][2] == ["13", "9", "2021"]
        assert record["front"]["article-meta"]["pub-date"][2].attributes == {"pub-type": "pmc-release"}
        assert record["front"]["article-meta"]["conference"] == []
        assert record["front"]["article-meta"]["supplementary-material"] == []
        assert len(record["front"]["article-meta"]["related-article"]) == 1
        assert record["front"]["article-meta"]["related-article"][0] == ""
        assert record["front"]["article-meta"]["related-article"][0].attributes == {
                "related-article-type": "corrected-article",
                "id": "d31e52",
                "ext-link-type": "doi",
                "http://www.w3.org/1999/xlink href": "10.1183/23120541.00193-2021",
            }
        assert record["front"]["article-meta"]["kwd-group"] == []
        assert record["front"]["article-meta"]["contrib-group"] == []
        assert record["front"]["article-meta"]["issue-sponsor"] == []
        assert record["front"]["article-meta"]["self-uri"] == []
        assert record["front"]["article-meta"]["product"] == []
        assert record["front"]["article-meta"]["issue"] == ["3"]
        assert record["front"]["article-meta"]["ext-link"] == []
        assert record["front"]["article-meta"]["support-group"] == []
        assert len(record["front"]["article-meta"]["article-id"]) == 4
        assert record["front"]["article-meta"]["article-id"][0] == "34527728"
        assert record["front"]["article-meta"]["article-id"][0].attributes == {"pub-id-type": "pmid"}
        assert record["front"]["article-meta"]["article-id"][1] == "8435807"
        assert record["front"]["article-meta"]["article-id"][1].attributes == {"pub-id-type": "pmc"}
        assert record["front"]["article-meta"]["article-id"][2] == "10.1183/23120541.50193-2021"
        assert record["front"]["article-meta"]["article-id"][2].attributes == {"pub-id-type": "doi"}
        assert record["front"]["article-meta"]["article-id"][3] == "50193-2021"
        assert record["front"]["article-meta"]["article-id"][3].attributes == {"pub-id-type": "publisher-id"}
        assert record["front"]["article-meta"]["issue-title-group"] == []
        assert record["front"]["article-meta"]["x"] == []
        assert record["front"]["article-meta"]["uri"] == []
        assert record["front"]["article-meta"]["email"] == []
        assert record["front"]["article-meta"]["volume-id"] == []
        assert record["front"]["article-meta"]["issue-id"] == []
        assert record["front"]["article-meta"]["trans-abstract"] == []
        assert record["front"]["article-meta"]["volume-issue-group"] == []
        assert record["front"]["article-meta"]["related-object"] == []
        assert record["front"]["article-meta"]["isbn"] == []
        assert record["front"]["article-meta"]["volume"] == ["7"]
        assert record["front"]["article-meta"]["aff-alternatives"] == []
        assert record["front"]["article-meta"]["article-version"] == "Version of Record"
        assert len(record["front"]["article-meta"]["article-version"].attributes) == 3
        assert record["front"]["article-meta"]["article-version"].attributes["vocab"] == "JAV"
        assert record["front"]["article-meta"]["article-version"].attributes[
                "vocab-identifier"
            ] == "http://www.niso.org/publications/rp/RP-8-2008.pdf"
        assert record["front"]["article-meta"]["article-version"].attributes[
                "article-version-type"
            ] == "VoR"
        assert len(record["front"]["article-meta"]["article-categories"]) == 3
        assert record["front"]["article-meta"]["article-categories"]["series-text"] == []
        assert len(record["front"]["article-meta"]["article-categories"]["subj-group"]) == 1
        assert len(record["front"]["article-meta"]["article-categories"]["subj-group"][0]) == 3
        assert record["front"]["article-meta"]["article-categories"]["subj-group"][0][
                "subject"
            ] == ["Author Correction"]
        assert record["front"]["article-meta"]["article-categories"]["subj-group"][0][
                "subj-group"
            ] == []
        assert record["front"]["article-meta"]["article-categories"]["subj-group"][0][
                "compound-subject"
            ] == []
        assert record["front"]["article-meta"]["article-categories"]["subj-group"][
                0
            ].attributes == {"subj-group-type": "heading"}
        assert len(record["front"]["article-meta"]["title-group"]) == 4
        assert record["front"]["article-meta"]["title-group"]["trans-title-group"] == []
        assert record["front"]["article-meta"]["title-group"]["alt-title"] == []
        assert record["front"]["article-meta"]["title-group"]["subtitle"] == []
        assert record["front"]["article-meta"]["title-group"]["article-title"] == '“Lung diffusing capacity for nitric oxide measured by two commercial devices: a randomised crossover comparison in healthy adults”. Thomas Radtke, Quintin de Groot, Sarah R. Haile, Marion Maggi, Connie C.W. Hsia and Holger Dressel. <italic toggle="yes">ERJ Open Res</italic> 2021; 7: 00193-2021.'
        assert record["front"]["article-meta"]["elocation-id"] == "50193-2021"
        assert len(record["front"]["article-meta"]["permissions"]) == 5
        assert record["front"]["article-meta"]["permissions"]["copyright-year"] == ["2021"]
        assert record["front"]["article-meta"]["permissions"]["copyright-holder"] == []
        assert len(record["front"]["article-meta"]["permissions"]["license"]) == 1
        assert len(record["front"]["article-meta"]["permissions"]["license"][0]) == 2
        assert record["front"]["article-meta"]["permissions"]["license"][0][0] == "https://creativecommons.org/licenses/by-nc/4.0/"
        assert record["front"]["article-meta"]["permissions"]["license"][0][0].attributes == {"specific-use": "textmining", "content-type": "ccbynclicense"}
        assert record["front"]["article-meta"]["permissions"]["license"][0][1] == 'This version is distributed under the terms of the Creative Commons Attribution Non-Commercial Licence 4.0. For commercial reproduction rights and permissions contact <ext-link ext-link-type="uri" http://www.w3.org/1999/xlink href="mailto:permissions@ersnet.org">permissions@ersnet.org</ext-link>'
        assert record["front"]["article-meta"]["permissions"]["copyright-statement"] == ["Copyright ©The authors 2021"]
        assert record["front"]["article-meta"]["permissions"]["ali:free_to_read"] == []
        assert record["front"]["glossary"] == []
        assert record["front"]["fn-group"] == []
        assert record["front"]["notes"] == []
        assert record["front"]["bio"] == []
        assert record["front"]["list"] == []
        assert record["front"]["def-list"] == []
        assert record["front"]["ack"] == []
        assert len(record["body"]) == 37

        assert record["body"]["table-wrap-group"] == []
        assert record["body"]["disp-formula"] == []
        assert record["body"]["answer-set"] == []
        assert record["body"]["graphic"] == []
        assert record["body"]["statement"] == []
        assert record["body"]["fig-group"] == []
        assert record["body"]["verse-group"] == []
        assert record["body"]["supplementary-material"] == []
        assert record["body"]["related-article"] == []
        assert record["body"]["code"] == []
        assert record["body"]["question"] == []
        assert record["body"]["preformat"] == []
        assert record["body"]["tex-math"] == []
        assert record["body"]["mml:math"] == []
        assert record["body"]["speech"] == []
        assert record["body"]["block-alternatives"] == []
        assert record["body"]["explanation"] == []
        assert record["body"]["array"] == []
        assert record["body"]["question-wrap-group"] == []
        assert record["body"]["alternatives"] == []
        assert record["body"]["media"] == []
        assert record["body"]["x"] == []
        assert record["body"]["sec"] == []
        assert record["body"]["address"] == []
        assert record["body"]["disp-quote"] == []
        assert record["body"]["table-wrap"] == []
        assert record["body"]["ack"] == []
        assert record["body"]["chem-struct-wrap"] == []
        assert record["body"]["related-object"] == []
        assert record["body"]["list"] == []
        assert record["body"]["def-list"] == []
        assert record["body"]["p"] == [
                "This article was originally published with an error in table 2. The upper 95% confidence limit of the per cent difference in the primary end-point (diffusing capacity of the lung for nitric oxide) was incorrectly given as 15.1% and has now been corrected to −15.1% in the published article.\n"
            ]
        assert record["body"]["fig"] == []
        assert record["body"]["answer"] == []
        assert record["body"]["boxed-text"] == []
        assert record["body"]["disp-formula-group"] == []
        assert record["body"]["question-wrap"] == []

    def test_taxonomy(self):
        # Access the Taxonomy database using efetch.
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db="taxonomy", id="9615,9685", retmode="xml")
        with open("Entrez/taxonomy.xml", "rb") as stream:
            record = Entrez.read(stream)
        # fmt: off
        assert len(record) == 2
        assert record[0]["TaxId"] == "9615"
        assert record[0]["ScientificName"] == "Canis lupus familiaris"
        assert record[0]["OtherNames"]["GenbankCommonName"] == "dog"
        assert len(record[0]["OtherNames"]["Synonym"]) == 5
        assert record[0]["OtherNames"]["Synonym"][0] == "Canis borealis"
        assert record[0]["OtherNames"]["Synonym"][1] == "Canis canis"
        assert record[0]["OtherNames"]["Synonym"][2] == "Canis domesticus"
        assert record[0]["OtherNames"]["Synonym"][3] == "Canis familiaris"
        assert record[0]["OtherNames"]["Synonym"][4] == "Canis lupus borealis"
        assert len(record[0]["OtherNames"]["CommonName"]) == 1
        assert record[0]["OtherNames"]["CommonName"][0] == "dogs"
        assert record[0]["OtherNames"]["Includes"][0] == "beagle dog"
        assert record[0]["OtherNames"]["Includes"][1] == "beagle dogs"
        assert record[0]["ParentTaxId"] == "9612"
        assert record[0]["Rank"] == "subspecies"
        assert record[0]["Division"] == "Mammals"
        assert record[0]["GeneticCode"]["GCId"] == "1"
        assert record[0]["GeneticCode"]["GCName"] == "Standard"
        assert record[0]["MitoGeneticCode"]["MGCId"] == "2"
        assert record[0]["MitoGeneticCode"]["MGCName"] == "Vertebrate Mitochondrial"
        assert record[0]["Lineage"] == "cellular organisms; Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Deuterostomia; Chordata; Craniata; Vertebrata; Gnathostomata; Teleostomi; Euteleostomi; Sarcopterygii; Dipnotetrapodomorpha; Tetrapoda; Amniota; Mammalia; Theria; Eutheria; Boreoeutheria; Laurasiatheria; Carnivora; Caniformia; Canidae; Canis; Canis lupus"

        assert record[0]["LineageEx"][0]["TaxId"] == "131567"
        assert record[0]["LineageEx"][0]["ScientificName"] == "cellular organisms"
        assert record[0]["LineageEx"][0]["Rank"] == "cellular root"
        assert record[0]["LineageEx"][1]["TaxId"] == "2759"
        assert record[0]["LineageEx"][1]["ScientificName"] == "Eukaryota"
        assert record[0]["LineageEx"][1]["Rank"] == "domain"
        assert record[0]["LineageEx"][2]["TaxId"] == "33154"
        assert record[0]["LineageEx"][2]["ScientificName"] == "Opisthokonta"
        assert record[0]["LineageEx"][2]["Rank"] == "clade"
        assert record[0]["LineageEx"][3]["TaxId"] == "33208"
        assert record[0]["LineageEx"][3]["ScientificName"] == "Metazoa"
        assert record[0]["LineageEx"][3]["Rank"] == "kingdom"
        assert record[0]["LineageEx"][4]["TaxId"] == "6072"
        assert record[0]["LineageEx"][4]["ScientificName"] == "Eumetazoa"
        assert record[0]["LineageEx"][4]["Rank"] == "clade"
        assert record[0]["LineageEx"][5]["TaxId"] == "33213"
        assert record[0]["LineageEx"][5]["ScientificName"] == "Bilateria"
        assert record[0]["LineageEx"][5]["Rank"] == "clade"
        assert record[0]["LineageEx"][6]["TaxId"] == "33511"
        assert record[0]["LineageEx"][6]["ScientificName"] == "Deuterostomia"
        assert record[0]["LineageEx"][6]["Rank"] == "clade"
        assert record[0]["LineageEx"][7]["TaxId"] == "7711"
        assert record[0]["LineageEx"][7]["ScientificName"] == "Chordata"
        assert record[0]["LineageEx"][7]["Rank"] == "phylum"
        assert record[0]["LineageEx"][8]["TaxId"] == "89593"
        assert record[0]["LineageEx"][8]["ScientificName"] == "Craniata"
        assert record[0]["LineageEx"][8]["Rank"] == "subphylum"
        assert record[0]["LineageEx"][9]["TaxId"] == "7742"
        assert record[0]["LineageEx"][9]["ScientificName"] == "Vertebrata"
        assert record[0]["LineageEx"][9]["Rank"] == "clade"
        assert record[0]["LineageEx"][10]["TaxId"] == "7776"
        assert record[0]["LineageEx"][10]["ScientificName"] == "Gnathostomata"
        assert record[0]["LineageEx"][10]["Rank"] == "clade"
        assert record[0]["LineageEx"][11]["TaxId"] == "117570"
        assert record[0]["LineageEx"][11]["ScientificName"] == "Teleostomi"
        assert record[0]["LineageEx"][11]["Rank"] == "clade"
        assert record[0]["LineageEx"][12]["TaxId"] == "117571"
        assert record[0]["LineageEx"][12]["ScientificName"] == "Euteleostomi"
        assert record[0]["LineageEx"][12]["Rank"] == "clade"
        assert record[0]["LineageEx"][13]["TaxId"] == "8287"
        assert record[0]["LineageEx"][13]["ScientificName"] == "Sarcopterygii"
        assert record[0]["LineageEx"][13]["Rank"] == "superclass"
        assert record[0]["LineageEx"][14]["TaxId"] == "1338369"
        assert record[0]["LineageEx"][14]["ScientificName"] == "Dipnotetrapodomorpha"
        assert record[0]["LineageEx"][14]["Rank"] == "clade"
        assert record[0]["LineageEx"][15]["TaxId"] == "32523"
        assert record[0]["LineageEx"][15]["ScientificName"] == "Tetrapoda"
        assert record[0]["LineageEx"][15]["Rank"] == "clade"
        assert record[0]["LineageEx"][16]["TaxId"] == "32524"
        assert record[0]["LineageEx"][16]["ScientificName"] == "Amniota"
        assert record[0]["LineageEx"][16]["Rank"] == "clade"
        assert record[0]["LineageEx"][17]["TaxId"] == "40674"
        assert record[0]["LineageEx"][17]["ScientificName"] == "Mammalia"
        assert record[0]["LineageEx"][17]["Rank"] == "class"
        assert record[0]["LineageEx"][18]["TaxId"] == "32525"
        assert record[0]["LineageEx"][18]["ScientificName"] == "Theria"
        assert record[0]["LineageEx"][18]["Rank"] == "clade"
        assert record[0]["LineageEx"][19]["TaxId"] == "9347"
        assert record[0]["LineageEx"][19]["ScientificName"] == "Eutheria"
        assert record[0]["LineageEx"][19]["Rank"] == "clade"
        assert record[0]["LineageEx"][20]["TaxId"] == "1437010"
        assert record[0]["LineageEx"][20]["ScientificName"] == "Boreoeutheria"
        assert record[0]["LineageEx"][20]["Rank"] == "clade"
        assert record[0]["LineageEx"][21]["TaxId"] == "314145"
        assert record[0]["LineageEx"][21]["ScientificName"] == "Laurasiatheria"
        assert record[0]["LineageEx"][21]["Rank"] == "superorder"
        assert record[0]["LineageEx"][22]["TaxId"] == "33554"
        assert record[0]["LineageEx"][22]["ScientificName"] == "Carnivora"
        assert record[0]["LineageEx"][22]["Rank"] == "order"
        assert record[0]["LineageEx"][23]["TaxId"] == "379584"
        assert record[0]["LineageEx"][23]["ScientificName"] == "Caniformia"
        assert record[0]["LineageEx"][23]["Rank"] == "suborder"
        assert record[0]["LineageEx"][24]["TaxId"] == "9608"
        assert record[0]["LineageEx"][24]["ScientificName"] == "Canidae"
        assert record[0]["LineageEx"][24]["Rank"] == "family"
        assert record[0]["LineageEx"][25]["TaxId"] == "9611"
        assert record[0]["LineageEx"][25]["ScientificName"] == "Canis"
        assert record[0]["LineageEx"][25]["Rank"] == "genus"
        assert record[0]["LineageEx"][26]["TaxId"] == "9612"
        assert record[0]["LineageEx"][26]["ScientificName"] == "Canis lupus"
        assert record[0]["LineageEx"][26]["Rank"] == "species"
        assert record[0]["CreateDate"] == "1995/02/27 09:24:00"
        assert record[0]["UpdateDate"] == "2024/02/09 13:11:20"
        assert record[0]["PubDate"] == "1993/04/27 01:00:00"
        assert record[1]["TaxId"] == "9685"
        assert record[1]["ScientificName"] == "Felis catus"
        assert record[1]["OtherNames"]["GenbankCommonName"] == "domestic cat"
        assert len(record[1]["OtherNames"]["Synonym"]) == 2
        assert record[1]["OtherNames"]["Synonym"][0] == "Felis domesticus"
        assert record[1]["OtherNames"]["Synonym"][1] == "Felis silvestris catus"
        assert len(record[1]["OtherNames"]["CommonName"]) == 2
        assert record[1]["OtherNames"]["CommonName"][0] == "cat"
        assert record[1]["OtherNames"]["CommonName"][1] == "cats"
        assert record[1]["OtherNames"]["Includes"][0] == "Korat cats"
        assert record[1]["ParentTaxId"] == "9682"
        assert record[1]["Rank"] == "species"
        assert record[1]["Division"] == "Mammals"
        assert record[1]["GeneticCode"]["GCId"] == "1"
        assert record[1]["GeneticCode"]["GCName"] == "Standard"
        assert record[1]["MitoGeneticCode"]["MGCId"] == "2"
        assert record[1]["MitoGeneticCode"]["MGCName"] == "Vertebrate Mitochondrial"
        assert record[1]["Lineage"] == "cellular organisms; Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Deuterostomia; Chordata; Craniata; Vertebrata; Gnathostomata; Teleostomi; Euteleostomi; Sarcopterygii; Dipnotetrapodomorpha; Tetrapoda; Amniota; Mammalia; Theria; Eutheria; Boreoeutheria; Laurasiatheria; Carnivora; Feliformia; Felidae; Felinae; Felis"

        assert record[1]["LineageEx"][0]["TaxId"] == "131567"
        assert record[1]["LineageEx"][0]["ScientificName"] == "cellular organisms"
        assert record[1]["LineageEx"][0]["Rank"] == "cellular root"
        assert record[1]["LineageEx"][1]["TaxId"] == "2759"
        assert record[1]["LineageEx"][1]["ScientificName"] == "Eukaryota"
        assert record[1]["LineageEx"][1]["Rank"] == "domain"
        assert record[1]["LineageEx"][2]["TaxId"] == "33154"
        assert record[1]["LineageEx"][2]["ScientificName"] == "Opisthokonta"
        assert record[1]["LineageEx"][2]["Rank"] == "clade"
        assert record[1]["LineageEx"][3]["TaxId"] == "33208"
        assert record[1]["LineageEx"][3]["ScientificName"] == "Metazoa"
        assert record[1]["LineageEx"][3]["Rank"] == "kingdom"
        assert record[1]["LineageEx"][4]["TaxId"] == "6072"
        assert record[1]["LineageEx"][4]["ScientificName"] == "Eumetazoa"
        assert record[1]["LineageEx"][4]["Rank"] == "clade"
        assert record[1]["LineageEx"][5]["TaxId"] == "33213"
        assert record[1]["LineageEx"][5]["ScientificName"] == "Bilateria"
        assert record[1]["LineageEx"][5]["Rank"] == "clade"
        assert record[1]["LineageEx"][6]["TaxId"] == "33511"
        assert record[1]["LineageEx"][6]["ScientificName"] == "Deuterostomia"
        assert record[1]["LineageEx"][6]["Rank"] == "clade"
        assert record[1]["LineageEx"][7]["TaxId"] == "7711"
        assert record[1]["LineageEx"][7]["ScientificName"] == "Chordata"
        assert record[1]["LineageEx"][7]["Rank"] == "phylum"
        assert record[1]["LineageEx"][8]["TaxId"] == "89593"
        assert record[1]["LineageEx"][8]["ScientificName"] == "Craniata"
        assert record[1]["LineageEx"][8]["Rank"] == "subphylum"
        assert record[1]["LineageEx"][9]["TaxId"] == "7742"
        assert record[1]["LineageEx"][9]["ScientificName"] == "Vertebrata"
        assert record[1]["LineageEx"][9]["Rank"] == "clade"
        assert record[1]["LineageEx"][10]["TaxId"] == "7776"
        assert record[1]["LineageEx"][10]["ScientificName"] == "Gnathostomata"
        assert record[1]["LineageEx"][10]["Rank"] == "clade"
        assert record[1]["LineageEx"][11]["TaxId"] == "117570"
        assert record[1]["LineageEx"][11]["ScientificName"] == "Teleostomi"
        assert record[1]["LineageEx"][11]["Rank"] == "clade"
        assert record[1]["LineageEx"][12]["TaxId"] == "117571"
        assert record[1]["LineageEx"][12]["ScientificName"] == "Euteleostomi"
        assert record[1]["LineageEx"][12]["Rank"] == "clade"
        assert record[1]["LineageEx"][13]["TaxId"] == "8287"
        assert record[1]["LineageEx"][13]["ScientificName"] == "Sarcopterygii"
        assert record[1]["LineageEx"][13]["Rank"] == "superclass"
        assert record[1]["LineageEx"][14]["TaxId"] == "1338369"
        assert record[1]["LineageEx"][14]["ScientificName"] == "Dipnotetrapodomorpha"
        assert record[1]["LineageEx"][14]["Rank"] == "clade"
        assert record[1]["LineageEx"][15]["TaxId"] == "32523"
        assert record[1]["LineageEx"][15]["ScientificName"] == "Tetrapoda"
        assert record[1]["LineageEx"][15]["Rank"] == "clade"
        assert record[1]["LineageEx"][16]["TaxId"] == "32524"
        assert record[1]["LineageEx"][16]["ScientificName"] == "Amniota"
        assert record[1]["LineageEx"][16]["Rank"] == "clade"
        assert record[1]["LineageEx"][17]["TaxId"] == "40674"
        assert record[1]["LineageEx"][17]["ScientificName"] == "Mammalia"
        assert record[1]["LineageEx"][17]["Rank"] == "class"
        assert record[1]["LineageEx"][18]["TaxId"] == "32525"
        assert record[1]["LineageEx"][18]["ScientificName"] == "Theria"
        assert record[1]["LineageEx"][18]["Rank"] == "clade"
        assert record[1]["LineageEx"][19]["TaxId"] == "9347"
        assert record[1]["LineageEx"][19]["ScientificName"] == "Eutheria"
        assert record[1]["LineageEx"][19]["Rank"] == "clade"
        assert record[1]["LineageEx"][20]["TaxId"] == "1437010"
        assert record[1]["LineageEx"][20]["ScientificName"] == "Boreoeutheria"
        assert record[1]["LineageEx"][20]["Rank"] == "clade"
        assert record[1]["LineageEx"][21]["TaxId"] == "314145"
        assert record[1]["LineageEx"][21]["ScientificName"] == "Laurasiatheria"
        assert record[1]["LineageEx"][21]["Rank"] == "superorder"
        assert record[1]["LineageEx"][22]["TaxId"] == "33554"
        assert record[1]["LineageEx"][22]["ScientificName"] == "Carnivora"
        assert record[1]["LineageEx"][22]["Rank"] == "order"
        assert record[1]["LineageEx"][23]["TaxId"] == "379583"
        assert record[1]["LineageEx"][23]["ScientificName"] == "Feliformia"
        assert record[1]["LineageEx"][23]["Rank"] == "suborder"
        assert record[1]["LineageEx"][24]["TaxId"] == "9681"
        assert record[1]["LineageEx"][24]["ScientificName"] == "Felidae"
        assert record[1]["LineageEx"][24]["Rank"] == "family"
        assert record[1]["LineageEx"][25]["TaxId"] == "338152"
        assert record[1]["LineageEx"][25]["ScientificName"] == "Felinae"
        assert record[1]["LineageEx"][25]["Rank"] == "subfamily"
        assert record[1]["LineageEx"][26]["TaxId"] == "9682"
        assert record[1]["LineageEx"][26]["ScientificName"] == "Felis"
        assert record[1]["LineageEx"][26]["Rank"] == "genus"
        assert record[1]["CreateDate"] == "1995/02/27 09:24:00"
        assert record[1]["UpdateDate"] == "2024/03/03 11:27:08"
        assert record[1]["PubDate"] == "1993/07/26 01:00:00"
        # fmt: on

    def test_nucleotide1(self):
        """Test parsing XML returned by EFetch, Nucleotide database (first test)."""
        # Access the nucleotide database using efetch.
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db='nucleotide', id=5, retmode='xml')
        with open("Entrez/nucleotide1.xml", "rb") as stream:
            record = Entrez.read(stream)
        # fmt: off
        assert record[0]["GBSeq_locus"] == "X60065"
        assert record[0]["GBSeq_length"] == "1136"
        assert record[0]["GBSeq_strandedness"] == "single"
        assert record[0]["GBSeq_moltype"] == "mRNA"
        assert record[0]["GBSeq_topology"] == "linear"
        assert record[0]["GBSeq_division"] == "MAM"
        assert record[0]["GBSeq_update-date"] == "26-JUL-2016"
        assert record[0]["GBSeq_create-date"] == "05-MAY-1992"
        assert record[0]["GBSeq_definition"] == "B.bovis beta-2-gpI mRNA for beta-2-glycoprotein I"
        assert record[0]["GBSeq_primary-accession"] == "X60065"
        assert record[0]["GBSeq_accession-version"] == "X60065.1"
        assert record[0]["GBSeq_other-seqids"][0] == "emb|X60065.1|"
        assert record[0]["GBSeq_other-seqids"][1] == "gi|5"
        assert record[0]["GBSeq_keywords"][0] == "beta-2 glycoprotein I"
        assert record[0]["GBSeq_source"] == "Bos taurus (domestic cattle)"
        assert record[0]["GBSeq_organism"] == "Bos taurus"
        assert record[0]["GBSeq_taxonomy"] == "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Laurasiatheria; Artiodactyla; Ruminantia; Pecora; Bovidae; Bovinae; Bos"
        assert record[0]["GBSeq_references"][0]["GBReference_reference"] == "1"
        assert record[0]["GBSeq_references"][0]["GBReference_authors"][0] == "Bendixen,E."
        assert record[0]["GBSeq_references"][0]["GBReference_authors"][1] == "Halkier,T."
        assert record[0]["GBSeq_references"][0]["GBReference_authors"][2] == "Magnusson,S."
        assert record[0]["GBSeq_references"][0]["GBReference_authors"][3] == "Sottrup-Jensen,L."
        assert record[0]["GBSeq_references"][0]["GBReference_authors"][4] == "Kristensen,T."
        assert record[0]["GBSeq_references"][0]["GBReference_title"] == "Complete primary structure of bovine beta 2-glycoprotein I: localization of the disulfide bridges"
        assert record[0]["GBSeq_references"][0]["GBReference_journal"] == "Biochemistry 31 (14), 3611-3617 (1992)"
        assert record[0]["GBSeq_references"][0]["GBReference_pubmed"] == "1567819"
        assert record[0]["GBSeq_references"][1]["GBReference_reference"] == "2"
        assert record[0]["GBSeq_references"][1]["GBReference_position"] == "1..1136"
        assert record[0]["GBSeq_references"][1]["GBReference_authors"][0] == "Kristensen,T."
        assert record[0]["GBSeq_references"][1]["GBReference_title"] == "Direct Submission"
        assert record[0]["GBSeq_references"][1]["GBReference_journal"] == "Submitted (11-JUN-1991) T. Kristensen, Dept of Mol Biology, University of Aarhus, C F Mollers Alle 130, DK-8000 Aarhus C, DENMARK"
        assert len(record[0]["GBSeq_feature-table"]) == 7
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_key"] == "source"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_location"] == "1..1136"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_intervals"][0]["GBInterval_from"] == "1"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_intervals"][0]["GBInterval_to"] == "1136"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_intervals"][0]["GBInterval_accession"] == "X60065.1"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][0]["GBQualifier_name"] == "organism"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][0]["GBQualifier_value"] == "Bos taurus"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][1]["GBQualifier_name"] == "mol_type"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][1]["GBQualifier_value"] == "mRNA"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][2]["GBQualifier_name"] == "db_xref"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][2]["GBQualifier_value"] == "taxon:9913"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][3]["GBQualifier_name"] == "clone"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][3]["GBQualifier_value"] == "pBB2I"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][4]["GBQualifier_name"] == "tissue_type"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][4]["GBQualifier_value"] == "liver"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_key"] == "gene"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_location"] == "<1..1136"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_intervals"][0]["GBInterval_from"] == "1"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_intervals"][0]["GBInterval_to"] == "1136"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_intervals"][0]["GBInterval_accession"] == "X60065.1"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_partial5"] == ""
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_partial5"].attributes["value"] == "true"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_quals"][0]["GBQualifier_name"] == "gene"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_quals"][0]["GBQualifier_value"] == "beta-2-gpI"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_key"] == "CDS"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_location"] == "<1..1029"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_intervals"][0]["GBInterval_from"] == "1"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_intervals"][0]["GBInterval_to"] == "1029"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_intervals"][0]["GBInterval_accession"] == "X60065.1"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_partial5"] == ""
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_partial5"].attributes["value"] == "true"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][0]["GBQualifier_name"] == "gene"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][0]["GBQualifier_value"] == "beta-2-gpI"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][1]["GBQualifier_name"] == "codon_start"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][1]["GBQualifier_value"] == "1"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][2]["GBQualifier_name"] == "transl_table"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][2]["GBQualifier_value"] == "1"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][3]["GBQualifier_name"] == "product"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][3]["GBQualifier_value"] == "beta-2-glycoprotein I"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][4]["GBQualifier_name"] == "protein_id"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][4]["GBQualifier_value"] == "CAA42669.1"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][5]["GBQualifier_name"] == "db_xref"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][5]["GBQualifier_value"] == "GOA:P17690"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][6]["GBQualifier_name"] == "db_xref"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][6]["GBQualifier_value"] == "InterPro:IPR000436"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][7]["GBQualifier_name"] == "db_xref"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][7]["GBQualifier_value"] == "InterPro:IPR015104"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][8]["GBQualifier_name"] == "db_xref"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][8]["GBQualifier_value"] == "InterPro:IPR016060"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_key"] == "sig_peptide"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_location"] == "<1..48"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_intervals"][0]["GBInterval_from"] == "1"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_intervals"][0]["GBInterval_to"] == "48"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_intervals"][0]["GBInterval_accession"] == "X60065.1"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_partial5"] == ""
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_partial5"].attributes["value"] == "true"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][0]["GBQualifier_name"] == "gene"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][0]["GBQualifier_value"] == "beta-2-gpI"
        assert record[0]["GBSeq_feature-table"][4]["GBFeature_key"] == "mat_peptide"
        assert record[0]["GBSeq_feature-table"][4]["GBFeature_location"] == "49..1026"
        assert record[0]["GBSeq_feature-table"][4]["GBFeature_intervals"][0]["GBInterval_from"] == "49"
        assert record[0]["GBSeq_feature-table"][4]["GBFeature_intervals"][0]["GBInterval_to"] == "1026"
        assert record[0]["GBSeq_feature-table"][4]["GBFeature_intervals"][0]["GBInterval_accession"] == "X60065.1"
        assert record[0]["GBSeq_feature-table"][4]["GBFeature_quals"][0]["GBQualifier_name"] == "gene"
        assert record[0]["GBSeq_feature-table"][4]["GBFeature_quals"][0]["GBQualifier_value"] == "beta-2-gpI"
        assert record[0]["GBSeq_feature-table"][4]["GBFeature_quals"][1]["GBQualifier_name"] == "product"
        assert record[0]["GBSeq_feature-table"][4]["GBFeature_quals"][1]["GBQualifier_value"] == "beta-2-glycoprotein I"
        assert record[0]["GBSeq_feature-table"][4]["GBFeature_quals"][2]["GBQualifier_name"] == "peptide"
        assert record[0]["GBSeq_feature-table"][4]["GBFeature_quals"][2]["GBQualifier_value"] == "GRTCPKPDELPFSTVVPLKRTYEPGEQIVFSCQPGYVSRGGIRRFTCPLTGLWPINTLKCMPRVCPFAGILENGTVRYTTFEYPNTISFSCHTGFYLKGASSAKCTEEGKWSPDLPVCAPITCPPPPIPKFASLSVYKPLAGNNSFYGSKAVFKCLPHHAMFGNDTVTCTEHGNWTQLPECREVRCPFPSRPDNGFVNHPANPVLYYKDTATFGCHETYSLDGPEEVECSKFGNWSAQPSCKASCKLSIKRATVIYEGERVAIQNKFKNGMLHGQKVSFFCKHKEKKCSYTEDAQCIDGTIEIPKCFKEHSSLAFWKTDASDVKPC"
        assert record[0]["GBSeq_feature-table"][5]["GBFeature_key"] == "regulatory"
        assert record[0]["GBSeq_feature-table"][5]["GBFeature_location"] == "1101..1106"
        assert record[0]["GBSeq_feature-table"][5]["GBFeature_intervals"][0]["GBInterval_from"] == "1101"
        assert record[0]["GBSeq_feature-table"][5]["GBFeature_intervals"][0]["GBInterval_to"] == "1106"
        assert record[0]["GBSeq_feature-table"][5]["GBFeature_intervals"][0]["GBInterval_accession"] == "X60065.1"
        assert record[0]["GBSeq_feature-table"][5]["GBFeature_quals"][0]["GBQualifier_name"] == "regulatory_class"
        assert record[0]["GBSeq_feature-table"][5]["GBFeature_quals"][0]["GBQualifier_value"] == "polyA_signal_sequence"
        assert record[0]["GBSeq_feature-table"][5]["GBFeature_quals"][1]["GBQualifier_name"] == "gene"
        assert record[0]["GBSeq_feature-table"][5]["GBFeature_quals"][1]["GBQualifier_value"] == "beta-2-gpI"
        assert record[0]["GBSeq_feature-table"][6]["GBFeature_key"] == "polyA_site"
        assert record[0]["GBSeq_feature-table"][6]["GBFeature_location"] == "1130"
        assert record[0]["GBSeq_feature-table"][6]["GBFeature_intervals"][0]["GBInterval_point"] == "1130"
        assert record[0]["GBSeq_feature-table"][6]["GBFeature_intervals"][0]["GBInterval_accession"] == "X60065.1"
        assert record[0]["GBSeq_feature-table"][6]["GBFeature_quals"][0]["GBQualifier_name"] == "gene"
        assert record[0]["GBSeq_feature-table"][6]["GBFeature_quals"][0]["GBQualifier_value"] == "beta-2-gpI"
        assert record[0]["GBSeq_sequence"] == "ccagcgctcgtcttgctgttggggtttctctgccacgttgctatcgcaggacgaacctgccccaagccagatgagctaccgttttccacggtggttccactgaaacggacctatgagcccggggagcagatagtcttctcctgccagccgggctacgtgtcccggggagggatccggcggtttacatgcccgctcacaggactctggcccatcaacacgctgaaatgcatgcccagagtatgtccttttgctgggatcttagaaaacggaacggtacgctatacaacgtttgagtatcccaacaccatcagcttttcttgccacacggggttttatctgaaaggagctagttctgcaaaatgcactgaggaagggaagtggagcccagaccttcctgtctgtgcccctataacctgccctccaccacccatacccaagtttgcaagtctcagcgtttacaagccgttggctgggaacaactccttctatggcagcaaggcagtctttaagtgcttgccacaccacgcgatgtttggaaatgacaccgttacctgcacggaacatgggaactggacgcagttgccagaatgcagggaagtaagatgcccattcccatcaagaccagacaatgggtttgtgaaccatcctgcaaatccagtgctctactataaggacaccgccacctttggctgccatgaaacgtattccttggatggaccggaagaagtagaatgcagcaaattcggaaactggtctgcacagccaagctgtaaagcatcttgtaagttatctattaaaagagctactgtgatatatgaaggagagagagtagctatccagaacaaatttaagaatggaatgctgcatggccaaaaggtttctttcttctgcaagcataaggaaaagaagtgcagctacacagaagatgctcagtgcatagacggcaccatcgagattcccaaatgcttcaaggagcacagttctttagctttctggaaaacggatgcatctgacgtaaaaccatgctaagctggttttcacactgaaaattaaatgtcatgcttatatgtgtctgtctgagaatctgatggaaacggaaaaataaagagactgaatttaccgtgtcaagaaaaaaa"
        # fmt: off

    def test_nucleotide2(self):
        """Test parsing XML returned by EFetch, Nucleotide database (second test)."""
        # Access the nucleotide database using efetch.
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db='nucleotide', id=5,
        #                       rettype='fasta', complexity=0, retmode='xml')
        with open("Entrez/nucleotide2.xml", "rb") as stream:
            record = Entrez.read(stream)
        assert record[0]["TSeq_seqtype"] == ""
        assert record[0]["TSeq_seqtype"].attributes["value"] == "nucleotide"
        assert record[0]["TSeq_accver"] == "X60065.1"
        assert record[0]["TSeq_taxid"] == "9913"
        assert record[0]["TSeq_orgname"] == "Bos taurus"
        assert record[0]["TSeq_defline"] == "B.bovis beta-2-gpI mRNA for beta-2-glycoprotein I"
        assert record[0]["TSeq_length"] == "1136"
        assert record[0]["TSeq_sequence"] == "CCAGCGCTCGTCTTGCTGTTGGGGTTTCTCTGCCACGTTGCTATCGCAGGACGAACCTGCCCCAAGCCAGATGAGCTACCGTTTTCCACGGTGGTTCCACTGAAACGGACCTATGAGCCCGGGGAGCAGATAGTCTTCTCCTGCCAGCCGGGCTACGTGTCCCGGGGAGGGATCCGGCGGTTTACATGCCCGCTCACAGGACTCTGGCCCATCAACACGCTGAAATGCATGCCCAGAGTATGTCCTTTTGCTGGGATCTTAGAAAACGGAACGGTACGCTATACAACGTTTGAGTATCCCAACACCATCAGCTTTTCTTGCCACACGGGGTTTTATCTGAAAGGAGCTAGTTCTGCAAAATGCACTGAGGAAGGGAAGTGGAGCCCAGACCTTCCTGTCTGTGCCCCTATAACCTGCCCTCCACCACCCATACCCAAGTTTGCAAGTCTCAGCGTTTACAAGCCGTTGGCTGGGAACAACTCCTTCTATGGCAGCAAGGCAGTCTTTAAGTGCTTGCCACACCACGCGATGTTTGGAAATGACACCGTTACCTGCACGGAACATGGGAACTGGACGCAGTTGCCAGAATGCAGGGAAGTAAGATGCCCATTCCCATCAAGACCAGACAATGGGTTTGTGAACCATCCTGCAAATCCAGTGCTCTACTATAAGGACACCGCCACCTTTGGCTGCCATGAAACGTATTCCTTGGATGGACCGGAAGAAGTAGAATGCAGCAAATTCGGAAACTGGTCTGCACAGCCAAGCTGTAAAGCATCTTGTAAGTTATCTATTAAAAGAGCTACTGTGATATATGAAGGAGAGAGAGTAGCTATCCAGAACAAATTTAAGAATGGAATGCTGCATGGCCAAAAGGTTTCTTTCTTCTGCAAGCATAAGGAAAAGAAGTGCAGCTACACAGAAGATGCTCAGTGCATAGACGGCACCATCGAGATTCCCAAATGCTTCAAGGAGCACAGTTCTTTAGCTTTCTGGAAAACGGATGCATCTGACGTAAAACCATGCTAAGCTGGTTTTCACACTGAAAATTAAATGTCATGCTTATATGTGTCTGTCTGAGAATCTGATGGAAACGGAAAAATAAAGAGACTGAATTTACCGTGTCAAGAAAAAAA"
        assert record[1]["TSeq_seqtype"] == ""
        assert record[1]["TSeq_seqtype"].attributes["value"] == "protein"
        assert record[1]["TSeq_accver"] == "CAA42669.1"
        assert record[1]["TSeq_taxid"] == "9913"
        assert record[1]["TSeq_orgname"] == "Bos taurus"
        assert record[1]["TSeq_defline"] == "beta-2-glycoprotein I, partial [Bos taurus]"
        assert record[1]["TSeq_length"] == "342"
        assert record[1]["TSeq_sequence"] == "PALVLLLGFLCHVAIAGRTCPKPDELPFSTVVPLKRTYEPGEQIVFSCQPGYVSRGGIRRFTCPLTGLWPINTLKCMPRVCPFAGILENGTVRYTTFEYPNTISFSCHTGFYLKGASSAKCTEEGKWSPDLPVCAPITCPPPPIPKFASLSVYKPLAGNNSFYGSKAVFKCLPHHAMFGNDTVTCTEHGNWTQLPECREVRCPFPSRPDNGFVNHPANPVLYYKDTATFGCHETYSLDGPEEVECSKFGNWSAQPSCKASCKLSIKRATVIYEGERVAIQNKFKNGMLHGQKVSFFCKHKEKKCSYTEDAQCIDGTIEIPKCFKEHSSLAFWKTDASDVKPC"

    def test_protein(self):
        """Test parsing XML returned by EFetch, Protein database."""
        # Access the protein database using efetch.
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db='protein', id=8, rettype='gp', retmode='xml')
        with open("Entrez/protein.xml", "rb") as stream:
            record = Entrez.read(stream)
        # fmt: off
        assert record[0]["GBSeq_locus"] == "CAA35997"
        assert record[0]["GBSeq_length"] == "100"
        assert record[0]["GBSeq_moltype"] == "AA"
        assert record[0]["GBSeq_topology"] == "linear"
        assert record[0]["GBSeq_division"] == "MAM"
        assert record[0]["GBSeq_update-date"] == "12-SEP-1993"
        assert record[0]["GBSeq_create-date"] == "03-APR-1990"
        assert record[0]["GBSeq_definition"] == "unnamed protein product [Bos taurus]"
        assert record[0]["GBSeq_primary-accession"] == "CAA35997"
        assert record[0]["GBSeq_accession-version"] == "CAA35997.1"
        assert record[0]["GBSeq_other-seqids"][0] == "emb|CAA35997.1|"
        assert record[0]["GBSeq_other-seqids"][1] == "gi|8"
        assert record[0]["GBSeq_source"] == "Bos taurus (domestic cattle)"
        assert record[0]["GBSeq_organism"] == "Bos taurus"
        assert record[0]["GBSeq_taxonomy"] == "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Laurasiatheria; Artiodactyla; Ruminantia; Pecora; Bovidae; Bovinae; Bos"
        assert record[0]["GBSeq_references"][0]["GBReference_reference"] == "1"
        assert record[0]["GBSeq_references"][0]["GBReference_position"] == "1..100"
        assert record[0]["GBSeq_references"][0]["GBReference_authors"][0] == "Kiefer,M.C."
        assert record[0]["GBSeq_references"][0]["GBReference_authors"][1] == "Saphire,A.C.S."
        assert record[0]["GBSeq_references"][0]["GBReference_authors"][2] == "Bauer,D.M."
        assert record[0]["GBSeq_references"][0]["GBReference_authors"][3] == "Barr,P.J."
        assert record[0]["GBSeq_references"][0]["GBReference_journal"] == "Unpublished"
        assert record[0]["GBSeq_references"][1]["GBReference_reference"] == "2"
        assert record[0]["GBSeq_references"][1]["GBReference_position"] == "1..100"
        assert record[0]["GBSeq_references"][1]["GBReference_authors"][0] == "Kiefer,M.C."
        assert record[0]["GBSeq_references"][1]["GBReference_title"] == "Direct Submission"
        assert record[0]["GBSeq_references"][1]["GBReference_journal"] == "Submitted (30-JAN-1990) Kiefer M.C., Chiron Corporation, 4560 Hortom St, Emeryville CA 94608-2916, U S A"
        assert record[0]["GBSeq_comment"] == "See <X15699> for Human sequence.~~Data kindly reviewed (08-MAY-1990) by Kiefer M.C."
        assert record[0]["GBSeq_source-db"] == "embl accession X51700.1"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_key"] == "source"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_location"] == "1..100"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_intervals"][0]["GBInterval_from"] == "1"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_intervals"][0]["GBInterval_to"] == "100"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_intervals"][0]["GBInterval_accession"] == "CAA35997.1"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][0]["GBQualifier_name"] == "organism"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][0]["GBQualifier_value"] == "Bos taurus"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][1]["GBQualifier_name"] == "db_xref"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][1]["GBQualifier_value"] == "taxon:9913"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][2]["GBQualifier_name"] == "clone"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][2]["GBQualifier_value"] == "bBGP-3"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][3]["GBQualifier_name"] == "tissue_type"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][3]["GBQualifier_value"] == "bone matrix"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][4]["GBQualifier_name"] == "clone_lib"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][4]["GBQualifier_value"] == "Zap-bb"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_key"] == "Protein"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_location"] == "1..100"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_intervals"][0]["GBInterval_from"] == "1"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_intervals"][0]["GBInterval_to"] == "100"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_intervals"][0]["GBInterval_accession"] == "CAA35997.1"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_quals"][0]["GBQualifier_name"] == "name"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_quals"][0]["GBQualifier_value"] == "unnamed protein product"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_key"] == "Region"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_location"] == "33..97"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_intervals"][0]["GBInterval_from"] == "33"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_intervals"][0]["GBInterval_to"] == "97"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_intervals"][0]["GBInterval_accession"] == "CAA35997.1"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][0]["GBQualifier_name"] == "region_name"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][0]["GBQualifier_value"] == "GLA"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][1]["GBQualifier_name"] == "note"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][1]["GBQualifier_value"] == "Domain containing Gla (gamma-carboxyglutamate) residues; smart00069"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][2]["GBQualifier_name"] == "db_xref"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][2]["GBQualifier_value"] == "CDD:214503"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_key"] == "CDS"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_location"] == "1..100"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_intervals"][0]["GBInterval_from"] == "1"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_intervals"][0]["GBInterval_to"] == "100"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_intervals"][0]["GBInterval_accession"] == "CAA35997.1"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][0]["GBQualifier_name"] == "coded_by"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][0]["GBQualifier_value"] == "X51700.1:28..330"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][1]["GBQualifier_name"] == "note"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][1]["GBQualifier_value"] == "bone Gla precursor (100 AA)"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][2]["GBQualifier_name"] == "transl_table"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][2]["GBQualifier_value"] == "1"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][3]["GBQualifier_name"] == "db_xref"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][3]["GBQualifier_value"] == "GOA:P02820"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][4]["GBQualifier_name"] == "db_xref"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][4]["GBQualifier_value"] == "InterPro:IPR000294"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][5]["GBQualifier_name"] == "db_xref"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][5]["GBQualifier_value"] == "InterPro:IPR002384"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][6]["GBQualifier_name"] == "db_xref"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][6]["GBQualifier_value"] == "PDB:1Q3M"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][7]["GBQualifier_name"] == "db_xref"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][7]["GBQualifier_value"] == "UniProtKB/Swiss-Prot:P02820"
        assert record[0]["GBSeq_sequence"] == "mrtpmllallalatlclagradakpgdaesgkgaafvskqegsevvkrlrryldhwlgapapypdplepkrevcelnpdcdeladhigfqeayrrfygpv"
        # fmt: on

    def test_efetch_schemas(self):
        """Test parsing XML using Schemas."""
        # To create the XML file,use
        # >>> Bio.Entrez.efetch("protein", id="783730874", rettype="ipg", retmode="xml")
        with open("Entrez/efetch_schemas.xml", "rb") as stream:
            records = Entrez.read(stream)
        assert len(records) == 1
        record = records["IPGReport"]
        assert len(record.attributes) == 2
        assert record.attributes["product_acc"] == "KJV04014.1"
        assert record.attributes["ipg"] == "79092155"
        assert len(record) == 3
        assert record["Product"] == ""
        assert record["Product"].attributes["kingdom"] == "Bacteria"
        assert record["Product"].attributes["slen"] == "513"
        assert record["Product"].attributes["name"] == "methylmalonate-semialdehyde dehydrogenase (CoA acylating)"
        assert record["Product"].attributes["org"] == "Rhodococcus sp. PML026"
        assert record["Product"].attributes["kingdom_taxid"] == "2"
        assert record["Product"].attributes["accver"] == "WP_045840896.1"
        assert record["Product"].attributes["taxid"] == "1356405"
        assert len(record["ProteinList"]) == 2
        protein = record["ProteinList"][0]
        assert protein.tag == "Protein"
        assert protein.attributes["accver"] == "KJV04014.1"
        assert protein.attributes["kingdom"] == "Bacteria"
        assert protein.attributes["kingdom_taxid"] == "2"
        assert protein.attributes["name"] == "methylmalonic acid semialdehyde dehydrogenase mmsa"
        assert protein.attributes["org"] == "Rhodococcus sp. PML026"
        assert protein.attributes["priority"] == "0"
        assert protein.attributes["source"] == "INSDC"
        assert protein.attributes["taxid"] == "1356405"
        assert len(protein) == 1
        assert protein["CDSList"].tag == "CDSList"
        assert protein["CDSList"].attributes == {}
        assert len(protein["CDSList"]) == 2
        assert protein["CDSList"][0] == ""
        assert protein["CDSList"][0].attributes["kingdom"] == "Bacteria"
        assert protein["CDSList"][0].attributes["assembly"] == "GCA_000963615.1"
        assert protein["CDSList"][0].attributes["start"] == "264437"
        assert protein["CDSList"][0].attributes["stop"] == "265978"
        assert protein["CDSList"][0].attributes["taxid"] == "1356405"
        assert protein["CDSList"][0].attributes["strain"] == "PML026"
        assert protein["CDSList"][0].attributes["org"] == "Rhodococcus sp. PML026"
        assert protein["CDSList"][0].attributes["kingdom_taxid"] == "2"
        assert protein["CDSList"][0].attributes["accver"] == "JZIS01000004.1"
        assert protein["CDSList"][0].attributes["strand"] == "-"
        assert protein["CDSList"][1] == ""
        assert protein["CDSList"][1].attributes["kingdom"] == "Bacteria"
        assert protein["CDSList"][1].attributes["assembly"] == "GCA_000963615.1"
        assert protein["CDSList"][1].attributes["start"] == "264437"
        assert protein["CDSList"][1].attributes["stop"] == "265978"
        assert protein["CDSList"][1].attributes["taxid"] == "1356405"
        assert protein["CDSList"][1].attributes["strain"] == "PML026"
        assert protein["CDSList"][1].attributes["org"] == "Rhodococcus sp. PML026"
        assert protein["CDSList"][1].attributes["kingdom_taxid"] == "2"
        assert protein["CDSList"][1].attributes["accver"] == "KQ031368.1"
        assert protein["CDSList"][1].attributes["strand"] == "-"
        protein = record["ProteinList"][1]
        assert protein.attributes["accver"] == "WP_045840896.1"
        assert protein.attributes["source"] == "RefSeq"
        assert protein.attributes["name"] == "methylmalonate-semialdehyde dehydrogenase (CoA acylating)"
        assert protein.attributes["taxid"] == "1356405"
        assert protein.attributes["org"] == "Rhodococcus sp. PML026"
        assert protein.attributes["kingdom_taxid"] == "2"
        assert protein.attributes["kingdom"] == "Bacteria"
        assert protein.attributes["priority"] == "1"
        assert len(protein) == 1
        assert protein["CDSList"].tag == "CDSList"
        assert protein["CDSList"].attributes == {}
        assert len(protein["CDSList"]) == 1
        assert protein["CDSList"][0].attributes["assembly"] == "GCF_000963615.1"
        assert protein["CDSList"][0].attributes["start"] == "264437"
        assert protein["CDSList"][0].attributes["stop"] == "265978"
        assert protein["CDSList"][0].attributes["taxid"] == "1356405"
        assert protein["CDSList"][0].attributes["strain"] == "PML026"
        assert protein["CDSList"][0].attributes["org"] == "Rhodococcus sp. PML026"
        assert protein["CDSList"][0].attributes["kingdom_taxid"] == "2"
        assert protein["CDSList"][0].attributes["accver"] == "NZ_KQ031368.1"
        assert protein["CDSList"][0].attributes["strand"] == "-"
        assert record["Statistics"] == ""
        assert record["Statistics"].attributes["assmb_count"] == "2"
        assert record["Statistics"].attributes["nuc_count"] == "3"
        assert record["Statistics"].attributes["prot_count"] == "2"

    def test_genbank(self):
        """Test error handling when presented with GenBank non-XML data."""
        # Access the nucleotide database using efetch, but return the data
        # in GenBank format.
        # To create the GenBank file, use
        # >>> Bio.Entrez.efetch(db='nucleotide', id='NT_019265', rettype='gb')
        from Bio.Entrez import Parser

        with open("GenBank/NT_019265.gb", "rb") as stream:
            with pytest.raises(Parser.NotXMLError):
                Entrez.read(stream)
        with open("GenBank/NT_019265.gb", "rb") as stream:
            iterator = Entrez.parse(stream)
            with pytest.raises(Parser.NotXMLError):
                next(iterator)

    def test_fasta(self):
        """Test error handling when presented with Fasta non-XML data."""
        from Bio.Entrez import Parser

        with open("Fasta/wisteria.nu", "rb") as stream:
            with pytest.raises(Parser.NotXMLError):
                Entrez.read(stream)
        with open("Fasta/wisteria.nu", "rb") as stream:
            iterator = Entrez.parse(stream)
            with pytest.raises(Parser.NotXMLError):
                next(iterator)

    def test_pubmed_html(self):
        """Test error handling when presented with HTML (so XML-like) data."""
        # To create the HTML file, use
        # >>> Bio.Entrez.efetch(db="pubmed", id="19304878")
        from Bio.Entrez import Parser

        with open("Entrez/pubmed3.html", "rb") as stream:
            with pytest.raises(Parser.NotXMLError):
                Entrez.read(stream)
        # Test if the error is also raised with Entrez.parse
        with open("Entrez/pubmed3.html", "rb") as stream:
            records = Entrez.parse(stream)
            with pytest.raises(Parser.NotXMLError):
                next(records)

    def test_xml_without_declaration(self):
        """Test error handling for a missing XML declaration."""
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db="journals",id="2830,6011,7473",retmode='xml')
        from Bio.Entrez import Parser

        with open("Entrez/journals.xml", "rb") as stream:
            with pytest.raises(Parser.NotXMLError):
                Entrez.read(stream)
        # Test if the error is also raised with Entrez.parse
        with open("Entrez/journals.xml", "rb") as stream:
            records = Entrez.parse(stream)
            with pytest.raises(Parser.NotXMLError):
                next(records)

    def test_xml_without_definition(self):
        """Test error handling for a missing DTD or XML Schema."""
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db="biosample", id="3502652", rettype="xml")
        with open("Entrez/biosample.xml", "rb") as stream:
            with pytest.raises(ValueError):
                Entrez.read(stream)
        # Test if the error is also raised with Entrez.parse
        with open("Entrez/biosample.xml", "rb") as stream:
            records = Entrez.parse(stream)
            with pytest.raises(ValueError):
                next(records)

    def test_truncated_xml(self):
        """Test error handling for a truncated XML declaration."""
        from io import BytesIO

        from Bio.Entrez.Parser import CorruptedXMLError

        truncated_xml = b"""<?xml version="1.0" encoding="UTF-8"  ?>
<!DOCTYPE GBSet PUBLIC "-//NCBI//NCBI GBSeq/EN" "https://www.ncbi.nlm.nih.gov/dtd/NCBI_GBSeq.dtd">
<GBSet>
  <GBSeq>

    <GBSeq_locus>
        """
        stream = BytesIO()
        stream.write(truncated_xml)
        stream.seek(0)
        records = Entrez.parse(stream)
        with pytest.raises(CorruptedXMLError):
            next(records)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
