# Copyright 2010-2013 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Testing Bio.TogoWS online code."""

import unittest
import pytest
from io import StringIO
from urllib.error import HTTPError


from Bio import Medline

# In order to check any sequences returned
from Bio import SeqIO

# We want to test these:
from Bio import TogoWS
from Bio.SeqUtils.CheckSum import seguid

pytestmark = pytest.mark.online


#####################################################################


class TogoFields(unittest.TestCase):
    def test_invalid_database(self):
        """Check asking for fields of invalid database fails."""
        with pytest.raises(IOError):
            TogoWS._get_fields("http://togows.dbcls.jp/entry/invalid?fields")

    def test_databases(self):
        """Check supported databases."""
        dbs = set(TogoWS._get_entry_dbs())
        expected = {
            "nuccore",
            "nucest",
            "nucgss",
            "nucleotide",
            "protein",
            "gene",
            "homologene",
            "snp",
            "mesh",
            "pubmed",  # 'embl',
            "uniprot",
            "uniparc",
            "uniref100",
            "uniref90",
            "uniref50",
            "ddbj",
            "dad",
            "pdb",
            "compound",
            "drug",
            "enzyme",
            "genes",
            "glycan",
            "orthology",
            "reaction",
            "module",
            "pathway",
        }
        assert dbs.issuperset(expected), f"Missing DB: {', '.join(sorted(expected.difference(dbs)))}"

    def test_pubmed(self):
        """Check supported fields for pubmed database."""
        fields = set(TogoWS._get_entry_fields("pubmed"))
        assert fields.issuperset(
                ["abstract", "au", "authors", "doi", "mesh", "so", "title"]
            ), fields

    def test_ncbi_protein(self):
        """Check supported fields for NCBI protein database."""
        fields = set(TogoWS._get_entry_fields("ncbi-protein"))
        assert fields.issuperset(
                [
                    "entry_id",
                    "length",
                    "strand",
                    "moltype",
                    "linearity",
                    "division",
                    "date",
                    "definition",
                    "accession",
                    "accessions",
                    "version",
                    "versions",
                    "acc_version",
                    "gi",
                    "keywords",
                    "organism",
                    "common_name",
                    "taxonomy",
                    "comment",
                    "seq",
                ]
            ), fields

    def test_ddbj(self):
        """Check supported fields for ddbj database."""
        fields = set(TogoWS._get_entry_fields("ddbj"))
        assert fields.issuperset(
                [
                    "entry_id",
                    "length",
                    "strand",
                    "moltype",
                    "linearity",
                    "division",
                    "date",
                    "definition",
                    "accession",
                    "accessions",
                    "version",
                    "versions",
                    "acc_version",
                    "gi",
                    "keywords",
                    "organism",
                    "common_name",
                    "taxonomy",
                    "comment",
                    "seq",
                ]
            ), fields

    def test_uniprot(self):
        """Check supported fields for uniprot database."""
        fields = set(TogoWS._get_entry_fields("uniprot"))
        assert fields.issuperset(["definition", "entry_id", "seq"]), fields

    def test_pdb(self):
        """Check supported fields for pdb database."""
        fields = set(TogoWS._get_entry_fields("pdb"))
        assert fields.issuperset(["accession", "chains", "keywords", "models"]), fields


class TogoEntry(unittest.TestCase):
    def test_pubmed_16381885(self):
        """Bio.TogoWS.entry("pubmed", "16381885")."""
        # Gives Medline plain text
        handle = TogoWS.entry("pubmed", "16381885")
        data = Medline.read(handle)
        handle.close()
        assert data["TI"] == "From genomics to chemical genomics: new developments in KEGG."
        assert data["AU"] == [
                "Kanehisa M",
                "Goto S",
                "Hattori M",
                "Aoki-Kinoshita KF",
                "Itoh M",
                "Kawashima S",
                "Katayama T",
                "Araki M",
                "Hirakawa M",
            ]

    def test_pubmed_16381885_ti(self):
        """Bio.TogoWS.entry("pubmed", "16381885", field="title")."""
        handle = TogoWS.entry("pubmed", "16381885", field="title")
        data = handle.read().strip()
        handle.close()
        assert data == "From genomics to chemical genomics: new developments in KEGG."

    def test_pubmed_16381885_title(self):
        """Bio.TogoWS.entry("pubmed", "16381885", field="title")."""
        handle = TogoWS.entry("pubmed", "16381885", field="title")
        data = handle.read().strip()
        handle.close()
        assert data == "From genomics to chemical genomics: new developments in KEGG."

    def test_pubmed_16381885_au(self):
        """Bio.TogoWS.entry("pubmed", "16381885", field="au")."""
        # Gives one name per line (i.e. \n separated), no dots
        handle = TogoWS.entry("pubmed", "16381885", field="au")
        data = handle.read().strip().split("\n")
        handle.close()
        assert data == [
                "Kanehisa M",
                "Goto S",
                "Hattori M",
                "Aoki-Kinoshita KF",
                "Itoh M",
                "Kawashima S",
                "Katayama T",
                "Araki M",
                "Hirakawa M",
            ]

    def test_pubmed_16381885_authors(self):
        """Bio.TogoWS.entry("pubmed", "16381885", field="authors")."""
        # Gives names tab separated (i.e. \t separated)
        handle = TogoWS.entry("pubmed", "16381885", field="authors")
        data = handle.read().strip().split("\t")
        handle.close()
        assert data == [
                "Kanehisa, M.",
                "Goto, S.",
                "Hattori, M.",
                "Aoki-Kinoshita, K. F.",
                "Itoh, M.",
                "Kawashima, S.",
                "Katayama, T.",
                "Araki, M.",
                "Hirakawa, M.",
            ]

    def test_pubmed_16381885_invalid_field(self):
        """Bio.TogoWS.entry("pubmed", "16381885", field="invalid_for_testing")."""
        with pytest.raises(ValueError):
            TogoWS.entry("pubmed", "16381885", field="invalid_for_testing")

    def test_pubmed_16381885_invalid_format(self):
        """Bio.TogoWS.entry("pubmed", "16381885", format="invalid_for_testing")."""
        with pytest.raises(ValueError):
            TogoWS.entry("pubmed", "16381885", format="invalid_for_testing")

    def test_pubmed_invalid_id(self):
        """Bio.TogoWS.entry("pubmed", "invalid_for_testing")."""
        with pytest.raises(IOError):
            TogoWS.entry("pubmed", "invalid_for_testing")

    def test_pubmed_16381885_and_19850725(self):
        """Bio.TogoWS.entry("pubmed", "16381885,19850725")."""
        handle = TogoWS.entry("pubmed", "16381885,19850725")
        records = list(Medline.parse(handle))
        handle.close()
        assert len(records) == 2
        assert records[0]["TI"] == "From genomics to chemical genomics: new developments in KEGG."
        assert records[0]["AU"] == [
                "Kanehisa M",
                "Goto S",
                "Hattori M",
                "Aoki-Kinoshita KF",
                "Itoh M",
                "Kawashima S",
                "Katayama T",
                "Araki M",
                "Hirakawa M",
            ]
        assert (records[1]["TI"] == "DDBJ launches a new archive database with "
            "analytical tools for next-generation sequence data.")
        assert records[1]["AU"] == [
                "Kaminuma E",
                "Mashima J",
                "Kodama Y",
                "Gojobori T",
                "Ogasawara O",
                "Okubo K",
                "Takagi T",
                "Nakamura Y",
            ]

    def test_pubmed_16381885_and_19850725_authors(self):
        """Bio.TogoWS.entry("pubmed", "16381885,19850725", field="authors")."""
        handle = TogoWS.entry("pubmed", "16381885,19850725", field="authors")
        # Little hack to remove blank lines...
        # names = handle.read().replace("\n\n", "\n").strip().split("\n")
        names = handle.read().strip().split("\n")
        handle.close()
        assert 2 == len(names)
        names1, names2 = names
        assert names1.split("\t") == [
                "Kanehisa, M.",
                "Goto, S.",
                "Hattori, M.",
                "Aoki-Kinoshita, K. F.",
                "Itoh, M.",
                "Kawashima, S.",
                "Katayama, T.",
                "Araki, M.",
                "Hirakawa, M.",
            ]
        assert names2.split("\t") == [
                "Kaminuma, E.",
                "Mashima, J.",
                "Kodama, Y.",
                "Gojobori, T.",
                "Ogasawara, O.",
                "Okubo, K.",
                "Takagi, T.",
                "Nakamura, Y.",
            ]

    def test_invalid_db(self):
        """Bio.TogoWS.entry("invalid_db", "invalid_id")."""
        with pytest.raises(ValueError):
            TogoWS.entry("invalid_db", "invalid_id")

    def test_ddbj_genbank_length(self):
        """Bio.TogoWS.entry("ddbj", "X52960", field="length")."""
        handle = TogoWS.entry("ddbj", "X52960", field="length")
        data = handle.read().strip()  # ignore trailing \n
        handle.close()
        assert data == "248"

    def test_ddbj_genbank(self):
        """Bio.TogoWS.entry("ddbj", "X52960")."""
        handle = TogoWS.entry("ddbj", "X52960")  # Returns "genbank" format
        record = SeqIO.read(handle, "gb")
        handle.close()
        assert record.id == "X52960.1"
        assert record.name == "X52960"
        assert len(record) == 248
        assert seguid(record.seq) == "Ktxz0HgMlhQmrKTuZpOxPZJ6zGU"

    def test_nucleotide_genbank_length(self):
        """Bio.TogoWS.entry("nucleotide", "X52960", field="length")."""
        handle = TogoWS.entry("nucleotide", "X52960", field="length")
        data = handle.read().strip()  # ignore trailing \n
        handle.close()
        assert data == "248"

    def test_nucleotide_genbank_seq(self):
        """Bio.TogoWS.entry("nucleotide", "X52960", field="seq")."""
        handle = TogoWS.entry("nucleotide", "X52960", field="seq")
        data = handle.read().strip()  # ignore trailing \n
        handle.close()
        assert seguid(data) == "Ktxz0HgMlhQmrKTuZpOxPZJ6zGU"

    def test_nucleotide_genbank_definition(self):
        """Bio.TogoWS.entry("nucleotide", "X52960", field="definition")."""
        handle = TogoWS.entry("nucleotide", "X52960", field="definition")
        data = handle.read().strip()  # ignore trailing \n
        handle.close()
        assert data == "Coleus blumei viroid 1 (CbVd) RNA."

    def test_nucleotide_genbank_accession(self):
        """Bio.TogoWS.entry("nucleotide", "X52960", field="accession")."""
        handle = TogoWS.entry("nucleotide", "X52960", field="accession")
        data = handle.read().strip()  # ignore trailing \n
        handle.close()
        assert data == "X52960"

    def test_nucleotide_genbank_version(self):
        """Bio.TogoWS.entry("nucleotide", "X52960", field="version")."""
        handle = TogoWS.entry("nucleotide", "X52960", field="version")
        data = handle.read().strip()  # ignore trailing \n
        handle.close()
        assert data == "1"

    def test_nucleotide_genbank_acc_version(self):
        """Bio.TogoWS.entry("nucleotide", "X52960", field="acc_version")."""
        handle = TogoWS.entry("nucleotide", "X52960", field="acc_version")
        data = handle.read().strip()  # ignore trailing \n
        handle.close()
        assert data == "X52960.1"

    def test_nucleotide_genbank_organism(self):
        """Bio.TogoWS.entry("nucleotide", "X52960", field="organism")."""
        handle = TogoWS.entry("nucleotide", "X52960", field="organism")
        data = handle.read().strip()  # ignore trailing \n
        handle.close()
        assert data == "Coleus blumei viroid 1"

    def test_ddbj_genbank_invalid_field(self):
        """Bio.TogoWS.entry("nucleotide", "X52960", field="invalid_for_testing")."""
        with pytest.raises(ValueError):
            TogoWS.entry("nucleotide", "X52960", field="invalid_for_testing")

    def test_nucleotide_invalid_format(self):
        """Bio.TogoWS.entry("nucleotide", "X52960", format="invalid_for_testing")."""
        with pytest.raises(ValueError):
            TogoWS.entry("nucleotide", "X52960", format="invalid_for_testing")

    def test_ddbj_gff3(self):
        """Bio.TogoWS.entry("ddbj", "X52960", format="gff")."""
        handle = TogoWS.entry("ddbj", "X52960", format="gff")
        data = handle.read()
        handle.close()
        assert data.startswith("##gff-version 3\nX52960\tDDBJ\t"), data

    def test_genbank_gff3(self):
        """Bio.TogoWS.entry("nucleotide", "X52960", format="gff")."""
        # Note - Using manual URL with genbank instead of nucleotide works
        handle = TogoWS.entry("nucleotide", "X52960", format="gff")
        data = handle.read()
        handle.close()
        assert data.startswith("##gff-version 3\nX52960\tGenbank\t"), data

    def test_ddbj_fasta(self):
        """Bio.TogoWS.entry("ddbj", "X52960", "fasta")."""
        handle = TogoWS.entry("ddbj", "X52960", "fasta")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        assert "X52960" in record.id
        assert "X52960" in record.name
        assert len(record) == 248
        assert seguid(record.seq) == "Ktxz0HgMlhQmrKTuZpOxPZJ6zGU"

    def test_nucleotide_fasta(self):
        """Bio.TogoWS.entry("nucleotide", "6273291", "fasta")."""
        handle = TogoWS.entry("nucleotide", "6273291", "fasta")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        # NCBI is phasing out GI numbers, so no longer true:
        # self.assertIn("6273291", record.id)
        # self.assertIn("6273291", record.name)
        assert "AF191665.1" in record.id
        assert "AF191665.1" in record.name
        assert len(record) == 902
        assert seguid(record.seq) == "bLhlq4mEFJOoS9PieOx4nhGnjAQ"

    def test_protein_fasta(self):
        """Bio.TogoWS.entry("protein", "16130152", "fasta")."""
        handle = TogoWS.entry("protein", "16130152", "fasta")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        # NCBI is phasing out GI numbers, so no longer true:
        # self.assertIn("16130152", record.id)
        # self.assertIn("16130152", record.name)
        assert "NP_416719.1" in record.id
        assert "NP_416719.1" in record.name
        assert " porin " in record.description
        assert len(record) == 367
        assert seguid(record.seq) == "fCjcjMFeGIrilHAn6h+yju267lg"


class TogoSearch(unittest.TestCase):
    """Search tests."""

    def test_bad_args_just_limit(self):
        """Reject Bio.TogoWS.search(...) with just limit."""
        with pytest.raises(ValueError):
            TogoWS.search("pubmed", "lung+cancer", limit=10)

    def test_bad_args_just_offset(self):
        """Reject Bio.TogoWS.search(...) with just offset."""
        with pytest.raises(ValueError):
            TogoWS.search("pubmed", "lung+cancer", offset=10)

    def test_bad_args_zero_limit(self):
        """Reject Bio.TogoWS.search(...) with zero limit."""
        with pytest.raises(ValueError):
            TogoWS.search("pubmed", "lung+cancer", offset=1, limit=0)

    def test_bad_args_zero_offset(self):
        """Reject Bio.TogoWS.search(...) with zero offset."""
        with pytest.raises(ValueError):
            TogoWS.search("pubmed", "lung+cancer", offset=0, limit=10)

    def test_bad_args_non_int_offset(self):
        """Reject Bio.TogoWS.search(...) with non-integer offset."""
        with pytest.raises(ValueError):
            TogoWS.search("pubmed", "lung+cancer", offset="test", limit=10)

    def test_bad_args_non_int_limit(self):
        """Reject Bio.TogoWS.search(...) with non-integer limit."""
        with pytest.raises(ValueError):
            TogoWS.search("pubmed", "lung+cancer", offset=1, limit="lots")

    def test_pubmed_search_togows(self):
        """Bio.TogoWS.search_iter("pubmed", "TogoWS") etc."""
        self.check("pubmed", "TogoWS", ["20472643"])

    def test_pubmed_search_bioruby(self):
        """Bio.TogoWS.search_iter("pubmed", "BioRuby") etc."""
        self.check(
            "pubmed",
            "BioRuby",
            ["22994508", "22399473", "20739307", "20015970", "14693808"],
        )

    def test_pubmed_search_porin(self):
        """Bio.TogoWS.search_iter("pubmed", "human porin") etc.

        Count was 357 at time of writing, this was chosen to
        be larger than the default chunk size for iteration,
        but still not too big to download the full list.
        """
        self.check("pubmed", "human porin", ["21189321", "21835183"])

    # TogoWS search for PDBj currently unavailable
    #    def test_pdb_search_porin(self):
    #        """Bio.TogoWS.search_iter("pdb", "porin") etc
    #
    #        Count was about 161 at time of writing.
    #        """
    #        self.check("pdb", "porin", ["2j1n", "2vqg", "3m8b", "2k0l"])

    def test_uniprot_search_lung_cancer(self):
        """Bio.TogoWS.search_iter("uniprot", "terminal+lung+cancer", limit=150) etc.

        Search count was 211 at time of writing, a bit large to
        download all the results in a unit test. Want to use a limit
        larger than the batch size (100) to ensure at least two
        batches.
        """
        self.check("uniprot", "terminal+lung+cancer", limit=150)

    def check(self, database, search_term, expected_matches=(), limit=None):
        if expected_matches and limit:
            raise ValueError("Bad test - TogoWS makes no promises about order")
        try:
            search_count = TogoWS.search_count(database, search_term)
        except HTTPError as err:
            raise ValueError(f"{err} from {err.url}") from None
        if expected_matches:
            assert search_count >= len(expected_matches)
        if search_count > 5000 and not limit:
            print("%i results, skipping" % search_count)
            return
        if limit:
            count = min(search_count, limit)
        else:
            count = search_count

        # Iteration should find everything... unless a limit is used
        search_iter = list(TogoWS.search_iter(database, search_term, limit))
        assert count == len(search_iter)
        for match in expected_matches:
            assert match in search_iter, f"Expected {match} in results"


class TogoConvert(unittest.TestCase):
    """Conversion tests."""

    def test_invalid_format(self):
        """Check convert file format checking."""
        with pytest.raises(ValueError):
            TogoWS.convert(StringIO("PLACEHOLDER"), "genbank", "invalid_for_testing")
        with pytest.raises(ValueError):
            TogoWS.convert(StringIO("PLACEHOLDER"), "invalid_for_testing", "fasta")

    def test_genbank_to_fasta(self):
        """Conversion of GenBank to FASTA."""
        filename = "GenBank/NC_005816.gb"
        old = SeqIO.read(filename, "gb")
        with open(filename) as handle:
            new = SeqIO.read(TogoWS.convert(handle, "genbank", "fasta"), "fasta")
        assert old.seq == new.seq


#    def test_genbank_to_embl(self):
#        """Conversion of GenBank to EMBL."""
#        filename = "GenBank/NC_005816.gb"
#        old = SeqIO.read(filename, "gb")
#        with open(filename) as handle:
#            new = SeqIO.read(TogoWS.convert(handle, "genbank", "embl"), "embl")
#        self.assertEqual(str(old.seq), str(new.seq))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
