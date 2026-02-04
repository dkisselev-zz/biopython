# Copyright 2009-2017 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""SeqFeature related tests for SeqRecord objects from Bio.SeqIO.

Initially this takes matched tests of GenBank and FASTA files from the NCBI
and confirms they are consistent using our different parsers.
"""

import unittest
import pytest

try:
    import numpy as np
except ImportError:
    np = None  # type: ignore

from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
from Bio.SeqFeature import AfterPosition
from Bio.SeqFeature import BeforePosition
from Bio.SeqFeature import ExactPosition
from Bio.SeqFeature import OneOfPosition
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import SimpleLocation
from Bio.SeqFeature import WithinPosition
from Bio.SeqRecord import SeqRecord


class SeqRecordCreation(unittest.TestCase):
    """Test basic creation of SeqRecords."""

    def test_annotations(self):
        """Pass in annotations to SeqRecords."""
        rec = SeqRecord(Seq("ACGT"), id="Test", name="Test", description="Test")
        assert rec.annotations == {}
        rec = SeqRecord(
            Seq("ACGT"),
            id="Test",
            name="Test",
            description="Test",
            annotations={"test": ["a test"]},
        )
        assert rec.annotations["test"] == ["a test"]

    def test_letter_annotations(self):
        """Pass in letter annotations to SeqRecords."""
        rec = SeqRecord(Seq("ACGT"), id="Test", name="Test", description="Test")
        assert rec.annotations == {}
        rec = SeqRecord(
            Seq("ACGT"),
            id="Test",
            name="Test",
            description="Test",
            letter_annotations={"test": [1, 2, 3, 4]},
        )
        assert rec.letter_annotations["test"] == [1, 2, 3, 4]
        # Now try modifying it to a bad value...
        try:
            rec.letter_annotations["bad"] = "abc"
            raise AssertionError("Adding a bad letter_annotation should fail!")
        except (TypeError, ValueError) as e:
            pass
        # Now try setting it afterwards to a bad value...
        rec = SeqRecord(Seq("ACGT"), id="Test", name="Test", description="Test")
        try:
            rec.letter_annotations = {"test": [1, 2, 3]}
            raise AssertionError("Changing to bad letter_annotations should fail!")
        except (TypeError, ValueError) as e:
            pass
        # Now try setting it at creation time to a bad value...
        try:
            rec = SeqRecord(
                Seq("ACGT"),
                id="Test",
                name="Test",
                description="Test",
                letter_annotations={"test": [1, 2, 3]},
            )
            raise AssertionError("Wrong length letter_annotations should fail!")
        except (TypeError, ValueError) as e:
            pass

    def test_replacing_seq(self):
        """Replacing .seq if .letter_annotation present."""
        rec = SeqRecord(
            Seq("ACGT"),
            id="Test",
            name="Test",
            description="Test",
            letter_annotations={"example": [1, 2, 3, 4]},
        )
        try:
            rec.seq = Seq("ACGTACGT")
            raise AssertionError("Changing .seq length with letter_annotations present should fail!")
        except ValueError as e:
            assert str(e) == "You must empty the letter annotations first!"
        # Check we can replace IF the length is the same
        assert rec.seq == "ACGT"
        assert rec.letter_annotations == {"example": [1, 2, 3, 4]}
        rec.seq = Seq("NNNN")
        assert rec.seq == "NNNN"
        assert rec.letter_annotations == {"example": [1, 2, 3, 4]}

    def test_valid_id(self):
        with pytest.raises(TypeError):
            SeqRecord(Seq("ACGT"), id={})

    def test_valid_name(self):
        with pytest.raises(TypeError):
            SeqRecord(Seq("ACGT"), name={})

    def test_valid_seq(self):
        with pytest.raises(TypeError):
            SeqRecord("ACGT")

    def test_valid_description(self):
        with pytest.raises(TypeError):
            SeqRecord(Seq("ACGT"), description={})

    def test_valid_dbxrefs(self):
        with pytest.raises(TypeError):
            SeqRecord(Seq("ACGT"), dbxrefs={})

    def test_valid_annotations(self):
        with pytest.raises(TypeError):
            SeqRecord(Seq("ACGT"), annotations=[])

    def test_valid_features(self):
        with pytest.raises(TypeError):
            SeqRecord(Seq("ACGT"), features={})

    def test_default_properties(self):
        seqobj = Seq("A")
        default__dict__ = {
            "_seq": seqobj,
            "id": "<unknown id>",
            "name": "<unknown name>",
            "description": "<unknown description>",
            "dbxrefs": [],
            "annotations": {},
            "_per_letter_annotations": None,
            "features": [],
        }
        bsr = SeqRecord(seqobj)
        assert bsr.__dict__ == default__dict__


class SeqRecordMethods(unittest.TestCase):
    """Test SeqRecord methods."""

    def setUp(self):
        f0 = SeqFeature(
            SimpleLocation(0, 26),
            type="source",
            qualifiers={"mol_type": ["fake protein"]},
        )
        f1 = SeqFeature(SimpleLocation(0, ExactPosition(10)))
        f2 = SeqFeature(
            SimpleLocation(WithinPosition(12, left=12, right=15), BeforePosition(22))
        )
        f3 = SeqFeature(
            SimpleLocation(
                AfterPosition(16),
                OneOfPosition(26, [ExactPosition(25), AfterPosition(26)]),
            )
        )
        self.record = SeqRecord(
            Seq("ABCDEFGHIJKLMNOPQRSTUVWZYX"),
            id="TestID",
            name="TestName",
            description="TestDescr",
            dbxrefs=["TestXRef"],
            annotations={"k": "v"},
            letter_annotations={"fake": "X" * 26},
            features=[f0, f1, f2, f3],
        )

    def test_iter(self):
        for amino in self.record:
            assert "A" == amino
            break

    def test_contains(self):
        assert Seq("ABC") in self.record

    def test_bytes(self):
        assert b"ABCDEFGHIJKLMNOPQRSTUVWZYX" == bytes(self.record)

    def test_str(self):
        expected = """
ID: TestID
Name: TestName
Description: TestDescr
Database cross-references: TestXRef
Number of features: 4
/k=v
Per letter annotation for: fake
Seq('ABCDEFGHIJKLMNOPQRSTUVWZYX')"""
        assert expected.lstrip() == str(self.record)

    def test_repr(self):
        expected = (
            "SeqRecord(seq=Seq('ABCDEFGHIJKLMNOPQRSTUVWZYX'), "
            "id='TestID', name='TestName', description='TestDescr', dbxrefs=['TestXRef'])"
        )
        assert expected == repr(self.record)

    def test_format(self):
        expected = ">TestID TestDescr\nABCDEFGHIJKLMNOPQRSTUVWZYX\n"
        assert expected == self.record.format("fasta")

    def test_format_str(self):
        expected = ">TestID TestDescr\nABCDEFGHIJKLMNOPQRSTUVWZYX\n"
        assert expected == f"{self.record:fasta}"

    def test_format_str_binary(self):
        with pytest.raises(ValueError, match="Binary format sff cannot be used with SeqRecord format method"):
            f"{self.record:sff}"

    def test_format_spaces(self):
        rec = SeqRecord(
            Seq("ABCDEFGHIJKLMNOPQRSTUVWZYX"),
            id="TestID",
            name="TestName",
            description="TestDescr",
        )
        rec.description = "TestDescr     with5spaces"
        expected = ">TestID TestDescr     with5spaces\nABCDEFGHIJKLMNOPQRSTUVWZYX\n"
        assert expected == rec.format("fasta")

    def test_count(self):
        assert self.record.count("HIJK") == 1
        with pytest.raises(TypeError):
            SeqRecord(Seq("AC777GT")).count(7)
        with pytest.raises(TypeError):
            SeqRecord(Seq("AC777GT")).count(None)

    def test_upper(self):
        assert "ABCDEFGHIJKLMNOPQRSTUVWZYX" == self.record.lower().upper().seq
        seqobj = Seq("A")
        default__dict__ = {
            "_seq": seqobj,
            "id": "<unknown id>",
            "name": "<unknown name>",
            "description": "<unknown description>",
            "dbxrefs": [],
            "annotations": {},
            "_per_letter_annotations": None,
            "features": [],
        }
        bsr = SeqRecord(seqobj)
        bsru = bsr.upper()
        assert bsru.__dict__ == default__dict__

    def test_lower(self):
        assert "abcdefghijklmnopqrstuvwzyx" == self.record.lower().seq

    def test_isupper(self):
        assert self.record.isupper()
        assert not self.record.lower().isupper()

    def test_islower(self):
        assert not self.record.islower()
        assert self.record.lower().islower()

    def test_slicing(self):
        assert "B" == self.record[1]
        assert "BC" == self.record[1:3].seq
        with pytest.raises(ValueError):
            c = self.record["a"].seq
        if np is not None:
            start, stop = np.array([1, 3])  # numpy integers
            assert "B" == self.record[start]
            assert "BC" == self.record[start:stop].seq

    def test_slice_variants(self):
        """Simple slices using different start/end values."""
        for start in list(range(-30, 30)) + [None]:
            for end in list(range(-30, 30)) + [None]:
                if start is None and end is None:
                    continue
                rec = self.record[start:end]
                seq = self.record.seq[start:end]
                seq_str = str(self.record.seq)[start:end]
                assert seq_str == str(seq)
                assert seq_str == str(rec.seq)
                assert "X" * len(seq_str) == rec.letter_annotations["fake"]

    def test_slice_simple(self):
        """Simple slice."""
        rec = self.record
        assert len(rec) == 26
        left = rec[:10]
        assert left.seq == rec.seq[:10]
        right = rec[-10:]
        assert right.seq == rec.seq[-10:]
        mid = rec[12:22]
        assert mid.seq == rec.seq[12:22]
        for sub in [left, right, mid]:
            assert len(sub) == 10
            assert sub.id == "TestID"
            assert sub.name == "TestName"
            assert sub.description == "TestDescr"
            assert sub.letter_annotations == {"fake": "X" * 10}
            assert sub.dbxrefs == []  # May change this...
            assert sub.annotations == {}  # May change this...
            assert len(sub.features) == 1
            # By construction, each feature matches the full sliced region:
            assert sub.features[0].extract(sub.seq) == sub.seq
            assert sub.features[0].extract(sub.seq) == sub.seq

    def test_slice_zero(self):
        """Zero slice."""
        rec = self.record
        assert len(rec) == 26
        assert len(rec[2:-2]) == 22
        assert len(rec[5:2]) == 0
        assert len(rec[5:2][2:-2]) == 0

    def test_add_simple(self):
        """Simple addition."""
        rec = self.record + self.record
        assert len(rec) == 52
        assert rec.id == "TestID"
        assert rec.name == "TestName"
        assert rec.description == "TestDescr"
        assert rec.dbxrefs == ["TestXRef"]
        assert rec.annotations == {"k": "v"}
        assert rec.letter_annotations == {"fake": "X" * 52}
        assert len(rec.features) == 2 * len(self.record.features)

    def test_add_seq(self):
        """Simple addition of Seq or string."""
        for other in [Seq("BIO"), "BIO"]:
            rec = self.record + other  # will use SeqRecord's __add__ method
            assert len(rec) == 26 + 3
            assert rec.seq == str(self.record.seq) + "BIO"
            assert rec.id == "TestID"
            assert rec.name == "TestName"
            assert rec.description == "TestDescr"
            assert rec.dbxrefs == ["TestXRef"]
            assert rec.annotations == {"k": "v"}
            assert rec.letter_annotations == {}
            assert len(rec.features) == len(self.record.features)
            assert rec.features[0].type == "source"
            assert rec.features[0].location.start == 0
            assert rec.features[0].location.end == 26  # not +3

    def test_add_seqrecord(self):
        """Simple left addition of SeqRecord from genbank file."""
        other = SeqIO.read("GenBank/dbsource_wrap.gb", "gb")
        other.dbxrefs = ["dummy"]
        rec = self.record + other
        assert len(rec) == len(self.record) + len(other)
        assert rec.seq == self.record.seq + other.seq
        assert rec.id == "<unknown id>"
        assert rec.name == "<unknown name>"
        assert rec.description == "<unknown description>"
        assert rec.dbxrefs == ["TestXRef", "dummy"]
        assert len(rec.annotations) == 0
        assert len(rec.letter_annotations) == 0
        assert len(rec.features) == len(self.record.features) + len(other.features)
        assert rec.features[0].type == "source"
        assert rec.features[0].location.start == 0
        assert rec.features[0].location.end == len(self.record)  # not +3
        i = len(self.record.features)
        assert rec.features[i].type == "source"
        assert rec.features[i].location.start == len(self.record)
        assert rec.features[i].location.end == len(rec)

    def test_add_seq_left(self):
        """Simple left addition of Seq or string."""
        for other in [Seq("BIO"), "BIO"]:
            rec = other + self.record  # will use SeqRecord's __radd__ method
            assert len(rec) == 26 + 3
            assert rec.seq == "BIO" + self.record.seq
            assert rec.id == "TestID"
            assert rec.name == "TestName"
            assert rec.description == "TestDescr"
            assert rec.dbxrefs == ["TestXRef"]
            assert rec.annotations == {"k": "v"}
            assert rec.letter_annotations == {}
            assert len(rec.features) == len(self.record.features)
            assert rec.features[0].type == "source"
            assert rec.features[0].location.start == 3
            assert rec.features[0].location.end == 26 + 3

    def test_slice_add_simple(self):
        """Simple slice and add."""
        for cut in range(27):
            rec = self.record[:cut] + self.record[cut:]
            assert rec.seq == self.record.seq
            assert len(rec) == 26
            assert rec.id == "TestID"
            assert rec.name == "TestName"
            assert rec.description == "TestDescr"
            assert rec.dbxrefs == []  # May change this...
            assert rec.annotations == {}  # May change this...
            assert rec.letter_annotations == {"fake": "X" * 26}
            assert len(rec.features) <= len(self.record.features)

    def test_slice_add_shift(self):
        """Simple slice and add to shift."""
        for cut in range(27):
            rec = self.record[cut:] + self.record[:cut]
            assert rec.seq == self.record.seq[cut:] + self.record.seq[:cut]
            assert len(rec) == 26
            assert rec.id == "TestID"
            assert rec.name == "TestName"
            assert rec.description == "TestDescr"
            assert rec.dbxrefs == []  # May change this...
            assert rec.annotations == {}  # May change this...
            assert rec.letter_annotations == {"fake": "X" * 26}
            assert len(rec.features) <= len(self.record.features)


class SeqRecordMethodsMore(unittest.TestCase):
    """Test SeqRecord methods cont."""

    # This class does not have a setUp defining self.record

    def test_reverse_complement_seq(self):
        s = SeqRecord(
            Seq("ACTG"),
            id="TestID",
            name="TestName",
            description="TestDescription",
            dbxrefs=["TestDbxrefs"],
            features=[SeqFeature(SimpleLocation(0, 3), type="Site")],
            annotations={"organism": "bombyx"},
            letter_annotations={"test": "abcd"},
        )
        rc = s.reverse_complement(
            id=True,
            name=True,
            description=True,
            dbxrefs=True,
            features=True,
            annotations=True,
            letter_annotations=True,
        )

        assert "CAGT" == rc.seq
        assert "TestID" == rc.id
        assert "TestID" == s.reverse_complement(id="TestID").id

        assert "TestName" == rc.name
        assert "TestName" == s.reverse_complement(name="TestName").name

        assert "TestDescription" == rc.description
        assert "TestDescription" == s.reverse_complement(description="TestDescription").description

        assert ["TestDbxrefs"] == rc.dbxrefs
        assert ["TestDbxrefs"] == s.reverse_complement(dbxrefs=["TestDbxrefs"]).dbxrefs

        assert "[SeqFeature(SimpleLocation(ExactPosition(1), ExactPosition(4)), type='Site')]" == repr(rc.features)
        rc2 = s.reverse_complement(
            features=[SeqFeature(SimpleLocation(1, 4), type="Site")]
        )
        assert "[SeqFeature(SimpleLocation(ExactPosition(1), ExactPosition(4)), type='Site')]" == repr(rc2.features)

        assert {"organism": "bombyx"} == rc.annotations
        assert {"organism": "bombyx"} == s.reverse_complement(annotations={"organism": "bombyx"}).annotations

        assert {"test": "dcba"} == rc.letter_annotations
        assert {"test": "abcd"} == s.reverse_complement(
                letter_annotations={"test": "abcd"}
            ).letter_annotations

    def test_reverse_complement_mutable_seq(self):
        s = SeqRecord(MutableSeq("ACTG"))
        assert "CAGT" == s.reverse_complement().seq

    def test_translate(self):
        s = SeqRecord(
            Seq("ATGGTGTAA"),
            id="TestID",
            name="TestName",
            description="TestDescription",
            dbxrefs=["TestDbxrefs"],
            features=[SeqFeature(SimpleLocation(0, 3), type="Site")],
            annotations={"organism": "bombyx"},
            letter_annotations={"test": "abcdefghi"},
        )

        t = s.translate()
        assert t.seq == "MV*"
        assert t.id == "<unknown id>"
        assert t.name == "<unknown name>"
        assert t.description == "<unknown description>"
        assert not t.dbxrefs
        assert not t.features
        assert t.annotations == {"molecule_type": "protein"}
        assert not t.letter_annotations

        t = s.translate(
            cds=True,
            id=True,
            name=True,
            description=True,
            dbxrefs=True,
            annotations=True,
        )
        assert t.seq == "MV"
        assert t.id == "TestID"
        assert t.name == "TestName"
        assert t.description == "TestDescription"
        assert t.dbxrefs == ["TestDbxrefs"]
        assert not t.features
        assert t.annotations == {"organism": "bombyx", "molecule_type": "protein"}
        assert not t.letter_annotations

    def test_no_side_effects(self):
        a = SeqRecord(Seq("AAA"))
        assert a._per_letter_annotations is None
        a.reverse_complement()
        assert a._per_letter_annotations is None

        a = SeqRecord(Seq("AAA"))
        assert a._per_letter_annotations is None
        a.translate()
        assert a._per_letter_annotations is None

    def test_lt_exception(self):
        def lt():
            return SeqRecord(Seq("A")) < SeqRecord(Seq("A"))

        with pytest.raises(NotImplementedError):
            lt()

    def test_le_exception(self):
        def le():
            return SeqRecord(Seq("A")) <= SeqRecord(Seq("A"))  # type: ignore

        with pytest.raises(NotImplementedError):
            le()

    def test_eq_exception(self):
        def equality():
            return SeqRecord(Seq("A")) == SeqRecord(Seq("A"))  # type: ignore

        with pytest.raises(NotImplementedError):
            equality()

    def test_ne_exception(self):
        def notequality():
            return SeqRecord(Seq("A")) != SeqRecord(Seq("A"))  # type: ignore

        with pytest.raises(NotImplementedError):
            notequality()

    def test_gt_exception(self):
        def gt():
            return SeqRecord(Seq("A")) > SeqRecord(Seq("A"))  # type: ignore

        with pytest.raises(NotImplementedError):
            gt()

    def test_ge_exception(self):
        def ge():
            return SeqRecord(Seq("A")) >= SeqRecord(Seq("A"))  # type: ignore

        with pytest.raises(NotImplementedError):
            ge()

    def test_hash_exception(self):
        def hash1():
            hash(SeqRecord(Seq("A")))

        with pytest.raises(TypeError):
            hash1()

        def hash2():
            SeqRecord(Seq("A")).__hash__()

        with pytest.raises(TypeError):
            hash2()


class TestTranslation(unittest.TestCase):
    def setUp(self):
        self.s = SeqRecord(
            Seq("ATGGTGTAA"),
            id="TestID",
            name="TestName",
            description="TestDescription",
            dbxrefs=["TestDbxrefs"],
            features=[SeqFeature(SimpleLocation(0, 3), type="Site")],
            annotations={"organism": "bombyx"},
            letter_annotations={"test": "abcdefghi"},
        )

    def test_defaults(self):
        t = self.s.translate()
        assert t.seq == "MV*"
        assert t.id == "<unknown id>"
        assert t.name == "<unknown name>"
        assert t.description == "<unknown description>"
        assert not t.dbxrefs
        assert not t.features
        assert t.annotations == {"molecule_type": "protein"}
        assert not t.letter_annotations

    def test_preserve(self):
        t = self.s.translate(
            cds=True,
            id=True,
            name=True,
            description=True,
            dbxrefs=True,
            annotations=True,
        )
        assert t.seq == "MV"
        assert t.id == "TestID"
        assert t.name == "TestName"
        assert t.description == "TestDescription"
        assert t.dbxrefs == ["TestDbxrefs"]
        assert not t.features
        assert t.annotations == {"organism": "bombyx", "molecule_type": "protein"}
        assert not t.letter_annotations

        # Should not preserve these
        with pytest.raises(TypeError):
            self.s.translate(features=True)
        with pytest.raises(TypeError):
            self.s.translate(letter_annotations=True)

    def test_new_annot(self):
        t = self.s.translate(
            1,
            to_stop=True,
            gap="-",
            id="Foo",
            name="Bar",
            description="Baz",
            dbxrefs=["Nope"],
            features=[SeqFeature(SimpleLocation(0, 3), type="Site")],
            annotations={"a": "team"},
            letter_annotations={"aa": ["Met", "Val"]},
        )
        assert t.seq == "MV"
        assert t.id == "Foo"
        assert t.name == "Bar"
        assert t.description == "Baz"
        assert t.dbxrefs == ["Nope"]
        assert len(t.features) == 1
        assert t.annotations == {"a": "team", "molecule_type": "protein"}
        assert t.letter_annotations == {"aa": ["Met", "Val"]}


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
