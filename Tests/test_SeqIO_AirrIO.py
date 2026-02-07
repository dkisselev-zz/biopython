# Copyright 2024 by the Biopython contributors. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Tests for Bio.SeqIO.AirrIO module."""

import unittest
from io import StringIO

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.AirrIO import AirrIterator
from Bio.SeqIO.AirrIO import AirrWriter
from Bio.SeqIO.AirrIO import _format_field
from Bio.SeqIO.AirrIO import _parse_field
from Bio.SeqIO.AirrIO import _STANDARD_AIRR_FIELDS
from Bio.SeqRecord import SeqRecord


class TestParseField(unittest.TestCase):
    """Tests for the _parse_field helper."""

    def test_empty_to_none(self):
        self.assertIsNone(_parse_field("v_call", ""))

    def test_bool_true(self):
        self.assertIs(_parse_field("productive", "T"), True)

    def test_bool_false(self):
        self.assertIs(_parse_field("productive", "F"), False)

    def test_bool_invalid(self):
        self.assertRaises(ValueError, _parse_field, "productive", "true")

    def test_int(self):
        self.assertEqual(_parse_field("duplicate_count", "5"), 5)

    def test_int_invalid(self):
        self.assertRaises(ValueError, _parse_field, "duplicate_count", "abc")

    def test_float(self):
        self.assertAlmostEqual(_parse_field("v_identity", "0.987"), 0.987)

    def test_float_invalid(self):
        self.assertRaises(ValueError, _parse_field, "v_identity", "abc")

    def test_string(self):
        self.assertEqual(_parse_field("v_call", "IGHV1-2*01"), "IGHV1-2*01")

    def test_stop_codon_bool(self):
        self.assertIs(_parse_field("stop_codon", "T"), True)


class TestFormatField(unittest.TestCase):
    """Tests for the _format_field helper."""

    def test_none_to_empty(self):
        self.assertEqual(_format_field("v_call", None), "")

    def test_bool_true(self):
        self.assertEqual(_format_field("productive", True), "T")

    def test_bool_false(self):
        self.assertEqual(_format_field("productive", False), "F")

    def test_string(self):
        self.assertEqual(_format_field("v_call", "IGHV1-2*01"), "IGHV1-2*01")

    def test_int(self):
        self.assertEqual(_format_field("duplicate_count", 5), "5")

    def test_float(self):
        self.assertEqual(_format_field("v_identity", 0.987), "0.987")


class TestAirrIterator(unittest.TestCase):
    """Tests for AirrIterator."""

    def test_basic_parse(self):
        data = "sequence_id\tsequence\tv_call\tproductive\n"
        data += "s1\tACGT\tIGHV1-2*01\tT\n"
        records = list(AirrIterator(StringIO(data)))
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].id, "s1")
        self.assertEqual(str(records[0].seq), "ACGT")
        self.assertIs(records[0].annotations["productive"], True)
        self.assertEqual(records[0].annotations["v_call"], "IGHV1-2*01")

    def test_empty_file(self):
        records = list(AirrIterator(StringIO("")))
        self.assertEqual(len(records), 0)

    def test_header_only(self):
        records = list(AirrIterator(StringIO("sequence_id\tsequence\n")))
        self.assertEqual(len(records), 0)

    def test_missing_sequence_id(self):
        data = "seq_id\tsequence\ns1\tACGT\n"
        self.assertRaises(ValueError, AirrIterator, StringIO(data))

    def test_missing_sequence(self):
        data = "sequence_id\tseq\ns1\tACGT\n"
        self.assertRaises(ValueError, AirrIterator, StringIO(data))

    def test_blank_lines_skipped(self):
        data = "sequence_id\tsequence\ns1\tACGT\n\ns2\tGCTA\n"
        records = list(AirrIterator(StringIO(data)))
        self.assertEqual(len(records), 2)
        self.assertEqual(records[1].id, "s2")

    def test_bom_stripped(self):
        data = "\ufeffsequence_id\tsequence\ns1\tACGT\n"
        records = list(AirrIterator(StringIO(data)))
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].id, "s1")

    def test_file_parse(self):
        with open("AIRR/example.airr") as handle:
            records = list(AirrIterator(handle))
        self.assertEqual(len(records), 3)
        self.assertEqual(records[0].id, "seq1")
        self.assertIs(records[0].annotations["productive"], True)
        self.assertEqual(records[0].annotations["duplicate_count"], 5)
        self.assertAlmostEqual(records[0].annotations["v_identity"], 0.987)

    def test_header_only_file(self):
        with open("AIRR/header_only.airr") as handle:
            records = list(AirrIterator(handle))
        self.assertEqual(len(records), 0)

    def test_extra_fields_raise(self):
        data = "sequence_id\tsequence\ns1\tACGT\textra\n"
        self.assertRaises(ValueError, list, AirrIterator(StringIO(data)))

    def test_missing_trailing_fields_padded(self):
        # Header has 3 fields but data has only 2 -> padded with None
        data = "sequence_id\tsequence\tv_call\ns1\tACGT\n"
        records = list(AirrIterator(StringIO(data)))
        self.assertEqual(len(records), 1)
        self.assertIsNone(records[0].annotations["v_call"])

    def test_description_empty(self):
        data = "sequence_id\tsequence\ns1\tACGT\n"
        records = list(AirrIterator(StringIO(data)))
        self.assertEqual(records[0].description, "")
        self.assertEqual(records[0].name, "s1")


class TestAirrWriter(unittest.TestCase):
    """Tests for AirrWriter."""

    def test_round_trip(self):
        data = "sequence_id\tsequence\tv_call\tproductive\n"
        data += "s1\tACGT\tIGHV1-2*01\tT\n"
        data += "s2\tGCTA\tIGHV3-30*01\tF\n"
        records = list(SeqIO.parse(StringIO(data), "airr"))
        out = StringIO()
        count = SeqIO.write(records, out, "airr")
        self.assertEqual(count, 2)
        out.seek(0)
        records2 = list(SeqIO.parse(out, "airr"))
        self.assertEqual(len(records2), 2)
        self.assertEqual(records2[0].id, "s1")
        self.assertIs(records2[0].annotations["productive"], True)
        self.assertEqual(records2[1].id, "s2")
        self.assertIs(records2[1].annotations["productive"], False)

    def test_empty_write(self):
        out = StringIO()
        count = SeqIO.write([], out, "airr")
        self.assertEqual(count, 0)
        self.assertEqual(out.getvalue(), "")

    def test_buffer_includes_sparse_keys(self):
        # Buffer strategy collects the union of keys across all buffered records
        r1 = SeqRecord(
            Seq("ACGT"), id="s1", name="s1", description="",
            annotations={"v_call": "IGHV1"},
        )
        r2 = SeqRecord(
            Seq("GCTA"), id="s2", name="s2", description="",
            annotations={"v_call": "IGHV2", "extra": "x"},
        )
        out = StringIO()
        # default buffer_size=100 via SeqIO.write -> both records buffered
        count = SeqIO.write([r1, r2], out, "airr")
        self.assertEqual(count, 2)
        lines = out.getvalue().splitlines()
        # Header is union of keys: v_call from r1, extra added from r2
        self.assertEqual(lines[0], "sequence_id\tsequence\tv_call\textra")
        # r1 gets an empty extra column
        self.assertEqual(lines[1], "s1\tACGT\tIGHV1\t")
        # r2 has both values
        self.assertEqual(lines[2], "s2\tGCTA\tIGHV2\tx")

    def test_post_buffer_keys_dropped(self):
        # Keys seen only after the buffer window are silently dropped
        r1 = SeqRecord(
            Seq("ACGT"), id="s1", name="s1", description="",
            annotations={"v_call": "IGHV1"},
        )
        r2 = SeqRecord(
            Seq("GCTA"), id="s2", name="s2", description="",
            annotations={"v_call": "IGHV2", "late_key": "x"},
        )
        out = StringIO()
        writer = AirrWriter(out, buffer_size=1)
        count = writer.write_records(iter([r1, r2]))
        self.assertEqual(count, 2)
        lines = out.getvalue().splitlines()
        # late_key was not seen during the single-record buffer window
        self.assertEqual(lines[0], "sequence_id\tsequence\tv_call")
        self.assertEqual(lines[2], "s2\tGCTA\tIGHV2")

    def test_buffer_larger_than_dataset(self):
        # buffer_size bigger than total records still works correctly
        records = [
            SeqRecord(
                Seq("ACGT"), id=f"s{i}", name=f"s{i}", description="",
                annotations={"v_call": f"IGHV{i}", "j_call": f"IGHJ{i}"},
            )
            for i in range(3)
        ]
        out = StringIO()
        writer = AirrWriter(out, buffer_size=100)
        count = writer.write_records(iter(records))
        self.assertEqual(count, 3)
        lines = out.getvalue().splitlines()
        self.assertEqual(len(lines), 4)  # 1 header + 3 data
        self.assertEqual(lines[0], "sequence_id\tsequence\tv_call\tj_call")

    def test_strict_streaming_standard_fields(self):
        # buffer_size=0 uses _STANDARD_AIRR_FIELDS as the header
        r = SeqRecord(
            Seq("ACGT"), id="s1", name="s1", description="",
            annotations={"v_call": "IGHV1"},
        )
        out = StringIO()
        writer = AirrWriter(out, buffer_size=0)
        count = writer.write_records(iter([r]))
        self.assertEqual(count, 1)
        cols = out.getvalue().splitlines()[0].split("\t")
        self.assertEqual(cols[0], "sequence_id")
        self.assertEqual(cols[1], "sequence")
        for field in _STANDARD_AIRR_FIELDS:
            self.assertIn(field, cols)

    def test_strict_streaming_empty(self):
        # buffer_size=0 with no records still writes nothing
        out = StringIO()
        writer = AirrWriter(out, buffer_size=0)
        count = writer.write_records(iter([]))
        self.assertEqual(count, 0)
        self.assertEqual(out.getvalue(), "")

    def test_explicit_fields_override_buffer(self):
        # Explicit fields take priority over buffer inference
        r1 = SeqRecord(
            Seq("ACGT"), id="s1", name="s1", description="",
            annotations={"v_call": "IGHV1", "extra": "y"},
        )
        out = StringIO()
        writer = AirrWriter(out, fields=["sequence_id", "sequence", "v_call"])
        count = writer.write_records(iter([r1]))
        self.assertEqual(count, 1)
        lines = out.getvalue().splitlines()
        self.assertEqual(lines[0], "sequence_id\tsequence\tv_call")
        # extra dropped because not in explicit schema
        self.assertEqual(lines[1], "s1\tACGT\tIGHV1")

    def test_none_annotation_writes_empty(self):
        r = SeqRecord(
            Seq("ACGT"), id="s1", name="s1", description="",
            annotations={"v_call": None}
        )
        out = StringIO()
        count = SeqIO.write([r], out, "airr")
        self.assertEqual(count, 1)
        lines = out.getvalue().splitlines()
        self.assertEqual(lines[0], "sequence_id\tsequence\tv_call")
        self.assertEqual(lines[1], "s1\tACGT\t")


class TestAirrIndex(unittest.TestCase):
    """Tests for SeqIO.index with AIRR format."""

    def test_index_keys(self):
        d = SeqIO.index("AIRR/example.airr", "airr")
        self.assertEqual(set(d.keys()), {"seq1", "seq2", "seq3"})
        d.close()

    def test_index_retrieve(self):
        d = SeqIO.index("AIRR/example.airr", "airr")
        record = d["seq1"]
        self.assertEqual(record.id, "seq1")
        self.assertIs(record.annotations["productive"], True)
        self.assertEqual(record.annotations["duplicate_count"], 5)
        d.close()

    def test_get_raw(self):
        d = SeqIO.index("AIRR/example.airr", "airr")
        raw = d.get_raw("seq1")
        self.assertIn(b"sequence_id", raw)
        self.assertIn(b"seq1", raw)
        self.assertIn(b"sequence", raw)
        d.close()


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
