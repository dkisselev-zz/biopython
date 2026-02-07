# Copyright 2024 by Biopython Contributors. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Tests for Bio.SeqIO.AirrIO (AIRR format I/O)."""

import os
import unittest
from io import StringIO

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class TestAirrReading(unittest.TestCase):
    """Tests for reading AIRR files."""

    def test_read_minimal(self):
        """Read minimal AIRR file with required fields only."""
        records = list(SeqIO.parse("Airr/minimal.tsv", "airr"))

        # Check we got 3 records
        self.assertEqual(len(records), 3)

        # Check first record
        rec = records[0]
        self.assertEqual(rec.id, "SEQ001")
        self.assertEqual(rec.name, "SEQ001")
        self.assertEqual(str(rec.seq), "ATGCGTATCGATCGCGATACGATTAGGCGGATGTAA")

        # Check AIRR annotations
        self.assertIn("airr", rec.annotations)
        airr = rec.annotations["airr"]
        self.assertEqual(airr["rev_comp"], False)
        self.assertEqual(airr["productive"], True)

        # Check second record
        rec = records[1]
        self.assertEqual(rec.id, "SEQ002")
        self.assertEqual(airr["productive"], True)

        # Check third record (non-productive)
        rec = records[2]
        self.assertEqual(rec.id, "SEQ003")
        self.assertEqual(rec.annotations["airr"]["productive"], False)

    def test_read_standard(self):
        """Read standard AIRR file with full fields."""
        records = list(SeqIO.parse("Airr/standard.tsv", "airr"))

        # Check we got 10 records
        self.assertEqual(len(records), 10)

        # Check first record has comprehensive fields
        rec = records[0]
        self.assertEqual(rec.id, "SEQ001")
        airr = rec.annotations["airr"]

        # Check gene calls
        self.assertEqual(airr["v_call"], "IGHV1-69*01")
        self.assertEqual(airr["d_call"], "IGHD3-10*01")
        self.assertEqual(airr["j_call"], "IGHJ4*02")

        # Check sequences
        self.assertEqual(
            airr["junction"],
            "TGTGCGAGAGATAGGCTCGGGACTAGCTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG",
        )
        self.assertEqual(airr["junction_aa"], "CARDRLGTSWGQGTLVTVSS")
        self.assertEqual(airr["cdr3_aa"], "CARDRLGTS")

        # Check integer fields
        self.assertEqual(airr["v_sequence_start"], 1)
        self.assertEqual(airr["v_sequence_end"], 87)
        self.assertEqual(airr["junction_length"], 65)

        # Check float fields
        self.assertAlmostEqual(airr["v_identity"], 0.97)
        self.assertAlmostEqual(airr["v_score"], 487.5)

        # Check boolean fields
        self.assertEqual(airr["rev_comp"], False)
        self.assertEqual(airr["productive"], True)
        self.assertEqual(airr["vj_in_frame"], True)

    def test_read_custom_columns(self):
        """Read AIRR file with custom non-standard columns."""
        records = list(SeqIO.parse("Airr/custom_columns.tsv", "airr"))

        # Check we got 5 records
        self.assertEqual(len(records), 5)

        # Check first record has custom fields
        rec = records[0]
        airr = rec.annotations["airr"]

        # Check standard fields
        self.assertEqual(airr["v_call"], "IGHV1-69*01")
        self.assertEqual(airr["j_call"], "IGHJ4*02")

        # Check custom fields are preserved as strings
        self.assertEqual(airr["custom_score"], "0.95")
        self.assertEqual(airr["custom_flag"], "pass")
        self.assertEqual(airr["sample_id"], "SAMPLE_A")

    def test_read_missing_fields(self):
        """Test handling of empty/missing field values."""
        records = list(SeqIO.parse("Airr/standard.tsv", "airr"))

        # SEQ008 has empty v_score
        rec = records[7]
        self.assertEqual(rec.id, "SEQ008")
        self.assertIsNone(rec.annotations["airr"]["v_score"])

        # SEQ009 has empty d_call
        rec = records[8]
        self.assertEqual(rec.id, "SEQ009")
        self.assertIsNone(rec.annotations["airr"]["d_call"])

    def test_field_order_preserved(self):
        """Test that field order is preserved for round-trip."""
        records = list(SeqIO.parse("Airr/standard.tsv", "airr"))
        rec = records[0]

        # Check field order is stored
        self.assertIn("_airr_field_order", rec.annotations)
        field_order = rec.annotations["_airr_field_order"]

        # Should start with sequence_id and sequence
        self.assertEqual(field_order[0], "sequence_id")
        self.assertEqual(field_order[1], "sequence")

    def test_malformed_no_header(self):
        """Test error handling for file without header."""
        data = "SEQ001\tATGCGT\n"
        handle = StringIO(data)

        with self.assertRaises(ValueError) as cm:
            list(SeqIO.parse(handle, "airr"))
        self.assertIn("sequence_id", str(cm.exception))

    def test_malformed_missing_sequence_id(self):
        """Test error handling for missing sequence_id column."""
        data = "sequence\trev_comp\nATGCGT\tF\n"
        handle = StringIO(data)

        with self.assertRaises(ValueError) as cm:
            list(SeqIO.parse(handle, "airr"))
        self.assertIn("sequence_id", str(cm.exception))

    def test_malformed_missing_sequence(self):
        """Test error handling for missing sequence column."""
        data = "sequence_id\trev_comp\nSEQ001\tF\n"
        handle = StringIO(data)

        with self.assertRaises(ValueError) as cm:
            list(SeqIO.parse(handle, "airr"))
        self.assertIn("sequence", str(cm.exception))

    def test_truncated_line_missing_required_fields(self):
        """Test error when truncated line is missing required fields."""
        # Header has 5 columns, but data line only has 1 (missing sequence)
        data = "sequence_id\tsequence\trev_comp\tproductive\tv_call\nSEQ001\n"
        handle = StringIO(data)

        with self.assertRaises(ValueError) as cm:
            list(SeqIO.parse(handle, "airr"))
        self.assertIn("Truncated line", str(cm.exception))
        self.assertIn("required fields", str(cm.exception))

    def test_truncated_line_significantly_missing(self):
        """Test error when line is missing >20% of columns."""
        # Header has 10 columns, data line has only 7 (30% missing)
        data = (
            "sequence_id\tsequence\trev_comp\tproductive\tv_call\tj_call\t"
            "junction\tjunction_aa\tcdr3\tcdr3_aa\n"
            "SEQ001\tATGCGT\tF\tT\tIGHV1-69\tIGHJ4\tTGTGCGAGA\n"
        )
        handle = StringIO(data)

        with self.assertRaises(ValueError) as cm:
            list(SeqIO.parse(handle, "airr"))
        self.assertIn("Truncated line", str(cm.exception))
        self.assertIn("30", str(cm.exception))  # 30% missing

    def test_mildly_truncated_line_warns(self):
        """Test warning when line is missing <20% of columns."""
        import warnings

        # Header has 10 columns, data line has 9 (10% missing)
        data = (
            "sequence_id\tsequence\trev_comp\tproductive\tv_call\tj_call\t"
            "junction\tjunction_aa\tcdr3\tcdr3_aa\n"
            "SEQ001\tATGCGT\tF\tT\tIGHV1-69\tIGHJ4\tTGTGCGAGA\tCARDRL\tTGCGAGA\n"
        )
        handle = StringIO(data)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            records = list(SeqIO.parse(handle, "airr"))

            # Should have parsed successfully
            self.assertEqual(len(records), 1)

            # Should have issued a warning
            self.assertEqual(len(w), 1)
            self.assertIn("Missing 1 trailing columns", str(w[0].message))

    def test_invalid_type_conversion_warns(self):
        """Test warnings for invalid type conversions."""
        import warnings

        # Create data with invalid types
        data = (
            "sequence_id\tsequence\tproductive\tv_sequence_start\tv_identity\n"
            "SEQ001\tATGCGT\tYES\tabc\tinvalid\n"  # Invalid bool, int, float
        )
        handle = StringIO(data)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            records = list(SeqIO.parse(handle, "airr"))

            # Should have parsed successfully
            self.assertEqual(len(records), 1)

            # Should have issued 3 warnings (one for each invalid conversion)
            self.assertEqual(len(w), 3)

            # Check warning messages
            warning_messages = [str(warning.message) for warning in w]
            self.assertTrue(
                any("Invalid boolean value" in msg for msg in warning_messages)
            )
            self.assertTrue(
                any(
                    "Cannot convert 'abc' to integer" in msg for msg in warning_messages
                )
            )
            self.assertTrue(
                any(
                    "Cannot convert 'invalid' to float" in msg
                    for msg in warning_messages
                )
            )

            # Values should be None after failed conversion
            airr = records[0].annotations["airr"]
            self.assertIsNone(airr["productive"])
            self.assertIsNone(airr["v_sequence_start"])
            self.assertIsNone(airr["v_identity"])


class TestAirrWriting(unittest.TestCase):
    """Tests for writing AIRR files."""

    def test_write_simple(self):
        """Write simple AIRR records."""
        # Create test records
        rec1 = SeqRecord(Seq("ATGCGT"), id="SEQ001", name="SEQ001", description="")
        rec1.annotations["airr"] = {
            "rev_comp": False,
            "productive": True,
            "v_call": "IGHV1-69*01",
        }
        rec1.annotations["_airr_field_order"] = [
            "sequence_id",
            "sequence",
            "rev_comp",
            "productive",
            "v_call",
        ]

        rec2 = SeqRecord(Seq("GCATGC"), id="SEQ002", name="SEQ002", description="")
        rec2.annotations["airr"] = {
            "rev_comp": True,
            "productive": False,
            "v_call": "IGHV1-2*02",
        }
        rec2.annotations["_airr_field_order"] = [
            "sequence_id",
            "sequence",
            "rev_comp",
            "productive",
            "v_call",
        ]

        # Write to string
        output = StringIO()
        count = SeqIO.write([rec1, rec2], output, "airr")
        self.assertEqual(count, 2)

        # Check output
        output.seek(0)
        lines = output.readlines()
        self.assertEqual(len(lines), 3)  # Header + 2 records

        # Check header
        self.assertEqual(
            lines[0].strip(), "sequence_id\tsequence\trev_comp\tproductive\tv_call"
        )

        # Check first record
        fields = lines[1].strip().split("\t")
        self.assertEqual(fields[0], "SEQ001")
        self.assertEqual(fields[1], "ATGCGT")
        self.assertEqual(fields[2], "F")
        self.assertEqual(fields[3], "T")
        self.assertEqual(fields[4], "IGHV1-69*01")

        # Check second record
        fields = lines[2].strip().split("\t")
        self.assertEqual(fields[0], "SEQ002")
        self.assertEqual(fields[2], "T")
        self.assertEqual(fields[3], "F")

    def test_write_none_values(self):
        """Test that None values are written as empty strings."""
        rec = SeqRecord(Seq("ATGCGT"), id="SEQ001", name="SEQ001", description="")
        rec.annotations["airr"] = {
            "v_call": "IGHV1-69*01",
            "d_call": None,
            "j_call": "IGHJ4*02",
            "v_score": None,
        }
        rec.annotations["_airr_field_order"] = [
            "sequence_id",
            "sequence",
            "v_call",
            "d_call",
            "j_call",
            "v_score",
        ]

        output = StringIO()
        SeqIO.write([rec], output, "airr")

        output.seek(0)
        lines = output.readlines()
        # Don't strip the line to preserve trailing tabs
        fields = lines[1].rstrip("\n\r").split("\t")

        # Check that None values are empty strings
        self.assertEqual(len(fields), 6)  # All 6 fields should be present
        self.assertEqual(fields[3], "")  # d_call
        self.assertEqual(fields[5], "")  # v_score

    def test_write_without_field_order(self):
        """Test writing without pre-defined field order."""
        rec = SeqRecord(Seq("ATGCGT"), id="SEQ001", name="SEQ001", description="")
        rec.annotations["airr"] = {
            "v_call": "IGHV1-69*01",
            "productive": True,
        }
        # No _airr_field_order

        output = StringIO()
        SeqIO.write([rec], output, "airr")

        output.seek(0)
        lines = output.readlines()

        # Should have header with sequence_id, sequence, then sorted others
        header = lines[0].strip().split("\t")
        self.assertEqual(header[0], "sequence_id")
        self.assertEqual(header[1], "sequence")
        self.assertIn("productive", header)
        self.assertIn("v_call", header)

    def test_explicit_fields(self):
        """Test writing with explicitly specified field schema."""
        from Bio.SeqIO.AirrIO import AirrWriter

        # Create records with different fields
        rec1 = SeqRecord(Seq("ATGCGT"), id="SEQ001", name="SEQ001", description="")
        rec1.annotations["airr"] = {
            "v_call": "IGHV1-69*01",
            "productive": True,
            "extra_field": "value1",
        }

        rec2 = SeqRecord(Seq("GCATGC"), id="SEQ002", name="SEQ002", description="")
        rec2.annotations["airr"] = {
            "v_call": "IGHV1-2*02",
            "productive": False,
            "different_field": "value2",
        }

        # Specify fields explicitly
        fields = ["sequence_id", "sequence", "v_call", "productive"]
        output = StringIO()
        writer = AirrWriter(output, fields=fields)
        writer.write_record(rec1)
        writer.write_record(rec2)

        # Check output
        output.seek(0)
        lines = output.readlines()

        # Should have exactly the specified fields
        header = lines[0].strip().split("\t")
        self.assertEqual(header, fields)

        # extra_field and different_field should not be in output
        self.assertNotIn("extra_field", header)
        self.assertNotIn("different_field", header)

    def test_explicit_fields_missing_required(self):
        """Test that missing required fields raises error."""
        from Bio.SeqIO.AirrIO import AirrWriter

        output = StringIO()

        # Missing sequence_id
        with self.assertRaises(ValueError) as cm:
            AirrWriter(output, fields=["sequence", "v_call"])
        self.assertIn("sequence_id", str(cm.exception))

        # Missing sequence
        with self.assertRaises(ValueError) as cm:
            AirrWriter(output, fields=["sequence_id", "v_call"])
        self.assertIn("sequence", str(cm.exception))

    def test_inferred_schema_warns(self):
        """Test that inferred schema issues warning."""
        import warnings

        from Bio.SeqIO.AirrIO import AirrWriter

        rec = SeqRecord(Seq("ATGCGT"), id="SEQ001", name="SEQ001", description="")
        rec.annotations["airr"] = {"v_call": "IGHV1-69*01"}

        output = StringIO()

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            writer = AirrWriter(output)  # No fields specified
            writer.write_record(rec)

            # Should have issued warning
            self.assertEqual(len(w), 1)
            self.assertIn("not explicitly specified", str(w[0].message))
            self.assertIn("Inferring field order", str(w[0].message))

    def test_field_mismatch_warns(self):
        """Test warning when record has fields not in schema."""
        import warnings

        from Bio.SeqIO.AirrIO import AirrWriter

        # First record defines schema
        rec1 = SeqRecord(Seq("ATGCGT"), id="SEQ001", name="SEQ001", description="")
        rec1.annotations["airr"] = {"v_call": "IGHV1-69*01"}

        # Second record has extra field
        rec2 = SeqRecord(Seq("GCATGC"), id="SEQ002", name="SEQ002", description="")
        rec2.annotations["airr"] = {
            "v_call": "IGHV1-2*02",
            "extra_field": "value",  # Not in first record
        }

        output = StringIO()

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            writer = AirrWriter(output)  # Schema inferred from first record
            writer.write_record(rec1)  # First warning: schema inference
            writer.write_record(rec2)  # Second warning: field mismatch

            # Should have two warnings
            self.assertEqual(len(w), 2)

            # First warning about schema inference
            self.assertIn("not explicitly specified", str(w[0].message))

            # Second warning about field mismatch
            self.assertIn("extra_field", str(w[1].message))
            self.assertIn("not in the output schema", str(w[1].message))

    def test_explicit_fields_no_warnings(self):
        """Test that explicit fields don't produce warnings."""
        import warnings

        from Bio.SeqIO.AirrIO import AirrWriter

        rec1 = SeqRecord(Seq("ATGCGT"), id="SEQ001", name="SEQ001", description="")
        rec1.annotations["airr"] = {"v_call": "IGHV1-69*01", "extra_field": "value1"}

        rec2 = SeqRecord(Seq("GCATGC"), id="SEQ002", name="SEQ002", description="")
        rec2.annotations["airr"] = {
            "v_call": "IGHV1-2*02",
            "different_field": "value2",
        }

        # Explicit fields
        fields = ["sequence_id", "sequence", "v_call"]
        output = StringIO()

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            writer = AirrWriter(output, fields=fields)
            writer.write_record(rec1)
            writer.write_record(rec2)

            # Should have NO warnings (extra fields silently dropped with explicit schema)
            self.assertEqual(len(w), 0)


class TestAirrRoundTrip(unittest.TestCase):
    """Tests for round-trip reading and writing."""

    def test_roundtrip_minimal(self):
        """Test round-trip with minimal file."""
        # Read original
        records = list(SeqIO.parse("Airr/minimal.tsv", "airr"))

        # Write to string
        output = StringIO()
        SeqIO.write(records, output, "airr")

        # Read back
        output.seek(0)
        records2 = list(SeqIO.parse(output, "airr"))

        # Check same number of records
        self.assertEqual(len(records), len(records2))

        # Check each record
        for rec1, rec2 in zip(records, records2):
            self.assertEqual(rec1.id, rec2.id)
            self.assertEqual(str(rec1.seq), str(rec2.seq))
            self.assertEqual(rec1.annotations["airr"], rec2.annotations["airr"])

    def test_roundtrip_standard(self):
        """Test round-trip with standard file."""
        records = list(SeqIO.parse("Airr/standard.tsv", "airr"))

        output = StringIO()
        SeqIO.write(records, output, "airr")

        output.seek(0)
        records2 = list(SeqIO.parse(output, "airr"))

        self.assertEqual(len(records), len(records2))

        for rec1, rec2 in zip(records, records2):
            self.assertEqual(rec1.id, rec2.id)
            self.assertEqual(str(rec1.seq), str(rec2.seq))
            # Compare AIRR annotations (should be identical)
            self.assertEqual(rec1.annotations["airr"], rec2.annotations["airr"])

    def test_roundtrip_preserves_field_order(self):
        """Test that field order is preserved through round-trip."""
        records = list(SeqIO.parse("Airr/custom_columns.tsv", "airr"))
        original_order = records[0].annotations["_airr_field_order"]

        output = StringIO()
        SeqIO.write(records, output, "airr")

        output.seek(0)
        records2 = list(SeqIO.parse(output, "airr"))

        # Field order should be preserved
        self.assertEqual(original_order, records2[0].annotations["_airr_field_order"])


class TestAirrIndexing(unittest.TestCase):
    """Tests for random access indexing."""

    def test_index_basic(self):
        """Test basic indexing functionality."""
        idx = SeqIO.index("Airr/minimal.tsv", "airr")

        # Check length
        self.assertEqual(len(idx), 3)

        # Check keys
        self.assertIn("SEQ001", idx)
        self.assertIn("SEQ002", idx)
        self.assertIn("SEQ003", idx)

        # Access by key
        rec = idx["SEQ001"]
        self.assertEqual(rec.id, "SEQ001")
        self.assertEqual(str(rec.seq), "ATGCGTATCGATCGCGATACGATTAGGCGGATGTAA")

        # Check annotations
        self.assertEqual(rec.annotations["airr"]["productive"], True)

        idx.close()

    def test_index_standard(self):
        """Test indexing with standard file."""
        idx = SeqIO.index("Airr/standard.tsv", "airr")

        self.assertEqual(len(idx), 10)

        # Random access
        rec = idx["SEQ005"]
        self.assertEqual(rec.id, "SEQ005")
        self.assertEqual(rec.annotations["airr"]["v_call"], "IGHV3-23*01")

        # Access another
        rec = idx["SEQ010"]
        self.assertEqual(rec.id, "SEQ010")

        idx.close()

    def test_index_large_file(self):
        """Test indexing with large file (1000+ sequences)."""
        idx = SeqIO.index("Airr/large.tsv", "airr")

        # Check we can index 1000 sequences
        self.assertEqual(len(idx), 1000)

        # Random access to various records
        rec = idx["SEQ0001"]
        self.assertEqual(rec.id, "SEQ0001")

        rec = idx["SEQ0500"]
        self.assertEqual(rec.id, "SEQ0500")

        rec = idx["SEQ1000"]
        self.assertEqual(rec.id, "SEQ1000")

        idx.close()

    def test_index_iteration(self):
        """Test iterating over index."""
        idx = SeqIO.index("Airr/minimal.tsv", "airr")

        ids = set(idx.keys())
        self.assertEqual(ids, {"SEQ001", "SEQ002", "SEQ003"})

        # Iterate over records
        count = 0
        for key in idx:
            rec = idx[key]
            self.assertEqual(rec.id, key)
            count += 1

        self.assertEqual(count, 3)

        idx.close()


class TestAirrTypeConversion(unittest.TestCase):
    """Tests for type conversion."""

    def test_boolean_conversion(self):
        """Test boolean field conversion."""
        records = list(SeqIO.parse("Airr/standard.tsv", "airr"))

        # Check True values
        self.assertIs(records[0].annotations["airr"]["productive"], True)
        self.assertIs(records[0].annotations["airr"]["rev_comp"], False)

        # Check False values
        self.assertIs(records[2].annotations["airr"]["productive"], False)

    def test_integer_conversion(self):
        """Test integer field conversion."""
        records = list(SeqIO.parse("Airr/standard.tsv", "airr"))
        airr = records[0].annotations["airr"]

        self.assertIsInstance(airr["v_sequence_start"], int)
        self.assertEqual(airr["v_sequence_start"], 1)
        self.assertIsInstance(airr["v_sequence_end"], int)
        self.assertEqual(airr["v_sequence_end"], 87)
        self.assertIsInstance(airr["junction_length"], int)
        self.assertEqual(airr["junction_length"], 65)

    def test_float_conversion(self):
        """Test float field conversion."""
        records = list(SeqIO.parse("Airr/standard.tsv", "airr"))
        airr = records[0].annotations["airr"]

        self.assertIsInstance(airr["v_identity"], float)
        self.assertAlmostEqual(airr["v_identity"], 0.97)
        self.assertIsInstance(airr["v_score"], float)
        self.assertAlmostEqual(airr["v_score"], 487.5)

    def test_cigar_string_preservation(self):
        """Test that CIGAR strings are kept as strings, not converted to float."""
        # Create record with CIGAR strings
        rec = SeqRecord(Seq("ATGCGT"), id="SEQ001", name="SEQ001", description="")
        rec.annotations["airr"] = {
            "v_call": "IGHV1-69*01",
            "v_cigar": "100M2D50M",
            "d_cigar": "50M",
            "j_cigar": "25M1I25M",
        }
        rec.annotations["_airr_field_order"] = [
            "sequence_id",
            "sequence",
            "v_call",
            "v_cigar",
            "d_cigar",
            "j_cigar",
        ]

        # Write and read back
        output = StringIO()
        SeqIO.write([rec], output, "airr")
        output.seek(0)
        records = list(SeqIO.parse(output, "airr"))

        # Verify CIGAR strings are preserved as strings
        airr = records[0].annotations["airr"]
        self.assertIsInstance(airr["v_cigar"], str)
        self.assertEqual(airr["v_cigar"], "100M2D50M")
        self.assertIsInstance(airr["d_cigar"], str)
        self.assertEqual(airr["d_cigar"], "50M")
        self.assertIsInstance(airr["j_cigar"], str)
        self.assertEqual(airr["j_cigar"], "25M1I25M")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
