# Copyright 2012 by Jeff Hussmann.  All rights reserved.
# Revisions copyright 2013-2016 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SffIO module."""

import re
import unittest
import pytest
from io import BytesIO

from Bio import SeqIO
from Bio.SeqIO.SffIO import _sff_do_slow_index
from Bio.SeqIO.SffIO import _sff_find_roche_index
from Bio.SeqIO.SffIO import _sff_read_roche_index
from Bio.SeqIO.SffIO import ReadRocheXmlManifest
from Bio.SeqIO.SffIO import SffIterator
from Bio.SeqIO.SffIO import SffWriter

# sffinfo E3MFGYR02_random_10_reads.sff | sed -n '/>\|Run Prefix\|Region\|XY/p'
test_data = """
>E3MFGYR02JWQ7T
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  3946_2103
>E3MFGYR02JA6IL
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  3700_3115
>E3MFGYR02JHD4H
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  3771_2095
>E3MFGYR02GFKUC
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2520_2738
>E3MFGYR02FTGED
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2268_2739
>E3MFGYR02FR9G7
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2255_0361
>E3MFGYR02GAZMS
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2468_1618
>E3MFGYR02HHZ8O
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2958_1574
>E3MFGYR02GPGB1
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2633_0607
>E3MFGYR02F7Z7G
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2434_1658"""


class TestUAN(unittest.TestCase):
    """Test annotations."""

    def setUp(self):
        self.records = list(SeqIO.parse("Roche/E3MFGYR02_random_10_reads.sff", "sff"))
        self.test_annotations = {}
        for line in test_data.splitlines():
            fields = re.split(r"\s+", line.strip())
            if ">" in line:
                current_name = fields[0].lstrip(">")
                self.test_annotations[current_name] = {}
            elif "Prefix" in line:
                time_list = [int(v) for v in fields[2].split("_")[1:-1]]
                self.test_annotations[current_name]["time"] = time_list
            elif "Region" in line:
                region = int(fields[-1])
                self.test_annotations[current_name]["region"] = region
            elif "XY" in line:
                x, y = (int(v) for v in fields[-1].split("_"))
                self.test_annotations[current_name]["coords"] = (x, y)

    def test_time(self):
        for record in self.records:
            assert record.annotations["time"] == self.test_annotations[record.name]["time"]

    def test_region(self):
        for record in self.records:
            assert record.annotations["region"] == self.test_annotations[record.name]["region"]

    def test_coords(self):
        for record in self.records:
            assert record.annotations["coords"] == self.test_annotations[record.name]["coords"]


class TestErrors(unittest.TestCase):
    with open("Roche/E3MFGYR02_random_10_reads.sff", "rb") as handle:
        good = handle.read()

    def test_empty(self):
        fh = BytesIO()
        with pytest.raises(ValueError) as cm:
            records = list(SeqIO.parse(fh, "sff"))
        assert str(cm.value) == "Empty file."

    def check_bad_header(self, header, msg):
        with pytest.raises(ValueError) as cm:
            records = list(SeqIO.parse(BytesIO(header), "sff"))
        err = str(cm.value)
        if isinstance(msg, (tuple, list)):
            assert err in msg, f"Unexpected error: {err}"
        else:
            assert err == msg

    def test_30bytes(self):
        self.check_bad_header(b"x" * 30, "File too small to hold a valid SFF header.")

    def test_31bytes(self):
        self.check_bad_header(
            b"x" * 31,
            (
                "SFF file did not start '.sff', but 'xxxx'",
                "SFF file did not start '.sff', but b'xxxx'",
            ),
        )

    def test_31bytes_index_header(self):
        self.check_bad_header(
            b".srt" + b"x" * 27, "Handle seems to be at SFF index block, not start"
        )

    def test_31bytes_bad_ver(self):
        self.check_bad_header(
            b".sff1.00" + b"x" * 23, "Unsupported SFF version in header, 49.46.48.48"
        )

    def test_31bytes_bad_flowgram(self):
        self.check_bad_header(
            b".sff\x00\x00\x00\x01" + b"x" * 23,
            "Flowgram format code 120 not supported",
        )

    def test_bad_index_offset(self):
        bad = self.good[:12] + b"\x00\x00\x00\x00" + self.good[16:]
        self.check_bad_header(bad, "Index offset 0 but index length 764")

    def test_bad_index_length(self):
        bad = self.good[:16] + b"\x00\x00\x00\x00" + self.good[20:]
        self.check_bad_header(bad, "Index offset 16824 but index length 0")

    def test_bad_index_eof(self):
        # Semi-random edit to the index offset value,
        bad = self.good[:13] + b"\x01" + self.good[14:]
        self.check_bad_header(
            bad,
            "Gap of 65536 bytes after final record end 16824, "
            "before 82360 where index starts?",
        )

    def test_no_index(self):
        # Does a lot of work to create a no-index SFF file
        # (in the process checking this bit of SffWriter works)
        records = list(SeqIO.parse(BytesIO(self.good), "sff"))
        with BytesIO() as handle:
            writer = SffWriter(handle, index=False)
            count = writer.write_file(records)
            assert count == len(records)
            handle.seek(0)
            new = list(SeqIO.parse(handle, "sff"))
            assert len(records) == len(new)
            for a, b in zip(records, new):
                assert a.id == b.id
            handle.seek(0)
            with pytest.raises(ValueError) as cm:
                values = _sff_find_roche_index(handle)
            err = str(cm.value)
            assert err == "No index present in this SFF file"

    def test_unknown_index(self):
        # TODO - Add SFF file with no index,
        # self.assertEqual(str(err), "No index present in this SFF file")
        with open("Roche/E3MFGYR02_alt_index_in_middle.sff", "rb") as handle:
            with pytest.raises(ValueError) as cm:
                values = _sff_find_roche_index(handle)
        assert str(cm.value) in (
                "Unknown magic number '.diy' in SFF index header:\n'.diy1.00'",
                "Unknown magic number b'.diy' in SFF index header:\nb'.diy1.00'",
            )

    def check_sff_read_roche_index(self, data, msg):
        handle = BytesIO(data)
        with pytest.raises(ValueError) as cm:
            index = list(_sff_read_roche_index(handle))
        assert str(cm.value) == msg

    def test_premature_end_of_index(self):
        self.check_sff_read_roche_index(self.good[:-50], "Premature end of file!")

    def test_index_name_no_null(self):
        assert self.good[17502:17503] == b"\x00"
        self.check_sff_read_roche_index(
            self.good[:17502] + b"x" + self.good[17503:],
            "Expected a null terminator to the read name.",
        )

    def test_index_mft_version(self):
        assert self.good[16824:16832] == b".mft1.00"
        self.check_sff_read_roche_index(
            self.good[:16828] + b"\x01\x02\x03\x04" + self.good[16832:],
            "Unsupported version in .mft index header, 1.2.3.4",
        )

    def test_index_mft_data_size(self):
        assert self.good[16824:16832] == b".mft1.00"
        self.check_sff_read_roche_index(
            self.good[:16836] + b"\x00\x00\x00\x00" + self.good[16840:],
            "Problem understanding .mft index header, 764 != 8 + 8 + 548 + 0",
        )

    def test_index_lengths(self):
        # Reduce the number of reads from 10 to 9 so index loading fails...
        assert self.good[20:24] == b"\x00\x00\x00\x0a"
        self.check_sff_read_roche_index(
            self.good[:20] + b"\x00\x00\x00\x09" + self.good[24:],
            "Problem with index length? 17568 vs 17588",
        )

    def test_no_manifest_xml(self):
        with open("Roche/E3MFGYR02_no_manifest.sff", "rb") as handle:
            with pytest.raises(ValueError) as cm:
                xml = ReadRocheXmlManifest(handle)
            assert str(cm.value) == "No XML manifest found"


class TestIndex(unittest.TestCase):
    """Test SFF index."""

    def test_manifest(self):
        filename = "Roche/E3MFGYR02_random_10_reads.sff"
        with open(filename, "rb") as handle:
            metadata = ReadRocheXmlManifest(handle)

    def test_both_ways(self):
        filename = "Roche/E3MFGYR02_random_10_reads.sff"
        with open(filename, "rb") as handle:
            index1 = sorted(_sff_read_roche_index(handle))
        with open(filename, "rb") as handle:
            index2 = sorted(_sff_do_slow_index(handle))
        assert index1 == index2
        with open(filename, "rb") as handle:
            assert len(index1) == len(list(SffIterator(handle)))
        with open(filename, "rb") as handle:
            assert len(index1) == len(list(SffIterator(BytesIO(handle.read()))))


class TestAlternativeIndexes(unittest.TestCase):
    filename = "Roche/E3MFGYR02_random_10_reads.sff"
    with open(filename, "rb") as handle:
        sff = list(SffIterator(handle))

    def check_same(self, new_sff):
        assert len(self.sff) == len(new_sff)
        for old, new in zip(self.sff, new_sff):
            assert old.id == new.id
            assert old.seq == new.seq

    def test_alt_index_at_end(self):
        with open("Roche/E3MFGYR02_alt_index_at_end.sff", "rb") as handle:
            sff2 = list(SffIterator(handle))
        self.check_same(sff2)

    def test_alt_index_at_start(self):
        with open("Roche/E3MFGYR02_alt_index_at_start.sff", "rb") as handle:
            sff2 = list(SffIterator(handle))
        self.check_same(sff2)

    def test_alt_index_in_middle(self):
        with open("Roche/E3MFGYR02_alt_index_in_middle.sff", "rb") as handle:
            sff2 = list(SffIterator(handle))
        self.check_same(sff2)

    def test_index_at_start(self):
        with open("Roche/E3MFGYR02_index_at_start.sff", "rb") as handle:
            sff2 = list(SffIterator(handle))
        self.check_same(sff2)

    def test_index_in_middle(self):
        with open("Roche/E3MFGYR02_index_in_middle.sff", "rb") as handle:
            sff2 = list(SffIterator(handle))
        self.check_same(sff2)

    def test_trim(self):
        with open(self.filename, "rb") as handle:
            sff_trim = list(SffIterator(handle, trim=True))
        assert len(self.sff) == len(sff_trim)
        for old, new in zip(self.sff, sff_trim):
            assert old.id == new.id


class TestConcatenated(unittest.TestCase):
    """Test concatenated SFF files."""

    def test_parses_gzipped_stream(self):
        import gzip

        count = 0
        with gzip.open("Roche/E3MFGYR02_random_10_reads.sff.gz", "rb") as fh:
            for record in SeqIO.parse(fh, "sff"):
                count += 1
        assert 10 == count

    def test_parse1(self):
        count = 0
        caught = False
        try:
            for record in SeqIO.parse("Roche/invalid_greek_E3MFGYR02.sff", "sff"):
                count += 1
        except ValueError as err:
            assert ("Additional data at end of SFF file, perhaps "
                "multiple SFF files concatenated? "
                "See offset 65296" in str(err)), err
            caught = True
        assert caught, "Didn't spot concatenation"
        assert count == 24

    def test_index1(self):
        with pytest.raises(ValueError) as cm:
            d = SeqIO.index("Roche/invalid_greek_E3MFGYR02.sff", "sff")
        err = str(cm.value)
        assert ("Additional data at end of SFF file, perhaps "
            "multiple SFF files concatenated? See offset 65296" in err)

    def test_parse2(self):
        count = 0
        caught = False
        try:
            for record in SeqIO.parse("Roche/invalid_paired_E3MFGYR02.sff", "sff"):
                count += 1
        except ValueError as err:
            assert ("Your SFF file is invalid, post index 5 byte "
                "null padding region ended '.sff' which could "
                "be the start of a concatenated SFF file? "
                "See offset 54371" in str(err)), err
            caught = True
        assert caught, "Didn't spot concatenation"
        assert count == 20

    def test_index2(self):
        with pytest.raises(ValueError) as cm:
            d = SeqIO.index("Roche/invalid_paired_E3MFGYR02.sff", "sff")
        assert ("Your SFF file is invalid, post index 5 byte "
            "null padding region ended '.sff' which could "
            "be the start of a concatenated SFF file? "
            "See offset 54371" in str(cm.value))


class TestSelf(unittest.TestCase):
    """These tests were originally defined in SffIO.py as self-tests."""

    def test_read(self):
        filename = "Roche/E3MFGYR02_random_10_reads.sff"
        with open(filename, "rb") as handle:
            sff = list(SffIterator(handle))
        with open(filename, "rb") as handle:
            sff_trim = list(SffIterator(handle, trim=True))

        filename = "Roche/E3MFGYR02_random_10_reads_no_trim.fasta"
        fasta_no_trim = list(SeqIO.parse(filename, "fasta"))
        filename = "Roche/E3MFGYR02_random_10_reads_no_trim.qual"
        qual_no_trim = list(SeqIO.parse(filename, "qual"))

        filename = "Roche/E3MFGYR02_random_10_reads.fasta"
        fasta_trim = list(SeqIO.parse(filename, "fasta"))
        filename = "Roche/E3MFGYR02_random_10_reads.qual"
        qual_trim = list(SeqIO.parse(filename, "qual"))

        for s, sT, f, q, fT, qT in zip(
            sff, sff_trim, fasta_no_trim, qual_no_trim, fasta_trim, qual_trim
        ):
            assert len({s.id, f.id, q.id}) == 1  # All values are the same
            assert s.seq == f.seq
            assert s.letter_annotations["phred_quality"] == q.letter_annotations["phred_quality"]
            assert len({s.id, sT.id, fT.id, qT.id}) == 1  # All values are the same
            assert sT.seq == fT.seq
            assert sT.letter_annotations["phred_quality"] == qT.letter_annotations["phred_quality"]

    def test_write(self):
        filename = "Roche/E3MFGYR02_random_10_reads.sff"
        with open(filename, "rb") as handle:
            metadata = ReadRocheXmlManifest(handle)
        with open(filename, "rb") as handle:
            sff = list(SffIterator(handle))
        b_handle = BytesIO()
        w = SffWriter(b_handle, xml=metadata)
        w.write_file(sff)  # list
        data = b_handle.getvalue()
        # And again with an iterator...
        handle = BytesIO()
        w = SffWriter(handle, xml=metadata)
        w.write_file(iter(sff))
        assert data == handle.getvalue()
        # Check 100% identical to the original:
        with open(filename, "rb") as handle:
            original = handle.read()
        assert len(data) == len(original)
        assert data == original
        del data

    def test_index(self):
        filename = "Roche/greek.sff"
        with open(filename, "rb") as a_handle, open(filename, "rb") as b_handle:
            index1 = sorted(_sff_read_roche_index(a_handle))
            index2 = sorted(_sff_do_slow_index(b_handle))
            assert index1 == index2

    def test_read_wrong(self):
        filename = "Roche/greek.sff"
        with open(filename, "rb") as handle:
            with pytest.raises(ValueError):
                ReadRocheXmlManifest(handle)

        with open(filename, "rb") as handle:
            for record in SffIterator(handle):
                pass

            def fileiter(handle):
                for record in SffIterator(handle):
                    i = record.id

            with pytest.raises(ValueError):
                fileiter(handle)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
