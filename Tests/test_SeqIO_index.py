# Copyright 2009-2017 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Unit tests for Bio.SeqIO.index(...) and index_db() functions."""

try:
    import sqlite3
except ImportError:
    # Try to run what tests we can in case sqlite3 was not installed
    sqlite3 = None

import gzip
import os
import tempfile
import threading
import unittest
import pytest
import warnings
from io import BytesIO
from io import StringIO
from pathlib import Path

from seq_tests_common import SeqRecordTestBaseClass
from test_SeqIO import SeqIOTestBaseClass

from Bio import BiopythonParserWarning
from Bio import SeqIO
from Bio.SeqIO._index import _FormatToRandomAccess
from Bio.SeqRecord import SeqRecord

CUR_DIR = os.getcwd()


if sqlite3:

    def raw_filenames(index_filename):
        """Open SQLite index and extract filenames (as is).

        Returns a 2-tuple, holding a list of strings, and the value
        of the meta_data.filenames_relative_to_index (or None).
        """
        con = sqlite3.dbapi2.connect(index_filename)

        filenames = [
            row[0]
            for row in con.execute(
                "SELECT name FROM file_data ORDER BY file_number;"
            ).fetchall()
        ]

        try:
            (filenames_relative_to_index,) = con.execute(
                "SELECT value FROM meta_data WHERE key=?;",
                ("filenames_relative_to_index",),
            ).fetchone()
            filenames_relative_to_index = filenames_relative_to_index.upper() == "TRUE"
        except TypeError:
            filenames_relative_to_index = None

        con.close()
        return filenames, filenames_relative_to_index

    class OldIndexTest(unittest.TestCase):
        """Testing a pre-built index (make sure cross platform etc).

        >>> from Bio import SeqIO
        >>> d = SeqIO.index_db("triple_sff.idx", ["E3MFGYR02_no_manifest.sff", "greek.sff", "paired.sff"], "sff")
        >>> len(d)
        54
        """

        def setUp(self):
            os.chdir(CUR_DIR)

        def tearDown(self):
            os.chdir(CUR_DIR)

        def test_old(self):
            """Load existing index with no options (from parent directory)."""
            d = SeqIO.index_db("Roche/triple_sff.idx")
            assert 54 == len(d)
            with pytest.raises(FileNotFoundError):
                d.get_raw("alpha")

        def test_pathobj(self):
            """Load existing index from a pathlib.Path object."""
            d = SeqIO.index_db(Path("Roche/triple_sff.idx"))
            assert 54 == len(d)

        def test_old_check_same_thread(self):
            """Setting check_same_thread to False doesn't raise an exception."""
            d = SeqIO.index_db("Roche/triple_sff_rel_paths.idx")

            def reader_thread():
                try:
                    d["alpha"]
                except sqlite3.ProgrammingError:
                    raise AssertionError("Raised sqlite3.ProgrammingError in violation of check_same_thread=False")

            reader = threading.Thread(target=reader_thread)
            reader.start()
            reader.join()

        def test_old_rel(self):
            """Load existing index (with relative paths) with no options (from parent directory)."""
            d = SeqIO.index_db("Roche/triple_sff_rel_paths.idx")
            assert 54 == len(d)
            assert 395 == len(d["alpha"])

        def test_old_contents(self):
            """Check actual filenames in existing indexes."""
            filenames, flag = raw_filenames("Roche/triple_sff.idx")
            assert flag is None
            assert filenames == ["E3MFGYR02_no_manifest.sff", "greek.sff", "paired.sff"]

            filenames, flag = raw_filenames("Roche/triple_sff_rel_paths.idx")
            assert flag
            assert filenames == ["E3MFGYR02_no_manifest.sff", "greek.sff", "paired.sff"]

        def test_old_same_dir(self):
            """Load existing index with no options (from same directory)."""
            os.chdir("Roche")
            d = SeqIO.index_db("triple_sff.idx")
            assert 54 == len(d)
            assert 395 == len(d["alpha"])

        def test_old_same_dir_rel(self):
            """Load existing index (with relative paths) with no options (from same directory)."""
            os.chdir("Roche")
            d = SeqIO.index_db("triple_sff_rel_paths.idx")
            assert 54 == len(d)
            assert 395 == len(d["alpha"])

        def test_old_format(self):
            """Load existing index with correct format."""
            d = SeqIO.index_db("Roche/triple_sff.idx", format="sff")
            assert 54 == len(d)

        def test_old_format_wrong(self):
            """Load existing index with wrong format."""
            with pytest.raises(ValueError):
                SeqIO.index_db("Roche/triple_sff.idx", format="fasta")

        def test_old_files(self):
            """Load existing index with correct files (from parent directory)."""
            d = SeqIO.index_db(
                "Roche/triple_sff.idx",
                ["E3MFGYR02_no_manifest.sff", "greek.sff", "paired.sff"],
            )
            assert 54 == len(d)
            with pytest.raises(FileNotFoundError):
                d.get_raw("alpha")

        def test_old_files_same_dir(self):
            """Load existing index with correct files (from same directory)."""
            os.chdir("Roche")
            d = SeqIO.index_db(
                "triple_sff.idx",
                ["E3MFGYR02_no_manifest.sff", "greek.sff", "paired.sff"],
            )
            assert 54 == len(d)
            assert 395 == len(d["alpha"])

        def test_old_files_wrong(self):
            """Load existing index with wrong files."""
            with pytest.raises(ValueError):
                SeqIO.index_db("Roche/triple_sff.idx", ["a.sff", "b.sff", "c.sff"])

        def test_old_files_wrong2(self):
            """Load existing index with wrong number of files."""
            with pytest.raises(ValueError):
                SeqIO.index_db("Roche/triple_sff.idx", ["E3MFGYR02_no_manifest.sff", "greek.sff"])

    class NewIndexTest(unittest.TestCase):
        """Check paths etc in newly built index."""

        def setUp(self):
            os.chdir(CUR_DIR)

        def tearDown(self):
            os.chdir(CUR_DIR)
            for i in ["temp.idx", "Roche/temp.idx"]:
                if os.path.isfile(i):
                    os.remove(i)

        def check(self, index_file, sff_files, expt_sff_files):
            if os.path.isfile(index_file):
                os.remove(index_file)
            # Build index...
            d = SeqIO.index_db(index_file, sff_files, "sff")
            assert 395 == len(d["alpha"])
            d._con.close()  # hack for PyPy
            d.close()
            assert [os.path.abspath(f) for f in sff_files] == [os.path.abspath(f) for f in d._filenames]

            # Now directly check the filenames inside the SQLite index:
            filenames, flag = raw_filenames(index_file)
            assert flag
            assert filenames == expt_sff_files

            # Load index...
            d = SeqIO.index_db(index_file, sff_files)
            assert 395 == len(d["alpha"])
            d._con.close()  # hack for PyPy
            d.close()
            assert [os.path.abspath(f) for f in sff_files] == d._filenames

            os.remove(index_file)

        def test_child_folder_rel(self):
            """Check relative links to child folder."""
            # Note we expect relative paths recorded with Unix slashes!
            expt_sff_files = [
                "Roche/E3MFGYR02_no_manifest.sff",
                "Roche/greek.sff",
                "Roche/paired.sff",
            ]

            self.check("temp.idx", expt_sff_files, expt_sff_files)
            # Here index is given as abs
            self.check(
                os.path.abspath("temp.idx"),
                [
                    "Roche/E3MFGYR02_no_manifest.sff",
                    os.path.abspath("Roche/greek.sff"),
                    "Roche/paired.sff",
                ],
                expt_sff_files,
            )
            # Here index is given as relative path
            self.check(
                "temp.idx",
                [
                    "Roche/E3MFGYR02_no_manifest.sff",
                    os.path.abspath("Roche/greek.sff"),
                    "Roche/paired.sff",
                ],
                expt_sff_files,
            )

        def test_same_folder(self):
            """Check relative links in same folder."""
            os.chdir("Roche")
            expt_sff_files = ["E3MFGYR02_no_manifest.sff", "greek.sff", "paired.sff"]

            # Here everything is relative,
            self.check("temp.idx", expt_sff_files, expt_sff_files)
            self.check(
                os.path.abspath("temp.idx"),
                [
                    "E3MFGYR02_no_manifest.sff",
                    os.path.abspath("greek.sff"),
                    "../Roche/paired.sff",
                ],
                expt_sff_files,
            )
            self.check(
                "temp.idx",
                [
                    "E3MFGYR02_no_manifest.sff",
                    os.path.abspath("greek.sff"),
                    "../Roche/paired.sff",
                ],
                expt_sff_files,
            )
            self.check(
                "../Roche/temp.idx",
                [
                    "E3MFGYR02_no_manifest.sff",
                    os.path.abspath("greek.sff"),
                    "../Roche/paired.sff",
                ],
                expt_sff_files,
            )

        def test_some_abs(self):
            """Check absolute filenames in index.

            Unless the repository and tests themselves are under the temp
            directory (as detected by ``tempfile``), we expect the index to
            use absolute filenames.
            """
            h, t = tempfile.mkstemp(prefix="index_test_", suffix=".idx")
            os.close(h)
            os.remove(t)

            abs_sff_files = [
                os.path.abspath("Roche/E3MFGYR02_no_manifest.sff"),
                os.path.abspath("Roche/greek.sff"),
                os.path.abspath(os.path.join("Roche", "paired.sff")),
            ]

            if os.getcwd().startswith(os.path.dirname(t)):
                # The tests are being run from within the temp directory,
                # e.g. index filename /tmp/index_test_XYZ.idx
                # and working directory of /tmp/biopython/Tests/
                # This means the indexing will use a RELATIVE path
                # e.g. biopython/Tests/Roche/E3MFGYR02_no_manifest.sff
                # not /tmp/biopython/Tests/Roche/E3MFGYR02_no_manifest.sff
                expt_sff_files = [
                    os.path.relpath(f, os.path.dirname(t)) for f in abs_sff_files
                ]
            else:
                expt_sff_files = abs_sff_files

            # Providing absolute paths...
            self.check(t, abs_sff_files, expt_sff_files)
            # Now try with mix of abs and relative paths...
            self.check(
                t,
                [
                    os.path.abspath("Roche/E3MFGYR02_no_manifest.sff"),
                    os.path.join("Roche", "greek.sff"),
                    os.path.abspath("Roche/paired.sff"),
                ],
                expt_sff_files,
            )


class IndexDictTests(SeqRecordTestBaseClass, SeqIOTestBaseClass):
    tests = [
        ("Ace/contig1.ace", "ace"),
        ("Ace/consed_sample.ace", "ace"),
        ("Ace/seq.cap.ace", "ace"),
        ("Quality/wrapping_original_sanger.fastq", "fastq"),
        ("Quality/example.fastq", "fastq"),  # Unix newlines
        ("Quality/example.fastq", "fastq-sanger"),
        ("Quality/example_dos.fastq", "fastq"),  # DOS/Windows newlines
        ("Quality/tricky.fastq", "fastq"),
        ("Quality/sanger_faked.fastq", "fastq-sanger"),
        ("Quality/solexa_faked.fastq", "fastq-solexa"),
        ("Quality/illumina_faked.fastq", "fastq-illumina"),
        ("Quality/zero_length.fastq", "fastq"),
        ("EMBL/epo_prt_selection.embl", "embl"),
        ("EMBL/U87107.embl", "embl"),
        ("EMBL/TRBG361.embl", "embl"),
        ("EMBL/kipo_prt_sample.embl", "embl"),
        ("EMBL/A04195.imgt", "embl"),  # Not a proper EMBL file, an IMGT file
        ("EMBL/A04195.imgt", "imgt"),
        ("EMBL/hla_3260_sample.imgt", "imgt"),
        ("EMBL/patents.embl", "embl"),
        ("EMBL/AAA03323.embl", "embl"),
        ("GenBank/NC_000932.faa", "fasta"),
        ("GenBank/NC_005816.faa", "fasta"),
        ("GenBank/NC_005816.tsv", "tab"),
        ("GenBank/NC_005816.ffn", "fasta"),
        ("GenBank/NC_005816.fna", "fasta"),
        ("GenBank/NC_005816.gb", "gb"),
        ("GenBank/cor6_6.gb", "genbank"),
        ("GenBank/empty_accession.gbk", "gb"),
        ("GenBank/empty_version.gbk", "gb"),
        ("IntelliGenetics/vpu_nucaligned.txt", "ig"),
        ("IntelliGenetics/TAT_mase_nuc.txt", "ig"),
        ("IntelliGenetics/VIF_mase-pro.txt", "ig"),
        ("Phd/phd1", "phd"),
        ("Phd/phd2", "phd"),
        ("Phd/phd_solexa", "phd"),
        ("Phd/phd_454", "phd"),
        ("NBRF/B_nuc.pir", "pir"),
        ("NBRF/Cw_prot.pir", "pir"),
        ("NBRF/clustalw.pir", "pir"),
        ("SwissProt/Q13454.txt", "swiss"),
        ("SwissProt/Q13639.txt", "swiss"),
        ("SwissProt/sp016", "swiss"),
        ("SwissProt/multi_ex.txt", "swiss"),
        ("SwissProt/multi_ex.xml", "uniprot-xml"),
        ("SwissProt/multi_ex.fasta", "fasta"),
        ("Roche/E3MFGYR02_random_10_reads.sff", "sff"),
        ("Roche/E3MFGYR02_random_10_reads.sff", "sff-trim"),
        ("Roche/E3MFGYR02_index_at_start.sff", "sff"),
        ("Roche/E3MFGYR02_index_in_middle.sff", "sff"),
        ("Roche/E3MFGYR02_alt_index_at_start.sff", "sff"),
        ("Roche/E3MFGYR02_alt_index_in_middle.sff", "sff"),
        ("Roche/E3MFGYR02_alt_index_at_end.sff", "sff"),
        ("Roche/E3MFGYR02_no_manifest.sff", "sff"),
        ("Roche/greek.sff", "sff"),
        ("Roche/greek.sff", "sff-trim"),
        ("Roche/paired.sff", "sff"),
        ("Roche/paired.sff", "sff-trim"),
    ]

    def setUp(self):
        os.chdir(CUR_DIR)
        h, self.index_tmp = tempfile.mkstemp("_idx.tmp")
        os.close(h)

    def tearDown(self):
        os.chdir(CUR_DIR)
        if os.path.isfile(self.index_tmp):
            os.remove(self.index_tmp)

    def check_dict_methods(self, rec_dict, keys, ids, msg):
        assert sorted(keys) == sorted(rec_dict.keys()), msg
        # This is redundant, I just want to make sure len works:
        assert len(keys) == len(rec_dict), msg
        # Make sure boolean evaluation works
        assert bool(keys) == bool(rec_dict), msg
        for key, id in zip(keys, ids):
            assert key in rec_dict, msg
            assert id == rec_dict[key].id, msg
            assert id == rec_dict.get(key).id, msg
        # Check non-existent keys,
        assert chr(0) not in keys, "Bad example in test"
        with pytest.raises(KeyError):
            rec = rec_dict[chr(0)]
        assert rec_dict.get(chr(0)) is None, msg
        assert rec_dict.get(chr(0), chr(1)) == chr(1), msg
        with pytest.raises(AttributeError):
            rec_dict.iteritems
        for key, rec in rec_dict.items():
            assert key in keys, msg
            assert isinstance(rec, SeqRecord), msg
            assert rec.id in ids, msg
        for rec in rec_dict.values():
            assert key in keys, msg
            assert isinstance(rec, SeqRecord), msg
            assert rec.id in ids, msg

    def simple_check(self, filename, fmt, comp):
        """Check indexing (without a key function)."""
        msg = f"Test failure parsing file {filename} with format {fmt}"
        if comp:
            mode = "r" + self.get_mode(fmt)
            with gzip.open(filename, mode) as handle:
                id_list = [rec.id for rec in SeqIO.parse(handle, fmt)]
        else:
            id_list = [rec.id for rec in SeqIO.parse(filename, fmt)]

        with warnings.catch_warnings():
            if "_alt_index_" in str(filename):
                # BiopythonParserWarning: Could not parse the SFF index:
                # Unknown magic number b'.diy' in SFF index header:
                # b'.diy1.00'
                warnings.simplefilter("ignore", BiopythonParserWarning)

            rec_dict = SeqIO.index(filename, fmt)
            self.check_dict_methods(rec_dict, id_list, id_list, msg=msg)
            rec_dict.close()

            if not sqlite3:
                return

            # In memory,
            # note here give filenames as list of strings
            rec_dict = SeqIO.index_db(":memory:", [filename], fmt)
            self.check_dict_methods(rec_dict, id_list, id_list, msg=msg)
            rec_dict.close()

            # check error conditions
            with pytest.raises(ValueError):
                SeqIO.index_db(":memory:", format="dummy")
            with pytest.raises(ValueError):
                SeqIO.index_db(":memory:", filenames=["dummy"])

            # Saving to file...
            index_tmp = self.index_tmp
            if os.path.isfile(index_tmp):
                os.remove(index_tmp)

            # To disk,
            # note here we give the filename as a single string
            # to confirm that works too.
            rec_dict = SeqIO.index_db(index_tmp, filename, fmt)
            self.check_dict_methods(rec_dict, id_list, id_list, msg=msg)
            rec_dict.close()
            rec_dict._con.close()  # hack for PyPy

            # Now reload it...
            rec_dict = SeqIO.index_db(index_tmp, [filename], fmt)
            self.check_dict_methods(rec_dict, id_list, id_list, msg=msg)
            rec_dict.close()
            rec_dict._con.close()  # hack for PyPy

            # Now reload without passing filenames and format
            # and switch directory to check  paths still work
            index_tmp = os.path.abspath(index_tmp)
            os.chdir(os.path.dirname(filename))
            try:
                rec_dict = SeqIO.index_db(index_tmp)
            finally:
                os.chdir(CUR_DIR)
            self.check_dict_methods(rec_dict, id_list, id_list, msg=msg)
            rec_dict.close()
            rec_dict._con.close()  # hack for PyPy

            os.remove(index_tmp)

    def add_prefix(self, key):
        """Sample key_function for testing index code."""
        return "id_" + key

    def key_check(self, filename, fmt, comp):
        """Check indexing with a key function."""
        msg = f"Test failure parsing file {filename} with format {fmt}"
        if comp:
            mode = "r" + self.get_mode(fmt)
            with gzip.open(filename, mode) as handle:
                id_list = [rec.id for rec in SeqIO.parse(handle, fmt)]
        else:
            id_list = [rec.id for rec in SeqIO.parse(filename, fmt)]

        key_list = [self.add_prefix(id) for id in id_list]

        with warnings.catch_warnings():
            if "_alt_index_" in str(filename):
                # BiopythonParserWarning: Could not parse the SFF index:
                # Unknown magic number b'.diy' in SFF index header:
                # b'.diy1.00'
                warnings.simplefilter("ignore", BiopythonParserWarning)

            rec_dict = SeqIO.index(filename, fmt, key_function=self.add_prefix)
            self.check_dict_methods(rec_dict, key_list, id_list, msg=msg)
            rec_dict.close()

            if not sqlite3:
                return

            # In memory,
            rec_dict = SeqIO.index_db(
                ":memory:", [filename], fmt, key_function=self.add_prefix
            )
            self.check_dict_methods(rec_dict, key_list, id_list, msg=msg)
            # check error conditions
            with pytest.raises(ValueError):
                SeqIO.index_db(":memory:", format="dummy", key_function=self.add_prefix)
            with pytest.raises(ValueError):
                SeqIO.index_db(
                    ":memory:", filenames=["dummy"], key_function=self.add_prefix
                )
            rec_dict.close()

            # Saving to file...
            index_tmp = None
            if isinstance(filename, str):
                index_tmp = filename + ".key.idx"
            elif isinstance(filename, Path):
                index_tmp = filename.with_suffix(".key.idx")

            if os.path.isfile(str(index_tmp)):
                os.remove(str(index_tmp))
            rec_dict = SeqIO.index_db(
                index_tmp, [filename], fmt, key_function=self.add_prefix
            )
            self.check_dict_methods(rec_dict, key_list, id_list, msg=msg)
            rec_dict.close()
            rec_dict._con.close()  # hack for PyPy

            # Now reload it...
            rec_dict = SeqIO.index_db(
                index_tmp, [filename], fmt, key_function=self.add_prefix
            )
            self.check_dict_methods(rec_dict, key_list, id_list, msg=msg)
            rec_dict.close()
            rec_dict._con.close()  # hack for PyPy

            # Now reload without passing filenames and format
            rec_dict = SeqIO.index_db(index_tmp, key_function=self.add_prefix)
            self.check_dict_methods(rec_dict, key_list, id_list, msg=msg)
            rec_dict.close()
            rec_dict._con.close()  # hack for PyPy
            os.remove(str(index_tmp))
            # Done

    def get_raw_check(self, filename, fmt, comp):
        # Also checking the key_function here
        msg = f"Test failure parsing file {filename} with format {fmt}"
        if comp:
            with gzip.open(filename, "rb") as handle:
                raw_file = handle.read()
            mode = "r" + self.get_mode(fmt)
            with gzip.open(filename, mode) as handle:
                id_list = [rec.id.lower() for rec in SeqIO.parse(handle, fmt)]
        else:
            with open(filename, "rb") as handle:
                raw_file = handle.read()
            id_list = [rec.id.lower() for rec in SeqIO.parse(filename, fmt)]

        if fmt in ["sff"]:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", BiopythonParserWarning)
                rec_dict = SeqIO.index(filename, fmt, key_function=str.lower)
                if sqlite3:
                    rec_dict_db = SeqIO.index_db(
                        ":memory:", filename, fmt, key_function=str.lower
                    )
        else:
            rec_dict = SeqIO.index(filename, fmt, key_function=str.lower)
            if sqlite3:
                rec_dict_db = SeqIO.index_db(
                    ":memory:", filename, fmt, key_function=str.lower
                )

        assert sorted(id_list) == sorted(rec_dict.keys()), msg
        if sqlite3:
            assert sorted(id_list) == sorted(rec_dict_db.keys()), msg
        for key in id_list:
            assert key in rec_dict, msg
            assert key == rec_dict[key].id.lower(), msg
            assert key == rec_dict.get(key).id.lower(), msg
            raw = rec_dict.get_raw(key)
            assert isinstance(raw, bytes), msg
            assert raw.strip(), msg
            assert raw in raw_file, msg

            if sqlite3:
                raw_db = rec_dict_db.get_raw(key)
                # Via index using format-specific get_raw which scans the file,
                # Via index_db in general using raw length found when indexing.
                assert raw == raw_db, msg

            rec1 = rec_dict[key]
            # Following isn't very elegant, but it lets me test the
            # __getitem__ SFF code is working.
            mode = self.get_mode(fmt)
            if mode == "b":
                handle = BytesIO(raw)
            elif mode == "t":
                handle = StringIO(raw.decode())
            else:
                raise RuntimeError(f"Unexpected mode {mode}")
            if fmt == "sff":
                rec_dict._proxy.trim = False
                rec_dict._proxy._offset = handle.tell()
                rec2 = SeqIO.SffIO.SffIterator._sff_read_seq_record(
                    rec_dict._proxy, handle
                )
            elif fmt == "sff-trim":
                rec_dict._proxy.trim = True
                rec_dict._proxy._offset = handle.tell()
                rec2 = SeqIO.SffIO.SffIterator._sff_read_seq_record(
                    rec_dict._proxy, handle
                )
            elif fmt == "uniprot-xml":
                assert raw.startswith(b"<entry "), msg
                assert raw.endswith(b"</entry>"), msg
                # Currently the __getitem__ method uses this
                # trick too, but we hope to fix that later
                raw = (
                    b"""<?xml version='1.0' encoding='UTF-8'?>
                <uniprot xmlns="http://uniprot.org/uniprot"
                xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
                xsi:schemaLocation="http://uniprot.org/uniprot
                http://www.uniprot.org/support/docs/uniprot.xsd">
                %s
                </uniprot>
                """
                    % raw
                )
                handle = BytesIO(raw)
                rec2 = SeqIO.read(handle, fmt)
            else:
                rec2 = SeqIO.read(handle, fmt)
            self.compare_record(rec1, rec2)
        rec_dict.close()
        del rec_dict

    if sqlite3:

        def test_alpha_fails_db(self):
            """Reject alphabet argument in Bio.SeqIO.index_db()."""
            # In historic usage, alphabet=... would be a Bio.Alphabet object.
            with pytest.raises(ValueError):
                SeqIO.index_db(":memory:", ["Fasta/dups.fasta"], "fasta", alphabet="XXX")

    def test_alpha_fails(self):
        """Reject alphabet argument in Bio.SeqIO.index()."""
        # In historic usage, alphabet=... would be a Bio.Alphabet object.
        with pytest.raises(ValueError):
            SeqIO.index("Fasta/dups.fasta", "fasta", alphabet="XXX")

    if sqlite3:

        def test_duplicates_index_db(self):
            """Index file with duplicate identifiers with Bio.SeqIO.index_db()."""
            with pytest.raises(ValueError):
                SeqIO.index_db(":memory:", ["Fasta/dups.fasta"], "fasta")

    def test_duplicates_index(self):
        """Index file with duplicate identifiers with Bio.SeqIO.index()."""
        with pytest.raises(ValueError):
            SeqIO.index("Fasta/dups.fasta", "fasta")

    def test_duplicates_to_dict(self):
        """Index file with duplicate identifiers with Bio.SeqIO.to_dict()."""
        with open("Fasta/dups.fasta") as handle:
            iterator = SeqIO.parse(handle, "fasta")
            with pytest.raises(ValueError):
                SeqIO.to_dict(iterator)

    def test_simple_checks(self):
        for filename1, fmt in self.tests:
            assert fmt in _FormatToRandomAccess
            tasks = [(filename1, None)]
            if os.path.isfile(filename1 + ".bgz"):
                tasks.append((filename1 + ".bgz", "bgzf"))
            for filename2, comp in tasks:
                self.simple_check(filename2, fmt, comp)

    def test_simple_checks_with_pathobj(self):
        for filename1, fmt in self.tests:
            assert fmt in _FormatToRandomAccess
            tasks = [(filename1, None)]
            if os.path.isfile(filename1 + ".bgz"):
                tasks.append((filename1 + ".bgz", "bgzf"))
            for filename2, comp in tasks:
                self.simple_check(Path(filename2), fmt, comp)

    def test_key_checks(self):
        for filename1, fmt in self.tests:
            assert fmt in _FormatToRandomAccess
            tasks = [(filename1, None)]
            if os.path.isfile(filename1 + ".bgz"):
                tasks.append((filename1 + ".bgz", "bgzf"))
            for filename2, comp in tasks:
                self.key_check(filename2, fmt, comp)

    def test_key_checks_with_pathobj(self):
        for filename1, fmt in self.tests:
            assert fmt in _FormatToRandomAccess
            tasks = [(filename1, None)]
            if os.path.isfile(filename1 + ".bgz"):
                tasks.append((filename1 + ".bgz", "bgzf"))
            for filename2, comp in tasks:
                self.key_check(Path(filename2), fmt, comp)

    def test_raw_checks(self):
        for filename1, fmt in self.tests:
            assert fmt in _FormatToRandomAccess
            tasks = [(filename1, None)]
            if os.path.isfile(filename1 + ".bgz"):
                tasks.append((filename1 + ".bgz", "bgzf"))
            for filename2, comp in tasks:
                self.get_raw_check(filename2, fmt, comp)

    def test_raw_checks_with_pathobj(self):
        for filename1, fmt in self.tests:
            assert fmt in _FormatToRandomAccess
            tasks = [(filename1, None)]
            if os.path.isfile(filename1 + ".bgz"):
                tasks.append((filename1 + ".bgz", "bgzf"))
            for filename2, comp in tasks:
                self.get_raw_check(Path(filename2), fmt, comp)


class IndexOrderingSingleFile(unittest.TestCase):
    f = "GenBank/NC_000932.faa"
    ids = [r.id for r in SeqIO.parse(f, "fasta")]

    def test_order_to_dict(self):
        """Check to_dict preserves order in indexed file."""
        d = SeqIO.to_dict(SeqIO.parse(self.f, "fasta"))
        assert self.ids == list(d)

    def test_order_index(self):
        """Check index preserves order in indexed file."""
        d = SeqIO.index(self.f, "fasta")
        assert self.ids == list(d)

    if sqlite3:

        def test_order_index_db(self):
            """Check index_db preserves ordering indexed file."""
            d = SeqIO.index_db(":memory:", [self.f], "fasta")
            assert self.ids == list(d)


if sqlite3:

    class IndexOrderingManyFiles(unittest.TestCase):
        def test_order_index_db(self):
            """Check index_db preserves order in multiple indexed files."""
            files = ["GenBank/NC_000932.faa", "GenBank/NC_005816.faa"]
            ids = []
            for f in files:
                ids.extend(r.id for r in SeqIO.parse(f, "fasta"))
            d = SeqIO.index_db(":memory:", files, "fasta")
            assert ids == list(d)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
