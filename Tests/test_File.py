# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.File module."""

import os.path
import shutil
import tempfile
import unittest
import pytest
from io import StringIO

from Bio import bgzf
from Bio import File


class RandomAccess(unittest.TestCase):
    """Random access tests."""

    def test_plain(self):
        """Test plain text file."""
        with File._open_for_random_access("Quality/example.fastq") as handle:
            assert "r" in handle.mode
            assert "b" in handle.mode

    def test_bgzf(self):
        """Test BGZF compressed file."""
        with File._open_for_random_access("Quality/example.fastq.bgz") as handle:
            assert isinstance(handle, bgzf.BgzfReader)

    def test_gzip(self):
        """Test gzip compressed file."""
        with pytest.raises(ValueError):
            File._open_for_random_access("Quality/example.fastq.gz")


class AsHandleTestCase(unittest.TestCase):
    """Tests for as_handle function."""

    def setUp(self):
        """Initialise temporary directory."""
        # Create a directory to work in
        self.temp_dir = tempfile.mkdtemp(prefix="biopython-test")

    def tearDown(self):
        """Remove temporary directory."""
        shutil.rmtree(self.temp_dir)

    def _path(self, *args):
        return os.path.join(self.temp_dir, *args)

    def test_handle(self):
        """Test as_handle with a file-like object argument."""
        p = self._path("test_file.fasta")
        with open(p, "wb") as fp:
            with File.as_handle(fp) as handle:
                assert fp == handle, "as_handle should return argument when given a file-like object"
                assert not handle.closed

            assert not handle.closed, "Exiting as_handle given a file-like object should not close the file"

    def test_string_path(self):
        """Test as_handle with a string path argument."""
        p = self._path("test_file.fasta")
        mode = "wb"
        with File.as_handle(p, mode=mode) as handle:
            assert p == handle.name
            assert mode == handle.mode
            assert not handle.closed
        assert handle.closed

    def test_path_object(self):
        """Test as_handle with a pathlib.Path object."""
        from pathlib import Path

        p = Path(self._path("test_file.fasta"))
        mode = "wb"
        with File.as_handle(p, mode=mode) as handle:
            assert str(p.absolute()) == handle.name
            assert mode == handle.mode
            assert not handle.closed
        assert handle.closed

    def test_custom_path_like_object(self):
        """Test as_handle with a custom path-like object."""

        class CustomPathLike:
            def __init__(self, path):
                self.path = path

            def __fspath__(self):
                return self.path

        p = CustomPathLike(self._path("test_file.fasta"))
        mode = "wb"
        with File.as_handle(p, mode=mode) as handle:
            assert p.path == handle.name
            assert mode == handle.mode
            assert not handle.closed
        assert handle.closed

    def test_stringio(self):
        """Testing passing StringIO handles."""
        s = StringIO()
        with File.as_handle(s) as handle:
            assert s is handle


class BaseClassTests(unittest.TestCase):
    """Tests for _IndexedSeqFileProxy base class."""

    def test_instance_exception(self):
        with pytest.raises(TypeError):
            File._IndexedSeqFileProxy()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
