# Copyright (C) 2011, 2019 by Brandon Invergo (b.invergo@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for PAML yn00 module."""

import glob
import os
import os.path
import unittest
import pytest

from Bio.Phylo.PAML import yn00


class ModTest(unittest.TestCase):
    align_dir = os.path.join("PAML", "Alignments")
    tree_dir = os.path.join("PAML", "Trees")
    ctl_dir = os.path.join("PAML", "Control_files")
    results_dir = os.path.join("PAML", "Results")
    working_dir = os.path.join("PAML", "yn00_test")

    align_file = os.path.join(align_dir, "alignment.phylip")
    dotname_align_file = os.path.join(align_dir, "alignment_dottednames.phylip")
    out_file = os.path.join(results_dir, "test.out")
    results_file = os.path.join(results_dir, "bad_results.out")
    bad_ctl_file1 = os.path.join(ctl_dir, "bad1.ctl")
    bad_ctl_file2 = os.path.join(ctl_dir, "bad2.ctl")
    ctl_file = os.path.join(ctl_dir, "yn00", "yn00.ctl")

    def tearDown(self):
        """Just in case yn00 creates some junk files, do a clean-up."""
        del_files = [self.out_file, "2YN.dN", "2YN.dS", "2YN.t", "rst", "rst1", "rub"]
        for filename in del_files:
            if os.path.exists(filename):
                os.remove(filename)
        if os.path.exists(self.working_dir):
            for filename in os.listdir(self.working_dir):
                filepath = os.path.join(self.working_dir, filename)
                os.remove(filepath)
            os.rmdir(self.working_dir)

    def setUp(self):
        self.yn00 = yn00.Yn00()

    def testAlignmentFileIsValid(self):
        with pytest.raises((AttributeError, TypeError, OSError)):
            yn00.Yn00(alignment=[])
        self.yn00.alignment = []
        self.yn00.out_file = self.out_file
        with pytest.raises((AttributeError, TypeError, OSError)):
            self.yn00.run()

    def testAlignmentExists(self):
        with pytest.raises((EnvironmentError, IOError)):
            yn00.Yn00(alignment="nonexistent")
        self.yn00.alignment = "nonexistent"
        self.yn00.out_file = self.out_file
        with pytest.raises(IOError):
            self.yn00.run()

    def testWorkingDirValid(self):
        self.yn00.alignment = self.align_file
        self.yn00.out_file = self.out_file
        self.yn00.working_dir = []
        with pytest.raises((AttributeError, TypeError, OSError)):
            self.yn00.run()

    def testOptionExists(self):
        with pytest.raises((AttributeError, KeyError)):
            self.yn00.set_options(xxxx=1)
        with pytest.raises((AttributeError, KeyError)):
            self.yn00.get_option("xxxx")

    def testAlignmentSpecified(self):
        self.yn00.out_file = self.out_file
        with pytest.raises((AttributeError, ValueError)):
            self.yn00.run()

    def testOutputFileSpecified(self):
        self.yn00.alignment = self.align_file
        with pytest.raises((AttributeError, ValueError)):
            self.yn00.run()

    # def testPamlErrorsCaught(self):
    #     self.yn00.alignment = self.align_file
    #     self.yn00.out_file = self.out_file
    #     self.assertRaises((EnvironmentError, PamlError),
    #                        self.yn00.run)

    def testCtlFileValidOnRun(self):
        self.yn00.alignment = self.align_file
        self.yn00.out_file = self.out_file
        with pytest.raises((AttributeError, TypeError, OSError)):
            self.yn00.run(ctl_file=[])

    def testCtlFileExistsOnRun(self):
        self.yn00.alignment = self.align_file
        self.yn00.out_file = self.out_file
        with pytest.raises(IOError):
            self.yn00.run(ctl_file="nonexistent")

    def testCtlFileValidOnRead(self):
        with pytest.raises((AttributeError, TypeError, OSError)):
            self.yn00.read_ctl_file([])
        with pytest.raises((AttributeError, KeyError)):
            self.yn00.read_ctl_file(self.bad_ctl_file1)
        with pytest.raises(AttributeError):
            self.yn00.read_ctl_file(self.bad_ctl_file2)
        target_options = {
            "verbose": 1,
            "icode": 0,
            "weighting": 0,
            "commonf3x4": 0,
            "ndata": 1,
        }
        self.yn00.read_ctl_file(self.ctl_file)
        assert self.yn00._options == target_options

    def testCtlFileExistsOnRead(self):
        with pytest.raises(IOError):
            self.yn00.read_ctl_file(ctl_file="nonexistent")

    def testResultsValid(self):
        with pytest.raises((AttributeError, TypeError, OSError)):
            yn00.read([])

    def testResultsExist(self):
        with pytest.raises((EnvironmentError, IOError)):
            yn00.read("nonexistent")

    def testResultsParsable(self):
        with pytest.raises(ValueError):
            yn00.read(self.results_file)

    def testParseAllVersions(self):
        pattern = os.path.join(self.results_dir, "yn00", "yn00-*")
        for results_file in glob.glob(pattern):
            results = yn00.read(results_file)
            assert len(results) == 5
            assert len(results["Homo_sapie"]) == 4
            assert len(results["Homo_sapie"]["Pan_troglo"]) == 5

    def testParseLongNames(self):
        pattern = os.path.join(self.results_dir, "yn00", "yn00_long-*")
        for results_file in glob.glob(pattern):
            results = yn00.read(results_file)
            # Expect seven taxa...
            assert len(results) == 7
            # ...each of which is compared to the other six.
            assert {len(v) for v in results.values()} == {6}
            # ...each of which has five measures.
            assert {len(v) for taxa in results.values() for v in taxa.values()} == {5}

    def testParseDottedNames(self):
        pattern = os.path.join(self.results_dir, "yn00", "yn00_dotted-*")
        for results_file in glob.glob(pattern):
            results = yn00.read(results_file)
            # Expect seven taxa...
            assert len(results) == 5
            # ...each of which is compared to the other six.
            assert {len(v) for v in results.values()} == {4}
            # ...each of which has five measures.
            assert {len(v) for taxa in results.values() for v in taxa.values()} == {5}
            assert len(results["Homo.sapie"]) == 4
            assert len(results["Homo.sapie"]["Pan.troglo"]) == 5

    def testParseDottedNumNames(self):
        pattern = os.path.join(self.results_dir, "yn00", "yn00_dottednum-*")
        for results_file in glob.glob(pattern):
            results = yn00.read(results_file)
            # Expect seven taxa...
            assert len(results) == 7
            # ...each of which is compared to the other six.
            assert {len(v) for v in results.values()} == {6}
            # ...each of which has five measures.
            assert {len(v) for taxa in results.values() for v in taxa.values()} == {5}


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
