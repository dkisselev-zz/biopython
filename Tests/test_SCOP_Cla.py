# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# Modifications Copyright 2010 Jeffrey Finkelstein. All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Unit test for Cla."""

import unittest
import pytest

from Bio.SCOP import Cla


class ClaTests(unittest.TestCase):
    def setUp(self):
        self.filename = "./SCOP/dir.cla.scop.txt_test"

    def testParse(self):
        """Test if all records in a CLA file are being read."""
        count = 0
        with open(self.filename) as f:
            records = Cla.parse(f)
            for record in records:
                count += 1
        assert count == 14

    def testStr(self):
        """Test if we can convert each record to a string correctly."""
        with open(self.filename) as f:
            for line in f:
                record = Cla.Record(line)
                # The SCOP Classification file format which can be found at
                # http://scop.mrc-lmb.cam.ac.uk/scop/release-notes.html states
                # that the list of classification hierarchy key-value pairs is
                # unordered, therefore we need only check that they are all
                # there, NOT that they are in the same order.
                # End of line is platform dependent. Strip it off
                expected_hierarchy = line.rstrip().split("\t")[5].split(",")
                expected_hierarchy = dict(
                    pair.split("=") for pair in expected_hierarchy
                )
                actual_hierarchy = str(record).rstrip().split("\t")[5].split(",")
                actual_hierarchy = dict(pair.split("=") for pair in actual_hierarchy)
                assert len(actual_hierarchy) == len(expected_hierarchy)
                for key, actual_value in actual_hierarchy.items():
                    assert actual_value == expected_hierarchy[key]

    def testError(self):
        """Test if a corrupt record raises the appropriate exception."""
        corruptRec = "49268\tsp\tb.1.2.1\t-\n"
        with pytest.raises(ValueError):
            Cla.Record(corruptRec)

    def testRecord(self):
        """Test one record in detail."""
        recLine = "d1dan.1\t1dan\tT:,U:91-106\tb.1.2.1\t21953\tcl=48724,cf=48725,sf=49265,fa=49266,dm=49267,sp=49268,px=21953"

        record = Cla.Record(recLine)
        assert record.sid == "d1dan.1"
        assert record.residues.pdbid == "1dan"
        assert record.residues.fragments == (("T", "", ""), ("U", "91", "106"))
        assert record.sccs == "b.1.2.1"
        assert record.sunid == 21953
        assert record.hierarchy == {
                "cl": 48724,
                "cf": 48725,
                "sf": 49265,
                "fa": 49266,
                "dm": 49267,
                "sp": 49268,
                "px": 21953,
            }

    def testIndex(self):
        """Test CLA file indexing."""
        index = Cla.Index(self.filename)

        assert len(index) == 14
        assert "d4hbia_" in index

        rec = index["d1hbia_"]
        assert rec.sunid == 14996


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
