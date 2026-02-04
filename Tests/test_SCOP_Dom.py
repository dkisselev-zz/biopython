# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Unit test for Dom.

This test requires the mini DOM file 'testDom.txt'
"""

import unittest
import pytest

from Bio.SCOP import Dom


class DomTests(unittest.TestCase):
    def setUp(self):
        self.filename = "./SCOP/testDom.txt"

    def testParse(self):
        """Test if all records in a DOM file are being read."""
        count = 0
        with open(self.filename) as f:
            for record in Dom.parse(f):
                count += 1
        assert count == 10

    def testStr(self):
        """Test if we can convert each record to a string correctly."""
        with open(self.filename) as f:
            for line in f:
                record = Dom.Record(line)
                # End of line is platform dependent. Strip it off
                assert str(record).rstrip() == line.rstrip()

    def testError(self):
        """Test if a corrupt record raises the appropriate exception."""
        corruptDom = "49xxx268\tsp\tb.1.2.1\t-\n"
        with pytest.raises(ValueError):
            Dom.Record(corruptDom)

    def testRecord(self):
        """Test one record in detail."""
        recLine = "d7hbib_\t7hbi\tb:\t1.001.001.001.001.001"

        rec = Dom.Record(recLine)
        assert rec.sid == "d7hbib_"
        assert rec.residues.pdbid == "7hbi"
        assert rec.residues.fragments == (("b", "", ""),)
        assert rec.hierarchy == "1.001.001.001.001.001"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
