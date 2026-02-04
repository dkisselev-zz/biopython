# Copyright 2009 by James Casbon.  All rights reserved.
# Revisions copyright 2009-2010 by Michiel de Hoon. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for parsing Compass output."""

import os
import unittest
import pytest

from Bio import Compass


class CompassTest(unittest.TestCase):
    def setUp(self):
        file_dir = os.path.join("Compass")
        self.test_files = [
            os.path.join(file_dir, "comtest1"),
            os.path.join(file_dir, "comtest2"),
        ]

    def testCompassScanAndConsume(self):
        with open(self.test_files[0]) as handle:
            com_record = Compass.read(handle)

        assert "60456.blo.gz.aln" == com_record.query
        assert "60456.blo.gz.aln" == com_record.hit
        assert 0.5 == pytest.approx(com_record.gap_threshold, abs=5e-8)

        assert 388 == com_record.query_length
        assert 386 == com_record.query_filtered_length
        assert 388 == com_record.hit_length
        assert 386 == com_record.hit_filtered_length

        assert 399 == com_record.query_nseqs
        assert 12.972 == pytest.approx(com_record.query_neffseqs, abs=5e-8)
        assert 399 == com_record.hit_nseqs
        assert 12.972 == pytest.approx(com_record.hit_neffseqs, abs=5e-8)

        assert 2759 == com_record.sw_score
        assert 0.0 == pytest.approx(com_record.evalue, abs=5e-8)

    def testCompassParser(self):
        with open(self.test_files[0]) as handle:
            com_record = Compass.read(handle)

        assert "60456.blo.gz.aln" == com_record.query

    def testCompassIteratorEasy(self):
        with open(self.test_files[0]) as handle:
            records = Compass.parse(handle)
            com_record = next(records)
        assert "60456.blo.gz.aln" == com_record.query
        with pytest.raises(StopIteration):
            next(records)

    def testCompassIteratorHard(self):
        with open(self.test_files[1]) as handle:
            records = Compass.parse(handle)

            com_record = next(records)
            assert "allscop//14982.blo.gz.aln" == com_record.hit
            assert 1.01e03 == pytest.approx(com_record.evalue, abs=5e-8)

            com_record = next(records)
            assert "allscop//14983.blo.gz.aln" == com_record.hit
            assert 1.01e03 == pytest.approx(com_record.evalue, abs=5e-8)

            com_record = next(records)
            assert "allscop//14984.blo.gz.aln" == com_record.hit
            assert 5.75e02 == pytest.approx(com_record.evalue, abs=5e-8)

    def testAlignmentParsingOne(self):
        with open(self.test_files[1]) as handle:
            records = Compass.parse(handle)

            com_record = next(records)
            assert 178 == com_record.query_start
            assert "KKDLEEIAD" == com_record.query_aln
            assert 9 == com_record.hit_start
            assert "QAAVQAVTA" == com_record.hit_aln
            assert "++ ++++++" == com_record.positives

            com_record = next(records)
            com_record = next(records)
            assert 371 == com_record.query_start
            assert "LEEAMDRMER~~~V" == com_record.query_aln
            assert 76 == com_record.hit_start
            assert "LQNFIDQLDNpddL" == com_record.hit_aln
            assert "+ ++++ + +   +" == com_record.positives

    def testAlignmentParsingTwo(self):
        with open(self.test_files[0]) as handle:
            records = Compass.parse(handle)
            com_record = next(records)
        assert 2 == com_record.query_start
        assert 2 == com_record.hit_start
        assert "LKERKL" == com_record.hit_aln[-6:]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
