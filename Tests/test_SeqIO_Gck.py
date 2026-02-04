# Copyright 2019 Damien Goutte-Gattat.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Tests for the SeqIO Gck module."""

import unittest
import pytest
from io import BytesIO

from Bio import SeqIO


class TestGckWithArtificialData(unittest.TestCase):
    def setUp(self):
        with open("Gck/artificial.gck", "rb") as f:
            self.buffer = f.read()

    def test_read(self):
        """Read an artificial sample file."""
        h = BytesIO(self.buffer)
        record = SeqIO.read(h, "gck")
        assert "ACGTACGTACGT" == record.seq
        assert "Sample construct" == record.description
        assert "linear" == record.annotations["topology"]
        assert 2 == len(record.features)

        assert 2 == record.features[0].location.start
        assert 6 == record.features[0].location.end
        assert 1 == record.features[0].location.strand
        assert "misc_feature" == record.features[0].type
        assert "FeatureA" == record.features[0].qualifiers["label"][0]

        assert 7 == record.features[1].location.start
        assert 11 == record.features[1].location.end
        assert -1 == record.features[1].location.strand
        assert "CDS" == record.features[1].type
        assert "FeatureB" == record.features[1].qualifiers["label"][0]

        h.close()

    def munge_buffer(self, position, value):
        mod_buffer = bytearray(self.buffer)
        if isinstance(value, list):
            mod_buffer[position : position + len(value) - 1] = value
        else:
            mod_buffer[position] = value
        return BytesIO(mod_buffer)

    def test_conflicting_lengths(self):
        """Read a file with incorrect length."""
        # Change the sequence length as indicated in the sequence packet
        h = self.munge_buffer(0x1C, [0x00, 0x00, 0x20, 0x15])
        with pytest.raises(ValueError, match="Conflicting sequence length values"):
            SeqIO.read(h, "gck")
        h.close()

        # Change the sequence length as indicated in the features packet
        h = self.munge_buffer(0x36, [0x00, 0x00, 0x20, 0x15])
        with pytest.raises(ValueError, match="Conflicting sequence length values"):
            SeqIO.read(h, "gck")
        h.close()

        # Change the number of features
        h = self.munge_buffer(0x3B, 0x30)
        with pytest.raises(ValueError, match="Features packet size inconsistent with number of features"):
            SeqIO.read(h, "gck")
        h.close()

        # Change the number of restriction sites
        h = self.munge_buffer(0x137, 0x30)
        with pytest.raises(ValueError, match="Sites packet size inconsistent with number of sites"):
            SeqIO.read(h, "gck")
        h.close()


class TestGckWithImproperHeader(unittest.TestCase):
    def test_read(self):
        """Read a file with an incomplete header."""
        stream = BytesIO(b"tiny")
        with pytest.raises(ValueError, match="Improper header, cannot read 24 bytes from stream"):
            SeqIO.read(stream, "gck")
        stream.close()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
