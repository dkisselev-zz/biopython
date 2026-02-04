# Copyright 2017 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Testing Bio.ExPASy online code."""

import unittest
import pytest


# We want to test these:
from Bio import ExPASy

# In order to check any records returned
from Bio.ExPASy import Prodoc
from Bio.ExPASy import Prosite
from Bio.ExPASy import ScanProsite

pytestmark = pytest.mark.online



class ExPASyOnlineTests(unittest.TestCase):
    """Test ExPASy online resources."""

    def test_prosite_raw(self):
        with ExPASy.get_prosite_raw("PS00001") as handle:
            record = Prosite.read(handle)
        assert record.accession == "PS00001"
        assert record.name == "ASN_GLYCOSYLATION"

    def test_prodoc_raw(self):
        with ExPASy.get_prosite_raw("PDOC00001") as handle:
            record = Prodoc.read(handle)
        assert record.accession == "PDOC00001"

    def test_prosite_html(self):
        with ExPASy.get_prosite_entry("PS00001") as handle:
            html = handle.read()
        assert handle.url == "https://prosite.expasy.org/cgi-bin/prosite/get-prosite-entry?PS00001"
        assert "<title>PROSITE - PS00001</title>" in html

    def test_prodoc_html(self):
        with ExPASy.get_prodoc_entry("PDOC00001") as handle:
            html = handle.read()
        assert handle.url == "https://prosite.expasy.org/cgi-bin/prosite/get-prodoc-entry?PDOC00001"
        assert "{PS00001; ASN_GLYCOSYLATION}" in html

    def test_scanprosite_swissprot(self):
        pattern = "P-x(2)-G-E-S-G(2)-[AS]"
        result = ScanProsite.scan(
            sig=pattern, mirror=ScanProsite.PROSITE_URL, output="xml", db="sp"
        )
        sequences = ScanProsite.read(result)
        assert len(sequences) > 0

    def test_scanprosite_trembl(self):
        pattern = "P-x(2)-G-E-S-G(2)-[AS]"
        result = ScanProsite.scan(
            sig=pattern, mirror=ScanProsite.PROSITE_URL, output="xml", db="tr"
        )
        sequences = ScanProsite.read(result)
        assert len(sequences) > 0

    def test_scanprosite_pdb(self):
        pattern = "P-x(2)-G-E-S-G(2)-[AS]"
        result = ScanProsite.scan(
            sig=pattern, mirror=ScanProsite.PROSITE_URL, output="xml", db="pdb"
        )
        sequences = ScanProsite.read(result)
        assert len(sequences) > 0

    def test_scanprosite_output_not_implemented(self):
        pattern = "P-x(2)-G-E-S-G(2)-[AS]"
        with pytest.raises(NotImplementedError):
            ScanProsite.scan(
                sig=pattern, mirror=ScanProsite.PROSITE_URL, output="txt", db="sp"
            )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
