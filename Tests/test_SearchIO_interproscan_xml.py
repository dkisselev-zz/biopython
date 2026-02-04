# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO InterproscanIO parsers."""

import os
import unittest
import pytest

from Bio.SearchIO import parse

# test case files are in the Blast directory
TEST_DIR = "InterProScan"
FMT = "interproscan-xml"


def get_file(filename):
    """Return the path of a test file."""
    return os.path.join(TEST_DIR, filename)


class InterproscanXmlCases(unittest.TestCase):
    def test_xml_001(self):
        xml_file = get_file("test_001.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = next(qresults)
        counter += 1

        assert "5.26-65.0" == qresult.version

        # test parsed values of qresult
        assert "AT5G23090.4" == qresult.id
        assert ("pacid=19665592 transcript=AT5G23090.4 locus=AT5G23090 "
            "ID=AT5G23090.4.TAIR10 annot-version=TAIR10" == qresult.description)
        assert 4 == len(qresult)

        hit = qresult[0]
        assert "PF00808" == hit.id
        assert "Histone-like transcription factor (CBF/NF-Y) and archaeal histone" == hit.description
        assert "PFAM" == hit.attributes["Target"]
        assert "31.0" == hit.attributes["Target version"]
        assert "hmmer3" == hit.attributes["Hit type"]
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert 76.7 == hsp.bitscore
        assert 1.1e-21 == hsp.evalue
        assert 13 == hsp.query_start
        assert 79 == hsp.query_end
        assert 0 == hsp.hit_start
        assert 65 == hsp.hit_end
        assert 66 == hsp.aln_span
        assert ("MDPMDIVGKSKEDASLPKATMTKIIKEMLPPDVRVARDAQDLLIECCVEFINLVSSESNDVCNKEDKRTIAPEH"
            "VLKALQVLGFGEYIEEVYAAYEQHKYETMDTQRSVKWNPGAQMTEEEAAAEQQRMFAEARARMNGGVSVPQPEH"
            "PETDQRSPQS" == hsp.query.seq)

        # parse last hit
        hit = qresult[-1]
        assert "SSF47113" == hit.id
        assert 1 == len(hit)
        assert "IPR:IPR009072" == hit.dbxrefs[0]
        assert "GO:0046982" == hit.dbxrefs[1]

        hsp = hit.hsps[0]
        assert 11 == hsp.query_start
        assert 141 == hsp.query_end


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
