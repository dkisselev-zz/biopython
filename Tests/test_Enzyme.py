# Copyright 1999 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""Tests for Enzyme module."""

import os
import unittest
import pytest
from io import StringIO

from Bio.ExPASy import Enzyme


class TestEnzyme(unittest.TestCase):
    def test_parse_zero(self):
        handle = StringIO("")
        records = list(Enzyme.parse(handle))
        assert len(records) == 0

    def test_parse_one(self):
        """Check parse function with one record."""
        with open("Enzymes/lipoprotein.txt") as handle:
            records = list(Enzyme.parse(handle))
        assert len(records) == 1
        assert records[0]["ID"] == "3.1.1.34"

    def test_parse_many(self):
        """Check parse function with multiple records."""
        data = ""
        for filename in [
            "Enzymes/lipoprotein.txt",
            "Enzymes/proline.txt",
            "Enzymes/valine.txt",
        ]:
            with open(filename) as handle:
                data += handle.read()
        handle = StringIO(data)
        records = list(Enzyme.parse(handle))
        assert len(records) == 3
        assert records[0]["ID"] == "3.1.1.34"
        assert records[1]["ID"] == "5.1.1.4"
        assert records[2]["ID"] == "4.1.1.14"

    def test_lipoprotein(self):
        """Parsing ENZYME record for lipoprotein lipase (3.1.1.34)."""
        filename = os.path.join("Enzymes", "lipoprotein.txt")
        with open(filename) as handle:
            record = Enzyme.read(handle)
        assert record["ID"] == "3.1.1.34"
        assert record["DE"] == "Lipoprotein lipase."
        assert len(record["AN"]) == 3
        assert record["AN"][0] == "Clearing factor lipase."
        assert record["AN"][1] == "Diacylglycerol lipase."
        assert record["AN"][2] == "Diglyceride lipase."
        assert record["CA"] == "Triacylglycerol + H(2)O = diacylglycerol + a carboxylate."
        assert record["CC"][0] == "Hydrolyzes triacylglycerols in chylomicrons and very low-density lipoproteins (VLDL)."
        assert record["CC"][1] == "Also hydrolyzes diacylglycerol."
        assert record["PR"] == ["PDOC00110"]
        assert record["DR"][0] == ["P11151", "LIPL_BOVIN"]
        assert record["DR"][1] == ["P11153", "LIPL_CAVPO"]
        assert record["DR"][2] == ["P11602", "LIPL_CHICK"]
        assert record["DR"][3] == ["P55031", "LIPL_FELCA"]
        assert record["DR"][4] == ["P06858", "LIPL_HUMAN"]
        assert record["DR"][5] == ["P11152", "LIPL_MOUSE"]
        assert record["DR"][6] == ["O46647", "LIPL_MUSVI"]
        assert record["DR"][7] == ["P49060", "LIPL_PAPAN"]
        assert record["DR"][8] == ["P49923", "LIPL_PIG"]
        assert record["DR"][9] == ["Q06000", "LIPL_RAT"]
        assert record["DR"][10] == ["Q29524", "LIPL_SHEEP"]
        assert str(record).startswith("ID: 3.1.1.34\nDE: Lipoprotein lipase.\n"), f"Did not expect:\n{record}"

    def test_proline(self):
        """Parsing ENZYME record for proline racemase (5.1.1.4)."""
        filename = os.path.join("Enzymes", "proline.txt")
        with open(filename) as handle:
            record = Enzyme.read(handle)
        assert record["ID"] == "5.1.1.4"
        assert record["DE"] == "Proline racemase."
        assert record["CA"] == "L-proline = D-proline."
        assert len(record["DR"]) == 9
        assert record["DR"][0] == ["Q17ZY4", "PRAC_CLOD6"]
        assert record["DR"][1] == ["A8DEZ8", "PRAC_CLODI"]
        assert record["DR"][2] == ["Q4DA80", "PRCMA_TRYCR"]
        assert record["DR"][3] == ["Q868H8", "PRCMB_TRYCR"]
        assert record["DR"][4] == ["Q3SX04", "PRCM_BOVIN"]
        assert record["DR"][5] == ["Q96EM0", "PRCM_HUMAN"]
        assert record["DR"][6] == ["Q9CXA2", "PRCM_MOUSE"]
        assert record["DR"][7] == ["Q5RC28", "PRCM_PONAB"]
        assert record["DR"][8] == ["Q66II5", "PRCM_XENTR"]
        assert str(record).startswith("ID: 5.1.1.4\nDE: Proline racemase.\n"), f"Did not expect:\n{record}"

    def test_valine(self):
        """Parsing ENZYME record for valine decarboxylase (4.1.1.14)."""
        filename = os.path.join("Enzymes", "valine.txt")
        with open(filename) as handle:
            record = Enzyme.read(handle)
        assert record["ID"] == "4.1.1.14"
        assert record["DE"] == "Valine decarboxylase."
        assert record["CA"] == "L-valine = 2-methylpropanamine + CO(2)."
        assert record["CF"] == "Pyridoxal 5'-phosphate."
        assert record["CC"] == ["Also acts on L-leucine."]
        assert len(record["DR"]) == 0
        assert str(record).startswith("ID: 4.1.1.14\nDE: Valine decarboxylase.\n"), f"Did not expect:\n{record}"

    def test_lactate(self):
        """Parsing ENZYME record for lactate racemase (5.1.2.1)."""
        filename = os.path.join("Enzymes", "lactate.txt")
        with open(filename) as handle:
            record = Enzyme.read(handle)
        assert record["ID"] == "5.1.2.1"
        assert record["DE"] == "Lactate racemase."
        assert len(record["AN"]) == 3
        assert record["AN"][0] == "Hydroxyacid racemase."
        assert record["AN"][1] == "Lactic acid racemase."
        assert record["AN"][2] == "Lacticoracemase."
        assert record["CA"] == "(S)-lactate = (R)-lactate."
        assert len(record["DR"]) == 0
        assert str(record).startswith("ID: 5.1.2.1\nDE: Lactate racemase.\n"), f"Did not expect:\n{record}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
