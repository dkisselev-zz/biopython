# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# Modifications Copyright 2010 Jeffrey Finkelstein. All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Unit test for Scop."""

import unittest
import pytest
from io import StringIO

from Bio.SCOP import cmp_sccs
from Bio.SCOP import parse_domain
from Bio.SCOP import Scop


class ScopTests(unittest.TestCase):
    def _compare_cla_lines(self, cla_line_1, cla_line_2):
        """Compare the two specified Cla lines for equality.

        The order of the key-value pairs in the sixth field of the lines does
        not matter. For more information, see
        http://scop.mrc-lmb.cam.ac.uk/scop/release-notes.html.
        """
        fields1 = cla_line_1.rstrip().split("\t")
        fields2 = cla_line_2.rstrip().split("\t")
        # compare the first five fields in a Cla line, which should be exactly
        # the same
        assert fields1[:5] == fields2[:5]
        # compare the hierarchy key-value pairs, which are unordered
        assert sorted(fields1[5].split(",")) == sorted(fields2[5].split(","))

    def testParse(self):
        with open("./SCOP/dir.cla.scop.txt_test") as f:
            cla = f.read()
        with open("./SCOP/dir.des.scop.txt_test") as f:
            des = f.read()
        with open("./SCOP/dir.hie.scop.txt_test") as f:
            hie = f.read()

        scop = Scop(StringIO(cla), StringIO(des), StringIO(hie))

        cla_out = StringIO()
        scop.write_cla(cla_out)
        lines = zip(cla.rstrip().split("\n"), cla_out.getvalue().rstrip().split("\n"))
        for expected_line, line in lines:
            self._compare_cla_lines(expected_line, line)

        des_out = StringIO()
        scop.write_des(des_out)
        assert des_out.getvalue() == des

        hie_out = StringIO()
        scop.write_hie(hie_out)
        assert hie_out.getvalue() == hie

        domain = scop.getDomainBySid("d1hbia_")
        assert domain.sunid == 14996

        domains = scop.getDomains()
        assert len(domains) == 14
        assert domains[4].sunid == 14988

        dom = scop.getNodeBySunid(-111)
        assert dom is None
        dom = scop.getDomainBySid("no such domain")
        assert dom is None

    def testSccsOrder(self):
        assert cmp_sccs("a.1.1.1", "a.1.1.1") == 0
        assert cmp_sccs("a.1.1.2", "a.1.1.1") == 1
        assert cmp_sccs("a.1.1.2", "a.1.1.11") == -1
        assert cmp_sccs("a.1.2.2", "a.1.1.11") == 1
        assert cmp_sccs("a.1.2.2", "a.5.1.11") == -1
        assert cmp_sccs("b.1.2.2", "a.5.1.11") == 1
        assert cmp_sccs("b.1.2.2", "b.1.2") == 1

    def testParseDomain(self):
        s = ">d1tpt_1 a.46.2.1 (1-70) Thymidine phosphorylase {Escherichia coli}"
        dom = parse_domain(s)

        assert dom.sid == "d1tpt_1"
        assert dom.sccs == "a.46.2.1"
        assert dom.residues.pdbid == "1tpt"
        assert dom.description == "Thymidine phosphorylase {Escherichia coli}"

        s2 = "d1tpt_1 a.46.2.1 (1tpt 1-70) Thymidine phosphorylase {E. coli}"
        assert s2 == str(parse_domain(s2))

        # Genetic domains (See Astral release notes)
        s3 = "g1cph.1 g.1.1.1 (1cph B:,A:) Insulin {Cow (Bos taurus)}"
        assert s3 == str(parse_domain(s3))

        s4 = "e1cph.1a g.1.1.1 (1cph A:) Insulin {Cow (Bos taurus)}"
        assert s4 == str(parse_domain(s4))

        # Raw Astral header
        s5 = ">e1cph.1a g.1.1.1 (A:) Insulin {Cow (Bos taurus)}"
        assert s4 == str(parse_domain(s5))

        with pytest.raises(ValueError):
            parse_domain("Totally wrong")

    def testConstructFromDirectory(self):
        scop = Scop(dir_path="SCOP", version="test")
        assert isinstance(scop, Scop)

        domain = scop.getDomainBySid("d1hbia_")
        assert domain.sunid == 14996

    def testGetAscendent(self):
        scop = Scop(dir_path="SCOP", version="test")
        domain = scop.getDomainBySid("d1hbia_")

        # get the fold
        fold = domain.getAscendent("cf")
        assert fold.sunid == 46457

        # get the superfamily
        sf = domain.getAscendent("superfamily")
        assert sf.sunid == 46458

        # px has no px ascendent
        px = domain.getAscendent("px")
        assert px is None

        # an sf has no px ascendent
        px2 = sf.getAscendent("px")
        assert px2 is None

    def test_get_descendents(self):
        """Test getDescendents method."""
        scop = Scop(dir_path="SCOP", version="test")
        fold = scop.getNodeBySunid(46457)

        # get px descendents
        domains = fold.getDescendents("px")
        assert len(domains) == 14
        for d in domains:
            assert d.type == "px"

        sfs = fold.getDescendents("superfamily")
        assert len(sfs) == 1
        for d in sfs:
            assert d.type == "sf"

        # cl has no cl descendent
        cl = fold.getDescendents("cl")
        assert cl == []


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
