# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""Testing code for Restriction enzyme classes of Biopython."""

import unittest
import pytest

from Bio import BiopythonWarning
from Bio.Restriction import AanI
from Bio.Restriction import Acc65I
from Bio.Restriction import AllEnzymes
from Bio.Restriction import Analysis
from Bio.Restriction import Asp718I
from Bio.Restriction import BamHI
from Bio.Restriction import BsmBI
from Bio.Restriction import CommOnly
from Bio.Restriction import EarI
from Bio.Restriction import EcoRI
from Bio.Restriction import EcoRV
from Bio.Restriction import FormattedSeq
from Bio.Restriction import KpnI
from Bio.Restriction import McrI
from Bio.Restriction import MluCI
from Bio.Restriction import NdeI
from Bio.Restriction import NonComm
from Bio.Restriction import Restriction
from Bio.Restriction import RestrictionBatch
from Bio.Restriction import SmaI
from Bio.Restriction import SnaI
from Bio.Restriction import SphI
from Bio.Restriction import BsaI
from Bio.Restriction import BsaXI
from Bio.Restriction import BspCNI
from Bio.Seq import MutableSeq
from Bio.Seq import Seq


class SequenceTesting(unittest.TestCase):
    """Tests for dealing with input."""

    def test_sequence_object(self):
        """Test if sequence must be a Seq or MutableSeq object."""
        with pytest.raises(TypeError):
            seq = FormattedSeq("GATC")
        seq = FormattedSeq(Seq("TAGC"))
        seq = FormattedSeq(MutableSeq("AGTC"))
        seq = FormattedSeq(seq)
        with pytest.raises(TypeError):
            EcoRI.search("GATC")
        EcoRI.search(Seq("ATGC"))
        EcoRI.search(MutableSeq("TCAG"))

    def test_non_allowed_characters(self):
        """Test if non-allowed characters raise a TypeError."""
        # Any letter is accepted, even if it's not a nucleotide
        FormattedSeq(Seq("ABCDEFGHIJKLMNOPQRSTUVWXYZ"))
        # Other characters are not accepted
        with pytest.raises(TypeError):
            FormattedSeq(Seq("GATCZE-"))

    def test_formatted_seq(self):
        """Test several methods of FormattedSeq."""
        assert str(FormattedSeq(Seq("GATC"))) == "FormattedSeq(Seq('GATC'), linear=True)"
        assert FormattedSeq(Seq("GATC")) != FormattedSeq(Seq("TAGC"))
        assert FormattedSeq(Seq("TAGC")) != Seq("TAGC")
        assert FormattedSeq(Seq("ATGC")) == FormattedSeq(Seq("ATGC"))
        linear_seq = FormattedSeq(Seq("T"))
        assert linear_seq.is_linear()
        linear_seq.circularise()
        assert not linear_seq.is_linear()
        linear_seq.linearise()
        circular_seq = linear_seq.to_circular()
        assert not circular_seq.is_linear()
        linear_seq = circular_seq.to_linear()
        assert linear_seq.is_linear()


class SimpleEnzyme(unittest.TestCase):
    """Tests for dealing with basic enzymes using the Restriction package."""

    def test_init(self):
        """Check for error during __init__."""
        with pytest.raises(ValueError) as ve:
            Restriction.OneCut("bla-me", (Restriction.RestrictionType,), {})
            assert "hyphen" in str(ve.value)

    def setUp(self):
        """Set up some sequences for later use."""
        base_seq = Seq("AAAA")
        self.ecosite_seq = base_seq + Seq(EcoRI.site) + base_seq
        self.smasite_seq = base_seq + Seq(SmaI.site) + base_seq
        self.kpnsite_seq = base_seq + Seq(KpnI.site) + base_seq

    def test_eco_cutting(self):
        """Test basic cutting with EcoRI (5'overhang)."""
        assert EcoRI.site == "GAATTC"
        assert EcoRI.cut_once()
        assert not EcoRI.is_blunt()
        assert EcoRI.is_5overhang()
        assert not EcoRI.is_3overhang()
        assert EcoRI.overhang() == "5' overhang"
        assert EcoRI.is_defined()
        assert not EcoRI.is_ambiguous()
        assert not EcoRI.is_unknown()
        assert EcoRI.is_palindromic()
        assert EcoRI.is_comm()
        assert "Thermo Fisher Scientific" in EcoRI.supplier_list()
        assert EcoRI.elucidate() == "G^AATT_C"
        assert EcoRI.search(self.ecosite_seq) == [6]
        assert EcoRI.characteristic() == (1, -1, None, None, "GAATTC")

        parts = EcoRI.catalyse(self.ecosite_seq)
        assert len(parts) == 2
        assert str(parts[1]) == "AATTCAAAA"
        parts = EcoRI.catalyze(self.ecosite_seq)
        assert len(parts) == 2

    def test_kpn_cutting(self):
        """Test basic cutting with KpnI (3'overhang)."""
        assert KpnI.is_3overhang()
        assert not KpnI.is_5overhang()
        assert not KpnI.is_blunt()
        assert KpnI.overhang() == "3' overhang"
        parts = KpnI.catalyse(self.kpnsite_seq)
        assert len(parts) == 2
        assert KpnI.catalyse(self.kpnsite_seq) == KpnI.catalyze(self.kpnsite_seq)

    def test_sma_cutting(self):
        """Test basic cutting with SmaI (blunt cutter)."""
        assert SmaI.is_blunt()
        assert not SmaI.is_3overhang()
        assert not SmaI.is_5overhang()
        assert SmaI.overhang() == "blunt"
        parts = SmaI.catalyse(self.smasite_seq)
        assert len(parts) == 2
        assert str(parts[1]) == "GGGAAAA"
        parts = SmaI.catalyze(self.smasite_seq)
        assert len(parts) == 2

    def test_ear_cutting(self):
        """Test basic cutting with EarI (ambiguous overhang)."""
        assert not EarI.is_palindromic()
        assert not EarI.is_defined()
        assert EarI.is_ambiguous()
        assert not EarI.is_unknown()
        assert EarI.elucidate() == "CTCTTCN^NNN_N"

    def test_sna_cutting(self):
        """Test basic cutting with SnaI (unknown)."""
        assert SnaI.elucidate() == "? GTATAC ?"
        assert not SnaI.is_defined()
        assert not SnaI.is_ambiguous()
        assert SnaI.is_unknown()
        assert not SnaI.is_comm()
        assert SnaI.suppliers() is None
        assert SnaI.supplier_list() == []
        with pytest.raises(TypeError):
            SnaI.buffers("no company")

    def test_circular_sequences(self):
        """Deal with cutting circular sequences."""
        parts = EcoRI.catalyse(self.ecosite_seq, linear=False)
        assert len(parts) == 1
        locations = EcoRI.search(parts[0], linear=False)
        assert locations == [1]

        parts = KpnI.catalyse(self.kpnsite_seq, linear=False)
        assert len(parts) == 1
        locations = KpnI.search(parts[0], linear=False)
        assert locations == [1]

        parts = SmaI.catalyse(self.smasite_seq, linear=False)
        assert len(parts) == 1
        locations = SmaI.search(parts[0], linear=False)
        assert locations == [1]

        assert EarI.search(FormattedSeq(Seq("CTCTTCAAAAA")), linear=False) == [8]
        assert SnaI.search(FormattedSeq(Seq("GTATACAAAAA")), linear=False) == [1]

    def test_shortcuts(self):
        """Check if '/' and '//' work as '.search' and '.catalyse'."""
        assert EcoRI / self.ecosite_seq == [6]
        assert self.ecosite_seq / EcoRI == [6]
        assert len(EcoRI // self.ecosite_seq) == 2
        assert len(self.ecosite_seq // EcoRI) == 2

    def test_cutting_border_positions(self):
        """Check if cutting after first and penultimate position works."""
        # Use EarI, cuts as follows: CTCTTCN^NNN_N
        # Only when the cut produces double stranded DNA on both outputs
        # it is returned.
        seq = Seq("CTCTTCA")
        assert EarI.search(seq) == []
        seq += "AAA"
        assert EarI.search(seq) == []
        seq += "A"
        assert EarI.search(seq) == [8]
        # Recognition site on reverse-complement strand
        seq = Seq("AAAAGAAGAG")
        assert EarI.search(seq) == []
        seq = "A" + seq
        assert EarI.search(seq) == [2]

        # Examples from https://github.com/biopython/biopython/issues/4604
        assert BsaI.search(Seq("GGTCTCATAAAA")) == [8]
        assert BsaI.search(Seq("GGTCTCATAAA")) == []
        assert BsaI.search(Seq("GGTCTCGT")) == []
        assert BsaI.search(Seq("GGTCTCGT").reverse_complement()) == []

        assert BsaXI.search(Seq("AAATAAAAAAAAAACAAAAACTCC")) == [5]
        assert BsaXI.search(Seq("AATAAAAAAAAAACAAAAACTCC")) == []

        assert BspCNI.search(Seq("CTCAGAAAAAAAAAT")) == [15]
        assert BspCNI.search(Seq("CTCAGAAAAAAAAA")) == []

    def test_recognition_site_on_both_strands(self):
        """Check if recognition sites on both strands are properly handled."""
        seq = Seq("CTCTTCGAAGAG")
        assert EarI.search(seq) == [3, 8]

    def test_overlapping_cut_sites(self):
        """Check if overlapping recognition sites are properly handled."""
        seq = Seq("CATGCACGCATGCATGCACGC")
        assert SphI.search(seq) == [13, 17]


class EnzymeComparison(unittest.TestCase):
    """Tests for comparing various enzymes."""

    def test_basic_isochizomers(self):
        """Test to be sure isochizomer and neoschizomers are as expected."""
        assert Acc65I.isoschizomers() == [Asp718I, KpnI]
        assert Acc65I.elucidate() == "G^GTAC_C"
        assert Asp718I.elucidate() == "G^GTAC_C"
        assert KpnI.elucidate() == "G_GTAC^C"
        assert Acc65I.is_isoschizomer(KpnI)
        assert not Acc65I.is_equischizomer(KpnI)
        assert Acc65I.is_neoschizomer(KpnI)
        assert Acc65I in Asp718I.equischizomers()
        assert KpnI in Asp718I.neoschizomers()
        assert KpnI in Acc65I.isoschizomers()

    def test_comparisons(self):
        """Test comparison operators between different enzymes."""
        # Comparison of iso- and neoschizomers
        assert Acc65I == Acc65I
        assert Acc65I != KpnI
        assert not Acc65I == Asp718I  # noqa: A500
        # self.assertNotEqual(Acc65I, Asp718I) it doesn't work as expected
        assert not Acc65I != Asp718I  # noqa: A500
        assert Acc65I != EcoRI
        assert Acc65I >> KpnI
        assert not Acc65I >> Asp718I

        # Compare length of recognition sites
        assert not EcoRI >= EcoRV
        assert EcoRV >= EcoRI
        with pytest.raises(NotImplementedError):
            EcoRV >= 3
        assert not EcoRI > EcoRV
        assert EcoRV > EcoRI
        with pytest.raises(NotImplementedError):
            EcoRV > 3
        assert EcoRI <= EcoRV
        assert not EcoRV <= EcoRI
        with pytest.raises(NotImplementedError):
            EcoRV <= 3
        assert EcoRI < EcoRV
        assert not EcoRV < EcoRI
        with pytest.raises(NotImplementedError):
            EcoRV < 3

        # Compare compatible overhangs
        assert Acc65I % Asp718I
        assert Acc65I % Acc65I
        assert not Acc65I % KpnI
        with pytest.raises(TypeError):
            Acc65I % "KpnI"
        assert SmaI % EcoRV
        assert EarI % EarI
        assert EcoRV in SmaI.compatible_end()
        assert Acc65I in Asp718I.compatible_end()


class RestrictionBatchPrintTest(unittest.TestCase):
    """Tests Restriction.Analysis printing functionality."""

    def createAnalysis(self, seq_str, batch_ary):
        """Restriction.Analysis creation helper method."""
        rb = Restriction.RestrictionBatch(batch_ary)
        seq = Seq(seq_str)
        return Restriction.Analysis(rb, seq)

    def assertAnalysisFormat(self, analysis, expected):
        """Test make_format.

        Test that the Restriction.Analysis make_format(print_that) matches
        some string.
        """
        dct = analysis.mapping
        ls, nc = [], []
        for k, v in dct.items():
            if v:
                ls.append((k, v))
            else:
                nc.append(k)
        result = analysis.make_format(ls, "", [], "")
        assert result.replace(" ", "") == expected.replace(" ", "")

    def test_make_format_map1(self):
        """Test that print_as('map'); print_that() correctly wraps round.

        1. With no marker.
        """
        analysis = self.createAnalysis(
            "CCAGTCTATAATTCG"
            + Restriction.BamHI.site
            + "GCGGCATCATACTCGAATATCGCGTGATGATACGTAGTAATTACGCATG",
            ["BamHI"],
        )
        analysis.print_as("map")
        expected = [
            "                17 BamHI",
            "                |                                           ",
            "CCAGTCTATAATTCGGGATCCGCGGCATCATACTCGAATATCGCGTGATGATACGTAGTA",
            "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||",
            "GGTCAGATATTAAGCCCTAGGCGCCGTAGTATGAGCTTATAGCGCACTACTATGCATCAT",
            "1                                                         60",
            "",
            "ATTACGCATG",
            "||||||||||",
            "TAATGCGTAC",
            "61                          70",
            "",
            "",
        ]
        self.assertAnalysisFormat(analysis, "\n".join(expected))

    def test_make_format_map2(self):
        """Test that print_as('map'); print_that() correctly wraps round.

        2. With marker.
        """
        analysis = self.createAnalysis(
            "CCAGTCTATAATTCG"
            + Restriction.BamHI.site
            + "GCGGCATCATACTCGA"
            + Restriction.BamHI.site
            + "ATATCGCGTGATGATA"
            + Restriction.NdeI.site
            + "CGTAGTAATTACGCATG",
            ["NdeI", "EcoRI", "BamHI", "BsmBI"],
        )
        analysis.print_as("map")
        expected = [
            "                17 BamHI",
            "                |                                           ",
            "                |                     39 BamHI",
            "                |                     |                     ",
            "CCAGTCTATAATTCGGGATCCGCGGCATCATACTCGAGGATCCATATCGCGTGATGATAC",
            "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||",
            "GGTCAGATATTAAGCCCTAGGCGCCGTAGTATGAGCTCCTAGGTATAGCGCACTACTATG",
            "1                                                         60",
            "",
            " 62 NdeI",
            " |                                                          ",
            "ATATGCGTAGTAATTACGCATG",
            "||||||||||||||||||||||",
            "TATACGCATCATTAATGCGTAC",
            "61                          82",
            "",
            "",
        ]
        self.assertAnalysisFormat(analysis, "\n".join(expected))

    def test_make_format_map3(self):
        """Test that print_as('map'); print_that() correctly wraps round.

        3. With marker restricted.
        """
        analysis = self.createAnalysis(
            "CCAGTCTATAATTCG"
            + Restriction.BamHI.site
            + "GCGGCATCATACTCGA"
            + Restriction.BamHI.site
            + "ATATCGCGTGATGATA"
            + Restriction.EcoRV.site
            + "CGTAGTAATTACGCATG",
            ["NdeI", "EcoRI", "BamHI", "BsmBI"],
        )
        analysis.print_as("map")
        expected = [
            "                17 BamHI",
            "                |                                           ",
            "                |                     39 BamHI",
            "                |                     |                     ",
            "CCAGTCTATAATTCGGGATCCGCGGCATCATACTCGAGGATCCATATCGCGTGATGATAG",
            "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||",
            "GGTCAGATATTAAGCCCTAGGCGCCGTAGTATGAGCTCCTAGGTATAGCGCACTACTATC",
            "1                                                         60",
            "",
            "ATATCCGTAGTAATTACGCATG",
            "||||||||||||||||||||||",
            "TATAGGCATCATTAATGCGTAC",
            "61                          82",
            "",
            "",
        ]
        self.assertAnalysisFormat(analysis, "\n".join(expected))

    def test_change(self):
        """Test that change() changes something."""
        seq = Seq(
            "CCAGTCTATAATTCG"
            + BamHI.site
            + "GCGGCATCATACTCGA"
            + BamHI.site
            + "ATATCGCGTGATGATA"
            + EcoRV.site
            + "CGTAGTAATTACGCATG"
        )
        batch = NdeI + EcoRI + BamHI + BsmBI
        analysis = Analysis(batch, seq)
        assert analysis.full()[BamHI] == [17, 39]
        batch = NdeI + EcoRI + BsmBI
        seq += NdeI.site
        analysis.change(sequence=seq)
        analysis.change(rb=batch)
        assert len(analysis.full()) == 3
        assert analysis.full()[NdeI] == [85]
        with pytest.raises(AttributeError):
            analysis.change(**{"NameWidth": 3, "KonsoleWidth": 40})  # Console


class RestrictionBatches(unittest.TestCase):
    """Tests for dealing with batches of restriction enzymes."""

    def test_creating_batch(self):
        """Creating and modifying a restriction batch."""
        batch = RestrictionBatch()
        assert batch.suppl_codes()["N"] == "New England Biolabs"
        assert batch.is_restriction(EcoRI)
        batch = RestrictionBatch([EcoRI])
        batch.add(KpnI)
        batch += EcoRV
        assert len(batch) == 3
        assert batch.elements() == ["EcoRI", "EcoRV", "KpnI"]
        # Problem with Python 3, as sequence of list may be different:
        # self.assertEqual(batch.as_string(), ['EcoRI', 'KpnI', 'EcoRV'])
        assert "EcoRI" in batch.as_string()

        # The usual way to test batch membership
        assert EcoRV in batch
        assert EcoRI in batch
        assert KpnI in batch
        assert SmaI not in batch
        # Syntax sugar for the above
        assert "EcoRV" in batch
        assert "SmaI" not in batch

        batch.get(EcoRV)
        with pytest.raises(ValueError):
            batch.get(SmaI)
        batch.get(SmaI, add=True)
        assert len(batch) == 4
        batch.remove(SmaI)
        batch.remove(EcoRV)
        assert len(batch) == 2

        assert EcoRV not in batch
        assert "EcoRV" not in batch

        # Creating a batch by addition of restriction enzymes
        new_batch = EcoRI + KpnI
        assert batch == new_batch
        # or by addition of a batch with an enzyme
        another_new_batch = new_batch + EcoRV
        new_batch += EcoRV
        assert another_new_batch == new_batch
        with pytest.raises(TypeError):
            EcoRI.__add__(1)

        # Create a batch with suppliers and other supplier related methods
        # These tests may be 'update sensitive' since company names and
        # products may change often...
        batch = RestrictionBatch((), ("S"))  # Sigma
        assert batch.current_suppliers() == ["Sigma Chemical Corporation"]
        assert EcoRI in batch
        assert AanI not in batch
        batch.add_supplier("B")  # Thermo Fisher Scientific
        assert AanI in batch

    def test_batch_analysis(self):
        """Sequence analysis with a restriction batch."""
        seq = Seq("AAAA" + EcoRV.site + "AAAA" + EcoRI.site + "AAAA")
        batch = RestrictionBatch([EcoRV, EcoRI])

        hits = batch.search(seq)
        assert hits[EcoRV] == [8]
        assert hits[EcoRI] == [16]

    def test_premade_batches(self):
        """Test content of premade batches CommOnly, NoComm, AllEnzymes."""
        assert len(AllEnzymes) == (len(CommOnly) + len(NonComm))
        assert len(AllEnzymes) > len(CommOnly) > len(NonComm)

    def test_search_premade_batches(self):
        """Test search with pre-made batches CommOnly, NoComm, AllEnzymes."""
        seq = Seq("ACCCGAATTCAAAACTGACTGATCGATCGTCGACTG")
        search = AllEnzymes.search(seq)
        assert search[MluCI] == [6]
        # Check if '/' operator works as 'search':
        search = CommOnly / seq
        assert search[MluCI] == [6]
        # Also in reverse order:
        search = seq / NonComm
        assert search[McrI] == [28]

    def test_analysis_restrictions(self):
        """Test Fancier restriction analysis."""
        new_seq = Seq("TTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAA")
        rb = RestrictionBatch([EcoRI, KpnI, EcoRV])
        ana = Analysis(rb, new_seq, linear=False)
        # Output only the result for enzymes which cut blunt:
        assert ana.blunt() == {EcoRV: []}
        assert ana.full() == {KpnI: [], EcoRV: [], EcoRI: [33]}
        # Output only the result for enzymes which have a site:
        assert ana.with_sites() == {EcoRI: [33]}
        # Output only the enzymes which have no site:
        assert ana.without_site() == {KpnI: [], EcoRV: []}
        assert ana.with_site_size([32]) == {}
        # Output only enzymes which produce 5' overhangs
        assert ana.overhang5() == {EcoRI: [33]}
        # Output only enzymes which produce 3' overhangs
        assert ana.overhang3() == {KpnI: []}
        # Output only enzymes which produce defined ends
        assert ana.defined() == {KpnI: [], EcoRV: [], EcoRI: [33]}
        # Output only enzymes hich cut N times
        assert ana.with_N_sites(2) == {}
        # The enzymes which cut between position x and y:
        with pytest.raises(TypeError):
            ana.only_between("t", 20)
        with pytest.raises(TypeError):
            ana.only_between(1, "t")
        assert ana.only_between(1, 20) == {}
        assert ana.only_between(20, 34) == {EcoRI: [33]}
        # Mix start/end order:
        assert ana.only_between(34, 20) == {EcoRI: [33]}
        assert ana.only_outside(20, 34) == {}
        with pytest.warns(BiopythonWarning):
            ana.with_name(["fake"])
        assert ana.with_name([EcoRI]) == {EcoRI: [33]}
        assert (ana._boundaries(1, 20)[:2]) == (1, 20)
        # Reverse order:
        assert (ana._boundaries(20, 1)[:2]) == (1, 20)
        # Fix negative start:
        assert (ana._boundaries(-1, 20)[:2]) == (20, 33)
        # Fix negative end:
        assert (ana._boundaries(1, -1)[:2]) == (1, 33)
        # Sites in- and outside of boundaries
        new_seq = Seq("GAATTCAAAAAAGAATTC")
        rb = RestrictionBatch([EcoRI])
        ana = Analysis(rb, new_seq)
        # Cut at least inside
        assert ana.between(1, 7) == {EcoRI: [2, 14]}
        # Cut at least inside and report only inside site
        assert ana.show_only_between(1, 7) == {EcoRI: [2]}
        # Cut at least outside
        assert ana.outside(1, 7) == {EcoRI: [2, 14]}
        # Don't cut within
        assert ana.do_not_cut(7, 12) == {EcoRI: [2, 14]}


class TestPrintOutputs(unittest.TestCase):
    """Class to test various print outputs."""

    import sys
    from io import StringIO

    def test_supplier(self):
        """Test output of supplier list for different enzyme types."""
        out = self.StringIO()
        self.sys.stdout = out
        EcoRI.suppliers()
        assert "Thermo Fisher Scientific" in out.getvalue()
        assert SnaI.suppliers() is None
        EcoRI.all_suppliers()  # Independent of enzyme, list of all suppliers
        assert "Agilent Technologies" in out.getvalue()
        batch = EcoRI + SnaI
        batch.show_codes()
        assert "N = New England Biolabs" in out.getvalue()
        self.sys.stdout = self.sys.__stdout__

    def test_print_that(self):
        """Test print_that function."""
        out = self.StringIO()
        self.sys.stdout = out
        my_batch = EcoRI + SmaI + KpnI
        my_seq = Seq("GAATTCCCGGGATATA")  # EcoRI and SmaI sites
        analysis = Analysis(my_batch, my_seq)
        analysis.print_that(None, title="My sequence\n\n", s1="Non Cutters\n\n")
        assert "My sequence" in out.getvalue()
        assert "Non Cutters" in out.getvalue()
        assert "2." in out.getvalue()
        self.sys.stdout = self.sys.__stdout__

    def test_str_method(self):
        """Test __str__ and __repr__ outputs."""
        batch = EcoRI + SmaI + KpnI
        assert str(batch) == "EcoRI+KpnI+SmaI"
        batch += Asp718I
        batch += SnaI
        assert str(batch) == "Asp718I+EcoRI...SmaI+SnaI"
        assert repr(batch) == "RestrictionBatch(['Asp718I', 'EcoRI', 'KpnI', 'SmaI', 'SnaI'])"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
