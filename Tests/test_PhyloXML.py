# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the PhyloXML and PhyloXMLIO modules."""

import os
import platform  # for Windows hack, see issue #3944
import sys  # for Windows hack
import tempfile
import unittest
import pytest
from itertools import chain

from Bio.Align import Alignment
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo import PhyloXML as PX
from Bio.Phylo import PhyloXMLIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Example PhyloXML files
EX_APAF = "PhyloXML/apaf.xml"
EX_BCL2 = "PhyloXML/bcl_2.xml"
EX_MADE = "PhyloXML/made_up.xml"
EX_PHYLO = "PhyloXML/phyloxml_examples.xml"
EX_DOLLO = "PhyloXML/o_tol_332_d_dollo.xml"

# Temporary file name for Writer tests below
DUMMY = tempfile.NamedTemporaryFile(delete=False).name


# ---------------------------------------------------------
# Parser tests


def _test_read_factory(source, count):
    """Generate a test method for read()ing the given source.

    The generated function reads an example file to produce a phyloXML object,
    then tests for existence of the root node, and counts the number of
    phylogenies under the root.
    """
    fname = os.path.basename(source)

    def test_read(self):
        phx = PhyloXMLIO.read(source)
        assert phx
        assert len(phx) == count[0]
        assert len(phx.other) == count[1]

    test_read.__doc__ = f"Read {fname} to produce a phyloXML object."
    return test_read


def _test_parse_factory(source, count):
    """Generate a test method for parse()ing the given source.

    The generated function extracts each phylogenetic tree using the parse()
    function and counts the total number of trees extracted.
    """
    fname = os.path.basename(source)

    def test_parse(self):
        trees = PhyloXMLIO.parse(source)
        assert len(list(trees)) == count

    test_parse.__doc__ = f"Parse the phylogenies in {fname}."
    return test_parse


def _test_shape_factory(source, shapes):
    """Generate a test method for checking tree shapes.

    Counts the branches at each level of branching in a phylogenetic tree, 3
    clades deep.
    """
    fname = os.path.basename(source)

    def test_shape(self):
        trees = PhyloXMLIO.parse(source)
        for tree, shape_expect in zip(trees, shapes):
            assert len(tree.clade) == len(shape_expect)
            for clade, sub_expect in zip(tree.clade, shape_expect):
                assert len(clade) == sub_expect[0]
                for subclade, len_expect in zip(clade, sub_expect[1]):
                    assert len(subclade) == len_expect

    test_shape.__doc__ = f"Check the branching structure of {fname}."
    return test_shape


class ParseTests(unittest.TestCase):
    """Tests for proper parsing of example phyloXML files."""

    test_read_apaf = _test_read_factory(EX_APAF, (1, 0))
    test_read_bcl2 = _test_read_factory(EX_BCL2, (1, 0))
    test_read_made = _test_read_factory(EX_MADE, (6, 0))
    test_read_phylo = _test_read_factory(EX_PHYLO, (14, 1))
    test_read_dollo = _test_read_factory(EX_DOLLO, (1, 0))

    test_parse_apaf = _test_parse_factory(EX_APAF, 1)
    test_parse_bcl2 = _test_parse_factory(EX_BCL2, 1)
    test_parse_made = _test_parse_factory(EX_MADE, 6)
    test_parse_phylo = _test_parse_factory(EX_PHYLO, 14)
    test_parse_dollo = _test_parse_factory(EX_DOLLO, 1)

    # lvl-2 clades, sub-clade counts, lvl-3 clades
    test_shape_apaf = _test_shape_factory(EX_APAF, (((2, (2, 2)), (2, (2, 2))),))
    test_shape_bcl2 = _test_shape_factory(EX_BCL2, (((2, (2, 2)), (2, (2, 2))),))
    test_shape_phylo = _test_shape_factory(
        EX_PHYLO,
        (
            ((2, (0, 0)), (0, ())),
            ((2, (0, 0)), (0, ())),
            ((2, (0, 0)), (0, ())),
            ((2, (0, 0)), (0, ())),
            ((2, (0, 0)), (0, ())),
            ((2, (0, 0)), (0, ())),
            ((2, (0, 0)), (0, ())),
            ((2, (0, 0)), (0, ())),
            ((2, (0, 0)), (0, ())),
            ((0, ()), (2, (0, 0))),
            ((3, (0, 0, 0)), (0, ())),
            ((2, (0, 0)), (0, ())),
            ((2, (0, 0)), (0, ())),
        ),
    )
    test_shape_dollo = _test_shape_factory(EX_DOLLO, (((2, (2, 2)), (2, (2, 2))),))


class TreeTests(unittest.TestCase):
    """Tests for instantiation and attributes of each complex type."""

    # ENH: also test check_str() regexps wherever they're used

    def test_Phyloxml(self):
        """Instantiation of Phyloxml objects."""
        phx = PhyloXMLIO.read(EX_PHYLO)
        assert isinstance(phx, PX.Phyloxml)
        for tree in phx:
            assert isinstance(tree, PX.Phylogeny)
        for otr in phx.other:
            assert isinstance(otr, PX.Other)

    def test_Other(self):
        """Instantiation of Other objects."""
        phx = PhyloXMLIO.read(EX_PHYLO)
        otr = phx.other[0]
        assert isinstance(otr, PX.Other)
        assert otr.tag == "alignment"
        assert otr.namespace == "http://example.org/align"
        assert len(otr.children) == 3
        for child, name, value in zip(
            otr,
            ("A", "B", "C"),
            (
                "acgtcgcggcccgtggaagtcctctcct",
                "aggtcgcggcctgtggaagtcctctcct",
                "taaatcgc--cccgtgg-agtccc-cct",
            ),
        ):
            assert child.tag == "seq"
            assert child.attributes["name"] == name
            assert child.value == value

    def test_Phylogeny(self):
        """Instantiation of Phylogeny objects."""
        trees = list(PhyloXMLIO.parse(EX_PHYLO))
        # Monitor lizards
        assert trees[9].name == "monitor lizards"
        assert trees[9].description == "a pylogeny of some monitor lizards"
        assert trees[9].rooted
        # Network (unrooted)
        assert trees[6].name == "network, node B is connected to TWO nodes: AB and C"
        assert not trees[6].rooted

    def test_Clade(self):
        """Instantiation of Clade objects."""
        tree = list(PhyloXMLIO.parse(EX_PHYLO))[6]
        clade_ab, clade_c = tree.clade.clades
        clade_a, clade_b = clade_ab.clades
        for clade, id_source, name, blen in zip(
            (clade_ab, clade_a, clade_b, clade_c),
            ("ab", "a", "b", "c"),
            ("AB", "A", "B", "C"),
            (0.06, 0.102, 0.23, 0.4),
        ):
            assert isinstance(clade, PX.Clade)
            assert clade.id_source == id_source
            assert clade.name == name
            assert clade.branch_length == pytest.approx(blen, abs=5e-8)

    def test_Annotation(self):
        """Instantiation of Annotation objects."""
        tree = list(PhyloXMLIO.parse(EX_PHYLO))[3]
        ann = tree.clade[1].sequences[0].annotations[0]
        assert isinstance(ann, PX.Annotation)
        assert ann.desc == "alcohol dehydrogenase"
        assert ann.confidence.value == pytest.approx(0.67, abs=5e-8)
        assert ann.confidence.type == "probability"

    def test_BinaryCharacters(self):
        """Instantiation of BinaryCharacters objects."""
        with open(EX_DOLLO) as handle:
            tree = next(PhyloXMLIO.parse(handle))
        bchars = tree.clade[0, 0].binary_characters
        assert isinstance(bchars, PX.BinaryCharacters)
        assert bchars.type == "parsimony inferred"
        for name, count, value in (
            ("gained", 2, ["Cofilin_ADF", "Gelsolin"]),
            ("lost", 0, []),
            ("present", 2, ["Cofilin_ADF", "Gelsolin"]),
            ("absent", None, []),
        ):
            assert getattr(bchars, name + "_count") == count
            assert getattr(bchars, name) == value

    # TODO: BranchColor -- see made_up.xml

    def test_CladeRelation(self):
        """Instantiation of CladeRelation objects."""
        tree = list(PhyloXMLIO.parse(EX_PHYLO))[6]
        crel = tree.clade_relations[0]
        assert isinstance(crel, PX.CladeRelation)
        assert crel.id_ref_0 == "b"
        assert crel.id_ref_1 == "c"
        assert crel.type == "network_connection"

    def test_Confidence(self):
        """Instantiation of Confidence objects."""
        with open(EX_MADE) as handle:
            tree = next(PhyloXMLIO.parse(handle))
        assert tree.name == "testing confidence"
        for conf, type, val in zip(
            tree.confidences, ("bootstrap", "probability"), (89.0, 0.71)
        ):
            assert isinstance(conf, PX.Confidence)
            assert conf.type == type
            assert conf.value == pytest.approx(val, abs=5e-8)
        assert tree.clade.name == "b"
        assert tree.clade.width == pytest.approx(0.2, abs=5e-8)
        for conf, val in zip(tree.clade[0].confidences, (0.9, 0.71)):
            assert isinstance(conf, PX.Confidence)
            assert conf.type == "probability"
            assert conf.value == pytest.approx(val, abs=5e-8)

    def test_Date(self):
        """Instantiation of Date objects."""
        tree = list(PhyloXMLIO.parse(EX_PHYLO))[11]
        silurian = tree.clade[0, 0].date
        devonian = tree.clade[0, 1].date
        ediacaran = tree.clade[1].date
        for date, desc, val in zip(
            (silurian, devonian, ediacaran),
            # (10, 20, 30), # range is deprecated
            ("Silurian", "Devonian", "Ediacaran"),
            (425, 320, 600),
        ):
            assert isinstance(date, PX.Date)
            assert date.unit == "mya"
            # self.assertAlmostEqual(date.range, rang)
            assert date.desc == desc
            assert date.value == pytest.approx(val, abs=5e-8)

    def test_Distribution(self):
        """Instantiation of Distribution objects.

        Also checks Point type and safe Unicode handling (?).
        """
        tree = list(PhyloXMLIO.parse(EX_PHYLO))[10]
        hirschweg = tree.clade[0, 0].distributions[0]
        nagoya = tree.clade[0, 1].distributions[0]
        eth_zurich = tree.clade[0, 2].distributions[0]
        san_diego = tree.clade[1].distributions[0]
        for dist, desc, lati, longi, alti in zip(
            (hirschweg, nagoya, eth_zurich, san_diego),
            (
                "Hirschweg, Winterthur, Switzerland",
                "Nagoya, Aichi, Japan",
                "ETH Zürich",
                "San Diego",
            ),
            (47.481277, 35.155904, 47.376334, 32.880933),
            (8.769303, 136.915863, 8.548108, -117.217543),
            (472, 10, 452, 104),
        ):
            assert isinstance(dist, PX.Distribution)
            assert dist.desc == desc
            point = dist.points[0]
            assert isinstance(point, PX.Point)
            assert point.geodetic_datum == "WGS84"
            assert point.lat == lati
            assert point.long == longi
            assert point.alt == alti

    def test_DomainArchitecture(self):
        """Instantiation of DomainArchitecture objects.

        Also checks ProteinDomain type.
        """
        with open(EX_APAF) as handle:
            tree = next(PhyloXMLIO.parse(handle))
        clade = tree.clade[0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        darch = clade.sequences[0].domain_architecture
        assert isinstance(darch, PX.DomainArchitecture)
        assert darch.length == 1249
        for domain, start, end, conf, value in zip(
            darch.domains,
            (6, 109, 605, 647, 689, 733, 872, 993, 1075, 1117, 1168),
            (90, 414, 643, 685, 729, 771, 910, 1031, 1113, 1155, 1204),
            (
                7.0e-26,
                7.2e-117,
                2.4e-6,
                1.1e-12,
                2.4e-7,
                4.7e-14,
                2.5e-8,
                4.6e-6,
                6.3e-7,
                1.4e-7,
                0.3,
            ),
            (
                "CARD",
                "NB-ARC",
                "WD40",
                "WD40",
                "WD40",
                "WD40",
                "WD40",
                "WD40",
                "WD40",
                "WD40",
                "WD40",
            ),
        ):
            assert isinstance(domain, PX.ProteinDomain)
            assert domain.start + 1 == start
            assert domain.end == end
            assert domain.confidence == pytest.approx(conf, abs=5e-8)
            assert domain.value == value

    def test_Events(self):
        """Instantiation of Events objects."""
        tree = list(PhyloXMLIO.parse(EX_PHYLO))[4]
        event_s = tree.clade.events
        assert isinstance(event_s, PX.Events)
        assert event_s.speciations == 1
        event_d = tree.clade[0].events
        assert isinstance(event_d, PX.Events)
        assert event_d.duplications == 1

    def test_Polygon(self):
        """Instantiation of Polygon objects."""
        tree = PhyloXMLIO.read(EX_MADE).phylogenies[1]
        assert tree.name == "testing polygon"
        dist = tree.clade[0].distributions[0]
        for poly in dist.polygons:
            assert isinstance(poly, PX.Polygon)
            assert len(poly.points) == 3
        assert dist.polygons[0].points[0].alt_unit == "m"
        for point, lati, longi, alti in zip(
            chain(dist.polygons[0].points, dist.polygons[1].points),
            (47.481277, 35.155904, 47.376334, 40.481277, 25.155904, 47.376334),
            (8.769303, 136.915863, 8.548108, 8.769303, 136.915863, 7.548108),
            (472, 10, 452, 42, 10, 452),
        ):
            assert isinstance(point, PX.Point)
            assert point.geodetic_datum == "WGS84"
            assert point.lat == lati
            assert point.long == longi
            assert point.alt == alti

    def test_Property(self):
        """Instantiation of Property objects."""
        tree = list(PhyloXMLIO.parse(EX_PHYLO))[8]
        for prop, id_ref, value in zip(
            tree.properties, ("id_a", "id_b", "id_c"), ("1200", "2300", "200")
        ):
            assert isinstance(prop, PX.Property)
            assert prop.id_ref == id_ref
            assert prop.datatype == "xsd:integer"
            assert prop.ref == "NOAA:depth"
            assert prop.applies_to == "node"
            assert prop.unit == "METRIC:m"
            assert prop.value == value

    def test_Reference(self):
        """Instantiation of Reference objects."""
        with open(EX_DOLLO) as handle:
            tree = next(PhyloXMLIO.parse(handle))
        reference = tree.clade[0, 0, 0, 0, 0, 0].references[0]
        assert isinstance(reference, PX.Reference)
        assert reference.doi == "10.1038/nature06614"
        assert reference.desc is None

    def test_Sequence(self):
        """Instantiation of Sequence objects.

        Also checks Accession and Annotation types.
        """
        trees = list(PhyloXMLIO.parse(EX_PHYLO))
        # Simple element with id_source
        seq0 = trees[4].clade[1].sequences[0]
        assert isinstance(seq0, PX.Sequence)
        assert seq0.id_source == "z"
        assert seq0.symbol == "ADHX"
        assert seq0.accession.source == "ncbi"
        assert seq0.accession.value == "Q17335"
        assert seq0.name == "alcohol dehydrogenase"
        assert seq0.annotations[0].ref == "InterPro:IPR002085"
        # More complete elements
        seq1 = trees[5].clade[0, 0].sequences[0]
        seq2 = trees[5].clade[0, 1].sequences[0]
        seq3 = trees[5].clade[1].sequences[0]
        for seq, sym, acc, name, mol_seq, ann_refs in zip(
            (seq1, seq2, seq3),
            ("ADHX", "RT4I1", "ADHB"),
            ("P81431", "Q54II4", "Q04945"),
            (
                "Alcohol dehydrogenase class-3",
                "Reticulon-4-interacting protein 1 homolog, mitochondrial precursor",
                "NADH-dependent butanol dehydrogenase B",
            ),
            (
                "TDATGKPIKCMAAIAWEAKKPLSIEEVEVAPPKSGEVRIKILHSGVCHTD",
                "MKGILLNGYGESLDLLEYKTDLPVPKPIKSQVLIKIHSTSINPLDNVMRK",
                "MVDFEYSIPTRIFFGKDKINVLGRELKKYGSKVLIVYGGGSIKRNGIYDK",
            ),
            (
                ("EC:1.1.1.1", "GO:0004022"),
                ("GO:0008270", "GO:0016491"),
                ("GO:0046872", "KEGG:Tetrachloroethene degradation"),
            ),
        ):
            assert isinstance(seq, PX.Sequence)
            assert seq.symbol == sym
            assert seq.accession.source == "UniProtKB"
            assert seq.accession.value == acc
            assert seq.name == name
            assert seq.mol_seq.value == mol_seq
            assert seq.annotations[0].ref == ann_refs[0]
            assert seq.annotations[1].ref == ann_refs[1]

    def test_SequenceRelation(self):
        """Instantiation of SequenceRelation objects."""
        tree = list(PhyloXMLIO.parse(EX_PHYLO))[4]
        for seqrel, id_ref_0, id_ref_1, type in zip(
            tree.sequence_relations,
            ("x", "x", "y"),
            ("y", "z", "z"),
            ("paralogy", "orthology", "orthology"),
        ):
            assert isinstance(seqrel, PX.SequenceRelation)
            assert seqrel.id_ref_0 == id_ref_0
            assert seqrel.id_ref_1 == id_ref_1
            assert seqrel.type == type

    def test_Taxonomy(self):
        """Instantiation of Taxonomy objects.

        Also checks Id type.
        """
        trees = list(PhyloXMLIO.parse(EX_PHYLO))
        # Octopus
        tax5 = trees[5].clade[0, 0].taxonomies[0]
        assert isinstance(tax5, PX.Taxonomy)
        assert tax5.id.value == "6645"
        assert tax5.id.provider == "NCBI"
        assert tax5.code == "OCTVU"
        assert tax5.scientific_name == "Octopus vulgaris"
        # Nile monitor
        tax9 = trees[9].clade[0].taxonomies[0]
        assert isinstance(tax9, PX.Taxonomy)
        assert tax9.id.value == "62046"
        assert tax9.id.provider == "NCBI"
        assert tax9.scientific_name == "Varanus niloticus"
        assert tax9.common_names[0] == "Nile monitor"
        assert tax9.rank == "species"

    def test_Uri(self):
        """Instantiation of Uri objects."""
        tree = list(PhyloXMLIO.parse(EX_PHYLO))[9]
        uri = tree.clade.taxonomies[0].uri
        assert isinstance(uri, PX.Uri)
        assert uri.desc == "EMBL REPTILE DATABASE"
        assert uri.value == "http://www.embl-heidelberg.de/~uetz/families/Varanidae.html"


# ---------------------------------------------------------
# Serialization tests


class WriterTests(unittest.TestCase):
    """Tests for serialization of objects to phyloXML format.

    Modifies the globally defined filenames in order to run the other parser
    tests on files (re)generated by PhyloXMLIO's own writer.
    """

    def _rewrite_and_call(self, orig_fname, test_cases):
        """Parse, rewrite and retest a phyloXML example file."""
        with open(orig_fname) as infile:
            phx = PhyloXMLIO.read(infile)
        with open(DUMMY, "w") as outfile:
            PhyloXMLIO.write(phx, outfile)
        for cls, tests in test_cases:
            inst = cls("setUp")
            for test in tests:
                if test == "test_Distribution" and platform.system() == "Windows":
                    continue  # Skip, see issue #3944
                getattr(inst, test)()

    def test_apaf(self):
        """Round-trip parsing and serialization of apaf.xml."""
        global EX_APAF
        orig_fname = EX_APAF
        try:
            EX_APAF = DUMMY
            self._rewrite_and_call(
                orig_fname,
                (
                    (
                        ParseTests,
                        ["test_read_apaf", "test_parse_apaf", "test_shape_apaf"],
                    ),
                    (TreeTests, ["test_DomainArchitecture"]),
                ),
            )
        finally:
            EX_APAF = orig_fname

    def test_bcl2(self):
        """Round-trip parsing and serialization of bcl_2.xml."""
        global EX_BCL2
        orig_fname = EX_BCL2
        try:
            EX_BCL2 = DUMMY
            self._rewrite_and_call(
                orig_fname,
                (
                    (
                        ParseTests,
                        ["test_read_bcl2", "test_parse_bcl2", "test_shape_bcl2"],
                    ),
                    (TreeTests, ["test_Confidence"]),
                ),
            )
        finally:
            EX_BCL2 = orig_fname

    def test_made(self):
        """Round-trip parsing and serialization of made_up.xml."""
        global EX_MADE
        orig_fname = EX_MADE
        try:
            EX_MADE = DUMMY
            self._rewrite_and_call(
                orig_fname,
                (
                    (ParseTests, ["test_read_made", "test_parse_made"]),
                    (TreeTests, ["test_Confidence", "test_Polygon"]),
                ),
            )
        finally:
            EX_MADE = orig_fname

    def test_phylo(self):
        """Round-trip parsing and serialization of phyloxml_examples.xml."""
        global EX_PHYLO
        orig_fname = EX_PHYLO
        try:
            EX_PHYLO = DUMMY
            self._rewrite_and_call(
                orig_fname,
                (
                    (
                        ParseTests,
                        ["test_read_phylo", "test_parse_phylo", "test_shape_phylo"],
                    ),
                    (
                        TreeTests,
                        [
                            "test_Phyloxml",
                            "test_Other",
                            "test_Phylogeny",
                            "test_Clade",
                            "test_Annotation",
                            "test_CladeRelation",
                            "test_Date",
                            "test_Distribution",
                            "test_Events",
                            "test_Property",
                            "test_Sequence",
                            "test_SequenceRelation",
                            "test_Taxonomy",
                            "test_Uri",
                        ],
                    ),
                ),
            )
        finally:
            EX_PHYLO = orig_fname

    def test_dollo(self):
        """Round-trip parsing and serialization of o_tol_332_d_dollo.xml."""
        global EX_DOLLO
        orig_fname = EX_DOLLO
        try:
            EX_DOLLO = DUMMY
            self._rewrite_and_call(
                orig_fname,
                (
                    (ParseTests, ["test_read_dollo", "test_parse_dollo"]),
                    (TreeTests, ["test_BinaryCharacters"]),
                ),
            )
        finally:
            EX_DOLLO = orig_fname


# ---------------------------------------------------------
# Method tests


class MethodTests(unittest.TestCase):
    """Tests for methods on specific classes/objects."""

    def setUp(self):
        self.phyloxml = PhyloXMLIO.read(EX_PHYLO)

    # Type conversions

    def test_clade_to_phylogeny(self):
        """Convert a Clade object to a new Phylogeny."""
        clade = self.phyloxml.phylogenies[0].clade[0]
        tree = clade.to_phylogeny(rooted=True)
        assert isinstance(tree, PX.Phylogeny)

    def test_phylogeny_to_phyloxml(self):
        """Convert a Phylogeny object to a new Phyloxml."""
        tree = self.phyloxml.phylogenies[0]
        doc = tree.to_phyloxml_container()
        assert isinstance(doc, PX.Phyloxml)

    def test_sequence_conversion(self):
        pseq = PX.Sequence(
            type="protein",
            # id_ref=None,
            # id_source=None,
            symbol="ADHX",
            accession=PX.Accession("P81431", source="UniProtKB"),
            name="Alcohol dehydrogenase class-3",
            # location=None,
            mol_seq=PX.MolSeq("TDATGKPIKCMAAIAWEAKKPLSIEEVEVAPPKSGEVRIKILHSGVCHTD"),
            uri=None,
            annotations=[
                PX.Annotation(ref="EC:1.1.1.1"),
                PX.Annotation(ref="GO:0004022"),
            ],
            domain_architecture=PX.DomainArchitecture(
                length=50,
                domains=[
                    PX.ProteinDomain(*args)
                    for args in (
                        # value, start, end, confidence
                        ("FOO", 0, 5, 7.0e-26),
                        ("BAR", 8, 13, 7.2e-117),
                        ("A-OK", 21, 34, 2.4e-06),
                        ("WD40", 40, 50, 0.3),
                    )
                ],
            ),
        )
        srec = pseq.to_seqrecord()
        # TODO: check seqrec-specific traits (see args)
        #   Seq(letters), id, name, description, features
        pseq2 = PX.Sequence.from_seqrecord(srec)
        # TODO: check the round-tripped attributes again

    def test_to_alignment(self):
        tree = self.phyloxml.phylogenies[0]
        aln = tree.to_alignment()
        assert isinstance(aln, MultipleSeqAlignment)
        assert len(aln) == 0
        # Add sequences to the terminals
        for tip, seqstr in zip(tree.get_terminals(), ("AA--TTA", "AA--TTG", "AACCTTC")):
            tip.sequences.append(
                PX.Sequence.from_seqrecord(
                    SeqRecord(Seq(seqstr), id=str(tip)), is_aligned=True
                )
            )
        # Check the alignment
        aln = tree.to_alignment()
        assert isinstance(aln, MultipleSeqAlignment)
        assert len(aln) == 3
        assert aln.get_alignment_length() == 7

    def test_alignment(self):
        tree = self.phyloxml.phylogenies[0]
        aln = tree.alignment
        assert isinstance(aln, Alignment)
        assert len(aln) == 0
        assert aln.shape == (0, 0)
        # Add sequences to the terminals
        for tip, seqstr in zip(tree.get_terminals(), ("AA--TTA", "AA--TTG", "AACCTTC")):
            tip.sequences.append(
                PX.Sequence.from_seqrecord(
                    SeqRecord(Seq(seqstr), id=str(tip)), is_aligned=True
                )
            )
        # Check the alignment
        aln = tree.alignment
        assert isinstance(aln, Alignment)
        assert aln.shape == (3, 7)
        assert aln.sequences[0].id == ":A"
        assert aln.sequences[1].id == ":B"
        assert aln.sequences[2].id == ":C"
        assert aln.sequences[0].seq == "AATTA"
        assert aln.sequences[1].seq == "AATTG"
        assert aln.sequences[2].seq == "AACCTTC"
        assert aln[0] == "AA--TTA"
        assert aln[1] == "AA--TTG"
        assert aln[2] == "AACCTTC"
        assert str(aln) == """\
:A                0 AA--TTA 5
:B                0 AA--TTG 5
:C                0 AACCTTC 7
"""

    # Syntax sugar

    def test_clade_getitem(self):
        """Clade.__getitem__: get sub-clades by extended indexing."""
        tree = self.phyloxml.phylogenies[3]
        assert tree.clade[0, 0] == tree.clade.clades[0].clades[0]
        assert tree.clade[0, 1] == tree.clade.clades[0].clades[1]
        assert tree.clade[1] == tree.clade.clades[1]
        assert len(tree.clade[:]) == len(tree.clade.clades)
        assert len(tree.clade[0, :]) == len(tree.clade.clades[0].clades)

    def test_phyloxml_getitem(self):
        """Phyloxml.__getitem__: get phylogenies by name or index."""
        assert self.phyloxml.phylogenies[9] is self.phyloxml[9]
        assert self.phyloxml["monitor lizards"] is self.phyloxml[9]
        assert len(self.phyloxml[:]) == len(self.phyloxml)

    def test_events(self):
        """Events: Mapping-type behavior."""
        evts = self.phyloxml.phylogenies[4].clade.events
        # Container behavior: __len__, __contains__
        assert len(evts) == 1
        assert "speciations" in evts
        assert "duplications" not in evts
        # Attribute access: __get/set/delitem__
        assert evts["speciations"] == 1
        with pytest.raises(KeyError):
            lambda k: evts[k]("duplications")  # noqa: E731
        evts["duplications"] = 3
        assert evts.duplications == 3
        assert len(evts) == 2
        del evts["speciations"]
        assert evts.speciations is None
        assert len(evts) == 1
        # Iteration: __iter__, keys, values, items
        assert list(iter(evts)) == ["duplications"]
        assert list(evts.keys()) == ["duplications"]
        assert list(evts.values()) == [3]
        assert list(evts.items()) == [("duplications", 3)]

    def test_singlular(self):
        """Clade, Phylogeny: Singular properties for plural attributes."""
        conf = PX.Confidence(0.9, "bootstrap")
        taxo = PX.Taxonomy(rank="genus")
        # Clade.taxonomy, Clade.confidence
        clade = PX.Clade(confidences=[conf], taxonomies=[taxo])
        assert clade.confidence.type == "bootstrap"
        assert clade.taxonomy.rank == "genus"
        # raise if len > 1
        clade.confidences.append(conf)
        with pytest.raises(AttributeError):
            getattr(clade, "confidence")
        clade.taxonomies.append(taxo)
        with pytest.raises(AttributeError):
            getattr(clade, "taxonomy")
        # None if []
        clade.confidences = []
        assert clade.confidence is None
        clade.taxonomies = []
        assert clade.taxonomy is None
        # Phylogeny.confidence
        tree = PX.Phylogeny(True, confidences=[conf])
        assert tree.confidence.type == "bootstrap"
        tree.confidences.append(conf)
        with pytest.raises(AttributeError):
            getattr(tree, "confidence")
        tree.confidences = []
        assert tree.confidence is None

    # Other methods

    def test_color_hex(self):
        """BranchColor: to_hex() method."""
        black = PX.BranchColor(0, 0, 0)
        assert black.to_hex() == "#000000"
        white = PX.BranchColor(255, 255, 255)
        assert white.to_hex() == "#ffffff"
        green = PX.BranchColor(14, 192, 113)
        assert green.to_hex() == "#0ec071"


# ---------------------------------------------------------


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
