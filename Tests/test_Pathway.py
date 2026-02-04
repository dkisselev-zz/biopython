# Copyright 2001 by Tarjei Mikkelsen.  All rights reserved.
# Revisions copyright 2018 by Maximilian Greil. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# python unittest framework
"""Tests for Pathway module."""

import unittest
import pytest

# modules to be tested
from Bio.Pathway import Reaction
from Bio.Pathway.Rep.Graph import Graph
from Bio.Pathway.Rep.MultiGraph import MultiGraph


class GraphTestCase(unittest.TestCase):
    def test_Equals(self):
        a = Graph(["a", "b", "c"])
        a.add_edge("a", "b", "label1")
        a.add_edge("b", "c", "label1")
        a.add_edge("b", "a", "label2")
        b = Graph(["a", "b", "c"])
        assert a != b, "equal to similar nodes, no edges"
        b.add_edge("a", "b", "label1")
        assert a != b, "equal to similar nodes, edge subset"
        b.add_edge("b", "c", "label1")
        b.add_edge("b", "a", "label2")
        assert a == b, "not equal to similar"
        c = Graph(["a", "b", "c"])
        c.add_edge("a", "b", "label2")
        c.add_edge("b", "c", "label2")
        c.add_edge("b", "a", "label1")
        assert a != c, "equal to similar with different labels"
        assert c != Graph(), "equal to empty graph"
        assert Graph() == Graph(), "empty graph not equal to self"

    def test_Nodes(self):
        a = Graph()
        assert a.nodes() == [], "default graph not empty"
        a.add_node("a")
        assert a.nodes() == ["a"], "one node not added"
        a.add_node("a")
        assert a.nodes() == ["a"], "duplicate node added"
        a.add_node("b")
        assert sorted(a.nodes()) == ["a", "b"], "second node not added"

    def test_Edges(self):
        a = Graph(["a", "b", "c", "d"])
        a.add_edge("a", "b", "label1")
        assert a.child_edges("a") == [("b", "label1")]
        a.add_edge("b", "a", "label2")
        assert a.parent_edges("a") == [("b", "label2")]
        a.add_edge("b", "c", "label3")
        assert a.parent_edges("c") == [("b", "label3")]
        assert sorted(a.children("b")) == ["a", "c"], "incorrect children"
        assert a.children("d") == [], "incorrect children for singleton"
        assert a.parents("a") == ["b"], "incorrect parents"

    def test_RemoveNode(self):
        a = Graph(["a", "b", "c", "d", "e"])
        a.add_edge("a", "e", "label1")
        a.add_edge("b", "e", "label1")
        a.add_edge("c", "e", "label2")
        a.add_edge("d", "e", "label3")
        a.add_edge("e", "d", "label4")
        a.add_edge("a", "b", "label5")
        a.remove_node("e")
        b = Graph(["a", "b", "c", "d"])
        b.add_edge("a", "b", "label5")
        assert a == b  # , "incorrect node removal")

    def testAdditionalFunctions(self):
        a = Graph(["a", "b", "c"])
        a.add_edge("a", "b", "label1")
        a.add_edge("b", "c", "label1")
        a.add_edge("b", "a", "label2")
        assert str(a) == "<Graph: 3 node(s), 3 edge(s), 2 unique label(s)>"
        assert repr(a) == "<Graph: ('a': ('b', 'label1'))('b': ('a', 'label2'),('c', 'label1'))('c': )>"
        assert a.edges("label1") == [("a", "b"), ("b", "c")]
        assert a.labels() == ["label1", "label2"]


class MultiGraphTestCase(unittest.TestCase):
    def test_Equals(self):
        a = MultiGraph(["a", "b", "c"])
        a.add_edge("a", "b", "label1")
        a.add_edge("b", "c", "label1")
        a.add_edge("b", "a", "label2")
        b = MultiGraph(["a", "b", "c"])
        assert a != b, "equal to similar nodes, no edges"
        b.add_edge("a", "b", "label1")
        assert a != b, "equal to similar nodes, edge subset"
        b.add_edge("b", "c", "label1")
        b.add_edge("b", "a", "label2")
        assert a == b, "not equal to similar"
        c = MultiGraph(["a", "b", "c"])
        c.add_edge("a", "b", "label2")
        c.add_edge("b", "c", "label2")
        c.add_edge("b", "a", "label1")
        assert a != c, "equal to similar with different labels"
        assert c != MultiGraph(), "equal to empty graph"
        assert MultiGraph() == MultiGraph(), "empty graph not equal to self"

    def test_Nodes(self):
        a = MultiGraph()
        assert a.nodes() == [], "default graph not empty"
        a.add_node("a")
        assert a.nodes() == ["a"], "one node not added"
        a.add_node("a")
        assert a.nodes() == ["a"], "duplicate node added"
        a.add_node("b")
        assert sorted(a.nodes()) == ["a", "b"], "second node not added"

    def test_Edges(self):
        a = MultiGraph(["a", "b", "c", "d"])
        a.add_edge("a", "b", "label1")
        assert a.child_edges("a") == [("b", "label1")]
        a.add_edge("a", "b", "label2")
        assert sorted(a.child_edges("a")) == [("b", "label1"), ("b", "label2")]
        a.add_edge("b", "a", "label2")
        assert a.parent_edges("a") == [("b", "label2")]
        a.add_edge("b", "c", "label3")
        assert a.parent_edges("c") == [("b", "label3")]
        children = a.children("b")
        children.sort()
        assert children == ["a", "c"], "incorrect children"
        assert a.children("d") == [], "incorrect children for singleton"
        assert a.parents("a") == ["b"], "incorrect parents"

    def test_RemoveNode(self):
        a = MultiGraph(["a", "b", "c", "d", "e"])
        a.add_edge("a", "e", "label1")
        assert repr(a) == "<MultiGraph: ('a': ('e', 'label1'))('b': )('c': )('d': )('e': )>"
        a.add_edge("b", "e", "label1")
        a.add_edge("c", "e", "label2")
        a.add_edge("d", "e", "label3")
        a.add_edge("e", "d", "label4")
        a.add_edge("a", "b", "label5")
        assert repr(a) == "<MultiGraph: ('a': ('b', 'label5'),('e', 'label1'))('b': ('e', 'label1'))('c': ('e', 'label2'))('d': ('e', 'label3'))('e': ('d', 'label4'))>"
        a.remove_node("e")
        assert repr(a) == "<MultiGraph: ('a': ('b', 'label5'))('b': )('c': )('d': )>"
        b = MultiGraph(["a", "b", "c", "d"])
        b.add_edge("a", "b", "label5")
        assert repr(b) == "<MultiGraph: ('a': ('b', 'label5'))('b': )('c': )('d': )>"
        assert repr(a) == repr(b)
        assert a == b  # , "incorrect node removal")

    def testAdditionalFunctions(self):
        a = MultiGraph(["a", "b", "c"])
        a.add_edge("a", "b", "label1")
        a.add_edge("b", "c", "label1")
        a.add_edge("b", "a", "label2")
        assert str(a) == "<MultiGraph: 3 node(s), 3 edge(s), 2 unique label(s)>"
        assert a.edges("label1") == [("a", "b"), ("b", "c")]
        assert a.labels() == ["label1", "label2"]


class ReactionTestCase(unittest.TestCase):
    def setUp(self):
        self.r_empty = Reaction()
        self.r_prod = Reaction({"a": 1})
        self.r_dest = Reaction({"a": -1})
        self.r_1 = Reaction({"a": -1, "b": 1})
        self.r_1i = Reaction({"a": -1, "b": 1, "c": 0})
        self.r_2 = Reaction({"b": -1, "c": 1})
        self.r_3 = Reaction({"a": -1, "d": 2})
        self.r_4 = Reaction({"c": -1, "d": -1, "a": 1, "e": 2})

    def test_eq(self):
        assert self.r_1 == self.r_1i  # , "not equal to similar")
        assert self.r_3 != self.r_4  # , "equal to different")

    def test_rev(self):
        assert self.r_empty.reverse() == self.r_empty, "empty reversed not empty"
        assert self.r_prod.reverse() == self.r_dest, "reversed reaction not equal to similar"
        assert self.r_4.reverse() == Reaction({"c": 1, "d": 1, "a": -1, "e": -2}), "reversed reaction not equal to similar"
        assert self.r_3.reverse().reverse() == self.r_3, "double reversal not identity"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
