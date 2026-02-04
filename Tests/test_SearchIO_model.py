# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO objects.

Tests the methods and behaviors of QueryResult, Hit, and HSP objects. All tests
are format-independent and are meant to check the fundamental behavior common
to all formats.

"""

import pickle
import unittest
import pytest
from copy import deepcopy
from io import BytesIO

from search_tests_common import SearchTestBaseClass

from Bio.Align import MultipleSeqAlignment
from Bio.SearchIO._model import Hit
from Bio.SearchIO._model import HSP
from Bio.SearchIO._model import HSPFragment
from Bio.SearchIO._model import QueryResult
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# mock HSPFragments
frag111 = HSPFragment("hit1", "query1", hit="ATGCGCAT", query="ATGCGCAT")
frag112 = HSPFragment("hit1", "query1", hit="ATG", query="GAT")
frag113 = HSPFragment("hit1", "query1", hit="ATTCG", query="AT-CG")
frag113b = HSPFragment("hit1", "query1", hit="ATTCG", query="AT-CG")
frag114 = HSPFragment("hit1", "query1", hit="AT", query="AT")
frag114b = HSPFragment("hit1", "query1", hit="ATCG", query="ATGG")
frag211 = HSPFragment("hit2", "query1", hit="GGGCCC", query="GGGCC-")
frag311 = HSPFragment("hit3", "query1", hit="GATG", query="GTTG")
frag312 = HSPFragment("hit3", "query1", hit="ATATAT", query="ATATAT")
frag411 = HSPFragment("hit4", "query1", hit="CC-ATG", query="CCCATG")
frag121 = HSPFragment("hit1", "query2", hit="GCGAG", query="GCGAC")
# mock HSPs
hsp111 = HSP([frag111])
hsp112 = HSP([frag112])
hsp113 = HSP([frag113, frag113b])
hsp114 = HSP([frag114, frag114b])
hsp211 = HSP([frag211])
hsp311 = HSP([frag311])
hsp312 = HSP([frag312])
hsp411 = HSP([frag411])
hsp121 = HSP([frag121])
# mock Hits
hit11 = Hit([hsp111, hsp112, hsp113, hsp114])
hit21 = Hit([hsp211])
hit31 = Hit([hsp311, hsp312])
hit41 = Hit([hsp411])
hit12 = Hit([hsp121])


class QueryResultCases(SearchTestBaseClass):
    def setUp(self):
        self.qresult = QueryResult([hit11, hit21, hit31], "query1")
        # set mock attributes
        self.qresult.seq_len = 1102
        self.qresult.target = "refseq_rna"

    def test_pickle(self):
        """Test pickling and unpickling of QueryResult."""
        buf = BytesIO()
        pickle.dump(self.qresult, buf)
        unp = pickle.loads(buf.getvalue())
        self.compare_search_obj(self.qresult, unp)

    def test_order(self):
        # added hits should be ordered
        assert self.qresult[0] == hit11
        assert self.qresult[2] == hit31
        # removal of second item should bump item #2 to #1
        del self.qresult["hit2"]
        assert self.qresult[0] == hit11
        assert self.qresult[1] == hit31

    def test_init_none(self):
        """Test QueryResult.__init__, no arguments."""
        qresult = QueryResult()
        assert qresult.id is None
        assert qresult.description is None

    def test_init_id_only(self):
        """Test QueryResult.__init__, with ID only."""
        qresult = QueryResult(id="query1")
        assert "query1" == qresult.id
        assert qresult.description is None

    def test_init_hits_only(self):
        """Test QueryResult.__init__, with hits only."""
        qresult = QueryResult([hit11, hit21, hit31])
        assert "query1" == qresult.id
        assert "<unknown description>" == qresult.description

    def test_repr(self):
        """Test QueryResult.__repr__."""
        assert "QueryResult(id='query1', 3 hits)" == repr(self.qresult)

    def test_iter(self):
        """Test QueryResult.__iter__."""
        # iteration should return hits contained
        for counter, hit in enumerate(self.qresult):
            assert hit in (hit11, hit21, hit31)
        assert 2 == counter

    def test_hits(self):
        """Test QueryResult.hits."""
        # hits should return hits contained in qresult
        hits = list(self.qresult.hits)
        assert [hit11, hit21, hit31] == hits

    def test_hit_keys(self):
        """Test QueryResult.hit_keys."""
        # hit_keys should return hit keys (which default to hit ids)
        hit_keys = list(self.qresult.hit_keys)
        assert ["hit1", "hit2", "hit3"] == hit_keys

    def test_items(self):
        """Test QueryResult.items."""
        # items should return tuples of hit key, hit object pair
        items = list(self.qresult.items)
        assert [("hit1", hit11), ("hit2", hit21), ("hit3", hit31)] == items

    def test_hsps(self):
        """Test QueryResult.hsps."""
        # hsps should return all hsps contained in qresult
        hsps = self.qresult.hsps
        assert [hsp111, hsp112, hsp113, hsp114, hsp211, hsp311, hsp312] == hsps

    def test_fragments(self):
        """Test QueryResult.fragments."""
        # fragments should return all fragments contained in qresult
        frags = self.qresult.fragments
        assert [
                frag111,
                frag112,
                frag113,
                frag113b,
                frag114,
                frag114b,
                frag211,
                frag311,
                frag312,
            ] == frags

    def test_contains(self):
        """Test QueryResult.__contains__."""
        # contains should work with hit ids or hit objects
        assert "hit1" in self.qresult
        assert hit21 in self.qresult
        assert "hit5" not in self.qresult
        assert hit41 not in self.qresult

    def test_contains_alt(self):
        """Test QueryResult.__contains__, with alternative IDs."""
        # contains should work with alternative hit IDs
        hit11._id_alt = ["alt1"]
        query = QueryResult([hit11])
        assert "alt1" in query
        hit11._id_alt = []

    def test_len(self):
        """Test QueryResult.__len__."""
        # len() should return the number of hits contained
        assert 3 == len(self.qresult)

    def test_bool(self):
        """Test QueryResult.__bool__."""
        # should return true only if the qresult has hits
        assert self.qresult
        blank_qresult = QueryResult()
        assert not blank_qresult

    def test_setitem_ok(self):
        """Test QueryResult.__setitem__."""
        # hit objects assignment should work with arbitrary string keys
        self.qresult["hit4"] = hit41
        assert [hit11, hit21, hit31, hit41] == list(self.qresult.hits)
        # and if the key already exist, the object should be overwritten
        self.qresult["hit4"] = hit11
        assert [hit11, hit21, hit31, hit11] == list(self.qresult.hits)

    def test_setitem_ok_alt(self):
        """Test QueryResult.__setitem__, checking alt hit IDs."""
        # hit objects assignment should make alt IDs visible
        hit11._id_alt = ["alt1", "alt11"]
        query = QueryResult()
        query["hit1"] = hit11
        assert hit11 == query["hit1"]
        assert hit11 == query["alt1"]
        assert hit11 == query["alt11"]
        assert hit11.id != "alt1"
        assert hit11.id != "alt11"
        hit11._id_alt = []

    def test_setitem_ok_alt_existing(self):
        """Test QueryResult.__setitem__, existing key."""
        # hit objects assignment on existing hits should also update alt IDs
        hit11._id_alt = ["alt1"]
        hit21._id_alt = ["alt2"]
        query = QueryResult()
        query["hit"] = hit11
        assert hit11 == query["hit"]
        assert hit11 == query["alt1"]
        query["hit"] = hit21
        assert hit21 == query["hit"]
        assert hit21 == query["alt2"]
        with pytest.raises(KeyError):
            query.__getitem__("alt1")
        hit11._id_alt = []
        hit21._id_alt = []

    def test_setitem_ok_alt_ok_promote(self):
        """Test QueryResult.__setitem__, previously alt ID."""
        # hit objects assignment with ID previously existing as alternative
        # should make the ID primary
        hit11._id_alt = ["alt1"]
        hit41._id_alt = ["alt4"]
        hit31._id_alt = ["alt3"]
        query = QueryResult([hit11, hit41])
        assert hit11 == query["alt1"]
        assert hit41 == query["alt4"]
        assert "alt1" not in query._items
        assert "alt1" in query._QueryResult__alt_hit_ids
        query["alt1"] = hit31
        assert hit31 == query["alt1"]
        assert hit41 == query["alt4"]
        assert "alt1" in query._items
        assert "alt1" not in query._QueryResult__alt_hit_ids
        hit11._id_alt = []
        hit41._id_alt = []
        hit31._id_alt = []

    def test_setitem_wrong_key_type(self):
        """Test QueryResult.__setitem__, wrong key type."""
        # item assignment should fail if the key is not string
        with pytest.raises(TypeError):
            self.qresult.__setitem__(0, hit41)
        with pytest.raises(TypeError):
            self.qresult.__setitem__(slice(0, 2), [hit41, hit31])

    def test_setitem_wrong_type(self):
        """Test QueryResult.__setitem__, wrong type."""
        # item assignment should fail if the object assigned is not a hit object
        with pytest.raises(TypeError):
            self.qresult.__setitem__("hit4", hsp111)
        with pytest.raises(TypeError):
            self.qresult.__setitem__("hit5", "hit5")

    def test_setitem_wrong_query_id(self):
        """Test QueryResult.__setitem__, wrong query ID."""
        # item assignment should fail if the hit object does not have the same
        # query id
        with pytest.raises(ValueError):
            self.qresult.__setitem__("hit4", hit12)

    def test_setitem_from_empty(self):
        """Test QueryResult.__setitem__, from empty container."""
        qresult = QueryResult()
        # initial desc and id is None
        assert qresult.id is None
        assert qresult.description is None
        # but changes to the first item's after append
        qresult.append(hit11)
        assert "query1" == qresult.id
        assert "<unknown description>" == qresult.description
        # and remains the same after popping the last item
        qresult.pop()
        assert "query1" == qresult.id
        assert "<unknown description>" == qresult.description

    def test_getitem_default_ok(self):
        """Test QueryResult.__getitem__."""
        # hits should be retrievable by their keys (default to id)
        assert hit21 == self.qresult["hit2"]
        assert hit11 == self.qresult["hit1"]

    def test_getitem_int_ok(self):
        """Test QueryResult.__getitem__, with integer."""
        # hits should be retrievable by their index
        assert hit21 == self.qresult[1]
        assert hit31 == self.qresult[-1]

    def test_getitem_slice_ok(self):
        """Test QueryResult.__getitem__, with slice."""
        # if the index is a slice object, a new qresult object with the same
        # instance attributes should be returned
        assert 1102 == self.qresult.seq_len
        assert "refseq_rna" == self.qresult.target
        new_qresult = self.qresult[1:]
        assert [hit21, hit31] == list(new_qresult.hits)
        assert 1102 == new_qresult.seq_len
        assert "refseq_rna" == new_qresult.target

    def test_getitm_slice_alt_ok(self):
        """Test QueryResult.__getitem__, with slice and alt IDs."""
        # slicing should be reflected in the alt IDs as well
        hit31._id_alt = ["alt3"]
        hit11._id_alt = ["alt1"]
        query = QueryResult([hit31, hit11])
        assert hit11 == query["hit1"]
        assert hit11 == query["alt1"]
        assert hit31 == query["hit3"]
        assert hit31 == query["alt3"]
        query = query[:1]
        assert hit31 == query["hit3"]
        assert hit31 == query["alt3"]
        with pytest.raises(KeyError):
            query.__getitem__("hit1")
        with pytest.raises(KeyError):
            query.__getitem__("alt1")
        hit31._id_alt = []
        hit11._id_alt = []

    def test_getitem_alt_ok(self):
        """Test QueryResult.__getitem__, single item with alternative ID."""
        hit11._id_alt = ["alt1"]
        query = QueryResult([hit11])
        assert hit11 == query["hit1"]
        assert hit11 == query["alt1"]
        assert hit11.id != "alt1"
        hit11._id_alt = []

    def test_delitem_string_ok(self):
        """Test QueryResult.__getitem__, with string."""
        # delitem should work with string index
        del self.qresult["hit1"]
        assert 2 == len(self.qresult)
        assert [hit21, hit31], list(self.qresult.hits)

    def test_delitem_int_ok(self):
        """Test QueryResult.__delitem__."""
        # delitem should work with int index
        del self.qresult[-1]
        assert 2 == len(self.qresult)
        assert [hit11, hit21] == list(self.qresult.hits)
        del self.qresult[0]
        assert 1 == len(self.qresult)
        assert [hit21], list(self.qresult.hits)

    def test_delitem_slice_ok(self):
        """Test QueryResult.__delitem__, with slice."""
        # delitem should work with slice objects
        del self.qresult[:-1]
        assert 1 == len(self.qresult)
        assert [hit31], self.qresult.hits

    def test_delitem_alt_ok(self):
        """Test QueryResult.__delitem__, with alt ID."""
        # delitem should work with alt IDs
        hit31._id_alt = ["alt3"]
        qresult = QueryResult([hit31, hit41])
        assert 2 == len(qresult)
        del qresult["alt3"]
        assert 1 == len(qresult)
        assert hit41 == qresult["hit4"]
        with pytest.raises(KeyError):
            qresult.__getitem__("alt3")
        hit31._id_alt = []

    def test_description_set(self):
        """Test QueryResult.description setter."""
        # setting the description should change the query seqrecord description
        # of the contained hsps, if they have an alignment
        # test for default value
        qresult = deepcopy(self.qresult)
        new_desc = "unicorn hox homolog"
        # test initial condition
        for hit in qresult:
            assert new_desc != hit.query_description
            for hsp in hit:
                assert new_desc != hsp.query_description
                for fragment in hsp:
                    assert new_desc != fragment.query_description
                    assert new_desc != fragment.query.description
        qresult.description = new_desc
        # test after setting
        for hit in qresult:
            assert new_desc == hit.query_description
            for hsp in hit:
                assert new_desc == hsp.query_description
                for fragment in hsp:
                    assert new_desc == fragment.query_description
                    assert new_desc == fragment.query.description

    def test_description_set_no_seqrecord(self):
        """Test QueryResult.description setter, without HSP SeqRecords."""
        frag1 = HSPFragment("hit1", "query")
        frag2 = HSPFragment("hit1", "query")
        frag3 = HSPFragment("hit2", "query")
        hit1 = Hit([HSP([x]) for x in [frag1, frag2]])
        hit2 = Hit([HSP([frag3])])
        qresult = QueryResult([hit1, hit2])
        # test initial condition
        for hit in qresult:
            for hsp in hit.hsps:
                assert getattr(hsp, "query") is None
        qresult.description = "unicorn hox homolog"
        # test after setting
        for hit in qresult:
            for hsp in hit.hsps:
                assert getattr(hsp, "query") is None

    def test_id_set(self):
        """Test QueryResult.id setter."""
        # setting an ID should change the query IDs of all contained Hit and HSPs
        qresult = deepcopy(self.qresult)
        assert "query1" == qresult.id
        for hit in qresult:
            assert "query1" == hit.query_id
            for hsp in hit:
                assert "query1" == hsp.query_id
                for fragment in hsp:
                    assert "query1" == fragment.query_id
                    assert "query1" == fragment.query.id
        qresult.id = "new_id"
        assert "new_id" == qresult.id
        for hit in qresult:
            assert "new_id" == hit.query_id
            for hsp in hit:
                assert "new_id" == hsp.query_id
                for fragment in hsp:
                    assert "new_id" == fragment.query_id
                    assert "new_id" == fragment.query.id

    def test_absorb_hit_does_not_exist(self):
        """Test QueryResult.absorb, hit does not exist."""
        # absorb should work like append when the hit does not exist
        assert [hit11, hit21, hit31] == list(self.qresult.hits)
        self.qresult.absorb(hit41)
        assert [hit11, hit21, hit31, hit41] == list(self.qresult.hits)
        assert ["hit1", "hit2", "hit3", "hit4"] == list(self.qresult.hit_keys)

    def test_absorb_hit_exists(self):
        """Test QueryResult.absorb, hit with the same ID exists."""
        # absorb should combine the hit's hsps if an existing one is present
        assert [hit11, hit21, hit31] == list(self.qresult.hits)
        assert 2 == len(self.qresult["hit3"])
        hit = Hit([HSP([HSPFragment("hit3", "query1")])])
        self.qresult.absorb(hit)
        assert [hit11, hit21, hit31] == list(self.qresult.hits)
        assert ["hit1", "hit2", "hit3"] == list(self.qresult.hit_keys)
        assert 3 == len(self.qresult["hit3"])
        # remove the mock hsp
        del self.qresult["hit3"][-1]

    def test_append_ok(self):
        """Test QueryResult.append."""
        # append should work with Hit objects
        assert [hit11, hit21, hit31] == list(self.qresult.hits)
        self.qresult.append(hit41)
        assert [hit11, hit21, hit31, hit41] == list(self.qresult.hits)
        assert ["hit1", "hit2", "hit3", "hit4"] == list(self.qresult.hit_keys)

    def test_append_custom_hit_key_function_ok(self):
        """Test QueryResult.append, with custom hit key function."""
        self.qresult._hit_key_function = lambda hit: hit.id + "_custom"  # noqa: E731
        # append should assign hit keys according to _hit_key_function
        assert ["hit1", "hit2", "hit3"] == list(self.qresult.hit_keys)
        self.qresult.append(hit41)
        assert ["hit1", "hit2", "hit3", "hit4_custom"] == list(self.qresult.hit_keys)

    def test_append_id_exists(self):
        """Test QueryResult.append, when ID exists."""
        # append should raise an error if hit_key already exists
        with pytest.raises(ValueError):
            self.qresult.append(hit11)

    def test_append_alt_id_exists(self):
        """Test QueryResult.append, when alt ID exists."""
        # append should raise an error if hit_key already exists as alt ID
        hit11._id_alt = ["alt"]
        hit21._id_alt = ["alt"]
        qresult = QueryResult([hit11])
        with pytest.raises(ValueError):
            qresult.append(hit21)
        hit11._id_alt = []
        hit21._id_alt = []

    def test_append_alt_id_exists_alt(self):
        """Test QueryResult.append, when alt ID exists as primary."""
        # append should raise an error if alt ID already exists as primary ID
        hit21._id_alt = ["hit1"]
        qresult = QueryResult([hit11])
        with pytest.raises(ValueError):
            qresult.append(hit21)
        hit21._id_alt = []

    def test_hit_filter(self):
        """Test QueryResult.hit_filter."""
        # hit_filter should return a new QueryResult object (shallow copy),
        assert [hit11, hit21, hit31] == list(self.qresult.hits)
        # filter func: min hit length == 2
        # this would filter out hit21, since it only has 1 HSP
        filter_func = lambda hit: len(hit) >= 2  # noqa: E731
        filtered = self.qresult.hit_filter(filter_func)
        assert [hit11, hit31] == list(filtered.hits)
        # make sure all remaining hits return True for the filter function
        assert all(filter_func(hit) for hit in filtered)
        assert 1102 == filtered.seq_len
        assert "refseq_rna" == filtered.target

    def test_hit_filter_no_func(self):
        """Test QueryResult.hit_filter, without arguments."""
        # when given no arguments, hit_filter should create a new object with
        # the same contents
        filtered = self.qresult.hit_filter()
        self.compare_search_obj(filtered, self.qresult)
        assert id(filtered) != id(self.qresult)
        assert 1102 == filtered.seq_len
        assert "refseq_rna" == filtered.target

    def test_hit_filter_no_filtered(self):
        """Test QueryResult.hit_filter, all hits filtered out."""
        # when the filter filters out all hits, hit_filter should return an
        # empty QueryResult object
        filter_func = lambda hit: len(hit) > 50  # noqa: E731
        filtered = self.qresult.hit_filter(filter_func)
        assert 0 == len(filtered)
        assert isinstance(filtered, QueryResult)
        assert 1102 == filtered.seq_len
        assert "refseq_rna" == filtered.target

    def test_hit_map(self):
        """Test QueryResult.hit_map."""
        # hit_map should apply the given function to all contained Hits
        # deepcopy the qresult since we'll change the objects within
        qresult = deepcopy(self.qresult)
        # map func: capitalize hit IDs

        def map_func(hit):
            hit.id = hit.id.upper()
            return hit

        # test before mapping
        assert "hit1" == qresult[0].id
        assert "hit2" == qresult[1].id
        assert "hit3" == qresult[2].id
        mapped = qresult.hit_map(map_func)
        assert "HIT1" == mapped[0].id
        assert "HIT2" == mapped[1].id
        assert "HIT3" == mapped[2].id
        # and make sure the attributes are transferred
        assert 1102 == mapped.seq_len
        assert "refseq_rna" == mapped.target

    def test_hit_map_no_func(self):
        """Test QueryResult.hit_map, without arguments."""
        # when given no arguments, hit_map should create a new object with
        # the same contents
        mapped = self.qresult.hit_map()
        self.compare_search_obj(mapped, self.qresult)
        assert id(mapped) != id(self.qresult)
        assert 1102 == mapped.seq_len
        assert "refseq_rna" == mapped.target

    def test_hsp_filter(self):
        """Test QueryResult.hsp_filter."""
        # hsp_filter should return a new QueryResult object (shallow copy)
        # and any empty hits should be discarded
        assert [hit11, hit21, hit31] == list(self.qresult.hits)
        # filter func: no '-' in hsp query sequence
        # this would filter out hsp113 and hsp211, effectively removing hit21
        filter_func = lambda hsp: "-" not in hsp.fragments[0].query  # noqa: E731
        filtered = self.qresult.hsp_filter(filter_func)
        assert "hit1" in filtered
        assert "hit2" not in filtered
        assert "hit3" in filtered
        # test hsps in hit11
        assert all(hsp in filtered["hit1"] for hsp in [hsp111, hsp112, hsp114])
        # test hsps in hit31
        assert all(hsp in filtered["hit3"] for hsp in [hsp311, hsp312])

    def test_hsp_filter_no_func(self):
        """Test QueryResult.hsp_filter, no arguments."""
        # when given no arguments, hsp_filter should create a new object with
        # the same contents
        filtered = self.qresult.hsp_filter()
        self.compare_search_obj(filtered, self.qresult)
        assert id(filtered) != id(self.qresult)
        assert 1102 == filtered.seq_len
        assert "refseq_rna" == filtered.target

    def test_hsp_filter_no_filtered(self):
        """Test QueryResult.hsp_filter, all hits filtered out."""
        # when the filter filters out all hits, hsp_filter should return an
        # empty QueryResult object
        filter_func = lambda hsp: len(hsp) > 50  # noqa: E731
        filtered = self.qresult.hsp_filter(filter_func)
        assert 0 == len(filtered)
        assert isinstance(filtered, QueryResult)
        assert 1102 == filtered.seq_len
        assert "refseq_rna" == filtered.target

    def test_hsp_map(self):
        """Test QueryResult.hsp_map."""
        # hsp_map should apply the given function to all contained HSPs
        # deepcopy the qresult since we'll change the objects within
        qresult = deepcopy(self.qresult)
        # apply mock attributes to hsp, for testing mapped hsp attributes
        for hit in qresult:
            for hsp in hit:
                setattr(hsp, "mock", 13)

        # map func: remove first letter of all HSP.aln
        def map_func(hsp):
            mapped_frags = [x[1:] for x in hsp]
            return HSP(mapped_frags)

        mapped = qresult.hsp_map(map_func)
        # make sure old hsp attributes is not transferred to mapped hsps
        for hit in mapped:
            for hsp in hit.hsps:
                assert not hasattr(hsp, "mock")
        # check hsps in hit1
        assert "TGCGCAT" == mapped["hit1"][0][0].hit.seq
        assert "TGCGCAT" == mapped["hit1"][0][0].query.seq
        assert "TG" == mapped["hit1"][1][0].hit.seq
        assert "AT" == mapped["hit1"][1][0].query.seq
        assert "TTCG" == mapped["hit1"][2][0].hit.seq
        assert "T-CG" == mapped["hit1"][2][0].query.seq
        assert "TTCG" == mapped["hit1"][2][1].hit.seq
        assert "T-CG" == mapped["hit1"][2][1].query.seq
        assert "T" == mapped["hit1"][3][0].hit.seq
        assert "T" == mapped["hit1"][3][0].query.seq
        assert "TCG" == mapped["hit1"][3][1].hit.seq
        assert "TGG" == mapped["hit1"][3][1].query.seq
        # check hsps in hit2
        assert "GGCCC" == mapped["hit2"][0][0].hit.seq
        assert "GGCC-" == mapped["hit2"][0][0].query.seq
        # check hsps in hit3
        assert "ATG" == mapped["hit3"][0][0].hit.seq
        assert "TTG" == mapped["hit3"][0][0].query.seq
        assert "TATAT" == mapped["hit3"][1][0].hit.seq
        assert "TATAT" == mapped["hit3"][1][0].query.seq
        # and make sure the attributes are transferred
        assert 1102 == mapped.seq_len
        assert "refseq_rna" == mapped.target

    def test_hsp_map_no_func(self):
        """Test QueryResult.hsp_map, without arguments."""
        # when given no arguments, hit_map should create a new object with
        # the same contents
        mapped = self.qresult.hsp_map()
        self.compare_search_obj(mapped, self.qresult)
        assert id(mapped) != id(self.qresult)
        assert 1102 == mapped.seq_len
        assert "refseq_rna" == mapped.target

    def test_pop_nonexistent_with_default(self):
        """Test QueryResult.pop with default for nonexistent key."""
        default = "An arbitrary default return value for this test only."
        nonexistent_key = "neither a standard nor alternative key"
        hit = self.qresult.pop(nonexistent_key, default)
        assert hit == default

    def test_pop_nonexistent_key(self):
        """Test QueryResult.pop with default for nonexistent key."""
        nonexistent_key = "neither a standard nor alternative key"
        with pytest.raises(KeyError):
            self.qresult.pop(nonexistent_key)

    def test_pop_ok(self):
        """Test QueryResult.pop."""
        assert 3 == len(self.qresult)
        hit = self.qresult.pop()
        assert hit == hit31
        assert [hit11, hit21] == list(self.qresult.hits)

    def test_pop_int_index_ok(self):
        """Test QueryResult.pop, with integer index."""
        # pop should work if given an int index
        assert 3 == len(self.qresult)
        hit = self.qresult.pop(1)
        assert hit == hit21
        assert [hit11, hit31] == list(self.qresult.hits)

    def test_pop_string_index_ok(self):
        """Test QueryResult.pop, with string index."""
        # pop should work if given a string index
        assert 3 == len(self.qresult)
        hit = self.qresult.pop("hit2")
        assert hit == hit21
        assert [hit11, hit31] == list(self.qresult.hits)

    def test_pop_string_alt_ok(self):
        """Test QueryResult.pop, with alternative ID."""
        # pop should work with alternative index
        hit11._id_alt = ["alt1"]
        hit21._id_alt = ["alt2"]
        qresult = QueryResult([hit11, hit21])
        hit = qresult.pop("alt1")
        assert hit == hit11
        assert [hit21] == list(qresult)
        assert "hit1" not in qresult
        hit11._id_alt = []
        hit21._id_alt = []

    def test_index(self):
        """Test QueryResult.index."""
        # index should accept hit objects or hit key strings
        assert 2 == self.qresult.index("hit3")
        assert 2 == self.qresult.index(hit31)

    def test_index_alt(self):
        """Test QueryResult.index, with alt ID."""
        # index should work with alt IDs
        hit11._id_alt = ["alt1"]
        qresult = QueryResult([hit21, hit11])
        assert 1 == qresult.index("alt1")
        hit11._id_alt = []

    def test_index_not_present(self):
        """Test QueryResult.index, when index is not present."""
        with pytest.raises(ValueError):
            self.qresult.index("hit4")
        with pytest.raises(ValueError):
            self.qresult.index(hit41)

    def test_sort_ok(self):
        """Test QueryResult.sort."""
        # sort without any arguments should keep the Hits in the same order
        assert [hit11, hit21, hit31] == list(self.qresult.hits)
        self.qresult.sort()
        assert [hit11, hit21, hit31] == list(self.qresult.hits)

    def test_sort_not_in_place_ok(self):
        """Test QueryResult.sort, not in place."""
        # sort without any arguments should keep the Hits in the same order
        assert [hit11, hit21, hit31] == list(self.qresult.hits)
        sorted_qresult = self.qresult.sort(in_place=False)
        assert [hit11, hit21, hit31] == list(sorted_qresult.hits)
        assert [hit11, hit21, hit31] == list(self.qresult.hits)

    def test_sort_reverse_ok(self):
        """Test QueryResult.sort, reverse."""
        # sorting with reverse=True should return a QueryResult with Hits reversed
        assert [hit11, hit21, hit31] == list(self.qresult.hits)
        self.qresult.sort(reverse=True)
        assert [hit31, hit21, hit11] == list(self.qresult.hits)

    def test_sort_reverse_not_in_place_ok(self):
        """Test QueryResult.sort, reverse, not in place."""
        # sorting with reverse=True should return a QueryResult with Hits reversed
        assert [hit11, hit21, hit31] == list(self.qresult.hits)
        sorted_qresult = self.qresult.sort(reverse=True, in_place=False)
        assert [hit31, hit21, hit11] == list(sorted_qresult.hits)
        assert [hit11, hit21, hit31] == list(self.qresult.hits)

    def test_sort_key_ok(self):
        """Test QueryResult.sort, with custom key."""
        # if custom key is given, sort using it
        key = lambda hit: len(hit)  # noqa: E731
        assert [hit11, hit21, hit31] == list(self.qresult.hits)
        self.qresult.sort(key=key)
        assert [hit21, hit31, hit11] == list(self.qresult.hits)

    def test_sort_key_not_in_place_ok(self):
        """Test QueryResult.sort, with custom key, not in place."""
        # if custom key is given, sort using it
        key = lambda hit: len(hit)  # noqa: E731
        assert [hit11, hit21, hit31] == list(self.qresult.hits)
        sorted_qresult = self.qresult.sort(key=key, in_place=False)
        assert [hit21, hit31, hit11] == list(sorted_qresult.hits)
        assert [hit11, hit21, hit31] == list(self.qresult.hits)


class HitCases(SearchTestBaseClass):
    def setUp(self):
        self.hit = Hit([hsp111, hsp112, hsp113])
        self.hit.evalue = 5e-10
        self.hit.name = "test"

    def test_pickle(self):
        """Test pickling and unpickling of Hit."""
        buf = BytesIO()
        pickle.dump(self.hit, buf)
        unp = pickle.loads(buf.getvalue())
        self.compare_search_obj(self.hit, unp)

    def test_init_none(self):
        """Test Hit.__init__, no arguments."""
        hit = Hit()
        assert hit.id is None
        assert hit.description is None
        assert hit.query_id is None
        assert hit.query_description is None

    def test_init_id_only(self):
        """Test Hit.__init__, with ID only."""
        hit = Hit(id="hit1")
        assert "hit1" == hit.id
        assert hit.description is None
        assert hit.query_id is None
        assert hit.query_description is None

    def test_init_hsps_only(self):
        """Test Hit.__init__, with hsps only."""
        hit = Hit([hsp111, hsp112, hsp113])
        assert "hit1" == hit.id
        assert "<unknown description>" == hit.description
        assert "query1" == hit.query_id  # set from the HSPs
        assert "<unknown description>" == hit.query_description

    def test_repr(self):
        """Test Hit.__repr__."""
        # test for cases with 1 or other alignment numbers
        assert "Hit(id='hit1', query_id='query1', 3 hsps)" == repr(self.hit)

    def test_hsps(self):
        """Test Hit.hsps."""
        # hsps should return the list of hsps contained
        assert [hsp111, hsp112, hsp113] == self.hit.hsps

    def test_fragments(self):
        """Test Hit.fragments."""
        # fragments should return the list of fragments in each hsps
        # as a flat list
        assert [frag111, frag112, frag113, frag113b] == self.hit.fragments

    def test_iter(self):
        """Test Hit.__iter__."""
        # iteration should return hsps contained
        for counter, hsp in enumerate(self.hit):
            assert hsp in [hsp111, hsp112, hsp113]
        assert 2 == counter

    def test_len(self):
        """Test Hit.__len__."""
        # len() on Hit objects should return how many hsps it has
        assert 3 == len(self.hit)

    def test_bool(self):
        """Test Hit.__bool__."""
        # bool() on Hit objects should return True only if hsps is filled
        # which is always true
        assert self.hit

    def test_setitem_single(self):
        """Test Hit.__setitem__, single item."""
        # test regular setitem overwrite
        self.hit[1] = hsp114
        assert self.hit.hsps == [hsp111, hsp114, hsp113]

    def test_item_multiple(self):
        """Test Hit.__setitem__, multiple items."""
        # test iterable setitem
        self.hit[:] = [hsp113, hsp112, hsp111]
        assert self.hit.hsps == [hsp113, hsp112, hsp111]

    def test_getitem_single(self):
        """Test Hit.__getitem__, single item."""
        # getitem using integer index should return a hsp object
        hsp1 = self.hit[0]
        assert hsp111 == hsp1
        hsp3 = self.hit[-1]
        assert hsp113 == hsp3

    def test_getitem_multiple(self):
        """Test Hit.__getitem__, multiple items."""
        # getitem using slices should return another hit object
        # with the hsps sliced accordingly, but other attributes preserved
        new_hit = self.hit[:2]
        assert 2 == len(new_hit)
        assert [hsp111, hsp112] == new_hit.hsps
        assert self.hit.id == new_hit.id
        assert self.hit.query_id == new_hit.query_id
        assert 5e-10 == new_hit.evalue
        assert "test" == new_hit.name

    def test_delitem(self):
        """Test Hit.__delitem__."""
        # test delitem
        del self.hit[0]
        assert 2 == len(self.hit)
        assert [hsp112, hsp113] == self.hit.hsps

    def test_validate_hsp_ok(self):
        """Test Hit._validate_hsp."""
        # validation should pass if item is an hsp object with matching
        # query and hit ids
        # if validation passes, None is returned
        assert self.hit._validate_hsp(hsp114) is None

    def test_validate_hsp_wrong_type(self):
        """Test Hit._validate_hsp, wrong type."""
        # validation should fail if item is not an hsp object
        with pytest.raises(TypeError):
            self.hit._validate_hsp(1)
        with pytest.raises(TypeError):
            self.hit._validate_hsp(Seq(""))

    def test_validate_hsp_wrong_query_id(self):
        """Test Hit._validate_hsp, wrong query ID."""
        # validation should fail if query id does not match
        with pytest.raises(ValueError):
            self.hit._validate_hsp(hsp211)

    def test_validate_hsp_wrong_hit_id(self):
        """Test Hit._validate_hsp, wrong hit ID."""
        # validation should vail if hit id does not match
        with pytest.raises(ValueError):
            self.hit._validate_hsp(hsp121)

    def test_desc_set(self):
        """Test Hit.description setter."""
        # setting the description should change the hit seqrecord description
        # of the contained hsps, if they have an alignment
        # test for default value
        hit = deepcopy(self.hit)
        new_desc = "unicorn hox homolog"
        # test initial condition
        for hsp in hit:
            assert new_desc != hsp.hit_description
            for fragment in hsp:
                assert new_desc != fragment.hit_description
                assert new_desc != fragment.hit.description
        hit.description = new_desc
        # test after setting
        for hsp in hit:
            assert new_desc == hsp.hit_description
            for fragment in hsp:
                assert new_desc == fragment.hit_description
                assert new_desc == fragment.hit.description

    def test_desc_set_no_seqrecord(self):
        """Test Hit.description setter, without HSP SeqRecords."""
        frag1 = HSPFragment("hit1", "query")
        frag2 = HSPFragment("hit1", "query")
        hit = Hit([HSP([x]) for x in [frag1, frag2]])
        new_desc = "unicorn hox homolog"
        # test initial condition
        assert hit.description == "<unknown description>"
        for hsp in hit:
            assert hsp.hit_description == "<unknown description>"
            for fragment in hsp:
                assert hsp.hit_description == "<unknown description>"
        hit.description = new_desc
        # test after setting
        assert hit.description == new_desc
        for hsp in hit:
            assert hsp.hit_description, new_desc
            for fragment in hsp:
                assert hsp.hit_description == new_desc

    def test_id_set(self):
        """Test Hit.id setter."""
        # setting an ID should change the query IDs of all contained HSPs
        hit = deepcopy(self.hit)
        assert "hit1" == hit.id
        for hsp in hit.hsps:
            assert "hit1" == hsp.hit_id
            for fragment in hsp:
                assert fragment.hit_id == "hit1"
                assert fragment.hit.id == "hit1"
        hit.id = "new_id"
        assert "new_id" == hit.id
        for hsp in hit.hsps:
            assert "new_id" == hsp.hit_id
            for fragment in hsp:
                assert fragment.hit_id == "new_id"
                assert fragment.hit.id == "new_id"

    def test_append(self):
        """Test Hit.append."""
        # append should add hits to the last position
        self.hit.append(hsp114)
        assert 4 == len(self.hit)
        assert hsp114 == self.hit[-1]

    def test_filter(self):
        """Test Hit.filter."""
        # filter should return a new QueryResult object (shallow copy),
        assert [hsp111, hsp112, hsp113] == self.hit.hsps
        # filter func: min hsp length == 4
        filter_func = lambda hsp: len(hsp[0]) >= 4  # noqa: E731
        filtered = self.hit.filter(filter_func)
        assert [hsp111, hsp113] == filtered.hsps
        # make sure all remaining hits return True for the filter function
        assert all(filter_func(hit) for hit in filtered)
        assert 5e-10 == filtered.evalue
        assert "test" == filtered.name

    def test_filter_no_func(self):
        """Test Hit.filter, without arguments."""
        # when given no arguments, filter should create a new object with
        # the same contents
        filtered = self.hit.filter()
        self.compare_search_obj(filtered, self.hit)
        assert id(filtered) != id(self.hit)
        assert 5e-10 == filtered.evalue
        assert "test" == filtered.name

    def test_filter_no_filtered(self):
        """Test Hit.hit_filter, all hits filtered out."""
        # when the filter filters out all hits, it should return None
        filter_func = lambda hsp: len(hsp[0]) > 50  # noqa: E731
        filtered = self.hit.filter(filter_func)
        assert filtered is None

    def test_index(self):
        """Test Hit.index."""
        # index should accept hsp objects
        assert 1 == self.hit.index(hsp112)

    def test_index_not_present(self):
        """Test Hit.index, when index is not present."""
        with pytest.raises(ValueError):
            self.hit.index(hsp114)

    def test_map(self):
        """Test Hit.hsp_map."""
        # map should apply the given function to all contained HSPs
        # deepcopy hit since we'll change the objects within
        hit = deepcopy(self.hit)
        # apply mock attributes to hsp, for testing mapped hsp attributes
        for hsp in hit:
            setattr(hsp, "mock", 13)

        # map func: remove first letter of all HSP.alignment
        def map_func(hsp):
            mapped_frags = [x[1:] for x in hsp]
            return HSP(mapped_frags)

        mapped = hit.map(map_func)
        # make sure old hsp attributes is not transferred to mapped hsps
        for hsp in mapped:
            assert not hasattr(hsp, "mock")
        # check hsps in hit1
        assert "TGCGCAT" == mapped[0][0].hit.seq
        assert "TGCGCAT" == mapped[0][0].query.seq
        assert "TG" == mapped[1][0].hit.seq
        assert "AT" == mapped[1][0].query.seq
        assert "TTCG" == mapped[2][0].hit.seq
        assert "T-CG" == mapped[2][0].query.seq
        assert "TTCG" == mapped[2][1].hit.seq
        assert "T-CG" == mapped[2][1].query.seq
        # and make sure the attributes are transferred
        assert 5e-10 == mapped.evalue
        assert "test" == mapped.name

    def test_hsp_map_no_func(self):
        """Test Hit.map, without arguments."""
        # when given no arguments, map should create a new object with
        # the same contents
        mapped = self.hit.map()
        self.compare_search_obj(mapped, self.hit)
        assert id(mapped) != id(self.hit)
        assert 5e-10 == mapped.evalue
        assert "test" == mapped.name

    def test_pop(self):
        """Test Hit.pop."""
        # pop should return the last item by default
        assert hsp113 == self.hit.pop()
        assert hsp111 == self.hit.pop(0)

    def test_sort(self):
        """Test Hit.sort."""
        assert [hsp111, hsp112, hsp113] == self.hit.hsps
        # sort by hsp length
        key = lambda batch_hsp: len(batch_hsp[0])  # noqa: E731
        self.hit.sort(key=key)
        assert [hsp112, hsp113, hsp111] == self.hit.hsps

    def test_sort_not_in_place(self):
        """Test Hit.sort, not in place."""
        assert [hsp111, hsp112, hsp113] == self.hit.hsps
        # sort by hsp length
        key = lambda hsp: len(hsp[0])  # noqa: E731
        sorted_hit = self.hit.sort(key=key, in_place=False)
        assert [hsp112, hsp113, hsp111] == sorted_hit.hsps
        assert [hsp111, hsp112, hsp113] == self.hit.hsps
        assert 5e-10 == sorted_hit.evalue
        assert "test" == sorted_hit.name


class HSPSingleFragmentCases(unittest.TestCase):
    def setUp(self):
        self.frag = HSPFragment("hit_id", "query_id", "ATCAGT", "AT-ACT")
        self.frag.query_start = 0
        self.frag.query_end = 6
        self.frag.hit_start = 15
        self.frag.hit_end = 20
        self.hsp = HSP([self.frag])

    def test_init_no_fragment(self):
        """Test HSP.__init__ without fragments."""
        with pytest.raises(ValueError):
            HSP([])

    def test_len(self):
        """Test HSP.__len__."""
        assert 1 == len(self.hsp)

    def test_fragment(self):
        """Test HSP.fragment property."""
        assert self.frag is self.hsp.fragment

    def test_is_fragmented(self):
        """Test HSP.is_fragmented property."""
        assert not self.hsp.is_fragmented

    def test_seq(self):
        """Test HSP sequence properties."""
        assert "ATCAGT" == self.hsp.hit.seq
        assert "AT-ACT" == self.hsp.query.seq

    def test_alignment(self):
        """Test HSP.alignment property."""
        aln = self.hsp.aln
        assert isinstance(aln, MultipleSeqAlignment)
        assert 2 == len(aln)
        assert "ATCAGT", aln[0].seq
        assert "AT-ACT", aln[1].seq

    def test_aln_span(self):
        """Test HSP.aln_span property."""
        assert 6 == self.hsp.aln_span

    def test_span(self):
        """Test HSP span properties."""
        assert 5 == self.hsp.hit_span
        assert 6 == self.hsp.query_span

    def test_range(self):
        """Test HSP range properties."""
        assert (15, 20) == self.hsp.hit_range
        assert (0, 6) == self.hsp.query_range

    def test_setters_readonly(self):
        """Test HSP read-only properties."""
        read_onlies = ("range", "span", "strand", "frame", "start", "end")
        for seq_type in ("query", "hit"):
            with pytest.raises(AttributeError):
                setattr(self.hsp, seq_type, "A")
            for attr in read_onlies:
                with pytest.raises(AttributeError):
                    setattr(self.hsp, f"{seq_type}_{attr}", 5)
        with pytest.raises(AttributeError):
            setattr(self.hsp, "aln", None)


class HSPMultipleFragmentCases(SearchTestBaseClass):
    def setUp(self):
        self.frag1 = HSPFragment("hit_id", "query_id", "ATCAGT", "AT-ACT")
        self.frag1.query_start = 0
        self.frag1.query_end = 6
        self.frag1.hit_start = 15
        self.frag1.hit_end = 20
        self.frag2 = HSPFragment("hit_id", "query_id", "GGG", "CCC")
        self.frag2.query_start = 10
        self.frag2.query_end = 13
        self.frag2.hit_start = 158
        self.frag2.hit_end = 161
        self.hsp = HSP([self.frag1, self.frag2])

    def test_pickle(self):
        """Test pickling and unpickling of HSP."""
        buf = BytesIO()
        pickle.dump(self.hsp, buf)
        unp = pickle.loads(buf.getvalue())
        self.compare_search_obj(self.hsp, unp)

    def test_len(self):
        """Test HSP.__len__."""
        assert 2 == len(self.hsp)

    def test_getitem(self):
        """Test HSP.__getitem__."""
        assert self.frag1 is self.hsp[0]
        assert self.frag2 is self.hsp[1]

    def test_setitem_single(self):
        """Test HSP.__setitem__, single item."""
        frag3 = HSPFragment("hit_id", "query_id", "AAA", "AAT")
        self.hsp[1] = frag3
        assert 2 == len(self.hsp)
        assert self.frag1 is self.hsp[0]
        assert frag3 is self.hsp[1]

    def test_setitem_multiple(self):
        """Test HSP.__setitem__, multiple items."""
        frag3 = HSPFragment("hit_id", "query_id", "AAA", "AAT")
        frag4 = HSPFragment("hit_id", "query_id", "GGG", "GAG")
        self.hsp[:2] = [frag3, frag4]
        assert 2 == len(self.hsp)
        assert frag3 is self.hsp[0]
        assert frag4 is self.hsp[1]

    def test_delitem(self):
        """Test HSP.__delitem__."""
        del self.hsp[0]
        assert 1 == len(self.hsp)
        assert self.frag2 is self.hsp[0]

    def test_contains(self):
        """Test HSP.__contains__."""
        frag3 = HSPFragment("hit_id", "query_id", "AAA", "AAT")
        assert self.frag1 in self.hsp
        assert frag3 not in self.hsp

    def test_fragments(self):
        """Test HSP.fragments property."""
        assert [self.frag1, self.frag2] == self.hsp.fragments

    def test_is_fragmented(self):
        """Test HSP.is_fragmented property."""
        assert self.hsp.is_fragmented

    def test_seqs(self):
        """Test HSP sequence properties."""
        assert ["ATCAGT", "GGG"] == [x.seq for x in self.hsp.hit_all]
        assert ["AT-ACT", "CCC"] == [x.seq for x in self.hsp.query_all]

    def test_id_desc_set(self):
        """Test HSP query and hit id and description setters."""
        for seq_type in ("query", "hit"):
            for attr in ("id", "description"):
                attr_name = f"{seq_type}_{attr}"
                value = getattr(self.hsp, attr_name)
                if attr == "id":
                    # because we happen to have the same value for
                    # IDs and the actual attribute name
                    assert value == attr_name
                    for fragment in self.hsp:
                        assert getattr(fragment, attr_name) == attr_name
                else:
                    assert value == "<unknown description>"
                    for fragment in self.hsp:
                        assert getattr(fragment, attr_name) == "<unknown description>"
                new_value = "new_" + value
                setattr(self.hsp, attr_name, new_value)
                assert getattr(self.hsp, attr_name) == new_value
                assert getattr(self.hsp, attr_name) != value
                for fragment in self.hsp:
                    assert getattr(fragment, attr_name) == new_value
                    assert getattr(fragment, attr_name) != value

    def test_molecule_type(self):
        """Test HSP.molecule_type getter."""
        assert self.hsp.molecule_type is None

    def test_molecule_type_set(self):
        """Test HSP.molecule_type setter."""
        # test initial values
        assert self.hsp.molecule_type is None
        for frag in self.hsp.fragments:
            assert frag.molecule_type is None
        self.hsp.molecule_type = "DNA"
        # test values after setting
        assert self.hsp.molecule_type == "DNA"
        for frag in self.hsp.fragments:
            assert frag.molecule_type == "DNA"

    def test_range(self):
        """Test HSP range properties."""
        # range on HSP with multiple fragment should give the
        # min start and max end coordinates
        assert (15, 161) == self.hsp.hit_range
        assert (0, 13) == self.hsp.query_range

    def test_ranges(self):
        """Test HSP ranges properties."""
        assert [(15, 20), (158, 161)] == self.hsp.hit_range_all
        assert [(0, 6), (10, 13)] == self.hsp.query_range_all

    def test_span(self):
        """Test HSP span properties."""
        # span is always end - start
        assert 146 == self.hsp.hit_span
        assert 13 == self.hsp.query_span

    def test_setters_readonly(self):
        """Test HSP read-only properties."""
        read_onlies = ("range_all", "strand_all", "frame_all")
        for seq_type in ("query", "hit"):
            for attr in read_onlies:
                with pytest.raises(AttributeError):
                    setattr(self.hsp, f"{seq_type}_{attr}", 5)
        with pytest.raises(AttributeError):
            setattr(self.hsp, "aln_all", None)
        with pytest.raises(AttributeError):
            setattr(self.hsp, "hit_all", None)
        with pytest.raises(AttributeError):
            setattr(self.hsp, "query_all", None)


class HSPFragmentWithoutSeqCases(unittest.TestCase):
    def setUp(self):
        self.fragment = HSPFragment("hit_id", "query_id")

    def test_init(self):
        """Test HSPFragment.__init__ attributes."""
        fragment = HSPFragment("hit_id", "query_id")
        for seq_type in ("query", "hit"):
            assert getattr(fragment, seq_type) is None
            for attr in ("strand", "frame", "start", "end"):
                attr_name = f"{seq_type}_{attr}"
                assert getattr(fragment, attr_name) is None
        assert fragment.aln is None
        assert fragment.molecule_type is None
        assert fragment.aln_annotation == {}

    def test_seqmodel(self):
        """Test HSPFragment sequence attributes, no alignments."""
        # all query, hit, and alignment objects should be None
        assert self.fragment.query is None
        assert self.fragment.hit is None
        assert self.fragment.aln is None

    def test_len(self):
        """Test HSPFragment.__len__, no alignments."""
        with pytest.raises(TypeError):
            len(self)
        # len is a shorthand for .aln_span, and it can be set manually
        self.fragment.aln_span = 5
        assert 5 == len(self.fragment)

    def test_repr(self):
        """Test HSPFragment.__repr__, no alignments."""
        # test for minimum repr
        assert "HSPFragment(hit_id='hit_id', query_id='query_id')" == repr(self.fragment)
        self.fragment.aln_span = 5
        assert "HSPFragment(hit_id='hit_id', query_id='query_id', 5 columns)" == repr(self.fragment)

    def test_getitem(self):
        """Test HSPFragment.__getitem__, no alignments."""
        # getitem not supported without alignment
        with pytest.raises(TypeError):
            self.fragment.__getitem__(0)
        with pytest.raises(TypeError):
            self.fragment.__getitem__(slice(0, 2))

    def test_getitem_only_query(self):
        """Test HSPFragment.__getitem__, only query."""
        # getitem should work if only query is present
        self.fragment.query = "AATCG"
        assert "ATCG" == self.fragment[1:].query.seq

    def test_getitem_only_hit(self):
        """Test HSPFragment.__getitem__, only hit."""
        # getitem should work if only query is present
        self.fragment.hit = "CATGC"
        assert "ATGC" == self.fragment[1:].hit.seq

    def test_iter(self):
        """Test HSP.__iter__, no alignments."""
        # iteration not supported
        with pytest.raises(TypeError):
            iter(self)


class HSPFragmentCases(SearchTestBaseClass):
    def setUp(self):
        self.fragment = HSPFragment(
            "hit_id", "query_id", "ATGCTAGCTACA", "ATG--AGCTAGG"
        )

    def test_pickle(self):
        """Test pickling and unpickling of HSPFragment."""
        buf = BytesIO()
        pickle.dump(self.fragment, buf)
        unp = pickle.loads(buf.getvalue())
        self.compare_search_obj(self.fragment, unp)

    def test_init_with_seqrecord(self):
        """Test HSPFragment.__init__, with SeqRecord."""
        # init should work with seqrecords
        hit_seq = SeqRecord(Seq("ATGCTAGCTACA"))
        query_seq = SeqRecord(Seq("ATG--AGCTAGG"))
        hsp = HSPFragment("hit_id", "query_id", hit_seq, query_seq)
        assert isinstance(hsp.query, SeqRecord)
        assert isinstance(hsp.hit, SeqRecord)
        assert isinstance(hsp.aln, MultipleSeqAlignment)

    def test_init_wrong_seqtypes(self):
        """Test HSPFragment.__init__, wrong sequence argument types."""
        # init should only work with string or seqrecords
        wrong_query = Seq("ATGC")
        wrong_hit = Seq("ATGC")
        with pytest.raises(TypeError):
            HSPFragment("hit_id", "query_id", wrong_hit, wrong_query)

    def test_seqmodel(self):
        """Test HSPFragment sequence attribute types and default values."""
        # check hit
        assert isinstance(self.fragment.hit, SeqRecord)
        assert "<unknown description>" == self.fragment.hit.description
        assert "aligned hit sequence" == self.fragment.hit.name
        assert self.fragment.hit.annotations["molecule_type"] is None
        # check query
        assert isinstance(self.fragment.query, SeqRecord)
        assert "<unknown description>" == self.fragment.query.description
        assert "aligned query sequence" == self.fragment.query.name
        assert self.fragment.query.annotations["molecule_type"] is None
        # check alignment
        assert isinstance(self.fragment.aln, MultipleSeqAlignment)
        with pytest.raises(AttributeError):
            self.fragment.aln.molecule_type

    def test_molecule_type_no_seq(self):
        """Test HSPFragment molecule_type property, query and hit sequences not present."""
        assert self.fragment.molecule_type is None
        self.fragment.molecule_type = "DNA"
        assert self.fragment.molecule_type == "DNA"

    def test_molecule_type_with_seq(self):
        """Test HSPFragment molecule_type property, query or hit sequences present."""
        assert self.fragment.molecule_type is None
        self.fragment._hit = SeqRecord(Seq("AAA"))
        self.fragment._query = SeqRecord(Seq("AAA"))
        self.fragment.molecule_type = "DNA"
        assert self.fragment.molecule_type == "DNA"
        assert self.fragment.hit.annotations["molecule_type"] == "DNA"
        assert self.fragment.query.annotations["molecule_type"] == "DNA"

    def test_seq_unequal_hit_query_len(self):
        """Test HSPFragment sequence setter with unequal hit and query lengths."""
        for seq_type in ("hit", "query"):
            opp_type = "query" if seq_type == "hit" else "hit"
            # reset values first
            fragment = HSPFragment("hit_id", "query_id")
            # and test it against the opposite
            setattr(fragment, seq_type, "ATGCACAACAGGA")
            with pytest.raises(ValueError):
                setattr(fragment, opp_type, "ATGCGA")

    def test_len(self):
        """Test HSPFragment.__len__."""
        # len should equal alignment column length
        assert 12 == len(self.fragment)

    def test_repr(self):
        """Test HSPFragment.__repr__."""
        # test for minimum repr
        assert "HSPFragment(hit_id='hit_id', query_id='query_id', 12 columns)" == repr(self.fragment)

    def test_getitem(self):
        """Test HSPFragment.__getitem__."""
        # getitem is supported when alignment is present
        sliced_fragment = self.fragment[:5]
        assert isinstance(sliced_fragment, HSPFragment)
        assert 5 == len(sliced_fragment)
        assert "ATGCT" == sliced_fragment.hit.seq
        assert "ATG--" == sliced_fragment.query.seq

    def test_getitem_attrs(self):
        """Test HSPFragment.__getitem__, with attributes."""
        # attributes from the original instance should not be present in the new
        # objects, except for query, hit, and alignment - related attributes
        setattr(self.fragment, "attr_original", 1000)
        setattr(self.fragment, "hit_description", "yeah")
        setattr(self.fragment, "hit_strand", 1)
        setattr(self.fragment, "query_frame", None)
        # test values prior to slicing
        assert 1000 == getattr(self.fragment, "attr_original")
        assert "yeah" == getattr(self.fragment, "hit_description")
        assert 1 == getattr(self.fragment, "hit_strand")
        assert getattr(self.fragment, "query_frame") is None
        new_hsp = self.fragment[:5]
        # test values after slicing
        assert not hasattr(new_hsp, "attr_original")
        assert 1000 == getattr(self.fragment, "attr_original")
        assert "yeah" == getattr(self.fragment, "hit_description")
        assert 1 == getattr(self.fragment, "hit_strand")
        assert getattr(self.fragment, "query_frame") is None

    def test_getitem_alignment_annot(self):
        """Test HSPFragment.__getitem__, with alignment annotation."""
        # the alignment is annotated, it should be sliced accordingly
        # and transferred to the new object
        setattr(self.fragment, "aln_annotation", {"test": "182718738172"})
        new_hsp = self.fragment[:5]
        assert "18271" == new_hsp.aln_annotation["test"]

    def test_default_attrs(self):
        """Test HSPFragment attributes' default values."""
        fragment = HSPFragment()
        assert "<unknown id>" == fragment.hit_id
        assert "<unknown id>" == fragment.query_id
        assert "<unknown description>" == fragment.hit_description
        assert "<unknown description>" == fragment.query_description
        assert fragment.hit is None
        assert fragment.query is None
        assert fragment.aln is None
        assert [] == fragment.hit_features
        assert [] == fragment.query_features
        assert fragment.hit_strand is None
        assert fragment.query_strand is None
        assert fragment.hit_frame is None
        assert fragment.query_frame is None

    def test_id_desc_set(self):
        """Test HSPFragment query and hit id and description setters."""
        for seq_type in ("query", "hit"):
            for attr in ("id", "description"):
                attr_name = f"{seq_type}_{attr}"
                value = getattr(self.fragment, attr_name)
                if attr == "id":
                    # because we happen to have the same value for
                    # IDs and the actual attribute name
                    assert value == attr_name
                else:
                    assert value == "<unknown description>"
                new_value = "new_" + value
                setattr(self.fragment, attr_name, new_value)
                assert getattr(self.fragment, attr_name) == new_value
                assert getattr(self.fragment, attr_name) != value

    def test_frame_set_ok(self):
        """Test HSPFragment query and hit frame setters."""
        attr = "frame"
        for seq_type in ("query", "hit"):
            attr_name = f"{seq_type}_{attr}"
            for value in (-3, -2, -1, 0, 1, 2, 3, None):
                setattr(self.fragment, attr_name, value)
                assert value == getattr(self.fragment, attr_name)

    def test_frame_set_error(self):
        """Test HSPFragment query and hit frame setters, invalid values."""
        attr = "frame"
        for seq_type in ("query", "hit"):
            func_name = f"_{seq_type}_{attr}_set"
            func = getattr(self.fragment, func_name)
            for value in ("3", "+3", "-2", "plus"):
                with pytest.raises(ValueError):
                    func(value)

    def test_strand_set_ok(self):
        """Test HSPFragment query and hit strand setters."""
        attr = "strand"
        for seq_type in ("query", "hit"):
            attr_name = f"{seq_type}_{attr}"
            for value in (-1, 0, 1, None):
                setattr(self.fragment, attr_name, value)
                assert value == getattr(self.fragment, attr_name)

    def test_strand_set_error(self):
        """Test HSPFragment query and hit strand setters, invalid values."""
        attr = "strand"
        for seq_type in ("query", "hit"):
            func_name = f"_{seq_type}_{attr}_set"
            func = getattr(self.fragment, func_name)
            for value in (3, "plus", "minus", "-", "+"):
                with pytest.raises(ValueError):
                    func(value)

    def test_strand_set_from_plus_frame(self):
        """Test HSPFragment query and hit strand getters, from plus frame."""
        for seq_type in ("query", "hit"):
            attr_name = f"{seq_type}_strand"
            assert getattr(self.fragment, attr_name) is None
            setattr(self.fragment, f"{seq_type}_frame", 3)
            assert 1 == getattr(self.fragment, attr_name)

    def test_strand_set_from_minus_frame(self):
        """Test HSPFragment query and hit strand getters, from minus frame."""
        for seq_type in ("query", "hit"):
            attr_name = f"{seq_type}_strand"
            assert getattr(self.fragment, attr_name) is None
            setattr(self.fragment, f"{seq_type}_frame", -2)
            assert -1 == getattr(self.fragment, attr_name)

    def test_strand_set_from_zero_frame(self):
        """Test HSPFragment query and hit strand getters, from zero frame."""
        for seq_type in ("query", "hit"):
            attr_name = f"{seq_type}_strand"
            assert getattr(self.fragment, attr_name) is None
            setattr(self.fragment, f"{seq_type}_frame", 0)
            assert 0 == getattr(self.fragment, attr_name)

    def test_coords_setters_getters(self):
        """Test HSPFragment query and hit coordinate-related setters and getters."""
        for seq_type in ("query", "hit"):
            attr_start = f"{seq_type}_{'start'}"
            attr_end = f"{seq_type}_{'end'}"
            setattr(self.fragment, attr_start, 9)
            setattr(self.fragment, attr_end, 99)
            # check for span value
            span = getattr(self.fragment, f"{seq_type}_span")
            assert 90 == span
            # and range as well
            range = getattr(self.fragment, f"{seq_type}_range")
            assert (9, 99) == range

    def test_coords_setters_readonly(self):
        """Test HSPFragment query and hit coordinate-related read-only getters."""
        read_onlies = ("range", "span")
        for seq_type in ("query", "hit"):
            for attr in read_onlies:
                with pytest.raises(AttributeError):
                    setattr(self.fragment, f"{seq_type}_{attr}", 5)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
