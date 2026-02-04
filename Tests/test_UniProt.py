"""Tests for the UniProt module."""

import unittest
import pytest
from itertools import islice


from Bio import UniProt

pytestmark = pytest.mark.online



class SearchTests(unittest.TestCase):
    def test_search_result_count(self):
        query = "(organism_id:2697049) AND (reviewed:true)"
        results = UniProt.search(query)
        search_result_count = len(results)

        assert isinstance(search_result_count, int)
        assert search_result_count >= 1
        assert search_result_count == len(list(results))

    def test_search_all_fields(self):
        query = "Insulin AND (reviewed:true)"
        results = list(islice(UniProt.search(query, batch_size=50), 50))
        assert len(results) == 50

        for result in results:
            assert isinstance(result, dict)
            assert len(result) > 0

    def test_search_with_fields(self):
        query = "Insulin AND (reviewed:true)"
        fields = ["id", "accession"]
        results = list(islice(UniProt.search(query, fields=fields, batch_size=50), 50))
        assert len(results) == 50
        expected_keys = {"entryType", "primaryAccession", "uniProtkbId"}

        for result in results:
            assert isinstance(result, dict)
            assert len(result) > 0
            assert set(result.keys()) == expected_keys

    def test_search_results_slicing(self):
        query = "Insulin AND (reviewed:true)"
        results = UniProt.search(query, batch_size=25)
        results_sliced = results[:50]
        assert 50 == len(results_sliced)
        results_list = list(islice(results, 50))
        assert 50 == len(results_list)
        assert results_sliced == results_list

        # We need to fetch the results again to reset the iterator
        results = UniProt.search(query, batch_size=25)
        results_sliced = results[49:30:-1]
        assert 19 == len(results_sliced)
        results_list = list(islice(results, 50))[49:30:-1]
        assert 19 == len(results_list)
        assert results_sliced == results_list

        results_sliced = results[0:0]
        assert 0 == len(results_sliced)

        # This query should return less than 500 results.
        query = "(organism_id:2697049) AND (reviewed:true)"
        results = UniProt.search(query)
        results_sliced = results[-10:-5]
        assert 5 == len(results_sliced)
        results_list = list(results)[-10:-5]
        assert 5 == len(results_list)
        assert results_sliced == results_list

        results = UniProt.search(query)
        results_sliced = results[-1_000_001:1_000_000]
        assert len(results) == len(results_sliced)
        results_list = list(results)[-1_000_001:1_000_000]
        assert len(results_list) == len(results)
        assert results_sliced == results_list

    def test_search_results_indexing(self):
        query = "(organism_id:2697049) AND (reviewed:true)"
        results = UniProt.search(query)
        assert isinstance(results[0], dict)
        assert isinstance(results[-1], dict)
        with pytest.raises(IndexError):
            _ = results[1_000_000]
        with pytest.raises(IndexError):
            _ = results[-1_000_000]
