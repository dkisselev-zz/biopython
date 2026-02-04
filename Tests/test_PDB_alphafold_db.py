"""Tests for the alphafold_db module."""

import unittest
import pytest

from Bio.PDB import alphafold_db

pytestmark = pytest.mark.online



class AlphafoldDBTests(unittest.TestCase):
    def test_get_predictions(self):
        predictions = alphafold_db.get_predictions("P00520")

        for prediction in predictions:
            assert isinstance(prediction, dict)
            assert len(prediction) > 0

    def test_get_mmcif_file_path_for(self):
        prediction = {
            "cifUrl": "https://alphafold.ebi.ac.uk/files/AF-P00520-F1-model_v4.cif",
        }
        file_path = alphafold_db._get_mmcif_file_path_for(prediction, "test")
        assert file_path == "test/AF-P00520-F1-model_v4.cif"
