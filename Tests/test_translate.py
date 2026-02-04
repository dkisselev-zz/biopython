# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""Tests of the transcription and translation methods of Seq objects."""

import unittest
import pytest

from Bio import Seq


class TestTranscriptionTranslation(unittest.TestCase):
    def test_transcription(self):
        s = "ATA"
        dna = Seq.Seq(s)
        rna = dna.transcribe()
        assert rna == "AUA"
        s = "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATATT"
        dna = Seq.Seq(s)
        rna = dna.transcribe()
        assert rna == "GAAAAUUCAUUUUCUUUGGACUUUCUCUGAAAUCCGAGUCCUAGGAAAGAUGCGUGAGAUUCUUCAUAUU"
        s = "GAAAAUUCAUUUUCUUUGGACUUUCUCUGAAAUCCGAGUCCUAGGAAAGAUGCGUGAGAUUCUUCAUAUU"
        rna = Seq.Seq(s)
        dna = rna.back_transcribe()
        assert dna == "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATATT"

    def test_translation(self):
        s = ""
        dna = Seq.Seq(s)
        protein = dna.translate(to_stop=True)
        assert protein == ""
        s = "TAA"
        dna = Seq.Seq(s)
        protein = dna.translate(to_stop=True)
        assert protein == ""
        s = "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCA"
        dna = Seq.Seq(s)
        protein = dna.translate(to_stop=True)
        assert protein == "ENSFSLDFL"
        s = "GAA"
        dna = Seq.Seq(s)
        protein = dna.translate(15, to_stop=True)
        assert protein == "E"
        s = "ATA"
        dna = Seq.Seq(s)
        protein = dna.translate("Vertebrate Mitochondrial", to_stop=True)
        assert protein == "M"
        s = "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATAT"
        dna = Seq.Seq(s)
        protein = dna.translate("SGC8", to_stop=True)
        assert protein == "ENSFSLDFLWNPSPSNDAWDSSY"

    def test_dna_rna_translation(self):
        s = "TCAAAAAGGTGCATCTAGATG"
        dna = Seq.Seq(s)
        protein = dna.translate(to_stop=True)
        assert protein == "SKRCI"
        gapped_protein = dna.translate()
        assert gapped_protein == "SKRCI*M"
        # The table used here has "AGG" as a stop codon:
        p2 = dna.translate(table=2, to_stop=True)
        assert p2 == "SK"
        p2 = dna.translate(table=2)
        assert p2 == "SK*CI*M"
        p2 = dna.translate(table=2, stop_symbol="+")
        assert p2 == "SK+CI+M"
        r = s.replace("T", "U")
        rna = Seq.Seq(r)
        protein = rna.translate(to_stop=True)
        assert protein == "SKRCI"
        gapped_protein = rna.translate()
        assert gapped_protein == "SKRCI*M"

    def test_ambiguous(self):
        s = "RATGATTARAATYTA"
        dna = Seq.Seq(s)
        protein = dna.translate("Vertebrate Mitochondrial")
        assert protein == "BD*NL"
        stop_protein = dna.translate("SGC1", to_stop=True)
        assert stop_protein == "BD"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
