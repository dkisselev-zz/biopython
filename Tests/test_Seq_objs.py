# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Unittests for the Seq objects."""

import array
import unittest
import pytest
import warnings

from Bio import BiopythonWarning
from Bio import SeqIO
from Bio.Data.CodonTable import CodonTable
from Bio.Data.CodonTable import TranslationError
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Data.IUPACData import ambiguous_rna_values
from Bio.Seq import _UndefinedSequenceData
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
from Bio.Seq import translate
from Bio.Seq import UndefinedSequenceError
from Bio.SeqRecord import SeqRecord

try:
    import numpy as np
except ImportError:
    np = None

# This is just the standard table with less stop codons
# (replaced with coding for O as an artificial example)

# Turn black code style off
# fmt: off

special_table = CodonTable(forward_table={
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "O",
    "TGT": "C", "TGC": "C", "TGA": "O", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"},
    start_codons=["TAA", "TAG", "TGA"],
    stop_codons=["TAG"],
)

Chilodonella_uncinata_table = CodonTable(forward_table={
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y",             "TAG": "Q",  # noqa: E241
    "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"},
    start_codons=["ATG"],
    stop_codons=["TAA"],
)

# Turn black code style on
# fmt: on


class StringMethodTests(unittest.TestCase):
    _examples = [
        # These are length 9, a multiple of 3 for translation tests:
        Seq("ACGTGGGGT"),
        Seq("ACGUGGGGU"),
        Seq("GG"),
        Seq("A"),
    ]
    for seq in _examples[:]:
        _examples.append(MutableSeq(seq))
    _start_end_values = [0, 1, 2, 1000, -1, -2, -999, None]
    if np is not None:
        # test with numpy integers (np.int32, np.int64 etc.)
        _start_end_values.extend(np.array([3, 5]))

    def _test_method(self, method_name, start_end=False):
        """Check this method matches the plain string's method."""
        assert isinstance(method_name, str)
        for example1 in self._examples:
            if not hasattr(example1, method_name):
                # e.g. MutableSeq does not support transcribe
                continue
            str1 = str(example1)

            for example2 in self._examples:
                if not hasattr(example2, method_name):
                    # e.g. MutableSeq does not support transcribe
                    continue
                str2 = str(example2)

                try:
                    i = getattr(example1, method_name)(str2)
                except ValueError:
                    i = ValueError
                try:
                    j = getattr(str1, method_name)(str2)
                except ValueError:
                    j = ValueError
                assert i == j, f"{example1!r}.{method_name}({str2!r})"
                try:
                    i = getattr(example1, method_name)(example2)
                except ValueError:
                    i = ValueError
                try:
                    j = getattr(str1, method_name)(str2)
                except ValueError:
                    j = ValueError
                assert i == j, f"{example1!r}.{method_name}({example2!r})"

                if start_end:
                    for start in self._start_end_values:
                        try:
                            i = getattr(example1, method_name)(str2, start)
                        except ValueError:
                            i = ValueError
                        try:
                            j = getattr(str1, method_name)(str2, start)
                        except ValueError:
                            j = ValueError
                        assert i == j, f"{example1!r}.{method_name}({str2!r}, {start})"

                        for end in self._start_end_values:
                            try:
                                i = getattr(example1, method_name)(str2, start, end)
                            except ValueError:
                                i = ValueError
                            try:
                                j = getattr(str1, method_name)(str2, start, end)
                            except ValueError:
                                j = ValueError
                            assert i == j, f"{example1!r}.{method_name}({str2!r}, {start}, {end})"

    def test_str_count(self):
        """Check matches the python string count method."""
        self._test_method("count", start_end=True)
        assert Seq("AC777GT").count("7") == 3
        with pytest.raises(TypeError):
            Seq("AC777GT").count(7)
        with pytest.raises(TypeError):
            Seq("AC777GT").count(None)

    def test_count_overlap(self):
        """Check count_overlap exception matches python string count method."""
        assert Seq("AC777GT").count("77") == 1
        assert Seq("AC777GT").count_overlap("77") == 2
        assert Seq("AC777GT").count_overlap("7") == 3
        with pytest.raises(TypeError):
            Seq("AC777GT").count_overlap(7)
        with pytest.raises(TypeError):
            Seq("AC777GT").count_overlap(None)

    def test_str_count_overlap_GG(self):
        """Check our count_overlap method using GG."""
        # Testing with self._examples
        expected = [
            3,
            3,
            1,
            0,  # Seq() Tests
            3,
            3,
            1,
            0,  # MutableSeq() Tests
        ]

        assert len(self._examples) == len(expected)

        for seq, exp in zip(self._examples, expected):
            # Using search term GG as a string
            assert seq.count_overlap("GG") == exp
            assert seq.count_overlap("G" * 5) == 0
            # Using search term GG as a Seq
            assert seq.count_overlap(Seq("GG")) == exp
            assert seq.count_overlap(Seq("G" * 5)) == 0
            # Using search term GG as a MutableSeq
            assert seq.count_overlap(MutableSeq("GG")) == exp
            assert seq.count_overlap(MutableSeq("G" * 5)) == 0

    def test_count_overlap_start_end_GG(self):
        """Check our count_overlap method using GG with variable ends and starts."""
        # Testing Seq() and MutableSeq() with variable start and end arguments
        start_end_exp = [
            (1, 7, 3),
            (3, None, 3),
            (3, 6, 2),
            (4, 6, 1),
            (4, -1, 2),
            (-5, None, 2),
            (-5, 7, 2),
            (7, -5, 0),
            (-100, None, 3),
            (None, 100, 3),
            (-100, 1000, 3),
        ]

        testing_seq = "GTAGGGGAG"

        for start, end, exp in start_end_exp:
            assert Seq(testing_seq).count_overlap("GG", start, end) == exp
            assert MutableSeq(testing_seq).count_overlap("GG", start, end) == exp

        # Testing Seq() and MutableSeq() with a more heterogeneous sequenece
        assert Seq("GGGTGGTAGGG").count_overlap("GG") == 5
        assert MutableSeq("GGGTGGTAGGG").count_overlap("GG") == 5
        assert Seq("GGGTGGTAGGG").count_overlap("GG", 2, 8) == 1
        assert MutableSeq("GGGTGGTAGGG").count_overlap("GG", 2, 8) == 1
        assert Seq("GGGTGGTAGGG").count_overlap("GG", -11, 6) == 3
        assert MutableSeq("GGGTGGTAGGG").count_overlap("GG", -11, 6) == 3
        assert Seq("GGGTGGTAGGG").count_overlap("GG", 7, 2) == 0
        assert MutableSeq("GGGTGGTAGGG").count_overlap("GG", 7, 2) == 0
        assert Seq("GGGTGGTAGGG").count_overlap("GG", -2, -10) == 0

    def test_str_count_overlap_NN(self):
        """Check our count_overlap method using NN."""
        # Testing with self._examples
        for seq in self._examples:
            # Using search term NN as a string
            assert seq.count_overlap("NN") == 0
            assert seq.count_overlap("N" * 13) == 0
            # Using search term NN as a Seq
            assert seq.count_overlap(Seq("NN")) == 0
            assert seq.count_overlap(Seq("N" * 13)) == 0

    def test_count_overlap_start_end_NN(self):
        """Check our count_overlap method using NN with variable ends and starts."""
        # Testing Seq() and MutableSeq() with variable start and end arguments
        start_end_exp = [
            (1, 7, 0),
            (3, None, 0),
            (3, 6, 0),
            (4, 6, 0),
            (4, -1, 0),
            (-5, None, 0),
            (-5, 7, 0),
            (7, -5, 0),
            (-100, None, 0),
            (None, 100, 0),
            (-100, 1000, 0),
        ]

        testing_seq = "GTAGGGGAG"

        for start, end, exp in start_end_exp:
            assert Seq(testing_seq).count_overlap("NN", start, end) == exp
            assert MutableSeq(testing_seq).count_overlap("NN", start, end) == exp

        # Testing Seq() and MutableSeq() with a more heterogeneous sequenece
        assert Seq("GGGTGGTAGGG").count_overlap("NN") == 0
        assert MutableSeq("GGGTGGTAGGG").count_overlap("NN") == 0
        assert Seq("GGGTGGTAGGG").count_overlap("NN", 2, 8) == 0
        assert MutableSeq("GGGTGGTAGGG").count_overlap("NN", 2, 8) == 0
        assert Seq("GGGTGGTAGGG").count_overlap("NN", -11, 6) == 0
        assert MutableSeq("GGGTGGTAGGG").count_overlap("NN", -11, 6) == 0
        assert Seq("GGGTGGTAGGG").count_overlap("NN", 7, 2) == 0
        assert MutableSeq("GGGTGGTAGGG").count_overlap("NN", 7, 2) == 0
        assert Seq("GGGTGGTAGGG").count_overlap("NN", -10, -2) == 0

    def test_str_find(self):
        """Check matches the python string find method."""
        self._test_method("find", start_end=True)
        assert Seq("AC7GT").find("7") == 2
        with pytest.raises(TypeError):
            Seq("AC7GT").find(7)
        with pytest.raises(TypeError):
            Seq("ACGT").find(None)

    def test_str_rfind(self):
        """Check matches the python string rfind method."""
        self._test_method("rfind", start_end=True)
        assert Seq("AC7GT").rfind("7") == 2
        with pytest.raises(TypeError):
            Seq("AC7GT").rfind(7)
        with pytest.raises(TypeError):
            Seq("ACGT").rfind(None)

    def test_str_index(self):
        """Check matches the python string index method."""
        self._test_method("index", start_end=True)
        assert Seq("AC7GT").index("7") == 2
        with pytest.raises(TypeError):
            Seq("AC7GT").index(7)
        with pytest.raises(TypeError):
            Seq("ACGT").index(None)
        assert MutableSeq("AC7GT").index("7") == 2
        with pytest.raises(TypeError):
            MutableSeq("AC7GT").index(7)
        with pytest.raises(TypeError):
            MutableSeq("ACGT").index(None)

    def test_str_rindex(self):
        """Check matches the python string rindex method."""
        self._test_method("rindex", start_end=True)
        assert Seq("AC7GT").rindex("7") == 2
        with pytest.raises(TypeError):
            Seq("AC7GT").rindex(7)
        with pytest.raises(TypeError):
            Seq("ACGT").rindex(None)
        assert MutableSeq("AC7GT").rindex("7") == 2
        with pytest.raises(TypeError):
            MutableSeq("AC7GT").rindex(7)
        with pytest.raises(TypeError):
            MutableSeq("ACGT").rindex(None)

    def test_str_startswith(self):
        """Check matches the python string startswith method."""
        self._test_method("startswith", start_end=True)
        assert "ABCDE".startswith(("ABE", "OBE", "ABC"))
        with pytest.raises(TypeError):
            Seq("ACGT").startswith(None)
        with pytest.raises(TypeError):
            MutableSeq("ACGT").startswith(None)

        # Now check with a tuple of sub sequences
        for example1 in self._examples:
            subs = tuple(
                example1[start : start + 2] for start in range(0, len(example1) - 2, 3)
            )
            subs_str = tuple(str(s) for s in subs)

            assert str(example1).startswith(subs_str) == example1.startswith(subs)
            assert str(example1).startswith(subs_str) == example1.startswith(subs_str)  # strings!
            assert str(example1).startswith(subs_str, 3) == example1.startswith(subs, 3)
            assert str(example1).startswith(subs_str, 2, 6) == example1.startswith(subs, 2, 6)

    def test_str_endswith(self):
        """Check matches the python string endswith method."""
        self._test_method("endswith", start_end=True)
        assert "ABCDE".endswith(("ABE", "OBE", "CDE"))
        with pytest.raises(TypeError):
            Seq("ACGT").endswith(None)

        # Now check with a tuple of sub sequences
        for example1 in self._examples:
            subs = tuple(
                example1[start : start + 2] for start in range(0, len(example1) - 2, 3)
            )
            subs_str = tuple(str(s) for s in subs)

            assert str(example1).endswith(subs_str) == example1.endswith(subs)
            assert str(example1).startswith(subs_str) == example1.startswith(subs_str)  # strings!
            assert str(example1).endswith(subs_str, 3) == example1.endswith(subs, 3)
            assert str(example1).endswith(subs_str, 2, 6) == example1.endswith(subs, 2, 6)

    def test_str_strip(self):
        """Check matches the python string strip method."""
        self._test_method("strip")
        s = Seq(" ACGT ")
        m = MutableSeq(" ACGT ")
        assert s.strip() == "ACGT"
        with pytest.raises(TypeError):
            s.strip(7)
        assert s == " ACGT "
        assert m.strip() == "ACGT"
        with pytest.raises(TypeError):
            m.strip(7)
        assert m == " ACGT "
        assert m.strip(inplace=True) == "ACGT"
        assert m == "ACGT"
        with pytest.raises(TypeError) as cm:
            s.strip("ACGT", inplace=True)
        assert str(cm.value) == "Sequence is immutable"

    def test_str_lstrip(self):
        """Check matches the python string lstrip method."""
        self._test_method("lstrip")
        s = Seq(" ACGT ")
        m = MutableSeq(" ACGT ")
        assert s.lstrip() == "ACGT "
        with pytest.raises(TypeError):
            s.lstrip(7)
        assert s == " ACGT "
        assert m.lstrip() == "ACGT "
        with pytest.raises(TypeError):
            m.lstrip(7)
        assert m == " ACGT "
        assert m.lstrip(inplace=True) == "ACGT "
        assert m == "ACGT "
        with pytest.raises(TypeError) as cm:
            s.lstrip("ACGT", inplace=True)
        assert str(cm.value) == "Sequence is immutable"

    def test_str_removeprefix(self):
        """Check matches the python string removeprefix method."""
        with pytest.raises(TypeError):
            Seq("ACGT").removeprefix(7)
        with pytest.raises(TypeError):
            Seq("ACGT").removeprefix(7.0)
        with pytest.raises(TypeError):
            Seq("ACGT").removeprefix(None)
        with pytest.raises(UndefinedSequenceError):
            Seq(data=None, length=1000).removeprefix(7.0)
        with pytest.raises(UndefinedSequenceError):
            Seq(data=None, length=1000).removeprefix("A")

    def test_str_removesuffix(self):
        """Check matches the python string removesuffix method."""
        with pytest.raises(TypeError):
            Seq("ACGT").removesuffix(7)
        with pytest.raises(TypeError):
            Seq("ACGT").removesuffix(7.0)
        with pytest.raises(TypeError):
            Seq("ACGT").removesuffix(None)
        with pytest.raises(UndefinedSequenceError):
            Seq(data=None, length=1000).removesuffix(7.0)
        with pytest.raises(UndefinedSequenceError):
            Seq(data=None, length=1000).removesuffix("A")

    def test_str_rstrip(self):
        """Check matches the python string rstrip method."""
        self._test_method("rstrip")
        s = Seq(" ACGT ")
        m = MutableSeq(" ACGT ")
        assert s.rstrip() == " ACGT"
        with pytest.raises(TypeError):
            s.rstrip(7)
        assert s == " ACGT "
        assert m.rstrip() == " ACGT"
        with pytest.raises(TypeError):
            m.rstrip(7)
        assert m == " ACGT "
        assert m.rstrip(inplace=True) == " ACGT"
        assert m == " ACGT"
        with pytest.raises(TypeError) as cm:
            s.rstrip("ACGT", inplace=True)
        assert str(cm.value) == "Sequence is immutable"

    def test_str_split(self):
        """Check matches the python string split method."""
        self._test_method("split")
        assert Seq("AC7GT").split("7") == "AC7GT".split("7")
        with pytest.raises(TypeError):
            Seq("AC7GT").split(7)
        assert MutableSeq("AC7GT").split("7") == "AC7GT".split("7")
        with pytest.raises(TypeError):
            MutableSeq("AC7GT").split(7)

    def test_str_rsplit(self):
        """Check matches the python string rsplit method."""
        self._test_method("rsplit")
        assert Seq("AC7GT").rsplit("7") == "AC7GT".rsplit("7")
        with pytest.raises(TypeError):
            Seq("AC7GT").rsplit(7)
        assert MutableSeq("AC7GT").rsplit("7") == "AC7GT".rsplit("7")
        with pytest.raises(TypeError):
            MutableSeq("AC7GT").rsplit(7)

    def test_str_length(self):
        """Check matches the python string __len__ method."""
        for example1 in self._examples:
            str1 = str(example1)
            assert len(example1) == len(str1)

    def test_str_upper(self):
        """Check matches the python string upper method."""
        for example1 in self._examples:
            str1 = str(example1)
            example1 = example1.upper()
            assert example1 == str1.upper()
        with pytest.raises(TypeError) as cm:
            Seq("abcd").upper(inplace=True)
        assert str(cm.value) == "Sequence is immutable"

    def test_str_lower(self):
        """Check matches the python string lower method."""
        for example1 in self._examples:
            str1 = str(example1)
            example1 = example1.lower()
            assert example1 == str1.lower()
        with pytest.raises(TypeError) as cm:
            Seq("ABCD").lower(inplace=True)
        assert str(cm.value) == "Sequence is immutable"

    def test_str_isupper(self):
        """Check matches the python string isupper method."""
        for example1 in self._examples:
            str1 = str(example1)
            if isinstance(example1, _UndefinedSequenceData):
                with pytest.raises(UndefinedSequenceError):
                    example1.isupper()
            else:
                example1 = example1.isupper()
            assert example1 == str1.isupper()

    def test_str_islower(self):
        """Check matches the python string islower method."""
        for example1 in self._examples:
            str1 = str(example1)
            if isinstance(example1, _UndefinedSequenceData):
                with pytest.raises(UndefinedSequenceError):
                    example1.islower()
            else:
                example1 = example1.islower()
            assert example1 == str1.islower()

    def test_str_replace(self):
        """Check matches the python string replace method."""
        s = Seq("AAGTACGT")
        m = MutableSeq("AAGTACGT")
        t = s.replace("AC", "XYZ")
        assert s == "AAGTACGT"
        assert t == "AAGTXYZGT"
        with pytest.raises(TypeError):
            s.replace("AC", "XYZ", True)
        t = m.replace("AC", "XYZ")
        assert m == "AAGTACGT"
        assert t == "AAGTXYZGT"
        t = m.replace("AC", "XYZ", inplace=True)
        assert m == "AAGTXYZGT"
        assert t == "AAGTXYZGT"
        u = Seq(None, length=20)
        t = u.replace("AT", "CG")
        assert repr(t) == "Seq(None, length=20)"
        with pytest.raises(UndefinedSequenceError) as cm:
            u.replace("AT", "ACGT")  # unequal length
        assert str(cm.value) == "Sequence content is undefined"
        records = SeqIO.parse("TwoBit/sequence.littleendian.2bit", "twobit")
        v = records["seq6"].seq  # ACGTacgtNNNNnn, lazy-loaded
        s = Seq("xyzACGTacgtNNNNnnXYZ")
        t = s.replace(v, "KLM")
        assert t == "xyzKLMXYZ"
        s = Seq("xyzKLMabcd")
        t = s.replace("KLM", v)
        assert t == "xyzACGTacgtNNNNnnabcd"

    def test_str_encode(self):
        """Check matches the python string encode method."""
        for example1 in self._examples:
            str1 = str(example1)
            assert bytes(example1) == str1.encode("ascii")

    def test_str_hash(self):
        for example1 in self._examples:
            if isinstance(example1, MutableSeq):
                continue
            with warnings.catch_warnings():
                # Silence change in behaviour warning
                warnings.simplefilter("ignore", BiopythonWarning)
                assert hash(str(example1)) == hash(example1), ("Hash mismatch, %r for %r vs %r for %r"
                    % (hash(str(example1)), id(example1), hash(example1), example1))

    def test_str_comparison(self):
        for example1 in self._examples:
            for example2 in self._examples:
                with warnings.catch_warnings():
                    assert str(example1) == str(example2) == example1 == example2, f"Checking {example1!r} == {example2!r}"
                    assert str(example1) != str(example2) == example1 != example2, f"Checking {example1!r} != {example2!r}"
                    assert str(example1) < str(example2) == example1 < example2, f"Checking {example1!r} < {example2!r}"
                    assert str(example1) <= str(example2) == example1 <= example2, f"Checking {example1!r} <= {example2!r}"
                    assert str(example1) > str(example2) == example1 > example2, f"Checking {example1!r} > {example2!r}"
                    assert str(example1) >= str(example2) == example1 >= example2, f"Checking {example1!r} >= {example2!r}"

    def test_str_getitem(self):
        """Check slicing and indexing works like a string."""
        for example1 in self._examples:
            str1 = str(example1)
            for i in self._start_end_values:
                if i is not None and abs(i) < len(example1):
                    assert example1[i] == str1[i]
                assert example1[:i] == str1[:i]
                assert example1[i:] == str1[i:]
                for j in self._start_end_values:
                    assert example1[i:j] == str1[i:j]
                    for step in range(-3, 4):
                        if step == 0:
                            with pytest.raises(ValueError) as cm:
                                example1[i:j:step]
                            assert str(cm.value) == "slice step cannot be zero"
                        else:
                            assert example1[i:j:step] == str1[i:j:step]
        u = Seq(None, length=0)
        assert u == ""

    def test_search(self):
        """Check the search method of Seq objects."""
        s = Seq("ACGTACGT")
        matches = s.search(["CGT", Seq("CG"), b"ACGT", bytearray(b"GTA")])
        assert next(matches) == (0, "ACGT")
        assert next(matches) == (1, "CGT")
        assert next(matches) == (1, "CG")
        assert next(matches) == (2, "GTA")
        assert next(matches) == (4, "ACGT")
        assert next(matches) == (5, "CGT")
        assert next(matches) == (5, "CG")
        with pytest.raises(StopIteration):
            next(matches)
        s = MutableSeq("ACGTACGT")
        matches = s.search(["CGT", Seq("CG"), b"ACGT", bytearray(b"GTA")])
        assert next(matches) == (0, "ACGT")
        assert next(matches) == (1, "CGT")
        assert next(matches) == (1, "CG")
        assert next(matches) == (2, "GTA")
        assert next(matches) == (4, "ACGT")
        assert next(matches) == (5, "CGT")
        assert next(matches) == (5, "CG")
        with pytest.raises(StopIteration):
            next(matches)

    def test_MutableSeq_setitem(self):
        """Check setting sequence contents of a MutableSeq object."""
        m = MutableSeq("ABCD")
        m[1] = "X"
        assert m == "AXCD"
        m[1:3] = MutableSeq("XY")
        assert m == "AXYD"
        m[1:3] = Seq("KL")
        assert m == "AKLD"
        m[1:3] = Seq("bc")
        assert m == "AbcD"
        with pytest.raises(TypeError) as cm:
            m[1:3] = 9
        assert str(cm.value) == "received unexpected type 'int'"

    def test_MutableSeq_extend(self):
        """Check extending a MutableSeq object."""
        m = MutableSeq("ABCD")
        m.extend(MutableSeq("xyz"))
        assert m == "ABCDxyz"
        m.extend(Seq("KLM"))
        assert m == "ABCDxyzKLM"
        m.extend("PQRST")
        assert m == "ABCDxyzKLMPQRST"
        with pytest.raises(TypeError) as cm:
            m.extend(5)
        assert str(cm.value) == "expected a string, Seq or MutableSeq"

    def test_tomutable(self):
        """Check creating a MutableSeq object."""
        for example1 in self._examples:
            mut = MutableSeq(example1)
            assert isinstance(mut, MutableSeq)
            assert mut == example1

    def test_toseq(self):
        """Check creating a Seq object."""
        for example1 in self._examples:
            seq = Seq(example1)
            assert isinstance(seq, Seq)
            assert seq == example1

    def test_the_complement(self):
        """Check obj.complement() method."""
        mapping = ""
        for example1 in self._examples:
            if isinstance(example1, MutableSeq):
                continue
            if "u" in example1 or "U" in example1:
                comp = example1.complement_rna()
            else:
                comp = example1.complement()
            str1 = str(example1)
            if "U" in str1 or "u" in str1:
                mapping = str.maketrans("ACGUacgu", "UGCAugca")
            else:
                # Default to DNA, e.g. complement("A") -> "T" not "U"
                mapping = str.maketrans("ACGTacgt", "TGCAtgca")
            assert str1.translate(mapping) == comp
        with pytest.raises(TypeError) as cm:
            Seq("ACGT").complement(inplace=True)
        assert str(cm.value) == "Sequence is immutable"
        with pytest.raises(TypeError) as cm:
            Seq("ACGT").complement_rna(inplace=True)
        assert str(cm.value) == "Sequence is immutable"
        assert repr(Seq(None, length=20).complement()) == "Seq(None, length=20)"
        assert repr(Seq(None, length=20).complement_rna()) == "Seq(None, length=20)"

    def test_the_reverse_complement(self):
        """Check obj.reverse_complement() method."""
        mapping = ""
        for example1 in self._examples:
            if isinstance(example1, MutableSeq):
                continue
            if "u" in example1 or "U" in example1:
                comp = example1.reverse_complement_rna()
            else:
                comp = example1.reverse_complement()
            str1 = str(example1)
            if "U" in str1 or "u" in str1:
                mapping = str.maketrans("ACGUacgu", "UGCAugca")
            else:
                # Defaults to DNA, so reverse_complement("A") --> "T" not "U"
                mapping = str.maketrans("ACGTacgt", "TGCAtgca")
            assert str1.translate(mapping)[::-1] == comp
        with pytest.raises(TypeError) as cm:
            Seq("ACGT").reverse_complement(inplace=True)
        assert str(cm.value) == "Sequence is immutable"
        with pytest.raises(TypeError) as cm:
            Seq("ACGT").reverse_complement_rna(inplace=True)
        assert str(cm.value) == "Sequence is immutable"
        assert repr(Seq(None, length=20).reverse_complement()) == "Seq(None, length=20)"
        assert repr(Seq(None, length=20).reverse_complement_rna()) == "Seq(None, length=20)"

    def test_the_transcription(self):
        """Check obj.transcribe() method."""
        mapping = ""
        for example1 in self._examples:
            if isinstance(example1, MutableSeq):
                continue
            tran = example1.transcribe()
            str1 = str(example1)
            if len(str1) % 3 != 0:
                # TODO - Check for or silence the expected warning?
                continue
            assert str1.replace("T", "U").replace("t", "u") == tran
        with pytest.raises(TypeError) as cm:
            Seq("ACGT").transcribe(inplace=True)
        assert str(cm.value) == "Sequence is immutable"
        assert repr(Seq(None, length=20).transcribe()) == "Seq(None, length=20)"

    def test_the_back_transcription(self):
        """Check obj.back_transcribe() method."""
        mapping = ""
        for example1 in self._examples:
            if isinstance(example1, MutableSeq):
                continue
            tran = example1.back_transcribe()
            str1 = str(example1)
            assert str1.replace("U", "T").replace("u", "t") == tran
        with pytest.raises(TypeError) as cm:
            Seq("ACGU").back_transcribe(inplace=True)
        assert str(cm.value) == "Sequence is immutable"
        assert repr(Seq(None, length=20).back_transcribe()) == "Seq(None, length=20)"

    def test_the_translate(self):
        """Check obj.translate() method."""
        mapping = ""
        for example1 in self._examples:
            if len(example1) % 3 != 0:
                # TODO - Check for or silence the expected warning?
                continue
            tran = example1.translate()
            # Try with positional vs named argument:
            assert example1.translate(11) == example1.translate(table=11)

            # TODO - check the actual translation, and all the optional args
        with pytest.raises(ValueError):
            Seq("ABCD").translate("XYZ")
        with pytest.warns(BiopythonWarning) as cm:
            Seq(None, length=20).translate()
        assert str(cm.warning) == "Partial codon, len(sequence) not a multiple of three. This may become an error in future."

    def test_the_translation_of_stops(self):
        """Check obj.translate() method with stop codons."""
        misc_stops = "TAATAGTGAAGAAGG"
        nuc = Seq(misc_stops)
        assert "***RR" == nuc.translate()
        assert "***RR" == nuc.translate(1)
        assert "***RR" == nuc.translate("SGC0")
        assert "**W**" == nuc.translate(table=2)
        assert "**WRR" == nuc.translate(table="Yeast Mitochondrial")
        assert "**WSS" == nuc.translate(table=5)
        assert "**WSS" == nuc.translate(table=9)
        assert "**CRR" == nuc.translate(table="Euplotid Nuclear")
        assert "***RR" == nuc.translate(table=11)
        assert "***RR" == nuc.translate(table="11")
        assert "***RR" == nuc.translate(table="Bacterial")
        assert "**GRR" == nuc.translate(table=25)
        assert "" == nuc.translate(to_stop=True)
        assert "O*ORR" == nuc.translate(table=special_table)
        assert "*QWRR" == nuc.translate(table=Chilodonella_uncinata_table)
        nuc = MutableSeq(misc_stops)
        assert "***RR" == nuc.translate()
        assert "***RR" == nuc.translate(1)
        assert "***RR" == nuc.translate("SGC0")
        assert "**W**" == nuc.translate(table=2)
        assert "**WRR" == nuc.translate(table="Yeast Mitochondrial")
        assert "**WSS" == nuc.translate(table=5)
        assert "**WSS" == nuc.translate(table=9)
        assert "**CRR" == nuc.translate(table="Euplotid Nuclear")
        assert "***RR" == nuc.translate(table=11)
        assert "***RR" == nuc.translate(table="11")
        assert "***RR" == nuc.translate(table="Bacterial")
        assert "**GRR" == nuc.translate(table=25)
        assert "" == nuc.translate(to_stop=True)
        assert "O*ORR" == nuc.translate(table=special_table)
        assert "*QWRR" == nuc.translate(table=Chilodonella_uncinata_table)
        # These test the Bio.Seq.translate() function - move these?:
        assert "*QWRR" == translate(str(nuc), table=Chilodonella_uncinata_table)
        assert "O*ORR" == translate(str(nuc), table=special_table)
        assert "" == translate(str(nuc), to_stop=True)
        assert "***RR" == translate(str(nuc), table="Bacterial")
        assert "***RR" == translate(str(nuc), table="11")
        assert "***RR" == translate(str(nuc), table=11)
        assert "**W**" == translate(str(nuc), table=2)
        assert Seq("TAT").translate() == "Y"
        assert Seq("TAR").translate() == "*"
        assert Seq("TAN").translate() == "X"
        assert Seq("NNN").translate() == "X"
        assert Seq("TAt").translate() == "Y"
        assert Seq("TaR").translate() == "*"
        assert Seq("TaN").translate() == "X"
        assert Seq("nnN").translate() == "X"
        assert Seq("tat").translate() == "Y"
        assert Seq("tar").translate() == "*"
        assert Seq("tan").translate() == "X"
        assert Seq("nnn").translate() == "X"

    def test_the_translation_of_invalid_codons(self):
        """Check obj.translate() method with invalid codons."""
        for codon in ["TA?", "N-N", "AC_", "Ac_"]:
            msg = f"Translating {codon} should fail"
            nuc = Seq(codon)
            with pytest.raises(TranslationError):
                nuc.translate()
            nuc = MutableSeq(codon)
            with pytest.raises(TranslationError):
                nuc.translate()

    def test_the_translation_of_ambig_codons(self):
        """Check obj.translate() method with ambiguous codons."""
        for ambig_values in [ambiguous_dna_values, ambiguous_rna_values]:
            ambig = set(ambig_values.keys())
            ambig.remove("X")
            for c1 in ambig:
                for c2 in ambig:
                    for c3 in ambig:
                        values = {
                            str(Seq(a + b + c).translate())
                            for a in ambig_values[c1]
                            for b in ambig_values[c2]
                            for c in ambig_values[c3]
                        }
                        t = Seq(c1 + c2 + c3).translate()
                        if t == "*":
                            assert values == set("*")
                        elif t == "X":
                            assert len(values) > 1, ("translate('%s') = '%s' not '%s'"
                                % (c1 + c2 + c3, t, ",".join(values)))
                        elif t == "Z":
                            assert values == set("EQ")
                        elif t == "B":
                            assert values == set("DN")
                        elif t == "J":
                            assert values == set("LI")
                        else:
                            assert values == set(t)
                        # TODO - Use the Bio.Data.IUPACData module for the
                        # ambiguous protein mappings?

    def test_Seq_init_error(self):
        """Check Seq __init__ raises the appropriate exceptions."""
        with pytest.raises(TypeError):
            Seq(("A", "C", "G", "T"))
        with pytest.raises(TypeError):
            Seq(["A", "C", "G", "T"])
        with pytest.raises(TypeError):
            Seq(1)
        with pytest.raises(TypeError):
            Seq(1.0)
        with pytest.raises(ValueError):
            Seq(None)

    def test_MutableSeq_init_error(self):
        """Check MutableSeq __init__ raises the appropriate exceptions."""
        with pytest.raises(TypeError):
            MutableSeq(("A", "C", "G", "T"))
        with pytest.raises(TypeError):
            MutableSeq(["A", "C", "G", "T"])
        with pytest.raises(TypeError):
            MutableSeq(1)
        with pytest.raises(TypeError):
            MutableSeq(1.0)
        with pytest.raises(TypeError):
            MutableSeq(array.array("i", [1, 2, 3, 4]))

    def test_join_Seq_TypeError(self):
        """Checks that a TypeError is thrown for all non-iterable types."""
        # No iterable types which contain non-accepted types either.

        spacer = Seq("NNNNN")
        with pytest.raises(TypeError):
            spacer.join(5)
        with pytest.raises(TypeError):
            spacer.join(SeqRecord(Seq("ATG")))
        with pytest.raises(TypeError):
            spacer.join(["ATG", "ATG", 5, "ATG"])

    def test_join_MutableSeq_TypeError_iter(self):
        """Checks that a TypeError is thrown for all non-iterable types."""
        # No iterable types which contain non-accepted types either.

        spacer = MutableSeq("MMMMM")
        with pytest.raises(TypeError):
            spacer.join(5)
        with pytest.raises(TypeError):
            spacer.join(["ATG", "ATG", 5, "ATG"])

    def test_join_Seq(self):
        """Checks if Seq join correctly concatenates sequence with the spacer."""
        spacer = Seq("NNNNN")
        assert "N" * 15 == spacer.join([Seq("NNNNN"), Seq("NNNNN")])

        spacer1 = Seq("")
        spacers = [spacer1, Seq("NNNNN"), Seq("GGG")]
        example_strings = ["ATG", "ATG", "ATG", "ATG"]
        example_strings_seqs = ["ATG", "ATG", Seq("ATG"), "ATG"]

        # strings with empty spacer
        str_concatenated = spacer1.join(example_strings)

        assert str_concatenated == "".join(example_strings)

        for spacer in spacers:
            seq_concatenated = spacer.join(example_strings_seqs)
            assert seq_concatenated == str(spacer).join(example_strings)
            # Now try single sequence arguments, should join the letters
            for target in example_strings + example_strings_seqs:
                assert str(spacer).join(str(target)) == str(spacer.join(target))

    def test_join_MutableSeq_mixed(self):
        """Check MutableSeq objects can be joined."""
        spacer = MutableSeq("NNNNN")
        assert "N" * 15 == spacer.join([MutableSeq("NNNNN"), MutableSeq("NNNNN")])
        with pytest.raises(TypeError):
            spacer.join([Seq("NNNNN"), MutableSeq("NNNNN")])()

    def test_join_Seq_with_file(self):
        """Checks if Seq join correctly concatenates sequence from a file with the spacer."""
        filename = "Fasta/f003.fa"
        seqlist = [record.seq for record in SeqIO.parse(filename, "fasta")]
        seqlist_as_strings = [str(_) for _ in seqlist]

        spacer = Seq("NNNNN")
        spacer1 = Seq("")
        # seq objects with spacer
        seq_concatenated = spacer.join(seqlist)
        # seq objects with empty spacer
        seq_concatenated1 = spacer1.join(seqlist)

        ref_data = ref_data1 = ""
        ref_data = str(spacer).join(seqlist_as_strings)
        ref_data1 = str(spacer1).join(seqlist_as_strings)

        assert seq_concatenated == ref_data
        assert seq_concatenated1 == ref_data1
        with pytest.raises(TypeError):
            spacer.join(SeqIO.parse(filename, "fasta"))

    def test_join_MutableSeq(self):
        """Checks if MutableSeq join correctly concatenates sequence with the spacer."""
        # Only expect it to take Seq objects and/or strings in an iterable!

        spacer1 = MutableSeq("")
        spacers = [spacer1, MutableSeq("NNNNN"), MutableSeq("GGG")]
        example_strings = ["ATG", "ATG", "ATG", "ATG"]
        example_strings_seqs = ["ATG", "ATG", Seq("ATG"), "ATG"]

        # strings with empty spacer
        str_concatenated = spacer1.join(example_strings)

        assert str_concatenated == "".join(example_strings)

        for spacer in spacers:
            seq_concatenated = spacer.join(example_strings_seqs)
            assert seq_concatenated == str(spacer).join(example_strings)

    def test_join_MutableSeq_with_file(self):
        """Checks if MutableSeq join correctly concatenates sequence from a file with the spacer."""
        filename = "Fasta/f003.fa"
        seqlist = [record.seq for record in SeqIO.parse(filename, "fasta")]
        seqlist_as_strings = [str(_) for _ in seqlist]

        spacer = MutableSeq("NNNNN")
        spacer1 = MutableSeq("")
        # seq objects with spacer
        seq_concatenated = spacer.join(seqlist)
        # seq objects with empty spacer
        seq_concatenated1 = spacer1.join(seqlist)

        ref_data = ref_data1 = ""
        ref_data = str(spacer).join(seqlist_as_strings)
        ref_data1 = str(spacer1).join(seqlist_as_strings)

        assert seq_concatenated == ref_data
        assert seq_concatenated1 == ref_data1
        with pytest.raises(TypeError):
            spacer.join(SeqIO.parse(filename, "fasta"))

    def test_equality(self):
        """Test equality when mixing types."""
        assert Seq("6") == "6"
        assert Seq("6") != 6
        assert Seq("") == ""
        assert Seq("") != None
        assert Seq("None") == "None"
        assert Seq("None") != None

        assert MutableSeq("6") == "6"
        assert MutableSeq("6") != 6
        assert MutableSeq("") == ""
        assert MutableSeq("") != None
        assert MutableSeq("None") == "None"
        assert MutableSeq("None") != None

    # TODO - Addition...


class FileBasedTests(unittest.TestCase):
    """Test Seq objects created from files by SeqIO."""

    def test_unknown_seq(self):
        """Test if feature extraction works properly for unknown sequences."""
        rec = SeqIO.read("GenBank/NT_019265.gb", "genbank")
        assert isinstance(rec.seq, Seq)
        with pytest.raises(UndefinedSequenceError):
            bytes(rec.seq)

        feature = rec.features[1]
        seq = feature.extract(rec.seq)
        assert isinstance(seq, Seq)
        assert len(seq) == len(feature)
        with pytest.raises(UndefinedSequenceError):
            bytes(seq)


class ComparisonTests(unittest.TestCase):
    """Test comparisons of Seq and MutableSeq."""

    def test_eq(self):
        s1 = "GCATGTATGT"
        s2 = "TTGATCAGTT"
        for seq1 in (Seq(s1), MutableSeq(s1)):
            for seq2 in (Seq(s2), MutableSeq(s2)):
                msg = f"{type(seq1)} vs {type(seq2)}"
                assert seq1[2:4] == seq2[3:5], msg
                assert seq1[2:4] == "AT", msg
                assert seq2[3:5] == "AT", msg
                assert seq1[2:4] == b"AT", msg
                assert seq2[3:5] == b"AT", msg
                assert seq1[2:4] == Seq("AT"), msg
                assert seq2[3:5] == Seq("AT"), msg
                assert seq1[2:4] == MutableSeq("AT"), msg
                assert seq2[3:5] == MutableSeq("AT"), msg
                assert seq2[3:5] == seq1[2:4], msg
                assert "AT" == seq1[2:4], msg
                assert "AT" == seq2[3:5], msg
                assert b"AT" == seq1[2:4], msg
                assert b"AT" == seq2[3:5], msg
                assert Seq("AT") == seq1[2:4], msg
                assert Seq("AT") == seq2[3:5], msg
                assert MutableSeq("AT") == seq1[2:4], msg
                assert MutableSeq("AT") == seq2[3:5], msg
                with pytest.raises(UndefinedSequenceError):
                    seq1 == Seq(None, len(seq1))
                with pytest.raises(UndefinedSequenceError):
                    seq2 == Seq(None, len(seq2))
                with pytest.raises(UndefinedSequenceError):
                    seq1 == Seq(None, 10)
                with pytest.raises(UndefinedSequenceError):
                    seq2 == Seq(None, 10)
                with pytest.raises(UndefinedSequenceError):
                    Seq(None, len(seq1)) == seq1
                with pytest.raises(UndefinedSequenceError):
                    Seq(None, len(seq2)) == seq2
                with pytest.raises(UndefinedSequenceError):
                    Seq(None, 10) == seq1
                with pytest.raises(UndefinedSequenceError):
                    Seq(None, 10) == seq2

    def test_ne(self):
        s1 = "GCATGTATGT"
        s2 = "TTGATCAGTT"
        for seq1 in (Seq(s1), MutableSeq(s1)):
            for seq2 in (Seq(s2), MutableSeq(s2)):
                msg = f"{type(seq1)} vs {type(seq2)}"
                assert seq1 != seq2, msg
                assert seq1 != "AT", msg
                assert seq2 != "AT", msg
                assert seq1 != b"AT", msg
                assert seq2 != b"AT", msg
                assert seq1 != Seq("AT"), msg
                assert seq2 != Seq("AT"), msg
                assert seq1 != MutableSeq("AT"), msg
                assert seq2 != MutableSeq("AT"), msg
                assert seq1[2:4] != "CG", msg
                assert seq2[3:5] != "CG", msg
                assert seq1[2:4] != b"CG", msg
                assert seq2[3:5] != b"CG", msg
                assert seq1[2:4] != Seq("CG"), msg
                assert seq2[3:5] != Seq("CG"), msg
                assert seq1[2:4] != MutableSeq("CG"), msg
                assert seq2[3:5] != MutableSeq("CG"), msg
                assert seq2 != seq1, msg
                assert "AT" != seq1, msg
                assert "AT" != seq2, msg
                assert b"AT" != seq1, msg
                assert b"AT" != seq2, msg
                assert Seq("AT") != seq1, msg
                assert Seq("AT") != seq2, msg
                assert MutableSeq("AT") != seq1, msg
                assert MutableSeq("AT") != seq2, msg
                assert "CG" != seq1[2:4], msg
                assert "CG" != seq2[3:5], msg
                assert b"CG" != seq1[2:4], msg
                assert b"CG" != seq2[3:5], msg
                assert Seq("CG") != seq1[2:4], msg
                assert Seq("CG") != seq2[3:5], msg
                assert MutableSeq("CG") != seq1[2:4], msg
                assert MutableSeq("CG") != seq2[3:5], msg
                with pytest.raises(UndefinedSequenceError):
                    seq1 != Seq(None, len(seq1))
                with pytest.raises(UndefinedSequenceError):
                    seq2 != Seq(None, len(seq2))
                with pytest.raises(UndefinedSequenceError):
                    seq1 != Seq(None, 10)
                with pytest.raises(UndefinedSequenceError):
                    seq2 != Seq(None, 10)
                with pytest.raises(UndefinedSequenceError):
                    Seq(None, len(seq1)) != seq1
                with pytest.raises(UndefinedSequenceError):
                    Seq(None, len(seq2)) != seq2
                with pytest.raises(UndefinedSequenceError):
                    Seq(None, 10) != seq1
                with pytest.raises(UndefinedSequenceError):
                    Seq(None, 10) != seq2

    def test_lt(self):
        s1 = "GCATGTATGT"
        s2 = "TTGATCAGTT"
        for seq1 in (Seq(s1), MutableSeq(s1)):
            for seq2 in (Seq(s2), MutableSeq(s2)):
                msg = f"{type(seq1)} vs {type(seq2)}"
                assert seq1 < seq2, msg
                assert "AA" < seq1, msg
                assert seq1 < "TT", msg
                assert "AA" < seq2, msg
                assert seq2 < "TTT", msg
                assert b"AA" < seq1, msg
                assert seq1 < b"TT", msg
                assert b"AA" < seq2, msg
                assert seq2 < b"TTT", msg
                assert Seq("AA") < seq1, msg
                assert seq1 < Seq("TT"), msg
                assert Seq("AA") < seq2, msg
                assert seq2 < Seq("TTT"), msg
                assert MutableSeq("AA") < seq1, msg
                assert seq1 < MutableSeq("TT"), msg
                assert MutableSeq("AA") < seq2, msg
                assert seq2 < MutableSeq("TTT"), msg
                with pytest.raises(UndefinedSequenceError):
                    seq1 < Seq(None, len(seq1))
                with pytest.raises(UndefinedSequenceError):
                    seq2 < Seq(None, len(seq2))
                with pytest.raises(UndefinedSequenceError):
                    seq1 < Seq(None, 10)
                with pytest.raises(UndefinedSequenceError):
                    seq2 < Seq(None, 10)
                assert "AA" < seq1[2:4], msg
                assert "AA" < seq2[3:5], msg
                assert b"AA" < seq1[2:4], msg
                assert b"AA" < seq2[3:5], msg
                assert seq2[3:5] < seq1[2:6], msg
                assert Seq("AA") < seq1[2:4], msg
                assert Seq("AA") < seq2[3:5], msg
                assert MutableSeq("AA") < seq1[2:4], msg
                assert MutableSeq("AA") < seq2[3:5], msg
                assert seq1[2:4] < "TT", msg
                assert seq2[3:5] < "TT", msg
                assert seq1[2:4] < b"TT", msg
                assert seq2[3:5] < b"TT", msg
                assert seq1[2:4] < Seq("TT"), msg
                assert seq2[3:5] < Seq("TT"), msg
                assert seq1[2:4] < MutableSeq("TT"), msg
                assert seq2[3:5] < MutableSeq("TT"), msg

    def test_le(self):
        s1 = "GCATGTATGT"
        s2 = "TTGATCAGTT"
        for seq1 in (Seq(s1), MutableSeq(s1)):
            for seq2 in (Seq(s2), MutableSeq(s2)):
                msg = f"{type(seq1)} vs {type(seq2)}"
                assert seq1 <= seq2, msg
                assert seq1 <= "TT", msg
                assert "TT" <= seq2, msg
                assert seq1 <= b"TT", msg
                assert "TT" <= seq2, msg
                assert seq1 <= Seq("TT"), msg
                assert "TT" <= seq2, msg
                assert seq1 <= MutableSeq("TT"), msg
                assert MutableSeq("TT") <= seq2, msg
                assert "AA" <= seq1, msg
                assert "AA" <= seq2, msg
                assert b"AA" <= seq1, msg
                assert b"AA" <= seq2, msg
                assert Seq("AA") <= seq1, msg
                assert Seq("AA") <= seq2, msg
                assert MutableSeq("AA") <= seq1, msg
                assert MutableSeq("AA") <= seq2, msg
                assert "GC" <= seq1, msg
                assert "GC" <= seq2, msg
                assert b"GC" <= seq1, msg
                assert b"GC" <= seq2, msg
                assert Seq("GC") <= seq1, msg
                assert Seq("GC") <= seq2, msg
                assert MutableSeq("GC") <= seq1, msg
                assert MutableSeq("GC") <= seq2, msg
                with pytest.raises(UndefinedSequenceError):
                    seq1 <= Seq(None, len(seq1))
                with pytest.raises(UndefinedSequenceError):
                    seq2 <= Seq(None, len(seq2))
                with pytest.raises(UndefinedSequenceError):
                    seq1 <= Seq(None, 10)
                with pytest.raises(UndefinedSequenceError):
                    seq2 <= Seq(None, 10)
                assert "AA" <= seq1[2:4], msg
                assert "AA" <= seq2[3:5], msg
                assert b"AA" <= seq1[2:4], msg
                assert b"AA" <= seq2[3:5], msg
                assert seq1[2:4] <= seq2[3:5], msg
                assert Seq("AA") <= seq1[2:4], msg
                assert Seq("AA") <= seq2[3:5], msg
                assert MutableSeq("AA") <= seq1[2:4], msg
                assert MutableSeq("AA") <= seq2[3:5], msg
                assert seq1[2:4] <= "TT", msg
                assert seq2[3:5] <= "TT", msg
                assert seq1[2:4] <= b"TT", msg
                assert seq2[3:5] <= b"TT", msg
                assert seq1[2:4] <= seq2[3:5], msg
                assert seq1[2:4] <= Seq("TT"), msg
                assert seq2[3:5] <= Seq("TT"), msg
                assert seq1[2:4] <= MutableSeq("TT"), msg
                assert seq2[3:5] <= MutableSeq("TT"), msg

    def test_gt(self):
        s1 = "GCATGTATGT"
        s2 = "TTGATCAGTT"
        for seq1 in (Seq(s1), MutableSeq(s1)):
            for seq2 in (Seq(s2), MutableSeq(s2)):
                msg = f"{type(seq1)} vs {type(seq2)}"
                assert seq2 > seq1, msg
                assert "TT" > seq1, msg
                assert seq2 > "TT", msg
                assert b"TT" > seq1, msg
                assert seq2 > b"TT", msg
                assert Seq("TT") > seq1, msg
                assert seq2 > Seq("TT"), msg
                assert MutableSeq("TT") > seq1, msg
                assert seq2 > MutableSeq("TT"), msg
                assert seq1 > "AA", msg
                assert seq2 > "AA", msg
                assert seq1 > b"AA", msg
                assert seq2 > b"AA", msg
                assert seq1 > Seq("AA"), msg
                assert seq2 > Seq("AA"), msg
                assert seq1 > MutableSeq("AA"), msg
                assert seq2 > MutableSeq("AA"), msg
                assert seq1 > "GC", msg
                assert seq2 > "GC", msg
                assert seq1 > b"GC", msg
                assert seq2 > b"GC", msg
                assert seq1 > Seq("GC"), msg
                assert seq2 > Seq("GC"), msg
                assert seq1 > MutableSeq("GC"), msg
                assert seq2 > MutableSeq("GC"), msg
                with pytest.raises(UndefinedSequenceError):
                    seq1 > Seq(None, len(seq1))
                with pytest.raises(UndefinedSequenceError):
                    seq2 > Seq(None, len(seq2))
                with pytest.raises(UndefinedSequenceError):
                    seq1 > Seq(None, 10)
                with pytest.raises(UndefinedSequenceError):
                    seq2 > Seq(None, 10)
                assert seq1[2:4] > "AA", msg
                assert seq2[3:5] > "AA", msg
                assert seq1[2:4] > b"AA", msg
                assert seq2[3:5] > b"AA", msg
                assert seq1[2:6] > seq2[3:5], msg
                assert seq1[2:4] > Seq("AA"), msg
                assert seq2[3:5] > Seq("AA"), msg
                assert seq1[2:4] > MutableSeq("AA"), msg
                assert seq2[3:5] > MutableSeq("AA"), msg
                assert "TT" > seq1[2:4], msg
                assert "TT" > seq2[3:5], msg
                assert b"TT" > seq1[2:4], msg
                assert b"TT" > seq2[3:5], msg
                assert Seq("TT") > seq1[2:4], msg
                assert Seq("TT") > seq2[3:5], msg
                assert MutableSeq("TT") > seq1[2:4], msg
                assert MutableSeq("TT") > seq2[3:5], msg

    def test_ge(self):
        s1 = "GCATGTATGT"
        s2 = "TTGATCAGTT"
        for seq1 in (Seq(s1), MutableSeq(s1)):
            for seq2 in (Seq(s2), MutableSeq(s2)):
                msg = f"{type(seq1)} vs {type(seq2)}"
                assert seq2 >= seq1, msg
                assert "TT" >= seq1, msg
                assert seq2 >= "TT", msg
                assert b"TT" >= seq1, msg
                assert seq2 >= b"TT", msg
                assert Seq("TT") >= seq1, msg
                assert seq2 >= Seq("TT"), msg
                assert MutableSeq("TT") >= seq1, msg
                assert seq2 >= MutableSeq("TT"), msg
                assert seq1 >= "AA", msg
                assert seq2 >= "AA", msg
                assert seq1 >= b"AA", msg
                assert seq2 >= b"AA", msg
                assert seq1 >= Seq("AA"), msg
                assert seq2 >= Seq("AA"), msg
                assert seq1 >= MutableSeq("AA"), msg
                assert seq2 >= MutableSeq("AA"), msg
                assert seq1 >= "GC", msg
                assert seq2 >= "GC", msg
                assert seq1 >= b"GC", msg
                assert seq2 >= b"GC", msg
                assert seq1 >= Seq("GC"), msg
                assert seq2 >= Seq("GC"), msg
                assert seq1 >= MutableSeq("GC"), msg
                assert seq2 >= MutableSeq("GC"), msg
                with pytest.raises(UndefinedSequenceError):
                    seq1 >= Seq(None, len(seq1))
                with pytest.raises(UndefinedSequenceError):
                    seq2 >= Seq(None, len(seq2))
                with pytest.raises(UndefinedSequenceError):
                    seq1 >= Seq(None, 10)
                with pytest.raises(UndefinedSequenceError):
                    seq2 >= Seq(None, 10)
                assert seq1[2:4] >= "AA", msg
                assert seq2[3:5] >= "AA", msg
                assert seq1[2:4] >= b"AA", msg
                assert seq2[3:5] >= b"AA", msg
                assert seq1[2:4] >= seq2[3:5], msg
                assert seq1[2:4] >= Seq("AA"), msg
                assert seq2[3:5] >= Seq("AA"), msg
                assert seq1[2:4] >= MutableSeq("AA"), msg
                assert seq2[3:5] >= MutableSeq("AA"), msg
                assert "TT" >= seq1[2:4], msg
                assert "TT" >= seq2[3:5], msg
                assert b"TT" >= seq1[2:4], msg
                assert b"TT" >= seq2[3:5], msg
                assert seq1[2:4] >= seq2[3:5], msg
                assert Seq("TT") >= seq1[2:4], msg
                assert Seq("TT") >= seq2[3:5], msg
                assert MutableSeq("TT") >= seq1[2:4], msg
                assert MutableSeq("TT") >= seq2[3:5], msg


class PartialSequenceTests(unittest.TestCase):
    """Test Seq objects with partially defined sequences."""

    def test_init(self):
        seq = Seq({5: "ACGT"}, length=20)
        assert repr(seq) == "Seq({5: 'ACGT'}, length=20)"
        with pytest.raises(ValueError) as cm:
            Seq({5: "ACGT"}, length=-10)
        assert str(cm.value) == "length must not be negative."
        with pytest.raises(ValueError) as cm:
            Seq({5: 1.5}, length=10)
        assert str(cm.value) == "Expected bytes-like objects or strings"
        with pytest.raises(ValueError) as cm:
            Seq({5: "PQRST", 8: "KLM"}, length=10)
        assert str(cm.value) == "Sequence data are overlapping."
        with pytest.raises(ValueError) as cm:
            Seq({5: "PQRST"}, length=8)
        assert str(cm.value) == "Provided sequence data extend beyond sequence length."

    def test_repr(self):
        seq = Seq({5: "ACGT", 14: "GGC"}, length=20)
        assert repr(seq) == "Seq({5: 'ACGT', 14: 'GGC'}, length=20)"
        seq = Seq({5: "ACGT" * 25, 140: "GGC"}, length=143)
        assert repr(seq) == "Seq({5: 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC...CGT', 140: 'GGC'}, length=143)"
        seq = Seq({5: "ACGT" * 25, 140: "GGC"}, length=150)
        assert repr(seq) == "Seq({5: 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC...CGT', 140: 'GGC'}, length=150)"
        seq = Seq({5: "ACGT" * 25, 140: "acgt" * 20}, length=250)
        assert repr(seq) == "Seq({5: 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC...CGT', 140: 'acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtac...cgt'}, length=250)"

    def test_getitem(self):
        #           1    1    2
        # 0    5    0    5    0
        # ?????ABCD?????EFG???
        seq = Seq({5: "ABCD", 14: "EFG"}, length=20)
        assert repr(seq) == "Seq({5: 'ABCD', 14: 'EFG'}, length=20)"
        # index
        with pytest.raises(UndefinedSequenceError) as cm:
            seq[0]
        assert str(cm.value) == "Sequence at position 0 is undefined"
        with pytest.raises(UndefinedSequenceError) as cm:
            seq[1]
        assert str(cm.value) == "Sequence at position 1 is undefined"
        assert seq[5] == "A"
        assert seq[6] == "B"
        assert seq[7] == "C"
        assert seq[8] == "D"
        with pytest.raises(UndefinedSequenceError) as cm:
            seq[10]
        assert str(cm.value) == "Sequence at position 10 is undefined"
        assert seq[14] == "E"
        assert seq[15] == "F"
        assert seq[16] == "G"
        with pytest.raises(UndefinedSequenceError) as cm:
            seq[17]
        assert str(cm.value) == "Sequence at position 17 is undefined"
        with pytest.raises(IndexError) as cm:
            seq[30]
        assert str(cm.value) == "sequence index out of range"
        # step = 0
        with pytest.raises(ValueError) as cm:
            s = seq[::0]
        assert str(cm.value) == "slice step cannot be zero"
        # step = 1, stop = +inf
        s = seq[:]  # ?????ABCD?????EFG???
        assert repr(s) == "Seq({5: 'ABCD', 14: 'EFG'}, length=20)"
        s = seq[0:]  # ?????ABCD?????EFG???
        assert repr(s) == "Seq({5: 'ABCD', 14: 'EFG'}, length=20)"
        s = seq[1:]  # ????ABCD?????EFG???
        assert repr(s) == "Seq({4: 'ABCD', 13: 'EFG'}, length=19)"
        s = seq[4:]  # ?ABCD?????EFG???
        assert repr(s) == "Seq({1: 'ABCD', 10: 'EFG'}, length=16)"
        s = seq[5:]  # ABCD?????EFG???
        assert repr(s) == "Seq({0: 'ABCD', 9: 'EFG'}, length=15)"
        s = seq[6:]  # BCD?????EFG???
        assert repr(s) == "Seq({0: 'BCD', 8: 'EFG'}, length=14)"
        s = seq[7:]  # CD?????EFG???
        assert repr(s) == "Seq({0: 'CD', 7: 'EFG'}, length=13)"
        s = seq[8:]  # D?????EFG???
        assert repr(s) == "Seq({0: 'D', 6: 'EFG'}, length=12)"
        s = seq[9:]  # ?????EFG???
        assert repr(s) == "Seq({5: 'EFG'}, length=11)"
        s = seq[10:]  # ????EFG???
        assert repr(s) == "Seq({4: 'EFG'}, length=10)"
        s = seq[13:]  # ?EFG???
        assert repr(s) == "Seq({1: 'EFG'}, length=7)"
        s = seq[14:]  # EFG???
        assert repr(s) == "Seq({0: 'EFG'}, length=6)"
        s = seq[15:]  # FG???
        assert repr(s) == "Seq({0: 'FG'}, length=5)"
        s = seq[16:]  # G???
        assert repr(s) == "Seq({0: 'G'}, length=4)"
        s = seq[17:]  # ???
        assert repr(s) == "Seq(None, length=3)"
        s = seq[18:]  # ??
        assert repr(s) == "Seq(None, length=2)"
        s = seq[19:]  # ?
        assert repr(s) == "Seq(None, length=1)"
        s = seq[20:]  # empty sequence
        assert repr(s) == "Seq('')"
        # step = 1, stop = 9
        s = seq[:9]  # ?????ABCD
        assert repr(s) == "Seq({5: 'ABCD'}, length=9)"
        s = seq[0:9]  # ?????ABCD
        assert repr(s) == "Seq({5: 'ABCD'}, length=9)"
        s = seq[1:9]  # ????ABCD
        assert repr(s) == "Seq({4: 'ABCD'}, length=8)"
        s = seq[4:9]  # ?ABCD
        assert repr(s) == "Seq({1: 'ABCD'}, length=5)"
        s = seq[5:9]  # ABCD
        assert s._data == b"ABCD"
        s = seq[6:9]  # BCD
        assert s._data == b"BCD"
        s = seq[7:9]  # CD
        assert s._data == b"CD"
        s = seq[8:9]  # D
        assert s._data == b"D"
        s = seq[9:9]  # empty sequence
        assert repr(s) == "Seq('')"
        s = seq[10:9]  # empty sequence
        assert repr(s) == "Seq('')"
        # step = 2, stop = +inf
        s = seq[::2]  # ???BD??EG?
        assert repr(s) == "Seq({3: 'BD', 7: 'EG'}, length=10)"
        s = seq[0::2]  # ???BD??EG?
        assert repr(s) == "Seq({3: 'BD', 7: 'EG'}, length=10)"
        s = seq[1::2]  # ??AC???F??
        assert repr(s) == "Seq({2: 'AC', 7: 'F'}, length=10)"
        s = seq[4::2]  # ?BD??EG?
        assert repr(s) == "Seq({1: 'BD', 5: 'EG'}, length=8)"
        s = seq[5::2]  # AC???F??
        assert repr(s) == "Seq({0: 'AC', 5: 'F'}, length=8)"
        s = seq[6::2]  # BD??EG?
        assert repr(s) == "Seq({0: 'BD', 4: 'EG'}, length=7)"
        s = seq[7::2]  # C???F??
        assert repr(s) == "Seq({0: 'C', 4: 'F'}, length=7)"
        s = seq[8::2]  # D??EG?
        assert repr(s) == "Seq({0: 'D', 3: 'EG'}, length=6)"
        s = seq[9::2]  # ???F??
        assert repr(s) == "Seq({3: 'F'}, length=6)"
        s = seq[10::2]  # ??EG?
        assert repr(s) == "Seq({2: 'EG'}, length=5)"
        s = seq[13::2]  # ?F??
        assert repr(s) == "Seq({1: 'F'}, length=4)"
        s = seq[14::2]  # EG?
        assert repr(s) == "Seq({0: 'EG'}, length=3)"
        s = seq[15::2]  # F??
        assert repr(s) == "Seq({0: 'F'}, length=3)"
        s = seq[16::2]  # G?
        assert repr(s) == "Seq({0: 'G'}, length=2)"
        s = seq[17::2]  # ??
        assert repr(s) == "Seq(None, length=2)"
        s = seq[18::2]  # ?
        assert repr(s) == "Seq(None, length=1)"
        s = seq[19::2]  # ?
        assert repr(s) == "Seq(None, length=1)"
        s = seq[20::2]  # empty sequence
        assert repr(s) == "Seq('')"
        # step = -1, start = None
        s = seq[::-1]  # ???GFE?????DCBA?????
        assert repr(s) == "Seq({3: 'GFE', 11: 'DCBA'}, length=20)"
        s = seq[:0:-1]  # ???GFE?????DCBA????
        assert repr(s) == "Seq({3: 'GFE', 11: 'DCBA'}, length=19)"
        s = seq[:1:-1]  # ???GFE?????DCBA???
        assert repr(s) == "Seq({3: 'GFE', 11: 'DCBA'}, length=18)"
        s = seq[:4:-1]  # ???GFE?????DCBA
        assert repr(s) == "Seq({3: 'GFE', 11: 'DCBA'}, length=15)"
        s = seq[:5:-1]  # ???GFE?????DCB
        assert repr(s) == "Seq({3: 'GFE', 11: 'DCB'}, length=14)"
        s = seq[:6:-1]  # ???GFE?????DC
        assert repr(s) == "Seq({3: 'GFE', 11: 'DC'}, length=13)"
        s = seq[:7:-1]  # ???GFE?????D
        assert repr(s) == "Seq({3: 'GFE', 11: 'D'}, length=12)"
        s = seq[:8:-1]  # ???GFE?????
        assert repr(s) == "Seq({3: 'GFE'}, length=11)"
        s = seq[:9:-1]  # ???GFE????
        assert repr(s) == "Seq({3: 'GFE'}, length=10)"
        s = seq[:10:-1]  # ???GFE???
        assert repr(s) == "Seq({3: 'GFE'}, length=9)"
        s = seq[:13:-1]  # ???GFE
        assert repr(s) == "Seq({3: 'GFE'}, length=6)"
        s = seq[:14:-1]  # ???GF
        assert repr(s) == "Seq({3: 'GF'}, length=5)"
        s = seq[:15:-1]  # ???G
        assert repr(s) == "Seq({3: 'G'}, length=4)"
        s = seq[:16:-1]  # ???
        assert repr(s) == "Seq(None, length=3)"
        s = seq[:17:-1]  # ??
        assert repr(s) == "Seq(None, length=2)"
        s = seq[:18:-1]  # ?
        assert repr(s) == "Seq(None, length=1)"
        s = seq[:19:-1]  # empty sequence
        assert repr(s) == "Seq('')"
        s = seq[:20:-1]  # empty sequence
        assert repr(s) == "Seq('')"
        # step = -1, stop = 9
        s = seq[9::-1]  # ?DCBA?????
        assert repr(s) == "Seq({1: 'DCBA'}, length=10)"
        s = seq[9:0:-1]  # ?DCBA????
        assert repr(s) == "Seq({1: 'DCBA'}, length=9)"
        s = seq[9:1:-1]  # ?DCBA???
        assert repr(s) == "Seq({1: 'DCBA'}, length=8)"
        s = seq[9:4:-1]  # ?DCBA
        assert repr(s) == "Seq({1: 'DCBA'}, length=5)"
        s = seq[9:5:-1]  # ?DCB
        assert repr(s) == "Seq({1: 'DCB'}, length=4)"
        s = seq[9:6:-1]  # ?DC
        assert repr(s) == "Seq({1: 'DC'}, length=3)"
        s = seq[9:7:-1]  # ?D
        assert repr(s) == "Seq({1: 'D'}, length=2)"
        s = seq[9:8:-1]  # ?
        assert repr(s) == "Seq(None, length=1)"
        s = seq[9:9:-1]  # empty sequence
        assert repr(s) == "Seq('')"
        s = seq[9:10:-1]  # empty sequence
        assert repr(s) == "Seq('')"
        # step = -2, stop = None
        s = seq[::-2]  # ??F???CA??
        assert repr(s) == "Seq({2: 'F', 6: 'CA'}, length=10)"
        s = seq[:0:-2]  # ??F???CA??
        assert repr(s) == "Seq({2: 'F', 6: 'CA'}, length=10)"
        s = seq[:1:-2]  # ??F???CA?
        assert repr(s) == "Seq({2: 'F', 6: 'CA'}, length=9)"
        s = seq[:4:-2]  # ??F???CA
        assert repr(s) == "Seq({2: 'F', 6: 'CA'}, length=8)"
        s = seq[:5:-2]  # ??F???C
        assert repr(s) == "Seq({2: 'F', 6: 'C'}, length=7)"
        s = seq[:6:-2]  # ??F???C
        assert repr(s) == "Seq({2: 'F', 6: 'C'}, length=7)"
        s = seq[:7:-2]  # ??F???
        assert repr(s) == "Seq({2: 'F'}, length=6)"
        s = seq[:8:-2]  # ??F???
        assert repr(s) == "Seq({2: 'F'}, length=6)"
        s = seq[:9:-2]  # ??F??
        assert repr(s) == "Seq({2: 'F'}, length=5)"
        s = seq[:10:-2]  # ??F??
        assert repr(s) == "Seq({2: 'F'}, length=5)"
        s = seq[:13:-2]  # ??F
        assert repr(s) == "Seq({2: 'F'}, length=3)"
        s = seq[:14:-2]  # ??F
        assert repr(s) == "Seq({2: 'F'}, length=3)"
        s = seq[:15:-2]  # ??
        assert repr(s) == "Seq(None, length=2)"
        s = seq[:16:-2]  # ??
        assert repr(s) == "Seq(None, length=2)"
        s = seq[:17:-2]  # ?
        assert repr(s) == "Seq(None, length=1)"
        s = seq[:18:-2]  # ?
        assert repr(s) == "Seq(None, length=1)"
        s = seq[:19:-2]  # empty sequence
        assert repr(s) == "Seq('')"
        s = seq[:20:-2]  # empty sequence
        assert repr(s) == "Seq('')"
        # test merging segments
        seq = Seq({5: "ABCD", 11: "EFGH"}, length=20)
        s = seq[5::2]
        assert repr(s) == "Seq({0: 'AC', 3: 'EG'}, length=8)"
        s = seq[5::3]
        assert repr(s) == "Seq({0: 'ADEH'}, length=5)"
        s = seq[5::4]
        assert repr(s) == "Seq({0: 'A', 2: 'G'}, length=4)"
        s = seq[5::5]
        assert repr(s) == "Seq({0: 'A'}, length=3)"
        s = seq[4:]
        assert repr(s) == "Seq({1: 'ABCD', 7: 'EFGH'}, length=16)"
        s = seq[4::2]
        assert repr(s) == "Seq({1: 'BD', 4: 'FH'}, length=8)"
        s = seq[4::3]
        assert repr(s) == "Seq({1: 'C', 3: 'G'}, length=6)"
        s = seq[4::4]
        assert repr(s) == "Seq({1: 'DF'}, length=4)"
        # test length 0
        assert Seq({}, length=0) == ""

    def test_addition(self):
        s1 = Seq("ABCD")
        s2 = Seq("EFG")
        u1 = Seq(None, length=7)
        u2 = Seq(None, length=9)
        p1 = Seq({3: "KLM", 11: "XYZ"}, length=17)
        p2 = Seq({0: "PQRST", 8: "HIJ"}, length=13)
        records = SeqIO.parse("TwoBit/sequence.littleendian.2bit", "twobit")
        t = records["seq6"].seq  # ACGTacgtNNNNnn, lazy-loaded
        assert s1 + s1 == Seq("ABCDABCD")
        assert s1 + s2 == Seq("ABCDEFG")
        assert repr(s1 + u1) == "Seq({0: 'ABCD'}, length=11)"
        assert repr(s1 + u2) == "Seq({0: 'ABCD'}, length=13)"
        assert repr(s1 + p1) == "Seq({0: 'ABCD', 7: 'KLM', 15: 'XYZ'}, length=21)"
        assert repr(s1 + p2) == "Seq({0: 'ABCDPQRST', 12: 'HIJ'}, length=17)"
        assert s1 + t == Seq("ABCDACGTacgtNNNNnn")
        assert s2 + s1 == Seq("EFGABCD")
        assert s2 + s2 == Seq("EFGEFG")
        assert repr(s2 + u1) == "Seq({0: 'EFG'}, length=10)"
        assert repr(s2 + u2) == "Seq({0: 'EFG'}, length=12)"
        assert repr(s2 + p1) == "Seq({0: 'EFG', 6: 'KLM', 14: 'XYZ'}, length=20)"
        assert repr(s2 + p2) == "Seq({0: 'EFGPQRST', 11: 'HIJ'}, length=16)"
        assert s2 + t == Seq("EFGACGTacgtNNNNnn")
        assert repr(u1 + s1) == "Seq({7: 'ABCD'}, length=11)"
        assert repr(u1 + s2) == "Seq({7: 'EFG'}, length=10)"
        assert repr(u1 + u1) == "Seq(None, length=14)"
        assert repr(u1 + u2) == "Seq(None, length=16)"
        assert repr(u1 + p1) == "Seq({10: 'KLM', 18: 'XYZ'}, length=24)"
        assert repr(u1 + p2) == "Seq({7: 'PQRST', 15: 'HIJ'}, length=20)"
        assert repr(u1 + t) == "Seq({7: 'ACGTacgtNNNNnn'}, length=21)"
        assert repr(u2 + s1) == "Seq({9: 'ABCD'}, length=13)"
        assert repr(u2 + s2) == "Seq({9: 'EFG'}, length=12)"
        assert repr(u2 + u1) == "Seq(None, length=16)"
        assert repr(u2 + u2) == "Seq(None, length=18)"
        assert repr(u2 + p1) == "Seq({12: 'KLM', 20: 'XYZ'}, length=26)"
        assert repr(u2 + p2) == "Seq({9: 'PQRST', 17: 'HIJ'}, length=22)"
        assert repr(u2 + t) == "Seq({9: 'ACGTacgtNNNNnn'}, length=23)"
        assert repr(p1 + s1) == "Seq({3: 'KLM', 11: 'XYZ', 17: 'ABCD'}, length=21)"
        assert repr(p1 + s2) == "Seq({3: 'KLM', 11: 'XYZ', 17: 'EFG'}, length=20)"
        assert repr(p1 + u1) == "Seq({3: 'KLM', 11: 'XYZ'}, length=24)"
        assert repr(p1 + u2) == "Seq({3: 'KLM', 11: 'XYZ'}, length=26)"
        assert repr(p1 + p1) == "Seq({3: 'KLM', 11: 'XYZ', 20: 'KLM', 28: 'XYZ'}, length=34)"
        assert repr(p1 + p2) == "Seq({3: 'KLM', 11: 'XYZ', 17: 'PQRST', 25: 'HIJ'}, length=30)"
        assert repr(p1 + t) == "Seq({3: 'KLM', 11: 'XYZ', 17: 'ACGTacgtNNNNnn'}, length=31)"
        assert repr(p2 + s1) == "Seq({0: 'PQRST', 8: 'HIJ', 13: 'ABCD'}, length=17)"
        assert repr(p2 + s2) == "Seq({0: 'PQRST', 8: 'HIJ', 13: 'EFG'}, length=16)"
        assert repr(p2 + u1) == "Seq({0: 'PQRST', 8: 'HIJ'}, length=20)"
        assert repr(p2 + u2) == "Seq({0: 'PQRST', 8: 'HIJ'}, length=22)"
        assert repr(p2 + p1) == "Seq({0: 'PQRST', 8: 'HIJ', 16: 'KLM', 24: 'XYZ'}, length=30)"
        assert repr(p2 + p2) == "Seq({0: 'PQRST', 8: 'HIJ', 13: 'PQRST', 21: 'HIJ'}, length=26)"
        assert repr(p2 + t) == "Seq({0: 'PQRST', 8: 'HIJ', 13: 'ACGTacgtNNNNnn'}, length=27)"
        assert t + s1 == Seq("ACGTacgtNNNNnnABCD")
        assert t + s2 == Seq("ACGTacgtNNNNnnEFG")
        assert repr(t + u1) == "Seq({0: 'ACGTacgtNNNNnn'}, length=21)"
        assert repr(t + u2) == "Seq({0: 'ACGTacgtNNNNnn'}, length=23)"
        assert repr(t + p1) == "Seq({0: 'ACGTacgtNNNNnn', 17: 'KLM', 25: 'XYZ'}, length=31)"
        assert repr(t + p2) == "Seq({0: 'ACGTacgtNNNNnnPQRST', 22: 'HIJ'}, length=27)"
        assert t + t == Seq("ACGTacgtNNNNnnACGTacgtNNNNnn")
        p1 = Seq({3: "KLM", 11: "XYZ"}, length=14)
        p2 = Seq({0: "PQRST", 8: "HIJ"}, length=11)
        assert repr(p1 + p2) == "Seq({3: 'KLM', 11: 'XYZPQRST', 22: 'HIJ'}, length=25)"
        assert repr(p2 + p1) == "Seq({0: 'PQRST', 8: 'HIJ', 14: 'KLM', 22: 'XYZ'}, length=25)"
        assert repr(p1 + s1) == "Seq({3: 'KLM', 11: 'XYZABCD'}, length=18)"
        assert repr(p1 + s2) == "Seq({3: 'KLM', 11: 'XYZEFG'}, length=17)"

    def test_multiplication(self):
        p1 = Seq({3: "KLM", 11: "XYZ"}, length=17)
        p2 = Seq({0: "PQRST", 8: "HIJ"}, length=11)
        assert repr(3 * p1) == "Seq({3: 'KLM', 11: 'XYZ', 20: 'KLM', 28: 'XYZ', 37: 'KLM', 45: 'XYZ'}, length=51)"
        assert repr(3 * p2) == "Seq({0: 'PQRST', 8: 'HIJPQRST', 19: 'HIJPQRST', 30: 'HIJ'}, length=33)"

    def test_lower_upper(self):
        u = Seq({3: "KLM", 11: "XYZ"}, length=17)
        l = Seq({0: "pqrst", 8: "hij"}, length=13)  # noqa: E741
        m = Seq({5: "ABCD", 10: "efgh"}, length=20)
        assert repr(u.upper()) == "Seq({3: 'KLM', 11: 'XYZ'}, length=17)"
        assert repr(u.lower()) == "Seq({3: 'klm', 11: 'xyz'}, length=17)"
        assert repr(l.upper()) == "Seq({0: 'PQRST', 8: 'HIJ'}, length=13)"
        assert repr(l.lower()) == "Seq({0: 'pqrst', 8: 'hij'}, length=13)"
        assert repr(m.upper()) == "Seq({5: 'ABCD', 10: 'EFGH'}, length=20)"
        assert repr(m.lower()) == "Seq({5: 'abcd', 10: 'efgh'}, length=20)"

    def test_complement(self):
        s = Seq({3: "AACC", 11: "CGT"}, length=20)
        u = Seq({3: "AACC", 11: "CGU"}, length=20)
        assert repr(s.complement()) == "Seq({3: 'TTGG', 11: 'GCA'}, length=20)"
        assert repr(u.complement()) == "Seq({3: 'TTGG', 11: 'GCA'}, length=20)"
        assert repr(s.reverse_complement()) == "Seq({6: 'ACG', 13: 'GGTT'}, length=20)"
        assert repr(u.reverse_complement()) == "Seq({6: 'ACG', 13: 'GGTT'}, length=20)"
        assert repr(s.complement_rna()) == "Seq({3: 'UUGG', 11: 'GCA'}, length=20)"
        assert repr(u.complement_rna()) == "Seq({3: 'UUGG', 11: 'GCA'}, length=20)"
        assert repr(s.reverse_complement_rna()) == "Seq({6: 'ACG', 13: 'GGUU'}, length=20)"
        assert repr(u.reverse_complement_rna()) == "Seq({6: 'ACG', 13: 'GGUU'}, length=20)"

    def test_replace(self):
        s = Seq({3: "AACC", 11: "CGT"}, length=20)
        assert repr(s.replace("A", "X")) == "Seq({3: 'XXCC', 11: 'CGT'}, length=20)"
        assert repr(s.replace("CC", "YY")) == "Seq({3: 'AAYY', 11: 'CGT'}, length=20)"
        with pytest.raises(UndefinedSequenceError):
            s.replace("A", "XX")

    def test_transcribe(self):
        s = Seq({3: "acgt", 11: "ACGT"}, length=20)
        u = s.transcribe()
        assert repr(u) == "Seq({3: 'acgu', 11: 'ACGU'}, length=20)"
        s = u.back_transcribe()
        assert repr(s) == "Seq({3: 'acgt', 11: 'ACGT'}, length=20)"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
