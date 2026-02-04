# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO BlastIO parsers."""

import os
import unittest
import pytest

from Bio.SearchIO import parse
from Bio.SearchIO.BlastIO.blast_tab import _LONG_SHORT_MAP as all_fields

# test case files are in the Blast directory
TEST_DIR = "Blast"
FMT = "blast-tab"


def get_file(filename):
    """Return the path of a test file."""
    return os.path.join(TEST_DIR, filename)


class BlastTabCases(unittest.TestCase):
    """Tests for the tab-separated BLAST parser."""

    def test_tab_2228_tblastn_001(self):
        """Test parsing TBLASTN 2.2.28+ tabular output (tab_2228_tblastn_001)."""
        tab_file = get_file("tab_2228_tblastn_001.txt")
        qresults = list(
            parse(
                tab_file, FMT, fields=["evalue", "sallseqid", "qseqid"], comments=True
            )
        )

        assert 1 == len(qresults)
        assert 10 == len(qresults[0].hits)
        # there is one hit with an alternative ID
        assert qresults[0]["gi|148227873|ref|NM_001095167.1|"] == qresults[0]["gi|55250552|gb|BC086280.1|"]

        # check some of the HSPs
        assert 0.0 == qresults[0][0][0].evalue
        assert 8e-173 == qresults[0][-1][0].evalue

    def test_tab_2228_tblastx_001(self):
        """Test parsing TBLASTX 2.2.28+ tabular output (tab_2228_tblastx_001)."""
        tab_file = get_file("tab_2228_tblastx_001.txt")
        qresults = list(
            parse(tab_file, FMT, fields=list(all_fields.values()), comments=True)
        )

        # this a single query, with 192 hits and 243 hsps
        assert 1 == len(qresults)
        assert 192 == len(qresults[0].hits)
        assert 243 == sum(len(x) for x in qresults[0])
        # there is one hit with an alternative ID
        assert qresults[0]["gi|31126987|gb|AY255526.2|"] == qresults[0]["gi|31342050|ref|NM_181083.2|"]

        # only checking the new fields in 2.2.28+
        hit = qresults[0][0]
        assert ["NM_001183135", "EF059095"] == hit.accession_all
        assert ["32630", "559292"] == hit.tax_ids
        assert ["N/A", "N/A"] == hit.sci_names
        assert ["N/A", "N/A"] == hit.com_names
        assert ["N/A"] == hit.blast_names
        assert ["N/A"] == hit.super_kingdoms
        assert "Saccharomyces cerevisiae S288c Mon2p (MON2), mRNA" == hit.title
        assert [
                "Saccharomyces cerevisiae S288c Mon2p (MON2), mRNA",
                "Synthetic construct Saccharomyces cerevisiae clone "
                "FLH203015.01X MON2, complete sequence",
            ] == hit.title_all
        assert "N/A" == hit.strand
        assert 100.0 == hit.query_coverage

        for hsp in hit[:4]:
            # shorthand ~ the values just happen to all be 99
            # in other cases, they may be different
            assert 99.0 == hsp.query_coverage
        assert 73.0 == hit[5].query_coverage
        assert 12.0 == hit[6].query_coverage

    def test_tab_2226_tblastn_001(self):
        """Test parsing TBLASTN 2.2.26+ tabular output (tab_2226_tblastn_001)."""
        xml_file = get_file("tab_2226_tblastn_001.txt")
        qresults = parse(xml_file, FMT)
        counter = 0

        # test first qresult
        qresult = next(qresults)
        counter += 1

        assert "gi|16080617|ref|NP_391444.1|" == qresult.id
        assert 3 == len(qresult)

        hit = qresult[0]
        assert "gi|145479850|ref|XM_001425911.1|" == hit.id
        assert "gi|16080617|ref|NP_391444.1|" == hit.query_id
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert "gi|145479850|ref|XM_001425911.1|" == hsp.hit_id
        assert "gi|16080617|ref|NP_391444.1|" == hsp.query_id
        assert 34.88 == hsp.ident_pct
        assert 43 == hsp.aln_span
        assert 28 == hsp.mismatch_num
        assert 0 == hsp.gapopen_num
        assert 30 == hsp.query_start
        assert 73 == hsp.query_end
        assert 1743 == hsp.hit_start
        assert 1872 == hsp.hit_end
        assert 1e-05 == hsp.evalue
        assert 34.7 == hsp.bitscore

        hit = qresult[-1]
        assert "gi|115975252|ref|XM_001180111.1|" == hit.id
        assert "gi|16080617|ref|NP_391444.1|" == hit.query_id
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert "gi|115975252|ref|XM_001180111.1|" == hsp.hit_id
        assert "gi|16080617|ref|NP_391444.1|" == hsp.query_id
        assert 33.90 == hsp.ident_pct
        assert 59 == hsp.aln_span
        assert 31 == hsp.mismatch_num
        assert 1 == hsp.gapopen_num
        assert 43 == hsp.query_start
        assert 94 == hsp.query_end
        assert 1056 == hsp.hit_start
        assert 1233 == hsp.hit_end
        assert 1e-04 == hsp.evalue
        assert 31.6 == hsp.bitscore

        # test last qresult
        qresult = next(qresults)
        counter += 1

        assert "gi|11464971:4-101" == qresult.id
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "gi|350596019|ref|XM_003360601.2|" == hit.id
        assert "gi|11464971:4-101" == hit.query_id
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert "gi|350596019|ref|XM_003360601.2|" == hsp.hit_id
        assert "gi|11464971:4-101" == hsp.query_id
        assert 95.92 == hsp.ident_pct
        assert 98 == hsp.aln_span
        assert 4 == hsp.mismatch_num
        assert 0 == hsp.gapopen_num
        assert 0 == hsp.query_start
        assert 98 == hsp.query_end
        assert 94 == hsp.hit_start
        assert 388 == hsp.hit_end
        assert 2e-67 == hsp.evalue
        assert 199 == hsp.bitscore

        hsp = hit.hsps[-1]
        assert "gi|350596019|ref|XM_003360601.2|" == hsp.hit_id
        assert "gi|11464971:4-101" == hsp.query_id
        assert 29.58 == hsp.ident_pct
        assert 71 == hsp.aln_span
        assert 46 == hsp.mismatch_num
        assert 2 == hsp.gapopen_num
        assert 29 == hsp.query_start
        assert 96 == hsp.query_end
        assert 541 == hsp.hit_start
        assert 754 == hsp.hit_end
        assert 4e-05 == hsp.evalue
        assert 32.7 == hsp.bitscore

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 2 == counter

    def test_tab_2226_tblastn_002(self):
        """Test parsing TBLASTN 2.2.26+ tabular output (tab_2226_tblastn_002)."""
        xml_file = get_file("tab_2226_tblastn_002.txt")
        qresults = parse(xml_file, FMT)

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)

    def test_tab_2226_tblastn_003(self):
        """Test parsing TBLASTN 2.2.26+ tabular output (tab_2226_tblastn_003)."""
        xml_file = get_file("tab_2226_tblastn_003.txt")
        qresults = parse(xml_file, FMT)
        counter = 0

        # test first qresult
        qresult = next(qresults)
        counter += 1

        assert "gi|16080617|ref|NP_391444.1|" == qresult.id
        assert 3 == len(qresult)

        hit = qresult[0]
        assert "gi|145479850|ref|XM_001425911.1|" == hit.id
        assert "gi|16080617|ref|NP_391444.1|" == hit.query_id
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert "gi|145479850|ref|XM_001425911.1|" == hsp.hit_id
        assert "gi|16080617|ref|NP_391444.1|" == hsp.query_id
        assert 34.88 == hsp.ident_pct
        assert 43 == hsp.aln_span
        assert 28 == hsp.mismatch_num
        assert 0 == hsp.gapopen_num
        assert 30 == hsp.query_start
        assert 73 == hsp.query_end
        assert 1743 == hsp.hit_start
        assert 1872 == hsp.hit_end
        assert 1e-05 == hsp.evalue
        assert 34.7 == hsp.bitscore

        hit = qresult[-1]
        assert "gi|115975252|ref|XM_001180111.1|" == hit.id
        assert "gi|16080617|ref|NP_391444.1|" == hit.query_id
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert "gi|115975252|ref|XM_001180111.1|" == hsp.hit_id
        assert "gi|16080617|ref|NP_391444.1|" == hsp.query_id
        assert 33.90 == hsp.ident_pct
        assert 59 == hsp.aln_span
        assert 31 == hsp.mismatch_num
        assert 1 == hsp.gapopen_num
        assert 43 == hsp.query_start
        assert 94 == hsp.query_end
        assert 1056 == hsp.hit_start
        assert 1233 == hsp.hit_end
        assert 1e-04 == hsp.evalue
        assert 31.6 == hsp.bitscore

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_tab_2226_tblastn_004(self):
        """Test parsing TBLASTN 2.2.26+ tabular output (tab_2226_tblastn_004)."""
        xml_file = get_file("tab_2226_tblastn_004.txt")
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1

        assert "gi|11464971:4-101" == qresult.id
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "gi|350596019|ref|XM_003360601.2|" == hit.id
        assert "gi|11464971:4-101" == hit.query_id
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert "gi|350596019|ref|XM_003360601.2|" == hsp.hit_id
        assert "gi|11464971:4-101" == hsp.query_id
        assert 95.92 == hsp.ident_pct
        assert 98 == hsp.aln_span
        assert 4 == hsp.mismatch_num
        assert 0 == hsp.gapopen_num
        assert 0 == hsp.query_start
        assert 98 == hsp.query_end
        assert 94 == hsp.hit_start
        assert 388 == hsp.hit_end
        assert 2e-67 == hsp.evalue
        assert 199 == hsp.bitscore

        hsp = hit.hsps[-1]
        assert "gi|350596019|ref|XM_003360601.2|" == hsp.hit_id
        assert "gi|11464971:4-101" == hsp.query_id
        assert 29.58 == hsp.ident_pct
        assert 71 == hsp.aln_span
        assert 46 == hsp.mismatch_num
        assert 2 == hsp.gapopen_num
        assert 29 == hsp.query_start
        assert 96 == hsp.query_end
        assert 541 == hsp.hit_start
        assert 754 == hsp.hit_end
        assert 4e-05 == hsp.evalue
        assert 32.7 == hsp.bitscore

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_tab_2226_tblastn_005(self):
        """Test parsing TBLASTN 2.2.26+ tabular output with comments (tab_2226_tblastn_005)."""
        xml_file = get_file("tab_2226_tblastn_005.txt")
        qresults = parse(xml_file, FMT, comments=True)
        counter = 0

        # test first qresult
        qresult = next(qresults)
        counter += 1

        assert "tblastn" == qresult.program
        assert "db/minirefseq_mrna" == qresult.target
        assert "2.2.26+" == qresult.version
        assert "random_s00" == qresult.id
        assert 0 == len(qresult)

        # test second qresult
        qresult = next(qresults)
        counter += 1

        assert "tblastn" == qresult.program
        assert "db/minirefseq_mrna" == qresult.target
        assert "2.2.26+" == qresult.version
        assert "gi|16080617|ref|NP_391444.1|" == qresult.id
        assert 3 == len(qresult)

        hit = qresult[0]
        assert "gi|145479850|ref|XM_001425911.1|" == hit.id
        assert "gi|16080617|ref|NP_391444.1|" == hit.query_id
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert "gi|145479850|ref|XM_001425911.1|" == hsp.hit_id
        assert "gi|16080617|ref|NP_391444.1|" == hsp.query_id
        assert 34.88 == hsp.ident_pct
        assert 43 == hsp.aln_span
        assert 28 == hsp.mismatch_num
        assert 0 == hsp.gapopen_num
        assert 30 == hsp.query_start
        assert 73 == hsp.query_end
        assert 1743 == hsp.hit_start
        assert 1872 == hsp.hit_end
        assert 1e-05 == hsp.evalue
        assert 34.7 == hsp.bitscore

        hit = qresult[-1]
        assert "gi|115975252|ref|XM_001180111.1|" == hit.id
        assert "gi|16080617|ref|NP_391444.1|" == hit.query_id
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert "gi|115975252|ref|XM_001180111.1|" == hsp.hit_id
        assert "gi|16080617|ref|NP_391444.1|" == hsp.query_id
        assert 33.90 == hsp.ident_pct
        assert 59 == hsp.aln_span
        assert 31 == hsp.mismatch_num
        assert 1 == hsp.gapopen_num
        assert 43 == hsp.query_start
        assert 94 == hsp.query_end
        assert 1056 == hsp.hit_start
        assert 1233 == hsp.hit_end
        assert 1e-04 == hsp.evalue
        assert 31.6 == hsp.bitscore

        # test last qresult
        qresult = next(qresults)
        counter += 1

        assert "tblastn" == qresult.program
        assert "db/minirefseq_mrna" == qresult.target
        assert "2.2.26+" == qresult.version
        assert "gi|11464971:4-101" == qresult.id
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "gi|350596019|ref|XM_003360601.2|" == hit.id
        assert "gi|11464971:4-101" == hit.query_id
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert "gi|350596019|ref|XM_003360601.2|" == hsp.hit_id
        assert "gi|11464971:4-101" == hsp.query_id
        assert 95.92 == hsp.ident_pct
        assert 98 == hsp.aln_span
        assert 4 == hsp.mismatch_num
        assert 0 == hsp.gapopen_num
        assert 0 == hsp.query_start
        assert 98 == hsp.query_end
        assert 94 == hsp.hit_start
        assert 388 == hsp.hit_end
        assert 2e-67 == hsp.evalue
        assert 199 == hsp.bitscore

        hsp = hit.hsps[-1]
        assert "gi|350596019|ref|XM_003360601.2|" == hsp.hit_id
        assert "gi|11464971:4-101" == hsp.query_id
        assert 29.58 == hsp.ident_pct
        assert 71 == hsp.aln_span
        assert 46 == hsp.mismatch_num
        assert 2 == hsp.gapopen_num
        assert 29 == hsp.query_start
        assert 96 == hsp.query_end
        assert 541 == hsp.hit_start
        assert 754 == hsp.hit_end
        assert 4e-05 == hsp.evalue
        assert 32.7 == hsp.bitscore

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 3 == counter

    def test_tab_2226_tblastn_005_comments_false(self):
        """Test parsing TBLASTN 2.2.26+ tabular output with comments (tab_2226_tblastn_005)."""
        tab_file = get_file("tab_2226_tblastn_005.txt")
        exc_msg = (
            "Encountered unexpected character '#' at the beginning of a line. "
            "Set comments=True if the file is a commented file."
        )
        qresults = parse(tab_file, FMT)
        with pytest.raises(ValueError):
            next(qresults)

    def test_tab_2226_tblastn_006(self):
        """Test parsing TBLASTN 2.2.26+ tabular output with comments (tab_2226_tblastn_006)."""
        xml_file = get_file("tab_2226_tblastn_006.txt")
        qresults = parse(xml_file, FMT, comments=True)
        counter = 0

        qresult = next(qresults)
        counter += 1

        assert "tblastn" == qresult.program
        assert "db/minirefseq_mrna" == qresult.target
        assert "2.2.26+" == qresult.version
        assert "random_s00" == qresult.id
        assert 0 == len(qresult)

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_tab_2226_tblastn_007(self):
        """Test parsing TBLASTN 2.2.26+ tabular output with comments (tab_2226_tblastn_007)."""
        xml_file = get_file("tab_2226_tblastn_007.txt")
        qresults = parse(xml_file, FMT, comments=True)
        counter = 0

        qresult = next(qresults)
        counter += 1

        assert "tblastn" == qresult.program
        assert "db/minirefseq_mrna" == qresult.target
        assert "2.2.26+" == qresult.version
        assert "gi|16080617|ref|NP_391444.1|" == qresult.id
        assert 3 == len(qresult)

        hit = qresult[0]
        assert "gi|145479850|ref|XM_001425911.1|" == hit.id
        assert "gi|16080617|ref|NP_391444.1|" == hit.query_id
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert "gi|145479850|ref|XM_001425911.1|" == hsp.hit_id
        assert "gi|16080617|ref|NP_391444.1|" == hsp.query_id
        assert 34.88 == hsp.ident_pct
        assert 43 == hsp.aln_span
        assert 28 == hsp.mismatch_num
        assert 0 == hsp.gapopen_num
        assert 30 == hsp.query_start
        assert 73 == hsp.query_end
        assert 1743 == hsp.hit_start
        assert 1872 == hsp.hit_end
        assert 1e-05 == hsp.evalue
        assert 34.7 == hsp.bitscore

        hit = qresult[-1]
        assert "gi|115975252|ref|XM_001180111.1|" == hit.id
        assert "gi|16080617|ref|NP_391444.1|" == hit.query_id
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert "gi|115975252|ref|XM_001180111.1|" == hsp.hit_id
        assert "gi|16080617|ref|NP_391444.1|" == hsp.query_id
        assert 33.90 == hsp.ident_pct
        assert 59 == hsp.aln_span
        assert 31 == hsp.mismatch_num
        assert 1 == hsp.gapopen_num
        assert 43 == hsp.query_start
        assert 94 == hsp.query_end
        assert 1056 == hsp.hit_start
        assert 1233 == hsp.hit_end
        assert 1e-04 == hsp.evalue
        assert 31.6 == hsp.bitscore

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_tab_2226_tblastn_008(self):
        """Test parsing TBLASTN 2.2.26+ tabular output with comments (tab_2226_tblastn_008)."""
        xml_file = get_file("tab_2226_tblastn_008.txt")
        qresults = parse(xml_file, FMT, comments=True)
        counter = 0

        qresult = next(qresults)
        counter += 1

        assert "tblastn" == qresult.program
        assert "db/minirefseq_mrna" == qresult.target
        assert "2.2.26+" == qresult.version
        assert "gi|11464971:4-101" == qresult.id
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "gi|350596019|ref|XM_003360601.2|" == hit.id
        assert "gi|11464971:4-101" == hit.query_id
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert "gi|350596019|ref|XM_003360601.2|" == hsp.hit_id
        assert "gi|11464971:4-101" == hsp.query_id
        assert 95.92 == hsp.ident_pct
        assert 98 == hsp.aln_span
        assert 4 == hsp.mismatch_num
        assert 0 == hsp.gapopen_num
        assert 0 == hsp.query_start
        assert 98 == hsp.query_end
        assert 94 == hsp.hit_start
        assert 388 == hsp.hit_end
        assert 2e-67 == hsp.evalue
        assert 199 == hsp.bitscore

        hsp = hit.hsps[-1]
        assert "gi|350596019|ref|XM_003360601.2|" == hsp.hit_id
        assert "gi|11464971:4-101" == hsp.query_id
        assert 29.58 == hsp.ident_pct
        assert 71 == hsp.aln_span
        assert 46 == hsp.mismatch_num
        assert 2 == hsp.gapopen_num
        assert 29 == hsp.query_start
        assert 96 == hsp.query_end
        assert 541 == hsp.hit_start
        assert 754 == hsp.hit_end
        assert 4e-05 == hsp.evalue
        assert 32.7 == hsp.bitscore

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter

    def test_tab_2226_tblastn_009(self):
        """Test parsing TBLASTN 2.2.26+ tabular output (tab_2226_tblastn_009)."""
        xml_file = get_file("tab_2226_tblastn_009.txt")
        qresults = parse(xml_file, FMT, fields=("qseqid", "sseqid"))
        counter = 0

        # test first qresult
        qresult = next(qresults)
        counter += 1

        assert "<unknown program>" == qresult.program
        assert "<unknown target>" == qresult.target
        assert "<unknown version>" == qresult.version
        assert "gi|16080617|ref|NP_391444.1|" == qresult.id
        assert 3 == len(qresult)

        hit = qresult[0]
        assert "gi|145479850|ref|XM_001425911.1|" == hit.id
        assert "gi|16080617|ref|NP_391444.1|" == hit.query_id
        assert 1 == len(hit)

        hsp = hit[0]
        assert "gi|145479850|ref|XM_001425911.1|" == hsp.hit_id
        assert "gi|16080617|ref|NP_391444.1|" == hsp.query_id

        hit = qresult[-1]
        assert "gi|115975252|ref|XM_001180111.1|" == hit.id
        assert "gi|16080617|ref|NP_391444.1|" == hit.query_id
        assert 1 == len(hit)

        hsp = hit[0]
        assert "gi|115975252|ref|XM_001180111.1|" == hsp.hit_id
        assert "gi|16080617|ref|NP_391444.1|" == hsp.query_id

        # test last qresult
        qresult = next(qresults)
        counter += 1

        assert "<unknown program>" == qresult.program
        assert "<unknown target>" == qresult.target
        assert "<unknown version>" == qresult.version
        assert "gi|11464971:4-101" == qresult.id
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "gi|350596019|ref|XM_003360601.2|" == hit.id
        assert "gi|11464971:4-101" == hit.query_id
        assert 2 == len(hit)

        hsp = hit[0]
        assert "gi|350596019|ref|XM_003360601.2|" == hsp.hit_id
        assert "gi|11464971:4-101" == hsp.query_id

        hsp = hit[-1]
        assert "gi|350596019|ref|XM_003360601.2|" == hsp.hit_id
        assert "gi|11464971:4-101" == hsp.query_id

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 2 == counter

    def test_tab_2226_tblastn_010(self):
        """Test parsing TBLASTN 2.2.26+ tabular output with comments (tab_2226_tblastn_010)."""
        xml_file = get_file("tab_2226_tblastn_010.txt")
        qresults = parse(xml_file, FMT, comments=True)
        counter = 0

        # test first qresult
        qresult = next(qresults)
        counter += 1

        assert "tblastn" == qresult.program
        assert "db/minirefseq_mrna" == qresult.target
        assert "2.2.26+" == qresult.version
        assert "random_s00" == qresult.id
        assert 0 == len(qresult)

        # test second qresult
        qresult = next(qresults)
        counter += 1

        assert "tblastn" == qresult.program
        assert "db/minirefseq_mrna" == qresult.target
        assert "2.2.26+" == qresult.version
        assert "gi|16080617|ref|NP_391444.1|" == qresult.id
        assert 3 == len(qresult)

        hit = qresult[0]
        assert "gi|145479850|ref|XM_001425911.1|" == hit.id
        assert "gi|16080617|ref|NP_391444.1|" == hit.query_id
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert "gi|145479850|ref|XM_001425911.1|" == hsp.hit_id
        assert "gi|16080617|ref|NP_391444.1|" == hsp.query_id
        assert 1e-05 == hsp.evalue
        assert 34.7 == hsp.bitscore

        hit = qresult[-1]
        assert "gi|115975252|ref|XM_001180111.1|" == hit.id
        assert "gi|16080617|ref|NP_391444.1|" == hit.query_id
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert "gi|115975252|ref|XM_001180111.1|" == hsp.hit_id
        assert "gi|16080617|ref|NP_391444.1|" == hsp.query_id
        assert 1e-04 == hsp.evalue
        assert 31.6 == hsp.bitscore

        # test last qresult
        qresult = next(qresults)
        counter += 1

        assert "tblastn" == qresult.program
        assert "db/minirefseq_mrna" == qresult.target
        assert "2.2.26+" == qresult.version
        assert "gi|11464971:4-101" == qresult.id
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "gi|350596019|ref|XM_003360601.2|" == hit.id
        assert "gi|11464971:4-101" == hit.query_id
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert "gi|350596019|ref|XM_003360601.2|" == hsp.hit_id
        assert "gi|11464971:4-101" == hsp.query_id
        assert 2e-67 == hsp.evalue
        assert 199 == hsp.bitscore

        hsp = hit.hsps[-1]
        assert "gi|350596019|ref|XM_003360601.2|" == hsp.hit_id
        assert "gi|11464971:4-101" == hsp.query_id
        assert 4e-05 == hsp.evalue
        assert 32.7 == hsp.bitscore

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 3 == counter

    def test_tab_2226_tblastn_011(self):
        """Test parsing TBLASTN 2.2.26+ tabular output with comments (tab_2226_tblastn_011)."""
        xml_file = get_file("tab_2226_tblastn_011.txt")
        qresults = parse(xml_file, FMT, comments=True)
        counter = 0

        # test first qresult
        qresult = next(qresults)
        counter += 1

        assert "tblastn" == qresult.program
        assert "db/minirefseq_mrna" == qresult.target
        assert "2.2.26+" == qresult.version
        assert "random_s00" == qresult.id
        assert 0 == len(qresult)

        # test second qresult
        qresult = next(qresults)
        counter += 1

        assert "tblastn" == qresult.program
        assert "db/minirefseq_mrna" == qresult.target
        assert "2.2.26+" == qresult.version
        assert "gi|16080617|ref|NP_391444.1|" == qresult.id
        assert "gi|16080617|ref|NP_391444.1|" == qresult.accession
        assert "gi|16080617|ref|NP_391444.1|" == qresult.accession_version
        assert "0" == qresult.gi
        assert 102 == qresult.seq_len
        assert 3 == len(qresult)

        hit = qresult[0]
        assert "gi|145479850|ref|XM_001425911.1|" == hit.id
        assert ["gi|145479850|ref|XM_001425911.1|"] == hit.id_all
        assert "gi|145479850|ref|XM_001425911.1|" == hit.accession
        assert "gi|145479850|ref|XM_001425911.1|" == hit.accession_version
        assert "0" == hit.gi
        assert "0" == hit.gi_all
        assert "gi|16080617|ref|NP_391444.1|" == hit.query_id
        assert 4632 == hit.seq_len
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert "gi|145479850|ref|XM_001425911.1|" == hsp.hit_id
        assert "gi|16080617|ref|NP_391444.1|" == hsp.query_id
        assert 34.88 == hsp.ident_pct
        assert 43 == hsp.aln_span
        assert 28 == hsp.mismatch_num
        assert 0 == hsp.gapopen_num
        assert 30 == hsp.query_start
        assert 73 == hsp.query_end
        assert 1743 == hsp.hit_start
        assert 1872 == hsp.hit_end
        assert 1e-05 == hsp.evalue
        assert 34.7 == hsp.bitscore
        assert "PDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLD" == hsp.query.seq
        assert "PKTATGTKKGTIIGLLSIHTILFILTSHALSLEVKEQT*KDID" == hsp.hit.seq
        assert 78 == hsp.bitscore_raw
        assert 15 == hsp.ident_num
        assert 26 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 60.47 == hsp.pos_pct
        assert 0 == hsp.query_frame
        assert 1 == hsp.hit_frame

        hit = qresult[-1]
        assert "gi|115975252|ref|XM_001180111.1|" == hit.id
        assert "gi|115975252|ref|XM_001180111.1|" == hit.accession
        assert "gi|115975252|ref|XM_001180111.1|" == hit.accession_version
        assert "gi|16080617|ref|NP_391444.1|" == hit.query_id
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert "gi|115975252|ref|XM_001180111.1|" == hsp.hit_id
        assert "gi|16080617|ref|NP_391444.1|" == hsp.query_id
        assert 33.90 == hsp.ident_pct
        assert 59 == hsp.aln_span
        assert 31 == hsp.mismatch_num
        assert 1 == hsp.gapopen_num
        assert 43 == hsp.query_start
        assert 94 == hsp.query_end
        assert 1056 == hsp.hit_start
        assert 1233 == hsp.hit_end
        assert 1e-04 == hsp.evalue
        assert 31.6 == hsp.bitscore
        assert "GLADTHTIEVTVDNEPVSLDITEESTSDLDKFNSG--------DKVTITYEKNDEGQLL" == hsp.query.seq
        assert "GLVPDHTLILPVGHYQSMLDLTEEVQTELDQFKSALRKYYLSKGKTCVIYERNFRTQHL" == hsp.hit.seq
        assert 70.0 == hsp.bitscore_raw
        assert 20 == hsp.ident_num
        assert 29 == hsp.pos_num
        assert 8 == hsp.gap_num
        assert 49.15 == hsp.pos_pct
        assert 0 == hsp.query_frame
        assert 1 == hsp.hit_frame

        # test last qresult
        qresult = next(qresults)
        counter += 1

        assert "tblastn" == qresult.program
        assert "db/minirefseq_mrna" == qresult.target
        assert "2.2.26+" == qresult.version
        assert "gi|11464971:4-101" == qresult.id
        assert "gi|11464971:4-101" == qresult.accession
        assert "gi|11464971:4-101" == qresult.accession_version
        assert "0" == qresult.gi
        assert 98 == qresult.seq_len
        assert 5 == len(qresult)

        hit = qresult[0]
        assert "gi|350596019|ref|XM_003360601.2|" == hit.id
        assert ["gi|350596019|ref|XM_003360601.2|"] == hit.id_all
        assert "gi|350596019|ref|XM_003360601.2|" == hit.accession
        assert "gi|350596019|ref|XM_003360601.2|" == hit.accession_version
        assert "0" == hit.gi
        assert "0" == hit.gi_all
        assert "gi|11464971:4-101" == hit.query_id
        assert 772 == hit.seq_len
        assert 2 == len(hit)

        hsp = hit.hsps[0]
        assert "gi|350596019|ref|XM_003360601.2|" == hsp.hit_id
        assert "gi|11464971:4-101" == hsp.query_id
        assert 95.92 == hsp.ident_pct
        assert 98 == hsp.aln_span
        assert 4 == hsp.mismatch_num
        assert 0 == hsp.gapopen_num
        assert 0 == hsp.query_start
        assert 98 == hsp.query_end
        assert 94 == hsp.hit_start
        assert 388 == hsp.hit_end
        assert 2e-67 == hsp.evalue
        assert 199 == hsp.bitscore
        assert "KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK" == hsp.query.seq
        assert "KRIREGYLVKKGSMFNTWKPMWVILLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDGWVRDIKKAIK" == hsp.hit.seq
        assert 506.0 == hsp.bitscore_raw
        assert 94 == hsp.ident_num
        assert 96 == hsp.pos_num
        assert 0 == hsp.gap_num
        assert 97.96 == hsp.pos_pct
        assert 0 == hsp.query_frame
        assert 2 == hsp.hit_frame

        hsp = hit.hsps[-1]
        assert "gi|350596019|ref|XM_003360601.2|" == hsp.hit_id
        assert "gi|11464971:4-101" == hsp.query_id
        assert 29.58 == hsp.ident_pct
        assert 71 == hsp.aln_span
        assert 46 == hsp.mismatch_num
        assert 2 == hsp.gapopen_num
        assert 29 == hsp.query_start
        assert 96 == hsp.query_end
        assert 541 == hsp.hit_start
        assert 754 == hsp.hit_end
        assert 4e-05 == hsp.evalue
        assert 32.7 == hsp.bitscore
        assert "IEFYKKKSDNSPKGMIPLKGSTLTS-PCQDFGKRMFVLK---ITTTKQQDHFFQAAFLEERDAWVRDIKKA" == hsp.query.seq
        assert "LHYYDPAGGEDPLGAIHLRGCVVTSVESNTDGKNGFLWERAXXITADEVHYFLQAANPKERTEWIKAIQVA" == hsp.hit.seq
        assert 73.0 == hsp.bitscore_raw
        assert 21 == hsp.ident_num
        assert 33 == hsp.pos_num
        assert 4 == hsp.gap_num
        assert 46.48 == hsp.pos_pct
        assert 0 == hsp.query_frame
        assert 2 == hsp.hit_frame

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 3 == counter

    def test_tab_2226_tblastn_012(self):
        """Test parsing TBLASTN 2.2.26+ tabular output with comments (tab_2226_tblastn_012)."""
        xml_file = get_file("tab_2226_tblastn_012.txt")
        qresults = parse(xml_file, FMT, comments=True)
        counter = 0

        # test first qresult
        qresult = next(qresults)
        counter += 1

        assert "tblastn" == qresult.program
        assert "refseq_rna" == qresult.target
        assert "2.2.26+" == qresult.version
        assert "random_s00" == qresult.id
        assert "X76FDCG9016" == qresult.rid
        assert 0 == len(qresult)

        # test second qresult
        qresult = next(qresults)
        counter += 1

        assert "tblastn" == qresult.program
        assert "refseq_rna" == qresult.target
        assert "2.2.26+" == qresult.version
        assert "gi|16080617|ref|NP_391444.1|" == qresult.id
        assert "X76FDCG9016" == qresult.rid
        assert 3 == len(qresult)

        # test last qresult
        qresult = next(qresults)
        counter += 1

        assert "tblastn" == qresult.program
        assert "refseq_rna" == qresult.target
        assert "2.2.26+" == qresult.version
        assert "gi|11464971:4-101" == qresult.id
        assert "X76FDCG9016" == qresult.rid
        assert 5 == len(qresult)

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 3 == counter

    def test_tab_2226_tblastn_013(self):
        """Test parsing TBLASTN 2.2.26+ tabular output (tab_2226_tblastn_013)."""
        xml_file = get_file("tab_2226_tblastn_013.txt")
        qresults = parse(xml_file, FMT, fields="qseq std sseq")
        counter = 0

        qresult = next(qresults)
        counter += 1

        assert "<unknown program>" == qresult.program
        assert "<unknown target>" == qresult.target
        assert "<unknown version>" == qresult.version
        assert 3 == len(qresult)

        hit = qresult[0]
        assert "gi|145479850|ref|XM_001425911.1|" == hit.id
        assert "gi|16080617|ref|NP_391444.1|" == hit.query_id
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert "gi|145479850|ref|XM_001425911.1|" == hsp.hit_id
        assert "gi|16080617|ref|NP_391444.1|" == hsp.query_id
        assert 34.88 == hsp.ident_pct
        assert 43 == hsp.aln_span
        assert 28 == hsp.mismatch_num
        assert 0 == hsp.gapopen_num
        assert 30 == hsp.query_start
        assert 73 == hsp.query_end
        assert 1743 == hsp.hit_start
        assert 1872 == hsp.hit_end
        assert 1e-05 == hsp.evalue
        assert 34.7 == hsp.bitscore
        assert "PDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLD" == hsp.query.seq
        assert "PKTATGTKKGTIIGLLSIHTILFILTSHALSLEVKEQT*KDID" == hsp.hit.seq

        hit = qresult[-1]
        assert "gi|115975252|ref|XM_001180111.1|" == hit.id
        assert "gi|16080617|ref|NP_391444.1|" == hit.query_id
        assert 1 == len(hit)

        hsp = hit.hsps[0]
        assert "gi|115975252|ref|XM_001180111.1|" == hsp.hit_id
        assert "gi|16080617|ref|NP_391444.1|" == hsp.query_id
        assert 33.90 == hsp.ident_pct
        assert 59 == hsp.aln_span
        assert 31 == hsp.mismatch_num
        assert 1 == hsp.gapopen_num
        assert 43 == hsp.query_start
        assert 94 == hsp.query_end
        assert 1056 == hsp.hit_start
        assert 1233 == hsp.hit_end
        assert 1e-04 == hsp.evalue
        assert 31.6 == hsp.bitscore
        assert "GLADTHTIEVTVDNEPVSLDITEESTSDLDKFNSG--------DKVTITYEKNDEGQLL" == hsp.query.seq
        assert "GLVPDHTLILPVGHYQSMLDLTEEVQTELDQFKSALRKYYLSKGKTCVIYERNFRTQHL" == hsp.hit.seq

        # check if we've finished iteration over qresults
        with pytest.raises(StopIteration):
            next(qresults)
        assert 1 == counter


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
