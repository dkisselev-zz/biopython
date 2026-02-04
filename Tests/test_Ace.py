# Copyright 2004 by Frank Kauff.  All rights reserved.
# Revisions copyright 2008-2013 by Peter Cock. All rights reserved.
# Revisions copyright 2009-2009 by Michiel de Hoon. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for Ace module."""

import pytest

from Bio.Sequencing import Ace


class TestAceTestOne:
    @pytest.fixture(autouse=True)
    def _setup(self):
        self.handle = open("Ace/contig1.ace")

        yield
        self.handle.close()

    def test_check_ACEParser(self):
        """Test to check that ACEParser can parse the whole file into one record."""
        record = Ace.read(self.handle)
        assert record.ncontigs == 2
        assert record.nreads == 16
        assert len(record.wa) == 1
        assert record.wa[0].tag_type == "phrap_params"
        assert record.wa[0].program == "phrap"
        assert record.wa[0].date == "040203:114710"
        assert record.wa[0].info == [
                "phrap 304_nuclsu.fasta.screen -new_ace -retain_duplicates",
                "phrap version 0.990329",
            ]
        assert len(record.contigs) == 2
        assert len(record.contigs[0].reads) == 2
        assert record.contigs[0].name == "Contig1"
        assert record.contigs[0].nbases == 856
        assert record.contigs[0].nreads == 2
        assert record.contigs[0].nsegments == 31
        assert record.contigs[0].uorc == "U"
        center = len(record.contigs[0].sequence) // 2
        assert record.contigs[0].sequence[:10] == "aatacgGGAT"
        assert record.contigs[0].sequence[center - 5 : center + 5] == "ACATCATCTG"
        assert record.contigs[0].sequence[-10:] == "cATCTAGtac"
        center = len(record.contigs[0].quality) // 2
        assert record.contigs[0].quality[:10] == [0, 0, 0, 0, 0, 0, 22, 23, 25, 28]
        assert record.contigs[0].quality[center - 5 : center + 5] == [90, 90, 90, 90, 90, 90, 90, 90, 90, 90]
        assert record.contigs[0].quality[-10:] == [15, 22, 30, 24, 28, 22, 21, 15, 19, 0]
        assert len(record.contigs[0].af) == 2
        assert record.contigs[0].af[1].name == "BL060c3-LR0R.b.ab1"
        assert record.contigs[0].af[1].coru == "U"
        assert record.contigs[0].af[1].padded_start == 1
        assert len(record.contigs[0].bs) == 31
        assert record.contigs[0].bs[15].name == "BL060c3-LR5.g.ab1"
        assert record.contigs[0].bs[15].padded_start == 434
        assert record.contigs[0].bs[15].padded_end == 438
        assert record.contigs[0].bs[30].name == "BL060c3-LR0R.b.ab1"
        assert record.contigs[0].bs[30].padded_start == 823
        assert record.contigs[0].bs[30].padded_end == 856
        assert len(record.contigs[0].ct) == 1
        assert record.contigs[0].ct[0].name == "Contig1"
        assert record.contigs[0].ct[0].tag_type == "repeat"
        assert record.contigs[0].ct[0].program == "phrap"
        assert record.contigs[0].ct[0].padded_start == 52
        assert record.contigs[0].ct[0].padded_end == 53
        assert record.contigs[0].ct[0].date == "555456:555432"
        assert record.contigs[0].ct[0].info == ["This is the first line of comment for c1", "and this the second for c1"]
        assert record.contigs[0].wa is None
        assert len(record.contigs[0].reads) == 2
        assert record.contigs[0].reads[0].rd.name == "BL060c3-LR5.g.ab1"
        assert record.contigs[0].reads[0].rd.padded_bases == 868
        assert record.contigs[0].reads[0].rd.info_items == 0
        assert record.contigs[0].reads[0].rd.read_tags == 0
        center = len(record.contigs[0].reads[0].rd.sequence) // 2
        assert record.contigs[0].reads[0].rd.sequence[:10] == "tagcgaggaa"
        assert record.contigs[0].reads[0].rd.sequence[center - 5 : center + 5] == "CCGAGGCCAA"
        assert record.contigs[0].reads[0].rd.sequence[-10:] == "gaaccatcag"
        assert record.contigs[0].reads[0].qa.qual_clipping_start == 80
        assert record.contigs[0].reads[0].qa.qual_clipping_end == 853
        assert record.contigs[0].reads[0].qa.align_clipping_start == 22
        assert record.contigs[0].reads[0].qa.align_clipping_end == 856
        assert record.contigs[0].reads[0].ds is None
        assert len(record.contigs[0].reads[0].rt) == 4
        assert record.contigs[0].reads[0].rt[0].name == "BL060c3-LR5.g.ab1"
        assert record.contigs[0].reads[0].rt[0].tag_type == "matchElsewhereHighQual"
        assert record.contigs[0].reads[0].rt[0].program == "phrap"
        assert record.contigs[0].reads[0].rt[0].padded_start == 590
        assert record.contigs[0].reads[0].rt[0].padded_end == 607
        assert record.contigs[0].reads[0].rt[0].date == "040217:110357"
        assert record.contigs[0].reads[0].rt[1].name == "BL060c3-LR5.g.ab1"
        assert record.contigs[0].reads[0].rt[1].tag_type == "matchElsewhereHighQual"
        assert record.contigs[0].reads[0].rt[1].program == "phrap"
        assert record.contigs[0].reads[0].rt[1].padded_start == 617
        assert record.contigs[0].reads[0].rt[1].padded_end == 631
        assert record.contigs[0].reads[0].rt[1].date == "040217:110357"
        assert record.contigs[0].reads[0].rt[2].name == "BL060c3-LR5.g.ab1"
        assert record.contigs[0].reads[0].rt[2].tag_type == "matchElsewhereHighQual"
        assert record.contigs[0].reads[0].rt[2].program == "phrap"
        assert record.contigs[0].reads[0].rt[2].padded_start == 617
        assert record.contigs[0].reads[0].rt[2].padded_end == 631
        assert record.contigs[0].reads[0].rt[2].date == "040217:110357"
        assert record.contigs[0].reads[0].rt[3].name == "BL060c3-LR5.g.ab1"
        assert record.contigs[0].reads[0].rt[3].tag_type == "matchElsewhereHighQual"
        assert record.contigs[0].reads[0].rt[3].program == "phrap"
        assert record.contigs[0].reads[0].rt[3].padded_start == 617
        assert record.contigs[0].reads[0].rt[3].padded_end == 631
        assert record.contigs[0].reads[0].rt[3].date == "040217:110357"

        assert len(record.contigs[0].reads[0].wr) == 1
        assert record.contigs[0].reads[0].wr[0].name == "BL060c3-LR5.g.ab1"
        assert record.contigs[0].reads[0].wr[0].aligned == "unaligned"
        assert record.contigs[0].reads[0].wr[0].program == "phrap"
        assert record.contigs[0].reads[0].wr[0].date == "040217:110357"
        assert record.contigs[0].reads[1].rd.name == "BL060c3-LR0R.b.ab1"
        assert record.contigs[0].reads[1].rd.padded_bases == 856
        assert record.contigs[0].reads[1].rd.info_items == 0
        assert record.contigs[0].reads[1].rd.read_tags == 0
        center = len(record.contigs[0].reads[1].rd.sequence) // 2
        assert record.contigs[0].reads[1].rd.sequence[:10] == "aatacgGGAT"
        assert record.contigs[0].reads[1].rd.sequence[center - 5 : center + 5] == "ACATCATCTG"
        assert record.contigs[0].reads[1].rd.sequence[-10:] == "cATCTAGtac"
        assert record.contigs[0].reads[1].qa.qual_clipping_start == 7
        assert record.contigs[0].reads[1].qa.qual_clipping_end == 778
        assert record.contigs[0].reads[1].qa.align_clipping_start == 1
        assert record.contigs[0].reads[1].qa.align_clipping_end == 856
        assert record.contigs[0].reads[1].ds is None
        assert record.contigs[0].reads[1].rt is None
        assert record.contigs[0].reads[1].wr is None

        assert len(record.contigs[1].reads) == 14
        assert record.contigs[1].name == "Contig2"
        assert record.contigs[1].nbases == 3296
        assert record.contigs[1].nreads == 14
        assert record.contigs[1].nsegments == 214
        assert record.contigs[1].uorc == "U"
        center = len(record.contigs[1].sequence) // 2
        assert record.contigs[1].sequence[:10] == "cacggatgat"
        assert record.contigs[1].sequence[center - 5 : center + 5] == "TTTGAATATT"
        assert record.contigs[1].sequence[-10:] == "Atccttgtag"
        center = len(record.contigs[1].quality) // 2
        assert record.contigs[1].quality[:10] == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        assert record.contigs[1].quality[center - 5 : center + 5] == [90, 90, 90, 90, 90, 90, 90, 90, 90, 90]
        assert record.contigs[1].quality[-10:] == [24, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        assert len(record.contigs[1].af) == 14
        assert record.contigs[1].af[7].name == "BL060-LR3R.b.ab1"
        assert record.contigs[1].af[7].coru == "C"
        assert record.contigs[1].af[7].padded_start == 1601
        assert record.contigs[1].af[13].name == "BL060c2-LR0R.b.ab1"
        assert record.contigs[1].af[13].coru == "C"
        assert record.contigs[1].af[13].padded_start == 2445
        assert len(record.contigs[1].bs) == 214
        assert record.contigs[1].bs[107].name == "BL060-c1-LR3R.b.ab1"
        assert record.contigs[1].bs[107].padded_start == 2286
        assert record.contigs[1].bs[107].padded_end == 2292
        assert record.contigs[1].bs[213].name == "BL060c2-LR0R.b.ab1"
        assert record.contigs[1].bs[213].padded_start == 3236
        assert record.contigs[1].bs[213].padded_end == 3296
        assert len(record.contigs[1].ct) == 1
        assert record.contigs[1].ct[0].name == "Contig2"
        assert record.contigs[1].ct[0].tag_type == "repeat"
        assert record.contigs[1].ct[0].program == "phrap"
        assert record.contigs[1].ct[0].padded_start == 42
        assert record.contigs[1].ct[0].padded_end == 43
        assert record.contigs[1].ct[0].date == "123456:765432"
        assert record.contigs[1].ct[0].info == ["This is the first line of comment for c2", "and this the second for c2"]
        assert len(record.contigs[1].wa) == 1
        assert record.contigs[1].wa[0].tag_type == "phrap_params"
        assert record.contigs[1].wa[0].program == "phrap"
        assert record.contigs[1].wa[0].date == "040203:114710"
        assert record.contigs[1].wa[0].info == [
                "phrap 304_nuclsu.fasta.screen -new_ace -retain_duplicates",
                "phrap version 0.990329",
            ]

        assert len(record.contigs[1].reads) == 14

        # Read 0
        assert record.contigs[1].reads[0].rd.name == "BL060-c1-LR12.g.ab1"
        assert record.contigs[1].reads[0].rd.padded_bases == 862
        assert record.contigs[1].reads[0].rd.info_items == 0
        assert record.contigs[1].reads[0].rd.read_tags == 0
        center = len(record.contigs[1].reads[0].rd.sequence) // 2
        assert record.contigs[1].reads[0].rd.sequence[:10] == "cacggatgat"
        assert record.contigs[1].reads[0].rd.sequence[center - 5 : center + 5] == "GTTCTCGTTG"
        assert record.contigs[1].reads[0].rd.sequence[-10:] == "CGTTTACCcg"
        assert record.contigs[1].reads[0].qa.qual_clipping_start == 81
        assert record.contigs[1].reads[0].qa.qual_clipping_end == 842
        assert record.contigs[1].reads[0].qa.align_clipping_start == 1
        assert record.contigs[1].reads[0].qa.align_clipping_end == 862
        assert record.contigs[1].reads[0].ds.chromat_file == "BL060-c1-LR12.g.ab1"
        assert record.contigs[1].reads[0].ds.phd_file == "BL060-c1-LR12.g.ab1.phd.1"
        assert record.contigs[1].reads[0].ds.time == "Tue Feb  3 11:01:16 2004"
        assert record.contigs[1].reads[0].ds.chem == "term"
        assert record.contigs[1].reads[0].ds.dye == "big"
        assert record.contigs[1].reads[0].ds.template == ""
        assert record.contigs[1].reads[0].ds.direction == ""
        assert record.contigs[1].reads[0].rt is None
        assert record.contigs[1].reads[0].wr is None

        # Read 1
        assert record.contigs[1].reads[1].rd.name == "BL060-c1-LR11.g.ab1"
        assert record.contigs[1].reads[1].rd.padded_bases == 880
        assert record.contigs[1].reads[1].rd.info_items == 0
        assert record.contigs[1].reads[1].rd.read_tags == 0
        center = len(record.contigs[1].reads[1].rd.sequence) // 2
        assert record.contigs[1].reads[1].rd.sequence[:10] == "ctttctgacC"
        assert record.contigs[1].reads[1].rd.sequence[center - 5 : center + 5] == "CTGTGGTTTC"
        assert record.contigs[1].reads[1].rd.sequence[-10:] == "cggagttacg"
        assert record.contigs[1].reads[1].qa.qual_clipping_start == 11
        assert record.contigs[1].reads[1].qa.qual_clipping_end == 807
        assert record.contigs[1].reads[1].qa.align_clipping_start == 8
        assert record.contigs[1].reads[1].qa.align_clipping_end == 880
        assert record.contigs[1].reads[1].ds.chromat_file == "BL060-c1-LR11.g.ab1"
        assert record.contigs[1].reads[1].ds.phd_file == "BL060-c1-LR11.g.ab1.phd.1"
        assert record.contigs[1].reads[1].ds.time == "Tue Feb  3 11:01:16 2004"
        assert record.contigs[1].reads[1].ds.chem == "term"
        assert record.contigs[1].reads[1].ds.dye == "big"
        assert record.contigs[1].reads[1].ds.template == ""
        assert record.contigs[1].reads[1].ds.direction == ""
        assert len(record.contigs[1].reads[1].rt) == 0
        assert record.contigs[1].reads[1].wr is None

        # Read 2
        assert record.contigs[1].reads[2].rd.name == "BL060-c1-LR9.g.ab1"
        assert record.contigs[1].reads[2].rd.padded_bases == 864
        assert record.contigs[1].reads[2].rd.info_items == 0
        assert record.contigs[1].reads[2].rd.read_tags == 0
        center = len(record.contigs[1].reads[2].rd.sequence) // 2
        assert record.contigs[1].reads[2].rd.sequence[:10] == "cacccaCTTT"
        assert record.contigs[1].reads[2].rd.sequence[center - 5 : center + 5] == "ACCAAACATT"
        assert record.contigs[1].reads[2].rd.sequence[-10:] == "GGTAGCACgc"
        assert record.contigs[1].reads[2].qa.qual_clipping_start == 7
        assert record.contigs[1].reads[2].qa.qual_clipping_end == 840
        assert record.contigs[1].reads[2].qa.align_clipping_start == 4
        assert record.contigs[1].reads[2].qa.align_clipping_end == 864
        assert record.contigs[1].reads[2].ds.chromat_file == "BL060-c1-LR9.g.ab1"
        assert record.contigs[1].reads[2].ds.phd_file == "BL060-c1-LR9.g.ab1.phd.1"
        assert record.contigs[1].reads[2].ds.time == "Tue Feb  3 11:01:16 2004"
        assert record.contigs[1].reads[2].ds.chem == "term"
        assert record.contigs[1].reads[2].ds.dye == "big"
        assert record.contigs[1].reads[2].ds.template == ""
        assert record.contigs[1].reads[2].ds.direction == ""
        assert record.contigs[1].reads[2].rt is None
        assert record.contigs[1].reads[2].wr is None

        # Read 3
        assert record.contigs[1].reads[3].rd.name == "BL060-c1-LR17R.b.ab1"
        assert record.contigs[1].reads[3].rd.padded_bases == 863
        assert record.contigs[1].reads[3].rd.info_items == 0
        assert record.contigs[1].reads[3].rd.read_tags == 0
        center = len(record.contigs[1].reads[3].rd.sequence) // 2
        assert record.contigs[1].reads[3].rd.sequence[:10] == "ctaattggcc"
        assert record.contigs[1].reads[3].rd.sequence[center - 5 : center + 5] == "GGAACCTTTC"
        assert record.contigs[1].reads[3].rd.sequence[-10:] == "CAACCTgact"
        assert record.contigs[1].reads[3].qa.qual_clipping_start == 63
        assert record.contigs[1].reads[3].qa.qual_clipping_end == 857
        assert record.contigs[1].reads[3].qa.align_clipping_start == 1
        assert record.contigs[1].reads[3].qa.align_clipping_end == 861
        assert record.contigs[1].reads[3].ds.chromat_file == "BL060-c1-LR17R.b.ab1"
        assert record.contigs[1].reads[3].ds.phd_file == "BL060-c1-LR17R.b.ab1.phd.1"
        assert record.contigs[1].reads[3].ds.time == "Tue Feb  3 11:01:16 2004"
        assert record.contigs[1].reads[3].ds.chem == "term"
        assert record.contigs[1].reads[3].ds.dye == "big"
        assert record.contigs[1].reads[3].ds.template == ""
        assert record.contigs[1].reads[3].ds.direction == ""
        assert record.contigs[1].reads[3].rt == []
        assert record.contigs[1].reads[3].wr is None

        # Read 4
        assert record.contigs[1].reads[4].rd.name == "BL060-LR8.5.g.ab1"
        assert record.contigs[1].reads[4].rd.padded_bases == 877
        assert record.contigs[1].reads[4].rd.info_items == 0
        assert record.contigs[1].reads[4].rd.read_tags == 0
        center = len(record.contigs[1].reads[4].rd.sequence) // 2
        assert record.contigs[1].reads[4].rd.sequence[:10] == "tgCTGCGGTT"
        assert record.contigs[1].reads[4].rd.sequence[center - 5 : center + 5] == "GGCAGTTTCA"
        assert record.contigs[1].reads[4].rd.sequence[-10:] == "tactcataaa"
        assert record.contigs[1].reads[4].qa.qual_clipping_start == 13
        assert record.contigs[1].reads[4].qa.qual_clipping_end == 729
        assert record.contigs[1].reads[4].qa.align_clipping_start == 1
        assert record.contigs[1].reads[4].qa.align_clipping_end == 877
        assert record.contigs[1].reads[4].ds.chromat_file == "BL060-LR8.5.g.ab1"
        assert record.contigs[1].reads[4].ds.phd_file == "BL060-LR8.5.g.ab1.phd.1"
        assert record.contigs[1].reads[4].ds.time == "Fri Nov 14 09:46:03 2003"
        assert record.contigs[1].reads[4].ds.chem == "term"
        assert record.contigs[1].reads[4].ds.dye == "big"
        assert record.contigs[1].reads[4].ds.template == ""
        assert record.contigs[1].reads[4].ds.direction == ""
        assert record.contigs[1].reads[4].rt is None
        assert record.contigs[1].reads[4].wr is None

        # Read 5
        assert record.contigs[1].reads[5].rd.name == "BL060-LR3R.b.ab1"
        assert record.contigs[1].reads[5].rd.padded_bases == 874
        assert record.contigs[1].reads[5].rd.info_items == 0
        assert record.contigs[1].reads[5].rd.read_tags == 0
        center = len(record.contigs[1].reads[5].rd.sequence) // 2
        assert record.contigs[1].reads[5].rd.sequence[:10] == "ctCTTAGGAT"
        assert record.contigs[1].reads[5].rd.sequence[center - 5 : center + 5] == "AACTCACATT"
        assert record.contigs[1].reads[5].rd.sequence[-10:] == "*CACCCAAac"
        assert record.contigs[1].reads[5].qa.qual_clipping_start == 65
        assert record.contigs[1].reads[5].qa.qual_clipping_end == 874
        assert record.contigs[1].reads[5].qa.align_clipping_start == 1
        assert record.contigs[1].reads[5].qa.align_clipping_end == 874
        assert record.contigs[1].reads[5].ds.chromat_file == "BL060-LR3R.b.ab1"
        assert record.contigs[1].reads[5].ds.phd_file == "BL060-LR3R.b.ab1.phd.1"
        assert record.contigs[1].reads[5].ds.time == "Fri Nov 14 09:46:03 2003"
        assert record.contigs[1].reads[5].ds.chem == "term"
        assert record.contigs[1].reads[5].ds.dye == "big"
        assert record.contigs[1].reads[5].ds.template == ""
        assert record.contigs[1].reads[5].ds.direction == ""
        assert record.contigs[1].reads[5].rt is None
        assert record.contigs[1].reads[5].wr is None

        # Read 6
        assert record.contigs[1].reads[6].rd.name == "BL060-c1-LR3R.b.ab1"
        assert record.contigs[1].reads[6].rd.padded_bases == 864
        assert record.contigs[1].reads[6].rd.info_items == 0
        assert record.contigs[1].reads[6].rd.read_tags == 0
        center = len(record.contigs[1].reads[6].rd.sequence) // 2
        assert record.contigs[1].reads[6].rd.sequence[:10] == "CCaTGTCCAA"
        assert record.contigs[1].reads[6].rd.sequence[center - 5 : center + 5] == "AAGGGTT*CA"
        assert record.contigs[1].reads[6].rd.sequence[-10:] == "ACACTCGCga"
        assert record.contigs[1].reads[6].qa.qual_clipping_start == 73
        assert record.contigs[1].reads[6].qa.qual_clipping_end == 862
        assert record.contigs[1].reads[6].qa.align_clipping_start == 1
        assert record.contigs[1].reads[6].qa.align_clipping_end == 863
        assert record.contigs[1].reads[6].ds.chromat_file == "BL060-c1-LR3R.b.ab1"
        assert record.contigs[1].reads[6].ds.phd_file == "BL060-c1-LR3R.b.ab1.phd.1"
        assert record.contigs[1].reads[6].ds.time == "Tue Feb  3 11:01:16 2004"
        assert record.contigs[1].reads[6].ds.chem == "term"
        assert record.contigs[1].reads[6].ds.dye == "big"
        assert record.contigs[1].reads[6].ds.template == ""
        assert record.contigs[1].reads[6].ds.direction == ""
        assert record.contigs[1].reads[6].rt is None
        assert record.contigs[1].reads[6].wr is None

        # Read 7
        assert record.contigs[1].reads[7].rd.name == "BL060-LR3R.b.ab1"
        assert record.contigs[1].reads[7].rd.padded_bases == 857
        assert record.contigs[1].reads[7].rd.info_items == 0
        assert record.contigs[1].reads[7].rd.read_tags == 0
        center = len(record.contigs[1].reads[7].rd.sequence) // 2
        assert record.contigs[1].reads[7].rd.sequence[:10] == "agaaagagga"
        assert record.contigs[1].reads[7].rd.sequence[center - 5 : center + 5] == "nnnannnnnn"
        assert record.contigs[1].reads[7].rd.sequence[-10:] == "gtctttgctc"
        assert record.contigs[1].reads[7].qa.qual_clipping_start == 548
        assert record.contigs[1].reads[7].qa.qual_clipping_end == 847
        assert record.contigs[1].reads[7].qa.align_clipping_start == 442
        assert record.contigs[1].reads[7].qa.align_clipping_end == 854
        assert record.contigs[1].reads[7].ds.chromat_file == "BL060-LR3R.b.ab1"
        assert record.contigs[1].reads[7].ds.phd_file == "BL060-LR3R.b.ab1.phd.1"
        assert record.contigs[1].reads[7].ds.time == "Fri Jan 16 09:01:10 2004"
        assert record.contigs[1].reads[7].ds.chem == "term"
        assert record.contigs[1].reads[7].ds.dye == "big"
        assert record.contigs[1].reads[7].ds.template == ""
        assert record.contigs[1].reads[7].ds.direction == ""
        assert record.contigs[1].reads[7].rt is None
        assert record.contigs[1].reads[7].wr is None

        # Read 8
        assert record.contigs[1].reads[8].rd.name == "BL060-c1-LR7.g.ab1"
        assert record.contigs[1].reads[8].rd.padded_bases == 878
        assert record.contigs[1].reads[8].rd.info_items == 0
        assert record.contigs[1].reads[8].rd.read_tags == 0
        center = len(record.contigs[1].reads[8].rd.sequence) // 2
        assert record.contigs[1].reads[8].rd.sequence[:10] == "agTttc*ctc"
        assert record.contigs[1].reads[8].rd.sequence[center - 5 : center + 5] == "TCATAAAACT"
        assert record.contigs[1].reads[8].rd.sequence[-10:] == "xxxxxxxxxx"
        assert record.contigs[1].reads[8].qa.qual_clipping_start == 20
        assert record.contigs[1].reads[8].qa.qual_clipping_end == 798
        assert record.contigs[1].reads[8].qa.align_clipping_start == 1
        assert record.contigs[1].reads[8].qa.align_clipping_end == 798
        assert record.contigs[1].reads[8].ds.chromat_file == "BL060-c1-LR7.g.ab1"
        assert record.contigs[1].reads[8].ds.phd_file == "BL060-c1-LR7.g.ab1.phd.1"
        assert record.contigs[1].reads[8].ds.time == "Tue Feb  3 11:01:16 2004"
        assert record.contigs[1].reads[8].ds.chem == "term"
        assert record.contigs[1].reads[8].ds.dye == "big"
        assert record.contigs[1].reads[8].ds.template == ""
        assert record.contigs[1].reads[8].ds.direction == ""
        assert record.contigs[1].reads[8].rt is None
        assert record.contigs[1].reads[8].wr is None

        # Read 9
        assert record.contigs[1].reads[9].rd.name == "BL060-LR7.g.ab1"
        assert record.contigs[1].reads[9].rd.padded_bases == 880
        assert record.contigs[1].reads[9].rd.info_items == 0
        assert record.contigs[1].reads[9].rd.read_tags == 0
        center = len(record.contigs[1].reads[9].rd.sequence) // 2
        assert record.contigs[1].reads[9].rd.sequence[:10] == "ggctaCGCCc"
        assert record.contigs[1].reads[9].rd.sequence[center - 5 : center + 5] == "ATTGAGTTTC"
        assert record.contigs[1].reads[9].rd.sequence[-10:] == "tggcgttgcg"
        assert record.contigs[1].reads[9].qa.qual_clipping_start == 14
        assert record.contigs[1].reads[9].qa.qual_clipping_end == 765
        assert record.contigs[1].reads[9].qa.align_clipping_start == 4
        assert record.contigs[1].reads[9].qa.align_clipping_end == 765
        assert record.contigs[1].reads[9].ds.chromat_file == "BL060-LR7.g.ab1"
        assert record.contigs[1].reads[9].ds.phd_file == "BL060-LR7.g.ab1.phd.1"
        assert record.contigs[1].reads[9].ds.time == "Fri Nov 14 09:46:03 2003"
        assert record.contigs[1].reads[9].ds.chem == "term"
        assert record.contigs[1].reads[9].ds.dye == "big"
        assert record.contigs[1].reads[9].ds.template == ""
        assert record.contigs[1].reads[9].ds.direction == ""
        assert record.contigs[1].reads[9].rt is None
        assert record.contigs[1].reads[9].wr is None

        # Read 10
        assert record.contigs[1].reads[10].rd.name == "BL060c5-LR5.g.ab1"
        assert record.contigs[1].reads[10].rd.padded_bases == 871
        assert record.contigs[1].reads[10].rd.info_items == 0
        assert record.contigs[1].reads[10].rd.read_tags == 0
        center = len(record.contigs[1].reads[10].rd.sequence) // 2
        assert record.contigs[1].reads[10].rd.sequence[:10] == "ggtTCGATTA"
        assert record.contigs[1].reads[10].rd.sequence[center - 5 : center + 5] == "ACCAATTGAC"
        assert record.contigs[1].reads[10].rd.sequence[-10:] == "ACCACCCatt"
        assert record.contigs[1].reads[10].qa.qual_clipping_start == 12
        assert record.contigs[1].reads[10].qa.qual_clipping_end == 767
        assert record.contigs[1].reads[10].qa.align_clipping_start == 1
        assert record.contigs[1].reads[10].qa.align_clipping_end == 871
        assert record.contigs[1].reads[10].ds.chromat_file == "BL060c5-LR5.g.ab1"
        assert record.contigs[1].reads[10].ds.phd_file == "BL060c5-LR5.g.ab1.phd.1"
        assert record.contigs[1].reads[10].ds.time == "Fri Nov 14 09:46:03 2003"
        assert record.contigs[1].reads[10].ds.chem == "term"
        assert record.contigs[1].reads[10].ds.dye == "big"
        assert record.contigs[1].reads[10].ds.template == ""
        assert record.contigs[1].reads[10].ds.direction == ""
        assert record.contigs[1].reads[10].rt is None
        assert record.contigs[1].reads[10].wr is None

        # Read 11
        assert record.contigs[1].reads[11].rd.name == "BL060c2-LR5.g.ab1"
        assert record.contigs[1].reads[11].rd.padded_bases == 839
        assert record.contigs[1].reads[11].rd.info_items == 0
        assert record.contigs[1].reads[11].rd.read_tags == 0
        center = len(record.contigs[1].reads[11].rd.sequence) // 2
        assert record.contigs[1].reads[11].rd.sequence[:10] == "ggttcatatg"
        assert record.contigs[1].reads[11].rd.sequence[center - 5 : center + 5] == "TAAAATCAGT"
        assert record.contigs[1].reads[11].rd.sequence[-10:] == "TCTTGCaata"
        assert record.contigs[1].reads[11].qa.qual_clipping_start == 11
        assert record.contigs[1].reads[11].qa.qual_clipping_end == 757
        assert record.contigs[1].reads[11].qa.align_clipping_start == 10
        assert record.contigs[1].reads[11].qa.align_clipping_end == 835
        assert record.contigs[1].reads[11].ds is None
        assert len(record.contigs[1].reads[11].rt) == 1
        assert record.contigs[1].reads[11].rt[0].name == "BL060c2-LR5.g.ab1"
        assert record.contigs[1].reads[11].rt[0].tag_type == "matchElsewhereHighQual"
        assert record.contigs[1].reads[11].rt[0].program == "phrap"
        assert record.contigs[1].reads[11].rt[0].padded_start == 617
        assert record.contigs[1].reads[11].rt[0].padded_end == 631
        assert record.contigs[1].reads[11].rt[0].date == "040217:110357"
        assert record.contigs[1].reads[11].wr is None

        # Read 12
        assert record.contigs[1].reads[12].rd.name == "BL060c5-LR0R.b.ab1"
        assert record.contigs[1].reads[12].rd.padded_bases == 855
        assert record.contigs[1].reads[12].rd.info_items == 0
        assert record.contigs[1].reads[12].rd.read_tags == 0
        center = len(record.contigs[1].reads[12].rd.sequence) // 2
        assert record.contigs[1].reads[12].rd.sequence[:10] == "cACTCGCGTA"
        assert record.contigs[1].reads[12].rd.sequence[center - 5 : center + 5] == "CTCGTAAAAT"
        assert record.contigs[1].reads[12].rd.sequence[-10:] == "aacccctgca"
        assert record.contigs[1].reads[12].qa.qual_clipping_start == 94
        assert record.contigs[1].reads[12].qa.qual_clipping_end == 835
        assert record.contigs[1].reads[12].qa.align_clipping_start == 1
        assert record.contigs[1].reads[12].qa.align_clipping_end == 847
        assert record.contigs[1].reads[12].ds.chromat_file == "BL060c5-LR0R.b.ab1"
        assert record.contigs[1].reads[12].ds.phd_file == "BL060c5-LR0R.b.ab1.phd.1"
        assert record.contigs[1].reads[12].ds.time == "Wed Nov 12 08:16:30 2003"
        assert record.contigs[1].reads[12].ds.chem == "term"
        assert record.contigs[1].reads[12].ds.dye == "big"
        assert record.contigs[1].reads[12].ds.template == ""
        assert record.contigs[1].reads[12].ds.direction == ""
        assert len(record.contigs[1].reads[12].rt) == 1
        assert record.contigs[1].reads[12].rt[0].name == "BL060c5-LR0R.b.ab1"
        assert record.contigs[1].reads[12].rt[0].tag_type == "matchElsewhereHighQual"
        assert record.contigs[1].reads[12].rt[0].program == "phrap"
        assert record.contigs[1].reads[12].rt[0].padded_start == 617
        assert record.contigs[1].reads[12].rt[0].padded_end == 631
        assert record.contigs[1].reads[12].rt[0].date == "040217:110357"
        assert record.contigs[1].reads[12].wr is None

        # Read 13
        assert record.contigs[1].reads[13].rd.name == "BL060c2-LR0R.b.ab1"
        assert record.contigs[1].reads[13].rd.padded_bases == 852
        assert record.contigs[1].reads[13].rd.info_items == 0
        assert record.contigs[1].reads[13].rd.read_tags == 0
        center = len(record.contigs[1].reads[13].rd.sequence) // 2
        assert record.contigs[1].reads[13].rd.sequence[:10] == "cgCGTa*tTG"
        assert record.contigs[1].reads[13].rd.sequence[center - 5 : center + 5] == "GTAAAATATT"
        assert record.contigs[1].reads[13].rd.sequence[-10:] == "Atccttgtag"
        assert record.contigs[1].reads[13].qa.qual_clipping_start == 33
        assert record.contigs[1].reads[13].qa.qual_clipping_end == 831
        assert record.contigs[1].reads[13].qa.align_clipping_start == 1
        assert record.contigs[1].reads[13].qa.align_clipping_end == 852
        assert record.contigs[1].reads[13].ds.chromat_file == "BL060c2-LR0R.b.ab1"
        assert record.contigs[1].reads[13].ds.phd_file == "BL060c2-LR0R.b.ab1.phd.1"
        assert record.contigs[1].reads[13].ds.time == "Wed Nov 12 08:16:29 2003"
        assert record.contigs[1].reads[13].ds.chem == "term"
        assert record.contigs[1].reads[13].ds.dye == "big"
        assert record.contigs[1].reads[13].ds.template == ""
        assert record.contigs[1].reads[13].ds.direction == ""
        assert record.contigs[1].reads[13].rt == []
        assert len(record.contigs[1].reads[13].wr) == 1
        assert record.contigs[1].reads[13].wr[0].name == "BL060c2-LR0R.b.ab1"
        assert record.contigs[1].reads[13].wr[0].aligned == "unaligned"
        assert record.contigs[1].reads[13].wr[0].program == "phrap"
        assert record.contigs[1].reads[13].wr[0].date == "040217:110357"

    def test_check_record_parser(self):
        """Test to check that contig parser parses each contig into a contig."""
        contigs = Ace.parse(self.handle)

        # First contig
        contig = next(contigs)
        assert len(contig.reads) == 2
        assert contig.name == "Contig1"
        assert contig.nbases == 856
        assert contig.nreads == 2
        assert contig.nsegments == 31
        assert contig.uorc == "U"
        center = len(contig.sequence) // 2
        assert contig.sequence[:10] == "aatacgGGAT"
        assert contig.sequence[center - 5 : center + 5] == "ACATCATCTG"
        assert contig.sequence[-10:] == "cATCTAGtac"
        center = len(contig.quality) // 2
        assert contig.quality[:10] == [0, 0, 0, 0, 0, 0, 22, 23, 25, 28]
        assert contig.quality[center - 5 : center + 5] == [90, 90, 90, 90, 90, 90, 90, 90, 90, 90]
        assert contig.quality[-10:] == [15, 22, 30, 24, 28, 22, 21, 15, 19, 0]
        assert len(contig.af) == 2
        assert contig.af[1].name == "BL060c3-LR0R.b.ab1"
        assert contig.af[1].coru == "U"
        assert contig.af[1].padded_start == 1
        assert len(contig.bs) == 31
        assert contig.bs[15].name == "BL060c3-LR5.g.ab1"
        assert contig.bs[15].padded_start == 434
        assert contig.bs[15].padded_end == 438
        assert contig.bs[30].name == "BL060c3-LR0R.b.ab1"
        assert contig.bs[30].padded_start == 823
        assert contig.bs[30].padded_end == 856
        assert contig.ct is None
        assert contig.wa is None
        assert len(contig.reads) == 2
        assert contig.reads[0].rd.name == "BL060c3-LR5.g.ab1"
        assert contig.reads[0].rd.padded_bases == 868
        assert contig.reads[0].rd.info_items == 0
        assert contig.reads[0].rd.read_tags == 0
        center = len(contig.reads[0].rd.sequence) // 2
        assert contig.reads[0].rd.sequence[:10] == "tagcgaggaa"
        assert contig.reads[0].rd.sequence[center - 5 : center + 5] == "CCGAGGCCAA"
        assert contig.reads[0].rd.sequence[-10:] == "gaaccatcag"
        assert contig.reads[0].qa.qual_clipping_start == 80
        assert contig.reads[0].qa.qual_clipping_end == 853
        assert contig.reads[0].qa.align_clipping_start == 22
        assert contig.reads[0].qa.align_clipping_end == 856
        assert contig.reads[0].ds is None
        assert len(contig.reads[0].rt) == 2
        assert contig.reads[0].rt[0].name == "BL060c3-LR5.g.ab1"
        assert contig.reads[0].rt[0].tag_type == "matchElsewhereHighQual"
        assert contig.reads[0].rt[0].program == "phrap"
        assert contig.reads[0].rt[0].padded_start == 590
        assert contig.reads[0].rt[0].padded_end == 607
        assert contig.reads[0].rt[0].date == "040217:110357"
        assert contig.reads[0].rt[1].name == "BL060c3-LR5.g.ab1"
        assert contig.reads[0].rt[1].tag_type == "matchElsewhereHighQual"
        assert contig.reads[0].rt[1].program == "phrap"
        assert contig.reads[0].rt[1].padded_start == 617
        assert contig.reads[0].rt[1].padded_end == 631
        assert contig.reads[0].rt[1].date == "040217:110357"

        assert len(contig.reads[0].wr) == 1
        assert contig.reads[0].wr[0].name == "BL060c3-LR5.g.ab1"
        assert contig.reads[0].wr[0].aligned == "unaligned"
        assert contig.reads[0].wr[0].program == "phrap"
        assert contig.reads[0].wr[0].date == "040217:110357"

        assert contig.reads[1].rd.name == "BL060c3-LR0R.b.ab1"
        assert contig.reads[1].rd.padded_bases == 856
        assert contig.reads[1].rd.info_items == 0
        assert contig.reads[1].rd.read_tags == 0
        center = len(contig.reads[1].rd.sequence) // 2
        assert contig.reads[1].rd.sequence[:10] == "aatacgGGAT"
        assert contig.reads[1].rd.sequence[center - 5 : center + 5] == "ACATCATCTG"
        assert contig.reads[1].rd.sequence[-10:] == "cATCTAGtac"
        assert contig.reads[1].qa.qual_clipping_start == 7
        assert contig.reads[1].qa.qual_clipping_end == 778
        assert contig.reads[1].qa.align_clipping_start == 1
        assert contig.reads[1].qa.align_clipping_end == 856
        assert contig.reads[1].ds is None
        assert contig.reads[1].rt is None
        assert contig.reads[1].wr is None

        # Second contig
        contig = next(contigs)
        assert len(contig.reads) == 14
        assert contig.name == "Contig2"
        assert contig.nbases == 3296
        assert contig.nreads == 14
        assert contig.nsegments == 214
        assert contig.uorc == "U"
        center = len(contig.sequence) // 2
        assert contig.sequence[:10] == "cacggatgat"
        assert contig.sequence[center - 5 : center + 5] == "TTTGAATATT"
        assert contig.sequence[-10:] == "Atccttgtag"
        center = len(contig.quality) // 2
        assert contig.quality[:10] == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        assert contig.quality[center - 5 : center + 5] == [90, 90, 90, 90, 90, 90, 90, 90, 90, 90]
        assert contig.quality[-10:] == [24, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        assert len(contig.af) == 14
        assert contig.af[7].name == "BL060-LR3R.b.ab1"
        assert contig.af[7].coru == "C"
        assert contig.af[7].padded_start == 1601
        assert contig.af[13].name == "BL060c2-LR0R.b.ab1"
        assert contig.af[13].coru == "C"
        assert contig.af[13].padded_start == 2445
        assert len(contig.bs) == 214
        assert contig.bs[107].name == "BL060-c1-LR3R.b.ab1"
        assert contig.bs[107].padded_start == 2286
        assert contig.bs[107].padded_end == 2292
        assert contig.bs[213].name == "BL060c2-LR0R.b.ab1"
        assert contig.bs[213].padded_start == 3236
        assert contig.bs[213].padded_end == 3296
        assert len(contig.ct) == 3
        assert contig.ct[0].name == "Contig2"
        assert contig.ct[0].tag_type == "repeat"
        assert contig.ct[0].program == "phrap"
        assert contig.ct[0].padded_start == 42
        assert contig.ct[0].padded_end == 43
        assert contig.ct[0].date == "123456:765432"
        assert contig.ct[0].info == ["This is the first line of comment for c2", "and this the second for c2"]
        assert contig.ct[1].name == "unrelated_Contig"
        assert contig.ct[1].tag_type == "repeat"
        assert contig.ct[1].program == "phrap"
        assert contig.ct[1].padded_start == 1142
        assert contig.ct[1].padded_end == 143
        assert contig.ct[1].date == "122226:722232"
        assert contig.ct[1].info == [
                "This is the first line of comment for the unrelated ct tag",
                "and this the second",
            ]

        assert contig.ct[2].name == "Contig1"
        assert contig.ct[2].tag_type == "repeat"
        assert contig.ct[2].program == "phrap"
        assert contig.ct[2].padded_start == 52
        assert contig.ct[2].padded_end == 53
        assert contig.ct[2].date == "555456:555432"
        assert contig.ct[2].info == ["This is the first line of comment for c1", "and this the second for c1"]

        assert len(contig.wa) == 1
        assert contig.wa[0].tag_type == "phrap_params"
        assert contig.wa[0].program == "phrap"
        assert contig.wa[0].date == "040203:114710"
        assert contig.wa[0].info == [
                "phrap 304_nuclsu.fasta.screen -new_ace -retain_duplicates",
                "phrap version 0.990329",
            ]

        assert len(contig.reads) == 14

        # Read 0
        assert contig.reads[0].rd.name == "BL060-c1-LR12.g.ab1"
        assert contig.reads[0].rd.padded_bases == 862
        assert contig.reads[0].rd.info_items == 0
        assert contig.reads[0].rd.read_tags == 0
        center = len(contig.reads[0].rd.sequence) // 2
        assert contig.reads[0].rd.sequence[:10] == "cacggatgat"
        assert contig.reads[0].rd.sequence[center - 5 : center + 5] == "GTTCTCGTTG"
        assert contig.reads[0].rd.sequence[-10:] == "CGTTTACCcg"
        assert contig.reads[0].qa.qual_clipping_start == 81
        assert contig.reads[0].qa.qual_clipping_end == 842
        assert contig.reads[0].qa.align_clipping_start == 1
        assert contig.reads[0].qa.align_clipping_end == 862
        assert contig.reads[0].ds.chromat_file == "BL060-c1-LR12.g.ab1"
        assert contig.reads[0].ds.phd_file == "BL060-c1-LR12.g.ab1.phd.1"
        assert contig.reads[0].ds.time == "Tue Feb  3 11:01:16 2004"
        assert contig.reads[0].ds.chem == "term"
        assert contig.reads[0].ds.dye == "big"
        assert contig.reads[0].ds.template == ""
        assert contig.reads[0].ds.direction == ""
        assert contig.reads[0].rt is None
        assert contig.reads[0].wr is None

        # Read 1
        assert contig.reads[1].rd.name == "BL060-c1-LR11.g.ab1"
        assert contig.reads[1].rd.padded_bases == 880
        assert contig.reads[1].rd.info_items == 0
        assert contig.reads[1].rd.read_tags == 0
        center = len(contig.reads[1].rd.sequence) // 2
        assert contig.reads[1].rd.sequence[:10] == "ctttctgacC"
        assert contig.reads[1].rd.sequence[center - 5 : center + 5] == "CTGTGGTTTC"
        assert contig.reads[1].rd.sequence[-10:] == "cggagttacg"
        assert contig.reads[1].qa.qual_clipping_start == 11
        assert contig.reads[1].qa.qual_clipping_end == 807
        assert contig.reads[1].qa.align_clipping_start == 8
        assert contig.reads[1].qa.align_clipping_end == 880
        assert contig.reads[1].ds.chromat_file == "BL060-c1-LR11.g.ab1"
        assert contig.reads[1].ds.phd_file == "BL060-c1-LR11.g.ab1.phd.1"
        assert contig.reads[1].ds.time == "Tue Feb  3 11:01:16 2004"
        assert contig.reads[1].ds.chem == "term"
        assert contig.reads[1].ds.dye == "big"
        assert contig.reads[1].ds.template == ""
        assert contig.reads[1].ds.direction == ""
        assert len(contig.reads[1].rt) == 1
        assert contig.reads[1].rt[0].name == "BL060c3-LR5.g.ab1"
        assert contig.reads[1].rt[0].tag_type == "matchElsewhereHighQual"
        assert contig.reads[1].rt[0].program == "phrap"
        assert contig.reads[1].rt[0].padded_start == 617
        assert contig.reads[1].rt[0].padded_end == 631
        assert contig.reads[1].rt[0].date == "040217:110357"
        assert contig.reads[1].wr is None

        # Read 2
        assert contig.reads[2].rd.name == "BL060-c1-LR9.g.ab1"
        assert contig.reads[2].rd.padded_bases == 864
        assert contig.reads[2].rd.info_items == 0
        assert contig.reads[2].rd.read_tags == 0
        center = len(contig.reads[2].rd.sequence) // 2
        assert contig.reads[2].rd.sequence[:10] == "cacccaCTTT"
        assert contig.reads[2].rd.sequence[center - 5 : center + 5] == "ACCAAACATT"
        assert contig.reads[2].rd.sequence[-10:] == "GGTAGCACgc"
        assert contig.reads[2].qa.qual_clipping_start == 7
        assert contig.reads[2].qa.qual_clipping_end == 840
        assert contig.reads[2].qa.align_clipping_start == 4
        assert contig.reads[2].qa.align_clipping_end == 864
        assert contig.reads[2].ds.chromat_file == "BL060-c1-LR9.g.ab1"
        assert contig.reads[2].ds.phd_file == "BL060-c1-LR9.g.ab1.phd.1"
        assert contig.reads[2].ds.time == "Tue Feb  3 11:01:16 2004"
        assert contig.reads[2].ds.chem == "term"
        assert contig.reads[2].ds.dye == "big"
        assert contig.reads[2].ds.template == ""
        assert contig.reads[2].ds.direction == ""
        assert contig.reads[2].rt is None
        assert contig.reads[2].wr is None

        # Read 3
        assert contig.reads[3].rd.name == "BL060-c1-LR17R.b.ab1"
        assert contig.reads[3].rd.padded_bases == 863
        assert contig.reads[3].rd.info_items == 0
        assert contig.reads[3].rd.read_tags == 0
        center = len(contig.reads[3].rd.sequence) // 2
        assert contig.reads[3].rd.sequence[:10] == "ctaattggcc"
        assert contig.reads[3].rd.sequence[center - 5 : center + 5] == "GGAACCTTTC"
        assert contig.reads[3].rd.sequence[-10:] == "CAACCTgact"
        assert contig.reads[3].qa.qual_clipping_start == 63
        assert contig.reads[3].qa.qual_clipping_end == 857
        assert contig.reads[3].qa.align_clipping_start == 1
        assert contig.reads[3].qa.align_clipping_end == 861
        assert contig.reads[3].ds.chromat_file == "BL060-c1-LR17R.b.ab1"
        assert contig.reads[3].ds.phd_file == "BL060-c1-LR17R.b.ab1.phd.1"
        assert contig.reads[3].ds.time == "Tue Feb  3 11:01:16 2004"
        assert contig.reads[3].ds.chem == "term"
        assert contig.reads[3].ds.dye == "big"
        assert contig.reads[3].ds.template == ""
        assert contig.reads[3].ds.direction == ""
        assert len(contig.reads[3].rt) == 1
        assert contig.reads[3].rt[0].name == "BL060c3-LR5.g.ab1"
        assert contig.reads[3].rt[0].tag_type == "matchElsewhereHighQual"
        assert contig.reads[3].rt[0].program == "phrap"
        assert contig.reads[3].rt[0].padded_start == 617
        assert contig.reads[3].rt[0].padded_end == 631
        assert contig.reads[3].rt[0].date == "040217:110357"
        assert contig.reads[3].wr is None

        # Read 4
        assert contig.reads[4].rd.name == "BL060-LR8.5.g.ab1"
        assert contig.reads[4].rd.padded_bases == 877
        assert contig.reads[4].rd.info_items == 0
        assert contig.reads[4].rd.read_tags == 0
        center = len(contig.reads[4].rd.sequence) // 2
        assert contig.reads[4].rd.sequence[:10] == "tgCTGCGGTT"
        assert contig.reads[4].rd.sequence[center - 5 : center + 5] == "GGCAGTTTCA"
        assert contig.reads[4].rd.sequence[-10:] == "tactcataaa"
        assert contig.reads[4].qa.qual_clipping_start == 13
        assert contig.reads[4].qa.qual_clipping_end == 729
        assert contig.reads[4].qa.align_clipping_start == 1
        assert contig.reads[4].qa.align_clipping_end == 877
        assert contig.reads[4].ds.chromat_file == "BL060-LR8.5.g.ab1"
        assert contig.reads[4].ds.phd_file == "BL060-LR8.5.g.ab1.phd.1"
        assert contig.reads[4].ds.time == "Fri Nov 14 09:46:03 2003"
        assert contig.reads[4].ds.chem == "term"
        assert contig.reads[4].ds.dye == "big"
        assert contig.reads[4].ds.template == ""
        assert contig.reads[4].ds.direction == ""
        assert contig.reads[4].rt is None
        assert contig.reads[4].wr is None

        # Read 5
        assert contig.reads[5].rd.name == "BL060-LR3R.b.ab1"
        assert contig.reads[5].rd.padded_bases == 874
        assert contig.reads[5].rd.info_items == 0
        assert contig.reads[5].rd.read_tags == 0
        center = len(contig.reads[5].rd.sequence) // 2
        assert contig.reads[5].rd.sequence[:10] == "ctCTTAGGAT"
        assert contig.reads[5].rd.sequence[center - 5 : center + 5] == "AACTCACATT"
        assert contig.reads[5].rd.sequence[-10:] == "*CACCCAAac"
        assert contig.reads[5].qa.qual_clipping_start == 65
        assert contig.reads[5].qa.qual_clipping_end == 874
        assert contig.reads[5].qa.align_clipping_start == 1
        assert contig.reads[5].qa.align_clipping_end == 874
        assert contig.reads[5].ds.chromat_file == "BL060-LR3R.b.ab1"
        assert contig.reads[5].ds.phd_file == "BL060-LR3R.b.ab1.phd.1"
        assert contig.reads[5].ds.time == "Fri Nov 14 09:46:03 2003"
        assert contig.reads[5].ds.chem == "term"
        assert contig.reads[5].ds.dye == "big"
        assert contig.reads[5].ds.template == ""
        assert contig.reads[5].ds.direction == ""
        assert contig.reads[5].rt is None
        assert contig.reads[5].wr is None

        # Read 6
        assert contig.reads[6].rd.name == "BL060-c1-LR3R.b.ab1"
        assert contig.reads[6].rd.padded_bases == 864
        assert contig.reads[6].rd.info_items == 0
        assert contig.reads[6].rd.read_tags == 0
        center = len(contig.reads[6].rd.sequence) // 2
        assert contig.reads[6].rd.sequence[:10] == "CCaTGTCCAA"
        assert contig.reads[6].rd.sequence[center - 5 : center + 5] == "AAGGGTT*CA"
        assert contig.reads[6].rd.sequence[-10:] == "ACACTCGCga"
        assert contig.reads[6].qa.qual_clipping_start == 73
        assert contig.reads[6].qa.qual_clipping_end == 862
        assert contig.reads[6].qa.align_clipping_start == 1
        assert contig.reads[6].qa.align_clipping_end == 863
        assert contig.reads[6].ds.chromat_file == "BL060-c1-LR3R.b.ab1"
        assert contig.reads[6].ds.phd_file == "BL060-c1-LR3R.b.ab1.phd.1"
        assert contig.reads[6].ds.time == "Tue Feb  3 11:01:16 2004"
        assert contig.reads[6].ds.chem == "term"
        assert contig.reads[6].ds.dye == "big"
        assert contig.reads[6].ds.template == ""
        assert contig.reads[6].ds.direction == ""
        assert contig.reads[6].rt is None
        assert contig.reads[6].wr is None

        # Read 7
        assert contig.reads[7].rd.name == "BL060-LR3R.b.ab1"
        assert contig.reads[7].rd.padded_bases == 857
        assert contig.reads[7].rd.info_items == 0
        assert contig.reads[7].rd.read_tags == 0
        center = len(contig.reads[7].rd.sequence) // 2
        assert contig.reads[7].rd.sequence[:10] == "agaaagagga"
        assert contig.reads[7].rd.sequence[center - 5 : center + 5] == "nnnannnnnn"
        assert contig.reads[7].rd.sequence[-10:] == "gtctttgctc"
        assert contig.reads[7].qa.qual_clipping_start == 548
        assert contig.reads[7].qa.qual_clipping_end == 847
        assert contig.reads[7].qa.align_clipping_start == 442
        assert contig.reads[7].qa.align_clipping_end == 854
        assert contig.reads[7].ds.chromat_file == "BL060-LR3R.b.ab1"
        assert contig.reads[7].ds.phd_file == "BL060-LR3R.b.ab1.phd.1"
        assert contig.reads[7].ds.time == "Fri Jan 16 09:01:10 2004"
        assert contig.reads[7].ds.chem == "term"
        assert contig.reads[7].ds.dye == "big"
        assert contig.reads[7].ds.template == ""
        assert contig.reads[7].ds.direction == ""
        assert contig.reads[7].rt is None
        assert contig.reads[7].wr is None

        # Read 8
        assert contig.reads[8].rd.name == "BL060-c1-LR7.g.ab1"
        assert contig.reads[8].rd.padded_bases == 878
        assert contig.reads[8].rd.info_items == 0
        assert contig.reads[8].rd.read_tags == 0
        center = len(contig.reads[8].rd.sequence) // 2
        assert contig.reads[8].rd.sequence[:10] == "agTttc*ctc"
        assert contig.reads[8].rd.sequence[center - 5 : center + 5] == "TCATAAAACT"
        assert contig.reads[8].rd.sequence[-10:] == "xxxxxxxxxx"
        assert contig.reads[8].qa.qual_clipping_start == 20
        assert contig.reads[8].qa.qual_clipping_end == 798
        assert contig.reads[8].qa.align_clipping_start == 1
        assert contig.reads[8].qa.align_clipping_end == 798
        assert contig.reads[8].ds.chromat_file == "BL060-c1-LR7.g.ab1"
        assert contig.reads[8].ds.phd_file == "BL060-c1-LR7.g.ab1.phd.1"
        assert contig.reads[8].ds.time == "Tue Feb  3 11:01:16 2004"
        assert contig.reads[8].ds.chem == "term"
        assert contig.reads[8].ds.dye == "big"
        assert contig.reads[8].ds.template == ""
        assert contig.reads[8].ds.direction == ""
        assert contig.reads[8].rt is None
        assert contig.reads[8].wr is None

        # Read 9
        assert contig.reads[9].rd.name == "BL060-LR7.g.ab1"
        assert contig.reads[9].rd.padded_bases == 880
        assert contig.reads[9].rd.info_items == 0
        assert contig.reads[9].rd.read_tags == 0
        center = len(contig.reads[9].rd.sequence) // 2
        assert contig.reads[9].rd.sequence[:10] == "ggctaCGCCc"
        assert contig.reads[9].rd.sequence[center - 5 : center + 5] == "ATTGAGTTTC"
        assert contig.reads[9].rd.sequence[-10:] == "tggcgttgcg"
        assert contig.reads[9].qa.qual_clipping_start == 14
        assert contig.reads[9].qa.qual_clipping_end == 765
        assert contig.reads[9].qa.align_clipping_start == 4
        assert contig.reads[9].qa.align_clipping_end == 765
        assert contig.reads[9].ds.chromat_file == "BL060-LR7.g.ab1"
        assert contig.reads[9].ds.phd_file == "BL060-LR7.g.ab1.phd.1"
        assert contig.reads[9].ds.time == "Fri Nov 14 09:46:03 2003"
        assert contig.reads[9].ds.chem == "term"
        assert contig.reads[9].ds.dye == "big"
        assert contig.reads[9].ds.template == ""
        assert contig.reads[9].ds.direction == ""
        assert contig.reads[9].rt is None
        assert contig.reads[9].wr is None

        # Read 10
        assert contig.reads[10].rd.name == "BL060c5-LR5.g.ab1"
        assert contig.reads[10].rd.padded_bases == 871
        assert contig.reads[10].rd.info_items == 0
        assert contig.reads[10].rd.read_tags == 0
        center = len(contig.reads[10].rd.sequence) // 2
        assert contig.reads[10].rd.sequence[:10] == "ggtTCGATTA"
        assert contig.reads[10].rd.sequence[center - 5 : center + 5] == "ACCAATTGAC"
        assert contig.reads[10].rd.sequence[-10:] == "ACCACCCatt"
        assert contig.reads[10].qa.qual_clipping_start == 12
        assert contig.reads[10].qa.qual_clipping_end == 767
        assert contig.reads[10].qa.align_clipping_start == 1
        assert contig.reads[10].qa.align_clipping_end == 871
        assert contig.reads[10].ds.chromat_file == "BL060c5-LR5.g.ab1"
        assert contig.reads[10].ds.phd_file == "BL060c5-LR5.g.ab1.phd.1"
        assert contig.reads[10].ds.time == "Fri Nov 14 09:46:03 2003"
        assert contig.reads[10].ds.chem == "term"
        assert contig.reads[10].ds.dye == "big"
        assert contig.reads[10].ds.template == ""
        assert contig.reads[10].ds.direction == ""
        assert contig.reads[10].rt is None
        assert contig.reads[10].wr is None

        # Read 11
        assert contig.reads[11].rd.name == "BL060c2-LR5.g.ab1"
        assert contig.reads[11].rd.padded_bases == 839
        assert contig.reads[11].rd.info_items == 0
        assert contig.reads[11].rd.read_tags == 0
        center = len(contig.reads[11].rd.sequence) // 2
        assert contig.reads[11].rd.sequence[:10] == "ggttcatatg"
        assert contig.reads[11].rd.sequence[center - 5 : center + 5] == "TAAAATCAGT"
        assert contig.reads[11].rd.sequence[-10:] == "TCTTGCaata"
        assert contig.reads[11].qa.qual_clipping_start == 11
        assert contig.reads[11].qa.qual_clipping_end == 757
        assert contig.reads[11].qa.align_clipping_start == 10
        assert contig.reads[11].qa.align_clipping_end == 835
        assert contig.reads[11].ds is None
        assert len(contig.reads[11].rt) == 1
        assert contig.reads[11].rt[0].name == "BL060c2-LR5.g.ab1"
        assert contig.reads[11].rt[0].tag_type == "matchElsewhereHighQual"
        assert contig.reads[11].rt[0].program == "phrap"
        assert contig.reads[11].rt[0].padded_start == 617
        assert contig.reads[11].rt[0].padded_end == 631
        assert contig.reads[11].rt[0].date == "040217:110357"
        assert contig.reads[11].wr is None

        # Read 12
        assert contig.reads[12].rd.name == "BL060c5-LR0R.b.ab1"
        assert contig.reads[12].rd.padded_bases == 855
        assert contig.reads[12].rd.info_items == 0
        assert contig.reads[12].rd.read_tags == 0
        center = len(contig.reads[12].rd.sequence) // 2
        assert contig.reads[12].rd.sequence[:10] == "cACTCGCGTA"
        assert contig.reads[12].rd.sequence[center - 5 : center + 5] == "CTCGTAAAAT"
        assert contig.reads[12].rd.sequence[-10:] == "aacccctgca"
        assert contig.reads[12].qa.qual_clipping_start == 94
        assert contig.reads[12].qa.qual_clipping_end == 835
        assert contig.reads[12].qa.align_clipping_start == 1
        assert contig.reads[12].qa.align_clipping_end == 847
        assert contig.reads[12].ds.chromat_file == "BL060c5-LR0R.b.ab1"
        assert contig.reads[12].ds.phd_file == "BL060c5-LR0R.b.ab1.phd.1"
        assert contig.reads[12].ds.time == "Wed Nov 12 08:16:30 2003"
        assert contig.reads[12].ds.chem == "term"
        assert contig.reads[12].ds.dye == "big"
        assert contig.reads[12].ds.template == ""
        assert contig.reads[12].ds.direction == ""
        assert contig.reads[12].rt is None
        assert contig.reads[12].wr is None

        # Read 13
        assert contig.reads[13].rd.name == "BL060c2-LR0R.b.ab1"
        assert contig.reads[13].rd.padded_bases == 852
        assert contig.reads[13].rd.info_items == 0
        assert contig.reads[13].rd.read_tags == 0
        center = len(contig.reads[13].rd.sequence) // 2
        assert contig.reads[13].rd.sequence[:10] == "cgCGTa*tTG"
        assert contig.reads[13].rd.sequence[center - 5 : center + 5] == "GTAAAATATT"
        assert contig.reads[13].rd.sequence[-10:] == "Atccttgtag"
        assert contig.reads[13].qa.qual_clipping_start == 33
        assert contig.reads[13].qa.qual_clipping_end == 831
        assert contig.reads[13].qa.align_clipping_start == 1
        assert contig.reads[13].qa.align_clipping_end == 852
        assert contig.reads[13].ds.chromat_file == "BL060c2-LR0R.b.ab1"
        assert contig.reads[13].ds.phd_file == "BL060c2-LR0R.b.ab1.phd.1"
        assert contig.reads[13].ds.time == "Wed Nov 12 08:16:29 2003"
        assert contig.reads[13].ds.chem == "term"
        assert contig.reads[13].ds.dye == "big"
        assert contig.reads[13].ds.template == ""
        assert contig.reads[13].ds.direction == ""
        assert len(contig.reads[13].rt) == 1
        assert contig.reads[13].rt[0].name == "BL060c5-LR0R.b.ab1"
        assert contig.reads[13].rt[0].tag_type == "matchElsewhereHighQual"
        assert contig.reads[13].rt[0].program == "phrap"
        assert contig.reads[13].rt[0].padded_start == 617
        assert contig.reads[13].rt[0].padded_end == 631
        assert contig.reads[13].rt[0].date == "040217:110357"
        assert len(contig.reads[13].wr) == 1
        assert contig.reads[13].wr[0].name == "BL060c2-LR0R.b.ab1"
        assert contig.reads[13].wr[0].aligned == "unaligned"
        assert contig.reads[13].wr[0].program == "phrap"
        assert contig.reads[13].wr[0].date == "040217:110357"

        # Make sure there are no more contigs
        with pytest.raises(StopIteration):
            next(contigs)


class TestAceTestTwo:
    """Test parsing example output from CAP3.

    The sample input file seq.cap.ace was downloaded from:
    http://genome.cs.mtu.edu/cap/data/seq.cap.ace
    """

    @pytest.fixture(autouse=True)
    def _setup(self):
        self.handle = open("Ace/seq.cap.ace")

        yield
        self.handle.close()

    def test_check_ACEParser(self):
        """Test to check that ACEParser can parse the whole file into one record."""
        record = Ace.read(self.handle)
        assert record.ncontigs == 1
        assert record.nreads == 6
        assert record.wa is None
        assert len(record.contigs) == 1

        assert len(record.contigs[0].reads) == 6
        assert record.contigs[0].name == "Contig1"
        assert record.contigs[0].nbases == 1222
        assert record.contigs[0].nreads == 6
        assert record.contigs[0].nsegments == 0
        assert record.contigs[0].uorc == "U"
        center = len(record.contigs[0].sequence) // 2
        assert record.contigs[0].sequence[:10] == "AGTTTTAGTT"
        assert record.contigs[0].sequence[center - 5 : center + 5] == "TGTGCGCGCA"
        assert record.contigs[0].sequence[-10:] == "ATATCACATT"
        center = len(record.contigs[0].quality) // 2
        assert record.contigs[0].quality[:10] == [61, 66, 67, 70, 71, 73, 73, 77, 77, 87]
        assert record.contigs[0].quality[center - 5 : center + 5] == [97, 97, 97, 97, 97, 97, 97, 97, 97, 97]
        assert record.contigs[0].quality[-10:] == [56, 51, 49, 41, 38, 39, 45, 44, 49, 46]
        assert len(record.contigs[0].af) == 6
        assert len(record.contigs[0].bs) == 0
        assert record.contigs[0].af[3].name == "R5"
        assert record.contigs[0].af[3].coru == "C"
        assert record.contigs[0].af[3].padded_start == 320
        assert record.contigs[0].af[5].name == "R6"
        assert record.contigs[0].af[5].coru == "C"
        assert record.contigs[0].af[5].padded_start == 517
        assert record.contigs[0].bs == []
        assert record.contigs[0].ct is None
        assert record.contigs[0].wa is None
        assert len(record.contigs[0].reads) == 6

        assert record.contigs[0].reads[0].rd.name == "R3"
        assert record.contigs[0].reads[0].rd.padded_bases == 919
        assert record.contigs[0].reads[0].rd.info_items == 0
        assert record.contigs[0].reads[0].rd.read_tags == 0
        center = len(record.contigs[0].reads[0].rd.sequence) // 2
        assert record.contigs[0].reads[0].rd.sequence[:10] == "NNNNNNNNNN"
        assert record.contigs[0].reads[0].rd.sequence[center - 5 : center + 5] == "ATGTGCGCTC"
        assert record.contigs[0].reads[0].rd.sequence[-10:] == "CAGCTCACCA"
        assert record.contigs[0].reads[0].qa.qual_clipping_start == 55
        assert record.contigs[0].reads[0].qa.qual_clipping_end == 916
        assert record.contigs[0].reads[0].qa.align_clipping_start == 55
        assert record.contigs[0].reads[0].qa.align_clipping_end == 916
        assert record.contigs[0].reads[0].ds.chromat_file == ""
        assert record.contigs[0].reads[0].ds.phd_file == ""
        assert record.contigs[0].reads[0].ds.time == ""
        assert record.contigs[0].reads[0].ds.chem == ""
        assert record.contigs[0].reads[0].ds.dye == ""
        assert record.contigs[0].reads[0].ds.template == ""
        assert record.contigs[0].reads[0].ds.direction == ""
        assert record.contigs[0].reads[0].rt is None
        assert record.contigs[0].reads[0].wr is None

        assert record.contigs[0].reads[1].rd.name == "R1"
        assert record.contigs[0].reads[1].rd.padded_bases == 864
        assert record.contigs[0].reads[1].rd.info_items == 0
        assert record.contigs[0].reads[1].rd.read_tags == 0
        center = len(record.contigs[0].reads[1].rd.sequence) // 2
        assert record.contigs[0].reads[1].rd.sequence[:10] == "AGCCGGTACC"
        assert record.contigs[0].reads[1].rd.sequence[center - 5 : center + 5] == "GGGATGGCAC"
        assert record.contigs[0].reads[1].rd.sequence[-10:] == "GGGCTGGGAG"
        assert record.contigs[0].reads[1].qa.qual_clipping_start == 12
        assert record.contigs[0].reads[1].qa.qual_clipping_end == 863
        assert record.contigs[0].reads[1].qa.align_clipping_start == 12
        assert record.contigs[0].reads[1].qa.align_clipping_end == 863
        assert record.contigs[0].reads[1].ds.chromat_file == ""
        assert record.contigs[0].reads[1].ds.phd_file == ""
        assert record.contigs[0].reads[1].ds.time == ""
        assert record.contigs[0].reads[1].ds.chem == ""
        assert record.contigs[0].reads[1].ds.dye == ""
        assert record.contigs[0].reads[1].ds.template == ""
        assert record.contigs[0].reads[1].ds.direction == ""
        assert record.contigs[0].reads[1].rt is None
        assert record.contigs[0].reads[1].wr is None

        assert record.contigs[0].reads[2].rd.name == "R2"
        assert record.contigs[0].reads[2].rd.padded_bases == 1026
        assert record.contigs[0].reads[2].rd.info_items == 0
        assert record.contigs[0].reads[2].rd.read_tags == 0
        center = len(record.contigs[0].reads[2].rd.sequence) // 2
        assert record.contigs[0].reads[2].rd.sequence[:10] == "NNNNNNNNNN"
        assert record.contigs[0].reads[2].rd.sequence[center - 5 : center + 5] == "GGATGCCTGG"
        assert record.contigs[0].reads[2].rd.sequence[-10:] == "GGTTGAGGCC"
        assert record.contigs[0].reads[2].qa.qual_clipping_start == 55
        assert record.contigs[0].reads[2].qa.qual_clipping_end == 1000
        assert record.contigs[0].reads[2].qa.align_clipping_start == 55
        assert record.contigs[0].reads[2].qa.align_clipping_end == 1000
        assert record.contigs[0].reads[2].ds.chromat_file == ""
        assert record.contigs[0].reads[2].ds.phd_file == ""
        assert record.contigs[0].reads[2].ds.time == ""
        assert record.contigs[0].reads[2].ds.chem == ""
        assert record.contigs[0].reads[2].ds.dye == ""
        assert record.contigs[0].reads[2].ds.template == ""
        assert record.contigs[0].reads[2].ds.direction == ""
        assert record.contigs[0].reads[2].rt is None
        assert record.contigs[0].reads[2].wr is None

        assert record.contigs[0].reads[3].rd.name == "R5"
        assert record.contigs[0].reads[3].rd.padded_bases == 925
        assert record.contigs[0].reads[3].rd.info_items == 0
        assert record.contigs[0].reads[3].rd.read_tags == 0
        center = len(record.contigs[0].reads[3].rd.sequence) // 2
        assert record.contigs[0].reads[3].rd.sequence[:10] == "NNNNNNNNNN"
        assert record.contigs[0].reads[3].rd.sequence[center - 5 : center + 5] == "CCTCCCTACA"
        assert record.contigs[0].reads[3].rd.sequence[-10:] == "GCCCCCGGNN"
        assert record.contigs[0].reads[3].qa.qual_clipping_start == 293
        assert record.contigs[0].reads[3].qa.qual_clipping_end == 874
        assert record.contigs[0].reads[3].qa.align_clipping_start == 293
        assert record.contigs[0].reads[3].qa.align_clipping_end == 874
        assert record.contigs[0].reads[3].ds.chromat_file == ""
        assert record.contigs[0].reads[3].ds.phd_file == ""
        assert record.contigs[0].reads[3].ds.time == ""
        assert record.contigs[0].reads[3].ds.chem == ""
        assert record.contigs[0].reads[3].ds.dye == ""
        assert record.contigs[0].reads[3].ds.template == ""
        assert record.contigs[0].reads[3].ds.direction == ""
        assert record.contigs[0].reads[3].rt is None
        assert record.contigs[0].reads[3].wr is None

        assert record.contigs[0].reads[4].rd.name == "R4"
        assert record.contigs[0].reads[4].rd.padded_bases == 816
        assert record.contigs[0].reads[4].rd.info_items == 0
        assert record.contigs[0].reads[4].rd.read_tags == 0
        center = len(record.contigs[0].reads[4].rd.sequence) // 2
        assert record.contigs[0].reads[4].rd.sequence[:10] == "CACTCAGCTC"
        assert record.contigs[0].reads[4].rd.sequence[center - 5 : center + 5] == "TCCAAAGGGT"
        assert record.contigs[0].reads[4].rd.sequence[-10:] == "AGCTGAATCG"
        assert record.contigs[0].reads[4].qa.qual_clipping_start == 1
        assert record.contigs[0].reads[4].qa.qual_clipping_end == 799
        assert record.contigs[0].reads[4].qa.align_clipping_start == 1
        assert record.contigs[0].reads[4].qa.align_clipping_end == 799
        assert record.contigs[0].reads[4].ds.chromat_file == ""
        assert record.contigs[0].reads[4].ds.phd_file == ""
        assert record.contigs[0].reads[4].ds.time == ""
        assert record.contigs[0].reads[4].ds.chem == ""
        assert record.contigs[0].reads[4].ds.dye == ""
        assert record.contigs[0].reads[4].ds.template == ""
        assert record.contigs[0].reads[4].ds.direction == ""
        assert record.contigs[0].reads[4].rt is None
        assert record.contigs[0].reads[4].wr is None

        assert record.contigs[0].reads[5].rd.name == "R6"
        assert record.contigs[0].reads[5].rd.padded_bases == 857
        assert record.contigs[0].reads[5].rd.info_items == 0
        assert record.contigs[0].reads[5].rd.read_tags == 0
        center = len(record.contigs[0].reads[5].rd.sequence) // 2
        assert record.contigs[0].reads[5].rd.sequence[:10] == "CCGGCAGTGA"
        assert record.contigs[0].reads[5].rd.sequence[center - 5 : center + 5] == "AAAAAAAACC"
        assert record.contigs[0].reads[5].rd.sequence[-10:] == "NNNNNNNNNN"
        assert record.contigs[0].reads[5].qa.qual_clipping_start == 24
        assert record.contigs[0].reads[5].qa.qual_clipping_end == 706
        assert record.contigs[0].reads[5].qa.align_clipping_start == 24
        assert record.contigs[0].reads[5].qa.align_clipping_end == 706
        assert record.contigs[0].reads[5].ds.chromat_file == ""
        assert record.contigs[0].reads[5].ds.phd_file == ""
        assert record.contigs[0].reads[5].ds.time == ""
        assert record.contigs[0].reads[5].ds.chem == ""
        assert record.contigs[0].reads[5].ds.dye == ""
        assert record.contigs[0].reads[5].ds.template == ""
        assert record.contigs[0].reads[5].ds.direction == ""
        assert record.contigs[0].reads[5].rt is None
        assert record.contigs[0].reads[5].wr is None

    def test_check_record_parser(self):
        """Test to check that record parser parses each contig into a record."""
        contigs = Ace.parse(self.handle)

        # First (and only) contig
        contig = next(contigs)

        assert len(contig.reads) == 6
        assert contig.name == "Contig1"
        assert contig.nbases == 1222
        assert contig.nreads == 6
        assert contig.nsegments == 0
        assert contig.uorc == "U"
        center = len(contig.sequence) // 2
        assert contig.sequence[:10] == "AGTTTTAGTT"
        assert contig.sequence[center - 5 : center + 5] == "TGTGCGCGCA"
        assert contig.sequence[-10:] == "ATATCACATT"
        center = len(contig.quality) // 2
        assert contig.quality[:10] == [61, 66, 67, 70, 71, 73, 73, 77, 77, 87]
        assert contig.quality[center - 5 : center + 5] == [97, 97, 97, 97, 97, 97, 97, 97, 97, 97]
        assert contig.quality[-10:] == [56, 51, 49, 41, 38, 39, 45, 44, 49, 46]
        assert len(contig.af) == 6
        assert len(contig.bs) == 0
        assert contig.af[3].name == "R5"
        assert contig.af[3].coru == "C"
        assert contig.af[3].padded_start == 320
        assert contig.af[5].name == "R6"
        assert contig.af[5].coru == "C"
        assert contig.af[5].padded_start == 517
        assert contig.bs == []
        assert contig.ct is None
        assert contig.wa is None
        assert len(contig.reads) == 6

        assert contig.reads[0].rd.name == "R3"
        assert contig.reads[0].rd.padded_bases == 919
        assert contig.reads[0].rd.info_items == 0
        assert contig.reads[0].rd.read_tags == 0
        center = len(contig.reads[0].rd.sequence) // 2
        assert contig.reads[0].rd.sequence[:10] == "NNNNNNNNNN"
        assert contig.reads[0].rd.sequence[center - 5 : center + 5] == "ATGTGCGCTC"
        assert contig.reads[0].rd.sequence[-10:] == "CAGCTCACCA"
        assert contig.reads[0].qa.qual_clipping_start == 55
        assert contig.reads[0].qa.qual_clipping_end == 916
        assert contig.reads[0].qa.align_clipping_start == 55
        assert contig.reads[0].qa.align_clipping_end == 916
        assert contig.reads[0].ds.chromat_file == ""
        assert contig.reads[0].ds.phd_file == ""
        assert contig.reads[0].ds.time == ""
        assert contig.reads[0].ds.chem == ""
        assert contig.reads[0].ds.dye == ""
        assert contig.reads[0].ds.template == ""
        assert contig.reads[0].ds.direction == ""
        assert contig.reads[0].rt is None
        assert contig.reads[0].wr is None

        assert contig.reads[1].rd.name == "R1"
        assert contig.reads[1].rd.padded_bases == 864
        assert contig.reads[1].rd.info_items == 0
        assert contig.reads[1].rd.read_tags == 0
        center = len(contig.reads[1].rd.sequence) // 2
        assert contig.reads[1].rd.sequence[:10] == "AGCCGGTACC"
        assert contig.reads[1].rd.sequence[center - 5 : center + 5] == "GGGATGGCAC"
        assert contig.reads[1].rd.sequence[-10:] == "GGGCTGGGAG"
        assert contig.reads[1].qa.qual_clipping_start == 12
        assert contig.reads[1].qa.qual_clipping_end == 863
        assert contig.reads[1].qa.align_clipping_start == 12
        assert contig.reads[1].qa.align_clipping_end == 863
        assert contig.reads[1].ds.chromat_file == ""
        assert contig.reads[1].ds.phd_file == ""
        assert contig.reads[1].ds.time == ""
        assert contig.reads[1].ds.chem == ""
        assert contig.reads[1].ds.dye == ""
        assert contig.reads[1].ds.template == ""
        assert contig.reads[1].ds.direction == ""
        assert contig.reads[1].rt is None
        assert contig.reads[1].wr is None

        assert contig.reads[2].rd.name == "R2"
        assert contig.reads[2].rd.padded_bases == 1026
        assert contig.reads[2].rd.info_items == 0
        assert contig.reads[2].rd.read_tags == 0
        center = len(contig.reads[2].rd.sequence) // 2
        assert contig.reads[2].rd.sequence[:10] == "NNNNNNNNNN"
        assert contig.reads[2].rd.sequence[center - 5 : center + 5] == "GGATGCCTGG"
        assert contig.reads[2].rd.sequence[-10:] == "GGTTGAGGCC"
        assert contig.reads[2].qa.qual_clipping_start == 55
        assert contig.reads[2].qa.qual_clipping_end == 1000
        assert contig.reads[2].qa.align_clipping_start == 55
        assert contig.reads[2].qa.align_clipping_end == 1000
        assert contig.reads[2].ds.chromat_file == ""
        assert contig.reads[2].ds.phd_file == ""
        assert contig.reads[2].ds.time == ""
        assert contig.reads[2].ds.chem == ""
        assert contig.reads[2].ds.dye == ""
        assert contig.reads[2].ds.template == ""
        assert contig.reads[2].ds.direction == ""
        assert contig.reads[2].rt is None
        assert contig.reads[2].wr is None

        assert contig.reads[3].rd.name == "R5"
        assert contig.reads[3].rd.padded_bases == 925
        assert contig.reads[3].rd.info_items == 0
        assert contig.reads[3].rd.read_tags == 0
        center = len(contig.reads[3].rd.sequence) // 2
        assert contig.reads[3].rd.sequence[:10] == "NNNNNNNNNN"
        assert contig.reads[3].rd.sequence[center - 5 : center + 5] == "CCTCCCTACA"
        assert contig.reads[3].rd.sequence[-10:] == "GCCCCCGGNN"
        assert contig.reads[3].qa.qual_clipping_start == 293
        assert contig.reads[3].qa.qual_clipping_end == 874
        assert contig.reads[3].qa.align_clipping_start == 293
        assert contig.reads[3].qa.align_clipping_end == 874
        assert contig.reads[3].ds.chromat_file == ""
        assert contig.reads[3].ds.phd_file == ""
        assert contig.reads[3].ds.time == ""
        assert contig.reads[3].ds.chem == ""
        assert contig.reads[3].ds.dye == ""
        assert contig.reads[3].ds.template == ""
        assert contig.reads[3].ds.direction == ""
        assert contig.reads[3].rt is None
        assert contig.reads[3].wr is None

        assert contig.reads[4].rd.name == "R4"
        assert contig.reads[4].rd.padded_bases == 816
        assert contig.reads[4].rd.info_items == 0
        assert contig.reads[4].rd.read_tags == 0
        center = len(contig.reads[4].rd.sequence) // 2
        assert contig.reads[4].rd.sequence[:10] == "CACTCAGCTC"
        assert contig.reads[4].rd.sequence[center - 5 : center + 5] == "TCCAAAGGGT"
        assert contig.reads[4].rd.sequence[-10:] == "AGCTGAATCG"
        assert contig.reads[4].qa.qual_clipping_start == 1
        assert contig.reads[4].qa.qual_clipping_end == 799
        assert contig.reads[4].qa.align_clipping_start == 1
        assert contig.reads[4].qa.align_clipping_end == 799
        assert contig.reads[4].ds.chromat_file == ""
        assert contig.reads[4].ds.phd_file == ""
        assert contig.reads[4].ds.time == ""
        assert contig.reads[4].ds.chem == ""
        assert contig.reads[4].ds.dye == ""
        assert contig.reads[4].ds.template == ""
        assert contig.reads[4].ds.direction == ""
        assert contig.reads[4].rt is None
        assert contig.reads[4].wr is None

        assert contig.reads[5].rd.name == "R6"
        assert contig.reads[5].rd.padded_bases == 857
        assert contig.reads[5].rd.info_items == 0
        assert contig.reads[5].rd.read_tags == 0
        center = len(contig.reads[5].rd.sequence) // 2
        assert contig.reads[5].rd.sequence[:10] == "CCGGCAGTGA"
        assert contig.reads[5].rd.sequence[center - 5 : center + 5] == "AAAAAAAACC"
        assert contig.reads[5].rd.sequence[-10:] == "NNNNNNNNNN"
        assert contig.reads[5].qa.qual_clipping_start == 24
        assert contig.reads[5].qa.qual_clipping_end == 706
        assert contig.reads[5].qa.align_clipping_start == 24
        assert contig.reads[5].qa.align_clipping_end == 706
        assert contig.reads[5].ds.chromat_file == ""
        assert contig.reads[5].ds.phd_file == ""
        assert contig.reads[5].ds.time == ""
        assert contig.reads[5].ds.chem == ""
        assert contig.reads[5].ds.dye == ""
        assert contig.reads[5].ds.template == ""
        assert contig.reads[5].ds.direction == ""
        assert contig.reads[5].rt is None
        assert contig.reads[5].wr is None

        # Make sure there are no more contigs
        with pytest.raises(StopIteration):
            next(contigs)


class TestAceTestThree:
    """Test parsing example ACE input file for CONSED.

    The sample input file was downloaded from:
    http://bozeman.mbt.washington.edu/consed/distributions/README.16.0.txt
    """

    @pytest.fixture(autouse=True)
    def _setup(self):
        self.handle = open("Ace/consed_sample.ace")

        yield
        self.handle.close()

    def test_check_ACEParser(self):
        """Test to check that ACEParser can parse the whole file into one record."""
        record = Ace.read(self.handle)
        assert record.ncontigs == 1
        assert record.nreads == 8
        assert len(record.wa) == 1
        assert record.wa[0].tag_type == "phrap_params"
        assert record.wa[0].program == "phrap"
        assert record.wa[0].date == "990621:161947"
        assert record.wa[0].info == [
                "/usr/local/genome/bin/phrap standard.fasta.screen -new_ace -view",
                "phrap version 0.990319",
            ]
        assert len(record.contigs) == 1

        assert len(record.contigs[0].reads) == 8
        assert record.contigs[0].name == "Contig1"
        assert record.contigs[0].nbases == 1475
        assert record.contigs[0].nreads == 8
        assert record.contigs[0].nsegments == 156
        assert record.contigs[0].uorc == "U"
        center = len(record.contigs[0].sequence) // 2
        assert record.contigs[0].sequence[:10] == "agccccgggc"
        assert record.contigs[0].sequence[center - 5 : center + 5] == "CTTCCCCAGG"
        assert record.contigs[0].sequence[-10:] == "gttgggtttg"
        center = len(record.contigs[0].quality) // 2
        assert record.contigs[0].quality[:10] == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        assert record.contigs[0].quality[center - 5 : center + 5] == [90, 90, 90, 90, 90, 90, 90, 90, 89, 89]
        assert record.contigs[0].quality[-10:] == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        assert len(record.contigs[0].af) == 8
        assert len(record.contigs[0].bs) == 156
        assert record.contigs[0].af[4].name == "K26-291s"
        assert record.contigs[0].af[4].coru == "U"
        assert record.contigs[0].af[4].padded_start == 828
        assert record.contigs[0].af[7].name == "K26-766c"
        assert record.contigs[0].af[7].coru == "C"
        assert record.contigs[0].af[7].padded_start == 408
        assert record.contigs[0].bs[78].name == "K26-394c"
        assert record.contigs[0].bs[78].padded_start == 987
        assert record.contigs[0].bs[78].padded_end == 987
        assert record.contigs[0].bs[155].name == "K26-822c"
        assert record.contigs[0].bs[155].padded_start == 1303
        assert record.contigs[0].bs[155].padded_end == 1475
        assert len(record.contigs[0].ct) == 3
        assert record.contigs[0].ct[0].name == "Contig1"
        assert record.contigs[0].ct[0].tag_type == "repeat"
        assert record.contigs[0].ct[0].program == "consed"
        assert record.contigs[0].ct[0].padded_start == 976
        assert record.contigs[0].ct[0].padded_end == 986
        assert record.contigs[0].ct[0].date == "971218:180623"
        assert record.contigs[0].ct[0].info == []
        assert record.contigs[0].ct[1].name == "Contig1"
        assert record.contigs[0].ct[1].tag_type == "comment"
        assert record.contigs[0].ct[1].program == "consed"
        assert record.contigs[0].ct[1].padded_start == 996
        assert record.contigs[0].ct[1].padded_end == 1007
        assert record.contigs[0].ct[1].date == "971218:180623"
        assert record.contigs[0].ct[1].info == ["This is line 1 of a comment", "There may be any number of lines"]
        assert record.contigs[0].ct[2].name == "Contig1"
        assert record.contigs[0].ct[2].tag_type == "oligo"
        assert record.contigs[0].ct[2].program == "consed"
        assert record.contigs[0].ct[2].padded_start == 963
        assert record.contigs[0].ct[2].padded_end == 987
        assert record.contigs[0].ct[2].date == "971218:180623"
        assert record.contigs[0].ct[2].info == ["standard.1 acataagacattctaaatttttact 50 U", "seq from clone"]
        assert len(record.contigs[0].wa) == 1
        assert record.contigs[0].wa[0].tag_type == "phrap_params"
        assert record.contigs[0].wa[0].program == "phrap"
        assert record.contigs[0].wa[0].date == "990621:161947"
        assert record.contigs[0].wa[0].info == [
                "/usr/local/genome/bin/phrap standard.fasta.screen -new_ace -view",
                "phrap version 0.990319",
            ]

        assert len(record.contigs[0].reads) == 8

        assert record.contigs[0].reads[0].rd.name == "K26-217c"
        assert record.contigs[0].reads[0].rd.padded_bases == 563
        assert record.contigs[0].reads[0].rd.info_items == 0
        assert record.contigs[0].reads[0].rd.read_tags == 0
        center = len(record.contigs[0].reads[0].rd.sequence) // 2
        assert record.contigs[0].reads[0].rd.sequence[:10] == "tcccCgtgag"
        assert record.contigs[0].reads[0].rd.sequence[center - 5 : center + 5] == "CTCCTGcctg"
        assert record.contigs[0].reads[0].rd.sequence[-10:] == "ggcccccctc"
        assert record.contigs[0].reads[0].qa.qual_clipping_start == 19
        assert record.contigs[0].reads[0].qa.qual_clipping_end == 349
        assert record.contigs[0].reads[0].qa.align_clipping_start == 19
        assert record.contigs[0].reads[0].qa.align_clipping_end == 424
        assert record.contigs[0].reads[0].ds.chromat_file == "K26-217c"
        assert record.contigs[0].reads[0].ds.phd_file == "K26-217c.phd.1"
        assert record.contigs[0].reads[0].ds.time == "Thu Sep 12 15:42:38 1996"
        assert record.contigs[0].reads[0].ds.chem == ""
        assert record.contigs[0].reads[0].ds.dye == ""
        assert record.contigs[0].reads[0].ds.template == ""
        assert record.contigs[0].reads[0].ds.direction == ""
        assert record.contigs[0].reads[0].rt is None
        assert record.contigs[0].reads[0].wr is None

        assert record.contigs[0].reads[1].rd.name == "K26-526t"
        assert record.contigs[0].reads[1].rd.padded_bases == 687
        assert record.contigs[0].reads[1].rd.info_items == 0
        assert record.contigs[0].reads[1].rd.read_tags == 0
        center = len(record.contigs[0].reads[1].rd.sequence) // 2
        assert record.contigs[0].reads[1].rd.sequence[:10] == "ccgtcctgag"
        assert record.contigs[0].reads[1].rd.sequence[center - 5 : center + 5] == "cacagcccT*"
        assert record.contigs[0].reads[1].rd.sequence[-10:] == "Ttttgtttta"
        assert record.contigs[0].reads[1].qa.qual_clipping_start == 12
        assert record.contigs[0].reads[1].qa.qual_clipping_end == 353
        assert record.contigs[0].reads[1].qa.align_clipping_start == 9
        assert record.contigs[0].reads[1].qa.align_clipping_end == 572
        assert record.contigs[0].reads[1].ds.chromat_file == "K26-526t"
        assert record.contigs[0].reads[1].ds.phd_file == "K26-526t.phd.1"
        assert record.contigs[0].reads[1].ds.time == "Thu Sep 12 15:42:33 1996"
        assert record.contigs[0].reads[1].ds.chem == ""
        assert record.contigs[0].reads[1].ds.dye == ""
        assert record.contigs[0].reads[1].ds.template == ""
        assert record.contigs[0].reads[1].ds.direction == ""
        assert record.contigs[0].reads[1].rt is None
        assert record.contigs[0].reads[1].wr is None

        assert record.contigs[0].reads[2].rd.name == "K26-961c"
        assert record.contigs[0].reads[2].rd.padded_bases == 517
        assert record.contigs[0].reads[2].rd.info_items == 0
        assert record.contigs[0].reads[2].rd.read_tags == 0
        center = len(record.contigs[0].reads[2].rd.sequence) // 2
        assert record.contigs[0].reads[2].rd.sequence[:10] == "aatattaccg"
        assert record.contigs[0].reads[2].rd.sequence[center - 5 : center + 5] == "CAGATGGGTT"
        assert record.contigs[0].reads[2].rd.sequence[-10:] == "ctattcaggg"
        assert record.contigs[0].reads[2].qa.qual_clipping_start == 20
        assert record.contigs[0].reads[2].qa.qual_clipping_end == 415
        assert record.contigs[0].reads[2].qa.align_clipping_start == 26
        assert record.contigs[0].reads[2].qa.align_clipping_end == 514
        assert record.contigs[0].reads[2].ds.chromat_file == "K26-961c"
        assert record.contigs[0].reads[2].ds.phd_file == "K26-961c.phd.1"
        assert record.contigs[0].reads[2].ds.time == "Thu Sep 12 15:42:37 1996"
        assert record.contigs[0].reads[2].ds.chem == ""
        assert record.contigs[0].reads[2].ds.dye == ""
        assert record.contigs[0].reads[2].ds.template == ""
        assert record.contigs[0].reads[2].ds.direction == ""
        assert record.contigs[0].reads[2].rt is None
        assert record.contigs[0].reads[2].wr is None

        assert record.contigs[0].reads[3].rd.name == "K26-394c"
        assert record.contigs[0].reads[3].rd.padded_bases == 628
        assert record.contigs[0].reads[3].rd.info_items == 0
        assert record.contigs[0].reads[3].rd.read_tags == 0
        center = len(record.contigs[0].reads[3].rd.sequence) // 2
        assert record.contigs[0].reads[3].rd.sequence[:10] == "ctgcgtatcg"
        assert record.contigs[0].reads[3].rd.sequence[center - 5 : center + 5] == "AGGATTGCTT"
        assert record.contigs[0].reads[3].rd.sequence[-10:] == "aaccctgggt"
        assert record.contigs[0].reads[3].qa.qual_clipping_start == 18
        assert record.contigs[0].reads[3].qa.qual_clipping_end == 368
        assert record.contigs[0].reads[3].qa.align_clipping_start == 11
        assert record.contigs[0].reads[3].qa.align_clipping_end == 502
        assert record.contigs[0].reads[3].ds.chromat_file == "K26-394c"
        assert record.contigs[0].reads[3].ds.phd_file == "K26-394c.phd.1"
        assert record.contigs[0].reads[3].ds.time == "Thu Sep 12 15:42:32 1996"
        assert record.contigs[0].reads[3].ds.chem == ""
        assert record.contigs[0].reads[3].ds.dye == ""
        assert record.contigs[0].reads[3].ds.template == ""
        assert record.contigs[0].reads[3].ds.direction == ""
        assert record.contigs[0].reads[3].rt is None
        assert record.contigs[0].reads[3].wr is None

        assert record.contigs[0].reads[4].rd.name == "K26-291s"
        assert record.contigs[0].reads[4].rd.padded_bases == 556
        assert record.contigs[0].reads[4].rd.info_items == 0
        assert record.contigs[0].reads[4].rd.read_tags == 0
        center = len(record.contigs[0].reads[4].rd.sequence) // 2
        assert record.contigs[0].reads[4].rd.sequence[:10] == "gaggatcgct"
        assert record.contigs[0].reads[4].rd.sequence[center - 5 : center + 5] == "GTgcgaggat"
        assert record.contigs[0].reads[4].rd.sequence[-10:] == "caggcagatg"
        assert record.contigs[0].reads[4].qa.qual_clipping_start == 11
        assert record.contigs[0].reads[4].qa.qual_clipping_end == 373
        assert record.contigs[0].reads[4].qa.align_clipping_start == 11
        assert record.contigs[0].reads[4].qa.align_clipping_end == 476
        assert record.contigs[0].reads[4].ds.chromat_file == "K26-291s"
        assert record.contigs[0].reads[4].ds.phd_file == "K26-291s.phd.1"
        assert record.contigs[0].reads[4].ds.time == "Thu Sep 12 15:42:31 1996"
        assert record.contigs[0].reads[4].ds.chem == ""
        assert record.contigs[0].reads[4].ds.dye == ""
        assert record.contigs[0].reads[4].ds.template == ""
        assert record.contigs[0].reads[4].ds.direction == ""
        assert record.contigs[0].reads[4].rt is None
        assert record.contigs[0].reads[4].wr is None

        assert record.contigs[0].reads[5].rd.name == "K26-822c"
        assert record.contigs[0].reads[5].rd.padded_bases == 593
        assert record.contigs[0].reads[5].rd.info_items == 0
        assert record.contigs[0].reads[5].rd.read_tags == 0
        center = len(record.contigs[0].reads[5].rd.sequence) // 2
        assert record.contigs[0].reads[5].rd.sequence[:10] == "ggggatccg*"
        assert record.contigs[0].reads[5].rd.sequence[center - 5 : center + 5] == "GCaAgacCCt"
        assert record.contigs[0].reads[5].rd.sequence[-10:] == "gttgggtttg"

        assert record.contigs[0].reads[5].qa.qual_clipping_start == 25
        assert record.contigs[0].reads[5].qa.qual_clipping_end == 333
        assert record.contigs[0].reads[5].qa.align_clipping_start == 16
        assert record.contigs[0].reads[5].qa.align_clipping_end == 593
        assert record.contigs[0].reads[5].ds.chromat_file == "K26-822c"
        assert record.contigs[0].reads[5].ds.phd_file == "K26-822c.phd.1"
        assert record.contigs[0].reads[5].ds.time == "Thu Sep 12 15:42:36 1996"
        assert record.contigs[0].reads[5].ds.chem == ""
        assert record.contigs[0].reads[5].ds.dye == ""
        assert record.contigs[0].reads[5].ds.template == ""
        assert record.contigs[0].reads[5].ds.direction == ""
        assert record.contigs[0].reads[5].rt is None
        assert record.contigs[0].reads[5].wr is None

        assert record.contigs[0].reads[6].rd.name == "K26-572c"
        assert record.contigs[0].reads[6].rd.padded_bases == 594
        assert record.contigs[0].reads[6].rd.info_items == 0
        assert record.contigs[0].reads[6].rd.read_tags == 0
        center = len(record.contigs[0].reads[6].rd.sequence) // 2
        assert record.contigs[0].reads[6].rd.sequence[:10] == "agccccgggc"
        assert record.contigs[0].reads[6].rd.sequence[center - 5 : center + 5] == "ggatcACATA"
        assert record.contigs[0].reads[6].rd.sequence[-10:] == "aatagtaaca"
        assert record.contigs[0].reads[6].qa.qual_clipping_start == 249
        assert record.contigs[0].reads[6].qa.qual_clipping_end == 584
        assert record.contigs[0].reads[6].qa.align_clipping_start == 1
        assert record.contigs[0].reads[6].qa.align_clipping_end == 586
        assert record.contigs[0].reads[6].ds.chromat_file == "K26-572c"
        assert record.contigs[0].reads[6].ds.phd_file == "K26-572c.phd.1"
        assert record.contigs[0].reads[6].ds.time == "Thu Sep 12 15:42:34 1996"
        assert record.contigs[0].reads[6].ds.chem == ""
        assert record.contigs[0].reads[6].ds.dye == ""
        assert record.contigs[0].reads[6].ds.template == ""
        assert record.contigs[0].reads[6].ds.direction == ""
        assert record.contigs[0].reads[6].rt is None
        assert record.contigs[0].reads[6].wr is None

        assert record.contigs[0].reads[7].rd.name == "K26-766c"
        assert record.contigs[0].reads[7].rd.padded_bases == 603
        assert record.contigs[0].reads[7].rd.info_items == 0
        assert record.contigs[0].reads[7].rd.read_tags == 0
        center = len(record.contigs[0].reads[7].rd.sequence) // 2
        assert record.contigs[0].reads[7].rd.sequence[:10] == "gaataattgg"
        assert record.contigs[0].reads[7].rd.sequence[center - 5 : center + 5] == "TggCCCATCT"
        assert record.contigs[0].reads[7].rd.sequence[-10:] == "gaaccacacg"
        assert record.contigs[0].reads[7].qa.qual_clipping_start == 240
        assert record.contigs[0].reads[7].qa.qual_clipping_end == 584
        assert record.contigs[0].reads[7].qa.align_clipping_start == 126
        assert record.contigs[0].reads[7].qa.align_clipping_end == 583
        assert record.contigs[0].reads[7].ds.chromat_file == "K26-766c"
        assert record.contigs[0].reads[7].ds.phd_file == "K26-766c.phd.1"
        assert record.contigs[0].reads[7].ds.time == "Thu Sep 12 15:42:35 1996"
        assert record.contigs[0].reads[7].ds.chem == ""
        assert record.contigs[0].reads[7].ds.dye == ""
        assert record.contigs[0].reads[7].ds.template == ""
        assert record.contigs[0].reads[7].ds.direction == ""
        assert record.contigs[0].reads[7].rt is None
        assert record.contigs[0].reads[7].wr is None

    def test_check_record_parser(self):
        """Test to check that record parser parses each contig into a record."""
        contigs = Ace.parse(self.handle)

        # First (and only) contig
        contig = next(contigs)

        assert len(contig.reads) == 8
        assert contig.name == "Contig1"
        assert contig.nbases == 1475
        assert contig.nreads == 8
        assert contig.nsegments == 156
        assert contig.uorc == "U"
        center = len(contig.sequence) // 2
        assert contig.sequence[:10] == "agccccgggc"
        assert contig.sequence[center - 5 : center + 5] == "CTTCCCCAGG"
        assert contig.sequence[-10:] == "gttgggtttg"
        center = len(contig.quality) // 2
        assert contig.quality[:10] == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        assert contig.quality[center - 5 : center + 5] == [90, 90, 90, 90, 90, 90, 90, 90, 89, 89]
        assert contig.quality[-10:] == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        assert len(contig.af) == 8
        assert len(contig.bs) == 156
        assert contig.af[4].name == "K26-291s"
        assert contig.af[4].coru == "U"
        assert contig.af[4].padded_start == 828
        assert contig.af[7].name == "K26-766c"
        assert contig.af[7].coru == "C"
        assert contig.af[7].padded_start == 408
        assert contig.bs[78].name == "K26-394c"
        assert contig.bs[78].padded_start == 987
        assert contig.bs[78].padded_end == 987
        assert contig.bs[155].name == "K26-822c"
        assert contig.bs[155].padded_start == 1303
        assert contig.bs[155].padded_end == 1475
        assert len(contig.ct) == 3
        assert contig.ct[0].name == "Contig1"
        assert contig.ct[0].tag_type == "repeat"
        assert contig.ct[0].program == "consed"
        assert contig.ct[0].padded_start == 976
        assert contig.ct[0].padded_end == 986
        assert contig.ct[0].date == "971218:180623"
        assert contig.ct[0].info == []
        assert contig.ct[1].name == "Contig1"
        assert contig.ct[1].tag_type == "comment"
        assert contig.ct[1].program == "consed"
        assert contig.ct[1].padded_start == 996
        assert contig.ct[1].padded_end == 1007
        assert contig.ct[1].date == "971218:180623"
        assert contig.ct[1].info == ["This is line 1 of a comment", "There may be any number of lines"]
        assert contig.ct[2].name == "Contig1"
        assert contig.ct[2].tag_type == "oligo"
        assert contig.ct[2].program == "consed"
        assert contig.ct[2].padded_start == 963
        assert contig.ct[2].padded_end == 987
        assert contig.ct[2].date == "971218:180623"
        assert contig.ct[2].info == ["standard.1 acataagacattctaaatttttact 50 U", "seq from clone"]
        assert len(contig.wa) == 1
        assert contig.wa[0].tag_type == "phrap_params"
        assert contig.wa[0].program == "phrap"
        assert contig.wa[0].date == "990621:161947"
        assert contig.wa[0].info == [
                "/usr/local/genome/bin/phrap standard.fasta.screen -new_ace -view",
                "phrap version 0.990319",
            ]

        assert len(contig.reads) == 8

        assert contig.reads[0].rd.name == "K26-217c"
        assert contig.reads[0].rd.padded_bases == 563
        assert contig.reads[0].rd.info_items == 0
        assert contig.reads[0].rd.read_tags == 0
        center = len(contig.reads[0].rd.sequence) // 2
        assert contig.reads[0].rd.sequence[:10] == "tcccCgtgag"
        assert contig.reads[0].rd.sequence[center - 5 : center + 5] == "CTCCTGcctg"
        assert contig.reads[0].rd.sequence[-10:] == "ggcccccctc"
        assert contig.reads[0].qa.qual_clipping_start == 19
        assert contig.reads[0].qa.qual_clipping_end == 349
        assert contig.reads[0].qa.align_clipping_start == 19
        assert contig.reads[0].qa.align_clipping_end == 424
        assert contig.reads[0].ds.chromat_file == "K26-217c"
        assert contig.reads[0].ds.phd_file == "K26-217c.phd.1"
        assert contig.reads[0].ds.time == "Thu Sep 12 15:42:38 1996"
        assert contig.reads[0].ds.chem == ""
        assert contig.reads[0].ds.dye == ""
        assert contig.reads[0].ds.template == ""
        assert contig.reads[0].ds.direction == ""
        assert contig.reads[0].rt is None
        assert contig.reads[0].wr is None

        assert contig.reads[1].rd.name == "K26-526t"
        assert contig.reads[1].rd.padded_bases == 687
        assert contig.reads[1].rd.info_items == 0
        assert contig.reads[1].rd.read_tags == 0
        center = len(contig.reads[1].rd.sequence) // 2
        assert contig.reads[1].rd.sequence[:10] == "ccgtcctgag"
        assert contig.reads[1].rd.sequence[center - 5 : center + 5] == "cacagcccT*"
        assert contig.reads[1].rd.sequence[-10:] == "Ttttgtttta"
        assert contig.reads[1].qa.qual_clipping_start == 12
        assert contig.reads[1].qa.qual_clipping_end == 353
        assert contig.reads[1].qa.align_clipping_start == 9
        assert contig.reads[1].qa.align_clipping_end == 572
        assert contig.reads[1].ds.chromat_file == "K26-526t"
        assert contig.reads[1].ds.phd_file == "K26-526t.phd.1"
        assert contig.reads[1].ds.time == "Thu Sep 12 15:42:33 1996"
        assert contig.reads[1].ds.chem == ""
        assert contig.reads[1].ds.dye == ""
        assert contig.reads[1].ds.template == ""
        assert contig.reads[1].ds.direction == ""
        assert contig.reads[1].rt is None
        assert contig.reads[1].wr is None

        assert contig.reads[2].rd.name == "K26-961c"
        assert contig.reads[2].rd.padded_bases == 517
        assert contig.reads[2].rd.info_items == 0
        assert contig.reads[2].rd.read_tags == 0
        center = len(contig.reads[2].rd.sequence) // 2
        assert contig.reads[2].rd.sequence[:10] == "aatattaccg"
        assert contig.reads[2].rd.sequence[center - 5 : center + 5] == "CAGATGGGTT"
        assert contig.reads[2].rd.sequence[-10:] == "ctattcaggg"
        assert contig.reads[2].qa.qual_clipping_start == 20
        assert contig.reads[2].qa.qual_clipping_end == 415
        assert contig.reads[2].qa.align_clipping_start == 26
        assert contig.reads[2].qa.align_clipping_end == 514
        assert contig.reads[2].ds.chromat_file == "K26-961c"
        assert contig.reads[2].ds.phd_file == "K26-961c.phd.1"
        assert contig.reads[2].ds.time == "Thu Sep 12 15:42:37 1996"
        assert contig.reads[2].ds.chem == ""
        assert contig.reads[2].ds.dye == ""
        assert contig.reads[2].ds.template == ""
        assert contig.reads[2].ds.direction == ""
        assert contig.reads[2].rt is None
        assert contig.reads[2].wr is None

        assert contig.reads[3].rd.name == "K26-394c"
        assert contig.reads[3].rd.padded_bases == 628
        assert contig.reads[3].rd.info_items == 0
        assert contig.reads[3].rd.read_tags == 0
        center = len(contig.reads[3].rd.sequence) // 2
        assert contig.reads[3].rd.sequence[:10] == "ctgcgtatcg"
        assert contig.reads[3].rd.sequence[center - 5 : center + 5] == "AGGATTGCTT"
        assert contig.reads[3].rd.sequence[-10:] == "aaccctgggt"
        assert contig.reads[3].qa.qual_clipping_start == 18
        assert contig.reads[3].qa.qual_clipping_end == 368
        assert contig.reads[3].qa.align_clipping_start == 11
        assert contig.reads[3].qa.align_clipping_end == 502
        assert contig.reads[3].ds.chromat_file == "K26-394c"
        assert contig.reads[3].ds.phd_file == "K26-394c.phd.1"
        assert contig.reads[3].ds.time == "Thu Sep 12 15:42:32 1996"
        assert contig.reads[3].ds.chem == ""
        assert contig.reads[3].ds.dye == ""
        assert contig.reads[3].ds.template == ""
        assert contig.reads[3].ds.direction == ""
        assert contig.reads[3].rt is None
        assert contig.reads[3].wr is None

        assert contig.reads[4].rd.name == "K26-291s"
        assert contig.reads[4].rd.padded_bases == 556
        assert contig.reads[4].rd.info_items == 0
        assert contig.reads[4].rd.read_tags == 0
        center = len(contig.reads[4].rd.sequence) // 2
        assert contig.reads[4].rd.sequence[:10] == "gaggatcgct"
        assert contig.reads[4].rd.sequence[center - 5 : center + 5] == "GTgcgaggat"
        assert contig.reads[4].rd.sequence[-10:] == "caggcagatg"
        assert contig.reads[4].qa.qual_clipping_start == 11
        assert contig.reads[4].qa.qual_clipping_end == 373
        assert contig.reads[4].qa.align_clipping_start == 11
        assert contig.reads[4].qa.align_clipping_end == 476
        assert contig.reads[4].ds.chromat_file == "K26-291s"
        assert contig.reads[4].ds.phd_file == "K26-291s.phd.1"
        assert contig.reads[4].ds.time == "Thu Sep 12 15:42:31 1996"
        assert contig.reads[4].ds.chem == ""
        assert contig.reads[4].ds.dye == ""
        assert contig.reads[4].ds.template == ""
        assert contig.reads[4].ds.direction == ""
        assert contig.reads[4].rt is None
        assert contig.reads[4].wr is None

        assert contig.reads[5].rd.name == "K26-822c"
        assert contig.reads[5].rd.padded_bases == 593
        assert contig.reads[5].rd.info_items == 0
        assert contig.reads[5].rd.read_tags == 0
        center = len(contig.reads[5].rd.sequence) // 2
        assert contig.reads[5].rd.sequence[:10] == "ggggatccg*"
        assert contig.reads[5].rd.sequence[center - 5 : center + 5] == "GCaAgacCCt"
        assert contig.reads[5].rd.sequence[-10:] == "gttgggtttg"

        assert contig.reads[5].qa.qual_clipping_start == 25
        assert contig.reads[5].qa.qual_clipping_end == 333
        assert contig.reads[5].qa.align_clipping_start == 16
        assert contig.reads[5].qa.align_clipping_end == 593
        assert contig.reads[5].ds.chromat_file == "K26-822c"
        assert contig.reads[5].ds.phd_file == "K26-822c.phd.1"
        assert contig.reads[5].ds.time == "Thu Sep 12 15:42:36 1996"
        assert contig.reads[5].ds.chem == ""
        assert contig.reads[5].ds.dye == ""
        assert contig.reads[5].ds.template == ""
        assert contig.reads[5].ds.direction == ""
        assert contig.reads[5].rt is None
        assert contig.reads[5].wr is None

        assert contig.reads[6].rd.name == "K26-572c"
        assert contig.reads[6].rd.padded_bases == 594
        assert contig.reads[6].rd.info_items == 0
        assert contig.reads[6].rd.read_tags == 0
        center = len(contig.reads[6].rd.sequence) // 2
        assert contig.reads[6].rd.sequence[:10] == "agccccgggc"
        assert contig.reads[6].rd.sequence[center - 5 : center + 5] == "ggatcACATA"
        assert contig.reads[6].rd.sequence[-10:] == "aatagtaaca"
        assert contig.reads[6].qa.qual_clipping_start == 249
        assert contig.reads[6].qa.qual_clipping_end == 584
        assert contig.reads[6].qa.align_clipping_start == 1
        assert contig.reads[6].qa.align_clipping_end == 586
        assert contig.reads[6].ds.chromat_file == "K26-572c"
        assert contig.reads[6].ds.phd_file == "K26-572c.phd.1"
        assert contig.reads[6].ds.time == "Thu Sep 12 15:42:34 1996"
        assert contig.reads[6].ds.chem == ""
        assert contig.reads[6].ds.dye == ""
        assert contig.reads[6].ds.template == ""
        assert contig.reads[6].ds.direction == ""
        assert contig.reads[6].rt is None
        assert contig.reads[6].wr is None

        assert contig.reads[7].rd.name == "K26-766c"
        assert contig.reads[7].rd.padded_bases == 603
        assert contig.reads[7].rd.info_items == 0
        assert contig.reads[7].rd.read_tags == 0
        center = len(contig.reads[7].rd.sequence) // 2
        assert contig.reads[7].rd.sequence[:10] == "gaataattgg"
        assert contig.reads[7].rd.sequence[center - 5 : center + 5] == "TggCCCATCT"
        assert contig.reads[7].rd.sequence[-10:] == "gaaccacacg"
        assert contig.reads[7].qa.qual_clipping_start == 240
        assert contig.reads[7].qa.qual_clipping_end == 584
        assert contig.reads[7].qa.align_clipping_start == 126
        assert contig.reads[7].qa.align_clipping_end == 583
        assert contig.reads[7].ds.chromat_file == "K26-766c"
        assert contig.reads[7].ds.phd_file == "K26-766c.phd.1"
        assert contig.reads[7].ds.time == "Thu Sep 12 15:42:35 1996"
        assert contig.reads[7].ds.chem == ""
        assert contig.reads[7].ds.dye == ""
        assert contig.reads[7].ds.template == ""
        assert contig.reads[7].ds.direction == ""
        assert contig.reads[7].rt is None
        assert contig.reads[7].wr is None

        # Make sure there are no more contigs
        with pytest.raises(StopIteration):
            next(contigs)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
