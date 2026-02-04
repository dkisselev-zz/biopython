# Copyright 2008 Michiel de Hoon
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Tests for Medline module."""

import unittest
import pytest

from Bio import Medline


class TestMedline(unittest.TestCase):
    def test_read(self):
        with open("Medline/pubmed_result1.txt") as handle:
            record = Medline.read(handle)
        assert record["PMID"] == "12230038"
        assert record["OWN"] == "NLM"
        assert record["STAT"] == "MEDLINE"
        assert record["DA"] == "20020916"
        assert record["DCOM"] == "20030606"
        assert record["LR"] == "20041117"
        assert record["PUBM"] == "Print"
        assert record["IS"] == "1467-5463 (Print)"
        assert record["VI"] == "3"
        assert record["IP"] == "3"
        assert record["DP"] == "2002 Sep"
        assert record["TI"] == "The Bio* toolkits--a brief overview."
        assert record["PG"] == "296-302"
        assert record["AB"] == "Bioinformatics research is often difficult to do with commercial software. The Open Source BioPerl, BioPython and Biojava projects provide toolkits with multiple functionality that make it easier to create customised pipelines or analysis. This review briefly compares the quirks of the underlying languages and the functionality, documentation, utility and relative advantages of the Bio counterparts, particularly from the point of view of the beginning biologist programmer."
        assert record["AD"] == ["tacg Informatics, Irvine, CA 92612, USA. hjm@tacgi.com"]
        assert record["FAU"] == ["Mangalam, Harry"]
        assert record["AU"] == ["Mangalam H"]
        assert record["LA"] == ["eng"]
        assert record["PT"] == ["Journal Article"]
        assert record["PL"] == "England"
        assert record["TA"] == "Brief Bioinform"
        assert record["JT"] == "Briefings in bioinformatics"
        assert record["JID"] == "100912837"
        assert record["SB"] == "IM"
        assert record["MH"] == [
                "*Computational Biology",
                "Computer Systems",
                "Humans",
                "Internet",
                "*Programming Languages",
                "*Software",
                "User-Computer Interface",
            ]
        assert record["EDAT"] == "2002/09/17 10:00"
        assert record["MHDA"] == "2003/06/07 05:00"
        assert record["PST"] == "ppublish"
        assert record["SO"] == "Brief Bioinform. 2002 Sep;3(3):296-302."

    def test_parse(self):
        with open("Medline/pubmed_result2.txt") as handle:
            records = Medline.parse(handle)
            record = next(records)
            assert record["PMID"] == "16403221"
            assert record["OWN"] == "NLM"
            assert record["STAT"] == "MEDLINE"
            assert record["DA"] == "20060220"
            assert record["DCOM"] == "20060314"
            assert record["PUBM"] == "Electronic"
            assert record["IS"] == "1471-2105 (Electronic)"
            assert record["VI"] == "7"
            assert record["DP"] == "2006"
            assert record["TI"] == "A high level interface to SCOP and ASTRAL implemented in python."
            assert record["PG"] == "10"
            assert record["AB"] == "BACKGROUND: Benchmarking algorithms in structural bioinformatics often involves the construction of datasets of proteins with given sequence and structural properties. The SCOP database is a manually curated structural classification which groups together proteins on the basis of structural similarity. The ASTRAL compendium provides non redundant subsets of SCOP domains on the basis of sequence similarity such that no two domains in a given subset share more than a defined degree of sequence similarity. Taken together these two resources provide a 'ground truth' for assessing structural bioinformatics algorithms. We present a small and easy to use API written in python to enable construction of datasets from these resources. RESULTS: We have designed a set of python modules to provide an abstraction of the SCOP and ASTRAL databases. The modules are designed to work as part of the Biopython distribution. Python users can now manipulate and use the SCOP hierarchy from within python programs, and use ASTRAL to return sequences of domains in SCOP, as well as clustered representations of SCOP from ASTRAL. CONCLUSION: The modules make the analysis and generation of datasets for use in structural genomics easier and more principled."
            assert record["AD"] == [
                    "Bioinformatics, Institute of Cell and Molecular Science, School of Medicine and Dentistry, Queen Mary, University of London, London EC1 6BQ, UK. j.a.casbon@qmul.ac.uk"
                ]
            assert record["FAU"] == ["Casbon, James A", "Crooks, Gavin E", "Saqi, Mansoor A S"]
            assert record["AU"] == ["Casbon JA", "Crooks GE", "Saqi MA"]
            assert record["LA"] == ["eng"]
            assert record["PT"] == ["Evaluation Studies", "Journal Article"]
            assert record["DEP"] == "20060110"
            assert record["PL"] == "England"
            assert record["TA"] == "BMC Bioinformatics"
            assert record["JT"] == "BMC bioinformatics"
            assert record["JID"] == "100965194"
            assert record["SB"] == "IM"
            assert record["MH"] == [
                    "*Database Management Systems",
                    "*Databases, Protein",
                    "Information Storage and Retrieval/*methods",
                    "Programming Languages",
                    "Sequence Alignment/*methods",
                    "Sequence Analysis, Protein/*methods",
                    "Sequence Homology, Amino Acid",
                    "*Software",
                    "*User-Computer Interface",
                ]
            assert record["PMC"] == "PMC1373603"
            assert record["EDAT"] == "2006/01/13 09:00"
            assert record["MHDA"] == "2006/03/15 09:00"
            assert record["PHST"] == [
                    "2005/06/17 [received]",
                    "2006/01/10 [accepted]",
                    "2006/01/10 [aheadofprint]",
                ]
            assert record["AID"] == ["1471-2105-7-10 [pii]", "10.1186/1471-2105-7-10 [doi]"]
            assert record["PST"] == "epublish"
            assert record["SO"] == "BMC Bioinformatics. 2006 Jan 10;7:10."
            record = next(records)
            assert record["PMID"] == "16377612"
            assert record["OWN"] == "NLM"
            assert record["STAT"] == "MEDLINE"
            assert record["DA"] == "20060223"
            assert record["DCOM"] == "20060418"
            assert record["LR"] == "20061115"
            assert record["PUBM"] == "Print-Electronic"
            assert record["IS"] == "1367-4803 (Print)"
            assert record["VI"] == "22"
            assert record["IP"] == "5"
            assert record["DP"] == "2006 Mar 1"
            assert record["TI"] == "GenomeDiagram: a python package for the visualization of large-scale genomic data."
            assert record["PG"] == "616-7"
            assert record["AB"] == "SUMMARY: We present GenomeDiagram, a flexible, open-source Python module for the visualization of large-scale genomic, comparative genomic and other data with reference to a single chromosome or other biological sequence. GenomeDiagram may be used to generate publication-quality vector graphics, rastered images and in-line streamed graphics for webpages. The package integrates with datatypes from the BioPython project, and is available for Windows, Linux and Mac OS X systems. AVAILABILITY: GenomeDiagram is freely available as source code (under GNU Public License) at http://bioinf.scri.ac.uk/lp/programs.html, and requires Python 2.3 or higher, and recent versions of the ReportLab and BioPython packages. SUPPLEMENTARY INFORMATION: A user manual, example code and images are available at http://bioinf.scri.ac.uk/lp/programs.html."
            assert record["AD"] == [
                    "Plant Pathogen Programme, Scottish Crop Research Institute, Invergowrie, Dundee DD2 5DA, Scotland, UK. lpritc@scri.ac.uk"
                ]
            assert record["FAU"] == [
                    "Pritchard, Leighton",
                    "White, Jennifer A",
                    "Birch, Paul R J",
                    "Toth, Ian K",
                ]
            assert record["AU"] == ["Pritchard L", "White JA", "Birch PR", "Toth IK"]
            assert record["LA"] == ["eng"]
            assert record["PT"] == ["Journal Article", "Research Support, Non-U.S. Gov't"]
            assert record["DEP"] == "20051223"
            assert record["PL"] == "England"
            assert record["TA"] == "Bioinformatics"
            assert record["JT"] == "Bioinformatics (Oxford, England)"
            assert record["JID"] == "9808944"
            assert record["SB"] == "IM"
            assert record["MH"] == [
                    "Chromosome Mapping/*methods",
                    "*Computer Graphics",
                    "*Database Management Systems",
                    "*Databases, Genetic",
                    "Information Storage and Retrieval/methods",
                    "*Programming Languages",
                    "*Software",
                    "*User-Computer Interface",
                ]
            assert record["EDAT"] == "2005/12/27 09:00"
            assert record["MHDA"] == "2006/04/19 09:00"
            assert record["PHST"] == ["2005/12/23 [aheadofprint]"]
            assert record["AID"] == ["btk021 [pii]", "10.1093/bioinformatics/btk021 [doi]"]
            assert record["PST"] == "ppublish"
            assert record["SO"] == "Bioinformatics. 2006 Mar 1;22(5):616-7. Epub 2005 Dec 23."
            record = next(records)
            assert record["PMID"] == "14871861"
            assert record["OWN"] == "NLM"
            assert record["STAT"] == "MEDLINE"
            assert record["DA"] == "20040611"
            assert record["DCOM"] == "20050104"
            assert record["LR"] == "20061115"
            assert record["PUBM"] == "Print-Electronic"
            assert record["IS"] == "1367-4803 (Print)"
            assert record["VI"] == "20"
            assert record["IP"] == "9"
            assert record["DP"] == "2004 Jun 12"
            assert record["TI"] == "Open source clustering software."
            assert record["PG"] == "1453-4"
            assert record["AB"] == "SUMMARY: We have implemented k-means clustering, hierarchical clustering and self-organizing maps in a single multipurpose open-source library of C routines, callable from other C and C++ programs. Using this library, we have created an improved version of Michael Eisen's well-known Cluster program for Windows, Mac OS X and Linux/Unix. In addition, we generated a Python and a Perl interface to the C Clustering Library, thereby combining the flexibility of a scripting language with the speed of C. AVAILABILITY: The C Clustering Library and the corresponding Python C extension module Pycluster were released under the Python License, while the Perl module Algorithm::Cluster was released under the Artistic License. The GUI code Cluster 3.0 for Windows, Macintosh and Linux/Unix, as well as the corresponding command-line program, were released under the same license as the original Cluster code. The complete source code is available at http://bonsai.ims.u-tokyo.ac.jp/mdehoon/software/cluster. Alternatively, Algorithm::Cluster can be downloaded from CPAN, while Pycluster is also available as part of the Biopython distribution."
            assert record["AD"] == [
                    "Human Genome Center, Institute of Medical Science, University of Tokyo, 4-6-1 Shirokanedai, Minato-ku, Tokyo, 108-8639 Japan. mdehoon@ims.u-tokyo.ac.jp"
                ]
            assert record["FAU"] == ["de Hoon, M J L", "Imoto, S", "Nolan, J", "Miyano, S"]
            assert record["AU"] == ["de Hoon MJ", "Imoto S", "Nolan J", "Miyano S"]
            assert record["LA"] == ["eng"]
            assert record["PT"] == [
                    "Comparative Study",
                    "Evaluation Studies",
                    "Journal Article",
                    "Validation Studies",
                ]
            assert record["DEP"] == "20040210"
            assert record["PL"] == "England"
            assert record["TA"] == "Bioinformatics"
            assert record["JT"] == "Bioinformatics (Oxford, England)"
            assert record["JID"] == "9808944"
            assert record["SB"] == "IM"
            assert record["MH"] == [
                    "*Algorithms",
                    "*Cluster Analysis",
                    "Gene Expression Profiling/*methods",
                    "Pattern Recognition, Automated/methods",
                    "*Programming Languages",
                    "Sequence Alignment/*methods",
                    "Sequence Analysis, DNA/*methods",
                    "*Software",
                ]
            assert record["EDAT"] == "2004/02/12 05:00"
            assert record["MHDA"] == "2005/01/05 09:00"
            assert record["PHST"] == ["2004/02/10 [aheadofprint]"]
            assert record["AID"] == ["10.1093/bioinformatics/bth078 [doi]", "bth078 [pii]"]
            assert record["PST"] == "ppublish"
            assert record["SO"] == "Bioinformatics. 2004 Jun 12;20(9):1453-4. Epub 2004 Feb 10."
            record = next(records)
            assert record["PMID"] == "14630660"
            assert record["OWN"] == "NLM"
            assert record["STAT"] == "MEDLINE"
            assert record["DA"] == "20031121"
            assert record["DCOM"] == "20040722"
            assert record["LR"] == "20061115"
            assert record["PUBM"] == "Print"
            assert record["IS"] == "1367-4803 (Print)"
            assert record["VI"] == "19"
            assert record["IP"] == "17"
            assert record["DP"] == "2003 Nov 22"
            assert record["TI"] == "PDB file parser and structure class implemented in Python."
            assert record["PG"] == "2308-10"
            assert record["AB"] == "The biopython project provides a set of bioinformatics tools implemented in Python. Recently, biopython was extended with a set of modules that deal with macromolecular structure. Biopython now contains a parser for PDB files that makes the atomic information available in an easy-to-use but powerful data structure. The parser and data structure deal with features that are often left out or handled inadequately by other packages, e.g. atom and residue disorder (if point mutants are present in the crystal), anisotropic B factors, multiple models and insertion codes. In addition, the parser performs some sanity checking to detect obvious errors. AVAILABILITY: The Biopython distribution (including source code and documentation) is freely available (under the Biopython license) from http://www.biopython.org"
            assert record["AD"] == [
                    "Department of Cellular and Molecular Interactions, Vlaams Interuniversitair Instituut voor Biotechnologie and Computational Modeling Lab, Department of Computer Science, Vrije Universiteit Brussel, Pleinlaan 2, 1050 Brussels, Belgium. thamelry@vub.ac.be"
                ]
            assert record["FAU"] == ["Hamelryck, Thomas", "Manderick, Bernard"]
            assert record["AU"] == ["Hamelryck T", "Manderick B"]
            assert record["LA"] == ["eng"]
            assert record["PT"] == [
                    "Comparative Study",
                    "Evaluation Studies",
                    "Journal Article",
                    "Research Support, Non-U.S. Gov't",
                    "Validation Studies",
                ]
            assert record["PL"] == "England"
            assert record["TA"] == "Bioinformatics"
            assert record["JT"] == "Bioinformatics (Oxford, England)"
            assert record["JID"] == "9808944"
            assert record["RN"] == ["0 (Macromolecular Substances)"]
            assert record["SB"] == "IM"
            assert record["MH"] == [
                    "Computer Simulation",
                    "Database Management Systems/*standards",
                    "*Databases, Protein",
                    "Information Storage and Retrieval/*methods/*standards",
                    "Macromolecular Substances",
                    "*Models, Molecular",
                    "*Programming Languages",
                    "Protein Conformation",
                    "*Software",
                ]
            assert record["EDAT"] == "2003/11/25 05:00"
            assert record["MHDA"] == "2004/07/23 05:00"
            assert record["PST"] == "ppublish"
            assert record["SO"] == "Bioinformatics. 2003 Nov 22;19(17):2308-10."
            with pytest.raises(StopIteration):
                next(records)

    def test_multiline_mesh(self):
        with open("Medline/pubmed_result3.txt") as handle:
            record = Medline.read(handle)
            assert record["PMID"] == "23039619"
        assert record["MH"] == [
                "Blood Circulation",
                "High-Intensity Focused Ultrasound Ablation/adverse effects/instrumentation/*methods",
                "Humans",
                "Models, Biological",
                "Sonication",
                "Temperature",
                "Time Factors",
                "Transducers",
            ]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
