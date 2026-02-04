# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""Tests the basic functionality of the GEO parsers."""

import unittest
import pytest

from Bio import Geo


class TestGeo(unittest.TestCase):
    def test_soft_ex_dual(self):
        path = "Geo/soft_ex_dual.txt"
        with open(path, encoding="latin") as handle:
            records = Geo.parse(handle)
            record = next(records)
            assert record.entity_type == "SAMPLE"
            assert record.entity_id == "Control Embyronic Stem Cell Replicate 1"
            assert len(record.entity_attributes) == 24
            assert record.entity_attributes["Sample_extract_protocol_ch2"] == "TriZol procedure"
            assert record.entity_attributes["Sample_hyb_protocol"] == "Oligoarray control targets and hybridization buffer (Agilent In Situ Hybridization Kit Plus) were added, and samples were applied to microarrays enclosed in Agilent SureHyb-enabled hybridization chambers. After hybridization, slides were washed sequentially with 6x SSC/0.005% Triton X-102 and 0.1x SSC/0.005% Triton X-102 before scanning. Slides were hybridized for 17 h at 60\xb0C in a rotating oven, and washed."
            assert record.entity_attributes["Sample_extract_protocol_ch1"] == "TriZol procedure"
            assert record.entity_attributes["Sample_platform_id"] == "GPL3759"
            assert record.entity_attributes["Sample_title"] == "Control Embyronic Stem Cell Replicate 1"
            assert record.entity_attributes["Sample_supplementary_file"] == "file1.gpr"
            assert record.entity_attributes["Sample_organism_ch2"] == "Mus musculus"
            assert record.entity_attributes["Sample_organism_ch1"] == "Mus musculus"
            assert record.entity_attributes["Sample_label_ch1"] == "Cy5"
            assert len(record.entity_attributes["Sample_scan_protocol"]) == 2
            assert record.entity_attributes["Sample_scan_protocol"][0] == "Scanned on an Agilent G2565AA scanner."
            assert record.entity_attributes["Sample_scan_protocol"][1] == "Images were quantified using Agilent Feature Extraction Software (version A.7.5)."
            assert record.entity_attributes["sample_table_begin"] == ""
            assert record.entity_attributes["Sample_label_protocol_ch2"] == "10 \xb5g of total RNA were primed with 2 \xb5l of 100 \xb5M T16N2 DNA primer at 70\xb0C for 10 min, then reversed transcribed at 42\xb0C for 1 h in the presence of 400 U SuperScript II RTase (Invitrogen), and 100 \xb5M each dATP, dTTP, dGTP, with 25 \xb5M dCTP, 25 \xb5M Cy5-labeled dCTP (NEN Life Science, Boston, MA), and RNase inhibitor (Invitrogen). RNA was then degraded with RNase A, and labeled cDNAs were purified using QIAquick PCR columns (Qiagen)."
            assert record.entity_attributes["Sample_label_protocol_ch1"] == "10 \xb5g of total RNA were primed with 2 \xb5l of 100 \xb5M T16N2 DNA primer at 70\xb0C for 10 min, then reversed transcribed at 42\xb0C for 1 h in the presence of 400 U SuperScript II RTase (Invitrogen), and 100 \xb5M each dATP, dTTP, dGTP, with 25 \xb5M dCTP, 25 \xb5M Cy5-labeled dCTP (NEN Life Science, Boston, MA), and RNase inhibitor (Invitrogen). RNA was then degraded with RNase A, and labeled cDNAs were purified using QIAquick PCR columns (Qiagen)."
            assert record.entity_attributes["Sample_data_processing"] == "LOWESS normalized, background subtracted VALUE data obtained from log of processed Red signal/processed Green signal."
            assert record.entity_attributes["sample_table_end"] == ""
            assert record.entity_attributes["Sample_label_ch2"] == "Cy3"
            assert record.entity_attributes["Sample_description"] == "Biological replicate 1 of 4. Control embryonic stem cells, untreated, harvested after several passages."
            assert record.entity_attributes["Sample_source_name_ch1"] == "Total RNA from murine ES-D3 embryonic stem cells labeled with Cyanine-5 (red)."
            assert record.entity_attributes["Sample_source_name_ch2"] == "Total RNA from pooled whole mouse embryos e17.5, labeled with Cyanine-3 (green)."
            assert record.entity_attributes["Sample_molecule_ch1"] == "total RNA"
            assert record.entity_attributes["Sample_molecule_ch2"] == "total RNA"
            assert record.entity_attributes["Sample_growth_protocol_ch1"] == "ES cells were kept in an undifferentiated, pluripotent state by using 1000 IU/ml leukemia inhibitory factor (LIF; Chemicon, ESGRO, ESG1107), and grown on top of murine embryonic fibroblasts feeder layer inactivated by 10 ug/ml of mitomycin C (Sigma, St. Louis). ES cells were cultured on 0.1% gelatin-coated plastic dishes in ES medium containing Dulbecco modified Eagle medium supplemented with 15% fetal calf serum, 0.1 mM beta-mercaptoethanol, 2 mM glutamine, and 0.1 mN non-essential amino acids."
            assert len(record.entity_attributes["Sample_characteristics_ch1"]) == 4
            assert record.entity_attributes["Sample_characteristics_ch1"][0] == "ES-D3 cell line (CRL-1934)"
            assert record.entity_attributes["Sample_characteristics_ch1"][1] == "Age: day 4"
            assert record.entity_attributes["Sample_characteristics_ch1"][2] == "Tissue: blastocytes"
            assert record.entity_attributes["Sample_characteristics_ch1"][3] == "Strain: 129/Sv mice"
            assert len(record.entity_attributes["Sample_characteristics_ch2"]) == 3
            assert record.entity_attributes["Sample_characteristics_ch2"][0] == "Strain: C57BL/6"
            assert record.entity_attributes["Sample_characteristics_ch2"][1] == "Age: e17.5 d"
            assert record.entity_attributes["Sample_characteristics_ch2"][2] == "Tissue: whole embryo"
            assert len(record.col_defs) == 6
            assert record.col_defs["ID_REF"] == ""
            assert record.col_defs["VALUE"] == "log(REDsignal/GREENsignal) per feature (processed signals used)."
            assert record.col_defs["gProcessedSignal"] == 'Dye-normalized signal after surrogate "algorithm," green "channel," used for computation of log ratio.'
            assert record.col_defs["LogRatioError"] == "error of the log ratio calculated according to the error model chosen."
            assert record.col_defs["PValueLogRatio"] == "Significance level of the Log Ratio computed for a feature."
            assert record.col_defs["rProcessedSignal"] == 'Dye-normalized signal after surrogate "algorithm," red "channel," used for computation of log ratio.'
            assert len(record.table_rows) == 21
            assert len(record.table_rows[0]) == 6
            assert record.table_rows[0][0] == "ID_REF"
            assert record.table_rows[0][1] == "VALUE"
            assert record.table_rows[0][2] == "LogRatioError"
            assert record.table_rows[0][3] == "PValueLogRatio"
            assert record.table_rows[0][4] == "gProcessedSignal"
            assert record.table_rows[0][5] == "rProcessedSignal"
            assert len(record.table_rows[1]) == 6
            assert record.table_rows[1][0] == "1"
            assert record.table_rows[1][1] == "-1.6274758"
            assert record.table_rows[1][2] == "1.36E-01"
            assert record.table_rows[1][3] == "6.41E-33"
            assert record.table_rows[1][4] == "9.13E+03"
            assert record.table_rows[1][5] == "2.15E+02"
            assert len(record.table_rows[2]) == 6
            assert record.table_rows[2][0] == "2"
            assert record.table_rows[2][1] == "0.1412248"
            assert record.table_rows[2][2] == "1.34E+00"
            assert record.table_rows[2][3] == "1.00E+00"
            assert record.table_rows[2][4] == "4.14E+01"
            assert record.table_rows[2][5] == "5.72E+01"
            assert len(record.table_rows[3]) == 6
            assert record.table_rows[3][0] == "3"
            assert record.table_rows[3][1] == "0.1827684"
            assert record.table_rows[3][2] == "5.19E-02"
            assert record.table_rows[3][3] == "4.33E-04"
            assert record.table_rows[3][4] == "5.13E+03"
            assert record.table_rows[3][5] == "7.81E+03"
            assert len(record.table_rows[4]) == 6
            assert record.table_rows[4][0] == "4"
            assert record.table_rows[4][1] == "-0.3932267"
            assert record.table_rows[4][2] == "6.08E-02"
            assert record.table_rows[4][3] == "1.02E-10"
            assert record.table_rows[4][4] == "4.65E+03"
            assert record.table_rows[4][5] == "1.88E+03"
            assert len(record.table_rows[5]) == 6
            assert record.table_rows[5][0] == "5"
            assert record.table_rows[5][1] == "-0.9865994"
            assert record.table_rows[5][2] == "1.05E-01"
            assert record.table_rows[5][3] == "6.32E-21"
            assert record.table_rows[5][4] == "2.91E+03"
            assert record.table_rows[5][5] == "3.01E+02"
            assert len(record.table_rows[6]) == 6
            assert record.table_rows[6][0] == "6"
            assert record.table_rows[6][1] == "0.0238812"
            assert record.table_rows[6][2] == "1.02E-01"
            assert record.table_rows[6][3] == "8.15E-01"
            assert record.table_rows[6][4] == "7.08E+02"
            assert record.table_rows[6][5] == "7.48E+02"
            assert len(record.table_rows[7]) == 6
            assert record.table_rows[7][0] == "7"
            assert record.table_rows[7][1] == "-1.4841822"
            assert record.table_rows[7][2] == "1.25E-01"
            assert record.table_rows[7][3] == "1.42E-32"
            assert record.table_rows[7][4] == "1.02E+04"
            assert record.table_rows[7][5] == "3.36E+02"
            assert len(record.table_rows[8]) == 6
            assert record.table_rows[8][0] == "8"
            assert record.table_rows[8][1] == "-1.8261356"
            assert record.table_rows[8][2] == "4.15E-01"
            assert record.table_rows[8][3] == "1.10E-05"
            assert record.table_rows[8][4] == "7.19E+02"
            assert record.table_rows[8][5] == "1.07E+01"
            assert len(record.table_rows[9]) == 6
            assert record.table_rows[9][0] == "9"
            assert record.table_rows[9][1] == "-1.0344779"
            assert record.table_rows[9][2] == "1.78E+00"
            assert record.table_rows[9][3] == "1.00E+00"
            assert record.table_rows[9][4] == "9.62E+01"
            assert record.table_rows[9][5] == "8.89E+00"
            assert len(record.table_rows[10]) == 6
            assert record.table_rows[10][0] == "10"
            assert record.table_rows[10][1] == "0.2405891"
            assert record.table_rows[10][2] == "3.09E-01"
            assert record.table_rows[10][3] == "4.36E-01"
            assert record.table_rows[10][4] == "1.61E+02"
            assert record.table_rows[10][5] == "2.80E+02"
            assert len(record.table_rows[11]) == 6
            assert record.table_rows[11][0] == "11"
            assert record.table_rows[11][1] == "0.3209366"
            assert record.table_rows[11][2] == "3.59E-01"
            assert record.table_rows[11][3] == "3.71E-01"
            assert record.table_rows[11][4] == "1.25E+02"
            assert record.table_rows[11][5] == "2.61E+02"
            assert len(record.table_rows[12]) == 6
            assert record.table_rows[12][0] == "12"
            assert record.table_rows[12][1] == "0.358304"
            assert record.table_rows[12][2] == "2.06E+00"
            assert record.table_rows[12][3] == "1.00E+00"
            assert record.table_rows[12][4] == "2.04E+01"
            assert record.table_rows[12][5] == "4.66E+01"
            assert len(record.table_rows[13]) == 6
            assert record.table_rows[13][0] == "13"
            assert record.table_rows[13][1] == "-0.0122072"
            assert record.table_rows[13][2] == "3.64E-01"
            assert record.table_rows[13][3] == "9.73E-01"
            assert record.table_rows[13][4] == "1.84E+02"
            assert record.table_rows[13][5] == "1.79E+02"
            assert len(record.table_rows[14]) == 6
            assert record.table_rows[14][0] == "14"
            assert record.table_rows[14][1] == "-1.5480396"
            assert record.table_rows[14][2] == "1.30E-01"
            assert record.table_rows[14][3] == "7.21E-33"
            assert record.table_rows[14][4] == "1.02E+04"
            assert record.table_rows[14][5] == "2.90E+02"
            assert len(record.table_rows[15]) == 6
            assert record.table_rows[15][0] == "15"
            assert record.table_rows[15][1] == "0.0073419"
            assert record.table_rows[15][2] == "2.98E-01"
            assert record.table_rows[15][3] == "9.80E-01"
            assert record.table_rows[15][4] == "2.21E+02"
            assert record.table_rows[15][5] == "2.25E+02"
            assert len(record.table_rows[16]) == 6
            assert record.table_rows[16][0] == "16"
            assert record.table_rows[16][1] == "-0.2267015"
            assert record.table_rows[16][2] == "9.44E-01"
            assert record.table_rows[16][3] == "8.10E-01"
            assert record.table_rows[16][4] == "8.90E+01"
            assert record.table_rows[16][5] == "5.28E+01"
            assert len(record.table_rows[17]) == 6
            assert record.table_rows[17][0] == "17"
            assert record.table_rows[17][1] == "-0.1484023"
            assert record.table_rows[17][2] == "8.01E-01"
            assert record.table_rows[17][3] == "8.53E-01"
            assert record.table_rows[17][4] == "9.65E+01"
            assert record.table_rows[17][5] == "6.86E+01"
            assert len(record.table_rows[18]) == 6
            assert record.table_rows[18][0] == "18"
            assert record.table_rows[18][1] == "-0.6122195"
            assert record.table_rows[18][2] == "1.28E-01"
            assert record.table_rows[18][3] == "1.69E-06"
            assert record.table_rows[18][4] == "1.12E+03"
            assert record.table_rows[18][5] == "2.73E+02"
            assert len(record.table_rows[19]) == 6
            assert record.table_rows[19][0] == "19"
            assert record.table_rows[19][1] == "0.0796905"
            assert record.table_rows[19][2] == "8.78E-02"
            assert record.table_rows[19][3] == "3.64E-01"
            assert record.table_rows[19][4] == "8.21E+02"
            assert record.table_rows[19][5] == "9.87E+02"
            assert len(record.table_rows[20]) == 6
            assert record.table_rows[20][0] == "20"
            assert record.table_rows[20][1] == "-0.084895"
            assert record.table_rows[20][2] == "9.38E-01"
            assert record.table_rows[20][3] == "9.28E-01"
            assert record.table_rows[20][4] == "7.68E+01"
            assert record.table_rows[20][5] == "6.32E+01"
            record = next(records)
            assert record.entity_type == "SAMPLE"
            assert record.entity_id == "Control Embyronic Stem Cell Replicate 2"
            assert len(record.entity_attributes) == 24
            assert record.entity_attributes["Sample_extract_protocol_ch2"] == "TriZol procedure"
            assert record.entity_attributes["Sample_hyb_protocol"] == "Oligoarray control targets and hybridization buffer (Agilent In Situ Hybridization Kit Plus) were added, and samples were applied to microarrays enclosed in Agilent SureHyb-enabled hybridization chambers. After hybridization, slides were washed sequentially with 6x SSC/0.005% Triton X-102 and 0.1x SSC/0.005% Triton X-102 before scanning. Slides were hybridized for 17 h at 60\xb0C in a rotating oven, and washed."
            assert record.entity_attributes["Sample_extract_protocol_ch1"] == "TriZol procedure"
            assert record.entity_attributes["Sample_platform_id"] == "GPL3759"
            assert record.entity_attributes["Sample_title"] == "Control Embyronic Stem Cell Replicate 2"
            assert record.entity_attributes["Sample_supplementary_file"] == "file2.gpr"
            assert record.entity_attributes["Sample_organism_ch2"] == "Mus musculus"
            assert record.entity_attributes["Sample_organism_ch1"] == "Mus musculus"
            assert record.entity_attributes["Sample_label_ch1"] == "Cy5"
            assert len(record.entity_attributes["Sample_scan_protocol"]) == 2
            assert record.entity_attributes["Sample_scan_protocol"][0] == "Scanned on an Agilent G2565AA scanner."
            assert record.entity_attributes["Sample_scan_protocol"][1] == "Images were quantified using Agilent Feature Extraction Software (version A.7.5)."
            assert record.entity_attributes["sample_table_begin"] == ""
            assert record.entity_attributes["Sample_label_protocol_ch2"] == "10 \xb5g of total RNA were primed with 2 \xb5l of 100 \xb5M T16N2 DNA primer at 70\xb0C for 10 min, then reversed transcribed at 42\xb0C for 1 h in the presence of 400 U SuperScript II RTase (Invitrogen), and 100 \xb5M each dATP, dTTP, dGTP, with 25 \xb5M dCTP, 25 \xb5M Cy5-labeled dCTP (NEN Life Science, Boston, MA), and RNase inhibitor (Invitrogen). RNA was then degraded with RNase A, and labeled cDNAs were purified using QIAquick PCR columns (Qiagen)."
            assert record.entity_attributes["Sample_label_protocol_ch1"] == "10 \xb5g of total RNA were primed with 2 \xb5l of 100 \xb5M T16N2 DNA primer at 70\xb0C for 10 min, then reversed transcribed at 42\xb0C for 1 h in the presence of 400 U SuperScript II RTase (Invitrogen), and 100 \xb5M each dATP, dTTP, dGTP, with 25 \xb5M dCTP, 25 \xb5M Cy5-labeled dCTP (NEN Life Science, Boston, MA), and RNase inhibitor (Invitrogen). RNA was then degraded with RNase A, and labeled cDNAs were purified using QIAquick PCR columns (Qiagen)."
            assert record.entity_attributes["Sample_data_processing"] == "LOWESS normalized, background subtracted VALUE data obtained from log of processed Red signal/processed Green signal."
            assert record.entity_attributes["sample_table_end"] == ""
            assert record.entity_attributes["Sample_label_ch2"] == "Cy3"
            assert record.entity_attributes["Sample_description"] == "Biological replicate 2 of 4. Control embryonic stem cells, untreated, harvested after several passages."
            assert record.entity_attributes["Sample_source_name_ch1"] == "Total RNA from murine ES-D3 embryonic stem cells labeled with Cyanine-5 (red)."
            assert record.entity_attributes["Sample_source_name_ch2"] == "Total RNA from pooled whole mouse embryos e17.5, labeled with Cyanine-3 (green)."
            assert record.entity_attributes["Sample_molecule_ch1"] == "total RNA"
            assert record.entity_attributes["Sample_molecule_ch2"] == "total RNA"
            assert record.entity_attributes["Sample_growth_protocol_ch1"] == "ES cells were kept in an undifferentiated, pluripotent state by using 1000 IU/ml leukemia inhibitory factor (LIF; Chemicon, ESGRO, ESG1107), and grown on top of murine embryonic fibroblasts feeder layer inactivated by 10 ug/ml of mitomycin C (Sigma, St. Louis). ES cells were cultured on 0.1% gelatin-coated plastic dishes in ES medium containing Dulbecco modified Eagle medium supplemented with 15% fetal calf serum, 0.1 mM beta-mercaptoethanol, 2 mM glutamine, and 0.1 mN non-essential amino acids."
            assert len(record.entity_attributes["Sample_characteristics_ch1"]) == 4
            assert record.entity_attributes["Sample_characteristics_ch1"][0] == "ES-D3 cell line (CRL-1934)"
            assert record.entity_attributes["Sample_characteristics_ch1"][1] == "Age: day 4"
            assert record.entity_attributes["Sample_characteristics_ch1"][2] == "Tissue: blastocytes"
            assert record.entity_attributes["Sample_characteristics_ch1"][3] == "Strain: 129/Sv mice"
            assert len(record.entity_attributes["Sample_characteristics_ch2"]) == 3
            assert record.entity_attributes["Sample_characteristics_ch2"][0] == "Strain: C57BL/6"
            assert record.entity_attributes["Sample_characteristics_ch2"][1] == "Age: e17.5 d"
            assert record.entity_attributes["Sample_characteristics_ch2"][2] == "Tissue: whole embryo"
            assert len(record.col_defs) == 6
            assert record.col_defs["ID_REF"] == ""
            assert record.col_defs["VALUE"] == "log(REDsignal/GREENsignal) per feature (processed signals used)."
            assert record.col_defs["gProcessedSignal"] == 'Dye-normalized signal after surrogate "algorithm," green "channel," used for computation of log ratio.'
            assert record.col_defs["LogRatioError"] == "error of the log ratio calculated according to the error model chosen."
            assert record.col_defs["PValueLogRatio"] == "Significance level of the Log Ratio computed for a feature."
            assert record.col_defs["rProcessedSignal"] == 'Dye-normalized signal after surrogate "algorithm," red "channel," used for computation of log ratio.'
            assert len(record.table_rows) == 21
            assert len(record.table_rows[0]) == 6
            assert record.table_rows[0][0] == "ID_REF"
            assert record.table_rows[0][1] == "VALUE"
            assert record.table_rows[0][2] == "LogRatioError"
            assert record.table_rows[0][3] == "PValueLogRatio"
            assert record.table_rows[0][4] == "gProcessedSignal"
            assert record.table_rows[0][5] == "rProcessedSignal"
            assert len(record.table_rows[1]) == 6
            assert record.table_rows[1][0] == "1"
            assert record.table_rows[1][1] == "-1.1697263"
            assert record.table_rows[1][2] == "1.23E-01"
            assert record.table_rows[1][3] == "2.14E-21"
            assert record.table_rows[1][4] == "3.17E+03"
            assert record.table_rows[1][5] == "2.14E+02"
            assert len(record.table_rows[2]) == 6
            assert record.table_rows[2][0] == "2"
            assert record.table_rows[2][1] == "-0.1111353"
            assert record.table_rows[2][2] == "1.63E+00"
            assert record.table_rows[2][3] == "9.46E-01"
            assert record.table_rows[2][4] == "5.43E+01"
            assert record.table_rows[2][5] == "4.20E+01"
            assert len(record.table_rows[3]) == 6
            assert record.table_rows[3][0] == "3"
            assert record.table_rows[3][1] == "0.1400597"
            assert record.table_rows[3][2] == "5.11E-02"
            assert record.table_rows[3][3] == "6.17E-03"
            assert record.table_rows[3][4] == "6.72E+03"
            assert record.table_rows[3][5] == "9.28E+03"
            assert len(record.table_rows[4]) == 6
            assert record.table_rows[4][0] == "4"
            assert record.table_rows[4][1] == "-0.4820633"
            assert record.table_rows[4][2] == "6.38E-02"
            assert record.table_rows[4][3] == "4.06E-14"
            assert record.table_rows[4][4] == "6.46E+03"
            assert record.table_rows[4][5] == "2.13E+03"
            assert len(record.table_rows[5]) == 6
            assert record.table_rows[5][0] == "5"
            assert record.table_rows[5][1] == "-1.2116196"
            assert record.table_rows[5][2] == "1.22E-01"
            assert record.table_rows[5][3] == "2.31E-23"
            assert record.table_rows[5][4] == "3.62E+03"
            assert record.table_rows[5][5] == "2.22E+02"
            assert len(record.table_rows[6]) == 6
            assert record.table_rows[6][0] == "6"
            assert record.table_rows[6][1] == "-0.0230528"
            assert record.table_rows[6][2] == "1.04E-01"
            assert record.table_rows[6][3] == "8.24E-01"
            assert record.table_rows[6][4] == "8.76E+02"
            assert record.table_rows[6][5] == "8.31E+02"
            assert len(record.table_rows[7]) == 6
            assert record.table_rows[7][0] == "7"
            assert record.table_rows[7][1] == "-1.1380152"
            assert record.table_rows[7][2] == "1.13E-01"
            assert record.table_rows[7][3] == "9.23E-24"
            assert record.table_rows[7][4] == "3.94E+03"
            assert record.table_rows[7][5] == "2.86E+02"
            assert len(record.table_rows[8]) == 6
            assert record.table_rows[8][0] == "8"
            assert record.table_rows[8][1] == "-1.834596"
            assert record.table_rows[8][2] == "5.40E-01"
            assert record.table_rows[8][3] == "6.74E-04"
            assert record.table_rows[8][4] == "6.44E+02"
            assert record.table_rows[8][5] == "9.43E+00"
            assert len(record.table_rows[9]) == 6
            assert record.table_rows[9][0] == "9"
            assert record.table_rows[9][1] == "-0.9747637"
            assert record.table_rows[9][2] == "2.14E+00"
            assert record.table_rows[9][3] == "1.00E+00"
            assert record.table_rows[9][4] == "9.17E+01"
            assert record.table_rows[9][5] == "9.72E+00"
            assert len(record.table_rows[10]) == 6
            assert record.table_rows[10][0] == "10"
            assert record.table_rows[10][1] == "0.3874005"
            assert record.table_rows[10][2] == "2.92E-01"
            assert record.table_rows[10][3] == "1.85E-01"
            assert record.table_rows[10][4] == "1.69E+02"
            assert record.table_rows[10][5] == "4.11E+02"
            assert len(record.table_rows[11]) == 6
            assert record.table_rows[11][0] == "11"
            assert record.table_rows[11][1] == "0.5340442"
            assert record.table_rows[11][2] == "3.29E-01"
            assert record.table_rows[11][3] == "1.04E-01"
            assert record.table_rows[11][4] == "1.23E+02"
            assert record.table_rows[11][5] == "4.20E+02"
            assert len(record.table_rows[12]) == 6
            assert record.table_rows[12][0] == "12"
            assert record.table_rows[12][1] == "0.3260696"
            assert record.table_rows[12][2] == "1.92E+00"
            assert record.table_rows[12][3] == "8.65E-01"
            assert record.table_rows[12][4] == "2.73E+01"
            assert record.table_rows[12][5] == "5.77E+01"
            assert len(record.table_rows[13]) == 6
            assert record.table_rows[13][0] == "13"
            assert record.table_rows[13][1] == "0.3010618"
            assert record.table_rows[13][2] == "2.84E-01"
            assert record.table_rows[13][3] == "2.90E-01"
            assert record.table_rows[13][4] == "1.93E+02"
            assert record.table_rows[13][5] == "3.87E+02"
            assert len(record.table_rows[14]) == 6
            assert record.table_rows[14][0] == "14"
            assert record.table_rows[14][1] == "-1.0760413"
            assert record.table_rows[14][2] == "1.08E-01"
            assert record.table_rows[14][3] == "1.63E-23"
            assert record.table_rows[14][4] == "4.06E+03"
            assert record.table_rows[14][5] == "3.41E+02"
            assert len(record.table_rows[15]) == 6
            assert record.table_rows[15][0] == "15"
            assert record.table_rows[15][1] == "-0.1167371"
            assert record.table_rows[15][2] == "3.87E-01"
            assert record.table_rows[15][3] == "7.63E-01"
            assert record.table_rows[15][4] == "2.32E+02"
            assert record.table_rows[15][5] == "1.77E+02"
            assert len(record.table_rows[16]) == 6
            assert record.table_rows[16][0] == "16"
            assert record.table_rows[16][1] == "-0.1936322"
            assert record.table_rows[16][2] == "9.44E-01"
            assert record.table_rows[16][3] == "8.38E-01"
            assert record.table_rows[16][4] == "1.02E+02"
            assert record.table_rows[16][5] == "6.56E+01"
            assert len(record.table_rows[17]) == 6
            assert record.table_rows[17][0] == "17"
            assert record.table_rows[17][1] == "-0.3275898"
            assert record.table_rows[17][2] == "7.87E-01"
            assert record.table_rows[17][3] == "6.77E-01"
            assert record.table_rows[17][4] == "1.41E+02"
            assert record.table_rows[17][5] == "6.65E+01"
            assert len(record.table_rows[18]) == 6
            assert record.table_rows[18][0] == "18"
            assert record.table_rows[18][1] == "-0.4805853"
            assert record.table_rows[18][2] == "1.14E-01"
            assert record.table_rows[18][3] == "2.41E-05"
            assert record.table_rows[18][4] == "1.34E+03"
            assert record.table_rows[18][5] == "4.42E+02"
            assert len(record.table_rows[19]) == 6
            assert record.table_rows[19][0] == "19"
            assert record.table_rows[19][1] == "0.1109524"
            assert record.table_rows[19][2] == "9.56E-02"
            assert record.table_rows[19][3] == "2.46E-01"
            assert record.table_rows[19][4] == "8.38E+02"
            assert record.table_rows[19][5] == "1.08E+03"
            assert len(record.table_rows[20]) == 6
            assert record.table_rows[20][0] == "20"
            assert record.table_rows[20][1] == "0.1677912"
            assert record.table_rows[20][2] == "6.51E-01"
            assert record.table_rows[20][3] == "7.97E-01"
            assert record.table_rows[20][4] == "9.84E+01"
            assert record.table_rows[20][5] == "1.45E+02"
            record = next(records)
            assert record.entity_type == "SAMPLE"
            assert record.entity_id == "Triple-Fusion Transfected Embryonic Stem Cells Replicate 1"
            assert len(record.entity_attributes) == 25
            assert record.entity_attributes["Sample_extract_protocol_ch2"] == "TriZol procedure"
            assert record.entity_attributes["Sample_hyb_protocol"] == "Oligoarray control targets and hybridization buffer (Agilent In Situ Hybridization Kit Plus) were added, and samples were applied to microarrays enclosed in Agilent SureHyb-enabled hybridization chambers. After hybridization, slides were washed sequentially with 6x SSC/0.005% Triton X-102 and 0.1x SSC/0.005% Triton X-102 before scanning. Slides were hybridized for 17 h at 60\xb0C in a rotating oven, and washed."
            assert record.entity_attributes["Sample_extract_protocol_ch1"] == "TriZol procedure"
            assert record.entity_attributes["Sample_platform_id"] == "GPL3759"
            assert record.entity_attributes["Sample_title"] == "Triple-Fusion Transfected Embryonic Stem Cells Replicate 1"
            assert record.entity_attributes["Sample_supplementary_file"] == "file3.gpr"
            assert record.entity_attributes["Sample_organism_ch2"] == "Mus musculus"
            assert record.entity_attributes["Sample_organism_ch1"] == "Mus musculus"
            assert record.entity_attributes["Sample_label_ch1"] == "Cy5"
            assert record.entity_attributes["Sample_scan_protocol"] == "Scanned on an Agilent G2565AA scanner."
            assert record.entity_attributes["sample_table_begin"] == ""
            assert record.entity_attributes["Sample_label_protocol_ch2"] == "10 \xb5g of total RNA were primed with 2 \xb5l of 100 \xb5M T16N2 DNA primer at 70\xb0C for 10 min, then reversed transcribed at 42\xb0C for 1 h in the presence of 400 U SuperScript II RTase (Invitrogen), and 100 \xb5M each dATP, dTTP, dGTP, with 25 \xb5M dCTP, 25 \xb5M Cy5-labeled dCTP (NEN Life Science, Boston, MA), and RNase inhibitor (Invitrogen). RNA was then degraded with RNase A, and labeled cDNAs were purified using QIAquick PCR columns (Qiagen)."
            assert record.entity_attributes["Sample_label_protocol_ch1"] == "10 \xb5g of total RNA were primed with 2 \xb5l of 100 \xb5M T16N2 DNA primer at 70\xb0C for 10 min, then reversed transcribed at 42\xb0C for 1 h in the presence of 400 U SuperScript II RTase (Invitrogen), and 100 \xb5M each dATP, dTTP, dGTP with 25 \xb5M dCTP, 25 \xb5M Cy5-labeled dCTP (NEN Life Science, Boston, MA), and RNase inhibitor (Invitrogen). RNA was then degraded with RNase A, and labeled cDNAs were purified using QIAquick PCR columns (Qiagen)."
            assert record.entity_attributes["Sample_data_processing"] == "LOWESS normalized, background subtracted VALUE data obtained from log of processed Red signal/processed Green signal."
            assert record.entity_attributes["sample_table_end"] == ""
            assert record.entity_attributes["Sample_label_ch2"] == "Cy3"
            assert record.entity_attributes["Sample_description"] == "Biological replicate 1 of 3. Stable triple-fusion-reporter-gene transfected embryonic stem cells, harvested after several passages."
            assert record.entity_attributes["Sample_source_name_ch1"] == "Total RNA from murine ES-D3 triple-transfected embryonic stem cells labeled with Cyanine-5 (red)."
            assert record.entity_attributes["Sample_source_name_ch2"] == "Total RNA from pooled whole mouse embryos e17.5, labeled with Cyanine-3 (green)."
            assert record.entity_attributes["Sample_molecule_ch1"] == "total RNA"
            assert record.entity_attributes["Sample_molecule_ch2"] == "total RNA"
            assert record.entity_attributes["Sample_growth_protocol_ch1"] == "ES cells were kept in an undifferentiated, pluripotent state by using 1000 IU/ml leukemia inhibitory factor (LIF; Chemicon, ESGRO, ESG1107), and grown on top of murine embryonic fibroblasts feeder layer inactivated by 10 ug/ml of mitomycin C (Sigma, St. Louis). ES cells were cultured on 0.1% gelatin-coated plastic dishes in ES medium containing Dulbecco modified Eagle medium supplemented with 15% fetal calf serum, 0.1 mM beta-mercaptoethanol, 2 mM glutamine, and 0.1 mN non-essential amino acids."
            assert len(record.entity_attributes["Sample_characteristics_ch1"]) == 5
            assert record.entity_attributes["Sample_characteristics_ch1"][0] == "ES-D3 cell line (CRL-1934)"
            assert record.entity_attributes["Sample_characteristics_ch1"][1] == "Transfected with pUb-fluc-mrfp-ttk triple fusion reporter gene."
            assert record.entity_attributes["Sample_characteristics_ch1"][2] == "Age: day 4"
            assert record.entity_attributes["Sample_characteristics_ch1"][3] == "Tissue: blastocytes"
            assert record.entity_attributes["Sample_characteristics_ch1"][4] == "Strain: 129/Sv mice"
            assert len(record.entity_attributes["Sample_characteristics_ch2"]) == 3
            assert record.entity_attributes["Sample_characteristics_ch2"][0] == "Strain: C57BL/6"
            assert record.entity_attributes["Sample_characteristics_ch2"][1] == "Age: e17.5 d"
            assert record.entity_attributes["Sample_characteristics_ch2"][2] == "Tissue: whole embryo"
            assert record.entity_attributes["Sample_treatment_protocol_ch1"] == "PCR amplification and standard cloning techniques were used to insert fluc and mrfp genes from plasmids pCDNA 3.1-CMV-fluc (Promega, Madison, WI) and pCDNA3.1-CMV-mrfp in frame with the ttk gene into the pCDNA3.1-truncated	sr39tk. This triple fusion (TF) reporter gene fragment (3.3 kbp) was released from the plasmid with Not1 and BamH1 restriction enzymes before blunt-end ligation into the multiple cloning site of lentiviral transfer vector, FUG, driven by the human ubiquitin-C promoter. Self-inactivating (SIN) lentivirus was prepared by transient transfection of 293T cells. Briefly, pFUG-TF containing the triple fusion reporter gene was co-transfected into 293T cells with HIV-1 packaging vector (?8.9) and vesicular stomatitis virus G glycoprotein-pseudotyped envelop vector (pVSVG). Lentivirus supernatant was concentrated by sediment centrifugation using a SW29 rotor at 50,000 x g for two hours. Concentrated virus was titered on 293T cells. Murine ES cells were transfected with LV-pUb-fluc-mrfp-ttk at a multiplicity of infection (MOI) of 10."
            assert len(record.col_defs) == 6
            assert record.col_defs["ID_REF"] == ""
            assert record.col_defs["VALUE"] == "log(REDsignal/GREENsignal) per feature (processed signals used)."
            assert record.col_defs["gProcessedSignal"] == 'Dye-normalized signal after surrogate "algorithm," green "channel," used for computation of log ratio.'
            assert record.col_defs["LogRatioError"] == "error of the log ratio calculated according to the error model chosen."
            assert record.col_defs["PValueLogRatio"] == "Significance level of the Log Ratio computed for a feature."
            assert record.col_defs["rProcessedSignal"] == 'Dye-normalized signal after surrogate "algorithm," red "channel," used for computation of log ratio.'
            assert len(record.table_rows) == 21
            assert len(record.table_rows[0]) == 6
            assert record.table_rows[0][0] == "ID_REF"
            assert record.table_rows[0][1] == "VALUE"
            assert record.table_rows[0][2] == "LogRatioError"
            assert record.table_rows[0][3] == "PValueLogRatio"
            assert record.table_rows[0][4] == "gProcessedSignal"
            assert record.table_rows[0][5] == "rProcessedSignal"
            assert len(record.table_rows[1]) == 6
            assert record.table_rows[1][0] == "1"
            assert record.table_rows[1][1] == "-0.7837546"
            assert record.table_rows[1][2] == "1.30E-01"
            assert record.table_rows[1][3] == "1.70E-09"
            assert record.table_rows[1][4] == "2.10E+03"
            assert record.table_rows[1][5] == "3.46E+02"
            assert len(record.table_rows[2]) == 6
            assert record.table_rows[2][0] == "2"
            assert record.table_rows[2][1] == "0.3797837"
            assert record.table_rows[2][2] == "1.15E+00"
            assert record.table_rows[2][3] == "7.41E-01"
            assert record.table_rows[2][4] == "5.59E+01"
            assert record.table_rows[2][5] == "1.34E+02"
            assert len(record.table_rows[3]) == 6
            assert record.table_rows[3][0] == "3"
            assert record.table_rows[3][1] == "0.2079269"
            assert record.table_rows[3][2] == "5.38E-02"
            assert record.table_rows[3][3] == "1.12E-04"
            assert record.table_rows[3][4] == "5.04E+03"
            assert record.table_rows[3][5] == "8.14E+03"
            assert len(record.table_rows[4]) == 6
            assert record.table_rows[4][0] == "4"
            assert record.table_rows[4][1] == "-0.4730291"
            assert record.table_rows[4][2] == "6.71E-02"
            assert record.table_rows[4][3] == "1.86E-12"
            assert record.table_rows[4][4] == "5.66E+03"
            assert record.table_rows[4][5] == "1.91E+03"
            assert len(record.table_rows[5]) == 6
            assert record.table_rows[5][0] == "5"
            assert record.table_rows[5][1] == "-0.9481128"
            assert record.table_rows[5][2] == "1.19E-01"
            assert record.table_rows[5][3] == "1.30E-15"
            assert record.table_rows[5][4] == "3.10E+03"
            assert record.table_rows[5][5] == "3.49E+02"
            assert len(record.table_rows[6]) == 6
            assert record.table_rows[6][0] == "6"
            assert record.table_rows[6][1] == "-0.0159867"
            assert record.table_rows[6][2] == "1.33E-01"
            assert record.table_rows[6][3] == "9.05E-01"
            assert record.table_rows[6][4] == "8.45E+02"
            assert record.table_rows[6][5] == "8.14E+02"
            assert len(record.table_rows[7]) == 6
            assert record.table_rows[7][0] == "7"
            assert record.table_rows[7][1] == "-0.819922"
            assert record.table_rows[7][2] == "1.14E-01"
            assert record.table_rows[7][3] == "7.01E-13"
            assert record.table_rows[7][4] == "2.75E+03"
            assert record.table_rows[7][5] == "4.16E+02"
            assert len(record.table_rows[8]) == 6
            assert record.table_rows[8][0] == "8"
            assert record.table_rows[8][1] == "-0.1559774"
            assert record.table_rows[8][2] == "9.16E-01"
            assert record.table_rows[8][3] == "8.65E-01"
            assert record.table_rows[8][4] == "1.34E+02"
            assert record.table_rows[8][5] == "9.34E+01"
            assert len(record.table_rows[9]) == 6
            assert record.table_rows[9][0] == "9"
            assert record.table_rows[9][1] == "0.145267"
            assert record.table_rows[9][2] == "3.90E+00"
            assert record.table_rows[9][3] == "1.00E+00"
            assert record.table_rows[9][4] == "2.22E+01"
            assert record.table_rows[9][5] == "3.10E+01"
            assert len(record.table_rows[10]) == 6
            assert record.table_rows[10][0] == "10"
            assert record.table_rows[10][1] == "0.3611211"
            assert record.table_rows[10][2] == "3.40E-01"
            assert record.table_rows[10][3] == "2.88E-01"
            assert record.table_rows[10][4] == "1.97E+02"
            assert record.table_rows[10][5] == "4.52E+02"
            assert len(record.table_rows[11]) == 6
            assert record.table_rows[11][0] == "11"
            assert record.table_rows[11][1] == "0.5092089"
            assert record.table_rows[11][2] == "4.39E-01"
            assert record.table_rows[11][3] == "2.46E-01"
            assert record.table_rows[11][4] == "1.24E+02"
            assert record.table_rows[11][5] == "4.01E+02"
            assert len(record.table_rows[12]) == 6
            assert record.table_rows[12][0] == "12"
            assert record.table_rows[12][1] == "0.3715387"
            assert record.table_rows[12][2] == "1.69E+00"
            assert record.table_rows[12][3] == "8.26E-01"
            assert record.table_rows[12][4] == "3.84E+01"
            assert record.table_rows[12][5] == "9.04E+01"
            assert len(record.table_rows[13]) == 6
            assert record.table_rows[13][0] == "13"
            assert record.table_rows[13][1] == "0.1734934"
            assert record.table_rows[13][2] == "3.57E-01"
            assert record.table_rows[13][3] == "6.27E-01"
            assert record.table_rows[13][4] == "2.37E+02"
            assert record.table_rows[13][5] == "3.53E+02"
            assert len(record.table_rows[14]) == 6
            assert record.table_rows[14][0] == "14"
            assert record.table_rows[14][1] == "-0.9340707"
            assert record.table_rows[14][2] == "1.20E-01"
            assert record.table_rows[14][3] == "6.90E-15"
            assert record.table_rows[14][4] == "2.96E+03"
            assert record.table_rows[14][5] == "3.45E+02"
            assert len(record.table_rows[15]) == 6
            assert record.table_rows[15][0] == "15"
            assert record.table_rows[15][1] == "-0.2956317"
            assert record.table_rows[15][2] == "5.78E-01"
            assert record.table_rows[15][3] == "6.09E-01"
            assert record.table_rows[15][4] == "2.46E+02"
            assert record.table_rows[15][5] == "1.25E+02"
            assert len(record.table_rows[16]) == 6
            assert record.table_rows[16][0] == "16"
            assert record.table_rows[16][1] == "-0.2321102"
            assert record.table_rows[16][2] == "1.22E+00"
            assert record.table_rows[16][3] == "8.49E-01"
            assert record.table_rows[16][4] == "1.09E+02"
            assert record.table_rows[16][5] == "6.37E+01"
            assert len(record.table_rows[17]) == 6
            assert record.table_rows[17][0] == "17"
            assert record.table_rows[17][1] == "-0.1603561"
            assert record.table_rows[17][2] == "1.16E+00"
            assert record.table_rows[17][3] == "8.90E-01"
            assert record.table_rows[17][4] == "1.06E+02"
            assert record.table_rows[17][5] == "7.34E+01"
            assert len(record.table_rows[18]) == 6
            assert record.table_rows[18][0] == "18"
            assert record.table_rows[18][1] == "-0.5063897"
            assert record.table_rows[18][2] == "1.63E-01"
            assert record.table_rows[18][3] == "1.95E-03"
            assert record.table_rows[18][4] == "1.15E+03"
            assert record.table_rows[18][5] == "3.58E+02"
            assert len(record.table_rows[19]) == 6
            assert record.table_rows[19][0] == "19"
            assert record.table_rows[19][1] == "0.1990761"
            assert record.table_rows[19][2] == "1.32E-01"
            assert record.table_rows[19][3] == "1.32E-01"
            assert record.table_rows[19][4] == "6.65E+02"
            assert record.table_rows[19][5] == "1.05E+03"
            assert len(record.table_rows[20]) == 6
            assert record.table_rows[20][0] == "20"
            assert record.table_rows[20][1] == "0.2985912"
            assert record.table_rows[20][2] == "8.89E-01"
            assert record.table_rows[20][3] == "7.37E-01"
            assert record.table_rows[20][4] == "8.06E+01"
            assert record.table_rows[20][5] == "1.60E+02"

    def test_soft_ex_affy(self):
        path = "Geo/soft_ex_affy.txt"
        with open(path, encoding="latin") as handle:
            records = Geo.parse(handle)
            record = next(records)
            assert record.entity_type == "SAMPLE"
            assert record.entity_id == "Drosophila_T0-1"
            assert len(record.entity_attributes) == 18
            assert record.entity_attributes["Sample_organism_ch1"] == "Drosophila melanogaster"
            assert record.entity_attributes["Sample_label_ch1"] == "biotin"
            assert record.entity_attributes["Sample_description"] == "Gene expression data from embryos younger than nuclear cycle 9, i.e. before zygotic genome activation."
            assert record.entity_attributes["Sample_table_end"] == ""
            assert record.entity_attributes["Sample_growth_protocol_ch1"] == "30 min egg collections of OreR and yw flies at 25C were aged at room temperature (RT) according to the different temporal classes T0-T4."
            assert len(record.entity_attributes["Sample_characteristics_ch1"]) == 2
            assert record.entity_attributes["Sample_characteristics_ch1"][0] == "Genotype: yellow white and Oregon R parents"
            assert record.entity_attributes["Sample_characteristics_ch1"][1] == "Age: embryos younger than nuclear cycle 9, i.e. before pole cells budding"
            assert record.entity_attributes["Sample_scan_protocol"] == "GeneChips were scanned using the Hewlett-Packard GeneArray Scanner G2500A."
            assert record.entity_attributes["Sample_hyb_protocol"] == "Following fragmentation, 10 microg of cRNA were hybridized for 16 hr at 45C on GeneChip Drosophila Genome Array. GeneChips were washed and stained in the Affymetrix Fluidics Station 400."
            assert record.entity_attributes["Sample_extract_protocol_ch1"] == "Trizol extraction of total RNA was performed according to the manufacturer's instructions."
            assert record.entity_attributes["Sample_source_name_ch1"] == "Drosophila embryos before nuclear cycle 9 (maternal transcripts)"
            assert record.entity_attributes["Sample_table_begin"] == ""
            assert record.entity_attributes["Sample_molecule_ch1"] == "total RNA"
            assert record.entity_attributes["Sample_label_protocol_ch1"] == "Biotinylated cRNA were prepared according to the standard Affymetrix protocol from 6 microg total RNA (Expression Analysis Technical Manual, 2001, Affymetrix)."
            assert record.entity_attributes["Sample_data_processing"] == "The data were analyzed with Microarray Suite version 5.0 (MAS 5.0) using Affymetrix default analysis settings and global scaling as normalization method. The trimmed mean target intensity of each array was arbitrarily set to 100."
            assert record.entity_attributes["Sample_treatment_protocol_ch1"] == "Embryos were dechorionated with 50% bleach, put on a cover slip and covered with Halocarbon oil 27 (Sigma). Embryos of the appropriate stage were manually selected under the dissecting scope. Selected embryos were transferred to a basket, rinsed with PBS with 0,7% NaCl, 0,04% triton-X100 and placed on ice in the Trizol solution (GibcoBRL)."
            assert record.entity_attributes["Sample_title"] == "embryo at T0, biological rep1"
            assert record.entity_attributes["Sample_supplementary_file"] == "Drosophila_T0-1.CEL"
            assert record.entity_attributes["Sample_platform_id"] == "GPL72"
            assert len(record.col_defs) == 4
            assert record.col_defs["DETECTION P-VALUE"] == "'detection p-value', p-value that indicates the significance level of the detection call"
            assert record.col_defs["ID_REF"] == ""
            assert record.col_defs["VALUE"] == "MAS5-calculated Signal intensity"
            assert record.col_defs["ABS_CALL"] == "the call in an absolute analysis that indicates if the transcript was present (P), absent (A), marginal (M), or no call (NC)"
            assert len(record.table_rows) == 22
            assert len(record.table_rows[0]) == 4
            assert record.table_rows[0][0] == "ID_REF"
            assert record.table_rows[0][1] == "VALUE"
            assert record.table_rows[0][2] == "ABS_CALL"
            assert record.table_rows[0][3] == "DETECTION P-VALUE"
            assert len(record.table_rows[1]) == 4
            assert record.table_rows[1][0] == "141200_at"
            assert record.table_rows[1][1] == "36.6"
            assert record.table_rows[1][2] == "A"
            assert record.table_rows[1][3] == "0.818657"
            assert len(record.table_rows[2]) == 4
            assert record.table_rows[2][0] == "141201_at"
            assert record.table_rows[2][1] == "41.5"
            assert record.table_rows[2][2] == "A"
            assert record.table_rows[2][3] == "0.703191"
            assert len(record.table_rows[3]) == 4
            assert record.table_rows[3][0] == "141202_at"
            assert record.table_rows[3][1] == "607.3"
            assert record.table_rows[3][2] == "P"
            assert record.table_rows[3][3] == "0.000944"
            assert len(record.table_rows[4]) == 4
            assert record.table_rows[4][0] == "141203_at"
            assert record.table_rows[4][1] == "1509.1"
            assert record.table_rows[4][2] == "P"
            assert record.table_rows[4][3] == "0.000762"
            assert len(record.table_rows[5]) == 4
            assert record.table_rows[5][0] == "141204_at"
            assert record.table_rows[5][1] == "837.3"
            assert record.table_rows[5][2] == "P"
            assert record.table_rows[5][3] == "0.000613"
            assert len(record.table_rows[6]) == 4
            assert record.table_rows[6][0] == "141205_at"
            assert record.table_rows[6][1] == "363.2"
            assert record.table_rows[6][2] == "P"
            assert record.table_rows[6][3] == "0.003815"
            assert len(record.table_rows[7]) == 4
            assert record.table_rows[7][0] == "141206_at"
            assert record.table_rows[7][1] == "1193.6"
            assert record.table_rows[7][2] == "P"
            assert record.table_rows[7][3] == "0.000491"
            assert len(record.table_rows[8]) == 4
            assert record.table_rows[8][0] == "141207_at"
            assert record.table_rows[8][1] == "346.6"
            assert record.table_rows[8][2] == "P"
            assert record.table_rows[8][3] == "0.001165"
            assert len(record.table_rows[9]) == 4
            assert record.table_rows[9][0] == "141208_at"
            assert record.table_rows[9][1] == "257.8"
            assert record.table_rows[9][2] == "P"
            assert record.table_rows[9][3] == "0.006575"
            assert len(record.table_rows[10]) == 4
            assert record.table_rows[10][0] == "141209_at"
            assert record.table_rows[10][1] == "337.1"
            assert record.table_rows[10][2] == "P"
            assert record.table_rows[10][3] == "0.002607"
            assert len(record.table_rows[11]) == 4
            assert record.table_rows[11][0] == "141210_at"
            assert record.table_rows[11][1] == "48"
            assert record.table_rows[11][2] == "A"
            assert record.table_rows[11][3] == "0.150145"
            assert len(record.table_rows[12]) == 4
            assert record.table_rows[12][0] == "141211_at"
            assert record.table_rows[12][1] == "130.7"
            assert record.table_rows[12][2] == "P"
            assert record.table_rows[12][3] == "0.005504"
            assert len(record.table_rows[13]) == 4
            assert record.table_rows[13][0] == "141212_at"
            assert record.table_rows[13][1] == "1454.3"
            assert record.table_rows[13][2] == "P"
            assert record.table_rows[13][3] == "0.000491"
            assert len(record.table_rows[14]) == 4
            assert record.table_rows[14][0] == "141213_at"
            assert record.table_rows[14][1] == "21.2"
            assert record.table_rows[14][2] == "A"
            assert record.table_rows[14][3] == "0.635055"
            assert len(record.table_rows[15]) == 4
            assert record.table_rows[15][0] == "142121_at"
            assert record.table_rows[15][1] == "133.7"
            assert record.table_rows[15][2] == "A"
            assert record.table_rows[15][3] == "0.889551"
            assert len(record.table_rows[16]) == 4
            assert record.table_rows[16][0] == "142122_at"
            assert record.table_rows[16][1] == "275.3"
            assert record.table_rows[16][2] == "A"
            assert record.table_rows[16][3] == "0.611218"
            assert len(record.table_rows[17]) == 4
            assert record.table_rows[17][0] == "142123_at"
            assert record.table_rows[17][1] == "307.6"
            assert record.table_rows[17][2] == "A"
            assert record.table_rows[17][3] == "0.611218"
            assert len(record.table_rows[18]) == 4
            assert record.table_rows[18][0] == "142124_at"
            assert record.table_rows[18][1] == "132.6"
            assert record.table_rows[18][2] == "A"
            assert record.table_rows[18][3] == "0.437646"
            assert len(record.table_rows[19]) == 4
            assert record.table_rows[19][0] == "142125_at"
            assert record.table_rows[19][1] == "195.8"
            assert record.table_rows[19][2] == "A"
            assert record.table_rows[19][3] == "0.110449"
            assert len(record.table_rows[20]) == 4
            assert record.table_rows[20][0] == "142126_at"
            assert record.table_rows[20][1] == "174.1"
            assert record.table_rows[20][2] == "A"
            assert record.table_rows[20][3] == "0.681117"
            assert len(record.table_rows[21]) == 4
            assert record.table_rows[21][0] == "142127_at"
            assert record.table_rows[21][1] == "316.3"
            assert record.table_rows[21][2] == "A"
            assert record.table_rows[21][3] == "0.65838"
            record = next(records)
            assert record.entity_type == "SAMPLE"
            assert record.entity_id == "Drosophila_T0-2"
            assert len(record.entity_attributes) == 18
            assert record.entity_attributes["Sample_organism_ch1"] == "Drosophila melanogaster"
            assert record.entity_attributes["Sample_label_ch1"] == "biotin"
            assert record.entity_attributes["Sample_description"] == "Gene expression data from embryos younger than nuclear cycle 9, i.e. before zygotic genome activation."
            assert record.entity_attributes["Sample_table_end"] == ""
            assert record.entity_attributes["Sample_growth_protocol_ch1"] == "30 min egg collections of OreR and yw flies at 25C were aged at room temperature (RT) according to the different temporal classes T0-T4."
            assert len(record.entity_attributes["Sample_characteristics_ch1"]) == 2
            assert record.entity_attributes["Sample_characteristics_ch1"][0] == "Genotype: yellow white and Oregon R parents"
            assert record.entity_attributes["Sample_characteristics_ch1"][1] == "Age: embryos younger than nuclear cycle 9, i.e. before pole cells budding"
            assert record.entity_attributes["Sample_scan_protocol"] == "GeneChips were scanned using the Hewlett-Packard GeneArray Scanner G2500A."
            assert record.entity_attributes["Sample_hyb_protocol"] == "Following fragmentation, 10 microg of cRNA were hybridized for 16 hr at 45C on GeneChip Drosophila Genome Array. GeneChips were washed and stained in the Affymetrix Fluidics Station 400."
            assert record.entity_attributes["Sample_extract_protocol_ch1"] == "Trizol extraction of total RNA was performed according to the manufacturer's instructions."
            assert record.entity_attributes["Sample_source_name_ch1"] == "Drosophila embryos before nuclear cycle 9 (maternal transcripts)"
            assert record.entity_attributes["Sample_table_begin"] == ""
            assert record.entity_attributes["Sample_molecule_ch1"] == "total RNA"
            assert record.entity_attributes["Sample_label_protocol_ch1"] == "Biotinylated cRNA were prepared according to the standard Affymetrix protocol from 6 microg total RNA (Expression Analysis Technical Manual, 2001, Affymetrix)."
            assert record.entity_attributes["Sample_data_processing"] == "The data were analyzed with Microarray Suite version 5.0 (MAS 5.0) using Affymetrix default analysis settings and global scaling as normalization method. The trimmed mean target intensity of each array was arbitrarily set to 100."
            assert record.entity_attributes["Sample_treatment_protocol_ch1"] == "Embryos were dechorionated with 50% bleach, put on a cover slip and covered with Halocarbon oil 27 (Sigma). Embryos of the appropriate stage were manually selected under the dissecting scope. Selected embryos were transferred to a basket, rinsed with PBS with 0,7% NaCl, 0,04% triton-X100 and placed on ice in the Trizol solution (GibcoBRL)."
            assert record.entity_attributes["Sample_title"] == "embryo at T0, biological rep2"
            assert record.entity_attributes["Sample_supplementary_file"] == "Drosophila_T0-2.CEL"
            assert record.entity_attributes["Sample_platform_id"] == "GPL72"
            assert len(record.col_defs) == 4
            assert record.col_defs["DETECTION P-VALUE"] == "'detection p-value', p-value that indicates the significance level of the detection call"
            assert record.col_defs["ID_REF"] == ""
            assert record.col_defs["VALUE"] == "MAS5-calculated Signal intensity"
            assert record.col_defs["ABS_CALL"] == "the call in an absolute analysis that indicates if the transcript was present (P), absent (A), marginal (M), or no call (NC)"
            assert len(record.table_rows) == 22
            assert len(record.table_rows[0]) == 4
            assert record.table_rows[0][0] == "ID_REF"
            assert record.table_rows[0][1] == "VALUE"
            assert record.table_rows[0][2] == "ABS_CALL"
            assert record.table_rows[0][3] == "DETECTION P-VALUE"
            assert len(record.table_rows[1]) == 4
            assert record.table_rows[1][0] == "141200_at"
            assert record.table_rows[1][1] == "70.3"
            assert record.table_rows[1][2] == "A"
            assert record.table_rows[1][3] == "0.216313"
            assert len(record.table_rows[2]) == 4
            assert record.table_rows[2][0] == "141201_at"
            assert record.table_rows[2][1] == "38"
            assert record.table_rows[2][2] == "A"
            assert record.table_rows[2][3] == "0.635055"
            assert len(record.table_rows[3]) == 4
            assert record.table_rows[3][0] == "141202_at"
            assert record.table_rows[3][1] == "831.8"
            assert record.table_rows[3][2] == "P"
            assert record.table_rows[3][3] == "0.000613"
            assert len(record.table_rows[4]) == 4
            assert record.table_rows[4][0] == "141203_at"
            assert record.table_rows[4][1] == "2215.5"
            assert record.table_rows[4][2] == "P"
            assert record.table_rows[4][3] == "0.000944"
            assert len(record.table_rows[5]) == 4
            assert record.table_rows[5][0] == "141204_at"
            assert record.table_rows[5][1] == "965.6"
            assert record.table_rows[5][2] == "P"
            assert record.table_rows[5][3] == "0.000491"
            assert len(record.table_rows[6]) == 4
            assert record.table_rows[6][0] == "141205_at"
            assert record.table_rows[6][1] == "383.2"
            assert record.table_rows[6][2] == "P"
            assert record.table_rows[6][3] == "0.001432"
            assert len(record.table_rows[7]) == 4
            assert record.table_rows[7][0] == "141206_at"
            assert record.table_rows[7][1] == "1195"
            assert record.table_rows[7][2] == "P"
            assert record.table_rows[7][3] == "0.000491"
            assert len(record.table_rows[8]) == 4
            assert record.table_rows[8][0] == "141207_at"
            assert record.table_rows[8][1] == "413.7"
            assert record.table_rows[8][2] == "P"
            assert record.table_rows[8][3] == "0.000613"
            assert len(record.table_rows[9]) == 4
            assert record.table_rows[9][0] == "141208_at"
            assert record.table_rows[9][1] == "447.3"
            assert record.table_rows[9][2] == "P"
            assert record.table_rows[9][3] == "0.000762"
            assert len(record.table_rows[10]) == 4
            assert record.table_rows[10][0] == "141209_at"
            assert record.table_rows[10][1] == "294.4"
            assert record.table_rows[10][2] == "P"
            assert record.table_rows[10][3] == "0.004591"
            assert len(record.table_rows[11]) == 4
            assert record.table_rows[11][0] == "141210_at"
            assert record.table_rows[11][1] == "81.7"
            assert record.table_rows[11][2] == "M"
            assert record.table_rows[11][3] == "0.054711"
            assert len(record.table_rows[12]) == 4
            assert record.table_rows[12][0] == "141211_at"
            assert record.table_rows[12][1] == "84.9"
            assert record.table_rows[12][2] == "P"
            assert record.table_rows[12][3] == "0.005504"
            assert len(record.table_rows[13]) == 4
            assert record.table_rows[13][0] == "141212_at"
            assert record.table_rows[13][1] == "1456.4"
            assert record.table_rows[13][2] == "P"
            assert record.table_rows[13][3] == "0.000491"
            assert len(record.table_rows[14]) == 4
            assert record.table_rows[14][0] == "141213_at"
            assert record.table_rows[14][1] == "37"
            assert record.table_rows[14][2] == "A"
            assert record.table_rows[14][3] == "0.122747"
            assert len(record.table_rows[15]) == 4
            assert record.table_rows[15][0] == "142121_at"
            assert record.table_rows[15][1] == "133.7"
            assert record.table_rows[15][2] == "A"
            assert record.table_rows[15][3] == "0.889551"
            assert len(record.table_rows[16]) == 4
            assert record.table_rows[16][0] == "142122_at"
            assert record.table_rows[16][1] == "275.3"
            assert record.table_rows[16][2] == "A"
            assert record.table_rows[16][3] == "0.611218"
            assert len(record.table_rows[17]) == 4
            assert record.table_rows[17][0] == "142123_at"
            assert record.table_rows[17][1] == "307.6"
            assert record.table_rows[17][2] == "A"
            assert record.table_rows[17][3] == "0.611218"
            assert len(record.table_rows[18]) == 4
            assert record.table_rows[18][0] == "142124_at"
            assert record.table_rows[18][1] == "132.6"
            assert record.table_rows[18][2] == "A"
            assert record.table_rows[18][3] == "0.437646"
            assert len(record.table_rows[19]) == 4
            assert record.table_rows[19][0] == "142125_at"
            assert record.table_rows[19][1] == "195.8"
            assert record.table_rows[19][2] == "A"
            assert record.table_rows[19][3] == "0.110449"
            assert len(record.table_rows[20]) == 4
            assert record.table_rows[20][0] == "142126_at"
            assert record.table_rows[20][1] == "174.1"
            assert record.table_rows[20][2] == "A"
            assert record.table_rows[20][3] == "0.681117"
            assert len(record.table_rows[21]) == 4
            assert record.table_rows[21][0] == "142127_at"
            assert record.table_rows[21][1] == "316.3"
            assert record.table_rows[21][2] == "A"
            assert record.table_rows[21][3] == "0.65838"
            record = next(records)
            assert record.entity_type == "SAMPLE"
            assert record.entity_id == "Drosophila_T1-1"
            assert len(record.entity_attributes) == 18
            assert record.entity_attributes["Sample_organism_ch1"] == "Drosophila melanogaster"
            assert record.entity_attributes["Sample_label_ch1"] == "biotin"
            assert record.entity_attributes["Sample_description"] == "Gene expression data from embryos in slow phase of cellularisation."
            assert record.entity_attributes["Sample_table_end"] == ""
            assert record.entity_attributes["Sample_growth_protocol_ch1"] == "30 min egg collections of OreR and yw flies at 25C were aged at room temperature (RT) according to the different temporal classes T0-T4."
            assert len(record.entity_attributes["Sample_characteristics_ch1"]) == 2
            assert record.entity_attributes["Sample_characteristics_ch1"][0] == "Genotype: yellow white and Oregon R parents"
            assert record.entity_attributes["Sample_characteristics_ch1"][1] == "Age: embryos in slow phase of cellularisation"
            assert record.entity_attributes["Sample_scan_protocol"] == "GeneChips were scanned using the Hewlett-Packard GeneArray Scanner G2500A."
            assert record.entity_attributes["Sample_hyb_protocol"] == "Following fragmentation, 10 microg of cRNA were hybridized for 16 hr at 45C on GeneChip Drosophila Genome Array. GeneChips were washed and stained in the Affymetrix Fluidics Station 400."
            assert record.entity_attributes["Sample_extract_protocol_ch1"] == "Trizol extraction of total RNA was performed according to the manufacturer's instructions."
            assert record.entity_attributes["Sample_source_name_ch1"] == "Drosophila embryos in slow phase of cellularisation"
            assert record.entity_attributes["Sample_table_begin"] == ""
            assert record.entity_attributes["Sample_molecule_ch1"] == "total RNA"
            assert record.entity_attributes["Sample_label_protocol_ch1"] == "Biotinylated cRNA were prepared according to the standard Affymetrix protocol from 6 microg total RNA (Expression Analysis Technical Manual, 2001, Affymetrix)."
            assert record.entity_attributes["Sample_data_processing"] == "The data were analyzed with Microarray Suite version 5.0 (MAS 5.0) using Affymetrix default analysis settings and global scaling as normalization method. The trimmed mean target intensity of each array was arbitrarily set to 100."
            assert record.entity_attributes["Sample_treatment_protocol_ch1"] == "Embryos were dechorionated with 50% bleach, put on a cover slip and covered with Halocarbon oil 27 (Sigma). Embryos of the appropriate stage were manually selected under the dissecting scope. Selected embryos were transferred to a basket, rinsed with PBS with 0,7% NaCl, 0,04% triton-X100 and placed on ice in the Trizol solution (GibcoBRL)."
            assert record.entity_attributes["Sample_title"] == "embryo at T1, biological rep1"
            assert record.entity_attributes["Sample_supplementary_file"] == "Drosophila_T1-1.CEL"
            assert record.entity_attributes["Sample_platform_id"] == "GPL72"
            assert len(record.col_defs) == 4
            assert record.col_defs["DETECTION P-VALUE"] == "'detection p-value', p-value that indicates the significance level of the detection call"
            assert record.col_defs["ID_REF"] == ""
            assert record.col_defs["VALUE"] == "MAS5-calculated Signal intensity"
            assert record.col_defs["ABS_CALL"] == "the call in an absolute analysis that indicates if the transcript was present (P), absent (A), marginal (M), or no call (NC)"
            assert len(record.table_rows) == 22
            assert len(record.table_rows[0]) == 4
            assert record.table_rows[0][0] == "ID_REF"
            assert record.table_rows[0][1] == "VALUE"
            assert record.table_rows[0][2] == "ABS_CALL"
            assert record.table_rows[0][3] == "DETECTION P-VALUE"
            assert len(record.table_rows[1]) == 4
            assert record.table_rows[1][0] == "141200_at"
            assert record.table_rows[1][1] == "20.8"
            assert record.table_rows[1][2] == "A"
            assert record.table_rows[1][3] == "0.801637"
            assert len(record.table_rows[2]) == 4
            assert record.table_rows[2][0] == "141201_at"
            assert record.table_rows[2][1] == "85.8"
            assert record.table_rows[2][2] == "A"
            assert record.table_rows[2][3] == "0.48748"
            assert len(record.table_rows[3]) == 4
            assert record.table_rows[3][0] == "141202_at"
            assert record.table_rows[3][1] == "704.8"
            assert record.table_rows[3][2] == "P"
            assert record.table_rows[3][3] == "0.000613"
            assert len(record.table_rows[4]) == 4
            assert record.table_rows[4][0] == "141203_at"
            assert record.table_rows[4][1] == "1036.6"
            assert record.table_rows[4][2] == "P"
            assert record.table_rows[4][3] == "0.000944"
            assert len(record.table_rows[5]) == 4
            assert record.table_rows[5][0] == "141204_at"
            assert record.table_rows[5][1] == "700.3"
            assert record.table_rows[5][2] == "P"
            assert record.table_rows[5][3] == "0.000491"
            assert len(record.table_rows[6]) == 4
            assert record.table_rows[6][0] == "141205_at"
            assert record.table_rows[6][1] == "462.4"
            assert record.table_rows[6][2] == "P"
            assert record.table_rows[6][3] == "0.003159"
            assert len(record.table_rows[7]) == 4
            assert record.table_rows[7][0] == "141206_at"
            assert record.table_rows[7][1] == "1301.9"
            assert record.table_rows[7][2] == "P"
            assert record.table_rows[7][3] == "0.000491"
            assert len(record.table_rows[8]) == 4
            assert record.table_rows[8][0] == "141207_at"
            assert record.table_rows[8][1] == "454.8"
            assert record.table_rows[8][2] == "P"
            assert record.table_rows[8][3] == "0.000944"
            assert len(record.table_rows[9]) == 4
            assert record.table_rows[9][0] == "141208_at"
            assert record.table_rows[9][1] == "438.6"
            assert record.table_rows[9][2] == "P"
            assert record.table_rows[9][3] == "0.000944"
            assert len(record.table_rows[10]) == 4
            assert record.table_rows[10][0] == "141209_at"
            assert record.table_rows[10][1] == "264.4"
            assert record.table_rows[10][2] == "P"
            assert record.table_rows[10][3] == "0.004591"
            assert len(record.table_rows[11]) == 4
            assert record.table_rows[11][0] == "141210_at"
            assert record.table_rows[11][1] == "65.6"
            assert record.table_rows[11][2] == "A"
            assert record.table_rows[11][3] == "0.150145"
            assert len(record.table_rows[12]) == 4
            assert record.table_rows[12][0] == "141211_at"
            assert record.table_rows[12][1] == "72.2"
            assert record.table_rows[12][2] == "A"
            assert record.table_rows[12][3] == "0.070073"
            assert len(record.table_rows[13]) == 4
            assert record.table_rows[13][0] == "141212_at"
            assert record.table_rows[13][1] == "1200"
            assert record.table_rows[13][2] == "P"
            assert record.table_rows[13][3] == "0.000491"
            assert len(record.table_rows[14]) == 4
            assert record.table_rows[14][0] == "141213_at"
            assert record.table_rows[14][1] == "13.7"
            assert record.table_rows[14][2] == "A"
            assert record.table_rows[14][3] == "0.635055"
            assert len(record.table_rows[15]) == 4
            assert record.table_rows[15][0] == "142121_at"
            assert record.table_rows[15][1] == "133.7"
            assert record.table_rows[15][2] == "A"
            assert record.table_rows[15][3] == "0.889551"
            assert len(record.table_rows[16]) == 4
            assert record.table_rows[16][0] == "142122_at"
            assert record.table_rows[16][1] == "275.3"
            assert record.table_rows[16][2] == "A"
            assert record.table_rows[16][3] == "0.611218"
            assert len(record.table_rows[17]) == 4
            assert record.table_rows[17][0] == "142123_at"
            assert record.table_rows[17][1] == "307.6"
            assert record.table_rows[17][2] == "A"
            assert record.table_rows[17][3] == "0.611218"
            assert len(record.table_rows[18]) == 4
            assert record.table_rows[18][0] == "142124_at"
            assert record.table_rows[18][1] == "132.6"
            assert record.table_rows[18][2] == "A"
            assert record.table_rows[18][3] == "0.437646"
            assert len(record.table_rows[19]) == 4
            assert record.table_rows[19][0] == "142125_at"
            assert record.table_rows[19][1] == "195.8"
            assert record.table_rows[19][2] == "A"
            assert record.table_rows[19][3] == "0.110449"
            assert len(record.table_rows[20]) == 4
            assert record.table_rows[20][0] == "142126_at"
            assert record.table_rows[20][1] == "174.1"
            assert record.table_rows[20][2] == "A"
            assert record.table_rows[20][3] == "0.681117"
            assert len(record.table_rows[21]) == 4
            assert record.table_rows[21][0] == "142127_at"
            assert record.table_rows[21][1] == "316.3"
            assert record.table_rows[21][2] == "A"
            assert record.table_rows[21][3] == "0.65838"
            record = next(records)
            assert record.entity_type == "SERIES"
            assert record.entity_id == "Dros_embryo_timecourse"
            assert len(record.entity_attributes) == 6
            assert len(record.entity_attributes["Series_sample_id"]) == 3
            assert record.entity_attributes["Series_sample_id"][0] == "Drosophila_T0-1"
            assert record.entity_attributes["Series_sample_id"][1] == "Drosophila_T0-2"
            assert record.entity_attributes["Series_sample_id"][2] == "Drosophila_T1-1"
            assert len(record.entity_attributes["Series_contributor"]) == 5
            assert record.entity_attributes["Series_contributor"][0] == "Jane,Doe"
            assert record.entity_attributes["Series_contributor"][1] == "John,A,Smith"
            assert record.entity_attributes["Series_contributor"][2] == "Hans,van Elton"
            assert record.entity_attributes["Series_contributor"][3] == "John,Smithers Jr"
            assert record.entity_attributes["Series_contributor"][4] == "Jie,D,Chen"
            assert len(record.entity_attributes["Series_summary"]) == 2
            assert record.entity_attributes["Series_summary"][0] == "Morphogenesis of epithelial tissues relies on the precise developmental control of cell polarity and architecture. In the early Drosophila embryo, the primary epithelium forms during cellularisation, following a tightly controlled genetic programme where specific sets of genes are up-regulated. Some of them, for instance, control membrane invagination between the nuclei anchored at the apical surface of the syncytium."
            assert record.entity_attributes["Series_summary"][1] == "We used microarrays to detail the global programme of gene expression underlying cellularisation and identified distinct classes of up-regulated genes during this process."
            assert record.entity_attributes["Series_type"] == "time course"
            assert record.entity_attributes["Series_title"] == "Expression data from early Drosophila embryo"
            assert record.entity_attributes["Series_overall_design"] == "Drosophila embryos were selected at successive stages of early development for RNA extraction and hybridization on Affymetrix microarrays. We sought to obtain homogeneous populations of embryos at each developmental stage in order to increase the temporal resolution of expression profiles. To that end, we hand-selected embryos according to morphological criteria at five time-points: before pole cell formation, i.e. before zygotic transcription (T0), during the slow phase (T1) and the fast phase (T2) of cellularisation and at the beginning (T3) and the end (T4) of gastrulation."
            assert len(record.col_defs) == 0
            assert len(record.table_rows) == 0

    def test_GSE16(self):
        path = "Geo/GSE16.txt"
        with open(path, encoding="latin") as handle:
            records = Geo.parse(handle)
            record = next(records)
            assert record.entity_type == "SAMPLE"
            assert record.entity_id == "GSM804"
            assert len(record.entity_attributes) == 18
            assert record.entity_attributes["Sample_pubmed_id"] == "11687795"
            assert record.entity_attributes["Sample_submitter_institute"] == "University of California San Francisco"
            assert len(record.entity_attributes["Sample_author"]) == 19
            assert record.entity_attributes["Sample_author"][0] == "Antoine,M,Snijders"
            assert record.entity_attributes["Sample_author"][1] == "Norma,,Nowak"
            assert record.entity_attributes["Sample_author"][2] == "Richard,,Segraves"
            assert record.entity_attributes["Sample_author"][3] == "Stephanie,,Blackwood"
            assert record.entity_attributes["Sample_author"][4] == "Nils,,Brown"
            assert record.entity_attributes["Sample_author"][5] == "Jeffery,,Conroy"
            assert record.entity_attributes["Sample_author"][6] == "Greg,,Hamilton"
            assert record.entity_attributes["Sample_author"][7] == "Anna,K,Hindle"
            assert record.entity_attributes["Sample_author"][8] == "Bing,,Huey"
            assert record.entity_attributes["Sample_author"][9] == "Karen,,Kimura"
            assert record.entity_attributes["Sample_author"][10] == "Sindy,,Law"
            assert record.entity_attributes["Sample_author"][11] == "Ken,,Myambo"
            assert record.entity_attributes["Sample_author"][12] == "Joel,,Palmer"
            assert record.entity_attributes["Sample_author"][13] == "Bauke,,Ylstra"
            assert record.entity_attributes["Sample_author"][14] == "Jingzhu,P,Yue"
            assert record.entity_attributes["Sample_author"][15] == "Joe,W,Gray"
            assert record.entity_attributes["Sample_author"][16] == "Ajay,N,Jain"
            assert record.entity_attributes["Sample_author"][17] == "Daniel,,Pinkel"
            assert record.entity_attributes["Sample_author"][18] == "Donna,G,Albertson"
            assert record.entity_attributes["Sample_submitter_phone"] == "415 502-8463"
            assert record.entity_attributes["Sample_submitter_department"] == "Comprehensive Cancer Center"
            assert len(record.entity_attributes["Sample_description"]) == 4
            assert record.entity_attributes["Sample_description"][0] == 'Coriell Cell Repositories cell line <a href="http://locus.umdnj.edu/nigms/nigms_cgi/display.cgi?GM05296">GM05296</a>.'
            assert record.entity_attributes["Sample_description"][1] == "Fibroblast cell line derived from a 1 month old female with multiple congenital malformations, dysmorphic features, intrauterine growth retardation, heart murmur, cleft palate, equinovarus deformity, microcephaly, coloboma of right iris, clinodactyly, reduced RBC catalase activity, and 1 copy of catalase gene."
            assert record.entity_attributes["Sample_description"][2] == "Chromosome abnormalities are present."
            assert record.entity_attributes["Sample_description"][3] == "Karyotype is 46,XX,-11,+der(11)inv ins(11;10)(11pter> 11p13::10q21>10q24::11p13>11qter)mat"
            assert record.entity_attributes["Sample_target_source2"] == "normal male reference genomic DNA"
            assert record.entity_attributes["Sample_target_source1"] == "Cell line GM05296"
            assert record.entity_attributes["Sample_submitter_name"] == "Donna,G,Albertson"
            assert record.entity_attributes["Sample_platform_id"] == "GPL28"
            assert record.entity_attributes["Sample_type"] == "dual channel genomic"
            assert record.entity_attributes["Sample_status"] == "Public on Feb 12 2002"
            assert record.entity_attributes["Sample_submitter_email"] == "albertson@cc.ucsf.edu"
            assert record.entity_attributes["Sample_title"] == "CGH_Albertson_GM05296-001218"
            assert record.entity_attributes["Sample_organism"] == "Homo sapiens"
            assert record.entity_attributes["Sample_series_id"] == "GSE16"
            assert record.entity_attributes["Sample_submission_date"] == "Jan 17 2002"
            assert record.entity_attributes["Sample_submitter_city"] == "San Francisco,CA,94143,USA"
            assert len(record.col_defs) == 5
            assert record.col_defs["NO_REPLICATES"] == "Number of replicate spot measurements"
            assert record.col_defs["LOG2STDDEV"] == "Standard deviation of VALUE"
            assert record.col_defs["ID_REF"] == "Unique row identifier, genome position order"
            assert record.col_defs["VALUE"] == "aka LOG2RATIO, mean of log base 2 of LINEAR_RATIO"
            assert record.col_defs["LINEAR_RATIO"] == "Mean of replicate Cy3/Cy5 ratios"
            assert len(record.table_rows) == 21
            assert len(record.table_rows[0]) == 5
            assert record.table_rows[0][0] == "ID_REF"
            assert record.table_rows[0][1] == "VALUE"
            assert record.table_rows[0][2] == "LINEAR_RATIO"
            assert record.table_rows[0][3] == "LOG2STDDEV"
            assert record.table_rows[0][4] == "NO_REPLICATES"
            assert len(record.table_rows[1]) == 5
            assert record.table_rows[1][0] == "1"
            assert record.table_rows[1][1] == ""
            assert record.table_rows[1][2] == "1.047765"
            assert record.table_rows[1][3] == "0.011853"
            assert record.table_rows[1][4] == "3"
            assert len(record.table_rows[2]) == 5
            assert record.table_rows[2][0] == "2"
            assert record.table_rows[2][1] == ""
            assert record.table_rows[2][2] == ""
            assert record.table_rows[2][3] == ""
            assert record.table_rows[2][4] == "0"
            assert len(record.table_rows[3]) == 5
            assert record.table_rows[3][0] == "3"
            assert record.table_rows[3][1] == "0.008824"
            assert record.table_rows[3][2] == "1.006135"
            assert record.table_rows[3][3] == "0.00143"
            assert record.table_rows[3][4] == "3"
            assert len(record.table_rows[4]) == 5
            assert record.table_rows[4][0] == "4"
            assert record.table_rows[4][1] == "-0.000894"
            assert record.table_rows[4][2] == "0.99938"
            assert record.table_rows[4][3] == "0.001454"
            assert record.table_rows[4][4] == "3"
            assert len(record.table_rows[5]) == 5
            assert record.table_rows[5][0] == "5"
            assert record.table_rows[5][1] == "0.075875"
            assert record.table_rows[5][2] == "1.054"
            assert record.table_rows[5][3] == "0.003077"
            assert record.table_rows[5][4] == "3"
            assert len(record.table_rows[6]) == 5
            assert record.table_rows[6][0] == "6"
            assert record.table_rows[6][1] == "0.017303"
            assert record.table_rows[6][2] == "1.012066"
            assert record.table_rows[6][3] == "0.005876"
            assert record.table_rows[6][4] == "2"
            assert len(record.table_rows[7]) == 5
            assert record.table_rows[7][0] == "7"
            assert record.table_rows[7][1] == "-0.006766"
            assert record.table_rows[7][2] == "0.995321"
            assert record.table_rows[7][3] == "0.013881"
            assert record.table_rows[7][4] == "3"
            assert len(record.table_rows[8]) == 5
            assert record.table_rows[8][0] == "8"
            assert record.table_rows[8][1] == "0.020755"
            assert record.table_rows[8][2] == "1.014491"
            assert record.table_rows[8][3] == "0.005506"
            assert record.table_rows[8][4] == "3"
            assert len(record.table_rows[9]) == 5
            assert record.table_rows[9][0] == "9"
            assert record.table_rows[9][1] == "-0.094938"
            assert record.table_rows[9][2] == "0.936313"
            assert record.table_rows[9][3] == "0.012662"
            assert record.table_rows[9][4] == "3"
            assert len(record.table_rows[10]) == 5
            assert record.table_rows[10][0] == "10"
            assert record.table_rows[10][1] == "-0.054527"
            assert record.table_rows[10][2] == "0.96291"
            assert record.table_rows[10][3] == "0.01073"
            assert record.table_rows[10][4] == "3"
            assert len(record.table_rows[11]) == 5
            assert record.table_rows[11][0] == "11"
            assert record.table_rows[11][1] == "-0.025057"
            assert record.table_rows[11][2] == "0.982782"
            assert record.table_rows[11][3] == "0.003855"
            assert record.table_rows[11][4] == "3"
            assert len(record.table_rows[12]) == 5
            assert record.table_rows[12][0] == "12"
            assert record.table_rows[12][1] == ""
            assert record.table_rows[12][2] == ""
            assert record.table_rows[12][3] == ""
            assert record.table_rows[12][4] == "0"
            assert len(record.table_rows[13]) == 5
            assert record.table_rows[13][0] == "13"
            assert record.table_rows[13][1] == "0.108454"
            assert record.table_rows[13][2] == "1.078072"
            assert record.table_rows[13][3] == "0.005196"
            assert record.table_rows[13][4] == "3"
            assert len(record.table_rows[14]) == 5
            assert record.table_rows[14][0] == "14"
            assert record.table_rows[14][1] == "0.078633"
            assert record.table_rows[14][2] == "1.056017"
            assert record.table_rows[14][3] == "0.009165"
            assert record.table_rows[14][4] == "3"
            assert len(record.table_rows[15]) == 5
            assert record.table_rows[15][0] == "15"
            assert record.table_rows[15][1] == "0.098571"
            assert record.table_rows[15][2] == "1.070712"
            assert record.table_rows[15][3] == "0.007834"
            assert record.table_rows[15][4] == "3"
            assert len(record.table_rows[16]) == 5
            assert record.table_rows[16][0] == "16"
            assert record.table_rows[16][1] == "0.044048"
            assert record.table_rows[16][2] == "1.031003"
            assert record.table_rows[16][3] == "0.013651"
            assert record.table_rows[16][4] == "3"
            assert len(record.table_rows[17]) == 5
            assert record.table_rows[17][0] == "17"
            assert record.table_rows[17][1] == "0.018039"
            assert record.table_rows[17][2] == "1.012582"
            assert record.table_rows[17][3] == "0.005471"
            assert record.table_rows[17][4] == "3"
            assert len(record.table_rows[18]) == 5
            assert record.table_rows[18][0] == "18"
            assert record.table_rows[18][1] == "-0.088807"
            assert record.table_rows[18][2] == "0.9403"
            assert record.table_rows[18][3] == "0.010571"
            assert record.table_rows[18][4] == "3"
            assert len(record.table_rows[19]) == 5
            assert record.table_rows[19][0] == "19"
            assert record.table_rows[19][1] == "0.016349"
            assert record.table_rows[19][2] == "1.011397"
            assert record.table_rows[19][3] == "0.007113"
            assert record.table_rows[19][4] == "3"
            assert len(record.table_rows[20]) == 5
            assert record.table_rows[20][0] == "20"
            assert record.table_rows[20][1] == "0.030977"
            assert record.table_rows[20][2] == "1.021704"
            assert record.table_rows[20][3] == "0.016798"
            assert record.table_rows[20][4] == "3"

    def test_soft_ex_platform(self):
        path = "Geo/soft_ex_platform.txt"
        with open(path, encoding="latin") as handle:
            records = Geo.parse(handle)
            record = next(records)
            assert record.entity_type == "PLATFORM"
            assert record.entity_id == "Murine 15K long oligo array version 2.0"
            assert len(record.entity_attributes) == 12
            assert record.entity_attributes["Platform_title"] == "Murine 15K long oligo array version 2.0"
            assert record.entity_attributes["Platform_web_link"] == "http://www.microarray.protocols.html"
            assert record.entity_attributes["platform_table_end"] == ""
            assert record.entity_attributes["Platform_support"] == "glass"
            assert record.entity_attributes["Platform_manufacturer"] == "Un. London microarray facility"
            assert record.entity_attributes["Platform_coating"] == "polysine"
            assert record.entity_attributes["Platform_technology"] == "spotted oligonucleotide"
            assert record.entity_attributes["platform_table_begin"] == ""
            assert len(record.entity_attributes["Platform_manufacture_protocol"]) == 12
            assert record.entity_attributes["Platform_manufacture_protocol"][0] == "1.  Oligos are arrayed in Greiner 384-well flat-bottom plates. Each well contains 600 pmol of 70-mer oligo."
            assert record.entity_attributes["Platform_manufacture_protocol"][1] == "2. Resuspend oligos in water to 20 uM and rearray 5 \xb5L into 384-well, Genetix polystyrene V-bottom plates (cat# X6004)."
            assert record.entity_attributes["Platform_manufacture_protocol"][2] == "3. Allow Genetix plates to dry through passive water evaporation in a protected environment (e.g., chemical hood)."
            assert record.entity_attributes["Platform_manufacture_protocol"][3] == "4. Before printing, add 5 \xb5L of 1X Printing Buffer to each well. This can be done the night before a print run is started."
            assert record.entity_attributes["Platform_manufacture_protocol"][4] == "5. Seal plates with Corning seals."
            assert record.entity_attributes["Platform_manufacture_protocol"][5] == "6. Incubate at 37\xb0C for 30 minutes to aid resuspension of DNA."
            assert record.entity_attributes["Platform_manufacture_protocol"][6] == "7. Shake plates near maximum rotational speed on flat-bed shaker for 1 minute."
            assert record.entity_attributes["Platform_manufacture_protocol"][7] == "8. Centrifuge plates at 2000 rpm for 3 minutes."
            assert record.entity_attributes["Platform_manufacture_protocol"][8] == "9. Remove seals and cover with plate lids. Place in appropriate location of plate cassette. This should be done with first plates just before print run is started to minimize evaporation time before printing. For second and third cassettes, wait until 30 minutes before next cassette is needed to begin centrifugation."
            assert record.entity_attributes["Platform_manufacture_protocol"][9] == "10. Make sure plates rest behind both holding clips in the cassettes. Push plates back into the cassettes as far as they will go, putting them in the proper position for the server arm."
            assert record.entity_attributes["Platform_manufacture_protocol"][10] == "11. After the print run is completed, allow plates to dry through passive evaporation in a protected environment."
            assert record.entity_attributes["Platform_manufacture_protocol"][11] == "12. For each subsequent preparation of these plates for a print run, add water to the wells instead of sodium phosphate buffer. The amount of water should be decreased by 0.25 \xb5L per print run, as this is the amount drawn up by the pin capillary during each dip."
            assert record.entity_attributes["Platform_organism"] == "Mus musculus"
            assert len(record.entity_attributes["Platform_contributor"]) == 5
            assert record.entity_attributes["Platform_contributor"][0] == "Jane,Doe"
            assert record.entity_attributes["Platform_contributor"][1] == "John,A,Smith"
            assert record.entity_attributes["Platform_contributor"][2] == "Hans,van Elton"
            assert record.entity_attributes["Platform_contributor"][3] == "John,Smithers Jr"
            assert record.entity_attributes["Platform_contributor"][4] == "Jie,D,Chen"
            assert record.entity_attributes["Platform_distribution"] == "non-commercial"
            assert len(record.col_defs) == 6
            assert record.col_defs["Gene_Desc"] == "Gene description"
            assert record.col_defs["SEQUENCE"] == "Probe sequence information"
            assert record.col_defs["Gene_Sym"] == "Gene symbols"
            assert record.col_defs["GB_ACC"] == 'GenBank accession number of sequence used to design oligonucleotide probe   LINK_PRE:"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Search&db=Nucleotide&term="'
            assert record.col_defs["SPOT_ID"] == "alternative identifier"
            assert record.col_defs["ID"] == ""
            assert len(record.table_rows) == 21
            assert len(record.table_rows[0]) == 6
            assert record.table_rows[0][0] == "ID"
            assert record.table_rows[0][1] == "GB_ACC"
            assert record.table_rows[0][2] == "Gene_Desc"
            assert record.table_rows[0][3] == "Gene_Sym"
            assert record.table_rows[0][4] == "SPOT_ID"
            assert record.table_rows[0][5] == "SEQUENCE"
            assert len(record.table_rows[1]) == 6
            assert record.table_rows[1][0] == "1"
            assert record.table_rows[1][1] == "U02079"
            assert record.table_rows[1][2] == "nuclear factor of activated T-cells, cytoplasmic 2"
            assert record.table_rows[1][3] == "Nfatc2"
            assert record.table_rows[1][4] == ""
            assert record.table_rows[1][5] == "ACCTGGATGACGCAGCCACTTCAGAAAGCTGGGTTGGGACAGAAAGGTATATAGAGAGAAAATTTTGGAA"
            assert len(record.table_rows[2]) == 6
            assert record.table_rows[2][0] == "2"
            assert record.table_rows[2][1] == "NM_008154"
            assert record.table_rows[2][2] == "G-protein coupled receptor 3"
            assert record.table_rows[2][3] == "Gpr3"
            assert record.table_rows[2][4] == ""
            assert record.table_rows[2][5] == "CTGTACAATGCTCTCACTTACTACTCAGAGACAACGGTAACTCGGACTTATGTGATGCTGGCCTTGGTGT"
            assert len(record.table_rows[3]) == 6
            assert record.table_rows[3][0] == "3"
            assert record.table_rows[3][1] == "AK015719"
            assert record.table_rows[3][2] == "tropomodulin 2"
            assert record.table_rows[3][3] == "Tmod2"
            assert record.table_rows[3][4] == ""
            assert record.table_rows[3][5] == "CACCAGGCTCAGTGCCTAGTATCGGCTTCACCTAGTGTGGTTACTCAGGGCACGCAGAGCTACAGAACAC"
            assert len(record.table_rows[4]) == 6
            assert record.table_rows[4][0] == "4"
            assert record.table_rows[4][1] == "AK003367"
            assert record.table_rows[4][2] == "mitochondrial ribosomal protein L15"
            assert record.table_rows[4][3] == "Mrpl15"
            assert record.table_rows[4][4] == ""
            assert record.table_rows[4][5] == "CAAGAAGTCTAGAAATTCTGTGCAAGCCTATTCCATTCTTTCTGCGGGGACAACCAATTCCGAAAAGAAT"
            assert len(record.table_rows[5]) == 6
            assert record.table_rows[5][0] == "5"
            assert record.table_rows[5][1] == "BC003333"
            assert record.table_rows[5][2] == "RIKEN cDNA 0610033I05 gene"
            assert record.table_rows[5][3] == "0610033I05Rik"
            assert record.table_rows[5][4] == ""
            assert record.table_rows[5][5] == "AGAACTGGGTGGCAGATATCCTAGAGTTTTGACCAACGTTCACAGCACACATATTGATCTTATAGGACCT"
            assert len(record.table_rows[6]) == 6
            assert record.table_rows[6][0] == "6"
            assert record.table_rows[6][1] == "NM_008462"
            assert record.table_rows[6][2] == "killer cell lectin-like receptor, subfamily A, member 2"
            assert record.table_rows[6][3] == "Klra2"
            assert record.table_rows[6][4] == ""
            assert record.table_rows[6][5] == "TGAATTGAAGTTCCTTAAATCCCAACTTCAAAGAAACACATACTGGATTTCACTGACACATCATAAAAGC"
            assert len(record.table_rows[7]) == 6
            assert record.table_rows[7][0] == "7"
            assert record.table_rows[7][1] == "NM_008029"
            assert record.table_rows[7][2] == "FMS-like tyrosine kinase 4"
            assert record.table_rows[7][3] == "Flt4"
            assert record.table_rows[7][4] == ""
            assert record.table_rows[7][5] == "GAGGTGCTGTGGGATGACCGCCGGGGCATGCGGGTGCCCACTCAACTGTTGCGCGATGCCCTGTACCTGC"
            assert len(record.table_rows[8]) == 6
            assert record.table_rows[8][0] == "8"
            assert record.table_rows[8][1] == "NM_054088"
            assert record.table_rows[8][2] == "adiponutrin"
            assert record.table_rows[8][3] == "Adpn"
            assert record.table_rows[8][4] == ""
            assert record.table_rows[8][5] == "GTCTGAGTTCCATTCCAAAGACGAAGTCGTGGATGCCCTGGTGTGTTCCTGCTTCATTCCCCTCTTCTCT"
            assert len(record.table_rows[9]) == 6
            assert record.table_rows[9][0] == "9"
            assert record.table_rows[9][1] == "NM_009750"
            assert record.table_rows[9][2] == "nerve growth factor receptor (TNFRSF16) associated protein 1"
            assert record.table_rows[9][3] == "Ngfrap1"
            assert record.table_rows[9][4] == ""
            assert record.table_rows[9][5] == "TACAGCTGAGAAATTGTCTACGCATCCTTATGGGGGAGCTGTCTAACCACCACGATCACCATGATGAATT"
            assert len(record.table_rows[10]) == 6
            assert record.table_rows[10][0] == "10"
            assert record.table_rows[10][1] == "AB045323"
            assert record.table_rows[10][2] == "DNA segment, Chr 8, ERATO Doi 594, expressed"
            assert record.table_rows[10][3] == "D8Ertd594e"
            assert record.table_rows[10][4] == ""
            assert record.table_rows[10][5] == "GATTCAGACTCGGGAGGAGCATCCCAACCTCTCCTTGAGGATAAAGGCCTGAGCGATTGCCCTGGGGAGC"
            assert len(record.table_rows[11]) == 6
            assert record.table_rows[11][0] == "11"
            assert record.table_rows[11][1] == "AK005789"
            assert record.table_rows[11][2] == "dynein, cytoplasmic, light chain 2B"
            assert record.table_rows[11][3] == "Dncl2b"
            assert record.table_rows[11][4] == ""
            assert record.table_rows[11][5] == "TGCAGAAGGCATTCCAATCCGAACAACCCTGGACAACTCCACAACGGTTCAGTATGCGGGTCTTCTCCAC"
            assert len(record.table_rows[12]) == 6
            assert record.table_rows[12][0] == "12"
            assert record.table_rows[12][1] == "NM_010517"
            assert record.table_rows[12][2] == "insulin-like growth factor binding protein 4"
            assert record.table_rows[12][3] == "Igfbp4"
            assert record.table_rows[12][4] == ""
            assert record.table_rows[12][5] == "GGAGAAGCTGGCGCGCTGCCGCCCCCCCGTGGGTTGCGAGGAGTTGGTGCGGGAGCCAGGCTGCGGTTGT"
            assert len(record.table_rows[13]) == 6
            assert record.table_rows[13][0] == "13"
            assert record.table_rows[13][1] == "AK010722"
            assert record.table_rows[13][2] == "RIKEN cDNA 2410075D05 gene"
            assert record.table_rows[13][3] == "2410075D05Rik"
            assert record.table_rows[13][4] == ""
            assert record.table_rows[13][5] == "GGAGCATCTGGAGTTCCGCTTACCGGAAATAAAGTCTTTACTATCGGTGATTGGAGGGCAGTTCACTAAC"
            assert len(record.table_rows[14]) == 6
            assert record.table_rows[14][0] == "14"
            assert record.table_rows[14][1] == "AK003755"
            assert record.table_rows[14][2] == "DNA segment, Chr 4, ERATO Doi 421, expressed"
            assert record.table_rows[14][3] == "D4Ertd421e"
            assert record.table_rows[14][4] == ""
            assert record.table_rows[14][5] == "AGCAAAGAGATCTCCCTCAGTGTGCCCATAGGTGGCGGTGCGAGCTTGCGGTTATTGGCCAGTGACTTGC"
            assert len(record.table_rows[15]) == 6
            assert record.table_rows[15][0] == "15"
            assert record.table_rows[15][1] == "BC003241"
            assert record.table_rows[15][2] == "cleavage stimulation factor, 3' pre-RNA, subunit 3"
            assert record.table_rows[15][3] == "Cstf3"
            assert record.table_rows[15][4] == ""
            assert record.table_rows[15][5] == "AAATTAGAAGAAAATCCATATGACCTTGATGCTTGGAGCATTCTCATTCGAGAGGCACAGAATCAACCTA"
            assert len(record.table_rows[16]) == 6
            assert record.table_rows[16][0] == "16"
            assert record.table_rows[16][1] == "AK004937"
            assert record.table_rows[16][2] == "RIKEN cDNA 1300007O09 gene"
            assert record.table_rows[16][3] == "1300007O09Rik"
            assert record.table_rows[16][4] == ""
            assert record.table_rows[16][5] == "CAGACACAAACCCTAGGTTGTATTGTAGACCGGAGTTTAAGCAGGCACTACCTGTCTGTCTTTTCTTCAT"
            assert len(record.table_rows[17]) == 6
            assert record.table_rows[17][0] == "17"
            assert record.table_rows[17][1] == "AK004524"
            assert record.table_rows[17][2] == "unnamed protein product; hypothetical SOCS domain"
            assert record.table_rows[17][3] == ""
            assert record.table_rows[17][4] == ""
            assert record.table_rows[17][5] == "CGGAGCCCTGCGCGCCCAGAGCCCCCTCCCACCCGCTTCCACCAAGTGCATGGAGCCAACATCCGCATGG"
            assert len(record.table_rows[18]) == 6
            assert record.table_rows[18][0] == "18"
            assert record.table_rows[18][1] == "NM_025999"
            assert record.table_rows[18][2] == "RIKEN cDNA 2610110L04 gene"
            assert record.table_rows[18][3] == "2610110L04Rik"
            assert record.table_rows[18][4] == ""
            assert record.table_rows[18][5] == "TGCATTGATAAATGGAGTGATCGACACAGGAACTGCCCCATTTGTCGCCTACAGATGACTGGAGCAAATG"
            assert len(record.table_rows[19]) == 6
            assert record.table_rows[19][0] == "19"
            assert record.table_rows[19][1] == ""
            assert record.table_rows[19][2] == ""
            assert record.table_rows[19][3] == ""
            assert record.table_rows[19][4] == "-- CONTROL"
            assert record.table_rows[19][5] == ""
            assert len(record.table_rows[20]) == 6
            assert record.table_rows[20][0] == "20"
            assert record.table_rows[20][1] == "NM_023120"
            assert record.table_rows[20][2] == "guanine nucleotide binding protein (G protein), beta polypeptide 1-like"
            assert record.table_rows[20][3] == "Gnb1l"
            assert record.table_rows[20][4] == ""
            assert record.table_rows[20][5] == "ACCGCCTGGTCCCAGATTTGTCCTCCGAGGCACACAGTCGGCTGTGAACACGCTCCATTTCTGCCCACCA"

    def test_GSM700(self):
        path = "Geo/GSM700.txt"
        with open(path, encoding="latin") as handle:
            records = Geo.parse(handle)
            record = next(records)
            assert record.entity_type == "SAMPLE"
            assert record.entity_id == "GSM700"
            assert len(record.entity_attributes) == 20
            assert record.entity_attributes["Sample_submitter_institute"] == "National Cancer Institute"
            assert record.entity_attributes["Sample_submitter_department"] == "Cancer Genome Anatomy Project"
            assert record.entity_attributes["Sample_submitter_web_link"] == "http://cgap.nci.nih.gov/"
            assert len(record.entity_attributes["Sample_description"]) == 14
            assert record.entity_attributes["Sample_description"][0] == "This library represents a Cancer Genome Anatomy Project library, which was either produced through CGAP funding, or donated to CGAP."
            assert record.entity_attributes["Sample_description"][1] == "The Cancer Genome Anatomy Project (CGAP: http://cgap.nci.nih.gov) is an interdisciplinary program established and administered by the National Cancer Institute (NCI: http://www.nci.nih.gov) to generate the information and technological tools needed to decipher the molecular anatomy of the cancer cell."
            assert record.entity_attributes["Sample_description"][2] == "Cell line grown under 1.5% oxygen conditions for 24 hours prior to harvesting in zinc option media with 10% RBS and harvested at passage 102. Library constructed in the laboratory of G. Riggins, M.D., Ph.D. (Duke University)."
            assert record.entity_attributes["Sample_description"][3] == "Organ: brain"
            assert record.entity_attributes["Sample_description"][4] == "Tissue_type: glioblastoma multiforme"
            assert record.entity_attributes["Sample_description"][5] == "Cell_line: H247"
            assert record.entity_attributes["Sample_description"][6] == "Lab host: DH10B"
            assert record.entity_attributes["Sample_description"][7] == "Vector: pZErO-1"
            assert record.entity_attributes["Sample_description"][8] == "Vector type: plasmid"
            assert record.entity_attributes["Sample_description"][9] == "R. Site 1: Sph1"
            assert record.entity_attributes["Sample_description"][10] == "R. Site 2: Sph1"
            assert record.entity_attributes["Sample_description"][11] == "Library treatment: non-normalized"
            assert record.entity_attributes["Sample_description"][12] == "Tissue description: Brain, Duke glioblastoma multiforme cell line, H247, grown under 1.5% oxygen conditions  for 24 hours prior to harvesting."
            assert record.entity_attributes["Sample_description"][13] == "Tissue"
            assert len(record.entity_attributes["Sample_author"]) == 2
            assert record.entity_attributes["Sample_author"][0] == "Gregory,J,Riggins"
            assert record.entity_attributes["Sample_author"][1] == "Robert,L,Strausberg"
            assert record.entity_attributes["Sample_web_link"] == "http://cgap.nci.nih.gov"
            assert record.entity_attributes["Sample_submitter_phone"] == "301-496-1550"
            assert record.entity_attributes["Sample_series_id"] == "GSE14"
            assert record.entity_attributes["Sample_tag_count"] == "72031"
            assert record.entity_attributes["Sample_type"] == "sage"
            assert record.entity_attributes["Sample_submitter_name"] == "Robert,L,Strausberg"
            assert record.entity_attributes["Sample_platform_id"] == "GPL4"
            assert record.entity_attributes["Sample_submitter_city"] == "Bethesda,MD,20892,USA"
            assert record.entity_attributes["Sample_status"] == "Public on Nov 28 2001"
            assert record.entity_attributes["Sample_anchor"] == "NlaIII"
            assert record.entity_attributes["Sample_title"] == "SAGE_Duke_H247_Hypoxia"
            assert record.entity_attributes["Sample_organism"] == "Homo sapiens"
            assert record.entity_attributes["Sample_target_source"] == "Brain, glioblastoma multiforme, cell-line H247"
            assert record.entity_attributes["Sample_submission_date"] == "Nov 28 2001"
            assert record.entity_attributes["Sample_submitter_email"] == "cgapbs-r@mail.nih.gov"
            assert len(record.col_defs) == 3
            assert record.col_defs["COUNT"] == "Absolute tag count"
            assert record.col_defs["TPM"] == "Tags per million, or (1000000*COUNT)/(Total tags)"
            assert record.col_defs["TAG"] == 'Ten base SAGE tag, LINK_PRE:"http://www.ncbi.nlm.nih.gov/SAGE/SAGEtag.cgi?tag="'
            assert len(record.table_rows) == 21
            assert len(record.table_rows[0]) == 3
            assert record.table_rows[0][0] == "TAG"
            assert record.table_rows[0][1] == "COUNT"
            assert record.table_rows[0][2] == "TPM"
            assert len(record.table_rows[1]) == 3
            assert record.table_rows[1][0] == "TCCAAATCGA"
            assert record.table_rows[1][1] == "520"
            assert record.table_rows[1][2] == "7219.11"
            assert len(record.table_rows[2]) == 3
            assert record.table_rows[2][0] == "TACCATCAAT"
            assert record.table_rows[2][1] == "434"
            assert record.table_rows[2][2] == "6025.18"
            assert len(record.table_rows[3]) == 3
            assert record.table_rows[3][0] == "TTGGGGTTTC"
            assert record.table_rows[3][1] == "389"
            assert record.table_rows[3][2] == "5400.45"
            assert len(record.table_rows[4]) == 3
            assert record.table_rows[4][0] == "CCCATCGTCC"
            assert record.table_rows[4][1] == "367"
            assert record.table_rows[4][2] == "5095.03"
            assert len(record.table_rows[5]) == 3
            assert record.table_rows[5][0] == "GTGAAACCCC"
            assert record.table_rows[5][1] == "365"
            assert record.table_rows[5][2] == "5067.26"
            assert len(record.table_rows[6]) == 3
            assert record.table_rows[6][0] == "GGGGAAATCG"
            assert record.table_rows[6][1] == "357"
            assert record.table_rows[6][2] == "4956.2"
            assert len(record.table_rows[7]) == 3
            assert record.table_rows[7][0] == "CCTGTAATCC"
            assert record.table_rows[7][1] == "346"
            assert record.table_rows[7][2] == "4803.49"
            assert len(record.table_rows[8]) == 3
            assert record.table_rows[8][0] == "TGATTTCACT"
            assert record.table_rows[8][1] == "334"
            assert record.table_rows[8][2] == "4636.89"
            assert len(record.table_rows[9]) == 3
            assert record.table_rows[9][0] == "TGTGTTGAGA"
            assert record.table_rows[9][1] == "315"
            assert record.table_rows[9][2] == "4373.12"
            assert len(record.table_rows[10]) == 3
            assert record.table_rows[10][0] == "GCCCCCAATA"
            assert record.table_rows[10][1] == "303"
            assert record.table_rows[10][2] == "4206.52"
            assert len(record.table_rows[11]) == 3
            assert record.table_rows[11][0] == "CTAAGACTTC"
            assert record.table_rows[11][1] == "279"
            assert record.table_rows[11][2] == "3873.33"
            assert len(record.table_rows[12]) == 3
            assert record.table_rows[12][0] == "GCGACCGTCA"
            assert record.table_rows[12][1] == "276"
            assert record.table_rows[12][2] == "3831.68"
            assert len(record.table_rows[13]) == 3
            assert record.table_rows[13][0] == "TTGGTCCTCT"
            assert record.table_rows[13][1] == "276"
            assert record.table_rows[13][2] == "3831.68"
            assert len(record.table_rows[14]) == 3
            assert record.table_rows[14][0] == "CCTAGCTGGA"
            assert record.table_rows[14][1] == "268"
            assert record.table_rows[14][2] == "3720.62"
            assert len(record.table_rows[15]) == 3
            assert record.table_rows[15][0] == "GATGAGGAGA"
            assert record.table_rows[15][1] == "251"
            assert record.table_rows[15][2] == "3484.61"
            assert len(record.table_rows[16]) == 3
            assert record.table_rows[16][0] == "ACTTTTTCAA"
            assert record.table_rows[16][1] == "244"
            assert record.table_rows[16][2] == "3387.43"
            assert len(record.table_rows[17]) == 3
            assert record.table_rows[17][0] == "CCACTGCACT"
            assert record.table_rows[17][1] == "223"
            assert record.table_rows[17][2] == "3095.89"
            assert len(record.table_rows[18]) == 3
            assert record.table_rows[18][0] == "GTGTGTTTGT"
            assert record.table_rows[18][1] == "223"
            assert record.table_rows[18][2] == "3095.89"
            assert len(record.table_rows[19]) == 3
            assert record.table_rows[19][0] == "GAAATACAGT"
            assert record.table_rows[19][1] == "218"
            assert record.table_rows[19][2] == "3026.47"
            assert len(record.table_rows[20]) == 3
            assert record.table_rows[20][0] == "GCTTTATTTG"
            assert record.table_rows[20][1] == "218"
            assert record.table_rows[20][2] == "3026.47"

    def test_GSM645(self):
        path = "Geo/GSM645.txt"
        with open(path, encoding="latin") as handle:
            records = Geo.parse(handle)
            record = next(records)
            assert record.entity_type == "SAMPLE"
            assert record.entity_id == "GSM645"
            assert len(record.entity_attributes) == 17
            assert record.entity_attributes["Sample_submitter_institute"] == "Max von Pettenkofer Institut"
            assert record.entity_attributes["Sample_submitter_department"] == "Bacteriology"
            assert len(record.entity_attributes["Sample_author"]) == 4
            assert record.entity_attributes["Sample_author"][0] == "Reinhard,,Hoffmann"
            assert record.entity_attributes["Sample_author"][1] == "Thomas,,Seidl"
            assert record.entity_attributes["Sample_author"][2] == "Ton,,Rolink"
            assert record.entity_attributes["Sample_author"][3] == "Fritz,,Melchers"
            assert record.entity_attributes["Sample_submitter_phone"] == "+49-89-5160-5424"
            assert record.entity_attributes["Sample_series_id"] == "GSE13"
            assert record.entity_attributes["Sample_description"] == "B220+CD25+sIg- Large Pre BII cells sorted out of mouse bone marrow, sort no. 8"
            assert record.entity_attributes["Sample_type"] == "single channel"
            assert record.entity_attributes["Sample_submitter_name"] == "Reinhard,,Hoffmann"
            assert record.entity_attributes["Sample_platform_id"] == "GPL22"
            assert record.entity_attributes["Sample_submitter_city"] == "Munich,80336,Germany"
            assert record.entity_attributes["Sample_status"] == "Public on Dec 17 2001"
            assert record.entity_attributes["Sample_submitter_email"] == "r_hoffmann@m3401.mpk.med.uni-muenchen.de"
            assert record.entity_attributes["Sample_title"] == "Large Pre-BII cells 8b"
            assert record.entity_attributes["Sample_organism"] == "Mus musculus"
            assert record.entity_attributes["Sample_target_source"] == "Large Pre-BII cells"
            assert record.entity_attributes["Sample_submission_date"] == "Nov 27 2001"
            assert record.entity_attributes["Sample_submitter_address"] == "Pettenkoferstr. 9a"
            assert len(record.col_defs) == 14
            assert record.col_defs["PAIRS"] == "number of probe set specific probe pairs on the array"
            assert record.col_defs["ABS_CALL"] == "Whether a probe set is present, marginal, or absent; see Affymetrix Literature"
            assert record.col_defs["PM Excess"] == "number of probe pairs  where PM/MM exceeds the ratio limit (10 by default)"
            assert record.col_defs["POSITIVE"] == "number of poisitive probe pairs"
            assert record.col_defs["MM Excess"] == "Number of probe peirs where MM/PM exceeds 1/ratio limit (10 by default)"
            assert record.col_defs["ID_REF"] == "Affymetrix Probe Set Identifier"
            assert record.col_defs["NEGATIVE"] == "number of negative probe pairs"
            assert record.col_defs["VALUE"] == "Average Difference Intensity"
            assert record.col_defs["POS_FRACTION"] == "Positive/Pairs Used"
            assert record.col_defs["Experiment Name"] == "Experiment Name"
            assert record.col_defs["POS/NEG"] == "Positive/Negative"
            assert record.col_defs["PAIRS_USED"] == ""
            assert record.col_defs["Log Avg"] == ""
            assert record.col_defs["PAIRS_IN_AVG"] == "Trimmed probe pair set"
            assert len(record.table_rows) == 21
            assert len(record.table_rows[0]) == 14
            assert record.table_rows[0][0] == "ID_REF"
            assert record.table_rows[0][1] == "Experiment Name"
            assert record.table_rows[0][2] == "POSITIVE"
            assert record.table_rows[0][3] == "NEGATIVE"
            assert record.table_rows[0][4] == "PAIRS"
            assert record.table_rows[0][5] == "PAIRS_USED"
            assert record.table_rows[0][6] == "PAIRS_IN_AVG"
            assert record.table_rows[0][7] == "POS_FRACTION"
            assert record.table_rows[0][8] == "Log Avg"
            assert record.table_rows[0][9] == "PM Excess"
            assert record.table_rows[0][10] == "MM Excess"
            assert record.table_rows[0][11] == "POS/NEG"
            assert record.table_rows[0][12] == "VALUE"
            assert record.table_rows[0][13] == "ABS_CALL"
            assert len(record.table_rows[1]) == 14
            assert record.table_rows[1][0] == "IL2_at"
            assert record.table_rows[1][1] == "RHMu8LarB"
            assert record.table_rows[1][2] == "4"
            assert record.table_rows[1][3] == "4"
            assert record.table_rows[1][4] == "19"
            assert record.table_rows[1][5] == "19"
            assert record.table_rows[1][6] == "19"
            assert record.table_rows[1][7] == "0.21"
            assert record.table_rows[1][8] == "-0.58"
            assert record.table_rows[1][9] == "0"
            assert record.table_rows[1][10] == "0"
            assert record.table_rows[1][11] == "1.0"
            assert record.table_rows[1][12] == "-78"
            assert record.table_rows[1][13] == "A"
            assert len(record.table_rows[2]) == 14
            assert record.table_rows[2][0] == "IL10_at"
            assert record.table_rows[2][1] == "RHMu8LarB"
            assert record.table_rows[2][2] == "7"
            assert record.table_rows[2][3] == "4"
            assert record.table_rows[2][4] == "20"
            assert record.table_rows[2][5] == "20"
            assert record.table_rows[2][6] == "18"
            assert record.table_rows[2][7] == "0.35"
            assert record.table_rows[2][8] == "1.87"
            assert record.table_rows[2][9] == "1"
            assert record.table_rows[2][10] == "0"
            assert record.table_rows[2][11] == "1.8"
            assert record.table_rows[2][12] == "161"
            assert record.table_rows[2][13] == "A"
            assert len(record.table_rows[3]) == 14
            assert record.table_rows[3][0] == "GMCSF_at"
            assert record.table_rows[3][1] == "RHMu8LarB"
            assert record.table_rows[3][2] == "4"
            assert record.table_rows[3][3] == "4"
            assert record.table_rows[3][4] == "20"
            assert record.table_rows[3][5] == "20"
            assert record.table_rows[3][6] == "19"
            assert record.table_rows[3][7] == "0.20"
            assert record.table_rows[3][8] == "0.39"
            assert record.table_rows[3][9] == "0"
            assert record.table_rows[3][10] == "0"
            assert record.table_rows[3][11] == "1.0"
            assert record.table_rows[3][12] == "-11"
            assert record.table_rows[3][13] == "A"
            assert len(record.table_rows[4]) == 14
            assert record.table_rows[4][0] == "TNFRII_at"
            assert record.table_rows[4][1] == "RHMu8LarB"
            assert record.table_rows[4][2] == "2"
            assert record.table_rows[4][3] == "2"
            assert record.table_rows[4][4] == "20"
            assert record.table_rows[4][5] == "20"
            assert record.table_rows[4][6] == "18"
            assert record.table_rows[4][7] == "0.10"
            assert record.table_rows[4][8] == "0.48"
            assert record.table_rows[4][9] == "0"
            assert record.table_rows[4][10] == "0"
            assert record.table_rows[4][11] == "1.0"
            assert record.table_rows[4][12] == "52"
            assert record.table_rows[4][13] == "A"
            assert len(record.table_rows[5]) == 14
            assert record.table_rows[5][0] == "MIP1-B_at"
            assert record.table_rows[5][1] == "RHMu8LarB"
            assert record.table_rows[5][2] == "6"
            assert record.table_rows[5][3] == "4"
            assert record.table_rows[5][4] == "20"
            assert record.table_rows[5][5] == "20"
            assert record.table_rows[5][6] == "19"
            assert record.table_rows[5][7] == "0.30"
            assert record.table_rows[5][8] == "0.43"
            assert record.table_rows[5][9] == "0"
            assert record.table_rows[5][10] == "0"
            assert record.table_rows[5][11] == "1.5"
            assert record.table_rows[5][12] == "373"
            assert record.table_rows[5][13] == "A"
            assert len(record.table_rows[6]) == 14
            assert record.table_rows[6][0] == "IL4_at"
            assert record.table_rows[6][1] == "RHMu8LarB"
            assert record.table_rows[6][2] == "3"
            assert record.table_rows[6][3] == "3"
            assert record.table_rows[6][4] == "20"
            assert record.table_rows[6][5] == "20"
            assert record.table_rows[6][6] == "19"
            assert record.table_rows[6][7] == "0.15"
            assert record.table_rows[6][8] == "0.29"
            assert record.table_rows[6][9] == "0"
            assert record.table_rows[6][10] == "0"
            assert record.table_rows[6][11] == "1.0"
            assert record.table_rows[6][12] == "27"
            assert record.table_rows[6][13] == "A"
            assert len(record.table_rows[7]) == 14
            assert record.table_rows[7][0] == "IL12_P40_at"
            assert record.table_rows[7][1] == "RHMu8LarB"
            assert record.table_rows[7][2] == "3"
            assert record.table_rows[7][3] == "5"
            assert record.table_rows[7][4] == "20"
            assert record.table_rows[7][5] == "20"
            assert record.table_rows[7][6] == "19"
            assert record.table_rows[7][7] == "0.15"
            assert record.table_rows[7][8] == "-0.22"
            assert record.table_rows[7][9] == "0"
            assert record.table_rows[7][10] == "0"
            assert record.table_rows[7][11] == "0.6"
            assert record.table_rows[7][12] == "-163"
            assert record.table_rows[7][13] == "A"
            assert len(record.table_rows[8]) == 14
            assert record.table_rows[8][0] == "TNFa_at"
            assert record.table_rows[8][1] == "RHMu8LarB"
            assert record.table_rows[8][2] == "3"
            assert record.table_rows[8][3] == "4"
            assert record.table_rows[8][4] == "20"
            assert record.table_rows[8][5] == "20"
            assert record.table_rows[8][6] == "20"
            assert record.table_rows[8][7] == "0.15"
            assert record.table_rows[8][8] == "-0.57"
            assert record.table_rows[8][9] == "1"
            assert record.table_rows[8][10] == "0"
            assert record.table_rows[8][11] == "0.8"
            assert record.table_rows[8][12] == "-95"
            assert record.table_rows[8][13] == "A"
            assert len(record.table_rows[9]) == 14
            assert record.table_rows[9][0] == "TCRa_at"
            assert record.table_rows[9][1] == "RHMu8LarB"
            assert record.table_rows[9][2] == "1"
            assert record.table_rows[9][3] == "4"
            assert record.table_rows[9][4] == "20"
            assert record.table_rows[9][5] == "20"
            assert record.table_rows[9][6] == "19"
            assert record.table_rows[9][7] == "0.05"
            assert record.table_rows[9][8] == "-0.50"
            assert record.table_rows[9][9] == "0"
            assert record.table_rows[9][10] == "0"
            assert record.table_rows[9][11] == "0.3"
            assert record.table_rows[9][12] == "-186"
            assert record.table_rows[9][13] == "A"
            assert len(record.table_rows[10]) == 14
            assert record.table_rows[10][0] == "AFFX-BioB-5_at"
            assert record.table_rows[10][1] == "RHMu8LarB"
            assert record.table_rows[10][2] == "0"
            assert record.table_rows[10][3] == "1"
            assert record.table_rows[10][4] == "20"
            assert record.table_rows[10][5] == "20"
            assert record.table_rows[10][6] == "19"
            assert record.table_rows[10][7] == "0.00"
            assert record.table_rows[10][8] == "0.35"
            assert record.table_rows[10][9] == "0"
            assert record.table_rows[10][10] == "0"
            assert record.table_rows[10][11] == "0.0"
            assert record.table_rows[10][12] == "120"
            assert record.table_rows[10][13] == "A"
            assert len(record.table_rows[11]) == 14
            assert record.table_rows[11][0] == "AFFX-BioB-M_at"
            assert record.table_rows[11][1] == "RHMu8LarB"
            assert record.table_rows[11][2] == "0"
            assert record.table_rows[11][3] == "1"
            assert record.table_rows[11][4] == "20"
            assert record.table_rows[11][5] == "20"
            assert record.table_rows[11][6] == "19"
            assert record.table_rows[11][7] == "0.00"
            assert record.table_rows[11][8] == "0.02"
            assert record.table_rows[11][9] == "0"
            assert record.table_rows[11][10] == "0"
            assert record.table_rows[11][11] == "0.0"
            assert record.table_rows[11][12] == "-13"
            assert record.table_rows[11][13] == "A"
            assert len(record.table_rows[12]) == 14
            assert record.table_rows[12][0] == "AFFX-BioB-3_at"
            assert record.table_rows[12][1] == "RHMu8LarB"
            assert record.table_rows[12][2] == "2"
            assert record.table_rows[12][3] == "0"
            assert record.table_rows[12][4] == "20"
            assert record.table_rows[12][5] == "20"
            assert record.table_rows[12][6] == "19"
            assert record.table_rows[12][7] == "0.10"
            assert record.table_rows[12][8] == "0.38"
            assert record.table_rows[12][9] == "0"
            assert record.table_rows[12][10] == "0"
            assert record.table_rows[12][11] == "Undef"
            assert record.table_rows[12][12] == "136"
            assert record.table_rows[12][13] == "A"
            assert len(record.table_rows[13]) == 14
            assert record.table_rows[13][0] == "AFFX-BioC-5_at"
            assert record.table_rows[13][1] == "RHMu8LarB"
            assert record.table_rows[13][2] == "9"
            assert record.table_rows[13][3] == "0"
            assert record.table_rows[13][4] == "20"
            assert record.table_rows[13][5] == "20"
            assert record.table_rows[13][6] == "20"
            assert record.table_rows[13][7] == "0.45"
            assert record.table_rows[13][8] == "1.33"
            assert record.table_rows[13][9] == "0"
            assert record.table_rows[13][10] == "0"
            assert record.table_rows[13][11] == "Undef"
            assert record.table_rows[13][12] == "606"
            assert record.table_rows[13][13] == "P"
            assert len(record.table_rows[14]) == 14
            assert record.table_rows[14][0] == "AFFX-BioC-3_at"
            assert record.table_rows[14][1] == "RHMu8LarB"
            assert record.table_rows[14][2] == "2"
            assert record.table_rows[14][3] == "0"
            assert record.table_rows[14][4] == "20"
            assert record.table_rows[14][5] == "20"
            assert record.table_rows[14][6] == "19"
            assert record.table_rows[14][7] == "0.10"
            assert record.table_rows[14][8] == "0.64"
            assert record.table_rows[14][9] == "0"
            assert record.table_rows[14][10] == "0"
            assert record.table_rows[14][11] == "Undef"
            assert record.table_rows[14][12] == "257"
            assert record.table_rows[14][13] == "A"
            assert len(record.table_rows[15]) == 14
            assert record.table_rows[15][0] == "AFFX-BioDn-5_at"
            assert record.table_rows[15][1] == "RHMu8LarB"
            assert record.table_rows[15][2] == "8"
            assert record.table_rows[15][3] == "0"
            assert record.table_rows[15][4] == "20"
            assert record.table_rows[15][5] == "20"
            assert record.table_rows[15][6] == "20"
            assert record.table_rows[15][7] == "0.40"
            assert record.table_rows[15][8] == "1.23"
            assert record.table_rows[15][9] == "0"
            assert record.table_rows[15][10] == "0"
            assert record.table_rows[15][11] == "Undef"
            assert record.table_rows[15][12] == "380"
            assert record.table_rows[15][13] == "P"
            assert len(record.table_rows[16]) == 14
            assert record.table_rows[16][0] == "AFFX-BioDn-3_at"
            assert record.table_rows[16][1] == "RHMu8LarB"
            assert record.table_rows[16][2] == "16"
            assert record.table_rows[16][3] == "0"
            assert record.table_rows[16][4] == "20"
            assert record.table_rows[16][5] == "20"
            assert record.table_rows[16][6] == "19"
            assert record.table_rows[16][7] == "0.80"
            assert record.table_rows[16][8] == "2.79"
            assert record.table_rows[16][9] == "0"
            assert record.table_rows[16][10] == "0"
            assert record.table_rows[16][11] == "Undef"
            assert record.table_rows[16][12] == "2764"
            assert record.table_rows[16][13] == "P"
            assert len(record.table_rows[17]) == 14
            assert record.table_rows[17][0] == "AFFX-CreX-5_at"
            assert record.table_rows[17][1] == "RHMu8LarB"
            assert record.table_rows[17][2] == "19"
            assert record.table_rows[17][3] == "0"
            assert record.table_rows[17][4] == "20"
            assert record.table_rows[17][5] == "20"
            assert record.table_rows[17][6] == "19"
            assert record.table_rows[17][7] == "0.95"
            assert record.table_rows[17][8] == "5.65"
            assert record.table_rows[17][9] == "0"
            assert record.table_rows[17][10] == "0"
            assert record.table_rows[17][11] == "Undef"
            assert record.table_rows[17][12] == "4391"
            assert record.table_rows[17][13] == "P"
            assert len(record.table_rows[18]) == 14
            assert record.table_rows[18][0] == "AFFX-CreX-3_at"
            assert record.table_rows[18][1] == "RHMu8LarB"
            assert record.table_rows[18][2] == "19"
            assert record.table_rows[18][3] == "0"
            assert record.table_rows[18][4] == "20"
            assert record.table_rows[18][5] == "20"
            assert record.table_rows[18][6] == "20"
            assert record.table_rows[18][7] == "0.95"
            assert record.table_rows[18][8] == "6.42"
            assert record.table_rows[18][9] == "2"
            assert record.table_rows[18][10] == "0"
            assert record.table_rows[18][11] == "Undef"
            assert record.table_rows[18][12] == "10787"
            assert record.table_rows[18][13] == "P"
            assert len(record.table_rows[19]) == 14
            assert record.table_rows[19][0] == "AFFX-BioB-5_st"
            assert record.table_rows[19][1] == "RHMu8LarB"
            assert record.table_rows[19][2] == "5"
            assert record.table_rows[19][3] == "3"
            assert record.table_rows[19][4] == "20"
            assert record.table_rows[19][5] == "20"
            assert record.table_rows[19][6] == "19"
            assert record.table_rows[19][7] == "0.25"
            assert record.table_rows[19][8] == "0.48"
            assert record.table_rows[19][9] == "0"
            assert record.table_rows[19][10] == "0"
            assert record.table_rows[19][11] == "1.7"
            assert record.table_rows[19][12] == "80"
            assert record.table_rows[19][13] == "A"
            assert len(record.table_rows[20]) == 14
            assert record.table_rows[20][0] == "AFFX-BioB-M_st"
            assert record.table_rows[20][1] == "RHMu8LarB"
            assert record.table_rows[20][2] == "2"
            assert record.table_rows[20][3] == "3"
            assert record.table_rows[20][4] == "20"
            assert record.table_rows[20][5] == "20"
            assert record.table_rows[20][6] == "17"
            assert record.table_rows[20][7] == "0.10"
            assert record.table_rows[20][8] == "0.16"
            assert record.table_rows[20][9] == "0"
            assert record.table_rows[20][10] == "0"
            assert record.table_rows[20][11] == "0.7"
            assert record.table_rows[20][12] == "24"
            assert record.table_rows[20][13] == "A"

    def test_soft_ex_series(self):
        path = "Geo/soft_ex_series.txt"
        with open(path, encoding="latin") as handle:
            records = Geo.parse(handle)
            record = next(records)
            assert record.entity_type == "SERIES"
            assert record.entity_id == "Bone_marrow_stromal_cells"
            assert len(record.entity_attributes) == 13
            assert record.entity_attributes["Series_variable_description_1"] == "HS-5"
            assert record.entity_attributes["Series_variable_sample_list_1"] == "GSM10001, GSM10002"
            assert len(record.entity_attributes["Series_sample_id"]) == 4
            assert record.entity_attributes["Series_sample_id"][0] == "GSM10001"
            assert record.entity_attributes["Series_sample_id"][1] == "GSM10002"
            assert record.entity_attributes["Series_sample_id"][2] == "GSM10003"
            assert record.entity_attributes["Series_sample_id"][3] == "GSM10004"
            assert record.entity_attributes["Series_variable_2"] == "cell line"
            assert record.entity_attributes["Series_pubmed_id"] == "123456789"
            assert len(record.entity_attributes["Series_contributor"]) == 5
            assert record.entity_attributes["Series_contributor"][0] == "Jane,Doe"
            assert record.entity_attributes["Series_contributor"][1] == "John,A,Smith"
            assert record.entity_attributes["Series_contributor"][2] == "Hans,van Elton"
            assert record.entity_attributes["Series_contributor"][3] == "John,Smithers Jr"
            assert record.entity_attributes["Series_contributor"][4] == "Jie,D,Chen"
            assert record.entity_attributes["Series_summary"] == "Two human stromal cell lines, HS-5 and HS-27a, represent functionally distinct components of the bone marrow microenvironment.1,2 HS-27a supports cobblestone area formation by early hematopoietic progenitors, whereas HS-5 secretes multiple cytokines that support the proliferation of committed progenitors. These cell lines have been distributed to research groups worldwide for use as a tool to understand interactions between hematopoietic cells and their microenvironment. We have used DNA microarray technology to characterize and compare the expression of over 17 000 genes in these cell lines. Gene expression differences in cytokines/chemokines, G-protein signaling molecules, and multiple extracellular matrix proteins add to the known protein and functional characterization of the lines, leading to new insight into the differences in their support function for hematopoietic progenitors."
            assert record.entity_attributes["Series_type"] == "Cell Line Comparison"
            assert record.entity_attributes["Series_variable_1"] == "cell line"
            assert record.entity_attributes["Series_variable_description_2"] == "HS-27a"
            assert record.entity_attributes["Series_title"] == "Profiling of the functionally distinct human bone marrow stromal cell lines HS-5 and HS-27a."
            assert record.entity_attributes["Series_variable_sample_list_2"] == "GSM10003, GSM10004"
            assert record.entity_attributes["Series_overall_design"] == "We analyzed 2 arrays for HS-5 cell line and 2 arrays for HS-27a cell line"
            assert len(record.col_defs) == 0
            assert len(record.table_rows) == 0

    def test_GSM691(self):
        path = "Geo/GSM691.txt"
        with open(path, encoding="latin") as handle:
            records = Geo.parse(handle)
            record = next(records)
            assert record.entity_type == "SAMPLE"
            assert record.entity_id == "GSM691"
            assert len(record.entity_attributes) == 20
            assert record.entity_attributes["Sample_submitter_institute"] == "National Cancer Institute"
            assert record.entity_attributes["Sample_submitter_department"] == "Cancer Genome Anatomy Project"
            assert record.entity_attributes["Sample_submitter_web_link"] == "http://cgap.nci.nih.gov/"
            assert len(record.entity_attributes["Sample_description"]) == 12
            assert record.entity_attributes["Sample_description"][0] == "This library represents a Cancer Genome Anatomy Project library, which was either produced through CGAP funding, or donated to CGAP."
            assert record.entity_attributes["Sample_description"][1] == "The Cancer Genome Anatomy Project (CGAP: http://cgap.nci.nih.gov) is an interdisciplinary program established and administered by the National Cancer Institute (NCI: http://www.nci.nih.gov) to generate the information and technological tools needed to decipher the molecular anatomy of the cancer cell."
            assert record.entity_attributes["Sample_description"][2] == "Library constructed by Riggins laboratory Tissue supplied by Jeffrey Marks, Ph.D."
            assert record.entity_attributes["Sample_description"][3] == "Organ: Breast"
            assert record.entity_attributes["Sample_description"][4] == "Tissue_type: normal epithelial organoids"
            assert record.entity_attributes["Sample_description"][5] == "Library treatment: non-normalized"
            assert record.entity_attributes["Sample_description"][6] == "Tissue description: Breast, Isolated normal epithelial organoids. Derived from a reduction mammoplasty."
            assert record.entity_attributes["Sample_description"][7] == "Tissue supplier: Jeffrey Marks, Ph.D."
            assert record.entity_attributes["Sample_description"][8] == "Sample type: Bulk"
            assert record.entity_attributes["Sample_description"][9] == "Producer: Riggins Laboratory"
            assert record.entity_attributes["Sample_description"][10] == "Clones generated to date: 768"
            assert record.entity_attributes["Sample_description"][11] == "Sequences generated to date: 572"
            assert len(record.entity_attributes["Sample_author"]) == 3
            assert record.entity_attributes["Sample_author"][0] == "Jeffrey,,Marks"
            assert record.entity_attributes["Sample_author"][1] == "Gregory,J,Riggins"
            assert record.entity_attributes["Sample_author"][2] == "Robert,L,Strausberg"
            assert record.entity_attributes["Sample_web_link"] == "http://cgap.nci.nih.gov"
            assert record.entity_attributes["Sample_submitter_phone"] == "301-496-1550"
            assert record.entity_attributes["Sample_series_id"] == "GSE14"
            assert record.entity_attributes["Sample_tag_count"] == "7165"
            assert record.entity_attributes["Sample_type"] == "sage"
            assert record.entity_attributes["Sample_submitter_name"] == "Robert,L,Strausberg"
            assert record.entity_attributes["Sample_platform_id"] == "GPL4"
            assert record.entity_attributes["Sample_submitter_city"] == "Bethesda,MD,20892,USA"
            assert record.entity_attributes["Sample_status"] == "Public on Nov 28 2001"
            assert record.entity_attributes["Sample_anchor"] == "NlaIII"
            assert record.entity_attributes["Sample_title"] == "SAGE_Duke_40N"
            assert record.entity_attributes["Sample_organism"] == "Homo sapiens"
            assert record.entity_attributes["Sample_target_source"] == "Breast, isolated normal epithelial organoids"
            assert record.entity_attributes["Sample_submission_date"] == "Nov 28 2001"
            assert record.entity_attributes["Sample_submitter_email"] == "cgapbs-r@mail.nih.gov"
            assert len(record.col_defs) == 3
            assert record.col_defs["COUNT"] == "Absolute tag count"
            assert record.col_defs["TPM"] == "Tags per million, or (1000000*COUNT)/(Total tags)"
            assert record.col_defs["TAG"] == 'Ten base SAGE tag, LINK_PRE:"http://www.ncbi.nlm.nih.gov/SAGE/SAGEtag.cgi?tag="'
            assert len(record.table_rows) == 21
            assert len(record.table_rows[0]) == 3
            assert record.table_rows[0][0] == "TAG"
            assert record.table_rows[0][1] == "COUNT"
            assert record.table_rows[0][2] == "TPM"
            assert len(record.table_rows[1]) == 3
            assert record.table_rows[1][0] == "TTGGGGTTTC"
            assert record.table_rows[1][1] == "202"
            assert record.table_rows[1][2] == "28192.6"
            assert len(record.table_rows[2]) == 3
            assert record.table_rows[2][0] == "TAGGTTGTCT"
            assert record.table_rows[2][1] == "129"
            assert record.table_rows[2][2] == "18004.2"
            assert len(record.table_rows[3]) == 3
            assert record.table_rows[3][0] == "GAGGGAGTTT"
            assert record.table_rows[3][1] == "109"
            assert record.table_rows[3][2] == "15212.8"
            assert len(record.table_rows[4]) == 3
            assert record.table_rows[4][0] == "TGCACGTTTT"
            assert record.table_rows[4][1] == "92"
            assert record.table_rows[4][2] == "12840.2"
            assert len(record.table_rows[5]) == 3
            assert record.table_rows[5][0] == "CTGGGTTAAT"
            assert record.table_rows[5][1] == "83"
            assert record.table_rows[5][2] == "11584.1"
            assert len(record.table_rows[6]) == 3
            assert record.table_rows[6][0] == "GTTGTGGTTA"
            assert record.table_rows[6][1] == "82"
            assert record.table_rows[6][2] == "11444.5"
            assert len(record.table_rows[7]) == 3
            assert record.table_rows[7][0] == "GATCCCAACT"
            assert record.table_rows[7][1] == "63"
            assert record.table_rows[7][2] == "8792.74"
            assert len(record.table_rows[8]) == 3
            assert record.table_rows[8][0] == "TGCAGTCACT"
            assert record.table_rows[8][1] == "59"
            assert record.table_rows[8][2] == "8234.47"
            assert len(record.table_rows[9]) == 3
            assert record.table_rows[9][0] == "GGATTTGGCC"
            assert record.table_rows[9][1] == "58"
            assert record.table_rows[9][2] == "8094.91"
            assert len(record.table_rows[10]) == 3
            assert record.table_rows[10][0] == "GGGCTGGGGT"
            assert record.table_rows[10][1] == "56"
            assert record.table_rows[10][2] == "7815.77"
            assert len(record.table_rows[11]) == 3
            assert record.table_rows[11][0] == "ATAATTCTTT"
            assert record.table_rows[11][1] == "44"
            assert record.table_rows[11][2] == "6140.96"
            assert len(record.table_rows[12]) == 3
            assert record.table_rows[12][0] == "CTTCCTTGCC"
            assert record.table_rows[12][1] == "42"
            assert record.table_rows[12][2] == "5861.83"
            assert len(record.table_rows[13]) == 3
            assert record.table_rows[13][0] == "TTGGTCCTCT"
            assert record.table_rows[13][1] == "40"
            assert record.table_rows[13][2] == "5582.69"
            assert len(record.table_rows[14]) == 3
            assert record.table_rows[14][0] == "GGCAAGCCCC"
            assert record.table_rows[14][1] == "36"
            assert record.table_rows[14][2] == "5024.42"
            assert len(record.table_rows[15]) == 3
            assert record.table_rows[15][0] == "AACTAAAAAA"
            assert record.table_rows[15][1] == "34"
            assert record.table_rows[15][2] == "4745.29"
            assert len(record.table_rows[16]) == 3
            assert record.table_rows[16][0] == "AGGGCTTCCA"
            assert record.table_rows[16][1] == "34"
            assert record.table_rows[16][2] == "4745.29"
            assert len(record.table_rows[17]) == 3
            assert record.table_rows[17][0] == "AGGCTACGGA"
            assert record.table_rows[17][1] == "33"
            assert record.table_rows[17][2] == "4605.72"
            assert len(record.table_rows[18]) == 3
            assert record.table_rows[18][0] == "GTGAAACCCC"
            assert record.table_rows[18][1] == "32"
            assert record.table_rows[18][2] == "4466.15"
            assert len(record.table_rows[19]) == 3
            assert record.table_rows[19][0] == "AACTAACAAA"
            assert record.table_rows[19][1] == "31"
            assert record.table_rows[19][2] == "4326.59"
            assert len(record.table_rows[20]) == 3
            assert record.table_rows[20][0] == "GAAAAATGGT"
            assert record.table_rows[20][1] == "30"
            assert record.table_rows[20][2] == "4187.02"

    def test_soft_ex_family(self):
        path = "Geo/soft_ex_family.txt"
        with open(path, encoding="latin") as handle:
            records = Geo.parse(handle)
            record = next(records)
            assert record.entity_type == "PLATFORM"
            assert record.entity_id == "Murine 15K long oligo array version 2.0"
            assert len(record.entity_attributes) == 12
            assert record.entity_attributes["Platform_title"] == "Murine 15K long oligo array version 2.0"
            assert record.entity_attributes["Platform_web_link"] == "http://www.microarray.protocols.html"
            assert record.entity_attributes["platform_table_end"] == ""
            assert record.entity_attributes["Platform_support"] == "glass"
            assert record.entity_attributes["Platform_manufacturer"] == "Un. London microarray facility"
            assert record.entity_attributes["Platform_coating"] == "polysine"
            assert record.entity_attributes["Platform_technology"] == "spotted oligonucleotide"
            assert record.entity_attributes["platform_table_begin"] == ""
            assert len(record.entity_attributes["Platform_manufacture_protocol"]) == 12
            assert record.entity_attributes["Platform_manufacture_protocol"][0] == "1.  Oligos are arrayed in Greiner 384-well flat-bottom plates. Each well contains 600 pmol of 70-mer oligo."
            assert record.entity_attributes["Platform_manufacture_protocol"][1] == "2. Resuspend oligos in water to 20 uM and rearray 5 \xb5L into 384-well, Genetix polystyrene V-bottom plates (cat# X6004)."
            assert record.entity_attributes["Platform_manufacture_protocol"][2] == "3. Allow Genetix plates to dry through passive water evaporation in a protected environment (e.g., chemical hood)."
            assert record.entity_attributes["Platform_manufacture_protocol"][3] == "4. Before printing, add 5 \xb5L of 1X Printing Buffer to each well. This can be done the night before a print run is started."
            assert record.entity_attributes["Platform_manufacture_protocol"][4] == "5. Seal plates with Corning seals."
            assert record.entity_attributes["Platform_manufacture_protocol"][5] == "6. Incubate at 37\xb0C for 30 minutes to aid resuspension of DNA."
            assert record.entity_attributes["Platform_manufacture_protocol"][6] == "7. Shake plates near maximum rotational speed on flat-bed shaker for 1 minute."
            assert record.entity_attributes["Platform_manufacture_protocol"][7] == "8. Centrifuge plates at 2000 rpm for 3 minutes."
            assert record.entity_attributes["Platform_manufacture_protocol"][8] == "9. Remove seals and cover with plate lids. Place in appropriate location of plate cassette. This should be done with first plates just before print run is started to minimize evaporation time before printing. For second and third cassettes, wait until 30 minutes before next cassette is needed to begin centrifugation."
            assert record.entity_attributes["Platform_manufacture_protocol"][9] == "10. Make sure plates rest behind both holding clips in the cassettes. Push plates back into the cassettes as far as they will go, putting them in the proper position for the server arm."
            assert record.entity_attributes["Platform_manufacture_protocol"][10] == "11. After the print run is completed, allow plates to dry through passive evaporation in a protected environment."
            assert record.entity_attributes["Platform_manufacture_protocol"][11] == "12. For each subsequent preparation of these plates for a print run, add water to the wells instead of sodium phosphate buffer. The amount of water should be decreased by 0.25 \xb5L per print run, as this is the amount drawn up by the pin capillary during each dip."
            assert record.entity_attributes["Platform_organism"] == "Mus musculus"
            assert len(record.entity_attributes["Platform_contributor"]) == 5
            assert record.entity_attributes["Platform_contributor"][0] == "Jane,Doe"
            assert record.entity_attributes["Platform_contributor"][1] == "John,A,Smith"
            assert record.entity_attributes["Platform_contributor"][2] == "Hans,van Elton"
            assert record.entity_attributes["Platform_contributor"][3] == "John,Smithers Jr"
            assert record.entity_attributes["Platform_contributor"][4] == "Jie,D,Chen"
            assert record.entity_attributes["Platform_distribution"] == "non-commercial"
            assert len(record.col_defs) == 6
            assert record.col_defs["Gene_Desc"] == "Gene description"
            assert record.col_defs["SEQUENCE"] == "Probe sequence information"
            assert record.col_defs["Gene_Sym"] == "Gene symbols"
            assert record.col_defs["GB_ACC"] == 'GenBank accession number of sequence used to design oligonucleotide probe   LINK_PRE:"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Search&db=Nucleotide&term="'
            assert record.col_defs["SPOT_ID"] == "alternative identifier"
            assert record.col_defs["ID"] == ""
            assert len(record.table_rows) == 21
            assert len(record.table_rows[0]) == 6
            assert record.table_rows[0][0] == "ID"
            assert record.table_rows[0][1] == "GB_ACC"
            assert record.table_rows[0][2] == "Gene_Desc"
            assert record.table_rows[0][3] == "Gene_Sym"
            assert record.table_rows[0][4] == "SPOT_ID"
            assert record.table_rows[0][5] == "SEQUENCE"
            assert len(record.table_rows[1]) == 6
            assert record.table_rows[1][0] == "1"
            assert record.table_rows[1][1] == "U02079"
            assert record.table_rows[1][2] == "nuclear factor of activated T-cells, cytoplasmic 2"
            assert record.table_rows[1][3] == "Nfatc2"
            assert record.table_rows[1][4] == ""
            assert record.table_rows[1][5] == "ACCTGGATGACGCAGCCACTTCAGAAAGCTGGGTTGGGACAGAAAGGTATATAGAGAGAAAATTTTGGAA"
            assert len(record.table_rows[2]) == 6
            assert record.table_rows[2][0] == "2"
            assert record.table_rows[2][1] == "NM_008154"
            assert record.table_rows[2][2] == "G-protein coupled receptor 3"
            assert record.table_rows[2][3] == "Gpr3"
            assert record.table_rows[2][4] == ""
            assert record.table_rows[2][5] == "CTGTACAATGCTCTCACTTACTACTCAGAGACAACGGTAACTCGGACTTATGTGATGCTGGCCTTGGTGT"
            assert len(record.table_rows[3]) == 6
            assert record.table_rows[3][0] == "3"
            assert record.table_rows[3][1] == "AK015719"
            assert record.table_rows[3][2] == "tropomodulin 2"
            assert record.table_rows[3][3] == "Tmod2"
            assert record.table_rows[3][4] == ""
            assert record.table_rows[3][5] == "CACCAGGCTCAGTGCCTAGTATCGGCTTCACCTAGTGTGGTTACTCAGGGCACGCAGAGCTACAGAACAC"
            assert len(record.table_rows[4]) == 6
            assert record.table_rows[4][0] == "4"
            assert record.table_rows[4][1] == "AK003367"
            assert record.table_rows[4][2] == "mitochondrial ribosomal protein L15"
            assert record.table_rows[4][3] == "Mrpl15"
            assert record.table_rows[4][4] == ""
            assert record.table_rows[4][5] == "CAAGAAGTCTAGAAATTCTGTGCAAGCCTATTCCATTCTTTCTGCGGGGACAACCAATTCCGAAAAGAAT"
            assert len(record.table_rows[5]) == 6
            assert record.table_rows[5][0] == "5"
            assert record.table_rows[5][1] == "BC003333"
            assert record.table_rows[5][2] == "RIKEN cDNA 0610033I05 gene"
            assert record.table_rows[5][3] == "0610033I05Rik"
            assert record.table_rows[5][4] == ""
            assert record.table_rows[5][5] == "AGAACTGGGTGGCAGATATCCTAGAGTTTTGACCAACGTTCACAGCACACATATTGATCTTATAGGACCT"
            assert len(record.table_rows[6]) == 6
            assert record.table_rows[6][0] == "6"
            assert record.table_rows[6][1] == "NM_008462"
            assert record.table_rows[6][2] == "killer cell lectin-like receptor, subfamily A, member 2"
            assert record.table_rows[6][3] == "Klra2"
            assert record.table_rows[6][4] == ""
            assert record.table_rows[6][5] == "TGAATTGAAGTTCCTTAAATCCCAACTTCAAAGAAACACATACTGGATTTCACTGACACATCATAAAAGC"
            assert len(record.table_rows[7]) == 6
            assert record.table_rows[7][0] == "7"
            assert record.table_rows[7][1] == "NM_008029"
            assert record.table_rows[7][2] == "FMS-like tyrosine kinase 4"
            assert record.table_rows[7][3] == "Flt4"
            assert record.table_rows[7][4] == ""
            assert record.table_rows[7][5] == "GAGGTGCTGTGGGATGACCGCCGGGGCATGCGGGTGCCCACTCAACTGTTGCGCGATGCCCTGTACCTGC"
            assert len(record.table_rows[8]) == 6
            assert record.table_rows[8][0] == "8"
            assert record.table_rows[8][1] == "NM_054088"
            assert record.table_rows[8][2] == "adiponutrin"
            assert record.table_rows[8][3] == "Adpn"
            assert record.table_rows[8][4] == ""
            assert record.table_rows[8][5] == "GTCTGAGTTCCATTCCAAAGACGAAGTCGTGGATGCCCTGGTGTGTTCCTGCTTCATTCCCCTCTTCTCT"
            assert len(record.table_rows[9]) == 6
            assert record.table_rows[9][0] == "9"
            assert record.table_rows[9][1] == "NM_009750"
            assert record.table_rows[9][2] == "nerve growth factor receptor (TNFRSF16) associated protein 1"
            assert record.table_rows[9][3] == "Ngfrap1"
            assert record.table_rows[9][4] == ""
            assert record.table_rows[9][5] == "TACAGCTGAGAAATTGTCTACGCATCCTTATGGGGGAGCTGTCTAACCACCACGATCACCATGATGAATT"
            assert len(record.table_rows[10]) == 6
            assert record.table_rows[10][0] == "10"
            assert record.table_rows[10][1] == "AB045323"
            assert record.table_rows[10][2] == "DNA segment, Chr 8, ERATO Doi 594, expressed"
            assert record.table_rows[10][3] == "D8Ertd594e"
            assert record.table_rows[10][4] == ""
            assert record.table_rows[10][5] == "GATTCAGACTCGGGAGGAGCATCCCAACCTCTCCTTGAGGATAAAGGCCTGAGCGATTGCCCTGGGGAGC"
            assert len(record.table_rows[11]) == 6
            assert record.table_rows[11][0] == "11"
            assert record.table_rows[11][1] == "AK005789"
            assert record.table_rows[11][2] == "dynein, cytoplasmic, light chain 2B"
            assert record.table_rows[11][3] == "Dncl2b"
            assert record.table_rows[11][4] == ""
            assert record.table_rows[11][5] == "TGCAGAAGGCATTCCAATCCGAACAACCCTGGACAACTCCACAACGGTTCAGTATGCGGGTCTTCTCCAC"
            assert len(record.table_rows[12]) == 6
            assert record.table_rows[12][0] == "12"
            assert record.table_rows[12][1] == "NM_010517"
            assert record.table_rows[12][2] == "insulin-like growth factor binding protein 4"
            assert record.table_rows[12][3] == "Igfbp4"
            assert record.table_rows[12][4] == ""
            assert record.table_rows[12][5] == "GGAGAAGCTGGCGCGCTGCCGCCCCCCCGTGGGTTGCGAGGAGTTGGTGCGGGAGCCAGGCTGCGGTTGT"
            assert len(record.table_rows[13]) == 6
            assert record.table_rows[13][0] == "13"
            assert record.table_rows[13][1] == "AK010722"
            assert record.table_rows[13][2] == "RIKEN cDNA 2410075D05 gene"
            assert record.table_rows[13][3] == "2410075D05Rik"
            assert record.table_rows[13][4] == ""
            assert record.table_rows[13][5] == "GGAGCATCTGGAGTTCCGCTTACCGGAAATAAAGTCTTTACTATCGGTGATTGGAGGGCAGTTCACTAAC"
            assert len(record.table_rows[14]) == 6
            assert record.table_rows[14][0] == "14"
            assert record.table_rows[14][1] == "AK003755"
            assert record.table_rows[14][2] == "DNA segment, Chr 4, ERATO Doi 421, expressed"
            assert record.table_rows[14][3] == "D4Ertd421e"
            assert record.table_rows[14][4] == ""
            assert record.table_rows[14][5] == "AGCAAAGAGATCTCCCTCAGTGTGCCCATAGGTGGCGGTGCGAGCTTGCGGTTATTGGCCAGTGACTTGC"
            assert len(record.table_rows[15]) == 6
            assert record.table_rows[15][0] == "15"
            assert record.table_rows[15][1] == "BC003241"
            assert record.table_rows[15][2] == "cleavage stimulation factor, 3' pre-RNA, subunit 3"
            assert record.table_rows[15][3] == "Cstf3"
            assert record.table_rows[15][4] == ""
            assert record.table_rows[15][5] == "AAATTAGAAGAAAATCCATATGACCTTGATGCTTGGAGCATTCTCATTCGAGAGGCACAGAATCAACCTA"
            assert len(record.table_rows[16]) == 6
            assert record.table_rows[16][0] == "16"
            assert record.table_rows[16][1] == "AK004937"
            assert record.table_rows[16][2] == "RIKEN cDNA 1300007O09 gene"
            assert record.table_rows[16][3] == "1300007O09Rik"
            assert record.table_rows[16][4] == ""
            assert record.table_rows[16][5] == "CAGACACAAACCCTAGGTTGTATTGTAGACCGGAGTTTAAGCAGGCACTACCTGTCTGTCTTTTCTTCAT"
            assert len(record.table_rows[17]) == 6
            assert record.table_rows[17][0] == "17"
            assert record.table_rows[17][1] == "AK004524"
            assert record.table_rows[17][2] == "unnamed protein product; hypothetical SOCS domain"
            assert record.table_rows[17][3] == ""
            assert record.table_rows[17][4] == ""
            assert record.table_rows[17][5] == "CGGAGCCCTGCGCGCCCAGAGCCCCCTCCCACCCGCTTCCACCAAGTGCATGGAGCCAACATCCGCATGG"
            assert len(record.table_rows[18]) == 6
            assert record.table_rows[18][0] == "18"
            assert record.table_rows[18][1] == "NM_025999"
            assert record.table_rows[18][2] == "RIKEN cDNA 2610110L04 gene"
            assert record.table_rows[18][3] == "2610110L04Rik"
            assert record.table_rows[18][4] == ""
            assert record.table_rows[18][5] == "TGCATTGATAAATGGAGTGATCGACACAGGAACTGCCCCATTTGTCGCCTACAGATGACTGGAGCAAATG"
            assert len(record.table_rows[19]) == 6
            assert record.table_rows[19][0] == "19"
            assert record.table_rows[19][1] == ""
            assert record.table_rows[19][2] == ""
            assert record.table_rows[19][3] == ""
            assert record.table_rows[19][4] == "-- CONTROL"
            assert record.table_rows[19][5] == ""
            assert len(record.table_rows[20]) == 6
            assert record.table_rows[20][0] == "20"
            assert record.table_rows[20][1] == "NM_023120"
            assert record.table_rows[20][2] == "guanine nucleotide binding protein (G protein), beta polypeptide 1-like"
            assert record.table_rows[20][3] == "Gnb1l"
            assert record.table_rows[20][4] == ""
            assert record.table_rows[20][5] == "ACCGCCTGGTCCCAGATTTGTCCTCCGAGGCACACAGTCGGCTGTGAACACGCTCCATTTCTGCCCACCA"
            record = next(records)
            assert record.entity_type == "SAMPLE"
            assert record.entity_id == "Control Embyronic Stem Cell Replicate 1"
            assert len(record.entity_attributes) == 24
            assert record.entity_attributes["Sample_extract_protocol_ch2"] == "TriZol procedure"
            assert record.entity_attributes["Sample_hyb_protocol"] == "Oligoarray control targets and hybridization buffer (Agilent In Situ Hybridization Kit Plus) were added, and samples were applied to microarrays enclosed in Agilent SureHyb-enabled hybridization chambers. After hybridization, slides were washed sequentially with 6x SSC/0.005% Triton X-102 and 0.1x SSC/0.005% Triton X-102 before scanning. Slides were hybridized for 17 h at 60\xb0C in a rotating oven, and washed."
            assert record.entity_attributes["Sample_extract_protocol_ch1"] == "TriZol procedure"
            assert record.entity_attributes["Sample_platform_id"] == "Murine 15K long oligo array version 2.0"
            assert record.entity_attributes["Sample_title"] == "Control Embyronic Stem Cell Replicate 1"
            assert record.entity_attributes["Sample_supplementary_file"] == "file1.gpr"
            assert record.entity_attributes["Sample_organism_ch2"] == "Mus musculus"
            assert record.entity_attributes["Sample_organism_ch1"] == "Mus musculus"
            assert record.entity_attributes["Sample_label_ch1"] == "Cy5"
            assert len(record.entity_attributes["Sample_scan_protocol"]) == 2
            assert record.entity_attributes["Sample_scan_protocol"][0] == "Scanned on an Agilent G2565AA scanner."
            assert record.entity_attributes["Sample_scan_protocol"][1] == "Images were quantified using Agilent Feature Extraction Software (version A.7.5)."
            assert record.entity_attributes["sample_table_begin"] == ""
            assert record.entity_attributes["Sample_label_protocol_ch2"] == "10 \xb5g of total RNA were primed with 2 \xb5l of 100 \xb5M T16N2 DNA primer at 70\xb0C for 10 min, then reversed transcribed at 42\xb0C for 1 h in the presence of 400 U SuperScript II RTase (Invitrogen), and 100 \xb5M each dATP, dTTP, dGTP, with 25 \xb5M dCTP, 25 \xb5M Cy5-labeled dCTP (NEN Life Science, Boston, MA), and RNase inhibitor (Invitrogen). RNA was then degraded with RNase A, and labeled cDNAs were purified using QIAquick PCR columns (Qiagen)."
            assert record.entity_attributes["Sample_label_protocol_ch1"] == "10 \xb5g of total RNA were primed with 2 \xb5l of 100 \xb5M T16N2 DNA primer at 70\xb0C for 10 min, then reversed transcribed at 42\xb0C for 1 h in the presence of 400 U SuperScript II RTase (Invitrogen), and 100 \xb5M each dATP, dTTP, dGTP, with 25 \xb5M dCTP, 25 \xb5M Cy5-labeled dCTP (NEN Life Science, Boston, MA), and RNase inhibitor (Invitrogen). RNA was then degraded with RNase A, and labeled cDNAs were purified using QIAquick PCR columns (Qiagen)."
            assert record.entity_attributes["Sample_data_processing"] == "LOWESS normalized, background subtracted VALUE data obtained from log of processed Red signal/processed Green signal."
            assert record.entity_attributes["sample_table_end"] == ""
            assert record.entity_attributes["Sample_label_ch2"] == "Cy3"
            assert record.entity_attributes["Sample_description"] == "Biological replicate 1 of 4. Control embryonic stem cells, untreated, harvested after several passages."
            assert record.entity_attributes["Sample_source_name_ch1"] == "Total RNA from murine ES-D3 embryonic stem cells labeled with Cyanine-5 (red)."
            assert record.entity_attributes["Sample_source_name_ch2"] == "Total RNA from pooled whole mouse embryos e17.5, labeled with Cyanine-3 (green)."
            assert record.entity_attributes["Sample_molecule_ch1"] == "total RNA"
            assert record.entity_attributes["Sample_molecule_ch2"] == "total RNA"
            assert record.entity_attributes["Sample_growth_protocol_ch1"] == "ES cells were kept in an undifferentiated, pluripotent state by using 1000 IU/ml leukemia inhibitory factor (LIF; Chemicon, ESGRO, ESG1107), and grown on top of murine embryonic fibroblasts feeder layer inactivated by 10 ug/ml of mitomycin C (Sigma, St. Louis). ES cells were cultured on 0.1% gelatin-coated plastic dishes in ES medium containing Dulbecco modified Eagle medium supplemented with 15% fetal calf serum, 0.1 mM beta-mercaptoethanol, 2 mM glutamine, and 0.1 mN non-essential amino acids."
            assert len(record.entity_attributes["Sample_characteristics_ch1"]) == 4
            assert record.entity_attributes["Sample_characteristics_ch1"][0] == "ES-D3 cell line (CRL-1934)"
            assert record.entity_attributes["Sample_characteristics_ch1"][1] == "Age: day 4"
            assert record.entity_attributes["Sample_characteristics_ch1"][2] == "Tissue: blastocytes"
            assert record.entity_attributes["Sample_characteristics_ch1"][3] == "Strain: 129/Sv mice"
            assert len(record.entity_attributes["Sample_characteristics_ch2"]) == 3
            assert record.entity_attributes["Sample_characteristics_ch2"][0] == "Strain: C57BL/6"
            assert record.entity_attributes["Sample_characteristics_ch2"][1] == "Age: e17.5 d"
            assert record.entity_attributes["Sample_characteristics_ch2"][2] == "Tissue: whole embryo"
            assert len(record.col_defs) == 6
            assert record.col_defs["ID_REF"] == ""
            assert record.col_defs["VALUE"] == "log(REDsignal/GREENsignal) per feature (processed signals used)."
            assert record.col_defs["gProcessedSignal"] == 'Dye-normalized signal after surrogate "algorithm," green "channel," used for computation of log ratio.'
            assert record.col_defs["LogRatioError"] == "error of the log ratio calculated according to the error model chosen."
            assert record.col_defs["PValueLogRatio"] == "Significance level of the Log Ratio computed for a feature."
            assert record.col_defs["rProcessedSignal"] == 'Dye-normalized signal after surrogate "algorithm," red "channel," used for computation of log ratio.'
            assert len(record.table_rows) == 21
            assert len(record.table_rows[0]) == 6
            assert record.table_rows[0][0] == "ID_REF"
            assert record.table_rows[0][1] == "VALUE"
            assert record.table_rows[0][2] == "LogRatioError"
            assert record.table_rows[0][3] == "PValueLogRatio"
            assert record.table_rows[0][4] == "gProcessedSignal"
            assert record.table_rows[0][5] == "rProcessedSignal"
            assert len(record.table_rows[1]) == 6
            assert record.table_rows[1][0] == "1"
            assert record.table_rows[1][1] == "-1.6274758"
            assert record.table_rows[1][2] == "1.36E-01"
            assert record.table_rows[1][3] == "6.41E-33"
            assert record.table_rows[1][4] == "9.13E+03"
            assert record.table_rows[1][5] == "2.15E+02"
            assert len(record.table_rows[2]) == 6
            assert record.table_rows[2][0] == "2"
            assert record.table_rows[2][1] == "0.1412248"
            assert record.table_rows[2][2] == "1.34E+00"
            assert record.table_rows[2][3] == "1.00E+00"
            assert record.table_rows[2][4] == "4.14E+01"
            assert record.table_rows[2][5] == "5.72E+01"
            assert len(record.table_rows[3]) == 6
            assert record.table_rows[3][0] == "3"
            assert record.table_rows[3][1] == "0.1827684"
            assert record.table_rows[3][2] == "5.19E-02"
            assert record.table_rows[3][3] == "4.33E-04"
            assert record.table_rows[3][4] == "5.13E+03"
            assert record.table_rows[3][5] == "7.81E+03"
            assert len(record.table_rows[4]) == 6
            assert record.table_rows[4][0] == "4"
            assert record.table_rows[4][1] == "-0.3932267"
            assert record.table_rows[4][2] == "6.08E-02"
            assert record.table_rows[4][3] == "1.02E-10"
            assert record.table_rows[4][4] == "4.65E+03"
            assert record.table_rows[4][5] == "1.88E+03"
            assert len(record.table_rows[5]) == 6
            assert record.table_rows[5][0] == "5"
            assert record.table_rows[5][1] == "-0.9865994"
            assert record.table_rows[5][2] == "1.05E-01"
            assert record.table_rows[5][3] == "6.32E-21"
            assert record.table_rows[5][4] == "2.91E+03"
            assert record.table_rows[5][5] == "3.01E+02"
            assert len(record.table_rows[6]) == 6
            assert record.table_rows[6][0] == "6"
            assert record.table_rows[6][1] == "0.0238812"
            assert record.table_rows[6][2] == "1.02E-01"
            assert record.table_rows[6][3] == "8.15E-01"
            assert record.table_rows[6][4] == "7.08E+02"
            assert record.table_rows[6][5] == "7.48E+02"
            assert len(record.table_rows[7]) == 6
            assert record.table_rows[7][0] == "7"
            assert record.table_rows[7][1] == "-1.4841822"
            assert record.table_rows[7][2] == "1.25E-01"
            assert record.table_rows[7][3] == "1.42E-32"
            assert record.table_rows[7][4] == "1.02E+04"
            assert record.table_rows[7][5] == "3.36E+02"
            assert len(record.table_rows[8]) == 6
            assert record.table_rows[8][0] == "8"
            assert record.table_rows[8][1] == "-1.8261356"
            assert record.table_rows[8][2] == "4.15E-01"
            assert record.table_rows[8][3] == "1.10E-05"
            assert record.table_rows[8][4] == "7.19E+02"
            assert record.table_rows[8][5] == "1.07E+01"
            assert len(record.table_rows[9]) == 6
            assert record.table_rows[9][0] == "9"
            assert record.table_rows[9][1] == "-1.0344779"
            assert record.table_rows[9][2] == "1.78E+00"
            assert record.table_rows[9][3] == "1.00E+00"
            assert record.table_rows[9][4] == "9.62E+01"
            assert record.table_rows[9][5] == "8.89E+00"
            assert len(record.table_rows[10]) == 6
            assert record.table_rows[10][0] == "10"
            assert record.table_rows[10][1] == "0.2405891"
            assert record.table_rows[10][2] == "3.09E-01"
            assert record.table_rows[10][3] == "4.36E-01"
            assert record.table_rows[10][4] == "1.61E+02"
            assert record.table_rows[10][5] == "2.80E+02"
            assert len(record.table_rows[11]) == 6
            assert record.table_rows[11][0] == "11"
            assert record.table_rows[11][1] == "0.3209366"
            assert record.table_rows[11][2] == "3.59E-01"
            assert record.table_rows[11][3] == "3.71E-01"
            assert record.table_rows[11][4] == "1.25E+02"
            assert record.table_rows[11][5] == "2.61E+02"
            assert len(record.table_rows[12]) == 6
            assert record.table_rows[12][0] == "12"
            assert record.table_rows[12][1] == "0.358304"
            assert record.table_rows[12][2] == "2.06E+00"
            assert record.table_rows[12][3] == "1.00E+00"
            assert record.table_rows[12][4] == "2.04E+01"
            assert record.table_rows[12][5] == "4.66E+01"
            assert len(record.table_rows[13]) == 6
            assert record.table_rows[13][0] == "13"
            assert record.table_rows[13][1] == "-0.0122072"
            assert record.table_rows[13][2] == "3.64E-01"
            assert record.table_rows[13][3] == "9.73E-01"
            assert record.table_rows[13][4] == "1.84E+02"
            assert record.table_rows[13][5] == "1.79E+02"
            assert len(record.table_rows[14]) == 6
            assert record.table_rows[14][0] == "14"
            assert record.table_rows[14][1] == "-1.5480396"
            assert record.table_rows[14][2] == "1.30E-01"
            assert record.table_rows[14][3] == "7.21E-33"
            assert record.table_rows[14][4] == "1.02E+04"
            assert record.table_rows[14][5] == "2.90E+02"
            assert len(record.table_rows[15]) == 6
            assert record.table_rows[15][0] == "15"
            assert record.table_rows[15][1] == "0.0073419"
            assert record.table_rows[15][2] == "2.98E-01"
            assert record.table_rows[15][3] == "9.80E-01"
            assert record.table_rows[15][4] == "2.21E+02"
            assert record.table_rows[15][5] == "2.25E+02"
            assert len(record.table_rows[16]) == 6
            assert record.table_rows[16][0] == "16"
            assert record.table_rows[16][1] == "-0.2267015"
            assert record.table_rows[16][2] == "9.44E-01"
            assert record.table_rows[16][3] == "8.10E-01"
            assert record.table_rows[16][4] == "8.90E+01"
            assert record.table_rows[16][5] == "5.28E+01"
            assert len(record.table_rows[17]) == 6
            assert record.table_rows[17][0] == "17"
            assert record.table_rows[17][1] == "-0.1484023"
            assert record.table_rows[17][2] == "8.01E-01"
            assert record.table_rows[17][3] == "8.53E-01"
            assert record.table_rows[17][4] == "9.65E+01"
            assert record.table_rows[17][5] == "6.86E+01"
            assert len(record.table_rows[18]) == 6
            assert record.table_rows[18][0] == "18"
            assert record.table_rows[18][1] == "-0.6122195"
            assert record.table_rows[18][2] == "1.28E-01"
            assert record.table_rows[18][3] == "1.69E-06"
            assert record.table_rows[18][4] == "1.12E+03"
            assert record.table_rows[18][5] == "2.73E+02"
            assert len(record.table_rows[19]) == 6
            assert record.table_rows[19][0] == "19"
            assert record.table_rows[19][1] == "0.0796905"
            assert record.table_rows[19][2] == "8.78E-02"
            assert record.table_rows[19][3] == "3.64E-01"
            assert record.table_rows[19][4] == "8.21E+02"
            assert record.table_rows[19][5] == "9.87E+02"
            assert len(record.table_rows[20]) == 6
            assert record.table_rows[20][0] == "20"
            assert record.table_rows[20][1] == "-0.084895"
            assert record.table_rows[20][2] == "9.38E-01"
            assert record.table_rows[20][3] == "9.28E-01"
            assert record.table_rows[20][4] == "7.68E+01"
            assert record.table_rows[20][5] == "6.32E+01"
            record = next(records)
            assert record.entity_type == "SAMPLE"
            assert record.entity_id == "Control Embyronic Stem Cell Replicate 2"
            assert len(record.entity_attributes) == 24
            assert record.entity_attributes["Sample_extract_protocol_ch2"] == "TriZol procedure"
            assert record.entity_attributes["Sample_hyb_protocol"] == "Oligoarray control targets and hybridization buffer (Agilent In Situ Hybridization Kit Plus) were added, and samples were applied to microarrays enclosed in Agilent SureHyb-enabled hybridization chambers. After hybridization, slides were washed sequentially with 6x SSC/0.005% Triton X-102 and 0.1x SSC/0.005% Triton X-102 before scanning. Slides were hybridized for 17 h at 60\xb0C in a rotating oven, and washed."
            assert record.entity_attributes["Sample_extract_protocol_ch1"] == "TriZol procedure"
            assert record.entity_attributes["Sample_platform_id"] == "Murine 15K long oligo array version 2.0"
            assert record.entity_attributes["Sample_title"] == "Control Embyronic Stem Cell Replicate 2"
            assert record.entity_attributes["Sample_supplementary_file"] == "file2.gpr"
            assert record.entity_attributes["Sample_organism_ch2"] == "Mus musculus"
            assert record.entity_attributes["Sample_organism_ch1"] == "Mus musculus"
            assert record.entity_attributes["Sample_label_ch1"] == "Cy5"
            assert len(record.entity_attributes["Sample_scan_protocol"]) == 2
            assert record.entity_attributes["Sample_scan_protocol"][0] == "Scanned on an Agilent G2565AA scanner."
            assert record.entity_attributes["Sample_scan_protocol"][1] == "Images were quantified using Agilent Feature Extraction Software (version A.7.5)."
            assert record.entity_attributes["sample_table_begin"] == ""
            assert record.entity_attributes["Sample_label_protocol_ch2"] == "10 \xb5g of total RNA were primed with 2 \xb5l of 100 \xb5M T16N2 DNA primer at 70\xb0C for 10 min, then reversed transcribed at 42\xb0C for 1 h in the presence of 400 U SuperScript II RTase (Invitrogen), and 100 \xb5M each dATP, dTTP, dGTP, with 25 \xb5M dCTP, 25 \xb5M Cy5-labeled dCTP (NEN Life Science, Boston, MA), and RNase inhibitor (Invitrogen). RNA was then degraded with RNase A, and labeled cDNAs were purified using QIAquick PCR columns (Qiagen)."
            assert record.entity_attributes["Sample_label_protocol_ch1"] == "10 \xb5g of total RNA were primed with 2 \xb5l of 100 \xb5M T16N2 DNA primer at 70\xb0C for 10 min, then reversed transcribed at 42\xb0C for 1 h in the presence of 400 U SuperScript II RTase (Invitrogen), and 100 \xb5M each dATP, dTTP, dGTP, with 25 \xb5M dCTP, 25 \xb5M Cy5-labeled dCTP (NEN Life Science, Boston, MA), and RNase inhibitor (Invitrogen). RNA was then degraded with RNase A, and labeled cDNAs were purified using QIAquick PCR columns (Qiagen)."
            assert record.entity_attributes["Sample_data_processing"] == "LOWESS normalized, background subtracted VALUE data obtained from log of processed Red signal/processed Green signal."
            assert record.entity_attributes["sample_table_end"] == ""
            assert record.entity_attributes["Sample_label_ch2"] == "Cy3"
            assert record.entity_attributes["Sample_description"] == "Biological replicate 2 of 4. Control embryonic stem cells, untreated, harvested after several passages."
            assert record.entity_attributes["Sample_source_name_ch1"] == "Total RNA from murine ES-D3 embryonic stem cells labeled with Cyanine-5 (red)."
            assert record.entity_attributes["Sample_source_name_ch2"] == "Total RNA from pooled whole mouse embryos e17.5, labeled with Cyanine-3 (green)."
            assert record.entity_attributes["Sample_molecule_ch1"] == "total RNA"
            assert record.entity_attributes["Sample_molecule_ch2"] == "total RNA"
            assert record.entity_attributes["Sample_growth_protocol_ch1"] == "ES cells were kept in an undifferentiated, pluripotent state by using 1000 IU/ml leukemia inhibitory factor (LIF; Chemicon, ESGRO, ESG1107), and grown on top of murine embryonic fibroblasts feeder layer inactivated by 10 ug/ml of mitomycin C (Sigma, St. Louis). ES cells were cultured on 0.1% gelatin-coated plastic dishes in ES medium containing Dulbecco modified Eagle medium supplemented with 15% fetal calf serum, 0.1 mM beta-mercaptoethanol, 2 mM glutamine, and 0.1 mN non-essential amino acids."
            assert len(record.entity_attributes["Sample_characteristics_ch1"]) == 4
            assert record.entity_attributes["Sample_characteristics_ch1"][0] == "ES-D3 cell line (CRL-1934)"
            assert record.entity_attributes["Sample_characteristics_ch1"][1] == "Age: day 4"
            assert record.entity_attributes["Sample_characteristics_ch1"][2] == "Tissue: blastocytes"
            assert record.entity_attributes["Sample_characteristics_ch1"][3] == "Strain: 129/Sv mice"
            assert len(record.entity_attributes["Sample_characteristics_ch2"]) == 3
            assert record.entity_attributes["Sample_characteristics_ch2"][0] == "Strain: C57BL/6"
            assert record.entity_attributes["Sample_characteristics_ch2"][1] == "Age: e17.5 d"
            assert record.entity_attributes["Sample_characteristics_ch2"][2] == "Tissue: whole embryo"
            assert len(record.col_defs) == 6
            assert record.col_defs["ID_REF"] == ""
            assert record.col_defs["VALUE"] == "log(REDsignal/GREENsignal) per feature (processed signals used)."
            assert record.col_defs["gProcessedSignal"] == 'Dye-normalized signal after surrogate "algorithm," green "channel," used for computation of log ratio.'
            assert record.col_defs["LogRatioError"] == "error of the log ratio calculated according to the error model chosen."
            assert record.col_defs["PValueLogRatio"] == "Significance level of the Log Ratio computed for a feature."
            assert record.col_defs["rProcessedSignal"] == 'Dye-normalized signal after surrogate "algorithm," red "channel," used for computation of log ratio.'
            assert len(record.table_rows) == 21
            assert len(record.table_rows[0]) == 6
            assert record.table_rows[0][0] == "ID_REF"
            assert record.table_rows[0][1] == "VALUE"
            assert record.table_rows[0][2] == "LogRatioError"
            assert record.table_rows[0][3] == "PValueLogRatio"
            assert record.table_rows[0][4] == "gProcessedSignal"
            assert record.table_rows[0][5] == "rProcessedSignal"
            assert len(record.table_rows[1]) == 6
            assert record.table_rows[1][0] == "1"
            assert record.table_rows[1][1] == "-1.1697263"
            assert record.table_rows[1][2] == "1.23E-01"
            assert record.table_rows[1][3] == "2.14E-21"
            assert record.table_rows[1][4] == "3.17E+03"
            assert record.table_rows[1][5] == "2.14E+02"
            assert len(record.table_rows[2]) == 6
            assert record.table_rows[2][0] == "2"
            assert record.table_rows[2][1] == "-0.1111353"
            assert record.table_rows[2][2] == "1.63E+00"
            assert record.table_rows[2][3] == "9.46E-01"
            assert record.table_rows[2][4] == "5.43E+01"
            assert record.table_rows[2][5] == "4.20E+01"
            assert len(record.table_rows[3]) == 6
            assert record.table_rows[3][0] == "3"
            assert record.table_rows[3][1] == "0.1400597"
            assert record.table_rows[3][2] == "5.11E-02"
            assert record.table_rows[3][3] == "6.17E-03"
            assert record.table_rows[3][4] == "6.72E+03"
            assert record.table_rows[3][5] == "9.28E+03"
            assert len(record.table_rows[4]) == 6
            assert record.table_rows[4][0] == "4"
            assert record.table_rows[4][1] == "-0.4820633"
            assert record.table_rows[4][2] == "6.38E-02"
            assert record.table_rows[4][3] == "4.06E-14"
            assert record.table_rows[4][4] == "6.46E+03"
            assert record.table_rows[4][5] == "2.13E+03"
            assert len(record.table_rows[5]) == 6
            assert record.table_rows[5][0] == "5"
            assert record.table_rows[5][1] == "-1.2116196"
            assert record.table_rows[5][2] == "1.22E-01"
            assert record.table_rows[5][3] == "2.31E-23"
            assert record.table_rows[5][4] == "3.62E+03"
            assert record.table_rows[5][5] == "2.22E+02"
            assert len(record.table_rows[6]) == 6
            assert record.table_rows[6][0] == "6"
            assert record.table_rows[6][1] == "-0.0230528"
            assert record.table_rows[6][2] == "1.04E-01"
            assert record.table_rows[6][3] == "8.24E-01"
            assert record.table_rows[6][4] == "8.76E+02"
            assert record.table_rows[6][5] == "8.31E+02"
            assert len(record.table_rows[7]) == 6
            assert record.table_rows[7][0] == "7"
            assert record.table_rows[7][1] == "-1.1380152"
            assert record.table_rows[7][2] == "1.13E-01"
            assert record.table_rows[7][3] == "9.23E-24"
            assert record.table_rows[7][4] == "3.94E+03"
            assert record.table_rows[7][5] == "2.86E+02"
            assert len(record.table_rows[8]) == 6
            assert record.table_rows[8][0] == "8"
            assert record.table_rows[8][1] == "-1.834596"
            assert record.table_rows[8][2] == "5.40E-01"
            assert record.table_rows[8][3] == "6.74E-04"
            assert record.table_rows[8][4] == "6.44E+02"
            assert record.table_rows[8][5] == "9.43E+00"
            assert len(record.table_rows[9]) == 6
            assert record.table_rows[9][0] == "9"
            assert record.table_rows[9][1] == "-0.9747637"
            assert record.table_rows[9][2] == "2.14E+00"
            assert record.table_rows[9][3] == "1.00E+00"
            assert record.table_rows[9][4] == "9.17E+01"
            assert record.table_rows[9][5] == "9.72E+00"
            assert len(record.table_rows[10]) == 6
            assert record.table_rows[10][0] == "10"
            assert record.table_rows[10][1] == "0.3874005"
            assert record.table_rows[10][2] == "2.92E-01"
            assert record.table_rows[10][3] == "1.85E-01"
            assert record.table_rows[10][4] == "1.69E+02"
            assert record.table_rows[10][5] == "4.11E+02"
            assert len(record.table_rows[11]) == 6
            assert record.table_rows[11][0] == "11"
            assert record.table_rows[11][1] == "0.5340442"
            assert record.table_rows[11][2] == "3.29E-01"
            assert record.table_rows[11][3] == "1.04E-01"
            assert record.table_rows[11][4] == "1.23E+02"
            assert record.table_rows[11][5] == "4.20E+02"
            assert len(record.table_rows[12]) == 6
            assert record.table_rows[12][0] == "12"
            assert record.table_rows[12][1] == "0.3260696"
            assert record.table_rows[12][2] == "1.92E+00"
            assert record.table_rows[12][3] == "8.65E-01"
            assert record.table_rows[12][4] == "2.73E+01"
            assert record.table_rows[12][5] == "5.77E+01"
            assert len(record.table_rows[13]) == 6
            assert record.table_rows[13][0] == "13"
            assert record.table_rows[13][1] == "0.3010618"
            assert record.table_rows[13][2] == "2.84E-01"
            assert record.table_rows[13][3] == "2.90E-01"
            assert record.table_rows[13][4] == "1.93E+02"
            assert record.table_rows[13][5] == "3.87E+02"
            assert len(record.table_rows[14]) == 6
            assert record.table_rows[14][0] == "14"
            assert record.table_rows[14][1] == "-1.0760413"
            assert record.table_rows[14][2] == "1.08E-01"
            assert record.table_rows[14][3] == "1.63E-23"
            assert record.table_rows[14][4] == "4.06E+03"
            assert record.table_rows[14][5] == "3.41E+02"
            assert len(record.table_rows[15]) == 6
            assert record.table_rows[15][0] == "15"
            assert record.table_rows[15][1] == "-0.1167371"
            assert record.table_rows[15][2] == "3.87E-01"
            assert record.table_rows[15][3] == "7.63E-01"
            assert record.table_rows[15][4] == "2.32E+02"
            assert record.table_rows[15][5] == "1.77E+02"
            assert len(record.table_rows[16]) == 6
            assert record.table_rows[16][0] == "16"
            assert record.table_rows[16][1] == "-0.1936322"
            assert record.table_rows[16][2] == "9.44E-01"
            assert record.table_rows[16][3] == "8.38E-01"
            assert record.table_rows[16][4] == "1.02E+02"
            assert record.table_rows[16][5] == "6.56E+01"
            assert len(record.table_rows[17]) == 6
            assert record.table_rows[17][0] == "17"
            assert record.table_rows[17][1] == "-0.3275898"
            assert record.table_rows[17][2] == "7.87E-01"
            assert record.table_rows[17][3] == "6.77E-01"
            assert record.table_rows[17][4] == "1.41E+02"
            assert record.table_rows[17][5] == "6.65E+01"
            assert len(record.table_rows[18]) == 6
            assert record.table_rows[18][0] == "18"
            assert record.table_rows[18][1] == "-0.4805853"
            assert record.table_rows[18][2] == "1.14E-01"
            assert record.table_rows[18][3] == "2.41E-05"
            assert record.table_rows[18][4] == "1.34E+03"
            assert record.table_rows[18][5] == "4.42E+02"
            assert len(record.table_rows[19]) == 6
            assert record.table_rows[19][0] == "19"
            assert record.table_rows[19][1] == "0.1109524"
            assert record.table_rows[19][2] == "9.56E-02"
            assert record.table_rows[19][3] == "2.46E-01"
            assert record.table_rows[19][4] == "8.38E+02"
            assert record.table_rows[19][5] == "1.08E+03"
            assert len(record.table_rows[20]) == 6
            assert record.table_rows[20][0] == "20"
            assert record.table_rows[20][1] == "0.1677912"
            assert record.table_rows[20][2] == "6.51E-01"
            assert record.table_rows[20][3] == "7.97E-01"
            assert record.table_rows[20][4] == "9.84E+01"
            assert record.table_rows[20][5] == "1.45E+02"
            record = next(records)
            assert record.entity_type == "SAMPLE"
            assert record.entity_id == "Triple-Fusion Transfected Embryonic Stem Cells Replicate 1"
            assert len(record.entity_attributes) == 25
            assert record.entity_attributes["Sample_extract_protocol_ch2"] == "TriZol procedure"
            assert record.entity_attributes["Sample_hyb_protocol"] == "Oligoarray control targets and hybridization buffer (Agilent In Situ Hybridization Kit Plus) were added, and samples were applied to microarrays enclosed in Agilent SureHyb-enabled hybridization chambers. After hybridization, slides were washed sequentially with 6x SSC/0.005% Triton X-102 and 0.1x SSC/0.005% Triton X-102 before scanning. Slides were hybridized for 17 h at 60\xb0C in a rotating oven, and washed."
            assert record.entity_attributes["Sample_extract_protocol_ch1"] == "TriZol procedure"
            assert record.entity_attributes["Sample_platform_id"] == "Murine 15K long oligo array version 2.0"
            assert record.entity_attributes["Sample_title"] == "Triple-Fusion Transfected Embryonic Stem Cells Replicate 1"
            assert record.entity_attributes["Sample_supplementary_file"] == "file3.gpr"
            assert record.entity_attributes["Sample_organism_ch2"] == "Mus musculus"
            assert record.entity_attributes["Sample_organism_ch1"] == "Mus musculus"
            assert record.entity_attributes["Sample_label_ch1"] == "Cy5"
            assert record.entity_attributes["Sample_scan_protocol"] == "Scanned on an Agilent G2565AA scanner."
            assert record.entity_attributes["sample_table_begin"] == ""
            assert record.entity_attributes["Sample_label_protocol_ch2"] == "10 \xb5g of total RNA were primed with 2 \xb5l of 100 \xb5M T16N2 DNA primer at 70\xb0C for 10 min, then reversed transcribed at 42\xb0C for 1 h in the presence of 400 U SuperScript II RTase (Invitrogen), and 100 \xb5M each dATP, dTTP, dGTP, with 25 \xb5M dCTP, 25 \xb5M Cy5-labeled dCTP (NEN Life Science, Boston, MA), and RNase inhibitor (Invitrogen). RNA was then degraded with RNase A, and labeled cDNAs were purified using QIAquick PCR columns (Qiagen)."
            assert record.entity_attributes["Sample_label_protocol_ch1"] == "10 \xb5g of total RNA were primed with 2 \xb5l of 100 \xb5M T16N2 DNA primer at 70\xb0C for 10 min, then reversed transcribed at 42\xb0C for 1 h in the presence of 400 U SuperScript II RTase (Invitrogen), and 100 \xb5M each dATP, dTTP, dGTP with 25 \xb5M dCTP, 25 \xb5M Cy5-labeled dCTP (NEN Life Science, Boston, MA), and RNase inhibitor (Invitrogen). RNA was then degraded with RNase A, and labeled cDNAs were purified using QIAquick PCR columns (Qiagen)."
            assert record.entity_attributes["Sample_data_processing"] == "LOWESS normalized, background subtracted VALUE data obtained from log of processed Red signal/processed Green signal."
            assert record.entity_attributes["sample_table_end"] == ""
            assert record.entity_attributes["Sample_label_ch2"] == "Cy3"
            assert record.entity_attributes["Sample_description"] == "Biological replicate 1 of 3. Stable triple-fusion-reporter-gene transfected embryonic stem cells, harvested after several passages."
            assert record.entity_attributes["Sample_source_name_ch1"] == "Total RNA from murine ES-D3 triple-transfected embryonic stem cells labeled with Cyanine-5 (red)."
            assert record.entity_attributes["Sample_source_name_ch2"] == "Total RNA from pooled whole mouse embryos e17.5, labeled with Cyanine-3 (green)."
            assert record.entity_attributes["Sample_molecule_ch1"] == "total RNA"
            assert record.entity_attributes["Sample_molecule_ch2"] == "total RNA"
            assert record.entity_attributes["Sample_growth_protocol_ch1"] == "ES cells were kept in an undifferentiated, pluripotent state by using 1000 IU/ml leukemia inhibitory factor (LIF; Chemicon, ESGRO, ESG1107), and grown on top of murine embryonic fibroblasts feeder layer inactivated by 10 ug/ml of mitomycin C (Sigma, St. Louis). ES cells were cultured on 0.1% gelatin-coated plastic dishes in ES medium containing Dulbecco modified Eagle medium supplemented with 15% fetal calf serum, 0.1 mM beta-mercaptoethanol, 2 mM glutamine, and 0.1 mN non-essential amino acids."
            assert len(record.entity_attributes["Sample_characteristics_ch1"]) == 5
            assert record.entity_attributes["Sample_characteristics_ch1"][0] == "ES-D3 cell line (CRL-1934)"
            assert record.entity_attributes["Sample_characteristics_ch1"][1] == "Transfected with pUb-fluc-mrfp-ttk triple fusion reporter gene."
            assert record.entity_attributes["Sample_characteristics_ch1"][2] == "Age: day 4"
            assert record.entity_attributes["Sample_characteristics_ch1"][3] == "Tissue: blastocytes"
            assert record.entity_attributes["Sample_characteristics_ch1"][4] == "Strain: 129/Sv mice"
            assert len(record.entity_attributes["Sample_characteristics_ch2"]) == 3
            assert record.entity_attributes["Sample_characteristics_ch2"][0] == "Strain: C57BL/6"
            assert record.entity_attributes["Sample_characteristics_ch2"][1] == "Age: e17.5 d"
            assert record.entity_attributes["Sample_characteristics_ch2"][2] == "Tissue: whole embryo"
            assert record.entity_attributes["Sample_treatment_protocol_ch1"] == "PCR amplification and standard cloning techniques were used to insert fluc and mrfp genes from plasmids pCDNA 3.1-CMV-fluc (Promega, Madison, WI) and pCDNA3.1-CMV-mrfp in frame with the ttk gene into the pCDNA3.1-truncated	sr39tk. This triple fusion (TF) reporter gene fragment (3.3 kbp) was released from the plasmid with Not1 and BamH1 restriction enzymes before blunt-end ligation into the multiple cloning site of lentiviral transfer vector, FUG, driven by the human ubiquitin-C promoter. Self-inactivating (SIN) lentivirus was prepared by transient transfection of 293T cells. Briefly, pFUG-TF containing the triple fusion reporter gene was co-transfected into 293T cells with HIV-1 packaging vector (?8.9) and vesicular stomatitis virus G glycoprotein-pseudotyped envelop vector (pVSVG). Lentivirus supernatant was concentrated by sediment centrifugation using a SW29 rotor at 50,000 x g for two hours. Concentrated virus was titered on 293T cells. Murine ES cells were transfected with LV-pUb-fluc-mrfp-ttk at a multiplicity of infection (MOI) of 10."
            assert len(record.col_defs) == 6
            assert record.col_defs["ID_REF"] == ""
            assert record.col_defs["VALUE"] == "log(REDsignal/GREENsignal) per feature (processed signals used)."
            assert record.col_defs["gProcessedSignal"] == 'Dye-normalized signal after surrogate "algorithm," green "channel," used for computation of log ratio.'
            assert record.col_defs["LogRatioError"] == "error of the log ratio calculated according to the error model chosen."
            assert record.col_defs["PValueLogRatio"] == "Significance level of the Log Ratio computed for a feature."
            assert record.col_defs["rProcessedSignal"] == 'Dye-normalized signal after surrogate "algorithm," red "channel," used for computation of log ratio.'
            assert len(record.table_rows) == 21
            assert len(record.table_rows[0]) == 6
            assert record.table_rows[0][0] == "ID_REF"
            assert record.table_rows[0][1] == "VALUE"
            assert record.table_rows[0][2] == "LogRatioError"
            assert record.table_rows[0][3] == "PValueLogRatio"
            assert record.table_rows[0][4] == "gProcessedSignal"
            assert record.table_rows[0][5] == "rProcessedSignal"
            assert len(record.table_rows[1]) == 6
            assert record.table_rows[1][0] == "1"
            assert record.table_rows[1][1] == "-0.7837546"
            assert record.table_rows[1][2] == "1.30E-01"
            assert record.table_rows[1][3] == "1.70E-09"
            assert record.table_rows[1][4] == "2.10E+03"
            assert record.table_rows[1][5] == "3.46E+02"
            assert len(record.table_rows[2]) == 6
            assert record.table_rows[2][0] == "2"
            assert record.table_rows[2][1] == "0.3797837"
            assert record.table_rows[2][2] == "1.15E+00"
            assert record.table_rows[2][3] == "7.41E-01"
            assert record.table_rows[2][4] == "5.59E+01"
            assert record.table_rows[2][5] == "1.34E+02"
            assert len(record.table_rows[3]) == 6
            assert record.table_rows[3][0] == "3"
            assert record.table_rows[3][1] == "0.2079269"
            assert record.table_rows[3][2] == "5.38E-02"
            assert record.table_rows[3][3] == "1.12E-04"
            assert record.table_rows[3][4] == "5.04E+03"
            assert record.table_rows[3][5] == "8.14E+03"
            assert len(record.table_rows[4]) == 6
            assert record.table_rows[4][0] == "4"
            assert record.table_rows[4][1] == "-0.4730291"
            assert record.table_rows[4][2] == "6.71E-02"
            assert record.table_rows[4][3] == "1.86E-12"
            assert record.table_rows[4][4] == "5.66E+03"
            assert record.table_rows[4][5] == "1.91E+03"
            assert len(record.table_rows[5]) == 6
            assert record.table_rows[5][0] == "5"
            assert record.table_rows[5][1] == "-0.9481128"
            assert record.table_rows[5][2] == "1.19E-01"
            assert record.table_rows[5][3] == "1.30E-15"
            assert record.table_rows[5][4] == "3.10E+03"
            assert record.table_rows[5][5] == "3.49E+02"
            assert len(record.table_rows[6]) == 6
            assert record.table_rows[6][0] == "6"
            assert record.table_rows[6][1] == "-0.0159867"
            assert record.table_rows[6][2] == "1.33E-01"
            assert record.table_rows[6][3] == "9.05E-01"
            assert record.table_rows[6][4] == "8.45E+02"
            assert record.table_rows[6][5] == "8.14E+02"
            assert len(record.table_rows[7]) == 6
            assert record.table_rows[7][0] == "7"
            assert record.table_rows[7][1] == "-0.819922"
            assert record.table_rows[7][2] == "1.14E-01"
            assert record.table_rows[7][3] == "7.01E-13"
            assert record.table_rows[7][4] == "2.75E+03"
            assert record.table_rows[7][5] == "4.16E+02"
            assert len(record.table_rows[8]) == 6
            assert record.table_rows[8][0] == "8"
            assert record.table_rows[8][1] == "-0.1559774"
            assert record.table_rows[8][2] == "9.16E-01"
            assert record.table_rows[8][3] == "8.65E-01"
            assert record.table_rows[8][4] == "1.34E+02"
            assert record.table_rows[8][5] == "9.34E+01"
            assert len(record.table_rows[9]) == 6
            assert record.table_rows[9][0] == "9"
            assert record.table_rows[9][1] == "0.145267"
            assert record.table_rows[9][2] == "3.90E+00"
            assert record.table_rows[9][3] == "1.00E+00"
            assert record.table_rows[9][4] == "2.22E+01"
            assert record.table_rows[9][5] == "3.10E+01"
            assert len(record.table_rows[10]) == 6
            assert record.table_rows[10][0] == "10"
            assert record.table_rows[10][1] == "0.3611211"
            assert record.table_rows[10][2] == "3.40E-01"
            assert record.table_rows[10][3] == "2.88E-01"
            assert record.table_rows[10][4] == "1.97E+02"
            assert record.table_rows[10][5] == "4.52E+02"
            assert len(record.table_rows[11]) == 6
            assert record.table_rows[11][0] == "11"
            assert record.table_rows[11][1] == "0.5092089"
            assert record.table_rows[11][2] == "4.39E-01"
            assert record.table_rows[11][3] == "2.46E-01"
            assert record.table_rows[11][4] == "1.24E+02"
            assert record.table_rows[11][5] == "4.01E+02"
            assert len(record.table_rows[12]) == 6
            assert record.table_rows[12][0] == "12"
            assert record.table_rows[12][1] == "0.3715387"
            assert record.table_rows[12][2] == "1.69E+00"
            assert record.table_rows[12][3] == "8.26E-01"
            assert record.table_rows[12][4] == "3.84E+01"
            assert record.table_rows[12][5] == "9.04E+01"
            assert len(record.table_rows[13]) == 6
            assert record.table_rows[13][0] == "13"
            assert record.table_rows[13][1] == "0.1734934"
            assert record.table_rows[13][2] == "3.57E-01"
            assert record.table_rows[13][3] == "6.27E-01"
            assert record.table_rows[13][4] == "2.37E+02"
            assert record.table_rows[13][5] == "3.53E+02"
            assert len(record.table_rows[14]) == 6
            assert record.table_rows[14][0] == "14"
            assert record.table_rows[14][1] == "-0.9340707"
            assert record.table_rows[14][2] == "1.20E-01"
            assert record.table_rows[14][3] == "6.90E-15"
            assert record.table_rows[14][4] == "2.96E+03"
            assert record.table_rows[14][5] == "3.45E+02"
            assert len(record.table_rows[15]) == 6
            assert record.table_rows[15][0] == "15"
            assert record.table_rows[15][1] == "-0.2956317"
            assert record.table_rows[15][2] == "5.78E-01"
            assert record.table_rows[15][3] == "6.09E-01"
            assert record.table_rows[15][4] == "2.46E+02"
            assert record.table_rows[15][5] == "1.25E+02"
            assert len(record.table_rows[16]) == 6
            assert record.table_rows[16][0] == "16"
            assert record.table_rows[16][1] == "-0.2321102"
            assert record.table_rows[16][2] == "1.22E+00"
            assert record.table_rows[16][3] == "8.49E-01"
            assert record.table_rows[16][4] == "1.09E+02"
            assert record.table_rows[16][5] == "6.37E+01"
            assert len(record.table_rows[17]) == 6
            assert record.table_rows[17][0] == "17"
            assert record.table_rows[17][1] == "-0.1603561"
            assert record.table_rows[17][2] == "1.16E+00"
            assert record.table_rows[17][3] == "8.90E-01"
            assert record.table_rows[17][4] == "1.06E+02"
            assert record.table_rows[17][5] == "7.34E+01"
            assert len(record.table_rows[18]) == 6
            assert record.table_rows[18][0] == "18"
            assert record.table_rows[18][1] == "-0.5063897"
            assert record.table_rows[18][2] == "1.63E-01"
            assert record.table_rows[18][3] == "1.95E-03"
            assert record.table_rows[18][4] == "1.15E+03"
            assert record.table_rows[18][5] == "3.58E+02"
            assert len(record.table_rows[19]) == 6
            assert record.table_rows[19][0] == "19"
            assert record.table_rows[19][1] == "0.1990761"
            assert record.table_rows[19][2] == "1.32E-01"
            assert record.table_rows[19][3] == "1.32E-01"
            assert record.table_rows[19][4] == "6.65E+02"
            assert record.table_rows[19][5] == "1.05E+03"
            assert len(record.table_rows[20]) == 6
            assert record.table_rows[20][0] == "20"
            assert record.table_rows[20][1] == "0.2985912"
            assert record.table_rows[20][2] == "8.89E-01"
            assert record.table_rows[20][3] == "7.37E-01"
            assert record.table_rows[20][4] == "8.06E+01"
            assert record.table_rows[20][5] == "1.60E+02"
            record = next(records)
            assert record.entity_type == "SERIES"
            assert record.entity_id == "Murine ES Cells"
            assert len(record.entity_attributes) == 7
            assert len(record.entity_attributes["Series_sample_id"]) == 3
            assert record.entity_attributes["Series_sample_id"][0] == "Control Embyronic Stem Cell Replicate 1"
            assert record.entity_attributes["Series_sample_id"][1] == "Control Embyronic Stem Cell Replicate 2"
            assert record.entity_attributes["Series_sample_id"][2] == "Triple-Fusion Transfected Embryonic Stem Cells Replicate 1"
            assert record.entity_attributes["Series_pubmed_id"] == "16390873"
            assert len(record.entity_attributes["Series_contributor"]) == 9
            assert record.entity_attributes["Series_contributor"][0] == "Joseph,C,Wu"
            assert record.entity_attributes["Series_contributor"][1] == "Joshua,M,Spin"
            assert record.entity_attributes["Series_contributor"][2] == "Feng,,Cao"
            assert record.entity_attributes["Series_contributor"][3] == "Shaun,,Lin"
            assert record.entity_attributes["Series_contributor"][4] == "Olivier,,Gheysens"
            assert record.entity_attributes["Series_contributor"][5] == "Ian,Y,Chen"
            assert record.entity_attributes["Series_contributor"][6] == "Anya,,Tsalenko"
            assert record.entity_attributes["Series_contributor"][7] == "Sanjiv,S,Ghambhir"
            assert record.entity_attributes["Series_contributor"][8] == "Thomas,,Quertermous"
            assert record.entity_attributes["Series_summary"] == "Transcriptional profiling of mouse embryonic stem cells comparing control untreated ES cells with ES cells transfected with a pUb-fluc-mrfp-ttk triple fusion reporter gene. The latter makes ES visualization possible by FACS and single ce"
            assert record.entity_attributes["Series_type"] == "Genetic modification"
            assert record.entity_attributes["Series_title"] == "Murine ES Cells: Control vs. Triple-Fusion Transfected"
            assert record.entity_attributes["Series_overall_design"] == "Two-condition experiment, ES vs. TF-ES cells. Biological replicates: 4 control, 3 transfected, independently grown and harvested. One replicate per array."
            assert len(record.col_defs) == 0
            assert len(record.table_rows) == 0

    def test_GSM804(self):
        path = "Geo/GSM804.txt"
        with open(path, encoding="latin") as handle:
            records = Geo.parse(handle)
            record = next(records)
            assert record.entity_type == "SAMPLE"
            assert record.entity_id == "GSM804"
            assert len(record.entity_attributes) == 18
            assert record.entity_attributes["Sample_pubmed_id"] == "11687795"
            assert record.entity_attributes["Sample_submitter_institute"] == "University of California San Francisco"
            assert len(record.entity_attributes["Sample_author"]) == 19
            assert record.entity_attributes["Sample_author"][0] == "Antoine,M,Snijders"
            assert record.entity_attributes["Sample_author"][1] == "Norma,,Nowak"
            assert record.entity_attributes["Sample_author"][2] == "Richard,,Segraves"
            assert record.entity_attributes["Sample_author"][3] == "Stephanie,,Blackwood"
            assert record.entity_attributes["Sample_author"][4] == "Nils,,Brown"
            assert record.entity_attributes["Sample_author"][5] == "Jeffery,,Conroy"
            assert record.entity_attributes["Sample_author"][6] == "Greg,,Hamilton"
            assert record.entity_attributes["Sample_author"][7] == "Anna,K,Hindle"
            assert record.entity_attributes["Sample_author"][8] == "Bing,,Huey"
            assert record.entity_attributes["Sample_author"][9] == "Karen,,Kimura"
            assert record.entity_attributes["Sample_author"][10] == "Sindy,,Law"
            assert record.entity_attributes["Sample_author"][11] == "Ken,,Myambo"
            assert record.entity_attributes["Sample_author"][12] == "Joel,,Palmer"
            assert record.entity_attributes["Sample_author"][13] == "Bauke,,Ylstra"
            assert record.entity_attributes["Sample_author"][14] == "Jingzhu,P,Yue"
            assert record.entity_attributes["Sample_author"][15] == "Joe,W,Gray"
            assert record.entity_attributes["Sample_author"][16] == "Ajay,N,Jain"
            assert record.entity_attributes["Sample_author"][17] == "Daniel,,Pinkel"
            assert record.entity_attributes["Sample_author"][18] == "Donna,G,Albertson"
            assert record.entity_attributes["Sample_submitter_phone"] == "415 502-8463"
            assert record.entity_attributes["Sample_submitter_department"] == "Comprehensive Cancer Center"
            assert len(record.entity_attributes["Sample_description"]) == 4
            assert record.entity_attributes["Sample_description"][0] == 'Coriell Cell Repositories cell line <a href="http://locus.umdnj.edu/nigms/nigms_cgi/display.cgi?GM05296">GM05296</a>.'
            assert record.entity_attributes["Sample_description"][1] == "Fibroblast cell line derived from a 1 month old female with multiple congenital malformations, dysmorphic features, intrauterine growth retardation, heart murmur, cleft palate, equinovarus deformity, microcephaly, coloboma of right iris, clinodactyly, reduced RBC catalase activity, and 1 copy of catalase gene."
            assert record.entity_attributes["Sample_description"][2] == "Chromosome abnormalities are present."
            assert record.entity_attributes["Sample_description"][3] == "Karyotype is 46,XX,-11,+der(11)inv ins(11;10)(11pter> 11p13::10q21>10q24::11p13>11qter)mat"
            assert record.entity_attributes["Sample_target_source2"] == "normal male reference genomic DNA"
            assert record.entity_attributes["Sample_target_source1"] == "Cell line GM05296"
            assert record.entity_attributes["Sample_submitter_name"] == "Donna,G,Albertson"
            assert record.entity_attributes["Sample_platform_id"] == "GPL28"
            assert record.entity_attributes["Sample_type"] == "dual channel genomic"
            assert record.entity_attributes["Sample_status"] == "Public on Feb 12 2002"
            assert record.entity_attributes["Sample_submitter_email"] == "albertson@cc.ucsf.edu"
            assert record.entity_attributes["Sample_title"] == "CGH_Albertson_GM05296-001218"
            assert record.entity_attributes["Sample_organism"] == "Homo sapiens"
            assert record.entity_attributes["Sample_series_id"] == "GSE16"
            assert record.entity_attributes["Sample_submission_date"] == "Jan 17 2002"
            assert record.entity_attributes["Sample_submitter_city"] == "San Francisco,CA,94143,USA"
            assert len(record.col_defs) == 5
            assert record.col_defs["NO_REPLICATES"] == "Number of replicate spot measurements"
            assert record.col_defs["LOG2STDDEV"] == "Standard deviation of VALUE"
            assert record.col_defs["ID_REF"] == "Unique row identifier, genome position order"
            assert record.col_defs["VALUE"] == "aka LOG2RATIO, mean of log base 2 of LINEAR_RATIO"
            assert record.col_defs["LINEAR_RATIO"] == "Mean of replicate Cy3/Cy5 ratios"
            assert len(record.table_rows) == 21
            assert len(record.table_rows[0]) == 5
            assert record.table_rows[0][0] == "ID_REF"
            assert record.table_rows[0][1] == "VALUE"
            assert record.table_rows[0][2] == "LINEAR_RATIO"
            assert record.table_rows[0][3] == "LOG2STDDEV"
            assert record.table_rows[0][4] == "NO_REPLICATES"
            assert len(record.table_rows[1]) == 5
            assert record.table_rows[1][0] == "1"
            assert record.table_rows[1][1] == ""
            assert record.table_rows[1][2] == "1.047765"
            assert record.table_rows[1][3] == "0.011853"
            assert record.table_rows[1][4] == "3"
            assert len(record.table_rows[2]) == 5
            assert record.table_rows[2][0] == "2"
            assert record.table_rows[2][1] == ""
            assert record.table_rows[2][2] == ""
            assert record.table_rows[2][3] == ""
            assert record.table_rows[2][4] == "0"
            assert len(record.table_rows[3]) == 5
            assert record.table_rows[3][0] == "3"
            assert record.table_rows[3][1] == "0.008824"
            assert record.table_rows[3][2] == "1.006135"
            assert record.table_rows[3][3] == "0.00143"
            assert record.table_rows[3][4] == "3"
            assert len(record.table_rows[4]) == 5
            assert record.table_rows[4][0] == "4"
            assert record.table_rows[4][1] == "-0.000894"
            assert record.table_rows[4][2] == "0.99938"
            assert record.table_rows[4][3] == "0.001454"
            assert record.table_rows[4][4] == "3"
            assert len(record.table_rows[5]) == 5
            assert record.table_rows[5][0] == "5"
            assert record.table_rows[5][1] == "0.075875"
            assert record.table_rows[5][2] == "1.054"
            assert record.table_rows[5][3] == "0.003077"
            assert record.table_rows[5][4] == "3"
            assert len(record.table_rows[6]) == 5
            assert record.table_rows[6][0] == "6"
            assert record.table_rows[6][1] == "0.017303"
            assert record.table_rows[6][2] == "1.012066"
            assert record.table_rows[6][3] == "0.005876"
            assert record.table_rows[6][4] == "2"
            assert len(record.table_rows[7]) == 5
            assert record.table_rows[7][0] == "7"
            assert record.table_rows[7][1] == "-0.006766"
            assert record.table_rows[7][2] == "0.995321"
            assert record.table_rows[7][3] == "0.013881"
            assert record.table_rows[7][4] == "3"
            assert len(record.table_rows[8]) == 5
            assert record.table_rows[8][0] == "8"
            assert record.table_rows[8][1] == "0.020755"
            assert record.table_rows[8][2] == "1.014491"
            assert record.table_rows[8][3] == "0.005506"
            assert record.table_rows[8][4] == "3"
            assert len(record.table_rows[9]) == 5
            assert record.table_rows[9][0] == "9"
            assert record.table_rows[9][1] == "-0.094938"
            assert record.table_rows[9][2] == "0.936313"
            assert record.table_rows[9][3] == "0.012662"
            assert record.table_rows[9][4] == "3"
            assert len(record.table_rows[10]) == 5
            assert record.table_rows[10][0] == "10"
            assert record.table_rows[10][1] == "-0.054527"
            assert record.table_rows[10][2] == "0.96291"
            assert record.table_rows[10][3] == "0.01073"
            assert record.table_rows[10][4] == "3"
            assert len(record.table_rows[11]) == 5
            assert record.table_rows[11][0] == "11"
            assert record.table_rows[11][1] == "-0.025057"
            assert record.table_rows[11][2] == "0.982782"
            assert record.table_rows[11][3] == "0.003855"
            assert record.table_rows[11][4] == "3"
            assert len(record.table_rows[12]) == 5
            assert record.table_rows[12][0] == "12"
            assert record.table_rows[12][1] == ""
            assert record.table_rows[12][2] == ""
            assert record.table_rows[12][3] == ""
            assert record.table_rows[12][4] == "0"
            assert len(record.table_rows[13]) == 5
            assert record.table_rows[13][0] == "13"
            assert record.table_rows[13][1] == "0.108454"
            assert record.table_rows[13][2] == "1.078072"
            assert record.table_rows[13][3] == "0.005196"
            assert record.table_rows[13][4] == "3"
            assert len(record.table_rows[14]) == 5
            assert record.table_rows[14][0] == "14"
            assert record.table_rows[14][1] == "0.078633"
            assert record.table_rows[14][2] == "1.056017"
            assert record.table_rows[14][3] == "0.009165"
            assert record.table_rows[14][4] == "3"
            assert len(record.table_rows[15]) == 5
            assert record.table_rows[15][0] == "15"
            assert record.table_rows[15][1] == "0.098571"
            assert record.table_rows[15][2] == "1.070712"
            assert record.table_rows[15][3] == "0.007834"
            assert record.table_rows[15][4] == "3"
            assert len(record.table_rows[16]) == 5
            assert record.table_rows[16][0] == "16"
            assert record.table_rows[16][1] == "0.044048"
            assert record.table_rows[16][2] == "1.031003"
            assert record.table_rows[16][3] == "0.013651"
            assert record.table_rows[16][4] == "3"
            assert len(record.table_rows[17]) == 5
            assert record.table_rows[17][0] == "17"
            assert record.table_rows[17][1] == "0.018039"
            assert record.table_rows[17][2] == "1.012582"
            assert record.table_rows[17][3] == "0.005471"
            assert record.table_rows[17][4] == "3"
            assert len(record.table_rows[18]) == 5
            assert record.table_rows[18][0] == "18"
            assert record.table_rows[18][1] == "-0.088807"
            assert record.table_rows[18][2] == "0.9403"
            assert record.table_rows[18][3] == "0.010571"
            assert record.table_rows[18][4] == "3"
            assert len(record.table_rows[19]) == 5
            assert record.table_rows[19][0] == "19"
            assert record.table_rows[19][1] == "0.016349"
            assert record.table_rows[19][2] == "1.011397"
            assert record.table_rows[19][3] == "0.007113"
            assert record.table_rows[19][4] == "3"
            assert len(record.table_rows[20]) == 5
            assert record.table_rows[20][0] == "20"
            assert record.table_rows[20][1] == "0.030977"
            assert record.table_rows[20][2] == "1.021704"
            assert record.table_rows[20][3] == "0.016798"
            assert record.table_rows[20][4] == "3"

    def test_soft_ex_affy_chp(self):
        path = "Geo/soft_ex_affy_chp.txt"
        with open(path, encoding="latin") as handle:
            records = Geo.parse(handle)
            record = next(records)
            assert record.entity_type == "SAMPLE"
            assert record.entity_id == "Drosophila_T0-1"
            assert len(record.entity_attributes) == 16
            assert record.entity_attributes["Sample_organism_ch1"] == "Drosophila melanogaster"
            assert record.entity_attributes["Sample_label_ch1"] == "biotin"
            assert record.entity_attributes["Sample_description"] == "Gene expression data from embryos younger than nuclear cycle 9, i.e. before zygotic genome activation."
            assert record.entity_attributes["Sample_growth_protocol_ch1"] == "30 min egg collections of OreR and yw flies at 25C were aged at room temperature (RT) according to the different temporal classes T0-T4."
            assert len(record.entity_attributes["Sample_characteristics_ch1"]) == 2
            assert record.entity_attributes["Sample_characteristics_ch1"][0] == "Genotype: yellow white and Oregon R parents"
            assert record.entity_attributes["Sample_characteristics_ch1"][1] == "Age: embryos younger than nuclear cycle 9, i.e. before pole cells budding"
            assert record.entity_attributes["Sample_scan_protocol"] == "GeneChips were scanned using the Hewlett-Packard GeneArray Scanner G2500A."
            assert record.entity_attributes["Sample_hyb_protocol"] == "Following fragmentation, 10 microg of cRNA were hybridized for 16 hr at 45C on GeneChip Drosophila Genome Array. GeneChips were washed and stained in the Affymetrix Fluidics Station 400."
            assert record.entity_attributes["Sample_extract_protocol_ch1"] == "Trizol extraction of total RNA was performed according to the manufacturer's instructions."
            assert record.entity_attributes["Sample_source_name_ch1"] == "Drosophila embryos before nuclear cycle 9 (maternal transcripts)"
            assert record.entity_attributes["Sample_table"] == "Drosophila_T0-1.CHP"
            assert record.entity_attributes["Sample_molecule_ch1"] == "total RNA"
            assert record.entity_attributes["Sample_label_protocol_ch1"] == "Biotinylated cRNA were prepared according to the standard Affymetrix protocol from 6 microg total RNA (Expression Analysis Technical Manual, 2001, Affymetrix)."
            assert record.entity_attributes["Sample_data_processing"] == "The data were analyzed with Microarray Suite version 5.0 (MAS 5.0) using Affymetrix default analysis settings and global scaling as normalization method. The trimmed mean target intensity of each array was arbitrarily set to 100."
            assert record.entity_attributes["Sample_treatment_protocol_ch1"] == "Embryos were dechorionated with 50% bleach, put on a cover slip and covered with Halocarbon oil 27 (Sigma). Embryos of the appropriate stage were manually selected under the dissecting scope. Selected embryos were transferred to a basket, rinsed with PBS with 0,7% NaCl, 0,04% triton-X100 and placed on ice in the Trizol solution (GibcoBRL)."
            assert record.entity_attributes["Sample_title"] == "embryo at T0, biological rep1"
            assert record.entity_attributes["Sample_supplementary_file"] == "Drosophila_T0-1.CEL"
            assert len(record.col_defs) == 0
            assert len(record.table_rows) == 0
            record = next(records)
            assert record.entity_type == "SAMPLE"
            assert record.entity_id == "Drosophila_T0-2"
            assert len(record.entity_attributes) == 16
            assert record.entity_attributes["Sample_organism_ch1"] == "Drosophila melanogaster"
            assert record.entity_attributes["Sample_label_ch1"] == "biotin"
            assert record.entity_attributes["Sample_description"] == "Gene expression data from embryos younger than nuclear cycle 9, i.e. before zygotic genome activation."
            assert record.entity_attributes["Sample_growth_protocol_ch1"] == "30 min egg collections of OreR and yw flies at 25C were aged at room temperature (RT) according to the different temporal classes T0-T4."
            assert len(record.entity_attributes["Sample_characteristics_ch1"]) == 2
            assert record.entity_attributes["Sample_characteristics_ch1"][0] == "Genotype: yellow white and Oregon R parents"
            assert record.entity_attributes["Sample_characteristics_ch1"][1] == "Age: embryos younger than nuclear cycle 9, i.e. before pole cells budding"
            assert record.entity_attributes["Sample_scan_protocol"] == "GeneChips were scanned using the Hewlett-Packard GeneArray Scanner G2500A."
            assert record.entity_attributes["Sample_hyb_protocol"] == "Following fragmentation, 10 microg of cRNA were hybridized for 16 hr at 45C on GeneChip Drosophila Genome Array. GeneChips were washed and stained in the Affymetrix Fluidics Station 400."
            assert record.entity_attributes["Sample_extract_protocol_ch1"] == "Trizol extraction of total RNA was performed according to the manufacturer's instructions."
            assert record.entity_attributes["Sample_source_name_ch1"] == "Drosophila embryos before nuclear cycle 9 (maternal transcripts)"
            assert record.entity_attributes["Sample_table"] == "Drosophila_T0-2.CHP"
            assert record.entity_attributes["Sample_molecule_ch1"] == "total RNA"
            assert record.entity_attributes["Sample_label_protocol_ch1"] == "Biotinylated cRNA were prepared according to the standard Affymetrix protocol from 6 microg total RNA (Expression Analysis Technical Manual, 2001, Affymetrix)."
            assert record.entity_attributes["Sample_data_processing"] == "The data were analyzed with Microarray Suite version 5.0 (MAS 5.0) using Affymetrix default analysis settings and global scaling as normalization method. The trimmed mean target intensity of each array was arbitrarily set to 100."
            assert record.entity_attributes["Sample_treatment_protocol_ch1"] == "Embryos were dechorionated with 50% bleach, put on a cover slip and covered with Halocarbon oil 27 (Sigma). Embryos of the appropriate stage were manually selected under the dissecting scope. Selected embryos were transferred to a basket, rinsed with PBS with 0,7% NaCl, 0,04% triton-X100 and placed on ice in the Trizol solution (GibcoBRL)."
            assert record.entity_attributes["Sample_title"] == "embryo at T0, biological rep2"
            assert record.entity_attributes["Sample_supplementary_file"] == "Drosophila_T0-2.CEL"
            assert len(record.col_defs) == 0
            assert len(record.table_rows) == 0
            record = next(records)
            assert record.entity_type == "SAMPLE"
            assert record.entity_id == "Drosophila_T1-1"
            assert len(record.entity_attributes) == 16
            assert record.entity_attributes["Sample_organism_ch1"] == "Drosophila melanogaster"
            assert record.entity_attributes["Sample_label_ch1"] == "biotin"
            assert record.entity_attributes["Sample_description"] == "Gene expression data from embryos in slow phase of cellularisation."
            assert record.entity_attributes["Sample_growth_protocol_ch1"] == "30 min egg collections of OreR and yw flies at 25C were aged at room temperature (RT) according to the different temporal classes T0-T4."
            assert len(record.entity_attributes["Sample_characteristics_ch1"]) == 2
            assert record.entity_attributes["Sample_characteristics_ch1"][0] == "Genotype: yellow white and Oregon R parents"
            assert record.entity_attributes["Sample_characteristics_ch1"][1] == "Age: embryos in slow phase of cellularisation"
            assert record.entity_attributes["Sample_scan_protocol"] == "GeneChips were scanned using the Hewlett-Packard GeneArray Scanner G2500A."
            assert record.entity_attributes["Sample_hyb_protocol"] == "Following fragmentation, 10 microg of cRNA were hybridized for 16 hr at 45C on GeneChip Drosophila Genome Array. GeneChips were washed and stained in the Affymetrix Fluidics Station 400."
            assert record.entity_attributes["Sample_extract_protocol_ch1"] == "Trizol extraction of total RNA was performed according to the manufacturer's instructions."
            assert record.entity_attributes["Sample_source_name_ch1"] == "Drosophila embryos in slow phase of cellularisation"
            assert record.entity_attributes["Sample_table"] == "Drosophila_T1-1.CHP"
            assert record.entity_attributes["Sample_molecule_ch1"] == "total RNA"
            assert record.entity_attributes["Sample_label_protocol_ch1"] == "Biotinylated cRNA were prepared according to the standard Affymetrix protocol from 6 microg total RNA (Expression Analysis Technical Manual, 2001, Affymetrix)."
            assert record.entity_attributes["Sample_data_processing"] == "The data were analyzed with Microarray Suite version 5.0 (MAS 5.0) using Affymetrix default analysis settings and global scaling as normalization method. The trimmed mean target intensity of each array was arbitrarily set to 100."
            assert record.entity_attributes["Sample_treatment_protocol_ch1"] == "Embryos were dechorionated with 50% bleach, put on a cover slip and covered with Halocarbon oil 27 (Sigma). Embryos of the appropriate stage were manually selected under the dissecting scope. Selected embryos were transferred to a basket, rinsed with PBS with 0,7% NaCl, 0,04% triton-X100 and placed on ice in the Trizol solution (GibcoBRL)."
            assert record.entity_attributes["Sample_title"] == "embryo at T1, biological rep1"
            assert record.entity_attributes["Sample_supplementary_file"] == "Drosophila_T1-1.CEL"
            assert len(record.col_defs) == 0
            assert len(record.table_rows) == 0
            record = next(records)
            assert record.entity_type == "SERIES"
            assert record.entity_id == "Dros_embryo_timecourse"
            assert len(record.entity_attributes) == 6
            assert len(record.entity_attributes["Series_sample_id"]) == 3
            assert record.entity_attributes["Series_sample_id"][0] == "Drosophila_T0-1"
            assert record.entity_attributes["Series_sample_id"][1] == "Drosophila_T0-2"
            assert record.entity_attributes["Series_sample_id"][2] == "Drosophila_T1-1"
            assert len(record.entity_attributes["Series_contributor"]) == 5
            assert record.entity_attributes["Series_contributor"][0] == "Jane,Doe"
            assert record.entity_attributes["Series_contributor"][1] == "John,A,Smith"
            assert record.entity_attributes["Series_contributor"][2] == "Hans,van Elton"
            assert record.entity_attributes["Series_contributor"][3] == "John,Smithers Jr"
            assert record.entity_attributes["Series_contributor"][4] == "Jie,D,Chen"
            assert len(record.entity_attributes["Series_summary"]) == 2
            assert record.entity_attributes["Series_summary"][0] == "Morphogenesis of epithelial tissues relies on the precise developmental control of cell polarity and architecture. In the early Drosophila embryo, the primary epithelium forms during cellularisation, following a tightly controlled genetic programme where specific sets of genes are up-regulated. Some of them, for instance, control membrane invagination between the nuclei anchored at the apical surface of the syncytium."
            assert record.entity_attributes["Series_summary"][1] == "We used microarrays to detail the global programme of gene expression underlying cellularisation and identified distinct classes of up-regulated genes during this process."
            assert record.entity_attributes["Series_type"] == "time course"
            assert record.entity_attributes["Series_title"] == "Expression data from early Drosophila embryo"
            assert record.entity_attributes["Series_overall_design"] == "Drosophila embryos were selected at successive stages of early development for RNA extraction and hybridization on Affymetrix microarrays. We sought to obtain homogeneous populations of embryos at each developmental stage in order to increase the temporal resolution of expression profiles. To that end, we hand-selected embryos according to morphological criteria at five time-points: before pole cell formation, i.e. before zygotic transcription (T0), during the slow phase (T1) and the fast phase (T2) of cellularisation and at the beginning (T3) and the end (T4) of gastrulation."
            assert len(record.col_defs) == 0
            assert len(record.table_rows) == 0

    def test_record_str(self):
        path = "Geo/GSM804.txt"
        with open(path, encoding="latin") as handle:
            records = Geo.parse(handle)
            record = next(records)
            assert str(record) == """\
GEO Type: SAMPLE
GEO Id: GSM804
Sample_author: Antoine,M,Snijders

Sample_author: Norma,,Nowak

Sample_author: Richard,,Segraves

Sample_author: Stephanie,,Blackwood

Sample_author: Nils,,Brown

Sample_author: Jeffery,,Conroy

Sample_author: Greg,,Hamilton

Sample_author: Anna,K,Hindle

Sample_author: Bing,,Huey

Sample_author: Karen,,Kimura

Sample_author: Sindy,,Law

Sample_author: Ken,,Myambo

Sample_author: Joel,,Palmer

Sample_author: Bauke,,Ylstra

Sample_author: Jingzhu,P,Yue

Sample_author: Joe,W,Gray

Sample_author: Ajay,N,Jain

Sample_author: Daniel,,Pinkel

Sample_author: Donna,G,Albertson

Sample_description: Coriell Cell Repositories cell line <a h
ref="http://locus.umdnj.edu/nigms/nigms_cgi/display.cgi?GM05296">GM05296</a>.

Sample_description: Fibroblast cell line derived from a 1 mo
nth old female with multiple congenital malformations, dysmorphic features, intr
auterine growth retardation, heart murmur, cleft palate, equinovarus deformity, \
\nmicrocephaly, coloboma of right iris, clinodactyly, reduced RBC catalase activit
y, and 1 copy of catalase gene.

Sample_description: Chromosome abnormalities are present.

Sample_description: Karyotype is 46,XX,-11,+der(11)inv ins(1
1;10)(11pter> 11p13::10q21>10q24::11p13>11qter)mat

Sample_organism: Homo sapiens

Sample_platform_id: GPL28

Sample_pubmed_id: 11687795

Sample_series_id: GSE16

Sample_status: Public on Feb 12 2002

Sample_submission_date: Jan 17 2002

Sample_submitter_city: San Francisco,CA,94143,USA

Sample_submitter_department: Comprehensive Cancer Center

Sample_submitter_email: albertson@cc.ucsf.edu

Sample_submitter_institute: University of California San Francisco

Sample_submitter_name: Donna,G,Albertson

Sample_submitter_phone: 415 502-8463

Sample_target_source1: Cell line GM05296

Sample_target_source2: normal male reference genomic DNA

Sample_title: CGH_Albertson_GM05296-001218

Sample_type: dual channel genomic

Column Header Definitions
    ID_REF: Unique row identifier, genome position o
    rder

    LINEAR_RATIO: Mean of replicate Cy3/Cy5 ratios

    LOG2STDDEV: Standard deviation of VALUE

    NO_REPLICATES: Number of replicate spot measurements

    VALUE: aka LOG2RATIO, mean of log base 2 of LIN
    EAR_RATIO

0: ID_REF	VALUE	LINEAR_RATIO	LOG2STDDEV	NO_REPLICATES\t
1: 1		1.047765	0.011853	3\t
2: 2				0\t
3: 3	0.008824	1.006135	0.00143	3\t
4: 4	-0.000894	0.99938	0.001454	3\t
5: 5	0.075875	1.054	0.003077	3\t
6: 6	0.017303	1.012066	0.005876	2\t
7: 7	-0.006766	0.995321	0.013881	3\t
8: 8	0.020755	1.014491	0.005506	3\t
9: 9	-0.094938	0.936313	0.012662	3\t
10: 10	-0.054527	0.96291	0.01073	3\t
11: 11	-0.025057	0.982782	0.003855	3\t
12: 12				0\t
13: 13	0.108454	1.078072	0.005196	3\t
14: 14	0.078633	1.056017	0.009165	3\t
15: 15	0.098571	1.070712	0.007834	3\t
16: 16	0.044048	1.031003	0.013651	3\t
17: 17	0.018039	1.012582	0.005471	3\t
18: 18	-0.088807	0.9403	0.010571	3\t
19: 19	0.016349	1.011397	0.007113	3\t
20: 20	0.030977	1.021704	0.016798	3\t
"""


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
