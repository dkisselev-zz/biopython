# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for Affy module."""

import os
import struct
import unittest
import pytest

pytest.importorskip("numpy")
import numpy.testing
from numpy import array
from Bio.Affy import CelFile


class AffyTest(unittest.TestCase):
    def setUp(self):
        self.affy3 = "Affy/affy_v3_example.CEL"
        self.affy4 = "Affy/affy_v4_example.CEL"
        self.affy4Bad = "Affy/affy_v4_bad_example.CEL"
        with open(self.affy4Bad, "wb") as f:
            self.writeExampleV4(f, bad=True)

    def tearDown(self):
        os.remove(self.affy4Bad)

    # tests the Affymetrix v3 parser
    def testAffy3(self):
        with open(self.affy3) as f:
            record = CelFile.read(f)
            assert len(record.DatHeader) > 0
            assert record.intensities.shape == (5, 5)
            assert record.intensities.shape == record.stdevs.shape
            assert record.intensities.shape == record.npix.shape
            assert record.ncols == 5
            assert record.nrows == 5
            assert record.version == 3
            assert record.GridCornerUL == (206, 129)
            assert record.GridCornerUR == (3570, 107)
            assert record.GridCornerLR == (3597, 3470)
            assert record.GridCornerLL == (234, 3492)
            assert record.DatHeader["filename"] == "1g_A9AF"
            assert record.DatHeader["CLS"] == 3684
            assert record.DatHeader["RWS"] == 3684
            assert record.DatHeader["XIN"] == 1
            assert record.DatHeader["YIN"] == 1
            assert record.DatHeader["VE"] == 30
            assert record.DatHeader["laser-power"] == pytest.approx(2.0, abs=5e-8)
            assert record.DatHeader["scan-date"] == "08/23/07"
            assert record.DatHeader["scan-time"] == "11:23:24"
            assert record.DatHeader["scanner-id"] == "50205880"
            assert record.DatHeader["scanner-type"] == "M10"
            assert record.DatHeader["array-type"] == "Tgondii_SNP1.1sq"
            assert record.DatHeader["filter-wavelength"] == 570
            assert record.DatHeader["arc-radius"] == pytest.approx(25356.509766, abs=5e-8)
            assert record.DatHeader["laser-spotsize"] == pytest.approx(3.5, abs=5e-8)
            assert record.DatHeader["pixel-size"] == pytest.approx(1.56, abs=5e-8)
            assert record.DatHeader["image-orientation"] == 6
            assert record.Algorithm == "Percentile"
            assert len(record.AlgorithmParameters) == 16
            assert record.AlgorithmParameters["Percentile"] == 75
            assert record.AlgorithmParameters["CellMargin"] == 2
            assert record.AlgorithmParameters["OutlierHigh"] == pytest.approx(1.500, abs=5e-8)
            assert record.AlgorithmParameters["OutlierLow"] == pytest.approx(1.004, abs=5e-8)
            assert record.AlgorithmParameters["AlgVersion"] == "6.0"
            assert record.AlgorithmParameters["FixedCellSize"] == True  # noqa: A502
            assert record.AlgorithmParameters["FullFeatureWidth"] == 7
            assert record.AlgorithmParameters["FullFeatureHeight"] == 7
            assert record.AlgorithmParameters["IgnoreOutliersInShiftRows"] == False  # noqa: A502
            assert record.AlgorithmParameters["FeatureExtraction"] == True  # noqa: A502
            assert record.AlgorithmParameters["PoolWidthExtenstion"] == 2
            assert record.AlgorithmParameters["PoolHeightExtension"] == 2
            assert record.AlgorithmParameters["UseSubgrids"] == False  # noqa: A502
            assert record.AlgorithmParameters["RandomizePixels"] == False  # noqa: A502
            assert record.AlgorithmParameters["ErrorBasis"] == "StdvMean"
            assert record.AlgorithmParameters["StdMult"] == pytest.approx(1.0, abs=5e-8)
            assert record.NumberCells == 25

            global message
            try:
                numpy.testing.assert_allclose(
                    record.intensities,
                    [
                        [234.0, 170.0, 22177.0, 164.0, 22104.0],
                        [188.0, 188.0, 21871.0, 168.0, 21883.0],
                        [188.0, 193.0, 21455.0, 198.0, 21300.0],
                        [188.0, 182.0, 21438.0, 188.0, 20945.0],
                        [193.0, 20370.0, 174.0, 20605.0, 168.0],
                    ],
                )
                message = None
            except AssertionError as err:
                message = str(err)
            if message is not None:
                raise AssertionError(message)
            try:
                numpy.testing.assert_allclose(
                    record.stdevs,
                    [
                        [24.0, 34.5, 2669.0, 19.7, 3661.2],
                        [29.8, 29.8, 2795.9, 67.9, 2792.4],
                        [29.8, 88.7, 2976.5, 62.0, 2914.5],
                        [29.8, 76.2, 2759.5, 49.2, 2762.0],
                        [38.8, 2611.8, 26.6, 2810.7, 24.1],
                    ],
                )
                message = None
            except AssertionError as err:
                message = str(err)
            if message is not None:
                raise AssertionError(message)
            try:
                numpy.testing.assert_array_equal(
                    record.npix,
                    [
                        [25, 25, 25, 25, 25],
                        [25, 25, 25, 25, 25],
                        [25, 25, 25, 25, 25],
                        [25, 25, 25, 25, 25],
                        [25, 25, 25, 25, 25],
                    ],
                )
                message = None
            except AssertionError as err:
                message = str(err)
            if message is not None:
                raise AssertionError(message)
            assert record.nmask == 3
            try:
                numpy.testing.assert_array_equal(
                    record.mask,
                    [
                        [False, False, False, False, False],
                        [False, False, False, True, True],
                        [False, False, False, False, True],
                        [False, False, False, False, False],
                        [False, False, False, False, False],
                    ],
                )
                message = None
            except AssertionError as err:
                message = str(err)
            if message is not None:
                raise AssertionError(message)
            assert record.noutliers == 3
            try:
                numpy.testing.assert_array_equal(
                    record.outliers,
                    [
                        [False, False, False, False, False],
                        [False, True, True, False, False],
                        [False, False, False, False, False],
                        [False, True, False, False, False],
                        [False, False, False, False, False],
                    ],
                )
                message = None
            except AssertionError as err:
                message = str(err)
            if message is not None:
                raise AssertionError(message)
            assert record.nmodified == 3
            try:
                numpy.testing.assert_allclose(
                    record.modified,
                    [
                        [0.0, 0.0, 0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0, 189.0, 220.0],
                        [0.0, 0.0, 0.0, 21775.0, 0.0],
                        [0.0, 0.0, 0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0, 0.0, 0.0],
                    ],
                )
                message = None
            except AssertionError as err:
                message = str(err)
            if message is not None:
                raise AssertionError(message)

    def testAffy4(self):
        with open(self.affy4, "rb") as f:
            record = CelFile.read(f)
        assert record.intensities.shape == (5, 5)
        assert record.intensities.shape == record.stdevs.shape
        assert record.intensities.shape == record.npix.shape
        assert record.ncols == 5
        assert record.nrows == 5
        global message
        try:
            numpy.testing.assert_allclose(
                record.intensities,
                [
                    [0.0, 1.0, 2.0, 3.0, 4.0],
                    [5.0, 6.0, 7.0, 8.0, 9.0],
                    [10.0, 11.0, 12.0, 13.0, 14.0],
                    [15.0, 16.0, 17.0, 18.0, 19.0],
                    [20.0, 21.0, 22.0, 23.0, 24.0],
                ],
            )
            message = None
        except AssertionError as err:
            message = str(err)
        if message is not None:
            raise AssertionError(message)
        try:
            numpy.testing.assert_allclose(
                record.stdevs,
                [
                    [0.0, -1.0, -2.0, -3.0, -4.0],
                    [-5.0, -6.0, -7.0, -8.0, -9.0],
                    [-10.0, -11.0, -12.0, -13.0, -14.0],
                    [-15.0, -16.0, -17.0, -18.0, -19.0],
                    [-20.0, -21.0, -22.0, -23.0, -24.0],
                ],
            )
            message = None
        except AssertionError as err:
            message = str(err)
        if message is not None:
            raise AssertionError(message)
        try:
            numpy.testing.assert_allclose(
                record.npix,
                [
                    [9, 9, 9, 9, 9],
                    [9, 9, 9, 9, 9],
                    [9, 9, 9, 9, 9],
                    [9, 9, 9, 9, 9],
                    [9, 9, 9, 9, 9],
                ],
            )
            message = None
        except AssertionError as err:
            message = str(err)
        if message is not None:
            raise AssertionError(message)
        assert len(record.AlgorithmParameters) == 329
        assert len(record.GridCornerUL) == 7
        assert record.AlgorithmParameters[-3:] == "169"

    def testAffyBadHeader(self):
        with pytest.raises(CelFile.ParserError):
            with open(self.affy4Bad, "rb") as f:
                record = CelFile.read(f)

    def testAffyWrongModeReadV3(self):
        with pytest.raises(ValueError):
            with open(self.affy3, "rb") as f:
                record = CelFile.read(f, version=3)

    def testAffyWrongModeReadV4(self):
        with pytest.raises(ValueError):
            with open(self.affy4) as f:
                record = CelFile.read(f, version=4)

    # Writes a small example Affymetrix V4 CEL File
    def writeExampleV4(self, f, bad=False):
        preHeaders = {
            "cellNo": 25,
            "columns": 5,
            "headerLen": 752,
            "magic": 64,
            "rows": 5,
            "version": 4,
        }
        goodH = {"Axis-invertX": b"0"}
        badH = {"Axis-invertX": b"1"}

        headers = {
            "Algorithm": b"Percentile",
            "AlgorithmParameters": b"Percentile:75;CellMargin:4;Outlie"
            b"rHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize"
            b":TRUE;FullFeatureWidth:7;FullFeatureHeight:7;IgnoreOutlie"
            b"rsInShiftRows:FALSE;FeatureExtraction:TRUE;PoolWidthExten"
            b"stion:1;PoolHeightExtension:1;UseSubgrids:FALSE;Randomize"
            b"Pixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000;NumDATS"
            b"ubgrids:169",
            "AxisInvertY": b"0",
            "Cols": b"5",
            "DatHeader": b"[0..65534]  20_10N:CLS=19420RWS=19420XIN=0"
            b" YIN=0  VE=30        2.0 05/25/05 23:19:07 50102310  M10 "
            b"  \x14  \x14 HuEx-1_0-st-v2.1sq \x14  \x14  \x14  \x14  "
            b"\x14570 \x14 25540.671875 \x14 3.500000 \x14 0.7000 \x14"
            b" 3",
            "GridCornerLL": b"518 18668",
            "GridCornerLR": b"18800 18825",
            "GridCornerUL": b"659 469",
            "GridCornerUR": b"18942 623",
            "OffsetX": b"0",
            "OffsetY": b"0",
            "Rows": b"5",
            "TotalX": b"2560",
            "TotalY": b"2560",
            "swapXY": b"0",
        }
        if not bad:
            headers.update(goodH)
        else:
            headers.update(badH)
        prePadding = b"this text doesn't matter and is ignored\x04"
        preHeadersOrder = ["magic", "version", "columns", "rows", "cellNo", "headerLen"]
        headersEncoded = struct.pack(
            "<" + "i" * len(preHeadersOrder),
            *(preHeaders[header] for header in preHeadersOrder),
        )

        f.write(headersEncoded)
        for header in headers:
            try:
                f.write(
                    bytes(header, encoding="utf-8") + b"=" + headers[header] + b"\n"
                )
            except TypeError:
                f.write(header + b"=" + headers[header] + b"\n")
        f.write(prePadding)
        f.write(b"\x00" * 15)
        for i in range(25):
            f.write(struct.pack("< f f h", i, -i, 9))  # intensity, sdev, pixel


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
