# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for Cluster module."""

import unittest
import pytest

np = pytest.importorskip("numpy")
class TestCluster(unittest.TestCase):
    module = "Bio.Cluster"

    def test_matrix_parse(self):
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster import treecluster
        elif TestCluster.module == "Pycluster":
            from Pycluster import treecluster

        # Normal matrix, no errors
        data1 = np.array(
            [
                [1.1, 1.2],
                [1.4, 1.3],
                [1.1, 1.5],
                [2.0, 1.5],
                [1.7, 1.9],
                [1.7, 1.9],
                [5.7, 5.9],
                [5.7, 5.9],
                [3.1, 3.3],
                [5.4, 5.3],
                [5.1, 5.5],
                [5.0, 5.5],
                [5.1, 5.2],
            ]
        )

        # Another normal matrix, no errors; written as a list
        data2 = [
            [1.1, 2.2, 3.3, 4.4, 5.5],
            [3.1, 3.2, 1.3, 2.4, 1.5],
            [4.1, 2.2, 0.3, 5.4, 0.5],
            [2.1, 2.0, 0.0, 5.0, 0.0],
        ]

        # Rows are not contiguous
        data3 = data1[::2, :]

        # Columns are not contiguous
        data4 = np.array(data2)[:, ::2]

        # Matrix using float32
        data5 = np.array(
            [
                [1.1, 2.2, 3.3, 4.4, 5.5],
                [3.1, 3.2, 1.3, 2.4, 1.5],
                [4.1, 2.2, 0.3, 5.4, 0.5],
                [2.1, 2.0, 0.0, 5.0, 0.0],
            ],
            np.float32,
        )

        # Matrix using int
        # fmt: off
        data6 = np.array(
            [
                [1, 2, 3, 4, 5],
                [3, 3, 1, 2, 1],
                [4, 2, 0, 5, 0],
                [2, 2, 0, 5, 0]
            ],
            np.int32,
        )
        # fmt: on
        try:
            treecluster(data1)
        except Exception:
            raise AssertionError("treecluster failed to accept matrix data1")

        try:
            treecluster(data2)
        except Exception:
            raise AssertionError("treecluster failed to accept matrix data2")

        try:
            treecluster(data3)
        except Exception:
            raise AssertionError("treecluster failed to accept matrix data3")

        try:
            treecluster(data4)
        except Exception:
            raise AssertionError("treecluster failed to accept matrix data4")

        try:
            treecluster(data5)
        except Exception:
            raise AssertionError("treecluster failed to accept matrix data5")

        try:
            treecluster(data6)
        except Exception:
            raise AssertionError("treecluster failed to accept matrix data6")

        # Ragged matrix
        data7 = [
            [91.1, 92.2, 93.3, 94.4, 95.5],
            [93.1, 93.2, 91.3, 92.4],
            [94.1, 92.2, 90.3],
            [12.1, 92.0, 90.0, 95.0, 90.0],
        ]

        # Matrix with bad cells
        data8 = [
            [7.1, 7.2, 7.3, 7.4, 7.5],
            [7.1, 7.2, 7.3, 7.4, "snoopy"],
            [7.1, 7.2, 7.3, None, None],
        ]

        # Matrix with a bad row
        # fmt: off
        data9 = [
            [23.1, 23.2, 23.3, 23.4, 23.5],
            None,
            [23.1, 23.0, 23.0, 23.0, 23.0]
        ]
        # fmt: on

        # Various references that don't point to matrices at all
        data10 = "snoopy"
        data11 = {"a": [[2.3, 1.2], [3.3, 5.6]]}
        data12 = []
        data13 = [None]

        # Array of incorrect rank
        data14 = np.array(
            [
                [[1.1, 1.2], [2.3, 1.2], [3.4, 1.6]],
                [[1.4, 1.3], [3.2, 4.5], [9.8, 4.9]],
                [[1.1, 1.5], [1.1, 2.3], [6.5, 0.4]],
            ]
        )

        # Array with non-numerical values
        data15 = np.array([["a", "b", "c"], ["e", "f", "g"]], "c")

        # Empty array
        data16 = np.array([[]], "d")

        with pytest.raises(ValueError):
            treecluster(data7)
        with pytest.raises(ValueError):
            treecluster(data8)
        with pytest.raises(ValueError):
            treecluster(data9)
        with pytest.raises(ValueError):
            treecluster(data10)
        with pytest.raises(TypeError):
            treecluster(data11)
        with pytest.raises(ValueError):
            treecluster(data12)
        with pytest.raises(ValueError):
            treecluster(data13)
        with pytest.raises(ValueError):
            treecluster(data14)
        with pytest.raises(ValueError):
            treecluster(data15)
        with pytest.raises(ValueError):
            treecluster(data16)

    def test_mask_parse(self):
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster import treecluster
        elif TestCluster.module == "Pycluster":
            from Pycluster import treecluster

        # data matrix
        data = np.array(
            [
                [1.1, 2.2, 3.3, 4.4, 5.5],
                [3.1, 3.2, 1.3, 2.4, 1.5],
                [4.1, 2.2, 0.3, 5.4, 0.5],
                [2.1, 2.0, 0.0, 5.0, 0.0],
            ]
        )

        # Normal mask, no errors
        # fmt: off
        mask1 = np.array(
            [
                [1, 1, 0, 1, 0],
                [1, 1, 1, 0, 0],
                [1, 1, 0, 1, 1],
                [1, 0, 1, 1, 0],
            ]
        )

        # Same mask, no errors; written as a list
        mask2 = [
            [1, 1, 0, 1, 0],
            [1, 1, 1, 0, 0],
            [1, 1, 0, 1, 1],
            [1, 0, 1, 1, 0],
        ]
        # fmt: on

        # Rows are not contiguous
        mask3 = np.array(
            [
                [1, 1, 0, 1, 0],
                [1, 1, 1, 0, 0],
                [1, 1, 1, 0, 0],
                [1, 1, 0, 1, 1],
                [1, 1, 1, 0, 0],
                [1, 1, 0, 1, 1],
                [1, 1, 0, 1, 1],
                [1, 0, 1, 1, 0],
            ]
        )
        mask3 = mask3[::2, :]

        # Columns are not contiguous
        mask4 = np.array(
            [
                [1, 1, 0, 1, 0, 1, 0, 0, 1, 1],
                [1, 1, 1, 0, 0, 1, 1, 0, 0, 1],
                [1, 1, 0, 1, 1, 1, 0, 1, 1, 0],
                [1, 0, 1, 1, 0, 1, 0, 0, 1, 1],
            ]
        )
        mask4 = mask4[:, ::2]

        # Matrix using int16
        # fmt: off
        mask5 = np.array(
            [
                [1, 1, 0, 1, 0],
                [1, 1, 1, 0, 0],
                [1, 1, 1, 0, 0],
                [1, 1, 0, 1, 1],
            ],
            np.int16,
        )
        # fmt: on

        # Matrix using float
        mask6 = np.array(
            [
                [1.0, 2.2, 3.1, 4.8, 5.1],
                [3.3, 3.3, 1.4, 2.4, 1.2],
                [4.1, 2.2, 0.6, 5.5, 0.6],
                [2.7, 2.5, 0.4, 5.7, 0.2],
            ],
            float,
        )
        try:
            treecluster(data, mask1)
        except Exception:
            raise AssertionError("treecluster failed to accept matrix mask1")

        try:
            treecluster(data, mask2)
        except Exception:
            raise AssertionError("treecluster failed to accept matrix mask2")

        try:
            treecluster(data, mask3)
        except Exception:
            raise AssertionError("treecluster failed to accept matrix mask3")

        try:
            treecluster(data, mask4)
        except Exception:
            raise AssertionError("treecluster failed to accept matrix mask4")

        try:
            treecluster(data, mask5)
        except Exception:
            raise AssertionError("treecluster failed to accept matrix mask5")

        try:
            treecluster(data, mask6)
        except Exception:
            raise AssertionError("treecluster failed to accept matrix mask6")

        # Ragged mask
        # fmt: off
        mask7 = [
            [1, 1, 0, 1],
            [1, 1, 1, 0, 0],
            [1, 1, 0, 1, 1],
            [1, 1, 0],
        ]
        # fmt: on

        # Mask with incorrect number of rows
        mask8 = np.array(
            [
                [1, 1, 0, 1, 0],
                [1, 1, 1, 0, 0],
                [1, 1, 0, 1, 1],
                [0, 1, 1, 0, 1],
                [1, 0, 1, 1, 0],
            ]
        )

        # Mask with incorrect number of columns
        mask9 = np.array(
            [
                [1, 1, 0, 1, 0, 1],
                [1, 1, 1, 0, 0, 0],
                [0, 1, 1, 0, 1, 1],
                [1, 0, 1, 1, 0, 1],
            ]
        )

        # Matrix with bad cells
        mask10 = [
            [1, 1, 0, 1, 0],
            [1, 1, 1, 0, "snoopy"],
            [1, 1, 0, 1, 1],
            [1, 0, 1, 1, 0],
        ]

        # Matrix with a bad row
        # fmt: off
        mask11 = [
            [1, 1, 0, 1, 0],
            None,
            [1, 1, 0, 1, 1],
            [1, 0, 1, 1, 0],
        ]
        # fmt: on

        # Array with non-numerical values
        mask12 = np.array([["a", "b", "c"], ["e", "f", "g"]], "c")

        # Empty arrays
        mask13 = np.array([[]], "d")
        mask14 = []

        # Array of incorrect rank
        mask15 = np.array(
            [
                [[1, 1], [0, 1], [1, 1]],
                [[1, 1], [0, 1], [1, 1]],
                [[1, 1], [1, 1], [1, 0]],
            ]
        )

        # References that cannot be converted to a matrix of int
        mask16 = "snoopy"
        mask17 = {"a": [[1, 0], [1, 1]]}
        mask18 = [None]

        with pytest.raises(ValueError):
            treecluster(data, mask7)
        with pytest.raises(ValueError):
            treecluster(data, mask8)
        with pytest.raises(ValueError):
            treecluster(data, mask9)
        with pytest.raises(ValueError):
            treecluster(data, mask10)
        with pytest.raises(ValueError):
            treecluster(data, mask11)
        with pytest.raises(ValueError):
            treecluster(data, mask12)
        with pytest.raises(ValueError):
            treecluster(data, mask13)
        with pytest.raises(ValueError):
            treecluster(data, mask14)
        with pytest.raises(ValueError):
            treecluster(data, mask15)
        with pytest.raises(ValueError):
            treecluster(data, mask16)
        with pytest.raises(TypeError):
            treecluster(data, mask17)
        with pytest.raises(TypeError):
            treecluster(data, mask18)

    def test_kcluster_arguments(self):
        # Test if incorrect arguments are caught by the C code
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster._cluster import clustercentroids
            from Bio.Cluster._cluster import kcluster
        elif TestCluster.module == "Pycluster":
            from Pycluster._cluster import clustercentroids
            from Pycluster._cluster import kcluster

        nclusters = 3
        weight = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
        data = np.array(
            [
                [1.1, 2.2, 3.3, 4.4, 5.5],
                [3.1, 3.2, 1.3, 2.4, 1.5],
                [4.1, 2.2, 0.3, 5.4, 0.5],
                [9.9, 2.0, 0.0, 5.0, 0.0],
            ]
        )
        # fmt: off
        mask = np.array(
            [
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
            ],
            np.int32,
        )
        # fmt: on
        nrows, ncols = data.shape
        clusterid = np.zeros(nrows, np.int32)

        message = "^data matrix is empty$"
        with pytest.raises(ValueError, match=message):
            kcluster(
                data[:0, :],
                nclusters=nclusters,
                mask=mask,
                weight=weight,
                transpose=False,
                npass=100,
                method="a",
                dist="e",
            )
        message = "^mask has incorrect rank 1 \\(expected 2\\)$"
        with pytest.raises(ValueError, match=message):
            kcluster(
                data,
                nclusters=nclusters,
                mask=np.zeros(3),
                weight=weight,
                transpose=False,
                npass=100,
                method="a",
                dist="e",
            )
        message = "^mask has incorrect dimensions 4 x 3 \\(expected 4 x 5\\)$"
        with pytest.raises(ValueError, match=message):
            kcluster(
                data,
                nclusters=nclusters,
                mask=np.zeros((4, 3), np.int32),
                weight=weight,
                transpose=False,
                npass=100,
                method="a",
                dist="e",
                clusterid=clusterid,
            )
        message = "^incorrect rank 2 \\(expected 1\\)$"
        with pytest.raises(ValueError, match=message):
            kcluster(
                data,
                nclusters=nclusters,
                mask=mask,
                weight=np.zeros((2, 2)),
                transpose=False,
                npass=100,
                method="a",
                dist="e",
            )
        message = "^weight has incorrect size 3 \\(expected 5\\)$"
        with pytest.raises(ValueError, match=message):
            kcluster(
                data,
                nclusters=nclusters,
                mask=mask,
                weight=np.zeros(3),
                transpose=False,
                npass=100,
                method="a",
                dist="e",
                clusterid=clusterid,
            )
        message = "^nclusters should be positive$"
        with pytest.raises(ValueError, match=message):
            kcluster(
                data,
                nclusters=-1,
                mask=mask,
                weight=weight,
                transpose=False,
                npass=100,
                method="a",
                dist="e",
                clusterid=clusterid,
            )
        message = "^more clusters than items to be clustered$"
        with pytest.raises(ValueError, match=message):
            kcluster(
                data,
                nclusters=1234,
                mask=mask,
                weight=weight,
                transpose=False,
                npass=100,
                method="a",
                dist="e",
                clusterid=clusterid,
            )
        message = "^incorrect size \\(3, expected 4\\)$"
        with pytest.raises(ValueError, match=message):
            kcluster(
                data,
                nclusters=nclusters,
                mask=mask,
                weight=weight,
                transpose=False,
                npass=0,
                method="a",
                dist="e",
                clusterid=clusterid[:3],
            )
        message = "^more clusters requested than found in clusterid$"
        with pytest.raises(ValueError, match=message):
            kcluster(
                data,
                nclusters=nclusters,
                mask=mask,
                weight=weight,
                transpose=False,
                npass=0,
                method="a",
                dist="e",
                clusterid=clusterid,
            )
        clusterid = np.array([0, -1, 2, 3], np.int32)
        message = "^negative cluster number found$"
        with pytest.raises(ValueError, match=message):
            kcluster(
                data,
                nclusters=nclusters,
                mask=mask,
                weight=weight,
                transpose=False,
                npass=0,
                method="a",
                dist="e",
                clusterid=clusterid,
            )
        clusterid = np.array([0, 0, 2, 3], np.int32)
        message = "^cluster 1 is empty$"
        with pytest.raises(ValueError, match=message):
            kcluster(
                data,
                nclusters=nclusters,
                mask=mask,
                weight=weight,
                transpose=False,
                npass=0,
                method="a",
                dist="e",
                clusterid=clusterid,
            )

    def test_kcluster(self):
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster import clustercentroids
            from Bio.Cluster import kcluster
        elif TestCluster.module == "Pycluster":
            from Pycluster import clustercentroids
            from Pycluster import kcluster

        nclusters = 3

        # First data set
        weight = np.array([1, 1, 1, 1, 1])
        data = np.array(
            [
                [1.1, 2.2, 3.3, 4.4, 5.5],
                [3.1, 3.2, 1.3, 2.4, 1.5],
                [4.1, 2.2, 0.3, 5.4, 0.5],
                [9.9, 2.0, 0.0, 5.0, 0.0],
            ]
        )
        # fmt: off
        mask = np.array(
            [
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
            ],
            int,
        )
        # fmt: on
        nrows, ncols = data.shape

        clusterid, error, nfound = kcluster(
            data,
            nclusters=nclusters,
            mask=mask,
            weight=weight,
            transpose=False,
            npass=100,
            method="a",
            dist="e",
        )
        assert len(clusterid) == len(data)

        correct = [0, 1, 1, 2]
        mapping = [clusterid[correct.index(i)] for i in range(nclusters)]
        for i in range(len(clusterid)):
            assert clusterid[i] == mapping[correct[i]]

        cdata, cmask = clustercentroids(
            data, mask=mask, clusterid=clusterid, method="a", transpose=False
        )

        assert cdata.shape == (nclusters, ncols)
        assert cmask.shape == (nclusters, ncols)
        for value in cmask.flat:
            assert value == 1

        correct = np.array(
            [
                [1.1, 2.2, 3.3, 4.4, 5.5],
                [3.6, 2.7, 0.8, 3.9, 1.0],
                [9.9, 2.0, 0.0, 5.0, 0.0],
            ]
        )
        for i in range(nclusters):
            for j in range(ncols):
                assert cdata[mapping[i], j] == pytest.approx(correct[i, j], abs=5e-8)

        # First data set, using transpose=True
        weight = np.array([1, 1, 1, 1])
        clusterid, error, nfound = kcluster(
            data,
            nclusters=nclusters,
            mask=mask,
            weight=weight,
            transpose=True,
            npass=100,
            method="a",
            dist="e",
        )
        assert len(clusterid) == ncols

        correct = [0, 1, 1, 2, 1]
        mapping = [clusterid[correct.index(i)] for i in range(nclusters)]
        for i in range(len(clusterid)):
            assert clusterid[i] == mapping[correct[i]]

        cdata, cmask = clustercentroids(
            data, mask=mask, clusterid=clusterid, method="a", transpose=True
        )

        assert cdata.shape == (nrows, nclusters)
        assert cmask.shape == (nrows, nclusters)
        for value in cmask.flat:
            assert value == 1

        correct = np.array(
            [
                [1.1, 3.6666666667, 4.4],
                [3.1, 2.0000000000, 2.4],
                [4.1, 1.0000000000, 5.4],
                [9.9, 0.6666666667, 5.0],
            ]
        )
        for i in range(nrows):
            for j in range(nclusters):
                assert cdata[i, mapping[j]] == pytest.approx(correct[i, j], abs=5e-8)

        # Second data set
        weight = np.array([1, 1])
        data = np.array(
            [
                [1.1, 1.2],
                [1.4, 1.3],
                [1.1, 1.5],
                [2.0, 1.5],
                [1.7, 1.9],
                [1.7, 1.9],
                [5.7, 5.9],
                [5.7, 5.9],
                [3.1, 3.3],
                [5.4, 5.3],
                [5.1, 5.5],
                [5.0, 5.5],
                [5.1, 5.2],
            ]
        )
        mask = np.array(
            [
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
            ],
            int,
        )
        nrows, ncols = data.shape

        clusterid, error, nfound = kcluster(
            data,
            nclusters=3,
            mask=mask,
            weight=weight,
            transpose=False,
            npass=100,
            method="a",
            dist="e",
        )
        assert len(clusterid) == len(data)

        correct = [0, 0, 0, 0, 0, 0, 1, 1, 2, 1, 1, 1, 1]
        mapping = [clusterid[correct.index(i)] for i in range(nclusters)]
        for i in range(len(clusterid)):
            assert clusterid[i] == mapping[correct[i]]

        cdata, cmask = clustercentroids(
            data, mask=mask, clusterid=clusterid, method="a", transpose=False
        )

        assert cdata.shape == (nclusters, ncols)
        assert cmask.shape == (nclusters, ncols)
        for value in cmask.flat:
            assert value == 1

        correct = np.array([[1.5000000, 1.55], [5.3333333, 5.55], [3.1000000, 3.30]])
        for i in range(nclusters):
            for j in range(ncols):
                assert cdata[mapping[i], j] == pytest.approx(correct[i, j], abs=5e-8)

    def test_clusterdistance_arguments(self):
        # Test if incorrect arguments are caught by the C code
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster._cluster import clusterdistance
        elif TestCluster.module == "Pycluster":
            from Pycluster._cluster import clusterdistance

        # First data set
        weight = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
        data = np.array(
            [
                [1.1, 2.2, 3.3, 4.4, 5.5],
                [3.1, 3.2, 1.3, 2.4, 1.5],
                [4.1, 2.2, 0.3, 5.4, 0.5],
                [9.9, 2.0, 0.0, 5.0, 0.0],
            ]
        )
        # fmt: off
        mask = np.array(
            [
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
            ],
            np.int32,
        )
        # fmt: on

        # Cluster assignments
        c1 = np.array([0], np.int32)
        c2 = np.array([1, 2], np.int32)
        c3 = np.array([3], np.int32)

        message = "^data is None$"
        with pytest.raises(RuntimeError, match=message):
            clusterdistance(
                data=None,
                mask=mask,
                weight=weight,
                index1=c1,
                index2=c2,
                dist="e",
                method="a",
                transpose=False,
            )
        message = "^data matrix has unexpected format.$"
        with pytest.raises(RuntimeError, match=message):
            clusterdistance(
                data=[None],
                mask=mask,
                weight=weight,
                index1=c1,
                index2=c2,
                dist="e",
                method="a",
                transpose=False,
            )
        message = "^data matrix has incorrect rank 1 \\(expected 2\\)$"
        with pytest.raises(RuntimeError, match=message):
            clusterdistance(
                data=np.zeros(3),
                mask=mask,
                weight=weight,
                index1=c1,
                index2=c2,
                dist="e",
                method="a",
                transpose=False,
            )
        message = "^data matrix has incorrect data type$"
        with pytest.raises(RuntimeError, match=message):
            clusterdistance(
                data=np.zeros((3, 3), dtype=np.int16),
                mask=mask,
                weight=weight,
                index1=c1,
                index2=c2,
                dist="e",
                method="a",
                transpose=False,
            )
        message = "^data matrix is empty$"
        with pytest.raises(ValueError, match=message):
            clusterdistance(
                data=data[:0],
                mask=mask,
                weight=weight,
                index1=c1,
                index2=c2,
                dist="e",
                method="a",
                transpose=False,
            )
        message = "^data is not contiguous$"
        with pytest.raises(RuntimeError, match=message):
            clusterdistance(
                data=data[:, ::2],
                mask=mask,
                weight=weight,
                index1=c1,
                index2=c2,
                dist="e",
                method="a",
                transpose=False,
            )
        message = "^mask has unexpected format.$"
        with pytest.raises(RuntimeError, match=message):
            clusterdistance(
                data=data,
                mask=[None],
                weight=weight,
                index1=c1,
                index2=c2,
                dist="e",
                method="a",
                transpose=False,
            )
        message = "^mask has incorrect rank 1 \\(expected 2\\)$"
        with pytest.raises(ValueError, match=message):
            clusterdistance(
                data=data,
                mask=np.zeros(3),
                weight=weight,
                index1=c1,
                index2=c2,
                dist="e",
                method="a",
                transpose=False,
            )
        message = "^mask has incorrect data type$"
        with pytest.raises(RuntimeError, match=message):
            clusterdistance(
                data=data,
                mask=np.ones((2, 2), dtype=np.int16),
                weight=weight,
                index1=c1,
                index2=c2,
                dist="e",
                method="a",
                transpose=False,
            )
        message = "^mask is not contiguous$"
        with pytest.raises(RuntimeError, match=message):
            clusterdistance(
                data=data,
                mask=mask[:, ::2],
                weight=weight,
                index1=c1,
                index2=c2,
                dist="e",
                method="a",
                transpose=False,
            )
        message = "^unexpected format.$"
        with pytest.raises(RuntimeError, match=message):
            clusterdistance(
                data=data,
                mask=mask,
                weight="nothing",
                index1=c1,
                index2=c2,
                dist="e",
                method="a",
                transpose=False,
            )
        message = "^incorrect rank 2 \\(expected 1\\)$"
        with pytest.raises(ValueError, match=message):
            clusterdistance(
                data=data,
                mask=mask,
                weight=np.zeros((2, 2)),
                index1=c1,
                index2=c2,
                dist="e",
                method="a",
                transpose=False,
            )
        message = "^array has incorrect data type$"
        with pytest.raises(RuntimeError, match=message):
            clusterdistance(
                data=data,
                mask=mask,
                weight=np.ones(3, dtype=np.int16),
                index1=c1,
                index2=c2,
                dist="e",
                method="a",
                transpose=False,
            )
        message = "^unexpected format.$"
        with pytest.raises(RuntimeError, match=message):
            clusterdistance(
                data=data,
                mask=mask,
                weight=weight,
                index1=None,
                index2=c2,
                dist="e",
                method="a",
                transpose=False,
            )
        message = "^incorrect rank 2 \\(expected 1\\)$"
        with pytest.raises(ValueError, match=message):
            clusterdistance(
                data=data,
                mask=mask,
                weight=weight,
                index1=np.zeros((2, 2)),
                index2=c2,
                dist="e",
                method="a",
                transpose=False,
            )
        message = "^argument has incorrect data type$"
        with pytest.raises(RuntimeError, match=message):
            clusterdistance(
                data=data,
                mask=mask,
                weight=weight,
                index1=np.zeros(2, np.int16),
                index2=c2,
                dist="e",
                method="a",
                transpose=False,
            )
        message = "^unexpected format.$"
        with pytest.raises(RuntimeError, match=message):
            clusterdistance(
                data=data,
                mask=mask,
                weight=weight,
                index1=c1,
                index2=None,
                dist="e",
                method="a",
                transpose=False,
            )
        message = "^incorrect rank 2 \\(expected 1\\)$"
        with pytest.raises(ValueError, match=message):
            clusterdistance(
                data=data,
                mask=mask,
                weight=weight,
                index1=c1,
                index2=np.zeros((2, 2)),
                dist="e",
                method="a",
                transpose=False,
            )
        message = "^argument has incorrect data type$"
        with pytest.raises(RuntimeError, match=message):
            clusterdistance(
                data=data,
                mask=mask,
                weight=weight,
                index1=c1,
                index2=np.zeros(2, np.int16),
                dist="e",
                method="a",
                transpose=False,
            )

    def test_clusterdistance(self):
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster import clusterdistance
        elif TestCluster.module == "Pycluster":
            from Pycluster import clusterdistance

        # First data set
        weight = np.array([1, 1, 1, 1, 1])
        data = np.array(
            [
                [1.1, 2.2, 3.3, 4.4, 5.5],
                [3.1, 3.2, 1.3, 2.4, 1.5],
                [4.1, 2.2, 0.3, 5.4, 0.5],
                [9.9, 2.0, 0.0, 5.0, 0.0],
            ]
        )
        # fmt: off
        mask = np.array(
            [
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
            ],
            int,
        )
        # fmt: on

        # Cluster assignments
        c1 = [0]
        c2 = [1, 2]
        c3 = [3]

        distance = clusterdistance(
            data,
            mask=mask,
            weight=weight,
            index1=c1,
            index2=c2,
            dist="e",
            method="a",
            transpose=False,
        )
        assert distance == pytest.approx(6.650, abs=0.0005)
        distance = clusterdistance(
            data,
            mask=mask,
            weight=weight,
            index1=c1,
            index2=c3,
            dist="e",
            method="a",
            transpose=False,
        )
        assert distance == pytest.approx(23.796, abs=0.0005)
        distance = clusterdistance(
            data,
            mask=mask,
            weight=weight,
            index1=c2,
            index2=c3,
            dist="e",
            method="a",
            transpose=False,
        )
        assert distance == pytest.approx(8.606, abs=0.0005)

        # First data set, using transpose=True
        weight = np.array([1, 1, 1, 1])

        # Cluster assignments
        c1 = [0, 2]
        c2 = [1, 4]
        c3 = [3]

        distance = clusterdistance(
            data,
            mask=mask,
            weight=weight,
            index1=c1,
            index2=c2,
            dist="e",
            method="a",
            transpose=True,
        )
        assert distance == pytest.approx(4.7675, abs=0.0005)
        distance = clusterdistance(
            data,
            mask=mask,
            weight=weight,
            index1=c1,
            index2=c3,
            dist="e",
            method="a",
            transpose=True,
        )
        assert distance == pytest.approx(3.780625, abs=0.0005)
        distance = clusterdistance(
            data,
            mask=mask,
            weight=weight,
            index1=c2,
            index2=c3,
            dist="e",
            method="a",
            transpose=True,
        )
        assert distance == pytest.approx(8.176875, abs=0.0005)

        # Second data set
        weight = np.array([1, 1])
        data = np.array(
            [
                [1.1, 1.2],
                [1.4, 1.3],
                [1.1, 1.5],
                [2.0, 1.5],
                [1.7, 1.9],
                [1.7, 1.9],
                [5.7, 5.9],
                [5.7, 5.9],
                [3.1, 3.3],
                [5.4, 5.3],
                [5.1, 5.5],
                [5.0, 5.5],
                [5.1, 5.2],
            ]
        )
        mask = np.array(
            [
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
            ],
            int,
        )

        # Cluster assignments
        c1 = [0, 1, 2, 3]
        c2 = [4, 5, 6, 7]
        c3 = [8]

        distance = clusterdistance(
            data,
            mask=mask,
            weight=weight,
            index1=c1,
            index2=c2,
            dist="e",
            method="a",
            transpose=False,
        )
        assert distance == pytest.approx(5.833, abs=0.0005)
        distance = clusterdistance(
            data,
            mask=mask,
            weight=weight,
            index1=c1,
            index2=c3,
            dist="e",
            method="a",
            transpose=False,
        )
        assert distance == pytest.approx(3.298, abs=0.0005)
        distance = clusterdistance(
            data,
            mask=mask,
            weight=weight,
            index1=c2,
            index2=c3,
            dist="e",
            method="a",
            transpose=False,
        )
        assert distance == pytest.approx(0.360, abs=0.0005)

    def test_treecluster_arguments(self):
        # Test if incorrect arguments are caught by the C code
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster._cluster import Tree
            from Bio.Cluster._cluster import treecluster
        elif TestCluster.module == "Pycluster":
            from Pycluster._cluster import Tree
            from Pycluster._cluster import treecluster
        weight = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
        data = np.array(
            [
                [1.1, 2.2, 3.3, 4.4, 5.5],
                [3.1, 3.2, 1.3, 2.4, 1.5],
                [4.1, 2.2, 0.3, 5.4, 0.5],
                [9.7, 2.0, 0.0, 5.0, 0.0],
            ]
        )
        # fmt: off
        mask = np.array(
            [
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
            ],
            np.int32,
        )
        # fmt: on

        message = "^argument 1 must be _cluster.Tree, not None$"
        with pytest.raises(TypeError, match=message):
            treecluster(
                None,
                data=data,
                mask=mask,
                weight=weight,
                transpose=False,
                method="a",
                dist="e",
                distancematrix=None,
            )
        tree = Tree()
        message = "^neither data nor distancematrix was given$"
        with pytest.raises(ValueError, match=message):
            treecluster(
                tree,
                data=None,
                mask=mask,
                weight=weight,
                transpose=False,
                method="a",
                dist="e",
                distancematrix=None,
            )
        message = "^data matrix has unexpected format.$"
        with pytest.raises(RuntimeError, match=message):
            treecluster(
                tree,
                data=[],
                mask=mask,
                weight=weight,
                transpose=False,
                method="a",
                dist="e",
                distancematrix=None,
            )
        message = "^data matrix has incorrect data type$"
        with pytest.raises(RuntimeError, match=message):
            treecluster(
                tree,
                data=np.zeros((3, 3), np.int32),
                mask=mask,
                weight=weight,
                transpose=False,
                method="a",
                dist="e",
                distancematrix=None,
            )
        message = "^data matrix has incorrect rank 1 \\(expected 2\\)$"
        with pytest.raises(RuntimeError, match=message):
            treecluster(
                tree,
                data=np.zeros(3),
                mask=mask,
                weight=weight,
                transpose=False,
                method="a",
                dist="e",
                distancematrix=None,
            )

    def test_tree_arguments(self):
        # Test if incorrect arguments are caught by the C code
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster._cluster import Node
            from Bio.Cluster._cluster import Tree
        elif TestCluster.module == "Pycluster":
            from Pycluster._cluster import Node
            from Pycluster._cluster import Tree

        nodes = [Node(1, 2, 0.2), Node(0, -1, 0.5), Node(3, -2, 0.6)]
        indices = np.zeros(4, np.int32)
        tree = Tree(nodes)
        message = "^requested number of clusters should be positive$"
        with pytest.raises(ValueError, match=message):
            tree.cut(indices, -5)
        message = "^more clusters requested than items available$"
        with pytest.raises(ValueError, match=message):
            tree.cut(indices, +5)
        message = "^unexpected format.$"
        with pytest.raises(RuntimeError, match=message):
            tree.sort(indices, "nothing")
        message = "^incorrect rank 2 \\(expected 1\\)$"
        with pytest.raises(ValueError, match=message):
            tree.sort(indices, np.zeros((5, 5)))
        message = "^order array has incorrect size 2 \\(expected 4\\)$"
        with pytest.raises(ValueError, match=message):
            tree.sort(indices, np.zeros(2))
        message = "^order array has incorrect size 6 \\(expected 4\\)$"
        with pytest.raises(ValueError, match=message):
            tree.sort(indices, np.zeros(6))

    def test_tree(self):
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster import Node
            from Bio.Cluster import Tree
        elif TestCluster.module == "Pycluster":
            from Pycluster import Node
            from Pycluster import Tree

        node = Node(2, 3)
        assert node.left == 2
        assert node.right == 3
        assert node.distance == pytest.approx(0.0, abs=0.0005)
        node.left = 6
        node.right = 2
        node.distance = 0.73
        assert node.left == 6
        assert node.right == 2
        assert node.distance == pytest.approx(0.73, abs=0.0005)
        nodes = [Node(1, 2, 0.2), Node(0, 3, 0.5), Node(-2, 4, 0.6), Node(-1, -3, 0.9)]
        try:
            tree = Tree(nodes)
        except Exception:
            raise AssertionError("failed to construct tree from nodes")
        nodes = [Node(1, 2, 0.2), Node(0, 2, 0.5)]
        with pytest.raises(ValueError):
            Tree(nodes)
        nodes = [Node(1, 2, 0.2), Node(0, -1, 0.5)]
        tree = Tree(nodes)
        assert tree[0].left == 1
        assert tree[0].right == 2
        assert tree[0].distance == pytest.approx(0.2, abs=5e-8)
        assert tree[1].left == 0
        assert tree[1].right == -1
        assert tree[1].distance == pytest.approx(0.5, abs=5e-8)
        tree = Tree([Node(1, 2, 0.1), Node(0, -1, 0.5), Node(-2, 3, 0.9)])
        nodes = tree[:]
        nodes[0] = Node(0, 1, 0.2)
        nodes[1].left = 2
        tree = Tree(nodes)
        assert tree[0].left == 0
        assert tree[0].right == 1
        assert tree[0].distance == pytest.approx(0.2, abs=5e-8)
        assert tree[1].left == 2
        assert tree[1].right == -1
        assert tree[1].distance == pytest.approx(0.5, abs=5e-8)
        assert tree[2].left == -2
        assert tree[2].right == 3
        assert tree[2].distance == pytest.approx(0.9, abs=5e-8)

    def test_treecluster(self):
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster import treecluster
        elif TestCluster.module == "Pycluster":
            from Pycluster import treecluster

        # First data set
        weight1 = [1, 1, 1, 1, 1]
        data1 = np.array(
            [
                [1.1, 2.2, 3.3, 4.4, 5.5],
                [3.1, 3.2, 1.3, 2.4, 1.5],
                [4.1, 2.2, 0.3, 5.4, 0.5],
                [9.7, 2.0, 0.0, 5.0, 0.0],
            ]
        )
        # fmt: off
        mask1 = np.array(
            [
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
            ],
            int,
        )
        # fmt: on

        # Pairwise average-linkage clustering
        tree = treecluster(
            data=data1,
            mask=mask1,
            weight=weight1,
            transpose=False,
            method="a",
            dist="e",
        )
        assert len(tree) == len(data1) - 1
        assert tree[0].left == 2
        assert tree[0].right == 1
        assert tree[0].distance == pytest.approx(2.600, abs=0.0005)
        assert tree[1].left == -1
        assert tree[1].right == 0
        assert tree[1].distance == pytest.approx(7.300, abs=0.0005)
        assert tree[2].left == 3
        assert tree[2].right == -2
        assert tree[2].distance == pytest.approx(13.540, abs=0.0005)
        indices = tree.cut(nclusters=1)
        assert len(indices) == len(data1)
        assert indices[0] == 0
        assert indices[1] == 0
        assert indices[2] == 0
        assert indices[3] == 0
        indices = tree.cut(nclusters=2)
        assert len(indices) == len(data1)
        assert indices[0] == 1
        assert indices[1] == 1
        assert indices[2] == 1
        assert indices[3] == 0
        indices = tree.cut(nclusters=3)
        assert len(indices) == len(data1)
        assert indices[0] == 2
        assert indices[1] == 1
        assert indices[2] == 1
        assert indices[3] == 0
        indices = tree.cut(nclusters=4)
        assert len(indices) == len(data1)
        assert indices[0] == 3
        assert indices[1] == 2
        assert indices[2] == 1
        assert indices[3] == 0
        indices = tree.sort([0, 1, 2, 3])
        assert len(indices) == len(data1)
        assert indices[0] == 0
        assert indices[1] == 1
        assert indices[2] == 2
        assert indices[3] == 3
        indices = tree.sort([0, 3, 2, 1])
        assert len(indices) == len(data1)
        assert indices[0] == 3
        assert indices[1] == 0
        assert indices[2] == 2
        assert indices[3] == 1

        # Pairwise single-linkage clustering
        tree = treecluster(
            data=data1,
            mask=mask1,
            weight=weight1,
            transpose=False,
            method="s",
            dist="e",
        )
        assert len(tree) == len(data1) - 1
        assert tree[0].left == 1
        assert tree[0].right == 2
        assert tree[0].distance == pytest.approx(2.600, abs=0.0005)
        assert tree[1].left == 0
        assert tree[1].right == -1
        assert tree[1].distance == pytest.approx(5.800, abs=0.0005)
        assert tree[2].left == -2
        assert tree[2].right == 3
        assert tree[2].distance == pytest.approx(6.380, abs=0.0005)
        indices = tree.cut(nclusters=1)
        assert len(indices) == len(data1)
        assert indices[0] == 0
        assert indices[1] == 0
        assert indices[2] == 0
        assert indices[3] == 0
        indices = tree.cut(nclusters=2)
        assert len(indices) == len(data1)
        assert indices[0] == 0
        assert indices[1] == 0
        assert indices[2] == 0
        assert indices[3] == 1
        indices = tree.cut(nclusters=3)
        assert len(indices) == len(data1)
        assert indices[0] == 0
        assert indices[1] == 1
        assert indices[2] == 1
        assert indices[3] == 2
        indices = tree.cut(nclusters=4)
        assert len(indices) == len(data1)
        assert indices[0] == 0
        assert indices[1] == 1
        assert indices[2] == 2
        assert indices[3] == 3
        indices = tree.sort([0, 1, 2, 3])
        assert len(indices) == len(data1)
        assert indices[0] == 0
        assert indices[1] == 1
        assert indices[2] == 2
        assert indices[3] == 3
        indices = tree.sort([0, 3, 2, 1])
        assert len(indices) == len(data1)
        assert indices[0] == 3
        assert indices[1] == 0
        assert indices[2] == 2
        assert indices[3] == 1

        # Pairwise centroid-linkage clustering
        tree = treecluster(
            data=data1,
            mask=mask1,
            weight=weight1,
            transpose=False,
            method="c",
            dist="e",
        )
        assert len(tree) == len(data1) - 1
        assert tree[0].left == 1
        assert tree[0].right == 2
        assert tree[0].distance == pytest.approx(2.600, abs=0.0005)
        assert tree[1].left == 0
        assert tree[1].right == -1
        assert tree[1].distance == pytest.approx(6.650, abs=0.0005)
        assert tree[2].left == -2
        assert tree[2].right == 3
        assert tree[2].distance == pytest.approx(11.629, abs=0.0005)
        indices = tree.sort([0, 1, 2, 3])
        assert len(indices) == len(data1)
        assert indices[0] == 0
        assert indices[1] == 1
        assert indices[2] == 2
        assert indices[3] == 3
        indices = tree.cut(nclusters=1)
        assert len(indices) == len(data1)
        assert indices[0] == 0
        assert indices[1] == 0
        assert indices[2] == 0
        assert indices[3] == 0
        indices = tree.cut(nclusters=2)
        assert len(indices) == len(data1)
        assert indices[0] == 0
        assert indices[1] == 0
        assert indices[2] == 0
        assert indices[3] == 1
        indices = tree.cut(nclusters=3)
        assert len(indices) == len(data1)
        assert indices[0] == 0
        assert indices[1] == 1
        assert indices[2] == 1
        assert indices[3] == 2
        indices = tree.cut(nclusters=4)
        assert len(indices) == len(data1)
        assert indices[0] == 0
        assert indices[1] == 1
        assert indices[2] == 2
        assert indices[3] == 3
        indices = tree.sort([0, 3, 2, 1])
        assert len(indices) == len(data1)
        assert indices[0] == 3
        assert indices[1] == 0
        assert indices[2] == 2
        assert indices[3] == 1

        # Pairwise maximum-linkage clustering
        tree = treecluster(
            data=data1,
            mask=mask1,
            weight=weight1,
            transpose=False,
            method="m",
            dist="e",
        )
        assert len(tree) == len(data1) - 1
        assert tree[0].left == 2
        assert tree[0].right == 1
        assert tree[0].distance == pytest.approx(2.600, abs=0.0005)
        assert tree[1].left == -1
        assert tree[1].right == 0
        assert tree[1].distance == pytest.approx(8.800, abs=0.0005)
        assert tree[2].left == 3
        assert tree[2].right == -2
        assert tree[2].distance == pytest.approx(23.100, abs=0.0005)
        indices = tree.cut(nclusters=1)
        assert len(indices) == len(data1)
        assert indices[0] == 0
        assert indices[1] == 0
        assert indices[2] == 0
        assert indices[3] == 0
        indices = tree.cut(nclusters=2)
        assert len(indices) == len(data1)
        assert indices[0] == 1
        assert indices[1] == 1
        assert indices[2] == 1
        assert indices[3] == 0
        indices = tree.cut(nclusters=3)
        assert len(indices) == len(data1)
        assert indices[0] == 2
        assert indices[1] == 1
        assert indices[2] == 1
        assert indices[3] == 0
        indices = tree.cut(nclusters=4)
        assert len(indices) == len(data1)
        assert indices[0] == 3
        assert indices[1] == 2
        assert indices[2] == 1
        assert indices[3] == 0
        indices = tree.sort([0, 1, 2, 3])
        assert len(indices) == len(data1)
        assert indices[0] == 0
        assert indices[1] == 1
        assert indices[2] == 2
        assert indices[3] == 3
        indices = tree.sort([0, 3, 2, 1])
        assert len(indices) == len(data1)
        assert indices[0] == 3
        assert indices[1] == 0
        assert indices[2] == 2
        assert indices[3] == 1

        # First data set, using transpose=True
        weight1 = [1, 1, 1, 1]
        data1 = np.array(
            [
                [1.1, 2.2, 3.3, 4.4, 5.5],
                [3.1, 3.2, 1.3, 2.4, 1.5],
                [4.1, 2.2, 0.3, 5.4, 0.5],
                [9.7, 2.0, 0.0, 5.0, 0.0],
            ]
        )
        # fmt: off
        mask1 = np.array(
            [
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
            ],
            int,
        )
        # fmt: on
        nrows, ncols = data1.shape

        # Pairwise average-linkage clustering
        tree = treecluster(
            data=data1, mask=mask1, weight=weight1, transpose=True, method="a", dist="e"
        )
        assert len(tree) == ncols - 1
        assert tree[0].left == 4
        assert tree[0].right == 2
        assert tree[0].distance == pytest.approx(1.230, abs=0.0005)
        assert tree[1].left == -1
        assert tree[1].right == 1
        assert tree[1].distance == pytest.approx(4.1375, abs=0.0005)
        assert tree[2].left == 3
        assert tree[2].right == 0
        assert tree[2].distance == pytest.approx(8.790, abs=0.0005)
        assert tree[3].left == -2
        assert tree[3].right == -3
        assert tree[3].distance == pytest.approx(18.2867, abs=0.0005)
        indices = tree.cut(nclusters=1)
        assert len(indices) == ncols
        assert indices[0] == 0
        assert indices[1] == 0
        assert indices[2] == 0
        assert indices[3] == 0
        assert indices[4] == 0
        indices = tree.cut(nclusters=2)
        assert len(indices) == ncols
        assert indices[0] == 1
        assert indices[1] == 0
        assert indices[2] == 0
        assert indices[3] == 1
        assert indices[4] == 0
        indices = tree.cut(nclusters=3)
        assert len(indices) == ncols
        assert indices[0] == 2
        assert indices[1] == 0
        assert indices[2] == 0
        assert indices[3] == 1
        assert indices[4] == 0
        indices = tree.cut(nclusters=4)
        assert len(indices) == ncols
        assert indices[0] == 3
        assert indices[1] == 1
        assert indices[2] == 0
        assert indices[3] == 2
        assert indices[4] == 0
        indices = tree.sort([0, 1, 2, 3, 4])
        assert len(indices) == ncols
        assert indices[0] == 0
        assert indices[1] == 3
        assert indices[2] == 1
        assert indices[3] == 2
        assert indices[4] == 4
        indices = tree.sort([0, 4, 3, 2, 1])
        assert len(indices) == ncols
        assert indices[0] == 0
        assert indices[1] == 3
        assert indices[2] == 4
        assert indices[3] == 2
        assert indices[4] == 1

        # Pairwise single-linkage clustering
        tree = treecluster(
            data=data1, mask=mask1, weight=weight1, transpose=True, method="s", dist="e"
        )
        assert len(tree) == ncols - 1
        assert tree[0].left == 2
        assert tree[0].right == 4
        assert tree[0].distance == pytest.approx(1.230, abs=0.0005)
        assert tree[1].left == 1
        assert tree[1].right == -1
        assert tree[1].distance == pytest.approx(3.1075, abs=0.0005)
        assert tree[2].left == 3
        assert tree[2].right == -2
        assert tree[2].distance == pytest.approx(6.180, abs=0.0005)
        assert tree[3].left == 0
        assert tree[3].right == -3
        assert tree[3].distance == pytest.approx(8.790, abs=0.0005)
        indices = tree.cut(nclusters=1)
        assert len(indices) == ncols
        assert indices[0] == 0
        assert indices[1] == 0
        assert indices[2] == 0
        assert indices[3] == 0
        assert indices[4] == 0
        indices = tree.cut(nclusters=2)
        assert len(indices) == ncols
        assert indices[0] == 0
        assert indices[1] == 1
        assert indices[2] == 1
        assert indices[3] == 1
        assert indices[4] == 1
        indices = tree.cut(nclusters=3)
        assert len(indices) == ncols
        assert indices[0] == 0
        assert indices[1] == 2
        assert indices[2] == 2
        assert indices[3] == 1
        assert indices[4] == 2
        indices = tree.cut(nclusters=4)
        assert len(indices) == ncols
        assert indices[0] == 0
        assert indices[1] == 2
        assert indices[2] == 3
        assert indices[3] == 1
        assert indices[4] == 3
        indices = tree.sort([0, 1, 2, 3, 4])
        assert len(indices) == ncols
        assert indices[0] == 0
        assert indices[1] == 1
        assert indices[2] == 2
        assert indices[3] == 4
        assert indices[4] == 3
        indices = tree.sort([0, 4, 3, 2, 1])
        assert len(indices) == ncols
        assert indices[0] == 0
        assert indices[1] == 3
        assert indices[2] == 4
        assert indices[3] == 2
        assert indices[4] == 1

        # Pairwise centroid-linkage clustering
        tree = treecluster(
            data=data1, mask=mask1, weight=weight1, transpose=True, method="c", dist="e"
        )
        assert len(tree) == ncols - 1
        assert tree[0].left == 2
        assert tree[0].right == 4
        assert tree[0].distance == pytest.approx(1.23, abs=0.0005)
        assert tree[1].left == 1
        assert tree[1].right == -1
        assert tree[1].distance == pytest.approx(3.83, abs=0.0005)
        assert tree[2].left == 0
        assert tree[2].right == 3
        assert tree[2].distance == pytest.approx(8.79, abs=0.0005)
        assert tree[3].left == -3
        assert tree[3].right == -2
        assert tree[3].distance == pytest.approx(15.0331, abs=0.0005)
        indices = tree.sort([0, 1, 2, 3, 4])
        assert len(indices) == ncols
        assert indices[0] == 0
        assert indices[1] == 3
        assert indices[2] == 1
        assert indices[3] == 2
        assert indices[4] == 4
        indices = tree.cut(nclusters=1)
        assert len(indices) == ncols
        assert indices[0] == 0
        assert indices[1] == 0
        assert indices[2] == 0
        assert indices[3] == 0
        assert indices[4] == 0
        indices = tree.cut(nclusters=2)
        assert len(indices) == ncols
        assert indices[0] == 0
        assert indices[1] == 1
        assert indices[2] == 1
        assert indices[3] == 0
        assert indices[4] == 1
        indices = tree.cut(nclusters=3)
        assert len(indices) == ncols
        assert indices[0] == 0
        assert indices[1] == 2
        assert indices[2] == 2
        assert indices[3] == 1
        assert indices[4] == 2
        indices = tree.cut(nclusters=4)
        assert len(indices) == ncols
        assert indices[0] == 0
        assert indices[1] == 2
        assert indices[2] == 3
        assert indices[3] == 1
        assert indices[4] == 3
        indices = tree.sort([0, 4, 3, 2, 1])
        assert len(indices) == ncols
        assert indices[0] == 0
        assert indices[1] == 3
        assert indices[2] == 4
        assert indices[3] == 2
        assert indices[4] == 1

        # Pairwise maximum-linkage clustering
        tree = treecluster(
            data=data1, mask=mask1, weight=weight1, transpose=True, method="m", dist="e"
        )
        assert len(tree) == ncols - 1
        assert tree[0].left == 4
        assert tree[0].right == 2
        assert tree[0].distance == pytest.approx(1.230, abs=0.0005)
        assert tree[1].left == -1
        assert tree[1].right == 1
        assert tree[1].distance == pytest.approx(5.1675, abs=0.0005)
        assert tree[2].left == 3
        assert tree[2].right == 0
        assert tree[2].distance == pytest.approx(8.790, abs=0.0005)
        assert tree[3].left == -2
        assert tree[3].right == -3
        assert tree[3].distance == pytest.approx(32.2425, abs=0.0005)
        indices = tree.cut(nclusters=1)
        assert len(indices) == ncols
        assert indices[0] == 0
        assert indices[1] == 0
        assert indices[2] == 0
        assert indices[3] == 0
        assert indices[4] == 0
        indices = tree.cut(nclusters=2)
        assert len(indices) == ncols
        assert indices[0] == 1
        assert indices[1] == 0
        assert indices[2] == 0
        assert indices[3] == 1
        assert indices[4] == 0
        indices = tree.cut(nclusters=3)
        assert len(indices) == ncols
        assert indices[0] == 2
        assert indices[1] == 0
        assert indices[2] == 0
        assert indices[3] == 1
        assert indices[4] == 0
        indices = tree.cut(nclusters=4)
        assert len(indices) == ncols
        assert indices[0] == 3
        assert indices[1] == 1
        assert indices[2] == 0
        assert indices[3] == 2
        assert indices[4] == 0
        indices = tree.sort([0, 1, 2, 3, 4])
        assert len(indices) == ncols
        assert indices[0] == 0
        assert indices[1] == 3
        assert indices[2] == 1
        assert indices[3] == 2
        assert indices[4] == 4
        indices = tree.sort([0, 4, 3, 2, 1])
        assert len(indices) == ncols
        assert indices[0] == 0
        assert indices[1] == 3
        assert indices[2] == 4
        assert indices[3] == 2
        assert indices[4] == 1

        # Second data set
        weight2 = [1, 1]
        data2 = np.array(
            [
                [0.8223, 0.9295],
                [1.4365, 1.3223],
                [1.1623, 1.5364],
                [2.1826, 1.1934],
                [1.7763, 1.9352],
                [1.7215, 1.9912],
                [2.1812, 5.9935],
                [5.3290, 5.9452],
                [3.1491, 3.3454],
                [5.1923, 5.3156],
                [4.7735, 5.4012],
                [5.1297, 5.5645],
                [5.3934, 5.1823],
            ]
        )
        mask2 = np.array(
            [
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
            ],
            int,
        )

        # Test second data set
        # Pairwise average-linkage clustering
        tree = treecluster(
            data=data2,
            mask=mask2,
            weight=weight2,
            transpose=False,
            method="a",
            dist="e",
        )
        assert len(tree) == len(data2) - 1
        assert tree[0].left == 5
        assert tree[0].right == 4
        assert tree[0].distance == pytest.approx(0.003, abs=0.0005)
        assert tree[1].left == 9
        assert tree[1].right == 12
        assert tree[1].distance == pytest.approx(0.029, abs=0.0005)
        assert tree[2].left == 2
        assert tree[2].right == 1
        assert tree[2].distance == pytest.approx(0.061, abs=0.0005)
        assert tree[3].left == 11
        assert tree[3].right == -2
        assert tree[3].distance == pytest.approx(0.070, abs=0.0005)
        assert tree[4].left == -4
        assert tree[4].right == 10
        assert tree[4].distance == pytest.approx(0.128, abs=0.0005)
        assert tree[5].left == 7
        assert tree[5].right == -5
        assert tree[5].distance == pytest.approx(0.224, abs=0.0005)
        assert tree[6].left == -3
        assert tree[6].right == 0
        assert tree[6].distance == pytest.approx(0.254, abs=0.0005)
        assert tree[7].left == -1
        assert tree[7].right == 3
        assert tree[7].distance == pytest.approx(0.391, abs=0.0005)
        assert tree[8].left == -8
        assert tree[8].right == -7
        assert tree[8].distance == pytest.approx(0.532, abs=0.0005)
        assert tree[9].left == 8
        assert tree[9].right == -9
        assert tree[9].distance == pytest.approx(3.234, abs=0.0005)
        assert tree[10].left == -6
        assert tree[10].right == 6
        assert tree[10].distance == pytest.approx(4.636, abs=0.0005)
        assert tree[11].left == -11
        assert tree[11].right == -10
        assert tree[11].distance == pytest.approx(12.741, abs=0.0005)
        indices = tree.cut(nclusters=1)
        assert len(indices) == len(data2)
        assert indices[0] == 0
        assert indices[1] == 0
        assert indices[2] == 0
        assert indices[3] == 0
        assert indices[4] == 0
        assert indices[5] == 0
        assert indices[6] == 0
        assert indices[7] == 0
        assert indices[8] == 0
        assert indices[9] == 0
        assert indices[10] == 0
        assert indices[11] == 0
        assert indices[12] == 0
        indices = tree.cut(nclusters=2)
        assert len(indices) == len(data2)
        assert indices[0] == 1
        assert indices[1] == 1
        assert indices[2] == 1
        assert indices[3] == 1
        assert indices[4] == 1
        assert indices[5] == 1
        assert indices[6] == 0
        assert indices[7] == 0
        assert indices[8] == 1
        assert indices[9] == 0
        assert indices[10] == 0
        assert indices[11] == 0
        assert indices[12] == 0
        indices = tree.cut(nclusters=3)
        assert len(indices) == len(data2)
        assert indices[0] == 2
        assert indices[1] == 2
        assert indices[2] == 2
        assert indices[3] == 2
        assert indices[4] == 2
        assert indices[5] == 2
        assert indices[6] == 1
        assert indices[7] == 0
        assert indices[8] == 2
        assert indices[9] == 0
        assert indices[10] == 0
        assert indices[11] == 0
        assert indices[12] == 0
        indices = tree.cut(nclusters=4)
        assert len(indices) == len(data2)
        assert indices[0] == 3
        assert indices[1] == 3
        assert indices[2] == 3
        assert indices[3] == 3
        assert indices[4] == 3
        assert indices[5] == 3
        assert indices[6] == 1
        assert indices[7] == 0
        assert indices[8] == 2
        assert indices[9] == 0
        assert indices[10] == 0
        assert indices[11] == 0
        assert indices[12] == 0
        indices = tree.cut(nclusters=5)
        assert len(indices) == len(data2)
        assert indices[0] == 4
        assert indices[1] == 4
        assert indices[2] == 4
        assert indices[3] == 3
        assert indices[4] == 3
        assert indices[5] == 3
        assert indices[6] == 1
        assert indices[7] == 0
        assert indices[8] == 2
        assert indices[9] == 0
        assert indices[10] == 0
        assert indices[11] == 0
        assert indices[12] == 0
        indices = tree.sort()
        assert len(indices) == len(data2)
        assert indices[0] == 7
        assert indices[1] == 11
        assert indices[2] == 9
        assert indices[3] == 12
        assert indices[4] == 10
        assert indices[5] == 6
        assert indices[6] == 8
        assert indices[7] == 5
        assert indices[8] == 4
        assert indices[9] == 3
        assert indices[10] == 2
        assert indices[11] == 1
        assert indices[12] == 0

        # Pairwise single-linkage clustering
        tree = treecluster(
            data=data2,
            mask=mask2,
            weight=weight2,
            transpose=False,
            method="s",
            dist="e",
        )
        assert len(tree) == len(data2) - 1
        assert tree[0].left == 4
        assert tree[0].right == 5
        assert tree[0].distance == pytest.approx(0.003, abs=0.0005)
        assert tree[1].left == 9
        assert tree[1].right == 12
        assert tree[1].distance == pytest.approx(0.029, abs=0.0005)
        assert tree[2].left == 11
        assert tree[2].right == -2
        assert tree[2].distance == pytest.approx(0.033, abs=0.0005)
        assert tree[3].left == 1
        assert tree[3].right == 2
        assert tree[3].distance == pytest.approx(0.061, abs=0.0005)
        assert tree[4].left == 10
        assert tree[4].right == -3
        assert tree[4].distance == pytest.approx(0.077, abs=0.0005)
        assert tree[5].left == 7
        assert tree[5].right == -5
        assert tree[5].distance == pytest.approx(0.092, abs=0.0005)
        assert tree[6].left == 0
        assert tree[6].right == -4
        assert tree[6].distance == pytest.approx(0.242, abs=0.0005)
        assert tree[7].left == -7
        assert tree[7].right == -1
        assert tree[7].distance == pytest.approx(0.246, abs=0.0005)
        assert tree[8].left == 3
        assert tree[8].right == -8
        assert tree[8].distance == pytest.approx(0.287, abs=0.0005)
        assert tree[9].left == -9
        assert tree[9].right == 8
        assert tree[9].distance == pytest.approx(1.936, abs=0.0005)
        assert tree[10].left == -10
        assert tree[10].right == -6
        assert tree[10].distance == pytest.approx(3.432, abs=0.0005)
        assert tree[11].left == 6
        assert tree[11].right == -11
        assert tree[11].distance == pytest.approx(3.535, abs=0.0005)
        indices = tree.sort()
        assert len(indices) == len(data2)
        assert indices[0] == 6
        assert indices[1] == 3
        assert indices[2] == 0
        assert indices[3] == 1
        assert indices[4] == 2
        assert indices[5] == 4
        assert indices[6] == 5
        assert indices[7] == 8
        assert indices[8] == 7
        assert indices[9] == 10
        assert indices[10] == 11
        assert indices[11] == 9
        assert indices[12] == 12

        # Pairwise centroid-linkage clustering
        tree = treecluster(
            data=data2,
            mask=mask2,
            weight=weight2,
            transpose=False,
            method="c",
            dist="e",
        )
        assert len(tree) == len(data2) - 1
        assert tree[0].left == 4
        assert tree[0].right == 5
        assert tree[0].distance == pytest.approx(0.003, abs=0.0005)
        assert tree[1].left == 12
        assert tree[1].right == 9
        assert tree[1].distance == pytest.approx(0.029, abs=0.0005)
        assert tree[2].left == 1
        assert tree[2].right == 2
        assert tree[2].distance == pytest.approx(0.061, abs=0.0005)
        assert tree[3].left == -2
        assert tree[3].right == 11
        assert tree[3].distance == pytest.approx(0.063, abs=0.0005)
        assert tree[4].left == 10
        assert tree[4].right == -4
        assert tree[4].distance == pytest.approx(0.109, abs=0.0005)
        assert tree[5].left == -5
        assert tree[5].right == 7
        assert tree[5].distance == pytest.approx(0.189, abs=0.0005)
        assert tree[6].left == 0
        assert tree[6].right == -3
        assert tree[6].distance == pytest.approx(0.239, abs=0.0005)
        assert tree[7].left == 3
        assert tree[7].right == -1
        assert tree[7].distance == pytest.approx(0.390, abs=0.0005)
        assert tree[8].left == -7
        assert tree[8].right == -8
        assert tree[8].distance == pytest.approx(0.382, abs=0.0005)
        assert tree[9].left == -9
        assert tree[9].right == 8
        assert tree[9].distance == pytest.approx(3.063, abs=0.0005)
        assert tree[10].left == 6
        assert tree[10].right == -6
        assert tree[10].distance == pytest.approx(4.578, abs=0.0005)
        assert tree[11].left == -10
        assert tree[11].right == -11
        assert tree[11].distance == pytest.approx(11.536, abs=0.0005)
        indices = tree.sort()
        assert len(indices) == len(data2)
        assert indices[0] == 0
        assert indices[1] == 1
        assert indices[2] == 2
        assert indices[3] == 3
        assert indices[4] == 4
        assert indices[5] == 5
        assert indices[6] == 8
        assert indices[7] == 6
        assert indices[8] == 10
        assert indices[9] == 12
        assert indices[10] == 9
        assert indices[11] == 11
        assert indices[12] == 7

        # Pairwise maximum-linkage clustering
        tree = treecluster(
            data=data2,
            mask=mask2,
            weight=weight2,
            transpose=False,
            method="m",
            dist="e",
        )
        assert len(tree) == len(data2) - 1
        assert tree[0].left == 5
        assert tree[0].right == 4
        assert tree[0].distance == pytest.approx(0.003, abs=0.0005)
        assert tree[1].left == 9
        assert tree[1].right == 12
        assert tree[1].distance == pytest.approx(0.029, abs=0.0005)
        assert tree[2].left == 2
        assert tree[2].right == 1
        assert tree[2].distance == pytest.approx(0.061, abs=0.0005)
        assert tree[3].left == 11
        assert tree[3].right == 10
        assert tree[3].distance == pytest.approx(0.077, abs=0.0005)
        assert tree[4].left == -2
        assert tree[4].right == -4
        assert tree[4].distance == pytest.approx(0.216, abs=0.0005)
        assert tree[5].left == -3
        assert tree[5].right == 0
        assert tree[5].distance == pytest.approx(0.266, abs=0.0005)
        assert tree[6].left == -5
        assert tree[6].right == 7
        assert tree[6].distance == pytest.approx(0.302, abs=0.0005)
        assert tree[7].left == -1
        assert tree[7].right == 3
        assert tree[7].distance == pytest.approx(0.425, abs=0.0005)
        assert tree[8].left == -8
        assert tree[8].right == -6
        assert tree[8].distance == pytest.approx(0.968, abs=0.0005)
        assert tree[9].left == 8
        assert tree[9].right == 6
        assert tree[9].distance == pytest.approx(3.975, abs=0.0005)
        assert tree[10].left == -10
        assert tree[10].right == -7
        assert tree[10].distance == pytest.approx(5.755, abs=0.0005)
        assert tree[11].left == -11
        assert tree[11].right == -9
        assert tree[11].distance == pytest.approx(22.734, abs=0.0005)
        indices = tree.sort()
        assert len(indices) == len(data2)
        assert indices[0] == 8
        assert indices[1] == 6
        assert indices[2] == 9
        assert indices[3] == 12
        assert indices[4] == 11
        assert indices[5] == 10
        assert indices[6] == 7
        assert indices[7] == 5
        assert indices[8] == 4
        assert indices[9] == 3
        assert indices[10] == 2
        assert indices[11] == 1
        assert indices[12] == 0

    def test_somcluster_arguments(self):
        # Test if incorrect arguments are caught by the C code
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster._cluster import somcluster
        elif TestCluster.module == "Pycluster":
            from Pycluster._cluster import somcluster

        weight = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
        data = np.array(
            [
                [1.1, 2.2, 3.3, 4.4, 5.5],
                [3.1, 3.2, 1.3, 2.4, 1.5],
                [4.1, 2.2, 0.3, 5.4, 0.5],
                [9.9, 2.0, 0.0, 5.0, 0.0],
            ]
        )
        # fmt: off
        mask = np.array(
            [
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
            ],
            np.int32,
        )
        # fmt: on
        nitems, ndata = data.shape
        nxgrid, nygrid = 10, 10
        clusterids = np.ones((nitems, 2), np.int32)
        celldata = np.zeros((nxgrid, nygrid, ndata), dtype="d")

        message = "^unexpected format.$"
        with pytest.raises(RuntimeError, match=message):
            somcluster(
                clusterids=None,
                celldata=celldata,
                data=data,
                mask=mask,
                weight=weight,
                transpose=False,
                inittau=0.02,
                niter=100,
                dist="e",
            )
        message = "^incorrect rank 1 \\(expected 2\\)$"
        with pytest.raises(ValueError, match=message):
            somcluster(
                clusterids=np.ones(nitems, np.int32),
                celldata=celldata,
                data=data,
                mask=mask,
                weight=weight,
                transpose=False,
                inittau=0.02,
                niter=100,
                dist="e",
            )
        message = "^argument has incorrect data type$"
        with pytest.raises(RuntimeError, match=message):
            somcluster(
                clusterids=np.ones((nitems, 2), np.int16),
                celldata=celldata,
                data=data,
                mask=mask,
                weight=weight,
                transpose=False,
                inittau=0.02,
                niter=100,
                dist="e",
            )
        message = "^array has 3 columns \\(expected 2\\)$"
        with pytest.raises(ValueError, match=message):
            somcluster(
                clusterids=np.ones((nitems, 3), np.int32),
                celldata=celldata,
                data=data,
                mask=mask,
                weight=weight,
                transpose=False,
                inittau=0.02,
                niter=100,
                dist="e",
            )
        message = "^celldata array has unexpected format.$"
        with pytest.raises(RuntimeError, match=message):
            somcluster(
                clusterids=clusterids,
                celldata=None,
                data=data,
                mask=mask,
                weight=weight,
                transpose=False,
                inittau=0.02,
                niter=100,
                dist="e",
            )
        message = "^celldata array has incorrect data type$"
        with pytest.raises(RuntimeError, match=message):
            somcluster(
                clusterids=clusterids,
                celldata=np.zeros((nxgrid, nygrid, ndata), dtype=np.int32),
                data=data,
                mask=mask,
                weight=weight,
                transpose=False,
                inittau=0.02,
                niter=100,
                dist="e",
            )
        message = "^data is None$"
        with pytest.raises(RuntimeError, match=message):
            somcluster(
                clusterids=clusterids,
                celldata=celldata,
                data=None,
                mask=mask,
                weight=weight,
                transpose=False,
                inittau=0.02,
                niter=100,
                dist="e",
            )
        message = "^data matrix has unexpected format.$"
        with pytest.raises(RuntimeError, match=message):
            somcluster(
                clusterids=clusterids,
                celldata=celldata,
                data=[None],
                mask=mask,
                weight=weight,
                transpose=False,
                inittau=0.02,
                niter=100,
                dist="e",
            )
        message = "^data matrix has incorrect data type$"
        with pytest.raises(RuntimeError, match=message):
            somcluster(
                clusterids=clusterids,
                celldata=celldata,
                data=np.zeros((4, 5), dtype=np.int16),
                mask=mask,
                weight=weight,
                transpose=False,
                inittau=0.02,
                niter=100,
                dist="e",
            )
        message = "^data matrix has incorrect rank 1 \\(expected 2\\)$"
        with pytest.raises(RuntimeError, match=message):
            somcluster(
                clusterids=clusterids,
                celldata=celldata,
                data=np.zeros(4),
                mask=mask,
                weight=weight,
                transpose=False,
                inittau=0.02,
                niter=100,
                dist="e",
            )
        message = "^data matrix is empty$"
        with pytest.raises(ValueError, match=message):
            somcluster(
                clusterids=clusterids,
                celldata=celldata,
                data=data[:0],
                mask=mask,
                weight=weight,
                transpose=False,
                inittau=0.02,
                niter=100,
                dist="e",
            )
        message = "^data is not contiguous$"
        with pytest.raises(RuntimeError, match=message):
            somcluster(
                clusterids=clusterids,
                celldata=celldata,
                data=data[:, ::2],
                mask=mask,
                weight=weight,
                transpose=False,
                inittau=0.02,
                niter=100,
                dist="e",
            )
        message = "^mask is None$"
        with pytest.raises(RuntimeError, match=message):
            somcluster(
                clusterids=clusterids,
                celldata=celldata,
                data=data,
                mask=None,
                weight=weight,
                transpose=False,
                inittau=0.02,
                niter=100,
                dist="e",
            )
        message = "^mask has unexpected format.$"
        with pytest.raises(RuntimeError, match=message):
            somcluster(
                clusterids=clusterids,
                celldata=celldata,
                data=data,
                mask=[None],
                weight=weight,
                transpose=False,
                inittau=0.02,
                niter=100,
                dist="e",
            )
        message = "^mask has incorrect rank 1 \\(expected 2\\)$"
        with pytest.raises(ValueError, match=message):
            somcluster(
                clusterids=clusterids,
                celldata=celldata,
                data=data,
                mask=np.array([1, 1, 1]),
                weight=weight,
                transpose=False,
                inittau=0.02,
                niter=100,
                dist="e",
            )
        message = "^mask has incorrect data type$"
        with pytest.raises(RuntimeError, match=message):
            somcluster(
                clusterids=clusterids,
                celldata=celldata,
                data=data,
                mask=np.array([[1, 1], [1, 1]], dtype=np.int16),
                weight=weight,
                transpose=False,
                inittau=0.02,
                niter=100,
                dist="e",
            )
        message = "^mask is not contiguous$"
        with pytest.raises(RuntimeError, match=message):
            somcluster(
                clusterids=clusterids,
                celldata=celldata,
                data=data,
                mask=mask[:, ::2],
                weight=weight,
                transpose=False,
                inittau=0.02,
                niter=100,
                dist="e",
            )
        message = "^unexpected format.$"
        with pytest.raises(RuntimeError, match=message):
            somcluster(
                clusterids=clusterids,
                celldata=celldata,
                data=data,
                mask=mask,
                weight=None,
                transpose=False,
                inittau=0.02,
                niter=100,
                dist="e",
            )
        message = "^incorrect rank 3 \\(expected 1\\)$"
        with pytest.raises(ValueError, match=message):
            somcluster(
                clusterids=clusterids,
                celldata=celldata,
                data=data,
                mask=mask,
                weight=np.zeros((2, 2, 2)),
                transpose=False,
                inittau=0.02,
                niter=100,
                dist="e",
            )
        message = "^array has incorrect data type$"
        with pytest.raises(RuntimeError, match=message):
            somcluster(
                clusterids=clusterids,
                celldata=celldata,
                data=data,
                mask=mask,
                weight=np.array([1, 1, 1], dtype=np.int16),
                transpose=False,
                inittau=0.02,
                niter=100,
                dist="e",
            )
        message = "^dist should be a string$"
        with pytest.raises(ValueError, match=message):
            somcluster(
                clusterids=clusterids,
                celldata=celldata,
                data=data,
                mask=mask,
                weight=weight,
                transpose=False,
                inittau=0.02,
                niter=100,
                dist=5,
            )
        message = "^dist should be a single character$"
        with pytest.raises(ValueError, match=message):
            somcluster(
                clusterids=clusterids,
                celldata=celldata,
                data=data,
                mask=mask,
                weight=weight,
                transpose=False,
                inittau=0.02,
                niter=100,
                dist="Pearson",
            )
        message = "^unknown dist function specified \\(should be one of 'ebcauxsk'\\)$"
        with pytest.raises(ValueError, match=message):
            somcluster(
                clusterids=clusterids,
                celldata=celldata,
                data=data,
                mask=mask,
                weight=weight,
                transpose=False,
                inittau=0.02,
                niter=100,
                dist="X",
            )

    def test_somcluster(self):
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster import somcluster
        elif TestCluster.module == "Pycluster":
            from Pycluster import somcluster

        # First data set
        weight = [1, 1, 1, 1, 1]
        data = np.array(
            [
                [1.1, 2.2, 3.3, 4.4, 5.5],
                [3.1, 3.2, 1.3, 2.4, 1.5],
                [4.1, 2.2, 0.3, 5.4, 0.5],
                [9.9, 2.0, 0.0, 5.0, 0.0],
            ]
        )
        # fmt: off
        mask = np.array(
            [
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
            ],
            int,
        )
        # fmt: on
        nrows, ncols = data.shape

        clusterid, celldata = somcluster(
            data=data,
            mask=mask,
            weight=weight,
            transpose=False,
            nxgrid=10,
            nygrid=10,
            inittau=0.02,
            niter=100,
            dist="e",
        )
        assert len(clusterid) == nrows
        assert len(clusterid[0]) == 2

        # First data set, using transpose=True
        weight = [1, 1, 1, 1]
        clusterid, celldata = somcluster(
            data=data,
            mask=mask,
            weight=weight,
            transpose=True,
            nxgrid=10,
            nygrid=10,
            inittau=0.02,
            niter=100,
            dist="e",
        )
        assert len(clusterid) == ncols
        assert len(clusterid[0]) == 2

        # Second data set
        weight = [1, 1]
        data = np.array(
            [
                [1.1, 1.2],
                [1.4, 1.3],
                [1.1, 1.5],
                [2.0, 1.5],
                [1.7, 1.9],
                [1.7, 1.9],
                [5.7, 5.9],
                [5.7, 5.9],
                [3.1, 3.3],
                [5.4, 5.3],
                [5.1, 5.5],
                [5.0, 5.5],
                [5.1, 5.2],
            ]
        )
        mask = np.array(
            [
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
            ],
            int,
        )

        clusterid, celldata = somcluster(
            data=data,
            mask=mask,
            weight=weight,
            transpose=False,
            nxgrid=10,
            nygrid=10,
            inittau=0.02,
            niter=100,
            dist="e",
        )
        assert len(clusterid) == len(data)
        assert len(clusterid[0]) == 2

    def test_distancematrix_arguments(self):
        # Test if incorrect arguments are caught by the C code
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster._cluster import distancematrix
        elif TestCluster.module == "Pycluster":
            from Pycluster._cluster import distancematrix

        data = np.array(
            [
                [2.2, 3.3, 4.4],
                [2.1, 1.4, 5.6],
                [7.8, 9.0, 1.2],
                [4.5, 2.3, 1.5],
                [4.2, 2.4, 1.9],
                [3.6, 3.1, 9.3],
                [2.3, 1.2, 3.9],
                [4.2, 9.6, 9.3],
                [1.7, 8.9, 1.1],
            ]
        )
        mask = np.array(
            [
                [1, 1, 1],
                [1, 1, 1],
                [0, 1, 1],
                [1, 1, 1],
                [1, 1, 1],
                [0, 1, 0],
                [1, 1, 1],
                [1, 0, 1],
                [1, 1, 1],
            ],
            np.int32,
        )
        weight = np.array([2.0, 1.0, 0.5])
        message = "^data matrix is empty$"
        with pytest.raises(ValueError, match=message):
            distancematrix(data[:0, :], mask=mask, weight=weight)
        message = "^mask has incorrect rank 1 \\(expected 2\\)$"
        with pytest.raises(ValueError, match=message):
            distancematrix(data, mask=np.zeros(3), weight=weight)
        message = "^mask has incorrect dimensions \\(4 x 3, expected 9 x 3\\)$"
        with pytest.raises(ValueError, match=message):
            distancematrix(
                data,
                mask=mask[:4, :],
                weight=weight,
                transpose=False,
                dist="c",
                distancematrix=[],
            )
        message = "^incorrect rank 2 \\(expected 1\\)$"
        with pytest.raises(ValueError, match=message):
            distancematrix(data, mask=mask, weight=np.zeros((2, 2)))
        message = "^weight has incorrect size 4 \\(expected 3\\)$"
        with pytest.raises(ValueError, match=message):
            distancematrix(
                data,
                mask=mask,
                weight=np.zeros(4),
                transpose=False,
                dist="c",
                distancematrix=[],
            )

    def test_kmedoids_arguments(self):
        # Test if incorrect arguments are caught by the C code
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster._cluster import distancematrix
            from Bio.Cluster._cluster import kmedoids
        elif TestCluster.module == "Pycluster":
            from Pycluster._cluster import distancematrix
            from Pycluster._cluster import kmedoids

        clusterid = np.zeros(10, np.int32)
        message = "^failed to parse row 0.$"
        with pytest.raises(RuntimeError, match=message):
            kmedoids([None])
        message = "^more clusters requested than items to be clustered$"
        with pytest.raises(ValueError, match=message):
            kmedoids([], nclusters=2, npass=1000, clusterid=clusterid)
        message = "^distance matrix is not square.$"
        with pytest.raises(ValueError, match=message):
            kmedoids(np.zeros((2, 3)), npass=1000)
        message = "^distance matrix has incorrect rank 3 \\(expected 1 or 2\\)$"
        with pytest.raises(ValueError, match=message):
            kmedoids(np.zeros((2, 3, 4)), npass=1000)

    def test_distancematrix_kmedoids(self):
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster import distancematrix
            from Bio.Cluster import kmedoids
        elif TestCluster.module == "Pycluster":
            from Pycluster import distancematrix
            from Pycluster import kmedoids

        # transpose=False
        data = np.array(
            [
                [2.2, 3.3, 4.4],
                [2.1, 1.4, 5.6],
                [7.8, 9.0, 1.2],
                [4.5, 2.3, 1.5],
                [4.2, 2.4, 1.9],
                [3.6, 3.1, 9.3],
                [2.3, 1.2, 3.9],
                [4.2, 9.6, 9.3],
                [1.7, 8.9, 1.1],
            ]
        )
        mask = np.array(
            [
                [1, 1, 1],
                [1, 1, 1],
                [0, 1, 1],
                [1, 1, 1],
                [1, 1, 1],
                [0, 1, 0],
                [1, 1, 1],
                [1, 0, 1],
                [1, 1, 1],
            ],
            int,
        )
        weight = np.array([2.0, 1.0, 0.5])
        matrix = distancematrix(data, mask=mask, weight=weight)

        assert matrix[1][0] == pytest.approx(1.243, abs=0.0005)

        assert matrix[2][0] == pytest.approx(25.073, abs=0.0005)
        assert matrix[2][1] == pytest.approx(44.960, abs=0.0005)

        assert matrix[3][0] == pytest.approx(4.510, abs=0.0005)
        assert matrix[3][1] == pytest.approx(5.924, abs=0.0005)
        assert matrix[3][2] == pytest.approx(29.957, abs=0.0005)

        assert matrix[4][0] == pytest.approx(3.410, abs=0.0005)
        assert matrix[4][1] == pytest.approx(4.761, abs=0.0005)
        assert matrix[4][2] == pytest.approx(29.203, abs=0.0005)
        assert matrix[4][3] == pytest.approx(0.077, abs=0.0005)

        assert matrix[5][0] == pytest.approx(0.040, abs=0.0005)
        assert matrix[5][1] == pytest.approx(2.890, abs=0.0005)
        assert matrix[5][2] == pytest.approx(34.810, abs=0.0005)
        assert matrix[5][3] == pytest.approx(0.640, abs=0.0005)
        assert matrix[5][4] == pytest.approx(0.490, abs=0.0005)

        assert matrix[6][0] == pytest.approx(1.301, abs=0.0005)
        assert matrix[6][1] == pytest.approx(0.447, abs=0.0005)
        assert matrix[6][2] == pytest.approx(42.990, abs=0.0005)
        assert matrix[6][3] == pytest.approx(3.934, abs=0.0005)
        assert matrix[6][4] == pytest.approx(3.046, abs=0.0005)
        assert matrix[6][5] == pytest.approx(3.610, abs=0.0005)

        assert matrix[7][0] == pytest.approx(8.002, abs=0.0005)
        assert matrix[7][1] == pytest.approx(6.266, abs=0.0005)
        assert matrix[7][2] == pytest.approx(65.610, abs=0.0005)
        assert matrix[7][3] == pytest.approx(12.240, abs=0.0005)
        assert matrix[7][4] == pytest.approx(10.952, abs=0.0005)
        assert matrix[7][5] == pytest.approx(0.000, abs=0.0005)
        assert matrix[7][6] == pytest.approx(8.720, abs=0.0005)

        assert matrix[8][0] == pytest.approx(10.659, abs=0.0005)
        assert matrix[8][1] == pytest.approx(19.056, abs=0.0005)
        assert matrix[8][2] == pytest.approx(0.010, abs=0.0005)
        assert matrix[8][3] == pytest.approx(16.949, abs=0.0005)
        assert matrix[8][4] == pytest.approx(15.734, abs=0.0005)
        assert matrix[8][5] == pytest.approx(33.640, abs=0.0005)
        assert matrix[8][6] == pytest.approx(18.266, abs=0.0005)
        assert matrix[8][7] == pytest.approx(18.448, abs=0.0005)

        clusterid, error, nfound = kmedoids(matrix, npass=1000)
        assert clusterid[0] == 5
        assert clusterid[1] == 5
        assert clusterid[2] == 2
        assert clusterid[3] == 5
        assert clusterid[4] == 5
        assert clusterid[5] == 5
        assert clusterid[6] == 5
        assert clusterid[7] == 5
        assert clusterid[8] == 2
        assert error == pytest.approx(7.680, abs=0.0005)

        # check if default weights can be used
        matrix = distancematrix(data, mask=mask)
        assert len(matrix) == 9
        for i in range(3):
            assert len(matrix[i]) == i

        assert matrix[1][0] == pytest.approx(1.687, abs=0.0005)

        assert matrix[2][0] == pytest.approx(21.365, abs=0.0005)
        assert matrix[2][1] == pytest.approx(38.560, abs=0.0005)

        assert matrix[3][0] == pytest.approx(4.900, abs=0.0005)
        assert matrix[3][1] == pytest.approx(7.793, abs=0.0005)
        assert matrix[3][2] == pytest.approx(22.490, abs=0.0005)

        assert matrix[4][0] == pytest.approx(3.687, abs=0.0005)
        assert matrix[4][1] == pytest.approx(6.367, abs=0.0005)
        assert matrix[4][2] == pytest.approx(22.025, abs=0.0005)
        assert matrix[4][3] == pytest.approx(0.087, abs=0.0005)

        assert matrix[5][0] == pytest.approx(0.040, abs=0.0005)
        assert matrix[5][1] == pytest.approx(2.890, abs=0.0005)
        assert matrix[5][2] == pytest.approx(34.810, abs=0.0005)
        assert matrix[5][3] == pytest.approx(0.640, abs=0.0005)
        assert matrix[5][4] == pytest.approx(0.490, abs=0.0005)

        assert matrix[6][0] == pytest.approx(1.557, abs=0.0005)
        assert matrix[6][1] == pytest.approx(0.990, abs=0.0005)
        assert matrix[6][2] == pytest.approx(34.065, abs=0.0005)
        assert matrix[6][3] == pytest.approx(3.937, abs=0.0005)
        assert matrix[6][4] == pytest.approx(3.017, abs=0.0005)
        assert matrix[6][5] == pytest.approx(3.610, abs=0.0005)

        assert matrix[7][0] == pytest.approx(14.005, abs=0.0005)
        assert matrix[7][1] == pytest.approx(9.050, abs=0.0005)
        assert matrix[7][2] == pytest.approx(65.610, abs=0.0005)
        assert matrix[7][3] == pytest.approx(30.465, abs=0.0005)
        assert matrix[7][4] == pytest.approx(27.380, abs=0.0005)
        assert matrix[7][5] == pytest.approx(0.000, abs=0.0005)
        assert matrix[7][6] == pytest.approx(16.385, abs=0.0005)

        assert matrix[8][0] == pytest.approx(14.167, abs=0.0005)
        assert matrix[8][1] == pytest.approx(25.553, abs=0.0005)
        assert matrix[8][2] == pytest.approx(0.010, abs=0.0005)
        assert matrix[8][3] == pytest.approx(17.187, abs=0.0005)
        assert matrix[8][4] == pytest.approx(16.380, abs=0.0005)
        assert matrix[8][5] == pytest.approx(33.640, abs=0.0005)
        assert matrix[8][6] == pytest.approx(22.497, abs=0.0005)
        assert matrix[8][7] == pytest.approx(36.745, abs=0.0005)

        # transpose=True
        weight = np.array([2.0, 1.0, 0.5, 0.1, 0.9, 3.0, 2.0, 1.5, 0.2])
        matrix = distancematrix(data, mask=mask, weight=weight, transpose=True)
        assert len(matrix) == 3
        for i in range(3):
            assert len(matrix[i]) == i

        assert matrix[1][0] == pytest.approx(3.080323, abs=0.0005)
        assert matrix[2][0] == pytest.approx(9.324416, abs=0.0005)
        assert matrix[2][1] == pytest.approx(11.569701, abs=0.0005)

        clusterid, error, nfound = kmedoids(matrix, npass=1000)
        assert clusterid[0] == 0
        assert clusterid[1] == 0
        assert clusterid[2] == 2
        assert error == pytest.approx(3.08032258, abs=0.0005)

        # check if default weights can be used
        matrix = distancematrix(data, mask=mask, transpose=True)
        assert len(matrix) == 3
        for i in range(3):
            assert len(matrix[i]) == i

        assert matrix[1][0] == pytest.approx(10.47166667, abs=0.0005)
        assert matrix[2][0] == pytest.approx(8.61571429, abs=0.0005)
        assert matrix[2][1] == pytest.approx(21.24428571, abs=0.0005)

    def test_pca_arguments(self):
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster._cluster import pca
        elif TestCluster.module == "Pycluster":
            from Pycluster._cluster import pca

        data = np.zeros((4, 2))
        columnmean = np.zeros(2)
        pc = np.zeros((2, 2), dtype="d")
        coordinates = np.zeros((4, 2), dtype="d")
        eigenvalues = np.zeros(2, dtype="d")

        message = "^data matrix has unexpected format.$"
        with pytest.raises(RuntimeError, match=message):
            pca([None], columnmean, coordinates, pc, eigenvalues)
        message = "^data matrix has incorrect rank 1 \\(expected 2\\)$"
        with pytest.raises(RuntimeError, match=message):
            pca(np.zeros(3), columnmean, coordinates, pc, eigenvalues)
        message = "^data matrix has incorrect data type$"
        with pytest.raises(RuntimeError, match=message):
            pca(
                np.zeros((3, 3), dtype=np.int16),
                columnmean,
                coordinates,
                pc,
                eigenvalues,
            )
        message = "^unexpected format.$"
        with pytest.raises(RuntimeError, match=message):
            pca(data, "nothing", coordinates, pc, eigenvalues)
        message = "^incorrect rank 2 \\(expected 1\\)$"
        with pytest.raises(ValueError, match=message):
            pca(data, np.zeros((2, 2)), coordinates, pc, eigenvalues)
        message = "^array has incorrect data type$"
        with pytest.raises(RuntimeError, match=message):
            pca(data, np.ones(3, dtype=np.int16), coordinates, pc, eigenvalues)
        message = "^data matrix has unexpected format.$"
        with pytest.raises(RuntimeError, match=message):
            pca(data, columnmean, [None], pc, eigenvalues)
        message = "^data matrix has incorrect rank 1 \\(expected 2\\)$"
        with pytest.raises(RuntimeError, match=message):
            pca(data, columnmean, np.zeros(3), pc, eigenvalues)
        message = "^data matrix has incorrect data type$"
        with pytest.raises(RuntimeError, match=message):
            pca(
                data,
                columnmean,
                np.zeros((3, 3), dtype=np.int16),
                pc,
                eigenvalues,
            )
        message = "^data matrix has unexpected format.$"
        with pytest.raises(RuntimeError, match=message):
            pca(data, columnmean, coordinates, [None], eigenvalues)
        message = "^data matrix has incorrect rank 1 \\(expected 2\\)$"
        with pytest.raises(RuntimeError, match=message):
            pca(data, columnmean, coordinates, np.zeros(3), eigenvalues)
        message = "^data matrix has incorrect data type$"
        with pytest.raises(RuntimeError, match=message):
            pca(
                data,
                columnmean,
                coordinates,
                np.zeros((3, 3), dtype=np.int16),
                eigenvalues,
            )
        message = "^unexpected format.$"
        with pytest.raises(RuntimeError, match=message):
            pca(data, columnmean, coordinates, pc, "nothing")
        message = "^incorrect rank 2 \\(expected 1\\)$"
        with pytest.raises(ValueError, match=message):
            pca(data, columnmean, coordinates, pc, np.zeros((2, 2)))
        message = "^array has incorrect data type$"
        with pytest.raises(RuntimeError, match=message):
            pca(data, columnmean, coordinates, pc, np.ones(3, dtype=np.int16))

    def test_pca(self):
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster import pca
        elif TestCluster.module == "Pycluster":
            from Pycluster import pca

        data = np.array(
            [
                [3.1, 1.2],
                [1.4, 1.3],
                [1.1, 1.5],
                [2.0, 1.5],
                [1.7, 1.9],
                [1.7, 1.9],
                [5.7, 5.9],
                [5.7, 5.9],
                [3.1, 3.3],
                [5.4, 5.3],
                [5.1, 5.5],
                [5.0, 5.5],
                [5.1, 5.2],
            ]
        )

        mean, coordinates, pc, eigenvalues = pca(data)
        assert mean[0] == pytest.approx(3.5461538461538464, abs=5e-8)
        assert mean[1] == pytest.approx(3.5307692307692311, abs=5e-8)
        assert coordinates[0, 0] == pytest.approx(2.0323189722653883, abs=5e-8)
        assert coordinates[0, 1] == pytest.approx(1.2252420399694917, abs=5e-8)
        assert coordinates[1, 0] == pytest.approx(3.0936985166252251, abs=5e-8)
        assert coordinates[1, 1] == pytest.approx(-0.10647619705157851, abs=5e-8)
        assert coordinates[2, 0] == pytest.approx(3.1453186907749426, abs=5e-8)
        assert coordinates[2, 1] == pytest.approx(-0.46331699855941139, abs=5e-8)
        assert coordinates[3, 0] == pytest.approx(2.5440202962223761, abs=5e-8)
        assert coordinates[3, 1] == pytest.approx(0.20633980959571077, abs=5e-8)
        assert coordinates[4, 0] == pytest.approx(2.4468278463376221, abs=5e-8)
        assert coordinates[4, 1] == pytest.approx(-0.28412285736824866, abs=5e-8)
        assert coordinates[5, 0] == pytest.approx(2.4468278463376221, abs=5e-8)
        assert coordinates[5, 1] == pytest.approx(-0.28412285736824866, abs=5e-8)
        assert coordinates[6, 0] == pytest.approx(-3.2018619434743254, abs=5e-8)
        assert coordinates[6, 1] == pytest.approx(0.019692314198662915, abs=5e-8)
        assert coordinates[7, 0] == pytest.approx(-3.2018619434743254, abs=5e-8)
        assert coordinates[7, 1] == pytest.approx(0.019692314198662915, abs=5e-8)
        assert coordinates[8, 0] == pytest.approx(0.46978641990344067, abs=5e-8)
        assert coordinates[8, 1] == pytest.approx(-0.17778754731982949, abs=5e-8)
        assert coordinates[9, 0] == pytest.approx(-2.5549912731867215, abs=5e-8)
        assert coordinates[9, 1] == pytest.approx(0.19733897451533403, abs=5e-8)
        assert coordinates[10, 0] == pytest.approx(-2.5033710990370044, abs=5e-8)
        assert coordinates[10, 1] == pytest.approx(-0.15950182699250004, abs=5e-8)
        assert coordinates[11, 0] == pytest.approx(-2.4365601663089413, abs=5e-8)
        assert coordinates[11, 1] == pytest.approx(-0.23390813900973562, abs=5e-8)
        assert coordinates[12, 0] == pytest.approx(-2.2801521629852974, abs=5e-8)
        assert coordinates[12, 1] == pytest.approx(0.0409309711916888, abs=5e-8)
        assert pc[0, 0] == pytest.approx(-0.66810932728062988, abs=5e-8)
        assert pc[0, 1] == pytest.approx(-0.74406312017235743, abs=5e-8)
        assert pc[1, 0] == pytest.approx(0.74406312017235743, abs=5e-8)
        assert pc[1, 1] == pytest.approx(-0.66810932728062988, abs=5e-8)
        assert eigenvalues[0] == pytest.approx(9.3110471246032844, abs=5e-8)
        assert eigenvalues[1] == pytest.approx(1.4437456297481428, abs=5e-8)

        data = np.array(
            [
                [2.3, 4.5, 1.2, 6.7, 5.3, 7.1],
                [1.3, 6.5, 2.2, 5.7, 6.2, 9.1],
                [3.2, 7.2, 3.2, 7.4, 7.3, 8.9],
                [4.2, 5.2, 9.2, 4.4, 6.3, 7.2],
            ]
        )
        mean, coordinates, pc, eigenvalues = pca(data)
        assert mean[0] == pytest.approx(2.7500, abs=5e-8)
        assert mean[1] == pytest.approx(5.8500, abs=5e-8)
        assert mean[2] == pytest.approx(3.9500, abs=5e-8)
        assert mean[3] == pytest.approx(6.0500, abs=5e-8)
        assert mean[4] == pytest.approx(6.2750, abs=5e-8)
        assert mean[5] == pytest.approx(8.0750, abs=5e-8)
        assert coordinates[0, 0] == pytest.approx(2.6460846688406905, abs=5e-8)
        assert coordinates[0, 1] == pytest.approx(-2.1421701432732418, abs=5e-8)
        assert coordinates[0, 2] == pytest.approx(-0.56620932754145858, abs=5e-8)
        assert coordinates[0, 3] == pytest.approx(0.0, abs=5e-8)
        assert coordinates[1, 0] == pytest.approx(2.0644120899917544, abs=5e-8)
        assert coordinates[1, 1] == pytest.approx(0.55542108669180323, abs=5e-8)
        assert coordinates[1, 2] == pytest.approx(1.4818772348457117, abs=5e-8)
        assert coordinates[1, 3] == pytest.approx(0.0, abs=5e-8)
        assert coordinates[2, 0] == pytest.approx(1.0686641862092987, abs=5e-8)
        assert coordinates[2, 1] == pytest.approx(1.9994412069101073, abs=5e-8)
        assert coordinates[2, 2] == pytest.approx(-1.000720598980291, abs=5e-8)
        assert coordinates[2, 3] == pytest.approx(0.0, abs=5e-8)
        assert coordinates[3, 0] == pytest.approx(-5.77916094504174, abs=5e-8)
        assert coordinates[3, 1] == pytest.approx(-0.41269215032867046, abs=5e-8)
        assert coordinates[3, 2] == pytest.approx(0.085052691676038017, abs=5e-8)
        assert coordinates[3, 3] == pytest.approx(0.0, abs=5e-8)
        assert pc[0, 0] == pytest.approx(-0.26379660005997291, abs=5e-8)
        assert pc[0, 1] == pytest.approx(0.064814972617134495, abs=5e-8)
        assert pc[0, 2] == pytest.approx(-0.91763310094893846, abs=5e-8)
        assert pc[0, 3] == pytest.approx(0.26145408875373249, abs=5e-8)
        assert pc[1, 0] == pytest.approx(0.05073770520434398, abs=5e-8)
        assert pc[1, 1] == pytest.approx(0.68616983388698793, abs=5e-8)
        assert pc[1, 2] == pytest.approx(0.13819106187213354, abs=5e-8)
        assert pc[1, 3] == pytest.approx(0.19782544121828985, abs=5e-8)
        assert pc[2, 0] == pytest.approx(-0.63000893660095947, abs=5e-8)
        assert pc[2, 1] == pytest.approx(0.091155993862151397, abs=5e-8)
        assert pc[2, 2] == pytest.approx(0.045630391256086845, abs=5e-8)
        assert pc[2, 3] == pytest.approx(-0.67456694780914772, abs=5e-8)
        # As the last eigenvalue is zero, the corresponding eigenvector is
        # strongly affected by roundoff error, and is not being tested here.
        # For PCA, this doesn't matter since all data have a zero coefficient
        # along this eigenvector.
        assert eigenvalues[0] == pytest.approx(6.7678878332578778, abs=5e-8)
        assert eigenvalues[1] == pytest.approx(3.0108911400291856, abs=5e-8)
        assert eigenvalues[2] == pytest.approx(1.8775592718563467, abs=5e-8)
        assert eigenvalues[3] == pytest.approx(0.0, abs=5e-8)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
