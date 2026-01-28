# Copyright 2024 The Biopython Contributors.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""SciPy-accelerated distance computation strategy (PRIVATE).

Uses ``scipy.spatial.distance.pdist`` for the identity model.  Scoring
matrix models delegate silently to ``NumpyDistanceStrategy`` because
SciPy has no built-in facility for weighted substitution-matrix scoring.
"""

try:
    import numpy as np
    from scipy.spatial.distance import pdist
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Please install SciPy if you want to use method='scipy' in "
        "DistanceCalculator.  See https://docs.scipy.org/"
    ) from None

from typing import Any

from Bio.Phylo.DistanceComputation._base import _DistanceStrategy
from Bio.Phylo.DistanceComputation._base import _validate_distance_matrix
from Bio.Phylo.DistanceComputation._numpy_strategy import NumpyDistanceStrategy


class ScipyDistanceStrategy(_DistanceStrategy):
    """Distance computation leveraging scipy.spatial.distance (PRIVATE).

    Identity model: encodes the alignment as an ``(n_seqs, n_cols)`` uint8
    matrix and calls ``pdist`` with a custom normalised-Hamming metric that
    respects ``skip_letters``.

    Scoring matrix models: delegates to ``NumpyDistanceStrategy`` because
    ``scipy.spatial.distance`` has no support for asymmetric or weighted
    character-level substitution matrices.  No warning is emitted; the
    NumPy path still delivers vectorised performance.
    """

    def __init__(self, n_jobs=None):
        """Initialise with optional worker count.

        Arguments:
         - n_jobs - Passed through to the numpy fallback strategy.
        """
        super().__init__(n_jobs=n_jobs)
        self._fallback = NumpyDistanceStrategy(n_jobs=n_jobs)

    def compute(
        self,
        sequences: list[str],
        names: list[str],
        scoring_matrix: Any,
        skip_letters: tuple[str, ...],
    ) -> list[list[float]]:
        """Compute distances, using pdist where applicable (PRIVATE).

        Arguments:
         - sequences - Aligned sequence strings (all same length).
         - names - Sequence identifiers (unused internally).
         - scoring_matrix - Substitution matrix or None for identity.
         - skip_letters - Characters to exclude from scoring.

        Returns lower-triangular list of lists.
        """
        if __debug__:
            if len(sequences) != len(names):
                raise ValueError(
                    f"len(sequences)={len(sequences)} != len(names)={len(names)}"
                )
            if sequences and len({len(s) for s in sequences}) != 1:
                raise ValueError("All sequences must have the same length")

        if scoring_matrix is not None:
            return self._fallback.compute(
                sequences, names, scoring_matrix, skip_letters
            )

        n_seqs = len(sequences)
        n_cols = len(sequences[0])

        seq_matrix = np.array(
            [np.frombuffer(s.encode("ascii"), dtype=np.uint8) for s in sequences],
            dtype=np.uint8,
        )

        skip_codes = (
            np.array([ord(c) for c in skip_letters], dtype=np.uint8)
            if skip_letters
            else np.array([], dtype=np.uint8)
        )

        def _hamming_norm(u, v):
            """Normalised identity distance respecting skip_letters (PRIVATE)."""
            if len(skip_codes) > 0:
                valid = ~np.isin(u, skip_codes) & ~np.isin(v, skip_codes)
                matches = int(np.sum(u[valid] == v[valid]))
            else:
                matches = int(np.sum(u == v))
            # Denominator is always n_cols, matching legacy identity model.
            return 1.0 - (matches / n_cols)

        condensed = pdist(seq_matrix, metric=_hamming_norm)

        # Convert condensed vector to lower-triangular list of lists.
        # pdist condensed index for pair (a, b) where a < b is:
        #   n*a - a*(a+1)//2 + b - a - 1
        # Our lower-triangular matrix[i][j] stores distance for i > j,
        # so the condensed pair is (j, i).
        matrix: list[list[float]] = [[0.0] * (row + 1) for row in range(n_seqs)]
        for i in range(1, n_seqs):
            for j in range(i):
                idx = n_seqs * j - j * (j + 1) // 2 + i - j - 1
                matrix[i][j] = float(condensed[idx])

        if __debug__:
            _validate_distance_matrix(matrix)

        return matrix
