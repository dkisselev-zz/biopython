# Copyright 2024 The Biopython Contributors.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""NumPy-accelerated distance computation strategy (PRIVATE).

Encodes sequences as integer arrays and leverages NumPy broadcasting
and a precomputed 256x256 scoring lookup table to vectorise the inner
distance loop.  Uses ``ThreadPoolExecutor`` for parallelism because
NumPy array operations release the GIL.
"""

try:
    import numpy as np
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Please install NumPy if you want to use method='numpy' in "
        "DistanceCalculator.  See https://www.numpy.org/"
    ) from None

from typing import Any

from Bio.Phylo.DistanceComputation._base import _DistanceStrategy
from Bio.Phylo.DistanceComputation._base import _validate_distance_matrix
from Bio.Phylo.DistanceComputation._parallel import _ParallelMixin


class NumpyDistanceStrategy(_DistanceStrategy, _ParallelMixin):
    """Vectorised pairwise distance via NumPy array operations (PRIVATE).

    Identity model: sequences are compared as uint8 arrays element-wise;
    the denominator is always the full alignment length (matching legacy
    behaviour where ``skip_letters`` defaults to the empty tuple).

    Scoring model: a 256x256 lookup table indexed by character ordinal
    is built once from the substitution matrix.  Per pair, scores and
    diagonal maxima are summed using fancy indexing over valid (non-skip)
    positions.

    Threading is used for parallelism because NumPy releases the GIL
    during array computations.
    """

    def compute(
        self,
        sequences: list[str],
        names: list[str],
        scoring_matrix: Any,
        skip_letters: tuple[str, ...],
    ) -> list[list[float]]:
        """Compute distances using NumPy vectorisation (PRIVATE).

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

        n_seqs = len(sequences)
        n_cols = len(sequences[0])

        # Encode each sequence as a uint8 array of character ordinals
        seq_arrays = [
            np.frombuffer(s.encode("ascii"), dtype=np.uint8) for s in sequences
        ]

        skip_codes = list({ord(c) for c in skip_letters}) if skip_letters else []

        if scoring_matrix is None:
            # --- Identity model ---
            # Denominator is ALWAYS n_cols (alignment length), matching legacy.
            def pair_dist(i: int, j: int) -> float:
                arr_i, arr_j = seq_arrays[i], seq_arrays[j]
                if skip_codes:
                    valid = ~np.isin(arr_i, skip_codes) & ~np.isin(arr_j, skip_codes)
                    matches = int(np.sum(arr_i[valid] == arr_j[valid]))
                else:
                    matches = int(np.sum(arr_i == arr_j))
                return 1.0 - (matches / n_cols)

        else:
            # --- Scoring model ---
            # Build 256x256 score table and 256-element diagonal table.
            table = np.zeros((256, 256), dtype=np.float64)
            diag = np.zeros(256, dtype=np.float64)
            for a in scoring_matrix.alphabet:
                code_a = ord(a)
                diag[code_a] = float(scoring_matrix[a, a])
                for b in scoring_matrix.alphabet:
                    table[code_a, ord(b)] = float(scoring_matrix[a, b])

            def pair_dist(i: int, j: int) -> float:
                arr_i, arr_j = seq_arrays[i], seq_arrays[j]
                if skip_codes:
                    valid = ~np.isin(arr_i, skip_codes) & ~np.isin(arr_j, skip_codes)
                    arr_i_v, arr_j_v = arr_i[valid], arr_j[valid]
                else:
                    arr_i_v, arr_j_v = arr_i, arr_j

                score = float(np.sum(table[arr_i_v, arr_j_v]))
                max_s1 = float(np.sum(diag[arr_i_v]))
                max_s2 = float(np.sum(diag[arr_j_v]))
                max_score = max(max_s1, max_s2)
                if max_score == 0:
                    return 1.0
                return 1.0 - (score / max_score)

        # Generate all lower-triangular pairs
        pairs = [(i, j) for i in range(n_seqs) for j in range(i)]
        results = self._compute_pairs(pair_dist, pairs, use_threads=True)

        # Assemble into lower-triangular list of lists
        matrix: list[list[float]] = [[0.0] * (row + 1) for row in range(n_seqs)]
        for i, j, dist in results:
            matrix[i][j] = dist

        if __debug__:
            _validate_distance_matrix(matrix)

        return matrix
