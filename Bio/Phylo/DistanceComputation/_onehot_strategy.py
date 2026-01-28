# Copyright 2024 The Biopython Contributors.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""One-hot encoding distance computation strategy (PRIVATE).

Each sequence is encoded as a binary ``(n_cols x alphabet_size)`` matrix.
Pairwise scoring is performed as a matrix product with the substitution
matrix, fully vectorising the inner loop.  Uses ``ThreadPoolExecutor``
because NumPy matrix multiplication releases the GIL.
"""

try:
    import numpy as np
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Please install NumPy if you want to use method='onehot' in "
        "DistanceCalculator.  See https://www.numpy.org/"
    ) from None

from typing import Any

from Bio.Phylo.DistanceComputation._base import _DistanceStrategy
from Bio.Phylo.DistanceComputation._base import _validate_distance_matrix
from Bio.Phylo.DistanceComputation._parallel import _ParallelMixin


class OneHotDistanceStrategy(_DistanceStrategy, _ParallelMixin):
    """One-hot matrix-multiply distance strategy (PRIVATE).

    Identity model: the substitution matrix is ``numpy.eye`` over the
    observed alphabet, so scoring reduces to counting matches.

    Scoring model: the substitution matrix is mapped onto the observed
    alphabet indices and applied via ``onehot[i] @ score_mat``.

    Diagonal normalisation mirrors the legacy behaviour: each sequence's
    self-score is computed from its own one-hot encoding and the diagonal
    of the substitution matrix.
    """

    def compute(
        self,
        sequences: list[str],
        names: list[str],
        scoring_matrix: Any,
        skip_letters: tuple[str, ...],
    ) -> list[list[float]]:
        """Compute distances via one-hot encoding and matrix multiply (PRIVATE).

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
        skip_set = set(skip_letters) if skip_letters else set()

        # Determine observed alphabet (excluding skip characters)
        all_chars: set = set()
        for s in sequences:
            all_chars.update(s)
        all_chars -= skip_set
        alphabet = sorted(all_chars)
        char_to_idx = {c: i for i, c in enumerate(alphabet)}
        alpha_size = len(alphabet)

        # Build alpha_size x alpha_size scoring matrix
        if scoring_matrix is None:
            score_mat = np.eye(alpha_size, dtype=np.float64)
        else:
            score_mat = np.zeros((alpha_size, alpha_size), dtype=np.float64)
            for i_a, a in enumerate(alphabet):
                for i_b, b in enumerate(alphabet):
                    try:
                        score_mat[i_a, i_b] = float(scoring_matrix[a, b])
                    except (IndexError, KeyError):
                        score_mat[i_a, i_b] = 0.0

        # One-hot encode: (n_seqs, n_cols, alpha_size)
        onehot = np.zeros((n_seqs, n_cols, alpha_size), dtype=np.float64)
        for s_idx, seq in enumerate(sequences):
            for col, ch in enumerate(seq):
                if ch in skip_set:
                    continue
                if ch in char_to_idx:
                    onehot[s_idx, col, char_to_idx[ch]] = 1.0

        # Precompute per-sequence diagonal scores for normalisation
        diag_vec = np.array(
            [score_mat[i, i] for i in range(alpha_size)], dtype=np.float64
        )
        # diag_scores[s] = sum over columns of (onehot[s] @ diag_vec)
        diag_scores = np.array(
            [float(np.sum(onehot[s] @ diag_vec)) for s in range(n_seqs)],
            dtype=np.float64,
        )

        def pair_dist(i: int, j: int) -> float:
            scored_i = onehot[i] @ score_mat  # (n_cols, alpha_size)
            score = float(np.sum(scored_i * onehot[j]))
            max_score = max(float(diag_scores[i]), float(diag_scores[j]))
            if max_score == 0:
                return 1.0
            return 1.0 - (score / max_score)

        pairs = [(i, j) for i in range(n_seqs) for j in range(i)]
        results = self._compute_pairs(pair_dist, pairs, use_threads=True)

        matrix: list[list[float]] = [[0.0] * (row + 1) for row in range(n_seqs)]
        for i, j, dist in results:
            matrix[i][j] = dist

        if __debug__:
            _validate_distance_matrix(matrix)

        return matrix
