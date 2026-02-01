# Copyright 2024 by Biopython contributors. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Distance calculation strategy implementations.

This module provides multiple implementations for computing pairwise
distance matrices from sequence alignments.
"""

from collections.abc import Sequence
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from numpy.typing import NDArray

from Bio.Phylo.DistanceMatrix._base import create_gap_mask
from Bio.Phylo.DistanceMatrix._base import register_strategy


# Module-level function for multiprocessing (must be picklable)
def _pairwise_identity(
    seq1: NDArray[np.bytes_],
    seq2: NDArray[np.bytes_],
    skip_bytes: tuple[bytes, ...],
) -> float:
    """Calculate identity distance between two sequences."""
    matches = 0
    valid = 0
    for c1, c2 in zip(seq1, seq2):
        if c1 in skip_bytes or c2 in skip_bytes:
            continue
        valid += 1
        if c1 == c2:
            matches += 1
    if valid == 0:
        return 1.0
    return 1.0 - (matches / valid)


def _pairwise_substitution(
    seq1: NDArray[np.bytes_],
    seq2: NDArray[np.bytes_],
    skip_bytes: tuple[bytes, ...],
    scoring_matrix: NDArray[np.float64],
    char_to_idx: dict,
) -> float:
    """Calculate substitution matrix distance between two sequences."""
    score = 0.0
    max_score1 = 0.0
    max_score2 = 0.0
    for c1, c2 in zip(seq1, seq2):
        if c1 in skip_bytes or c2 in skip_bytes:
            continue
        try:
            idx1 = char_to_idx[c1]
            idx2 = char_to_idx[c2]
        except KeyError:
            continue
        score += scoring_matrix[idx1, idx2]
        max_score1 += scoring_matrix[idx1, idx1]
        max_score2 += scoring_matrix[idx2, idx2]
    max_score = max(max_score1, max_score2)
    if max_score == 0:
        return 1.0
    return 1.0 - (score / max_score)


@register_strategy("python")
class PythonStrategy:
    """Pure Python implementation (legacy compatibility).

    Uses nested loops with optional multiprocessing for parallelism.
    This strategy serves as the baseline for correctness verification.
    """

    def compute(
        self,
        alignment: NDArray[np.bytes_],
        scoring_matrix: NDArray[np.float64] | None = None,
        alphabet: str | None = None,
        skip_letters: Sequence[str] = ("-", "*"),
        n_jobs: int = 1,
    ) -> NDArray[np.float64]:
        """Compute pairwise distances using pure Python loops."""
        n_seqs = alignment.shape[0]
        n_pairs = n_seqs * (n_seqs - 1) // 2
        skip_bytes = tuple(s.encode() for s in skip_letters)

        if scoring_matrix is not None and alphabet is not None:
            char_to_idx = {c.encode(): i for i, c in enumerate(alphabet)}

        distances = np.empty(n_pairs, dtype=np.float64)
        idx = 0

        if n_jobs == 1:
            # Serial execution
            for i in range(n_seqs):
                for j in range(i):
                    if scoring_matrix is None:
                        distances[idx] = _pairwise_identity(
                            alignment[i], alignment[j], skip_bytes
                        )
                    else:
                        distances[idx] = _pairwise_substitution(
                            alignment[i],
                            alignment[j],
                            skip_bytes,
                            scoring_matrix,
                            char_to_idx,
                        )
                    idx += 1
        else:
            # Parallel execution with ProcessPoolExecutor
            pairs = [(i, j) for i in range(n_seqs) for j in range(i)]
            chunk_size = max(1, len(pairs) // (n_jobs * 4))

            if scoring_matrix is None:
                identity_args = [
                    (alignment[i], alignment[j], skip_bytes) for i, j in pairs
                ]
                with ProcessPoolExecutor(max_workers=n_jobs) as executor:
                    results = list(
                        executor.map(
                            lambda x: _pairwise_identity(*x),
                            identity_args,
                            chunksize=chunk_size,
                        )
                    )
            else:
                subst_args = [
                    (
                        alignment[i],
                        alignment[j],
                        skip_bytes,
                        scoring_matrix,
                        char_to_idx,
                    )
                    for i, j in pairs
                ]
                with ProcessPoolExecutor(max_workers=n_jobs) as executor:
                    results = list(
                        executor.map(
                            lambda x: _pairwise_substitution(*x),
                            subst_args,
                            chunksize=chunk_size,
                        )
                    )
            distances = np.array(results, dtype=np.float64)

        return distances


@register_strategy("numpy")
class NumPyStrategy:
    """NumPy vectorized implementation using 2D broadcasting.

    Computes all pairwise distances using array operations.
    Threading controlled via BLAS environment variables.
    """

    def compute(
        self,
        alignment: NDArray[np.bytes_],
        scoring_matrix: NDArray[np.float64] | None = None,
        alphabet: str | None = None,
        skip_letters: Sequence[str] = ("-", "*"),
        n_jobs: int = 1,
    ) -> NDArray[np.float64]:
        """Compute pairwise distances using NumPy vectorization."""
        # Create gap mask
        gap_mask = create_gap_mask(alignment, skip_letters)

        if scoring_matrix is None:
            return self._identity_vectorized(alignment, gap_mask)
        else:
            assert alphabet is not None, "alphabet required for substitution model"
            return self._substitution_vectorized(
                alignment, gap_mask, scoring_matrix, alphabet
            )

    def _identity_vectorized(
        self, alignment: NDArray[np.bytes_], gap_mask: NDArray[np.bool_]
    ) -> NDArray[np.float64]:
        """Vectorized identity distance computation."""
        n_seqs, align_len = alignment.shape

        # Expand to 3D for pairwise comparison: (n, 1, L) vs (1, n, L)
        seq_i = alignment[:, np.newaxis, :]  # (n, 1, L)
        seq_j = alignment[np.newaxis, :, :]  # (1, n, L)

        # Masks for gaps
        gap_i = gap_mask[:, np.newaxis, :]  # (n, 1, L)
        gap_j = gap_mask[np.newaxis, :, :]  # (1, n, L)

        # Valid positions: neither is a gap
        valid = ~(gap_i | gap_j)  # (n, n, L)

        # Matches: same character at valid positions
        matches = (seq_i == seq_j) & valid  # (n, n, L)

        # Sum along alignment axis
        n_matches = matches.sum(axis=2).astype(np.float64)  # (n, n)
        n_valid = valid.sum(axis=2).astype(np.float64)  # (n, n)

        # Distance = 1 - (matches / valid), handle division by zero
        with np.errstate(divide="ignore", invalid="ignore"):
            distances = 1.0 - (n_matches / n_valid)
            distances = np.where(n_valid == 0, 1.0, distances)

        # Extract lower triangle as condensed form
        return distances[np.tril_indices(n_seqs, k=-1)]

    def _substitution_vectorized(
        self,
        alignment: NDArray[np.bytes_],
        gap_mask: NDArray[np.bool_],
        scoring_matrix: NDArray[np.float64],
        alphabet: str,
    ) -> NDArray[np.float64]:
        """Vectorized substitution matrix distance computation."""
        n_seqs, align_len = alignment.shape

        # Create character to index mapping
        char_to_idx = {c.encode(): i for i, c in enumerate(alphabet)}
        unknown_idx = -1

        # Convert alignment to integer indices
        idx_alignment = np.zeros(alignment.shape, dtype=np.int32)
        for i in range(n_seqs):
            for j in range(align_len):
                char = alignment[i, j]
                idx_alignment[i, j] = char_to_idx.get(char, unknown_idx)

        # Mark unknown characters as gaps
        unknown_mask = idx_alignment == unknown_idx
        gap_mask = gap_mask | unknown_mask

        # Get scores using advanced indexing
        # For each position, we need scores[i, j] for all pairs
        idx_i = idx_alignment[:, np.newaxis, :]  # (n, 1, L)
        idx_j = idx_alignment[np.newaxis, :, :]  # (1, n, L)

        # Clip indices to valid range for scoring matrix lookup
        idx_i_clipped = np.clip(idx_i, 0, scoring_matrix.shape[0] - 1)
        idx_j_clipped = np.clip(idx_j, 0, scoring_matrix.shape[1] - 1)

        # Get pairwise scores at each position
        scores = scoring_matrix[idx_i_clipped, idx_j_clipped]  # (n, n, L)
        max_scores_i = scoring_matrix[idx_i_clipped, idx_i_clipped]  # (n, n, L)
        max_scores_j = scoring_matrix[idx_j_clipped, idx_j_clipped]  # (n, n, L)

        # Apply gap mask
        gap_i = gap_mask[:, np.newaxis, :]
        gap_j = gap_mask[np.newaxis, :, :]
        valid = ~(gap_i | gap_j)

        scores = np.where(valid, scores, 0.0)
        max_scores_i = np.where(valid, max_scores_i, 0.0)
        max_scores_j = np.where(valid, max_scores_j, 0.0)

        # Sum along alignment axis
        total_score = scores.sum(axis=2)
        total_max = np.maximum(max_scores_i.sum(axis=2), max_scores_j.sum(axis=2))

        # Compute distance
        with np.errstate(divide="ignore", invalid="ignore"):
            distances = 1.0 - (total_score / total_max)
            distances = np.where(total_max == 0, 1.0, distances)

        return distances[np.tril_indices(n_seqs, k=-1)]


@register_strategy("scipy")
class SciPyStrategy:
    """SciPy implementation using scipy.spatial.distance.pdist.

    Provides optimized distance computation with custom metric support.
    """

    def compute(
        self,
        alignment: NDArray[np.bytes_],
        scoring_matrix: NDArray[np.float64] | None = None,
        alphabet: str | None = None,
        skip_letters: Sequence[str] = ("-", "*"),
        n_jobs: int = 1,
    ) -> NDArray[np.float64]:
        """Compute pairwise distances using scipy.spatial.distance.pdist."""
        try:
            from scipy.spatial.distance import pdist
            from scipy.spatial.distance import squareform
        except ImportError:
            raise ImportError(
                "SciPy is required for the 'scipy' strategy. "
                "Install with: pip install scipy"
            )

        skip_bytes = tuple(s.encode() for s in skip_letters)
        n_seqs = alignment.shape[0]

        if scoring_matrix is None:
            # Use custom identity metric with gap handling
            def identity_metric(u, v):
                return _pairwise_identity(u, v, skip_bytes)

            condensed_upper = pdist(alignment, metric=identity_metric)
        else:
            assert alphabet is not None, "alphabet required for substitution model"
            char_to_idx = {c.encode(): i for i, c in enumerate(alphabet)}

            def substitution_metric(u, v):
                return _pairwise_substitution(
                    u, v, skip_bytes, scoring_matrix, char_to_idx
                )

            condensed_upper = pdist(alignment, metric=substitution_metric)

        # pdist returns upper-triangular order (i < j), but we need
        # lower-triangular order (i > j) to match other strategies.
        # Convert via squareform and extract lower triangle.
        square = squareform(condensed_upper)
        return square[np.tril_indices(n_seqs, k=-1)]


@register_strategy("onehot")
class OneHotStrategy:
    """One-hot encoding with matrix multiplication.

    Converts sequences to one-hot encoding and uses matrix operations
    for fast distance computation. Particularly efficient for identity
    distance with many sequences.
    """

    def compute(
        self,
        alignment: NDArray[np.bytes_],
        scoring_matrix: NDArray[np.float64] | None = None,
        alphabet: str | None = None,
        skip_letters: Sequence[str] = ("-", "*"),
        n_jobs: int = 1,
    ) -> NDArray[np.float64]:
        """Compute pairwise distances using one-hot encoding."""
        n_seqs, align_len = alignment.shape

        # Determine alphabet if not provided
        if alphabet is None:
            unique_chars = set(alignment.flatten())
            skip_bytes = {s.encode() for s in skip_letters}
            alphabet = "".join(
                sorted(c.decode() for c in unique_chars if c not in skip_bytes)
            )

        # Create one-hot encoding
        onehot, valid_mask = self._to_onehot(alignment, alphabet, skip_letters)

        if scoring_matrix is None:
            return self._identity_onehot(onehot, valid_mask)
        else:
            return self._substitution_onehot(onehot, valid_mask, scoring_matrix)

    def _to_onehot(
        self,
        alignment: NDArray[np.bytes_],
        alphabet: str,
        skip_letters: Sequence[str],
    ) -> tuple[NDArray[np.float32], NDArray[np.bool_]]:
        """Convert alignment to one-hot encoding.

        Returns (onehot, valid_mask) where valid_mask indicates non-gap positions.
        """
        n_seqs, align_len = alignment.shape
        n_chars = len(alphabet)

        char_to_idx = {c.encode(): i for i, c in enumerate(alphabet)}
        skip_bytes = {s.encode() for s in skip_letters}

        onehot = np.zeros((n_seqs, align_len, n_chars), dtype=np.float32)
        valid_mask = np.ones((n_seqs, align_len), dtype=bool)

        for i in range(n_seqs):
            for j in range(align_len):
                char = alignment[i, j]
                if char in skip_bytes:
                    valid_mask[i, j] = False
                elif char in char_to_idx:
                    onehot[i, j, char_to_idx[char]] = 1.0
                else:
                    valid_mask[i, j] = False

        return onehot, valid_mask

    def _identity_onehot(
        self,
        onehot: NDArray[np.float32],
        valid_mask: NDArray[np.bool_],
    ) -> NDArray[np.float64]:
        """Identity distance using one-hot dot products."""
        n_seqs, align_len, n_chars = onehot.shape

        # Flatten to (n_seqs, align_len * n_chars)
        flat = onehot.reshape(n_seqs, -1)

        # Dot product gives number of matches
        similarity = flat @ flat.T  # (n_seqs, n_seqs)

        # Count valid positions for each pair
        valid_i = valid_mask[:, np.newaxis, :]  # (n, 1, L)
        valid_j = valid_mask[np.newaxis, :, :]  # (1, n, L)
        valid_both = (valid_i & valid_j).sum(axis=2).astype(np.float64)  # (n, n)

        with np.errstate(divide="ignore", invalid="ignore"):
            distances = 1.0 - (similarity / valid_both)
            distances = np.where(valid_both == 0, 1.0, distances)

        return distances[np.tril_indices(n_seqs, k=-1)]

    def _substitution_onehot(
        self,
        onehot: NDArray[np.float32],
        valid_mask: NDArray[np.bool_],
        scoring_matrix: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        """Substitution distance using one-hot and scoring matrix."""
        n_seqs, align_len, n_chars = onehot.shape

        # Ensure scoring matrix is float64
        scoring_matrix = scoring_matrix.astype(np.float64)

        # Compute valid mask for each pair: positions where BOTH sequences are valid
        # valid_mask: (n, L) -> valid_both: (n, n, L)
        valid_i = valid_mask[:, np.newaxis, :]  # (n, 1, L)
        valid_j = valid_mask[np.newaxis, :, :]  # (1, n, L)
        valid_both = valid_i & valid_j  # (n, n, L)

        # Transform onehot by scoring matrix: onehot @ scoring_matrix
        # Shape: (n, L, n_chars) @ (n_chars, n_chars) -> (n, L, n_chars)
        transformed = np.einsum("ijk,kl->ijl", onehot, scoring_matrix)

        # Compute per-position scores for all pairs
        # score_per_pos[i, j, k] = dot(transformed[i, k, :], onehot[j, k, :])
        score_per_pos = np.einsum("ikc,jkc->ijk", transformed, onehot)  # (n, n, L)

        # Mask invalid positions and sum
        score_per_pos = np.where(valid_both, score_per_pos, 0.0)
        scores = score_per_pos.sum(axis=2)  # (n, n)

        # Max scores per position: diagonal element for each character
        diag = np.diag(scoring_matrix)
        # max_score_per_pos[seq, pos] = diag[char_idx] for the character at that position
        max_score_per_pos = np.einsum("ijk,k->ij", onehot, diag)  # (n, L)

        # For each pair (i, j), compute max_score_i and max_score_j
        # considering only positions where BOTH are valid
        max_i_per_pos = max_score_per_pos[:, np.newaxis, :]  # (n, 1, L)
        max_j_per_pos = max_score_per_pos[np.newaxis, :, :]  # (1, n, L)

        # Mask and sum for each pair
        max_i_masked = np.where(valid_both, max_i_per_pos, 0.0)
        max_j_masked = np.where(valid_both, max_j_per_pos, 0.0)
        max_scores_i = max_i_masked.sum(axis=2)  # (n, n)
        max_scores_j = max_j_masked.sum(axis=2)  # (n, n)

        # Take maximum of the two max_scores for each pair
        max_score_pairs = np.maximum(max_scores_i, max_scores_j)  # (n, n)

        # Compute distances
        with np.errstate(divide="ignore", invalid="ignore"):
            distances = 1.0 - (scores / max_score_pairs)
            distances = np.where(max_score_pairs == 0, 1.0, distances)

        return distances[np.tril_indices(n_seqs, k=-1)]
