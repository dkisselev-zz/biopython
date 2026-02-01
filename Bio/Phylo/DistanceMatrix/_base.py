# Copyright 2024 by Biopython contributors. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Base classes and utilities for distance matrix computation.

This module provides the Protocol for distance calculation strategies
and the registry for registering new strategies.
"""

from collections.abc import Callable
from collections.abc import Sequence
from typing import Protocol

import numpy as np
from numpy.typing import NDArray


class DistanceStrategy(Protocol):
    """Protocol for distance calculation strategies.

    Implementations must provide a compute() method that accepts an alignment
    array and optional parameters, returning a condensed distance vector.

    The condensed form is a 1D array of length n*(n-1)/2 containing the
    lower-triangular elements of the distance matrix in row-major order.
    """

    def compute(
        self,
        alignment: NDArray[np.bytes_],
        scoring_matrix: NDArray[np.float64] | None = None,
        alphabet: str | None = None,
        skip_letters: Sequence[str] = ("-", "*"),
        n_jobs: int = 1,
    ) -> NDArray[np.float64]:
        """Compute pairwise distances for all sequence pairs.

        :param alignment: 2D array of shape (n_sequences, alignment_length),
            dtype S1 (single-byte characters).
        :param scoring_matrix: Optional 2D scoring matrix for substitution models.
            If None, identity model is used.
        :param alphabet: Alphabet string mapping indices to characters in the
            scoring matrix. Required if scoring_matrix is provided.
        :param skip_letters: Characters to skip in distance calculation (gaps).
        :param n_jobs: Number of parallel workers (1 = serial).
        :returns: Condensed distance matrix (1D array, length n*(n-1)/2).
        """
        ...


# Strategy registry
_STRATEGY_REGISTRY: dict[str, type[DistanceStrategy]] = {}


def register_strategy(name: str) -> Callable[[type], type]:
    """Decorator to register a distance strategy.

    :param name: The name to register the strategy under.
    :returns: Decorator function.

    Example
    -------
    >>> @register_strategy("custom")
    ... class CustomStrategy:
    ...     def compute(self, alignment, **kwargs):
    ...         ...

    """

    def decorator(cls: type) -> type:
        _STRATEGY_REGISTRY[name] = cls
        return cls

    return decorator


def get_strategy(name: str) -> DistanceStrategy:
    """Get a strategy instance by name.

    :param name: The registered name of the strategy.
    :returns: An instance of the strategy.
    :raises ValueError: If the strategy name is not registered.
    """
    if name not in _STRATEGY_REGISTRY:
        available = ", ".join(sorted(_STRATEGY_REGISTRY.keys()))
        raise ValueError(f"Unknown strategy '{name}'. Available: {available}")
    return _STRATEGY_REGISTRY[name]()


def list_strategies() -> list:
    """Return a list of all registered strategy names.

    :returns: List of registered strategy names.
    """
    return sorted(_STRATEGY_REGISTRY.keys())


def create_gap_mask(
    alignment: NDArray[np.bytes_], skip_letters: Sequence[str]
) -> NDArray[np.bool_]:
    """Create boolean mask where True indicates a gap/skip character.

    :param alignment: 2D array of shape (n_sequences, alignment_length), dtype S1.
    :param skip_letters: Characters to skip (e.g., ['-', '*']).
    :returns: Boolean mask of same shape as alignment.
    """
    if not skip_letters:
        return np.zeros(alignment.shape, dtype=bool)
    # Convert skip_letters to bytes array for vectorized set-membership check
    skip_bytes = np.array([s.encode() for s in skip_letters], dtype="S1")
    return np.isin(alignment, skip_bytes)


def condensed_to_squareform(
    condensed: NDArray[np.float64], n: int
) -> NDArray[np.float64]:
    """Convert condensed distance vector to square distance matrix.

    :param condensed: 1D condensed distance array of length n*(n-1)/2.
    :param n: Number of sequences (matrix dimension).
    :returns: 2D symmetric distance matrix of shape (n, n).
    """
    matrix = np.zeros((n, n), dtype=np.float64)
    # Get lower triangular indices (excluding diagonal)
    row_idx, col_idx = np.tril_indices(n, k=-1)
    # Fill lower triangle
    matrix[row_idx, col_idx] = condensed
    # Fill upper triangle (symmetric)
    matrix[col_idx, row_idx] = condensed
    return matrix
