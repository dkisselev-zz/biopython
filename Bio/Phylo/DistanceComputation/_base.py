# Copyright 2024 The Biopython Contributors.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Abstract base class for distance computation strategies (PRIVATE)."""

import abc
from typing import Any


def _validate_distance_matrix(matrix: list[list[float]]) -> None:
    """Verify lower-triangular shape and zero diagonals (PRIVATE).

    Arguments:
     - matrix - Lower-triangular list of lists as returned by
       ``_DistanceStrategy.compute``.

    Raises ValueError if any row has the wrong length or any diagonal
    entry is not exactly ``0.0``.  This function is intended for use in
    unit tests and debug builds (``python -O`` disables it via the
    ``__debug__`` guard in ``_DistanceStrategy.compute``).
    """
    for i, row in enumerate(matrix):
        expected_len = i + 1
        if len(row) != expected_len:
            raise ValueError(
                f"Row {i} has {len(row)} entries; expected {expected_len} "
                f"(lower-triangular format requires row i to have i + 1 entries)"
            )
        if row[i] != 0.0:
            raise ValueError(
                f"Diagonal entry matrix[{i}][{i}] is {row[i]}; must be 0.0"
            )


class _DistanceStrategy(abc.ABC):
    """Base class for all distance computation backends (PRIVATE).

    Subclasses implement ``compute``, which accepts pre-extracted sequence
    strings and scoring configuration, and returns a lower-triangular
    distance matrix as a list of lists.  The caller (``DistanceCalculator``)
    wraps the result into a ``DistanceMatrix`` object.

    Contract requirements enforced on the inputs to ``compute``:

    - ``len(sequences) == len(names)`` — every sequence must have a
      corresponding identifier.
    - All strings in ``sequences`` must have the same length (i.e. the
      caller must supply a gapless-or-gapped *aligned* MSA).

    These preconditions are validated automatically when Python is run
    without the ``-O`` flag (i.e. while ``__debug__`` is ``True``).
    """

    def __init__(self, n_jobs: int | None = None) -> None:
        """Initialise the strategy.

        Arguments:
         - n_jobs - Number of parallel workers.  None triggers
           auto-detection in ``_ParallelMixin``.
        """
        self._n_jobs = n_jobs

    @abc.abstractmethod
    def compute(
        self,
        sequences: list[str],
        names: list[str],
        scoring_matrix: Any,
        skip_letters: tuple[str, ...],
    ) -> list[list[float]]:
        """Compute pairwise distances and return lower-triangular matrix (PRIVATE).

        Arguments:
         - sequences - Aligned sequence strings.  All entries must have
           the same length (columns).  ``len(sequences)`` must equal
           ``len(names)``.
         - names - Corresponding sequence identifiers.  Must be the same
           length as ``sequences``; order must match.  The values are not
           used by the computation itself but are required so that the
           caller can map rows back to sequence identifiers.
         - scoring_matrix - A substitution matrix (``Bio.Align.substitution_matrices``
           Array) or ``None`` for the identity model.
         - skip_letters - Tuple of single-character strings to ignore at
           any column where either sequence contains one of them.  Pass
           an empty tuple ``()`` to score every column.

        Returns a list of lists in lower-triangular format where row *i*
        contains *i + 1* entries and the diagonal entry ``matrix[i][i]``
        is ``0.0``.

        The returned matrix is validated for correct shape and zero
        diagonals whenever ``__debug__`` is ``True`` (the default;
        disabled only when the interpreter is invoked with ``-O``).
        """
        ...
