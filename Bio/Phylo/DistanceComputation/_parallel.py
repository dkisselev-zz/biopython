# Copyright 2024 The Biopython Contributors.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Parallel execution helpers for distance strategies (PRIVATE)."""

import os
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import ThreadPoolExecutor
from collections.abc import Callable
from typing import Optional


def _auto_n_jobs() -> int:
    """Return a sensible default worker count (PRIVATE).

    Capped at 8 to avoid over-subscription on shared machines.
    """
    return min(os.cpu_count() or 1, 8)


class _ParallelMixin:
    """Mixin providing parallel pair-wise computation (PRIVATE).

    Subclasses that inherit this mixin gain access to ``_compute_pairs``,
    which distributes a pairwise distance function across workers using
    either threading (for GIL-releasing backends such as NumPy) or
    multiprocessing (for CPU-bound pure-Python routines).

    Sequential execution is used automatically when fewer than 50 pairs
    are requested or when ``n_jobs`` resolves to 1, avoiding the overhead
    of executor creation for small matrices.
    """

    _n_jobs: int | None

    def _compute_pairs(
        self,
        pair_func: Callable[[int, int], float],
        pairs: list[tuple[int, int]],
        use_threads: bool = True,
    ) -> list[tuple[int, int, float]]:
        """Execute pair_func over pairs, optionally in parallel (PRIVATE).

        Arguments:
         - pair_func - Callable ``(i, j) -> float`` returning the pairwise
           distance for sequence indices *i* and *j*.
         - pairs - List of ``(i, j)`` index tuples to compute.
         - use_threads - ``True`` selects ``ThreadPoolExecutor`` (appropriate
           when the inner computation releases the GIL, e.g. NumPy array
           operations).  ``False`` selects ``ProcessPoolExecutor`` for
           CPU-bound Python code.

        Returns a list of ``(i, j, distance)`` triples.
        """
        n_jobs = self._n_jobs if self._n_jobs is not None else _auto_n_jobs()

        if n_jobs <= 1 or len(pairs) < 50:
            return [(i, j, pair_func(i, j)) for i, j in pairs]

        executor_cls = ThreadPoolExecutor if use_threads else ProcessPoolExecutor
        with executor_cls(max_workers=n_jobs) as executor:
            futures = {executor.submit(pair_func, i, j): (i, j) for i, j in pairs}
            return [(futures[f][0], futures[f][1], f.result()) for f in futures]
