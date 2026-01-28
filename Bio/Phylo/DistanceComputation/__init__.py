# Copyright 2024 The Biopython Contributors.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Pluggable distance computation strategies for DistanceCalculator.

This sub-package provides alternative backends for pairwise sequence
distance calculation.  Strategies are resolved by name through the
``_resolve_strategy`` factory and consumed exclusively by
``DistanceCalculator`` in ``Bio.Phylo.TreeConstruction``.

Available method names: ``"numpy"``, ``"scipy"``, ``"onehot"``.
The ``"python"`` method is handled by the legacy loop in
``DistanceCalculator`` itself and does not use this package.

See Also: the Phylo_ wiki and the chapter_phylo_ tutorial.

.. _Phylo: http://biopython.org/wiki/Phylo
.. _chapter_phylo: https://biopython.org/docs/latest/Tutorial/chapter_phylo.html
"""

import os

# Registry of accelerated backends supported by this factory.
# "python" is intentionally absent — it is handled by the legacy loop
# inside DistanceCalculator and never reaches _resolve_strategy.
# Each entry is a callable ``() -> type`` that performs the lazy import
# and returns the strategy class.  This keeps the top-level import cost
# near zero and makes adding a new backend a one-line registry entry.
_STRATEGY_REGISTRY: dict[str, str] = {
    "numpy": "Bio.Phylo.DistanceComputation._numpy_strategy.NumpyDistanceStrategy",
    "scipy": "Bio.Phylo.DistanceComputation._scipy_strategy.ScipyDistanceStrategy",
    "onehot": "Bio.Phylo.DistanceComputation._onehot_strategy.OneHotDistanceStrategy",
}

# Public frozenset used by TreeConstruction for the full method list
# (including "python", which that module handles itself).
_VALID_METHODS = frozenset(("python",) + tuple(_STRATEGY_REGISTRY))


def _parse_n_jobs(n_jobs: int | None) -> int | None:
    """Resolve the effective worker count from argument and env override (PRIVATE).

    Arguments:
     - n_jobs - Caller-supplied worker count, or None for auto-detect.

    The ``BIOPYTHON_DIST_JOBS`` environment variable overrides *n_jobs*
    when set.  The resolved value must be a positive integer or None.

    Raises ValueError if the environment variable is non-numeric or
    resolves to a value less than 1.
    """
    env_jobs = os.environ.get("BIOPYTHON_DIST_JOBS")
    if env_jobs is not None:
        try:
            n_jobs = int(env_jobs)
        except ValueError:
            raise ValueError(
                f"BIOPYTHON_DIST_JOBS must be a positive integer, got '{env_jobs}'"
            ) from None
        if n_jobs < 1:
            raise ValueError(
                f"BIOPYTHON_DIST_JOBS must be a positive integer, got {n_jobs}"
            )
    return n_jobs


def _resolve_strategy(method, n_jobs=None):
    """Return an instantiated strategy for the given method name (PRIVATE).

    Arguments:
     - method - One of the keys in ``_STRATEGY_REGISTRY``
       (currently ``"numpy"``, ``"scipy"``, ``"onehot"``).
     - n_jobs - Number of parallel workers.  None means auto-detect.
       Overridden by the ``BIOPYTHON_DIST_JOBS`` environment variable
       when that variable is set.  Must be a positive integer or None.

    Raises ValueError for unrecognised method names or invalid *n_jobs*.
    """
    n_jobs = _parse_n_jobs(n_jobs)

    if method not in _STRATEGY_REGISTRY:
        raise ValueError(
            f"Unknown method '{method}'. "
            f"Accelerated methods: {tuple(_STRATEGY_REGISTRY)}"
        )

    # Lazy-import the class via the dotted path in the registry.
    module_path, class_name = _STRATEGY_REGISTRY[method].rsplit(".", 1)
    import importlib

    module = importlib.import_module(module_path)
    strategy_cls = getattr(module, class_name)
    return strategy_cls(n_jobs=n_jobs)
