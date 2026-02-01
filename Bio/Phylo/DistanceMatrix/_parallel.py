# Copyright 2024 by Biopython contributors. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Parallelism utilities for distance matrix computation.

This module provides utilities for controlling parallel execution across
different distance calculation strategies.
"""

import os

BIOPYTHON_NUM_THREADS_ENV = "BIOPYTHON_NUM_THREADS"


def get_n_jobs(n_jobs: int | None = None) -> int:
    """Determine number of parallel jobs.

    The priority order is:
    1. Explicit n_jobs argument
    2. BIOPYTHON_NUM_THREADS environment variable
    3. Default: 1 (serial execution)

    :param n_jobs: Number of parallel workers.
        - n_jobs=1: Serial execution (default)
        - n_jobs=-1: Use all available CPUs
        - n_jobs>1: Use specified number of workers
    :returns: Number of parallel workers to use.
    :raises ValueError: If n_jobs is 0.
    """
    if n_jobs is not None:
        if n_jobs == 0:
            raise ValueError("n_jobs cannot be 0")
        if n_jobs < 0:
            return os.cpu_count() or 1
        return n_jobs

    env_value = os.environ.get(BIOPYTHON_NUM_THREADS_ENV)
    if env_value:
        try:
            return max(1, int(env_value))
        except ValueError:
            pass
    return 1


def configure_blas_threads(n_threads: int) -> bool:
    """Configure BLAS/OpenMP threading for NumPy operations.

    Attempts to limit BLAS threading at runtime using threadpoolctl if
    available. Falls back to setting environment variables, which only
    take effect if set before NumPy/BLAS libraries are first imported.

    :param n_threads: Number of threads to use for BLAS operations.
    :returns: True if runtime thread control was successful (threadpoolctl),
        False if only environment variables were set (may not take effect).
    """
    # Prioritize threadpoolctl for runtime control (works even after NumPy import)
    try:
        import threadpoolctl

        threadpoolctl.threadpool_limits(limits=n_threads, user_api="blas")
        return True
    except ImportError:
        pass

    # Fallback: set environment variables for BLAS threading.
    # NOTE: These only take effect if set BEFORE NumPy/BLAS is imported.
    # In a running process where NumPy is already loaded, these have no effect.
    # They are useful for child processes spawned via multiprocessing.
    os.environ["OMP_NUM_THREADS"] = str(n_threads)
    os.environ["MKL_NUM_THREADS"] = str(n_threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(n_threads)
    return False
