"""Private module: pluggable distance-computation strategies for DistanceCalculator.

This module provides alternative backends for computing pairwise sequence
distances used by ``Bio.Phylo.TreeConstruction.DistanceCalculator``.  All
backends produce numerically identical (or near-identical) results to the
original pure-Python pair loop (``method="python"``); vectorised backends
trade memory for speed.

Available methods::

    python  – reference implementation; the only one supporting ``n_jobs > 1``
    numpy   – chunked NumPy broadcasting (identity) / position-loop vectorised
              (scoring)
    scipy   – SciPy ``pdist`` with C-optimised hamming (identity) or callable
              (scoring)
    onehot  – one-hot encoding + BLAS GEMM (``X @ X.T``)

The public entry-point is ``METHOD_REGISTRY``, a dict mapping method names to
strategy classes.  ``_resolve_n_jobs`` handles the ``BIOPYTHON_DISTANCE_JOBS``
environment variable.
"""

import os
from multiprocessing import Pool

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _resolve_n_jobs(n_jobs):
    """Return effective worker count, honouring BIOPYTHON_DISTANCE_JOBS env var.

    ``n_jobs=-1`` is converted to ``os.cpu_count()``.  The environment
    variable, when set to a valid integer, overrides the *n_jobs* argument.
    """
    env = os.environ.get("BIOPYTHON_DISTANCE_JOBS")
    if env is not None:
        try:
            n_jobs = int(env)
        except ValueError:
            pass
    if n_jobs == -1:
        return os.cpu_count() or 1
    return max(1, n_jobs)


def _scoring_matrix_to_dict(scoring_matrix):
    """Convert a substitution_matrices.Array to a plain dict for pickle safety.

    ``multiprocessing`` workers receive arguments via pickle; the Array
    subclass's custom *alphabet* attribute does not survive the default
    ndarray ``__reduce__``.  A ``{(char_a, char_b): float}`` dict is
    universally picklable and fast to look up.
    """
    sm_dict = {}
    for a in scoring_matrix.alphabet:
        for b in scoring_matrix.alphabet:
            sm_dict[(a, b)] = float(scoring_matrix[a, b])
    return sm_dict


def _pairwise_python(seq1, seq2, scoring_matrix, skip_letters):
    """Compute pairwise distance between two aligned sequences.

    This is a **module-level** function so that ``multiprocessing.Pool``
    can pickle it.

    Parameters
    ----------
    seq1, seq2 : str
        Aligned sequences of equal length.
    scoring_matrix : dict or None
        ``{(char_a, char_b): score}`` dict, or ``None`` for the identity model.
    skip_letters : tuple of str
        Characters to skip (e.g. ``("-", "*")``).

    Returns
    -------
    float
        Pairwise distance in [0, 1].
    """
    if scoring_matrix is None:
        score = sum(1 for a, b in zip(seq1, seq2) if a == b)
        max_score = len(seq1)
        if max_score == 0:
            return 1
        return 1 - score / max_score

    score = 0
    max_score1 = 0
    max_score2 = 0
    for c1, c2 in zip(seq1, seq2):
        if c1 in skip_letters or c2 in skip_letters:
            continue
        try:
            s = scoring_matrix[(c1, c2)]
        except KeyError:
            if (c1, c1) not in scoring_matrix:
                raise ValueError(f"Bad letter '{c1}' in sequence") from None
            raise ValueError(f"Bad letter '{c2}' in sequence") from None
        score += s
        max_score1 += scoring_matrix[(c1, c1)]
        max_score2 += scoring_matrix[(c2, c2)]
    max_score = max(max_score1, max_score2)
    if max_score == 0:
        return 1
    return 1 - score / max_score


def _extract_lower_tri(matrix, n):
    """Extract lower-triangular entries in row-major order.

    Returns values for pairs ``(1,0), (2,0), (2,1), (3,0), …`` as a list.
    Requires NumPy.
    """
    import numpy as np

    return matrix[np.tril_indices(n, k=-1)].tolist()


def _pdist_to_lower_tri(condensed, n):
    """Reorder scipy pdist condensed vector to lower-tri row-major format.

    pdist output order (upper-tri row-major): ``(0,1), (0,2), …``
    Our output order (lower-tri row-major): ``(1,0), (2,0), (2,1), …``

    Uses vectorised index arithmetic; no N×N intermediate matrix is created.
    """
    import numpy as np

    rows, cols = np.tril_indices(n, k=-1)
    # pdist condensed index for pair (i, j) where i < j:
    #   i * (2n - i - 1) // 2 + (j - i - 1)
    # Here cols < rows, so the pdist pair is (cols, rows).
    indices = cols * (2 * n - cols - 1) // 2 + (rows - cols - 1)
    return condensed[indices].tolist()


# ---------------------------------------------------------------------------
# Base class
# ---------------------------------------------------------------------------


class _DistanceMethod:
    """Base class for distance-computation strategies (PRIVATE).

    Subclasses implement ``compute``, which receives pre-extracted sequence
    strings and returns distances in lower-triangular row-major order:
    ``(1,0), (2,0), (2,1), (3,0), (3,1), (3,2), …``
    """

    def compute(self, sequences, scoring_matrix, skip_letters, n_jobs=1):
        """Return a flat list of pairwise distances.

        Parameters
        ----------
        sequences : list of str
            Aligned sequence strings, all the same length.
        scoring_matrix : substitution_matrices.Array or None
            ``None`` selects the identity model.
        skip_letters : tuple of str
            Characters ignored in the scoring-matrix model.
        n_jobs : int
            Parallelism hint (honoured only by ``_PythonMethod``).

        Returns
        -------
        list of float
            Lower-triangular distances in row-major order.
        """
        raise NotImplementedError


# ---------------------------------------------------------------------------
# 1. Pure-Python reference (supports n_jobs > 1)
# ---------------------------------------------------------------------------


class _PythonMethod(_DistanceMethod):
    """Reference implementation: pure-Python pair loop.

    This is the only method that honours *n_jobs*.  When ``n_jobs > 1``,
    pairs are distributed via ``multiprocessing.Pool.starmap``.
    """

    def compute(self, sequences, scoring_matrix, skip_letters, n_jobs=1):
        n = len(sequences)
        sm_dict = (
            _scoring_matrix_to_dict(scoring_matrix)
            if scoring_matrix is not None
            else None
        )

        pairs = [
            (sequences[i], sequences[j], sm_dict, skip_letters)
            for i in range(1, n)
            for j in range(i)
        ]

        if n_jobs == 1:
            return [_pairwise_python(*p) for p in pairs]
        with Pool(processes=n_jobs) as pool:
            return pool.starmap(_pairwise_python, pairs)


# ---------------------------------------------------------------------------
# 2. NumPy vectorised
# ---------------------------------------------------------------------------


class _NumpyMethod(_DistanceMethod):
    """Vectorised implementation using NumPy broadcasting.

    Identity model: chunked row-broadcast over a character array, with a
    512 MB memory budget to avoid materialising the full ``(N, N, L)`` tensor.

    Scoring model: iterates over alignment positions (L loop) but vectorises
    the N×N pair dimension at each step.

    *n_jobs* is ignored; NumPy/BLAS releases the GIL internally.
    """

    _MEMORY_BUDGET = 512 * 1024 * 1024  # 512 MB

    def compute(self, sequences, scoring_matrix, skip_letters, n_jobs=1):
        import numpy as np  # noqa: F811

        n = len(sequences)
        L = len(sequences[0])

        if scoring_matrix is None:
            return self._identity(sequences, n, L)
        return self._scoring(sequences, scoring_matrix, skip_letters, n, L)

    def _identity(self, sequences, n, L):
        import numpy as np

        # Single-pass ingestion: join → encode → frombuffer → reshape.
        # Comparison on uint8 ordinals is identical to comparison on characters.
        arr = np.frombuffer(
            "".join(sequences).encode("ascii"), dtype=np.uint8
        ).reshape(n, L)

        # Adaptive chunk size so chunk * N * L bytes stays under budget
        chunk = max(1, self._MEMORY_BUDGET // (n * L))
        matches = np.zeros((n, n), dtype=np.float64)

        for start in range(0, n, chunk):
            end = min(start + chunk, n)
            matches[start:end, :] = (arr[start:end, None, :] == arr[None, :, :]).sum(
                axis=2
            )

        distances = 1.0 - matches / L
        return _extract_lower_tri(distances, n)

    def _scoring(self, sequences, scoring_matrix, skip_letters, n, L):
        import numpy as np

        alphabet = list(scoring_matrix.alphabet)
        sentinel = len(alphabet)

        # 128-entry LUT: alphabet → index, skip_letters → sentinel, else → -1.
        # Vectorised lookup replaces both validation and index-array construction.
        lut = np.full(128, -1, dtype=np.intp)
        for i, c in enumerate(alphabet):
            lut[ord(c)] = i
        for c in skip_letters:
            lut[ord(c)] = sentinel

        # Single-pass ingestion + vectorised lookup
        ord_arr = np.frombuffer(
            "".join(sequences).encode("ascii"), dtype=np.uint8
        ).reshape(n, L)
        idx_arr = lut[ord_arr]  # (N, L)

        # Vectorised validation: -1 means unrecognised character
        bad_mask = idx_arr == -1
        if bad_mask.any():
            pos = int(np.argmax(bad_mask))
            row, col = divmod(pos, L)
            raise ValueError(f"Bad letter '{sequences[row][col]}' in sequence")

        # Extended scoring matrix: sentinel row/column are zeros
        sm_size = len(alphabet) + 1
        sm = np.zeros((sm_size, sm_size), dtype=np.float64)
        for i, a in enumerate(alphabet):
            for j, b in enumerate(alphabet):
                sm[i, j] = float(scoring_matrix[a, b])

        score = np.zeros((n, n), dtype=np.float64)
        max_score1 = np.zeros((n, n), dtype=np.float64)
        max_score2 = np.zeros((n, n), dtype=np.float64)

        for k in range(L):
            chars_k = idx_arr[:, k]  # (N,)
            valid_k = chars_k != sentinel  # (N,) bool
            pair_valid = valid_k[:, None] & valid_k[None, :]  # (N, N)

            score += sm[chars_k[:, None], chars_k[None, :]] * pair_valid
            self_scores = sm[chars_k, chars_k]  # (N,)
            max_score1 += self_scores[:, None] * pair_valid
            max_score2 += self_scores[None, :] * pair_valid

        max_score = np.maximum(max_score1, max_score2)
        distances = np.where(max_score > 0, 1.0 - score / max_score, 1.0)
        return _extract_lower_tri(distances, n)


# ---------------------------------------------------------------------------
# 3. SciPy pdist
# ---------------------------------------------------------------------------


class _ScipyMethod(_DistanceMethod):
    """Implementation using SciPy's ``pdist`` (identity) and chunked NumPy
    broadcasting (scoring).

    Identity model: delegates to ``pdist(metric='hamming')``, a C-optimised
    path.  Note: hamming computes ``mismatch_count / L`` while the Python
    reference computes ``1 - match_count / L``.  These are mathematically
    identical but may differ at the last floating-point bit.

    Scoring model: chunked row-broadcast over the L-loop, materialising only
    ``(chunk, N)`` accumulators at a time.  Peak memory is bounded by
    ``_MEMORY_BUDGET``; no per-pair Python callable is involved.

    *n_jobs* is ignored; ``pdist`` / NumPy are internally optimised.
    """

    _MEMORY_BUDGET = 512 * 1024 * 1024  # 512 MB

    def compute(self, sequences, scoring_matrix, skip_letters, n_jobs=1):
        n = len(sequences)
        if scoring_matrix is None:
            return self._identity(sequences, n)
        return self._scoring(sequences, scoring_matrix, skip_letters, n)

    def _identity(self, sequences, n):
        import numpy as np
        from scipy.spatial.distance import pdist

        arr = np.array([[ord(c) for c in seq] for seq in sequences], dtype=np.int32)
        condensed = pdist(arr, metric="hamming")
        return _pdist_to_lower_tri(condensed, n)

    def _scoring(self, sequences, scoring_matrix, skip_letters, n):
        import numpy as np

        alphabet = list(scoring_matrix.alphabet)
        sentinel = len(alphabet)
        L = len(sequences[0])

        # 128-entry LUT: alphabet → index, skip_letters → sentinel, else → -1.
        # Vectorised lookup replaces both validation and index-array construction.
        lut = np.full(128, -1, dtype=np.intp)
        for i, c in enumerate(alphabet):
            lut[ord(c)] = i
        for c in skip_letters:
            lut[ord(c)] = sentinel

        # Single-pass ingestion + vectorised lookup
        ord_arr = np.frombuffer(
            "".join(sequences).encode("ascii"), dtype=np.uint8
        ).reshape(n, L)
        idx_arr = lut[ord_arr]  # (N, L)

        # Vectorised validation: -1 means unrecognised character
        bad_mask = idx_arr == -1
        if bad_mask.any():
            pos = int(np.argmax(bad_mask))
            row, col = divmod(pos, L)
            raise ValueError(f"Bad letter '{sequences[row][col]}' in sequence")

        # Extended scoring matrix: sentinel row/column are zeros
        sm_size = len(alphabet) + 1
        sm = np.zeros((sm_size, sm_size), dtype=np.float64)
        for i, a in enumerate(alphabet):
            for j, b in enumerate(alphabet):
                sm[i, j] = float(scoring_matrix[a, b])

        # Chunked row-broadcast: each chunk materialises (chunk, N)
        # accumulators rather than a full (N, N) matrix.  Chunk size is chosen
        # so that the three persistent accumulators plus the largest L-loop
        # temporary stay within the memory budget.
        bytes_per_row = n * 8 * 4  # score, max_score1, max_score2 + peak temp
        chunk_size = max(1, self._MEMORY_BUDGET // bytes_per_row)
        chunk_size = min(chunk_size, n)

        col_indices = np.arange(n)  # reused for lower-tri mask each chunk
        flat = []

        for i_start in range(0, n, chunk_size):
            i_end = min(i_start + chunk_size, n)
            chunk_n = i_end - i_start

            score = np.zeros((chunk_n, n), dtype=np.float64)
            max_score1 = np.zeros((chunk_n, n), dtype=np.float64)
            max_score2 = np.zeros((chunk_n, n), dtype=np.float64)

            for k in range(L):
                row_chars = idx_arr[i_start:i_end, k]  # (chunk_n,)
                col_chars = idx_arr[:, k]  # (N,)

                row_valid = row_chars != sentinel  # (chunk_n,) bool
                col_valid = col_chars != sentinel  # (N,) bool
                pair_valid = row_valid[:, None] & col_valid[None, :]  # (chunk_n, N)

                score += sm[row_chars[:, None], col_chars[None, :]] * pair_valid
                row_self = sm[row_chars, row_chars]  # (chunk_n,)
                col_self = sm[col_chars, col_chars]  # (N,)
                max_score1 += row_self[:, None] * pair_valid
                max_score2 += col_self[None, :] * pair_valid

            max_score = np.maximum(max_score1, max_score2)
            dists = np.where(max_score > 0, 1.0 - score / max_score, 1.0)

            # Extract lower-triangular entries in row-major order.
            # For global row i we need columns [0, i).  The boolean mask
            # selects exactly those entries; dists[mask] returns them in
            # row-major order, which is the required flat layout.
            row_globals = np.arange(i_start, i_end)[:, None]  # (chunk_n, 1)
            mask = col_indices[None, :] < row_globals  # (chunk_n, N)
            flat.extend(dists[mask].tolist())

        return flat


# ---------------------------------------------------------------------------
# 4. One-hot + GEMM
# ---------------------------------------------------------------------------


class _OneHotMethod(_DistanceMethod):
    """One-hot encoding with block-wise BLAS matrix multiplication.

    Identity model: row-blocks of one-hot vectors are generated on the fly
    and paired with column-blocks via GEMM (``X_i @ X_j.T``).  Only two
    blocks of shape ``(chunk, L*|A|)`` coexist at any time; peak memory is
    bounded by ``_MEMORY_BUDGET``.  The alphabet is derived from the sequences
    themselves, so gaps are included and gap–gap matches score 1 (correct for
    identity semantics).

    Scoring model: each row block is scoring-transformed (``X_s = X @ S``)
    then paired with raw column-block one-hot via GEMM.  ``max_score`` is
    accumulated block-wise via ``self_scores @ valid.T``.  The entire
    computation is free of explicit position loops.

    *n_jobs* is ignored; BLAS GEMM is a single optimised call per block.
    """

    _MEMORY_BUDGET = 512 * 1024 * 1024  # 512 MB

    def compute(self, sequences, scoring_matrix, skip_letters, n_jobs=1):
        n = len(sequences)
        L = len(sequences[0])

        if scoring_matrix is None:
            return self._identity(sequences, n, L)
        return self._scoring(sequences, scoring_matrix, skip_letters, n, L)

    def _identity(self, sequences, n, L):
        import numpy as np

        # Single-pass ingestion
        ord_arr = np.frombuffer(
            "".join(sequences).encode("ascii"), dtype=np.uint8
        ).reshape(n, L)

        # Derive alphabet from unique byte values present
        unique_ords = np.unique(ord_arr)
        A = len(unique_ords)

        # LUT: byte value → alphabet index
        lut = np.zeros(128, dtype=np.intp)
        for i, o in enumerate(unique_ords):
            lut[int(o)] = i

        # Chunk size: two one-hot blocks (X_i, X_j) of shape (chunk, L*A)
        # coexist during off-diagonal GEMM; the (chunk, N) row accumulator
        # is small by comparison.
        chunk = max(1, self._MEMORY_BUDGET // (2 * L * A * 8))
        chunk = min(chunk, n)

        # Reusable position-offset array for one-hot column indices
        pos_offsets = np.arange(L).reshape(1, L) * A  # (1, L)

        def _onehot_block(start, end):
            """Vectorised one-hot for sequences[start:end]."""
            m = end - start
            idx = lut[ord_arr[start:end]]  # (m, L)
            rows = np.repeat(np.arange(m), L)
            cols = (pos_offsets + idx).ravel()
            X = np.zeros((m, L * A), dtype=np.float64)
            X[rows, cols] = 1.0
            return X

        col_indices = np.arange(n)
        flat = []

        for i_start in range(0, n, chunk):
            i_end = min(i_start + chunk, n)
            chunk_i = i_end - i_start

            X_i = _onehot_block(i_start, i_end)  # (chunk_i, L*A)
            sim_row = np.zeros((chunk_i, n), dtype=np.float64)

            for j_start in range(0, i_end, chunk):
                j_end = min(j_start + chunk, i_end)

                if j_start == i_start:
                    # Diagonal block: reuse X_i (identity GEMM is symmetric)
                    sim_row[:, j_start:j_end] = X_i @ X_i.T
                else:
                    X_j = _onehot_block(j_start, j_end)
                    sim_row[:, j_start:j_end] = X_i @ X_j.T
                    del X_j

            del X_i

            # Convert to distances and extract lower-triangular entries.
            # mask[local_i, j] is True iff j < (i_start + local_i).
            distances_row = 1.0 - sim_row / L
            row_globals = np.arange(i_start, i_end)[:, None]  # (chunk_i, 1)
            mask = col_indices[None, :] < row_globals  # (chunk_i, N)
            flat.extend(distances_row[mask].tolist())

        return flat

    def _scoring(self, sequences, scoring_matrix, skip_letters, n, L):
        import numpy as np

        alphabet = list(scoring_matrix.alphabet)
        A = len(alphabet)
        sentinel = A  # index reserved for skip_letters

        # LUT: alphabet → index, skip_letters → sentinel, else → -1
        lut = np.full(128, -1, dtype=np.intp)
        for i, c in enumerate(alphabet):
            lut[ord(c)] = i
        for c in skip_letters:
            lut[ord(c)] = sentinel

        # Single-pass ingestion + vectorised lookup
        ord_arr = np.frombuffer(
            "".join(sequences).encode("ascii"), dtype=np.uint8
        ).reshape(n, L)
        idx_arr = lut[ord_arr]  # (N, L)

        # Vectorised validation: -1 means unrecognised character
        bad_mask = idx_arr == -1
        if bad_mask.any():
            pos = int(np.argmax(bad_mask))
            row, col = divmod(pos, L)
            raise ValueError(f"Bad letter '{sequences[row][col]}' in sequence")

        # Scoring matrix as numpy array
        S = np.zeros((A, A), dtype=np.float64)
        for i, a in enumerate(alphabet):
            for j, b in enumerate(alphabet):
                S[i, j] = float(scoring_matrix[a, b])
        diag_S = np.diag(S)  # (A,)

        # Chunk size: during off-diagonal GEMM, X_s_i (chunk, L*A) and
        # X_j (chunk, L*A) coexist.  X_i is freed after X_s_i is produced,
        # so peak one-hot memory is 2 × chunk × L × A × 8.
        chunk = max(1, self._MEMORY_BUDGET // (2 * L * A * 8))
        chunk = min(chunk, n)

        def _encode_block(start, end):
            """One-hot + validity mask + self-scores for sequences[start:end]."""
            m = end - start
            idx = idx_arr[start:end]  # (m, L) — already validated
            valid_mask = idx != sentinel  # (m, L)
            idx_safe = np.where(valid_mask, idx, 0)

            rows = np.repeat(np.arange(m), L)
            poss = np.tile(np.arange(L), m)
            flat_valid = valid_mask.ravel()

            X = np.zeros((m, L, A), dtype=np.float64)
            X[rows[flat_valid], poss[flat_valid], idx_safe.ravel()[flat_valid]] = 1.0

            valid = valid_mask.astype(np.float64)  # (m, L)
            self_scores = np.where(valid_mask, diag_S[idx_safe], 0.0)  # (m, L)
            return X, valid, self_scores

        col_indices = np.arange(n)
        flat = []

        for i_start in range(0, n, chunk):
            i_end = min(i_start + chunk, n)
            chunk_i = i_end - i_start

            # Encode row block, produce scoring-transformed one-hot, then
            # free the raw X_i — only X_s_i is needed on the row side.
            X_i, valid_i, ss_i = _encode_block(i_start, i_end)
            X_s_i = np.matmul(X_i, S)  # (chunk_i, L, A)
            del X_i  # frees chunk × L × A; X_s_i is independent (matmul output)
            X_s_i_flat = X_s_i.reshape(chunk_i, L * A)

            sim_row = np.zeros((chunk_i, n), dtype=np.float64)
            ms1_row = np.zeros((chunk_i, n), dtype=np.float64)
            ms2_row = np.zeros((chunk_i, n), dtype=np.float64)

            for j_start in range(0, i_end, chunk):
                j_end = min(j_start + chunk, i_end)

                # Column block is always raw one-hot: the GEMM is
                # X_s_i @ X_j.T (transformed rows × raw columns).
                # The diagonal block cannot reuse X_s_i, so it re-encodes.
                X_j, valid_j, ss_j = _encode_block(j_start, j_end)
                X_j_flat = X_j.reshape(j_end - j_start, L * A)
                del X_j  # flat view keeps the data alive for the GEMM

                sim_row[:, j_start:j_end] = X_s_i_flat @ X_j_flat.T
                ms1_row[:, j_start:j_end] = ss_i @ valid_j.T
                ms2_row[:, j_start:j_end] = valid_i @ ss_j.T

                del X_j_flat, valid_j, ss_j

            del X_s_i, X_s_i_flat, valid_i, ss_i

            max_score = np.maximum(ms1_row, ms2_row)
            dists_row = np.where(max_score > 0, 1.0 - sim_row / max_score, 1.0)

            # Extract lower-triangular entries in row-major order.
            row_globals = np.arange(i_start, i_end)[:, None]  # (chunk_i, 1)
            mask = col_indices[None, :] < row_globals  # (chunk_i, N)
            flat.extend(dists_row[mask].tolist())

        return flat


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

METHOD_REGISTRY = {
    "python": _PythonMethod,
    "numpy": _NumpyMethod,
    "scipy": _ScipyMethod,
    "onehot": _OneHotMethod,
}
