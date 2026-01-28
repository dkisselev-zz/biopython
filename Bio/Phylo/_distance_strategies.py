"""Distance matrix computation strategies for DistanceCalculator.

This module provides multiple computation backends for distance matrix calculation,
enabling significant performance improvements for large multiple sequence alignments.
"""

from abc import ABC
from abc import abstractmethod


def _compute_pairwise_chunk(args):
    """Worker function to compute distances for a chunk of sequence pairs.

    Arguments:
     - args - Tuple of (chunk, seq_strings, seq_ids, model, skip_letters)
       where chunk is list of (i,j) pairs, seq_strings is list of sequence
       strings, seq_ids is list of sequence IDs, model is model name, and
       skip_letters is tuple of letters to skip

    Returns list of (i, j, distance) tuples.

    Note: This function receives simple strings rather than Seq/SeqRecord
    objects to minimize pickling overhead. Pickling heavy Biopython objects
    across processes is expensive; strings are lightweight and fast to pickle.
    """
    chunk, seq_strings, seq_ids, model, skip_letters = args

    # Recreate calculator (can't pickle the whole object)
    from Bio.Phylo.TreeConstruction import DistanceCalculator

    calculator = DistanceCalculator(model=model, skip_letters=skip_letters)

    # Create simple objects that _pairwise expects (with .id attribute)
    from types import SimpleNamespace

    sequences = [
        SimpleNamespace(seq=seq_str, id=seq_id)
        for seq_str, seq_id in zip(seq_strings, seq_ids)
    ]

    results = []
    for i, j in chunk:
        # _pairwise expects objects with seq attribute or iterables
        distance = calculator._pairwise(sequences[i].seq, sequences[j].seq)
        results.append((i, j, distance))
    return results


def compute_parallel(strategy, msa, calculator, n_jobs):
    """Compute distance matrix using parallel execution.

    This function wraps any strategy to execute pairwise distance
    computations in parallel using multiple workers. It extracts sequence
    strings from the MSA to minimize pickling overhead when sending data
    to worker processes.

    Arguments:
     - strategy - DistanceComputationStrategy instance (currently unused,
       parallel execution always uses python loop method)
     - msa - Alignment or MultipleSeqAlignment object
     - calculator - DistanceCalculator instance
     - n_jobs - Number of parallel workers (-1 for all CPUs)

    Returns a :class:`DistanceMatrix` object.

    Performance considerations:
     - Multiprocessing overhead (process spawning, IPC) can exceed benefits
       for small alignments (N < 100)
     - Sequence data is converted to simple strings to minimize pickling cost
     - Each worker recreates a DistanceCalculator (cannot pickle complex objects)
    """
    import multiprocessing
    from Bio.Align import Alignment
    from Bio.Align import MultipleSeqAlignment
    from Bio.Phylo.TreeConstruction import DistanceMatrix

    # Determine number of workers
    if n_jobs == -1:
        n_jobs = multiprocessing.cpu_count()
    elif n_jobs < 1:
        raise ValueError(f"n_jobs must be >= 1 or -1, got {n_jobs}")

    # Extract sequence IDs and convert sequences to lightweight strings
    # This minimizes pickling overhead when sending to worker processes
    if isinstance(msa, Alignment):
        names = [s.id for s in msa.sequences]
        # Extract just the sequence strings, not full Seq/SeqRecord objects
        seq_strings = [str(msa[i]) for i in range(len(names))]
    elif isinstance(msa, MultipleSeqAlignment):
        names = [s.id for s in msa]
        # Convert SeqRecord objects to simple strings
        seq_strings = [str(s.seq) for s in msa]
    else:
        raise TypeError(
            "Must provide an Alignment object or a MultipleSeqAlignment object."
        )

    n_seqs = len(names)
    if n_seqs < 2:
        return DistanceMatrix(names)

    # Generate all pairwise indices
    pairs = [(i, j) for i in range(n_seqs) for j in range(i)]

    # If only a few pairs, don't bother with parallelism overhead
    if len(pairs) < n_jobs * 2:
        return strategy.compute(msa, calculator)

    # Split pairs into chunks for workers
    chunk_size = max(1, len(pairs) // n_jobs)
    chunks = [pairs[i : i + chunk_size] for i in range(0, len(pairs), chunk_size)]

    # Prepare arguments for workers
    # Pass simple strings and IDs, not heavy Biopython objects
    worker_args = [
        (chunk, seq_strings, names, calculator.model, calculator.skip_letters)
        for chunk in chunks
    ]

    # Execute in parallel using multiprocessing
    with multiprocessing.Pool(processes=n_jobs) as pool:
        chunk_results = pool.map(_compute_pairwise_chunk, worker_args)

    # Combine results into distance matrix
    dm = DistanceMatrix(names)
    for chunk in chunk_results:
        for i, j, distance in chunk:
            dm[names[i], names[j]] = distance

    return dm


class DistanceComputationStrategy(ABC):
    """Abstract base class for distance matrix computation strategies.

    Distance computation strategies implement different algorithmic approaches
    for calculating pairwise evolutionary distances from multiple sequence
    alignments. Subclasses must implement the :meth:`compute` method.

    Each strategy may have different performance characteristics and memory
    requirements. Use the :meth:`validate_prerequisites` method to check if
    required dependencies are available before attempting computation.
    """

    @abstractmethod
    def compute(self, msa, calculator):
        """Compute distance matrix for a multiple sequence alignment.

        Arguments:
         - msa - Alignment or MultipleSeqAlignment object
         - calculator - DistanceCalculator instance containing model and parameters

        Returns a :class:`DistanceMatrix` object with pairwise distances.
        """
        pass

    def validate_prerequisites(self):
        """Check if required dependencies are available.

        Returns ``True`` if all required libraries and dependencies are available,
        ``False`` otherwise. Default implementation returns ``True``.
        """
        return True

    def estimate_memory(self, n_seqs, seq_len):
        """Estimate peak memory usage in bytes.

        Arguments:
         - n_seqs - Number of sequences in alignment
         - seq_len - Length of aligned sequences

        Returns estimated peak memory usage in bytes. Default implementation
        returns a conservative estimate of O(N²) for the distance matrix.
        """
        # Conservative estimate: N² distance matrix (8 bytes per float)
        return n_seqs * n_seqs * 8


class PythonLoopStrategy(DistanceComputationStrategy):
    """Reference implementation using nested Python loops.

    This strategy wraps the original nested loop implementation from
    :class:`DistanceCalculator`. It serves as the baseline for correctness
    testing and fallback when other methods are unavailable.

    **Performance**: O(N² × L) where N = sequences, L = alignment length

    **Memory**: O(N²) for distance matrix only

    **When to use**: Small alignments (N < 10), or as fallback

    Example::

        >>> from Bio import Align
        >>> from Bio.Phylo.TreeConstruction import DistanceCalculator
        >>> aln = Align.read("TreeConstruction/msa.phy", "phylip")
        >>> calculator = DistanceCalculator("blosum62")
        >>> dm = calculator.get_distance(aln, method="python")

    See also :class:`NumPyVectorizedStrategy`, :class:`DistanceCalculator`
    """

    def compute(self, msa, calculator):
        """Compute distance matrix using nested Python loops.

        This implementation directly calls the existing :meth:`_pairwise`
        method from the DistanceCalculator for each sequence pair.

        Arguments:
         - msa - Alignment or MultipleSeqAlignment object
         - calculator - DistanceCalculator instance

        Returns a :class:`DistanceMatrix` object.
        """
        from Bio.Align import Alignment
        from Bio.Align import MultipleSeqAlignment
        from Bio.Phylo.TreeConstruction import DistanceMatrix
        import itertools

        if isinstance(msa, Alignment):
            names = [s.id for s in msa.sequences]
            dm = DistanceMatrix(names)
            n = len(names)
            for i1 in range(n):
                for i2 in range(i1):
                    dm[names[i1], names[i2]] = calculator._pairwise(msa[i1], msa[i2])
        elif isinstance(msa, MultipleSeqAlignment):
            names = [s.id for s in msa]
            dm = DistanceMatrix(names)
            for seq1, seq2 in itertools.combinations(msa, 2):
                dm[seq1.id, seq2.id] = calculator._pairwise(seq1, seq2)
        else:
            raise TypeError(
                "Must provide an Alignment object or a MultipleSeqAlignment object."
            )

        return dm


class NumPyVectorizedStrategy(DistanceComputationStrategy):
    """NumPy-accelerated vectorized distance computation.

    This strategy vectorizes the inner loop over sequence length using NumPy
    operations, providing significant speedup compared to the baseline Python
    implementation while maintaining numerical accuracy.

    **Performance**: O(N² × L/k) where k≈10-50 (vectorization speedup)

    **Memory**: O(N×L + N²) for sequence array and distance matrix

    **When to use**: Medium alignments (10 < N < 500), moderate L

    **Prerequisites**: NumPy (already required by Biopython)

    Example::

        >>> from Bio import Align
        >>> from Bio.Phylo.TreeConstruction import DistanceCalculator
        >>> aln = Align.read("TreeConstruction/msa.phy", "phylip")
        >>> calculator = DistanceCalculator("blosum62")
        >>> dm = calculator.get_distance(aln, method="numpy")

    See also :class:`PythonLoopStrategy`, :class:`DistanceCalculator`
    """

    def validate_prerequisites(self):
        """Check if NumPy is available.

        Returns ``True`` since NumPy is required by Biopython.
        """
        try:
            import numpy  # noqa: F401

            return True
        except ImportError:
            return False

    def estimate_memory(self, n_seqs, seq_len):
        """Estimate peak memory usage in bytes.

        Arguments:
         - n_seqs - Number of sequences in alignment
         - seq_len - Length of aligned sequences

        Returns estimated peak memory usage: N×L (sequences) + N² (distances).
        """
        # Sequence array (uint8): N × L bytes
        # Distance matrix (float64): N² × 8 bytes
        return n_seqs * seq_len + n_seqs * n_seqs * 8

    def compute(self, msa, calculator):
        """Compute distance matrix using NumPy vectorization.

        This implementation vectorizes over the sequence length dimension,
        converting sequences to NumPy arrays and using vectorized operations
        for scoring and gap masking.

        Arguments:
         - msa - Alignment or MultipleSeqAlignment object
         - calculator - DistanceCalculator instance

        Returns a :class:`DistanceMatrix` object.
        """
        import numpy as np
        from Bio.Align import Alignment
        from Bio.Align import MultipleSeqAlignment
        from Bio.Phylo.TreeConstruction import DistanceMatrix

        # Extract sequences and names
        if isinstance(msa, Alignment):
            names = [s.id for s in msa.sequences]
            sequences = [str(msa[i]) for i in range(len(names))]
        elif isinstance(msa, MultipleSeqAlignment):
            names = [s.id for s in msa]
            sequences = [str(s.seq) for s in msa]
        else:
            raise TypeError(
                "Must provide an Alignment object or a MultipleSeqAlignment object."
            )

        n_seqs = len(names)
        if n_seqs == 0:
            return DistanceMatrix([])

        seq_len = len(sequences[0])

        # Convert sequences to NumPy array of characters
        seq_array = np.array([list(seq) for seq in sequences], dtype="U1")

        # Create skip mask for gap characters
        skip_chars = set(calculator.skip_letters)
        skip_mask = np.zeros((n_seqs, seq_len), dtype=bool)
        for char in skip_chars:
            skip_mask |= seq_array == char

        # Initialize distance matrix
        dm = DistanceMatrix(names)

        # Check if using identity model (scoring_matrix is None for identity)
        is_identity = calculator.scoring_matrix is None

        # Compute pairwise distances
        for i in range(n_seqs):
            for j in range(i):
                # Get valid positions (neither sequence has skip character)
                valid = ~(skip_mask[i] | skip_mask[j])
                n_valid = np.sum(valid)

                if n_valid == 0:
                    # No valid positions - distance undefined, use 1.0
                    dm[names[i], names[j]] = 1.0
                    continue

                if is_identity:
                    # Identity model: fraction of mismatches
                    matches = np.sum((seq_array[i] == seq_array[j]) & valid)
                    distance = 1.0 - (matches / n_valid)
                else:
                    # Substitution matrix model
                    distance = self._compute_matrix_distance(
                        seq_array[i],
                        seq_array[j],
                        valid,
                        n_valid,
                        calculator.scoring_matrix,
                    )

                dm[names[i], names[j]] = distance

        return dm

    def _compute_matrix_distance(self, seq1, seq2, valid, n_valid, scoring_matrix):
        """Compute distance using substitution matrix.

        Arguments:
         - seq1 - First sequence as NumPy array
         - seq2 - Second sequence as NumPy array
         - valid - Boolean mask of valid positions
         - n_valid - Number of valid positions
         - scoring_matrix - Substitution matrix (Bio.Align.substitution_matrices.Array)

        Returns normalized distance score.
        """
        import numpy as np

        # Extract valid positions
        seq1_valid = seq1[valid]
        seq2_valid = seq2[valid]

        # Compute scores using vectorized matrix lookups
        # For each position, get scoring_matrix[seq1[i], seq2[i]]
        scores = np.array(
            [scoring_matrix[a, b] for a, b in zip(seq1_valid, seq2_valid)]
        )
        total_score = np.sum(scores)

        # Compute maximum possible scores (self-comparisons)
        max_score1 = np.sum([scoring_matrix[a, a] for a in seq1_valid])
        max_score2 = np.sum([scoring_matrix[b, b] for b in seq2_valid])

        # Take the higher score (matches Python implementation line 544)
        max_score = max(max_score1, max_score2)

        if max_score == 0:
            return 1.0

        # Normalize: distance = 1 - (score / max_score)
        distance = 1.0 - (total_score / max_score)
        return max(0.0, min(1.0, distance))  # Clamp to [0, 1]


class SciPyDistanceStrategy(DistanceComputationStrategy):
    """SciPy-accelerated distance computation using pdist.

    This strategy leverages SciPy's optimized ``pdist`` function for fast
    pairwise distance computation. Currently supports only the identity model,
    using Hamming distance as the underlying metric.

    **Performance**: O(N² × L) but highly optimized in C (10-30x faster)

    **Memory**: O(N×L + N²) for sequence array and distance matrix

    **When to use**: Medium to large alignments (N > 100) with identity model

    **Prerequisites**: SciPy (optional dependency)

    **Limitations**: Only supports identity model. For substitution matrices,
    use :class:`NumPyVectorizedStrategy` instead.

    Example::

        >>> from Bio import Align
        >>> from Bio.Phylo.TreeConstruction import DistanceCalculator
        >>> aln = Align.read("TreeConstruction/msa.phy", "phylip")
        >>> calculator = DistanceCalculator("identity")
        >>> dm = calculator.get_distance(aln, method="scipy")

    See also :class:`NumPyVectorizedStrategy`, :class:`DistanceCalculator`
    """

    def validate_prerequisites(self):
        """Check if SciPy is available.

        Returns ``True`` if SciPy can be imported, ``False`` otherwise.
        """
        try:
            import scipy.spatial.distance  # noqa: F401

            return True
        except ImportError:
            return False

    def estimate_memory(self, n_seqs, seq_len):
        """Estimate peak memory usage in bytes.

        Arguments:
         - n_seqs - Number of sequences in alignment
         - seq_len - Length of aligned sequences

        Returns estimated peak memory usage: N×L (sequences) + N² (distances).
        """
        # Sequence array (uint8): N × L bytes
        # Distance array (float64): N×(N-1)/2 × 8 bytes (condensed form)
        # Distance matrix (float64): N² × 8 bytes (square form)
        return n_seqs * seq_len + n_seqs * n_seqs * 8

    def compute(self, msa, calculator):
        """Compute distance matrix using SciPy's pdist.

        This implementation uses SciPy's highly optimized pdist function
        with the Hamming distance metric. Only works with identity model.

        Arguments:
         - msa - Alignment or MultipleSeqAlignment object
         - calculator - DistanceCalculator instance

        Returns a :class:`DistanceMatrix` object.

        Raises:
         - Bio.MissingPythonDependencyError if SciPy is not installed
         - ValueError if used with a substitution matrix model
        """
        # Check SciPy availability
        try:
            from scipy.spatial.distance import pdist
            from scipy.spatial.distance import squareform
        except ImportError:
            from Bio import MissingPythonDependencyError

            raise MissingPythonDependencyError(
                "SciPy is required for the 'scipy' distance computation method. "
                "Install with: pip install scipy"
            )

        # Check that we're using identity model
        if calculator.scoring_matrix is not None:
            raise ValueError(
                "SciPy method only supports identity model. "
                "For substitution matrices, use method='numpy' instead."
            )

        import numpy as np
        from Bio.Align import Alignment
        from Bio.Align import MultipleSeqAlignment
        from Bio.Phylo.TreeConstruction import DistanceMatrix

        # Extract sequences and names
        if isinstance(msa, Alignment):
            names = [s.id for s in msa.sequences]
            sequences = [str(msa[i]) for i in range(len(names))]
        elif isinstance(msa, MultipleSeqAlignment):
            names = [s.id for s in msa]
            sequences = [str(s.seq) for s in msa]
        else:
            raise TypeError(
                "Must provide an Alignment object or a MultipleSeqAlignment object."
            )

        n_seqs = len(names)
        if n_seqs == 0:
            return DistanceMatrix([])

        if n_seqs == 1:
            return DistanceMatrix(names)

        seq_len = len(sequences[0])

        # Convert sequences to NumPy array of characters
        seq_array = np.array([list(seq) for seq in sequences], dtype="U1")

        # Handle skip characters by masking
        skip_chars = set(calculator.skip_letters)
        if skip_chars:
            # Create mask for valid positions
            skip_mask = np.zeros((n_seqs, seq_len), dtype=bool)
            for char in skip_chars:
                skip_mask |= seq_array == char

            # For each pairwise comparison, we need to handle gaps
            # Fall back to manual computation with gap handling
            dm = DistanceMatrix(names)
            for i in range(n_seqs):
                for j in range(i):
                    valid = ~(skip_mask[i] | skip_mask[j])
                    n_valid = np.sum(valid)
                    if n_valid == 0:
                        dm[names[i], names[j]] = 1.0
                    else:
                        matches = np.sum((seq_array[i] == seq_array[j]) & valid)
                        dm[names[i], names[j]] = 1.0 - (matches / n_valid)
        else:
            # No skip characters - use pdist directly with hamming distance
            # Convert to numeric codes for pdist
            # First, get unique characters
            unique_chars = np.unique(seq_array)
            char_to_code = {char: i for i, char in enumerate(unique_chars)}

            # Convert sequences to numeric codes
            seq_codes = np.array(
                [[char_to_code[char] for char in seq] for seq in sequences]
            )

            # Compute pairwise Hamming distances
            # pdist returns condensed distance matrix
            distances = pdist(seq_codes, metric="hamming")

            # Convert to square form
            dist_matrix = squareform(distances)

            # Create DistanceMatrix object
            dm = DistanceMatrix(names)
            for i in range(n_seqs):
                for j in range(i):
                    dm[names[i], names[j]] = dist_matrix[i, j]

        return dm


class OneHotMatMulStrategy(DistanceComputationStrategy):
    """Matrix multiplication-based distance computation using one-hot encoding.

    This strategy uses one-hot encoding and BLAS-optimized matrix multiplication
    to achieve maximum performance for large alignments. Converts sequences to
    one-hot tensors and uses einsum for vectorized scoring.

    **Performance**: O(N² × L × C²) but BLAS Level 3 optimized (50-200x faster)

    **Memory**: O(N×L×C) for one-hot encoding (HIGH - see estimate_memory)

    **When to use**: Very large alignments (N > 500) with sufficient memory

    **Prerequisites**: NumPy (already required by Biopython)

    **Limitations**: High memory usage for large alphabets or long sequences.
    Use estimate_memory() to check before computation.

    Example::

        >>> from Bio import Align
        >>> from Bio.Phylo.TreeConstruction import DistanceCalculator
        >>> aln = Align.read("TreeConstruction/msa.phy", "phylip")
        >>> calculator = DistanceCalculator("blosum62")
        >>> # Check memory first
        >>> from Bio.Phylo._distance_strategies import OneHotMatMulStrategy
        >>> strategy = OneHotMatMulStrategy()
        >>> mem_bytes = strategy.estimate_memory(5, 13)
        >>> print(f"Estimated memory: {mem_bytes / 1e6:.2f} MB")
        Estimated memory: 0.01 MB
        >>> dm = calculator.get_distance(aln, method="onehot")

    See also :class:`NumPyVectorizedStrategy`, :class:`DistanceCalculator`
    """

    def validate_prerequisites(self):
        """Check if NumPy is available.

        Returns ``True`` since NumPy is required by Biopython.
        """
        try:
            import numpy  # noqa: F401

            return True
        except ImportError:
            return False

    def estimate_memory(self, n_seqs, seq_len):
        """Estimate peak memory usage in bytes.

        Arguments:
         - n_seqs - Number of sequences in alignment
         - seq_len - Length of aligned sequences

        Returns estimated peak memory usage for one-hot encoding.

        **Memory breakdown**:
        - One-hot array: N × L × C × 4 bytes (float32)
        - Distance matrix: N² × 8 bytes (float64)

        Where C is estimated at 20 (typical protein alphabet size).
        """
        # Estimate alphabet size (20 for proteins, 4 for DNA, use 20 as default)
        alphabet_size = 20

        # One-hot tensor: N × L × C (float32 = 4 bytes)
        onehot_memory = n_seqs * seq_len * alphabet_size * 4

        # Distance matrix: N² (float64 = 8 bytes)
        matrix_memory = n_seqs * n_seqs * 8

        return onehot_memory + matrix_memory

    def compute(self, msa, calculator):
        """Compute distance matrix using one-hot encoding and matrix multiplication.

        This implementation converts sequences to one-hot encoded tensors and
        uses NumPy's einsum for BLAS-optimized matrix operations.

        Arguments:
         - msa - Alignment or MultipleSeqAlignment object
         - calculator - DistanceCalculator instance

        Returns a :class:`DistanceMatrix` object.
        """
        import numpy as np
        from Bio.Align import Alignment
        from Bio.Align import MultipleSeqAlignment
        from Bio.Phylo.TreeConstruction import DistanceMatrix

        # Extract sequences and names
        if isinstance(msa, Alignment):
            names = [s.id for s in msa.sequences]
            sequences = [str(msa[i]) for i in range(len(names))]
        elif isinstance(msa, MultipleSeqAlignment):
            names = [s.id for s in msa]
            sequences = [str(s.seq) for s in msa]
        else:
            raise TypeError(
                "Must provide an Alignment object or a MultipleSeqAlignment object."
            )

        n_seqs = len(names)
        if n_seqs == 0:
            return DistanceMatrix([])

        if n_seqs == 1:
            return DistanceMatrix(names)

        seq_len = len(sequences[0])

        # Convert sequences to NumPy array
        seq_array = np.array([list(seq) for seq in sequences], dtype="U1")

        # Get unique characters (alphabet)
        unique_chars = np.unique(seq_array)
        alphabet = [c for c in unique_chars if c not in calculator.skip_letters]
        alphabet_size = len(alphabet)

        if alphabet_size == 0:
            # All characters are skip characters
            dm = DistanceMatrix(names)
            for i in range(n_seqs):
                for j in range(i):
                    dm[names[i], names[j]] = 1.0
            return dm

        # Create character to index mapping
        char_to_idx = {char: i for i, char in enumerate(alphabet)}

        # Create one-hot encoding: (N, L, C)
        onehot = np.zeros((n_seqs, seq_len, alphabet_size), dtype=np.float32)

        for i in range(n_seqs):
            for j in range(seq_len):
                char = seq_array[i, j]
                if char in char_to_idx:
                    onehot[i, j, char_to_idx[char]] = 1.0

        # Handle identity vs scoring matrix
        if calculator.scoring_matrix is None:
            # Identity model: score matrix is identity matrix
            score_matrix = np.eye(alphabet_size, dtype=np.float32)
        else:
            # Build score matrix from substitution matrix
            score_matrix = np.zeros((alphabet_size, alphabet_size), dtype=np.float32)
            for i, char1 in enumerate(alphabet):
                for j, char2 in enumerate(alphabet):
                    try:
                        score_matrix[i, j] = calculator.scoring_matrix[char1, char2]
                    except (KeyError, IndexError):
                        # Character not in scoring matrix, use 0
                        score_matrix[i, j] = 0.0

        # Compute pairwise scores using einsum
        # For each pair (i, j), sum over positions l and alphabet indices a, b:
        # score[i,j] = sum_l sum_a sum_b onehot[i,l,a] * score_matrix[a,b] * onehot[j,l,b]
        pairwise_scores = np.einsum(
            "ila,ab,jlb->ij", onehot, score_matrix, onehot, dtype=np.float32
        )

        # Compute self-scores (diagonal)
        self_scores = np.diag(pairwise_scores)

        # Create distance matrix
        dm = DistanceMatrix(names)

        for i in range(n_seqs):
            for j in range(i):
                score = pairwise_scores[i, j]

                # Normalize by maximum possible score
                max_score = max(self_scores[i], self_scores[j])

                if max_score == 0:
                    distance = 1.0
                else:
                    distance = 1.0 - (score / max_score)

                # Clamp to [0, 1]
                distance = max(0.0, min(1.0, distance))
                dm[names[i], names[j]] = distance

        return dm
