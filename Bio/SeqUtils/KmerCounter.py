#!/usr/bin/env python
"""K-mer counting with optional Rust acceleration.

This module provides efficient k-mer counting for DNA sequences using
either a Rust-accelerated implementation or pure Python fallback.

Classes
-------
KmerCounter
    Main class for k-mer counting with configurable parameters.

Functions
---------
count_kmers
    Convenience function for one-shot k-mer counting.
streaming_count_kmers
    Stream large FASTA/FASTQ files for memory-efficient counting.

Examples
--------
Basic k-mer counting:

>>> from Bio.SeqUtils.KmerCounter import count_kmers
>>> counts = count_kmers("ACGTACGT", k=3)
>>> counts["ACG"]
2
>>> counts["CGT"]
2

Canonical k-mers (treat k-mer and reverse complement as same):

>>> from Bio.SeqUtils.KmerCounter import KmerCounter
>>> counter = KmerCounter(k=3, canonical=True)
>>> counts = counter.count("ACGTACGT")
>>> # ACG and CGT are reverse complements, counted together

Working with Bio.Seq objects:

>>> from Bio.Seq import Seq
>>> counts = count_kmers(Seq("ATCGATCGATCG"), k=4)
>>> "ATCG" in counts
True

Streaming large files:

>>> from Bio.SeqUtils.KmerCounter import streaming_count_kmers
>>> counts = streaming_count_kmers("large_genome.fasta", k=21)  # doctest: +SKIP

Notes
-----
The Rust implementation is 20-50x faster than pure Python for typical use cases,
but requires Rust toolchain for installation. If Rust extension is not available,
the module automatically falls back to a pure Python implementation.
"""

import warnings
from collections import Counter

from Bio import BiopythonWarning

# Try to import Rust extension
try:
    from Bio.SeqUtils._kmer_counter_rust import (
        KmerCounter as RustKmerCounter,
        count_kmers_rust,
    )

    RUST_AVAILABLE = True
except ImportError:
    RUST_AVAILABLE = False
    RustKmerCounter = None
    count_kmers_rust = None
    warnings.warn(
        "Rust k-mer counter extension not available. "
        "Using pure Python fallback (will be slower for large datasets). "
        "To enable Rust acceleration, install with: pip install setuptools-rust && pip install -e .",
        BiopythonWarning,
    )


class KmerCounter:
    """Count k-mers in DNA sequences with optional Rust acceleration.

    This class provides an interface for counting k-mers (subsequences of length k)
    in DNA sequences. For k < 32, it can use a Rust implementation with 2-bit packing
    and parallel processing for significant speedup over pure Python.

    Parameters
    ----------
    k : int
        K-mer length. Must be > 0. For Rust acceleration, must be < 32.
    canonical : bool, optional
        If True, treat a k-mer and its reverse complement as the same entity.
        The lexicographically smaller sequence is used as the canonical form.
        Default is False.
    use_rust : bool, optional
        Use Rust implementation if available. Default is True.
        Set to False to force use of Python fallback.

    Attributes
    ----------
    k : int
        The k-mer length.
    canonical : bool
        Whether using canonical k-mers.
    using_rust : bool
        Whether Rust implementation is actually being used.

    Examples
    --------
    Basic usage:

    >>> counter = KmerCounter(k=3)
    >>> counts = counter.count("ACGTACGT")
    >>> counts["ACG"]
    2

    Canonical mode:

    >>> counter = KmerCounter(k=3, canonical=True)
    >>> counts = counter.count("ACGT")
    >>> # ACG and CGT (reverse complements) counted together

    Multiple sequences:

    >>> counter = KmerCounter(k=4)
    >>> seqs = ["ACGTACGT", "TGCATGCA"]
    >>> counts = counter.count(seqs)
    >>> "ACGT" in counts
    True

    Notes
    -----
    K-mers containing ambiguous bases (N, R, Y, etc.) are automatically excluded
    from counting. Sequences are case-insensitive.

    See Also
    --------
    count_kmers : Convenience function for simple k-mer counting
    streaming_count_kmers : Memory-efficient counting from large files
    """

    def __init__(self, k, canonical=False, use_rust=True):
        """Initialize K-mer counter.

        Parameters
        ----------
        k : int
            K-mer length
        canonical : bool, optional
            Use canonical k-mers
        use_rust : bool, optional
            Use Rust if available
        """
        if k <= 0:
            raise ValueError(f"k must be positive, got {k}")

        self.k = k
        self.canonical = canonical
        self.use_rust = use_rust and RUST_AVAILABLE and k < 32
        self.using_rust = self.use_rust

        if k >= 32 and use_rust and RUST_AVAILABLE:
            warnings.warn(
                f"k={k} is too large for Rust 2-bit packing (max k=31). "
                "Falling back to Python implementation.",
                BiopythonWarning,
            )

    def count(self, sequences):
        """Count k-mers in sequences.

        Parameters
        ----------
        sequences : str, Seq, SeqRecord, or iterable
            Input sequence(s). Can be:
            - Single string
            - Bio.Seq.Seq object
            - Bio.SeqRecord.SeqRecord object
            - Iterable of any of the above

        Returns
        -------
        dict
            Dictionary mapping k-mer strings to their counts.

        Examples
        --------
        >>> counter = KmerCounter(k=3)
        >>> counter.count("ACGTACGT")  # doctest: +ELLIPSIS
        {...}

        >>> from Bio.Seq import Seq
        >>> counter.count(Seq("ACGT"))  # doctest: +ELLIPSIS
        {...}

        >>> counter.count(["ACGT", "TGCA"])  # doctest: +ELLIPSIS
        {...}

        Notes
        -----
        Empty sequences or sequences shorter than k return empty dictionaries.
        K-mers containing N or other ambiguous bases are excluded.
        """
        # Normalize input to list of strings
        seq_list = self._normalize_sequences(sequences)

        if self.using_rust:
            # Use stateful Rust KmerCounter class
            rust_counter = RustKmerCounter(self.k, self.canonical)
            rust_counter.update(seq_list)
            return rust_counter.get_counts()
        else:
            return self._count_python(seq_list)

    def _normalize_sequences(self, sequences):
        """Convert various input types to list of uppercase strings.

        Parameters
        ----------
        sequences : str, Seq, SeqRecord, or iterable
            Input sequences

        Returns
        -------
        list of str
            List of uppercase DNA sequences

        Raises
        ------
        TypeError
            If sequence type is not supported
        """
        # Handle single string
        if isinstance(sequences, str):
            return [sequences.upper()]

        # Handle Seq object (has upper method)
        if hasattr(sequences, "upper") and callable(getattr(sequences, "upper")):
            return [str(sequences).upper()]

        # Handle SeqRecord
        if hasattr(sequences, "seq"):
            return [str(sequences.seq).upper()]

        # Handle iterable of sequences
        result = []
        try:
            for seq in sequences:
                if isinstance(seq, str):
                    result.append(seq.upper())
                elif hasattr(seq, "seq"):  # SeqRecord
                    result.append(str(seq.seq).upper())
                elif hasattr(seq, "upper"):  # Seq
                    result.append(str(seq).upper())
                else:
                    raise TypeError(f"Unsupported sequence type: {type(seq)}")
        except TypeError:
            # sequences is not iterable
            raise TypeError(
                f"Sequences must be str, Seq, SeqRecord, or iterable of these. Got: {type(sequences)}"
            )

        return result

    def _count_python(self, seq_list):
        """Pure Python k-mer counting fallback.

        Parameters
        ----------
        seq_list : list of str
            List of uppercase DNA sequences

        Returns
        -------
        dict
            K-mer counts
        """
        from Bio.Seq import reverse_complement

        counts = Counter()

        for seq in seq_list:
            # Skip k-mers with ambiguous bases
            for i in range(len(seq) - self.k + 1):
                kmer = seq[i : i + self.k]

                # Skip if contains N or other ambiguous bases
                if "N" in kmer or not all(c in "ACGT" for c in kmer):
                    continue

                if self.canonical:
                    rc = str(reverse_complement(kmer))
                    kmer = min(kmer, rc)  # Lexicographically smaller

                counts[kmer] += 1

        return dict(counts)


def count_kmers(sequences, k, canonical=False):
    """Count k-mers in sequences (convenience function).

    This is a convenience wrapper around KmerCounter for simple use cases.

    Parameters
    ----------
    sequences : str, Seq, SeqRecord, or iterable
        Input sequence(s)
    k : int
        K-mer length
    canonical : bool, optional
        Use canonical k-mers (treat k-mer and reverse complement as same).
        Default is False.

    Returns
    -------
    dict
        K-mer counts

    Examples
    --------
    >>> counts = count_kmers("ACGTACGT", k=3)
    >>> counts["ACG"]
    2

    >>> counts = count_kmers("ACGTACGT", k=3, canonical=True)
    >>> # ACG and CGT counted together

    >>> from Bio.Seq import Seq
    >>> counts = count_kmers(Seq("ATCGATCG"), k=4)
    >>> "ATCG" in counts
    True

    See Also
    --------
    KmerCounter : Class interface with more control
    streaming_count_kmers : Memory-efficient counting from files
    """
    counter = KmerCounter(k=k, canonical=canonical)
    return counter.count(sequences)


def streaming_count_kmers(file_path, k, format="fasta", canonical=False):
    """Count k-mers from large sequence file using streaming I/O.

    Uses Bio.SeqIO to stream sequences from file, avoiding the need to load
    the entire file into memory. Processes sequences in batches for efficiency.

    Parameters
    ----------
    file_path : str
        Path to sequence file
    k : int
        K-mer length
    format : str, optional
        SeqIO format (default: "fasta"). Common formats include:
        "fasta", "fastq", "genbank", "embl", etc.
    canonical : bool, optional
        Use canonical k-mers (treat k-mer and reverse complement as same).
        Default is False.

    Returns
    -------
    dict
        K-mer counts

    Examples
    --------
    >>> counts = streaming_count_kmers("sequences.fasta", k=21)  # doctest: +SKIP

    >>> counts = streaming_count_kmers("reads.fastq", k=31, format="fastq")  # doctest: +SKIP

    >>> counts = streaming_count_kmers("genome.fasta", k=21, canonical=True)  # doctest: +SKIP

    Notes
    -----
    This function is designed for files too large to fit in memory. It processes
    sequences in batches of 1000 records at a time, accumulating k-mer counts
    incrementally.

    For small files that fit in memory, using count_kmers() with Bio.SeqIO.parse()
    may be more straightforward.

    See Also
    --------
    count_kmers : Simple k-mer counting for in-memory sequences
    KmerCounter : Class interface for more control
    Bio.SeqIO : Sequence file I/O
    """
    from Bio import SeqIO

    counter = KmerCounter(k=k, canonical=canonical)
    total_counts = Counter()

    # Process in batches for efficiency
    batch_size = 1000
    batch = []

    for record in SeqIO.parse(file_path, format):
        batch.append(record)

        if len(batch) >= batch_size:
            batch_counts = counter.count(batch)
            total_counts.update(batch_counts)
            batch = []

    # Process remaining records
    if batch:
        batch_counts = counter.count(batch)
        total_counts.update(batch_counts)

    return dict(total_counts)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
