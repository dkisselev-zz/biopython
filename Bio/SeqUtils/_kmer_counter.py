# Copyright 2025 by the Biopython Contributors.
# All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""K-mer counting with an optional Rust acceleration layer.

This module provides :class:`KmerCounter`, a high-speed k-mer counter backed
by a Rust extension (PyO3 + Rayon) when available, falling back to a
pure-Python implementation otherwise.

Example usage::

    >>> from Bio.SeqUtils import KmerCounter
    >>> counter = KmerCounter(k=3)
    >>> counts = counter.count("ATCGATCG")
    >>> sorted(counts.items())
    [('ATC', 2), ('CGA', 1), ('GAT', 1), ('TCG', 2)]
"""

import contextlib
import itertools
import os
import re
from collections import Counter
from typing import Dict, Iterable, Iterator, Optional, Union

import Bio.SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ---------------------------------------------------------------------------
# Optional Rust backend
# ---------------------------------------------------------------------------
try:
    from Bio.SeqUtils._kmer_rust import count_kmers as _rust_count_kmers

    _USE_RUST = True
except ImportError:
    _USE_RUST = False

# ---------------------------------------------------------------------------
# Batch size for streaming file reads
# ---------------------------------------------------------------------------
_BATCH_SIZE = 10_000

# ---------------------------------------------------------------------------
# Pure-Python helpers
# ---------------------------------------------------------------------------
_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")
_NON_ACGT = re.compile(r"[^ACGT]+")


def _reverse_complement(seq: str) -> str:
    """Return the reverse complement of an ACGT string.

    :param seq: A DNA sequence containing only A, C, G, T (case-insensitive).
    :type seq: str
    :returns: The reverse complement.
    :rtype: str
    """
    return seq.translate(_COMPLEMENT)[::-1]


def _count_kmers_python(seq_bytes: bytes, k: int, use_canonical: bool) -> Dict[str, int]:
    """Count k-mers in *seq_bytes* using pure Python.

    The sequence is uppercased and split on any run of non-ACGT characters;
    k-mers are extracted from each resulting segment by string slicing.  A
    non-ACGT byte therefore acts as a window break — identical behaviour to
    the Rust path.

    :param seq_bytes: A single sequence as ASCII bytes.
    :type seq_bytes: bytes
    :param k: K-mer length.
    :type k: int
    :param use_canonical: If ``True``, replace each k-mer with the
        lexicographically smaller of itself and its reverse complement.
    :type use_canonical: bool
    :returns: Mapping from k-mer string to count.
    :rtype: dict[str, int]
    """
    counts: Counter = Counter()
    seq_str = seq_bytes.decode("ascii").upper()
    for segment in _NON_ACGT.split(seq_str):
        for i in range(len(segment) - k + 1):
            kmer = segment[i : i + k]
            if use_canonical:
                kmer = min(kmer, _reverse_complement(kmer))
            counts[kmer] += 1
    return dict(counts)


# ---------------------------------------------------------------------------
# Input normalisation
# ---------------------------------------------------------------------------
def _normalise_input(sequences: Union[str, Seq, SeqRecord, Iterable]) -> Iterator[bytes]:
    """Yield byte-strings from diverse input types, one per sequence.

    No intermediate list is materialised: each input element is converted and
    yielded immediately, so memory usage is proportional to the size of a
    single sequence rather than the entire input.

    Accepted types
    --------------
    * ``str``              — uppercased and ASCII-encoded.
    * :class:`Bio.Seq.Seq`           — converted via ``str()``, then uppercased and encoded.
    * :class:`Bio.SeqRecord.SeqRecord` — the ``.seq`` attribute is normalised.
    * An iterable of any of the above types.

    :param sequences: One or more nucleotide sequences.
    :type sequences: str | Seq | SeqRecord | Iterable
    :returns: An iterator of byte-strings ready for counting.
    :rtype: Iterator[bytes]
    :raises TypeError: If an element cannot be converted.
    """
    if isinstance(sequences, str):
        yield sequences.upper().encode("ascii")
        return
    if isinstance(sequences, Seq):
        yield str(sequences).upper().encode("ascii")
        return
    if isinstance(sequences, SeqRecord):
        yield from _normalise_input(sequences.seq)
        return
    # Iterable fallback — recurse on each element and flatten.
    for item in sequences:
        yield from _normalise_input(item)


# ---------------------------------------------------------------------------
# Public class
# ---------------------------------------------------------------------------
class KmerCounter:
    """High-speed k-mer counter with an optional Rust back-end.

    When the compiled Rust extension ``Bio.SeqUtils._kmer_rust`` is available,
    counting is performed in parallel using Rayon.  Otherwise a pure-Python
    fallback is used transparently.

    :param k: K-mer length.  Must satisfy ``1 <= k <= 32``.
    :type k: int
    :param canonical: If ``True``, each k-mer is canonicalised by taking the
        lexicographically smaller of itself and its reverse complement.  This
        collapses forward/reverse pairs into a single entry.
    :type canonical: bool
    :param threads: Number of worker threads for the Rust back-end.
        Defaults to :func:`os.cpu_count` when ``None``.  Ignored when the
        Python fallback is in use.
    :type threads: int | None

    :raises ValueError: If *k* is outside the range 1–32.

    Example::

        >>> from Bio.SeqUtils import KmerCounter
        >>> c = KmerCounter(k=3)
        >>> sorted(c.count("ATCGATCG").items())
        [('ATC', 2), ('CGA', 1), ('GAT', 1), ('TCG', 2)]
    """

    def __init__(self, k: int, canonical: bool = False, threads: Optional[int] = None):
        if not (1 <= k <= 32):
            raise ValueError("k must be between 1 and 32 inclusive")
        self.k = k
        self.canonical = canonical
        self.threads = threads if threads is not None else (os.cpu_count() or 1)

    def count(self, sequences: Union[str, Seq, SeqRecord, Iterable]) -> Dict[str, int]:
        """Count k-mers in one or more sequences.

        :param sequences: A single sequence (``str``, :class:`~Bio.Seq.Seq`,
            or :class:`~Bio.SeqRecord.SeqRecord`) or an iterable thereof.
        :type sequences: str | Seq | SeqRecord | Iterable
        :returns: Mapping from k-mer string to occurrence count.  K-mers are
            upper-case ACGT strings of length *k*.
        :rtype: dict[str, int]

        Example::

            >>> from Bio.SeqUtils import KmerCounter
            >>> c = KmerCounter(k=2)
            >>> sorted(c.count("ACGT").items())
            [('AC', 1), ('CG', 1), ('GT', 1)]
        """
        byte_iter = _normalise_input(sequences)
        merged: Dict[str, int] = {}
        if _USE_RUST:
            # Consume the generator in _BATCH_SIZE chunks.  Each chunk is
            # materialised only long enough to hand off to Rust; peak memory
            # is bounded to one chunk regardless of total input length.
            while True:
                batch = list(itertools.islice(byte_iter, _BATCH_SIZE))
                if not batch:
                    break
                for kmer, cnt in _rust_count_kmers(
                    batch, self.k, self.canonical, self.threads
                ).items():
                    merged[kmer] = merged.get(kmer, 0) + cnt
            return merged
        # Pure-Python path: count each sequence and merge.
        for seq_bytes in byte_iter:
            for kmer, cnt in _count_kmers_python(seq_bytes, self.k, self.canonical).items():
                merged[kmer] = merged.get(kmer, 0) + cnt
        return merged

    def count_file(
        self, handle, format: str = "fasta"  # noqa: A002
    ) -> Dict[str, int]:
        """Stream a sequence file and count k-mers without loading it all into memory.

        Records are read in batches of 10 000 via :func:`Bio.SeqIO.parse`.

        :param handle: An open file handle or a filename string.
        :type handle: IO | str
        :param format: A SeqIO-recognised format string (e.g. ``"fasta"``,
            ``"fastq"``).
        :type format: str
        :returns: Merged k-mer counts across all records in the file.
        :rtype: dict[str, int]

        Example::

            >>> import tempfile, os
            >>> from Bio.SeqUtils import KmerCounter
            >>> with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fh:
            ...     _ = fh.write(">seq1\\nATCGATCG\\n")
            ...     name = fh.name
            >>> c = KmerCounter(k=3)
            >>> sorted(c.count_file(name).items())
            [('ATC', 2), ('CGA', 1), ('GAT', 1), ('TCG', 2)]
            >>> os.unlink(name)
        """
        ctx = open(handle) if isinstance(handle, str) else contextlib.nullcontext(handle)
        with ctx as fh:
            merged: Dict[str, int] = {}
            batch: list = []
            for record in Bio.SeqIO.parse(fh, format):
                batch.append(record)
                if len(batch) >= _BATCH_SIZE:
                    for kmer, cnt in self.count(batch).items():
                        merged[kmer] = merged.get(kmer, 0) + cnt
                    batch = []
            if batch:
                for kmer, cnt in self.count(batch).items():
                    merged[kmer] = merged.get(kmer, 0) + cnt
            return merged
