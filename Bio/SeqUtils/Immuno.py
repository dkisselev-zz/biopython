# Copyright 2024 by the Biopython contributors. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Utilities for immune repertoire (AIRR) data analysis.

Functions included:

1. translate_junction     - translate a nucleotide junction/CDR3 sequence
2. cdr3_distance          - compute Hamming or Levenshtein distance between CDR3s
3. group_clonotypes       - group SeqRecords by (V gene, J gene, CDR3)
4. gene_usage             - count V/D/J gene usage across records

Examples
--------
>>> from Bio.SeqUtils.Immuno import translate_junction
>>> translate_junction("TGTGCAAGAGATCGACTG")
'CARDRL'
>>> translate_junction("TGTGCAAGACTGGACTAC")
'CARLDY'

"""

from collections import Counter

from Bio.Seq import Seq


def translate_junction(junction_nt, productive=True):
    """Translate a nucleotide junction sequence to amino acids.

    :param junction_nt: Nucleotide junction/CDR3 sequence.
    :type junction_nt: str
    :param productive: If ``True`` (default), return ``None`` for sequences
        whose length is not a multiple of 3 or which contain an internal
        stop codon.  If ``False``, truncate to the nearest codon boundary
        and translate without stop-codon checking.
    :type productive: bool
    :returns: Amino acid string, or ``None`` if translation is not possible.
    :rtype: str | None

    Examples
    --------
    >>> translate_junction("TGTGCAAGAGATCGACTG")
    'CARDRL'
    >>> translate_junction("TGTGCAAGACTGGACTAC")
    'CARLDY'
    >>> translate_junction("TGTGCA") is None
    False
    >>> translate_junction("TGTG") is None
    True
    >>> translate_junction("TGTG", productive=False)
    'C'

    """
    if not junction_nt:
        return None
    if len(junction_nt) % 3 != 0:
        if productive:
            return None
        # Truncate to nearest codon boundary
        junction_nt = junction_nt[: (len(junction_nt) // 3) * 3]
        if not junction_nt:
            return None
    aa = str(Seq(junction_nt).translate())
    if productive and "*" in aa:
        return None
    return aa


def cdr3_distance(seq1, seq2, metric="hamming"):
    """Compute the distance between two CDR3 amino acid sequences.

    :param seq1: First CDR3 sequence.
    :type seq1: str
    :param seq2: Second CDR3 sequence.
    :type seq2: str
    :param metric: Distance metric to use.  Either ``"hamming"`` (requires
        equal-length sequences) or ``"levenshtein"`` (allows different
        lengths, pure-Python O(m*n) DP).
    :type metric: str
    :returns: Integer distance.
    :rtype: int
    :raises ValueError: If ``metric="hamming"`` and sequences have different
        lengths, or if an unknown metric is requested.

    Examples
    --------
    >>> cdr3_distance("CARDRL", "CARLDY")
    3
    >>> cdr3_distance("CARDRL", "CARDRLL", metric="levenshtein")
    1
    >>> cdr3_distance("CARDRL", "CARDRL")
    0

    """
    if metric == "hamming":
        if len(seq1) != len(seq2):
            raise ValueError(
                f"Hamming distance requires equal-length sequences, "
                f"got {len(seq1)} and {len(seq2)}"
            )
        return sum(a != b for a, b in zip(seq1, seq2))
    elif metric == "levenshtein":
        m, n = len(seq1), len(seq2)
        # Two-row DP
        prev = list(range(n + 1))
        for i in range(1, m + 1):
            curr = [i] + [0] * n
            for j in range(1, n + 1):
                if seq1[i - 1] == seq2[j - 1]:
                    curr[j] = prev[j - 1]
                else:
                    curr[j] = 1 + min(prev[j], curr[j - 1], prev[j - 1])
            prev = curr
        return prev[n]
    else:
        raise ValueError(f"Unknown metric {metric!r}; use 'hamming' or 'levenshtein'")


def group_clonotypes(
    records, v_field="v_call", j_field="j_call", cdr3_field="junction_aa"
):
    """Group immune receptor records into clonotypes.

    A clonotype is defined by a unique combination of V gene, J gene, and
    CDR3 amino acid sequence.  Multi-allele calls (e.g. ``"IGHV1-2*01,IGHV1-2*02"``)
    are normalised by taking the first entry.  Allele suffixes (``*01``) are
    stripped.

    :param records: Iterable of :class:`Bio.SeqRecord.SeqRecord` objects.
    :param v_field: Annotation key for the V gene call (default ``"v_call"``).
    :type v_field: str
    :param j_field: Annotation key for the J gene call (default ``"j_call"``).
    :type j_field: str
    :param cdr3_field: Annotation key for the CDR3/junction AA
        (default ``"junction_aa"``).
    :type cdr3_field: str
    :returns: Dictionary mapping ``(v_gene, j_gene, cdr3_aa)`` tuples to lists
        of SeqRecords.
    :rtype: dict

    Examples
    --------
    >>> from Bio.SeqUtils.Immuno import group_clonotypes
    >>> # (see tutorial for full example with AIRR data)

    """
    clonotypes = {}
    for record in records:
        v_raw = record.annotations.get(v_field) or ""
        j_raw = record.annotations.get(j_field) or ""
        cdr3 = record.annotations.get(cdr3_field) or ""
        # Normalise multi-allele: take first entry
        v_call = v_raw.split(",")[0].strip()
        j_call = j_raw.split(",")[0].strip()
        # Strip allele suffix
        v_gene = v_call.split("*")[0]
        j_gene = j_call.split("*")[0]
        key = (v_gene, j_gene, cdr3)
        if key in clonotypes:
            clonotypes[key].append(record)
        else:
            clonotypes[key] = [record]
    return clonotypes


def gene_usage(records, gene="v_call", strip_allele=True):
    """Tally V, D, or J gene usage across a set of records.

    Multi-allele calls (comma-separated) are split and each allele counted
    separately.  Empty or ``None`` calls are silently skipped.

    :param records: Iterable of :class:`Bio.SeqRecord.SeqRecord` objects.
    :param gene: Annotation key to tally (default ``"v_call"``).
    :type gene: str
    :param strip_allele: Strip ``*allele`` suffixes before counting
        (default ``True``).
    :type strip_allele: bool
    :returns: :class:`collections.Counter` mapping gene names to counts.
    :rtype: collections.Counter

    Examples
    --------
    >>> from Bio.SeqUtils.Immuno import gene_usage
    >>> # (see tutorial for full example with AIRR data)

    """
    counter = Counter()
    for record in records:
        call = record.annotations.get(gene)
        if not call:
            continue
        for entry in call.split(","):
            entry = entry.strip()
            if not entry:
                continue
            if strip_allele:
                entry = entry.split("*")[0]
            counter[entry] += 1
    return counter


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
