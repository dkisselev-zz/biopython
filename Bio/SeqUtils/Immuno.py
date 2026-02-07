# Copyright 2024 by Biopython Contributors. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Utilities for immune repertoire sequence analysis.

This module provides functions for analyzing adaptive immune receptor
repertoire (AIRR) data, including:

- CDR3 translation
- Sequence distance metrics (Hamming, Levenshtein)
- Clonotype grouping
- Gene usage counting

Examples
--------
Translate a CDR3 nucleotide sequence:

>>> from Bio.SeqUtils.Immuno import translate_cdr3
>>> translate_cdr3("TGTGCGAGAGATAGGCTC")
'CARDRL'
>>> translate_cdr3("TGTGCGAGAGATAG")  # Not divisible by 3
>>>

Efficient translation with SeqRecord:

>>> from Bio import SeqIO
>>> from Bio.SeqUtils.Immuno import translate_cdr3
>>> for record in SeqIO.parse("Tests/Airr/minimal.tsv", "airr"):
...     # Automatically extracts junction and checks productive status
...     aa_seq = translate_cdr3(record)
...     if aa_seq:
...         print(f"{record.id}: {aa_seq}")
>>>

Calculate CDR3 distance (unified interface):

>>> from Bio.SeqUtils.Immuno import cdr3_distance
>>> cdr3_distance("CARD", "CAR")  # Levenshtein (default)
1
>>> cdr3_distance("CARD", "CARE", method="hamming")  # Hamming distance
1
>>> cdr3_distance("CARDRLGTS", "CARDYYGTS")  # Edit distance
3

Calculate specific distance metrics:

>>> from Bio.SeqUtils.Immuno import hamming_distance, levenshtein_distance
>>> hamming_distance("CARD", "CARE")  # Substitutions only
1
>>> levenshtein_distance("CARD", "CAR")  # Insertions, deletions, substitutions
1

Group sequences by clonotype (nucleotide-based, default):

>>> from Bio import SeqIO
>>> from Bio.SeqUtils.Immuno import group_clonotypes
>>> records = list(SeqIO.parse("Tests/Airr/standard.tsv", "airr"))
>>> clonotypes = group_clonotypes(records)  # Uses nucleotide junction
>>> print(f"Found {len(clonotypes)} unique clonotypes")

Functional grouping by amino acid sequences:

>>> clonotypes_aa = group_clonotypes(records, cdr3_field="junction_aa")
>>> print(f"Found {len(clonotypes_aa)} functional groups")

Count V gene usage:

>>> from Bio.SeqUtils.Immuno import count_gene_usage
>>> v_usage = count_gene_usage(records, gene_field="v_call")
>>> for gene, count in list(v_usage.items())[:3]:
...     print(f"{gene}: {count}")

"""

import warnings
from collections import defaultdict

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def translate_cdr3(
    nucleotide_seq, codon_table=1, productive=None, junction_field="junction"
):
    """Translate nucleotide junction/CDR3 sequence to amino acids.

    This function can accept either a nucleotide sequence string or a SeqRecord
    with AIRR annotations. It performs validation and returns None for
    non-productive sequences (frameshifts, stop codons).

    Parameters
    ----------
    nucleotide_seq : str or SeqRecord
        Nucleotide sequence to translate. Can be:
        - str: Direct nucleotide sequence
        - SeqRecord: Record with AIRR annotations, will extract sequence from
          annotations["airr"][junction_field]
    codon_table : int, optional
        Codon table to use for translation (default: 1 = standard genetic code)
    productive : bool or None, optional
        If explicitly False, returns None immediately without translation.
        Useful for skipping known non-productive sequences from AIRR annotations.
        If None (default), performs translation and validation.
        If True, proceeds with translation (no special behavior).
    junction_field : str, optional
        Field name to extract junction sequence from SeqRecord AIRR annotations
        (default: "junction" for nucleotide sequences).
        Use "cdr3" for CDR3 sequences if available in your data.

    Returns
    -------
    str or None
        Amino acid sequence, or None if:
        - productive parameter is explicitly False
        - Length not divisible by 3 (frameshift)
        - Contains stop codons (non-productive)
        - Invalid nucleotides
        - SeqRecord missing junction field

    Examples
    --------
    Basic usage with nucleotide string:

    >>> translate_cdr3("TGTGCGAGAGATAGGCTC")
    'CARDRL'
    >>> translate_cdr3("TGTGCGAGAGATAG")  # Not divisible by 3
    >>> translate_cdr3("TGTGCGAGATAGTAA")  # Contains stop codon
    >>> translate_cdr3("")
    ''

    Skip translation if known to be non-productive:

    >>> translate_cdr3("TGTGCGAGATAGTAA", productive=False)

    With SeqRecord from AIRR file:

    >>> from Bio import SeqIO
    >>> record = next(SeqIO.parse("Tests/Airr/minimal.tsv", "airr"))
    >>> translate_cdr3(record)  # Extracts junction from annotations
    >>> # Or with explicit productive check
    >>> productive = record.annotations["airr"].get("productive")
    >>> translate_cdr3(record, productive=productive)

    """
    # Early return if explicitly marked as non-productive
    if productive is False:
        return None

    # Handle SeqRecord input
    if isinstance(nucleotide_seq, SeqRecord):
        # Extract junction sequence from AIRR annotations
        airr = nucleotide_seq.annotations.get("airr", {})

        # If productive not specified, try to get from record
        if productive is None:
            productive = airr.get("productive")
            # If record says non-productive, return None immediately
            if productive is False:
                return None

        # Extract junction sequence
        junction = airr.get(junction_field)
        if not junction:
            # No junction sequence available
            warnings.warn(
                f"SeqRecord '{nucleotide_seq.id}' has no '{junction_field}' "
                f"field in AIRR annotations",
                UserWarning,
            )
            return None

        nucleotide_seq = junction

    # Handle empty sequences
    if not nucleotide_seq:
        return ""

    # Check if length is divisible by 3 (frameshift check)
    if len(nucleotide_seq) % 3 != 0:
        return None

    try:
        # Translate sequence
        aa_seq = str(Seq(nucleotide_seq).translate(table=codon_table))

        # Check for stop codons (non-productive)
        if "*" in aa_seq:
            return None

        return aa_seq

    except Exception:
        # Invalid nucleotides or other translation errors
        return None


def hamming_distance(seq1, seq2):
    """Calculate Hamming distance between two sequences.

    Hamming distance is the number of positions at which the corresponding
    characters are different. Sequences must be the same length.

    Parameters
    ----------
    seq1 : str
        First sequence
    seq2 : str
        Second sequence

    Returns
    -------
    int
        Number of mismatched positions

    Raises
    ------
    ValueError
        If sequences have different lengths

    Examples
    --------
    >>> hamming_distance("CARD", "CARE")
    1
    >>> hamming_distance("CARD", "CARD")
    0
    >>> hamming_distance("CARD", "RARE")
    2

    """
    if len(seq1) != len(seq2):
        raise ValueError(
            f"Sequences must be same length (got {len(seq1)} and {len(seq2)})"
        )

    return sum(c1.upper() != c2.upper() for c1, c2 in zip(seq1, seq2))


def levenshtein_distance(seq1, seq2):
    """Calculate Levenshtein (edit) distance between two sequences.

    Levenshtein distance is the minimum number of single-character edits
    (insertions, deletions, or substitutions) required to change one
    sequence into the other.

    Uses dynamic programming with O(min(m,n)) space optimization.

    Parameters
    ----------
    seq1 : str
        First sequence
    seq2 : str
        Second sequence

    Returns
    -------
    int
        Edit distance between sequences

    Examples
    --------
    >>> levenshtein_distance("CARD", "CAR")
    1
    >>> levenshtein_distance("CARD", "CARE")
    1
    >>> levenshtein_distance("CARD", "HARD")
    1
    >>> levenshtein_distance("CARD", "CARD")
    0

    """
    # Handle empty sequences
    if not seq1:
        return len(seq2)
    if not seq2:
        return len(seq1)

    # Ensure seq1 is shorter for space efficiency
    if len(seq1) > len(seq2):
        seq1, seq2 = seq2, seq1

    # Convert to uppercase for case-insensitive comparison
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    # Initialize previous row (base case: distance from empty string)
    prev_row = list(range(len(seq2) + 1))

    # Dynamic programming: compute each row
    for i, c1 in enumerate(seq1, 1):
        curr_row = [i]  # First column: distance from empty string
        for j, c2 in enumerate(seq2, 1):
            if c1 == c2:
                # Characters match: no additional cost
                cost = prev_row[j - 1]
            else:
                # Characters don't match: min of insert, delete, substitute
                cost = 1 + min(
                    prev_row[j],  # Deletion
                    curr_row[j - 1],  # Insertion
                    prev_row[j - 1],  # Substitution
                )
            curr_row.append(cost)
        prev_row = curr_row

    return prev_row[-1]


def cdr3_distance(seq1, seq2, method="levenshtein"):
    """Calculate distance between two CDR3 sequences.

    Unified interface for computing distance between CDR3 amino acid or
    nucleotide sequences using different distance metrics.

    Parameters
    ----------
    seq1 : str
        First CDR3 sequence
    seq2 : str
        Second CDR3 sequence
    method : str, optional
        Distance metric to use:
        - "levenshtein" (default): Edit distance (insertions, deletions, substitutions)
        - "hamming": Hamming distance (substitutions only, requires same length)

    Returns
    -------
    int
        Distance between sequences

    Raises
    ------
    ValueError
        If method is "hamming" and sequences have different lengths, or
        if method is not recognized

    Examples
    --------
    >>> from Bio.SeqUtils.Immuno import cdr3_distance

    Levenshtein distance (default):

    >>> cdr3_distance("CARD", "CAR")
    1
    >>> cdr3_distance("CARD", "CARE")
    1
    >>> cdr3_distance("CARDRLGTS", "CARDYYGTS")
    3

    Hamming distance:

    >>> cdr3_distance("CARD", "CARE", method="hamming")
    1
    >>> cdr3_distance("CARD", "CARD", method="hamming")
    0

    Case insensitive:

    >>> cdr3_distance("CARD", "card")
    0

    """
    method = method.lower()

    if method == "hamming":
        return hamming_distance(seq1, seq2)
    elif method == "levenshtein":
        return levenshtein_distance(seq1, seq2)
    else:
        raise ValueError(
            f"Unknown distance method: '{method}'. "
            "Supported methods: 'hamming', 'levenshtein'"
        )


def group_clonotypes(
    records,
    v_field="v_call",
    j_field="j_call",
    cdr3_field="junction",
    strip_alleles=True,
):
    """Group sequences into clonotypes by V/J genes and CDR3.

    Clonotypes are groups of B or T cell receptors that likely originated
    from the same ancestor cell. Sequences are grouped by V gene, J gene,
    and CDR3 sequence.

    **Clonotyping Strategy:**

    - **Strict clonotyping** (default): Uses nucleotide sequences (``junction``)
      for precise clonotype definition. This avoids grouping sequences with
      identical amino acid translations but different nucleotide sequences
      (due to codon degeneracy).

    - **Functional grouping**: Use amino acid sequences (``junction_aa``)
      by setting ``cdr3_field="junction_aa"``. This groups sequences with
      the same functional CDR3 region regardless of silent mutations, useful
      for identifying convergent recombination.

    The default nucleotide-based approach is recommended for strict clonotype
    identification in repertoire analysis, as it maintains the highest resolution
    and accurately reflects clonal relationships.

    Parameters
    ----------
    records : iterable of SeqRecord
        Sequence records with AIRR annotations
    v_field : str, optional
        Name of V gene field in annotations["airr"] (default: "v_call")
    j_field : str, optional
        Name of J gene field in annotations["airr"] (default: "j_call")
    cdr3_field : str, optional
        Name of CDR3 field in annotations["airr"].
        Use "junction" for nucleotide-based strict clonotyping (default),
        or "junction_aa" for amino acid-based functional grouping.
        Can also use "cdr3" or "cdr3_aa" depending on data availability.
    strip_alleles : bool, optional
        If True, strip allele designations (IGHV1-69*01 → IGHV1-69)
        (default: True)

    Returns
    -------
    dict
        Dictionary mapping (v_gene, j_gene, cdr3_sequence) tuples to lists of
        SeqRecord objects. Keys with None values indicate missing data.

    Examples
    --------
    Strict clonotyping with nucleotide sequences (default):

    >>> from Bio import SeqIO
    >>> from Bio.SeqUtils.Immuno import group_clonotypes
    >>> records = list(SeqIO.parse("Tests/Airr/standard.tsv", "airr"))
    >>> clonotypes = group_clonotypes(records)  # Uses nucleotide junction
    >>> print(f"Found {len(clonotypes)} unique clonotypes")

    Functional grouping with amino acid sequences:

    >>> clonotypes_aa = group_clonotypes(records, cdr3_field="junction_aa")
    >>> print(f"Found {len(clonotypes_aa)} functional groups")
    >>> # May be fewer groups due to codon degeneracy

    Find largest clonotype:

    >>> largest = max(clonotypes.items(), key=lambda x: len(x[1]))
    >>> v_gene, j_gene, cdr3 = largest[0]
    >>> print(f"Largest clonotype: {v_gene}/{j_gene}/{cdr3[:20]}... ({len(largest[1])} sequences)")

    """
    clonotypes = defaultdict(list)

    for record in records:
        airr = record.annotations.get("airr", {})

        # Extract gene calls and CDR3
        v_gene = airr.get(v_field)
        j_gene = airr.get(j_field)
        cdr3 = airr.get(cdr3_field)

        # Strip allele designations if requested
        if strip_alleles:
            if v_gene and "*" in v_gene:
                v_gene = v_gene.split("*")[0]
            if j_gene and "*" in j_gene:
                j_gene = j_gene.split("*")[0]

        # Group by (V, J, CDR3) tuple
        clonotypes[(v_gene, j_gene, cdr3)].append(record)

    return dict(clonotypes)


def count_gene_usage(records, gene_field="v_call", strip_alleles=True):
    """Count V/D/J gene frequencies in a set of sequences.

    Parameters
    ----------
    records : iterable of SeqRecord
        Sequence records with AIRR annotations
    gene_field : str, optional
        Name of gene field in annotations["airr"]
        (e.g., "v_call", "d_call", "j_call") (default: "v_call")
    strip_alleles : bool, optional
        If True, strip allele designations (IGHV1-69*01 → IGHV1-69)
        (default: True)

    Returns
    -------
    dict
        Dictionary mapping gene names to counts, sorted by count (descending).
        None key indicates sequences with missing gene assignments.

    Examples
    --------
    >>> from Bio import SeqIO
    >>> from Bio.SeqUtils.Immuno import count_gene_usage
    >>> records = list(SeqIO.parse("Tests/Airr/standard.tsv", "airr"))
    >>> # Count V gene usage
    >>> v_usage = count_gene_usage(records, gene_field="v_call")
    >>> for gene, count in list(v_usage.items())[:5]:
    ...     print(f"{gene}: {count}")
    >>> # Count J gene usage
    >>> j_usage = count_gene_usage(records, gene_field="j_call")

    """
    gene_counts = defaultdict(int)

    for record in records:
        airr = record.annotations.get("airr", {})
        gene = airr.get(gene_field)

        # Strip allele designation if requested
        if strip_alleles and gene and "*" in gene:
            gene = gene.split("*")[0]

        gene_counts[gene] += 1

    # Sort by count (descending)
    return dict(sorted(gene_counts.items(), key=lambda x: x[1], reverse=True))


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest(verbose=0)
