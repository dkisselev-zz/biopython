# Copyright 2024 by Biopython contributors. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Fast distance matrix computation for phylogenetic analysis.

This module provides modern, high-performance implementations for computing
pairwise distance matrices from sequence alignments. Multiple computation
strategies are available with different performance characteristics.

Example
-------
>>> from Bio import Align
>>> from Bio.Phylo.DistanceMatrix import FastDistanceCalculator
>>> alignment = Align.read("example.phy", "phylip")
>>> calc = FastDistanceCalculator("identity", method="numpy")
>>> dm = calc.get_distance(alignment)

Available strategies:
    - ``python``: Pure Python loops (baseline, supports multiprocessing)
    - ``numpy``: NumPy vectorized operations (fast, memory efficient)
    - ``scipy``: SciPy pdist with custom metrics (requires scipy)
    - ``onehot``: One-hot encoding with matrix multiplication (very fast for identity)

"""

from collections.abc import Sequence
from typing import Union

import numpy as np

from Bio.Align import Alignment
from Bio.Align import MultipleSeqAlignment
from Bio.Align import substitution_matrices
from Bio.Phylo.TreeConstruction import DistanceMatrix

from Bio.Phylo.DistanceMatrix._base import condensed_to_squareform
from Bio.Phylo.DistanceMatrix._base import get_strategy
from Bio.Phylo.DistanceMatrix._base import list_strategies
from Bio.Phylo.DistanceMatrix._base import register_strategy
from Bio.Phylo.DistanceMatrix._parallel import get_n_jobs

# Import strategies to register them
from Bio.Phylo.DistanceMatrix import _strategies  # noqa: F401


class FastDistanceCalculator:
    """Modern distance calculator with pluggable strategies.

    This class provides the same interface as DistanceCalculator but with
    support for multiple computation backends and parallelism.

    :param model: Distance model name ('identity', 'blosum62', etc.).
    :param method: Computation strategy ('python', 'numpy', 'scipy', 'onehot').
    :param skip_letters: Characters to skip in distance calculation.
    :param n_jobs: Number of parallel workers (1=serial, -1=all CPUs).

    Examples
    --------
    >>> from Bio.Phylo.DistanceMatrix import FastDistanceCalculator
    >>> from Bio import Align
    >>> aln = Align.read("TreeConstruction/msa.phy", "phylip")
    >>> calc = FastDistanceCalculator("identity", method="numpy")
    >>> dm = calc.get_distance(aln)

    """

    # Class-level model lists (populated on first access)
    _dna_models: list[str] | None = None
    _protein_models: list[str] | None = None

    @classmethod
    def _load_models(cls) -> None:
        """Load available substitution matrix models."""
        if cls._dna_models is not None:
            return

        cls._dna_models = []
        cls._protein_models = []
        protein_alphabet = set("ACDEFGHIKLMNPQRSTVWY")

        names = substitution_matrices.load()
        for name in names:
            matrix = substitution_matrices.load(name)
            if name == "NUC.4.4":
                model_name = "blastn"
            else:
                model_name = name.lower()

            if protein_alphabet.issubset(set(matrix.alphabet)):
                cls._protein_models.append(model_name)
            else:
                cls._dna_models.append(model_name)

    @property
    def dna_models(self) -> list[str]:
        """List of available DNA distance models."""
        self._load_models()
        assert self._dna_models is not None
        return self._dna_models

    @property
    def protein_models(self) -> list[str]:
        """List of available protein distance models."""
        self._load_models()
        assert self._protein_models is not None
        return self._protein_models

    @property
    def models(self) -> list[str]:
        """List of all available distance models."""
        self._load_models()
        assert self._dna_models is not None
        assert self._protein_models is not None
        return ["identity"] + self._dna_models + self._protein_models

    def __init__(
        self,
        model: str = "identity",
        method: str = "numpy",
        skip_letters: Sequence[str] | None = None,
        n_jobs: int | None = None,
    ):
        """Initialize with a distance model and computation method."""
        self.model = model
        self.method = method
        self.n_jobs = get_n_jobs(n_jobs)

        # Set up skip_letters (same logic as original DistanceCalculator)
        if skip_letters is not None:
            self.skip_letters = tuple(skip_letters)
        elif model == "identity":
            self.skip_letters = ()
        else:
            self.skip_letters = ("-", "*")

        # Load scoring matrix
        if model == "identity":
            self.scoring_matrix: np.ndarray | None = None
            self.alphabet: str | None = None
        else:
            self._load_models()
            all_models = self.models
            if model not in all_models:
                raise ValueError(
                    f"Model not supported. Available models: {', '.join(all_models)}"
                )
            if model == "blastn":
                name = "NUC.4.4"
            else:
                name = model.upper()
            matrix = substitution_matrices.load(name)
            self.scoring_matrix = np.array(matrix, dtype=np.float64)
            self.alphabet = matrix.alphabet

        # Get strategy
        self._strategy = get_strategy(method)

    def get_distance(self, msa: Alignment | MultipleSeqAlignment) -> DistanceMatrix:
        """Return a DistanceMatrix for an Alignment or MultipleSeqAlignment.

        :param msa: Alignment or MultipleSeqAlignment object representing a
            DNA or protein multiple sequence alignment.
        :returns: DistanceMatrix with pairwise distances.
        """
        # Convert to numpy array
        if isinstance(msa, Alignment):
            names = [s.id for s in msa.sequences]
            alignment_array = np.array(msa)  # Uses __array__()
        elif isinstance(msa, MultipleSeqAlignment):
            names = [s.id for s in msa]
            # Vectorized conversion: concatenate all sequences into single buffer
            n_seqs = len(msa)
            seq_len = msa.get_alignment_length()
            # Join all sequences into one byte string and convert in single operation
            all_bytes = b"".join(str(s.seq).encode() for s in msa)
            alignment_array = np.frombuffer(all_bytes, dtype="S1").reshape(
                n_seqs, seq_len
            )
        else:
            raise TypeError(
                "Must provide an Alignment object or a MultipleSeqAlignment object."
            )

        # Compute distances
        condensed = self._strategy.compute(
            alignment_array,
            scoring_matrix=self.scoring_matrix,
            alphabet=self.alphabet,
            skip_letters=self.skip_letters,
            n_jobs=self.n_jobs,
        )

        # Convert condensed to DistanceMatrix
        return self._condensed_to_matrix(names, condensed)

    def _condensed_to_matrix(
        self, names: list, condensed: np.ndarray
    ) -> DistanceMatrix:
        """Convert condensed distance vector to DistanceMatrix."""
        n = len(names)

        # Create square matrix using NumPy indexing
        square = np.zeros((n, n), dtype=np.float64)
        tril_idx = np.tril_indices(n, k=-1)
        square[tril_idx] = condensed

        # Convert to list-of-lists (lower triangular with diagonal)
        # DistanceMatrix expects matrix[i] to have i+1 elements
        matrix = [square[i, : i + 1].tolist() for i in range(n)]

        return DistanceMatrix(names, matrix)


__all__ = [
    "FastDistanceCalculator",
    "get_strategy",
    "list_strategies",
    "register_strategy",
    "condensed_to_squareform",
]
