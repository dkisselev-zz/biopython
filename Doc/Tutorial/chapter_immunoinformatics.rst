.. _`chapter:immunoinformatics`:

Immunoinformatics with Biopython
=================================

This chapter covers immune receptor repertoire analysis using Biopython.
The AIRR (Adaptive Immune Receptor Repertoire) Community defines a standard
tab-separated values format for representing immunoglobulin and T-cell
receptor rearrangement data.  Biopython supports reading and writing AIRR
files through ``Bio.SeqIO`` and provides utilities for CDR3 analysis,
clonotype grouping, and gene usage tallying through ``Bio.SeqUtils.Immuno``.

.. _`sec:airr_io`:

Reading and Writing AIRR Files
-------------------------------

AIRR rearrangement TSV files have a header row followed by one rearrangement
per line.  The mandatory ``sequence_id`` and ``sequence`` columns map to
:attr:`Bio.SeqRecord.SeqRecord.id` and the sequence data; all other columns
are stored in :attr:`~Bio.SeqRecord.SeqRecord.annotations`.

.. doctest examples/immunoinformatics

.. code:: pycon

   >>> from Bio import SeqIO
   >>> records = list(SeqIO.parse("example.airr", "airr"))
   >>> len(records)
   5
   >>> records[0].id
   'heavy_seq1'
   >>> records[0].annotations["v_call"]
   'IGHV1-2*01'
   >>> records[0].annotations["productive"]
   True

Writing works symmetrically:

.. cont-doctest

.. code:: pycon

   >>> from io import StringIO
   >>> out = StringIO()
   >>> count = SeqIO.write(records, out, "airr")
   >>> count
   5

.. _`sec:airr_filtering`:

Filtering Productive Sequences
--------------------------------

The ``productive`` field indicates whether a rearrangement results in a
functional receptor.  Filtering on this flag is a common first step in
repertoire analysis:

.. cont-doctest

.. code:: pycon

   >>> productive = [r for r in records if r.annotations.get("productive")]
   >>> len(productive)
   4

.. _`sec:airr_translate`:

Translating Junction Sequences
--------------------------------

The junction (CDR3 region including anchoring residues) is stored as a
nucleotide string.  :func:`Bio.SeqUtils.Immuno.translate_junction` translates
it to amino acids, returning ``None`` for non-productive or out-of-frame
sequences:

.. doctest examples/immunoinformatics

.. code:: pycon

   >>> from Bio.SeqUtils.Immuno import translate_junction
   >>> translate_junction("TGTGCGAGAGATCGGCTG")
   'CARDRL'
   >>> translate_junction("TGTGCAAGACTGGACTAC")
   'CARLDY'

.. _`sec:airr_distance`:

CDR3 Distance Calculation
--------------------------

:func:`Bio.SeqUtils.Immuno.cdr3_distance` computes pairwise distance between
CDR3 sequences using Hamming (equal-length) or Levenshtein (variable-length)
metrics:

.. cont-doctest

.. code:: pycon

   >>> from Bio.SeqUtils.Immuno import cdr3_distance
   >>> cdr3_distance("CARDRL", "CARLDY", metric="hamming")
   3
   >>> cdr3_distance("CARDRL", "CARDRLL", metric="levenshtein")
   1

.. _`sec:airr_clonotypes`:

Grouping into Clonotypes
-------------------------

A clonotype is conventionally defined by a unique combination of V gene,
J gene, and CDR3 sequence.  :func:`Bio.SeqUtils.Immuno.group_clonotypes`
returns a dictionary keyed by ``(v_gene, j_gene, cdr3_aa)`` tuples:

.. doctest examples/immunoinformatics

.. code:: pycon

   >>> from Bio import SeqIO
   >>> from Bio.SeqUtils.Immuno import group_clonotypes
   >>> records = list(SeqIO.parse("example.airr", "airr"))
   >>> productive = [r for r in records if r.annotations.get("productive")]
   >>> clones = group_clonotypes(productive)
   >>> len(clones)
   4

.. _`sec:airr_gene_usage`:

Gene Usage Analysis
--------------------

Tallying V, D, or J gene usage is a standard repertoire characterization
step.  :func:`Bio.SeqUtils.Immuno.gene_usage` returns a
:class:`collections.Counter`:

.. cont-doctest

.. code:: pycon

   >>> from Bio.SeqUtils.Immuno import gene_usage
   >>> usage = gene_usage(productive, gene="v_call")
   >>> sorted(usage.items())
   [('IGHV1-2', 1), ('IGHV3-23', 1), ('IGKV1-5', 1), ('IGKV3-20', 1)]
