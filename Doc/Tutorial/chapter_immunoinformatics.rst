.. _`chapter:immunoinformatics`:

Immunoinformatics and Adaptive Immune Receptors
================================================

This chapter covers Biopython's support for analyzing adaptive immune receptor
repertoire (AIRR) sequencing data. This includes:

- Reading and writing AIRR format files with ``Bio.SeqIO``
- Filtering and analyzing immune repertoire sequences
- Clonotype grouping and gene usage analysis with ``Bio.SeqUtils.Immuno``
- CDR3 sequence translation and distance metrics

The AIRR format is a standardized tab-separated values (TSV) format developed by
the AIRR Community for representing adaptive immune receptor repertoire sequencing
data from B-cell and T-cell receptors. For more information about the AIRR data
standards, see https://docs.airr-community.org/.

.. _`sec:airr-parsing`:

Parsing AIRR Files
------------------

AIRR files can be parsed using ``Bio.SeqIO.parse()`` just like any other sequence
format. The format string is ``"airr"``.

Basic AIRR File Reading
~~~~~~~~~~~~~~~~~~~~~~~~

Here's a simple example of reading an AIRR file:

.. code:: python

   from Bio import SeqIO

   # Parse AIRR file
   for record in SeqIO.parse("repertoire.tsv", "airr"):
       print(f"Sequence ID: {record.id}")
       print(f"Sequence length: {len(record.seq)} bp")

       # Access AIRR-specific annotations
       airr = record.annotations["airr"]
       print(f"V gene: {airr.get('v_call')}")
       print(f"J gene: {airr.get('j_call')}")
       print(f"Productive: {airr.get('productive')}")
       print()

The ``record.annotations["airr"]`` dictionary contains all AIRR-specific fields
from the file. Common fields include:

- ``v_call``, ``d_call``, ``j_call``: V, D, and J gene assignments
- ``junction``, ``junction_aa``: CDR3 junction sequence (nucleotide and amino acid)
- ``productive``: Boolean indicating if sequence is productive (in-frame, no stop codons)
- ``rev_comp``: Boolean indicating if sequence is reverse complement
- ``v_identity``, ``v_score``: V gene alignment identity and score

Type Conversion
~~~~~~~~~~~~~~~

AIRR fields are automatically converted to appropriate Python types:

.. code:: python

   from Bio import SeqIO

   record = next(SeqIO.parse("repertoire.tsv", "airr"))
   airr = record.annotations["airr"]

   # Boolean fields: T/F in file -> True/False in Python
   is_productive = airr["productive"]  # True or False
   print(f"Productive: {is_productive} (type: {type(is_productive).__name__})")

   # Integer fields: numeric strings -> int
   junction_len = airr["junction_length"]  # int
   v_start = airr["v_sequence_start"]  # int

   # Float fields: numeric strings -> float
   v_identity = airr["v_identity"]  # float (0.0 to 1.0)
   v_score = airr["v_score"]  # float

   # String fields: remain as strings
   v_gene = airr["v_call"]  # str (e.g., "IGHV1-69*01")

Empty or missing fields are represented as ``None``.

.. _`sec:airr-filtering`:

Filtering Productive Sequences
-------------------------------

A common first step in repertoire analysis is filtering for productive sequences.
Productive sequences are in-frame and lack stop codons, representing functional
immune receptors.

Filtering During Parsing
~~~~~~~~~~~~~~~~~~~~~~~~~

You can filter sequences as you parse them:

.. code:: python

   from Bio import SeqIO

   # Read only productive sequences
   productive_records = []

   for record in SeqIO.parse("repertoire.tsv", "airr"):
       airr = record.annotations["airr"]

       # Filter for productive sequences
       if airr.get("productive") is True:
           productive_records.append(record)

   print(f"Found {len(productive_records)} productive sequences")

You can combine multiple filtering criteria:

.. code:: python

   from Bio import SeqIO

   # More stringent filtering
   high_quality_records = []

   for record in SeqIO.parse("repertoire.tsv", "airr"):
       airr = record.annotations["airr"]

       # Filter for:
       # - Productive sequences
       # - High V gene identity (>95%)
       # - Complete VDJ rearrangement
       if (airr.get("productive") is True and
           airr.get("v_identity") and airr["v_identity"] > 0.95 and
           airr.get("complete_vdj") is True):
           high_quality_records.append(record)

   print(f"Found {len(high_quality_records)} high-quality productive sequences")

Writing Filtered Sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can write filtered sequences to a new AIRR file:

.. code:: python

   from Bio import SeqIO

   # Read and filter
   records = SeqIO.parse("repertoire.tsv", "airr")
   productive = (rec for rec in records
                 if rec.annotations["airr"].get("productive") is True)

   # Write to new file
   count = SeqIO.write(productive, "productive_only.tsv", "airr")
   print(f"Wrote {count} productive sequences")

.. _`sec:clonotype-grouping`:

Grouping Clonotypes
-------------------

Clonotypes are groups of B or T cell receptors that likely originated from the
same ancestor cell. The ``Bio.SeqUtils.Immuno`` module provides functions for
clonotype analysis.

Basic Clonotype Grouping
~~~~~~~~~~~~~~~~~~~~~~~~~

Clonotypes are typically defined by:

1. **V gene** - which V gene segment is used
2. **J gene** - which J gene segment is used
3. **CDR3 sequence** - the hypervariable junction region

.. code:: python

   from Bio import SeqIO
   from Bio.SeqUtils.Immuno import group_clonotypes

   # Parse AIRR file
   records = list(SeqIO.parse("repertoire.tsv", "airr"))

   # Group into clonotypes (nucleotide-based by default)
   clonotypes = group_clonotypes(records)

   print(f"Found {len(clonotypes)} unique clonotypes")

   # Examine clonotypes
   for (v_gene, j_gene, cdr3), sequences in clonotypes.items():
       if len(sequences) > 1:  # Expanded clonotypes
           print(f"Clonotype: {v_gene}/{j_gene}")
           print(f"  CDR3: {cdr3[:30]}...")  # First 30 bp
           print(f"  Size: {len(sequences)} sequences")

The default behavior uses nucleotide sequences (``junction`` field) for strict
clonotype identification. This distinguishes sequences with identical amino acid
translations but different nucleotide sequences due to codon degeneracy.

Nucleotide vs. Amino Acid Grouping
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can choose between two clonotyping strategies:

**Strict clonotyping** (default) uses nucleotide sequences:

.. code:: python

   from Bio import SeqIO
   from Bio.SeqUtils.Immuno import group_clonotypes

   records = list(SeqIO.parse("repertoire.tsv", "airr"))

   # Nucleotide-based: distinguishes synonymous codons
   nt_clonotypes = group_clonotypes(records, cdr3_field="junction")
   print(f"Nucleotide clonotypes: {len(nt_clonotypes)}")

**Functional grouping** uses amino acid sequences:

.. code:: python

   # Amino acid-based: groups convergent recombination
   aa_clonotypes = group_clonotypes(records, cdr3_field="junction_aa")
   print(f"Amino acid clonotypes: {len(aa_clonotypes)}")

The amino acid approach typically yields fewer clonotypes because multiple
nucleotide sequences can encode the same amino acid sequence:

.. code:: python

   # Compare the two approaches
   print(f"Nucleotide clonotypes: {len(nt_clonotypes)}")
   print(f"Amino acid clonotypes: {len(aa_clonotypes)}")
   print(f"Difference: {len(nt_clonotypes) - len(aa_clonotypes)}")

Finding Expanded Clonotypes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Expanded clonotypes (those with many sequences) often represent immune responses:

.. code:: python

   from Bio import SeqIO
   from Bio.SeqUtils.Immuno import group_clonotypes

   records = list(SeqIO.parse("repertoire.tsv", "airr"))
   clonotypes = group_clonotypes(records)

   # Find largest clonotypes
   sorted_clonotypes = sorted(clonotypes.items(),
                              key=lambda x: len(x[1]),
                              reverse=True)

   print("Top 10 expanded clonotypes:")
   for i, ((v_gene, j_gene, cdr3), sequences) in enumerate(sorted_clonotypes[:10], 1):
       print(f"{i}. {v_gene}/{j_gene}: {len(sequences)} sequences")

       # Show amino acid translation if available
       first_seq = sequences[0]
       cdr3_aa = first_seq.annotations["airr"].get("junction_aa")
       if cdr3_aa:
           print(f"   CDR3 (AA): {cdr3_aa}")

Allele Stripping
~~~~~~~~~~~~~~~~

By default, allele designations are stripped from gene names (e.g.,
``IGHV1-69*01`` becomes ``IGHV1-69``). You can disable this:

.. code:: python

   # Keep allele information
   clonotypes_with_alleles = group_clonotypes(records, strip_alleles=False)

   # Now "IGHV1-69*01" and "IGHV1-69*02" are separate groups

.. _`sec:gene-usage`:

Calculating Gene Usage Statistics
----------------------------------

The ``count_gene_usage()`` function tabulates how often each V, D, or J gene
is used in a repertoire.

V Gene Usage
~~~~~~~~~~~~

.. code:: python

   from Bio import SeqIO
   from Bio.SeqUtils.Immuno import count_gene_usage

   records = list(SeqIO.parse("repertoire.tsv", "airr"))

   # Count V gene usage
   v_usage = count_gene_usage(records, gene_field="v_call")

   print("Top 10 V genes:")
   for gene, count in list(v_usage.items())[:10]:
       print(f"{gene}: {count} sequences")

The results are automatically sorted by frequency (most common first).

J Gene Usage
~~~~~~~~~~~~

Similarly for J genes:

.. code:: python

   # Count J gene usage
   j_usage = count_gene_usage(records, gene_field="j_call")

   print("J gene usage:")
   for gene, count in j_usage.items():
       print(f"{gene}: {count}")

Calculating Frequencies
~~~~~~~~~~~~~~~~~~~~~~~

To get frequencies instead of counts:

.. code:: python

   from Bio import SeqIO
   from Bio.SeqUtils.Immuno import count_gene_usage

   records = list(SeqIO.parse("repertoire.tsv", "airr"))
   v_usage = count_gene_usage(records, gene_field="v_call")

   # Calculate total sequences
   total = sum(v_usage.values())

   print("V gene frequencies:")
   for gene, count in list(v_usage.items())[:10]:
       frequency = count / total * 100
       print(f"{gene}: {count} ({frequency:.1f}%)")

Filtering by Clonotype
~~~~~~~~~~~~~~~~~~~~~~~

You can combine gene usage with clonotype analysis:

.. code:: python

   from Bio import SeqIO
   from Bio.SeqUtils.Immuno import group_clonotypes, count_gene_usage

   records = list(SeqIO.parse("repertoire.tsv", "airr"))

   # Find expanded clonotypes (size >= 5)
   clonotypes = group_clonotypes(records)
   expanded_sequences = []

   for clonotype, sequences in clonotypes.items():
       if len(sequences) >= 5:
           expanded_sequences.extend(sequences)

   print(f"Found {len(expanded_sequences)} sequences in expanded clonotypes")

   # Gene usage within expanded clonotypes only
   v_usage_expanded = count_gene_usage(expanded_sequences, gene_field="v_call")

   print("V gene usage in expanded clonotypes:")
   for gene, count in list(v_usage_expanded.items())[:5]:
       print(f"{gene}: {count}")

.. _`sec:cdr3-analysis`:

CDR3 Sequence Analysis
----------------------

The CDR3 (Complementarity Determining Region 3) is the most variable part of
an immune receptor and determines antigen specificity.

Translating CDR3 Sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``translate_cdr3()`` function translates nucleotide CDR3 sequences to
amino acids with validation:

.. code:: python

   from Bio.SeqUtils.Immuno import translate_cdr3

   # Valid translation
   nt_seq = "TGTGCGAGAGATAGGCTC"
   aa_seq = translate_cdr3(nt_seq)
   print(f"Translation: {aa_seq}")  # CARDRL

   # Invalid sequences return None
   invalid = "TGTGCGAGAGATAG"  # Not divisible by 3
   result = translate_cdr3(invalid)
   print(f"Frameshift: {result}")  # None

   stop_codon = "TGTGCGAGATAGTAA"  # Contains stop codon
   result = translate_cdr3(stop_codon)
   print(f"Stop codon: {result}")  # None

This is useful for validating productive sequences:

.. code:: python

   from Bio import SeqIO
   from Bio.SeqUtils.Immuno import translate_cdr3

   validated_records = []

   for record in SeqIO.parse("repertoire.tsv", "airr"):
       airr = record.annotations["airr"]
       junction_nt = airr.get("junction")

       if junction_nt:
           # Validate translation
           junction_aa = translate_cdr3(junction_nt)

           if junction_aa is not None:
               # Valid productive sequence
               validated_records.append(record)

   print(f"Validated {len(validated_records)} productive sequences")

CDR3 Distance Metrics
~~~~~~~~~~~~~~~~~~~~~

The ``cdr3_distance()`` function calculates similarity between CDR3 sequences
using different distance metrics.

**Levenshtein distance** (default) counts insertions, deletions, and substitutions:

.. code:: python

   from Bio.SeqUtils.Immuno import cdr3_distance

   # Edit distance (default: Levenshtein)
   dist = cdr3_distance("CARDRLGTS", "CARDYYGTS")
   print(f"Edit distance: {dist}")  # 2 (R->Y, L->Y)

   # With different lengths
   dist = cdr3_distance("CARD", "CAR")
   print(f"Edit distance: {dist}")  # 1 (deletion of D)

**Hamming distance** counts only substitutions (sequences must be same length):

.. code:: python

   # Hamming distance (substitutions only)
   dist = cdr3_distance("CARD", "CARE", method="hamming")
   print(f"Hamming distance: {dist}")  # 1 (D->E)

   # Different lengths raise error
   try:
       dist = cdr3_distance("CARD", "CAR", method="hamming")
   except ValueError as e:
       print(f"Error: {e}")

Finding Similar CDR3 Sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can use distance metrics to find similar CDR3 sequences:

.. code:: python

   from Bio import SeqIO
   from Bio.SeqUtils.Immuno import cdr3_distance

   records = list(SeqIO.parse("repertoire.tsv", "airr"))

   # Pick a reference CDR3
   ref_record = records[0]
   ref_cdr3 = ref_record.annotations["airr"]["junction_aa"]

   # Find sequences within edit distance of 2
   similar_sequences = []

   for record in records[1:]:
       cdr3 = record.annotations["airr"].get("junction_aa")

       if cdr3 and len(cdr3) > 0:
           distance = cdr3_distance(ref_cdr3, cdr3)

           if distance <= 2:
               similar_sequences.append((record, distance))

   print(f"Found {len(similar_sequences)} similar sequences")
   for record, dist in similar_sequences[:5]:
       print(f"  {record.id}: distance {dist}")

.. _`sec:airr-writing`:

Writing AIRR Files
------------------

You can write sequence records to AIRR format using ``SeqIO.write()``.

Basic Writing
~~~~~~~~~~~~~

.. code:: python

   from Bio import SeqIO

   # Read, process, and write
   records = SeqIO.parse("input.tsv", "airr")

   # Filter for productive sequences
   productive = (rec for rec in records
                 if rec.annotations["airr"].get("productive") is True)

   # Write to new file
   count = SeqIO.write(productive, "output.tsv", "airr")
   print(f"Wrote {count} records")

Explicit Field Schema
~~~~~~~~~~~~~~~~~~~~~

For best results, specify which fields to include explicitly:

.. code:: python

   from Bio import SeqIO
   from Bio.SeqIO import AirrIO

   records = SeqIO.parse("input.tsv", "airr")

   # Specify fields explicitly
   fields = [
       "sequence_id", "sequence", "rev_comp", "productive",
       "v_call", "d_call", "j_call",
       "junction", "junction_aa", "junction_length",
       "v_identity", "v_score"
   ]

   with open("output.tsv", "w") as handle:
       writer = AirrIO.AirrWriter(handle, fields=fields)
       writer.write_file(records)

This ensures consistent output even if records have different fields.

.. _`sec:airr-random-access`:

Random Access and Indexing
---------------------------

For large AIRR files, you can use indexing for efficient random access:

.. code:: python

   from Bio import SeqIO

   # Build index (only reads file once)
   idx = SeqIO.index("large_repertoire.tsv", "airr")

   print(f"Indexed {len(idx)} sequences")

   # Fast random access by sequence_id
   record = idx["SEQ123456"]
   print(f"Retrieved: {record.id}")
   print(f"V gene: {record.annotations['airr']['v_call']}")

   # Can iterate through all records
   for seq_id in idx:
       # Access only needed records
       if seq_id.startswith("PRODUCTIVE"):
           record = idx[seq_id]
           # Process record...

   idx.close()

This is much more memory-efficient than loading all sequences at once.

.. _`sec:airr-complete-workflow`:

Complete Analysis Workflow
---------------------------

Here's a complete example that ties everything together:

.. code:: python

   from Bio import SeqIO
   from Bio.SeqUtils.Immuno import (group_clonotypes, count_gene_usage,
                                     translate_cdr3, cdr3_distance)

   # 1. Parse AIRR file and filter for productive sequences
   print("Step 1: Loading and filtering sequences...")
   all_records = list(SeqIO.parse("repertoire.tsv", "airr"))

   productive_records = [
       rec for rec in all_records
       if rec.annotations["airr"].get("productive") is True
   ]

   print(f"  Total sequences: {len(all_records)}")
   print(f"  Productive sequences: {len(productive_records)}")

   # 2. Group into clonotypes
   print("\nStep 2: Grouping clonotypes...")
   clonotypes = group_clonotypes(productive_records)
   print(f"  Unique clonotypes: {len(clonotypes)}")

   # Find expanded clonotypes
   expanded = [(key, seqs) for key, seqs in clonotypes.items()
               if len(seqs) >= 5]
   print(f"  Expanded clonotypes (>=5): {len(expanded)}")

   # Show top 3
   sorted_expanded = sorted(expanded, key=lambda x: len(x[1]), reverse=True)
   print("\n  Top 3 expanded clonotypes:")
   for i, ((v, j, cdr3), seqs) in enumerate(sorted_expanded[:3], 1):
       cdr3_aa = seqs[0].annotations["airr"].get("junction_aa", "N/A")
       print(f"    {i}. {v}/{j}: {len(seqs)} sequences")
       print(f"       CDR3 (AA): {cdr3_aa}")

   # 3. Calculate gene usage statistics
   print("\nStep 3: Gene usage analysis...")
   v_usage = count_gene_usage(productive_records, gene_field="v_call")
   j_usage = count_gene_usage(productive_records, gene_field="j_call")

   total_seqs = len(productive_records)
   print(f"  Top 5 V genes:")
   for gene, count in list(v_usage.items())[:5]:
       freq = count / total_seqs * 100
       print(f"    {gene}: {count} ({freq:.1f}%)")

   print(f"\n  Top 5 J genes:")
   for gene, count in list(j_usage.items())[:5]:
       freq = count / total_seqs * 100
       print(f"    {gene}: {count} ({freq:.1f}%)")

   # 4. CDR3 diversity analysis
   print("\nStep 4: CDR3 diversity...")

   # Get all CDR3 amino acid sequences
   cdr3_sequences = []
   for record in productive_records:
       junction_aa = record.annotations["airr"].get("junction_aa")
       if junction_aa:
           cdr3_sequences.append(junction_aa)

   # Calculate average CDR3 length
   avg_length = sum(len(seq) for seq in cdr3_sequences) / len(cdr3_sequences)
   print(f"  Average CDR3 length: {avg_length:.1f} amino acids")

   # Find most common CDR3 length
   from collections import Counter
   length_counts = Counter(len(seq) for seq in cdr3_sequences)
   most_common_length = length_counts.most_common(1)[0]
   print(f"  Most common length: {most_common_length[0]} AA ({most_common_length[1]} sequences)")

   # 5. Write filtered results
   print("\nStep 5: Writing output...")

   # Write expanded clonotypes to file
   expanded_sequences = []
   for (v, j, cdr3), seqs in expanded:
       expanded_sequences.extend(seqs)

   count = SeqIO.write(expanded_sequences, "expanded_clonotypes.tsv", "airr")
   print(f"  Wrote {count} sequences from expanded clonotypes")

   print("\nAnalysis complete!")

This workflow demonstrates the typical steps in immune repertoire analysis:
loading data, quality filtering, clonotype identification, gene usage analysis,
and result output.

.. _`sec:airr-summary`:

Summary
-------

This chapter covered:

- Parsing AIRR files with ``Bio.SeqIO.parse()``
- Filtering productive sequences using annotations
- Grouping clonotypes with ``group_clonotypes()``
- Calculating gene usage statistics with ``count_gene_usage()``
- CDR3 translation and distance metrics
- Writing AIRR files and random access indexing

For more information:

- AIRR Community standards: https://docs.airr-community.org/
- ``Bio.SeqIO`` documentation: Chapter :ref:`chapter:seqio`
- ``Bio.SeqUtils.Immuno`` API documentation: :py:mod:`Bio.SeqUtils.Immuno`
