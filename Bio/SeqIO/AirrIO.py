# Copyright 2024 by the Biopython contributors. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SeqIO support for AIRR Rearrangement TSV format.

You are expected to use this module via the Bio.SeqIO functions.

The AIRR (Adaptive Immune Receptor Repertoire) format is a tab-separated
values file for representing immune receptor rearrangement data as defined
by the AIRR Community (https://docs.airr-community.org/en/stable/).

Each row represents one rearrangement/read, with required columns
``sequence_id`` and ``sequence``.  Additional annotation columns (such as
``v_call``, ``j_call``, ``junction``, ``productive``, etc.) are stored in
:attr:`Bio.SeqRecord.SeqRecord.annotations`.

Boolean fields use ``"T"`` / ``"F"`` encoding.  Integer and float fields are
converted automatically.

Examples
--------
Reading an AIRR file:

>>> from Bio import SeqIO
>>> records = list(SeqIO.parse("AIRR/example.airr", "airr"))
>>> len(records)
3
>>> records[0].id
'seq1'
>>> records[0].annotations["productive"]
True

"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .Interfaces import _clean
from .Interfaces import SequenceIterator
from .Interfaces import SequenceWriter

_BOOL_FIELDS = frozenset(
    {
        "productive",
        "vj_in_frame",
        "stop_codon",
        "complete_vdj",
        "rev_comp",
    }
)

_INT_FIELDS = frozenset(
    {
        "duplicate_count",
        "v_sequence_start",
        "v_sequence_end",
        "d_sequence_start",
        "d_sequence_end",
        "j_sequence_start",
        "j_sequence_end",
        "cdr3_start",
        "cdr3_end",
        "v_germline_start",
        "v_germline_end",
        "d_germline_start",
        "d_germline_end",
        "j_germline_start",
        "j_germline_end",
        "np1_length",
        "np2_length",
        "v_score",
        "d_score",
        "j_score",
        "junction_length",
    }
)

_FLOAT_FIELDS = frozenset(
    {
        "v_identity",
        "d_identity",
        "j_identity",
        "v_support",
        "d_support",
        "j_support",
        "v_evalue",
        "d_evalue",
        "j_evalue",
    }
)

# Ordered superset of well-known AIRR Rearrangement fields used as the
# fallback schema when strict one-pass streaming is requested
# (AirrWriter buffer_size=0).  Follows the AIRR schema field order:
# gene calls → alignment strings → junction/CDR3 → quality flags →
# coordinates → statistical scores → counts.
_STANDARD_AIRR_FIELDS = (
    "locus",
    "v_call",
    "d_call",
    "j_call",
    "c_call",
    "sequence_alignment",
    "germline_alignment",
    "junction",
    "junction_aa",
    "np1",
    "np2",
    "productive",
    "vj_in_frame",
    "stop_codon",
    "complete_vdj",
    "rev_comp",
    "v_sequence_start",
    "v_sequence_end",
    "d_sequence_start",
    "d_sequence_end",
    "j_sequence_start",
    "j_sequence_end",
    "cdr3_start",
    "cdr3_end",
    "v_germline_start",
    "v_germline_end",
    "d_germline_start",
    "d_germline_end",
    "j_germline_start",
    "j_germline_end",
    "np1_length",
    "np2_length",
    "junction_length",
    "v_score",
    "d_score",
    "j_score",
    "v_identity",
    "d_identity",
    "j_identity",
    "v_support",
    "d_support",
    "j_support",
    "v_evalue",
    "d_evalue",
    "j_evalue",
    "duplicate_count",
)


def _parse_field(name, value):
    """Parse a single TSV field value according to the AIRR schema (PRIVATE).

    :param name: Column name.
    :type name: str
    :param value: Raw string value from the TSV cell.
    :type value: str
    :returns: Typed Python value.
    :rtype: bool | int | float | str | None
    :raises ValueError: If a boolean field contains an unexpected token,
        or an int/float field cannot be parsed.
    """
    if value == "":
        return None
    if name in _BOOL_FIELDS:
        if value == "T":
            return True
        elif value == "F":
            return False
        else:
            raise ValueError(
                f"Boolean field {name!r} expected 'T' or 'F', got {value!r}"
            )
    if name in _INT_FIELDS:
        try:
            return int(value)
        except ValueError:
            raise ValueError(
                f"Integer field {name!r} could not be parsed: {value!r}"
            ) from None
    if name in _FLOAT_FIELDS:
        try:
            return float(value)
        except ValueError:
            raise ValueError(
                f"Float field {name!r} could not be parsed: {value!r}"
            ) from None
    return value


def _format_field(name, value):
    """Format a single annotation value for TSV output (PRIVATE).

    :param name: Column name.
    :type name: str
    :param value: Python value to format.
    :returns: String suitable for writing to a TSV cell.
    :rtype: str
    """
    if value is None:
        return ""
    if name in _BOOL_FIELDS:
        return "T" if value else "F"
    return str(value)


class AirrIterator(SequenceIterator):
    """Parser for AIRR Rearrangement TSV files.

    Each line (after the header) is returned as a :class:`Bio.SeqRecord.SeqRecord`
    where the ``sequence`` column becomes the sequence data, ``sequence_id``
    becomes the record id/name, and all other columns are stored in
    :attr:`~Bio.SeqRecord.SeqRecord.annotations`.
    """

    modes = "t"

    def __init__(self, source):
        """Iterate over AIRR TSV records as SeqRecord objects.

        :param source: File-like object opened in text mode, or a path to a file.

        The file must have a header line with at least ``sequence_id`` and
        ``sequence`` columns separated by tabs.  Any blank lines in the body
        are silently skipped.  A UTF-8 BOM (``\\ufeff``) on the header line
        is stripped automatically.

        Examples
        --------
        >>> from Bio.SeqIO.AirrIO import AirrIterator
        >>> with open("AIRR/example.airr") as handle:
        ...     records = list(AirrIterator(handle))
        >>> len(records)
        3

        """
        super().__init__(source, fmt="AIRR")
        header = self.stream.readline()
        if not header:
            self._fields = None
            return
        # Strip UTF-8 BOM if present
        if header.startswith("\ufeff"):
            header = header[1:]
        fields = header.rstrip("\r\n").split("\t")
        if "sequence_id" not in fields:
            raise ValueError("AIRR file missing required 'sequence_id' column")
        if "sequence" not in fields:
            raise ValueError("AIRR file missing required 'sequence' column")
        self._fields = fields
        self._seq_id_idx = fields.index("sequence_id")
        self._seq_idx = fields.index("sequence")

    def __next__(self):
        """Return the next SeqRecord parsed from the AIRR TSV stream."""
        if self._fields is None:
            raise StopIteration
        for line in self.stream:
            if not line.strip():
                continue
            parts = line.rstrip("\r\n").split("\t")
            n_fields = len(self._fields)
            if len(parts) < n_fields:
                # Pad missing trailing fields with empty string (lenient)
                parts.extend([""] * (n_fields - len(parts)))
            elif len(parts) > n_fields:
                raise ValueError(
                    f"Data line has more fields ({len(parts)}) than header "
                    f"({n_fields}): {line!r}"
                )
            seq_id = parts[self._seq_id_idx]
            seq_str = parts[self._seq_idx]
            annotations = {}
            for i, field in enumerate(self._fields):
                if field in ("sequence_id", "sequence"):
                    continue
                annotations[field] = _parse_field(field, parts[i])
            return SeqRecord(
                Seq(seq_str),
                id=seq_id,
                name=seq_id,
                description="",
                annotations=annotations,
            )
        raise StopIteration


class AirrWriter(SequenceWriter):
    """Write AIRR Rearrangement TSV files.

    The column schema is determined by one of three strategies (in priority
    order):

    1. **Explicit schema** — pass a *fields* list to the constructor.
    2. **Buffer inference** — the first *buffer_size* records are inspected
       and the union of their annotation keys becomes the header
       (default, ``buffer_size=100``).  Records beyond the buffer window
       that contain keys absent from the inferred schema will have those
       values silently dropped.
    3. **Standard-field fallback** — when ``buffer_size=0``, the predefined
       :data:`_STANDARD_AIRR_FIELDS` superset is used so that all well-known
       AIRR columns always appear in the header even in a strict one-pass
       streaming context.

    This class is not intended to be used directly.  Instead, please use
    the top-level :func:`Bio.SeqIO.write` function with ``format="airr"``.
    """

    modes = "t"

    def __init__(self, target, fields=None, buffer_size=100):
        """Initialise the AIRR writer.

        :param target: File-like object opened in text mode, or a path.
        :param fields: Optional explicit list of column names (must start
            with ``sequence_id`` and ``sequence``).  When provided,
            *buffer_size* is ignored.
        :type fields: list[str] | None
        :param buffer_size: Number of records to buffer for schema inference
            when *fields* is not given.  The union of annotation keys across
            all buffered records becomes the header.  Set to ``0`` to disable
            buffering and fall back to :data:`_STANDARD_AIRR_FIELDS` instead
            (useful when strict one-pass streaming is required and the full
            key set is not known in advance).
        :type buffer_size: int
        """
        super().__init__(target)
        self._explicit_fields = fields
        self._fields = None
        self._buffer_size = buffer_size

    def write_records(self, records):
        """Write all records and return the count.

        The header is written before the first data row using the schema
        strategy selected at construction time.

        :param records: Iterable of :class:`Bio.SeqRecord.SeqRecord` objects.
        :returns: Number of records written.
        :rtype: int
        """
        records = iter(records)

        # Strategy 1: explicit field list supplied by caller
        if self._explicit_fields is not None:
            try:
                first = next(records)
            except StopIteration:
                return 0
            self._fields = self._explicit_fields
            self.handle.write("\t".join(self._fields) + "\n")
            self.write_record(first)
            count = 1
            for record in records:
                self.write_record(record)
                count += 1
            return count

        # Strategy 2: buffer a batch and infer the union of annotation keys
        if self._buffer_size > 0:
            buffer = []
            seen = set()
            key_order = []
            for record in records:
                buffer.append(record)
                for key in record.annotations:
                    if key not in seen:
                        seen.add(key)
                        key_order.append(key)
                if len(buffer) >= self._buffer_size:
                    break
            if not buffer:
                return 0
            self._fields = ["sequence_id", "sequence"] + key_order
            self.handle.write("\t".join(self._fields) + "\n")
            count = 0
            for record in buffer:
                self.write_record(record)
                count += 1
            for record in records:
                self.write_record(record)
                count += 1
            return count

        # Strategy 3: strict streaming — fall back to standard AIRR superset
        try:
            first = next(records)
        except StopIteration:
            return 0
        self._fields = ["sequence_id", "sequence"] + list(_STANDARD_AIRR_FIELDS)
        self.handle.write("\t".join(self._fields) + "\n")
        self.write_record(first)
        count = 1
        for record in records:
            self.write_record(record)
            count += 1
        return count

    def write_record(self, record):
        """Write a single AIRR data line.

        :param record: A :class:`Bio.SeqRecord.SeqRecord` object.
        """
        values = [_clean(record.id), str(record.seq)]
        for field in self._fields[2:]:
            raw = _format_field(field, record.annotations.get(field))
            # Sanitise: tabs would corrupt the TSV; newlines would create
            # spurious rows when the file is read back.
            raw = raw.replace("\t", " ").replace("\n", " ").replace("\r", " ")
            values.append(raw)
        self.handle.write("\t".join(values) + "\n")


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest(verbose=0)
