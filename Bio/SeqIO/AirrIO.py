# Copyright 2024 by Biopython Contributors. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SeqIO support for the AIRR (Adaptive Immune Receptor Repertoire) format.

You are expected to use this module via the Bio.SeqIO functions.

The AIRR format is a tab-separated values (TSV) file format developed by the
AIRR Community for representing adaptive immune receptor repertoire sequencing
data. The format includes a header line defining field names, followed by data
lines with tab-separated values.

For more information, see:
    https://docs.airr-community.org/

Examples
--------
Reading AIRR files:

>>> from Bio import SeqIO
>>> for record in SeqIO.parse("Tests/Airr/minimal.tsv", "airr"):
...     print(f"{record.id}: {len(record)} bp")
...     v_gene = record.annotations["airr"].get("v_call")
...     if v_gene:
...         print(f"  V gene: {v_gene}")

Writing AIRR files:

>>> records = SeqIO.parse("Tests/Airr/minimal.tsv", "airr")
>>> count = SeqIO.write(records, "output.tsv", "airr")
>>> print(f"Wrote {count} records")

Random access with indexing:

>>> idx = SeqIO.index("Tests/Airr/minimal.tsv", "airr")
>>> record = idx["SEQ001"]
>>> print(record.annotations["airr"].get("productive"))
>>> idx.close()

"""

import warnings

from Bio import BiopythonParserWarning
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .Interfaces import SequenceIterator
from .Interfaces import SequenceWriter

# AIRR standard field definitions (for type conversion)
# Based on AIRR Community standards v1.4+
BOOLEAN_FIELDS = {
    "rev_comp",
    "productive",
    "vj_in_frame",
    "stop_codon",
    "complete_vdj",
}

INTEGER_FIELDS = {
    "v_sequence_start",
    "v_sequence_end",
    "v_germline_start",
    "v_germline_end",
    "d_sequence_start",
    "d_sequence_end",
    "d_germline_start",
    "d_germline_end",
    "j_sequence_start",
    "j_sequence_end",
    "j_germline_start",
    "j_germline_end",
    "junction_length",
    "np1_length",
    "np2_length",
    "n1_length",
    "n2_length",
    "p3v_length",
    "p5d_length",
    "p3d_length",
    "p5j_length",
    "duplicate_count",
    "consensus_count",
}

FLOAT_FIELDS = {
    "v_score",
    "d_score",
    "j_score",
    "v_identity",
    "v_evalue",
    "d_identity",
    "d_evalue",
    "j_identity",
    "j_evalue",
}


def _convert_bool(value, field_name, record_id=None):
    """Convert string to boolean (T/F -> True/False).

    Parameters
    ----------
    value : str
        String value from TSV file
    field_name : str
        Name of the AIRR field (for error messages)
    record_id : str, optional
        Record identifier for better error messages

    Returns
    -------
    bool or None
        True, False, or None if empty/invalid

    """
    if not value or value.strip() == "":
        return None

    value = value.strip()
    if value == "T":
        return True
    elif value == "F":
        return False
    else:
        record_info = f" in record '{record_id}'" if record_id else ""
        warnings.warn(
            f"Invalid boolean value '{value}' for field '{field_name}'{record_info}. "
            f"Expected 'T' or 'F'. Value will be treated as None.",
            BiopythonParserWarning,
        )
        return None


def _convert_int(value, field_name, record_id=None):
    """Convert string to integer.

    Parameters
    ----------
    value : str
        String value from TSV file
    field_name : str
        Name of the AIRR field (for error messages)
    record_id : str, optional
        Record identifier for better error messages

    Returns
    -------
    int or None
        Integer value or None if empty/invalid

    """
    if not value or value.strip() == "":
        return None

    try:
        return int(value.strip())
    except ValueError:
        record_info = f" in record '{record_id}'" if record_id else ""
        warnings.warn(
            f"Cannot convert '{value}' to integer for field '{field_name}'{record_info}. "
            f"Value will be treated as None.",
            BiopythonParserWarning,
        )
        return None


def _convert_float(value, field_name, record_id=None):
    """Convert string to float.

    Parameters
    ----------
    value : str
        String value from TSV file
    field_name : str
        Name of the AIRR field (for error messages)
    record_id : str, optional
        Record identifier for better error messages

    Returns
    -------
    float or None
        Float value or None if empty/invalid

    """
    if not value or value.strip() == "":
        return None

    try:
        return float(value.strip())
    except ValueError:
        record_info = f" in record '{record_id}'" if record_id else ""
        warnings.warn(
            f"Cannot convert '{value}' to float for field '{field_name}'{record_info}. "
            f"Value will be treated as None.",
            BiopythonParserWarning,
        )
        return None


def _convert_str(value, field_name, record_id=None):
    """Convert string to string (passthrough with whitespace stripping).

    Parameters
    ----------
    value : str
        String value from TSV file
    field_name : str
        Name of the AIRR field (unused, for signature compatibility)
    record_id : str, optional
        Record identifier (unused, for signature compatibility)

    Returns
    -------
    str or None
        String value or None if empty

    """
    if not value or value.strip() == "":
        return None
    return value.strip()


def _convert_field_value(field_name, value, record_id=None):
    """Convert AIRR field value from string to appropriate Python type.

    Parameters
    ----------
    field_name : str
        Name of the AIRR field
    value : str
        String value from TSV file
    record_id : str, optional
        Record identifier for better error messages

    Returns
    -------
    bool, int, float, str, or None
        Converted value with appropriate type, or None if empty/missing

    Warnings
    --------
    Issues BiopythonParserWarning when type conversion fails for a non-empty value

    """
    # Handle empty/missing values
    if not value or value.strip() == "":
        return None

    value = value.strip()

    # Boolean fields: T/F -> True/False
    if field_name in BOOLEAN_FIELDS:
        if value == "T":
            return True
        elif value == "F":
            return False
        else:
            # Warn about invalid boolean value
            record_info = f" in record '{record_id}'" if record_id else ""
            warnings.warn(
                f"Invalid boolean value '{value}' for field '{field_name}'{record_info}. "
                f"Expected 'T' or 'F'. Value will be treated as None.",
                BiopythonParserWarning,
            )
            return None

    # Integer fields
    if field_name in INTEGER_FIELDS:
        try:
            return int(value)
        except ValueError:
            # Warn about invalid integer value
            record_info = f" in record '{record_id}'" if record_id else ""
            warnings.warn(
                f"Cannot convert '{value}' to integer for field '{field_name}'{record_info}. "
                f"Value will be treated as None.",
                BiopythonParserWarning,
            )
            return None

    # Float fields
    if field_name in FLOAT_FIELDS:
        try:
            return float(value)
        except ValueError:
            # Warn about invalid float value
            record_info = f" in record '{record_id}'" if record_id else ""
            warnings.warn(
                f"Cannot convert '{value}' to float for field '{field_name}'{record_info}. "
                f"Value will be treated as None.",
                BiopythonParserWarning,
            )
            return None

    # Default: return as string
    return value


def _format_field_value(value):
    """Convert Python value to AIRR format string.

    Parameters
    ----------
    value : bool, int, float, str, or None
        Python value to convert

    Returns
    -------
    str
        Formatted string for AIRR TSV output

    """
    if value is None:
        return ""
    if isinstance(value, bool):
        return "T" if value else "F"
    return str(value)


class AirrIterator(SequenceIterator):
    """Parser for AIRR (Adaptive Immune Receptor Repertoire) TSV files.

    The AIRR format is a tab-separated values file with:
    - Header line containing field names
    - Data lines with tab-separated values
    - Required fields: sequence_id, sequence

    Each record is returned as a SeqRecord with:
    - record.id and record.name set to sequence_id
    - record.seq set to the sequence
    - record.annotations["airr"] containing all AIRR fields as a dictionary
    - record.annotations["_airr_field_order"] preserving field order for round-trip
    - record.annotations["molecule_type"] set to "DNA"

    """

    modes = "t"

    def __init__(self, source):
        """Initialize AIRR format parser.

        Parameters
        ----------
        source : file-like object or path
            File opened in text mode, or path to file

        """
        super().__init__(source, fmt="AIRR")

        # Read and parse header line
        header_line = self.stream.readline()
        if not header_line:
            raise ValueError("Empty file or missing header line")

        self._field_names = header_line.rstrip("\n\r").split("\t")

        # Validate required fields
        if "sequence_id" not in self._field_names:
            raise ValueError("AIRR file must contain 'sequence_id' column")
        if "sequence" not in self._field_names:
            raise ValueError("AIRR file must contain 'sequence' column")

        # Pre-compute field indices for O(1) lookup
        self._field_indices = {name: i for i, name in enumerate(self._field_names)}
        self._seq_id_index = self._field_indices["sequence_id"]
        self._seq_index = self._field_indices["sequence"]

        # Pre-compute converter functions for each column (performance optimization)
        # This avoids repeated set lookups during parsing
        self._converters = []
        for field_name in self._field_names:
            if field_name in BOOLEAN_FIELDS:
                self._converters.append(_convert_bool)
            elif field_name in INTEGER_FIELDS:
                self._converters.append(_convert_int)
            elif field_name in FLOAT_FIELDS:
                self._converters.append(_convert_float)
            else:
                self._converters.append(_convert_str)

    def __next__(self):
        """Parse next record from AIRR file.

        Returns
        -------
        SeqRecord
            Next sequence record from file

        Raises
        ------
        StopIteration
            When end of file is reached
        ValueError
            When data line has significantly fewer columns than header,
            indicating a truncated or corrupt line

        """
        for line_num, line in enumerate(
            self.stream, start=2
        ):  # Start at 2 (after header)
            # Skip blank lines
            if not line.strip():
                continue

            # Split on tabs
            fields = line.rstrip("\n\r").split("\t")

            # Check for truncated lines (missing > 20% of columns)
            expected_columns = len(self._field_names)
            actual_columns = len(fields)
            missing_columns = expected_columns - actual_columns

            if missing_columns > 0:
                # Calculate percentage of missing columns
                missing_percentage = (missing_columns / expected_columns) * 100

                # Check if required fields are missing
                if (
                    actual_columns <= self._seq_id_index
                    or actual_columns <= self._seq_index
                ):
                    raise ValueError(
                        f"Truncated line at position {line_num}: "
                        f"line has {actual_columns} columns but required fields "
                        f"(sequence_id at column {self._seq_id_index + 1}, "
                        f"sequence at column {self._seq_index + 1}) are missing. "
                        f"Expected {expected_columns} columns."
                    )

                # Check if significantly truncated (> 20% missing)
                if missing_percentage > 20:
                    raise ValueError(
                        f"Truncated line at position {line_num}: "
                        f"line has {actual_columns} columns but header defines {expected_columns} "
                        f"({missing_columns} columns missing, {missing_percentage:.1f}%). "
                        f"This likely indicates a corrupt or truncated file."
                    )

                # If only moderately truncated (< 20%), warn and pad
                if missing_columns > 0:
                    warnings.warn(
                        f"Line {line_num} has {actual_columns} columns but header defines "
                        f"{expected_columns}. Missing {missing_columns} trailing columns "
                        f"will be treated as empty.",
                        BiopythonParserWarning,
                    )

            # Pad with empty strings if needed
            while len(fields) < expected_columns:
                fields.append("")

            # Extract sequence_id and sequence
            seq_id = fields[self._seq_id_index].strip()
            seq = fields[self._seq_index].strip()

            # Build annotations dictionary with type conversion
            # Use pre-calculated converters for performance
            airr_annotations = {}
            for i, field_name in enumerate(self._field_names):
                # Skip sequence_id and sequence (stored separately)
                if field_name in ("sequence_id", "sequence"):
                    continue

                # Get field value and convert using pre-calculated converter
                raw_value = fields[i] if i < len(fields) else ""
                converted_value = self._converters[i](raw_value, field_name, seq_id)
                airr_annotations[field_name] = converted_value

            # Create SeqRecord
            record = SeqRecord(
                Seq(seq),
                id=seq_id,
                name=seq_id,
                description="",
            )

            # Store AIRR annotations
            record.annotations["airr"] = airr_annotations
            record.annotations["_airr_field_order"] = self._field_names.copy()
            record.annotations["molecule_type"] = "DNA"

            return record

        raise StopIteration


class AirrWriter(SequenceWriter):
    """Writer for AIRR (Adaptive Immune Receptor Repertoire) TSV files.

    Writes SeqRecord objects to AIRR format TSV with:
    - Header line on first record
    - Tab-separated values
    - Type conversion (bool -> T/F, None -> empty string)

    **Field Schema Determination:**

    The output schema (column order and which fields to include) can be
    specified explicitly via the ``fields`` parameter in ``__init__``, or
    inferred from the first record written.

    **Recommended:** Always specify ``fields`` explicitly to ensure consistent
    output and prevent data loss. If fields are not specified, the writer will
    infer the schema from the first record, and any unique fields in subsequent
    records will be silently dropped (a warning is issued).

    Examples
    --------
    Explicit field specification (recommended):

    >>> from Bio import SeqIO  # doctest: +SKIP
    >>> fields = ["sequence_id", "sequence", "rev_comp", "productive", "v_call", "j_call"]  # doctest: +SKIP
    >>> with open("output.tsv", "w") as handle:  # doctest: +SKIP
    ...     writer = SeqIO.AirrIO.AirrWriter(handle, fields=fields)  # doctest: +SKIP
    ...     writer.write_file(records)  # noqa: F821  # doctest: +SKIP

    Automatic field inference (issues warning):

    >>> with open("output.tsv", "w") as handle:  # doctest: +SKIP
    ...     writer = SeqIO.AirrIO.AirrWriter(handle)  # Fields inferred  # doctest: +SKIP
    ...     writer.write_file(records)  # noqa: F821  # doctest: +SKIP

    """

    modes = "t"

    def __init__(self, target, fields=None):
        """Initialize AIRR format writer.

        Parameters
        ----------
        target : file-like object or path
            File opened in text mode, or path to file
        fields : list of str, optional
            Explicit list of field names to include in output, in desired order.
            Must include at minimum "sequence_id" and "sequence".
            If not provided, fields are inferred from the first record written,
            and a warning is issued about potential data loss.

        Raises
        ------
        ValueError
            If fields is provided but missing required fields "sequence_id" or "sequence"

        Warnings
        --------
        Issues BiopythonWarning if fields is not specified (schema inferred from
        first record, subsequent unique fields may be lost)

        Examples
        --------
        >>> from io import StringIO  # doctest: +SKIP
        >>> handle = StringIO()  # doctest: +SKIP
        >>> writer = AirrWriter(handle, fields=["sequence_id", "sequence", "v_call"])  # doctest: +SKIP
        >>> writer = AirrWriter(handle)  # Fields inferred, warning issued  # doctest: +SKIP

        """
        super().__init__(target)
        self._header_written = False
        self._schema_warning_issued = False

        # Validate and store explicit field order
        if fields is not None:
            if not isinstance(fields, (list, tuple)):
                raise TypeError(f"fields must be a list or tuple, got {type(fields)}")

            if "sequence_id" not in fields:
                raise ValueError(
                    "fields must include 'sequence_id' as a required field"
                )
            if "sequence" not in fields:
                raise ValueError("fields must include 'sequence' as a required field")

            self._field_order = list(fields)  # Make a copy
        else:
            # Fields will be inferred from first record
            self._field_order = None

    def write_record(self, record):
        """Write a single record to AIRR format.

        Parameters
        ----------
        record : SeqRecord
            Record to write

        Warnings
        --------
        Issues BiopythonWarning if:
        - Schema was not explicitly provided and is being inferred from first record
        - Current record has fields not in the schema (data will be lost)

        """
        # Determine field order from first record if not explicitly provided
        if not self._header_written:
            if self._field_order is None:
                # Schema not explicitly provided - infer from first record
                # Issue warning about potential data loss
                warnings.warn(
                    "AIRR output schema not explicitly specified. Inferring field "
                    "order from first record. Any unique fields in subsequent records "
                    "will be silently dropped. To avoid data loss, specify fields "
                    "explicitly: AirrWriter(handle, fields=[...]).",
                    BiopythonParserWarning,
                )
                self._schema_warning_issued = True

                # Try to get field order from record annotations
                if "_airr_field_order" in record.annotations:
                    self._field_order = record.annotations["_airr_field_order"]
                else:
                    # Default field order: sequence_id, sequence, then sorted others
                    airr_dict = record.annotations.get("airr", {})
                    other_fields = sorted(airr_dict.keys())
                    self._field_order = ["sequence_id", "sequence"] + other_fields

            # Write header
            self.handle.write("\t".join(self._field_order) + "\n")
            self._header_written = True

        # Check if current record has fields not in schema (potential data loss)
        if self._schema_warning_issued:
            airr_dict = record.annotations.get("airr", {})
            # Get fields in this record that are not in schema
            schema_fields = set(self._field_order) - {"sequence_id", "sequence"}
            record_fields = set(airr_dict.keys())
            missing_fields = record_fields - schema_fields

            if missing_fields:
                warnings.warn(
                    f"Record '{record.id}' has fields {sorted(missing_fields)} "
                    f"that are not in the output schema. These fields will be dropped. "
                    f"To include all fields, specify the complete schema explicitly.",
                    BiopythonParserWarning,
                )

        # Extract values for each field
        airr_dict = record.annotations.get("airr", {})
        values = []

        for field_name in self._field_order:
            if field_name == "sequence_id":
                value = record.id
            elif field_name == "sequence":
                value = str(record.seq)
            else:
                value = airr_dict.get(field_name)

            # Format value (convert types to strings)
            formatted_value = _format_field_value(value)
            values.append(formatted_value)

        # Write data line
        self.handle.write("\t".join(values) + "\n")


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest(verbose=0)
