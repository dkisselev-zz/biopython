# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is **Biopython**, a major open-source Python library for computational molecular biology and bioinformatics. Version 1.87.dev0, targeting Python 3.10+. The sole required runtime dependency is NumPy; many optional dependencies exist (matplotlib, networkx, reportlab, scipy, etc.).

## Local Development Environment

| Property | Value |
|---|---|
| OS | macOS (Darwin 24.6.0, Xcode CLI tools) |
| Architecture | x86_64 (Intel), 64-bit |
| CPU cores | 8 |
| RAM | 16 GB |
| Disk (free) | ~170 GB |
| Python | 3.12.10 (virtualenv at `.venv`) |
| C compiler | Apple clang 17.0.0 (via Xcode command-line tools) |
| GPU | None — PyTorch is installed but CUDA is not available |

### Installed Optional Dependencies

Available: `numpy 1.26.4`, `scipy 1.13.1`, `matplotlib 3.10.8`, `networkx 2.8.8`

Not installed: `reportlab`, `igraph`, `rdflib`, `psycopg2`, `mysqlclient`, `mmtf-python`

### Known Limitations on This Platform

- **No GPU/CUDA**: PyTorch is present but without CUDA. Any code path that requires GPU acceleration will fall back to CPU or be skipped.
- **No BioSQL database backends**: Neither PostgreSQL (`psycopg2`) nor MySQL connectors are installed. BioSQL-related tests will be skipped automatically.
- **No ReportLab**: `Bio.Graphics` tests requiring ReportLab will be skipped.
- **No igraph/rdflib**: Phylo tests that depend on these graph libraries will be skipped.
- **NumPy pinned at 1.26.x**: This is below NumPy 2.x. Code that uses NumPy 2-specific APIs (e.g., `np.strings`) will not work in this environment without an upgrade.
- **macOS x86_64 platform tag**: Wheel builds here target `macosx-10.13-x86_64`. If you need arm64 (Apple Silicon) artifacts, they must be built on a different machine or via CI.
- **Python 3.12 vs CI matrix**: CI tests 3.10–3.14 and PyPy. Local testing covers only 3.12. Syntax or runtime differences in 3.10 or 3.13/3.14 won't be caught locally.

## Key Commands

### Development Setup
```bash
pip install -e .                    # Editable install (builds C extensions)
pip install pre-commit && pre-commit install  # Set up pre-commit hooks
```

### Running Tests
```bash
# All tests (from repo root):
python setup.py test

# All tests, skip internet-dependent ones (preferred for local dev):
python setup.py test --offline

# Run a single test file:
cd Tests && python run_tests.py test_SeqIO

# Run with unittest directly (for specific test classes/methods):
cd Tests && python -m unittest test_SeqIO.TestSeqIO -v

# Run doctests only:
cd Tests && python run_tests.py doctest

# With coverage:
cd Tests && coverage run --source Bio,BioSQL run_tests.py --offline && coverage xml
```

### Linting & Style
```bash
# Run all pre-commit checks (black, ruff, flake8, mypy, etc.):
pre-commit run --all-files

# Individual tools:
black --check --target-version py310 Bio/
ruff check Bio/
flake8 Bio/
mypy Bio/
```

### Documentation
```bash
pip install -r .circleci/requirements-sphinx.txt
make -C Doc html     # Build HTML docs
make -C Doc latexpdf # Build PDF docs
```

## Architecture

### Source Layout
- `Bio/` — Main package (~50 subpackages). Each subpackage maps to a bioinformatics domain.
- `BioSQL/` — Database backend for storing/retrieving sequence data (PostgreSQL, MySQL, SQLite).
- `Tests/` — Test suite with 200+ `test_*.py` files and test data subdirectories.
- `Doc/` — Sphinx documentation with tutorials and API docs.
- `Scripts/` — Standalone utility scripts.

### Key Bio Subpackages (grouped by function)

**Core data types** — `Bio.Seq` (sequences), `Bio.SeqRecord` (annotated sequences), `Bio.SeqFeature` (feature annotations). Nearly everything else builds on these.

**I/O layer** — `Bio.SeqIO` (26+ sequence formats: FASTA, GenBank, EMBL, etc.), `Bio.AlignIO` (alignment formats), `Bio.SearchIO` (search result formats: BLAST, HMMER, etc.). These use a unified `read()`/`parse()`/`write()` interface pattern.

**Alignment & comparison** — `Bio.Align` (pairwise and multiple alignment), `Bio.Align.substitution_matrices` (BLOSUM, PAM), `Bio.codonalign`, `Bio.motifs`.

**Phylogenetics** — `Bio.Phylo` (tree I/O and manipulation), `Bio.Nexus` (Nexus format), `Bio.Phylo.PAML` (PAML wrapper).

**Structure** — `Bio.PDB` (Protein Data Bank, 3D structure handling — largest subpackage, ~47 modules). `Bio.SVDSuperimposer` for structure alignment.

**Database access** — `Bio.Entrez` (NCBI), `Bio.ExPASy`/`Bio.UniProt`/`Bio.SwissProt` (protein databases), `Bio.KEGG`, `Bio.Geo`.

**Analysis** — `Bio.Restriction` (restriction enzymes), `Bio.Cluster` (clustering with C extension), `Bio.HMM`, `Bio.PopGen`.

**Graphics** — `Bio.Graphics` (requires ReportLab), `Bio.Phylo` plotting (requires matplotlib).

### C Extensions
Performance-critical code is implemented as C extensions (compiled at install time). These live in source files alongside their Python counterparts (e.g., `Bio/Align/_Align.c`). The `setup.py` defines all Extension objects.

### Test Conventions
- Framework: `unittest` (not pytest).
- Test files: `Tests/test_*.py`.
- Test data: subdirectories under `Tests/` named by module (e.g., `Tests/GenBank/`, `Tests/PDB/`).
- Tests reference data files with relative paths assuming `Tests/` is the working directory.
- Internet-dependent tests check for connectivity and skip gracefully via `MissingExternalDependencyError`.
- Tests for missing optional packages (MySQL, etc.) are skipped, not failed.

### CI
- **GitHub Actions** (`.github/workflows/ci.yml`): style checks, builds (wheel + sdist), tests on Linux/macOS/Windows across Python 3.10–3.14 and PyPy. All CI tests run `--offline`.
- **AppVeyor**: Windows builds with MySQL/PostgreSQL for BioSQL testing.
- **CircleCI**: Documentation build and deployment.
- Pre-commit CI runs on changed files only; full checks are optional locally.

---

## Contribution Requirements

### Licensing
Every new file must include the standard copyright header and dual-license acknowledgment:

```python
# Copyright YYYY by Your Name.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
```

In your commit message or pull request body, explicitly state agreement to dual licensing under both the "Biopython License Agreement" and the "3-Clause BSD License".

### Code Style Rules (enforced by pre-commit)

| Tool | What it checks |
|---|---|
| **black** | Canonical formatting, target Python 3.10 |
| **ruff** | Linting (bugbear, comprehensions, docstring presence, import sorting, upgrade suggestions) |
| **flake8** | PEP 8 + pyflakes + RST docstring validity |
| **mypy** | Type correctness on `Bio/` and `BioSQL/` |

Key rules that are *not* enforced automatically but must be followed:

- **Imports**: one per line, sorted alphabetically (ruff isort with `force-single-line=true`, `order-by-type=false`). Example:
  ```python
  from typing import Optional
  from typing import Union

  from Bio import BiopythonWarning
  from Bio.Seq import Seq
  ```
- **Line length**: E501 is intentionally ignored — long lines are acceptable when they aid readability.
- **Module naming**: Existing modules use mixed case (e.g., `SeqIO`, `SeqRecord`). Match the existing convention for the package you're contributing to; do not rename modules to lowercase.
- **Internal symbols**: Mark private/internal functions and classes by appending `(PRIVATE)` to the first line of the docstring. Example: `"""Do internal bookkeeping (PRIVATE)."""`
- **Experimental code**: New code that is alpha/beta quality should emit `BiopythonExperimentalWarning` (defined in `Bio/__init__.py`). This signals to users that the API may change.
- **Deprecation**: When removing or replacing functionality, emit `BiopythonDeprecationWarning` first across at least one release cycle. Document the change in `DEPRECATED.rst`.
- **Type annotations**: Preferred on new public APIs, but the codebase is not fully typed. `Bio/PDB/` has relaxed mypy rules due to its size and complexity.
- **Optional dependency guards**: Wrap imports of optional packages in try/except and raise `MissingPythonDependencyError` with an informative message if the package is absent.

### Testing Requirements

- Every new feature or bug fix must include a test. The project policy (per CONTRIBUTING.rst) is that untested changes will not be accepted.
- Use `unittest.TestCase` subclasses — do not introduce pytest-specific fixtures or markers.
- Place test data files in the appropriate subdirectory under `Tests/` (e.g., `Tests/PDB/` for PDB test data). Reference them with paths relative to `Tests/`.
- Tests that require network access must be guarded so they are skipped when run with `--offline`.
- Tests that require missing optional packages must raise `MissingExternalDependencyError` or `MissingPythonDependencyError` to be reported as skipped, not failed.

---

## Documentation Style Requirements

### Docstrings

All public modules, classes, and functions must have docstrings. The project follows **PEP 257** with **reStructuredText (RST)** markup for rendering in Sphinx-generated API docs.

**Module docstrings** should provide an overview and link to the relevant wiki page and tutorial chapter:

```python
"""Short summary of what this module does.

See also the ModuleName_ wiki and the chapter in our Tutorial_.

.. _ModuleName: http://biopython.org/wiki/ModuleName
.. _Tutorial: https://biopython.org/docs/latest/Tutorial/index.html

"""
```

**Class and function docstrings** use a plain-text argument style (not Sphinx `:param:` directives), indented under the summary line:

```python
def translate(
    self, table="Standard", stop_symbol="*", to_stop=False, cds=False, gap="-"
):
    """Turn a nucleotide sequence into a protein sequence.

    This method will translate DNA or RNA sequences. It should not
    be used on protein sequences as any result will be biologically
    meaningless.

    Arguments:
     - table - Which codon table to use?  This can be either a name
       (string), an NCBI identifier (integer), or a CodonTable
       object (useful for non-standard genetic codes).  This
       defaults to the "Standard" table.
     - stop_symbol - What single character to use to represent any
       in-frame stop codons?  This defaults to an asterisk, "*".
     - to_stop - Boolean, defaults to False meaning do a full
       translation continuing on past any stop codons (and using
       the stop_symbol to represent these).  If True, the
       translation will terminate at the first in-frame stop codon
       (no stop symbol is added to the protein sequence).
     - cds - Boolean, defaults to False.  If True, ...
     - gap - Single character used to represent gaps in the input ...

    For example:

    >>> from Bio.Seq import Seq
    >>> coding_dna = Seq("ATGAAACCCGGGTTTTAA")
    >>> coding_dna.translate()
    Seq('MKPGF*')

    """
```

Key conventions:
- First line is a short imperative summary (fits on one line).
- Arguments listed under an `Arguments:` header with ` - name - description` format.
- Include inline `>>>` doctests where they illustrate usage concisely. These are run as part of the test suite (`python run_tests.py doctest`).
- Use RST inline markup: ``code``, *emphasis*, **strong**, and hyperlink targets where helpful.
- Avoid Sphinx-specific directives (`:param:`, `:type:`, `:returns:`, `:rtype:`) — the project does not use them.

### Tutorial and Documentation Files

- All `.rst` files are validated by `rstcheck` and `doc8` via pre-commit.
- Code blocks in RST docs are formatted by `blacken-docs`.
- DOI references must use the canonical `https://doi.org/` form (enforced by a pre-commit hook; `doi:`, `dx.doi.org`, and `http://doi.org` are rejected).
- Tutorial code examples in `Doc/Tutorial/` are tested via `Tests/test_Tutorial.py` — any code block you add there will be executed during the test run.
- Contributor names in `CONTRIB.rst` must be kept in alphabetical order (enforced by pre-commit).
