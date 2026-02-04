# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is **Biopython** (v1.87.dev0), a major Python library for computational molecular biology. It contains 69 packages under `Bio/` and `BioSQL/`, with 12 C extensions for performance-critical paths. Only dependency is NumPy; many optional dependencies exist for specific submodules.

> **CRITICAL — run before any test or Python code in this repo:**
> A sibling copy of this repository exists at `../model_b` (same `Bio/` package structure).
> Without an editable install the Python import system can silently resolve `Bio` to the
> wrong tree, producing confusing failures or quietly testing stale code.
> **Always** run the install step below first in every session, and always run tests
> from within this directory.
> ```bash
> pip install -e .
> ```

## Commands

### Install (development mode — must run first every session)
```bash
pip install -e .
```
Compiles the C extensions and pins the `Bio` package to **this** directory in the
`mercor_3` virtualenv. Required before running any tests or importing `Bio`.
Because `../model_b` is an identical repo, skipping this step risks importing from
the wrong tree with no error.

### Run all tests
```bash
python setup.py test
python setup.py test --offline   # skip tests requiring internet
```
Alternatively, run directly from the Tests directory:
```bash
cd Tests && python run_tests.py
cd Tests && python run_tests.py --offline
```

### Run a single test
The runner must be invoked **from the `Tests/` directory** and the name must match an actual `test_*.py` file (without the `.py` extension). There is no `test_Seq.py`; the sequence-object tests live in `test_Seq_objs.py`.
```bash
cd Tests && python run_tests.py test_SeqUtils       # .py extension optional
cd Tests && python run_tests.py test_SeqIO
cd Tests && python run_tests.py test_Seq_objs
```

### Run doctests only
```bash
cd Tests && python run_tests.py doctest
```

### Linting / style checks
`pre-commit` is **not installed** in the active environment by default. Install it first:
```bash
pip install pre-commit
pre-commit install          # activate git hooks
pre-commit run --all-files  # run all hooks manually
```
The pre-commit config runs: **black** (formatting), **ruff** + **flake8** (linting), **mypy** (type checking on Bio/ and BioSQL/), **rstcheck** + **doc8** (RST docs), **codespell** (spelling).

### Build docs
```bash
cd Doc && make html
```
Requires sphinx, sphinx_rtd_theme, numpydoc, and biopython installed.

## Architecture

### Package layout
- **`Bio/`** — all Biopython modules (69 packages). The main namespace.
- **`BioSQL/`** — database persistence layer for biological sequences (PostgreSQL, MySQL).
- **`Tests/`** — unittest-based test suite with a custom runner (`run_tests.py`). Tests are `test_*.py` files; test data lives alongside them in the same directory.
- **`Doc/`** — Sphinx documentation, including a Tutorial with per-chapter subdirectories.
- **`Scripts/`** — standalone utility scripts.

### C extensions
Defined in `setup.py` under `EXTENSIONS`. Performance-critical modules include:
- `Bio.Align._aligncore`, `_pairwisealigner`, `_codonaligner` — sequence alignment internals
- `Bio.cpairwise2` — pairwise alignment (legacy C module)
- `Bio.PDB.kdtrees`, `Bio.PDB.ccealign`, `Bio.PDB._bcif_helper` — PDB structure analysis
- `Bio.Cluster._cluster` — clustering algorithms
- `Bio.SeqIO._twoBitIO` — 2bit genome format reader
- `Bio.Nexus.cnexus`, `Bio.motifs._pwm` — Nexus parsing, motif scoring

### Key architectural patterns
- **Optional-dependency modules**: Many submodules (Bio.Graphics, Bio.Cluster, Bio.Phylo plotting, etc.) gracefully raise `Bio.MissingPythonDependencyError` when their optional dependency is absent. Tests use this to skip appropriately.
- **Experimental code**: Modules under active development use `Bio.BiopythonExperimentalWarning` to signal instability.
- **Deprecated module pattern**: `Bio.Alphabet` was fully removed; importing it raises an informative `ImportError` guiding users to the replacement.
- **Version is in `Bio/__init__.py`** as `__version__`.

### Major submodule groups
| Area | Key modules |
|------|-------------|
| Sequences & I/O | `Bio.Seq`, `Bio.SeqIO`, `Bio.SeqRecord`, `Bio.SeqUtils` |
| Alignment | `Bio.Align`, `Bio.AlignIO`, `Bio.codonalign`, `Bio.pairwise2` |
| Structure (PDB) | `Bio.PDB` (large sub-package with mmCIF, mmtf, kdtrees) |
| NCBI interfaces | `Bio.Entrez`, `Bio.Blast`, `Bio.SearchIO` |
| Phylogenetics | `Bio.Phylo`, `Bio.Nexus` |
| Motifs & restriction | `Bio.Motifs`, `Bio.Restriction` |
| ML / stats | `Bio.Cluster`, `Bio.kNN`, `Bio.MarkovModel`, `Bio.NaiveBayes` |
| Databases | `BioSQL` |

## Testing conventions
- The test runner is **not pytest** — it is a custom `unittest`-based runner in `Tests/run_tests.py`.
- Tests that need the network are tagged; use `--offline` to exclude them.
- Tests that need optional dependencies (NumPy, matplotlib, etc.) skip themselves via `MissingPythonDependencyError` or `unittest.skipIf`.
- CI runs with `AddressSanitizer` (`-fsanitize=address,undefined`) on Linux to catch memory bugs in the C extensions.

## Style rules to follow
- **Black** for formatting (target Python 3.10+). Run before committing.
- **Existing module names are not lower_case** — do not rename them; this is a known, intentional exception to PEP 8.
- **Docstrings use reStructuredText (RST)** markup, not Google or NumPy style.
- `.flake8` and `ruff` configs are in the repo root; both ignore `E501` (line length is Black's job).
- PRs must include dual-license agreement in the commit message (Biopython License + 3-clause BSD).

## Local Development Environment

### Hardware & OS
| Property | Value |
|----------|-------|
| Machine | MacBookPro15,2 |
| CPU | Quad-Core Intel Core i7 @ 2.7 GHz (x86_64) |
| Cores | 8 (logical) |
| RAM | 16 GB |
| GPU | None (Intel integrated graphics only — no discrete GPU, no Metal compute, no CUDA) |
| OS | macOS Sequoia (Darwin 24.6.0) |

### Python environment
Managed by **pyenv**. The active environment is `mercor_3`, a virtualenv built on **Python 3.12.6**. Multiple other pyenv envs exist (`jhu`, `jhu-new`, `mercor_2`) — do not switch into them; stay on `mercor_3`.

`setuptools` is **not** pre-installed in `mercor_3`. The test runner (`Tests/run_tests.py`) imports it at the top level, so it must be present:
```bash
pip install setuptools      # required before running tests
```

**Model-A / Model-B collision risk.** Both `model_a/` and `model_b/` are checked-out copies of the same Biopython repo sharing the same virtualenv. `pip install -e .` writes a `.egg-link` / editable `.pth` entry that points the interpreter to whichever copy was installed last. If you run code without reinstalling first, `import Bio` may silently resolve to `model_b`. The safe sequence at the start of every session is:
```bash
pip install -e .            # pins Bio → model_a in the active env
```
Run this before every `python` or `run_tests.py` invocation.

### Installed optional dependencies (verified working)
| Package | Version | Used by |
|---------|---------|---------|
| numpy | 2.4.2 | Core requirement; almost all of Bio/ |
| scipy | 1.17.0 | Bio.Cluster, Bio.PDB, statistical tests |
| matplotlib | 3.10.8 | Bio.Graphics, Bio.Phylo plotting |
| reportlab | 4.4.9 | Bio.Graphics PDF/PostScript output |
| networkx | 3.6.1 | Bio.Phylo graph algorithms |
| rdflib | 7.5.0 | Bio.Phylo CDAO parser |
| igraph | 1.0.0 | Bio.Phylo (alternative graph backend) |
| mmtf-python | 1.1.3 | Bio.PDB MMTF binary format |

### Not installed / unavailable
- **psycopg2 / mysqlclient** — no PostgreSQL or MySQL servers available; BioSQL tests will skip.
- **pygraphviz / pydot** — not installed; Phylo graph-layout tests that depend on Graphviz will skip.
- **torch / tensorflow / jax** — none present; no GPU-accelerated ML frameworks available.
- **pre-commit** — not pre-installed; `pip install pre-commit` before using hooks.

### Toolchains

**C/C++ compiler:** Apple Clang 17.0.0 (`/usr/bin/clang`). Note that `/usr/bin/gcc` is actually clang (Xcode shim). All 12 Biopython C extensions compile and import successfully under this toolchain.

**Rust:** rustc 1.92.0 and cargo 1.92.0 are installed via Homebrew. The active toolchain is `stable-x86_64-apple-darwin`. Only the `x86_64-apple-darwin` target is available (no `aarch64` / ARM cross-compilation set up).

**jellyfish (k-mer counting):** jellyfish 2.3.1 is installed at `/usr/local/bin/jellyfish` and fully functional. Verified end-to-end: `count`, `histo`, and `dump` all work. Canonical mode (`-C`) is supported. Typical workflow:
```bash
jellyfish count -m 21 -s 100M -t 4 -C input.fa -o out.jf
jellyfish histo out.jf
jellyfish dump -c out.jf          # tab-separated k-mer counts
```
Limitation: single-threaded performance is the practical ceiling here (4 logical cores on one physical Quad-Core i7; no hyper-threading benefit beyond 4 threads).

### Key compatibility notes
- The environment is **x86_64 only**. Wheels and compiled extensions will be the `x86_64` variants even though some cached wheels reference `macosx_14_0` (the pip resolver picks compatible wheels).
- `numpy 2.4.2` is installed. CI explicitly excludes `numpy==2.1.0` (`ci-dependencies.txt`). If a test failure looks ABI-related, check for numpy version mismatches.
- Biopython emits a `BiopythonWarning` when imported from inside the source tree (`setup.py` is detected in the working directory). This is harmless during development but means tests should ideally be run from `Tests/` rather than the repo root.

## Documentation Requirements

All new or modified public-facing functionality must be accompanied by documentation updates:

1. **Docstrings are mandatory** on all public classes, methods, and functions. Use **reStructuredText (RST)** markup — not Google or NumPy docstring style. The Sphinx `autodoc` system pulls these into the API reference.
2. **Tutorial chapters** live in `Doc/Tutorial/`. If the feature is user-facing and not obvious, add or update the relevant chapter. Each chapter directory contains a `.tex` source and accompanying `.py` example scripts.
3. **RST files** (README.rst, NEWS.rst, DEPRECATED.rst) are validated by `rstcheck` and `doc8` in pre-commit. Run them before pushing.
4. **Code examples in docs** are formatted by `blacken-docs` (a pre-commit hook). They must be valid, runnable Python.
5. **NEWS.rst** — new features and bug fixes should be noted here for the next release. Follow the existing entry format (version heading, bullet points with issue/PR links).
6. **Doc build verification**: after touching anything under `Doc/`, build locally with `cd Doc && make html` and check for warnings. CI gates the docs job on all tests passing first.

## Contribution Requirements

These are the baseline rules that all code changes (bug fixes, features, refactors) must satisfy before merging:

1. **Tests are required.** Every new feature or bug fix must include a corresponding test in `Tests/`. The project does not accept untested contributions. Tests must pass with `--offline` if they do not genuinely require network access.
2. **Dual-license agreement.** Each commit or PR must explicitly state agreement to dual-licensing under the Biopython License Agreement and the 3-clause BSD License. This is a legal requirement stated in `CONTRIBUTING.rst`.
3. **Style must pass CI.** The GitHub Actions `style` job runs pre-commit on changed files. Locally: run `black`, `flake8`, and `ruff` on any file you touch. Do not leave trailing whitespace or missing end-of-file newlines (checked by pre-commit hooks).
4. **Type annotations on new code in Bio/ and BioSQL/.** `mypy` is run as a pre-commit hook with strict settings (see `.mypy.ini`). New public APIs should be fully annotated. Partially-typed modules (mostly in `Bio.PDB`) have per-module `ignore_errors` in `.mypy.ini` — do not expand that list.
5. **Do not break the offline test suite.** Many CI matrix entries run with minimal dependencies. If your code requires an optional package, guard the import and skip the test with `MissingPythonDependencyError` or `unittest.skipIf`, following the pattern already used throughout the codebase.
6. **Module naming.** Existing module names in Bio/ are intentionally not all lower_case. Do not rename them. New modules should follow PEP 8 (lower_case_with_underscores).
7. **C extension changes require extra care.** CI runs AddressSanitizer (`-fsanitize=address,undefined`) on Linux. Any change to code under `Bio/` that touches the `.c` or `.pyx` sources must not introduce memory leaks or undefined behaviour — the sanitiser build will catch them.
8. **CI must be green.** All GitHub Actions checks (style, build, test on Linux/macOS/Windows/PyPy, docs) must pass. CodeCov coverage checks are advisory but declining coverage on changed lines will be flagged in review.
