"""Pytest configuration and fixtures for the Biopython test suite.

Replaces every implicit behaviour of the legacy run_tests.py runner:
  - --offline flag (poisons urlopen + requires_internet)
  - skip conversion for MissingXxxDependencyError at collection time
  - per-module doctest parametrisation
  - Tutorial RST doctest collection
  - CWD and LANG safety fixtures
"""

from __future__ import annotations

import doctest
import os
import sys
import types
import warnings
from collections.abc import Generator, Iterator
from pathlib import Path
from typing import Any, NoReturn, TextIO

import pytest

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
TESTS_DIR = os.path.abspath(os.path.dirname(__file__))

# Ensure Tests/ is on sys.path so helpers like requires_internet are importable
if TESTS_DIR not in sys.path:
    sys.path.insert(0, TESTS_DIR)

# Modules that must never be collected for doctests
ALWAYS_EXCLUDE_DOCTEST = {"Bio.Alphabet"}

# Modules excluded only under --offline
ONLINE_DOCTEST_MODULES = {
    "Bio.Entrez",
    "Bio.ExPASy",
    "Bio.ExPASy.cellosaurus",
    "Bio.TogoWS",
    "Bio.UniProt",
}

# Modules that require numpy (ported verbatim from run_tests.py)
NUMPY_DEPENDENT_MODULES = {
    "Bio.Affy.CelFile",
    "Bio.Align",
    "Bio.Align.substitution_matrices",
    "Bio.Cluster",
    "Bio.kNN",
    "Bio.LogisticRegression",
    "Bio.MarkovModel",
    "Bio.MaxEntropy",
    "Bio.NaiveBayes",
    "Bio.PDB.Chain",
    "Bio.PDB.Dice",
    "Bio.PDB.HSExposure",
    "Bio.PDB.MMCIF2Dict",
    "Bio.PDB.MMCIFParser",
    "Bio.PDB.mmtf.DefaultParser",
    "Bio.PDB.mmtf",
    "Bio.PDB.Model",
    "Bio.PDB.NACCESS",
    "Bio.PDB.NeighborSearch",
    "Bio.PDB.parse_pdb_header",
    "Bio.PDB.PDBExceptions",
    "Bio.PDB.PDBList",
    "Bio.PDB.PDBParser",
    "Bio.PDB.Polypeptide",
    "Bio.PDB.PSEA",
    "Bio.PDB.Residue",
    "Bio.PDB.Selection",
    "Bio.PDB.StructureAlignment",
    "Bio.PDB.StructureBuilder",
    "Bio.PDB.Structure",
    "Bio.PDB.Superimposer",
    "Bio.PDB.Vector",
    "Bio.phenotype",
    "Bio.phenotype.parse",
    "Bio.phenotype.phen_micro",
    "Bio.phenotype.pm_fitting",
    "Bio.SeqIO.PdbIO",
    "Bio.SVDSuperimposer",
}


# ---------------------------------------------------------------------------
# Command-line option and markers
# ---------------------------------------------------------------------------
def pytest_addoption(parser: pytest.Parser) -> None:
    """Register the --offline command-line flag."""
    parser.addoption(
        "--offline",
        action="store_true",
        default=False,
        help="Skip tests which require internet access.",
    )


def pytest_configure(config: pytest.Config) -> None:
    """Register markers and apply --offline side-effects."""
    config.addinivalue_line(
        "markers", "online: test requires internet access; skipped under --offline"
    )
    if config.getoption("--offline", default=False):
        # Poison urllib so any accidental network access is caught immediately
        import urllib.request

        def _dummy_urlopen(url: str) -> NoReturn:
            raise RuntimeError(
                "Internal test suite error, attempting to use internet "
                "despite --offline setting"
            )

        urllib.request.urlopen = _dummy_urlopen

        # Poison the requires_internet helper used by many legacy test files
        import requires_internet

        requires_internet.check.available = False


def pytest_collection_modifyitems(
    config: pytest.Config, items: list[pytest.Item]
) -> None:
    """Skip @pytest.mark.online tests when --offline is given."""
    if not config.getoption("--offline"):
        return
    skip_online = pytest.mark.skip(reason="skipped under --offline")
    for item in items:
        if "online" in item.keywords:
            item.add_marker(skip_online)


class _SafeModule(pytest.Module):
    """Module collector that converts Bio dependency errors into clean skips.

    When a test module raises MissingPythonDependencyError,
    MissingExternalDependencyError, or an ImportError mentioning Bio.Alphabet
    at import time, pytest.skip() is called so the module is reported as
    skipped rather than as a collection error.
    """

    def _getobj(self) -> types.ModuleType:
        try:
            return super()._getobj()
        except Exception as exc:
            # Two cases:
            # 1) MissingExternalDependencyError (NOT an ImportError subclass) –
            #    propagates as-is; type(exc).__qualname__ identifies it.
            # 2) MissingPythonDependencyError (IS an ImportError subclass) –
            #    pytest wraps it in CollectError; the original lives in __cause__
            #    and the type name appears in str(exc).
            cause = exc.__cause__ if exc.__cause__ else exc
            cause_type = type(cause).__qualname__
            if "MissingPythonDependencyError" in cause_type:
                pytest.skip(str(cause))
            elif "MissingExternalDependencyError" in cause_type:
                pytest.skip(str(cause))
            elif "Bio.Alphabet" in str(cause) and "ImportError" in cause_type:
                pytest.skip("Bio.Alphabet has been removed")
            # Fallback: scan the full string (covers double-wrapped cases)
            msg = str(exc)
            if "MissingPythonDependencyError" in msg:
                pytest.skip(_extract_skip_reason(msg))
            elif "MissingExternalDependencyError" in msg:
                pytest.skip(_extract_skip_reason(msg))
            raise


def _extract_skip_reason(msg: str) -> str:
    """Pull the skip reason from a CollectError traceback string."""
    for marker in (
        "Bio.MissingPythonDependencyError: ",
        "Bio.MissingExternalDependencyError: ",
    ):
        idx = msg.find(marker)
        if idx != -1:
            return msg[idx + len(marker) :].split("\n")[0].strip()
    return "Missing dependency"


@pytest.hookimpl(tryfirst=True)
def pytest_pycollect_makemodule(
    module_path: Path, parent: pytest.Collector
) -> _SafeModule:
    """Use SafeModule for all test files to handle dependency-skip patterns."""
    return _SafeModule.from_parent(parent, path=module_path)


# ---------------------------------------------------------------------------
# Module-level doctest parametrisation  (feeds test_doctests.py)
# ---------------------------------------------------------------------------
def pytest_generate_tests(metafunc: pytest.Metafunc) -> None:
    """Parametrise test_module_doctest with the Bio.* module list."""
    if "bio_module_name" not in metafunc.fixturenames:
        return
    offline = metafunc.config.getoption("--offline", default=False)
    modules = _find_bio_modules()
    excluded = _excluded_doctest_modules(offline)
    metafunc.parametrize("bio_module_name", sorted(modules - excluded))


def _find_bio_modules() -> set[str]:
    """Return set of importable Bio.* / BioSQL.* module names.

    Walks the directory tree under Bio/ and BioSQL/ using only pkgutil —
    no setuptools required.
    """
    from pkgutil import iter_modules

    root = os.path.join(TESTS_DIR, "..")
    modules: set[str] = set()

    def _walk(pkg_name: str) -> None:
        modules.add(pkg_name)
        pkg_path = os.path.join(root, pkg_name.replace(".", os.sep))
        for info in iter_modules([pkg_path]):
            full_name = pkg_name + "." + info.name
            modules.add(full_name)
            if info.ispkg:
                _walk(full_name)

    for top in ("Bio", "BioSQL"):
        top_path = os.path.join(root, top)
        if os.path.isdir(top_path) and os.path.exists(
            os.path.join(top_path, "__init__.py")
        ):
            _walk(top)

    return modules


def _excluded_doctest_modules(offline: bool) -> set[str]:
    """Build the exclusion set: Alphabet + numpy-dep + sqlite3-dep + online."""
    excluded = set(ALWAYS_EXCLUDE_DOCTEST)

    try:
        import numpy  # noqa: F401
    except ImportError:
        excluded.update(NUMPY_DEPENDENT_MODULES)

    try:
        import sqlite3  # noqa: F401
    except ImportError:
        excluded.add("Bio.SeqIO")
        excluded.add("Bio.SearchIO")

    if offline:
        excluded.update(ONLINE_DOCTEST_MODULES)

    return excluded


# ---------------------------------------------------------------------------
# Autouse fixtures – safety nets that every test_*.py test inherits
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _cwd_safety() -> Generator[None, None, None]:
    """Guarantee CWD == Tests/ around each test; fail if a test changed it."""
    os.chdir(TESTS_DIR)
    yield
    current = os.path.abspath(".")
    os.chdir(TESTS_DIR)  # Restore first so the next test isn't affected
    if current != TESTS_DIR:
        pytest.fail(
            f"Test changed the current working directory to {current}. "
            "Restore it in a fixture or tearDown.",
            pytrace=False,
        )


@pytest.fixture(autouse=True)
def _lang_reset() -> Generator[None, None, None]:
    """Restore LANG env var after each test (some tests tweak it)."""
    original = os.environ.get("LANG", "C")
    yield
    os.environ["LANG"] = original


# ---------------------------------------------------------------------------
# Tutorial RST doctest collector
# ---------------------------------------------------------------------------
def pytest_collect_file(
    parent: pytest.Collector, file_path: Path
) -> TutorialDoctestFile | None:
    """Collect Tutorial chapter_*.rst files as doctest suites."""
    if (
        file_path.suffix == ".rst"
        and file_path.name.startswith("chapter_")
        and "Tutorial" in file_path.parts
    ):
        return TutorialDoctestFile.from_parent(parent, path=file_path)


class TutorialDoctestFile(pytest.File):
    """One pytest File node per Tutorial chapter_*.rst."""

    def collect(
        self,
    ) -> Iterator[TutorialDoctestItem | _TutorialParseErrorItem]:
        tutorial_base = str(self.path.parent.parent)  # Doc/ (parent of Tutorial/)
        try:
            doctests = _extract_tutorial_doctests(str(self.path))
        except ValueError as exc:
            yield _TutorialParseErrorItem.from_parent(
                self,
                name=f"PARSE_ERROR_{self.path.name}",
                error_message=str(exc),
            )
            return
        for item in doctests:
            yield TutorialDoctestItem.from_parent(
                self,
                name=item["name"],
                callobj=item,
                tutorial_base=tutorial_base,
            )


class _TutorialParseErrorItem(pytest.Item):
    """Placeholder that surfaces an RST-parsing ValueError as a single failure.

    Yielded by TutorialDoctestFile.collect when _extract_tutorial_doctests
    raises, so the rest of the suite continues to collect normally.
    """

    def __init__(
        self,
        name: str,
        parent: pytest.Collector,
        error_message: str,
        **kwargs: Any,
    ) -> None:
        super().__init__(name, parent, **kwargs)
        self._error_message = error_message

    def runtest(self) -> None:
        pytest.fail(self._error_message, pytrace=False)

    def repr_failure(self, excinfo: pytest.ExceptionInfo[BaseException]) -> str:
        return str(excinfo.value)

    def reportinfo(self) -> tuple[Path, int, str]:
        return self.path, 0, self.name


class TutorialDoctestItem(pytest.Item):
    """Single doctest block extracted from a Tutorial RST file."""

    def __init__(
        self,
        name: str,
        parent: pytest.Collector,
        callobj: dict[str, str | list[str]],
        tutorial_base: str,
        **kwargs: Any,
    ) -> None:
        super().__init__(name, parent, **kwargs)
        self._callobj = callobj
        self._tutorial_base = tutorial_base
        self._original_cwd: str | None = None

    def setup(self) -> None:
        """Save CWD and skip if any lib:/internet dependency is missing."""
        self._original_cwd = os.getcwd()
        deps = self._callobj["deps"]
        missing = []
        for dep in deps:
            if dep == "internet":
                if self.config.getoption("--offline", default=False):
                    missing.append("internet")
            else:
                assert dep.startswith("lib:"), dep
                lib = dep[4:]
                try:
                    __import__(lib)
                except ImportError:
                    missing.append(lib)
        if missing:
            pytest.skip(f"Missing: {', '.join(missing)}")

    def runtest(self) -> None:
        """Execute the doctest block, optionally chdir-ing first."""
        example = self._callobj["example"]
        folder = self._callobj["folder"]
        if folder:
            workdir = os.path.join(self._tutorial_base, folder)
            code = f">>> import os\n>>> os.chdir({workdir!r})\n{example}"
        else:
            code = example

        with warnings.catch_warnings():
            from Bio import BiopythonDeprecationWarning
            from Bio import BiopythonExperimentalWarning

            warnings.simplefilter("ignore", BiopythonDeprecationWarning)
            warnings.simplefilter("ignore", BiopythonExperimentalWarning)

            parser = doctest.DocTestParser()
            globs = {}
            test = parser.get_doctest(code, globs, self.name, str(self.path), 0)
            runner = doctest.DocTestRunner(optionflags=doctest.ELLIPSIS)
            results = runner.run(test)

        if results.failed:
            from io import StringIO

            out = StringIO()
            runner2 = doctest.DocTestRunner(
                optionflags=doctest.ELLIPSIS, verbose=True
            )
            parser2 = doctest.DocTestParser()
            test2 = parser2.get_doctest(code, {}, self.name, str(self.path), 0)
            runner2.run(test2, out=out.write)
            pytest.fail(out.getvalue())

    def teardown(self) -> None:
        """Restore CWD and clean up tutorial output files."""
        if self._original_cwd:
            os.chdir(self._original_cwd)
        # Files created by chapter_phylo doctests (paths relative to Doc/)
        for fname in ("examples/tree1.nwk", "examples/other_trees.xml"):
            path = os.path.join(self._tutorial_base, fname)
            if os.path.exists(path):
                os.remove(path)
        # Files created by chapter_cluster doctests (paths relative to Tests/)
        tests_dir = os.path.join(self._tutorial_base, "..", "Tests")
        for fname in (
            "Cluster/cyano_result.atr",
            "Cluster/cyano_result.cdt",
            "Cluster/cyano_result.gtr",
            "Cluster/cyano_result_K_A2.kag",
            "Cluster/cyano_result_K_G5.kgg",
            "Cluster/cyano_result_K_G5_A2.cdt",
        ):
            path = os.path.join(tests_dir, fname)
            if os.path.exists(path):
                os.remove(path)

    def repr_failure(self, excinfo: pytest.ExceptionInfo[BaseException]) -> str:
        """Short failure representation."""
        return str(excinfo.value)

    def reportinfo(self) -> tuple[Path, int, str]:
        return self.path, 0, self.name


# ---------------------------------------------------------------------------
# RST parsing helpers (ported from test_Tutorial.py)
# ---------------------------------------------------------------------------
def _extract(handle: TextIO) -> list[str]:
    """Read one code block after a doctest/cont-doctest directive.

    Faithful port of test_Tutorial._extract().
    """
    line = handle.readline()
    if line != "\n":
        raise ValueError(
            "Any '.. doctest' or '.. cont-doctest' line should "
            "be followed by an empty line"
        )

    line = handle.readline()
    if line.lstrip() != ".. code:: pycon\n":
        raise ValueError(
            "Any '.. doctest' or '.. cont-doctest' line should "
            r"be followed by '\n.. code:: pycon\n\n"
        )

    line = handle.readline()
    if line != "\n":
        raise ValueError(
            "Any '.. doctest' or '.. cont-doctest' line should "
            r"be followed by '\n.. code:: pycon\n\n"
        )

    lines = []
    while True:
        line = handle.readline()
        if not line:
            if lines:
                break
            else:
                raise ValueError("Didn't find lines!")
        elif line == "\n":
            break
        else:
            lines.append(line)
    return lines


def _extract_tutorial_doctests(rst_filename: str) -> list[dict[str, str | list[str]]]:
    """Scan one Tutorial chapter RST and return doctest dicts.

    Faithful port of test_Tutorial.extract_doctests().
    Each dict has keys: name, example, folder, deps.
    """
    base_name = os.path.splitext(os.path.basename(rst_filename))[0]
    name = None
    deps = []
    folder = ""
    results = []

    with open(rst_filename, encoding="utf8") as handle:
        line_number = 0
        lines = []
        while True:
            line = handle.readline()
            line_number += 1
            if not line:
                # End of file
                break
            elif line.lstrip().startswith(".. cont-doctest"):
                x = _extract(handle)
                lines.extend(x)
                line_number += len(x) + 2
            elif line.lstrip().startswith(".. doctest"):
                if lines:
                    if not lines[0].lstrip().startswith(">>> "):
                        raise ValueError(
                            f"Should start with '>>> ' (indented), not {lines[0]!r}"
                        )
                    results.append(
                        {
                            "name": name,
                            "example": "".join(lines),
                            "folder": folder,
                            "deps": deps,
                        }
                    )
                    lines = []
                deps = [x.strip() for x in line.split()[2:]]
                if deps:
                    folder = deps[0]
                    deps = deps[1:]
                else:
                    folder = ""
                name = "test_%s_line_%05i" % (base_name, line_number)
                x = _extract(handle)
                lines.extend(x)
                line_number += len(x) + 2
    if lines:
        if not name:
            raise ValueError(f"Unanchored doctest in {rst_filename}: {lines}")
        if not lines[0].lstrip().startswith(">>> "):
            raise ValueError(f"Should start '>>> ' not {lines[0]!r}")
        results.append(
            {
                "name": name,
                "example": "".join(lines),
                "folder": folder,
                "deps": deps,
            }
        )

    return results
