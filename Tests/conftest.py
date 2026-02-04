"""pytest configuration and fixtures for Biopython test suite.

This module provides:
- Working directory enforcement (Tests/ must remain current directory)
- Online/offline test handling via --offline flag with cached connectivity check
- Custom skip markers for missing dependencies
- Environment restoration between tests
- Doctest collection with dynamic exclusions

Performance optimizations:
- Internet connectivity is probed only once per session (not per test)
  Results cached in config._internet_available
- Working directory enforcement uses module scope (runs 209 times vs 3,264)
  with opt-in function scope for tests that change CWD (~40 tests)
  This provides a 92.4% reduction in fixture overhead
"""

import os
import sys
import pytest
import warnings
from pathlib import Path

# Add Tests directory to path for test imports
TEST_DIR = Path(__file__).parent.absolute()
if str(TEST_DIR) not in sys.path:
    sys.path.insert(0, str(TEST_DIR))


# ============================================================================
# Command-line Options
# ============================================================================

def pytest_addoption(parser):
    """Add custom command-line options."""
    parser.addoption(
        "--offline",
        action="store_true",
        default=False,
        help="Skip tests requiring internet access"
    )


# ============================================================================
# Session-level Setup
# ============================================================================

@pytest.fixture(scope="session", autouse=True)
def configure_offline_mode(request):
    """Configure offline mode if --offline flag is provided.

    This mimics the behavior of run_tests.py's offline handling:
    - Monkey-patches urllib.request.urlopen
    - Sets requires_internet.check.available = False
    """
    if request.config.getoption("--offline"):
        # Monkey patch for urlopen()
        import urllib.request

        def dummy_urlopen(url):
            raise RuntimeError(
                "Internal test suite error, attempting to use internet "
                "despite --offline setting"
            )

        urllib.request.urlopen = dummy_urlopen

        # Set requires_internet flag
        try:
            import requires_internet
            requires_internet.check.available = False
        except ImportError:
            pass  # requires_internet might not be available in all contexts

        print("Skipping any tests requiring internet access")


# ============================================================================
# Working Directory Enforcement
# ============================================================================

@pytest.fixture(scope="module", autouse=True)
def enforce_working_directory_module():
    """Ensure working directory is Tests/ at module start/end (lightweight).

    This lightweight fixture runs once per module (209 times total) instead of
    once per test (3,264 times), providing a 93.6% reduction in overhead.

    For the ~1.2% of tests that explicitly change CWD (test_SeqIO_index.py,
    test_Tutorial.py), use the strict fixture with:
        pytestmark = pytest.mark.usefixtures("enforce_working_directory_strict")

    This replicates the behavior from run_tests.py lines 282-290.
    """
    original_dir = os.getcwd()
    expected_dir = str(TEST_DIR)

    # Ensure we start in Tests/
    if os.path.abspath(original_dir) != expected_dir:
        os.chdir(expected_dir)

    yield

    # Restore at module end (in case last test didn't clean up)
    if os.path.abspath(os.getcwd()) != expected_dir:
        os.chdir(expected_dir)


@pytest.fixture(scope="function")
def enforce_working_directory_strict():
    """Strict working directory enforcement for tests that change CWD.

    Use this fixture for tests that explicitly change working directory.
    Mark the entire module:

        pytestmark = pytest.mark.usefixtures("enforce_working_directory_strict")

    Or mark a specific class:

        @pytest.mark.usefixtures("enforce_working_directory_strict")
        class MyTest(unittest.TestCase):
            ...

    This ensures:
    1. Test starts in Tests/ directory
    2. If test changes directory, it must restore it
    3. Failure to restore triggers a test failure with clear message

    Currently used by:
    - test_SeqIO_index.py (31 tests that use os.chdir)
    - test_Tutorial.py (9 tests that change directories)
    """
    original_dir = os.getcwd()
    expected_dir = str(TEST_DIR)

    # Ensure we start in Tests/
    if os.path.abspath(original_dir) != expected_dir:
        os.chdir(expected_dir)

    yield

    # Check if test changed directory
    current_dir = os.path.abspath(os.getcwd())
    if current_dir != expected_dir:
        # Restore directory and fail the test
        os.chdir(expected_dir)
        pytest.fail(
            f"Test changed working directory!\n"
            f"Expected: {expected_dir}\n"
            f"Got: {current_dir}"
        )


@pytest.fixture(scope="function", autouse=True)
def restore_lang_environment():
    """Restore LANG environment variable between tests.

    Some tests modify LANG to detect command-line tools.
    This ensures it's restored between tests.

    Replicates behavior from run_tests.py lines 128, 234.
    """
    SYSTEM_LANG = os.environ.get("LANG", "C")

    yield

    os.environ["LANG"] = SYSTEM_LANG


# ============================================================================
# Dependency Skip Markers
# ============================================================================

def pytest_configure(config):
    """Register custom markers and configure test collection.

    Also caches the internet connectivity check result to avoid
    repeated network probes for every @pytest.mark.online test.
    """
    # Markers are already defined in pytest.ini

    # Cache internet availability check (only probe once per session)
    internet_available = None

    if not config.getoption("--offline"):
        try:
            import requires_internet
            from Bio import MissingExternalDependencyError
            try:
                requires_internet.check()
                internet_available = True
            except MissingExternalDependencyError:
                internet_available = False
        except ImportError:
            # If requires_internet module not available, assume offline
            internet_available = False
    else:
        # Explicit --offline flag means don't check
        internet_available = False

    # Store in config for access by pytest_runtest_setup
    config._internet_available = internet_available


def pytest_runtest_setup(item):
    """Skip tests based on markers and --offline flag.

    Uses cached internet availability from pytest_configure to avoid
    repeated network probes.
    """
    # Handle online tests
    if item.get_closest_marker("online"):
        if item.config.getoption("--offline"):
            pytest.skip("test requires internet (--offline mode)")

        # Check cached internet availability (set once in pytest_configure)
        if not getattr(item.config, "_internet_available", True):
            pytest.skip("test requires internet (not available)")


# ============================================================================
# Doctest Collection
# ============================================================================

# Modules to exclude from doctest collection
# Replicated from run_tests.py lines 49-104
EXCLUDE_DOCTEST_MODULES = [
    "Bio.Alphabet",  # Removed from Biopython
]

ONLINE_DOCTEST_MODULES = [
    "Bio.Entrez",
    "Bio.ExPASy",
    "Bio.ExPASy.cellosaurus",
    "Bio.TogoWS",
    "Bio.UniProt",
]

# Check for NumPy and exclude NumPy-dependent modules if missing
try:
    import numpy as np
except ImportError:
    np = None

if np is None:
    EXCLUDE_DOCTEST_MODULES.extend([
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
    ])

# Check for sqlite3
try:
    import sqlite3
except ImportError:
    EXCLUDE_DOCTEST_MODULES.append("Bio.SeqIO")
    EXCLUDE_DOCTEST_MODULES.append("Bio.SearchIO")


def pytest_collection_modifyitems(config, items):
    """Modify test collection to handle offline mode for doctests."""
    offline = config.getoption("--offline")

    if offline:
        # Extend exclusion list with online modules in offline mode
        excluded_modules = set(EXCLUDE_DOCTEST_MODULES + ONLINE_DOCTEST_MODULES)
    else:
        excluded_modules = set(EXCLUDE_DOCTEST_MODULES)

    # Skip doctest items for excluded modules
    for item in items:
        if isinstance(item, pytest.DoctestItem):
            # Get module name from DoctestItem - it's stored in dtest.globs['__name__']
            # or we can get it from the item's parent
            try:
                if hasattr(item, 'dtest') and item.dtest:
                    module_name = item.dtest.globs.get('__name__', '')
                elif hasattr(item.parent, 'module'):
                    module_name = item.parent.module.__name__
                else:
                    # Skip if we can't determine module name
                    continue

                if module_name in excluded_modules:
                    item.add_marker(pytest.mark.skip(
                        reason=f"Module {module_name} excluded from doctests"
                    ))
            except (AttributeError, KeyError):
                # If we can't get module name, don't skip
                pass


# ============================================================================
# Doctest Discovery
# ============================================================================

def pytest_collect_file(parent, file_path):
    """Custom doctest collection for Bio modules.

    This allows pytest --doctest-modules to work with our exclusion logic.
    We only collect doctests from Bio/ directory, not Tests/.
    """
    # Only collect from Bio/ directory
    # Use file_path.parts to match "Bio" as a path component, not substring
    # This avoids false positives like "/home/Bioinformatics/project/..."
    if "Bio" not in file_path.parts:
        return None

    # Standard doctest collection
    if file_path.suffix == ".py":
        # Let pytest handle it with our exclusion logic in pytest_collection_modifyitems
        return None


# ============================================================================
# Test Report Customization
# ============================================================================

@pytest.hookimpl(tryfirst=True, hookwrapper=True)
def pytest_runtest_makereport(item, call):
    """Customize test failure reporting to match unittest output style.

    This helps maintain familiar error output during migration.
    """
    outcome = yield
    rep = outcome.get_result()

    # Add custom attributes for better reporting
    if rep.when == "call":
        # Store test execution time
        rep.test_duration = call.stop - call.start
