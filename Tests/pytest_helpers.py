"""Pytest helpers for Biopython test suite.

Provides decorators and utilities for handling dependencies
and skipping tests with missing requirements.

Performance optimizations:
- _check_import and _check_external_tool are cached with @lru_cache
- Subprocess calls for tool availability only run once per tool/flag combination
"""

import pytest
import subprocess
from functools import wraps, lru_cache


@lru_cache(maxsize=None)
def _check_import(module_name):
    """Check if a Python module can be imported.

    Cached to avoid repeated import attempts for the same module.
    """
    try:
        __import__(module_name)
        return True
    except ImportError:
        return False


@lru_cache(maxsize=None)
def _check_external_tool(tool_name, version_flag="--version", timeout=5):
    """Check if external command-line tool is available.

    Cached to avoid repeated subprocess calls for the same tool.
    Each unique combination of (tool_name, version_flag, timeout) is
    checked only once and the result is cached for the entire session.
    """
    try:
        subprocess.check_output(
            [tool_name, version_flag],
            stderr=subprocess.DEVNULL,
            timeout=timeout
        )
        return True
    except (FileNotFoundError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return False


# ============================================================================
# Python Dependency Decorators
# ============================================================================

def requires_numpy(func):
    """Skip test if NumPy is not available."""
    return pytest.mark.skipif(
        not _check_import("numpy"),
        reason="Install NumPy to run this test"
    )(func)


def requires_scipy(func):
    """Skip test if SciPy is not available."""
    return pytest.mark.skipif(
        not _check_import("scipy"),
        reason="Install SciPy to run this test"
    )(func)


def requires_matplotlib(func):
    """Skip test if matplotlib is not available."""
    return pytest.mark.skipif(
        not _check_import("matplotlib"),
        reason="Install matplotlib to run this test"
    )(func)


def requires_networkx(func):
    """Skip test if networkx is not available."""
    return pytest.mark.skipif(
        not _check_import("networkx"),
        reason="Install networkx to run this test"
    )(func)


def requires_igraph(func):
    """Skip test if igraph is not available."""
    return pytest.mark.skipif(
        not _check_import("igraph"),
        reason="Install igraph to run this test"
    )(func)


def requires_rdflib(func):
    """Skip test if rdflib is not available."""
    return pytest.mark.skipif(
        not _check_import("rdflib"),
        reason="Install rdflib to run this test"
    )(func)


def requires_reportlab(func):
    """Skip test if reportlab is not available."""
    return pytest.mark.skipif(
        not _check_import("reportlab"),
        reason="Install reportlab to run this test"
    )(func)


def requires_mmtf(func):
    """Skip test if mmtf-python is not available."""
    return pytest.mark.skipif(
        not _check_import("mmtf"),
        reason="Install mmtf-python to run this test"
    )(func)


# ============================================================================
# External Tool Decorators
# ============================================================================

def requires_external_tool(tool_name, version_flag="--version"):
    """Decorator to skip test if external tool is not available.

    Args:
        tool_name: Name of the command-line tool
        version_flag: Flag to check tool availability (default: --version)

    Example:
        @requires_external_tool("dssp")
        def test_dssp():
            ...
    """
    tool_available = _check_external_tool(tool_name, version_flag)

    def decorator(func):
        return pytest.mark.skipif(
            not tool_available,
            reason=f"{tool_name} not available"
        )(func)

    return decorator


# ============================================================================
# Pre-computed Tool Availability
# ============================================================================

# Pre-compute availability of common tools at module import time
# This avoids repeated subprocess calls during test collection
DSSP_AVAILABLE = _check_external_tool("dssp")
NACCESS_AVAILABLE = _check_external_tool("naccess", "-q")
PSEA_AVAILABLE = _check_external_tool("psea", "-h")
EMBOSS_AVAILABLE = _check_external_tool("needle", "-version")
CLUSTALW_AVAILABLE = _check_external_tool("clustalw2", "-version")
MUSCLE_AVAILABLE = _check_external_tool("muscle", "-version")
MAFFT_AVAILABLE = _check_external_tool("mafft", "--version")
PRANK_AVAILABLE = _check_external_tool("prank")
PROBCONS_AVAILABLE = _check_external_tool("probcons")
TCoffee_AVAILABLE = _check_external_tool("t_coffee", "-version")
MSAPROBS_AVAILABLE = _check_external_tool("msaprobs")
DIALIGN2_AVAILABLE = _check_external_tool("dialign2-2")


# ============================================================================
# Export Commonly Used Skip Markers
# ============================================================================

# External tools
requires_dssp = pytest.mark.skipif(not DSSP_AVAILABLE, reason="DSSP not available")
requires_naccess = pytest.mark.skipif(not NACCESS_AVAILABLE, reason="NACCESS not available")
requires_psea = pytest.mark.skipif(not PSEA_AVAILABLE, reason="PSEA not available")
requires_emboss = pytest.mark.skipif(not EMBOSS_AVAILABLE, reason="EMBOSS not available")
requires_clustalw = pytest.mark.skipif(not CLUSTALW_AVAILABLE, reason="ClustalW not available")
requires_muscle = pytest.mark.skipif(not MUSCLE_AVAILABLE, reason="MUSCLE not available")
requires_mafft = pytest.mark.skipif(not MAFFT_AVAILABLE, reason="MAFFT not available")
requires_prank = pytest.mark.skipif(not PRANK_AVAILABLE, reason="PRANK not available")
requires_probcons = pytest.mark.skipif(not PROBCONS_AVAILABLE, reason="PROBCONS not available")
requires_tcoffee = pytest.mark.skipif(not TCoffee_AVAILABLE, reason="T-Coffee not available")
requires_msaprobs = pytest.mark.skipif(not MSAPROBS_AVAILABLE, reason="MSAProbs not available")
requires_dialign = pytest.mark.skipif(not DIALIGN2_AVAILABLE, reason="DIALIGN2 not available")


# ============================================================================
# Utility Functions
# ============================================================================

def check_tool_version(tool_name, min_version=None):
    """Check if a tool is available and optionally check minimum version.

    Args:
        tool_name: Name of the command-line tool
        min_version: Minimum required version (tuple), e.g., (2, 0, 0)

    Returns:
        tuple: (available: bool, version: tuple or None)
    """
    try:
        output = subprocess.check_output(
            [tool_name, "--version"],
            stderr=subprocess.STDOUT,
            text=True,
            timeout=5
        )

        # Try to parse version from output
        # This is a simple heuristic and may need adjustment per tool
        import re
        version_match = re.search(r'(\d+)\.(\d+)(?:\.(\d+))?', output)
        if version_match:
            version = tuple(int(x) if x else 0 for x in version_match.groups())
            if min_version and version < min_version:
                return (False, version)
            return (True, version)
        return (True, None)
    except (FileNotFoundError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return (False, None)
