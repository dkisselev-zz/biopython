"""Run doctests from all Bio.* and BioSQL.* source modules."""
import doctest
import importlib
import unittest

import pytest

from Bio import MissingExternalDependencyError
from Bio import MissingPythonDependencyError


def test_module_doctest(bio_module_name):
    """Doctest one Bio/BioSQL module.  Parametrised via conftest.pytest_generate_tests."""
    try:
        module = importlib.import_module(bio_module_name)
    except (MissingPythonDependencyError, MissingExternalDependencyError) as exc:
        pytest.skip(str(exc))
    except ImportError as exc:
        pytest.skip(f"ImportError: {exc}")

    suite = doctest.DocTestSuite(module, optionflags=doctest.ELLIPSIS)
    if suite.countTestCases() == 0:
        pytest.skip(f"No doctests in {bio_module_name}")

    result = unittest.TestResult()
    suite.run(result)

    if result.errors:
        pytest.fail(
            f"Doctest errors in {bio_module_name}:\n"
            + "\n".join(tb for _, tb in result.errors)
        )
    if result.failures:
        pytest.fail(
            f"Doctest failures in {bio_module_name}:\n"
            + "\n".join(tb for _, tb in result.failures)
        )
