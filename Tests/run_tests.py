#!/usr/bin/env python
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Backward-compatibility shim: forwards all invocations to pytest.

Supported legacy invocations (all still work unchanged):
    python run_tests.py --offline
    python run_tests.py test_Seq
    python run_tests.py test_Seq.py
    python run_tests.py doctest
    coverage run --source Bio,BioSQL run_tests.py --offline
"""
import sys

import pytest


def main():
    args = sys.argv[1:]
    pytest_args = []

    # --- flags ---
    if "--offline" in args:
        args.remove("--offline")
        pytest_args.append("--offline")
    if "-v" in args or "--verbose" in args:
        args = [a for a in args if a not in ("-v", "--verbose")]
        pytest_args.append("-v")
    if "--help" in args:
        pytest_args.append("--help")
        args = [a for a in args if a != "--help"]

    # --- positional arguments (test names / "doctest") ---
    for arg in args:
        if arg == "doctest":
            # Module doctests are collected automatically via test_doctests.py;
            # "doctest" alone (no other args) is a no-op here.  If other test
            # names were also given, they are handled below.
            continue
        name = arg[:-3] if arg.endswith(".py") else arg
        if name.startswith("test_"):
            pytest_args.append(name + ".py")
        else:
            pytest_args.append(arg)

    sys.exit(pytest.main(pytest_args))


if __name__ == "__main__":
    main()
