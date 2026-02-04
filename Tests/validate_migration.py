#!/usr/bin/env python
"""Validation script to compare unittest and pytest test runs.

This script runs the test suite with both unittest (via run_tests.py) and
pytest, then compares the results to ensure the migration is complete and
accurate.

Usage:
    python validate_migration.py
    python validate_migration.py --verbose
"""

import argparse
import subprocess
import sys
import time
from pathlib import Path


def run_command(cmd, description):
    """Run a command and return timing and result information."""
    print(f"\n{'='*70}")
    print(f"Running: {description}")
    print(f"Command: {' '.join(cmd)}")
    print(f"{'='*70}")

    start_time = time.time()
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600  # 10 minute timeout
        )
        elapsed = time.time() - start_time
        return {
            'success': result.returncode == 0,
            'returncode': result.returncode,
            'stdout': result.stdout,
            'stderr': result.stderr,
            'elapsed': elapsed
        }
    except subprocess.TimeoutExpired:
        elapsed = time.time() - start_time
        return {
            'success': False,
            'returncode': -1,
            'stdout': '',
            'stderr': 'Timeout after 10 minutes',
            'elapsed': elapsed
        }


def parse_unittest_output(output):
    """Parse unittest output to extract test counts."""
    lines = output.split('\n')
    info = {
        'total': 0,
        'passed': 0,
        'skipped': 0,
        'failed': 0,
        'errors': 0
    }

    for line in lines:
        # Look for the final result line like:
        # "Ran 3264 tests in 120.456s"
        if line.startswith('Ran '):
            parts = line.split()
            if len(parts) >= 2:
                try:
                    info['total'] = int(parts[1])
                except ValueError:
                    pass

        # Look for status line like: "OK (skipped=25)"
        if line.startswith('OK'):
            info['passed'] = info['total']
            if 'skipped=' in line:
                try:
                    skip_str = line.split('skipped=')[1].split(')')[0]
                    info['skipped'] = int(skip_str)
                    info['passed'] -= info['skipped']
                except (ValueError, IndexError):
                    pass

        # Look for failure line
        if line.startswith('FAILED'):
            if 'failures=' in line:
                try:
                    fail_str = line.split('failures=')[1].split(',')[0].split(')')[0]
                    info['failed'] = int(fail_str)
                except (ValueError, IndexError):
                    pass
            if 'errors=' in line:
                try:
                    err_str = line.split('errors=')[1].split(',')[0].split(')')[0]
                    info['errors'] = int(err_str)
                except (ValueError, IndexError):
                    pass

    return info


def parse_pytest_output(output):
    """Parse pytest output to extract test counts."""
    lines = output.split('\n')
    info = {
        'total': 0,
        'passed': 0,
        'skipped': 0,
        'failed': 0,
        'errors': 0
    }

    for line in lines:
        # Look for summary line like:
        # "=== 3264 passed, 25 skipped in 120.45s ==="
        if 'passed' in line and 'in ' in line and '=' in line:
            parts = line.split()
            for i, part in enumerate(parts):
                if part == 'passed,' or part == 'passed':
                    try:
                        info['passed'] = int(parts[i-1])
                    except (ValueError, IndexError):
                        pass
                elif 'skipped' in part:
                    try:
                        info['skipped'] = int(parts[i-1])
                    except (ValueError, IndexError):
                        pass
                elif 'failed' in part:
                    try:
                        info['failed'] = int(parts[i-1])
                    except (ValueError, IndexError):
                        pass
                elif 'error' in part:
                    try:
                        info['errors'] = int(parts[i-1])
                    except (ValueError, IndexError):
                        pass

    info['total'] = info['passed'] + info['skipped'] + info['failed'] + info['errors']
    return info


def validate():
    """Run validation comparing unittest and pytest."""
    print("\n" + "="*70)
    print("Biopython Test Migration Validation")
    print("="*70)

    # Check we're in Tests directory
    if not Path('run_tests.py').exists():
        print("ERROR: Must run from Tests/ directory")
        print("Current directory:", Path.cwd())
        return False

    # Run unittest suite
    print("\n1. Running unittest suite (run_tests.py)...")
    unittest_result = run_command(
        ['python', 'run_tests.py', '--offline'],
        "Unittest test suite"
    )

    if not unittest_result['success']:
        print(f"WARNING: Unittest suite failed with return code {unittest_result['returncode']}")
        if unittest_result['stderr']:
            print("Stderr:", unittest_result['stderr'][:500])

    unittest_info = parse_unittest_output(unittest_result['stdout'])

    print(f"\nUnittest Results:")
    print(f"  Total tests:  {unittest_info['total']}")
    print(f"  Passed:       {unittest_info['passed']}")
    print(f"  Skipped:      {unittest_info['skipped']}")
    print(f"  Failed:       {unittest_info['failed']}")
    print(f"  Errors:       {unittest_info['errors']}")
    print(f"  Time:         {unittest_result['elapsed']:.2f}s")

    # Run pytest suite
    print("\n2. Running pytest suite...")
    pytest_result = run_command(
        ['pytest', '--offline', '-v'],
        "Pytest test suite"
    )

    if not pytest_result['success']:
        print(f"WARNING: Pytest suite failed with return code {pytest_result['returncode']}")
        if pytest_result['stderr']:
            print("Stderr:", pytest_result['stderr'][:500])

    pytest_info = parse_pytest_output(pytest_result['stdout'])

    print(f"\nPytest Results:")
    print(f"  Total tests:  {pytest_info['total']}")
    print(f"  Passed:       {pytest_info['passed']}")
    print(f"  Skipped:      {pytest_info['skipped']}")
    print(f"  Failed:       {pytest_info['failed']}")
    print(f"  Errors:       {pytest_info['errors']}")
    print(f"  Time:         {pytest_result['elapsed']:.2f}s")

    # Compare results
    print(f"\n{'='*70}")
    print("Comparison Analysis")
    print(f"{'='*70}")

    validation_passed = True

    # Compare test counts (±2% tolerance)
    if unittest_info['total'] > 0:
        total_diff_pct = abs(unittest_info['total'] - pytest_info['total']) / unittest_info['total'] * 100
        print(f"\nTest Count Comparison:")
        print(f"  Unittest: {unittest_info['total']}")
        print(f"  Pytest:   {pytest_info['total']}")
        print(f"  Diff:     {pytest_info['total'] - unittest_info['total']} ({total_diff_pct:.2f}%)")

        if total_diff_pct > 2.0:
            print(f"  ⚠️  WARNING: Test count differs by more than 2%")
            validation_passed = False
        else:
            print(f"  ✅ Test counts match within tolerance")

    # Compare skip counts (±10% tolerance)
    if unittest_info['skipped'] > 0:
        skip_diff_pct = abs(unittest_info['skipped'] - pytest_info['skipped']) / unittest_info['skipped'] * 100
        print(f"\nSkip Count Comparison:")
        print(f"  Unittest: {unittest_info['skipped']}")
        print(f"  Pytest:   {pytest_info['skipped']}")
        print(f"  Diff:     {pytest_info['skipped'] - unittest_info['skipped']} ({skip_diff_pct:.2f}%)")

        if skip_diff_pct > 10.0:
            print(f"  ⚠️  WARNING: Skip count differs by more than 10%")
            validation_passed = False
        else:
            print(f"  ✅ Skip counts match within tolerance")

    # Compare failure counts
    print(f"\nFailure Comparison:")
    print(f"  Unittest failures: {unittest_info['failed']}")
    print(f"  Pytest failures:   {pytest_info['failed']}")

    if unittest_info['failed'] == 0 and pytest_info['failed'] == 0:
        print(f"  ✅ Both suites passed all tests")
    elif unittest_info['failed'] == pytest_info['failed']:
        print(f"  ⚠️  Same number of failures in both suites")
    else:
        print(f"  ⚠️  Different failure counts")
        validation_passed = False

    # Performance comparison
    print(f"\nPerformance:")
    speedup = unittest_result['elapsed'] / pytest_result['elapsed'] if pytest_result['elapsed'] > 0 else 1.0
    print(f"  Unittest time: {unittest_result['elapsed']:.2f}s")
    print(f"  Pytest time:   {pytest_result['elapsed']:.2f}s")
    if speedup > 1.0:
        print(f"  Pytest is {speedup:.2f}x faster")
    elif speedup < 1.0:
        print(f"  Unittest is {1/speedup:.2f}x faster")
    else:
        print(f"  Similar performance")

    # Final verdict
    print(f"\n{'='*70}")
    if validation_passed:
        print("✅ VALIDATION PASSED")
        print("Migration is complete and both test runners produce similar results.")
    else:
        print("⚠️  VALIDATION WARNINGS")
        print("Some metrics differ beyond tolerance. Review output above.")
    print(f"{'='*70}\n")

    return validation_passed


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Show verbose output')
    args = parser.parse_args()

    try:
        success = validate()
        sys.exit(0 if success else 1)
    except KeyboardInterrupt:
        print("\n\nValidation interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n\nERROR during validation: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
