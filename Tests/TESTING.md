# Testing Guide for Biopython

This guide covers running tests in Biopython using both pytest (recommended) and the legacy unittest runner.

## Quick Start

```bash
cd Tests

# Run all tests with pytest (recommended)
pytest --offline

# Run specific test file
pytest test_SeqIO.py -v

# Run tests matching a pattern
pytest -k "test_SeqIO" --offline
```

## Prerequisites

Make sure you've installed Biopython in development mode before running tests:

```bash
# From repository root
pip install -e .

# Verify correct installation
python -c "import Bio; print(Bio.__file__)"
# Should show: .../model_a/Bio/__init__.py
```

**Important**: Always run `pip install -e .` after pulling changes to avoid import errors.

## Running Tests with Pytest (Recommended)

### Basic Usage

```bash
cd Tests

# All offline tests (recommended for local development)
pytest --offline

# All tests including online tests
pytest

# Verbose output
pytest --offline -v

# Very verbose (show test names and output)
pytest --offline -vv

# Quiet mode (minimal output)
pytest --offline -q
```

### Running Specific Tests

```bash
# Single test file
pytest test_SeqIO.py --offline

# Multiple files
pytest test_SeqIO.py test_Align.py --offline

# Test by name pattern
pytest -k "test_fasta" --offline

# Test by marker
pytest -m "not online" --offline

# Specific test function
pytest test_SeqIO.py::TestSeqIO::test_fasta --offline
```

### Parallel Execution

Run tests in parallel for faster execution:

```bash
# Auto-detect CPU count
pytest --offline -n auto

# Specific number of workers
pytest --offline -n 4

# Distribute by file (faster for many small test files)
pytest --offline -n auto --dist loadfile
```

### Test Markers

Biopython uses custom pytest markers to categorize tests:

```bash
# List all markers
pytest --markers

# Skip online tests
pytest -m "not online" --offline

# Run only numpy tests
pytest -m numpy

# Exclude slow tests
pytest -m "not slow"
```

Available markers:
- `online` - Tests requiring internet access
- `numpy` - Tests requiring NumPy
- `scipy` - Tests requiring SciPy
- `matplotlib` - Tests requiring matplotlib
- `external_tool` - Tests requiring command-line tools (DSSP, NACCESS, etc.)
- `slow` - Slow-running tests

### Coverage

Generate coverage reports:

```bash
# Run with coverage
pytest --offline --cov=Bio --cov=BioSQL --cov-report=html

# View report
open htmlcov/index.html
```

### Doctest Integration

Run doctests from Bio modules:

```bash
# Include doctests from Bio modules
pytest --doctest-modules --offline

# Specific module doctests
pytest ../Bio/Seq.py --doctest-modules

# Tutorial doctests
pytest test_Tutorial.py --offline
```

## Running Tests with Unittest (Legacy)

The original test runner still works for backward compatibility:

```bash
cd Tests

# All offline tests
python run_tests.py --offline

# All tests including online
python run_tests.py

# Specific test module
python run_tests.py test_SeqIO

# Single test file directly
python test_SeqIO.py

# Via setup.py
cd ..
python setup.py test --offline
```

## Test Organization

### Test Discovery

Pytest automatically discovers tests using these patterns:
- Files: `test_*.py` or `*_test.py`
- Classes: `Test*` (inheriting from `unittest.TestCase`)
- Functions: `test_*`

### Test Structure

```python
import unittest
import pytest

# Mark test requiring internet
@pytest.mark.online
class TestOnlineFeature(unittest.TestCase):
    def test_something(self):
        # Test code
        pass

# Mark test requiring NumPy
@pytest.mark.numpy
class TestNumPyFeature(unittest.TestCase):
    def test_array_operation(self):
        # Test code
        pass
```

### Module-level Skips

For tests requiring external dependencies:

```python
# At module level
import pytest

# Simple dependency check
pytest.importorskip("numpy", reason="NumPy required")

# External tool check (keep try/except for unittest compatibility)
try:
    import some_optional_module
except ImportError:
    from Bio import MissingExternalDependencyError
    raise MissingExternalDependencyError(
        "Install some_optional_module for these tests"
    ) from None
```

## Offline Mode

Use `--offline` to skip tests requiring internet:

```bash
# Pytest
pytest --offline

# Unittest
python run_tests.py --offline
```

Tests marked with `@pytest.mark.online` are automatically skipped in offline mode.

## Debugging Failed Tests

### Show Test Output

```bash
# Show print statements
pytest test_SeqIO.py -s

# Show detailed output on failure
pytest test_SeqIO.py -vv

# Stop on first failure
pytest test_SeqIO.py -x

# Drop into debugger on failure
pytest test_SeqIO.py --pdb
```

### Rerun Failed Tests

```bash
# Run only failed tests from last run
pytest --offline --lf

# Run failed tests first, then others
pytest --offline --ff
```

### Verbose Error Messages

```bash
# Show full diff for assertion failures
pytest test_SeqIO.py --offline -vv --tb=long

# Show only short traceback
pytest test_SeqIO.py --offline --tb=short

# Show only failing line
pytest test_SeqIO.py --offline --tb=line
```

## Writing New Tests

### Basic Test Structure

```python
import unittest
import pytest

class TestNewFeature(unittest.TestCase):
    """Tests for new feature."""

    def setUp(self):
        """Set up test fixtures."""
        # Run before each test method
        pass

    def tearDown(self):
        """Clean up after tests."""
        # Run after each test method
        pass

    def test_basic_functionality(self):
        """Test basic functionality."""
        result = some_function()
        self.assertEqual(result, expected_value)

    @pytest.mark.online
    def test_online_feature(self):
        """Test requiring internet."""
        # Will be skipped with --offline
        pass
```

### Using Pytest Fixtures

While we maintain unittest.TestCase for compatibility, you can use pytest fixtures:

```python
import pytest

@pytest.fixture
def sample_sequence():
    """Provide a sample sequence for tests."""
    return "ACGTACGT"

def test_with_fixture(sample_sequence):
    """Test using pytest fixture."""
    assert len(sample_sequence) == 8
```

### Parametrized Tests

```python
import pytest

@pytest.mark.parametrize("input,expected", [
    ("ACGT", 4),
    ("ACGTACGT", 8),
    ("", 0),
])
def test_sequence_length(input, expected):
    """Test sequence length calculation."""
    assert len(input) == expected
```

## Common Issues

### Import Errors

If you see import errors about Bio modules:

```bash
# Reinstall in development mode
pip install -e .

# Verify correct installation
python -c "import Bio; print(Bio.__file__)"
```

### Working Directory Issues

Always run tests from the `Tests/` directory:

```bash
cd Tests
pytest --offline
```

pytest automatically enforces this via `conftest.py`.

### Skipped Tests

If many tests are skipped:

```bash
# See why tests were skipped
pytest --offline -v -rs

# Run without --offline to include online tests
pytest -v
```

### Collection Errors

If you see collection errors for dependency-related tests, this is expected:

```bash
# These are expected (25 files need optional dependencies)
pytest --offline --collect-only
# 3264 tests collected, 25 errors
```

The errors are from tests requiring:
- Database drivers (MySQLdb, psycopg2, mysql-connector)
- Graphics libraries (reportlab, matplotlib)
- Scientific libraries (scipy, networkx, igraph)
- External tools (PAML, DSSP, NACCESS, PSEA)

## CI/CD Integration

### GitHub Actions

Tests run automatically on push/PR via `.github/workflows/ci.yml`:

```yaml
- name: Run test suite
  run: |
    cd Tests
    pip install -e ..
    coverage run -m pytest --offline -v
    coverage xml
```

### Local CI Simulation

Simulate CI environment locally:

```bash
cd Tests
pip install -e ..
coverage run -m pytest --offline -v
coverage xml
```

## Performance Tips

1. **Use parallel execution** for faster runs:
   ```bash
   pytest --offline -n auto
   ```

2. **Run subset of tests** during development:
   ```bash
   pytest test_SeqIO.py --offline
   ```

3. **Skip slow tests** for quick validation:
   ```bash
   pytest -m "not slow" --offline
   ```

4. **Use --lf** to rerun only failed tests:
   ```bash
   pytest --offline --lf
   ```

## Getting Help

```bash
# Show pytest help
pytest --help

# Show available markers
pytest --markers

# Show available fixtures
pytest --fixtures

# Show test collection
pytest --offline --collect-only
```

## Additional Resources

- Pytest documentation: https://docs.pytest.org/
- Biopython testing docs: http://biopython.org/wiki/Contributing
- Migration plan: See `Tests/migrate_to_pytest.py` for analysis tools
