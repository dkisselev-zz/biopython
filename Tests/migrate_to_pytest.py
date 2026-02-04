#!/usr/bin/env python
"""Semi-automated migration helper for converting tests to pytest.

This script provides utilities to:
1. Analyze test files for migration opportunities
2. Detect dependency checks (MissingExternalDependencyError)
3. Detect internet usage (requires_internet)
4. Generate migration reports
5. Identify test complexity (base classes, setUp/tearDown, etc.)

Usage:
    python migrate_to_pytest.py --analyze test_Seq.py
    python migrate_to_pytest.py --report-all
    python migrate_to_pytest.py --list-simple
"""

import ast
import sys
import argparse
from pathlib import Path
from typing import List, Tuple, Dict, Set


class MigrationAnalyzer(ast.NodeVisitor):
    """Analyze test file for migration opportunities."""

    def __init__(self):
        self.needs_pytest_import = False
        self.has_missing_dependency_check = False
        self.has_requires_internet = False
        self.online_test_methods = []
        self.external_tool_checks = []
        self.import_statements = []
        self.base_classes = set()
        self.test_classes = []
        self.test_methods = []
        self.has_setup = False
        self.has_teardown = False
        self.has_setupclass = False
        self.has_teardownclass = False

    def visit_Import(self, node):
        for alias in node.names:
            self.import_statements.append(alias.name)
            if alias.name == "requires_internet":
                self.has_requires_internet = True
        self.generic_visit(node)

    def visit_ImportFrom(self, node):
        if node.module == "Bio" and any(
            alias.name in ("MissingExternalDependencyError", "MissingPythonDependencyError")
            for alias in node.names
        ):
            self.has_missing_dependency_check = True
        self.generic_visit(node)

    def visit_ClassDef(self, node):
        # Check if it's a test class
        if node.name.startswith("Test") or "TestCase" in [b.id if isinstance(b, ast.Name) else str(b) for b in node.bases]:
            self.test_classes.append(node.name)

            # Check for base classes
            for base in node.bases:
                if isinstance(base, ast.Name):
                    self.base_classes.add(base.id)
                elif isinstance(base, ast.Attribute):
                    self.base_classes.add(f"{base.value.id if isinstance(base.value, ast.Name) else ''}.{base.attr}")

        self.generic_visit(node)

    def visit_FunctionDef(self, node):
        # Check if it's a test method
        if node.name.startswith("test_"):
            self.test_methods.append(node.name)

        # Check for setUp/tearDown
        if node.name == "setUp":
            self.has_setup = True
        elif node.name == "tearDown":
            self.has_teardown = True
        elif node.name == "setUpClass":
            self.has_setupclass = True
        elif node.name == "tearDownClass":
            self.has_teardownclass = True

        # Check if test calls requires_internet.check()
        for child in ast.walk(node):
            if isinstance(child, ast.Call):
                if isinstance(child.func, ast.Attribute):
                    if (isinstance(child.func.value, ast.Name) and
                        child.func.value.id == "requires_internet" and
                        child.func.attr == "check"):
                        self.online_test_methods.append(node.name)

        self.generic_visit(node)

    def visit_Raise(self, node):
        """Detect raises of MissingExternalDependencyError."""
        if isinstance(node.exc, ast.Call):
            if isinstance(node.exc.func, ast.Name):
                if node.exc.func.id == "MissingExternalDependencyError":
                    self.has_missing_dependency_check = True
        self.generic_visit(node)


def analyze_file(filepath: Path) -> Tuple[MigrationAnalyzer, str]:
    """Analyze a test file for migration opportunities."""
    with open(filepath) as f:
        content = f.read()

    try:
        tree = ast.parse(content)
        analyzer = MigrationAnalyzer()
        analyzer.visit(tree)
        return analyzer, content
    except SyntaxError as e:
        print(f"Syntax error in {filepath}: {e}")
        return None, content


def calculate_complexity_score(analyzer: MigrationAnalyzer) -> int:
    """Calculate complexity score for a test file.

    Lower score = simpler to migrate
    Higher score = more complex
    """
    score = 0

    # Base classes add complexity
    if "unittest.TestCase" in analyzer.base_classes:
        score += 0  # Standard, no extra complexity
    else:
        for base in analyzer.base_classes:
            if base not in ["unittest.TestCase", "TestCase"]:
                score += 10  # Custom base class

    # setUp/tearDown adds moderate complexity
    if analyzer.has_setup or analyzer.has_teardown:
        score += 2

    # setUpClass/tearDownClass adds more complexity
    if analyzer.has_setupclass or analyzer.has_teardownclass:
        score += 5

    # External dependencies add complexity
    if analyzer.has_missing_dependency_check:
        score += 3

    # Internet tests are moderate complexity
    if analyzer.has_requires_internet:
        score += 2

    return score


def generate_migration_report(filepath: Path) -> str:
    """Generate a migration report for a test file."""
    result = analyze_file(filepath)
    if result[0] is None:
        return f"ERROR: Could not parse {filepath.name}"

    analyzer, content = result

    lines = []
    lines.append(f"\n{'='*70}")
    lines.append(f"Migration Report: {filepath.name}")
    lines.append(f"{'='*70}")

    complexity = calculate_complexity_score(analyzer)
    lines.append(f"\nComplexity Score: {complexity} ({'Simple' if complexity < 5 else 'Moderate' if complexity < 15 else 'Complex'})")

    lines.append(f"\nStatistics:")
    lines.append(f"  Test classes: {len(analyzer.test_classes)}")
    lines.append(f"  Test methods: {len(analyzer.test_methods)}")

    if analyzer.base_classes:
        lines.append(f"\nBase Classes:")
        for base in analyzer.base_classes:
            lines.append(f"  - {base}")

    if analyzer.has_setup or analyzer.has_teardown or analyzer.has_setupclass or analyzer.has_teardownclass:
        lines.append(f"\nFixtures:")
        if analyzer.has_setup:
            lines.append("  ✓ Has setUp()")
        if analyzer.has_teardown:
            lines.append("  ✓ Has tearDown()")
        if analyzer.has_setupclass:
            lines.append("  ✓ Has setUpClass()")
        if analyzer.has_teardownclass:
            lines.append("  ✓ Has tearDownClass()")

    if analyzer.has_missing_dependency_check:
        lines.append("\n⚠️  Has dependency checks (consider pytest.importorskip)")

    if analyzer.has_requires_internet:
        lines.append(f"\n⚠️  Has internet dependency (add @pytest.mark.online)")
        if analyzer.online_test_methods:
            lines.append(f"   Methods: {', '.join(analyzer.online_test_methods)}")

    lines.append("\nMigration Steps:")
    steps = []

    if not analyzer.has_missing_dependency_check and not analyzer.has_requires_internet:
        steps.append("1. No changes needed - already compatible!")
    else:
        step_num = 1
        if analyzer.has_requires_internet:
            steps.append(f"{step_num}. Add @pytest.mark.online to test classes")
            step_num += 1
        if analyzer.has_missing_dependency_check:
            steps.append(f"{step_num}. Convert module-level dependency checks to pytest.importorskip")
            step_num += 1
        steps.append(f"{step_num}. Test with: pytest {filepath.name} -v")
        step_num += 1
        steps.append(f"{step_num}. Verify with unittest: python {filepath.name}")

    for step in steps:
        lines.append(f"  {step}")

    return "\n".join(lines)


def list_simple_tests(test_dir: Path) -> List[Tuple[Path, int]]:
    """List test files sorted by complexity (simplest first)."""
    test_files = sorted(test_dir.glob("test_*.py"))

    results = []
    for filepath in test_files:
        result = analyze_file(filepath)
        if result[0] is None:
            continue

        analyzer, _ = result
        complexity = calculate_complexity_score(analyzer)
        results.append((filepath, complexity))

    # Sort by complexity (simplest first)
    results.sort(key=lambda x: x[1])
    return results


def report_all_tests(test_dir: Path):
    """Generate a summary report for all test files."""
    results = list_simple_tests(test_dir)

    print(f"\n{'='*70}")
    print(f"Test Files Migration Priority Report")
    print(f"{'='*70}")
    print(f"\nTotal test files: {len(results)}")

    # Group by complexity
    simple = [r for r in results if r[1] < 5]
    moderate = [r for r in results if 5 <= r[1] < 15]
    complex_tests = [r for r in results if r[1] >= 15]

    print(f"\nSimple tests (complexity < 5): {len(simple)}")
    print(f"Moderate tests (5 <= complexity < 15): {len(moderate)}")
    print(f"Complex tests (complexity >= 15): {len(complex_tests)}")

    print(f"\n{'='*70}")
    print("Top 20 Simplest Tests (Priority for Initial Migration)")
    print(f"{'='*70}")

    for i, (filepath, complexity) in enumerate(results[:20], 1):
        print(f"{i:2d}. {filepath.name:40s} (complexity: {complexity:2d})")

    print(f"\n{'='*70}")
    print("Most Complex Tests (Migrate Last)")
    print(f"{'='*70}")

    for i, (filepath, complexity) in enumerate(results[-10:], 1):
        print(f"{i:2d}. {filepath.name:40s} (complexity: {complexity:2d})")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analyze", metavar="FILE", help="Analyze a specific test file")
    parser.add_argument("--report-all", action="store_true", help="Generate summary report for all tests")
    parser.add_argument("--list-simple", action="store_true", help="List simplest tests first")

    args = parser.parse_args()

    test_dir = Path(__file__).parent

    if args.analyze:
        filepath = test_dir / args.analyze
        if not filepath.exists():
            print(f"Error: {filepath} not found")
            sys.exit(1)
        print(generate_migration_report(filepath))

    elif args.report_all:
        report_all_tests(test_dir)

    elif args.list_simple:
        results = list_simple_tests(test_dir)
        print(f"\nTest files sorted by complexity (simplest first):\n")
        for filepath, complexity in results[:30]:
            status = "✓ Simple" if complexity < 5 else "◐ Moderate" if complexity < 15 else "⚠ Complex"
            print(f"{status:10s} [{complexity:2d}] {filepath.name}")

    else:
        parser.print_help()


if __name__ == "__main__":
    main()
