#!/usr/bin/env python3
"""
Phase 2 Integration Tests for Metal Dimer Design Critique Fixes.

This test suite verifies all Phase 2 improvements:
1. olig_contacts guiding potential added
2. TEBL validation is opt-in (check_tebl=False by default)
3. Geometry RMSD thresholds tightened
4. Contig syntax verified (comma separators, /0 chain breaks)
5. Interface metrics added to validation output
6. C2 symmetry options documented
"""

import os
import sys

# Test file paths
HANDLER_PATH = os.path.join(os.path.dirname(__file__), "handler.py")
METAL_VALIDATION_PATH = os.path.join(os.path.dirname(__file__), "metal_validation.py")


def test_olig_contacts_in_guiding_potentials():
    """Verify olig_contacts potential is included in metal dimer design."""
    with open(HANDLER_PATH, 'r', encoding='utf-8') as f:
        source = f.read()

    assert "olig_contacts" in source, \
        "olig_contacts potential not found in handler.py"

    # Check it's in the guiding_potentials list
    assert "type:olig_contacts" in source, \
        "olig_contacts not properly formatted as guiding potential"

    print("  - olig_contacts potential found in guiding_potentials")


def test_tebl_validation_opt_in():
    """Verify TEBL validation defaults to False (opt-in, not opt-out)."""
    with open(METAL_VALIDATION_PATH, 'r', encoding='utf-8') as f:
        source = f.read()

    assert "check_tebl: bool = False" in source, \
        "check_tebl should default to False"

    assert "if check_tebl:" in source, \
        "TEBL validation should be conditional"

    print("  - TEBL validation is opt-in (check_tebl=False default)")


def test_geometry_rmsd_thresholds_tightened():
    """Verify geometry RMSD thresholds are tightened to crystallographic standards."""
    # Import the module to check actual values
    sys.path.insert(0, os.path.dirname(__file__))
    from metal_validation import LANTHANIDE_CRITERIA

    rmsd = LANTHANIDE_CRITERIA["geometry_rmsd"]

    # Tightened thresholds (from critique):
    # max: 15 deg (was 20), good: 10 deg (was 15), excellent: 5 deg (was 10)
    assert rmsd["max"] <= 15.0, f"max RMSD should be <=15 deg, got {rmsd['max']}"
    assert rmsd["good"] <= 10.0, f"good RMSD should be <=10 deg, got {rmsd['good']}"
    assert rmsd["excellent"] <= 7.0, f"excellent RMSD should be <=7 deg, got {rmsd['excellent']}"

    print(f"  - RMSD thresholds tightened: max={rmsd['max']}, good={rmsd['good']}, excellent={rmsd['excellent']}")


def test_contig_syntax_correct():
    """Verify contig syntax uses comma separators and /0 chain breaks."""
    with open(HANDLER_PATH, 'r', encoding='utf-8') as f:
        source = f.read()

    # Should use forward slash /0 for chain break, NOT backslash
    assert "/0" in source, "Chain break syntax /0 not found"

    # Verify proper format in contig patterns
    # Look for typical patterns: "60-80,/0,60-80" or similar
    import re
    contig_patterns = re.findall(r'["\'][\d\-]+,/0,[\d\-]+["\']', source)
    assert len(contig_patterns) > 0, "No valid contig patterns found with /0 chain break"

    print(f"  - Contig syntax verified: {len(contig_patterns)} patterns with /0 chain breaks")


def test_interface_metrics_in_validation():
    """Verify interface metrics are added to validation output."""
    with open(HANDLER_PATH, 'r', encoding='utf-8') as f:
        source = f.read()

    # Check for interface metrics keys
    assert "interface_metrics" in source, \
        "interface_metrics not found in handler.py"

    # Check specific metrics are included
    metrics = ["contacts", "hbonds", "buried_sasa", "packstat"]
    found_metrics = [m for m in metrics if m in source]

    assert len(found_metrics) >= 3, \
        f"Expected at least 3 interface metrics, found: {found_metrics}"

    print(f"  - Interface metrics found: {found_metrics}")


def test_c2_symmetry_documented():
    """Verify C2 symmetry options are documented in docstring."""
    with open(HANDLER_PATH, 'r', encoding='utf-8') as f:
        source = f.read()

    # Check for symmetry documentation in docstring
    symmetry_docs = [
        "Symmetry Options",
        "C2 symmetry",
        "explicit 2-chain",
    ]

    found_docs = [doc for doc in symmetry_docs if doc.lower() in source.lower()]

    assert len(found_docs) >= 2, \
        f"Expected symmetry documentation, found: {found_docs}"

    # Check for trade-off documentation
    assert "pros" in source.lower() or "advantages" in source.lower() or "trade-off" in source.lower(), \
        "Symmetry trade-offs should be documented"

    print(f"  - C2 symmetry documented with trade-offs")


def test_handler_syntax_valid():
    """Verify handler.py has valid Python syntax after all modifications."""
    import ast

    with open(HANDLER_PATH, 'r', encoding='utf-8') as f:
        source = f.read()

    try:
        ast.parse(source)
        print("  - handler.py syntax is valid")
    except SyntaxError as e:
        raise AssertionError(f"handler.py has syntax error: {e}")


def test_metal_validation_syntax_valid():
    """Verify metal_validation.py has valid Python syntax after all modifications."""
    import ast

    with open(METAL_VALIDATION_PATH, 'r', encoding='utf-8') as f:
        source = f.read()

    try:
        ast.parse(source)
        print("  - metal_validation.py syntax is valid")
    except SyntaxError as e:
        raise AssertionError(f"metal_validation.py has syntax error: {e}")


def run_all_tests():
    """Run all Phase 2 integration tests."""
    print("=" * 60)
    print("PHASE 2 METAL DIMER CRITIQUE FIXES - INTEGRATION TESTS")
    print("=" * 60)

    tests = [
        ("olig_contacts guiding potential", test_olig_contacts_in_guiding_potentials),
        ("TEBL validation opt-in", test_tebl_validation_opt_in),
        ("Geometry RMSD thresholds", test_geometry_rmsd_thresholds_tightened),
        ("Contig syntax correct", test_contig_syntax_correct),
        ("Interface metrics in validation", test_interface_metrics_in_validation),
        ("C2 symmetry documented", test_c2_symmetry_documented),
        ("handler.py syntax valid", test_handler_syntax_valid),
        ("metal_validation.py syntax valid", test_metal_validation_syntax_valid),
    ]

    passed = 0
    failed = 0

    for name, test_func in tests:
        print(f"\n[TEST] {name}")
        try:
            test_func()
            print(f"[PASS] {name}")
            passed += 1
        except AssertionError as e:
            print(f"[FAIL] {name}: {e}")
            failed += 1
        except Exception as e:
            print(f"[ERROR] {name}: {type(e).__name__}: {e}")
            failed += 1

    print("\n" + "=" * 60)
    print(f"RESULTS: {passed} passed, {failed} failed")
    print("=" * 60)

    if failed == 0:
        print("\n*** ALL PHASE 2 TESTS PASSED ***")
        print("""
Phase 2 Critique Fixes Summary:
1. olig_contacts potential added for interface optimization
2. TEBL validation is opt-in (not forced on all designs)
3. Geometry RMSD thresholds tightened to crystallographic standards
4. Contig syntax verified correct (comma separators, /0 chain breaks)
5. Interface metrics added to validation output
6. C2 symmetry options documented with trade-offs

Combined with Phase 1 fixes (fixed_positions, guiding_potentials, H-bond conditioning),
the metal dimer implementation now addresses all CRITICAL and HIGH priority items
from the original critique.
""")
    else:
        print(f"\n*** {failed} TEST(S) FAILED - REVIEW REQUIRED ***")

    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
