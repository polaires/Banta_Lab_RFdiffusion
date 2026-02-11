#!/usr/bin/env python3
"""Tests for validation options."""

import ast
import os

METAL_VALIDATION_PATH = os.path.join(os.path.dirname(__file__), "metal_validation.py")


def get_function_params(filepath: str, func_name: str) -> list:
    """Extract parameter names from a function definition using AST."""
    with open(filepath, 'r', encoding='utf-8') as f:
        source = f.read()

    tree = ast.parse(source)

    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == func_name:
            params = []
            for arg in node.args.args:
                params.append(arg.arg)
            for arg in node.args.kwonlyargs:
                params.append(arg.arg)
            return params

    raise ValueError(f"Function {func_name} not found in {filepath}")


def test_validate_lanthanide_has_check_tebl_default_false():
    """Verify check_tebl defaults to False for general validation."""
    with open(METAL_VALIDATION_PATH, 'r', encoding='utf-8') as f:
        source = f.read()

    # check_tebl should default to False
    assert "check_tebl: bool = False" in source, \
        "check_tebl should default to False"


def test_tebl_validation_is_conditional():
    """Verify TEBL validation only runs when check_tebl=True."""
    with open(METAL_VALIDATION_PATH, 'r', encoding='utf-8') as f:
        source = f.read()

    # TEBL validation should be inside if check_tebl:
    assert "if check_tebl:" in source, \
        "TEBL validation should be conditional on check_tebl flag"


def test_geometry_rmsd_thresholds_tightened():
    """Verify geometry RMSD thresholds are tightened per critique."""
    import sys
    sys.path.insert(0, os.path.dirname(__file__))

    from metal_validation import LANTHANIDE_CRITERIA

    rmsd = LANTHANIDE_CRITERIA["geometry_rmsd"]

    # Tightened thresholds (from critique):
    # max: 15 deg (was 20), good: 10 deg (was 15), excellent: 5 deg (was 10)
    assert rmsd["max"] <= 15.0, f"max RMSD should be <=15 deg, got {rmsd['max']}"
    assert rmsd["good"] <= 10.0, f"good RMSD should be <=10 deg, got {rmsd['good']}"
    assert rmsd["excellent"] <= 7.0, f"excellent RMSD should be <=7 deg, got {rmsd['excellent']}"


if __name__ == "__main__":
    print("=== Validation Options Tests ===")

    try:
        test_validate_lanthanide_has_check_tebl_default_false()
        print("[PASS] test_validate_lanthanide_has_check_tebl_default_false")
    except AssertionError as e:
        print(f"[FAIL] test_validate_lanthanide_has_check_tebl_default_false: {e}")

    try:
        test_tebl_validation_is_conditional()
        print("[PASS] test_tebl_validation_is_conditional")
    except AssertionError as e:
        print(f"[FAIL] test_tebl_validation_is_conditional: {e}")

    try:
        test_geometry_rmsd_thresholds_tightened()
        print("[PASS] test_geometry_rmsd_thresholds_tightened")
    except AssertionError as e:
        print(f"[FAIL] test_geometry_rmsd_thresholds_tightened: {e}")

    print("\nTest run complete!")
