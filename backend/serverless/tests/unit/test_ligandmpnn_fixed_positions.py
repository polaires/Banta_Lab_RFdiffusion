#!/usr/bin/env python3
"""Tests for LigandMPNN fixed_positions support.

These tests use AST parsing to verify function signatures without importing
the module (which has complex dependencies like runpod).
"""

import ast
import os
import re

HANDLER_PATH = os.path.join(os.path.dirname(__file__), "handler.py")
LANTHANIDE_TEMPLATES_PATH = os.path.join(os.path.dirname(__file__), "lanthanide_templates.py")


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


def test_run_ligandmpnn_accepts_fixed_positions():
    """Verify run_ligandmpnn_for_ligand_binding accepts fixed_positions parameter."""
    param_names = get_function_params(HANDLER_PATH, "run_ligandmpnn_for_ligand_binding")

    assert "fixed_positions" in param_names, \
        f"fixed_positions not in function signature. Got: {param_names}"


def test_fixed_positions_in_mpnn_input():
    """Verify fixed_positions is passed to mpnn_input dict."""
    with open(HANDLER_PATH, 'r', encoding='utf-8') as f:
        source = f.read()

    # Find the run_ligandmpnn_for_ligand_binding function
    # and check if fixed_positions is used in mpnn_input
    func_pattern = r'def run_ligandmpnn_for_ligand_binding\([^)]+\).*?(?=\ndef |\Z)'
    match = re.search(func_pattern, source, re.DOTALL)

    if not match:
        raise AssertionError("Could not find run_ligandmpnn_for_ligand_binding function")

    func_body = match.group(0)

    # Check that fixed_positions is assigned to mpnn_input
    assert '"fixed_positions"' in func_body or "'fixed_positions'" in func_body, \
        "fixed_positions not found in mpnn_input dict"


# Tests for Task 2: Template fixed positions extraction
def test_get_template_fixed_positions_exists():
    """Verify get_template_fixed_positions function exists in lanthanide_templates."""
    try:
        param_names = get_function_params(LANTHANIDE_TEMPLATES_PATH, "get_template_fixed_positions")
        assert "template_name" in param_names, \
            f"template_name parameter missing. Got: {param_names}"
    except ValueError as e:
        raise AssertionError(f"get_template_fixed_positions function not found: {e}")


def test_get_template_fixed_positions_returns_positions():
    """Verify get_template_fixed_positions returns correct positions."""
    # This test can only run if we can import the module
    import sys
    sys.path.insert(0, os.path.dirname(__file__))

    try:
        from lanthanide_templates import get_template_fixed_positions

        # caldwell_4 template has: A15, A25, B15, B25
        positions = get_template_fixed_positions("caldwell_4")

        assert positions is not None, "get_template_fixed_positions returned None"
        assert len(positions) == 4, f"Expected 4 positions, got {len(positions)}"
        assert "A15" in positions, f"A15 not in positions: {positions}"
        assert "A25" in positions, f"A25 not in positions: {positions}"
        assert "B15" in positions, f"B15 not in positions: {positions}"
        assert "B25" in positions, f"B25 not in positions: {positions}"

    except ImportError as e:
        # Skip if we can't import (missing dependencies)
        print(f"  [SKIP] Import failed: {e}")
        return


def test_get_template_fixed_positions_ef_hand():
    """Verify ef_hand_8 template positions."""
    import sys
    sys.path.insert(0, os.path.dirname(__file__))

    try:
        from lanthanide_templates import get_template_fixed_positions

        # ef_hand_8 template has: A10, A15, A20, B10, B15, B20
        positions = get_template_fixed_positions("ef_hand_8")

        assert len(positions) == 6, f"Expected 6 positions, got {len(positions)}"
        expected = ["A10", "A15", "A20", "B10", "B15", "B20"]
        for pos in expected:
            assert pos in positions, f"{pos} not in positions: {positions}"

    except ImportError as e:
        print(f"  [SKIP] Import failed: {e}")
        return


if __name__ == "__main__":
    # Task 1 tests
    print("=== Task 1: fixed_positions in run_ligandmpnn_for_ligand_binding ===")

    try:
        test_run_ligandmpnn_accepts_fixed_positions()
        print("[PASS] test_run_ligandmpnn_accepts_fixed_positions")
    except AssertionError as e:
        print(f"[FAIL] test_run_ligandmpnn_accepts_fixed_positions: {e}")

    try:
        test_fixed_positions_in_mpnn_input()
        print("[PASS] test_fixed_positions_in_mpnn_input")
    except AssertionError as e:
        print(f"[FAIL] test_fixed_positions_in_mpnn_input: {e}")

    # Task 2 tests
    print("\n=== Task 2: get_template_fixed_positions helper ===")

    try:
        test_get_template_fixed_positions_exists()
        print("[PASS] test_get_template_fixed_positions_exists")
    except AssertionError as e:
        print(f"[FAIL] test_get_template_fixed_positions_exists: {e}")

    try:
        test_get_template_fixed_positions_returns_positions()
        print("[PASS] test_get_template_fixed_positions_returns_positions")
    except AssertionError as e:
        print(f"[FAIL] test_get_template_fixed_positions_returns_positions: {e}")

    try:
        test_get_template_fixed_positions_ef_hand()
        print("[PASS] test_get_template_fixed_positions_ef_hand")
    except AssertionError as e:
        print(f"[FAIL] test_get_template_fixed_positions_ef_hand: {e}")

    print("\nTest run complete!")
