#!/usr/bin/env python3
"""Tests for RFD3 guiding_potentials support.

Uses AST parsing to verify function signatures without importing
the module (which has complex dependencies).
"""

import ast
import os

INFERENCE_UTILS_PATH = os.path.join(os.path.dirname(__file__), "inference_utils.py")


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


def test_run_rfd3_inference_accepts_guiding_potentials():
    """Verify run_rfd3_inference accepts guiding_potentials parameter."""
    param_names = get_function_params(INFERENCE_UTILS_PATH, "run_rfd3_inference")

    assert "guiding_potentials" in param_names, \
        f"guiding_potentials not in function signature. Got: {param_names}"
    assert "guide_scale" in param_names, \
        f"guide_scale not in function signature. Got: {param_names}"


def test_guiding_potentials_in_spec():
    """Verify guiding_potentials is used in the spec dict."""
    with open(INFERENCE_UTILS_PATH, 'r', encoding='utf-8') as f:
        source = f.read()

    # Check that guiding_potentials is referenced in the function body
    assert '"guiding_potentials"' in source or "'guiding_potentials'" in source, \
        "guiding_potentials not found in inference_utils.py"


HANDLER_PATH = os.path.join(os.path.dirname(__file__), "handler.py")


def test_olig_contacts_in_metal_dimer():
    """Verify olig_contacts potential is used in metal dimer design."""
    with open(HANDLER_PATH, 'r', encoding='utf-8') as f:
        source = f.read()

    # Check that olig_contacts is in guiding_potentials
    assert "olig_contacts" in source, \
        "olig_contacts potential not found in handler.py"


if __name__ == "__main__":
    print("=== Guiding Potentials Tests ===")

    try:
        test_run_rfd3_inference_accepts_guiding_potentials()
        print("[PASS] test_run_rfd3_inference_accepts_guiding_potentials")
    except AssertionError as e:
        print(f"[FAIL] test_run_rfd3_inference_accepts_guiding_potentials: {e}")

    try:
        test_guiding_potentials_in_spec()
        print("[PASS] test_guiding_potentials_in_spec")
    except AssertionError as e:
        print(f"[FAIL] test_guiding_potentials_in_spec: {e}")

    try:
        test_olig_contacts_in_metal_dimer()
        print("[PASS] test_olig_contacts_in_metal_dimer")
    except AssertionError as e:
        print(f"[FAIL] test_olig_contacts_in_metal_dimer: {e}")

    print("\nTest run complete!")
