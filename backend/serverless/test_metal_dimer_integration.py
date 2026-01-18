#!/usr/bin/env python3
"""Integration test for metal dimer fixes.

Verifies all changes from the fix-interface-metal-dimer implementation plan.
"""

import ast
import os
import sys

# Setup path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

HANDLER_PATH = os.path.join(os.path.dirname(__file__), "handler.py")
LANTHANIDE_TEMPLATES_PATH = os.path.join(os.path.dirname(__file__), "lanthanide_templates.py")
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


def test_all_imports():
    """Verify all modified modules have valid syntax."""
    print("Testing module syntax...")

    # Test handler.py
    with open(HANDLER_PATH, 'r', encoding='utf-8') as f:
        ast.parse(f.read())
    print("  [OK] handler.py syntax valid")

    # Test lanthanide_templates.py
    with open(LANTHANIDE_TEMPLATES_PATH, 'r', encoding='utf-8') as f:
        ast.parse(f.read())
    print("  [OK] lanthanide_templates.py syntax valid")

    # Test inference_utils.py
    with open(INFERENCE_UTILS_PATH, 'r', encoding='utf-8') as f:
        ast.parse(f.read())
    print("  [OK] inference_utils.py syntax valid")

    return True


def test_fixed_positions_extraction():
    """Verify template fixed positions work end-to-end."""
    from lanthanide_templates import get_template_fixed_positions, TEMPLATE_LIBRARY

    print("\nTesting fixed position extraction...")

    templates_tested = 0
    for template_name in ["caldwell_4", "ef_hand_8", "lanm_mixed", "high_coord_9"]:
        if template_name not in TEMPLATE_LIBRARY:
            print(f"  [SKIP] {template_name} not in TEMPLATE_LIBRARY")
            continue

        positions = get_template_fixed_positions(template_name)
        template = TEMPLATE_LIBRARY[template_name]
        expected_count = len(template["residues"])

        assert len(positions) == expected_count, \
            f"{template_name}: expected {expected_count} positions, got {len(positions)}"
        print(f"  [OK] {template_name}: {positions}")
        templates_tested += 1

    assert templates_tested >= 2, f"Only tested {templates_tested} templates"
    return True


def test_function_signatures():
    """Verify all function signatures have new parameters."""
    print("\nTesting function signatures...")

    # Check run_ligandmpnn_for_ligand_binding
    params = get_function_params(HANDLER_PATH, "run_ligandmpnn_for_ligand_binding")
    assert "fixed_positions" in params, "fixed_positions missing from MPNN function"
    print("  [OK] run_ligandmpnn_for_ligand_binding has fixed_positions")

    # Check run_rfd3_inference
    params = get_function_params(INFERENCE_UTILS_PATH, "run_rfd3_inference")
    assert "guiding_potentials" in params, "guiding_potentials missing from RFD3 function"
    assert "guide_scale" in params, "guide_scale missing from RFD3 function"
    assert "guide_decay" in params, "guide_decay missing from RFD3 function"
    print("  [OK] run_rfd3_inference has guiding_potentials, guide_scale, guide_decay")

    # Check get_template_fixed_positions exists
    params = get_function_params(LANTHANIDE_TEMPLATES_PATH, "get_template_fixed_positions")
    assert "template_name" in params, "template_name missing from get_template_fixed_positions"
    print("  [OK] get_template_fixed_positions exists with template_name param")

    return True


def test_handler_integration():
    """Verify handler.py has the integration code."""
    print("\nTesting handler.py integration code...")

    with open(HANDLER_PATH, 'r', encoding='utf-8') as f:
        source = f.read()

    # Check fixed_positions is passed to LigandMPNN
    assert "fixed_positions=coord_fixed_positions" in source, \
        "fixed_positions not passed to run_ligandmpnn_for_ligand_binding"
    print("  [OK] fixed_positions passed to LigandMPNN call")

    # Check get_template_fixed_positions is imported
    assert "get_template_fixed_positions" in source, \
        "get_template_fixed_positions not imported in handler.py"
    print("  [OK] get_template_fixed_positions imported")

    # Check H-bond conditioning is added
    assert "select_hbond_acceptor" in source, \
        "select_hbond_acceptor not found in handler.py"
    print("  [OK] H-bond conditioning (select_hbond_acceptor) present")

    # Check guiding_potentials is in rfd3_input
    assert "guiding_potentials" in source, \
        "guiding_potentials not found in handler.py"
    print("  [OK] guiding_potentials present")

    # Check substrate_contacts potential is used
    assert "substrate_contacts" in source, \
        "substrate_contacts potential not found in handler.py"
    print("  [OK] substrate_contacts potential configured")

    return True


def test_inference_utils_integration():
    """Verify inference_utils.py has guiding potentials code."""
    print("\nTesting inference_utils.py integration code...")

    with open(INFERENCE_UTILS_PATH, 'r', encoding='utf-8') as f:
        source = f.read()

    # Check guiding_potentials is added to spec
    assert 'spec["guiding_potentials"]' in source, \
        "guiding_potentials not added to spec dict"
    print("  [OK] guiding_potentials added to spec dict")

    # Check guide_scale is added to spec
    assert 'spec["guide_scale"]' in source, \
        "guide_scale not added to spec dict"
    print("  [OK] guide_scale added to spec dict")

    return True


if __name__ == "__main__":
    print("=" * 60)
    print("METAL DIMER FIX INTEGRATION TESTS")
    print("=" * 60)

    all_passed = True

    try:
        test_all_imports()
    except Exception as e:
        print(f"  [FAIL] Import test failed: {e}")
        all_passed = False

    try:
        test_fixed_positions_extraction()
    except Exception as e:
        print(f"  [FAIL] Fixed positions test failed: {e}")
        all_passed = False

    try:
        test_function_signatures()
    except Exception as e:
        print(f"  [FAIL] Signature test failed: {e}")
        all_passed = False

    try:
        test_handler_integration()
    except Exception as e:
        print(f"  [FAIL] Handler integration test failed: {e}")
        all_passed = False

    try:
        test_inference_utils_integration()
    except Exception as e:
        print(f"  [FAIL] Inference utils integration test failed: {e}")
        all_passed = False

    print("\n" + "=" * 60)
    if all_passed:
        print("ALL INTEGRATION TESTS PASSED")
    else:
        print("SOME TESTS FAILED")
    print("=" * 60)

    sys.exit(0 if all_passed else 1)
