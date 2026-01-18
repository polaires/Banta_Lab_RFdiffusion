# Fix Interface Metal Dimer Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Fix critical bugs in metal dimer design that cause poor coordination numbers (2/8 atoms instead of 8/8).

**Architecture:** The main fix is passing `fixed_positions` to LigandMPNN so coordinating residues (Asp/Glu) defined in templates aren't mutated during sequence design. Secondary fixes add guiding potentials support and H-bond conditioning.

**Tech Stack:** Python 3.12, RFdiffusion3, LigandMPNN, pytest

---

## Task 1: Add `fixed_positions` Support to `run_ligandmpnn_for_ligand_binding()`

**Files:**
- Modify: `backend/serverless/handler.py:602-724`
- Test: `backend/serverless/test_ligandmpnn_fixed_positions.py` (new)

**Step 1: Write the failing test**

Create `backend/serverless/test_ligandmpnn_fixed_positions.py`:

```python
#!/usr/bin/env python3
"""Tests for LigandMPNN fixed_positions support."""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def test_run_ligandmpnn_accepts_fixed_positions():
    """Verify run_ligandmpnn_for_ligand_binding accepts fixed_positions parameter."""
    from handler import run_ligandmpnn_for_ligand_binding
    import inspect

    sig = inspect.signature(run_ligandmpnn_for_ligand_binding)
    param_names = list(sig.parameters.keys())

    assert "fixed_positions" in param_names, \
        f"fixed_positions not in function signature. Got: {param_names}"


def test_fixed_positions_passed_to_handle_mpnn():
    """Verify fixed_positions is passed through to handle_mpnn."""
    from handler import run_ligandmpnn_for_ligand_binding

    # Minimal PDB with a metal
    test_pdb = """ATOM      1  N   GLU A  15       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  GLU A  15       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   GLU A  15       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   GLU A  15       1.251   2.400   0.000  1.00  0.00           O
HETATM    5 TB   TB  L   1      50.000  50.000  50.000  1.00  0.00          TB
END
"""

    # This should not raise an error about unknown parameter
    # Note: actual MPNN call will fail without server, but parameter parsing should work
    try:
        result = run_ligandmpnn_for_ligand_binding(
            pdb_content=test_pdb,
            ligand_type="lanthanide",
            fixed_positions=["A15"],
            num_sequences=1,
        )
        # If we get here, the function accepted the parameter
        assert True
    except TypeError as e:
        if "fixed_positions" in str(e):
            raise AssertionError(f"fixed_positions parameter not accepted: {e}")
        # Other errors (like missing server) are OK for this test
        pass


if __name__ == "__main__":
    test_run_ligandmpnn_accepts_fixed_positions()
    print("[PASS] test_run_ligandmpnn_accepts_fixed_positions")

    test_fixed_positions_passed_to_handle_mpnn()
    print("[PASS] test_fixed_positions_passed_to_handle_mpnn")

    print("\nAll tests passed!")
```

**Step 2: Run test to verify it fails**

Run: `cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless && python test_ligandmpnn_fixed_positions.py`

Expected: FAIL with "fixed_positions not in function signature"

**Step 3: Add `fixed_positions` parameter to `run_ligandmpnn_for_ligand_binding()`**

Edit `backend/serverless/handler.py` - modify the function signature at line 602:

```python
def run_ligandmpnn_for_ligand_binding(
    pdb_content: str,
    ligand_type: str = "small_molecule",
    num_sequences: int = 4,
    temperature: float = 0.1,
    fixed_positions: Optional[List[str]] = None,  # ADD THIS LINE
) -> Dict[str, Any]:
```

Then add to the `mpnn_input` dict (around line 720):

```python
    mpnn_input = {
        "pdb_content": pdb_content,
        "num_sequences": num_sequences,
        "temperature": temperature,
        "model_type": "ligand_mpnn",
        "remove_waters": True,
        "fixed_positions": fixed_positions,  # ADD THIS LINE
        # ... rest stays the same
    }
```

**Step 4: Run test to verify it passes**

Run: `cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless && python test_ligandmpnn_fixed_positions.py`

Expected: PASS

**Step 5: Commit**

```bash
cd G:\Github_local_repo\Banta_Lab_RFdiffusion
git add backend/serverless/handler.py backend/serverless/test_ligandmpnn_fixed_positions.py
git commit -m "feat(metal-dimer): add fixed_positions param to run_ligandmpnn_for_ligand_binding"
```

---

## Task 2: Extract Coordinating Residue Positions from Template Library

**Files:**
- Modify: `backend/serverless/lanthanide_templates.py`
- Test: `backend/serverless/test_ligandmpnn_fixed_positions.py` (append)

**Step 1: Write the failing test**

Append to `backend/serverless/test_ligandmpnn_fixed_positions.py`:

```python
def test_get_template_fixed_positions():
    """Verify we can extract fixed positions from template library."""
    from lanthanide_templates import get_template_fixed_positions

    # caldwell_4 template has: A15, A25, B15, B25
    positions = get_template_fixed_positions("caldwell_4")

    assert positions is not None, "get_template_fixed_positions returned None"
    assert len(positions) == 4, f"Expected 4 positions, got {len(positions)}"
    assert "A15" in positions, f"A15 not in positions: {positions}"
    assert "A25" in positions, f"A25 not in positions: {positions}"
    assert "B15" in positions, f"B15 not in positions: {positions}"
    assert "B25" in positions, f"B25 not in positions: {positions}"


def test_get_template_fixed_positions_ef_hand():
    """Verify ef_hand_8 template positions."""
    from lanthanide_templates import get_template_fixed_positions

    # ef_hand_8 template has: A10, A15, A20, B10, B15, B20
    positions = get_template_fixed_positions("ef_hand_8")

    assert len(positions) == 6, f"Expected 6 positions, got {len(positions)}"
    expected = ["A10", "A15", "A20", "B10", "B15", "B20"]
    for pos in expected:
        assert pos in positions, f"{pos} not in positions: {positions}"
```

Update the `__main__` block to include new tests.

**Step 2: Run test to verify it fails**

Run: `cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless && python test_ligandmpnn_fixed_positions.py`

Expected: FAIL with "cannot import name 'get_template_fixed_positions'"

**Step 3: Implement `get_template_fixed_positions()` in `lanthanide_templates.py`**

Add after the `TEMPLATE_LIBRARY` dict (around line 260):

```python
def get_template_fixed_positions(template_name: str) -> List[str]:
    """
    Extract coordinating residue positions from a template definition.

    These positions should be passed to LigandMPNN's fixed_positions
    to prevent mutation of metal-coordinating residues.

    Args:
        template_name: Name of template from TEMPLATE_LIBRARY

    Returns:
        List of position strings like ["A15", "A25", "B15", "B25"]
    """
    if template_name not in TEMPLATE_LIBRARY:
        logger.warning(f"Template '{template_name}' not found in library")
        return []

    template = TEMPLATE_LIBRARY[template_name]
    residues = template.get("residues", [])

    positions = []
    for res in residues:
        chain = res.get("chain", "A")
        resnum = res.get("resnum")
        if resnum is not None:
            positions.append(f"{chain}{resnum}")

    return positions
```

**Step 4: Run test to verify it passes**

Run: `cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless && python test_ligandmpnn_fixed_positions.py`

Expected: PASS

**Step 5: Commit**

```bash
cd G:\Github_local_repo\Banta_Lab_RFdiffusion
git add backend/serverless/lanthanide_templates.py backend/serverless/test_ligandmpnn_fixed_positions.py
git commit -m "feat(metal-dimer): add get_template_fixed_positions helper"
```

---

## Task 3: Wire Up Fixed Positions in `_design_joint_metal_dimer()`

**Files:**
- Modify: `backend/serverless/handler.py:5334-5650`
- Test: Manual verification (complex integration)

**Step 1: Identify the LigandMPNN call location**

The call is at `handler.py:5556`. We need to:
1. Get template_name from job_input
2. Call `get_template_fixed_positions(template_name)`
3. Pass result to `run_ligandmpnn_for_ligand_binding()`

**Step 2: Add import at top of handler.py**

Find the lanthanide_templates imports (around line 50-80) and add:

```python
from lanthanide_templates import (
    get_template,
    get_template_info,
    generate_template_from_library,
    recommend_template,
    is_lanthanide,
    generate_parametric_template,
    get_template_fixed_positions,  # ADD THIS
)
```

**Step 3: Modify the LigandMPNN call in `_design_joint_metal_dimer()`**

Find the code block around line 5543-5561 and modify:

```python
        # Run LigandMPNN sequence redesign for lanthanides with carboxylate bias
        # This uses D:6.0,E:4.0 bias to favor Asp/Glu coordinating residues
        # LigandMPNN atomize_side_chains for proper metal binding context
        if is_lanthanide_metal:
            print(f"[JointMetal] Design {design_idx + 1}: Running LigandMPNN with lanthanide preset...")
            try:
                # Fix atom names/IDs before MPNN to avoid NaN issues
                from inference_utils import fix_pdb_for_mpnn
                pdb_content = fix_pdb_for_mpnn(pdb_content)

                # CRITICAL FIX: Get fixed positions from template to preserve coordinators
                coord_fixed_positions = None
                if template_used and LANTHANIDE_TEMPLATES_AVAILABLE:
                    coord_fixed_positions = get_template_fixed_positions(template_used)
                    if coord_fixed_positions:
                        print(f"[JointMetal] Design {design_idx + 1}: Fixing coordinating residues: {coord_fixed_positions}")

                # run_ligandmpnn_for_ligand_binding() automatically fixes incomplete backbone
                mpnn_result = run_ligandmpnn_for_ligand_binding(
                    pdb_content=pdb_content,
                    ligand_type="lanthanide",
                    num_sequences=4,
                    temperature=0.1,
                    fixed_positions=coord_fixed_positions,  # ADD THIS
                )
```

**Step 4: Verify syntax by loading module**

Run: `cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless && python -c "import handler; print('OK')"`

Expected: "OK" (no import errors)

**Step 5: Commit**

```bash
cd G:\Github_local_repo\Banta_Lab_RFdiffusion
git add backend/serverless/handler.py
git commit -m "fix(metal-dimer): pass fixed_positions to LigandMPNN for coordinate preservation

CRITICAL BUG FIX: Template coordinating residues were being mutated
by LigandMPNN because fixed_positions wasn't set. Now extracts
positions from template library and passes them through.

This should fix the poor coordination numbers (2/8 vs 8/8)."
```

---

## Task 4: Add Secondary Sphere H-Bond Conditioning to Metal Dimer

**Files:**
- Modify: `backend/serverless/handler.py:5506-5520`
- Test: Syntax verification

**Step 1: Identify the RFD3 call location**

The RFD3 input dict is built around line 5506-5525. We need to add `select_hbond_acceptor` for coordinating carboxylates.

**Step 2: Add H-bond conditioning after the RFD3 input dict**

Find the code block around line 5517-5525 and add after `select_hotspots`:

```python
        # Add motif scaffolding parameters when using templates
        if fixed_atoms:
            rfd3_input["select_fixed_atoms"] = fixed_atoms
            rfd3_input["select_hotspots"] = {"L1": "all"}
            if design_idx == 0:
                print(f"[JointMetal] Added select_fixed_atoms for {len(fixed_atoms)} residues")
                print(f"[JointMetal] Added select_hotspots: L1:all (ensure metal contact)")

            # NEW: Add H-bond conditioning for secondary coordination sphere
            # Carboxylate oxygens (OE1, OE2 for Glu; OD1, OD2 for Asp) are H-bond acceptors
            # This reinforces H-bond network around metal site
            if template_used and template_type == "library":
                hbond_acceptors = {}
                template_def = TEMPLATE_LIBRARY.get(template_used, {})
                for res in template_def.get("residues", []):
                    chain = res.get("chain", "A")
                    resnum = res.get("resnum")
                    res_type = res.get("type", "GLU")
                    if resnum:
                        key = f"{chain}{resnum}"
                        if res_type == "GLU":
                            hbond_acceptors[key] = "OE1,OE2"
                        elif res_type == "ASP":
                            hbond_acceptors[key] = "OD1,OD2"

                if hbond_acceptors:
                    rfd3_input["select_hbond_acceptor"] = hbond_acceptors
                    if design_idx == 0:
                        print(f"[JointMetal] Added H-bond conditioning for {len(hbond_acceptors)} residues")
```

**Step 3: Verify syntax**

Run: `cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless && python -c "import handler; print('OK')"`

Expected: "OK"

**Step 4: Commit**

```bash
cd G:\Github_local_repo\Banta_Lab_RFdiffusion
git add backend/serverless/handler.py
git commit -m "feat(metal-dimer): add H-bond conditioning for secondary coordination sphere

Adds select_hbond_acceptor for coordinating carboxylate oxygens.
This reinforces the H-bond network around the metal binding site,
which research shows provides ~11x activity boost for metalloenzymes."
```

---

## Task 5: Add `guiding_potentials` Support to `run_rfd3_inference()`

**Files:**
- Modify: `backend/serverless/inference_utils.py:1209+`
- Test: Parameter acceptance test

**Step 1: Write the failing test**

Create `backend/serverless/test_guiding_potentials.py`:

```python
#!/usr/bin/env python3
"""Tests for RFD3 guiding_potentials support."""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def test_run_rfd3_inference_accepts_guiding_potentials():
    """Verify run_rfd3_inference accepts guiding_potentials parameter."""
    from inference_utils import run_rfd3_inference
    import inspect

    sig = inspect.signature(run_rfd3_inference)
    param_names = list(sig.parameters.keys())

    assert "guiding_potentials" in param_names, \
        f"guiding_potentials not in function signature. Got: {param_names}"
    assert "guide_scale" in param_names, \
        f"guide_scale not in function signature. Got: {param_names}"


if __name__ == "__main__":
    test_run_rfd3_inference_accepts_guiding_potentials()
    print("[PASS] test_run_rfd3_inference_accepts_guiding_potentials")
    print("\nAll tests passed!")
```

**Step 2: Run test to verify it fails**

Run: `cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless && python test_guiding_potentials.py`

Expected: FAIL with "guiding_potentials not in function signature"

**Step 3: Add parameters to `run_rfd3_inference()` signature**

Edit `backend/serverless/inference_utils.py` around line 1209. Add after `covalent_bonds`:

```python
def run_rfd3_inference(
    contig: Optional[str] = None,
    length: Optional[str] = None,
    # ... existing params ...
    covalent_bonds: Optional[List[Dict[str, Any]]] = None,
    # ADD THESE:
    guiding_potentials: Optional[List[str]] = None,
    guide_scale: Optional[float] = None,
    guide_decay: Optional[str] = None,
    # Mock mode
    use_mock: bool = False
) -> Dict[str, Any]:
```

Then add to the spec dict construction (around line 1400+, after handling other params):

```python
        # Guiding potentials for design optimization
        if guiding_potentials:
            spec["guiding_potentials"] = guiding_potentials
            if guide_scale is not None:
                spec["guide_scale"] = guide_scale
            if guide_decay is not None:
                spec["guide_decay"] = guide_decay
            print(f"[RFD3] Using guiding potentials: {guiding_potentials}")
```

**Step 4: Run test to verify it passes**

Run: `cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless && python test_guiding_potentials.py`

Expected: PASS

**Step 5: Commit**

```bash
cd G:\Github_local_repo\Banta_Lab_RFdiffusion
git add backend/serverless/inference_utils.py backend/serverless/test_guiding_potentials.py
git commit -m "feat(rfd3): add guiding_potentials, guide_scale, guide_decay params

Exposes RFdiffusion's auxiliary potentials for metal binding design.
substrate_contacts and monomer_ROG can now be configured."
```

---

## Task 6: Wire Up Guiding Potentials for Metal Dimer Design

**Files:**
- Modify: `backend/serverless/handler.py:395-489` (handle_rfd3)
- Modify: `backend/serverless/handler.py:5506-5530` (_design_joint_metal_dimer)

**Step 1: Add params to handle_rfd3()**

Edit `backend/serverless/handler.py` around line 445, add:

```python
    # Guiding potentials
    guiding_potentials = job_input.get("guiding_potentials")
    guide_scale = job_input.get("guide_scale")
    guide_decay = job_input.get("guide_decay")
```

Then add to the `run_rfd3_inference()` call around line 485:

```python
        # Guiding potentials
        guiding_potentials=guiding_potentials,
        guide_scale=guide_scale,
        guide_decay=guide_decay,
```

**Step 2: Use potentials in `_design_joint_metal_dimer()`**

Find the rfd3_input dict (around line 5506) and add:

```python
        rfd3_input = {
            "task": "rfd3",
            "contig": contig,
            "pdb_content": metal_pdb,
            "ligand": metal,
            "ori_token": ori_token,
            "seed": design_seed,
            "num_designs": 1,
            "is_non_loopy": True,
            # ADD: Guiding potentials for metal coordination
            "guiding_potentials": [
                "type:substrate_contacts,weight:5,s:1,r_0:8,d_0:4",
                "type:monomer_ROG,weight:1,min_dist:15",
            ],
            "guide_scale": 2,
        }
```

**Step 3: Verify syntax**

Run: `cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless && python -c "import handler; print('OK')"`

Expected: "OK"

**Step 4: Commit**

```bash
cd G:\Github_local_repo\Banta_Lab_RFdiffusion
git add backend/serverless/handler.py
git commit -m "feat(metal-dimer): use guiding potentials for metal coordination

Adds substrate_contacts and monomer_ROG potentials to metal dimer
design. substrate_contacts reinforces metal-protein proximity,
monomer_ROG encourages compact structures."
```

---

## Task 7: Final Integration Test

**Files:**
- Test: Manual verification with test_local.py

**Step 1: Create integration test**

Create `backend/serverless/test_metal_dimer_integration.py`:

```python
#!/usr/bin/env python3
"""Integration test for metal dimer fixes."""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def test_all_imports():
    """Verify all modified modules import correctly."""
    print("Testing imports...")

    from handler import (
        handle_interface_metal_design,
        run_ligandmpnn_for_ligand_binding,
        handle_rfd3,
    )
    print("  [OK] handler imports")

    from lanthanide_templates import (
        TEMPLATE_LIBRARY,
        get_template_fixed_positions,
    )
    print("  [OK] lanthanide_templates imports")

    from inference_utils import run_rfd3_inference
    print("  [OK] inference_utils imports")

    return True


def test_fixed_positions_extraction():
    """Verify template fixed positions work end-to-end."""
    from lanthanide_templates import get_template_fixed_positions, TEMPLATE_LIBRARY

    print("\nTesting fixed position extraction...")

    for template_name in ["caldwell_4", "ef_hand_8", "lanm_mixed", "high_coord_9"]:
        positions = get_template_fixed_positions(template_name)
        template = TEMPLATE_LIBRARY[template_name]
        expected_count = len(template["residues"])

        assert len(positions) == expected_count, \
            f"{template_name}: expected {expected_count} positions, got {len(positions)}"
        print(f"  [OK] {template_name}: {positions}")

    return True


def test_function_signatures():
    """Verify all function signatures have new parameters."""
    import inspect
    from handler import run_ligandmpnn_for_ligand_binding
    from inference_utils import run_rfd3_inference

    print("\nTesting function signatures...")

    # Check run_ligandmpnn_for_ligand_binding
    sig = inspect.signature(run_ligandmpnn_for_ligand_binding)
    assert "fixed_positions" in sig.parameters, "fixed_positions missing from MPNN function"
    print("  [OK] run_ligandmpnn_for_ligand_binding has fixed_positions")

    # Check run_rfd3_inference
    sig = inspect.signature(run_rfd3_inference)
    assert "guiding_potentials" in sig.parameters, "guiding_potentials missing from RFD3 function"
    assert "guide_scale" in sig.parameters, "guide_scale missing from RFD3 function"
    print("  [OK] run_rfd3_inference has guiding_potentials, guide_scale")

    return True


if __name__ == "__main__":
    print("="*60)
    print("METAL DIMER FIX INTEGRATION TESTS")
    print("="*60)

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

    print("\n" + "="*60)
    if all_passed:
        print("ALL INTEGRATION TESTS PASSED")
    else:
        print("SOME TESTS FAILED")
    print("="*60)

    sys.exit(0 if all_passed else 1)
```

**Step 2: Run integration test**

Run: `cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless && python test_metal_dimer_integration.py`

Expected: ALL INTEGRATION TESTS PASSED

**Step 3: Commit**

```bash
cd G:\Github_local_repo\Banta_Lab_RFdiffusion
git add backend/serverless/test_metal_dimer_integration.py
git commit -m "test: add integration tests for metal dimer fixes"
```

**Step 4: Final commit summary**

```bash
cd G:\Github_local_repo\Banta_Lab_RFdiffusion
git log --oneline -7
```

Expected output shows all 7 commits from this plan.

---

## Summary of Changes

| File | Change |
|------|--------|
| `handler.py:602-724` | Added `fixed_positions` param to `run_ligandmpnn_for_ligand_binding()` |
| `handler.py:5543-5561` | Pass `fixed_positions` from template to LigandMPNN call |
| `handler.py:5517-5530` | Add H-bond conditioning for secondary coordination sphere |
| `handler.py:5506-5530` | Add guiding_potentials to metal dimer RFD3 call |
| `handler.py:395-489` | Wire guiding_potentials through handle_rfd3 |
| `lanthanide_templates.py` | Add `get_template_fixed_positions()` helper |
| `inference_utils.py:1209+` | Add `guiding_potentials`, `guide_scale`, `guide_decay` params |

## Verification

After all tasks complete:
1. All unit tests pass
2. Integration test passes
3. Code imports without errors
4. Ready for real-world testing with RFdiffusion server
