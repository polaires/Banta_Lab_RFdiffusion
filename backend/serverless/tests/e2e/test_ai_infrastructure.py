"""
Test script for AI Infrastructure Enhancement
Tests the original query that prompted this enhancement.
"""
import sys
import json

# Test query
TEST_QUERY = "scaffold the active pocket pqq-ca of 4cvb and keep all interacting residues and make protein more stable, and i want it able to preserve its ADHD functionality."

def test_nl_parser():
    """Test NL parser with typo correction and enzyme detection."""
    print("=" * 70)
    print("TEST 1: NL Parser - Typo Correction & Enzyme Detection")
    print("=" * 70)
    print(f"\nQuery: {TEST_QUERY}\n")

    from nl_design_parser import SimpleFallbackParser

    parser = SimpleFallbackParser()
    intent = parser.parse(TEST_QUERY)

    print("Parsed DesignIntent:")
    print(f"  - source_pdb_id: {intent.source_pdb_id}")
    print(f"  - ligand_name: {intent.ligand_name}")
    print(f"  - metal_type: {intent.metal_type}")
    print(f"  - is_scaffolding: {intent.is_scaffolding}")
    print(f"  - stability_focus: {intent.stability_focus}")
    print(f"  - include_all_contacts: {intent.include_all_contacts}")
    print(f"  - enzyme_class: {intent.enzyme_class}")
    print(f"  - enzyme_class_confidence: {intent.enzyme_class_confidence}")
    print(f"  - preserve_function: {intent.preserve_function}")
    print(f"  - typo_corrections: {intent.typo_corrections}")
    print(f"  - corrected_query: {intent.corrected_query[:80]}..." if intent.corrected_query else "  - corrected_query: (none)")

    # Verify key detections
    checks = []

    if intent.source_pdb_id and intent.source_pdb_id.upper() == "4CVB":
        checks.append("[PASS] PDB ID detected: 4CVB")
    else:
        checks.append(f"[FAIL] PDB ID not detected (got: {intent.source_pdb_id})")

    if intent.ligand_name and "pqq" in intent.ligand_name.lower():
        checks.append("[PASS] Ligand detected: PQQ")
    else:
        checks.append(f"[FAIL] Ligand not detected (got: {intent.ligand_name})")

    if intent.metal_type and intent.metal_type.upper() == "CA":
        checks.append("[PASS] Metal detected: CA")
    else:
        checks.append(f"[FAIL] Metal not detected (got: {intent.metal_type})")

    if intent.stability_focus:
        checks.append("[PASS] Stability focus detected")
    else:
        checks.append("[FAIL] Stability focus not detected")

    if intent.enzyme_class:
        checks.append(f"[PASS] Enzyme class detected: {intent.enzyme_class}")
    else:
        checks.append("[FAIL] Enzyme class not detected")

    if intent.preserve_function:
        checks.append("[PASS] Preserve function intent detected")
    else:
        checks.append("[FAIL] Preserve function not detected")

    if intent.typo_corrections and any("adhd" in t.lower() for t in intent.typo_corrections):
        checks.append("[PASS] Typo correction applied: ADHD -> ADH")
    elif intent.typo_corrections:
        checks.append(f"[INFO] Typo corrections: {intent.typo_corrections}")
    else:
        checks.append("[INFO] No typo corrections needed or applied")

    print("\nVerification:")
    for check in checks:
        print(f"  {check}")

    return intent


def test_stability_profile(intent):
    """Test stability profile selection."""
    print("\n" + "=" * 70)
    print("TEST 2: Stability Profile Selection")
    print("=" * 70)

    from design_rules import select_stability_profile, STABILITY_PROFILES

    print("\nAvailable profiles:")
    for name, profile in STABILITY_PROFILES.items():
        print(f"  - {name}: cfg_scale={profile.cfg_scale}, step_scale={profile.step_scale}, "
              f"gamma_0={profile.gamma_0}, num_timesteps={profile.num_timesteps}")

    stability = select_stability_profile(intent, intent.design_goal or "binding")

    print(f"\nSelected profile for stability_focus={intent.stability_focus}:")
    print(f"  - name: {stability.name}")
    print(f"  - cfg_scale: {stability.cfg_scale}")
    print(f"  - step_scale: {stability.step_scale}")
    print(f"  - gamma_0: {stability.gamma_0}")
    print(f"  - num_timesteps: {stability.num_timesteps}")
    print(f"  - plddt_threshold: {stability.plddt_threshold}")

    # Verify stability-focused profile is selected
    if intent.stability_focus and stability.name != "balanced":
        print(f"\n[PASS] Stability-focused profile selected: {stability.name}")
    elif not intent.stability_focus and stability.name == "balanced":
        print("\n[PASS] Balanced profile selected (stability_focus=False)")
    else:
        print(f"\n[INFO] Profile: {stability.name}")

    return stability


def test_enzyme_detection():
    """Test enzyme class database and detection."""
    print("\n" + "=" * 70)
    print("TEST 3: Enzyme Class Detection")
    print("=" * 70)

    from enzyme_chemistry import detect_enzyme_class, get_preservation_requirements, ENZYME_CLASS_DATABASE

    print(f"\nEnzyme database has {len(ENZYME_CLASS_DATABASE)} classes")

    # Test detection with various queries
    test_queries = [
        ("ADHD functionality", "pqq"),  # Typo
        ("ADH functionality", "pqq"),   # Correct
        ("dehydrogenase", None),
        ("quinoprotein", "pqq"),
        ("preserve catalytic activity", None),
    ]

    print("\nDetection tests:")
    for query, ligand in test_queries:
        enzyme_class, confidence = detect_enzyme_class(query, ligand)
        print(f"  Query: '{query}' (ligand={ligand})")
        print(f"    -> enzyme_class={enzyme_class}, confidence={confidence:.2f}")

    # Get preservation requirements for quinoprotein (PQQ-dependent)
    print("\nPreservation requirements for 'quinoprotein':")
    reqs = get_preservation_requirements("quinoprotein")
    if reqs:
        for key, value in reqs.items():
            print(f"  - {key}: {value}")
    else:
        print("  (no requirements found)")

    return True


def test_rfd3_config():
    """Test RFD3 config generation with new parameters."""
    print("\n" + "=" * 70)
    print("TEST 4: RFD3 Config Generation")
    print("=" * 70)

    from rfd3_config_generator import RFD3Config

    # Create a config with new stability parameters
    config = RFD3Config(
        contig="X1,L1,/0,100-130",
        cfg_scale=3.0,
        step_scale=1.2,
        gamma_0=0.4,
        num_timesteps=250,
        select_exposed={"L1": "O1,O2,C3"},
        select_buried={"L1": "C4,C5,C6"},
    )

    api_params = config.to_api_params()

    print("\nRFD3 Config API params:")
    for key, value in api_params.items():
        if key != "pdb_content":
            print(f"  - {key}: {value}")

    # Verify new parameters are included
    checks = []
    if "gamma_0" in api_params:
        checks.append(f"[PASS] gamma_0 included: {api_params['gamma_0']}")
    else:
        checks.append("[FAIL] gamma_0 not in API params")

    if "select_exposed" in api_params:
        checks.append(f"[PASS] select_exposed included: {api_params['select_exposed']}")
    else:
        checks.append("[FAIL] select_exposed not in API params")

    if api_params.get("num_timesteps") == 250:
        checks.append("[PASS] num_timesteps=250 (non-default)")

    print("\nVerification:")
    for check in checks:
        print(f"  {check}")

    return True


def test_unified_analyzer():
    """Test unified analyzer enzyme preservation methods."""
    print("\n" + "=" * 70)
    print("TEST 5: Unified Analyzer - Enzyme Preservation")
    print("=" * 70)

    from unified_analyzer import UnifiedDesignAnalyzer

    analyzer = UnifiedDesignAnalyzer()

    # Check if new methods exist
    methods = [
        "_analyze_active_site_rmsd",
        "_analyze_substrate_access",
        "_analyze_enzyme_preservation",
    ]

    print("\nNew analysis methods:")
    for method in methods:
        if hasattr(analyzer, method):
            print(f"  [PASS] {method}")
        else:
            print(f"  [FAIL] {method} not found")

    # Test with mock PDB content
    mock_pdb = """ATOM      1  CA  ALA A  45      10.000  10.000  10.000  1.00  0.00           C
ATOM      2  CA  GLU A  78      15.000  15.000  15.000  1.00  0.00           C
ATOM      3  CA  HIS A 123      20.000  20.000  20.000  1.00  0.00           C
HETATM    4  O1  PQQ L   1      12.000  12.000  12.000  1.00  0.00           O
HETATM    5  O2  PQQ L   1      13.000  13.000  13.000  1.00  0.00           O
END
"""

    # Test substrate access analysis
    print("\nTesting _analyze_substrate_access:")
    result = analyzer._analyze_substrate_access(
        mock_pdb,
        {"L": "O1,O2,O3"}  # O3 doesn't exist
    )
    print(f"  Status: {result.status}")
    if result.metrics:
        print(f"  Found: {result.metrics.get('substrate_channel', {}).get('atoms_found', 0)}")
        print(f"  Expected: {result.metrics.get('substrate_channel', {}).get('atoms_expected', 0)}")

    # Test enzyme preservation analysis
    print("\nTesting _analyze_enzyme_preservation:")
    result = analyzer._analyze_enzyme_preservation(
        designed_pdb=mock_pdb,
        enzyme_class="quinoprotein",
        exposed_atoms={"L": "O1,O2"},
    )
    print(f"  Status: {result.status}")
    if result.metrics:
        print(f"  Enzyme class: {result.metrics.get('enzyme_class')}")
        print(f"  Checks performed: {result.metrics.get('checks_performed', [])}")
        print(f"  Passed: {result.metrics.get('passed', False)}")

    return True


def main():
    """Run all tests."""
    print("\n" + "=" * 70)
    print("AI INFRASTRUCTURE ENHANCEMENT TEST")
    print("Testing with original query that prompted this enhancement")
    print("=" * 70)

    try:
        # Test 1: NL Parser
        intent = test_nl_parser()

        # Test 2: Stability Profile
        stability = test_stability_profile(intent)

        # Test 3: Enzyme Detection
        test_enzyme_detection()

        # Test 4: RFD3 Config
        test_rfd3_config()

        # Test 5: Unified Analyzer
        test_unified_analyzer()

        print("\n" + "=" * 70)
        print("ALL TESTS COMPLETED")
        print("=" * 70)

        # Summary
        print("\nSummary for original query:")
        print(f"  - PDB: {intent.source_pdb_id}")
        print(f"  - Ligand: {intent.ligand_name}")
        print(f"  - Metal: {intent.metal_type}")
        print(f"  - Stability: {stability.name} (cfg={stability.cfg_scale}, step={stability.step_scale})")
        print(f"  - Enzyme: {intent.enzyme_class}")
        print(f"  - Preserve Function: {intent.preserve_function}")

    except Exception as e:
        print(f"\n[FAIL] TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
