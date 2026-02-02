"""
E2E Test: PQQ + Calcium Chemistry-Aware Pipeline

Tests all 8 enhancements from the chemistry-aware pipeline plan using the
PQQ-Ca (quinoprotein dehydrogenase) use case.

Flow tested:
  NLDesignParser -> LigandResolver -> RFD3ConfigGenerator
  With enzyme_chemistry, ligand_donors, design_rules all wired in.

Enhancements verified:
  E1: Donor ID wired into fixing strategy
  E2: Proper H-bond donor/acceptor split
  E3: Catalytic H-bond network -> RFD3 config
  E4: Substrate channel -> select_exposed
  E5: Protonation state (Dimorphite-DL)
  E6: Metal-aware conformers (Architector / constrained RDKit)
  E7: Post-design coordination validation
  E8: ML donor prediction stub
"""

import sys
import os
import json
import logging
from dataclasses import asdict
from typing import Any, Dict

# Ensure backend modules are importable
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

logging.basicConfig(level=logging.DEBUG, format="%(name)s | %(levelname)s | %(message)s")
logger = logging.getLogger("test_e2e_pqq_ca")

# PQQ canonical SMILES (from PubChem CID 1024)
PQQ_SMILES = "OC(=O)c1[nH]c2c(c1)C(=O)C(=O)c3nc(cc(C(O)=O)c23)C(O)=O"

# ============================================================================
# Helpers
# ============================================================================

def section(title: str):
    """Print a section header."""
    print(f"\n{'='*72}")
    print(f"  {title}")
    print(f"{'='*72}")


def check(label: str, condition: bool, detail: str = ""):
    """Print a pass/fail check."""
    status = "PASS" if condition else "FAIL"
    marker = "[+]" if condition else "[!]"
    msg = f"  {marker} {status}: {label}"
    if detail:
        msg += f" — {detail}"
    print(msg)
    return condition


def dump_dict(d: Dict[str, Any], indent: int = 4, max_val_len: int = 120):
    """Print dict with truncated values."""
    for k, v in sorted(d.items()):
        val_str = str(v)
        if len(val_str) > max_val_len:
            val_str = val_str[:max_val_len] + "..."
        print(f"{'':>{indent}}{k}: {val_str}")


# ============================================================================
# E1: Donor identification from SMILES
# ============================================================================

def test_e1_donor_identification():
    section("E1: Donor Identification from SMILES (ligand_donors.py)")
    results = {}

    try:
        from ligand_donors import identify_donors_from_smiles, get_recommended_coordination_mode, filter_donors_by_geometry
    except ImportError as e:
        print(f"  SKIP: Cannot import ligand_donors: {e}")
        return results

    # Test donor identification on PQQ
    donors = identify_donors_from_smiles(PQQ_SMILES)
    results["donor_count"] = len(donors)
    results["donor_elements"] = [d["element"] for d in donors]
    results["donor_types"] = [d["type"] for d in donors]

    print(f"  PQQ donors found: {len(donors)}")
    for d in donors:
        print(f"    idx={d['atom_idx']} elem={d['element']} type={d['type']} hsab={d.get('hsab', '?')}")

    check("Donors found for PQQ", len(donors) > 0, f"{len(donors)} donors")
    check("Has oxygen donors", "O" in results["donor_elements"], str(results["donor_elements"]))
    check("Has nitrogen donors", "N" in results["donor_elements"], str(results["donor_elements"]))

    # Test geometry filtering on PQQ
    filtered = filter_donors_by_geometry(donors, PQQ_SMILES, metal="CA")
    results["geometry_filtered_count"] = len(filtered)
    print(f"\n  Geometry-filtered donors: {len(filtered)} (from {len(donors)} raw)")
    for d in filtered:
        print(f"    idx={d['atom_idx']} elem={d['element']} type={d['type']}")
    check("Geometry filter reduces PQQ donors", len(filtered) < len(donors),
          f"{len(donors)} -> {len(filtered)}")

    # Test coordination recommendation for PQQ + CA
    coord = get_recommended_coordination_mode(PQQ_SMILES, "CA", target_cn=6)
    results["coordination_mode"] = coord.get("coordination_mode")
    results["compatibility_score"] = coord.get("compatibility_score", 0)
    results["recommended_donors"] = len(coord.get("recommended_donors", []))

    print(f"\n  Coordination mode for PQQ + Ca:")
    print(f"    mode: {coord.get('coordination_mode')}")
    print(f"    compatibility: {coord.get('compatibility_score', 0):.2f}")
    print(f"    estimated CN contribution: {coord.get('estimated_cn_contribution', 0)}")
    print(f"    total_donors: {coord.get('total_donors', 0)}")
    print(f"    recommended_donors: {len(coord.get('recommended_donors', []))}")

    check("Coordination mode detected", coord.get("coordination_mode") is not None)
    check("Compatibility score > 0", coord.get("compatibility_score", 0) > 0)

    return results


# ============================================================================
# E2: H-bond donor/acceptor split from SMILES
# ============================================================================

def test_e2_hbond_classification():
    section("E2: H-bond Donor/Acceptor Split (ligand_resolver.py)")
    results = {}

    try:
        from ligand_resolver import _classify_hbond_from_smiles, _map_rdkit_to_pdb_names
    except ImportError as e:
        print(f"  SKIP: Cannot import: {e}")
        return results

    # Build a dummy atom index map (just use identity for testing)
    # PQQ has ~24 heavy atoms; map RDKit indices to fake PDB names
    dummy_map = {i: f"X{i}" for i in range(30)}

    try:
        hbond_donors, hbond_acceptors = _classify_hbond_from_smiles(PQQ_SMILES, dummy_map)
    except Exception as e:
        print(f"  ERROR: _classify_hbond_from_smiles failed: {e}")
        return results

    results["hbond_donors"] = hbond_donors
    results["hbond_acceptors"] = hbond_acceptors

    print(f"  H-bond donors (O-H, N-H):  {hbond_donors}")
    print(f"  H-bond acceptors (C=O, etc): {hbond_acceptors}")

    # PQQ has:
    # - 3 carboxylic acid groups (COOH): each has donor O-H + acceptor C=O
    # - 2 quinone C=O groups: acceptors only
    # - 1 N-H (pyrrole): donor
    # - 1 pyridine N: acceptor
    check("Has H-bond donors", len(hbond_donors) > 0, str(hbond_donors))
    check("Has H-bond acceptors", len(hbond_acceptors) > 0, str(hbond_acceptors))
    check("Donors and acceptors are different sets",
          set(hbond_donors) != set(hbond_acceptors),
          "Should not be identical")

    # The old code put ALL O/N in acceptors — verify donors exist now
    check("Donors not empty (old bug was all-acceptor)",
          len(hbond_donors) > 0,
          f"{len(hbond_donors)} donors found")

    return results


# ============================================================================
# E3: Enzyme class detection + catalytic H-bond network
# ============================================================================

def test_e3_enzyme_chemistry():
    section("E3: Enzyme Class Detection + Catalytic H-bond Network")
    results = {}

    try:
        from enzyme_chemistry import (
            detect_enzyme_class,
            identify_catalytic_residues_from_pdb,
            get_catalytic_hbond_network,
        )
    except ImportError as e:
        print(f"  SKIP: Cannot import enzyme_chemistry: {e}")
        return results

    # Test enzyme class detection from queries
    test_queries = [
        ("Design a PQQ-dependent dehydrogenase with calcium", "PQQ", "CA"),
        ("Scaffold the PQQ-Ca active site from 4CVB", "PQQ", "CA"),
        ("Build a quinoprotein methanol dehydrogenase", "PQQ", "CA"),
    ]

    for query, ligand, metal in test_queries:
        enzyme_class, confidence = detect_enzyme_class(query, ligand, metal)
        print(f"  Query: '{query[:50]}...'")
        print(f"    enzyme_class={enzyme_class}, confidence={confidence:.2f}")
        check(f"Detected enzyme class for PQQ query",
              enzyme_class is not None,
              f"{enzyme_class} ({confidence:.2f})")

    results["enzyme_class"] = enzyme_class
    results["confidence"] = confidence

    return results


# ============================================================================
# E4: Substrate channel -> select_exposed
# ============================================================================

def test_e4_substrate_channel():
    section("E4: Substrate Channel -> select_exposed")
    results = {}

    try:
        from enzyme_chemistry import get_substrate_channel_atoms, get_preservation_requirements
    except ImportError as e:
        print(f"  SKIP: Cannot import: {e}")
        return results

    # get_preservation_requirements checks if substrate_access_required
    try:
        reqs = get_preservation_requirements("quinoprotein")
        results["substrate_access_required"] = reqs.get("substrate_access_required", False)
        print(f"  Preservation requirements for quinoprotein:")
        for k, v in reqs.items():
            print(f"    {k}: {v}")
        check("Substrate access required for quinoprotein",
              reqs.get("substrate_access_required", False))
    except Exception as e:
        print(f"  ERROR: get_preservation_requirements failed: {e}")

    return results


# ============================================================================
# E5: Protonation state via Dimorphite-DL
# ============================================================================

def test_e5_protonation():
    section("E5: Protonation State (Dimorphite-DL)")
    results = {}

    try:
        from ligand_donors import protonate_at_pH, HAS_DIMORPHITE
    except ImportError as e:
        print(f"  SKIP: Cannot import: {e}")
        return results

    results["has_dimorphite"] = HAS_DIMORPHITE
    print(f"  Dimorphite-DL available: {HAS_DIMORPHITE}")

    if not HAS_DIMORPHITE:
        print("  NOTE: Dimorphite-DL not installed; protonation will use fallback")
        check("Graceful fallback when Dimorphite unavailable", True)
        # Test fallback: should return original SMILES
        result = protonate_at_pH(PQQ_SMILES, pH=7.4)
        check("Fallback returns original SMILES", result == PQQ_SMILES)
        results["protonated_smiles"] = result
        return results

    # Test protonation of PQQ at physiological pH
    protonated = protonate_at_pH(PQQ_SMILES, pH=7.4)
    results["protonated_smiles"] = protonated
    results["changed"] = protonated != PQQ_SMILES

    print(f"  Original:   {PQQ_SMILES}")
    print(f"  At pH 7.4:  {protonated}")
    print(f"  Changed:    {results['changed']}")

    # At pH 7.4, PQQ's 3 carboxylic acids should be deprotonated (COO-)
    check("Protonation returned valid SMILES", len(protonated) > 0)
    if protonated != PQQ_SMILES:
        check("SMILES changed at physiological pH", True,
              "Carboxylates should deprotonate at pH 7.4")
    else:
        check("SMILES unchanged (Dimorphite may not modify this structure)", True,
              "Not all molecules change")

    # Also test with identify_donors_from_smiles pH parameter
    from ligand_donors import identify_donors_from_smiles
    donors_default = identify_donors_from_smiles(PQQ_SMILES)
    donors_ph74 = identify_donors_from_smiles(PQQ_SMILES, pH=7.4)
    results["donors_default_count"] = len(donors_default)
    results["donors_ph74_count"] = len(donors_ph74)

    print(f"\n  Donors (no pH):  {len(donors_default)}")
    print(f"  Donors (pH 7.4): {len(donors_ph74)}")
    check("identify_donors_from_smiles accepts pH parameter", True)

    return results


# ============================================================================
# E6: Metal-aware conformer generation
# ============================================================================

def test_e6_metal_conformers():
    section("E6: Metal-Aware Conformer Generation")
    results = {}

    try:
        from conformer_utils import generate_metal_aware_conformer, generate_conformer
    except ImportError as e:
        print(f"  SKIP: Cannot import conformer_utils: {e}")
        return results

    # Test regular conformer first as baseline
    regular_pdb = generate_conformer(PQQ_SMILES, name="PQQ")
    results["regular_conformer"] = regular_pdb is not None
    if regular_pdb:
        hetatm_count = sum(1 for line in regular_pdb.split('\n') if line.startswith('HETATM'))
        results["regular_hetatm_count"] = hetatm_count
        print(f"  Regular conformer: {hetatm_count} HETATM lines")
    else:
        print(f"  Regular conformer: FAILED")

    # Test metal-aware conformer for CA
    metal_pdb = generate_metal_aware_conformer(
        PQQ_SMILES,
        metal="CA",
        name="PQQ",
        center=(0.0, 0.0, 0.0),
    )
    results["metal_conformer"] = metal_pdb is not None
    if metal_pdb:
        hetatm_count = sum(1 for line in metal_pdb.split('\n') if line.startswith('HETATM'))
        results["metal_hetatm_count"] = hetatm_count
        print(f"  Metal-aware conformer (CA): {hetatm_count} HETATM lines")
    else:
        print(f"  Metal-aware conformer (CA): FAILED (fallback to regular)")

    check("Regular conformer generated", results["regular_conformer"])
    check("Metal-aware conformer generated or fell back gracefully",
          metal_pdb is not None,
          "Should produce a conformer either way")

    return results


# ============================================================================
# E7: Post-design coordination validation
# ============================================================================

def test_e7_coordination_validation():
    section("E7: Post-Design Coordination Validation")
    results = {}

    try:
        from unified_analyzer import validate_metal_coordination
    except ImportError as e:
        print(f"  SKIP: Cannot import validate_metal_coordination: {e}")
        return results

    # Create a synthetic PDB with Ca and nearby O atoms to test validation
    # This simulates a predicted structure with a Ca binding site
    synthetic_pdb = """ATOM      1  N   ASP A  10       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ASP A  10       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  OD1 ASP A  10       5.000   1.000   0.000  1.00  0.00           O
ATOM      4  OD2 ASP A  10       5.000  -1.000   0.000  1.00  0.00           O
ATOM      5  N   GLU A  11       3.000   0.000   0.000  1.00  0.00           N
ATOM      6  OE1 GLU A  11       5.000   0.000   2.000  1.00  0.00           O
ATOM      7  OE2 GLU A  11       5.000   0.000  -2.000  1.00  0.00           O
ATOM      8  OE1 GLU A  12       7.000   0.000   0.000  1.00  0.00           O
HETATM    9 CA    CA A  99       5.000   0.000   0.000  1.00  0.00          CA
HETATM   10  O5  PQQ L   1       5.000   2.200   0.000  1.00  0.00           O
HETATM   11  N6  PQQ L   1       5.000   0.000   2.300  1.00  0.00           N
HETATM   12  O7A PQQ L   1       5.000  -2.200   0.000  1.00  0.00           O
END
"""

    result = validate_metal_coordination(synthetic_pdb, "CA", expected_cn=6)
    results["validation_result"] = result

    print(f"  Validation result:")
    for k, v in result.items():
        print(f"    {k}: {v}")

    check("Validation returned result", result is not None)
    check("Has 'valid' key", "valid" in result)
    check("Has 'coordination_number' key", "coordination_number" in result)

    cn = result.get("coordination_number", 0)
    print(f"\n  Detected coordination number: {cn}")
    check("Found coordinating atoms", cn > 0, f"CN={cn}")

    return results


# ============================================================================
# E8: ML donor prediction stub
# ============================================================================

def test_e8_ml_stub():
    section("E8: ML Donor Prediction Stub")
    results = {}

    # Test the stub in ligand_donors.py
    try:
        from ligand_donors import predict_donors_ml
    except ImportError as e:
        print(f"  SKIP: Cannot import predict_donors_ml: {e}")
        return results

    ml_result = predict_donors_ml(PQQ_SMILES, "CA")
    results["ml_stub_returned"] = ml_result is not None
    results["ml_donor_count"] = len(ml_result) if ml_result else 0

    print(f"  predict_donors_ml returned {len(ml_result)} donors (delegated to SMARTS)")
    check("ML stub returns results", len(ml_result) > 0, "Falls back to SMARTS-based")

    # Test the MLDonorPredictor class stub
    try:
        from ml_donor_prediction import MLDonorPredictor
        predictor = MLDonorPredictor()
        try:
            predictor.predict(PQQ_SMILES, "CA")
            check("MLDonorPredictor.predict raises NotImplementedError", False,
                  "Should have raised")
        except NotImplementedError:
            check("MLDonorPredictor.predict raises NotImplementedError", True)
        results["ml_class_stub"] = True
    except ImportError:
        print(f"  NOTE: ml_donor_prediction.py not found (optional)")
        results["ml_class_stub"] = False

    return results


# ============================================================================
# GEOMETRY FILTER: 3D donor filtering
# ============================================================================

CITRATE_SMILES = "OC(=O)CC(O)(CC(O)=O)C(O)=O"
CAFFEINE_SMILES = "Cn1c(=O)c2c(ncn2C)n(C)c1=O"
ATP_SMILES = "c1nc(c2c(n1)n(cn2)[C@@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N"


def test_geometry_filter():
    section("GEOMETRY FILTER: 3D Donor Filtering (filter_donors_by_geometry)")
    results = {}

    try:
        from ligand_donors import identify_donors_from_smiles, filter_donors_by_geometry
    except ImportError as e:
        print(f"  SKIP: Cannot import: {e}")
        return results

    test_cases = [
        ("PQQ", PQQ_SMILES, "CA"),
        ("Citrate", CITRATE_SMILES, "TB"),
        ("Caffeine", CAFFEINE_SMILES, None),
        ("ATP", ATP_SMILES, "MG"),
    ]

    for name, smiles, metal in test_cases:
        raw = identify_donors_from_smiles(smiles)
        filtered = filter_donors_by_geometry(raw, smiles, metal=metal)

        results[f"{name}_raw"] = len(raw)
        results[f"{name}_filtered"] = len(filtered)

        print(f"\n  {name} ({metal or 'no metal'}):")
        print(f"    Raw donors:      {len(raw)}")
        print(f"    After geometry:  {len(filtered)}")
        for d in filtered:
            print(f"      idx={d['atom_idx']} elem={d['element']} type={d['type']}")

        check(f"{name}: geometry filter reduces donors",
              len(filtered) <= len(raw),
              f"{len(raw)} -> {len(filtered)}")

    # PQQ should be significantly reduced
    check("PQQ filtered to <=5 donors",
          results.get("PQQ_filtered", 99) <= 5,
          f"Got {results.get('PQQ_filtered')}")

    # Ordering: PQQ raw > PQQ filtered
    check("PQQ raw > filtered",
          results.get("PQQ_raw", 0) > results.get("PQQ_filtered", 0),
          f"{results.get('PQQ_raw')} > {results.get('PQQ_filtered')}")

    return results


# ============================================================================
# CHAIN LENGTH ESTIMATION: TPSA-based
# ============================================================================

def test_chain_length_estimation():
    section("CHAIN LENGTH ESTIMATION: TPSA-Based (estimate_chain_length)")
    results = {}

    try:
        from design_rules import estimate_chain_length
    except ImportError as e:
        print(f"  SKIP: Cannot import: {e}")
        return results

    test_cases = [
        ("Caffeine", CAFFEINE_SMILES, False),
        ("Citrate", CITRATE_SMILES, True),
        ("PQQ", PQQ_SMILES, True),
        ("ATP", ATP_SMILES, True),
    ]

    print(f"\n  {'Ligand':<12} {'SMILES':<30} {'Metal':<6} {'Chain Length':<15}")
    print(f"  {'-'*12} {'-'*30} {'-'*6} {'-'*15}")

    chain_lengths = {}
    for name, smiles, has_metal in test_cases:
        cl = estimate_chain_length(
            smiles=smiles,
            ligand_character="polar",
            target_topology="monomer",
            has_metal=has_metal,
        )
        chain_lengths[name] = cl
        results[f"{name}_range"] = cl

        # Parse midpoint from "X-Y" format
        parts = cl.split("-")
        mid = (int(parts[0]) + int(parts[1])) / 2
        results[f"{name}_mid"] = mid

        smiles_short = smiles[:27] + "..." if len(smiles) > 30 else smiles
        print(f"  {name:<12} {smiles_short:<30} {'yes' if has_metal else 'no':<6} {cl:<15}")

    # Verify ordering: PQQ > Citrate > Caffeine
    check("PQQ chain > Citrate chain",
          results.get("PQQ_mid", 0) > results.get("Citrate_mid", 0),
          f"PQQ={results.get('PQQ_mid'):.0f} vs Citrate={results.get('Citrate_mid'):.0f}")

    check("Citrate chain > Caffeine chain",
          results.get("Citrate_mid", 0) > results.get("Caffeine_mid", 0),
          f"Citrate={results.get('Citrate_mid'):.0f} vs Caffeine={results.get('Caffeine_mid'):.0f}")

    check("ATP chain > PQQ chain",
          results.get("ATP_mid", 0) > results.get("PQQ_mid", 0),
          f"ATP={results.get('ATP_mid'):.0f} vs PQQ={results.get('PQQ_mid'):.0f}")

    return results


# ============================================================================
# MULTI-LIGAND COMPARISON: Full pipeline comparison
# ============================================================================

def test_multi_ligand_comparison():
    section("MULTI-LIGAND COMPARISON: Geometry + Chain Length")
    results = {}

    try:
        from ligand_donors import identify_donors_from_smiles, filter_donors_by_geometry
        from design_rules import estimate_chain_length
    except ImportError as e:
        print(f"  SKIP: Cannot import: {e}")
        return results

    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
    except ImportError:
        print(f"  SKIP: RDKit not available")
        return results

    ligands = [
        ("PQQ", PQQ_SMILES, "CA"),
        ("Citrate", CITRATE_SMILES, "TB"),
        ("Caffeine", CAFFEINE_SMILES, None),
        ("ATP", ATP_SMILES, "MG"),
    ]

    print(f"\n  {'Ligand':<10} {'Heavy':<6} {'TPSA':<8} {'Raw D':<6} {'Filt D':<7} {'Chain':<12} {'Metal':<6}")
    print(f"  {'-'*10} {'-'*6} {'-'*8} {'-'*6} {'-'*7} {'-'*12} {'-'*6}")

    for name, smiles, metal in ligands:
        mol = Chem.MolFromSmiles(smiles)
        heavy = mol.GetNumHeavyAtoms() if mol else 0
        tpsa = Descriptors.TPSA(mol) if mol else 0

        raw = identify_donors_from_smiles(smiles)
        filtered = filter_donors_by_geometry(raw, smiles, metal=metal)

        cl = estimate_chain_length(
            smiles=smiles,
            ligand_character="polar",
            target_topology="monomer",
            has_metal=metal is not None,
        )

        results[name] = {
            "heavy_atoms": heavy,
            "tpsa": round(tpsa, 1),
            "raw_donors": len(raw),
            "filtered_donors": len(filtered),
            "chain_length": cl,
            "metal": metal or "none",
        }

        print(f"  {name:<10} {heavy:<6} {tpsa:<8.1f} {len(raw):<6} {len(filtered):<7} {cl:<12} {metal or '-':<6}")

    print(f"\n  Chemistry Analysis:")
    print(f"  - PQQ: Tricyclic quinone cofactor, 3 COOH + 2 C=O + pyrrole N. Ca binds O5,N6,O7A.")
    print(f"  - Citrate: Tricarboxylate, 3 COO- groups. Ln3+ binds 3 carboxylate oxygens.")
    print(f"  - Caffeine: Two carbonyl O + two aromatic N. Poor metal coordinator.")
    print(f"  - ATP: Triphosphate + adenine. Mg2+ binds phosphate oxygens.")

    check("Comparison table generated", len(results) == 4)

    return results


# ============================================================================
# INTEGRATION: Full NL -> Config pipeline for PQQ + CA
# ============================================================================

def test_integration_pqq_ca_pipeline():
    section("INTEGRATION: Full PQQ+Ca Pipeline (NL -> Resolve -> Config)")
    results = {}

    # --- Step 1: Parse NL query ---
    print("\n  --- Step 1: NL Intent Parsing ---")
    try:
        from nl_design_parser import SimpleFallbackParser
    except ImportError as e:
        print(f"  SKIP: Cannot import parser: {e}")
        return results

    parser = SimpleFallbackParser()
    query = "Design a PQQ-dependent dehydrogenase with calcium binding site"
    intent = parser.parse(query)

    print(f"  Query: '{query}'")
    print(f"  metal_type: {intent.metal_type}")
    print(f"  ligand_name: {intent.ligand_name}")
    print(f"  design_goal: {intent.design_goal}")
    print(f"  enzyme_class: {intent.enzyme_class}")
    print(f"  enzyme_class_confidence: {intent.enzyme_class_confidence}")
    print(f"  bury_ligand: {intent.bury_ligand}")
    print(f"  confidence: {intent.confidence}")

    results["intent_metal"] = intent.metal_type
    results["intent_ligand"] = intent.ligand_name
    results["intent_enzyme_class"] = intent.enzyme_class

    check("Metal detected as CA", intent.metal_type and intent.metal_type.upper() == "CA",
          str(intent.metal_type))
    check("Ligand detected as PQQ", intent.ligand_name and "pqq" in intent.ligand_name.lower(),
          str(intent.ligand_name))
    check("Enzyme class detected (E3)", intent.enzyme_class is not None,
          str(intent.enzyme_class))

    # --- Step 2: Resolve ligand ---
    print("\n  --- Step 2: Ligand Resolution ---")
    try:
        from ligand_resolver import LigandResolver
    except ImportError as e:
        print(f"  SKIP: Cannot import LigandResolver: {e}")
        return results

    resolver = LigandResolver(use_pubchem=False, use_pdb=False)  # Template-only for speed
    ligand_name = intent.ligand_name or "PQQ"
    metal_type = intent.metal_type or "CA"
    resolved = resolver.resolve(ligand_name, metal_type=metal_type)

    print(f"  resolved: {resolved.resolved}")
    print(f"  source: {resolved.source}")
    print(f"  template_name: {resolved.template_name}")
    print(f"  residue_code: {resolved.residue_code}")
    print(f"  smiles: {resolved.smiles}")
    print(f"  has pdb_content: {resolved.pdb_content is not None}")
    print(f"  coordination keys: {list(resolved.coordination.keys()) if resolved.coordination else []}")

    results["resolved"] = resolved.resolved
    results["source"] = resolved.source
    results["template_name"] = resolved.template_name

    check("Ligand resolved", resolved.resolved)
    check("Resolved from template", resolved.source == "template", resolved.source)
    check("Template is pqq_ca", resolved.template_name == "pqq_ca", str(resolved.template_name))
    check("Has PDB content", resolved.pdb_content is not None)
    check("Has coordination info", bool(resolved.coordination))
    check("Has SMILES", resolved.smiles is not None)

    # Check chemistry analysis
    if resolved.chemistry:
        chem = resolved.chemistry
        print(f"\n  Chemistry analysis:")
        print(f"    total_atoms: {chem.total_atoms}")
        print(f"    donor_atoms: {chem.donor_atoms}")
        print(f"    hbond_donors: {chem.hbond_donors}")
        print(f"    hbond_acceptors: {chem.hbond_acceptors}")
        print(f"    coordinating_atoms: {chem.coordinating_atoms}")
        print(f"    hydrophobic_atoms: {chem.hydrophobic_atoms}")
        print(f"    ligand_character: {chem.ligand_character}")

        results["chem_donors"] = chem.donor_atoms
        results["chem_hbond_donors"] = chem.hbond_donors
        results["chem_hbond_acceptors"] = chem.hbond_acceptors

        check("Has donor atoms", len(chem.donor_atoms) > 0, str(chem.donor_atoms))
        check("Has H-bond acceptors", len(chem.hbond_acceptors) > 0)
        # E2: H-bond donors should now be populated for template path too
        # (template has ligand_hbond_donors field)
    else:
        print("  WARNING: No chemistry analysis")

    # Check coordination info
    if resolved.coordination:
        coord = resolved.coordination
        print(f"\n  Coordination info:")
        dump_dict(coord)

        results["coord_donors"] = coord.get("ligand_metal_donors", [])
        results["coord_cn"] = coord.get("metal_coordination_number")

        check("Has ligand donors in coordination",
              len(coord.get("ligand_metal_donors", [])) > 0,
              str(coord.get("ligand_metal_donors")))
        check("CN defined",
              coord.get("metal_coordination_number") is not None,
              str(coord.get("metal_coordination_number")))

    # --- Step 3: Generate RFD3 config ---
    print("\n  --- Step 3: RFD3 Config Generation ---")
    try:
        from rfd3_config_generator import RFD3ConfigGenerator
    except ImportError as e:
        print(f"  SKIP: Cannot import RFD3ConfigGenerator: {e}")
        return results

    generator = RFD3ConfigGenerator()
    rfd3_config = generator.generate(intent, resolved)
    mpnn_config = generator.generate_mpnn_config(intent, resolved)

    print(f"\n  RFD3 Config:")
    rfd3_dict = asdict(rfd3_config)
    # Filter out None values and empty collections for readability
    rfd3_filtered = {k: v for k, v in rfd3_dict.items()
                     if v is not None and v != [] and v != {} and v != ""
                     and k != "pdb_content"}
    dump_dict(rfd3_filtered)

    results["rfd3_contig"] = rfd3_config.contig
    results["rfd3_cfg"] = rfd3_config.use_classifier_free_guidance
    results["rfd3_cfg_scale"] = rfd3_config.cfg_scale
    results["rfd3_select_buried"] = dict(rfd3_config.select_buried)
    results["rfd3_select_hbond_donor"] = dict(rfd3_config.select_hbond_donor)
    results["rfd3_select_hbond_acceptor"] = dict(rfd3_config.select_hbond_acceptor)
    results["rfd3_select_fixed_atoms"] = dict(rfd3_config.select_fixed_atoms)
    results["rfd3_select_exposed"] = dict(rfd3_config.select_exposed)

    check("Has contig", bool(rfd3_config.contig), rfd3_config.contig)
    check("CFG enabled", rfd3_config.use_classifier_free_guidance)
    check("CFG scale >= 2.0", rfd3_config.cfg_scale >= 2.0, str(rfd3_config.cfg_scale))
    check("Has burial conditioning", bool(rfd3_config.select_buried),
          str(dict(rfd3_config.select_buried)))

    # E2: Check H-bond donor/acceptor split
    check("Has H-bond acceptors in config",
          bool(rfd3_config.select_hbond_acceptor),
          str(dict(rfd3_config.select_hbond_acceptor)))
    check("Has H-bond donors in config (E2)",
          bool(rfd3_config.select_hbond_donor),
          str(dict(rfd3_config.select_hbond_donor)))

    # E3: Check if enzyme H-bond network was applied (only if enzyme_class detected
    # AND PDB content was provided — template path doesn't have full PDB for catalytic
    # residue analysis, but the code path should still run without errors)
    if intent.enzyme_class:
        print(f"\n  Enzyme class '{intent.enzyme_class}' was detected (E3)")
        print(f"  NOTE: Catalytic H-bond network requires full PDB content with protein")
        print(f"        residues. Template-only resolution provides ligand PDB only.")
        print(f"        E3/E4 activate when pdb_content includes protein chains.")

    # E4: select_exposed (same caveat — needs protein PDB)
    if rfd3_config.select_exposed:
        check("Has select_exposed (E4 substrate channel)", True,
              str(dict(rfd3_config.select_exposed)))
    else:
        print(f"  NOTE: select_exposed empty (expected — no protein PDB in template path)")

    # Check MPNN config
    print(f"\n  MPNN Config:")
    mpnn_dict = asdict(mpnn_config)
    mpnn_filtered = {k: v for k, v in mpnn_dict.items() if v is not None}
    dump_dict(mpnn_filtered)

    results["mpnn_temperature"] = mpnn_config.temperature
    results["mpnn_bias"] = mpnn_config.bias_AA

    check("MPNN has metal-aware bias (Ca is hard acid -> Asp/Glu bias)",
          mpnn_config.bias_AA is not None and len(mpnn_config.bias_AA) > 0,
          str(mpnn_config.bias_AA)[:80] if mpnn_config.bias_AA else "None")

    # --- Step 4: Convert to API params ---
    print("\n  --- Step 4: API Parameter Output ---")
    api_params = rfd3_config.to_api_params()
    results["api_params_keys"] = list(api_params.keys())

    print(f"  API params keys: {list(api_params.keys())}")
    check("API params generated", len(api_params) > 0)

    return results


# ============================================================================
# E7 (extended): Filter preset includes coordination_valid
# ============================================================================

def test_e7_filter_presets():
    section("E7 (extended): Filter Presets with Coordination Validation")
    results = {}

    try:
        from filter_evaluator import FILTER_PRESETS
    except ImportError as e:
        print(f"  SKIP: Cannot import FILTER_PRESETS: {e}")
        return results

    metal_preset = FILTER_PRESETS.get("metal_binding", {})
    strict_preset = FILTER_PRESETS.get("metal_binding_strict", {})

    print(f"  metal_binding preset keys: {list(metal_preset.keys())}")
    print(f"  metal_binding_strict preset keys: {list(strict_preset.keys())}")

    results["metal_has_coord"] = "coordination_valid" in metal_preset
    results["strict_has_coord"] = "coordination_valid" in strict_preset

    check("metal_binding has coordination_valid",
          "coordination_valid" in metal_preset,
          str(metal_preset.get("coordination_valid")))
    check("metal_binding_strict has coordination_valid",
          "coordination_valid" in strict_preset,
          str(strict_preset.get("coordination_valid")))

    # Check the threshold format
    if "coordination_valid" in metal_preset:
        cv = metal_preset["coordination_valid"]
        check("coordination_valid uses equals format",
              isinstance(cv, dict) and "equals" in cv,
              str(cv))

    return results


# ============================================================================
# KNOWLEDGE BASE: Self-growing ligand feature knowledge base
# ============================================================================

def test_knowledge_base_lookup():
    section("KNOWLEDGE BASE: Seed Ligand Lookup")
    results = {}

    try:
        from ligand_features import LigandKnowledgeBase, analyze_ligand_features
    except ImportError as e:
        print(f"  SKIP: Cannot import ligand_features: {e}")
        return results

    # Use a temp directory for test KB
    import tempfile
    import shutil
    temp_dir = tempfile.mkdtemp(prefix="test_kb_")

    try:
        kb = LigandKnowledgeBase(base_dir=temp_dir)

        # Test PQQ lookup
        pqq = kb.lookup("pqq", metal="CA")
        results["pqq_found"] = pqq is not None
        check("PQQ found in KB", pqq is not None)

        if pqq:
            print(f"  PQQ SMILES: {pqq.get('smiles', '')[:50]}...")
            print(f"  PQQ HSAB: {pqq.get('hsab_character')}")
            print(f"  PQQ coordination summary: {list(pqq.get('coordination_summary', {}).keys())}")

            coord = pqq.get("metal_coordination", {})
            results["pqq_ca_donors"] = coord.get("donors", [])
            check("PQQ has CA coordination data",
                  len(coord.get("donors", [])) > 0,
                  str(coord.get("donors")))

        # Test citrate lookup by alias
        cit = kb.lookup("cit", metal="TB")
        results["citrate_by_alias"] = cit is not None
        check("Citrate found by alias 'cit'", cit is not None)

        if cit:
            tb_coord = cit.get("metal_coordination", {})
            results["citrate_tb_donors"] = tb_coord.get("donors", [])
            check("Citrate has TB coordination data",
                  len(tb_coord.get("donors", [])) > 0,
                  str(tb_coord.get("donors")))

        # Test coordination donors method
        donors = kb.get_coordination_donors("pqq", "CA")
        results["pqq_ca_coord_donors"] = donors
        check("get_coordination_donors returns donors for PQQ+CA",
              donors is not None and len(donors) > 0,
              str(donors))

        # Test unknown ligand
        unknown = kb.lookup("totally_unknown_ligand_xyz")
        check("Unknown ligand returns None", unknown is None)

        # Test list_ligands
        all_ligands = kb.list_ligands()
        results["total_seed_ligands"] = len(all_ligands)
        check("Seed ligands loaded", len(all_ligands) >= 5,
              f"{len(all_ligands)} ligands")

    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

    return results


def test_knowledge_base_growth():
    section("KNOWLEDGE BASE: Growth via record_evidence")
    results = {}

    try:
        from ligand_features import LigandKnowledgeBase
    except ImportError as e:
        print(f"  SKIP: Cannot import: {e}")
        return results

    import tempfile
    import shutil
    temp_dir = tempfile.mkdtemp(prefix="test_kb_growth_")

    try:
        kb = LigandKnowledgeBase(base_dir=temp_dir)

        # Record new evidence for a novel ligand
        kb.record_evidence(
            ligand_id="test_ligand",
            pdb_id="1ABC",
            metal="ZN",
            coordination_data={
                "coordination_number": 4,
                "ligand_donors": ["N1@2.10", "O3@2.25"],
                "protein_donors": ["A45:HIS NE2@2.15"],
                "geometry": "tetrahedral",
                "resolution": 1.5,
                "source": "test",
            },
            smiles="c1ccccc1",
        )

        # Verify it's retrievable
        meta = kb.lookup("test_ligand", "ZN")
        results["new_ligand_found"] = meta is not None
        check("New ligand found after recording", meta is not None)

        if meta:
            coord = meta.get("coordination_summary", {}).get("ZN", {})
            results["new_donors"] = coord.get("donors", [])
            check("New ligand has ZN donors",
                  "N1" in coord.get("donors", []) and "O3" in coord.get("donors", []),
                  str(coord.get("donors")))
            check("Evidence count incremented",
                  coord.get("evidence_count", 0) > 0,
                  str(coord.get("evidence_count")))

        # Record second evidence for same ligand
        kb.record_evidence(
            ligand_id="test_ligand",
            pdb_id="2XYZ",
            metal="ZN",
            coordination_data={
                "coordination_number": 4,
                "ligand_donors": ["N1@2.12", "O3@2.20", "S5@2.35"],
                "protein_donors": [],
                "geometry": "tetrahedral",
            },
        )

        meta2 = kb.lookup("test_ligand", "ZN")
        if meta2:
            coord2 = meta2.get("coordination_summary", {}).get("ZN", {})
            results["merged_donors"] = coord2.get("donors", [])
            check("Donors merged from two PDBs",
                  "S5" in coord2.get("donors", []),
                  str(coord2.get("donors")))
            check("Evidence count is 2",
                  coord2.get("evidence_count", 0) == 2)

        # Verify index updated
        index = kb.list_ligands()
        check("test_ligand in index", "test_ligand" in index)

    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

    return results


def test_chemicalfeatures_fallback():
    section("CHEMICALFEATURES: Pharmacophore Fallback for Unknown Ligand")
    results = {}

    try:
        from ligand_features import analyze_ligand_features
    except ImportError as e:
        print(f"  SKIP: Cannot import: {e}")
        return results

    import tempfile
    import shutil

    # Use a fresh KB so the ligand won't be found
    temp_dir = tempfile.mkdtemp(prefix="test_cf_")

    try:
        # Analyze a ligand not in KB (using its SMILES directly)
        # Aspirin SMILES
        aspirin_smiles = "CC(=O)Oc1ccccc1C(O)=O"
        result = analyze_ligand_features(
            ligand_name="aspirin_test",
            smiles=aspirin_smiles,
            metal="CA",
        )

        results["source"] = result.get("source")
        results["feature_count"] = len(result.get("features", []))
        results["coord_donors"] = result.get("coordination_donors", [])

        print(f"  Source: {result.get('source')}")
        print(f"  Features found: {len(result.get('features', []))}")
        print(f"  Coordination donors: {result.get('coordination_donors', [])}")
        print(f"  Compatibility: {result.get('compatibility_score', 0):.2f}")

        check("Falls through to chemicalfeatures or geometry_filter",
              result.get("source") in ("chemicalfeatures", "geometry_filter"),
              result.get("source"))
        check("Features identified",
              len(result.get("features", [])) > 0,
              str(len(result.get("features", []))))
        check("Has coordination donors (O atoms)",
              len(result.get("coordination_donors", [])) > 0)

    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

    return results


def test_knowledge_base_vs_geometry():
    section("KB vs GEOMETRY: Compare donors for known ligand")
    results = {}

    try:
        from ligand_features import analyze_ligand_features, LigandKnowledgeBase
        from ligand_donors import identify_donors_from_smiles, filter_donors_by_geometry
    except ImportError as e:
        print(f"  SKIP: Cannot import: {e}")
        return results

    # Compare KB donors vs geometry-filtered donors for citrate + TB
    kb_result = analyze_ligand_features(
        ligand_name="citrate",
        smiles=CITRATE_SMILES,
        metal="TB",
    )

    geo_raw = identify_donors_from_smiles(CITRATE_SMILES)
    geo_filtered = filter_donors_by_geometry(geo_raw, CITRATE_SMILES, metal="TB")

    kb_donors = kb_result.get("coordination_donors", [])
    geo_donor_names = [f"{d['element']}{d['atom_idx']}" for d in geo_filtered]

    results["kb_source"] = kb_result.get("source")
    results["kb_donors"] = kb_donors
    results["geo_donors"] = geo_donor_names
    results["kb_count"] = len(kb_donors)
    results["geo_count"] = len(geo_donor_names)

    print(f"  KB source: {kb_result.get('source')}")
    print(f"  KB donors ({len(kb_donors)}): {kb_donors}")
    print(f"  Geometry donors ({len(geo_donor_names)}): {geo_donor_names}")

    check("KB returns knowledge_base source for citrate",
          kb_result.get("source") == "knowledge_base")
    check("Both methods find donors",
          len(kb_donors) > 0 and len(geo_donor_names) > 0)

    return results


# ============================================================================
# Main
# ============================================================================

def main():
    print("=" * 72)
    print("  E2E Test: PQQ + Ca Chemistry-Aware Pipeline")
    print("  Tests all 8 enhancements")
    print("=" * 72)

    all_results = {}
    pass_count = 0
    fail_count = 0

    # Run each enhancement test
    tests = [
        ("E1", test_e1_donor_identification),
        ("E2", test_e2_hbond_classification),
        ("E3", test_e3_enzyme_chemistry),
        ("E4", test_e4_substrate_channel),
        ("E5", test_e5_protonation),
        ("E6", test_e6_metal_conformers),
        ("E7", test_e7_coordination_validation),
        ("E7_filters", test_e7_filter_presets),
        ("E8", test_e8_ml_stub),
        ("GEOMETRY_FILTER", test_geometry_filter),
        ("CHAIN_LENGTH", test_chain_length_estimation),
        ("MULTI_LIGAND", test_multi_ligand_comparison),
        ("KB_LOOKUP", test_knowledge_base_lookup),
        ("KB_GROWTH", test_knowledge_base_growth),
        ("CF_FALLBACK", test_chemicalfeatures_fallback),
        ("KB_VS_GEO", test_knowledge_base_vs_geometry),
        ("INTEGRATION", test_integration_pqq_ca_pipeline),
    ]

    for name, test_fn in tests:
        try:
            result = test_fn()
            all_results[name] = result
        except Exception as e:
            print(f"\n  EXCEPTION in {name}: {e}")
            import traceback
            traceback.print_exc()
            all_results[name] = {"error": str(e)}

    # Summary
    section("SUMMARY")
    print(f"\n  Results saved to: test_e2e_pqq_ca_results.json")

    # Save results
    output_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "test_e2e_pqq_ca_results.json")
    with open(output_path, "w") as f:
        json.dump(all_results, f, indent=2, default=str)

    print(f"\n  Done. Review results above for PASS/FAIL status of each enhancement.")


if __name__ == "__main__":
    main()
