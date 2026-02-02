"""
Focused chemistry validation for ligand_features.py

Tests real-world ligands to verify:
1. Seed knowledge base has chemically correct donors
2. ChemicalFeatures pharmacophore perception identifies sensible features
3. HSAB compatibility scoring is reasonable
4. analyze_ligand_features unified entry point works end-to-end
"""

import sys
import os
import json
import tempfile
import shutil

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def separator(title):
    print(f"\n{'='*60}")
    print(f"  {title}")
    print(f"{'='*60}")


def check(label, condition, detail=None):
    status = "PASS" if condition else "FAIL"
    print(f"  [{status}] {label}")
    if detail:
        print(f"         {detail}")
    return condition


# ============================================================================
# Test 1: Seed Knowledge Base Chemistry Validation
# ============================================================================
def test_seed_kb_chemistry():
    separator("Seed Knowledge Base - Chemistry Validation")

    from ligand_features import LigandKnowledgeBase, SEED_LIGANDS

    temp_dir = tempfile.mkdtemp(prefix="test_chem_")
    try:
        kb = LigandKnowledgeBase(base_dir=temp_dir)

        print(f"\n  Seed ligands: {list(SEED_LIGANDS.keys())}")
        all_pass = True

        # --- PQQ (Pyrroloquinoline quinone) ---
        print(f"\n  --- PQQ ---")
        pqq = kb.lookup("pqq", metal="CA")
        if not check("PQQ found", pqq is not None):
            return False

        pqq_donors = pqq.get("coordination_summary", {}).get("CA", {}).get("donors", [])
        print(f"  CA donors: {pqq_donors}")
        print(f"  HSAB character: {pqq.get('hsab_character')}")
        # PQQ coordinates Ca2+ via O5 (quinone), N6 (pyridine N), O7 (carboxylate)
        # The ONO tridentate pocket is well-characterized in PQQ-Ca crystal structures
        # N1 in our atom naming = N6 in PDB convention (pyridine ring nitrogen)
        donor_elements = set(d[0] for d in pqq_donors)
        all_pass &= check("PQQ-CA has O and N donors (ONO tridentate pocket)",
                          "O" in donor_elements and "N" in donor_elements,
                          f"donors={pqq_donors}, elements={donor_elements}")
        all_pass &= check("PQQ-CA denticity 2-3 (typical for PQQ)",
                          2 <= len(pqq_donors) <= 4,
                          f"denticity={len(pqq_donors)}")
        all_pass &= check("PQQ is hard ligand (matches Ca2+ hard acid)",
                          pqq.get("hsab_character") == "hard")

        # --- Citrate ---
        print(f"\n  --- CITRATE ---")
        cit = kb.lookup("citrate", metal="TB")
        if not check("Citrate found", cit is not None):
            return False

        cit_tb = cit.get("coordination_summary", {}).get("TB", {}).get("donors", [])
        cit_ca = cit.get("coordination_summary", {}).get("CA", {}).get("donors", [])
        print(f"  TB donors: {cit_tb}")
        print(f"  CA donors: {cit_ca}")
        print(f"  HSAB: {cit.get('hsab_character')}")
        # Citrate: tricarboxylate, coordinates via carboxylate O atoms + hydroxyl O
        all_pass &= check("Citrate-TB has O donors",
                          all(d.startswith("O") for d in cit_tb),
                          f"donors={cit_tb}")
        all_pass &= check("Citrate can be tridentate with Tb",
                          len(cit_tb) >= 3,
                          f"denticity={len(cit_tb)}")
        all_pass &= check("Citrate is hard ligand",
                          cit.get("hsab_character") == "hard")

        # --- ATP ---
        print(f"\n  --- ATP ---")
        atp = kb.lookup("atp", metal="MG")
        if not check("ATP found", atp is not None):
            return False

        atp_mg = atp.get("coordination_summary", {}).get("MG", {}).get("donors", [])
        print(f"  MG donors: {atp_mg}")
        print(f"  HSAB: {atp.get('hsab_character')}")
        # ATP coordinates Mg2+ via phosphate oxygens
        # Key donors: O from beta and gamma phosphates
        all_pass &= check("ATP-MG has donors",
                          len(atp_mg) > 0,
                          f"donors={atp_mg}")
        all_pass &= check("ATP-MG donors are O atoms (phosphate)",
                          all(d.startswith("O") for d in atp_mg),
                          f"donors={atp_mg}")

        # --- EDTA ---
        print(f"\n  --- EDTA ---")
        edta = kb.lookup("edta", metal="CA")
        if not check("EDTA found", edta is not None):
            return False

        edta_ca = edta.get("coordination_summary", {}).get("CA", {}).get("donors", [])
        print(f"  CA donors: {edta_ca}")
        print(f"  HSAB: {edta.get('hsab_character')}")
        # EDTA is a hexadentate chelator: 4 carboxylate O + 2 amine N
        all_pass &= check("EDTA-CA has 4-6 donors (hexadentate chelator)",
                          4 <= len(edta_ca) <= 6,
                          f"denticity={len(edta_ca)}")
        edta_elements = set(d[0] for d in edta_ca)
        all_pass &= check("EDTA-CA has both N and O donors",
                          "N" in edta_elements and "O" in edta_elements,
                          f"elements={edta_elements}, donors={edta_ca}")

        # --- Heme ---
        print(f"\n  --- HEME ---")
        heme = kb.lookup("heme", metal="FE")
        if not check("Heme found", heme is not None):
            return False

        heme_fe = heme.get("coordination_summary", {}).get("FE", {}).get("donors", [])
        print(f"  FE donors: {heme_fe}")
        print(f"  HSAB: {heme.get('hsab_character')}")
        # Heme: Fe is coordinated by 4 pyrrole nitrogens
        all_pass &= check("Heme-FE has 4 N donors (porphyrin)",
                          len(heme_fe) == 4,
                          f"denticity={len(heme_fe)}")
        all_pass &= check("Heme-FE donors are all N atoms",
                          all(d.startswith("N") for d in heme_fe),
                          f"donors={heme_fe}")

        # --- Alias lookups ---
        print(f"\n  --- ALIAS LOOKUPS ---")
        cit_alias = kb.lookup("cit")
        all_pass &= check("'cit' resolves to citrate", cit_alias is not None)

        pqq_alias = kb.lookup("methoxatin")
        all_pass &= check("'methoxatin' resolves to PQQ", pqq_alias is not None)

        return all_pass

    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


# ============================================================================
# Test 2: ChemicalFeatures Pharmacophore Perception
# ============================================================================
def test_chemicalfeatures_real_ligands():
    separator("ChemicalFeatures - Real Ligand Pharmacophores")

    from ligand_features import get_chemicalfeatures_pharmacophore

    all_pass = True

    test_cases = [
        {
            "name": "Aspirin",
            "smiles": "CC(=O)Oc1ccccc1C(O)=O",
            "expect_donors": True,   # carboxylic acid OH
            "expect_acceptors": True, # carbonyls, ether O
            "expect_aromatic": True,  # benzene ring
            "min_features": 3,
        },
        {
            "name": "Ethanol",
            "smiles": "CCO",
            "expect_donors": True,    # OH
            "expect_acceptors": True,  # OH (also acceptor)
            "expect_aromatic": False,
            "min_features": 1,
        },
        {
            "name": "Benzene",
            "smiles": "c1ccccc1",
            "expect_donors": False,
            "expect_acceptors": False,
            "expect_aromatic": True,
            "min_features": 1,
        },
        {
            "name": "Imidazole (His sidechain analog)",
            "smiles": "c1c[nH]cn1",
            "expect_donors": True,    # NH
            "expect_acceptors": True,  # lone pair N
            "expect_aromatic": True,
            "min_features": 2,
        },
        {
            "name": "Acetic acid",
            "smiles": "CC(O)=O",
            "expect_donors": True,
            "expect_acceptors": True,
            "expect_aromatic": False,
            "min_features": 2,
        },
        {
            "name": "Pyridine",
            "smiles": "c1ccncc1",
            "expect_donors": False,    # no NH
            "expect_acceptors": True,   # N lone pair
            "expect_aromatic": True,
            "min_features": 2,
        },
    ]

    for tc in test_cases:
        print(f"\n  --- {tc['name']} ({tc['smiles']}) ---")
        features = get_chemicalfeatures_pharmacophore(tc["smiles"])

        types_found = set(f["type"] for f in features)
        elements_found = set(f["element"] for f in features)
        print(f"  Feature types: {types_found}")
        print(f"  Elements: {elements_found}")
        print(f"  Total features: {len(features)}")
        for f in features:
            print(f"    {f['atom_name']} ({f['element']}): {f['type']}"
                  + (f" coords={[round(c,1) for c in f['coords']]}" if f.get("coords") else ""))

        all_pass &= check(f"{tc['name']}: enough features",
                          len(features) >= tc["min_features"],
                          f"got {len(features)}, need >= {tc['min_features']}")

        if tc["expect_donors"]:
            all_pass &= check(f"{tc['name']}: has donor features",
                              "donor" in types_found,
                              f"types={types_found}")
        if tc["expect_acceptors"]:
            all_pass &= check(f"{tc['name']}: has acceptor features",
                              "acceptor" in types_found,
                              f"types={types_found}")
        if tc["expect_aromatic"]:
            all_pass &= check(f"{tc['name']}: has aromatic features",
                              "aromatic" in types_found,
                              f"types={types_found}")

    return all_pass


# ============================================================================
# Test 3: analyze_ligand_features End-to-End
# ============================================================================
def test_analyze_e2e():
    separator("analyze_ligand_features - End-to-End")

    from ligand_features import analyze_ligand_features

    all_pass = True
    temp_dir = tempfile.mkdtemp(prefix="test_e2e_")

    try:
        # Case 1: Known ligand (KB hit)
        print(f"\n  --- Case 1: PQQ + CA (should hit KB) ---")
        result = analyze_ligand_features("pqq", metal="CA")
        print(f"  Source: {result['source']}")
        print(f"  Coordination donors: {result['coordination_donors']}")
        print(f"  Max denticity: {result['max_denticity']}")
        print(f"  Compatibility: {result['compatibility_score']:.2f}")
        print(f"  Evidence count: {result['evidence_count']}")
        print(f"  Coordination mode: {result['coordination_mode']}")
        print(f"  Notes: {result['notes'][:80]}...")

        all_pass &= check("PQQ source is knowledge_base",
                          result["source"] == "knowledge_base")
        all_pass &= check("PQQ has coordination donors",
                          len(result["coordination_donors"]) > 0)
        all_pass &= check("PQQ compatibility > 0.5 (hard-hard match)",
                          result["compatibility_score"] > 0.5,
                          f"score={result['compatibility_score']:.2f}")

        # Case 2: Citrate + Tb
        print(f"\n  --- Case 2: Citrate + TB (should hit KB) ---")
        result2 = analyze_ligand_features("citrate",
                                          smiles="OC(=O)CC(O)(CC(O)=O)C(O)=O",
                                          metal="TB")
        print(f"  Source: {result2['source']}")
        print(f"  Coordination donors: {result2['coordination_donors']}")
        print(f"  Max denticity: {result2['max_denticity']}")
        print(f"  Compatibility: {result2['compatibility_score']:.2f}")
        print(f"  Features: {len(result2['features'])} total")

        all_pass &= check("Citrate source is knowledge_base",
                          result2["source"] == "knowledge_base")
        all_pass &= check("Citrate tridentate for Tb",
                          result2["max_denticity"] >= 3,
                          f"denticity={result2['max_denticity']}")

        # Case 3: Unknown ligand with SMILES (should fall to ChemicalFeatures)
        print(f"\n  --- Case 3: Ibuprofen (not in KB, should use ChemicalFeatures) ---")
        ibu_smiles = "CC(C)Cc1ccc(cc1)C(C)C(O)=O"
        result3 = analyze_ligand_features("ibuprofen",
                                          smiles=ibu_smiles,
                                          metal="CA")
        print(f"  Source: {result3['source']}")
        print(f"  Coordination donors: {result3['coordination_donors']}")
        print(f"  Features: {len(result3['features'])}")
        for f in result3["features"]:
            print(f"    {f['atom_name']} ({f['element']}): {f['type']}"
                  + (" [COORD]" if f["is_coordination_donor"] else ""))

        all_pass &= check("Ibuprofen not from KB",
                          result3["source"] in ("chemicalfeatures", "geometry_filter"),
                          f"source={result3['source']}")
        all_pass &= check("Ibuprofen has features",
                          len(result3["features"]) > 0)
        # Ibuprofen has a carboxylic acid - should find O donors
        donor_elements = set(f["element"] for f in result3["features"]
                           if f["is_coordination_donor"])
        all_pass &= check("Ibuprofen carboxylate O identified as donor",
                          "O" in donor_elements,
                          f"donor elements={donor_elements}")

        # Case 4: Unknown ligand, no SMILES (should skip gracefully)
        print(f"\n  --- Case 4: Unknown ligand, no SMILES ---")
        result4 = analyze_ligand_features("totally_fake_ligand_xyz")
        print(f"  Source: {result4['source']}")
        print(f"  Features: {len(result4['features'])}")
        all_pass &= check("Unknown ligand returns gracefully",
                          result4["source"] in ("unknown", "geometry_filter", "chemicalfeatures") or
                          len(result4["features"]) == 0)

        # Case 5: HSAB compatibility check
        print(f"\n  --- Case 5: HSAB Compatibility Scoring ---")
        # Hard acid (Ca2+) + hard base (citrate) = good match
        hard_hard = analyze_ligand_features("citrate", metal="CA")
        # Soft acid-like (hypothetical) - Cu2+ is borderline
        borderline = analyze_ligand_features("citrate", metal="CU")
        print(f"  Citrate + CA compatibility: {hard_hard['compatibility_score']:.2f}")
        print(f"  Citrate + CU compatibility: {borderline['compatibility_score']:.2f}")
        all_pass &= check("Hard-hard (citrate+CA) >= borderline (citrate+CU)",
                          hard_hard["compatibility_score"] >= borderline["compatibility_score"],
                          f"CA={hard_hard['compatibility_score']:.2f} vs CU={borderline['compatibility_score']:.2f}")

    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

    return all_pass


# ============================================================================
# Test 4: KB Growth and Evidence Merging
# ============================================================================
def test_kb_growth_chemistry():
    separator("KB Growth - Evidence Merging Chemistry")

    from ligand_features import LigandKnowledgeBase

    temp_dir = tempfile.mkdtemp(prefix="test_growth_")
    all_pass = True

    try:
        kb = LigandKnowledgeBase(base_dir=temp_dir)

        # Simulate scaffold_search discovering a ligand in two PDB structures
        print(f"\n  Recording evidence from PDB 1ABC (ZN + novel_chelator)...")
        kb.record_evidence(
            ligand_id="novel_chelator",
            pdb_id="1ABC",
            metal="ZN",
            coordination_data={
                "coordination_number": 4,
                "ligand_donors": ["N1@2.05", "N3@2.10", "O5@2.20"],
                "protein_donors": ["A50:HIS NE2@2.08", "A55:CYS SG@2.30"],
                "geometry": "tetrahedral",
                "resolution": 1.8,
            },
            smiles="c1ncc(O)cn1",
        )

        print(f"  Recording second evidence from PDB 2XYZ (ZN + novel_chelator)...")
        kb.record_evidence(
            ligand_id="novel_chelator",
            pdb_id="2XYZ",
            metal="ZN",
            coordination_data={
                "coordination_number": 4,
                "ligand_donors": ["N1@2.08", "N3@2.12"],  # Only 2 donors this time
                "protein_donors": ["B30:HIS NE2@2.10", "B35:CYS SG@2.28", "B40:HIS NE2@2.15"],
                "geometry": "tetrahedral",
                "resolution": 2.0,
            },
        )

        meta = kb.lookup("novel_chelator", metal="ZN")
        all_pass &= check("Novel chelator found after recording", meta is not None)

        if meta:
            coord = meta.get("coordination_summary", {}).get("ZN", {})
            donors = coord.get("donors", [])
            print(f"  Merged donors: {donors}")
            print(f"  Evidence count: {coord.get('evidence_count')}")
            print(f"  Max denticity: {coord.get('denticity')}")

            # N1 appears in both PDBs, N3 in both, O5 only in first
            all_pass &= check("N1 in merged donors (appears in both PDBs)",
                              "N1" in donors, str(donors))
            all_pass &= check("N3 in merged donors (appears in both PDBs)",
                              "N3" in donors, str(donors))
            all_pass &= check("O5 in merged donors (appears in one PDB)",
                              "O5" in donors, str(donors))
            all_pass &= check("Evidence count is 2",
                              coord.get("evidence_count") == 2)

        # Verify PDB evidence files exist
        evidence_dir = os.path.join(temp_dir, "ligands", "novel_chelator", "pdb_evidence")
        if os.path.isdir(evidence_dir):
            files = os.listdir(evidence_dir)
            print(f"  Evidence files: {files}")
            all_pass &= check("Two evidence files created",
                              len(files) == 2, str(files))

    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

    return all_pass


# ============================================================================
# Test 5: Real-world ligand chemistry validation
# ============================================================================
def test_real_world_chemistry():
    separator("Real-World Ligand Chemistry")

    from ligand_features import analyze_ligand_features

    all_pass = True

    # Gluconate + Ca - common in biochemistry
    print(f"\n  --- Gluconate + CA ---")
    gluconate_smiles = "OCC(O)C(O)C(O)C(O)C(=O)[O-]"
    r = analyze_ligand_features("gluconate", smiles=gluconate_smiles, metal="CA")
    print(f"  Source: {r['source']}, Features: {len(r['features'])}")
    print(f"  Coord donors: {r['coordination_donors']}")
    print(f"  Compatibility: {r['compatibility_score']:.2f}")
    all_pass &= check("Gluconate has features",
                      len(r["features"]) > 0)
    all_pass &= check("Gluconate has O donors for Ca",
                      any(f["element"] == "O" and f["is_coordination_donor"]
                          for f in r["features"]))

    # Histidine sidechain analog (imidazole) + Zn
    print(f"\n  --- Imidazole + ZN ---")
    r2 = analyze_ligand_features("imidazole", smiles="c1c[nH]cn1", metal="ZN")
    print(f"  Source: {r2['source']}, Features: {len(r2['features'])}")
    print(f"  Coord donors: {r2['coordination_donors']}")
    all_pass &= check("Imidazole has N donor for Zn",
                      any(f["element"] == "N" and f["is_coordination_donor"]
                          for f in r2["features"]),
                      f"donors={r2['coordination_donors']}")

    # Catechol + Fe (classic chelator)
    print(f"\n  --- Catechol + FE ---")
    r3 = analyze_ligand_features("catechol", smiles="Oc1ccccc1O", metal="FE")
    print(f"  Source: {r3['source']}, Features: {len(r3['features'])}")
    print(f"  Coord donors: {r3['coordination_donors']}")
    all_pass &= check("Catechol has O donors",
                      any(f["element"] == "O" for f in r3["features"]))

    # Bipyridine + Ru (classic N-N chelator)
    print(f"\n  --- 2,2'-Bipyridine + RU ---")
    bipy_smiles = "c1ccnc(-c2ccccn2)c1"
    r4 = analyze_ligand_features("bipyridine", smiles=bipy_smiles, metal="RU")
    print(f"  Source: {r4['source']}, Features: {len(r4['features'])}")
    print(f"  Coord donors: {r4['coordination_donors']}")
    all_pass &= check("Bipyridine has N donors",
                      any(f["element"] == "N" for f in r4["features"]),
                      f"features={[(f['atom_name'], f['type']) for f in r4['features']]}")

    # Thiolate + soft metal (Au)
    print(f"\n  --- Methanethiol + AU (soft-soft) ---")
    r5 = analyze_ligand_features("methanethiol", smiles="CS", metal="AU")
    print(f"  Source: {r5['source']}, Features: {len(r5['features'])}")
    print(f"  Coord donors: {r5['coordination_donors']}")
    # S should be identified as a donor, especially for soft Au
    all_pass &= check("Thiol S identified as feature",
                      any(f["element"] == "S" for f in r5["features"]),
                      f"features={[(f['atom_name'], f['element'], f['type']) for f in r5['features']]}")

    return all_pass


# ============================================================================
# MAIN
# ============================================================================
if __name__ == "__main__":
    print("=" * 60)
    print("  Ligand Features Chemistry Validation")
    print("=" * 60)

    tests = [
        ("Seed KB Chemistry", test_seed_kb_chemistry),
        ("ChemicalFeatures Perception", test_chemicalfeatures_real_ligands),
        ("E2E analyze_ligand_features", test_analyze_e2e),
        ("KB Growth & Merging", test_kb_growth_chemistry),
        ("Real-World Chemistry", test_real_world_chemistry),
    ]

    results = {}
    for name, fn in tests:
        try:
            passed = fn()
            results[name] = "PASS" if passed else "FAIL"
        except Exception as e:
            print(f"\n  EXCEPTION: {e}")
            import traceback
            traceback.print_exc()
            results[name] = f"ERROR: {e}"

    separator("FINAL SUMMARY")
    for name, status in results.items():
        icon = "OK" if status == "PASS" else "XX"
        print(f"  [{icon}] {name}: {status}")
