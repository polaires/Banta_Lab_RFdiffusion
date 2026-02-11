#!/usr/bin/env python3
"""
Test Crystal Scaffold Approach for Metal-Ligand Complex Dimer Design

This script tests the enhanced template-based design approach that uses
crystal structure coordinates for metal-ligand complexes.

Usage:
    # Direct Python test (requires dependencies)
    python test_crystal_scaffold_dimer.py --direct

    # API test (requires Docker running on localhost:8000)
    python test_crystal_scaffold_dimer.py --api

    # Both tests
    python test_crystal_scaffold_dimer.py --all
"""

import argparse
import json
import sys
import tempfile
from pathlib import Path


def test_crystal_scaffold_direct():
    """
    Test crystal scaffold approach directly (without API).

    This tests the new design_dimer_from_crystal_structure() function.
    """
    print("\n" + "=" * 70)
    print("Testing Crystal Scaffold Approach (Direct)")
    print("=" * 70)

    try:
        from handler import design_dimer_from_crystal_structure
    except ImportError as e:
        print(f"[ERROR] Could not import handler: {e}")
        print("[INFO] Make sure you're running from the serverless directory")
        return False

    # Test 1: Citrate-Terbium from 6MI5 (Lanmodulin)
    print("\n[Test 1] Citrate-Terbium from PDB 6MI5")
    print("-" * 50)

    result = design_dimer_from_crystal_structure(
        pdb_id="6MI5",
        metal="TB",
        ligand="CIT",
        chain_length="50-70",
        num_designs=1,
        seed=42,
    )

    if result.get("status") == "completed":
        designs = result.get("result", {}).get("designs", [])
        print(f"[SUCCESS] Generated {len(designs)} design(s)")
        for i, d in enumerate(designs):
            print(f"  Design {i+1}: coordination={d.get('metrics', {}).get('coordination', 'N/A')}")

            # Save design for analysis
            pdb_content = d.get("pdb_content", "")
            if pdb_content:
                output_path = Path(tempfile.gettempdir()) / f"crystal_scaffold_test_{i+1}.pdb"
                with open(output_path, "w") as f:
                    f.write(pdb_content)
                print(f"  Saved to: {output_path}")
        return True
    else:
        print(f"[FAILED] {result.get('error', 'Unknown error')}")
        return False


def test_crystal_scaffold_api(base_url: str = "http://localhost:8000"):
    """
    Test crystal scaffold approach via API.

    This tests the interface_metal_ligand_design task with approach="crystal".
    """
    print("\n" + "=" * 70)
    print(f"Testing Crystal Scaffold Approach (API: {base_url})")
    print("=" * 70)

    try:
        import requests
    except ImportError:
        print("[ERROR] requests library not installed")
        print("[INFO] Run: pip install requests")
        return False

    # First check if API is available
    try:
        health = requests.post(
            f"{base_url}/runsync",
            json={"input": {"task": "health"}},
            timeout=5,
        )
        if health.status_code != 200:
            print(f"[ERROR] API health check failed: {health.status_code}")
            return False
        print("[OK] API is running")
    except requests.exceptions.ConnectionError:
        print(f"[ERROR] Could not connect to {base_url}")
        print("[INFO] Make sure Docker is running:")
        print("  cd backend/serverless")
        print("  docker-compose -f docker-compose.local.yml up --build")
        return False

    # Test: Crystal scaffold approach for Citrate-Terbium
    print("\n[Test] Crystal Scaffold Design via API")
    print("-" * 50)

    request_payload = {
        "input": {
            "task": "interface_metal_ligand_design",
            "template_name": "citrate_tb",
            "approach": "crystal",  # Use crystal scaffold approach
            "contig_str": "50-70,/0,50-70",
            "num_designs": 1,
            "seed": 42,
            "validate_coordination": True,
        }
    }

    print(f"[INFO] Request: {json.dumps(request_payload, indent=2)}")

    try:
        response = requests.post(
            f"{base_url}/runsync",
            json=request_payload,
            timeout=300,  # 5 minute timeout for design
        )

        if response.status_code != 200:
            print(f"[ERROR] API returned status {response.status_code}")
            print(f"[INFO] Response: {response.text[:500]}")
            return False

        result = response.json()
        status = result.get("status", "")

        if status == "COMPLETED":
            output = result.get("output", {})
            designs = output.get("result", {}).get("designs", [])
            print(f"[SUCCESS] Generated {len(designs)} design(s)")

            for i, d in enumerate(designs):
                coord = d.get("metrics", {}).get("coordination", "N/A")
                source = d.get("scaffold_source", "N/A")
                print(f"  Design {i+1}: coordination={coord}, source={source}")

                # Save design
                pdb_content = d.get("pdb_content", "")
                if pdb_content:
                    output_path = Path(tempfile.gettempdir()) / f"crystal_api_test_{i+1}.pdb"
                    with open(output_path, "w") as f:
                        f.write(pdb_content)
                    print(f"  Saved to: {output_path}")

            # Return the result for analyze-design
            return result
        else:
            error = result.get("output", {}).get("error", "Unknown error")
            print(f"[FAILED] Status: {status}, Error: {error}")
            return False

    except requests.exceptions.Timeout:
        print("[ERROR] Request timed out after 5 minutes")
        return False
    except Exception as e:
        print(f"[ERROR] {e}")
        return False


def test_scaffold_extraction():
    """
    Test that scaffold extraction includes HETATM records.
    """
    print("\n" + "=" * 70)
    print("Testing Scaffold Extraction (HETATM included)")
    print("=" * 70)

    try:
        from handler import _extract_coordination_scaffold
        from metal_site_fetcher import find_metal_ligand_active_site, _fetch_pdb_content
    except ImportError as e:
        print(f"[ERROR] Could not import: {e}")
        return False

    pdb_id = "1W6S"  # PQQ-Ca from MDH
    metal = "CA"
    ligand = "PQQ"

    print(f"\n[Test] Extracting scaffold from PDB {pdb_id} ({metal}+{ligand})")
    print("-" * 50)

    # Get site info
    site_info = find_metal_ligand_active_site(pdb_id, metal, ligand)
    if not site_info:
        print("[FAILED] Could not find active site")
        return False

    print(f"[OK] Found active site at {site_info['metal_chain']}{site_info['metal_resnum']}")

    # Build protein_coords
    amino_acids = {"GLU", "ASP", "HIS", "CYS", "ASN", "TYR", "MET", "SER", "THR", "GLN"}
    protein_coords = []
    for atom in site_info.get("coordinating_atoms", []):
        res_name = atom.get("res_name", "")
        if res_name in amino_acids:
            protein_coords.append({
                "residue": res_name,
                "chain": atom.get("chain_id", ""),
                "resnum": atom.get("res_seq", 0),
                "atom": atom.get("atom_name", ""),
                "distance": atom.get("distance", 0),
            })

    # Fetch PDB and extract scaffold
    pdb_content = _fetch_pdb_content(pdb_id)
    if not pdb_content:
        print("[FAILED] Could not fetch PDB content")
        return False

    scaffold = _extract_coordination_scaffold(
        pdb_content=pdb_content,
        metal=metal,
        ligand=ligand,
        protein_coords=protein_coords,
        metal_chain=site_info['metal_chain'],
        metal_resnum=site_info['metal_resnum'],
    )

    if not scaffold:
        print("[FAILED] Could not extract scaffold")
        return False

    scaffold_pdb = scaffold.get("scaffold_pdb", "")

    # Check for HETATM
    hetatm_count = sum(1 for line in scaffold_pdb.split('\n') if line.startswith('HETATM'))
    atom_count = sum(1 for line in scaffold_pdb.split('\n') if line.startswith('ATOM'))

    print(f"[INFO] Scaffold has {atom_count} ATOM records, {hetatm_count} HETATM records")

    if hetatm_count > 0:
        print("[SUCCESS] HETATM records included in scaffold")

        # Check for original_hetatm_pdb
        if scaffold.get("original_hetatm_pdb"):
            print("[SUCCESS] original_hetatm_pdb stored for geometry restoration")
        else:
            print("[WARNING] original_hetatm_pdb not stored")

        return True
    else:
        print("[FAILED] No HETATM records in scaffold")
        return False


def main():
    parser = argparse.ArgumentParser(description="Test Crystal Scaffold Dimer Design")
    parser.add_argument("--direct", action="store_true", help="Run direct Python test")
    parser.add_argument("--api", action="store_true", help="Run API test")
    parser.add_argument("--extraction", action="store_true", help="Test scaffold extraction")
    parser.add_argument("--all", action="store_true", help="Run all tests")
    parser.add_argument("--url", default="http://localhost:8000", help="API base URL")

    args = parser.parse_args()

    # Default to extraction test if no args
    if not any([args.direct, args.api, args.extraction, args.all]):
        args.extraction = True

    results = {}

    if args.extraction or args.all:
        results["extraction"] = test_scaffold_extraction()

    if args.direct or args.all:
        results["direct"] = test_crystal_scaffold_direct()

    if args.api or args.all:
        results["api"] = test_crystal_scaffold_api(args.url)

    # Summary
    print("\n" + "=" * 70)
    print("Test Summary")
    print("=" * 70)
    for test_name, passed in results.items():
        status = "PASSED" if passed else "FAILED"
        print(f"  {test_name}: {status}")

    # Return non-zero if any test failed
    if all(results.values()):
        print("\n[ALL TESTS PASSED]")
        return 0
    else:
        print("\n[SOME TESTS FAILED]")
        return 1


if __name__ == "__main__":
    sys.exit(main())
