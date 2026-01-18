#!/usr/bin/env python3
"""
Master Script for Testing Four Topology-Correct Dimer Design Approaches

This script coordinates the execution of all four approaches for designing
azobenzene dimers with ligand at the interface:
1. Approach 1: Asymmetric Monomer Design (Offset Origin)
2. Approach 2: Binder-to-Binder Design (Sequential)
3. Approach 3: Template-Based Design
4. Approach 4: Domain-Based Cleavage (Better Selection)

Usage:
    python run_all_approaches.py --approach 1  # Run specific approach
    python run_all_approaches.py --all         # Run all approaches
    python run_all_approaches.py --analyze     # Analyze existing outputs
"""

import os
import sys
import json
import subprocess
import argparse
import requests
import time
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from datetime import datetime

# Add parent directories to path for imports
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "backend" / "serverless"))

try:
    from topology_validation import validate_dimer_topology, filter_valid_designs
    TOPOLOGY_VALIDATION_AVAILABLE = True
except ImportError:
    TOPOLOGY_VALIDATION_AVAILABLE = False
    print("Warning: topology_validation module not available")

# Configuration paths
SETUP_DIR = SCRIPT_DIR / "setup"
CONFIGS_DIR = SETUP_DIR / "configs"
OUTPUTS_DIR = SCRIPT_DIR / "outputs"

# RFD3 API endpoint (Docker container at localhost:8000)
RFD3_API_URL = "http://localhost:8000/runsync"

# Azobenzene SMILES
AZOBENZENE_SMILES = "c1ccc(cc1)N=Nc2ccccc2"

# RDKit atom naming for azobenzene:
# Ring 1: C1, C2, C3, C4, C5, C6 (first phenyl, attached to N7)
# Azo bridge: N7, N8
# Ring 2: C9, C10, C11, C12, C13, C14 (second phenyl, attached to N8)
AZOBENZENE_RING1 = "C1,C2,C3,C4,C5,C6"
AZOBENZENE_RING2 = "C9,C10,C11,C12,C13,C14"
AZOBENZENE_AZO_BRIDGE = "N7,N8"


def log(message: str, level: str = "INFO"):
    """Print timestamped log message."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] [{level}] {message}")


def check_api_health() -> bool:
    """Check if the RFD3 API is available."""
    try:
        response = requests.post(
            RFD3_API_URL,
            json={"input": {"task": "health"}},
            timeout=10
        )
        return response.status_code == 200
    except Exception as e:
        log(f"API health check failed: {e}", "ERROR")
        return False


def run_rfd3_design(
    config_path: str,
    output_dir: str,
    n_batches: int = 2,
    batch_size: int = 5,
    use_symmetry: bool = False,
    dry_run: bool = False
) -> bool:
    """
    Run RFD3 design via HTTP API.

    Args:
        config_path: Path to JSON config file
        output_dir: Output directory for designs
        n_batches: Number of batches
        batch_size: Designs per batch
        use_symmetry: Whether to use symmetry sampler
        dry_run: If True, print request but don't execute

    Returns:
        True if successful, False otherwise
    """
    # Load config
    try:
        with open(config_path, 'r') as f:
            config = json.load(f)
    except Exception as e:
        log(f"Failed to load config {config_path}: {e}", "ERROR")
        return False

    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Process each design in the config
    total_designs = 0
    success_count = 0

    for design_name, design_config in config.items():
        # Skip comment fields
        if design_name.startswith('_'):
            continue

        log(f"Running design: {design_name}")

        # Build API request
        api_request = {
            "input": {
                "task": "rfd3",
                "contig": design_config.get("contig", "60-80"),
                "num_designs": batch_size * n_batches,
            }
        }

        # Add optional parameters
        if "ligand" in design_config:
            api_request["input"]["ligand"] = design_config["ligand"]

        if "input" in design_config:
            # Read input PDB and include as string
            input_pdb_path = PROJECT_ROOT / design_config["input"]
            if input_pdb_path.exists():
                with open(input_pdb_path, 'r') as f:
                    api_request["input"]["input_pdb"] = f.read()

        if "ori_token" in design_config:
            api_request["input"]["ori_token"] = design_config["ori_token"]

        if "select_exposed" in design_config:
            api_request["input"]["select_exposed"] = design_config["select_exposed"]

        if "select_buried" in design_config:
            api_request["input"]["select_buried"] = design_config["select_buried"]

        if "select_partially_buried" in design_config:
            api_request["input"]["select_partially_buried"] = design_config["select_partially_buried"]

        if "select_hotspots" in design_config:
            api_request["input"]["select_hotspots"] = design_config["select_hotspots"]

        if "partial_t" in design_config:
            api_request["input"]["partial_t"] = design_config["partial_t"]

        if use_symmetry:
            api_request["input"]["symmetry"] = "C2"

        log(f"API Request: {json.dumps(api_request, indent=2)[:500]}...")

        if dry_run:
            log("DRY RUN - request not sent", "DEBUG")
            total_designs += 1
            success_count += 1
            continue

        try:
            # Send request to API
            response = requests.post(
                RFD3_API_URL,
                json=api_request,
                timeout=3600  # 1 hour timeout
            )

            if response.status_code != 200:
                log(f"API request failed: {response.status_code} - {response.text}", "ERROR")
                total_designs += 1
                continue

            result = response.json()

            # Check for errors in response
            if "error" in result:
                log(f"Design error: {result['error']}", "ERROR")
                total_designs += 1
                continue

            # Save output PDBs
            output = result.get("output", {})
            pdb_outputs = output.get("pdb_outputs", [])

            for i, pdb_content in enumerate(pdb_outputs):
                output_path = Path(output_dir) / f"{design_name}_{i}.pdb"
                with open(output_path, 'w') as f:
                    f.write(pdb_content)
                log(f"Saved: {output_path}")

            total_designs += 1
            success_count += 1
            log(f"Design {design_name} completed: {len(pdb_outputs)} structures")

        except requests.Timeout:
            log(f"Design {design_name} timed out", "ERROR")
            total_designs += 1
        except Exception as e:
            log(f"Design {design_name} error: {e}", "ERROR")
            total_designs += 1

    log(f"Completed {success_count}/{total_designs} designs")
    return success_count > 0


def apply_c2_symmetry(input_pdb: str, output_pdb: str, axis: str = "x") -> bool:
    """
    Apply C2 symmetry rotation to create dimer from monomer.

    For azobenzene where N=N bond is along X axis:
    - Rotate 180 degrees around X axis (the C2 symmetry axis of the ligand)
    - This transforms: y -> -y, z -> -z, x -> x
    - Chain A stays on +Y side, Chain B goes to -Y side

    Args:
        input_pdb: Path to monomer PDB
        output_pdb: Path to output dimer PDB
        axis: Rotation axis - "x", "y", or "z"

    Returns:
        True if successful
    """
    try:
        import numpy as np

        # Read input PDB
        with open(input_pdb, 'r') as f:
            lines = f.readlines()

        # Parse atoms
        atoms_a = []
        ligand_lines = []

        for line in lines:
            if line.startswith('ATOM'):
                atoms_a.append(line)
            elif line.startswith('HETATM'):
                ligand_lines.append(line)

        # C2 rotation functions for different axes
        def rotate_c2_x(x, y, z):
            """180 degrees around X axis: y -> -y, z -> -z"""
            return x, -y, -z

        def rotate_c2_y(x, y, z):
            """180 degrees around Y axis: x -> -x, z -> -z"""
            return -x, y, -z

        def rotate_c2_z(x, y, z):
            """180 degrees around Z axis: x -> -x, y -> -y"""
            return -x, -y, z

        # Select rotation function based on axis
        rotate_func = {
            "x": rotate_c2_x,
            "y": rotate_c2_y,
            "z": rotate_c2_z,
        }.get(axis.lower(), rotate_c2_x)

        # Generate chain B by rotating chain A
        atoms_b = []
        for line in atoms_a:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            new_x, new_y, new_z = rotate_func(x, y, z)

            # Create new line with chain B and rotated coordinates
            new_line = line[:21] + 'B' + line[22:30]
            new_line += f"{new_x:8.3f}{new_y:8.3f}{new_z:8.3f}"
            new_line += line[54:]

            atoms_b.append(new_line)

        # Write output
        with open(output_pdb, 'w') as f:
            f.write(f"REMARK  C2 symmetric dimer generated from monomer\n")
            f.write(f"REMARK  Chain A: original, Chain B: C2 rotated around {axis.upper()} axis\n")

            # Write chain A
            for line in atoms_a:
                f.write(line)

            f.write("TER\n")

            # Write chain B
            for line in atoms_b:
                f.write(line)

            f.write("TER\n")

            # Write ligand (keep in place - at interface)
            for line in ligand_lines:
                f.write(line)

            f.write("END\n")

        log(f"Created C2 dimer (axis={axis}): {output_pdb}")
        return True

    except Exception as e:
        log(f"Failed to apply C2 symmetry: {e}", "ERROR")
        return False


def run_interface_ligand_design(
    approach: str,
    chain_length: str,
    num_designs: int,
    ori_offset: Optional[List[float]] = None,
    exposed_atoms: Optional[str] = None,
    side: str = "left",
    design_name: str = "design",
    chain_a_pdb: Optional[str] = None,
    dry_run: bool = False
) -> Dict[str, Any]:
    """
    Run interface_ligand_design API task with SMILES input.

    This is the correct way to design with RASA conditioning.
    The API handles ligand preparation and atom naming internally.

    Args:
        approach: "asymmetric", "sequential", or "full"
        chain_length: e.g., "60-80"
        num_designs: Number of designs
        ori_offset: [x, y, z] offset from ligand (optional)
        exposed_atoms: Atoms to keep exposed (e.g., "C7,C8,C9")
        side: "left" or "right" for asymmetric approach
        design_name: Name prefix for output files
        chain_a_pdb: PDB content for sequential approach (required if approach="sequential")
        dry_run: Print request without executing

    Returns:
        Result dict with status, designs, and metrics
    """
    # Build API request
    api_request = {
        "input": {
            "task": "interface_ligand_design",
            "ligand_smiles": AZOBENZENE_SMILES,
            "approach": approach,
            "chain_length": chain_length,
            "num_designs": num_designs,
            "side": side,
        }
    }

    # Optional parameters
    if ori_offset:
        api_request["input"]["ori_offset"] = ori_offset
        api_request["input"]["use_ori_token"] = True

    if exposed_atoms:
        api_request["input"]["exposed_atoms"] = exposed_atoms

    if chain_a_pdb and approach == "sequential":
        api_request["input"]["chain_a_pdb"] = chain_a_pdb

    log(f"API Request: {json.dumps(api_request, indent=2)[:800]}...")

    if dry_run:
        log("DRY RUN - request not sent", "DEBUG")
        return {"status": "dry_run", "designs": []}

    try:
        response = requests.post(
            RFD3_API_URL,
            json=api_request,
            timeout=3600  # 1 hour timeout
        )

        if response.status_code != 200:
            log(f"API request failed: {response.status_code} - {response.text[:500]}", "ERROR")
            return {"status": "failed", "error": response.text}

        result = response.json()

        # RunPod API wraps output in "output" key
        if "output" in result:
            result = result["output"]

        # Check for errors
        if result.get("status") == "failed":
            log(f"Design failed: {result.get('error')}", "ERROR")
            return result

        return result

    except requests.Timeout:
        log(f"Design {design_name} timed out", "ERROR")
        return {"status": "failed", "error": "timeout"}
    except Exception as e:
        log(f"Design {design_name} error: {e}", "ERROR")
        return {"status": "failed", "error": str(e)}


def run_approach_1(
    n_batches: int = 2,
    batch_size: int = 5,
    dry_run: bool = False
) -> Dict[str, Any]:
    """
    Run Approach 1: Asymmetric Monomer Design with ori_token offset.

    Uses interface_ligand_design API with SMILES input.
    Designs one-sided binders, then applies C2 symmetry post-processing.
    """
    log("=" * 60)
    log("APPROACH 1: Asymmetric Monomer Design (Offset Origin)")
    log("=" * 60)

    results = {
        "approach": 1,
        "name": "offset_origin",
        "status": "running",
        "designs": [],
        "dimers": []
    }

    output_dir = OUTPUTS_DIR / "approach_1_offset"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Test different configurations
    # Note: ori_token doesn't reliably position proteins on opposite sides
    # Use select_exposed as the primary mechanism for one-sided binding
    #
    # Key insight: We need the protein to bind ONLY one half of the ligand
    # so that C2 rotation creates a true dimer with separable chains
    offset_configs = [
        # Expose Ring 2 (C9-C14) - protein binds Ring 1 side
        {"name": "expose_ring2_60_80", "offset": None, "length": "60-80",
         "exposed": AZOBENZENE_RING2, "buried": None},
        {"name": "expose_ring2_80_100", "offset": None, "length": "80-100",
         "exposed": AZOBENZENE_RING2, "buried": None},
        # Expose Ring 1 + azo - protein binds Ring 2 only
        {"name": "expose_ring1_azo_60_80", "offset": None, "length": "60-80",
         "exposed": f"{AZOBENZENE_RING1},{AZOBENZENE_AZO_BRIDGE}", "buried": None},
        # Shorter chains for more compact binding
        {"name": "expose_ring2_40_60", "offset": None, "length": "40-60",
         "exposed": AZOBENZENE_RING2, "buried": None},
    ]

    num_designs_each = batch_size * n_batches
    all_designs = []

    for config in offset_configs:
        log(f"Running design: {config['name']}")

        result = run_interface_ligand_design(
            approach="asymmetric",
            chain_length=config["length"],
            num_designs=num_designs_each,
            ori_offset=config.get("offset"),
            exposed_atoms=config.get("exposed", AZOBENZENE_RING2),
            side="left",  # Bind from left (Ring 1 side)
            design_name=config["name"],
            dry_run=dry_run
        )

        if result.get("status") == "completed":
            # Extract designs from result
            api_result = result.get("result", {})
            designs = api_result.get("designs", [])

            for i, design in enumerate(designs):
                pdb_content = design.get("pdb_content", "")
                if pdb_content:
                    output_path = output_dir / f"{config['name']}_{i}.pdb"
                    with open(output_path, 'w') as f:
                        f.write(pdb_content)
                    all_designs.append({
                        "path": str(output_path),
                        "config": config["name"],
                        "metrics": design.get("metrics", {})
                    })
                    log(f"Saved: {output_path}")

            log(f"Design {config['name']} completed: {len(designs)} structures")
        else:
            log(f"Design {config['name']} failed: {result.get('error', 'unknown')}", "ERROR")

    results["designs"] = all_designs

    # Apply C2 symmetry to create dimers
    if not dry_run and all_designs:
        dimer_dir = output_dir / "dimers"
        dimer_dir.mkdir(exist_ok=True)

        for design in all_designs:
            monomer_pdb = Path(design["path"])
            dimer_pdb = dimer_dir / f"{monomer_pdb.stem}_dimer.pdb"
            if apply_c2_symmetry(str(monomer_pdb), str(dimer_pdb)):
                results["dimers"].append(str(dimer_pdb))

    results["status"] = "completed" if all_designs else "failed"
    log(f"Completed {len(all_designs)} designs")
    return results


def run_approach_2(
    n_batches: int = 2,
    batch_size: int = 5,
    dry_run: bool = False
) -> Dict[str, Any]:
    """
    Run Approach 2: Binder-to-Binder Design (Full Workflow).

    Uses interface_ligand_design API with "full" approach which:
    1. Designs asymmetric Chain A (one side of ligand)
    2. Designs complementary Chain B against Chain A + ligand
    """
    log("=" * 60)
    log("APPROACH 2: Binder-to-Binder Design (Full Workflow)")
    log("=" * 60)

    results = {
        "approach": 2,
        "name": "binder_to_binder",
        "status": "running",
        "designs": [],
        "dimers": []
    }

    output_dir = OUTPUTS_DIR / "approach_2_sequential"
    output_dir.mkdir(parents=True, exist_ok=True)

    num_designs = batch_size * n_batches
    all_designs = []

    # Use "full" approach which handles both chains automatically
    log("Running full dimer design workflow...")

    result = run_interface_ligand_design(
        approach="full",
        chain_length="50-70",  # Shorter chains for sequential design
        num_designs=num_designs,
        ori_offset=[12.0, 0.0, 0.0],  # Standard offset
        exposed_atoms=AZOBENZENE_RING2,  # Keep Ring 2 exposed (C9-C14)
        side="left",  # Start binding from left
        design_name="binder_full",
        dry_run=dry_run
    )

    if result.get("status") == "completed":
        # Extract designs from result
        api_result = result.get("result", {})
        designs = api_result.get("designs", [])

        for i, design in enumerate(designs):
            pdb_content = design.get("pdb_content", "")
            if pdb_content:
                output_path = output_dir / f"dimer_full_{i}.pdb"
                with open(output_path, 'w') as f:
                    f.write(pdb_content)
                all_designs.append({
                    "path": str(output_path),
                    "metrics": design.get("metrics", {})
                })
                log(f"Saved: {output_path}")

        log(f"Full workflow completed: {len(designs)} dimers")
    else:
        log(f"Full workflow failed: {result.get('error', 'unknown')}", "ERROR")

    results["designs"] = all_designs
    results["dimers"] = [d["path"] for d in all_designs]
    results["status"] = "completed" if all_designs else "failed"

    log(f"Completed {len(all_designs)} designs")
    return results


def run_approach_3(
    n_batches: int = 2,
    batch_size: int = 5,
    dry_run: bool = False
) -> Dict[str, Any]:
    """
    Run Approach 3: Template-Based Design.

    Uses natural dimers from PDB as templates, redesigns binding pocket.
    """
    log("=" * 60)
    log("APPROACH 3: Template-Based Design")
    log("=" * 60)

    results = {
        "approach": 3,
        "name": "template_based",
        "status": "running",
        "designs": []
    }

    config_path = CONFIGS_DIR / "approach3_template.json"
    output_dir = OUTPUTS_DIR / "approach_3_template"

    # Check if template files exist
    templates_dir = SETUP_DIR / "templates"
    if not templates_dir.exists():
        log("Templates directory not found. Creating placeholder...", "WARNING")
        templates_dir.mkdir(parents=True, exist_ok=True)

        # Create README for templates
        readme = templates_dir / "README.md"
        readme.write_text("""# Template PDBs for Approach 3

Download and prepare these template dimers:

1. **2M0Z** - Photoswitchable PDZ domain with azobenzene derivative
   - Download from RCSB PDB
   - Replace ligand with azobenzene
   - Save as 2m0z_with_azo.pdb

2. **7LMN** - Immunoglobulin light chain homodimer
   - Download from RCSB PDB
   - Add azobenzene at interface
   - Save as 7lmn_with_azo.pdb

3. **4DM6** - RARb LBD homodimer
   - Download from RCSB PDB
   - Replace TTNPB with azobenzene
   - Save as 4dm6_with_azo.pdb
""")

        results["status"] = "templates_needed"
        results["message"] = "Template PDBs need to be downloaded and prepared"
        return results

    # Run RFD3 design
    success = run_rfd3_design(
        str(config_path),
        str(output_dir),
        n_batches=n_batches,
        batch_size=batch_size,
        dry_run=dry_run
    )

    if not success:
        results["status"] = "failed"
        return results

    if not dry_run:
        design_pdbs = list(output_dir.glob("*.pdb"))
        results["designs"] = [str(p) for p in design_pdbs]

    results["status"] = "completed"
    return results


def run_approach_4(
    n_batches: int = 2,
    batch_size: int = 5,
    dry_run: bool = False
) -> Dict[str, Any]:
    """
    Run Approach 4: Domain-Based Cleavage with Better Selection.

    Designs monomer around ligand, then finds valid cleavage sites
    using topology validation.
    """
    log("=" * 60)
    log("APPROACH 4: Domain-Based Cleavage (Better Selection)")
    log("=" * 60)

    results = {
        "approach": 4,
        "name": "better_cleavage",
        "status": "running",
        "monomers": [],
        "dimers": []
    }

    # First, design monomers around the ligand (full wrap, not offset)
    monomer_config = {
        "cleavable_monomer": {
            "input": "experiments/azobenzene_dimer/setup/azobenzene_oriented.pdb",
            "contig": "80-120",
            "ligand": "AZO",
            "select_partially_buried": {"AZO": "ALL"}
        }
    }

    config_path = CONFIGS_DIR / "approach4_cleavable.json"
    with open(config_path, 'w') as f:
        json.dump(monomer_config, f, indent=2)

    output_dir = OUTPUTS_DIR / "approach_4_cleavage"

    success = run_rfd3_design(
        str(config_path),
        str(output_dir),
        n_batches=n_batches,
        batch_size=batch_size,
        dry_run=dry_run
    )

    if not success:
        results["status"] = "failed"
        return results

    if not dry_run:
        monomer_pdbs = list(output_dir.glob("*.pdb"))
        results["monomers"] = [str(p) for p in monomer_pdbs]

        # TODO: Apply cleavage algorithm with topology validation
        # This requires the enhanced cleavage_utils.py
        log("Cleavage analysis requires enhanced cleavage_utils.py", "WARNING")

    results["status"] = "completed"
    return results


def analyze_all_outputs() -> Dict[str, Any]:
    """
    Analyze all design outputs from all approaches.

    Runs GNINA scoring and topology validation on all designs.
    """
    log("=" * 60)
    log("ANALYZING ALL DESIGN OUTPUTS")
    log("=" * 60)

    results = {
        "timestamp": datetime.now().isoformat(),
        "approaches": {}
    }

    # Collect all PDB files from each approach
    approach_dirs = {
        1: OUTPUTS_DIR / "approach_1_offset" / "dimers",
        2: OUTPUTS_DIR / "approach_2_sequential" / "step2",
        3: OUTPUTS_DIR / "approach_3_template",
        4: OUTPUTS_DIR / "approach_4_cleavage"
    }

    for approach_num, approach_dir in approach_dirs.items():
        approach_results = {
            "directory": str(approach_dir),
            "designs": [],
            "valid_count": 0,
            "total_count": 0
        }

        if not approach_dir.exists():
            approach_results["status"] = "no_outputs"
            results["approaches"][approach_num] = approach_results
            continue

        pdb_files = list(approach_dir.glob("*.pdb"))
        approach_results["total_count"] = len(pdb_files)

        if TOPOLOGY_VALIDATION_AVAILABLE:
            valid_paths, validations = filter_valid_designs(
                [str(p) for p in pdb_files]
            )

            approach_results["valid_count"] = len(valid_paths)
            approach_results["designs"] = validations
        else:
            approach_results["status"] = "topology_validation_not_available"
            approach_results["designs"] = [str(p) for p in pdb_files]

        results["approaches"][approach_num] = approach_results

    # Summary
    total_valid = sum(
        r.get("valid_count", 0) for r in results["approaches"].values()
    )
    total_designs = sum(
        r.get("total_count", 0) for r in results["approaches"].values()
    )

    results["summary"] = {
        "total_designs": total_designs,
        "valid_designs": total_valid,
        "success_rate": total_valid / total_designs if total_designs > 0 else 0
    }

    log(f"Analysis complete: {total_valid}/{total_designs} designs have valid topology")

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Run topology-correct dimer design approaches"
    )
    parser.add_argument(
        "--approach", "-a",
        type=int,
        choices=[1, 2, 3, 4],
        help="Run specific approach (1-4)"
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Run all approaches"
    )
    parser.add_argument(
        "--analyze",
        action="store_true",
        help="Analyze existing outputs"
    )
    parser.add_argument(
        "--n-batches",
        type=int,
        default=2,
        help="Number of batches per design (default: 2)"
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=5,
        help="Designs per batch (default: 5)"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing"
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Output JSON file for results"
    )

    args = parser.parse_args()

    # Ensure output directories exist
    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)

    results = {"timestamp": datetime.now().isoformat()}

    if args.analyze:
        results = analyze_all_outputs()

    elif args.all:
        log("Running all four approaches...")
        results["approach_1"] = run_approach_1(
            args.n_batches, args.batch_size, args.dry_run
        )
        results["approach_2"] = run_approach_2(
            args.n_batches, args.batch_size, args.dry_run
        )
        results["approach_3"] = run_approach_3(
            args.n_batches, args.batch_size, args.dry_run
        )
        results["approach_4"] = run_approach_4(
            args.n_batches, args.batch_size, args.dry_run
        )

        # Run analysis on completed designs
        if not args.dry_run:
            results["analysis"] = analyze_all_outputs()

    elif args.approach:
        approach_funcs = {
            1: run_approach_1,
            2: run_approach_2,
            3: run_approach_3,
            4: run_approach_4
        }
        results[f"approach_{args.approach}"] = approach_funcs[args.approach](
            args.n_batches, args.batch_size, args.dry_run
        )

    else:
        parser.print_help()
        return

    # Save results
    if args.output:
        output_path = Path(args.output)
    else:
        output_path = OUTPUTS_DIR / f"results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"

    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)

    log(f"Results saved to: {output_path}")

    # Print summary
    log("=" * 60)
    log("SUMMARY")
    log("=" * 60)
    print(json.dumps(results, indent=2, default=str))


if __name__ == "__main__":
    main()
