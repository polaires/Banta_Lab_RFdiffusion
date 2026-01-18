#!/usr/bin/env python3
"""
Test script for lanthanide metal binding template improvements.

Tests the following template/donor combinations:
1. EF-hand with ASP (recommended for tight coordination)
2. EF-hand with GLU (original behavior baseline)
3. EF-hand with MIXED (5 Asp + 3 Glu)
4. C4-symmetric with ASP
5. No template (de novo baseline)

Usage:
    python test_lanthanide_templates.py --backend-url http://localhost:8000
    python test_lanthanide_templates.py --config ef_hand_asp  # Run single config
    python test_lanthanide_templates.py --all  # Run all configs
"""

import argparse
import json
import os
import sys
from datetime import datetime
from pathlib import Path
import requests
from typing import Dict, Any, List, Optional

# Configuration
SCRIPT_DIR = Path(__file__).parent
CONFIG_DIR = SCRIPT_DIR / "setup" / "configs"
OUTPUT_DIR = SCRIPT_DIR / "outputs"

# Available test configurations
TEST_CONFIGS = [
    "ef_hand_asp",
    "ef_hand_glu",
    "ef_hand_mixed",
    "c4_asp",
    "no_template",
]


def load_config(config_name: str) -> Dict[str, Any]:
    """Load a test configuration from JSON file."""
    config_path = CONFIG_DIR / f"{config_name}.json"
    if not config_path.exists():
        raise FileNotFoundError(f"Config not found: {config_path}")
    with open(config_path, "r") as f:
        return json.load(f)


def run_test(config_name: str, backend_url: str, timeout: int = 600) -> Dict[str, Any]:
    """Run a single test configuration against the backend."""
    config = load_config(config_name)
    print(f"\n{'='*60}")
    print(f"Running: {config['name']}")
    print(f"Description: {config['description']}")
    print(f"{'='*60}")

    # Make API request
    url = f"{backend_url}/runsync"
    payload = {"input": config["input"]}

    print(f"Request URL: {url}")
    print(f"Input: {json.dumps(config['input'], indent=2)}")

    try:
        response = requests.post(url, json=payload, timeout=timeout)
        response.raise_for_status()
        result = response.json()
    except requests.exceptions.Timeout:
        return {
            "config_name": config_name,
            "status": "timeout",
            "error": f"Request timed out after {timeout}s"
        }
    except Exception as e:
        return {
            "config_name": config_name,
            "status": "error",
            "error": str(e)
        }

    return {
        "config_name": config_name,
        "config": config,
        "result": result,
        "status": result.get("status", "unknown"),
    }


def analyze_result(test_result: Dict[str, Any]) -> Dict[str, Any]:
    """Analyze test result and extract key metrics."""
    if test_result["status"] != "completed":
        return {
            "success": False,
            "error": test_result.get("error", "Unknown error"),
        }

    result_data = test_result.get("result", {}).get("result", {})
    designs = result_data.get("designs", [])

    if not designs:
        return {
            "success": False,
            "error": "No designs returned",
        }

    # Extract metrics from all designs
    metrics_list = []
    for design in designs:
        metrics = design.get("metrics", {})
        coord_analysis = design.get("coordination_analysis", {})

        # Count carboxylate donors (D and E residues)
        donor_residues = coord_analysis.get("donor_residues", [])
        carboxylate_count = sum(
            1 for r in donor_residues
            if r.get("resname") in ["ASP", "GLU"]
        )

        metrics_list.append({
            "coordination_number": metrics.get("coordination_number", 0),
            "carboxylate_donors": carboxylate_count,
            "chain_a_donors": metrics.get("chain_a_donors", 0),
            "chain_b_donors": metrics.get("chain_b_donors", 0),
            "average_distance": metrics.get("average_distance"),
            "sequence_identity": metrics.get("sequence_identity"),
            "donor_residues": donor_residues,
        })

    # Compute averages
    avg_coord = sum(m["coordination_number"] for m in metrics_list) / len(metrics_list)
    avg_carboxylate = sum(m["carboxylate_donors"] for m in metrics_list) / len(metrics_list)

    return {
        "success": True,
        "num_designs": len(designs),
        "average_coordination": round(avg_coord, 2),
        "average_carboxylate_donors": round(avg_carboxylate, 2),
        "per_design_metrics": metrics_list,
    }


def save_results(
    config_name: str,
    test_result: Dict[str, Any],
    analysis: Dict[str, Any],
    output_dir: Path
) -> Path:
    """Save test results to files."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = output_dir / f"{config_name}_{timestamp}"
    run_dir.mkdir(parents=True, exist_ok=True)

    # Save raw result
    with open(run_dir / "result.json", "w") as f:
        json.dump(test_result, f, indent=2)

    # Save analysis
    with open(run_dir / "analysis.json", "w") as f:
        json.dump(analysis, f, indent=2)

    # Save PDB files
    if test_result["status"] == "completed":
        result_data = test_result.get("result", {}).get("result", {})
        designs = result_data.get("designs", [])
        for i, design in enumerate(designs):
            pdb_content = design.get("pdb_content")
            if pdb_content:
                with open(run_dir / f"design_{i+1}.pdb", "w") as f:
                    f.write(pdb_content)

    return run_dir


def print_summary(results: List[Dict[str, Any]]):
    """Print summary of all test results."""
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"{'Config':<20} {'Status':<12} {'Coord':<8} {'Carbox':<8} {'Notes'}")
    print("-"*80)

    for r in results:
        config_name = r["config_name"]
        status = r["status"]
        analysis = r.get("analysis", {})

        if analysis.get("success"):
            coord = f"{analysis['average_coordination']:.1f}"
            carbox = f"{analysis['average_carboxylate_donors']:.1f}"
            notes = f"{analysis['num_designs']} designs"
        else:
            coord = "-"
            carbox = "-"
            notes = analysis.get("error", "Error")[:30]

        print(f"{config_name:<20} {status:<12} {coord:<8} {carbox:<8} {notes}")

    print("="*80)


def main():
    parser = argparse.ArgumentParser(description="Test lanthanide template improvements")
    parser.add_argument(
        "--backend-url",
        default="http://localhost:8000",
        help="Backend API URL (default: http://localhost:8000)"
    )
    parser.add_argument(
        "--config",
        choices=TEST_CONFIGS,
        help="Run a specific configuration"
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Run all test configurations"
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=600,
        help="Request timeout in seconds (default: 600)"
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List available configurations"
    )

    args = parser.parse_args()

    if args.list:
        print("Available test configurations:")
        for config_name in TEST_CONFIGS:
            config = load_config(config_name)
            print(f"  {config_name}: {config['name']}")
        return 0

    if not args.config and not args.all:
        print("Please specify --config <name> or --all")
        print("Use --list to see available configurations")
        return 1

    # Determine which configs to run
    configs_to_run = TEST_CONFIGS if args.all else [args.config]

    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Run tests
    results = []
    for config_name in configs_to_run:
        test_result = run_test(config_name, args.backend_url, args.timeout)
        analysis = analyze_result(test_result)

        # Save results
        output_path = save_results(config_name, test_result, analysis, OUTPUT_DIR)
        print(f"Results saved to: {output_path}")

        results.append({
            "config_name": config_name,
            "status": test_result["status"],
            "analysis": analysis,
            "output_path": str(output_path),
        })

    # Print summary
    print_summary(results)

    return 0


if __name__ == "__main__":
    sys.exit(main())
