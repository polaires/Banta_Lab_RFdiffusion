#!/usr/bin/env python3
"""
Run a single round of Ln-citrate scaffold design experiments.

Usage:
    python run_round.py --round 1 --n-designs 8

Requires Docker running locally (use docker-wsl skill).
"""

import argparse
import json
import os
import sys
from datetime import datetime
from pathlib import Path

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "backend" / "serverless"))

import requests

API_URL = "http://localhost:8000/runsync"


def load_configs(round_num: int) -> dict:
    """Load configs for a specific round."""
    config_path = Path(__file__).parent.parent / "configs" / f"round_{round_num:02d}" / f"r{round_num}_all_configs.json"

    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path) as f:
        configs = json.load(f)

    # Remove metadata key
    configs.pop("_meta", None)

    return configs


def run_rfd3_design(config: dict, n_designs: int = 8) -> dict:
    """
    Run RFD3 design via local Docker API.

    Args:
        config: RFD3 configuration dict
        n_designs: Number of designs to generate

    Returns:
        API response dict with designs
    """
    # Remove description field (not part of RFD3 spec)
    config = {k: v for k, v in config.items() if not k.startswith("_")}

    payload = {
        "input": {
            "task": "rfd3",
            "config": config,
            "n_designs": n_designs
        }
    }

    try:
        response = requests.post(API_URL, json=payload, timeout=600)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        return {"error": str(e), "status": "failed"}


def save_designs(results: dict, config_name: str, round_num: int, output_dir: Path):
    """Save generated PDB designs to output directory."""
    if results.get("status") == "failed" or "error" in results:
        print(f"  [ERROR] {config_name}: {results.get('error', 'Unknown error')}")
        return []

    saved_files = []
    designs = results.get("output", {}).get("designs", [])

    if not designs:
        designs = results.get("designs", [])

    for i, design in enumerate(designs):
        pdb_content = design.get("pdb_content", design.get("pdb", ""))
        if not pdb_content:
            continue

        filename = f"r{round_num}_{config_name}_{i:03d}.pdb"
        filepath = output_dir / filename

        with open(filepath, "w") as f:
            f.write(pdb_content)

        saved_files.append(filepath)
        print(f"  Saved: {filename}")

    return saved_files


def main():
    parser = argparse.ArgumentParser(description="Run Ln-citrate scaffold design round")
    parser.add_argument("--round", "-r", type=int, required=True, help="Round number")
    parser.add_argument("--n-designs", "-n", type=int, default=8, help="Designs per config")
    parser.add_argument("--config", "-c", type=str, default=None, help="Run specific config only")
    parser.add_argument("--dry-run", action="store_true", help="Print configs without running")
    args = parser.parse_args()

    # Load configs
    print(f"\n{'='*60}")
    print(f"Ln-Citrate Scaffold Design - Round {args.round}")
    print(f"{'='*60}\n")

    try:
        configs = load_configs(args.round)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)

    if args.config:
        if args.config not in configs:
            print(f"Error: Config '{args.config}' not found. Available: {list(configs.keys())}")
            sys.exit(1)
        configs = {args.config: configs[args.config]}

    print(f"Configs to run: {list(configs.keys())}")
    print(f"Designs per config: {args.n_designs}")
    print(f"Total expected: {len(configs) * args.n_designs}\n")

    if args.dry_run:
        print("DRY RUN - Configs that would be run:")
        for name, config in configs.items():
            print(f"\n{name}:")
            print(json.dumps(config, indent=2))
        return

    # Check API health
    print("Checking API health...")
    try:
        health = requests.post(API_URL, json={"input": {"task": "health"}}, timeout=30)
        if health.status_code != 200:
            print("Warning: API health check failed. Is Docker running?")
    except requests.exceptions.RequestException:
        print("Error: Cannot connect to API. Start Docker first (use docker-wsl skill).")
        sys.exit(1)

    # Setup output directory
    output_dir = Path(__file__).parent.parent / "outputs" / f"round_{args.round:02d}"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Run each config
    all_results = {}
    all_files = []

    for config_name, config in configs.items():
        print(f"\nRunning: {config_name}")
        print(f"  Description: {config.get('_description', 'N/A')}")

        results = run_rfd3_design(config, args.n_designs)
        saved = save_designs(results, config_name, args.round, output_dir)

        all_results[config_name] = {
            "n_generated": len(saved),
            "files": [str(f) for f in saved],
            "timestamp": datetime.now().isoformat()
        }
        all_files.extend(saved)

        print(f"  Generated: {len(saved)} designs")

    # Save run summary
    summary_path = output_dir / f"run_summary_r{args.round}.json"
    with open(summary_path, "w") as f:
        json.dump({
            "round": args.round,
            "timestamp": datetime.now().isoformat(),
            "n_designs_requested": args.n_designs,
            "total_generated": len(all_files),
            "configs": all_results
        }, f, indent=2)

    print(f"\n{'='*60}")
    print(f"Round {args.round} Complete")
    print(f"{'='*60}")
    print(f"Total designs generated: {len(all_files)}")
    print(f"Output directory: {output_dir}")
    print(f"Summary: {summary_path}")
    print(f"\nNext: Run analyze_round.py --round {args.round}")


if __name__ == "__main__":
    main()
