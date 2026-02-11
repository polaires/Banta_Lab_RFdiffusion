#!/usr/bin/env python3
"""
CLI for Unified Design Analyzer

Usage:
    python -m analyze_design_cli output.pdb --params params.json --metal TB
    python -m analyze_design_cli output.pdb --ligand mol.sdf --session my_session
"""
import argparse
import json
import os
import sys

from unified_analyzer import UnifiedDesignAnalyzer
from design_history import DesignHistoryManager
from filter_evaluator import FilterEvaluator
from lesson_detector import LessonDetector


def main():
    parser = argparse.ArgumentParser(
        description="Analyze protein design and save to history"
    )
    parser.add_argument("pdb_file", help="Path to PDB file to analyze")
    parser.add_argument("--params", help="Path to design parameters JSON")
    parser.add_argument("--ligand", help="Path to ligand SDF file")
    parser.add_argument("--metal", help="Metal type code (e.g., TB, ZN). Auto-detected if not specified.")
    parser.add_argument("--metal-chain", default=None, help="Metal chain ID (auto-detected if not specified)")
    parser.add_argument("--metal-resnum", type=int, default=None, help="Metal residue number (auto-detected if not specified)")
    parser.add_argument("--session", help="Session name (creates new if not exists)")
    parser.add_argument(
        "--history-dir",
        default="experiments/design_history",
        help="Path to design history directory"
    )
    parser.add_argument("--no-save", action="store_true", help="Don't save to history")
    parser.add_argument("--json", action="store_true", help="Output JSON only")

    args = parser.parse_args()

    # Read PDB file
    with open(args.pdb_file, "r") as f:
        pdb_content = f.read()

    # Read params if provided
    params = {}
    if args.params:
        with open(args.params, "r") as f:
            params = json.load(f)

    # Read ligand if provided
    ligand_sdf = None
    if args.ligand:
        with open(args.ligand, "r") as f:
            ligand_sdf = f.read()

    # Initialize analyzer
    analyzer = UnifiedDesignAnalyzer()

    # Run analysis
    metrics = analyzer.analyze(
        pdb_content=pdb_content,
        design_params=params,
        pdb_path=args.pdb_file,
        ligand_sdf=ligand_sdf,
        metal_type=args.metal,
        metal_chain=args.metal_chain,
        metal_resnum=args.metal_resnum,
    )

    # Evaluate against filter presets
    presets_dir = os.path.join(args.history_dir, "filter_presets")
    if os.path.exists(presets_dir):
        evaluator = FilterEvaluator(presets_dir)
        metrics["filter_results"] = evaluator.evaluate_all_presets(metrics)

    # Save to history if requested
    if not args.no_save:
        history = DesignHistoryManager(args.history_dir)

        # Get or create session
        session_name = args.session or "cli_session"
        session = history.start_session(session_name)

        # Save run
        run_id = history.save_run(
            session=session,
            params=params,
            outputs={"pdb": pdb_content},
            metrics=metrics,
        )

        # Check for lesson triggers
        index = history.load_index()
        detector = LessonDetector()

        # Build history for detector
        design_history = []
        for design in index["designs"][-20:]:
            metrics_path = os.path.join(
                args.history_dir, "runs", design["run_id"],
                "analysis", "metrics.json"
            )
            if os.path.exists(metrics_path):
                with open(metrics_path) as f:
                    design_history.append(json.load(f))

        trigger = detector.check_triggers(metrics, design_history)

        if not args.json:
            print(f"\nSaved to: {run_id}")
            if trigger:
                print(f"\n*** LESSON TRIGGER: {trigger.trigger_type} ***")
                print(f"    {trigger.description}")

    # Output
    if args.json:
        print(json.dumps(metrics, indent=2))
    else:
        print(f"\nDesign Analysis: {metrics['design_id']}")
        print(f"Type: {metrics['design_type']}")

        # Show auto-detected information
        if "auto_detected" in metrics:
            print("\nAuto-detected:")
            auto = metrics["auto_detected"]
            if "metal" in auto:
                m = auto["metal"]
                print(f"  Metal: {m['type']} (chain {m['chain']}, resnum {m['resnum']})")
            if "ligand" in auto:
                lig = auto["ligand"]
                print(f"  Ligand: {lig['name']} (chain {lig['chain']}, resnum {lig['resnum']}, {lig['atom_count']} atoms)")

        print("\nAnalyses:")
        for name, data in metrics["analyses"].items():
            status = data["status"]
            if status == "success":
                print(f"  {name}: {status}")
                for k, v in data.get("metrics", {}).items():
                    if not isinstance(v, (list, dict)):
                        print(f"    - {k}: {v}")
            else:
                print(f"  {name}: {status}")
                if data.get("reason"):
                    print(f"    ({data['reason']})")

        if "filter_results" in metrics:
            print("\nFilter Results:")
            for preset, result in metrics["filter_results"].items():
                status = "PASS" if result.get("pass") else "FAIL"
                print(f"  {preset}: {status}")


if __name__ == "__main__":
    main()
