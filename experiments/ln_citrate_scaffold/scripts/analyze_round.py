#!/usr/bin/env python3
"""
Analyze designs from a round of Ln-citrate scaffold experiments.

Usage:
    python analyze_round.py --round 1 --save

Computes:
- Coordination number and geometry
- Metal burial (SASA)
- Citrate contacts
- Secondary structure content
- Tier assignment (S/A/B/C/F)
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "backend" / "serverless"))

try:
    from unified_analyzer import UnifiedDesignAnalyzer
    ANALYZER_AVAILABLE = True
except ImportError:
    ANALYZER_AVAILABLE = False
    print("Warning: UnifiedDesignAnalyzer not available. Using fallback metrics.")


# Tier thresholds for Ln-citrate designs
TIER_THRESHOLDS = {
    "S": {"cn_min": 8, "geom_max": 0.8, "sasa_max": 2, "cit_min": 3},
    "A": {"cn_min": 8, "geom_max": 1.2, "sasa_max": 5, "cit_min": 3},
    "B": {"cn_min": 7, "geom_max": 1.5, "sasa_max": 10, "cit_min": 2},
    "C": {"cn_min": 6, "geom_max": 2.0, "sasa_max": 20, "cit_min": 1},
}


def assign_tier(cn: int, geom_rmsd: float, metal_sasa: float, cit_contacts: int) -> str:
    """Assign quality tier based on metrics."""
    for tier in ["S", "A", "B", "C"]:
        t = TIER_THRESHOLDS[tier]
        if (cn >= t["cn_min"] and
            geom_rmsd <= t["geom_max"] and
            metal_sasa <= t["sasa_max"] and
            cit_contacts >= t["cit_min"]):
            return tier
    return "F"


def compute_score(cn: int, geom_rmsd: float, metal_sasa: float, cit_contacts: int, ss_content: float = 0.5) -> float:
    """Compute composite score (0-100)."""
    score = 0.0
    score += min(30, cn * 3.75)                     # CN: max 30
    score += max(0, 25 - geom_rmsd * 16.67)         # Geometry: max 25
    score += max(0, 20 - metal_sasa * 1.0)          # Burial: max 20
    score += min(15, cit_contacts * 5)              # Citrate: max 15
    score += min(10, ss_content * 10)               # SS: max 10
    return round(score, 2)


def analyze_pdb_fallback(pdb_path: Path) -> Dict[str, Any]:
    """
    Fallback analysis when UnifiedDesignAnalyzer is not available.
    Returns placeholder metrics.
    """
    return {
        "file": pdb_path.name,
        "coordination_number": 0,
        "geometry_rmsd": 999.0,
        "metal_sasa": 999.0,
        "citrate_contacts": 0,
        "ss_content": 0.5,
        "tier": "F",
        "score": 0.0,
        "error": "UnifiedDesignAnalyzer not available"
    }


def analyze_pdb(pdb_path: Path, analyzer=None) -> Dict[str, Any]:
    """
    Analyze a single PDB file.

    Returns:
        Dict with metrics, tier, and score
    """
    if analyzer is None:
        return analyze_pdb_fallback(pdb_path)

    try:
        result = analyzer.analyze(str(pdb_path))

        # Extract metrics
        cn = result.get("coordination", {}).get("coordination_number", 0)
        geom_rmsd = result.get("coordination", {}).get("geometry_rmsd", 999.0)
        metal_sasa = result.get("burial", {}).get("metal_sasa", 999.0)
        cit_contacts = result.get("ligand_interface", {}).get("contacts_to_metal", 0)
        ss = result.get("secondary_structure", {})
        ss_content = ss.get("helix", 0) + ss.get("sheet", 0)

        # Assign tier and score
        tier = assign_tier(cn, geom_rmsd, metal_sasa, cit_contacts)
        score = compute_score(cn, geom_rmsd, metal_sasa, cit_contacts, ss_content)

        return {
            "file": pdb_path.name,
            "coordination_number": cn,
            "geometry_rmsd": round(geom_rmsd, 3),
            "metal_sasa": round(metal_sasa, 2),
            "citrate_contacts": cit_contacts,
            "ss_content": round(ss_content, 3),
            "tier": tier,
            "score": score,
            "raw_result": result
        }

    except Exception as e:
        return {
            "file": pdb_path.name,
            "coordination_number": 0,
            "geometry_rmsd": 999.0,
            "metal_sasa": 999.0,
            "citrate_contacts": 0,
            "ss_content": 0.0,
            "tier": "F",
            "score": 0.0,
            "error": str(e)
        }


def compute_statistics(results: List[Dict]) -> Dict[str, Any]:
    """Compute aggregate statistics for a set of results."""
    if not results:
        return {}

    # Filter out errors
    valid = [r for r in results if "error" not in r or r.get("score", 0) > 0]

    if not valid:
        return {"error": "No valid results"}

    # Tier distribution
    tier_counts = {"S": 0, "A": 0, "B": 0, "C": 0, "F": 0}
    for r in valid:
        tier_counts[r.get("tier", "F")] += 1

    # Metrics
    cns = [r["coordination_number"] for r in valid]
    geoms = [r["geometry_rmsd"] for r in valid if r["geometry_rmsd"] < 900]
    sasas = [r["metal_sasa"] for r in valid if r["metal_sasa"] < 900]
    scores = [r["score"] for r in valid]

    # Pass rate (Tier B or better)
    n_pass = tier_counts["S"] + tier_counts["A"] + tier_counts["B"]
    pass_rate = n_pass / len(valid) if valid else 0

    return {
        "n_total": len(results),
        "n_valid": len(valid),
        "n_errors": len(results) - len(valid),
        "tier_distribution": tier_counts,
        "pass_rate": round(pass_rate, 3),
        "metrics": {
            "cn": {
                "mean": round(np.mean(cns), 2) if cns else 0,
                "std": round(np.std(cns), 2) if cns else 0,
                "min": min(cns) if cns else 0,
                "max": max(cns) if cns else 0
            },
            "geometry_rmsd": {
                "mean": round(np.mean(geoms), 3) if geoms else 0,
                "std": round(np.std(geoms), 3) if geoms else 0,
                "min": round(min(geoms), 3) if geoms else 0,
                "max": round(max(geoms), 3) if geoms else 0
            },
            "metal_sasa": {
                "mean": round(np.mean(sasas), 2) if sasas else 0,
                "std": round(np.std(sasas), 2) if sasas else 0,
                "min": round(min(sasas), 2) if sasas else 0,
                "max": round(max(sasas), 2) if sasas else 0
            },
            "score": {
                "mean": round(np.mean(scores), 2) if scores else 0,
                "std": round(np.std(scores), 2) if scores else 0,
                "min": round(min(scores), 2) if scores else 0,
                "max": round(max(scores), 2) if scores else 0
            }
        },
        "best_design": max(valid, key=lambda x: x["score"])["file"] if valid else None,
        "best_score": max(scores) if scores else 0
    }


def generate_judgment(stats_by_config: Dict[str, Dict], round_num: int) -> Dict[str, Any]:
    """
    Generate judgment and recommendations based on round results.
    """
    # Find best performing config
    best_config = None
    best_pass_rate = 0

    for config_name, stats in stats_by_config.items():
        if stats.get("pass_rate", 0) > best_pass_rate:
            best_pass_rate = stats["pass_rate"]
            best_config = config_name

    # Overall pass rate
    all_results = []
    for stats in stats_by_config.values():
        tier_dist = stats.get("tier_distribution", {})
        for tier, count in tier_dist.items():
            all_results.extend([tier] * count)

    overall_pass = sum(1 for t in all_results if t in ["S", "A", "B"]) / len(all_results) if all_results else 0

    # Determine action
    if overall_pass < 0.10:
        action = "REVISE_FUNDAMENTALLY"
        recommendations = [
            "Re-examine input structure preparation",
            "Check if Ln-citrate geometry is physically reasonable",
            "Consider alternative coordinating residue arrangements",
            "Verify metal and citrate are correctly positioned in input"
        ]
    elif overall_pass < 0.30:
        action = "ADJUST_KEY_PARAMETERS"
        recommendations = [
            f"Focus on best config: {best_config}",
            "Vary partial_t (try ±3Å from current)",
            "Adjust select_fixed_atoms (try BKBN vs TIP)",
            "Try different contig lengths"
        ]
    elif overall_pass < 0.50:
        action = "FINE_TUNE"
        recommendations = [
            f"Continue with {best_config} strategy",
            "Add H-bond conditioning (select_hbond_acceptor for citrate)",
            "Refine RASA selections",
            "Increase n_designs for better statistics"
        ]
    else:
        action = "PROCEED_TO_MPNN"
        recommendations = [
            "Backbone quality sufficient for sequence design",
            "Select top 5-10 backbones for LigandMPNN",
            "Use bias_AA: 'E:3.0,D:3.0,N:2.4,Q:2.4,C:-5.0,A:-2.0'",
            "Fix coordinating residue positions"
        ]

    return {
        "round": round_num,
        "overall_pass_rate": round(overall_pass, 3),
        "best_config": best_config,
        "best_config_pass_rate": round(best_pass_rate, 3),
        "action": action,
        "recommendations": recommendations
    }


def main():
    parser = argparse.ArgumentParser(description="Analyze Ln-citrate scaffold design round")
    parser.add_argument("--round", "-r", type=int, required=True, help="Round number")
    parser.add_argument("--save", "-s", action="store_true", help="Save analysis to file")
    parser.add_argument("--verbose", "-v", action="store_true", help="Print detailed results")
    args = parser.parse_args()

    print(f"\n{'='*60}")
    print(f"Ln-Citrate Scaffold Analysis - Round {args.round}")
    print(f"{'='*60}\n")

    # Find output directory
    output_dir = Path(__file__).parent.parent / "outputs" / f"round_{args.round:02d}"
    if not output_dir.exists():
        print(f"Error: Output directory not found: {output_dir}")
        print("Run run_round.py first.")
        sys.exit(1)

    # Find all PDB files
    pdb_files = list(output_dir.glob("*.pdb"))
    if not pdb_files:
        print(f"Error: No PDB files found in {output_dir}")
        sys.exit(1)

    print(f"Found {len(pdb_files)} PDB files to analyze\n")

    # Initialize analyzer
    analyzer = None
    if ANALYZER_AVAILABLE:
        try:
            analyzer = UnifiedDesignAnalyzer(metal="TB")
            print("Using UnifiedDesignAnalyzer")
        except Exception as e:
            print(f"Warning: Could not initialize analyzer: {e}")
    else:
        print("Warning: Using fallback analysis (limited metrics)")

    # Analyze all files
    all_results = []
    results_by_config = {}

    for pdb_file in sorted(pdb_files):
        result = analyze_pdb(pdb_file, analyzer)
        all_results.append(result)

        # Group by config
        # Expected format: r1_configname_001.pdb
        parts = pdb_file.stem.split("_")
        if len(parts) >= 2:
            config_name = "_".join(parts[:-1])  # Everything except index
            if config_name not in results_by_config:
                results_by_config[config_name] = []
            results_by_config[config_name].append(result)

        if args.verbose:
            print(f"  {result['file']}: Tier={result['tier']}, Score={result['score']:.1f}, CN={result['coordination_number']}")

    # Compute statistics by config
    stats_by_config = {}
    for config_name, results in results_by_config.items():
        stats_by_config[config_name] = compute_statistics(results)

    # Compute overall statistics
    overall_stats = compute_statistics(all_results)

    # Generate judgment
    judgment = generate_judgment(stats_by_config, args.round)

    # Print summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}\n")

    print("Overall Statistics:")
    print(f"  Total designs: {overall_stats.get('n_total', 0)}")
    print(f"  Valid designs: {overall_stats.get('n_valid', 0)}")
    print(f"  Pass rate: {overall_stats.get('pass_rate', 0)*100:.1f}%")
    print(f"  Best score: {overall_stats.get('best_score', 0):.1f}")
    print(f"  Best design: {overall_stats.get('best_design', 'N/A')}")

    print(f"\nTier Distribution:")
    tier_dist = overall_stats.get("tier_distribution", {})
    for tier in ["S", "A", "B", "C", "F"]:
        count = tier_dist.get(tier, 0)
        pct = count / overall_stats.get("n_valid", 1) * 100
        bar = "#" * int(pct / 5)
        print(f"  {tier}: {count:3d} ({pct:5.1f}%) {bar}")

    print(f"\nBy Config:")
    for config_name, stats in stats_by_config.items():
        print(f"  {config_name}:")
        print(f"    Pass rate: {stats.get('pass_rate', 0)*100:.1f}%")
        print(f"    Mean score: {stats.get('metrics', {}).get('score', {}).get('mean', 0):.1f}")
        print(f"    Mean CN: {stats.get('metrics', {}).get('cn', {}).get('mean', 0):.1f}")

    print(f"\n{'='*60}")
    print("JUDGMENT")
    print(f"{'='*60}\n")

    print(f"Action: {judgment['action']}")
    print(f"Best config: {judgment['best_config']} ({judgment['best_config_pass_rate']*100:.1f}% pass)")
    print(f"\nRecommendations:")
    for i, rec in enumerate(judgment['recommendations'], 1):
        print(f"  {i}. {rec}")

    # Save results
    if args.save:
        analysis_dir = Path(__file__).parent.parent / "analysis"
        analysis_dir.mkdir(exist_ok=True)

        analysis_file = analysis_dir / f"round_{args.round:02d}_analysis.json"
        with open(analysis_file, "w") as f:
            json.dump({
                "round": args.round,
                "timestamp": datetime.now().isoformat(),
                "overall_stats": overall_stats,
                "stats_by_config": stats_by_config,
                "judgment": judgment,
                "individual_results": all_results
            }, f, indent=2, default=str)

        print(f"\nAnalysis saved to: {analysis_file}")

        # Update cumulative lessons
        lessons_file = analysis_dir / "cumulative_lessons.md"
        with open(lessons_file, "a") as f:
            f.write(f"\n\n## Round {args.round} ({datetime.now().strftime('%Y-%m-%d')})\n\n")
            f.write(f"- Pass rate: {overall_stats.get('pass_rate', 0)*100:.1f}%\n")
            f.write(f"- Best config: {judgment['best_config']}\n")
            f.write(f"- Action: {judgment['action']}\n")
            f.write(f"- Recommendations:\n")
            for rec in judgment['recommendations']:
                f.write(f"  - {rec}\n")

        print(f"Lessons updated: {lessons_file}")

    print(f"\nNext steps based on judgment:")
    if judgment['action'] == "PROCEED_TO_MPNN":
        print("  1. Select top backbones from outputs/")
        print("  2. Run LigandMPNN with lanthanide bias")
        print("  3. Validate with AF2/PyRosetta")
    else:
        print(f"  1. Review recommendations above")
        print(f"  2. Update configs/round_{args.round+1:02d}/")
        print(f"  3. Run: python run_round.py --round {args.round+1}")


if __name__ == "__main__":
    main()
