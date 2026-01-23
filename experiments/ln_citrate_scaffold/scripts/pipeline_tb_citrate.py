"""
Automated TB-Citrate Binding Protein Design Pipeline

Workflow:
1. Parameter sweep (10 trials per config)
2. RF3 validation & filtering
3. Rank configs by success rate
4. Production run with best config (1000 designs)

Usage:
    python pipeline_tb_citrate.py --mode sweep      # Parameter optimization
    python pipeline_tb_citrate.py --mode production # Production with best config
"""

import argparse
import requests
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Any, Optional
from dataclasses import dataclass, asdict
import statistics

API_URL = "http://localhost:8000/runsync"

# ============================================================================
# CONFIGURATION
# ============================================================================

@dataclass
class DesignConfig:
    """Configuration for a design run."""
    name: str
    contig: str
    cfg_scale: float
    select_buried: Dict[str, str]
    select_hbond_acceptor: Dict[str, str]
    select_hbond_donor: Dict[str, str]
    temperature: float = 0.2
    num_seqs_per_backbone: int = 8

# Parameter sweep configurations
SWEEP_CONFIGS = [
    DesignConfig(
        name="small_low_cfg",
        contig="100-120",
        cfg_scale=1.5,
        select_buried={"X1": "all"},
        select_hbond_acceptor={"L1": "O1,O2,O3,O4,O5,O6"},
        select_hbond_donor={"L1": "O7"},
    ),
    DesignConfig(
        name="small_med_cfg",
        contig="100-120",
        cfg_scale=2.0,
        select_buried={"X1": "all"},
        select_hbond_acceptor={"L1": "O1,O2,O3,O4,O5,O6"},
        select_hbond_donor={"L1": "O7"},
    ),
    DesignConfig(
        name="small_high_cfg",
        contig="100-120",
        cfg_scale=2.5,
        select_buried={"X1": "all"},
        select_hbond_acceptor={"L1": "O1,O2,O3,O4,O5,O6"},
        select_hbond_donor={"L1": "O7"},
    ),
    DesignConfig(
        name="medium_low_cfg",
        contig="110-130",
        cfg_scale=1.5,
        select_buried={"X1": "all"},
        select_hbond_acceptor={"L1": "O1,O2,O3,O4,O5,O6"},
        select_hbond_donor={"L1": "O7"},
    ),
    DesignConfig(
        name="medium_med_cfg",
        contig="110-130",
        cfg_scale=2.0,
        select_buried={"X1": "all"},
        select_hbond_acceptor={"L1": "O1,O2,O3,O4,O5,O6"},
        select_hbond_donor={"L1": "O7"},
    ),
    DesignConfig(
        name="medium_high_cfg",
        contig="110-130",
        cfg_scale=2.5,
        select_buried={"X1": "all"},
        select_hbond_acceptor={"L1": "O1,O2,O3,O4,O5,O6"},
        select_hbond_donor={"L1": "O7"},
    ),
    DesignConfig(
        name="large_low_cfg",
        contig="130-150",
        cfg_scale=1.5,
        select_buried={"X1": "all"},
        select_hbond_acceptor={"L1": "O1,O2,O3,O4,O5,O6"},
        select_hbond_donor={"L1": "O7"},
    ),
    DesignConfig(
        name="large_med_cfg",
        contig="130-150",
        cfg_scale=2.0,
        select_buried={"X1": "all"},
        select_hbond_acceptor={"L1": "O1,O2,O3,O4,O5,O6"},
        select_hbond_donor={"L1": "O7"},
    ),
    DesignConfig(
        name="large_high_cfg",
        contig="130-150",
        cfg_scale=2.5,
        select_buried={"X1": "all"},
        select_hbond_acceptor={"L1": "O1,O2,O3,O4,O5,O6"},
        select_hbond_donor={"L1": "O7"},
    ),
]

# Filter thresholds
FILTER_STRICT = {"plddt": 0.80, "ptm": 0.80, "pae": 5.0}
FILTER_RELAXED = {"plddt": 0.75, "ptm": 0.70, "pae": 10.0}

# ============================================================================
# PIPELINE STAGES
# ============================================================================

def run_rfd3_scaffolding(pdb_content: str, config: DesignConfig, num_designs: int = 1) -> List[str]:
    """Stage 1: Generate scaffold backbones."""
    payload = {
        "input": {
            "task": "rfd3",
            "pdb_content": pdb_content,
            "contig": config.contig,
            "ligand": "CIT,TB",
            "num_designs": num_designs,
            "select_fixed_atoms": {"X1": "all", "L1": "all"},
            "select_buried": config.select_buried,
            "select_hbond_acceptor": config.select_hbond_acceptor,
            "select_hbond_donor": config.select_hbond_donor,
            "use_classifier_free_guidance": True,
            "cfg_scale": config.cfg_scale,
            "infer_ori_strategy": "com",
        }
    }

    response = requests.post(API_URL, json=payload, timeout=1200)
    result = response.json()

    if result.get("status", "").upper() == "COMPLETED":
        designs = result.get("output", {}).get("result", {}).get("designs", [])
        return [d.get("content", d) if isinstance(d, dict) else d for d in designs]
    return []


def run_mpnn_sequence_design(pdb_content: str, config: DesignConfig) -> List[str]:
    """Stage 2: Generate sequences for backbone."""
    payload = {
        "input": {
            "task": "mpnn",
            "pdb_content": pdb_content,
            "num_seqs": config.num_seqs_per_backbone,
            "ligand_mpnn_use_atom_context": 1,
            "temperature": config.temperature,
        }
    }

    response = requests.post(API_URL, json=payload, timeout=300)
    result = response.json()

    sequences = []
    if result.get("status", "").upper() == "COMPLETED":
        seq_data = result.get("output", {}).get("result", {}).get("sequences", [])
        for item in seq_data:
            if isinstance(item, dict) and "content" in item:
                # Parse FASTA content
                lines = item["content"].strip().split("\n")
                for i in range(0, len(lines), 2):
                    if i + 1 < len(lines) and lines[i].startswith(">"):
                        sequences.append(lines[i + 1])
            elif isinstance(item, str):
                sequences.append(item)
    return sequences


def run_rf3_validation(sequence: str, name: str) -> Dict[str, float]:
    """Stage 3: Validate sequence with RF3."""
    payload = {
        "input": {
            "task": "rf3",
            "sequence": sequence,
            "name": name
        }
    }

    response = requests.post(API_URL, json=payload, timeout=300)
    result = response.json()

    if result.get("status", "").upper() == "COMPLETED":
        confidences = result.get("output", {}).get("result", {}).get("confidences", {})
        return {
            "plddt": confidences.get("mean_plddt", 0),
            "ptm": confidences.get("ptm", 0),
            "pae": confidences.get("overall_pae", 999),
        }
    return {"plddt": 0, "ptm": 0, "pae": 999}


def analyze_sequence(sequence: str) -> Dict[str, Any]:
    """Stage 4: Analyze sequence composition."""
    length = len(sequence)
    e_count = sequence.count('E')
    d_count = sequence.count('D')
    a_count = sequence.count('A')

    carbox_pct = 100 * (e_count + d_count) / length
    ala_pct = 100 * a_count / length

    # Binding potential score
    score = 0
    if carbox_pct >= 15:
        score += 40
    elif carbox_pct >= 10:
        score += 25
    if sequence.count('H') >= 2:
        score += 20
    if (sequence.count('S') + sequence.count('T')) / length >= 0.10:
        score += 20
    if ala_pct < 25:
        score += 20

    return {
        "length": length,
        "carboxylate_pct": round(carbox_pct, 1),
        "alanine_pct": round(ala_pct, 1),
        "binding_score": score,
        "ala_flag": ala_pct > 25,
    }


def passes_filter(metrics: Dict, thresholds: Dict) -> bool:
    """Check if design passes filter thresholds."""
    return (
        metrics.get("plddt", 0) >= thresholds["plddt"] and
        metrics.get("ptm", 0) >= thresholds["ptm"] and
        metrics.get("pae", 999) <= thresholds["pae"]
    )


# ============================================================================
# PIPELINE ORCHESTRATION
# ============================================================================

def run_parameter_sweep(
    input_pdb: str,
    output_dir: Path,
    trials_per_config: int = 10,
    seqs_per_backbone: int = 4,
) -> Dict[str, Any]:
    """
    Run parameter sweep to find best configuration.

    For each config:
    1. Generate `trials_per_config` backbones
    2. Generate `seqs_per_backbone` sequences per backbone
    3. Validate all with RF3
    4. Calculate success rate and average metrics
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    results = {}

    for config in SWEEP_CONFIGS:
        print(f"\n{'='*60}")
        print(f"Testing: {config.name}")
        print(f"  Contig: {config.contig}, CFG: {config.cfg_scale}")
        print(f"{'='*60}")

        config.num_seqs_per_backbone = seqs_per_backbone

        all_metrics = []
        passed_strict = 0
        passed_relaxed = 0

        # Generate backbones
        backbones = run_rfd3_scaffolding(input_pdb, config, num_designs=trials_per_config)
        print(f"  Generated {len(backbones)} backbones")

        for i, backbone in enumerate(backbones):
            # Generate sequences
            sequences = run_mpnn_sequence_design(backbone, config)
            print(f"    Backbone {i}: {len(sequences)} sequences")

            for j, seq in enumerate(sequences):
                name = f"{config.name}_{i}_{j}"

                # RF3 validation
                rf3_metrics = run_rf3_validation(seq, name)
                seq_analysis = analyze_sequence(seq)

                metrics = {**rf3_metrics, **seq_analysis, "name": name, "sequence": seq}
                all_metrics.append(metrics)

                if passes_filter(rf3_metrics, FILTER_STRICT):
                    passed_strict += 1
                if passes_filter(rf3_metrics, FILTER_RELAXED):
                    passed_relaxed += 1

        # Calculate summary statistics
        if all_metrics:
            plddts = [m["plddt"] for m in all_metrics if m["plddt"] > 0]
            ptms = [m["ptm"] for m in all_metrics if m["ptm"] > 0]

            results[config.name] = {
                "config": asdict(config),
                "total_designs": len(all_metrics),
                "passed_strict": passed_strict,
                "passed_relaxed": passed_relaxed,
                "strict_rate": passed_strict / len(all_metrics) if all_metrics else 0,
                "relaxed_rate": passed_relaxed / len(all_metrics) if all_metrics else 0,
                "avg_plddt": statistics.mean(plddts) if plddts else 0,
                "max_plddt": max(plddts) if plddts else 0,
                "avg_ptm": statistics.mean(ptms) if ptms else 0,
                "max_ptm": max(ptms) if ptms else 0,
                "all_metrics": all_metrics,
            }

            print(f"  Results: {passed_strict}/{len(all_metrics)} strict, {passed_relaxed}/{len(all_metrics)} relaxed")
            print(f"  Avg pLDDT: {results[config.name]['avg_plddt']:.3f}, Max: {results[config.name]['max_plddt']:.3f}")

    # Save results
    results_file = output_dir / "sweep_results.json"
    with open(results_file, 'w') as f:
        # Remove sequences from saved results to reduce size
        save_results = {}
        for name, data in results.items():
            save_data = {k: v for k, v in data.items() if k != "all_metrics"}
            save_data["best_designs"] = sorted(
                data["all_metrics"],
                key=lambda x: x["plddt"],
                reverse=True
            )[:5]
            save_results[name] = save_data
        json.dump(save_results, f, indent=2)

    return results


def find_best_config(results: Dict) -> str:
    """Find best configuration based on sweep results."""
    # Rank by: strict_rate * avg_plddt * avg_ptm
    scores = {}
    for name, data in results.items():
        score = (
            data["strict_rate"] * 0.4 +
            data["avg_plddt"] * 0.3 +
            data["avg_ptm"] * 0.3
        )
        scores[name] = score

    best = max(scores, key=scores.get)
    return best


def run_production(
    input_pdb: str,
    output_dir: Path,
    config: DesignConfig,
    num_designs: int = 1000,
    batch_size: int = 10,
) -> List[Dict]:
    """
    Production run with best configuration.

    Generates designs in batches, validates, and saves passing designs.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    all_passing = []
    total_generated = 0

    print(f"\n{'='*60}")
    print(f"PRODUCTION RUN: {config.name}")
    print(f"Target: {num_designs} designs")
    print(f"{'='*60}")

    while total_generated < num_designs:
        batch_num = total_generated // batch_size

        # Generate batch of backbones
        backbones = run_rfd3_scaffolding(input_pdb, config, num_designs=batch_size)

        for i, backbone in enumerate(backbones):
            sequences = run_mpnn_sequence_design(backbone, config)

            for j, seq in enumerate(sequences):
                name = f"prod_{batch_num}_{i}_{j}"
                total_generated += 1

                # Validate
                rf3_metrics = run_rf3_validation(seq, name)
                seq_analysis = analyze_sequence(seq)

                if passes_filter(rf3_metrics, FILTER_STRICT):
                    design = {
                        "name": name,
                        "sequence": seq,
                        **rf3_metrics,
                        **seq_analysis,
                    }
                    all_passing.append(design)

                    # Save passing design
                    fasta_file = output_dir / f"{name}.fasta"
                    with open(fasta_file, 'w') as f:
                        f.write(f">{name}_pLDDT{rf3_metrics['plddt']:.2f}\n{seq}\n")

                if total_generated % 100 == 0:
                    print(f"  Progress: {total_generated}/{num_designs}, Passing: {len(all_passing)}")

    # Save summary
    summary = {
        "config": asdict(config),
        "total_generated": total_generated,
        "total_passing": len(all_passing),
        "pass_rate": len(all_passing) / total_generated if total_generated > 0 else 0,
        "best_designs": sorted(all_passing, key=lambda x: x["plddt"], reverse=True)[:20],
    }

    with open(output_dir / "production_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)

    # Save all passing sequences
    with open(output_dir / "all_passing.fasta", 'w') as f:
        for d in sorted(all_passing, key=lambda x: x["plddt"], reverse=True):
            f.write(f">{d['name']}_pLDDT{d['plddt']:.2f}_pTM{d['ptm']:.2f}\n")
            f.write(f"{d['sequence']}\n")

    print(f"\n{'='*60}")
    print(f"PRODUCTION COMPLETE")
    print(f"  Total generated: {total_generated}")
    print(f"  Total passing: {len(all_passing)} ({100*len(all_passing)/total_generated:.1f}%)")
    print(f"{'='*60}")

    return all_passing


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="TB-Citrate Design Pipeline")
    parser.add_argument("--mode", choices=["sweep", "production"], required=True)
    parser.add_argument("--input", type=str, default=None, help="Input PDB path")
    parser.add_argument("--output", type=str, default=None, help="Output directory")
    parser.add_argument("--trials", type=int, default=10, help="Trials per config (sweep mode)")
    parser.add_argument("--num-designs", type=int, default=1000, help="Number of designs (production mode)")
    parser.add_argument("--config", type=str, default=None, help="Config name for production")
    args = parser.parse_args()

    # Defaults
    base_dir = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\experiments\ln_citrate_scaffold")
    input_pdb = args.input or str(base_dir / "inputs" / "tb_citrate_motif_scaffold.pdb")
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    with open(input_pdb, 'r') as f:
        pdb_content = f.read()

    if args.mode == "sweep":
        output_dir = Path(args.output) if args.output else base_dir / "outputs" / f"sweep_{timestamp}"
        results = run_parameter_sweep(pdb_content, output_dir, trials_per_config=args.trials)

        # Find and report best config
        best = find_best_config(results)
        print(f"\n{'='*60}")
        print(f"BEST CONFIGURATION: {best}")
        print(f"  Strict pass rate: {results[best]['strict_rate']*100:.1f}%")
        print(f"  Avg pLDDT: {results[best]['avg_plddt']:.3f}")
        print(f"  Max pLDDT: {results[best]['max_plddt']:.3f}")
        print(f"{'='*60}")

    elif args.mode == "production":
        # Find config
        config_name = args.config or "large_high_cfg"  # Default to best from our experiments
        config = next((c for c in SWEEP_CONFIGS if c.name == config_name), SWEEP_CONFIGS[-1])

        output_dir = Path(args.output) if args.output else base_dir / "outputs" / f"production_{timestamp}"
        run_production(pdb_content, output_dir, config, num_designs=args.num_designs)


if __name__ == "__main__":
    main()
