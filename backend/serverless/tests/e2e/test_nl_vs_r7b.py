"""
E2E Test: NL Pipeline vs Round 7b -Tb-Citrate Design Comparison

Runs both design approaches end-to-end against the local Docker API,
generating equivalent numbers of designs, and produces a structured
comparison of quality metrics.

Usage:
    cd backend/serverless
    python test_nl_vs_r7b.py
"""

import json
import os
import re
import statistics
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

import requests

# ============== Configuration ==============

API_URL = "http://localhost:8000/runsync"

TIMEOUTS = {
    "rfd3": 1800,
    "metal_binding_design": 1800,
    "mpnn": 600,
    "rf3": 600,
    "analyze_design": 120,
}

# Input PDB for R7b arm
R7B_INPUT_PDB = os.path.join(
    os.path.dirname(__file__),
    "..", "..", "experiments", "ln_citrate_scaffold", "inputs",
    "tb_citrate_motif_scaffold.pdb",
)

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "output_nl_vs_r7b")

# Citrate SMILES for NL pipeline RF3 calls
CITRATE_SMILES = "OC(=O)CC(O)(CC(O)=O)C(O)=O"


# ============== Data Classes ==============

@dataclass
class DesignResult:
    seq_id: str
    sequence: str
    backbone_id: str
    rf3_pdb: Optional[str] = None
    ptm: Optional[float] = None
    plddt: Optional[float] = None
    pae: Optional[float] = None
    coordination_number: Optional[float] = None
    geometry_rmsd: Optional[float] = None
    alanine_pct: Optional[float] = None
    total_residues: Optional[int] = None
    filter_passed: Optional[bool] = None
    failed_filters: Optional[List[str]] = None


@dataclass
class ArmResult:
    name: str
    backbones: List[Dict[str, Any]] = field(default_factory=list)
    sequences: List[str] = field(default_factory=list)
    designs: List[DesignResult] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)
    timings: Dict[str, float] = field(default_factory=dict)


# ============== API Helper ==============

def call_api(task: str, params: Dict[str, Any], timeout: Optional[int] = None,
             retries: int = 2) -> Dict[str, Any]:
    """Call the local Docker API with retry logic."""
    timeout = timeout or TIMEOUTS.get(task, 300)
    payload = {"input": {"task": task, **params}}

    for attempt in range(retries + 1):
        try:
            resp = requests.post(API_URL, json=payload, timeout=timeout)
            result = resp.json()

            # RunPod wraps in output
            if "output" in result:
                inner = result["output"]
                if isinstance(inner, dict):
                    return inner
            return result

        except requests.exceptions.Timeout:
            if attempt < retries:
                print(f"  [Retry {attempt+1}] Timeout on {task}, retrying...")
                time.sleep(5)
            else:
                return {"status": "failed", "error": f"Timeout after {timeout}s ({retries+1} attempts)"}
        except Exception as e:
            if attempt < retries:
                print(f"  [Retry {attempt+1}] Error on {task}: {e}")
                time.sleep(5)
            else:
                return {"status": "failed", "error": str(e)}

    return {"status": "failed", "error": "Exhausted retries"}


# ============== Sequence Parsing ==============

def parse_sequences(mpnn_result: Dict[str, Any]) -> List[str]:
    """Extract sequences from MPNN FASTA result."""
    sequences_data = mpnn_result.get("result", {}).get("sequences", [])
    if not sequences_data:
        return []

    fasta_content = ""
    for item in sequences_data:
        if isinstance(item, dict) and "content" in item:
            fasta_content += item["content"]
        elif isinstance(item, str):
            fasta_content += item

    seqs = []
    for line in fasta_content.strip().split("\n"):
        line = line.strip()
        if line and not line.startswith(">"):
            seqs.append(line)
    return seqs


# ============== RF3 Confidence Extraction ==============

def extract_rf3_metrics(rf3_result: Dict[str, Any]) -> Dict[str, Any]:
    """Extract pTM, pLDDT, PAE from RF3 result."""
    metrics = {"ptm": None, "plddt": None, "pae": None, "pdb_content": None}
    result = rf3_result.get("result", {})

    # Get PDB content
    predictions = result.get("predictions", [])
    if predictions:
        metrics["pdb_content"] = predictions[0].get("content")

    # Get confidences
    conf = result.get("confidences") or {}
    summary = conf.get("summary_confidences", conf)

    metrics["ptm"] = summary.get("ptm")
    metrics["plddt"] = summary.get("overall_plddt")
    metrics["pae"] = summary.get("overall_pae")

    return metrics


# ============== Analysis Extraction ==============

def extract_analysis_metrics(analysis_result: Dict[str, Any]) -> Dict[str, Any]:
    """Extract analysis metrics from the 'analyze' task result.

    The 'analyze' task (handle_analyze) returns binding site info, not
    UnifiedDesignAnalyzer output. We extract coordination numbers from
    binding_sites data.
    """
    out = {
        "coordination_number": None,
        "geometry_rmsd": None,
        "alanine_pct": None,
        "total_residues": None,
        "filter_passed": None,
        "failed_filters": [],
    }
    result = analysis_result.get("result", {})

    # Extract coordination from binding_sites (analyze task output)
    binding_sites = result.get("binding_sites", [])
    for site in binding_sites:
        # Find Tb site
        if site.get("ligand", {}).get("name") in ("TB",):
            residues = site.get("coordinating_residues", [])
            out["coordination_number"] = len(residues)
            break

    return out


def compute_alanine_pct(sequence: str) -> float:
    """Compute alanine percentage from a sequence string."""
    if not sequence:
        return 0.0
    return round(sequence.count("A") / len(sequence) * 100, 1)


# ============== RFD3 Functions ==============

def run_r7b_rfd3(pdb_content: str, contig: str, num_designs: int) -> List[Dict[str, Any]]:
    """Run RFD3 with Round 7b expert-tuned parameters."""
    params = {
        "pdb_content": pdb_content,
        "contig": contig,
        "ligand": "CIT,TB",
        "select_fixed_atoms": {"X1": "all", "L1": "all"},
        "select_buried": {"X1": "all"},
        "select_hbond_acceptor": {"L1": "O1,O2,O3,O4,O5,O6"},
        "select_hbond_donor": {"L1": "O7"},
        "use_classifier_free_guidance": True,
        "cfg_scale": 2.5,
        "infer_ori_strategy": "com",
        "num_designs": num_designs,
    }
    print(f"  [R7b RFD3] contig={contig}, num_designs={num_designs}")
    result = call_api("rfd3", params)

    if result.get("status") != "completed":
        print(f"  [R7b RFD3] FAILED: {result.get('error', 'unknown')}")
        return []

    designs = result.get("result", {}).get("designs", [])
    print(f"  [R7b RFD3] Got {len(designs)} backbones")
    return designs


def run_nl_rfd3(contig: str, num_designs: int) -> List[Dict[str, Any]]:
    """Run RFD3 via metal_binding_design single mode (NL pipeline path)."""
    params = {
        "mode": "single",
        "metal": "TB",
        "ligand": "CIT",
        "contig": contig,
        "cfg_scale": 2.0,
        "num_designs": num_designs,
        "bury_ligand": True,
    }
    print(f"  [NL RFD3] contig={contig}, num_designs={num_designs}")
    result = call_api("metal_binding_design", params)

    if result.get("status") != "completed":
        print(f"  [NL RFD3] FAILED: {result.get('error', 'unknown')}")
        return []

    designs = result.get("result", {}).get("designs", [])
    print(f"  [NL RFD3] Got {len(designs)} backbones")
    return designs


# ============== Arm Runner ==============

def run_arm(
    name: str,
    rfd3_fn: Callable,
    mpnn_params: Dict[str, Any],
    rf3_extra_params: Dict[str, Any],
) -> ArmResult:
    """Run a full design arm: RFD3 ->MPNN ->RF3 ->Analyze."""
    arm = ArmResult(name=name)
    t0 = time.time()

    # ---- Step 1: RFD3 (2 configs ->5 backbones) ----
    print(f"\n{'='*60}")
    print(f"  [{name}] Step 1: RFD3 backbone generation")
    print(f"{'='*60}")

    t_rfd3 = time.time()
    small_designs = rfd3_fn(contig="110-130", num_designs=3)
    medium_designs = rfd3_fn(contig="130-150", num_designs=2)

    all_backbones = []
    for i, d in enumerate(small_designs):
        content = d.get("content", d) if isinstance(d, dict) else d
        all_backbones.append({"id": f"small_{i:03d}", "content": content})
    for i, d in enumerate(medium_designs):
        content = d.get("content", d) if isinstance(d, dict) else d
        all_backbones.append({"id": f"medium_{i:03d}", "content": content})

    arm.backbones = all_backbones
    arm.timings["rfd3"] = time.time() - t_rfd3
    print(f"  [{name}] RFD3 complete: {len(all_backbones)} backbones "
          f"({arm.timings['rfd3']:.0f}s)")

    if not all_backbones:
        arm.errors.append("No backbones generated")
        return arm

    # ---- Step 2: MPNN (per backbone ->4 seqs each ->20 total) ----
    print(f"\n{'='*60}")
    print(f"  [{name}] Step 2: MPNN sequence design")
    print(f"{'='*60}")

    t_mpnn = time.time()
    all_seqs = []  # (backbone_id, sequence)

    for bb in all_backbones:
        params = {
            "pdb_content": bb["content"],
            "num_sequences": mpnn_params.get("num_seqs", 4),
            "temperature": mpnn_params.get("temperature", 0.1),
            "ligand_mpnn_use_atom_context": mpnn_params.get("ligand_mpnn_use_atom_context", 1),
        }
        # Optional params
        if mpnn_params.get("bias_AA"):
            params["bias_AA"] = mpnn_params["bias_AA"]
        if mpnn_params.get("omit_AA"):
            params["omit_AA"] = mpnn_params["omit_AA"]
        if mpnn_params.get("pack_side_chains"):
            params["pack_side_chains"] = True

        print(f"  [{name}] MPNN for backbone {bb['id']}...")
        result = call_api("mpnn", params)

        if result.get("status") != "completed":
            err = f"MPNN failed for {bb['id']}: {result.get('error', 'unknown')}"
            print(f"  [{name}] {err}")
            arm.errors.append(err)
            continue

        seqs = parse_sequences(result)
        print(f"  [{name}] Got {len(seqs)} sequences from {bb['id']}")
        for seq in seqs:
            all_seqs.append((bb["id"], seq))

    arm.sequences = [s for _, s in all_seqs]
    arm.timings["mpnn"] = time.time() - t_mpnn
    print(f"  [{name}] MPNN complete: {len(all_seqs)} sequences "
          f"({arm.timings['mpnn']:.0f}s)")

    if not all_seqs:
        arm.errors.append("No sequences designed")
        return arm

    # ---- Step 3: RF3 (per sequence ->predicted PDBs) ----
    print(f"\n{'='*60}")
    print(f"  [{name}] Step 3: RF3 structure prediction")
    print(f"{'='*60}")

    t_rf3 = time.time()
    rf3_results = []  # (backbone_id, seq, rf3_metrics)

    for i, (bb_id, seq) in enumerate(all_seqs):
        params = {"sequence": seq, "name": f"{name}_seq_{i+1:03d}"}
        params.update(rf3_extra_params)

        print(f"  [{name}] RF3 {i+1}/{len(all_seqs)} (backbone {bb_id})...")
        result = call_api("rf3", params)

        if result.get("status") != "completed":
            err = f"RF3 failed for seq {i+1}: {result.get('error', 'unknown')}"
            print(f"  [{name}] {err}")
            arm.errors.append(err)
            rf3_results.append((bb_id, seq, None))
            continue

        metrics = extract_rf3_metrics(result)
        print(f"  [{name}] RF3 seq {i+1}: pTM={metrics['ptm']}, "
              f"pLDDT={metrics['plddt']}, PAE={metrics['pae']}")
        rf3_results.append((bb_id, seq, metrics))

    arm.timings["rf3"] = time.time() - t_rf3
    completed_rf3 = sum(1 for _, _, m in rf3_results if m is not None)
    print(f"  [{name}] RF3 complete: {completed_rf3}/{len(rf3_results)} "
          f"({arm.timings['rf3']:.0f}s)")

    # ---- Step 4: Analyze (per RF3 output) ----
    print(f"\n{'='*60}")
    print(f"  [{name}] Step 4: Design analysis")
    print(f"{'='*60}")

    t_analyze = time.time()

    for i, (bb_id, seq, rf3_metrics) in enumerate(rf3_results):
        seq_id = f"seq_{i+1:03d}"
        design = DesignResult(
            seq_id=seq_id,
            sequence=seq,
            backbone_id=bb_id,
        )

        if rf3_metrics is None:
            arm.designs.append(design)
            continue

        design.rf3_pdb = rf3_metrics.get("pdb_content")
        design.ptm = rf3_metrics.get("ptm")
        design.plddt = rf3_metrics.get("plddt")
        design.pae = rf3_metrics.get("pae")

        # Compute sequence-derived metrics
        design.alanine_pct = compute_alanine_pct(seq)
        design.total_residues = len(seq)

        # Simple pass/fail based on pTM threshold
        if design.ptm is not None:
            design.filter_passed = design.ptm >= 0.6
            design.failed_filters = [] if design.filter_passed else ["ptm < 0.6"]

        # Run structural analysis for coordination info
        if design.rf3_pdb:
            print(f"  [{name}] Analyzing {seq_id}...")
            analysis = call_api("analyze", {
                "pdb_content": design.rf3_pdb,
                "target_ligands": ["TB"],
            })

            if analysis.get("status") == "completed":
                am = extract_analysis_metrics(analysis)
                design.coordination_number = am["coordination_number"]
            else:
                err = f"Analysis failed for {seq_id}: {analysis.get('error', 'unknown')}"
                print(f"  [{name}] {err}")
                arm.errors.append(err)

        arm.designs.append(design)

    arm.timings["analyze"] = time.time() - t_analyze
    arm.timings["total"] = time.time() - t0

    print(f"\n  [{name}] ARM COMPLETE -{len(arm.designs)} designs, "
          f"{arm.timings['total']:.0f}s total")

    return arm


# ============== Comparison ==============

def safe_mean(values: List[Optional[float]]) -> Optional[float]:
    """Mean of non-None values, or None if empty."""
    valid = [v for v in values if v is not None]
    return round(statistics.mean(valid), 4) if valid else None


def compare(r7b: ArmResult, nl: ArmResult) -> Dict[str, Any]:
    """Produce structured comparison of two arms."""
    def arm_stats(arm: ArmResult) -> Dict[str, Any]:
        designs = arm.designs
        rf3_done = [d for d in designs if d.ptm is not None]
        passed = [d for d in designs if d.filter_passed is True]

        return {
            "backbones": len(arm.backbones),
            "sequences": len(arm.sequences),
            "rf3_completed": len(rf3_done),
            "filter_passed": len(passed),
            "pass_rate": round(len(passed) / len(rf3_done) * 100, 1) if rf3_done else 0,
            "avg_ptm": safe_mean([d.ptm for d in designs]),
            "avg_plddt": safe_mean([d.plddt for d in designs]),
            "avg_pae": safe_mean([d.pae for d in designs]),
            "avg_coordination": safe_mean([d.coordination_number for d in designs]),
            "avg_geometry_rmsd": safe_mean([d.geometry_rmsd for d in designs]),
            "avg_alanine_pct": safe_mean([d.alanine_pct for d in designs]),
            "best_ptm_design": None,
            "best_ptm_value": None,
            "timings": arm.timings,
            "errors": arm.errors,
        }

    r7b_stats = arm_stats(r7b)
    nl_stats = arm_stats(nl)

    # Find best pTM across both
    all_designs = [(r7b.name, d) for d in r7b.designs] + [(nl.name, d) for d in nl.designs]
    best = max(
        [(name, d) for name, d in all_designs if d.ptm is not None],
        key=lambda x: x[1].ptm,
        default=(None, None),
    )
    best_name, best_design = best

    # Update best in each arm
    for arm, stats in [(r7b, r7b_stats), (nl, nl_stats)]:
        arm_best = max(
            [d for d in arm.designs if d.ptm is not None],
            key=lambda d: d.ptm,
            default=None,
        )
        if arm_best:
            stats["best_ptm_design"] = arm_best.seq_id
            stats["best_ptm_value"] = arm_best.ptm

    return {
        "r7b": r7b_stats,
        "nl": nl_stats,
        "overall_best": {
            "arm": best_name,
            "design": best_design.seq_id if best_design else None,
            "ptm": best_design.ptm if best_design else None,
        },
        "parameter_differences": {
            "input_pdb": ["Hand-crafted motif + 5 motif residues", "Auto template (metal+ligand only)"],
            "cfg_scale": [2.5, 2.0],
            "burial": ["Metal only", "Metal + ligand"],
            "hbond_acceptor": ["O1-O6", "O1-O7 (all O as acceptors)"],
            "hbond_donor": ["O7", "None"],
            "mpnn_temp": [0.2, 0.1],
            "mpnn_bias": ["None", "A:-2.0, omit C"],
            "rf3_context": ["Seq only", "Seq + ligand + metal"],
        },
    }


# ============== Output ==============

def save_outputs(arm: ArmResult, base_dir: str):
    """Save arm outputs to disk."""
    arm_dir = os.path.join(base_dir, arm.name)
    os.makedirs(arm_dir, exist_ok=True)

    # Save backbones
    for bb in arm.backbones:
        path = os.path.join(arm_dir, f"backbone_{bb['id']}.pdb")
        with open(path, "w") as f:
            f.write(bb["content"])

    # Save sequences as FASTA
    fasta_path = os.path.join(arm_dir, "sequences.fasta")
    with open(fasta_path, "w") as f:
        for i, seq in enumerate(arm.sequences):
            f.write(f">seq_{i+1:03d}\n{seq}\n")

    # Save RF3 PDBs
    for d in arm.designs:
        if d.rf3_pdb:
            path = os.path.join(arm_dir, f"rf3_{d.seq_id}.pdb")
            with open(path, "w") as f:
                f.write(d.rf3_pdb)

    # Save arm summary
    summary = {
        "name": arm.name,
        "num_backbones": len(arm.backbones),
        "num_sequences": len(arm.sequences),
        "num_designs": len(arm.designs),
        "timings": arm.timings,
        "errors": arm.errors,
        "designs": [
            {
                "seq_id": d.seq_id,
                "backbone_id": d.backbone_id,
                "sequence": d.sequence,
                "ptm": d.ptm,
                "plddt": d.plddt,
                "pae": d.pae,
                "coordination_number": d.coordination_number,
                "geometry_rmsd": d.geometry_rmsd,
                "alanine_pct": d.alanine_pct,
                "total_residues": d.total_residues,
                "filter_passed": d.filter_passed,
                "failed_filters": d.failed_filters,
            }
            for d in arm.designs
        ],
    }
    with open(os.path.join(arm_dir, "arm_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)


def format_report(comparison: Dict[str, Any]) -> str:
    """Format a human-readable comparison report."""
    r = comparison["r7b"]
    n = comparison["nl"]
    best = comparison["overall_best"]

    def fmt(val, decimals=3):
        if val is None:
            return "N/A"
        if isinstance(val, float):
            return f"{val:.{decimals}f}"
        return str(val)

    lines = [
        "=" * 64,
        "        NL PIPELINE vs ROUND 7b -COMPARISON",
        "=" * 64,
        "",
        f"{'':24s} {'R7b':>12s}   {'NL Pipeline':>12s}",
        f"{'Backbones generated:':24s} {r['backbones']:>12d}   {n['backbones']:>12d}",
        f"{'Sequences designed:':24s} {r['sequences']:>12d}   {n['sequences']:>12d}",
        f"{'RF3 completed:':24s} {r['rf3_completed']:>12d}   {n['rf3_completed']:>12d}",
        f"{'Filter passed:':24s} {r['filter_passed']:>12d}   {n['filter_passed']:>12d}",
        f"{'Pass rate:':24s} {fmt(r['pass_rate'], 1) + '%':>12s}   {fmt(n['pass_rate'], 1) + '%':>12s}",
        "",
        f"{'Avg pTM:':24s} {fmt(r['avg_ptm']):>12s}   {fmt(n['avg_ptm']):>12s}",
        f"{'Avg pLDDT:':24s} {fmt(r['avg_plddt']):>12s}   {fmt(n['avg_plddt']):>12s}",
        f"{'Avg PAE:':24s} {fmt(r['avg_pae'], 2):>12s}   {fmt(n['avg_pae'], 2):>12s}",
        f"{'Avg coordination #:':24s} {fmt(r['avg_coordination'], 1):>12s}   {fmt(n['avg_coordination'], 1):>12s}",
        f"{'Avg geometry RMSD:':24s} {fmt(r['avg_geometry_rmsd'], 2):>12s}   {fmt(n['avg_geometry_rmsd'], 2):>12s}",
        f"{'Avg alanine %:':24s} {fmt(r['avg_alanine_pct'], 1) + '%' if r['avg_alanine_pct'] is not None else 'N/A':>12s}"
        f"   {fmt(n['avg_alanine_pct'], 1) + '%' if n['avg_alanine_pct'] is not None else 'N/A':>12s}",
        "",
        f"Best pTM (R7b):          {r['best_ptm_design'] or 'N/A'} ({fmt(r['best_ptm_value'])})",
        f"Best pTM (NL):           {n['best_ptm_design'] or 'N/A'} ({fmt(n['best_ptm_value'])})",
        f"Overall best:            {best['arm']} {best['design'] or 'N/A'} ({fmt(best['ptm'])})",
        "",
        "Timings (seconds):",
        f"  R7b:  RFD3={r['timings'].get('rfd3', 0):.0f}  MPNN={r['timings'].get('mpnn', 0):.0f}"
        f"  RF3={r['timings'].get('rf3', 0):.0f}  Analyze={r['timings'].get('analyze', 0):.0f}"
        f"  Total={r['timings'].get('total', 0):.0f}",
        f"  NL:   RFD3={n['timings'].get('rfd3', 0):.0f}  MPNN={n['timings'].get('mpnn', 0):.0f}"
        f"  RF3={n['timings'].get('rf3', 0):.0f}  Analyze={n['timings'].get('analyze', 0):.0f}"
        f"  Total={n['timings'].get('total', 0):.0f}",
        "",
        "Parameter Differences:",
    ]

    for key, (r_val, n_val) in comparison["parameter_differences"].items():
        lines.append(f"  {key:20s}  R7b: {r_val}")
        lines.append(f"  {'':20s}  NL:  {n_val}")

    if r["errors"] or n["errors"]:
        lines.append("")
        lines.append("Errors:")
        for e in r["errors"]:
            lines.append(f"  [R7b] {e}")
        for e in n["errors"]:
            lines.append(f"  [NL]  {e}")

    lines.append("=" * 64)
    return "\n".join(lines)


# ============== Main ==============

def main():
    print("=" * 64)
    print("  NL PIPELINE vs ROUND 7b -Tb-Citrate E2E Comparison")
    print("=" * 64)

    # Health check
    print("\nChecking API health...")
    health = call_api("health", {}, timeout=30)
    if health.get("status") not in ("completed", "ok"):
        print(f"API health check failed: {health}")
        print("Make sure Docker is running (use docker-wsl skill)")
        sys.exit(1)
    print("API is healthy.\n")

    # Load R7b input PDB
    r7b_pdb_path = os.path.normpath(R7B_INPUT_PDB)
    if not os.path.exists(r7b_pdb_path):
        print(f"ERROR: R7b input PDB not found: {r7b_pdb_path}")
        sys.exit(1)

    with open(r7b_pdb_path, "r") as f:
        r7b_pdb_content = f.read()
    print(f"Loaded R7b input PDB: {len(r7b_pdb_content)} chars")

    # ---- Run R7b Arm ----
    print("\n" + "#" * 64)
    print("#  ARM 1: ROUND 7b (Expert-Tuned)")
    print("#" * 64)

    r7b_result = run_arm(
        name="r7b",
        rfd3_fn=lambda contig, num_designs: run_r7b_rfd3(r7b_pdb_content, contig, num_designs),
        mpnn_params={
            "num_seqs": 4,
            "temperature": 0.2,
            "ligand_mpnn_use_atom_context": 1,
        },
        rf3_extra_params={},  # Sequence only -no metal, no ligand
    )

    # ---- Run NL Arm ----
    print("\n" + "#" * 64)
    print("#  ARM 2: NL PIPELINE (Automated Defaults)")
    print("#" * 64)

    nl_result = run_arm(
        name="nl",
        rfd3_fn=lambda contig, num_designs: run_nl_rfd3(contig, num_designs),
        mpnn_params={
            "num_seqs": 4,
            "temperature": 0.1,
            "bias_AA": "A:-2.0",
            "omit_AA": "C",
            "pack_side_chains": True,
            "ligand_mpnn_use_atom_context": 1,
        },
        rf3_extra_params={
            "ligand_smiles": CITRATE_SMILES,
            "metal": "TB",
        },
    )

    # ---- Compare ----
    print("\n" + "#" * 64)
    print("#  COMPARISON")
    print("#" * 64)

    comparison = compare(r7b_result, nl_result)

    # ---- Save outputs ----
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    save_outputs(r7b_result, OUTPUT_DIR)
    save_outputs(nl_result, OUTPUT_DIR)

    # Save comparison JSON
    comp_path = os.path.join(OUTPUT_DIR, "comparison.json")
    with open(comp_path, "w") as f:
        json.dump(comparison, f, indent=2)
    print(f"\nSaved comparison.json ->{comp_path}")

    # Save report
    report = format_report(comparison)
    report_path = os.path.join(OUTPUT_DIR, "report.txt")
    with open(report_path, "w") as f:
        f.write(report)
    print(f"Saved report.txt ->{report_path}")

    # Print report
    print("\n" + report)


if __name__ == "__main__":
    main()
