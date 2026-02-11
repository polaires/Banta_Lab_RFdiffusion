#!/usr/bin/env python3
"""
E2E Test: H-bond Network Scaffolding for 3C9H Tb-Citrate Design

Tests the H-bond-only motif (NO metal coordinators) approach where:
1. Citrate H-bond network is FIXED (for selectivity): 69-70, 214-216, 238
2. Metal coordination is NOT fixed - RFD3 designs for Tb's preferred CN=8-9
3. Uses unindex parameter for multi-island scaffolding

Usage:
    cd backend/serverless
    python test_hbond_scaffolding.py --pilot     # Quick test (5 backbones)
    python test_hbond_scaffolding.py --batch 1   # Production batch (50 backbones)
"""

import argparse
import json
import os
import re
import statistics
import sys
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import requests

# ============== Configuration ==============

# API endpoint - use environment variable or default to localhost
# Set RUNPOD_API_URL for production, RUNPOD_API_KEY for auth
API_URL = os.environ.get("RUNPOD_API_URL", "http://localhost:8000/runsync")
API_KEY = os.environ.get("RUNPOD_API_KEY", "")

TIMEOUTS = {
    "rfd3": 3600,       # Scaffolding can be slower
    "mpnn": 600,
    "rf3": 600,
    "analyze_design": 300,
}

# H-bond motif PDB (created for US-001)
MOTIF_PDB_PATH = Path(__file__).parent.parent.parent / "experiments" / "ln_citrate_scaffold" / "scaffolding_3c9h" / "inputs" / "3c9h_hbond_motif.pdb"

OUTPUT_DIR = Path(__file__).parent.parent.parent / "experiments" / "ln_citrate_scaffold" / "scaffolding_3c9h" / "outputs"

# Citrate SMILES for RF3 ligand-aware prediction
CITRATE_SMILES = "OC(=O)CC(O)(CC(O)=O)C(O)=O"

# Best parameters from de novo exploration (progress.txt)
BEST_PARAMS = {
    "cfg_scale": 1.5,
    "contig": "100-130",
    "num_timesteps": 200,
    "step_scale": 1.0,
    "mpnn_temperature": 0.05,
    "mpnn_bias_AA": "A:-2.0",
    "mpnn_omit_AA": "C",
}


# ============== Data Classes ==============

@dataclass
class DesignResult:
    seq_id: str
    sequence: str
    backbone_id: str
    backbone_ptm: Optional[float] = None
    backbone_plddt: Optional[float] = None
    rf3_pdb: Optional[str] = None
    ptm: Optional[float] = None
    plddt: Optional[float] = None
    pae: Optional[float] = None
    coordination_number: Optional[int] = None
    geometry_rmsd: Optional[float] = None
    alanine_pct: Optional[float] = None
    total_residues: Optional[int] = None
    filter_passed: Optional[bool] = None
    failed_filters: Optional[List[str]] = None
    # H-bond specific metrics
    binding_site_plddt: Optional[float] = None
    citrate_hbond_count: Optional[int] = None


@dataclass
class BatchResult:
    name: str
    backbones: List[Dict[str, Any]] = field(default_factory=list)
    designs: List[DesignResult] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)
    timings: Dict[str, float] = field(default_factory=dict)


# ============== API Helper ==============

def call_api(task: str, params: Dict[str, Any], timeout: Optional[int] = None,
             retries: int = 2) -> Dict[str, Any]:
    """Call the API with retry logic. Supports local Docker or RunPod."""
    timeout = timeout or TIMEOUTS.get(task, 300)
    payload = {"input": {"task": task, **params}}

    headers = {"Content-Type": "application/json"}
    if API_KEY:
        headers["Authorization"] = f"Bearer {API_KEY}"

    for attempt in range(retries + 1):
        try:
            resp = requests.post(API_URL, json=payload, headers=headers, timeout=timeout)
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


# ============== RFD3 Scaffolding ==============

def run_scaffold_rfd3(motif_pdb: str, num_designs: int, seed: int = 42,
                      use_hotspots: bool = False) -> List[Dict[str, Any]]:
    """Run RFD3 with H-bond motif scaffolding.

    Key: Uses `unindex` for multi-island motif (3 segments).
    Does NOT fix metal coordinators - lets RFD3 design for Tb's CN=8-9.

    If use_hotspots=True, adds metal as hotspot to force protein contact.
    """
    params = {
        "pdb_content": motif_pdb,
        "contig": BEST_PARAMS["contig"],
        # Multi-island scaffolding - original chain A numbering
        "unindex": "A69-70,A214-216,A238",
        # Fix only protein H-bond residues (NOT metal coordination)
        "select_fixed_atoms": {
            "A69": "all", "A70": "all",
            "A214": "all", "A215": "all", "A216": "all",
            "A238": "all",
            "X1": "all",  # Fix metal position
            "L1": "all",  # Fix ligand position
        },
        # Tell RFD3 about metal and ligand
        "ligand": "CIT,TB",
        # RASA conditioning - bury metal aggressively to force coordination
        "select_buried": {"X1": "all", "L1": "O2,O5,O7"},
        # H-bond conditioning for citrate selectivity
        "select_hbond_acceptor": {"L1": "O1,O3,O4,O6"},
        # CRITICAL: Position diffusion cloud around metal/ligand center of mass
        # Without this, protein generates at random position and folds away from metal
        "infer_ori_strategy": "com",
        # CFG for conditioning - higher for scaffolding
        "use_classifier_free_guidance": True,
        "cfg_scale": 2.5,  # Higher CFG for better conditioning
        "num_timesteps": BEST_PARAMS["num_timesteps"],
        "step_scale": BEST_PARAMS["step_scale"],
        "num_designs": num_designs,
        "seed": seed,
    }

    # Add hotspots to force protein-metal contact
    if use_hotspots:
        params["hotspots"] = ["X1", "L1"]
        print(f"  [Scaffold RFD3] Using hotspots to force metal contact")

    print(f"  [Scaffold RFD3] contig={params['contig']}, unindex={params['unindex']}, "
          f"num_designs={num_designs}")
    result = call_api("rfd3", params)

    if result.get("status") != "completed":
        error = result.get("error", "unknown")
        print(f"  [Scaffold RFD3] FAILED: {error}")
        return []

    designs = result.get("result", {}).get("designs", [])
    print(f"  [Scaffold RFD3] Got {len(designs)} backbones")
    return designs


# ============== Tiered Sequence Strategy ==============

def get_seq_count_for_ptm(ptm: float) -> int:
    """Return number of sequences based on backbone pTM quality."""
    if ptm is None or ptm < 0.65:
        return 0  # Skip low quality backbones
    elif ptm < 0.75:
        return 2  # Low quality
    elif ptm < 0.85:
        return 4  # Medium quality
    else:
        return 6  # High quality


# ============== Analysis ==============

def compute_alanine_pct(sequence: str) -> float:
    """Compute alanine percentage from a sequence string."""
    if not sequence:
        return 0.0
    return round(sequence.count("A") / len(sequence) * 100, 1)


def analyze_design_metal(pdb_content: str) -> Dict[str, Any]:
    """Analyze metal coordination using analyze_design task.

    Uses same parameters as production_runner.py.
    """
    params = {
        "pdb_content": pdb_content,
        "metal_type": "TB",
        "ligand_name": "CIT",
        "design_type": "metal",
        "filter_tier": "standard",
    }
    result = call_api("analyze_design", params)

    if result.get("status") != "completed":
        print(f"    [Analyze] FAILED: {result.get('error', 'unknown')}")
        return {"coordination_number": None, "geometry_rmsd": None, "filter_passed": None}

    analysis = result.get("result", {})
    metrics = analysis.get("metrics", {})

    # Extract metrics from proper location (production_runner pattern)
    cn = metrics.get("coordination_number")
    geo_rmsd = metrics.get("geometry_rmsd")
    filter_passed = analysis.get("filter_passed")
    failed_filters = [
        f"{f.get('metric')}={f.get('value')}" for f in analysis.get("failed_filters", [])
    ]

    return {
        "coordination_number": cn,
        "geometry_rmsd": geo_rmsd,
        "filter_passed": filter_passed,
        "failed_filters": failed_filters,
    }


# ============== Main Pipeline ==============

def run_scaffold_batch(
    motif_pdb: str,
    num_backbones: int,
    batch_name: str,
    seed: int = 42,
    use_hotspots: bool = False,
) -> BatchResult:
    """Run full scaffolding batch: RFD3 -> MPNN -> RF3 -> Analyze."""

    result = BatchResult(name=batch_name)
    t0 = time.time()

    # ---- Step 1: RFD3 Scaffolding ----
    print(f"\n{'='*60}")
    print(f"  [{batch_name}] Step 1: RFD3 scaffolding ({num_backbones} backbones)")
    print(f"{'='*60}")

    t_rfd3 = time.time()
    designs = run_scaffold_rfd3(motif_pdb, num_designs=num_backbones, seed=seed,
                                use_hotspots=use_hotspots)

    backbones = []
    for i, d in enumerate(designs):
        content = d.get("content", d) if isinstance(d, dict) else d
        # Extract backbone confidence if available
        ptm = d.get("ptm") if isinstance(d, dict) else None
        plddt = d.get("plddt") if isinstance(d, dict) else None
        backbones.append({
            "id": f"b{i:03d}",
            "content": content,
            "ptm": ptm,
            "plddt": plddt,
        })

    result.backbones = backbones
    result.timings["rfd3"] = time.time() - t_rfd3
    print(f"  [{batch_name}] RFD3 complete: {len(backbones)} backbones "
          f"({result.timings['rfd3']:.0f}s)")

    if not backbones:
        result.errors.append("No backbones generated")
        return result

    # ---- Step 2: MPNN (tiered sequences per backbone) ----
    print(f"\n{'='*60}")
    print(f"  [{batch_name}] Step 2: MPNN sequence design (tiered)")
    print(f"{'='*60}")

    t_mpnn = time.time()
    all_seqs = []  # (backbone_id, backbone_ptm, backbone_plddt, sequence)

    for bb in backbones:
        # Determine sequence count based on pTM
        bb_ptm = bb.get("ptm")
        num_seqs = get_seq_count_for_ptm(bb_ptm) if bb_ptm else 4  # Default 4 if no pTM

        if num_seqs == 0:
            print(f"  [{batch_name}] Skipping backbone {bb['id']} (pTM={bb_ptm:.3f} < 0.65)")
            continue

        params = {
            "pdb_content": bb["content"],
            "num_sequences": num_seqs,
            "temperature": BEST_PARAMS["mpnn_temperature"],
            "bias_AA": BEST_PARAMS["mpnn_bias_AA"],
            "omit_AA": BEST_PARAMS["mpnn_omit_AA"],
            "ligand_mpnn_use_atom_context": 1,
            "pack_side_chains": True,
        }

        print(f"  [{batch_name}] MPNN {bb['id']} (pTM={bb_ptm or 'N/A'}) -> {num_seqs} seqs...")
        mpnn_result = call_api("mpnn", params)

        if mpnn_result.get("status") != "completed":
            err = f"MPNN failed for {bb['id']}: {mpnn_result.get('error', 'unknown')}"
            print(f"  [{batch_name}] {err}")
            result.errors.append(err)
            continue

        seqs = parse_sequences(mpnn_result)
        print(f"  [{batch_name}] Got {len(seqs)} sequences from {bb['id']}")
        for seq in seqs:
            all_seqs.append((bb["id"], bb.get("ptm"), bb.get("plddt"), seq))

    result.timings["mpnn"] = time.time() - t_mpnn
    print(f"  [{batch_name}] MPNN complete: {len(all_seqs)} sequences "
          f"({result.timings['mpnn']:.0f}s)")

    if not all_seqs:
        result.errors.append("No sequences designed")
        return result

    # ---- Step 3: RF3 (per sequence) ----
    print(f"\n{'='*60}")
    print(f"  [{batch_name}] Step 3: RF3 structure prediction")
    print(f"{'='*60}")

    t_rf3 = time.time()
    rf3_results = []

    for i, (bb_id, bb_ptm, bb_plddt, seq) in enumerate(all_seqs):
        params = {
            "sequence": seq,
            "name": f"{batch_name}_seq_{i:04d}",
            "ligand_smiles": CITRATE_SMILES,
            "metal": "TB",
        }

        print(f"  [{batch_name}] RF3 {i+1}/{len(all_seqs)} (backbone {bb_id})...")
        rf3_result = call_api("rf3", params)

        if rf3_result.get("status") != "completed":
            err = f"RF3 failed for seq {i+1}: {rf3_result.get('error', 'unknown')}"
            print(f"  [{batch_name}] {err}")
            result.errors.append(err)
            rf3_results.append((bb_id, bb_ptm, bb_plddt, seq, None))
            continue

        metrics = extract_rf3_metrics(rf3_result)
        print(f"  [{batch_name}] RF3 seq {i+1}: pTM={metrics['ptm']}, "
              f"pLDDT={metrics['plddt']}, PAE={metrics['pae']}")
        rf3_results.append((bb_id, bb_ptm, bb_plddt, seq, metrics))

    result.timings["rf3"] = time.time() - t_rf3
    completed_rf3 = sum(1 for _, _, _, _, m in rf3_results if m is not None)
    print(f"  [{batch_name}] RF3 complete: {completed_rf3}/{len(rf3_results)} "
          f"({result.timings['rf3']:.0f}s)")

    # ---- Step 4: Analyze with metal filter ----
    print(f"\n{'='*60}")
    print(f"  [{batch_name}] Step 4: Metal binding analysis")
    print(f"{'='*60}")

    t_analyze = time.time()

    for i, (bb_id, bb_ptm, bb_plddt, seq, rf3_metrics) in enumerate(rf3_results):
        seq_id = f"seq_{i:04d}"
        design = DesignResult(
            seq_id=seq_id,
            sequence=seq,
            backbone_id=bb_id,
            backbone_ptm=bb_ptm,
            backbone_plddt=bb_plddt,
        )

        if rf3_metrics is None:
            result.designs.append(design)
            continue

        design.rf3_pdb = rf3_metrics.get("pdb_content")
        design.ptm = rf3_metrics.get("ptm")
        design.plddt = rf3_metrics.get("plddt")
        design.pae = rf3_metrics.get("pae")

        # Compute sequence-derived metrics
        design.alanine_pct = compute_alanine_pct(seq)
        design.total_residues = len(seq)

        # Analyze metal coordination
        if design.rf3_pdb:
            analysis = analyze_design_metal(design.rf3_pdb)
            design.coordination_number = analysis.get("coordination_number")
            design.geometry_rmsd = analysis.get("geometry_rmsd")
            design.filter_passed = analysis.get("filter_passed")
            design.failed_filters = analysis.get("failed_filters")

        result.designs.append(design)

    result.timings["analyze"] = time.time() - t_analyze
    result.timings["total"] = time.time() - t0

    # Summary
    passing = [d for d in result.designs if d.filter_passed]
    print(f"\n{'='*60}")
    print(f"  [{batch_name}] SUMMARY")
    print(f"{'='*60}")
    print(f"  Backbones: {len(backbones)}")
    print(f"  Sequences: {len(result.designs)}")
    print(f"  Passing metal_tb_standard: {len(passing)} ({100*len(passing)/len(result.designs):.1f}%)")

    # CN distribution
    cn_values = [d.coordination_number for d in result.designs if d.coordination_number is not None]
    if cn_values:
        print(f"  CN range: {min(cn_values)}-{max(cn_values)}, avg: {statistics.mean(cn_values):.1f}")
        cn_ge_6 = sum(1 for cn in cn_values if cn >= 6)
        cn_ge_7 = sum(1 for cn in cn_values if cn >= 7)
        cn_ge_8 = sum(1 for cn in cn_values if cn >= 8)
        print(f"  CN>=6: {cn_ge_6}, CN>=7: {cn_ge_7}, CN>=8: {cn_ge_8}")

    print(f"  Total time: {result.timings['total']:.0f}s")

    return result


def save_results(result: BatchResult, output_dir: Path):
    """Save batch results to JSON and individual PDBs."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save summary JSON
    summary = {
        "name": result.name,
        "num_backbones": len(result.backbones),
        "num_designs": len(result.designs),
        "num_passing": len([d for d in result.designs if d.filter_passed]),
        "timings": result.timings,
        "errors": result.errors,
        "designs": [asdict(d) for d in result.designs],
    }

    # Remove PDB content from summary (too large)
    for d in summary["designs"]:
        d.pop("rf3_pdb", None)

    summary_path = output_dir / f"{result.name}_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Saved summary to {summary_path}")

    # Save passing PDBs
    pdb_dir = output_dir / f"{result.name}_pdbs"
    pdb_dir.mkdir(exist_ok=True)

    for d in result.designs:
        if d.filter_passed and d.rf3_pdb:
            pdb_path = pdb_dir / f"{d.seq_id}.pdb"
            with open(pdb_path, "w") as f:
                f.write(d.rf3_pdb)

    num_saved = len([d for d in result.designs if d.filter_passed and d.rf3_pdb])
    print(f"Saved {num_saved} passing PDBs to {pdb_dir}")


# ============== Main ==============

def main():
    parser = argparse.ArgumentParser(description="H-bond scaffolding E2E test")
    parser.add_argument("--pilot", action="store_true", help="Run pilot (5 backbones)")
    parser.add_argument("--batch", type=int, help="Run production batch N (50 backbones)")
    parser.add_argument("--subbatch", type=str, help="Run sub-batch N.M (10 backbones), e.g., 2.1 for batch 2 sub 1")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args()

    # Load motif PDB
    if not MOTIF_PDB_PATH.exists():
        print(f"ERROR: Motif PDB not found at {MOTIF_PDB_PATH}")
        print("Run US-001 first to create the H-bond motif.")
        sys.exit(1)

    with open(MOTIF_PDB_PATH) as f:
        motif_pdb = f.read()

    print(f"Loaded motif PDB: {len(motif_pdb)} chars")
    print(f"Output directory: {OUTPUT_DIR}")

    if args.pilot:
        # Quick test with 5 backbones - try with hotspots
        result = run_scaffold_batch(
            motif_pdb=motif_pdb,
            num_backbones=5,
            batch_name="pilot_hotspots",
            seed=args.seed,
            use_hotspots=True,  # Force metal contact
        )
        save_results(result, OUTPUT_DIR)

    elif args.subbatch is not None:
        # Sub-batch (10 backbones) to avoid hitting same endpoint
        # Format: N.M where N is batch number, M is sub-batch (1-5)
        parts = args.subbatch.split(".")
        batch_num = int(parts[0])
        sub_num = int(parts[1]) if len(parts) > 1 else 1

        # Each sub-batch gets unique seed: base + batch*1000 + sub*100
        sub_seed = args.seed + (batch_num - 1) * 1000 + (sub_num - 1) * 100
        batch_name = f"batch_{batch_num:02d}_sub{sub_num}"

        print(f"Running sub-batch {batch_num}.{sub_num} (10 backbones, seed={sub_seed})")
        result = run_scaffold_batch(
            motif_pdb=motif_pdb,
            num_backbones=10,
            batch_name=batch_name,
            seed=sub_seed,
            use_hotspots=True,
        )
        save_results(result, OUTPUT_DIR)

    elif args.batch is not None:
        # Production batch (50 backbones) with hotspots
        batch_seed = args.seed + (args.batch - 1) * 1000
        result = run_scaffold_batch(
            motif_pdb=motif_pdb,
            num_backbones=50,
            batch_name=f"batch_{args.batch:02d}",
            seed=batch_seed,
            use_hotspots=True,  # Use hotspots for metal contact
        )
        save_results(result, OUTPUT_DIR)

    else:
        parser.print_help()
        print("\nExamples:")
        print("  python test_hbond_scaffolding.py --pilot          # Quick test (5 backbones)")
        print("  python test_hbond_scaffolding.py --batch 1        # Full batch (50 backbones)")
        print("  python test_hbond_scaffolding.py --subbatch 2.1   # Sub-batch (10 backbones)")
        print("\nFor 400 designs (40 sub-batches of 10):")
        print("  Run --subbatch 1.1 through 8.5 (8 batches × 5 subs each)")


if __name__ == "__main__":
    main()
