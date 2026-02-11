"""
E2E Test: Ln-Citrate Production Pipeline

Mirrors the frontend NL pipeline exactly, adding scout filter,
stability filtering, metal+ligand coordination analysis, design
history persistence, and lesson detection.

Pipeline:
  1. RFD3 backbone generation (metal_binding_design single)
  2. Scout filter (pre-filter + optional MPNN+RF3 per backbone)
  3. MPNN sequence design (with metal bias)
  4. RF3 validation (with citrate SMILES + TB metal context)
  5. Stability filter (pTM, pLDDT, alanine %)
  6. Analyze design (UnifiedDesignAnalyzer + FILTER_PRESETS)
  7. Save history (non-fatal)
  8. Check lessons (non-fatal)
  9. Report — summary table, production readiness, top designs

Usage:
    cd backend/serverless
    python test_production_pipeline.py
"""

import json
import os
import statistics
import sys
import time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import requests

# ============== Configuration ==============

API_URL = "http://localhost:8000/runsync"

# Scale controls — adjust for production runs
NUM_DESIGNS = 5              # Backbones per RFD3 call (production: 20-50)
NUM_SEQS_PER_BACKBONE = 4    # Sequences per backbone
CONTIG = "110-150"           # Protein length range

# Scout thresholds (matches frontend defaults)
SCOUT_PTM_THRESHOLD = 0.6
SCOUT_PLDDT_THRESHOLD = 0.65

# Stability thresholds (sequence-level pre-filter)
STABILITY_PTM_MIN = 0.6
STABILITY_PLDDT_MIN = 0.65
STABILITY_ALA_MAX = 25.0     # Percent

# Metal/ligand context
METAL = "TB"
LIGAND = "CIT"
CITRATE_SMILES = "OC(=O)CC(O)(CC(O)=O)C(O)=O"

# API timeouts (seconds)
TIMEOUTS = {
    "metal_binding_design": 1800,
    "scout_filter": 1200,
    "mpnn": 600,
    "rf3": 600,
    "analyze_design": 120,
    "save_design_history": 60,
    "check_lessons": 60,
}

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "output_production_pipeline")


# ============== Data Classes ==============

@dataclass
class BackboneResult:
    id: str
    pdb_content: str
    pre_filter_passed: bool = True
    pre_filter_checks: dict = field(default_factory=dict)
    scout_ptm: Optional[float] = None
    scout_plddt: Optional[float] = None
    scout_passed: Optional[bool] = None


@dataclass
class DesignResult:
    seq_id: str
    sequence: str
    backbone_id: str
    # RF3 metrics
    rf3_pdb: Optional[str] = None
    ptm: Optional[float] = None
    plddt: Optional[float] = None
    pae: Optional[float] = None
    # Stability
    alanine_pct: Optional[float] = None
    total_residues: Optional[int] = None
    stability_passed: bool = False
    # Coordination (from analyze_design)
    coordination_number: Optional[float] = None
    geometry_rmsd: Optional[float] = None
    filter_preset: Optional[str] = None
    filter_passed: Optional[bool] = None
    failed_filters: list = field(default_factory=list)
    # Ligand analysis
    ligand_contacts: Optional[int] = None
    hsab_score: Optional[float] = None


@dataclass
class PipelineResult:
    backbones_generated: int = 0
    backbones_after_scout: int = 0
    sequences_designed: int = 0
    sequences_after_stability: int = 0
    designs_analyzed: int = 0
    designs_passed_metal_binding: int = 0
    designs_passed_strict: int = 0
    # Collections
    backbones: List[BackboneResult] = field(default_factory=list)
    designs: List[DesignResult] = field(default_factory=list)
    timings: Dict[str, float] = field(default_factory=dict)
    errors: List[str] = field(default_factory=list)
    # History/lessons
    history_saved: bool = False
    history_run_id: Optional[str] = None
    lessons_detected: list = field(default_factory=list)


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


# ============== Utility Functions ==============

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


def extract_rf3_metrics(rf3_result: Dict[str, Any]) -> Dict[str, Any]:
    """Extract pTM, pLDDT, PAE from RF3 result."""
    metrics = {"ptm": None, "plddt": None, "pae": None, "pdb_content": None}
    result = rf3_result.get("result", {})

    predictions = result.get("predictions", [])
    if predictions:
        metrics["pdb_content"] = predictions[0].get("content")

    conf = result.get("confidences") or {}
    summary = conf.get("summary_confidences", conf)

    metrics["ptm"] = summary.get("ptm")
    metrics["plddt"] = summary.get("overall_plddt")
    metrics["pae"] = summary.get("overall_pae")

    return metrics


def compute_alanine_pct(sequence: str) -> float:
    """Compute alanine percentage from a sequence string."""
    if not sequence:
        return 0.0
    return round(sequence.count("A") / len(sequence) * 100, 1)


def safe_mean(values: List[Optional[float]]) -> Optional[float]:
    """Mean of non-None values, or None if empty."""
    valid = [v for v in values if v is not None]
    return round(statistics.mean(valid), 4) if valid else None


def composite_score(design: DesignResult) -> float:
    """Compute composite ranking score for a design."""
    ptm = design.ptm or 0
    plddt = design.plddt or 0
    geo_rmsd = design.geometry_rmsd if design.geometry_rmsd is not None else 2.0
    return ptm * 0.4 + plddt * 0.3 + max(0, (1 - geo_rmsd / 2)) * 0.3


# ============== Pipeline Steps ==============

def step1_rfd3_backbone_generation(result: PipelineResult) -> PipelineResult:
    """Step 1: RFD3 backbone generation via metal_binding_design single mode."""
    print(f"\n{'='*64}")
    print(f"  Step 1: RFD3 Backbone Generation")
    print(f"  contig={CONTIG}, num_designs={NUM_DESIGNS}")
    print(f"{'='*64}")

    t0 = time.time()

    params = {
        "mode": "single",
        "metal": METAL,
        "ligand": LIGAND,
        "contig": CONTIG,
        "cfg_scale": 2.0,
        "num_designs": NUM_DESIGNS,
        "bury_ligand": True,
    }

    api_result = call_api("metal_binding_design", params)

    if api_result.get("status") != "completed":
        err = f"RFD3 failed: {api_result.get('error', 'unknown')}"
        print(f"  ERROR: {err}")
        result.errors.append(err)
        result.timings["rfd3"] = time.time() - t0
        return result

    designs = api_result.get("result", {}).get("designs", [])
    print(f"  Got {len(designs)} backbones")

    for i, d in enumerate(designs):
        content = d.get("content", d) if isinstance(d, dict) else d
        bb = BackboneResult(id=f"bb_{i:03d}", pdb_content=content)
        result.backbones.append(bb)

    result.backbones_generated = len(result.backbones)
    result.timings["rfd3"] = time.time() - t0
    print(f"  Complete: {result.backbones_generated} backbones ({result.timings['rfd3']:.0f}s)")

    return result


def step2_scout_filter(result: PipelineResult) -> PipelineResult:
    """Step 2: Scout filter — pre-filter + optional MPNN+RF3 per backbone.

    For metal_binding_design single mode, we use pre_filter_only=True
    (matches frontend logic where metal backbones skip GPU scout).
    """
    print(f"\n{'='*64}")
    print(f"  Step 2: Scout Filter")
    print(f"  {result.backbones_generated} backbones, "
          f"pTM>={SCOUT_PTM_THRESHOLD}, pLDDT>={SCOUT_PLDDT_THRESHOLD}")
    print(f"{'='*64}")

    if not result.backbones:
        print("  SKIP: No backbones to filter")
        return result

    t0 = time.time()

    backbone_pdbs = [bb.pdb_content for bb in result.backbones]

    params = {
        "backbone_pdbs": backbone_pdbs,
        "ptm_threshold": SCOUT_PTM_THRESHOLD,
        "plddt_threshold": SCOUT_PLDDT_THRESHOLD,
        "target_metal": METAL,
        "ligand_smiles": CITRATE_SMILES,
        "ligand_name": LIGAND,
        "pre_filter_only": True,
    }

    api_result = call_api("scout_filter", params)

    if api_result.get("status") != "completed":
        err = f"Scout filter failed: {api_result.get('error', 'unknown')}"
        print(f"  WARNING: {err}")
        print("  Proceeding with all backbones (scout is non-blocking)")
        result.errors.append(err)
        result.backbones_after_scout = result.backbones_generated
        result.timings["scout"] = time.time() - t0
        return result

    scout_data = api_result.get("result", {})
    scout_results = scout_data.get("scout_results", [])
    pre_filter_results = scout_data.get("pre_filter_results", [])

    # Update backbone results with scout data
    for i, bb in enumerate(result.backbones):
        # Pre-filter results
        if i < len(pre_filter_results):
            pf = pre_filter_results[i]
            bb.pre_filter_passed = pf.get("passed", True)
            bb.pre_filter_checks = pf.get("checks", {})

        # Scout results (if GPU scout ran)
        if i < len(scout_results):
            sr = scout_results[i]
            bb.scout_ptm = sr.get("ptm")
            bb.scout_plddt = sr.get("plddt")
            bb.scout_passed = sr.get("passed", True)

    # Filter to passing backbones
    passing = [bb for bb in result.backbones if bb.pre_filter_passed and
               (bb.scout_passed is None or bb.scout_passed)]
    failed = result.backbones_generated - len(passing)

    result.backbones_after_scout = len(passing)
    result.timings["scout"] = time.time() - t0

    print(f"  Pre-filter: {scout_data.get('pre_filter_passed', '?')} passed, "
          f"{scout_data.get('pre_filter_failed', '?')} failed")
    if not scout_data.get("pre_filter_only", True):
        print(f"  GPU scout: {scout_data.get('passing_count', '?')} passed")
    print(f"  Result: {len(passing)}/{result.backbones_generated} backbones pass "
          f"({result.timings['scout']:.0f}s)")

    # Replace backbones list with only passing ones
    result.backbones = passing

    return result


def step3_mpnn_sequence_design(result: PipelineResult) -> PipelineResult:
    """Step 3: MPNN sequence design with metal bias on passing backbones."""
    print(f"\n{'='*64}")
    print(f"  Step 3: MPNN Sequence Design")
    print(f"  {len(result.backbones)} backbones x {NUM_SEQS_PER_BACKBONE} seqs each")
    print(f"{'='*64}")

    if not result.backbones:
        print("  SKIP: No backbones for MPNN")
        return result

    t0 = time.time()
    all_seqs = []  # (backbone_id, sequence)

    for bb in result.backbones:
        params = {
            "pdb_content": bb.pdb_content,
            "num_sequences": NUM_SEQS_PER_BACKBONE,
            "temperature": 0.1,
            "ligand_mpnn_use_atom_context": 1,
            "bias_AA": "A:-2.0",
            "omit_AA": "C",
            "pack_side_chains": True,
        }

        print(f"  MPNN for {bb.id}...")
        api_result = call_api("mpnn", params)

        if api_result.get("status") != "completed":
            err = f"MPNN failed for {bb.id}: {api_result.get('error', 'unknown')}"
            print(f"  ERROR: {err}")
            result.errors.append(err)
            continue

        seqs = parse_sequences(api_result)
        print(f"    Got {len(seqs)} sequences")
        for seq in seqs:
            all_seqs.append((bb.id, seq))

    result.sequences_designed = len(all_seqs)
    result.timings["mpnn"] = time.time() - t0
    print(f"  Complete: {result.sequences_designed} sequences ({result.timings['mpnn']:.0f}s)")

    # Store sequences temporarily for next step (as partial DesignResult)
    for i, (bb_id, seq) in enumerate(all_seqs):
        design = DesignResult(
            seq_id=f"seq_{i+1:03d}",
            sequence=seq,
            backbone_id=bb_id,
            alanine_pct=compute_alanine_pct(seq),
            total_residues=len(seq),
        )
        result.designs.append(design)

    return result


def step4_rf3_validation(result: PipelineResult) -> PipelineResult:
    """Step 4: RF3 structure prediction with citrate + TB context."""
    print(f"\n{'='*64}")
    print(f"  Step 4: RF3 Validation")
    print(f"  {len(result.designs)} sequences")
    print(f"{'='*64}")

    if not result.designs:
        print("  SKIP: No sequences for RF3")
        return result

    t0 = time.time()
    completed = 0

    for i, design in enumerate(result.designs):
        params = {
            "sequence": design.sequence,
            "name": f"prod_{design.seq_id}",
            "ligand_smiles": CITRATE_SMILES,
            "metal": METAL,
        }

        print(f"  RF3 {i+1}/{len(result.designs)} ({design.seq_id}, "
              f"backbone {design.backbone_id})...")
        api_result = call_api("rf3", params)

        if api_result.get("status") != "completed":
            err = f"RF3 failed for {design.seq_id}: {api_result.get('error', 'unknown')}"
            print(f"    ERROR: {err}")
            result.errors.append(err)
            continue

        metrics = extract_rf3_metrics(api_result)
        design.rf3_pdb = metrics.get("pdb_content")
        design.ptm = metrics.get("ptm")
        design.plddt = metrics.get("plddt")
        design.pae = metrics.get("pae")
        completed += 1

        print(f"    pTM={design.ptm}, pLDDT={design.plddt}, PAE={design.pae}")

    result.timings["rf3"] = time.time() - t0
    print(f"  Complete: {completed}/{len(result.designs)} validated ({result.timings['rf3']:.0f}s)")

    return result


def step5_stability_filter(result: PipelineResult) -> PipelineResult:
    """Step 5: Stability filter — remove weak sequences before expensive analysis."""
    print(f"\n{'='*64}")
    print(f"  Step 5: Stability Filter")
    print(f"  pTM>={STABILITY_PTM_MIN}, pLDDT>={STABILITY_PLDDT_MIN}, "
          f"ala%<={STABILITY_ALA_MAX}")
    print(f"{'='*64}")

    before = len(result.designs)
    passing = []

    for design in result.designs:
        # Must have RF3 results
        if design.ptm is None:
            print(f"  {design.seq_id}: FAIL (no RF3 result)")
            continue

        failed = []
        if design.ptm < STABILITY_PTM_MIN:
            failed.append(f"pTM={design.ptm:.3f}<{STABILITY_PTM_MIN}")
        if design.plddt is not None and design.plddt < STABILITY_PLDDT_MIN:
            failed.append(f"pLDDT={design.plddt:.3f}<{STABILITY_PLDDT_MIN}")
        if design.alanine_pct is not None and design.alanine_pct > STABILITY_ALA_MAX:
            failed.append(f"ala%={design.alanine_pct:.1f}>{STABILITY_ALA_MAX}")

        if failed:
            print(f"  {design.seq_id}: FAIL ({', '.join(failed)})")
            continue

        design.stability_passed = True
        passing.append(design)
        print(f"  {design.seq_id}: PASS (pTM={design.ptm:.3f}, "
              f"pLDDT={design.plddt:.3f}, ala%={design.alanine_pct:.1f})")

    result.designs = passing
    result.sequences_after_stability = len(passing)
    print(f"\n  Result: {len(passing)}/{before} pass stability filter")

    return result


def step6_analyze_design(result: PipelineResult) -> PipelineResult:
    """Step 6: Full design analysis via analyze_design task.

    Uses UnifiedDesignAnalyzer + FILTER_PRESETS for coordination geometry,
    contacts, and metal_binding / metal_binding_strict preset evaluation.
    """
    print(f"\n{'='*64}")
    print(f"  Step 6: Design Analysis (UnifiedDesignAnalyzer)")
    print(f"  {len(result.designs)} designs, preset=metal_binding")
    print(f"{'='*64}")

    if not result.designs:
        print("  SKIP: No designs to analyze")
        return result

    t0 = time.time()
    metal_pass = 0
    strict_pass = 0

    for design in result.designs:
        if not design.rf3_pdb:
            print(f"  {design.seq_id}: SKIP (no RF3 PDB)")
            continue

        # Run analyze_design with metal context
        params = {
            "pdb_content": design.rf3_pdb,
            "metal_type": METAL,
            "ligand_name": LIGAND,
            "design_type": "metal",
        }

        print(f"  Analyzing {design.seq_id}...")
        api_result = call_api("analyze_design", params)

        if api_result.get("status") != "completed":
            err = f"Analysis failed for {design.seq_id}: {api_result.get('error', 'unknown')}"
            print(f"    ERROR: {err}")
            result.errors.append(err)
            continue

        analysis = api_result.get("result", {})
        metrics = analysis.get("metrics", {})

        # Extract coordination metrics
        design.coordination_number = metrics.get("coordination_number")
        design.geometry_rmsd = metrics.get("geometry_rmsd")
        design.ligand_contacts = metrics.get("total_contacts")
        design.hsab_score = metrics.get("hsab_score")

        # Filter evaluation — metal_binding preset
        design.filter_preset = analysis.get("filter_preset", "metal_binding")
        design.filter_passed = analysis.get("filter_passed")
        design.failed_filters = [
            f"{f.get('metric')}={f.get('value')}" for f in analysis.get("failed_filters", [])
        ]

        result.designs_analyzed += 1

        if design.filter_passed:
            metal_pass += 1

        # Also check strict preset manually
        cn = design.coordination_number or 0
        geo = design.geometry_rmsd if design.geometry_rmsd is not None else 999
        plddt = design.plddt or 0
        if cn >= 8 and geo <= 0.8 and plddt >= 0.85:
            strict_pass += 1

        print(f"    CN={design.coordination_number}, geo_rmsd={design.geometry_rmsd}, "
              f"contacts={design.ligand_contacts}, "
              f"filter={'PASS' if design.filter_passed else 'FAIL'}"
              f"{' (' + ', '.join(design.failed_filters) + ')' if design.failed_filters else ''}")

    result.designs_passed_metal_binding = metal_pass
    result.designs_passed_strict = strict_pass
    result.timings["analyze"] = time.time() - t0
    print(f"\n  Complete: {result.designs_analyzed} analyzed, "
          f"{metal_pass} pass metal_binding, "
          f"{strict_pass} pass metal_binding_strict "
          f"({result.timings['analyze']:.0f}s)")

    return result


def step7_save_history(result: PipelineResult) -> PipelineResult:
    """Step 7: Save design history (non-fatal)."""
    print(f"\n{'='*64}")
    print(f"  Step 7: Save Design History")
    print(f"{'='*64}")

    t0 = time.time()

    # Compile design metrics for history
    design_metrics = {
        "backbones_generated": result.backbones_generated,
        "backbones_after_scout": result.backbones_after_scout,
        "sequences_designed": result.sequences_designed,
        "sequences_after_stability": result.sequences_after_stability,
        "designs_analyzed": result.designs_analyzed,
        "designs_passed_metal_binding": result.designs_passed_metal_binding,
        "designs_passed_strict": result.designs_passed_strict,
        "avg_ptm": safe_mean([d.ptm for d in result.designs]),
        "avg_plddt": safe_mean([d.plddt for d in result.designs]),
        "avg_coordination": safe_mean([d.coordination_number for d in result.designs]),
        "avg_geometry_rmsd": safe_mean([d.geometry_rmsd for d in result.designs]),
        "top_designs": [
            {
                "seq_id": d.seq_id,
                "ptm": d.ptm,
                "plddt": d.plddt,
                "coordination_number": d.coordination_number,
                "geometry_rmsd": d.geometry_rmsd,
                "filter_passed": d.filter_passed,
            }
            for d in sorted(result.designs, key=composite_score, reverse=True)[:5]
        ],
    }

    params = {
        "session_name": "production_pipeline_test",
        "design_params": {
            "metal": METAL,
            "ligand": LIGAND,
            "contig": CONTIG,
            "num_designs": NUM_DESIGNS,
            "num_seqs_per_backbone": NUM_SEQS_PER_BACKBONE,
            "citrate_smiles": CITRATE_SMILES,
        },
        "design_outputs": {
            "num_backbones": result.backbones_generated,
            "num_sequences": result.sequences_designed,
            "num_analyzed": result.designs_analyzed,
        },
        "design_metrics": design_metrics,
    }

    api_result = call_api("save_design_history", params)

    if api_result.get("status") == "completed":
        run_id = api_result.get("result", {}).get("run_id", "unknown")
        result.history_saved = True
        result.history_run_id = run_id
        print(f"  Saved: run_id={run_id}")
    else:
        err = f"Save history failed: {api_result.get('error', 'unknown')}"
        print(f"  WARNING: {err} (non-fatal)")
        result.errors.append(err)

    result.timings["save_history"] = time.time() - t0
    return result


def step8_check_lessons(result: PipelineResult) -> PipelineResult:
    """Step 8: Check for lesson triggers (non-fatal)."""
    print(f"\n{'='*64}")
    print(f"  Step 8: Check Lessons")
    print(f"{'='*64}")

    t0 = time.time()

    params = {
        "result": {
            "metal": METAL,
            "ligand": LIGAND,
            "avg_ptm": safe_mean([d.ptm for d in result.designs]),
            "avg_plddt": safe_mean([d.plddt for d in result.designs]),
            "pass_rate": (result.designs_passed_metal_binding / result.designs_analyzed * 100
                         if result.designs_analyzed > 0 else 0),
            "num_designs": result.designs_analyzed,
        },
    }

    api_result = call_api("check_lessons", params)

    if api_result.get("status") == "completed":
        lesson_data = api_result.get("result", {})
        if lesson_data.get("trigger_detected"):
            trigger = lesson_data.get("trigger", {})
            result.lessons_detected.append(trigger)
            print(f"  TRIGGER: {trigger.get('type', '?')} — {trigger.get('description', '?')}")
        else:
            print(f"  No lesson triggers detected "
                  f"(analyzed {lesson_data.get('history_count', '?')} recent designs)")
    else:
        err = f"Check lessons failed: {api_result.get('error', 'unknown')}"
        print(f"  WARNING: {err} (non-fatal)")
        result.errors.append(err)

    result.timings["check_lessons"] = time.time() - t0
    return result


# ============== Report ==============

def format_report(result: PipelineResult) -> str:
    """Format a human-readable production pipeline report."""

    def fmt(val, decimals=3):
        if val is None:
            return "N/A"
        if isinstance(val, float):
            return f"{val:.{decimals}f}"
        return str(val)

    lines = [
        "=" * 72,
        "      Ln-CITRATE PRODUCTION PIPELINE REPORT",
        "=" * 72,
        "",
        "PIPELINE FUNNEL",
        f"  Backbones generated:     {result.backbones_generated}",
        f"  After scout filter:      {result.backbones_after_scout}",
        f"  Sequences designed:      {result.sequences_designed}",
        f"  After stability filter:  {result.sequences_after_stability}",
        f"  Designs analyzed:        {result.designs_analyzed}",
        f"  Pass metal_binding:      {result.designs_passed_metal_binding}",
        f"  Pass metal_binding_strict: {result.designs_passed_strict}",
        "",
    ]

    # Production readiness
    if result.designs_analyzed > 0:
        pass_rate = result.designs_passed_metal_binding / result.designs_analyzed * 100
        lines.append(f"PRODUCTION READINESS: {pass_rate:.1f}% pass metal_binding preset")
        if pass_rate >= 50:
            lines.append("  -> Ready for scale-up (>=50% pass rate)")
        elif pass_rate >= 20:
            lines.append("  -> Marginal — consider parameter tuning before scale-up")
        else:
            lines.append("  -> NOT ready — significant optimization needed")
    else:
        lines.append("PRODUCTION READINESS: No designs analyzed")
    lines.append("")

    # Per-design table
    if result.designs:
        lines.append("DESIGN DETAILS")
        header = (f"  {'ID':>8s}  {'Backbone':>8s}  {'pTM':>6s}  {'pLDDT':>6s}  "
                  f"{'PAE':>6s}  {'Ala%':>5s}  {'CN':>4s}  {'geoRMSD':>7s}  "
                  f"{'Filter':>6s}  {'Score':>6s}")
        lines.append(header)
        lines.append("  " + "-" * (len(header) - 2))

        ranked = sorted(result.designs, key=composite_score, reverse=True)
        for d in ranked:
            score = composite_score(d)
            filt = "PASS" if d.filter_passed else "FAIL"
            lines.append(
                f"  {d.seq_id:>8s}  {d.backbone_id:>8s}  {fmt(d.ptm):>6s}  "
                f"{fmt(d.plddt):>6s}  {fmt(d.pae, 2):>6s}  "
                f"{fmt(d.alanine_pct, 1):>5s}  "
                f"{fmt(d.coordination_number, 0):>4s}  "
                f"{fmt(d.geometry_rmsd, 2):>7s}  {filt:>6s}  {score:.3f}"
            )
        lines.append("")

    # Top 5
    if result.designs:
        lines.append("TOP 5 DESIGNS (by composite score)")
        ranked = sorted(result.designs, key=composite_score, reverse=True)[:5]
        for i, d in enumerate(ranked):
            score = composite_score(d)
            lines.append(
                f"  #{i+1}  {d.seq_id} (backbone {d.backbone_id}) "
                f"score={score:.3f}  pTM={fmt(d.ptm)}  CN={fmt(d.coordination_number, 0)}  "
                f"geoRMSD={fmt(d.geometry_rmsd, 2)}"
            )
        lines.append("")

    # Aggregate stats
    lines.append("AGGREGATE METRICS")
    lines.append(f"  Avg pTM:            {fmt(safe_mean([d.ptm for d in result.designs]))}")
    lines.append(f"  Avg pLDDT:          {fmt(safe_mean([d.plddt for d in result.designs]))}")
    lines.append(f"  Avg PAE:            {fmt(safe_mean([d.pae for d in result.designs]), 2)}")
    lines.append(f"  Avg Alanine %:      {fmt(safe_mean([d.alanine_pct for d in result.designs]), 1)}")
    lines.append(f"  Avg Coordination #: {fmt(safe_mean([d.coordination_number for d in result.designs]), 1)}")
    lines.append(f"  Avg Geometry RMSD:  {fmt(safe_mean([d.geometry_rmsd for d in result.designs]), 2)}")
    lines.append("")

    # Timings
    lines.append("TIMING BREAKDOWN (seconds)")
    total = sum(result.timings.values())
    for step_name, t in result.timings.items():
        pct = t / total * 100 if total > 0 else 0
        lines.append(f"  {step_name:20s}  {t:>7.0f}s  ({pct:.1f}%)")
    lines.append(f"  {'TOTAL':20s}  {total:>7.0f}s")
    lines.append("")

    # History/lessons
    lines.append("HISTORY & LESSONS")
    lines.append(f"  History saved:    {'Yes (run_id=' + result.history_run_id + ')' if result.history_saved else 'No'}")
    if result.lessons_detected:
        for lesson in result.lessons_detected:
            lines.append(f"  Lesson trigger:   [{lesson.get('type', '?')}] {lesson.get('description', '?')}")
    else:
        lines.append("  Lesson triggers:  None")
    lines.append("")

    # Errors
    if result.errors:
        lines.append("ERRORS")
        for e in result.errors:
            lines.append(f"  - {e}")
        lines.append("")

    lines.append("=" * 72)
    return "\n".join(lines)


# ============== Output ==============

def save_outputs(result: PipelineResult):
    """Save all pipeline outputs to disk."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Save backbones
    bb_dir = os.path.join(OUTPUT_DIR, "backbones")
    os.makedirs(bb_dir, exist_ok=True)
    for bb in result.backbones:
        path = os.path.join(bb_dir, f"{bb.id}.pdb")
        with open(path, "w") as f:
            f.write(bb.pdb_content)

    # Save sequences as FASTA
    fasta_path = os.path.join(OUTPUT_DIR, "sequences.fasta")
    with open(fasta_path, "w") as f:
        for d in result.designs:
            f.write(f">{d.seq_id} backbone={d.backbone_id}\n{d.sequence}\n")

    # Save RF3 PDBs
    rf3_dir = os.path.join(OUTPUT_DIR, "rf3_predictions")
    os.makedirs(rf3_dir, exist_ok=True)
    for d in result.designs:
        if d.rf3_pdb:
            path = os.path.join(rf3_dir, f"rf3_{d.seq_id}.pdb")
            with open(path, "w") as f:
                f.write(d.rf3_pdb)

    # Save pipeline summary JSON
    summary = {
        "config": {
            "metal": METAL,
            "ligand": LIGAND,
            "contig": CONTIG,
            "num_designs": NUM_DESIGNS,
            "num_seqs_per_backbone": NUM_SEQS_PER_BACKBONE,
            "scout_ptm_threshold": SCOUT_PTM_THRESHOLD,
            "scout_plddt_threshold": SCOUT_PLDDT_THRESHOLD,
            "stability_ptm_min": STABILITY_PTM_MIN,
            "stability_plddt_min": STABILITY_PLDDT_MIN,
            "stability_ala_max": STABILITY_ALA_MAX,
        },
        "funnel": {
            "backbones_generated": result.backbones_generated,
            "backbones_after_scout": result.backbones_after_scout,
            "sequences_designed": result.sequences_designed,
            "sequences_after_stability": result.sequences_after_stability,
            "designs_analyzed": result.designs_analyzed,
            "designs_passed_metal_binding": result.designs_passed_metal_binding,
            "designs_passed_strict": result.designs_passed_strict,
        },
        "timings": result.timings,
        "errors": result.errors,
        "history_saved": result.history_saved,
        "history_run_id": result.history_run_id,
        "lessons_detected": result.lessons_detected,
        "designs": [
            {
                "seq_id": d.seq_id,
                "backbone_id": d.backbone_id,
                "sequence": d.sequence,
                "ptm": d.ptm,
                "plddt": d.plddt,
                "pae": d.pae,
                "alanine_pct": d.alanine_pct,
                "total_residues": d.total_residues,
                "stability_passed": d.stability_passed,
                "coordination_number": d.coordination_number,
                "geometry_rmsd": d.geometry_rmsd,
                "filter_preset": d.filter_preset,
                "filter_passed": d.filter_passed,
                "failed_filters": d.failed_filters,
                "ligand_contacts": d.ligand_contacts,
                "hsab_score": d.hsab_score,
                "composite_score": composite_score(d),
            }
            for d in result.designs
        ],
    }

    summary_path = os.path.join(OUTPUT_DIR, "pipeline_summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  Saved pipeline_summary.json -> {summary_path}")

    return summary_path


# ============== Main ==============

def main():
    print("=" * 72)
    print("  Ln-CITRATE PRODUCTION PIPELINE E2E TEST")
    print(f"  Metal={METAL}, Ligand={LIGAND}, Contig={CONTIG}")
    print(f"  Designs={NUM_DESIGNS}, Seqs/backbone={NUM_SEQS_PER_BACKBONE}")
    print("=" * 72)

    # Health check
    print("\nChecking API health...")
    health = call_api("health", {}, timeout=30)
    if health.get("status") not in ("completed", "ok"):
        print(f"API health check failed: {health}")
        print("Make sure Docker is running (use docker-wsl skill)")
        sys.exit(1)
    print("API is healthy.\n")

    result = PipelineResult()
    pipeline_t0 = time.time()

    # Run pipeline steps sequentially
    result = step1_rfd3_backbone_generation(result)
    result = step2_scout_filter(result)
    result = step3_mpnn_sequence_design(result)
    result = step4_rf3_validation(result)
    result = step5_stability_filter(result)
    result = step6_analyze_design(result)
    result = step7_save_history(result)
    result = step8_check_lessons(result)

    result.timings["pipeline_total"] = time.time() - pipeline_t0

    # Step 9: Report
    print(f"\n{'#'*72}")
    print(f"#  Step 9: Report")
    print(f"{'#'*72}")

    # Save outputs
    save_outputs(result)

    # Generate and display report
    report = format_report(result)
    report_path = os.path.join(OUTPUT_DIR, "report.txt")
    with open(report_path, "w") as f:
        f.write(report)
    print(f"  Saved report.txt -> {report_path}")

    print("\n" + report)


if __name__ == "__main__":
    main()
