#!/usr/bin/env python3
"""AF3 cross-validation for Tb-citrate design candidates.

Generates AF3 input JSONs for dual validation (Ca-citrate + Tb-citrate),
runs AF3 inference, and extracts confidence metrics.

Supports tiered validation:
  Phase 1 (screen):   all 173 designs, 1 seed, Ca+Tb only, noMSA
  Phase 2 (validate): hits from Phase 1, 3 seeds, Ca+Tb+apo, noMSA
  Phase 3 (final):    top candidates, 5 seeds, Ca+Tb+apo, noMSA

Usage:
    # Phase 1: Screen all 173 designs with 1 seed, noMSA
    python af3_validate.py prepare --all --seeds 1 --noMSA --variant ca_tb

    # Phase 1: Analyze screening results
    python af3_validate.py analyze --phase screening

    # Phase 2: Validate hits with 3 seeds
    python af3_validate.py prepare --from-phase1 --seeds 3 --noMSA

    # Phase 3: Final top candidates with 5 seeds
    python af3_validate.py prepare --top 10 --seeds 5 --noMSA

    # Run everything end-to-end
    python af3_validate.py run --all --seeds 1 --noMSA --variant ca_tb

    # Legacy: top N with full MSA
    python af3_validate.py prepare --top 10 --seeds 5
"""

import argparse
import csv
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional


# Paths
PROJ_ROOT = Path(__file__).resolve().parent
PROMISING_DIR = PROJ_ROOT / "final" / "promising"
ANALYSIS_DIR = PROJ_ROOT / "analysis"
AF3_RESULTS_DIR = PROJ_ROOT / "analysis" / "af3_validation"

# WSL paths — AF3 runs in WSL, so input/output use WSL home
WSL_HOME = "/home/polaire"
WSL_AF3_INPUT = f"{WSL_HOME}/af_input/tb_citrate_validation"
WSL_AF3_OUTPUT = f"{WSL_HOME}/af_output"

# Windows-accessible paths via UNC for reading/writing WSL files
WSL_UNC = "//wsl.localhost/Ubuntu"
AF3_INPUT_DIR = Path(f"{WSL_UNC}{WSL_HOME}/af_input/tb_citrate_validation")
AF3_OUTPUT_DIR = Path(f"{WSL_UNC}{WSL_HOME}/af_output")

# Thresholds
THRESHOLDS = {
    "ca_cit": {
        "iptm_min": 0.8,
        "ptm_min": 0.7,
        "plddt_min": 70,
        "note": "Strict - Ca well-represented in AF3 training (16K structures)",
    },
    "tb_cit": {
        "iptm_min": 0.6,
        "ptm_min": 0.5,
        "plddt_min": 60,
        "note": "Relaxed - Tb sparse in AF3 training (62 PDB structures)",
    },
}


def load_designs(top_n: Optional[int] = None) -> list[dict]:
    """Load passing designs by composite score from all batch JSONs.

    Args:
        top_n: If set, return only top N designs. If None, return all passing.
    """
    all_designs = []

    for batch_json in sorted(ANALYSIS_DIR.glob("production_batch_*.json")):
        with open(batch_json) as f:
            data = json.load(f)
        batch_num = data.get("batch", 0)
        for d in data.get("designs", []):
            if d.get("filter_passed"):
                d["batch"] = batch_num
                all_designs.append(d)

    # Sort by composite score descending
    all_designs.sort(key=lambda x: x.get("composite_score", 0), reverse=True)

    if top_n is not None:
        return all_designs[:top_n]
    return all_designs


def load_phase1_hits() -> list[dict]:
    """Load designs that passed Phase 1 screening.

    Uses phase2_candidates.json if available (curated selection),
    otherwise falls back to Ca>=0.8 OR Tb>=0.8 from phase1 results.
    """
    # Prefer curated candidate list
    candidates_path = AF3_RESULTS_DIR / "phase2_candidates.json"
    if candidates_path.exists():
        with open(candidates_path) as f:
            data = json.load(f)
        hits = data.get("designs", [])
        if hits:
            print(f"  Loaded {len(hits)} Phase 2 candidates from {candidates_path.name}")
            return hits

    phase1_path = AF3_RESULTS_DIR / "phase1_screening.json"
    if not phase1_path.exists():
        print(f"No Phase 1 results at {phase1_path}")
        print("Run Phase 1 screening first.")
        return []

    with open(phase1_path) as f:
        phase1 = json.load(f)

    # Strict selection: Ca iPTM >= 0.8 OR Tb iPTM >= 0.8
    hits = []
    for d in phase1.get("designs", []):
        ca_iptm = d.get("ca_cit_iptm") or 0
        tb_iptm = d.get("tb_cit_iptm") or 0
        if ca_iptm >= 0.8 or tb_iptm >= 0.8:
            hits.append(d)

    if not hits:
        # Fallback: relax to iPTM >= 0.6
        print("  No hits at strict thresholds. Relaxing to iPTM >= 0.6...")
        for d in phase1.get("designs", []):
            ca_iptm = d.get("ca_cit_iptm") or 0
            tb_iptm = d.get("tb_cit_iptm") or 0
            if ca_iptm >= 0.6 or tb_iptm >= 0.6:
                hits.append(d)

    return hits


def load_phase2_hits() -> list[dict]:
    """Load designs that passed Phase 2 validation."""
    phase2_path = AF3_RESULTS_DIR / "phase2_validation.json"
    if not phase2_path.exists():
        print(f"No Phase 2 results at {phase2_path}")
        return []

    with open(phase2_path) as f:
        phase2 = json.load(f)

    return [d for d in phase2.get("designs", []) if d.get("phase2_pass")]


def make_af3_input(
    name: str,
    sequence: str,
    metal_ccd: str = "CA",
    ligand_ccd: str = "CIT",
    seeds: list[int] = None,
    no_msa: bool = False,
) -> dict:
    """Create AF3 input JSON for a design with metal + ligand.

    Args:
        name: Job name
        sequence: Protein sequence
        metal_ccd: Metal CCD code ("CA" for calcium, "TB" for terbium)
        ligand_ccd: Ligand CCD code ("CIT" for citric acid)
        seeds: Random seeds for sampling
        no_msa: If True, add unpairedMsa field for --run_data_pipeline=false mode
    """
    if seeds is None:
        seeds = list(range(1, 6))

    protein_entry = {
        "id": "A",
        "sequence": sequence,
    }
    if no_msa:
        # AF3 requires all three fields when --run_data_pipeline=false:
        # unpairedMsa, pairedMsa, and templates
        a3m_query = f">query\n{sequence}\n"
        protein_entry["unpairedMsa"] = a3m_query
        protein_entry["pairedMsa"] = a3m_query
        protein_entry["templates"] = []

    return {
        "name": name,
        "sequences": [
            {"protein": protein_entry},
            {"ligand": {"id": "B", "ccdCodes": [metal_ccd]}},
            {"ligand": {"id": "C", "ccdCodes": [ligand_ccd]}},
        ],
        "modelSeeds": seeds,
        "dialect": "alphafold3",
        "version": 2,
    }


def make_apo_input(
    name: str,
    sequence: str,
    seeds: list[int] = None,
    no_msa: bool = False,
) -> dict:
    """Create AF3 input JSON for protein-only (apo) control."""
    if seeds is None:
        seeds = list(range(1, 6))

    protein_entry = {
        "id": "A",
        "sequence": sequence,
    }
    if no_msa:
        a3m_query = f">query\n{sequence}\n"
        protein_entry["unpairedMsa"] = a3m_query
        protein_entry["pairedMsa"] = a3m_query
        protein_entry["templates"] = []

    return {
        "name": name,
        "sequences": [{"protein": protein_entry}],
        "modelSeeds": seeds,
        "dialect": "alphafold3",
        "version": 2,
    }


def prepare_inputs(
    top_n: Optional[int] = None,
    seeds: int = 5,
    no_msa: bool = False,
    variant: str = "all",
    from_phase1: bool = False,
    from_phase2: bool = False,
    phase: str = "screening",
):
    """Generate AF3 input JSONs for designs.

    Args:
        top_n: Number of top designs (None = all passing)
        seeds: Number of seeds per job
        no_msa: Use noMSA mode (--run_data_pipeline=false)
        variant: "all" (Ca+Tb+apo), "ca_tb" (Ca+Tb only), "ca_cit", "tb_cit", "apo"
        from_phase1: Load hits from Phase 1 results
        from_phase2: Load hits from Phase 2 results
        phase: Phase label for manifest ("screening", "validation", "final")
    """
    # Load designs based on source
    if from_phase1:
        hits = load_phase1_hits()
        if not hits:
            return
        # Reconstruct full design info from batch JSONs
        all_designs = load_designs()
        design_map = {f"b{d['batch']:02d}_{d['seq_id']}": d for d in all_designs}
        designs = []
        for h in hits:
            key = h.get("design", "")
            if key in design_map:
                designs.append(design_map[key])
        if not designs:
            print("Could not match Phase 1 hits to batch designs")
            return
        print(f"\nLoaded {len(designs)} Phase 1 hits for Phase 2 validation")
    elif from_phase2:
        hits = load_phase2_hits()
        if not hits:
            return
        all_designs = load_designs()
        design_map = {f"b{d['batch']:02d}_{d['seq_id']}": d for d in all_designs}
        designs = []
        for h in hits:
            key = h.get("design", "")
            if key in design_map:
                designs.append(design_map[key])
        if not designs:
            print("Could not match Phase 2 hits to batch designs")
            return
        if top_n:
            designs = designs[:top_n]
        print(f"\nLoaded {len(designs)} Phase 2 hits for Phase 3 final")
    else:
        designs = load_designs(top_n)

    if not designs:
        print("No passing designs found!")
        return

    # Determine which variants to generate
    if variant == "ca_tb":
        variant_list = ["ca_cit", "tb_cit"]
    elif variant == "all":
        variant_list = ["ca_cit", "tb_cit", "apo"]
    else:
        variant_list = [variant]

    # Use phase-specific subdirectory
    phase_input_dir = AF3_INPUT_DIR / phase
    phase_input_dir.mkdir(parents=True, exist_ok=True)

    seed_list = list(range(1, seeds + 1))

    print(f"\nPreparing AF3 inputs — Phase: {phase}")
    print(f"  Designs: {len(designs)}")
    print(f"  Variants: {variant_list}")
    print(f"  Seeds: {seed_list}")
    print(f"  noMSA: {no_msa}")
    print(f"  Output: {phase_input_dir}\n")

    for i, design in enumerate(designs):
        seq_id = design["seq_id"]
        batch = design["batch"]
        sequence = design["sequence"]
        score = design.get("composite_score", 0)
        prefix = f"b{batch:02d}_{seq_id}"

        print(f"  {i+1:>3}. {prefix} (score={score:.3f}, pTM={design['ptm']:.3f}, "
              f"pLDDT={design['plddt']:.3f}, CN={design['coordination_number']})")

        for var in variant_list:
            if var == "apo":
                inp = make_apo_input(
                    name=f"{prefix}_apo",
                    sequence=sequence,
                    seeds=seed_list,
                    no_msa=no_msa,
                )
            else:
                metal = "CA" if var == "ca_cit" else "TB"
                inp = make_af3_input(
                    name=f"{prefix}_{var}",
                    sequence=sequence,
                    metal_ccd=metal,
                    ligand_ccd="CIT",
                    seeds=seed_list,
                    no_msa=no_msa,
                )
            out_path = phase_input_dir / f"{prefix}_{var}.json"
            with open(out_path, "w") as f:
                json.dump(inp, f, indent=2)

    # Write batch runner script (runs in WSL — use WSL paths)
    runner_path = phase_input_dir / "run_all.sh"
    af3_extra_flags = " --run_data_pipeline=false" if no_msa else ""
    wsl_input_dir = f"{WSL_AF3_INPUT}/{phase}"

    with open(runner_path, "w", newline="\n") as f:
        f.write("#!/bin/bash\n")
        f.write(f"# AF3 batch runner — Phase: {phase}\n")
        f.write(f"# {len(designs)} designs x {len(variant_list)} variants x {seeds} seeds\n")
        f.write(f"# noMSA: {no_msa}\n\n")
        f.write("set -e\n\n")
        f.write("TOTAL=0\nPASSED=0\nFAILED=0\n")
        f.write(f'START_TIME=$(date +%s)\n\n')

        for i, design in enumerate(designs):
            seq_id = design["seq_id"]
            batch = design["batch"]
            prefix = f"b{batch:02d}_{seq_id}"
            for var in variant_list:
                job_name = f"{prefix}_{var}"
                wsl_input_file = f"{wsl_input_dir}/{job_name}.json"
                # Check if output already exists (skip completed)
                f.write(f'# [{i+1}/{len(designs)}] {job_name}\n')
                f.write(f'if [ -d "$HOME/af_output/{job_name}" ] && ls "$HOME/af_output/{job_name}"/*_model.cif >/dev/null 2>&1; then\n')
                f.write(f'    echo "[SKIP] {job_name} (done)"\n')
                f.write(f'else\n')
                f.write(f'    echo "[RUN] {job_name} ({i+1}/{len(designs)})"\n')
                f.write(f'    if ~/run_alphafold3.sh "{wsl_input_file}"{af3_extra_flags}; then\n')
                f.write(f'        echo "  OK"\n')
                f.write(f'        PASSED=$((PASSED+1))\n')
                f.write(f'    else\n')
                f.write(f'        echo "  FAILED"\n')
                f.write(f'        FAILED=$((FAILED+1))\n')
                f.write(f'    fi\n')
                f.write(f'fi\n')
                f.write(f'TOTAL=$((TOTAL+1))\n\n')

        f.write('END_TIME=$(date +%s)\n')
        f.write('ELAPSED=$((END_TIME - START_TIME))\n')
        f.write('echo ""\n')
        f.write('echo "========================================"\n')
        f.write(f'echo "  AF3 {phase.title()} Complete"\n')
        f.write('echo "  Total: $TOTAL, Passed: $PASSED, Failed: $FAILED"\n')
        f.write('echo "  Time: ${ELAPSED}s"\n')
        f.write('echo "========================================"\n')

    os.chmod(runner_path, 0o755)

    total_jobs = len(designs) * len(variant_list)
    print(f"\n  Generated {total_jobs} AF3 input JSONs")
    print(f"  Batch runner: {runner_path}")
    print(f"\n  To run: bash {runner_path}")

    # Save manifest
    manifest = {
        "phase": phase,
        "no_msa": no_msa,
        "designs": [
            {
                "rank": i + 1,
                "batch": d["batch"],
                "seq_id": d["seq_id"],
                "sequence": d["sequence"],
                "rf3_ptm": d["ptm"],
                "rf3_plddt": d["plddt"],
                "rf3_cn": d["coordination_number"],
                "rf3_lc": d.get("ligand_contacts"),
                "composite_score": d.get("composite_score", 0),
            }
            for i, d in enumerate(designs)
        ],
        "af3_seeds": seed_list,
        "variants": variant_list,
        "thresholds": THRESHOLDS,
    }
    manifest_path = phase_input_dir / "manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    # Also save/overwrite top-level manifest for analyze compatibility
    top_manifest_path = AF3_INPUT_DIR / "manifest.json"
    with open(top_manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"  Manifest: {manifest_path}")


def parse_af3_output(output_dir: Path, job_name: str) -> Optional[dict]:
    """Parse AF3 output for a single job, extracting confidence metrics."""
    job_dir = output_dir / job_name
    if not job_dir.exists():
        return None

    results = []

    # AF3 outputs multiple seeds as separate mmCIF files
    for cif_file in sorted(job_dir.glob("*_model.cif")):
        summary_file = cif_file.with_name(
            cif_file.name.replace("_model.cif", "_summary_confidences.json")
        )
        if summary_file.exists():
            with open(summary_file) as f:
                conf = json.load(f)
            results.append({
                "seed_file": cif_file.name,
                "ptm": conf.get("ptm", None),
                "iptm": conf.get("iptm", None),
                "ranking_score": conf.get("ranking_score", None),
                "fraction_disordered": conf.get("fraction_disordered", None),
                "has_clash": conf.get("has_clash", None),
            })

    if not results:
        # Try alternative output format
        for conf_file in sorted(job_dir.glob("*confidence*.json")):
            with open(conf_file) as f:
                conf = json.load(f)
            results.append({
                "file": conf_file.name,
                **{k: conf.get(k) for k in ["ptm", "iptm", "ranking_score",
                                              "fraction_disordered", "has_clash"]},
            })

    if not results:
        return None

    best = max(results, key=lambda x: x.get("ranking_score") or x.get("ptm") or 0)

    # Compute stats across seeds
    iptms = [r["iptm"] for r in results if r.get("iptm") is not None]
    ptms = [r["ptm"] for r in results if r.get("ptm") is not None]

    return {
        "job_name": job_name,
        "n_models": len(results),
        "best_ptm": best.get("ptm"),
        "best_iptm": best.get("iptm"),
        "best_ranking_score": best.get("ranking_score"),
        "mean_iptm": sum(iptms) / len(iptms) if iptms else None,
        "std_iptm": (sum((x - sum(iptms)/len(iptms))**2 for x in iptms) / len(iptms))**0.5 if len(iptms) > 1 else 0.0,
        "mean_ptm": sum(ptms) / len(ptms) if ptms else None,
        "all_models": results,
    }


def analyze_results(phase: str = "screening", save_as: Optional[str] = None):
    """Analyze AF3 outputs and compare to RF3 metrics.

    Args:
        phase: Which phase to analyze ("screening", "validation", "final")
        save_as: Override output filename (default: phase-based)
    """
    # Try phase-specific manifest (in meta/ subdir or directly), fall back to top-level
    manifest_path = AF3_INPUT_DIR / phase / "meta" / "manifest.json"
    if not manifest_path.exists():
        manifest_path = AF3_INPUT_DIR / phase / "manifest.json"
    if not manifest_path.exists():
        manifest_path = AF3_INPUT_DIR / "manifest.json"
    if not manifest_path.exists():
        print(f"No manifest found. Run 'prepare' first.")
        return

    with open(manifest_path) as f:
        manifest = json.load(f)

    AF3_RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    variants = manifest.get("variants", ["ca_cit", "tb_cit", "apo"])
    ca_thresholds = manifest["thresholds"]["ca_cit"]
    tb_thresholds = manifest["thresholds"]["tb_cit"]

    print("\n" + "=" * 70)
    print(f"  AF3 Cross-Validation Results — Phase: {phase}")
    print("=" * 70)

    results = []
    for design in manifest["designs"]:
        rank = design["rank"]
        batch = design["batch"]
        seq_id = design["seq_id"]
        prefix = f"b{batch:02d}_{seq_id}"

        row = {
            "rank": rank,
            "design": prefix,
            "sequence": design["sequence"],
            "rf3_ptm": design["rf3_ptm"],
            "rf3_plddt": design["rf3_plddt"],
            "rf3_cn": design["rf3_cn"],
            "composite_score": design.get("composite_score", 0),
        }

        for var in ["ca_cit", "tb_cit", "apo"]:
            job_name = f"{prefix}_{var}"
            af3 = parse_af3_output(AF3_OUTPUT_DIR, job_name)
            if af3:
                row[f"{var}_ptm"] = af3["best_ptm"]
                row[f"{var}_iptm"] = af3["best_iptm"]
                row[f"{var}_ranking"] = af3["best_ranking_score"]
                row[f"{var}_models"] = af3["n_models"]
                row[f"{var}_mean_iptm"] = af3.get("mean_iptm")
                row[f"{var}_std_iptm"] = af3.get("std_iptm")
            else:
                for suffix in ["_ptm", "_iptm", "_ranking", "_models", "_mean_iptm", "_std_iptm"]:
                    row[f"{var}{suffix}"] = None if suffix != "_models" else 0

        # Evaluate pass/fail
        ca_pass = (
            row.get("ca_cit_iptm") is not None
            and row["ca_cit_iptm"] >= ca_thresholds["iptm_min"]
            and (row.get("ca_cit_ptm") or 0) >= ca_thresholds["ptm_min"]
        )
        tb_pass = (
            row.get("tb_cit_iptm") is not None
            and row["tb_cit_iptm"] >= tb_thresholds["iptm_min"]
            and (row.get("tb_cit_ptm") or 0) >= tb_thresholds["ptm_min"]
        )
        row["ca_pass"] = ca_pass
        row["tb_pass"] = tb_pass

        # Confidence delta: holo vs apo
        if row.get("ca_cit_ptm") and row.get("apo_ptm"):
            row["ca_delta_ptm"] = row["ca_cit_ptm"] - row["apo_ptm"]
        else:
            row["ca_delta_ptm"] = None

        # Phase 2 pass: consistent across seeds + passes threshold
        if phase == "validation" and row.get("ca_cit_std_iptm") is not None:
            row["phase2_pass"] = (
                (ca_pass or tb_pass)
                and (row.get("ca_cit_std_iptm") or 0) < 0.05
            )
        else:
            row["phase2_pass"] = ca_pass or tb_pass

        results.append(row)

    # Print results table
    has_apo = "apo" in variants
    print(f"\n{'Rank':>4} {'Design':<18} {'RF3':>6} {'RF3':>6} {'CN':>3} "
          f"{'Ca iPTM':>8} {'Ca pTM':>7} {'Ca':>3} "
          f"{'Tb iPTM':>8} {'Tb pTM':>7} {'Tb':>3}", end="")
    if has_apo:
        print(f" {'Apo pTM':>8} {'dPTM':>6}", end="")
    print()

    print(f"{'':>4} {'':>18} {'pTM':>6} {'pLDDT':>6} {'':>3} "
          f"{'':>8} {'':>7} {'':>3} "
          f"{'':>8} {'':>7} {'':>3}", end="")
    if has_apo:
        print(f" {'':>8} {'Ca-Apo':>6}", end="")
    print()
    print("-" * (110 if has_apo else 94))

    ca_passed = tb_passed = both_passed = 0

    for r in results:
        ca_flag = "Y" if r["ca_pass"] else "N" if r.get("ca_cit_iptm") is not None else "?"
        tb_flag = "Y" if r["tb_pass"] else "N" if r.get("tb_cit_iptm") is not None else "?"

        def fmt(v, w=8):
            return f"{v:.3f}" if v is not None else " " * (w - 3) + "---"

        line = (f"{r['rank']:>4} {r['design']:<18} {r['rf3_ptm']:.3f} {r['rf3_plddt']:.3f} "
                f"{r['rf3_cn']:>3} "
                f"{fmt(r.get('ca_cit_iptm')):>8} {fmt(r.get('ca_cit_ptm'), 7):>7} {ca_flag:>3} "
                f"{fmt(r.get('tb_cit_iptm')):>8} {fmt(r.get('tb_cit_ptm'), 7):>7} {tb_flag:>3}")
        if has_apo:
            delta = f"{r['ca_delta_ptm']:+.3f}" if r.get("ca_delta_ptm") is not None else "   ---"
            line += f" {fmt(r.get('apo_ptm')):>8} {delta:>6}"
        print(line)

        if r["ca_pass"]:
            ca_passed += 1
        if r["tb_pass"]:
            tb_passed += 1
        if r["ca_pass"] and r["tb_pass"]:
            both_passed += 1

    print("-" * (110 if has_apo else 94))
    n = len(results)
    n_run = sum(1 for r in results if r.get("ca_cit_iptm") is not None)
    print(f"\nCompleted: {n_run}/{n} designs")
    print(f"Ca-citrate pass (iPTM>={ca_thresholds['iptm_min']}): {ca_passed}/{n_run}")
    print(f"Tb-citrate pass (iPTM>={tb_thresholds['iptm_min']}): {tb_passed}/{n_run}")
    print(f"Both pass: {both_passed}/{n_run}")

    # Interpretation
    print("\nInterpretation:")
    if n_run == 0:
        print("  No AF3 results found yet. Run AF3 first.")
    elif ca_passed > 0 and tb_passed == 0:
        print("  Ca passes but Tb fails -> Scaffolds are sound, AF3 lacks Tb training data")
    elif ca_passed > 0 and tb_passed > 0:
        print("  Both pass -> Strong candidates for experimental testing!")
    elif ca_passed == 0:
        print("  Ca fails -> Scaffold design may need revision, or relax thresholds")

    # Confidence delta analysis
    deltas = [r["ca_delta_ptm"] for r in results if r.get("ca_delta_ptm") is not None]
    if deltas:
        avg_delta = sum(deltas) / len(deltas)
        print(f"\n  Avg pTM delta (Ca-holo vs apo): {avg_delta:+.3f}")
        if avg_delta > 0.05:
            print("  Positive delta = AF3 recognizes ligand binding -> good sign")
        elif avg_delta < -0.05:
            print("  Negative delta = ligand may destabilize prediction")
        else:
            print("  Near-zero delta = ligand doesn't significantly affect confidence")

    # RF3 vs AF3 correlation
    pairs = [(r["rf3_ptm"], r.get("ca_cit_iptm"))
             for r in results if r.get("ca_cit_iptm") is not None]
    if len(pairs) >= 5:
        rf3_vals = [p[0] for p in pairs]
        af3_vals = [p[1] for p in pairs]
        mean_rf3 = sum(rf3_vals) / len(rf3_vals)
        mean_af3 = sum(af3_vals) / len(af3_vals)
        cov = sum((r - mean_rf3) * (a - mean_af3) for r, a in zip(rf3_vals, af3_vals)) / len(pairs)
        std_rf3 = (sum((r - mean_rf3)**2 for r in rf3_vals) / len(rf3_vals))**0.5
        std_af3 = (sum((a - mean_af3)**2 for a in af3_vals) / len(af3_vals))**0.5
        if std_rf3 > 0 and std_af3 > 0:
            corr = cov / (std_rf3 * std_af3)
            print(f"\n  RF3 pTM vs AF3 Ca iPTM correlation: r={corr:.3f} (n={len(pairs)})")

    # Save results
    output_name = save_as or f"phase{'1' if phase == 'screening' else '2' if phase == 'validation' else '3'}_{phase}.json"
    results_path = AF3_RESULTS_DIR / output_name
    # Strip sequences from saved results to keep file size manageable
    save_results = []
    for r in results:
        sr = {k: v for k, v in r.items() if k != "sequence"}
        save_results.append(sr)

    with open(results_path, "w") as f:
        json.dump({
            "phase": phase,
            "designs": save_results,
            "summary": {
                "total": n,
                "completed": n_run,
                "ca_passed": ca_passed,
                "tb_passed": tb_passed,
                "both_passed": both_passed,
                "ca_thresholds": ca_thresholds,
                "tb_thresholds": tb_thresholds,
            },
        }, f, indent=2)
    print(f"\n  Results saved: {results_path}")

    # Also save as the canonical phase file
    phase_map = {"screening": "phase1_screening.json",
                 "validation": "phase2_validation.json",
                 "final": "phase3_final.json"}
    if phase in phase_map:
        canonical = AF3_RESULTS_DIR / phase_map[phase]
        with open(canonical, "w") as f:
            json.dump({
                "phase": phase,
                "designs": save_results,
                "summary": {
                    "total": n,
                    "completed": n_run,
                    "ca_passed": ca_passed,
                    "tb_passed": tb_passed,
                    "both_passed": both_passed,
                    "ca_thresholds": ca_thresholds,
                    "tb_thresholds": tb_thresholds,
                },
            }, f, indent=2)

    return results


def run_af3(
    top_n: Optional[int] = None,
    seeds: int = 5,
    variant: str = "all",
    no_msa: bool = False,
    from_phase1: bool = False,
    from_phase2: bool = False,
    phase: str = "screening",
):
    """Run full AF3 validation: prepare + execute + analyze."""
    # Check if Docker container is using GPU
    try:
        result = subprocess.run(
            ["docker", "ps", "--format", "{{.Names}}"],
            capture_output=True, text=True, timeout=5,
        )
        if "serverless-rfdiffusion" in result.stdout:
            print("WARNING: RFD3 Docker container is running!")
            print("Stop it first: docker compose -f docker-compose.local.yml down")
            print("AF3 needs the GPU exclusively.")
            return
    except Exception:
        pass

    prepare_inputs(top_n, seeds, no_msa, variant, from_phase1, from_phase2, phase)

    print("\n" + "=" * 50)
    print(f"  Running AF3 Inference — Phase: {phase}")
    print("=" * 50)

    phase_input_dir = AF3_INPUT_DIR / phase
    manifest_path = phase_input_dir / "manifest.json"
    if not manifest_path.exists():
        manifest_path = AF3_INPUT_DIR / "manifest.json"

    with open(manifest_path) as f:
        manifest = json.load(f)

    variants = manifest.get("variants", ["ca_cit", "tb_cit", "apo"])
    af3_extra_args = ["--run_data_pipeline=false"] if no_msa else []

    total = 0
    completed = 0
    failed = 0
    start_time = time.time()

    for i, design in enumerate(manifest["designs"]):
        batch = design["batch"]
        seq_id = design["seq_id"]
        prefix = f"b{batch:02d}_{seq_id}"

        for var in variants:
            job_name = f"{prefix}_{var}"
            input_path = phase_input_dir / f"{job_name}.json"
            total += 1

            # Skip if already completed
            output_check = AF3_OUTPUT_DIR / job_name
            if output_check.exists() and any(output_check.glob("*_model.cif")):
                print(f"  [SKIP] {job_name} (already completed)")
                completed += 1
                continue

            print(f"\n  [RUN] {job_name} ({i+1}/{len(manifest['designs'])})...")
            try:
                cmd = [str(Path.home() / "run_alphafold3.sh"), str(input_path)] + af3_extra_args
                result = subprocess.run(
                    cmd,
                    capture_output=True, text=True, timeout=3600,
                )
                if result.returncode == 0:
                    print(f"  [OK] {job_name}")
                    completed += 1
                else:
                    print(f"  [FAIL] {job_name}: {result.stderr[-200:]}")
                    failed += 1
            except subprocess.TimeoutExpired:
                print(f"  [TIMEOUT] {job_name}")
                failed += 1
            except Exception as e:
                print(f"  [ERROR] {job_name}: {e}")
                failed += 1

    elapsed = time.time() - start_time
    print(f"\n  AF3 {phase}: {completed}/{total} completed, {failed} failed ({elapsed:.0f}s)")

    analyze_results(phase)


def main():
    parser = argparse.ArgumentParser(
        description="AF3 cross-validation for Tb-citrate designs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Phase 1: Screen all 173 with 1 seed, noMSA
  python af3_validate.py prepare --all --seeds 1 --noMSA --variant ca_tb
  python af3_validate.py run --all --seeds 1 --noMSA --variant ca_tb

  # Phase 2: Validate Phase 1 hits with 3 seeds
  python af3_validate.py run --from-phase1 --seeds 3 --noMSA

  # Phase 3: Final top 10 with 5 seeds
  python af3_validate.py run --top 10 --seeds 5 --noMSA --phase final

  # Analyze results
  python af3_validate.py analyze --phase screening
  python af3_validate.py analyze --phase validation
  python af3_validate.py analyze --phase final
""",
    )
    sub = parser.add_subparsers(dest="command")

    # Common args helper
    def add_common_args(p):
        p.add_argument("--top", type=int, default=None, help="Number of top designs (default: 10 unless --all)")
        p.add_argument("--all", action="store_true", help="Use all 173 passing designs")
        p.add_argument("--seeds", type=int, default=5, help="Number of AF3 seeds")
        p.add_argument("--noMSA", action="store_true", help="Skip MSA search (fast mode)")
        p.add_argument("--variant", default="all",
                       choices=["all", "ca_tb", "ca_cit", "tb_cit", "apo"],
                       help="Which variants to generate")
        p.add_argument("--from-phase1", action="store_true", help="Load hits from Phase 1")
        p.add_argument("--from-phase2", action="store_true", help="Load hits from Phase 2")
        p.add_argument("--phase", default="screening",
                       choices=["screening", "validation", "final"],
                       help="Phase label")

    # prepare
    p_prep = sub.add_parser("prepare", help="Generate AF3 input JSONs")
    add_common_args(p_prep)

    # analyze
    p_analyze = sub.add_parser("analyze", help="Parse AF3 outputs and compare to RF3")
    p_analyze.add_argument("--phase", default="screening",
                           choices=["screening", "validation", "final"],
                           help="Phase to analyze")

    # run
    p_run = sub.add_parser("run", help="Full pipeline: prepare + run AF3 + analyze")
    add_common_args(p_run)

    args = parser.parse_args()

    # Default top_n logic: 10 unless --all
    if hasattr(args, "top") and hasattr(args, "all"):
        if args.all:
            top_n = None  # all designs
        elif args.top is not None:
            top_n = args.top
        else:
            top_n = 10  # default
    else:
        top_n = 10

    if args.command == "prepare":
        prepare_inputs(
            top_n=top_n,
            seeds=args.seeds,
            no_msa=args.noMSA,
            variant=args.variant,
            from_phase1=args.from_phase1,
            from_phase2=args.from_phase2,
            phase=args.phase,
        )
    elif args.command == "analyze":
        analyze_results(phase=args.phase)
    elif args.command == "run":
        run_af3(
            top_n=top_n,
            seeds=args.seeds,
            variant=args.variant,
            no_msa=args.noMSA,
            from_phase1=args.from_phase1,
            from_phase2=args.from_phase2,
            phase=args.phase,
        )
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
