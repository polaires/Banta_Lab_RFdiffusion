"""
Production pipeline for enzyme design: RFD3 + LigandMPNN + RF3 validation.

End-to-end automated pipeline with:
- Configurable parameters via CLI
- 8-stage pipeline with checkpointing
- Optimized LigandMPNN settings for large cofactors
- Both ESMFold and RF3 validation
- Comprehensive output structure

Usage:
    python production_run.py \
        --pdb-id 4CVB --ligand PQQ --metal CA \
        --num-backbones 50 --num-seqs 8 \
        --output-dir ./output_production_01

    # Resume from checkpoint:
    python production_run.py --resume ./output_production_01

    # Parameter sweep mode:
    python production_run.py --sweep-config tier1_sweep.json --output-dir ./sweep_tier1
"""
import argparse
import asyncio
import json
import math
import os
import sys
import time
from datetime import datetime
from pathlib import Path

# PQQ canonical SMILES (PubChem CID 1024)
PQQ_SMILES = "C1=C(C2=C(C(=O)C(=O)C3=C2NC(=C3)C(=O)O)N=C1C(=O)O)C(=O)O"

# Default API URL
BASE_URL = os.getenv("DESIGN_API_URL", "http://localhost:8000/runsync")


###############################################################################
# CLI
###############################################################################

def build_parser():
    parser = argparse.ArgumentParser(
        description="Production enzyme design pipeline (RFD3 + LigandMPNN + RF3 validation)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Required
    parser.add_argument("--pdb-id", default="4CVB", help="Source PDB ID")
    parser.add_argument("--ligand", default="PQQ", help="Ligand 3-letter code")
    parser.add_argument("--metal", default="CA", help="Metal element symbol")

    # Pipeline options
    parser.add_argument("--enzyme-class", default="oxidoreductase", help="Enzyme class for motif selection")
    parser.add_argument("--num-backbones", type=int, default=50, help="Number of RFD3 backbone designs")
    parser.add_argument("--num-seqs", type=int, default=8, help="Sequences per backbone (LigandMPNN)")
    parser.add_argument("--scout-mode", action="store_true", default=True, help="Enable scout validation (generate 1 seq first, skip bad backbones)")
    parser.add_argument("--no-scout", dest="scout_mode", action="store_false", help="Disable scout validation")
    parser.add_argument("--scout-ptm-threshold", type=float, default=0.6, help="Minimum pTM for scout sequence to pass")
    parser.add_argument("--run-name", default=None, help="Run name for output directory naming")

    # RFD3 parameters
    parser.add_argument("--step-scale", type=float, default=1.5, help="RFD3 step scale")
    parser.add_argument("--cfg-scale", type=float, default=2.0, help="RFD3 classifier-free guidance scale")
    parser.add_argument("--scaffold-size", default=None, help="Scaffold size range (e.g. '150-200'). Auto-calculated if not set.")
    parser.add_argument("--num-timesteps", type=int, default=200, help="RFD3 diffusion timesteps")
    parser.add_argument("--gamma-0", type=float, default=None, help="RFD3 gamma_0 (None=default)")
    parser.add_argument("--enable-hbond", action="store_true", default=True, help="Enable H-bond conditioning for ligand")
    parser.add_argument("--no-hbond", dest="enable_hbond", action="store_false", help="Disable H-bond conditioning")

    # LigandMPNN parameters
    parser.add_argument("--mpnn-temperature", type=float, default=0.1, help="LigandMPNN sampling temperature")
    parser.add_argument("--mpnn-noise-level", default="010", help="LigandMPNN model noise level (005, 010, 020)")
    parser.add_argument("--mpnn-ligand-cutoff", type=float, default=10.0, help="LigandMPNN ligand cutoff for scoring (A)")
    parser.add_argument("--mpnn-num-packs", type=int, default=6, help="Number of packing iterations per design")
    parser.add_argument("--mpnn-repack-everything", action="store_true", default=True, help="Repack all residues")

    # Validation
    parser.add_argument("--validate-rf3", action="store_true", default=True, help="Run RF3 validation")
    parser.add_argument("--no-validate-rf3", dest="validate_rf3", action="store_false")
    parser.add_argument("--ligand-smiles", default=None, help="SMILES for ligand-aware RF3 (auto-detected for PQQ)")

    # Output
    parser.add_argument("--output-dir", default=None, help="Output directory (auto-generated if not set)")
    parser.add_argument("--resume", default=None, help="Resume from checkpoint directory")

    # Sweep mode
    parser.add_argument("--sweep-config", default=None, help="JSON file with parameter sweep configuration")

    return parser


###############################################################################
# CHECKPOINT HELPERS
###############################################################################

STAGES = [
    "scaffold",           # 1: Scaffold extraction
    "backbones",          # 2: RFD3 backbone generation
    "clash_filter",       # 3: Clash filtering
    "continuity",         # 4: Continuity check
    "scout_mpnn",         # 5a: Scout sequence generation (1 per backbone)
    "scout_validation",   # 5b: Scout RF3 validation (filter bad backbones)
    "mpnn",               # 6: LigandMPNN sequence design (remaining seqs for passing backbones)
    "rf3_validation",     # 7: RF3 validation
    "filtering",          # 8: Filtering + ranking
    "summary",            # 9: Summary generation
]


def save_checkpoint(output_dir, stage, data):
    """Save checkpoint for a pipeline stage."""
    cp_file = Path(output_dir) / f"checkpoint_{stage}.json"
    with open(cp_file, "w") as f:
        json.dump(data, f, indent=2, default=str)


def load_checkpoint(output_dir, stage):
    """Load checkpoint for a pipeline stage. Returns None if not found."""
    cp_file = Path(output_dir) / f"checkpoint_{stage}.json"
    if cp_file.exists():
        with open(cp_file) as f:
            return json.load(f)
    return None


def get_resume_stage(output_dir):
    """Find the last completed stage for resuming."""
    for stage in reversed(STAGES):
        if load_checkpoint(output_dir, stage) is not None:
            return stage
    return None


###############################################################################
# PIPELINE STAGES
###############################################################################

async def stage_scaffold(args, output_dir):
    """Stage 1: Extract scaffold from PDB."""
    print("\n" + "=" * 70)
    print("STAGE 1: Scaffold Extraction")
    print("=" * 70)

    from scaffolding_workflow import ScaffoldingWorkflow

    workflow = ScaffoldingWorkflow()
    scaffold_result = await workflow.run(
        pdb_id=args.pdb_id,
        ligand_code=args.ligand,
        metal=args.metal,
        chain_length=args.scaffold_size,
        include_all_ligand_contacts=True,
        fixed_atom_type="ALL",
        use_minimal_motif=True,
        enzyme_class=args.enzyme_class,
    )

    if not scaffold_result or not scaffold_result.success:
        print(f"[FAIL] Scaffold extraction failed: {getattr(scaffold_result, 'error_message', 'unknown')}")
        return None

    print(f"  Design approach: {scaffold_result.source_info.get('design_approach')}")
    print(f"  Catalytic residues: {len(scaffold_result.source_info.get('catalytic_residue_ids', []))}")
    print(f"  Fixed atoms: {len(scaffold_result.fixed_atoms)} entries")
    print(f"  Length: {scaffold_result.length}")
    print(f"  Unindex: {scaffold_result.unindex}")
    print(f"  Ligand codes: {scaffold_result.ligand_codes}")
    print(f"  H-bond acceptors: {scaffold_result.hbond_acceptors}")
    print(f"  H-bond donors: {scaffold_result.hbond_donors}")

    # Save checkpoint (serialize the result)
    checkpoint = {
        "motif_pdb": scaffold_result.motif_pdb,
        "length": scaffold_result.length,
        "unindex": scaffold_result.unindex,
        "ligand_codes": scaffold_result.ligand_codes,
        "fixed_atoms": scaffold_result.fixed_atoms,
        "rasa_targets": scaffold_result.rasa_targets,
        "hbond_acceptors": scaffold_result.hbond_acceptors,
        "hbond_donors": scaffold_result.hbond_donors,
        "source_info": scaffold_result.source_info,
        "catalytic_residue_ids": scaffold_result.source_info.get("catalytic_residue_ids", []),
    }
    save_checkpoint(output_dir, "scaffold", checkpoint)

    return checkpoint


def stage_backbones(args, output_dir, scaffold_data):
    """Stage 2: Generate backbone designs with RFD3."""
    print("\n" + "=" * 70)
    print("STAGE 2: RFD3 Backbone Generation")
    print("=" * 70)

    from api_client import run_rfd3_api
    from run_50_designs import _merge_hetatm_to_backbone

    rfd3_params = {
        "length": args.scaffold_size or scaffold_data["length"],
        "unindex": scaffold_data["unindex"],
        "ligand": scaffold_data["ligand_codes"],
        "pdb_content": scaffold_data["motif_pdb"],
        "select_fixed_atoms": scaffold_data["fixed_atoms"],
        "select_buried": scaffold_data["rasa_targets"] or None,
        "num_designs": args.num_backbones,
        "is_non_loopy": True,
        "cfg_scale": args.cfg_scale,
        "step_scale": args.step_scale,
        "num_timesteps": args.num_timesteps,
        "use_classifier_free_guidance": True,
    }

    # Gamma_0 override
    if args.gamma_0 is not None:
        rfd3_params["gamma_0"] = args.gamma_0

    # H-bond conditioning
    if args.enable_hbond:
        if scaffold_data.get("hbond_acceptors"):
            rfd3_params["select_hbond_acceptor"] = scaffold_data["hbond_acceptors"]
            print(f"  H-bond acceptors: {scaffold_data['hbond_acceptors']}")
        if scaffold_data.get("hbond_donors"):
            rfd3_params["select_hbond_donor"] = scaffold_data["hbond_donors"]
            print(f"  H-bond donors: {scaffold_data['hbond_donors']}")
    else:
        print("  H-bond conditioning: DISABLED")

    print(f"  Generating {args.num_backbones} backbones...")
    print(f"  cfg_scale={args.cfg_scale}, step_scale={args.step_scale}, timesteps={args.num_timesteps}")

    # Batch RFD3 to avoid VRAM/timeout issues on consumer GPUs
    # RTX 3090 (24GB) can handle ~10 designs per batch for 170+ residue proteins
    BATCH_SIZE = 10
    total_designs = rfd3_params.pop("num_designs")
    raw_pdbs = []
    start = time.time()

    num_batches = math.ceil(total_designs / BATCH_SIZE)
    for batch_idx in range(num_batches):
        batch_n = min(BATCH_SIZE, total_designs - batch_idx * BATCH_SIZE)
        print(f"    Batch {batch_idx+1}/{num_batches}: generating {batch_n} backbones...")
        rfd3_result = run_rfd3_api(**rfd3_params, num_designs=batch_n, timeout=1800)

        if rfd3_result["status"] != "completed":
            print(f"[FAIL] RFD3 batch {batch_idx+1} failed: {rfd3_result.get('error')}")
            if raw_pdbs:
                print(f"  Continuing with {len(raw_pdbs)} backbones from previous batches")
                break
            return None

        designs = rfd3_result["result"].get("designs", [])
        batch_pdbs = [d.get("content", "") for d in designs if d.get("content")]
        raw_pdbs.extend(batch_pdbs)
        print(f"    Batch {batch_idx+1}: got {len(batch_pdbs)} backbones ({len(raw_pdbs)} total)")

    elapsed = time.time() - start
    print(f"  Generated: {len(raw_pdbs)} backbones in {elapsed:.1f}s")

    # Merge HETATM if needed
    merged_pdbs = []
    for pdb in raw_pdbs:
        has_hetatm = any(line.startswith("HETATM") for line in pdb.split("\n"))
        if has_hetatm:
            merged_pdbs.append(pdb)
        else:
            merged_pdbs.append(_merge_hetatm_to_backbone(pdb, scaffold_data["motif_pdb"], []))

    # Save PDBs
    bb_dir = Path(output_dir) / "backbones"
    bb_dir.mkdir(parents=True, exist_ok=True)
    for i, pdb in enumerate(merged_pdbs):
        with open(bb_dir / f"backbone_{i+1:03d}.pdb", "w") as f:
            f.write(pdb)

    checkpoint = {
        "num_generated": len(merged_pdbs),
        "rfd3_time": elapsed,
        "backbone_files": [f"backbones/backbone_{i+1:03d}.pdb" for i in range(len(merged_pdbs))],
    }
    save_checkpoint(output_dir, "backbones", checkpoint)

    return merged_pdbs


def stage_clash_filter(args, output_dir, backbone_pdbs):
    """Stage 3: Filter designs with ligand clashes."""
    print("\n" + "=" * 70)
    print("STAGE 3: Clash Filtering")
    print("=" * 70)

    from run_50_designs import _check_ligand_clashes

    ligand_code = args.ligand.upper()
    clean_pdbs = []
    clash_count = 0

    for i, pdb in enumerate(backbone_pdbs):
        clash = _check_ligand_clashes(pdb, ligand_code, clash_threshold=2.0)
        if clash["has_clashes"]:
            clash_count += 1
            if clash_count <= 5:
                print(f"  [CLASH] Backbone {i+1}: {clash['clash_count']} hard clashes, "
                      f"min dist {clash['min_distance']:.2f}A")
        else:
            clean_pdbs.append(pdb)

    pct = 100 * len(clean_pdbs) / max(len(backbone_pdbs), 1)
    print(f"  Passed: {len(clean_pdbs)}/{len(backbone_pdbs)} ({pct:.0f}%)")
    print(f"  Rejected: {clash_count} designs with hard clashes")

    if not clean_pdbs:
        print("[FAIL] All designs had ligand clashes.")
        return None

    # Save filtered PDBs
    filt_dir = Path(output_dir) / "filtered"
    filt_dir.mkdir(parents=True, exist_ok=True)
    for i, pdb in enumerate(clean_pdbs):
        with open(filt_dir / f"clean_{i+1:03d}.pdb", "w") as f:
            f.write(pdb)

    checkpoint = {
        "num_input": len(backbone_pdbs),
        "num_passed": len(clean_pdbs),
        "num_clashed": clash_count,
        "pass_rate": len(clean_pdbs) / max(len(backbone_pdbs), 1),
    }
    save_checkpoint(output_dir, "clash_filter", checkpoint)

    return clean_pdbs


def stage_continuity(args, output_dir, backbone_pdbs):
    """Stage 4: Check backbone continuity (chain breaks)."""
    print("\n" + "=" * 70)
    print("STAGE 4: Backbone Continuity Check")
    print("=" * 70)

    from scaffolding_workflow import check_backbone_continuity

    continuous_pdbs = []
    break_count = 0

    for i, pdb in enumerate(backbone_pdbs):
        cont = check_backbone_continuity(pdb)
        if cont["continuous"]:
            continuous_pdbs.append(pdb)
        else:
            break_count += 1
            if break_count <= 5:
                print(f"  [BREAK] Backbone {i+1}: {cont['num_breaks']} breaks (max CA-CA={cont['max_observed']}A)")

    if continuous_pdbs:
        pct = 100 * len(continuous_pdbs) / max(len(backbone_pdbs), 1)
        print(f"  Continuous: {len(continuous_pdbs)}/{len(backbone_pdbs)} ({pct:.0f}%)")
    else:
        print(f"  [WARN] No continuous backbones. Using all {len(backbone_pdbs)} designs.")
        continuous_pdbs = backbone_pdbs

    checkpoint = {
        "num_input": len(backbone_pdbs),
        "num_continuous": len(continuous_pdbs),
        "num_breaks": break_count,
    }
    save_checkpoint(output_dir, "continuity", checkpoint)

    return continuous_pdbs


def stage_scout_mpnn(args, output_dir, backbone_pdbs, scaffold_data):
    """Stage 5a: Generate 1 scout sequence per backbone for early filtering."""
    print("\n" + "=" * 70)
    print("STAGE 5a: Scout Sequence Generation (1 per backbone)")
    print("=" * 70)

    from api_client import run_mpnn_api

    catalytic_ids = scaffold_data.get("catalytic_residue_ids", [])

    has_metal = bool(args.metal)
    has_ligand = bool(args.ligand)
    if has_metal and has_ligand:
        mpnn_bias = "H:0.5,D:0.5,E:0.5,W:0.3,Y:0.3"
        mpnn_omit = None
    elif has_metal:
        mpnn_bias = "H:1.0,C:1.0,D:0.5,E:0.5"
        mpnn_omit = None
    elif has_ligand:
        mpnn_bias = "W:0.5,Y:0.5,F:0.3"
        mpnn_omit = "C"
    else:
        mpnn_bias = None
        mpnn_omit = "C"

    print(f"  Generating 1 scout sequence per backbone ({len(backbone_pdbs)} backbones)")
    print(f"  Temperature: {args.mpnn_temperature}")

    start = time.time()
    scout_sequences = {}  # backbone_idx -> sequence string

    for i, pdb in enumerate(backbone_pdbs):
        result = run_mpnn_api(
            pdb_content=pdb,
            num_seqs=1,
            temperature=args.mpnn_temperature,
            fixed_positions=catalytic_ids if catalytic_ids else None,
            bias_AA=mpnn_bias,
            omit_AA=mpnn_omit,
            pack_side_chains=True,
            pack_with_ligand_context=True,
            number_of_packs_per_design=args.mpnn_num_packs,
            ligand_cutoff_for_score=args.mpnn_ligand_cutoff,
            use_side_chain_context=bool(catalytic_ids),
            save_stats=True,
        )

        if result["status"] == "completed":
            seqs = result["result"].get("sequences", [])
            for seq_data in seqs:
                content = seq_data.get("content", "")
                for line in content.split("\n"):
                    if line and not line.startswith(">"):
                        scout_sequences[i] = line.strip()
                        break
                if i in scout_sequences:
                    break

        if (i + 1) % 10 == 0 or i + 1 == len(backbone_pdbs):
            print(f"    Processed {i+1}/{len(backbone_pdbs)} backbones...")

    elapsed = time.time() - start
    print(f"  Generated {len(scout_sequences)} scout sequences in {elapsed:.1f}s")

    checkpoint = {
        "scout_sequences": scout_sequences,
        "scout_time": elapsed,
    }
    save_checkpoint(output_dir, "scout_mpnn", checkpoint)

    return scout_sequences


def stage_scout_validation(args, output_dir, backbone_pdbs, scout_sequences):
    """Stage 5b: Validate scout sequences with RF3 to filter bad backbones."""
    print("\n" + "=" * 70)
    print("STAGE 5b: Scout RF3 Validation")
    print("=" * 70)

    from run_50_designs import run_rf3_validation

    threshold = args.scout_ptm_threshold
    print(f"  Scout pTM threshold: {threshold}")

    start = time.time()
    scout_results = {}  # backbone_idx -> {passed, ptm, plddt, sequence}

    for bb_idx, seq in scout_sequences.items():
        bb_idx_int = int(bb_idx)  # JSON keys may be strings after checkpoint load
        if bb_idx_int >= len(backbone_pdbs):
            continue

        backbone_pdb = backbone_pdbs[bb_idx_int]
        name = f"scout_bb{bb_idx_int + 1}"

        rf3_result = run_rf3_validation(seq, backbone_pdb, name)

        if "error" in rf3_result:
            print(f"  BB {bb_idx_int + 1}: RF3 ERROR — {rf3_result['error']}")
            scout_results[bb_idx_int] = {
                "passed": False,
                "ptm": None,
                "plddt": None,
                "sequence": seq,
                "error": rf3_result["error"],
            }
        else:
            ptm = rf3_result["rf3_ptm"]
            plddt = rf3_result["rf3_plddt"]
            rmsd = rf3_result.get("rf3_rmsd")
            passed = ptm is not None and ptm > threshold

            status = "PASS" if passed else "SKIP"
            rmsd_s = f", RMSD={rmsd:.2f}A" if rmsd is not None else ""
            print(f"  BB {bb_idx_int + 1}: [{status}] pTM={ptm:.3f}, pLDDT={plddt:.3f}{rmsd_s}")

            scout_results[bb_idx_int] = {
                "passed": passed,
                "ptm": ptm,
                "plddt": plddt,
                "rmsd": rmsd,
                "sequence": seq,
            }

    elapsed = time.time() - start

    # Summary
    num_passed = sum(1 for r in scout_results.values() if r["passed"])
    num_total = len(scout_results)
    num_skipped = num_total - num_passed
    print(f"\n  Scout validation: {num_passed}/{num_total} backbones passed "
          f"({num_skipped} will be skipped)")
    print(f"  Scout validation time: {elapsed:.1f}s")

    checkpoint = {
        "scout_results": {str(k): v for k, v in scout_results.items()},
        "num_passed": num_passed,
        "num_total": num_total,
        "scout_validation_time": elapsed,
        "threshold": threshold,
    }
    save_checkpoint(output_dir, "scout_validation", checkpoint)

    return scout_results


def stage_mpnn(args, output_dir, backbone_pdbs, scaffold_data, scout_results=None):
    """Stage 6: LigandMPNN sequence design with optimized parameters.

    When scout_results is provided, only processes backbones that passed scout
    validation and generates (num_seqs - 1) additional sequences per backbone,
    prepending the scout sequence to maintain the total at num_seqs.
    """
    stage_label = "STAGE 6" if scout_results is not None else "STAGE 5"
    print("\n" + "=" * 70)
    print(f"{stage_label}: LigandMPNN Sequence Design")
    print("=" * 70)

    from api_client import run_mpnn_api

    catalytic_ids = scaffold_data.get("catalytic_residue_ids", [])

    # Bias based on ligand/metal presence
    has_metal = bool(args.metal)
    has_ligand = bool(args.ligand)
    if has_metal and has_ligand:
        mpnn_bias = "H:0.5,D:0.5,E:0.5,W:0.3,Y:0.3"
        mpnn_omit = None
    elif has_metal:
        mpnn_bias = "H:1.0,C:1.0,D:0.5,E:0.5"
        mpnn_omit = None
    elif has_ligand:
        mpnn_bias = "W:0.5,Y:0.5,F:0.3"
        mpnn_omit = "C"
    else:
        mpnn_bias = None
        mpnn_omit = "C"

    # Determine which backbones to process and how many seqs to generate
    if scout_results is not None:
        passing_indices = sorted(
            int(idx) for idx, r in scout_results.items() if r.get("passed")
        )
        seqs_to_generate = max(1, args.num_seqs - 1)  # Scout provides 1
        print(f"  Scout mode: processing {len(passing_indices)}/{len(backbone_pdbs)} passing backbones")
        print(f"  Generating {seqs_to_generate} additional sequences per backbone (+ 1 scout = {args.num_seqs} total)")
    else:
        passing_indices = list(range(len(backbone_pdbs)))
        seqs_to_generate = args.num_seqs

    print(f"  Temperature: {args.mpnn_temperature}")
    print(f"  Ligand cutoff: {args.mpnn_ligand_cutoff}A")
    print(f"  Packs per design: {args.mpnn_num_packs}")
    print(f"  Fixed catalytic residues: {len(catalytic_ids)}")
    print(f"  Bias: {mpnn_bias}")

    start = time.time()
    all_sequences = []  # flat list: [bb0_seq1, bb0_seq2, ..., bb1_seq1, ...]
    bb_sequence_map = {}  # backbone_idx -> list of sequences (for FASTA output)
    mpnn_stats = []

    for count, bb_idx in enumerate(passing_indices):
        pdb = backbone_pdbs[bb_idx]
        bb_seqs = []

        # Prepend scout sequence if available
        if scout_results is not None:
            scout_data = scout_results.get(bb_idx) or scout_results.get(str(bb_idx))
            if scout_data and scout_data.get("sequence"):
                bb_seqs.append(scout_data["sequence"])

        # Generate additional sequences
        result = run_mpnn_api(
            pdb_content=pdb,
            num_seqs=seqs_to_generate,
            temperature=args.mpnn_temperature,
            fixed_positions=catalytic_ids if catalytic_ids else None,
            bias_AA=mpnn_bias,
            omit_AA=mpnn_omit,
            pack_side_chains=True,
            pack_with_ligand_context=True,
            number_of_packs_per_design=args.mpnn_num_packs,
            ligand_cutoff_for_score=args.mpnn_ligand_cutoff,
            use_side_chain_context=bool(catalytic_ids),
            save_stats=True,
        )

        if result["status"] == "completed":
            seqs = result["result"].get("sequences", [])
            for seq_data in seqs:
                content = seq_data.get("content", "")
                for line in content.split("\n"):
                    if line and not line.startswith(">"):
                        bb_seqs.append(line.strip())

            stats = result["result"].get("stats", {})
            if stats:
                mpnn_stats.append({"backbone": bb_idx + 1, "stats": stats})

        bb_sequence_map[bb_idx] = bb_seqs
        all_sequences.extend(bb_seqs)

        processed = count + 1
        if processed % 10 == 0 or processed == len(passing_indices):
            print(f"    Processed {processed}/{len(passing_indices)} backbones "
                  f"({len(all_sequences)} seqs so far)...")

    elapsed = time.time() - start
    print(f"  Generated: {len(all_sequences)} sequences in {elapsed:.1f}s")

    # Save sequences as FASTA
    seq_file = Path(output_dir) / "sequences.fasta"
    with open(seq_file, "w") as f:
        for bb_idx in sorted(bb_sequence_map.keys()):
            for sq_idx, seq in enumerate(bb_sequence_map[bb_idx]):
                label = "scout" if scout_results is not None and sq_idx == 0 else f"seq{sq_idx + 1}"
                f.write(f">bb{bb_idx + 1}_{label}\n{seq}\n")

    # Save MPNN stats
    if mpnn_stats:
        stats_file = Path(output_dir) / "mpnn_stats.json"
        with open(stats_file, "w") as f:
            json.dump(mpnn_stats, f, indent=2)

    checkpoint = {
        "num_sequences": len(all_sequences),
        "mpnn_time": elapsed,
        "temperature": args.mpnn_temperature,
        "ligand_cutoff": args.mpnn_ligand_cutoff,
        "scout_mode": scout_results is not None,
        "passing_backbone_indices": [int(i) for i in passing_indices],
        "bb_sequence_counts": {str(k): len(v) for k, v in bb_sequence_map.items()},
    }
    save_checkpoint(output_dir, "mpnn", checkpoint)

    return all_sequences


def stage_rf3_validation(args, output_dir, backbone_pdbs, all_sequences):
    """Stage 6: RF3 structure prediction + Kabsch RMSD validation."""
    if not args.validate_rf3:
        print("\n  [SKIP] RF3 validation disabled")
        return []

    from run_50_designs import run_validation_stage

    # Determine ligand SMILES
    ligand_smiles = args.ligand_smiles
    if ligand_smiles is None and args.ligand and args.ligand.upper() == "PQQ":
        ligand_smiles = PQQ_SMILES

    results = run_validation_stage(
        backbone_pdbs=backbone_pdbs,
        all_sequences=all_sequences,
        num_sequences_per_design=args.num_seqs,
        ligand_smiles=ligand_smiles,
        output_dir=Path(output_dir) / "validation",
    )

    # Save checkpoint
    save_results = []
    for r in results:
        sr = dict(r)
        sr.pop("cif_content", None)
        if "sequence" in sr and len(sr["sequence"]) > 20:
            sr["sequence"] = sr["sequence"][:20] + "..."
        save_results.append(sr)

    checkpoint = {"results": save_results}
    save_checkpoint(output_dir, "rf3_validation", checkpoint)

    return results


def stage_filtering(args, output_dir, validation_results):
    """Stage 7: Filter and rank final results."""
    print("\n" + "=" * 70)
    print("STAGE 7: Filtering & Ranking")
    print("=" * 70)

    if not validation_results:
        print("  No validation results to filter")
        return []

    # Score each design
    scored = []
    for r in validation_results:
        score = 0.0
        rmsd = r.get("rf3_rmsd")
        plddt = r.get("rf3_plddt")
        ptm = r.get("rf3_ptm")
        lig_ptm = r.get("rf3_lig_ptm")

        if rmsd is not None and rmsd < 2.0:
            score += 40 * (1.0 - rmsd / 2.0)  # 0-40 points
        if plddt is not None:
            score += 30 * min(plddt, 1.0)  # 0-30 points
        if ptm is not None:
            score += 15 * min(ptm, 1.0)  # 0-15 points
        if lig_ptm is not None:
            score += 15 * min(lig_ptm, 1.0)  # 0-15 points

        r["composite_score"] = round(score, 2)
        scored.append(r)

    # Sort by composite score
    scored.sort(key=lambda x: x.get("composite_score", 0), reverse=True)

    # Print top 10
    print(f"\n  Top designs (by composite score):")
    print(f"  {'Rank':>4} {'BB':>3} {'Seq':>3} | {'Score':>5} {'RMSD':>7} {'pLDDT':>6} {'pTM':>5} {'Pass':>4}")
    for rank, r in enumerate(scored[:10], 1):
        rmsd_s = f"{r['rf3_rmsd']:.2f}A" if r.get("rf3_rmsd") is not None else "N/A"
        plddt_s = f"{r['rf3_plddt']:.3f}" if r.get("rf3_plddt") is not None else "N/A"
        ptm_s = f"{r['rf3_ptm']:.3f}" if r.get("rf3_ptm") is not None else "N/A"
        pass_s = "Y" if r.get("passed_any") else "N"
        print(f"  {rank:>4} {r['backbone']:>3} {r['seq']:>3} | "
              f"{r['composite_score']:>5.1f} {rmsd_s:>7} {plddt_s:>6} {ptm_s:>5} {pass_s:>4}")

    # Save final results
    final_file = Path(output_dir) / "final_results.json"
    save_results = []
    for r in scored:
        sr = dict(r)
        sr.pop("cif_content", None)
        if "sequence" in sr and len(sr["sequence"]) > 20:
            sr["sequence_preview"] = sr["sequence"][:20] + "..."
        save_results.append(sr)
    with open(final_file, "w") as f:
        json.dump(save_results, f, indent=2)

    checkpoint = {
        "num_scored": len(scored),
        "num_passed_rf3": sum(1 for r in scored if r.get("passed_rf3")),
        "num_passed_any": sum(1 for r in scored if r.get("passed_any")),
        "top_score": scored[0]["composite_score"] if scored else 0,
    }
    save_checkpoint(output_dir, "filtering", checkpoint)

    return scored


def stage_summary(args, output_dir, scaffold_data, backbone_pdbs, all_sequences, validation_results,
                   scout_stats=None):
    """Stage 9: Generate human-readable summary."""
    print("\n" + "=" * 70)
    print("STAGE 9: Summary Generation")
    print("=" * 70)

    n_total = len(validation_results) if validation_results else 0
    n_passed_rf3 = sum(1 for r in validation_results if r.get("passed_rf3")) if validation_results else 0
    n_passed_any = sum(1 for r in validation_results if r.get("passed_any")) if validation_results else 0

    rmsds = [r["rf3_rmsd"] for r in validation_results if r.get("rf3_rmsd") is not None] if validation_results else []
    plddts = [r["rf3_plddt"] for r in validation_results if r.get("rf3_plddt") is not None] if validation_results else []
    ptms = [r["rf3_ptm"] for r in validation_results if r.get("rf3_ptm") is not None] if validation_results else []

    summary_lines = [
        f"Production Run Summary",
        f"=" * 50,
        f"Date: {datetime.now().isoformat()}",
        f"Run name: {args.run_name or 'unnamed'}",
        f"",
        f"Source: PDB {args.pdb_id}, Ligand {args.ligand}, Metal {args.metal}",
        f"Enzyme class: {args.enzyme_class}",
        f"",
        f"=== RFD3 Parameters ===",
        f"  cfg_scale: {args.cfg_scale}",
        f"  step_scale: {args.step_scale}",
        f"  num_timesteps: {args.num_timesteps}",
        f"  scaffold_size: {args.scaffold_size or scaffold_data.get('length', 'auto')}",
        f"  H-bond conditioning: {'enabled' if args.enable_hbond else 'disabled'}",
        f"",
        f"=== LigandMPNN Parameters ===",
        f"  temperature: {args.mpnn_temperature}",
        f"  ligand_cutoff: {args.mpnn_ligand_cutoff}A",
        f"  num_packs: {args.mpnn_num_packs}",
        f"",
        f"=== Pipeline Results ===",
        f"  Backbones requested: {args.num_backbones}",
        f"  Backbones after clash filter: {len(backbone_pdbs)}",
    ]

    # Add scout statistics if available
    if scout_stats:
        summary_lines += [
            f"  Scout mode: enabled (pTM threshold={scout_stats.get('threshold', 'N/A')})",
            f"  Backbones scouted: {scout_stats.get('num_total', 'N/A')}",
            f"  Backbones passed scout: {scout_stats.get('num_passed', 'N/A')}",
            f"  Backbones skipped: {scout_stats.get('num_total', 0) - scout_stats.get('num_passed', 0)}",
        ]
    elif hasattr(args, 'scout_mode') and not args.scout_mode:
        summary_lines.append(f"  Scout mode: disabled")

    summary_lines += [
        f"  Total sequences: {len(all_sequences)}",
        f"  Sequences validated: {n_total}",
        f"",
        f"=== Validation Metrics ===",
        f"  Passed RF3 (RMSD<2A + pLDDT>0.7): {n_passed_rf3}/{n_total} ({100*n_passed_rf3/max(n_total,1):.1f}%)",
        f"  Passed ANY criterion: {n_passed_any}/{n_total} ({100*n_passed_any/max(n_total,1):.1f}%)",
    ]

    if rmsds:
        summary_lines += [
            f"",
            f"  RMSD: min={min(rmsds):.2f}A, max={max(rmsds):.2f}A, mean={sum(rmsds)/len(rmsds):.2f}A",
        ]
    if plddts:
        summary_lines += [
            f"  pLDDT: min={min(plddts):.3f}, max={max(plddts):.3f}, mean={sum(plddts)/len(plddts):.3f}",
        ]
    if ptms:
        summary_lines += [
            f"  pTM: min={min(ptms):.3f}, max={max(ptms):.3f}, mean={sum(ptms)/len(ptms):.3f}",
        ]

    summary_text = "\n".join(summary_lines)
    print(summary_text)

    # Save summary
    summary_file = Path(output_dir) / "summary.txt"
    with open(summary_file, "w") as f:
        f.write(summary_text)

    # Save config
    config = {
        "pdb_id": args.pdb_id,
        "ligand": args.ligand,
        "metal": args.metal,
        "enzyme_class": args.enzyme_class,
        "num_backbones": args.num_backbones,
        "num_seqs": args.num_seqs,
        "step_scale": args.step_scale,
        "cfg_scale": args.cfg_scale,
        "scaffold_size": args.scaffold_size,
        "num_timesteps": args.num_timesteps,
        "gamma_0": args.gamma_0,
        "enable_hbond": args.enable_hbond,
        "mpnn_temperature": args.mpnn_temperature,
        "mpnn_ligand_cutoff": args.mpnn_ligand_cutoff,
        "mpnn_num_packs": args.mpnn_num_packs,
        "mpnn_noise_level": args.mpnn_noise_level,
        "validate_rf3": args.validate_rf3,
        "scout_mode": args.scout_mode,
        "scout_ptm_threshold": args.scout_ptm_threshold,
        "run_name": args.run_name,
        "timestamp": datetime.now().isoformat(),
    }
    config_file = Path(output_dir) / "config.json"
    with open(config_file, "w") as f:
        json.dump(config, f, indent=2)

    save_checkpoint(output_dir, "summary", {"completed": True})
    print(f"\n  Summary saved to {summary_file}")
    print(f"  Config saved to {config_file}")

    return summary_text


###############################################################################
# MAIN PIPELINE
###############################################################################

async def run_pipeline(args):
    """Run the full production pipeline."""
    # Setup output directory
    if args.output_dir:
        output_dir = args.output_dir
    else:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        name = args.run_name or f"{args.pdb_id}_{args.ligand}_{args.metal}"
        output_dir = f"output_{name}_{ts}"

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print(f"PRODUCTION PIPELINE — {args.pdb_id} / {args.ligand} / {args.metal}")
    print(f"Output: {output_dir}")
    print("=" * 70)

    total_start = time.time()

    # Check for resume
    resume_stage = None
    if args.resume:
        output_dir = Path(args.resume)
        resume_stage = get_resume_stage(output_dir)
        if resume_stage:
            print(f"\n  Resuming from stage: {resume_stage}")

    # Stage 1: Scaffold
    scaffold_data = load_checkpoint(output_dir, "scaffold") if resume_stage else None
    if scaffold_data is None:
        scaffold_data = await stage_scaffold(args, output_dir)
        if scaffold_data is None:
            return 1

    # Stage 2: Backbones
    backbone_pdbs = None
    bb_checkpoint = load_checkpoint(output_dir, "backbones") if resume_stage else None
    if bb_checkpoint:
        # Reload PDBs from files
        backbone_pdbs = []
        for f in bb_checkpoint.get("backbone_files", []):
            pdb_path = output_dir / f
            if pdb_path.exists():
                backbone_pdbs.append(pdb_path.read_text())
        print(f"\n  [RESUME] Loaded {len(backbone_pdbs)} backbones from checkpoint")
    else:
        backbone_pdbs = stage_backbones(args, output_dir, scaffold_data)
        if backbone_pdbs is None:
            return 1

    # Stage 3: Clash filtering
    clash_cp = load_checkpoint(output_dir, "clash_filter") if resume_stage else None
    if clash_cp is None:
        backbone_pdbs = stage_clash_filter(args, output_dir, backbone_pdbs)
        if backbone_pdbs is None:
            return 1
    else:
        # Reload filtered PDBs
        filt_dir = output_dir / "filtered"
        if filt_dir.exists():
            backbone_pdbs = []
            for f in sorted(filt_dir.glob("clean_*.pdb")):
                backbone_pdbs.append(f.read_text())
            print(f"\n  [RESUME] Loaded {len(backbone_pdbs)} filtered backbones")

    # Stage 4: Continuity
    cont_cp = load_checkpoint(output_dir, "continuity") if resume_stage else None
    if cont_cp is None:
        backbone_pdbs = stage_continuity(args, output_dir, backbone_pdbs)

    # Stage 5a/5b: Scout (optional)
    scout_results = None
    scout_stats = None
    backbone_pdbs_for_validation = backbone_pdbs  # default: validate all
    if args.scout_mode and args.num_seqs > 1:
        # Stage 5a: Scout MPNN
        scout_mpnn_cp = load_checkpoint(output_dir, "scout_mpnn") if resume_stage else None
        if scout_mpnn_cp:
            scout_sequences = scout_mpnn_cp.get("scout_sequences", {})
            print(f"\n  [RESUME] Loaded {len(scout_sequences)} scout sequences from checkpoint")
        else:
            scout_sequences = stage_scout_mpnn(args, output_dir, backbone_pdbs, scaffold_data)

        # Stage 5b: Scout validation
        scout_val_cp = load_checkpoint(output_dir, "scout_validation") if resume_stage else None
        if scout_val_cp:
            scout_results = scout_val_cp.get("scout_results", {})
            scout_stats = scout_val_cp
            num_passed = scout_val_cp.get("num_passed", 0)
            num_total = scout_val_cp.get("num_total", 0)
            print(f"\n  [RESUME] Loaded scout validation: {num_passed}/{num_total} passed")
        else:
            scout_results = stage_scout_validation(args, output_dir, backbone_pdbs, scout_sequences)
            scout_stats = {
                "num_passed": sum(1 for r in scout_results.values() if r.get("passed")),
                "num_total": len(scout_results),
                "threshold": args.scout_ptm_threshold,
            }

        # Filter backbone list to only passing backbones for downstream stages
        passing_indices = sorted(
            int(idx) for idx, r in scout_results.items() if r.get("passed")
        )
        if not passing_indices:
            print("\n  [WARN] No backbones passed scout validation. Running all backbones without scout filtering.")
            scout_results = None
            scout_stats["fallback"] = True
        else:
            # For validation stage, only pass the backbone PDBs that passed scout
            backbone_pdbs_for_validation = [backbone_pdbs[i] for i in passing_indices]

    # Stage 6: MPNN
    all_sequences = None
    mpnn_cp = load_checkpoint(output_dir, "mpnn") if resume_stage else None
    if mpnn_cp:
        # Reload sequences from FASTA
        seq_file = output_dir / "sequences.fasta"
        if seq_file.exists():
            all_sequences = []
            for line in seq_file.read_text().splitlines():
                if line and not line.startswith(">"):
                    all_sequences.append(line.strip())
            print(f"\n  [RESUME] Loaded {len(all_sequences)} sequences from checkpoint")
    if all_sequences is None:
        all_sequences = stage_mpnn(args, output_dir, backbone_pdbs, scaffold_data,
                                   scout_results=scout_results)

    # Stage 7: RF3 validation
    # backbone_pdbs_for_validation is set to only passing backbones when scout mode is active
    validation_backbone_pdbs = backbone_pdbs_for_validation

    validation_results = None
    rf3_cp = load_checkpoint(output_dir, "rf3_validation") if resume_stage else None
    if rf3_cp:
        validation_results = rf3_cp.get("results", [])
        print(f"\n  [RESUME] Loaded {len(validation_results)} validation results")
    if validation_results is None:
        validation_results = stage_rf3_validation(args, output_dir, validation_backbone_pdbs, all_sequences)

    # Stage 8: Filtering & ranking
    scored = stage_filtering(args, output_dir, validation_results)

    # Stage 9: Summary
    stage_summary(args, output_dir, scaffold_data, backbone_pdbs, all_sequences, validation_results,
                  scout_stats=scout_stats)

    total_time = time.time() - total_start
    print(f"\n{'=' * 70}")
    print(f"PIPELINE COMPLETE — Total time: {total_time:.1f}s ({total_time/60:.1f} min)")
    print(f"Output: {output_dir}")
    print(f"{'=' * 70}")

    return 0


###############################################################################
# SWEEP MODE
###############################################################################

async def run_sweep(args):
    """Run parameter sweep from JSON config."""
    config_path = Path(args.sweep_config)
    if not config_path.exists():
        print(f"[FAIL] Sweep config not found: {config_path}")
        return 1

    with open(config_path) as f:
        sweep_config = json.load(f)

    base = sweep_config.get("base", {})
    grid = sweep_config.get("grid", {})

    # Build all parameter combinations
    param_names = list(grid.keys())
    param_values = list(grid.values())

    # Cartesian product
    combos = [{}]
    for name, values in zip(param_names, param_values):
        new_combos = []
        for combo in combos:
            for val in values:
                new_combo = dict(combo)
                new_combo[name] = val
                new_combos.append(new_combo)
        combos = new_combos

    print(f"Sweep: {len(combos)} configurations from {len(param_names)} parameters")

    output_base = Path(args.output_dir or "sweep_output")
    output_base.mkdir(parents=True, exist_ok=True)

    results_summary = []

    for i, combo in enumerate(combos):
        combo_name = "_".join(f"{k}{v}" for k, v in combo.items())
        combo_dir = output_base / f"config_{i+1:02d}_{combo_name}"

        print(f"\n{'#' * 70}")
        print(f"SWEEP CONFIG {i+1}/{len(combos)}: {combo}")
        print(f"{'#' * 70}")

        # Build args for this combo
        sweep_args = argparse.Namespace(**vars(args))
        # Apply base params
        for k, v in base.items():
            k_attr = k.replace("-", "_")
            setattr(sweep_args, k_attr, v)
        # Apply combo params
        for k, v in combo.items():
            k_attr = k.replace("-", "_")
            setattr(sweep_args, k_attr, v)

        sweep_args.output_dir = str(combo_dir)
        sweep_args.run_name = combo_name
        sweep_args.sweep_config = None  # Don't recurse

        rc = await run_pipeline(sweep_args)

        # Collect result summary
        filtering_cp = load_checkpoint(combo_dir, "filtering")
        results_summary.append({
            "config_index": i + 1,
            "params": combo,
            "output_dir": str(combo_dir),
            "exit_code": rc,
            "num_passed_rf3": filtering_cp.get("num_passed_rf3", 0) if filtering_cp else 0,
            "num_passed_any": filtering_cp.get("num_passed_any", 0) if filtering_cp else 0,
            "top_score": filtering_cp.get("top_score", 0) if filtering_cp else 0,
        })

    # Save sweep summary
    sweep_summary_file = output_base / "sweep_summary.json"
    with open(sweep_summary_file, "w") as f:
        json.dump(results_summary, f, indent=2)

    # Print sweep results
    print(f"\n{'=' * 70}")
    print("SWEEP RESULTS")
    print(f"{'=' * 70}")
    print(f"{'Config':>6} | {'Passed RF3':>10} {'Passed Any':>10} {'Top Score':>9} | Params")
    for r in sorted(results_summary, key=lambda x: x["top_score"], reverse=True):
        print(f"  {r['config_index']:>4} | {r['num_passed_rf3']:>10} {r['num_passed_any']:>10} "
              f"{r['top_score']:>9.1f} | {r['params']}")

    print(f"\nSweep summary saved to {sweep_summary_file}")
    return 0


###############################################################################
# ENTRY POINT
###############################################################################

def main():
    parser = build_parser()
    args = parser.parse_args()

    if args.sweep_config:
        return asyncio.run(run_sweep(args))
    else:
        return asyncio.run(run_pipeline(args))


if __name__ == "__main__":
    sys.exit(main())
