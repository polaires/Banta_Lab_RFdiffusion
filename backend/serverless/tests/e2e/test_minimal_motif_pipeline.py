"""
Test minimal motif pipeline with 25 designs at proper 140-170 size.
Runs in batches to avoid API timeout.
"""
import sys
import time
import json
import math
import builtins
from datetime import datetime
from pathlib import Path

# Force unbuffered output
_orig_print = builtins.print
def print(*args, **kwargs):
    kwargs.setdefault('flush', True)
    _orig_print(*args, **kwargs)

sys.path.insert(0, r"G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless")

from scaffolding_workflow import ScaffoldingWorkflow, check_backbone_continuity
from api_client import run_rfd3_api, run_mpnn_api
from run_50_designs import _merge_hetatm_to_backbone, _check_ligand_clashes
import asyncio


def main():
    NUM_DESIGNS = 10
    NUM_SEQS = 4
    BATCH_SIZE = 2  # 2 designs per API call (140-170 res + CFG is slow)

    print("=" * 70)
    print(f"MINIMAL MOTIF PIPELINE TEST - {NUM_DESIGNS} DESIGNS")
    print("=" * 70)
    print(f"Batch size: {BATCH_SIZE} designs per API call")
    print(f"Sequences per backbone: {NUM_SEQS}")
    print()

    # Step 1: Extract scaffold with minimal motif
    print("Step 1: Extracting scaffold with minimal motif...")
    workflow = ScaffoldingWorkflow()
    scaffold_result = asyncio.run(workflow.run(
        pdb_id="4CVB",
        metal="CA",
        ligand_code="PQQ",
        include_all_ligand_contacts=True,
        use_minimal_motif=True,
    ))

    if not scaffold_result or not scaffold_result.success:
        print(f"[FAIL] Scaffold extraction failed")
        return

    print(f"  Length: {scaffold_result.length}")
    print(f"  Unindex: {scaffold_result.unindex}")
    print(f"  Fixed atoms: {len(scaffold_result.fixed_atoms)} entries")
    for k, v in scaffold_result.fixed_atoms.items():
        print(f"    {k}: {v}")

    if scaffold_result.motif_result:
        mr = scaffold_result.motif_result
        print(f"\n  Motif: {len(mr.motif_residues)} residues (Tier 1-3)")
        print(f"  Context: {len(mr.context_residues)} excluded (Tier 4)")
        print(f"  Fixed atoms: {mr.total_fixed_atoms}")
        for r in mr.motif_residues:
            print(f"    T{r.tier} {r.chain}{r.resnum} {r.resname} {r.fix_type} atoms={r.functional_atoms}")

    # Step 2: Run RFD3 in batches
    print(f"\nStep 2: Running RFD3 ({NUM_DESIGNS} designs in batches of {BATCH_SIZE})...")
    all_raw_pdbs = []

    num_batches = math.ceil(NUM_DESIGNS / BATCH_SIZE)
    for batch_idx in range(num_batches):
        batch_num = min(BATCH_SIZE, NUM_DESIGNS - batch_idx * BATCH_SIZE)
        print(f"\n  Batch {batch_idx+1}/{num_batches}: {batch_num} designs...")

        # Build burial targets: bury both ligand and metal
        burial_targets = {}
        if scaffold_result.rasa_targets:
            burial_targets.update(scaffold_result.rasa_targets)
        burial_targets["PQQ"] = "ALL"  # Bury ligand to create proper pocket

        rfd3_params = {
            "length": scaffold_result.length,
            "unindex": scaffold_result.unindex,
            "ligand": scaffold_result.ligand_codes,
            "pdb_content": scaffold_result.motif_pdb,
            "select_fixed_atoms": scaffold_result.fixed_atoms,
            "select_buried": burial_targets,  # Bury ligand+metal
            "num_designs": batch_num,
            "is_non_loopy": True,
            "cfg_scale": 2.0,  # Standard, needed for RASA + H-bond conditioning
            "use_classifier_free_guidance": True,
        }

        # H-bond conditioning: tell RFD3 that PQQ oxygens are H-bond acceptors
        # This guides RFD3 to place H-bond donors near them, not backbone through them
        if scaffold_result.hbond_acceptors:
            rfd3_params["select_hbond_acceptor"] = scaffold_result.hbond_acceptors
            if batch_idx == 0:
                print(f"    H-bond acceptors: {scaffold_result.hbond_acceptors}")
        if scaffold_result.hbond_donors:
            rfd3_params["select_hbond_donor"] = scaffold_result.hbond_donors
            if batch_idx == 0:
                print(f"    H-bond donors: {scaffold_result.hbond_donors}")

        start = time.time()
        rfd3_result = run_rfd3_api(**rfd3_params, timeout=1800)
        elapsed = time.time() - start

        if rfd3_result["status"] != "completed":
            print(f"    [FAIL] Batch {batch_idx+1}: {rfd3_result.get('error')}")
            continue

        designs = rfd3_result["result"].get("designs", [])
        pdbs = [d.get("content", "") for d in designs if d.get("content")]
        all_raw_pdbs.extend(pdbs)
        print(f"    Got {len(pdbs)} designs in {elapsed:.1f}s")

    print(f"\n  Total raw designs: {len(all_raw_pdbs)}")

    if not all_raw_pdbs:
        print("[FAIL] No designs generated")
        return

    # Step 3: Check design sizes
    print(f"\nStep 3: Design sizes...")
    sizes = []
    for pdb in all_raw_pdbs:
        n_res = len(set(
            line[22:27].strip()
            for line in pdb.split('\n')
            if line.startswith('ATOM') and line[12:16].strip() == 'CA'
        ))
        sizes.append(n_res)
    print(f"  Sizes: min={min(sizes)}, max={max(sizes)}, avg={sum(sizes)/len(sizes):.0f}")

    # Step 4: Merge HETATM
    print(f"\nStep 4: Merging HETATM records...")
    merged_pdbs = []
    for pdb in all_raw_pdbs:
        has_hetatm = any(line.startswith('HETATM') for line in pdb.split('\n'))
        if has_hetatm:
            merged_pdbs.append(pdb)
        else:
            merged = _merge_hetatm_to_backbone(pdb, scaffold_result.motif_pdb, [])
            merged_pdbs.append(merged)
    print(f"  Merged: {len(merged_pdbs)} designs")

    # Step 5: Filter clashes (tiered: hard clashes rejected, polar contacts tolerated)
    print(f"\nStep 5: Filtering ligand clashes (tiered filter)...")
    clean_pdbs = []
    clash_count = 0
    soft_total = 0
    for i, pdb in enumerate(merged_pdbs):
        clash = _check_ligand_clashes(pdb, "PQQ", clash_threshold=2.0)
        soft = clash.get('soft_contacts', 0)
        soft_total += soft
        if clash['has_clashes']:
            clash_count += 1
            if clash_count <= 5:
                print(f"  [CLASH] Design {i+1}: {clash['clash_count']} hard clashes, "
                      f"{soft} soft contacts, min dist {clash['min_distance']:.2f}A")
                for c in clash.get('clashes', []):
                    print(f"          {c['protein']} <-> {c['ligand']} = {c['distance']:.2f}A")
        else:
            clean_pdbs.append(pdb)
            if soft > 0:
                print(f"  [PASS+] Design {i+1}: {soft} soft polar contacts (tolerated), min dist {clash['min_distance']:.2f}A")
            else:
                print(f"  [PASS]  Design {i+1}: clean, min dist {clash['min_distance']:.2f}A")
    print(f"  Passed: {len(clean_pdbs)}/{len(merged_pdbs)} ({100*len(clean_pdbs)/max(len(merged_pdbs),1):.0f}%)")
    print(f"  Soft polar contacts (tolerated): {soft_total} across all designs")

    if not clean_pdbs:
        print("[FAIL] All designs clashed. Using all anyway for analysis.")
        clean_pdbs = merged_pdbs

    # Step 6: Check backbone continuity
    print(f"\nStep 6: Backbone continuity...")
    continuous_pdbs = []
    break_details = []
    for i, pdb in enumerate(clean_pdbs):
        cont = check_backbone_continuity(pdb)
        if cont["continuous"]:
            continuous_pdbs.append(pdb)
        else:
            break_details.append(cont)
            if len(break_details) <= 3:
                print(f"  [BREAK] Design {i+1}: {cont['num_breaks']} breaks (max CA-CA={cont['max_observed']:.1f}A)")

    cont_rate = len(continuous_pdbs) / max(len(clean_pdbs), 1) * 100
    print(f"  Continuous: {len(continuous_pdbs)}/{len(clean_pdbs)} ({cont_rate:.0f}%)")

    # Use continuous if available, else fall back
    final_pdbs = continuous_pdbs if continuous_pdbs else clean_pdbs
    print(f"  Using {len(final_pdbs)} designs for MPNN")

    # Step 7: Run MPNN
    print(f"\nStep 7: Running LigandMPNN ({NUM_SEQS} seqs/backbone)...")
    catalytic_ids = scaffold_result.source_info.get("catalytic_residue_ids", [])
    print(f"  Fixed positions: {catalytic_ids}")

    all_sequences = []
    start = time.time()
    for i, pdb in enumerate(final_pdbs):
        mpnn_result = run_mpnn_api(
            pdb_content=pdb,
            num_seqs=NUM_SEQS,
            temperature=0.1,
            fixed_positions=catalytic_ids if catalytic_ids else None,
            bias_AA="H:0.5,D:0.5,E:0.5,W:0.3,Y:0.3",
            omit_AA=None,
            pack_side_chains=True,
            pack_with_ligand_context=True,
            ligand_cutoff_for_score=5.0,
            use_side_chain_context=True,
            save_stats=True,
        )
        if mpnn_result["status"] == "completed":
            seqs = mpnn_result["result"].get("sequences", [])
            for seq_data in seqs:
                content = seq_data.get("content", "")
                for line in content.split('\n'):
                    if line and not line.startswith('>'):
                        all_sequences.append(line.strip())
        if (i + 1) % 5 == 0:
            print(f"    Processed {i+1}/{len(final_pdbs)} backbones...")

    mpnn_time = time.time() - start
    print(f"  Generated: {len(all_sequences)} sequences in {mpnn_time:.1f}s")

    # Step 8: Summary
    print("\n" + "=" * 70)
    print("RESULTS SUMMARY")
    print("=" * 70)
    print(f"  Scaffold length: {scaffold_result.length}")
    print(f"  Motif residues: {len(scaffold_result.motif_result.motif_residues) if scaffold_result.motif_result else 'N/A'}")
    print(f"  Design sizes: {min(sizes)}-{max(sizes)} (avg {sum(sizes)/len(sizes):.0f})")
    print(f"  Raw designs: {len(all_raw_pdbs)}")
    print(f"  After clash filter: {len(clean_pdbs)}")
    print(f"  Continuous backbones: {len(continuous_pdbs)}/{len(clean_pdbs)} ({cont_rate:.0f}%)")
    print(f"  Final backbones: {len(final_pdbs)}")
    print(f"  Total sequences: {len(all_sequences)}")

    # Save for analysis
    output_dir = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless\output_minimal_motif_25")
    output_dir.mkdir(exist_ok=True)

    # Save backbone PDBs
    for i, pdb in enumerate(final_pdbs):
        (output_dir / f"backbone_{i+1}.pdb").write_text(pdb)

    # Save sequences
    with open(output_dir / "sequences.fasta", "w") as f:
        for i, seq in enumerate(all_sequences):
            f.write(f">design_{i+1}\n{seq}\n")

    # Save summary
    summary = {
        "timestamp": datetime.now().isoformat(),
        "scaffold_length": scaffold_result.length,
        "motif_residues": len(scaffold_result.motif_result.motif_residues) if scaffold_result.motif_result else 0,
        "design_sizes": sizes,
        "raw_designs": len(all_raw_pdbs),
        "clash_free": len(clean_pdbs),
        "continuous": len(continuous_pdbs),
        "continuity_rate": cont_rate,
        "final_backbones": len(final_pdbs),
        "total_sequences": len(all_sequences),
        "fixed_atoms": {k: v for k, v in scaffold_result.fixed_atoms.items()},
    }
    with open(output_dir / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n  Output saved to: {output_dir}")
    print("=" * 70)


if __name__ == "__main__":
    main()
