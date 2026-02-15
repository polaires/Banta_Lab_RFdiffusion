#!/usr/bin/env python3
import json
from collections import defaultdict
from pathlib import Path

ANALYSIS_DIR = Path(__file__).parent
OUTPUTS_DIR = ANALYSIS_DIR.parent / "outputs"
BACKEND_DIR = ANALYSIS_DIR.parent.parent.parent / "backend" / "serverless"
PREFILTER_CN = 6
PREFILTER_LIGAND_CONTACTS = 3


def load_all_designs():
    all_designs = []
    for i in range(1, 6):
        fpath = ANALYSIS_DIR / f"exploration_round_0{i}.json"
        if not fpath.exists():
            print(f"  [SKIP] {fpath.name} not found")
            continue
        data = json.load(open(fpath))
        designs = data.get("designs", [])
        for d in designs:
            d["_round"] = i
        all_designs.extend(designs)
        print(f"  Round {i}: {len(designs)} designs loaded")
    return all_designs


def group_by_backbone(designs):
    backbones = defaultdict(lambda: {
        "designs": [], "cn_values": [], "lc_values": [],
        "plddt_values": [], "ptm_values": [],
        "any_passed": False, "pass_count": 0, "rounds": set(),
    })
    for d in designs:
        bid = d.get("backbone_id", "unknown")
        b = backbones[bid]
        b["designs"].append(d)
        b["rounds"].add(d["_round"])
        cn = d.get("coordination_number")
        lc = d.get("ligand_contacts")
        plddt = d.get("plddt")
        ptm = d.get("ptm")
        if cn is not None: b["cn_values"].append(cn)
        if lc is not None: b["lc_values"].append(lc)
        if plddt is not None: b["plddt_values"].append(plddt)
        if ptm is not None: b["ptm_values"].append(ptm)
        if d.get("filter_passed", False):
            b["any_passed"] = True
            b["pass_count"] += 1
    return dict(backbones)


def stats(values):
    if not values:
        return None, None, None
    return min(values), max(values), sum(values) / len(values)


def check_backbone_pdbs():
    print()
    print("=" * 70)
    print("BACKBONE PDB FILE CHECK")
    print("=" * 70)
    if OUTPUTS_DIR.exists():
        subdirs = sorted([d for d in OUTPUTS_DIR.iterdir() if d.is_dir()])
        print(f"Experiment outputs ({OUTPUTS_DIR}):")
        for sd in subdirs:
            pdbs = list(sd.glob("*.pdb"))
            print(f"  {sd.name}: {len(pdbs)} PDB files")
    else:
        print(f"  Outputs dir not found: {OUTPUTS_DIR}")
    print("Backend exploration outputs:")
    for i in range(1, 6):
        odir = BACKEND_DIR / f"output_exploration_r{i}"
        if odir.exists():
            pdbs = list(odir.glob("*.pdb"))
            print(f"  output_exploration_r{i}: {len(pdbs)} PDB files")
        else:
            print(f"  output_exploration_r{i}: not found")


def main():
    sep = "=" * 70
    sep2 = "=" * 60

    print(sep)
    print("PRE-FILTER VALIDATION ANALYSIS")
    print(f"Proposed thresholds: CN >= {PREFILTER_CN}, ligand_contacts >= {PREFILTER_LIGAND_CONTACTS}")
    print(sep)

    print("Loading exploration round data...")
    designs = load_all_designs()
    print(f"Total designs loaded: {len(designs)}")

    backbones = group_by_backbone(designs)
    promising = {k: v for k, v in backbones.items() if v["any_passed"]}
    non_promising = {k: v for k, v in backbones.items() if not v["any_passed"]}

    print(f"Total unique backbones: {len(backbones)}")
    print(f"  Promising (any seq passed metal_binding): {len(promising)}")
    print(f"  Non-promising (no seq passed): {len(non_promising)}")

    for label, group in [("PROMISING", promising), ("NON-PROMISING", non_promising)]:
        print(f"\n{sep2}")
        print(f"  {label} BACKBONES ({len(group)})")
        print(sep2)
        if not group:
            print("  (none)")
            continue
        all_cn, all_lc, all_plddt, all_ptm = [], [], [], []
        for v in group.values():
            all_cn.extend(v["cn_values"])
            all_lc.extend(v["lc_values"])
            all_plddt.extend(v["plddt_values"])
            all_ptm.extend(v["ptm_values"])
        cn_s = stats(all_cn)
        lc_s = stats(all_lc)
        pl_s = stats(all_plddt)
        pt_s = stats(all_ptm)
        total_d = sum(len(v["designs"]) for v in group.values())
        print(f"  Designs in group: {total_d}")
        if cn_s[2] is not None:
            print(f"  Coordination Number: min={cn_s[0]}, max={cn_s[1]}, avg={cn_s[2]:.2f}")
        if lc_s[2] is not None:
            print(f"  Ligand Contacts:     min={lc_s[0]}, max={lc_s[1]}, avg={lc_s[2]:.2f}")
        if pl_s[2] is not None:
            print(f"  pLDDT:               min={pl_s[0]:.3f}, max={pl_s[1]:.3f}, avg={pl_s[2]:.3f}")
        if pt_s[2] is not None:
            print(f"  pTM:                 min={pt_s[0]:.3f}, max={pt_s[1]:.3f}, avg={pt_s[2]:.3f}")

    print(f"\n{sep}")
    print("KEY QUESTION: Would pre-filter reject promising backbones?")
    print(f"Pre-filter: CN >= {PREFILTER_CN} AND ligand_contacts >= {PREFILTER_LIGAND_CONTACTS}")
    print(sep)

    rejected_promising = []
    passed_promising = []
    for bid, v in promising.items():
        best_cn = max(v["cn_values"]) if v["cn_values"] else 0
        best_lc = max(v["lc_values"]) if v["lc_values"] else 0
        if (best_cn >= PREFILTER_CN) and (best_lc >= PREFILTER_LIGAND_CONTACTS):
            passed_promising.append((bid, v, best_cn, best_lc))
        else:
            rejected_promising.append((bid, v, best_cn, best_lc))

    rejected_promising_avg = []
    for bid, v in promising.items():
        avg_cn = sum(v["cn_values"]) / len(v["cn_values"]) if v["cn_values"] else 0
        avg_lc = sum(v["lc_values"]) / len(v["lc_values"]) if v["lc_values"] else 0
        if not ((avg_cn >= PREFILTER_CN) and (avg_lc >= PREFILTER_LIGAND_CONTACTS)):
            rejected_promising_avg.append((bid, v, avg_cn, avg_lc))

    print(f"Using BEST (max) metrics per backbone:")
    print(f"  Promising backbones that PASS pre-filter:  {len(passed_promising)}")
    print(f"  Promising backbones that FAIL pre-filter:  {len(rejected_promising)}")

    if rejected_promising:
        print(f"\n  *** WARNING: {len(rejected_promising)} promising backbone(s) would be REJECTED! ***")
        for bid, v, cn, lc in rejected_promising:
            n_pass = v["pass_count"]
            n_total = len(v["designs"])
            avg_plddt = sum(v["plddt_values"]) / len(v["plddt_values"]) if v["plddt_values"] else 0
            avg_ptm = sum(v["ptm_values"]) / len(v["ptm_values"]) if v["ptm_values"] else 0
            best_composite = max((d.get("composite_score", 0) or 0) for d in v["designs"])
            rds = sorted(v["rounds"])
            print(f"\n    Backbone: {bid}")
            print(f"      Best CN: {cn}, Best LC: {lc}")
            print(f"      Seqs: {n_total} total, {n_pass} passed metal_binding")
            print(f"      Avg pLDDT: {avg_plddt:.3f}, Avg pTM: {avg_ptm:.3f}")
            print(f"      Best composite: {best_composite:.4f}")
            print(f"      Rounds: {rds}")
            fail_reasons = []
            if cn < PREFILTER_CN:
                fail_reasons.append(f"CN {cn} < {PREFILTER_CN}")
            if lc < PREFILTER_LIGAND_CONTACTS:
                fail_reasons.append(f"LC {lc} < {PREFILTER_LIGAND_CONTACTS}")
            reason_str = ", ".join(fail_reasons)
            print(f"      Reason: {reason_str}")
    else:
        print("  All promising backbones PASS the pre-filter (using best metrics).")

    print(f"Using AVERAGE metrics per backbone:")
    print(f"  Promising backbones that FAIL pre-filter:  {len(rejected_promising_avg)}")
    if rejected_promising_avg:
        for bid, v, cn, lc in rejected_promising_avg:
            print(f"    {bid}: avg CN={cn:.1f}, avg LC={lc:.1f}")

    print(f"\n{sep2}")
    print("PRE-FILTER EFFECTIVENESS ON NON-PROMISING BACKBONES")
    print(sep2)

    rejected_np = 0
    passed_np = 0
    for bid, v in non_promising.items():
        best_cn = max(v["cn_values"]) if v["cn_values"] else 0
        best_lc = max(v["lc_values"]) if v["lc_values"] else 0
        if (best_cn >= PREFILTER_CN) and (best_lc >= PREFILTER_LIGAND_CONTACTS):
            passed_np += 1
        else:
            rejected_np += 1

    print(f"  Non-promising correctly rejected: {rejected_np} / {len(non_promising)}")
    print(f"  Non-promising that would pass (wasted compute): {passed_np} / {len(non_promising)}")
    if len(non_promising) > 0:
        rej_pct = rejected_np / len(non_promising) * 100
        print(f"  Rejection rate: {rej_pct:.1f}%")

    print(f"\n{sep2}")
    print("COORDINATION NUMBER DISTRIBUTION BY BACKBONE STATUS")
    print(sep2)
    for label, group in [("Promising", promising), ("Non-promising", non_promising)]:
        cn_dist = defaultdict(int)
        for v in group.values():
            for cn in v["cn_values"]:
                cn_dist[cn] += 1
        if cn_dist:
            print(f"\n  {label}:")
            for cn in sorted(cn_dist.keys()):
                bar = "#" * min(cn_dist[cn], 50)
                marker = " <-- below threshold" if cn < PREFILTER_CN else ""
                print(f"    CN={cn:2d}: {cn_dist[cn]:4d} designs {bar}{marker}")

    print(f"\n{sep2}")
    print("LIGAND CONTACTS DISTRIBUTION BY BACKBONE STATUS")
    print(sep2)
    for label, group in [("Promising", promising), ("Non-promising", non_promising)]:
        lc_dist = defaultdict(int)
        for v in group.values():
            for lc in v["lc_values"]:
                lc_dist[lc] += 1
        if lc_dist:
            print(f"\n  {label}:")
            for lc in sorted(lc_dist.keys()):
                bar = "#" * min(lc_dist[lc], 50)
                marker = " <-- below threshold" if lc < PREFILTER_LIGAND_CONTACTS else ""
                print(f"    LC={lc:2d}: {lc_dist[lc]:4d} designs {bar}{marker}")

    print(f"\n{sep}")
    print("SUMMARY")
    print(sep)
    print(f"Total backbones: {len(backbones)}")
    print(f"Promising: {len(promising)}, Non-promising: {len(non_promising)}")
    print(f"Pre-filter (CN>={PREFILTER_CN}, LC>={PREFILTER_LIGAND_CONTACTS}):")
    print(f"  Promising rejected (false negatives): {len(rejected_promising)}")
    print(f"  Non-promising rejected (true negatives): {rejected_np}")
    print(f"  Non-promising passed (false positives): {passed_np}")
    safe = len(rejected_promising) == 0
    verdict = "no" if safe else str(len(rejected_promising))
    status = "SAFE" if safe else "NOT SAFE"
    print(f"Conclusion: Pre-filter is {status} - {verdict} promising backbone(s) would be lost.")

    check_backbone_pdbs()


if __name__ == "__main__":
    main()
