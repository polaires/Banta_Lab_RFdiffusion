"""
Generalized Metal-Ligand Binding Protein Design Pipeline

Supports different metal-ligand complexes by configuring:
1. Metal properties (coordination number, preferred donors)
2. Ligand properties (H-bond acceptors/donors, coordinating atoms)
3. Design parameters (size, CFG, filters)

Usage:
    # TB-Citrate (default)
    python pipeline_metal_ligand.py --metal TB --ligand CIT --mode sweep

    # Calcium-EDTA
    python pipeline_metal_ligand.py --metal CA --ligand EDT --mode sweep

    # Zinc-Histidine motif
    python pipeline_metal_ligand.py --metal ZN --ligand HIS --mode sweep

    # Custom complex from config file
    python pipeline_metal_ligand.py --config my_complex.json --mode sweep
"""

import argparse
import requests
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Any, Optional
from dataclasses import dataclass, asdict, field
import statistics

API_URL = "http://localhost:8000/runsync"

# ============================================================================
# METAL-LIGAND TEMPLATES
# ============================================================================

@dataclass
class MetalProperties:
    """Properties of a metal ion for binding site design."""
    code: str                          # PDB code (TB, CA, ZN, etc.)
    name: str                          # Full name
    charge: int                        # Formal charge
    coordination_number: int           # Target CN
    preferred_donors: List[str]        # O, N, S preferences
    ionic_radius: float                # Angstroms
    hsab_class: str                    # hard/soft/borderline

    def to_dict(self):
        return asdict(self)


@dataclass
class LigandProperties:
    """Properties of a ligand for binding site design."""
    code: str                          # PDB code (CIT, EDT, PQQ, etc.)
    name: str                          # Full name
    coordinating_atoms: List[str]      # Atoms that coordinate metal
    hbond_acceptors: List[str]         # O/N atoms for protein H-bonds
    hbond_donors: List[str]            # OH/NH atoms for protein H-bonds
    smiles: str                        # For AF3 submission
    metal_donors: int                  # Number of atoms coordinating metal

    def to_dict(self):
        return asdict(self)


@dataclass
class MetalLigandComplex:
    """Complete metal-ligand complex configuration."""
    metal: MetalProperties
    ligand: LigandProperties
    pdb_path: str                      # Input PDB with complex
    metal_chain: str = "X"             # Chain ID for metal
    ligand_chain: str = "L"            # Chain ID for ligand
    metal_resnum: int = 1
    ligand_resnum: int = 1

    # Design parameters
    contig_ranges: List[str] = field(default_factory=lambda: ["100-120", "110-130", "130-150"])
    cfg_scales: List[float] = field(default_factory=lambda: [1.5, 2.0, 2.5])

    # Filter thresholds (can be adjusted per complex)
    filter_plddt: float = 0.80
    filter_ptm: float = 0.80
    filter_pae: float = 5.0

    def to_dict(self):
        return {
            "metal": self.metal.to_dict(),
            "ligand": self.ligand.to_dict(),
            "pdb_path": self.pdb_path,
            "metal_chain": self.metal_chain,
            "ligand_chain": self.ligand_chain,
            "contig_ranges": self.contig_ranges,
            "cfg_scales": self.cfg_scales,
            "filters": {
                "plddt": self.filter_plddt,
                "ptm": self.filter_ptm,
                "pae": self.filter_pae,
            }
        }


# ============================================================================
# PRE-DEFINED METAL-LIGAND TEMPLATES
# ============================================================================

METAL_TEMPLATES = {
    # Lanthanides (hard acids, prefer O donors)
    "TB": MetalProperties(
        code="TB", name="Terbium", charge=3,
        coordination_number=9, preferred_donors=["O"],
        ionic_radius=1.04, hsab_class="hard"
    ),
    "GD": MetalProperties(
        code="GD", name="Gadolinium", charge=3,
        coordination_number=9, preferred_donors=["O"],
        ionic_radius=1.05, hsab_class="hard"
    ),
    "EU": MetalProperties(
        code="EU", name="Europium", charge=3,
        coordination_number=9, preferred_donors=["O"],
        ionic_radius=1.07, hsab_class="hard"
    ),
    "LA": MetalProperties(
        code="LA", name="Lanthanum", charge=3,
        coordination_number=9, preferred_donors=["O"],
        ionic_radius=1.16, hsab_class="hard"
    ),

    # Alkaline earth (hard acids)
    "CA": MetalProperties(
        code="CA", name="Calcium", charge=2,
        coordination_number=7, preferred_donors=["O"],
        ionic_radius=1.00, hsab_class="hard"
    ),
    "MG": MetalProperties(
        code="MG", name="Magnesium", charge=2,
        coordination_number=6, preferred_donors=["O"],
        ionic_radius=0.72, hsab_class="hard"
    ),

    # Transition metals (borderline/soft)
    "ZN": MetalProperties(
        code="ZN", name="Zinc", charge=2,
        coordination_number=4, preferred_donors=["N", "S", "O"],
        ionic_radius=0.74, hsab_class="borderline"
    ),
    "FE": MetalProperties(
        code="FE", name="Iron", charge=2,
        coordination_number=6, preferred_donors=["N", "O", "S"],
        ionic_radius=0.78, hsab_class="borderline"
    ),
    "FE3": MetalProperties(
        code="FE", name="Iron(III)", charge=3,
        coordination_number=6, preferred_donors=["O", "N"],
        ionic_radius=0.65, hsab_class="hard"
    ),
    "CU": MetalProperties(
        code="CU", name="Copper", charge=2,
        coordination_number=4, preferred_donors=["N", "S", "O"],
        ionic_radius=0.73, hsab_class="borderline"
    ),
    "MN": MetalProperties(
        code="MN", name="Manganese", charge=2,
        coordination_number=6, preferred_donors=["O", "N"],
        ionic_radius=0.83, hsab_class="borderline"
    ),
    "CO": MetalProperties(
        code="CO", name="Cobalt", charge=2,
        coordination_number=6, preferred_donors=["N", "O", "S"],
        ionic_radius=0.75, hsab_class="borderline"
    ),
    "NI": MetalProperties(
        code="NI", name="Nickel", charge=2,
        coordination_number=6, preferred_donors=["N", "S", "O"],
        ionic_radius=0.69, hsab_class="borderline"
    ),
}

LIGAND_TEMPLATES = {
    # Carboxylate ligands
    "CIT": LigandProperties(
        code="CIT", name="Citrate",
        coordinating_atoms=["O1", "O3", "O5", "O7"],
        hbond_acceptors=["O1", "O2", "O3", "O4", "O5", "O6"],
        hbond_donors=["O7"],
        smiles="OC(CC(O)(C([O-])=O)CC([O-])=O)([O-])=O",
        metal_donors=4
    ),
    "EDT": LigandProperties(
        code="EDT", name="EDTA",
        coordinating_atoms=["O1", "O3", "O5", "O7", "N1", "N2"],
        hbond_acceptors=["O1", "O2", "O3", "O4", "O5", "O6", "O7", "O8"],
        hbond_donors=[],
        smiles="C(CN(CC(=O)[O-])CC(=O)[O-])N(CC(=O)[O-])CC(=O)[O-]",
        metal_donors=6
    ),
    "NTA": LigandProperties(
        code="NTA", name="Nitrilotriacetate",
        coordinating_atoms=["O1", "O3", "O5", "N1"],
        hbond_acceptors=["O1", "O2", "O3", "O4", "O5", "O6"],
        hbond_donors=[],
        smiles="C(C(=O)[O-])N(CC(=O)[O-])CC(=O)[O-]",
        metal_donors=4
    ),

    # Cofactors
    "PQQ": LigandProperties(
        code="PQQ", name="Pyrroloquinoline quinone",
        coordinating_atoms=["O1", "O2", "O3"],
        hbond_acceptors=["O1", "O2", "O3", "O4", "O5"],
        hbond_donors=["N1"],
        smiles="C1=C2C(=CC(=O)C(=O)C2=NC3=C1C(=O)C(=CC3=O)C(=O)O)C(=O)O",
        metal_donors=3
    ),
    "HEM": LigandProperties(
        code="HEM", name="Heme",
        coordinating_atoms=["NA", "NB", "NC", "ND"],
        hbond_acceptors=["O1A", "O2A", "O1D", "O2D"],
        hbond_donors=[],
        smiles="",  # Complex - use PDB
        metal_donors=4
    ),

    # Simple ligands
    "ACE": LigandProperties(
        code="ACE", name="Acetate",
        coordinating_atoms=["O1", "O2"],
        hbond_acceptors=["O1", "O2"],
        hbond_donors=[],
        smiles="CC([O-])=O",
        metal_donors=2
    ),
    "SO4": LigandProperties(
        code="SO4", name="Sulfate",
        coordinating_atoms=["O1", "O2", "O3", "O4"],
        hbond_acceptors=["O1", "O2", "O3", "O4"],
        hbond_donors=[],
        smiles="[O-]S([O-])(=O)=O",
        metal_donors=2
    ),
    "PO4": LigandProperties(
        code="PO4", name="Phosphate",
        coordinating_atoms=["O1", "O2", "O3", "O4"],
        hbond_acceptors=["O1", "O2", "O3", "O4"],
        hbond_donors=[],
        smiles="[O-]P([O-])([O-])=O",
        metal_donors=2
    ),
}


# ============================================================================
# PIPELINE FUNCTIONS
# ============================================================================

def create_rfd3_payload(pdb_content: str, complex_config: MetalLigandComplex,
                        contig: str, cfg_scale: float, num_designs: int = 1) -> Dict:
    """Create RFD3 payload based on metal-ligand complex configuration."""

    metal_key = f"{complex_config.metal_chain}{complex_config.metal_resnum}"
    ligand_key = f"{complex_config.ligand_chain}{complex_config.ligand_resnum}"

    # Build H-bond conditioning from ligand properties
    hbond_acceptors = ",".join(complex_config.ligand.hbond_acceptors)
    hbond_donors = ",".join(complex_config.ligand.hbond_donors)

    payload = {
        "input": {
            "task": "rfd3",
            "pdb_content": pdb_content,
            "contig": contig,
            "ligand": f"{complex_config.ligand.code},{complex_config.metal.code}",
            "num_designs": num_designs,
            "select_fixed_atoms": {
                metal_key: "all",
                ligand_key: "all",
            },
            "select_buried": {metal_key: "all"},
            "use_classifier_free_guidance": True,
            "cfg_scale": cfg_scale,
            "infer_ori_strategy": "com",
        }
    }

    # Add H-bond conditioning if ligand has acceptors/donors
    if hbond_acceptors:
        payload["input"]["select_hbond_acceptor"] = {ligand_key: hbond_acceptors}
    if hbond_donors:
        payload["input"]["select_hbond_donor"] = {ligand_key: hbond_donors}

    return payload


def run_rfd3_scaffolding(pdb_content: str, complex_config: MetalLigandComplex,
                         contig: str, cfg_scale: float, num_designs: int = 1) -> List[str]:
    """Stage 1: Generate scaffold backbones."""
    payload = create_rfd3_payload(pdb_content, complex_config, contig, cfg_scale, num_designs)

    try:
        response = requests.post(API_URL, json=payload, timeout=1200)
        result = response.json()

        if result.get("status", "").upper() == "COMPLETED":
            designs = result.get("output", {}).get("result", {}).get("designs", [])
            return [d.get("content", d) if isinstance(d, dict) else d for d in designs]
    except Exception as e:
        print(f"    RFD3 error: {e}")
    return []


def run_mpnn_sequence_design(pdb_content: str, num_seqs: int = 8, temperature: float = 0.2) -> List[str]:
    """Stage 2: Generate sequences for backbone."""
    payload = {
        "input": {
            "task": "mpnn",
            "pdb_content": pdb_content,
            "num_seqs": num_seqs,
            "ligand_mpnn_use_atom_context": 1,
            "temperature": temperature,
        }
    }

    sequences = []
    try:
        response = requests.post(API_URL, json=payload, timeout=300)
        result = response.json()

        if result.get("status", "").upper() == "COMPLETED":
            seq_data = result.get("output", {}).get("result", {}).get("sequences", [])
            for item in seq_data:
                if isinstance(item, dict) and "content" in item:
                    lines = item["content"].strip().split("\n")
                    for i in range(0, len(lines), 2):
                        if i + 1 < len(lines) and lines[i].startswith(">"):
                            sequences.append(lines[i + 1])
                elif isinstance(item, str):
                    sequences.append(item)
    except Exception as e:
        print(f"    MPNN error: {e}")
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

    try:
        response = requests.post(API_URL, json=payload, timeout=300)
        result = response.json()

        if result.get("status", "").upper() == "COMPLETED":
            confidences = result.get("output", {}).get("result", {}).get("confidences", {})
            return {
                "plddt": confidences.get("mean_plddt", 0),
                "ptm": confidences.get("ptm", 0),
                "pae": confidences.get("overall_pae", 999),
            }
    except Exception as e:
        print(f"    RF3 error: {e}")
    return {"plddt": 0, "ptm": 0, "pae": 999}


def analyze_sequence_for_metal(sequence: str, metal: MetalProperties) -> Dict[str, Any]:
    """Analyze sequence composition based on metal preferences."""
    length = len(sequence)

    # Count potential coordinating residues based on metal preferences
    donor_counts = {}

    if "O" in metal.preferred_donors:
        donor_counts["carboxylate"] = sequence.count('E') + sequence.count('D')
        donor_counts["hydroxyl"] = sequence.count('S') + sequence.count('T')

    if "N" in metal.preferred_donors:
        donor_counts["imidazole"] = sequence.count('H')
        donor_counts["amine"] = sequence.count('K') + sequence.count('R')

    if "S" in metal.preferred_donors:
        donor_counts["thiol"] = sequence.count('C')
        donor_counts["thioether"] = sequence.count('M')

    # Calculate binding potential score
    total_donors = sum(donor_counts.values())
    target_donors = metal.coordination_number - 4  # Assume ligand provides ~4

    score = min(100, int(100 * total_donors / max(target_donors * 2, 1)))

    # Alanine check
    ala_pct = 100 * sequence.count('A') / length

    return {
        "length": length,
        "donor_counts": donor_counts,
        "total_donors": total_donors,
        "binding_score": score,
        "alanine_pct": round(ala_pct, 1),
        "ala_flag": ala_pct > 25,
    }


def passes_filter(metrics: Dict, complex_config: MetalLigandComplex) -> bool:
    """Check if design passes filter thresholds."""
    return (
        metrics.get("plddt", 0) >= complex_config.filter_plddt and
        metrics.get("ptm", 0) >= complex_config.filter_ptm and
        metrics.get("pae", 999) <= complex_config.filter_pae
    )


# ============================================================================
# MAIN PIPELINE
# ============================================================================

def run_parameter_sweep(complex_config: MetalLigandComplex, output_dir: Path,
                        trials_per_config: int = 10, seqs_per_backbone: int = 4) -> Dict:
    """Run parameter sweep to find best configuration."""
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(complex_config.pdb_path, 'r') as f:
        pdb_content = f.read()

    results = {}

    for contig in complex_config.contig_ranges:
        for cfg_scale in complex_config.cfg_scales:
            config_name = f"{contig}_cfg{cfg_scale}"

            print(f"\n{'='*60}")
            print(f"Testing: {config_name}")
            print(f"  Metal: {complex_config.metal.name}, Ligand: {complex_config.ligand.name}")
            print(f"{'='*60}")

            all_metrics = []
            passed = 0

            # Generate backbones
            backbones = run_rfd3_scaffolding(
                pdb_content, complex_config, contig, cfg_scale, num_designs=trials_per_config
            )
            print(f"  Generated {len(backbones)} backbones")

            for i, backbone in enumerate(backbones):
                sequences = run_mpnn_sequence_design(backbone, num_seqs=seqs_per_backbone)

                for j, seq in enumerate(sequences):
                    name = f"{config_name}_{i}_{j}"
                    rf3_metrics = run_rf3_validation(seq, name)
                    seq_analysis = analyze_sequence_for_metal(seq, complex_config.metal)

                    metrics = {**rf3_metrics, **seq_analysis, "name": name, "sequence": seq}
                    all_metrics.append(metrics)

                    if passes_filter(rf3_metrics, complex_config):
                        passed += 1

            if all_metrics:
                plddts = [m["plddt"] for m in all_metrics if m["plddt"] > 0]
                ptms = [m["ptm"] for m in all_metrics if m["ptm"] > 0]

                results[config_name] = {
                    "contig": contig,
                    "cfg_scale": cfg_scale,
                    "total_designs": len(all_metrics),
                    "passed": passed,
                    "pass_rate": passed / len(all_metrics) if all_metrics else 0,
                    "avg_plddt": statistics.mean(plddts) if plddts else 0,
                    "max_plddt": max(plddts) if plddts else 0,
                    "avg_ptm": statistics.mean(ptms) if ptms else 0,
                    "best_designs": sorted(all_metrics, key=lambda x: x["plddt"], reverse=True)[:5],
                }

                print(f"  Pass rate: {passed}/{len(all_metrics)} ({100*passed/len(all_metrics):.1f}%)")
                print(f"  Avg pLDDT: {results[config_name]['avg_plddt']:.3f}")

    # Save results
    save_results = {
        "complex": complex_config.to_dict(),
        "sweep_results": results,
        "best_config": max(results, key=lambda k: results[k]["pass_rate"]) if results else None,
    }

    with open(output_dir / "sweep_results.json", 'w') as f:
        json.dump(save_results, f, indent=2, default=str)

    return results


def run_production(complex_config: MetalLigandComplex, output_dir: Path,
                   contig: str, cfg_scale: float, num_designs: int = 1000,
                   batch_size: int = 10) -> List[Dict]:
    """Production run with specified configuration."""
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(complex_config.pdb_path, 'r') as f:
        pdb_content = f.read()

    all_passing = []
    total_generated = 0

    print(f"\n{'='*60}")
    print(f"PRODUCTION: {complex_config.metal.name}-{complex_config.ligand.name}")
    print(f"Config: {contig}, CFG {cfg_scale}")
    print(f"Target: {num_designs} designs")
    print(f"{'='*60}")

    while total_generated < num_designs:
        backbones = run_rfd3_scaffolding(
            pdb_content, complex_config, contig, cfg_scale, num_designs=batch_size
        )

        for i, backbone in enumerate(backbones):
            sequences = run_mpnn_sequence_design(backbone, num_seqs=8)

            for j, seq in enumerate(sequences):
                total_generated += 1
                name = f"prod_{total_generated}"

                rf3_metrics = run_rf3_validation(seq, name)

                if passes_filter(rf3_metrics, complex_config):
                    seq_analysis = analyze_sequence_for_metal(seq, complex_config.metal)
                    design = {"name": name, "sequence": seq, **rf3_metrics, **seq_analysis}
                    all_passing.append(design)

                if total_generated % 100 == 0:
                    print(f"  Progress: {total_generated}/{num_designs}, Passing: {len(all_passing)}")

    # Save results
    with open(output_dir / "all_passing.fasta", 'w') as f:
        for d in sorted(all_passing, key=lambda x: x["plddt"], reverse=True):
            f.write(f">{d['name']}_pLDDT{d['plddt']:.2f}\n{d['sequence']}\n")

    summary = {
        "complex": complex_config.to_dict(),
        "config": {"contig": contig, "cfg_scale": cfg_scale},
        "total_generated": total_generated,
        "total_passing": len(all_passing),
        "pass_rate": len(all_passing) / total_generated if total_generated > 0 else 0,
        "best_designs": sorted(all_passing, key=lambda x: x["plddt"], reverse=True)[:20],
    }

    with open(output_dir / "production_summary.json", 'w') as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\nComplete: {len(all_passing)}/{total_generated} passing ({100*len(all_passing)/total_generated:.1f}%)")

    return all_passing


# ============================================================================
# CLI
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="Metal-Ligand Binding Protein Design Pipeline")
    parser.add_argument("--mode", choices=["sweep", "production"], required=True)
    parser.add_argument("--metal", type=str, default="TB", help="Metal code (TB, CA, ZN, etc.)")
    parser.add_argument("--ligand", type=str, default="CIT", help="Ligand code (CIT, EDT, etc.)")
    parser.add_argument("--input", type=str, required=True, help="Input PDB with metal-ligand complex")
    parser.add_argument("--output", type=str, default=None, help="Output directory")
    parser.add_argument("--config", type=str, default=None, help="JSON config file (overrides --metal/--ligand)")
    parser.add_argument("--trials", type=int, default=10, help="Trials per config (sweep)")
    parser.add_argument("--num-designs", type=int, default=1000, help="Designs (production)")
    parser.add_argument("--contig", type=str, default="130-150", help="Contig for production")
    parser.add_argument("--cfg", type=float, default=2.5, help="CFG scale for production")
    args = parser.parse_args()

    # Build complex configuration
    if args.config:
        with open(args.config, 'r') as f:
            config_data = json.load(f)
        metal = MetalProperties(**config_data["metal"])
        ligand = LigandProperties(**config_data["ligand"])
        complex_config = MetalLigandComplex(
            metal=metal, ligand=ligand, pdb_path=args.input,
            **{k: v for k, v in config_data.items() if k not in ["metal", "ligand"]}
        )
    else:
        metal = METAL_TEMPLATES.get(args.metal.upper())
        ligand = LIGAND_TEMPLATES.get(args.ligand.upper())

        if not metal:
            print(f"Unknown metal: {args.metal}. Available: {list(METAL_TEMPLATES.keys())}")
            return
        if not ligand:
            print(f"Unknown ligand: {args.ligand}. Available: {list(LIGAND_TEMPLATES.keys())}")
            return

        complex_config = MetalLigandComplex(metal=metal, ligand=ligand, pdb_path=args.input)

    # Set output directory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    if args.output:
        output_dir = Path(args.output)
    else:
        output_dir = Path(args.input).parent / f"{args.mode}_{complex_config.metal.code}_{complex_config.ligand.code}_{timestamp}"

    # Run pipeline
    if args.mode == "sweep":
        results = run_parameter_sweep(complex_config, output_dir, trials_per_config=args.trials)

        if results:
            best = max(results, key=lambda k: results[k]["pass_rate"])
            print(f"\n{'='*60}")
            print(f"BEST CONFIG: {best}")
            print(f"  Pass rate: {results[best]['pass_rate']*100:.1f}%")
            print(f"  Avg pLDDT: {results[best]['avg_plddt']:.3f}")
            print(f"{'='*60}")

    elif args.mode == "production":
        run_production(complex_config, output_dir, args.contig, args.cfg, num_designs=args.num_designs)


if __name__ == "__main__":
    main()
