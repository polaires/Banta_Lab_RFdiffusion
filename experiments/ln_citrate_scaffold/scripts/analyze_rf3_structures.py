"""
Analyze RF3 predicted structures for metal-binding potential.

RF3 outputs CIF files - we analyze these for:
1. Potential metal coordination sites (E, D, H)
2. Secondary structure content
3. Fold compactness
"""

from pathlib import Path
import re

OUTPUT_DIR = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\experiments\ln_citrate_scaffold\outputs\round_07_scaffold")

# RF3 validation results
rf3_results = {
    "r7_design_1": {"plddt": 0.8646, "ptm": 0.8913, "pae": 2.34},
    "r7_design_7": {"plddt": 0.8276, "ptm": 0.7157, "pae": 5.58},
    "r7_design_4": {"plddt": 0.7977, "ptm": 0.5400, "pae": 12.36},
    "r7_design_8": {"plddt": 0.8139, "ptm": 0.5932, "pae": 10.16},
}

sequences = {
    "r7_design_1": "TLRLTITFAPGDLKVTFFDAETGEKLGTYVGRDAIIAANNELRDAGVWHTEVVARADAPEKAVRDSVPHRTLSIEQIAPNTIVAVVETADPAAFVAEETAEVAALGGSLTYEVL",
    "r7_design_4": "RLRVTITYAPGEKLVRFFDAETELLGTYVGMDAILAANADAAAGKWITEEVALADRRPEKAVRDSVPHRILSLERIAPNTRAVVETDDPAAFAAEETAEVAAMGGSMTWEVL",
    "r7_design_7": "TLRVTITWAPGEKLVTFFDAETKEKLGTYVGRDAILAAKNELTAAGKWHTEEVALADRKPEKAVRESVPHKILSLEQVAPDTVRAVIETDDPDALIAAETAAVAAMGGSMTAERL",
    "r7_design_8": "RLRLTITYAPGDLLVRFFNAETGELLGTFVGRDAILAANEELAAGVWHTEEVARADRPEDAVAETTPHILSLEQVAPNTVVAEVDTADPDALVAEETAAVAALGGSLTYERL",
}


def find_carboxylate_clusters(sequence, window=7):
    """
    Find clusters of carboxylate residues (E, D) that could coordinate metals.

    Clusters are regions with multiple E/D within a window.
    Good metal binding sites often have 3+ carboxylates within ~7 residues.
    """
    clusters = []

    for i in range(len(sequence) - window + 1):
        window_seq = sequence[i:i + window]
        e_count = window_seq.count('E')
        d_count = window_seq.count('D')
        total = e_count + d_count

        if total >= 3:  # At least 3 carboxylates in window
            clusters.append({
                "start": i + 1,  # 1-indexed
                "end": i + window,
                "sequence": window_seq,
                "E_count": e_count,
                "D_count": d_count,
                "total": total,
            })

    # Merge overlapping clusters
    if not clusters:
        return []

    merged = [clusters[0]]
    for cluster in clusters[1:]:
        if cluster["start"] <= merged[-1]["end"]:
            # Overlapping - extend the previous cluster
            merged[-1]["end"] = max(merged[-1]["end"], cluster["end"])
            merged[-1]["total"] = max(merged[-1]["total"], cluster["total"])
        else:
            merged.append(cluster)

    return merged


def analyze_binding_potential(sequence):
    """
    Analyze sequence for metal binding potential.

    Lanthanides prefer:
    - Carboxylate oxygens (E, D) - primary
    - Backbone carbonyls
    - Ser/Thr hydroxyls (less common)
    """
    length = len(sequence)

    # Count binding residues
    e_count = sequence.count('E')
    d_count = sequence.count('D')
    h_count = sequence.count('H')
    s_count = sequence.count('S')
    t_count = sequence.count('T')

    # Find E/D positions
    e_positions = [i + 1 for i, aa in enumerate(sequence) if aa == 'E']
    d_positions = [i + 1 for i, aa in enumerate(sequence) if aa == 'D']

    # Find carboxylate clusters
    clusters = find_carboxylate_clusters(sequence)

    # Score binding potential (0-100)
    score = 0

    # Carboxylate content (max 40)
    carbox_ratio = (e_count + d_count) / length
    if carbox_ratio >= 0.15:
        score += 40
    elif carbox_ratio >= 0.10:
        score += 30
    elif carbox_ratio >= 0.05:
        score += 20

    # Carboxylate clusters (max 30)
    if len(clusters) >= 2:
        score += 30
    elif len(clusters) >= 1:
        score += 20

    # His for additional coordination (max 15)
    if h_count >= 2:
        score += 15
    elif h_count >= 1:
        score += 10

    # Ser/Thr for H-bond support (max 15)
    st_ratio = (s_count + t_count) / length
    if st_ratio >= 0.10:
        score += 15
    elif st_ratio >= 0.05:
        score += 10

    return {
        "binding_potential_score": score,
        "carboxylate_count": e_count + d_count,
        "carboxylate_percent": round(100 * (e_count + d_count) / length, 1),
        "glu_positions": e_positions,
        "asp_positions": d_positions,
        "his_count": h_count,
        "carboxylate_clusters": clusters,
        "cluster_count": len(clusters),
    }


def parse_cif_coordinates(cif_content):
    """Extract CA coordinates from CIF file."""
    ca_coords = []

    # CIF atom site loop
    in_atom_site = False
    coord_indices = {}

    for line in cif_content.split('\n'):
        if '_atom_site.Cartn_x' in line:
            in_atom_site = True
            continue

        if in_atom_site:
            if line.startswith('_atom_site.'):
                # Parse column definition
                field = line.split('.')[1].strip()
                idx = len(coord_indices)
                coord_indices[field] = idx
            elif line.startswith('ATOM') or (line and line[0].isupper()):
                parts = line.split()
                if len(parts) > 5:
                    try:
                        atom_name = parts[coord_indices.get('label_atom_id', 3)]
                        if atom_name == 'CA':
                            x = float(parts[coord_indices.get('Cartn_x', 10)])
                            y = float(parts[coord_indices.get('Cartn_y', 11)])
                            z = float(parts[coord_indices.get('Cartn_z', 12)])
                            ca_coords.append((x, y, z))
                    except (ValueError, IndexError, KeyError):
                        continue
            elif line.startswith('#') or line.startswith('_'):
                in_atom_site = False

    return ca_coords


def calculate_radius_of_gyration(coords):
    """Calculate radius of gyration from CA coordinates."""
    if not coords:
        return None

    # Calculate centroid
    n = len(coords)
    cx = sum(c[0] for c in coords) / n
    cy = sum(c[1] for c in coords) / n
    cz = sum(c[2] for c in coords) / n

    # Calculate RoG
    sq_distances = []
    for x, y, z in coords:
        sq_dist = (x - cx)**2 + (y - cy)**2 + (z - cz)**2
        sq_distances.append(sq_dist)

    rg = (sum(sq_distances) / n) ** 0.5
    return round(rg, 2)


print("=" * 70)
print("RF3 STRUCTURE ANALYSIS")
print("=" * 70)

print("\n[1] BINDING POTENTIAL ANALYSIS")
print("-" * 60)

binding_results = {}
for name, seq in sequences.items():
    analysis = analyze_binding_potential(seq)
    binding_results[name] = analysis

    print(f"\n{name}:")
    print(f"  Binding potential score: {analysis['binding_potential_score']}/100")
    print(f"  Carboxylates (E+D): {analysis['carboxylate_count']} ({analysis['carboxylate_percent']}%)")
    print(f"  Histidines: {analysis['his_count']}")
    print(f"  Carboxylate clusters: {analysis['cluster_count']}")

    if analysis['carboxylate_clusters']:
        for cluster in analysis['carboxylate_clusters'][:2]:  # Show max 2
            print(f"    Region {cluster['start']}-{cluster['end']}: {cluster['total']} E/D")

print("\n[2] STRUCTURAL COMPACTNESS")
print("-" * 60)

for name in sequences.keys():
    cif_path = OUTPUT_DIR / f"{name}_rf3.cif"
    if cif_path.exists():
        with open(cif_path, 'r') as f:
            cif_content = f.read()

        coords = parse_cif_coordinates(cif_content)
        rg = calculate_radius_of_gyration(coords)

        if rg:
            # Expected RoG for globular protein: ~0.395 * N^0.6 (Flory)
            n = len(sequences[name])
            expected_rg = 0.395 * (n ** 0.6)
            compactness = expected_rg / rg if rg > 0 else 0

            print(f"\n{name}:")
            print(f"  Radius of gyration: {rg} A")
            print(f"  Expected (Flory): {expected_rg:.1f} A")
            print(f"  Compactness ratio: {compactness:.2f}")
            if compactness > 0.9:
                print(f"  [GOOD] Compact globular fold")
            elif compactness > 0.7:
                print(f"  [OK] Moderately compact")
            else:
                print(f"  [WARN] Extended or disordered")
        else:
            print(f"\n{name}: Could not parse CA coordinates")
    else:
        print(f"\n{name}: CIF file not found")

print("\n[3] COMBINED RANKING")
print("-" * 60)

# Combine RF3 metrics with binding potential
rankings = []
for name in sequences.keys():
    combined = {
        "name": name,
        "plddt": rf3_results[name]["plddt"],
        "ptm": rf3_results[name]["ptm"],
        "pae": rf3_results[name]["pae"],
        "binding_score": binding_results[name]["binding_potential_score"],
        "carboxylate_count": binding_results[name]["carboxylate_count"],
    }

    # Combined score: RF3 confidence * binding potential
    combined["final_score"] = (
        combined["plddt"] * 0.4 +
        combined["ptm"] * 0.3 +
        (100 - combined["pae"]) / 100 * 0.1 +
        combined["binding_score"] / 100 * 0.2
    )
    rankings.append(combined)

rankings.sort(key=lambda x: x["final_score"], reverse=True)

print(f"{'Rank':<5} {'Design':<15} {'pLDDT':>7} {'pTM':>7} {'Bind':>6} {'Final':>7}")
print("-" * 55)
for i, r in enumerate(rankings, 1):
    print(f"{i:<5} {r['name']:<15} {r['plddt']:>7.3f} {r['ptm']:>7.3f} {r['binding_score']:>6} {r['final_score']:>7.3f}")

print("\n[4] FINAL RECOMMENDATION")
print("-" * 60)

best = rankings[0]
print(f"BEST DESIGN FOR AF3: {best['name']}")
print(f"  Final Score: {best['final_score']:.4f}")
print(f"  RF3 pLDDT: {best['plddt']:.4f} (fold confidence)")
print(f"  RF3 pTM: {best['ptm']:.4f} (template matching)")
print(f"  Binding Potential: {best['binding_score']}/100")
print(f"  Carboxylates: {best['carboxylate_count']} E+D residues")

print("\nFor AF3 submission:")
print(f"  Sequence: {sequences[best['name']]}")

print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
