"""
Check the actual distances between protein and Dy-TriNOx.
"""

from pathlib import Path
import math

OUTPUT_DIR = Path("G:/Github_local_repo/Banta_Lab_RFdiffusion/experiments/Dy_TriNOx_scaffold/outputs")

def parse_pdb(pdb_path):
    atoms = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21]
                res_num = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms.append({
                    "name": atom_name,
                    "res_name": res_name,
                    "chain": chain,
                    "res_num": res_num,
                    "x": x, "y": y, "z": z
                })
    return atoms

def distance(a1, a2):
    return math.sqrt(
        (a1["x"] - a2["x"])**2 +
        (a1["y"] - a2["y"])**2 +
        (a1["z"] - a2["z"])**2
    )

def main():
    # Check one design
    pdb_file = OUTPUT_DIR / "r1b_medium_000.pdb"
    atoms = parse_pdb(pdb_file)

    # Find Dy
    dy = None
    for a in atoms:
        if a["res_name"] == "DY" or a["name"] == "DY":
            dy = a
            break

    if not dy:
        print("No Dy found!")
        return

    print(f"Dy position: ({dy['x']:.2f}, {dy['y']:.2f}, {dy['z']:.2f})")

    # Find ligand atoms
    ligand_atoms = [a for a in atoms if a["res_name"] == "UNL"]
    print(f"\nLigand (UNL) atoms: {len(ligand_atoms)}")
    if ligand_atoms:
        ligand_center = (
            sum(a["x"] for a in ligand_atoms) / len(ligand_atoms),
            sum(a["y"] for a in ligand_atoms) / len(ligand_atoms),
            sum(a["z"] for a in ligand_atoms) / len(ligand_atoms)
        )
        print(f"Ligand center: ({ligand_center[0]:.2f}, {ligand_center[1]:.2f}, {ligand_center[2]:.2f})")

    # Find protein atoms
    protein_atoms = [a for a in atoms if a["chain"] not in ["X", "L"] and a["res_name"] not in ["DY", "UNL"]]
    print(f"\nProtein atoms: {len(protein_atoms)}")

    if not protein_atoms:
        print("No protein atoms found!")
        return

    # Find closest protein atom to Dy
    min_dist = float('inf')
    closest_atom = None
    for a in protein_atoms:
        d = distance(dy, a)
        if d < min_dist:
            min_dist = d
            closest_atom = a

    print(f"\nClosest protein atom to Dy:")
    print(f"  {closest_atom['res_name']}{closest_atom['res_num']} {closest_atom['name']}")
    print(f"  Distance: {min_dist:.2f} A")
    print(f"  Position: ({closest_atom['x']:.2f}, {closest_atom['y']:.2f}, {closest_atom['z']:.2f})")

    # Find protein center
    protein_center = (
        sum(a["x"] for a in protein_atoms) / len(protein_atoms),
        sum(a["y"] for a in protein_atoms) / len(protein_atoms),
        sum(a["z"] for a in protein_atoms) / len(protein_atoms)
    )
    print(f"\nProtein center of mass: ({protein_center[0]:.2f}, {protein_center[1]:.2f}, {protein_center[2]:.2f})")

    # Distance from Dy to protein center
    dy_to_protein_center = math.sqrt(
        (dy["x"] - protein_center[0])**2 +
        (dy["y"] - protein_center[1])**2 +
        (dy["z"] - protein_center[2])**2
    )
    print(f"Distance from Dy to protein center: {dy_to_protein_center:.2f} A")

    # Check distances of all protein atoms to Dy
    distances = [distance(dy, a) for a in protein_atoms]
    print(f"\nProtein-Dy distance distribution:")
    print(f"  Min: {min(distances):.2f} A")
    print(f"  Max: {max(distances):.2f} A")
    print(f"  Mean: {sum(distances)/len(distances):.2f} A")

    # Count atoms at different distance shells
    within_5 = sum(1 for d in distances if d <= 5)
    within_8 = sum(1 for d in distances if d <= 8)
    within_10 = sum(1 for d in distances if d <= 10)
    print(f"\n  Atoms within 5A: {within_5}")
    print(f"  Atoms within 8A: {within_8}")
    print(f"  Atoms within 10A: {within_10}")

    # Check the open coordination site (apical position)
    # For TriNOx, the open site is typically above the Dy, in +Y direction
    print(f"\n--- APICAL SITE ANALYSIS ---")
    # The Dy has coordination from TriNOx in the XZ plane
    # The open apical site should be above (+Y from input structure)

    # Check if any protein atoms are in the apical direction
    apical_atoms = []
    for a in protein_atoms:
        # Check if atom is roughly above Dy (higher Y, within reasonable XZ)
        dy_to_atom_y = a["y"] - dy["y"]
        dy_to_atom_xz = math.sqrt((a["x"] - dy["x"])**2 + (a["z"] - dy["z"])**2)

        # If atom is above Dy (positive Y direction) and within 10A XZ
        if dy_to_atom_y > 0 and dy_to_atom_xz < 10:
            d = distance(dy, a)
            apical_atoms.append((a, d, dy_to_atom_y))

    if apical_atoms:
        apical_atoms.sort(key=lambda x: x[1])  # Sort by distance
        print(f"Protein atoms above Dy (apical direction):")
        for a, d, y_off in apical_atoms[:10]:
            print(f"  {a['res_name']}{a['res_num']} {a['name']}: {d:.2f}A (Y offset: {y_off:.2f}A)")
    else:
        print("No protein atoms found in apical direction!")
        print("The protein may be positioned to the side of the complex, not above it.")

if __name__ == "__main__":
    main()
