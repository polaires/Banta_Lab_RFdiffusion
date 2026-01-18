#!/usr/bin/env python3
"""Test parametric template coordination."""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
BACKEND_DIR = SCRIPT_DIR.parent.parent / "backend" / "serverless"
sys.path.insert(0, str(BACKEND_DIR))

from lanthanide_templates import generate_parametric_template
from validate_coordination import analyze_coordination

print("=== Testing Parametric Template ===")
print()

# Test with default settings (CN=8, 50% bidentate)
templates = generate_parametric_template(
    metal="TB",
    coordination_number=8,
    num_waters=0,
    bidentate_fraction=0.5,
    randomize=False,
    seed=42,
)

pdb_content = templates[0]

# Save for inspection
with open(SCRIPT_DIR / "outputs" / "debug_parametric.pdb", "w") as f:
    f.write(pdb_content)
print("Saved to: outputs/debug_parametric.pdb")
print()

analysis = analyze_coordination(pdb_content, max_coord_distance=3.5)

print(f"Coordination Number: {analysis.get('coordination_number', 0)}")
print(f"Average Distance: {analysis.get('average_distance', 'N/A')} A")
print()

# Print individual oxygen distances
print("Individual oxygen distances:")
for donor in analysis.get("coordinating_donors", []):
    print(f"  {donor['resname']} {donor['chain']}{donor['resnum']} {donor['atom_name']}: {donor['distance']:.2f} A")

print()
print("Expected: CN around 8, distances around 2.5A")
