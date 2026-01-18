"""
Quick validation script for symmetric approach designs.
Run from experiments/azobenzene_dimer directory.
"""

import sys
import json
from pathlib import Path

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "backend" / "serverless"))

from topology_validation import validate_dimer_topology

def main():
    symmetric_dir = Path(__file__).parent / "outputs" / "approach_symmetric"

    # Find all dimer PDB files
    dimer_files = list(symmetric_dir.glob("*_dimer.pdb"))

    print(f"Found {len(dimer_files)} dimer designs to validate\n")

    results = []
    valid_count = 0

    for pdb_path in sorted(dimer_files):
        result = validate_dimer_topology(str(pdb_path), ligand_resname='UNL')
        result['name'] = pdb_path.stem.replace('_dimer', '')
        results.append(result)

        status = "VALID" if result.get('valid', False) else "INVALID"
        if result.get('valid', False):
            valid_count += 1

        contacts_a = result.get('contacts_chain_a', 0)
        contacts_b = result.get('contacts_chain_b', 0)
        hull = result.get('hull_overlap', 1.0)
        sep = result.get('separable', False)

        print(f"{pdb_path.stem}: {status}")
        print(f"  Hull overlap: {hull:.1%}, Separable: {sep}")
        print(f"  Contacts A: {contacts_a}, B: {contacts_b}")
        print()

    print("=" * 60)
    print(f"SUMMARY: {valid_count}/{len(results)} designs VALID ({100*valid_count/len(results):.0f}%)")
    print("=" * 60)

    # Sort by total contacts for ranking
    valid_results = [r for r in results if r.get('valid', False)]
    valid_results.sort(
        key=lambda x: x.get('contacts_chain_a', 0) + x.get('contacts_chain_b', 0),
        reverse=True
    )

    print("\nTOP 5 VALID DESIGNS (by contacts):")
    for i, r in enumerate(valid_results[:5], 1):
        total = r.get('contacts_chain_a', 0) + r.get('contacts_chain_b', 0)
        print(f"  {i}. {r['name']}: {total} contacts (A:{r.get('contacts_chain_a')}, B:{r.get('contacts_chain_b')})")

    # Save results to JSON
    output_file = Path(__file__).parent / "outputs" / "symmetric_validation_results.json"
    with open(output_file, 'w') as f:
        json.dump({
            'total_designs': len(results),
            'valid_count': valid_count,
            'valid_percentage': 100 * valid_count / len(results),
            'results': results
        }, f, indent=2, default=str)

    print(f"\nResults saved to: {output_file}")

if __name__ == '__main__':
    main()
