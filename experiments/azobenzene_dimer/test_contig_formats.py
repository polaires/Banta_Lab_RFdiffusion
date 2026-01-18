#!/usr/bin/env python3
"""Test different contig formats for binder design in RFD3."""

import requests
import urllib.request

API_URL = "http://localhost:8000/runsync"

# Fetch rubredoxin
url = "https://files.rcsb.org/download/1RB9.pdb"
response = urllib.request.urlopen(url)
pdb_lines = []
for line in response.read().decode('utf-8').split('\n'):
    if line.startswith('ATOM'):
        chain = line[21]
        if chain == 'A':
            pdb_lines.append(line)
pdb_lines.append('END')
target_pdb = '\n'.join(pdb_lines)

print(f"Target PDB: {len(pdb_lines)} lines")

# Count residues in target
res_nums = set()
for line in pdb_lines:
    if line.startswith('ATOM'):
        res_nums.add(line[22:26].strip())
print(f"Target residues: {len(res_nums)} (range: {min(int(r) for r in res_nums)}-{max(int(r) for r in res_nums)})")

# Test different contig formats
contig_formats = [
    "A2-52/0 40-60",         # Original format
    "A2-52 40-60",           # Without /0
    "A2-52/0 50",            # Fixed length binder
    "50/A2-52/0 50",         # Binder before and after
    "A2-52/50",              # Different separator
]

for contig in contig_formats:
    print(f"\n{'='*60}")
    print(f"Testing contig: {contig}")
    print('='*60)

    request = {
        "input": {
            "task": "rfd3",
            "contig": contig,
            "pdb_content": target_pdb,
            "num_designs": 1,
        }
    }

    try:
        response = requests.post(API_URL, json=request, timeout=120)
        result = response.json()

        if "output" in result:
            result = result["output"]

        if result.get("status") == "completed":
            designs = result.get("result", {}).get("designs", [])
            if designs:
                pdb_content = designs[0].get("content") or designs[0].get("pdb_content")
                if pdb_content:
                    # Count chains and residues
                    chains = {}
                    for line in pdb_content.split('\n'):
                        if line.startswith('ATOM'):
                            chain = line[21]
                            res_num = line[22:26].strip()
                            if chain not in chains:
                                chains[chain] = set()
                            chains[chain].add(res_num)

                    total_residues = sum(len(r) for r in chains.values())
                    print(f"Status: completed")
                    print(f"Total residues: {total_residues}")
                    for chain, residues in sorted(chains.items()):
                        min_res = min(int(r) for r in residues)
                        max_res = max(int(r) for r in residues)
                        print(f"  Chain {chain}: {len(residues)} residues (range: {min_res}-{max_res})")

                    # Check if binder was generated (total > target)
                    if total_residues > 51:
                        print(f"✓ BINDER GENERATED: {total_residues - 51} new residues")
                    else:
                        print(f"✗ NO BINDER: same as target")
                else:
                    print("No PDB content")
            else:
                print("No designs returned")
        else:
            print(f"Error: {result.get('error')}")
    except Exception as e:
        print(f"Exception: {e}")

print("\n" + "="*60)
print("DONE")
