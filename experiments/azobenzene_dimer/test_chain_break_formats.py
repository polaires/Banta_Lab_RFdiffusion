#!/usr/bin/env python3
"""Test different chain break formats for RFD3."""

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

print(f"Target: 51 residues (2-52)")

# Test different chain break formats
formats = [
    "A2-52,/0,40-60",       # Forward slash with commas
    "A2-52/0 40-60",         # Forward slash with space (original)
    "A2-52,0,40-60",         # Just zero with commas
    "A2-52,,40-60",          # Double comma
]

for contig in formats:
    print(f"\n{'='*60}")
    print(f"Testing contig: '{contig}'")
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
                    chains = {}
                    for line in pdb_content.split('\n'):
                        if line.startswith('ATOM'):
                            chain = line[21]
                            res_num = line[22:26].strip()
                            if chain not in chains:
                                chains[chain] = set()
                            chains[chain].add(res_num)

                    total = sum(len(r) for r in chains.values())
                    print(f"Status: completed, Total residues: {total}")
                    for chain, residues in sorted(chains.items()):
                        print(f"  Chain {chain}: {len(residues)} residues")

                    if total > 60:
                        print("SUCCESS: Binder generated!")
                    else:
                        print("NO BINDER (same as target)")
                else:
                    print("No PDB content")
            else:
                print("No designs")
        else:
            error = result.get('error', 'Unknown error')
            # Truncate error for readability
            if len(error) > 100:
                error = error[:100] + "..."
            print(f"Error: {error}")
    except Exception as e:
        print(f"Exception: {e}")

print("\n" + "="*60)
print("DONE")
