#!/usr/bin/env python3
"""Test RFD3 binder design with separate contig and length parameters."""

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

# Count residues in target
res_nums = set()
for line in pdb_lines:
    if line.startswith('ATOM'):
        res_nums.add(line[22:26].strip())
print(f"Target: {len(res_nums)} residues (range: {min(int(r) for r in res_nums)}-{max(int(r) for r in res_nums)})")

# Test with separate contig and length
print("\n" + "="*60)
print("Testing RFD3 with contig='A2-52' and length='40-60'")
print("="*60)

request = {
    "input": {
        "task": "rfd3",
        "contig": "A2-52",
        "length": "40-60",
        "pdb_content": target_pdb,
        "num_designs": 1,
    }
}

try:
    response = requests.post(API_URL, json=request, timeout=180)
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

                # Check if binder was generated
                if total_residues > 51:
                    print(f"\nSUCCESS: Binder generated with {total_residues - 51} new residues")
                else:
                    print(f"\nFAILURE: No binder generated (same as target)")

                # Save for inspection
                with open('binder_test_output.pdb', 'w') as f:
                    f.write(pdb_content)
                print(f"Saved to binder_test_output.pdb")
            else:
                print("No PDB content")
        else:
            print("No designs returned")
    else:
        print(f"Error: {result.get('error')}")
except Exception as e:
    print(f"Exception: {e}")
