#!/usr/bin/env python3
"""Check what chain labels RFD3 outputs for binder design."""

import json
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
target_pdb = '\n'.join(pdb_lines)  # Full PDB

print(f"Target PDB: {len(pdb_lines)} lines")

# Test RFD3
request = {
    "input": {
        "task": "rfd3",
        "contig": "A2-52/0 40-60",
        "pdb_content": target_pdb,
        "num_designs": 1,
    }
}

print("Calling RFD3...")
response = requests.post(API_URL, json=request, timeout=120)
result = response.json()

if "output" in result:
    result = result["output"]

print(f"Status: {result.get('status')}")

if result.get("status") == "completed":
    designs = result.get("result", {}).get("designs", [])
    if designs:
        pdb_content = designs[0].get("content") or designs[0].get("pdb_content")
        if pdb_content:
            # Check chain labels
            chains = set()
            for line in pdb_content.split('\n'):
                if line.startswith('ATOM'):
                    chain = line[21]
                    chains.add(chain)

            print(f"Output chains: {sorted(chains)}")

            # Count residues per chain
            for chain in sorted(chains):
                res_nums = set()
                for line in pdb_content.split('\n'):
                    if line.startswith('ATOM') and line[21] == chain:
                        res_nums.add(line[22:26].strip())
                print(f"  Chain {chain}: {len(res_nums)} residues")

            # Save sample
            with open('/tmp/rfd3_sample.pdb', 'w') as f:
                f.write(pdb_content)
            print("Saved to /tmp/rfd3_sample.pdb")
        else:
            print("No PDB content in design")
    else:
        print("No designs returned")
else:
    print(f"Error: {result.get('error')}")
