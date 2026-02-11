"""Test metal binding API with Round7b parameters."""
import requests
import json

API_URL = "http://localhost:8000/runsync"

# Read the Round7b input PDB
with open(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\experiments\ln_citrate_scaffold\inputs\tb_citrate_motif_scaffold.pdb", 'r') as f:
    pdb_content = f.read()

print("=" * 70)
print("Testing RFD3 with Round7b parameters")
print("=" * 70)

# Test 1: Direct RFD3 call (like Round7b does)
payload = {
    "input": {
        "task": "rfd3",
        "pdb_content": pdb_content,
        "contig": "110-130",
        "ligand": "CIT,TB",
        "select_fixed_atoms": {
            "X1": "all",  # TB
            "L1": "all",  # Citrate
        },
        "select_buried": {"X1": "all"},
        "select_hbond_acceptor": {"L1": "O1,O2,O3,O4,O5,O6"},
        "use_classifier_free_guidance": True,
        "cfg_scale": 2.5,
        "infer_ori_strategy": "com",
        "num_designs": 1,
    }
}

print("\nPayload:")
print(json.dumps({k: v if k != 'pdb_content' else f'<{len(v)} chars>' for k, v in payload['input'].items()}, indent=2))

try:
    response = requests.post(API_URL, json=payload, timeout=600)
    result = response.json()

    print(f"\nStatus code: {response.status_code}")
    print(f"Response status: {result.get('status')}")

    if result.get('status') == 'COMPLETED':
        output = result.get('output', {})
        res = output.get('result', {})
        designs = res.get('designs', res.get('pdbs', []))
        print(f"Generated {len(designs)} designs")

        if designs:
            design = designs[0]
            content = design.get('content', design) if isinstance(design, dict) else design
            print(f"First design: {len(content)} chars")
            print(content[:500] if content else "No content")
    else:
        print(f"Error: {result}")

except Exception as e:
    print(f"Exception: {e}")
