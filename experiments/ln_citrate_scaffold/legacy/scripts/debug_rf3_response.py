"""
Debug RF3 response structure.
"""

import requests
import json

API_URL = "http://localhost:8000/runsync"

seq = "TLRLTITFAPGDLKVTFFDAETGEKLGTYVGRDAIIAANNELRDAGVWHTEVVARADAPEKAVRDSVPHRTLSIEQIAPNTIVAVVETADPAAFVAEETAEVAALGGSLTYEVL"

payload = {
    "input": {
        "task": "rf3",
        "sequence": seq,
        "name": "debug_test"
    }
}

print("Sending RF3 request...")
response = requests.post(API_URL, json=payload, timeout=300)
result = response.json()

# Print full structure without the large content fields
def summarize_structure(obj, depth=0, max_depth=4):
    indent = "  " * depth
    if depth > max_depth:
        return f"{indent}..."

    if isinstance(obj, dict):
        lines = [f"{indent}{{"]
        for k, v in obj.items():
            if k in ['content', 'pdb_content', 'cif_content', 'plddt'] and isinstance(v, str) and len(v) > 100:
                lines.append(f"{indent}  '{k}': <string len={len(v)}>")
            elif k == 'plddt' and isinstance(v, list):
                lines.append(f"{indent}  '{k}': <list of {len(v)} floats, mean={sum(v)/len(v):.4f}>")
            else:
                val_str = summarize_structure(v, depth + 1, max_depth)
                if '\n' in val_str:
                    lines.append(f"{indent}  '{k}':")
                    lines.append(val_str)
                else:
                    lines.append(f"{indent}  '{k}': {val_str}")
        lines.append(f"{indent}}}")
        return '\n'.join(lines)
    elif isinstance(obj, list):
        if len(obj) == 0:
            return "[]"
        elif len(obj) > 5:
            return f"[{len(obj)} items]"
        else:
            lines = [f"{indent}["]
            for item in obj:
                lines.append(summarize_structure(item, depth + 1, max_depth))
            lines.append(f"{indent}]")
            return '\n'.join(lines)
    elif isinstance(obj, str) and len(obj) > 100:
        return f"<string len={len(obj)}>"
    else:
        return repr(obj)

print("\n=== RESPONSE STRUCTURE ===")
print(summarize_structure(result))

# Also print the actual keys at each level
print("\n=== KEY PATHS ===")
print(f"result.keys(): {list(result.keys())}")
if 'output' in result:
    print(f"result['output'].keys(): {list(result.get('output', {}).keys())}")
    if 'result' in result.get('output', {}):
        res = result['output']['result']
        print(f"result['output']['result'].keys(): {list(res.keys())}")
        if 'predictions' in res:
            pred = res['predictions']
            print(f"predictions length: {len(pred)}")
            if pred:
                print(f"predictions[0].keys(): {list(pred[0].keys())}")
                # Print actual metric values
                p0 = pred[0]
                for key in ['mean_plddt', 'overall_plddt', 'ptm', 'iptm', 'overall_pae', 'overall_pde']:
                    if key in p0:
                        print(f"  {key}: {p0[key]}")
