"""
Validate top v3 designs with LigandMPNN + RF3 structure prediction.

Top designs based on TriNOx-protein contact analysis:
1. v3_nterm_004 - Score 85 (excellent), 85% hydrophobic coverage, 8 aromatic stacking
2. v3_large_004 - Score 75 (good), 94% hydrophobic coverage, 13 aromatic stacking
3. v3_nterm_003 - Score 75 (good), 76% hydrophobic coverage, 6 aromatic stacking
4. v3_cterm_002 - Score 75 (good), 79% hydrophobic coverage, 7 aromatic stacking
"""

import sys
import json
import time
from pathlib import Path

# Add serverless path for imports
serverless_path = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless")
sys.path.insert(0, str(serverless_path))

try:
    import requests
except ImportError:
    print("ERROR: requests module not found. Install with: pip install requests")
    sys.exit(1)

API_URL = "http://localhost:8000/runsync"
OUTPUT_DIR = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\experiments\Dy_TriNOx_scaffold\outputs_v3")
VALIDATION_DIR = OUTPUT_DIR / "validation"
VALIDATION_DIR.mkdir(exist_ok=True)

# Top designs to validate
TOP_DESIGNS = [
    "v3_nterm_004",
    "v3_large_004",
    "v3_nterm_003",
    "v3_cterm_002",
]

TIMEOUT = 600  # 10 minutes


def call_api(payload, timeout=TIMEOUT):
    """Call the local Docker API."""
    try:
        response = requests.post(
            API_URL,
            json={"input": payload},
            timeout=timeout
        )
        return response.json()
    except requests.exceptions.Timeout:
        return {"error": "API timeout"}
    except Exception as e:
        return {"error": str(e)}


def run_ligandmpnn(pdb_content, design_name, num_seqs=4):
    """Run LigandMPNN on the design."""
    print(f"  Running LigandMPNN ({num_seqs} sequences)...")

    payload = {
        "task": "mpnn",
        "pdb_content": pdb_content,
        "num_sequences": num_seqs,  # Correct parameter name
        "temperature": 0.2,
        "model_type": "ligand_mpnn",  # Use LigandMPNN for ligand context
    }

    result = call_api(payload)

    if "error" in result:
        print(f"    ERROR: {result['error']}")
        return None

    output = result.get("output", {})

    # Check for nested result structure
    inner_result = output.get("result", output)

    if output.get("status") != "completed":
        print(f"    FAILED: {output.get('error', output)}")
        return None

    # Sequences can be in FASTA format (from python API) or list (from CLI)
    raw_sequences = inner_result.get("sequences", [])

    # Parse FASTA format if needed
    sequences = []
    for item in raw_sequences:
        if isinstance(item, dict) and "content" in item:
            # FASTA format: parse the content
            fasta_content = item["content"]
            for block in fasta_content.split(">"):
                if not block.strip():
                    continue
                lines = block.strip().split("\n")
                if len(lines) >= 2:
                    header = lines[0]
                    seq = "".join(lines[1:]).replace(" ", "")
                    # Extract score from header if present
                    score = 0.0
                    if "score=" in header.lower():
                        try:
                            score = float(header.split("score=")[1].split()[0].rstrip(","))
                        except:
                            pass
                    sequences.append({"sequence": seq, "score": score, "header": header})
        elif isinstance(item, dict) and "sequence" in item:
            sequences.append(item)
        elif isinstance(item, str):
            sequences.append({"sequence": item, "score": 0.0})

    if sequences:
        print(f"    Generated {len(sequences)} sequences")
        for i, s in enumerate(sequences[:2]):
            print(f"      Seq {i+1}: len={len(s.get('sequence', ''))}")
        return sequences
    else:
        print(f"    No sequences found in output")
        return None


def run_rf3_validation(sequence, ligand_pdb, design_name, seq_idx):
    """Run RF3/AF3 structure prediction to validate the sequence."""
    print(f"    Validating sequence {seq_idx+1}...")

    # RF3 expects fasta + ligand
    payload = {
        "task": "rf3",
        "sequence": sequence,
        "ligand_pdb": ligand_pdb,
        "num_models": 1,
    }

    result = call_api(payload, timeout=900)  # 15 min for structure prediction

    if "error" in result:
        print(f"      ERROR: {result['error']}")
        return None

    output = result.get("output", {})

    # Debug: print the output structure
    print(f"      Output keys: {list(output.keys())}")
    if output.get("status"):
        print(f"      Status: {output.get('status')}")

    # Check different success indicators
    if output.get("status") == "completed":
        inner_result = output.get("result", output)

        # Try to get metrics from different locations
        plddt = 0
        ptm = 0
        pdb_content = ""

        # Check for 'confidences' key (RF3 format)
        if "confidences" in inner_result:
            confidences = inner_result["confidences"]
            if isinstance(confidences, dict):
                # RF3 returns: mean_plddt, overall_plddt, ptm, iptm, ranking_score
                plddt = confidences.get("mean_plddt") or confidences.get("overall_plddt") or confidences.get("plddt", 0)
                ptm = confidences.get("ptm") or confidences.get("pTM", 0)
            elif isinstance(confidences, list) and len(confidences) > 0:
                c = confidences[0]
                if isinstance(c, dict):
                    plddt = c.get("mean_plddt") or c.get("overall_plddt") or c.get("plddt", 0)
                    ptm = c.get("ptm") or c.get("pTM", 0)

        # Check for 'predictions' key
        if "predictions" in inner_result:
            predictions = inner_result["predictions"]
            if isinstance(predictions, list) and len(predictions) > 0:
                pred = predictions[0]
                if isinstance(pred, dict):
                    plddt = plddt or pred.get("plddt", pred.get("avg_plddt", 0))
                    ptm = ptm or pred.get("ptm", pred.get("pTM", 0))
                    pdb_content = pred.get("pdb_content", pred.get("structure", ""))

        # Direct keys
        if not plddt:
            plddt = inner_result.get("plddt") or inner_result.get("avg_plddt", 0)
        if not ptm:
            ptm = inner_result.get("ptm") or inner_result.get("pTM", 0)
        if not pdb_content:
            pdb_content = inner_result.get("pdb_content", "")

        print(f"      pLDDT: {plddt:.1f}, pTM: {ptm:.3f}")
        return {
            "plddt": plddt,
            "ptm": ptm,
            "pdb_content": pdb_content,
        }
    elif output.get("success"):
        plddt = output.get("plddt", 0)
        ptm = output.get("ptm", 0)
        print(f"      pLDDT: {plddt:.1f}, pTM: {ptm:.3f}")
        return output
    else:
        print(f"      FAILED: {output.get('error', output)}")
        return None


def extract_ligand_pdb(pdb_content):
    """Extract just the ligand (DY + UNL) from the PDB."""
    ligand_lines = []
    for line in pdb_content.split('\n'):
        if line.startswith("HETATM"):
            res_name = line[17:20].strip()
            if res_name in ["DY", "UNL"]:
                ligand_lines.append(line)
    return '\n'.join(ligand_lines)


def validate_design(design_name):
    """Full validation pipeline for a single design."""
    print(f"\n{'='*60}")
    print(f"Validating: {design_name}")
    print('='*60)

    pdb_path = OUTPUT_DIR / f"{design_name}.pdb"
    if not pdb_path.exists():
        print(f"  ERROR: PDB file not found: {pdb_path}")
        return None

    with open(pdb_path, 'r') as f:
        pdb_content = f.read()

    # Extract ligand for RF3 validation
    ligand_pdb = extract_ligand_pdb(pdb_content)

    # Step 1: LigandMPNN sequence design
    sequences = run_ligandmpnn(pdb_content, design_name, num_seqs=4)
    if not sequences:
        return None

    # Step 2: Validate top sequences with RF3
    results = []
    for i, seq_data in enumerate(sequences[:2]):  # Validate top 2 sequences
        if isinstance(seq_data, str):
            seq = seq_data
            score = 0.0
        else:
            seq = seq_data.get("sequence", "")
            score = seq_data.get("score", 0)
        print(f"  Sequence {i+1}: score={score:.2f}, len={len(seq)}")

        rf3_result = run_rf3_validation(seq, ligand_pdb, design_name, i)
        if rf3_result:
            results.append({
                "sequence": seq,
                "mpnn_score": score,
                "plddt": rf3_result.get("plddt", 0),
                "ptm": rf3_result.get("ptm", 0),
                "pdb_content": rf3_result.get("pdb_content", ""),
            })

    return {
        "design": design_name,
        "sequences": len(sequences),
        "validated": results,
    }


def main():
    print("="*70)
    print("TOP V3 DESIGN VALIDATION")
    print("LigandMPNN + RF3 Structure Prediction")
    print("="*70)

    # Check API health
    print("\nChecking API health...")
    health = call_api({"task": "health"}, timeout=30)
    if "error" in health:
        print(f"ERROR: API not responding: {health['error']}")
        print("Start the Docker container first!")
        return
    print("API is healthy")

    all_results = []

    for design_name in TOP_DESIGNS:
        result = validate_design(design_name)
        if result:
            all_results.append(result)

            # Save intermediate results
            output_file = VALIDATION_DIR / f"{design_name}_validation.json"
            with open(output_file, 'w') as f:
                json.dump(result, f, indent=2)
            print(f"  Saved: {output_file}")

    # Summary
    print("\n" + "="*70)
    print("VALIDATION SUMMARY")
    print("="*70)

    print(f"\n{'Design':<20} {'Seqs':>6} {'Validated':>10} {'Best pLDDT':>12} {'Best pTM':>10}")
    print("-"*65)

    for r in all_results:
        validated = r.get("validated", [])
        best_plddt = max([v.get("plddt", 0) for v in validated]) if validated else 0
        best_ptm = max([v.get("ptm", 0) for v in validated]) if validated else 0
        print(f"{r['design']:<20} {r['sequences']:>6} {len(validated):>10} {best_plddt:>12.1f} {best_ptm:>10.3f}")

    # Save full results
    summary_file = VALIDATION_DIR / "validation_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nFull results saved to: {summary_file}")


if __name__ == "__main__":
    main()
