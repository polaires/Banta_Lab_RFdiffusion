"""
Run LigandMPNN on R7 scaffold design for sequence design.
Uses ligand-aware design without strong amino acid bias.
"""

import requests
import json
from pathlib import Path

API_URL = "http://localhost:8000/runsync"
BASE_DIR = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\experiments\ln_citrate_scaffold")
SCAFFOLD_PDB = BASE_DIR / "outputs" / "round_07_scaffold" / "r7_single_000.pdb"
OUTPUT_DIR = BASE_DIR / "outputs" / "round_07_scaffold"

def load_pdb_content(pdb_path: Path) -> str:
    with open(pdb_path, 'r') as f:
        return f.read()

pdb_content = load_pdb_content(SCAFFOLD_PDB)

# LigandMPNN config - based on lessons learned:
# - Use ligand context (let model see citrate and TB)
# - Don't use strong amino acid bias (kills foldability)
# - Don't omit cysteine for metal binding
config = {
    "task": "mpnn",
    "pdb_content": pdb_content,
    "model_type": "ligand_mpnn",
    "num_sequences": 8,
    "temperature": 0.1,
    "ligand_mpnn_use_atom_context": 1,
    "ligand_mpnn_use_side_chain_context": 1,
    # No strong bias - let ligand context guide the design
    "omit_AA": "C",  # Omit cysteine only
    "save_stats": True,
}

print("Running LigandMPNN on R7 scaffold...")
print(f"Input: {SCAFFOLD_PDB.name}")

try:
    response = requests.post(
        API_URL,
        json={"input": config},
        timeout=300
    )
    result = response.json()

    print(f"\nStatus: {result.get('status')}")

    output = result.get("output", {})
    res = output.get("result", {})

    sequences = res.get("sequences", [])
    print(f"Number of sequences: {len(sequences)}")

    # Save sequences to FASTA
    fasta_path = OUTPUT_DIR / "r7_single_000_seqs.fasta"
    with open(fasta_path, 'w') as f:
        for i, seq_info in enumerate(sequences):
            if isinstance(seq_info, dict):
                seq = seq_info.get("sequence", "")
                score = seq_info.get("score", 0)
                confidence = seq_info.get("overall_confidence", 0)
                ligand_conf = seq_info.get("ligand_confidence", 0)
                f.write(f">r7_seq{i:02d}_score{score:.2f}_conf{confidence:.2f}_ligconf{ligand_conf:.2f}\n")
                f.write(f"{seq}\n")
            else:
                f.write(f">r7_seq{i:02d}\n")
                f.write(f"{seq_info}\n")

    print(f"Saved: {fasta_path}")

    # Print sequence summaries
    print("\nSequence Summaries:")
    for i, seq_info in enumerate(sequences):
        if isinstance(seq_info, dict):
            seq = seq_info.get("sequence", "")
            # Count amino acids
            ala_count = seq.count('A')
            glu_count = seq.count('E')
            asp_count = seq.count('D')
            total = len(seq)
            print(f"  seq{i:02d}: A={ala_count} ({100*ala_count/total:.0f}%), E+D={glu_count+asp_count} ({100*(glu_count+asp_count)/total:.0f}%)")
        else:
            seq = seq_info
            ala_count = seq.count('A')
            glu_count = seq.count('E')
            asp_count = seq.count('D')
            total = len(seq)
            print(f"  seq{i:02d}: len={total}, A={ala_count} ({100*ala_count/total:.0f}%), E+D={glu_count+asp_count} ({100*(glu_count+asp_count)/total:.0f}%)")

except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
