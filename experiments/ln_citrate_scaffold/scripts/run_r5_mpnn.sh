#!/bin/bash
# Round 5 LigandMPNN: Sequence design on best backbones

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
EXPERIMENT_DIR="$(dirname "$SCRIPT_DIR")"
INPUT_DIR="$EXPERIMENT_DIR/outputs/round_05"
OUTPUT_DIR="$EXPERIMENT_DIR/outputs/round_05_mpnn"

API_URL="http://localhost:8000/runsync"
NUM_SEQS=8

mkdir -p "$OUTPUT_DIR"

echo "============================================================"
echo "Round 5 LigandMPNN: Sequence Design"
echo "============================================================"
echo "Strategy: Metal-aware LigandMPNN, no global bias, omit C"
echo ""

# Check API
echo "Checking API..."
if ! curl -s -X POST $API_URL -H "Content-Type: application/json" -d '{"input": {"task": "health"}}' | grep -q '"healthy":true'; then
    echo "ERROR: API not healthy"
    exit 1
fi
echo "API healthy"
echo ""

# Best backbones from Round 5 analysis (top CN designs)
BEST_BACKBONES=(
    "r5_a_both_fixed_000"   # CN=5: GLU, ASP, HIS
    "r5_b_hbond_011"        # CN=5: GLU, ASN, GLN
    "r5_a_both_fixed_001"   # CN=4: GLU, HIS
    "r5_a_both_fixed_006"   # CN=3: THR, ASN, HIS
    "r5_a_both_fixed_013"   # CN=3: ASP, THR
    "r5_a_both_fixed_014"   # CN=3: ASP, ASN
)

run_mpnn() {
    local BACKBONE=$1
    local PDB_FILE="$INPUT_DIR/${BACKBONE}.pdb"

    if [ ! -f "$PDB_FILE" ]; then
        echo "  ERROR: $PDB_FILE not found"
        return 1
    fi

    PDB_CONTENT=$(cat "$PDB_FILE" | python3 -c "import sys,json; print(json.dumps(sys.stdin.read()))")

    PAYLOAD=$(cat <<EOFPAYLOAD
{
  "input": {
    "task": "mpnn",
    "pdb_content": $PDB_CONTENT,
    "num_seqs": $NUM_SEQS,
    "model_type": "ligand_mpnn",
    "omit_AA": "C",
    "ligand_mpnn_use_atom_context": 1,
    "ligand_mpnn_use_side_chain_context": 1,
    "temperature": 0.1,
    "pack_side_chains": 1
  }
}
EOFPAYLOAD
)

    RESULT=$(curl -s -X POST $API_URL \
        -H "Content-Type: application/json" \
        -d "$PAYLOAD" \
        --max-time 300)

    # Extract FASTA content
    FASTA=$(echo "$RESULT" | python3 -c "
import sys, json
try:
    data = json.load(sys.stdin)
    seqs = data.get('output', {}).get('result', {}).get('sequences', [])
    if seqs and isinstance(seqs[0], dict):
        print(seqs[0].get('content', ''))
except Exception as e:
    print('')
" 2>/dev/null)

    if [ -n "$FASTA" ] && [ ${#FASTA} -gt 10 ]; then
        # Save combined FASTA
        echo "$FASTA" > "$OUTPUT_DIR/${BACKBONE}_seqs.fasta"

        # Parse and save individual sequences
        echo "$FASTA" | python3 -c "
import sys
backbone = '$BACKBONE'
lines = sys.stdin.read().strip().split('\n')
seq_idx = 0
current_header = ''
current_seq = []

for line in lines:
    if line.startswith('>'):
        if current_seq:
            # Count D/E
            seq = ''.join(current_seq)
            d_count = seq.count('D')
            e_count = seq.count('E')
            # Save individual file
            fname = f'$OUTPUT_DIR/{backbone}_seq{seq_idx:02d}.fasta'
            with open(fname, 'w') as f:
                f.write(f'>{backbone}_seq{seq_idx:02d} D={d_count} E={e_count}\\n')
                f.write(seq + '\\n')
            print(f'  seq{seq_idx:02d}: D={d_count}, E={e_count}')
            seq_idx += 1
        current_header = line
        current_seq = []
    else:
        current_seq.append(line.strip())

# Don't forget the last sequence
if current_seq:
    seq = ''.join(current_seq)
    d_count = seq.count('D')
    e_count = seq.count('E')
    fname = f'$OUTPUT_DIR/{backbone}_seq{seq_idx:02d}.fasta'
    with open(fname, 'w') as f:
        f.write(f'>{backbone}_seq{seq_idx:02d} D={d_count} E={e_count}\\n')
        f.write(seq + '\\n')
    print(f'  seq{seq_idx:02d}: D={d_count}, E={e_count}')
"
        return 0
    else
        echo "  ERROR: No sequences returned"
        return 1
    fi
}

# Process each backbone
for BACKBONE in "${BEST_BACKBONES[@]}"; do
    echo "=================================================="
    echo "Backbone: $BACKBONE"
    echo "  Generating $NUM_SEQS sequences..."

    START=$(date +%s)
    run_mpnn "$BACKBONE"
    END=$(date +%s)

    echo "  Done in $((END-START))s"
    echo ""
done

# Summary
echo "============================================================"
echo "LigandMPNN Complete"
echo "============================================================"
N_FASTA=$(ls "$OUTPUT_DIR"/*_seq*.fasta 2>/dev/null | wc -l)
echo "Total sequences: $N_FASTA"
echo "Output: $OUTPUT_DIR"
echo ""
echo "Next: RF3 validation"
