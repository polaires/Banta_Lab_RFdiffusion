#!/bin/bash
# Round 5 RF3: Structure prediction on LigandMPNN sequences

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
EXPERIMENT_DIR="$(dirname "$SCRIPT_DIR")"
INPUT_DIR="$EXPERIMENT_DIR/outputs/round_05_mpnn"
OUTPUT_DIR="$EXPERIMENT_DIR/outputs/round_05_rf3"

API_URL="http://localhost:8000/runsync"
CITRATE_SMILES="C(C(=O)[O-])(CC(=O)[O-])(CC(=O)[O-])O"

mkdir -p "$OUTPUT_DIR"

echo "============================================================"
echo "Round 5 RF3: Structure Prediction"
echo "============================================================"
echo "Validating LigandMPNN sequences with RoseTTAFold3"
echo ""

# Check API
echo "Checking API..."
if ! curl -s -X POST $API_URL -H "Content-Type: application/json" -d '{"input": {"task": "health"}}' | grep -q '"healthy":true'; then
    echo "ERROR: API not healthy"
    exit 1
fi
echo "API healthy"
echo ""

# Process only seq00 and seq01 from each backbone (for speed)
SEQS_PER_BACKBONE=2

run_rf3() {
    local FASTA_FILE=$1
    local NAME=$(basename "$FASTA_FILE" .fasta)

    # Read sequence from FASTA
    SEQUENCE=$(tail -1 "$FASTA_FILE" | tr -d '\n\r')

    if [ -z "$SEQUENCE" ] || [ ${#SEQUENCE} -lt 10 ]; then
        echo "  ERROR: Invalid sequence"
        return 1
    fi

    SEQ_JSON=$(echo "$SEQUENCE" | python3 -c "import sys,json; print(json.dumps(sys.stdin.read().strip()))")

    PAYLOAD=$(cat <<EOFPAYLOAD
{
  "input": {
    "task": "rf3",
    "sequence": $SEQ_JSON,
    "ligand_smiles": "$CITRATE_SMILES",
    "metal": "TB"
  }
}
EOFPAYLOAD
)

    RESULT=$(curl -s -X POST $API_URL \
        -H "Content-Type: application/json" \
        -d "$PAYLOAD" \
        --max-time 600)

    # Extract PDB content (RF3 returns mmCIF in predictions[0]["content"])
    PDB=$(echo "$RESULT" | python3 -c "
import sys, json
try:
    data = json.load(sys.stdin)
    result = data.get('output', {}).get('result', {})
    predictions = result.get('predictions', [])
    if predictions and isinstance(predictions[0], dict):
        pdb = predictions[0].get('content', '')
        print(pdb)
    else:
        pdb = result.get('pdb_content') or result.get('pdb') or ''
        print(pdb)
except:
    print('')
" 2>/dev/null)

    if [ -n "$PDB" ] && [ ${#PDB} -gt 100 ]; then
        echo "$PDB" > "$OUTPUT_DIR/${NAME}_rf3.pdb"
        echo "OK (${#SEQUENCE} res)"
        return 0
    else
        echo "FAIL"
        return 1
    fi
}

# Find all FASTA files and process first N per backbone
BACKBONES=$(ls "$INPUT_DIR"/*_seq00.fasta 2>/dev/null | sed 's/_seq00.fasta//' | xargs -n1 basename)

for BACKBONE in $BACKBONES; do
    echo "=================================================="
    echo "Backbone: $BACKBONE"

    for i in $(seq 0 $((SEQS_PER_BACKBONE - 1))); do
        FASTA="$INPUT_DIR/${BACKBONE}_seq$(printf '%02d' $i).fasta"
        if [ -f "$FASTA" ]; then
            printf "  seq%02d... " $i
            START=$(date +%s)
            run_rf3 "$FASTA"
            END=$(date +%s)
            echo "    ($((END-START))s)"
        fi
    done
    echo ""
done

# Summary
echo "============================================================"
echo "RF3 Validation Complete"
echo "============================================================"
N_PDB=$(ls "$OUTPUT_DIR"/*.pdb 2>/dev/null | wc -l)
echo "Total predictions: $N_PDB"
echo "Output: $OUTPUT_DIR"
echo ""
echo "Next: Analyze predictions for coordination"
