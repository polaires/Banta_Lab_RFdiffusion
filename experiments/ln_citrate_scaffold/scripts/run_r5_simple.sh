#!/bin/bash
# Round 5: Simple RFD3 runner with FIXED TB-citrate complex

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
EXPERIMENT_DIR="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="$EXPERIMENT_DIR/outputs/round_05"
INPUT_PDB="$EXPERIMENT_DIR/inputs/citrate_ln_only.pdb"

API_URL="http://localhost:8000/runsync"
DESIGNS_PER_CONFIG=16

mkdir -p "$OUTPUT_DIR"

echo "============================================================"
echo "Round 5: FIXED TB-Citrate Complex"
echo "============================================================"
echo "TB-CIT coordination FIXED - both atoms fixed in place"
echo ""

# Check API
echo "Checking API..."
if ! curl -s -X POST $API_URL -H "Content-Type: application/json" -d '{"input": {"task": "health"}}' | grep -q '"healthy":true'; then
    echo "ERROR: API not healthy"
    exit 1
fi
echo "API healthy"
echo ""

# Prepare PDB content
PDB_JSON=$(cat "$INPUT_PDB" | python3 -c "import sys,json; print(json.dumps(sys.stdin.read()))")

run_design() {
    local CONFIG=$1
    local CONTIG=$2
    local IDX=$3
    local EXTRA=$4

    PAYLOAD=$(cat <<EOFPAYLOAD
{
  "input": {
    "task": "rfd3",
    "pdb_content": $PDB_JSON,
    "contig": "$CONTIG",
    "ligand": "CIT,TB",
    "select_fixed_atoms": {"TB": "ALL", "CIT": "ALL"},
    "select_buried": {"TB": "ALL", "CIT": "ALL"},
    "infer_ori_strategy": "com",
    "use_classifier_free_guidance": true,
    "cfg_scale": 2.0
    $EXTRA
  }
}
EOFPAYLOAD
)

    RESULT=$(curl -s -X POST $API_URL \
        -H "Content-Type: application/json" \
        -d "$PAYLOAD" \
        --max-time 300)

    # Extract PDB
    PDB_OUT=$(echo "$RESULT" | python3 -c "
import sys, json
try:
    data = json.load(sys.stdin)
    designs = data.get('output', {}).get('result', {}).get('designs', [])
    if designs:
        print(designs[0].get('content', ''))
except:
    pass
" 2>/dev/null)

    if [ -n "$PDB_OUT" ] && [ ${#PDB_OUT} -gt 100 ]; then
        OUT_FILE="$OUTPUT_DIR/${CONFIG}_$(printf '%03d' $IDX).pdb"
        echo "$PDB_OUT" > "$OUT_FILE"
        echo "OK"
        return 0
    else
        echo "FAIL"
        return 1
    fi
}

# Config A: Basic fixed + buried
echo "Config: r5_a_both_fixed (80-100 residues)"
for i in $(seq 0 $((DESIGNS_PER_CONFIG - 1))); do
    printf "  Design %d/%d... " $((i+1)) $DESIGNS_PER_CONFIG
    run_design "r5_a_both_fixed" "80-100" $i ""
done
echo ""

# Config B: With H-bond conditioning
echo "Config: r5_b_hbond (90-110 residues)"
for i in $(seq 0 $((DESIGNS_PER_CONFIG - 1))); do
    printf "  Design %d/%d... " $((i+1)) $DESIGNS_PER_CONFIG
    run_design "r5_b_hbond" "90-110" $i ', "select_hbond_acceptor": {"CIT": "O1,O2,O3,O4,O5,O6,O7"}'
done
echo ""

# Summary
echo "============================================================"
echo "Round 5 Complete"
echo "============================================================"
N_DESIGNS=$(ls "$OUTPUT_DIR"/*.pdb 2>/dev/null | wc -l)
echo "Total designs: $N_DESIGNS"
echo "Output: $OUTPUT_DIR"
echo ""
echo "Next: Verify TB-CIT distances are preserved"
