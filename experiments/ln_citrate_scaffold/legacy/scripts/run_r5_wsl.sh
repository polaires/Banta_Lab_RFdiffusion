#!/bin/bash
# Round 5: RFD3 with FIXED TB-citrate complex
# Run from WSL to access Docker API

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
EXPERIMENT_DIR="$(dirname "$SCRIPT_DIR")"
CONFIG_FILE="$EXPERIMENT_DIR/configs/round_05/r5_configs.json"
OUTPUT_DIR="$EXPERIMENT_DIR/outputs/round_05"
INPUT_DIR="$EXPERIMENT_DIR/inputs"

API_URL="http://localhost:8000/runsync"
DESIGNS_PER_CONFIG=16

mkdir -p "$OUTPUT_DIR"

echo "============================================================"
echo "Ln-Citrate Scaffold - Round 5: FIXED TB-Citrate Complex"
echo "============================================================"
echo ""
echo "CRITICAL FIX:"
echo "  - Round 2b had CIT NOT fixed -> drifted 17A from TB"
echo "  - Round 5 fixes BOTH TB AND CIT atoms"
echo "  - This maintains coordination geometry (TB-O: 2.35-2.5A)"
echo ""

# Check API health
echo "Checking API health..."
HEALTH=$(curl -s -X POST $API_URL -H "Content-Type: application/json" -d '{"input": {"task": "health"}}')
if echo "$HEALTH" | grep -q '"healthy":true'; then
    echo "  API healthy"
else
    echo "  ERROR: API not healthy"
    echo "$HEALTH"
    exit 1
fi

# Read input PDB
INPUT_PDB="$INPUT_DIR/citrate_ln_only.pdb"
PDB_CONTENT=$(cat "$INPUT_PDB" | python3 -c "import sys,json; print(json.dumps(sys.stdin.read()))")

echo ""
echo "Input: $INPUT_PDB"
echo "Output: $OUTPUT_DIR"
echo ""

# Run config: r5_a_both_fixed_buried
run_config() {
    local CONFIG_NAME=$1
    local CONTIG=$2
    local LIGAND=$3
    local SELECT_FIXED_ATOMS=$4
    local SELECT_BURIED=$5
    local EXTRA_OPTS=$6

    echo "=================================================="
    echo "Config: $CONFIG_NAME"
    echo "  Contig: $CONTIG"
    echo "  Fixed atoms: $SELECT_FIXED_ATOMS"
    echo ""

    for i in $(seq 0 $((DESIGNS_PER_CONFIG - 1))); do
        printf "  Design %d/%d... " $((i+1)) $DESIGNS_PER_CONFIG

        # Build payload
        PAYLOAD=$(cat <<EOF
{
  "input": {
    "task": "rfd3",
    "pdb_content": $PDB_CONTENT,
    "contig": "$CONTIG",
    "ligand": "$LIGAND",
    "select_fixed_atoms": $SELECT_FIXED_ATOMS,
    "select_buried": $SELECT_BURIED,
    "infer_ori_strategy": "com",
    "use_classifier_free_guidance": true,
    "cfg_scale": 2.0,
    "is_non_loopy": true,
    "plddt_enhanced": true
    $EXTRA_OPTS
  }
}
EOF
)

        # Run API call
        START_TIME=$(date +%s)
        RESULT=$(curl -s -X POST $API_URL \
            -H "Content-Type: application/json" \
            -d "$PAYLOAD" \
            --max-time 600)
        END_TIME=$(date +%s)
        ELAPSED=$((END_TIME - START_TIME))

        # Extract PDB content and save
        PDB_OUT=$(echo "$RESULT" | python3 -c "
import sys, json
try:
    data = json.load(sys.stdin)
    result = data.get('output', {}).get('result', {})
    pdb = result.get('pdb_content') or result.get('pdb') or result.get('structure')
    if not pdb and 'predictions' in result:
        pdb = result['predictions'][0].get('pdb_content')
    print(pdb or '')
except:
    print('')
")

        if [ -n "$PDB_OUT" ] && [ "$PDB_OUT" != "None" ]; then
            OUT_FILE="$OUTPUT_DIR/${CONFIG_NAME}_$(printf '%03d' $i).pdb"
            echo "$PDB_OUT" > "$OUT_FILE"
            echo "saved (${ELAPSED}s)"
        else
            echo "FAILED (${ELAPSED}s)"
            # Debug: show error
            echo "$RESULT" | python3 -c "import sys,json; d=json.load(sys.stdin); print('  Error:', d.get('error', d.get('output',{}).get('error','unknown')))" 2>/dev/null || true
        fi
    done
    echo ""
}

# Config r5_a: Both TB and CIT fixed + buried
# NOTE: Use component names (TB, CIT) not chain IDs (X, L)
run_config "r5_a_both_fixed_buried" \
    "80-100" \
    "CIT,TB" \
    '{"TB": "ALL", "CIT": "ALL"}' \
    '{"TB": "ALL", "CIT": "ALL"}' \
    ""

# Config r5_b: Both fixed + H-bond conditioning
run_config "r5_b_both_fixed_hbond" \
    "90-110" \
    "CIT,TB" \
    '{"TB": "ALL", "CIT": "ALL"}' \
    '{"TB": "ALL"}' \
    ', "select_partially_buried": {"CIT": "O1,O2,O3,O4,O5,O6,O7"}, "select_hbond_acceptor": {"CIT": "O1,O2,O3,O4,O5,O6,O7"}'

echo "============================================================"
echo "Round 5 Complete"
echo "============================================================"
echo ""
ls -la "$OUTPUT_DIR"/*.pdb 2>/dev/null | wc -l | xargs -I{} echo "Total designs: {}"
echo ""
echo "Next: Analyze backbones for TB-citrate distance"
