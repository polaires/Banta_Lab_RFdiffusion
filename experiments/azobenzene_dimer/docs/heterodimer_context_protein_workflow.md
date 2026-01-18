# Heterodimer Design: Context Protein Workflow

## Overview

This document describes the **Context Protein** approach for designing true heterodimers with a small molecule ligand at the protein-protein interface. This approach fixes the coordinate frame problem that causes independent chain designs to fail.

## The Problem

When designing heterodimers independently:
1. Chain A is designed with ligand at origin (0,0,0)
2. Chain B is designed with ligand at origin (0,0,0) - but in a DIFFERENT coordinate frame
3. Combining them places Chain B far from Chain A's ligand

**Result:** Chain B has 0 contacts with ligand, not a true interface dimer.

## The Solution: Context Protein Design

Design Chain B WITH Chain A + ligand as input context:

```
Step 1: Design Chain A with ligand (standard RFD3)
Step 2: Design Chain B using Chain A + ligand as FIXED INPUT
        RFD3 designs Chain B in the same coordinate frame as Chain A
Step 3: Output already has both chains around the same ligand
```

## RFD3 Contig Syntax

```
contig = "A1-{chain_a_length},/0,{min_len}-{max_len}"
```

- `A1-{len}` = Keep Chain A residues 1 to len FIXED
- `/0` = Chain break (new chain starts)
- `{min}-{max}` = Design new chain with this length range

## API Usage

### Asymmetric RASA Approach
```python
{
    "task": "interface_ligand_design",
    "ligand_smiles": "c1ccc(/N=N\\c2ccccc2)cc1",  # Cis-azobenzene
    "approach": "asymmetric_rasa",
    "chain_length": "50-70",
    "num_designs": 5,
    "seed": 100,
}
```

### Induced Dimerization Approach
```python
{
    "task": "interface_ligand_design",
    "ligand_smiles": "c1ccc(/N=N\\c2ccccc2)cc1",
    "approach": "induced",
    "chain_length": "50-70",
    "num_designs": 5,
    "seed": 100,
}
```

### Joint Approach (Reference)
```python
{
    "task": "interface_ligand_design",
    "ligand_smiles": "c1ccc(/N=N\\c2ccccc2)cc1",
    "approach": "joint",
    "chain_length": "50-70",
    "num_designs": 5,
    "seed": 100,
}
```

## Implementation Details

### Key Code Changes (handler.py)

**Chain B Design with Context:**
```python
# Count residues in Chain A
chain_a_len = _count_chain_residues(pdb_a, "A")

# Build context contig
context_contig = f"A1-{chain_a_len},/0,{min_len}-{max_len}"

# Design Chain B with Chain A + ligand as context
chain_b_input = {
    "task": "rfd3",
    "pdb_content": pdb_a,        # Chain A + ligand as INPUT
    "contig": context_contig,     # Fix A, design new B
    "ligand": "UNL",              # Use existing ligand
    "ori_token": [0.0, 0.0, 0.0], # Center on ligand
    "select_buried": {"UNL": buried_atoms},
    "select_exposed": {"UNL": exposed_atoms},
}
```

### Critical Parameters

| Parameter | Purpose |
|-----------|---------|
| `pdb_content` | Provides Chain A + ligand as context |
| `contig` | Fixes Chain A, specifies new chain length |
| `ligand: "UNL"` | Uses existing ligand from input PDB |
| `ori_token: [0,0,0]` | Centers new chain design on ligand |

## Test Results (2026-01-13)

### 5 Designs Per Approach

| Approach | Avg Affinity | Best Affinity | Success Rate |
|----------|-------------|---------------|--------------|
| **Joint** | -1.33 | **-6.38** | 100% |
| **Asymmetric RASA** | 1.47 | **-6.44** | 100% |
| **Induced** | -1.30 | **-6.73** | 100% |

### Validation Metrics

All successful designs have:
- Both chains with contacts (5+ each)
- No steric clashes
- Chains on opposite sides (angle > 120°)
- True heterodimer (sequence identity < 70%)

## Validation Checklist

After generating designs:

1. **Check contacts:** Both chains should have ≥5 contacts with ligand
2. **Check clashes:** `has_clashes: false` in results
3. **Check affinity:** Negative affinity (< 0 kcal/mol)
4. **Check positioning:** Angle between chains > 120° (opposite sides)
5. **Visual verification:** Open in PyMOL, confirm ligand is at interface

## Analysis Tool

Use `analyze_heterodimer.py` to verify designs:

```bash
cd experiments/azobenzene_dimer
python analyze_heterodimer.py
```

Output shows:
- Chain distances from ligand
- Contact counts per chain
- Angle between chains (dot product analysis)
- PASS/FAIL status for positioning

## Troubleshooting

| Issue | Cause | Solution |
|-------|-------|----------|
| Chain B at 25Å from ligand | Old translation code | Restart Docker to pick up context protein fix |
| `contacts_b: 0` | Chain B not near ligand | Verify `pdb_content` includes Chain A + ligand |
| Positive affinity | Clashes or poor fit | Run more designs, check `ori_token` is set |
| Chains on same side | RFD3 placement variance | Generate more designs, filter by angle |

## File Locations

- **Handler code:** `backend/serverless/handler.py`
  - `_design_asymmetric_rasa_heterodimer()` ~line 3850
  - `_design_induced_heterodimer()` ~line 4060
- **Test script:** `experiments/azobenzene_dimer/test_heterodimer_approaches.py`
- **Analysis tool:** `experiments/azobenzene_dimer/analyze_heterodimer.py`
- **Output:** `experiments/azobenzene_dimer/outputs/heterodimer_test_*/`
