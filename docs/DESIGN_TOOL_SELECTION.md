# Protein Design Tool Selection Guide

## Quick Reference Card

| I want to design... | Backbone Tool | Sequence Tool | Relaxation | Validation |
|---------------------|---------------|---------------|------------|------------|
| Metal-binding protein | RFD3 | LigandMPNN | Optional | AF2 + Metal validation |
| Small molecule binder | RFD3 | LigandMPNN | Optional | AF2 + Docking |
| DNA/RNA binder | RFD3 | LigandMPNN | FastRelax | AF3 |
| Enzyme with active site | RFD3 (motif) | LigandMPNN | FastRelax | AF2 + Docking |
| Protein-protein binder | RFDiffusion | ProteinMPNN | FastRelax | AF2-Multimer |
| Symmetric oligomer | RFDiffusion | ProteinMPNN | FastRelax | AF2-Multimer |
| Symmetric ligand binder | RFD3 + symmetry | LigandMPNN | Optional | AF2-Multimer |
| Metal-mediated dimer | RFD3 + symmetry | LigandMPNN | Optional | AF2-Multimer |
| General scaffold | RFDiffusion | ProteinMPNN | Optional | AF2 |
| Peptide binder | BindCraft | ProteinMPNN | Built-in | Built-in |

## The One Rule

**If your design involves ANY non-protein atoms (small molecules, metals, nucleotides), use LigandMPNN.**

The performance difference is substantial:
- Metals: 77.5% vs 40.6% sequence recovery
- Small molecules: 63.3% vs 50.5%
- Nucleotides: 50.5% vs 34.0%

## Using the Unified Design API

```python
import requests

API_URL = "http://localhost:8000/runsync"

# Automatic tool selection based on parameters
response = requests.post(API_URL, json={
    "input": {
        "task": "design",
        "metal": "ZN",           # Triggers LigandMPNN + metal presets
        "pdb_content": pdb,
        "num_sequences": 8,
    }
})

# Explicit design type
response = requests.post(API_URL, json={
    "input": {
        "task": "design",
        "design_type": "METAL_MEDIATED_DIMER",
        "metal": "ZN",
        "pdb_content": pdb,
    }
})
```

## Design Types

The system supports these design types:

| Design Type | Description | Sequence Tool |
|-------------|-------------|---------------|
| `METAL_BINDING` | Metal ion coordination (Zn, Fe, Cu, etc.) | LigandMPNN |
| `SMALL_MOLECULE_BINDER` | Drug-like molecule binding | LigandMPNN |
| `DNA_RNA_BINDER` | Nucleic acid binding | LigandMPNN |
| `PROTEIN_PROTEIN_BINDER` | Protein-protein interface | ProteinMPNN |
| `ENZYME_ACTIVE_SITE` | Catalytic site scaffolding | LigandMPNN |
| `SYMMETRIC_OLIGOMER` | Homo-oligomeric assemblies | ProteinMPNN |
| `SYMMETRIC_LIGAND_BINDER` | Symmetric binder with ligand | LigandMPNN |
| `METAL_MEDIATED_DIMER` | Dimer with metal at interface | LigandMPNN |
| `GENERAL_SCAFFOLD` | Novel folds, no specific function | ProteinMPNN |
| `PEPTIDE_BINDER` | Peptide binding (BindCraft-style) | ProteinMPNN |

## Metal-Specific Presets

| Metal | bias_AA | omit_AA | Coordinating Residues |
|-------|---------|---------|----------------------|
| Zinc | A:-2.0,H:2.0,C:1.5,E:1.0,D:1.0 | - | HIS, CYS, GLU, ASP |
| Lanthanide | A:-2.0,E:2.5,D:2.5,N:1.0,Q:1.0,C:-2.0 | C | GLU, ASP, ASN, GLN |
| Iron | A:-2.0,H:2.0,E:1.5,D:1.5,C:1.0 | - | HIS, CYS, GLU, ASP, MET |
| Copper | A:-2.0,H:2.0,C:1.5,M:1.5,E:1.0,D:1.0 | - | HIS, CYS, MET, GLU, ASP |
| Calcium | A:-2.0,E:2.0,D:2.0,N:1.5,Q:1.0 | C | GLU, ASP, ASN |

## The Alanine Problem

Without proper biasing, LigandMPNN can generate 60%+ alanine sequences ("AAAA" patterns). The solution:

```python
bias_AA = "A:-2.0,H:2.0,E:1.0,D:1.0"  # Penalize alanine, favor coordinating residues
```

This reduces alanine from 60% to <10%.

## When FastRelax is Needed

**Use FastRelax:**
- After ProteinMPNN (no sidechain packing)
- Protein-protein interfaces
- RFDiffusion backbones (may have strain)

**Skip FastRelax:**
- Using LigandMPNN with pack_side_chains=True
- Metal coordination sites (may disrupt geometry)
- Quick screening

## Complete Workflow: Metal-Binding Dimer

```
1. BACKBONE GENERATION (RFD3)
   - symmetry: {type: C2}
   - ligand: "ZN"
   - select_partially_buried: {ZN: ALL}
   - ori_token: [0, 0, 0]

2. INTERFACE ANALYSIS
   - Auto-detect coordinating residues
   - Set fixed_positions for A63, B39

3. SEQUENCE DESIGN (LigandMPNN)
   - model_type: "ligand_mpnn"
   - bias_AA: "A:-2.0,H:2.0,E:1.0,D:1.0"
   - fixed_positions: ["A63", "B39"]
   - temperature: 0.1
   - pack_side_chains: True

4. VALIDATION
   - AF2-Multimer: Check dimer folds
   - Metal geometry: Verify coordination distances
   - pLDDT > 70, pAE < 10
```

## API Parameters Reference

| Parameter | Type | Description |
|-----------|------|-------------|
| `task` | str | Must be "design" for unified endpoint |
| `pdb_content` | str | Required. Input PDB content |
| `metal` | str | Metal code (ZN, FE, CU, CA, TB, etc.) |
| `ligand` | str | Ligand CCD code |
| `target_pdb` | str | Target protein for binder design |
| `symmetry` | str | Symmetry type (C2, C3, D2, etc.) |
| `design_type` | str | Force specific design type |
| `temperature` | float | Override temperature (default: 0.1) |
| `num_sequences` | int | Number of sequences (default: 8) |
| `bias_AA` | str | Override AA bias |
| `fixed_positions` | list | Positions to keep fixed |

## Response Format

```json
{
    "status": "completed",
    "design_type": "METAL_BINDING",
    "workflow": {
        "backbone_tool": "rfd3",
        "sequence_tool": "ligand_mpnn",
        "relaxation": "optional",
        "validation": "af2",
        "description": "Metal-binding protein design",
        "rationale": "LigandMPNN 77.5% vs ProteinMPNN 40.6% for metals..."
    },
    "sequences": [...],
    "config_used": {
        "sequence_tool": "ligand_mpnn",
        "temperature": 0.1,
        "bias_AA": "A:-2.0,H:2.0,C:1.5,E:1.0,D:1.0",
        "fixed_positions": ["A63", "B39"]
    }
}
```

## Further Reading

- [LigandMPNN Paper](https://www.nature.com/articles/s41592-024-02530-w) - Nature Methods 2025
- [RFD3 Documentation](./RFD3_DESIGN_INSTRUCTION.md)
- [BindCraft Source](https://github.com/martinpacesa/BindCraft)
