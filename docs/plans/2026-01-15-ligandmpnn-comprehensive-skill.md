# LigandMPNN Comprehensive Skill Document

> **Purpose**: Complete reference for LigandMPNN protein sequence design based on official documentation from dauparas/LigandMPNN. Use this skill when working on any sequence design task.

---

## 1. Core Concepts

### What is LigandMPNN?
LigandMPNN is a deep learning model for protein sequence design that extends ProteinMPNN with **ligand context awareness**. It can design sequences while considering:
- Small molecules
- Metals
- Nucleotides (DNA/RNA)
- Other non-protein atoms

### Key Features
- **Ligand-aware design**: Considers atomic context from bound ligands
- **Side chain packing**: Built-in conformational optimization
- **Symmetry support**: Homo-oligomer and custom symmetric constraints
- **Multiple model types**: ProteinMPNN, LigandMPNN, SolubleMPNN, Membrane variants
- **Batch processing**: Design multiple sequences and process multiple PDBs

### Model Architecture
- Node features: 128
- Edge features: 128
- Hidden dim: 128
- Encoder layers: 3
- Decoder layers: 3
- k_neighbors: 32 (LigandMPNN) or 48 (ProteinMPNN)

---

## 2. Available Model Types

### ProteinMPNN (Original)
**Purpose**: General protein sequence design without ligand context.

**Checkpoints** (noise levels in Å):
| Checkpoint | Noise Level |
|------------|-------------|
| `proteinmpnn_v_48_002.pt` | 0.02 Å |
| `proteinmpnn_v_48_010.pt` | 0.10 Å |
| `proteinmpnn_v_48_020.pt` | 0.20 Å (default) |
| `proteinmpnn_v_48_030.pt` | 0.30 Å |

### LigandMPNN
**Purpose**: Ligand-aware protein design with atomic context.

**Checkpoints** (num_edges=32, atom_context_num=25):
| Checkpoint | Noise Level |
|------------|-------------|
| `ligandmpnn_v_32_005_25.pt` | 0.05 Å |
| `ligandmpnn_v_32_010_25.pt` | 0.10 Å (default) |
| `ligandmpnn_v_32_020_25.pt` | 0.20 Å |
| `ligandmpnn_v_32_030_25.pt` | 0.30 Å |

### SolubleMPNN
**Purpose**: Trained exclusively on soluble proteins (excludes membrane proteins).

**Checkpoints**:
- `solublempnn_v_48_002.pt`
- `solublempnn_v_48_010.pt`
- `solublempnn_v_48_020.pt` (default)
- `solublempnn_v_48_030.pt`

### Membrane Models
**Purpose**: Design for transmembrane proteins.

**Global Label Model** (`global_label_membrane_mpnn`):
- Binary classification: membrane (1) vs soluble (0)
- Checkpoint: `global_label_membrane_mpnn_v_48_020.pt`

**Per-Residue Label Model** (`per_residue_label_membrane_mpnn`):
- Per-residue labels: buried (2), interface (1), other (0)
- Checkpoint: `per_residue_label_membrane_mpnn_v_48_020.pt`

### Side Chain Packing Model
**Purpose**: Conformational optimization of side chains.
- Checkpoint: `ligandmpnn_sc_v_32_002_16.pt`

---

## 3. Command-Line Arguments Reference

### Basic Parameters
```bash
--seed 111                    # Random seed (0 = random)
--pdb_path ./input.pdb        # Input PDB file
--out_folder ./output/        # Output directory
--temperature 0.1             # Sampling temperature (default: 0.1)
--verbose 1                   # Verbosity level (0 = silent)
```

### Model Selection
```bash
--model_type protein_mpnn                    # Original ProteinMPNN
--model_type ligand_mpnn                     # Ligand-aware (default for ligand)
--model_type soluble_mpnn                    # Soluble proteins only
--model_type global_label_membrane_mpnn      # Membrane with global label
--model_type per_residue_label_membrane_mpnn # Membrane with per-residue labels
```

### Residue Control
```bash
# Fix specific residues (keep original sequence)
--fixed_residues "A12 A13 A14 B2 B25"

# Only redesign specific residues (fix everything else)
--redesigned_residues "A12 A13 A14 B2 B25"

# Only design specific chains
--chains_to_design "A,B,C"

# Only parse specific chains from PDB
--parse_these_chains_only "A,B"
```

### Amino Acid Biasing
```bash
# Global bias (positive = favor, negative = disfavor)
--bias_AA "W:3.0,P:3.0,C:-3.0,A:-1.0"

# Per-residue bias from JSON file
--bias_AA_per_residue ./bias.json
# JSON format: {"A12": {"G": -0.3, "C": -2.0, "H": 0.8}}

# Global amino acid omission
--omit_AA "CDFGHILMNPQRSTVWY"  # Only allow AEK

# Per-residue omission from JSON
--omit_AA_per_residue ./omit.json
# JSON format: {"A12": "PG", "A13": "QST"}
```

### Symmetry Constraints
```bash
# Define symmetric residue groups (| separates groups, , separates residues)
--symmetry_residues "C1,C2,C3|C4,C5|C6,C7"

# Weighting for symmetric positions
--symmetry_weights "0.33,0.33,0.33|0.5,0.5|0.5,0.5"

# Automatic homo-oligomer symmetry (same sequence for identical chains)
--homo_oligomer 1
```

### Batch Processing
```bash
--batch_size 3           # Sequences per batch
--number_of_batches 5    # Total batches (total = batch_size × number_of_batches)
```

### LigandMPNN-Specific
```bash
# Enable/disable ligand atom context
--ligand_mpnn_use_atom_context 1      # Use ligand context (default)
--ligand_mpnn_use_atom_context 0      # Ignore ligand atoms

# Include side chain atoms of fixed residues as context
--ligand_mpnn_use_side_chain_context 1

# Distance cutoff for "ligand-proximal" residue classification
--ligand_mpnn_cutoff_for_score 8.0    # Default: 8.0 Å
```

### Membrane Model Parameters
```bash
# Global label (for global_label_membrane_mpnn)
--global_transmembrane_label 1   # 1 = membrane, 0 = soluble

# Per-residue labels (for per_residue_label_membrane_mpnn)
--transmembrane_buried "A1 A2 A3 A11"      # Hydrophobic core
--transmembrane_interface "A4 A5 A6 A22"   # Lipid-water interface
```

### Side Chain Packing
```bash
--pack_side_chains 1                # Enable side chain packing
--number_of_packs_per_design 4      # Samples per sequence (0 = fast single)
--pack_with_ligand_context 1        # Include ligand in packing
--repack_everything 1               # Repack all residues (0 = only redesigned)
--sc_num_denoising_steps 3          # Denoising iterations
--sc_num_samples 16                 # Mixture samples for likelihood
```

### Output Control
```bash
--file_ending "_xyz"           # Custom suffix for output files
--zero_indexed 1               # Start numbering from 0 (default: 1)
--fasta_seq_separation ":"     # Character between chains in FASTA
--save_stats 1                 # Export design statistics (.pt file)
```

### Multi-PDB Processing
```bash
--pdb_path_multi ./pdb_list.json              # JSON with PDB paths
--fixed_residues_multi ./fixed.json           # Per-PDB fixed residues
--redesigned_residues_multi ./redesign.json   # Per-PDB redesign targets
--omit_AA_per_residue_multi ./omit.json       # Per-PDB omissions
--bias_AA_per_residue_multi ./bias.json       # Per-PDB biases
```

---

## 4. Residue Specification Format

### Basic Format
```
ChainID + ResidueNumber + InsertionCode(optional)

Examples:
A12      # Chain A, residue 12
B25      # Chain B, residue 25
B82A     # Chain B, residue 82, insertion code A
C1-C10   # Chain C, residues 1-10 (range notation in some contexts)
```

### JSON Format Examples

**Per-residue bias:**
```json
{
  "A12": {"G": -0.3, "C": -2.0, "H": 0.8},
  "A13": {"G": -1.3},
  "B5": {"W": 2.0, "Y": 1.5}
}
```

**Per-residue omission:**
```json
{
  "A12": "PG",
  "A13": "QST",
  "B5": "C"
}
```

**Multi-PDB fixed residues:**
```json
{
  "./inputs/protein1.pdb": "A12 A13 A14 B2",
  "./inputs/protein2.pdb": "C1 C2 C3"
}
```

---

## 5. Output Files

### Directory Structure
```
out_folder/
├── seqs/          # FASTA files with designed sequences
├── backbones/     # PDB files with backbone + new sequence names
├── packed/        # Full PDB files with packed side chains
└── stats/         # Design statistics (.pt files if save_stats=1)
```

### FASTA Output Format
```
>name, T=0.1, seed=111, num_res=50, num_ligand_res=10, use_ligand_context=True, ...
ORIGINAL_SEQUENCE
>name, id=1, T=0.1, seed=111, overall_confidence=0.8523, ligand_confidence=0.7891, seq_rec=0.4500
DESIGNED_SEQUENCE_1
>name, id=2, T=0.1, seed=111, overall_confidence=0.8234, ligand_confidence=0.7654, seq_rec=0.4200
DESIGNED_SEQUENCE_2
```

### Confidence Metrics
| Metric | Description |
|--------|-------------|
| `overall_confidence` | Average confidence across redesigned residues (0-1) |
| `ligand_confidence` | Confidence for residues within `cutoff_for_score` of ligand |
| `seq_rec` | Sequence recovery vs original (redesigned positions only) |

### Statistics File (.pt)
```python
{
    "generated_sequences": tensor,      # [num_designs, seq_len]
    "sampling_probs": tensor,           # Softmax probabilities
    "log_probs": tensor,                # Log probabilities
    "decoding_order": tensor,           # Order residues were designed
    "native_sequence": tensor,          # Original sequence
    "mask": tensor,                     # Valid residue mask
    "chain_mask": tensor,               # Designable residue mask
    "seed": int,
    "temperature": float
}
```

---

## 6. Common Workflows

### 6.1 Basic Protein Design
```bash
python run.py \
    --model_type protein_mpnn \
    --seed 111 \
    --pdb_path ./protein.pdb \
    --out_folder ./output/ \
    --temperature 0.1 \
    --batch_size 4 \
    --number_of_batches 8
```

### 6.2 Ligand-Binding Site Design
```bash
python run.py \
    --model_type ligand_mpnn \
    --seed 111 \
    --pdb_path ./protein_with_ligand.pdb \
    --out_folder ./output/ \
    --ligand_mpnn_use_atom_context 1 \
    --ligand_mpnn_cutoff_for_score 6.0 \
    --temperature 0.1 \
    --batch_size 4 \
    --number_of_batches 8
```

### 6.3 Metal-Binding Site Design
```bash
python run.py \
    --model_type ligand_mpnn \
    --seed 111 \
    --pdb_path ./protein_with_metal.pdb \
    --fixed_residues "A10 A15 A20 A25" \
    --bias_AA "D:5.0,E:4.0,H:3.0,C:3.0" \
    --out_folder ./output/ \
    --ligand_mpnn_use_atom_context 1 \
    --temperature 0.1
```

### 6.4 Homo-Oligomer Design
```bash
python run.py \
    --model_type ligand_mpnn \
    --seed 111 \
    --pdb_path ./homo_tetramer.pdb \
    --homo_oligomer 1 \
    --out_folder ./output/ \
    --batch_size 4 \
    --number_of_batches 4
```

### 6.5 Partial Redesign with Fixed Binding Site
```bash
python run.py \
    --model_type ligand_mpnn \
    --seed 111 \
    --pdb_path ./enzyme.pdb \
    --fixed_residues "A45 A67 A89 A123" \
    --ligand_mpnn_use_side_chain_context 1 \
    --out_folder ./output/
```

### 6.6 Design with Side Chain Packing
```bash
python run.py \
    --model_type ligand_mpnn \
    --seed 111 \
    --pdb_path ./protein.pdb \
    --out_folder ./output/ \
    --pack_side_chains 1 \
    --number_of_packs_per_design 4 \
    --pack_with_ligand_context 1
```

### 6.7 Membrane Protein Design
```bash
# Global label approach
python run.py \
    --model_type global_label_membrane_mpnn \
    --seed 111 \
    --pdb_path ./membrane_protein.pdb \
    --global_transmembrane_label 1 \
    --out_folder ./output/

# Per-residue label approach
python run.py \
    --model_type per_residue_label_membrane_mpnn \
    --seed 111 \
    --pdb_path ./membrane_protein.pdb \
    --transmembrane_buried "A1 A2 A3 A4 A5" \
    --transmembrane_interface "A6 A7 A8" \
    --out_folder ./output/
```

### 6.8 Batch Processing Multiple PDBs
```bash
python run.py \
    --pdb_path_multi ./pdb_list.json \
    --fixed_residues_multi ./fixed_residues.json \
    --seed 111 \
    --out_folder ./output/
```

---

## 7. Amino Acid Biasing Guide

### Common Bias Presets

**Metal binding (Zn, Fe, Cu):**
```bash
--bias_AA "H:5.0,C:4.0,D:3.0,E:3.0" --omit_AA "P,G"
```

**Lanthanide binding (Tb, Gd, Eu):**
```bash
--bias_AA "D:6.0,E:4.0,N:1.0,Q:1.0" --omit_AA "C,H"
```

**Small molecule binding (aromatic):**
```bash
--bias_AA "W:3.0,Y:2.0,F:2.0,H:1.5"
```

**Nucleotide binding:**
```bash
--bias_AA "R:4.0,K:3.0,N:2.0,Q:2.0,S:1.5,T:1.5"
```

**Increase solubility:**
```bash
--bias_AA "K:2.0,R:2.0,E:2.0,D:2.0" --omit_AA "C,M,W"
```

**Hydrophobic core:**
```bash
--bias_AA "L:2.0,I:2.0,V:2.0,F:1.5,A:1.0"
```

### Bias Value Interpretation
- **Positive values**: Increase probability (log-odds shift)
- **Negative values**: Decrease probability
- **~3.0**: Strong preference (~20x more likely)
- **~1.0**: Mild preference (~2.7x more likely)
- **~-3.0**: Strong disfavor (~20x less likely)

---

## 8. Temperature Guide

| Temperature | Effect | Use Case |
|-------------|--------|----------|
| 0.05 | Very deterministic | Minimal diversity, highest confidence |
| 0.1 | Low diversity (default) | Standard design, good confidence |
| 0.2 | Moderate diversity | Exploring sequence space |
| 0.5 | High diversity | Maximum exploration |
| 1.0 | Very high diversity | Random-like sampling |

**Recommendation**: Start with 0.1, increase to 0.2-0.3 if more diversity needed.

---

## 9. Best Practices

### Input Preparation
1. **Clean PDB**: Remove waters, alternate conformations unless needed
2. **Include ligands**: HETATM records for metals, small molecules
3. **Chain IDs**: Ensure unique chain identifiers
4. **Occupancy**: Set `--parse_atoms_with_zero_occupancy 1` if needed

### Fixed Residues
1. **Catalytic sites**: Always fix active site residues
2. **Binding contacts**: Fix residues making key ligand contacts
3. **Disulfides**: Fix cysteine pairs forming disulfide bonds
4. **Use side chain context**: Enable `--ligand_mpnn_use_side_chain_context 1`

### Sequence Diversity
1. **Batch size**: Use larger batches for efficiency
2. **Temperature**: Increase for more diversity
3. **Multiple runs**: Different seeds for independent sampling

### Quality Control
1. **Check confidence**: Higher overall_confidence = better
2. **Ligand confidence**: Critical for binding site designs
3. **Sequence recovery**: Lower = more novel design

---

## 10. Scoring Sequences

LigandMPNN can score existing sequences without designing new ones.

### Autoregressive Scoring
Conditional probability: p(AA_i | backbone, AA_1...AA_{i-1})

```bash
python score.py \
    --pdb_path ./protein.pdb \
    --model_type ligand_mpnn \
    --autoregressive_score 1
```

### Single AA Scoring
Marginal probability: p(AA_i | backbone, all other AAs)

```bash
python score.py \
    --pdb_path ./protein.pdb \
    --model_type ligand_mpnn \
    --single_aa_score 1
```

---

## 11. Quick Reference Card

```
MODEL TYPES:
  protein_mpnn                    → Original ProteinMPNN
  ligand_mpnn                     → Ligand-aware design
  soluble_mpnn                    → Soluble proteins only
  global_label_membrane_mpnn      → Membrane (global label)
  per_residue_label_membrane_mpnn → Membrane (per-residue)

RESIDUE FORMAT:
  A12                             → Chain A, residue 12
  B82A                            → Chain B, residue 82, insertion A
  "A12 A13 B2"                    → Space-separated list

COMMON PARAMETERS:
  --seed 111                      → Reproducibility
  --temperature 0.1               → Sampling diversity
  --batch_size 4                  → Sequences per batch
  --number_of_batches 8           → Total batches

LIGAND CONTEXT:
  --ligand_mpnn_use_atom_context 1      → Enable (default)
  --ligand_mpnn_use_side_chain_context 1 → Include fixed SC atoms
  --ligand_mpnn_cutoff_for_score 8.0    → Distance cutoff (Å)

RESIDUE CONTROL:
  --fixed_residues "A10 A15"      → Keep these positions
  --redesigned_residues "A10 A15" → Only design these
  --chains_to_design "A,B"        → Design specific chains

BIASING:
  --bias_AA "H:5.0,D:3.0"         → Favor H, D globally
  --omit_AA "C,M"                 → Never use C, M
  --homo_oligomer 1               → Symmetric homo-oligomer

SIDE CHAIN PACKING:
  --pack_side_chains 1            → Enable packing
  --number_of_packs_per_design 4  → Samples per design
  --pack_with_ligand_context 1    → Include ligand

OUTPUT:
  seqs/        → FASTA files
  backbones/   → PDB with sequences
  packed/      → Full PDB with side chains
  stats/       → Design statistics
```

---

## 12. Troubleshooting

| Issue | Cause | Solution |
|-------|-------|----------|
| No ligand atoms parsed | Missing HETATM records | Verify ligand in PDB file |
| Low confidence scores | High temperature or poor structure | Lower temperature, check input |
| Sequence recovery too high | Temperature too low | Increase to 0.2-0.3 |
| Missing chains in output | Wrong chain parsing | Use `--parse_these_chains_only` |
| Side chains clash | Packing without context | Enable `--pack_with_ligand_context 1` |
| Memory error | Batch too large | Reduce `--batch_size` |

---

## 13. Official Resources

| Resource | URL |
|----------|-----|
| GitHub Repository | https://github.com/dauparas/LigandMPNN |
| Original ProteinMPNN Paper | Science 2022 |
| LigandMPNN Paper | bioRxiv 2023 |
| Model Weights | Via `get_model_params.sh` script |

---

*Document generated: 2026-01-15*
*Source: Official LigandMPNN repository (dauparas/LigandMPNN)*
