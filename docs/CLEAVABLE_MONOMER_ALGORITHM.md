# Cleavable Monomer Algorithm for Interface Ligand Design

## Concept

Instead of designing a dimer with ligand at interface (which RFD3 struggles with), 
we design a **monomer with ligand in the CENTER**, then **cleave** the protein at 
a loop region to create a dimer where the ligand is naturally at the interface.

```
DESIGN PHASE:                          CLEAVAGE PHASE:

     ┌──────────────┐                   Chain A    Chain B
     │   N-term     │                   ────┐      ┌────
     │      ↓       │                       │      │
     │  ┌───────┐   │                   ┌───┘      └───┐
     │  │  LIG  │   │       ──►         │    [LIG]    │
     │  └───────┘   │                   └───┐      ┌───┘
     │      ↓       │                       │      │
     │   C-term     │                   ────┘      └────
     └──────────────┘                   
                                        Ligand at interface!
   Complete binding pocket              Pocket remains intact
```

## Why This Works

1. **RFD3 designs complete pockets**: Monomer wraps around ligand properly
2. **Binding contacts preserved**: Cleavage doesn't remove contacts, just separates chains
3. **Natural interface**: The cleavage point IS the interface
4. **Strong affinity**: Pocket was designed whole, not as two halves

---

## Algorithm Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  INPUT                                                                      │
│  - ligand_smiles: SMILES string                                             │
│  - target_length: 80-120 residues                                           │
│  - dimer_type: "homo" or "hetero"                                           │
│  - min_chain_length: 25 residues (minimum per chain after cleavage)         │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 1: DESIGN MONOMER WITH CENTRAL LIGAND                                 │
│                                                                             │
│  RFD3 Config:                                                               │
│    - contig: "80-120"                                                       │
│    - ligand: "UNL"                                                          │
│    - select_buried: {UNL: "ALL"}  ← Force ligand to center                  │
│    - ori_token: [0, 0, 0]         ← Center generation on ligand             │
│                                                                             │
│  Output: Monomer with ligand fully buried in center                         │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 2: IDENTIFY CANDIDATE CLEAVAGE SITES                                  │
│                                                                             │
│  For each residue i in protein:                                             │
│    1. Calculate distance from Cα to ligand centroid                         │
│    2. Identify residues where backbone "crosses" near ligand                │
│    3. Check secondary structure (prefer loops, avoid helix/sheet)           │
│    4. Check that both resulting chains have ligand contacts                 │
│                                                                             │
│  Candidate criteria:                                                        │
│    - Distance to ligand: 4-8 Å (close but not clashing)                     │
│    - Secondary structure: loop/coil (not helix/sheet)                       │
│    - Chain A contacts: ≥3 atoms within 5Å of ligand                         │
│    - Chain B contacts: ≥3 atoms within 5Å of ligand                         │
│    - Min chain length: ≥25 residues each                                    │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 3: SCORE CLEAVAGE SITES                                               │
│                                                                             │
│  For each candidate site i:                                                 │
│                                                                             │
│    contacts_A = ligand contacts from residues 1 to i-1                      │
│    contacts_B = ligand contacts from residues i+1 to N                      │
│    symmetry = min(contacts_A, contacts_B) / max(contacts_A, contacts_B)     │
│    ss_penalty = penalty if cutting through helix/sheet                      │
│    length_balance = 1 - abs(i - N/2) / (N/2)  # Prefer middle cuts          │
│                                                                             │
│    score = (contacts_A + contacts_B) * symmetry * length_balance - ss_penalty│
│                                                                             │
│  Select site with highest score                                             │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 4: PERFORM CLEAVAGE                                                   │
│                                                                             │
│  Option A - Simple cleavage (hetero-dimer):                                 │
│    Chain A = residues 1 to i-1                                              │
│    Chain B = residues i+1 to N                                              │
│    (residue i is removed - it's the "linker")                               │
│                                                                             │
│  Option B - Loop removal (cleaner interface):                               │
│    Identify loop boundaries around cleavage site                            │
│    Remove entire loop (3-7 residues typically)                              │
│    Chain A = residues 1 to loop_start-1                                     │
│    Chain B = residues loop_end+1 to N                                       │
│                                                                             │
│  Output: Two-chain PDB with ligand                                          │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 5: VALIDATE DIMER                                                     │
│                                                                             │
│  Structural:                                                                │
│    - No clashes between chains (min distance > 2.5Å)                        │
│    - Both chains have secondary structure (not just loops)                  │
│    - Chain termini not sterically clashing with ligand                      │
│                                                                             │
│  Binding:                                                                   │
│    - GNINA affinity < -5 kcal/mol                                           │
│    - Contacts: ≥3 per chain                                                 │
│    - Contact symmetry > 50%                                                 │
│                                                                             │
│  Optional - Terminus capping:                                               │
│    - Add ACE cap to new N-termini                                           │
│    - Add NME cap to new C-termini                                           │
│    - Prevents charge-charge artifacts                                       │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 6 (OPTIONAL): HOMO-DIMER GENERATION                                   │
│                                                                             │
│  If user wants homo-dimer:                                                  │
│    1. Take the LARGER chain from cleavage                                   │
│    2. Apply C2 symmetry rotation                                            │
│    3. Check that symmetric copy doesn't clash with original                 │
│    4. Ligand must be on the C2 axis                                         │
│                                                                             │
│  Note: True homo-dimer requires ligand with C2 symmetry (like azobenzene)   │
└─────────────────────────────────────────────────────────────────────────────┘

---

## Detailed Algorithms

### Algorithm 1: Find Cleavage Sites

```python
import numpy as np
from Bio.PDB import DSSP
import biotite.structure as struc

def find_cleavage_sites(structure, ligand_atoms, min_chain_length=25):
    """
    Find candidate cleavage sites where the backbone crosses near the ligand.
    
    A good cleavage site:
    1. Is in a loop (not helix/sheet)
    2. Is close to the ligand (backbone crosses the binding pocket)
    3. Results in both chains having ligand contacts
    4. Results in both chains having sufficient length
    """
    
    protein_atoms = structure[structure.atom_name == "CA"]
    n_residues = len(protein_atoms)
    ligand_centroid = ligand_atoms.coord.mean(axis=0)
    
    # Calculate secondary structure
    dssp = compute_dssp(structure)  # Returns 'H', 'E', 'C' for each residue
    
    candidates = []
    
    for i in range(min_chain_length, n_residues - min_chain_length):
        # Check 1: Is this residue in a loop?
        if dssp[i] != 'C':  # Not coil/loop
            continue
        
        # Check 2: Is this residue close to ligand?
        ca_pos = protein_atoms.coord[i]
        dist_to_ligand = np.linalg.norm(ca_pos - ligand_centroid)
        if dist_to_ligand > 12.0:  # Too far from ligand
            continue
        
        # Check 3: Do both resulting chains contact the ligand?
        chain_a_atoms = structure[structure.res_id <= i]
        chain_b_atoms = structure[structure.res_id > i]
        
        contacts_a = count_contacts(chain_a_atoms, ligand_atoms, cutoff=5.0)
        contacts_b = count_contacts(chain_b_atoms, ligand_atoms, cutoff=5.0)
        
        if contacts_a < 3 or contacts_b < 3:
            continue
        
        # This is a valid candidate
        candidates.append({
            'residue_index': i,
            'residue_id': protein_atoms.res_id[i],
            'distance_to_ligand': dist_to_ligand,
            'contacts_a': contacts_a,
            'contacts_b': contacts_b,
            'secondary_structure': dssp[i],
            'chain_a_length': i,
            'chain_b_length': n_residues - i
        })
    
    return candidates


def score_cleavage_site(candidate, n_residues):
    """
    Score a cleavage site. Higher is better.
    """
    contacts_a = candidate['contacts_a']
    contacts_b = candidate['contacts_b']
    i = candidate['residue_index']
    
    # Contact score: more contacts = better
    contact_score = contacts_a + contacts_b
    
    # Symmetry score: balanced contacts = better
    symmetry = min(contacts_a, contacts_b) / max(contacts_a, contacts_b)
    
    # Length balance: prefer cuts near middle
    length_balance = 1 - abs(i - n_residues/2) / (n_residues/2)
    
    # Distance penalty: very close to ligand might cause issues
    dist = candidate['distance_to_ligand']
    if dist < 5.0:
        distance_penalty = (5.0 - dist) * 2
    else:
        distance_penalty = 0
    
    score = contact_score * symmetry * (0.5 + 0.5 * length_balance) - distance_penalty
    
    return score
```

### Algorithm 2: Perform Cleavage

```python
def cleave_protein(structure, cleavage_residue_id, ligand_atoms, 
                   remove_loop=True, loop_window=3):
    """
    Cleave protein at specified residue to create two chains.
    
    Parameters:
    - structure: protein + ligand structure
    - cleavage_residue_id: residue ID where to cut
    - ligand_atoms: ligand atom array
    - remove_loop: if True, remove entire loop around cleavage site
    - loop_window: residues to check on each side for loop removal
    """
    
    protein = structure[struc.filter_amino_acids(structure)]
    
    if remove_loop:
        # Find loop boundaries
        dssp = compute_dssp(structure)
        
        # Expand cleavage to full loop
        loop_start = cleavage_residue_id
        loop_end = cleavage_residue_id
        
        # Expand backward
        for j in range(cleavage_residue_id - 1, 
                       max(0, cleavage_residue_id - loop_window), -1):
            if dssp[j] == 'C':
                loop_start = j
            else:
                break
        
        # Expand forward
        for j in range(cleavage_residue_id + 1,
                       min(len(dssp), cleavage_residue_id + loop_window)):
            if dssp[j] == 'C':
                loop_end = j
            else:
                break
    else:
        loop_start = cleavage_residue_id
        loop_end = cleavage_residue_id
    
    # Create Chain A (before loop)
    chain_a = protein[protein.res_id < loop_start].copy()
    chain_a.chain_id[:] = 'A'
    
    # Create Chain B (after loop)  
    chain_b = protein[protein.res_id > loop_end].copy()
    chain_b.chain_id[:] = 'B'
    # Renumber residues starting from 1
    chain_b.res_id = chain_b.res_id - loop_end
    
    # Combine
    dimer = chain_a + chain_b + ligand_atoms
    
    return dimer, {
        'loop_start': loop_start,
        'loop_end': loop_end,
        'residues_removed': loop_end - loop_start + 1,
        'chain_a_length': len(set(chain_a.res_id)),
        'chain_b_length': len(set(chain_b.res_id))
    }
```

### Algorithm 3: Homo-dimer from Hetero-dimer

```python
def create_homo_dimer(structure, ligand_atoms, chain_to_use='larger'):
    """
    Create a homo-dimer by taking one chain and applying C2 symmetry.
    
    Requirements:
    - Ligand must be on the C2 axis (e.g., azobenzene with N=N along Z)
    - Selected chain must not clash with its symmetric copy
    """
    
    protein = structure[struc.filter_amino_acids(structure)]
    chain_a = protein[protein.chain_id == 'A']
    chain_b = protein[protein.chain_id == 'B']
    
    # Select which chain to duplicate
    if chain_to_use == 'larger':
        template = chain_a if len(chain_a) > len(chain_b) else chain_b
    elif chain_to_use == 'A':
        template = chain_a
    else:
        template = chain_b
    
    template.chain_id[:] = 'A'
    
    # Apply C2 rotation (180° around Z-axis through ligand centroid)
    ligand_center = ligand_atoms.coord.mean(axis=0)
    
    chain_b_homo = template.copy()
    chain_b_homo.chain_id[:] = 'B'
    
    # C2 rotation: reflect through the ligand center in XY plane
    chain_b_homo.coord[:, 0] = 2 * ligand_center[0] - chain_b_homo.coord[:, 0]
    chain_b_homo.coord[:, 1] = 2 * ligand_center[1] - chain_b_homo.coord[:, 1]
    # Z unchanged for C2 around Z-axis
    
    # Check for clashes
    min_dist = compute_min_distance(template, chain_b_homo)
    if min_dist < 2.5:
        return None, "Chains clash after C2 symmetry"
    
    homo_dimer = template + chain_b_homo + ligand_atoms
    
    return homo_dimer, {
        'min_interchain_distance': min_dist,
        'chain_length': len(set(template.res_id))
    }
```

---

## Comparison to Previous Approaches

| Aspect | Post-Processing Symmetry | Cleavable Monomer |
|--------|-------------------------|-------------------|
| Binding pocket | Two incomplete halves | One complete pocket |
| Design intent | Designed as monomer, forced into dimer | Designed as monomer, cleaved into dimer |
| Contact distribution | Often asymmetric | Both chains contact same pocket |
| Expected affinity | -2 to -4 kcal/mol | -5 to -8 kcal/mol (predicted) |
| Pocket integrity | Compromised | Preserved |

---

## Potential Issues and Solutions

### Issue 1: No Good Cleavage Site Found

**Cause:** Backbone doesn't cross near ligand; all crossings are in helices/sheets.

**Solutions:**
1. Generate more monomer designs (different folds)
2. Use longer chains (more likely to have loops near ligand)
3. Try `is_non_loopy: false` to encourage more loops

### Issue 2: Very Asymmetric Cleavage

**Cause:** Best cleavage site near N or C terminus.

**Solutions:**
1. Add `length_balance` weight to scoring
2. Accept asymmetry (hetero-dimer is fine)
3. For homo-dimer: only use designs with central cleavage

### Issue 3: Chains Drift Apart After Cleavage

**Cause:** Without covalent linkage, chains may not stay associated.

**Solutions:**
1. This is actually DESIRED for dimerization binder applications
2. For stable complex: run brief MD equilibration
3. Validate with AF3 that chains + ligand form stable complex

### Issue 4: New Termini Clash with Ligand

**Cause:** Cleavage creates N/C termini pointing at ligand.

**Solutions:**
1. Remove larger loop (more residues)
2. Add terminus caps (ACE/NME)
3. Brief energy minimization

---

## Validation Protocol

```python
def validate_cleaved_dimer(dimer, ligand_name="UNL"):
    """Full validation of cleaved dimer."""
    
    results = {}
    
    # 1. Structural integrity
    chain_a = dimer[dimer.chain_id == 'A']
    chain_b = dimer[dimer.chain_id == 'B']
    ligand = dimer[dimer.res_name == ligand_name]
    
    results['chain_a_length'] = len(set(chain_a.res_id))
    results['chain_b_length'] = len(set(chain_b.res_id))
    results['length_pass'] = min(results['chain_a_length'], 
                                  results['chain_b_length']) >= 20
    
    # 2. No clashes
    min_dist = compute_min_distance(chain_a, chain_b)
    results['min_interchain_distance'] = min_dist
    results['clash_pass'] = min_dist > 2.5
    
    # 3. Both chains contact ligand
    contacts_a = count_contacts(chain_a, ligand, cutoff=4.0)
    contacts_b = count_contacts(chain_b, ligand, cutoff=4.0)
    results['contacts_a'] = contacts_a
    results['contacts_b'] = contacts_b
    results['contact_pass'] = contacts_a >= 3 and contacts_b >= 3
    
    # 4. Contact symmetry
    symmetry = min(contacts_a, contacts_b) / max(contacts_a, contacts_b)
    results['contact_symmetry'] = symmetry
    results['symmetry_pass'] = symmetry > 0.4
    
    # 5. GNINA affinity
    affinity = run_gnina(dimer, ligand_name)
    results['gnina_affinity'] = affinity
    results['affinity_pass'] = affinity < -5.0
    
    # 6. Overall
    results['overall_pass'] = all([
        results['length_pass'],
        results['clash_pass'],
        results['contact_pass'],
        results['affinity_pass']
    ])
    
    return results
```

---

## Example Usage

```python
# Step 1: Design monomer with buried ligand
config = {
    "azo_central": {
        "input": "azobenzene.pdb",
        "contig": "80-100",
        "ligand": "UNL",
        "select_buried": {"UNL": "ALL"},  # Force to center
        "ori_token": [0.0, 0.0, 0.0]
    }
}

monomer = run_rfd3(config)

# Step 2: Find cleavage sites
candidates = find_cleavage_sites(monomer, ligand_atoms)
best_site = max(candidates, key=lambda c: score_cleavage_site(c, n_residues))

# Step 3: Cleave
dimer, info = cleave_protein(monomer, best_site['residue_id'], ligand_atoms)

# Step 4: Validate
results = validate_cleaved_dimer(dimer)
print(f"GNINA affinity: {results['gnina_affinity']:.2f} kcal/mol")
print(f"Contacts: A={results['contacts_a']}, B={results['contacts_b']}")

# Step 5 (optional): Create homo-dimer
if user_wants_homo_dimer:
    homo_dimer, homo_info = create_homo_dimer(dimer, ligand_atoms)
```

---

## Summary

| Step | Action | Key Consideration |
|------|--------|-------------------|
| 1 | Design monomer | Use `select_buried` to center ligand |
| 2 | Find cleavage sites | Must be in loop, both chains contact ligand |
| 3 | Score sites | Balance contacts, symmetry, and length |
| 4 | Cleave | Remove full loop for clean interface |
| 5 | Validate | GNINA affinity must be < -5 kcal/mol |
| 6 | Homo-dimer (optional) | Requires C2-symmetric ligand |
