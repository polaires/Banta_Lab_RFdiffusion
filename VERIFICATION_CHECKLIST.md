# Binder Pipeline Verification Checklist

## Pre-requisites

- [ ] Docker Desktop running
- [ ] Node.js and npm/pnpm installed
- [ ] GPU available (for full ESMFold testing)

---

## Part 1: Docker Container Build & Test

### 1.1 Build the Docker Container

```bash
cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless
docker build -t banta-rfdiffusion:latest .
```

**Expected**: Build completes without errors

### 1.2 Verify New Python Files Exist in Container

```bash
docker run --rm banta-rfdiffusion:latest ls -la /app/
```

**Expected files**:
- [ ] `handler.py` (updated with protocol presets)
- [ ] `binding_analysis.py` (updated with SC, hydrophobicity, unsaturated H-bonds)
- [ ] `esmfold_utils.py` (NEW - ESMFold validation)
- [ ] `hotspot_detection.py`

### 1.3 Run Container with GPU

```bash
docker run --gpus all -p 8080:8080 banta-rfdiffusion:latest
```

**Expected**: Container starts, logs show RFDiffusion model loading

---

## Part 2: Backend API Tests

### 2.1 Test Basic Health Check

```bash
curl http://localhost:8080/health
```

**Expected**: `{"status": "ok"}`

### 2.2 Test Hotspot Detection (Minimal)

```bash
curl -X POST http://localhost:8080/detect_hotspots \
  -H "Content-Type: application/json" \
  -d '{"pdb_id": "1BRF", "target_chain": "A"}'
```

**Expected response should include**:
- [ ] `status: "completed"`
- [ ] `hotspots` array with residues
- [ ] `hotspot_centers` clusters

### 2.3 Test Binder Design with Protocol Preset

```bash
curl -X POST http://localhost:8080/protein_binder_design \
  -H "Content-Type: application/json" \
  -d '{
    "pdb_id": "1BRF",
    "target_chain": "A",
    "protocol": "miniprotein_default",
    "num_designs": 2,
    "validate_structure": false
  }'
```

**Expected response should include**:
- [ ] `status: "completed"`
- [ ] `designs` array with 0-2 designs
- [ ] Each design has `binder_sequence`, `pdb_content`
- [ ] `statistics` object with counts

### 2.4 Test New Metrics in Binder Design

Run a single design and verify new metrics:

```bash
curl -X POST http://localhost:8080/protein_binder_design \
  -H "Content-Type: application/json" \
  -d '{
    "pdb_id": "1BRF",
    "target_chain": "A",
    "protocol": "miniprotein_default",
    "num_designs": 4,
    "quality_threshold": "relaxed"
  }'
```

**Verify each passing design includes**:
- [ ] `shape_complementarity` (float, typically 0.3-0.8)
- [ ] `surface_hydrophobicity` (float, < 0.37 for standard)
- [ ] `unsaturated_hbonds` (int, < 6 for standard)
- [ ] `interface_residue_count` (int, > 6 for standard)

### 2.5 Test ESMFold Validation (GPU Required)

```bash
curl -X POST http://localhost:8080/protein_binder_design \
  -H "Content-Type: application/json" \
  -d '{
    "pdb_id": "1BRF",
    "target_chain": "A",
    "protocol": "miniprotein_default",
    "num_designs": 2,
    "validate_structure": true
  }'
```

**Verify designs include**:
- [ ] `esmfold_plddt` (float, 0-1)
- [ ] `esmfold_rmsd` (float, Angstroms)
- [ ] `esmfold_validation_passed` (bool)

### 2.6 Test Protocol Presets

Test each protocol preset returns appropriate binder lengths:

| Protocol | Expected Length Range |
|----------|----------------------|
| `miniprotein_default` | 60-100 |
| `miniprotein_hardtarget` | 80-120 |
| `peptide_default` | 15-30 |
| `peptide_helical` | 15-25 |
| `large_binder` | 100-150 |

### 2.7 Test Quality Thresholds

Run with different thresholds and verify filter behavior:

```bash
# Relaxed - more designs pass
curl -X POST ... -d '{"quality_threshold": "relaxed", ...}'

# Standard - balanced filtering
curl -X POST ... -d '{"quality_threshold": "standard", ...}'

# Strict - fewer designs pass
curl -X POST ... -d '{"quality_threshold": "strict", ...}'
```

**Expected**: Pass rate decreases as strictness increases

---

## Part 3: Frontend Tests

### 3.1 Build Frontend

```bash
cd G:\Github_local_repo\Banta_Lab_RFdiffusion\frontend
pnpm install
pnpm build
```

**Expected**: Build completes without TypeScript errors

### 3.2 Run Frontend Dev Server

```bash
pnpm dev
```

**Expected**: Server starts on http://localhost:3000

### 3.3 Test Binder Interview Flow

1. [ ] Navigate to binder design page
2. [ ] Click "Design a Binder" or start interview
3. [ ] **Step 1**: Select target (Rubredoxin)
4. [ ] **Step 2**: Select protocol type (NEW STEP!)
   - [ ] Verify 5 protocol options displayed
   - [ ] Select "Miniprotein (Recommended)"
5. [ ] **Step 3**: Select binder size
6. [ ] **Step 4**: Select hotspot method
7. [ ] **Step 5**: Select number of designs
8. [ ] **Step 6**: Select quality threshold
9. [ ] **Step 7**: Select priority
10. [ ] Click "Review Design"

### 3.4 Verify Protocol Selection Works

After completing interview, check that preferences include:
- [ ] `protocol: "miniprotein_default"` (or selected value)
- [ ] `protocolLabel: "Miniprotein (Standard)"`
- [ ] `validateStructure: true`

### 3.5 Test API Integration

If backend is running:
1. [ ] Complete interview and submit
2. [ ] Verify job starts successfully
3. [ ] Check job result includes new metrics
4. [ ] Verify design cards show shape complementarity

---

## Part 4: End-to-End Workflow Test

### 4.1 Full Workflow Test

1. [ ] Start Docker container with GPU
2. [ ] Start frontend dev server
3. [ ] Complete binder interview selecting:
   - Target: Rubredoxin (1BRF)
   - Protocol: Miniprotein (Recommended)
   - Length: Medium (60-80)
   - Hotspots: Auto-detect
   - Designs: 4
   - Quality: Standard
   - Priority: Balanced

4. [ ] Submit design job
5. [ ] Wait for completion (5-10 min)
6. [ ] Verify results:
   - [ ] At least 1 design passes
   - [ ] Design shows shape complementarity score
   - [ ] Design shows ESMFold pLDDT (if enabled)
   - [ ] 3D viewer displays binder structure

### 4.2 Error Handling Test

1. [ ] Submit with invalid PDB ID - should return error
2. [ ] Submit with no hotspots found - should handle gracefully
3. [ ] Cancel job mid-execution - should stop cleanly

---

## Part 5: Regression Tests

### 5.1 Hotspot Detection Still Works

```bash
curl -X POST http://localhost:8080/detect_hotspots \
  -d '{"pdb_id": "1BRF", "target_chain": "A"}'
```

- [ ] Returns hotspots with SASA scores
- [ ] Returns clustered hotspot centers

### 5.2 Basic Binder Design Still Works

Without new features:

```bash
curl -X POST http://localhost:8080/protein_binder_design \
  -d '{
    "pdb_id": "1BRF",
    "target_chain": "A",
    "binder_length": "60-80",
    "num_designs": 2
  }'
```

- [ ] Returns designs without protocol preset
- [ ] Backward compatible with old API calls

---

## Troubleshooting

### Docker GPU Issues

```bash
# Check NVIDIA Docker runtime
docker run --gpus all nvidia/cuda:11.8-base nvidia-smi
```

### ESMFold Memory Issues

If OOM errors occur:
- Reduce sequence length
- Set `validate_structure: false` to skip ESMFold
- Use smaller binders (peptide protocols)

### Frontend TypeScript Errors

```bash
# Check for type errors
pnpm tsc --noEmit
```

### API Connection Issues

```bash
# Check container logs
docker logs <container_id>

# Check if port is accessible
curl http://localhost:8080/health
```

---

## Sign-off

| Test Section | Passed | Notes |
|--------------|--------|-------|
| Docker Build | [ ] | |
| Backend API | [ ] | |
| Frontend Build | [ ] | |
| Interview Flow | [ ] | |
| End-to-End | [ ] | |
| Regression | [ ] | |

**Tested by**: _______________
**Date**: _______________
**Version**: Post-BindCraft improvements (Phase 1-5)
