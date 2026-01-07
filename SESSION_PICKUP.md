# Session Pickup Guide - Foundry Protein Design

## Quick Pod Restart (Copy-Paste Ready)

### Step 1: Deploy New RunPod
1. Go to [RunPod Console](https://console.runpod.io/deploy)
2. Select **A40** GPU ($0.40/hr) or any GPU with 24GB+ VRAM
3. **IMPORTANT**: Set HTTP Ports to `8888, 8000`
4. Attach network volume `foundry-checkpoints` if available
5. Deploy On-Demand

### Step 2: Start Backend (Run in Jupyter)
Open JupyterLab and run these cells:

```python
# Cell 1: Install dependencies
!pip install fastapi uvicorn python-multipart -q

# Cell 2: Download latest backend from GitHub
!curl -o main.py https://raw.githubusercontent.com/polaires/Banta_Lab_RFdiffusion/main/backend/main.py

# Cell 3: Start server
import subprocess
subprocess.Popen(["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"])

# Cell 4: Verify
!sleep 2 && curl -s http://localhost:8000/health | python -m json.tool
```

### Step 3: Connect Frontend
Update `frontend/.env.local`:
```
NEXT_PUBLIC_API_URL=https://{YOUR_POD_ID}-8000.proxy.runpod.net
```
Then restart frontend: `npm run dev`

---

## Current State Summary (Updated Jan 7, 2026)

### Completed Features
- [x] **Unified Backend** (`backend/main.py`) - Auto-detects Foundry, falls back to mock
- [x] **Enhanced Health Endpoint** - Shows mode, GPU name/memory, model status
- [x] **Notification Toast System** - Guides users through workflow
- [x] **Panel Data Flow** - "Use Latest Design" button in MPNN panel
- [x] **Molstar 3D Viewer** - Protein structure visualization
- [x] **Field Name Compatibility** - Accepts both `contig` and `contigs`
- [x] **ContigBuilder.tsx** - Visual UI for building contig specifications
- [x] **WorkflowStepper.tsx** - Horizontal progress indicator showing RFD3 → RF3 → MPNN

### In Progress / Next Steps
- [ ] **Complete RFD3 Checkpoint Download** - RF3 downloaded successfully, RFD3 still needed
- [ ] **Test Real Inference End-to-End** - Once all checkpoints are installed

---

## For Real Model Inference

### Official Foundry Documentation
https://github.com/RosettaCommons/foundry/tree/production

### Setup Steps (Verified Working)

**Step 1: Create Python 3.12 venv**
```python
!apt update && apt install -y python3.12 python3.12-venv python3.12-dev -qq
!python3.12 -m venv /workspace/foundry_env
!/workspace/foundry_env/bin/pip install --upgrade pip -q
```

**Step 2: Install rc-foundry**
```python
!/workspace/foundry_env/bin/pip install "rc-foundry[all]" -q
```

**Step 3: Download Checkpoints (CRITICAL)**
```python
# Download all base models at once
!/workspace/foundry_env/bin/foundry install base-models --checkpoint-dir /workspace/checkpoints

# OR download individually
!/workspace/foundry_env/bin/foundry install rfd3 --checkpoint-dir /workspace/checkpoints
!/workspace/foundry_env/bin/foundry install rf3 --checkpoint-dir /workspace/checkpoints
!/workspace/foundry_env/bin/foundry install proteinmpnn --checkpoint-dir /workspace/checkpoints
```

**Verify installation:**
```bash
/workspace/foundry_env/bin/foundry list-installed
```

**Step 4: Start Backend**
```python
import subprocess, os
env = os.environ.copy()
env["FOUNDRY_CHECKPOINT_DIRS"] = "/workspace/checkpoints"
subprocess.Popen(
    ["/workspace/foundry_env/bin/python", "-m", "uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"],
    env=env
)
```

---

## Key Files

| File | Purpose |
|------|---------|
| `backend/main.py` | Unified FastAPI backend (mock + real mode) |
| `frontend/src/components/ConnectionStatus.tsx` | Backend status display |
| `frontend/src/components/RFD3Panel.tsx` | RFdiffusion design panel with ContigBuilder |
| `frontend/src/components/ContigBuilder.tsx` | Visual segment builder for contig strings |
| `frontend/src/components/WorkflowStepper.tsx` | Horizontal workflow progress indicator |
| `frontend/src/components/MPNNPanel.tsx` | ProteinMPNN panel with "Use Latest Design" |
| `frontend/src/components/ProteinViewer.tsx` | Molstar 3D viewer |
| `frontend/src/components/NotificationToast.tsx` | Workflow notifications |
| `frontend/src/lib/store.ts` | Zustand state management |

---

## Prompt for Next Session

```
Continue working on Banta_Lab_RFdiffusion:

Current state:
- Backend unified with mock/real mode auto-detection (working)
- Frontend complete with ContigBuilder, WorkflowStepper, notifications, 3D viewer
- rc-foundry installed in Python 3.12 venv
- RF3 checkpoint downloaded successfully
- RFD3 checkpoint still needs to be downloaded

Next priorities:
1. Run: /workspace/foundry_env/bin/foundry install rfd3 --checkpoint-dir /workspace/checkpoints
2. Test real RFD3 inference from frontend
3. If working, test full workflow: RFD3 → MPNN

See RUNPOD_QUICK_SETUP.md for detailed setup instructions.
Official Foundry docs: https://github.com/RosettaCommons/foundry/tree/production
```
