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

## Current State Summary

### Completed Features
- [x] **Unified Backend** (`backend/main.py`) - Auto-detects Foundry, falls back to mock
- [x] **Enhanced Health Endpoint** - Shows mode, GPU name/memory, model status
- [x] **Notification Toast System** - Guides users through workflow
- [x] **Panel Data Flow** - "Use Latest Design" button in MPNN panel
- [x] **Molstar 3D Viewer** - Protein structure visualization
- [x] **Field Name Compatibility** - Accepts both `contig` and `contigs`

### Pending Features
- [ ] **ContigBuilder.tsx** - Visual UI for building contig specifications
- [ ] **WorkflowStepper.tsx** - Horizontal progress indicator in header
- [ ] **Real Foundry Integration** - Requires Python 3.12+ pod for `rc-foundry[all]`

---

## For Real Model Inference

The `rc-foundry` package requires Python 3.12+. Current RunPod PyTorch templates use Python 3.10.

**Options:**
1. Use a custom Docker image with Python 3.12
2. Wait for RosettaCommons to release a Python 3.10 compatible version
3. Install Python 3.12 manually on the pod

Once available, install with:
```bash
pip install "rc-foundry[all]"
export FOUNDRY_CHECKPOINT_DIRS=/workspace/checkpoints
```

---

## Key Files

| File | Purpose |
|------|---------|
| `backend/main.py` | Unified FastAPI backend (mock + real mode) |
| `frontend/src/components/ConnectionStatus.tsx` | Backend status display |
| `frontend/src/components/RFD3Panel.tsx` | RFdiffusion design panel |
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
- Frontend has notifications, 3D viewer, panel data flow (working)
- rc-foundry requires Python 3.12+ (blocked on RunPod template)

Next priorities:
1. Create ContigBuilder.tsx - Visual UI to help beginners build contig strings
2. Create WorkflowStepper.tsx - Progress indicator: Design → Validate → Sequences
3. (Optional) Find Python 3.12 RunPod solution for real Foundry models

The plan file at ~/.claude/plans/velvet-snuggling-pony.md has full details.
```
