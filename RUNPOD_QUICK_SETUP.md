# RunPod Quick Setup Guide

## Overview
This guide enables fast GPU backend setup on RunPod for the Foundry Protein Design frontend.

## Prerequisites
- RunPod account with credits
- Network volume `foundry-checkpoints` (optional but recommended - persists data between pod restarts)

---

## One-Cell Setup (Recommended)

### Step 1: Deploy Pod
1. Go to [RunPod Console](https://console.runpod.io/deploy)
2. Select **A40** GPU ($0.40/hr) or **RTX 4090** or any GPU with 24GB+ VRAM
3. Click **Edit** on the template section
4. Set **HTTP Ports**: `8888, 8000` (IMPORTANT: add port 8000!)
5. Optionally attach network volume `foundry-checkpoints`
6. Click **Set Overrides** → **Deploy On-Demand**

### Step 2: Run Single Setup Cell
1. Wait for pod to show "Ready" status
2. Click **Jupyter Lab** link
3. Create new Python notebook and run this **ONE CELL**:

```python
# ============================================================
# FOUNDRY PROTEIN DESIGN - COMPLETE SETUP (Single Cell)
# ============================================================
# This cell handles everything: dependencies, models, and server
# Run this after ANY pod restart to get a working backend
# Takes ~10-15 minutes on first run (downloading checkpoints)
# Takes ~30 seconds on subsequent runs (if using network volume)
# ============================================================

import subprocess
import os
import sys
import time

print("="*60)
print("FOUNDRY PROTEIN DESIGN - AUTOMATED SETUP")
print("="*60)

# Configuration
WORKSPACE = "/workspace"
VENV_DIR = f"{WORKSPACE}/foundry_env"
CHECKPOINT_DIR = f"{WORKSPACE}/checkpoints"
MAIN_PY_URL = "https://raw.githubusercontent.com/polaires/Banta_Lab_RFdiffusion/main/backend/main.py"

def run(cmd, desc="", check=True):
    """Run a command with status output"""
    print(f"\n[STEP] {desc or cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0 and check:
        print(f"[WARN] {result.stderr[:500] if result.stderr else 'Command had non-zero exit'}")
    return result.returncode == 0

def is_server_running():
    """Check if server is already running"""
    try:
        import urllib.request
        urllib.request.urlopen("http://localhost:8000/health", timeout=2)
        return True
    except:
        return False

# Check if server is already running
if is_server_running():
    print("\n[OK] Server is already running!")
    print("[OK] Visit: http://localhost:8000/health")
    print("[OK] Or use your RunPod proxy URL")
else:
    # Step 1: Install system dependencies
    print("\n" + "="*60)
    print("STEP 1/6: Installing system dependencies...")
    print("="*60)
    run("apt-get update -qq && apt-get install -y -qq python3.12 python3.12-venv python3.12-dev curl > /dev/null 2>&1",
        "Installing Python 3.12 and curl")

    # Step 2: Create virtual environment (if not exists)
    print("\n" + "="*60)
    print("STEP 2/6: Setting up Python virtual environment...")
    print("="*60)
    if not os.path.exists(f"{VENV_DIR}/bin/python"):
        run(f"python3.12 -m venv {VENV_DIR}", "Creating virtual environment")
        run(f"{VENV_DIR}/bin/pip install --upgrade pip -q", "Upgrading pip")
    else:
        print("[OK] Virtual environment already exists")

    # Step 3: Install rc-foundry (if not installed)
    print("\n" + "="*60)
    print("STEP 3/6: Installing rc-foundry...")
    print("="*60)
    foundry_installed = os.path.exists(f"{VENV_DIR}/bin/foundry")
    if not foundry_installed:
        run(f'{VENV_DIR}/bin/pip install "rc-foundry[all]" -q', "Installing rc-foundry (this may take a few minutes)")
    else:
        print("[OK] rc-foundry already installed")

    # Step 4: Install FastAPI dependencies
    print("\n" + "="*60)
    print("STEP 4/6: Installing FastAPI dependencies...")
    print("="*60)
    run(f"{VENV_DIR}/bin/pip install fastapi uvicorn python-multipart biotite -q", "Installing FastAPI, uvicorn, biotite")

    # Step 5: Download model checkpoints (if not present)
    print("\n" + "="*60)
    print("STEP 5/6: Downloading model checkpoints...")
    print("="*60)
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)

    # Check if checkpoints exist
    rfd3_exists = os.path.exists(f"{CHECKPOINT_DIR}/rfd3") and len(os.listdir(f"{CHECKPOINT_DIR}/rfd3")) > 0

    if not rfd3_exists:
        print("[INFO] Downloading checkpoints (~10GB). This may take 5-10 minutes...")
        run(f"{VENV_DIR}/bin/foundry install base-models --checkpoint-dir {CHECKPOINT_DIR}",
            "Downloading RFD3, RF3, and MPNN checkpoints")
    else:
        print("[OK] Checkpoints already downloaded")
        # List installed checkpoints
        run(f"{VENV_DIR}/bin/foundry list-installed --checkpoint-dir {CHECKPOINT_DIR}",
            "Listing installed checkpoints", check=False)

    # Step 6: Download and start server
    print("\n" + "="*60)
    print("STEP 6/6: Starting FastAPI server...")
    print("="*60)

    # Download latest main.py
    run(f"curl -sL {MAIN_PY_URL} -o {WORKSPACE}/main.py", "Downloading latest main.py from GitHub")

    # Start server with proper environment
    env = os.environ.copy()
    env["FOUNDRY_CHECKPOINT_DIRS"] = CHECKPOINT_DIR

    # Kill any existing server on port 8000
    subprocess.run("fuser -k 8000/tcp 2>/dev/null", shell=True, capture_output=True)
    time.sleep(1)

    # Start server
    subprocess.Popen(
        [f"{VENV_DIR}/bin/python", "-m", "uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"],
        env=env,
        cwd=WORKSPACE,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )

    # Wait for server to start
    print("[INFO] Waiting for server to start...")
    for i in range(30):
        time.sleep(1)
        if is_server_running():
            break

    # Verify server is running
    print("\n" + "="*60)
    print("SETUP COMPLETE!")
    print("="*60)

    if is_server_running():
        import urllib.request
        import json
        try:
            response = urllib.request.urlopen("http://localhost:8000/health", timeout=5)
            health = json.loads(response.read().decode())
            print(f"\n[SUCCESS] Server is running!")
            print(f"  Mode: {health.get('mode', 'unknown')}")
            print(f"  GPU: {health.get('gpu_name', 'N/A')} ({health.get('gpu_memory_gb', 0):.0f}GB)")
            print(f"  Status: {health.get('status', 'unknown')}")
        except Exception as e:
            print(f"[OK] Server started (health check: {e})")
    else:
        print("[WARN] Server may still be starting. Check in a few seconds.")

    print("\n" + "-"*60)
    print("NEXT STEPS:")
    print("-"*60)
    print("1. Copy your API URL from RunPod console:")
    print("   https://<POD_ID>-8000.proxy.runpod.net")
    print("")
    print("2. Enter this URL in the Frontend:")
    print("   - Open the Foundry Protein Design web app")
    print("   - Paste the URL in 'Backend Connection' field")
    print("   - Click 'Connect'")
    print("")
    print("3. Test locally (optional):")
    print('   !curl -s http://localhost:8000/health | python -m json.tool')
    print("="*60)
```

### Step 3: Connect Frontend
1. Copy your pod's API URL from RunPod console: `https://{POD_ID}-8000.proxy.runpod.net`
2. Open the Foundry Protein Design web app
3. Paste the URL in the **Backend Connection** field at the top
4. Click **Connect**

That's it! The frontend will automatically connect and you can start designing proteins.

---

## Verify Setup (Optional)

Run this cell to check server status:

```python
!curl -s http://localhost:8000/health | python -m json.tool
```

Expected output:
```json
{
    "status": "healthy",
    "mode": "real",
    "gpu_available": true,
    "gpu_name": "NVIDIA A40",
    "gpu_memory_gb": 45.0,
    "models": {
        "rfd3": {"available": true, "checkpoint_exists": true, "checkpoint_size_gb": 3.5},
        "rf3": {"available": true, "checkpoint_exists": true, "checkpoint_size_gb": 2.8},
        "proteinmpnn": {"available": true, "checkpoint_exists": true, "checkpoint_size_gb": 0.5}
    }
}
```

Test RFD3 design:
```python
!curl -s -X POST http://localhost:8000/api/rfd3/design \
  -H "Content-Type: application/json" \
  -d '{"contig": "100"}' | python -m json.tool
```

---

## After Pod Restart

Just run the **same single cell** again! It will:
- Skip already-installed packages
- Skip already-downloaded checkpoints
- Start the server immediately

If using a network volume, restart takes ~30 seconds.

---

## Troubleshooting

### Port 8000 not accessible externally
**Cause:** Port 8000 wasn't added to HTTP services when creating the pod.
**Fix:** Terminate pod and redeploy with `8888, 8000` in HTTP Ports field.

### "Connection refused" or CORS errors
**Cause:** Server not running or wrong URL.
**Fix:**
1. Re-run the setup cell
2. Verify URL matches pod ID exactly (check RunPod console)

### Server starts but shows "mock mode"
**Cause:** Checkpoints not downloaded or wrong path.
**Fix:**
```python
# Check if checkpoints exist
!ls -la /workspace/checkpoints/
# Re-download if needed
!/workspace/foundry_env/bin/foundry install base-models --checkpoint-dir /workspace/checkpoints
```

### Out of GPU memory
**Cause:** A40 has 45GB, some operations need more.
**Fix:** Try RTX 4090 or A100 if available.

---

## Cost Management

| Resource | Cost |
|----------|------|
| A40 GPU running | $0.40/hr |
| Network volume (50GB) | $0.005/hr |
| Pod stopped | $0 (only volume cost) |

**Tip:** Stop pod when not in use. Network volume persists all files including checkpoints.

---

## Pod Configuration Summary

| Setting | Value |
|---------|-------|
| GPU | A40 (48GB VRAM) recommended |
| Template | Runpod Pytorch 2.4.0 |
| HTTP Ports | `8888, 8000` |
| Network Volume | `foundry-checkpoints` (optional) |
| Mount Path | `/workspace` |

---

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/health` | GET | Health check with GPU, mode, and model status |
| `/api/rfd3/design` | POST | RFdiffusion3 structure design |
| `/api/rf3/predict` | POST | RosettaFold3 structure prediction |
| `/api/mpnn/design` | POST | ProteinMPNN sequence design |
| `/api/jobs/{job_id}` | GET | Get job status and results |
| `/api/jobs` | GET | List all jobs |
| `/api/jobs/{job_id}` | DELETE | Delete a job |
