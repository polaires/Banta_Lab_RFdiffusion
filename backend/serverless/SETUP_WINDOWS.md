# Complete Windows Setup Guide for Local RFdiffusion Testing

Your system:
- GPU: NVIDIA RTX 3090 (24GB) ✅
- Driver: 581.63 ✅
- WSL2: Ubuntu installed ✅
- Docker: Not installed ❌

## Step 1: Install Docker Desktop

1. **Download Docker Desktop:**
   - Go to: https://www.docker.com/products/docker-desktop/
   - Click "Download for Windows"

2. **Run the installer:**
   - Double-click `Docker Desktop Installer.exe`
   - ✅ Check "Use WSL 2 instead of Hyper-V"
   - ✅ Check "Add shortcut to desktop"
   - Click "Ok" and wait for installation

3. **Restart your computer** when prompted

4. **Start Docker Desktop:**
   - Launch from Start Menu or desktop shortcut
   - Accept the license agreement
   - Skip the tutorial/sign-in (optional)

5. **Verify installation** (in PowerShell):
   ```powershell
   docker --version
   docker run hello-world
   ```

## Step 2: Configure Docker for WSL2 + GPU

1. **Open Docker Desktop Settings** (gear icon)

2. **General:**
   - ✅ "Use the WSL 2 based engine" (should be checked)

3. **Resources → WSL Integration:**
   - ✅ Enable integration with your Ubuntu distro
   - Click "Apply & Restart"

4. **Verify GPU access:**
   ```powershell
   docker run --rm --gpus all nvidia/cuda:12.4.1-base-ubuntu22.04 nvidia-smi
   ```

   You should see your RTX 3090 listed.

## Step 3: Download Model Checkpoints

The models are ~15GB total. Choose one method:

### Option A: Download via RunPod (Recommended if you have checkpoints there)

If your RunPod Network Volume already has checkpoints, we can copy them.

### Option B: Download locally via Foundry CLI

```powershell
# In PowerShell:
cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless
mkdir checkpoints

# Install Foundry (requires Python 3.12)
pip install rc-foundry[all]

# Set checkpoint directory and download
$env:FOUNDRY_CHECKPOINT_DIRS = "G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless\checkpoints"
foundry install base-models
```

### Option C: Let container download on first run

Start without checkpoints - it will run in "mock mode" for testing the API:
```powershell
cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless
mkdir checkpoints
docker-compose -f docker-compose.local.yml up --build
```

Then trigger download from inside container:
```powershell
curl http://localhost:8000/runsync -X POST -H "Content-Type: application/json" -d "{\"input\": {\"task\": \"download_checkpoints\"}}"
```

## Step 4: Build and Run Container

```powershell
cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless

# Build and start (first time takes ~10-15 minutes)
docker-compose -f docker-compose.local.yml up --build
```

You should see:
```
[Handler] Starting RunPod serverless handler...
[Handler] Initializing worker...
[Handler] GPU: {'available': True, 'name': 'NVIDIA GeForce RTX 3090', ...}
[Handler] Foundry available: True
...
[Handler] Worker initialization complete
```

## Step 5: Test

In another PowerShell window:

```powershell
cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless

# Health check
python test_local.py health

# Test RFD3 (will take ~30-60s first time as models load)
python test_local.py rfd3
```

## Troubleshooting

### "error during connect: This error may indicate that the docker daemon is not running"
→ Start Docker Desktop from the Start Menu

### "docker: Error response from daemon: could not select device driver"
→ GPU not accessible. Run:
```powershell
docker run --rm --gpus all nvidia/cuda:12.4.1-base-ubuntu22.04 nvidia-smi
```
If this fails, reinstall NVIDIA drivers or Docker Desktop.

### "CUDA out of memory"
→ Close other GPU applications (games, other ML models)

### Container exits immediately
→ Check logs: `docker-compose -f docker-compose.local.yml logs`

### Slow build
→ First build downloads ~10GB of base image + dependencies. Subsequent builds are cached.
