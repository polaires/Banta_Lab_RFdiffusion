#!/usr/bin/env python3
"""
RunPod FastAPI Startup Script
Run this single file in Jupyter to start the backend:

    %run /workspace/startup.py

Or from terminal:
    python /workspace/startup.py
"""

import subprocess
import sys
import time

def main():
    print("=" * 50)
    print("Foundry Backend Startup")
    print("=" * 50)

    # Step 1: Install dependencies
    print("\n[1/3] Installing dependencies...")
    subprocess.run([sys.executable, "-m", "pip", "install", "-q",
                   "fastapi", "uvicorn", "python-multipart"])
    print("      Done!")

    # Step 2: Start server
    print("\n[2/3] Starting FastAPI server on port 8000...")
    proc = subprocess.Popen(
        ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"],
        cwd="/workspace"
    )
    time.sleep(2)  # Wait for server to start

    # Step 3: Verify
    print("\n[3/3] Verifying server...")
    result = subprocess.run(
        ["curl", "-s", "http://localhost:8000/health"],
        capture_output=True, text=True
    )

    if "healthy" in result.stdout:
        print("      Server is running!")
        print(f"      Response: {result.stdout}")
        print("\n" + "=" * 50)
        print("SUCCESS! Backend is ready.")
        print("=" * 50)
        print("\nYour API URL:")

        # Try to get pod ID from environment or hostname
        import os
        pod_id = os.environ.get("RUNPOD_POD_ID", "<POD_ID>")
        print(f"  https://{pod_id}-8000.proxy.runpod.net")
        print("\nUpdate frontend/.env.local with this URL")
    else:
        print("      WARNING: Server may not be running correctly")
        print(f"      Response: {result.stdout}")
        print(f"      Error: {result.stderr}")

    return proc

if __name__ == "__main__":
    main()
