"""Utility handlers for checkpoint management and file operations.

Tasks: download_checkpoints, delete_file
"""

import os
import subprocess
import traceback
from typing import Dict, Any

# Read checkpoint dir from env (set by handler.py at startup)
CHECKPOINT_DIR = os.environ.get(
    "FOUNDRY_CHECKPOINT_DIRS",
    os.environ.get("CHECKPOINT_DIR", "/runpod-volume/checkpoints"),
)


def handle_download_checkpoints(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Download model checkpoints to the Network Volume using Foundry CLI.

    This task downloads the required model checkpoints (RFD3, RF3, MPNN) to the
    Network Volume at /runpod-volume/checkpoints/. This only needs to be run once.

    Input:
        force: bool - Force re-download even if checkpoints exist (default: False)

    Returns:
        status: "completed" | "failed"
        result: {
            checkpoints_dir: str,
            before: {files, total_size_gb},
            after: {files, total_size_gb},
            downloaded: bool,
            message: str
        }
    """
    import glob as glob_module

    force = job_input.get("force", False)

    def get_checkpoint_info(directory: str) -> Dict[str, Any]:
        """Get info about checkpoint files in directory"""
        files = []
        total_size = 0

        if os.path.exists(directory):
            # Look for checkpoint files
            patterns = ["**/*.ckpt", "**/*.pt", "**/*.pth", "**/*.safetensors"]
            for pattern in patterns:
                for filepath in glob_module.glob(os.path.join(directory, pattern), recursive=True):
                    size = os.path.getsize(filepath)
                    files.append({
                        "name": os.path.basename(filepath),
                        "path": filepath,
                        "size_mb": round(size / (1024 * 1024), 2)
                    })
                    total_size += size

        return {
            "files": files,
            "file_count": len(files),
            "total_size_gb": round(total_size / (1024 * 1024 * 1024), 2)
        }

    try:
        print(f"[Handler] Checking checkpoints at: {CHECKPOINT_DIR}")

        # Ensure checkpoint directory exists
        os.makedirs(CHECKPOINT_DIR, exist_ok=True)

        # Get current state
        before_info = get_checkpoint_info(CHECKPOINT_DIR)
        print(f"[Handler] Before: {before_info['file_count']} files, {before_info['total_size_gb']} GB")

        # Check if we already have checkpoints
        has_checkpoints = before_info["total_size_gb"] > 1.0  # At least 1GB of checkpoints

        if has_checkpoints and not force:
            return {
                "status": "completed",
                "result": {
                    "checkpoints_dir": CHECKPOINT_DIR,
                    "before": before_info,
                    "after": before_info,
                    "downloaded": False,
                    "message": f"Checkpoints already exist ({before_info['total_size_gb']} GB). Use force=True to re-download."
                }
            }

        # If force is enabled, delete any 0-byte (corrupted) checkpoint files first
        if force:
            for file_info in before_info.get("files", []):
                if file_info.get("size_mb", 0) == 0:
                    filepath = file_info.get("path")
                    if filepath and os.path.exists(filepath):
                        print(f"[Handler] Deleting corrupted 0-byte file: {filepath}")
                        os.remove(filepath)

        # Download checkpoints using Foundry CLI
        cmd = ["foundry", "install", "base-models"]
        if force:
            cmd.append("--force")
        print(f"[Handler] Downloading checkpoints using Foundry CLI...")
        print(f"[Handler] Running: {' '.join(cmd)}")

        # Set environment for Foundry
        env = os.environ.copy()
        env["FOUNDRY_CHECKPOINT_DIRS"] = CHECKPOINT_DIR

        # Run foundry install
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            env=env,
            timeout=3600  # 1 hour timeout for downloads
        )

        print(f"[Handler] Foundry stdout: {result.stdout}")
        if result.stderr:
            print(f"[Handler] Foundry stderr: {result.stderr}")

        if result.returncode != 0:
            # Try alternative: foundry download
            print("[Handler] Trying alternative: foundry download...")
            result2 = subprocess.run(
                ["foundry", "download", "--all"],
                capture_output=True,
                text=True,
                env=env,
                timeout=3600
            )
            print(f"[Handler] Alternative stdout: {result2.stdout}")
            if result2.stderr:
                print(f"[Handler] Alternative stderr: {result2.stderr}")

        # Get updated state
        after_info = get_checkpoint_info(CHECKPOINT_DIR)
        print(f"[Handler] After: {after_info['file_count']} files, {after_info['total_size_gb']} GB")

        downloaded_size = after_info["total_size_gb"] - before_info["total_size_gb"]

        return {
            "status": "completed",
            "result": {
                "checkpoints_dir": CHECKPOINT_DIR,
                "before": before_info,
                "after": after_info,
                "downloaded": True,
                "downloaded_gb": round(downloaded_size, 2),
                "message": f"Downloaded {round(downloaded_size, 2)} GB of checkpoints. Total: {after_info['total_size_gb']} GB",
                "foundry_output": result.stdout[:2000] if result.stdout else None
            }
        }

    except subprocess.TimeoutExpired:
        return {
            "status": "failed",
            "error": "Checkpoint download timed out after 1 hour"
        }
    except Exception as e:
        print(f"[Handler] Download error: {e}")
        traceback.print_exc()
        return {
            "status": "failed",
            "error": str(e),
            "traceback": traceback.format_exc()
        }


def handle_delete_file(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Delete a specific file from the Network Volume.
    Only allows deleting files within the checkpoint directory for safety.

    Input:
        filename: str - Name of file to delete (e.g., "rfd3_latest.ckpt")
    """
    filename = job_input.get("filename")
    if not filename:
        return {"status": "failed", "error": "Missing 'filename' parameter"}

    # Safety: only allow deleting from checkpoint directory
    filepath = os.path.join(CHECKPOINT_DIR, filename)

    # Prevent path traversal
    if not os.path.abspath(filepath).startswith(os.path.abspath(CHECKPOINT_DIR)):
        return {"status": "failed", "error": "Invalid filename - path traversal not allowed"}

    if not os.path.exists(filepath):
        return {
            "status": "completed",
            "result": {
                "deleted": False,
                "message": f"File does not exist: {filepath}"
            }
        }

    try:
        file_size = os.path.getsize(filepath)
        os.remove(filepath)
        return {
            "status": "completed",
            "result": {
                "deleted": True,
                "filepath": filepath,
                "size_bytes": file_size,
                "message": f"Successfully deleted {filepath} ({file_size} bytes)"
            }
        }
    except Exception as e:
        return {
            "status": "failed",
            "error": f"Failed to delete {filepath}: {str(e)}"
        }
