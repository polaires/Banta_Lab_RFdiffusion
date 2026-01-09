"""
RunPod Serverless Handler for Foundry Protein Design

This handler processes inference requests for RFD3, RF3, and MPNN models
in a serverless environment. Models are loaded once at worker startup.

Usage:
    - Deploy to RunPod Serverless with network volume for checkpoints
    - Local testing: python handler.py --rp_serve_api
"""

import runpod
import os
import sys
import traceback
from typing import Dict, Any, Optional

# Import inference utilities
from inference_utils import (
    run_rfd3_inference,
    run_rf3_inference,
    run_mpnn_inference,
    calculate_rmsd,
    get_gpu_info,
    check_foundry_available,
)

# ============== Configuration ==============

CHECKPOINT_DIR = os.environ.get(
    "CHECKPOINT_DIR",
    os.environ.get("FOUNDRY_CHECKPOINT_DIRS", "/runpod-volume/checkpoints")
)
os.environ["FOUNDRY_CHECKPOINT_DIRS"] = CHECKPOINT_DIR

# Global model state
MODELS_LOADED = False
FOUNDRY_AVAILABLE = False
GPU_INFO: Dict[str, Any] = {}


# ============== Model Loading ==============

def load_models():
    """Load models at worker startup (called once per worker lifecycle)"""
    global MODELS_LOADED, FOUNDRY_AVAILABLE, GPU_INFO

    if MODELS_LOADED:
        print("[Handler] Models already loaded, skipping")
        return

    print(f"[Handler] Initializing worker...")
    print(f"[Handler] Checkpoint directory: {CHECKPOINT_DIR}")
    print(f"[Handler] Python version: {sys.version}")

    # Check GPU
    GPU_INFO = get_gpu_info()
    print(f"[Handler] GPU: {GPU_INFO}")

    # Check Foundry availability
    FOUNDRY_AVAILABLE = check_foundry_available()
    print(f"[Handler] Foundry available: {FOUNDRY_AVAILABLE}")

    if FOUNDRY_AVAILABLE:
        # Pre-import modules to warm up (optional, reduces first-request latency)
        try:
            print("[Handler] Pre-warming RFD3...")
            from rfd3.engine import RFD3InferenceEngine
            print("[Handler] RFD3 module loaded")
        except ImportError as e:
            print(f"[Handler] RFD3 not available: {e}")

        try:
            print("[Handler] Pre-warming RF3...")
            from rf3.inference_engines.rf3 import RF3InferenceEngine
            print("[Handler] RF3 module loaded")
        except ImportError as e:
            print(f"[Handler] RF3 not available: {e}")

        try:
            print("[Handler] Pre-warming MPNN...")
            from mpnn.inference_engines.mpnn import MPNNInferenceEngine
            print("[Handler] MPNN module loaded")
        except ImportError as e:
            print(f"[Handler] MPNN not available: {e}")

    MODELS_LOADED = True
    print("[Handler] Worker initialization complete")


# Load models at module import time (worker startup)
load_models()


# ============== Handler ==============

def handler(job: Dict[str, Any]) -> Dict[str, Any]:
    """
    Process a single inference request.

    Input format:
    {
        "input": {
            "task": "rfd3" | "rf3" | "mpnn" | "rmsd" | "health",
            // Task-specific parameters
        }
    }

    Output format:
    {
        "status": "completed" | "failed",
        "result": {...} | None,
        "error": str | None
    }
    """
    try:
        job_input = job.get("input", {})
        task = job_input.get("task", "").lower()

        print(f"[Handler] Processing task: {task}")

        if task == "health":
            return handle_health()
        elif task == "rfd3":
            return handle_rfd3(job_input)
        elif task == "rf3":
            return handle_rf3(job_input)
        elif task == "mpnn":
            return handle_mpnn(job_input)
        elif task == "rmsd":
            return handle_rmsd(job_input)
        elif task == "download_checkpoints":
            return handle_download_checkpoints(job_input)
        elif task == "delete_file":
            return handle_delete_file(job_input)
        else:
            return {
                "status": "failed",
                "error": f"Unknown task: {task}. Valid tasks: health, rfd3, rf3, mpnn, rmsd, download_checkpoints, delete_file"
            }

    except Exception as e:
        print(f"[Handler] Error: {e}")
        traceback.print_exc()

        # Collect detailed error context
        error_context = {
            "task": task,
            "input_keys": list(job_input.keys()) if job_input else [],
            "gpu_info": GPU_INFO,
            "foundry_available": FOUNDRY_AVAILABLE,
            "checkpoint_dir": CHECKPOINT_DIR,
        }

        # Try to get GPU memory status
        try:
            import subprocess
            mem_result = subprocess.run(
                ["nvidia-smi", "--query-gpu=memory.used,memory.total", "--format=csv,noheader,nounits"],
                capture_output=True, text=True, timeout=5
            )
            if mem_result.returncode == 0:
                parts = mem_result.stdout.strip().split(", ")
                error_context["gpu_memory_used_mb"] = int(parts[0]) if parts else None
                error_context["gpu_memory_total_mb"] = int(parts[1]) if len(parts) > 1 else None
        except Exception:
            pass

        return {
            "status": "failed",
            "error": str(e),
            "error_type": type(e).__name__,
            "traceback": traceback.format_exc(),
            "context": error_context
        }


def handle_health() -> Dict[str, Any]:
    """Return health status"""
    return {
        "status": "completed",
        "result": {
            "healthy": True,
            "mode": "real" if FOUNDRY_AVAILABLE else "mock",
            "gpu_available": GPU_INFO.get("available", False),
            "gpu_name": GPU_INFO.get("name"),
            "gpu_memory_gb": GPU_INFO.get("memory_gb"),
            "checkpoint_dir": CHECKPOINT_DIR,
        }
    }


def handle_rfd3(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """Handle RFD3 design request"""
    # Core parameters
    contig = job_input.get("contig") or job_input.get("contigs")
    length = job_input.get("length")

    # Either contig or length is required
    if not contig and not length:
        return {"status": "failed", "error": "Missing 'contig' or 'length' parameter"}

    num_designs = job_input.get("num_designs", 1)
    seed = job_input.get("seed")
    pdb_content = job_input.get("pdb_content")

    # Quality settings
    num_timesteps = job_input.get("num_timesteps")
    step_scale = job_input.get("step_scale")
    gamma_0 = job_input.get("gamma_0")
    is_non_loopy = job_input.get("is_non_loopy")

    # Partial diffusion (refinement)
    partial_t = job_input.get("partial_t")

    # Symmetry
    symmetry = job_input.get("symmetry")

    # Small molecule / enzyme design
    ligand = job_input.get("ligand")
    select_fixed_atoms = job_input.get("select_fixed_atoms")
    unindex = job_input.get("unindex")

    # RASA conditioning (binding pocket design)
    select_buried = job_input.get("select_buried")
    select_exposed = job_input.get("select_exposed")
    select_partially_buried = job_input.get("select_partially_buried")

    # Protein binder design
    hotspots = job_input.get("hotspots")
    infer_ori_strategy = job_input.get("infer_ori_strategy")

    # Nucleic acid binder design
    na_chains = job_input.get("na_chains")
    ori_token = job_input.get("ori_token")
    select_hbond_donor = job_input.get("select_hbond_donor")
    select_hbond_acceptor = job_input.get("select_hbond_acceptor")

    result = run_rfd3_inference(
        contig=contig,
        length=length,
        num_designs=num_designs,
        seed=seed,
        pdb_content=pdb_content,
        # Quality
        num_timesteps=num_timesteps,
        step_scale=step_scale,
        gamma_0=gamma_0,
        is_non_loopy=is_non_loopy,
        # Refinement
        partial_t=partial_t,
        # Symmetry
        symmetry=symmetry,
        # Small molecule / enzyme
        ligand=ligand,
        select_fixed_atoms=select_fixed_atoms,
        unindex=unindex,
        # RASA conditioning
        select_buried=select_buried,
        select_exposed=select_exposed,
        select_partially_buried=select_partially_buried,
        # Protein binder
        hotspots=hotspots,
        infer_ori_strategy=infer_ori_strategy,
        # Nucleic acid binder
        na_chains=na_chains,
        ori_token=ori_token,
        select_hbond_donor=select_hbond_donor,
        select_hbond_acceptor=select_hbond_acceptor,
        # Mock mode
        use_mock=not FOUNDRY_AVAILABLE
    )

    return result


def handle_rf3(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """Handle RF3 prediction request"""
    sequence = job_input.get("sequence")
    if not sequence:
        return {"status": "failed", "error": "Missing 'sequence' parameter"}

    name = job_input.get("name", "prediction")
    pdb_content = job_input.get("pdb_content")

    result = run_rf3_inference(
        sequence=sequence,
        name=name,
        pdb_content=pdb_content,
        use_mock=not FOUNDRY_AVAILABLE
    )

    return result


def handle_mpnn(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """Handle MPNN sequence design request"""
    pdb_content = job_input.get("pdb_content")
    if not pdb_content:
        return {"status": "failed", "error": "Missing 'pdb_content' parameter"}

    num_sequences = job_input.get("num_sequences", 8)
    temperature = job_input.get("temperature", 0.1)
    model_type = job_input.get("model_type", "ligand_mpnn")
    remove_waters = job_input.get("remove_waters", True)

    result = run_mpnn_inference(
        pdb_content=pdb_content,
        num_sequences=num_sequences,
        temperature=temperature,
        model_type=model_type,
        remove_waters=remove_waters,
        use_mock=not FOUNDRY_AVAILABLE
    )

    return result


def handle_rmsd(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """Handle RMSD calculation request"""
    pdb1 = job_input.get("pdb_content_1")
    pdb2 = job_input.get("pdb_content_2")

    if not pdb1 or not pdb2:
        return {"status": "failed", "error": "Missing pdb_content_1 or pdb_content_2"}

    backbone_only = job_input.get("backbone_only", True)

    result = calculate_rmsd(pdb1, pdb2, backbone_only)

    if "error" in result and result.get("rmsd") is None:
        return {"status": "failed", "error": result["error"]}

    return {"status": "completed", "result": result}


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
    import subprocess
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


# ============== RunPod Entry Point ==============

if __name__ == "__main__":
    print("[Handler] Starting RunPod serverless handler...")
    runpod.serverless.start({"handler": handler})
