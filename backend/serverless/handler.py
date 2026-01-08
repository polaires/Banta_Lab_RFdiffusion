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
        else:
            return {
                "status": "failed",
                "error": f"Unknown task: {task}. Valid tasks: health, rfd3, rf3, mpnn, rmsd"
            }

    except Exception as e:
        print(f"[Handler] Error: {e}")
        traceback.print_exc()
        return {
            "status": "failed",
            "error": str(e),
            "traceback": traceback.format_exc()
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
    contig = job_input.get("contig") or job_input.get("contigs")
    if not contig:
        return {"status": "failed", "error": "Missing 'contig' parameter"}

    num_designs = job_input.get("num_designs", 1)
    seed = job_input.get("seed")
    pdb_content = job_input.get("pdb_content")

    result = run_rfd3_inference(
        contig=contig,
        num_designs=num_designs,
        seed=seed,
        pdb_content=pdb_content,
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


# ============== RunPod Entry Point ==============

if __name__ == "__main__":
    print("[Handler] Starting RunPod serverless handler...")
    runpod.serverless.start({"handler": handler})
