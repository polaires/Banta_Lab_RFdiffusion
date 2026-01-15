"""
Foundry API Backend for RunPod
Provides REST API endpoints for RFD3, RF3, and ProteinMPNN inference
"""

from fastapi import FastAPI, HTTPException, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Optional, List, Dict, Any
import subprocess
import tempfile
import os
import sys
import json
import uuid
from datetime import datetime

# Add serverless directory to path for shared module access
serverless_dir = os.path.join(os.path.dirname(__file__), '..', 'serverless')
if serverless_dir not in sys.path:
    sys.path.insert(0, serverless_dir)

# Import shared interaction analysis module (if available)
try:
    from shared.interaction_analysis import (
        analyze_all_interactions,
        format_for_frontend,
        format_for_ai_assistant,
        generate_recommendations,
    )
    INTERACTION_ANALYSIS_AVAILABLE = True
except ImportError:
    INTERACTION_ANALYSIS_AVAILABLE = False
    print("[API] Warning: shared.interaction_analysis not available")

app = FastAPI(
    title="Foundry Protein Design API",
    description="API for RFdiffusion3, RosettaFold3, and ProteinMPNN inference",
    version="1.0.0"
)

# CORS configuration for Vercel frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:3000",
        "http://localhost:5173",
        "https://*.vercel.app",
        "*"  # Remove in production
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# In-memory job storage (use Redis in production)
jobs: Dict[str, Dict[str, Any]] = {}


# ============ Request/Response Models ============

class RFD3Request(BaseModel):
    """RFdiffusion3 design request"""
    contigs: str  # e.g., "A1-50/0 50-100"
    num_designs: int = 1
    config: Optional[Dict[str, Any]] = None


class RF3Request(BaseModel):
    """RosettaFold3 prediction request"""
    sequence: str
    config: Optional[Dict[str, Any]] = None


class ProteinMPNNRequest(BaseModel):
    """ProteinMPNN/LigandMPNN sequence design request"""
    pdb_content: str
    num_sequences: int = 8
    temperature: float = 0.1
    model_type: str = "ligand_mpnn"  # "ligand_mpnn" or "protein_mpnn"
    remove_waters: bool = True
    config: Optional[Dict[str, Any]] = None


class JobResponse(BaseModel):
    """Job submission response"""
    job_id: str
    status: str
    message: str


class JobStatus(BaseModel):
    """Job status response"""
    job_id: str
    status: str  # "pending", "running", "completed", "failed"
    created_at: str
    completed_at: Optional[str] = None
    result: Optional[Dict[str, Any]] = None
    error: Optional[str] = None


class HealthResponse(BaseModel):
    """Health check response"""
    status: str
    gpu_available: bool
    models_loaded: List[str]


class InteractionAnalysisRequest(BaseModel):
    """Protein-ligand interaction analysis request"""
    pdb_content: str
    ligand_name: str = "UNL"
    include_visualization: bool = True
    include_recommendations: bool = True
    ligand_has_aromatics: bool = False


class InteractionAnalysisResponse(BaseModel):
    """Interaction analysis response"""
    status: str
    interactions: Dict[str, Any]
    key_residues: List[str]
    recommendations: Optional[List[str]] = None
    visualization: Optional[Dict[str, Any]] = None
    analysis_method: str
    ai_summary: Optional[str] = None


# ============ Helper Functions ============

def check_gpu_available() -> bool:
    """Check if CUDA GPU is available"""
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
            capture_output=True, text=True, timeout=10
        )
        return result.returncode == 0
    except Exception:
        return False


def get_loaded_models() -> List[str]:
    """Get list of available models"""
    checkpoint_dir = os.environ.get("FOUNDRY_CHECKPOINT_DIRS", "/workspace/checkpoints")
    models = []

    if os.path.exists(os.path.join(checkpoint_dir, "rfd3")):
        models.append("rfd3")
    if os.path.exists(os.path.join(checkpoint_dir, "rf3")):
        models.append("rf3")
    if os.path.exists(os.path.join(checkpoint_dir, "mpnn")):
        models.append("proteinmpnn")

    # If no checkpoints found, assume all models available (they'll be downloaded on first use)
    if not models:
        models = ["rfd3", "rf3", "proteinmpnn"]

    return models


# ============ Background Tasks ============

def run_rfd3_job(job_id: str, request: RFD3Request):
    """Run RFD3 design in background"""
    jobs[job_id]["status"] = "running"

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create input configuration
            input_config = {
                "contigs": [request.contigs],
                **(request.config or {})
            }
            input_path = os.path.join(tmpdir, "input.json")
            with open(input_path, "w") as f:
                json.dump(input_config, f)

            out_dir = os.path.join(tmpdir, "output")
            os.makedirs(out_dir, exist_ok=True)

            # Run RFD3
            cmd = f"rfd3 design out_dir={out_dir} inputs={input_path}"
            result = subprocess.run(
                cmd, shell=True, capture_output=True, text=True, timeout=3600
            )

            if result.returncode == 0:
                # Collect output PDB files
                outputs = []
                for filename in os.listdir(out_dir):
                    if filename.endswith(".pdb"):
                        filepath = os.path.join(out_dir, filename)
                        with open(filepath) as pdb_file:
                            outputs.append({
                                "filename": filename,
                                "content": pdb_file.read()
                            })

                jobs[job_id]["status"] = "completed"
                jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
                jobs[job_id]["result"] = {
                    "designs": outputs,
                    "stdout": result.stdout
                }
            else:
                jobs[job_id]["status"] = "failed"
                jobs[job_id]["error"] = result.stderr or "Unknown error"

    except subprocess.TimeoutExpired:
        jobs[job_id]["status"] = "failed"
        jobs[job_id]["error"] = "Job timed out (max 1 hour)"
    except Exception as e:
        jobs[job_id]["status"] = "failed"
        jobs[job_id]["error"] = str(e)


def run_rf3_job(job_id: str, request: RF3Request):
    """Run RF3 prediction in background"""
    jobs[job_id]["status"] = "running"

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            # Write sequence to FASTA file
            fasta_path = os.path.join(tmpdir, "input.fasta")
            with open(fasta_path, "w") as f:
                f.write(f">query\n{request.sequence}\n")

            out_dir = os.path.join(tmpdir, "output")
            os.makedirs(out_dir, exist_ok=True)

            # Run RF3
            cmd = f"rf3 predict out_dir={out_dir} inputs={fasta_path}"
            result = subprocess.run(
                cmd, shell=True, capture_output=True, text=True, timeout=3600
            )

            if result.returncode == 0:
                # Collect output
                outputs = []
                for filename in os.listdir(out_dir):
                    if filename.endswith(".pdb") or filename.endswith(".cif"):
                        filepath = os.path.join(out_dir, filename)
                        with open(filepath) as out_file:
                            outputs.append({
                                "filename": filename,
                                "content": out_file.read()
                            })

                jobs[job_id]["status"] = "completed"
                jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
                jobs[job_id]["result"] = {
                    "predictions": outputs,
                    "stdout": result.stdout
                }
            else:
                jobs[job_id]["status"] = "failed"
                jobs[job_id]["error"] = result.stderr or "Unknown error"

    except Exception as e:
        jobs[job_id]["status"] = "failed"
        jobs[job_id]["error"] = str(e)


def run_mpnn_job(job_id: str, request: ProteinMPNNRequest):
    """Run ProteinMPNN/LigandMPNN in background"""
    jobs[job_id]["status"] = "running"

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            # Write PDB file
            pdb_path = os.path.join(tmpdir, "input.pdb")
            with open(pdb_path, "w") as f:
                f.write(request.pdb_content)

            out_dir = os.path.join(tmpdir, "output")
            os.makedirs(out_dir, exist_ok=True)

            # Get model type from request
            model_type = request.model_type or 'ligand_mpnn'
            remove_waters = request.remove_waters

            # Run MPNN using Foundry CLI with correct arguments
            # CLI: mpnn --structure_path --out_directory --model_type --batch_size --temperature
            cmd_parts = [
                "mpnn",
                f"--structure_path {pdb_path}",
                f"--out_directory {out_dir}",
                "--name design",
                f"--model_type {model_type}",
                f"--batch_size {request.num_sequences}",
                f"--temperature {request.temperature}",
                "--is_legacy_weights True",
                "--write_fasta True",
                "--write_structures False",
            ]

            # Add remove_waters
            if remove_waters:
                cmd_parts.append("--remove_waters True")

            cmd = " ".join(cmd_parts)
            print(f"[MPNN] Running: {cmd}")

            result = subprocess.run(
                cmd, shell=True, capture_output=True, text=True, timeout=1800
            )

            if result.returncode == 0:
                # Collect designed sequences
                sequences = []
                for filename in os.listdir(out_dir):
                    if filename.endswith(".fa") or filename.endswith(".fasta"):
                        filepath = os.path.join(out_dir, filename)
                        with open(filepath) as seq_file:
                            sequences.append({
                                "filename": filename,
                                "content": seq_file.read()
                            })

                # Also check subdirectories
                if not sequences:
                    for root, dirs, files in os.walk(out_dir):
                        for filename in files:
                            if filename.endswith((".fa", ".fasta")):
                                with open(os.path.join(root, filename)) as f:
                                    sequences.append({"filename": filename, "content": f.read()})

                if sequences:
                    jobs[job_id]["status"] = "completed"
                    jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
                    jobs[job_id]["result"] = {
                        "sequences": sequences,
                        "model_type": model_type,
                        "stdout": result.stdout
                    }
                else:
                    jobs[job_id]["status"] = "failed"
                    jobs[job_id]["error"] = f"No sequences generated. stdout: {result.stdout}"
            else:
                jobs[job_id]["status"] = "failed"
                jobs[job_id]["error"] = result.stderr or result.stdout or "Unknown error"

    except Exception as e:
        import traceback
        jobs[job_id]["status"] = "failed"
        jobs[job_id]["error"] = f"{str(e)}\n{traceback.format_exc()}"


# ============ API Endpoints ============

@app.get("/", response_model=HealthResponse)
async def health_check():
    """Health check endpoint"""
    return HealthResponse(
        status="healthy",
        gpu_available=check_gpu_available(),
        models_loaded=get_loaded_models()
    )


@app.get("/health", response_model=HealthResponse)
async def health():
    """Alias for health check"""
    return await health_check()


@app.post("/api/rfd3/design", response_model=JobResponse)
async def submit_rfd3_design(request: RFD3Request, background_tasks: BackgroundTasks):
    """Submit RFdiffusion3 design job"""
    job_id = str(uuid.uuid4())
    jobs[job_id] = {
        "status": "pending",
        "created_at": datetime.utcnow().isoformat(),
        "type": "rfd3"
    }

    background_tasks.add_task(run_rfd3_job, job_id, request)

    return JobResponse(
        job_id=job_id,
        status="pending",
        message="RFD3 design job submitted successfully"
    )


@app.post("/api/rf3/predict", response_model=JobResponse)
async def submit_rf3_prediction(request: RF3Request, background_tasks: BackgroundTasks):
    """Submit RosettaFold3 prediction job"""
    job_id = str(uuid.uuid4())
    jobs[job_id] = {
        "status": "pending",
        "created_at": datetime.utcnow().isoformat(),
        "type": "rf3"
    }

    background_tasks.add_task(run_rf3_job, job_id, request)

    return JobResponse(
        job_id=job_id,
        status="pending",
        message="RF3 prediction job submitted successfully"
    )


@app.post("/api/mpnn/design", response_model=JobResponse)
async def submit_mpnn_design(request: ProteinMPNNRequest, background_tasks: BackgroundTasks):
    """Submit ProteinMPNN sequence design job"""
    job_id = str(uuid.uuid4())
    jobs[job_id] = {
        "status": "pending",
        "created_at": datetime.utcnow().isoformat(),
        "type": "mpnn"
    }

    background_tasks.add_task(run_mpnn_job, job_id, request)

    return JobResponse(
        job_id=job_id,
        status="pending",
        message="ProteinMPNN design job submitted successfully"
    )


@app.get("/api/jobs/{job_id}", response_model=JobStatus)
async def get_job_status(job_id: str):
    """Get job status and results"""
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")

    job = jobs[job_id]
    return JobStatus(
        job_id=job_id,
        status=job["status"],
        created_at=job["created_at"],
        completed_at=job.get("completed_at"),
        result=job.get("result"),
        error=job.get("error")
    )


@app.get("/api/jobs")
async def list_jobs():
    """List all jobs"""
    return {
        job_id: {
            "status": job["status"],
            "type": job.get("type"),
            "created_at": job["created_at"]
        }
        for job_id, job in jobs.items()
    }


@app.delete("/api/jobs/{job_id}")
async def delete_job(job_id: str):
    """Delete a job from memory"""
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")

    del jobs[job_id]
    return {"message": "Job deleted"}


@app.post("/api/analyze/interactions", response_model=InteractionAnalysisResponse)
async def analyze_interactions(request: InteractionAnalysisRequest):
    """
    Analyze protein-ligand interactions using PLIP.

    Detects multiple interaction types:
    - Hydrogen bonds (with D-H-A angle validation)
    - Hydrophobic contacts
    - Pi-stacking (face-to-face and edge-to-face)
    - Salt bridges
    - Halogen bonds

    Returns interaction data formatted for frontend visualization.
    """
    if not INTERACTION_ANALYSIS_AVAILABLE:
        raise HTTPException(
            status_code=503,
            detail="Interaction analysis module not available. Install with: pip install plip"
        )

    try:
        # Run comprehensive interaction analysis
        summary = analyze_all_interactions(
            pdb_content=request.pdb_content,
            ligand_name=request.ligand_name,
            include_visualization_data=request.include_visualization,
        )

        if summary.status == "error":
            raise HTTPException(
                status_code=500,
                detail=f"Analysis failed: {summary.error}"
            )

        # Format response
        frontend_data = format_for_frontend(summary)

        response = InteractionAnalysisResponse(
            status="completed",
            interactions=frontend_data["interactions"],
            key_residues=summary.key_residues,
            analysis_method=summary.analysis_method,
            visualization=frontend_data.get("visualization") if request.include_visualization else None,
        )

        # Add recommendations if requested
        if request.include_recommendations:
            response.recommendations = generate_recommendations(
                summary,
                ligand_has_aromatics=request.ligand_has_aromatics
            )

        # Add AI summary
        response.ai_summary = format_for_ai_assistant(summary)

        return response

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Interaction analysis failed: {str(e)}"
        )


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
