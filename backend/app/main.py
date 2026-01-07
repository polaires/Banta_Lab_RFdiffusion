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
import json
import uuid
from datetime import datetime

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
    """ProteinMPNN sequence design request"""
    pdb_content: str
    num_sequences: int = 8
    temperature: float = 0.1
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
    """Run ProteinMPNN in background"""
    jobs[job_id]["status"] = "running"

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            # Write PDB file
            pdb_path = os.path.join(tmpdir, "input.pdb")
            with open(pdb_path, "w") as f:
                f.write(request.pdb_content)

            out_dir = os.path.join(tmpdir, "output")
            os.makedirs(out_dir, exist_ok=True)

            # Run ProteinMPNN
            cmd = (
                f"proteinmpnn "
                f"--pdb {pdb_path} "
                f"--out_dir {out_dir} "
                f"--num_seq_per_target {request.num_sequences} "
                f"--sampling_temp {request.temperature}"
            )
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

                jobs[job_id]["status"] = "completed"
                jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
                jobs[job_id]["result"] = {
                    "sequences": sequences,
                    "stdout": result.stdout
                }
            else:
                jobs[job_id]["status"] = "failed"
                jobs[job_id]["error"] = result.stderr or "Unknown error"

    except Exception as e:
        jobs[job_id]["status"] = "failed"
        jobs[job_id]["error"] = str(e)


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


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
