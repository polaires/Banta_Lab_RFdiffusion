"""
Foundry Protein Design API - Unified Backend
FastAPI backend for RFdiffusion3, RosettaFold3, and ProteinMPNN

Features:
- Auto-detects Foundry CLI availability
- Falls back to mock mode when models unavailable
- Accepts both 'contig' and 'contigs' field names for compatibility

Deploy to RunPod:
1. Copy this file to /workspace/main.py on your RunPod pod
2. Run: uvicorn main:app --host 0.0.0.0 --port 8000
"""

from fastapi import FastAPI, HTTPException, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, model_validator
from typing import Optional, List, Dict, Any
import subprocess
import tempfile
import uuid
import time
import os
import json
import asyncio
from datetime import datetime
from contextlib import asynccontextmanager

# ============== Configuration ==============

# MOCK_MODE: "auto" (detect), "true" (force mock), "false" (force real)
MOCK_MODE = os.environ.get("MOCK_MODE", "auto").lower()
CHECKPOINT_DIR = os.environ.get("FOUNDRY_CHECKPOINT_DIRS", "/workspace/checkpoints")

# Global state
FOUNDRY_AVAILABLE = False
GPU_INFO: Dict[str, Any] = {}


def check_foundry_available() -> bool:
    """Check if Foundry CLI tools are installed"""
    # rc-foundry uses 'mpnn' not 'proteinmpnn'
    for cmd in ["rfd3", "rf3", "mpnn"]:
        try:
            result = subprocess.run([cmd, "--help"], capture_output=True, timeout=10)
            if result.returncode == 0:
                return True
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
    return False


def get_gpu_info() -> Dict[str, Any]:
    """Get GPU information from nvidia-smi"""
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=name,memory.total", "--format=csv,noheader,nounits"],
            capture_output=True, text=True, timeout=10
        )
        if result.returncode == 0:
            parts = result.stdout.strip().split(", ")
            return {
                "available": True,
                "name": parts[0] if parts else "Unknown",
                "memory_gb": float(parts[1]) / 1024 if len(parts) > 1 else 0,
            }
    except Exception:
        pass
    return {"available": False, "name": None, "memory_gb": None}


def get_model_status() -> Dict[str, Dict[str, Any]]:
    """Get detailed status of each model"""
    status = {}
    for model, subdir in [("rfd3", "rfd3"), ("rf3", "rf3"), ("proteinmpnn", "mpnn")]:
        path = os.path.join(CHECKPOINT_DIR, subdir)
        exists = os.path.exists(path)
        size_gb = 0
        if exists:
            try:
                total = sum(
                    os.path.getsize(os.path.join(dp, f))
                    for dp, _, filenames in os.walk(path)
                    for f in filenames
                )
                size_gb = total / (1024**3)
            except Exception:
                pass
        status[model] = {
            "available": FOUNDRY_AVAILABLE,
            "checkpoint_exists": exists,
            "checkpoint_size_gb": round(size_gb, 2),
        }
    return status


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Startup and shutdown events"""
    global FOUNDRY_AVAILABLE, GPU_INFO

    # Determine mode
    if MOCK_MODE == "true":
        FOUNDRY_AVAILABLE = False
        print("[STARTUP] Mock mode forced via MOCK_MODE=true")
    elif MOCK_MODE == "false":
        FOUNDRY_AVAILABLE = True
        print("[STARTUP] Real mode forced via MOCK_MODE=false")
    else:
        FOUNDRY_AVAILABLE = check_foundry_available()
        print(f"[STARTUP] Auto-detected Foundry available: {FOUNDRY_AVAILABLE}")

    GPU_INFO = get_gpu_info()
    print(f"[STARTUP] GPU: {GPU_INFO}")
    print(f"[STARTUP] Checkpoint dir: {CHECKPOINT_DIR}")

    yield


app = FastAPI(
    title="Foundry Protein Design API",
    description="Backend API for RFdiffusion3, RosettaFold3, and ProteinMPNN",
    version="1.0.0",
    lifespan=lifespan,
)

# CORS - allow all origins for development
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# In-memory job storage (use Redis/DB in production)
jobs: Dict[str, Dict[str, Any]] = {}

# ============== Models ==============


class HealthResponse(BaseModel):
    status: str
    mode: str  # "real" or "mock"
    gpu_available: bool
    gpu_name: Optional[str] = None
    gpu_memory_gb: Optional[float] = None
    models: Dict[str, Dict[str, Any]]


class RFD3Request(BaseModel):
    """RFdiffusion3 design request - accepts both field names"""
    contig: Optional[str] = None
    contigs: Optional[str] = None  # Alias for backward compat
    num_designs: int = 1
    pdb_content: Optional[str] = None
    config: Optional[Dict[str, Any]] = None

    @model_validator(mode="after")
    def normalize_contig(self):
        """Normalize contig/contigs to internal format"""
        if not self.contig and not self.contigs:
            raise ValueError("Either 'contig' or 'contigs' must be provided")
        if not self.contig and self.contigs:
            self.contig = self.contigs
        return self


class RF3Request(BaseModel):
    """RosettaFold3 prediction request"""
    sequence: str
    name: Optional[str] = "prediction"
    config: Optional[Dict[str, Any]] = None


class MPNNRequest(BaseModel):
    """ProteinMPNN sequence design request"""
    pdb_content: str
    num_sequences: int = 8
    temperature: float = 0.1
    config: Optional[Dict[str, Any]] = None


class JobResponse(BaseModel):
    job_id: str
    status: str
    message: str


class JobStatus(BaseModel):
    job_id: str
    status: str  # pending, running, completed, failed
    created_at: str
    completed_at: Optional[str] = None
    result: Optional[Dict[str, Any]] = None
    error: Optional[str] = None


# ============== Mock Implementations ==============


def generate_mock_pdb(length: int = 100) -> str:
    """Generate a mock PDB file for testing"""
    lines = ["HEADER    MOCK PROTEIN STRUCTURE"]
    for i in range(1, length + 1):
        x, y, z = i * 0.5, i * 0.3, i * 0.2
        lines.append(
            f"ATOM  {i:5d}  CA  ALA A{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
        )
    lines.append("END")
    return "\n".join(lines)


def generate_mock_fasta(num_seq: int = 8, length: int = 100) -> str:
    """Generate mock FASTA sequences for testing"""
    import random
    aa = "ACDEFGHIKLMNPQRSTVWY"
    lines = []
    for i in range(num_seq):
        seq = "".join(random.choices(aa, k=length))
        lines.append(f">design_{i+1}")
        lines.append(seq)
    return "\n".join(lines)


# ============== Job Processing ==============


async def process_rfd3_job(job_id: str, request: RFD3Request):
    """Process RFD3 job - uses real Foundry or mock"""
    jobs[job_id]["status"] = "running"

    if not FOUNDRY_AVAILABLE:
        await process_rfd3_mock(job_id, request)
    else:
        await process_rfd3_real(job_id, request)


async def process_rfd3_mock(job_id: str, request: RFD3Request):
    """Mock implementation for testing without Foundry"""
    await asyncio.sleep(3)  # Simulate processing

    try:
        length = int(request.contig) if request.contig and request.contig.isdigit() else 100
    except Exception:
        length = 100

    jobs[job_id]["status"] = "completed"
    jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
    jobs[job_id]["result"] = {
        "designs": [{"filename": "design_1.pdb", "content": generate_mock_pdb(length)}],
        "mode": "mock",
    }


async def process_rfd3_real(job_id: str, request: RFD3Request):
    """Real Foundry RFD3 implementation"""
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create input configuration
            input_config = {"contigs": [request.contig], **(request.config or {})}
            input_path = os.path.join(tmpdir, "input.json")
            with open(input_path, "w") as f:
                json.dump(input_config, f)

            out_dir = os.path.join(tmpdir, "output")
            os.makedirs(out_dir, exist_ok=True)

            # Run RFD3 with 1-hour timeout
            cmd = f"rfd3 design out_dir={out_dir} inputs={input_path}"
            result = subprocess.run(
                cmd, shell=True, capture_output=True, text=True, timeout=3600
            )

            if result.returncode == 0:
                outputs = []
                for filename in os.listdir(out_dir):
                    if filename.endswith(".pdb"):
                        with open(os.path.join(out_dir, filename)) as f:
                            outputs.append({"filename": filename, "content": f.read()})

                jobs[job_id]["status"] = "completed"
                jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
                jobs[job_id]["result"] = {
                    "designs": outputs,
                    "stdout": result.stdout,
                    "mode": "real",
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


async def process_rf3_job(job_id: str, request: RF3Request):
    """Process RF3 job - uses real Foundry or mock"""
    jobs[job_id]["status"] = "running"

    if not FOUNDRY_AVAILABLE:
        await process_rf3_mock(job_id, request)
    else:
        await process_rf3_real(job_id, request)


async def process_rf3_mock(job_id: str, request: RF3Request):
    """Mock implementation for RF3"""
    await asyncio.sleep(3)

    jobs[job_id]["status"] = "completed"
    jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
    jobs[job_id]["result"] = {
        "predictions": [{
            "filename": f"{request.name}.pdb",
            "content": generate_mock_pdb(len(request.sequence)),
        }],
        "mode": "mock",
    }


async def process_rf3_real(job_id: str, request: RF3Request):
    """Real Foundry RF3 implementation"""
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = os.path.join(tmpdir, "input.fasta")
            with open(fasta_path, "w") as f:
                f.write(f">query\n{request.sequence}\n")

            out_dir = os.path.join(tmpdir, "output")
            os.makedirs(out_dir, exist_ok=True)

            cmd = f"rf3 predict out_dir={out_dir} inputs={fasta_path}"
            result = subprocess.run(
                cmd, shell=True, capture_output=True, text=True, timeout=3600
            )

            if result.returncode == 0:
                outputs = []
                for filename in os.listdir(out_dir):
                    if filename.endswith((".pdb", ".cif")):
                        with open(os.path.join(out_dir, filename)) as f:
                            outputs.append({"filename": filename, "content": f.read()})

                jobs[job_id]["status"] = "completed"
                jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
                jobs[job_id]["result"] = {
                    "predictions": outputs,
                    "stdout": result.stdout,
                    "mode": "real",
                }
            else:
                jobs[job_id]["status"] = "failed"
                jobs[job_id]["error"] = result.stderr or "Unknown error"

    except Exception as e:
        jobs[job_id]["status"] = "failed"
        jobs[job_id]["error"] = str(e)


async def process_mpnn_job(job_id: str, request: MPNNRequest):
    """Process MPNN job - uses real Foundry or mock"""
    jobs[job_id]["status"] = "running"

    if not FOUNDRY_AVAILABLE:
        await process_mpnn_mock(job_id, request)
    else:
        await process_mpnn_real(job_id, request)


async def process_mpnn_mock(job_id: str, request: MPNNRequest):
    """Mock implementation for MPNN"""
    await asyncio.sleep(2)

    jobs[job_id]["status"] = "completed"
    jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
    jobs[job_id]["result"] = {
        "sequences": [{
            "filename": "sequences.fasta",
            "content": generate_mock_fasta(request.num_sequences, 100),
        }],
        "mode": "mock",
    }


async def process_mpnn_real(job_id: str, request: MPNNRequest):
    """Real Foundry MPNN implementation"""
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            pdb_path = os.path.join(tmpdir, "input.pdb")
            with open(pdb_path, "w") as f:
                f.write(request.pdb_content)

            out_dir = os.path.join(tmpdir, "output")
            os.makedirs(out_dir, exist_ok=True)

            # rc-foundry uses 'mpnn' command
            cmd = (
                f"mpnn "
                f"--pdb {pdb_path} "
                f"--out_dir {out_dir} "
                f"--num_seq_per_target {request.num_sequences} "
                f"--sampling_temp {request.temperature}"
            )
            result = subprocess.run(
                cmd, shell=True, capture_output=True, text=True, timeout=1800
            )

            if result.returncode == 0:
                sequences = []
                for filename in os.listdir(out_dir):
                    if filename.endswith((".fa", ".fasta")):
                        with open(os.path.join(out_dir, filename)) as f:
                            sequences.append({"filename": filename, "content": f.read()})

                jobs[job_id]["status"] = "completed"
                jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
                jobs[job_id]["result"] = {
                    "sequences": sequences,
                    "stdout": result.stdout,
                    "mode": "real",
                }
            else:
                jobs[job_id]["status"] = "failed"
                jobs[job_id]["error"] = result.stderr or "Unknown error"

    except Exception as e:
        jobs[job_id]["status"] = "failed"
        jobs[job_id]["error"] = str(e)


# ============== Endpoints ==============


@app.get("/", response_model=HealthResponse)
@app.get("/health", response_model=HealthResponse)
async def health():
    """Health check endpoint with detailed model status"""
    return HealthResponse(
        status="healthy",
        mode="real" if FOUNDRY_AVAILABLE else "mock",
        gpu_available=GPU_INFO.get("available", False),
        gpu_name=GPU_INFO.get("name"),
        gpu_memory_gb=GPU_INFO.get("memory_gb"),
        models=get_model_status(),
    )


@app.post("/api/rfd3/design", response_model=JobResponse)
async def submit_rfd3_design(request: RFD3Request, background_tasks: BackgroundTasks):
    """Submit RFdiffusion3 design job"""
    job_id = str(uuid.uuid4())

    jobs[job_id] = {
        "job_id": job_id,
        "status": "pending",
        "created_at": datetime.utcnow().isoformat(),
        "type": "rfd3",
        "request": {"contig": request.contig, "num_designs": request.num_designs},
    }

    background_tasks.add_task(process_rfd3_job, job_id, request)

    return JobResponse(
        job_id=job_id,
        status="pending",
        message="RFD3 design job submitted successfully",
    )


@app.post("/api/rf3/predict", response_model=JobResponse)
async def submit_rf3_prediction(request: RF3Request, background_tasks: BackgroundTasks):
    """Submit RosettaFold3 prediction job"""
    job_id = str(uuid.uuid4())

    jobs[job_id] = {
        "job_id": job_id,
        "status": "pending",
        "created_at": datetime.utcnow().isoformat(),
        "type": "rf3",
        "request": {"sequence_length": len(request.sequence)},
    }

    background_tasks.add_task(process_rf3_job, job_id, request)

    return JobResponse(
        job_id=job_id,
        status="pending",
        message="RF3 prediction job submitted successfully",
    )


@app.post("/api/mpnn/design", response_model=JobResponse)
async def submit_mpnn_design(request: MPNNRequest, background_tasks: BackgroundTasks):
    """Submit ProteinMPNN sequence design job"""
    job_id = str(uuid.uuid4())

    jobs[job_id] = {
        "job_id": job_id,
        "status": "pending",
        "created_at": datetime.utcnow().isoformat(),
        "type": "mpnn",
        "request": {"num_sequences": request.num_sequences, "temperature": request.temperature},
    }

    background_tasks.add_task(process_mpnn_job, job_id, request)

    return JobResponse(
        job_id=job_id,
        status="pending",
        message="ProteinMPNN design job submitted successfully",
    )


@app.get("/api/jobs/{job_id}", response_model=JobStatus)
async def get_job_status(job_id: str):
    """Get job status by ID"""
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")

    job = jobs[job_id]
    return JobStatus(
        job_id=job["job_id"],
        status=job["status"],
        created_at=job["created_at"],
        completed_at=job.get("completed_at"),
        result=job.get("result"),
        error=job.get("error"),
    )


@app.get("/api/jobs")
async def list_jobs():
    """List all jobs"""
    return {
        job_id: {
            "status": job["status"],
            "type": job.get("type"),
            "created_at": job["created_at"],
        }
        for job_id, job in jobs.items()
    }


@app.delete("/api/jobs/{job_id}")
async def delete_job(job_id: str):
    """Delete a job"""
    if job_id in jobs:
        del jobs[job_id]
    return {"status": "deleted"}


# ============== Main ==============

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
