"""
Foundry Protein Design API - Unified Backend (v2.0)
FastAPI backend for RFdiffusion3, RosettaFold3, and ProteinMPNN

Features:
- Python API integration for all models (matching IPD official pipeline)
- RF3 confidence metrics (pLDDT, PAE, pTM, ranking_score)
- MPNN model type selection (ligand_mpnn/protein_mpnn)
- RMSD validation endpoint for designability assessment
- Seed control for reproducible generation
- CIF export support
- Falls back to mock mode when models unavailable

Deploy to RunPod:
1. Copy this file to /workspace/main.py on your RunPod pod
2. Run: uvicorn main:app --host 0.0.0.0 --port 8000
"""

from fastapi import FastAPI, HTTPException, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, model_validator, Field
from typing import Optional, List, Dict, Any, Literal
import subprocess
import tempfile
import uuid
import os
import json
import asyncio
from datetime import datetime
from contextlib import asynccontextmanager

# ============== Configuration ==============

MOCK_MODE = os.environ.get("MOCK_MODE", "auto").lower()
CHECKPOINT_DIR = os.environ.get("FOUNDRY_CHECKPOINT_DIRS", "/workspace/checkpoints")

import sys
VENV_BIN_DIR = os.path.dirname(sys.executable)

KNOWN_VENV_PATHS = [
    VENV_BIN_DIR,
    "/workspace/foundry_env/bin",
    os.path.expanduser("~/foundry_env/bin"),
]

CLI_PREFIX = ""
for venv_path in KNOWN_VENV_PATHS:
    if os.path.exists(os.path.join(venv_path, "rfd3")):
        CLI_PREFIX = venv_path
        break

# Global state
FOUNDRY_AVAILABLE = False
GPU_INFO: Dict[str, Any] = {}

# In-memory job storage
jobs: Dict[str, Dict[str, Any]] = {}

# Store atom arrays for cross-panel data flow
stored_structures: Dict[str, Any] = {}


def get_cli_path(cmd: str) -> str:
    if CLI_PREFIX:
        return os.path.join(CLI_PREFIX, cmd)
    return cmd


def check_foundry_available() -> bool:
    """Check if Foundry Python modules are importable"""
    try:
        from rfd3.engine import RFD3InferenceEngine
        print("[STARTUP] Found rfd3.engine module")
        return True
    except ImportError:
        pass

    # Fallback: check CLI tools
    for cmd in ["rfd3", "rf3", "mpnn"]:
        try:
            cli_path = get_cli_path(cmd)
            result = subprocess.run([cli_path, "--help"], capture_output=True, timeout=10)
            if result.returncode == 0:
                print(f"[STARTUP] Found CLI tool: {cli_path}")
                return True
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
    return False


def get_gpu_info() -> Dict[str, Any]:
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
    global FOUNDRY_AVAILABLE, GPU_INFO

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
    print(f"[STARTUP] CLI prefix: {CLI_PREFIX or '(system PATH)'}")

    yield


app = FastAPI(
    title="Foundry Protein Design API",
    description="Backend API for RFdiffusion3, RosettaFold3, and ProteinMPNN (v2.0 - Python API)",
    version="2.0.0",
    lifespan=lifespan,
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# ============== Request/Response Models ==============

class HealthResponse(BaseModel):
    status: str
    mode: str
    gpu_available: bool
    gpu_name: Optional[str] = None
    gpu_memory_gb: Optional[float] = None
    models: Dict[str, Dict[str, Any]]


class RFD3Request(BaseModel):
    """RFdiffusion3 design request with seed control"""
    contig: Optional[str] = None
    contigs: Optional[str] = None
    num_designs: int = Field(default=1, ge=1, le=10)
    pdb_content: Optional[str] = None
    seed: Optional[int] = None  # For reproducibility
    config: Optional[Dict[str, Any]] = None

    @model_validator(mode="after")
    def normalize_contig(self):
        if not self.contig and not self.contigs:
            raise ValueError("Either 'contig' or 'contigs' must be provided")
        if not self.contig and self.contigs:
            self.contig = self.contigs
        return self


class RF3Request(BaseModel):
    """RosettaFold3 prediction request"""
    sequence: str
    name: Optional[str] = "prediction"
    pdb_content: Optional[str] = None  # For structure-based input
    config: Optional[Dict[str, Any]] = None


class MPNNRequest(BaseModel):
    """ProteinMPNN sequence design request with model type selection"""
    pdb_content: str
    num_sequences: int = Field(default=8, ge=1, le=100)
    temperature: float = Field(default=0.1, ge=0.01, le=2.0)
    model_type: Literal["ligand_mpnn", "protein_mpnn"] = "ligand_mpnn"
    remove_waters: bool = True
    config: Optional[Dict[str, Any]] = None


class RMSDRequest(BaseModel):
    """RMSD validation request"""
    pdb_content_1: str  # Reference structure (e.g., RFD3 output)
    pdb_content_2: str  # Comparison structure (e.g., RF3 refolded)
    backbone_only: bool = True


class ExportRequest(BaseModel):
    """Structure export request"""
    pdb_content: str
    format: Literal["pdb", "cif"] = "cif"
    include_confidences: bool = False
    confidences: Optional[Dict[str, Any]] = None


class JobResponse(BaseModel):
    job_id: str
    status: str
    message: str


class JobStatus(BaseModel):
    job_id: str
    status: str
    created_at: str
    completed_at: Optional[str] = None
    result: Optional[Dict[str, Any]] = None
    error: Optional[str] = None


# ============== Utility Functions ==============

def atom_array_to_pdb(atom_array) -> str:
    """Convert biotite AtomArray to PDB string"""
    from biotite.structure.io.pdb import PDBFile
    import io

    pdb_file = PDBFile()
    pdb_file.set_structure(atom_array)
    buf = io.StringIO()
    pdb_file.write(buf)
    return buf.getvalue()


def atom_array_to_cif(atom_array) -> str:
    """Convert biotite AtomArray to CIF string"""
    try:
        from atomworks.io.utils.io_utils import to_cif_file
        import tempfile
        import os

        with tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as f:
            temp_path = f.name

        to_cif_file(atom_array, temp_path)
        with open(temp_path) as f:
            content = f.read()
        os.unlink(temp_path)
        return content
    except ImportError:
        # Fallback to biotite CIF writer
        from biotite.structure.io.pdbx import PDBxFile, set_structure
        import io

        pdbx_file = PDBxFile()
        set_structure(pdbx_file, atom_array)
        buf = io.StringIO()
        pdbx_file.write(buf)
        return buf.getvalue()


def pdb_to_atom_array(pdb_content: str):
    """Convert PDB string to biotite AtomArray"""
    from biotite.structure.io.pdb import PDBFile
    import io

    pdb_file = PDBFile.read(io.StringIO(pdb_content))
    return pdb_file.get_structure(model=1)


def extract_sequence_from_atom_array(atom_array) -> str:
    """Extract amino acid sequence from AtomArray"""
    from biotite.structure import get_residue_starts
    from biotite.sequence import ProteinSequence

    res_starts = get_residue_starts(atom_array)
    seq = ''.join(
        ProteinSequence.convert_letter_3to1(res_name)
        for res_name in atom_array.res_name[res_starts]
        if res_name in ProteinSequence._dict_3to1
    )
    return seq


# ============== Mock Implementations ==============

def generate_mock_pdb(length: int = 100) -> str:
    lines = ["HEADER    MOCK PROTEIN STRUCTURE"]
    for i in range(1, length + 1):
        x, y, z = i * 0.5, i * 0.3, i * 0.2
        lines.append(
            f"ATOM  {i:5d}  CA  ALA A{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
        )
    lines.append("END")
    return "\n".join(lines)


def generate_mock_fasta(num_seq: int = 8, length: int = 100) -> str:
    import random
    aa = "ACDEFGHIKLMNPQRSTVWY"
    lines = []
    for i in range(num_seq):
        seq = "".join(random.choices(aa, k=length))
        lines.append(f">design_{i+1}")
        lines.append(seq)
    return "\n".join(lines)


def generate_mock_confidences(length: int = 100) -> Dict[str, Any]:
    """Generate mock confidence metrics"""
    import random
    return {
        "summary_confidences": {
            "overall_plddt": round(random.uniform(0.7, 0.95), 3),
            "overall_pae": round(random.uniform(2.0, 8.0), 2),
            "overall_pde": round(random.uniform(0.1, 0.3), 3),
            "ptm": round(random.uniform(0.6, 0.9), 3),
            "iptm": None,
            "ranking_score": round(random.uniform(0.5, 0.85), 3),
            "has_clash": False,
        },
        "per_residue_plddt": [round(random.uniform(0.6, 1.0), 2) for _ in range(length)],
    }


# ============== RFD3 Job Processing ==============

async def process_rfd3_job(job_id: str, request: RFD3Request):
    jobs[job_id]["status"] = "running"

    if not FOUNDRY_AVAILABLE:
        await process_rfd3_mock(job_id, request)
    else:
        await process_rfd3_real(job_id, request)


async def process_rfd3_mock(job_id: str, request: RFD3Request):
    await asyncio.sleep(3)

    try:
        length = int(request.contig) if request.contig and request.contig.isdigit() else 100
    except Exception:
        length = 100

    pdb_content = generate_mock_pdb(length)

    jobs[job_id]["status"] = "completed"
    jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
    jobs[job_id]["result"] = {
        "designs": [{"filename": "design_1.pdb", "content": pdb_content}],
        "mode": "mock",
    }

    # Store for data flow
    stored_structures[job_id] = {"pdb": pdb_content, "type": "rfd3"}


async def process_rfd3_real(job_id: str, request: RFD3Request):
    try:
        from rfd3.engine import RFD3InferenceConfig, RFD3InferenceEngine

        # Set seed for reproducibility
        if request.seed is not None:
            from lightning.fabric import seed_everything
            seed_everything(request.seed)
            print(f"[RFD3] Set seed to {request.seed}")

        contig_str = request.contig or "100"
        is_simple_length = contig_str.replace("-", "").isdigit()

        if is_simple_length:
            if "-" in contig_str:
                spec = {"length": contig_str}
            else:
                spec = {"length": int(contig_str)}
        else:
            spec = {"contig": contig_str}

        config = RFD3InferenceConfig(
            specification=spec,
            diffusion_batch_size=request.num_designs,
        )

        model = RFD3InferenceEngine(**config)
        outputs_dict = model.run(
            inputs=None,
            out_dir=None,
            n_batches=1,
        )

        outputs = []
        first_atom_array = None

        if outputs_dict:
            for key, result_list in outputs_dict.items():
                items = result_list if isinstance(result_list, list) else [result_list]

                for idx, item in enumerate(items):
                    pdb_content = None
                    cif_content = None
                    filename = f"{key}_{idx}.pdb" if len(items) > 1 else f"{key}.pdb"

                    if hasattr(item, 'atom_array'):
                        try:
                            atom_array = item.atom_array
                            if first_atom_array is None:
                                first_atom_array = atom_array

                            pdb_content = atom_array_to_pdb(atom_array)
                            try:
                                cif_content = atom_array_to_cif(atom_array)
                            except Exception:
                                pass
                        except Exception as e:
                            print(f"Error converting atom_array: {e}")

                    if pdb_content:
                        output_item = {"filename": filename, "content": pdb_content}
                        if cif_content:
                            output_item["cif_content"] = cif_content
                        outputs.append(output_item)

        jobs[job_id]["status"] = "completed"
        jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
        jobs[job_id]["result"] = {
            "designs": outputs,
            "mode": "real",
            "seed": request.seed,
        }

        # Store first structure for data flow
        if outputs:
            stored_structures[job_id] = {
                "pdb": outputs[0]["content"],
                "cif": outputs[0].get("cif_content"),
                "type": "rfd3",
            }

    except ImportError as e:
        jobs[job_id]["status"] = "failed"
        jobs[job_id]["error"] = f"Foundry not installed: {e}"
    except Exception as e:
        import traceback
        jobs[job_id]["status"] = "failed"
        jobs[job_id]["error"] = f"{str(e)}\n{traceback.format_exc()}"


# ============== RF3 Job Processing (Python API with Confidences) ==============

async def process_rf3_job(job_id: str, request: RF3Request):
    jobs[job_id]["status"] = "running"

    if not FOUNDRY_AVAILABLE:
        await process_rf3_mock(job_id, request)
    else:
        await process_rf3_real(job_id, request)


async def process_rf3_mock(job_id: str, request: RF3Request):
    await asyncio.sleep(3)

    length = len(request.sequence)
    pdb_content = generate_mock_pdb(length)
    confidences = generate_mock_confidences(length)

    jobs[job_id]["status"] = "completed"
    jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
    jobs[job_id]["result"] = {
        "predictions": [{
            "filename": f"{request.name}.pdb",
            "content": pdb_content,
        }],
        "confidences": confidences,
        "mode": "mock",
    }

    stored_structures[job_id] = {"pdb": pdb_content, "type": "rf3"}


async def process_rf3_real(job_id: str, request: RF3Request):
    """Real RF3 implementation using Python API with confidence metrics"""
    try:
        from rf3.inference_engines.rf3 import RF3InferenceEngine
        from rf3.utils.inference import InferenceInput

        # For sequence-only input, use CLI approach (Python API needs atom_array)
        if not request.pdb_content:
            print("[RF3] Sequence-only input, using CLI approach")
            await process_rf3_cli(job_id, request)
            return

        # Initialize RF3 engine for structure input
        inference_engine = RF3InferenceEngine(ckpt_path='rf3', verbose=False)

        # Structure-based input
        atom_array = pdb_to_atom_array(request.pdb_content)
        input_data = InferenceInput.from_atom_array(
            atom_array,
            example_id=request.name or "prediction"
        )

        # Run inference
        rf3_outputs = inference_engine.run(inputs=input_data)

        # Process outputs
        outputs = []
        confidences = None
        first_pdb = None

        example_id = request.name or "prediction"
        if example_id in rf3_outputs:
            results = rf3_outputs[example_id]

            for idx, rf3_output in enumerate(results):
                # Extract structure
                if hasattr(rf3_output, 'atom_array'):
                    pdb_content = atom_array_to_pdb(rf3_output.atom_array)
                    if first_pdb is None:
                        first_pdb = pdb_content

                    output_item = {
                        "filename": f"{example_id}_{idx}.pdb" if len(results) > 1 else f"{example_id}.pdb",
                        "content": pdb_content,
                    }

                    # Try to add CIF
                    try:
                        output_item["cif_content"] = atom_array_to_cif(rf3_output.atom_array)
                    except Exception:
                        pass

                    outputs.append(output_item)

                # Extract confidence metrics (from first output)
                if idx == 0 and confidences is None:
                    confidences = {}

                    # Summary confidences
                    if hasattr(rf3_output, 'summary_confidences') and rf3_output.summary_confidences:
                        sc = rf3_output.summary_confidences
                        confidences["summary_confidences"] = {
                            "overall_plddt": float(sc.get('overall_plddt', 0)),
                            "overall_pae": float(sc.get('overall_pae', 0)),
                            "overall_pde": float(sc.get('overall_pde', 0)) if 'overall_pde' in sc else None,
                            "ptm": float(sc.get('ptm', 0)),
                            "iptm": float(sc.get('iptm')) if sc.get('iptm') is not None else None,
                            "ranking_score": float(sc.get('ranking_score', 0)),
                            "has_clash": bool(sc.get('has_clash', False)),
                        }

                    # Per-residue confidences
                    if hasattr(rf3_output, 'confidences') and rf3_output.confidences:
                        conf = rf3_output.confidences
                        if 'atom_plddts' in conf:
                            # Convert to per-residue by averaging
                            atom_plddts = list(conf['atom_plddts'])
                            confidences["per_residue_plddt"] = [round(float(p), 3) for p in atom_plddts[:500]]  # Limit size

                        if 'pae' in conf:
                            # PAE matrix - limit size for JSON
                            pae = conf['pae']
                            pae_list = [[round(float(x), 2) for x in row[:100]] for row in pae[:100]]
                            confidences["pae_matrix"] = pae_list

        jobs[job_id]["status"] = "completed"
        jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
        jobs[job_id]["result"] = {
            "predictions": outputs,
            "confidences": confidences,
            "mode": "real",
        }

        if first_pdb:
            stored_structures[job_id] = {"pdb": first_pdb, "type": "rf3", "confidences": confidences}

    except ImportError as e:
        # Fallback to CLI if Python API not available
        print(f"[RF3] Python API not available ({e}), falling back to CLI")
        await process_rf3_cli(job_id, request)
    except Exception as e:
        import traceback
        jobs[job_id]["status"] = "failed"
        jobs[job_id]["error"] = f"{str(e)}\n{traceback.format_exc()}"


async def process_rf3_cli(job_id: str, request: RF3Request):
    """Fallback CLI implementation for RF3"""
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = os.path.join(tmpdir, "input.fasta")
            with open(fasta_path, "w") as f:
                f.write(f">query\n{request.sequence}\n")

            out_dir = os.path.join(tmpdir, "output")
            os.makedirs(out_dir, exist_ok=True)

            # Create .project-root file in /workspace for rootutils
            # IMPORTANT: Must be in /workspace, NOT in tmpdir, otherwise RF3 looks for
            # checkpoints in the wrong location
            workspace_root = "/workspace/.project-root"
            if not os.path.exists(workspace_root):
                try:
                    with open(workspace_root, "w") as f:
                        f.write("")
                    print(f"[RF3] Created {workspace_root}")
                except Exception as e:
                    print(f"[RF3] Could not create {workspace_root}: {e}")

            # Set up environment with checkpoint directory
            env = os.environ.copy()
            env["FOUNDRY_CHECKPOINT_DIRS"] = "/workspace/checkpoints"

            rf3_cli = get_cli_path("rf3")
            cmd = f"{rf3_cli} predict out_dir={out_dir} inputs={fasta_path}"
            print(f"[RF3] Running: {cmd}")
            print(f"[RF3] FOUNDRY_CHECKPOINT_DIRS={env.get('FOUNDRY_CHECKPOINT_DIRS')}")
            result = subprocess.run(
                cmd, shell=True, capture_output=True, text=True, timeout=3600,
                cwd="/workspace", env=env
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
                    "confidences": None,  # Not available from CLI
                    "mode": "real (cli fallback)",
                }
            else:
                jobs[job_id]["status"] = "failed"
                jobs[job_id]["error"] = result.stderr or "Unknown error"

    except Exception as e:
        jobs[job_id]["status"] = "failed"
        jobs[job_id]["error"] = str(e)


# ============== MPNN Job Processing (Python API with model type) ==============

async def process_mpnn_job(job_id: str, request: MPNNRequest):
    jobs[job_id]["status"] = "running"

    if not FOUNDRY_AVAILABLE:
        await process_mpnn_mock(job_id, request)
    else:
        await process_mpnn_real(job_id, request)


async def process_mpnn_mock(job_id: str, request: MPNNRequest):
    await asyncio.sleep(2)

    # Parse PDB to get approximate length
    lines = request.pdb_content.split('\n')
    ca_count = sum(1 for l in lines if l.startswith('ATOM') and ' CA ' in l)
    length = ca_count if ca_count > 0 else 100

    jobs[job_id]["status"] = "completed"
    jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
    jobs[job_id]["result"] = {
        "sequences": [{
            "filename": "sequences.fasta",
            "content": generate_mock_fasta(request.num_sequences, length),
        }],
        "model_type": request.model_type,
        "mode": "mock",
    }


async def process_mpnn_real(job_id: str, request: MPNNRequest):
    """Real MPNN implementation using Python API"""
    try:
        from mpnn.inference_engines.mpnn import MPNNInferenceEngine
        from biotite.structure import get_residue_starts
        from biotite.sequence import ProteinSequence

        # Parse input PDB to atom array
        atom_array = pdb_to_atom_array(request.pdb_content)

        # Configure MPNN engine (per IPD official example)
        engine_config = {
            "model_type": request.model_type,  # "ligand_mpnn" or "protein_mpnn"
            "is_legacy_weights": True,  # Required for current models
            "out_directory": None,  # Return in memory
            "write_structures": False,
            "write_fasta": False,
        }

        # Per-input config
        input_configs = [{
            "batch_size": request.num_sequences,
            "remove_waters": request.remove_waters,
        }]

        # Run MPNN
        model = MPNNInferenceEngine(**engine_config)
        mpnn_outputs = model.run(input_dicts=input_configs, atom_arrays=[atom_array])

        # Extract sequences from outputs
        sequences = []
        for i, item in enumerate(mpnn_outputs):
            if hasattr(item, 'atom_array'):
                # Extract sequence from atom array
                res_starts = get_residue_starts(item.atom_array)
                seq = ''.join(
                    ProteinSequence.convert_letter_3to1(res_name)
                    for res_name in item.atom_array.res_name[res_starts]
                    if res_name in ProteinSequence._dict_3to1
                )
                sequences.append(f">design_{i+1}\n{seq}")

        fasta_content = "\n".join(sequences)

        jobs[job_id]["status"] = "completed"
        jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
        jobs[job_id]["result"] = {
            "sequences": [{
                "filename": "sequences.fasta",
                "content": fasta_content,
            }],
            "num_sequences": len(sequences),
            "model_type": request.model_type,
            "mode": "real",
        }

    except ImportError as e:
        print(f"[MPNN] Python API not available ({e}), falling back to CLI")
        await process_mpnn_cli(job_id, request)
    except Exception as e:
        # Try CLI fallback on any Python API error
        import traceback
        print(f"[MPNN] Python API failed ({e}), falling back to CLI\n{traceback.format_exc()}")
        try:
            await process_mpnn_cli(job_id, request)
        except Exception as cli_error:
            jobs[job_id]["status"] = "failed"
            jobs[job_id]["error"] = f"Python API: {str(e)}\nCLI fallback: {str(cli_error)}"


async def process_mpnn_cli(job_id: str, request: MPNNRequest):
    """Fallback CLI implementation for MPNN"""
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            pdb_path = os.path.join(tmpdir, "input.pdb")
            with open(pdb_path, "w") as f:
                f.write(request.pdb_content)

            out_dir = os.path.join(tmpdir, "output")
            os.makedirs(out_dir, exist_ok=True)

            mpnn_cli = get_cli_path("mpnn")

            # Build command with correct MPNN CLI arguments
            # CLI uses: --structure_path, --out_directory, --batch_size, --temperature
            # --model_type (ligand_mpnn or protein_mpnn), --is_legacy_weights True
            cmd_parts = [
                mpnn_cli,
                f"--structure_path {pdb_path}",
                f"--out_directory {out_dir}",
                "--name design",
                f"--model_type {request.model_type}",
                f"--batch_size {request.num_sequences}",
                f"--temperature {request.temperature}",
                "--is_legacy_weights True",
                "--write_fasta True",
                "--write_structures False",
            ]

            # Add remove_waters if specified (use True/False/None format)
            if request.remove_waters:
                cmd_parts.append("--remove_waters True")
            else:
                cmd_parts.append("--remove_waters None")

            cmd = " ".join(cmd_parts)
            print(f"[MPNN CLI] Running: {cmd}")

            result = subprocess.run(
                cmd, shell=True, capture_output=True, text=True, timeout=1800
            )

            if result.returncode == 0:
                sequences = []
                # MPNN outputs FASTA files
                for filename in sorted(os.listdir(out_dir)):
                    if filename.endswith((".fa", ".fasta")):
                        with open(os.path.join(out_dir, filename)) as f:
                            sequences.append({"filename": filename, "content": f.read()})

                # Also check for sequences in subdirectories
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
                        "model_type": request.model_type,
                        "mode": "real (cli)",
                    }
                else:
                    jobs[job_id]["status"] = "failed"
                    jobs[job_id]["error"] = f"No output sequences found. stdout: {result.stdout}"
            else:
                jobs[job_id]["status"] = "failed"
                jobs[job_id]["error"] = result.stderr or result.stdout or "Unknown error"

    except Exception as e:
        import traceback
        jobs[job_id]["status"] = "failed"
        jobs[job_id]["error"] = f"{str(e)}\n{traceback.format_exc()}"


# ============== RMSD Validation ==============

def calculate_rmsd(pdb1: str, pdb2: str, backbone_only: bool = True) -> Dict[str, Any]:
    """Calculate RMSD between two structures"""
    try:
        from biotite.structure import rmsd, superimpose
        import numpy as np

        # Backbone atom names
        BACKBONE_ATOMS = ['N', 'CA', 'C', 'O']

        # Parse structures
        aa1 = pdb_to_atom_array(pdb1)
        aa2 = pdb_to_atom_array(pdb2)

        if backbone_only:
            # Filter to backbone atoms
            mask1 = np.isin(aa1.atom_name, BACKBONE_ATOMS)
            mask2 = np.isin(aa2.atom_name, BACKBONE_ATOMS)
            aa1 = aa1[mask1]
            aa2 = aa2[mask2]

        # Check if structures have same number of atoms
        if len(aa1) != len(aa2):
            # Try to align by residue
            min_len = min(len(aa1), len(aa2))
            aa1 = aa1[:min_len]
            aa2 = aa2[:min_len]

        # Superimpose and calculate RMSD
        aa2_fitted, _ = superimpose(aa1, aa2)
        rmsd_value = float(rmsd(aa1, aa2_fitted))

        # Interpret RMSD
        if rmsd_value < 1.0:
            interpretation = "Excellent"
            description = "Very high designability - structure is highly likely to fold as designed"
        elif rmsd_value < 2.0:
            interpretation = "Good"
            description = "Good designability - structure will likely fold correctly"
        elif rmsd_value < 3.0:
            interpretation = "Moderate"
            description = "Moderate designability - some structural deviation expected"
        else:
            interpretation = "Poor"
            description = "Low designability - significant structural deviation"

        return {
            "rmsd": round(rmsd_value, 3),
            "interpretation": interpretation,
            "description": description,
            "backbone_only": backbone_only,
            "num_atoms_compared": len(aa1),
        }

    except Exception as e:
        return {
            "error": str(e),
            "rmsd": None,
        }


# ============== Endpoints ==============

@app.get("/", response_model=HealthResponse)
@app.get("/health", response_model=HealthResponse)
async def health():
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
    job_id = str(uuid.uuid4())

    jobs[job_id] = {
        "job_id": job_id,
        "status": "pending",
        "created_at": datetime.utcnow().isoformat(),
        "type": "rfd3",
        "request": {"contig": request.contig, "num_designs": request.num_designs, "seed": request.seed},
    }

    background_tasks.add_task(process_rfd3_job, job_id, request)

    return JobResponse(
        job_id=job_id,
        status="pending",
        message="RFD3 design job submitted successfully",
    )


@app.post("/api/rf3/predict", response_model=JobResponse)
async def submit_rf3_prediction(request: RF3Request, background_tasks: BackgroundTasks):
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
    job_id = str(uuid.uuid4())

    jobs[job_id] = {
        "job_id": job_id,
        "status": "pending",
        "created_at": datetime.utcnow().isoformat(),
        "type": "mpnn",
        "request": {
            "num_sequences": request.num_sequences,
            "temperature": request.temperature,
            "model_type": request.model_type,
        },
    }

    background_tasks.add_task(process_mpnn_job, job_id, request)

    return JobResponse(
        job_id=job_id,
        status="pending",
        message="ProteinMPNN design job submitted successfully",
    )


@app.post("/api/validate/rmsd")
async def validate_rmsd(request: RMSDRequest):
    """Calculate RMSD between two structures for designability validation"""
    result = calculate_rmsd(
        request.pdb_content_1,
        request.pdb_content_2,
        request.backbone_only
    )

    if "error" in result:
        raise HTTPException(status_code=400, detail=result["error"])

    return result


@app.post("/api/export/structure")
async def export_structure(request: ExportRequest):
    """Export structure in different formats"""
    try:
        atom_array = pdb_to_atom_array(request.pdb_content)

        if request.format == "cif":
            content = atom_array_to_cif(atom_array)
            filename = "structure.cif"
        else:
            content = request.pdb_content
            filename = "structure.pdb"

        result = {
            "filename": filename,
            "content": content,
            "format": request.format,
        }

        if request.include_confidences and request.confidences:
            result["confidences_json"] = json.dumps(request.confidences, indent=2)

        return result

    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/stored/{job_id}")
async def get_stored_structure(job_id: str):
    """Get stored structure from a completed job for cross-panel data flow"""
    if job_id not in stored_structures:
        raise HTTPException(status_code=404, detail="Structure not found")

    return stored_structures[job_id]


@app.get("/api/jobs/{job_id}", response_model=JobStatus)
async def get_job_status(job_id: str):
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
    if job_id in jobs:
        del jobs[job_id]
    if job_id in stored_structures:
        del stored_structures[job_id]
    return {"status": "deleted"}


# ============== Main ==============

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
