"""
E2E Test: NL Pipeline — PQQ + Ca Protein Design

Replicates the EXACT frontend natural-language pipeline flow against
the local Docker API (localhost:8000). Mirrors natural-language.ts step
by step, including:
  1. ai_parse         — NL intent parsing (→ metal, ligand, contig)
  2. scaffold_search   — RCSB PDB scaffold discovery (→ motif PDB)
  3. metal_binding_design (single mode) — Backbone generation
  4. scout_filter      — Pre-filter metal backbones (CPU only)
  5. mpnn              — Sequence design with HSAB metal config
  6. rf3               — Ligand-aware structure validation
  7. analyze_design    — Coordination analysis + filter evaluation
  8. save_design_history — Persist to design history (non-fatal)
  9. check_lessons     — Detect failure patterns (non-fatal)

Response format note:
  Docker/RunPod returns: { id, status: "COMPLETED", output: { status: "completed", result: {...} } }
  This script unwraps the RunPod envelope to match what the frontend sees after proxy.

Usage:
  cd backend/serverless
  python test_nl_pqq_ca.py
"""
import builtins
import json
import math
import sys
import time
import requests
from datetime import datetime
from pathlib import Path

# Force unbuffered output
_orig_print = builtins.print
def print(*args, **kwargs):
    kwargs.setdefault('flush', True)
    _orig_print(*args, **kwargs)

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
BASE_URL = "http://localhost:8000"
API_URL = f"{BASE_URL}/runsync"

PQQ_SMILES = "OC(=O)c1[nH]c2c(c1)C(=O)C(=O)c3nc(cc(C(O)=O)c23)C(O)=O"

NUM_BACKBONES = 5
NUM_SEQS_PER_BACKBONE = 4
DEFAULT_CONTIG = "70-100"  # frontend fallback for metal+ligand with no AI chain lengths

OUTPUT_DIR = Path(__file__).parent / "output_nl_pqq_ca"

TIMEOUTS = {
    "ai_parse": 60,
    "scaffold_search": 120,
    "metal_binding_design": 1800,
    "scout_filter": 600,
    "mpnn": 600,
    "rf3": 600,
    "analyze_design": 120,
    "save_design_history": 30,
    "check_lessons": 30,
}


# ---------------------------------------------------------------------------
# HTTP helper — unwraps RunPod response envelope
# ---------------------------------------------------------------------------
def call_api(task: str, params: dict, timeout: int | None = None) -> dict:
    """Call the local Docker API and unwrap the RunPod response envelope.

    Docker returns: { id, status: "COMPLETED", output: { status: "completed", result: {...} } }
    This function returns the inner output: { status: "completed", result: {...} }
    matching what the frontend sees after proxy unwrapping.
    """
    if timeout is None:
        timeout = TIMEOUTS.get(task, 600)

    payload = {"input": {"task": task, **params}}

    for attempt in range(1, 4):
        try:
            resp = requests.post(
                API_URL,
                json=payload,
                headers={"Content-Type": "application/json"},
                timeout=timeout,
            )
            raw = resp.json()

            # Unwrap RunPod envelope: { id, status, output: { status, result, error } }
            outer_status = str(raw.get("status", "")).upper()
            if outer_status == "FAILED" and "output" not in raw:
                # RunPod-level failure (no inner output)
                return {"status": "failed", "error": raw.get("error", "RunPod FAILED (no output)")}
            if "output" in raw:
                inner = raw["output"]
                # Normalize: RunPod outer status is "COMPLETED"/"FAILED",
                # inner status is "completed"/"failed"/"error"
                return inner

            # Direct handler response (no RunPod wrapper)
            return raw

        except requests.exceptions.ConnectionError:
            if attempt < 3:
                print(f"  [RETRY] Connection refused (attempt {attempt}/3), waiting 5s...")
                time.sleep(5)
            else:
                return {"status": "failed", "error": "Connection refused after 3 attempts. Is Docker running?"}
        except requests.exceptions.Timeout:
            return {"status": "failed", "error": f"Timeout after {timeout}s"}
        except Exception as e:
            return {"status": "failed", "error": str(e)}


def check_health() -> bool:
    """Quick health check."""
    result = call_api("health", {}, timeout=10)
    return result.get("status") != "failed"


# ---------------------------------------------------------------------------
# Parsing helpers — match frontend's parseMpnnSequences / extractConfidences
# ---------------------------------------------------------------------------
def parse_sequences(mpnn_result: dict) -> list[str]:
    """Extract sequences from MPNN result (matches frontend parseMpnnSequences).

    MPNN returns { sequences: [{ content: ">header\\nSEQ\\n>header2\\nSEQ2", ... }] }
    """
    sequences = []
    for seq_entry in mpnn_result.get("sequences", []):
        content = seq_entry.get("content", "") or seq_entry.get("sequence", "")
        if ">" in content:
            # Multi-sequence FASTA
            current_seq = ""
            for line in content.split("\n"):
                line = line.strip()
                if line.startswith(">"):
                    if current_seq:
                        sequences.append(current_seq)
                    current_seq = ""
                elif line:
                    current_seq += line
            if current_seq:
                sequences.append(current_seq)
        elif content.strip():
            sequences.append(content.strip())
    return sequences


def extract_rf3_pdb_content(result: dict) -> str:
    """Extract PDB content from RF3 result (matches frontend extractRf3PdbContent).

    RF3 returns { predictions: [{ content: "ATOM..." }] } or { pdb_content: "..." }
    """
    predictions = result.get("predictions", [])
    if predictions and isinstance(predictions, list):
        pred = predictions[0]
        pdb = pred.get("content") or pred.get("pdb_content", "")
        if pdb:
            return pdb

    # Fallback: direct fields
    return result.get("pdb_content", "") or result.get("content", "")


def extract_rf3_confidences(result: dict) -> dict:
    """Extract confidence metrics from RF3 result (matches frontend extractConfidences).

    RF3 returns confidences in varying formats:
    - { confidences: { summary_confidences: { ptm, overall_plddt, overall_pae } } }
    - { predictions: [{ confidences: { summary_confidences: { ... } } }] }
    - { confidences: { overall_plddt, chain_ptm: [...], overall_pae } }
    """
    # Try per-prediction confidences first
    predictions = result.get("predictions", [])
    pred_conf = predictions[0].get("confidences", {}) if predictions else {}

    # Then top-level confidences
    top_conf = result.get("confidences", {})
    confidences = pred_conf or top_conf

    if not confidences:
        return {"ptm": None, "plddt": None, "pae": None}

    # Try summary_confidences (nested format), then fall back to direct fields
    summary = confidences.get("summary_confidences", {})
    src = summary if summary else confidences

    # pTM: try ptm, then chain_ptm[0]
    ptm = src.get("ptm")
    if ptm is None:
        chain_ptm = confidences.get("chain_ptm") or src.get("chain_ptm")
        if isinstance(chain_ptm, list) and chain_ptm:
            ptm = chain_ptm[0]

    plddt = src.get("overall_plddt") or src.get("mean_plddt")
    pae = src.get("overall_pae")

    return {"ptm": ptm, "plddt": plddt, "pae": pae}


def find_metal_coordinating_residues(pdb_content: str, cutoff: float = 3.5) -> list[str]:
    """Find residues coordinating a metal ion (matches frontend findMetalCoordinatingResidues).

    Returns list of fixed position strings like ["A1", "A5", "A10"].
    """
    metal_atoms = []
    lines = pdb_content.split("\n")

    for line in lines:
        if not line.startswith("HETATM"):
            continue
        res_name = line[17:20].strip()
        # Metal residues: 1-2 uppercase letters (TB, MG, CA, ZN, FE, etc.)
        import re
        if not re.match(r'^[A-Z]{1,2}$', res_name):
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            metal_atoms.append((x, y, z))
        except (ValueError, IndexError):
            pass

    if not metal_atoms:
        return []

    coord_residues = set()
    cutoff_sq = cutoff * cutoff

    for line in lines:
        if not line.startswith("ATOM"):
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except (ValueError, IndexError):
            continue

        for mx, my, mz in metal_atoms:
            dx, dy, dz = x - mx, y - my, z - mz
            if dx * dx + dy * dy + dz * dz <= cutoff_sq:
                try:
                    res_num = int(line[22:26].strip())
                    coord_residues.add(res_num)
                except ValueError:
                    pass
                break

    return [f"A{r}" for r in sorted(coord_residues)]


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------
def save_outputs(results: dict, output_dir: Path):
    """Save all pipeline outputs to disk."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Parse intent
    with open(output_dir / "parse_intent.json", "w") as f:
        json.dump(results.get("parse_intent", {}), f, indent=2)

    # Scaffold search
    with open(output_dir / "scaffold_search.json", "w") as f:
        # Strip large PDB content from scaffold search for json
        scaffold_data = results.get("scaffold_search", {})
        scaffold_slim = {k: v for k, v in scaffold_data.items()
                        if k not in ("candidate_pdbs", "best_pdb_content")}
        json.dump(scaffold_slim, f, indent=2)

    # Backbone PDBs
    for i, pdb in enumerate(results.get("backbones", [])):
        (output_dir / f"backbone_{i:03d}.pdb").write_text(pdb)

    # Sequences FASTA
    seqs = results.get("sequences", [])
    if seqs:
        with open(output_dir / "sequences.fasta", "w") as f:
            for i, seq in enumerate(seqs):
                f.write(f">design_{i+1:03d}\n{seq}\n")

    # RF3 PDBs
    for i, entry in enumerate(results.get("rf3_results", [])):
        pdb = entry.get("pdb_content", "")
        if pdb:
            (output_dir / f"rf3_seq_{i+1:03d}.pdb").write_text(pdb)

    # Summary JSON
    with open(output_dir / "summary.json", "w") as f:
        summary = _build_summary(results)
        json.dump(summary, f, indent=2)

    # Human-readable report
    report = format_report(results)
    (output_dir / "report.txt").write_text(report)


def _build_summary(results: dict) -> dict:
    """Build a structured summary without bulky PDB strings."""
    rf3_summaries = []
    for entry in results.get("rf3_results", []):
        rf3_summaries.append({
            "backbone_idx": entry.get("backbone_idx"),
            "seq_idx": entry.get("seq_idx"),
            "ptm": entry.get("ptm"),
            "plddt": entry.get("plddt"),
            "pae": entry.get("pae"),
        })

    analysis_summaries = []
    for entry in results.get("analysis_results", []):
        analysis_summaries.append({
            "backbone_idx": entry.get("backbone_idx"),
            "seq_idx": entry.get("seq_idx"),
            "coordination_number": entry.get("coordination_number"),
            "ligand_coordination": entry.get("ligand_coordination"),
            "protein_coordination": entry.get("protein_coordination"),
            "protein_donors": entry.get("protein_donors"),
            "ligand_donors": entry.get("ligand_donors"),
            "alanine_pct": entry.get("alanine_pct"),
            "aromatic_pct": entry.get("aromatic_pct"),
            "total_residues": entry.get("total_residues"),
            "filter_passed": entry.get("filter_passed"),
            "failed_filters": entry.get("failed_filters"),
        })

    # Scaffold search summary (without PDB blobs)
    scaffold_raw = results.get("scaffold_search", {})
    scaffold_summary = {}
    if scaffold_raw.get("status") != "failed":
        scaffold_summary = {
            "recommended_action": scaffold_raw.get("recommended_action"),
            "num_pdb_hits": scaffold_raw.get("num_pdb_hits"),
            "num_validated": scaffold_raw.get("num_validated"),
            "best_candidate": scaffold_raw.get("best_candidate"),
            "num_candidates": len(scaffold_raw.get("candidates", [])),
        }

    return {
        "timestamp": datetime.now().isoformat(),
        "query": "Designing a protein to bind PQQ and Ca",
        "parse_intent": results.get("parse_intent", {}),
        "scaffold_search": scaffold_summary,
        "num_backbones": len(results.get("backbones", [])),
        "num_sequences": len(results.get("sequences", [])),
        "contig_used": results.get("contig_used"),
        "rf3_results": rf3_summaries,
        "analysis_results": analysis_summaries,
        "scout_filter": results.get("scout_filter", {}),
        "save_history": results.get("save_history", {}),
        "check_lessons": results.get("check_lessons", {}),
        "timings": results.get("timings", {}),
    }


def format_report(results: dict) -> str:
    """Generate a human-readable report."""
    lines = []
    lines.append("=" * 70)
    lines.append("NL Pipeline E2E Report: PQQ + Ca Protein Design")
    lines.append(f"Timestamp: {datetime.now().isoformat()}")
    lines.append("=" * 70)

    # Intent
    intent = results.get("parse_intent", {}).get("result", results.get("parse_intent", {}))
    lines.append("\n--- Step 1: Parse Intent ---")
    lines.append(f"  Metal: {intent.get('metal_type', 'N/A')}")
    lines.append(f"  Ligand: {intent.get('ligand_name', 'N/A')}")
    lines.append(f"  Chain length: {intent.get('chain_length_min', '?')}-{intent.get('chain_length_max', '?')}")
    lines.append(f"  Design type: metal (mapped from intent)")
    lines.append(f"  Contig used: {results.get('contig_used', 'N/A')}")

    # Scaffold search
    scaffold = results.get("scaffold_search", {})
    candidates = scaffold.get("candidates", [])
    lines.append(f"\n--- Step 2: Scaffold Search ---")
    lines.append(f"  Recommended action: {scaffold.get('recommended_action', 'N/A')}")
    lines.append(f"  Candidates found: {len(candidates)}")
    for c in candidates[:5]:
        lines.append(f"    {c.get('pdb_id', '?')} - score: {c.get('total_score', '?')}")

    # Backbones
    backbones = results.get("backbones", [])
    lines.append(f"\n--- Step 3: Backbone Generation (metal_binding_design single) ---")
    lines.append(f"  Generated: {len(backbones)} backbones")

    # Scout filter
    scout = results.get("scout_filter", {})
    if scout:
        lines.append(f"\n--- Step 4: Scout Pre-filter (CPU only, metal mode) ---")
        lines.append(f"  Pre-filter passed: {scout.get('pre_filter_passed', 'N/A')}/{scout.get('original_count', 'N/A')}")

    # Sequences
    seqs = results.get("sequences", [])
    lines.append(f"\n--- Step 5: MPNN Sequence Design (HSAB hard-acid config for CA) ---")
    lines.append(f"  Total sequences: {len(seqs)}")
    lines.append(f"  bias_AA: A:-2.0, omit_AA: C")
    lines.append(f"  pack_side_chains: true, pack_with_ligand_context: true")

    # RF3
    rf3_list = results.get("rf3_results", [])
    lines.append(f"\n--- Step 6: RF3 Validation (ligand-aware) ---")
    lines.append(f"  Validated: {len(rf3_list)} sequences")
    lines.append(f"  ligand_smiles: {PQQ_SMILES[:50]}...")

    ptms = [r["ptm"] for r in rf3_list if r.get("ptm") is not None]
    plddts = [r["plddt"] for r in rf3_list if r.get("plddt") is not None]
    if ptms:
        lines.append(f"  pTM:   min={min(ptms):.3f}  max={max(ptms):.3f}  avg={sum(ptms)/len(ptms):.3f}")
    if plddts:
        lines.append(f"  pLDDT: min={min(plddts):.1f}  max={max(plddts):.1f}  avg={sum(plddts)/len(plddts):.1f}")

    # Analysis
    analysis_list = results.get("analysis_results", [])
    lines.append(f"\n--- Step 7: Coordination Analysis (UnifiedDesignAnalyzer) ---")
    lines.append(f"  Analyzed: {len(analysis_list)} designs")

    passed = [a for a in analysis_list if a.get("filter_passed")]
    lines.append(f"  Passed filters: {len(passed)}/{len(analysis_list)}")

    for a in analysis_list:
        tag = "PASS" if a.get("filter_passed") else "FAIL"
        lines.append(
            f"  [{tag}] bb{a.get('backbone_idx', '?')}/seq{a.get('seq_idx', '?')}: "
            f"coord={a.get('coordination_number', '?')} "
            f"(prot={a.get('protein_coordination', '?')}, lig={a.get('ligand_coordination', '?')}) "
            f"pTM={a.get('ptm', '?')} pLDDT={a.get('plddt', '?')}"
        )
        if a.get("failed_filters"):
            lines.append(f"         failed: {a['failed_filters']}")

    # History & lessons
    history = results.get("save_history", {})
    lessons = results.get("check_lessons", {})
    lines.append(f"\n--- Step 8-9: History & Lessons (automatic, non-fatal) ---")
    lines.append(f"  History: {history.get('status', 'N/A')} run_id={history.get('run_id', 'N/A')}")
    lines.append(f"  Lessons: trigger={lessons.get('trigger_detected', 'N/A')}")

    # Timings
    timings = results.get("timings", {})
    lines.append(f"\n--- Timings ---")
    for step, t in timings.items():
        lines.append(f"  {step}: {t:.1f}s")
    total = sum(timings.values())
    lines.append(f"  TOTAL: {total:.1f}s")

    lines.append("\n" + "=" * 70)
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main pipeline — mirrors frontend natural-language.ts step by step
# ---------------------------------------------------------------------------
def run_pipeline() -> dict:
    """Execute the full NL pipeline, mirroring the exact frontend flow."""

    results = {
        "parse_intent": {},
        "scaffold_search": {},
        "backbones": [],
        "sequences": [],
        "rf3_results": [],
        "analysis_results": [],
        "scout_filter": {},
        "save_history": {},
        "check_lessons": {},
        "contig_used": None,
        "timings": {},
    }

    # ------------------------------------------------------------------
    # Health check
    # ------------------------------------------------------------------
    print("Checking Docker API health...")
    if not check_health():
        print("[FAIL] Docker API not reachable at", API_URL)
        print("       Start Docker with: docker-wsl skill or docker compose up")
        sys.exit(1)
    print("[OK] API is healthy\n")

    # ==================================================================
    # Step 1: Parse Intent (frontend: parse_intent step)
    # Frontend calls: api.parseIntent(prompt) → POST runsync { task: "ai_parse", query }
    # ==================================================================
    print("=" * 60)
    print("STEP 1: Parse Intent (ai_parse)")
    print("=" * 60)

    t0 = time.time()
    parse_result = call_api("ai_parse", {
        "query": "Designing a protein to bind PQQ and Ca",
    })
    results["timings"]["parse_intent"] = time.time() - t0

    if parse_result.get("status") == "failed" or parse_result.get("status") == "error":
        print(f"[FAIL] ai_parse: {parse_result.get('error')}")
        print("       Continuing with defaults (keyword fallback)...")
        # Frontend keywordFallbackParse: detects "ca" → metal=CA, "pqq" → ligand_name
        intent = {}
        contig = DEFAULT_CONTIG
        metal = "CA"
        ligand_name = "pqq"
        ligand_smiles = ""
        results["parse_intent"] = {"status": "fallback", "result": {
            "metal_type": "CA", "ligand_name": "pqq", "design_type": "metal"
        }}
    else:
        intent = parse_result.get("result", parse_result)
        results["parse_intent"] = intent
        metal = intent.get("metal_type", "CA")
        ligand_name = intent.get("ligand_name", "pqq")
        ligand_smiles = intent.get("ligand_smiles", "")
        chain_min = intent.get("chain_length_min")
        chain_max = intent.get("chain_length_max")

        # Frontend configure step logic:
        # if AI returns chain_length_min/max → use them
        # else if metal+ligand → "70-100"
        if chain_min and chain_max:
            contig = f"{chain_min}-{chain_max}"
        else:
            contig = DEFAULT_CONTIG

        print(f"  Metal: {metal}")
        print(f"  Ligand: {ligand_name}")
        print(f"  Contig: {contig}")
        print(f"  Parser type: {intent.get('parser_type', 'N/A')}")
        print(f"  Confidence: {intent.get('confidence', 'N/A')}")
        print(f"  Elapsed: {results['timings']['parse_intent']:.1f}s")

    results["contig_used"] = contig

    # Frontend mapDesignType: if metal_type → "metal"
    design_type = "metal"
    print()

    # ==================================================================
    # Step 2: Resolve Structure (frontend: resolve_structure step)
    # For metal+ligand with no PDB input, frontend returns:
    #   "No input PDB — backend will generate {metal} binding template"
    # This is a pass-through step for our case.
    # ==================================================================
    print("=" * 60)
    print("STEP 2: Resolve Structure (skipped — no input PDB)")
    print("=" * 60)
    print("  Metal-ligand template will be generated by backend")
    print()

    # ==================================================================
    # Step 3: Scaffold Search (frontend: scaffold_search_nl step)
    # Frontend calls: api.searchScaffold({ metal, ligand_name, fetch_pdb, limit })
    # ==================================================================
    print("=" * 60)
    print("STEP 3: Scaffold Search")
    print("=" * 60)

    t0 = time.time()
    scaffold_raw = call_api("scaffold_search", {
        "metal": metal,
        "ligand_name": ligand_name,
        "fetch_pdb": True,
        "limit": 10,
    })
    results["timings"]["scaffold_search"] = time.time() - t0

    motif_pdb = ""
    scaffold_pdb_id = None

    if scaffold_raw.get("status") in ("failed", "error"):
        print(f"[WARN] scaffold_search failed: {scaffold_raw.get('error')}")
        print("       Proceeding without scaffold (auto-template mode)")
        results["scaffold_search"] = {"status": "failed", "error": scaffold_raw.get("error")}
    else:
        sr = scaffold_raw.get("result", scaffold_raw)
        results["scaffold_search"] = sr
        candidates = sr.get("candidates", [])
        recommended = sr.get("recommended_action", "de_novo")
        print(f"  Recommended action: {recommended}")
        print(f"  Found {len(candidates)} scaffold candidates")
        for c in candidates[:5]:
            print(f"    {c.get('pdb_id', '?')}: score={c.get('total_score', '?')}, "
                  f"coord={c.get('coordination_number', '?')}")

        # Frontend extracts PDB from candidate_pdbs dict (matches searchScaffold response)
        candidate_pdbs = sr.get("candidate_pdbs", {})
        if recommended == "scaffold" and candidates:
            best = candidates[0]
            scaffold_pdb_id = best.get("pdb_id")
            # Try candidate_pdbs[pdb_id] first, then best_pdb_content fallback
            motif_pdb = candidate_pdbs.get(scaffold_pdb_id, "") or sr.get("best_pdb_content", "")
            if motif_pdb:
                print(f"  Using scaffold: {scaffold_pdb_id}")
            else:
                print(f"  Scaffold {scaffold_pdb_id} recommended but no PDB content fetched")

    print(f"  Elapsed: {results['timings']['scaffold_search']:.1f}s")
    print()

    # ==================================================================
    # Step 4: Configure (frontend: configure step)
    # Merges intent + scaffold + user params → design config
    # This is done inline here, not a separate API call.
    # ==================================================================
    # Frontend configure step determines:
    # - contig: AI chain lengths or metal+ligand default "70-100"
    # - bury_ligand: from AI parser (default true)
    # - ligand code: maps name → 3-letter code (pqq → PQQ)
    ligand_code_map = {
        "citrate": "CIT", "pqq": "PQQ", "atp": "ATP", "nad": "NAD", "heme": "HEM",
    }
    ligand_code = ligand_code_map.get(ligand_name.lower(), ligand_name.upper())
    bury_ligand = intent.get("bury_ligand", True) if intent else True

    # ==================================================================
    # Step 5: Backbone Generation (frontend: rfd3_nl step — metal single mode)
    # Frontend calls: api.submitMetalSingleDesign({...})
    # which sends: POST runsync { task: "metal_binding_design", mode: "single", ... }
    # ==================================================================
    print("=" * 60)
    print(f"STEP 5: Backbone Generation — metal_binding_design single ({NUM_BACKBONES} backbones)")
    print("=" * 60)

    backbone_params = {
        "mode": "single",
        "metal": metal,
        "ligand": ligand_code,
        "contig": contig,
        "cfg_scale": 2.0,
        "num_designs": NUM_BACKBONES,
        "num_timesteps": 200,
        "step_scale": 1.5,   # frontend default
        "bury_ligand": bury_ligand,
        "seed": 42,          # frontend hardcodes seed=42
    }
    if motif_pdb:
        backbone_params["motif_pdb"] = motif_pdb

    print(f"  Contig: {contig}")
    print(f"  Metal: {metal}, Ligand: {ligand_code}")
    print(f"  cfg_scale: 2.0, step_scale: 1.5, seed: 42")
    print(f"  bury_ligand: {bury_ligand}")
    print(f"  Motif PDB: {'from scaffold ' + (scaffold_pdb_id or '') if motif_pdb else 'auto-template'}")

    t0 = time.time()
    bb_result = call_api("metal_binding_design", backbone_params)
    results["timings"]["backbone_generation"] = time.time() - t0

    # Debug: print raw response keys to diagnose extraction issues
    print(f"  [DEBUG] bb_result keys: {list(bb_result.keys())}")
    print(f"  [DEBUG] bb_result status: {bb_result.get('status')}")
    if "result" in bb_result:
        inner_debug = bb_result["result"]
        if isinstance(inner_debug, dict):
            print(f"  [DEBUG] result keys: {list(inner_debug.keys())}")
            if "designs" in inner_debug:
                print(f"  [DEBUG] designs count: {len(inner_debug['designs'])}")
            if "backbones" in inner_debug:
                print(f"  [DEBUG] backbones count: {len(inner_debug['backbones'])}")
        else:
            print(f"  [DEBUG] result type: {type(inner_debug).__name__}")

    if bb_result.get("status") in ("failed", "error"):
        print(f"[FAIL] metal_binding_design: {bb_result.get('error')}")
        print("       Cannot continue without backbones.")
        save_outputs(results, OUTPUT_DIR)
        return results

    # Frontend: singleResult.result.designs[].content
    inner = bb_result.get("result", bb_result)
    designs = inner.get("designs", [])
    if not designs:
        designs = inner.get("backbones", [])

    backbones = []
    for d in designs:
        pdb = d if isinstance(d, str) else d.get("content", d.get("pdb_content", ""))
        if pdb:
            backbones.append(pdb)

    results["backbones"] = backbones
    print(f"  Generated: {len(backbones)} backbones")
    if inner.get("errors"):
        print(f"  Errors: {inner['errors']}")
    print(f"  Elapsed: {results['timings']['backbone_generation']:.1f}s")

    if not backbones:
        print("[FAIL] No backbones produced.")
        save_outputs(results, OUTPUT_DIR)
        return results

    print()

    # ==================================================================
    # Step 5.5: Scout Filter (frontend: scout_filter_nl step — metal pre-filter)
    # Frontend: metal_single_mode → pre_filter_only=true (CPU checks, no GPU)
    # Skipped if ≤1 backbone. Does NOT re-emit pdbOutputs.
    # ==================================================================
    print("=" * 60)
    print(f"STEP 5.5: Scout Pre-filter ({len(backbones)} backbones)")
    print("=" * 60)

    if len(backbones) <= 1:
        print("  Skipped: single backbone — pre-filter not needed")
        results["scout_filter"] = {"skipped": True, "reason": "Only 1 backbone"}
    else:
        t0 = time.time()
        scout_resp = call_api("scout_filter", {
            "backbone_pdbs": backbones,
            "ptm_threshold": 0.6,
            "plddt_threshold": 0.65,
            "target_metal": metal,
            "ligand_name": ligand_name,
            "pre_filter_only": True,  # CPU only for metal single mode
        })
        results["timings"]["scout_filter"] = time.time() - t0

        if scout_resp.get("status") in ("failed", "error"):
            print(f"[WARN] Scout filter failed (non-fatal): {scout_resp.get('error')}")
            results["scout_filter"] = {"status": "failed", "error": scout_resp.get("error")}
        else:
            sr = scout_resp.get("result", scout_resp)
            results["scout_filter"] = sr
            pf_passed = sr.get("pre_filter_passed", len(backbones))
            pf_failed = sr.get("pre_filter_failed", 0)
            print(f"  Pre-filter: {pf_passed}/{len(backbones)} passed, {pf_failed} eliminated")
            # Note: frontend does NOT re-emit pdbOutputs here.
            # MPNN finds backbones from rfd3_nl step directly.
            # We keep all backbones for MPNN (pre-filter is advisory in single mode).
        print(f"  Elapsed: {results['timings'].get('scout_filter', 0):.1f}s")

    print()

    # ==================================================================
    # Step 6: MPNN Sequence Design (frontend: mpnn_nl step)
    # Frontend uses getMetalMpnnConfig('CA') → HSAB hard-acid config:
    #   bias_AA: "A:-2.0", omit_AA: "C"
    #   pack_side_chains: true, pack_with_ligand_context: true,
    #   use_side_chain_context: true
    # Plus auto-detected fixed_positions for metal-coordinating residues.
    # ==================================================================
    print("=" * 60)
    print(f"STEP 6: MPNN Sequence Design ({NUM_SEQS_PER_BACKBONE} seqs/backbone)")
    print("=" * 60)
    print("  HSAB config for CA (hard acid):")
    print("    bias_AA=A:-2.0, omit_AA=C")
    print("    pack_side_chains=true, pack_with_ligand_context=true")
    print("    use_side_chain_context=true")

    all_sequences = []
    seq_backbone_map = []  # Track which backbone each seq came from

    t0 = time.time()
    for bb_idx, pdb in enumerate(backbones):
        print(f"  Backbone {bb_idx + 1}/{len(backbones)}...", end=" ")

        # Frontend: findMetalCoordinatingResidues() → fixed_positions
        fixed_positions = find_metal_coordinating_residues(pdb)
        if fixed_positions:
            print(f"fixing {len(fixed_positions)} coord residues...", end=" ")

        mpnn_params = {
            "pdb_content": pdb,
            "num_sequences": NUM_SEQS_PER_BACKBONE,
            "temperature": 0.1,
            # HSAB hard-acid config (matches getMetalMpnnConfig('CA'))
            "bias_AA": "A:-2.0",
            "omit_AA": "C",
            "pack_side_chains": True,
            "pack_with_ligand_context": True,
            "use_side_chain_context": True,
        }
        if fixed_positions:
            mpnn_params["fixed_positions"] = fixed_positions

        mpnn_result = call_api("mpnn", mpnn_params)

        if mpnn_result.get("status") in ("failed", "error"):
            print(f"FAIL: {mpnn_result.get('error', 'unknown')}")
            continue

        # Frontend: parseMpnnSequences handles FASTA content in sequences[].content
        inner = mpnn_result.get("result", mpnn_result)
        seqs = parse_sequences(inner)
        print(f"got {len(seqs)} sequences")

        for seq in seqs:
            all_sequences.append(seq)
            seq_backbone_map.append(bb_idx)

    results["timings"]["mpnn"] = time.time() - t0
    results["sequences"] = all_sequences
    print(f"  Total sequences: {len(all_sequences)}")
    print(f"  Elapsed: {results['timings']['mpnn']:.1f}s")

    if not all_sequences:
        print("[FAIL] No sequences produced.")
        save_outputs(results, OUTPUT_DIR)
        return results

    print()

    # ==================================================================
    # Step 7: RF3 Validation (frontend: rf3_nl step)
    # Frontend calls: api.submitRF3Prediction({ sequence, ligand_smiles, metal })
    # Extracts PDB from predictions[0].content, confidences from
    # predictions[0].confidences.summary_confidences or chain_ptm[0]
    # ==================================================================
    print("=" * 60)
    print(f"STEP 7: RF3 Validation ({len(all_sequences)} sequences)")
    print("=" * 60)

    # Frontend gets ligand_smiles from intent result data
    # For PQQ, the LIGAND_SMILES lookup table in natural-language.ts doesn't have PQQ,
    # so ligand_smiles comes from AI parser or is empty.
    # We use PQQ_SMILES directly since we know the ligand.
    rf3_ligand_smiles = ligand_smiles if ligand_smiles else PQQ_SMILES

    rf3_results = []

    t0 = time.time()
    for seq_idx, seq in enumerate(all_sequences):
        bb_idx = seq_backbone_map[seq_idx]
        print(f"  Seq {seq_idx + 1}/{len(all_sequences)} (bb{bb_idx})...", end=" ")

        rf3_resp = call_api("rf3", {
            "sequence": seq,
            "name": f"Design {bb_idx + 1} - Seq {(seq_idx % NUM_SEQS_PER_BACKBONE) + 1}",
            "ligand_smiles": rf3_ligand_smiles,
            "metal": metal,
        })

        if rf3_resp.get("status") in ("failed", "error"):
            print(f"FAIL: {rf3_resp.get('error', 'unknown')}")
            rf3_results.append({
                "backbone_idx": bb_idx,
                "seq_idx": seq_idx,
                "ptm": None,
                "plddt": None,
                "pae": None,
                "pdb_content": "",
                "error": rf3_resp.get("error"),
            })
            continue

        # Unwrap: frontend extractRf3PdbContent + extractConfidences
        inner = rf3_resp.get("result", rf3_resp)
        pdb_content = extract_rf3_pdb_content(inner)
        confidences = extract_rf3_confidences(inner)

        entry = {
            "backbone_idx": bb_idx,
            "seq_idx": seq_idx,
            "ptm": confidences.get("ptm"),
            "plddt": confidences.get("plddt"),
            "pae": confidences.get("pae"),
            "pdb_content": pdb_content,
        }
        rf3_results.append(entry)

        ptm = entry["ptm"]
        plddt = entry["plddt"]
        ptm_str = f"{ptm:.3f}" if ptm is not None else "N/A"
        plddt_str = f"{plddt:.1f}" if plddt is not None else "N/A"
        has_pdb = "PDB" if pdb_content else "no-PDB"
        print(f"pTM={ptm_str}  pLDDT={plddt_str}  [{has_pdb}]")

    results["timings"]["rf3"] = time.time() - t0
    results["rf3_results"] = rf3_results

    valid_count = sum(1 for r in rf3_results if r.get("pdb_content"))
    print(f"  Validated: {valid_count}/{len(all_sequences)} (got PDB)")
    print(f"  Elapsed: {results['timings']['rf3']:.1f}s")
    print()

    # ==================================================================
    # Step 8: Analysis (frontend: analysis step)
    # Frontend calls: api.analyzeDesign({ pdb_content, metal_type, ligand_name, design_type })
    # Only analyzes RF3 outputs (pdbOutputs with id starting "rf3-")
    # ==================================================================
    print("=" * 60)
    print("STEP 8: Coordination Analysis (analyze_design)")
    print("=" * 60)

    analysis_results = []

    t0 = time.time()
    for rf3_entry in rf3_results:
        pdb_content = rf3_entry.get("pdb_content", "")
        if not pdb_content:
            continue

        bb_idx = rf3_entry["backbone_idx"]
        seq_idx = rf3_entry["seq_idx"]
        print(f"  Analyzing bb{bb_idx}/seq{seq_idx}...", end=" ")

        analysis_resp = call_api("analyze_design", {
            "pdb_content": pdb_content,
            "metal_type": metal,
            "ligand_name": ligand_name,
            "design_type": design_type,
            "design_params": {},
        })

        entry = {
            "backbone_idx": bb_idx,
            "seq_idx": seq_idx,
            "ptm": rf3_entry.get("ptm"),
            "plddt": rf3_entry.get("plddt"),
            "pae": rf3_entry.get("pae"),
        }

        if analysis_resp.get("status") in ("failed", "error"):
            print(f"FAIL: {analysis_resp.get('error', 'unknown')}")
            entry["error"] = analysis_resp.get("error")
        else:
            a = analysis_resp.get("result", analysis_resp)
            # Frontend reads: a.metrics.*, a.filter_passed, a.failed_filters
            metrics = a.get("metrics", a)
            entry.update({
                "coordination_number": metrics.get("coordination_number"),
                "ligand_coordination": metrics.get("ligand_coordination"),
                "protein_coordination": metrics.get("protein_coordination"),
                "protein_donors": metrics.get("protein_donors"),
                "ligand_donors": metrics.get("ligand_donors"),
                "alanine_pct": metrics.get("alanine_pct"),
                "aromatic_pct": metrics.get("aromatic_pct"),
                "total_residues": metrics.get("total_residues"),
                "filter_passed": a.get("filter_passed"),
                "filter_preset": a.get("filter_preset"),
                "failed_filters": [
                    f"{f['metric']}={f.get('value', '?')}"
                    for f in a.get("failed_filters", [])
                ],
            })
            coord = entry.get("coordination_number", "?")
            passed = "PASS" if entry.get("filter_passed") else "FAIL"
            print(f"[{passed}] coord={coord} preset={entry.get('filter_preset', '?')}")

        analysis_results.append(entry)

    results["timings"]["analysis"] = time.time() - t0
    results["analysis_results"] = analysis_results

    passed_count = sum(1 for a in analysis_results if a.get("filter_passed"))
    print(f"  Passed: {passed_count}/{len(analysis_results)}")
    print(f"  Elapsed: {results['timings']['analysis']:.1f}s")
    print()

    # ==================================================================
    # Step 9: Save Design History (frontend: save_history_nl step)
    # Automatic, non-fatal. Frontend calls: api.saveDesignHistory({...})
    # ==================================================================
    print("=" * 60)
    print("STEP 9: Save Design History (non-fatal)")
    print("=" * 60)

    t0 = time.time()
    try:
        history_resp = call_api("save_design_history", {
            "session_name": f"nl_{design_type}",
            "design_params": {
                "design_type": design_type,
                "target_metal": metal,
                "ligand_name": ligand_name,
                "user_prompt": "Designing a protein to bind PQQ and Ca",
                "contig": contig,
            },
            "design_outputs": {
                "backbone_count": len(backbones),
                "sequence_count": len(all_sequences),
            },
            "design_metrics": {
                "design_type": design_type,
                "total_analyzed": len(analysis_results),
                "top_plddt": max((r.get("plddt") or 0 for r in rf3_results), default=0),
            },
        })
        results["timings"]["save_history"] = time.time() - t0

        if history_resp.get("status") in ("failed", "error"):
            print(f"  [WARN] Save failed (non-fatal): {history_resp.get('error')}")
            results["save_history"] = {"status": "failed", "error": history_resp.get("error")}
        else:
            inner = history_resp.get("result", history_resp)
            run_id = inner.get("run_id", "unknown")
            session_id = inner.get("session_id", "unknown")
            print(f"  Saved: run_id={run_id}, session_id={session_id}")
            results["save_history"] = {"status": "completed", "run_id": run_id, "session_id": session_id}
    except Exception as e:
        print(f"  [WARN] Save failed (non-fatal): {e}")
        results["save_history"] = {"status": "failed", "error": str(e)}
        results["timings"]["save_history"] = time.time() - t0

    print()

    # ==================================================================
    # Step 10: Check Lessons (frontend: check_lessons_nl step)
    # Automatic, non-fatal. Frontend calls: api.checkLessons({ result: metricsResult })
    # ==================================================================
    print("=" * 60)
    print("STEP 10: Check Lessons (non-fatal)")
    print("=" * 60)

    t0 = time.time()
    try:
        top_plddt = max((r.get("plddt") or 0 for r in rf3_results), default=0)
        lessons_resp = call_api("check_lessons", {
            "result": {
                "design_type": design_type,
                "total_analyzed": len(analysis_results),
                "top_plddt": top_plddt,
                "outcome": "failure" if top_plddt < 0.5 else None,
            },
        })
        results["timings"]["check_lessons"] = time.time() - t0

        if lessons_resp.get("status") in ("failed", "error"):
            print(f"  [WARN] Check failed (non-fatal): {lessons_resp.get('error')}")
            results["check_lessons"] = {"status": "failed", "error": lessons_resp.get("error")}
        else:
            inner = lessons_resp.get("result", lessons_resp)
            trigger = inner.get("trigger_detected", False)
            history_count = inner.get("history_count", 0)
            print(f"  Trigger detected: {trigger}")
            print(f"  History count: {history_count}")
            if trigger and inner.get("trigger"):
                t = inner["trigger"]
                print(f"  Type: {t.get('type')}")
                print(f"  Description: {t.get('description')}")
            results["check_lessons"] = {
                "status": "completed",
                "trigger_detected": trigger,
                "history_count": history_count,
                "trigger": inner.get("trigger"),
            }
    except Exception as e:
        print(f"  [WARN] Check failed (non-fatal): {e}")
        results["check_lessons"] = {"status": "failed", "error": str(e)}
        results["timings"]["check_lessons"] = time.time() - t0

    return results


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("NL Pipeline E2E Test: PQQ + Ca Protein Design")
    print(f"Started: {datetime.now().isoformat()}")
    print(f"API: {API_URL}")
    print(f"Output: {OUTPUT_DIR}")
    print()
    print("This test mirrors the EXACT frontend flow (natural-language.ts)")
    print("Steps: parse -> scaffold_search -> configure -> rfd3(metal single)")
    print("       -> scout(pre-filter) -> mpnn(HSAB) -> rf3(ligand) -> analysis")
    print("       -> save_history -> check_lessons")
    print("=" * 70)
    print()

    results = run_pipeline()

    # Save everything
    print("\n" + "=" * 60)
    print("Saving outputs...")
    save_outputs(results, OUTPUT_DIR)
    print(f"  Saved to: {OUTPUT_DIR}")

    # Print report
    report = format_report(results)
    print("\n" + report)

    print(f"\nFinished: {datetime.now().isoformat()}")


if __name__ == "__main__":
    main()
