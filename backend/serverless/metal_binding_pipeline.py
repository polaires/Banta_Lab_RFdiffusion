"""
Metal Binding Site Design Pipeline - Round 7b Style

Implements iterative parameter optimization for MONOMER scaffold design
around metal-ligand binding sites.

Key difference from interface_metal_ligand_design:
- MONOMER scaffold (single chain) vs dimer (two chains)
- Parameter sweep capability (3 sizes × 3 CFG scales)
- Quality tier filtering (S/A/B/C/F)
- Production mode with best config
"""

from dataclasses import dataclass, field
from typing import Dict, Any, List, Optional, Tuple
from enum import Enum
import time
import uuid


class QualityTier(str, Enum):
    """Design quality tiers based on structural metrics."""
    S = "S"  # Superior: CN≥8, RMSD<0.8, SASA<2, pLDDT≥0.85
    A = "A"  # Excellent: CN≥8, RMSD<1.2, SASA<5, pLDDT≥0.80
    B = "B"  # Good: CN≥7, RMSD<1.5, pLDDT≥0.75
    C = "C"  # Acceptable: CN≥6, RMSD<2.0, pLDDT≥0.70
    F = "F"  # Failed: Below thresholds


# Quality thresholds based on BindCraft + ln_citrate.json analysis
QUALITY_THRESHOLDS = {
    "S": {
        "coordination_number": 8,
        "geometry_rmsd": 0.8,      # Å - strict geometric match
        "metal_sasa": 2,           # Å² - well buried metal
        "ligand_contacts": 3,      # contacts with ligand
        "plddt": 0.85,             # from BindCraft stringent
        "ptm": 0.80,
    },
    "A": {
        "coordination_number": 8,
        "geometry_rmsd": 1.2,
        "metal_sasa": 5,
        "ligand_contacts": 3,
        "plddt": 0.80,             # from BindCraft standard
        "ptm": 0.75,
    },
    "B": {
        "coordination_number": 7,
        "geometry_rmsd": 1.5,
        "metal_sasa": 10,
        "ligand_contacts": 2,
        "plddt": 0.75,
        "ptm": 0.70,
    },
    "C": {
        "coordination_number": 6,
        "geometry_rmsd": 2.0,
        "metal_sasa": 20,
        "ligand_contacts": 1,
        "plddt": 0.70,             # minimum viable
        "ptm": 0.65,
    },
}

# Proceed to MPNN criteria
PROCEED_CRITERIA = {
    "min_tier": "B",
    "min_pass_rate": 0.30,
    "min_tier_a_count": 3,
}

# Default contig ranges for MONOMER design (larger than dimers)
MONOMER_CONTIG_RANGES = {
    "small": "100-120",
    "medium": "110-130",
    "large": "130-150",
}

# CFG scale values for sweep
CFG_SCALE_VALUES = [2.5, 3.0, 3.5]


@dataclass
class MetalBindingConfig:
    """Single configuration for parameter sweep or single-run."""
    name: str
    contig_range: str      # "100-120", "110-130", "130-150"
    cfg_scale: float       # 1.5, 2.0, 2.5
    metal: str
    ligand: Optional[str] = None
    num_designs: int = 10
    # Optional overrides from configure step (educated defaults)
    num_timesteps: Optional[int] = None   # Default 200 in build_rfd3_params
    step_scale: Optional[float] = None    # Default 1.5 in build_rfd3_params
    gamma_0: Optional[float] = None       # RFD3 gamma_0 param
    ligand_fix_atoms: Optional[str] = None  # Partial ligand fixing: "O2,O5,O7" or None for "all"

    def to_dict(self) -> Dict[str, Any]:
        d = {
            "name": self.name,
            "contig_range": self.contig_range,
            "cfg_scale": self.cfg_scale,
            "metal": self.metal,
            "ligand": self.ligand,
            "num_designs": self.num_designs,
        }
        if self.num_timesteps is not None:
            d["num_timesteps"] = self.num_timesteps
        if self.step_scale is not None:
            d["step_scale"] = self.step_scale
        if self.gamma_0 is not None:
            d["gamma_0"] = self.gamma_0
        return d


@dataclass
class DesignResult:
    """Result from a single design run."""
    name: str
    config_name: str
    pdb_content: str
    sequence: str
    plddt: float
    ptm: float
    pae: Optional[float] = None
    coordination_number: Optional[int] = None
    geometry_rmsd: Optional[float] = None
    metal_sasa: Optional[float] = None
    tier: str = "F"
    status: str = "pending"  # "pass", "review", "fail"
    seed: Optional[int] = None

    def to_dict(self, include_pdb: bool = True) -> Dict[str, Any]:
        d = {
            "name": self.name,
            "config_name": self.config_name,
            "sequence": self.sequence,
            "plddt": self.plddt,
            "ptm": self.ptm,
            "pae": self.pae,
            "coordination_number": self.coordination_number,
            "geometry_rmsd": self.geometry_rmsd,
            "metal_sasa": self.metal_sasa,
            "tier": self.tier,
            "status": self.status,
            "seed": self.seed,
        }
        if include_pdb and self.pdb_content:
            d["pdb_content"] = self.pdb_content
        return d


@dataclass
class SweepProgress:
    """Progress tracking for parameter sweep."""
    session_id: str
    mode: str  # "sweep" or "production"
    status: str  # "running", "completed", "cancelled", "failed"
    current_config: int
    total_configs: int
    current_design: int
    designs_per_config: int
    total_generated: int
    total_passing: int
    total_review: int
    total_failed: int
    pass_rate: float
    best_design: Optional[Dict[str, Any]] = None
    results: List[DesignResult] = field(default_factory=list)
    config_rankings: List[Dict[str, Any]] = field(default_factory=list)
    error: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "session_id": self.session_id,
            "mode": self.mode,
            "status": self.status,
            "current_config": self.current_config,
            "total_configs": self.total_configs,
            "current_design": self.current_design,
            "designs_per_config": self.designs_per_config,
            "total_generated": self.total_generated,
            "total_passing": self.total_passing,
            "total_review": self.total_review,
            "total_failed": self.total_failed,
            "pass_rate": self.pass_rate,
            "best_design": self.best_design,
            "config_rankings": self.config_rankings,
            "error": self.error,
        }


def generate_sweep_configs(
    metal: str,
    ligand: Optional[str] = None,
    sizes: Optional[List[str]] = None,
    cfg_scales: Optional[List[float]] = None,
    num_designs: int = 10,
) -> List[MetalBindingConfig]:
    """
    Generate sweep configurations (3 sizes × 3 CFG scales = 9 configs).

    Args:
        metal: Metal code (TB, EU, CA, ZN, etc.)
        ligand: Optional ligand code (CIT, PQQ, etc.)
        sizes: List of sizes to test ["small", "medium", "large"]
        cfg_scales: List of CFG scales to test [1.5, 2.0, 2.5]
        num_designs: Designs per config

    Returns:
        List of MetalBindingConfig objects
    """
    sizes = sizes or ["small", "medium", "large"]
    cfg_scales = cfg_scales or CFG_SCALE_VALUES

    configs = []
    for size_name in sizes:
        contig = MONOMER_CONTIG_RANGES.get(size_name, "110-130")
        for cfg in cfg_scales:
            cfg_label = "low" if cfg == 1.5 else "mid" if cfg == 2.0 else "high"
            configs.append(MetalBindingConfig(
                name=f"{size_name}_cfg{cfg_label}",
                contig_range=contig,
                cfg_scale=cfg,
                metal=metal,
                ligand=ligand,
                num_designs=num_designs,
            ))
    return configs


def assign_quality_tier(metrics: Dict[str, Any]) -> str:
    """
    Assign S/A/B/C/F quality tier based on design metrics.

    When coordination geometry metrics (coordination_number, geometry_rmsd) are
    available, they are used alongside pLDDT/pTM for full evaluation.
    When only pLDDT/pTM are available (e.g., sweep mode before metal validation),
    tier is assigned based on confidence metrics alone.

    Args:
        metrics: Dict with plddt, ptm, and optionally coordination_number, geometry_rmsd, etc.

    Returns:
        Tier string: "S", "A", "B", "C", or "F"
    """
    plddt = metrics.get("plddt", 0)
    ptm = metrics.get("ptm", 0)
    cn = metrics.get("coordination_number")
    rmsd = metrics.get("geometry_rmsd")

    has_metal_metrics = cn is not None and rmsd is not None

    for tier in ["S", "A", "B", "C"]:
        thresh = QUALITY_THRESHOLDS[tier]
        # Check confidence metrics (always required)
        if plddt < thresh["plddt"] or ptm < thresh.get("ptm", 0):
            continue
        # Check metal coordination metrics only if available
        if has_metal_metrics:
            if cn < thresh["coordination_number"] or rmsd > thresh["geometry_rmsd"]:
                continue
        return tier
    return "F"


def determine_status(tier: str, plddt: float, ptm: float) -> str:
    """
    Determine pass/review/fail status from tier and metrics.

    Args:
        tier: Quality tier (S/A/B/C/F)
        plddt: Per-residue pLDDT score
        ptm: pTM score

    Returns:
        Status string: "pass", "review", or "fail"
    """
    if tier in ["S", "A"]:
        return "pass"
    elif tier == "B":
        return "review"
    else:
        return "fail"


def rank_configs_by_performance(
    results: List[DesignResult],
    configs: List[MetalBindingConfig],
) -> List[Dict[str, Any]]:
    """
    Rank sweep configs by pass rate and average quality metrics.

    Args:
        results: List of design results
        configs: List of configs used

    Returns:
        Sorted list of config rankings
    """
    config_stats = {}

    for config in configs:
        config_results = [r for r in results if r.config_name == config.name]
        if not config_results:
            continue

        passing = [r for r in config_results if r.status == "pass"]
        review = [r for r in config_results if r.status == "review"]

        avg_plddt = sum(r.plddt for r in config_results) / len(config_results)
        avg_ptm = sum(r.ptm for r in config_results) / len(config_results)
        pass_rate = len(passing) / len(config_results)

        # Tier distribution
        tier_counts = {"S": 0, "A": 0, "B": 0, "C": 0, "F": 0}
        for r in config_results:
            tier_counts[r.tier] = tier_counts.get(r.tier, 0) + 1

        config_stats[config.name] = {
            "config_name": config.name,
            "contig_range": config.contig_range,
            "cfg_scale": config.cfg_scale,
            "total_designs": len(config_results),
            "passing": len(passing),
            "review": len(review),
            "failed": len(config_results) - len(passing) - len(review),
            "pass_rate": pass_rate,
            "avg_plddt": avg_plddt,
            "avg_ptm": avg_ptm,
            "tier_distribution": tier_counts,
        }

    # Sort by: pass_rate (desc), then avg_plddt (desc)
    rankings = sorted(
        config_stats.values(),
        key=lambda x: (x["pass_rate"], x["avg_plddt"]),
        reverse=True,
    )

    for i, r in enumerate(rankings):
        r["rank"] = i + 1

    return rankings


def build_rfd3_params(
    config: MetalBindingConfig,
    motif_pdb: str,
    seed: int,
    fixed_atoms: Optional[Dict[str, str]] = None,
    buried_atoms: Optional[Dict[str, str]] = None,
    hbond_acceptors: Optional[Dict[str, str]] = None,
) -> Dict[str, Any]:
    """
    Build RFD3 parameters for MONOMER metal binding design.

    Key difference from dimer design:
    - Single contig range (e.g., "110-130") not "60-80,/0,60-80"
    - No chain separator in contig

    Note: RFD3 automatically assigns chain IDs to non-polymer residues during
    parsing. We cannot use hardcoded chain names like "X1" or "L1" because
    the actual chain IDs depend on the input structure. Instead, we rely on
    the `ligand` parameter to identify metal/ligand residues.

    Args:
        config: Sweep configuration
        motif_pdb: PDB content with metal-ligand complex
        seed: Random seed for reproducibility
        fixed_atoms: Dict of residue/ligand to fixed atoms (optional, requires actual chain IDs)
        buried_atoms: Dict for RASA buried conditioning (optional, requires actual chain IDs)
        hbond_acceptors: Dict for H-bond acceptor conditioning (optional, requires actual chain IDs)

    Returns:
        RFD3 parameters dict
    """
    # Build ligand string - tells RFD3 which residue names are non-protein
    ligand_str = config.metal
    if config.ligand:
        ligand_str = f"{config.metal},{config.ligand}"

    # Default fixed atoms: fix metal and ligand completely
    # Use RESIDUE NAMES (e.g., "TB", "CIT") not chain+resnum ("X1", "L1")
    # Foundry SDK resolves components by residue name from the atom array
    #
    # Partial ligand fixing: config.ligand_fix_atoms can specify which atoms
    # to fix (e.g., "O2,O5,O7" for coordinating oxygens only). Non-listed
    # atoms co-diffuse with the protein, allowing optimal H-bond contacts.
    if fixed_atoms is None:
        fixed_atoms = {
            config.metal.upper(): "all",  # Metal center: always fully fixed
        }
        if config.ligand:
            ligand_fix = getattr(config, "ligand_fix_atoms", None) or "all"
            fixed_atoms[config.ligand.upper()] = ligand_fix

    # Default buried atoms: bury metal AND ligand for pocket formation
    # Both must be buried so RFD3 wraps protein around the entire complex
    if buried_atoms is None:
        buried_atoms = {config.metal.upper(): "all"}
        if config.ligand:
            buried_atoms[config.ligand.upper()] = "all"

    # Default H-bond acceptors — extracted dynamically from ligand atoms in PDB
    # IMPORTANT: Exclude coordinating atoms (those in ligand_fix_atoms) from H-bond
    # acceptors. Their lone pairs are occupied by metal coordination bonds, so
    # conditioning H-bonds on them conflicts with the metal binding geometry.
    if hbond_acceptors is None and config.ligand:
        ligand_fix = getattr(config, "ligand_fix_atoms", None)
        coordinating_atoms = set()
        if ligand_fix and ligand_fix != "all":
            coordinating_atoms = {a.strip() for a in ligand_fix.split(",")}
        hbond_acceptors = get_ligand_hbond_atoms(
            config.ligand, motif_pdb=motif_pdb,
            exclude_atoms=coordinating_atoms,
        )

    params = {
        "pdb_content": motif_pdb,
        "contig": config.contig_range,  # MONOMER: single range like "110-130"
        "ligand": ligand_str,
        "select_fixed_atoms": fixed_atoms,
        "select_buried": buried_atoms,
        "use_classifier_free_guidance": True,
        "cfg_scale": config.cfg_scale,
        "num_designs": 1,  # One at a time for tracking
        "seed": seed,
        "num_timesteps": config.num_timesteps or 200,
        "step_scale": config.step_scale or 1.5,
        "infer_ori_strategy": "com",  # Center on fixed atoms
    }

    if hbond_acceptors:
        params["select_hbond_acceptor"] = hbond_acceptors

    # Optional gamma_0 from configure step
    if config.gamma_0 is not None:
        params["gamma_0"] = config.gamma_0

    # Scaffold mode: when fixed_atoms includes protein residues (keys like "A1", "A2"),
    # add unindex so RFD3 knows these are spatially fixed but sequence-flexible
    # Match "A" followed by digits only (not ligand codes like "ATP", "ADP")
    import re
    protein_fixed = [k for k in fixed_atoms if re.match(r'^A\d+$', k)]
    if protein_fixed:
        params["unindex"] = ",".join(protein_fixed)

    return params


def get_ligand_hbond_atoms(
    ligand: str,
    motif_pdb: Optional[str] = None,
    exclude_atoms: Optional[set] = None,
) -> Optional[Dict[str, str]]:
    """
    Get H-bond acceptor atoms for a ligand by extracting O/N atoms from the PDB.

    Dynamically identifies H-bond acceptor atoms (oxygen, nitrogen) from the ligand's
    HETATM records in the motif PDB. Falls back to RCSB CCD lookup if no PDB provided.

    Args:
        ligand: Ligand code (CIT, PQQ, etc.)
        motif_pdb: PDB content containing the ligand (optional)
        exclude_atoms: Set of atom names to exclude (e.g. coordinating atoms
            whose lone pairs are occupied by metal bonds)

    Returns:
        Dict for select_hbond_acceptor parameter, or None
    """
    ligand_upper = ligand.upper()
    exclude = exclude_atoms or set()

    # Extract O and N atoms from ligand HETATM records in the PDB
    if motif_pdb:
        acceptor_atoms = []
        for line in motif_pdb.split("\n"):
            if not line.startswith("HETATM"):
                continue
            resname = line[17:20].strip()
            if resname != ligand_upper:
                continue
            atom_name = line[12:16].strip()
            element = line[76:78].strip() if len(line) >= 78 else ""
            # H-bond acceptors are O and N atoms
            if element in ("O", "N") or atom_name.startswith("O") or atom_name.startswith("N"):
                if atom_name not in exclude:
                    acceptor_atoms.append(atom_name)
                else:
                    print(f"[HBond] Excluding {atom_name} from acceptors "
                          f"(metal-coordinating)")

        if acceptor_atoms:
            atoms_str = ",".join(sorted(set(acceptor_atoms)))
            return {ligand_upper: atoms_str}

    # Fallback: try to get atoms from RCSB CCD
    try:
        import requests
        resp = requests.get(
            f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_upper}",
            timeout=5,
        )
        if resp.status_code == 200:
            data = resp.json()
            # Extract atom names that are O or N from the ideal coordinates
            atoms = data.get("rcsb_chem_comp_descriptor", {})
            # CCD doesn't directly list atoms, but we can check the formula
            formula = data.get("chem_comp", {}).get("formula", "")
            if "O" in formula or "N" in formula:
                # Has potential H-bond acceptors, but we can't get specific atom names
                # from CCD descriptor alone — return None to skip conditioning
                print(f"[HBond] {ligand_upper} has O/N atoms per CCD formula ({formula}) "
                      f"but no PDB to extract atom names — skipping H-bond conditioning")
    except Exception:
        pass

    return None


def build_mpnn_params(
    pdb_content: str,
    metal: str,
    ligand: Optional[str] = None,
    num_seqs: int = 4,
    temperature: float = 0.1,
    hsab_bias: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Build LigandMPNN parameters for sequence design.

    Uses educated defaults from ln_citrate experiments:
    - temperature=0.1 (conservative, preserves coordination geometry)
    - No bias for hard Lewis acids (atom context handles it naturally)
    - omit_AA='C' for hard metals (Cys incompatible with hard Lewis acids)
    - pack_side_chains=True (optimizes chi angles around metal)
    - use_side_chain_context=True (uses fixed residue sidechains as context)

    Args:
        pdb_content: Backbone PDB from RFD3
        metal: Metal code for context
        ligand: Optional ligand code
        num_seqs: Number of sequences to generate
        temperature: Sampling temperature (0.1 recommended for metal sites)
        hsab_bias: HSAB-compliant amino acid bias string (None for hard acids)

    Returns:
        MPNN parameters dict compatible with run_mpnn_inference()
    """
    hsab_class = classify_hsab(metal)

    params: Dict[str, Any] = {
        "pdb_content": pdb_content,
        "num_sequences": num_seqs,
        "model_type": "ligand_mpnn",
        "temperature": temperature,
        "pack_side_chains": True,
        "use_side_chain_context": True,
    }

    # Omit cysteine for hard Lewis acids (thiolate incompatible)
    if hsab_class == "hard":
        params["omit_AA"] = "C"

    # Only apply AA bias for soft/borderline metals
    # Hard acids: atom context naturally places carboxylates at coordination sites
    if hsab_bias:
        params["bias_AA"] = hsab_bias

    return params


# Session storage for tracking active pipelines
_active_sessions: Dict[str, SweepProgress] = {}


def create_session(mode: str, total_configs: int, designs_per_config: int) -> str:
    """Create a new pipeline session."""
    session_id = str(uuid.uuid4())[:8]
    _active_sessions[session_id] = SweepProgress(
        session_id=session_id,
        mode=mode,
        status="running",
        current_config=0,
        total_configs=total_configs,
        current_design=0,
        designs_per_config=designs_per_config,
        total_generated=0,
        total_passing=0,
        total_review=0,
        total_failed=0,
        pass_rate=0.0,
    )
    return session_id


def get_session(session_id: str) -> Optional[SweepProgress]:
    """Get session by ID."""
    return _active_sessions.get(session_id)


def update_session(session_id: str, updates: Dict[str, Any]) -> None:
    """Update session progress."""
    session = _active_sessions.get(session_id)
    if session:
        for key, value in updates.items():
            if hasattr(session, key):
                setattr(session, key, value)


def cancel_session(session_id: str) -> bool:
    """Cancel a running session."""
    session = _active_sessions.get(session_id)
    if session and session.status == "running":
        session.status = "cancelled"
        return True
    return False


def delete_session(session_id: str) -> bool:
    """Delete a session."""
    if session_id in _active_sessions:
        del _active_sessions[session_id]
        return True
    return False


# HSAB classification for metals
HSAB_CLASSIFICATION = {
    # Hard Lewis acids (lanthanides, alkaline earth)
    "TB": "hard", "EU": "hard", "GD": "hard", "CA": "hard", "MG": "hard",
    "LA": "hard", "CE": "hard", "SM": "hard", "YB": "hard", "DY": "hard",
    # Borderline Lewis acids (transition metals)
    "ZN": "borderline", "FE": "borderline", "MN": "borderline",
    "CO": "borderline", "NI": "borderline",
    # Soft Lewis acids
    "CU": "soft", "AG": "soft",
}

# HSAB amino acid biases for common metals
# Hard Lewis acids: NO bias — LigandMPNN atom context handles coordination
# naturally (ln_citrate R6 lesson: bias-free + atom_context → 0.81A RMSD,
# vs E:3.0,D:3.0 bias → 61% charged → RMSD 21-25A non-foldable)
HSAB_BIASES = {
    # Hard Lewis acids: no bias needed (atom context is sufficient)
    "TB": None, "EU": None, "GD": None, "CA": None, "MG": None,
    "LA": None, "CE": None, "SM": None, "YB": None, "DY": None,
    # Borderline Lewis acids (Zn2+, Fe2+): moderate bias
    "ZN": "H:3.0,C:2.0,E:2.0,D:2.0,N:1.5,Q:1.5",
    "FE": "H:2.5,C:2.0,E:2.5,D:2.5,N:1.5",
    # Soft Lewis acids (Cu+): strong S/N donor bias
    "CU": "C:3.0,M:2.5,H:2.0,E:1.0,D:1.0",
}

# Cache for resolved ligand SMILES (avoids repeated HTTP lookups)
_ligand_smiles_cache: Dict[str, Optional[str]] = {}


def _fetch_ccd_smiles(ligand_code: str) -> Optional[str]:
    """Fetch SMILES from RCSB Chemical Component Dictionary (CCD).

    Uses the RCSB REST API to look up canonical SMILES for any PDB ligand code.
    Falls back gracefully on network errors.
    """
    try:
        import requests
        url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_code.upper()}"
        resp = requests.get(url, timeout=10)
        if resp.status_code == 200:
            data = resp.json()
            # CCD stores SMILES in rcsb_chem_comp_descriptor
            descriptors = data.get("rcsb_chem_comp_descriptor", [])
            for desc in descriptors:
                if desc.get("type") == "SMILES_CANONICAL":
                    return desc.get("descriptor")
                if desc.get("type") == "SMILES":
                    return desc.get("descriptor")
        print(f"[MetalBinding] CCD lookup for {ligand_code}: HTTP {resp.status_code}")
    except Exception as e:
        print(f"[MetalBinding] CCD lookup failed for {ligand_code}: {e}")
    return None


def classify_hsab(metal: str) -> str:
    """Classify a metal by HSAB theory."""
    return HSAB_CLASSIFICATION.get(metal.upper(), "borderline")


def get_hsab_bias(metal: str) -> Optional[str]:
    """Get HSAB-compliant amino acid bias for a metal.
    Returns None for hard Lewis acids (atom context is sufficient)."""
    return HSAB_BIASES.get(metal.upper())


def get_ligand_smiles(ligand_code: str) -> Optional[str]:
    """Get SMILES string for a ligand code (for RF3 ligand-aware prediction).

    Resolution order:
    1. In-memory cache (from previous lookups)
    2. LigandResolver (templates + PDB structures)
    3. RCSB CCD HTTP lookup (canonical SMILES for any PDB ligand)
    """
    if not ligand_code:
        return None

    code = ligand_code.upper()

    # Check cache first
    if code in _ligand_smiles_cache:
        return _ligand_smiles_cache[code]

    # Try the existing LigandResolver (templates, PDB, PubChem)
    try:
        from ligand_resolver import resolve_ligand
        resolved = resolve_ligand(code)
        if resolved.resolved and resolved.smiles:
            _ligand_smiles_cache[code] = resolved.smiles
            print(f"[MetalBinding] Resolved {code} SMILES via {resolved.source}")
            return resolved.smiles
    except Exception as e:
        print(f"[MetalBinding] LigandResolver failed for {code}: {e}")

    # Fall back to RCSB CCD database
    smiles = _fetch_ccd_smiles(code)
    if smiles:
        _ligand_smiles_cache[code] = smiles
        print(f"[MetalBinding] Resolved {code} SMILES via CCD")
        return smiles

    _ligand_smiles_cache[code] = None
    print(f"[MetalBinding] Could not resolve SMILES for {code}")
    return None
