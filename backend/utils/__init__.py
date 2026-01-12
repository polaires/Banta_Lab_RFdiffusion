"""
Backend Utilities for AI-Guided Protein Engineering

This package contains analysis tools for:
- PDB fetching from RCSB
- SASA (Solvent Accessible Surface Area) calculation
- Metal coordination geometry analysis
- Secondary structure assignment (DSSP)
- ESM-2 sequence analysis
"""

from .pdb_fetch import fetch_pdb, fetch_cif, validate_pdb_id
from .coordination import (
    analyze_coordination_geometry,
    get_coordinating_atoms,
    classify_donor_type,
    calculate_coordination_number,
)
from .sasa import calculate_sasa, classify_burial
from .dssp import assign_secondary_structure

__all__ = [
    # PDB fetching
    "fetch_pdb",
    "fetch_cif",
    "validate_pdb_id",
    # Coordination analysis
    "analyze_coordination_geometry",
    "get_coordinating_atoms",
    "classify_donor_type",
    "calculate_coordination_number",
    # SASA
    "calculate_sasa",
    "classify_burial",
    # Secondary structure
    "assign_secondary_structure",
]
