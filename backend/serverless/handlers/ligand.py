"""Ligand resolution and feature analysis handlers.

Tasks: resolve_ligand, analyze_ligand_features
"""

import traceback
from typing import Dict, Any


def handle_resolve_ligand(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Resolve a ligand name to SMILES, residue code, and chemistry info.

    Uses the LigandResolver priority chain:
    1. Pre-built templates (curated metal-ligand complexes)
    2. PDB experimental structures
    3. PubChem SMILES + isomeric SMILES table
    4. Calculated fallback

    Input:
        ligand_name: str - Ligand name (e.g., "citrate", "azobenzene", "caffeine")
        metal_type: str - Optional metal symbol (e.g., "TB", "CA")
        isomer_spec: str - Optional isomer ("cis", "trans")

    Returns:
        {
            "status": "completed",
            "result": {
                "success": true/false,
                "smiles": "...",
                "residue_code": "LIG",
                "source": "template|pdb|pubchem|calculated",
                "ligand_fixing_strategy": "fix_all|fix_coordination|diffuse_all",
                "warnings": [...]
            }
        }
    """
    ligand_name = job_input.get("ligand_name", "")
    if not ligand_name:
        return {"status": "failed", "error": "Must provide 'ligand_name' parameter"}

    metal_type = job_input.get("metal_type")
    isomer_spec = job_input.get("isomer_spec")

    try:
        from ligand_resolver import LigandResolver
        resolver = LigandResolver()
        resolved = resolver.resolve(
            ligand_name=ligand_name,
            metal_type=metal_type,
            isomer_spec=isomer_spec,
        )

        result = {
            "success": resolved.resolved,
            "smiles": resolved.smiles or "",
            "residue_code": resolved.residue_code,
            "source": resolved.source,
            "name": resolved.name,
            "warnings": resolved.warnings,
        }

        # Include coordination info if available
        if resolved.coordination:
            result["ligand_fixing_strategy"] = resolved.coordination.get(
                "recommended_fixing_strategy", "fix_all"
            )
            result["coordination_donors"] = resolved.coordination.get(
                "ligand_metal_donors", []
            )

        return {"status": "completed", "result": result}
    except Exception as e:
        print(f"[Resolve Ligand] Error: {e}")
        traceback.print_exc()
        return {"status": "failed", "error": f"Ligand resolution error: {str(e)}"}


def handle_analyze_ligand_features(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Analyze ligand features using three-layer system:
    1. Knowledge base (cached PDB crystal structure data)
    2. ChemicalFeatures (RDKit pharmacophore perception)
    3. Geometry filter (SMARTS-based fallback)

    Input:
        ligand_name: str - Ligand name (e.g., "citrate", "pqq")
        smiles: str - Optional SMILES string
        metal_type: str - Optional metal symbol (e.g., "TB", "CA")
        pdb_content: str - Optional PDB content for 3D analysis
        record_evidence: dict - Optional PDB evidence to cache

    Returns:
        {
            "status": "completed",
            "result": {
                "ligand_name": "...",
                "smiles": "...",
                "source": "knowledge_base|chemicalfeatures|geometry_filter",
                "features": [...],
                "coordination_donors": [...],
                ...
            }
        }
    """
    ligand_name = job_input.get("ligand_name", "")
    if not ligand_name:
        return {"status": "failed", "error": "Must provide 'ligand_name' parameter"}

    smiles = job_input.get("smiles")
    metal_type = job_input.get("metal_type")
    pdb_content = job_input.get("pdb_content")
    record_evidence = job_input.get("record_evidence")

    try:
        from ligand_features import analyze_ligand_features
        result = analyze_ligand_features(
            ligand_name=ligand_name,
            smiles=smiles,
            metal=metal_type,
            pdb_content=pdb_content,
            record_evidence_data=record_evidence,
        )
        return {"status": "completed", "result": result}
    except Exception as e:
        print(f"[Analyze Ligand Features] Error: {e}")
        traceback.print_exc()
        return {"status": "failed", "error": f"Ligand feature analysis error: {str(e)}"}
