# ligand_features.py
"""
Self-Growing Ligand Knowledge Base + ChemicalFeatures Fallback

Three-layer ligand feature identification:
1. Knowledge Base: Self-growing JSON store caching coordination data from PDB crystal structures
2. ChemicalFeatures: RDKit pharmacophore perception for any ligand with SMILES
3. Geometry Filter: Existing filter_donors_by_geometry() as final fallback

Follows DesignHistoryManager storage pattern for persistence.
"""

import json
import logging
import os
from datetime import datetime
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

try:
    from rdkit import RDConfig
    from rdkit.Chem import ChemicalFeatures
    HAS_CHEMICAL_FEATURES = True
except ImportError:
    HAS_CHEMICAL_FEATURES = False


# =============================================================================
# SEED DATA — Bootstrap knowledge base with well-characterized ligands
# =============================================================================

SEED_LIGANDS: Dict[str, Dict[str, Any]] = {
    "pqq": {
        "ligand_id": "pqq",
        "smiles": "OC(=O)c1cc2c(O)c3[nH]c(C(O)=O)c(=O)c3c(C(O)=O)c2[nH]1",
        "aliases": ["pyrroloquinoline quinone", "methoxatin"],
        "pdb_codes": ["PQQ"],
        "heavy_atom_count": 22,
        "hsab_character": "hard",
        "pharmacophore_features": [
            {"atom_idx": 1, "atom_name": "O1", "element": "O", "type": "acceptor"},
            {"atom_idx": 2, "atom_name": "O2", "element": "O", "type": "donor"},
            {"atom_idx": 8, "atom_name": "O5", "element": "O", "type": "donor"},
            {"atom_idx": 11, "atom_name": "N1", "element": "N", "type": "donor"},
            {"atom_idx": 14, "atom_name": "O7", "element": "O", "type": "acceptor"},
            {"atom_idx": 15, "atom_name": "O8", "element": "O", "type": "donor"},
            {"atom_idx": 18, "atom_name": "O9", "element": "O", "type": "acceptor"},
            {"atom_idx": 19, "atom_name": "O10", "element": "O", "type": "donor"},
            {"atom_idx": 21, "atom_name": "N2", "element": "N", "type": "donor"},
        ],
        "coordination_summary": {
            "CA": {"donors": ["O5", "N1", "O7"], "denticity": 3, "evidence_count": 2},
        },
        "notes": "Redox cofactor, binds Ca2+ in quinoprotein dehydrogenases",
    },
    "citrate": {
        "ligand_id": "citrate",
        "smiles": "OC(=O)CC(O)(CC(O)=O)C(O)=O",
        "aliases": ["cit", "citric acid", "FLC"],
        "pdb_codes": ["CIT", "FLC"],
        "heavy_atom_count": 13,
        "hsab_character": "hard",
        "pharmacophore_features": [
            {"atom_idx": 1, "atom_name": "O1", "element": "O", "type": "acceptor"},
            {"atom_idx": 2, "atom_name": "O2", "element": "O", "type": "donor"},
            {"atom_idx": 5, "atom_name": "O3", "element": "O", "type": "donor"},
            {"atom_idx": 8, "atom_name": "O5", "element": "O", "type": "acceptor"},
            {"atom_idx": 9, "atom_name": "O6", "element": "O", "type": "donor"},
            {"atom_idx": 11, "atom_name": "O7", "element": "O", "type": "acceptor"},
            {"atom_idx": 12, "atom_name": "O8", "element": "O", "type": "donor"},
        ],
        "coordination_summary": {
            "TB": {"donors": ["O2", "O5", "O7"], "denticity": 3, "evidence_count": 1},
            "CA": {"donors": ["O2", "O5"], "denticity": 2, "evidence_count": 1},
            "LA": {"donors": ["O2", "O5", "O7"], "denticity": 3, "evidence_count": 1},
        },
        "notes": "Tricarboxylate chelator, common lanthanide binder",
    },
    "atp": {
        "ligand_id": "atp",
        "smiles": "c1nc(c2c(n1)n(cn2)[C@@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N",
        "aliases": ["adenosine triphosphate"],
        "pdb_codes": ["ATP"],
        "heavy_atom_count": 31,
        "hsab_character": "hard",
        "pharmacophore_features": [
            {"atom_idx": 0, "atom_name": "N1", "element": "N", "type": "acceptor"},
            {"atom_idx": 6, "atom_name": "N6", "element": "N", "type": "donor"},
        ],
        "coordination_summary": {
            "MG": {"donors": ["O2B", "O2G"], "denticity": 2, "evidence_count": 3},
        },
        "notes": "Triphosphate chelates Mg2+ through beta/gamma phosphate oxygens",
    },
    "edta": {
        "ligand_id": "edta",
        "smiles": "OC(=O)CN(CCN(CC(O)=O)CC(O)=O)CC(O)=O",
        "aliases": ["ethylenediaminetetraacetic acid"],
        "pdb_codes": ["EDT"],
        "heavy_atom_count": 18,
        "hsab_character": "hard",
        "pharmacophore_features": [
            {"atom_idx": 1, "atom_name": "O1", "element": "O", "type": "acceptor"},
            {"atom_idx": 2, "atom_name": "O2", "element": "O", "type": "donor"},
            {"atom_idx": 4, "atom_name": "N1", "element": "N", "type": "donor"},
            {"atom_idx": 7, "atom_name": "N2", "element": "N", "type": "donor"},
            {"atom_idx": 10, "atom_name": "O5", "element": "O", "type": "acceptor"},
            {"atom_idx": 11, "atom_name": "O6", "element": "O", "type": "donor"},
            {"atom_idx": 14, "atom_name": "O7", "element": "O", "type": "acceptor"},
            {"atom_idx": 15, "atom_name": "O8", "element": "O", "type": "donor"},
        ],
        "coordination_summary": {
            "CA": {"donors": ["O2", "N1", "N2", "O6", "O8"], "denticity": 5, "evidence_count": 1},
            "FE": {"donors": ["O2", "N1", "N2", "O6", "O8", "O1"], "denticity": 6, "evidence_count": 1},
        },
        "notes": "Hexadentate chelator, binds most transition metals and lanthanides",
    },
    "heme": {
        "ligand_id": "heme",
        "smiles": "",
        "aliases": ["hem", "protoporphyrin IX", "heme b"],
        "pdb_codes": ["HEM", "HEB"],
        "heavy_atom_count": 43,
        "hsab_character": "borderline",
        "pharmacophore_features": [
            {"atom_idx": 0, "atom_name": "NA", "element": "N", "type": "donor"},
            {"atom_idx": 1, "atom_name": "NB", "element": "N", "type": "donor"},
            {"atom_idx": 2, "atom_name": "NC", "element": "N", "type": "donor"},
            {"atom_idx": 3, "atom_name": "ND", "element": "N", "type": "donor"},
        ],
        "coordination_summary": {
            "FE": {"donors": ["NA", "NB", "NC", "ND"], "denticity": 4, "evidence_count": 100},
        },
        "notes": "Porphyrin macrocycle, tetradentate equatorial Fe coordination",
    },
}


# =============================================================================
# KNOWLEDGE BASE CLASS
# =============================================================================

class LigandKnowledgeBase:
    """Self-growing ligand feature knowledge base.
    Follows DesignHistoryManager storage pattern."""

    def __init__(self, base_dir: Optional[str] = None):
        if base_dir is None:
            base_dir = os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                "ligand_knowledge_base",
            )
        self.base_dir = base_dir
        self._ensure_structure()

    def _ensure_structure(self):
        """Create directory structure and seed data on first use."""
        index_path = os.path.join(self.base_dir, "index.json")
        if os.path.exists(index_path):
            return

        os.makedirs(os.path.join(self.base_dir, "ligands"), exist_ok=True)
        # Write seed data
        index = {
            "version": "1.0.0",
            "created": datetime.now().isoformat(),
            "ligands": {},
        }

        for ligand_id, seed in SEED_LIGANDS.items():
            ligand_dir = os.path.join(self.base_dir, "ligands", ligand_id)
            os.makedirs(os.path.join(ligand_dir, "pdb_evidence"), exist_ok=True)

            # Write meta.json
            meta = {
                "ligand_id": ligand_id,
                "smiles": seed["smiles"],
                "aliases": seed.get("aliases", []),
                "pdb_codes": seed.get("pdb_codes", []),
                "heavy_atom_count": seed.get("heavy_atom_count", 0),
                "hsab_character": seed.get("hsab_character", "unknown"),
                "pharmacophore_features": seed.get("pharmacophore_features", []),
                "coordination_summary": seed.get("coordination_summary", {}),
                "notes": seed.get("notes", ""),
            }
            with open(os.path.join(ligand_dir, "meta.json"), "w") as f:
                json.dump(meta, f, indent=2)

            # Update index
            index["ligands"][ligand_id] = {
                "smiles": seed["smiles"],
                "aliases": seed.get("aliases", []),
                "pdb_codes": seed.get("pdb_codes", []),
                "evidence_count": sum(
                    v.get("evidence_count", 0)
                    for v in seed.get("coordination_summary", {}).values()
                ),
                "metals_observed": list(seed.get("coordination_summary", {}).keys()),
                "last_updated": datetime.now().isoformat(),
            }

        with open(index_path, "w") as f:
            json.dump(index, f, indent=2)

        logger.info(f"Initialized ligand knowledge base at {self.base_dir} with {len(SEED_LIGANDS)} seed ligands")

    def _load_index(self) -> Dict[str, Any]:
        index_path = os.path.join(self.base_dir, "index.json")
        if not os.path.exists(index_path):
            return {"version": "1.0.0", "created": datetime.now().isoformat(), "ligands": {}}
        with open(index_path, "r") as f:
            return json.load(f)

    def _save_index(self, index: Dict[str, Any]):
        index_path = os.path.join(self.base_dir, "index.json")
        with open(index_path, "w") as f:
            json.dump(index, f, indent=2)

    def _normalize_name(self, name: str) -> Optional[str]:
        """Normalize ligand name to KB ID. Checks exact match, aliases, and PDB codes."""
        name_lower = name.lower().strip()
        index = self._load_index()

        # Direct match
        if name_lower in index.get("ligands", {}):
            return name_lower

        # Check aliases and PDB codes
        for lid, entry in index.get("ligands", {}).items():
            aliases = [a.lower() for a in entry.get("aliases", [])]
            pdb_codes = [c.lower() for c in entry.get("pdb_codes", [])]
            if name_lower in aliases or name_lower in pdb_codes:
                return lid

        return None

    def lookup(self, ligand_name: str, metal: Optional[str] = None) -> Optional[Dict[str, Any]]:
        """Check knowledge base for cached ligand features.

        Returns meta.json content enriched with metal-specific coordination data.
        """
        lid = self._normalize_name(ligand_name)
        if lid is None:
            return None

        meta_path = os.path.join(self.base_dir, "ligands", lid, "meta.json")
        if not os.path.exists(meta_path):
            return None

        with open(meta_path, "r") as f:
            meta = json.load(f)

        # Enrich with metal-specific data
        if metal:
            metal_upper = metal.upper()
            coord = meta.get("coordination_summary", {}).get(metal_upper)
            if coord:
                meta["metal_coordination"] = coord

        return meta

    def get_coordination_donors(self, ligand_id: str, metal: str) -> Optional[List[str]]:
        """Get consensus coordination donors for a specific metal."""
        meta = self.lookup(ligand_id, metal)
        if meta is None:
            return None

        metal_upper = metal.upper()
        coord = meta.get("coordination_summary", {}).get(metal_upper)
        if coord:
            return coord.get("donors", [])

        # Check PDB evidence files for consensus
        evidence_dir = os.path.join(self.base_dir, "ligands", ligand_id, "pdb_evidence")
        if not os.path.exists(evidence_dir):
            return None

        donor_counts: Dict[str, int] = {}
        evidence_count = 0

        for fname in os.listdir(evidence_dir):
            if not fname.endswith(".json"):
                continue
            # Filter to matching metal
            parts = fname.replace(".json", "").split("_")
            if len(parts) >= 2 and parts[-1].upper() != metal_upper:
                continue

            try:
                with open(os.path.join(evidence_dir, fname), "r") as f:
                    ev = json.load(f)
                for donor_str in ev.get("ligand_donors", []):
                    atom_name = donor_str.split("@")[0] if "@" in donor_str else donor_str
                    donor_counts[atom_name] = donor_counts.get(atom_name, 0) + 1
                evidence_count += 1
            except (json.JSONDecodeError, IOError):
                continue

        if evidence_count == 0:
            return None

        # Consensus: atoms appearing in >50% of evidence
        threshold = evidence_count * 0.5
        consensus = [atom for atom, count in donor_counts.items() if count > threshold]
        return consensus if consensus else None

    def get_coordination_donors_by_smiles(self, smiles: str, metal: str) -> Optional[List[str]]:
        """Get coordination donors by SMILES match (for ligand_donors.py integration)."""
        index = self._load_index()
        for lid, entry in index.get("ligands", {}).items():
            if entry.get("smiles") == smiles:
                return self.get_coordination_donors(lid, metal)
        return None

    def record_evidence(
        self,
        ligand_id: str,
        pdb_id: str,
        metal: str,
        coordination_data: Dict[str, Any],
        smiles: Optional[str] = None,
    ):
        """Cache coordination data from a PDB structure.

        Called after scaffold_search extracts coordination info.
        """
        ligand_id = ligand_id.lower().strip()
        metal_upper = metal.upper()
        pdb_upper = pdb_id.upper()

        ligand_dir = os.path.join(self.base_dir, "ligands", ligand_id)
        evidence_dir = os.path.join(ligand_dir, "pdb_evidence")
        os.makedirs(evidence_dir, exist_ok=True)

        # Save PDB evidence
        evidence = {
            "pdb_id": pdb_upper,
            "metal": metal_upper,
            "resolution": coordination_data.get("resolution"),
            "coordination_number": coordination_data.get("coordination_number", 0),
            "ligand_donors": coordination_data.get("ligand_donors", []),
            "protein_donors": coordination_data.get("protein_donors", []),
            "geometry": coordination_data.get("geometry", "unknown"),
            "date_extracted": datetime.now().isoformat(),
            "source": coordination_data.get("source", "scaffold_search"),
        }
        evidence_path = os.path.join(evidence_dir, f"{pdb_upper}_{metal_upper}.json")
        with open(evidence_path, "w") as f:
            json.dump(evidence, f, indent=2)

        # Update or create meta.json
        meta_path = os.path.join(ligand_dir, "meta.json")
        if os.path.exists(meta_path):
            with open(meta_path, "r") as f:
                meta = json.load(f)
        else:
            meta = {
                "ligand_id": ligand_id,
                "smiles": smiles or "",
                "aliases": [],
                "pdb_codes": [],
                "heavy_atom_count": 0,
                "hsab_character": "unknown",
                "pharmacophore_features": [],
                "coordination_summary": {},
                "notes": "",
            }

        if smiles and not meta.get("smiles"):
            meta["smiles"] = smiles

        # Update coordination summary
        coord_summary = meta.get("coordination_summary", {})
        if metal_upper not in coord_summary:
            coord_summary[metal_upper] = {"donors": [], "denticity": 0, "evidence_count": 0}

        # Parse donor atom names from evidence
        new_donors = []
        for d in coordination_data.get("ligand_donors", []):
            atom_name = d.split("@")[0] if "@" in d else d
            new_donors.append(atom_name)

        # Merge donors (keep unique)
        existing_donors = set(coord_summary[metal_upper].get("donors", []))
        existing_donors.update(new_donors)
        coord_summary[metal_upper]["donors"] = sorted(existing_donors)
        coord_summary[metal_upper]["denticity"] = len(existing_donors)
        coord_summary[metal_upper]["evidence_count"] = coord_summary[metal_upper].get("evidence_count", 0) + 1

        meta["coordination_summary"] = coord_summary

        with open(meta_path, "w") as f:
            json.dump(meta, f, indent=2)

        # Update index
        index = self._load_index()
        if ligand_id not in index.get("ligands", {}):
            index.setdefault("ligands", {})[ligand_id] = {
                "smiles": smiles or "",
                "aliases": [],
                "pdb_codes": [],
                "evidence_count": 0,
                "metals_observed": [],
                "last_updated": datetime.now().isoformat(),
            }

        entry = index["ligands"][ligand_id]
        entry["evidence_count"] = entry.get("evidence_count", 0) + 1
        if metal_upper not in entry.get("metals_observed", []):
            entry.setdefault("metals_observed", []).append(metal_upper)
        entry["last_updated"] = datetime.now().isoformat()
        if smiles and not entry.get("smiles"):
            entry["smiles"] = smiles

        self._save_index(index)

        logger.info(f"Recorded evidence for {ligand_id} + {metal_upper} from {pdb_upper}")

    def list_ligands(self) -> Dict[str, Any]:
        """List all known ligands."""
        return self._load_index().get("ligands", {})


# =============================================================================
# CHEMICALFEATURES PHARMACOPHORE PERCEPTION
# =============================================================================

def get_chemicalfeatures_pharmacophore(smiles: str) -> List[Dict[str, Any]]:
    """RDKit ChemicalFeatures pharmacophore perception.

    Uses BaseFeatures.fdef for donor/acceptor/aromatic/hydrophobic detection.

    Args:
        smiles: SMILES string of the ligand

    Returns:
        List of feature dicts: {atom_idx, atom_name, element, type, coords}
    """
    if not HAS_RDKIT:
        return []

    if not smiles:
        return []

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.warning(f"Failed to parse SMILES for pharmacophore: {smiles}")
        return []

    mol = Chem.AddHs(mol)

    # Generate 3D conformer for coordinates
    coords_available = False
    try:
        params = AllChem.ETKDGv3()
        result = AllChem.EmbedMolecule(mol, params)
        if result == -1:
            params.useRandomCoords = True
            result = AllChem.EmbedMolecule(mol, params)
        if result != -1:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
            coords_available = True
    except Exception:
        pass

    features = []

    if HAS_CHEMICAL_FEATURES:
        # Use RDKit's ChemicalFeatures factory
        try:
            fdef_path = os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
            factory = ChemicalFeatures.BuildFeatureFactory(fdef_path)
            feats = factory.GetFeaturesForMol(mol)

            for feat in feats:
                feat_type = feat.GetType()
                atom_indices = feat.GetAtomIds()

                # Map feature family to our types
                family = feat.GetFamily()
                if family == "Donor":
                    ftype = "donor"
                elif family == "Acceptor":
                    ftype = "acceptor"
                elif family == "Aromatic":
                    ftype = "aromatic"
                elif family == "Hydrophobe":
                    ftype = "hydrophobic"
                else:
                    continue

                for aidx in atom_indices:
                    atom = mol.GetAtomWithIdx(aidx)
                    if atom.GetSymbol() == "H":
                        continue

                    coords = None
                    if coords_available:
                        conf = mol.GetConformer()
                        pos = conf.GetAtomPosition(aidx)
                        coords = [round(pos.x, 3), round(pos.y, 3), round(pos.z, 3)]

                    features.append({
                        "atom_idx": aidx,
                        "atom_name": f"{atom.GetSymbol()}{aidx}",
                        "element": atom.GetSymbol(),
                        "type": ftype,
                        "coords": coords,
                    })
        except Exception as e:
            logger.warning(f"ChemicalFeatures failed: {e}, falling back to SMARTS")
            features = []

    # Fallback: manual donor/acceptor detection from atom types
    if not features:
        for atom in mol.GetAtoms():
            elem = atom.GetSymbol()
            aidx = atom.GetIdx()
            if elem in ("O", "N", "S", "P"):
                # Determine donor/acceptor heuristically
                total_h = atom.GetTotalNumHs()
                is_donor = total_h > 0 or (elem == "N" and atom.GetIsAromatic())
                is_acceptor = elem in ("O", "N", "S")

                coords = None
                if coords_available:
                    conf = mol.GetConformer()
                    pos = conf.GetAtomPosition(aidx)
                    coords = [round(pos.x, 3), round(pos.y, 3), round(pos.z, 3)]

                if is_donor:
                    features.append({
                        "atom_idx": aidx,
                        "atom_name": f"{elem}{aidx}",
                        "element": elem,
                        "type": "donor",
                        "coords": coords,
                    })
                if is_acceptor and not is_donor:
                    features.append({
                        "atom_idx": aidx,
                        "atom_name": f"{elem}{aidx}",
                        "element": elem,
                        "type": "acceptor",
                        "coords": coords,
                    })

    # Deduplicate by atom_idx + type
    seen = set()
    deduped = []
    for f in features:
        key = (f["atom_idx"], f["type"])
        if key not in seen:
            seen.add(key)
            deduped.append(f)

    return deduped


# =============================================================================
# UNIFIED ENTRY POINT
# =============================================================================

def analyze_ligand_features(
    ligand_name: str,
    smiles: Optional[str] = None,
    metal: Optional[str] = None,
    pdb_content: Optional[str] = None,
    record_evidence_data: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Unified entry point for ligand feature analysis.

    Tries layers in order:
    1. Knowledge base lookup
    2. ChemicalFeatures pharmacophore perception
    3. Geometry-filtered donors (existing fallback)

    Args:
        ligand_name: Ligand name or PDB code
        smiles: Optional SMILES string
        metal: Optional metal symbol
        pdb_content: Optional PDB content for 3D analysis
        record_evidence_data: Optional evidence to record (from scaffold search)

    Returns:
        Unified result dict with features, coordination donors, and metadata
    """
    kb = LigandKnowledgeBase()

    # Record evidence if provided (from scaffold search)
    if record_evidence_data:
        try:
            kb.record_evidence(
                ligand_id=ligand_name,
                pdb_id=record_evidence_data.get("pdb_id", ""),
                metal=record_evidence_data.get("metal", metal or ""),
                coordination_data=record_evidence_data,
                smiles=smiles,
            )
        except Exception as e:
            logger.warning(f"Failed to record evidence: {e}")

    result = {
        "ligand_name": ligand_name,
        "smiles": smiles or "",
        "source": "unknown",
        "metal": metal,
        "features": [],
        "coordination_donors": [],
        "max_denticity": 0,
        "evidence_count": 0,
        "compatibility_score": 0.0,
        "coordination_mode": "unknown",
        "notes": "",
        "pdb_evidence": [],
        "ligand_pdb_content": pdb_content,
    }

    # Layer 1: Knowledge base
    kb_data = kb.lookup(ligand_name, metal)
    if kb_data:
        logger.info(f"Knowledge base hit for {ligand_name}")
        result["source"] = "knowledge_base"
        result["notes"] = kb_data.get("notes", "")

        # Build features from KB pharmacophore data
        kb_features = kb_data.get("pharmacophore_features", [])
        coord_donors = []
        metal_coord = kb_data.get("metal_coordination") or {}
        kb_donor_names = set(metal_coord.get("donors", []))

        for f in kb_features:
            is_coord = f.get("atom_name", "") in kb_donor_names
            feature = {
                "atom_idx": f.get("atom_idx", 0),
                "atom_name": f.get("atom_name", ""),
                "element": f.get("element", ""),
                "type": f.get("type", "donor"),
                "is_coordination_donor": is_coord,
                "coords": f.get("coords"),
                "hsab": kb_data.get("hsab_character"),
                "enabled": True,
            }
            result["features"].append(feature)
            if is_coord:
                coord_donors.append(f.get("atom_name", ""))

        result["coordination_donors"] = coord_donors
        result["max_denticity"] = metal_coord.get("denticity", len(coord_donors))
        result["evidence_count"] = metal_coord.get("evidence_count", 0)
        result["smiles"] = kb_data.get("smiles") or smiles or ""

        # List PDB evidence
        lid = kb._normalize_name(ligand_name)
        if lid:
            evidence_dir = os.path.join(kb.base_dir, "ligands", lid, "pdb_evidence")
            if os.path.exists(evidence_dir):
                result["pdb_evidence"] = [
                    f.replace(".json", "") for f in os.listdir(evidence_dir) if f.endswith(".json")
                ]

        # Compute compatibility score
        if metal:
            result["compatibility_score"] = _compute_compatibility(kb_data, metal)
            if coord_donors:
                result["coordination_mode"] = _get_coordination_mode(len(coord_donors))
            else:
                # KB has the ligand but no specific coordination data for this metal.
                # Mark donors as best-guess from pharmacophore (O/N/S donor/acceptor atoms).
                for feat in result["features"]:
                    if feat["type"] in ("donor", "acceptor") and feat["element"] in ("O", "N", "S"):
                        feat["is_coordination_donor"] = True
                        coord_donors.append(feat["atom_name"])
                result["coordination_donors"] = coord_donors
                result["max_denticity"] = len(coord_donors)
                result["coordination_mode"] = _get_coordination_mode(len(coord_donors))
                result["notes"] = (result["notes"] + " No crystal evidence for this metal; "
                                   "donors inferred from pharmacophore features.").strip()

        return result

    # Layer 2: ChemicalFeatures
    if smiles:
        logger.info(f"Trying ChemicalFeatures for {ligand_name}")
        cf_features = get_chemicalfeatures_pharmacophore(smiles)
        if cf_features:
            result["source"] = "chemicalfeatures"
            result["smiles"] = smiles

            for f in cf_features:
                feature = {
                    "atom_idx": f["atom_idx"],
                    "atom_name": f["atom_name"],
                    "element": f["element"],
                    "type": f["type"],
                    "is_coordination_donor": f["type"] in ("donor", "acceptor") and f["element"] in ("O", "N", "S"),
                    "coords": f.get("coords"),
                    "hsab": _element_to_hsab(f["element"]),
                    "enabled": True,
                }
                result["features"].append(feature)

            # Identify likely coordination donors (O, N, S that are donors/acceptors)
            coord_donors = [
                f["atom_name"]
                for f in result["features"]
                if f["is_coordination_donor"]
            ]
            result["coordination_donors"] = coord_donors
            result["max_denticity"] = len(coord_donors)

            if metal:
                result["compatibility_score"] = _compute_compatibility_from_features(result["features"], metal)
                result["coordination_mode"] = _get_coordination_mode(len(coord_donors))

            return result

    # Layer 3: Geometry-filtered donors (existing fallback)
    if smiles:
        logger.info(f"Falling back to geometry filter for {ligand_name}")
        try:
            from ligand_donors import identify_donors_from_smiles, filter_donors_by_geometry, score_ligand_metal_compatibility

            raw_donors = identify_donors_from_smiles(smiles)
            donors = filter_donors_by_geometry(raw_donors, smiles, metal=metal) if metal else raw_donors

            result["source"] = "geometry_filter"
            result["smiles"] = smiles

            for d in donors:
                feature = {
                    "atom_idx": d["atom_idx"],
                    "atom_name": f"{d['element']}{d['atom_idx']}",
                    "element": d["element"],
                    "type": d["type"],
                    "is_coordination_donor": True,
                    "coords": d.get("coords"),
                    "hsab": d.get("hsab"),
                    "enabled": True,
                }
                result["features"].append(feature)

            result["coordination_donors"] = [f["atom_name"] for f in result["features"]]
            result["max_denticity"] = len(result["coordination_donors"])

            if metal:
                result["compatibility_score"] = score_ligand_metal_compatibility(smiles, metal)
                result["coordination_mode"] = _get_coordination_mode(len(donors))

        except ImportError:
            logger.warning("ligand_donors module not available for fallback")

    return result


# =============================================================================
# HELPERS
# =============================================================================

def _element_to_hsab(element: str) -> str:
    """Simple element to HSAB classification."""
    mapping = {
        "O": "hard",
        "N": "borderline",
        "S": "soft",
        "P": "hard",
    }
    return mapping.get(element, "unknown")


def _get_coordination_mode(num_donors: int) -> str:
    if num_donors >= 3:
        return "polydentate"
    elif num_donors == 2:
        return "bidentate"
    elif num_donors == 1:
        return "monodentate"
    return "unknown"


def _compute_compatibility(kb_data: Dict[str, Any], metal: str) -> float:
    """Compute compatibility score from KB data."""
    hsab_char = kb_data.get("hsab_character", "unknown")
    try:
        from metal_chemistry import get_hsab_class, METAL_DATABASE
        metal_upper = metal.upper()
        if metal_upper in METAL_DATABASE:
            ox = METAL_DATABASE[metal_upper]["default_oxidation"]
            metal_hsab = get_hsab_class(metal, ox)
        else:
            metal_hsab = "borderline"
    except (ImportError, ValueError):
        metal_hsab = "borderline"

    compatibility_matrix = {
        "hard": {"hard": 1.0, "borderline": 0.7, "soft": 0.1},
        "borderline": {"hard": 0.7, "borderline": 1.0, "soft": 0.7},
        "soft": {"hard": 0.1, "borderline": 0.7, "soft": 1.0},
    }
    return compatibility_matrix.get(metal_hsab, {}).get(hsab_char, 0.5)


def _compute_compatibility_from_features(features: List[Dict[str, Any]], metal: str) -> float:
    """Compute compatibility from feature list."""
    if not features:
        return 0.0

    try:
        from metal_chemistry import get_hsab_class, METAL_DATABASE
        metal_upper = metal.upper()
        if metal_upper in METAL_DATABASE:
            ox = METAL_DATABASE[metal_upper]["default_oxidation"]
            metal_hsab = get_hsab_class(metal, ox)
        else:
            metal_hsab = "borderline"
    except (ImportError, ValueError):
        metal_hsab = "borderline"

    compatibility_matrix = {
        "hard": {"hard": 1.0, "borderline": 0.7, "soft": 0.1},
        "borderline": {"hard": 0.7, "borderline": 1.0, "soft": 0.7},
        "soft": {"hard": 0.1, "borderline": 0.7, "soft": 1.0},
    }

    coord_features = [f for f in features if f.get("is_coordination_donor")]
    if not coord_features:
        return 0.5

    total = sum(
        compatibility_matrix.get(metal_hsab, {}).get(f.get("hsab", "borderline"), 0.5)
        for f in coord_features
    )
    return total / len(coord_features)
