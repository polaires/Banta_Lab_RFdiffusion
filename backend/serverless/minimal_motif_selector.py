"""
Minimal Motif Selector — Evidence-Based Residue Selection for RFD3.

Combines evidence from M-CSA, PLIP/distance-based interactions, enzyme chemistry DB,
and optionally ConSurf conservation to classify pocket residues into tiers:

| Tier | Criteria                      | Fix Type | Expected Count |
|------|-------------------------------|----------|----------------|
| 1    | M-CSA primary catalytic role  | TIP      | 2-4            |
| 2    | >=2 strong sources agree      | TIP      | 2-4            |
| 3    | 1 strong source / ConSurf 3   | CA       | 0-5            |
| 4    | Proximity only                | NONE     | remainder       |

Only Tier 1-3 enter the motif (<=12 residues total). Tier 4 are excluded.
"""

import logging
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)


# =============================================================================
# TIP ATOM LOOKUP — functional group atoms per amino acid
# =============================================================================

TIP_ATOMS: Dict[str, List[str]] = {
    "ASP": ["OD1", "OD2", "CG"],
    "GLU": ["OE1", "OE2", "CD"],
    "HIS": ["ND1", "NE2", "CE1", "CD2", "CG"],
    "LYS": ["NZ", "CE"],
    "ARG": ["NH1", "NH2", "NE", "CZ"],
    "SER": ["OG", "CB"],
    "THR": ["OG1", "CB", "CG2"],
    "CYS": ["SG", "CB"],
    "TYR": ["OH", "CZ", "CE1", "CE2"],
    "ASN": ["OD1", "ND2", "CG"],
    "GLN": ["OE1", "NE2", "CD"],
    "TRP": ["NE1", "CE2", "CD1", "CZ2", "CH2", "CZ3", "CE3", "CD2", "CG"],
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "ALA": ["CB"],
    "VAL": ["CG1", "CG2", "CB"],
    "LEU": ["CD1", "CD2", "CG"],
    "ILE": ["CD1", "CG1", "CG2"],
    "MET": ["SD", "CE", "CG"],
    "PRO": ["CG", "CD"],
    "GLY": [],
}

# M-CSA roles considered "primary catalytic"
MCSA_PRIMARY_ROLES = {
    "nucleophile",
    "proton donor",
    "proton acceptor",
    "metal ligand",
    "electrostatic stabiliser",
    "steric role",
    "activator",
}


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class MotifResidue:
    """A single residue classified by the motif selector."""
    chain: str
    resnum: int
    resname: str          # 3-letter code
    tier: int             # 1-4
    fix_type: str         # "TIP", "CA", "BKBN", "NONE"
    evidence: List[str]   # e.g. ["mcsa:nucleophile", "plip:hbond", "consurf:grade2"]
    functional_atoms: List[str]  # Specific atoms to fix
    confidence: float = 0.0


@dataclass
class MinimalMotifResult:
    """Result from MinimalMotifSelector.select()."""
    motif_residues: List[MotifResidue]     # Tier 1-3
    context_residues: List[MotifResidue]   # Tier 4 (not fixed)
    total_fixed_atoms: int = 0
    evidence_summary: Dict[str, int] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)

    def to_fixed_atoms_dict(self) -> Dict[str, str]:
        """Convert to RFD3 select_fixed_atoms format.

        Returns dict like {"A1": "OD1,OD2,CG", "A2": "NE2,ND1", "A3": "CA"}.
        Keys use the renumbered residue IDs (A1, A2, ...) based on position
        in motif_residues list.
        """
        result = {}
        for idx, res in enumerate(self.motif_residues, start=1):
            res_id = f"A{idx}"
            if res.fix_type == "TIP" and res.functional_atoms:
                result[res_id] = ",".join(res.functional_atoms)
            elif res.fix_type == "CA":
                result[res_id] = "CA"
            elif res.fix_type == "BKBN":
                result[res_id] = "N,CA,C,O"
            # fix_type == "NONE" -> not added
        return result

    def get_motif_residue_keys(self) -> Set[Tuple[str, int]]:
        """Get set of (chain, resnum) for motif residues only."""
        return {(r.chain, r.resnum) for r in self.motif_residues}


# =============================================================================
# MAIN CLASS
# =============================================================================

class MinimalMotifSelector:
    """
    Selects a minimal motif from pocket residues using multi-source evidence.

    Usage:
        selector = MinimalMotifSelector()
        result = await selector.select(
            pdb_id="4CVB", pdb_content=pdb_str,
            ligand_code="PQQ", metal_type="CA",
            enzyme_class="quinoprotein", chain="A",
            pocket_residues=[{"chain": "A", "resnum": 356, "res_name": "ASP"}, ...]
        )
        fixed_atoms = result.to_fixed_atoms_dict()
    """

    def __init__(
        self,
        max_tier1_2: int = 7,
        max_tier3: int = 5,
        max_total: int = 12,
        use_conservation: bool = False,
        conservation_timeout: float = 60.0,
    ):
        self.max_tier1_2 = max_tier1_2
        self.max_tier3 = max_tier3
        self.max_total = max_total
        self.use_conservation = use_conservation
        self.conservation_timeout = conservation_timeout

    async def select(
        self,
        pdb_id: str,
        pdb_content: str,
        ligand_code: Optional[str] = None,
        metal_type: Optional[str] = None,
        enzyme_class: Optional[str] = None,
        chain: str = "A",
        pocket_residues: Optional[List[Dict]] = None,
    ) -> MinimalMotifResult:
        """
        Main entry: gather evidence -> merge -> classify -> enforce limits.

        Args:
            pdb_id: 4-char PDB ID
            pdb_content: Full PDB file content
            ligand_code: Ligand 3-letter code (e.g. "PQQ")
            metal_type: Metal element (e.g. "CA")
            enzyme_class: Key from enzyme_chemistry.ENZYME_CLASS_DATABASE
            chain: Protein chain to analyze
            pocket_residues: Pre-identified pocket residues (list of dicts with
                            chain, resnum, res_name keys)

        Returns:
            MinimalMotifResult with classified residues and fixed_atoms mapping
        """
        warnings = []

        # Build pocket residue set for filtering
        pocket_keys = set()
        pocket_resname = {}
        if pocket_residues:
            for r in pocket_residues:
                c = r.get("chain", chain)
                rn = r.get("resnum", 0)
                pocket_keys.add((c, rn))
                pocket_resname[(c, rn)] = r.get("res_name", "UNK")

        # Gather evidence from all sources (in parallel where possible)
        mcsa_ev = await self._gather_mcsa_evidence(pdb_id)
        plip_ev = self._gather_plip_evidence(pdb_content, ligand_code)
        chem_ev = self._gather_enzyme_chemistry_evidence(
            pdb_content, enzyme_class, metal_type, ligand_code,
            pocket_residues=pocket_residues,
        )

        cons_ev: Dict = {}
        if self.use_conservation:
            try:
                cons_ev = await self._gather_conservation_evidence(pdb_content, chain)
            except Exception as e:
                warnings.append(f"Conservation analysis failed: {e}")

        # Log evidence counts
        logger.info(
            f"Evidence gathered: mcsa={len(mcsa_ev)}, plip={len(plip_ev)}, "
            f"chem={len(chem_ev)}, cons={len(cons_ev)}"
        )

        # Merge and classify
        all_residues = self._merge_evidence(
            mcsa_ev, plip_ev, chem_ev, cons_ev,
            pocket_keys, pocket_resname, chain
        )

        # Enforce limits
        motif, context = self._enforce_limits(all_residues)

        # If all sources failed and we got nothing, use proximity fallback
        if not motif and pocket_residues:
            warnings.append("All evidence sources empty; using proximity fallback")
            motif, context = self._proximity_fallback(pocket_residues, chain)

        # Ligand-proximity upgrade: ensure exposed ligand atoms have nearby
        # TIP-fixed residues for physical shielding. Upgrades CA-only or
        # promotes context residues near unshielded ligand atoms.
        if ligand_code and pdb_content and motif:
            motif, context, promo_warns = self._upgrade_ligand_shielding(
                pdb_content, ligand_code, motif, context, chain
            )
            warnings.extend(promo_warns)

        # Count fixed atoms
        total_fixed = 0
        for r in motif:
            if r.fix_type == "TIP":
                total_fixed += len(r.functional_atoms)
            elif r.fix_type == "CA":
                total_fixed += 1
            elif r.fix_type == "BKBN":
                total_fixed += 4

        # Build evidence summary
        ev_summary: Dict[str, int] = {}
        for r in motif + context:
            for ev in r.evidence:
                src = ev.split(":")[0]
                ev_summary[src] = ev_summary.get(src, 0) + 1

        logger.info(
            f"Motif selection: {len(motif)} motif residues, "
            f"{len(context)} context, {total_fixed} fixed atoms"
        )

        return MinimalMotifResult(
            motif_residues=motif,
            context_residues=context,
            total_fixed_atoms=total_fixed,
            evidence_summary=ev_summary,
            warnings=warnings,
        )

    # =========================================================================
    # EVIDENCE GATHERING
    # =========================================================================

    async def _gather_mcsa_evidence(self, pdb_id: str) -> Dict[Tuple[str, int], Dict]:
        """Query M-CSA for catalytic residues."""
        result: Dict[Tuple[str, int], Dict] = {}
        try:
            from mcsa_client_lite import query_mcsa_by_pdb
            residues = await query_mcsa_by_pdb(pdb_id)
            for r in residues:
                key = (r["chain"], r["resnum"])
                role = (r.get("role") or "").lower()
                # Determine if primary catalytic
                is_primary = any(pr in role for pr in MCSA_PRIMARY_ROLES)
                result[key] = {
                    "resname": r["resname"],
                    "evidence_tags": [f"mcsa:{role}" if role else "mcsa:unknown"],
                    "confidence": r.get("confidence", 1.0),
                    "is_mcsa_primary": is_primary,
                }
        except Exception as e:
            logger.warning(f"M-CSA evidence gathering failed: {e}")
        return result

    def _gather_plip_evidence(
        self, pdb_content: str, ligand_code: Optional[str]
    ) -> Dict[Tuple[str, int], Dict]:
        """Gather PLIP/distance-based interaction evidence."""
        result: Dict[Tuple[str, int], Dict] = {}
        if not ligand_code:
            return result

        try:
            from shared.interaction_analysis import analyze_all_interactions
            summary = analyze_all_interactions(
                pdb_content, ligand_name=ligand_code.upper(),
                include_visualization_data=False,
            )

            if summary.status != "completed":
                logger.warning(f"PLIP analysis status: {summary.status}")
                return result

            # Process H-bonds (strong evidence)
            for hb in summary.hydrogen_bonds:
                key = self._parse_residue_key(hb.protein_chain, hb.protein_residue)
                if key:
                    if key not in result:
                        result[key] = {
                            "resname": self._extract_resname(hb.protein_residue),
                            "evidence_tags": [],
                            "confidence": 0.8,
                        }
                    result[key]["evidence_tags"].append("plip:hbond")

            # Process salt bridges (strong evidence)
            for sb in summary.salt_bridges:
                key = self._parse_residue_key(sb.protein_chain, sb.protein_residue)
                if key:
                    if key not in result:
                        result[key] = {
                            "resname": self._extract_resname(sb.protein_residue),
                            "evidence_tags": [],
                            "confidence": 0.8,
                        }
                    result[key]["evidence_tags"].append("plip:salt_bridge")

            # Process pi-stacking (moderate evidence)
            for ps in summary.pi_stacking:
                key = self._parse_residue_key(ps.protein_chain, ps.protein_residue)
                if key:
                    if key not in result:
                        result[key] = {
                            "resname": self._extract_resname(ps.protein_residue),
                            "evidence_tags": [],
                            "confidence": 0.6,
                        }
                    result[key]["evidence_tags"].append("plip:pi_stack")

            # Process hydrophobic contacts (weak evidence)
            for hc in summary.hydrophobic_contacts:
                key = self._parse_residue_key(hc.protein_chain, hc.protein_residue)
                if key:
                    if key not in result:
                        result[key] = {
                            "resname": self._extract_resname(hc.protein_residue),
                            "evidence_tags": [],
                            "confidence": 0.4,
                        }
                    result[key]["evidence_tags"].append("plip:hydrophobic")

        except Exception as e:
            logger.warning(f"PLIP evidence gathering failed: {e}")

        return result

    def _gather_enzyme_chemistry_evidence(
        self,
        pdb_content: str,
        enzyme_class: Optional[str],
        metal_type: Optional[str],
        ligand_code: Optional[str],
        pocket_residues: Optional[List[Dict]] = None,
    ) -> Dict[Tuple[str, int], Dict]:
        """Gather enzyme chemistry DB evidence.

        Uses two approaches:
        1. identify_catalytic_residues_from_pdb (CA-distance to HETATM)
        2. Cross-reference pocket residues against essential residue types
           (catches residues whose CA is far but sidechain is close)
        """
        result: Dict[Tuple[str, int], Dict] = {}

        # Try to detect enzyme class if not provided
        if not enzyme_class:
            try:
                from enzyme_chemistry import detect_enzyme_class
                detected, conf = detect_enzyme_class(
                    "", ligand_name=ligand_code, metal_type=metal_type
                )
                if detected and conf >= 0.5:
                    enzyme_class = detected
                    logger.info(f"Auto-detected enzyme class: {enzyme_class} (conf={conf:.2f})")
            except ImportError:
                pass

        if not enzyme_class:
            return result

        try:
            from enzyme_chemistry import (
                identify_catalytic_residues_from_pdb,
                ENZYME_CLASS_DATABASE,
                AA_CODES_REV,
            )

            # Approach 1: CA-distance based identification
            catalytic = identify_catalytic_residues_from_pdb(
                pdb_content, enzyme_class,
                metal_type=metal_type, ligand_name=ligand_code,
            )
            for r in catalytic:
                key = (r["chain"], r["resnum"])
                role = r.get("role", "unknown")
                conf = r.get("confidence", 0.6)
                tag = f"enzyme_chem:{role}"
                is_essential = conf >= 0.8

                if key not in result:
                    result[key] = {
                        "resname": r["resname"],
                        "evidence_tags": [],
                        "confidence": conf,
                        "is_enzyme_essential": is_essential,
                    }
                result[key]["evidence_tags"].append(tag)

            # Approach 2: Cross-reference pocket residues with essential types
            # This catches residues whose sidechain is close but CA is far
            profile = ENZYME_CLASS_DATABASE.get(enzyme_class)
            if profile and pocket_residues:
                essential_1letter = set(profile.essential_residues)
                for r in pocket_residues:
                    c = r.get("chain", "A")
                    rn = r.get("resnum", 0)
                    resname = r.get("res_name", "UNK")
                    key = (c, rn)

                    if key in result:
                        continue  # Already found by CA-distance approach

                    res_1letter = AA_CODES_REV.get(resname, "X")
                    if res_1letter in essential_1letter:
                        result[key] = {
                            "resname": resname,
                            "evidence_tags": [f"enzyme_chem:essential_type_{res_1letter}"],
                            "confidence": 0.85,
                            "is_enzyme_essential": True,
                        }
                        logger.debug(
                            f"Pocket residue {c}{rn} {resname} matches "
                            f"essential type for {enzyme_class}"
                        )

        except Exception as e:
            logger.warning(f"Enzyme chemistry evidence failed: {e}")

        return result

    async def _gather_conservation_evidence(
        self, pdb_content: str, chain: str
    ) -> Dict[Tuple[str, int], Dict]:
        """Gather conservation evidence (placeholder for ConSurf integration)."""
        # ConSurf integration is optional and often unavailable in serverless.
        # This is a placeholder that returns empty.
        # A real implementation would query ConSurf API or run rate4site.
        return {}

    # =========================================================================
    # MERGING & CLASSIFICATION
    # =========================================================================

    def _merge_evidence(
        self,
        mcsa_ev: Dict[Tuple[str, int], Dict],
        plip_ev: Dict[Tuple[str, int], Dict],
        chem_ev: Dict[Tuple[str, int], Dict],
        cons_ev: Dict[Tuple[str, int], Dict],
        pocket_keys: Set[Tuple[str, int]],
        pocket_resname: Dict[Tuple[str, int], str],
        chain: str,
    ) -> List[MotifResidue]:
        """Merge all evidence sources and classify each residue."""

        # Collect all residue keys from all sources
        all_keys: Set[Tuple[str, int]] = set()
        all_keys.update(mcsa_ev.keys())
        all_keys.update(plip_ev.keys())
        all_keys.update(chem_ev.keys())
        all_keys.update(cons_ev.keys())
        # Also include pocket residues (they may have no evidence -> Tier 4)
        all_keys.update(pocket_keys)

        residues: List[MotifResidue] = []

        for key in all_keys:
            c, rn = key
            # Collect all evidence tags
            all_tags: List[str] = []
            resname = "UNK"
            has_mcsa_primary = False

            if key in mcsa_ev:
                all_tags.extend(mcsa_ev[key]["evidence_tags"])
                resname = mcsa_ev[key].get("resname", resname)
                has_mcsa_primary = mcsa_ev[key].get("is_mcsa_primary", False)

            if key in plip_ev:
                all_tags.extend(plip_ev[key]["evidence_tags"])
                if resname == "UNK":
                    resname = plip_ev[key].get("resname", resname)

            if key in chem_ev:
                all_tags.extend(chem_ev[key]["evidence_tags"])
                if resname == "UNK":
                    resname = chem_ev[key].get("resname", resname)

            if key in cons_ev:
                all_tags.extend(cons_ev[key]["evidence_tags"])

            # Fall back to pocket_resname
            if resname == "UNK" and key in pocket_resname:
                resname = pocket_resname[key]

            # Count strong sources
            strong_count = self._count_strong_sources(all_tags, mcsa_ev.get(key), chem_ev.get(key))

            # Classify tier
            tier = self._classify_tier(all_tags, strong_count, has_mcsa_primary)

            # Map tier to fix_type
            if tier <= 2:
                fix_type = "TIP"
            elif tier == 3:
                fix_type = "CA"
            else:
                fix_type = "NONE"

            # Get functional atoms
            functional_atoms = TIP_ATOMS.get(resname, []) if fix_type == "TIP" else []
            if fix_type == "CA":
                functional_atoms = ["CA"]

            # Compute confidence from evidence
            confidence = min(1.0, 0.2 * strong_count + 0.1 * len(all_tags))
            if has_mcsa_primary:
                confidence = max(confidence, 0.95)

            residues.append(MotifResidue(
                chain=c,
                resnum=rn,
                resname=resname,
                tier=tier,
                fix_type=fix_type,
                evidence=all_tags,
                functional_atoms=functional_atoms,
                confidence=confidence,
            ))

        # Sort by tier (ascending), then confidence (descending)
        residues.sort(key=lambda r: (r.tier, -r.confidence, r.resnum))
        return residues

    def _count_strong_sources(
        self,
        tags: List[str],
        mcsa_data: Optional[Dict],
        chem_data: Optional[Dict],
    ) -> int:
        """Count number of 'strong' evidence sources."""
        strong = 0
        has_mcsa = any(t.startswith("mcsa:") for t in tags)
        has_plip_hbond = any("plip:hbond" in t for t in tags)
        has_plip_salt = any("plip:salt_bridge" in t for t in tags)
        has_chem_essential = chem_data.get("is_enzyme_essential", False) if chem_data else False
        has_cons_high = any("consurf:grade1" in t or "consurf:grade2" in t for t in tags)

        if has_mcsa:
            strong += 1
        if has_plip_hbond or has_plip_salt:
            strong += 1
        if has_chem_essential:
            strong += 1
        if has_cons_high:
            strong += 1

        return strong

    def _classify_tier(
        self,
        evidence_tags: List[str],
        strong_count: int,
        has_mcsa_primary: bool,
    ) -> int:
        """Classify a residue into tier 1-4."""
        # Tier 1: M-CSA primary catalytic role
        if has_mcsa_primary:
            return 1

        # Tier 2: >= 2 strong sources agree
        if strong_count >= 2:
            return 2

        # Tier 3: 1 strong source, or ConSurf grade 3, or PLIP hydrophobic
        if strong_count >= 1:
            return 3
        if any("consurf:grade3" in t for t in evidence_tags):
            return 3
        if any("plip:hydrophobic" in t for t in evidence_tags):
            return 3

        # Tier 4: everything else
        return 4

    def _enforce_limits(
        self, residues: List[MotifResidue]
    ) -> Tuple[List[MotifResidue], List[MotifResidue]]:
        """Enforce max limits. Excess residues get demoted to context."""
        tier1_2 = [r for r in residues if r.tier <= 2]
        tier3 = [r for r in residues if r.tier == 3]
        tier4 = [r for r in residues if r.tier >= 4]

        # Cap Tier 1+2
        if len(tier1_2) > self.max_tier1_2:
            # Demote lowest-confidence Tier 1+2 to Tier 3
            excess = tier1_2[self.max_tier1_2:]
            for r in excess:
                r.tier = 3
                r.fix_type = "CA"
                r.functional_atoms = ["CA"]
            tier3 = excess + tier3
            tier1_2 = tier1_2[:self.max_tier1_2]

        # Cap Tier 3
        if len(tier3) > self.max_tier3:
            excess = tier3[self.max_tier3:]
            for r in excess:
                r.tier = 4
                r.fix_type = "NONE"
                r.functional_atoms = []
            tier4 = excess + tier4
            tier3 = tier3[:self.max_tier3]

        motif = tier1_2 + tier3

        # Cap total
        if len(motif) > self.max_total:
            excess = motif[self.max_total:]
            for r in excess:
                r.tier = 4
                r.fix_type = "NONE"
                r.functional_atoms = []
            tier4 = excess + tier4
            motif = motif[:self.max_total]

        return motif, tier4

    def _upgrade_ligand_shielding(
        self,
        pdb_content: str,
        ligand_code: str,
        motif: List[MotifResidue],
        context: List[MotifResidue],
        chain: str,
        shield_radius: float = 5.0,
    ) -> Tuple[List[MotifResidue], List[MotifResidue], List[str]]:
        """
        Ensure exposed ligand atoms have nearby TIP-fixed residues.

        For each ligand heavy atom, check if any TIP-fixed motif residue
        has a sidechain atom within shield_radius. If not, the atom is
        "exposed" — upgrade the nearest CA-only motif residue to TIP,
        or promote the nearest context residue to Tier 3 TIP.

        This prevents backbone-through-ligand clashes by ensuring physical
        shielding around all parts of the ligand.
        """
        warnings: List[str] = []

        # Parse ligand atom positions
        ligand_atoms = []  # [(atom_name, x, y, z), ...]
        for line in pdb_content.split('\n'):
            if line.startswith('HETATM'):
                resname = line[17:20].strip()
                if resname == ligand_code:
                    atom_name = line[12:16].strip()
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        ligand_atoms.append((atom_name, x, y, z))
                    except (ValueError, IndexError):
                        continue

        if not ligand_atoms:
            return motif, context, warnings

        # Parse protein atom positions (sidechain atoms only, not backbone)
        backbone_names = {"N", "CA", "C", "O"}
        # Map: (chain, resnum) -> [(atom_name, x, y, z), ...]
        residue_atoms: Dict[Tuple[str, int], List[Tuple[str, float, float, float]]] = {}
        for line in pdb_content.split('\n'):
            if line.startswith('ATOM'):
                atom_name = line[12:16].strip()
                ch = line[21].strip() or chain
                try:
                    resnum = int(line[22:26].strip())
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except (ValueError, IndexError):
                    continue
                key = (ch, resnum)
                if key not in residue_atoms:
                    residue_atoms[key] = []
                residue_atoms[key].append((atom_name, x, y, z))

        # Build set of TIP-fixed residue keys
        tip_keys = {(r.chain, r.resnum) for r in motif if r.fix_type == "TIP"}

        # For each ligand atom, check if any TIP-fixed residue has a
        # sidechain atom within shield_radius
        exposed_ligand_atoms = []  # [(atom_name, x, y, z)]
        for la_name, lx, ly, lz in ligand_atoms:
            shielded = False
            for rkey in tip_keys:
                atoms = residue_atoms.get(rkey, [])
                for a_name, ax, ay, az in atoms:
                    if a_name in backbone_names:
                        continue
                    dist = math.sqrt((lx - ax)**2 + (ly - ay)**2 + (lz - az)**2)
                    if dist < shield_radius:
                        shielded = True
                        break
                if shielded:
                    break
            if not shielded:
                exposed_ligand_atoms.append((la_name, lx, ly, lz))

        if not exposed_ligand_atoms:
            return motif, context, warnings

        exposed_names = [a[0] for a in exposed_ligand_atoms]
        logger.info(
            f"Ligand shielding: {len(exposed_ligand_atoms)}/{len(ligand_atoms)} "
            f"atoms exposed (no TIP residue within {shield_radius}A): "
            f"{exposed_names[:10]}"
        )

        # Find which residues are closest to exposed ligand atoms
        # Check CA-only motif residues first, then context residues
        ca_only = [r for r in motif if r.fix_type == "CA"]
        upgradable = ca_only + context

        upgraded_keys: Set[Tuple[str, int]] = set()

        for la_name, lx, ly, lz in exposed_ligand_atoms:
            # Find closest upgradable residue (by any sidechain atom)
            best_res = None
            best_dist = float('inf')

            for r in upgradable:
                rkey = (r.chain, r.resnum)
                if rkey in upgraded_keys:
                    continue
                atoms = residue_atoms.get(rkey, [])
                for a_name, ax, ay, az in atoms:
                    if a_name in backbone_names:
                        continue
                    dist = math.sqrt((lx - ax)**2 + (ly - ay)**2 + (lz - az)**2)
                    if dist < best_dist:
                        best_dist = dist
                        best_res = r

            if best_res and best_dist < shield_radius + 2.0:
                rkey = (best_res.chain, best_res.resnum)
                upgraded_keys.add(rkey)

                tip_atoms = TIP_ATOMS.get(best_res.resname, ["CB"])
                was_context = best_res.fix_type == "NONE"

                best_res.fix_type = "TIP"
                best_res.functional_atoms = tip_atoms
                if best_res.tier > 3:
                    best_res.tier = 3
                if "shielding:ligand_proximity" not in best_res.evidence:
                    best_res.evidence.append("shielding:ligand_proximity")

                action = "promoted from context" if was_context else "upgraded CA→TIP"
                warnings.append(
                    f"Shielding upgrade: {best_res.chain}{best_res.resnum} "
                    f"{best_res.resname} {action} (nearest exposed ligand atom "
                    f"{la_name} at {best_dist:.1f}A)"
                )
                logger.info(warnings[-1])

        # Rebuild motif/context lists after promotions
        if upgraded_keys:
            new_motif = []
            new_context = []
            all_res = motif + context
            for r in all_res:
                if r.tier <= 3:
                    new_motif.append(r)
                else:
                    new_context.append(r)
            motif = new_motif
            context = new_context

        return motif, context, warnings

    def _proximity_fallback(
        self,
        pocket_residues: List[Dict],
        chain: str,
    ) -> Tuple[List[MotifResidue], List[MotifResidue]]:
        """
        Fallback when all evidence sources fail.
        Pick 5 closest pocket residues as Tier 3 (CA only).
        """
        # Sort by distance if available, otherwise just take first 5
        sorted_res = sorted(
            pocket_residues,
            key=lambda r: r.get("min_distance", r.get("resnum", 0))
        )

        motif: List[MotifResidue] = []
        context: List[MotifResidue] = []

        for i, r in enumerate(sorted_res):
            c = r.get("chain", chain)
            rn = r.get("resnum", 0)
            resname = r.get("res_name", "UNK")

            if i < 5:
                motif.append(MotifResidue(
                    chain=c, resnum=rn, resname=resname,
                    tier=3, fix_type="CA",
                    evidence=["proximity:fallback"],
                    functional_atoms=["CA"],
                    confidence=0.3,
                ))
            else:
                context.append(MotifResidue(
                    chain=c, resnum=rn, resname=resname,
                    tier=4, fix_type="NONE",
                    evidence=["proximity:excluded"],
                    functional_atoms=[],
                    confidence=0.1,
                ))

        return motif, context

    # =========================================================================
    # HELPERS
    # =========================================================================

    @staticmethod
    def _parse_residue_key(
        chain: str, residue_str: str
    ) -> Optional[Tuple[str, int]]:
        """Parse 'ASP356' or 'A:ASP356' into (chain, resnum).

        Forces plain Python str for chain to avoid numpy.str_ hash mismatches.
        """
        import re
        # Try "RESNAME123" format
        match = re.match(r'([A-Z]{3})(\d+)', str(residue_str))
        if match:
            return (str(chain), int(match.group(2)))
        return None

    @staticmethod
    def _extract_resname(residue_str: str) -> str:
        """Extract 3-letter code from 'ASP356' etc."""
        import re
        match = re.match(r'([A-Z]{3})\d+', residue_str)
        return match.group(1) if match else "UNK"
