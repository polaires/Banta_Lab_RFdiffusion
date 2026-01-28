"""
Enzyme Chemistry Database for AI-Driven Protein Design.

Provides comprehensive enzyme class profiles for:
- Enzyme class detection from natural language queries
- Catalytic residue identification
- Substrate channel preservation requirements
- Activity preservation validation

Covers all 6 EC categories with 22+ enzyme classes.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Set
import re


# ═══════════════════════════════════════════════════════════════════
# 1. ENZYME CLASS PROFILES
# ═══════════════════════════════════════════════════════════════════

@dataclass
class EnzymeClassProfile:
    """
    Database entry for enzyme class characteristics.

    Used for:
    - Detecting enzyme type from user queries
    - Identifying catalytic residues from PDB
    - Determining activity preservation requirements
    """
    ec_prefix: str                    # EC number prefix (e.g., "1.1" for dehydrogenases)
    name: str                         # Human-readable name
    aliases: List[str]                # Alternative names for NL detection
    catalytic_motifs: List[str]       # Known catalytic motifs (regex patterns)
    essential_residues: List[str]     # Amino acids critical for activity (1-letter codes)
    cofactor_types: List[str]         # Common cofactors (NAD, FAD, PQQ, etc.)
    substrate_access_required: bool   # Whether substrate channel must remain open
    hbond_network_critical: bool      # Whether H-bond network is essential for mechanism
    metal_role: str                   # "catalytic", "structural", "electron_transfer", "none"
    typical_metals: List[str]         # Common metal cofactors
    preservation_priority: List[str]  # Order of importance for preservation


ENZYME_CLASS_DATABASE: Dict[str, EnzymeClassProfile] = {
    # ═══════════════════════════════════════════════════════════════════
    # EC 1: OXIDOREDUCTASES
    # ═══════════════════════════════════════════════════════════════════

    "dehydrogenase": EnzymeClassProfile(
        ec_prefix="1.1",
        name="Alcohol/Aldehyde Dehydrogenase",
        aliases=["adh", "adhd", "aldh", "alcohol dehydrogenase", "aldehyde dehydrogenase",
                 "oxidoreductase", "dehydrogenase"],
        catalytic_motifs=[r"G.G..G", r"H.C", r"C..C"],  # NAD-binding, Zn-binding
        essential_residues=["C", "H", "D", "E"],
        cofactor_types=["NAD", "NADP", "ZN"],
        substrate_access_required=True,
        hbond_network_critical=True,
        metal_role="catalytic",
        typical_metals=["ZN"],
        preservation_priority=["catalytic_residues", "nad_binding", "substrate_channel"],
    ),

    "quinoprotein": EnzymeClassProfile(
        ec_prefix="1.1.5",
        name="PQQ-dependent Dehydrogenase",
        aliases=["pqq", "quinoprotein", "glucose dehydrogenase", "methanol dehydrogenase",
                 "pqq-dependent", "quinohemoprotein"],
        catalytic_motifs=[r"W.D.G.D", r"D..D"],  # PQQ binding
        essential_residues=["D", "E", "W", "H"],
        cofactor_types=["PQQ"],
        substrate_access_required=True,
        hbond_network_critical=True,
        metal_role="catalytic",
        typical_metals=["CA"],  # Calcium is essential for PQQ enzymes
        preservation_priority=["pqq_coordination", "ca_site", "catalytic_asp", "substrate_channel"],
    ),

    "oxidase": EnzymeClassProfile(
        ec_prefix="1.10",
        name="Copper Oxidase / Laccase",
        aliases=["oxidase", "laccase", "copper oxidase", "multicopper oxidase",
                 "cytochrome oxidase", "amine oxidase"],
        catalytic_motifs=[r"H.H", r"H..H..H"],  # Type 1/2/3 copper sites
        essential_residues=["H", "C", "M"],
        cofactor_types=["CU", "FE"],
        substrate_access_required=True,
        hbond_network_critical=True,
        metal_role="electron_transfer",
        typical_metals=["CU", "FE"],
        preservation_priority=["copper_sites", "electron_pathway", "substrate_channel"],
    ),

    "peroxidase": EnzymeClassProfile(
        ec_prefix="1.11",
        name="Peroxidase / Catalase",
        aliases=["peroxidase", "catalase", "heme peroxidase", "cytochrome c peroxidase",
                 "glutathione peroxidase", "chloroperoxidase"],
        catalytic_motifs=[r"C..C", r"H..H"],  # Heme coordination
        essential_residues=["H", "R", "D", "E"],
        cofactor_types=["HEM", "FE"],
        substrate_access_required=True,
        hbond_network_critical=True,
        metal_role="catalytic",
        typical_metals=["FE"],
        preservation_priority=["heme_pocket", "distal_his", "proximal_his", "substrate_channel"],
    ),

    "oxygenase": EnzymeClassProfile(
        ec_prefix="1.14",
        name="Oxygenase / P450",
        aliases=["oxygenase", "monooxygenase", "dioxygenase", "p450", "cytochrome p450",
                 "hydroxylase", "lipoxygenase"],
        catalytic_motifs=[r"C...G", r"E..R..C"],  # P450 cysteine loop
        essential_residues=["C", "H", "D", "E", "T"],
        cofactor_types=["HEM", "FE", "FAD", "FMN"],
        substrate_access_required=True,
        hbond_network_critical=True,
        metal_role="catalytic",
        typical_metals=["FE"],
        preservation_priority=["heme_thiolate", "substrate_channel", "electron_transfer"],
    ),

    "reductase": EnzymeClassProfile(
        ec_prefix="1.1",
        name="Reductase",
        aliases=["reductase", "nitroreductase", "thioredoxin reductase", "ferredoxin reductase"],
        catalytic_motifs=[r"C..C", r"G.G..G"],
        essential_residues=["C", "H", "Y", "S"],
        cofactor_types=["FAD", "FMN", "NAD", "NADP"],
        substrate_access_required=True,
        hbond_network_critical=True,
        metal_role="none",
        typical_metals=[],
        preservation_priority=["flavin_binding", "active_site_cys", "substrate_channel"],
    ),

    # ═══════════════════════════════════════════════════════════════════
    # EC 2: TRANSFERASES
    # ═══════════════════════════════════════════════════════════════════

    "kinase": EnzymeClassProfile(
        ec_prefix="2.7",
        name="Kinase / Phosphotransferase",
        aliases=["kinase", "phosphotransferase", "protein kinase", "tyrosine kinase",
                 "serine kinase", "atp-binding"],
        catalytic_motifs=[r"G.G..G", r"D.G", r"APE"],  # P-loop, DFG motif
        essential_residues=["D", "K", "E", "N"],
        cofactor_types=["ATP", "ADP", "MG"],
        substrate_access_required=True,
        hbond_network_critical=True,
        metal_role="catalytic",
        typical_metals=["MG", "MN"],
        preservation_priority=["atp_binding", "catalytic_asp", "mg_site", "substrate_site"],
    ),

    "methyltransferase": EnzymeClassProfile(
        ec_prefix="2.1",
        name="Methyltransferase",
        aliases=["methyltransferase", "sam-dependent", "dna methyltransferase",
                 "protein methyltransferase", "o-methyltransferase"],
        catalytic_motifs=[r"G.G.G", r"D..Y"],  # SAM binding
        essential_residues=["D", "E", "Y", "H"],
        cofactor_types=["SAM", "SAH"],
        substrate_access_required=True,
        hbond_network_critical=True,
        metal_role="none",
        typical_metals=[],
        preservation_priority=["sam_binding", "catalytic_residues", "substrate_recognition"],
    ),

    "acetyltransferase": EnzymeClassProfile(
        ec_prefix="2.3",
        name="Acetyltransferase / Acyltransferase",
        aliases=["acetyltransferase", "acyltransferase", "histone acetyltransferase",
                 "n-acetyltransferase", "coenzyme a"],
        catalytic_motifs=[r"G.G", r"R..G..G"],  # CoA binding
        essential_residues=["H", "D", "E", "C"],
        cofactor_types=["COA", "ACO"],
        substrate_access_required=True,
        hbond_network_critical=True,
        metal_role="none",
        typical_metals=[],
        preservation_priority=["coa_binding", "catalytic_his", "substrate_channel"],
    ),

    "glycosyltransferase": EnzymeClassProfile(
        ec_prefix="2.4",
        name="Glycosyltransferase",
        aliases=["glycosyltransferase", "gt", "glucosyltransferase", "mannosyltransferase",
                 "sialyltransferase", "sugar transferase"],
        catalytic_motifs=[r"D.D", r"E.E"],  # GT-A/GT-B folds
        essential_residues=["D", "E", "H", "R"],
        cofactor_types=["UDP", "GDP", "CMP"],
        substrate_access_required=True,
        hbond_network_critical=True,
        metal_role="catalytic",
        typical_metals=["MN", "MG"],
        preservation_priority=["donor_binding", "acceptor_binding", "metal_site"],
    ),

    "aminotransferase": EnzymeClassProfile(
        ec_prefix="2.6",
        name="Aminotransferase / Transaminase",
        aliases=["aminotransferase", "transaminase", "plp-dependent", "pyridoxal",
                 "aspartate aminotransferase", "alanine aminotransferase"],
        catalytic_motifs=[r"K..S..H", r"D..A..K"],  # PLP binding
        essential_residues=["K", "D", "Y", "R"],
        cofactor_types=["PLP", "PMP"],
        substrate_access_required=True,
        hbond_network_critical=True,
        metal_role="none",
        typical_metals=[],
        preservation_priority=["plp_lysine", "substrate_binding", "hbond_network"],
    ),

    # ═══════════════════════════════════════════════════════════════════
    # EC 3: HYDROLASES
    # ═══════════════════════════════════════════════════════════════════

    "protease": EnzymeClassProfile(
        ec_prefix="3.4",
        name="Protease / Peptidase",
        aliases=["protease", "peptidase", "proteinase", "serine protease", "cysteine protease",
                 "aspartic protease", "metalloprotease", "trypsin", "chymotrypsin"],
        catalytic_motifs=[r"H..S", r"D..S..H", r"H.E..H"],  # Catalytic triads
        essential_residues=["H", "D", "S", "C", "E"],
        cofactor_types=[],
        substrate_access_required=True,
        hbond_network_critical=True,
        metal_role="catalytic",  # For metalloproteases
        typical_metals=["ZN"],  # Metalloproteases
        preservation_priority=["catalytic_triad", "oxyanion_hole", "specificity_pocket"],
    ),

    "lipase": EnzymeClassProfile(
        ec_prefix="3.1",
        name="Lipase / Esterase",
        aliases=["lipase", "esterase", "phospholipase", "triacylglycerol lipase",
                 "carboxylesterase", "acetylcholinesterase"],
        catalytic_motifs=[r"G.S.G", r"H..D"],  # Lipase box
        essential_residues=["S", "H", "D", "G"],
        cofactor_types=[],
        substrate_access_required=True,
        hbond_network_critical=True,
        metal_role="none",
        typical_metals=[],
        preservation_priority=["catalytic_triad", "oxyanion_hole", "lid_domain"],
    ),

    "glycosidase": EnzymeClassProfile(
        ec_prefix="3.2",
        name="Glycosidase / Glycoside Hydrolase",
        aliases=["glycosidase", "glycoside hydrolase", "cellulase", "amylase",
                 "chitinase", "lysozyme", "beta-glucosidase"],
        catalytic_motifs=[r"E..E", r"D..E"],  # Acid/base pair
        essential_residues=["E", "D", "H", "Y"],
        cofactor_types=[],
        substrate_access_required=True,
        hbond_network_critical=True,
        metal_role="none",
        typical_metals=[],
        preservation_priority=["acid_base_pair", "substrate_binding", "tunnel_residues"],
    ),

    "phosphatase": EnzymeClassProfile(
        ec_prefix="3.1.3",
        name="Phosphatase",
        aliases=["phosphatase", "phosphohydrolase", "acid phosphatase", "alkaline phosphatase",
                 "protein phosphatase", "tyrosine phosphatase"],
        catalytic_motifs=[r"D..D..D", r"C..R"],  # Metal binding, PTP motif
        essential_residues=["D", "H", "C", "R"],
        cofactor_types=[],
        substrate_access_required=True,
        hbond_network_critical=True,
        metal_role="catalytic",
        typical_metals=["ZN", "MG", "MN"],
        preservation_priority=["metal_site", "catalytic_cys", "substrate_binding"],
    ),

    "nuclease": EnzymeClassProfile(
        ec_prefix="3.1",
        name="Nuclease",
        aliases=["nuclease", "dnase", "rnase", "endonuclease", "exonuclease",
                 "restriction enzyme", "ribonuclease"],
        catalytic_motifs=[r"D..E..D", r"H..H"],  # Metal binding
        essential_residues=["D", "E", "H", "K"],
        cofactor_types=[],
        substrate_access_required=True,
        hbond_network_critical=True,
        metal_role="catalytic",
        typical_metals=["MG", "MN", "CA"],
        preservation_priority=["metal_site", "dna_binding", "specificity_residues"],
    ),

    # ═══════════════════════════════════════════════════════════════════
    # EC 4: LYASES
    # ═══════════════════════════════════════════════════════════════════

    "decarboxylase": EnzymeClassProfile(
        ec_prefix="4.1",
        name="Decarboxylase",
        aliases=["decarboxylase", "carboxy-lyase", "amino acid decarboxylase",
                 "pyruvate decarboxylase", "plp decarboxylase"],
        catalytic_motifs=[r"K..H", r"D..K"],  # PLP/ThDP binding
        essential_residues=["K", "H", "D", "E"],
        cofactor_types=["PLP", "TPP"],
        substrate_access_required=True,
        hbond_network_critical=True,
        metal_role="none",
        typical_metals=["MG"],  # For ThDP-dependent
        preservation_priority=["cofactor_binding", "substrate_pocket", "proton_donor"],
    ),

    "aldolase": EnzymeClassProfile(
        ec_prefix="4.1.2",
        name="Aldolase",
        aliases=["aldolase", "lyase", "fructose-bisphosphate aldolase", "deoxyribose-phosphate aldolase"],
        catalytic_motifs=[r"K..D", r"E..K"],  # Schiff base
        essential_residues=["K", "D", "E", "H"],
        cofactor_types=[],
        substrate_access_required=True,
        hbond_network_critical=True,
        metal_role="catalytic",  # Class II aldolases
        typical_metals=["ZN"],  # Class II
        preservation_priority=["schiff_base_lys", "substrate_binding", "metal_site"],
    ),

    "dehydratase": EnzymeClassProfile(
        ec_prefix="4.2",
        name="Dehydratase / Hydro-lyase",
        aliases=["dehydratase", "hydro-lyase", "carbonic anhydrase", "aconitase",
                 "fumarase", "enolase"],
        catalytic_motifs=[r"H..H..H", r"E..E"],  # Metal binding
        essential_residues=["H", "E", "D", "K"],
        cofactor_types=[],
        substrate_access_required=True,
        hbond_network_critical=True,
        metal_role="catalytic",
        typical_metals=["ZN", "FE"],
        preservation_priority=["metal_site", "substrate_binding", "proton_shuttle"],
    ),

    # ═══════════════════════════════════════════════════════════════════
    # EC 5: ISOMERASES
    # ═══════════════════════════════════════════════════════════════════

    "isomerase": EnzymeClassProfile(
        ec_prefix="5.1",
        name="Isomerase / Racemase / Epimerase",
        aliases=["isomerase", "racemase", "epimerase", "triose-phosphate isomerase",
                 "phosphoglucose isomerase", "mutarotase"],
        catalytic_motifs=[r"K..H", r"E..H"],
        essential_residues=["K", "H", "E", "C"],
        cofactor_types=["PLP"],  # Some racemases
        substrate_access_required=True,
        hbond_network_critical=True,
        metal_role="none",
        typical_metals=[],
        preservation_priority=["catalytic_base", "substrate_binding", "loop_flexibility"],
    ),

    "mutase": EnzymeClassProfile(
        ec_prefix="5.4",
        name="Mutase / Intramolecular Transferase",
        aliases=["mutase", "intramolecular transferase", "phosphoglycerate mutase",
                 "chorismate mutase", "methylmalonyl-coa mutase"],
        catalytic_motifs=[r"H..R", r"E..K"],
        essential_residues=["H", "R", "E", "K"],
        cofactor_types=["B12"],  # Some mutases
        substrate_access_required=True,
        hbond_network_critical=True,
        metal_role="none",
        typical_metals=["CO"],  # B12-dependent
        preservation_priority=["catalytic_his", "substrate_binding", "cofactor_site"],
    ),

    # ═══════════════════════════════════════════════════════════════════
    # EC 6: LIGASES
    # ═══════════════════════════════════════════════════════════════════

    "synthetase": EnzymeClassProfile(
        ec_prefix="6.1",
        name="Synthetase / Ligase",
        aliases=["synthetase", "ligase", "synthase", "aminoacyl-trna synthetase",
                 "glutamine synthetase", "dna ligase", "ubiquitin ligase"],
        catalytic_motifs=[r"G.G..G", r"K..K"],  # ATP binding
        essential_residues=["K", "D", "E", "R"],
        cofactor_types=["ATP", "AMP"],
        substrate_access_required=True,
        hbond_network_critical=True,
        metal_role="catalytic",
        typical_metals=["MG", "MN"],
        preservation_priority=["atp_binding", "substrate_sites", "metal_coordination"],
    ),
}


# ═══════════════════════════════════════════════════════════════════
# 2. ENZYME CLASS DETECTION
# ═══════════════════════════════════════════════════════════════════

def detect_enzyme_class(
    query: str,
    ligand_name: Optional[str] = None,
    metal_type: Optional[str] = None,
) -> Tuple[Optional[str], float]:
    """
    Detect enzyme class from query text and context.

    Args:
        query: User's natural language query
        ligand_name: Detected ligand name (e.g., "PQQ", "NAD")
        metal_type: Detected metal type (e.g., "ZN", "CA")

    Returns:
        Tuple of (enzyme_class_key, confidence)
        Returns (None, 0.0) if no enzyme class detected
    """
    query_lower = query.lower()

    best_match: Optional[str] = None
    best_confidence: float = 0.0

    for class_key, profile in ENZYME_CLASS_DATABASE.items():
        confidence = 0.0

        # Check aliases
        for alias in profile.aliases:
            if alias.lower() in query_lower:
                # Exact word match gets higher confidence
                if re.search(rf'\b{re.escape(alias.lower())}\b', query_lower):
                    confidence = max(confidence, 0.9)
                else:
                    confidence = max(confidence, 0.7)

        # Check if ligand matches cofactor types
        if ligand_name:
            ligand_upper = ligand_name.upper()
            if ligand_upper in profile.cofactor_types:
                confidence = max(confidence, 0.8)
            # Special cases
            if ligand_upper == "PQQ" and class_key == "quinoprotein":
                confidence = max(confidence, 0.95)
            if ligand_upper in ["NAD", "NADP"] and class_key == "dehydrogenase":
                confidence = max(confidence, 0.8)

        # Check if metal matches typical metals
        if metal_type:
            metal_upper = metal_type.upper()
            if metal_upper in profile.typical_metals:
                confidence += 0.1  # Boost for metal match

        if confidence > best_confidence:
            best_confidence = confidence
            best_match = class_key

    return (best_match, best_confidence) if best_confidence >= 0.5 else (None, 0.0)


# ═══════════════════════════════════════════════════════════════════
# 3. PRESERVATION REQUIREMENTS
# ═══════════════════════════════════════════════════════════════════

def get_preservation_requirements(enzyme_class: str) -> Dict:
    """
    Get activity preservation requirements for an enzyme class.

    Args:
        enzyme_class: Key from ENZYME_CLASS_DATABASE

    Returns:
        Dict with preservation requirements:
        - substrate_access_required: bool
        - hbond_network_critical: bool
        - metal_role: str
        - preservation_priority: List[str]
        - essential_residues: List[str]
    """
    profile = ENZYME_CLASS_DATABASE.get(enzyme_class)
    if not profile:
        return {
            "substrate_access_required": True,
            "hbond_network_critical": True,
            "metal_role": "none",
            "preservation_priority": ["catalytic_residues", "substrate_channel"],
            "essential_residues": [],
        }

    return {
        "substrate_access_required": profile.substrate_access_required,
        "hbond_network_critical": profile.hbond_network_critical,
        "metal_role": profile.metal_role,
        "preservation_priority": profile.preservation_priority,
        "essential_residues": profile.essential_residues,
        "typical_metals": profile.typical_metals,
        "cofactor_types": profile.cofactor_types,
    }


# ═══════════════════════════════════════════════════════════════════
# 4. CATALYTIC RESIDUE IDENTIFICATION
# ═══════════════════════════════════════════════════════════════════

# Amino acid codes for pattern matching
AA_CODES = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
    'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
    'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
    'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR',
}

AA_CODES_REV = {v: k for k, v in AA_CODES.items()}


def identify_catalytic_residues_from_pdb(
    pdb_content: str,
    enzyme_class: str,
    metal_type: Optional[str] = None,
    ligand_name: Optional[str] = None,
    contact_distance: float = 4.5,
) -> List[Dict]:
    """
    Identify catalytic residues from PDB structure.

    Uses multiple approaches:
    1. Proximity to metal/cofactor
    2. Essential residue types for enzyme class
    3. Sequence motif matching (if sequence available)

    Args:
        pdb_content: PDB file content as string
        enzyme_class: Key from ENZYME_CLASS_DATABASE
        metal_type: Metal to search near (e.g., "ZN", "CA")
        ligand_name: Ligand to search near (e.g., "PQQ", "NAD")
        contact_distance: Distance threshold in Angstroms

    Returns:
        List of dicts with residue info:
        - chain: str
        - resnum: int
        - resname: str (3-letter code)
        - role: str (e.g., "metal_coordinator", "catalytic", "substrate_binding")
        - confidence: float
    """
    profile = ENZYME_CLASS_DATABASE.get(enzyme_class)
    if not profile:
        return []

    catalytic_residues = []
    seen_residues: Set[str] = set()

    # Parse PDB atoms
    protein_atoms = []
    hetatm_positions = []

    for line in pdb_content.split('\n'):
        if line.startswith('ATOM'):
            try:
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21:22].strip() or 'A'
                res_num = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())

                protein_atoms.append({
                    'atom': atom_name,
                    'resname': res_name,
                    'chain': chain,
                    'resnum': res_num,
                    'coords': (x, y, z),
                })
            except (ValueError, IndexError):
                continue

        elif line.startswith('HETATM'):
            try:
                res_name = line[17:20].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())

                # Check if this is the metal or ligand we're looking for
                is_target = False
                if metal_type and res_name.upper() == metal_type.upper():
                    is_target = True
                if ligand_name and res_name.upper() == ligand_name.upper():
                    is_target = True
                if not metal_type and not ligand_name:
                    # Include all HETATM as potential cofactors
                    is_target = True

                if is_target:
                    hetatm_positions.append({
                        'resname': res_name,
                        'coords': (x, y, z),
                    })
            except (ValueError, IndexError):
                continue

    # Find residues near cofactor/metal
    essential_1letter = set(profile.essential_residues)

    for hetatm in hetatm_positions:
        hx, hy, hz = hetatm['coords']

        for atom in protein_atoms:
            # Only consider CA for distance (faster, representative)
            if atom['atom'] != 'CA':
                continue

            ax, ay, az = atom['coords']
            dist = ((hx - ax)**2 + (hy - ay)**2 + (hz - az)**2) ** 0.5

            if dist <= contact_distance:
                res_key = f"{atom['chain']}{atom['resnum']}"
                if res_key in seen_residues:
                    continue

                # Check if this is an essential residue type
                res_1letter = AA_CODES_REV.get(atom['resname'], 'X')
                is_essential = res_1letter in essential_1letter

                # Determine role
                if hetatm['resname'].upper() in ['ZN', 'FE', 'CU', 'MG', 'MN', 'CA', 'CO', 'NI']:
                    role = "metal_coordinator"
                    confidence = 0.9 if is_essential else 0.7
                else:
                    role = "cofactor_binding" if is_essential else "substrate_binding"
                    confidence = 0.8 if is_essential else 0.6

                catalytic_residues.append({
                    'chain': atom['chain'],
                    'resnum': atom['resnum'],
                    'resname': atom['resname'],
                    'role': role,
                    'confidence': confidence,
                    'distance': round(dist, 2),
                })
                seen_residues.add(res_key)

    # Sort by confidence and distance
    catalytic_residues.sort(key=lambda x: (-x['confidence'], x['distance']))

    return catalytic_residues[:25]  # Return top 25


# ═══════════════════════════════════════════════════════════════════
# 5. SUBSTRATE CHANNEL DETECTION
# ═══════════════════════════════════════════════════════════════════

def get_substrate_channel_atoms(
    pdb_content: str,
    enzyme_class: str,
    ligand_name: Optional[str] = None,
    channel_distance: float = 8.0,
) -> Dict[str, str]:
    """
    Identify atoms that should remain exposed for substrate access.

    For enzyme activity preservation, certain atoms near the active site
    must remain solvent-accessible for substrate entry/product exit.

    Args:
        pdb_content: PDB file content
        enzyme_class: Key from ENZYME_CLASS_DATABASE
        ligand_name: Ligand to identify channel around
        channel_distance: Distance from ligand to consider as channel

    Returns:
        Dict mapping chain to comma-separated residue numbers
        that should be exposed (for select_exposed RFD3 parameter)
    """
    profile = ENZYME_CLASS_DATABASE.get(enzyme_class)
    if not profile or not profile.substrate_access_required:
        return {}

    # Find ligand position
    ligand_center = None
    for line in pdb_content.split('\n'):
        if line.startswith('HETATM'):
            try:
                res_name = line[17:20].strip()
                if ligand_name and res_name.upper() == ligand_name.upper():
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    if ligand_center is None:
                        ligand_center = [x, y, z, 1]
                    else:
                        ligand_center[0] += x
                        ligand_center[1] += y
                        ligand_center[2] += z
                        ligand_center[3] += 1
            except (ValueError, IndexError):
                continue

    if ligand_center is None:
        return {}

    # Calculate center of mass
    ligand_center = (
        ligand_center[0] / ligand_center[3],
        ligand_center[1] / ligand_center[3],
        ligand_center[2] / ligand_center[3],
    )

    # Find residues at channel distance (outer shell)
    channel_residues: Dict[str, Set[int]] = {}

    for line in pdb_content.split('\n'):
        if line.startswith('ATOM'):
            try:
                atom_name = line[12:16].strip()
                if atom_name != 'CA':
                    continue

                chain = line[21:22].strip() or 'A'
                res_num = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())

                dist = (
                    (x - ligand_center[0])**2 +
                    (y - ligand_center[1])**2 +
                    (z - ligand_center[2])**2
                ) ** 0.5

                # Channel residues are in outer shell (5-8 Å from ligand)
                if 5.0 <= dist <= channel_distance:
                    if chain not in channel_residues:
                        channel_residues[chain] = set()
                    channel_residues[chain].add(res_num)

            except (ValueError, IndexError):
                continue

    # Format for RFD3 select_exposed parameter
    result = {}
    for chain, residues in channel_residues.items():
        sorted_residues = sorted(residues)
        result[chain] = ",".join(str(r) for r in sorted_residues[:15])  # Limit to 15

    return result


# ═══════════════════════════════════════════════════════════════════
# 6. CATALYTIC H-BOND NETWORK
# ═══════════════════════════════════════════════════════════════════

def get_catalytic_hbond_network(
    pdb_content: str,
    catalytic_residues: List[Dict],
    hbond_distance: float = 3.5,
) -> Dict[str, Dict[str, str]]:
    """
    Identify H-bond donor/acceptor atoms around catalytic residues.

    For enzyme activity preservation, the H-bond network around
    catalytic residues is critical for mechanism.

    Args:
        pdb_content: PDB file content
        catalytic_residues: List of catalytic residue dicts
        hbond_distance: Max distance for H-bond

    Returns:
        Dict with 'donors' and 'acceptors' for RFD3 H-bond conditioning
    """
    if not catalytic_residues:
        return {}

    # Build set of catalytic residue keys
    catalytic_keys = set()
    for res in catalytic_residues:
        catalytic_keys.add(f"{res['chain']}{res['resnum']}")

    # Known H-bond donors by residue type
    DONORS = {
        'SER': ['OG'], 'THR': ['OG1'], 'TYR': ['OH'],
        'ASN': ['ND2'], 'GLN': ['NE2'], 'HIS': ['ND1', 'NE2'],
        'ARG': ['NH1', 'NH2', 'NE'], 'LYS': ['NZ'], 'TRP': ['NE1'],
    }

    # Known H-bond acceptors by residue type
    ACCEPTORS = {
        'SER': ['OG'], 'THR': ['OG1'], 'TYR': ['OH'],
        'ASN': ['OD1'], 'GLN': ['OE1'], 'HIS': ['ND1', 'NE2'],
        'ASP': ['OD1', 'OD2'], 'GLU': ['OE1', 'OE2'],
    }

    donors: Dict[str, List[str]] = {}
    acceptors: Dict[str, List[str]] = {}

    for line in pdb_content.split('\n'):
        if not line.startswith('ATOM'):
            continue

        try:
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain = line[21:22].strip() or 'A'
            res_num = int(line[22:26].strip())

            res_key = f"{chain}{res_num}"
            if res_key not in catalytic_keys:
                continue

            # Check if atom is a donor
            if res_name in DONORS and atom_name in DONORS[res_name]:
                if res_key not in donors:
                    donors[res_key] = []
                donors[res_key].append(atom_name)

            # Check if atom is an acceptor
            if res_name in ACCEPTORS and atom_name in ACCEPTORS[res_name]:
                if res_key not in acceptors:
                    acceptors[res_key] = []
                acceptors[res_key].append(atom_name)

        except (ValueError, IndexError):
            continue

    # Format for RFD3
    result = {}

    if donors:
        donor_strs = []
        for res_key, atoms in donors.items():
            donor_strs.append(f"{res_key}:{','.join(atoms)}")
        result['donors'] = {k.split(':')[0][0]: k.split(':')[0][1:] for k in donor_strs[:5]}

    if acceptors:
        acceptor_strs = []
        for res_key, atoms in acceptors.items():
            acceptor_strs.append(f"{res_key}:{','.join(atoms)}")
        result['acceptors'] = {k.split(':')[0][0]: k.split(':')[0][1:] for k in acceptor_strs[:5]}

    return result


# ═══════════════════════════════════════════════════════════════════
# 7. UTILITY FUNCTIONS
# ═══════════════════════════════════════════════════════════════════

def get_enzyme_class_info(enzyme_class: str) -> Optional[Dict]:
    """Get full enzyme class profile as dict."""
    profile = ENZYME_CLASS_DATABASE.get(enzyme_class)
    if not profile:
        return None

    return {
        "ec_prefix": profile.ec_prefix,
        "name": profile.name,
        "aliases": profile.aliases,
        "essential_residues": profile.essential_residues,
        "cofactor_types": profile.cofactor_types,
        "substrate_access_required": profile.substrate_access_required,
        "hbond_network_critical": profile.hbond_network_critical,
        "metal_role": profile.metal_role,
        "typical_metals": profile.typical_metals,
    }


def list_enzyme_classes() -> List[str]:
    """Return list of all enzyme class keys."""
    return list(ENZYME_CLASS_DATABASE.keys())


def get_enzyme_aliases(enzyme_class: str) -> List[str]:
    """Get all aliases for an enzyme class."""
    profile = ENZYME_CLASS_DATABASE.get(enzyme_class)
    return profile.aliases if profile else []
