/**
 * Ligand Contact Analysis for MolStar
 * Extracted and adapted from Protein_engineering_tools/src/components/ProteinViewer.tsx
 *
 * Provides:
 * - Ligand detection (non-polymer, non-water, non-metal molecules)
 * - Protein-ligand contact identification
 * - Interaction type classification (H-bond, hydrophobic, salt bridge, pi-stacking)
 * - Binding site classification (functional vs crystal artifact)
 */

import { METAL_ELEMENTS } from './metalAnalysis';

// Types
export type InteractionType = 'hydrogen_bond' | 'hydrophobic' | 'salt_bridge' | 'pi_stacking' | 'other';

export type PharmacophoreType = 'donor' | 'acceptor' | 'aromatic' | 'hydrophobic' | 'positive' | 'negative';

export type QualityRating = 'EXCELLENT' | 'GOOD' | 'MODERATE' | 'POOR';

export interface BinderQualityResult {
  score: number;      // 0-100
  rating: QualityRating;
  breakdown: {
    hbondScore: number;       // H-bond contribution (max 30)
    hydrophobicScore: number; // Hydrophobic contribution (max 25)
    piStackScore: number;     // Pi-stacking contribution (max 15)
    saltBridgeScore: number;  // Salt bridge contribution (max 15)
    burialnessScore: number;  // Contact density (max 15)
  };
  suggestions: string[];
}

// Shell analysis cutoffs
const PRIMARY_SHELL_CUTOFF = 4.0;   // Direct contacts
const SECONDARY_SHELL_CUTOFF = 6.0;  // Extended environment

export interface ShellAnalysis {
  primary: LigandContact[];   // < 4A: direct binding
  secondary: LigandContact[]; // 4-6A: extended shell
  stats: {
    primaryCount: number;
    secondaryCount: number;
    totalCount: number;
    avgPrimaryDistance: number;
    avgSecondaryDistance: number;
  };
}

export interface PharmacophoreFeature {
  type: PharmacophoreType;
  residue: string;
  chain: string;
  atom: string;
  distance: number;
  position?: [number, number, number];
}

// Pi-stacking geometry types
export interface RingGeometry {
  centroid: [number, number, number];
  normal: [number, number, number];
}

export interface PiStackingResult {
  type: 'parallel' | 't-shaped' | 'offset-parallel' | 'none';
  distance: number;
  angle: number;
  offset: number;
  isStacking: boolean;
  ligandRingIdx?: number;
  proteinResidue?: string;
  proteinChain?: string;
}

// Aromatic ring atom patterns for known residues
export const AROMATIC_PATTERNS: Record<string, string[][]> = {
  'PHE': [['CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2']],
  'TYR': [['CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2']],
  'TRP': [['CG', 'CD1', 'NE1', 'CE2', 'CD2'], ['CE2', 'CD2', 'CE3', 'CZ3', 'CH2', 'CZ2']],
  'HIS': [['CG', 'ND1', 'CE1', 'NE2', 'CD2']],
};

/**
 * Calculate ring geometry (centroid and normal vector) from atom positions
 * Requires at least 5 atoms to form an aromatic ring
 */
export function calculateRingGeometry(
  positions: [number, number, number][]
): RingGeometry | null {
  if (positions.length < 5) return null;

  // Calculate centroid
  const centroid: [number, number, number] = [0, 0, 0];
  for (const pos of positions) {
    centroid[0] += pos[0];
    centroid[1] += pos[1];
    centroid[2] += pos[2];
  }
  centroid[0] /= positions.length;
  centroid[1] /= positions.length;
  centroid[2] /= positions.length;

  // Calculate normal using first 3 points (cross product of two edge vectors)
  const v1: [number, number, number] = [
    positions[1][0] - positions[0][0],
    positions[1][1] - positions[0][1],
    positions[1][2] - positions[0][2],
  ];
  const v2: [number, number, number] = [
    positions[2][0] - positions[0][0],
    positions[2][1] - positions[0][1],
    positions[2][2] - positions[0][2],
  ];

  // Cross product
  const normal: [number, number, number] = [
    v1[1] * v2[2] - v1[2] * v2[1],
    v1[2] * v2[0] - v1[0] * v2[2],
    v1[0] * v2[1] - v1[1] * v2[0],
  ];

  // Normalize
  const length = Math.sqrt(normal[0] ** 2 + normal[1] ** 2 + normal[2] ** 2);
  if (length < 0.0001) return null;

  normal[0] /= length;
  normal[1] /= length;
  normal[2] /= length;

  return { centroid, normal };
}

/**
 * Analyze pi-stacking interaction between two aromatic rings
 * Returns stacking type (parallel, t-shaped, offset-parallel, or none)
 */
export function analyzePiStacking(
  ring1: RingGeometry,
  ring2: RingGeometry
): PiStackingResult {
  // Distance between centroids
  const dx = ring2.centroid[0] - ring1.centroid[0];
  const dy = ring2.centroid[1] - ring1.centroid[1];
  const dz = ring2.centroid[2] - ring1.centroid[2];
  const distance = Math.sqrt(dx * dx + dy * dy + dz * dz);

  // Angle between normals (dot product)
  const dotProduct = Math.abs(
    ring1.normal[0] * ring2.normal[0] +
    ring1.normal[1] * ring2.normal[1] +
    ring1.normal[2] * ring2.normal[2]
  );
  const angle = Math.acos(Math.min(1, dotProduct)) * 180 / Math.PI;

  // Calculate offset (horizontal displacement from stacking axis)
  const projection =
    dx * ring1.normal[0] + dy * ring1.normal[1] + dz * ring1.normal[2];
  const offset = Math.sqrt(Math.max(0, distance * distance - projection * projection));

  // Classify stacking type
  let type: PiStackingResult['type'] = 'none';

  if (distance <= 4.5) {
    if (angle < 30) {
      // Rings are roughly parallel
      type = offset < 2.0 ? 'parallel' : 'offset-parallel';
    } else if (angle > 60) {
      // Rings are roughly perpendicular
      type = 't-shaped';
    }
  }

  return {
    type,
    distance,
    angle,
    offset,
    isStacking: type !== 'none',
  };
}

export interface LigandContact {
  residue: string;
  chain: string;
  atom: string;
  distance: number;
  interactionType: InteractionType;
}

export interface LigandData {
  name: string;
  info: string;
  atoms: number;
  chainId: string;
  resSeq: number;
  contacts: LigandContact[];
  proteinContactCount: number;
  waterContactCount: number;
  bindingSiteType: 'functional' | 'crystal_artifact' | 'uncertain';
  bindingSiteReason: string;
  piStackingResults?: PiStackingResult[];
}

export interface LigandAnalysisResult {
  ligandCount: number;
  ligandDetails: LigandData[];
}

// Common crystallization additives and buffers to flag
export const CRYSTAL_ADDITIVES = new Set([
  'SO4', 'PO4', 'GOL', 'EDO', 'PEG', 'MPD', 'DMS', 'ACT', 'CIT', 'TRS',
  'BME', 'MES', 'EPE', 'IMD', 'SCN', 'NO3', 'CL', 'BR', 'IOD', 'F',
  'BU3', 'PG4', '1PE', 'P6G', 'PGE', 'ARS'
]);

// H-bond capable atoms
const H_BOND_DONORS = new Set(['N', 'O', 'S']);
const H_BOND_ACCEPTORS = new Set(['N', 'O', 'S', 'F']);

// Hydrophobic atoms
const HYDROPHOBIC_ATOMS = new Set(['C']);

// Charged residues for salt bridges
const POSITIVE_RESIDUES = new Set(['ARG', 'LYS', 'HIS']);
const NEGATIVE_RESIDUES = new Set(['ASP', 'GLU']);

// Standard amino acids
const STANDARD_AMINO_ACIDS = new Set([
  'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
  'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
]);

// Distance calculation helper
function distance(
  a: [number, number, number],
  b: [number, number, number]
): number {
  const dx = a[0] - b[0];
  const dy = a[1] - b[1];
  const dz = a[2] - b[2];
  return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

/**
 * Determine interaction type based on atom properties and distance
 */
function determineInteractionType(
  atomElement: string,
  resName: string,
  dist: number
): InteractionType {
  // H-bond: donor/acceptor atoms within 3.5Å
  if (dist <= 3.5 && (H_BOND_DONORS.has(atomElement) || H_BOND_ACCEPTORS.has(atomElement))) {
    return 'hydrogen_bond';
  }

  // Hydrophobic: carbon atoms within 4.0Å
  if (dist <= 4.0 && HYDROPHOBIC_ATOMS.has(atomElement)) {
    return 'hydrophobic';
  }

  // Salt bridge: charged residues within 4.0Å
  if (dist <= 4.0 && (POSITIVE_RESIDUES.has(resName) || NEGATIVE_RESIDUES.has(resName))) {
    return 'salt_bridge';
  }

  // Aromatic residues for pi-stacking
  const aromaticResidues = new Set(['PHE', 'TYR', 'TRP', 'HIS']);
  if (dist <= 4.5 && aromaticResidues.has(resName)) {
    return 'pi_stacking';
  }

  return 'other';
}

/**
 * Classify binding site as functional vs crystal artifact
 */
function classifyLigandBindingSite(
  ligandName: string,
  uniqueProteinResidues: number
): { type: 'functional' | 'crystal_artifact' | 'uncertain'; reason: string } {
  const isCrystalAdditive = CRYSTAL_ADDITIVES.has(ligandName);

  if (isCrystalAdditive && uniqueProteinResidues < 3) {
    return {
      type: 'crystal_artifact',
      reason: `Known crystallization additive (${ligandName}) with few contacts`
    };
  } else if (uniqueProteinResidues >= 5) {
    return {
      type: 'functional',
      reason: `${uniqueProteinResidues} protein residue contacts - well-defined binding pocket`
    };
  } else if (uniqueProteinResidues >= 3) {
    return {
      type: isCrystalAdditive ? 'uncertain' : 'functional',
      reason: `${uniqueProteinResidues} protein contacts${isCrystalAdditive ? ' (common additive - verify)' : ''}`
    };
  } else if (uniqueProteinResidues <= 1) {
    return {
      type: 'crystal_artifact',
      reason: `Minimal protein contacts (${uniqueProteinResidues}) - likely surface-bound additive`
    };
  } else {
    return {
      type: 'uncertain',
      reason: `${uniqueProteinResidues} protein contacts - needs manual inspection`
    };
  }
}

/**
 * Parse PDB content and find ligands with their protein contacts
 * This is a standalone function that works without Mol* structure
 */
export function findLigandContactsFromPDB(
  pdbContent: string,
  contactRadius: number = 4.0
): LigandAnalysisResult {
  const lines = pdbContent.split('\n');

  // Parse all atoms
  interface ParsedAtom {
    serial: number;
    name: string;
    resName: string;
    chainId: string;
    resSeq: number;
    x: number;
    y: number;
    z: number;
    element: string;
    isHETATM: boolean;
  }

  const atoms: ParsedAtom[] = [];

  for (const line of lines) {
    if (line.startsWith('ATOM') || line.startsWith('HETATM')) {
      const isHETATM = line.startsWith('HETATM');
      try {
        atoms.push({
          serial: parseInt(line.slice(6, 11).trim()),
          name: line.slice(12, 16).trim(),
          resName: line.slice(17, 20).trim(),
          chainId: line.slice(21, 22).trim(),
          resSeq: parseInt(line.slice(22, 26).trim()),
          x: parseFloat(line.slice(30, 38).trim()),
          y: parseFloat(line.slice(38, 46).trim()),
          z: parseFloat(line.slice(46, 54).trim()),
          element: line.slice(76, 78).trim().toUpperCase() || line.slice(12, 14).trim().replace(/[0-9]/g, '').toUpperCase(),
          isHETATM
        });
      } catch {
        // Skip malformed lines
      }
    }
  }

  // Identify ligand residues (HETATM, not water, not metal)
  const ligandResidues = new Map<string, {
    name: string;
    chainId: string;
    resSeq: number;
    atoms: ParsedAtom[];
  }>();

  for (const atom of atoms) {
    // Check if this is a ligand (HETATM, not water, not metal, not standard amino acid)
    const isMetal = METAL_ELEMENTS.has(atom.element) || METAL_ELEMENTS.has(atom.resName);
    const isWater = atom.resName === 'HOH' || atom.resName === 'WAT';
    const isStandardAA = STANDARD_AMINO_ACIDS.has(atom.resName);

    if (atom.isHETATM && !isWater && !isMetal && !isStandardAA) {
      const key = `${atom.chainId}:${atom.resSeq}:${atom.resName}`;

      if (!ligandResidues.has(key)) {
        ligandResidues.set(key, {
          name: atom.resName,
          chainId: atom.chainId,
          resSeq: atom.resSeq,
          atoms: []
        });
      }
      ligandResidues.get(key)!.atoms.push(atom);
    }
  }

  // Identify protein atoms
  const proteinAtoms = atoms.filter(a =>
    !a.isHETATM || STANDARD_AMINO_ACIDS.has(a.resName)
  );

  // Identify water atoms
  const waterAtoms = atoms.filter(a =>
    a.resName === 'HOH' || a.resName === 'WAT'
  );

  const ligandDetails: LigandData[] = [];

  // Find contacts for each ligand
  for (const [ligKey, ligand] of ligandResidues) {
    const contacts: LigandContact[] = [];
    const contactedResidues = new Set<string>();
    let proteinContactCount = 0;
    let waterContactCount = 0;

    // Get ligand atom positions
    const ligandPositions: [number, number, number][] = ligand.atoms.map(a => [a.x, a.y, a.z]);

    // Check protein atoms for contacts
    for (const proteinAtom of proteinAtoms) {
      const proteinPos: [number, number, number] = [proteinAtom.x, proteinAtom.y, proteinAtom.z];

      // Find minimum distance to any ligand atom
      let minDist = Infinity;
      for (const ligPos of ligandPositions) {
        const dist = distance(proteinPos, ligPos);
        if (dist < minDist) minDist = dist;
      }

      if (minDist <= contactRadius) {
        proteinContactCount++;
        const resKey = `${proteinAtom.chainId}:${proteinAtom.resSeq}:${proteinAtom.resName}`;

        // Only add unique residue contacts
        if (!contactedResidues.has(resKey)) {
          contactedResidues.add(resKey);

          const interactionType = determineInteractionType(
            proteinAtom.element,
            proteinAtom.resName,
            minDist
          );

          contacts.push({
            residue: `${proteinAtom.resName}${proteinAtom.resSeq}`,
            chain: proteinAtom.chainId,
            atom: proteinAtom.name,
            distance: minDist,
            interactionType
          });
        }
      }
    }

    // Check water atoms for contacts
    for (const waterAtom of waterAtoms) {
      const waterPos: [number, number, number] = [waterAtom.x, waterAtom.y, waterAtom.z];

      let minDist = Infinity;
      for (const ligPos of ligandPositions) {
        const dist = distance(waterPos, ligPos);
        if (dist < minDist) minDist = dist;
      }

      if (minDist <= contactRadius) {
        waterContactCount++;
      }
    }

    // Sort contacts by distance
    contacts.sort((a, b) => a.distance - b.distance);

    // Classify binding site
    const classification = classifyLigandBindingSite(ligand.name, contacts.length);

    ligandDetails.push({
      name: ligand.name,
      info: `${ligand.name} (Chain ${ligand.chainId}, ${ligand.resSeq})`,
      atoms: ligand.atoms.length,
      chainId: ligand.chainId,
      resSeq: ligand.resSeq,
      contacts,
      proteinContactCount,
      waterContactCount,
      bindingSiteType: classification.type,
      bindingSiteReason: classification.reason
    });
  }

  return {
    ligandCount: ligandDetails.length,
    ligandDetails
  };
}

/**
 * Get interaction type color for visualization
 */
export function getInteractionColor(type: InteractionType): string {
  switch (type) {
    case 'hydrogen_bond':
      return '#3B82F6'; // blue
    case 'hydrophobic':
      return '#22C55E'; // green
    case 'salt_bridge':
      return '#EF4444'; // red
    case 'pi_stacking':
      return '#A855F7'; // purple
    case 'other':
    default:
      return '#6B7280'; // gray
  }
}

/**
 * Get interaction type label for display
 */
export function getInteractionLabel(type: InteractionType): string {
  switch (type) {
    case 'hydrogen_bond':
      return 'H-bond';
    case 'hydrophobic':
      return 'Hydrophobic';
    case 'salt_bridge':
      return 'Salt bridge';
    case 'pi_stacking':
      return 'π-stacking';
    case 'other':
    default:
      return 'Contact';
  }
}

// Aromatic residues for pharmacophore classification
const AROMATIC_RESIDUES = new Set(['PHE', 'TYR', 'TRP', 'HIS']);

/**
 * Extract residue name from a contact residue string (e.g., "PHE101" -> "PHE")
 */
function extractResidueName(residue: string): string {
  // Remove numbers from the residue string to get the residue name
  return residue.replace(/[0-9]/g, '');
}

/**
 * Get the element from an atom name (e.g., "NE2" -> "N", "OD1" -> "O", "CD1" -> "C")
 */
function getAtomElement(atomName: string): string {
  // Atom names in PDB format have element as first character(s)
  // Common patterns: N, O, C, S, CA, CB, CG, CD, CE, NE, OD, OE, etc.
  const firstChar = atomName.charAt(0);
  return firstChar;
}

/**
 * Extract pharmacophore features from ligand contacts
 *
 * Classifies each contact based on interaction type and residue/atom properties:
 * - Aromatic: PHE, TYR, TRP, HIS residues or pi_stacking interactions
 * - Positive: ARG, LYS, HIS with salt_bridge
 * - Negative: ASP, GLU with salt_bridge
 * - Donor: N atoms in hydrogen_bond
 * - Acceptor: O atoms in hydrogen_bond
 * - Hydrophobic: C atoms or hydrophobic interactions
 */
export function extractPharmacophoreFeatures(
  contacts: LigandContact[]
): PharmacophoreFeature[] {
  const features: PharmacophoreFeature[] = [];

  for (const contact of contacts) {
    const resName = extractResidueName(contact.residue);
    const atomElement = getAtomElement(contact.atom);
    let featureType: PharmacophoreType;

    // Priority 1: Pi-stacking interactions are always aromatic
    if (contact.interactionType === 'pi_stacking') {
      featureType = 'aromatic';
    }
    // Priority 2: Salt bridge interactions - check for charged residues
    else if (contact.interactionType === 'salt_bridge') {
      if (POSITIVE_RESIDUES.has(resName)) {
        featureType = 'positive';
      } else if (NEGATIVE_RESIDUES.has(resName)) {
        featureType = 'negative';
      } else {
        // Default to positive/negative based on atom type if residue not recognized
        featureType = atomElement === 'N' ? 'positive' : 'negative';
      }
    }
    // Priority 3: Aromatic residues (even without pi_stacking interaction type)
    else if (AROMATIC_RESIDUES.has(resName)) {
      featureType = 'aromatic';
    }
    // Priority 4: Hydrogen bonds - classify by atom type
    else if (contact.interactionType === 'hydrogen_bond') {
      if (atomElement === 'N') {
        featureType = 'donor';
      } else if (atomElement === 'O') {
        featureType = 'acceptor';
      } else if (atomElement === 'S') {
        // Sulfur can be both donor and acceptor, default to acceptor
        featureType = 'acceptor';
      } else {
        // For other atoms in H-bonds, default to acceptor
        featureType = 'acceptor';
      }
    }
    // Priority 5: Hydrophobic interactions or carbon atoms
    else if (contact.interactionType === 'hydrophobic' || atomElement === 'C') {
      featureType = 'hydrophobic';
    }
    // Default: classify based on atom element
    else {
      if (atomElement === 'N') {
        featureType = 'donor';
      } else if (atomElement === 'O') {
        featureType = 'acceptor';
      } else {
        featureType = 'hydrophobic';
      }
    }

    features.push({
      type: featureType,
      residue: contact.residue,
      chain: contact.chain,
      atom: contact.atom,
      distance: contact.distance,
    });
  }

  return features;
}

/**
 * Classify contacts into primary and secondary shells based on distance
 *
 * Primary shell (<4Å): Direct binding contacts
 * Secondary shell (4-6Å): Extended environment contacts
 * Contacts beyond 6Å are excluded
 */
export function classifyShellResidues(contacts: LigandContact[]): ShellAnalysis {
  const primary: LigandContact[] = [];
  const secondary: LigandContact[] = [];

  for (const contact of contacts) {
    if (contact.distance < PRIMARY_SHELL_CUTOFF) {
      primary.push(contact);
    } else if (contact.distance >= PRIMARY_SHELL_CUTOFF && contact.distance <= SECONDARY_SHELL_CUTOFF) {
      secondary.push(contact);
    }
    // Contacts beyond SECONDARY_SHELL_CUTOFF are excluded
  }

  // Calculate average distances
  const avgPrimaryDistance = primary.length > 0
    ? primary.reduce((sum, c) => sum + c.distance, 0) / primary.length
    : 0;

  const avgSecondaryDistance = secondary.length > 0
    ? secondary.reduce((sum, c) => sum + c.distance, 0) / secondary.length
    : 0;

  return {
    primary,
    secondary,
    stats: {
      primaryCount: primary.length,
      secondaryCount: secondary.length,
      totalCount: primary.length + secondary.length,
      avgPrimaryDistance,
      avgSecondaryDistance,
    },
  };
}

/**
 * Calculate binder quality score based on interaction diversity and density
 *
 * Scores each category with diminishing returns:
 * - H-bonds: max 30 points (3+ H-bonds = max)
 * - Hydrophobic: max 25 points (5+ contacts = max)
 * - Pi-stacking: max 15 points (2+ stacks = max)
 * - Salt bridges: max 15 points (2+ bridges = max)
 * - Burialness: max 15 points (10+ total contacts = max)
 *
 * Rating thresholds:
 * - EXCELLENT: >= 70
 * - GOOD: >= 50
 * - MODERATE: >= 30
 * - POOR: < 30
 */
export function calculateBinderQualityScore(
  contacts: LigandContact[]
): BinderQualityResult {
  // Count interactions by type
  let hbondCount = 0;
  let hydrophobicCount = 0;
  let piStackCount = 0;
  let saltBridgeCount = 0;

  for (const contact of contacts) {
    switch (contact.interactionType) {
      case 'hydrogen_bond':
        hbondCount++;
        break;
      case 'hydrophobic':
        hydrophobicCount++;
        break;
      case 'pi_stacking':
        piStackCount++;
        break;
      case 'salt_bridge':
        saltBridgeCount++;
        break;
      // 'other' type contacts still count toward burialness
    }
  }

  const totalContacts = contacts.length;

  // Calculate scores with diminishing returns (capped)
  const hbondScore = Math.min(30, hbondCount * 10);           // 3+ H-bonds = max
  const hydrophobicScore = Math.min(25, hydrophobicCount * 5); // 5+ contacts = max
  const piStackScore = Math.min(15, piStackCount * 7.5);       // 2+ stacks = max
  const saltBridgeScore = Math.min(15, saltBridgeCount * 7.5); // 2+ bridges = max
  const burialnessScore = Math.min(15, totalContacts * 1.5);   // 10+ contacts = max

  // Total score (max 100)
  const score = hbondScore + hydrophobicScore + piStackScore + saltBridgeScore + burialnessScore;

  // Determine rating
  let rating: QualityRating;
  if (score >= 70) {
    rating = 'EXCELLENT';
  } else if (score >= 50) {
    rating = 'GOOD';
  } else if (score >= 30) {
    rating = 'MODERATE';
  } else {
    rating = 'POOR';
  }

  // Generate suggestions for improvement
  const suggestions: string[] = [];

  if (hbondCount < 3) {
    suggestions.push(`Add ${3 - hbondCount} more H-bond${3 - hbondCount === 1 ? '' : 's'} for improved specificity`);
  }
  if (hydrophobicCount < 5) {
    suggestions.push(`Increase hydrophobic contacts (currently ${hydrophobicCount}, target 5+) for better binding affinity`);
  }
  if (piStackCount < 2) {
    suggestions.push(`Consider adding aromatic residues for pi-stacking interactions`);
  }
  if (saltBridgeCount < 1) {
    suggestions.push(`Add charged residues (Arg/Lys/Asp/Glu) for salt bridge formation`);
  }
  if (totalContacts < 10) {
    suggestions.push(`Binding pocket may be too shallow - aim for 10+ total contacts`);
  }

  return {
    score,
    rating,
    breakdown: {
      hbondScore,
      hydrophobicScore,
      piStackScore,
      saltBridgeScore,
      burialnessScore,
    },
    suggestions,
  };
}
