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

export interface PharmacophoreFeature {
  type: PharmacophoreType;
  residue: string;
  chain: string;
  atom: string;
  distance: number;
  position?: [number, number, number];
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
