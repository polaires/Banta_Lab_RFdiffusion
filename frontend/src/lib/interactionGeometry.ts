/**
 * Interaction Line Geometry Computation
 * Computes 3D line endpoints for visualizing protein-ligand interactions
 */

import type { InteractionType } from './ligandAnalysis';

export interface InteractionLineContact {
  ligandChain: string;
  ligandResSeq: number;
  ligandAtom: string;
  proteinResidue: string;
  proteinChain: string;
  proteinAtom: string;
  interactionType: InteractionType;
  distance: number;
}

export interface InteractionLine {
  start: [number, number, number];
  end: [number, number, number];
  type: InteractionType;
  color: number;
  dashed: boolean;
  label?: string;
}

// Interaction type colors (matching ligandAnalysis)
const INTERACTION_COLORS: Record<InteractionType, number> = {
  hydrogen_bond: 0x3B82F6,  // Blue
  hydrophobic: 0x22C55E,    // Green
  salt_bridge: 0xEF4444,    // Red
  pi_stacking: 0xA855F7,    // Purple
  other: 0x6B7280,          // Gray
};

// Dashed lines for weaker interactions
const DASHED_TYPES: Set<InteractionType> = new Set(['hydrophobic', 'other']);

interface ParsedAtom {
  name: string;
  resName: string;
  chainId: string;
  resSeq: number;
  x: number;
  y: number;
  z: number;
}

/**
 * Parse PDB to extract atom positions
 */
function parseAtomPositions(pdbContent: string): Map<string, ParsedAtom> {
  const atoms = new Map<string, ParsedAtom>();
  const lines = pdbContent.split('\n');

  for (const line of lines) {
    if (!line.startsWith('ATOM') && !line.startsWith('HETATM')) continue;

    try {
      const atom: ParsedAtom = {
        name: line.slice(12, 16).trim(),
        resName: line.slice(17, 20).trim(),
        chainId: line.slice(21, 22).trim(),
        resSeq: parseInt(line.slice(22, 26).trim()),
        x: parseFloat(line.slice(30, 38).trim()),
        y: parseFloat(line.slice(38, 46).trim()),
        z: parseFloat(line.slice(46, 54).trim()),
      };

      // Create unique key for lookup
      const key = `${atom.chainId}:${atom.resSeq}:${atom.name}`;
      atoms.set(key, atom);
    } catch {
      // Skip malformed lines
    }
  }

  return atoms;
}

/**
 * Extract residue number from residue string (e.g., "SER10" -> 10)
 */
function extractResSeq(residue: string): number {
  const match = residue.match(/\d+/);
  return match ? parseInt(match[0]) : 0;
}

/**
 * Compute 3D line segments for interaction visualization
 */
export function computeInteractionLines(
  pdbContent: string,
  contacts: InteractionLineContact[]
): InteractionLine[] {
  if (!pdbContent || contacts.length === 0) {
    return [];
  }

  const atoms = parseAtomPositions(pdbContent);
  const lines: InteractionLine[] = [];

  for (const contact of contacts) {
    // Find ligand atom position
    const ligandKey = `${contact.ligandChain}:${contact.ligandResSeq}:${contact.ligandAtom}`;
    const ligandAtom = atoms.get(ligandKey);

    // Find protein atom position
    const proteinResSeq = extractResSeq(contact.proteinResidue);
    const proteinKey = `${contact.proteinChain}:${proteinResSeq}:${contact.proteinAtom}`;
    const proteinAtom = atoms.get(proteinKey);

    if (!ligandAtom || !proteinAtom) {
      continue; // Skip if atoms not found
    }

    lines.push({
      start: [ligandAtom.x, ligandAtom.y, ligandAtom.z],
      end: [proteinAtom.x, proteinAtom.y, proteinAtom.z],
      type: contact.interactionType,
      color: INTERACTION_COLORS[contact.interactionType],
      dashed: DASHED_TYPES.has(contact.interactionType),
      label: `${contact.proteinResidue} ${contact.proteinAtom}`,
    });
  }

  return lines;
}
