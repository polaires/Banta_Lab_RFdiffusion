/**
 * Enzyme Analysis Utilities for H-Bond and RASA Conditioning
 *
 * Analyzes PDB files to detect:
 * - Metal ions and their coordination environment
 * - Ligands and their atom classifications
 * - Suggested buried/exposed atoms for RASA conditioning
 * - Suggested H-bond acceptor/donor atoms
 *
 * Used by EnzymeForm for auto-detection with manual override capability
 */

import { findMetalCoordinationFromPDB, METAL_ELEMENTS } from './metalAnalysis';
import type { EnzymeAnalysisResult } from './store';

// Expected coordination numbers for lanthanides vs other metals
const LANTHANIDE_ELEMENTS = new Set([
  'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU'
]);

// Typical coordination numbers by metal type
export const METAL_COORDINATION_INFO: Record<string, { typicalCN: number; range: [number, number]; donors: string[] }> = {
  // Alkali metals
  'NA': { typicalCN: 6, range: [4, 8], donors: ['O', 'N'] },
  'K': { typicalCN: 6, range: [4, 8], donors: ['O', 'N'] },
  // Alkaline earth
  'MG': { typicalCN: 6, range: [4, 6], donors: ['O', 'N'] },
  'CA': { typicalCN: 7, range: [6, 8], donors: ['O', 'N'] },
  // Transition metals
  'ZN': { typicalCN: 4, range: [4, 6], donors: ['N', 'S', 'O'] },
  'FE': { typicalCN: 6, range: [4, 6], donors: ['N', 'S', 'O'] },
  'MN': { typicalCN: 6, range: [4, 6], donors: ['O', 'N'] },
  'CO': { typicalCN: 6, range: [4, 6], donors: ['N', 'S', 'O'] },
  'NI': { typicalCN: 6, range: [4, 6], donors: ['N', 'S', 'O'] },
  'CU': { typicalCN: 4, range: [4, 6], donors: ['N', 'S', 'O'] },
  // Lanthanides - high coordination
  'TB': { typicalCN: 9, range: [8, 9], donors: ['O'] },
  'EU': { typicalCN: 9, range: [8, 9], donors: ['O'] },
  'GD': { typicalCN: 9, range: [8, 9], donors: ['O'] },
  'LA': { typicalCN: 9, range: [9, 10], donors: ['O'] },
  'LU': { typicalCN: 8, range: [8, 9], donors: ['O'] },
  'YB': { typicalCN: 8, range: [8, 9], donors: ['O'] },
  // Default
  'DEFAULT': { typicalCN: 6, range: [4, 8], donors: ['O', 'N'] }
};

// Common ligand atoms that can accept H-bonds (O and N atoms)
const HBOND_ACCEPTOR_ELEMENTS = new Set(['O', 'N', 'S']);

// Common organic ligands that should be analyzed
const KNOWN_LIGANDS = new Set([
  // Carboxylic acids
  'CIT', 'ACE', 'FMT', 'TAR', 'SUC', 'MAL',
  // Nucleotides
  'ATP', 'ADP', 'AMP', 'GTP', 'GDP', 'NAD', 'FAD',
  // Cofactors
  'HEM', 'SAM', 'PLP', 'COA',
  // Sugars
  'GLC', 'GAL', 'MAN', 'FUC',
  // Common substrates
  'EDO', 'GOL', 'PEG', 'MPD'
]);

// Atom names that typically indicate terminal/exposed positions
const TERMINAL_ATOM_PATTERNS = [
  /^O[1-9H]$/,      // O1, O2, O3H, etc. (carboxyl/hydroxyl termini)
  /^N[1-9H]$/,      // N1, N2, etc.
  /^C[1-9]$/,       // C1, C2, etc. (non-coordination carbons)
];

// Atom names that typically indicate coordination/buried positions
const COORDINATION_ATOM_PATTERNS = [
  /^O[0-9]*G$/,     // OG, O1G, O2G (gamma oxygens, often coordinate)
  /^O[ABCD]$/,      // OA, OB, OC, OD (alpha/beta oxygens)
  /^OE[12]$/,       // OE1, OE2 (epsilon carboxyl oxygens)
  /^OD[12]$/,       // OD1, OD2 (delta carboxyl oxygens)
];

/**
 * Parse PDB content and extract atom information
 */
interface PDBAtom {
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

function parsePDBAtoms(pdbContent: string): PDBAtom[] {
  const atoms: PDBAtom[] = [];
  const lines = pdbContent.split('\n');

  for (const line of lines) {
    if (line.startsWith('ATOM') || line.startsWith('HETATM')) {
      try {
        atoms.push({
          serial: parseInt(line.slice(6, 11).trim()),
          name: line.slice(12, 16).trim(),
          resName: line.slice(17, 20).trim(),
          chainId: line.slice(21, 22).trim() || 'A',
          resSeq: parseInt(line.slice(22, 26).trim()),
          x: parseFloat(line.slice(30, 38).trim()),
          y: parseFloat(line.slice(38, 46).trim()),
          z: parseFloat(line.slice(46, 54).trim()),
          element: line.slice(76, 78).trim().toUpperCase() || line.slice(12, 14).trim().replace(/[0-9]/g, '').toUpperCase(),
          isHETATM: line.startsWith('HETATM')
        });
      } catch {
        // Skip malformed lines
      }
    }
  }

  return atoms;
}

/**
 * Calculate distance between two atoms
 */
function atomDistance(a1: PDBAtom, a2: PDBAtom): number {
  const dx = a1.x - a2.x;
  const dy = a1.y - a2.y;
  const dz = a1.z - a2.z;
  return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

/**
 * Detect ligands from PDB content
 * Returns non-protein, non-water, non-metal HETATM residues
 */
function detectLigands(atoms: PDBAtom[]): Map<string, PDBAtom[]> {
  const ligandMap = new Map<string, PDBAtom[]>();

  // Group HETATM atoms by residue
  for (const atom of atoms) {
    if (!atom.isHETATM) continue;

    // Skip water and metals
    if (atom.resName === 'HOH' || atom.resName === 'WAT') continue;
    if (METAL_ELEMENTS.has(atom.element)) continue;

    // Create unique key for residue
    const key = `${atom.chainId}_${atom.resName}_${atom.resSeq}`;

    if (!ligandMap.has(key)) {
      ligandMap.set(key, []);
    }
    ligandMap.get(key)!.push(atom);
  }

  return ligandMap;
}

/**
 * Classify ligand atoms as coordination vs terminal
 * based on proximity to metals and atom naming patterns
 */
function classifyLigandAtoms(
  ligandAtoms: PDBAtom[],
  metalPositions: Array<{ x: number; y: number; z: number; element: string }>
): {
  buried: string[];      // Near metal, should be buried
  exposed: string[];     // Terminal, can be exposed
  hbondAcceptors: string[];  // O/N atoms that can accept H-bonds
} {
  const buried: string[] = [];
  const exposed: string[] = [];
  const hbondAcceptors: string[] = [];

  const METAL_DISTANCE_THRESHOLD = 3.5; // Å - atoms within this range of metal are coordinating

  for (const atom of ligandAtoms) {
    // Check if atom can be H-bond acceptor (O, N, S)
    if (HBOND_ACCEPTOR_ELEMENTS.has(atom.element)) {
      hbondAcceptors.push(atom.name);
    }

    // Check distance to all metals
    let nearMetal = false;
    for (const metal of metalPositions) {
      const dist = Math.sqrt(
        (atom.x - metal.x) ** 2 +
        (atom.y - metal.y) ** 2 +
        (atom.z - metal.z) ** 2
      );
      if (dist <= METAL_DISTANCE_THRESHOLD && HBOND_ACCEPTOR_ELEMENTS.has(atom.element)) {
        nearMetal = true;
        break;
      }
    }

    if (nearMetal) {
      // Atom is coordinating metal - should be buried
      buried.push(atom.name);
    } else if (HBOND_ACCEPTOR_ELEMENTS.has(atom.element)) {
      // Non-coordinating O/N/S - potentially exposed for substrate access
      exposed.push(atom.name);
    }

    // Also check atom naming patterns for edge cases
    if (!nearMetal && buried.indexOf(atom.name) === -1) {
      for (const pattern of COORDINATION_ATOM_PATTERNS) {
        if (pattern.test(atom.name) && HBOND_ACCEPTOR_ELEMENTS.has(atom.element)) {
          // Pattern suggests coordination but not near metal - might be in a different conformation
          // Don't auto-add, let user decide
          break;
        }
      }
    }
  }

  return {
    buried: [...new Set(buried)],      // Remove duplicates
    exposed: [...new Set(exposed)],
    hbondAcceptors: [...new Set(hbondAcceptors)]
  };
}

/**
 * Get coordination number mismatch warning for metal replacement
 */
export function getCoordinationWarning(
  currentMetal: string,
  targetMetal: string,
  currentCN: number
): { warning: string; severity: 'info' | 'warning' | 'error' } | null {
  const currentInfo = METAL_COORDINATION_INFO[currentMetal] || METAL_COORDINATION_INFO['DEFAULT'];
  const targetInfo = METAL_COORDINATION_INFO[targetMetal] || METAL_COORDINATION_INFO['DEFAULT'];

  const cnDiff = targetInfo.typicalCN - currentCN;

  if (cnDiff === 0) {
    return null; // Perfect match
  }

  if (cnDiff > 0) {
    // Target metal needs MORE coordination
    if (LANTHANIDE_ELEMENTS.has(targetMetal) && !LANTHANIDE_ELEMENTS.has(currentMetal)) {
      return {
        warning: `${targetMetal} requires CN ${targetInfo.typicalCN} (current: ${currentCN}). ` +
                `RFD3 will need to design ${cnDiff} additional Asp/Glu residues.`,
        severity: cnDiff > 3 ? 'error' : 'warning'
      };
    }
    return {
      warning: `${targetMetal} typically has CN ${targetInfo.typicalCN} (current: ${currentCN}). ` +
              `Consider "Explore" mode to find additional coordinators.`,
      severity: 'info'
    };
  }

  // Target metal needs LESS coordination
  return {
    warning: `${targetMetal} typically has CN ${targetInfo.typicalCN} (current: ${currentCN}). ` +
            `Some existing coordinators may not be used.`,
    severity: 'info'
  };
}

/**
 * Check if a metal is a lanthanide (rare earth)
 */
export function isLanthanide(element: string): boolean {
  return LANTHANIDE_ELEMENTS.has(element.toUpperCase());
}

/**
 * Get recommended donor residue types for a metal
 */
export function getRecommendedDonors(element: string): string[] {
  const info = METAL_COORDINATION_INFO[element.toUpperCase()] || METAL_COORDINATION_INFO['DEFAULT'];

  const donorResidues: string[] = [];

  if (info.donors.includes('O')) {
    donorResidues.push('Asp', 'Glu', 'Asn', 'Gln', 'Ser', 'Thr', 'Tyr');
  }
  if (info.donors.includes('N')) {
    donorResidues.push('His');
  }
  if (info.donors.includes('S')) {
    donorResidues.push('Cys', 'Met');
  }

  return donorResidues;
}

/**
 * Main analysis function: Analyze enzyme structure for conditioning
 *
 * @param pdbContent - Raw PDB file content
 * @returns EnzymeAnalysisResult with metals and ligands with suggestions
 */
export function analyzeEnzymeStructure(pdbContent: string): EnzymeAnalysisResult {
  // Parse atoms
  const atoms = parsePDBAtoms(pdbContent);

  // Find metals using existing metalAnalysis
  const metalCoordinations = findMetalCoordinationFromPDB(pdbContent, 3.5);

  // Detect ligands
  const ligandMap = detectLigands(atoms);

  // Convert metal coordinations to simpler format for EnzymeAnalysisResult
  const metals = metalCoordinations.map(mc => ({
    element: mc.element,
    chain: mc.chainId,
    residueNumber: mc.resSeq,
    coordinatingAtoms: mc.coordinating
      .filter(c => !c.isWater)  // Exclude water from coordination list
      .map(c => ({
        chain: c.chain,
        residue: parseInt(c.residue.replace(/[A-Z]/g, '')),
        atom: c.atom,
        distance: c.distance
      })),
    coordinationNumber: mc.coordinating.filter(c => !c.isWater).length
  }));

  // Get metal positions for ligand classification
  const metalPositions = metalCoordinations.map(mc => ({
    x: mc.pos[0],
    y: mc.pos[1],
    z: mc.pos[2],
    element: mc.element
  }));

  // Analyze each ligand
  const ligands: EnzymeAnalysisResult['ligands'] = [];

  for (const [key, ligandAtoms] of ligandMap) {
    const [chainId, resName, resSeqStr] = key.split('_');
    const resSeq = parseInt(resSeqStr);

    // Get unique atoms with their elements
    const atomsInfo = ligandAtoms.map(a => ({
      name: a.name,
      element: a.element
    }));

    // Remove duplicate atom names (keep first occurrence)
    const uniqueAtoms = atomsInfo.filter(
      (a, i, arr) => arr.findIndex(x => x.name === a.name) === i
    );

    // Classify atoms
    const classification = classifyLigandAtoms(ligandAtoms, metalPositions);

    ligands.push({
      name: resName,
      chain: chainId,
      residueNumber: resSeq,
      atoms: uniqueAtoms,
      suggestedBuried: classification.buried,
      suggestedExposed: classification.exposed,
      suggestedHBondAcceptors: classification.hbondAcceptors
    });
  }

  return {
    metals,
    ligands
  };
}

/**
 * Replace metal in PDB content
 * Used for metal replacement workflow (e.g., MG → TB)
 *
 * PDB HETATM format (columns are 1-indexed):
 * - Columns 1-6: Record type (HETATM)
 * - Columns 7-11: Atom serial number
 * - Columns 13-16: Atom name (RIGHT-justified for 1-2 char names)
 * - Column 17: Alternate location indicator
 * - Columns 18-20: Residue name (RIGHT-justified)
 * - Column 22: Chain identifier
 * - Columns 23-26: Residue sequence number
 * - Columns 77-78: Element symbol (RIGHT-justified)
 *
 * For metal ions, atom name should match the element symbol.
 */
export function replaceMetalInPdb(
  pdbContent: string,
  sourceMetal: string,
  targetMetal: string
): string {
  const sourceUpper = sourceMetal.toUpperCase();
  const targetUpper = targetMetal.toUpperCase();

  // Replace metal element in HETATM lines
  const lines = pdbContent.split('\n');
  const modifiedLines = lines.map(line => {
    if (line.startsWith('HETATM')) {
      const atomName = line.slice(12, 16).trim();
      const resName = line.slice(17, 20).trim();
      const element = line.slice(76, 78).trim().toUpperCase();

      // Check if this is the metal to replace (by residue name, element, or atom name)
      if (resName === sourceUpper || element === sourceUpper || atomName === sourceUpper) {
        // Ensure line is long enough
        let modifiedLine = line.padEnd(80);

        // Replace atom name (columns 13-16, 0-indexed: 12-16)
        // Metal atom names are typically right-justified in a 4-char field
        const paddedAtomName = targetUpper.padStart(2).padEnd(4);
        modifiedLine = modifiedLine.slice(0, 12) + paddedAtomName + modifiedLine.slice(16);

        // Replace residue name (columns 18-20, 0-indexed: 17-20)
        const paddedResName = targetUpper.padStart(3);
        modifiedLine = modifiedLine.slice(0, 17) + paddedResName + modifiedLine.slice(20);

        // Replace element symbol (columns 77-78, 0-indexed: 76-78)
        const paddedElement = targetUpper.padStart(2);
        modifiedLine = modifiedLine.slice(0, 76) + paddedElement + modifiedLine.slice(78);

        return modifiedLine.trimEnd();
      }
    }
    return line;
  });

  return modifiedLines.join('\n');
}

/**
 * Build select_buried parameter for RFD3 request
 */
export function buildBuriedSelection(
  selections: Record<string, string>
): Record<string, string> | undefined {
  const filtered: Record<string, string> = {};

  for (const [key, atoms] of Object.entries(selections)) {
    if (atoms && atoms.trim()) {
      filtered[key] = atoms;
    }
  }

  return Object.keys(filtered).length > 0 ? filtered : undefined;
}

/**
 * Build select_exposed parameter for RFD3 request
 */
export function buildExposedSelection(
  selections: Record<string, string>
): Record<string, string> | undefined {
  return buildBuriedSelection(selections); // Same format
}

/**
 * Build select_hbond_acceptor parameter for RFD3 request
 */
export function buildHBondAcceptorSelection(
  selections: Record<string, string>
): Record<string, string> | undefined {
  return buildBuriedSelection(selections); // Same format
}

/**
 * Build select_hbond_donor parameter for RFD3 request
 */
export function buildHBondDonorSelection(
  selections: Record<string, string>
): Record<string, string> | undefined {
  return buildBuriedSelection(selections); // Same format
}

/**
 * Validate that buried and exposed selections don't overlap
 */
export function validateBuriedExposedSelections(
  buried: Record<string, string>,
  exposed: Record<string, string>
): { valid: boolean; errors: string[] } {
  const errors: string[] = [];

  for (const key of Object.keys(buried)) {
    if (exposed[key]) {
      const buriedAtoms = new Set(buried[key].split(',').map(s => s.trim()));
      const exposedAtoms = exposed[key].split(',').map(s => s.trim());

      for (const atom of exposedAtoms) {
        if (buriedAtoms.has(atom)) {
          errors.push(`${key}:${atom} cannot be both buried and exposed`);
        }
      }
    }
  }

  return {
    valid: errors.length === 0,
    errors
  };
}

/**
 * Suggest protein donor atoms for H-bond conditioning
 * Based on amino acid type and common H-bond donor patterns
 */
export function suggestProteinDonors(
  aminoAcid: string
): Array<{ atom: string; description: string }> {
  const aa = aminoAcid.toUpperCase();

  const donors: Array<{ atom: string; description: string }> = [
    { atom: 'N', description: 'Backbone amide NH' }
  ];

  switch (aa) {
    case 'SER':
      donors.push({ atom: 'OG', description: 'Serine hydroxyl' });
      break;
    case 'THR':
      donors.push({ atom: 'OG1', description: 'Threonine hydroxyl' });
      break;
    case 'TYR':
      donors.push({ atom: 'OH', description: 'Tyrosine hydroxyl' });
      break;
    case 'ASN':
      donors.push({ atom: 'ND2', description: 'Asparagine amide NH2' });
      break;
    case 'GLN':
      donors.push({ atom: 'NE2', description: 'Glutamine amide NH2' });
      break;
    case 'HIS':
      donors.push({ atom: 'ND1', description: 'Histidine ND1' });
      donors.push({ atom: 'NE2', description: 'Histidine NE2' });
      break;
    case 'TRP':
      donors.push({ atom: 'NE1', description: 'Tryptophan indole NH' });
      break;
    case 'ARG':
      donors.push({ atom: 'NH1', description: 'Arginine guanidinium NH' });
      donors.push({ atom: 'NH2', description: 'Arginine guanidinium NH' });
      donors.push({ atom: 'NE', description: 'Arginine epsilon NH' });
      break;
    case 'LYS':
      donors.push({ atom: 'NZ', description: 'Lysine amine NH3+' });
      break;
    case 'CYS':
      donors.push({ atom: 'SG', description: 'Cysteine thiol (weak donor)' });
      break;
  }

  return donors;
}

// ============================================================================
// METAL REPLACEMENT PRESETS
// ============================================================================

/**
 * Preset types for metal replacement workflow
 * These optimize all parameters based on the coordination number difference
 */
export type MetalReplacementPreset =
  | 'keep_coordination'      // Same CN metals (Mg→Ca, Zn→Co)
  | 'expand_coordination'    // Small→Large CN (Mg→Tb, Ca→Eu)
  | 'reduce_coordination'    // Large→Small CN (rare)
  | 'lanthanide_design'      // Optimized for lanthanide binding
  | 'custom';                // Manual control

/**
 * Preset configuration that defines what parameters to set
 */
export interface PresetConfig {
  name: string;
  description: string;
  /** Whether to fix original metal-coordinating residues */
  fixMetalCoordinators: boolean;
  /** Whether to fix ligand position */
  fixLigandPosition: boolean;
  /** Coordination mode for RFD3 */
  coordinationMode: 'keep' | 'explore' | 'hybrid';
  /** Whether to enable RASA conditioning */
  enableRASA: boolean;
  /** Whether to bury the metal */
  buryMetal: boolean;
  /** Whether to enable H-bond conditioning */
  enableHBonds: boolean;
  /** Recommended design length multiplier (relative to original) */
  lengthMultiplier: [number, number];
  /** Number of extra coordinating residues needed */
  extraCoordinatorsNeeded: number;
  /** Warning message to display */
  warning?: string;
  /** Tips for the user */
  tips: string[];
}

/**
 * Preset definitions with full configuration
 */
export const METAL_REPLACEMENT_PRESETS: Record<MetalReplacementPreset, PresetConfig> = {
  keep_coordination: {
    name: 'Keep Coordination',
    description: 'For metals with similar coordination (Mg→Ca, Zn→Co, Fe→Mn)',
    fixMetalCoordinators: true,
    fixLigandPosition: true,
    coordinationMode: 'keep',
    enableRASA: true,
    buryMetal: true,
    enableHBonds: true,
    lengthMultiplier: [0.9, 1.1],
    extraCoordinatorsNeeded: 0,
    tips: [
      'Original coordinating residues will be preserved',
      'Good for similar-sized metal substitutions',
      'Minimal scaffold redesign needed'
    ]
  },
  expand_coordination: {
    name: 'Expand Coordination',
    description: 'For replacing small metals with larger CN requirements',
    fixMetalCoordinators: false,  // Don't fix - need new coordinators
    fixLigandPosition: true,
    coordinationMode: 'explore',
    enableRASA: true,
    buryMetal: true,
    enableHBonds: true,
    lengthMultiplier: [1.0, 1.3],
    extraCoordinatorsNeeded: 3,  // Will be calculated dynamically
    warning: 'Original coordinating residues will NOT be fixed to allow RFD3 to design new coordinators',
    tips: [
      'RFD3 will design additional Asp/Glu residues for coordination',
      'Generate multiple designs and filter by coordination quality',
      'Consider running PyRosetta analysis on outputs'
    ]
  },
  reduce_coordination: {
    name: 'Reduce Coordination',
    description: 'For replacing large CN metals with smaller requirements',
    fixMetalCoordinators: false,
    fixLigandPosition: true,
    coordinationMode: 'explore',
    enableRASA: true,
    buryMetal: true,
    enableHBonds: false,  // Less H-bonds needed
    lengthMultiplier: [0.8, 1.0],
    extraCoordinatorsNeeded: 0,
    warning: 'Some original coordinators may become redundant',
    tips: [
      'Fewer coordinating residues needed',
      'May result in more open active site',
      'Consider which ligand atoms to keep coordinated'
    ]
  },
  lanthanide_design: {
    name: 'Lanthanide Design',
    description: 'Optimized for Tb, Eu, Gd, La and other lanthanides (CN=8-9)',
    fixMetalCoordinators: false,  // Critical: need new coordinators
    fixLigandPosition: true,
    coordinationMode: 'explore',
    enableRASA: true,
    buryMetal: true,
    enableHBonds: true,
    lengthMultiplier: [1.1, 1.4],
    extraCoordinatorsNeeded: 3,
    warning: 'Lanthanides require CN=8-9. RFD3 will design 2-3 additional Asp/Glu coordinators.',
    tips: [
      'Lanthanides prefer oxygen donors (Asp/Glu carboxylates)',
      'Original MG/CA coordinators are too few - will be redesigned',
      'Use longer scaffolds to accommodate larger coordination sphere',
      'H-bond conditioning helps position coordinating residues'
    ]
  },
  custom: {
    name: 'Custom',
    description: 'Manual control over all parameters',
    fixMetalCoordinators: true,
    fixLigandPosition: true,
    coordinationMode: 'explore',
    enableRASA: false,
    buryMetal: false,
    enableHBonds: false,
    lengthMultiplier: [1.0, 1.0],
    extraCoordinatorsNeeded: 0,
    tips: [
      'Configure each parameter manually',
      'Use for specialized workflows',
      'Refer to RFD3 documentation for parameter details'
    ]
  }
};

/**
 * Calculate coordination number difference between two metals
 */
export function calculateCNDifference(sourceMetal: string, targetMetal: string): number {
  const sourceInfo = METAL_COORDINATION_INFO[sourceMetal.toUpperCase()] || METAL_COORDINATION_INFO['DEFAULT'];
  const targetInfo = METAL_COORDINATION_INFO[targetMetal.toUpperCase()] || METAL_COORDINATION_INFO['DEFAULT'];
  return targetInfo.typicalCN - sourceInfo.typicalCN;
}

/**
 * Get coordination info for a metal
 */
export function getMetalCoordinationInfo(metal: string): { typicalCN: number; range: [number, number]; donors: string[] } {
  return METAL_COORDINATION_INFO[metal.toUpperCase()] || METAL_COORDINATION_INFO['DEFAULT'];
}

// Note: isLanthanide is already defined above (line ~269), don't duplicate

/**
 * Recommend the best preset based on source and target metals
 */
export function recommendPreset(sourceMetal: string, targetMetal: string): {
  preset: MetalReplacementPreset;
  reason: string;
  cnDifference: number;
  extraCoordinatorsNeeded: number;
} {
  const source = sourceMetal.toUpperCase();
  const target = targetMetal.toUpperCase();
  const cnDiff = calculateCNDifference(source, target);
  const targetIsLanthanide = isLanthanide(target);

  // Calculate extra coordinators needed
  const sourceInfo = getMetalCoordinationInfo(source);
  const targetInfo = getMetalCoordinationInfo(target);
  const extraNeeded = Math.max(0, targetInfo.typicalCN - sourceInfo.typicalCN);

  // Special case: target is lanthanide
  if (targetIsLanthanide) {
    return {
      preset: 'lanthanide_design',
      reason: `${target} is a lanthanide requiring CN=${targetInfo.typicalCN}. ${source} has CN=${sourceInfo.typicalCN}. Need ${extraNeeded} additional coordinators.`,
      cnDifference: cnDiff,
      extraCoordinatorsNeeded: extraNeeded
    };
  }

  // CN increases significantly (need more coordinators)
  if (cnDiff >= 2) {
    return {
      preset: 'expand_coordination',
      reason: `Target ${target} needs ${cnDiff} more coordinators than ${source} (CN: ${sourceInfo.typicalCN}→${targetInfo.typicalCN})`,
      cnDifference: cnDiff,
      extraCoordinatorsNeeded: extraNeeded
    };
  }

  // CN decreases significantly
  if (cnDiff <= -2) {
    return {
      preset: 'reduce_coordination',
      reason: `Target ${target} needs ${Math.abs(cnDiff)} fewer coordinators than ${source} (CN: ${sourceInfo.typicalCN}→${targetInfo.typicalCN})`,
      cnDifference: cnDiff,
      extraCoordinatorsNeeded: 0
    };
  }

  // Similar CN - keep coordination
  return {
    preset: 'keep_coordination',
    reason: `${source} and ${target} have similar coordination requirements (CN: ${sourceInfo.typicalCN}→${targetInfo.typicalCN})`,
    cnDifference: cnDiff,
    extraCoordinatorsNeeded: 0
  };
}

/**
 * Get preset configuration with dynamic values based on metals
 */
export function getPresetConfig(
  preset: MetalReplacementPreset,
  sourceMetal?: string,
  targetMetal?: string
): PresetConfig {
  const baseConfig = { ...METAL_REPLACEMENT_PRESETS[preset] };

  // Calculate dynamic values if metals provided
  if (sourceMetal && targetMetal && preset !== 'custom') {
    const sourceInfo = getMetalCoordinationInfo(sourceMetal);
    const targetInfo = getMetalCoordinationInfo(targetMetal);
    baseConfig.extraCoordinatorsNeeded = Math.max(0, targetInfo.typicalCN - sourceInfo.typicalCN);

    // Update warning with specific values
    if (preset === 'lanthanide_design' || preset === 'expand_coordination') {
      baseConfig.warning = `${targetMetal.toUpperCase()} requires CN=${targetInfo.typicalCN} (current: ${sourceInfo.typicalCN}). RFD3 will design ${baseConfig.extraCoordinatorsNeeded} additional coordinating residues.`;
    }
  }

  return baseConfig;
}

/**
 * Common target metals for replacement, organized by category
 */
export const TARGET_METAL_OPTIONS = {
  'Alkaline Earth': [
    { symbol: 'MG', name: 'Magnesium', cn: 6 },
    { symbol: 'CA', name: 'Calcium', cn: 7 },
  ],
  'Transition Metals': [
    { symbol: 'ZN', name: 'Zinc', cn: 4 },
    { symbol: 'FE', name: 'Iron', cn: 6 },
    { symbol: 'MN', name: 'Manganese', cn: 6 },
    { symbol: 'CO', name: 'Cobalt', cn: 6 },
    { symbol: 'NI', name: 'Nickel', cn: 6 },
    { symbol: 'CU', name: 'Copper', cn: 4 },
  ],
  'Lanthanides': [
    { symbol: 'TB', name: 'Terbium', cn: 9 },
    { symbol: 'EU', name: 'Europium', cn: 9 },
    { symbol: 'GD', name: 'Gadolinium', cn: 9 },
    { symbol: 'LA', name: 'Lanthanum', cn: 9 },
    { symbol: 'YB', name: 'Ytterbium', cn: 8 },
    { symbol: 'LU', name: 'Lutetium', cn: 8 },
  ]
};

/**
 * Filter catalytic residues based on preset
 * For expand/lanthanide presets, remove metal-coordinating residues to let RFD3 redesign them
 */
export function filterCatalyticResiduesForPreset(
  allResidues: Array<{ chain: string; residue: number; name: string; role?: string }>,
  preset: MetalReplacementPreset,
  metalCoordinators: Array<{ chain: string; residue: number }>
): Array<{ chain: string; residue: number; name: string; role?: string; excluded?: boolean; excludeReason?: string }> {
  const config = METAL_REPLACEMENT_PRESETS[preset];

  if (config.fixMetalCoordinators) {
    // Keep all residues
    return allResidues.map(r => ({ ...r, excluded: false }));
  }

  // Mark metal coordinators as excluded for expand/lanthanide presets
  const coordinatorSet = new Set(
    metalCoordinators.map(c => `${c.chain}${c.residue}`)
  );

  return allResidues.map(r => {
    const key = `${r.chain}${r.residue}`;
    const isCoordinator = coordinatorSet.has(key);
    return {
      ...r,
      excluded: isCoordinator,
      excludeReason: isCoordinator ? 'Metal coordinator - will be redesigned for new CN' : undefined
    };
  });
}
