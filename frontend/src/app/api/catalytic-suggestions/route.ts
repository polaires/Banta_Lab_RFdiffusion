/**
 * Catalytic residue suggestions API endpoint.
 * Uses local binding pocket detection to find residues near ligands.
 * M-CSA is used for known enzymes as a secondary source.
 */

import { NextRequest, NextResponse } from 'next/server';

interface SuggestionRequest {
  pdb_id: string | null;
  pdb_content: string;
  /** Distance cutoff for organic ligands only (default: 5.0 Å). Metals always use tight 3.5 Å. */
  cutoff?: number;
}

interface CatalyticResidue {
  chain: string;
  residue: number;
  name: string;
  role: string | null;
  confidence: number;
  source: 'mcsa' | 'local';
  ligandCode?: string;  // The ligand/metal this residue is associated with
}

interface SuggestionResponse {
  source: 'mcsa' | 'local' | 'none';
  residues: CatalyticResidue[];
}

const MCSA_API_BASE = 'https://www.ebi.ac.uk/thornton-srv/m-csa/api';

// Default distance cutoff for binding pocket detection around organic ligands (Angstroms)
const DEFAULT_LIGAND_CUTOFF = 5.0;

// Tight cutoff for direct metal coordination (Angstroms)
// Metal-ligand bonds are typically 2.0-2.7 Å; 3.5 Å captures direct coordinators
const METAL_COORDINATION_CUTOFF = 3.5;

// Common solvent/ion residues to exclude from ligand detection
const EXCLUDED_LIGANDS = new Set([
  'HOH', 'WAT', 'H2O', 'DOD', 'DIS',  // Water
  'NA', 'CL', 'K', 'CA', 'ZN', 'MG', 'FE', 'CU', 'MN', 'CO', 'NI',  // Common ions (we want to find residues AROUND these)
  'SO4', 'PO4', 'GOL', 'EDO', 'PEG', 'DMS', 'ACT', 'ACE', 'BME',  // Common crystallization additives
]);

// Metal ions that we want to find binding residues around
const METAL_IONS = new Set(['MG', 'ZN', 'FE', 'CU', 'MN', 'CO', 'NI', 'CA']);

interface Atom {
  x: number;
  y: number;
  z: number;
  resName: string;
  chainId: string;
  resSeq: number;
  isHetatm: boolean;
}

/**
 * Parse PDB content to extract atom coordinates.
 */
function parsePdbAtoms(pdbContent: string): { proteinAtoms: Atom[]; ligandAtoms: Atom[]; metalAtoms: Atom[] } {
  const proteinAtoms: Atom[] = [];
  const ligandAtoms: Atom[] = [];
  const metalAtoms: Atom[] = [];
  const lines = pdbContent.split('\n');

  for (const line of lines) {
    const recordType = line.slice(0, 6).trim();
    if (recordType !== 'ATOM' && recordType !== 'HETATM') continue;

    const x = parseFloat(line.slice(30, 38));
    const y = parseFloat(line.slice(38, 46));
    const z = parseFloat(line.slice(46, 54));
    const resName = line.slice(17, 20).trim();
    const chainId = line.slice(21, 22).trim() || 'A';
    const resSeq = parseInt(line.slice(22, 26).trim(), 10);

    if (isNaN(x) || isNaN(y) || isNaN(z) || isNaN(resSeq)) continue;

    const atom: Atom = { x, y, z, resName, chainId, resSeq, isHetatm: recordType === 'HETATM' };

    if (recordType === 'ATOM') {
      // Standard amino acid
      proteinAtoms.push(atom);
    } else {
      // HETATM - could be ligand, metal, or solvent
      if (METAL_IONS.has(resName)) {
        metalAtoms.push(atom);
      } else if (!EXCLUDED_LIGANDS.has(resName)) {
        ligandAtoms.push(atom);
      }
    }
  }

  return { proteinAtoms, ligandAtoms, metalAtoms };
}

/**
 * Calculate distance between two atoms.
 */
function distance(a1: Atom, a2: Atom): number {
  const dx = a1.x - a2.x;
  const dy = a1.y - a2.y;
  const dz = a1.z - a2.z;
  return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

/**
 * Find binding pocket residues near ligands and metals.
 * Uses SMART cutoffs: tight 3.5 Å for direct metal coordination, broader cutoff for organic ligands.
 * Only analyzes the first protein chain to avoid duplicates from multiple biological assemblies.
 * Tracks which ligand/metal each residue is associated with for filtering.
 * @param pdbContent - PDB file content
 * @param ligandCutoff - Distance cutoff for organic ligands (default: 5.0 Å)
 */
function findBindingPocketResidues(pdbContent: string, ligandCutoff: number = DEFAULT_LIGAND_CUTOFF): CatalyticResidue[] {
  const { proteinAtoms, ligandAtoms, metalAtoms } = parsePdbAtoms(pdbContent);

  if (proteinAtoms.length === 0) {
    console.log('[Local] No protein atoms found');
    return [];
  }

  // Find unique protein chains and use only the first one (alphabetically)
  const proteinChains = [...new Set(proteinAtoms.map(a => a.chainId))].sort();
  const targetChain = proteinChains[0];

  console.log(`[Local] Found protein chains: ${proteinChains.join(', ')}. Using chain ${targetChain}`);

  // Filter protein atoms to target chain only
  const chainProteinAtoms = proteinAtoms.filter(a => a.chainId === targetChain);

  // Combine ligands and metals as targets
  const allTargets = [...ligandAtoms, ...metalAtoms];

  if (allTargets.length === 0) {
    console.log('[Local] No ligands/metals found');
    return [];
  }

  // Find ligands/metals that are closest to our target chain (within 10Å of any protein atom)
  // This associates ligands with the correct protein chain
  const LIGAND_ASSOCIATION_CUTOFF = 10.0;
  const chainTargets: Atom[] = [];

  for (const target of allTargets) {
    let minDistToChain = Infinity;
    for (const proteinAtom of chainProteinAtoms) {
      const dist = distance(target, proteinAtom);
      if (dist < minDistToChain) {
        minDistToChain = dist;
      }
    }
    // Only include ligands that are close to our target chain
    if (minDistToChain <= LIGAND_ASSOCIATION_CUTOFF) {
      chainTargets.push(target);
    }
  }

  if (chainTargets.length === 0) {
    console.log(`[Local] No ligands/metals associated with chain ${targetChain}`);
    return [];
  }

  // Separate chain-associated metals from organic ligands
  const chainMetals = chainTargets.filter(a => METAL_IONS.has(a.resName));
  const chainLigands = chainTargets.filter(a => !METAL_IONS.has(a.resName));

  // Log what we found
  const metalNames = [...new Set(chainMetals.map(a => a.resName))];
  const ligandNames = [...new Set(chainLigands.map(a => a.resName))];
  console.log(`[Local] Chain ${targetChain} - Metals: ${metalNames.join(', ') || 'none'}, Ligands: ${ligandNames.join(', ') || 'none'}`);
  console.log(`[Local] Using cutoffs - Metals: ${METAL_COORDINATION_CUTOFF}Å (direct coordination), Ligands: ${ligandCutoff}Å`);

  // Find unique residues with SMART cutoffs:
  // - Metals: tight 3.5 Å for direct coordinators only
  // - Ligands: broader cutoff for binding pocket
  const residueMap = new Map<string, {
    chain: string;
    residue: number;
    name: string;
    minDist: number;
    ligandCode: string;
    isMetal: boolean;  // Track if this is a metal coordinator
  }>();

  // First pass: find direct metal coordinators (tight cutoff)
  for (const proteinAtom of chainProteinAtoms) {
    for (const metalAtom of chainMetals) {
      const dist = distance(proteinAtom, metalAtom);
      if (dist <= METAL_COORDINATION_CUTOFF) {
        const key = `${proteinAtom.chainId}:${proteinAtom.resSeq}`;
        const existing = residueMap.get(key);
        if (!existing || dist < existing.minDist) {
          residueMap.set(key, {
            chain: proteinAtom.chainId,
            residue: proteinAtom.resSeq,
            name: proteinAtom.resName,
            minDist: dist,
            ligandCode: metalAtom.resName,
            isMetal: true,
          });
        }
      }
    }
  }

  // Second pass: find ligand binding pocket residues (broader cutoff)
  for (const proteinAtom of chainProteinAtoms) {
    for (const ligandAtom of chainLigands) {
      const dist = distance(proteinAtom, ligandAtom);
      if (dist <= ligandCutoff) {
        const key = `${proteinAtom.chainId}:${proteinAtom.resSeq}`;
        const existing = residueMap.get(key);
        // Only add if not already a metal coordinator, or if this is closer
        if (!existing || (!existing.isMetal && dist < existing.minDist)) {
          // Don't overwrite metal coordinators with ligand contacts
          if (!existing?.isMetal) {
            residueMap.set(key, {
              chain: proteinAtom.chainId,
              residue: proteinAtom.resSeq,
              name: proteinAtom.resName,
              minDist: dist,
              ligandCode: ligandAtom.resName,
              isMetal: false,
            });
          }
        }
      }
    }
  }

  // Convert to array and sort by distance (closest first)
  const residues: CatalyticResidue[] = Array.from(residueMap.values())
    .sort((a, b) => a.minDist - b.minDist)
    .slice(0, 25)
    .map(r => {
      // Calculate confidence based on appropriate cutoff
      const effectiveCutoff = r.isMetal ? METAL_COORDINATION_CUTOFF : ligandCutoff;
      return {
        chain: r.chain,
        residue: r.residue,
        name: r.name,
        role: r.isMetal ? 'metal coordinator' : 'binding pocket',
        confidence: Math.max(0.5, 1.0 - (r.minDist / effectiveCutoff)),
        source: 'local' as const,
        ligandCode: r.ligandCode,
      };
    });

  console.log(`[Local] Found ${residues.length} binding pocket residues for chain ${targetChain}`);
  return residues;
}

// Primary catalytic roles that indicate true active site residues
const PRIMARY_CATALYTIC_ROLES = [
  'proton donor',
  'proton acceptor',
  'nucleophile',
  'electrophile',
  'nucleofuge',
  'electrofuge',
  'covalently attached',
  'activator',
  'metal ligand',
];

/**
 * Check if a role string contains any primary catalytic roles.
 */
function hasPrimaryCatalyticRole(roleStr: string | null): boolean {
  if (!roleStr) return false;
  const roleLower = roleStr.toLowerCase();
  return PRIMARY_CATALYTIC_ROLES.some(role => roleLower.includes(role));
}

/**
 * Query M-CSA for catalytic residues by PDB ID.
 * Filters to only include residues from the queried PDB with primary catalytic roles.
 */
async function queryMcsa(pdbId: string): Promise<CatalyticResidue[]> {
  if (!pdbId || pdbId.length !== 4) {
    return [];
  }

  const pdbIdUpper = pdbId.toUpperCase();
  const url = `${MCSA_API_BASE}/residues/?pdb_id=${pdbIdUpper}&format=json`;

  try {
    const response = await fetch(url, {
      headers: { 'Accept': 'application/json' },
      signal: AbortSignal.timeout(10000),
    });

    if (!response.ok) {
      return [];
    }

    const data = await response.json();

    if (!data || !Array.isArray(data)) {
      return [];
    }

    const residues: CatalyticResidue[] = [];
    const seen = new Set<string>();

    for (const entry of data) {
      // residue_chains is an array - look for chain matching our PDB ID
      const chains = entry.residue_chains;
      if (!chains || !Array.isArray(chains) || chains.length === 0) {
        continue;
      }

      // Find chain info that matches our specific PDB ID
      // M-CSA returns data from reference structures which may differ from queried PDB
      const matchingChain = chains.find(
        (c: { pdb_id?: string }) => c.pdb_id?.toUpperCase() === pdbIdUpper
      );

      // If no matching PDB, skip this residue (it's from a different reference structure)
      if (!matchingChain) {
        continue;
      }

      const chain = matchingChain.chain_name || 'A';
      const resid = matchingChain.resid;
      const code = matchingChain.code || 'UNK';
      const role = entry.roles_summary || null;

      if (resid === null || resid === undefined) {
        continue;
      }

      // Only include residues with primary catalytic roles
      if (!hasPrimaryCatalyticRole(role)) {
        continue;
      }

      const key = `${chain}${resid}`;
      if (seen.has(key)) {
        continue;
      }
      seen.add(key);

      residues.push({
        chain,
        residue: Number(resid),
        name: code,
        role: role || null,
        confidence: 1.0,
        source: 'mcsa',
      });
    }

    // Limit to top 25 residues (increased from 15 for better coverage)
    return residues.slice(0, 25);
  } catch (error) {
    console.error('[M-CSA] Query failed:', error);
    return [];
  }
}

/**
 * Extract PDB ID from HEADER line of PDB file.
 */
function extractPdbIdFromContent(pdbContent: string): string | null {
  if (!pdbContent) {
    return null;
  }

  const lines = pdbContent.split('\n').slice(0, 20);
  for (const line of lines) {
    if (line.startsWith('HEADER')) {
      // PDB ID is typically at positions 62-66
      if (line.length >= 66) {
        const pdbId = line.slice(62, 66).trim();
        if (pdbId.length === 4 && /^[A-Za-z0-9]{4}$/.test(pdbId)) {
          return pdbId.toUpperCase();
        }
      }
      // Also try regex for flexibility
      const match = line.match(/\b([0-9][A-Za-z0-9]{3})\s*$/);
      if (match) {
        return match[1].toUpperCase();
      }
    }
  }

  return null;
}

export async function POST(request: NextRequest): Promise<NextResponse<SuggestionResponse>> {
  try {
    const body: SuggestionRequest = await request.json();
    const { pdb_id, pdb_content, cutoff } = body;

    // Cutoff only applies to organic ligands; metals always use tight 3.5 Å for direct coordination
    const ligandCutoff = cutoff ?? DEFAULT_LIGAND_CUTOFF;
    console.log(`[catalytic-suggestions] Processing request (metals: ${METAL_COORDINATION_CUTOFF}Å, ligands: ${ligandCutoff}Å)...`);

    // Primary method: Local binding pocket detection with smart cutoffs
    // - Metals: tight 3.5 Å for direct coordinators only
    // - Ligands: broader cutoff for binding pocket
    const localResidues = findBindingPocketResidues(pdb_content, ligandCutoff);

    if (localResidues.length > 0) {
      console.log(`[catalytic-suggestions] Found ${localResidues.length} binding pocket residues via local detection`);
      return NextResponse.json({
        source: 'local',
        residues: localResidues,
      });
    }

    // Secondary: Try M-CSA for known enzymes (if we have a PDB ID)
    let pdbId = pdb_id;
    if (!pdbId) {
      pdbId = extractPdbIdFromContent(pdb_content);
    }

    if (pdbId) {
      console.log(`[catalytic-suggestions] Trying M-CSA for PDB ID: ${pdbId}`);
      const mcsaResidues = await queryMcsa(pdbId);

      if (mcsaResidues.length > 0) {
        console.log(`[catalytic-suggestions] Found ${mcsaResidues.length} catalytic residues from M-CSA`);
        return NextResponse.json({
          source: 'mcsa',
          residues: mcsaResidues,
        });
      }
    }

    // No suggestions found
    console.log('[catalytic-suggestions] No binding pocket residues found');
    return NextResponse.json({
      source: 'none',
      residues: [],
    });
  } catch (error) {
    console.error('[catalytic-suggestions] Error:', error);
    return NextResponse.json({
      source: 'none',
      residues: [],
    });
  }
}
