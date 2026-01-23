/**
 * API client for catalytic residue detection.
 * Fetches suggestions from local Next.js API.
 * Primary: Local binding pocket detection (finds residues near ligands/metals)
 * Secondary: M-CSA for known enzymes
 */

import type { CatalyticSuggestion } from './store';

interface SuggestionResponse {
  source: 'mcsa' | 'local' | 'none';
  residues: CatalyticSuggestion[];
}

/**
 * Fetch catalytic residue suggestions for a PDB structure.
 * Uses SMART detection: metals get tight 3.5 Å cutoff for direct coordinators,
 * organic ligands use the provided cutoff for binding pocket residues.
 *
 * @param pdbContent - Full PDB file content
 * @param pdbId - Optional PDB ID if known (speeds up M-CSA lookup)
 * @param _backendUrl - Unused, kept for API compatibility
 * @param cutoff - Distance cutoff for organic ligands only (default: 5.0 Å). Metals always use 3.5 Å.
 */
export async function fetchCatalyticSuggestions(
  pdbContent: string,
  pdbId: string | null,
  _backendUrl: string,
  cutoff?: number
): Promise<SuggestionResponse> {
  // Use local Next.js API route instead of Python backend
  const response = await fetch('/api/catalytic-suggestions', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      pdb_id: pdbId,
      pdb_content: pdbContent,
      cutoff,
    }),
  });

  if (!response.ok) {
    throw new Error(`Failed to fetch suggestions: ${response.statusText}`);
  }

  return response.json();
}

/**
 * Extract PDB ID from structure content.
 * Checks HEADER line for 4-character PDB code.
 */
export function extractPdbId(pdbContent: string): string | null {
  if (!pdbContent) return null;

  const lines = pdbContent.split('\n').slice(0, 20);
  for (const line of lines) {
    if (line.startsWith('HEADER')) {
      // PDB ID at positions 62-66
      const pdbId = line.slice(62, 66).trim();
      if (pdbId.length === 4 && /^[0-9A-Z]{4}$/i.test(pdbId)) {
        return pdbId.toUpperCase();
      }
      // Try regex fallback
      const match = line.match(/\b([0-9][A-Za-z0-9]{3})\s*$/);
      if (match) {
        return match[1].toUpperCase();
      }
    }
  }

  return null;
}

/**
 * Get list of protein chains in PDB content.
 */
export function getProteinChains(pdbContent: string): string[] {
  const chains = new Set<string>();
  const lines = pdbContent.split('\n');

  for (const line of lines) {
    if (line.startsWith('ATOM  ')) {
      const chainId = line.charAt(21);
      if (chainId && chainId !== ' ') {
        chains.add(chainId);
      }
    }
  }

  return [...chains].sort();
}

/**
 * Filter PDB content to only include one biological assembly (single chain and its ligands).
 * This is necessary for structures with multiple copies to avoid duplicate ligands.
 *
 * @param pdbContent - Full PDB file content
 * @param targetChain - Chain to keep (defaults to first chain alphabetically)
 * @returns Filtered PDB content with only the target chain and associated ligands
 */
export function filterPdbToSingleAssembly(pdbContent: string, targetChain?: string): string {
  const lines = pdbContent.split('\n');

  // Find all protein chains
  const proteinChains = new Set<string>();
  for (const line of lines) {
    if (line.startsWith('ATOM  ')) {
      const chainId = line.charAt(21);
      if (chainId && chainId !== ' ') {
        proteinChains.add(chainId);
      }
    }
  }

  // If only one chain or no chains, return original content
  if (proteinChains.size <= 1) {
    return pdbContent;
  }

  // Use target chain or first alphabetically
  const chain = targetChain || [...proteinChains].sort()[0];
  console.log(`[PDB Filter] Multiple chains found: ${[...proteinChains].join(', ')}. Filtering to chain ${chain}`);

  // Parse protein atom coordinates for distance calculation
  interface Coord { x: number; y: number; z: number }
  const chainAtoms: Coord[] = [];

  for (const line of lines) {
    if (line.startsWith('ATOM  ') && line.charAt(21) === chain) {
      const x = parseFloat(line.slice(30, 38));
      const y = parseFloat(line.slice(38, 46));
      const z = parseFloat(line.slice(46, 54));
      if (!isNaN(x) && !isNaN(y) && !isNaN(z)) {
        chainAtoms.push({ x, y, z });
      }
    }
  }

  // Function to check if a HETATM is close to our target chain
  const isCloseToChain = (line: string): boolean => {
    const x = parseFloat(line.slice(30, 38));
    const y = parseFloat(line.slice(38, 46));
    const z = parseFloat(line.slice(46, 54));
    if (isNaN(x) || isNaN(y) || isNaN(z)) return false;

    // Check if within 15Å of any chain atom (generous cutoff for ligand association)
    const cutoff = 15.0;
    const cutoffSq = cutoff * cutoff;

    for (const atom of chainAtoms) {
      const dx = x - atom.x;
      const dy = y - atom.y;
      const dz = z - atom.z;
      if (dx * dx + dy * dy + dz * dz <= cutoffSq) {
        return true;
      }
    }
    return false;
  };

  // Track which HETATM residues (resSeq only, for ligands without chain ID) are close to our chain
  const closeHetatmResidues = new Set<string>();
  for (const line of lines) {
    if (line.startsWith('HETATM')) {
      const hetatmChain = line.charAt(21);
      // Only track distance for ligands without explicit chain ID
      if ((!hetatmChain || hetatmChain === ' ') && isCloseToChain(line)) {
        const resSeq = line.slice(22, 26).trim();
        closeHetatmResidues.add(resSeq);
      }
    }
  }

  // Filter lines
  const filteredLines: string[] = [];
  for (const line of lines) {
    const recordType = line.slice(0, 6).trim();

    if (recordType === 'ATOM') {
      // Only keep ATOM records from target chain
      if (line.charAt(21) === chain) {
        filteredLines.push(line);
      }
    } else if (recordType === 'HETATM') {
      const hetatmChain = line.charAt(21);
      const resSeq = line.slice(22, 26).trim();

      // Keep HETATM if:
      // 1. Chain ID matches target chain, OR
      // 2. No chain ID and is close to target chain (distance-based)
      if (hetatmChain === chain) {
        // Explicit chain match - keep it
        filteredLines.push(line);
      } else if ((!hetatmChain || hetatmChain === ' ') && closeHetatmResidues.has(resSeq)) {
        // No chain ID but close to target chain - keep it
        filteredLines.push(line);
      }
      // Otherwise skip (belongs to different chain/assembly)
    } else if (recordType === 'TER') {
      // Keep TER only for our chain
      if (line.charAt(21) === chain || line.charAt(21) === ' ') {
        filteredLines.push(line);
      }
    } else if (recordType === 'CONECT') {
      // Skip CONECT records - they reference atom serial numbers that may not exist
      // after filtering and cause parser errors (IndexError in biotite)
    } else {
      // Keep all other records (HEADER, REMARK, etc.)
      filteredLines.push(line);
    }
  }

  console.log(`[PDB Filter] Filtered from ${lines.length} to ${filteredLines.length} lines`);
  return filteredLines.join('\n');
}
