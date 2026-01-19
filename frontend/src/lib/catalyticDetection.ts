/**
 * API client for catalytic residue detection.
 * Fetches suggestions from backend which queries M-CSA and P2Rank.
 */

import type { CatalyticSuggestion } from './store';

interface SuggestionResponse {
  source: 'mcsa' | 'p2rank' | 'none';
  residues: CatalyticSuggestion[];
}

/**
 * Fetch catalytic residue suggestions for a PDB structure.
 *
 * @param pdbContent - Full PDB file content
 * @param pdbId - Optional PDB ID if known (speeds up M-CSA lookup)
 * @param backendUrl - Backend API URL
 */
export async function fetchCatalyticSuggestions(
  pdbContent: string,
  pdbId: string | null,
  backendUrl: string
): Promise<SuggestionResponse> {
  const response = await fetch(`${backendUrl}/api/catalytic-suggestions`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      pdb_id: pdbId,
      pdb_content: pdbContent,
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
