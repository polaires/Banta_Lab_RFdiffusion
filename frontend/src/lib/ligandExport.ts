/**
 * Ligand Analysis Export Utilities
 *
 * Provides export functionality for ligand analysis results:
 * - JSON export (complete data)
 * - CSV export (contacts table)
 * - RFdiffusion hotspot configuration export
 */

import type { LigandAnalysisResult, LigandData, LigandContact } from './ligandAnalysis';

/**
 * Export complete analysis as JSON
 *
 * Returns a formatted JSON string with all ligand analysis data
 */
export function exportToJSON(data: LigandAnalysisResult): string {
  return JSON.stringify(data, null, 2);
}

/**
 * Export contacts as CSV
 *
 * Headers: Ligand,Chain,ResSeq,ContactResidue,ContactChain,ContactAtom,Distance,InteractionType,BindingSiteType
 * One row per contact, not per ligand
 */
export function exportToCSV(data: LigandAnalysisResult): string {
  const headers = [
    'Ligand',
    'Chain',
    'ResSeq',
    'ContactResidue',
    'ContactChain',
    'ContactAtom',
    'Distance',
    'InteractionType',
    'BindingSiteType',
  ];

  const rows: string[] = [headers.join(',')];

  for (const ligand of data.ligandDetails) {
    for (const contact of ligand.contacts) {
      const row = [
        ligand.name,
        ligand.chainId,
        ligand.resSeq.toString(),
        contact.residue,
        contact.chain,
        contact.atom,
        contact.distance.toFixed(2),
        contact.interactionType,
        ligand.bindingSiteType,
      ];
      rows.push(row.join(','));
    }
  }

  return rows.join('\n');
}

/**
 * Extract residue number from a contact residue string
 *
 * Examples:
 * - "ASP45" -> 45
 * - "PHE101" -> 101
 * - "LYS9" -> 9
 */
function extractResidueNumber(residue: string): number | null {
  const match = residue.match(/(\d+)$/);
  if (match) {
    return parseInt(match[1], 10);
  }
  return null;
}

/**
 * Format residue as RFdiffusion hotspot format
 *
 * Examples:
 * - Chain "A", residue "ASP45" -> "A45"
 * - Chain "B", residue "PHE101" -> "B101"
 */
function formatHotspotResidue(chain: string, residue: string): string | null {
  const resNum = extractResidueNumber(residue);
  if (resNum === null) {
    return null;
  }
  return `${chain}${resNum}`;
}

/**
 * Export as RFdiffusion hotspot configuration
 *
 * Generates a configuration snippet for RFdiffusion with hotspot residues
 * Only includes contacts from functional binding sites
 * Format: hotspot_residues: [A45, A72, B89, ...]
 */
export function exportToRFdiffusionHotspots(data: LigandAnalysisResult): string {
  const hotspots = new Set<string>();

  // Only include functional binding sites
  for (const ligand of data.ligandDetails) {
    if (ligand.bindingSiteType !== 'functional') {
      continue;
    }

    for (const contact of ligand.contacts) {
      const hotspot = formatHotspotResidue(contact.chain, contact.residue);
      if (hotspot) {
        hotspots.add(hotspot);
      }
    }
  }

  // Sort hotspots by chain then by residue number for consistent output
  const sortedHotspots = Array.from(hotspots).sort((a, b) => {
    const chainA = a.charAt(0);
    const chainB = b.charAt(0);
    if (chainA !== chainB) {
      return chainA.localeCompare(chainB);
    }
    const numA = parseInt(a.slice(1), 10);
    const numB = parseInt(b.slice(1), 10);
    return numA - numB;
  });

  const hotspotsArray = sortedHotspots.join(', ');

  const lines = [
    '# RFdiffusion Hotspot Configuration',
    '# Generated from ligand contact analysis',
    '#',
    '# These residues form functional binding sites with ligands.',
    '# Use these as hotspot_residues to guide binder design toward',
    '# these contact positions.',
    '#',
    '# Usage in RFdiffusion:',
    '#   --potentials.guide_potentials "type:custom_hotspots,pdb_path:<target.pdb>,hotspot_residues:[' + hotspotsArray + ']"',
    '#',
    `# Total hotspot residues: ${sortedHotspots.length}`,
    '',
    `hotspot_residues: [${hotspotsArray}]`,
  ];

  return lines.join('\n');
}

/**
 * Trigger browser download
 *
 * Creates a blob from the content and triggers a download with the specified filename
 */
export function downloadFile(content: string, filename: string, mimeType: string): void {
  const blob = new Blob([content], { type: mimeType });
  const url = URL.createObjectURL(blob);

  const link = document.createElement('a');
  link.href = url;
  link.download = filename;

  // Append to body, click, and remove
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);

  // Clean up the URL object
  URL.revokeObjectURL(url);
}

/**
 * Download ligand analysis as JSON
 */
export function downloadJSON(data: LigandAnalysisResult, filename: string = 'ligand-analysis.json'): void {
  const content = exportToJSON(data);
  downloadFile(content, filename, 'application/json');
}

/**
 * Download ligand contacts as CSV
 */
export function downloadCSV(data: LigandAnalysisResult, filename: string = 'ligand-contacts.csv'): void {
  const content = exportToCSV(data);
  downloadFile(content, filename, 'text/csv');
}

/**
 * Download RFdiffusion hotspot configuration
 */
export function downloadRFdiffusionHotspots(
  data: LigandAnalysisResult,
  filename: string = 'rfdiffusion-hotspots.txt'
): void {
  const content = exportToRFdiffusionHotspots(data);
  downloadFile(content, filename, 'text/plain');
}
