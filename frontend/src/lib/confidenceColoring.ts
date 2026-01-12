/**
 * Confidence-based Coloring for MolStar
 *
 * Provides utilities for coloring protein structures by per-residue pLDDT
 * using the AlphaFold color scheme.
 */

// AlphaFold-style pLDDT color scheme (matching ConfidenceMetrics.tsx)
export const PLDDT_COLORS = {
  veryHigh: '#0053d6',  // pLDDT >= 90 - blue
  high: '#65cbf3',       // pLDDT 70-90 - cyan
  medium: '#ffdb13',     // pLDDT 50-70 - yellow
  low: '#ff7d45'         // pLDDT < 50 - orange
} as const;

/**
 * Get color for a pLDDT value (0-1 scale)
 * @param value pLDDT value between 0 and 1
 * @returns hex color string
 */
export function getPLDDTColor(value: number): string {
  if (value >= 0.9) return PLDDT_COLORS.veryHigh;
  if (value >= 0.7) return PLDDT_COLORS.high;
  if (value >= 0.5) return PLDDT_COLORS.medium;
  return PLDDT_COLORS.low;
}

/**
 * Get color for a pLDDT value (0-100 scale)
 * @param value pLDDT value between 0 and 100
 * @returns hex color string
 */
export function getPLDDTColor100(value: number): string {
  return getPLDDTColor(value / 100);
}

/**
 * Get label for a pLDDT value
 * @param value pLDDT value between 0 and 1
 * @returns confidence level label
 */
export function getPLDDTLabel(value: number): string {
  if (value >= 0.9) return 'Very High';
  if (value >= 0.7) return 'High';
  if (value >= 0.5) return 'Medium';
  return 'Low';
}

/**
 * Convert hex color to RGB tuple
 * @param hex hex color string (e.g., '#0053d6')
 * @returns RGB tuple [r, g, b] with values 0-255
 */
export function hexToRgb(hex: string): [number, number, number] {
  const result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
  if (!result) return [128, 128, 128];
  return [
    parseInt(result[1], 16),
    parseInt(result[2], 16),
    parseInt(result[3], 16)
  ];
}

/**
 * Convert RGB to Mol* Color (0xRRGGBB format)
 * @param r red component 0-255
 * @param g green component 0-255
 * @param b blue component 0-255
 * @returns Mol* color number
 */
export function rgbToMolstarColor(r: number, g: number, b: number): number {
  return (r << 16) | (g << 8) | b;
}

/**
 * Get pLDDT color as Mol* Color number
 * @param value pLDDT value between 0 and 1
 * @returns Mol* color number
 */
export function getPLDDTMolstarColor(value: number): number {
  const hex = getPLDDTColor(value);
  const [r, g, b] = hexToRgb(hex);
  return rgbToMolstarColor(r, g, b);
}

// Pre-computed Mol* colors for each confidence level
export const PLDDT_MOLSTAR_COLORS = {
  veryHigh: rgbToMolstarColor(0, 83, 214),  // #0053d6
  high: rgbToMolstarColor(101, 203, 243),    // #65cbf3
  medium: rgbToMolstarColor(255, 219, 19),   // #ffdb13
  low: rgbToMolstarColor(255, 125, 69)       // #ff7d45
} as const;

/**
 * Generate a color map from per-residue pLDDT values
 * @param perResiduePLDDT array of pLDDT values (0-1 scale)
 * @returns Map of residue index to hex color
 */
export function generatePLDDTColorMap(perResiduePLDDT: number[]): Map<number, string> {
  const colorMap = new Map<number, string>();
  for (let i = 0; i < perResiduePLDDT.length; i++) {
    colorMap.set(i + 1, getPLDDTColor(perResiduePLDDT[i])); // 1-indexed residues
  }
  return colorMap;
}

/**
 * Statistics about confidence distribution
 */
export interface ConfidenceStats {
  mean: number;
  min: number;
  max: number;
  veryHighCount: number;  // >= 90
  highCount: number;       // 70-90
  mediumCount: number;     // 50-70
  lowCount: number;        // < 50
  veryHighPercent: number;
  highPercent: number;
  mediumPercent: number;
  lowPercent: number;
}

/**
 * Calculate statistics from per-residue pLDDT values
 * @param perResiduePLDDT array of pLDDT values (0-1 scale)
 * @returns confidence statistics
 */
export function calculateConfidenceStats(perResiduePLDDT: number[]): ConfidenceStats {
  if (perResiduePLDDT.length === 0) {
    return {
      mean: 0, min: 0, max: 0,
      veryHighCount: 0, highCount: 0, mediumCount: 0, lowCount: 0,
      veryHighPercent: 0, highPercent: 0, mediumPercent: 0, lowPercent: 0
    };
  }

  const n = perResiduePLDDT.length;
  let sum = 0;
  let min = 1;
  let max = 0;
  let veryHighCount = 0;
  let highCount = 0;
  let mediumCount = 0;
  let lowCount = 0;

  for (const value of perResiduePLDDT) {
    sum += value;
    if (value < min) min = value;
    if (value > max) max = value;

    if (value >= 0.9) veryHighCount++;
    else if (value >= 0.7) highCount++;
    else if (value >= 0.5) mediumCount++;
    else lowCount++;
  }

  return {
    mean: sum / n,
    min,
    max,
    veryHighCount,
    highCount,
    mediumCount,
    lowCount,
    veryHighPercent: (veryHighCount / n) * 100,
    highPercent: (highCount / n) * 100,
    mediumPercent: (mediumCount / n) * 100,
    lowPercent: (lowCount / n) * 100
  };
}

/**
 * Color legend data for UI display
 */
export const PLDDT_LEGEND = [
  { label: '>90', color: PLDDT_COLORS.veryHigh, description: 'Very High Confidence' },
  { label: '70-90', color: PLDDT_COLORS.high, description: 'High Confidence' },
  { label: '50-70', color: PLDDT_COLORS.medium, description: 'Medium Confidence' },
  { label: '<50', color: PLDDT_COLORS.low, description: 'Low Confidence' }
] as const;

/**
 * RMSD-based deviation coloring for structure comparison
 */
export const RMSD_COLORS = {
  identical: '#22C55E',   // < 0.5A - green
  similar: '#3B82F6',     // 0.5-1.0A - blue
  moderate: '#F59E0B',    // 1.0-2.0A - amber
  different: '#EF4444'    // > 2.0A - red
} as const;

/**
 * Get color for per-residue RMSD deviation
 * @param rmsd RMSD value in Angstroms
 * @returns hex color string
 */
export function getRMSDDeviationColor(rmsd: number): string {
  if (rmsd < 0.5) return RMSD_COLORS.identical;
  if (rmsd < 1.0) return RMSD_COLORS.similar;
  if (rmsd < 2.0) return RMSD_COLORS.moderate;
  return RMSD_COLORS.different;
}

/**
 * Color legend for RMSD deviation
 */
export const RMSD_LEGEND = [
  { label: '<0.5A', color: RMSD_COLORS.identical, description: 'Nearly Identical' },
  { label: '0.5-1.0A', color: RMSD_COLORS.similar, description: 'Similar' },
  { label: '1.0-2.0A', color: RMSD_COLORS.moderate, description: 'Moderate Deviation' },
  { label: '>2.0A', color: RMSD_COLORS.different, description: 'Significant Deviation' }
] as const;
