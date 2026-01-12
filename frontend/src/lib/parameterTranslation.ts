/**
 * Parameter Translation Logic
 * Converts user-friendly interview preferences to RFD3 technical parameters
 */

import type { RFD3Request, UserPreferences, MetalBindingAnalysis } from './api';

// Metal-specific configuration
const METAL_CONFIG: Record<string, {
  typical_coordination: number[];
  preferred_donors: string[];
  bond_distance_range: [number, number];
  geometry: string;
}> = {
  TB: {
    typical_coordination: [8, 9],
    preferred_donors: ['O_carboxylate', 'O_carbonyl', 'O_water'],
    bond_distance_range: [2.3, 2.5],
    geometry: 'square_antiprism',
  },
  GD: {
    typical_coordination: [8, 9],
    preferred_donors: ['O_carboxylate', 'O_carbonyl', 'O_water'],
    bond_distance_range: [2.3, 2.5],
    geometry: 'square_antiprism',
  },
  EU: {
    typical_coordination: [8, 9],
    preferred_donors: ['O_carboxylate', 'O_carbonyl'],
    bond_distance_range: [2.3, 2.6],
    geometry: 'square_antiprism',
  },
  CA: {
    typical_coordination: [6, 8],
    preferred_donors: ['O_carboxylate', 'O_carbonyl', 'O_water'],
    bond_distance_range: [2.3, 2.6],
    geometry: 'octahedral',
  },
  ZN: {
    typical_coordination: [4, 6],
    preferred_donors: ['N_imidazole', 'S_thiolate', 'O_carboxylate'],
    bond_distance_range: [2.0, 2.3],
    geometry: 'tetrahedral',
  },
  FE: {
    typical_coordination: [4, 6],
    preferred_donors: ['S_thiolate', 'N_imidazole', 'O_carboxylate'],
    bond_distance_range: [2.0, 2.4],
    geometry: 'tetrahedral',
  },
};

// Structure info for contig generation
export interface StructureInfo {
  chains: string[];
  num_residues: number;
  num_atoms?: number;
}

/**
 * Translate user preferences from interview mode to RFD3 parameters
 */
export function translatePreferencesToParams(
  preferences: UserPreferences,
  analysis?: MetalBindingAnalysis | null,
  structureInfo?: StructureInfo | null
): Partial<RFD3Request> {
  const metalConfig = METAL_CONFIG[preferences.targetMetal] || METAL_CONFIG['TB'];

  const params: Partial<RFD3Request> = {
    ligand: preferences.targetMetal,
    num_designs: preferences.numDesigns,
  };

  // Generate contig for partial diffusion
  // For input structures, we need to specify the chain and residue range
  if (structureInfo && structureInfo.chains.length > 0 && structureInfo.num_residues > 0) {
    // Use the first chain and full residue range
    const chain = structureInfo.chains[0];
    params.contig = `${chain}1-${structureInfo.num_residues}`;
  } else {
    // Fallback: if no structure info, use a default length for de novo design
    params.length = '50-60';
  }

  // Design aggressiveness → partial_t and other noise parameters
  switch (preferences.designAggressiveness) {
    case 'conservative':
      params.partial_t = 8.0;
      params.num_timesteps = 150;
      params.step_scale = 1.2;
      params.gamma_0 = 0.8;
      break;
    case 'moderate':
      params.partial_t = 12.0;
      params.num_timesteps = 200;
      params.step_scale = 1.5;
      params.gamma_0 = 0.6;
      break;
    case 'aggressive':
      params.partial_t = 18.0;
      params.num_timesteps = 250;
      params.step_scale = 2.0;
      params.gamma_0 = 0.4;
      break;
  }

  // Coordination preference → additional constraints
  const isLanthanide = ['TB', 'GD', 'EU', 'LA', 'CE', 'ND'].includes(preferences.targetMetal);

  switch (preferences.coordinationPreference) {
    case 'tetrahedral':
      // Lower partial_t for compact sites
      params.partial_t = Math.min(params.partial_t || 12, 10);
      params.is_non_loopy = true;
      break;
    case 'octahedral':
      // Moderate settings
      params.is_non_loopy = true;
      break;
    case 'high':
      // Higher partial_t for larger pocket
      params.partial_t = Math.max(params.partial_t || 12, 15);
      params.step_scale = Math.max(params.step_scale || 1.5, 1.8);
      break;
    case 'auto':
      // Use metal-specific defaults
      if (isLanthanide) {
        params.partial_t = Math.max(params.partial_t || 12, 14);
        params.step_scale = 1.6;
      }
      break;
  }

  // Priority → quality settings optimization
  switch (preferences.priority) {
    case 'stability':
      params.num_timesteps = Math.max(params.num_timesteps || 200, 250);
      params.gamma_0 = Math.max(params.gamma_0 || 0.6, 0.7);
      params.is_non_loopy = true;
      break;
    case 'binding':
      params.step_scale = Math.max(params.step_scale || 1.5, 1.8);
      // Would add hotspot constraints if we have analysis data
      break;
    case 'expression':
      params.gamma_0 = Math.max(params.gamma_0 || 0.6, 0.8);
      params.is_non_loopy = true;
      break;
    case 'balanced':
      // Keep default balanced settings
      break;
  }

  // If we have analysis data, use it to inform parameters
  if (analysis?.success && analysis.coordination) {
    // Get coordinating residue positions for unindex
    const coordAtoms = analysis.coordination.coordinating_atoms;
    if (coordAtoms && coordAtoms.length > 0) {
      const positions = coordAtoms.map(atom =>
        `${atom.chain}${atom.residue_number}`
      );
      const uniquePositions = [...new Set(positions)];
      params.unindex = uniquePositions.join(',');

      // Set backbone-only for coordinating positions
      params.select_fixed_atoms = {};
      uniquePositions.forEach(pos => {
        params.select_fixed_atoms![pos] = 'BKBN';
      });
    }
  }

  // Add seed for reproducibility
  params.seed = Math.floor(Math.random() * 1000000);

  return params;
}

/**
 * Get human-readable explanation of the generated parameters
 */
export function explainParameters(params: Partial<RFD3Request>, preferences: UserPreferences): string[] {
  const explanations: string[] = [];

  explanations.push(`Target metal: ${preferences.targetMetalLabel}`);

  if (params.partial_t) {
    const level = params.partial_t < 10 ? 'conservative' : params.partial_t < 15 ? 'moderate' : 'aggressive';
    explanations.push(`Design approach: ${level} (partial diffusion level ${params.partial_t})`);
  }

  if (params.num_designs) {
    explanations.push(`Generating ${params.num_designs} candidate structures`);
  }

  if (params.contig) {
    explanations.push(`Using input structure (${params.contig})`);
  } else if (params.length) {
    explanations.push(`De novo design (${params.length} residues)`);
  }

  if (params.unindex) {
    const positions = params.unindex.split(',');
    explanations.push(`Redesigning ${positions.length} positions around the binding site`);
  }

  if (params.is_non_loopy) {
    explanations.push('Optimizing for clean secondary structures');
  }

  return explanations;
}

/**
 * Validate parameters are within safe ranges
 */
export function validateParameters(params: Partial<RFD3Request>): { valid: boolean; errors: string[] } {
  const errors: string[] = [];

  if (params.partial_t !== undefined) {
    if (params.partial_t < 1 || params.partial_t > 30) {
      errors.push('partial_t must be between 1 and 30');
    }
  }

  if (params.num_timesteps !== undefined) {
    if (params.num_timesteps < 50 || params.num_timesteps > 500) {
      errors.push('num_timesteps must be between 50 and 500');
    }
  }

  if (params.step_scale !== undefined) {
    if (params.step_scale < 0.5 || params.step_scale > 3.0) {
      errors.push('step_scale must be between 0.5 and 3.0');
    }
  }

  if (params.gamma_0 !== undefined) {
    if (params.gamma_0 < 0.1 || params.gamma_0 > 1.0) {
      errors.push('gamma_0 must be between 0.1 and 1.0');
    }
  }

  if (params.num_designs !== undefined) {
    if (params.num_designs < 1 || params.num_designs > 20) {
      errors.push('num_designs must be between 1 and 20');
    }
  }

  return {
    valid: errors.length === 0,
    errors,
  };
}
