/**
 * Metal-aware LigandMPNN configuration.
 *
 * Single source of truth for HSAB-based bias, omit, and context flags.
 * Callers spread the returned object directly into a ProteinMPNNRequest —
 * no manual flag assembly needed elsewhere.
 *
 * Bias values aligned with the LigandMPNN reference skill and
 * Nature Methods 2025 best-practice recommendations.
 */

import type { ProteinMPNNRequest } from '@/lib/api';

/** Fields this utility may set on a ProteinMPNNRequest. */
export type MetalMpnnConfig = Pick<
  ProteinMPNNRequest,
  | 'bias_AA'
  | 'omit_AA'
  | 'pack_side_chains'
  | 'pack_with_ligand_context'
  | 'use_side_chain_context'
  | 'fixed_positions'
>;

type HsabClass = 'hard' | 'soft' | 'borderline';

const HSAB_CLASSIFICATION: Record<string, HsabClass> = {
  // Hard acids – prefer O-donors (Asp/Glu/Asn/Gln)
  TB: 'hard', EU: 'hard', GD: 'hard', CA: 'hard', MG: 'hard',
  LA: 'hard', CE: 'hard', SM: 'hard', YB: 'hard', DY: 'hard',
  // Soft acids – prefer S/N-donors (Cys/His/Met)
  CU: 'soft', AG: 'soft',
  // Borderline – mixed donor preference
  ZN: 'borderline', FE: 'borderline', MN: 'borderline',
  CO: 'borderline', NI: 'borderline',
};

/**
 * Per-class bias & omit rules.
 *
 * Hard acids: NO bias — atom context naturally places carboxylates at
 * coordination sites (ln_citrate experiments: no bias + atom_context = 0.81Å
 * RMSD vs E:3.0,D:3.0 bias = 21-25Å RMSD / non-foldable).
 *
 * Omit rules are intentionally conservative:
 *  - Hard acids omit only C (thiolate incompatible with hard Lewis acids).
 *  - Soft / borderline omit nothing — His, Gly, Pro are all valid in metal
 *    binding motifs (EF-hand loops depend on Gly; His coordinates even
 *    harder metals in mixed-donor sites).
 */
const CONFIG_BY_CLASS: Record<HsabClass, Pick<MetalMpnnConfig, 'bias_AA' | 'omit_AA'>> = {
  hard: {
    // Lanthanide / hard-acid: NO bias — atom context handles coordination
    // Omit C only (thiolate incompatible)
    omit_AA: 'C',
  },
  soft: {
    // Soft acid (Cu, Ag): S/N donors
    bias_AA: 'C:4.0,H:3.0,M:2.0,A:-1.5',
  },
  borderline: {
    // Transition metals (Zn, Fe, Mn, Co, Ni): mixed O/N/S donors
    bias_AA: 'H:3.0,D:2.0,E:2.0,C:1.5,A:-1.5',
  },
};

/**
 * Return a complete LigandMPNN configuration for a given target metal.
 *
 * The result can be spread directly into a `ProteinMPNNRequest`:
 * ```ts
 * api.submitMPNNDesign({ pdb_content, ...getMetalMpnnConfig('TB') })
 * ```
 *
 * Falls back to borderline (mixed) config for unknown metals.
 */
export function getMetalMpnnConfig(targetMetal: string): MetalMpnnConfig {
  const key = targetMetal.toUpperCase().trim();
  const hsab = HSAB_CLASSIFICATION[key] ?? 'borderline';
  const { bias_AA, omit_AA } = CONFIG_BY_CLASS[hsab];

  return {
    bias_AA,
    omit_AA,
    pack_side_chains: true,
    pack_with_ligand_context: true,
    use_side_chain_context: true,
  };
}

/** Expose HSAB class for display / logging purposes. */
export function getHsabClass(targetMetal: string): HsabClass {
  const key = targetMetal.toUpperCase().trim();
  return HSAB_CLASSIFICATION[key] ?? 'borderline';
}
