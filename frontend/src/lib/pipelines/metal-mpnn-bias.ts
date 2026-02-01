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
 * Global bias: A:-2.0 alanine suppression for ALL metal classes.
 * At low temperature (T=0.1) MPNN defaults to overproducing alanine (~60%).
 * The penalty reduces Ala to <10%.
 *
 * Coordination-specific biases (D/E/H/C positive boosts) are NOT applied
 * globally — they would push the entire protein toward those residues.
 * Instead, atom context + fixed_positions handle the binding site:
 *  - use_side_chain_context=true gives MPNN awareness of HETATM atoms
 *  - pack_with_ligand_context=true scores sidechains against ligand contacts
 *  - fixed_positions lock metal-coordinating residues (set in shared-steps.ts)
 *
 * Omit rules are intentionally conservative:
 *  - Hard acids omit only C (thiolate incompatible with hard Lewis acids).
 *  - Soft / borderline omit nothing — His, Gly, Pro are all valid in metal
 *    binding motifs (EF-hand loops depend on Gly; His coordinates even
 *    harder metals in mixed-donor sites).
 */
const CONFIG_BY_CLASS: Record<HsabClass, Pick<MetalMpnnConfig, 'bias_AA' | 'omit_AA'>> = {
  hard: {
    // Lanthanide / hard-acid: global alanine suppression only.
    // Atom context naturally places carboxylates at coordination sites
    // (ln_citrate: no donor bias + atom_context = 0.81Å RMSD).
    bias_AA: 'A:-2.0',
    omit_AA: 'C',
  },
  soft: {
    // Soft acid (Cu, Ag): global alanine suppression only.
    // Atom context + fixed_positions handle S/N-donor placement at binding site.
    bias_AA: 'A:-2.0',
  },
  borderline: {
    // Transition metals (Zn, Fe, Mn, Co, Ni): global alanine suppression only.
    // Atom context + fixed_positions handle mixed O/N/S-donor placement.
    bias_AA: 'A:-2.0',
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
