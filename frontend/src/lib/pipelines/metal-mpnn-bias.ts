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

/**
 * Cofactor MPNN bias source: the backend's analyze_ligand_features endpoint
 * returns mpnn_bias_adjustments from the cofactor registry. The frontend
 * should pass that data through rather than duplicating it here.
 *
 * This is kept ONLY as a last-resort fallback when the backend response
 * is unavailable (e.g., offline mode or skipped coordination analysis).
 */
const COFACTOR_BIAS_FALLBACK: Record<string, string> = {
  PQQ: 'W:1.5,F:1.0,Y:1.0',
  HEM: 'H:2.0',
  FAD: 'W:1.0,Y:1.0',
  ATP: 'R:1.0,K:0.5',
};

/**
 * Merge two bias_AA strings, summing values for duplicate AAs.
 *
 * @example mergeBiasStrings('A:-2.0', 'W:1.5,F:1.0') => 'A:-2.0,F:1.0,W:1.5'
 */
function mergeBiasStrings(a: string | undefined, b: string): string {
  const parsed: Record<string, number> = {};
  for (const s of [a, b]) {
    if (!s) continue;
    for (const part of s.split(',')) {
      const [aa, val] = part.split(':');
      if (aa && val) {
        parsed[aa.trim()] = (parsed[aa.trim()] ?? 0) + parseFloat(val);
      }
    }
  }
  return Object.entries(parsed)
    .sort(([a], [b]) => a.localeCompare(b))
    .map(([aa, val]) => `${aa}:${val}`)
    .join(',');
}

/**
 * Return LigandMPNN config that merges HSAB metal bias with cofactor-specific
 * aromatic/coordination bias.
 *
 * @param targetMetal - Metal symbol (e.g., 'CA', 'TB')
 * @param ligandCode - Optional 3-letter ligand code (e.g., 'PQQ')
 * @param cofactorBias - Optional bias dict from backend (e.g., { W: 1.5, F: 1.0 }).
 *                       When provided, this takes priority over any local fallback.
 *
 * Falls back to plain HSAB config when no cofactor bias is available.
 */
export function getCofactorMpnnConfig(
  targetMetal: string,
  ligandCode?: string,
  cofactorBias?: Record<string, number>,
): MetalMpnnConfig {
  const base = getMetalMpnnConfig(targetMetal);
  if (!ligandCode) return base;

  // Primary: use backend-provided bias from analyze_ligand_features
  if (cofactorBias && Object.keys(cofactorBias).length > 0) {
    const biasStr = Object.entries(cofactorBias)
      .map(([aa, val]) => `${aa}:${val}`)
      .join(',');
    return { ...base, bias_AA: mergeBiasStrings(base.bias_AA, biasStr) };
  }

  // Fallback: local table (only when backend data unavailable)
  const extra = COFACTOR_BIAS_FALLBACK[ligandCode.toUpperCase()];
  if (!extra) return base;
  return { ...base, bias_AA: mergeBiasStrings(base.bias_AA, extra) };
}
