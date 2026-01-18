/**
 * Reusable MolScript Expressions
 * Common selection patterns for structure queries
 */

import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder';

/**
 * Select water molecules (HOH or WAT)
 */
export const WATER_EXPRESSION = MS.struct.generator.atomGroups({
  'residue-test': MS.core.logic.or([
    MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_comp_id(), 'HOH']),
    MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_comp_id(), 'WAT'])
  ])
});

/**
 * Select polymer entities (protein, DNA, RNA)
 */
export const POLYMER_EXPRESSION = MS.struct.generator.atomGroups({
  'entity-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.entityType(), 'polymer'])
});

/**
 * Select non-polymer entities (ligands, ions, etc.)
 */
export const NON_POLYMER_EXPRESSION = MS.struct.generator.atomGroups({
  'entity-test': MS.core.rel.neq([MS.struct.atomProperty.macromolecular.entityType(), 'polymer'])
});

/**
 * Create expression for a specific residue by chain and sequence number
 */
export function residueExpression(
  resName: string,
  chainId: string,
  resSeq: number
) {
  return MS.struct.generator.atomGroups({
    'residue-test': MS.core.logic.and([
      MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_comp_id(), resName]),
      MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_asym_id(), chainId]),
      MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_seq_id(), resSeq]),
    ])
  });
}

/**
 * Create surrounding sphere expression
 */
export function surroundingsExpression(
  centerExpression: ReturnType<typeof MS.struct.generator.atomGroups>,
  radius: number,
  asWholeResidues = true
) {
  return MS.struct.modifier.includeSurroundings({
    0: centerExpression,
    radius,
    'as-whole-residues': asWholeResidues
  });
}

/**
 * Exclude expression B from expression A
 */
export function exceptExpression(
  baseExpression: ReturnType<typeof MS.struct.generator.atomGroups>,
  excludeExpression: ReturnType<typeof MS.struct.generator.atomGroups>
) {
  return MS.struct.modifier.exceptBy({
    0: baseExpression,
    by: excludeExpression
  });
}

/**
 * Intersect expression A with expression B
 */
export function intersectExpression(
  baseExpression: ReturnType<typeof MS.struct.generator.atomGroups>,
  filterExpression: ReturnType<typeof MS.struct.generator.atomGroups>
) {
  return MS.struct.modifier.intersectBy({
    0: baseExpression,
    by: filterExpression
  });
}
