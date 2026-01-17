import { z } from 'zod';
import { QUALITY_PRESET_OPTIONS } from './metal-form';
// Re-export QualityPreset for consumers of this module
export type { QualityPreset } from './metal-form';

// Design approaches
export const LIGAND_APPROACHES = ['joint', 'asymmetric_rasa', 'induced', 'asymmetric'] as const;
export type LigandApproach = typeof LIGAND_APPROACHES[number];

// Side options for asymmetric approach
export const BINDING_SIDES = ['left', 'right'] as const;
export type BindingSide = typeof BINDING_SIDES[number];

// Quality presets imported from metal-form.ts to avoid duplicate exports

// Common ligands with SMILES
export const COMMON_LIGANDS_SMILES = [
  { id: 'azobenzene', name: 'Azobenzene', smiles: 'c1ccc(cc1)N=Nc2ccccc2', description: 'Light-switchable dye' },
  { id: 'caffeine', name: 'Caffeine', smiles: 'Cn1cnc2c1c(=O)n(c(=O)n2C)C', description: 'Stimulant alkaloid' },
  { id: 'benzene', name: 'Benzene', smiles: 'c1ccccc1', description: 'Simple aromatic ring' },
  { id: 'naphthalene', name: 'Naphthalene', smiles: 'c1ccc2ccccc2c1', description: 'Fused bicyclic aromatic' },
  { id: 'biphenyl', name: 'Biphenyl', smiles: 'c1ccc(cc1)c2ccccc2', description: 'Two connected phenyl rings' },
] as const;

// Main schema
export const ligandFormSchema = z.object({
  // Ligand Selection
  ligandSmiles: z.string().min(1, 'SMILES string is required'),
  selectedPreset: z.string(),
  customSmiles: z.boolean(),

  // Design Approach
  approach: z.enum(LIGAND_APPROACHES),
  side: z.enum(BINDING_SIDES),

  // Design Parameters
  chainLength: z.string().min(1, 'Chain length is required'),
  numDesigns: z.number().min(1).max(10),
  seed: z.string().optional(),

  // Quality Settings
  qualityPreset: z.enum(QUALITY_PRESET_OPTIONS),
  numTimesteps: z.number().min(10).max(500),
  stepScale: z.number().min(0.1).max(3.0),
  gamma0: z.number().min(0).max(1),

  // Advanced Options
  useOriToken: z.boolean(),
  oriOffset: z.string().optional(),
});

export type LigandFormValues = z.infer<typeof ligandFormSchema>;

export const ligandFormDefaults: LigandFormValues = {
  ligandSmiles: 'c1ccc(cc1)N=Nc2ccccc2', // Azobenzene
  selectedPreset: 'azobenzene',
  customSmiles: false,
  approach: 'joint',
  side: 'left',
  chainLength: '60-80',
  numDesigns: 3,
  seed: '',
  qualityPreset: 'Balanced',
  numTimesteps: 200,
  stepScale: 1.5,
  gamma0: 0.6,
  useOriToken: false,
  oriOffset: '12.0, 0.0, 0.0',
};

// Quality preset values are shared from metal-form.ts
// Import QUALITY_PRESET_VALUES from metal-form when needed
