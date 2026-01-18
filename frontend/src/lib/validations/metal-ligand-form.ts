import { z } from 'zod';
import { QUALITY_PRESET_OPTIONS } from './metal-form';
export type { QualityPreset } from './metal-form';

// Metal-ligand complex templates available in the backend
export const METAL_LIGAND_TEMPLATES = [
  {
    id: 'citrate_tb',
    name: 'Citrate-Terbium',
    description: 'Green luminescent biosensor',
    metal: 'TB',
    ligand: 'CIT',
    smiles: 'OC(CC(=O)[O-])(CC(=O)[O-])C(=O)[O-]',
    coordination: 9,
    ligandDonors: 3,
    proteinSitesNeeded: 6,
    luminescence: '545nm green',
  },
  {
    id: 'citrate_eu',
    name: 'Citrate-Europium',
    description: 'Red luminescent biosensor',
    metal: 'EU',
    ligand: 'CIT',
    smiles: 'OC(CC(=O)[O-])(CC(=O)[O-])C(=O)[O-]',
    coordination: 9,
    ligandDonors: 3,
    proteinSitesNeeded: 6,
    luminescence: '615nm red',
  },
  {
    id: 'citrate_gd',
    name: 'Citrate-Gadolinium',
    description: 'MRI contrast agent',
    metal: 'GD',
    ligand: 'CIT',
    smiles: 'OC(CC(=O)[O-])(CC(=O)[O-])C(=O)[O-]',
    coordination: 9,
    ligandDonors: 3,
    proteinSitesNeeded: 6,
    luminescence: 'UV/MRI',
  },
  {
    id: 'pqq_ca',
    name: 'PQQ-Calcium',
    description: 'Quinoprotein cofactor',
    metal: 'CA',
    ligand: 'PQQ',
    smiles: 'OC(=O)c1[nH]c2c(c1)C(=O)C(=O)c3nc(cc(C(O)=O)c23)C(O)=O',
    coordination: 6,
    ligandDonors: 3,
    proteinSitesNeeded: 3,
    luminescence: null,
  },
] as const;

export type MetalLigandTemplateId = typeof METAL_LIGAND_TEMPLATES[number]['id'];

// Design approaches for metal-ligand dimer
export const METAL_LIGAND_APPROACHES = ['joint', 'sequential'] as const;
export type MetalLigandApproach = typeof METAL_LIGAND_APPROACHES[number];

// Main schema
export const metalLigandFormSchema = z.object({
  // Template or Custom
  templateId: z.string(),
  useCustomComplex: z.boolean(),

  // Custom complex fields (when useCustomComplex is true)
  customSmiles: z.string().optional(),
  customMetal: z.string().optional(),
  customLigandName: z.string().optional(),
  customPdb: z.string().optional(),

  // Design Approach
  approach: z.enum(METAL_LIGAND_APPROACHES),

  // Design Parameters
  chainLength: z.string().min(1, 'Chain length is required'),
  numDesigns: z.number().min(1).max(10),
  seed: z.string().optional(),

  // Quality Settings
  qualityPreset: z.enum(QUALITY_PRESET_OPTIONS),
  numTimesteps: z.number().min(10).max(500),
  stepScale: z.number().min(0.1).max(3.0),
  gamma0: z.number().min(0).max(1),

  // Coordination Options
  targetCoordination: z.number().min(4).max(12).optional(),
  chainADonors: z.number().min(1).max(6),
  chainBDonors: z.number().min(1).max(6),

  // Validation Options
  validateCoordination: z.boolean(),
  runPipelineValidation: z.boolean(),

  // HSAB Bias Options
  useHsabBias: z.boolean(),
  customBiasAA: z.string().optional(),
});

export type MetalLigandFormValues = z.infer<typeof metalLigandFormSchema>;

export const metalLigandFormDefaults: MetalLigandFormValues = {
  templateId: 'citrate_tb',
  useCustomComplex: false,
  customSmiles: '',
  customMetal: 'TB',
  customLigandName: 'LIG',
  customPdb: '',
  approach: 'joint',
  chainLength: '60-80',
  numDesigns: 3,
  seed: '',
  qualityPreset: 'Balanced',
  numTimesteps: 200,
  stepScale: 1.5,
  gamma0: 0.6,
  targetCoordination: 9,
  chainADonors: 3,
  chainBDonors: 3,
  validateCoordination: true,
  runPipelineValidation: true,
  useHsabBias: true,
  customBiasAA: 'E:3.0,D:3.0,N:2.4,Q:2.4,H:0.5,C:-5.0',
};

// Helper to get template by ID
export function getTemplateById(id: string) {
  return METAL_LIGAND_TEMPLATES.find(t => t.id === id);
}

// Helper to get HSAB bias for a metal
export function getHsabBiasForMetal(metal: string): string {
  const upperMetal = metal.toUpperCase();

  // Lanthanides (hard acids) - favor hard bases (O donors)
  if (['TB', 'EU', 'GD', 'LA', 'CE', 'SM', 'YB'].includes(upperMetal)) {
    return 'E:3.0,D:3.0,N:2.4,Q:2.4,H:0.5,C:-5.0';
  }

  // Transition metals - varying preferences
  switch (upperMetal) {
    case 'ZN':
      return 'H:5.0,C:4.0,D:3.0,E:3.0';
    case 'FE':
      return 'H:5.0,C:3.0,D:2.0,E:2.0,Y:1.0';
    case 'CU':
      return 'H:5.0,C:4.0,M:3.0,D:2.0,E:2.0';
    case 'CA':
    case 'MG':
      return 'E:4.0,D:4.0,N:2.0,Q:2.0';
    default:
      return 'E:3.0,D:3.0,H:2.0,N:1.0';
  }
}
