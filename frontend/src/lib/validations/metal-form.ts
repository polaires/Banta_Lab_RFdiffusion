import { z } from 'zod';

// Metal ion identifiers
export const METAL_IDS = [
  'ZN', 'FE', 'CU', 'MN', 'CA', 'MG', 'TB', 'GD', 'EU'
] as const;

export type MetalId = typeof METAL_IDS[number];

// Design approaches for metal dimer
export const METAL_APPROACHES = [
  'joint_metal',
  'asymmetric_metal',
  'induced_metal',
  'bridging_metal',
  'redox_switch',
  'lanthanide_sensor',
] as const;

export type MetalApproach = typeof METAL_APPROACHES[number];

// Quality presets
export const QUALITY_PRESET_OPTIONS = ['Quick', 'Balanced', 'High Quality', 'Binder Optimized', 'Custom'] as const;
export type QualityPreset = typeof QUALITY_PRESET_OPTIONS[number];

// Geometry options
export const GEOMETRY_OPTIONS = [
  'auto',
  'tetrahedral',
  'square_planar',
  'trigonal_bipyramidal',
  'square_pyramidal',
  'octahedral',
  'pentagonal_bipyramidal',
  'capped_octahedral',
  'square_antiprismatic',
  'dodecahedral',
  'tricapped_trigonal_prismatic',
] as const;

export type Geometry = typeof GEOMETRY_OPTIONS[number];

// Donor types
export const DONOR_TYPES = ['His', 'Cys', 'Asp', 'Glu', 'Met', 'Tyr', 'Asn', 'Water', 'Backbone O'] as const;
export type DonorType = typeof DONOR_TYPES[number];

// Template library options
export const TEMPLATE_LIBRARY_OPTIONS = ['caldwell_4', 'ef_hand_8', 'lanm_mixed', 'high_coord_9', 'legacy'] as const;
export type TemplateLibraryId = typeof TEMPLATE_LIBRARY_OPTIONS[number];

// Legacy template types
export const LEGACY_TEMPLATE_TYPES = ['ef_hand', 'c4_symmetric', 'none'] as const;
export type LegacyTemplateType = typeof LEGACY_TEMPLATE_TYPES[number];

// Legacy donor residues
export const LEGACY_DONOR_RESIDUES = ['ASP', 'GLU', 'MIXED'] as const;
export type LegacyDonorResidue = typeof LEGACY_DONOR_RESIDUES[number];

// Main schema
export const metalFormSchema = z.object({
  // Metal Selection
  metal: z.enum(METAL_IDS),
  approach: z.enum(METAL_APPROACHES),

  // Coordination Split
  chainADonors: z.number().min(1).max(6),
  chainBDonors: z.number().min(1).max(6),
  chainADonorTypes: z.array(z.enum(DONOR_TYPES)),
  chainBDonorTypes: z.array(z.enum(DONOR_TYPES)),

  // Design Parameters
  chainLength: z.string().min(1, 'Chain length is required'),
  numDesigns: z.number().min(1).max(20),
  seed: z.string().optional(),

  // Quality Settings
  qualityPreset: z.enum(QUALITY_PRESET_OPTIONS),
  numTimesteps: z.number().min(10).max(500),
  stepScale: z.number().min(0.1).max(3.0),
  gamma0: z.number().min(0).max(1),

  // Advanced Options
  targetGeometry: z.enum(GEOMETRY_OPTIONS),
  includeWaters: z.boolean(),
  secondMetal: z.string().optional(),

  // Lanthanide Template Options (mode selection)
  useParametricMode: z.boolean(),
  useLegacyMode: z.boolean(),

  // Template Library Mode
  templateName: z.enum(TEMPLATE_LIBRARY_OPTIONS),

  // Parametric Mode
  paramCoordinationNumber: z.number().min(6).max(10),
  paramNumWaters: z.number().min(0).max(4),
  paramBidentateFraction: z.number().min(0).max(1),

  // Legacy Mode
  templateType: z.enum(LEGACY_TEMPLATE_TYPES),
  donorResidue: z.enum(LEGACY_DONOR_RESIDUES),
  useMotifScaffolding: z.boolean(),

  // Luminescence (for TB, EU)
  addTrpAntenna: z.boolean(),
  validateCoordination: z.boolean(),
});

export type MetalFormValues = z.infer<typeof metalFormSchema>;

export const metalFormDefaults: MetalFormValues = {
  metal: 'ZN',
  approach: 'joint_metal',
  chainADonors: 2,
  chainBDonors: 2,
  chainADonorTypes: ['His', 'His'],
  chainBDonorTypes: ['Cys', 'Cys'],
  chainLength: '60-80',
  numDesigns: 5,
  seed: '',
  qualityPreset: 'Balanced',
  numTimesteps: 200,
  stepScale: 1.5,
  gamma0: 0.6,
  targetGeometry: 'auto',
  includeWaters: false,
  secondMetal: '',
  useParametricMode: false,
  useLegacyMode: false,
  templateName: 'caldwell_4',
  paramCoordinationNumber: 8,
  paramNumWaters: 0,
  paramBidentateFraction: 0.5,
  templateType: 'ef_hand',
  donorResidue: 'ASP',
  useMotifScaffolding: true,
  addTrpAntenna: false,
  validateCoordination: true,
};

// Quality preset values
export const QUALITY_PRESET_VALUES = {
  Quick: { numTimesteps: 50, stepScale: 1.8, gamma0: 0.5 },
  Balanced: { numTimesteps: 200, stepScale: 1.5, gamma0: 0.6 },
  'High Quality': { numTimesteps: 500, stepScale: 1.2, gamma0: 0.7 },
  'Binder Optimized': { numTimesteps: 200, stepScale: 3.0, gamma0: 0.2 },
  Custom: { numTimesteps: 200, stepScale: 1.5, gamma0: 0.6 },
} as const;
