import { z } from 'zod';

export const metalFormSchema = z.object({
  // Metal Selection
  metal: z.enum([
    'Cu', 'Zn', 'Fe', 'Mn', 'Co', 'Ni', 'Ca', 'Mg',
    'Na', 'K', 'Cd', 'Hg', 'Pb', 'Mo', 'W'
  ]),
  coordination: z.number().min(2).max(8),
  geometry: z.enum([
    'tetrahedral',
    'square_planar',
    'octahedral',
    'trigonal_bipyramidal',
    'square_pyramidal',
    'linear',
    'trigonal_planar'
  ]).optional(),

  // Binding Site
  bindingSiteResidues: z.string().optional(),
  distanceConstraint: z.number().min(1.5).max(4.0).default(2.3),

  // Quality Settings
  qualityPreset: z.enum(['Quick', 'Balanced', 'High', 'Custom']).default('Balanced'),
  numTimesteps: z.number().min(10).max(500).default(50),
  stepScale: z.number().min(0.1).max(2.0).default(1.0),
  gamma0: z.number().min(0).max(1).default(0.0),

  // Output
  numDesigns: z.number().min(1).max(10).default(4),

  // Advanced
  usePotentials: z.boolean().default(true),
  potentialScale: z.number().min(0.1).max(10).default(1.0),
});

export type MetalFormValues = z.infer<typeof metalFormSchema>;

export const metalFormDefaults: MetalFormValues = {
  metal: 'Cu',
  coordination: 4,
  geometry: 'tetrahedral',
  bindingSiteResidues: '',
  distanceConstraint: 2.3,
  qualityPreset: 'Balanced',
  numTimesteps: 50,
  stepScale: 1.0,
  gamma0: 0.0,
  numDesigns: 4,
  usePotentials: true,
  potentialScale: 1.0,
};
