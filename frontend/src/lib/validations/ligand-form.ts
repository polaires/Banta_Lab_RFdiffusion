import { z } from 'zod';

export const ligandFormSchema = z.object({
  // Ligand Selection
  ligand: z.string().min(1, 'Ligand is required'),
  customSmiles: z.string().optional(),

  // Binding Site
  bindingPocket: z.enum(['auto', 'manual']).default('auto'),
  pocketResidues: z.string().optional(),
  pocketRadius: z.number().min(3).max(15).default(6),

  // Design Parameters
  scaffoldLength: z.object({
    min: z.number().min(30).max(500),
    max: z.number().min(30).max(500),
  }).refine(data => data.max >= data.min, {
    message: 'Max must be greater than or equal to min',
  }),

  // Quality Settings
  qualityPreset: z.enum(['Quick', 'Balanced', 'High', 'Custom']).default('Balanced'),
  numTimesteps: z.number().min(10).max(500).default(50),
  stepScale: z.number().min(0.1).max(2.0).default(1.0),
  gamma0: z.number().min(0).max(1).default(0.0),

  // Output
  numDesigns: z.number().min(1).max(10).default(4),

  // Advanced
  useConstraints: z.boolean().default(true),
  constraintWeight: z.number().min(0.1).max(10).default(1.0),
});

export type LigandFormValues = z.infer<typeof ligandFormSchema>;

export const ligandFormDefaults: LigandFormValues = {
  ligand: 'ATP',
  customSmiles: '',
  bindingPocket: 'auto',
  pocketResidues: '',
  pocketRadius: 6,
  scaffoldLength: { min: 80, max: 150 },
  qualityPreset: 'Balanced',
  numTimesteps: 50,
  stepScale: 1.0,
  gamma0: 0.0,
  numDesigns: 4,
  useConstraints: true,
  constraintWeight: 1.0,
};
