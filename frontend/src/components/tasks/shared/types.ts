/**
 * Shared types for RFD3 task forms
 */

// Quality presets for diffusion sampling
export const QUALITY_PRESETS = {
  Quick: {
    num_timesteps: 50,
    step_scale: 1.8,
    gamma_0: 0.5,
    description: 'Fast generation, good for exploration',
  },
  Balanced: {
    num_timesteps: 200,
    step_scale: 1.5,
    gamma_0: 0.6,
    description: 'Default settings, balanced quality',
  },
  'High Quality': {
    num_timesteps: 500,
    step_scale: 1.2,
    gamma_0: 0.7,
    description: 'Slower but higher quality designs',
  },
  'Binder Optimized': {
    num_timesteps: 200,
    step_scale: 3.0,
    gamma_0: 0.2,
    description: 'Recommended for protein binder design',
  },
} as const;

// Symmetry options
export const SYMMETRY_OPTIONS = [
  // Cyclic
  { id: 'C2', label: 'C2 (Dimer)', subunits: 2, category: 'Cyclic' },
  { id: 'C3', label: 'C3 (Trimer)', subunits: 3, category: 'Cyclic' },
  { id: 'C4', label: 'C4 (Tetramer)', subunits: 4, category: 'Cyclic' },
  { id: 'C5', label: 'C5 (Pentamer)', subunits: 5, category: 'Cyclic' },
  { id: 'C6', label: 'C6 (Hexamer)', subunits: 6, category: 'Cyclic' },
  // Dihedral
  { id: 'D2', label: 'D2', subunits: 4, category: 'Dihedral' },
  { id: 'D3', label: 'D3', subunits: 6, category: 'Dihedral' },
  { id: 'D4', label: 'D4', subunits: 8, category: 'Dihedral' },
  // Platonic
  { id: 'T', label: 'Tetrahedral', subunits: 12, category: 'Platonic' },
  { id: 'O', label: 'Octahedral', subunits: 24, category: 'Platonic' },
  { id: 'I', label: 'Icosahedral', subunits: 60, category: 'Platonic' },
] as const;

// Common ligands
// Codes are from PDB Chemical Component Dictionary: https://www.rcsb.org/chemical-component
export const COMMON_LIGANDS = [
  // Nucleotides
  { id: 'ATP', label: 'ATP', category: 'Nucleotide' },
  { id: 'ADP', label: 'ADP', category: 'Nucleotide' },
  { id: 'GTP', label: 'GTP', category: 'Nucleotide' },
  // Cofactors
  { id: 'NAD', label: 'NAD+', category: 'Cofactor' },
  { id: 'FAD', label: 'FAD', category: 'Cofactor' },
  { id: 'HEM', label: 'Heme', category: 'Cofactor' },
  { id: 'SAM', label: 'SAM', category: 'Cofactor' },
  { id: 'PLP', label: 'PLP', category: 'Cofactor' },
  // Common metals
  { id: 'ZN', label: 'Zinc', category: 'Metal' },
  { id: 'MG', label: 'Magnesium', category: 'Metal' },
  { id: 'CA', label: 'Calcium', category: 'Metal' },
  { id: 'FE', label: 'Iron', category: 'Metal' },
  { id: 'MN', label: 'Manganese', category: 'Metal' },
  { id: 'CO', label: 'Cobalt', category: 'Metal' },
  { id: 'CU', label: 'Copper', category: 'Metal' },
  { id: 'NI', label: 'Nickel', category: 'Metal' },
  // Lanthanides (rare earth metals)
  { id: 'LA', label: 'Lanthanum', category: 'Lanthanide' },
  { id: 'GD', label: 'Gadolinium', category: 'Lanthanide' },
  { id: 'EU', label: 'Europium', category: 'Lanthanide' },
  { id: 'TB', label: 'Terbium', category: 'Lanthanide' },
  { id: 'YB', label: 'Ytterbium', category: 'Lanthanide' },
  { id: 'LU', label: 'Lutetium', category: 'Lanthanide' },
] as const;

// Atom selection shortcuts
export const ATOM_SELECTION_OPTIONS = [
  { value: 'ALL', label: 'All atoms', description: 'All atoms in residue' },
  { value: 'BKBN', label: 'Backbone', description: 'N, CA, C, O atoms' },
  { value: 'TIP', label: 'Tip atom', description: 'Functional tip only' },
  { value: '', label: 'None (unfixed)', description: 'Allow diffusion' },
] as const;

// Base request type
export interface BaseRFD3Request {
  contig?: string;
  length?: string;
  num_designs?: number;
  seed?: number;
  pdb_content?: string;
  is_non_loopy?: boolean;
}

// Extended request type with all options
export interface RFD3Request extends BaseRFD3Request {
  // Diffusion parameters
  num_timesteps?: number;
  step_scale?: number;
  gamma_0?: number;

  // Motif specifications
  select_hotspots?: Record<string, string>;
  select_fixed_atoms?: Record<string, string>;
  select_unfixed_sequence?: Record<string, string>;

  // Ligand binding
  ligand?: string;
  select_buried?: Record<string, string>;
  select_exposed?: Record<string, string>;
  select_partially_buried?: Record<string, string>;

  // Nucleic acid binders
  ori_token?: [number, number, number];
  select_hbond_donor?: Record<string, string>;
  select_hbond_acceptor?: Record<string, string>;

  // Enzyme design
  unindex?: string;

  // Symmetry
  symmetry?: {
    id: string;
    is_unsym_motif?: string;
    is_symmetric_motif?: boolean;
  };

  // Partial diffusion
  partial_t?: number;

  // Binder design
  hotspots?: string[];  // Array of hotspot residues (e.g., ["A15", "A20"])
  infer_ori_strategy?: 'com' | 'hotspots';

  // Nucleic acid chains
  na_chains?: string;
}

// HealthResponse type (matches store)
export interface HealthResponse {
  status: string;
  gpu_available?: boolean;
  gpu_name?: string;
}

// Form submission handler type
export interface TaskFormProps {
  onSubmit: (request: RFD3Request) => Promise<void>;
  isSubmitting: boolean;
  health: HealthResponse | null;
}
