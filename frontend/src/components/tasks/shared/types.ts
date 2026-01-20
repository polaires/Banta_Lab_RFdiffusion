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

  // Classifier-Free Guidance (CFG)
  // Use when applying H-bond or RASA conditioning to improve adherence
  use_classifier_free_guidance?: boolean;
  cfg_scale?: number;  // Recommended: 2.0 when CFG enabled

  // Motif specifications
  select_hotspots?: Record<string, string>;
  select_fixed_atoms?: Record<string, string>;
  select_unfixed_sequence?: Record<string, string>;

  // Ligand binding
  ligand?: string;
  select_buried?: Record<string, string>;
  select_exposed?: Record<string, string>;
  select_partially_buried?: Record<string, string>;

  // H-bond conditioning
  select_hbond_donor?: Record<string, string>;
  select_hbond_acceptor?: Record<string, string>;

  // Origin/center-of-mass positioning
  ori_token?: [number, number, number];  // Explicit XYZ coordinates

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

  // Covalent modifications (for enzyme design with covalent inhibitors/substrates)
  covalent_bonds?: Array<{
    protein: { chain: string; res_name: string; res_num: number; atom_name: string };
    ligand: { chain: string; res_name: string; res_num: number; atom_name: string };
  }>;

  // Interface ligand design (separable dimer)
  task?: string;  // Task type for specialized handlers (e.g., 'interface_ligand_design', 'interface_metal_design', 'protein_binder_design')
  approach?:
    // Interface ligand approaches
    | 'asymmetric' | 'sequential' | 'full' | 'joint' | 'asymmetric_rasa' | 'induced' | 'symmetric'
    // Interface metal approaches
    | 'joint_metal' | 'asymmetric_metal' | 'induced_metal' | 'bridging_metal' | 'redox_switch' | 'lanthanide_sensor';  // Design approach
  ligand_smiles?: string;  // SMILES string for ligand
  chain_length?: string;  // Chain length range (e.g., "60-80")
  side?: 'left' | 'right';  // Which side of ligand to bind
  ori_offset?: number[];  // Offset from ligand center [x, y, z]
  use_ori_token?: boolean;  // Whether to use ori_token positioning

  // Interface metal ligand design
  template_name?: string;  // Template name for predefined metal complexes
  complex_pdb?: string;  // Custom PDB content for metal complex
  ligand_name?: string;  // Custom ligand name
  validate_coordination?: boolean;  // Validate coordination geometry
  bias_AA?: string;  // HSAB bias amino acid string

  // Interface metal design (metal-coordinated dimer)
  metal?: string;  // Metal ion code (ZN, FE, CA, TB, etc.)
  coordination_split?: [number, number];  // [chain_a_donors, chain_b_donors]
  coordination_geometry?: string;  // tetrahedral, octahedral, square_planar, etc.
  chain_a_donors?: string[];  // Donor types for chain A (His, Cys, Asp, etc.)
  chain_b_donors?: string[];  // Donor types for chain B
  allow_waters?: boolean;  // Allow water molecules in coordination sphere
  bridging_metal?: string;  // Second metal for bridging designs
  enable_luminescence?: boolean;  // Enable luminescence for lanthanide designs

  // Extended metal design fields
  metal_config?: {
    coordination_split: [number, number];
    geometry: string;
    chain_a_donors: string[];
    chain_b_donors: string[];
    include_waters: boolean;
    second_metal?: string;
  };
  parametric?: {
    coordination_number: number;
    num_waters: number;
    bidentate_fraction: number;
  };
  template_type?: string;  // Template type for legacy mode
  donor_residue?: string;  // Donor residue for legacy mode
  use_motif_scaffolding?: boolean;  // Use motif scaffolding
  add_trp_antenna?: boolean;  // Add tryptophan antenna for lanthanides

  // Protein binder design specific fields
  target_pdb?: string;  // PDB content for target protein
  binder_length?: string;  // Binder length range (e.g., "60-80")
  quality_threshold?: 'relaxed' | 'standard' | 'strict';  // Quality filter threshold
  protocol?: string;  // Protocol preset (e.g., 'miniprotein_default', 'peptide_default')
  validate_structure?: boolean;  // Enable ESMFold structure validation
  auto_hotspots?: boolean;  // Auto-detect hotspots using SASA analysis
  filter_wrap_around?: boolean;  // Filter out wrap-around binders using Rg ratio
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

// ====== PLIP Interaction Analysis Types ======

// Filter mode for success criteria display
export type FilterMode = 'heterodimer' | 'rfdiffusion3';

// PLIP interaction profile from ligand analysis
export interface InteractionProfile {
  hydrogen_bonds: number;
  hydrophobic_contacts: number;
  pi_stacking: number;
  salt_bridges: number;
  halogen_bonds?: number;
  total: number;
}

// Heterodimer-specific metrics (custom workflow)
export interface HeterodimerMetrics {
  affinity?: number;           // GNINA: < -5 = good (kcal/mol)
  contacts_a?: number;         // ≥ 5 = good
  contacts_b?: number;         // ≥ 5 = good
  sequence_identity?: number;  // < 70% = true heterodimer
  anti_homo_score?: number;    // > 60 = good (0-100 scale)
  n7_hbonds?: number;          // ≥ 1 = azobenzene N7 (azo) satisfied
  n8_hbonds?: number;          // ≥ 1 = azobenzene N8 (azo) satisfied
  has_clashes?: boolean;       // false = good
  is_heterodimer?: boolean;    // true = success
}

// Standard RFdiffusion3 metrics (Baker Lab reference)
export interface RFD3Metrics {
  pae_interaction?: number;  // < 10 = good
  plddt?: number;            // > 80 = good, > 90 = excellent
  ddg?: number;              // < -40 = good (Rosetta kcal/mol)
  total_hbonds?: number;     // > 11 = good
  shape_complementarity?: number; // > 0.65 = excellent
}

// Combined design metrics
export interface DesignMetrics extends HeterodimerMetrics, RFD3Metrics {}

// Design result with all possible fields
export interface DesignResult {
  pdb_content: string;
  cif_content?: string;
  sequence?: string;
  metrics?: DesignMetrics;
  interactions?: InteractionProfile;
  key_binding_residues?: string[];
  recommendations?: string[];
  analysis_method?: 'plip' | 'distance_based';
  interaction_summary?: string;
}

// Pipeline stage types
export type PipelineStage = 'idle' | 'backbone' | 'ligandmpnn' | 'validation' | 'analysis' | 'complete' | 'error';

// Stage details for pipeline progress
export interface StageDetails {
  backbone?: { progress?: number; message?: string };
  ligandmpnn?: { sequences?: number; message?: string };
  validation?: { message?: string };
  analysis?: { message?: string };
}
