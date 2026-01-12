/**
 * Demo mode case study data for the AI Design Assistant
 * Provides simulated data when backend is not connected
 */

import type { MetalBindingAnalysis, DesignEvaluation } from './api';

// Types for ligand design
export interface LigandAnalysis {
  success: boolean;
  ligand: {
    name: string;
    smiles: string;
    num_atoms: number;
    num_heavy_atoms: number;
    molecular_weight: number;
    num_rings: number;
    num_rotatable_bonds: number;
    num_hbd: number;
    num_hba: number;
    logp: number;
    symmetry: 'symmetric' | 'asymmetric';
    functional_groups: string[];
  };
  binding_site: {
    type: 'interface' | 'buried' | 'surface';
    approach: string;
    description: string;
    recommended_chain_length: string;
  };
  topology: {
    separable: boolean;
    description: string;
    ring_threading_risk: 'low' | 'medium' | 'high';
  };
  atom_analysis: {
    face_a_atoms: string[];
    face_b_atoms: string[];
    core_atoms: string[];
    contact_potential: number;
  };
  recommendations: string[];
}

export interface LigandEvaluation {
  success: boolean;
  approach: 'full' | 'asymmetric' | 'sequential';
  chain_a?: {
    contacts: number;
    exposed_atoms: string;
    affinity: number;
    hbonds?: number;
    hydrophobic_contacts?: number;
  };
  chain_b?: {
    contacts: number;
    exposed_atoms: string;
    affinity: number;
    hbonds?: number;
    hydrophobic_contacts?: number;
  };
  dimer?: {
    affinity: number;
    contacts_a: number;
    contacts_b: number;
    has_clashes: boolean;
    separable: boolean;
    interface_area?: number;
    energy_breakdown?: {
      vdw: number;
      electrostatic: number;
      hbond: number;
      desolvation: number;
    };
    shape_complementarity?: number;
  };
  gnina_score?: {
    cnn_score: number;
    cnn_affinity: number;
  };
  overall_pass: boolean;
}

export interface DesignResult {
  id: string;
  rank: number;
  pdbContent: string;
  metrics: {
    affinity?: number;
    contacts_a?: number;
    contacts_b?: number;
    has_clashes: boolean;
    separable?: boolean;
    interface_area?: number;
  };
  chain_a_metrics?: {
    contacts: number;
    exposed_atoms: string[];
    affinity: number;
  };
  chain_b_metrics?: {
    contacts: number;
    exposed_atoms: string[];
    affinity: number;
  };
}

export interface BinderDesign {
  rank: number;
  binder_sequence: string;
  mpnn_score: number;
  esm_perplexity: number;
  esm_confidence: number;
  interface_contacts: number;
  interface_hbonds: number;
  buried_sasa: number;
  packstat: number;
  pdb_content?: string;
  shape_complementarity?: number;
  surface_hydrophobicity?: number;
  unsaturated_hbonds?: number;
  interface_residue_count?: number;
  esmfold_plddt?: number;
  esmfold_rmsd?: number;
  esmfold_validation_passed?: boolean;
}

export interface BinderStatistics {
  generated: number;
  mpnn_designed: number;
  esm_passed: number;
  relaxed: number;
  interface_analyzed: number;
  passed_filters: number;
  returned: number;
}

export interface BinderEvaluation {
  success: boolean;
  statistics: BinderStatistics;
  best_design?: {
    rank: number;
    binder_sequence: string;
    mpnn_score: number;
    esm_perplexity: number;
    esm_confidence: number;
    interface_contacts: number;
    interface_hbonds: number;
    buried_sasa: number;
    packstat: number;
  };
  overall_pass: boolean;
}

// Demo PDB content for binder visualization
const DEMO_BINDER_PDB = `HEADER    DEMO BINDER DESIGN
TITLE     SIMULATED PROTEIN BINDER COMPLEX
REMARK    This is a simplified demo structure for visualization
REMARK    Chain A = Target (Rubredoxin fragment)
REMARK    Chain B = Designed Binder
ATOM      1  N   MET A   1       0.000   0.000   0.000  1.00 50.00           N
ATOM      2  CA  MET A   1       1.458   0.000   0.000  1.00 50.00           C
ATOM      3  C   MET A   1       2.009   1.420   0.000  1.00 50.00           C
ATOM      4  O   MET A   1       1.246   2.390   0.000  1.00 50.00           O
ATOM      5  N   LYS A   2       3.300   1.554   0.000  1.00 50.00           N
ATOM      6  CA  LYS A   2       3.933   2.868   0.000  1.00 50.00           C
ATOM      7  C   LYS A   2       5.451   2.764   0.000  1.00 50.00           C
ATOM      8  O   LYS A   2       6.034   1.684   0.000  1.00 50.00           O
ATOM      9  N   CYS A   3       6.089   3.928   0.000  1.00 50.00           N
ATOM     10  CA  CYS A   3       7.536   4.022   0.000  1.00 50.00           C
ATOM     11  C   CYS A   3       8.169   2.644   0.000  1.00 50.00           C
ATOM     12  O   CYS A   3       7.451   1.636   0.000  1.00 50.00           O
ATOM     13  SG  CYS A   3       8.060   5.000   1.500  1.00 50.00           S
ATOM     14  N   GLU A   4       9.480   2.591   0.000  1.00 50.00           N
ATOM     15  CA  GLU A   4      10.224   1.341   0.000  1.00 50.00           C
ATOM     16  C   GLU A   4      11.729   1.541   0.000  1.00 50.00           C
ATOM     17  O   GLU A   4      12.203   2.680   0.000  1.00 50.00           O
TER
ATOM    101  N   MET B   1      15.000  10.000   5.000  1.00 60.00           N
ATOM    102  CA  MET B   1      16.458  10.000   5.000  1.00 60.00           C
ATOM    103  C   MET B   1      17.009  11.420   5.000  1.00 60.00           C
ATOM    104  O   MET B   1      16.246  12.390   5.000  1.00 60.00           O
ATOM    105  N   LYS B   2      18.300  11.554   5.000  1.00 60.00           N
ATOM    106  CA  LYS B   2      18.933  12.868   5.000  1.00 60.00           C
ATOM    107  C   LYS B   2      20.451  12.764   5.000  1.00 60.00           C
ATOM    108  O   LYS B   2      21.034  11.684   5.000  1.00 60.00           O
ATOM    109  N   PHE B   3      21.089  13.928   5.000  1.00 60.00           N
ATOM    110  CA  PHE B   3      22.536  14.022   5.000  1.00 60.00           C
ATOM    111  C   PHE B   3      23.169  12.644   5.000  1.00 60.00           C
ATOM    112  O   PHE B   3      22.451  11.636   5.000  1.00 60.00           O
ATOM    113  N   LEU B   4      24.480  12.591   5.000  1.00 60.00           N
ATOM    114  CA  LEU B   4      25.224  11.341   5.000  1.00 60.00           C
ATOM    115  C   LEU B   4      26.729  11.541   5.000  1.00 60.00           C
ATOM    116  O   LEU B   4      27.203  12.680   5.000  1.00 60.00           O
ATOM    117  N   ILE B   5      27.400  10.400   5.000  1.00 60.00           N
ATOM    118  CA  ILE B   5      28.850  10.300   5.000  1.00 60.00           C
ATOM    119  C   ILE B   5      29.400   8.880   5.000  1.00 60.00           C
ATOM    120  O   ILE B   5      28.650   7.900   5.000  1.00 60.00           O
TER
END
`;

// Metal binding case study (Rubredoxin)
export const METAL_CASE_STUDY = {
  pdbInfo: {
    pdb_id: '1BRF',
    title: 'RUBREDOXIN FROM CLOSTRIDIUM PASTEURIANUM',
    chains: ['A'],
    num_atoms: 424,
    num_residues: 54,
    metals: [{ chain: 'A', residue: 'FE', res_num: '55', element: 'FE' }],
  },
  metalAnalysis: {
    success: true,
    metal: {
      element: 'FE',
      chain: 'A',
      residue: 'FE',
      resnum: 55,
      position: [0, 0, 0],
    },
    coordination: {
      number: 4,
      geometry: 'tetrahedral',
      geometry_rmsd: 0.15,
      coordinating_atoms: [
        { chain: 'A', residue: 'CYS', residue_number: 5, atom: 'SG', element: 'S', distance: 2.28, donor_type: 'S_thiolate' },
        { chain: 'A', residue: 'CYS', residue_number: 8, atom: 'SG', element: 'S', distance: 2.31, donor_type: 'S_thiolate' },
        { chain: 'A', residue: 'CYS', residue_number: 38, atom: 'SG', element: 'S', distance: 2.29, donor_type: 'S_thiolate' },
        { chain: 'A', residue: 'CYS', residue_number: 41, atom: 'SG', element: 'S', distance: 2.30, donor_type: 'S_thiolate' },
      ],
    },
    donor_analysis: {
      types: { 'S_thiolate': 4 },
      dominant_type: 'S_thiolate',
      lanthanide_compatible: false,
    },
    bond_analysis: {
      average_distance: 2.30,
      min_distance: 2.28,
      max_distance: 2.31,
      distances: [2.28, 2.31, 2.29, 2.30],
    },
    suggestions: [
      'Lanthanides prefer O-donors over S-donors. Consider replacing Cys with Asp/Glu.',
    ],
  } as MetalBindingAnalysis,
  simulatedEvaluation: {
    success: true,
    coordination_number: 8,
    avg_bond_distance: 2.42,
    geometry_type: 'square_antiprism',
    geometry_rmsd: 0.28,
    oxygen_donors: 7,
    nitrogen_donors: 1,
    sulfur_donors: 0,
    criteria_passed: 4,
    criteria_total: 4,
    overall_pass: true,
    suggestions: [],
  } as DesignEvaluation,
};

// Azobenzene interface dimer case study
export const AZOBENZENE_CASE_STUDY = {
  pdbInfo: {
    pdb_id: 'AZOB',
    title: 'Azobenzene Interface Dimer Design',
    description: 'Design separable protein dimers with azobenzene at the interface',
    ligand: {
      name: 'Azobenzene',
      smiles: 'c1ccc(cc1)N=Nc2ccccc2',
      description: 'Light-switchable diazo compound with two phenyl rings',
    },
  },
  ligandAnalysis: {
    success: true,
    ligand: {
      name: 'Azobenzene',
      smiles: 'c1ccc(cc1)N=Nc2ccccc2',
      num_atoms: 24,
      num_heavy_atoms: 14,
      molecular_weight: 182.22,
      num_rings: 2,
      num_rotatable_bonds: 1,
      num_hbd: 0,
      num_hba: 2,
      logp: 3.4,
      symmetry: 'symmetric' as const,
      functional_groups: ['Azo', 'Phenyl', 'Phenyl'],
    },
    binding_site: {
      type: 'interface' as const,
      approach: 'asymmetric_sequential',
      description: 'Ligand sits at dimer interface - chains A and B contact opposite sides',
      recommended_chain_length: '60-80 residues',
    },
    topology: {
      separable: true,
      description: 'Chains can physically separate (unlike buried ligand designs)',
      ring_threading_risk: 'low' as const,
    },
    atom_analysis: {
      face_a_atoms: ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'N1'],
      face_b_atoms: ['C8', 'C9', 'C10', 'C11', 'C12', 'C13', 'N2'],
      core_atoms: ['N1', 'N2'],
      contact_potential: 14.5,
    },
    recommendations: [
      'Symmetric ligand ideal for balanced dimer contacts',
      'Low ring threading risk - safe for interface design',
      'Consider hydrophobic residues (Phe, Leu) for phenyl contacts',
    ],
  } as LigandAnalysis,
  simulatedDesigns: [
    {
      id: 'design-1',
      rank: 1,
      pdbContent: 'SIMULATED_PDB_DESIGN_1',
      metrics: {
        affinity: -4.2,
        contacts_a: 8,
        contacts_b: 6,
        has_clashes: false,
        separable: true,
        interface_area: 850,
      },
      chain_a_metrics: {
        contacts: 8,
        exposed_atoms: ['C1', 'C2', 'C3', 'C4', 'C5', 'C6'],
        affinity: -2.1,
      },
      chain_b_metrics: {
        contacts: 6,
        exposed_atoms: ['C8', 'C9', 'C10', 'C11', 'C12', 'C13'],
        affinity: -1.8,
      },
    },
    {
      id: 'design-2',
      rank: 2,
      pdbContent: 'SIMULATED_PDB_DESIGN_2',
      metrics: {
        affinity: -3.8,
        contacts_a: 7,
        contacts_b: 7,
        has_clashes: false,
        separable: true,
        interface_area: 920,
      },
      chain_a_metrics: {
        contacts: 7,
        exposed_atoms: ['C1', 'C2', 'C3', 'C5', 'C6'],
        affinity: -1.9,
      },
      chain_b_metrics: {
        contacts: 7,
        exposed_atoms: ['C8', 'C9', 'C10', 'C12', 'C13'],
        affinity: -1.9,
      },
    },
    {
      id: 'design-3',
      rank: 3,
      pdbContent: 'SIMULATED_PDB_DESIGN_3',
      metrics: {
        affinity: -3.5,
        contacts_a: 6,
        contacts_b: 5,
        has_clashes: false,
        separable: true,
        interface_area: 780,
      },
      chain_a_metrics: {
        contacts: 6,
        exposed_atoms: ['C1', 'C2', 'C4', 'C5', 'C6'],
        affinity: -1.8,
      },
      chain_b_metrics: {
        contacts: 5,
        exposed_atoms: ['C9', 'C10', 'C11', 'C12'],
        affinity: -1.5,
      },
    },
    {
      id: 'design-4',
      rank: 4,
      pdbContent: 'SIMULATED_PDB_DESIGN_4',
      metrics: {
        affinity: -2.9,
        contacts_a: 5,
        contacts_b: 4,
        has_clashes: true,
        separable: true,
        interface_area: 650,
      },
      chain_a_metrics: {
        contacts: 5,
        exposed_atoms: ['C1', 'C3', 'C4', 'C5'],
        affinity: -1.5,
      },
      chain_b_metrics: {
        contacts: 4,
        exposed_atoms: ['C9', 'C10', 'C12'],
        affinity: -1.2,
      },
    },
    {
      id: 'design-5',
      rank: 5,
      pdbContent: 'SIMULATED_PDB_DESIGN_5',
      metrics: {
        affinity: -2.5,
        contacts_a: 4,
        contacts_b: 5,
        has_clashes: false,
        separable: false,
        interface_area: 580,
      },
      chain_a_metrics: {
        contacts: 4,
        exposed_atoms: ['C2', 'C3', 'C5'],
        affinity: -1.2,
      },
      chain_b_metrics: {
        contacts: 5,
        exposed_atoms: ['C8', 'C10', 'C11', 'C12'],
        affinity: -1.4,
      },
    },
  ] as DesignResult[],
  simulatedEvaluation: {
    success: true,
    approach: 'full' as const,
    chain_a: {
      contacts: 8,
      exposed_atoms: 'C1,C2,C3,C4,C5,C6',
      affinity: -2.1,
      hbonds: 2,
      hydrophobic_contacts: 6,
    },
    chain_b: {
      contacts: 6,
      exposed_atoms: 'C8,C9,C10,C11,C12,C13',
      affinity: -1.8,
      hbonds: 1,
      hydrophobic_contacts: 5,
    },
    dimer: {
      affinity: -4.2,
      contacts_a: 8,
      contacts_b: 6,
      has_clashes: false,
      separable: true,
      interface_area: 850,
      energy_breakdown: {
        vdw: -2.8,
        electrostatic: -0.5,
        hbond: -0.6,
        desolvation: -0.3,
      },
      shape_complementarity: 0.72,
    },
    gnina_score: {
      cnn_score: 0.85,
      cnn_affinity: -5.2,
    },
    overall_pass: true,
  } as LigandEvaluation,
};

// Protein binder design case study
export const BINDER_CASE_STUDY = {
  pdbInfo: {
    pdb_id: 'BIND',
    title: 'Protein Binder Design Pipeline',
    description: 'Design high-affinity protein binders with multi-stage validation',
    target: {
      name: 'Rubredoxin',
      pdb_id: '1BRF',
      chains: ['A'],
      residues: 54,
    },
  },
  simulatedStatistics: {
    generated: 4,
    mpnn_designed: 4,
    esm_passed: 4,
    relaxed: 4,
    interface_analyzed: 4,
    passed_filters: 4,
    returned: 2,
  } as BinderStatistics,
  simulatedDesigns: [
    {
      rank: 1,
      binder_sequence: 'MKFLILLFNILCSGACMDELEEEELKKLEEEEKQKLEEEEKQKLEEEKK',
      mpnn_score: 0.85,
      esm_perplexity: 4.2,
      esm_confidence: 0.89,
      interface_contacts: 6842,
      interface_hbonds: 14,
      buried_sasa: 1250.5,
      packstat: 0.72,
      pdb_content: DEMO_BINDER_PDB,
    },
    {
      rank: 2,
      binder_sequence: 'MKVFLILNILCSGACMDELEEEELKKLEEEKQKLEEEEEKQKLEEEEKK',
      mpnn_score: 0.82,
      esm_perplexity: 5.1,
      esm_confidence: 0.84,
      interface_contacts: 5921,
      interface_hbonds: 11,
      buried_sasa: 1180.2,
      packstat: 0.68,
      pdb_content: DEMO_BINDER_PDB,
    },
  ] as BinderDesign[],
  simulatedEvaluation: {
    success: true,
    statistics: {
      generated: 4,
      mpnn_designed: 4,
      esm_passed: 4,
      relaxed: 4,
      interface_analyzed: 4,
      passed_filters: 4,
      returned: 2,
    },
    best_design: {
      rank: 1,
      binder_sequence: 'MKFLILLFNILCSGACMDELEEEELKKLEEEEKQKLEEEEKQKLEEEKK',
      mpnn_score: 0.85,
      esm_perplexity: 4.2,
      esm_confidence: 0.89,
      interface_contacts: 6842,
      interface_hbonds: 14,
      buried_sasa: 1250.5,
      packstat: 0.72,
    },
    overall_pass: true,
  } as BinderEvaluation,
};

/**
 * Get demo data for a specific workflow type
 */
export function getDemoDataForWorkflow(
  workflowType: 'metal' | 'ligand' | 'binder'
): {
  statistics?: BinderStatistics;
  designs?: BinderDesign[] | DesignResult[];
  evaluation?: DesignEvaluation | LigandEvaluation | BinderEvaluation;
  analysis?: MetalBindingAnalysis | LigandAnalysis;
} {
  switch (workflowType) {
    case 'metal':
      return {
        evaluation: METAL_CASE_STUDY.simulatedEvaluation,
        analysis: METAL_CASE_STUDY.metalAnalysis,
      };
    case 'ligand':
      return {
        designs: AZOBENZENE_CASE_STUDY.simulatedDesigns,
        evaluation: AZOBENZENE_CASE_STUDY.simulatedEvaluation,
        analysis: AZOBENZENE_CASE_STUDY.ligandAnalysis,
      };
    case 'binder':
      return {
        statistics: BINDER_CASE_STUDY.simulatedStatistics,
        designs: BINDER_CASE_STUDY.simulatedDesigns,
        evaluation: BINDER_CASE_STUDY.simulatedEvaluation,
      };
    default:
      return {};
  }
}

/**
 * Check if a design result is from demo mode (simulated content)
 */
export function isSimulatedContent(pdbContent: string | undefined): boolean {
  if (!pdbContent) return true;
  return pdbContent.startsWith('SIMULATED_');
}

/**
 * Count passing designs in a ligand design result set
 */
export function countPassingLigandDesigns(designs: DesignResult[]): number {
  return designs.filter(
    (d) => !d.metrics.has_clashes && d.metrics.separable !== false
  ).length;
}
