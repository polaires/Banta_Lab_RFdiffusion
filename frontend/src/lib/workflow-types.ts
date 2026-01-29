/**
 * Workflow Types — TypeScript interfaces for the modular workflow system.
 *
 * Matches the backend JSON workflow spec format so the AI assistant
 * can compose workflows that the frontend can display and track.
 */

/** A single step in a workflow spec */
export interface WorkflowStepSpec {
  module: string;
  params?: Record<string, unknown>;
}

/** Strategy configuration for iterative workflows */
export interface WorkflowStrategy {
  type: 'scout' | 'sweep';
  ptm_threshold?: number;
  plddt_threshold?: number;
  grid?: Record<string, unknown[]>;
  designs_per_config?: number;
}

/** Complete workflow specification (matches backend JSON spec) */
export interface WorkflowSpec {
  name: string;
  backend?: 'auto' | 'in_process' | 'http';
  params: Record<string, unknown>;
  steps: WorkflowStepSpec[];
  strategy?: WorkflowStrategy;
  checkpoint_dir?: string;
}

/** Per-step result in a running workflow */
export interface WorkflowStepResult {
  status: 'completed' | 'failed';
  summary: string;
  timing?: number;
  metrics?: Record<string, number>;
}

/** Progress shape for workflow_run jobs (returned by polling) */
export interface WorkflowProgress {
  status: 'running' | 'completed' | 'failed';
  workflow_name: string;
  current_step: string;
  current_step_index: number;
  total_steps: number;
  step_results: Record<string, WorkflowStepResult>;
  result?: WorkflowResult;
}

/** Final result from a completed workflow */
export interface WorkflowResult {
  params: Record<string, unknown>;
  num_backbones: number;
  num_sequences: number;
  num_predictions: number;
  backbone_pdbs?: string[];
  sequences?: Array<{
    sequence: string;
    score: number;
    confidence: number;
    backbone_index: number;
  }>;
  predictions?: Array<{
    plddt: number;
    ptm: number;
    rmsd: number | null;
    passed: boolean;
    predictor: string;
    backbone_index: number;
  }>;
  analysis?: Record<string, unknown>;
  report?: string;
}

/** Map module names to human-readable display names */
export const MODULE_DISPLAY_NAMES: Record<string, string> = {
  intent_parser: 'Parse Intent',
  design_configurator: 'Configure Design',
  scaffolder: 'Extract Scaffold',
  backbone_generator: 'Generate Backbones',
  sequence_designer: 'Design Sequences',
  structure_predictor: 'Predict Structures',
  analyzer: 'Analyze Results',
  reporter: 'Generate Report',
  scout_filter: 'Scout Filter',
  sweep_configurator: 'Sweep Config',
};

/** Map module names to short descriptions */
export const MODULE_DESCRIPTIONS: Record<string, string> = {
  intent_parser: 'Parse natural language into design parameters',
  design_configurator: 'Generate RFD3 and MPNN configurations',
  scaffolder: 'Extract motifs from PDB structures',
  backbone_generator: 'Generate protein backbone structures via RFD3',
  sequence_designer: 'Design amino acid sequences via LigandMPNN',
  structure_predictor: 'Predict 3D structures via RF3/ESMFold',
  analyzer: 'Analyze and filter design candidates',
  reporter: 'Generate summary report with recommendations',
};
