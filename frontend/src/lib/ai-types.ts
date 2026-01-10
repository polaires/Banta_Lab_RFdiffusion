/**
 * AI-Driven Protein Engineering Types
 *
 * TypeScript types matching the backend AI engine for:
 * - Task detection and routing
 * - Dynamic form configuration
 * - Parameter generation and refinement
 */

// Supported design tasks (RFdiffusion + ESM3)
export type TaskType =
  | 'de_novo'
  | 'binder'
  | 'metal_redesign'
  | 'enzyme'
  | 'refinement'
  | 'symmetric'
  | 'scaffold'
  // ESM3 tasks
  | 'esm3_score'     // Score sequences for quality
  | 'esm3_generate'  // Generate sequences with function conditioning
  | 'esm3_embed'     // Get sequence embeddings
  | 'unknown';

// Form field types for dynamic UI generation
export type FormFieldType =
  | 'text'
  | 'number'
  | 'select'
  | 'range'
  | 'structure'
  | 'residue_picker';

// Option for select fields
export interface SelectOption {
  value: string;
  label: string;
}

// Range configuration for slider fields
export interface RangeConfig {
  min: number;
  max: number;
  step: number;
}

// Configuration for a dynamic form field
export interface FormFieldConfig {
  id: string;
  label: string;
  type: FormFieldType;
  required: boolean;
  default_value?: unknown;
  suggested_value?: unknown;
  options?: SelectOption[];
  range_config?: RangeConfig;
  help_text: string;
  ai_reasoning?: string;
}

// Response from AI analysis endpoint
export interface AIAnalyzeResponse {
  success: boolean;
  task_type: TaskType;
  params: Record<string, unknown>;
  reasoning: string;
  confidence: number;
  clarifying_questions: string[];
  form_config: FormFieldConfig[];
  error?: string;
}

// Request to AI analysis endpoint
export interface AIAnalyzeRequest {
  user_input: string;
  pdb_content?: string;
  pdb_id?: string;
  structure_info?: {
    chains: string[];
    num_residues: number;
    num_atoms?: number;
  };
  conversation_history?: Array<{
    role: 'user' | 'assistant';
    content: string;
  }>;
}

// Request to refine parameters
export interface AIRefineRequest {
  current_params: Record<string, unknown>;
  user_feedback: string;
  task_type: TaskType;
}

// Response from refine endpoint
export interface AIRefineResponse {
  success: boolean;
  params: Record<string, unknown>;
  reasoning: string;
  confidence: number;
}

// Chat message in AI conversation
export interface AIChatMessage {
  id: string;
  role: 'user' | 'assistant' | 'system';
  content: string;
  timestamp: number;
  data?: {
    type: 'analysis' | 'params' | 'question' | 'result';
    payload: AIAnalyzeResponse | Record<string, unknown>;
  };
}

// Workflow state for AI-driven design
export interface AIWorkflowState {
  phase: AIWorkflowPhase;
  taskType: TaskType | null;
  suggestedParams: Record<string, unknown> | null;
  formConfig: FormFieldConfig[] | null;
  clarifyingQuestions: string[];
  messages: AIChatMessage[];
  confidence: number;
  isProcessing: boolean;
  error: string | null;
}

// Workflow phases
export type AIWorkflowPhase =
  | 'initial'           // Waiting for user input
  | 'input_structure'   // Waiting for structure upload/PDB code
  | 'analyzing'         // AI analyzing request
  | 'showing_form'      // Smart form displayed with suggestions
  | 'refining'          // User chatting to refine params
  | 'confirming'        // Final confirmation before run
  | 'running'           // Job executing
  | 'evaluating'        // AI analyzing results
  | 'complete'          // Done, showing results
  | 'iterating';        // User requested adjustments

// Task-specific presets
export interface TaskPreset {
  id: TaskType;
  label: string;
  description: string;
  icon: string;
  requiresStructure: boolean;
  defaultParams: Record<string, unknown>;
}

// All available task presets
export const TASK_PRESETS: TaskPreset[] = [
  {
    id: 'de_novo',
    label: 'De Novo Design',
    description: 'Create a new protein from scratch',
    icon: 'sparkles',
    requiresStructure: false,
    defaultParams: { length: '80-100', num_designs: 10 },
  },
  {
    id: 'metal_redesign',
    label: 'Metal Binding',
    description: 'Redesign or create metal binding sites',
    icon: 'atom',
    requiresStructure: true,
    defaultParams: { partial_t: 15, ligand: 'TB' },
  },
  {
    id: 'binder',
    label: 'Binder Design',
    description: 'Design a protein to bind a target',
    icon: 'link',
    requiresStructure: true,
    defaultParams: { step_scale: 3.0, gamma_0: 0.2 },
  },
  {
    id: 'refinement',
    label: 'Structure Refinement',
    description: 'Improve an existing protein structure',
    icon: 'tool',
    requiresStructure: true,
    defaultParams: { partial_t: 6, gamma_0: 0.8 },
  },
  {
    id: 'symmetric',
    label: 'Symmetric Oligomer',
    description: 'Design homo-oligomeric complexes',
    icon: 'shape',
    requiresStructure: false,
    defaultParams: { symmetry: { id: 'C3' }, length: '60-80' },
  },
  {
    id: 'enzyme',
    label: 'Enzyme Design',
    description: 'Design or modify catalytic sites',
    icon: 'flask',
    requiresStructure: true,
    defaultParams: { partial_t: 12 },
  },
  // ESM3 tasks
  {
    id: 'esm3_score',
    label: 'Score Sequences',
    description: 'Evaluate sequence quality using ESM3',
    icon: 'score',
    requiresStructure: false,
    defaultParams: {},
  },
  {
    id: 'esm3_generate',
    label: 'Generate Sequence',
    description: 'Generate sequences with function conditioning',
    icon: 'auto_awesome',
    requiresStructure: false,
    defaultParams: { num_sequences: 4, temperature: 0.7 },
  },
];

// Helper to get task preset by ID
export function getTaskPreset(taskType: TaskType): TaskPreset | undefined {
  return TASK_PRESETS.find((p) => p.id === taskType);
}

// Helper to format confidence as percentage
export function formatConfidence(confidence: number): string {
  return `${Math.round(confidence * 100)}%`;
}

// Helper to get confidence color
export function getConfidenceColor(confidence: number): string {
  if (confidence >= 0.8) return 'text-green-600';
  if (confidence >= 0.6) return 'text-yellow-600';
  return 'text-red-600';
}

// Metal options for selection
export const METAL_OPTIONS: SelectOption[] = [
  { value: 'TB', label: 'Terbium (TB) - Luminescent' },
  { value: 'GD', label: 'Gadolinium (GD) - MRI contrast' },
  { value: 'EU', label: 'Europium (EU) - Red emission' },
  { value: 'CA', label: 'Calcium (CA) - Signaling' },
  { value: 'ZN', label: 'Zinc (ZN) - Catalytic' },
  { value: 'FE', label: 'Iron (FE) - Redox' },
  { value: 'CU', label: 'Copper (CU) - Electron transfer' },
  { value: 'MN', label: 'Manganese (MN) - Catalysis' },
  { value: 'MG', label: 'Magnesium (MG) - ATP binding' },
];

// Symmetry options for selection
export const SYMMETRY_OPTIONS: SelectOption[] = [
  { value: 'C2', label: 'C2 - Dimer (2 subunits)' },
  { value: 'C3', label: 'C3 - Trimer (3 subunits)' },
  { value: 'C4', label: 'C4 - Tetramer (4 subunits)' },
  { value: 'C5', label: 'C5 - Pentamer (5 subunits)' },
  { value: 'C6', label: 'C6 - Hexamer (6 subunits)' },
  { value: 'D2', label: 'D2 - Dihedral (4 subunits)' },
  { value: 'D3', label: 'D3 - Dihedral (6 subunits)' },
];
