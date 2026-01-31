import type { LucideIcon } from 'lucide-react';

// ---- Step Result ----

/** Output item from a step that produced PDB structures */
export interface PdbOutput {
  id: string;
  label: string;
  pdbContent: string;
  /** Confidence or quality metrics */
  metrics?: Record<string, number | string>;
  /** Extracted sequence */
  sequence?: string;
}

/** Output item from a step that produced sequences */
export interface SequenceOutput {
  id: string;
  sequence: string;
  score?: number;
  label?: string;
  metrics?: Record<string, number | string>;
}

/** Generic output from any pipeline step */
export interface StepResult {
  /** Unique ID for this result set */
  id: string;
  /** Human-readable summary */
  summary: string;
  /** PDB structures produced */
  pdbOutputs?: PdbOutput[];
  /** Sequences produced */
  sequences?: SequenceOutput[];
  /** Arbitrary analysis data for custom result previews */
  data?: Record<string, unknown>;
}

// ---- Step Status ----

export type StepStatus =
  | 'pending'     // Not yet started
  | 'ready'       // Dependencies met, can run
  | 'running'     // Currently executing
  | 'paused'      // Completed, waiting for user review
  | 'completed'   // Done, user approved or auto-continued
  | 'failed'      // Error occurred
  | 'skipped';    // User chose to skip

// ---- Step Execution ----

/** Context passed into every step's execute function */
export interface StepExecutionContext {
  /** Results from ALL previous steps, keyed by step ID */
  previousResults: Record<string, StepResult>;
  /** IDs of outputs the user selected to carry forward from the prior step */
  selectedItems: string[];
  /** User-edited parameters for THIS step */
  params: Record<string, unknown>;
  /** Initial pipeline parameters */
  initialParams: Record<string, unknown>;
  /** Signal for cancellation */
  abortSignal: AbortSignal;
  /** Callback for granular progress within a step (0-100) */
  onProgress: (percent: number, message?: string) => void;
  /** Callback when a backend job is created */
  onJobCreated?: (jobId: string, type: 'rfd3' | 'rf3' | 'mpnn' | 'workflow') => void;
}

// ---- Parameter Schema ----

/** Descriptor for a single editable parameter */
export interface StepParameterField {
  id: string;
  label: string;
  type: 'number' | 'text' | 'select' | 'slider' | 'boolean' | 'hidden';
  required: boolean;
  defaultValue?: unknown;
  options?: Array<{ value: string; label: string }>;
  range?: { min: number; max: number; step: number };
  helpText?: string;
}

// ---- Step Definition ----

/** Static definition of a single pipeline step */
export interface PipelineStepDefinition {
  /** Unique ID within the pipeline (e.g., 'rfd3_backbone') */
  id: string;
  /** Display name */
  name: string;
  /** Short description */
  description: string;
  /** Lucide icon */
  icon: LucideIcon;
  /** Whether the user MUST review results before continuing */
  requiresReview: boolean;
  /** Whether the user can select a subset of outputs to carry forward */
  supportsSelection: boolean;
  /** Whether this step can be skipped */
  optional: boolean;
  /** Default parameters */
  defaultParams: Record<string, unknown>;
  /** Parameter schema for the auto-generated editor */
  parameterSchema: StepParameterField[];
  /** Execute the step. Throws on failure. */
  execute: (ctx: StepExecutionContext) => Promise<StepResult>;
  /** Optional custom result preview component */
  ResultPreview?: React.ComponentType<{
    result: StepResult;
    onSelectDesign?: (pdbContent: string) => void;
    /** Selection props passed when supportsSelection is true and step is paused */
    selectedIds?: string[];
    onSelectionChange?: (ids: string[]) => void;
    isPaused?: boolean;
  }>;
}

// ---- Pipeline Definition ----

/** Complete static definition of a pipeline */
export interface PipelineDefinition {
  /** Unique pipeline ID */
  id: string;
  /** Display name */
  name: string;
  /** Description */
  description: string;
  /** Icon */
  icon: LucideIcon;
  /** Ordered array of steps */
  steps: PipelineStepDefinition[];
  /** Parameters needed before the first step */
  initialParams: StepParameterField[];
  /** Whether this pipeline requires a backend connection */
  requiresBackend: boolean;
}

// ---- Runtime State ----

/** Runtime state for one step */
export interface StepRuntimeState {
  stepId: string;
  status: StepStatus;
  /** Progress within the step (0-100) */
  progress: number;
  /** Current sub-status message */
  progressMessage?: string;
  /** Result once completed */
  result?: StepResult;
  /** Error message if failed */
  error?: string;
  /** IDs of outputs selected by user to carry forward */
  selectedOutputIds: string[];
  /** User-modified parameters for this step */
  params: Record<string, unknown>;
  /** Timestamps */
  startedAt?: number;
  completedAt?: number;
}

/** Overall pipeline status */
export type PipelineStatus = 'idle' | 'running' | 'paused' | 'completed' | 'failed' | 'cancelled';

/** Full runtime state of a pipeline execution */
export interface PipelineRuntimeState {
  /** Which pipeline definition is running */
  pipelineId: string;
  /** Unique execution session ID */
  sessionId: string;
  /** Overall status */
  status: PipelineStatus;
  /** Index of the currently active step */
  activeStepIndex: number;
  /** Per-step runtime states */
  steps: StepRuntimeState[];
  /** Initial parameters provided before pipeline start */
  initialParams: Record<string, unknown>;
  /** Timestamps */
  startedAt?: number;
  completedAt?: number;
}
