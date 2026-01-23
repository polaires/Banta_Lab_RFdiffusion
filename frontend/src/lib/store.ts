/**
 * Zustand store for application state (v2.2)
 * Supports cross-panel data flow, confidence metrics, enhanced viewer state,
 * and versioned persistence with migration support
 */

import { create } from 'zustand';
import { persist, createJSONStorage, StateStorage } from 'zustand/middleware';

// Store version for migration handling
const STORE_VERSION = 2;

// Maximum number of jobs to keep in history
const MAX_JOB_HISTORY = 100;
import type { JobStatus, HealthResponse, ConfidenceMetrics, RMSDResult, MetalBindingAnalysis, UserPreferences, DesignEvaluation } from './api';
import type { MetalCoordination } from './metalAnalysis';
import type { LigandAnalysisResult, PharmacophoreFeature } from './ligandAnalysis';
import type { MetalReplacementPreset, DetectedDonor } from './enzymeAnalysis';

// Viewer mode for different visualization states
export type ViewerMode = 'default' | 'metal' | 'ligand' | 'confidence' | 'comparison';

// Viewer representation styles
export type RepresentationStyle = 'cartoon' | 'ball-and-stick' | 'spacefill' | 'surface';

// Viewer color schemes
export type ColorScheme = 'default' | 'chain' | 'residue-type' | 'secondary-structure' | 'confidence' | 'hydrophobicity';

// Catalytic residue suggestion from API
export interface CatalyticSuggestion {
  chain: string;
  residue: number;
  name: string;
  role?: string;
  confidence: number;
  source: 'mcsa' | 'local';
  ligandCode?: string;  // The ligand/metal this residue is associated with (for filtering)
}

// Bottom panel mode
export type BottomPanelMode = 'suggestions' | 'metal-analysis' | 'ligand-analysis' | 'none';

// Enzyme analysis result from auto-detection (for H-bond/RASA conditioning)
export interface EnzymeAnalysisResult {
  metals: Array<{
    element: string;
    chain: string;
    residueNumber: number;
    coordinatingAtoms: Array<{ chain: string; residue: number; atom: string; distance: number }>;
    coordinationNumber: number;
  }>;
  ligands: Array<{
    name: string;
    chain: string;
    residueNumber: number;
    atoms: Array<{ name: string; element: string }>;
    suggestedBuried: string[];      // Coordination atoms (near metal)
    suggestedExposed: string[];     // Entry/exit atoms (terminal)
    suggestedHBondAcceptors: string[];  // O/N atoms that can accept H-bonds
    suggestedProteinDonors: DetectedDonor[];  // Auto-detected protein donors near ligand
  }>;
}

// Coordination mode for metal replacement workflow
export type CoordinationMode = 'keep' | 'explore' | 'hybrid';

interface ErrorContext {
  task?: string;
  input_keys?: string[];
  gpu_info?: {
    available?: boolean;
    name?: string;
    memory_gb?: number;
  };
  gpu_memory_used_mb?: number;
  gpu_memory_total_mb?: number;
  foundry_available?: boolean;
  checkpoint_dir?: string;
}

interface Job {
  id: string;
  type: 'rfd3' | 'rf3' | 'mpnn';
  status: 'pending' | 'running' | 'completed' | 'failed';
  createdAt: string;
  completedAt?: string;
  result?: JobStatus['result'];
  error?: string;
  errorType?: string;
  traceback?: string;
  errorContext?: ErrorContext;
}

// Tab type for navigation
export type TabId = 'ai' | 'task' | 'rfd3' | 'rf3' | 'mpnn' | 'jobs';

// AI Design Assistant state for case study demo
export interface AIConversationMessage {
  role: 'user' | 'assistant' | 'system';
  content: string;
  timestamp: number;
  data?: {
    type: 'analysis' | 'recommendation' | 'evaluation' | 'pdb_info';
    payload: Record<string, unknown>;
  };
}

// Workflow phase for interview mode
export type WorkflowPhase =
  | 'idle'
  | 'structure_input'
  | 'fetching'
  | 'analyzing'
  | 'interview'
  | 'confirming'
  | 'running'
  | 'evaluating'
  | 'complete'
  | 'error';

export interface AICaseStudyState {
  pdbId: string | null;
  pdbContent: string | null;
  targetMetal: string | null;
  analysisResult: MetalBindingAnalysis | null;
  recommendation: Record<string, unknown> | null;
  conversation: AIConversationMessage[];
  isProcessing: boolean;
  currentStep: 'idle' | 'fetching' | 'analyzing' | 'recommending' | 'complete';
  // Interview mode state
  workflowPhase: WorkflowPhase;
  userPreferences: UserPreferences | null;
  pendingJobId: string | null;
  evaluationResult: DesignEvaluation | null;
  jobProgress: number;
}

// Notification types for workflow guidance
export interface Notification {
  id: string;
  type: 'success' | 'error' | 'info';
  title: string;
  message: string;
  action?: {
    label: string;
    tab: TabId;
  };
  createdAt: number;
}

// Stored structure for cross-panel data flow
export interface StoredDesign {
  jobId: string;
  pdbContent: string;
  cifContent?: string;
  source: 'rfd3' | 'rf3' | 'mpnn';
  sequence?: string;  // Extracted from PDB or MPNN output
  confidences?: ConfidenceMetrics;
  timestamp: number;
}

interface AppState {
  // Backend connection
  backendUrl: string;
  setBackendUrl: (url: string) => void;

  // Health status
  health: HealthResponse | null;
  setHealth: (health: HealthResponse | null) => void;

  // Jobs
  jobs: Job[];
  addJob: (job: Job) => void;
  updateJob: (id: string, updates: Partial<Job>) => void;
  removeJob: (id: string) => void;
  cleanupExpiredJobs: () => void;

  // Selected structure for visualization
  selectedPdb: string | null;
  setSelectedPdb: (pdb: string | null) => void;

  // Active tab
  activeTab: TabId;
  setActiveTab: (tab: TabId) => void;

  // Selected design task (persisted between task and rfd3 tabs)
  selectedDesignTask: string | null;
  setSelectedDesignTask: (task: string | null) => void;

  // Notifications for workflow guidance
  notifications: Notification[];
  addNotification: (notification: Omit<Notification, 'id' | 'createdAt'>) => void;
  dismissNotification: (id: string) => void;
  clearNotifications: () => void;

  // Panel data flow - latest design for MPNN input
  latestDesignPdb: string | null;
  setLatestDesignPdb: (pdb: string | null) => void;
  lastCompletedJobType: 'rfd3' | 'rf3' | 'mpnn' | null;
  setLastCompletedJobType: (type: 'rfd3' | 'rf3' | 'mpnn' | null) => void;

  // Enhanced cross-panel data flow (v2.0)
  latestRfd3Design: StoredDesign | null;
  setLatestRfd3Design: (design: StoredDesign | null) => void;
  latestRf3Prediction: StoredDesign | null;
  setLatestRf3Prediction: (design: StoredDesign | null) => void;

  // Confidence metrics from RF3
  latestConfidences: ConfidenceMetrics | null;
  setLatestConfidences: (confidences: ConfidenceMetrics | null) => void;

  // RMSD validation result
  latestRmsdResult: RMSDResult | null;
  setLatestRmsdResult: (result: RMSDResult | null) => void;

  // UI state
  viewerCollapsed: boolean;
  setViewerCollapsed: (collapsed: boolean) => void;
  connectionModalOpen: boolean;
  setConnectionModalOpen: (open: boolean) => void;

  // Enhanced viewer state (v2.1)
  viewerMode: ViewerMode;
  setViewerMode: (mode: ViewerMode) => void;
  representationStyle: RepresentationStyle;
  setRepresentationStyle: (style: RepresentationStyle) => void;
  colorScheme: ColorScheme;
  setColorScheme: (scheme: ColorScheme) => void;

  // Metal coordination analysis results
  metalCoordination: MetalCoordination[] | null;
  setMetalCoordination: (data: MetalCoordination[] | null) => void;

  // Ligand analysis results
  ligandData: LigandAnalysisResult | null;
  setLigandData: (data: LigandAnalysisResult | null) => void;

  // Pharmacophore visualization state
  pharmacophoreFeatures: PharmacophoreFeature[] | null;
  setPharmacophoreFeatures: (features: PharmacophoreFeature[] | null) => void;
  showPharmacophores3D: boolean;
  setShowPharmacophores3D: (show: boolean) => void;

  // Structure comparison state
  comparisonEnabled: boolean;
  setComparisonEnabled: (enabled: boolean) => void;
  comparisonMode: 'overlay' | 'side-by-side';
  setComparisonMode: (mode: 'overlay' | 'side-by-side') => void;
  referenceStructure: { pdb: string; label: string; source: 'rfd3' | 'rf3' | 'upload' } | null;
  setReferenceStructure: (struct: { pdb: string; label: string; source: 'rfd3' | 'rf3' | 'upload' } | null) => void;

  // Focused binding site indices
  focusedMetalIndex: number | null;
  setFocusedMetalIndex: (index: number | null) => void;
  focusedLigandIndex: number | null;
  setFocusedLigandIndex: (index: number | null) => void;

  // Focus mode visualization settings
  focusSettings: {
    coordinationRadius: number;  // Metal coordination sphere radius (default 3.0)
    bindingPocketRadius: number; // Ligand binding pocket radius (default 5.0)
    showWaters: boolean;         // Show coordinating waters in focus view
    showInteractionLines: boolean; // Show H-bond, salt bridge lines
    showPharmacophores: boolean;  // Show pharmacophore spheres during focus
    ligandCarbonColor: number;    // Ligand carbon color (default green 0x50C878)
  };
  setFocusSettings: (settings: Partial<AppState['focusSettings']>) => void;

  // Analysis loading state
  analysisLoading: boolean;
  setAnalysisLoading: (loading: boolean) => void;

  // AI Design Assistant state
  aiCaseStudy: AICaseStudyState;
  setAiCaseStudy: (state: Partial<AICaseStudyState>) => void;
  addAiMessage: (message: Omit<AIConversationMessage, 'timestamp'>) => void;
  clearAiConversation: () => void;

  // Manual mode toggle (sidebar)
  manualMode: boolean;
  setManualMode: (enabled: boolean) => void;

  // Catalytic residue suggestions
  catalyticSuggestions: CatalyticSuggestion[];
  suggestionsSource: 'mcsa' | 'local' | 'none';
  suggestionsLoading: boolean;
  suggestionsError: string | null;
  bottomPanelMode: BottomPanelMode;
  /** Currently hovered catalytic suggestion for 3D highlight */
  hoveredCatalyticSuggestion: CatalyticSuggestion | null;

  // Catalytic suggestions actions
  setCatalyticSuggestions: (suggestions: CatalyticSuggestion[], source: 'mcsa' | 'local' | 'none') => void;
  setSuggestionsLoading: (loading: boolean) => void;
  setSuggestionsError: (error: string | null) => void;
  setBottomPanelMode: (mode: BottomPanelMode) => void;
  setHoveredCatalyticSuggestion: (suggestion: CatalyticSuggestion | null) => void;
  clearSuggestions: () => void;

  // Enzyme form shared state (for cross-component communication)
  enzymeCatalyticResidues: Array<{ chain: string; residue: number; name: string }>;
  enzymeFixedAtomTypes: Record<string, string>;
  /** Current ligand codes for filtering catalytic suggestions */
  enzymeLigandCodes: string;
  addEnzymeCatalyticResidue: (chain: string, residue: number, name: string, atomType: string) => void;
  removeEnzymeCatalyticResidue: (chain: string, residue: number) => void;
  clearEnzymeCatalyticResidues: () => void;
  /** Batch set catalytic residues with their atom types */
  setEnzymeCatalyticResidues: (residues: Array<{ chain: string; residue: number; name: string; atomType: string }>) => void;
  /** Set ligand codes for filtering catalytic suggestions */
  setEnzymeLigandCodes: (codes: string) => void;
  /** Set fixed atom types for all residues */
  setEnzymeFixedAtomTypes: (types: Record<string, string>) => void;

  // Enzyme analysis state (auto-detection results)
  enzymeAnalysis: EnzymeAnalysisResult | null;
  enzymeAnalysisLoading: boolean;
  setEnzymeAnalysis: (analysis: EnzymeAnalysisResult | null) => void;
  setEnzymeAnalysisLoading: (loading: boolean) => void;

  // Metal replacement workflow
  metalReplacementEnabled: boolean;
  sourceMetal: string | null;
  targetMetal: string | null;
  coordinationMode: CoordinationMode;
  fixLigandPosition: boolean;
  metalReplacementPreset: MetalReplacementPreset;
  /** Residues excluded from fixing due to preset (metal coordinators for expand/lanthanide) */
  excludedCatalyticResidues: Set<string>;
  setMetalReplacementEnabled: (enabled: boolean) => void;
  setSourceMetal: (metal: string | null) => void;
  setTargetMetal: (metal: string | null) => void;
  setCoordinationMode: (mode: CoordinationMode) => void;
  setFixLigandPosition: (fix: boolean) => void;
  setMetalReplacementPreset: (preset: MetalReplacementPreset) => void;
  /** Apply a preset - configures all parameters based on preset config */
  applyMetalReplacementPreset: (preset: MetalReplacementPreset) => void;

  // RASA conditioning (burial/exposure)
  selectedBuriedAtoms: Record<string, string>;    // {"CIT": "O3,O4,O5", "TB": "ALL"}
  selectedExposedAtoms: Record<string, string>;   // {"CIT": "O1,O2,O6,O7"}
  buriedOverridden: boolean;
  exposedOverridden: boolean;
  setBuriedAtoms: (key: string, atoms: string) => void;
  setExposedAtoms: (key: string, atoms: string) => void;
  removeBuriedAtoms: (key: string) => void;
  removeExposedAtoms: (key: string) => void;
  clearRASAConditioning: () => void;

  // H-Bond conditioning
  selectedHBondAcceptors: Record<string, string>; // {"CIT": "O1,O2,O3,O4,O5,O6,O7"}
  selectedHBondDonors: Record<string, string>;    // {"A193": "N", "A195": "N"}
  hbondOverridden: boolean;
  setHBondAcceptors: (key: string, atoms: string) => void;
  setHBondDonors: (key: string, atoms: string) => void;
  removeHBondAcceptors: (key: string) => void;
  removeHBondDonors: (key: string) => void;
  clearHBondConditioning: () => void;

  // Apply auto-suggestions from analysis
  applyEnzymeSuggestions: () => void;
  clearAllEnzymeConditioning: () => void;

  // Pipeline state for parameter sweeps and production runs
  pipelineState: PipelineState;
  setPipelineMode: (mode: PipelineMode) => void;
  startPipeline: (sessionId: string, mode: PipelineMode, totalConfigs: number, designsPerConfig: number) => void;
  updatePipelineProgress: (progress: Partial<PipelineProgress>) => void;
  setPipelineResults: (results: PipelineDesign[]) => void;
  setPipelineFilters: (filters: Partial<PipelineFilters>) => void;
  setSweepConfigs: (configs: SweepConfig[]) => void;
  resetPipeline: () => void;
  cancelPipeline: () => void;
}

// Pipeline types
export type PipelineMode = 'single' | 'sweep' | 'production';

export interface PipelineFilters {
  plddt: number;
  ptm: number;
  pae: number;
}

export interface SweepConfig {
  name: string;
  contigSize: 'small' | 'medium' | 'large';
  contigRange: string;
  cfgScale: number;
  numDesigns: number;
}

export interface PipelineDesign {
  name: string;
  config: string;
  sequence: string;
  plddt: number;
  ptm: number;
  pae: number;
  status: 'pass' | 'review' | 'fail';
  timestamp?: string;
}

export interface PipelineProgress {
  currentConfig: number;
  totalConfigs: number;
  currentDesign: number;
  designsPerConfig: number;
  totalGenerated: number;
  totalPassing: number;
  totalReview: number;
  totalFailed: number;
  passRate: number;
  bestDesign: {
    name: string;
    plddt: number;
    ptm: number;
    pae: number;
  } | null;
}

export interface PipelineState {
  mode: PipelineMode;
  isRunning: boolean;
  sessionId: string | null;
  progress: PipelineProgress;
  results: PipelineDesign[];
  filters: PipelineFilters;
  sweepConfigs: SweepConfig[];
}

// Default backend URL from environment or fallback
const DEFAULT_BACKEND_URL = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000';

export const useStore = create<AppState>()(
  persist(
    (set) => ({
  // Backend connection - persisted via Zustand middleware
  backendUrl: DEFAULT_BACKEND_URL,
  setBackendUrl: (url) => set({ backendUrl: url }),

  // Health status
  health: null,
  setHealth: (health) => set({ health }),

  // Jobs
  jobs: [],
  addJob: (job) => set((state) => {
    // Add new job at the start, then enforce max limit
    const newJobs = [job, ...state.jobs];
    // Keep only the most recent MAX_JOB_HISTORY jobs
    if (newJobs.length > MAX_JOB_HISTORY) {
      return { jobs: newJobs.slice(0, MAX_JOB_HISTORY) };
    }
    return { jobs: newJobs };
  }),
  updateJob: (id, updates) => set((state) => ({
    jobs: state.jobs.map((job) => job.id === id ? { ...job, ...updates } : job)
  })),
  removeJob: (id) => set((state) => ({
    jobs: state.jobs.filter((job) => job.id !== id)
  })),
  cleanupExpiredJobs: () => set((state) => {
    const now = Date.now();
    const TWENTY_FOUR_HOURS = 24 * 60 * 60 * 1000;
    return {
      jobs: state.jobs.filter((job) => {
        // Keep all non-failed jobs
        if (job.status !== 'failed') return true;
        // Remove failed jobs older than 24 hours
        const createdAt = new Date(job.createdAt).getTime();
        const expiresAt = createdAt + TWENTY_FOUR_HOURS;
        return now < expiresAt;
      }),
    };
  }),

  // Selected structure
  selectedPdb: null,
  setSelectedPdb: (pdb) => set({ selectedPdb: pdb }),

  // Active tab - start at AI assistant
  activeTab: 'ai',
  setActiveTab: (tab) => set({ activeTab: tab }),

  // Selected design task
  selectedDesignTask: null,
  setSelectedDesignTask: (task) => set({ selectedDesignTask: task }),

  // Notifications
  notifications: [],
  addNotification: (notification) => set((state) => ({
    notifications: [
      ...state.notifications,
      {
        ...notification,
        id: `notif-${Date.now()}-${Math.random().toString(36).slice(2, 9)}`,
        createdAt: Date.now(),
      },
    ],
  })),
  dismissNotification: (id) => set((state) => ({
    notifications: state.notifications.filter((n) => n.id !== id),
  })),
  clearNotifications: () => set({ notifications: [] }),

  // Panel data flow
  latestDesignPdb: null,
  setLatestDesignPdb: (pdb) => set({ latestDesignPdb: pdb }),
  lastCompletedJobType: null,
  setLastCompletedJobType: (type) => set({ lastCompletedJobType: type }),

  // Enhanced cross-panel data flow (v2.0)
  latestRfd3Design: null,
  setLatestRfd3Design: (design) => set({ latestRfd3Design: design }),
  latestRf3Prediction: null,
  setLatestRf3Prediction: (design) => set({ latestRf3Prediction: design }),

  // Confidence metrics
  latestConfidences: null,
  setLatestConfidences: (confidences) => set({ latestConfidences: confidences }),

  // RMSD validation
  latestRmsdResult: null,
  setLatestRmsdResult: (result) => set({ latestRmsdResult: result }),

  // UI state
  viewerCollapsed: false,
  setViewerCollapsed: (collapsed) => set({ viewerCollapsed: collapsed }),
  connectionModalOpen: false,
  setConnectionModalOpen: (open) => set({ connectionModalOpen: open }),

  // Enhanced viewer state (v2.1)
  viewerMode: 'default',
  setViewerMode: (mode) => set({ viewerMode: mode }),
  representationStyle: 'cartoon',
  setRepresentationStyle: (style) => set({ representationStyle: style }),
  colorScheme: 'default',
  setColorScheme: (scheme) => set({ colorScheme: scheme }),

  // Metal coordination analysis
  metalCoordination: null,
  setMetalCoordination: (data) => set({ metalCoordination: data }),

  // Ligand analysis
  ligandData: null,
  setLigandData: (data) => set({ ligandData: data }),

  // Pharmacophore visualization state
  pharmacophoreFeatures: null,
  setPharmacophoreFeatures: (features) => set({ pharmacophoreFeatures: features }),
  showPharmacophores3D: false,
  setShowPharmacophores3D: (show) => set({ showPharmacophores3D: show }),

  // Structure comparison
  comparisonEnabled: false,
  setComparisonEnabled: (enabled) => set({ comparisonEnabled: enabled }),
  comparisonMode: 'overlay',
  setComparisonMode: (mode) => set({ comparisonMode: mode }),
  referenceStructure: null,
  setReferenceStructure: (struct) => set({ referenceStructure: struct }),

  // Focused binding site indices
  focusedMetalIndex: null,
  setFocusedMetalIndex: (index) => set({ focusedMetalIndex: index }),
  focusedLigandIndex: null,
  setFocusedLigandIndex: (index) => set({ focusedLigandIndex: index }),

  // Focus mode visualization settings
  focusSettings: {
    coordinationRadius: 3.0,
    bindingPocketRadius: 5.0,
    showWaters: false,
    showInteractionLines: true,
    showPharmacophores: false,
    ligandCarbonColor: 0x50C878,
  },
  setFocusSettings: (settings) => set((state) => ({
    focusSettings: { ...state.focusSettings, ...settings },
  })),

  // Analysis loading state
  analysisLoading: false,
  setAnalysisLoading: (loading) => set({ analysisLoading: loading }),

  // AI Design Assistant state
  aiCaseStudy: {
    pdbId: null,
    pdbContent: null,
    targetMetal: null,
    analysisResult: null,
    recommendation: null,
    conversation: [],
    isProcessing: false,
    currentStep: 'idle',
    // Interview mode state
    workflowPhase: 'idle',
    userPreferences: null,
    pendingJobId: null,
    evaluationResult: null,
    jobProgress: 0,
  },
  setAiCaseStudy: (updates) => set((state) => ({
    aiCaseStudy: { ...state.aiCaseStudy, ...updates },
  })),
  addAiMessage: (message) => set((state) => ({
    aiCaseStudy: {
      ...state.aiCaseStudy,
      conversation: [
        ...state.aiCaseStudy.conversation,
        { ...message, timestamp: Date.now() },
      ],
    },
  })),
  clearAiConversation: () => set((state) => ({
    aiCaseStudy: {
      ...state.aiCaseStudy,
      conversation: [],
      pdbId: null,
      pdbContent: null,
      targetMetal: null,
      analysisResult: null,
      recommendation: null,
      currentStep: 'idle',
      // Reset interview state
      workflowPhase: 'idle',
      userPreferences: null,
      pendingJobId: null,
      evaluationResult: null,
      jobProgress: 0,
    },
  })),

  // Manual mode toggle (sidebar)
  manualMode: false,
  setManualMode: (enabled) => set({ manualMode: enabled }),

  // Catalytic residue suggestions
  catalyticSuggestions: [],
  suggestionsSource: 'none',
  suggestionsLoading: false,
  suggestionsError: null,
  bottomPanelMode: 'none',
  hoveredCatalyticSuggestion: null,

  setCatalyticSuggestions: (suggestions, source) => set({
    catalyticSuggestions: suggestions,
    suggestionsSource: source,
    suggestionsError: null,
    // Auto-show panel if suggestions found
    bottomPanelMode: suggestions.length > 0 ? 'suggestions' : 'none',
  }),
  setSuggestionsLoading: (loading) => set({ suggestionsLoading: loading }),
  setSuggestionsError: (error) => set({ suggestionsError: error, suggestionsLoading: false }),
  setBottomPanelMode: (mode) => set({ bottomPanelMode: mode }),
  setHoveredCatalyticSuggestion: (suggestion) => set({ hoveredCatalyticSuggestion: suggestion }),
  clearSuggestions: () => set({
    catalyticSuggestions: [],
    suggestionsSource: 'none',
    suggestionsError: null,
    bottomPanelMode: 'none',
  }),

  // Enzyme form shared state
  enzymeCatalyticResidues: [],
  enzymeFixedAtomTypes: {},
  enzymeLigandCodes: '',
  addEnzymeCatalyticResidue: (chain, residue, name, atomType) => set((state) => {
    const key = `${chain}${residue}`;
    if (state.enzymeCatalyticResidues.some((r) => r.chain === chain && r.residue === residue)) {
      return state;
    }
    return {
      enzymeCatalyticResidues: [...state.enzymeCatalyticResidues, { chain, residue, name }],
      enzymeFixedAtomTypes: { ...state.enzymeFixedAtomTypes, [key]: atomType },
    };
  }),
  removeEnzymeCatalyticResidue: (chain, residue) => set((state) => {
    const key = `${chain}${residue}`;
    const { [key]: _, ...restAtomTypes } = state.enzymeFixedAtomTypes;
    return {
      enzymeCatalyticResidues: state.enzymeCatalyticResidues.filter(
        (r) => !(r.chain === chain && r.residue === residue)
      ),
      enzymeFixedAtomTypes: restAtomTypes,
    };
  }),
  clearEnzymeCatalyticResidues: () => set({
    enzymeCatalyticResidues: [],
    enzymeFixedAtomTypes: {},
  }),
  setEnzymeCatalyticResidues: (residues) => set(() => {
    const newResidues: Array<{ chain: string; residue: number; name: string }> = [];
    const newAtomTypes: Record<string, string> = {};

    for (const r of residues) {
      // Avoid duplicates
      if (!newResidues.some(existing => existing.chain === r.chain && existing.residue === r.residue)) {
        newResidues.push({ chain: r.chain, residue: r.residue, name: r.name });
        newAtomTypes[`${r.chain}${r.residue}`] = r.atomType;
      }
    }

    return {
      enzymeCatalyticResidues: newResidues,
      enzymeFixedAtomTypes: newAtomTypes,
    };
  }),
  setEnzymeLigandCodes: (codes) => set({ enzymeLigandCodes: codes }),
  setEnzymeFixedAtomTypes: (types) => set({ enzymeFixedAtomTypes: types }),

  // Enzyme analysis state
  enzymeAnalysis: null,
  enzymeAnalysisLoading: false,
  setEnzymeAnalysis: (analysis) => set({ enzymeAnalysis: analysis }),
  setEnzymeAnalysisLoading: (loading) => set({ enzymeAnalysisLoading: loading }),

  // Metal replacement workflow
  metalReplacementEnabled: false,
  sourceMetal: null,
  targetMetal: null,
  coordinationMode: 'explore',
  fixLigandPosition: true,
  metalReplacementPreset: 'custom' as MetalReplacementPreset,
  excludedCatalyticResidues: new Set<string>(),
  setMetalReplacementEnabled: (enabled) => set({ metalReplacementEnabled: enabled }),
  setSourceMetal: (metal) => set({ sourceMetal: metal }),
  setTargetMetal: (metal) => set({ targetMetal: metal }),
  setCoordinationMode: (mode) => set({ coordinationMode: mode }),
  setFixLigandPosition: (fix) => set({ fixLigandPosition: fix }),
  setMetalReplacementPreset: (preset) => set({ metalReplacementPreset: preset }),

  // Apply a metal replacement preset - configures ALL parameters including catalytic residues
  applyMetalReplacementPreset: (preset) => set((state) => {
    // Import preset config dynamically to avoid circular deps
    const { METAL_REPLACEMENT_PRESETS, getPresetConfig } = require('./enzymeAnalysis');
    const config = getPresetConfig(preset, state.sourceMetal, state.targetMetal);

    // Determine which residues to exclude (metal coordinators for expand/lanthanide presets)
    const excludedResidues = new Set<string>();
    const metalCoordinators: Array<{ chain: string; residue: number; name: string; atom: string }> = [];

    if (state.enzymeAnalysis) {
      // Collect metal coordinators
      for (const metal of state.enzymeAnalysis.metals) {
        for (const coord of metal.coordinatingAtoms) {
          metalCoordinators.push({
            chain: coord.chain,
            residue: coord.residue,
            name: '', // Will be filled from suggestions or defaults
            atom: coord.atom,
          });
          // Mark as excluded if preset doesn't fix metal coordinators
          if (!config.fixMetalCoordinators) {
            excludedResidues.add(`${coord.chain}${coord.residue}`);
          }
        }
      }
    }

    // Build catalytic residues from suggestions and metal coordinators
    const newCatalyticResidues: Array<{ chain: string; residue: number; name: string }> = [];
    const newAtomTypes: Record<string, string> = {};

    // Add residues from catalytic suggestions (MCSA/local analysis)
    for (const suggestion of state.catalyticSuggestions) {
      const key = `${suggestion.chain}${suggestion.residue}`;
      // Skip if excluded (metal coordinators for expand/lanthanide)
      if (excludedResidues.has(key)) continue;

      if (!newCatalyticResidues.some(r => r.chain === suggestion.chain && r.residue === suggestion.residue)) {
        newCatalyticResidues.push({
          chain: suggestion.chain,
          residue: suggestion.residue,
          name: suggestion.name,
        });
        // Use 'ALL' to fix entire residue (valid RFD3 values: BKBN, ALL, TIP, '')
        newAtomTypes[key] = 'ALL';
      }
    }

    // Add metal coordinators if preset fixes them
    if (config.fixMetalCoordinators) {
      for (const coord of metalCoordinators) {
        const key = `${coord.chain}${coord.residue}`;
        if (!newCatalyticResidues.some(r => r.chain === coord.chain && r.residue === coord.residue)) {
          // Try to find name from suggestions
          const suggestion = state.catalyticSuggestions.find(
            s => s.chain === coord.chain && s.residue === coord.residue
          );
          newCatalyticResidues.push({
            chain: coord.chain,
            residue: coord.residue,
            name: suggestion?.name || 'UNK',
          });
          // Metal coordinators should fix entire residue (valid RFD3 values: BKBN, ALL, TIP, '')
          newAtomTypes[key] = 'ALL';
        }
      }
    }

    // Build buried atoms if preset enables RASA and burial
    const newBuried: Record<string, string> = {};
    if (config.buryMetal && state.targetMetal) {
      newBuried[state.targetMetal] = 'ALL';
    }
    // Also bury ligand coordination atoms if analysis available
    if (config.enableRASA && state.enzymeAnalysis) {
      for (const ligand of state.enzymeAnalysis.ligands) {
        if (ligand.suggestedBuried.length > 0) {
          newBuried[ligand.name] = ligand.suggestedBuried.join(',');
        }
      }
    }

    // Build exposed atoms
    const newExposed: Record<string, string> = {};
    if (config.enableRASA && state.enzymeAnalysis) {
      for (const ligand of state.enzymeAnalysis.ligands) {
        if (ligand.suggestedExposed.length > 0) {
          newExposed[ligand.name] = ligand.suggestedExposed.join(',');
        }
      }
    }

    // Build H-bond acceptors
    const newAcceptors: Record<string, string> = {};
    if (config.enableHBonds && state.enzymeAnalysis) {
      for (const ligand of state.enzymeAnalysis.ligands) {
        if (ligand.suggestedHBondAcceptors.length > 0) {
          newAcceptors[ligand.name] = ligand.suggestedHBondAcceptors.join(',');
        }
      }
    }

    return {
      metalReplacementPreset: preset,
      coordinationMode: config.coordinationMode,
      fixLigandPosition: config.fixLigandPosition,
      excludedCatalyticResidues: excludedResidues,
      // Set catalytic residues from combined suggestions + coordinators
      enzymeCatalyticResidues: newCatalyticResidues,
      enzymeFixedAtomTypes: newAtomTypes,
      // Apply RASA if enabled by preset
      selectedBuriedAtoms: config.enableRASA ? newBuried : state.selectedBuriedAtoms,
      selectedExposedAtoms: config.enableRASA ? newExposed : state.selectedExposedAtoms,
      // Apply H-bonds if enabled by preset
      selectedHBondAcceptors: config.enableHBonds ? newAcceptors : state.selectedHBondAcceptors,
      // Reset override flags
      buriedOverridden: false,
      exposedOverridden: false,
      hbondOverridden: false,
    };
  }),

  // RASA conditioning
  selectedBuriedAtoms: {},
  selectedExposedAtoms: {},
  buriedOverridden: false,
  exposedOverridden: false,
  setBuriedAtoms: (key, atoms) => set((state) => ({
    selectedBuriedAtoms: { ...state.selectedBuriedAtoms, [key]: atoms },
    buriedOverridden: true,
  })),
  setExposedAtoms: (key, atoms) => set((state) => ({
    selectedExposedAtoms: { ...state.selectedExposedAtoms, [key]: atoms },
    exposedOverridden: true,
  })),
  removeBuriedAtoms: (key) => set((state) => {
    const { [key]: _, ...rest } = state.selectedBuriedAtoms;
    return { selectedBuriedAtoms: rest };
  }),
  removeExposedAtoms: (key) => set((state) => {
    const { [key]: _, ...rest } = state.selectedExposedAtoms;
    return { selectedExposedAtoms: rest };
  }),
  clearRASAConditioning: () => set({
    selectedBuriedAtoms: {},
    selectedExposedAtoms: {},
    buriedOverridden: false,
    exposedOverridden: false,
  }),

  // H-Bond conditioning
  selectedHBondAcceptors: {},
  selectedHBondDonors: {},
  hbondOverridden: false,
  setHBondAcceptors: (key, atoms) => set((state) => ({
    selectedHBondAcceptors: { ...state.selectedHBondAcceptors, [key]: atoms },
    hbondOverridden: true,
  })),
  setHBondDonors: (key, atoms) => set((state) => ({
    selectedHBondDonors: { ...state.selectedHBondDonors, [key]: atoms },
    hbondOverridden: true,
  })),
  removeHBondAcceptors: (key) => set((state) => {
    const { [key]: _, ...rest } = state.selectedHBondAcceptors;
    return { selectedHBondAcceptors: rest };
  }),
  removeHBondDonors: (key) => set((state) => {
    const { [key]: _, ...rest } = state.selectedHBondDonors;
    return { selectedHBondDonors: rest };
  }),
  clearHBondConditioning: () => set({
    selectedHBondAcceptors: {},
    selectedHBondDonors: {},
    hbondOverridden: false,
  }),

  // Apply auto-suggestions from analysis
  // NOTE: This function ALWAYS applies suggestions regardless of override flags,
  // because it's called when user explicitly clicks "Apply Suggestions" button.
  // Override flags are only used to prevent auto-application on PDB upload.
  applyEnzymeSuggestions: () => set((state) => {
    if (!state.enzymeAnalysis) return state;

    const newBuried: Record<string, string> = {};
    const newExposed: Record<string, string> = {};
    const newAcceptors: Record<string, string> = {};

    // Apply metal burial - use TARGET metal if replacement is enabled, otherwise source
    for (const metal of state.enzymeAnalysis.metals) {
      const metalToUse = (state.metalReplacementEnabled && state.targetMetal)
        ? state.targetMetal
        : metal.element;
      newBuried[metalToUse] = 'ALL';
    }

    // Apply ligand suggestions (always apply - user clicked the button)
    for (const ligand of state.enzymeAnalysis.ligands) {
      if (ligand.suggestedBuried.length > 0) {
        newBuried[ligand.name] = ligand.suggestedBuried.join(',');
      }
      if (ligand.suggestedExposed.length > 0) {
        newExposed[ligand.name] = ligand.suggestedExposed.join(',');
      }
      if (ligand.suggestedHBondAcceptors.length > 0) {
        newAcceptors[ligand.name] = ligand.suggestedHBondAcceptors.join(',');
      }
    }

    // Always apply suggestions and reset override flags
    return {
      selectedBuriedAtoms: newBuried,
      selectedExposedAtoms: newExposed,
      selectedHBondAcceptors: newAcceptors,
      // Reset override flags so auto-suggestion works on next PDB upload
      buriedOverridden: false,
      exposedOverridden: false,
      hbondOverridden: false,
    };
  }),

  clearAllEnzymeConditioning: () => set({
    enzymeAnalysis: null,
    metalReplacementEnabled: false,
    sourceMetal: null,
    targetMetal: null,
    coordinationMode: 'explore',
    fixLigandPosition: true,
    metalReplacementPreset: 'custom' as MetalReplacementPreset,
    excludedCatalyticResidues: new Set<string>(),
    selectedBuriedAtoms: {},
    selectedExposedAtoms: {},
    buriedOverridden: false,
    exposedOverridden: false,
    selectedHBondAcceptors: {},
    selectedHBondDonors: {},
    hbondOverridden: false,
  }),

  // Pipeline state for parameter sweeps and production runs
  pipelineState: {
    mode: 'single' as PipelineMode,
    isRunning: false,
    sessionId: null,
    progress: {
      currentConfig: 0,
      totalConfigs: 0,
      currentDesign: 0,
      designsPerConfig: 10,
      totalGenerated: 0,
      totalPassing: 0,
      totalReview: 0,
      totalFailed: 0,
      passRate: 0,
      bestDesign: null,
    },
    results: [],
    filters: {
      plddt: 0.80,
      ptm: 0.80,
      pae: 5.0,
    },
    sweepConfigs: [],
  },

  setPipelineMode: (mode) => set((state) => ({
    pipelineState: { ...state.pipelineState, mode },
  })),

  startPipeline: (sessionId, mode, totalConfigs, designsPerConfig) => set((state) => ({
    pipelineState: {
      ...state.pipelineState,
      mode,
      isRunning: true,
      sessionId,
      progress: {
        ...state.pipelineState.progress,
        currentConfig: 0,
        totalConfigs,
        currentDesign: 0,
        designsPerConfig,
        totalGenerated: 0,
        totalPassing: 0,
        totalReview: 0,
        totalFailed: 0,
        passRate: 0,
        bestDesign: null,
      },
      results: [],
    },
  })),

  updatePipelineProgress: (progress) => set((state) => ({
    pipelineState: {
      ...state.pipelineState,
      progress: { ...state.pipelineState.progress, ...progress },
    },
  })),

  setPipelineResults: (results) => set((state) => ({
    pipelineState: {
      ...state.pipelineState,
      isRunning: false,
      results,
    },
  })),

  setPipelineFilters: (filters) => set((state) => ({
    pipelineState: {
      ...state.pipelineState,
      filters: { ...state.pipelineState.filters, ...filters },
    },
  })),

  setSweepConfigs: (configs) => set((state) => ({
    pipelineState: {
      ...state.pipelineState,
      sweepConfigs: configs,
    },
  })),

  resetPipeline: () => set((state) => ({
    pipelineState: {
      mode: 'single' as PipelineMode,
      isRunning: false,
      sessionId: null,
      progress: {
        currentConfig: 0,
        totalConfigs: 0,
        currentDesign: 0,
        designsPerConfig: 10,
        totalGenerated: 0,
        totalPassing: 0,
        totalReview: 0,
        totalFailed: 0,
        passRate: 0,
        bestDesign: null,
      },
      results: [],
      filters: state.pipelineState.filters, // Keep filters
      sweepConfigs: [], // Reset configs
    },
  })),

  cancelPipeline: () => set((state) => ({
    pipelineState: {
      ...state.pipelineState,
      isRunning: false,
    },
  })),
}),
    {
      name: 'rfd3-design-history',
      version: STORE_VERSION,
      storage: createJSONStorage(() => localStorage),
      partialize: (state) => ({
        // Only persist jobs and backend URL - not transient UI state
        jobs: state.jobs,
        backendUrl: state.backendUrl,
      }),
      migrate: (persistedState: unknown, version: number) => {
        // Handle migrations from older versions
        const state = persistedState as Partial<AppState>;

        if (version < 2) {
          // v1 -> v2: Ensure jobs array exists and has required fields
          if (state.jobs) {
            state.jobs = state.jobs.map((job: Job) => ({
              ...job,
              // Ensure errorContext exists for older jobs
              errorContext: job.errorContext || undefined,
            }));
          }
        }

        // Enforce max job limit on load
        if (state.jobs && state.jobs.length > MAX_JOB_HISTORY) {
          state.jobs = state.jobs.slice(0, MAX_JOB_HISTORY);
        }

        return state as AppState;
      },
    }
  )
);
