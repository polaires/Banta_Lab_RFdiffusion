/**
 * Zustand store for application state (v2.1)
 * Supports cross-panel data flow, confidence metrics, and enhanced viewer state
 */

import { create } from 'zustand';
import { persist, createJSONStorage } from 'zustand/middleware';
import type { JobStatus, HealthResponse, ConfidenceMetrics, RMSDResult, MetalBindingAnalysis, UserPreferences, DesignEvaluation } from './api';
import type { MetalCoordination } from './metalAnalysis';
import type { LigandAnalysisResult } from './ligandAnalysis';

// Viewer mode for different visualization states
export type ViewerMode = 'default' | 'metal' | 'ligand' | 'confidence' | 'comparison';

// Viewer representation styles
export type RepresentationStyle = 'cartoon' | 'ball-and-stick' | 'spacefill' | 'surface';

// Viewer color schemes
export type ColorScheme = 'default' | 'chain' | 'residue-type' | 'secondary-structure' | 'confidence' | 'hydrophobicity';

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

  // Structure comparison state
  comparisonEnabled: boolean;
  setComparisonEnabled: (enabled: boolean) => void;
  referenceStructure: { pdb: string; label: string; source: 'rfd3' | 'rf3' | 'upload' } | null;
  setReferenceStructure: (struct: { pdb: string; label: string; source: 'rfd3' | 'rf3' | 'upload' } | null) => void;

  // Focused binding site indices
  focusedMetalIndex: number | null;
  setFocusedMetalIndex: (index: number | null) => void;
  focusedLigandIndex: number | null;
  setFocusedLigandIndex: (index: number | null) => void;

  // Analysis loading state
  analysisLoading: boolean;
  setAnalysisLoading: (loading: boolean) => void;

  // AI Design Assistant state
  aiCaseStudy: AICaseStudyState;
  setAiCaseStudy: (state: Partial<AICaseStudyState>) => void;
  addAiMessage: (message: Omit<AIConversationMessage, 'timestamp'>) => void;
  clearAiConversation: () => void;
}

// Helper to get initial backend URL (localStorage > env > default)
const getInitialBackendUrl = (): string => {
  if (typeof window !== 'undefined') {
    const stored = localStorage.getItem('foundry-backend-url');
    if (stored) return stored;
  }
  return process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000';
};

export const useStore = create<AppState>((set) => ({
  // Backend connection - persisted to localStorage
  backendUrl: getInitialBackendUrl(),
  setBackendUrl: (url) => {
    // Persist to localStorage
    if (typeof window !== 'undefined') {
      localStorage.setItem('foundry-backend-url', url);
    }
    set({ backendUrl: url });
  },

  // Health status
  health: null,
  setHealth: (health) => set({ health }),

  // Jobs
  jobs: [],
  addJob: (job) => set((state) => ({ jobs: [job, ...state.jobs] })),
  updateJob: (id, updates) => set((state) => ({
    jobs: state.jobs.map((job) => job.id === id ? { ...job, ...updates } : job)
  })),
  removeJob: (id) => set((state) => ({
    jobs: state.jobs.filter((job) => job.id !== id)
  })),

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

  // Structure comparison
  comparisonEnabled: false,
  setComparisonEnabled: (enabled) => set({ comparisonEnabled: enabled }),
  referenceStructure: null,
  setReferenceStructure: (struct) => set({ referenceStructure: struct }),

  // Focused binding site indices
  focusedMetalIndex: null,
  setFocusedMetalIndex: (index) => set({ focusedMetalIndex: index }),
  focusedLigandIndex: null,
  setFocusedLigandIndex: (index) => set({ focusedLigandIndex: index }),

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
}));
