/**
 * Zustand store for application state (v2.0)
 * Supports cross-panel data flow and confidence metrics
 */

import { create } from 'zustand';
import { persist, createJSONStorage } from 'zustand/middleware';
import type { JobStatus, HealthResponse, ConfidenceMetrics, RMSDResult } from './api';

interface Job {
  id: string;
  type: 'rfd3' | 'rf3' | 'mpnn';
  status: 'pending' | 'running' | 'completed' | 'failed';
  createdAt: string;
  completedAt?: string;
  result?: JobStatus['result'];
  error?: string;
}

// Notification types for workflow guidance
export interface Notification {
  id: string;
  type: 'success' | 'error' | 'info';
  title: string;
  message: string;
  action?: {
    label: string;
    tab: 'rfd3' | 'rf3' | 'mpnn' | 'jobs';
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
  activeTab: 'rfd3' | 'rf3' | 'mpnn' | 'jobs';
  setActiveTab: (tab: 'rfd3' | 'rf3' | 'mpnn' | 'jobs') => void;

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

  // Active tab
  activeTab: 'rfd3',
  setActiveTab: (tab) => set({ activeTab: tab }),

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
}));
