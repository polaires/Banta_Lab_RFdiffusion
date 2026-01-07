/**
 * Zustand store for application state
 */

import { create } from 'zustand';
import type { JobStatus, HealthResponse } from './api';

interface Job {
  id: string;
  type: 'rfd3' | 'rf3' | 'mpnn';
  status: 'pending' | 'running' | 'completed' | 'failed';
  createdAt: string;
  completedAt?: string;
  result?: JobStatus['result'];
  error?: string;
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
}

export const useStore = create<AppState>((set) => ({
  // Backend connection
  backendUrl: process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000',
  setBackendUrl: (url) => set({ backendUrl: url }),

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
}));
