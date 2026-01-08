/**
 * API client for Foundry backend on RunPod (v2.0)
 * Supports confidence metrics, RMSD validation, and cross-panel data flow
 *
 * Modes:
 * - "traditional": Direct connection to RunPod Pod backend
 * - "serverless": Uses Vercel Edge Function proxy to RunPod Serverless
 */

import {
  saveJob as supabaseSaveJob,
  updateJob as supabaseUpdateJob,
  getJobs as supabaseGetJobs,
  isSupabaseConfigured,
} from './supabase';

const API_BASE_URL = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000';

// Check if serverless mode is available
const isServerlessAvailable = () => {
  return !!(
    process.env.NEXT_PUBLIC_RUNPOD_SERVERLESS === 'true' ||
    (typeof window !== 'undefined' && (window as any).__RUNPOD_SERVERLESS__)
  );
};

export type ApiMode = 'traditional' | 'serverless';

// ============== Request Types ==============

export interface RFD3Request {
  contig: string;
  num_designs?: number;
  pdb_content?: string;
  seed?: number;  // For reproducibility
}

export interface RF3Request {
  sequence: string;
  name?: string;
  pdb_content?: string;  // For structure-based input
  config?: Record<string, unknown>;
}

export interface ProteinMPNNRequest {
  pdb_content: string;
  num_sequences?: number;
  temperature?: number;
  model_type?: 'ligand_mpnn' | 'protein_mpnn';
  remove_waters?: boolean;
  config?: Record<string, unknown>;
}

export interface RMSDRequest {
  pdb_content_1: string;  // Reference (e.g., RFD3 output)
  pdb_content_2: string;  // Comparison (e.g., RF3 refolded)
  backbone_only?: boolean;
}

export interface ExportRequest {
  pdb_content: string;
  format: 'pdb' | 'cif';
  include_confidences?: boolean;
  confidences?: ConfidenceMetrics;
}

// ============== Response Types ==============

export interface JobResponse {
  job_id: string;
  status: string;
  message: string;
}

export interface SummaryConfidences {
  overall_plddt: number;
  overall_pae: number;
  overall_pde?: number | null;
  ptm: number;
  iptm?: number | null;
  ranking_score: number;
  has_clash: boolean;
}

export interface ConfidenceMetrics {
  summary_confidences?: SummaryConfidences;
  per_residue_plddt?: number[];
  pae_matrix?: number[][];
}

export interface DesignOutput {
  filename: string;
  content: string;
  cif_content?: string;  // CIF format output
}

export interface JobStatus {
  job_id: string;
  status: 'pending' | 'running' | 'completed' | 'failed';
  created_at: string;
  completed_at?: string;
  result?: {
    designs?: DesignOutput[];
    predictions?: DesignOutput[];
    sequences?: Array<{ filename: string; content: string }>;
    confidences?: ConfidenceMetrics;
    model_type?: string;
    seed?: number;
    mode?: string;
    stdout?: string;
  };
  error?: string;
}

export interface ModelStatus {
  available: boolean;
  checkpoint_exists: boolean;
  checkpoint_size_gb: number;
}

export interface HealthResponse {
  status: string;
  mode: 'real' | 'mock';
  gpu_available: boolean;
  gpu_name?: string;
  gpu_memory_gb?: number;
  models: Record<string, ModelStatus>;
  // Legacy fields for backward compatibility
  models_loaded?: string[];
}

export interface RMSDResult {
  rmsd: number;
  interpretation: 'Excellent' | 'Good' | 'Moderate' | 'Poor';
  description: string;
  backbone_only: boolean;
  num_atoms_compared: number;
}

export interface ExportResult {
  filename: string;
  content: string;
  format: 'pdb' | 'cif';
  confidences_json?: string;
}

export interface StoredStructure {
  pdb: string;
  cif?: string;
  type: 'rfd3' | 'rf3' | 'mpnn';
  confidences?: ConfidenceMetrics;
}

// ============== API Client ==============

class FoundryAPI {
  private baseUrl: string;
  private mode: ApiMode;

  constructor(baseUrl: string = API_BASE_URL) {
    this.baseUrl = baseUrl;
    // Auto-detect mode: serverless if RUNPOD_SERVERLESS env is set
    this.mode = isServerlessAvailable() ? 'serverless' : 'traditional';
  }

  setBaseUrl(url: string) {
    this.baseUrl = url;
  }

  getBaseUrl(): string {
    return this.baseUrl;
  }

  setMode(mode: ApiMode) {
    this.mode = mode;
  }

  getMode(): ApiMode {
    return this.mode;
  }

  // Get the appropriate endpoint URL based on mode
  private getEndpoint(path: string): string {
    if (this.mode === 'serverless') {
      // Use Vercel Edge Function proxy
      return `/api/runpod/${path.replace(/^\/?(api\/)?/, '')}`;
    }
    // Traditional mode: direct to backend
    return `${this.baseUrl}/${path.replace(/^\//, '')}`;
  }

  // Health check
  async checkHealth(): Promise<HealthResponse> {
    const endpoint = this.mode === 'serverless'
      ? '/api/runpod/health'
      : `${this.baseUrl}/health`;

    const response = await fetch(endpoint);
    if (!response.ok) {
      throw new Error('Backend not available');
    }
    return response.json();
  }

  // RFD3 Design
  async submitRFD3Design(request: RFD3Request): Promise<JobResponse> {
    const endpoint = this.getEndpoint('api/rfd3/design');
    const response = await fetch(endpoint, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(request),
    });
    if (!response.ok) {
      throw new Error(`Failed to submit RFD3 job: ${response.statusText}`);
    }
    const result = await response.json();

    // Save to Supabase in serverless mode
    if (this.mode === 'serverless' && isSupabaseConfigured()) {
      await supabaseSaveJob({
        runpod_id: result.job_id,
        type: 'rfd3',
        request: { contig: request.contig, num_designs: request.num_designs, seed: request.seed },
      });
    }

    return result;
  }

  // RF3 Prediction
  async submitRF3Prediction(request: RF3Request): Promise<JobResponse> {
    const endpoint = this.getEndpoint('api/rf3/predict');
    const response = await fetch(endpoint, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(request),
    });
    if (!response.ok) {
      throw new Error(`Failed to submit RF3 job: ${response.statusText}`);
    }
    const result = await response.json();

    // Save to Supabase in serverless mode
    if (this.mode === 'serverless' && isSupabaseConfigured()) {
      await supabaseSaveJob({
        runpod_id: result.job_id,
        type: 'rf3',
        request: { sequence_length: request.sequence?.length },
      });
    }

    return result;
  }

  // MPNN Design
  async submitMPNNDesign(request: ProteinMPNNRequest): Promise<JobResponse> {
    const endpoint = this.getEndpoint('api/mpnn/design');
    const response = await fetch(endpoint, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(request),
    });
    if (!response.ok) {
      throw new Error(`Failed to submit MPNN job: ${response.statusText}`);
    }
    const result = await response.json();

    // Save to Supabase in serverless mode
    if (this.mode === 'serverless' && isSupabaseConfigured()) {
      await supabaseSaveJob({
        runpod_id: result.job_id,
        type: 'mpnn',
        request: {
          num_sequences: request.num_sequences,
          temperature: request.temperature,
          model_type: request.model_type,
        },
      });
    }

    return result;
  }

  // RMSD Validation
  async calculateRMSD(request: RMSDRequest): Promise<RMSDResult> {
    const endpoint = this.getEndpoint('api/validate/rmsd');
    const response = await fetch(endpoint, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(request),
    });
    if (!response.ok) {
      const error = await response.json().catch(() => ({ detail: response.statusText }));
      throw new Error(error.detail || 'Failed to calculate RMSD');
    }
    return response.json();
  }

  // Structure Export (not available in serverless mode)
  async exportStructure(request: ExportRequest): Promise<ExportResult> {
    if (this.mode === 'serverless') {
      // In serverless mode, export locally
      return {
        filename: `structure.${request.format}`,
        content: request.pdb_content,
        format: request.format,
        confidences_json: request.confidences ? JSON.stringify(request.confidences, null, 2) : undefined,
      };
    }

    const response = await fetch(`${this.baseUrl}/api/export/structure`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(request),
    });
    if (!response.ok) {
      const error = await response.json().catch(() => ({ detail: response.statusText }));
      throw new Error(error.detail || 'Failed to export structure');
    }
    return response.json();
  }

  // Get stored structure from job (for cross-panel data flow)
  async getStoredStructure(jobId: string): Promise<StoredStructure> {
    if (this.mode === 'serverless') {
      // In serverless mode, structures are stored in frontend state
      throw new Error('Stored structures not available in serverless mode');
    }

    const response = await fetch(`${this.baseUrl}/api/stored/${jobId}`);
    if (!response.ok) {
      throw new Error('Structure not found');
    }
    return response.json();
  }

  // Job Management
  async getJobStatus(jobId: string): Promise<JobStatus> {
    const endpoint = this.mode === 'serverless'
      ? `/api/runpod/jobs/${jobId}`
      : `${this.baseUrl}/api/jobs/${jobId}`;

    const response = await fetch(endpoint);
    if (!response.ok) {
      throw new Error(`Failed to get job status: ${response.statusText}`);
    }
    const status = await response.json();

    // Update Supabase in serverless mode
    if (this.mode === 'serverless' && isSupabaseConfigured()) {
      await supabaseUpdateJob(jobId, {
        status: status.status,
        result: status.result,
        completed_at: status.completed_at,
      });
    }

    return status;
  }

  async listJobs(): Promise<Record<string, { status: string; type: string; created_at: string }>> {
    if (this.mode === 'serverless' && isSupabaseConfigured()) {
      // In serverless mode, list from Supabase
      const jobs = await supabaseGetJobs({ limit: 50 });
      const result: Record<string, { status: string; type: string; created_at: string }> = {};
      for (const job of jobs) {
        result[job.runpod_id] = {
          status: job.status,
          type: job.type,
          created_at: job.created_at,
        };
      }
      return result;
    }

    const response = await fetch(`${this.baseUrl}/api/jobs`);
    if (!response.ok) {
      throw new Error(`Failed to list jobs: ${response.statusText}`);
    }
    return response.json();
  }

  async deleteJob(jobId: string): Promise<void> {
    if (this.mode === 'serverless') {
      // In serverless mode, jobs are deleted from Supabase
      // (RunPod auto-deletes after 1 hour)
      return;
    }

    const response = await fetch(`${this.baseUrl}/api/jobs/${jobId}`, {
      method: 'DELETE',
    });
    if (!response.ok) {
      throw new Error(`Failed to delete job: ${response.statusText}`);
    }
  }

  // Poll job status until completion
  async waitForJob(
    jobId: string,
    onStatusChange?: (status: JobStatus) => void,
    pollInterval: number = 2000
  ): Promise<JobStatus> {
    while (true) {
      const status = await this.getJobStatus(jobId);
      if (onStatusChange) {
        onStatusChange(status);
      }

      if (status.status === 'completed' || status.status === 'failed') {
        return status;
      }

      await new Promise(resolve => setTimeout(resolve, pollInterval));
    }
  }
}

export const api = new FoundryAPI();
export default api;
