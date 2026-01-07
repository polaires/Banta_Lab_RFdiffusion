/**
 * API client for Foundry backend on RunPod (v2.0)
 * Supports confidence metrics, RMSD validation, and cross-panel data flow
 */

const API_BASE_URL = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000';

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

  constructor(baseUrl: string = API_BASE_URL) {
    this.baseUrl = baseUrl;
  }

  setBaseUrl(url: string) {
    this.baseUrl = url;
  }

  getBaseUrl(): string {
    return this.baseUrl;
  }

  // Health check
  async checkHealth(): Promise<HealthResponse> {
    const response = await fetch(`${this.baseUrl}/health`);
    if (!response.ok) {
      throw new Error('Backend not available');
    }
    return response.json();
  }

  // RFD3 Design
  async submitRFD3Design(request: RFD3Request): Promise<JobResponse> {
    const response = await fetch(`${this.baseUrl}/api/rfd3/design`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(request),
    });
    if (!response.ok) {
      throw new Error(`Failed to submit RFD3 job: ${response.statusText}`);
    }
    return response.json();
  }

  // RF3 Prediction
  async submitRF3Prediction(request: RF3Request): Promise<JobResponse> {
    const response = await fetch(`${this.baseUrl}/api/rf3/predict`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(request),
    });
    if (!response.ok) {
      throw new Error(`Failed to submit RF3 job: ${response.statusText}`);
    }
    return response.json();
  }

  // MPNN Design
  async submitMPNNDesign(request: ProteinMPNNRequest): Promise<JobResponse> {
    const response = await fetch(`${this.baseUrl}/api/mpnn/design`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(request),
    });
    if (!response.ok) {
      throw new Error(`Failed to submit MPNN job: ${response.statusText}`);
    }
    return response.json();
  }

  // RMSD Validation
  async calculateRMSD(request: RMSDRequest): Promise<RMSDResult> {
    const response = await fetch(`${this.baseUrl}/api/validate/rmsd`, {
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

  // Structure Export
  async exportStructure(request: ExportRequest): Promise<ExportResult> {
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
    const response = await fetch(`${this.baseUrl}/api/stored/${jobId}`);
    if (!response.ok) {
      throw new Error('Structure not found');
    }
    return response.json();
  }

  // Job Management
  async getJobStatus(jobId: string): Promise<JobStatus> {
    const response = await fetch(`${this.baseUrl}/api/jobs/${jobId}`);
    if (!response.ok) {
      throw new Error(`Failed to get job status: ${response.statusText}`);
    }
    return response.json();
  }

  async listJobs(): Promise<Record<string, { status: string; type: string; created_at: string }>> {
    const response = await fetch(`${this.baseUrl}/api/jobs`);
    if (!response.ok) {
      throw new Error(`Failed to list jobs: ${response.statusText}`);
    }
    return response.json();
  }

  async deleteJob(jobId: string): Promise<void> {
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
