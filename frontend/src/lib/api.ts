/**
 * API client for Foundry backend on RunPod
 */

const API_BASE_URL = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000';

export interface RFD3Request {
  contig: string;
  num_designs?: number;
  pdb_content?: string;
}

export interface RF3Request {
  sequence: string;
  config?: Record<string, unknown>;
}

export interface ProteinMPNNRequest {
  pdb_content: string;
  num_sequences?: number;
  temperature?: number;
  config?: Record<string, unknown>;
}

export interface JobResponse {
  job_id: string;
  status: string;
  message: string;
}

export interface JobStatus {
  job_id: string;
  status: 'pending' | 'running' | 'completed' | 'failed';
  created_at: string;
  completed_at?: string;
  result?: {
    designs?: Array<{ filename: string; content: string }>;
    predictions?: Array<{ filename: string; content: string }>;
    sequences?: Array<{ filename: string; content: string }>;
    stdout?: string;
  };
  error?: string;
}

export interface HealthResponse {
  status: string;
  gpu_available: boolean;
  models_loaded: string[];
}

class FoundryAPI {
  private baseUrl: string;

  constructor(baseUrl: string = API_BASE_URL) {
    this.baseUrl = baseUrl;
  }

  setBaseUrl(url: string) {
    this.baseUrl = url;
  }

  async checkHealth(): Promise<HealthResponse> {
    const response = await fetch(`${this.baseUrl}/health`);
    if (!response.ok) {
      throw new Error('Backend not available');
    }
    return response.json();
  }

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
