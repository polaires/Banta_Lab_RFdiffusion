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

// Phase 2: Symmetry configuration for oligomeric design
export interface SymmetryConfig {
  id: string;  // C2, C3, C4, C5, C6, D2, D3, D4, T, O, I
  is_unsym_motif?: string;  // Chains to exclude from symmetry
  is_symmetric_motif?: boolean;  // Default true
}

// Phase 2: Classifier-Free Guidance configuration
export interface CFGConfig {
  enabled: boolean;
  scale?: number;  // 0.5-3.0, default 1.5
  t_max?: number;  // 0.0-1.0, default 0.8
  features?: string[];  // active_donor, active_acceptor, ref_atomwise_rasa
}

export interface RFD3Request {
  // Core parameters - either contig OR length is required
  contig?: string;
  length?: string;  // Can be "100" or "80-120" range
  num_designs?: number;
  pdb_content?: string;
  seed?: number;  // For reproducibility

  // Quality/sampling parameters
  num_timesteps?: number;  // 50-500, default 200
  step_scale?: number;     // 0.5-3.0, default 1.5 (higher = less diverse, more designable)
  noise_scale?: number;    // 0.9-1.1, default 1.003
  gamma_0?: number;        // 0.1-1.0, default 0.6 (lower = more designable, less diverse)
  is_non_loopy?: boolean;  // Cleaner secondary structures

  // Partial diffusion (refinement)
  partial_t?: number;  // Noise level in Angstroms (5-20 recommended)

  // Symmetry design - for oligomeric proteins
  symmetry?: SymmetryConfig;

  // Protein binder design
  hotspots?: string[];  // Hotspot residues (e.g., ["A15", "A20", "A25"])
  infer_ori_strategy?: string;  // "hotspots" recommended for binders

  // Small molecule / enzyme design
  ligand?: string;  // Chemical component ID (e.g., "ATP", "NAD", "ZN")
  unindex?: string;  // Catalytic residues with inferred positions
  select_fixed_atoms?: Record<string, string>;  // Fixed atoms on residues/ligands

  // RASA conditioning (binding pocket design)
  select_buried?: Record<string, string>;  // Atoms to bury in binding pocket
  select_exposed?: Record<string, string>;  // Atoms to keep solvent accessible
  select_partially_buried?: Record<string, string>;  // Atoms at surface edge

  // Nucleic acid binder design
  na_chains?: string;  // NA chain IDs (e.g., "A,B")
  ori_token?: number[];  // Origin position [x, y, z]
  select_hbond_donor?: Record<string, string>;  // H-bond donor conditioning
  select_hbond_acceptor?: Record<string, string>;  // H-bond acceptor conditioning

  // Interface ligand design (separable dimer)
  task?: string;  // Task type for specialized handlers (e.g., 'interface_ligand_design', 'protein_binder_design')
  approach?: 'asymmetric' | 'sequential' | 'full';  // Design approach
  ligand_smiles?: string;  // SMILES string for ligand
  chain_length?: string;  // Chain length range (e.g., "60-80")
  side?: 'left' | 'right';  // Which side of ligand to bind
  ori_offset?: number[];  // Offset from ligand center [x, y, z]
  use_ori_token?: boolean;  // Whether to use ori_token positioning

  // Protein binder design (BindCraft-inspired pipeline)
  target_pdb?: string;  // PDB content for target protein
  binder_length?: string;  // Binder length range (e.g., "60-80")
  quality_threshold?: 'relaxed' | 'standard' | 'strict';  // Quality filter threshold
  protocol?: string;  // Protocol preset (e.g., 'miniprotein_default', 'peptide_default')
  validate_structure?: boolean;  // Enable ESMFold structure validation
  auto_hotspots?: boolean;  // Auto-detect hotspots using SASA analysis
  filter_wrap_around?: boolean;  // Filter out wrap-around binders using Rg ratio

  // Legacy/deprecated
  select_hotspots?: Record<string, string>;  // Old hotspot format
  cfg?: CFGConfig;  // Classifier-Free Guidance
}

export interface RF3Request {
  sequence: string;
  name?: string;
  pdb_content?: string;  // For structure-based input
  msa_content?: string;  // MSA file content (.a3m or .fasta format)
  config?: Record<string, unknown>;
}

export interface ProteinMPNNRequest {
  pdb_content: string;
  num_sequences?: number;
  temperature?: number;
  model_type?: 'ligand_mpnn' | 'protein_mpnn';
  remove_waters?: boolean;
  fixed_positions?: string[];  // e.g., ["A35", "A36", "B35"]

  // Advanced LigandMPNN parameters (Nature Methods 2025)
  pack_side_chains?: boolean;  // Enable sidechain packing
  pack_with_ligand_context?: boolean;  // Include ligand when packing
  number_of_packs_per_design?: number;  // Packing samples per sequence
  bias_AA?: string;  // e.g., "W:3.0,Y:2.0,C:-5.0"
  omit_AA?: string;  // e.g., "C" to omit cysteine
  model_noise_level?: '005' | '010' | '020' | '030';  // Model noise level
  ligand_cutoff_for_score?: number;  // Angstroms (default 8.0)
  use_side_chain_context?: boolean;  // Use fixed sidechains as context
  save_stats?: boolean;  // Return confidence metrics

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

export interface ErrorContext {
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
  error_type?: string;
  traceback?: string;
  context?: ErrorContext;
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

// Structure Analysis Types
export interface LigandInfo {
  name: string;
  chain: string;
  res_id: number;
  num_atoms: number;
  center: number[];
  atom_names: string[];
}

export interface CoordinatingResidue {
  chain: string;
  res_id: number;
  res_name: string;
  atoms: Array<{ name: string; distance: number }>;
  min_distance: number;
}

export interface BindingSite {
  ligand: LigandInfo;
  coordinating_residues: CoordinatingResidue[];
  likely_coordinators: string[];
  coordination_number: number;
}

export interface DesignSuggestion {
  type: string;
  description: string;
  rfd3_params: {
    ligand?: string;
    partial_t?: number;
    unindex?: string;
    is_non_loopy?: boolean;
    num_timesteps?: number;
  };
  notes: string[];
}

export interface AnalysisResult {
  num_residues: number;
  ligands: LigandInfo[];
  binding_sites: BindingSite[];
  suggestions: DesignSuggestion[];
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

// ============== AI Assistant Types ==============

export interface FetchPdbRequest {
  pdb_id: string;
}

export interface FetchPdbResponse {
  success: boolean;
  content: string;
  source: string;
  pdb_id: string;
  info?: {
    title: string;
    chains: string[];
    num_atoms: number;
    num_residues: number;
    metals: Array<{
      element: string;
      chain: string;
      residue: string;
      res_num: string;
    }>;
  };
  error?: string;
}

export interface MetalBindingAnalysisRequest {
  pdb_content?: string;
  pdb_id?: string;
  metal_chain?: string;
  metal_residue: string;
  metal_resnum: number;
  distance_cutoff?: number;
}

export interface CoordinatingAtom {
  chain: string;
  residue: string;
  residue_number: number;
  atom: string;
  element: string;
  distance: number;
  donor_type: string;
}

export interface MetalBindingAnalysis {
  success: boolean;
  metal: {
    element: string;
    chain: string;
    residue: string;
    resnum: number;
    position: number[];
  };
  coordination: {
    number: number;
    geometry: string;
    geometry_rmsd: number;
    coordinating_atoms: CoordinatingAtom[];
  };
  donor_analysis: {
    types: Record<string, number>;
    dominant_type: string;
    lanthanide_compatible: boolean;
  };
  bond_analysis: {
    average_distance: number;
    min_distance: number;
    max_distance: number;
    distances: number[];
  };
  suggestions?: string[];
  error?: string;
}

export interface ParameterRecommendationRequest {
  pdb_content?: string;
  pdb_id?: string;
  metal_chain?: string;
  metal_residue?: string;
  metal_resnum?: number;
  target_metal: string;
  user_description?: string;
}

export interface ParameterRecommendation {
  success: boolean;
  strategy: string;
  reasoning: string[];
  coordination_analysis: {
    current_coordination: number;
    target_coordination_range: number[];
    delta: number;
    recommendation: string;
  };
  parameters: {
    ligand: string;
    partial_t: number;
    unindex?: string;
    select_fixed_atoms?: Record<string, string>;
    num_designs: number;
    num_timesteps: number;
    step_scale: number;
    gamma_0: number;
  };
  evaluation_criteria: {
    target_coordination: number[];
    target_distance_range: number[];
    preferred_geometry: string;
    preferred_donors: string[];
  };
  error?: string;
}

export interface EvaluateDesignRequest {
  pdb_content: string;
  target_metal: string;
  metal_chain?: string;
  metal_resnum?: number;
  expected_coordination?: number[];
  expected_donors?: string[];
}

// Hotspot detection for protein binder design
export interface DetectHotspotsRequest {
  target_pdb: string;
  target_chain?: string;
  method?: 'exposed' | 'exposed_clustered' | 'patch' | 'interface_like';
  max_hotspots?: number;
  prefer_hydrophobic?: boolean;
  min_relative_sasa?: number;
}

export interface HotspotResidue {
  residue: string;  // e.g., "A25"
  restype: string;  // e.g., "LEU"
  relative_sasa: number;
  property: 'hydrophobic' | 'polar' | 'charged' | 'special';
  coords?: [number, number, number];
}

export interface DetectHotspotsResponse {
  status: 'completed' | 'failed';
  hotspots: string[];  // ["A25", "A30", ...]
  method: string;
  residue_details: HotspotResidue[];
  total_exposed: number;
  cluster_center?: { x: number; y: number; z: number };
  patch_area?: number;
  error?: string;
}

export interface DesignEvaluation {
  success: boolean;
  coordination_number: number;
  avg_bond_distance: number;
  geometry_type: string;
  geometry_rmsd: number;
  oxygen_donors: number;
  nitrogen_donors: number;
  sulfur_donors: number;
  criteria_passed: number;
  criteria_total: number;
  overall_pass: boolean;
  suggestions?: string[];
  error?: string;
}

// User preferences from interview mode
export interface UserPreferences {
  targetMetal: string;
  targetMetalLabel: string;
  designAggressiveness: 'conservative' | 'moderate' | 'aggressive';
  aggressivenessLabel: string;
  coordinationPreference: 'tetrahedral' | 'octahedral' | 'high' | 'auto';
  coordinationLabel: string;
  numDesigns: number;
  priority: 'stability' | 'binding' | 'expression' | 'balanced';
  priorityLabel: string;
}

// ============== ESM3 Types ==============

export interface ESM3ScoreRequest {
  sequences: string[];
}

export interface ESM3ScoreResult {
  scores: Array<{
    sequence: string;
    length: number;
    perplexity: number;
    score: number;
    error?: string;
  }>;
}

export interface ESM3GenerateRequest {
  prompt?: string;
  functions?: string[];  // e.g., ["zinc-binding", "hydrolase"]
  num_sequences?: number;
  temperature?: number;
  max_length?: number;
}

export interface ESM3GenerateResult {
  sequences: string[];
  num_generated: number;
  temperature: number;
  max_length: number;
  function_keywords?: string[];
}

export interface ESM3EmbedRequest {
  sequences: string[];
}

export interface ESM3EmbedResult {
  embeddings: Array<{
    sequence: string;
    length: number;
    embedding_dim: number;
    per_residue: number[][];
    global: number[];
    error?: string;
  }>;
  num_processed: number;
}

// ============== API Client ==============

class FoundryAPI {
  private baseUrl: string;
  private mode: ApiMode;
  // Cache for sync job results in traditional mode (jobs complete immediately)
  private syncJobResults: Map<string, any> = new Map();

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
    // Traditional mode: use Next.js proxy to avoid CORS issues
    const cleanPath = path.replace(/^\/?(api\/)?/, '');
    return `/api/traditional/${cleanPath}?url=${encodeURIComponent(this.baseUrl)}`;
  }

  // Health check
  async checkHealth(): Promise<HealthResponse> {
    if (this.mode === 'serverless') {
      // Serverless mode uses the Next.js API route
      const response = await fetch('/api/runpod/health');
      if (!response.ok) {
        throw new Error('Backend not available');
      }
      return response.json();
    } else {
      // Traditional mode uses Next.js proxy to avoid CORS issues
      const response = await fetch(`/api/traditional/health?url=${encodeURIComponent(this.baseUrl)}`);
      if (!response.ok) {
        throw new Error('Backend not available');
      }
      return response.json();
    }
  }

  // RFD3 Design
  async submitRFD3Design(request: RFD3Request): Promise<JobResponse & { result?: any; error?: string; syncCompleted?: boolean }> {
    if (this.mode === 'serverless') {
      // Serverless mode uses the Next.js API route
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
      if (isSupabaseConfigured()) {
        await supabaseSaveJob({
          runpod_id: result.job_id,
          type: 'rfd3',
          request: {
            contig: request.contig,
            length: request.length,
            num_designs: request.num_designs,
            seed: request.seed,
            num_timesteps: request.num_timesteps,
            step_scale: request.step_scale,
            gamma_0: request.gamma_0,
            is_non_loopy: request.is_non_loopy,
            partial_t: request.partial_t,
            symmetry: request.symmetry,
            hotspots: request.hotspots,
            ligand: request.ligand,
            unindex: request.unindex,
          },
        });
      }

      return result;
    } else {
      // Traditional mode uses Next.js proxy to avoid CORS issues
      // runsync is SYNCHRONOUS - result comes back in same request
      const proxyUrl = `/api/traditional/runsync?url=${encodeURIComponent(this.baseUrl)}`;
      console.log('[API] Submitting RFD3 design via traditional mode:', proxyUrl);

      let response: Response;
      try {
        response = await fetch(proxyUrl, {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ input: { task: 'rfd3', ...request } }),
        });
      } catch (fetchError) {
        console.error('[API] Fetch error:', fetchError);
        throw new Error(`Network error connecting to backend. Check if backend is running at ${this.baseUrl}`);
      }

      if (!response.ok) {
        const errorText = await response.text().catch(() => 'No error details');
        console.error('[API] Backend error:', response.status, errorText);
        throw new Error(`Backend error (${response.status}): ${errorText}`);
      }

      const result = await response.json();
      console.log('[API] RFD3 result:', result.status, result.id);

      // RunPod runsync returns result directly in output
      // Format: { id, status: "COMPLETED"|"FAILED", output: { status: "completed"|"failed"|"error", ... } }
      // The outer status is the RunPod wrapper status
      // The inner status (output.status) is the handler's actual status
      const output = result.output || {};
      const innerStatus = output.status;

      // Check both outer RunPod status and inner handler status
      // Handler returns "completed", "failed", or "error"
      const isSuccess = result.status === 'COMPLETED' &&
                       innerStatus === 'completed' &&
                       !output.error;

      console.log('[API] Inner status:', innerStatus, 'isSuccess:', isSuccess);

      return {
        job_id: result.id || 'sync-' + Date.now(),
        status: isSuccess ? 'completed' : 'failed',
        message: isSuccess ? 'Design completed' : (output.error || 'Design failed'),
        // Return the inner result, not the full output wrapper
        // output = { status: "completed", result: { designs: [...] } }
        // We want just output.result = { designs: [...] }
        result: isSuccess ? output.result : undefined,
        error: output.error || (innerStatus === 'error' || innerStatus === 'failed' ? 'Handler returned error status' : undefined),
        syncCompleted: true,  // Signal that result is already available
      };
    }
  }

  // RF3 Prediction
  async submitRF3Prediction(request: RF3Request): Promise<JobResponse> {
    if (this.mode === 'serverless') {
      // Serverless mode uses the Next.js API route
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
      if (isSupabaseConfigured()) {
        await supabaseSaveJob({
          runpod_id: result.job_id,
          type: 'rf3',
          request: { sequence_length: request.sequence?.length },
        });
      }

      return result;
    }

    // Traditional mode uses Next.js proxy to avoid CORS issues
    const proxyUrl = `/api/traditional/runsync?url=${encodeURIComponent(this.baseUrl)}`;
    const response = await fetch(proxyUrl, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ input: { task: 'rf3', ...request } }),
    });

    if (!response.ok) {
      throw new Error(`Failed to submit RF3 job: ${response.statusText}`);
    }

    const result = await response.json();

    // Parse RunPod response format
    const output = result.output || {};
    const jobId = result.id || 'sync-' + Date.now();
    const jobResult = output.result;

    // Cache the result for getJobStatus() to return
    this.syncJobResults.set(jobId, jobResult);

    return {
      job_id: jobId,
      status: output.status === 'completed' ? 'completed' : 'failed',
      message: output.status === 'completed' ? 'RF3 prediction completed' : (output.error || 'RF3 prediction failed'),
      result: jobResult,
      error: output.error,
      syncCompleted: true,
    };
  }

  // MPNN Design
  async submitMPNNDesign(request: ProteinMPNNRequest): Promise<JobResponse> {
    if (this.mode === 'serverless') {
      // Serverless mode uses the Next.js API route
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
      if (isSupabaseConfigured()) {
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

    // Traditional mode uses Next.js proxy to avoid CORS issues
    const proxyUrl = `/api/traditional/runsync?url=${encodeURIComponent(this.baseUrl)}`;
    const response = await fetch(proxyUrl, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ input: { task: 'mpnn', ...request } }),
    });

    if (!response.ok) {
      throw new Error(`Failed to submit MPNN job: ${response.statusText}`);
    }

    const result = await response.json();

    // Parse RunPod response format
    const output = result.output || {};
    const jobId = result.id || 'sync-' + Date.now();
    const jobResult = output.result;

    // Cache the result for getJobStatus() to return
    this.syncJobResults.set(jobId, jobResult);

    return {
      job_id: jobId,
      status: output.status === 'completed' ? 'completed' : 'failed',
      message: output.status === 'completed' ? 'MPNN design completed' : (output.error || 'MPNN design failed'),
      result: jobResult,
      error: output.error,
      syncCompleted: true,
    };
  }

  // RMSD Validation
  async calculateRMSD(request: RMSDRequest): Promise<RMSDResult> {
    if (this.mode === 'serverless') {
      // Serverless mode uses the Next.js API route
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

    // Traditional mode uses Next.js proxy to avoid CORS issues
    const proxyUrl = `/api/traditional/runsync?url=${encodeURIComponent(this.baseUrl)}`;
    const response = await fetch(proxyUrl, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ input: { task: 'rmsd', ...request } }),
    });

    if (!response.ok) {
      const error = await response.json().catch(() => ({ detail: response.statusText }));
      throw new Error(error.detail || 'Failed to calculate RMSD');
    }

    const result = await response.json();

    // Parse RunPod response format
    const output = result.output || {};
    if (output.status === 'completed' && output.result) {
      return output.result;
    } else {
      throw new Error(output.error || 'RMSD calculation failed');
    }
  }

  // Structure Analysis - AI-assisted binding site analysis
  async analyzeStructure(pdbContent: string, targetLigands?: string[]): Promise<AnalysisResult> {
    if (this.mode === 'serverless') {
      // Serverless mode - not yet implemented
      throw new Error('Structure analysis not yet available in serverless mode');
    }

    // Traditional mode uses Next.js proxy
    const response = await fetch(`/api/traditional/runsync?url=${encodeURIComponent(this.baseUrl)}`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        input: {
          task: 'analyze',
          pdb_content: pdbContent,
          target_ligands: targetLigands,
        }
      }),
    });

    if (!response.ok) {
      throw new Error(`Failed to analyze structure: ${response.statusText}`);
    }

    const result = await response.json();

    if (result.output?.status === 'failed') {
      throw new Error(result.output.error || 'Analysis failed');
    }

    return result.output?.result || result.output;
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

  // ============== AI Assistant Methods ==============

  // Fetch PDB from RCSB via backend or directly from RCSB
  async fetchPdb(pdbId: string): Promise<FetchPdbResponse> {
    // Try to fetch directly from RCSB (works in browser)
    const pdbIdUpper = pdbId.toUpperCase();
    const rcsbUrl = `https://files.rcsb.org/download/${pdbIdUpper}.pdb`;

    try {
      const response = await fetch(rcsbUrl);
      if (!response.ok) {
        throw new Error(`PDB ${pdbIdUpper} not found`);
      }

      const content = await response.text();

      // Parse basic info from PDB content
      const info = this.parsePdbInfo(content);

      return {
        success: true,
        content,
        source: 'RCSB PDB',
        pdb_id: pdbIdUpper,
        info,
      };
    } catch (error) {
      throw new Error(`Failed to fetch PDB ${pdbIdUpper}: ${error instanceof Error ? error.message : 'Unknown error'}`);
    }
  }

  // Parse basic info from PDB content
  private parsePdbInfo(pdbContent: string): FetchPdbResponse['info'] {
    const lines = pdbContent.split('\n');
    let title = '';
    const chains = new Set<string>();
    const metals: Array<{ element: string; chain: string; residue: string; res_num: string }> = [];
    let numAtoms = 0;
    const residues = new Set<string>();

    // Common metal elements
    const metalElements = ['FE', 'ZN', 'CA', 'MG', 'MN', 'CU', 'CO', 'NI', 'TB', 'GD', 'EU', 'LA', 'CE', 'ND'];

    for (const line of lines) {
      // Get title
      if (line.startsWith('TITLE')) {
        title += line.slice(10).trim() + ' ';
      }

      // Get atoms
      if (line.startsWith('ATOM') || line.startsWith('HETATM')) {
        numAtoms++;
        const chain = line.slice(21, 22).trim() || 'A';
        chains.add(chain);

        const resName = line.slice(17, 20).trim();
        const resNum = line.slice(22, 26).trim();
        residues.add(`${chain}:${resNum}`);

        // Check for metals
        if (line.startsWith('HETATM')) {
          const element = line.slice(76, 78).trim() || resName.slice(0, 2);
          if (metalElements.includes(element.toUpperCase()) || metalElements.includes(resName.toUpperCase())) {
            // Avoid duplicates
            const metalKey = `${chain}:${resName}:${resNum}`;
            if (!metals.find(m => `${m.chain}:${m.residue}:${m.res_num}` === metalKey)) {
              metals.push({
                element: element.toUpperCase() || resName.toUpperCase(),
                chain,
                residue: resName,
                res_num: resNum,
              });
            }
          }
        }
      }
    }

    return {
      title: title.trim(),
      chains: Array.from(chains),
      num_atoms: numAtoms,
      num_residues: residues.size,
      metals,
    };
  }

  // Analyze metal binding site
  async analyzeMetalBinding(request: MetalBindingAnalysisRequest): Promise<MetalBindingAnalysis> {
    if (this.mode === 'serverless') {
      throw new Error('Metal binding analysis not available in serverless mode');
    }

    const response = await fetch(`/api/traditional/runsync?url=${encodeURIComponent(this.baseUrl)}`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        input: {
          task: 'analyze_metal_binding',
          ...request,
        }
      }),
    });

    if (!response.ok) {
      throw new Error(`Failed to analyze metal binding: ${response.statusText}`);
    }

    const result = await response.json();
    if (result.output?.status === 'failed') {
      throw new Error(result.output.error || 'Metal binding analysis failed');
    }

    return result.output?.result || result.output;
  }

  // Get AI parameter recommendation
  async getParameterRecommendation(request: ParameterRecommendationRequest): Promise<ParameterRecommendation> {
    if (this.mode === 'serverless') {
      throw new Error('Parameter recommendation not available in serverless mode');
    }

    const response = await fetch(`/api/traditional/runsync?url=${encodeURIComponent(this.baseUrl)}`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        input: {
          task: 'recommend_parameters',
          ...request,
        }
      }),
    });

    if (!response.ok) {
      throw new Error(`Failed to get recommendations: ${response.statusText}`);
    }

    const result = await response.json();
    if (result.output?.status === 'failed') {
      throw new Error(result.output.error || 'Parameter recommendation failed');
    }

    return result.output?.result || result.output;
  }

  // Evaluate design output
  async evaluateDesign(request: EvaluateDesignRequest): Promise<DesignEvaluation> {
    if (this.mode === 'serverless') {
      throw new Error('Design evaluation not available in serverless mode');
    }

    const response = await fetch(`/api/traditional/runsync?url=${encodeURIComponent(this.baseUrl)}`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        input: {
          task: 'evaluate_design',
          ...request,
        }
      }),
    });

    if (!response.ok) {
      throw new Error(`Failed to evaluate design: ${response.statusText}`);
    }

    const result = await response.json();
    if (result.output?.status === 'failed') {
      throw new Error(result.output.error || 'Design evaluation failed');
    }

    return result.output?.result || result.output;
  }

  // ============== Hotspot Detection ==============

  /**
   * Detect optimal binding hotspots for protein binder design.
   * Uses SASA analysis and spatial clustering to find binding patches.
   */
  async detectHotspots(request: DetectHotspotsRequest): Promise<DetectHotspotsResponse> {
    console.log('[API] Detecting hotspots for target chain:', request.target_chain || 'A');

    // Traditional mode uses Next.js proxy
    const response = await fetch(`/api/traditional/runsync?url=${encodeURIComponent(this.baseUrl)}`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        input: {
          task: 'detect_hotspots',
          ...request,
        }
      }),
    });

    if (!response.ok) {
      throw new Error(`Failed to detect hotspots: ${response.statusText}`);
    }

    const result = await response.json();
    console.log('[API] Hotspot detection result:', result.output?.status);

    if (result.output?.status === 'failed') {
      throw new Error(result.output.error || 'Hotspot detection failed');
    }

    return result.output || result;
  }

  // Job Management
  async getJobStatus(jobId: string): Promise<JobStatus> {
    if (this.mode === 'serverless') {
      // Serverless mode: poll RunPod API via proxy
      const endpoint = `/api/runpod/jobs/${jobId}`;
      const response = await fetch(endpoint);
      if (!response.ok) {
        throw new Error(`Failed to get job status: ${response.statusText}`);
      }
      const status = await response.json();

      // Update Supabase in serverless mode
      if (isSupabaseConfigured()) {
        await supabaseUpdateJob(jobId, {
          status: status.status,
          result: status.result,
          completed_at: status.completed_at,
        });
      }

      return status;
    }

    // Traditional mode: jobs are synchronous (runsync), no polling needed
    // Return cached result if available
    const cachedResult = this.syncJobResults.get(jobId);
    if (cachedResult !== undefined) {
      // Clean up cache after returning (one-time use)
      this.syncJobResults.delete(jobId);
      return {
        status: 'completed',
        result: cachedResult,
      } as JobStatus;
    }

    // For any job ID without cached result, return completed with null
    // (result was already consumed or job didn't cache)
    return {
      status: 'completed',
      result: null,
    } as JobStatus;
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

    // Traditional mode: no job tracking, return empty
    return {};
  }

  async deleteJob(jobId: string): Promise<void> {
    if (this.mode === 'serverless') {
      // In serverless mode, jobs are deleted from Supabase
      // (RunPod auto-deletes after 1 hour)
      return;
    }

    // Traditional mode: no job tracking, no-op
    return;
  }

  // Poll job status until completion
  async waitForJob(
    jobId: string,
    onStatusChange?: (status: JobStatus) => void,
    pollInterval: number = 2000
  ): Promise<JobStatus> {
    console.log(`[API] Starting to poll job ${jobId}...`);
    const startTime = Date.now();
    let pollCount = 0;

    while (true) {
      pollCount++;
      const status = await this.getJobStatus(jobId);

      // Log every 10th poll or status changes
      if (pollCount % 10 === 1) {
        const elapsed = Math.round((Date.now() - startTime) / 1000);
        console.log(`[API] Job ${jobId} status: ${status.status} (poll #${pollCount}, ${elapsed}s elapsed)`);
      }

      if (onStatusChange) {
        onStatusChange(status);
      }

      if (status.status === 'completed' || status.status === 'failed') {
        const elapsed = Math.round((Date.now() - startTime) / 1000);
        console.log(`[API] Job ${jobId} finished with status: ${status.status} after ${elapsed}s (${pollCount} polls)`);
        return status;
      }

      await new Promise(resolve => setTimeout(resolve, pollInterval));
    }
  }

  // ============== ESM3 Methods ==============

  /**
   * Score protein sequences using ESM3 perplexity.
   * Lower perplexity = higher quality/more natural sequence.
   */
  async esm3ScoreSequences(sequences: string[]): Promise<ESM3ScoreResult> {
    if (this.mode === 'serverless') {
      const response = await fetch('/api/runpod/esm3/score', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ sequences }),
      });
      if (!response.ok) {
        throw new Error(`ESM3 scoring failed: ${response.statusText}`);
      }
      const result = await response.json();
      if (result.status === 'failed') {
        throw new Error(result.error || 'ESM3 scoring failed');
      }
      return result.result;
    }

    // Traditional mode
    const response = await fetch(`/api/traditional/runsync?url=${encodeURIComponent(this.baseUrl)}`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        input: { task: 'esm3_score', sequences }
      }),
    });

    if (!response.ok) {
      throw new Error(`ESM3 scoring failed: ${response.statusText}`);
    }

    const result = await response.json();
    if (result.output?.status === 'failed') {
      throw new Error(result.output.error || 'ESM3 scoring failed');
    }
    return result.output?.result || result.output;
  }

  /**
   * Generate protein sequences using ESM3 with optional function conditioning.
   */
  async esm3GenerateSequence(params: ESM3GenerateRequest): Promise<ESM3GenerateResult> {
    if (this.mode === 'serverless') {
      const response = await fetch('/api/runpod/esm3/generate', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(params),
      });
      if (!response.ok) {
        throw new Error(`ESM3 generation failed: ${response.statusText}`);
      }
      const result = await response.json();
      if (result.status === 'failed') {
        throw new Error(result.error || 'ESM3 generation failed');
      }
      return result.result;
    }

    // Traditional mode
    const response = await fetch(`/api/traditional/runsync?url=${encodeURIComponent(this.baseUrl)}`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        input: { task: 'esm3_generate', ...params }
      }),
    });

    if (!response.ok) {
      throw new Error(`ESM3 generation failed: ${response.statusText}`);
    }

    const result = await response.json();
    if (result.output?.status === 'failed') {
      throw new Error(result.output.error || 'ESM3 generation failed');
    }
    return result.output?.result || result.output;
  }

  /**
   * Get ESM3 embeddings for protein sequences.
   * Useful for sequence similarity and clustering.
   */
  async esm3GetEmbeddings(sequences: string[]): Promise<ESM3EmbedResult> {
    if (this.mode === 'serverless') {
      const response = await fetch('/api/runpod/esm3/embed', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ sequences }),
      });
      if (!response.ok) {
        throw new Error(`ESM3 embedding failed: ${response.statusText}`);
      }
      const result = await response.json();
      if (result.status === 'failed') {
        throw new Error(result.error || 'ESM3 embedding failed');
      }
      return result.result;
    }

    // Traditional mode
    const response = await fetch(`/api/traditional/runsync?url=${encodeURIComponent(this.baseUrl)}`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        input: { task: 'esm3_embed', sequences }
      }),
    });

    if (!response.ok) {
      throw new Error(`ESM3 embedding failed: ${response.statusText}`);
    }

    const result = await response.json();
    if (result.output?.status === 'failed') {
      throw new Error(result.output.error || 'ESM3 embedding failed');
    }
    return result.output?.result || result.output;
  }

  // ============== AI-Driven Analysis ==============

  /**
   * Analyze user request with AI to determine task type and parameters.
   * Uses hybrid handbook + LLM approach.
   */
  async aiAnalyze(request: {
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
  }): Promise<{
    success: boolean;
    task_type: string;
    params: Record<string, unknown>;
    reasoning: string;
    confidence: number;
    clarifying_questions: string[];
    form_config: Array<{
      id: string;
      label: string;
      type: string;
      required: boolean;
      suggested_value?: unknown;
      options?: Array<{ value: string; label: string }>;
      range_config?: { min: number; max: number; step: number };
      help_text: string;
      ai_reasoning?: string;
    }>;
    error?: string;
  }> {
    console.log('[API] AI analyze request:', request.user_input);

    const response = await fetch(`${this.baseUrl}/api/ai/analyze`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(request),
    });

    if (!response.ok) {
      const error = await response.text();
      throw new Error(`AI analysis failed: ${error}`);
    }

    const result = await response.json();
    console.log('[API] AI analyze result:', result);
    return result;
  }

  /**
   * Refine AI-suggested parameters based on user feedback.
   */
  async aiRefine(request: {
    current_params: Record<string, unknown>;
    user_feedback: string;
    task_type: string;
  }): Promise<{
    success: boolean;
    params: Record<string, unknown>;
    reasoning: string;
    confidence: number;
  }> {
    const response = await fetch(`${this.baseUrl}/api/ai/refine`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(request),
    });

    if (!response.ok) {
      const error = await response.text();
      throw new Error(`AI refinement failed: ${error}`);
    }

    return response.json();
  }
}

export const api = new FoundryAPI();
export default api;
