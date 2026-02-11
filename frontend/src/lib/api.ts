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
import { getCurrentUser } from './auth';

/**
 * Get the current authenticated user's ID, or null if not logged in.
 */
async function getCurrentUserId(): Promise<string | null> {
  try {
    const user = await getCurrentUser();
    return user?.id ?? null;
  } catch {
    return null;
  }
}

const API_BASE_URL = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000';

// Detect whether to use serverless mode.
// Vercel deployments always use RunPod serverless; local dev uses traditional (localhost backend).
const isServerlessAvailable = () => {
  // Explicit env override takes priority
  if (process.env.NEXT_PUBLIC_RUNPOD_SERVERLESS === 'true') return true;
  if (process.env.NEXT_PUBLIC_RUNPOD_SERVERLESS === 'false') return false;
  // Runtime override
  if (typeof window !== 'undefined' && (window as any).__RUNPOD_SERVERLESS__) return true;
  // Auto-detect: if running in browser and NOT on localhost, use serverless
  // (Vercel, preview deploys, custom domains all need serverless — only localhost has a backend)
  if (typeof window !== 'undefined') {
    const host = window.location.hostname;
    if (host !== 'localhost' && host !== '127.0.0.1') return true;
  }
  // Server-side: Vercel sets VERCEL=1 automatically
  if (process.env.VERCEL === '1') return true;
  // Default: local development → traditional mode
  return false;
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
  ligand_smiles?: string; // SMILES string for ligand-aware prediction
  metal?: string;         // Metal element code (TB, ZN, CA, etc.) for metal-binding validation
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
  result?: unknown;
  error?: string;
  syncCompleted?: boolean;
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
  type: 'rfd3' | 'rf3' | 'mpnn' | 'workflow';
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

export interface ConservationGradeAPI {
  position: number;
  residue: string;
  grade: number;         // 1-9 ConSurf grade
  score: number;
  confidence_interval: [number, number];
  data_quality: string;
}

export interface ConservationAnalysisResponse {
  status: string;
  grades: ConservationGradeAPI[];
  highly_conserved_positions: number[];
  conserved_positions?: number[];
  variable_positions: number[];
  msa_depth: number;
  method: string;
  average_conservation?: number;
  reliable?: boolean;
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

  // Get the runsync endpoint URL based on current mode
  private getRunsyncEndpoint(): string {
    if (this.mode === 'serverless') {
      return '/api/runpod/runsync';
    }
    return `/api/traditional/runsync?url=${encodeURIComponent(this.baseUrl)}`;
  }

  /**
   * Submit a task to RunPod asynchronously and poll for the result.
   * Avoids Vercel function timeout on Hobby plan by splitting the
   * blocking runsync call into: async submit (fast) + client-side polling (fast per request).
   * Falls back to runsync in traditional (local) mode where there is no timeout issue.
   */
  private async runpodAsync(input: Record<string, unknown>): Promise<any> {
    if (this.mode !== 'serverless') {
      // Traditional mode: runsync works fine (no Vercel timeout)
      const response = await fetch(this.getRunsyncEndpoint(), {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ input }),
      });
      if (!response.ok) {
        const errText = await response.text().catch(() => response.statusText);
        throw new Error(`Backend error (${response.status}): ${errText}`);
      }
      return response.json();
    }

    // Serverless mode: async submit + poll
    const submitRes = await fetch('/api/runpod/async', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ input }),
    });
    if (!submitRes.ok) {
      const errText = await submitRes.text().catch(() => submitRes.statusText);
      throw new Error(`Async submit failed (${submitRes.status}): ${errText}`);
    }
    const { job_id } = await submitRes.json();
    console.log('[API] Async job submitted:', job_id, '| task:', input.task);

    // Poll for result
    const maxPollTime = 5 * 60 * 1000; // 5 minutes
    const pollInterval = 2000;
    const start = Date.now();

    while (Date.now() - start < maxPollTime) {
      await new Promise(r => setTimeout(r, pollInterval));

      const statusRes = await fetch(`/api/runpod/jobs/${job_id}`);
      if (!statusRes.ok) {
        console.warn('[API] Poll error for', job_id, statusRes.status);
        continue;
      }

      const status = await statusRes.json();
      if (status.status === 'completed' || status.status === 'failed') {
        // Re-wrap into the same shape runsync returns: { output: { status, result, error } }
        return {
          id: job_id,
          status: status.status === 'completed' ? 'COMPLETED' : 'FAILED',
          output: {
            status: status.status,
            result: status.result,
            error: status.error,
          },
        };
      }
    }

    throw new Error(`Job ${job_id} timed out after ${Math.round(maxPollTime / 60000)} minutes`);
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
        const userId = await getCurrentUserId();
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
          user_id: userId,
        });
      }

      return result;
    } else {
      // Traditional mode uses Next.js proxy to avoid CORS issues
      // runsync is SYNCHRONOUS - result comes back in same request
      const proxyUrl = this.getRunsyncEndpoint();
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
        const userId = await getCurrentUserId();
        await supabaseSaveJob({
          runpod_id: result.job_id,
          type: 'rf3',
          request: { sequence_length: request.sequence?.length },
          user_id: userId,
        });
      }

      return result;
    }

    // Traditional mode uses Next.js proxy to avoid CORS issues
    const proxyUrl = this.getRunsyncEndpoint();
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
        const userId = await getCurrentUserId();
        await supabaseSaveJob({
          runpod_id: result.job_id,
          type: 'mpnn',
          request: {
            num_sequences: request.num_sequences,
            temperature: request.temperature,
            model_type: request.model_type,
          },
          user_id: userId,
        });
      }

      return result;
    }

    // Traditional mode uses Next.js proxy to avoid CORS issues
    const proxyUrl = this.getRunsyncEndpoint();
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

  // Workflow execution via modular pipeline
  async submitWorkflow(spec: import('./workflow-types').WorkflowSpec): Promise<JobResponse & { result?: any; error?: string; syncCompleted?: boolean }> {
    const jobId = `workflow-${Date.now()}-${Math.random().toString(36).slice(2, 8)}`;

    // Traditional mode: synchronous call through RunPod
    const proxyUrl = this.getRunsyncEndpoint();
    const response = await fetch(proxyUrl, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ input: { task: 'workflow_run', ...spec } }),
    });

    if (!response.ok) {
      const error = await response.json().catch(() => ({ detail: response.statusText }));
      return {
        job_id: jobId,
        status: 'failed',
        message: error.detail || 'Workflow submission failed',
        error: error.detail,
      };
    }

    const result = await response.json();
    const output = result.output || {};

    const jobResult = {
      type: 'workflow' as const,
      data: output.result || output,
    };
    this.syncJobResults.set(jobId, jobResult);

    return {
      job_id: jobId,
      status: output.status === 'completed' ? 'completed' : 'failed',
      message: output.status === 'completed' ? 'Workflow completed' : (output.error || 'Workflow failed'),
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
    const proxyUrl = this.getRunsyncEndpoint();
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
    const response = await fetch(this.getRunsyncEndpoint(), {
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

  // Fetch PDB from RCSB via our proxy (to avoid CORS issues)
  async fetchPdb(pdbId: string): Promise<FetchPdbResponse> {
    const pdbIdUpper = pdbId.toUpperCase();
    const proxyUrl = `/api/rcsb/${pdbIdUpper}`;

    try {
      const response = await fetch(proxyUrl);
      if (!response.ok) {
        const errorData = await response.json().catch(() => ({}));
        throw new Error(errorData.error || `PDB ${pdbIdUpper} not found`);
      }

      const content = await response.text();
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

    const response = await fetch(this.getRunsyncEndpoint(), {
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

    const response = await fetch(this.getRunsyncEndpoint(), {
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

    const response = await fetch(this.getRunsyncEndpoint(), {
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

  // ============== Unified Design Analysis ==============

  /**
   * Run UnifiedDesignAnalyzer + FILTER_PRESETS on a single design PDB.
   * Returns flat metrics, filter pass/fail, and full nested analysis data.
   */
  async analyzeDesign(request: {
    pdb_content: string;
    backbone_pdb?: string;
    metal_type?: string;
    ligand_name?: string;
    design_type?: string;
    filter_tier?: string;
    design_params?: Record<string, unknown>;
    // RF3/AF3 ligand-interface confidence (bioRxiv 2025.09.18.676967v2)
    rf3_confidences?: Record<string, unknown>;
    iptm?: number;
    min_chain_pair_pae?: number;
  }): Promise<{
    design_type: string;
    metrics: Record<string, number>;
    filter_preset: string;
    filter_passed: boolean;
    filter_tier?: string;
    chemistry_context?: {
      metal: string;
      hsab_class: string;
      formal_cn_range: [number, number];
      tier: string;
    };
    failed_filters: Array<{ metric: string; value: number; threshold: Record<string, number> }>;
    auto_detected: Record<string, unknown>;
    analyses: Record<string, unknown>;
  }> {
    const response = await fetch(this.getRunsyncEndpoint(), {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        input: {
          task: 'analyze_design',
          ...request,
        }
      }),
    });

    if (!response.ok) {
      throw new Error(`Design analysis failed: ${response.statusText}`);
    }

    const result = await response.json();
    if (result.output?.status === 'failed') {
      throw new Error(result.output.error || 'Design analysis failed');
    }

    return result.output?.result || result.result;
  }

  // ============== Hotspot Detection ==============

  /**
   * Detect optimal binding hotspots for protein binder design.
   * Uses SASA analysis and spatial clustering to find binding patches.
   */
  async detectHotspots(request: DetectHotspotsRequest): Promise<DetectHotspotsResponse> {
    console.log('[API] Detecting hotspots for target chain:', request.target_chain || 'A');

    // Traditional mode uses Next.js proxy
    const response = await fetch(this.getRunsyncEndpoint(), {
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

  /**
   * Analyze evolutionary conservation using ConSurf methodology.
   * Returns per-residue conservation grades (1-9 scale).
   */
  async analyzeConservation(request: { pdb_content: string; chain?: string }): Promise<ConservationAnalysisResponse> {
    console.log('[API] Analyzing conservation for chain:', request.chain || 'A');

    const response = await fetch(this.getRunsyncEndpoint(), {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        input: {
          task: 'analyze_conservation',
          ...request,
        }
      }),
    });

    if (!response.ok) {
      throw new Error(`Failed to analyze conservation: ${response.statusText}`);
    }

    const result = await response.json();
    console.log('[API] Conservation analysis result:', result.output?.status);

    if (result.output?.status === 'error' || result.output?.status === 'failed') {
      throw new Error(result.output.error || 'Conservation analysis failed');
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
        job_id: jobId,
        status: 'completed',
        created_at: new Date().toISOString(),
        result: cachedResult,
      } as JobStatus;
    }

    // For any job ID without cached result, return completed with null
    // (result was already consumed or job didn't cache)
    return {
      job_id: jobId,
      status: 'completed',
      created_at: new Date().toISOString(),
      result: undefined,
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

  // Poll job status until completion with retry and exponential backoff
  async waitForJob(
    jobId: string,
    onStatusChange?: (status: JobStatus) => void,
    pollInterval: number = 2000
  ): Promise<JobStatus> {
    console.log(`[API] Starting to poll job ${jobId}...`);
    const startTime = Date.now();
    let pollCount = 0;
    let consecutiveErrors = 0;
    const MAX_CONSECUTIVE_ERRORS = 5;
    const MAX_POLL_DURATION = 30 * 60 * 1000; // 30 minutes

    while (true) {
      // Timeout guard
      const elapsed = Date.now() - startTime;
      if (elapsed > MAX_POLL_DURATION) {
        throw new Error(`Job polling timed out after ${Math.round(elapsed / 60000)} minutes`);
      }

      pollCount++;

      try {
        const status = await this.getJobStatus(jobId);
        consecutiveErrors = 0; // Reset on success

        // Log every 10th poll
        if (pollCount % 10 === 1) {
          const elapsedSec = Math.round(elapsed / 1000);
          console.log(`[API] Job ${jobId} status: ${status.status} (poll #${pollCount}, ${elapsedSec}s elapsed)`);
        }

        if (onStatusChange) {
          onStatusChange(status);
        }

        if (status.status === 'completed' || status.status === 'failed') {
          const elapsedSec = Math.round(elapsed / 1000);
          console.log(`[API] Job ${jobId} finished with status: ${status.status} after ${elapsedSec}s (${pollCount} polls)`);
          return status;
        }
      } catch (err) {
        consecutiveErrors++;
        const errMsg = err instanceof Error ? err.message : String(err);
        console.warn(`[API] Poll error for job ${jobId} (${consecutiveErrors}/${MAX_CONSECUTIVE_ERRORS}): ${errMsg}`);

        if (consecutiveErrors >= MAX_CONSECUTIVE_ERRORS) {
          throw new Error(`Lost connection to backend after ${MAX_CONSECUTIVE_ERRORS} consecutive poll failures: ${errMsg}`);
        }
      }

      // Exponential backoff on errors: 2s, 4s, 8s, 16s, 32s (capped)
      const backoff = Math.min(pollInterval * Math.pow(2, consecutiveErrors), 32000);
      await new Promise(resolve => setTimeout(resolve, backoff));
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
    const response = await fetch(this.getRunsyncEndpoint(), {
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
    const response = await fetch(this.getRunsyncEndpoint(), {
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
    const response = await fetch(this.getRunsyncEndpoint(), {
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

  // ============== Pipeline Methods ==============

  /**
   * Start a pipeline design run (sweep or production mode).
   */
  async startPipeline(request: {
    mode: 'sweep' | 'production';
    metal: string;
    ligand: string;
    motif_pdb: string;
    sweep_configs?: Array<{
      name: string;
      contig_size: string;
      contig_range: string;
      cfg_scale: number;
      num_designs: number;
    }>;
    production_config?: {
      contig_size: string;
      contig_range: string;
      cfg_scale: number;
    };
    num_designs?: number;
    filters?: {
      plddt: number;
      ptm: number;
      pae: number;
    };
    designs_per_config?: number;
  }): Promise<{
    session_id: string;
    status: string;
    mode: string;
    total_configs: number;
    designs_per_config?: number;
    num_designs?: number;
  }> {
    console.log('[API] Starting pipeline:', request.mode);

    const response = await fetch(this.getRunsyncEndpoint(), {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        input: { task: 'pipeline_design', ...request }
      }),
    });

    if (!response.ok) {
      throw new Error(`Failed to start pipeline: ${response.statusText}`);
    }

    const result = await response.json();
    if (result.output?.status === 'failed') {
      throw new Error(result.output.error || 'Pipeline start failed');
    }

    return result.output?.result || result.output;
  }

  /**
   * Poll pipeline status for progress updates.
   */
  async getPipelineStatus(sessionId: string): Promise<{
    session_id: string;
    status: 'running' | 'completed' | 'cancelled' | 'failed' | 'not_found';
    mode: string;
    current_config: number;
    total_configs: number;
    current_design: number;
    designs_per_config: number;
    total_generated: number;
    total_passing: number;
    total_review: number;
    total_failed: number;
    pass_rate: number;
    best_design: {
      name: string;
      plddt: number;
      ptm: number;
      pae: number;
    } | null;
    results?: Array<{
      name: string;
      config_name: string;
      plddt: number;
      ptm: number;
      pae: number;
      status: string;
    }>;
    summary?: {
      total_generated: number;
      total_passing: number;
      pass_rate: number;
    };
  }> {
    const response = await fetch(this.getRunsyncEndpoint(), {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        input: { task: 'pipeline_status', session_id: sessionId }
      }),
    });

    if (!response.ok) {
      throw new Error(`Failed to get pipeline status: ${response.statusText}`);
    }

    const result = await response.json();
    if (result.output?.status === 'failed') {
      throw new Error(result.output.error || 'Pipeline status failed');
    }

    return result.output?.result || result.output;
  }

  /**
   * Poll pipeline status at regular intervals until completion.
   */
  async pollPipelineStatus(
    sessionId: string,
    onProgress?: (status: Awaited<ReturnType<typeof this.getPipelineStatus>>) => void,
    pollInterval: number = 2500,
  ): Promise<Awaited<ReturnType<typeof this.getPipelineStatus>>> {
    console.log(`[API] Starting to poll pipeline ${sessionId}...`);
    const startTime = Date.now();
    let pollCount = 0;

    while (true) {
      pollCount++;
      const status = await this.getPipelineStatus(sessionId);

      // Log every 10th poll
      if (pollCount % 10 === 1) {
        const elapsed = Math.round((Date.now() - startTime) / 1000);
        console.log(`[API] Pipeline ${sessionId} status: ${status.status} (poll #${pollCount}, ${elapsed}s elapsed)`);
      }

      if (onProgress) {
        onProgress(status);
      }

      if (status.status === 'completed' || status.status === 'cancelled' || status.status === 'failed') {
        const elapsed = Math.round((Date.now() - startTime) / 1000);
        console.log(`[API] Pipeline ${sessionId} finished with status: ${status.status} after ${elapsed}s`);
        return status;
      }

      await new Promise(resolve => setTimeout(resolve, pollInterval));
    }
  }

  /**
   * Cancel a running pipeline.
   */
  async cancelPipeline(sessionId: string): Promise<{ status: string; session_id: string }> {
    console.log(`[API] Cancelling pipeline ${sessionId}`);

    const response = await fetch(this.getRunsyncEndpoint(), {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        input: { task: 'pipeline_cancel', session_id: sessionId }
      }),
    });

    if (!response.ok) {
      throw new Error(`Failed to cancel pipeline: ${response.statusText}`);
    }

    const result = await response.json();
    return result.output?.result || result.output;
  }

  /**
   * Export pipeline results as FASTA.
   */
  async exportPipelineFasta(
    sessionId: string,
    includeReview: boolean = false,
  ): Promise<{ status: string; filename: string; content: string; num_sequences: number }> {
    console.log(`[API] Exporting pipeline ${sessionId} to FASTA`);

    const response = await fetch(this.getRunsyncEndpoint(), {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        input: {
          task: 'pipeline_export',
          session_id: sessionId,
          include_review: includeReview,
        }
      }),
    });

    if (!response.ok) {
      throw new Error(`Failed to export pipeline: ${response.statusText}`);
    }

    const result = await response.json();
    if (result.output?.status === 'failed') {
      throw new Error(result.output.error || 'Pipeline export failed');
    }

    return result.output?.result || result.output;
  }

  // ============== Metal Binding Design Methods (Round 7b Style) ==============

  /**
   * Analyze a PDB structure for metal binding site design potential.
   */
  async analyzeMetalBindingSite(request: {
    pdb_content?: string;
    pdb_id?: string;
  }): Promise<{
    status: string;
    result: {
      metals: Array<{ element: string; chain: string; res_num: string }>;
      ligands: Array<{ name: string; chain: string; res_num: string }>;
      num_metals: number;
      num_ligands: number;
      suitable_for_design: boolean;
      recommended_workflow: string | null;
      suggestions: string[];
      recommended_hsab_bias?: string;
    };
  }> {
    console.log('[API] Analyzing metal binding site');

    const response = await fetch(this.getRunsyncEndpoint(), {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        input: {
          task: 'metal_binding_design',
          mode: 'analyze',
          ...request,
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

    return result.output || result;
  }

  /**
   * Run parameter sweep for metal binding design (9 configs: 3 sizes x 3 CFG scales).
   */
  async startMetalBindingSweep(request: {
    motif_pdb: string;
    metal: string;
    ligand?: string;
    sizes?: string[];
    cfg_scales?: number[];
    designs_per_config?: number;
    filters?: { plddt: number; ptm: number; pae: number };
    seed?: number;
  }): Promise<{
    status: string;
    result: {
      session_id: string;
      mode: string;
      total_configs: number;
      designs_per_config: number;
      total_generated: number;
      total_passing: number;
      total_review: number;
      pass_rate: number;
      best_design: {
        name: string;
        plddt: number;
        ptm: number;
        tier: string;
      } | null;
      config_rankings: Array<{
        rank: number;
        config_name: string;
        contig_range: string;
        cfg_scale: number;
        pass_rate: number;
        avg_plddt: number;
        avg_ptm: number;
        tier_distribution: Record<string, number>;
      }>;
      results: Array<{
        name: string;
        config_name: string;
        sequence: string;
        plddt: number;
        ptm: number;
        pae: number;
        tier: string;
        status: string;
        pdb_content?: string;
        seed?: number;
        coordination_number?: number;
        geometry_rmsd?: number;
      }>;
      config_errors?: Record<string, string[]> | null;
    };
  }> {
    console.log('[API] Starting metal binding sweep');

    // Sweep can take 5-15 minutes depending on config count; use extended timeout
    const controller = new AbortController();
    const timeoutId = setTimeout(() => controller.abort(), 15 * 60 * 1000);

    const response = await fetch(this.getRunsyncEndpoint(), {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        input: {
          task: 'metal_binding_design',
          mode: 'sweep',
          ...request,
        }
      }),
      signal: controller.signal,
    });

    clearTimeout(timeoutId);

    if (!response.ok) {
      throw new Error(`Failed to start sweep: ${response.statusText}`);
    }

    const result = await response.json();
    if (result.output?.status === 'failed') {
      throw new Error(result.output.error || 'Sweep failed');
    }

    return result.output || result;
  }

  /**
   * Run single metal binding RFD3 design — backbone only, no MPNN/RF3.
   * Returns backbone PDBs for downstream MPNN/RF3 steps to process separately.
   */
  async submitMetalSingleDesign(request: {
    motif_pdb: string;
    metal: string;
    ligand?: string;
    contig?: string;
    cfg_scale?: number;
    num_designs?: number;
    num_timesteps?: number;
    step_scale?: number;
    gamma_0?: number;
    bury_ligand?: boolean;
    design_goal?: string;
    seed?: number;
  }): Promise<{
    status: string;
    result: {
      designs: Array<{
        content: string;
        name: string;
        seed?: number;
      }>;
      total_generated: number;
      total_requested: number;
      errors?: string[] | null;
      metal: string;
      ligand?: string;
    };
  }> {
    console.log('[API] Starting metal single design');

    const controller = new AbortController();
    const timeoutId = setTimeout(() => controller.abort(), 10 * 60 * 1000);

    const response = await fetch(this.getRunsyncEndpoint(), {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        input: {
          task: 'metal_binding_design',
          mode: 'single',
          ...request,
        }
      }),
      signal: controller.signal,
    });

    clearTimeout(timeoutId);

    if (!response.ok) {
      throw new Error(`Failed to start metal single design: ${response.statusText}`);
    }

    const result = await response.json();
    if (result.output?.status === 'failed') {
      throw new Error(result.output.error || 'Metal single design failed');
    }

    return result.output || result;
  }

  /**
   * Search RCSB PDB for scaffold candidates containing a metal-ligand complex.
   * Reusable across any pipeline that needs to find existing structures.
   */
  async searchScaffold(request: {
    metal: string;
    ligand_name?: string;
    ligand_code?: string;
    ligand_smiles?: string;
    resolution_max?: number;
    limit?: number;
    fetch_pdb?: boolean;
  }): Promise<{
    searched: boolean;
    query_metal: string;
    query_ligand: string;
    ligand_code: string;
    num_pdb_hits: number;
    num_validated: number;
    candidates: Array<{
      pdb_id: string;
      source_metal: string;
      target_metal: string;
      ligand_code: string;
      needs_substitution: boolean;
      coordination_number: number;
      total_score: number;
    }>;
    best_candidate: {
      pdb_id: string;
      source_metal: string;
      target_metal: string;
      ligand_code: string;
      needs_substitution: boolean;
      coordination_number: number;
      total_score: number;
    } | null;
    recommended_action: 'scaffold' | 'de_novo';
    reason: string;
    best_pdb_content?: string;
    candidate_pdbs?: Record<string, string | null>;
  }> {
    const endpoint = this.getRunsyncEndpoint();
    console.log('[API] Searching for scaffold candidates:', request.metal, request.ligand_name || request.ligand_code, '| mode:', this.mode, '| endpoint:', endpoint);

    const response = await fetch(endpoint, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        input: {
          task: 'scaffold_search',
          ...request,
        }
      }),
    });

    if (!response.ok) {
      const errorBody = await response.text().catch(() => '');
      console.error('[API] Scaffold search HTTP error:', response.status, response.statusText, errorBody);
      throw new Error(`Scaffold search failed (${response.status}): ${errorBody || response.statusText}`);
    }

    const result = await response.json();
    console.log('[API] Scaffold search response:', JSON.stringify(result).substring(0, 500));

    if (result.output?.status === 'failed') {
      console.error('[API] Scaffold search backend error:', result.output.error);
      throw new Error(result.output.error || 'Scaffold search failed');
    }

    // Handle RunPod timeout (runsync returned before completion)
    if (result.status === 'IN_QUEUE' || result.status === 'IN_PROGRESS') {
      console.warn('[API] Scaffold search timed out on RunPod runsync, job:', result.id);
      throw new Error(`Scaffold search timed out (job ${result.id}). RunPod runsync limit exceeded.`);
    }

    return result.output?.result || result.result;
  }

  /**
   * Resolve a ligand name to SMILES, residue code, and source via backend.
   * Uses template library, PDB, PubChem, and isomeric SMILES table.
   */
  async resolveLigand(request: {
    ligand_name: string;
    metal_type?: string;
    isomer_spec?: string;
  }): Promise<{
    success: boolean;
    smiles: string;
    residue_code: string;
    source: string;
    name: string;
    warnings: string[];
    ligand_fixing_strategy?: string;
    coordination_donors?: string[];
  }> {
    console.log('[API] Resolving ligand:', request.ligand_name, '| metal:', request.metal_type || 'none');

    const result = await this.runpodAsync({
      task: 'resolve_ligand',
      ...request,
    });

    if (result.output?.status === 'failed') {
      throw new Error(result.output.error || 'Ligand resolution failed');
    }

    return result.output?.result || result.result;
  }

  /**
   * Analyze ligand features using three-layer system:
   * 1. Knowledge base (cached PDB crystal structure data)
   * 2. ChemicalFeatures (RDKit pharmacophore perception)
   * 3. Geometry filter (SMARTS-based fallback)
   */
  async analyzeLigandFeatures(request: {
    ligand_name: string;
    smiles?: string;
    metal_type?: string;
    pdb_content?: string;
    record_evidence?: Record<string, unknown>;
  }): Promise<{
    ligand_name: string;
    smiles: string;
    source: 'knowledge_base' | 'chemicalfeatures' | 'geometry_filter' | 'unknown';
    metal: string | null;
    features: Array<{
      atom_idx: number;
      atom_name: string;
      element: string;
      type: 'donor' | 'acceptor' | 'aromatic' | 'hydrophobic';
      is_coordination_donor: boolean;
      coords: [number, number, number] | null;
      hsab: string | null;
      enabled: boolean;
    }>;
    coordination_donors: string[];
    max_denticity: number;
    evidence_count: number;
    compatibility_score: number;
    coordination_mode: string;
    notes: string;
    pdb_evidence: string[];
    ligand_pdb_content: string | null;
  }> {
    console.log('[API] Analyzing ligand features:', request.ligand_name, '| metal:', request.metal_type || 'none');

    const result = await this.runpodAsync({
      task: 'analyze_ligand_features',
      ...request,
    });

    console.log('[API] analyzeLigandFeatures response:', JSON.stringify(result).slice(0, 300));

    // Check for top-level failure (traditional backend format)
    if (result.status === 'FAILED' || result.status === 'failed') {
      console.error('[API] analyzeLigandFeatures backend error:', result.error);
      throw new Error(result.error || 'Ligand feature analysis failed');
    }

    // Check for nested failure (RunPod serverless format)
    if (result.output?.status === 'failed') {
      console.error('[API] analyzeLigandFeatures backend error:', result.output.error);
      throw new Error(result.output.error || 'Ligand feature analysis failed');
    }

    const data = result.output?.result || result.result;
    if (!data) {
      console.warn('[API] analyzeLigandFeatures: unexpected response shape:', JSON.stringify(result).slice(0, 300));
      throw new Error('Backend returned unexpected response format for ligand feature analysis');
    }
    return data;
  }

  /**
   * Parse a natural language design query using the backend AI parser (Claude API).
   * Falls back to keyword-based parsing if no API key is configured.
   */
  async parseIntent(query: string): Promise<{
    metal_type?: string;
    ligand_name?: string;
    isomer_specification?: string;
    design_goal: string;
    target_topology: string;
    chain_length_min: number;
    chain_length_max: number;
    design_mode: string;
    source_pdb_id?: string;
    pocket_description?: string;
    include_all_contacts: boolean;
    stability_focus: boolean;
    enzyme_class?: string;
    enzyme_class_confidence: number;
    preserve_function: boolean;
    bury_ligand: boolean;
    confidence: number;
    raw_query: string;
    corrected_query: string;
    warnings: string[];
    suggestions: string[];
    parsed_entities: Record<string, unknown>;
    typo_corrections: string[];
    parser_type: 'ai' | 'fallback';
  }> {
    // Call the Vercel edge function directly — no RunPod / GPU needed
    const endpoint = '/api/ai-parse';
    console.log('[API] Parsing intent with AI:', query.substring(0, 80), '| endpoint:', endpoint);

    const response = await fetch(endpoint, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        input: { query },
      }),
    });

    if (!response.ok) {
      const errorBody = await response.text().catch(() => '');
      console.error('[API] AI parse HTTP error:', response.status, response.statusText, errorBody);
      throw new Error(`AI parse failed (${response.status}): ${errorBody || response.statusText}`);
    }

    const result = await response.json();
    console.log('[API] AI parse response status:', result.output?.status, '| has result:', !!result.output?.result);

    if (result.error) {
      console.error('[API] AI parse error in response:', result.error);
      throw new Error(`AI parse failed: ${result.error}`);
    }

    if (result.output?.status === 'failed') {
      console.error('[API] AI parse output failed:', result.output.error);
      throw new Error(result.output.error || 'AI parse failed');
    }

    return result.output?.result || result.result;
  }

  /**
   * Run production mode for metal binding design with best or specified config.
   */
  async startMetalBindingProduction(request: {
    motif_pdb: string;
    metal: string;
    ligand?: string;
    session_id?: string;
    config?: { contig_range: string; cfg_scale: number };
    num_designs?: number;
  }): Promise<{
    status: string;
    result: {
      session_id: string;
      mode: string;
      config: { contig_range: string; cfg_scale: number };
      num_designs: number;
      message: string;
    };
  }> {
    console.log('[API] Starting metal binding production');

    const response = await fetch(this.getRunsyncEndpoint(), {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        input: {
          task: 'metal_binding_design',
          mode: 'production',
          ...request,
        }
      }),
    });

    if (!response.ok) {
      throw new Error(`Failed to start production: ${response.statusText}`);
    }

    const result = await response.json();
    if (result.output?.status === 'failed') {
      throw new Error(result.output.error || 'Production start failed');
    }

    return result.output || result;
  }

  // ============== Scout Filter, Design History, Lesson Detection ==============

  /**
   * Filter backbones by generating 1 scout sequence each and validating.
   * Returns only passing backbones with scout metrics.
   */
  async scoutFilter(request: {
    backbone_pdbs: string[];
    ptm_threshold?: number;
    plddt_threshold?: number;
    target_metal?: string;
    ligand_smiles?: string;
    ligand_name?: string;
    pre_filter_only?: boolean;
  }): Promise<{
    filtered_pdbs: string[];
    original_count: number;
    passing_count: number;
    scout_results: Array<{
      backbone_index: number;
      ptm: number;
      plddt: number;
      passed: boolean;
      sequence: string;
      pre_filter_failed?: boolean;
    }>;
    ptm_threshold: number;
    plddt_threshold: number;
    pre_filter_results?: Array<{
      passed: boolean;
      checks: Record<string, {
        passed: boolean;
        reasons?: string[];
        chain_breaks?: Array<{ chain: string; res_i: number; res_j: number; distance: number }>;
        coordination_number?: number;
        geometry_rmsd?: number;
        ligand_name?: string;
        atom_count?: number;
        error?: string;
      }>;
      failed_checks: string[];
      skipped_checks: string[];
    }>;
    pre_filter_passed?: number;
    pre_filter_failed?: number;
  }> {
    console.log('[API] Running scout filter on', request.backbone_pdbs.length, 'backbones');

    const response = await fetch(this.getRunsyncEndpoint(), {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        input: { task: 'scout_filter', ...request },
      }),
    });

    if (!response.ok) {
      throw new Error(`Scout filter failed: ${response.statusText}`);
    }

    const result = await response.json();
    if (result.output?.status === 'failed') {
      throw new Error(result.output.error || 'Scout filter failed');
    }

    return result.output?.result || result.result;
  }

  /**
   * Save a completed pipeline run to persistent design history.
   */
  async saveDesignHistory(request: {
    session_name: string;
    design_params: Record<string, unknown>;
    design_outputs: Record<string, unknown>;
    design_metrics: Record<string, unknown>;
  }): Promise<{ run_id: string; session_id: string }> {
    console.log('[API] Saving design history:', request.session_name);

    const response = await fetch(this.getRunsyncEndpoint(), {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        input: { task: 'save_design_history', ...request },
      }),
    });

    if (!response.ok) {
      throw new Error(`Save design history failed: ${response.statusText}`);
    }

    const result = await response.json();
    const output = result.output || {};
    if (output.status === 'failed') {
      throw new Error(output.error || 'Save design history failed');
    }

    const inner = output.result || result.result || {};
    return {
      run_id: inner.run_id || `local-${Date.now()}`,
      session_id: inner.session_id || 'unknown',
    };
  }

  /**
   * Check if this design run triggers any lesson patterns.
   */
  async checkLessons(request: {
    result: Record<string, unknown>;
  }): Promise<{
    trigger_detected: boolean;
    trigger?: {
      type: 'failure_pattern' | 'breakthrough' | 'improvement';
      description: string;
      relevant_designs: string[];
      metrics_involved: string[];
    };
    history_count: number;
  }> {
    console.log('[API] Checking for lesson triggers');

    const response = await fetch(this.getRunsyncEndpoint(), {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        input: { task: 'check_lessons', ...request },
      }),
    });

    if (!response.ok) {
      throw new Error(`Check lessons failed: ${response.statusText}`);
    }

    const result = await response.json();
    if (result.output?.status === 'failed') {
      throw new Error(result.output.error || 'Check lessons failed');
    }

    return result.output?.result || result.result;
  }

  /**
   * Get status of a metal binding design session.
   */
  async getMetalBindingStatus(sessionId: string): Promise<{
    status: string;
    result: {
      session_id: string;
      mode: string;
      status: string;
      current_config: number;
      total_configs: number;
      current_design: number;
      designs_per_config: number;
      total_generated: number;
      total_passing: number;
      total_review: number;
      total_failed: number;
      pass_rate: number;
      best_design: { name: string; plddt: number } | null;
      config_rankings: Array<{ rank: number; config_name: string; pass_rate: number }>;
    };
  }> {
    const response = await fetch(this.getRunsyncEndpoint(), {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        input: {
          task: 'metal_binding_design',
          mode: 'status',
          session_id: sessionId,
        }
      }),
    });

    if (!response.ok) {
      throw new Error(`Failed to get status: ${response.statusText}`);
    }

    const result = await response.json();
    return result.output || result;
  }
}

export const api = new FoundryAPI();
export default api;
