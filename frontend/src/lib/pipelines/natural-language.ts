/**
 * Natural Language Pipeline
 * Steps: parse_intent → resolve_structure → scaffold_search → configure → structure_gen → mpnn_sequence → validation → analysis
 *
 * This pipeline handles free-form natural language requests by parsing user intent,
 * resolving the target structure, and configuring the appropriate design parameters.
 */

import { MessageSquare, Search, Settings, Building2, Shield, Dna, BarChart3 } from 'lucide-react';
import api from '@/lib/api';
import { extractPdbContent } from '@/lib/workflowHandlers';
import type { PipelineDefinition, PdbOutput, SequenceOutput } from '@/lib/pipeline-types';
import {
  createMpnnStep,
  createRf3Step,
  createScaffoldSearchStep,
  createScoutFilterStep,
  createSaveHistoryStep,
  createCheckLessonsStep,
} from './shared-steps';
import { IntentResultPreview } from '@/components/pipeline/IntentResultPreview';
import { SweepResultPreview } from '@/components/pipeline/SweepResultPreview';

// ============================================================================
// AI Parse Helpers
// ============================================================================

/** SMILES lookup table for known ligands */
const LIGAND_SMILES: Record<string, string> = {
  citrate: 'OC(=O)CC(O)(CC(O)=O)C(O)=O',
  atp: 'c1nc(c2c(n1)n(cn2)[C@@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N',
  glucose: 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O',
  azobenzene: 'c1ccc(/N=N/c2ccccc2)cc1',
};

/** Map AI DesignIntent fields to the pipeline's design_type enum */
function mapDesignType(ai: { metal_type?: string; design_mode?: string; target_topology?: string; design_goal?: string }): string {
  if (ai.metal_type) return 'metal';
  if (ai.design_mode === 'scaffold') return 'scaffold';
  if (ai.target_topology?.includes('dimer')) return 'ligand';
  if (ai.design_goal === 'binding') return 'binder';
  return 'general';
}

/** Lookup SMILES for a ligand name */
function lookupSmiles(name?: string): string {
  if (!name) return '';
  return LIGAND_SMILES[name.toLowerCase()] || '';
}

/** Build a human-readable description from AI parsed intent */
function buildDescription(ai: { metal_type?: string; ligand_name?: string; design_goal?: string; enzyme_class?: string; design_mode?: string }): string {
  const parts: string[] = [];
  if (ai.enzyme_class) {
    parts.push(`${ai.enzyme_class} design`);
  } else if (ai.design_goal && ai.design_goal !== 'binding') {
    parts.push(`${ai.design_goal} design`);
  }
  if (ai.metal_type && ai.ligand_name) {
    parts.push(`Metal-ligand binding: ${ai.metal_type} + ${ai.ligand_name}`);
  } else if (ai.metal_type) {
    parts.push(`Metal binding: ${ai.metal_type}`);
  } else if (ai.ligand_name) {
    parts.push(`Ligand binding: ${ai.ligand_name}`);
  }
  if (ai.design_mode === 'scaffold') {
    parts.push('(scaffold mode)');
  }
  return parts.length > 0 ? parts.join(' — ') : 'General protein design';
}

/** Keyword-based fallback parser (original logic, used when backend AI is unavailable) */
function keywordFallbackParse(prompt: string, pdbId?: string): Record<string, unknown> {
  const promptLower = prompt.toLowerCase();
  let designType = 'general';
  let targetMetal = '';
  let ligandName = '';
  let ligandSmiles = '';
  let description = 'General protein design';

  // 1) Detect metals
  const metalNameMap: Record<string, string> = {
    zinc: 'ZN', iron: 'FE', copper: 'CU', manganese: 'MN',
    calcium: 'CA', magnesium: 'MG', cobalt: 'CO', nickel: 'NI',
    terbium: 'TB', lanthanide: 'TB', cadmium: 'CD', lead: 'PB',
    chromium: 'CR', molybdenum: 'MO', vanadium: 'V',
    europium: 'EU', gadolinium: 'GD', lanthanum: 'LA', cerium: 'CE',
    neodymium: 'ND',
  };
  const metalSymbolMap: Record<string, string> = {
    zn: 'ZN', fe: 'FE', cu: 'CU', mn: 'MN', ca: 'CA', mg: 'MG',
    co: 'CO', ni: 'NI', tb: 'TB', cd: 'CD', pb: 'PB', cr: 'CR',
    mo: 'MO', eu: 'EU', gd: 'GD', la: 'LA', ce: 'CE', nd: 'ND',
    v: 'V',
  };

  for (const [keyword, symbol] of Object.entries(metalNameMap)) {
    if (promptLower.includes(keyword)) {
      designType = 'metal';
      targetMetal = symbol;
      description = `Metal binding design targeting ${keyword}`;
      break;
    }
  }
  if (!targetMetal) {
    for (const [sym, symbol] of Object.entries(metalSymbolMap)) {
      const regex = new RegExp(`\\b${sym}\\b`, 'i');
      if (regex.test(promptLower)) {
        designType = 'metal';
        targetMetal = symbol;
        description = `Metal binding design targeting ${symbol}`;
        break;
      }
    }
  }
  if (!targetMetal && promptLower.includes('metal')) {
    designType = 'metal';
    targetMetal = 'ZN';
    description = 'Metal binding design (default: zinc)';
  }

  // 2) Detect ligands
  const ligandMap: Record<string, { name: string; smiles: string }> = {
    citrate:  { name: 'citrate',  smiles: 'OC(=O)CC(O)(CC(O)=O)C(O)=O' },
    citric:   { name: 'citrate',  smiles: 'OC(=O)CC(O)(CC(O)=O)C(O)=O' },
    atp:      { name: 'ATP',      smiles: 'c1nc(c2c(n1)n(cn2)[C@@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N' },
    nad:      { name: 'NAD+',     smiles: '' },
    glucose:  { name: 'glucose',  smiles: 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O' },
    azobenzene: { name: 'azobenzene', smiles: 'c1ccc(/N=N/c2ccccc2)cc1' },
    heme:     { name: 'heme',     smiles: '' },
  };

  for (const [keyword, info] of Object.entries(ligandMap)) {
    if (promptLower.includes(keyword)) {
      ligandName = info.name;
      ligandSmiles = info.smiles;
      break;
    }
  }

  if (targetMetal && ligandName) {
    description = `Metal-ligand binding design: ${targetMetal} + ${ligandName}`;
  }

  // 3) Detect binder requests
  if (designType === 'general' && (promptLower.includes('bind') || promptLower.includes('binder'))) {
    designType = 'binder';
    description = ligandName ? `Protein binder design for ${ligandName}` : 'Protein binder design';
  }

  // 4) Detect ligand-dimer/interface requests
  if (promptLower.includes('dimer') || promptLower.includes('interface')) {
    designType = 'ligand';
    description = 'Ligand-mediated dimer design';
  }

  // 5) Detect scaffold requests
  if (promptLower.includes('scaffold') || promptLower.includes('sweep')) {
    designType = 'scaffold';
    description = 'Scaffold sweep design';
  }

  return {
    design_type: designType,
    target_metal: targetMetal,
    ligand_name: ligandName,
    ligand_smiles: ligandSmiles,
    pdb_id: pdbId || 'Not specified',
    user_prompt: prompt,
    parsed_description: description,
    ai_parsed: false,
  };
}

// ============================================================================

export const naturalLanguagePipeline: PipelineDefinition = {
  id: 'natural-language',
  name: 'Natural Language Design',
  description: 'Design proteins from natural language descriptions',
  icon: MessageSquare,
  requiresBackend: true,

  initialParams: [
    { id: 'user_prompt', label: 'Design Request', type: 'text', required: true, helpText: 'Natural language description of the desired protein design' },
    { id: 'pdb_content', label: 'PDB Content', type: 'hidden', required: false },
    { id: 'pdb_id', label: 'PDB ID', type: 'hidden', required: false },
    // Resolved parameters (populated by earlier steps)
    { id: 'design_type', label: 'Design Type', type: 'hidden', required: false },
    { id: 'target_metal', label: 'Target Metal', type: 'hidden', required: false },
    { id: 'contig', label: 'Contig', type: 'hidden', required: false },
    { id: 'ligand', label: 'Ligand', type: 'hidden', required: false },
  ],

  steps: [
    // Step 1: Parse user intent
    {
      id: 'parse_intent',
      name: 'Parse Intent',
      description: 'Analyze the natural language request to determine design parameters',
      icon: MessageSquare,
      requiresReview: true,
      supportsSelection: false,
      optional: false,
      defaultParams: {},
      parameterSchema: [],
      ResultPreview: IntentResultPreview,

      async execute(ctx) {
        const { initialParams } = ctx;
        const prompt = initialParams.user_prompt as string || '';
        const pdbId = initialParams.pdb_id as string;

        ctx.onProgress(10, 'Analyzing design request with AI...');

        let intentData: Record<string, unknown>;
        let usedAI = false;

        try {
          // Try backend AI parser (Claude API or keyword fallback on server)
          const ai = await api.parseIntent(prompt);
          usedAI = true;
          intentData = {
            design_type: mapDesignType(ai),
            target_metal: ai.metal_type || '',
            ligand_name: ai.ligand_name || '',
            ligand_smiles: lookupSmiles(ai.ligand_name),
            pdb_id: ai.source_pdb_id || pdbId || 'Not specified',
            user_prompt: prompt,
            parsed_description: buildDescription(ai),
            // AI-specific fields
            design_goal: ai.design_goal,
            target_topology: ai.target_topology,
            chain_length_min: ai.chain_length_min,
            chain_length_max: ai.chain_length_max,
            design_mode: ai.design_mode,
            enzyme_class: ai.enzyme_class,
            enzyme_class_confidence: ai.enzyme_class_confidence,
            preserve_function: ai.preserve_function,
            confidence: ai.confidence,
            warnings: ai.warnings,
            suggestions: ai.suggestions,
            parsed_entities: ai.parsed_entities,
            typo_corrections: ai.typo_corrections,
            parser_type: ai.parser_type,
            bury_ligand: ai.bury_ligand ?? true,
            ai_parsed: true,
          };
        } catch (err) {
          // Fallback: local keyword-based parsing (works without backend)
          console.warn('[parse_intent] AI parse unavailable, using keyword fallback:', err);
          intentData = keywordFallbackParse(prompt, pdbId);
        }

        ctx.onProgress(100, usedAI ? 'AI parsing complete' : 'Intent parsed (fallback)');

        return {
          id: `intent-${Date.now()}`,
          summary: intentData.parsed_description as string,
          data: intentData,
        };
      },
    },

    // Step 2: Resolve structure
    {
      id: 'resolve_structure',
      name: 'Resolve Structure',
      description: 'Fetch or validate the target structure',
      icon: Search,
      requiresReview: true,
      supportsSelection: false,
      optional: true,
      defaultParams: {},
      parameterSchema: [],

      async execute(ctx) {
        const { initialParams, previousResults } = ctx;

        let pdbContent = initialParams.pdb_content as string;
        const pdbId = initialParams.pdb_id as string;

        if (!pdbContent && pdbId) {
          ctx.onProgress(20, `Fetching PDB ${pdbId}...`);
          const fetchResult = await api.fetchPdb(pdbId);
          pdbContent = fetchResult.content;

          ctx.onProgress(60, 'Analyzing structure...');

          return {
            id: `resolve-${Date.now()}`,
            summary: `Fetched ${pdbId}: ${fetchResult.info?.title || 'Unknown'}`,
            pdbOutputs: [{
              id: 'input-structure',
              label: pdbId,
              pdbContent,
            }],
            data: {
              pdb_id: pdbId,
              title: fetchResult.info?.title,
              chains: fetchResult.info?.chains?.join(', '),
              residues: fetchResult.info?.num_residues,
              metals: fetchResult.info?.metals?.map(m => `${m.element} (${m.chain}:${m.res_num})`).join(', ') || 'None',
            },
          };
        }

        if (pdbContent) {
          ctx.onProgress(100, 'Structure ready');
          return {
            id: `resolve-${Date.now()}`,
            summary: 'Input structure loaded',
            pdbOutputs: [{
              id: 'input-structure',
              label: pdbId || 'Uploaded structure',
              pdbContent,
            }],
            data: { source: pdbId ? 'RCSB PDB' : 'User upload' },
          };
        }

        // No PDB provided — check if this is a metal binding task
        // (template will be generated by backend or found via scaffold search)
        const intentResult = Object.values(previousResults).find(r => r.data?.design_type);
        const targetMetal = intentResult?.data?.target_metal as string;

        if (targetMetal) {
          ctx.onProgress(100, 'No input structure — metal-ligand template will be generated');
          return {
            id: `resolve-${Date.now()}`,
            summary: `No input PDB — backend will generate ${targetMetal} binding template`,
            data: { source: 'backend_template', note: 'Metal-ligand template generated by backend from coordination chemistry.' },
          };
        }

        ctx.onProgress(100, 'No input structure — using de novo design mode');
        return {
          id: `resolve-${Date.now()}`,
          summary: 'De novo design (no input structure)',
          data: { source: 'de_novo', note: 'No PDB provided. RFD3 will generate from contig specification.' },
        };
      },
    },

    // Step 3: Scaffold search — auto-discover PDB scaffolds for metal+ligand
    createScaffoldSearchStep({
      id: 'scaffold_search_nl',
      description: 'Search RCSB PDB for structures with the target metal-ligand complex (auto-skipped if not applicable)',
    }),

    // Step 4: Configure design parameters
    {
      id: 'configure',
      name: 'Configure Parameters',
      description: 'Set up design parameters based on parsed intent',
      icon: Settings,
      requiresReview: true,
      supportsSelection: false,
      optional: false,
      defaultParams: {},
      parameterSchema: [
        { id: 'num_designs', label: 'Number of Designs', type: 'slider', required: false, defaultValue: 4, range: { min: 1, max: 10, step: 1 } },
        { id: 'num_timesteps', label: 'Timesteps', type: 'slider', required: false, defaultValue: 200, range: { min: 50, max: 500, step: 50 } },
      ],

      async execute(ctx) {
        const { initialParams, previousResults, params } = ctx;

        ctx.onProgress(50, 'Configuring parameters...');

        // Gather intent from parse step
        const intentResult = Object.values(previousResults).find(r => r.data?.design_type);
        const designType = intentResult?.data?.design_type as string || initialParams.design_type as string || 'general';
        const targetMetal = intentResult?.data?.target_metal as string || initialParams.target_metal as string || '';

        // Get PDB from resolve step
        const structureResult = Object.values(previousResults).find(r => r.pdbOutputs && r.pdbOutputs.length > 0);
        const pdbContent = structureResult?.pdbOutputs?.[0]?.pdbContent;

        let configData: Record<string, unknown> = {
          design_type: designType,
          num_designs: params.num_designs ?? 4,
          num_timesteps: params.num_timesteps ?? 200,
        };

        // Try to get parameter recommendations from backend
        if (pdbContent && designType === 'metal' && targetMetal) {
          try {
            ctx.onProgress(60, 'Getting parameter recommendations...');
            const rec = await api.getParameterRecommendation({
              pdb_content: pdbContent,
              target_metal: targetMetal,
            });

            configData = {
              ...configData,
              strategy: rec.strategy,
              reasoning: rec.reasoning?.join('; '),
              ligand: rec.parameters?.ligand,
              partial_t: rec.parameters?.partial_t,
              unindex: rec.parameters?.unindex,
              gamma_0: rec.parameters?.gamma_0,
              step_scale: rec.parameters?.step_scale,
              target_metal: targetMetal,
            };
          } catch {
            // Use defaults
            configData.target_metal = targetMetal;
          }
        }

        // Gather ligand info from parse step
        const ligandName = intentResult?.data?.ligand_name as string || '';
        const ligandSmiles = intentResult?.data?.ligand_smiles as string || '';
        if (ligandName) configData.ligand_name = ligandName;
        if (ligandSmiles) configData.ligand_smiles = ligandSmiles;

        // Bury ligand setting from AI parser (default: true)
        configData.bury_ligand = intentResult?.data?.bury_ligand ?? true;

        // De novo mode: set design-type-aware contig when no PDB is available
        const resolveResult = Object.values(previousResults).find(r => r.data?.source === 'de_novo');
        if (resolveResult || !pdbContent) {
          configData.de_novo = true;

          if (!configData.contig) {
            // Use AI-recommended chain lengths if available
            if (intentResult?.data?.ai_parsed && intentResult.data.chain_length_min && intentResult.data.chain_length_max) {
              const min = intentResult.data.chain_length_min as number;
              const max = intentResult.data.chain_length_max as number;
              configData.contig = `${min}-${max}`;
            } else if (designType === 'metal' && ligandName) {
              // Metal + ligand binding pocket: ~80 residues, enough for a 4-helix bundle with binding site
              configData.contig = '70-100';
            } else if (designType === 'metal') {
              // Metal-only binding: compact fold
              configData.contig = '50-80';
            } else if (designType === 'binder') {
              // Binder: medium-sized protein
              configData.contig = '80-120';
            } else {
              // General de novo
              configData.contig = '60-100';
            }
          }

          if (targetMetal) {
            configData.target_metal = targetMetal;
            // RFD3 ligand field: specify metal element for coordination
            if (!configData.ligand) {
              configData.ligand = targetMetal;
            }
          }
        }

        ctx.onProgress(100, 'Configuration complete');

        return {
          id: `config-${Date.now()}`,
          summary: `Configured for ${designType} design${targetMetal ? ` (${targetMetal})` : ''}${ligandName ? ` + ${ligandName}` : ''}`,
          data: configData,
        };
      },
    },

    // Step 5: Backbone generation — single-run for metal, raw RFD3 for general
    {
      id: 'rfd3_nl',
      name: 'Backbone Generation',
      description: 'Generate backbone structures with RFD3 (metal-aware conditioning when applicable)',
      icon: Building2,
      requiresReview: true,
      supportsSelection: true,
      optional: false,
      ResultPreview: SweepResultPreview,
      defaultParams: {
        num_designs: 4,
        num_timesteps: 200,
      },
      parameterSchema: [
        { id: 'num_designs', label: 'Number of Designs', type: 'slider', required: false, defaultValue: 4, range: { min: 1, max: 10, step: 1 }, helpText: 'Total designs to generate' },
        { id: 'num_timesteps', label: 'Timesteps', type: 'slider', required: false, defaultValue: 200, range: { min: 50, max: 500, step: 50 } },
        { id: 'step_scale', label: 'Step Scale', type: 'slider', required: false, defaultValue: 1.5, range: { min: 0.5, max: 3, step: 0.1 }, helpText: 'Higher = more designable, less diverse' },
        { id: 'gamma_0', label: 'Gamma', type: 'slider', required: false, defaultValue: 0.6, range: { min: 0.1, max: 1, step: 0.05 }, helpText: 'Lower = more designable' },
        { id: 'use_sweep', label: 'Parameter Sweep', type: 'boolean', required: false, defaultValue: false, helpText: 'Run multi-config sweep (trial-and-error) instead of educated single-run' },
      ],

      async execute(ctx) {
        const { previousResults, params } = ctx;

        // Merge data from ALL previous steps (later steps override earlier)
        let mergedConfig: Record<string, unknown> = {};
        for (const result of Object.values(previousResults)) {
          if (result.data && typeof result.data === 'object') {
            mergedConfig = { ...mergedConfig, ...result.data as Record<string, unknown> };
          }
        }

        const structureResult = Object.values(previousResults).find(r => r.pdbOutputs && r.pdbOutputs.length > 0);
        const pdbContent = structureResult?.pdbOutputs?.[0]?.pdbContent;

        const designType = mergedConfig.design_type as string || 'general';
        const targetMetal = mergedConfig.target_metal as string || '';
        const ligandName = mergedConfig.ligand_name as string || '';

        // --- Metal path: single-run (default) or sweep (opt-in) ---
        if (designType === 'metal' && targetMetal) {
          // Check if scaffold search found a PDB to use as motif
          const scaffoldResult = Object.values(previousResults).find(
            r => r.data?.searched && r.data?.recommended_action === 'scaffold'
          );
          let scaffoldPdb: string | undefined;
          let scaffoldId: string | undefined;
          if (scaffoldResult?.pdbOutputs) {
            const selected = ctx.selectedItems.length > 0
              ? scaffoldResult.pdbOutputs.find(p => ctx.selectedItems.includes(p.id))
              : undefined;
            const chosen = selected || scaffoldResult.pdbOutputs[0];
            scaffoldPdb = chosen?.pdbContent;
            scaffoldId = chosen?.id?.replace('scaffold-', '') || (scaffoldResult.data?.scaffold_pdb_id as string);
          }

          // Map ligand_name to backend ligand code
          const ligandCodeMap: Record<string, string> = {
            citrate: 'CIT', pqq: 'PQQ', atp: 'ATP', nad: 'NAD', heme: 'HEM',
          };
          const ligandCode = ligandCodeMap[ligandName.toLowerCase()] || undefined;
          const motifPdb = scaffoldPdb || pdbContent || '';
          const numDesigns = mergedConfig.num_designs as number ?? params.num_designs as number ?? 4;
          const useSweep = params.use_sweep === true || mergedConfig.use_sweep === true;

          if (useSweep) {
            // --- SWEEP MODE (opt-in): trial-and-error across multiple configs ---
            ctx.onProgress(5, `Starting metal binding sweep (${targetMetal}${ligandName ? ` + ${ligandName}` : ''})...`);

            const allSizes = ['small', 'medium'];
            const allCfgScales = [2.5, 3.0, 3.5];
            const totalConfigs = allSizes.length * allCfgScales.length;

            let sizes = allSizes;
            let cfgScales = allCfgScales;
            let designsPerConfig = 1;

            if (numDesigns <= 3) {
              sizes = ['medium'];
              cfgScales = allCfgScales.slice(0, Math.max(1, numDesigns));
            } else if (numDesigns < totalConfigs) {
              sizes = allSizes;
              cfgScales = allCfgScales.slice(0, Math.ceil(numDesigns / sizes.length));
            } else {
              designsPerConfig = Math.ceil(numDesigns / totalConfigs);
            }

            console.log(`[rfd3_nl] SWEEP: ${sizes.length} sizes × ${cfgScales.length} cfgs × ${designsPerConfig}/config`);

            const sweepResult = await api.startMetalBindingSweep({
              motif_pdb: motifPdb,
              metal: targetMetal,
              ligand: ligandCode,
              sizes,
              cfg_scales: cfgScales,
              designs_per_config: designsPerConfig,
            });

            ctx.onProgress(90, 'Processing sweep results...');

            if (sweepResult.status === 'failed' || !sweepResult.result) {
              throw new Error((sweepResult as any).error || 'Metal binding sweep failed');
            }

            const sr = sweepResult.result;
            const pdbOutputs: PdbOutput[] = [];
            const sequences: SequenceOutput[] = [];

            if (sr.results && Array.isArray(sr.results)) {
              for (const design of sr.results) {
                const pdb = (design as any).pdb_content || (design as any).content;
                if (pdb) {
                  pdbOutputs.push({
                    id: `sweep-${design.name}`,
                    label: `${design.name} (${design.tier})`,
                    pdbContent: pdb,
                    metrics: {
                      pLDDT: design.plddt, pTM: design.ptm, pAE: design.pae,
                      tier: design.tier, config: design.config_name,
                    },
                  });
                }
                if (design.status !== 'fail' && design.sequence) {
                  sequences.push({
                    id: `sweep-${design.name}`,
                    sequence: design.sequence,
                    score: design.plddt,
                    label: `${design.name} (${design.tier}) pLDDT=${design.plddt.toFixed(2)} pTM=${design.ptm.toFixed(2)}`,
                  });
                }
              }
            }

            ctx.onProgress(100, 'Sweep complete');

            const reviewCount = sr.total_review ?? 0;
            const summaryParts = [`Sweep: ${sr.total_generated} designs`];
            if (sr.total_passing > 0) summaryParts.push(`${sr.total_passing} passing`);
            if (reviewCount > 0) summaryParts.push(`${reviewCount} review`);
            if (sr.best_design) summaryParts.push(`best: ${sr.best_design.tier} tier`);

            const sweepAllFailed = sr.total_passing === 0 && reviewCount === 0;
            if (sweepAllFailed) {
              const allDesigns = (sr.results || []) as Array<{ ptm: number }>;
              const bestPtm = allDesigns.length > 0 ? Math.max(...allDesigns.map(d => d.ptm)) : 0;
              summaryParts.length = 0;
              summaryParts.push(`All ${sr.total_generated} designs failed`, `best pTM=${bestPtm.toFixed(3)}`);
            }

            const resultData: Record<string, unknown> = {
              session_id: sr.session_id,
              total_generated: sr.total_generated,
              total_passing: sr.total_passing,
              total_review: sr.total_review,
              pass_rate: sr.pass_rate,
              best_design: sr.best_design,
              config_rankings: sr.config_rankings,
              config_errors: sr.config_errors || null,
            };

            if (sweepAllFailed) {
              resultData.sweep_all_failed = true;
              if (sr.results && Array.isArray(sr.results)) {
                resultData.failed_designs = (sr.results as Array<Record<string, unknown>>).map(d => ({
                  name: d.name, config_name: d.config_name,
                  plddt: d.plddt, ptm: d.ptm, pae: d.pae,
                  tier: d.tier, status: d.status,
                  pdb_content: d.pdb_content || d.content,
                }));
              }
            }

            return {
              id: `nl-metal-${Date.now()}`,
              summary: summaryParts.join(', '),
              pdbOutputs,
              sequences,
              data: resultData,
            };
          }

          // --- SINGLE MODE (default): educated RFD3 backbone generation ---
          // Returns backbone PDBs only; MPNN and RF3 run as separate pipeline steps
          if (scaffoldPdb && scaffoldId) {
            ctx.onProgress(5, `Using scaffold ${scaffoldId} for ${targetMetal}${ligandName ? ` + ${ligandName}` : ''} backbone generation...`);
          } else {
            ctx.onProgress(5, `Generating ${targetMetal}${ligandName ? ` + ${ligandName}` : ''} backbones...`);
          }

          const contig = mergedConfig.contig as string || '110-130';
          const cfgScale = params.step_scale !== undefined
            ? undefined  // don't double-specify
            : (mergedConfig.cfg_scale as number ?? 2.0);

          console.log(`[rfd3_nl] SINGLE: ${numDesigns} designs, contig=${contig}, cfg=${cfgScale}`);

          const buryLigand = mergedConfig.bury_ligand !== false; // default true
          const singleResult = await api.submitMetalSingleDesign({
            motif_pdb: motifPdb,
            metal: targetMetal,
            ligand: ligandCode,
            contig,
            cfg_scale: cfgScale ?? 2.0,
            num_designs: numDesigns,
            num_timesteps: params.num_timesteps as number ?? mergedConfig.num_timesteps as number ?? 200,
            step_scale: params.step_scale as number ?? mergedConfig.step_scale as number ?? 1.5,
            gamma_0: params.gamma_0 as number ?? mergedConfig.gamma_0 as number ?? undefined,
            bury_ligand: buryLigand,
            seed: 42,
          });

          ctx.onProgress(90, 'Processing backbone results...');

          if (singleResult.status === 'failed' || !singleResult.result) {
            throw new Error((singleResult as any).error || 'Metal backbone generation failed');
          }

          const designs = singleResult.result.designs;
          const pdbOutputs: PdbOutput[] = [];

          for (const [idx, d] of designs.entries()) {
            pdbOutputs.push({
              id: `metal-design-${idx + 1}`,
              label: `Design ${idx + 1}`,
              pdbContent: d.content,
            });
          }

          ctx.onProgress(100, 'Backbone generation complete');

          return {
            id: `nl-metal-${Date.now()}`,
            summary: `Generated ${pdbOutputs.length} backbone(s) for ${targetMetal}${ligandName ? ` + ${ligandName}` : ''}`,
            pdbOutputs,
            data: {
              metal_single_mode: true,
              target_metal: targetMetal,
              ligand_name: ligandName,
              ligand_code: ligandCode,
            },
          };
        }

        // --- General RFD3: raw backbone generation ---
        ctx.onProgress(5, 'Submitting RFD3 design job...');

        const request: Record<string, unknown> = {
          num_designs: params.num_designs ?? mergedConfig.num_designs ?? 4,
          num_timesteps: params.num_timesteps ?? mergedConfig.num_timesteps ?? 200,
          step_scale: params.step_scale ?? mergedConfig.step_scale ?? 1.5,
          gamma_0: params.gamma_0 ?? mergedConfig.gamma_0 ?? 0.6,
        };

        if (pdbContent) request.pdb_content = pdbContent;
        if (mergedConfig.ligand) request.ligand = mergedConfig.ligand;
        if (mergedConfig.unindex) request.unindex = mergedConfig.unindex;
        if (mergedConfig.partial_t) request.partial_t = mergedConfig.partial_t;
        if (mergedConfig.contig) request.contig = mergedConfig.contig;

        const response = await api.submitRFD3Design(request as any);
        ctx.onJobCreated?.(response.job_id, 'rfd3');
        ctx.onProgress(15, 'Generating structures...');

        let result: Record<string, unknown>;

        if (response.syncCompleted) {
          if (response.status !== 'completed' || !response.result) {
            throw new Error(response.error || 'Design failed');
          }
          result = response.result as Record<string, unknown>;
        } else {
          const polled = await api.waitForJob(response.job_id, (s) => {
            if (s.status === 'running') ctx.onProgress(50, 'Generating...');
          });
          if (polled.status !== 'completed') throw new Error(polled.error || 'Design failed');
          result = polled.result as Record<string, unknown>;
        }

        const designs = result.designs as Array<Record<string, unknown>> | undefined;
        const pdbOutputs: PdbOutput[] = [];

        if (Array.isArray(designs)) {
          designs.forEach((d, idx) => {
            const pdb = (d.content as string) || (d.pdb_content as string);
            if (pdb) {
              pdbOutputs.push({
                id: `nl-design-${idx + 1}`,
                label: `Design ${idx + 1}`,
                pdbContent: pdb,
              });
            }
          });
        }

        if (pdbOutputs.length === 0) {
          const pdb = extractPdbContent(result);
          if (pdb) pdbOutputs.push({ id: 'nl-design-1', label: 'Design 1', pdbContent: pdb });
        }

        ctx.onProgress(100, 'Structure generation complete');

        return {
          id: `nl-rfd3-${Date.now()}`,
          summary: `Generated ${pdbOutputs.length} backbone design(s)`,
          pdbOutputs,
        };
      },
    },

    // Step 5.5: Scout filter (optional — validates 1 seq per backbone and filters)
    createScoutFilterStep({ id: 'scout_filter_nl' }),

    // Step 6: MPNN sequence design
    createMpnnStep({ id: 'mpnn_nl' }),

    // Step 7: RF3 validation
    createRf3Step({ id: 'rf3_nl' }),

    // Step 8: Final analysis
    {
      id: 'analysis',
      name: 'Final Analysis',
      description: 'Summarize and rank all results',
      icon: BarChart3,
      requiresReview: true,
      supportsSelection: false,
      optional: false,
      defaultParams: {},
      parameterSchema: [],

      async execute(ctx) {
        const { previousResults, selectedItems, initialParams } = ctx;

        const allPdbs: PdbOutput[] = [];
        for (const r of Object.values(previousResults)) {
          if (r.pdbOutputs) allPdbs.push(...r.pdbOutputs);
        }
        // Filter to only validated structures (from RF3 step)
        const validated = allPdbs.filter(p => p.id.startsWith('rf3-') || p.id.startsWith('validated-'));
        const toAnalyze = validated.length > 0 ? validated : allPdbs.slice(-8);

        ctx.onProgress(30, 'Computing final metrics...');

        // Try metal evaluation if applicable
        const configResult = Object.values(previousResults).find(r => r.data?.design_type);
        const designType = configResult?.data?.design_type as string;
        const targetMetal = configResult?.data?.target_metal as string || initialParams.target_metal as string;

        const analyzed: PdbOutput[] = [];

        for (let i = 0; i < toAnalyze.length; i++) {
          const pdb = toAnalyze[i];
          ctx.onProgress(30 + (i / toAnalyze.length) * 60, `Analyzing ${pdb.label}...`);

          const metrics = { ...pdb.metrics };

          if (designType === 'metal' && targetMetal) {
            try {
              const evaluation = await api.evaluateDesign({
                pdb_content: pdb.pdbContent,
                target_metal: targetMetal,
              });
              metrics.coordination = evaluation.coordination_number;
              metrics.geometry = evaluation.geometry_type;
              metrics.pass = evaluation.overall_pass ? 'Yes' : 'No';
            } catch {
              // Non-fatal
            }
          }

          analyzed.push({ ...pdb, metrics });
        }

        // Sort by pLDDT if available
        analyzed.sort((a, b) => {
          const aScore = (a.metrics?.pLDDT as number) ?? 0;
          const bScore = (b.metrics?.pLDDT as number) ?? 0;
          return bScore - aScore;
        });

        ctx.onProgress(100, 'Analysis complete');

        return {
          id: `analysis-${Date.now()}`,
          summary: `Final analysis: ${analyzed.length} designs evaluated`,
          pdbOutputs: analyzed,
          data: {
            total_analyzed: analyzed.length,
            design_type: designType,
            top_score: analyzed[0]?.metrics?.pLDDT,
          },
        };
      },
    },

    // Step 9: Save to design history (automatic, no user interaction)
    createSaveHistoryStep({ id: 'save_history_nl' }),

    // Step 10: Check for lesson triggers (shows results if detected)
    createCheckLessonsStep({ id: 'check_lessons_nl' }),
  ],
};
