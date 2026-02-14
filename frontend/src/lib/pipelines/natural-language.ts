/**
 * Natural Language Pipeline
 * Steps: parse_intent → resolve_structure → ligand_features → scaffold_search → configure → structure_gen → mpnn_sequence → validation → analysis
 *
 * This pipeline handles free-form natural language requests by parsing user intent,
 * resolving the target structure, and configuring the appropriate design parameters.
 */

import { MessageSquare, Search, Settings, Building2, Shield, Dna, BarChart3, FlaskConical, Atom } from 'lucide-react';
import api from '@/lib/api';
import { extractPdbContent } from '@/lib/workflowHandlers';
import type { PipelineDefinition, PdbOutput, SequenceOutput } from '@/lib/pipeline-types';
import type { LigandCoordinationFeature, ProteinDonor } from '@/lib/store';
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
import { CoordinationAnalysisPreview } from '@/components/pipeline/CoordinationAnalysisPreview';

// ============================================================================
// AI Parse Helpers
// ============================================================================

/** Fallback SMILES table (used when backend is unreachable) */
const LIGAND_SMILES_FALLBACK: Record<string, string> = {
  citrate: 'OC(=O)CC(O)(CC(O)=O)C(O)=O',
  pqq: 'OC(=O)c1[nH]c2c(c1)C(=O)C(=O)c3nc(cc(C(O)=O)c23)C(O)=O',
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

/** Resolve ligand SMILES via backend API, with sync fallback */
async function resolveLigandSmiles(
  name?: string,
  metalType?: string,
  isomerSpec?: string,
): Promise<{ smiles: string; source: string }> {
  if (!name) return { smiles: '', source: 'none' };

  try {
    const result = await api.resolveLigand({
      ligand_name: name,
      metal_type: metalType,
      isomer_spec: isomerSpec,
    });
    if (result.success && result.smiles) {
      return { smiles: result.smiles, source: result.source };
    }
  } catch (err) {
    console.warn('[resolveLigandSmiles] Backend resolution failed, using fallback:', err);
  }

  // Sync fallback for offline/error
  const fallback = LIGAND_SMILES_FALLBACK[name.toLowerCase()] || '';
  return { smiles: fallback, source: fallback ? 'fallback' : 'none' };
}

/** Sync-only lookup for immediate use (no network) */
function lookupSmilesFallback(name?: string): string {
  if (!name) return '';
  return LIGAND_SMILES_FALLBACK[name.toLowerCase()] || '';
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

  // 2) Detect ligands — known ligands with SMILES, then dynamic detection for unknowns
  const ligandMap: Record<string, { name: string; smiles: string }> = {
    'pyrroloquinoline quinone': { name: 'PQQ', smiles: 'OC(=O)c1[nH]c2c(c1)C(=O)C(=O)c3nc(cc(C(O)=O)c23)C(O)=O' },
    pyrroloquinoline: { name: 'PQQ', smiles: 'OC(=O)c1[nH]c2c(c1)C(=O)C(=O)c3nc(cc(C(O)=O)c23)C(O)=O' },
    'citric acid': { name: 'citrate', smiles: 'OC(=O)CC(O)(CC(O)=O)C(O)=O' },
    citrate:  { name: 'citrate',  smiles: 'OC(=O)CC(O)(CC(O)=O)C(O)=O' },
    citric:   { name: 'citrate',  smiles: 'OC(=O)CC(O)(CC(O)=O)C(O)=O' },
    pqq:      { name: 'PQQ',      smiles: 'OC(=O)c1[nH]c2c(c1)C(=O)C(=O)c3nc(cc(C(O)=O)c23)C(O)=O' },
    atp:      { name: 'ATP',      smiles: 'c1nc(c2c(n1)n(cn2)[C@@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N' },
    adp:      { name: 'ADP',      smiles: '' },
    gtp:      { name: 'GTP',      smiles: '' },
    nad:      { name: 'NAD+',     smiles: '' },
    nadh:     { name: 'NAD+',     smiles: '' },
    nadp:     { name: 'NADP',     smiles: '' },
    fad:      { name: 'FAD',      smiles: '' },
    fmn:      { name: 'FMN',      smiles: '' },
    glucose:  { name: 'glucose',  smiles: 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O' },
    azobenzene: { name: 'azobenzene', smiles: 'c1ccc(/N=N/c2ccccc2)cc1' },
    heme:     { name: 'heme',     smiles: '' },
    haem:     { name: 'heme',     smiles: '' },
  };

  // Check known ligands first (longest match first via sorted keys)
  const sortedLigandKeys = Object.keys(ligandMap).sort((a, b) => b.length - a.length);
  for (const keyword of sortedLigandKeys) {
    // Use word boundary matching for short keywords to avoid false matches
    if (keyword.length <= 3) {
      const regex = new RegExp(`\\b${keyword}\\b`, 'i');
      if (regex.test(promptLower)) {
        ligandName = ligandMap[keyword].name;
        ligandSmiles = ligandMap[keyword].smiles;
        break;
      }
    } else if (promptLower.includes(keyword)) {
      ligandName = ligandMap[keyword].name;
      ligandSmiles = ligandMap[keyword].smiles;
      break;
    }
  }

  // Dynamic ligand detection: catch unknown ligand-like terms (2-5 uppercase letters
  // or short words near "bind" context) that aren't metals or common English words
  if (!ligandName) {
    const allMetals = new Set([
      ...Object.values(metalNameMap),
      ...Object.values(metalSymbolMap),
      ...Object.keys(metalNameMap),
      ...Object.keys(metalSymbolMap),
    ]);
    const stopWords = new Set([
      'the', 'and', 'for', 'with', 'that', 'this', 'from', 'bind', 'binding',
      'design', 'create', 'make', 'want', 'protein', 'metal', 'ligand',
      'a', 'an', 'to', 'of', 'in', 'on', 'is', 'it', 'not', 'my', 'i',
    ]);

    // Look for 2-5 letter uppercase abbreviation-like words (PQQ, ATP, NAD, etc.)
    const abbreviationMatch = prompt.match(/\b([A-Z][A-Z0-9]{1,4})\b/g);
    if (abbreviationMatch) {
      for (const abbr of abbreviationMatch) {
        const abbrUpper = abbr.toUpperCase();
        if (!allMetals.has(abbrUpper) && !allMetals.has(abbr.toLowerCase())) {
          ligandName = abbrUpper;
          break;
        }
      }
    }

    // Also look for unrecognized words adjacent to "bind" that could be ligand names
    if (!ligandName) {
      const bindContext = promptLower.match(/bind(?:ing|s)?\s+(\w+)/);
      if (bindContext) {
        const candidate = bindContext[1];
        if (
          !stopWords.has(candidate) &&
          !allMetals.has(candidate) &&
          !allMetals.has(candidate.toUpperCase()) &&
          !Object.keys(metalNameMap).includes(candidate)
        ) {
          ligandName = candidate.toUpperCase();
        }
      }
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
  name: 'Protein Design Studio',
  description: 'Describe what you want to build and we\'ll design it',
  icon: MessageSquare,
  requiresBackend: true,

  initialParams: [
    { id: 'user_prompt', label: 'Design Request', type: 'text', required: true, helpText: 'What protein would you like to design?' },
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
      name: 'Understanding Your Design',
      description: 'Reading your request and figuring out what to build',
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

          // Resolve ligand SMILES via backend (async, with isomer support)
          const isomerSpec = ai.isomer_specification;
          const ligandResolution = await resolveLigandSmiles(
            ai.ligand_name,
            ai.metal_type,
            isomerSpec,
          );

          intentData = {
            design_type: mapDesignType(ai),
            target_metal: ai.metal_type || '',
            ligand_name: ai.ligand_name || '',
            ligand_smiles: ligandResolution.smiles,
            ligand_smiles_source: ligandResolution.source,
            isomer_specification: isomerSpec || null,
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
          const fallback = keywordFallbackParse(prompt, pdbId);
          // Try async ligand resolution even for fallback path
          const fallbackLigand = fallback.ligand_name as string || '';
          if (fallbackLigand) {
            try {
              const resolution = await resolveLigandSmiles(fallbackLigand, fallback.target_metal as string);
              fallback.ligand_smiles = resolution.smiles;
              fallback.ligand_smiles_source = resolution.source;
            } catch {
              // Keep existing fallback smiles
            }
          }
          intentData = fallback;
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
      name: 'Finding Starting Point',
      description: 'Looking up or loading your starting structure',
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
    // MOVED UP: Now runs before coordination analysis so scaffold data is available
    createScaffoldSearchStep({
      id: 'scaffold_search_nl',
      name: 'Searching Nature\'s Library',
      description: 'Looking for existing structures in the Protein Data Bank that match your design',
    }),

    // Step 4: Coordination Analysis — combined protein + ligand donor identification
    {
      id: 'coordination_analysis_nl',
      name: 'Analyzing Binding Sites',
      description: 'Identifying how your metal and ligand interact at the molecular level',
      icon: Atom,
      requiresReview: true,
      supportsSelection: false,
      optional: true,
      defaultParams: {},
      parameterSchema: [],
      ResultPreview: CoordinationAnalysisPreview,

      async execute(ctx) {
        const { previousResults, selectedItems } = ctx;

        // Get intent data
        const intentResult = Object.values(previousResults).find(r => r.data?.design_type);
        const ligandName = intentResult?.data?.ligand_name as string || '';
        const ligandSmiles = intentResult?.data?.ligand_smiles as string || '';
        const targetMetal = intentResult?.data?.target_metal as string || '';

        // Skip if no metal/ligand
        if (!targetMetal && !ligandName) {
          return {
            id: `coordination-${Date.now()}`,
            summary: 'Skipped: no metal-ligand specified',
            data: { skipped: true, reason: 'No metal or ligand in design intent' },
          };
        }

        // Get scaffold results (now available from prior step)
        const scaffoldResult = previousResults['scaffold_search_nl'];
        const candidates = scaffoldResult?.data?.candidates as Array<Record<string, unknown>> | undefined;

        // Check if user selected a specific scaffold (selectedItems from context: ['scaffold-4CVB'])
        const scaffoldSelections = selectedItems.filter(id => id.startsWith('scaffold-'));
        let selectedCandidate: Record<string, unknown> | undefined;
        let scaffoldPdb: string | undefined;

        if (scaffoldSelections.length > 0 && candidates) {
          // User selected a specific scaffold - use that one
          const selectedPdbId = scaffoldSelections[0].replace('scaffold-', '');
          selectedCandidate = candidates.find(c => c.pdb_id === selectedPdbId);
          const pdbOutput = scaffoldResult?.pdbOutputs?.find(
            p => p.id === scaffoldSelections[0]
          );
          scaffoldPdb = pdbOutput?.pdbContent;
          console.log('[coordination_analysis] Using user-selected scaffold:', selectedPdbId);
        } else {
          // No selection - use best_candidate (first/default)
          selectedCandidate = scaffoldResult?.data?.best_candidate as Record<string, unknown> | undefined;
          scaffoldPdb = scaffoldResult?.pdbOutputs?.[0]?.pdbContent;
        }

        // Get the actual PDB ligand code from scaffold (e.g., "CIT" for citrate)
        // This is critical for Molstar selection - it needs the PDB residue code, not user's query string
        const scaffoldLigandCode = (selectedCandidate?.ligand_code as string) ||
                                    (scaffoldResult?.data?.ligand_code as string) || '';

        console.log('[coordination_analysis] scaffoldResult:', scaffoldResult ? 'present' : 'missing');
        console.log('[coordination_analysis] selectedCandidate:', selectedCandidate?.pdb_id);
        console.log('[coordination_analysis] scaffoldPdb available:', scaffoldPdb ? 'yes' : 'no');
        console.log('[coordination_analysis] scaffoldLigandCode:', scaffoldLigandCode);
        console.log('[coordination_analysis] protein_donors from scaffold:', selectedCandidate?.protein_donors);
        console.log('[coordination_analysis] ligand_donors from scaffold:', selectedCandidate?.ligand_donors);

        ctx.onProgress(20, 'Analyzing ligand features...');

        // Part 1: Get ligand features from backend
        let ligandFeatures: Record<string, unknown> | null = null;
        if (ligandName || ligandSmiles) {
          try {
            ligandFeatures = await api.analyzeLigandFeatures({
              ligand_name: ligandName,
              smiles: ligandSmiles || undefined,
              metal_type: targetMetal || undefined,
            }) as Record<string, unknown>;
          } catch (err) {
            console.warn('[coordination_analysis] Ligand features failed:', err);
          }
        }

        ctx.onProgress(60, 'Merging coordination data...');

        // Part 2: Parse protein donors from scaffold search result
        const proteinDonors: ProteinDonor[] = [];
        const scaffoldProteinDonors = selectedCandidate?.protein_donors as string[] | undefined;

        console.log('[coordination_analysis] Raw protein_donors:', scaffoldProteinDonors);

        if (scaffoldProteinDonors && Array.isArray(scaffoldProteinDonors)) {
          for (const donorStr of scaffoldProteinDonors) {
            // Try multiple patterns to handle format variations
            // Pattern 1: "A26:SER O@2.29A" or "A26:SER OG@2.29"
            let match = donorStr.match(/([A-Z])(\d+):(\w+)\s+(\w+)@([\d.]+)/);

            // Pattern 2: "A:26:SER:O@2.29" (colon-separated)
            if (!match) {
              match = donorStr.match(/([A-Z]):(\d+):(\w+):(\w+)@([\d.]+)/);
            }

            // Pattern 3: More lenient - just grab key parts
            if (!match) {
              match = donorStr.match(/([A-Z]+)[\s:]+(\d+)[\s:]+(\w{3})\s+(\w+)@([\d.]+)/);
            }

            if (match) {
              proteinDonors.push({
                chain: match[1],
                residue: parseInt(match[2], 10),
                name: match[3],
                atomName: match[4],
                distance: parseFloat(match[5]),
                enabled: true,
                source: 'scaffold',
              });
              console.log('[coordination_analysis] Parsed donor:', match[0]);
            } else {
              console.warn('[coordination_analysis] Failed to parse donor string:', donorStr);
            }
          }
        }

        console.log('[coordination_analysis] Final proteinDonors:', proteinDonors.length);

        // Part 3: Build ligand donors from scaffold data (primary source for visualization)
        // and enrich with KB features (for HSAB, type info)
        const ligandDonors: LigandCoordinationFeature[] = [];

        // Use scaffold ligand code for Molstar selection (PDB residue code like "CIT"),
        // falling back to user's ligand name if no scaffold data
        const pdbLigandCode = scaffoldLigandCode || ligandName.toUpperCase().slice(0, 3);

        // Parse scaffold ligand donors - these have the ACTUAL atom names from PDB
        const scaffoldLigandDonorsRaw = selectedCandidate?.ligand_donors as string[] | undefined;
        const scaffoldLigandAtoms: Array<{ atomName: string; distance: number }> = [];
        if (scaffoldLigandDonorsRaw) {
          for (const donorStr of scaffoldLigandDonorsRaw) {
            // Format: "O1@2.29" or "O1@2.29Å" or "O7@2.15"
            const match = donorStr.match(/(\w+)@([\d.]+)/);
            if (match) {
              scaffoldLigandAtoms.push({
                atomName: match[1],
                distance: parseFloat(match[2]),
              });
            }
          }
        }

        // Get KB features for enrichment (HSAB, type info)
        const features = ligandFeatures?.features as Array<Record<string, unknown>> | undefined;
        const kbByElement: Record<string, Record<string, unknown>> = {};
        if (features) {
          for (const f of features) {
            if (f.is_coordination_donor && f.element) {
              const elem = (f.element as string).toUpperCase();
              // Store first KB entry for each element (for HSAB/type lookup)
              if (!kbByElement[elem]) {
                kbByElement[elem] = f;
              }
            }
          }
        }

        console.log('[coordination_analysis] scaffoldLigandAtoms:', scaffoldLigandAtoms);
        console.log('[coordination_analysis] kbByElement:', Object.keys(kbByElement));

        // Build ligand donors from scaffold atoms (prioritized for correct visualization)
        if (scaffoldLigandAtoms.length > 0) {
          for (const { atomName, distance } of scaffoldLigandAtoms) {
            // Extract element from atom name (e.g., "O1" -> "O", "N2" -> "N")
            const element = atomName.replace(/\d+/g, '').toUpperCase();
            // Look up KB data for this element
            const kbData = kbByElement[element];

            ligandDonors.push({
              atom_idx: 0,
              atom_name: atomName, // Use actual PDB atom name
              element,
              type: (kbData?.type as string) || 'donor',
              is_coordination_donor: true,
              coords: null,
              hsab: (kbData?.hsab as string) || null,
              enabled: true,
              ligand_name: pdbLigandCode, // Use PDB residue code for Molstar selection
              distance,
            });
          }
        } else if (features) {
          // No scaffold data - fall back to KB features (won't work for visualization
          // but at least shows something in the UI)
          for (const f of features) {
            if (f.is_coordination_donor) {
              const atomName = f.atom_name as string;
              ligandDonors.push({
                atom_idx: f.atom_idx as number || 0,
                atom_name: atomName,
                element: f.element as string || '',
                type: f.type as string || 'donor',
                is_coordination_donor: true,
                coords: f.coords as [number, number, number] | null,
                hsab: f.hsab as string | null,
                enabled: true,
                ligand_name: pdbLigandCode, // Use PDB code for Molstar selection
              });
            }
          }
        }

        // Build PDB outputs (pass scaffold PDB + ligand PDB)
        const pdbOutputs: PdbOutput[] = [];
        if (scaffoldPdb) {
          pdbOutputs.push({
            id: 'coordination-scaffold',
            label: (selectedCandidate?.pdb_id as string) || 'Scaffold',
            pdbContent: scaffoldPdb,
          });
        }
        const ligandPdbContent = ligandFeatures?.ligand_pdb_content as string | undefined;
        if (ligandPdbContent) {
          pdbOutputs.push({
            id: 'coordination-ligand',
            label: ligandName || 'Ligand',
            pdbContent: ligandPdbContent,
          });
        }

        ctx.onProgress(100, 'Coordination analysis complete');

        const coordNum = (selectedCandidate?.coordination_number as number) || ligandDonors.length;
        const summary = proteinDonors.length > 0
          ? `${targetMetal}: ${proteinDonors.length} protein + ${ligandDonors.length} ligand donors`
          : `${ligandName}: ${ligandDonors.length} coordination donors`;

        return {
          id: `coordination-${Date.now()}`,
          summary,
          pdbOutputs,
          data: {
            metal: targetMetal,
            ligand_name: ligandName,
            coordination_number: coordNum,
            scaffold_pdb_id: selectedCandidate?.pdb_id,
            protein_donors: proteinDonors,
            ligand_donors: ligandDonors,
            ligand_features_raw: ligandFeatures,
            has_scaffold_data: !!selectedCandidate,
            has_ligand_kb_data: ligandFeatures?.source === 'knowledge_base',
          },
        };
      },
    },

    // Step 5: Configure design parameters
    {
      id: 'configure',
      name: 'Tuning the Recipe',
      description: 'Setting up the design parameters based on your requirements',
      icon: Settings,
      requiresReview: true,
      supportsSelection: false,
      optional: false,
      defaultParams: {},
      parameterSchema: [
        { id: 'num_designs', label: 'Number of Designs', type: 'slider', required: false, defaultValue: 4, range: { min: 1, max: 200, step: 1 }, helpText: 'More designs = better coverage but longer runtime' },
        { id: 'num_timesteps', label: 'Design Quality', type: 'slider', required: false, defaultValue: 200, range: { min: 50, max: 500, step: 50 }, helpText: 'More steps = finer-grained structures (default is good for most designs)' },
        {
          id: 'design_goal', label: 'Design Goal', type: 'select' as const, required: false, defaultValue: 'binding',
          options: [
            { value: 'binding', label: 'Tight binding' },
            { value: 'catalysis', label: 'Catalytic activity' },
            { value: 'sensing', label: 'Sensing / detection' },
            { value: 'structural', label: 'Structural role' },
          ],
          helpText: 'Shapes how the binding site is built — enclosed for binding, open for catalysis or sensing',
        },
        {
          id: 'filter_tier', label: 'Quality Bar', type: 'select' as const, required: false, defaultValue: 'standard',
          options: [
            { value: 'relaxed', label: 'Exploratory' },
            { value: 'standard', label: 'Standard' },
            { value: 'strict', label: 'Strict' },
          ],
          helpText: 'Controls which designs pass inspection',
        },
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

        // Design goal from wizard or configure step params → backend params
        const designGoal = params.design_goal as string || initialParams.design_goal as string || '';
        if (designGoal) {
          configData.design_goal = designGoal;
          if (designGoal === 'catalysis') {
            configData.bury_ligand = false;
          } else if (designGoal === 'sensing') {
            configData.bury_ligand = false;
          } else if (designGoal === 'structural') {
            configData.bury_ligand = true;
            configData.stability_focus = true;
          }
          // 'binding' → default (bury_ligand=true)
        }

        // Filter tier from wizard or configure step params → passed to analysis step
        const filterTier = params.filter_tier as string || initialParams.filter_tier as string || '';
        if (filterTier) {
          configData.filter_tier = filterTier;
        }

        // Wire in coordination analysis results (from coordination_analysis_nl step)
        const coordResult = previousResults['coordination_analysis_nl'];
        const coordData = coordResult?.data as Record<string, unknown> | undefined;

        if (coordData && !coordData.skipped) {
          // Check for user overrides
          const userOverrides = coordData.user_overrides as {
            active_protein_donors?: string[];
            active_ligand_donors?: string[];
            modified: boolean;
          } | undefined;

          // Ligand coordination donors
          if (userOverrides?.modified && userOverrides.active_ligand_donors) {
            configData.ligand_coordination_donors = userOverrides.active_ligand_donors;
            configData.recommended_fixing_strategy = 'fix_coordination';
          } else {
            // Use all enabled ligand donors from analysis
            const ligandDonors = coordData.ligand_donors as Array<{ atom_name: string; enabled: boolean }> | undefined;
            if (ligandDonors) {
              const activeLigandDonors = ligandDonors.filter(d => d.enabled).map(d => d.atom_name);
              if (activeLigandDonors.length > 0) {
                configData.ligand_coordination_donors = activeLigandDonors;
              }
            }
          }

          // Protein catalytic residues (from scaffold search)
          if (userOverrides?.modified && userOverrides.active_protein_donors) {
            configData.fixed_catalytic_residues = userOverrides.active_protein_donors;
            configData.recommended_fixing_strategy = 'fix_coordination';
          } else {
            const proteinDonors = coordData.protein_donors as Array<{ chain: string; residue: number; enabled: boolean }> | undefined;
            if (proteinDonors) {
              const activeProtein = proteinDonors.filter(d => d.enabled).map(d => `${d.chain}${d.residue}`);
              if (activeProtein.length > 0) {
                configData.fixed_catalytic_residues = activeProtein;
                configData.recommended_fixing_strategy = 'fix_coordination';
              }
            }
          }
        }

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
            } else if (designType === 'metal' && ligandName && coordData?.recommended_contig_min) {
              // Cofactor-aware contig from backend (scales with ligand heavy atom count)
              const cMin = coordData.recommended_contig_min as number;
              const cMax = coordData.recommended_contig_max as number;
              configData.contig = `${cMin}-${cMax}`;
              console.log(`[configure] Cofactor-aware contig: ${cMin}-${cMax} (from ligand analysis)`);
            } else if (designType === 'metal' && ligandName) {
              // Fallback for metal+ligand when no coordination analysis ran
              configData.contig = '90-120';
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

    // Step 6: Backbone generation — single-run for metal, raw RFD3 for general
    {
      id: 'rfd3_nl',
      name: 'Building Structures',
      description: 'Generating protein backbone structures that match your design',
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
        { id: 'num_designs', label: 'Number of Designs', type: 'slider', required: false, defaultValue: 4, range: { min: 1, max: 200, step: 1 }, helpText: 'More designs = better coverage but longer runtime (per config if sweep enabled)' },
        { id: 'num_timesteps', label: 'Design Quality', type: 'slider', required: false, defaultValue: 200, range: { min: 50, max: 500, step: 50 }, helpText: 'More steps = finer-grained structures (default is good for most designs)' },
        { id: 'step_scale', label: 'Step Scale', type: 'slider', required: false, defaultValue: 1.5, range: { min: 0.5, max: 3, step: 0.1 }, helpText: 'Higher = more realistic folds, lower = more variety' },
        { id: 'gamma_0', label: 'Gamma', type: 'slider', required: false, defaultValue: 0.6, range: { min: 0.1, max: 1, step: 0.05 }, helpText: 'Lower values favor more realistic folds' },
        { id: 'use_sweep', label: 'Parameter Sweep', type: 'boolean', required: false, defaultValue: false, helpText: 'Try multiple parameter combinations instead of a single educated run' },
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
            ctx.onProgress(5, `Generating ${numDesigns} backbone(s) from scaffold ${scaffoldId} for ${targetMetal}${ligandName ? ` + ${ligandName}` : ''}...`);
          } else {
            ctx.onProgress(5, `Generating ${numDesigns} ${targetMetal}${ligandName ? ` + ${ligandName}` : ''} backbone(s)...`);
          }

          const contig = mergedConfig.contig as string || '110-130';
          const cfgScale = params.step_scale !== undefined
            ? undefined  // don't double-specify
            : (mergedConfig.cfg_scale as number ?? 2.0);

          console.log(`[rfd3_nl] SINGLE: ${numDesigns} designs, contig=${contig}, cfg=${cfgScale}`);

          const buryLigand = mergedConfig.bury_ligand !== false; // default true
          ctx.onProgress(10, `Running RFdiffusion3: ${numDesigns} design(s), contig=${contig}...`);

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
            design_goal: mergedConfig.design_goal as string || 'binding',
            seed: 42,
          });

          ctx.onProgress(90, `Processing ${numDesigns} backbone result(s)...`);

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
        // Priority: mergedConfig (from configure step) → step params → default
        const numDesigns = (mergedConfig.num_designs ?? params.num_designs ?? 4) as number;
        const numTimesteps = (mergedConfig.num_timesteps ?? params.num_timesteps ?? 200) as number;
        ctx.onProgress(5, `Submitting RFD3 job: ${numDesigns} design(s), ${numTimesteps} timesteps...`);

        const request: Record<string, unknown> = {
          num_designs: numDesigns,
          num_timesteps: numTimesteps,
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
        ctx.onProgress(15, `Generating ${numDesigns} backbone(s) via RFdiffusion3...`);

        let result: Record<string, unknown>;

        if (response.syncCompleted) {
          if (response.status !== 'completed' || !response.result) {
            throw new Error(response.error || 'Design failed');
          }
          result = response.result as Record<string, unknown>;
        } else {
          let pollCount = 0;
          const polled = await api.waitForJob(response.job_id, (s) => {
            if (s.status === 'running') {
              pollCount++;
              // Simulate progress during generation (15% → 90%)
              const estimatedProgress = Math.min(90, 15 + pollCount * 5);
              ctx.onProgress(estimatedProgress, `Generating ${numDesigns} backbone(s)... (running)`);
            }
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

    // Step 6.5: Scout filter (optional — validates 1 seq per backbone and filters)
    createScoutFilterStep({
      id: 'scout_filter_nl',
      name: 'Finding the Strongest',
      description: 'Quick-testing each backbone to keep only the most promising ones',
    }),

    // Step 7: MPNN sequence design
    createMpnnStep({
      id: 'mpnn_nl',
      name: 'Writing Genetic Code',
      description: 'Designing amino acid sequences that fold into your structures',
    }),

    // Step 8: RF3 validation
    createRf3Step({
      id: 'rf3_nl',
      name: 'Testing Your Designs',
      description: 'Predicting how each sequence folds to check if the design works',
    }),

    // Step 9: Final analysis — UnifiedDesignAnalyzer + FILTER_PRESETS
    {
      id: 'analysis',
      name: 'Evaluating Results',
      description: 'Scoring each design and filtering for the best candidates',
      icon: BarChart3,
      requiresReview: true,
      supportsSelection: false,
      optional: false,
      defaultParams: {
        filter_tier: 'standard',
      },
      parameterSchema: [
        {
          id: 'filter_tier',
          label: 'Quality Standard',
          type: 'select' as const,
          required: false,
          defaultValue: 'standard',
          options: [
            { value: 'relaxed', label: 'Exploratory (cast a wide net)' },
            { value: 'standard', label: 'Standard' },
            { value: 'strict', label: 'High confidence only' },
          ],
          helpText: 'Controls which designs pass inspection. Stricter = fewer but higher-confidence results.',
        },
      ],

      async execute(ctx) {
        const { previousResults, initialParams, params } = ctx;

        const allPdbs: PdbOutput[] = [];
        for (const r of Object.values(previousResults)) {
          if (r.pdbOutputs) allPdbs.push(...r.pdbOutputs);
        }
        // Filter to only validated structures (from RF3 step)
        const validated = allPdbs.filter(p => p.id.startsWith('rf3-') || p.id.startsWith('validated-'));
        const toAnalyze = validated.length > 0 ? validated : allPdbs.slice(-8);

        // Find a backbone PDB with HETATM records (for metal/ligand coordination analysis).
        // RF3 outputs protein-only PDB, so we need the original backbone PDB to detect metal sites.
        const backbonePdbs = allPdbs.filter(p => p.id.startsWith('rfd3-') || p.id.startsWith('backbone-'));
        const backbonePdb = backbonePdbs.find(p => p.pdbContent?.includes('HETATM'))?.pdbContent;

        // Gather context from previous steps
        const configResult = Object.values(previousResults).find(r => r.data?.design_type);

        // filter_tier: step params (user override at review gate) → configure step data → default
        const filterTier = (params.filter_tier as string)
          || (configResult?.data?.filter_tier as string)
          || 'standard';
        ctx.onProgress(10, `Analyzing ${toAnalyze.length} designs (${filterTier} tier)...`);
        const designType = configResult?.data?.design_type as string || 'general';
        const targetMetal = configResult?.data?.target_metal as string || initialParams.target_metal as string || undefined;
        const ligandName = configResult?.data?.ligand_name as string || undefined;

        const analyzed: PdbOutput[] = [];
        let passCount = 0;
        let filterPreset = 'scout_relaxed';

        for (let i = 0; i < toAnalyze.length; i++) {
          const pdb = toAnalyze[i];
          const pct = 10 + Math.round((i / toAnalyze.length) * 80);
          ctx.onProgress(pct, `Analyzing ${pdb.label} (${i + 1}/${toAnalyze.length})...`);

          const metrics = { ...pdb.metrics };

          try {
            // Forward RF3 ligand-interface confidence metrics for AF3 paper filtering
            // (bioRxiv 2025.09.18.676967v2: iPTM > 0.8, min chain-pair PAE < 1.5)
            const iptm = pdb.metrics?.iPTM as number | undefined;
            const minChainPairPae = pdb.metrics?.min_chain_pair_pae as number | undefined;

            const result = await api.analyzeDesign({
              pdb_content: pdb.pdbContent,
              backbone_pdb: backbonePdb,
              metal_type: targetMetal,
              ligand_name: ligandName,
              design_type: designType,
              filter_tier: filterTier,
              design_params: {},
              iptm: iptm,
              min_chain_pair_pae: minChainPairPae,
            });

            // Enrich metrics with analyzer output
            if (result.metrics.plddt !== undefined) metrics.pLDDT_analysis = result.metrics.plddt;
            if (result.metrics.coordination_number !== undefined) metrics.coordination = result.metrics.coordination_number;
            if (result.metrics.geometry_rmsd !== undefined) metrics.geometry_rmsd = result.metrics.geometry_rmsd;
            if (result.metrics.alanine_pct !== undefined) metrics.ala_pct = result.metrics.alanine_pct;
            if (result.metrics.aromatic_pct !== undefined) metrics.aromatic_pct = result.metrics.aromatic_pct;
            if (result.metrics.total_residues !== undefined) metrics.residues = result.metrics.total_residues;
            if (result.metrics.ligand_contacts !== undefined) metrics.ligand_contacts = result.metrics.ligand_contacts;
            if (result.metrics.protein_coordination !== undefined) metrics.protein_coord = result.metrics.protein_coordination;
            // AF3 paper ligand-interface metrics (shown in analysis results)
            if (result.metrics.iptm !== undefined) metrics.iPTM_filter = result.metrics.iptm;
            if (result.metrics.min_chain_pair_pae !== undefined) metrics.chain_pair_PAE = result.metrics.min_chain_pair_pae;

            metrics.filter_passed = result.filter_passed ? 'Yes' : 'No';
            metrics.filter_preset = result.filter_preset;
            filterPreset = result.filter_preset;

            if (result.filter_passed) passCount++;

            if (result.failed_filters.length > 0) {
              metrics.failed_filters = result.failed_filters
                .map(f => `${f.metric}=${typeof f.value === 'number' ? f.value.toFixed(2) : f.value}`)
                .join(', ');
            }
          } catch (err) {
            // Non-fatal: fall back to existing metrics
            console.warn(`[analysis] analyzeDesign failed for ${pdb.label}:`, err);
            metrics.filter_passed = 'N/A';
          }

          analyzed.push({ ...pdb, metrics });
        }

        // Sort: passing designs first, then by pLDDT descending
        analyzed.sort((a, b) => {
          const aPass = a.metrics?.filter_passed === 'Yes' ? 1 : 0;
          const bPass = b.metrics?.filter_passed === 'Yes' ? 1 : 0;
          if (aPass !== bPass) return bPass - aPass;
          const aScore = (a.metrics?.pLDDT as number) ?? (a.metrics?.pLDDT_analysis as number) ?? 0;
          const bScore = (b.metrics?.pLDDT as number) ?? (b.metrics?.pLDDT_analysis as number) ?? 0;
          return bScore - aScore;
        });

        ctx.onProgress(100, 'Analysis complete');

        const summary = `${passCount}/${toAnalyze.length} designs passed ${filterPreset} filter`;

        return {
          id: `analysis-${Date.now()}`,
          summary,
          pdbOutputs: analyzed,
          data: {
            total_analyzed: analyzed.length,
            designs_passed: passCount,
            designs_failed: analyzed.length - passCount,
            filter_preset: filterPreset,
            filter_tier: filterTier,
            design_type: designType,
            top_score: analyzed[0]?.metrics?.pLDDT ?? analyzed[0]?.metrics?.pLDDT_analysis,
          },
        };
      },
    },

    // Step 10: Save to design history (automatic, no user interaction)
    createSaveHistoryStep({
      id: 'save_history_nl',
      name: 'Saving Your Work',
      description: 'Storing your results so you can come back to them later',
    }),

    // Step 11: Check for lesson triggers (shows results if detected)
    createCheckLessonsStep({
      id: 'check_lessons_nl',
      name: 'Learning What Works',
      description: 'Looking for patterns across your designs to improve future runs',
    }),
  ],
};
