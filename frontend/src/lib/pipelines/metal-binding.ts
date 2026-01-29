/**
 * Metal Binding Redesign Pipeline
 * Steps: analyze → rfd3_backbone → mpnn_sequence → rf3_validation → evaluate
 */

import { Atom } from 'lucide-react';
import api from '@/lib/api';
import type { PipelineDefinition, StepResult, StepExecutionContext } from '@/lib/pipeline-types';
import { getMetalMpnnConfig } from './metal-mpnn-bias';
import {
  createAnalyzeStep,
  createRfd3Step,
  createMpnnStep,
  createRf3Step,
  createEvaluateStep,
} from './shared-steps';

export const metalBindingPipeline: PipelineDefinition = {
  id: 'metal-binding',
  name: 'Metal Binding Redesign',
  description: 'Redesign a metal binding site with optimized coordination geometry',
  icon: Atom,
  requiresBackend: true,

  initialParams: [
    { id: 'pdb_content', label: 'PDB Content', type: 'hidden', required: true },
    { id: 'target_metal', label: 'Target Metal', type: 'text', required: true, defaultValue: 'ZN' },
    { id: 'metal_chain', label: 'Metal Chain', type: 'text', required: false },
    { id: 'metal_resnum', label: 'Metal Residue Number', type: 'number', required: false },
    { id: 'contig', label: 'Contig String', type: 'hidden', required: false },
    { id: 'ligand', label: 'Ligand', type: 'hidden', required: false },
    { id: 'unindex', label: 'Unindex', type: 'hidden', required: false },
    { id: 'partial_t', label: 'Partial T', type: 'hidden', required: false },
  ],

  steps: [
    // Step 1: Analyze the input structure
    createAnalyzeStep({
      id: 'analyze_metal',
      name: 'Metal Site Analysis',
      description: 'Analyze the metal binding site geometry and coordination',

      async execute(ctx) {
        const { initialParams } = ctx;
        const pdbContent = initialParams.pdb_content as string;
        const metalResidue = initialParams.target_metal as string || 'ZN';
        const metalChain = initialParams.metal_chain as string;
        const metalResnum = initialParams.metal_resnum as number;

        if (!pdbContent) {
          throw new Error('No PDB content provided');
        }

        ctx.onProgress(10, 'Analyzing metal binding site...');

        let analysisData: Record<string, unknown> = {};

        // Try metal-specific analysis first
        if (metalResnum) {
          try {
            const metalAnalysis = await api.analyzeMetalBinding({
              pdb_content: pdbContent,
              metal_residue: metalResidue,
              metal_chain: metalChain,
              metal_resnum: metalResnum,
            });

            ctx.onProgress(60, 'Metal binding analysis complete');

            analysisData = {
              coordination_number: metalAnalysis.coordination?.number,
              geometry: metalAnalysis.coordination?.geometry,
              geometry_rmsd: metalAnalysis.coordination?.geometry_rmsd,
              avg_bond_distance: metalAnalysis.bond_analysis?.average_distance,
              dominant_donor: metalAnalysis.donor_analysis?.dominant_type,
              coordinating_residues: metalAnalysis.coordination?.coordinating_atoms
                ?.map(a => `${a.chain}${a.residue_number} ${a.residue} (${a.atom})`)
                .join(', '),
            };
          } catch {
            // Fall back to general analysis
          }
        }

        // Also run general structure analysis
        ctx.onProgress(70, 'Running general structure analysis...');
        let rfd3Params: Record<string, unknown> = {};
        try {
          const general = await api.analyzeStructure(pdbContent);
          analysisData.num_residues = general.num_residues;
          analysisData.num_ligands = general.ligands?.length ?? 0;
          analysisData.num_binding_sites = general.binding_sites?.length ?? 0;

          // Extract RFD3-specific params from the first suggestion
          if (general.suggestions?.length > 0) {
            const suggestion = general.suggestions[0];
            rfd3Params = { ...suggestion.rfd3_params };
            analysisData.suggestion_type = suggestion.type;
            analysisData.suggestion_description = suggestion.description;
          }

          // Build a contig from the structure if not in suggestions
          if (!rfd3Params.contig && !rfd3Params.length && general.num_residues) {
            // Use full chain as contig for partial diffusion redesign
            const chains = pdbContent.match(/^ATOM.{17}(.)/gm);
            const firstChain = chains?.[0]?.charAt(21) || 'A';
            rfd3Params.contig = `${firstChain}1-${general.num_residues}/0 50-50`;
          }
        } catch {
          // Non-fatal
        }

        ctx.onProgress(100, 'Analysis complete');

        // Unified MPNN config from HSAB classification (single source of truth)
        const mpnnParams = getMetalMpnnConfig(metalResidue);

        return {
          id: `metal-analysis-${Date.now()}`,
          summary: analysisData.coordination_number
            ? `${metalResidue} site: ${analysisData.geometry} coordination (${analysisData.coordination_number}-fold)`
            : `Analyzed structure with ${analysisData.num_residues ?? '?'} residues`,
          data: {
            ...analysisData,
            rfd3_params: rfd3Params,
            mpnn_params: mpnnParams,
          },
        };
      },
    }),

    // Step 2: RFD3 backbone generation
    createRfd3Step({
      id: 'rfd3_metal',
      name: 'Backbone Generation',
      description: 'Generate redesigned backbone structures around the metal site',
      defaultParams: {
        num_designs: 4,
        num_timesteps: 200,
        step_scale: 1.5,
        gamma_0: 0.6,
        is_non_loopy: true,
      },
    }),

    // Step 3: MPNN sequence design
    createMpnnStep({
      id: 'mpnn_metal',
      defaultParams: {
        temperature: 0.1,
        num_sequences: 4,
        model_type: 'ligand_mpnn',
      },
    }),

    // Step 4: RF3 validation
    createRf3Step({
      id: 'rf3_metal',
    }),

    // Step 5: Evaluation
    createEvaluateStep({
      id: 'evaluate_metal',
      name: 'Metal Site Evaluation',
      description: 'Evaluate coordination geometry and binding quality',
    }),
  ],
};
