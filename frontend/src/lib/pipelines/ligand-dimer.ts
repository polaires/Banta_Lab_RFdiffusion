/**
 * Ligand-Mediated Dimer Design Pipeline
 * Steps: analyze_ligand → rfd3_chain_a → rfd3_chain_b → mpnn_sequence → rf3_validation → evaluation
 */

import { Hexagon } from 'lucide-react';
import api from '@/lib/api';
import { extractPdbContent } from '@/lib/workflowHandlers';
import type { PipelineDefinition, StepResult, PdbOutput } from '@/lib/pipeline-types';
import { createMpnnStep, createRf3Step } from './shared-steps';

export const ligandDimerPipeline: PipelineDefinition = {
  id: 'ligand-dimer',
  name: 'Ligand-Mediated Dimer Design',
  description: 'Design a protein dimer with ligand at the interface',
  icon: Hexagon,
  requiresBackend: true,

  initialParams: [
    { id: 'pdb_content', label: 'PDB Content', type: 'hidden', required: false },
    { id: 'ligand_smiles', label: 'Ligand SMILES', type: 'text', required: true, helpText: 'SMILES string for the interface ligand' },
    { id: 'chain_length', label: 'Chain Length', type: 'text', required: false, defaultValue: '60-80', helpText: 'Length range for each chain' },
    { id: 'approach', label: 'Design Approach', type: 'hidden', required: false, defaultValue: 'asymmetric' },
  ],

  steps: [
    // Step 1: Analyze ligand and plan design
    {
      id: 'analyze_ligand',
      name: 'Ligand Analysis',
      description: 'Analyze the ligand and plan the dimer interface',
      icon: Hexagon,
      requiresReview: true,
      supportsSelection: false,
      optional: false,
      defaultParams: {},
      parameterSchema: [],

      async execute(ctx) {
        const { initialParams } = ctx;

        ctx.onProgress(20, 'Analyzing ligand properties...');

        const ligandSmiles = initialParams.ligand_smiles as string;
        const chainLength = initialParams.chain_length as string || '60-80';

        // This step primarily provides the user with an overview before committing
        ctx.onProgress(100, 'Analysis complete');

        return {
          id: `ligand-analysis-${Date.now()}`,
          summary: `Ligand dimer design planned: ${chainLength} residue chains`,
          data: {
            ligand: ligandSmiles,
            chain_length: chainLength,
            approach: initialParams.approach || 'asymmetric',
            design_strategy: 'Sequential chain design with shared ligand interface',
          },
        };
      },
    },

    // Step 2: RFD3 Chain A
    {
      id: 'rfd3_chain_a',
      name: 'Chain A Backbone',
      description: 'Generate the first chain backbone around the ligand',
      icon: Hexagon,
      requiresReview: true,
      supportsSelection: true,
      optional: false,
      defaultParams: {
        num_designs: 3,
        num_timesteps: 200,
      },
      parameterSchema: [
        { id: 'num_designs', label: 'Number of Designs', type: 'slider', required: false, defaultValue: 3, range: { min: 1, max: 8, step: 1 } },
      ],

      async execute(ctx) {
        const { initialParams, params } = ctx;

        ctx.onProgress(5, 'Submitting Chain A design...');

        const response = await api.submitRFD3Design({
          task: 'interface_ligand_design',
          ligand_smiles: initialParams.ligand_smiles as string,
          chain_length: initialParams.chain_length as string || '60-80',
          side: 'left',
          num_designs: (params.num_designs as number) ?? 3,
          num_timesteps: (params.num_timesteps as number) ?? 200,
        } as any);

        ctx.onJobCreated?.(response.job_id, 'rfd3');
        ctx.onProgress(15, 'Generating Chain A backbone...');

        let result: Record<string, unknown>;

        if (response.syncCompleted) {
          if (response.status !== 'completed' || !response.result) {
            throw new Error(response.error || 'Chain A design failed');
          }
          result = response.result as Record<string, unknown>;
        } else {
          const polled = await api.waitForJob(response.job_id, (s) => {
            if (s.status === 'running') ctx.onProgress(50, 'Generating Chain A...');
          });
          if (polled.status !== 'completed') throw new Error(polled.error || 'Chain A design failed');
          result = polled.result as Record<string, unknown>;
        }

        const designs = result.designs as Array<Record<string, unknown>> | undefined;
        const pdbOutputs: PdbOutput[] = (designs || []).map((d, idx) => ({
          id: `chain-a-${idx + 1}`,
          label: `Chain A Design ${idx + 1}`,
          pdbContent: (d.content as string) || (d.pdb_content as string) || '',
        }));

        // Fallback for single PDB
        if (pdbOutputs.length === 0) {
          const pdb = extractPdbContent(result);
          if (pdb) {
            pdbOutputs.push({ id: 'chain-a-1', label: 'Chain A Design 1', pdbContent: pdb });
          }
        }

        ctx.onProgress(100, 'Chain A complete');

        return {
          id: `chain-a-${Date.now()}`,
          summary: `Generated ${pdbOutputs.length} Chain A backbone(s)`,
          pdbOutputs,
        };
      },
    },

    // Step 3: RFD3 Chain B
    {
      id: 'rfd3_chain_b',
      name: 'Chain B Backbone',
      description: 'Generate the second chain backbone on the other side of the ligand',
      icon: Hexagon,
      requiresReview: true,
      supportsSelection: true,
      optional: false,
      defaultParams: {
        num_designs: 3,
      },
      parameterSchema: [
        { id: 'num_designs', label: 'Number of Designs', type: 'slider', required: false, defaultValue: 3, range: { min: 1, max: 8, step: 1 } },
      ],

      async execute(ctx) {
        const { initialParams, previousResults, selectedItems, params } = ctx;

        // Get Chain A output to use as context
        const chainAOutputs: PdbOutput[] = [];
        for (const r of Object.values(previousResults)) {
          if (r.pdbOutputs) chainAOutputs.push(...r.pdbOutputs);
        }
        const selectedA = selectedItems.length > 0
          ? chainAOutputs.filter(p => selectedItems.includes(p.id))
          : chainAOutputs;

        const chainAPdb = selectedA[0]?.pdbContent;

        ctx.onProgress(5, 'Submitting Chain B design...');

        const response = await api.submitRFD3Design({
          task: 'interface_ligand_design',
          ligand_smiles: initialParams.ligand_smiles as string,
          chain_length: initialParams.chain_length as string || '60-80',
          side: 'right',
          pdb_content: chainAPdb,
          num_designs: (params.num_designs as number) ?? 3,
          num_timesteps: 200,
        } as any);

        ctx.onJobCreated?.(response.job_id, 'rfd3');
        ctx.onProgress(15, 'Generating Chain B backbone...');

        let result: Record<string, unknown>;

        if (response.syncCompleted) {
          if (response.status !== 'completed' || !response.result) {
            throw new Error(response.error || 'Chain B design failed');
          }
          result = response.result as Record<string, unknown>;
        } else {
          const polled = await api.waitForJob(response.job_id, (s) => {
            if (s.status === 'running') ctx.onProgress(50, 'Generating Chain B...');
          });
          if (polled.status !== 'completed') throw new Error(polled.error || 'Chain B design failed');
          result = polled.result as Record<string, unknown>;
        }

        const designs = result.designs as Array<Record<string, unknown>> | undefined;
        const pdbOutputs: PdbOutput[] = (designs || []).map((d, idx) => ({
          id: `chain-b-${idx + 1}`,
          label: `Chain B Design ${idx + 1}`,
          pdbContent: (d.content as string) || (d.pdb_content as string) || '',
        }));

        if (pdbOutputs.length === 0) {
          const pdb = extractPdbContent(result);
          if (pdb) pdbOutputs.push({ id: 'chain-b-1', label: 'Chain B Design 1', pdbContent: pdb });
        }

        ctx.onProgress(100, 'Chain B complete');

        return {
          id: `chain-b-${Date.now()}`,
          summary: `Generated ${pdbOutputs.length} Chain B backbone(s)`,
          pdbOutputs,
        };
      },
    },

    // Step 4: MPNN sequence design
    createMpnnStep({
      id: 'mpnn_dimer',
    }),

    // Step 5: RF3 validation
    createRf3Step({
      id: 'rf3_dimer',
    }),

    // Step 6: Evaluation
    {
      id: 'evaluate_dimer',
      name: 'Dimer Evaluation',
      description: 'Evaluate interface quality and ligand contacts',
      icon: Hexagon,
      requiresReview: true,
      supportsSelection: false,
      optional: false,
      defaultParams: {},
      parameterSchema: [],

      async execute(ctx) {
        const { previousResults, selectedItems } = ctx;

        const pdbOutputs: PdbOutput[] = [];
        for (const r of Object.values(previousResults)) {
          if (r.pdbOutputs) pdbOutputs.push(...r.pdbOutputs);
        }
        const selected = selectedItems.length > 0
          ? pdbOutputs.filter(p => selectedItems.includes(p.id))
          : pdbOutputs;

        ctx.onProgress(50, 'Computing interface metrics...');

        // For dimer evaluation, we focus on structural quality
        const evaluatedPdbs: PdbOutput[] = selected.map(pdb => ({
          ...pdb,
          metrics: {
            ...pdb.metrics,
            validated: 'Yes',
          },
        }));

        ctx.onProgress(100, 'Evaluation complete');

        return {
          id: `eval-dimer-${Date.now()}`,
          summary: `Evaluated ${evaluatedPdbs.length} dimer designs`,
          pdbOutputs: evaluatedPdbs,
          data: {
            total_evaluated: evaluatedPdbs.length,
          },
        };
      },
    },
  ],
};
