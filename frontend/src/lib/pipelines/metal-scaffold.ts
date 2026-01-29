/**
 * Metal Scaffold Sweep Pipeline
 * Steps: rfd3_backbone (sweep) → mpnn_sequence → rf3_validation → tier_assignment
 */

import { Grid3x3 } from 'lucide-react';
import api from '@/lib/api';
import { extractPdbContent } from '@/lib/workflowHandlers';
import type { PipelineDefinition, PdbOutput, SequenceOutput } from '@/lib/pipeline-types';
import { createMpnnStep, createRf3Step } from './shared-steps';

export const metalScaffoldPipeline: PipelineDefinition = {
  id: 'metal-scaffold',
  name: 'Metal Scaffold Sweep',
  description: 'Sweep across scaffold parameters to find optimal metal binding designs',
  icon: Grid3x3,
  requiresBackend: true,

  initialParams: [
    { id: 'pdb_content', label: 'PDB Content', type: 'hidden', required: true },
    { id: 'target_metal', label: 'Target Metal', type: 'text', required: true, defaultValue: 'ZN' },
    { id: 'contig', label: 'Contig', type: 'hidden', required: false },
    { id: 'ligand', label: 'Ligand', type: 'hidden', required: false },
    { id: 'unindex', label: 'Unindex', type: 'hidden', required: false },
    { id: 'partial_t', label: 'Partial T', type: 'hidden', required: false },
  ],

  steps: [
    // Step 1: RFD3 backbone sweep (iterates over multiple configs)
    {
      id: 'rfd3_sweep',
      name: 'Backbone Sweep',
      description: 'Run RFD3 with multiple parameter configurations',
      icon: Grid3x3,
      requiresReview: true,
      supportsSelection: true,
      optional: false,
      defaultParams: {
        configs_count: 3,
        designs_per_config: 2,
        num_timesteps: 200,
      },
      parameterSchema: [
        { id: 'configs_count', label: 'Parameter Configs', type: 'slider', required: false, defaultValue: 3, range: { min: 1, max: 9, step: 1 }, helpText: 'Number of parameter variations to try' },
        { id: 'designs_per_config', label: 'Designs per Config', type: 'slider', required: false, defaultValue: 2, range: { min: 1, max: 4, step: 1 } },
        { id: 'num_timesteps', label: 'Timesteps', type: 'slider', required: false, defaultValue: 200, range: { min: 50, max: 500, step: 50 } },
      ],

      async execute(ctx) {
        const { initialParams, params } = ctx;

        const configsCount = (params.configs_count as number) ?? 3;
        const designsPerConfig = (params.designs_per_config as number) ?? 2;
        const numTimesteps = (params.num_timesteps as number) ?? 200;

        // Generate sweep configurations
        const gammaValues = [0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.4, 0.55];
        const stepScaleValues = [1.0, 1.5, 2.0, 2.5, 1.2, 1.8, 0.8, 1.3, 1.7];

        const allOutputs: PdbOutput[] = [];

        for (let c = 0; c < Math.min(configsCount, gammaValues.length); c++) {
          const gamma = gammaValues[c];
          const stepScale = stepScaleValues[c];
          const configLabel = `\u03B3=${gamma} s=${stepScale}`;

          const progressBase = (c / configsCount) * 90;
          ctx.onProgress(progressBase + 2, `Config ${c + 1}/${configsCount}: ${configLabel}`);

          try {
            const response = await api.submitRFD3Design({
              ...initialParams as any,
              num_designs: designsPerConfig,
              num_timesteps: numTimesteps,
              gamma_0: gamma,
              step_scale: stepScale,
              is_non_loopy: true,
            });

            ctx.onJobCreated?.(response.job_id, 'rfd3');

            let result: Record<string, unknown>;

            if (response.syncCompleted) {
              if (response.status !== 'completed' || !response.result) continue;
              result = response.result as Record<string, unknown>;
            } else {
              const polled = await api.waitForJob(response.job_id, (s) => {
                if (s.status === 'running') {
                  ctx.onProgress(progressBase + (90 / configsCount) * 0.5, `${configLabel}: generating...`);
                }
              });
              if (polled.status !== 'completed') continue;
              result = polled.result as Record<string, unknown>;
            }

            const designs = result.designs as Array<Record<string, unknown>> | undefined;
            if (Array.isArray(designs)) {
              designs.forEach((d, idx) => {
                const pdb = (d.content as string) || (d.pdb_content as string);
                if (pdb) {
                  allOutputs.push({
                    id: `sweep-${c}-${idx}`,
                    label: `${configLabel} #${idx + 1}`,
                    pdbContent: pdb,
                    metrics: { gamma: gamma, step_scale: stepScale },
                  });
                }
              });
            } else {
              const pdb = extractPdbContent(result);
              if (pdb) {
                allOutputs.push({
                  id: `sweep-${c}-0`,
                  label: `${configLabel} #1`,
                  pdbContent: pdb,
                  metrics: { gamma: gamma, step_scale: stepScale },
                });
              }
            }
          } catch {
            // Continue with next config on failure
          }

          ctx.onProgress(progressBase + (90 / configsCount), `Config ${c + 1} done`);
        }

        ctx.onProgress(100, 'Sweep complete');

        return {
          id: `sweep-${Date.now()}`,
          summary: `Generated ${allOutputs.length} designs across ${configsCount} parameter configurations`,
          pdbOutputs: allOutputs,
          data: {
            total_designs: allOutputs.length,
            configs_tried: configsCount,
            designs_per_config: designsPerConfig,
          },
        };
      },
    },

    // Step 2: MPNN sequence design
    createMpnnStep({
      id: 'mpnn_scaffold',
      defaultParams: {
        temperature: 0.1,
        num_sequences: 2,
        model_type: 'ligand_mpnn',
      },
    }),

    // Step 3: RF3 validation
    createRf3Step({
      id: 'rf3_scaffold',
    }),

    // Step 4: Tier assignment
    {
      id: 'tier_assignment',
      name: 'Tier Assignment',
      description: 'Assign quality tiers based on validation metrics',
      icon: Grid3x3,
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
        const pdbs = selectedItems.length > 0
          ? allPdbs.filter(p => selectedItems.includes(p.id))
          : allPdbs;

        ctx.onProgress(20, 'Evaluating designs...');

        // Try to evaluate each design for metal coordination
        const tieredPdbs: PdbOutput[] = [];

        for (let i = 0; i < pdbs.length; i++) {
          const pdb = pdbs[i];
          ctx.onProgress(20 + (i / pdbs.length) * 70, `Evaluating ${pdb.label}...`);

          let tier = 'C';
          const metrics = { ...pdb.metrics };

          try {
            const evaluation = await api.evaluateDesign({
              pdb_content: pdb.pdbContent,
              target_metal: (initialParams.target_metal as string) || 'ZN',
            });

            metrics.coordination = evaluation.coordination_number;
            metrics.geometry = evaluation.geometry_type;
            metrics.bond_dist = evaluation.avg_bond_distance;

            // Tier assignment based on criteria
            if (evaluation.overall_pass && evaluation.criteria_passed >= evaluation.criteria_total * 0.8) {
              tier = 'A';
            } else if (evaluation.overall_pass || evaluation.criteria_passed >= evaluation.criteria_total * 0.5) {
              tier = 'B';
            }

            metrics.tier = tier;
            metrics.pass = evaluation.overall_pass ? 'Yes' : 'No';
          } catch {
            metrics.tier = 'C';
            metrics.pass = 'N/A';
          }

          tieredPdbs.push({ ...pdb, metrics });
        }

        // Sort by tier
        tieredPdbs.sort((a, b) => {
          const tierOrder: Record<string, number> = { A: 0, B: 1, C: 2 };
          return (tierOrder[a.metrics?.tier as string] ?? 2) - (tierOrder[b.metrics?.tier as string] ?? 2);
        });

        const tierA = tieredPdbs.filter(p => p.metrics?.tier === 'A').length;
        const tierB = tieredPdbs.filter(p => p.metrics?.tier === 'B').length;
        const tierC = tieredPdbs.filter(p => p.metrics?.tier === 'C').length;

        ctx.onProgress(100, 'Tier assignment complete');

        return {
          id: `tiers-${Date.now()}`,
          summary: `Tier A: ${tierA}, Tier B: ${tierB}, Tier C: ${tierC} (${tieredPdbs.length} total)`,
          pdbOutputs: tieredPdbs,
          data: {
            tier_a: tierA,
            tier_b: tierB,
            tier_c: tierC,
            total: tieredPdbs.length,
          },
        };
      },
    },
  ],
};
