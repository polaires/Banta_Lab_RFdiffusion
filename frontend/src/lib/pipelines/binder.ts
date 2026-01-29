/**
 * Protein Binder Design Pipeline
 * Steps: hotspot_selection → rfd3_backbone → mpnn_sequence → esm3_scoring → interface_analysis → filtering
 */

import { Target } from 'lucide-react';
import api from '@/lib/api';
import { extractPdbContent, parseBinderDesigns } from '@/lib/workflowHandlers';
import type { PipelineDefinition, PdbOutput, SequenceOutput } from '@/lib/pipeline-types';

export const binderPipeline: PipelineDefinition = {
  id: 'binder',
  name: 'Protein Binder Design',
  description: 'Design a de novo protein binder against a target protein',
  icon: Target,
  requiresBackend: true,

  initialParams: [
    { id: 'pdb_content', label: 'Target PDB', type: 'hidden', required: true },
    { id: 'target_chain', label: 'Target Chain', type: 'text', required: false, defaultValue: 'A' },
    { id: 'binder_length', label: 'Binder Length', type: 'text', required: false, defaultValue: '60-80', helpText: 'Length range (e.g., 60-80)' },
    { id: 'quality_threshold', label: 'Quality', type: 'hidden', required: false, defaultValue: 'standard' },
  ],

  steps: [
    // Step 1: Hotspot detection
    {
      id: 'hotspot_selection',
      name: 'Hotspot Detection',
      description: 'Identify optimal binding hotspots on the target surface',
      icon: Target,
      requiresReview: true,
      supportsSelection: false,
      optional: false,
      defaultParams: {
        method: 'exposed_clustered',
        max_hotspots: 10,
      },
      parameterSchema: [
        {
          id: 'method', label: 'Detection Method', type: 'select', required: false, defaultValue: 'exposed_clustered',
          options: [
            { value: 'exposed_clustered', label: 'Exposed Clustered' },
            { value: 'exposed', label: 'All Exposed' },
            { value: 'patch', label: 'Surface Patch' },
          ],
        },
        { id: 'max_hotspots', label: 'Max Hotspots', type: 'slider', required: false, defaultValue: 10, range: { min: 3, max: 20, step: 1 } },
      ],

      async execute(ctx) {
        const { initialParams, params } = ctx;

        ctx.onProgress(10, 'Detecting surface hotspots...');

        const result = await api.detectHotspots({
          target_pdb: initialParams.pdb_content as string,
          target_chain: initialParams.target_chain as string || 'A',
          method: (params.method as string) as any || 'exposed_clustered',
          max_hotspots: (params.max_hotspots as number) ?? 10,
        });

        ctx.onProgress(100, 'Hotspot detection complete');

        return {
          id: `hotspots-${Date.now()}`,
          summary: `Found ${result.hotspots?.length ?? 0} hotspot residues using ${result.method}`,
          data: {
            hotspots: result.hotspots?.join(', ') || 'None',
            total_exposed: result.total_exposed,
            method: result.method,
            patch_area: result.patch_area,
            residue_details: result.residue_details?.slice(0, 5).map(
              r => `${r.residue} (${r.restype}, SASA: ${r.relative_sasa?.toFixed(2)})`
            ).join(', '),
          },
        };
      },
    },

    // Step 2: RFD3 backbone generation for binder
    {
      id: 'rfd3_binder',
      name: 'Binder Backbone Design',
      description: 'Generate binder backbone structures using RFdiffusion3',
      icon: Target,
      requiresReview: true,
      supportsSelection: true,
      optional: false,
      defaultParams: {
        num_designs: 4,
        num_timesteps: 200,
      },
      parameterSchema: [
        { id: 'num_designs', label: 'Number of Designs', type: 'slider', required: false, defaultValue: 4, range: { min: 1, max: 10, step: 1 } },
        { id: 'num_timesteps', label: 'Timesteps', type: 'slider', required: false, defaultValue: 200, range: { min: 50, max: 500, step: 50 } },
      ],

      async execute(ctx) {
        const { initialParams, previousResults, params } = ctx;

        // Extract hotspots from analysis step
        const hotspotData = Object.values(previousResults).find(r => r.data?.hotspots);
        const hotspotsStr = hotspotData?.data?.hotspots as string || '';
        const hotspots = hotspotsStr ? hotspotsStr.split(', ').filter(Boolean) : undefined;

        ctx.onProgress(5, 'Submitting binder design...');

        const response = await api.submitRFD3Design({
          task: 'protein_binder_design',
          target_pdb: initialParams.pdb_content as string,
          binder_length: initialParams.binder_length as string || '60-80',
          hotspots,
          num_designs: (params.num_designs as number) ?? 4,
          num_timesteps: (params.num_timesteps as number) ?? 200,
          quality_threshold: (initialParams.quality_threshold as string) || 'standard',
        } as any);

        ctx.onJobCreated?.(response.job_id, 'rfd3');
        ctx.onProgress(15, 'Generating binder backbones...');

        let result: Record<string, unknown>;

        if (response.syncCompleted) {
          if (response.status !== 'completed' || !response.result) {
            throw new Error(response.error || 'Binder design failed');
          }
          result = response.result as Record<string, unknown>;
        } else {
          const polled = await api.waitForJob(response.job_id, (s) => {
            if (s.status === 'running') ctx.onProgress(50, 'Generating binder structures...');
          });
          if (polled.status !== 'completed') throw new Error(polled.error || 'Binder design failed');
          result = polled.result as Record<string, unknown>;
        }

        // Parse binder designs
        const binderDesigns = parseBinderDesigns(result);
        const pdbOutputs: PdbOutput[] = binderDesigns
          .filter(d => d.pdb_content)
          .map((d, idx) => ({
            id: `binder-${idx + 1}`,
            label: `Binder ${d.rank}`,
            pdbContent: d.pdb_content!,
            metrics: {
              ...(d.mpnn_score !== undefined ? { mpnn: d.mpnn_score } : {}),
              ...(d.interface_contacts !== undefined ? { contacts: d.interface_contacts } : {}),
              ...(d.buried_sasa !== undefined ? { bSASA: d.buried_sasa } : {}),
            },
          }));

        // Fallback for generic designs array
        if (pdbOutputs.length === 0) {
          const designs = result.designs as Array<Record<string, unknown>> | undefined;
          if (Array.isArray(designs)) {
            designs.forEach((d, idx) => {
              const pdb = (d.content as string) || (d.pdb_content as string);
              if (pdb) {
                pdbOutputs.push({
                  id: `binder-${idx + 1}`,
                  label: `Binder ${idx + 1}`,
                  pdbContent: pdb,
                });
              }
            });
          }
        }

        ctx.onProgress(100, 'Binder backbone generation complete');

        return {
          id: `binder-rfd3-${Date.now()}`,
          summary: `Generated ${pdbOutputs.length} binder designs`,
          pdbOutputs,
          data: {
            num_designs: pdbOutputs.length,
            hotspots_used: hotspots?.length ?? 0,
          },
        };
      },
    },

    // Step 3: MPNN sequence design
    {
      id: 'mpnn_binder',
      name: 'Binder Sequence Design',
      description: 'Design sequences for binder backbones',
      icon: Target,
      requiresReview: true,
      supportsSelection: true,
      optional: false,
      defaultParams: {
        temperature: 0.1,
        num_sequences: 4,
        model_type: 'protein_mpnn',
      },
      parameterSchema: [
        { id: 'temperature', label: 'Temperature', type: 'slider', required: false, defaultValue: 0.1, range: { min: 0.01, max: 1.0, step: 0.01 } },
        { id: 'num_sequences', label: 'Sequences per Design', type: 'slider', required: false, defaultValue: 4, range: { min: 1, max: 16, step: 1 } },
      ],

      async execute(ctx) {
        const { previousResults, selectedItems, params } = ctx;

        const allPdbs: PdbOutput[] = [];
        for (const r of Object.values(previousResults)) {
          if (r.pdbOutputs) allPdbs.push(...r.pdbOutputs);
        }
        const pdbs = selectedItems.length > 0
          ? allPdbs.filter(p => selectedItems.includes(p.id))
          : allPdbs;

        if (pdbs.length === 0) throw new Error('No binder structures for sequence design');

        const allSequences: SequenceOutput[] = [];

        for (let i = 0; i < pdbs.length; i++) {
          const pdb = pdbs[i];
          ctx.onProgress((i / pdbs.length) * 85 + 5, `Designing sequences for ${pdb.label}...`);

          const response = await api.submitMPNNDesign({
            pdb_content: pdb.pdbContent,
            temperature: (params.temperature as number) ?? 0.1,
            num_sequences: (params.num_sequences as number) ?? 4,
            model_type: (params.model_type as 'protein_mpnn') ?? 'protein_mpnn',
          });

          ctx.onJobCreated?.(response.job_id, 'mpnn');

          let jobResult: Record<string, unknown>;
          if (response.syncCompleted) {
            if (response.status !== 'completed') throw new Error(response.error || 'MPNN failed');
            jobResult = response.result as Record<string, unknown>;
          } else {
            const polled = await api.waitForJob(response.job_id);
            if (polled.status !== 'completed') throw new Error(polled.error || 'MPNN failed');
            jobResult = polled.result as Record<string, unknown>;
          }

          const sequences = jobResult.sequences as Array<Record<string, unknown>> | undefined;
          if (Array.isArray(sequences)) {
            sequences.forEach((s, idx) => {
              allSequences.push({
                id: `binder-seq-${i}-${idx}`,
                sequence: (s.content as string) || (s.sequence as string) || '',
                score: s.score as number | undefined,
                label: `${pdb.label} Seq ${idx + 1}`,
              });
            });
          }
        }

        ctx.onProgress(100, 'Sequence design complete');

        return {
          id: `binder-mpnn-${Date.now()}`,
          summary: `Designed ${allSequences.length} sequences for ${pdbs.length} binder(s)`,
          sequences: allSequences,
        };
      },
    },

    // Step 4: ESM3 scoring
    {
      id: 'esm3_scoring',
      name: 'ESM3 Sequence Scoring',
      description: 'Score sequences using ESM3 perplexity',
      icon: Target,
      requiresReview: true,
      supportsSelection: true,
      optional: true,
      defaultParams: {},
      parameterSchema: [],

      async execute(ctx) {
        const { previousResults, selectedItems } = ctx;

        const allSeqs: SequenceOutput[] = [];
        for (const r of Object.values(previousResults)) {
          if (r.sequences) allSeqs.push(...r.sequences);
        }
        const seqs = selectedItems.length > 0
          ? allSeqs.filter(s => selectedItems.includes(s.id))
          : allSeqs;

        if (seqs.length === 0) throw new Error('No sequences for scoring');

        ctx.onProgress(10, 'Scoring sequences with ESM3...');

        const sequences = seqs.map(s => s.sequence);

        try {
          const result = await api.esm3ScoreSequences(sequences);

          ctx.onProgress(100, 'Scoring complete');

          const scoredSeqs: SequenceOutput[] = seqs.map((seq, idx) => {
            const scoreData = result.scores?.[idx];
            return {
              ...seq,
              score: scoreData?.perplexity,
              metrics: {
                perplexity: scoreData?.perplexity ?? 0,
                esm_score: scoreData?.score ?? 0,
              },
            };
          });

          // Sort by perplexity (lower is better)
          scoredSeqs.sort((a, b) => (a.score ?? Infinity) - (b.score ?? Infinity));

          return {
            id: `esm3-${Date.now()}`,
            summary: `Scored ${scoredSeqs.length} sequences (best perplexity: ${scoredSeqs[0]?.score?.toFixed(2) ?? 'N/A'})`,
            sequences: scoredSeqs,
            data: {
              total_scored: scoredSeqs.length,
              best_perplexity: scoredSeqs[0]?.score,
              worst_perplexity: scoredSeqs[scoredSeqs.length - 1]?.score,
            },
          };
        } catch {
          // ESM3 scoring is optional; return sequences as-is
          ctx.onProgress(100, 'ESM3 scoring skipped (service unavailable)');

          return {
            id: `esm3-${Date.now()}`,
            summary: `ESM3 scoring unavailable; ${seqs.length} sequences passed through`,
            sequences: seqs,
            data: { scoring_skipped: true },
          };
        }
      },
    },

    // Step 5: Interface analysis (via RF3 validation)
    {
      id: 'interface_analysis',
      name: 'Structure Validation',
      description: 'Validate binder structures via RF3 prediction',
      icon: Target,
      requiresReview: true,
      supportsSelection: true,
      optional: false,
      defaultParams: {},
      parameterSchema: [],

      async execute(ctx) {
        const { previousResults, selectedItems } = ctx;

        const allSeqs: SequenceOutput[] = [];
        for (const r of Object.values(previousResults)) {
          if (r.sequences) allSeqs.push(...r.sequences);
        }
        const seqs = selectedItems.length > 0
          ? allSeqs.filter(s => selectedItems.includes(s.id))
          : allSeqs;

        // Take top sequences (limit to avoid excessive computation)
        const topSeqs = seqs.slice(0, 8);
        const pdbOutputs: PdbOutput[] = [];

        for (let i = 0; i < topSeqs.length; i++) {
          const seq = topSeqs[i];
          ctx.onProgress((i / topSeqs.length) * 90 + 5, `Validating ${seq.label ?? `seq ${i + 1}`}...`);

          try {
            const response = await api.submitRF3Prediction({
              sequence: seq.sequence,
              name: seq.label ?? seq.id,
            });

            ctx.onJobCreated?.(response.job_id, 'rf3');

            let jobResult: Record<string, unknown>;
            if ((response as any).syncCompleted) {
              jobResult = response.result as Record<string, unknown>;
            } else {
              const polled = await api.waitForJob(response.job_id);
              if (polled.status !== 'completed') continue;
              jobResult = polled.result as Record<string, unknown>;
            }

            const predictions = jobResult?.predictions as Array<Record<string, unknown>> | undefined;
            const pdb = predictions?.[0]?.content as string || extractPdbContent(jobResult);
            const confidences = jobResult?.confidences as Record<string, any> | undefined;
            const summary = confidences?.summary_confidences;

            if (pdb) {
              pdbOutputs.push({
                id: `validated-${seq.id}`,
                label: seq.label ?? `Validated ${i + 1}`,
                pdbContent: pdb,
                sequence: seq.sequence,
                metrics: {
                  ...(summary?.overall_plddt ? { pLDDT: summary.overall_plddt } : {}),
                  ...(summary?.ptm ? { pTM: summary.ptm } : {}),
                  ...(seq.score ? { perplexity: seq.score } : {}),
                },
              });
            }
          } catch {
            // Continue with other sequences
          }
        }

        ctx.onProgress(100, 'Validation complete');

        return {
          id: `interface-${Date.now()}`,
          summary: `Validated ${pdbOutputs.length}/${topSeqs.length} binder structures`,
          pdbOutputs,
          data: {
            validated: pdbOutputs.length,
            attempted: topSeqs.length,
          },
        };
      },
    },

    // Step 6: Final filtering
    {
      id: 'filtering',
      name: 'Final Ranking',
      description: 'Rank and filter the final binder designs',
      icon: Target,
      requiresReview: true,
      supportsSelection: false,
      optional: false,
      defaultParams: {},
      parameterSchema: [],

      async execute(ctx) {
        const { previousResults, selectedItems } = ctx;

        const allPdbs: PdbOutput[] = [];
        for (const r of Object.values(previousResults)) {
          if (r.pdbOutputs) allPdbs.push(...r.pdbOutputs);
        }
        const pdbs = selectedItems.length > 0
          ? allPdbs.filter(p => selectedItems.includes(p.id))
          : allPdbs;

        ctx.onProgress(50, 'Ranking designs...');

        // Sort by pLDDT descending
        const ranked = [...pdbs].sort((a, b) => {
          const aScore = (a.metrics?.pLDDT as number) ?? 0;
          const bScore = (b.metrics?.pLDDT as number) ?? 0;
          return bScore - aScore;
        });

        // Add rank to metrics
        const rankedPdbs: PdbOutput[] = ranked.map((pdb, idx) => ({
          ...pdb,
          metrics: { ...pdb.metrics, rank: idx + 1 },
        }));

        ctx.onProgress(100, 'Ranking complete');

        const topPlddt = rankedPdbs[0]?.metrics?.pLDDT;

        return {
          id: `ranking-${Date.now()}`,
          summary: `Ranked ${rankedPdbs.length} binder designs by structural quality`,
          pdbOutputs: rankedPdbs,
          data: {
            total_ranked: rankedPdbs.length,
            top_plddt: topPlddt,
          },
        };
      },
    },
  ],
};
