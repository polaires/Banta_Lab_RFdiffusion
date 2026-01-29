/**
 * Pipeline registry — single import point for all pipeline definitions.
 */

import type { PipelineDefinition } from '@/lib/pipeline-types';
import { metalBindingPipeline } from './metal-binding';
import { ligandDimerPipeline } from './ligand-dimer';
import { binderPipeline } from './binder';
import { metalScaffoldPipeline } from './metal-scaffold';
import { naturalLanguagePipeline } from './natural-language';

/** All available pipeline definitions, keyed by ID. */
export const pipelines: Record<string, PipelineDefinition> = {
  'metal-binding': metalBindingPipeline,
  'ligand-dimer': ligandDimerPipeline,
  'binder': binderPipeline,
  'metal-scaffold': metalScaffoldPipeline,
  'natural-language': naturalLanguagePipeline,
};

/** Get a pipeline definition by ID. */
export function getPipeline(id: string): PipelineDefinition | undefined {
  return pipelines[id];
}

/** Get all pipeline definitions as an array. */
export function getAllPipelines(): PipelineDefinition[] {
  return Object.values(pipelines);
}

// Re-export individual pipelines for direct imports
export {
  metalBindingPipeline,
  ligandDimerPipeline,
  binderPipeline,
  metalScaffoldPipeline,
  naturalLanguagePipeline,
};
