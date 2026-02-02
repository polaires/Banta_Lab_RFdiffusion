/**
 * Chat & Project data model types for the Claude-style chat interface.
 */

import type { StepResult } from './pipeline-types';

// ---- Chat Messages ----

export type ChatAttachment =
  | { type: 'pipeline_progress'; pipelineId: string }
  | { type: 'step_result'; stepId: string; result: StepResult }
  | { type: 'design_gallery'; designs: Array<{ id: string; name: string; pdbContent?: string; metrics?: Record<string, number | string> }> }
  | { type: 'evaluation_card'; data: Record<string, unknown> }
  | { type: 'interview'; interviewType: string }
  | { type: 'quick_start' }
  | { type: 'preference_summary'; preferences: Record<string, unknown>; designType: string };

export interface ChatMessage {
  id: string;
  role: 'user' | 'assistant' | 'system';
  content: string;
  timestamp: number;
  attachment?: ChatAttachment;
}

// ---- Project ----

export type ProjectStatus = 'active' | 'completed' | 'failed' | 'archived';

export interface ProjectContext {
  pdbId: string | null;
  pdbContent: string | null;
  targetMetal: string | null;
  nlQuery: string | null;
  designType: string | null; // 'metal_binding' | 'ligand_dimer' | 'binder' | 'nl' | etc.
  userPreferences: Record<string, unknown> | null;
  analysisResult: Record<string, unknown> | null;
  pipelineResults: Record<string, StepResult>;
}

export interface Project {
  id: string;
  name: string;
  status: ProjectStatus;
  createdAt: number;
  updatedAt: number;
  messages: ChatMessage[];
  activePipelineId: string | null;
  context: ProjectContext;
  jobIds: string[];
}

// ---- Helpers ----

let _counter = 0;

export function generateId(prefix = 'msg'): string {
  _counter++;
  return `${prefix}-${Date.now()}-${_counter}-${Math.random().toString(36).slice(2, 7)}`;
}

export function createEmptyContext(): ProjectContext {
  return {
    pdbId: null,
    pdbContent: null,
    targetMetal: null,
    nlQuery: null,
    designType: null,
    userPreferences: null,
    analysisResult: null,
    pipelineResults: {},
  };
}

export function createProject(name?: string): Project {
  const id = generateId('proj');
  return {
    id,
    name: name || 'New Design',
    status: 'active',
    createdAt: Date.now(),
    updatedAt: Date.now(),
    messages: [],
    activePipelineId: null,
    context: createEmptyContext(),
    jobIds: [],
  };
}

/** Generate a project name from the first user message */
export function deriveProjectName(content: string): string {
  // Take first ~40 chars, trim to word boundary
  const trimmed = content.trim().slice(0, 50);
  const lastSpace = trimmed.lastIndexOf(' ');
  const name = lastSpace > 20 ? trimmed.slice(0, lastSpace) : trimmed;
  return name + (content.length > 50 ? '...' : '');
}
