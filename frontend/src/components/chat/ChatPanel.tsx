'use client';

import { useState, useCallback, useRef, useEffect, useMemo } from 'react';
import {
  Sparkles,
  ChevronDown,
  ChevronUp,
  Eye,
  Trash2,
} from 'lucide-react';
import { useStore } from '@/lib/store';
import { cn } from '@/lib/utils';
import type { ChatMessage as ChatMessageType, ChatAttachment } from '@/lib/chat-types';
import type { StepResult } from '@/lib/pipeline-types';
import { ChatThread } from './ChatThread';
import { ChatInput } from './ChatInput';
import { PipelineCard } from './PipelineCard';
import { InitiativePanel } from './InitiativePanel';
import { GuidedDesignWizard } from './GuidedDesignWizard';
import type { WizardResult } from './GuidedDesignWizard';

// ---- Outputs Summary ----

function OutputsSummary({
  messages,
  onSelectDesign,
  onDeleteProject,
}: {
  messages: ChatMessageType[];
  onSelectDesign: (pdbContent: string) => void;
  onDeleteProject?: () => void;
}) {
  const [expanded, setExpanded] = useState(false);

  const designs: Array<{ id: string; name: string; pdbContent?: string; metrics?: Record<string, number | string> }> = [];
  const stepSummaries: Array<{ stepId: string; summary: string }> = [];
  let designType: string | null = null;

  for (const msg of messages) {
    if (!msg.attachment) continue;
    if (msg.attachment.type === 'step_result') {
      stepSummaries.push({ stepId: msg.attachment.stepId, summary: msg.attachment.result.summary });
      if (msg.attachment.result.pdbOutputs) {
        for (const pdb of msg.attachment.result.pdbOutputs) {
          designs.push({ id: pdb.id, name: pdb.label || pdb.id, pdbContent: pdb.pdbContent });
        }
      }
    }
    if (msg.attachment.type === 'design_gallery') {
      for (const d of msg.attachment.designs) {
        if (!designs.some(existing => existing.id === d.id)) {
          designs.push(d);
        }
      }
    }
    if (msg.attachment.type === 'preference_summary') {
      designType = msg.attachment.designType;
    }
  }

  if (stepSummaries.length === 0 && designs.length === 0) return null;

  return (
    <div className="border-b border-border">
      <button
        onClick={() => setExpanded(!expanded)}
        className="w-full flex items-center gap-2 px-6 py-2 text-xs text-muted-foreground hover:bg-muted/50 transition-colors"
      >
        <Eye className="h-3.5 w-3.5" />
        <span className="font-medium">
          {designs.length > 0
            ? `${designs.length} design${designs.length !== 1 ? 's' : ''} generated`
            : `${stepSummaries.length} step${stepSummaries.length !== 1 ? 's' : ''} completed`}
        </span>
        {designType && (
          <span className="bg-primary/10 text-primary px-1.5 py-0.5 rounded text-[10px] font-medium">
            {designType.replace(/_/g, ' ')}
          </span>
        )}
        <span className="flex-1" />
        {onDeleteProject && (
          <span
            role="button"
            onClick={(e) => { e.stopPropagation(); onDeleteProject(); }}
            className="p-1 hover:text-destructive rounded"
            title="Delete project"
          >
            <Trash2 className="h-3 w-3" />
          </span>
        )}
        {expanded ? <ChevronUp className="h-3.5 w-3.5" /> : <ChevronDown className="h-3.5 w-3.5" />}
      </button>

      {expanded && (
        <div className="px-6 pb-3 space-y-2">
          {stepSummaries.length > 0 && (
            <div className="space-y-1">
              <div className="text-[10px] font-medium text-muted-foreground uppercase tracking-wider">Pipeline Steps</div>
              {stepSummaries.map((s, i) => (
                <div key={i} className="flex items-center gap-2 text-xs text-muted-foreground">
                  <span className="font-mono text-[10px] bg-muted px-1 rounded">{s.stepId}</span>
                  <span className="truncate">{s.summary}</span>
                </div>
              ))}
            </div>
          )}

          {designs.length > 0 && (
            <div className="space-y-1">
              <div className="text-[10px] font-medium text-muted-foreground uppercase tracking-wider">
                Designs ({designs.length})
              </div>
              <div className="max-h-40 overflow-auto space-y-0.5">
                {designs.map((d) => (
                  <button
                    key={d.id}
                    onClick={() => d.pdbContent && onSelectDesign(d.pdbContent)}
                    disabled={!d.pdbContent}
                    className={cn(
                      'w-full text-left px-2 py-1 rounded text-xs transition-colors',
                      d.pdbContent
                        ? 'hover:bg-muted cursor-pointer'
                        : 'opacity-50 cursor-not-allowed'
                    )}
                  >
                    <span className="font-medium text-foreground">{d.name}</span>
                    {d.metrics && (
                      <span className="text-muted-foreground ml-2">
                        {Object.entries(d.metrics).slice(0, 3).map(([k, v]) => `${k}: ${v}`).join(' | ')}
                      </span>
                    )}
                  </button>
                ))}
              </div>
            </div>
          )}
        </div>
      )}
    </div>
  );
}

// ---- Main ChatPanel ----

interface ChatPanelProps {
  className?: string;
}

export function ChatPanel({ className }: ChatPanelProps) {
  const {
    activeProjectId,
    projects,
    addProjectMessage,
    updateProjectContext,
    setProjectPipeline,
    updateProjectStatus,
    addProjectJobId,
    setSelectedPdb,
    createProject,
    setActiveProject,
    deleteProject,
  } = useStore();

  const [initialMsgCount] = useState(() => {
    const proj = useStore.getState().projects.find((p) => p.id === useStore.getState().activeProjectId);
    return proj?.messages.length ?? 0;
  });
  const prevMsgCountRef = useRef(initialMsgCount);

  const [showWizard, setShowWizard] = useState(false);

  const project = useMemo(
    () => projects.find((p) => p.id === activeProjectId) ?? null,
    [projects, activeProjectId]
  );

  const messages = project?.messages ?? [];
  const activePipelineId = project?.activePipelineId ?? null;

  useEffect(() => {
    prevMsgCountRef.current = messages.length;
  }, [messages.length]);

  // ---- Helpers ----

  const addMsg = useCallback(
    (role: 'user' | 'assistant' | 'system', content: string, attachment?: ChatAttachment) => {
      if (!activeProjectId) return;
      addProjectMessage(activeProjectId, { role, content, attachment });
    },
    [activeProjectId, addProjectMessage]
  );

  // ---- Pipeline lifecycle (defined before handleSend so ref is available) ----

  const startPipelineRef = useRef<(projId: string, pipelineId: string, params: Record<string, unknown>) => void>(() => {});

  const startPipeline = useCallback(
    (projId: string, pipelineId: string, params: Record<string, unknown>) => {
      setProjectPipeline(projId, pipelineId);
      updateProjectContext(projId, {
        nlQuery: typeof params.user_prompt === 'string' ? params.user_prompt : null,
        designType: pipelineId === 'natural-language' ? 'nl' : pipelineId,
        userPreferences: params as Record<string, unknown>,
      });
      updateProjectStatus(projId, 'active');
    },
    [setProjectPipeline, updateProjectContext, updateProjectStatus]
  );

  startPipelineRef.current = startPipeline;

  // ---- Handle user input ----

  const handleSend = useCallback(
    (text: string) => {
      let projId = activeProjectId;
      if (!projId) {
        projId = createProject();
      }

      addProjectMessage(projId, { role: 'user', content: text });
      setShowWizard(false);

      if (!activePipelineId) {
        addProjectMessage(projId, {
          role: 'assistant',
          content: 'Starting the design pipeline for your request...',
        });
        startPipelineRef.current(projId, 'natural-language', { user_prompt: text });
      }
    },
    [activeProjectId, addProjectMessage, createProject, activePipelineId]
  );

  // ---- Wizard completion ----

  const handleWizardComplete = useCallback(
    (result: WizardResult) => {
      let projId = activeProjectId;
      if (!projId) {
        projId = createProject();
      }

      addProjectMessage(projId, {
        role: 'user',
        content: result.user_prompt,
      });

      const prefDisplay: Record<string, string> = {};
      if (result.target_metal) prefDisplay['Metal'] = result.target_metal;
      if (result.ligand_name) prefDisplay['Ligand'] = result.ligand_name;
      if (result.design_goal) prefDisplay['Goal'] = result.design_goal;
      if (result.filter_tier) prefDisplay['Quality'] = result.filter_tier;
      if (result.num_designs) prefDisplay['Designs'] = String(result.num_designs);

      addProjectMessage(projId, {
        role: 'assistant',
        content: 'Design configured. Starting the pipeline...',
        attachment: {
          type: 'preference_summary',
          designType: result.design_type,
          preferences: prefDisplay,
        },
      });

      updateProjectContext(projId, {
        targetMetal: result.target_metal || null,
        designType: result.design_type,
        userPreferences: result as unknown as Record<string, unknown>,
      });

      startPipelineRef.current(projId, 'natural-language', {
        user_prompt: result.user_prompt,
        target_metal: result.target_metal,
        ligand: result.ligand_name,
        num_designs: result.num_designs,
        use_sweep: result.use_sweep,
        design_goal: result.design_goal,
        filter_tier: result.filter_tier,
        symmetry: result.symmetry,
        target_pdb_id: result.target_pdb_id,
        protocol: result.protocol,
      });

      setShowWizard(false);
    },
    [activeProjectId, addProjectMessage, createProject, updateProjectContext]
  );

  const handleWizardCancel = useCallback(() => {
    setShowWizard(false);
  }, []);

  const handleStepComplete = useCallback(
    (stepId: string, result: StepResult) => {
      if (!activeProjectId) return;
      addMsg('assistant', `Step "${stepId}" completed: ${result.summary}`);
      updateProjectContext(activeProjectId, {
        pipelineResults: {
          ...(project?.context.pipelineResults ?? {}),
          [stepId]: result,
        },
      });
    },
    [activeProjectId, addMsg, updateProjectContext, project]
  );

  const handlePipelineComplete = useCallback(
    (results: Record<string, StepResult>) => {
      if (!activeProjectId) return;
      setProjectPipeline(activeProjectId, null);
      updateProjectStatus(activeProjectId, 'completed');
      addMsg('assistant', 'Pipeline completed successfully! Review the final results above.');
    },
    [activeProjectId, setProjectPipeline, updateProjectStatus, addMsg]
  );

  const handlePipelineFailed = useCallback(
    (stepId: string, error: string) => {
      if (!activeProjectId) return;
      updateProjectStatus(activeProjectId, 'failed');
      addMsg('assistant', `Pipeline failed at step "${stepId}": ${error}`);
    },
    [activeProjectId, updateProjectStatus, addMsg]
  );

  const handlePipelineCancel = useCallback(() => {
    if (!activeProjectId) return;
    setProjectPipeline(activeProjectId, null);
    updateProjectStatus(activeProjectId, 'active');
    addMsg('assistant', 'Pipeline cancelled.');
  }, [activeProjectId, setProjectPipeline, updateProjectStatus, addMsg]);

  const handleSelectDesign = useCallback(
    (pdbContent: string) => {
      setSelectedPdb(pdbContent);
    },
    [setSelectedPdb]
  );

  // ---- Build pipeline initial params from project context ----

  const pipelineInitialParams = useMemo(() => {
    if (!project) return {};
    const ctx = project.context;
    const params: Record<string, unknown> = {};
    if (ctx.nlQuery) params.user_prompt = ctx.nlQuery;
    if (ctx.pdbId) params.pdb_id = ctx.pdbId;
    if (ctx.pdbContent) params.pdb_content = ctx.pdbContent;
    if (ctx.targetMetal) params.target_metal = ctx.targetMetal;
    if (ctx.userPreferences) {
      Object.assign(params, ctx.userPreferences);
    }
    return params;
  }, [project]);

  // ---- Render ----

  const isEmpty = messages.length === 0 && !activePipelineId;

  // ---- Welcome state: Claude-style centered layout ----
  if (isEmpty && !showWizard) {
    return (
      <div className="h-full flex flex-col">
        <InitiativePanel
          onSuggestionClick={handleSend}
          onGuidedDesign={() => setShowWizard(true)}
          onSend={handleSend}
        />
      </div>
    );
  }

  // ---- Wizard state: header + centered wizard + bottom input ----
  if (showWizard && !activePipelineId) {
    return (
      <div className="h-full flex flex-col">
        {/* Minimal header */}
        <div className="border-b border-border px-6 py-3 shrink-0">
          <div className="flex items-center gap-3">
            <span className="bg-gradient-to-r from-primary to-primary/80 text-primary-foreground text-[11px] font-bold px-3 py-1 rounded-full uppercase tracking-wide flex items-center gap-1.5">
              <Sparkles className="h-4 w-4" />
              AI
            </span>
            <h2 className="text-lg font-bold text-foreground">Guided Design</h2>
          </div>
        </div>

        {/* Wizard centered */}
        <div className="flex-1 overflow-y-auto min-h-0">
          <div className="max-w-2xl mx-auto p-6">
            <GuidedDesignWizard
              onComplete={handleWizardComplete}
              onCancel={handleWizardCancel}
            />
          </div>
        </div>

        {/* Bottom input */}
        <div className="border-t border-border p-4 shrink-0">
          <ChatInput
            onSend={handleSend}
            placeholder="Or type your design request here to skip the wizard..."
          />
          <p className="mt-2 text-xs text-muted-foreground text-center">
            Complete the guided design above, or type to go directly
          </p>
        </div>
      </div>
    );
  }

  // ---- Active state: header + chat thread + pipeline + bottom input ----
  return (
    <div className="h-full flex flex-col">
      {/* Header */}
      <div className="border-b border-border px-6 py-3 shrink-0">
        <div className="flex items-center gap-3">
          <span className="bg-gradient-to-r from-primary to-primary/80 text-primary-foreground text-[11px] font-bold px-3 py-1 rounded-full uppercase tracking-wide flex items-center gap-1.5">
            <Sparkles className="h-4 w-4" />
            AI
          </span>
          <h2 className="text-lg font-bold text-foreground">
            {project?.name ?? 'AI Design Assistant'}
          </h2>
        </div>
      </div>

      {/* Outputs summary bar */}
      <OutputsSummary
        messages={messages}
        onSelectDesign={handleSelectDesign}
        onDeleteProject={activeProjectId ? () => deleteProject(activeProjectId) : undefined}
      />

      {/* Chat Thread */}
      <ChatThread
        messages={messages}
        onSelectDesign={handleSelectDesign}
        previousMessageCount={prevMsgCountRef.current}
      >
        {/* Active pipeline */}
        {activePipelineId && (
          <div className="pl-11">
            <PipelineCard
              key={`${activeProjectId}-${activePipelineId}`}
              pipelineId={activePipelineId}
              initialParams={pipelineInitialParams}
              onSelectDesign={handleSelectDesign}
              onStepComplete={handleStepComplete}
              onPipelineComplete={handlePipelineComplete}
              onPipelineFailed={handlePipelineFailed}
              onCancel={handlePipelineCancel}
            />
          </div>
        )}
      </ChatThread>

      {/* Input area */}
      <div className="border-t border-border p-4 shrink-0">
        <ChatInput
          onSend={handleSend}
          placeholder={
            activePipelineId
              ? 'Pipeline running — type a follow-up question...'
              : 'Describe your protein design (e.g., "Design a protein to bind citrate with terbium")...'
          }
        />
        <p className="mt-2 text-xs text-muted-foreground text-center">
          {activePipelineId
            ? 'Pipeline running — review each step above'
            : 'Natural language input or PDB code (e.g. 1BRF)'}
        </p>
      </div>
    </div>
  );
}
