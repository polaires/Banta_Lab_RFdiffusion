'use client';

import { useState, useCallback, useRef, useEffect, useMemo } from 'react';
import {
  Sparkles,
  Play,
  Gem,
  FlaskConical,
  Link,
  type LucideIcon,
} from 'lucide-react';
import { useStore } from '@/lib/store';
import type { ChatMessage as ChatMessageType, ChatAttachment } from '@/lib/chat-types';
import type { StepResult } from '@/lib/pipeline-types';
import type { UserPreferences } from '@/lib/api';
import { ChatThread } from './ChatThread';
import { ChatInput } from './ChatInput';
import { PipelineCard } from './PipelineCard';

// Interview components
import { InterviewMode } from '@/components/ai/InterviewMode';
import { MetalScaffoldInterviewMode } from '@/components/ai/MetalScaffoldInterviewMode';
import type { MetalScaffoldPreferences } from '@/components/ai/MetalScaffoldInterviewMode';
import { LigandInterviewMode } from '@/components/ai/LigandInterviewMode';
import type { LigandPreferences } from '@/components/ai/LigandInterviewMode';
import { BinderInterviewMode } from '@/components/ai/BinderInterviewMode';
import type { BinderPreferences } from '@/components/ai/BinderInterviewMode';

// ---- Design type detection from natural language ----

type InterviewType = 'metal_binding' | 'metal_scaffold' | 'ligand_dimer' | 'binder' | null;

/** Detect design type from user input to route to appropriate interview */
function detectDesignType(text: string): InterviewType {
  const lower = text.toLowerCase();

  // Metal scaffold: keywords like "scaffold", "de novo metal", "monomer scaffold"
  if (
    (lower.includes('scaffold') && (lower.includes('metal') || lower.includes('terbium') || lower.includes('lanthanide') || lower.includes('citrate'))) ||
    (lower.includes('de novo') && (lower.includes('metal') || lower.includes('binding')))
  ) {
    return 'metal_scaffold';
  }

  // Metal binding redesign: existing PDB + metal swap
  if (
    (lower.includes('redesign') || lower.includes('convert') || lower.includes('swap') || lower.includes('replace')) &&
    (lower.includes('metal') || lower.includes('iron') || lower.includes('zinc') || lower.includes('terbium') ||
     lower.includes('calcium') || lower.includes('europium') || lower.includes('gadolinium'))
  ) {
    return 'metal_binding';
  }

  // Ligand dimer: azobenzene, dimer, interface ligand
  if (
    lower.includes('azobenzene') || lower.includes('dimer') ||
    (lower.includes('interface') && lower.includes('ligand'))
  ) {
    return 'ligand_dimer';
  }

  // Binder design: binder, binding protein
  if (
    lower.includes('binder') ||
    (lower.includes('bind') && lower.includes('protein') && !lower.includes('metal'))
  ) {
    return 'binder';
  }

  return null;
}

// Icon mapping for quick-start cards
const QUICK_START_ICONS: Record<string, LucideIcon> = {
  token: Gem,
  science: FlaskConical,
  link: Link,
};

// ---- Quick Start Card ----

function QuickStartCard({
  icon,
  title,
  description,
  buttonText,
  onClick,
}: {
  icon: string;
  title: string;
  description: string;
  buttonText: string;
  onClick: () => void;
}) {
  const IconComponent = QUICK_START_ICONS[icon];
  return (
    <div className="p-5 bg-muted rounded-2xl border border-border">
      <div className="flex items-start gap-3">
        <div className="p-2.5 bg-card rounded-xl shadow-sm">
          {IconComponent && <IconComponent className="h-5 w-5 text-primary" />}
        </div>
        <div className="flex-1">
          <h3 className="font-bold text-foreground mb-1 text-sm">{title}</h3>
          <p className="text-xs text-muted-foreground mb-3">{description}</p>
          <button
            onClick={onClick}
            className="inline-flex items-center gap-1.5 px-3 py-1.5 bg-primary hover:bg-primary/90 text-primary-foreground rounded-lg font-medium text-xs transition-all"
          >
            <Play className="h-4 w-4" />
            {buttonText}
          </button>
        </div>
      </div>
    </div>
  );
}

// ---- Welcome message ----

function WelcomeSection({ onQuickStart }: { onQuickStart: (query: string) => void }) {
  return (
    <div className="space-y-4">
      {/* Welcome assistant message */}
      <div className="flex gap-3">
        <div className="flex-shrink-0 w-8 h-8 rounded-full bg-gradient-to-br from-primary to-primary/80 flex items-center justify-center">
          <Sparkles className="h-4 w-4 text-primary-foreground" />
        </div>
        <div className="flex-1">
          <div className="bg-muted rounded-2xl rounded-tl-sm px-4 py-3 text-sm text-foreground leading-relaxed">
            Welcome to the AI Design Assistant. Describe what you want to design in plain English, or pick a quick-start below.
          </div>
        </div>
      </div>

      {/* Quick start cards */}
      <div className="grid grid-cols-1 md:grid-cols-2 gap-3 pl-11">
        <QuickStartCard
          icon="token"
          title="Metal Binding Redesign"
          description="Convert iron-binding Rubredoxin (1BRF) to bind terbium."
          buttonText="Try 1BRF"
          onClick={() => onQuickStart('Redesign 1BRF to bind terbium instead of iron')}
        />
        <QuickStartCard
          icon="token"
          title="Metal Scaffold Design"
          description="Design monomer scaffolds with parameter sweep optimization."
          buttonText="Try METAL"
          onClick={() => onQuickStart('Design a metal binding scaffold for terbium with citrate')}
        />
        <QuickStartCard
          icon="science"
          title="Interface Dimer Design"
          description="Design separable dimer with azobenzene at the interface."
          buttonText="Try AZOB"
          onClick={() => onQuickStart('Design an azobenzene interface dimer')}
        />
        <QuickStartCard
          icon="link"
          title="Protein Binder Design"
          description="Design high-affinity binders with ESM-3 validation."
          buttonText="Try BIND"
          onClick={() => onQuickStart('Design a protein binder')}
        />
      </div>
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
  } = useStore();

  // Track previous message count for typing animation
  const prevMsgCountRef = useRef(0);

  // Interview state — which interview is currently active
  const [pendingInterview, setPendingInterview] = useState<InterviewType>(null);
  const interviewProjectIdRef = useRef<string | null>(null);
  const interviewQueryRef = useRef<string>('');

  // Get the active project
  const project = useMemo(
    () => projects.find((p) => p.id === activeProjectId) ?? null,
    [projects, activeProjectId]
  );

  const messages = project?.messages ?? [];
  const activePipelineId = project?.activePipelineId ?? null;

  // Update previous message count after render
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

  // ---- Handle user input ----

  const handleSend = useCallback(
    (text: string) => {
      let projId = activeProjectId;

      // If no active project, create one
      if (!projId) {
        projId = createProject();
      }

      // Add user message
      addProjectMessage(projId, { role: 'user', content: text });

      // Check if it looks like a PDB code (exactly 4 alphanum chars)
      const isPdbCode = /^[A-Za-z0-9]{4}$/.test(text.trim());

      if (isPdbCode) {
        // PDB code — add assistant message and start NL pipeline with resolve hint
        addProjectMessage(projId, {
          role: 'assistant',
          content: `Working with PDB ${text.trim().toUpperCase()}. Starting the design pipeline...`,
        });
        startPipeline(projId, 'natural-language', { user_prompt: text, pdb_id: text.trim().toUpperCase() });
        return;
      }

      // Detect design type for interview routing
      const designType = detectDesignType(text);

      if (designType && !activePipelineId) {
        // Show interview for specific design type
        const interviewLabels: Record<string, string> = {
          metal_binding: 'metal binding redesign',
          metal_scaffold: 'metal scaffold design',
          ligand_dimer: 'interface dimer design',
          binder: 'protein binder design',
        };
        addProjectMessage(projId, {
          role: 'assistant',
          content: `I'll guide you through ${interviewLabels[designType]} setup. Please configure your preferences below.`,
        });
        interviewProjectIdRef.current = projId;
        interviewQueryRef.current = text;
        setPendingInterview(designType);
      } else {
        // Generic NL query or follow-up during pipeline — start NL pipeline
        if (!activePipelineId) {
          addProjectMessage(projId, {
            role: 'assistant',
            content: 'Starting the design pipeline for your request...',
          });
          startPipeline(projId, 'natural-language', { user_prompt: text });
        } else {
          // Follow-up message during active pipeline — just add as user message
          // The pipeline is already running
        }
      }
    },
    [activeProjectId, addProjectMessage, createProject, activePipelineId]
  );

  // ---- Interview completion handlers ----

  const handleMetalBindingComplete = useCallback(
    (prefs: UserPreferences) => {
      const projId = interviewProjectIdRef.current;
      if (!projId) return;

      // Save preferences to project context
      updateProjectContext(projId, {
        targetMetal: prefs.targetMetal,
        designType: 'metal_binding',
        userPreferences: prefs as unknown as Record<string, unknown>,
      });

      // Add preference summary message
      addProjectMessage(projId, {
        role: 'assistant',
        content: 'Design preferences configured. Starting the pipeline...',
        attachment: {
          type: 'preference_summary',
          designType: 'metal_binding',
          preferences: {
            'Target Metal': prefs.targetMetalLabel,
            'Aggressiveness': prefs.aggressivenessLabel,
            'Coordination': prefs.coordinationLabel,
            'Designs': String(prefs.numDesigns),
            'Priority': prefs.priorityLabel,
          },
        },
      });

      // Start pipeline with preferences
      startPipelineRef.current(projId, 'natural-language', {
        ...prefs,
        user_prompt: interviewQueryRef.current,
        target_metal: prefs.targetMetal,
        num_designs: prefs.numDesigns,
      });

      setPendingInterview(null);
      interviewProjectIdRef.current = null;
    },
    [addProjectMessage, updateProjectContext]
  );

  const handleMetalScaffoldComplete = useCallback(
    (prefs: MetalScaffoldPreferences) => {
      const projId = interviewProjectIdRef.current;
      if (!projId) return;

      updateProjectContext(projId, {
        targetMetal: prefs.metal,
        designType: 'metal_scaffold',
        userPreferences: prefs as unknown as Record<string, unknown>,
      });

      addProjectMessage(projId, {
        role: 'assistant',
        content: `Scaffold design configured: ${prefs.numDesigns} total designs (${prefs.designsPerConfig}/config, ${prefs.optimizationModeLabel} mode). Starting pipeline...`,
        attachment: {
          type: 'preference_summary',
          designType: 'metal_scaffold',
          preferences: {
            'Metal': prefs.metalLabel,
            'Ligand': prefs.ligandLabel,
            'Scaffold Size': prefs.scaffoldSizeLabel,
            'Mode': prefs.optimizationModeLabel,
            'Designs/Config': String(prefs.designsPerConfig),
            'Total Designs': String(prefs.numDesigns),
          },
        },
      });

      startPipelineRef.current(projId, 'natural-language', {
        ...prefs,
        user_prompt: interviewQueryRef.current,
        target_metal: prefs.metal,
        ligand: prefs.ligand,
        contig_range: prefs.contigRange,
        num_designs: prefs.designsPerConfig,
        optimization_mode: prefs.optimizationMode,
        use_sweep: prefs.optimizationMode === 'sweep' || prefs.optimizationMode === 'production',
        designs_per_config: prefs.designsPerConfig,
      });

      setPendingInterview(null);
      interviewProjectIdRef.current = null;
    },
    [addProjectMessage, updateProjectContext]
  );

  const handleLigandComplete = useCallback(
    (prefs: LigandPreferences) => {
      const projId = interviewProjectIdRef.current;
      if (!projId) return;

      updateProjectContext(projId, {
        designType: 'ligand_dimer',
        userPreferences: prefs as unknown as Record<string, unknown>,
      });

      addProjectMessage(projId, {
        role: 'assistant',
        content: 'Dimer design configured. Starting the pipeline...',
        attachment: {
          type: 'preference_summary',
          designType: 'ligand_dimer',
          preferences: {
            'Ligand': prefs.ligandLabel,
            'Chain Length': prefs.chainLengthLabel,
            'Designs': String(prefs.numDesigns),
            'Priority': prefs.priorityLabel,
          },
        },
      });

      startPipelineRef.current(projId, 'natural-language', {
        ...prefs,
        user_prompt: interviewQueryRef.current,
        ligand_smiles: prefs.ligandSmiles,
        chain_length: prefs.chainLength,
        num_designs: prefs.numDesigns,
      });

      setPendingInterview(null);
      interviewProjectIdRef.current = null;
    },
    [addProjectMessage, updateProjectContext]
  );

  const handleBinderComplete = useCallback(
    (prefs: BinderPreferences) => {
      const projId = interviewProjectIdRef.current;
      if (!projId) return;

      updateProjectContext(projId, {
        pdbId: prefs.targetPdbId ?? null,
        designType: 'binder',
        userPreferences: prefs as unknown as Record<string, unknown>,
      });

      addProjectMessage(projId, {
        role: 'assistant',
        content: 'Binder design configured. Starting the pipeline...',
        attachment: {
          type: 'preference_summary',
          designType: 'binder',
          preferences: {
            'Target': prefs.targetLabel,
            'Protocol': prefs.protocolLabel,
            'Binder Size': prefs.binderLengthLabel,
            'Hotspots': prefs.hotspotMethodLabel,
            'Designs': String(prefs.numDesigns),
            'Quality': prefs.qualityThresholdLabel,
            'Priority': prefs.priorityLabel,
          },
        },
      });

      startPipelineRef.current(projId, 'natural-language', {
        ...prefs,
        user_prompt: interviewQueryRef.current,
        target_pdb_id: prefs.targetPdbId,
        binder_length: prefs.binderLength,
        num_designs: prefs.numDesigns,
        protocol: prefs.protocol,
      });

      setPendingInterview(null);
      interviewProjectIdRef.current = null;
    },
    [addProjectMessage, updateProjectContext]
  );

  const handleInterviewCancel = useCallback(() => {
    const projId = interviewProjectIdRef.current;
    if (projId) {
      addProjectMessage(projId, {
        role: 'assistant',
        content: 'Interview cancelled. You can describe your design in plain text and I\'ll handle the parameters automatically.',
      });
    }
    setPendingInterview(null);
    interviewProjectIdRef.current = null;
  }, [addProjectMessage]);

  // ---- Pipeline lifecycle ----

  // Use ref for startPipeline so interview handlers always see the latest
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

  // ---- Quick start ----

  const handleQuickStart = useCallback(
    (query: string) => {
      handleSend(query);
    },
    [handleSend]
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

  return (
    <div className="h-full flex flex-col">
      {/* Header */}
      <div className="border-b border-border px-6 py-4 shrink-0">
        <div className="flex items-center gap-3">
          <span className="bg-gradient-to-r from-primary to-primary/80 text-primary-foreground text-[11px] font-bold px-3 py-1 rounded-full uppercase tracking-wide flex items-center gap-1.5">
            <Sparkles className="h-4 w-4" />
            AI
          </span>
          <h2 className="text-lg font-bold text-foreground">
            {project?.name ?? 'AI Design Assistant'}
          </h2>
        </div>
        <p className="text-muted-foreground text-sm mt-1">
          Describe what you want to design in plain English. I&apos;ll handle the technical parameters.
        </p>
      </div>

      {/* Chat Thread */}
      <ChatThread
        messages={messages}
        onSelectDesign={handleSelectDesign}
        previousMessageCount={prevMsgCountRef.current}
      >
        {/* Welcome + quick-start when empty */}
        {isEmpty && !pendingInterview && <WelcomeSection onQuickStart={handleQuickStart} />}

        {/* Active interview */}
        {pendingInterview && (
          <div className="pl-11">
            {pendingInterview === 'metal_binding' && (
              <InterviewMode
                onComplete={handleMetalBindingComplete}
                onCancel={handleInterviewCancel}
              />
            )}
            {pendingInterview === 'metal_scaffold' && (
              <MetalScaffoldInterviewMode
                onComplete={handleMetalScaffoldComplete}
                onCancel={handleInterviewCancel}
              />
            )}
            {pendingInterview === 'ligand_dimer' && (
              <LigandInterviewMode
                onComplete={handleLigandComplete}
                onCancel={handleInterviewCancel}
              />
            )}
            {pendingInterview === 'binder' && (
              <BinderInterviewMode
                onComplete={handleBinderComplete}
                onCancel={handleInterviewCancel}
              />
            )}
          </div>
        )}

        {/* Active pipeline */}
        {activePipelineId && (
          <div className="pl-11">
            <PipelineCard
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

      {/* Input area — always enabled */}
      <div className="border-t border-border p-4 shrink-0">
        <ChatInput
          onSend={handleSend}
          placeholder={
            pendingInterview
              ? 'Complete the interview above, or type to skip and use automatic parameters...'
              : activePipelineId
              ? 'Pipeline running — type a follow-up question...'
              : 'Describe your protein design (e.g., "Design a protein to bind citrate with terbium")...'
          }
        />
        <p className="mt-2 text-xs text-muted-foreground text-center">
          {pendingInterview
            ? 'Configure your design preferences above'
            : activePipelineId
            ? 'Pipeline running — review each step above'
            : 'Natural language input or PDB code (e.g. 1BRF)'}
        </p>
      </div>
    </div>
  );
}
