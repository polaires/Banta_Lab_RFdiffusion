'use client';

import { useState, useRef, useEffect } from 'react';
import dynamic from 'next/dynamic';
import {
  Sparkles,
  User,
  Info,
  FlaskConical,
  Loader2,
  Play,
  Link,
  Box,
  RefreshCw,
  Gem,
  Search,
  Send,
  ListOrdered,
  CheckCircle,
  type LucideIcon,
} from 'lucide-react';

// Icon mapping for dynamic icons in QuickStartCard
const QUICK_START_ICONS: Record<string, LucideIcon> = {
  token: Gem,
  science: FlaskConical,
  link: Link,
};
import { useStore, type WorkflowPhase } from '@/lib/store';
import api, { type UserPreferences, type MetalBindingAnalysis, type DesignEvaluation } from '@/lib/api';
import { translatePreferencesToParams, explainParameters, type StructureInfo } from '@/lib/parameterTranslation';
import {
  useDesignWorkflow,
  createEvaluationFromDesign,
  METAL_CASE_STUDY,
  AZOBENZENE_CASE_STUDY,
  BINDER_CASE_STUDY,
} from '@/lib/useDesignWorkflow';
import { delay } from '@/lib/workflowHandlers';
import {
  InterviewMode,
  PreferenceSummaryCard,
  EvaluationCard,
  JobProgressCard,
  LigandInterviewMode,
  LigandAnalysisCard,
  LigandEvaluationCard,
  DesignGallery,
  BinderInterviewMode,
  BinderEvaluationCard,
  MetalScaffoldInterviewMode,
  MetalScaffoldDesignGallery,
  AIDesignPipelineWorkflow,
  type LigandPreferences,
  type LigandAnalysis,
  type LigandEvaluation,
  type DesignResult,
  type BinderPreferences,
  type BinderEvaluation,
  type MetalScaffoldPreferences,
  type MetalScaffoldDesignResult,
} from './ai';
import { useAIDesign } from '@/hooks/useAIDesign';
import {
  BinderResultsPanel,
  PipelineFunnel,
  HotspotSelector,
  type BinderDesign,
  type BinderStatistics,
} from './binder';

// Dynamic import of ProteinViewer to avoid SSR issues with Molstar
const ProteinViewer = dynamic(
  () => import('@/components/ProteinViewer').then((mod) => mod.ProteinViewer),
  {
    ssr: false,
    loading: () => (
      <div className="w-full h-[300px] bg-muted rounded-lg flex items-center justify-center">
        <div className="flex items-center gap-2 text-muted-foreground">
          <div className="w-5 h-5 border-2 border-muted-foreground/40 border-t-transparent rounded-full animate-spin" />
          Loading viewer...
        </div>
      </div>
    ),
  }
);

/**
 * AI Design Assistant Panel
 *
 * Demonstrates AI-guided protein engineering workflow with interview mode.
 * Users answer simple questions, and the system translates to RFD3 parameters.
 */

// Typing animation hook
function useTypingAnimation(text: string, speed: number = 20, enabled: boolean = true) {
  const [displayedText, setDisplayedText] = useState('');
  const [isComplete, setIsComplete] = useState(false);

  useEffect(() => {
    if (!enabled) {
      setDisplayedText(text);
      setIsComplete(true);
      return;
    }

    setDisplayedText('');
    setIsComplete(false);
    let index = 0;

    const timer = setInterval(() => {
      if (index < text.length) {
        setDisplayedText(text.slice(0, index + 1));
        index++;
      } else {
        setIsComplete(true);
        clearInterval(timer);
      }
    }, speed);

    return () => clearInterval(timer);
  }, [text, speed, enabled]);

  return { displayedText, isComplete };
}

// Message components
function AIMessage({
  content,
  isNew = false,
  children
}: {
  content: string;
  isNew?: boolean;
  children?: React.ReactNode;
}) {
  const { displayedText, isComplete } = useTypingAnimation(content, 15, isNew);

  return (
    <div className="flex gap-3">
      <div className="flex-shrink-0 w-8 h-8 rounded-full bg-gradient-to-br from-primary to-primary/80 flex items-center justify-center">
        <Sparkles className="h-4 w-4 text-primary-foreground" />
      </div>
      <div className="flex-1 space-y-3">
        <div className="bg-muted rounded-2xl rounded-tl-sm px-4 py-3 text-sm text-foreground leading-relaxed whitespace-pre-wrap">
          {displayedText}
          {!isComplete && <span className="animate-pulse">|</span>}
        </div>
        {isComplete && children}
      </div>
    </div>
  );
}

function UserMessage({ content }: { content: string }) {
  return (
    <div className="flex gap-3 justify-end">
      <div className="bg-primary text-primary-foreground rounded-2xl rounded-tr-sm px-4 py-3 text-sm max-w-[80%]">
        {content}
      </div>
      <div className="flex-shrink-0 w-8 h-8 rounded-full bg-muted flex items-center justify-center">
        <User className="h-4 w-4 text-muted-foreground" />
      </div>
    </div>
  );
}

function SystemMessage({ content }: { content: string }) {
  return (
    <div className="flex justify-center">
      <div className="bg-muted text-muted-foreground rounded-full px-4 py-1.5 text-xs flex items-center gap-2">
        <Info className="h-4 w-4" />
        {content}
      </div>
    </div>
  );
}

// Preference Summary Card for Ligand Design
function LigandPreferenceSummaryCard({
  preferences,
  onConfirm,
  onEdit,
  isRunning,
}: {
  preferences: LigandPreferences;
  onConfirm: () => void;
  onEdit: () => void;
  isRunning: boolean;
}) {
  return (
    <div className="bg-muted rounded-2xl p-6 border border-border shadow-sm">
      <div className="flex items-center gap-2 mb-4">
        <div className="w-8 h-8 rounded-lg bg-primary flex items-center justify-center">
          <FlaskConical className="h-4 w-4 text-primary-foreground" />
        </div>
        <h4 className="font-semibold text-foreground">Interface Dimer Configuration</h4>
      </div>

      <div className="space-y-2 mb-6">
        <PreferenceRow label="Ligand" value={preferences.ligandLabel} />
        <PreferenceRow label="Chain Length" value={preferences.chainLengthLabel} />
        <PreferenceRow label="Number of Designs" value={String(preferences.numDesigns)} />
        <PreferenceRow label="Priority" value={preferences.priorityLabel} isLast />
      </div>

      <div className="flex gap-3">
        <button
          onClick={onEdit}
          disabled={isRunning}
          className="flex-1 px-4 py-2.5 text-primary bg-card border border-border rounded-xl font-medium text-sm hover:bg-muted transition-all disabled:opacity-50"
        >
          Edit
        </button>
        <button
          onClick={onConfirm}
          disabled={isRunning}
          className="flex-1 px-4 py-2.5 bg-primary text-primary-foreground rounded-xl font-medium text-sm hover:bg-primary/90 transition-all disabled:opacity-50 flex items-center justify-center gap-2"
        >
          {isRunning ? (
            <>
              <Loader2 className="h-4 w-4 animate-spin" />
              Running...
            </>
          ) : (
            <>
              <Play className="h-4 w-4" />
              Run Design
            </>
          )}
        </button>
      </div>
    </div>
  );
}

// Preference Summary Card for Binder Design
function BinderPreferenceSummaryCard({
  preferences,
  onConfirm,
  onEdit,
  isRunning,
}: {
  preferences: BinderPreferences;
  onConfirm: () => void;
  onEdit: () => void;
  isRunning: boolean;
}) {
  return (
    <div className="bg-muted rounded-2xl p-6 border border-border shadow-sm">
      <div className="flex items-center gap-2 mb-4">
        <div className="w-8 h-8 rounded-lg bg-primary flex items-center justify-center">
          <Link className="h-4 w-4 text-primary-foreground" />
        </div>
        <h4 className="font-semibold text-foreground">Protein Binder Configuration</h4>
      </div>

      <div className="space-y-2 mb-6">
        <PreferenceRow label="Target Protein" value={preferences.targetLabel} />
        <PreferenceRow label="Binder Size" value={preferences.binderLengthLabel} />
        <PreferenceRow label="Number of Designs" value={String(preferences.numDesigns)} />
        <PreferenceRow label="Quality Threshold" value={preferences.qualityThresholdLabel} />
        <PreferenceRow label="Priority" value={preferences.priorityLabel} isLast />
      </div>

      <div className="flex gap-3">
        <button
          onClick={onEdit}
          disabled={isRunning}
          className="flex-1 px-4 py-2.5 text-primary bg-card border border-border rounded-xl font-medium text-sm hover:bg-muted transition-all disabled:opacity-50"
        >
          Edit
        </button>
        <button
          onClick={onConfirm}
          disabled={isRunning}
          className="flex-1 px-4 py-2.5 bg-primary text-primary-foreground rounded-xl font-medium text-sm hover:bg-primary/90 transition-all disabled:opacity-50 flex items-center justify-center gap-2"
        >
          {isRunning ? (
            <>
              <Loader2 className="h-4 w-4 animate-spin" />
              Running...
            </>
          ) : (
            <>
              <Play className="h-4 w-4" />
              Run Design
            </>
          )}
        </button>
      </div>
    </div>
  );
}

// Preference Summary Card for Metal Scaffold Design
function MetalScaffoldPreferenceSummaryCard({
  preferences,
  onConfirm,
  onEdit,
  isRunning,
}: {
  preferences: MetalScaffoldPreferences;
  onConfirm: () => void;
  onEdit: () => void;
  isRunning: boolean;
}) {
  return (
    <div className="bg-muted rounded-2xl p-6 border border-border shadow-sm">
      <div className="flex items-center gap-2 mb-4">
        <div className="w-8 h-8 rounded-lg bg-primary flex items-center justify-center">
          <Gem className="h-4 w-4 text-primary-foreground" />
        </div>
        <h4 className="font-semibold text-foreground">Metal Scaffold Configuration</h4>
      </div>

      <div className="space-y-2 mb-6">
        <PreferenceRow label="Target Metal" value={preferences.metalLabel} />
        <PreferenceRow label="Ligand" value={preferences.ligandLabel} />
        <PreferenceRow label="Scaffold Size" value={`${preferences.scaffoldSizeLabel} (${preferences.contigRange} aa)`} />
        <PreferenceRow label="Optimization" value={preferences.optimizationModeLabel} />
        <PreferenceRow label="Total Designs" value={String(preferences.numDesigns)} isLast />
      </div>

      {preferences.optimizationMode === 'sweep' && (
        <div className="mb-4 p-3 bg-card rounded-lg border border-border">
          <div className="text-xs text-muted-foreground mb-1">Parameter Sweep Details</div>
          <div className="text-sm text-foreground">
            9 configs (3 sizes × 3 CFG) × {preferences.designsPerConfig} designs = {preferences.numDesigns} total
          </div>
        </div>
      )}

      <div className="flex gap-3">
        <button
          onClick={onEdit}
          disabled={isRunning}
          className="flex-1 px-4 py-2.5 text-primary bg-card border border-border rounded-xl font-medium text-sm hover:bg-muted transition-all disabled:opacity-50"
        >
          Edit
        </button>
        <button
          onClick={onConfirm}
          disabled={isRunning}
          className="flex-1 px-4 py-2.5 bg-primary text-primary-foreground rounded-xl font-medium text-sm hover:bg-primary/90 transition-all disabled:opacity-50 flex items-center justify-center gap-2"
        >
          {isRunning ? (
            <>
              <Loader2 className="h-4 w-4 animate-spin" />
              Running...
            </>
          ) : (
            <>
              <Play className="h-4 w-4" />
              Run Design
            </>
          )}
        </button>
      </div>
    </div>
  );
}

// Shared preference row component
function PreferenceRow({ label, value, isLast = false }: { label: string; value: string; isLast?: boolean }) {
  return (
    <div className={`flex justify-between items-center py-2 ${!isLast ? 'border-b border-border' : ''}`}>
      <span className="text-sm text-muted-foreground">{label}</span>
      <span className="font-medium text-foreground">{value}</span>
    </div>
  );
}

// Structure Viewer section
function StructureViewerSection({
  pdbContent,
  title,
  badge,
  colorLegend,
  emptyMessage,
}: {
  pdbContent: string | null;
  title: string;
  badge?: string;
  colorLegend?: React.ReactNode;
  emptyMessage?: string;
}) {
  return (
    <div className="bg-card rounded-xl border border-border overflow-hidden">
      <div className="bg-muted px-4 py-2 border-b border-border flex items-center justify-between">
        <div className="flex items-center gap-2">
          <Box className="h-5 w-5 text-primary" />
          <h4 className="font-semibold text-foreground text-sm">{title}</h4>
          {badge && (
            <span className="text-xs bg-green-100 text-green-700 px-2 py-0.5 rounded-full">
              {badge}
            </span>
          )}
        </div>
        {colorLegend}
      </div>
      <div className="relative">
        {pdbContent ? (
          <ProteinViewer pdbContent={pdbContent} className="h-[350px]" />
        ) : (
          <div className="h-[350px] bg-gradient-to-br from-muted to-background flex flex-col items-center justify-center text-muted-foreground">
            <Box className="h-10 w-10 mb-2 text-muted-foreground/60" />
            <p className="text-sm font-medium">{emptyMessage || 'No Structure Available'}</p>
            <p className="text-xs text-muted-foreground/60 mt-1">Connect to backend to view actual designs</p>
          </div>
        )}
      </div>
    </div>
  );
}

export function AIDesignAssistantPanel() {
  const {
    aiCaseStudy,
    setAiCaseStudy,
    addAiMessage,
    clearAiConversation,
    setSelectedPdb,
    setLatestRfd3Design,
    backendUrl,
  } = useStore();

  // Use the consolidated workflow hook
  const { state: workflowState, isBackendConnected, executeMetal, executeLigand, executeBinder } = useDesignWorkflow();

  const [pdbInput, setPdbInput] = useState('');
  const [followUpInput, setFollowUpInput] = useState('');
  const [userPreferences, setUserPreferences] = useState<UserPreferences | null>(null);
  const [analysisResult, setAnalysisResult] = useState<MetalBindingAnalysis | null>(null);
  const [evaluationResult, setEvaluationResult] = useState<DesignEvaluation | null>(null);
  const [structureInfo, setStructureInfo] = useState<StructureInfo | null>(null);
  const messagesEndRef = useRef<HTMLDivElement>(null);
  const inputRef = useRef<HTMLInputElement>(null);

  // Ligand interface design state
  const [isLigandDesign, setIsLigandDesign] = useState(false);
  const [ligandAnalysisResult, setLigandAnalysisResult] = useState<LigandAnalysis | null>(null);
  const [ligandPreferences, setLigandPreferences] = useState<LigandPreferences | null>(null);
  const [ligandEvaluationResult, setLigandEvaluationResult] = useState<LigandEvaluation | null>(null);

  // Protein binder design state
  const [isBinderDesign, setIsBinderDesign] = useState(false);
  const [binderPreferences, setBinderPreferences] = useState<BinderPreferences | null>(null);
  const [binderEvaluationResult, setBinderEvaluationResult] = useState<BinderEvaluation | null>(null);
  const [binderDesigns, setBinderDesigns] = useState<BinderDesign[]>([]);
  const [binderStatistics, setBinderStatistics] = useState<BinderStatistics | null>(null);
  const [selectedBinderDesign, setSelectedBinderDesign] = useState<BinderDesign | null>(null);
  const [showHotspotSelector, setShowHotspotSelector] = useState(false);
  const [hotspotPdbContent, setHotspotPdbContent] = useState<string | null>(null);

  // Metal scaffold design state (Round 7b monomer workflow)
  const [isMetalScaffoldDesign, setIsMetalScaffoldDesign] = useState(false);
  const [metalScaffoldPreferences, setMetalScaffoldPreferences] = useState<MetalScaffoldPreferences | null>(null);
  const [metalScaffoldSweepResults, setMetalScaffoldSweepResults] = useState<{
    configs: Array<{
      name: string;
      passRate: number;
      avgPlddt: number;
      avgPtm: number;
      tierCounts: Record<string, number>;
    }>;
    bestConfig: string | null;
    totalDesigns: number;
    passingDesigns: number;
  } | null>(null);
  const [metalScaffoldDesigns, setMetalScaffoldDesigns] = useState<MetalScaffoldDesignResult[]>([]);
  const [selectedMetalDesignId, setSelectedMetalDesignId] = useState<string | null>(null);

  // Multiple design results state
  const [designResults, setDesignResults] = useState<DesignResult[]>([]);
  const [selectedDesignId, setSelectedDesignId] = useState<string | null>(null);

  // Viewer state
  const [viewerPdbContent, setViewerPdbContent] = useState<string | null>(null);

  // Natural Language AI Design state
  const [isNLDesign, setIsNLDesign] = useState(false);

  // Determine API URL based on backend URL
  // If using local Docker (localhost:8000), use traditional proxy
  // If using RunPod cloud, use runpod proxy
  const isLocalBackend = backendUrl.includes('localhost') || backendUrl.includes('127.0.0.1');
  const aiDesignApiUrl = isLocalBackend
    ? `/api/traditional/runsync?url=${encodeURIComponent(backendUrl)}`
    : '/api/runpod/runsync';

  const {
    stage: nlStage,
    stageInfo: nlStageInfo,
    result: nlResult,
    error: nlError,
    isRunning: nlIsRunning,
    runAIDesign,
    reset: resetNLDesign
  } = useAIDesign({
    apiUrl: aiDesignApiUrl,
    onStageChange: (stage) => {
      // Update conversation with stage changes
      if (stage === 'parsing') {
        addAiMessage({ role: 'assistant', content: 'Understanding your design request...' });
      } else if (stage === 'resolving') {
        addAiMessage({ role: 'assistant', content: 'Resolving ligand structure and chemistry...' });
      } else if (stage === 'backbone') {
        addAiMessage({ role: 'assistant', content: 'Generating backbone structures with RFD3...' });
      } else if (stage === 'sequence') {
        addAiMessage({ role: 'assistant', content: 'Designing sequences with LigandMPNN...' });
      } else if (stage === 'validation') {
        addAiMessage({ role: 'assistant', content: 'Validating structures with ESMFold...' });
      }
    }
  });

  const workflowPhase = (aiCaseStudy.workflowPhase || 'idle') as WorkflowPhase;
  const isRunningJob = workflowState.isRunning;

  // Auto-scroll to bottom when new messages arrive
  useEffect(() => {
    messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' });
  }, [aiCaseStudy.conversation, workflowPhase]);

  // Check if input is a natural language query vs PDB code
  const isNaturalLanguageQuery = (input: string): boolean => {
    const trimmed = input.trim();
    // PDB codes are exactly 4 alphanumeric chars, or special demo codes
    if (/^[A-Za-z0-9]{4}$/.test(trimmed)) return false;
    if (['AZOB', 'BIND', 'METAL'].includes(trimmed.toUpperCase())) return false;
    // If longer than 6 chars or contains spaces/keywords, it's NL
    if (trimmed.length > 6) return true;
    const nlKeywords = ['design', 'create', 'make', 'build', 'bind', 'protein', 'with', 'for', 'want'];
    return nlKeywords.some(kw => trimmed.toLowerCase().includes(kw));
  };

  // Handle natural language design query
  const handleNLDesign = async (query: string) => {
    setIsNLDesign(true);
    setPdbInput('');
    addAiMessage({ role: 'user', content: query });
    setAiCaseStudy({ workflowPhase: 'running', isProcessing: true });

    const result = await runAIDesign(query, {
      numDesigns: 4,
      numSequences: 8,
      validate: true
    });

    if (result?.success) {
      setAiCaseStudy({ workflowPhase: 'complete', isProcessing: false });
      addAiMessage({
        role: 'assistant',
        content: `Design complete!\n\n**Results:**\n- Backbones generated: ${result.num_backbones}\n- Sequences designed: ${result.num_sequences}\n- Pass rate: ${(result.pass_rate * 100).toFixed(0)}%\n${result.best_plddt ? `- Best pLDDT: ${result.best_plddt.toFixed(2)}` : ''}\n\n${result.recommendations?.length ? '**Recommendations:**\n' + result.recommendations.map(r => `- ${r}`).join('\n') : ''}`
      });
      if (result.best_sequence_pdb) {
        setSelectedPdb(result.best_sequence_pdb);
      }
    } else {
      setAiCaseStudy({ workflowPhase: 'idle', isProcessing: false });
      addAiMessage({
        role: 'assistant',
        content: `Design failed: ${nlError || 'Unknown error'}. Please try again with a different query.`
      });
    }
  };

  // Handle structure input (PDB code or natural language)
  const handleStructureInput = async () => {
    if (!pdbInput.trim()) return;

    const input = pdbInput.trim();

    // Check if it's a natural language query
    if (isNaturalLanguageQuery(input)) {
      await handleNLDesign(input);
      return;
    }

    const pdbId = input.toUpperCase();
    setPdbInput('');

    addAiMessage({ role: 'user', content: `I want to work with ${pdbId}` });
    setAiCaseStudy({ workflowPhase: 'fetching', isProcessing: true });

    await delay(500);
    addAiMessage({ role: 'assistant', content: `Fetching ${pdbId} from the Protein Data Bank...` });
    await delay(1500);

    // Handle special demo codes
    if (pdbId === 'AZOB') {
      await handleAzobenzeneDemo();
      return;
    }

    if (pdbId === 'BIND') {
      await handleBinderDemo();
      return;
    }

    if (pdbId === 'METAL') {
      await handleMetalScaffoldDemo();
      return;
    }

    // Handle metal binding demos and real PDB codes
    await handleMetalBindingInput(pdbId);
  };

  // Handle azobenzene interface dimer demo
  async function handleAzobenzeneDemo() {
    setIsLigandDesign(true);
    setAiCaseStudy({ pdbId: 'AZOB', workflowPhase: 'analyzing' });

    addAiMessage({
      role: 'assistant',
      content: `Great! You want to design an azobenzene interface dimer.\n\n- Ligand: Azobenzene (c1ccc(cc1)N=Nc2ccccc2)\n- Design Type: Interface Dimer\n- Approach: Two-step asymmetric-sequential\n\nAnalyzing ligand topology...`,
    });

    await delay(2000);

    const analysis = AZOBENZENE_CASE_STUDY.ligandAnalysis;
    setLigandAnalysisResult(analysis);
    setAiCaseStudy({ workflowPhase: 'interview', isProcessing: false });

    addAiMessage({
      role: 'assistant',
      content: `Ligand analysis complete!\n\n- Molecular Weight: ${analysis.ligand.molecular_weight} Da\n- Heavy Atoms: ${analysis.ligand.num_heavy_atoms}\n- Topology: ${analysis.topology.separable ? 'Separable' : 'Entangled'}\n\nThis ligand is ideal for interface dimer design. Let me ask you a few questions about your design preferences...`,
    });
  }

  // Handle protein binder design demo
  async function handleBinderDemo() {
    setIsBinderDesign(true);
    setAiCaseStudy({ pdbId: 'BIND', workflowPhase: 'analyzing' });

    addAiMessage({
      role: 'assistant',
      content: `Great! You want to design a protein binder.\n\n- Pipeline: Multi-stage validation with ESM-3 + MPNN\n- Design Type: De novo protein binder\n- Validation: Interface contacts, H-bonds, packing\n\nAnalyzing pipeline capabilities...`,
    });

    await delay(2000);
    setAiCaseStudy({ workflowPhase: 'interview', isProcessing: false });

    addAiMessage({
      role: 'assistant',
      content: `Pipeline ready!\n\n- Structure Generation: RFdiffusion 3\n- Sequence Design: ProteinMPNN\n- Sequence Scoring: ESM-3 confidence\n- Interface Analysis: Contact counting, H-bonds\n\nThis pipeline achieves ~50% binder success rate. Let me ask you a few questions about your design preferences...`,
    });
  }

  // Handle metal scaffold design (Round 7b monomer workflow)
  async function handleMetalScaffoldDemo() {
    setIsMetalScaffoldDesign(true);
    setAiCaseStudy({ pdbId: 'METAL', workflowPhase: 'analyzing' });

    addAiMessage({
      role: 'assistant',
      content: `Excellent choice! You want to design a metal binding scaffold.\n\n- Design Type: **Monomer** (single chain) scaffold\n- Pipeline: RFD3 → LigandMPNN → RF3 validation\n- Optimization: Parameter sweep (3 sizes × 3 CFG scales)\n- Filtering: Quality tiers (S/A/B/C/F)\n\nPreparing the design workflow...`,
    });

    await delay(2000);
    setAiCaseStudy({ workflowPhase: 'interview', isProcessing: false });

    addAiMessage({
      role: 'assistant',
      content: `Round 7b-style workflow ready!\n\n**Key Features:**\n- **Monomer focus**: Single-chain scaffolds (100-150 residues)\n- **Parameter sweep**: Test 9 configurations automatically\n- **Quality tiers**: S-tier (best) → F-tier (fail)\n- **HSAB compliance**: Hard acid donors for lanthanides\n\nLet me ask you a few questions to configure your metal binding scaffold design...`,
    });
  }

  // Handle metal binding input (demo or real)
  async function handleMetalBindingInput(pdbId: string) {
    if (!isBackendConnected) {
      // Demo mode
      if (pdbId === '1BRF') {
        await handleMetalDemo();
      } else {
        addAiMessage({
          role: 'assistant',
          content: `In demo mode, you can try:\n\n- "1BRF" - Metal binding site redesign (Iron to Terbium)\n- "AZOB" - Azobenzene interface dimer design\n- "BIND" - Protein binder design\n\nFor other structures, please connect to the backend.`,
        });
        setAiCaseStudy({ workflowPhase: 'structure_input', isProcessing: false });
      }
      return;
    }

    // Backend connected - fetch from RCSB
    try {
      const pdbResult = await api.fetchPdb(pdbId);
      if (!pdbResult.content) {
        throw new Error('No PDB content received');
      }

      setAiCaseStudy({ pdbId, pdbContent: pdbResult.content, workflowPhase: 'analyzing' });

      const info = pdbResult.info;
      if (info) {
        setStructureInfo({
          chains: info.chains,
          num_residues: info.num_residues,
          num_atoms: info.num_atoms,
        });
      }

      addAiMessage({
        role: 'assistant',
        content: `Got it! ${info?.title || pdbId}\n\n- ${info?.num_residues || '?'} residues\n- ${info?.num_atoms || '?'} atoms\n${info?.metals?.length ? `- Metal found: ${info.metals[0].element}` : '- No metal detected'}\n\nNow analyzing the structure...`,
      });

      await delay(1500);
      await analyzeMetalBinding(pdbResult);
    } catch (error) {
      addAiMessage({
        role: 'assistant',
        content: `Error: ${error instanceof Error ? error.message : 'Unknown error'}\n\nPlease check the PDB code and try again, or use the demo with "1BRF".`,
      });
      setAiCaseStudy({ workflowPhase: 'structure_input', isProcessing: false });
    }
  }

  // Handle metal binding demo (1BRF)
  async function handleMetalDemo() {
    const caseStudy = METAL_CASE_STUDY;
    setAiCaseStudy({
      pdbId: '1BRF',
      pdbContent: 'simulated_content',
      workflowPhase: 'analyzing',
    });

    setStructureInfo({
      chains: caseStudy.pdbInfo.chains,
      num_residues: caseStudy.pdbInfo.num_residues,
      num_atoms: caseStudy.pdbInfo.num_atoms,
    });

    addAiMessage({
      role: 'assistant',
      content: `Got it! ${caseStudy.pdbInfo.title}\n\n- ${caseStudy.pdbInfo.num_residues} residues\n- Metal found: Iron (Fe)\n\nNow let me analyze the metal binding site...`,
    });

    await delay(2000);

    setAnalysisResult(caseStudy.metalAnalysis);
    setAiCaseStudy({
      analysisResult: caseStudy.metalAnalysis,
      workflowPhase: 'interview',
      isProcessing: false,
    });

    const analysis = caseStudy.metalAnalysis;
    addAiMessage({
      role: 'assistant',
      content: `I found a ${analysis.metal.element} binding site with ${analysis.coordination.number} coordinating atoms.\n\nThe current geometry is ${analysis.coordination.geometry} with primarily ${analysis.donor_analysis.dominant_type} donors.\n\nNow let me ask you a few questions to understand how you'd like to redesign this site...`,
    });
  }

  // Analyze metal binding from PDB result
  async function analyzeMetalBinding(pdbResult: Awaited<ReturnType<typeof api.fetchPdb>>) {
    const metalInfo = pdbResult.info?.metals?.[0];

    if (!metalInfo) {
      addAiMessage({
        role: 'assistant',
        content: `I couldn't find any metal ions in this structure. The AI assistant works best with metal-binding proteins.\n\nYou can still proceed to design a new metal binding site. Let me ask you some questions...`,
      });

      const minimalAnalysis: MetalBindingAnalysis = {
        success: true,
        metal: { element: 'NONE', chain: 'A', residue: 'UNK', resnum: 0, position: [0, 0, 0] },
        coordination: { number: 0, geometry: 'none', geometry_rmsd: 0, coordinating_atoms: [] },
        donor_analysis: { types: {}, dominant_type: 'none', lanthanide_compatible: false },
        bond_analysis: { average_distance: 0, min_distance: 0, max_distance: 0, distances: [] },
      };
      setAnalysisResult(minimalAnalysis);
      setAiCaseStudy({ analysisResult: minimalAnalysis, workflowPhase: 'interview', isProcessing: false });
      return;
    }

    // Create frontend-based analysis
    const frontendAnalysis: MetalBindingAnalysis = {
      success: true,
      metal: {
        element: metalInfo.element,
        chain: metalInfo.chain,
        residue: metalInfo.residue,
        resnum: parseInt(metalInfo.res_num),
        position: [0, 0, 0],
      },
      coordination: {
        number: 4,
        geometry: 'tetrahedral',
        geometry_rmsd: 0.2,
        coordinating_atoms: [],
      },
      donor_analysis: {
        types: { 'unknown': 4 },
        dominant_type: 'unknown',
        lanthanide_compatible: false,
      },
      bond_analysis: {
        average_distance: 2.3,
        min_distance: 2.1,
        max_distance: 2.5,
        distances: [],
      },
      suggestions: ['Detailed analysis requires backend connection with analysis endpoints.'],
    };

    setAnalysisResult(frontendAnalysis);
    setAiCaseStudy({ analysisResult: frontendAnalysis, workflowPhase: 'interview', isProcessing: false });

    addAiMessage({
      role: 'assistant',
      content: `I found a ${metalInfo.element} binding site in chain ${metalInfo.chain}.\n\nNow let me ask you a few questions to understand how you'd like to redesign this site...`,
    });
  }

  // Handle interview completion (metal binding)
  const handleInterviewComplete = async (prefs: UserPreferences) => {
    setUserPreferences(prefs);
    setAiCaseStudy({
      userPreferences: prefs,
      targetMetal: prefs.targetMetal,
      workflowPhase: 'confirming',
    });

    addAiMessage({
      role: 'assistant',
      content: `Great choices! Here's a summary of your preferences:\n\n- Target Metal: ${prefs.targetMetalLabel}\n- Design Approach: ${prefs.aggressivenessLabel}\n- Coordination: ${prefs.coordinationLabel}\n- Variants: ${prefs.numDesigns} designs\n- Priority: ${prefs.priorityLabel}\n\nReview and click "Run Design" when ready!`,
    });
  };

  // Handle ligand interview completion
  const handleLigandInterviewComplete = async (prefs: LigandPreferences) => {
    setLigandPreferences(prefs);
    setAiCaseStudy({ workflowPhase: 'confirming' });

    addAiMessage({
      role: 'assistant',
      content: `Excellent! Here's your interface dimer design configuration:\n\n- Ligand: ${prefs.ligandLabel} (${prefs.ligandSmiles})\n- Chain Length: ${prefs.chainLengthLabel}\n- Number of Designs: ${prefs.numDesigns}\n- Priority: ${prefs.priorityLabel}\n\nThis will use the two-step workflow:\n1. Design Chain A (one-sided binding)\n2. Design Chain B (complementary binding)\n\nClick "Run Design" when ready!`,
    });
  };

  // Handle metal scaffold interview completion
  const handleMetalScaffoldInterviewComplete = async (prefs: MetalScaffoldPreferences) => {
    setMetalScaffoldPreferences(prefs);
    setAiCaseStudy({ workflowPhase: 'confirming' });

    const sweepInfo = prefs.optimizationMode === 'sweep'
      ? `\n\n**Parameter Sweep:**\n- 3 scaffold sizes × 3 CFG scales = 9 configurations\n- ${prefs.designsPerConfig} designs per config = ${prefs.numDesigns} total designs`
      : prefs.optimizationMode === 'quick'
      ? `\n\n**Quick Exploration:**\n- Single configuration: ${prefs.scaffoldSizeLabel}\n- ${prefs.numDesigns} designs total`
      : `\n\n**Production Run:**\n- Use best config from previous sweep\n- ${prefs.numDesigns} designs total`;

    addAiMessage({
      role: 'assistant',
      content: `Here's your metal scaffold design configuration:\n\n- **Metal**: ${prefs.metalLabel}\n- **Ligand**: ${prefs.ligandLabel}\n- **Scaffold Size**: ${prefs.scaffoldSizeLabel} (${prefs.contigRange} residues)\n- **Optimization**: ${prefs.optimizationModeLabel}${sweepInfo}\n\nThe pipeline will:\n1. Generate backbone scaffolds (RFD3 with RASA + H-bond conditioning)\n2. Design sequences (LigandMPNN with HSAB biases)\n3. Validate structures (RF3)\n4. Filter by quality tiers (S/A/B/C/F)\n\nClick "Run Design" when ready!`,
    });
  };

  // Handle binder interview completion
  const handleBinderInterviewComplete = async (prefs: BinderPreferences) => {
    setBinderPreferences(prefs);

    // If manual hotspot selection, fetch PDB and show selector
    if (prefs.hotspotMethod === 'manual' && prefs.targetPdbId) {
      addAiMessage({
        role: 'assistant',
        content: `Let me fetch the ${prefs.targetLabel} structure so you can select binding hotspots...`,
      });

      try {
        const pdbResult = await api.fetchPdb(prefs.targetPdbId);
        if (pdbResult.content) {
          setHotspotPdbContent(pdbResult.content);
          setShowHotspotSelector(true);

          addAiMessage({
            role: 'assistant',
            content: `Structure loaded! Use the 3D viewer to select residues where the binder should attach.\n\n**Tip:** Click "Detect" to let AI suggest optimal hotspots based on surface accessibility, then adjust as needed.`,
          });
        } else {
          throw new Error('Could not fetch PDB content');
        }
      } catch (err) {
        console.error('[AIPanel] Failed to fetch PDB for hotspot selection:', err);
        addAiMessage({
          role: 'assistant',
          content: `Could not fetch the target structure. Switching to auto-detect mode for hotspots.`,
        });
        const updatedPrefs = { ...prefs, hotspotMethod: 'auto' as const };
        setBinderPreferences(updatedPrefs);
        proceedToConfirmation(updatedPrefs);
      }
    } else {
      proceedToConfirmation(prefs);
    }
  };

  // Helper to proceed to confirmation phase
  function proceedToConfirmation(prefs: BinderPreferences) {
    setAiCaseStudy({ workflowPhase: 'confirming' });

    const hotspotInfo = prefs.hotspotMethod === 'manual' && prefs.manualHotspots?.length
      ? `\n- Hotspots: ${prefs.manualHotspots.join(', ')}`
      : prefs.hotspotMethod === 'auto'
      ? '\n- Hotspots: Auto-detect'
      : '';

    addAiMessage({
      role: 'assistant',
      content: `Excellent! Here's your protein binder design configuration:\n\n- Target: ${prefs.targetLabel}\n- Binder Size: ${prefs.binderLengthLabel}\n- Number of Designs: ${prefs.numDesigns}\n- Quality Threshold: ${prefs.qualityThresholdLabel}\n- Priority: ${prefs.priorityLabel}${hotspotInfo}\n\nThe pipeline will:\n1. Generate backbone structures (RFdiffusion)\n2. Design sequences (ProteinMPNN)\n3. Score sequences (ESM-3)\n4. Analyze interfaces\n5. Filter and rank results\n\nClick "Run Design" when ready!`,
    });
  }

  // Handle hotspot selection
  const handleHotspotConfirm = (hotspots: string[]) => {
    setShowHotspotSelector(false);

    if (binderPreferences) {
      const updatedPrefs = {
        ...binderPreferences,
        manualHotspots: hotspots,
        targetPdbContent: hotspotPdbContent || undefined,
      };
      setBinderPreferences(updatedPrefs);
      proceedToConfirmation(updatedPrefs);
    }
  };

  const handleHotspotCancel = () => {
    setShowHotspotSelector(false);

    if (binderPreferences) {
      const updatedPrefs = { ...binderPreferences, hotspotMethod: 'auto' as const };
      setBinderPreferences(updatedPrefs);
      proceedToConfirmation(updatedPrefs);
      addAiMessage({ role: 'assistant', content: `Switched to auto-detect mode for hotspots.` });
    }
  };

  // Run metal binding design
  const handleRunDesign = async () => {
    if (!userPreferences) return;

    setAiCaseStudy({ workflowPhase: 'running' });

    const params = translatePreferencesToParams(userPreferences, analysisResult, structureInfo);
    const explanations = explainParameters(params, userPreferences);

    addAiMessage({
      role: 'assistant',
      content: `Starting the design process with your preferences:\n\n${explanations.map(e => `- ${e}`).join('\n')}\n\nThis may take a few minutes...`,
    });

    const result = await executeMetal(params, aiCaseStudy.pdbContent || '', {
      onJobCreated: (jobId) => setAiCaseStudy({ pendingJobId: jobId }),
      onProgress: (progress) => setAiCaseStudy({ jobProgress: progress }),
      onError: (error) => addAiMessage({ role: 'assistant', content: error }),
    });

    if (result) {
      if (result.pdbContent) {
        setSelectedPdb(result.pdbContent);
        setViewerPdbContent(result.pdbContent);
        setLatestRfd3Design({
          jobId: aiCaseStudy.pendingJobId || 'demo',
          pdbContent: result.pdbContent,
          source: 'rfd3',
          timestamp: Date.now(),
        });
      }

      setEvaluationResult(result.evaluation);
      setAiCaseStudy({
        evaluationResult: result.evaluation,
        workflowPhase: 'complete',
      });

      const passCount = result.evaluation.criteria_passed;
      const message = result.evaluation.overall_pass
        ? `Excellent results! The design meets ${passCount}/${result.evaluation.criteria_total} target criteria for ${userPreferences.targetMetalLabel} binding. The coordination geometry looks ideal.`
        : `The design meets ${passCount}/${result.evaluation.criteria_total} criteria for ${userPreferences.targetMetalLabel} binding. Consider adjusting settings.`;

      addAiMessage({ role: 'assistant', content: message });
    }
  };

  // Run ligand design
  const handleRunLigandDesign = async () => {
    if (!ligandPreferences) return;

    setAiCaseStudy({ workflowPhase: 'running' });

    addAiMessage({
      role: 'assistant',
      content: `Starting two-step interface dimer design:\n\n1. Designing asymmetric Chain A (one-sided ${ligandPreferences.ligandLabel} binding)...\n2. Then designing complementary Chain B...\n\nThis may take a few minutes...`,
    });

    const result = await executeLigand(
      {
        ligand_smiles: ligandPreferences.ligandSmiles,
        chain_length: ligandPreferences.chainLength,
        num_designs: ligandPreferences.numDesigns,
        side: 'left',
      },
      {
        onJobCreated: (jobId) => setAiCaseStudy({ pendingJobId: jobId }),
        onProgress: (progress) => setAiCaseStudy({ jobProgress: progress }),
        onStage: (stage) => addAiMessage({ role: 'system', content: stage }),
        onError: (error) => addAiMessage({ role: 'assistant', content: error }),
      }
    );

    if (result) {
      setDesignResults(result.designs);
      setSelectedDesignId(result.designs[0]?.id || null);

      if (result.pdbContent && !result.pdbContent.startsWith('SIMULATED_')) {
        setViewerPdbContent(result.pdbContent);
        setSelectedPdb(result.pdbContent);
        setLatestRfd3Design({
          jobId: aiCaseStudy.pendingJobId || 'demo',
          pdbContent: result.pdbContent,
          source: 'rfd3',
          timestamp: Date.now(),
        });
      }

      setLigandEvaluationResult(result.evaluation);
      setAiCaseStudy({ workflowPhase: 'complete' });

      const passingDesigns = result.designs.filter(
        d => !d.metrics.has_clashes && d.metrics.separable !== false
      ).length;

      const dimer = result.evaluation.dimer;
      if (dimer) {
        addAiMessage({
          role: 'assistant',
          content: `Design complete! Generated ${result.designs.length} designs (${passingDesigns} pass all criteria).\n\n**Best Design (#1):**\n- Combined Affinity: ${dimer.affinity.toFixed(2)} kcal/mol\n- Chain A Contacts: ${dimer.contacts_a}\n- Chain B Contacts: ${dimer.contacts_b}\n- Separable: ${dimer.separable ? 'Yes' : 'No'}\n- Clashes: ${dimer.has_clashes ? 'Detected' : 'None'}\n\nUse the **Design Gallery** below to browse and compare all designs.`,
        });
      }
    }
  };

  // Run binder design
  const handleRunBinderDesign = async () => {
    if (!binderPreferences) return;

    setAiCaseStudy({ workflowPhase: 'running' });

    addAiMessage({
      role: 'assistant',
      content: `Starting protein binder design pipeline:\n\n1. Generating ${binderPreferences.numDesigns} backbone structures...\n2. Running ProteinMPNN sequence design...\n3. Scoring with ESM-3...\n4. Analyzing interfaces...\n5. Filtering by ${binderPreferences.qualityThresholdLabel} threshold...\n\nThis may take a few minutes...`,
    });

    // Get target PDB content
    let targetPdb: string;
    if (binderPreferences.targetPdbContent) {
      targetPdb = binderPreferences.targetPdbContent;
    } else if (binderPreferences.targetPdbId) {
      try {
        addAiMessage({ role: 'system', content: `Fetching target structure ${binderPreferences.targetPdbId}...` });
        const pdbResult = await api.fetchPdb(binderPreferences.targetPdbId);
        if (!pdbResult.content) {
          throw new Error(`Could not fetch PDB content for ${binderPreferences.targetPdbId}`);
        }
        targetPdb = pdbResult.content;
      } catch (err) {
        addAiMessage({
          role: 'assistant',
          content: `Error fetching target: ${err instanceof Error ? err.message : 'Unknown error'}`,
        });
        setAiCaseStudy({ workflowPhase: 'error' });
        return;
      }
    } else {
      addAiMessage({ role: 'assistant', content: 'No target PDB specified.' });
      setAiCaseStudy({ workflowPhase: 'error' });
      return;
    }

    // Build hotspot settings
    const hotspotSettings: Record<string, unknown> = {};
    if (binderPreferences.hotspotMethod === 'auto') {
      hotspotSettings.auto_hotspots = true;
      hotspotSettings.filter_wrap_around = true;
    } else if (binderPreferences.hotspotMethod === 'manual' && binderPreferences.manualHotspots?.length) {
      hotspotSettings.hotspots = binderPreferences.manualHotspots;
      hotspotSettings.auto_hotspots = false;
      hotspotSettings.filter_wrap_around = true;
    } else {
      hotspotSettings.auto_hotspots = false;
      hotspotSettings.filter_wrap_around = false;
    }

    const result = await executeBinder(
      {
        target_pdb: targetPdb,
        binder_length: binderPreferences.binderLength,
        num_designs: binderPreferences.numDesigns,
        quality_threshold: binderPreferences.qualityThreshold,
        protocol: binderPreferences.protocol,
        validate_structure: binderPreferences.validateStructure ?? true,
        ...hotspotSettings,
      },
      {
        onJobCreated: (jobId) => setAiCaseStudy({ pendingJobId: jobId }),
        onProgress: (progress) => setAiCaseStudy({ jobProgress: progress }),
        onStage: (stage) => addAiMessage({ role: 'system', content: stage }),
        onError: (error) => addAiMessage({ role: 'assistant', content: error }),
      }
    );

    if (result) {
      setBinderStatistics(result.statistics);
      setBinderDesigns(result.designs);
      setBinderEvaluationResult(result.evaluation);

      if (result.designs.length > 0) {
        setSelectedBinderDesign(result.designs[0]);
        if (result.designs[0].pdb_content) {
          setViewerPdbContent(result.designs[0].pdb_content);
        }
      }

      setAiCaseStudy({ workflowPhase: 'complete' });

      const best = result.designs[0];
      if (best) {
        addAiMessage({
          role: 'assistant',
          content: `Design complete! Generated ${result.statistics.generated} designs, ${result.statistics.returned} passed filters.\n\n**Best Design (#${best.rank}):**\n- ESM Confidence: ${best.esm_confidence ? `${(best.esm_confidence * 100).toFixed(0)}%` : 'N/A'}\n- ESM Perplexity: ${best.esm_perplexity?.toFixed(1) || 'N/A'}\n- Interface Contacts: ${best.interface_contacts || 'N/A'}\n- H-Bonds: ${best.interface_hbonds || 'N/A'}\n- Packstat: ${best.packstat?.toFixed(2) || 'N/A'}\n\nUse the results panel below to explore all designs.`,
        });
      } else {
        addAiMessage({
          role: 'assistant',
          content: `Design complete! Generated ${result.statistics.generated} designs, but none passed all filters.\n\nTry adjusting parameters or using a different target region.`,
        });
      }
    }
  };

  // Helper function to generate random amino acid sequence
  const generateRandomSequence = (length: number): string => {
    const aminoAcids = 'ACDEFGHIKLMNPQRSTVWY';
    let sequence = '';
    for (let i = 0; i < length; i++) {
      sequence += aminoAcids[Math.floor(Math.random() * aminoAcids.length)];
    }
    return sequence;
  };

  // Helper function to generate simulated individual designs for demo mode
  const generateSimulatedDesigns = (
    configs: Array<{ name: string; tierCounts: Record<string, number> }>,
    designsPerConfig: number
  ): MetalScaffoldDesignResult[] => {
    const designs: MetalScaffoldDesignResult[] = [];
    const tiers: Array<'S' | 'A' | 'B' | 'C' | 'F'> = ['S', 'A', 'B', 'C', 'F'];

    for (const config of configs) {
      // Generate designs based on tier distribution
      let designIndex = 0;
      for (const tier of tiers) {
        const count = config.tierCounts[tier] || 0;
        for (let i = 0; i < count; i++) {
          // Generate realistic metrics based on tier
          const tierMetrics = {
            S: { plddtRange: [0.88, 0.95], ptmRange: [0.82, 0.92], paeRange: [2.0, 4.0] },
            A: { plddtRange: [0.82, 0.88], ptmRange: [0.75, 0.82], paeRange: [4.0, 6.0] },
            B: { plddtRange: [0.75, 0.82], ptmRange: [0.68, 0.75], paeRange: [6.0, 8.0] },
            C: { plddtRange: [0.68, 0.75], ptmRange: [0.60, 0.68], paeRange: [8.0, 10.0] },
            F: { plddtRange: [0.55, 0.68], ptmRange: [0.45, 0.60], paeRange: [10.0, 15.0] },
          };

          const metrics = tierMetrics[tier];
          const plddt = metrics.plddtRange[0] + Math.random() * (metrics.plddtRange[1] - metrics.plddtRange[0]);
          const ptm = metrics.ptmRange[0] + Math.random() * (metrics.ptmRange[1] - metrics.ptmRange[0]);
          const pae = metrics.paeRange[0] + Math.random() * (metrics.paeRange[1] - metrics.paeRange[0]);

          // Generate scaffold length based on config name
          const lengthMap: Record<string, [number, number]> = {
            small: [100, 120],
            medium: [110, 130],
            large: [130, 150],
          };
          const sizeMatch = config.name.match(/^(small|medium|large)/);
          const sizeKey = sizeMatch ? sizeMatch[1] : 'medium';
          const [minLen, maxLen] = lengthMap[sizeKey] || [110, 130];
          const seqLength = Math.floor(minLen + Math.random() * (maxLen - minLen));

          designs.push({
            id: `${config.name}-${designIndex}`,
            name: `${config.name}_design_${designIndex}`,
            configName: config.name,
            sequence: generateRandomSequence(seqLength),
            tier,
            plddt,
            ptm,
            pae,
            status: ['S', 'A', 'B'].includes(tier) ? 'passed' : 'failed',
          });
          designIndex++;
        }
      }
    }

    // Sort by tier (S first) then by pLDDT descending
    const tierOrder: Record<string, number> = { S: 0, A: 1, B: 2, C: 3, F: 4 };
    return designs.sort((a, b) => {
      const tierDiff = tierOrder[a.tier] - tierOrder[b.tier];
      if (tierDiff !== 0) return tierDiff;
      return b.plddt - a.plddt;
    });
  };

  // Export FASTA for metal scaffold designs
  const handleExportMetalDesignsFasta = (designs: MetalScaffoldDesignResult[]) => {
    const passingDesigns = designs.filter(d => d.status === 'passed');
    if (passingDesigns.length === 0) return;

    const fastaContent = passingDesigns
      .map(d => `>${d.name}|tier=${d.tier}|plddt=${d.plddt.toFixed(2)}|ptm=${d.ptm.toFixed(2)}|pae=${d.pae.toFixed(1)}\n${d.sequence}`)
      .join('\n\n');

    const blob = new Blob([fastaContent], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `metal_scaffold_designs_${Date.now()}.fasta`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };

  // Handle metal design selection from gallery
  const handleSelectMetalDesign = (design: MetalScaffoldDesignResult) => {
    setSelectedMetalDesignId(design.id);
    if (design.pdbContent) {
      setViewerPdbContent(design.pdbContent);
    }
  };

  // Run metal scaffold design
  const handleRunMetalScaffoldDesign = async () => {
    if (!metalScaffoldPreferences) return;

    setAiCaseStudy({ workflowPhase: 'running' });

    const modeDescription = metalScaffoldPreferences.optimizationMode === 'sweep'
      ? `Starting parameter sweep with 9 configurations (3 sizes × 3 CFG scales)...\nGenerating ${metalScaffoldPreferences.numDesigns} total designs...`
      : metalScaffoldPreferences.optimizationMode === 'quick'
      ? `Starting quick exploration with single configuration...\nGenerating ${metalScaffoldPreferences.numDesigns} designs...`
      : `Starting production run with optimized configuration...\nGenerating ${metalScaffoldPreferences.numDesigns} designs...`;

    addAiMessage({
      role: 'assistant',
      content: `${modeDescription}\n\nPipeline stages:\n1. RFD3 backbone generation (RASA conditioning for metal burial)\n2. LigandMPNN sequence design (HSAB-biased amino acids)\n3. RF3 structure validation\n4. Quality tier assignment (S/A/B/C/F)\n\nThis may take several minutes...`,
    });

    try {
      // Demo mode - show simulated results
      if (!isBackendConnected) {
        const demoSessionId = `demo-${Date.now()}`;
        setAiCaseStudy({ pendingJobId: demoSessionId });

        await delay(3000);
        addAiMessage({ role: 'system', content: 'Stage 1/4: Generating backbones...' });
        setAiCaseStudy({ jobProgress: 25 });

        await delay(2000);
        addAiMessage({ role: 'system', content: 'Stage 2/4: Running LigandMPNN...' });
        setAiCaseStudy({ jobProgress: 50 });

        await delay(2000);
        addAiMessage({ role: 'system', content: 'Stage 3/4: Validating with RF3...' });
        setAiCaseStudy({ jobProgress: 75 });

        await delay(2000);
        addAiMessage({ role: 'system', content: 'Stage 4/4: Filtering and ranking...' });
        setAiCaseStudy({ jobProgress: 100 });

        // Simulated sweep results
        const simulatedResults = {
          configs: [
            { name: 'medium_cfg2.0', passRate: 0.67, avgPlddt: 0.84, avgPtm: 0.78, tierCounts: { S: 2, A: 4, B: 2, C: 1, F: 1 } },
            { name: 'large_cfg2.0', passRate: 0.56, avgPlddt: 0.82, avgPtm: 0.76, tierCounts: { S: 1, A: 3, B: 3, C: 2, F: 1 } },
            { name: 'small_cfg2.5', passRate: 0.44, avgPlddt: 0.81, avgPtm: 0.75, tierCounts: { S: 1, A: 2, B: 3, C: 3, F: 1 } },
            { name: 'medium_cfg1.5', passRate: 0.44, avgPlddt: 0.80, avgPtm: 0.74, tierCounts: { S: 0, A: 3, B: 4, C: 2, F: 1 } },
            { name: 'small_cfg2.0', passRate: 0.33, avgPlddt: 0.79, avgPtm: 0.73, tierCounts: { S: 0, A: 2, B: 4, C: 3, F: 1 } },
          ],
          bestConfig: 'medium_cfg2.0',
          totalDesigns: metalScaffoldPreferences.numDesigns,
          passingDesigns: Math.floor(metalScaffoldPreferences.numDesigns * 0.45),
        };

        // Generate individual designs for the gallery
        const individualDesigns = generateSimulatedDesigns(
          simulatedResults.configs,
          metalScaffoldPreferences.designsPerConfig
        );

        setMetalScaffoldSweepResults(simulatedResults);
        setMetalScaffoldDesigns(individualDesigns);
        setAiCaseStudy({ workflowPhase: 'complete' });

        const best = simulatedResults.configs[0];
        const totalPassing = individualDesigns.filter(d => d.status === 'passed').length;
        const sTierCount = individualDesigns.filter(d => d.tier === 'S').length;
        const aTierCount = individualDesigns.filter(d => d.tier === 'A').length;

        addAiMessage({
          role: 'assistant',
          content: `Design complete! Parameter sweep finished with ${individualDesigns.length} total designs.\n\n**Best Configuration: ${best.name}**\n- Pass Rate: ${(best.passRate * 100).toFixed(0)}%\n- Avg pLDDT: ${best.avgPlddt.toFixed(2)}\n- Avg pTM: ${best.avgPtm.toFixed(2)}\n- Tier Distribution: ${best.tierCounts.S}× S, ${best.tierCounts.A}× A, ${best.tierCounts.B}× B\n\n**Summary:**\n- ${totalPassing} designs passed (tier B or better)\n- ${sTierCount + aTierCount} high-quality (S/A tier) designs\n\nBrowse the **Individual Design Results** below to view, compare, and export passing designs.`,
        });
        return;
      }

      // Real backend mode - fetch reference structure and run sweep
      // Generate session ID upfront for progress tracking
      const sessionId = `sweep-${Date.now().toString(36)}-${Math.random().toString(36).slice(2, 8)}`;
      setAiCaseStudy({ pendingJobId: sessionId });

      addAiMessage({
        role: 'assistant',
        content: `Backend connected! Fetching reference Tb-citrate structure (3C9H)...`,
      });

      // Reference PDB IDs for common metal-ligand combinations
      const REFERENCE_PDBS: Record<string, string> = {
        'TB-CIT': '3C9H',  // Terbium-citrate reference
        'EU-CIT': '3C9H',  // Europium uses similar coordination
        'GD-CIT': '3C9H',  // Gadolinium uses similar coordination
        'CA-CIT': '1L3J',  // Calcium-citrate
        'ZN-CIT': '1CA2',  // Zinc - carbonic anhydrase
      };

      // Get reference PDB based on metal-ligand selection
      const metalLigandKey = `${metalScaffoldPreferences.metal}-${metalScaffoldPreferences.ligand}`;
      const referencePdbId = REFERENCE_PDBS[metalLigandKey] || '3C9H';

      // Progress simulation for long-running sweep
      // Each design takes ~90s (RFD3: 50s, MPNN: 10s, RF3: 30s)
      const totalDesigns = metalScaffoldPreferences.numDesigns;
      const estimatedTimePerDesign = 90; // seconds
      const totalEstimatedTime = totalDesigns * estimatedTimePerDesign;
      let progressInterval: NodeJS.Timeout | null = null;
      let currentProgress = 5;
      let currentDesign = 0;
      const startTime = Date.now();

      // Start progress simulation
      const startProgressSimulation = () => {
        // Set initial stage
        setAiCaseStudy({ currentStage: 'Initializing parameter sweep...' });

        progressInterval = setInterval(() => {
          const elapsed = (Date.now() - startTime) / 1000;
          const estimatedDesignsComplete = Math.floor(elapsed / estimatedTimePerDesign);

          if (estimatedDesignsComplete > currentDesign) {
            currentDesign = estimatedDesignsComplete;
            const configNum = Math.floor(currentDesign / metalScaffoldPreferences.designsPerConfig) + 1;
            const totalConfigs = metalScaffoldPreferences.optimizationMode === 'sweep' ? 9 : 1;
            const designInConfig = (currentDesign % metalScaffoldPreferences.designsPerConfig) + 1;

            // Update progress message
            if (configNum <= totalConfigs) {
              const stageMessage = `Config ${configNum}/${totalConfigs}: Design ${designInConfig}/${metalScaffoldPreferences.designsPerConfig}`;
              setAiCaseStudy({ currentStage: stageMessage });
              addAiMessage({
                role: 'system',
                content: `${stageMessage} complete...`
              });
            }
          }

          // Calculate progress (10% for setup, 85% for designs, 5% for filtering)
          currentProgress = Math.min(95, 10 + (elapsed / totalEstimatedTime) * 85);
          setAiCaseStudy({ jobProgress: Math.round(currentProgress) });
        }, 5000); // Update every 5 seconds
      };

      try {
        // Fetch the reference structure from RCSB
        setAiCaseStudy({ currentStage: `Fetching ${referencePdbId} from RCSB...` });
        addAiMessage({ role: 'system', content: `Fetching reference structure ${referencePdbId} from RCSB...` });
        const pdbResult = await api.fetchPdb(referencePdbId);
        const motifPdbContent = pdbResult.content;

        if (!motifPdbContent) {
          throw new Error(`Failed to fetch reference structure ${referencePdbId}`);
        }

        setAiCaseStudy({ currentStage: 'Starting parameter sweep...', jobProgress: 10 });
        addAiMessage({ role: 'system', content: `Reference structure loaded. Starting parameter sweep...` });

        // Start progress simulation
        startProgressSimulation();

        // Call the actual backend API
        const response = await api.startMetalBindingSweep({
          motif_pdb: motifPdbContent,
          metal: metalScaffoldPreferences.metal,
          ligand: metalScaffoldPreferences.ligand !== 'none' ? metalScaffoldPreferences.ligand : undefined,
          designs_per_config: metalScaffoldPreferences.designsPerConfig,
          filters: { plddt: 0.70, ptm: 0.60, pae: 10.0 },
        });

        // Stop progress simulation
        if (progressInterval) {
          clearInterval(progressInterval);
          progressInterval = null;
        }

        if (response.status === 'completed' && response.result) {
          const result = response.result;

          // Update with real session ID from backend
          if (result.session_id) {
            setAiCaseStudy({ pendingJobId: result.session_id });
          }

          setAiCaseStudy({ currentStage: 'Filtering and ranking designs...', jobProgress: 95 });
          addAiMessage({ role: 'system', content: `Sweep complete! Filtering and ranking ${result.total_generated} designs...` });

          // Convert backend result to frontend format
          // API returns: config_rankings with tier_distribution
          const backendResults = {
            configs: result.config_rankings?.map((cfg) => ({
              name: cfg.config_name,
              passRate: cfg.pass_rate,
              avgPlddt: cfg.avg_plddt,
              avgPtm: cfg.avg_ptm,
              tierCounts: {
                S: cfg.tier_distribution?.['S'] || 0,
                A: cfg.tier_distribution?.['A'] || 0,
                B: cfg.tier_distribution?.['B'] || 0,
                C: cfg.tier_distribution?.['C'] || 0,
                F: cfg.tier_distribution?.['F'] || 0,
              },
            })) || [],
            bestConfig: result.config_rankings?.[0]?.config_name || 'medium_cfg2.0',
            totalDesigns: result.total_generated || metalScaffoldPreferences.numDesigns,
            passingDesigns: result.total_passing || 0,
          };

          setMetalScaffoldSweepResults(backendResults);
          setAiCaseStudy({ workflowPhase: 'complete', jobProgress: 100, currentStage: null });

          const best = backendResults.configs[0] || { name: 'unknown', passRate: 0, avgPlddt: 0, avgPtm: 0, tierCounts: { S: 0, A: 0, B: 0, C: 0, F: 0 } };
          addAiMessage({
            role: 'assistant',
            content: `Design complete! Parameter sweep finished with ${backendResults.totalDesigns} total designs.\n\n**Best Configuration: ${best.name}**\n- Pass Rate: ${(best.passRate * 100).toFixed(0)}%\n- Avg pLDDT: ${best.avgPlddt.toFixed(2)}\n- Avg pTM: ${best.avgPtm.toFixed(2)}\n- Tier Distribution: ${best.tierCounts.S}× S, ${best.tierCounts.A}× A, ${best.tierCounts.B}× B\n\n**Summary:**\n- ${backendResults.passingDesigns} designs passed (tier B or better)`,
          });
        } else {
          throw new Error(`Backend returned status: ${response.status}`);
        }
      } catch (backendError) {
        // Stop progress simulation on error
        if (progressInterval) {
          clearInterval(progressInterval);
          progressInterval = null;
        }

        console.error('[AIPanel] Backend sweep error:', backendError);
        addAiMessage({
          role: 'assistant',
          content: `Backend error: ${backendError instanceof Error ? backendError.message : 'Unknown error'}\n\nFalling back to demo mode...`,
        });

        // Fallback to demo results
        await delay(2000);
        const simulatedResults = {
          configs: [
            { name: 'medium_cfg2.0', passRate: 0.67, avgPlddt: 0.84, avgPtm: 0.78, tierCounts: { S: 2, A: 4, B: 2, C: 1, F: 1 } },
            { name: 'large_cfg2.0', passRate: 0.56, avgPlddt: 0.82, avgPtm: 0.76, tierCounts: { S: 1, A: 3, B: 3, C: 2, F: 1 } },
            { name: 'small_cfg2.5', passRate: 0.44, avgPlddt: 0.81, avgPtm: 0.75, tierCounts: { S: 1, A: 2, B: 3, C: 3, F: 1 } },
          ],
          bestConfig: 'medium_cfg2.0',
          totalDesigns: metalScaffoldPreferences.numDesigns,
          passingDesigns: Math.floor(metalScaffoldPreferences.numDesigns * 0.45),
        };

        setMetalScaffoldSweepResults(simulatedResults);
        setAiCaseStudy({ workflowPhase: 'complete', currentStage: null });

        const best = simulatedResults.configs[0];
        addAiMessage({
          role: 'assistant',
          content: `(Demo mode) Design complete! Parameter sweep finished with ${simulatedResults.totalDesigns} total designs.\n\n**Best Configuration: ${best.name}**\n- Pass Rate: ${(best.passRate * 100).toFixed(0)}%\n- Avg pLDDT: ${best.avgPlddt.toFixed(2)}\n- Avg pTM: ${best.avgPtm.toFixed(2)}\n- Tier Distribution: ${best.tierCounts.S}× S, ${best.tierCounts.A}× A, ${best.tierCounts.B}× B`,
        });
      }
    } catch (error) {
      console.error('[AIPanel] Metal scaffold design error:', error);
      addAiMessage({
        role: 'assistant',
        content: `Error running metal scaffold design: ${error instanceof Error ? error.message : 'Unknown error'}\n\nPlease try again or check the backend connection.`,
      });
      setAiCaseStudy({ workflowPhase: 'error', pendingJobId: undefined, currentStage: null });
    }
  };

  // Handle design selection from gallery
  const handleSelectDesign = (design: DesignResult) => {
    setSelectedDesignId(design.id);

    if (design.pdbContent && !design.pdbContent.startsWith('SIMULATED_')) {
      setViewerPdbContent(design.pdbContent);
      setSelectedPdb(design.pdbContent);
    }

    setLigandEvaluationResult(createEvaluationFromDesign(design));

    addAiMessage({
      role: 'assistant',
      content: `Showing Design #${design.rank}:\n\n- Affinity: ${design.metrics.affinity?.toFixed(2) ?? 'N/A'} kcal/mol\n- Contacts: ${design.metrics.contacts_a ?? 0} (A) / ${design.metrics.contacts_b ?? 0} (B)\n- Status: ${!design.metrics.has_clashes && design.metrics.separable !== false ? 'PASS' : 'FAIL'}`,
    });
  };

  // Handle edit/retry/start new
  const handleEditPreferences = () => setAiCaseStudy({ workflowPhase: 'interview' });

  const handleRetry = () => {
    setEvaluationResult(null);
    setLigandEvaluationResult(null);
    setDesignResults([]);
    setSelectedDesignId(null);
    setViewerPdbContent(null);
    setAiCaseStudy({ workflowPhase: 'interview', evaluationResult: null });
    addAiMessage({ role: 'assistant', content: `Let's try again with different settings. What would you like to change?` });
  };

  const handleStartNew = () => {
    clearAiConversation();
    setUserPreferences(null);
    setAnalysisResult(null);
    setEvaluationResult(null);
    setPdbInput('');
    setIsLigandDesign(false);
    setLigandAnalysisResult(null);
    setLigandPreferences(null);
    setLigandEvaluationResult(null);
    setIsBinderDesign(false);
    setBinderPreferences(null);
    setBinderEvaluationResult(null);
    setBinderDesigns([]);
    setBinderStatistics(null);
    setSelectedBinderDesign(null);
    setIsMetalScaffoldDesign(false);
    setMetalScaffoldPreferences(null);
    setMetalScaffoldSweepResults(null);
    setMetalScaffoldDesigns([]);
    setSelectedMetalDesignId(null);
    setDesignResults([]);
    setSelectedDesignId(null);
    setViewerPdbContent(null);
  };

  // Quick start demos
  const handleQuickStart = async () => {
    setPdbInput('1BRF');
    await delay(100);
    handleStructureInput();
  };

  return (
    <div className="h-full flex flex-col">
      {/* Header */}
      <div className="border-b border-border pb-6 mb-6">
        <div className="flex items-center justify-between mb-2">
          <div className="flex items-center gap-3">
            <span className="bg-gradient-to-r from-primary to-primary/80 text-primary-foreground text-[11px] font-bold px-3 py-1 rounded-full uppercase tracking-wide flex items-center gap-1.5">
              <Sparkles className="h-4 w-4" />
              Phase 0
            </span>
            <h2 className="text-xl font-bold text-foreground">AI Design Assistant</h2>
          </div>
          {workflowPhase !== 'idle' && (
            <button
              onClick={handleStartNew}
              className="text-sm text-muted-foreground hover:text-foreground flex items-center gap-1"
            >
              <RefreshCw className="h-4 w-4" />
              Start New
            </button>
          )}
        </div>
        <p className="text-muted-foreground text-sm leading-relaxed max-w-3xl pl-1">
          {workflowPhase === 'idle'
            ? "Describe what you want to design in plain English, like 'Design a protein to bind citrate with terbium'. I'll handle the rest."
            : "Answer simple questions about your design goals. I'll handle the technical parameters for you."}
        </p>

        {!isBackendConnected && (
          <div className="mt-3 flex items-center gap-2 text-xs text-amber-600 bg-amber-50 px-3 py-1.5 rounded-lg w-fit">
            <Info className="h-4 w-4" />
            Demo mode - backend not connected
          </div>
        )}
      </div>

      {/* Quick Start Section */}
      {workflowPhase === 'idle' && aiCaseStudy.conversation.length === 0 && (
        <div className="mb-6 grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4">
          <QuickStartCard
            icon="token"
            title="Metal Binding Redesign"
            description="Convert iron-binding Rubredoxin (1BRF) to bind terbium."
            buttonText="Try 1BRF"
            onClick={handleQuickStart}
          />
          <QuickStartCard
            icon="token"
            title="Metal Scaffold Design"
            description="Design monomer scaffolds with parameter sweep optimization."
            buttonText="Try METAL"
            onClick={() => { setPdbInput('METAL'); delay(100).then(handleStructureInput); }}
          />
          <QuickStartCard
            icon="science"
            title="Interface Dimer Design"
            description="Design separable dimer with azobenzene at the interface."
            buttonText="Try AZOB"
            onClick={() => { setPdbInput('AZOB'); delay(100).then(handleStructureInput); }}
          />
          <QuickStartCard
            icon="link"
            title="Protein Binder Design"
            description="Design high-affinity binders with ESM-3 validation."
            buttonText="Try BIND"
            onClick={() => { setPdbInput('BIND'); delay(100).then(handleStructureInput); }}
          />
        </div>
      )}

      {/* Conversation Area */}
      <div className="flex-1 overflow-y-auto space-y-4 pb-4 min-h-[300px]">
        {aiCaseStudy.conversation.map((message, index) => {
          if (message.role === 'user') return <UserMessage key={index} content={message.content} />;
          if (message.role === 'system') return <SystemMessage key={index} content={message.content} />;
          return (
            <AIMessage
              key={index}
              content={message.content}
              isNew={index === aiCaseStudy.conversation.length - 1 && aiCaseStudy.isProcessing}
            />
          );
        })}

        {/* Ligand Analysis Card */}
        {isLigandDesign && ligandAnalysisResult && workflowPhase === 'interview' && (
          <div className="pl-11 mb-4">
            <LigandAnalysisCard analysis={ligandAnalysisResult} />
          </div>
        )}

        {/* Interview Modes */}
        {workflowPhase === 'interview' && !isLigandDesign && !isBinderDesign && !isMetalScaffoldDesign && (
          <div className="pl-11">
            <InterviewMode onComplete={handleInterviewComplete} onCancel={handleStartNew} />
          </div>
        )}

        {workflowPhase === 'interview' && isLigandDesign && (
          <div className="pl-11">
            <LigandInterviewMode onComplete={handleLigandInterviewComplete} onCancel={handleStartNew} />
          </div>
        )}

        {workflowPhase === 'interview' && isBinderDesign && (
          <div className="pl-11">
            <BinderInterviewMode onComplete={handleBinderInterviewComplete} onCancel={handleStartNew} />
          </div>
        )}

        {workflowPhase === 'interview' && isMetalScaffoldDesign && (
          <div className="pl-11">
            <MetalScaffoldInterviewMode onComplete={handleMetalScaffoldInterviewComplete} onCancel={handleStartNew} />
          </div>
        )}

        {/* Hotspot Selector */}
        {showHotspotSelector && hotspotPdbContent && (
          <div className="pl-11">
            <HotspotSelector
              pdbContent={hotspotPdbContent}
              targetChain="A"
              onConfirm={handleHotspotConfirm}
              onCancel={handleHotspotCancel}
              maxHotspots={3}
            />
          </div>
        )}

        {/* Preference Summaries */}
        {workflowPhase === 'confirming' && userPreferences && !isLigandDesign && !isBinderDesign && !isMetalScaffoldDesign && (
          <div className="pl-11">
            <PreferenceSummaryCard
              preferences={userPreferences}
              onConfirm={handleRunDesign}
              onEdit={handleEditPreferences}
              isRunning={isRunningJob}
            />
          </div>
        )}

        {workflowPhase === 'confirming' && ligandPreferences && isLigandDesign && (
          <div className="pl-11">
            <LigandPreferenceSummaryCard
              preferences={ligandPreferences}
              onConfirm={handleRunLigandDesign}
              onEdit={handleEditPreferences}
              isRunning={isRunningJob}
            />
          </div>
        )}

        {workflowPhase === 'confirming' && binderPreferences && isBinderDesign && (
          <div className="pl-11">
            <BinderPreferenceSummaryCard
              preferences={binderPreferences}
              onConfirm={handleRunBinderDesign}
              onEdit={handleEditPreferences}
              isRunning={isRunningJob}
            />
          </div>
        )}

        {workflowPhase === 'confirming' && metalScaffoldPreferences && isMetalScaffoldDesign && (
          <div className="pl-11">
            <MetalScaffoldPreferenceSummaryCard
              preferences={metalScaffoldPreferences}
              onConfirm={handleRunMetalScaffoldDesign}
              onEdit={handleEditPreferences}
              isRunning={isRunningJob}
            />
          </div>
        )}

        {/* Job Progress */}
        {workflowPhase === 'running' && (
          <div className="pl-11">
            <JobProgressCard
              jobId={aiCaseStudy.pendingJobId || 'sweep-pending'}
              status="running"
              progress={aiCaseStudy.jobProgress || workflowState.progress}
              message={aiCaseStudy.currentStage || workflowState.currentStage || 'Processing...'}
            />
          </div>
        )}

        {/* Results - Metal Binding */}
        {workflowPhase === 'complete' && evaluationResult && userPreferences && !isLigandDesign && !isBinderDesign && !isMetalScaffoldDesign && (
          <div className="pl-11 space-y-4">
            <EvaluationCard
              evaluation={evaluationResult}
              targetMetal={userPreferences.targetMetal}
              onRetry={handleRetry}
            />
            {viewerPdbContent && (
              <StructureViewerSection
                pdbContent={viewerPdbContent}
                title="Designed Structure"
                badge={`${userPreferences.targetMetalLabel} binding site`}
              />
            )}
          </div>
        )}

        {/* Results - Ligand Design */}
        {workflowPhase === 'complete' && ligandEvaluationResult && isLigandDesign && (
          <div className="pl-11 space-y-4">
            <StructureViewerSection
              pdbContent={viewerPdbContent}
              title="Structure Viewer"
              badge={selectedDesignId ? `Design #${designResults.find(d => d.id === selectedDesignId)?.rank || '?'}` : undefined}
              colorLegend={
                <div className="flex items-center gap-1 text-xs text-muted-foreground">
                  <span className="w-3 h-3 rounded-full bg-purple-400" />
                  <span>Chain A</span>
                  <span className="ml-2 w-3 h-3 rounded-full bg-cyan-400" />
                  <span>Chain B</span>
                </div>
              }
              emptyMessage="Demo Mode - No Structure Available"
            />

            {designResults.length > 0 && (
              <DesignGallery
                designs={designResults}
                selectedDesignId={selectedDesignId}
                onSelectDesign={handleSelectDesign}
                designType="ligand_interface"
              />
            )}

            <LigandEvaluationCard evaluation={ligandEvaluationResult} expanded={true} />

            <button
              onClick={handleRetry}
              className="w-full px-4 py-2.5 text-primary bg-card border border-border rounded-xl font-medium text-sm hover:bg-muted transition-all flex items-center justify-center gap-2"
            >
              <RefreshCw className="h-4 w-4" />
              Try Different Settings
            </button>
          </div>
        )}

        {/* Results - Binder Design */}
        {workflowPhase === 'complete' && binderEvaluationResult && isBinderDesign && (
          <div className="pl-11 space-y-4">
            {binderStatistics && <PipelineFunnel statistics={binderStatistics} />}

            <BinderEvaluationCard evaluation={binderEvaluationResult} expanded={true} />

            {/* Design list */}
            {binderDesigns.length > 1 && (
              <div className="bg-card rounded-xl border border-border overflow-hidden">
                <div className="bg-muted px-4 py-3 border-b border-border">
                  <div className="flex items-center gap-2">
                    <ListOrdered className="h-5 w-5 text-primary" />
                    <h4 className="font-semibold text-foreground text-sm">All Designs ({binderDesigns.length})</h4>
                  </div>
                </div>
                <div className="p-4 space-y-2">
                  {binderDesigns.map((design) => (
                    <button
                      key={design.rank}
                      onClick={() => setSelectedBinderDesign(design)}
                      className={`w-full text-left p-3 rounded-lg border transition-all ${
                        selectedBinderDesign?.rank === design.rank
                          ? 'border-primary bg-primary/5'
                          : 'border-border hover:border-primary/50 hover:bg-muted'
                      }`}
                    >
                      <div className="flex items-center justify-between">
                        <div className="flex items-center gap-3">
                          <span className={`w-6 h-6 rounded-full flex items-center justify-center text-xs font-bold ${
                            design.rank === 1 ? 'bg-amber-100 text-amber-700' : 'bg-muted text-muted-foreground'
                          }`}>
                            {design.rank}
                          </span>
                          <div>
                            <div className="text-sm font-medium text-foreground">Design #{design.rank}</div>
                            <div className="text-xs text-muted-foreground">
                              ESM: {design.esm_confidence ? `${(design.esm_confidence * 100).toFixed(0)}%` : 'N/A'} |
                              Contacts: {design.interface_contacts || 'N/A'} |
                              H-bonds: {design.interface_hbonds || 'N/A'}
                            </div>
                          </div>
                        </div>
                        {selectedBinderDesign?.rank === design.rank && (
                          <CheckCircle className="h-5 w-5 text-primary" />
                        )}
                      </div>
                    </button>
                  ))}
                </div>
              </div>
            )}

            {/* Structure viewer for selected binder */}
            {selectedBinderDesign?.pdb_content && (
              <StructureViewerSection
                pdbContent={selectedBinderDesign.pdb_content}
                title="Structure Viewer"
                badge={`Design #${selectedBinderDesign.rank}`}
                colorLegend={
                  <div className="flex items-center gap-1 text-xs text-muted-foreground">
                    <span className="w-3 h-3 rounded-full bg-muted-foreground/60" />
                    <span>Target</span>
                    <span className="ml-2 w-3 h-3 rounded-full bg-violet-400" />
                    <span>Binder</span>
                  </div>
                }
              />
            )}

            <button
              onClick={handleRetry}
              className="w-full px-4 py-2.5 text-primary bg-card border border-border rounded-xl font-medium text-sm hover:bg-muted transition-all flex items-center justify-center gap-2"
            >
              <RefreshCw className="h-4 w-4" />
              Try Different Settings
            </button>
          </div>
        )}

        {/* Results - Metal Scaffold Design */}
        {workflowPhase === 'complete' && metalScaffoldSweepResults && isMetalScaffoldDesign && (
          <div className="pl-11 space-y-4">
            {/* Sweep Results Summary */}
            <div className="bg-card rounded-xl border border-border overflow-hidden">
              <div className="bg-muted px-4 py-3 border-b border-border">
                <div className="flex items-center justify-between">
                  <div className="flex items-center gap-2">
                    <Gem className="h-5 w-5 text-primary" />
                    <h4 className="font-semibold text-foreground text-sm">Parameter Sweep Results</h4>
                  </div>
                  <div className="flex items-center gap-2 text-xs">
                    <span className="text-muted-foreground">
                      {metalScaffoldSweepResults.passingDesigns}/{metalScaffoldSweepResults.totalDesigns} passed
                    </span>
                    <span className="px-2 py-0.5 rounded-full bg-green-100 text-green-700 font-medium">
                      {((metalScaffoldSweepResults.passingDesigns / metalScaffoldSweepResults.totalDesigns) * 100).toFixed(0)}%
                    </span>
                  </div>
                </div>
              </div>

              <div className="p-4 space-y-3">
                {metalScaffoldSweepResults.configs.map((config, idx) => (
                  <div
                    key={config.name}
                    className={`p-3 rounded-lg border transition-all ${
                      config.name === metalScaffoldSweepResults.bestConfig
                        ? 'border-primary bg-primary/5'
                        : 'border-border'
                    }`}
                  >
                    <div className="flex items-center justify-between mb-2">
                      <div className="flex items-center gap-2">
                        <span className={`w-6 h-6 rounded-full flex items-center justify-center text-xs font-bold ${
                          idx === 0 ? 'bg-amber-100 text-amber-700' : 'bg-muted text-muted-foreground'
                        }`}>
                          {idx + 1}
                        </span>
                        <span className="font-medium text-foreground">{config.name}</span>
                        {config.name === metalScaffoldSweepResults.bestConfig && (
                          <span className="text-xs px-2 py-0.5 rounded-full bg-primary text-primary-foreground">
                            Best
                          </span>
                        )}
                      </div>
                      <span className="text-sm font-medium text-foreground">
                        {(config.passRate * 100).toFixed(0)}% pass
                      </span>
                    </div>

                    <div className="grid grid-cols-3 gap-4 text-xs">
                      <div>
                        <span className="text-muted-foreground">Avg pLDDT:</span>
                        <span className="ml-1 font-medium text-foreground">{config.avgPlddt.toFixed(2)}</span>
                      </div>
                      <div>
                        <span className="text-muted-foreground">Avg pTM:</span>
                        <span className="ml-1 font-medium text-foreground">{config.avgPtm.toFixed(2)}</span>
                      </div>
                      <div className="flex items-center gap-1">
                        <span className="w-2 h-2 rounded-full bg-green-500" title="S-tier" />
                        <span>{config.tierCounts.S || 0}</span>
                        <span className="w-2 h-2 rounded-full bg-blue-500 ml-1" title="A-tier" />
                        <span>{config.tierCounts.A || 0}</span>
                        <span className="w-2 h-2 rounded-full bg-yellow-500 ml-1" title="B-tier" />
                        <span>{config.tierCounts.B || 0}</span>
                      </div>
                    </div>
                  </div>
                ))}
              </div>
            </div>

            {/* Structure Viewer for selected design */}
            <StructureViewerSection
              pdbContent={viewerPdbContent}
              title="Structure Viewer"
              badge={
                selectedMetalDesignId
                  ? `${metalScaffoldDesigns.find(d => d.id === selectedMetalDesignId)?.name || 'Design'}`
                  : undefined
              }
              colorLegend={
                <div className="flex items-center gap-1 text-xs text-muted-foreground">
                  <span className="w-3 h-3 rounded-full bg-primary/60" />
                  <span>Scaffold</span>
                  <span className="ml-2 w-3 h-3 rounded-full bg-amber-400" />
                  <span>Metal Site</span>
                </div>
              }
              emptyMessage="Select a design to view structure (demo mode - no PDB available)"
            />

            {/* Individual Design Results Gallery */}
            {metalScaffoldDesigns.length > 0 && (
              <MetalScaffoldDesignGallery
                designs={metalScaffoldDesigns}
                selectedDesignId={selectedMetalDesignId}
                onSelectDesign={handleSelectMetalDesign}
                onExportFasta={handleExportMetalDesignsFasta}
              />
            )}

            {/* Quality Tier Legend */}
            <div className="bg-muted rounded-xl p-4 border border-border">
              <h5 className="font-medium text-foreground text-sm mb-2">Quality Tier Legend</h5>
              <div className="grid grid-cols-5 gap-2 text-xs">
                <div className="flex items-center gap-1">
                  <span className="w-3 h-3 rounded-full bg-green-500" />
                  <span className="text-foreground font-medium">S</span>
                  <span className="text-muted-foreground">Excellent</span>
                </div>
                <div className="flex items-center gap-1">
                  <span className="w-3 h-3 rounded-full bg-blue-500" />
                  <span className="text-foreground font-medium">A</span>
                  <span className="text-muted-foreground">Good</span>
                </div>
                <div className="flex items-center gap-1">
                  <span className="w-3 h-3 rounded-full bg-yellow-500" />
                  <span className="text-foreground font-medium">B</span>
                  <span className="text-muted-foreground">Pass</span>
                </div>
                <div className="flex items-center gap-1">
                  <span className="w-3 h-3 rounded-full bg-orange-500" />
                  <span className="text-foreground font-medium">C</span>
                  <span className="text-muted-foreground">Marginal</span>
                </div>
                <div className="flex items-center gap-1">
                  <span className="w-3 h-3 rounded-full bg-red-500" />
                  <span className="text-foreground font-medium">F</span>
                  <span className="text-muted-foreground">Fail</span>
                </div>
              </div>
            </div>

            <button
              onClick={handleRetry}
              className="w-full px-4 py-2.5 text-primary bg-card border border-border rounded-xl font-medium text-sm hover:bg-muted transition-all flex items-center justify-center gap-2"
            >
              <RefreshCw className="h-4 w-4" />
              Try Different Settings
            </button>
          </div>
        )}

        {/* Loading indicator */}
        {aiCaseStudy.isProcessing && workflowPhase !== 'interview' && (
          <div className="flex gap-3">
            <div className="flex-shrink-0 w-8 h-8 rounded-full bg-gradient-to-br from-primary to-primary/80 flex items-center justify-center">
              <Loader2 className="h-4 w-4 text-primary-foreground animate-spin" />
            </div>
            <div className="flex items-center gap-2 text-sm text-muted-foreground">
              <span>
                {workflowPhase === 'fetching' && 'Fetching structure from RCSB...'}
                {workflowPhase === 'analyzing' && !isLigandDesign && !isBinderDesign && !isMetalScaffoldDesign && 'Analyzing metal binding site...'}
                {workflowPhase === 'analyzing' && isLigandDesign && 'Analyzing ligand topology...'}
                {workflowPhase === 'analyzing' && isBinderDesign && 'Preparing binder design pipeline...'}
                {workflowPhase === 'analyzing' && isMetalScaffoldDesign && 'Preparing metal scaffold design pipeline...'}
                {workflowPhase === 'evaluating' && !isLigandDesign && !isBinderDesign && !isMetalScaffoldDesign && 'Evaluating design results...'}
                {workflowPhase === 'evaluating' && isLigandDesign && 'Evaluating dimer metrics...'}
                {workflowPhase === 'evaluating' && isBinderDesign && 'Evaluating binder designs...'}
                {workflowPhase === 'evaluating' && isMetalScaffoldDesign && 'Evaluating metal scaffold designs...'}
              </span>
            </div>
          </div>
        )}

        <div ref={messagesEndRef} />
      </div>

      {/* NL Design Pipeline Progress */}
      {isNLDesign && nlIsRunning && (
        <div className="mb-4">
          <AIDesignPipelineWorkflow
            currentStage={nlStage}
            stageInfo={nlStageInfo}
            error={nlError}
            query={pdbInput}
          />
        </div>
      )}

      {/* Input Area */}
      <div className="border-t border-border pt-4 mt-auto">
        {(workflowPhase === 'idle' || workflowPhase === 'structure_input') ? (
          <div className="flex gap-3">
            <div className="flex-1 relative">
              <textarea
                value={pdbInput || ''}
                onChange={(e) => setPdbInput(e.target.value)}
                onKeyDown={(e) => {
                  if (e.key === 'Enter' && !e.shiftKey) {
                    e.preventDefault();
                    handleStructureInput();
                  }
                }}
                placeholder="Describe your protein design (e.g., 'Design a protein to bind citrate with terbium') or enter a PDB code"
                className="w-full px-4 py-3 bg-muted rounded-xl border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 focus:outline-none text-foreground text-sm transition-all resize-none min-h-[48px] max-h-[120px]"
                disabled={aiCaseStudy.isProcessing || nlIsRunning}
                rows={1}
                style={{ height: 'auto', minHeight: '48px' }}
                onInput={(e) => {
                  const target = e.target as HTMLTextAreaElement;
                  target.style.height = 'auto';
                  target.style.height = Math.min(target.scrollHeight, 120) + 'px';
                }}
              />
            </div>
            <button
              onClick={handleStructureInput}
              disabled={!pdbInput.trim() || aiCaseStudy.isProcessing || nlIsRunning}
              className="px-6 py-3 bg-primary hover:bg-primary/90 disabled:bg-muted disabled:cursor-not-allowed text-primary-foreground rounded-xl font-medium text-sm transition-all flex items-center gap-2 self-end"
            >
              {nlIsRunning ? (
                <Loader2 className="h-5 w-5 animate-spin" />
              ) : (
                <Send className="h-5 w-5" />
              )}
              {nlIsRunning ? 'Running...' : 'Send'}
            </button>
          </div>
        ) : (
          <div className="flex gap-3">
            <input
              type="text"
              value={followUpInput}
              onChange={(e) => setFollowUpInput(e.target.value)}
              placeholder={workflowPhase === 'complete' ? "Ask a follow-up question..." : "Type a message..."}
              className="flex-1 px-4 py-3 bg-muted rounded-xl border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 focus:outline-none text-foreground text-sm transition-all"
              disabled={aiCaseStudy.isProcessing || workflowPhase === 'interview' || workflowPhase === 'running'}
            />
            <button
              disabled={aiCaseStudy.isProcessing || workflowPhase === 'interview' || workflowPhase === 'running' || !followUpInput.trim()}
              className="px-4 py-3 bg-primary hover:bg-primary/90 disabled:bg-muted disabled:cursor-not-allowed text-primary-foreground rounded-xl font-medium text-sm transition-all flex items-center gap-2"
            >
              <Send className="h-5 w-5" />
            </button>
          </div>
        )}
        <p className="mt-2 text-xs text-muted-foreground text-center">
          {workflowPhase === 'idle' && "Describe your design in natural language, or enter a PDB code (1BRF, METAL, AZOB, BIND)"}
          {workflowPhase === 'interview' && "Answer the questions above to configure your design"}
          {workflowPhase === 'confirming' && "Review your preferences and run the design"}
          {workflowPhase === 'running' && "Design in progress..."}
          {workflowPhase === 'complete' && "Design complete! You can start a new design or continue exploring."}
        </p>
      </div>
    </div>
  );
}

// Quick start card component
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
