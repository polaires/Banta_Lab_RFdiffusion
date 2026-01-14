'use client';

import { useState, useRef, useEffect } from 'react';
import dynamic from 'next/dynamic';
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
  type LigandPreferences,
  type LigandAnalysis,
  type LigandEvaluation,
  type DesignResult,
  type BinderPreferences,
  type BinderEvaluation,
} from './ai';
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
      <div className="w-full h-[300px] bg-slate-100 rounded-lg flex items-center justify-center">
        <div className="flex items-center gap-2 text-slate-500">
          <div className="w-5 h-5 border-2 border-slate-400 border-t-transparent rounded-full animate-spin" />
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
      <div className="flex-shrink-0 w-8 h-8 rounded-full bg-gradient-to-br from-violet-500 to-purple-600 flex items-center justify-center">
        <span className="material-symbols-outlined text-white text-sm">auto_awesome</span>
      </div>
      <div className="flex-1 space-y-3">
        <div className="bg-slate-50 rounded-2xl rounded-tl-sm px-4 py-3 text-sm text-slate-700 leading-relaxed whitespace-pre-wrap">
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
      <div className="bg-blue-600 text-white rounded-2xl rounded-tr-sm px-4 py-3 text-sm max-w-[80%]">
        {content}
      </div>
      <div className="flex-shrink-0 w-8 h-8 rounded-full bg-slate-200 flex items-center justify-center">
        <span className="material-symbols-outlined text-slate-600 text-sm">person</span>
      </div>
    </div>
  );
}

function SystemMessage({ content }: { content: string }) {
  return (
    <div className="flex justify-center">
      <div className="bg-slate-100 text-slate-600 rounded-full px-4 py-1.5 text-xs flex items-center gap-2">
        <span className="material-symbols-outlined text-sm">info</span>
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
    <div className="bg-gradient-to-br from-purple-50 to-pink-50 rounded-2xl p-6 border border-purple-200 shadow-sm">
      <div className="flex items-center gap-2 mb-4">
        <div className="w-8 h-8 rounded-lg bg-gradient-to-br from-purple-500 to-pink-500 flex items-center justify-center">
          <span className="material-symbols-outlined text-white text-sm">science</span>
        </div>
        <h4 className="font-semibold text-slate-900">Interface Dimer Configuration</h4>
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
          className="flex-1 px-4 py-2.5 text-purple-600 bg-white border border-purple-200 rounded-xl font-medium text-sm hover:bg-purple-50 transition-all disabled:opacity-50"
        >
          Edit
        </button>
        <button
          onClick={onConfirm}
          disabled={isRunning}
          className="flex-1 px-4 py-2.5 bg-gradient-to-r from-purple-600 to-pink-600 text-white rounded-xl font-medium text-sm hover:from-purple-700 hover:to-pink-700 transition-all disabled:opacity-50 flex items-center justify-center gap-2"
        >
          {isRunning ? (
            <>
              <span className="material-symbols-outlined text-sm animate-spin">progress_activity</span>
              Running...
            </>
          ) : (
            <>
              <span className="material-symbols-outlined text-sm">play_arrow</span>
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
    <div className="bg-gradient-to-br from-teal-50 to-emerald-50 rounded-2xl p-6 border border-teal-200 shadow-sm">
      <div className="flex items-center gap-2 mb-4">
        <div className="w-8 h-8 rounded-lg bg-gradient-to-br from-teal-500 to-emerald-500 flex items-center justify-center">
          <span className="material-symbols-outlined text-white text-sm">link</span>
        </div>
        <h4 className="font-semibold text-slate-900">Protein Binder Configuration</h4>
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
          className="flex-1 px-4 py-2.5 text-teal-600 bg-white border border-teal-200 rounded-xl font-medium text-sm hover:bg-teal-50 transition-all disabled:opacity-50"
        >
          Edit
        </button>
        <button
          onClick={onConfirm}
          disabled={isRunning}
          className="flex-1 px-4 py-2.5 bg-gradient-to-r from-teal-600 to-emerald-600 text-white rounded-xl font-medium text-sm hover:from-teal-700 hover:to-emerald-700 transition-all disabled:opacity-50 flex items-center justify-center gap-2"
        >
          {isRunning ? (
            <>
              <span className="material-symbols-outlined text-sm animate-spin">progress_activity</span>
              Running...
            </>
          ) : (
            <>
              <span className="material-symbols-outlined text-sm">play_arrow</span>
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
    <div className={`flex justify-between items-center py-2 ${!isLast ? 'border-b border-slate-100' : ''}`}>
      <span className="text-sm text-slate-600">{label}</span>
      <span className="font-medium text-slate-900">{value}</span>
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
    <div className="bg-white rounded-xl border border-slate-200 overflow-hidden">
      <div className="bg-slate-50 px-4 py-2 border-b border-slate-200 flex items-center justify-between">
        <div className="flex items-center gap-2">
          <span className="material-symbols-outlined text-purple-600">view_in_ar</span>
          <h4 className="font-semibold text-slate-900 text-sm">{title}</h4>
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
          <div className="h-[350px] bg-gradient-to-br from-slate-100 to-slate-50 flex flex-col items-center justify-center text-slate-500">
            <span className="material-symbols-outlined text-4xl mb-2 text-slate-400">view_in_ar</span>
            <p className="text-sm font-medium">{emptyMessage || 'No Structure Available'}</p>
            <p className="text-xs text-slate-400 mt-1">Connect to backend to view actual designs</p>
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

  // Multiple design results state
  const [designResults, setDesignResults] = useState<DesignResult[]>([]);
  const [selectedDesignId, setSelectedDesignId] = useState<string | null>(null);

  // Viewer state
  const [viewerPdbContent, setViewerPdbContent] = useState<string | null>(null);

  const workflowPhase = (aiCaseStudy.workflowPhase || 'idle') as WorkflowPhase;
  const isRunningJob = workflowState.isRunning;

  // Auto-scroll to bottom when new messages arrive
  useEffect(() => {
    messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' });
  }, [aiCaseStudy.conversation, workflowPhase]);

  // Handle structure input (PDB code)
  const handleStructureInput = async () => {
    if (!pdbInput.trim()) return;

    const pdbId = pdbInput.trim().toUpperCase();
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
      <div className="border-b border-slate-100 pb-6 mb-6">
        <div className="flex items-center justify-between mb-2">
          <div className="flex items-center gap-3">
            <span className="bg-gradient-to-r from-violet-500 to-purple-600 text-white text-[11px] font-bold px-3 py-1 rounded-full uppercase tracking-wide flex items-center gap-1.5">
              <span className="material-symbols-outlined text-sm">auto_awesome</span>
              Phase 0
            </span>
            <h2 className="text-xl font-bold text-slate-900">AI Design Assistant</h2>
          </div>
          {workflowPhase !== 'idle' && (
            <button
              onClick={handleStartNew}
              className="text-sm text-slate-500 hover:text-slate-700 flex items-center gap-1"
            >
              <span className="material-symbols-outlined text-sm">refresh</span>
              Start New
            </button>
          )}
        </div>
        <p className="text-slate-500 text-sm leading-relaxed max-w-3xl pl-1">
          {workflowPhase === 'idle'
            ? "Enter a PDB code to begin. I'll analyze the structure and guide you through designing a new metal binding site."
            : "Answer simple questions about your design goals. I'll handle the technical parameters for you."}
        </p>

        {!isBackendConnected && (
          <div className="mt-3 flex items-center gap-2 text-xs text-amber-600 bg-amber-50 px-3 py-1.5 rounded-lg w-fit">
            <span className="material-symbols-outlined text-sm">info</span>
            Demo mode - backend not connected
          </div>
        )}
      </div>

      {/* Quick Start Section */}
      {workflowPhase === 'idle' && aiCaseStudy.conversation.length === 0 && (
        <div className="mb-6 grid grid-cols-1 md:grid-cols-3 gap-4">
          <QuickStartCard
            icon="token"
            title="Metal Binding Redesign"
            description="Convert iron-binding Rubredoxin (1BRF) to bind terbium."
            buttonText="Try 1BRF"
            onClick={handleQuickStart}
            gradient="from-violet-50 to-purple-50"
            borderColor="violet-100"
            iconColor="violet-600"
            buttonGradient="bg-violet-600 hover:bg-violet-700"
          />
          <QuickStartCard
            icon="science"
            title="Interface Dimer Design"
            description="Design separable dimer with azobenzene at the interface."
            buttonText="Try AZOB"
            onClick={() => { setPdbInput('AZOB'); delay(100).then(handleStructureInput); }}
            gradient="from-purple-50 to-pink-50"
            borderColor="purple-200"
            iconColor="pink-600"
            buttonGradient="bg-gradient-to-r from-purple-600 to-pink-600 hover:from-purple-700 hover:to-pink-700"
          />
          <QuickStartCard
            icon="link"
            title="Protein Binder Design"
            description="Design high-affinity binders with ESM-3 validation."
            buttonText="Try BIND"
            onClick={() => { setPdbInput('BIND'); delay(100).then(handleStructureInput); }}
            gradient="from-teal-50 to-emerald-50"
            borderColor="teal-200"
            iconColor="teal-600"
            buttonGradient="bg-gradient-to-r from-teal-600 to-emerald-600 hover:from-teal-700 hover:to-emerald-700"
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
        {workflowPhase === 'interview' && !isLigandDesign && !isBinderDesign && (
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
        {workflowPhase === 'confirming' && userPreferences && !isLigandDesign && !isBinderDesign && (
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

        {/* Job Progress */}
        {workflowPhase === 'running' && (
          <div className="pl-11">
            <JobProgressCard
              jobId={aiCaseStudy.pendingJobId || 'demo-job'}
              status="running"
              progress={aiCaseStudy.jobProgress || workflowState.progress}
              message={workflowState.currentStage || 'Processing...'}
            />
          </div>
        )}

        {/* Results - Metal Binding */}
        {workflowPhase === 'complete' && evaluationResult && userPreferences && !isLigandDesign && !isBinderDesign && (
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
                <div className="flex items-center gap-1 text-xs text-slate-500">
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
              className="w-full px-4 py-2.5 text-purple-600 bg-white border border-purple-200 rounded-xl font-medium text-sm hover:bg-purple-50 transition-all flex items-center justify-center gap-2"
            >
              <span className="material-symbols-outlined text-sm">refresh</span>
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
              <div className="bg-white rounded-xl border border-slate-200 overflow-hidden">
                <div className="bg-gradient-to-r from-teal-50 to-emerald-50 px-4 py-3 border-b border-teal-100">
                  <div className="flex items-center gap-2">
                    <span className="material-symbols-outlined text-teal-600">format_list_numbered</span>
                    <h4 className="font-semibold text-slate-900 text-sm">All Designs ({binderDesigns.length})</h4>
                  </div>
                </div>
                <div className="p-4 space-y-2">
                  {binderDesigns.map((design) => (
                    <button
                      key={design.rank}
                      onClick={() => setSelectedBinderDesign(design)}
                      className={`w-full text-left p-3 rounded-lg border transition-all ${
                        selectedBinderDesign?.rank === design.rank
                          ? 'border-teal-500 bg-teal-50'
                          : 'border-slate-200 hover:border-teal-300 hover:bg-teal-50/50'
                      }`}
                    >
                      <div className="flex items-center justify-between">
                        <div className="flex items-center gap-3">
                          <span className={`w-6 h-6 rounded-full flex items-center justify-center text-xs font-bold ${
                            design.rank === 1 ? 'bg-amber-100 text-amber-700' : 'bg-slate-100 text-slate-600'
                          }`}>
                            {design.rank}
                          </span>
                          <div>
                            <div className="text-sm font-medium text-slate-900">Design #{design.rank}</div>
                            <div className="text-xs text-slate-500">
                              ESM: {design.esm_confidence ? `${(design.esm_confidence * 100).toFixed(0)}%` : 'N/A'} |
                              Contacts: {design.interface_contacts || 'N/A'} |
                              H-bonds: {design.interface_hbonds || 'N/A'}
                            </div>
                          </div>
                        </div>
                        {selectedBinderDesign?.rank === design.rank && (
                          <span className="material-symbols-outlined text-teal-600">check_circle</span>
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
                  <div className="flex items-center gap-1 text-xs text-slate-500">
                    <span className="w-3 h-3 rounded-full bg-slate-400" />
                    <span>Target</span>
                    <span className="ml-2 w-3 h-3 rounded-full bg-violet-400" />
                    <span>Binder</span>
                  </div>
                }
              />
            )}

            <button
              onClick={handleRetry}
              className="w-full px-4 py-2.5 text-teal-600 bg-white border border-teal-200 rounded-xl font-medium text-sm hover:bg-teal-50 transition-all flex items-center justify-center gap-2"
            >
              <span className="material-symbols-outlined text-sm">refresh</span>
              Try Different Settings
            </button>
          </div>
        )}

        {/* Loading indicator */}
        {aiCaseStudy.isProcessing && workflowPhase !== 'interview' && (
          <div className="flex gap-3">
            <div className="flex-shrink-0 w-8 h-8 rounded-full bg-gradient-to-br from-violet-500 to-purple-600 flex items-center justify-center">
              <span className="material-symbols-outlined text-white text-sm animate-spin">progress_activity</span>
            </div>
            <div className="flex items-center gap-2 text-sm text-slate-500">
              <span>
                {workflowPhase === 'fetching' && 'Fetching structure from RCSB...'}
                {workflowPhase === 'analyzing' && !isLigandDesign && !isBinderDesign && 'Analyzing metal binding site...'}
                {workflowPhase === 'analyzing' && isLigandDesign && 'Analyzing ligand topology...'}
                {workflowPhase === 'analyzing' && isBinderDesign && 'Preparing binder design pipeline...'}
                {workflowPhase === 'evaluating' && !isLigandDesign && !isBinderDesign && 'Evaluating design results...'}
                {workflowPhase === 'evaluating' && isLigandDesign && 'Evaluating dimer metrics...'}
                {workflowPhase === 'evaluating' && isBinderDesign && 'Evaluating binder designs...'}
              </span>
            </div>
          </div>
        )}

        <div ref={messagesEndRef} />
      </div>

      {/* Input Area */}
      <div className="border-t border-slate-100 pt-4 mt-auto">
        {(workflowPhase === 'idle' || workflowPhase === 'structure_input') ? (
          <div className="flex gap-3">
            <input
              ref={inputRef}
              type="text"
              value={pdbInput || ''}
              onChange={(e) => setPdbInput(e.target.value.toUpperCase())}
              onKeyDown={(e) => e.key === 'Enter' && handleStructureInput()}
              placeholder="Enter PDB code (e.g., 1BRF) or AZOB for demo"
              className="flex-1 px-4 py-3 bg-slate-50 rounded-xl border border-slate-200 focus:border-violet-500 focus:ring-2 focus:ring-violet-500/20 focus:outline-none text-slate-900 text-sm transition-all font-mono uppercase"
              disabled={aiCaseStudy.isProcessing}
              maxLength={4}
            />
            <button
              onClick={handleStructureInput}
              disabled={!pdbInput.trim() || aiCaseStudy.isProcessing}
              className="px-6 py-3 bg-violet-600 hover:bg-violet-700 disabled:bg-slate-300 disabled:cursor-not-allowed text-white rounded-xl font-medium text-sm transition-all flex items-center gap-2"
            >
              <span className="material-symbols-outlined text-lg">search</span>
              Analyze
            </button>
          </div>
        ) : (
          <div className="flex gap-3">
            <input
              type="text"
              value={followUpInput}
              onChange={(e) => setFollowUpInput(e.target.value)}
              placeholder={workflowPhase === 'complete' ? "Ask a follow-up question..." : "Type a message..."}
              className="flex-1 px-4 py-3 bg-slate-50 rounded-xl border border-slate-200 focus:border-violet-500 focus:ring-2 focus:ring-violet-500/20 focus:outline-none text-slate-900 text-sm transition-all"
              disabled={aiCaseStudy.isProcessing || workflowPhase === 'interview' || workflowPhase === 'running'}
            />
            <button
              disabled={aiCaseStudy.isProcessing || workflowPhase === 'interview' || workflowPhase === 'running' || !followUpInput.trim()}
              className="px-4 py-3 bg-violet-600 hover:bg-violet-700 disabled:bg-slate-300 disabled:cursor-not-allowed text-white rounded-xl font-medium text-sm transition-all flex items-center gap-2"
            >
              <span className="material-symbols-outlined text-lg">send</span>
            </button>
          </div>
        )}
        <p className="mt-2 text-xs text-slate-400 text-center">
          {workflowPhase === 'idle' && "Enter a 4-letter PDB code to start"}
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
  gradient,
  borderColor,
  iconColor,
  buttonGradient,
}: {
  icon: string;
  title: string;
  description: string;
  buttonText: string;
  onClick: () => void;
  gradient: string;
  borderColor: string;
  iconColor: string;
  buttonGradient: string;
}) {
  return (
    <div className={`p-5 bg-gradient-to-br ${gradient} rounded-2xl border border-${borderColor}`}>
      <div className="flex items-start gap-3">
        <div className="p-2.5 bg-white rounded-xl shadow-sm">
          <span className={`material-symbols-outlined text-xl text-${iconColor}`}>{icon}</span>
        </div>
        <div className="flex-1">
          <h3 className="font-bold text-slate-900 mb-1 text-sm">{title}</h3>
          <p className="text-xs text-slate-600 mb-3">{description}</p>
          <button
            onClick={onClick}
            className={`inline-flex items-center gap-1.5 px-3 py-1.5 ${buttonGradient} text-white rounded-lg font-medium text-xs transition-all`}
          >
            <span className="material-symbols-outlined text-sm">play_arrow</span>
            {buttonText}
          </button>
        </div>
      </div>
    </div>
  );
}
