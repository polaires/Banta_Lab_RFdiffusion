'use client';

import { useEffect, useState } from 'react';
import {
  MessageSquare,
  FlaskConical,
  Building2,
  Type,
  ShieldCheck,
  FileText,
  CheckCircle,
  AlertCircle,
  Check,
  X,
  Loader2,
  Beaker,
  Atom,
  ArrowRight,
  ChevronDown,
  ChevronUp
} from 'lucide-react';
import { cn } from '@/lib/utils';

// Pipeline stages for AI Design
export type AIDesignStage =
  | 'idle'
  | 'parsing'
  | 'resolving'
  | 'configuring'
  | 'backbone'
  | 'sequence'
  | 'validation'
  | 'analysis'
  | 'complete'
  | 'error';

interface StageInfo {
  metal?: string;
  ligand?: string;
  confidence?: number;
  numDesigns?: number;
  numSequences?: number;
  passRate?: number;
  bestRmsd?: number;
  bestPlddt?: number;
  message?: string;
}

interface AIDesignPipelineWorkflowProps {
  currentStage: AIDesignStage;
  stageInfo?: StageInfo;
  error?: string;
  query?: string;
  showDetails?: boolean;
}

// Define all pipeline stages
const PIPELINE_STAGES = [
  {
    id: 'parsing',
    name: 'Parse Query',
    Icon: MessageSquare,
    description: 'Understanding your request',
    phase: 0,
  },
  {
    id: 'resolving',
    name: 'Resolve Ligand',
    Icon: Beaker,
    description: 'Loading structure & chemistry',
    phase: 1,
  },
  {
    id: 'configuring',
    name: 'Configure',
    Icon: FlaskConical,
    description: 'Setting up RFD3 parameters',
    phase: 1,
  },
  {
    id: 'backbone',
    name: 'Backbone',
    Icon: Building2,
    description: 'RFD3 generating structure',
    phase: 2,
  },
  {
    id: 'sequence',
    name: 'Sequence',
    Icon: Type,
    description: 'LigandMPNN design',
    phase: 3,
  },
  {
    id: 'validation',
    name: 'Validate',
    Icon: ShieldCheck,
    description: 'ESMFold + composition',
    phase: 4,
  },
  {
    id: 'analysis',
    name: 'Analyze',
    Icon: Atom,
    description: 'Contact & metrics',
    phase: 4,
  },
];

// Phase descriptions matching the user's diagram
const PHASES = [
  { id: 0, name: 'Input', description: 'Natural language parsing' },
  { id: 1, name: 'Chemistry', description: 'Ligand analysis & configuration' },
  { id: 2, name: 'Backbone', description: 'RFD3 structure generation' },
  { id: 3, name: 'Sequence', description: 'LigandMPNN design' },
  { id: 4, name: 'Validation', description: 'Structure & metrics' },
];

export function AIDesignPipelineWorkflow({
  currentStage,
  stageInfo,
  error,
  query,
  showDetails = true
}: AIDesignPipelineWorkflowProps) {
  const [expanded, setExpanded] = useState(true);
  const [animatedStage, setAnimatedStage] = useState<string | null>(null);

  // Get current stage index
  const currentIndex = PIPELINE_STAGES.findIndex(s => s.id === currentStage);
  const isComplete = currentStage === 'complete';
  const isError = currentStage === 'error';
  const isIdle = currentStage === 'idle';

  // Animate current stage
  useEffect(() => {
    if (!isIdle && !isComplete && !isError) {
      setAnimatedStage(currentStage);
    } else {
      setAnimatedStage(null);
    }
  }, [currentStage, isIdle, isComplete, isError]);

  const getStageStatus = (stageId: string, index: number) => {
    if (isError) return index <= currentIndex ? 'error' : 'pending';
    if (isComplete) return 'complete';
    if (index < currentIndex) return 'complete';
    if (index === currentIndex) return 'active';
    return 'pending';
  };

  const getStageStyles = (status: string) => {
    switch (status) {
      case 'complete':
        return {
          icon: 'bg-emerald-100 text-emerald-600 dark:bg-emerald-900/30 dark:text-emerald-400',
          line: 'bg-emerald-500',
          text: 'text-foreground',
          desc: 'text-emerald-600 dark:text-emerald-400'
        };
      case 'active':
        return {
          icon: 'bg-blue-100 text-blue-600 ring-2 ring-blue-400/50 ring-offset-2 dark:bg-blue-900/30 dark:text-blue-400',
          line: 'bg-muted',
          text: 'text-foreground font-medium',
          desc: 'text-blue-600 dark:text-blue-400'
        };
      case 'error':
        return {
          icon: 'bg-red-100 text-red-600 dark:bg-red-900/30 dark:text-red-400',
          line: 'bg-red-300',
          text: 'text-red-700 dark:text-red-400',
          desc: 'text-red-500'
        };
      default:
        return {
          icon: 'bg-muted text-muted-foreground',
          line: 'bg-muted',
          text: 'text-muted-foreground',
          desc: 'text-muted-foreground'
        };
    }
  };

  // Get current phase
  const currentPhase = isComplete ? 4 :
    isError ? -1 :
    isIdle ? -1 :
    PIPELINE_STAGES.find(s => s.id === currentStage)?.phase ?? 0;

  return (
    <div className="w-full bg-card border border-border rounded-xl p-4 shadow-sm">
      {/* Header */}
      <div
        className="flex items-center justify-between cursor-pointer"
        onClick={() => setExpanded(!expanded)}
      >
        <div className="flex items-center gap-3">
          <div className="p-2 rounded-lg bg-gradient-to-br from-blue-500 to-purple-600">
            <FlaskConical className="w-5 h-5 text-white" />
          </div>
          <div>
            <h3 className="font-semibold text-foreground">AI Design Pipeline</h3>
            {query && (
              <p className="text-sm text-muted-foreground truncate max-w-[300px]">
                "{query}"
              </p>
            )}
          </div>
        </div>

        <div className="flex items-center gap-3">
          {isComplete && (
            <span className="text-sm text-emerald-600 dark:text-emerald-400 flex items-center gap-1">
              <CheckCircle className="w-4 h-4" />
              Complete
            </span>
          )}
          {isError && (
            <span className="text-sm text-red-600 flex items-center gap-1">
              <AlertCircle className="w-4 h-4" />
              Failed
            </span>
          )}
          {!isComplete && !isError && !isIdle && (
            <span className="text-sm text-blue-600 dark:text-blue-400 flex items-center gap-1">
              <Loader2 className="w-4 h-4 animate-spin" />
              Running
            </span>
          )}
          {expanded ? (
            <ChevronUp className="w-5 h-5 text-muted-foreground" />
          ) : (
            <ChevronDown className="w-5 h-5 text-muted-foreground" />
          )}
        </div>
      </div>

      {expanded && (
        <>
          {/* Phase Progress Bar */}
          <div className="mt-4 mb-6">
            <div className="flex justify-between mb-2">
              {PHASES.map((phase, index) => (
                <div
                  key={phase.id}
                  className={cn(
                    "text-xs font-medium",
                    index <= currentPhase ? "text-foreground" : "text-muted-foreground"
                  )}
                >
                  {phase.name}
                </div>
              ))}
            </div>
            <div className="h-2 bg-muted rounded-full overflow-hidden">
              <div
                className="h-full bg-gradient-to-r from-blue-500 to-emerald-500 transition-all duration-500 rounded-full"
                style={{
                  width: isComplete ? '100%' : `${Math.max(0, (currentPhase + 1) / PHASES.length * 100)}%`
                }}
              />
            </div>
          </div>

          {/* Pipeline Stages Grid */}
          <div className="grid grid-cols-7 gap-2">
            {PIPELINE_STAGES.map((stage, index) => {
              const status = getStageStatus(stage.id, index);
              const styles = getStageStyles(status);
              const isAnimating = animatedStage === stage.id;
              const StageIcon = stage.Icon;

              return (
                <div key={stage.id} className="flex flex-col items-center">
                  {/* Icon Circle */}
                  <div
                    className={cn(
                      "w-10 h-10 rounded-full flex items-center justify-center transition-all duration-300",
                      styles.icon,
                      isAnimating && "animate-pulse"
                    )}
                  >
                    {status === 'complete' ? (
                      <Check className="w-4 h-4" />
                    ) : status === 'error' && index === currentIndex ? (
                      <X className="w-4 h-4" />
                    ) : isAnimating ? (
                      <Loader2 className="w-4 h-4 animate-spin" />
                    ) : (
                      <StageIcon className="w-4 h-4" />
                    )}
                  </div>

                  {/* Stage Name */}
                  <span className={cn("mt-1 text-xs text-center", styles.text)}>
                    {stage.name}
                  </span>

                  {/* Arrow connector (except last) */}
                  {index < PIPELINE_STAGES.length - 1 && (
                    <ArrowRight
                      className={cn(
                        "absolute w-3 h-3 -right-2.5 top-3",
                        index < currentIndex ? "text-emerald-500" : "text-muted-foreground/30"
                      )}
                      style={{ display: 'none' }} // Hidden for grid layout
                    />
                  )}
                </div>
              );
            })}
          </div>

          {/* Stage Details */}
          {showDetails && stageInfo && (
            <div className="mt-4 pt-4 border-t border-border space-y-2">
              {/* Parsed Intent */}
              {(stageInfo.metal || stageInfo.ligand) && (
                <div className="flex items-center gap-4 text-sm">
                  <span className="text-muted-foreground">Detected:</span>
                  {stageInfo.metal && (
                    <span className="px-2 py-0.5 rounded-full bg-purple-100 text-purple-700 dark:bg-purple-900/30 dark:text-purple-300 text-xs font-medium">
                      Metal: {stageInfo.metal}
                    </span>
                  )}
                  {stageInfo.ligand && (
                    <span className="px-2 py-0.5 rounded-full bg-blue-100 text-blue-700 dark:bg-blue-900/30 dark:text-blue-300 text-xs font-medium">
                      Ligand: {stageInfo.ligand}
                    </span>
                  )}
                </div>
              )}

              {/* Design Counts */}
              {(stageInfo.numDesigns || stageInfo.numSequences) && (
                <div className="flex items-center gap-4 text-sm">
                  <span className="text-muted-foreground">Generated:</span>
                  {stageInfo.numDesigns && (
                    <span className="text-foreground">
                      {stageInfo.numDesigns} backbones
                    </span>
                  )}
                  {stageInfo.numSequences && (
                    <span className="text-foreground">
                      {stageInfo.numSequences} sequences
                    </span>
                  )}
                </div>
              )}

              {/* Validation Results */}
              {stageInfo.passRate !== undefined && (
                <div className="flex items-center gap-4 text-sm">
                  <span className="text-muted-foreground">Validation:</span>
                  <span className={cn(
                    "font-medium",
                    stageInfo.passRate >= 0.2 ? "text-emerald-600" : "text-amber-600"
                  )}>
                    {(stageInfo.passRate * 100).toFixed(0)}% pass rate
                  </span>
                  {stageInfo.bestRmsd && (
                    <span className="text-foreground">
                      Best RMSD: {stageInfo.bestRmsd.toFixed(2)}Å
                    </span>
                  )}
                  {stageInfo.bestPlddt && (
                    <span className="text-foreground">
                      pLDDT: {stageInfo.bestPlddt.toFixed(2)}
                    </span>
                  )}
                </div>
              )}

              {/* Current Status Message */}
              {stageInfo.message && (
                <div className="flex items-center gap-2 text-sm text-muted-foreground">
                  <Loader2 className="w-3 h-3 animate-spin" />
                  {stageInfo.message}
                </div>
              )}
            </div>
          )}

          {/* Error Message */}
          {isError && error && (
            <div className="mt-4 p-3 bg-red-50 dark:bg-red-900/20 border border-red-200 dark:border-red-800 rounded-lg">
              <div className="flex items-start gap-2">
                <AlertCircle className="w-4 h-4 text-red-500 flex-shrink-0 mt-0.5" />
                <p className="text-sm text-red-700 dark:text-red-300">{error}</p>
              </div>
            </div>
          )}
        </>
      )}
    </div>
  );
}

// Compact version for inline display
export function AIDesignPipelineStatus({
  currentStage,
  stageInfo
}: Pick<AIDesignPipelineWorkflowProps, 'currentStage' | 'stageInfo'>) {
  const currentIndex = PIPELINE_STAGES.findIndex(s => s.id === currentStage);
  const isComplete = currentStage === 'complete';
  const isError = currentStage === 'error';

  const progress = isComplete ? 100 :
    isError ? 0 :
    Math.max(0, (currentIndex + 1) / PIPELINE_STAGES.length * 100);

  return (
    <div className="flex items-center gap-3">
      <div className="flex-1 h-1.5 bg-muted rounded-full overflow-hidden">
        <div
          className={cn(
            "h-full transition-all duration-300 rounded-full",
            isError ? "bg-red-500" : "bg-gradient-to-r from-blue-500 to-emerald-500"
          )}
          style={{ width: `${progress}%` }}
        />
      </div>
      <span className="text-xs text-muted-foreground whitespace-nowrap">
        {isComplete ? 'Complete' :
         isError ? 'Failed' :
         PIPELINE_STAGES[currentIndex]?.name || 'Ready'}
      </span>
    </div>
  );
}

export default AIDesignPipelineWorkflow;
