'use client';

import { useEffect, useState } from 'react';
import { GitBranch, CheckCircle, AlertCircle, Check, X, Loader2, Building2, Type, ShieldCheck, FlaskConical } from 'lucide-react';

export type PipelineStage = 'idle' | 'backbone' | 'ligandmpnn' | 'validation' | 'analysis' | 'complete' | 'error';

interface StageDetails {
  backbone?: { progress?: number; message?: string };
  ligandmpnn?: { sequences?: number; message?: string };
  validation?: { message?: string };
  analysis?: { message?: string };
}

interface PipelineProgressProps {
  currentStage: PipelineStage;
  stageDetails?: StageDetails;
  error?: string;
}

const STAGES = [
  {
    id: 'backbone',
    name: 'Backbone',
    Icon: Building2,
    description: 'RFD3 generating structure',
    color: 'purple'
  },
  {
    id: 'ligandmpnn',
    name: 'LigandMPNN',
    Icon: Type,
    description: 'Sequence design',
    color: 'blue'
  },
  {
    id: 'validation',
    name: 'Validation',
    Icon: ShieldCheck,
    description: 'ESM + structure checks',
    color: 'teal'
  },
  {
    id: 'analysis',
    name: 'Analysis',
    Icon: FlaskConical,
    description: 'PLIP interaction profiling',
    color: 'emerald'
  },
];

export function PipelineProgress({ currentStage, stageDetails, error }: PipelineProgressProps) {
  const [animatedStage, setAnimatedStage] = useState<string | null>(null);

  // Get index of current stage
  const currentIndex = STAGES.findIndex(s => s.id === currentStage);
  const isComplete = currentStage === 'complete';
  const isError = currentStage === 'error';

  // Animate current stage
  useEffect(() => {
    if (currentStage !== 'idle' && currentStage !== 'complete' && currentStage !== 'error') {
      setAnimatedStage(currentStage);
    } else {
      setAnimatedStage(null);
    }
  }, [currentStage]);

  const getStageStatus = (stageId: string, index: number) => {
    if (isError) return index <= currentIndex ? 'error' : 'pending';
    if (isComplete) return 'complete';
    if (index < currentIndex) return 'complete';
    if (index === currentIndex) return 'active';
    return 'pending';
  };

  const getStageStyles = (status: string, color: string) => {
    switch (status) {
      case 'complete':
        return {
          icon: 'bg-emerald-100 text-emerald-600',
          line: 'bg-emerald-400',
          text: 'text-slate-900',
          desc: 'text-emerald-600'
        };
      case 'active':
        return {
          icon: `bg-${color}-100 text-${color}-600 ring-2 ring-${color}-400 ring-offset-2`,
          line: 'bg-slate-200',
          text: 'text-slate-900 font-medium',
          desc: `text-${color}-600`
        };
      case 'error':
        return {
          icon: 'bg-red-100 text-red-600',
          line: 'bg-red-200',
          text: 'text-red-700',
          desc: 'text-red-500'
        };
      default:
        return {
          icon: 'bg-slate-100 text-slate-400',
          line: 'bg-slate-200',
          text: 'text-slate-400',
          desc: 'text-slate-400'
        };
    }
  };

  return (
    <div className="w-full py-6">
      {/* Header */}
      <div className="flex items-center gap-2 mb-6">
        <GitBranch className="w-5 h-5 text-purple-600" />
        <span className="font-medium text-slate-700">Design Pipeline</span>
        {isComplete && (
          <span className="ml-auto text-sm text-emerald-600 flex items-center gap-1">
            <CheckCircle className="w-4 h-4" />
            Complete
          </span>
        )}
        {isError && (
          <span className="ml-auto text-sm text-red-600 flex items-center gap-1">
            <AlertCircle className="w-4 h-4" />
            Failed
          </span>
        )}
      </div>

      {/* Pipeline Steps */}
      <div className="relative">
        {/* Connection Lines */}
        <div className="absolute top-6 left-6 right-6 h-0.5 bg-slate-200 -z-10" />

        {/* Progress Line */}
        <div
          className="absolute top-6 left-6 h-0.5 bg-gradient-to-r from-purple-500 via-blue-500 to-emerald-500 transition-all duration-500 -z-10"
          style={{
            width: isComplete
              ? 'calc(100% - 48px)'
              : `calc(${Math.max(0, currentIndex) / (STAGES.length - 1) * 100}% - 24px)`
          }}
        />

        {/* Stages */}
        <div className="flex justify-between">
          {STAGES.map((stage, index) => {
            const status = getStageStatus(stage.id, index);
            const styles = getStageStyles(status, stage.color);
            const isAnimating = animatedStage === stage.id;
            const details = stageDetails?.[stage.id as keyof StageDetails];
            const StageIcon = stage.Icon;

            return (
              <div key={stage.id} className="flex flex-col items-center w-24">
                {/* Icon Circle */}
                <div
                  className={`
                    w-12 h-12 rounded-full flex items-center justify-center
                    transition-all duration-300
                    ${styles.icon}
                    ${isAnimating ? 'animate-pulse' : ''}
                  `}
                >
                  {status === 'complete' ? (
                    <Check className="w-5 h-5" />
                  ) : status === 'error' && index === currentIndex ? (
                    <X className="w-5 h-5" />
                  ) : isAnimating ? (
                    <Loader2 className="w-5 h-5 animate-spin" />
                  ) : (
                    <StageIcon className="w-5 h-5" />
                  )}
                </div>

                {/* Stage Name */}
                <span className={`mt-2 text-sm font-medium ${styles.text}`}>
                  {stage.name}
                </span>

                {/* Description or Status */}
                <span className={`text-xs text-center ${styles.desc}`}>
                  {status === 'active' && details?.message ? (
                    details.message
                  ) : status === 'complete' ? (
                    'Done'
                  ) : (
                    stage.description
                  )}
                </span>
              </div>
            );
          })}
        </div>
      </div>

      {/* Error Message */}
      {isError && error && (
        <div className="mt-4 p-3 bg-red-50 border border-red-200 rounded-lg">
          <div className="flex items-start gap-2">
            <AlertCircle className="w-4 h-4 text-red-500 flex-shrink-0 mt-0.5" />
            <p className="text-sm text-red-700">{error}</p>
          </div>
        </div>
      )}

      {/* Current Stage Details */}
      {currentStage !== 'idle' && currentStage !== 'complete' && currentStage !== 'error' && (
        <div className="mt-4 p-3 bg-slate-50 rounded-lg border border-slate-200">
          <div className="flex items-center gap-2">
            <Loader2 className="w-4 h-4 text-slate-400 animate-spin" />
            <span className="text-sm text-slate-600">
              {stageDetails?.[currentStage as keyof StageDetails]?.message ||
               `Running ${STAGES.find(s => s.id === currentStage)?.name}...`}
            </span>
          </div>
        </div>
      )}
    </div>
  );
}

export default PipelineProgress;
