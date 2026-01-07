'use client';

import { useStore } from '@/lib/store';
import { Dna, Atom, FlaskConical, Check, ChevronRight } from 'lucide-react';

interface WorkflowStep {
  id: 'rfd3' | 'rf3' | 'mpnn';
  label: string;
  shortLabel: string;
  description: string;
  icon: typeof Dna;
  color: string;
  bgColor: string;
  optional?: boolean;
}

const WORKFLOW_STEPS: WorkflowStep[] = [
  {
    id: 'rfd3',
    label: 'Design Structure',
    shortLabel: 'RFD3',
    description: 'Generate protein backbone',
    icon: Dna,
    color: 'text-blue-400',
    bgColor: 'bg-blue-600',
  },
  {
    id: 'rf3',
    label: 'Validate Fold',
    shortLabel: 'RF3',
    description: 'Predict & validate structure',
    icon: Atom,
    color: 'text-green-400',
    bgColor: 'bg-green-600',
    optional: true,
  },
  {
    id: 'mpnn',
    label: 'Design Sequences',
    shortLabel: 'MPNN',
    description: 'Generate amino acid sequences',
    icon: FlaskConical,
    color: 'text-purple-400',
    bgColor: 'bg-purple-600',
  },
];

export function WorkflowStepper() {
  const { jobs, activeTab, setActiveTab, lastCompletedJobType } = useStore();

  // Determine step status based on completed jobs
  const getStepStatus = (stepId: 'rfd3' | 'rf3' | 'mpnn'): 'completed' | 'current' | 'upcoming' | 'suggested' => {
    const hasCompletedJob = jobs.some(
      (job) => job.type === stepId && job.status === 'completed'
    );

    if (hasCompletedJob) return 'completed';
    if (activeTab === stepId) return 'current';

    // Suggest next step based on workflow
    if (lastCompletedJobType === 'rfd3' && stepId === 'mpnn') return 'suggested';
    if (lastCompletedJobType === 'rf3' && stepId === 'mpnn') return 'suggested';
    if (lastCompletedJobType === 'mpnn' && stepId === 'rf3') return 'suggested';

    return 'upcoming';
  };

  return (
    <div className="bg-gray-800/50 border border-gray-700 rounded-lg p-4">
      <div className="flex items-center justify-between">
        {WORKFLOW_STEPS.map((step, index) => {
          const status = getStepStatus(step.id);
          const isLast = index === WORKFLOW_STEPS.length - 1;

          return (
            <div key={step.id} className="flex items-center flex-1">
              {/* Step */}
              <button
                onClick={() => setActiveTab(step.id)}
                className={`flex items-center gap-3 px-4 py-2 rounded-lg transition group ${
                  status === 'current'
                    ? 'bg-gray-700 ring-2 ring-offset-2 ring-offset-gray-800 ring-gray-500'
                    : status === 'completed'
                    ? 'bg-gray-700/50 hover:bg-gray-700'
                    : status === 'suggested'
                    ? 'bg-gray-700/30 hover:bg-gray-700 ring-1 ring-blue-500/50'
                    : 'hover:bg-gray-700/50'
                }`}
              >
                {/* Step Icon/Number */}
                <div
                  className={`w-8 h-8 rounded-full flex items-center justify-center transition ${
                    status === 'completed'
                      ? 'bg-green-600 text-white'
                      : status === 'current'
                      ? step.bgColor + ' text-white'
                      : status === 'suggested'
                      ? 'bg-gray-600 text-white ring-1 ring-blue-400'
                      : 'bg-gray-700 text-gray-400'
                  }`}
                >
                  {status === 'completed' ? (
                    <Check className="w-4 h-4" />
                  ) : (
                    <step.icon className="w-4 h-4" />
                  )}
                </div>

                {/* Step Label */}
                <div className="text-left">
                  <div className="flex items-center gap-2">
                    <span
                      className={`text-sm font-medium ${
                        status === 'current' || status === 'completed'
                          ? 'text-white'
                          : status === 'suggested'
                          ? 'text-blue-300'
                          : 'text-gray-400'
                      }`}
                    >
                      {step.label}
                    </span>
                    {step.optional && (
                      <span className="text-xs text-gray-500 bg-gray-700 px-1.5 py-0.5 rounded">
                        optional
                      </span>
                    )}
                    {status === 'suggested' && (
                      <span className="text-xs text-blue-400 bg-blue-900/30 px-1.5 py-0.5 rounded">
                        next
                      </span>
                    )}
                  </div>
                  <p className="text-xs text-gray-500 hidden sm:block">
                    {step.description}
                  </p>
                </div>
              </button>

              {/* Connector */}
              {!isLast && (
                <div className="flex-1 flex items-center justify-center px-2">
                  <div
                    className={`h-0.5 w-full ${
                      getStepStatus(WORKFLOW_STEPS[index + 1].id) === 'completed' ||
                      status === 'completed'
                        ? 'bg-green-600'
                        : 'bg-gray-700'
                    }`}
                  />
                  <ChevronRight
                    className={`w-4 h-4 flex-shrink-0 ${
                      status === 'completed' ? 'text-green-500' : 'text-gray-600'
                    }`}
                  />
                </div>
              )}
            </div>
          );
        })}
      </div>

      {/* Workflow hint */}
      {!jobs.some((j) => j.status === 'completed') && (
        <p className="text-xs text-gray-500 mt-3 text-center">
          Start with RFD3 to design a protein backbone, then optionally validate with RF3, and design sequences with MPNN
        </p>
      )}
    </div>
  );
}
