'use client';

import { useStore } from '@/lib/store';
import { Dna, FlaskConical, Atom, History, Server, Check } from 'lucide-react';

const workflowSteps = [
  { id: 'rfd3' as const, label: 'Design Structure', shortLabel: 'RFD3', icon: Dna },
  { id: 'mpnn' as const, label: 'Design Sequences', shortLabel: 'MPNN', icon: FlaskConical },
  { id: 'rf3' as const, label: 'Validate Fold', shortLabel: 'RF3', icon: Atom },
];

export function Header() {
  const {
    activeTab,
    setActiveTab,
    health,
    jobs,
    setConnectionModalOpen,
    latestRfd3Design,
    lastCompletedJobType,
  } = useStore();

  const isConnected = health?.status === 'healthy';
  const pendingJobs = jobs.filter(j => j.status === 'pending' || j.status === 'running').length;

  // Determine step completion status based on workflow progress
  const getStepStatus = (stepId: string, index: number) => {
    if (activeTab === stepId) return 'active';

    // Check if step has been completed
    if (stepId === 'rfd3' && latestRfd3Design) return 'completed';
    if (stepId === 'mpnn' && lastCompletedJobType === 'mpnn') return 'completed';
    if (stepId === 'rf3' && lastCompletedJobType === 'rf3') return 'completed';

    return 'upcoming';
  };

  return (
    <header className="sticky top-0 z-50 border-b border-gray-200 bg-white/95 backdrop-blur-sm shadow-sm">
      <div className="max-w-7xl mx-auto px-4 py-3 flex items-center justify-between">
        {/* Logo & Title (Left) */}
        <div className="flex items-center gap-3 min-w-[180px]">
          <Dna className="w-7 h-7 text-blue-600" />
          <div>
            <h1 className="text-lg font-bold leading-tight text-gray-900">Foundry Protein Design</h1>
            <p className="text-xs text-gray-500 hidden sm:block">IPD Design Pipeline</p>
          </div>
        </div>

        {/* Workflow Stepper (Center) */}
        <nav className="hidden md:flex items-center gap-1">
          {workflowSteps.map((step, index) => {
            const status = getStepStatus(step.id, index);
            const StepIcon = step.icon;

            return (
              <div key={step.id} className="flex items-center">
                {/* Connector line (before step, except first) */}
                {index > 0 && (
                  <div
                    className={`w-8 h-0.5 ${
                      status === 'completed' || getStepStatus(workflowSteps[index - 1].id, index - 1) === 'completed'
                        ? 'bg-green-500'
                        : 'bg-gray-300'
                    }`}
                  />
                )}

                {/* Step button */}
                <button
                  onClick={() => setActiveTab(step.id)}
                  className={`flex items-center gap-2 px-3 py-2 rounded-lg transition-all ${
                    status === 'active'
                      ? 'bg-blue-600 text-white'
                      : status === 'completed'
                      ? 'bg-green-100 text-green-700 hover:bg-green-200'
                      : 'text-gray-500 hover:bg-gray-100 hover:text-gray-700'
                  }`}
                  title={step.label}
                >
                  {/* Step circle/icon */}
                  <div className={`w-6 h-6 rounded-full flex items-center justify-center text-xs font-bold ${
                    status === 'active'
                      ? 'bg-white/20'
                      : status === 'completed'
                      ? 'bg-green-500 text-white'
                      : 'bg-gray-200 text-gray-600'
                  }`}>
                    {status === 'completed' ? (
                      <Check className="w-3.5 h-3.5" />
                    ) : (
                      index + 1
                    )}
                  </div>

                  {/* Step label */}
                  <div className="text-left hidden lg:block">
                    <div className="text-xs font-medium">{step.shortLabel}</div>
                    <div className="text-[10px] opacity-70">{step.label}</div>
                  </div>
                </button>
              </div>
            );
          })}
        </nav>

        {/* Mobile step indicator */}
        <div className="md:hidden flex items-center gap-2">
          {workflowSteps.map((step, index) => {
            const status = getStepStatus(step.id, index);
            return (
              <button
                key={step.id}
                onClick={() => setActiveTab(step.id)}
                className={`w-8 h-8 rounded-full flex items-center justify-center text-xs font-bold transition ${
                  status === 'active'
                    ? 'bg-blue-600 text-white'
                    : status === 'completed'
                    ? 'bg-green-500 text-white'
                    : 'bg-gray-200 text-gray-600'
                }`}
              >
                {status === 'completed' ? <Check className="w-4 h-4" /> : index + 1}
              </button>
            );
          })}
        </div>

        {/* Icon buttons (Right) */}
        <div className="flex items-center gap-2 min-w-[100px] justify-end">
          {/* Jobs button */}
          <button
            onClick={() => setActiveTab('jobs')}
            className={`relative p-2 rounded-lg transition ${
              activeTab === 'jobs'
                ? 'bg-gray-200 text-gray-900'
                : 'text-gray-500 hover:bg-gray-100 hover:text-gray-700'
            }`}
            title="Job History"
          >
            <History className="w-5 h-5" />
            {pendingJobs > 0 && (
              <span className="absolute -top-1 -right-1 w-4 h-4 bg-blue-500 text-white text-[10px] font-bold rounded-full flex items-center justify-center">
                {pendingJobs}
              </span>
            )}
          </button>

          {/* Connection button */}
          <button
            onClick={() => setConnectionModalOpen(true)}
            className="relative p-2 rounded-lg text-gray-500 hover:bg-gray-100 hover:text-gray-700 transition"
            title={isConnected ? 'Connected to backend' : 'Connect to backend'}
          >
            <Server className="w-5 h-5" />
            {/* Status dot */}
            <span
              className={`absolute -top-0.5 -right-0.5 w-2.5 h-2.5 rounded-full border-2 border-white ${
                isConnected ? 'bg-green-500' : 'bg-red-500'
              }`}
            />
          </button>
        </div>
      </div>
    </header>
  );
}
