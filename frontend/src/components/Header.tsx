'use client';

import { useStore, TabId } from '@/lib/store';
import { UserMenu } from './auth/UserMenu';

const workflowSteps = [
  { id: 'task' as const, label: 'Design Task', step: 0 },
  { id: 'rfd3' as const, label: 'RFdiffusion3', step: 1 },
  { id: 'mpnn' as const, label: 'MPNN', step: 2 },
  { id: 'rf3' as const, label: 'Validate', step: 3 },
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
    selectedDesignTask,
  } = useStore();

  const isConnected = health?.status === 'healthy';
  const pendingJobs = jobs.filter(j => j.status === 'pending' || j.status === 'running').length;

  const getStepStatus = (stepId: string) => {
    if (activeTab === stepId) return 'active';
    if (stepId === 'task' && selectedDesignTask) return 'completed';
    if (stepId === 'rfd3' && latestRfd3Design) return 'completed';
    if (stepId === 'mpnn' && lastCompletedJobType === 'mpnn') return 'completed';
    if (stepId === 'rf3' && lastCompletedJobType === 'rf3') return 'completed';
    return 'upcoming';
  };

  return (
    <header className="bg-white/80 backdrop-blur-md border-b border-slate-200 sticky top-0 z-50">
      <div className="max-w-6xl mx-auto px-6 h-16 flex items-center justify-between">
        {/* Logo & Brand */}
        <div className="flex items-center gap-3">
          <div className="text-blue-600 bg-blue-50 rounded-lg p-1.5 flex items-center justify-center">
            <span className="material-symbols-outlined text-[24px]">biotech</span>
          </div>
          <div>
            <h1 className="font-bold text-sm text-slate-900 leading-tight">Foundry Protein Design</h1>
            <p className="text-[10px] text-slate-500 font-medium uppercase tracking-wide">IPD Design Pipeline</p>
          </div>
        </div>

        {/* Workflow Stepper */}
        <nav className="hidden md:flex items-center space-x-1">
          {workflowSteps.map((step, index) => {
            const status = getStepStatus(step.id);
            const isActive = status === 'active';

            return (
              <div key={step.id} className="flex items-center">
                {index > 0 && <div className="w-8 h-px bg-slate-200" />}

                <button
                  onClick={() => setActiveTab(step.id)}
                  className={`flex items-center gap-2 px-3 py-1.5 rounded-full transition-all ${
                    isActive
                      ? 'bg-white border border-blue-200 shadow-sm'
                      : 'opacity-50 grayscale hover:opacity-75'
                  }`}
                >
                  <span className={`w-5 h-5 flex items-center justify-center rounded-full text-[10px] font-bold ${
                    isActive
                      ? 'bg-blue-600 text-white'
                      : 'bg-slate-200 text-slate-600'
                  }`}>
                    {step.step}
                  </span>
                  <span className={`text-xs font-semibold ${
                    isActive ? 'text-blue-900' : 'text-slate-600'
                  }`}>
                    {step.label}
                  </span>
                </button>
              </div>
            );
          })}
        </nav>

        {/* Mobile step indicator */}
        <div className="md:hidden flex items-center gap-1">
          {workflowSteps.map((step) => {
            const status = getStepStatus(step.id);
            const isActive = status === 'active';
            return (
              <button
                key={step.id}
                onClick={() => setActiveTab(step.id)}
                className={`w-7 h-7 rounded-full flex items-center justify-center text-xs font-bold transition ${
                  isActive
                    ? 'bg-blue-600 text-white'
                    : 'bg-slate-200 text-slate-600'
                }`}
              >
                {step.step}
              </button>
            );
          })}
        </div>

        {/* Action buttons */}
        <div className="flex items-center gap-2">
          {/* History button */}
          <button
            onClick={() => setActiveTab('jobs')}
            className={`p-2 rounded-full transition-all relative ${
              activeTab === 'jobs'
                ? 'text-slate-700 bg-slate-100'
                : 'text-slate-400 hover:text-slate-700 hover:bg-slate-100'
            }`}
            title="Job History"
          >
            <span className="material-symbols-outlined text-xl">history</span>
            {pendingJobs > 0 && (
              <span className="absolute -top-0.5 -right-0.5 w-4 h-4 bg-blue-500 text-white text-[9px] font-bold rounded-full flex items-center justify-center border border-white">
                {pendingJobs}
              </span>
            )}
          </button>

          <div className="h-4 w-px bg-slate-200 mx-1" />

          {/* Connection status */}
          <button
            onClick={() => setConnectionModalOpen(true)}
            className="p-2 rounded-full text-slate-400 hover:text-slate-700 hover:bg-slate-100 transition-all relative"
            title={isConnected ? 'Connected to backend' : 'Connect to backend'}
          >
            <span className="material-symbols-outlined text-xl">dns</span>
            <span className={`absolute top-2 right-2 w-2 h-2 rounded-full border border-white ${
              isConnected ? 'bg-emerald-500' : 'bg-red-500'
            }`} />
          </button>

          <div className="h-4 w-px bg-slate-200 mx-1" />

          {/* User menu */}
          <UserMenu />
        </div>
      </div>
    </header>
  );
}
