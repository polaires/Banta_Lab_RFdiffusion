'use client';

import { useStore, TabId } from '@/lib/store';
import { UserMenu } from './auth/UserMenu';
import { FlaskConical, Sparkles, History, Server, type LucideIcon } from 'lucide-react';

const workflowSteps: { id: 'ai' | 'task' | 'rfd3' | 'mpnn' | 'rf3'; label: string; step: number; Icon?: LucideIcon; gradient?: boolean }[] = [
  { id: 'ai', label: 'AI Assistant', step: 0, Icon: Sparkles, gradient: true },
  { id: 'task', label: 'Design Task', step: 1 },
  { id: 'rfd3', label: 'RFdiffusion3', step: 2 },
  { id: 'mpnn', label: 'MPNN', step: 3 },
  { id: 'rf3', label: 'Validate', step: 4 },
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

  const aiCaseStudy = useStore((state) => state.aiCaseStudy);

  const getStepStatus = (stepId: string) => {
    if (activeTab === stepId) return 'active';
    if (stepId === 'ai' && aiCaseStudy?.currentStep === 'complete') return 'completed';
    if (stepId === 'task' && selectedDesignTask) return 'completed';
    if (stepId === 'rfd3' && latestRfd3Design) return 'completed';
    if (stepId === 'mpnn' && lastCompletedJobType === 'mpnn') return 'completed';
    if (stepId === 'rf3' && lastCompletedJobType === 'rf3') return 'completed';
    return 'upcoming';
  };

  return (
    <header className="bg-white/80 backdrop-blur-md border-b border-border sticky top-0 z-50">
      <div className="max-w-6xl mx-auto px-6 h-16 flex items-center justify-between">
        {/* Logo & Brand */}
        <div className="flex items-center gap-3">
          <div className="text-primary bg-primary/10 rounded-lg p-1.5 flex items-center justify-center">
            <FlaskConical className="w-6 h-6" />
          </div>
          <div>
            <h1 className="font-bold text-sm text-foreground leading-tight">Foundry Protein Design</h1>
            <p className="text-[10px] text-muted-foreground font-medium uppercase tracking-wide">IPD Design Pipeline</p>
          </div>
        </div>

        {/* Workflow Stepper */}
        <nav className="hidden md:flex items-center space-x-1">
          {workflowSteps.map((step, index) => {
            const status = getStepStatus(step.id);
            const isActive = status === 'active';

            return (
              <div key={step.id} className="flex items-center">
                {index > 0 && <div className="w-8 h-px bg-muted" />}

                <button
                  onClick={() => setActiveTab(step.id)}
                  className={`flex items-center gap-2 px-3 py-1.5 rounded-full transition-all ${
                    isActive
                      ? 'bg-muted border border-border shadow-sm'
                      : 'opacity-50 grayscale hover:opacity-75'
                  }`}
                >
                  {step.Icon ? (
                    <span className={`w-5 h-5 flex items-center justify-center rounded-full ${
                      isActive
                        ? 'bg-primary text-primary-foreground'
                        : 'bg-muted text-muted-foreground'
                    }`}>
                      <step.Icon className="w-3 h-3" />
                    </span>
                  ) : (
                    <span className={`w-5 h-5 flex items-center justify-center rounded-full text-[10px] font-bold ${
                      isActive
                        ? 'bg-primary text-primary-foreground'
                        : 'bg-muted text-muted-foreground'
                    }`}>
                      {step.step}
                    </span>
                  )}
                  <span className={`text-xs font-semibold ${
                    isActive
                      ? 'text-foreground'
                      : 'text-muted-foreground'
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
                    ? 'bg-primary text-primary-foreground'
                    : 'bg-muted text-muted-foreground'
                }`}
              >
                {step.Icon ? (
                  <step.Icon className="w-3.5 h-3.5" />
                ) : (
                  step.step
                )}
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
                ? 'text-foreground bg-muted'
                : 'text-muted-foreground hover:text-foreground hover:bg-muted'
            }`}
            title="Job History"
          >
            <History className="w-5 h-5" />
            {pendingJobs > 0 && (
              <span className="absolute -top-0.5 -right-0.5 w-4 h-4 bg-primary text-primary-foreground text-[9px] font-bold rounded-full flex items-center justify-center border border-white">
                {pendingJobs}
              </span>
            )}
          </button>

          <div className="h-4 w-px bg-muted mx-1" />

          {/* Connection status */}
          <button
            onClick={() => setConnectionModalOpen(true)}
            className="p-2 rounded-full text-muted-foreground hover:text-foreground hover:bg-muted transition-all relative"
            title={isConnected ? 'Connected to backend' : 'Connect to backend'}
          >
            <Server className="w-5 h-5" />
            <span className={`absolute top-2 right-2 w-2 h-2 rounded-full border border-white ${
              isConnected ? 'bg-emerald-500' : 'bg-red-500'
            }`} />
          </button>

          <div className="h-4 w-px bg-muted mx-1" />

          {/* User menu */}
          <UserMenu />
        </div>
      </div>
    </header>
  );
}
