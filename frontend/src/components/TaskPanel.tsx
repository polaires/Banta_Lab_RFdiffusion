'use client';

import { useStore } from '@/lib/store';
import { TaskSelector, DesignTask } from './TaskSelector';

export function TaskPanel() {
  const { setActiveTab, setSelectedDesignTask } = useStore();

  const handleTaskSelect = (task: DesignTask) => {
    setSelectedDesignTask(task);
    setActiveTab('rfd3');
  };

  return (
    <div className="p-8 space-y-8">
      {/* Header */}
      <div className="border-b border-slate-100 pb-6">
        <div className="flex items-center gap-3 mb-2">
          <span className="bg-slate-100 text-slate-600 border border-slate-200 text-[11px] font-bold px-2.5 py-0.5 rounded-full uppercase tracking-wide">
            Step 0
          </span>
          <h2 className="text-xl font-bold text-slate-900">Choose Design Task</h2>
        </div>
        <p className="text-slate-500 text-sm leading-relaxed max-w-3xl pl-1">
          Select the type of protein design you want to perform. Each task has tailored options and presets.
        </p>
      </div>

      {/* Task Selector */}
      <TaskSelector onTaskSelect={handleTaskSelect} />
    </div>
  );
}
