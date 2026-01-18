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
      <div className="border-b border-border pb-6">
        <div className="flex items-center gap-3 mb-2">
          <span className="bg-muted text-muted-foreground border border-border text-[11px] font-bold px-2.5 py-0.5 rounded-full uppercase tracking-wide">
            Step 0
          </span>
          <h2 className="text-xl font-bold text-foreground">Choose Design Task</h2>
        </div>
        <p className="text-muted-foreground text-sm leading-relaxed max-w-3xl pl-1">
          Select the type of protein design you want to perform. Each task has tailored options and presets.
        </p>
      </div>

      {/* Task Selector */}
      <TaskSelector onTaskSelect={handleTaskSelect} />
    </div>
  );
}
