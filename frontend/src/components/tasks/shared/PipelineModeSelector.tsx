'use client';

import { cn } from '@/lib/utils';
import { Beaker, Layers, Factory } from 'lucide-react';
import type { PipelineMode } from '@/lib/store';

interface PipelineModeSelectorProps {
  mode: PipelineMode;
  onModeChange: (mode: PipelineMode) => void;
  disabled?: boolean;
}

const MODES: Array<{
  id: PipelineMode;
  name: string;
  description: string;
  icon: React.ComponentType<{ className?: string }>;
}> = [
  {
    id: 'single',
    name: 'Single Design',
    description: 'Generate a few designs with one configuration',
    icon: Beaker,
  },
  {
    id: 'sweep',
    name: 'Parameter Sweep',
    description: 'Explore 9 configurations (3 sizes x 3 CFG values)',
    icon: Layers,
  },
  {
    id: 'production',
    name: 'Production Run',
    description: 'Generate many designs with optimal configuration',
    icon: Factory,
  },
];

export function PipelineModeSelector({
  mode,
  onModeChange,
  disabled = false,
}: PipelineModeSelectorProps) {
  return (
    <div className="flex gap-1 p-1 bg-muted rounded-lg">
      {MODES.map((modeOption) => {
        const Icon = modeOption.icon;
        const isSelected = mode === modeOption.id;

        return (
          <button
            key={modeOption.id}
            type="button"
            onClick={() => onModeChange(modeOption.id)}
            disabled={disabled}
            className={cn(
              "flex-1 flex items-center gap-2 px-3 py-2 rounded-md text-sm font-medium transition-all",
              isSelected
                ? "bg-background text-foreground shadow-sm"
                : "text-muted-foreground hover:text-foreground hover:bg-background/50",
              disabled && "opacity-50 cursor-not-allowed"
            )}
            title={modeOption.description}
          >
            <Icon className={cn(
              "w-4 h-4",
              isSelected ? "text-primary" : "text-muted-foreground"
            )} />
            <span className="hidden sm:inline">{modeOption.name}</span>
          </button>
        );
      })}
    </div>
  );
}

export function PipelineModeDescription({ mode }: { mode: PipelineMode }) {
  const modeInfo = MODES.find((m) => m.id === mode);
  if (!modeInfo) return null;

  return (
    <p className="text-xs text-muted-foreground mt-2">
      {modeInfo.description}
    </p>
  );
}
