'use client';

import { ClipboardList, Gem, SlidersHorizontal, Network, Copy, Star, Pencil, Play } from 'lucide-react';
import type { UserPreferences } from '@/lib/api';

interface PreferenceSummaryCardProps {
  preferences: UserPreferences;
  onConfirm: () => void;
  onEdit: () => void;
  isRunning?: boolean;
}

export function PreferenceSummaryCard({
  preferences,
  onConfirm,
  onEdit,
  isRunning = false,
}: PreferenceSummaryCardProps) {
  return (
    <div className="bg-card rounded-xl p-5 border border-border shadow-sm">
      <h4 className="font-semibold text-foreground mb-4 flex items-center gap-2">
        <ClipboardList className="h-5 w-5 text-primary" />
        Your Design Preferences
      </h4>

      <div className="space-y-3 text-sm">
        <div className="flex justify-between items-center py-2 border-b border-border">
          <span className="text-muted-foreground flex items-center gap-2">
            <Gem className="h-4 w-4 text-muted-foreground" />
            Target Metal
          </span>
          <span className="font-medium text-foreground">{preferences.targetMetalLabel}</span>
        </div>
        <div className="flex justify-between items-center py-2 border-b border-border">
          <span className="text-muted-foreground flex items-center gap-2">
            <SlidersHorizontal className="h-4 w-4 text-muted-foreground" />
            Design Approach
          </span>
          <span className="font-medium text-foreground">{preferences.aggressivenessLabel}</span>
        </div>
        <div className="flex justify-between items-center py-2 border-b border-border">
          <span className="text-muted-foreground flex items-center gap-2">
            <Network className="h-4 w-4 text-muted-foreground" />
            Coordination
          </span>
          <span className="font-medium text-foreground">{preferences.coordinationLabel}</span>
        </div>
        <div className="flex justify-between items-center py-2 border-b border-border">
          <span className="text-muted-foreground flex items-center gap-2">
            <Copy className="h-4 w-4 text-muted-foreground" />
            Variants
          </span>
          <span className="font-medium text-foreground">{preferences.numDesigns} designs</span>
        </div>
        <div className="flex justify-between items-center py-2">
          <span className="text-muted-foreground flex items-center gap-2">
            <Star className="h-4 w-4 text-muted-foreground" />
            Priority
          </span>
          <span className="font-medium text-foreground">{preferences.priorityLabel}</span>
        </div>
      </div>

      <div className="flex gap-3 mt-6">
        <button
          onClick={onEdit}
          disabled={isRunning}
          className={`flex-1 py-2.5 text-primary border border-border rounded-lg font-medium transition-all ${
            isRunning
              ? 'opacity-50 cursor-not-allowed'
              : 'hover:bg-muted hover:border-primary/50'
          }`}
        >
          <span className="flex items-center justify-center gap-2">
            <Pencil className="h-4 w-4" />
            Edit
          </span>
        </button>
        <button
          onClick={onConfirm}
          disabled={isRunning}
          className={`flex-1 py-2.5 rounded-lg font-medium transition-all flex items-center justify-center gap-2 ${
            isRunning
              ? 'bg-primary/70 text-primary-foreground cursor-not-allowed'
              : 'bg-primary text-primary-foreground hover:bg-primary/90 shadow-sm'
          }`}
        >
          {isRunning ? (
            <>
              <svg className="animate-spin h-4 w-4" viewBox="0 0 24 24">
                <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" fill="none" />
                <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z" />
              </svg>
              Running...
            </>
          ) : (
            <>
              <Play className="h-4 w-4" />
              Run Design
            </>
          )}
        </button>
      </div>
    </div>
  );
}
