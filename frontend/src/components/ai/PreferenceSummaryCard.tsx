'use client';

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
    <div className="bg-white rounded-xl p-5 border border-slate-200 shadow-sm">
      <h4 className="font-semibold text-slate-900 mb-4 flex items-center gap-2">
        <span className="material-symbols-outlined text-violet-600">checklist</span>
        Your Design Preferences
      </h4>

      <div className="space-y-3 text-sm">
        <div className="flex justify-between items-center py-2 border-b border-slate-100">
          <span className="text-slate-600 flex items-center gap-2">
            <span className="material-symbols-outlined text-sm text-slate-400">diamond</span>
            Target Metal
          </span>
          <span className="font-medium text-slate-900">{preferences.targetMetalLabel}</span>
        </div>
        <div className="flex justify-between items-center py-2 border-b border-slate-100">
          <span className="text-slate-600 flex items-center gap-2">
            <span className="material-symbols-outlined text-sm text-slate-400">tune</span>
            Design Approach
          </span>
          <span className="font-medium text-slate-900">{preferences.aggressivenessLabel}</span>
        </div>
        <div className="flex justify-between items-center py-2 border-b border-slate-100">
          <span className="text-slate-600 flex items-center gap-2">
            <span className="material-symbols-outlined text-sm text-slate-400">hub</span>
            Coordination
          </span>
          <span className="font-medium text-slate-900">{preferences.coordinationLabel}</span>
        </div>
        <div className="flex justify-between items-center py-2 border-b border-slate-100">
          <span className="text-slate-600 flex items-center gap-2">
            <span className="material-symbols-outlined text-sm text-slate-400">content_copy</span>
            Variants
          </span>
          <span className="font-medium text-slate-900">{preferences.numDesigns} designs</span>
        </div>
        <div className="flex justify-between items-center py-2">
          <span className="text-slate-600 flex items-center gap-2">
            <span className="material-symbols-outlined text-sm text-slate-400">star</span>
            Priority
          </span>
          <span className="font-medium text-slate-900">{preferences.priorityLabel}</span>
        </div>
      </div>

      <div className="flex gap-3 mt-6">
        <button
          onClick={onEdit}
          disabled={isRunning}
          className={`flex-1 py-2.5 text-violet-600 border border-violet-200 rounded-lg font-medium transition-all ${
            isRunning
              ? 'opacity-50 cursor-not-allowed'
              : 'hover:bg-violet-50 hover:border-violet-300'
          }`}
        >
          <span className="flex items-center justify-center gap-2">
            <span className="material-symbols-outlined text-sm">edit</span>
            Edit
          </span>
        </button>
        <button
          onClick={onConfirm}
          disabled={isRunning}
          className={`flex-1 py-2.5 rounded-lg font-medium transition-all flex items-center justify-center gap-2 ${
            isRunning
              ? 'bg-violet-400 text-white cursor-not-allowed'
              : 'bg-gradient-to-r from-violet-600 to-purple-600 text-white hover:from-violet-700 hover:to-purple-700 shadow-sm'
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
              <span className="material-symbols-outlined text-sm">play_arrow</span>
              Run Design
            </>
          )}
        </button>
      </div>
    </div>
  );
}
