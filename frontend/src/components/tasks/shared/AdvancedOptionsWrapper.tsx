'use client';

import { useState } from 'react';

interface AdvancedOptionsWrapperProps {
  children: React.ReactNode;
  defaultExpanded?: boolean;
  title?: string;
  className?: string;
}

export function AdvancedOptionsWrapper({
  children,
  defaultExpanded = false,
  title = 'Advanced Options',
  className = '',
}: AdvancedOptionsWrapperProps) {
  const [expanded, setExpanded] = useState(defaultExpanded);

  return (
    <div className={`border border-slate-200 rounded-xl overflow-hidden ${className}`}>
      {/* Toggle header */}
      <button
        type="button"
        onClick={() => setExpanded(!expanded)}
        className="w-full px-4 py-3 flex items-center justify-between bg-slate-50 hover:bg-slate-100 transition-colors"
      >
        <div className="flex items-center gap-2">
          <span className="material-symbols-outlined text-slate-500 text-lg">
            tune
          </span>
          <span className="text-sm font-medium text-slate-700">{title}</span>
        </div>
        <span
          className={`material-symbols-outlined text-slate-400 transition-transform duration-200 ${
            expanded ? 'rotate-180' : ''
          }`}
        >
          expand_more
        </span>
      </button>

      {/* Expandable content */}
      <div
        className={`grid transition-all duration-200 ease-in-out ${
          expanded ? 'grid-rows-[1fr]' : 'grid-rows-[0fr]'
        }`}
      >
        <div className="overflow-hidden">
          <div className="px-4 py-4 space-y-4 border-t border-slate-100">
            {children}
          </div>
        </div>
      </div>
    </div>
  );
}

export default AdvancedOptionsWrapper;
