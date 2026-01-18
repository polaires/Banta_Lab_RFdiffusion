'use client';

import { useState } from 'react';
import { SlidersHorizontal, ChevronDown } from 'lucide-react';

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
    <div className={`border border-border rounded-xl overflow-hidden ${className}`}>
      {/* Toggle header */}
      <button
        type="button"
        onClick={() => setExpanded(!expanded)}
        className="w-full px-4 py-3 flex items-center justify-between bg-muted hover:bg-muted/80 transition-colors"
      >
        <div className="flex items-center gap-2">
          <SlidersHorizontal className="h-5 w-5 text-muted-foreground" />
          <span className="text-sm font-medium text-foreground">{title}</span>
        </div>
        <ChevronDown
          className={`h-5 w-5 text-muted-foreground transition-transform duration-200 ${
            expanded ? 'rotate-180' : ''
          }`}
        />
      </button>

      {/* Expandable content */}
      <div
        className={`grid transition-all duration-200 ease-in-out ${
          expanded ? 'grid-rows-[1fr]' : 'grid-rows-[0fr]'
        }`}
      >
        <div className="overflow-hidden">
          <div className="px-4 py-4 space-y-4 border-t border-border">
            {children}
          </div>
        </div>
      </div>
    </div>
  );
}

export default AdvancedOptionsWrapper;
