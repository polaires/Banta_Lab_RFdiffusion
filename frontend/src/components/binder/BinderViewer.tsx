'use client';

import { useState, useRef, useEffect, useCallback } from 'react';
import dynamic from 'next/dynamic';

// Types for Mol* - loaded dynamically
type PluginUIContext = any;

interface BinderViewerProps {
  pdbContent: string | null;
  targetChain?: string;
  binderChain?: string;
  interfaceResidues?: {
    target: number[];
    binder: number[];
  };
  highlightedResidues?: { chain: string; residue: number }[];
  className?: string;
  colorMode?: 'chain' | 'interface' | 'confidence';
  onResidueClick?: (chain: string, residue: number) => void;
}

// Dynamically import Molstar with SSR disabled
const BinderViewerInner = dynamic(
  () => import('./BinderViewerInner').then((mod) => mod.BinderViewerInner),
  {
    ssr: false,
    loading: () => (
      <div className="relative bg-slate-100 rounded-lg overflow-hidden" style={{ minHeight: '300px' }}>
        <div className="absolute inset-0 flex items-center justify-center">
          <div className="flex items-center gap-2 text-slate-600">
            <div className="w-5 h-5 border-2 border-slate-400 border-t-transparent rounded-full animate-spin" />
            Loading viewer...
          </div>
        </div>
      </div>
    ),
  }
);

export function BinderViewer(props: BinderViewerProps) {
  return <BinderViewerInner {...props} />;
}
