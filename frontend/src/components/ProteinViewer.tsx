'use client';

import { useEffect, useRef } from 'react';

interface ProteinViewerProps {
  pdbContent: string | null;
  className?: string;
}

export function ProteinViewer({ pdbContent, className = '' }: ProteinViewerProps) {
  const containerRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    if (!pdbContent || !containerRef.current) return;

    // For now, display PDB as text
    // TODO: Integrate Mol* or 3Dmol.js for 3D visualization
    const pre = document.createElement('pre');
    pre.className = 'text-xs font-mono overflow-auto h-full p-4 bg-gray-900 text-green-400';
    pre.textContent = pdbContent.slice(0, 5000) + (pdbContent.length > 5000 ? '\n...' : '');

    containerRef.current.innerHTML = '';
    containerRef.current.appendChild(pre);
  }, [pdbContent]);

  if (!pdbContent) {
    return (
      <div className={`flex items-center justify-center bg-gray-800 rounded-lg ${className}`}>
        <p className="text-gray-400">No structure to display</p>
      </div>
    );
  }

  return (
    <div
      ref={containerRef}
      className={`bg-gray-900 rounded-lg overflow-hidden ${className}`}
    />
  );
}
