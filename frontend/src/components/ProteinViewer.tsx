'use client';

import dynamic from 'next/dynamic';
import type { MetalCoordination } from '@/lib/metalAnalysis';
import type { LigandData } from '@/lib/ligandAnalysis';

// Import Molstar CSS
import 'molstar/lib/mol-plugin-ui/skin/light.scss';

interface ProteinViewerProps {
  pdbContent: string | null;
  className?: string;
  // Focus props
  focusedMetalIndex?: number | null;
  focusedLigandIndex?: number | null;
  metalCoordination?: MetalCoordination[] | null;
  ligandData?: { ligandDetails: LigandData[] } | null;
}

// Dynamically import the viewer with SSR disabled to avoid Turbopack bundling issues
const ProteinViewerClient = dynamic(
  () => import('./ProteinViewerClient').then((mod) => mod.ProteinViewerClient),
  {
    ssr: false,
    loading: () => (
      <div className="relative bg-gray-100 rounded-lg overflow-hidden" style={{ minHeight: '300px' }}>
        <div className="absolute inset-0 flex items-center justify-center">
          <div className="flex items-center gap-2 text-gray-600">
            <div className="w-5 h-5 border-2 border-gray-400 border-t-transparent rounded-full animate-spin" />
            Loading viewer...
          </div>
        </div>
      </div>
    ),
  }
);

export function ProteinViewer({
  pdbContent,
  className = '',
  focusedMetalIndex,
  focusedLigandIndex,
  metalCoordination,
  ligandData,
}: ProteinViewerProps) {
  return (
    <ProteinViewerClient
      pdbContent={pdbContent}
      className={className}
      focusedMetalIndex={focusedMetalIndex}
      focusedLigandIndex={focusedLigandIndex}
      metalCoordination={metalCoordination}
      ligandData={ligandData}
    />
  );
}
