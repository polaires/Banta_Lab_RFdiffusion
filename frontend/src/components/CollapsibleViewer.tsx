'use client';

import dynamic from 'next/dynamic';
import { useStore } from '@/lib/store';
import { ChevronDown, ChevronUp, Eye } from 'lucide-react';

// Dynamic import for Molstar viewer (SSR incompatible)
const ProteinViewer = dynamic(
  () => import('@/components/ProteinViewer').then((mod) => mod.ProteinViewer),
  {
    ssr: false,
    loading: () => (
      <div className="h-80 flex items-center justify-center bg-gray-800 rounded-lg">
        <div className="flex items-center gap-2 text-gray-400">
          <div className="w-5 h-5 border-2 border-gray-400 border-t-transparent rounded-full animate-spin" />
          Loading 3D viewer...
        </div>
      </div>
    ),
  }
);

export function CollapsibleViewer() {
  const { selectedPdb, viewerCollapsed, setViewerCollapsed } = useStore();

  return (
    <div className="bg-white rounded-lg shadow-sm border border-gray-200 overflow-hidden mt-6">
      {/* Header */}
      <button
        onClick={() => setViewerCollapsed(!viewerCollapsed)}
        className="w-full flex items-center justify-between px-4 py-3 text-left hover:bg-gray-50 transition"
      >
        <div className="flex items-center gap-2">
          <Eye className="w-5 h-5 text-blue-600" />
          <h3 className="font-semibold text-gray-900">Structure Viewer</h3>
          {selectedPdb && (
            <span className="text-xs text-green-700 bg-green-100 px-2 py-0.5 rounded-full">
              Structure loaded
            </span>
          )}
          {!selectedPdb && (
            <span className="text-xs text-gray-500">No structure</span>
          )}
        </div>
        {viewerCollapsed ? (
          <ChevronDown className="w-5 h-5 text-gray-500" />
        ) : (
          <ChevronUp className="w-5 h-5 text-gray-500" />
        )}
      </button>

      {/* Viewer content */}
      {!viewerCollapsed && (
        <div className="px-4 pb-4">
          <ProteinViewer pdbContent={selectedPdb} className="h-80" />
        </div>
      )}
    </div>
  );
}
