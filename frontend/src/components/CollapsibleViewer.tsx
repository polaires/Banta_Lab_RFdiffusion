'use client';

import dynamic from 'next/dynamic';
import { useStore } from '@/lib/store';

// Dynamic import for Molstar viewer (SSR incompatible)
const ProteinViewer = dynamic(
  () => import('@/components/ProteinViewer').then((mod) => mod.ProteinViewer),
  {
    ssr: false,
    loading: () => (
      <div className="h-80 flex items-center justify-center bg-slate-900 rounded-xl">
        <div className="flex items-center gap-3 text-slate-400">
          <div className="w-5 h-5 border-2 border-slate-400 border-t-transparent rounded-full animate-spin" />
          <span className="text-sm font-medium">Loading 3D viewer...</span>
        </div>
      </div>
    ),
  }
);

export function CollapsibleViewer() {
  const { selectedPdb, viewerCollapsed, setViewerCollapsed } = useStore();

  return (
    <section className="bg-white rounded-2xl shadow-card overflow-hidden mt-8">
      {/* Header */}
      <button
        onClick={() => setViewerCollapsed(!viewerCollapsed)}
        className="w-full px-6 py-4 border-b border-slate-100 flex justify-between items-center bg-slate-50/50 hover:bg-slate-50 transition-colors"
      >
        <div className="flex items-center gap-3">
          <span className="material-symbols-outlined text-blue-600 text-xl">visibility</span>
          <h3 className="text-xs font-bold text-slate-700 uppercase tracking-widest">Structure Viewer</h3>
          {selectedPdb ? (
            <span className="text-[10px] font-semibold text-emerald-700 bg-emerald-50 px-2.5 py-1 rounded-full border border-emerald-200">
              Structure loaded
            </span>
          ) : (
            <span className="text-[10px] font-medium text-slate-500 bg-white px-2 py-0.5 rounded-full border border-slate-200 shadow-sm">
              No structure
            </span>
          )}
        </div>
        <span className="material-symbols-outlined text-slate-400 hover:text-slate-600 transition-colors">
          {viewerCollapsed ? 'expand_more' : 'expand_less'}
        </span>
      </button>

      {/* Viewer content */}
      {!viewerCollapsed && (
        <div className="p-6">
          {selectedPdb ? (
            <ProteinViewer pdbContent={selectedPdb} className="h-96 rounded-xl overflow-hidden" />
          ) : (
            <div className="h-64 flex flex-col items-center justify-center bg-slate-50 rounded-xl border-2 border-dashed border-slate-200">
              <span className="material-symbols-outlined text-4xl text-slate-300 mb-3">view_in_ar</span>
              <p className="text-sm font-medium text-slate-500">No structure to display</p>
              <p className="text-xs text-slate-400 mt-1">Run a design job to visualize structures</p>
            </div>
          )}
        </div>
      )}
    </section>
  );
}
