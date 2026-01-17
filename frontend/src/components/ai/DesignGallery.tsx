'use client';

import { useState } from 'react';
import { LayoutGrid, AlertTriangle, Link, CheckCircle, type LucideIcon } from 'lucide-react';
import { type DesignResult } from '@/lib/demoData';

// Icon mapping for status icons
type StatusIconInfo = { IconComponent: LucideIcon; color: string; bg: string };

const STATUS_ICON_MAP = {
  warning: { IconComponent: AlertTriangle, color: 'text-red-500', bg: 'bg-red-100' },
  link: { IconComponent: Link, color: 'text-amber-500', bg: 'bg-amber-100' },
  check_circle: { IconComponent: CheckCircle, color: 'text-green-500', bg: 'bg-green-100' },
} as const;

// Re-export for convenience
export type { DesignResult };

interface DesignGalleryProps {
  designs: DesignResult[];
  selectedDesignId: string | null;
  onSelectDesign: (design: DesignResult) => void;
  designType: 'ligand_interface' | 'metal_binding' | 'standard';
}

export function DesignGallery({
  designs,
  selectedDesignId,
  onSelectDesign,
  designType,
}: DesignGalleryProps) {
  const [sortBy, setSortBy] = useState<'rank' | 'affinity' | 'contacts'>('rank');

  // Sort designs based on selected criteria
  const sortedDesigns = [...designs].sort((a, b) => {
    if (sortBy === 'affinity') {
      const affinityA = a.metrics.affinity ?? 0;
      const affinityB = b.metrics.affinity ?? 0;
      return affinityA - affinityB; // Lower (more negative) is better
    }
    if (sortBy === 'contacts') {
      const contactsA = a.metrics.total_contacts ?? (a.metrics.contacts_a ?? 0) + (a.metrics.contacts_b ?? 0);
      const contactsB = b.metrics.total_contacts ?? (b.metrics.contacts_a ?? 0) + (b.metrics.contacts_b ?? 0);
      return contactsB - contactsA; // Higher is better
    }
    return a.rank - b.rank;
  });

  const getAffinityColor = (affinity: number | undefined) => {
    if (affinity === undefined) return 'text-slate-500';
    if (affinity < -4) return 'text-green-600';
    if (affinity < -2) return 'text-amber-600';
    return 'text-red-600';
  };

  const getStatusIcon = (design: DesignResult): StatusIconInfo => {
    const hasClashes = design.metrics.has_clashes;
    const separable = design.metrics.separable;

    if (hasClashes) return STATUS_ICON_MAP.warning;
    if (separable === false) return STATUS_ICON_MAP.link;
    return STATUS_ICON_MAP.check_circle;
  };

  return (
    <div className="bg-white rounded-xl border border-slate-200 overflow-hidden">
      {/* Header */}
      <div className="bg-slate-50 px-4 py-3 border-b border-slate-200 flex items-center justify-between">
        <div className="flex items-center gap-2">
          <LayoutGrid className="h-5 w-5 text-purple-600" />
          <h4 className="font-semibold text-slate-900 text-sm">Design Gallery</h4>
          <span className="text-xs text-slate-500 bg-slate-200 px-2 py-0.5 rounded-full">
            {designs.length} designs
          </span>
        </div>

        {/* Sort Controls */}
        <div className="flex items-center gap-2">
          <span className="text-xs text-slate-500">Sort by:</span>
          <select
            value={sortBy}
            onChange={(e) => setSortBy(e.target.value as 'rank' | 'affinity' | 'contacts')}
            className="text-xs border border-slate-200 rounded px-2 py-1 bg-white"
          >
            <option value="rank">Rank</option>
            <option value="affinity">Affinity</option>
            <option value="contacts">Contacts</option>
          </select>
        </div>
      </div>

      {/* Design Grid */}
      <div className="p-3 grid grid-cols-1 gap-2 max-h-[400px] overflow-y-auto">
        {sortedDesigns.map((design) => {
          const isSelected = design.id === selectedDesignId;
          const status = getStatusIcon(design);

          return (
            <button
              key={design.id}
              onClick={() => onSelectDesign(design)}
              className={`w-full text-left p-3 rounded-lg border-2 transition-all ${
                isSelected
                  ? 'border-purple-500 bg-purple-50'
                  : 'border-slate-200 bg-white hover:border-purple-300 hover:bg-purple-50/50'
              }`}
            >
              <div className="flex items-start justify-between">
                {/* Left: Rank and Status */}
                <div className="flex items-center gap-3">
                  <div className={`w-8 h-8 rounded-lg flex items-center justify-center font-bold text-sm ${
                    isSelected ? 'bg-purple-500 text-white' : 'bg-slate-200 text-slate-600'
                  }`}>
                    #{design.rank}
                  </div>

                  {/* Status Badge */}
                  <div className={`p-1.5 rounded-lg ${status.bg}`}>
                    <status.IconComponent className={`h-4 w-4 ${status.color}`} />
                  </div>
                </div>

                {/* Right: Key Metrics */}
                <div className="flex items-center gap-4 text-xs">
                  {/* Affinity */}
                  {design.metrics.affinity !== undefined && (
                    <div className="text-right">
                      <div className="text-slate-400 text-[10px] uppercase">Affinity</div>
                      <div className={`font-bold ${getAffinityColor(design.metrics.affinity)}`}>
                        {design.metrics.affinity.toFixed(2)} kcal/mol
                      </div>
                    </div>
                  )}

                  {/* Contacts */}
                  <div className="text-right">
                    <div className="text-slate-400 text-[10px] uppercase">Contacts</div>
                    <div className="font-medium text-slate-700">
                      {designType === 'ligand_interface' ? (
                        <>
                          <span className="text-purple-600">{design.metrics.contacts_a ?? 0}</span>
                          <span className="text-slate-400"> / </span>
                          <span className="text-cyan-600">{design.metrics.contacts_b ?? 0}</span>
                        </>
                      ) : (
                        design.metrics.total_contacts ?? ((design.metrics.contacts_a ?? 0) + (design.metrics.contacts_b ?? 0))
                      )}
                    </div>
                  </div>

                  {/* Interface Area (if available) */}
                  {design.metrics.interface_area && (
                    <div className="text-right">
                      <div className="text-slate-400 text-[10px] uppercase">Interface</div>
                      <div className="font-medium text-slate-700">
                        {design.metrics.interface_area.toFixed(0)} Å²
                      </div>
                    </div>
                  )}
                </div>
              </div>

              {/* Expanded details for selected design */}
              {isSelected && (
                <div className="mt-3 pt-3 border-t border-purple-200 grid grid-cols-2 gap-3">
                  {/* Chain A Metrics */}
                  {design.chain_a_metrics && (
                    <div className="bg-purple-100/50 rounded-lg p-2">
                      <div className="text-[10px] text-purple-700 font-medium uppercase mb-1">Chain A</div>
                      <div className="space-y-1 text-xs">
                        <div className="flex justify-between">
                          <span className="text-slate-500">Contacts:</span>
                          <span className="font-medium">{design.chain_a_metrics.contacts}</span>
                        </div>
                        {design.chain_a_metrics.affinity !== undefined && (
                          <div className="flex justify-between">
                            <span className="text-slate-500">Affinity:</span>
                            <span className="font-medium">{design.chain_a_metrics.affinity.toFixed(2)}</span>
                          </div>
                        )}
                        {design.chain_a_metrics.exposed_atoms.length > 0 && (
                          <div className="flex justify-between">
                            <span className="text-slate-500">Exposed:</span>
                            <span className="font-medium text-[10px]">
                              {design.chain_a_metrics.exposed_atoms.slice(0, 3).join(', ')}
                              {design.chain_a_metrics.exposed_atoms.length > 3 && '...'}
                            </span>
                          </div>
                        )}
                      </div>
                    </div>
                  )}

                  {/* Chain B Metrics */}
                  {design.chain_b_metrics && (
                    <div className="bg-cyan-100/50 rounded-lg p-2">
                      <div className="text-[10px] text-cyan-700 font-medium uppercase mb-1">Chain B</div>
                      <div className="space-y-1 text-xs">
                        <div className="flex justify-between">
                          <span className="text-slate-500">Contacts:</span>
                          <span className="font-medium">{design.chain_b_metrics.contacts}</span>
                        </div>
                        {design.chain_b_metrics.affinity !== undefined && (
                          <div className="flex justify-between">
                            <span className="text-slate-500">Affinity:</span>
                            <span className="font-medium">{design.chain_b_metrics.affinity.toFixed(2)}</span>
                          </div>
                        )}
                        {design.chain_b_metrics.exposed_atoms.length > 0 && (
                          <div className="flex justify-between">
                            <span className="text-slate-500">Exposed:</span>
                            <span className="font-medium text-[10px]">
                              {design.chain_b_metrics.exposed_atoms.slice(0, 3).join(', ')}
                              {design.chain_b_metrics.exposed_atoms.length > 3 && '...'}
                            </span>
                          </div>
                        )}
                      </div>
                    </div>
                  )}
                </div>
              )}

              {/* Additional Flags */}
              <div className="mt-2 flex items-center gap-2">
                {design.metrics.separable !== undefined && (
                  <span className={`text-[10px] px-2 py-0.5 rounded-full ${
                    design.metrics.separable
                      ? 'bg-green-100 text-green-700'
                      : 'bg-red-100 text-red-700'
                  }`}>
                    {design.metrics.separable ? 'Separable' : 'Entangled'}
                  </span>
                )}
                {design.metrics.has_clashes && (
                  <span className="text-[10px] px-2 py-0.5 rounded-full bg-red-100 text-red-700">
                    Clashes
                  </span>
                )}
                {design.metrics.plddt !== undefined && (
                  <span className="text-[10px] px-2 py-0.5 rounded-full bg-blue-100 text-blue-700">
                    pLDDT: {design.metrics.plddt.toFixed(1)}
                  </span>
                )}
              </div>
            </button>
          );
        })}
      </div>

      {/* Footer with summary stats */}
      <div className="bg-slate-50 px-4 py-2 border-t border-slate-200 flex items-center justify-between text-xs text-slate-500">
        <span>
          {designs.filter(d => !d.metrics.has_clashes && d.metrics.separable !== false).length} of {designs.length} pass criteria
        </span>
        {selectedDesignId && (
          <span className="text-purple-600 font-medium">
            Design #{designs.find(d => d.id === selectedDesignId)?.rank} selected
          </span>
        )}
      </div>
    </div>
  );
}
