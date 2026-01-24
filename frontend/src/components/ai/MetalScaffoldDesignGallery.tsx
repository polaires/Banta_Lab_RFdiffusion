'use client';

import { useState, useMemo } from 'react';
import {
  LayoutGrid,
  Download,
  Copy,
  CheckCircle,
  Trash2,
  Filter,
} from 'lucide-react';
import type { MetalScaffoldDesignResult } from './index';

// Tier color configuration matching shadcn/ui patterns
const TIER_COLORS = {
  S: { bg: 'bg-green-100', text: 'text-green-700', border: 'border-green-300', badge: 'bg-green-500' },
  A: { bg: 'bg-blue-100', text: 'text-blue-700', border: 'border-blue-300', badge: 'bg-blue-500' },
  B: { bg: 'bg-yellow-100', text: 'text-yellow-700', border: 'border-yellow-300', badge: 'bg-yellow-500' },
  C: { bg: 'bg-orange-100', text: 'text-orange-700', border: 'border-orange-300', badge: 'bg-orange-500' },
  F: { bg: 'bg-red-100', text: 'text-red-700', border: 'border-red-300', badge: 'bg-red-500' },
} as const;

// Tier sort order (S is best)
const TIER_ORDER = { S: 0, A: 1, B: 2, C: 3, F: 4 };

type FilterMode = 'all' | 'passed' | 'trashed';
type SortMode = 'tier' | 'plddt' | 'ptm' | 'config';

interface MetalScaffoldDesignGalleryProps {
  designs: MetalScaffoldDesignResult[];
  selectedDesignId: string | null;
  onSelectDesign: (design: MetalScaffoldDesignResult) => void;
  onExportFasta: (designs: MetalScaffoldDesignResult[]) => void;
}

export function MetalScaffoldDesignGallery({
  designs,
  selectedDesignId,
  onSelectDesign,
  onExportFasta,
}: MetalScaffoldDesignGalleryProps) {
  const [filterMode, setFilterMode] = useState<FilterMode>('passed');
  const [sortBy, setSortBy] = useState<SortMode>('tier');
  const [copiedId, setCopiedId] = useState<string | null>(null);

  // Count statistics
  const stats = useMemo(() => {
    const passed = designs.filter((d) => d.status === 'passed').length;
    const trashed = designs.filter((d) => d.status === 'failed').length;
    return { total: designs.length, passed, trashed };
  }, [designs]);

  // Filter designs based on mode
  const filteredDesigns = useMemo(() => {
    switch (filterMode) {
      case 'passed':
        return designs.filter((d) => d.status === 'passed');
      case 'trashed':
        return designs.filter((d) => d.status === 'failed');
      default:
        return designs;
    }
  }, [designs, filterMode]);

  // Sort designs
  const sortedDesigns = useMemo(() => {
    return [...filteredDesigns].sort((a, b) => {
      switch (sortBy) {
        case 'tier':
          return TIER_ORDER[a.tier] - TIER_ORDER[b.tier];
        case 'plddt':
          return b.plddt - a.plddt; // Higher is better
        case 'ptm':
          return b.ptm - a.ptm; // Higher is better
        case 'config':
          return a.configName.localeCompare(b.configName);
        default:
          return 0;
      }
    });
  }, [filteredDesigns, sortBy]);

  // Copy sequence to clipboard
  const handleCopySequence = async (design: MetalScaffoldDesignResult) => {
    try {
      await navigator.clipboard.writeText(design.sequence);
      setCopiedId(design.id);
      setTimeout(() => setCopiedId(null), 2000);
    } catch (err) {
      console.error('Failed to copy sequence:', err);
    }
  };

  // Truncate sequence for display
  const truncateSequence = (seq: string, maxLength: number = 20) => {
    if (seq.length <= maxLength) return seq;
    return `${seq.substring(0, maxLength)}...`;
  };

  // Get designs for export (only passed designs)
  const passingDesigns = designs.filter((d) => d.status === 'passed');

  return (
    <div className="bg-card rounded-xl border border-border overflow-hidden">
      {/* Header */}
      <div className="bg-muted/50 px-4 py-3 border-b border-border">
        <div className="flex items-center justify-between mb-3">
          <div className="flex items-center gap-2">
            <LayoutGrid className="h-5 w-5 text-primary" />
            <h4 className="font-semibold text-foreground text-sm">Individual Design Results</h4>
            <span className="text-xs text-muted-foreground bg-muted px-2 py-0.5 rounded-full">
              {stats.passed}/{stats.total} passed
            </span>
          </div>

          {/* Sort Controls */}
          <div className="flex items-center gap-2">
            <span className="text-xs text-muted-foreground">Sort by:</span>
            <select
              value={sortBy}
              onChange={(e) => setSortBy(e.target.value as SortMode)}
              className="text-xs border border-border rounded px-2 py-1 bg-card"
            >
              <option value="tier">Tier</option>
              <option value="plddt">pLDDT</option>
              <option value="ptm">pTM</option>
              <option value="config">Config</option>
            </select>
          </div>
        </div>

        {/* Filter Tabs */}
        <div className="flex gap-1">
          <button
            onClick={() => setFilterMode('all')}
            className={`px-3 py-1.5 rounded-md text-xs font-medium transition-colors ${
              filterMode === 'all'
                ? 'bg-primary text-primary-foreground'
                : 'bg-muted text-muted-foreground hover:bg-muted/80'
            }`}
          >
            <Filter className="h-3 w-3 inline-block mr-1" />
            All ({stats.total})
          </button>
          <button
            onClick={() => setFilterMode('passed')}
            className={`px-3 py-1.5 rounded-md text-xs font-medium transition-colors ${
              filterMode === 'passed'
                ? 'bg-green-500 text-white'
                : 'bg-green-100 text-green-700 hover:bg-green-200'
            }`}
          >
            <CheckCircle className="h-3 w-3 inline-block mr-1" />
            Passed ({stats.passed})
          </button>
          <button
            onClick={() => setFilterMode('trashed')}
            className={`px-3 py-1.5 rounded-md text-xs font-medium transition-colors ${
              filterMode === 'trashed'
                ? 'bg-red-500 text-white'
                : 'bg-red-100 text-red-700 hover:bg-red-200'
            }`}
          >
            <Trash2 className="h-3 w-3 inline-block mr-1" />
            Trashed ({stats.trashed})
          </button>
        </div>
      </div>

      {/* Design Grid */}
      <div className="p-3 grid grid-cols-1 md:grid-cols-2 gap-2 max-h-[450px] overflow-y-auto">
        {sortedDesigns.length === 0 ? (
          <div className="col-span-full text-center py-8 text-muted-foreground">
            No designs match the current filter.
          </div>
        ) : (
          sortedDesigns.map((design) => {
            const isSelected = design.id === selectedDesignId;
            const isTrashed = design.status === 'failed';
            const tierColors = TIER_COLORS[design.tier];

            return (
              <button
                key={design.id}
                onClick={() => onSelectDesign(design)}
                className={`w-full text-left p-3 rounded-lg border-2 transition-all ${
                  isSelected
                    ? 'border-primary bg-primary/5 shadow-sm'
                    : isTrashed
                    ? 'border-border bg-muted/50 opacity-60 hover:opacity-80'
                    : 'border-border bg-card hover:border-primary/50 hover:bg-muted/30'
                }`}
              >
                <div className="flex items-start justify-between gap-2">
                  {/* Left: Tier Badge + Name */}
                  <div className="flex items-center gap-2 min-w-0">
                    {/* Tier Badge */}
                    <div
                      className={`w-8 h-8 rounded-lg flex items-center justify-center font-bold text-sm text-white ${tierColors.badge}`}
                    >
                      {design.tier}
                    </div>
                    <div className="min-w-0">
                      <div
                        className={`font-medium text-sm truncate ${
                          isTrashed ? 'line-through text-muted-foreground' : 'text-foreground'
                        }`}
                      >
                        {design.name}
                      </div>
                      <div className="text-xs text-muted-foreground truncate">
                        {design.configName}
                      </div>
                    </div>
                  </div>

                  {/* Right: Metrics */}
                  <div className="flex items-center gap-3 text-xs shrink-0">
                    <div className="text-right">
                      <div className="text-muted-foreground text-[10px] uppercase">pLDDT</div>
                      <div
                        className={`font-bold ${
                          design.plddt >= 0.8
                            ? 'text-green-600'
                            : design.plddt >= 0.7
                            ? 'text-yellow-600'
                            : 'text-red-600'
                        }`}
                      >
                        {(design.plddt * 100).toFixed(0)}
                      </div>
                    </div>
                    <div className="text-right">
                      <div className="text-muted-foreground text-[10px] uppercase">pTM</div>
                      <div
                        className={`font-bold ${
                          design.ptm >= 0.7
                            ? 'text-green-600'
                            : design.ptm >= 0.6
                            ? 'text-yellow-600'
                            : 'text-red-600'
                        }`}
                      >
                        {design.ptm.toFixed(2)}
                      </div>
                    </div>
                    <div className="text-right">
                      <div className="text-muted-foreground text-[10px] uppercase">pAE</div>
                      <div
                        className={`font-bold ${
                          design.pae <= 5 ? 'text-green-600' : design.pae <= 8 ? 'text-yellow-600' : 'text-red-600'
                        }`}
                      >
                        {design.pae.toFixed(1)}
                      </div>
                    </div>
                  </div>
                </div>

                {/* Sequence preview */}
                <div className="mt-2 flex items-center gap-2">
                  <code
                    className={`text-[10px] px-2 py-0.5 rounded flex-1 truncate ${
                      isTrashed ? 'bg-muted text-muted-foreground' : 'bg-muted text-foreground'
                    }`}
                  >
                    {truncateSequence(design.sequence, 30)}
                  </code>
                  <span
                    role="button"
                    tabIndex={0}
                    onClick={(e) => {
                      e.stopPropagation();
                      handleCopySequence(design);
                    }}
                    onKeyDown={(e) => {
                      if (e.key === 'Enter' || e.key === ' ') {
                        e.stopPropagation();
                        handleCopySequence(design);
                      }
                    }}
                    className={`p-1 rounded hover:bg-muted/80 transition-colors cursor-pointer ${
                      copiedId === design.id ? 'text-green-500' : 'text-muted-foreground'
                    }`}
                    title="Copy sequence"
                  >
                    {copiedId === design.id ? (
                      <CheckCircle className="h-3.5 w-3.5" />
                    ) : (
                      <Copy className="h-3.5 w-3.5" />
                    )}
                  </span>
                </div>

                {/* Expanded details for selected design */}
                {isSelected && (
                  <div className="mt-3 pt-3 border-t border-border">
                    <div className="text-xs text-muted-foreground mb-1">Full Sequence:</div>
                    <code className="text-[10px] text-foreground break-all block bg-muted p-2 rounded">
                      {design.sequence}
                    </code>
                    <div className="text-xs text-muted-foreground mt-2">
                      Length: {design.sequence.length} aa
                    </div>
                  </div>
                )}
              </button>
            );
          })
        )}
      </div>

      {/* Footer with export */}
      <div className="bg-muted/50 px-4 py-3 border-t border-border flex items-center justify-between">
        <span className="text-xs text-muted-foreground">
          {filterMode === 'passed' && sortedDesigns.length > 0
            ? `Showing ${sortedDesigns.length} passing designs`
            : filterMode === 'trashed'
            ? `Showing ${sortedDesigns.length} trashed designs`
            : `Showing ${sortedDesigns.length} of ${stats.total} designs`}
        </span>

        <button
          onClick={() => onExportFasta(passingDesigns)}
          disabled={passingDesigns.length === 0}
          className={`flex items-center gap-2 px-3 py-1.5 rounded-md text-xs font-medium transition-colors ${
            passingDesigns.length > 0
              ? 'bg-primary text-primary-foreground hover:bg-primary/90'
              : 'bg-muted text-muted-foreground cursor-not-allowed'
          }`}
        >
          <Download className="h-3.5 w-3.5" />
          Export Passed FASTA ({passingDesigns.length})
        </button>
      </div>
    </div>
  );
}
