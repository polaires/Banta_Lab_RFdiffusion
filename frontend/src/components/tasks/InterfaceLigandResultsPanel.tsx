'use client';

import { useState, useMemo } from 'react';
import {
  AlertCircle, Loader2, AlertTriangle, Network, SlidersHorizontal,
  GitBranch, Link, Type, ShieldCheck, Layers, View, Download
} from 'lucide-react';
import { SuccessCriteriaCard, type FilterMode, type DesignMetrics } from '../binder/SuccessCriteriaCard';

// Design result from the backend
export interface InterfaceLigandDesign {
  id: string;
  rank: number;
  pdb_content: string;
  sequence_a?: string;
  sequence_b?: string;
  metrics: DesignMetrics;
  approach?: string;
  status: 'passed' | 'failed' | 'pending';
  pipeline_stage?: 'backbone' | 'ligandmpnn' | 'validation' | 'complete';
}

// Job result structure
export interface InterfaceLigandJobResult {
  status: 'completed' | 'failed' | 'running';
  approach: string;
  designs: InterfaceLigandDesign[];
  best_design_idx?: number;
  error?: string;
}

interface InterfaceLigandResultsPanelProps {
  result: InterfaceLigandJobResult;
  onSelectDesign: (design: InterfaceLigandDesign) => void;
  onRunLigandMPNN?: (designs: InterfaceLigandDesign[]) => Promise<void>;
  onRunValidation?: (designs: InterfaceLigandDesign[]) => Promise<void>;
  onDownloadPdb?: (design: InterfaceLigandDesign) => void;
  selectedDesignId?: string;
  isProcessing?: boolean;
}

// Approach display names - simplified without colors
const APPROACH_LABELS: Record<string, { name: string; Icon: React.ComponentType<{ className?: string }> }> = {
  joint: { name: 'Joint Heterodimer', Icon: Network },
  asymmetric_rasa: { name: 'Asymmetric RASA', Icon: SlidersHorizontal },
  induced: { name: 'Induced Dimerization', Icon: GitBranch },
  asymmetric: { name: 'Single Chain', Icon: Link },
};

export function InterfaceLigandResultsPanel({
  result,
  onSelectDesign,
  onRunLigandMPNN,
  onRunValidation,
  onDownloadPdb,
  selectedDesignId,
  isProcessing = false,
}: InterfaceLigandResultsPanelProps) {
  const [sortBy, setSortBy] = useState<'rank' | 'affinity' | 'contacts' | 'identity'>('rank');
  const [filterStatus, setFilterStatus] = useState<'all' | 'passed' | 'failed'>('all');
  const [filterMode, setFilterMode] = useState<FilterMode>('heterodimer');
  const [selectedDesigns, setSelectedDesigns] = useState<Set<string>>(new Set());

  // Filter and sort designs
  const filteredDesigns = useMemo(() => {
    let designs = [...result.designs];

    // Filter by status
    if (filterStatus !== 'all') {
      designs = designs.filter(d => d.status === filterStatus);
    }

    // Sort
    designs.sort((a, b) => {
      switch (sortBy) {
        case 'affinity':
          return (a.metrics.affinity ?? 0) - (b.metrics.affinity ?? 0);
        case 'contacts':
          const contactsA = (a.metrics.contacts_a ?? 0) + (a.metrics.contacts_b ?? 0);
          const contactsB = (b.metrics.contacts_a ?? 0) + (b.metrics.contacts_b ?? 0);
          return contactsB - contactsA;
        case 'identity':
          return (a.metrics.sequence_identity ?? 100) - (b.metrics.sequence_identity ?? 100);
        default:
          return a.rank - b.rank;
      }
    });

    return designs;
  }, [result.designs, filterStatus, sortBy]);

  const approachInfo = APPROACH_LABELS[result.approach] || APPROACH_LABELS.joint;
  const passedCount = result.designs.filter(d => d.status === 'passed').length;
  const failedCount = result.designs.filter(d => d.status === 'failed').length;

  const toggleSelectDesign = (designId: string) => {
    const newSelected = new Set(selectedDesigns);
    if (newSelected.has(designId)) {
      newSelected.delete(designId);
    } else {
      newSelected.add(designId);
    }
    setSelectedDesigns(newSelected);
  };

  const selectAllPassed = () => {
    const passed = result.designs.filter(d => d.status === 'passed').map(d => d.id);
    setSelectedDesigns(new Set(passed));
  };

  const clearSelection = () => {
    setSelectedDesigns(new Set());
  };

  const getStatusBadge = (status: string) => {
    switch (status) {
      case 'passed':
        return <span className="px-2 py-0.5 text-xs font-medium rounded-md border border-border bg-background text-foreground">Passed</span>;
      case 'failed':
        return <span className="px-2 py-0.5 text-xs font-medium rounded-md border border-destructive/30 bg-destructive/5 text-destructive">Failed</span>;
      default:
        return <span className="px-2 py-0.5 text-xs font-medium rounded-md border border-border bg-muted text-muted-foreground">Pending</span>;
    }
  };

  const getPipelineBadge = (stage?: string) => {
    const stages = {
      backbone: 'Backbone',
      ligandmpnn: 'MPNN',
      validation: 'Validated',
      complete: 'Complete',
    };
    const label = stages[stage as keyof typeof stages] || stages.backbone;
    return <span className="px-2 py-0.5 text-xs font-medium rounded-md border border-border bg-muted/50 text-muted-foreground">{label}</span>;
  };

  // Error state
  if (result.status === 'failed') {
    return (
      <div className="rounded-lg border border-destructive/20 bg-destructive/5 p-4">
        <div className="flex items-center gap-3">
          <AlertCircle className="w-5 h-5 text-destructive" />
          <div>
            <h4 className="font-medium text-foreground">Design Failed</h4>
            <p className="text-sm text-muted-foreground mt-1">{result.error || 'An error occurred during design.'}</p>
          </div>
        </div>
      </div>
    );
  }

  // Running state
  if (result.status === 'running') {
    return (
      <div className="rounded-lg border border-border bg-muted/30 p-4">
        <div className="flex items-center gap-3">
          <Loader2 className="w-5 h-5 text-muted-foreground animate-spin" />
          <div>
            <h4 className="font-medium text-foreground">Designing...</h4>
            <p className="text-sm text-muted-foreground mt-1">Running {approachInfo.name} pipeline.</p>
          </div>
        </div>
      </div>
    );
  }

  // No designs
  if (result.designs.length === 0) {
    return (
      <div className="rounded-lg border border-border bg-muted/30 p-4">
        <div className="flex items-center gap-3">
          <AlertTriangle className="w-5 h-5 text-muted-foreground" />
          <div>
            <h4 className="font-medium text-foreground">No Designs Generated</h4>
            <p className="text-sm text-muted-foreground mt-1">
              The {approachInfo.name} approach did not produce any designs. Try adjusting parameters.
            </p>
          </div>
        </div>
      </div>
    );
  }

  return (
    <div className="space-y-4">
      {/* Header Summary */}
      <div className="rounded-lg border border-border bg-card p-4">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-3">
            <div className="w-10 h-10 rounded-lg bg-muted flex items-center justify-center">
              <approachInfo.Icon className="w-5 h-5 text-foreground" />
            </div>
            <div>
              <h3 className="font-medium text-foreground">{approachInfo.name} Results</h3>
              <p className="text-sm text-muted-foreground">
                {result.designs.length} design{result.designs.length !== 1 ? 's' : ''} generated
              </p>
            </div>
          </div>
          <div className="flex items-center gap-6">
            <div className="text-right">
              <div className="text-2xl font-semibold text-foreground">{passedCount}</div>
              <div className="text-xs text-muted-foreground">Passed</div>
            </div>
            <div className="text-right">
              <div className="text-2xl font-semibold text-muted-foreground">{failedCount}</div>
              <div className="text-xs text-muted-foreground">Failed</div>
            </div>
          </div>
        </div>
      </div>

      {/* Pipeline Controls */}
      {(onRunLigandMPNN || onRunValidation) && (
        <div className="rounded-lg border border-border bg-card p-4">
          <div className="flex items-center justify-between mb-3">
            <div className="flex items-center gap-2">
              <GitBranch className="w-4 h-4 text-muted-foreground" />
              <h4 className="font-medium text-foreground text-sm">Pipeline Actions</h4>
            </div>
            <div className="flex items-center gap-3 text-xs">
              <span className="text-muted-foreground">{selectedDesigns.size} selected</span>
              <button
                onClick={selectAllPassed}
                className="text-foreground underline-offset-4 hover:underline"
              >
                Select all passed
              </button>
              <button
                onClick={clearSelection}
                className="text-muted-foreground underline-offset-4 hover:underline"
              >
                Clear
              </button>
            </div>
          </div>

          <div className="flex gap-2">
            {onRunLigandMPNN && (
              <button
                onClick={() => {
                  const selected = result.designs.filter(d => selectedDesigns.has(d.id));
                  onRunLigandMPNN(selected.length > 0 ? selected : result.designs.filter(d => d.status === 'passed'));
                }}
                disabled={isProcessing}
                className={`flex-1 py-2 px-4 rounded-md font-medium text-sm flex items-center justify-center gap-2 transition-colors border ${
                  isProcessing
                    ? 'bg-muted text-muted-foreground border-border cursor-not-allowed'
                    : 'bg-background text-foreground border-border hover:bg-muted'
                }`}
              >
                <Type className="w-4 h-4" />
                Run LigandMPNN
              </button>
            )}
            {onRunValidation && (
              <button
                onClick={() => {
                  const selected = result.designs.filter(d => selectedDesigns.has(d.id));
                  onRunValidation(selected.length > 0 ? selected : result.designs.filter(d => d.status === 'passed'));
                }}
                disabled={isProcessing}
                className={`flex-1 py-2 px-4 rounded-md font-medium text-sm flex items-center justify-center gap-2 transition-colors border ${
                  isProcessing
                    ? 'bg-muted text-muted-foreground border-border cursor-not-allowed'
                    : 'bg-background text-foreground border-border hover:bg-muted'
                }`}
              >
                <ShieldCheck className="w-4 h-4" />
                Validate Structure
              </button>
            )}
          </div>
        </div>
      )}

      {/* Filter and Sort Controls */}
      <div className="flex items-center justify-between">
        <div className="flex items-center gap-1">
          {/* Status filter */}
          {(['all', 'passed', 'failed'] as const).map((status) => (
            <button
              key={status}
              onClick={() => setFilterStatus(status)}
              className={`px-3 py-1.5 text-xs font-medium rounded-md transition-colors ${
                filterStatus === status
                  ? 'bg-foreground text-background'
                  : 'text-muted-foreground hover:text-foreground hover:bg-muted'
              }`}
            >
              {status === 'all' ? 'All' : status === 'passed' ? `Passed (${passedCount})` : `Failed (${failedCount})`}
            </button>
          ))}
        </div>

        {/* Sort dropdown */}
        <div className="flex items-center gap-2">
          <span className="text-xs text-muted-foreground">Sort by:</span>
          <select
            value={sortBy}
            onChange={(e) => setSortBy(e.target.value as typeof sortBy)}
            className="text-xs border border-border rounded-md px-2 py-1.5 bg-background text-foreground"
          >
            <option value="rank">Rank</option>
            <option value="affinity">Affinity</option>
            <option value="contacts">Contacts</option>
            <option value="identity">Identity</option>
          </select>
        </div>
      </div>

      {/* Design Gallery */}
      <div className="rounded-lg border border-border bg-card overflow-hidden">
        <div className="px-4 py-3 border-b border-border flex items-center justify-between">
          <div className="flex items-center gap-2">
            <Layers className="w-4 h-4 text-muted-foreground" />
            <h4 className="font-medium text-foreground text-sm">Design Gallery</h4>
            <span className="text-xs text-muted-foreground">
              {filteredDesigns.length} of {result.designs.length}
            </span>
          </div>
        </div>

        <div className="divide-y divide-border max-h-[500px] overflow-y-auto">
          {filteredDesigns.map((design) => {
            const isSelected = design.id === selectedDesignId;
            const isChecked = selectedDesigns.has(design.id);

            return (
              <div
                key={design.id}
                className={`p-4 transition-colors ${
                  isSelected
                    ? 'bg-muted/50'
                    : 'bg-card hover:bg-muted/30'
                }`}
              >
                {/* Header row: Checkbox, Rank, Badges, Actions */}
                <div className="flex items-center justify-between mb-3">
                  <div className="flex items-center gap-3">
                    <input
                      type="checkbox"
                      checked={isChecked}
                      onChange={() => toggleSelectDesign(design.id)}
                      className="w-4 h-4 rounded border-border"
                    />
                    <button
                      onClick={() => onSelectDesign(design)}
                      className={`w-8 h-8 rounded-md flex items-center justify-center font-medium text-sm flex-shrink-0 transition-colors ${
                        isSelected ? 'bg-foreground text-background' : 'bg-muted text-muted-foreground hover:bg-muted/80'
                      }`}
                    >
                      #{design.rank}
                    </button>
                    {getStatusBadge(design.status)}
                    {getPipelineBadge(design.pipeline_stage)}
                  </div>
                  <div className="flex items-center gap-1">
                    <button
                      onClick={() => onSelectDesign(design)}
                      className="p-2 rounded-md hover:bg-muted transition-colors"
                      title="View in 3D"
                    >
                      <View className="w-4 h-4 text-muted-foreground" />
                    </button>
                    {onDownloadPdb && (
                      <button
                        onClick={() => onDownloadPdb(design)}
                        className="p-2 rounded-md hover:bg-muted transition-colors"
                        title="Download PDB"
                      >
                        <Download className="w-4 h-4 text-muted-foreground" />
                      </button>
                    )}
                  </div>
                </div>

                {/* Metrics row: Compact bordered cards */}
                <div className="flex flex-wrap gap-2">
                  <div className="px-3 py-2 rounded-md border border-border bg-background min-w-[100px]">
                    <div className="text-sm font-semibold text-foreground">
                      {design.metrics.affinity?.toFixed(2) ?? 'N/A'}
                    </div>
                    <div className="text-[10px] text-muted-foreground">Affinity (kcal/mol)</div>
                  </div>
                  <div className="px-3 py-2 rounded-md border border-border bg-background min-w-[80px]">
                    <div className="text-sm font-semibold text-foreground">
                      {design.metrics.contacts_a ?? 0} / {design.metrics.contacts_b ?? 0}
                    </div>
                    <div className="text-[10px] text-muted-foreground">Contacts (A/B)</div>
                  </div>
                  <div className="px-3 py-2 rounded-md border border-border bg-background min-w-[70px]">
                    <div className="text-sm font-semibold text-foreground">
                      {design.metrics.sequence_identity?.toFixed(1) ?? 'N/A'}%
                    </div>
                    <div className="text-[10px] text-muted-foreground">Identity</div>
                  </div>
                  <div className="px-3 py-2 rounded-md border border-border bg-background min-w-[70px]">
                    <div className="text-sm font-semibold text-foreground">
                      {design.metrics.anti_homo_score?.toFixed(0) ?? 'N/A'}
                    </div>
                    <div className="text-[10px] text-muted-foreground">Anti-Homo</div>
                  </div>
                  {(design.metrics.n7_hbonds !== undefined || design.metrics.n8_hbonds !== undefined) && (
                    <div className="px-3 py-2 rounded-md border border-border bg-background min-w-[90px]">
                      <div className="text-sm font-semibold text-foreground">
                        {design.metrics.n7_hbonds ?? 0} / {design.metrics.n8_hbonds ?? 0}
                      </div>
                      <div className="text-[10px] text-muted-foreground">H-bonds (N7/N8)</div>
                    </div>
                  )}
                </div>
              </div>
            );
          })}
        </div>
      </div>

      {/* Selected Design Metrics */}
      {selectedDesignId && (
        <div>
          {(() => {
            const design = result.designs.find(d => d.id === selectedDesignId);
            if (!design) return null;
            return (
              <SuccessCriteriaCard
                metrics={design.metrics}
                filterMode={filterMode}
                onFilterModeChange={setFilterMode}
                showToggle
              />
            );
          })()}
        </div>
      )}
    </div>
  );
}
