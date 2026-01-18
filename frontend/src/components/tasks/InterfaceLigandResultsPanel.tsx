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

// Approach display names
const APPROACH_LABELS: Record<string, { name: string; Icon: React.ComponentType<{ className?: string }>; color: string }> = {
  joint: { name: 'Joint Heterodimer', Icon: Network, color: 'purple' },
  asymmetric_rasa: { name: 'Asymmetric RASA', Icon: SlidersHorizontal, color: 'blue' },
  induced: { name: 'Induced Dimerization', Icon: GitBranch, color: 'cyan' },
  asymmetric: { name: 'Single Chain', Icon: Link, color: 'amber' },
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
        return <span className="px-2 py-0.5 text-xs font-medium rounded-full bg-green-100 text-green-700">Passed</span>;
      case 'failed':
        return <span className="px-2 py-0.5 text-xs font-medium rounded-full bg-red-100 text-red-700">Failed</span>;
      default:
        return <span className="px-2 py-0.5 text-xs font-medium rounded-full bg-muted text-foreground">Pending</span>;
    }
  };

  const getPipelineBadge = (stage?: string) => {
    const stages = {
      backbone: { label: 'Backbone', color: 'bg-blue-100 text-blue-700' },
      ligandmpnn: { label: 'MPNN', color: 'bg-purple-100 text-purple-700' },
      validation: { label: 'Validated', color: 'bg-cyan-100 text-cyan-700' },
      complete: { label: 'Complete', color: 'bg-green-100 text-green-700' },
    };
    const info = stages[stage as keyof typeof stages] || stages.backbone;
    return <span className={`px-2 py-0.5 text-xs font-medium rounded-full ${info.color}`}>{info.label}</span>;
  };

  // Error state
  if (result.status === 'failed') {
    return (
      <div className="bg-red-50 rounded-xl border border-red-200 p-6">
        <div className="flex items-center gap-3">
          <AlertCircle className="w-6 h-6 text-red-500" />
          <div>
            <h4 className="font-semibold text-red-900">Design Failed</h4>
            <p className="text-sm text-red-700 mt-1">{result.error || 'An error occurred during design.'}</p>
          </div>
        </div>
      </div>
    );
  }

  // Running state
  if (result.status === 'running') {
    return (
      <div className="bg-blue-50 rounded-xl border border-blue-200 p-6">
        <div className="flex items-center gap-3">
          <Loader2 className="w-6 h-6 text-blue-500 animate-spin" />
          <div>
            <h4 className="font-semibold text-blue-900">Designing...</h4>
            <p className="text-sm text-blue-700 mt-1">Running {approachInfo.name} pipeline. This may take a few minutes.</p>
          </div>
        </div>
      </div>
    );
  }

  // No designs
  if (result.designs.length === 0) {
    return (
      <div className="bg-amber-50 rounded-xl border border-amber-200 p-6">
        <div className="flex items-center gap-3">
          <AlertTriangle className="w-6 h-6 text-amber-500" />
          <div>
            <h4 className="font-semibold text-amber-900">No Designs Generated</h4>
            <p className="text-sm text-amber-700 mt-1">
              The {approachInfo.name} approach did not produce any designs. Try adjusting parameters.
            </p>
          </div>
        </div>
      </div>
    );
  }

  return (
    <div className="space-y-6">
      {/* Header Summary */}
      <div className={`bg-gradient-to-r from-${approachInfo.color}-50 to-${approachInfo.color}-100 rounded-xl border border-${approachInfo.color}-200 p-4`}>
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-3">
            <div className={`w-10 h-10 rounded-lg bg-${approachInfo.color}-500 flex items-center justify-center`}>
              <approachInfo.Icon className="w-5 h-5 text-white" />
            </div>
            <div>
              <h3 className="font-semibold text-foreground">{approachInfo.name} Results</h3>
              <p className="text-sm text-muted-foreground">
                {result.designs.length} design{result.designs.length !== 1 ? 's' : ''} generated
              </p>
            </div>
          </div>
          <div className="flex items-center gap-4">
            <div className="text-right">
              <div className="text-2xl font-bold text-green-600">{passedCount}</div>
              <div className="text-xs text-muted-foreground">Passed</div>
            </div>
            <div className="text-right">
              <div className="text-2xl font-bold text-red-500">{failedCount}</div>
              <div className="text-xs text-muted-foreground">Failed</div>
            </div>
          </div>
        </div>
      </div>

      {/* Pipeline Controls */}
      {(onRunLigandMPNN || onRunValidation) && (
        <div className="bg-card rounded-xl border border-border p-4">
          <div className="flex items-center justify-between mb-3">
            <div className="flex items-center gap-2">
              <GitBranch className="w-5 h-5 text-muted-foreground" />
              <h4 className="font-semibold text-foreground text-sm">Pipeline Actions</h4>
            </div>
            <div className="flex items-center gap-2 text-xs text-muted-foreground">
              <span>{selectedDesigns.size} selected</span>
              <button
                onClick={selectAllPassed}
                className="text-purple-600 hover:underline"
              >
                Select all passed
              </button>
              <button
                onClick={clearSelection}
                className="text-muted-foreground hover:underline"
              >
                Clear
              </button>
            </div>
          </div>

          <div className="flex gap-3">
            {onRunLigandMPNN && (
              <button
                onClick={() => {
                  const selected = result.designs.filter(d => selectedDesigns.has(d.id));
                  onRunLigandMPNN(selected.length > 0 ? selected : result.designs.filter(d => d.status === 'passed'));
                }}
                disabled={isProcessing}
                className={`flex-1 py-2.5 px-4 rounded-lg font-medium text-sm flex items-center justify-center gap-2 transition-all ${
                  isProcessing
                    ? 'bg-muted text-muted-foreground cursor-not-allowed'
                    : 'bg-purple-100 text-purple-700 hover:bg-purple-200'
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
                className={`flex-1 py-2.5 px-4 rounded-lg font-medium text-sm flex items-center justify-center gap-2 transition-all ${
                  isProcessing
                    ? 'bg-muted text-muted-foreground cursor-not-allowed'
                    : 'bg-cyan-100 text-cyan-700 hover:bg-cyan-200'
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
        <div className="flex items-center gap-2">
          {/* Status filter */}
          <div className="flex rounded-lg border border-border overflow-hidden">
            {(['all', 'passed', 'failed'] as const).map((status) => (
              <button
                key={status}
                onClick={() => setFilterStatus(status)}
                className={`px-3 py-1.5 text-xs font-medium transition-colors ${
                  filterStatus === status
                    ? 'bg-purple-100 text-purple-700'
                    : 'bg-card text-muted-foreground hover:bg-muted/50'
                }`}
              >
                {status === 'all' ? 'All' : status === 'passed' ? `Passed (${passedCount})` : `Failed (${failedCount})`}
              </button>
            ))}
          </div>
        </div>

        {/* Sort dropdown */}
        <div className="flex items-center gap-2">
          <span className="text-xs text-muted-foreground">Sort by:</span>
          <select
            value={sortBy}
            onChange={(e) => setSortBy(e.target.value as typeof sortBy)}
            className="text-xs border border-border rounded-lg px-2 py-1.5 bg-card"
          >
            <option value="rank">Rank</option>
            <option value="affinity">Affinity</option>
            <option value="contacts">Contacts</option>
            <option value="identity">Identity</option>
          </select>
        </div>
      </div>

      {/* Design Gallery */}
      <div className="bg-card rounded-xl border border-border overflow-hidden">
        <div className="bg-muted/50 px-4 py-3 border-b border-border flex items-center justify-between">
          <div className="flex items-center gap-2">
            <Layers className="w-5 h-5 text-purple-600" />
            <h4 className="font-semibold text-foreground text-sm">Design Gallery</h4>
            <span className="text-xs text-muted-foreground bg-muted px-2 py-0.5 rounded-full">
              {filteredDesigns.length} of {result.designs.length}
            </span>
          </div>
        </div>

        <div className="p-3 space-y-2 max-h-[500px] overflow-y-auto">
          {filteredDesigns.map((design) => {
            const isSelected = design.id === selectedDesignId;
            const isChecked = selectedDesigns.has(design.id);
            const totalContacts = (design.metrics.contacts_a ?? 0) + (design.metrics.contacts_b ?? 0);

            return (
              <div
                key={design.id}
                className={`p-4 rounded-lg border-2 transition-all ${
                  isSelected
                    ? 'border-purple-500 bg-purple-50'
                    : 'border-border bg-card hover:border-purple-300 hover:bg-purple-50/30'
                }`}
              >
                <div className="flex items-start gap-3">
                  {/* Checkbox */}
                  <input
                    type="checkbox"
                    checked={isChecked}
                    onChange={() => toggleSelectDesign(design.id)}
                    className="mt-1 w-4 h-4 rounded border-border text-purple-600"
                  />

                  {/* Rank badge */}
                  <button
                    onClick={() => onSelectDesign(design)}
                    className={`w-8 h-8 rounded-lg flex items-center justify-center font-bold text-sm flex-shrink-0 ${
                      isSelected ? 'bg-purple-500 text-white' : 'bg-muted text-muted-foreground hover:bg-purple-200'
                    }`}
                  >
                    #{design.rank}
                  </button>

                  {/* Main content */}
                  <div className="flex-1 min-w-0">
                    <div className="flex items-center gap-2 mb-2">
                      {getStatusBadge(design.status)}
                      {getPipelineBadge(design.pipeline_stage)}
                    </div>

                    {/* Metrics grid */}
                    <div className="grid grid-cols-4 gap-3 text-xs">
                      <div>
                        <div className="text-muted-foreground text-[10px] uppercase">Affinity</div>
                        <div className={`font-bold ${
                          (design.metrics.affinity ?? 0) < -5 ? 'text-green-600' :
                          (design.metrics.affinity ?? 0) < -2 ? 'text-amber-600' : 'text-red-600'
                        }`}>
                          {design.metrics.affinity?.toFixed(2) ?? 'N/A'} kcal/mol
                        </div>
                      </div>
                      <div>
                        <div className="text-muted-foreground text-[10px] uppercase">Contacts</div>
                        <div className="font-medium text-foreground">
                          <span className="text-purple-600">{design.metrics.contacts_a ?? 0}</span>
                          <span className="text-muted-foreground"> / </span>
                          <span className="text-cyan-600">{design.metrics.contacts_b ?? 0}</span>
                        </div>
                      </div>
                      <div>
                        <div className="text-muted-foreground text-[10px] uppercase">Identity</div>
                        <div className={`font-medium ${
                          (design.metrics.sequence_identity ?? 100) < 50 ? 'text-green-600' :
                          (design.metrics.sequence_identity ?? 100) < 70 ? 'text-amber-600' : 'text-red-600'
                        }`}>
                          {design.metrics.sequence_identity?.toFixed(1) ?? 'N/A'}%
                        </div>
                      </div>
                      <div>
                        <div className="text-muted-foreground text-[10px] uppercase">Anti-Homo</div>
                        <div className={`font-medium ${
                          (design.metrics.anti_homo_score ?? 0) > 60 ? 'text-green-600' :
                          (design.metrics.anti_homo_score ?? 0) > 40 ? 'text-amber-600' : 'text-red-600'
                        }`}>
                          {design.metrics.anti_homo_score?.toFixed(0) ?? 'N/A'}/100
                        </div>
                      </div>
                    </div>

                    {/* H-bond info */}
                    {(design.metrics.n7_hbonds !== undefined || design.metrics.n8_hbonds !== undefined) && (
                      <div className="mt-2 flex items-center gap-3 text-xs">
                        <span className="text-muted-foreground">H-bonds:</span>
                        <span className={`font-medium ${(design.metrics.n7_hbonds ?? 0) >= 1 ? 'text-green-600' : 'text-red-600'}`}>
                          N7: {design.metrics.n7_hbonds ?? 0}
                        </span>
                        <span className={`font-medium ${(design.metrics.n8_hbonds ?? 0) >= 1 ? 'text-green-600' : 'text-red-600'}`}>
                          N8: {design.metrics.n8_hbonds ?? 0}
                        </span>
                      </div>
                    )}
                  </div>

                  {/* Actions */}
                  <div className="flex items-center gap-1">
                    <button
                      onClick={() => onSelectDesign(design)}
                      className="p-2 rounded-lg hover:bg-muted transition-colors"
                      title="View in 3D"
                    >
                      <View className="w-4 h-4 text-muted-foreground" />
                    </button>
                    {onDownloadPdb && (
                      <button
                        onClick={() => onDownloadPdb(design)}
                        className="p-2 rounded-lg hover:bg-muted transition-colors"
                        title="Download PDB"
                      >
                        <Download className="w-4 h-4 text-muted-foreground" />
                      </button>
                    )}
                  </div>
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
