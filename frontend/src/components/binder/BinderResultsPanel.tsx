'use client';

import { useState } from 'react';
import {
  AlertCircle,
  Loader2,
  AlertTriangle,
  CheckCircle2,
  LayoutGrid,
  Copy,
  Check,
  Download,
  BadgeCheck,
  Dna,
} from 'lucide-react';
import { PipelineFunnel } from './PipelineFunnel';
import { InterfaceMetrics } from './InterfaceMetrics';

// Design result from the protein_binder_design API
export interface BinderDesign {
  binder_sequence: string;
  mpnn_score: number;
  esm_perplexity: number;
  esm_confidence: number;
  interface_contacts: number;
  interface_hbonds: number;
  buried_sasa: number;
  packstat: number;
  rank: number;
  pdb_content?: string;
  // New BindCraft-style metrics
  shape_complementarity?: number;
  surface_hydrophobicity?: number;
  unsaturated_hbonds?: number;
  interface_residue_count?: number;
  // ESMFold validation metrics
  esmfold_plddt?: number;
  esmfold_rmsd?: number;
  esmfold_validation_passed?: boolean;
}

export interface BinderStatistics {
  generated: number;
  mpnn_designed: number;
  esm_passed: number;
  relaxed: number;
  interface_analyzed: number;
  passed_filters: number;
  returned: number;
}

interface BinderResultsPanelProps {
  status: 'completed' | 'error' | 'processing';
  statistics: BinderStatistics;
  designs: BinderDesign[];
  onDesignSelect?: (design: BinderDesign) => void;
  onDownloadPdb?: (design: BinderDesign) => void;
  onCopySequence?: (sequence: string) => void;
}

export function BinderResultsPanel({
  status,
  statistics,
  designs,
  onDesignSelect,
  onDownloadPdb,
  onCopySequence,
}: BinderResultsPanelProps) {
  const [selectedDesignIndex, setSelectedDesignIndex] = useState(0);
  const [copiedIndex, setCopiedIndex] = useState<number | null>(null);

  const selectedDesign = designs[selectedDesignIndex] || null;

  const handleCopySequence = (sequence: string, index: number) => {
    navigator.clipboard.writeText(sequence);
    setCopiedIndex(index);
    setTimeout(() => setCopiedIndex(null), 2000);
    onCopySequence?.(sequence);
  };

  const handleSelectDesign = (design: BinderDesign, index: number) => {
    setSelectedDesignIndex(index);
    onDesignSelect?.(design);
  };

  if (status === 'error') {
    return (
      <div className="bg-red-50 rounded-xl border border-red-200 p-6">
        <div className="flex items-center gap-3">
          <AlertCircle className="h-6 w-6 text-red-500" />
          <div>
            <h4 className="font-semibold text-red-900">Pipeline Error</h4>
            <p className="text-sm text-red-700">
              The binder design pipeline encountered an error. Please try again.
            </p>
          </div>
        </div>
      </div>
    );
  }

  if (status === 'processing') {
    return (
      <div className="bg-slate-50 rounded-xl border border-slate-200 p-6">
        <div className="flex items-center gap-3">
          <Loader2 className="h-6 w-6 text-blue-500 animate-spin" />
          <div>
            <h4 className="font-semibold text-slate-900">Designing Binders...</h4>
            <p className="text-sm text-slate-600">
              Running multi-stage validation pipeline. This may take a few minutes.
            </p>
          </div>
        </div>
      </div>
    );
  }

  if (designs.length === 0) {
    return (
      <div className="bg-amber-50 rounded-xl border border-amber-200 p-6">
        <div className="flex items-center gap-3">
          <AlertTriangle className="h-6 w-6 text-amber-500" />
          <div>
            <h4 className="font-semibold text-amber-900">No Designs Passed Filters</h4>
            <p className="text-sm text-amber-700">
              Try adjusting parameters or using a different target region.
            </p>
          </div>
        </div>
        {/* Still show the funnel to explain what happened */}
        <div className="mt-4">
          <PipelineFunnel statistics={statistics} />
        </div>
      </div>
    );
  }

  return (
    <div className="space-y-6">
      {/* Header with Summary */}
      <div className="bg-gradient-to-r from-green-50 to-emerald-50 rounded-xl border border-green-200 p-4">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-3">
            <div className="w-10 h-10 rounded-lg bg-green-500 flex items-center justify-center">
              <CheckCircle2 className="h-6 w-6 text-white" />
            </div>
            <div>
              <h3 className="font-semibold text-slate-900">
                {designs.length} High-Quality Binder{designs.length !== 1 ? 's' : ''} Designed
              </h3>
              <p className="text-sm text-slate-600">
                {statistics.returned} of {statistics.generated} designs passed all filters
              </p>
            </div>
          </div>
          <div className="text-right">
            <div className="text-2xl font-bold text-green-600">
              {((statistics.returned / statistics.generated) * 100).toFixed(0)}%
            </div>
            <div className="text-xs text-slate-500">Pass Rate</div>
          </div>
        </div>
      </div>

      {/* Pipeline Funnel */}
      <PipelineFunnel statistics={statistics} />

      {/* Design Gallery */}
      <div className="bg-white rounded-xl border border-slate-200 overflow-hidden">
        <div className="bg-slate-50 px-4 py-3 border-b border-slate-200">
          <div className="flex items-center gap-2">
            <LayoutGrid className="h-5 w-5 text-violet-600" />
            <h4 className="font-semibold text-slate-900 text-sm">Design Gallery</h4>
            <span className="text-xs text-slate-500 bg-slate-200 px-2 py-0.5 rounded-full">
              {designs.length} designs
            </span>
          </div>
        </div>

        <div className="p-4">
          {/* Design Cards */}
          <div className="grid grid-cols-1 sm:grid-cols-2 gap-3 mb-4">
            {designs.map((design, index) => (
              <button
                key={index}
                onClick={() => handleSelectDesign(design, index)}
                className={`text-left p-4 rounded-lg border-2 transition-all ${
                  selectedDesignIndex === index
                    ? 'border-violet-500 bg-violet-50'
                    : 'border-slate-200 bg-white hover:border-violet-300'
                }`}
              >
                <div className="flex items-center justify-between mb-2">
                  <div className={`w-8 h-8 rounded-lg flex items-center justify-center font-bold text-sm ${
                    selectedDesignIndex === index
                      ? 'bg-violet-500 text-white'
                      : 'bg-slate-200 text-slate-600'
                  }`}>
                    #{design.rank}
                  </div>
                  <div className="flex items-center gap-2">
                    <button
                      onClick={(e) => {
                        e.stopPropagation();
                        handleCopySequence(design.binder_sequence, index);
                      }}
                      className="p-1.5 rounded hover:bg-slate-100 transition-colors"
                      title="Copy sequence"
                    >
                      {copiedIndex === index ? (
                        <Check className="h-4 w-4 text-slate-500" />
                      ) : (
                        <Copy className="h-4 w-4 text-slate-500" />
                      )}
                    </button>
                    {design.pdb_content && onDownloadPdb && (
                      <button
                        onClick={(e) => {
                          e.stopPropagation();
                          onDownloadPdb(design);
                        }}
                        className="p-1.5 rounded hover:bg-slate-100 transition-colors"
                        title="Download PDB"
                      >
                        <Download className="h-4 w-4 text-slate-500" />
                      </button>
                    )}
                  </div>
                </div>

                {/* Quick Metrics */}
                <div className="grid grid-cols-2 gap-2 text-xs">
                  <div>
                    <span className="text-slate-400">Contacts:</span>
                    <span className="ml-1 font-medium text-slate-700">
                      {design.interface_contacts.toLocaleString()}
                    </span>
                  </div>
                  <div>
                    <span className="text-slate-400">H-bonds:</span>
                    <span className="ml-1 font-medium text-slate-700">{design.interface_hbonds}</span>
                  </div>
                  <div>
                    <span className="text-slate-400">ESM:</span>
                    <span className={`ml-1 font-medium ${
                      design.esm_perplexity < 5 ? 'text-green-600' :
                      design.esm_perplexity < 8 ? 'text-blue-600' : 'text-amber-600'
                    }`}>
                      {design.esm_perplexity.toFixed(2)}
                    </span>
                  </div>
                  {design.shape_complementarity != null ? (
                    <div>
                      <span className="text-slate-400">SC:</span>
                      <span className={`ml-1 font-medium ${
                        design.shape_complementarity >= 0.6 ? 'text-green-600' :
                        design.shape_complementarity >= 0.5 ? 'text-blue-600' : 'text-amber-600'
                      }`}>
                        {design.shape_complementarity.toFixed(2)}
                      </span>
                    </div>
                  ) : (
                    <div>
                      <span className="text-slate-400">Length:</span>
                      <span className="ml-1 font-medium text-slate-700">
                        {design.binder_sequence.length} aa
                      </span>
                    </div>
                  )}
                </div>
                {/* ESMFold validation badge if available */}
                {design.esmfold_plddt != null && (
                  <div className="mt-2 flex items-center gap-2">
                    <span className={`inline-flex items-center gap-1 px-2 py-0.5 rounded text-xs font-medium ${
                      design.esmfold_validation_passed
                        ? 'bg-green-100 text-green-700'
                        : 'bg-amber-100 text-amber-700'
                    }`}>
                      {design.esmfold_validation_passed ? (
                        <BadgeCheck className="h-3 w-3" />
                      ) : (
                        <AlertTriangle className="h-3 w-3" />
                      )}
                      pLDDT: {(design.esmfold_plddt * 100).toFixed(0)}%
                    </span>
                    {design.esmfold_rmsd != null && (
                      <span className="text-xs text-slate-500">
                        RMSD: {design.esmfold_rmsd.toFixed(1)}Å
                      </span>
                    )}
                  </div>
                )}

                {/* Sequence Preview */}
                <div className="mt-2 p-2 bg-slate-100 rounded text-xs font-mono text-slate-600 truncate">
                  {design.binder_sequence.substring(0, 40)}...
                </div>
              </button>
            ))}
          </div>
        </div>
      </div>

      {/* Selected Design Details */}
      {selectedDesign && (
        <div className="space-y-4">
          {/* Interface Metrics */}
          <InterfaceMetrics
            metrics={{
              interface_contacts: selectedDesign.interface_contacts,
              interface_hbonds: selectedDesign.interface_hbonds,
              buried_sasa: selectedDesign.buried_sasa,
              packstat: selectedDesign.packstat,
              shape_complementarity: selectedDesign.shape_complementarity,
              surface_hydrophobicity: selectedDesign.surface_hydrophobicity,
              unsaturated_hbonds: selectedDesign.unsaturated_hbonds,
              interface_residue_count: selectedDesign.interface_residue_count,
            }}
            esm={{
              perplexity: selectedDesign.esm_perplexity,
              confidence: selectedDesign.esm_confidence,
            }}
            mpnn_score={selectedDesign.mpnn_score}
            esmfold={selectedDesign.esmfold_plddt != null ? {
              plddt: selectedDesign.esmfold_plddt,
              rmsd: selectedDesign.esmfold_rmsd,
              validation_passed: selectedDesign.esmfold_validation_passed ?? false,
            } : undefined}
          />

          {/* Full Sequence */}
          <div className="bg-white rounded-xl border border-slate-200 overflow-hidden">
            <div className="bg-slate-50 px-4 py-3 border-b border-slate-200 flex items-center justify-between">
              <div className="flex items-center gap-2">
                <Dna className="h-5 w-5 text-slate-600" />
                <h4 className="font-semibold text-slate-900 text-sm">Binder Sequence</h4>
                <span className="text-xs text-slate-500">
                  {selectedDesign.binder_sequence.length} residues
                </span>
              </div>
              <button
                onClick={() => handleCopySequence(selectedDesign.binder_sequence, selectedDesignIndex)}
                className="flex items-center gap-1 px-3 py-1 text-xs font-medium text-violet-600 hover:bg-violet-50 rounded-lg transition-colors"
              >
                {copiedIndex === selectedDesignIndex ? (
                  <Check className="h-4 w-4" />
                ) : (
                  <Copy className="h-4 w-4" />
                )}
                {copiedIndex === selectedDesignIndex ? 'Copied!' : 'Copy'}
              </button>
            </div>
            <div className="p-4">
              <div className="font-mono text-sm text-slate-700 break-all leading-relaxed bg-slate-50 p-3 rounded-lg">
                {selectedDesign.binder_sequence}
              </div>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}
