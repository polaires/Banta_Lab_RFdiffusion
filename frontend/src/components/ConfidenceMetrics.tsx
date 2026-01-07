'use client';

import { useMemo } from 'react';
import type { ConfidenceMetrics, RMSDResult } from '@/lib/api';
import { CheckCircle, AlertCircle, XCircle, Info, Activity, Target, Layers } from 'lucide-react';

interface ConfidenceMetricsDisplayProps {
  confidences: ConfidenceMetrics | null;
  rmsdResult?: RMSDResult | null;
  compact?: boolean;
}

// Color scale for pLDDT (AlphaFold style)
function getPLDDTColor(value: number): string {
  if (value >= 0.9) return '#0053d6';  // Very high - blue
  if (value >= 0.7) return '#65cbf3';  // High - cyan
  if (value >= 0.5) return '#ffdb13';  // Medium - yellow
  return '#ff7d45';  // Low - orange
}

function getPLDDTLabel(value: number): string {
  if (value >= 0.9) return 'Very High';
  if (value >= 0.7) return 'High';
  if (value >= 0.5) return 'Medium';
  return 'Low';
}

function getRMSDIcon(interpretation: string) {
  switch (interpretation) {
    case 'Excellent':
      return <CheckCircle className="w-5 h-5 text-green-400" />;
    case 'Good':
      return <CheckCircle className="w-5 h-5 text-blue-400" />;
    case 'Moderate':
      return <AlertCircle className="w-5 h-5 text-yellow-400" />;
    default:
      return <XCircle className="w-5 h-5 text-red-400" />;
  }
}

function getRMSDColor(interpretation: string): string {
  switch (interpretation) {
    case 'Excellent': return 'text-green-400';
    case 'Good': return 'text-blue-400';
    case 'Moderate': return 'text-yellow-400';
    default: return 'text-red-400';
  }
}

export function ConfidenceMetricsDisplay({ confidences, rmsdResult, compact = false }: ConfidenceMetricsDisplayProps) {
  if (!confidences && !rmsdResult) return null;

  const summary = confidences?.summary_confidences;
  const perResiduePLDDT = confidences?.per_residue_plddt;

  // Calculate average pLDDT from per-residue if summary not available
  const avgPLDDT = useMemo(() => {
    if (summary?.overall_plddt) return summary.overall_plddt;
    if (perResiduePLDDT && perResiduePLDDT.length > 0) {
      return perResiduePLDDT.reduce((a, b) => a + b, 0) / perResiduePLDDT.length;
    }
    return null;
  }, [summary, perResiduePLDDT]);

  if (compact) {
    return (
      <div className="flex flex-wrap gap-3 text-sm">
        {avgPLDDT !== null && (
          <div className="flex items-center gap-1.5">
            <div
              className="w-3 h-3 rounded-full"
              style={{ backgroundColor: getPLDDTColor(avgPLDDT) }}
            />
            <span className="text-gray-400">pLDDT:</span>
            <span className="font-medium">{(avgPLDDT * 100).toFixed(1)}</span>
          </div>
        )}
        {summary?.ptm && (
          <div className="flex items-center gap-1.5">
            <Target className="w-3 h-3 text-gray-400" />
            <span className="text-gray-400">pTM:</span>
            <span className="font-medium">{(summary.ptm * 100).toFixed(1)}</span>
          </div>
        )}
        {rmsdResult && (
          <div className="flex items-center gap-1.5">
            {getRMSDIcon(rmsdResult.interpretation)}
            <span className="text-gray-400">RMSD:</span>
            <span className={`font-medium ${getRMSDColor(rmsdResult.interpretation)}`}>
              {rmsdResult.rmsd.toFixed(2)}A
            </span>
          </div>
        )}
      </div>
    );
  }

  return (
    <div className="space-y-4 p-4 bg-gray-800/50 rounded-lg border border-gray-700">
      <h3 className="font-medium flex items-center gap-2">
        <Activity className="w-4 h-4" />
        Confidence Metrics
      </h3>

      {/* Summary Metrics Grid */}
      {summary && (
        <div className="grid grid-cols-2 md:grid-cols-3 gap-3">
          {/* Overall pLDDT */}
          <div className="p-3 bg-gray-700/50 rounded">
            <div className="flex items-center gap-2 mb-1">
              <div
                className="w-3 h-3 rounded-full"
                style={{ backgroundColor: getPLDDTColor(summary.overall_plddt) }}
              />
              <span className="text-xs text-gray-400">pLDDT</span>
            </div>
            <div className="text-xl font-bold">
              {(summary.overall_plddt * 100).toFixed(1)}
            </div>
            <div className="text-xs text-gray-500">
              {getPLDDTLabel(summary.overall_plddt)}
            </div>
          </div>

          {/* pTM Score */}
          <div className="p-3 bg-gray-700/50 rounded">
            <div className="flex items-center gap-2 mb-1">
              <Target className="w-3 h-3 text-gray-400" />
              <span className="text-xs text-gray-400">pTM</span>
            </div>
            <div className="text-xl font-bold">
              {(summary.ptm * 100).toFixed(1)}
            </div>
            <div className="text-xs text-gray-500">
              Template Modeling
            </div>
          </div>

          {/* PAE */}
          <div className="p-3 bg-gray-700/50 rounded">
            <div className="flex items-center gap-2 mb-1">
              <Layers className="w-3 h-3 text-gray-400" />
              <span className="text-xs text-gray-400">PAE</span>
            </div>
            <div className="text-xl font-bold">
              {summary.overall_pae.toFixed(1)}A
            </div>
            <div className="text-xs text-gray-500">
              Predicted Aligned Error
            </div>
          </div>

          {/* Ranking Score */}
          {summary.ranking_score !== undefined && (
            <div className="p-3 bg-gray-700/50 rounded">
              <div className="flex items-center gap-2 mb-1">
                <CheckCircle className="w-3 h-3 text-gray-400" />
                <span className="text-xs text-gray-400">Ranking</span>
              </div>
              <div className="text-xl font-bold">
                {(summary.ranking_score * 100).toFixed(1)}
              </div>
              <div className="text-xs text-gray-500">
                Overall Quality
              </div>
            </div>
          )}

          {/* ipTM (if available) */}
          {summary.iptm !== null && summary.iptm !== undefined && (
            <div className="p-3 bg-gray-700/50 rounded">
              <div className="flex items-center gap-2 mb-1">
                <Target className="w-3 h-3 text-gray-400" />
                <span className="text-xs text-gray-400">ipTM</span>
              </div>
              <div className="text-xl font-bold">
                {(summary.iptm * 100).toFixed(1)}
              </div>
              <div className="text-xs text-gray-500">
                Interface pTM
              </div>
            </div>
          )}

          {/* Clash indicator */}
          {summary.has_clash !== undefined && (
            <div className="p-3 bg-gray-700/50 rounded">
              <div className="flex items-center gap-2 mb-1">
                {summary.has_clash ? (
                  <XCircle className="w-3 h-3 text-red-400" />
                ) : (
                  <CheckCircle className="w-3 h-3 text-green-400" />
                )}
                <span className="text-xs text-gray-400">Clashes</span>
              </div>
              <div className={`text-xl font-bold ${summary.has_clash ? 'text-red-400' : 'text-green-400'}`}>
                {summary.has_clash ? 'Yes' : 'None'}
              </div>
            </div>
          )}
        </div>
      )}

      {/* RMSD Result */}
      {rmsdResult && (
        <div className="p-3 bg-gray-700/50 rounded border border-gray-600">
          <div className="flex items-center justify-between mb-2">
            <div className="flex items-center gap-2">
              {getRMSDIcon(rmsdResult.interpretation)}
              <span className="font-medium">RMSD Validation</span>
            </div>
            <span className={`text-lg font-bold ${getRMSDColor(rmsdResult.interpretation)}`}>
              {rmsdResult.rmsd.toFixed(3)} A
            </span>
          </div>
          <div className={`text-sm ${getRMSDColor(rmsdResult.interpretation)}`}>
            {rmsdResult.interpretation} Designability
          </div>
          <p className="text-xs text-gray-400 mt-1">
            {rmsdResult.description}
          </p>
          <p className="text-xs text-gray-500 mt-1">
            {rmsdResult.num_atoms_compared} backbone atoms compared
          </p>
        </div>
      )}

      {/* Per-residue pLDDT visualization */}
      {perResiduePLDDT && perResiduePLDDT.length > 0 && (
        <div>
          <div className="flex items-center justify-between mb-2">
            <span className="text-xs text-gray-400">Per-residue pLDDT</span>
            <div className="flex items-center gap-2 text-xs">
              <span className="flex items-center gap-1">
                <div className="w-2 h-2 rounded" style={{ backgroundColor: '#ff7d45' }} />
                &lt;50
              </span>
              <span className="flex items-center gap-1">
                <div className="w-2 h-2 rounded" style={{ backgroundColor: '#ffdb13' }} />
                50-70
              </span>
              <span className="flex items-center gap-1">
                <div className="w-2 h-2 rounded" style={{ backgroundColor: '#65cbf3' }} />
                70-90
              </span>
              <span className="flex items-center gap-1">
                <div className="w-2 h-2 rounded" style={{ backgroundColor: '#0053d6' }} />
                &gt;90
              </span>
            </div>
          </div>
          <div className="flex h-6 rounded overflow-hidden">
            {perResiduePLDDT.slice(0, 200).map((value, i) => (
              <div
                key={i}
                className="flex-1 min-w-[2px]"
                style={{ backgroundColor: getPLDDTColor(value) }}
                title={`Residue ${i + 1}: ${(value * 100).toFixed(1)}`}
              />
            ))}
          </div>
          <div className="flex justify-between text-xs text-gray-500 mt-1">
            <span>1</span>
            <span>{Math.min(perResiduePLDDT.length, 200)}</span>
          </div>
        </div>
      )}

      {/* Help text */}
      <div className="flex items-start gap-2 text-xs text-gray-500">
        <Info className="w-3 h-3 mt-0.5 flex-shrink-0" />
        <div>
          <strong>pLDDT</strong>: Per-residue confidence (0-100, higher = better).
          <strong className="ml-2">pTM</strong>: Template modeling score.
          <strong className="ml-2">PAE</strong>: Predicted alignment error (lower = better).
        </div>
      </div>
    </div>
  );
}

export default ConfidenceMetricsDisplay;
