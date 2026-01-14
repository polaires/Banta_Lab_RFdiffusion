'use client';

// Interaction profile data from PLIP analysis
export interface InteractionProfile {
  hbonds: number;
  hydrophobic: number;
  pi_stacking: number;
  salt_bridges: number;
  halogen_bonds?: number;
  total: number;
}

// Binder evaluation result interface
export interface BinderEvaluation {
  success: boolean;
  statistics?: {
    generated: number;
    rg_filtered?: number;  // Elongated binders filtered out
    mpnn_designed: number;
    esm_passed: number;
    relaxed: number;
    interface_analyzed: number;
    passed_filters: number;
    returned: number;
  };
  best_design?: {
    rank: number;
    binder_sequence: string;
    mpnn_score?: number;
    esm_perplexity?: number;
    esm_confidence?: number;
    interface_contacts?: number;
    interface_hbonds?: number;
    buried_sasa?: number;
    packstat?: number;
    rg_ratio?: number;  // Radius of gyration ratio
    rg_assessment?: string;  // "compact", "acceptable", or "elongated"
  };
  overall_pass: boolean;
  error?: string;
  // Hotspot information
  hotspots_used?: string[];
  hotspots_auto_detected?: boolean;
  // Interaction analysis (from PLIP)
  interactions?: InteractionProfile;
  key_binding_residues?: string[];
  recommendations?: string[];
}

interface BinderEvaluationCardProps {
  evaluation: BinderEvaluation;
  expanded?: boolean;
  onViewDetails?: () => void;
}

export function BinderEvaluationCard({ evaluation, expanded = false, onViewDetails }: BinderEvaluationCardProps) {
  const stats = evaluation.statistics;
  const best = evaluation.best_design;

  // Helper function to get quality indicator
  const getQualityBadge = (value: number | undefined, thresholds: { good: number; medium: number }, higherIsBetter: boolean = true) => {
    if (value === undefined) return null;

    const isGood = higherIsBetter ? value >= thresholds.good : value <= thresholds.good;
    const isMedium = higherIsBetter ? value >= thresholds.medium : value <= thresholds.medium;

    if (isGood) {
      return <span className="text-xs px-2 py-0.5 bg-emerald-100 text-emerald-700 rounded-full">Excellent</span>;
    } else if (isMedium) {
      return <span className="text-xs px-2 py-0.5 bg-amber-100 text-amber-700 rounded-full">Good</span>;
    } else {
      return <span className="text-xs px-2 py-0.5 bg-red-100 text-red-700 rounded-full">Low</span>;
    }
  };

  return (
    <div className="bg-white rounded-xl border border-slate-200 overflow-hidden shadow-sm">
      {/* Header */}
      <div className={`px-4 py-3 border-b ${evaluation.overall_pass ? 'bg-emerald-50 border-emerald-100' : 'bg-red-50 border-red-100'}`}>
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-2">
            <span className={`material-symbols-outlined ${evaluation.overall_pass ? 'text-emerald-600' : 'text-red-600'}`}>
              {evaluation.overall_pass ? 'check_circle' : 'cancel'}
            </span>
            <h4 className="font-semibold text-slate-900 text-sm">
              {evaluation.overall_pass ? 'Binder Design Successful' : 'Design Needs Improvement'}
            </h4>
          </div>
          {evaluation.overall_pass && stats && (
            <span className="text-xs bg-emerald-100 text-emerald-700 px-2 py-0.5 rounded-full">
              {stats.returned}/{stats.generated} passed
            </span>
          )}
        </div>
      </div>

      {/* Pipeline Funnel Summary */}
      {stats && (
        <div className="px-4 py-4 border-b border-slate-100">
          <h5 className="text-xs font-medium text-slate-500 uppercase tracking-wider mb-3">Pipeline Results</h5>
          <div className="grid grid-cols-4 gap-2 text-center">
            <div className="p-2 bg-slate-50 rounded-lg">
              <div className="text-lg font-bold text-slate-900">{stats.generated}</div>
              <div className="text-xs text-slate-500">Generated</div>
            </div>
            <div className="p-2 bg-blue-50 rounded-lg">
              <div className="text-lg font-bold text-blue-600">{stats.mpnn_designed}</div>
              <div className="text-xs text-slate-500">MPNN</div>
            </div>
            <div className="p-2 bg-violet-50 rounded-lg">
              <div className="text-lg font-bold text-violet-600">{stats.esm_passed}</div>
              <div className="text-xs text-slate-500">ESM Pass</div>
            </div>
            <div className="p-2 bg-emerald-50 rounded-lg">
              <div className="text-lg font-bold text-emerald-600">{stats.returned}</div>
              <div className="text-xs text-slate-500">Returned</div>
            </div>
          </div>
        </div>
      )}

      {/* Best Design Metrics */}
      {best && (
        <div className="px-4 py-4">
          <h5 className="text-xs font-medium text-slate-500 uppercase tracking-wider mb-3">
            Best Design (#{best.rank})
          </h5>

          <div className="space-y-3">
            {/* ESM Metrics */}
            <div className="flex items-center justify-between py-2 border-b border-slate-100">
              <div className="flex items-center gap-2">
                <span className="material-symbols-outlined text-violet-500 text-lg">psychology</span>
                <span className="text-sm text-slate-600">ESM Confidence</span>
              </div>
              <div className="flex items-center gap-2">
                <span className="font-medium text-slate-900">
                  {best.esm_confidence != null ? `${(best.esm_confidence * 100).toFixed(0)}%` : 'N/A'}
                </span>
                {getQualityBadge(best.esm_confidence, { good: 0.8, medium: 0.6 })}
              </div>
            </div>

            {best.esm_perplexity != null && (
              <div className="flex items-center justify-between py-2 border-b border-slate-100">
                <div className="flex items-center gap-2">
                  <span className="material-symbols-outlined text-violet-500 text-lg">analytics</span>
                  <span className="text-sm text-slate-600">ESM Perplexity</span>
                </div>
                <div className="flex items-center gap-2">
                  <span className="font-medium text-slate-900">{best.esm_perplexity.toFixed(2)}</span>
                  {getQualityBadge(best.esm_perplexity, { good: 5.0, medium: 8.0 }, false)}
                </div>
              </div>
            )}

            {/* Interface Metrics */}
            {best.interface_contacts != null && (
              <div className="flex items-center justify-between py-2 border-b border-slate-100">
                <div className="flex items-center gap-2">
                  <span className="material-symbols-outlined text-teal-500 text-lg">link</span>
                  <span className="text-sm text-slate-600">Interface Contacts</span>
                </div>
                <div className="flex items-center gap-2">
                  <span className="font-medium text-slate-900">{best.interface_contacts}</span>
                  {getQualityBadge(best.interface_contacts, { good: 5000, medium: 2000 })}
                </div>
              </div>
            )}

            {best.interface_hbonds != null && (
              <div className="flex items-center justify-between py-2 border-b border-slate-100">
                <div className="flex items-center gap-2">
                  <span className="material-symbols-outlined text-blue-500 text-lg">humidity_percentage</span>
                  <span className="text-sm text-slate-600">H-Bonds</span>
                </div>
                <div className="flex items-center gap-2">
                  <span className="font-medium text-slate-900">{best.interface_hbonds}</span>
                  {getQualityBadge(best.interface_hbonds, { good: 10, medium: 5 })}
                </div>
              </div>
            )}

            {best.buried_sasa != null && (
              <div className="flex items-center justify-between py-2 border-b border-slate-100">
                <div className="flex items-center gap-2">
                  <span className="material-symbols-outlined text-amber-500 text-lg">layers</span>
                  <span className="text-sm text-slate-600">Buried SASA</span>
                </div>
                <span className="font-medium text-slate-900">{best.buried_sasa.toFixed(0)} Å²</span>
              </div>
            )}

            {best.packstat != null && (
              <div className="flex items-center justify-between py-2 border-b border-slate-100">
                <div className="flex items-center gap-2">
                  <span className="material-symbols-outlined text-emerald-500 text-lg">check_box</span>
                  <span className="text-sm text-slate-600">Packstat</span>
                </div>
                <div className="flex items-center gap-2">
                  <span className="font-medium text-slate-900">{best.packstat.toFixed(3)}</span>
                  {getQualityBadge(best.packstat, { good: 0.6, medium: 0.5 })}
                </div>
              </div>
            )}

            {best.rg_ratio != null && (
              <div className="flex items-center justify-between py-2">
                <div className="flex items-center gap-2">
                  <span className="material-symbols-outlined text-indigo-500 text-lg">circle</span>
                  <span className="text-sm text-slate-600">Compactness (Rg)</span>
                </div>
                <div className="flex items-center gap-2">
                  <span className="font-medium text-slate-900">{best.rg_ratio.toFixed(2)}</span>
                  {getQualityBadge(best.rg_ratio, { good: 1.2, medium: 1.5 }, false)}
                </div>
              </div>
            )}
          </div>

          {/* Interaction Profile (from PLIP analysis) */}
          {evaluation.interactions && expanded && (
            <div className="mt-4 p-3 bg-gradient-to-br from-blue-50 to-purple-50 rounded-lg border border-blue-100">
              <div className="flex items-center gap-2 mb-3">
                <span className="material-symbols-outlined text-blue-600 text-sm">science</span>
                <span className="text-xs font-medium text-blue-700 uppercase">Interaction Profile</span>
              </div>
              <div className="grid grid-cols-5 gap-2 text-center">
                <div className="p-2 bg-white/60 rounded-lg">
                  <div className="text-lg font-bold text-blue-600">{evaluation.interactions.hbonds}</div>
                  <div className="text-xs text-slate-500">H-bonds</div>
                </div>
                <div className="p-2 bg-white/60 rounded-lg">
                  <div className="text-lg font-bold text-yellow-600">{evaluation.interactions.hydrophobic}</div>
                  <div className="text-xs text-slate-500">Hydrophobic</div>
                </div>
                <div className="p-2 bg-white/60 rounded-lg">
                  <div className="text-lg font-bold text-purple-600">{evaluation.interactions.pi_stacking}</div>
                  <div className="text-xs text-slate-500">π-Stack</div>
                </div>
                <div className="p-2 bg-white/60 rounded-lg">
                  <div className="text-lg font-bold text-red-600">{evaluation.interactions.salt_bridges}</div>
                  <div className="text-xs text-slate-500">Salt Bridge</div>
                </div>
                <div className="p-2 bg-white/60 rounded-lg">
                  <div className="text-lg font-bold text-slate-800">{evaluation.interactions.total}</div>
                  <div className="text-xs text-slate-500">Total</div>
                </div>
              </div>
              {/* Key Binding Residues */}
              {evaluation.key_binding_residues && evaluation.key_binding_residues.length > 0 && (
                <div className="mt-3 pt-3 border-t border-blue-100">
                  <div className="text-xs text-slate-500 mb-1">Key Binding Residues</div>
                  <div className="flex flex-wrap gap-1">
                    {evaluation.key_binding_residues.slice(0, 8).map((residue) => (
                      <span key={residue} className="text-xs px-2 py-0.5 bg-emerald-100 text-emerald-700 rounded">
                        {residue}
                      </span>
                    ))}
                    {evaluation.key_binding_residues.length > 8 && (
                      <span className="text-xs px-2 py-0.5 text-slate-400">
                        +{evaluation.key_binding_residues.length - 8} more
                      </span>
                    )}
                  </div>
                </div>
              )}
            </div>
          )}

          {/* AI Recommendations */}
          {evaluation.recommendations && evaluation.recommendations.length > 0 && expanded && (
            <div className="mt-4 p-3 bg-amber-50 rounded-lg border border-amber-100">
              <div className="flex items-center gap-2 mb-2">
                <span className="material-symbols-outlined text-amber-600 text-sm">lightbulb</span>
                <span className="text-xs font-medium text-amber-700 uppercase">Suggestions for Improvement</span>
              </div>
              <ul className="space-y-1">
                {evaluation.recommendations.map((rec, idx) => (
                  <li key={idx} className="text-sm text-amber-800 flex items-start gap-2">
                    <span className="text-amber-500 mt-0.5">•</span>
                    {rec}
                  </li>
                ))}
              </ul>
            </div>
          )}

          {/* Hotspot Info */}
          {evaluation.hotspots_used && evaluation.hotspots_used.length > 0 && expanded && (
            <div className="mt-4 p-3 bg-indigo-50 rounded-lg">
              <div className="flex items-center gap-2 mb-2">
                <span className="material-symbols-outlined text-indigo-500 text-sm">my_location</span>
                <span className="text-xs font-medium text-indigo-700 uppercase">
                  {evaluation.hotspots_auto_detected ? 'Auto-Detected Hotspots' : 'Manual Hotspots'}
                </span>
              </div>
              <div className="flex flex-wrap gap-1">
                {evaluation.hotspots_used.map((hotspot) => (
                  <span key={hotspot} className="text-xs px-2 py-0.5 bg-indigo-100 text-indigo-700 rounded">
                    {hotspot}
                  </span>
                ))}
              </div>
            </div>
          )}

          {/* Sequence */}
          {expanded && best.binder_sequence && (
            <div className="mt-4 p-3 bg-slate-50 rounded-lg">
              <div className="flex items-center justify-between mb-2">
                <span className="text-xs font-medium text-slate-500 uppercase">Binder Sequence</span>
                <button
                  onClick={() => navigator.clipboard.writeText(best.binder_sequence)}
                  className="text-xs text-teal-600 hover:text-teal-700 flex items-center gap-1"
                >
                  <span className="material-symbols-outlined text-sm">content_copy</span>
                  Copy
                </button>
              </div>
              <code className="text-xs text-slate-700 font-mono break-all">
                {best.binder_sequence}
              </code>
            </div>
          )}
        </div>
      )}

      {/* Error State */}
      {evaluation.error && (
        <div className="px-4 py-4 bg-red-50">
          <div className="flex items-start gap-2">
            <span className="material-symbols-outlined text-red-500">error</span>
            <div>
              <p className="text-sm font-medium text-red-800">Design Failed</p>
              <p className="text-xs text-red-600 mt-1">{evaluation.error}</p>
            </div>
          </div>
        </div>
      )}

      {/* View Details Button */}
      {!expanded && onViewDetails && (
        <div className="px-4 py-3 bg-slate-50 border-t border-slate-100">
          <button
            onClick={onViewDetails}
            className="w-full text-sm text-teal-600 hover:text-teal-700 font-medium flex items-center justify-center gap-1"
          >
            View Full Results
            <span className="material-symbols-outlined text-sm">arrow_forward</span>
          </button>
        </div>
      )}
    </div>
  );
}
