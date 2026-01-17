'use client';

import { useState } from 'react';
import {
  ClipboardCheck,
  Zap,
  Link2,
  Fingerprint,
  Shield,
  Droplets,
  Network,
  Brain,
  FlaskConical,
  Box,
  HelpCircle,
  CheckCircle2,
  Info,
  type LucideIcon,
} from 'lucide-react';

export type FilterMode = 'heterodimer' | 'rfdiffusion3';

// Heterodimer-specific metrics (custom - default)
export interface HeterodimerMetrics {
  affinity?: number;           // GNINA: < -5 = good (kcal/mol)
  contacts_a?: number;         // ≥ 5 = good
  contacts_b?: number;         // ≥ 5 = good
  sequence_identity?: number;  // < 70% = true heterodimer
  anti_homo_score?: number;    // > 60 = good (0-100 scale)
  n7_hbonds?: number;          // ≥ 1 = azobenzene N7 (azo) satisfied
  n8_hbonds?: number;          // ≥ 1 = azobenzene N8 (azo) satisfied
  has_clashes?: boolean;       // false = good
  is_heterodimer?: boolean;    // true = success
}

// Standard RFdiffusion3 metrics (Baker Lab reference)
export interface RFD3Metrics {
  pae_interaction?: number;  // < 10 = good
  plddt?: number;            // > 80 = good, > 90 = excellent
  ddg?: number;              // < -40 = good (Rosetta kcal/mol)
  total_hbonds?: number;     // > 11 = good
  shape_complementarity?: number; // > 0.65 = excellent
}

export interface DesignMetrics extends HeterodimerMetrics, RFD3Metrics {}

interface SuccessCriteriaCardProps {
  metrics: DesignMetrics;
  filterMode?: FilterMode;
  onFilterModeChange?: (mode: FilterMode) => void;
  showToggle?: boolean;
}

// Thresholds for both modes
const THRESHOLDS = {
  heterodimer: {
    affinity: { excellent: -7, good: -5, poor: -2, unit: 'kcal/mol', lowerBetter: true },
    contacts_a: { excellent: 10, good: 5, poor: 2, unit: '', lowerBetter: false },
    contacts_b: { excellent: 10, good: 5, poor: 2, unit: '', lowerBetter: false },
    sequence_identity: { excellent: 40, good: 70, poor: 90, unit: '%', lowerBetter: true },
    anti_homo_score: { excellent: 80, good: 60, poor: 40, unit: '/100', lowerBetter: false },
    n7_hbonds: { excellent: 2, good: 1, poor: 0, unit: '', lowerBetter: false },
    n8_hbonds: { excellent: 2, good: 1, poor: 0, unit: '', lowerBetter: false },
  },
  rfdiffusion3: {
    pae_interaction: { excellent: 5, good: 10, poor: 15, unit: '', lowerBetter: true },
    plddt: { excellent: 90, good: 80, poor: 70, unit: '', lowerBetter: false },
    ddg: { excellent: -50, good: -40, poor: -20, unit: 'kcal/mol', lowerBetter: true },
    total_hbonds: { excellent: 15, good: 11, poor: 5, unit: '', lowerBetter: false },
    shape_complementarity: { excellent: 0.75, good: 0.65, poor: 0.5, unit: '', lowerBetter: false },
  },
};

const METRIC_LABELS: Record<string, { name: string; Icon: LucideIcon }> = {
  // Heterodimer
  affinity: { name: 'GNINA Affinity', Icon: Zap },
  contacts_a: { name: 'Chain A Contacts', Icon: Link2 },
  contacts_b: { name: 'Chain B Contacts', Icon: Link2 },
  sequence_identity: { name: 'Sequence Identity', Icon: Fingerprint },
  anti_homo_score: { name: 'Anti-Homo Score', Icon: Shield },
  n7_hbonds: { name: 'N7 H-bonds', Icon: Droplets },
  n8_hbonds: { name: 'N8 H-bonds', Icon: Droplets },
  // RFdiffusion3
  pae_interaction: { name: 'pAE (Interface)', Icon: Network },
  plddt: { name: 'pLDDT', Icon: Brain },
  ddg: { name: 'Rosetta ddG', Icon: FlaskConical },
  total_hbonds: { name: 'Total H-bonds', Icon: Droplets },
  shape_complementarity: { name: 'Shape Comp.', Icon: Box },
  // Fallback for unknown metrics
  help: { name: 'Unknown', Icon: HelpCircle },
};

type MetricStatus = 'excellent' | 'good' | 'moderate' | 'poor' | 'unknown';

// Type for threshold values
type ThresholdConfig = {
  excellent: number;
  good: number;
  poor: number;
  unit: string;
  lowerBetter: boolean;
};

function getMetricStatus(
  metricKey: string,
  value: number | undefined,
  mode: FilterMode
): MetricStatus {
  if (value === undefined || value === null) return 'unknown';

  const modeThresholds = THRESHOLDS[mode] as Record<string, ThresholdConfig>;
  const thresholds = modeThresholds?.[metricKey];
  if (!thresholds) return 'unknown';

  const { excellent, good, poor, lowerBetter } = thresholds;

  if (lowerBetter) {
    if (value <= excellent) return 'excellent';
    if (value <= good) return 'good';
    if (value <= poor) return 'moderate';
    return 'poor';
  } else {
    if (value >= excellent) return 'excellent';
    if (value >= good) return 'good';
    if (value >= poor) return 'moderate';
    return 'poor';
  }
}

const STATUS_STYLES: Record<MetricStatus, { bg: string; text: string; badge: string }> = {
  excellent: { bg: 'bg-emerald-50', text: 'text-emerald-700', badge: 'bg-emerald-100 text-emerald-700' },
  good: { bg: 'bg-blue-50', text: 'text-blue-700', badge: 'bg-blue-100 text-blue-700' },
  moderate: { bg: 'bg-amber-50', text: 'text-amber-700', badge: 'bg-amber-100 text-amber-700' },
  poor: { bg: 'bg-red-50', text: 'text-red-700', badge: 'bg-red-100 text-red-700' },
  unknown: { bg: 'bg-slate-50', text: 'text-slate-500', badge: 'bg-slate-100 text-slate-500' },
};

export function SuccessCriteriaCard({
  metrics,
  filterMode: externalMode,
  onFilterModeChange,
  showToggle = true,
}: SuccessCriteriaCardProps) {
  const [internalMode, setInternalMode] = useState<FilterMode>('heterodimer');
  const filterMode = externalMode ?? internalMode;

  const handleModeChange = (mode: FilterMode) => {
    setInternalMode(mode);
    onFilterModeChange?.(mode);
  };

  const renderMetric = (key: string, value: number | boolean | undefined) => {
    if (typeof value === 'boolean') {
      const status: MetricStatus = key === 'has_clashes' ? (value ? 'poor' : 'excellent') : (value ? 'excellent' : 'poor');
      const styles = STATUS_STYLES[status];
      const label = METRIC_LABELS[key] || METRIC_LABELS.help;
      const IconComponent = label.Icon;

      return (
        <div key={key} className={`p-3 rounded-lg ${styles.bg}`}>
          <div className="flex items-center gap-2 mb-1">
            <IconComponent className={`h-4 w-4 ${styles.text}`} />
            <span className="text-xs text-slate-500">{label.name}</span>
          </div>
          <div className="flex items-center justify-between">
            <span className={`font-semibold ${styles.text}`}>
              {key === 'has_clashes' ? (value ? 'Yes' : 'None') : (value ? 'Yes' : 'No')}
            </span>
            <span className={`text-xs px-2 py-0.5 rounded-full ${styles.badge}`}>
              {status === 'excellent' ? 'Good' : status === 'poor' ? 'Issue' : status}
            </span>
          </div>
        </div>
      );
    }

    const numValue = value as number | undefined;
    const status = getMetricStatus(key, numValue, filterMode);
    const styles = STATUS_STYLES[status];
    const label = METRIC_LABELS[key] || METRIC_LABELS.help;
    const IconComponent = label.Icon;
    // Get threshold for this metric
    const modeThresholds = THRESHOLDS[filterMode] as Record<string, ThresholdConfig>;
    const threshold = modeThresholds?.[key];
    const unit = threshold?.unit || '';

    return (
      <div key={key} className={`p-3 rounded-lg ${styles.bg}`}>
        <div className="flex items-center gap-2 mb-1">
          <IconComponent className={`h-4 w-4 ${styles.text}`} />
          <span className="text-xs text-slate-500">{label.name}</span>
        </div>
        <div className="flex items-center justify-between">
          <span className={`font-semibold ${styles.text}`}>
            {numValue !== undefined ? `${numValue.toFixed(2)}${unit}` : 'N/A'}
          </span>
          <span className={`text-xs px-2 py-0.5 rounded-full ${styles.badge}`}>
            {status === 'unknown' ? 'N/A' : status.charAt(0).toUpperCase() + status.slice(1)}
          </span>
        </div>
      </div>
    );
  };

  // Get metrics to display based on mode
  const heterodimerKeys: (keyof HeterodimerMetrics)[] = [
    'affinity', 'contacts_a', 'contacts_b', 'sequence_identity',
    'anti_homo_score', 'n7_hbonds', 'n8_hbonds'
  ];
  const rfd3Keys: (keyof RFD3Metrics)[] = [
    'pae_interaction', 'plddt', 'ddg', 'total_hbonds', 'shape_complementarity'
  ];

  const displayKeys = filterMode === 'heterodimer' ? heterodimerKeys : rfd3Keys;

  // Calculate overall status
  const calculateOverallStatus = (): { passed: number; total: number; status: MetricStatus } => {
    let passed = 0;
    let total = 0;

    displayKeys.forEach((key) => {
      const value = metrics[key];
      if (value !== undefined && typeof value === 'number') {
        total++;
        const status = getMetricStatus(key, value, filterMode);
        if (status === 'excellent' || status === 'good') passed++;
      }
    });

    // Add boolean checks for heterodimer mode
    if (filterMode === 'heterodimer') {
      if (metrics.has_clashes !== undefined) {
        total++;
        if (!metrics.has_clashes) passed++;
      }
      if (metrics.is_heterodimer !== undefined) {
        total++;
        if (metrics.is_heterodimer) passed++;
      }
    }

    const ratio = total > 0 ? passed / total : 0;
    let status: MetricStatus = 'poor';
    if (ratio >= 0.8) status = 'excellent';
    else if (ratio >= 0.6) status = 'good';
    else if (ratio >= 0.4) status = 'moderate';

    return { passed, total, status };
  };

  const overall = calculateOverallStatus();
  const overallStyles = STATUS_STYLES[overall.status];

  return (
    <div className="bg-white rounded-xl border border-slate-200 overflow-hidden shadow-sm">
      {/* Header with Toggle */}
      <div className="px-4 py-3 border-b border-slate-100 bg-slate-50">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-2">
            <ClipboardCheck className="h-5 w-5 text-slate-600" />
            <h4 className="font-semibold text-slate-900 text-sm">Success Criteria</h4>
          </div>

          {showToggle && (
            <div className="flex items-center gap-1 bg-white rounded-lg p-1 border border-slate-200">
              <button
                onClick={() => handleModeChange('heterodimer')}
                className={`px-3 py-1 text-xs font-medium rounded-md transition-all ${
                  filterMode === 'heterodimer'
                    ? 'bg-purple-100 text-purple-700'
                    : 'text-slate-500 hover:text-slate-700'
                }`}
              >
                Heterodimer
              </button>
              <button
                onClick={() => handleModeChange('rfdiffusion3')}
                className={`px-3 py-1 text-xs font-medium rounded-md transition-all ${
                  filterMode === 'rfdiffusion3'
                    ? 'bg-blue-100 text-blue-700'
                    : 'text-slate-500 hover:text-slate-700'
                }`}
              >
                RFD3 (Baker)
              </button>
            </div>
          )}
        </div>
      </div>

      {/* Overall Status */}
      <div className={`px-4 py-3 ${overallStyles.bg} border-b border-slate-100`}>
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-2">
            {overall.status === 'excellent' || overall.status === 'good' ? (
              <CheckCircle2 className={`h-5 w-5 ${overallStyles.text}`} />
            ) : (
              <Info className={`h-5 w-5 ${overallStyles.text}`} />
            )}
            <span className={`font-medium ${overallStyles.text}`}>
              {overall.passed}/{overall.total} criteria passed
            </span>
          </div>
          <span className={`text-xs px-2 py-0.5 rounded-full ${overallStyles.badge}`}>
            {overall.status.charAt(0).toUpperCase() + overall.status.slice(1)}
          </span>
        </div>
      </div>

      {/* Metrics Grid */}
      <div className="p-4">
        <div className="grid grid-cols-2 sm:grid-cols-3 lg:grid-cols-4 gap-3">
          {displayKeys.map((key) => renderMetric(key, metrics[key]))}

          {/* Boolean metrics for heterodimer mode */}
          {filterMode === 'heterodimer' && (
            <>
              {metrics.has_clashes !== undefined && renderMetric('has_clashes', metrics.has_clashes)}
              {metrics.is_heterodimer !== undefined && renderMetric('is_heterodimer', metrics.is_heterodimer)}
            </>
          )}
        </div>
      </div>

      {/* Mode Description */}
      <div className="px-4 py-2 bg-slate-50 border-t border-slate-100">
        <p className="text-xs text-slate-500">
          {filterMode === 'heterodimer' ? (
            <>
              <span className="font-medium">Heterodimer mode:</span> Custom thresholds for GNINA docking,
              anti-homodimerization validation, and azobenzene-specific H-bond tracking.
            </>
          ) : (
            <>
              <span className="font-medium">RFdiffusion3 mode:</span> Standard Baker Lab thresholds using
              AlphaFold2 pAE, pLDDT, and Rosetta ddG metrics.
            </>
          )}
        </p>
      </div>
    </div>
  );
}

export default SuccessCriteriaCard;
