'use client';

import { Network, Ruler, Shapes, CircleDot, CheckCircle, AlertCircle, ThumbsUp, Lightbulb, RefreshCw, BadgeCheck, Clock, type LucideIcon } from 'lucide-react';
import type { DesignEvaluation } from '@/lib/api';

// Icon mapping for criteria icons
const CRITERIA_ICONS: Record<string, LucideIcon> = {
  hub: Network,
  straighten: Ruler,
  category: Shapes,
  scatter_plot: CircleDot,
};

interface EvaluationCriteria {
  label: string;
  value: string;
  target: string;
  pass: boolean;
  icon: string;
}

interface EvaluationCardProps {
  evaluation: DesignEvaluation;
  targetMetal: string;
  onRetry?: () => void;
}

export function EvaluationCard({ evaluation, targetMetal, onRetry }: EvaluationCardProps) {
  // Define target ranges based on metal type
  const isLanthanide = ['TB', 'GD', 'EU', 'LA', 'CE', 'ND'].includes(targetMetal.toUpperCase());
  const targetCoordRange = isLanthanide ? [8, 9] : [4, 6];
  const targetDistRange = isLanthanide ? [2.3, 2.5] : [2.0, 2.4];

  const criteria: EvaluationCriteria[] = [
    {
      label: 'Coordination Number',
      value: String(evaluation.coordination_number),
      target: `${targetCoordRange[0]}-${targetCoordRange[1]}`,
      pass: evaluation.coordination_number >= targetCoordRange[0] && evaluation.coordination_number <= targetCoordRange[1],
      icon: 'hub',
    },
    {
      label: 'Bond Distance',
      value: `${evaluation.avg_bond_distance.toFixed(2)}A`,
      target: `${targetDistRange[0]}-${targetDistRange[1]}A`,
      pass: evaluation.avg_bond_distance >= targetDistRange[0] && evaluation.avg_bond_distance <= targetDistRange[1],
      icon: 'straighten',
    },
    {
      label: 'Geometry',
      value: evaluation.geometry_type,
      target: isLanthanide ? 'Square antiprism' : 'Tetrahedral/Octahedral',
      pass: evaluation.geometry_rmsd < 0.5,
      icon: 'category',
    },
    {
      label: 'Donor Types',
      value: isLanthanide
        ? `${evaluation.oxygen_donors} O-donors`
        : `${evaluation.oxygen_donors}O, ${evaluation.nitrogen_donors}N, ${evaluation.sulfur_donors}S`,
      target: isLanthanide ? '6+ oxygen' : 'Appropriate donors',
      pass: isLanthanide ? evaluation.oxygen_donors >= 6 : true,
      icon: 'scatter_plot',
    },
  ];

  const passCount = criteria.filter(c => c.pass).length;
  const overallPass = passCount >= 3;

  return (
    <div className={`rounded-xl p-5 border shadow-sm ${
      overallPass
        ? 'bg-success/10 border-success/20'
        : 'bg-warning/10 border-warning/20'
    }`}>
      {/* Header with score */}
      <div className="flex items-center justify-between mb-4">
        <h4 className="font-semibold text-foreground flex items-center gap-2">
          {overallPass
            ? <BadgeCheck className="h-5 w-5 text-success" />
            : <Clock className="h-5 w-5 text-warning" />
          }
          Design Evaluation
        </h4>
        <span className={`px-3 py-1 rounded-full text-xs font-bold ${
          overallPass
            ? 'bg-success text-success-foreground'
            : 'bg-warning text-warning-foreground'
        }`}>
          {passCount}/{criteria.length} PASS
        </span>
      </div>

      {/* Criteria checklist */}
      <div className="space-y-2 mb-4">
        {criteria.map((c, i) => (
          <div
            key={i}
            className={`flex items-center justify-between p-3 rounded-lg ${
              c.pass ? 'bg-white/60' : 'bg-white/80'
            }`}
          >
            <div className="flex items-center gap-3">
              {(() => {
                const IconComponent = CRITERIA_ICONS[c.icon];
                return IconComponent ? <IconComponent className={`h-4 w-4 ${c.pass ? 'text-success' : 'text-warning'}`} /> : null;
              })()}
              <div>
                <div className="text-sm font-medium text-foreground">{c.label}</div>
                <div className="text-xs text-muted-foreground">Target: {c.target}</div>
              </div>
            </div>
            <div className="flex items-center gap-2">
              <span className="font-mono text-sm text-foreground">{c.value}</span>
              {c.pass
                ? <CheckCircle className="h-5 w-5 text-success" />
                : <AlertCircle className="h-5 w-5 text-warning" />
              }
            </div>
          </div>
        ))}
      </div>

      {/* Overall assessment */}
      <div className={`p-3 rounded-lg text-sm ${
        overallPass
          ? 'bg-success/20 text-success'
          : 'bg-warning/20 text-warning'
      }`}>
        <div className="flex items-start gap-2">
          {overallPass ? <ThumbsUp className="h-4 w-4 mt-0.5" /> : <Lightbulb className="h-4 w-4 mt-0.5" />}
          <div>
            {overallPass ? (
              <p>
                <strong>Excellent results!</strong> This design meets the target criteria for {targetMetal} binding.
                Consider proceeding with sequence design (MPNN) and validation (RF3).
              </p>
            ) : (
              <p>
                <strong>Room for improvement.</strong> The design shows promise but may benefit from
                adjustments. Consider using a more aggressive design approach or modifying the coordination preferences.
              </p>
            )}
          </div>
        </div>
      </div>

      {/* Suggestions if not all pass */}
      {!overallPass && evaluation.suggestions && evaluation.suggestions.length > 0 && (
        <div className="mt-3 p-3 bg-white/60 rounded-lg">
          <div className="text-xs font-semibold text-muted-foreground mb-2 flex items-center gap-1">
            <Lightbulb className="h-4 w-4" />
            Suggestions
          </div>
          <ul className="text-xs text-muted-foreground space-y-1">
            {evaluation.suggestions.map((s, i) => (
              <li key={i} className="flex items-start gap-2">
                <span className="text-warning">•</span>
                <span>{s}</span>
              </li>
            ))}
          </ul>
        </div>
      )}

      {/* Retry button */}
      {onRetry && !overallPass && (
        <button
          onClick={onRetry}
          className="mt-4 w-full py-2.5 bg-white border border-warning/30 text-warning rounded-lg font-medium hover:bg-warning/10 transition-all flex items-center justify-center gap-2"
        >
          <RefreshCw className="h-4 w-4" />
          Try Different Settings
        </button>
      )}
    </div>
  );
}
