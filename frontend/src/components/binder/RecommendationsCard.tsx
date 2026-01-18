'use client';

import {
  Lightbulb,
  Droplets,
  Cylinder,
  Hexagon,
  Zap,
  type LucideIcon,
} from 'lucide-react';
import { InteractionProfile } from './InterfaceMetrics';

interface RecommendationsCardProps {
  recommendations: string[];
  interactions?: InteractionProfile;
}

// Categorize recommendations by type
const categorizeRecommendation = (rec: string): 'polar' | 'hydrophobic' | 'aromatic' | 'salt_bridge' | 'general' => {
  const recLower = rec.toLowerCase();
  if (recLower.includes('polar') || recLower.includes('gln') || recLower.includes('asn') || recLower.includes('ser') || recLower.includes('h-bond') || recLower.includes('hydrogen')) {
    return 'polar';
  }
  if (recLower.includes('hydrophobic') || recLower.includes('leu') || recLower.includes('ile') || recLower.includes('val') || recLower.includes('ala') || recLower.includes('met')) {
    return 'hydrophobic';
  }
  if (recLower.includes('aromatic') || recLower.includes('phe') || recLower.includes('tyr') || recLower.includes('trp') || recLower.includes('pi-stack') || recLower.includes('pi stack')) {
    return 'aromatic';
  }
  if (recLower.includes('salt') || recLower.includes('charged') || recLower.includes('asp') || recLower.includes('glu') || recLower.includes('lys') || recLower.includes('arg')) {
    return 'salt_bridge';
  }
  return 'general';
};

const CATEGORY_STYLES: Record<string, {
  Icon: LucideIcon;
  bg: string;
  border: string;
  text: string;
  badge: string;
}> = {
  polar: {
    Icon: Droplets,
    bg: 'bg-blue-50',
    border: 'border-blue-100',
    text: 'text-blue-700',
    badge: 'bg-blue-100 text-blue-700',
  },
  hydrophobic: {
    Icon: Cylinder,
    bg: 'bg-muted/50',
    border: 'border-border',
    text: 'text-foreground',
    badge: 'bg-muted text-foreground',
  },
  aromatic: {
    Icon: Hexagon,
    bg: 'bg-muted/50',
    border: 'border-border',
    text: 'text-foreground',
    badge: 'bg-muted text-foreground',
  },
  salt_bridge: {
    Icon: Zap,
    bg: 'bg-red-50',
    border: 'border-red-100',
    text: 'text-red-700',
    badge: 'bg-red-100 text-red-700',
  },
  general: {
    Icon: Lightbulb,
    bg: 'bg-muted/50',
    border: 'border-border',
    text: 'text-foreground',
    badge: 'bg-muted text-foreground',
  },
};

export function RecommendationsCard({ recommendations, interactions }: RecommendationsCardProps) {
  if (!recommendations || recommendations.length === 0) {
    return null;
  }

  // Categorize all recommendations
  const categorizedRecs = recommendations.map(rec => ({
    text: rec,
    category: categorizeRecommendation(rec),
  }));

  // Sort by category priority
  const categoryOrder: Array<'polar' | 'hydrophobic' | 'aromatic' | 'salt_bridge' | 'general'> = [
    'polar', 'hydrophobic', 'aromatic', 'salt_bridge', 'general'
  ];
  categorizedRecs.sort((a, b) => categoryOrder.indexOf(a.category) - categoryOrder.indexOf(b.category));

  return (
    <div className="bg-card rounded-xl border border-border overflow-hidden shadow-sm">
      {/* Header */}
      <div className="px-4 py-3 border-b border-border bg-muted">
        <div className="flex items-center gap-2">
          <Lightbulb className="h-5 w-5 text-primary" />
          <h4 className="font-semibold text-foreground text-sm">Design Recommendations</h4>
          <span className="ml-auto text-xs px-2 py-0.5 rounded-full bg-primary/10 text-primary">
            {recommendations.length} suggestions
          </span>
        </div>
        <p className="text-xs text-muted-foreground mt-1">
          AI-generated suggestions based on PLIP interaction analysis
        </p>
      </div>

      {/* Interaction Context (if available) */}
      {interactions && (
        <div className="px-4 py-2 bg-muted/50 border-b border-border">
          <div className="flex items-center gap-4 text-xs text-muted-foreground">
            <span>Current interactions:</span>
            <span className="font-medium text-foreground">{interactions.hydrogen_bonds} H-bonds</span>
            <span className="font-medium text-foreground">{interactions.hydrophobic_contacts} hydrophobic</span>
            <span className="font-medium text-foreground">{interactions.pi_stacking} pi-stack</span>
            <span className="font-medium text-foreground">{interactions.salt_bridges} salt bridges</span>
          </div>
        </div>
      )}

      {/* Recommendations List */}
      <div className="p-4 space-y-3">
        {categorizedRecs.map((rec, idx) => {
          const styles = CATEGORY_STYLES[rec.category];
          return (
            <div
              key={idx}
              className={`p-3 rounded-lg ${styles.bg} border ${styles.border}`}
            >
              <div className="flex items-start gap-3">
                <styles.Icon className={`h-5 w-5 ${styles.text} mt-0.5`} />
                <div className="flex-1">
                  <p className={`text-sm ${styles.text}`}>{rec.text}</p>
                </div>
                <span className={`text-xs px-2 py-0.5 rounded-full ${styles.badge} capitalize whitespace-nowrap`}>
                  {rec.category.replace('_', ' ')}
                </span>
              </div>
            </div>
          );
        })}
      </div>

      {/* Footer */}
      <div className="px-4 py-2 bg-muted/50 border-t border-border">
        <p className="text-xs text-muted-foreground text-center">
          Suggestions are based on interaction analysis. Consider protein stability when making modifications.
        </p>
      </div>
    </div>
  );
}

export default RecommendationsCard;
