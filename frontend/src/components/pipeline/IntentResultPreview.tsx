'use client';

import { Badge } from '@/components/ui/badge';
import {
  AlertTriangle,
  Lightbulb,
  Brain,
  Atom,
  FlaskConical,
  Target,
  Ruler,
  Zap,
} from 'lucide-react';
import { cn } from '@/lib/utils';
import type { StepResult } from '@/lib/pipeline-types';

interface IntentResultPreviewProps {
  result: StepResult;
  onSelectDesign?: (pdbContent: string) => void;
}

export function IntentResultPreview({ result }: IntentResultPreviewProps) {
  const data = result.data || {};
  const aiParsed = data.ai_parsed as boolean;
  const confidence = (data.confidence as number) || 0;
  const parserType = data.parser_type as string;
  const designType = data.design_type as string || 'general';
  const targetMetal = data.target_metal as string;
  const ligandName = data.ligand_name as string;
  const designGoal = data.design_goal as string;
  const topology = data.target_topology as string;
  const enzymeClass = data.enzyme_class as string;
  const preserveFunction = data.preserve_function as boolean;
  const chainMin = data.chain_length_min as number;
  const chainMax = data.chain_length_max as number;
  const designMode = data.design_mode as string;
  const pdbId = data.pdb_id as string;
  const warnings = (data.warnings as string[]) || [];
  const suggestions = (data.suggestions as string[]) || [];
  const typoCorrections = (data.typo_corrections as string[]) || [];

  const confidenceColor = confidence >= 0.7
    ? 'bg-green-100 text-green-800 dark:bg-green-900 dark:text-green-200'
    : confidence >= 0.4
      ? 'bg-yellow-100 text-yellow-800 dark:bg-yellow-900 dark:text-yellow-200'
      : 'bg-red-100 text-red-800 dark:bg-red-900 dark:text-red-200';

  const designTypeLabel: Record<string, string> = {
    metal: 'Metal Binding',
    binder: 'Protein Binder',
    ligand: 'Ligand-Mediated Dimer',
    scaffold: 'Scaffold Design',
    general: 'General Design',
  };

  return (
    <div className="space-y-3">
      {/* Parser type + confidence */}
      <div className="flex items-center gap-2 flex-wrap">
        {aiParsed && (
          <Badge variant="outline" className="gap-1">
            <Brain className="h-3 w-3" />
            {parserType === 'ai' ? 'AI Parsed' : 'Server Parsed'}
          </Badge>
        )}
        {!aiParsed && (
          <Badge variant="outline" className="gap-1 bg-orange-50 dark:bg-orange-950">
            <Zap className="h-3 w-3" />
            Local Fallback
          </Badge>
        )}
        {aiParsed && confidence > 0 && (
          <Badge className={cn('gap-1', confidenceColor)}>
            Confidence: {(confidence * 100).toFixed(0)}%
          </Badge>
        )}
        <Badge variant="secondary">
          {designTypeLabel[designType] || designType}
        </Badge>
      </div>

      {/* Detected entities */}
      <div className="grid grid-cols-2 gap-2 text-sm">
        {targetMetal && (
          <div className="flex items-center gap-1.5">
            <Atom className="h-3.5 w-3.5 text-blue-500" />
            <span className="text-muted-foreground">Metal:</span>
            <span className="font-medium">{targetMetal}</span>
          </div>
        )}
        {ligandName && (
          <div className="flex items-center gap-1.5">
            <FlaskConical className="h-3.5 w-3.5 text-purple-500" />
            <span className="text-muted-foreground">Ligand:</span>
            <span className="font-medium">{ligandName}</span>
          </div>
        )}
        {designGoal && designGoal !== 'binding' && (
          <div className="flex items-center gap-1.5">
            <Target className="h-3.5 w-3.5 text-green-500" />
            <span className="text-muted-foreground">Goal:</span>
            <span className="font-medium">{designGoal}</span>
          </div>
        )}
        {topology && topology !== 'monomer' && (
          <div className="flex items-center gap-1.5">
            <span className="text-muted-foreground">Topology:</span>
            <span className="font-medium">{topology}</span>
          </div>
        )}
        {enzymeClass && (
          <div className="flex items-center gap-1.5">
            <FlaskConical className="h-3.5 w-3.5 text-amber-500" />
            <span className="text-muted-foreground">Enzyme:</span>
            <span className="font-medium">{enzymeClass}</span>
          </div>
        )}
        {chainMin && chainMax && (
          <div className="flex items-center gap-1.5">
            <Ruler className="h-3.5 w-3.5 text-gray-500" />
            <span className="text-muted-foreground">Length:</span>
            <span className="font-medium">{chainMin}-{chainMax} residues</span>
          </div>
        )}
        {designMode === 'scaffold' && pdbId && pdbId !== 'Not specified' && (
          <div className="flex items-center gap-1.5">
            <span className="text-muted-foreground">Scaffold:</span>
            <span className="font-medium">{pdbId}</span>
          </div>
        )}
        {preserveFunction && (
          <div className="flex items-center gap-1.5 text-amber-600 dark:text-amber-400">
            <span className="font-medium">Preserve enzymatic function</span>
          </div>
        )}
      </div>

      {/* Typo corrections */}
      {typoCorrections.length > 0 && (
        <div className="text-xs text-muted-foreground italic">
          Auto-corrected: {typoCorrections.join(', ')}
        </div>
      )}

      {/* Warnings */}
      {warnings.length > 0 && (
        <div className="space-y-1">
          {warnings.map((w, i) => (
            <div key={i} className="flex items-start gap-1.5 text-xs text-yellow-700 dark:text-yellow-400">
              <AlertTriangle className="h-3.5 w-3.5 mt-0.5 shrink-0" />
              <span>{w}</span>
            </div>
          ))}
        </div>
      )}

      {/* Suggestions */}
      {suggestions.length > 0 && (
        <div className="space-y-1">
          {suggestions.map((s, i) => (
            <div key={i} className="flex items-start gap-1.5 text-xs text-blue-700 dark:text-blue-400">
              <Lightbulb className="h-3.5 w-3.5 mt-0.5 shrink-0" />
              <span>{s}</span>
            </div>
          ))}
        </div>
      )}
    </div>
  );
}
