'use client';

import { useState, useCallback } from 'react';
import { Send, Sparkles, RefreshCw, Download, Copy, Check, ChevronDown, ChevronUp } from 'lucide-react';
import { Button } from '@/components/ui/button';
import { Textarea } from '@/components/ui/textarea';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Slider } from '@/components/ui/slider';
import { Label } from '@/components/ui/label';
import { cn } from '@/lib/utils';
import { AIDesignPipelineWorkflow } from './AIDesignPipelineWorkflow';
import { useAIDesign } from '@/hooks/useAIDesign';

interface AIDesignPanelProps {
  claudeApiKey?: string;
  onDesignComplete?: (result: any) => void;
  onStructureLoad?: (pdbContent: string) => void;
}

const EXAMPLE_QUERIES = [
  "Design a protein to bind citrate with terbium",
  "Create a zinc finger DNA-binding protein",
  "Design a PQQ-binding dehydrogenase with calcium",
  "Make a luminescent biosensor with europium",
];

export function AIDesignPanel({
  claudeApiKey,
  onDesignComplete,
  onStructureLoad
}: AIDesignPanelProps) {
  const [query, setQuery] = useState('');
  const [numDesigns, setNumDesigns] = useState(4);
  const [numSequences, setNumSequences] = useState(8);
  const [showAdvanced, setShowAdvanced] = useState(false);
  const [copied, setCopied] = useState(false);

  const {
    stage,
    stageInfo,
    result,
    error,
    isRunning,
    runAIDesign,
    reset
  } = useAIDesign({ claudeApiKey });

  const handleSubmit = useCallback(async () => {
    if (!query.trim() || isRunning) return;

    const result = await runAIDesign(query, {
      numDesigns,
      numSequences,
      validate: true
    });

    if (result?.success) {
      onDesignComplete?.(result);
      if (result.best_sequence_pdb) {
        onStructureLoad?.(result.best_sequence_pdb);
      }
    }
  }, [query, numDesigns, numSequences, isRunning, runAIDesign, onDesignComplete, onStructureLoad]);

  const handleExampleClick = (example: string) => {
    setQuery(example);
  };

  const handleCopySequence = () => {
    if (result?.best_sequence?.sequence) {
      navigator.clipboard.writeText(result.best_sequence.sequence);
      setCopied(true);
      setTimeout(() => setCopied(false), 2000);
    }
  };

  const handleDownloadPDB = () => {
    if (result?.best_sequence_pdb) {
      const blob = new Blob([result.best_sequence_pdb], { type: 'chemical/x-pdb' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = 'ai_design_best.pdb';
      a.click();
      URL.revokeObjectURL(url);
    }
  };

  return (
    <div className="space-y-4">
      {/* Header */}
      <div className="flex items-center gap-2">
        <Sparkles className="w-5 h-5 text-purple-500" />
        <h2 className="text-lg font-semibold">AI Design Assistant</h2>
      </div>

      {/* Input Card */}
      <Card>
        <CardHeader className="pb-3">
          <CardTitle className="text-base">Describe Your Design</CardTitle>
          <CardDescription>
            Tell me what protein you want to design using natural language
          </CardDescription>
        </CardHeader>
        <CardContent className="space-y-4">
          {/* Query Input */}
          <div className="relative">
            <Textarea
              placeholder="e.g., Design a protein to bind citrate with terbium for luminescent sensing"
              value={query}
              onChange={(e) => setQuery(e.target.value)}
              className="min-h-[100px] pr-12 resize-none"
              disabled={isRunning}
            />
            <Button
              size="icon"
              className="absolute bottom-2 right-2"
              onClick={handleSubmit}
              disabled={!query.trim() || isRunning}
            >
              {isRunning ? (
                <RefreshCw className="w-4 h-4 animate-spin" />
              ) : (
                <Send className="w-4 h-4" />
              )}
            </Button>
          </div>

          {/* Example Queries */}
          <div className="flex flex-wrap gap-2">
            {EXAMPLE_QUERIES.map((example, i) => (
              <button
                key={i}
                onClick={() => handleExampleClick(example)}
                className={cn(
                  "text-xs px-2 py-1 rounded-full border transition-colors",
                  "hover:bg-muted hover:border-primary/30",
                  query === example && "bg-primary/10 border-primary/50"
                )}
                disabled={isRunning}
              >
                {example.length > 40 ? example.slice(0, 40) + '...' : example}
              </button>
            ))}
          </div>

          {/* Advanced Options */}
          <div>
            <button
              onClick={() => setShowAdvanced(!showAdvanced)}
              className="flex items-center gap-1 text-sm text-muted-foreground hover:text-foreground"
            >
              {showAdvanced ? <ChevronUp className="w-4 h-4" /> : <ChevronDown className="w-4 h-4" />}
              Advanced Options
            </button>

            {showAdvanced && (
              <div className="mt-3 space-y-4 p-3 bg-muted/50 rounded-lg">
                <div className="space-y-2">
                  <Label className="text-sm">
                    Number of Backbone Designs: {numDesigns}
                  </Label>
                  <Slider
                    value={[numDesigns]}
                    onValueChange={([v]) => setNumDesigns(v)}
                    min={1}
                    max={10}
                    step={1}
                    disabled={isRunning}
                  />
                </div>
                <div className="space-y-2">
                  <Label className="text-sm">
                    Sequences per Backbone: {numSequences}
                  </Label>
                  <Slider
                    value={[numSequences]}
                    onValueChange={([v]) => setNumSequences(v)}
                    min={4}
                    max={16}
                    step={4}
                    disabled={isRunning}
                  />
                </div>
              </div>
            )}
          </div>
        </CardContent>
      </Card>

      {/* Pipeline Progress */}
      {stage !== 'idle' && (
        <AIDesignPipelineWorkflow
          currentStage={stage}
          stageInfo={stageInfo}
          error={error}
          query={query}
        />
      )}

      {/* Results */}
      {result?.success && (
        <Card>
          <CardHeader className="pb-3">
            <CardTitle className="text-base flex items-center justify-between">
              <span>Design Results</span>
              <div className="flex gap-2">
                <Button
                  variant="outline"
                  size="sm"
                  onClick={handleCopySequence}
                  disabled={!result.best_sequence}
                >
                  {copied ? (
                    <Check className="w-4 h-4 mr-1" />
                  ) : (
                    <Copy className="w-4 h-4 mr-1" />
                  )}
                  Copy Sequence
                </Button>
                <Button
                  variant="outline"
                  size="sm"
                  onClick={handleDownloadPDB}
                  disabled={!result.best_sequence_pdb}
                >
                  <Download className="w-4 h-4 mr-1" />
                  Download PDB
                </Button>
              </div>
            </CardTitle>
          </CardHeader>
          <CardContent className="space-y-4">
            {/* Summary Stats */}
            <div className="grid grid-cols-2 md:grid-cols-4 gap-3">
              <div className="p-3 bg-muted/50 rounded-lg text-center">
                <div className="text-2xl font-bold text-foreground">
                  {result.num_backbones}
                </div>
                <div className="text-xs text-muted-foreground">Backbones</div>
              </div>
              <div className="p-3 bg-muted/50 rounded-lg text-center">
                <div className="text-2xl font-bold text-foreground">
                  {result.num_sequences}
                </div>
                <div className="text-xs text-muted-foreground">Sequences</div>
              </div>
              <div className="p-3 bg-muted/50 rounded-lg text-center">
                <div className={cn(
                  "text-2xl font-bold",
                  result.pass_rate >= 0.2 ? "text-emerald-600" : "text-amber-600"
                )}>
                  {(result.pass_rate * 100).toFixed(0)}%
                </div>
                <div className="text-xs text-muted-foreground">Pass Rate</div>
              </div>
              <div className="p-3 bg-muted/50 rounded-lg text-center">
                <div className="text-2xl font-bold text-foreground">
                  {result.total_time?.toFixed(1)}s
                </div>
                <div className="text-xs text-muted-foreground">Total Time</div>
              </div>
            </div>

            {/* Best Sequence */}
            {result.best_sequence && (
              <div className="space-y-2">
                <Label className="text-sm font-medium">Best Sequence</Label>
                <div className="p-3 bg-muted rounded-lg font-mono text-xs break-all">
                  {result.best_sequence.sequence}
                </div>
                <div className="flex gap-4 text-sm text-muted-foreground">
                  {result.best_rmsd && (
                    <span>RMSD: <span className="text-foreground font-medium">{result.best_rmsd.toFixed(2)}Å</span></span>
                  )}
                  {result.best_plddt && (
                    <span>pLDDT: <span className="text-foreground font-medium">{result.best_plddt.toFixed(2)}</span></span>
                  )}
                </div>
              </div>
            )}

            {/* Recommendations */}
            {result.recommendations && result.recommendations.length > 0 && (
              <div className="space-y-2">
                <Label className="text-sm font-medium">Recommendations</Label>
                <ul className="space-y-1">
                  {result.recommendations.map((rec, i) => (
                    <li key={i} className="text-sm text-muted-foreground flex items-start gap-2">
                      <span className="text-primary">•</span>
                      {rec}
                    </li>
                  ))}
                </ul>
              </div>
            )}

            {/* New Design Button */}
            <Button
              variant="outline"
              className="w-full"
              onClick={() => {
                reset();
                setQuery('');
              }}
            >
              <RefreshCw className="w-4 h-4 mr-2" />
              Start New Design
            </Button>
          </CardContent>
        </Card>
      )}
    </div>
  );
}

export default AIDesignPanel;
