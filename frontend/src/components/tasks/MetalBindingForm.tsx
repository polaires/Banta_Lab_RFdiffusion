'use client';

import { useState, useEffect } from 'react';
import {
  Sparkles,
  Info,
  Diamond,
  Hexagon,
  Layers,
  Settings,
  Ruler,
  Upload,
  CheckCircle,
  ChevronDown,
  ChevronRight,
  Atom,
  FlaskConical,
  Play,
  FileDown,
  XCircle,
} from 'lucide-react';

import { cn } from '@/lib/utils';
import { useStore, type PipelineMode, type SweepConfig } from '@/lib/store';
import { api } from '@/lib/api';

// Shared components
import { SweepConfigForm, generateDefaultSweepConfigs } from './SweepConfigForm';
import { FilterThresholdEditor } from './FilterThresholdEditor';
import { PipelineProgressPanel } from './PipelineProgressPanel';
import { PipelineResultsView } from './PipelineResultsView';

import { Button } from '@/components/ui/button';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Textarea } from '@/components/ui/textarea';
import { Badge } from '@/components/ui/badge';
import { Collapsible, CollapsibleContent, CollapsibleTrigger } from '@/components/ui/collapsible';

// Interview step configuration
interface InterviewOption {
  id: string;
  label: string;
  description: string;
}

interface InterviewStep {
  id: string;
  question: string;
  icon: React.ElementType;
  options: InterviewOption[];
}

const METAL_BINDING_INTERVIEW_STEPS: InterviewStep[] = [
  {
    id: 'target_metal',
    question: 'What metal do you want to bind?',
    icon: Diamond,
    options: [
      { id: 'TB', label: 'Terbium (Tb3+)', description: 'Lanthanide, CN=8-9, luminescent' },
      { id: 'EU', label: 'Europium (Eu3+)', description: 'Lanthanide, red fluorescence' },
      { id: 'CA', label: 'Calcium (Ca2+)', description: 'CN=6-8, signaling' },
      { id: 'ZN', label: 'Zinc (Zn2+)', description: 'CN=4-6, catalytic' },
    ],
  },
  {
    id: 'ligand_type',
    question: 'What ligand is coordinating the metal?',
    icon: Hexagon,
    options: [
      { id: 'CIT', label: 'Citrate', description: 'Tridentate, for lanthanides' },
      { id: 'PQQ', label: 'PQQ', description: 'Cofactor for dehydrogenases' },
      { id: 'none', label: 'None (metal only)', description: 'Design around bare metal' },
      { id: 'custom', label: 'Custom', description: 'Provide your own PDB' },
    ],
  },
  {
    id: 'scaffold_type',
    question: 'What scaffold type?',
    icon: Layers,
    options: [
      { id: 'monomer', label: 'Monomer (Recommended)', description: 'Single chain scaffold, 100-150 aa' },
      { id: 'small_monomer', label: 'Small Monomer', description: 'Compact single chain, 80-100 aa' },
    ],
  },
  {
    id: 'optimization_mode',
    question: 'How should we optimize?',
    icon: Settings,
    options: [
      { id: 'sweep', label: 'Parameter Sweep (Recommended)', description: 'Test 9 configs, find best' },
      { id: 'quick', label: 'Quick (10 designs)', description: 'Single config, fast exploration' },
      { id: 'production', label: 'Production (100+)', description: 'After finding best config' },
    ],
  },
];

// Quality tier display
const TIER_COLORS = {
  S: 'bg-yellow-500 text-black',
  A: 'bg-green-500 text-white',
  B: 'bg-blue-500 text-white',
  C: 'bg-orange-500 text-white',
  F: 'bg-red-500 text-white',
};

// Phase type
type FormPhase = 'input' | 'interview' | 'configure' | 'running' | 'results';

interface MetalBindingFormProps {
  health: boolean;
}

export function MetalBindingForm({ health }: MetalBindingFormProps) {
  // Phase state
  const [phase, setPhase] = useState<FormPhase>('input');
  const [currentStep, setCurrentStep] = useState(0);

  // Interview answers
  const [answers, setAnswers] = useState<Record<string, string>>({});

  // PDB input
  const [pdbInput, setPdbInput] = useState('');
  const [pdbCode, setPdbCode] = useState('');
  const [isLoadingPdb, setIsLoadingPdb] = useState(false);
  const [pdbError, setPdbError] = useState<string | null>(null);

  // Analysis result
  const [analysis, setAnalysis] = useState<{
    metals: Array<{ element: string; chain: string; res_num: string }>;
    ligands: Array<{ name: string; chain: string; res_num: string }>;
    suggestions: string[];
    recommended_hsab_bias?: string;
  } | null>(null);

  // Pipeline state from store
  const {
    pipelineState,
    setPipelineMode,
    startPipeline,
    updatePipelineProgress,
    setPipelineResults,
    setPipelineFilters,
    setSweepConfigs,
    resetPipeline,
    cancelPipeline,
  } = useStore();

  // Local state
  const [designsPerConfig, setDesignsPerConfig] = useState(10);
  const [productionDesigns, setProductionDesigns] = useState(100);

  // Initialize sweep configs when entering sweep mode
  useEffect(() => {
    if (pipelineState.mode === 'sweep' && pipelineState.sweepConfigs.length === 0) {
      // Use larger contig ranges for monomer design
      const monomerConfigs = generateMonomerSweepConfigs(designsPerConfig);
      setSweepConfigs(monomerConfigs);
    }
  }, [pipelineState.mode]);

  // Generate monomer-appropriate sweep configs
  function generateMonomerSweepConfigs(numDesigns: number): SweepConfig[] {
    const sizes: Array<'small' | 'medium' | 'large'> = ['small', 'medium', 'large'];
    const cfgValues = [1.5, 2.0, 2.5];

    // Monomer-specific ranges (larger than dimer)
    const monomerRanges: Record<'small' | 'medium' | 'large', string> = {
      small: '100-120',
      medium: '110-130',
      large: '130-150',
    };

    const configs: SweepConfig[] = [];
    for (const size of sizes) {
      for (const cfg of cfgValues) {
        const cfgName = cfg === 1.5 ? 'low_cfg' : cfg === 2.0 ? 'mid_cfg' : 'high_cfg';
        configs.push({
          name: `${size}_${cfgName}`,
          contigSize: size,
          contigRange: monomerRanges[size],
          cfgScale: cfg,
          numDesigns: numDesigns,
        });
      }
    }
    return configs;
  }

  // Handle PDB code lookup
  const handlePdbLookup = async () => {
    if (!pdbCode.trim()) return;

    setIsLoadingPdb(true);
    setPdbError(null);

    try {
      const response = await api.fetchPdb(pdbCode);
      if (response.success) {
        setPdbInput(response.content);

        // Run analysis
        const analysisResult = await runAnalysis(response.content);
        if (analysisResult) {
          setPhase('interview');
        }
      } else {
        setPdbError(response.error || 'Failed to fetch PDB');
      }
    } catch (error) {
      setPdbError(error instanceof Error ? error.message : 'Failed to fetch PDB');
    } finally {
      setIsLoadingPdb(false);
    }
  };

  // Run structure analysis
  const runAnalysis = async (pdbContent: string): Promise<boolean> => {
    try {
      // Call backend analysis
      const response = await fetch('/api/traditional/runsync?url=' + encodeURIComponent(useStore.getState().backendUrl), {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          input: {
            task: 'metal_binding_design',
            mode: 'analyze',
            pdb_content: pdbContent,
          },
        }),
      });

      const result = await response.json();
      if (result.output?.status === 'completed' && result.output?.result) {
        setAnalysis(result.output.result);

        // Auto-select metal if found
        if (result.output.result.metals?.length > 0) {
          setAnswers((prev) => ({
            ...prev,
            target_metal: result.output.result.metals[0].element,
          }));
        }

        // Auto-select ligand if found
        if (result.output.result.ligands?.length > 0) {
          setAnswers((prev) => ({
            ...prev,
            ligand_type: result.output.result.ligands[0].name,
          }));
        }

        return true;
      }
    } catch (error) {
      console.error('Analysis failed:', error);
    }
    return false;
  };

  // Handle interview answer
  const handleAnswer = (stepId: string, optionId: string) => {
    setAnswers((prev) => ({ ...prev, [stepId]: optionId }));

    // Move to next step or configure phase
    if (currentStep < METAL_BINDING_INTERVIEW_STEPS.length - 1) {
      setCurrentStep((prev) => prev + 1);
    } else {
      setPhase('configure');
      // Set pipeline mode based on answer
      const mode = answers.optimization_mode || optionId;
      if (mode === 'sweep') {
        setPipelineMode('sweep');
      } else if (mode === 'production') {
        setPipelineMode('production');
      } else {
        setPipelineMode('single');
      }
    }
  };

  // Start the pipeline
  const handleStartPipeline = async () => {
    const metal = answers.target_metal || 'TB';
    const ligand = answers.ligand_type !== 'none' && answers.ligand_type !== 'custom'
      ? answers.ligand_type
      : undefined;

    setPhase('running');

    try {
      if (pipelineState.mode === 'sweep') {
        const backendConfigs = pipelineState.sweepConfigs.map((c) => ({
          name: c.name,
          contig_size: c.contigSize,
          contig_range: c.contigRange,
          cfg_scale: c.cfgScale,
          num_designs: c.numDesigns,
        }));

        const result = await fetch('/api/traditional/runsync?url=' + encodeURIComponent(useStore.getState().backendUrl), {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            input: {
              task: 'metal_binding_design',
              mode: 'sweep',
              motif_pdb: pdbInput,
              metal,
              ligand,
              sweep_configs: backendConfigs,
              designs_per_config: designsPerConfig,
              filters: pipelineState.filters,
            },
          }),
        });

        const data = await result.json();
        if (!result.ok) {
          const errMsg = data.error || `Backend returned ${result.status}`;
          console.error('Pipeline API error:', errMsg);
          throw new Error(errMsg);
        }
        if (data.output?.status === 'failed') {
          const errMsg = data.output.error || 'Backend returned failure';
          console.error('Pipeline execution error:', errMsg);
          throw new Error(errMsg);
        }
        if (data.output?.status === 'completed') {
          const pipelineResult = data.output.result;

          startPipeline(
            pipelineResult.session_id,
            'sweep',
            pipelineResult.total_configs,
            pipelineResult.designs_per_config || designsPerConfig
          );

          // Update with final results
          updatePipelineProgress({
            currentConfig: pipelineResult.total_configs,
            totalConfigs: pipelineResult.total_configs,
            currentDesign: pipelineResult.designs_per_config,
            designsPerConfig: pipelineResult.designs_per_config,
            totalGenerated: pipelineResult.total_generated,
            totalPassing: pipelineResult.total_passing,
            totalReview: pipelineResult.total_review,
            totalFailed: pipelineResult.total_generated - pipelineResult.total_passing - pipelineResult.total_review,
            passRate: pipelineResult.pass_rate,
            bestDesign: pipelineResult.best_design,
          });

          if (pipelineResult.results) {
            setPipelineResults(
              pipelineResult.results.map((r: any) => ({
                name: r.name,
                config: r.config_name,
                sequence: r.sequence || '',
                plddt: r.plddt,
                ptm: r.ptm,
                pae: r.pae,
                status: r.status as 'pass' | 'review' | 'fail',
              }))
            );
          }

          setPhase('results');
        }
      }
    } catch (error) {
      console.error('Pipeline failed:', error);
      setPhase('configure');
    }
  };

  // Export FASTA
  const handleExportFasta = async (selectedIds: string[], includeReview: boolean) => {
    // Create FASTA from results
    const passingResults = pipelineState.results.filter(
      (r) => r.status === 'pass' || (includeReview && r.status === 'review')
    );

    const fastaContent = passingResults
      .map((r) => `>${r.name} pLDDT=${r.plddt.toFixed(2)} pTM=${r.ptm.toFixed(2)}\n${r.sequence}`)
      .join('\n\n');

    // Download
    const blob = new Blob([fastaContent], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `metal_binding_designs_${Date.now()}.fasta`;
    a.click();
    URL.revokeObjectURL(url);
  };

  // Reset workflow
  const handleReset = () => {
    setPhase('input');
    setCurrentStep(0);
    setAnswers({});
    setPdbInput('');
    setPdbCode('');
    setAnalysis(null);
    resetPipeline();
  };

  // Render based on phase
  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex items-center gap-3 pb-4 border-b border-border">
        <div className="w-10 h-10 rounded-lg bg-primary flex items-center justify-center">
          <Diamond className="w-5 h-5 text-primary-foreground" />
        </div>
        <div>
          <h2 className="font-semibold text-foreground">Metal Binding Site Design</h2>
          <p className="text-sm text-muted-foreground">
            Design monomer scaffolds around metal-ligand complexes (Round 7b style)
          </p>
        </div>
      </div>

      {/* Phase: Input */}
      {phase === 'input' && (
        <Card>
          <CardHeader>
            <CardTitle className="text-sm flex items-center gap-2">
              <Upload className="w-4 h-4" />
              Structure Input
            </CardTitle>
            <CardDescription>
              Provide a PDB with your metal-ligand complex or fetch from RCSB
            </CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            {/* PDB code lookup */}
            <div className="space-y-2">
              <Label>PDB Code</Label>
              <div className="flex gap-2">
                <Input
                  value={pdbCode}
                  onChange={(e) => setPdbCode(e.target.value.toUpperCase())}
                  placeholder="e.g., 3C9H"
                  className="font-mono"
                />
                <Button
                  onClick={handlePdbLookup}
                  disabled={!pdbCode.trim() || isLoadingPdb}
                >
                  {isLoadingPdb ? 'Loading...' : 'Fetch'}
                </Button>
              </div>
              {pdbError && (
                <p className="text-sm text-destructive">{pdbError}</p>
              )}
            </div>

            <div className="relative">
              <div className="absolute inset-0 flex items-center">
                <span className="w-full border-t" />
              </div>
              <div className="relative flex justify-center text-xs uppercase">
                <span className="bg-card px-2 text-muted-foreground">Or paste PDB</span>
              </div>
            </div>

            {/* PDB content */}
            <div className="space-y-2">
              <Label>PDB Content</Label>
              <Textarea
                value={pdbInput}
                onChange={(e) => setPdbInput(e.target.value)}
                placeholder="HETATM    1  C1  CIT L   1    ..."
                className="font-mono h-40"
              />
            </div>

            <Button
              onClick={async () => {
                if (pdbInput.trim()) {
                  const success = await runAnalysis(pdbInput);
                  if (success) setPhase('interview');
                }
              }}
              disabled={!pdbInput.trim()}
              className="w-full"
            >
              <ChevronRight className="w-4 h-4 mr-2" />
              Analyze Structure
            </Button>
          </CardContent>
        </Card>
      )}

      {/* Phase: Interview */}
      {phase === 'interview' && (
        <div className="space-y-4">
          {/* Analysis summary */}
          {analysis && (
            <Card className="bg-muted/50">
              <CardContent className="pt-4">
                <div className="flex items-start gap-3">
                  <Info className="w-5 h-5 text-primary mt-0.5 flex-shrink-0" />
                  <div className="text-sm">
                    <p className="font-medium mb-1">Structure Analysis</p>
                    <p className="text-muted-foreground">
                      Found {analysis.metals.length} metal(s)
                      {analysis.metals.length > 0 && ` (${analysis.metals.map((m) => m.element).join(', ')})`}
                      {analysis.ligands.length > 0 && `, ${analysis.ligands.length} ligand(s) (${analysis.ligands.map((l) => l.name).join(', ')})`}
                    </p>
                    {analysis.suggestions.map((s, i) => (
                      <p key={i} className="text-muted-foreground text-xs mt-1">{s}</p>
                    ))}
                  </div>
                </div>
              </CardContent>
            </Card>
          )}

          {/* Interview question */}
          <Card>
            <CardHeader>
              <CardTitle className="text-sm flex items-center gap-2">
                {(() => {
                  const Icon = METAL_BINDING_INTERVIEW_STEPS[currentStep].icon;
                  return <Icon className="w-4 h-4" />;
                })()}
                Step {currentStep + 1} of {METAL_BINDING_INTERVIEW_STEPS.length}
              </CardTitle>
              <CardDescription className="text-base font-medium text-foreground">
                {METAL_BINDING_INTERVIEW_STEPS[currentStep].question}
              </CardDescription>
            </CardHeader>
            <CardContent>
              <div className="grid grid-cols-1 sm:grid-cols-2 gap-3">
                {METAL_BINDING_INTERVIEW_STEPS[currentStep].options.map((option) => (
                  <button
                    key={option.id}
                    onClick={() => handleAnswer(METAL_BINDING_INTERVIEW_STEPS[currentStep].id, option.id)}
                    className={cn(
                      'p-4 rounded-lg border text-left transition-all',
                      answers[METAL_BINDING_INTERVIEW_STEPS[currentStep].id] === option.id
                        ? 'border-primary bg-primary/5 ring-1 ring-primary/20'
                        : 'border-border hover:border-muted-foreground/30 hover:bg-muted/50'
                    )}
                  >
                    <div className="font-medium text-foreground text-sm">{option.label}</div>
                    <p className="text-xs text-muted-foreground mt-1">{option.description}</p>
                  </button>
                ))}
              </div>
            </CardContent>
          </Card>

          {/* Progress dots */}
          <div className="flex justify-center gap-2">
            {METAL_BINDING_INTERVIEW_STEPS.map((_, i) => (
              <div
                key={i}
                className={cn(
                  'w-2 h-2 rounded-full transition-colors',
                  i < currentStep ? 'bg-primary' : i === currentStep ? 'bg-primary/50' : 'bg-muted'
                )}
              />
            ))}
          </div>
        </div>
      )}

      {/* Phase: Configure */}
      {phase === 'configure' && (
        <div className="space-y-4">
          {/* Configuration summary */}
          <Card>
            <CardHeader>
              <CardTitle className="text-sm flex items-center gap-2">
                <CheckCircle className="w-4 h-4 text-primary" />
                Configuration Summary
              </CardTitle>
            </CardHeader>
            <CardContent>
              <div className="grid grid-cols-2 gap-4 text-sm">
                <div>
                  <span className="text-muted-foreground">Metal:</span>
                  <span className="ml-2 font-medium">{answers.target_metal || 'Not selected'}</span>
                </div>
                <div>
                  <span className="text-muted-foreground">Ligand:</span>
                  <span className="ml-2 font-medium">
                    {answers.ligand_type === 'none' ? 'None' : answers.ligand_type || 'Not selected'}
                  </span>
                </div>
                <div>
                  <span className="text-muted-foreground">Scaffold:</span>
                  <span className="ml-2 font-medium">{answers.scaffold_type || 'Monomer'}</span>
                </div>
                <div>
                  <span className="text-muted-foreground">Mode:</span>
                  <span className="ml-2 font-medium">
                    {answers.optimization_mode === 'sweep' ? 'Parameter Sweep' :
                     answers.optimization_mode === 'production' ? 'Production' : 'Quick'}
                  </span>
                </div>
              </div>
            </CardContent>
          </Card>

          {/* Sweep configuration */}
          {pipelineState.mode === 'sweep' && (
            <SweepConfigForm
              configs={pipelineState.sweepConfigs}
              onConfigsChange={setSweepConfigs}
              designsPerConfig={designsPerConfig}
              onDesignsPerConfigChange={setDesignsPerConfig}
              disabled={pipelineState.isRunning}
            />
          )}

          {/* Production configuration */}
          {pipelineState.mode === 'production' && (
            <Card>
              <CardHeader>
                <CardTitle className="text-sm flex items-center gap-2">
                  <Layers className="w-4 h-4" />
                  Production Configuration
                </CardTitle>
              </CardHeader>
              <CardContent className="space-y-4">
                <div className="space-y-2">
                  <Label>Number of Designs</Label>
                  <Input
                    type="number"
                    min={10}
                    max={1000}
                    value={productionDesigns}
                    onChange={(e) => setProductionDesigns(parseInt(e.target.value) || 100)}
                  />
                </div>
              </CardContent>
            </Card>
          )}

          {/* Filter thresholds */}
          <FilterThresholdEditor
            filters={pipelineState.filters}
            onFiltersChange={setPipelineFilters}
            disabled={pipelineState.isRunning}
          />

          {/* Action buttons */}
          <div className="flex gap-2">
            <Button variant="outline" onClick={handleReset}>
              Start Over
            </Button>
            <Button onClick={handleStartPipeline} className="flex-1" disabled={!health}>
              <Play className="w-4 h-4 mr-2" />
              {pipelineState.mode === 'sweep'
                ? `Start Sweep (${pipelineState.sweepConfigs.length} configs × ${designsPerConfig} designs)`
                : `Start Production (${productionDesigns} designs)`}
            </Button>
          </div>
        </div>
      )}

      {/* Phase: Running */}
      {phase === 'running' && (
        <div className="space-y-4">
          <PipelineProgressPanel
            progress={pipelineState.progress}
            mode={pipelineState.mode}
            isRunning={pipelineState.isRunning}
            onCancel={() => {
              cancelPipeline();
              setPhase('configure');
            }}
          />
        </div>
      )}

      {/* Phase: Results */}
      {phase === 'results' && (
        <div className="space-y-4">
          {/* Summary */}
          <Card className="bg-muted/50">
            <CardContent className="pt-4">
              <div className="flex items-center gap-4">
                <div className="w-12 h-12 rounded-full bg-primary/20 flex items-center justify-center">
                  <CheckCircle className="w-6 h-6 text-primary" />
                </div>
                <div>
                  <h3 className="font-semibold">Sweep Complete</h3>
                  <p className="text-sm text-muted-foreground">
                    {pipelineState.progress.totalPassing} passing designs from {pipelineState.progress.totalGenerated} total
                    ({(pipelineState.progress.passRate * 100).toFixed(1)}% pass rate)
                  </p>
                </div>
              </div>
            </CardContent>
          </Card>

          {/* Results table */}
          <PipelineResultsView
            results={pipelineState.results}
            onExportFasta={handleExportFasta}
          />

          {/* Action buttons */}
          <div className="flex gap-2">
            <Button variant="outline" onClick={handleReset}>
              New Design
            </Button>
            <Button
              onClick={() => setPhase('configure')}
              variant="outline"
            >
              Run Production with Best Config
            </Button>
          </div>
        </div>
      )}

      {/* Backend status */}
      {!health && (
        <p className="text-center text-sm text-muted-foreground">
          Backend service unavailable
        </p>
      )}
    </div>
  );
}

export default MetalBindingForm;
