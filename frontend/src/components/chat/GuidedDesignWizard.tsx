'use client';

import { useState, useMemo } from 'react';
import {
  Gem,
  FlaskConical,
  Target,
  Gauge,
  Layers,
  Network,
  Shield,
  CheckCircle,
  ArrowLeft,
  ArrowRight,
  Rocket,
  type LucideIcon,
} from 'lucide-react';

// ---- Types ----

interface WizardOption {
  id: string;
  label: string;
  description: string;
}

interface WizardStep {
  id: string;
  question: string;
  description: string;
  options: WizardOption[];
  icon: LucideIcon;
}

export interface WizardResult {
  design_type: string;
  target_metal?: string;
  ligand_name?: string;
  design_goal?: string;
  filter_tier?: string;
  num_designs?: number;
  use_sweep?: boolean;
  symmetry?: string;
  target_pdb_id?: string;
  protocol?: string;
  user_prompt: string;
}

interface GuidedDesignWizardProps {
  knownParams?: Record<string, string>;
  onComplete: (result: WizardResult) => void;
  onCancel: () => void;
}

// ---- Step definitions per design type ----

const DESIGN_TYPE_STEP: WizardStep = {
  id: 'design_type',
  question: 'What kind of protein would you like to design?',
  description: 'This determines the pipeline configuration and available parameters.',
  icon: Layers,
  options: [
    { id: 'metal', label: 'Metal-Ligand Binding', description: 'Design a protein that coordinates a metal ion, optionally with an organic ligand' },
    { id: 'scaffold', label: 'Scaffold Design', description: 'Search for and build on existing metal-binding scaffold structures from PDB' },
    { id: 'ligand', label: 'Interface Dimer', description: 'Design two protein chains connected by a small molecule at the interface' },
    { id: 'binder', label: 'Protein Binder', description: 'Design a new protein that binds to an existing target protein' },
    { id: 'general', label: 'Something else', description: 'Describe your design goal in your own words' },
  ],
};

const METAL_STEP: WizardStep = {
  id: 'target_metal',
  question: 'Which metal ion should the protein coordinate?',
  description: 'The metal determines coordination geometry and residue preferences.',
  icon: Gem,
  options: [
    { id: 'TB', label: 'Terbium (Tb)', description: 'Lanthanide, 8-9 coordination, luminescent properties' },
    { id: 'CA', label: 'Calcium (Ca)', description: 'Common signaling ion, 6-8 coordination' },
    { id: 'ZN', label: 'Zinc (Zn)', description: 'Catalytic metal, tetrahedral or octahedral' },
    { id: 'EU', label: 'Europium (Eu)', description: 'Lanthanide with strong red fluorescence' },
    { id: 'GD', label: 'Gadolinium (Gd)', description: 'Lanthanide, MRI contrast applications' },
  ],
};

const LIGAND_STEP: WizardStep = {
  id: 'ligand_name',
  question: 'What ligand coordinates the metal?',
  description: 'An organic ligand can bridge between the protein and metal, or be the main binding target.',
  icon: FlaskConical,
  options: [
    { id: 'citrate', label: 'Citrate', description: 'Tricarboxylate chelator, common in biology' },
    { id: 'PQQ', label: 'PQQ', description: 'Pyrroloquinoline quinone, redox cofactor' },
    { id: 'ATP', label: 'ATP', description: 'Adenosine triphosphate, universal energy carrier' },
    { id: 'none', label: 'No ligand', description: 'Direct protein-metal coordination only' },
  ],
};

const DESIGN_GOAL_STEP: WizardStep = {
  id: 'design_goal',
  question: "What's the biological purpose of this design?",
  description: 'This affects how the binding site is shaped — buried for tight binding, accessible for catalysis.',
  icon: Target,
  options: [
    { id: 'binding', label: 'Tight binding', description: 'Maximize coordination number, bury the metal site for stability' },
    { id: 'catalysis', label: 'Enzymatic catalysis', description: 'Accessible active site with substrate channel' },
    { id: 'sensing', label: 'Sensing / detection', description: 'Exposed metal site for signal readout or fluorescence' },
    { id: 'structural', label: 'Structural role', description: 'Metal as part of the protein scaffold, maximize fold stability' },
  ],
};

const SCALE_STEP: WizardStep = {
  id: 'scale',
  question: 'How thorough should the design run be?',
  description: 'More designs give better coverage but take longer to compute.',
  icon: Gauge,
  options: [
    { id: 'quick', label: 'Quick exploration', description: '5 designs, fast results for initial testing' },
    { id: 'standard', label: 'Standard run', description: '10 designs, good coverage for most goals' },
    { id: 'sweep', label: 'Parameter sweep', description: '50+ designs across multiple configurations' },
    { id: 'production', label: 'Production', description: '200+ designs, maximum coverage' },
  ],
};

const QUALITY_STEP: WizardStep = {
  id: 'filter_tier',
  question: 'How strict should quality filtering be?',
  description: 'Stricter filters reject more designs but ensure higher confidence in results.',
  icon: Shield,
  options: [
    { id: 'relaxed', label: 'Exploratory', description: 'Lenient thresholds — keep more designs for analysis' },
    { id: 'standard', label: 'Standard', description: 'Balanced quality bar for most design goals' },
    { id: 'strict', label: 'Strict', description: 'High confidence only — suitable for production candidates' },
  ],
};

// Dimer-specific
const DIMER_LIGAND_STEP: WizardStep = {
  id: 'ligand_name',
  question: 'What molecule bridges the two chains?',
  description: 'The interface ligand mediates the dimer interaction.',
  icon: FlaskConical,
  options: [
    { id: 'azobenzene', label: 'Azobenzene', description: 'Light-switchable linker for reversible dimerization' },
    { id: 'citrate', label: 'Citrate', description: 'Tricarboxylate for stable bridging' },
    { id: 'ATP', label: 'ATP', description: 'Nucleotide-mediated interface' },
  ],
};

const SYMMETRY_STEP: WizardStep = {
  id: 'symmetry',
  question: 'Should the dimer chains be symmetric?',
  description: 'Symmetric dimers have identical chains; asymmetric allows different folds.',
  icon: Network,
  options: [
    { id: 'C2', label: 'Symmetric (C2)', description: 'Both chains identical — simpler to express and characterize' },
    { id: 'asymmetric', label: 'Asymmetric', description: 'Different chains — more design freedom' },
  ],
};

// Binder-specific
const BINDER_TARGET_STEP: WizardStep = {
  id: 'target_pdb_id',
  question: 'What is the target protein?',
  description: 'Enter a PDB code for the protein you want to bind.',
  icon: Target,
  options: [
    { id: 'custom', label: 'Enter PDB code', description: 'I\'ll type a PDB ID in the text box after this wizard' },
    { id: '1BRF', label: '1BRF (Rubredoxin)', description: 'Small iron-binding protein, good test target' },
    { id: '6LZG', label: '6LZG (SARS-CoV-2 RBD)', description: 'Receptor binding domain, therapeutic target' },
  ],
};

const BINDER_PROTOCOL_STEP: WizardStep = {
  id: 'protocol',
  question: 'Which design protocol?',
  description: 'Different protocols trade off speed vs design quality.',
  icon: Rocket,
  options: [
    { id: 'standard', label: 'Standard RFdiffusion', description: 'Fast backbone generation with hotspot targeting' },
    { id: 'partial_diffusion', label: 'Partial diffusion', description: 'Start from existing scaffold, refine for binding' },
  ],
};

// ---- Helper: build steps for a design type, skipping known params ----

function buildSteps(designType: string, known: Record<string, string>): WizardStep[] {
  const steps: WizardStep[] = [];

  if (designType === 'metal' || designType === 'scaffold') {
    if (!known.target_metal) steps.push(METAL_STEP);
    if (!known.ligand_name) steps.push(LIGAND_STEP);
    if (!known.design_goal) steps.push(DESIGN_GOAL_STEP);
    if (!known.scale) steps.push(SCALE_STEP);
    if (!known.filter_tier) steps.push(QUALITY_STEP);
  } else if (designType === 'ligand') {
    if (!known.ligand_name) steps.push(DIMER_LIGAND_STEP);
    if (!known.symmetry) steps.push(SYMMETRY_STEP);
    if (!known.scale) steps.push(SCALE_STEP);
    if (!known.filter_tier) steps.push(QUALITY_STEP);
  } else if (designType === 'binder') {
    if (!known.target_pdb_id) steps.push(BINDER_TARGET_STEP);
    if (!known.protocol) steps.push(BINDER_PROTOCOL_STEP);
    if (!known.scale) steps.push(SCALE_STEP);
    if (!known.filter_tier) steps.push(QUALITY_STEP);
  } else {
    // general — just scale + quality
    if (!known.scale) steps.push(SCALE_STEP);
    if (!known.filter_tier) steps.push(QUALITY_STEP);
  }

  return steps;
}

// ---- Scale → num_designs + use_sweep mapping ----

function scaleToParams(scale: string): { num_designs: number; use_sweep: boolean } {
  switch (scale) {
    case 'quick': return { num_designs: 5, use_sweep: false };
    case 'standard': return { num_designs: 10, use_sweep: false };
    case 'sweep': return { num_designs: 50, use_sweep: true };
    case 'production': return { num_designs: 200, use_sweep: true };
    default: return { num_designs: 10, use_sweep: false };
  }
}

// ---- Label lookups for review ----

const METAL_LABELS: Record<string, string> = {
  TB: 'Terbium (Tb)', CA: 'Calcium (Ca)', ZN: 'Zinc (Zn)',
  EU: 'Europium (Eu)', GD: 'Gadolinium (Gd)',
};

const DESIGN_TYPE_LABELS: Record<string, string> = {
  metal: 'Metal-Ligand Binding', scaffold: 'Scaffold Design',
  ligand: 'Interface Dimer', binder: 'Protein Binder', general: 'General Design',
};

const GOAL_LABELS: Record<string, string> = {
  binding: 'Tight binding', catalysis: 'Enzymatic catalysis',
  sensing: 'Sensing / detection', structural: 'Structural role',
};

const SCALE_LABELS: Record<string, string> = {
  quick: 'Quick (5)', standard: 'Standard (10)',
  sweep: 'Sweep (50+)', production: 'Production (200+)',
};

const TIER_LABELS: Record<string, string> = {
  relaxed: 'Exploratory', standard: 'Standard', strict: 'Strict',
};

// ---- Component ----

export function GuidedDesignWizard({ knownParams = {}, onComplete, onCancel }: GuidedDesignWizardProps) {
  // Phase: 'design_type' selection → 'questions' → 'review'
  const [phase, setPhase] = useState<'design_type' | 'questions' | 'review'>(
    knownParams.design_type ? 'questions' : 'design_type'
  );
  const [designType, setDesignType] = useState(knownParams.design_type || '');
  const [currentStep, setCurrentStep] = useState(0);
  const [answers, setAnswers] = useState<Record<string, string>>({ ...knownParams });
  const [selectedOption, setSelectedOption] = useState<string | null>(null);
  // For 'general' type, allow free-text input
  const [freeText, setFreeText] = useState('');

  // Build follow-up steps based on design type + known params
  const followUpSteps = useMemo(
    () => (designType ? buildSteps(designType, knownParams) : []),
    [designType, knownParams]
  );

  // ---- Design type selection phase ----

  if (phase === 'design_type') {
    const step = DESIGN_TYPE_STEP;
    return (
      <div className="bg-muted rounded-2xl p-6 border border-border shadow-sm">
        <div className="flex items-center gap-2 text-primary text-sm mb-4">
          <step.icon className="h-5 w-5" />
          <span className="font-medium">Guided Design</span>
        </div>

        <h3 className="text-xl font-semibold text-foreground mb-2">{step.question}</h3>
        <p className="text-muted-foreground text-sm mb-6">{step.description}</p>

        <div className="space-y-3 mb-6">
          {step.options.map((opt) => {
            const isSelected = selectedOption === opt.id;
            return (
              <button
                key={opt.id}
                onClick={() => setSelectedOption(opt.id)}
                className={`w-full text-left p-4 rounded-xl border-2 transition-all ${
                  isSelected
                    ? 'border-primary bg-primary/5 shadow-sm'
                    : 'border-border hover:border-primary/50 hover:bg-muted/50'
                }`}
              >
                <div className="flex items-center justify-between">
                  <div>
                    <div className="font-medium text-foreground">{opt.label}</div>
                    <div className="text-sm text-muted-foreground">{opt.description}</div>
                  </div>
                  {isSelected && <CheckCircle className="h-6 w-6 text-primary flex-shrink-0" />}
                </div>
              </button>
            );
          })}
        </div>

        {/* Show free-text input when "Something else" is selected */}
        {selectedOption === 'general' && (
          <div className="mb-6">
            <textarea
              value={freeText}
              onChange={(e) => setFreeText(e.target.value)}
              placeholder="Describe what you'd like to design..."
              className="w-full px-4 py-3 bg-card rounded-xl border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 focus:outline-none text-foreground text-sm resize-none min-h-[80px]"
              rows={3}
            />
          </div>
        )}

        <div className="flex items-center justify-between">
          <button
            onClick={onCancel}
            className="text-muted-foreground text-sm font-medium hover:text-foreground"
          >
            Cancel
          </button>
          <button
            onClick={() => {
              if (!selectedOption) return;
              setDesignType(selectedOption);
              setAnswers((prev) => ({ ...prev, design_type: selectedOption }));
              if (selectedOption === 'general' && freeText.trim()) {
                // "Something else" with text → skip questions, go to complete
                handleFinalComplete({ ...knownParams, ...answers, design_type: 'general' }, freeText.trim());
                return;
              }
              setSelectedOption(null);
              setCurrentStep(0);
              setPhase('questions');
            }}
            disabled={!selectedOption || (selectedOption === 'general' && !freeText.trim())}
            className={`px-6 py-2.5 rounded-xl font-medium transition-all flex items-center gap-2 ${
              selectedOption && !(selectedOption === 'general' && !freeText.trim())
                ? 'bg-primary text-primary-foreground hover:bg-primary/90 shadow-sm'
                : 'bg-muted text-muted-foreground cursor-not-allowed'
            }`}
          >
            {selectedOption === 'general' ? 'Start Design' : 'Next'}
            <ArrowRight className="h-4 w-4" />
          </button>
        </div>
      </div>
    );
  }

  // ---- Questions phase ----

  if (phase === 'questions' && followUpSteps.length > 0 && currentStep < followUpSteps.length) {
    const step = followUpSteps[currentStep];
    const progress = ((currentStep + 1) / (followUpSteps.length + 1)) * 100; // +1 for review
    const isLast = currentStep === followUpSteps.length - 1;

    return (
      <div className="bg-muted rounded-2xl p-6 border border-border shadow-sm">
        {/* Progress bar */}
        <div className="h-1.5 bg-muted-foreground/20 rounded-full mb-6 overflow-hidden">
          <div
            className="h-full bg-primary rounded-full transition-all duration-500"
            style={{ width: `${progress}%` }}
          />
        </div>

        <div className="flex items-center gap-2 text-primary text-sm mb-4">
          <step.icon className="h-5 w-5" />
          <span className="font-medium">
            Step {currentStep + 1} of {followUpSteps.length}
            <span className="text-muted-foreground ml-2">
              {DESIGN_TYPE_LABELS[designType] || designType}
            </span>
          </span>
        </div>

        <h3 className="text-xl font-semibold text-foreground mb-2">{step.question}</h3>
        <p className="text-muted-foreground text-sm mb-6">{step.description}</p>

        <div className="space-y-3 mb-6">
          {step.options.map((opt) => {
            const isSelected = selectedOption === opt.id;
            return (
              <button
                key={opt.id}
                onClick={() => setSelectedOption(opt.id)}
                className={`w-full text-left p-4 rounded-xl border-2 transition-all ${
                  isSelected
                    ? 'border-primary bg-primary/5 shadow-sm'
                    : 'border-border hover:border-primary/50 hover:bg-muted/50'
                }`}
              >
                <div className="flex items-center justify-between">
                  <div>
                    <div className="font-medium text-foreground">{opt.label}</div>
                    <div className="text-sm text-muted-foreground">{opt.description}</div>
                  </div>
                  {isSelected && <CheckCircle className="h-6 w-6 text-primary flex-shrink-0" />}
                </div>
              </button>
            );
          })}
        </div>

        <div className="flex items-center justify-between">
          <button
            onClick={() => {
              if (currentStep > 0) {
                setCurrentStep((s) => s - 1);
                setSelectedOption(answers[followUpSteps[currentStep - 1].id] || null);
              } else {
                // Go back to design type selection
                setPhase('design_type');
                setSelectedOption(designType);
              }
            }}
            className="text-primary text-sm font-medium hover:text-primary/80 flex items-center gap-1"
          >
            <ArrowLeft className="h-4 w-4" />
            Back
          </button>
          <button
            onClick={() => {
              if (!selectedOption) return;
              const newAnswers = { ...answers, [step.id]: selectedOption };
              setAnswers(newAnswers);
              if (isLast) {
                setPhase('review');
              } else {
                setCurrentStep((s) => s + 1);
              }
              setSelectedOption(null);
            }}
            disabled={!selectedOption}
            className={`px-6 py-2.5 rounded-xl font-medium transition-all flex items-center gap-2 ${
              selectedOption
                ? 'bg-primary text-primary-foreground hover:bg-primary/90 shadow-sm'
                : 'bg-muted text-muted-foreground cursor-not-allowed'
            }`}
          >
            {isLast ? 'Review' : 'Next'}
            {isLast ? <CheckCircle className="h-4 w-4" /> : <ArrowRight className="h-4 w-4" />}
          </button>
        </div>
      </div>
    );
  }

  // ---- Review phase ----

  const allAnswers = { ...knownParams, ...answers };

  // Build a human-readable summary
  const summaryItems: Array<{ label: string; value: string }> = [];
  summaryItems.push({ label: 'Design Type', value: DESIGN_TYPE_LABELS[allAnswers.design_type] || allAnswers.design_type || 'General' });

  if (allAnswers.target_metal) {
    summaryItems.push({ label: 'Metal', value: METAL_LABELS[allAnswers.target_metal] || allAnswers.target_metal });
  }
  if (allAnswers.ligand_name && allAnswers.ligand_name !== 'none') {
    summaryItems.push({ label: 'Ligand', value: allAnswers.ligand_name });
  }
  if (allAnswers.design_goal) {
    summaryItems.push({ label: 'Goal', value: GOAL_LABELS[allAnswers.design_goal] || allAnswers.design_goal });
  }
  if (allAnswers.symmetry) {
    summaryItems.push({ label: 'Symmetry', value: allAnswers.symmetry === 'C2' ? 'Symmetric (C2)' : 'Asymmetric' });
  }
  if (allAnswers.target_pdb_id && allAnswers.target_pdb_id !== 'custom') {
    summaryItems.push({ label: 'Target', value: allAnswers.target_pdb_id });
  }
  if (allAnswers.protocol) {
    summaryItems.push({ label: 'Protocol', value: allAnswers.protocol });
  }
  if (allAnswers.scale) {
    summaryItems.push({ label: 'Scale', value: SCALE_LABELS[allAnswers.scale] || allAnswers.scale });
  }
  if (allAnswers.filter_tier) {
    summaryItems.push({ label: 'Quality', value: TIER_LABELS[allAnswers.filter_tier] || allAnswers.filter_tier });
  }

  return (
    <div className="bg-muted rounded-2xl p-6 border border-border shadow-sm">
      <div className="h-1.5 bg-muted-foreground/20 rounded-full mb-6 overflow-hidden">
        <div className="h-full bg-primary rounded-full transition-all duration-500" style={{ width: '100%' }} />
      </div>

      <div className="flex items-center gap-2 text-primary text-sm mb-4">
        <Rocket className="h-5 w-5" />
        <span className="font-medium">Review & Run</span>
      </div>

      <h3 className="text-xl font-semibold text-foreground mb-2">Ready to start?</h3>
      <p className="text-muted-foreground text-sm mb-6">
        Review your design parameters below, then hit Start to begin the pipeline.
      </p>

      <div className="bg-card rounded-xl border border-border p-4 mb-6 space-y-2">
        {summaryItems.map((item) => (
          <div key={item.label} className="flex items-center justify-between text-sm">
            <span className="text-muted-foreground">{item.label}</span>
            <span className="font-medium text-foreground">{item.value}</span>
          </div>
        ))}
      </div>

      <div className="flex items-center justify-between">
        <button
          onClick={() => {
            if (followUpSteps.length > 0) {
              setPhase('questions');
              setCurrentStep(followUpSteps.length - 1);
              setSelectedOption(answers[followUpSteps[followUpSteps.length - 1].id] || null);
            } else {
              setPhase('design_type');
              setSelectedOption(designType);
            }
          }}
          className="text-primary text-sm font-medium hover:text-primary/80 flex items-center gap-1"
        >
          <ArrowLeft className="h-4 w-4" />
          Back
        </button>
        <button
          onClick={() => handleFinalComplete(allAnswers)}
          className="px-6 py-2.5 rounded-xl font-medium bg-primary text-primary-foreground hover:bg-primary/90 shadow-sm transition-all flex items-center gap-2"
        >
          <Rocket className="h-4 w-4" />
          Start Design
        </button>
      </div>
    </div>
  );

  // ---- Build result and call onComplete ----

  function handleFinalComplete(finalAnswers: Record<string, string>, freeTextPrompt?: string) {
    const scale = scaleToParams(finalAnswers.scale || 'standard');
    const dt = finalAnswers.design_type || 'general';

    // Build a natural language prompt from the wizard answers
    const promptParts: string[] = [];
    if (freeTextPrompt) {
      promptParts.push(freeTextPrompt);
    } else {
      if (dt === 'metal' || dt === 'scaffold') {
        const metalLabel = METAL_LABELS[finalAnswers.target_metal] || finalAnswers.target_metal || '';
        const ligand = finalAnswers.ligand_name !== 'none' ? finalAnswers.ligand_name : '';
        const goal = GOAL_LABELS[finalAnswers.design_goal] || '';
        promptParts.push(`Design a ${metalLabel.split(' (')[0]?.toLowerCase() || 'metal'}-binding protein`);
        if (ligand) promptParts.push(`with ${ligand}`);
        if (dt === 'scaffold') promptParts.push('using scaffold search');
        if (goal) promptParts.push(`for ${goal.toLowerCase()}`);
      } else if (dt === 'ligand') {
        const ligand = finalAnswers.ligand_name || 'ligand';
        promptParts.push(`Design a ${ligand} interface dimer`);
        if (finalAnswers.symmetry === 'C2') promptParts.push('with C2 symmetry');
      } else if (dt === 'binder') {
        const target = finalAnswers.target_pdb_id || '';
        promptParts.push('Design a protein binder');
        if (target && target !== 'custom') promptParts.push(`for ${target}`);
      } else {
        promptParts.push('Design a protein');
      }
    }

    const result: WizardResult = {
      design_type: dt,
      user_prompt: promptParts.join(' '),
      num_designs: scale.num_designs,
      use_sweep: scale.use_sweep,
    };

    if (finalAnswers.target_metal) result.target_metal = finalAnswers.target_metal;
    if (finalAnswers.ligand_name && finalAnswers.ligand_name !== 'none') result.ligand_name = finalAnswers.ligand_name;
    if (finalAnswers.design_goal) result.design_goal = finalAnswers.design_goal;
    if (finalAnswers.filter_tier) result.filter_tier = finalAnswers.filter_tier;
    if (finalAnswers.symmetry) result.symmetry = finalAnswers.symmetry;
    if (finalAnswers.target_pdb_id && finalAnswers.target_pdb_id !== 'custom') result.target_pdb_id = finalAnswers.target_pdb_id;
    if (finalAnswers.protocol) result.protocol = finalAnswers.protocol;

    onComplete(result);
  }
}
