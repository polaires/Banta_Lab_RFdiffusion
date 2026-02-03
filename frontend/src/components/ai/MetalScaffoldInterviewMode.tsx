'use client';

import { useState } from 'react';
import {
  Gem,
  Hexagon,
  Ruler,
  Layers,
  Copy,
  CheckCircle,
  ArrowLeft,
  ArrowRight,
  ClipboardList,
  type LucideIcon,
} from 'lucide-react';

// Icon mapping for step icons
const STEP_ICONS: Record<string, LucideIcon> = {
  diamond: Gem,
  hexagon: Hexagon,
  ruler: Ruler,
  layers: Layers,
  copy: Copy,
};

// Metal Scaffold Preferences interface
export interface MetalScaffoldPreferences {
  // Metal configuration
  metal: string;
  metalLabel: string;
  // Ligand configuration
  ligand: string;
  ligandLabel: string;
  // Scaffold size
  scaffoldSize: 'small' | 'medium' | 'large' | 'auto';
  scaffoldSizeLabel: string;
  contigRange: string;
  // Optimization mode
  optimizationMode: 'sweep' | 'quick' | 'production';
  optimizationModeLabel: string;
  // Number of designs
  numDesigns: number;
  designsPerConfig: number;
}

interface InterviewOption {
  id: string;
  label: string;
  description: string;
  contigRange?: string;
}

interface InterviewStep {
  id: string;
  question: string;
  description: string;
  options: InterviewOption[];
  icon: string;
}

const METAL_SCAFFOLD_INTERVIEW_STEPS: InterviewStep[] = [
  {
    id: 'metal',
    question: 'What metal do you want to bind?',
    description: 'Select your target metal ion for the new scaffold',
    icon: 'diamond',
    options: [
      { id: 'TB', label: 'Terbium (Tb3+)', description: 'Lanthanide, CN=8-9, luminescent properties' },
      { id: 'EU', label: 'Europium (Eu3+)', description: 'Lanthanide, red fluorescence for imaging' },
      { id: 'GD', label: 'Gadolinium (Gd3+)', description: 'Lanthanide, MRI contrast applications' },
      { id: 'CA', label: 'Calcium (Ca2+)', description: 'CN=6-8, signaling and structural roles' },
      { id: 'ZN', label: 'Zinc (Zn2+)', description: 'CN=4-6, catalytic metalloenzymes' },
    ],
  },
  {
    id: 'ligand',
    question: 'What ligand is coordinating the metal?',
    description: 'The organic molecule that helps coordinate the metal',
    icon: 'hexagon',
    options: [
      { id: 'CIT', label: 'Citrate', description: 'Tridentate carboxylate, ideal for lanthanides' },
      { id: 'PQQ', label: 'PQQ', description: 'Pyrroloquinoline quinone, redox cofactor' },
      { id: 'none', label: 'None (metal only)', description: 'Design scaffold around bare metal ion' },
    ],
  },
  {
    id: 'scaffold_size',
    question: 'What scaffold size do you prefer?',
    description: 'Larger scaffolds have more room but are harder to express',
    icon: 'ruler',
    options: [
      { id: 'small', label: 'Small (100-120 aa)', description: 'Compact, easier expression', contigRange: '100-120' },
      { id: 'medium', label: 'Medium (110-130 aa)', description: 'Balanced size and stability', contigRange: '110-130' },
      { id: 'large', label: 'Large (130-150 aa)', description: 'More room for binding pocket', contigRange: '130-150' },
      { id: 'auto', label: 'Let sweep decide', description: 'Test all sizes in parameter sweep', contigRange: '100-150' },
    ],
  },
  {
    id: 'optimization_mode',
    question: 'How should we optimize the design?',
    description: 'Choose between quick exploration or thorough optimization',
    icon: 'layers',
    options: [
      { id: 'sweep', label: 'Parameter Sweep (Recommended)', description: 'Test 9 configs (3 sizes × 3 CFG values), find best' },
      { id: 'quick', label: 'Quick Exploration', description: 'Single config, 10 designs, fast results' },
      { id: 'production', label: 'Production Run', description: 'Use best config from previous sweep, 100+ designs' },
    ],
  },
  {
    id: 'num_designs',
    question: 'How many designs per configuration?',
    description: 'More designs increase chances of finding optimal solutions',
    icon: 'copy',
    options: [
      { id: '5', label: '5 designs', description: 'Quick test (for sweep: 45 total)' },
      { id: '10', label: '10 designs (Recommended)', description: 'Balanced (for sweep: 90 total)' },
      { id: '25', label: '25 designs', description: 'Thorough (for sweep: 225 total)' },
      { id: '50', label: '50 designs (Production)', description: 'Full production run (for sweep: 450 total)' },
    ],
  },
];

// Labels for display
const METAL_LABELS: Record<string, string> = {
  TB: 'Terbium (Tb3+)',
  EU: 'Europium (Eu3+)',
  GD: 'Gadolinium (Gd3+)',
  CA: 'Calcium (Ca2+)',
  ZN: 'Zinc (Zn2+)',
};

const LIGAND_LABELS: Record<string, string> = {
  CIT: 'Citrate',
  PQQ: 'PQQ',
  none: 'None (metal only)',
};

const SCAFFOLD_SIZE_LABELS: Record<string, string> = {
  small: 'Small (100-120 aa)',
  medium: 'Medium (110-130 aa)',
  large: 'Large (130-150 aa)',
  auto: 'Auto (sweep all sizes)',
};

const CONTIG_RANGES: Record<string, string> = {
  small: '100-120',
  medium: '110-130',
  large: '130-150',
  auto: '100-150',
};

const OPTIMIZATION_MODE_LABELS: Record<string, string> = {
  sweep: 'Parameter Sweep',
  quick: 'Quick Exploration',
  production: 'Production Run',
};

interface MetalScaffoldInterviewModeProps {
  onComplete: (prefs: MetalScaffoldPreferences) => void;
  onCancel?: () => void;
}

export function MetalScaffoldInterviewMode({ onComplete, onCancel }: MetalScaffoldInterviewModeProps) {
  const [currentStep, setCurrentStep] = useState(0);
  const [answers, setAnswers] = useState<Record<string, string>>({});
  const [selectedOption, setSelectedOption] = useState<string | null>(null);

  const step = METAL_SCAFFOLD_INTERVIEW_STEPS[currentStep];
  const progress = ((currentStep + 1) / METAL_SCAFFOLD_INTERVIEW_STEPS.length) * 100;
  const isLastStep = currentStep === METAL_SCAFFOLD_INTERVIEW_STEPS.length - 1;

  const handleSelect = (optionId: string) => {
    setSelectedOption(optionId);
  };

  const handleNext = () => {
    if (!selectedOption) return;

    const newAnswers = { ...answers, [step.id]: selectedOption };
    setAnswers(newAnswers);

    if (isLastStep) {
      // Convert answers to MetalScaffoldPreferences
      const scaffoldSize = newAnswers.scaffold_size as 'small' | 'medium' | 'large' | 'auto';
      const optimizationMode = newAnswers.optimization_mode as 'sweep' | 'quick' | 'production';
      const numDesigns = parseInt(newAnswers.num_designs);

      const preferences: MetalScaffoldPreferences = {
        metal: newAnswers.metal,
        metalLabel: METAL_LABELS[newAnswers.metal] || newAnswers.metal,
        ligand: newAnswers.ligand,
        ligandLabel: LIGAND_LABELS[newAnswers.ligand] || newAnswers.ligand,
        scaffoldSize,
        scaffoldSizeLabel: SCAFFOLD_SIZE_LABELS[scaffoldSize],
        contigRange: CONTIG_RANGES[scaffoldSize],
        optimizationMode,
        optimizationModeLabel: OPTIMIZATION_MODE_LABELS[optimizationMode],
        numDesigns: optimizationMode === 'sweep' ? numDesigns * 9 : numDesigns,
        designsPerConfig: numDesigns,
      };
      onComplete(preferences);
    } else {
      setCurrentStep((prev) => prev + 1);
      setSelectedOption(null);
    }
  };

  const handleBack = () => {
    if (currentStep > 0) {
      setCurrentStep((prev) => prev - 1);
      setSelectedOption(answers[METAL_SCAFFOLD_INTERVIEW_STEPS[currentStep - 1].id] || null);
    }
  };

  return (
    <div className="bg-muted rounded-2xl p-6 border border-border shadow-sm">
      {/* Progress bar */}
      <div className="h-1.5 bg-muted-foreground/20 rounded-full mb-6 overflow-hidden">
        <div
          className="h-full bg-primary rounded-full transition-all duration-500"
          style={{ width: `${progress}%` }}
        />
      </div>

      {/* Step indicator */}
      <div className="flex items-center gap-2 text-primary text-sm mb-4">
        {(() => {
          const IconComponent = STEP_ICONS[step.icon];
          return IconComponent ? <IconComponent className="h-5 w-5" /> : null;
        })()}
        <span className="font-medium">
          Step {currentStep + 1} of {METAL_SCAFFOLD_INTERVIEW_STEPS.length}
        </span>
      </div>

      {/* Question */}
      <h3 className="text-xl font-semibold text-foreground mb-2">{step.question}</h3>
      <p className="text-muted-foreground text-sm mb-6">{step.description}</p>

      {/* Options */}
      <div className="space-y-3 mb-6">
        {step.options.map((option) => {
          const isSelected = selectedOption === option.id;
          return (
            <button
              key={option.id}
              onClick={() => handleSelect(option.id)}
              className={`w-full text-left p-4 rounded-xl border-2 transition-all ${
                isSelected
                  ? 'border-primary bg-primary/5 shadow-sm'
                  : 'border-border hover:border-primary/50 hover:bg-muted/50'
              }`}
            >
              <div className="flex items-center justify-between">
                <div>
                  <div className={`font-medium ${isSelected ? 'text-foreground' : 'text-foreground'}`}>
                    {option.label}
                  </div>
                  <div className={`text-sm ${isSelected ? 'text-muted-foreground' : 'text-muted-foreground'}`}>
                    {option.description}
                  </div>
                  {option.contigRange && (
                    <code
                      className={`text-xs mt-1 inline-block px-2 py-0.5 rounded ${
                        isSelected ? 'bg-primary/10 text-primary' : 'bg-muted text-muted-foreground'
                      }`}
                    >
                      {option.contigRange} residues
                    </code>
                  )}
                </div>
                {isSelected && <CheckCircle className="h-6 w-6 text-primary" />}
              </div>
            </button>
          );
        })}
      </div>

      {/* Navigation buttons */}
      <div className="flex items-center justify-between">
        <div>
          {currentStep > 0 ? (
            <button
              onClick={handleBack}
              className="text-primary text-sm font-medium hover:text-primary/80 flex items-center gap-1"
            >
              <ArrowLeft className="h-4 w-4" />
              Back
            </button>
          ) : onCancel ? (
            <button
              onClick={onCancel}
              className="text-muted-foreground text-sm font-medium hover:text-foreground"
            >
              Cancel
            </button>
          ) : null}
        </div>
        <button
          onClick={handleNext}
          disabled={!selectedOption}
          className={`px-6 py-2.5 rounded-xl font-medium transition-all flex items-center gap-2 ${
            selectedOption
              ? 'bg-primary text-primary-foreground hover:bg-primary/90 shadow-sm'
              : 'bg-muted text-muted-foreground cursor-not-allowed'
          }`}
        >
          {isLastStep ? 'Review Configuration' : 'Next'}
          {isLastStep ? <ClipboardList className="h-4 w-4" /> : <ArrowRight className="h-4 w-4" />}
        </button>
      </div>
    </div>
  );
}
