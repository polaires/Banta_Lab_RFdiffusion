'use client';

import { useState } from 'react';
import { FlaskConical, Ruler, Copy, Star, CheckCircle, ArrowLeft, ArrowRight, ClipboardList, type LucideIcon } from 'lucide-react';

// Icon mapping for dynamic icons
const STEP_ICONS: Record<string, LucideIcon> = {
  science: FlaskConical,
  straighten: Ruler,
  content_copy: Copy,
  star: Star,
};

// Ligand preferences interface
export interface LigandPreferences {
  ligandType: string;
  ligandLabel: string;
  ligandSmiles: string;
  chainLength: string;
  chainLengthLabel: string;
  numDesigns: number;
  priority: 'affinity' | 'separability' | 'balanced';
  priorityLabel: string;
}

interface InterviewOption {
  id: string;
  label: string;
  description: string;
  smiles?: string;
}

interface InterviewStep {
  id: string;
  question: string;
  description: string;
  options: InterviewOption[];
  icon: string;
}

const LIGAND_INTERVIEW_STEPS: InterviewStep[] = [
  {
    id: 'ligand_type',
    question: 'Which ligand do you want at the interface?',
    description: 'Select the molecule that will sit between the two protein chains',
    icon: 'science',
    options: [
      { id: 'azobenzene', label: 'Azobenzene', description: 'Light-switchable, symmetric diazo compound', smiles: 'c1ccc(cc1)N=Nc2ccccc2' },
      { id: 'caffeine', label: 'Caffeine', description: 'Common stimulant alkaloid', smiles: 'Cn1cnc2c1c(=O)n(c(=O)n2C)C' },
      { id: 'biphenyl', label: 'Biphenyl', description: 'Two connected phenyl rings', smiles: 'c1ccc(cc1)c2ccccc2' },
      { id: 'naphthalene', label: 'Naphthalene', description: 'Fused bicyclic aromatic', smiles: 'c1ccc2ccccc2c1' },
    ]
  },
  {
    id: 'chain_length',
    question: 'What chain length do you want?',
    description: 'Length of each protein chain in residues',
    icon: 'straighten',
    options: [
      { id: '40-60', label: 'Short (40-60)', description: 'Smaller proteins, faster design' },
      { id: '60-80', label: 'Medium (60-80)', description: 'Recommended for most ligands' },
      { id: '80-100', label: 'Long (80-100)', description: 'Larger binding surface' },
    ]
  },
  {
    id: 'num_designs',
    question: 'How many designs to generate?',
    description: 'More designs = better chance of finding optimal solution',
    icon: 'content_copy',
    options: [
      { id: '3', label: '3 designs', description: 'Quick exploration (~5 min)' },
      { id: '5', label: '5 designs (Recommended)', description: 'Balanced approach (~10 min)' },
      { id: '10', label: '10 designs', description: 'Thorough search (~20 min)' },
    ]
  },
  {
    id: 'priority',
    question: 'What is your priority?',
    description: 'Optimize for your specific goal',
    icon: 'star',
    options: [
      { id: 'affinity', label: 'Binding Affinity', description: 'Maximize ligand binding strength (kcal/mol)' },
      { id: 'separability', label: 'Separability', description: 'Ensure chains can physically dissociate' },
      { id: 'balanced', label: 'Balanced (Recommended)', description: 'Good affinity + guaranteed separability' },
    ]
  },
];

// Labels for display
const LIGAND_LABELS: Record<string, string> = {
  'azobenzene': 'Azobenzene',
  'caffeine': 'Caffeine',
  'biphenyl': 'Biphenyl',
  'naphthalene': 'Naphthalene',
};

const LIGAND_SMILES: Record<string, string> = {
  'azobenzene': 'c1ccc(cc1)N=Nc2ccccc2',
  'caffeine': 'Cn1cnc2c1c(=O)n(c(=O)n2C)C',
  'biphenyl': 'c1ccc(cc1)c2ccccc2',
  'naphthalene': 'c1ccc2ccccc2c1',
};

const CHAIN_LENGTH_LABELS: Record<string, string> = {
  '40-60': 'Short (40-60 residues)',
  '60-80': 'Medium (60-80 residues)',
  '80-100': 'Long (80-100 residues)',
};

const PRIORITY_LABELS: Record<string, string> = {
  'affinity': 'Binding Affinity',
  'separability': 'Separability',
  'balanced': 'Balanced',
};

interface LigandInterviewModeProps {
  onComplete: (prefs: LigandPreferences) => void;
  onCancel?: () => void;
}

export function LigandInterviewMode({ onComplete, onCancel }: LigandInterviewModeProps) {
  const [currentStep, setCurrentStep] = useState(0);
  const [answers, setAnswers] = useState<Record<string, string>>({});
  const [selectedOption, setSelectedOption] = useState<string | null>(null);

  const step = LIGAND_INTERVIEW_STEPS[currentStep];
  const progress = ((currentStep + 1) / LIGAND_INTERVIEW_STEPS.length) * 100;
  const isLastStep = currentStep === LIGAND_INTERVIEW_STEPS.length - 1;

  const handleSelect = (optionId: string) => {
    setSelectedOption(optionId);
  };

  const handleNext = () => {
    if (!selectedOption) return;

    const newAnswers = { ...answers, [step.id]: selectedOption };
    setAnswers(newAnswers);

    if (isLastStep) {
      // Convert answers to LigandPreferences
      const preferences: LigandPreferences = {
        ligandType: newAnswers.ligand_type,
        ligandLabel: LIGAND_LABELS[newAnswers.ligand_type] || newAnswers.ligand_type,
        ligandSmiles: LIGAND_SMILES[newAnswers.ligand_type] || newAnswers.ligand_type,
        chainLength: newAnswers.chain_length,
        chainLengthLabel: CHAIN_LENGTH_LABELS[newAnswers.chain_length],
        numDesigns: parseInt(newAnswers.num_designs),
        priority: newAnswers.priority as 'affinity' | 'separability' | 'balanced',
        priorityLabel: PRIORITY_LABELS[newAnswers.priority],
      };
      onComplete(preferences);
    } else {
      setCurrentStep(prev => prev + 1);
      setSelectedOption(null);
    }
  };

  const handleBack = () => {
    if (currentStep > 0) {
      setCurrentStep(prev => prev - 1);
      setSelectedOption(answers[LIGAND_INTERVIEW_STEPS[currentStep - 1].id] || null);
    }
  };

  return (
    <div className="bg-gradient-to-br from-purple-50 to-pink-50 rounded-2xl p-6 border border-purple-200 shadow-sm">
      {/* Progress bar */}
      <div className="h-1.5 bg-purple-100 rounded-full mb-6 overflow-hidden">
        <div
          className="h-full bg-gradient-to-r from-purple-500 to-pink-500 rounded-full transition-all duration-500"
          style={{ width: `${progress}%` }}
        />
      </div>

      {/* Step indicator */}
      <div className="flex items-center gap-2 text-purple-600 text-sm mb-4">
        {(() => {
          const IconComponent = STEP_ICONS[step.icon];
          return IconComponent ? <IconComponent className="h-5 w-5" /> : null;
        })()}
        <span className="font-medium">Step {currentStep + 1} of {LIGAND_INTERVIEW_STEPS.length}</span>
      </div>

      {/* Question */}
      <h3 className="text-xl font-semibold text-slate-900 mb-2">{step.question}</h3>
      <p className="text-slate-600 text-sm mb-6">{step.description}</p>

      {/* Options */}
      <div className="space-y-3 mb-6">
        {step.options.map(option => {
          const isSelected = selectedOption === option.id;
          return (
            <button
              key={option.id}
              onClick={() => handleSelect(option.id)}
              className={`w-full text-left p-4 rounded-xl border-2 transition-all ${
                isSelected
                  ? 'border-purple-500 bg-purple-50 shadow-sm'
                  : 'border-slate-200 hover:border-purple-300 hover:bg-purple-50/50'
              }`}
            >
              <div className="flex items-center justify-between">
                <div>
                  <div className={`font-medium ${isSelected ? 'text-purple-900' : 'text-slate-900'}`}>
                    {option.label}
                  </div>
                  <div className={`text-sm ${isSelected ? 'text-purple-700' : 'text-slate-500'}`}>
                    {option.description}
                  </div>
                  {option.smiles && (
                    <code className={`text-xs mt-1 inline-block px-2 py-0.5 rounded ${
                      isSelected ? 'bg-purple-100 text-purple-700' : 'bg-slate-100 text-slate-500'
                    }`}>
                      {option.smiles}
                    </code>
                  )}
                </div>
                {isSelected && (
                  <CheckCircle className="h-6 w-6 text-purple-600" />
                )}
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
              className="text-purple-600 text-sm font-medium hover:text-purple-800 flex items-center gap-1"
            >
              <ArrowLeft className="h-4 w-4" />
              Back
            </button>
          ) : onCancel ? (
            <button
              onClick={onCancel}
              className="text-slate-500 text-sm font-medium hover:text-slate-700"
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
              ? 'bg-gradient-to-r from-purple-600 to-pink-600 text-white hover:from-purple-700 hover:to-pink-700 shadow-sm'
              : 'bg-slate-200 text-slate-400 cursor-not-allowed'
          }`}
        >
          {isLastStep ? 'Review Design' : 'Next'}
          {isLastStep ? <ClipboardList className="h-4 w-4" /> : <ArrowRight className="h-4 w-4" />}
        </button>
      </div>
    </div>
  );
}
