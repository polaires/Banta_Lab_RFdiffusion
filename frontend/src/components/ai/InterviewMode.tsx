'use client';

import { useState } from 'react';
import { Gem, SlidersHorizontal, Network, Copy, Star, CheckCircle, ArrowLeft, ArrowRight, ClipboardList, type LucideIcon } from 'lucide-react';
import type { UserPreferences } from '@/lib/api';

// Icon mapping for step icons
const STEP_ICONS: Record<string, LucideIcon> = {
  diamond: Gem,
  tune: SlidersHorizontal,
  hub: Network,
  content_copy: Copy,
  star: Star,
};

interface InterviewOption {
  id: string;
  label: string;
  description: string;
}

interface InterviewStep {
  id: string;
  question: string;
  description: string;
  options: InterviewOption[];
  icon: string;
}

const INTERVIEW_STEPS: InterviewStep[] = [
  {
    id: 'target_metal',
    question: 'What metal do you want to bind?',
    description: 'Select your target metal ion for the redesigned binding site',
    icon: 'diamond',
    options: [
      { id: 'TB', label: 'Terbium (Tb)', description: 'Lanthanide with 8-9 coordination, luminescent properties' },
      { id: 'GD', label: 'Gadolinium (Gd)', description: 'Lanthanide, used in MRI contrast agents' },
      { id: 'EU', label: 'Europium (Eu)', description: 'Lanthanide with strong red fluorescence' },
      { id: 'CA', label: 'Calcium (Ca)', description: 'Common signaling ion, 6-8 coordination' },
      { id: 'ZN', label: 'Zinc (Zn)', description: 'Catalytic metal, 4-6 coordination' },
    ]
  },
  {
    id: 'design_aggressiveness',
    question: 'How much structural change is acceptable?',
    description: 'More aggressive designs allow larger modifications but may affect stability',
    icon: 'tune',
    options: [
      { id: 'conservative', label: 'Conservative', description: 'Minimal changes, preserve most of the original structure' },
      { id: 'moderate', label: 'Moderate (Recommended)', description: 'Balance between innovation and stability' },
      { id: 'aggressive', label: 'Aggressive', description: 'Allow significant redesign for optimal metal binding' }
    ]
  },
  {
    id: 'coordination_preference',
    question: 'What coordination geometry do you prefer?',
    description: 'Higher coordination numbers are typical for lanthanides',
    icon: 'hub',
    options: [
      { id: 'tetrahedral', label: 'Tetrahedral (4)', description: 'Common for Zn, Fe - compact site' },
      { id: 'octahedral', label: 'Octahedral (6)', description: 'Common for transition metals' },
      { id: 'high', label: 'High Coordination (8-9)', description: 'Ideal for lanthanides - larger binding pocket' },
      { id: 'auto', label: 'Let AI decide', description: 'Optimal geometry based on target metal' }
    ]
  },
  {
    id: 'num_designs',
    question: 'How many design variants to generate?',
    description: 'More variants increase chances of finding an optimal solution',
    icon: 'content_copy',
    options: [
      { id: '3', label: '3 designs', description: 'Quick exploration' },
      { id: '5', label: '5 designs (Recommended)', description: 'Balanced approach' },
      { id: '10', label: '10 designs', description: 'Thorough search' }
    ]
  },
  {
    id: 'priority',
    question: 'What is your top priority?',
    description: 'This helps the AI optimize for your specific needs',
    icon: 'star',
    options: [
      { id: 'stability', label: 'Protein Stability', description: 'Maximize fold confidence (pLDDT)' },
      { id: 'binding', label: 'Metal Binding Affinity', description: 'Optimize coordination geometry' },
      { id: 'expression', label: 'Expression/Solubility', description: 'Favor soluble residues' },
      { id: 'balanced', label: 'Balanced (Recommended)', description: 'Equal weight to all factors' }
    ]
  }
];

// Labels for display
const METAL_LABELS: Record<string, string> = {
  'TB': 'Terbium (Tb)',
  'GD': 'Gadolinium (Gd)',
  'EU': 'Europium (Eu)',
  'CA': 'Calcium (Ca)',
  'ZN': 'Zinc (Zn)',
};

const AGGRESSIVENESS_LABELS: Record<string, string> = {
  'conservative': 'Conservative',
  'moderate': 'Moderate',
  'aggressive': 'Aggressive',
};

const COORDINATION_LABELS: Record<string, string> = {
  'tetrahedral': 'Tetrahedral (4)',
  'octahedral': 'Octahedral (6)',
  'high': 'High Coordination (8-9)',
  'auto': 'AI-selected',
};

const PRIORITY_LABELS: Record<string, string> = {
  'stability': 'Protein Stability',
  'binding': 'Metal Binding',
  'expression': 'Expression/Solubility',
  'balanced': 'Balanced',
};

interface InterviewModeProps {
  onComplete: (prefs: UserPreferences) => void;
  onCancel?: () => void;
}

export function InterviewMode({ onComplete, onCancel }: InterviewModeProps) {
  const [currentStep, setCurrentStep] = useState(0);
  const [answers, setAnswers] = useState<Record<string, string>>({});
  const [selectedOption, setSelectedOption] = useState<string | null>(null);

  const step = INTERVIEW_STEPS[currentStep];
  const progress = ((currentStep + 1) / INTERVIEW_STEPS.length) * 100;
  const isLastStep = currentStep === INTERVIEW_STEPS.length - 1;

  const handleSelect = (optionId: string) => {
    setSelectedOption(optionId);
  };

  const handleNext = () => {
    if (!selectedOption) return;

    const newAnswers = { ...answers, [step.id]: selectedOption };
    setAnswers(newAnswers);

    if (isLastStep) {
      // Convert answers to UserPreferences
      const preferences: UserPreferences = {
        targetMetal: newAnswers.target_metal,
        targetMetalLabel: METAL_LABELS[newAnswers.target_metal] || newAnswers.target_metal,
        designAggressiveness: newAnswers.design_aggressiveness as 'conservative' | 'moderate' | 'aggressive',
        aggressivenessLabel: AGGRESSIVENESS_LABELS[newAnswers.design_aggressiveness],
        coordinationPreference: newAnswers.coordination_preference as 'tetrahedral' | 'octahedral' | 'high' | 'auto',
        coordinationLabel: COORDINATION_LABELS[newAnswers.coordination_preference],
        numDesigns: parseInt(newAnswers.num_designs),
        priority: newAnswers.priority as 'stability' | 'binding' | 'expression' | 'balanced',
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
      setSelectedOption(answers[INTERVIEW_STEPS[currentStep - 1].id] || null);
    }
  };

  return (
    <div className="bg-gradient-to-br from-violet-50 to-purple-50 rounded-2xl p-6 border border-violet-200 shadow-sm">
      {/* Progress bar */}
      <div className="h-1.5 bg-violet-100 rounded-full mb-6 overflow-hidden">
        <div
          className="h-full bg-gradient-to-r from-violet-500 to-purple-600 rounded-full transition-all duration-500"
          style={{ width: `${progress}%` }}
        />
      </div>

      {/* Step indicator */}
      <div className="flex items-center gap-2 text-violet-600 text-sm mb-4">
        {(() => {
          const IconComponent = STEP_ICONS[step.icon];
          return IconComponent ? <IconComponent className="h-5 w-5" /> : null;
        })()}
        <span className="font-medium">Step {currentStep + 1} of {INTERVIEW_STEPS.length}</span>
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
                  ? 'border-violet-500 bg-violet-50 shadow-sm'
                  : 'border-slate-200 hover:border-violet-300 hover:bg-violet-50/50'
              }`}
            >
              <div className="flex items-center justify-between">
                <div>
                  <div className={`font-medium ${isSelected ? 'text-violet-900' : 'text-slate-900'}`}>
                    {option.label}
                  </div>
                  <div className={`text-sm ${isSelected ? 'text-violet-700' : 'text-slate-500'}`}>
                    {option.description}
                  </div>
                </div>
                {isSelected && (
                  <CheckCircle className="h-6 w-6 text-violet-600" />
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
              className="text-violet-600 text-sm font-medium hover:text-violet-800 flex items-center gap-1"
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
              ? 'bg-gradient-to-r from-violet-600 to-purple-600 text-white hover:from-violet-700 hover:to-purple-700 shadow-sm'
              : 'bg-slate-200 text-slate-400 cursor-not-allowed'
          }`}
        >
          {isLastStep ? 'Review Preferences' : 'Next'}
          {isLastStep ? <ClipboardList className="h-4 w-4" /> : <ArrowRight className="h-4 w-4" />}
        </button>
      </div>
    </div>
  );
}
