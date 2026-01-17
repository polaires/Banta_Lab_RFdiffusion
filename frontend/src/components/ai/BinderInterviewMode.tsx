'use client';

import { useState } from 'react';
import { Target, FlaskConical, Ruler, Crosshair, Copy, BadgeCheck, Star, CheckCircle, ArrowLeft, ArrowRight, ClipboardList, type LucideIcon } from 'lucide-react';

// Icon mapping for step icons
const STEP_ICONS: Record<string, LucideIcon> = {
  target: Target,
  science: FlaskConical,
  straighten: Ruler,
  my_location: Crosshair,
  content_copy: Copy,
  verified: BadgeCheck,
  star: Star,
};

// Protocol types matching backend PROTOCOL_PRESETS
export type ProtocolType =
  | 'miniprotein_default'
  | 'miniprotein_hardtarget'
  | 'peptide_default'
  | 'peptide_helical'
  | 'large_binder';

// Binder preferences interface
export interface BinderPreferences {
  targetType: string;
  targetLabel: string;
  targetPdbId?: string;
  targetPdbContent?: string;  // PDB content for manual hotspot selection
  protocol: ProtocolType;
  protocolLabel: string;
  binderLength: string;
  binderLengthLabel: string;
  hotspotMethod: 'auto' | 'manual' | 'none';
  hotspotMethodLabel: string;
  manualHotspots?: string[];  // Selected hotspots when using manual method
  numDesigns: number;
  qualityThreshold: 'relaxed' | 'standard' | 'strict';
  qualityThresholdLabel: string;
  priority: 'affinity' | 'confidence' | 'balanced';
  priorityLabel: string;
  validateStructure?: boolean;  // Enable ESMFold validation
}

interface InterviewOption {
  id: string;
  label: string;
  description: string;
  pdbId?: string;
}

interface InterviewStep {
  id: string;
  question: string;
  description: string;
  options: InterviewOption[];
  icon: string;
}

const BINDER_INTERVIEW_STEPS: InterviewStep[] = [
  {
    id: 'target_type',
    question: 'What target protein do you want to bind?',
    description: 'Select the protein target for your binder design',
    icon: 'target',
    options: [
      { id: 'rubredoxin', label: 'Rubredoxin (1BRF)', description: 'Small iron-sulfur protein, ideal for testing', pdbId: '1BRF' },
      { id: 'gfp', label: 'GFP (1EMA)', description: 'Green fluorescent protein barrel', pdbId: '1EMA' },
      { id: 'lysozyme', label: 'Lysozyme (1LYZ)', description: 'Classic model protein for binding studies', pdbId: '1LYZ' },
      { id: 'custom', label: 'Upload Custom', description: 'Use your own target structure' },
    ]
  },
  {
    id: 'protocol',
    question: 'What type of binder do you want to design?',
    description: 'Protocol presets optimize settings for different binder types (BindCraft-inspired)',
    icon: 'science',
    options: [
      { id: 'miniprotein_default', label: 'Miniprotein (Recommended)', description: 'Standard 60-100 residue binder, balanced settings' },
      { id: 'miniprotein_hardtarget', label: 'Miniprotein (Hard Target)', description: 'For difficult targets - more sampling, relaxed filters' },
      { id: 'peptide_default', label: 'Peptide', description: 'Short 15-30 residue peptide binder' },
      { id: 'peptide_helical', label: 'Helical Peptide', description: 'Alpha-helical peptide, good for coiled-coil interfaces' },
      { id: 'large_binder', label: 'Large Binder', description: 'Large 100-150 residue binder for maximum interface' },
    ]
  },
  {
    id: 'binder_length',
    question: 'What binder size do you want?',
    description: 'Length of the designed binder protein',
    icon: 'straighten',
    options: [
      { id: '40-60', label: 'Small (40-60)', description: 'Compact binder, faster design' },
      { id: '60-80', label: 'Medium (60-80)', description: 'Recommended balance of size and affinity' },
      { id: '80-100', label: 'Large (80-100)', description: 'Larger interface, potentially higher affinity' },
    ]
  },
  {
    id: 'hotspot_method',
    question: 'How should I select the binding site?',
    description: 'Hotspots guide where the binder will attach to your target',
    icon: 'my_location',
    options: [
      { id: 'auto', label: 'Auto-detect (Recommended)', description: 'AI finds the best binding surface using SASA analysis' },
      { id: 'manual', label: 'Select manually', description: 'Choose specific residues on the 3D structure' },
      { id: 'none', label: 'No hotspots', description: 'Maximize overall contact (may result in wrap-around binding)' },
    ]
  },
  {
    id: 'num_designs',
    question: 'How many designs to generate?',
    description: 'More designs give better chance of finding high-affinity binders',
    icon: 'content_copy',
    options: [
      { id: '3', label: '3 designs', description: 'Quick exploration (~3 min)' },
      { id: '5', label: '5 designs (Recommended)', description: 'Good coverage (~5 min)' },
      { id: '10', label: '10 designs', description: 'Thorough search (~10 min)' },
    ]
  },
  {
    id: 'quality_threshold',
    question: 'What quality level do you need?',
    description: 'Higher thresholds filter more aggressively',
    icon: 'verified',
    options: [
      { id: 'relaxed', label: 'Relaxed', description: 'More designs pass, suitable for exploration' },
      { id: 'standard', label: 'Standard (Recommended)', description: 'Balanced filtering, good quality' },
      { id: 'strict', label: 'Strict', description: 'Only top-quality designs pass' },
    ]
  },
  {
    id: 'priority',
    question: 'What is your priority?',
    description: 'Optimize for your specific goal',
    icon: 'star',
    options: [
      { id: 'affinity', label: 'Binding Affinity', description: 'Maximize binding strength (GNINA score)' },
      { id: 'confidence', label: 'Sequence Confidence', description: 'Prioritize high-confidence sequences (ESM score)' },
      { id: 'balanced', label: 'Balanced (Recommended)', description: 'Good affinity + high confidence' },
    ]
  },
];

// Labels for display
const TARGET_LABELS: Record<string, string> = {
  'rubredoxin': 'Rubredoxin (1BRF)',
  'gfp': 'GFP (1EMA)',
  'lysozyme': 'Lysozyme (1LYZ)',
  'custom': 'Custom Target',
};

const TARGET_PDB_IDS: Record<string, string> = {
  'rubredoxin': '1BRF',
  'gfp': '1EMA',
  'lysozyme': '1LYZ',
};

const BINDER_LENGTH_LABELS: Record<string, string> = {
  '40-60': 'Small (40-60 residues)',
  '60-80': 'Medium (60-80 residues)',
  '80-100': 'Large (80-100 residues)',
};

const HOTSPOT_METHOD_LABELS: Record<string, string> = {
  'auto': 'Auto-detect',
  'manual': 'Manual selection',
  'none': 'No hotspots',
};

const QUALITY_LABELS: Record<string, string> = {
  'relaxed': 'Relaxed',
  'standard': 'Standard',
  'strict': 'Strict',
};

const PRIORITY_LABELS: Record<string, string> = {
  'affinity': 'Binding Affinity',
  'confidence': 'Sequence Confidence',
  'balanced': 'Balanced',
};

const PROTOCOL_LABELS: Record<string, string> = {
  'miniprotein_default': 'Miniprotein (Standard)',
  'miniprotein_hardtarget': 'Miniprotein (Hard Target)',
  'peptide_default': 'Peptide',
  'peptide_helical': 'Helical Peptide',
  'large_binder': 'Large Binder',
};

// Protocol-specific defaults that can inform UI
const PROTOCOL_DEFAULTS: Record<string, { binderLength: string; numDesigns: string }> = {
  'miniprotein_default': { binderLength: '60-80', numDesigns: '5' },
  'miniprotein_hardtarget': { binderLength: '80-100', numDesigns: '10' },
  'peptide_default': { binderLength: '40-60', numDesigns: '5' },
  'peptide_helical': { binderLength: '40-60', numDesigns: '5' },
  'large_binder': { binderLength: '80-100', numDesigns: '5' },
};

interface BinderInterviewModeProps {
  onComplete: (prefs: BinderPreferences) => void;
  onCancel?: () => void;
}

export function BinderInterviewMode({ onComplete, onCancel }: BinderInterviewModeProps) {
  const [currentStep, setCurrentStep] = useState(0);
  const [answers, setAnswers] = useState<Record<string, string>>({});
  const [selectedOption, setSelectedOption] = useState<string | null>(null);

  const step = BINDER_INTERVIEW_STEPS[currentStep];
  const progress = ((currentStep + 1) / BINDER_INTERVIEW_STEPS.length) * 100;
  const isLastStep = currentStep === BINDER_INTERVIEW_STEPS.length - 1;

  const handleSelect = (optionId: string) => {
    setSelectedOption(optionId);
  };

  const handleNext = () => {
    if (!selectedOption) return;

    const newAnswers = { ...answers, [step.id]: selectedOption };
    setAnswers(newAnswers);

    if (isLastStep) {
      // Convert answers to BinderPreferences
      const preferences: BinderPreferences = {
        targetType: newAnswers.target_type,
        targetLabel: TARGET_LABELS[newAnswers.target_type] || newAnswers.target_type,
        targetPdbId: TARGET_PDB_IDS[newAnswers.target_type],
        protocol: newAnswers.protocol as ProtocolType,
        protocolLabel: PROTOCOL_LABELS[newAnswers.protocol],
        binderLength: newAnswers.binder_length,
        binderLengthLabel: BINDER_LENGTH_LABELS[newAnswers.binder_length],
        hotspotMethod: newAnswers.hotspot_method as 'auto' | 'manual' | 'none',
        hotspotMethodLabel: HOTSPOT_METHOD_LABELS[newAnswers.hotspot_method],
        numDesigns: parseInt(newAnswers.num_designs),
        qualityThreshold: newAnswers.quality_threshold as 'relaxed' | 'standard' | 'strict',
        qualityThresholdLabel: QUALITY_LABELS[newAnswers.quality_threshold],
        priority: newAnswers.priority as 'affinity' | 'confidence' | 'balanced',
        priorityLabel: PRIORITY_LABELS[newAnswers.priority],
        validateStructure: true,  // Enable ESMFold validation by default
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
      setSelectedOption(answers[BINDER_INTERVIEW_STEPS[currentStep - 1].id] || null);
    }
  };

  return (
    <div className="bg-gradient-to-br from-teal-50 to-emerald-50 rounded-2xl p-6 border border-teal-200 shadow-sm">
      {/* Progress bar */}
      <div className="h-1.5 bg-teal-100 rounded-full mb-6 overflow-hidden">
        <div
          className="h-full bg-gradient-to-r from-teal-500 to-emerald-500 rounded-full transition-all duration-500"
          style={{ width: `${progress}%` }}
        />
      </div>

      {/* Step indicator */}
      <div className="flex items-center gap-2 text-teal-600 text-sm mb-4">
        {(() => {
          const IconComponent = STEP_ICONS[step.icon];
          return IconComponent ? <IconComponent className="h-5 w-5" /> : null;
        })()}
        <span className="font-medium">Step {currentStep + 1} of {BINDER_INTERVIEW_STEPS.length}</span>
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
                  ? 'border-teal-500 bg-teal-50 shadow-sm'
                  : 'border-slate-200 hover:border-teal-300 hover:bg-teal-50/50'
              }`}
            >
              <div className="flex items-center justify-between">
                <div>
                  <div className={`font-medium ${isSelected ? 'text-teal-900' : 'text-slate-900'}`}>
                    {option.label}
                  </div>
                  <div className={`text-sm ${isSelected ? 'text-teal-700' : 'text-slate-500'}`}>
                    {option.description}
                  </div>
                  {option.pdbId && (
                    <code className={`text-xs mt-1 inline-block px-2 py-0.5 rounded ${
                      isSelected ? 'bg-teal-100 text-teal-700' : 'bg-slate-100 text-slate-500'
                    }`}>
                      PDB: {option.pdbId}
                    </code>
                  )}
                </div>
                {isSelected && (
                  <CheckCircle className="h-6 w-6 text-teal-600" />
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
              className="text-teal-600 text-sm font-medium hover:text-teal-800 flex items-center gap-1"
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
              ? 'bg-gradient-to-r from-teal-600 to-emerald-600 text-white hover:from-teal-700 hover:to-emerald-700 shadow-sm'
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
