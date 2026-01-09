'use client';

import { useState } from 'react';
import { FormSection, FormField, FormRow } from './shared/FormSection';
import { PdbUploader } from './shared/PdbUploader';
import { QUALITY_PRESETS, RFD3Request, TaskFormProps } from './shared/types';

type QualityPreset = keyof typeof QUALITY_PRESETS;

const NOISE_LEVELS = [
  { value: 5, label: 'Very Low', description: 'Minor local adjustments' },
  { value: 10, label: 'Low', description: 'Small structural changes' },
  { value: 15, label: 'Medium', description: 'Moderate redesign' },
  { value: 20, label: 'High', description: 'Major structural changes' },
];

export function RefinementForm({ onSubmit, isSubmitting, health }: TaskFormProps) {
  // Required
  const [pdbContent, setPdbContent] = useState<string | null>(null);
  const [pdbFileName, setPdbFileName] = useState<string | null>(null);
  const [partialT, setPartialT] = useState(10);

  // Optional
  const [contig, setContig] = useState('');
  const [qualityPreset, setQualityPreset] = useState<QualityPreset>('Balanced');
  const [numDesigns, setNumDesigns] = useState(1);
  const [seed, setSeed] = useState<string>('');

  const handleSubmit = async () => {
    const preset = QUALITY_PRESETS[qualityPreset];

    const request: RFD3Request = {
      pdb_content: pdbContent || undefined,
      contig: contig.trim() || undefined,
      num_designs: numDesigns,
      partial_t: partialT,
      num_timesteps: preset.num_timesteps,
      step_scale: preset.step_scale,
      gamma_0: preset.gamma_0,
    };

    if (seed) {
      request.seed = parseInt(seed, 10);
    }

    await onSubmit(request);
  };

  const isValid = pdbContent !== null;

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex items-center gap-3 pb-4 border-b border-slate-200">
        <div className="w-10 h-10 rounded-lg bg-slate-200 flex items-center justify-center">
          <span className="material-symbols-outlined text-slate-600">tune</span>
        </div>
        <div>
          <h2 className="font-semibold text-slate-900">Structure Refinement</h2>
          <p className="text-sm text-slate-500">Refine an existing structure using partial diffusion</p>
        </div>
      </div>

      {/* Info Banner */}
      <div className="p-4 rounded-xl bg-slate-100 border border-slate-200">
        <div className="flex items-start gap-3">
          <span className="material-symbols-outlined text-slate-600 text-xl">lightbulb</span>
          <div className="text-sm text-slate-700">
            <strong>Partial diffusion</strong> adds controlled noise to your structure and then
            denoises it, allowing the model to explore nearby conformations while preserving
            the overall fold. Higher noise levels allow more dramatic changes.
          </div>
        </div>
      </div>

      {/* Input PDB - Required */}
      <FormSection
        title="Input Structure"
        description="Upload the PDB file you want to refine"
        required
      >
        <PdbUploader
          label="Input PDB"
          required
          value={pdbContent}
          fileName={pdbFileName}
          onChange={(content, name) => {
            setPdbContent(content);
            setPdbFileName(name);
          }}
        />
      </FormSection>

      {/* Noise Level - Required */}
      <FormSection
        title="Noise Level"
        description="How much structural change to allow during refinement"
        required
      >
        <div className="space-y-4">
          {/* Preset buttons */}
          <div className="grid grid-cols-2 sm:grid-cols-4 gap-2">
            {NOISE_LEVELS.map((level) => (
              <button
                key={level.value}
                onClick={() => setPartialT(level.value)}
                className={`p-3 rounded-xl border-2 text-left transition-all ${
                  partialT === level.value
                    ? 'border-blue-400 bg-blue-50'
                    : 'border-slate-200 hover:border-slate-300'
                }`}
              >
                <div className="font-medium text-sm text-slate-900">{level.label}</div>
                <div className="text-xs text-slate-500 mt-0.5">{level.description}</div>
              </button>
            ))}
          </div>

          {/* Custom slider */}
          <div className="px-2">
            <div className="flex justify-between text-xs text-slate-500 mb-2">
              <span>Less Change</span>
              <span className="font-medium text-slate-700">{partialT} Å</span>
              <span>More Change</span>
            </div>
            <input
              type="range"
              min={1}
              max={25}
              value={partialT}
              onChange={(e) => setPartialT(parseInt(e.target.value, 10))}
              className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-blue-600"
            />
          </div>

          {/* Visual indicator */}
          <div className="flex gap-1">
            {[...Array(25)].map((_, i) => (
              <div
                key={i}
                className={`flex-1 h-2 rounded-full transition-all ${
                  i < partialT
                    ? i < 5
                      ? 'bg-green-400'
                      : i < 10
                        ? 'bg-yellow-400'
                        : i < 15
                          ? 'bg-orange-400'
                          : 'bg-red-400'
                    : 'bg-slate-200'
                }`}
              />
            ))}
          </div>

          <div className="p-3 rounded-xl bg-slate-50 border border-slate-200 text-sm text-slate-600">
            <strong>Tip:</strong> Start with lower noise (5-10 Å) for subtle refinements.
            Use higher noise (15-20 Å) only if you want significant structural changes.
          </div>
        </div>
      </FormSection>

      {/* Optional: Contig specification */}
      <FormSection
        title="Region Selection (Optional)"
        description="Specify which regions to refine using contig syntax"
      >
        <FormField
          label="Contig String"
          hint="Leave empty to refine entire structure. Use contig syntax for partial refinement."
        >
          <input
            type="text"
            value={contig}
            onChange={(e) => setContig(e.target.value)}
            placeholder="e.g., A1-50,0,A60-100 (diffuse middle region)"
            className="w-full px-4 py-2.5 rounded-xl border border-slate-200 focus:border-blue-400 focus:ring-2 focus:ring-blue-100 outline-none transition-all font-mono text-sm"
          />
        </FormField>
      </FormSection>

      {/* Generation Settings */}
      <FormSection title="Generation Settings">
        <FormRow>
          <FormField label="# Designs">
            <input
              type="number"
              value={numDesigns}
              onChange={(e) => setNumDesigns(Math.max(1, parseInt(e.target.value) || 1))}
              min={1}
              max={10}
              className="w-full px-4 py-2.5 rounded-xl border border-slate-200 focus:border-blue-400 focus:ring-2 focus:ring-blue-100 outline-none transition-all"
            />
          </FormField>
          <FormField label="Random Seed" hint="For reproducible results">
            <input
              type="number"
              value={seed}
              onChange={(e) => setSeed(e.target.value)}
              placeholder="Optional"
              className="w-full px-4 py-2.5 rounded-xl border border-slate-200 focus:border-blue-400 focus:ring-2 focus:ring-blue-100 outline-none transition-all"
            />
          </FormField>
        </FormRow>
      </FormSection>

      {/* Quality Settings */}
      <FormSection
        title="Quality Settings"
        description="Higher quality takes longer but produces better refinements"
      >
        <div className="grid grid-cols-2 sm:grid-cols-4 gap-2">
          {(Object.keys(QUALITY_PRESETS) as QualityPreset[]).filter(
            (preset) => preset !== 'Binder Optimized'
          ).map((preset) => (
            <button
              key={preset}
              onClick={() => setQualityPreset(preset)}
              className={`p-3 rounded-xl border-2 text-left transition-all ${
                qualityPreset === preset
                  ? 'border-blue-400 bg-blue-50'
                  : 'border-slate-200 hover:border-slate-300'
              }`}
            >
              <div className="font-medium text-sm text-slate-900">{preset}</div>
              <div className="text-xs text-slate-500 mt-0.5">
                {QUALITY_PRESETS[preset].description}
              </div>
            </button>
          ))}
        </div>
      </FormSection>

      {/* Submit Button */}
      <div className="pt-4 border-t border-slate-200">
        <button
          onClick={handleSubmit}
          disabled={!isValid || isSubmitting || !health}
          className={`w-full py-3 px-6 rounded-xl font-semibold text-white transition-all flex items-center justify-center gap-2 ${
            isValid && !isSubmitting && !!health
              ? 'bg-blue-600 hover:bg-blue-700 shadow-lg shadow-blue-600/20'
              : 'bg-slate-300 cursor-not-allowed'
          }`}
        >
          {isSubmitting ? (
            <>
              <span className="material-symbols-outlined animate-spin">progress_activity</span>
              Submitting...
            </>
          ) : (
            <>
              <span className="material-symbols-outlined">rocket_launch</span>
              Refine {numDesigns} Structure{numDesigns > 1 ? 's' : ''}
            </>
          )}
        </button>
        {!health && (
          <p className="text-center text-sm text-amber-600 mt-2">
            Backend service unavailable. Please check connection.
          </p>
        )}
      </div>
    </div>
  );
}
