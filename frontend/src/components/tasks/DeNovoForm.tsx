'use client';

import { useState } from 'react';
import { FormSection, FormField, FormRow } from './shared/FormSection';
import {
  QUALITY_PRESETS,
  SYMMETRY_OPTIONS,
  RFD3Request,
  TaskFormProps,
} from './shared/types';

type QualityPreset = keyof typeof QUALITY_PRESETS;

export function DeNovoForm({ onSubmit, isSubmitting, health }: TaskFormProps) {
  // Required
  const [length, setLength] = useState('100');

  // Optional
  const [isNonLoopy, setIsNonLoopy] = useState(true);
  const [qualityPreset, setQualityPreset] = useState<QualityPreset>('Balanced');
  const [symmetry, setSymmetry] = useState<string>('');
  const [numDesigns, setNumDesigns] = useState(1);
  const [seed, setSeed] = useState<string>('');

  const handleSubmit = async () => {
    const preset = QUALITY_PRESETS[qualityPreset];

    const request: RFD3Request = {
      length,
      num_designs: numDesigns,
      is_non_loopy: isNonLoopy,
      num_timesteps: preset.num_timesteps,
      step_scale: preset.step_scale,
      gamma_0: preset.gamma_0,
    };

    if (symmetry) {
      request.symmetry = { id: symmetry };
    }

    if (seed) {
      request.seed = parseInt(seed, 10);
    }

    await onSubmit(request);
  };

  const isValid = length.trim() !== '';

  // Group symmetry options by category
  const symmetryByCategory = SYMMETRY_OPTIONS.reduce((acc, opt) => {
    if (!acc[opt.category]) acc[opt.category] = [];
    acc[opt.category].push(opt);
    return acc;
  }, {} as Record<string, typeof SYMMETRY_OPTIONS[number][]>);

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex items-center gap-3 pb-4 border-b border-slate-200">
        <div className="w-9 h-9 rounded-lg bg-slate-100 flex items-center justify-center">
          <span className="material-symbols-outlined text-slate-600">add_circle</span>
        </div>
        <div>
          <h2 className="font-semibold text-slate-900">De Novo Protein Design</h2>
          <p className="text-sm text-slate-500">Generate a new protein structure from scratch</p>
        </div>
      </div>

      {/* Length Section - Required */}
      <FormSection
        title="Protein Length"
        description="Specify the number of residues. Use a range (e.g., 80-120) for length sampling."
        required
      >
        <div className="flex gap-4 items-start">
          <FormField label="Length" required className="flex-1">
            <input
              type="text"
              value={length}
              onChange={(e) => setLength(e.target.value)}
              placeholder="100 or 80-120"
              className="w-full px-3 py-2 rounded-lg border border-slate-200 focus:border-blue-400 focus:ring-2 focus:ring-blue-100 outline-none transition-all text-sm"
            />
          </FormField>
          <FormField label="# Designs" className="w-28">
            <input
              type="number"
              value={numDesigns}
              onChange={(e) => setNumDesigns(Math.max(1, parseInt(e.target.value) || 1))}
              min={1}
              max={10}
              className="w-full px-3 py-2 rounded-lg border border-slate-200 focus:border-blue-400 focus:ring-2 focus:ring-blue-100 outline-none transition-all text-sm"
            />
          </FormField>
        </div>
      </FormSection>

      {/* Quality Preset */}
      <FormSection
        title="Quality Settings"
        description="Higher quality takes longer but produces better designs"
      >
        <div className="grid grid-cols-2 sm:grid-cols-4 gap-2">
          {(Object.keys(QUALITY_PRESETS) as QualityPreset[]).map((preset) => (
            <button
              key={preset}
              onClick={() => setQualityPreset(preset)}
              className={`p-3 rounded-lg border text-left transition-all ${
                qualityPreset === preset
                  ? 'border-blue-400 bg-white shadow-sm'
                  : 'border-slate-200 hover:border-slate-300 bg-white'
              }`}
            >
              <div className={`font-medium text-sm ${qualityPreset === preset ? 'text-blue-700' : 'text-slate-800'}`}>
                {preset}
              </div>
              <div className="text-xs text-slate-500 mt-0.5">
                {QUALITY_PRESETS[preset].description}
              </div>
            </button>
          ))}
        </div>
      </FormSection>

      {/* Structure Options */}
      <FormSection
        title="Structure Options"
        description="Control the structural properties of generated proteins"
      >
        <div className="space-y-4">
          {/* Non-loopy toggle */}
          <label className="flex items-center gap-3 p-3 rounded-lg bg-slate-50 hover:bg-slate-100 cursor-pointer transition-colors">
            <input
              type="checkbox"
              checked={isNonLoopy}
              onChange={(e) => setIsNonLoopy(e.target.checked)}
              className="w-4 h-4 rounded border-slate-300 text-blue-600 focus:ring-blue-500"
            />
            <div>
              <div className="font-medium text-sm text-slate-800">Non-loopy Mode</div>
              <div className="text-xs text-slate-500">
                Produces cleaner secondary structures (recommended)
              </div>
            </div>
          </label>

          {/* Symmetry */}
          <FormField label="Symmetry (Optional)" hint="Create homo-oligomeric assemblies">
            <select
              value={symmetry}
              onChange={(e) => setSymmetry(e.target.value)}
              className="w-full px-3 py-2 rounded-lg border border-slate-200 focus:border-blue-400 focus:ring-2 focus:ring-blue-100 outline-none transition-all bg-white text-sm"
            >
              <option value="">No symmetry (monomer)</option>
              {Object.entries(symmetryByCategory).map(([category, options]) => (
                <optgroup key={category} label={category}>
                  {options.map((opt) => (
                    <option key={opt.id} value={opt.id}>
                      {opt.label} ({opt.subunits} subunits)
                    </option>
                  ))}
                </optgroup>
              ))}
            </select>
          </FormField>

          {symmetry && (
            <div className="p-3 rounded-lg bg-slate-50 border border-slate-200">
              <div className="flex items-start gap-2">
                <span className="material-symbols-outlined text-slate-500 text-base">info</span>
                <div className="text-sm text-slate-600">
                  <strong>{symmetry}</strong> symmetry selected. The specified length is per subunit.
                  Total protein size will be{' '}
                  <strong>
                    {SYMMETRY_OPTIONS.find((s) => s.id === symmetry)?.subunits || 1}x
                  </strong>{' '}
                  the specified length.
                </div>
              </div>
            </div>
          )}
        </div>
      </FormSection>

      {/* Advanced Options */}
      <FormSection title="Advanced" description="Additional options for fine-tuning">
        <FormRow>
          <FormField label="Random Seed" hint="For reproducible results">
            <input
              type="number"
              value={seed}
              onChange={(e) => setSeed(e.target.value)}
              placeholder="Optional"
              className="w-full px-3 py-2 rounded-lg border border-slate-200 focus:border-blue-400 focus:ring-2 focus:ring-blue-100 outline-none transition-all text-sm"
            />
          </FormField>
        </FormRow>
      </FormSection>

      {/* Submit Button */}
      <div className="pt-4 border-t border-slate-200">
        <button
          onClick={handleSubmit}
          disabled={!isValid || isSubmitting || !health}
          className={`w-full py-3 px-6 rounded-lg font-semibold text-white transition-all flex items-center justify-center gap-2 ${
            isValid && !isSubmitting && !!health
              ? 'bg-blue-600 hover:bg-blue-700 shadow-lg shadow-blue-600/20'
              : 'bg-slate-300 cursor-not-allowed'
          }`}
        >
          {isSubmitting ? (
            <>
              <span className="material-symbols-outlined animate-spin text-lg">progress_activity</span>
              Submitting...
            </>
          ) : (
            <>
              <span className="material-symbols-outlined text-lg">play_circle</span>
              Generate {numDesigns} Design{numDesigns > 1 ? 's' : ''}
            </>
          )}
        </button>
        {!health && (
          <p className="text-center text-sm text-slate-500 mt-2">
            Connect to backend to enable design
          </p>
        )}
      </div>
    </div>
  );
}
