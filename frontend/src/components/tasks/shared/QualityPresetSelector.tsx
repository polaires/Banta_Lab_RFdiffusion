'use client';

import { QUALITY_PRESETS } from './types';

export type QualityPreset = keyof typeof QUALITY_PRESETS | 'Custom';

export interface QualityParams {
  num_timesteps: number;
  step_scale: number;
  gamma_0: number;
}

interface QualityPresetSelectorProps {
  value: QualityPreset;
  onChange: (preset: QualityPreset, params: QualityParams) => void;
  showDescription?: boolean;
  className?: string;
}

export function QualityPresetSelector({
  value,
  onChange,
  showDescription = true,
  className = '',
}: QualityPresetSelectorProps) {
  const presetKeys = Object.keys(QUALITY_PRESETS) as (keyof typeof QUALITY_PRESETS)[];
  const currentPreset = value !== 'Custom' ? QUALITY_PRESETS[value] : null;

  const handleChange = (e: React.ChangeEvent<HTMLSelectElement>) => {
    const newValue = e.target.value as QualityPreset;
    if (newValue === 'Custom') {
      // Keep current params when switching to custom
      onChange(newValue, currentPreset || QUALITY_PRESETS.Balanced);
    } else {
      const preset = QUALITY_PRESETS[newValue];
      onChange(newValue, {
        num_timesteps: preset.num_timesteps,
        step_scale: preset.step_scale,
        gamma_0: preset.gamma_0,
      });
    }
  };

  return (
    <div className={`space-y-2 ${className}`}>
      <select
        value={value}
        onChange={handleChange}
        className="w-full px-3 py-2 border border-slate-200 rounded-lg text-sm focus:ring-2 focus:ring-blue-500 focus:border-transparent bg-white"
      >
        {presetKeys.map((key) => (
          <option key={key} value={key}>
            {key}
          </option>
        ))}
        <option value="Custom">Custom</option>
      </select>

      {showDescription && currentPreset && (
        <p className="text-xs text-slate-500 flex items-center gap-1.5">
          <span className="material-symbols-outlined text-sm text-slate-400">info</span>
          {currentPreset.description}
        </p>
      )}

      {showDescription && value !== 'Custom' && currentPreset && (
        <div className="flex gap-3 text-xs text-slate-400">
          <span>Steps: {currentPreset.num_timesteps}</span>
          <span>Scale: {currentPreset.step_scale}</span>
          <span>Gamma: {currentPreset.gamma_0}</span>
        </div>
      )}
    </div>
  );
}

export default QualityPresetSelector;
