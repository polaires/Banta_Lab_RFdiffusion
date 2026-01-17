'use client';

import { useState, useMemo } from 'react';
import { FlaskConical, AlertTriangle } from 'lucide-react';
import { COMMON_LIGANDS } from './types';

interface LigandSelectorProps {
  value: string;
  onChange: (ligandId: string) => void;
  showCategories?: boolean;
  allowCustom?: boolean;
  customPlaceholder?: string;
  label?: string;
  className?: string;
}

export function LigandSelector({
  value,
  onChange,
  showCategories = true,
  allowCustom = true,
  customPlaceholder = 'Enter PDB ligand code',
  label,
  className = '',
}: LigandSelectorProps) {
  const [isCustom, setIsCustom] = useState(false);
  const [customValue, setCustomValue] = useState('');

  // Group ligands by category
  const groupedLigands = useMemo(() => {
    const groups: Record<string, typeof COMMON_LIGANDS[number][]> = {};
    COMMON_LIGANDS.forEach((ligand) => {
      if (!groups[ligand.category]) {
        groups[ligand.category] = [];
      }
      groups[ligand.category].push(ligand);
    });
    return groups;
  }, []);

  // Check if current value is a preset
  const isPreset = COMMON_LIGANDS.some((l) => l.id === value);
  const selectedLigand = COMMON_LIGANDS.find((l) => l.id === value);

  const handlePresetChange = (e: React.ChangeEvent<HTMLSelectElement>) => {
    const newValue = e.target.value;
    if (newValue === '__custom__') {
      setIsCustom(true);
    } else {
      setIsCustom(false);
      onChange(newValue);
    }
  };

  const handleCustomChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const newValue = e.target.value.toUpperCase();
    setCustomValue(newValue);
    onChange(newValue);
  };

  return (
    <div className={`space-y-2 ${className}`}>
      {label && (
        <label className="block text-sm font-medium text-slate-700">{label}</label>
      )}

      {!isCustom ? (
        <select
          value={isPreset ? value : '__custom__'}
          onChange={handlePresetChange}
          className="w-full px-3 py-2 border border-slate-200 rounded-lg text-sm focus:ring-2 focus:ring-blue-500 focus:border-transparent bg-white"
        >
          <option value="" disabled>
            Select a ligand...
          </option>

          {showCategories ? (
            // Grouped by category
            Object.entries(groupedLigands).map(([category, ligands]) => (
              <optgroup key={category} label={category}>
                {ligands.map((ligand) => (
                  <option key={ligand.id} value={ligand.id}>
                    {ligand.label} ({ligand.id})
                  </option>
                ))}
              </optgroup>
            ))
          ) : (
            // Flat list
            COMMON_LIGANDS.map((ligand) => (
              <option key={ligand.id} value={ligand.id}>
                {ligand.label} ({ligand.id})
              </option>
            ))
          )}

          {allowCustom && (
            <option value="__custom__">Enter custom code...</option>
          )}
        </select>
      ) : (
        <div className="flex gap-2">
          <input
            type="text"
            value={customValue}
            onChange={handleCustomChange}
            placeholder={customPlaceholder}
            maxLength={10}
            className="flex-1 px-3 py-2 border border-slate-200 rounded-lg text-sm focus:ring-2 focus:ring-blue-500 focus:border-transparent font-mono uppercase"
          />
          <button
            type="button"
            onClick={() => {
              setIsCustom(false);
              setCustomValue('');
              onChange('');
            }}
            className="px-3 py-2 text-slate-600 hover:text-slate-800 border border-slate-200 rounded-lg text-sm transition-colors"
          >
            Cancel
          </button>
        </div>
      )}

      {/* Show selected ligand info */}
      {selectedLigand && (
        <div className="flex items-center gap-2 text-xs text-slate-500">
          <FlaskConical className="h-4 w-4 text-slate-400" />
          <span>{selectedLigand.category}</span>
          <span className="text-slate-300">|</span>
          <a
            href={`https://www.rcsb.org/ligand/${selectedLigand.id}`}
            target="_blank"
            rel="noopener noreferrer"
            className="text-blue-500 hover:underline"
          >
            View in PDB
          </a>
        </div>
      )}

      {isCustom && customValue && (
        <div className="flex items-center gap-2 text-xs text-slate-500">
          <AlertTriangle className="h-4 w-4 text-amber-500" />
          <span>Custom code - verify it exists in the PDB</span>
          <a
            href={`https://www.rcsb.org/ligand/${customValue}`}
            target="_blank"
            rel="noopener noreferrer"
            className="text-blue-500 hover:underline"
          >
            Check
          </a>
        </div>
      )}
    </div>
  );
}

export default LigandSelector;
