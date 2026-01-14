'use client';

import { useState, useCallback } from 'react';

interface ResidueSelectorProps {
  value: string[];
  onChange: (residues: string[]) => void;
  label: string;
  description?: string;
  placeholder?: string;
  chainFilter?: string;
  className?: string;
}

// Validates residue format: A15, B20, etc.
function isValidResidue(residue: string): boolean {
  return /^[A-Z]\d+$/.test(residue);
}

// Parse residue string to get chain and number
function parseResidue(residue: string): { chain: string; num: number } | null {
  const match = residue.match(/^([A-Z])(\d+)$/);
  if (!match) return null;
  return { chain: match[1], num: parseInt(match[2], 10) };
}

export function ResidueSelector({
  value,
  onChange,
  label,
  description,
  placeholder = 'e.g., A15, A20, B30',
  chainFilter,
  className = '',
}: ResidueSelectorProps) {
  const [inputValue, setInputValue] = useState('');
  const [error, setError] = useState<string | null>(null);

  const addResidue = useCallback(() => {
    const trimmed = inputValue.trim().toUpperCase();
    if (!trimmed) return;

    // Split by comma or space
    const residues = trimmed.split(/[,\s]+/).filter(Boolean);
    const newResidues: string[] = [];
    const errors: string[] = [];

    for (const res of residues) {
      if (!isValidResidue(res)) {
        errors.push(`"${res}" - invalid format`);
        continue;
      }

      // Check chain filter if provided
      if (chainFilter) {
        const parsed = parseResidue(res);
        if (parsed && parsed.chain !== chainFilter) {
          errors.push(`"${res}" - must be chain ${chainFilter}`);
          continue;
        }
      }

      // Check for duplicates
      if (value.includes(res) || newResidues.includes(res)) {
        continue; // Skip duplicates silently
      }

      newResidues.push(res);
    }

    if (newResidues.length > 0) {
      onChange([...value, ...newResidues]);
      setInputValue('');
      setError(null);
    }

    if (errors.length > 0) {
      setError(errors.join(', '));
    }
  }, [inputValue, value, onChange, chainFilter]);

  const removeResidue = useCallback((residue: string) => {
    onChange(value.filter((r) => r !== residue));
  }, [value, onChange]);

  const handleKeyDown = (e: React.KeyboardEvent) => {
    if (e.key === 'Enter' || e.key === ',') {
      e.preventDefault();
      addResidue();
    }
  };

  // Group residues by chain for display
  const groupedByChain = value.reduce((acc, res) => {
    const parsed = parseResidue(res);
    if (parsed) {
      if (!acc[parsed.chain]) acc[parsed.chain] = [];
      acc[parsed.chain].push(res);
    }
    return acc;
  }, {} as Record<string, string[]>);

  return (
    <div className={`space-y-2 ${className}`}>
      <div>
        <label className="block text-sm font-medium text-slate-700">{label}</label>
        {description && (
          <p className="text-xs text-slate-500 mt-0.5">{description}</p>
        )}
      </div>

      {/* Input field */}
      <div className="flex gap-2">
        <input
          type="text"
          value={inputValue}
          onChange={(e) => {
            setInputValue(e.target.value);
            setError(null);
          }}
          onKeyDown={handleKeyDown}
          placeholder={placeholder}
          className={`flex-1 px-3 py-2 border rounded-lg text-sm focus:ring-2 focus:ring-blue-500 focus:border-transparent ${
            error ? 'border-red-300' : 'border-slate-200'
          }`}
        />
        <button
          type="button"
          onClick={addResidue}
          disabled={!inputValue.trim()}
          className="px-3 py-2 bg-blue-600 text-white text-sm font-medium rounded-lg hover:bg-blue-700 disabled:bg-slate-300 disabled:cursor-not-allowed transition-colors"
        >
          Add
        </button>
      </div>

      {error && (
        <p className="text-xs text-red-500">{error}</p>
      )}

      {/* Selected residues as chips */}
      {value.length > 0 && (
        <div className="space-y-2">
          {Object.entries(groupedByChain).map(([chain, residues]) => (
            <div key={chain} className="flex flex-wrap gap-1.5 items-center">
              <span className="text-xs font-medium text-slate-500 mr-1">
                Chain {chain}:
              </span>
              {residues.map((res) => (
                <span
                  key={res}
                  className="inline-flex items-center gap-1 px-2 py-1 bg-blue-50 text-blue-700 text-xs font-medium rounded-full border border-blue-100"
                >
                  {res}
                  <button
                    type="button"
                    onClick={() => removeResidue(res)}
                    className="hover:text-blue-900 transition-colors"
                  >
                    <span className="material-symbols-outlined text-sm">close</span>
                  </button>
                </span>
              ))}
            </div>
          ))}
        </div>
      )}

      {value.length === 0 && (
        <p className="text-xs text-slate-400 italic">No residues selected</p>
      )}
    </div>
  );
}

export default ResidueSelector;
