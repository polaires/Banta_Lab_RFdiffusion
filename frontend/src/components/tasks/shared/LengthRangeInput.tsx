'use client';

import { useState, useEffect } from 'react';
import { AlertCircle } from 'lucide-react';

interface LengthRangeInputProps {
  value: string;
  onChange: (value: string) => void;
  label?: string;
  min?: number;
  max?: number;
  placeholder?: string;
  hint?: string;
  className?: string;
}

// Validates length range format: "60-80" or "100"
function validateRange(value: string, min?: number, max?: number): string | null {
  if (!value.trim()) return null;

  // Check format
  const rangeMatch = value.match(/^(\d+)(?:-(\d+))?$/);
  if (!rangeMatch) {
    return 'Format: number or min-max (e.g., "60-80")';
  }

  const start = parseInt(rangeMatch[1], 10);
  const end = rangeMatch[2] ? parseInt(rangeMatch[2], 10) : start;

  // Validate range order
  if (end < start) {
    return 'End must be greater than or equal to start';
  }

  // Validate bounds
  if (min !== undefined && start < min) {
    return `Minimum length is ${min}`;
  }
  if (max !== undefined && end > max) {
    return `Maximum length is ${max}`;
  }

  return null;
}

export function LengthRangeInput({
  value,
  onChange,
  label,
  min = 20,
  max = 500,
  placeholder = '60-80',
  hint,
  className = '',
}: LengthRangeInputProps) {
  const [error, setError] = useState<string | null>(null);
  const [isFocused, setIsFocused] = useState(false);

  useEffect(() => {
    if (!isFocused && value) {
      setError(validateRange(value, min, max));
    }
  }, [value, min, max, isFocused]);

  const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const newValue = e.target.value;
    onChange(newValue);
    // Clear error while typing
    if (error) setError(null);
  };

  const handleBlur = () => {
    setIsFocused(false);
    if (value) {
      setError(validateRange(value, min, max));
    }
  };

  return (
    <div className={`space-y-1.5 ${className}`}>
      {label && (
        <label className="block text-sm font-medium text-foreground">
          {label}
        </label>
      )}
      <div className="relative">
        <input
          type="text"
          value={value}
          onChange={handleChange}
          onFocus={() => setIsFocused(true)}
          onBlur={handleBlur}
          placeholder={placeholder}
          className={`w-full px-3 py-2 border rounded-lg text-sm focus:ring-2 focus:ring-ring focus:border-transparent bg-card ${
            error ? 'border-red-300 bg-red-50' : 'border-border'
          }`}
        />
        {value && !error && (
          <span className="absolute right-3 top-1/2 -translate-y-1/2 text-xs text-muted-foreground">
            residues
          </span>
        )}
      </div>
      {error && (
        <p className="text-xs text-red-500 flex items-center gap-1">
          <AlertCircle className="h-4 w-4" />
          {error}
        </p>
      )}
      {!error && hint && (
        <p className="text-xs text-muted-foreground">{hint}</p>
      )}
    </div>
  );
}

export default LengthRangeInput;
