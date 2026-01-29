'use client';

import { useState, useCallback } from 'react';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Slider } from '@/components/ui/slider';
import { Switch } from '@/components/ui/switch';
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from '@/components/ui/select';
import type { StepParameterField } from '@/lib/pipeline-types';

interface StepParameterEditorProps {
  schema: StepParameterField[];
  params: Record<string, unknown>;
  onChange: (params: Record<string, unknown>) => void;
}

// FIX #14: Validate a field value against its schema
function validateField(field: StepParameterField, value: unknown): string | null {
  if (field.required && (value === undefined || value === null || value === '')) {
    return `${field.label} is required`;
  }

  if (field.type === 'number' || field.type === 'slider') {
    const num = Number(value);
    if (value !== undefined && value !== '' && isNaN(num)) {
      return `${field.label} must be a number`;
    }
    if (field.range) {
      if (num < field.range.min) return `Min: ${field.range.min}`;
      if (num > field.range.max) return `Max: ${field.range.max}`;
    }
  }

  return null;
}

export function StepParameterEditor({ schema, params, onChange }: StepParameterEditorProps) {
  const visibleFields = schema.filter(f => f.type !== 'hidden');
  // FIX #14: Track which fields have been touched (show errors only after interaction)
  const [touched, setTouched] = useState<Set<string>>(new Set());

  if (visibleFields.length === 0) return null;

  const updateParam = (id: string, value: unknown) => {
    setTouched(prev => new Set(prev).add(id));
    onChange({ ...params, [id]: value });
  };

  const markTouched = (id: string) => {
    setTouched(prev => new Set(prev).add(id));
  };

  return (
    <div className="space-y-3">
      {visibleFields.map(field => {
        const value = params[field.id] ?? field.defaultValue;
        const error = touched.has(field.id) ? validateField(field, value) : null;

        return (
          <div key={field.id} className="space-y-1.5">
            <div className="flex items-center justify-between">
              <Label htmlFor={field.id} className="text-xs font-medium">
                {field.label}
                {field.required && <span className="text-destructive ml-0.5">*</span>}
              </Label>
              {field.type === 'slider' && field.range && (
                <span className="text-xs text-muted-foreground font-mono">
                  {String(value ?? field.range.min)}
                </span>
              )}
            </div>

            {field.type === 'number' && (
              <Input
                id={field.id}
                type="number"
                value={String(value ?? '')}
                onChange={e => {
                  const raw = e.target.value;
                  // Allow empty for clearing, otherwise parse
                  updateParam(field.id, raw === '' ? '' : Number(raw));
                }}
                onBlur={() => markTouched(field.id)}
                min={field.range?.min}
                max={field.range?.max}
                step={field.range?.step}
                className={`h-8 text-sm ${error ? 'border-destructive' : ''}`}
              />
            )}

            {field.type === 'text' && (
              <Input
                id={field.id}
                type="text"
                value={String(value ?? '')}
                onChange={e => updateParam(field.id, e.target.value)}
                onBlur={() => markTouched(field.id)}
                className={`h-8 text-sm ${error ? 'border-destructive' : ''}`}
              />
            )}

            {field.type === 'select' && field.options && (
              <Select
                value={String(value ?? '')}
                onValueChange={v => updateParam(field.id, v)}
              >
                <SelectTrigger className={`h-8 text-sm ${error ? 'border-destructive' : ''}`}>
                  <SelectValue />
                </SelectTrigger>
                <SelectContent>
                  {field.options.map(opt => (
                    <SelectItem key={opt.value} value={opt.value}>
                      {opt.label}
                    </SelectItem>
                  ))}
                </SelectContent>
              </Select>
            )}

            {field.type === 'slider' && field.range && (
              <Slider
                id={field.id}
                value={[Number(value ?? field.range.min)]}
                min={field.range.min}
                max={field.range.max}
                step={field.range.step}
                onValueChange={([v]) => updateParam(field.id, v)}
                className="py-1"
              />
            )}

            {field.type === 'boolean' && (
              <Switch
                id={field.id}
                checked={Boolean(value ?? false)}
                onCheckedChange={v => updateParam(field.id, v)}
              />
            )}

            {/* Validation error */}
            {error && (
              <p className="text-[10px] text-destructive">{error}</p>
            )}

            {/* Help text — only show when no error */}
            {!error && field.helpText && (
              <p className="text-[10px] text-muted-foreground">{field.helpText}</p>
            )}
          </div>
        );
      })}
    </div>
  );
}
