'use client';

import { cn } from '@/lib/utils';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Label } from '@/components/ui/label';
import { Slider } from '@/components/ui/slider';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { Filter, AlertTriangle } from 'lucide-react';
import type { PipelineFilters } from '@/lib/store';

interface FilterThresholdEditorProps {
  filters: PipelineFilters;
  onFiltersChange: (filters: Partial<PipelineFilters>) => void;
  disabled?: boolean;
}

// Preset definitions
const FILTER_PRESETS = {
  strict: {
    label: 'Strict',
    description: 'High confidence only',
    values: { plddt: 0.80, ptm: 0.80, pae: 5.0 },
  },
  relaxed: {
    label: 'Relaxed',
    description: 'Include borderline designs',
    values: { plddt: 0.70, ptm: 0.65, pae: 10.0 },
  },
  custom: {
    label: 'Custom',
    description: 'Your own thresholds',
    values: null,
  },
} as const;

type PresetKey = keyof typeof FILTER_PRESETS;

function getActivePreset(filters: PipelineFilters): PresetKey {
  if (
    filters.plddt === FILTER_PRESETS.strict.values.plddt &&
    filters.ptm === FILTER_PRESETS.strict.values.ptm &&
    filters.pae === FILTER_PRESETS.strict.values.pae
  ) {
    return 'strict';
  }
  if (
    filters.plddt === FILTER_PRESETS.relaxed.values.plddt &&
    filters.ptm === FILTER_PRESETS.relaxed.values.ptm &&
    filters.pae === FILTER_PRESETS.relaxed.values.pae
  ) {
    return 'relaxed';
  }
  return 'custom';
}

export function FilterThresholdEditor({
  filters,
  onFiltersChange,
  disabled = false,
}: FilterThresholdEditorProps) {
  const activePreset = getActivePreset(filters);

  const applyPreset = (preset: PresetKey) => {
    if (preset === 'custom') return;
    const values = FILTER_PRESETS[preset].values;
    if (values) {
      onFiltersChange(values);
    }
  };

  return (
    <Card>
      <CardHeader>
        <CardTitle className="text-sm flex items-center gap-2">
          <Filter className="w-4 h-4" />
          Filter Thresholds
          <Badge variant="outline" className="ml-auto capitalize">
            {activePreset}
          </Badge>
        </CardTitle>
        <CardDescription>
          Designs below these thresholds will be filtered out
        </CardDescription>
      </CardHeader>
      <CardContent className="space-y-6">
        {/* Preset Buttons */}
        <div className="flex gap-2">
          {(Object.keys(FILTER_PRESETS) as PresetKey[]).map((key) => {
            const preset = FILTER_PRESETS[key];
            const isActive = activePreset === key;

            return (
              <Button
                key={key}
                type="button"
                variant={isActive ? 'default' : 'outline'}
                size="sm"
                onClick={() => applyPreset(key)}
                disabled={disabled || key === 'custom'}
                className="flex-1"
              >
                {preset.label}
              </Button>
            );
          })}
        </div>

        {/* pLDDT Slider */}
        <div className="space-y-3">
          <div className="flex items-center justify-between">
            <Label className="flex items-center gap-2">
              pLDDT
              <span className="text-xs text-muted-foreground font-normal">
                (predicted LDDT)
              </span>
            </Label>
            <span className="text-sm font-medium">
              {filters.plddt >= 1 ? filters.plddt.toFixed(0) : (filters.plddt * 100).toFixed(0)}%
            </span>
          </div>
          <Slider
            value={[filters.plddt]}
            onValueChange={([value]) => onFiltersChange({ plddt: value })}
            min={0.50}
            max={0.95}
            step={0.01}
            disabled={disabled}
          />
          <div className="flex justify-between text-xs text-muted-foreground">
            <span>50%</span>
            <span className="flex items-center gap-1">
              <span className={cn(
                "inline-block w-2 h-2 rounded-full",
                filters.plddt >= 0.80 ? "bg-success" : filters.plddt >= 0.70 ? "bg-warning" : "bg-destructive"
              )} />
              {filters.plddt >= 0.80 ? 'Confident' : filters.plddt >= 0.70 ? 'Marginal' : 'Low'}
            </span>
            <span>95%</span>
          </div>
        </div>

        {/* pTM Slider */}
        <div className="space-y-3">
          <div className="flex items-center justify-between">
            <Label className="flex items-center gap-2">
              pTM
              <span className="text-xs text-muted-foreground font-normal">
                (predicted TM-score)
              </span>
            </Label>
            <span className="text-sm font-medium">
              {filters.ptm >= 1 ? filters.ptm.toFixed(0) : (filters.ptm * 100).toFixed(0)}%
            </span>
          </div>
          <Slider
            value={[filters.ptm]}
            onValueChange={([value]) => onFiltersChange({ ptm: value })}
            min={0.50}
            max={0.95}
            step={0.01}
            disabled={disabled}
          />
          <div className="flex justify-between text-xs text-muted-foreground">
            <span>50%</span>
            <span className="flex items-center gap-1">
              <span className={cn(
                "inline-block w-2 h-2 rounded-full",
                filters.ptm >= 0.80 ? "bg-success" : filters.ptm >= 0.65 ? "bg-warning" : "bg-destructive"
              )} />
              {filters.ptm >= 0.80 ? 'Good match' : filters.ptm >= 0.65 ? 'Borderline' : 'Poor'}
            </span>
            <span>95%</span>
          </div>
        </div>

        {/* PAE Slider */}
        <div className="space-y-3">
          <div className="flex items-center justify-between">
            <Label className="flex items-center gap-2">
              PAE
              <span className="text-xs text-muted-foreground font-normal">
                (predicted aligned error)
              </span>
            </Label>
            <span className="text-sm font-medium">
              {filters.pae.toFixed(1)} A
            </span>
          </div>
          <Slider
            value={[filters.pae]}
            onValueChange={([value]) => onFiltersChange({ pae: value })}
            min={3.0}
            max={15.0}
            step={0.5}
            disabled={disabled}
          />
          <div className="flex justify-between text-xs text-muted-foreground">
            <span>3 A (strict)</span>
            <span className="flex items-center gap-1">
              <span className={cn(
                "inline-block w-2 h-2 rounded-full",
                filters.pae <= 5.0 ? "bg-success" : filters.pae <= 10.0 ? "bg-warning" : "bg-destructive"
              )} />
              {filters.pae <= 5.0 ? 'Low error' : filters.pae <= 10.0 ? 'Moderate' : 'High'}
            </span>
            <span>15 A (relaxed)</span>
          </div>
        </div>

        {/* Warning for very relaxed settings */}
        {(filters.plddt < 0.70 || filters.ptm < 0.65 || filters.pae > 10.0) && (
          <div className="flex items-start gap-2 p-3 rounded-lg bg-warning/10 border border-warning/20 text-warning">
            <AlertTriangle className="w-4 h-4 mt-0.5 flex-shrink-0" />
            <p className="text-xs">
              Very relaxed thresholds may include low-quality designs that are unlikely to fold correctly.
            </p>
          </div>
        )}
      </CardContent>
    </Card>
  );
}
