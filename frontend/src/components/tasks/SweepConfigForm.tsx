'use client';

import { useState } from 'react';
import { cn } from '@/lib/utils';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Checkbox } from '@/components/ui/checkbox';
import { Badge } from '@/components/ui/badge';
import { Slider } from '@/components/ui/slider';
import { Settings2, Layers, Zap } from 'lucide-react';
import type { SweepConfig } from '@/lib/store';

interface SweepConfigFormProps {
  configs: SweepConfig[];
  onConfigsChange: (configs: SweepConfig[]) => void;
  designsPerConfig: number;
  onDesignsPerConfigChange: (count: number) => void;
  disabled?: boolean;
}

// Default contig ranges for each size
const SIZE_RANGES: Record<'small' | 'medium' | 'large', string> = {
  small: '50-70',
  medium: '70-90',
  large: '90-120',
};

// CFG scale values
const CFG_VALUES = [1.5, 2.0, 2.5];

// Generate default sweep configs
export function generateDefaultSweepConfigs(designsPerConfig: number = 10): SweepConfig[] {
  const configs: SweepConfig[] = [];

  for (const size of ['small', 'medium', 'large'] as const) {
    for (const cfg of CFG_VALUES) {
      const cfgName = cfg === 1.5 ? 'low_cfg' : cfg === 2.0 ? 'mid_cfg' : 'high_cfg';
      configs.push({
        name: `${size}_${cfgName}`,
        contigSize: size,
        contigRange: SIZE_RANGES[size],
        cfgScale: cfg,
        numDesigns: designsPerConfig,
      });
    }
  }

  return configs;
}

export function SweepConfigForm({
  configs,
  onConfigsChange,
  designsPerConfig,
  onDesignsPerConfigChange,
  disabled = false,
}: SweepConfigFormProps) {
  const [customRanges, setCustomRanges] = useState<Record<string, string>>({
    small: SIZE_RANGES.small,
    medium: SIZE_RANGES.medium,
    large: SIZE_RANGES.large,
  });

  const [selectedSizes, setSelectedSizes] = useState<Set<'small' | 'medium' | 'large'>>(
    new Set(['small', 'medium', 'large'])
  );

  const [selectedCfgValues, setSelectedCfgValues] = useState<Set<number>>(
    new Set(CFG_VALUES)
  );

  // Toggle size selection
  const toggleSize = (size: 'small' | 'medium' | 'large') => {
    const newSizes = new Set(selectedSizes);
    if (newSizes.has(size)) {
      newSizes.delete(size);
    } else {
      newSizes.add(size);
    }
    setSelectedSizes(newSizes);
    updateConfigs(newSizes, selectedCfgValues, customRanges, designsPerConfig);
  };

  // Toggle CFG value selection
  const toggleCfg = (cfg: number) => {
    const newCfgValues = new Set(selectedCfgValues);
    if (newCfgValues.has(cfg)) {
      newCfgValues.delete(cfg);
    } else {
      newCfgValues.add(cfg);
    }
    setSelectedCfgValues(newCfgValues);
    updateConfigs(selectedSizes, newCfgValues, customRanges, designsPerConfig);
  };

  // Update contig range for a size
  const updateRange = (size: 'small' | 'medium' | 'large', range: string) => {
    const newRanges = { ...customRanges, [size]: range };
    setCustomRanges(newRanges);
    updateConfigs(selectedSizes, selectedCfgValues, newRanges, designsPerConfig);
  };

  // Regenerate configs based on selections
  const updateConfigs = (
    sizes: Set<'small' | 'medium' | 'large'>,
    cfgValues: Set<number>,
    ranges: Record<string, string>,
    numDesigns: number
  ) => {
    const newConfigs: SweepConfig[] = [];

    for (const size of ['small', 'medium', 'large'] as const) {
      if (!sizes.has(size)) continue;

      for (const cfg of CFG_VALUES) {
        if (!cfgValues.has(cfg)) continue;

        const cfgName = cfg === 1.5 ? 'low_cfg' : cfg === 2.0 ? 'mid_cfg' : 'high_cfg';
        newConfigs.push({
          name: `${size}_${cfgName}`,
          contigSize: size,
          contigRange: ranges[size] || SIZE_RANGES[size],
          cfgScale: cfg,
          numDesigns: numDesigns,
        });
      }
    }

    onConfigsChange(newConfigs);
  };

  // Update designs per config
  const handleDesignsChange = (count: number) => {
    onDesignsPerConfigChange(count);
    const updatedConfigs = configs.map((c) => ({ ...c, numDesigns: count }));
    onConfigsChange(updatedConfigs);
  };

  const totalConfigs = selectedSizes.size * selectedCfgValues.size;
  const totalDesigns = totalConfigs * designsPerConfig;

  return (
    <Card>
      <CardHeader>
        <CardTitle className="text-sm flex items-center gap-2">
          <Settings2 className="w-4 h-4" />
          Sweep Configuration
          <Badge variant="secondary" className="ml-auto">
            {totalConfigs} configs
          </Badge>
        </CardTitle>
        <CardDescription>
          Configure parameter sweep across contig sizes and CFG values
        </CardDescription>
      </CardHeader>
      <CardContent className="space-y-6">
        {/* Contig Sizes */}
        <div className="space-y-3">
          <Label className="flex items-center gap-2">
            <Layers className="w-4 h-4 text-muted-foreground" />
            Contig Sizes
          </Label>
          <div className="grid grid-cols-3 gap-3">
            {(['small', 'medium', 'large'] as const).map((size) => (
              <div key={size} className="space-y-2">
                <div className="flex items-center space-x-2">
                  <Checkbox
                    id={`size-${size}`}
                    checked={selectedSizes.has(size)}
                    onCheckedChange={() => toggleSize(size)}
                    disabled={disabled}
                  />
                  <Label
                    htmlFor={`size-${size}`}
                    className="text-sm font-medium capitalize"
                  >
                    {size}
                  </Label>
                </div>
                <Input
                  value={customRanges[size]}
                  onChange={(e) => updateRange(size, e.target.value)}
                  placeholder={SIZE_RANGES[size]}
                  disabled={disabled || !selectedSizes.has(size)}
                  className="h-8 text-sm font-mono"
                />
              </div>
            ))}
          </div>
        </div>

        {/* CFG Scale Values */}
        <div className="space-y-3">
          <Label className="flex items-center gap-2">
            <Zap className="w-4 h-4 text-muted-foreground" />
            CFG Scale Values
          </Label>
          <div className="flex gap-3">
            {CFG_VALUES.map((cfg) => (
              <button
                key={cfg}
                type="button"
                onClick={() => toggleCfg(cfg)}
                disabled={disabled}
                className={cn(
                  "flex-1 px-4 py-2 rounded-md border text-sm font-medium transition-all",
                  selectedCfgValues.has(cfg)
                    ? "border-primary bg-primary/10 text-primary"
                    : "border-border text-muted-foreground hover:border-muted-foreground/50",
                  disabled && "opacity-50 cursor-not-allowed"
                )}
              >
                {cfg.toFixed(1)}
                <div className="text-xs text-muted-foreground mt-1">
                  {cfg === 1.5 ? 'Low' : cfg === 2.0 ? 'Medium' : 'High'}
                </div>
              </button>
            ))}
          </div>
        </div>

        {/* Designs Per Config */}
        <div className="space-y-3">
          <div className="flex items-center justify-between">
            <Label>Designs per Configuration</Label>
            <span className="text-sm font-medium">{designsPerConfig}</span>
          </div>
          <Slider
            value={[designsPerConfig]}
            onValueChange={([value]) => handleDesignsChange(value)}
            min={1}
            max={50}
            step={1}
            disabled={disabled}
          />
          <div className="flex justify-between text-xs text-muted-foreground">
            <span>1</span>
            <span>50</span>
          </div>
        </div>

        {/* Summary */}
        <div className="p-3 rounded-lg bg-muted/50 border text-sm">
          <div className="flex justify-between items-center">
            <span className="text-muted-foreground">Total Configurations:</span>
            <span className="font-medium">{totalConfigs}</span>
          </div>
          <div className="flex justify-between items-center mt-1">
            <span className="text-muted-foreground">Total Designs:</span>
            <span className="font-medium">{totalDesigns}</span>
          </div>
        </div>
      </CardContent>
    </Card>
  );
}
