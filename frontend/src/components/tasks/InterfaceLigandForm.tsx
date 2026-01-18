'use client';

import { useForm } from 'react-hook-form';
import { zodResolver } from '@hookform/resolvers/zod';
import {
  Link2,
  Info,
  Network,
  SlidersHorizontal,
  GitBranch,
  Link,
  ChevronDown,
  CheckCircle,
} from 'lucide-react';

import { cn } from '@/lib/utils';
import {
  ligandFormSchema,
  ligandFormDefaults,
  COMMON_LIGANDS_SMILES,
  type LigandFormValues,
  type LigandApproach,
  type QualityPreset,
} from '@/lib/validations/ligand-form';
import { QUALITY_PRESET_VALUES } from '@/lib/validations/metal-form';
import { RFD3Request, TaskFormProps } from './shared/types';

import { Button } from '@/components/ui/button';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select';
import { Checkbox } from '@/components/ui/checkbox';
import { Badge } from '@/components/ui/badge';
import { Collapsible, CollapsibleContent, CollapsibleTrigger } from '@/components/ui/collapsible';
import {
  Form,
  FormControl,
  FormDescription,
  FormField,
  FormItem,
  FormLabel,
  FormMessage,
} from '@/components/ui/form';

// Design approaches
const APPROACHES: Array<{
  id: LigandApproach;
  name: string;
  description: string;
  details: string;
  icon: typeof Network;
  recommended?: boolean;
}> = [
  {
    id: 'joint',
    name: 'Joint Heterodimer',
    description: 'Both chains simultaneously',
    details: 'Multi-chain RFD3 diffusion. Both chains co-evolve around the ligand.',
    icon: Network,
    recommended: true,
  },
  {
    id: 'asymmetric_rasa',
    name: 'Asymmetric RASA',
    description: 'RASA-conditioned dimer',
    details: 'Complementary burial/exposure for each chain.',
    icon: SlidersHorizontal,
  },
  {
    id: 'induced',
    name: 'Induced Dimerization',
    description: 'Context protein workflow',
    details: 'Design Chain A first, then Chain B using A + ligand as fixed context.',
    icon: GitBranch,
  },
  {
    id: 'asymmetric',
    name: 'Single Chain',
    description: 'One-sided binder only',
    details: 'Design a protein binding one side of ligand.',
    icon: Link,
  },
];

// Expected Results thresholds
const EXPECTED_RESULTS = {
  excellent: {
    affinity: '< -7 kcal/mol',
    contacts: '>= 10 per chain',
    identity: '< 40%',
    antiHomo: '> 80/100',
  },
  good: {
    affinity: '-5 to -7 kcal/mol',
    contacts: '5-10 per chain',
    identity: '< 70%',
    antiHomo: '> 60/100',
  },
  minimum: {
    affinity: '< -2 kcal/mol',
    contacts: '>= 2 per chain',
    identity: '< 90%',
    hbonds: 'At least 1 total',
  },
};

export function InterfaceLigandForm({ onSubmit, isSubmitting, health }: TaskFormProps) {
  const form = useForm<LigandFormValues>({
    resolver: zodResolver(ligandFormSchema),
    defaultValues: ligandFormDefaults,
  });

  const watchApproach = form.watch('approach');
  const watchSelectedPreset = form.watch('selectedPreset');
  const watchCustomSmiles = form.watch('customSmiles');
  const watchQualityPreset = form.watch('qualityPreset');
  const watchUseOriToken = form.watch('useOriToken');

  // Quality preset change
  const handleQualityPresetChange = (preset: QualityPreset) => {
    form.setValue('qualityPreset', preset);
    if (preset !== 'Custom') {
      const values = QUALITY_PRESET_VALUES[preset];
      form.setValue('numTimesteps', values.numTimesteps);
      form.setValue('stepScale', values.stepScale);
      form.setValue('gamma0', values.gamma0);
    }
  };

  // Preset selection
  const handlePresetSelect = (presetId: string) => {
    form.setValue('selectedPreset', presetId);
    const preset = COMMON_LIGANDS_SMILES.find(l => l.id === presetId);
    if (preset) {
      form.setValue('ligandSmiles', preset.smiles);
      form.setValue('customSmiles', false);
    }
  };

  // Form submission
  const handleSubmit = async (data: LigandFormValues) => {
    // Parse ori_offset if provided
    let parsedOriOffset: number[] | undefined;
    if (data.useOriToken && data.oriOffset) {
      parsedOriOffset = data.oriOffset.split(',').map(v => parseFloat(v.trim()));
      if (parsedOriOffset.length !== 3 || parsedOriOffset.some(isNaN)) {
        parsedOriOffset = undefined;
      }
    }

    const request: RFD3Request = {
      task: 'interface_ligand_design',
      approach: data.approach,
      ligand_smiles: data.ligandSmiles,
      chain_length: data.chainLength,
      num_designs: data.numDesigns,
      side: data.side,
      num_timesteps: data.numTimesteps,
      step_scale: data.stepScale,
      gamma_0: data.gamma0,
    };

    if (data.seed) {
      request.seed = parseInt(data.seed, 10);
    }

    if (parsedOriOffset) {
      request.ori_offset = parsedOriOffset;
      request.use_ori_token = true;
    }

    await onSubmit(request);
  };

  const selectedLigandName = COMMON_LIGANDS_SMILES.find(l => l.id === watchSelectedPreset)?.name || 'Custom';

  return (
    <Form {...form}>
      <form onSubmit={form.handleSubmit(handleSubmit)} className="space-y-6">
        {/* Header */}
        <div className="flex items-center gap-3 pb-4 border-b border-border">
          <div className="w-10 h-10 rounded-lg bg-primary flex items-center justify-center">
            <Link2 className="w-5 h-5 text-primary-foreground" />
          </div>
          <div>
            <h2 className="font-semibold text-foreground">Interface Ligand Dimer Design</h2>
            <p className="text-sm text-muted-foreground">Design protein dimers with ligand at the interface (separable topology)</p>
          </div>
        </div>

        {/* Info Banner */}
        <Card className="bg-muted border-border">
          <CardContent className="flex items-start gap-3 pt-4">
            <Info className="w-5 h-5 text-primary mt-0.5 flex-shrink-0" />
            <div className="text-sm text-foreground">
              <p className="font-medium mb-1">Separable Dimer Design</p>
              <p className="text-muted-foreground">
                Unlike buried ligand designs, this creates dimers where chains A and B can physically separate.
                The ligand sits at the interface between the two chains.
              </p>
            </div>
          </CardContent>
        </Card>

        {/* Ligand Selection */}
        <Card>
          <CardHeader>
            <CardTitle className="text-sm flex items-center gap-2">
              Ligand
              <Badge variant="outline" className="text-xs">Required</Badge>
            </CardTitle>
            <CardDescription>Select a ligand or enter custom SMILES</CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            {/* Preset buttons */}
            <div className="grid grid-cols-2 sm:grid-cols-3 gap-2">
              {COMMON_LIGANDS_SMILES.map((lig) => (
                <button
                  key={lig.id}
                  type="button"
                  onClick={() => handlePresetSelect(lig.id)}
                  className={cn(
                    "p-3 rounded-lg border text-left transition-all",
                    watchSelectedPreset === lig.id && !watchCustomSmiles
                      ? "border-primary bg-primary/5 ring-1 ring-primary/20"
                      : "border-border hover:border-muted-foreground/30 hover:bg-muted/50"
                  )}
                >
                  <div className="font-medium text-sm text-foreground">{lig.name}</div>
                  <div className="text-xs text-muted-foreground mt-0.5">{lig.description}</div>
                </button>
              ))}
              <button
                type="button"
                onClick={() => form.setValue('customSmiles', true)}
                className={cn(
                  "p-3 rounded-lg border text-left transition-all",
                  watchCustomSmiles
                    ? "border-primary bg-primary/5 ring-1 ring-primary/20"
                    : "border-border hover:border-muted-foreground/30 hover:bg-muted/50"
                )}
              >
                <div className="font-medium text-sm text-foreground">Custom SMILES</div>
                <div className="text-xs text-muted-foreground mt-0.5">Enter your own molecule</div>
              </button>
            </div>

            {/* SMILES input */}
            <FormField
              control={form.control}
              name="ligandSmiles"
              render={({ field }) => (
                <FormItem>
                  <FormLabel className="text-xs">SMILES String</FormLabel>
                  <FormControl>
                    <Input
                      {...field}
                      placeholder="e.g., c1ccc(cc1)N=Nc2ccccc2"
                      className="font-mono"
                      onChange={(e) => {
                        field.onChange(e);
                        form.setValue('customSmiles', true);
                      }}
                    />
                  </FormControl>
                  <FormMessage />
                </FormItem>
              )}
            />
          </CardContent>
        </Card>

        {/* Design Approach */}
        <Card>
          <CardHeader>
            <CardTitle className="text-sm flex items-center gap-2">
              Design Approach
              <Badge variant="outline" className="text-xs">Required</Badge>
            </CardTitle>
            <CardDescription>Choose how to design the heterodimer</CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-3">
              {APPROACHES.map((app) => {
                const Icon = app.icon;
                return (
                  <button
                    key={app.id}
                    type="button"
                    onClick={() => form.setValue('approach', app.id)}
                    className={cn(
                      "p-4 rounded-lg border text-left transition-all relative",
                      watchApproach === app.id
                        ? "border-primary bg-primary/5 ring-1 ring-primary/20"
                        : "border-border hover:border-muted-foreground/30 hover:bg-muted/50"
                    )}
                  >
                    {app.recommended && (
                      <Badge className="absolute -top-2 -right-2">
                        Recommended
                      </Badge>
                    )}
                    <div className="flex items-center gap-2 mb-2">
                      <Icon className={cn(
                        "w-4 h-4",
                        watchApproach === app.id ? "text-primary" : "text-muted-foreground"
                      )} />
                      <span className="font-medium text-foreground text-sm">{app.name}</span>
                    </div>
                    <p className="text-xs text-muted-foreground">{app.description}</p>
                  </button>
                );
              })}
            </div>

            {/* Side selection for asymmetric */}
            {watchApproach === 'asymmetric' && (
              <div className="p-3 rounded-lg bg-muted/50 border">
                <Label className="text-xs text-muted-foreground">Binding Side</Label>
                <FormField
                  control={form.control}
                  name="side"
                  render={({ field }) => (
                    <FormItem>
                      <div className="mt-2 flex gap-3">
                        <Button
                          type="button"
                          variant={field.value === 'left' ? 'default' : 'outline'}
                          className="flex-1"
                          onClick={() => field.onChange('left')}
                        >
                          Left side (expose right)
                        </Button>
                        <Button
                          type="button"
                          variant={field.value === 'right' ? 'default' : 'outline'}
                          className="flex-1"
                          onClick={() => field.onChange('right')}
                        >
                          Right side (expose left)
                        </Button>
                      </div>
                      <FormDescription className="text-xs mt-2">
                        For azobenzene: left = first phenyl ring, right = second phenyl ring
                      </FormDescription>
                    </FormItem>
                  )}
                />
              </div>
            )}
          </CardContent>
        </Card>

        {/* Design Parameters */}
        <Card>
          <CardHeader>
            <CardTitle className="text-sm flex items-center gap-2">
              Design Parameters
              <Badge variant="outline" className="text-xs">Required</Badge>
            </CardTitle>
          </CardHeader>
          <CardContent className="space-y-4">
            <div className="grid grid-cols-1 sm:grid-cols-3 gap-4">
              <FormField
                control={form.control}
                name="chainLength"
                render={({ field }) => (
                  <FormItem>
                    <FormLabel>Chain Length</FormLabel>
                    <FormControl>
                      <Input placeholder="60-80" {...field} />
                    </FormControl>
                    <FormDescription className="text-xs">Per chain length range</FormDescription>
                    <FormMessage />
                  </FormItem>
                )}
              />

              <FormField
                control={form.control}
                name="numDesigns"
                render={({ field }) => (
                  <FormItem>
                    <FormLabel># Designs</FormLabel>
                    <FormControl>
                      <Input
                        type="number"
                        min={1}
                        max={10}
                        {...field}
                        onChange={(e) => field.onChange(parseInt(e.target.value) || 1)}
                      />
                    </FormControl>
                    <FormMessage />
                  </FormItem>
                )}
              />

              <FormField
                control={form.control}
                name="seed"
                render={({ field }) => (
                  <FormItem>
                    <FormLabel>Seed</FormLabel>
                    <FormControl>
                      <Input placeholder="Random" {...field} />
                    </FormControl>
                    <FormDescription className="text-xs">Optional</FormDescription>
                    <FormMessage />
                  </FormItem>
                )}
              />
            </div>

            {/* Quality Preset */}
            <div className="space-y-3">
              <FormField
                control={form.control}
                name="qualityPreset"
                render={({ field }) => (
                  <FormItem>
                    <FormLabel>Quality Preset</FormLabel>
                    <Select value={field.value} onValueChange={(v) => handleQualityPresetChange(v as QualityPreset)}>
                      <FormControl>
                        <SelectTrigger>
                          <SelectValue />
                        </SelectTrigger>
                      </FormControl>
                      <SelectContent>
                        <SelectItem value="Quick">Quick</SelectItem>
                        <SelectItem value="Balanced">Balanced</SelectItem>
                        <SelectItem value="High Quality">High Quality</SelectItem>
                        <SelectItem value="Binder Optimized">Binder Optimized</SelectItem>
                        <SelectItem value="Custom">Custom</SelectItem>
                      </SelectContent>
                    </Select>
                    <FormMessage />
                  </FormItem>
                )}
              />

              {watchQualityPreset !== 'Custom' && (
                <div className="flex gap-4 text-xs text-muted-foreground">
                  <span>Steps: <span className="font-medium text-foreground">{form.watch('numTimesteps')}</span></span>
                  <span>Scale: <span className="font-medium text-foreground">{form.watch('stepScale')}</span></span>
                  <span>Gamma: <span className="font-medium text-foreground">{form.watch('gamma0')}</span></span>
                </div>
              )}
            </div>
          </CardContent>
        </Card>

        {/* Advanced Options */}
        <Collapsible>
          <CollapsibleTrigger asChild>
            <Button variant="outline" className="w-full justify-between">
              <div className="flex items-center gap-2">
                <SlidersHorizontal className="w-4 h-4" />
                Advanced Options
              </div>
              <ChevronDown className="w-4 h-4" />
            </Button>
          </CollapsibleTrigger>
          <CollapsibleContent className="mt-3 space-y-4 p-4 border rounded-lg">
            <FormField
              control={form.control}
              name="useOriToken"
              render={({ field }) => (
                <FormItem className="flex flex-row items-start space-x-3 space-y-0">
                  <FormControl>
                    <Checkbox checked={field.value} onCheckedChange={field.onChange} />
                  </FormControl>
                  <div className="space-y-1 leading-none">
                    <FormLabel>Use ori_token offset (experimental)</FormLabel>
                  </div>
                </FormItem>
              )}
            />

            {watchUseOriToken && (
              <FormField
                control={form.control}
                name="oriOffset"
                render={({ field }) => (
                  <FormItem>
                    <FormLabel className="text-xs">Offset from ligand (x, y, z)</FormLabel>
                    <FormControl>
                      <Input
                        {...field}
                        placeholder="12.0, 0.0, 0.0"
                        className="font-mono"
                      />
                    </FormControl>
                    <FormDescription className="text-xs">
                      Positions protein origin away from ligand center. Use cautiously.
                    </FormDescription>
                    <FormMessage />
                  </FormItem>
                )}
              />
            )}
          </CollapsibleContent>
        </Collapsible>

        {/* Submit Button */}
        <div className="pt-4 border-t">
          <Button
            type="submit"
            disabled={!form.formState.isValid || isSubmitting || !health}
            className="w-full"
            size="lg"
          >
            {isSubmitting ? (
              <>
                <span className="animate-spin mr-2">
                  <Link2 className="w-4 h-4" />
                </span>
                Designing {watchApproach === 'asymmetric' ? 'Binder' : 'Dimer'}...
              </>
            ) : (
              <>
                <Link2 className="w-4 h-4 mr-2" />
                Design {selectedLigandName} {watchApproach === 'asymmetric' ? 'Binder' : 'Dimer'}
              </>
            )}
          </Button>
          {!health && (
            <p className="text-center text-sm text-muted-foreground mt-2">
              Backend service unavailable
            </p>
          )}

          {/* Success Criteria */}
          <Card className="mt-4 bg-muted/30">
            <CardContent className="pt-4">
              <div className="flex items-center gap-2 mb-3">
                <CheckCircle className="w-4 h-4 text-muted-foreground" />
                <p className="text-xs text-muted-foreground font-semibold uppercase tracking-wide">Success Criteria</p>
              </div>

              <div className="space-y-3">
                {/* Excellent */}
                <div className="p-2 rounded-lg bg-muted border border-border">
                  <div className="flex items-center gap-2 mb-1">
                    <span className="w-2 h-2 rounded-full bg-primary" />
                    <span className="text-xs font-medium text-foreground">Excellent</span>
                  </div>
                  <div className="grid grid-cols-2 gap-x-4 gap-y-0.5 text-xs text-muted-foreground ml-4">
                    <div>GNINA: {EXPECTED_RESULTS.excellent.affinity}</div>
                    <div>Contacts: {EXPECTED_RESULTS.excellent.contacts}</div>
                    <div>Identity: {EXPECTED_RESULTS.excellent.identity}</div>
                    <div>Anti-homo: {EXPECTED_RESULTS.excellent.antiHomo}</div>
                  </div>
                </div>

                {/* Good */}
                <div className="p-2 rounded-lg bg-muted/70 border border-border">
                  <div className="flex items-center gap-2 mb-1">
                    <span className="w-2 h-2 rounded-full bg-muted-foreground" />
                    <span className="text-xs font-medium text-foreground">Good</span>
                  </div>
                  <div className="grid grid-cols-2 gap-x-4 gap-y-0.5 text-xs text-muted-foreground ml-4">
                    <div>GNINA: {EXPECTED_RESULTS.good.affinity}</div>
                    <div>Contacts: {EXPECTED_RESULTS.good.contacts}</div>
                    <div>Identity: {EXPECTED_RESULTS.good.identity}</div>
                    <div>Anti-homo: {EXPECTED_RESULTS.good.antiHomo}</div>
                  </div>
                </div>

                {/* Minimum */}
                <div className="p-2 rounded-lg bg-muted/50 border border-border">
                  <div className="flex items-center gap-2 mb-1">
                    <span className="w-2 h-2 rounded-full bg-muted-foreground/60" />
                    <span className="text-xs font-medium text-foreground">Minimum Viable</span>
                  </div>
                  <div className="grid grid-cols-2 gap-x-4 gap-y-0.5 text-xs text-muted-foreground ml-4">
                    <div>GNINA: {EXPECTED_RESULTS.minimum.affinity}</div>
                    <div>Contacts: {EXPECTED_RESULTS.minimum.contacts}</div>
                    <div>Identity: {EXPECTED_RESULTS.minimum.identity}</div>
                    <div>H-bonds: {EXPECTED_RESULTS.minimum.hbonds}</div>
                  </div>
                </div>
              </div>

              <p className="text-xs text-muted-foreground mt-3 text-center">
                Pipeline: Backbone → LigandMPNN → Validation → PLIP Analysis
              </p>
            </CardContent>
          </Card>
        </div>
      </form>
    </Form>
  );
}
