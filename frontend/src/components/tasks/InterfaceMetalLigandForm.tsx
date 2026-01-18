'use client';

import { useForm } from 'react-hook-form';
import { zodResolver } from '@hookform/resolvers/zod';
import {
  Sparkles,
  Info,
  Network,
  ChevronDown,
  CheckCircle,
  Upload,
  Atom,
  FlaskConical,
} from 'lucide-react';

import { cn } from '@/lib/utils';
import {
  metalLigandFormSchema,
  metalLigandFormDefaults,
  METAL_LIGAND_TEMPLATES,
  METAL_LIGAND_APPROACHES,
  getTemplateById,
  getHsabBiasForMetal,
  type MetalLigandFormValues,
  type MetalLigandApproach,
} from '@/lib/validations/metal-ligand-form';
import { QUALITY_PRESET_VALUES, type QualityPreset } from '@/lib/validations/metal-form';
import { RFD3Request, TaskFormProps } from './shared/types';

import { Button } from '@/components/ui/button';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select';
import { Checkbox } from '@/components/ui/checkbox';
import { Badge } from '@/components/ui/badge';
import { Textarea } from '@/components/ui/textarea';
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
  id: MetalLigandApproach;
  name: string;
  description: string;
  details: string;
  recommended?: boolean;
}> = [
  {
    id: 'joint',
    name: 'Joint Design',
    description: 'Both chains simultaneously',
    details: 'RFD3 generates both protein chains around the metal-ligand complex in parallel.',
    recommended: true,
  },
  {
    id: 'sequential',
    name: 'Sequential Design',
    description: 'Chain A first, then B',
    details: 'Design Chain A first, then design Chain B using A + complex as fixed context.',
  },
];

// Expected Results thresholds
const EXPECTED_RESULTS = {
  excellent: {
    coordination: '>= 6 total donors',
    clashes: '0 clashes',
    chainBalance: '3+3 split',
    geometry: 'Correct geometry',
  },
  good: {
    coordination: '>= 4 total donors',
    clashes: '<= 3 clashes',
    chainBalance: '2+2 or better',
    geometry: 'Near correct',
  },
  minimum: {
    coordination: '>= 2 total donors',
    clashes: '<= 6 clashes',
    chainBalance: '1+1 minimum',
    hbonds: 'Ligand H-bonds present',
  },
};

export function InterfaceMetalLigandForm({ onSubmit, isSubmitting, health }: TaskFormProps) {
  const form = useForm<MetalLigandFormValues>({
    resolver: zodResolver(metalLigandFormSchema),
    defaultValues: metalLigandFormDefaults,
  });

  const watchTemplateId = form.watch('templateId');
  const watchUseCustomComplex = form.watch('useCustomComplex');
  const watchApproach = form.watch('approach');
  const watchQualityPreset = form.watch('qualityPreset');
  const watchUseHsabBias = form.watch('useHsabBias');

  const selectedTemplate = getTemplateById(watchTemplateId);

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

  // Template selection
  const handleTemplateSelect = (templateId: string) => {
    form.setValue('templateId', templateId);
    form.setValue('useCustomComplex', false);

    const template = getTemplateById(templateId);
    if (template) {
      // Update coordination split based on template
      const perChain = Math.ceil(template.proteinSitesNeeded / 2);
      form.setValue('chainADonors', perChain);
      form.setValue('chainBDonors', template.proteinSitesNeeded - perChain);
      form.setValue('targetCoordination', template.coordination);

      // Update HSAB bias
      if (form.watch('useHsabBias')) {
        form.setValue('customBiasAA', getHsabBiasForMetal(template.metal));
      }
    }
  };

  // Form submission
  const handleSubmit = async (data: MetalLigandFormValues) => {
    const template = getTemplateById(data.templateId);

    const request: RFD3Request = {
      task: 'interface_metal_ligand_design',
      approach: data.approach,
      contig_str: data.chainLength + ',/0,' + data.chainLength,
      num_designs: data.numDesigns,
      num_timesteps: data.numTimesteps,
      step_scale: data.stepScale,
      gamma_0: data.gamma0,
      validate_coordination: data.validateCoordination,
    };

    // Add template or custom complex info
    if (data.useCustomComplex) {
      if (data.customPdb) {
        request.complex_pdb = data.customPdb;
      } else {
        request.ligand_smiles = data.customSmiles;
        request.metal = data.customMetal;
        request.ligand_name = data.customLigandName;
      }
    } else {
      request.template_name = data.templateId;
    }

    // Add coordination parameters
    request.chain_a_donors = data.chainADonors;
    request.chain_b_donors = data.chainBDonors;

    // Add HSAB bias if enabled
    if (data.useHsabBias && data.customBiasAA) {
      request.bias_AA = data.customBiasAA;
    }

    // Add seed if provided
    if (data.seed) {
      request.seed = parseInt(data.seed, 10);
    }

    await onSubmit(request);
  };

  return (
    <Form {...form}>
      <form onSubmit={form.handleSubmit(handleSubmit)} className="space-y-6">
        {/* Header */}
        <div className="flex items-center gap-3 pb-4 border-b border-border">
          <div className="w-10 h-10 rounded-lg bg-primary flex items-center justify-center">
            <Sparkles className="w-5 h-5 text-primary-foreground" />
          </div>
          <div>
            <h2 className="font-semibold text-foreground">Metal-Ligand Complex Dimer Design</h2>
            <p className="text-sm text-muted-foreground">Design homodimers that bind metal-ligand complexes at the interface</p>
          </div>
        </div>

        {/* Info Banner */}
        <Card className="bg-muted border-border">
          <CardContent className="flex items-start gap-3 pt-4">
            <Info className="w-5 h-5 text-primary mt-0.5 flex-shrink-0" />
            <div className="text-sm text-foreground">
              <p className="font-medium mb-1">Metal-Ligand Complex Interface Design</p>
              <p className="text-muted-foreground">
                This designs protein homodimers where the metal-ligand complex (e.g., citrate-Tb) sits at the interface.
                Both chains coordinate the metal, with the ligand providing additional coordination sites.
                Uses HSAB chemistry for optimal binding.
              </p>
            </div>
          </CardContent>
        </Card>

        {/* Complex Template Selection */}
        <Card>
          <CardHeader>
            <CardTitle className="text-sm flex items-center gap-2">
              <Atom className="w-4 h-4" />
              Metal-Ligand Complex
              <Badge variant="outline" className="text-xs">Required</Badge>
            </CardTitle>
            <CardDescription>Select a pre-defined complex or provide custom PDB/SMILES</CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            {/* Template buttons */}
            <div className="grid grid-cols-1 sm:grid-cols-2 gap-3">
              {METAL_LIGAND_TEMPLATES.map((template) => (
                <button
                  key={template.id}
                  type="button"
                  onClick={() => handleTemplateSelect(template.id)}
                  className={cn(
                    "p-4 rounded-lg border text-left transition-all",
                    watchTemplateId === template.id && !watchUseCustomComplex
                      ? "border-primary bg-primary/5 ring-1 ring-primary/20"
                      : "border-border hover:border-muted-foreground/30 hover:bg-muted/50"
                  )}
                >
                  <div className="flex items-center justify-between mb-2">
                    <span className="font-medium text-foreground">{template.name}</span>
                    <Badge variant="secondary" className="text-xs">
                      {template.metal}
                    </Badge>
                  </div>
                  <p className="text-xs text-muted-foreground mb-2">{template.description}</p>
                  <div className="flex flex-wrap gap-2 text-xs">
                    <span className="px-2 py-0.5 rounded bg-muted text-muted-foreground">
                      CN: {template.coordination}
                    </span>
                    <span className="px-2 py-0.5 rounded bg-muted text-muted-foreground">
                      Ligand: {template.ligandDonors} donors
                    </span>
                    <span className="px-2 py-0.5 rounded bg-muted text-muted-foreground">
                      Protein: {template.proteinSitesNeeded} needed
                    </span>
                  </div>
                  {template.luminescence && (
                    <div className="mt-2 flex items-center gap-1 text-xs text-primary">
                      <Sparkles className="w-3 h-3" />
                      {template.luminescence}
                    </div>
                  )}
                </button>
              ))}

              {/* Custom option */}
              <button
                type="button"
                onClick={() => form.setValue('useCustomComplex', true)}
                className={cn(
                  "p-4 rounded-lg border text-left transition-all",
                  watchUseCustomComplex
                    ? "border-primary bg-primary/5 ring-1 ring-primary/20"
                    : "border-border hover:border-muted-foreground/30 hover:bg-muted/50"
                )}
              >
                <div className="flex items-center gap-2 mb-2">
                  <Upload className="w-4 h-4 text-muted-foreground" />
                  <span className="font-medium text-foreground">Custom Complex</span>
                </div>
                <p className="text-xs text-muted-foreground">
                  Provide your own PDB or SMILES for the metal-ligand complex
                </p>
              </button>
            </div>

            {/* Custom complex fields */}
            {watchUseCustomComplex && (
              <div className="space-y-4 p-4 rounded-lg bg-muted/50 border">
                <div className="grid grid-cols-1 sm:grid-cols-2 gap-4">
                  <FormField
                    control={form.control}
                    name="customMetal"
                    render={({ field }) => (
                      <FormItem>
                        <FormLabel>Metal Code</FormLabel>
                        <FormControl>
                          <Input placeholder="TB, EU, CA..." {...field} />
                        </FormControl>
                        <FormMessage />
                      </FormItem>
                    )}
                  />

                  <FormField
                    control={form.control}
                    name="customLigandName"
                    render={({ field }) => (
                      <FormItem>
                        <FormLabel>Ligand Name</FormLabel>
                        <FormControl>
                          <Input placeholder="CIT, PQQ, LIG..." {...field} />
                        </FormControl>
                        <FormDescription className="text-xs">3-letter residue code</FormDescription>
                        <FormMessage />
                      </FormItem>
                    )}
                  />
                </div>

                <FormField
                  control={form.control}
                  name="customSmiles"
                  render={({ field }) => (
                    <FormItem>
                      <FormLabel>Ligand SMILES</FormLabel>
                      <FormControl>
                        <Input
                          {...field}
                          placeholder="OC(CC(=O)[O-])(CC(=O)[O-])C(=O)[O-]"
                          className="font-mono"
                        />
                      </FormControl>
                      <FormDescription className="text-xs">
                        SMILES string for the ligand (without metal)
                      </FormDescription>
                      <FormMessage />
                    </FormItem>
                  )}
                />

                <div className="relative">
                  <div className="absolute inset-0 flex items-center">
                    <span className="w-full border-t" />
                  </div>
                  <div className="relative flex justify-center text-xs uppercase">
                    <span className="bg-muted px-2 text-muted-foreground">Or upload PDB</span>
                  </div>
                </div>

                <FormField
                  control={form.control}
                  name="customPdb"
                  render={({ field }) => (
                    <FormItem>
                      <FormLabel>Complex PDB Content</FormLabel>
                      <FormControl>
                        <Textarea
                          {...field}
                          placeholder="HETATM    1  C1  CIT L   1    ..."
                          className="font-mono h-32"
                        />
                      </FormControl>
                      <FormDescription className="text-xs">
                        PDB content with both HETATM records for ligand and metal
                      </FormDescription>
                      <FormMessage />
                    </FormItem>
                  )}
                />
              </div>
            )}

            {/* Selected template info */}
            {!watchUseCustomComplex && selectedTemplate && (
              <div className="p-3 rounded-lg bg-muted/30 border text-sm">
                <div className="flex items-center gap-2 mb-2">
                  <FlaskConical className="w-4 h-4 text-muted-foreground" />
                  <span className="font-medium">{selectedTemplate.name} Coordination</span>
                </div>
                <div className="grid grid-cols-2 gap-2 text-xs text-muted-foreground">
                  <div>Metal: <span className="font-medium text-foreground">{selectedTemplate.metal}</span></div>
                  <div>Ligand: <span className="font-medium text-foreground">{selectedTemplate.ligand}</span></div>
                  <div>Total CN: <span className="font-medium text-foreground">{selectedTemplate.coordination}</span></div>
                  <div>From ligand: <span className="font-medium text-foreground">{selectedTemplate.ligandDonors}</span></div>
                  <div className="col-span-2">Protein must provide: <span className="font-medium text-foreground">{selectedTemplate.proteinSitesNeeded} donors</span></div>
                </div>
              </div>
            )}
          </CardContent>
        </Card>

        {/* Design Approach */}
        <Card>
          <CardHeader>
            <CardTitle className="text-sm flex items-center gap-2">
              Design Approach
              <Badge variant="outline" className="text-xs">Required</Badge>
            </CardTitle>
            <CardDescription>Choose how to design the homodimer</CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            <div className="grid grid-cols-1 sm:grid-cols-2 gap-3">
              {APPROACHES.map((app) => (
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
                    <Network className={cn(
                      "w-4 h-4",
                      watchApproach === app.id ? "text-primary" : "text-muted-foreground"
                    )} />
                    <span className="font-medium text-foreground text-sm">{app.name}</span>
                  </div>
                  <p className="text-xs text-muted-foreground">{app.description}</p>
                  <p className="text-xs text-muted-foreground/70 mt-1">{app.details}</p>
                </button>
              ))}
            </div>
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

            {/* Coordination Split */}
            <div className="grid grid-cols-1 sm:grid-cols-2 gap-4">
              <FormField
                control={form.control}
                name="chainADonors"
                render={({ field }) => (
                  <FormItem>
                    <FormLabel>Chain A Donors</FormLabel>
                    <FormControl>
                      <Input
                        type="number"
                        min={1}
                        max={6}
                        {...field}
                        onChange={(e) => field.onChange(parseInt(e.target.value) || 1)}
                      />
                    </FormControl>
                    <FormDescription className="text-xs">Target coordination from Chain A</FormDescription>
                    <FormMessage />
                  </FormItem>
                )}
              />

              <FormField
                control={form.control}
                name="chainBDonors"
                render={({ field }) => (
                  <FormItem>
                    <FormLabel>Chain B Donors</FormLabel>
                    <FormControl>
                      <Input
                        type="number"
                        min={1}
                        max={6}
                        {...field}
                        onChange={(e) => field.onChange(parseInt(e.target.value) || 1)}
                      />
                    </FormControl>
                    <FormDescription className="text-xs">Target coordination from Chain B</FormDescription>
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
                <Atom className="w-4 h-4" />
                Advanced Options
              </div>
              <ChevronDown className="w-4 h-4" />
            </Button>
          </CollapsibleTrigger>
          <CollapsibleContent className="mt-3 space-y-4 p-4 border rounded-lg">
            {/* HSAB Bias */}
            <FormField
              control={form.control}
              name="useHsabBias"
              render={({ field }) => (
                <FormItem className="flex flex-row items-start space-x-3 space-y-0">
                  <FormControl>
                    <Checkbox checked={field.value} onCheckedChange={field.onChange} />
                  </FormControl>
                  <div className="space-y-1 leading-none">
                    <FormLabel>Use HSAB Chemistry Bias</FormLabel>
                    <FormDescription className="text-xs">
                      Bias LigandMPNN toward amino acids optimal for this metal (recommended)
                    </FormDescription>
                  </div>
                </FormItem>
              )}
            />

            {watchUseHsabBias && (
              <FormField
                control={form.control}
                name="customBiasAA"
                render={({ field }) => (
                  <FormItem>
                    <FormLabel className="text-xs">Amino Acid Bias</FormLabel>
                    <FormControl>
                      <Input
                        {...field}
                        placeholder="E:3.0,D:3.0,N:2.4,Q:2.4,H:0.5,C:-5.0"
                        className="font-mono"
                      />
                    </FormControl>
                    <FormDescription className="text-xs">
                      Format: AA:weight (positive = favor, negative = avoid)
                    </FormDescription>
                    <FormMessage />
                  </FormItem>
                )}
              />
            )}

            {/* Validation Options */}
            <div className="space-y-3">
              <FormField
                control={form.control}
                name="validateCoordination"
                render={({ field }) => (
                  <FormItem className="flex flex-row items-start space-x-3 space-y-0">
                    <FormControl>
                      <Checkbox checked={field.value} onCheckedChange={field.onChange} />
                    </FormControl>
                    <div className="space-y-1 leading-none">
                      <FormLabel>Validate Coordination</FormLabel>
                      <FormDescription className="text-xs">
                        Check metal coordination after design
                      </FormDescription>
                    </div>
                  </FormItem>
                )}
              />

              <FormField
                control={form.control}
                name="runPipelineValidation"
                render={({ field }) => (
                  <FormItem className="flex flex-row items-start space-x-3 space-y-0">
                    <FormControl>
                      <Checkbox checked={field.value} onCheckedChange={field.onChange} />
                    </FormControl>
                    <div className="space-y-1 leading-none">
                      <FormLabel>Run Pipeline Validation</FormLabel>
                      <FormDescription className="text-xs">
                        Full validation including clash detection and geometry analysis
                      </FormDescription>
                    </div>
                  </FormItem>
                )}
              />
            </div>

            {/* Custom quality settings */}
            {watchQualityPreset === 'Custom' && (
              <div className="grid grid-cols-3 gap-4">
                <FormField
                  control={form.control}
                  name="numTimesteps"
                  render={({ field }) => (
                    <FormItem>
                      <FormLabel className="text-xs">Timesteps</FormLabel>
                      <FormControl>
                        <Input
                          type="number"
                          min={10}
                          max={500}
                          {...field}
                          onChange={(e) => field.onChange(parseInt(e.target.value) || 200)}
                        />
                      </FormControl>
                    </FormItem>
                  )}
                />

                <FormField
                  control={form.control}
                  name="stepScale"
                  render={({ field }) => (
                    <FormItem>
                      <FormLabel className="text-xs">Step Scale</FormLabel>
                      <FormControl>
                        <Input
                          type="number"
                          step={0.1}
                          min={0.1}
                          max={3.0}
                          {...field}
                          onChange={(e) => field.onChange(parseFloat(e.target.value) || 1.5)}
                        />
                      </FormControl>
                    </FormItem>
                  )}
                />

                <FormField
                  control={form.control}
                  name="gamma0"
                  render={({ field }) => (
                    <FormItem>
                      <FormLabel className="text-xs">Gamma</FormLabel>
                      <FormControl>
                        <Input
                          type="number"
                          step={0.1}
                          min={0}
                          max={1}
                          {...field}
                          onChange={(e) => field.onChange(parseFloat(e.target.value) || 0.6)}
                        />
                      </FormControl>
                    </FormItem>
                  )}
                />
              </div>
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
                  <Sparkles className="w-4 h-4" />
                </span>
                Designing Metal-Ligand Dimer...
              </>
            ) : (
              <>
                <Sparkles className="w-4 h-4 mr-2" />
                Design {selectedTemplate?.name || 'Custom'} Dimer
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
                    <div>Coordination: {EXPECTED_RESULTS.excellent.coordination}</div>
                    <div>Clashes: {EXPECTED_RESULTS.excellent.clashes}</div>
                    <div>Balance: {EXPECTED_RESULTS.excellent.chainBalance}</div>
                    <div>Geometry: {EXPECTED_RESULTS.excellent.geometry}</div>
                  </div>
                </div>

                {/* Good */}
                <div className="p-2 rounded-lg bg-muted/70 border border-border">
                  <div className="flex items-center gap-2 mb-1">
                    <span className="w-2 h-2 rounded-full bg-muted-foreground" />
                    <span className="text-xs font-medium text-foreground">Good</span>
                  </div>
                  <div className="grid grid-cols-2 gap-x-4 gap-y-0.5 text-xs text-muted-foreground ml-4">
                    <div>Coordination: {EXPECTED_RESULTS.good.coordination}</div>
                    <div>Clashes: {EXPECTED_RESULTS.good.clashes}</div>
                    <div>Balance: {EXPECTED_RESULTS.good.chainBalance}</div>
                    <div>Geometry: {EXPECTED_RESULTS.good.geometry}</div>
                  </div>
                </div>

                {/* Minimum */}
                <div className="p-2 rounded-lg bg-muted/50 border border-border">
                  <div className="flex items-center gap-2 mb-1">
                    <span className="w-2 h-2 rounded-full bg-muted-foreground/60" />
                    <span className="text-xs font-medium text-foreground">Minimum Viable</span>
                  </div>
                  <div className="grid grid-cols-2 gap-x-4 gap-y-0.5 text-xs text-muted-foreground ml-4">
                    <div>Coordination: {EXPECTED_RESULTS.minimum.coordination}</div>
                    <div>Clashes: {EXPECTED_RESULTS.minimum.clashes}</div>
                    <div>Balance: {EXPECTED_RESULTS.minimum.chainBalance}</div>
                    <div>H-bonds: {EXPECTED_RESULTS.minimum.hbonds}</div>
                  </div>
                </div>
              </div>

              <p className="text-xs text-muted-foreground mt-3 text-center">
                Pipeline: RFD3 Backbone → LigandMPNN + HSAB → FastRelax → Validation
              </p>
            </CardContent>
          </Card>
        </div>
      </form>
    </Form>
  );
}
