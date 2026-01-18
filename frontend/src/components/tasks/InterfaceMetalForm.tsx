'use client';

import { useMemo } from 'react';
import { useForm } from 'react-hook-form';
import { zodResolver } from '@hookform/resolvers/zod';
import {
  Diamond,
  FlaskConical,
  Info,
  Network,
  SlidersHorizontal,
  Link,
  Zap,
  GitBranch,
  ChevronDown,
  Minus,
  Plus,
  Star,
  Hand,
  Droplet,
  Expand,
  Settings2,
  Sparkles,
  CheckCircle,
} from 'lucide-react';

import { cn } from '@/lib/utils';
import {
  metalFormSchema,
  metalFormDefaults,
  QUALITY_PRESET_VALUES,
  type MetalFormValues,
  type MetalId,
  type MetalApproach,
  type QualityPreset,
  type TemplateLibraryId,
} from '@/lib/validations/metal-form';
import { RFD3Request, TaskFormProps } from './shared/types';

import { Button } from '@/components/ui/button';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select';
import { Slider } from '@/components/ui/slider';
import { Switch } from '@/components/ui/switch';
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

// Metal ion profiles with coordination chemistry
const METAL_PROFILES = {
  ZN: {
    name: 'Zinc (Zn2+)',
    category: 'Transition',
    coordination: [4, 5, 6],
    preferredCoord: 4,
    geometry: 'tetrahedral',
    donorLabels: ['His', 'Cys', 'Asp', 'Glu'],
    distance: [2.0, 2.3],
    description: 'Most common catalytic metal. Tetrahedral preferred.',
  },
  FE: {
    name: 'Iron (Fe2+)',
    category: 'Transition',
    coordination: [4, 5, 6],
    preferredCoord: 6,
    geometry: 'octahedral',
    donorLabels: ['His', 'Cys', 'Asp', 'Glu', 'Tyr'],
    distance: [1.9, 2.2],
    description: 'Redox-active. Found in heme and non-heme sites.',
  },
  CU: {
    name: 'Copper (Cu2+)',
    category: 'Transition',
    coordination: [4, 5, 6],
    preferredCoord: 4,
    geometry: 'square_planar',
    donorLabels: ['His', 'Cys', 'Met'],
    distance: [1.9, 2.1],
    description: 'Electron transfer. Square planar or distorted.',
  },
  MN: {
    name: 'Manganese (Mn2+)',
    category: 'Transition',
    coordination: [6],
    preferredCoord: 6,
    geometry: 'octahedral',
    donorLabels: ['His', 'Asp', 'Glu', 'Water'],
    distance: [2.1, 2.3],
    description: 'Catalytic. Often in oxygen-evolving complexes.',
  },
  CA: {
    name: 'Calcium (Ca2+)',
    category: 'Alkaline Earth',
    coordination: [6, 7, 8],
    preferredCoord: 7,
    geometry: 'pentagonal_bipyramidal',
    donorLabels: ['Asp', 'Glu', 'Asn', 'Backbone O', 'Water'],
    distance: [2.3, 2.5],
    description: 'Structural/signaling. Flexible coordination.',
  },
  MG: {
    name: 'Magnesium (Mg2+)',
    category: 'Alkaline Earth',
    coordination: [6],
    preferredCoord: 6,
    geometry: 'octahedral',
    donorLabels: ['Asp', 'Glu', 'Water'],
    distance: [2.0, 2.2],
    description: 'Catalytic. Strict octahedral geometry.',
  },
  TB: {
    name: 'Terbium (Tb3+)',
    category: 'Lanthanide',
    coordination: [8, 9],
    preferredCoord: 9,
    geometry: 'tricapped_trigonal_prismatic',
    donorLabels: ['Asp', 'Glu', 'Asn', 'Backbone O', 'Water'],
    distance: [2.3, 2.5],
    description: 'Luminescent (green). Antenna effect with Trp.',
    special: 'luminescent' as const,
  },
  GD: {
    name: 'Gadolinium (Gd3+)',
    category: 'Lanthanide',
    coordination: [8, 9],
    preferredCoord: 9,
    geometry: 'tricapped_trigonal_prismatic',
    donorLabels: ['Asp', 'Glu', 'Asn', 'Backbone O', 'Water'],
    distance: [2.3, 2.5],
    description: 'Paramagnetic. Used in MRI contrast agents.',
    special: 'paramagnetic' as const,
  },
  EU: {
    name: 'Europium (Eu3+)',
    category: 'Lanthanide',
    coordination: [8, 9],
    preferredCoord: 9,
    geometry: 'tricapped_trigonal_prismatic',
    donorLabels: ['Asp', 'Glu', 'Asn', 'Backbone O', 'Water'],
    distance: [2.3, 2.5],
    description: 'Luminescent (red). Long-lived emission.',
    special: 'luminescent' as const,
  },
} as const;

// Template library options
const TEMPLATE_LIBRARY_OPTIONS = [
  {
    id: 'caldwell_4' as const,
    name: 'Caldwell (4 Glu)',
    description: '4 bidentate Glu, CN=8 - highest affinity (subfemtomolar)',
    coordinationNumber: 8,
    geometry: 'square_antiprism',
    recommendedFor: ['TB', 'EU'],
    icon: Star,
  },
  {
    id: 'ef_hand_8' as const,
    name: 'EF-Hand Mixed',
    description: '4 mono Asp + 2 bi Glu, CN=8 - balanced design',
    coordinationNumber: 8,
    geometry: 'distorted_square_antiprism',
    recommendedFor: ['GD'],
    icon: Hand,
  },
  {
    id: 'lanm_mixed' as const,
    name: 'Lanmodulin',
    description: '3 bi Asp + 2 waters, CN=9 - natural-like',
    coordinationNumber: 9,
    geometry: 'tricapped_trigonal_prism',
    recommendedFor: ['TB', 'GD'],
    icon: Droplet,
  },
  {
    id: 'high_coord_9' as const,
    name: 'High Coordination',
    description: '4 bi Glu + 1 water, CN=9 - for large lanthanides',
    coordinationNumber: 9,
    geometry: 'tricapped_trigonal_prism',
    recommendedFor: ['CA'],
    icon: Expand,
  },
];

// Design approaches
const APPROACHES: Array<{
  id: MetalApproach;
  name: string;
  description: string;
  details: string;
  icon: typeof Network;
  recommended?: boolean;
  advanced?: boolean;
}> = [
  {
    id: 'joint_metal',
    name: 'Joint Metal Dimer',
    description: 'Both chains co-evolve around metal',
    details: 'Simultaneous design of two chains sharing coordination.',
    icon: Network,
    recommended: true,
  },
  {
    id: 'asymmetric_metal',
    name: 'Asymmetric Metal Dimer',
    description: 'Different donor types per chain',
    details: 'Chain A provides some donor types, Chain B provides others.',
    icon: SlidersHorizontal,
  },
  {
    id: 'induced_metal',
    name: 'Metal-Induced Dimerization',
    description: 'Metal binding drives assembly',
    details: 'Monomers only become stable when dimerized.',
    icon: GitBranch,
  },
  {
    id: 'bridging_metal',
    name: 'Bridging Metals',
    description: 'Multiple metals at interface',
    details: 'Two or more metal ions bridging chains.',
    icon: Link,
    advanced: true,
  },
  {
    id: 'redox_switch',
    name: 'Redox-Switchable',
    description: 'Oxidation state controls dimer',
    details: 'Different affinity based on metal oxidation state.',
    icon: Zap,
    advanced: true,
  },
];

// Geometry options based on coordination number
const GEOMETRY_OPTIONS_BY_COORD: Record<number, string[]> = {
  4: ['tetrahedral', 'square_planar'],
  5: ['trigonal_bipyramidal', 'square_pyramidal'],
  6: ['octahedral'],
  7: ['pentagonal_bipyramidal', 'capped_octahedral'],
  8: ['square_antiprismatic', 'dodecahedral'],
  9: ['tricapped_trigonal_prismatic'],
};

export function InterfaceMetalForm({ onSubmit, isSubmitting, health }: TaskFormProps) {
  const form = useForm<MetalFormValues>({
    resolver: zodResolver(metalFormSchema),
    defaultValues: metalFormDefaults,
  });

  const watchMetal = form.watch('metal');
  const watchApproach = form.watch('approach');
  const watchChainADonors = form.watch('chainADonors');
  const watchChainBDonors = form.watch('chainBDonors');
  const watchQualityPreset = form.watch('qualityPreset');
  const watchUseParametricMode = form.watch('useParametricMode');
  const watchUseLegacyMode = form.watch('useLegacyMode');
  const watchTemplateName = form.watch('templateName');
  const watchParamCoordinationNumber = form.watch('paramCoordinationNumber');
  const watchParamNumWaters = form.watch('paramNumWaters');

  const metalProfile = METAL_PROFILES[watchMetal];
  const totalCoordination = watchChainADonors + watchChainBDonors;
  const geometryOptions = GEOMETRY_OPTIONS_BY_COORD[totalCoordination] || [];
  const isLanthanide = metalProfile.category === 'Lanthanide';
  const isLuminescent = 'special' in metalProfile && metalProfile.special === 'luminescent';

  // Update coordination when metal changes
  const handleMetalChange = (metalId: MetalId) => {
    form.setValue('metal', metalId);
    const profile = METAL_PROFILES[metalId];
    const coord = profile.preferredCoord;
    const split = Math.floor(coord / 2);
    form.setValue('chainADonors', split);
    form.setValue('chainBDonors', coord - split);
    form.setValue('targetGeometry', 'auto');

    // Reset TEBL options when switching away from luminescent metals
    const isLum = 'special' in profile && profile.special === 'luminescent';
    if (!isLum) {
      form.setValue('addTrpAntenna', false);
    }
  };

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

  // Mode switches for lanthanide templates
  const switchToTemplateLibrary = () => {
    form.setValue('useParametricMode', false);
    form.setValue('useLegacyMode', false);
  };

  const switchToParametric = () => {
    form.setValue('useParametricMode', true);
    form.setValue('useLegacyMode', false);
  };

  const switchToLegacy = () => {
    form.setValue('useParametricMode', false);
    form.setValue('useLegacyMode', true);
  };

  // Template change syncs coordination
  const handleTemplateChange = (templateId: TemplateLibraryId) => {
    form.setValue('templateName', templateId);
    const template = TEMPLATE_LIBRARY_OPTIONS.find(t => t.id === templateId);
    if (template) {
      const cn = template.coordinationNumber;
      const split = Math.floor(cn / 2);
      form.setValue('chainADonors', split);
      form.setValue('chainBDonors', cn - split);
    }
  };

  // Form submission
  const handleSubmit = async (data: MetalFormValues) => {
    const request: RFD3Request = {
      task: 'interface_metal_design',
      approach: data.approach,
      ligand: data.metal,
      metal: data.metal,
      chain_length: data.chainLength,
      num_designs: data.numDesigns,
      num_timesteps: data.numTimesteps,
      step_scale: data.stepScale,
      gamma_0: data.gamma0,
    };

    request.metal_config = {
      coordination_split: [data.chainADonors, data.chainBDonors],
      geometry: data.targetGeometry === 'auto' ? metalProfile.geometry : data.targetGeometry,
      chain_a_donors: data.chainADonorTypes,
      chain_b_donors: data.chainBDonorTypes,
      include_waters: data.includeWaters,
    };

    // Lanthanide-specific
    if (isLanthanide) {
      if (data.useParametricMode) {
        request.parametric = {
          coordination_number: data.paramCoordinationNumber,
          num_waters: data.paramNumWaters,
          bidentate_fraction: data.paramBidentateFraction,
        };
      } else if (data.useLegacyMode || data.templateName === 'legacy') {
        request.template_type = data.templateType;
        request.donor_residue = data.donorResidue;
        request.use_motif_scaffolding = data.useMotifScaffolding;
      } else {
        request.template_name = data.templateName;
      }
      request.add_trp_antenna = data.addTrpAntenna;
      request.validate_coordination = data.validateCoordination;
    }

    if (data.approach === 'bridging_metal' && data.secondMetal) {
      request.metal_config.second_metal = data.secondMetal;
    }

    if (data.seed) {
      request.seed = parseInt(data.seed, 10);
    }

    await onSubmit(request);
  };

  const isValidCoordination = (metalProfile.coordination as readonly number[]).includes(totalCoordination);

  return (
    <Form {...form}>
      <form onSubmit={form.handleSubmit(handleSubmit)} className="space-y-6">
        {/* Header */}
        <div className="flex items-center gap-3 pb-4 border-b border-border">
          <div className="w-10 h-10 rounded-lg bg-primary flex items-center justify-center">
            <Diamond className="w-5 h-5 text-primary-foreground" />
          </div>
          <div>
            <h2 className="font-semibold text-foreground">Interface Metal Dimer Design</h2>
            <p className="text-sm text-muted-foreground">Design protein dimers with metal coordination at the interface</p>
          </div>
        </div>

        {/* Info Banner */}
        <Card className="bg-muted border-border">
          <CardContent className="flex items-start gap-3 pt-4">
            <FlaskConical className="w-5 h-5 text-primary mt-0.5 flex-shrink-0" />
            <div className="text-sm text-foreground">
              <p className="font-medium mb-1">Metal-Mediated Dimerization</p>
              <p className="text-muted-foreground">
                Unlike organic ligands, metals have fixed coordination geometries. Each chain contributes
                donor atoms (His, Cys, Asp, Glu) that together complete the metal&apos;s coordination sphere.
              </p>
            </div>
          </CardContent>
        </Card>

        {/* Metal Selection */}
        <Card>
          <CardHeader>
            <CardTitle className="text-sm flex items-center gap-2">
              Metal Ion
              <Badge variant="outline" className="text-xs">Required</Badge>
            </CardTitle>
            <CardDescription>Select the metal for interface coordination</CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            {/* Category tabs */}
            {['Transition', 'Alkaline Earth', 'Lanthanide'].map((category) => (
              <div key={category} className="space-y-2">
                <span className="text-xs font-medium text-muted-foreground uppercase tracking-wide">{category}</span>
                <div className="flex gap-2 flex-wrap">
                  {Object.entries(METAL_PROFILES)
                    .filter(([, profile]) => profile.category === category)
                    .map(([id, profile]) => (
                      <Button
                        key={id}
                        type="button"
                        variant={watchMetal === id ? 'default' : 'outline'}
                        size="sm"
                        onClick={() => handleMetalChange(id as MetalId)}
                      >
                        <span className="font-bold">{id}</span>
                        {'special' in profile && profile.special === 'luminescent' && (
                          <Sparkles className="w-3 h-3 ml-1" />
                        )}
                      </Button>
                    ))}
                </div>
              </div>
            ))}

            {/* Selected metal info */}
            <div className="p-3 rounded-lg bg-muted/50 border">
              <div className="flex items-center gap-2 mb-2">
                <span className="text-lg font-bold text-primary">{watchMetal}</span>
                <span className="text-sm text-muted-foreground">{metalProfile.name}</span>
              </div>
              <p className="text-xs text-muted-foreground mb-2">{metalProfile.description}</p>
              <div className="grid grid-cols-3 gap-2 text-xs">
                <div>
                  <span className="text-muted-foreground">Coordination: </span>
                  <span className="font-medium">{metalProfile.coordination.join('-')}</span>
                </div>
                <div>
                  <span className="text-muted-foreground">Geometry: </span>
                  <span className="font-medium capitalize">{metalProfile.geometry.replace('_', ' ')}</span>
                </div>
                <div>
                  <span className="text-muted-foreground">Distance: </span>
                  <span className="font-medium">{metalProfile.distance[0]}-{metalProfile.distance[1]} A</span>
                </div>
              </div>
              <div className="mt-2 flex flex-wrap gap-1">
                <span className="text-xs text-muted-foreground">Donors: </span>
                {metalProfile.donorLabels.map((donor) => (
                  <Badge key={donor} variant="secondary" className="text-xs">
                    {donor}
                  </Badge>
                ))}
              </div>
            </div>
          </CardContent>
        </Card>

        {/* Design Approach */}
        <Card>
          <CardHeader>
            <CardTitle className="text-sm flex items-center gap-2">
              Design Approach
              <Badge variant="outline" className="text-xs">Required</Badge>
            </CardTitle>
            <CardDescription>How to coordinate the metal between chains</CardDescription>
          </CardHeader>
          <CardContent>
            <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-3">
              {APPROACHES.filter(a => !a.advanced).map((app) => {
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
          </CardContent>
        </Card>

        {/* Coordination Split */}
        <Card>
          <CardHeader>
            <CardTitle className="text-sm">Coordination Split</CardTitle>
            <CardDescription>How many donor atoms from each chain</CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            <div className="flex items-center gap-6 justify-center">
              <div className="text-center">
                <Label className="text-xs text-muted-foreground mb-1 block">Chain A Donors</Label>
                <div className="flex items-center justify-center gap-2">
                  <Button
                    type="button"
                    variant="outline"
                    size="icon"
                    className="h-8 w-8"
                    onClick={() => form.setValue('chainADonors', Math.max(1, watchChainADonors - 1))}
                  >
                    <Minus className="h-3 w-3" />
                  </Button>
                  <span className="w-12 h-12 rounded-xl bg-muted border-2 border-primary text-primary flex items-center justify-center text-2xl font-bold">
                    {watchChainADonors}
                  </span>
                  <Button
                    type="button"
                    variant="outline"
                    size="icon"
                    className="h-8 w-8"
                    onClick={() => form.setValue('chainADonors', watchChainADonors + 1)}
                  >
                    <Plus className="h-3 w-3" />
                  </Button>
                </div>
              </div>

              <div className="flex flex-col items-center">
                <span className="text-2xl text-muted-foreground/50">+</span>
                <Diamond className="w-5 h-5 text-primary" />
              </div>

              <div className="text-center">
                <Label className="text-xs text-muted-foreground mb-1 block">Chain B Donors</Label>
                <div className="flex items-center justify-center gap-2">
                  <Button
                    type="button"
                    variant="outline"
                    size="icon"
                    className="h-8 w-8"
                    onClick={() => form.setValue('chainBDonors', Math.max(1, watchChainBDonors - 1))}
                  >
                    <Minus className="h-3 w-3" />
                  </Button>
                  <span className="w-12 h-12 rounded-xl bg-muted border-2 border-muted-foreground text-muted-foreground flex items-center justify-center text-2xl font-bold">
                    {watchChainBDonors}
                  </span>
                  <Button
                    type="button"
                    variant="outline"
                    size="icon"
                    className="h-8 w-8"
                    onClick={() => form.setValue('chainBDonors', watchChainBDonors + 1)}
                  >
                    <Plus className="h-3 w-3" />
                  </Button>
                </div>
              </div>

              <div className="flex flex-col items-center">
                <span className="text-2xl text-muted-foreground/50">=</span>
              </div>

              <div className="text-center">
                <Label className="text-xs text-muted-foreground mb-1 block">Total</Label>
                <span className={cn(
                  "w-12 h-12 rounded-xl flex items-center justify-center text-2xl font-bold",
                  isValidCoordination
                    ? "bg-emerald-100 dark:bg-emerald-900/30 text-emerald-700 dark:text-emerald-300"
                    : "bg-red-100 dark:bg-red-900/30 text-red-700 dark:text-red-300"
                )}>
                  {totalCoordination}
                </span>
              </div>
            </div>

            {/* Validation message */}
            {!isValidCoordination && (
              <div className="p-2 rounded-lg bg-amber-50 dark:bg-amber-950/20 border border-amber-200 dark:border-amber-800 text-xs text-amber-700 dark:text-amber-300">
                <Info className="w-3 h-3 inline mr-1" />
                {watchMetal} typically has {metalProfile.coordination.join(' or ')}-coordinate geometry.
                Total {totalCoordination} may not be optimal.
              </div>
            )}

            {/* Geometry selector */}
            {geometryOptions.length > 0 && (
              <FormField
                control={form.control}
                name="targetGeometry"
                render={({ field }) => (
                  <FormItem>
                    <FormLabel className="text-xs">Target Geometry</FormLabel>
                    <div className="flex gap-2 flex-wrap">
                      <Button
                        type="button"
                        variant={field.value === 'auto' ? 'default' : 'secondary'}
                        size="sm"
                        onClick={() => field.onChange('auto')}
                      >
                        Auto ({metalProfile.geometry.replace('_', ' ')})
                      </Button>
                      {geometryOptions.map((geom) => (
                        <Button
                          key={geom}
                          type="button"
                          variant={field.value === geom ? 'default' : 'secondary'}
                          size="sm"
                          onClick={() => field.onChange(geom)}
                        >
                          {geom.replace('_', ' ')}
                        </Button>
                      ))}
                    </div>
                  </FormItem>
                )}
              />
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
                        max={20}
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

        {/* Lanthanide-Specific Options */}
        {isLanthanide && (
          <Card>
            <CardHeader>
              <CardTitle className="text-sm">Lanthanide Design Options</CardTitle>
              <CardDescription>Chemically realistic templates for optimal coordination (CN=8-9)</CardDescription>
            </CardHeader>
            <CardContent className="space-y-4">
              {/* Mode Selection */}
              <div className="flex gap-1 p-1 bg-muted rounded-lg">
                <Button
                  type="button"
                  variant={!watchUseParametricMode && !watchUseLegacyMode ? 'secondary' : 'ghost'}
                  size="sm"
                  className="flex-1"
                  onClick={switchToTemplateLibrary}
                >
                  Template Library
                </Button>
                <Button
                  type="button"
                  variant={watchUseParametricMode ? 'secondary' : 'ghost'}
                  size="sm"
                  className="flex-1"
                  onClick={switchToParametric}
                >
                  Parametric
                </Button>
                <Button
                  type="button"
                  variant={watchUseLegacyMode ? 'secondary' : 'ghost'}
                  size="sm"
                  className="flex-1"
                  onClick={switchToLegacy}
                >
                  Legacy
                </Button>
              </div>

              {/* Template Library */}
              {!watchUseParametricMode && !watchUseLegacyMode && (
                <div className="space-y-2">
                  <Label className="text-xs flex items-center gap-2">
                    Select Template
                    <Badge variant="secondary" className="text-xs">New!</Badge>
                  </Label>
                  <div className="grid grid-cols-1 sm:grid-cols-2 gap-3">
                    {TEMPLATE_LIBRARY_OPTIONS.map((template) => {
                      const Icon = template.icon;
                      const isRecommended = template.recommendedFor.includes(watchMetal);
                      return (
                        <button
                          key={template.id}
                          type="button"
                          onClick={() => handleTemplateChange(template.id)}
                          className={cn(
                            "p-3 rounded-lg border text-left transition-all",
                            watchTemplateName === template.id
                              ? "border-primary bg-primary/5 ring-1 ring-primary/20"
                              : "border-border hover:border-muted-foreground/30 hover:bg-muted/50"
                          )}
                        >
                          <div className="flex items-center gap-2 mb-1">
                            <Icon className={cn(
                              "w-4 h-4",
                              watchTemplateName === template.id ? "text-primary" : "text-muted-foreground"
                            )} />
                            <span className="font-medium text-sm">{template.name}</span>
                            {isRecommended && (
                              <Badge variant="secondary" className="text-xs">
                                Best for {watchMetal}
                              </Badge>
                            )}
                          </div>
                          <p className="text-xs text-muted-foreground">{template.description}</p>
                          <div className="flex items-center gap-2 mt-1">
                            <span className="text-xs text-primary">CN={template.coordinationNumber}</span>
                            <span className="text-xs text-muted-foreground">{template.geometry}</span>
                          </div>
                        </button>
                      );
                    })}
                  </div>
                </div>
              )}

              {/* Parametric Mode */}
              {watchUseParametricMode && (
                <div className="space-y-4 p-4 rounded-lg bg-muted border border-border">
                  <div className="flex items-center gap-2 mb-3">
                    <Settings2 className="w-4 h-4 text-primary" />
                    <span className="text-sm font-medium text-foreground">Custom Coordination Parameters</span>
                  </div>

                  <FormField
                    control={form.control}
                    name="paramCoordinationNumber"
                    render={({ field }) => (
                      <FormItem>
                        <FormLabel className="text-xs">
                          Coordination Number (CN): <span className="font-bold">{field.value}</span>
                        </FormLabel>
                        <FormControl>
                          <Slider
                            min={6}
                            max={10}
                            step={1}
                            value={[field.value]}
                            onValueChange={(v) => field.onChange(v[0])}
                          />
                        </FormControl>
                        <div className="flex justify-between text-xs text-muted-foreground">
                          <span>CN=6 (Octahedral)</span>
                          <span>CN=8 (SAP)</span>
                          <span>CN=10 (Rare)</span>
                        </div>
                      </FormItem>
                    )}
                  />

                  <FormField
                    control={form.control}
                    name="paramNumWaters"
                    render={({ field }) => (
                      <FormItem>
                        <FormLabel className="text-xs">
                          Water Molecules: <span className="font-bold">{field.value}</span>
                          <span className="text-muted-foreground ml-2">(Protein donors: {watchParamCoordinationNumber - field.value})</span>
                        </FormLabel>
                        <FormControl>
                          <Slider
                            min={0}
                            max={4}
                            step={1}
                            value={[field.value]}
                            onValueChange={(v) => field.onChange(v[0])}
                          />
                        </FormControl>
                      </FormItem>
                    )}
                  />

                  <FormField
                    control={form.control}
                    name="paramBidentateFraction"
                    render={({ field }) => (
                      <FormItem>
                        <FormLabel className="text-xs">
                          Bidentate Fraction: <span className="font-bold">{Math.round(field.value * 100)}%</span>
                        </FormLabel>
                        <FormControl>
                          <Slider
                            min={0}
                            max={1}
                            step={0.25}
                            value={[field.value]}
                            onValueChange={(v) => field.onChange(v[0])}
                          />
                        </FormControl>
                        <div className="flex justify-between text-xs text-muted-foreground">
                          <span>All mono</span>
                          <span>Mixed</span>
                          <span>All bi</span>
                        </div>
                      </FormItem>
                    )}
                  />
                </div>
              )}

              {/* Legacy Mode */}
              {watchUseLegacyMode && (
                <div className="space-y-4">
                  <div className="p-2 rounded-lg bg-amber-50 dark:bg-amber-950/20 border border-amber-200 dark:border-amber-800">
                    <p className="text-xs text-amber-700 dark:text-amber-300">
                      <span className="font-medium">Legacy Mode:</span> Original 8-residue templates.
                      For best results, use the new Template Library instead.
                    </p>
                  </div>

                  <FormField
                    control={form.control}
                    name="templateType"
                    render={({ field }) => (
                      <FormItem>
                        <FormLabel className="text-xs">Coordination Template (Legacy)</FormLabel>
                        <Select value={field.value} onValueChange={field.onChange}>
                          <FormControl>
                            <SelectTrigger>
                              <SelectValue />
                            </SelectTrigger>
                          </FormControl>
                          <SelectContent>
                            <SelectItem value="ef_hand">EF-Hand (8 residues)</SelectItem>
                            <SelectItem value="c4_symmetric">C4-Symmetric (4 + 4 waters)</SelectItem>
                            <SelectItem value="none">No Template (De novo)</SelectItem>
                          </SelectContent>
                        </Select>
                      </FormItem>
                    )}
                  />

                  <FormField
                    control={form.control}
                    name="donorResidue"
                    render={({ field }) => (
                      <FormItem>
                        <FormLabel className="text-xs">Coordinating Residue Type</FormLabel>
                        <Select value={field.value} onValueChange={field.onChange}>
                          <FormControl>
                            <SelectTrigger>
                              <SelectValue />
                            </SelectTrigger>
                          </FormControl>
                          <SelectContent>
                            <SelectItem value="ASP">ASP (shorter)</SelectItem>
                            <SelectItem value="GLU">GLU (longer)</SelectItem>
                            <SelectItem value="MIXED">MIXED (5 Asp + 3 Glu)</SelectItem>
                          </SelectContent>
                        </Select>
                      </FormItem>
                    )}
                  />

                  <FormField
                    control={form.control}
                    name="useMotifScaffolding"
                    render={({ field }) => (
                      <FormItem className="flex flex-row items-center justify-between rounded-lg border p-3">
                        <div className="space-y-0.5">
                          <FormLabel className="text-sm">Use Motif Scaffolding</FormLabel>
                          <FormDescription className="text-xs">Preserves template residues</FormDescription>
                        </div>
                        <FormControl>
                          <Switch checked={field.value} onCheckedChange={field.onChange} />
                        </FormControl>
                      </FormItem>
                    )}
                  />
                </div>
              )}

              {/* TEBL Options */}
              {isLuminescent && (
                <div className="p-4 rounded-lg bg-muted border border-border">
                  <div className="flex items-center gap-2 mb-3">
                    <Sparkles className="w-4 h-4 text-primary" />
                    <span className="text-sm font-medium text-foreground">TEBL Luminescence Assay</span>
                  </div>
                  <p className="text-xs text-muted-foreground mb-3">
                    Tryptophan-Enhanced Terbium Luminescence for binding validation.
                    Based on Caldwell et al. 2020 (subfemtomolar affinity achieved).
                  </p>
                  <FormField
                    control={form.control}
                    name="addTrpAntenna"
                    render={({ field }) => (
                      <FormItem className="flex flex-row items-start space-x-3 space-y-0">
                        <FormControl>
                          <Checkbox checked={field.value} onCheckedChange={field.onChange} />
                        </FormControl>
                        <div className="space-y-1 leading-none">
                          <FormLabel className="text-sm text-foreground">Add Trp antenna residue</FormLabel>
                          <FormDescription className="text-xs text-muted-foreground">
                            Position Trp within 4-6 A for energy transfer
                          </FormDescription>
                        </div>
                      </FormItem>
                    )}
                  />
                </div>
              )}

              {/* Validation toggle */}
              <FormField
                control={form.control}
                name="validateCoordination"
                render={({ field }) => (
                  <FormItem className="flex flex-row items-start space-x-3 space-y-0">
                    <FormControl>
                      <Checkbox checked={field.value} onCheckedChange={field.onChange} />
                    </FormControl>
                    <div className="space-y-1 leading-none">
                      <FormLabel className="text-sm">Run coordination validation</FormLabel>
                      <FormDescription className="text-xs">
                        Check geometry, bond distances, and donor types
                      </FormDescription>
                    </div>
                  </FormItem>
                )}
              />

              {/* Expected results info */}
              <div className="p-3 rounded-lg bg-muted/50 border">
                <div className="flex items-start gap-2">
                  <Info className="w-4 h-4 text-muted-foreground mt-0.5" />
                  <div className="text-xs text-muted-foreground">
                    <p className="font-medium mb-1 text-foreground">Expected lanthanide coordination:</p>
                    <ul className="list-disc list-inside space-y-0.5">
                      <li>8-9 coordinating atoms (carboxylate oxygens preferred)</li>
                      <li>Metal-O distances: 2.3-2.5 A</li>
                      <li>Both chains should contribute ≥3 donors each</li>
                      {form.watch('addTrpAntenna') && <li>Trp-metal distance: 4-5.5 A for TEBL</li>}
                    </ul>
                  </div>
                </div>
              </div>
            </CardContent>
          </Card>
        )}

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
              name="includeWaters"
              render={({ field }) => (
                <FormItem className="flex flex-row items-start space-x-3 space-y-0">
                  <FormControl>
                    <Checkbox checked={field.value} onCheckedChange={field.onChange} />
                  </FormControl>
                  <div className="space-y-1 leading-none">
                    <FormLabel>Include water molecules in coordination sphere</FormLabel>
                  </div>
                </FormItem>
              )}
            />

            {watchApproach === 'bridging_metal' && (
              <FormField
                control={form.control}
                name="secondMetal"
                render={({ field }) => (
                  <FormItem>
                    <FormLabel>Second Metal (for bridging)</FormLabel>
                    <Select value={field.value || ''} onValueChange={field.onChange}>
                      <FormControl>
                        <SelectTrigger>
                          <SelectValue placeholder={`Same as primary (${watchMetal})`} />
                        </SelectTrigger>
                      </FormControl>
                      <SelectContent>
                        <SelectItem value="">Same as primary ({watchMetal})</SelectItem>
                        {Object.entries(METAL_PROFILES).map(([id, profile]) => (
                          <SelectItem key={id} value={id}>{profile.name}</SelectItem>
                        ))}
                      </SelectContent>
                    </Select>
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
                  <Diamond className="w-4 h-4" />
                </span>
                Designing Metal Dimer...
              </>
            ) : (
              <>
                <Diamond className="w-4 h-4 mr-2" />
                Design {watchMetal} Metal Dimer
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
                <p className="text-xs text-muted-foreground font-semibold uppercase tracking-wide">Metal Dimer Success Criteria</p>
              </div>

              <div className="space-y-2">
                {[
                  `Coordination number: ${metalProfile.coordination.join(' or ')}`,
                  `Bond distance: ${metalProfile.distance[0]}-${metalProfile.distance[1]} A`,
                  'Geometry RMSD: < 0.5 A from ideal',
                  'Both chains contributing donors',
                ].map((criterion, i) => (
                  <div key={i} className="flex items-center gap-2 text-xs">
                    <span className="w-2 h-2 rounded-full bg-emerald-500" />
                    <span className="text-muted-foreground">{criterion}</span>
                  </div>
                ))}
                <div className="flex items-center gap-2 text-xs">
                  <span className="w-2 h-2 rounded-full bg-blue-500" />
                  <span className="text-muted-foreground">Sequence identity: &lt; 70% (for heterodimers)</span>
                </div>
              </div>

              <p className="text-xs text-muted-foreground mt-3 text-center">
                Pipeline: Backbone → LigandMPNN (metal-aware) → ESMFold Validation → Metal Site Analysis
              </p>
            </CardContent>
          </Card>
        </div>
      </form>
    </Form>
  );
}
