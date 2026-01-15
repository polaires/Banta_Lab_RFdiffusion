'use client';

import { useState, useMemo } from 'react';
import { FormSection, FormField } from './shared/FormSection';
import { LengthRangeInput } from './shared/LengthRangeInput';
import { QualityPresetSelector, QualityPreset, QualityParams } from './shared/QualityPresetSelector';
import { AdvancedOptionsWrapper } from './shared/AdvancedOptionsWrapper';
import { QUALITY_PRESETS, RFD3Request, TaskFormProps } from './shared/types';

// Metal ion profiles with coordination chemistry
const METAL_PROFILES = {
  // Transition metals
  ZN: {
    name: 'Zinc (Zn²⁺)',
    category: 'Transition',
    coordination: [4, 5, 6],
    preferredCoord: 4,
    geometry: 'tetrahedral',
    donors: ['His-NE2', 'Cys-SG', 'Asp-OD', 'Glu-OE'],
    donorLabels: ['His', 'Cys', 'Asp', 'Glu'],
    distance: [2.0, 2.3],
    description: 'Most common catalytic metal. Tetrahedral preferred.',
  },
  FE: {
    name: 'Iron (Fe²⁺)',
    category: 'Transition',
    coordination: [4, 5, 6],
    preferredCoord: 6,
    geometry: 'octahedral',
    donors: ['His-NE2', 'Cys-SG', 'Asp-OD', 'Glu-OE', 'Tyr-OH'],
    donorLabels: ['His', 'Cys', 'Asp', 'Glu', 'Tyr'],
    distance: [1.9, 2.2],
    description: 'Redox-active. Found in heme and non-heme sites.',
  },
  CU: {
    name: 'Copper (Cu²⁺)',
    category: 'Transition',
    coordination: [4, 5, 6],
    preferredCoord: 4,
    geometry: 'square_planar',
    donors: ['His-NE2', 'Cys-SG', 'Met-SD'],
    donorLabels: ['His', 'Cys', 'Met'],
    distance: [1.9, 2.1],
    description: 'Electron transfer. Square planar or distorted.',
  },
  MN: {
    name: 'Manganese (Mn²⁺)',
    category: 'Transition',
    coordination: [6],
    preferredCoord: 6,
    geometry: 'octahedral',
    donors: ['His-NE2', 'Asp-OD', 'Glu-OE', 'H2O'],
    donorLabels: ['His', 'Asp', 'Glu', 'Water'],
    distance: [2.1, 2.3],
    description: 'Catalytic. Often in oxygen-evolving complexes.',
  },
  // Alkaline earth
  CA: {
    name: 'Calcium (Ca²⁺)',
    category: 'Alkaline Earth',
    coordination: [6, 7, 8],
    preferredCoord: 7,
    geometry: 'pentagonal_bipyramidal',
    donors: ['Asp-OD', 'Glu-OE', 'Asn-OD', 'backbone-O', 'H2O'],
    donorLabels: ['Asp', 'Glu', 'Asn', 'Backbone O', 'Water'],
    distance: [2.3, 2.5],
    description: 'Structural/signaling. Flexible coordination.',
  },
  MG: {
    name: 'Magnesium (Mg²⁺)',
    category: 'Alkaline Earth',
    coordination: [6],
    preferredCoord: 6,
    geometry: 'octahedral',
    donors: ['Asp-OD', 'Glu-OE', 'H2O'],
    donorLabels: ['Asp', 'Glu', 'Water'],
    distance: [2.0, 2.2],
    description: 'Catalytic. Strict octahedral geometry.',
  },
  // Lanthanides
  TB: {
    name: 'Terbium (Tb³⁺)',
    category: 'Lanthanide',
    coordination: [8, 9],
    preferredCoord: 9,
    geometry: 'tricapped_trigonal_prismatic',
    donors: ['Asp-OD', 'Glu-OE', 'Asn-OD', 'backbone-O', 'H2O'],
    donorLabels: ['Asp', 'Glu', 'Asn', 'Backbone O', 'Water'],
    distance: [2.3, 2.5],
    description: 'Luminescent (green). Antenna effect with Trp.',
    special: 'luminescent',
  },
  GD: {
    name: 'Gadolinium (Gd³⁺)',
    category: 'Lanthanide',
    coordination: [8, 9],
    preferredCoord: 9,
    geometry: 'tricapped_trigonal_prismatic',
    donors: ['Asp-OD', 'Glu-OE', 'Asn-OD', 'backbone-O', 'H2O'],
    donorLabels: ['Asp', 'Glu', 'Asn', 'Backbone O', 'Water'],
    distance: [2.3, 2.5],
    description: 'Paramagnetic. Used in MRI contrast agents.',
    special: 'paramagnetic',
  },
  EU: {
    name: 'Europium (Eu³⁺)',
    category: 'Lanthanide',
    coordination: [8, 9],
    preferredCoord: 9,
    geometry: 'tricapped_trigonal_prismatic',
    donors: ['Asp-OD', 'Glu-OE', 'Asn-OD', 'backbone-O', 'H2O'],
    donorLabels: ['Asp', 'Glu', 'Asn', 'Backbone O', 'Water'],
    distance: [2.3, 2.5],
    description: 'Luminescent (red). Long-lived emission.',
    special: 'luminescent',
  },
} as const;

type MetalId = keyof typeof METAL_PROFILES;

// NEW: Template Library Options (chemically realistic coordination)
const TEMPLATE_LIBRARY_OPTIONS = [
  {
    id: 'caldwell_4',
    name: 'Caldwell (4 Glu)',
    description: '4 bidentate Glu, CN=8 - highest affinity (subfemtomolar)',
    coordinationNumber: 8,
    geometry: 'square_antiprism',
    recommendedFor: ['TB', 'EU'],
    icon: 'star',
    recommended: true,
    deprecated: false,
  },
  {
    id: 'ef_hand_8',
    name: 'EF-Hand Mixed',
    description: '4 mono Asp + 2 bi Glu, CN=8 - balanced design',
    coordinationNumber: 8,
    geometry: 'distorted_square_antiprism',
    recommendedFor: ['GD', 'YB'],
    icon: 'hand_gesture',
    deprecated: false,
  },
  {
    id: 'lanm_mixed',
    name: 'Lanmodulin',
    description: '3 bi Asp + 2 waters, CN=9 - natural-like',
    coordinationNumber: 9,
    geometry: 'tricapped_trigonal_prism',
    recommendedFor: ['TB', 'GD'],
    icon: 'water_drop',
    deprecated: false,
  },
  {
    id: 'high_coord_9',
    name: 'High Coordination',
    description: '4 bi Glu + 1 water, CN=9 - for large lanthanides',
    coordinationNumber: 9,
    geometry: 'tricapped_trigonal_prism',
    recommendedFor: ['LA', 'CE'],
    icon: 'expand',
    deprecated: false,
  },
  {
    id: 'legacy',
    name: 'Legacy (Custom)',
    description: 'Original 8-residue template with custom donor selection',
    coordinationNumber: 16,
    geometry: 'various',
    recommendedFor: [],
    icon: 'tune',
    deprecated: true,
  },
] as const;

type TemplateLibraryId = typeof TEMPLATE_LIBRARY_OPTIONS[number]['id'];

// Metal approach types
type MetalApproach = 'joint_metal' | 'asymmetric_metal' | 'induced_metal' | 'bridging_metal' | 'redox_switch' | 'lanthanide_sensor';

// Design approaches
const APPROACHES: Array<{
  id: MetalApproach;
  name: string;
  description: string;
  details: string;
  icon: string;
  recommended?: boolean;
  advanced?: boolean;
}> = [
  {
    id: 'joint_metal',
    name: 'Joint Metal Dimer',
    description: 'Both chains co-evolve around metal',
    details: 'Simultaneous design of two chains sharing coordination. Each chain contributes donors. Best for symmetric heterodimers.',
    icon: 'hub',
    recommended: true,
  },
  {
    id: 'asymmetric_metal',
    name: 'Asymmetric Metal Dimer',
    description: 'Different donor types per chain',
    details: 'Chain A provides some donor types (e.g., His), Chain B provides others (e.g., Cys). Creates distinct coordination chemistry.',
    icon: 'tune',
  },
  {
    id: 'induced_metal',
    name: 'Metal-Induced Dimerization',
    description: 'Metal binding drives assembly',
    details: 'Monomers with incomplete coordination only become stable when dimerized. Mimics natural metalloenzyme assembly.',
    icon: 'account_tree',
  },
  {
    id: 'bridging_metal',
    name: 'Bridging Metals',
    description: 'Multiple metals at interface',
    details: 'Two or more metal ions bridging chains. Like dinuclear enzyme active sites (Fe-Fe, Cu-Cu).',
    icon: 'link',
    advanced: true,
  },
  {
    id: 'redox_switch',
    name: 'Redox-Switchable',
    description: 'Oxidation state controls dimer',
    details: 'Design dimers with different affinity based on metal oxidation state. Electrochemical control of assembly.',
    icon: 'bolt',
    advanced: true,
  },
];

// Geometry options based on coordination number
const GEOMETRY_OPTIONS: Record<number, string[]> = {
  4: ['tetrahedral', 'square_planar'],
  5: ['trigonal_bipyramidal', 'square_pyramidal'],
  6: ['octahedral'],
  7: ['pentagonal_bipyramidal', 'capped_octahedral'],
  8: ['square_antiprismatic', 'dodecahedral'],
  9: ['tricapped_trigonal_prismatic'],
};

export function InterfaceMetalForm({ onSubmit, isSubmitting, health }: TaskFormProps) {
  // Metal selection
  const [selectedMetal, setSelectedMetal] = useState<MetalId>('ZN');
  const [approach, setApproach] = useState<MetalApproach>('joint_metal');

  // Coordination split
  const [chainADonors, setChainADonors] = useState(2);
  const [chainBDonors, setChainBDonors] = useState(2);

  // Donor type preferences per chain
  const [chainADonorTypes] = useState<string[]>(['His', 'His']);
  const [chainBDonorTypes] = useState<string[]>(['Cys', 'Cys']);

  // Design parameters
  const [chainLength, setChainLength] = useState('60-80');
  const [numDesigns, setNumDesigns] = useState(5);
  const [seed, setSeed] = useState<string>('');
  const [qualityPreset, setQualityPreset] = useState<QualityPreset>('Balanced');
  const [qualityParams, setQualityParams] = useState<QualityParams>(QUALITY_PRESETS.Balanced);

  // Advanced options
  const [targetGeometry, setTargetGeometry] = useState<string>('auto');
  const [includeWaters, setIncludeWaters] = useState(false);
  const [secondMetal, setSecondMetal] = useState<MetalId | ''>('');

  // Lanthanide-specific options (based on Caldwell et al. 2020)
  // NEW: Template library system with chemically realistic coordination
  const [templateName, setTemplateName] = useState<TemplateLibraryId>('caldwell_4');
  const [useLegacyMode, setUseLegacyMode] = useState(false);
  // Legacy options (shown only when legacy mode is enabled)
  const [templateType, setTemplateType] = useState<'ef_hand' | 'c4_symmetric' | 'none'>('ef_hand');
  const [donorResidue, setDonorResidue] = useState<'ASP' | 'GLU' | 'MIXED'>('ASP');
  const [useMotifScaffolding, setUseMotifScaffolding] = useState(true);
  const [addTrpAntenna, setAddTrpAntenna] = useState(false);
  // NEW: Parametric mode options
  const [useParametricMode, setUseParametricMode] = useState(false);
  const [paramCoordinationNumber, setParamCoordinationNumber] = useState(8);
  const [paramNumWaters, setParamNumWaters] = useState(0);
  const [paramBidentateFraction, setParamBidentateFraction] = useState(0.5);
  const [validateCoordination, setValidateCoordination] = useState(true);

  const metalProfile = METAL_PROFILES[selectedMetal];
  const totalCoordination = chainADonors + chainBDonors;
  const geometryOptions = GEOMETRY_OPTIONS[totalCoordination] || [];

  // Check if current metal is a lanthanide
  const isLanthanide = useMemo(() => {
    return metalProfile.category === 'Lanthanide';
  }, [metalProfile]);

  // Update coordination when metal changes
  const handleMetalChange = (metalId: MetalId) => {
    setSelectedMetal(metalId);
    const profile = METAL_PROFILES[metalId];
    const coord = profile.preferredCoord;
    const split = Math.floor(coord / 2);
    setChainADonors(split);
    setChainBDonors(coord - split);
    setTargetGeometry('auto');
  };

  const handleQualityChange = (preset: QualityPreset, params: QualityParams) => {
    setQualityPreset(preset);
    setQualityParams(params);
  };

  const handleSubmit = async () => {
    // Build the request
    const request: RFD3Request = {
      task: 'interface_metal_design',
      approach,
      ligand: selectedMetal, // Metal is specified as ligand
      metal: selectedMetal,  // Also pass as metal parameter
      chain_length: chainLength,
      num_designs: numDesigns,
      num_timesteps: qualityParams.num_timesteps,
      step_scale: qualityParams.step_scale,
      gamma_0: qualityParams.gamma_0,
    };

    // Add metal-specific parameters
    // @ts-expect-error - extended fields for metal design - extended fields for metal design
    request.metal_config = {
      coordination_split: [chainADonors, chainBDonors],
      geometry: targetGeometry === 'auto' ? metalProfile.geometry : targetGeometry,
      chain_a_donors: chainADonorTypes,
      chain_b_donors: chainBDonorTypes,
      include_waters: includeWaters,
    };

    // Add lanthanide-specific parameters
    if (isLanthanide) {
      // NEW: Template library system (chemically realistic coordination)
      if (useParametricMode) {
        // Parametric mode - custom coordination
        // @ts-expect-error - extended fields for metal design
        request.parametric = {
          coordination_number: paramCoordinationNumber,
          num_waters: paramNumWaters,
          bidentate_fraction: paramBidentateFraction,
        };
      } else if (useLegacyMode || templateName === 'legacy') {
        // Legacy mode - backward compatibility
        // @ts-expect-error - extended fields for metal design
        request.template_type = templateType;
        // @ts-expect-error - extended fields for metal design
        request.donor_residue = donorResidue;
        // @ts-expect-error - extended fields for metal design
        request.use_motif_scaffolding = useMotifScaffolding;
      } else {
        // NEW: Template library (recommended)
        // @ts-expect-error - extended fields for metal design
        request.template_name = templateName;
      }
      // @ts-expect-error - extended fields for metal design
      request.add_trp_antenna = addTrpAntenna;
      // @ts-expect-error - extended fields for metal design
      request.validate_coordination = validateCoordination;
      // @ts-expect-error - extended fields for metal design
      request.include_waters = includeWaters;
    }

    if (approach === 'bridging_metal' && secondMetal) {
      // @ts-expect-error - extended fields for metal design
      request.metal_config.second_metal = secondMetal;
    }

    if (seed) {
      request.seed = parseInt(seed, 10);
    }

    await onSubmit(request);
  };

  const isValid = chainLength.trim() !== '' && totalCoordination >= 4;

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex items-center gap-3 pb-4 border-b border-slate-200">
        <div className="w-10 h-10 rounded-lg bg-gradient-to-br from-amber-500 to-orange-500 flex items-center justify-center">
          <span className="material-symbols-outlined text-white">diamond</span>
        </div>
        <div>
          <h2 className="font-semibold text-slate-900">Interface Metal Dimer Design</h2>
          <p className="text-sm text-slate-500">Design protein dimers with metal coordination at the interface</p>
        </div>
      </div>

      {/* Info Banner */}
      <div className="p-4 rounded-xl bg-gradient-to-r from-amber-50 to-orange-50 border border-amber-100">
        <div className="flex items-start gap-3">
          <span className="material-symbols-outlined text-amber-600 mt-0.5">science</span>
          <div className="text-sm text-amber-900">
            <p className="font-medium mb-1">Metal-Mediated Dimerization</p>
            <p className="text-amber-700">
              Unlike organic ligands, metals have fixed coordination geometries. Each chain contributes
              donor atoms (His, Cys, Asp, Glu) that together complete the metal&apos;s coordination sphere.
            </p>
          </div>
        </div>
      </div>

      {/* Metal Selection */}
      <FormSection title="Metal Ion" description="Select the metal for interface coordination" required>
        <div className="space-y-4">
          {/* Category tabs */}
          <div className="flex gap-2 flex-wrap">
            {['Transition', 'Alkaline Earth', 'Lanthanide'].map((category) => (
              <div key={category} className="space-y-2">
                <span className="text-xs font-medium text-slate-500 uppercase tracking-wide">{category}</span>
                <div className="flex gap-2 flex-wrap">
                  {Object.entries(METAL_PROFILES)
                    .filter(([, profile]) => profile.category === category)
                    .map(([id, profile]) => (
                      <button
                        key={id}
                        onClick={() => handleMetalChange(id as MetalId)}
                        className={`px-4 py-2 rounded-lg border text-sm font-medium transition-all ${
                          selectedMetal === id
                            ? 'border-amber-400 bg-amber-50 text-amber-700 ring-1 ring-amber-200'
                            : 'border-slate-200 text-slate-600 hover:border-slate-300 hover:bg-slate-50'
                        }`}
                      >
                        <span className="font-bold">{id}</span>
                        <span className="text-slate-400 ml-1">
                          {'special' in profile && profile.special === 'luminescent' && '✨'}
                          {'special' in profile && profile.special === 'paramagnetic' && '🧲'}
                        </span>
                      </button>
                    ))}
                </div>
              </div>
            ))}
          </div>

          {/* Selected metal info */}
          <div className="p-3 rounded-lg bg-slate-50 border border-slate-200">
            <div className="flex items-center gap-2 mb-2">
              <span className="text-lg font-bold text-amber-600">{selectedMetal}</span>
              <span className="text-sm text-slate-600">{metalProfile.name}</span>
            </div>
            <p className="text-xs text-slate-500 mb-2">{metalProfile.description}</p>
            <div className="grid grid-cols-3 gap-2 text-xs">
              <div>
                <span className="text-slate-400">Coordination: </span>
                <span className="font-medium">{metalProfile.coordination.join('-')}</span>
              </div>
              <div>
                <span className="text-slate-400">Geometry: </span>
                <span className="font-medium capitalize">{metalProfile.geometry.replace('_', ' ')}</span>
              </div>
              <div>
                <span className="text-slate-400">Distance: </span>
                <span className="font-medium">{metalProfile.distance[0]}-{metalProfile.distance[1]} Å</span>
              </div>
            </div>
            <div className="mt-2 flex flex-wrap gap-1">
              <span className="text-xs text-slate-400">Donors: </span>
              {metalProfile.donorLabels.map((donor) => (
                <span key={donor} className="px-1.5 py-0.5 bg-slate-200 text-slate-600 rounded text-xs">
                  {donor}
                </span>
              ))}
            </div>
          </div>
        </div>
      </FormSection>

      {/* Design Approach */}
      <FormSection title="Design Approach" description="How to coordinate the metal between chains" required>
        <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-3">
          {APPROACHES.filter(a => !a.advanced || approach === a.id).map((app) => (
            <button
              key={app.id}
              onClick={() => setApproach(app.id)}
              className={`p-4 rounded-xl border text-left transition-all relative ${
                approach === app.id
                  ? 'border-amber-400 bg-amber-50 ring-1 ring-amber-200'
                  : 'border-slate-200 hover:border-slate-300 hover:bg-slate-50'
              }`}
            >
              {app.recommended && (
                <span className="absolute -top-2 -right-2 px-2 py-0.5 bg-emerald-500 text-white text-xs font-medium rounded-full">
                  Recommended
                </span>
              )}
              <div className="flex items-center gap-2 mb-2">
                <span className={`material-symbols-outlined text-lg ${
                  approach === app.id ? 'text-amber-600' : 'text-slate-400'
                }`}>
                  {app.icon}
                </span>
                <span className="font-medium text-slate-900">{app.name}</span>
              </div>
              <p className="text-sm text-slate-600">{app.description}</p>
              <p className="text-xs text-slate-400 mt-1">{app.details}</p>
            </button>
          ))}
        </div>

        {/* Show advanced approaches toggle */}
        <button
          onClick={() => setApproach(approach === 'bridging_metal' ? 'joint_metal' : 'bridging_metal')}
          className="text-xs text-amber-600 hover:underline mt-2"
        >
          {APPROACHES.some(a => a.advanced && approach !== a.id) ? 'Show advanced approaches' : 'Hide advanced'}
        </button>
      </FormSection>

      {/* Coordination Split */}
      <FormSection title="Coordination Split" description="How many donor atoms from each chain">
        <div className="space-y-4">
          <div className="flex items-center gap-6">
            <div className="flex-1 text-center">
              <label className="text-xs text-slate-500 mb-1 block">Chain A Donors</label>
              <div className="flex items-center justify-center gap-2">
                <button
                  onClick={() => setChainADonors(Math.max(1, chainADonors - 1))}
                  className="w-8 h-8 rounded-lg border border-slate-200 hover:bg-slate-100"
                >
                  -
                </button>
                <span className="w-12 h-12 rounded-xl bg-purple-100 text-purple-700 flex items-center justify-center text-2xl font-bold">
                  {chainADonors}
                </span>
                <button
                  onClick={() => setChainADonors(chainADonors + 1)}
                  className="w-8 h-8 rounded-lg border border-slate-200 hover:bg-slate-100"
                >
                  +
                </button>
              </div>
            </div>

            <div className="flex flex-col items-center">
              <span className="text-2xl text-slate-300">+</span>
              <span className="material-symbols-outlined text-amber-500">diamond</span>
            </div>

            <div className="flex-1 text-center">
              <label className="text-xs text-slate-500 mb-1 block">Chain B Donors</label>
              <div className="flex items-center justify-center gap-2">
                <button
                  onClick={() => setChainBDonors(Math.max(1, chainBDonors - 1))}
                  className="w-8 h-8 rounded-lg border border-slate-200 hover:bg-slate-100"
                >
                  -
                </button>
                <span className="w-12 h-12 rounded-xl bg-cyan-100 text-cyan-700 flex items-center justify-center text-2xl font-bold">
                  {chainBDonors}
                </span>
                <button
                  onClick={() => setChainBDonors(chainBDonors + 1)}
                  className="w-8 h-8 rounded-lg border border-slate-200 hover:bg-slate-100"
                >
                  +
                </button>
              </div>
            </div>

            <div className="flex flex-col items-center">
              <span className="text-2xl text-slate-300">=</span>
            </div>

            <div className="text-center">
              <label className="text-xs text-slate-500 mb-1 block">Total</label>
              <span className={`w-12 h-12 rounded-xl flex items-center justify-center text-2xl font-bold ${
                (metalProfile.coordination as readonly number[]).includes(totalCoordination)
                  ? 'bg-emerald-100 text-emerald-700'
                  : 'bg-red-100 text-red-700'
              }`}>
                {totalCoordination}
              </span>
            </div>
          </div>

          {/* Validation message */}
          {!(metalProfile.coordination as readonly number[]).includes(totalCoordination) && (
            <div className="p-2 rounded-lg bg-amber-50 border border-amber-200 text-xs text-amber-700">
              <span className="material-symbols-outlined text-sm align-middle mr-1">warning</span>
              {selectedMetal} typically has {metalProfile.coordination.join(' or ')}-coordinate geometry.
              Total {totalCoordination} may not be optimal.
            </div>
          )}

          {/* Geometry selector */}
          {geometryOptions.length > 0 && (
            <div>
              <label className="text-xs text-slate-500 block mb-1">Target Geometry</label>
              <div className="flex gap-2 flex-wrap">
                <button
                  onClick={() => setTargetGeometry('auto')}
                  className={`px-3 py-1.5 rounded-lg text-xs font-medium transition-all ${
                    targetGeometry === 'auto'
                      ? 'bg-amber-100 text-amber-700'
                      : 'bg-slate-100 text-slate-600 hover:bg-slate-200'
                  }`}
                >
                  Auto ({metalProfile.geometry.replace('_', ' ')})
                </button>
                {geometryOptions.map((geom) => (
                  <button
                    key={geom}
                    onClick={() => setTargetGeometry(geom)}
                    className={`px-3 py-1.5 rounded-lg text-xs font-medium capitalize transition-all ${
                      targetGeometry === geom
                        ? 'bg-amber-100 text-amber-700'
                        : 'bg-slate-100 text-slate-600 hover:bg-slate-200'
                    }`}
                  >
                    {geom.replace('_', ' ')}
                  </button>
                ))}
              </div>
            </div>
          )}
        </div>
      </FormSection>

      {/* Design Parameters */}
      <FormSection title="Design Parameters" required>
        <div className="space-y-4">
          <div className="flex gap-4 items-start">
            <div className="flex-1">
              <LengthRangeInput
                value={chainLength}
                onChange={setChainLength}
                label="Chain Length"
                placeholder="60-80"
                hint="Per chain length range"
              />
            </div>
            <FormField label="# Designs" className="w-24">
              <input
                type="number"
                value={numDesigns}
                onChange={(e) => setNumDesigns(Math.max(1, parseInt(e.target.value) || 1))}
                min={1}
                max={20}
                className="w-full px-3 py-2 rounded-lg border border-slate-200 focus:border-amber-400 focus:ring-2 focus:ring-amber-100 outline-none text-sm"
              />
            </FormField>
            <FormField label="Seed" hint="Optional" className="w-24">
              <input
                type="number"
                value={seed}
                onChange={(e) => setSeed(e.target.value)}
                placeholder="Random"
                className="w-full px-3 py-2 rounded-lg border border-slate-200 focus:border-amber-400 focus:ring-2 focus:ring-amber-100 outline-none text-sm"
              />
            </FormField>
          </div>

          <QualityPresetSelector
            value={qualityPreset}
            onChange={handleQualityChange}
            showDescription
          />
        </div>
      </FormSection>

      {/* Lanthanide-Specific Options - appears only for lanthanides */}
      {isLanthanide && (
        <FormSection
          title="Lanthanide Design Options"
          description="Chemically realistic templates for optimal coordination (CN=8-9)"
        >
          <div className="space-y-4">
            {/* Mode Selection - Library vs Parametric vs Legacy */}
            <div className="flex gap-2 p-1 bg-slate-100 rounded-lg">
              <button
                onClick={() => { setUseParametricMode(false); setUseLegacyMode(false); }}
                className={`flex-1 px-3 py-2 rounded-md text-sm font-medium transition-all ${
                  !useParametricMode && !useLegacyMode
                    ? 'bg-white text-purple-700 shadow-sm'
                    : 'text-slate-600 hover:text-slate-800'
                }`}
              >
                Template Library
              </button>
              <button
                onClick={() => { setUseParametricMode(true); setUseLegacyMode(false); }}
                className={`flex-1 px-3 py-2 rounded-md text-sm font-medium transition-all ${
                  useParametricMode
                    ? 'bg-white text-purple-700 shadow-sm'
                    : 'text-slate-600 hover:text-slate-800'
                }`}
              >
                Parametric
              </button>
              <button
                onClick={() => { setUseParametricMode(false); setUseLegacyMode(true); }}
                className={`flex-1 px-3 py-2 rounded-md text-sm font-medium transition-all ${
                  useLegacyMode
                    ? 'bg-white text-purple-700 shadow-sm'
                    : 'text-slate-600 hover:text-slate-800'
                }`}
              >
                Legacy
              </button>
            </div>

            {/* NEW: Template Library Selection */}
            {!useParametricMode && !useLegacyMode && (
              <div>
                <label className="text-xs text-slate-500 block mb-2">
                  Select Template
                  <span className="ml-2 px-1.5 py-0.5 bg-emerald-100 text-emerald-700 text-xs rounded">New!</span>
                </label>
                <div className="grid grid-cols-1 sm:grid-cols-2 gap-3">
                  {TEMPLATE_LIBRARY_OPTIONS.filter(t => !t.deprecated).map((template) => {
                    const isRecommended = (template.recommendedFor as readonly string[]).includes(selectedMetal);
                    return (
                      <button
                        key={template.id}
                        onClick={() => setTemplateName(template.id)}
                        className={`p-3 rounded-lg border text-left transition-all ${
                          templateName === template.id
                            ? 'border-purple-400 bg-purple-50 ring-1 ring-purple-200'
                            : 'border-slate-200 hover:border-slate-300 hover:bg-slate-50'
                        }`}
                      >
                        <div className="flex items-center gap-2 mb-1">
                          <span className={`material-symbols-outlined text-sm ${
                            templateName === template.id ? 'text-purple-600' : 'text-slate-400'
                          }`}>{template.icon}</span>
                          <span className="font-medium text-sm">{template.name}</span>
                          {isRecommended && (
                            <span className="px-1.5 py-0.5 bg-emerald-100 text-emerald-700 text-xs rounded">
                              Best for {selectedMetal}
                            </span>
                          )}
                        </div>
                        <p className="text-xs text-slate-500">{template.description}</p>
                        <div className="flex items-center gap-2 mt-1">
                          <span className="text-xs text-purple-600">CN={template.coordinationNumber}</span>
                          <span className="text-xs text-slate-400">•</span>
                          <span className="text-xs text-slate-500">{template.geometry}</span>
                        </div>
                      </button>
                    );
                  })}
                </div>
              </div>
            )}

            {/* NEW: Parametric Mode Options */}
            {useParametricMode && (
              <div className="space-y-4 p-4 rounded-lg bg-gradient-to-r from-blue-50 to-indigo-50 border border-blue-200">
                <div className="flex items-center gap-2 mb-3">
                  <span className="material-symbols-outlined text-blue-600">tune</span>
                  <span className="text-sm font-medium text-blue-800">Custom Coordination Parameters</span>
                </div>

                {/* Coordination Number Slider */}
                <div>
                  <label className="text-xs text-slate-600 block mb-2">
                    Coordination Number (CN): <span className="font-bold">{paramCoordinationNumber}</span>
                  </label>
                  <input
                    type="range"
                    min="6"
                    max="10"
                    value={paramCoordinationNumber}
                    onChange={(e) => setParamCoordinationNumber(parseInt(e.target.value))}
                    className="w-full h-2 bg-blue-200 rounded-lg appearance-none cursor-pointer"
                  />
                  <div className="flex justify-between text-xs text-slate-400 mt-1">
                    <span>CN=6 (Octahedral)</span>
                    <span>CN=8 (SAP)</span>
                    <span>CN=10 (Rare)</span>
                  </div>
                </div>

                {/* Waters Slider */}
                <div>
                  <label className="text-xs text-slate-600 block mb-2">
                    Water Molecules: <span className="font-bold">{paramNumWaters}</span>
                    <span className="text-slate-400 ml-2">(Protein donors: {paramCoordinationNumber - paramNumWaters})</span>
                  </label>
                  <input
                    type="range"
                    min="0"
                    max="4"
                    value={paramNumWaters}
                    onChange={(e) => setParamNumWaters(parseInt(e.target.value))}
                    className="w-full h-2 bg-blue-200 rounded-lg appearance-none cursor-pointer"
                  />
                </div>

                {/* Bidentate Fraction Slider */}
                <div>
                  <label className="text-xs text-slate-600 block mb-2">
                    Bidentate Fraction: <span className="font-bold">{Math.round(paramBidentateFraction * 100)}%</span>
                  </label>
                  <input
                    type="range"
                    min="0"
                    max="1"
                    step="0.25"
                    value={paramBidentateFraction}
                    onChange={(e) => setParamBidentateFraction(parseFloat(e.target.value))}
                    className="w-full h-2 bg-blue-200 rounded-lg appearance-none cursor-pointer"
                  />
                  <div className="flex justify-between text-xs text-slate-400 mt-1">
                    <span>All mono</span>
                    <span>Mixed</span>
                    <span>All bi</span>
                  </div>
                </div>
              </div>
            )}

            {/* LEGACY: Old Template Selection */}
            {useLegacyMode && (
              <>
                <div className="p-2 rounded-lg bg-amber-50 border border-amber-200">
                  <p className="text-xs text-amber-700">
                    <span className="font-medium">Legacy Mode:</span> Original 8-residue templates.
                    For best results, use the new Template Library instead.
                  </p>
                </div>
                <div>
                  <label className="text-xs text-slate-500 block mb-2">Coordination Template (Legacy)</label>
                  <div className="grid grid-cols-1 sm:grid-cols-3 gap-3">
                    <button
                      onClick={() => setTemplateType('ef_hand')}
                      className={`p-3 rounded-lg border text-left transition-all ${
                        templateType === 'ef_hand'
                          ? 'border-purple-400 bg-purple-50 ring-1 ring-purple-200'
                          : 'border-slate-200 hover:border-slate-300 hover:bg-slate-50'
                      }`}
                    >
                      <div className="flex items-center gap-2 mb-1">
                        <span className={`material-symbols-outlined text-sm ${
                          templateType === 'ef_hand' ? 'text-purple-600' : 'text-slate-400'
                        }`}>hand_gesture</span>
                        <span className="font-medium text-sm">EF-Hand</span>
                      </div>
                      <p className="text-xs text-slate-500">8 residues (4 per chain)</p>
                    </button>

                    <button
                      onClick={() => setTemplateType('c4_symmetric')}
                      className={`p-3 rounded-lg border text-left transition-all ${
                        templateType === 'c4_symmetric'
                          ? 'border-purple-400 bg-purple-50 ring-1 ring-purple-200'
                          : 'border-slate-200 hover:border-slate-300 hover:bg-slate-50'
                      }`}
                    >
                      <div className="flex items-center gap-2 mb-1">
                        <span className={`material-symbols-outlined text-sm ${
                          templateType === 'c4_symmetric' ? 'text-purple-600' : 'text-slate-400'
                        }`}>symmetry</span>
                        <span className="font-medium text-sm">C4-Symmetric</span>
                      </div>
                      <p className="text-xs text-slate-500">4 residues + 4 waters</p>
                    </button>

                    <button
                      onClick={() => setTemplateType('none')}
                      className={`p-3 rounded-lg border text-left transition-all ${
                        templateType === 'none'
                          ? 'border-purple-400 bg-purple-50 ring-1 ring-purple-200'
                          : 'border-slate-200 hover:border-slate-300 hover:bg-slate-50'
                      }`}
                    >
                      <div className="flex items-center gap-2 mb-1">
                        <span className={`material-symbols-outlined text-sm ${
                          templateType === 'none' ? 'text-purple-600' : 'text-slate-400'
                        }`}>edit_off</span>
                        <span className="font-medium text-sm">No Template</span>
                      </div>
                      <p className="text-xs text-slate-500">De novo design</p>
                    </button>
                  </div>
                </div>

                {/* Legacy Donor Residue Selection */}
                {templateType !== 'none' && (
                  <div>
                    <label className="text-xs text-slate-500 block mb-2">Coordinating Residue Type</label>
                    <div className="grid grid-cols-1 sm:grid-cols-3 gap-3">
                      <button
                        onClick={() => setDonorResidue('ASP')}
                        className={`p-3 rounded-lg border text-left transition-all ${
                          donorResidue === 'ASP'
                            ? 'border-emerald-400 bg-emerald-50 ring-1 ring-emerald-200'
                            : 'border-slate-200 hover:border-slate-300 hover:bg-slate-50'
                        }`}
                      >
                        <span className="font-bold text-sm">ASP</span>
                        <p className="text-xs text-slate-500">Shorter side chain</p>
                      </button>
                      <button
                        onClick={() => setDonorResidue('GLU')}
                        className={`p-3 rounded-lg border text-left transition-all ${
                          donorResidue === 'GLU'
                            ? 'border-emerald-400 bg-emerald-50 ring-1 ring-emerald-200'
                            : 'border-slate-200 hover:border-slate-300 hover:bg-slate-50'
                        }`}
                      >
                        <span className="font-bold text-sm">GLU</span>
                        <p className="text-xs text-slate-500">Longer side chain</p>
                      </button>
                      <button
                        onClick={() => setDonorResidue('MIXED')}
                        className={`p-3 rounded-lg border text-left transition-all ${
                          donorResidue === 'MIXED'
                            ? 'border-emerald-400 bg-emerald-50 ring-1 ring-emerald-200'
                            : 'border-slate-200 hover:border-slate-300 hover:bg-slate-50'
                        }`}
                      >
                        <span className="font-bold text-sm">MIXED</span>
                        <p className="text-xs text-slate-500">5 Asp + 3 Glu</p>
                      </button>
                    </div>
                  </div>
                )}

                {/* Legacy Motif Scaffolding Toggle */}
                {templateType !== 'none' && (
                  <div className="p-3 rounded-lg bg-gradient-to-r from-emerald-50 to-teal-50 border border-emerald-200">
                    <label className="flex items-center gap-3 cursor-pointer">
                      <input
                        type="checkbox"
                        checked={useMotifScaffolding}
                        onChange={(e) => setUseMotifScaffolding(e.target.checked)}
                        className="w-4 h-4 rounded border-emerald-300 text-emerald-600"
                      />
                      <div>
                        <span className="text-sm text-emerald-800 font-medium">Use Motif Scaffolding</span>
                        <p className="text-xs text-emerald-600">Preserves template residues</p>
                      </div>
                    </label>
                  </div>
                )}
              </>
            )}

            {/* TEBL Options - only for luminescent lanthanides */}
            {'special' in metalProfile && metalProfile.special === 'luminescent' && (
              <div className="p-4 rounded-lg bg-gradient-to-r from-purple-50 to-pink-50 border border-purple-200">
                <div className="flex items-center gap-2 mb-3">
                  <span className="material-symbols-outlined text-purple-600">auto_awesome</span>
                  <span className="text-sm font-medium text-purple-800">TEBL Luminescence Assay</span>
                </div>
                <p className="text-xs text-purple-700 mb-3">
                  Tryptophan-Enhanced Terbium Luminescence for binding validation.
                  Based on Caldwell et al. 2020 (subfemtomolar affinity achieved).
                </p>
                <label className="flex items-center gap-3 cursor-pointer">
                  <input
                    type="checkbox"
                    checked={addTrpAntenna}
                    onChange={(e) => setAddTrpAntenna(e.target.checked)}
                    className="w-4 h-4 rounded border-purple-300 text-purple-600"
                  />
                  <div>
                    <span className="text-sm text-purple-800 font-medium">Add Trp antenna residue</span>
                    <p className="text-xs text-purple-600">Position Trp within 4-6 Å for energy transfer</p>
                  </div>
                </label>
              </div>
            )}

            {/* Validation toggle */}
            <label className="flex items-center gap-3 cursor-pointer">
              <input
                type="checkbox"
                checked={validateCoordination}
                onChange={(e) => setValidateCoordination(e.target.checked)}
                className="w-4 h-4 rounded border-slate-300 text-purple-600"
              />
              <div>
                <span className="text-sm text-slate-700">Run coordination validation</span>
                <p className="text-xs text-slate-500">Check geometry, bond distances, and donor types using coordination.py</p>
              </div>
            </label>

            {/* Info about expected results */}
            <div className="p-3 rounded-lg bg-slate-50 border border-slate-200">
              <div className="flex items-start gap-2">
                <span className="material-symbols-outlined text-slate-500 text-sm mt-0.5">info</span>
                <div className="text-xs text-slate-600">
                  <p className="font-medium mb-1">Expected lanthanide coordination:</p>
                  <ul className="list-disc list-inside space-y-0.5">
                    <li>8-9 coordinating atoms (carboxylate oxygens preferred)</li>
                    <li>Metal-O distances: 2.3-2.5 Å</li>
                    <li>Both chains should contribute ≥3 donors each</li>
                    {addTrpAntenna && <li>Trp-metal distance: 4-5.5 Å for TEBL</li>}
                  </ul>
                </div>
              </div>
            </div>
          </div>
        </FormSection>
      )}

      {/* Advanced Options */}
      <AdvancedOptionsWrapper title="Advanced Options">
        <div className="space-y-4">
          <label className="flex items-center gap-3 cursor-pointer">
            <input
              type="checkbox"
              checked={includeWaters}
              onChange={(e) => setIncludeWaters(e.target.checked)}
              className="w-4 h-4 rounded border-slate-300 text-amber-600"
            />
            <span className="text-sm text-slate-600">Include water molecules in coordination sphere</span>
          </label>

          {approach === 'bridging_metal' && (
            <div>
              <label className="text-xs text-slate-500 block mb-1">Second Metal (for bridging)</label>
              <select
                value={secondMetal}
                onChange={(e) => setSecondMetal(e.target.value as MetalId)}
                className="w-full px-3 py-2 rounded-lg border border-slate-200 text-sm"
              >
                <option value="">Same as primary ({selectedMetal})</option>
                {Object.entries(METAL_PROFILES).map(([id, profile]) => (
                  <option key={id} value={id}>{profile.name}</option>
                ))}
              </select>
            </div>
          )}
        </div>
      </AdvancedOptionsWrapper>

      {/* Submit Button */}
      <div className="pt-4 border-t border-slate-200">
        <button
          onClick={handleSubmit}
          disabled={!isValid || isSubmitting || !health}
          className={`w-full py-3.5 px-6 rounded-xl font-semibold text-white transition-all flex items-center justify-center gap-2 ${
            isValid && !isSubmitting && !!health
              ? 'bg-gradient-to-r from-amber-500 to-orange-500 hover:from-amber-600 hover:to-orange-600 shadow-lg shadow-amber-500/20'
              : 'bg-slate-300 cursor-not-allowed'
          }`}
        >
          {isSubmitting ? (
            <>
              <span className="material-symbols-outlined animate-spin">progress_activity</span>
              Designing Metal Dimer...
            </>
          ) : (
            <>
              <span className="material-symbols-outlined">diamond</span>
              Design {selectedMetal} Metal Dimer
            </>
          )}
        </button>
        {!health && (
          <p className="text-center text-sm text-slate-500 mt-2">
            Backend service unavailable
          </p>
        )}

        {/* Success Criteria */}
        <div className="mt-4 p-4 rounded-xl bg-gradient-to-br from-slate-50 to-slate-100 border border-slate-200">
          <div className="flex items-center gap-2 mb-3">
            <span className="material-symbols-outlined text-slate-500 text-sm">checklist</span>
            <p className="text-xs text-slate-600 font-semibold uppercase tracking-wide">Metal Dimer Success Criteria</p>
          </div>

          <div className="space-y-2">
            <div className="flex items-center gap-2 text-xs">
              <span className="w-2 h-2 rounded-full bg-emerald-500" />
              <span className="text-slate-600">Coordination number: {metalProfile.coordination.join(' or ')}</span>
            </div>
            <div className="flex items-center gap-2 text-xs">
              <span className="w-2 h-2 rounded-full bg-emerald-500" />
              <span className="text-slate-600">Bond distance: {metalProfile.distance[0]}-{metalProfile.distance[1]} Å</span>
            </div>
            <div className="flex items-center gap-2 text-xs">
              <span className="w-2 h-2 rounded-full bg-emerald-500" />
              <span className="text-slate-600">Geometry RMSD: {"<"} 0.5 Å from ideal</span>
            </div>
            <div className="flex items-center gap-2 text-xs">
              <span className="w-2 h-2 rounded-full bg-emerald-500" />
              <span className="text-slate-600">Both chains contributing donors</span>
            </div>
            <div className="flex items-center gap-2 text-xs">
              <span className="w-2 h-2 rounded-full bg-blue-500" />
              <span className="text-slate-600">Sequence identity: {"<"} 70% (for heterodimers)</span>
            </div>
          </div>

          <p className="text-xs text-slate-500 mt-3 text-center">
            Pipeline: Backbone → LigandMPNN (metal-aware) → ESMFold Validation → Metal Site Analysis
          </p>
        </div>
      </div>
    </div>
  );
}
