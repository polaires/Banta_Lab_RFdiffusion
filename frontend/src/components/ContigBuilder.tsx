'use client';

import { useState } from 'react';
import { Wrench, HelpCircle, Zap, GripVertical, Trash2, Plus } from 'lucide-react';

type SegmentType = 'denovo' | 'fixed' | 'linker';

interface Segment {
  id: string;
  type: SegmentType;
  // For de novo: length range
  minLength?: number;
  maxLength?: number;
  // For fixed: chain and residue range
  chain?: string;
  startResidue?: number;
  endResidue?: number;
  // For linker: gap specification
  gapMin?: number;
  gapMax?: number;
}

const SEGMENT_PRESETS = {
  denovo: { label: 'De Novo Region', description: 'Design new residues from scratch', color: 'bg-blue-600', lightBg: 'bg-blue-50', border: 'border-blue-200' },
  fixed: { label: 'Fixed Region', description: 'Keep existing residues from input PDB', color: 'bg-emerald-600', lightBg: 'bg-emerald-50', border: 'border-emerald-200' },
  linker: { label: 'Flexible Linker', description: 'Variable-length connection between segments', color: 'bg-violet-600', lightBg: 'bg-violet-50', border: 'border-violet-200' },
};

// Binder design presets
const BINDER_PRESETS = [
  {
    name: 'Simple Binder',
    description: 'Design a binder to target chain A',
    segments: [
      { id: '1', type: 'fixed' as SegmentType, chain: 'A', startResidue: 1, endResidue: 100 },
      { id: '2', type: 'linker' as SegmentType, gapMin: 0, gapMax: 0 },
      { id: '3', type: 'denovo' as SegmentType, minLength: 60, maxLength: 100 },
    ],
  },
  {
    name: 'Hotspot Binder',
    description: 'Design binder targeting specific residues (hotspots)',
    segments: [
      { id: '1', type: 'fixed' as SegmentType, chain: 'A', startResidue: 25, endResidue: 75 },
      { id: '2', type: 'linker' as SegmentType, gapMin: 0, gapMax: 0 },
      { id: '3', type: 'denovo' as SegmentType, minLength: 80, maxLength: 120 },
    ],
  },
  {
    name: 'Scaffold Inpainting',
    description: 'Design scaffold connecting two fixed motifs',
    segments: [
      { id: '1', type: 'fixed' as SegmentType, chain: 'A', startResidue: 1, endResidue: 20 },
      { id: '2', type: 'linker' as SegmentType, gapMin: 0, gapMax: 0 },
      { id: '3', type: 'denovo' as SegmentType, minLength: 30, maxLength: 50 },
      { id: '4', type: 'linker' as SegmentType, gapMin: 0, gapMax: 0 },
      { id: '5', type: 'fixed' as SegmentType, chain: 'A', startResidue: 80, endResidue: 100 },
    ],
  },
  {
    name: 'Symmetric Scaffold',
    description: 'Design symmetric homo-oligomer scaffold',
    segments: [
      { id: '1', type: 'denovo' as SegmentType, minLength: 80, maxLength: 80 },
    ],
  },
];

interface ContigBuilderProps {
  onContigChange: (contig: string) => void;
  initialContig?: string;
}

export function ContigBuilder({ onContigChange, initialContig }: ContigBuilderProps) {
  const [segments, setSegments] = useState<Segment[]>([
    { id: '1', type: 'denovo', minLength: 100, maxLength: 100 }
  ]);
  const [showHelp, setShowHelp] = useState(false);
  const [showPresets, setShowPresets] = useState(false);

  const generateId = () => Math.random().toString(36).substring(2, 9);

  const applyPreset = (preset: typeof BINDER_PRESETS[0]) => {
    const newSegments = preset.segments.map(seg => ({
      ...seg,
      id: generateId(),
    }));
    setSegments(newSegments);
    setShowPresets(false);
  };

  const addSegment = (type: SegmentType) => {
    const newSegment: Segment = {
      id: generateId(),
      type,
      ...(type === 'denovo' && { minLength: 50, maxLength: 50 }),
      ...(type === 'fixed' && { chain: 'A', startResidue: 1, endResidue: 50 }),
      ...(type === 'linker' && { gapMin: 0, gapMax: 0 }),
    };
    setSegments([...segments, newSegment]);
  };

  const removeSegment = (id: string) => {
    if (segments.length > 1) {
      setSegments(segments.filter(s => s.id !== id));
    }
  };

  const updateSegment = (id: string, updates: Partial<Segment>) => {
    setSegments(segments.map(s => s.id === id ? { ...s, ...updates } : s));
  };

  // Build contig string from segments
  const buildContigString = (): string => {
    return segments.map(seg => {
      switch (seg.type) {
        case 'denovo':
          if (seg.minLength === seg.maxLength) {
            return `${seg.minLength}`;
          }
          return `${seg.minLength}-${seg.maxLength}`;
        case 'fixed':
          return `${seg.chain}${seg.startResidue}-${seg.endResidue}`;
        case 'linker':
          if (seg.gapMin === seg.gapMax) {
            return `/${seg.gapMin}`;
          }
          return `/${seg.gapMin}-${seg.gapMax}`;
        default:
          return '';
      }
    }).filter(Boolean).join(' ');
  };

  const contigString = buildContigString();

  const applyContig = () => {
    onContigChange(contigString);
  };

  return (
    <div className="bg-card rounded-xl p-5 space-y-4 border border-border shadow-sm">
      <div className="flex items-center justify-between">
        <div className="flex items-center gap-2">
          <Wrench className="w-5 h-5 text-blue-600" />
          <h3 className="font-semibold text-foreground">Contig Builder</h3>
        </div>
        <button
          onClick={() => setShowHelp(!showHelp)}
          className="p-1.5 hover:bg-muted rounded-lg transition-colors"
          title="Help"
        >
          <HelpCircle className="w-5 h-5 text-muted-foreground" />
        </button>
      </div>

      {showHelp && (
        <div className="bg-muted p-4 rounded-lg text-sm text-muted-foreground space-y-2 border border-border">
          <p><strong className="text-foreground">De Novo Region:</strong> Generates new amino acids. Specify length or range.</p>
          <p><strong className="text-foreground">Fixed Region:</strong> Keeps residues from input PDB. Specify chain and residue numbers.</p>
          <p><strong className="text-foreground">Flexible Linker:</strong> Gap between segments. Use 0 for direct connection.</p>
          <p className="pt-2 border-t border-border mt-2">
            <strong className="text-foreground">Binder Design:</strong> Upload target PDB, add Fixed region (target), Linker (/0), De Novo region (binder).
          </p>
        </div>
      )}

      {/* Design Presets */}
      <div className="space-y-2">
        <button
          onClick={() => setShowPresets(!showPresets)}
          className="flex items-center gap-2 text-sm text-blue-600 hover:text-blue-700 font-medium transition-colors"
        >
          <Zap className="w-4 h-4" />
          {showPresets ? 'Hide Design Presets' : 'Show Design Presets'}
        </button>

        {showPresets && (
          <div className="grid grid-cols-2 gap-2">
            {BINDER_PRESETS.map((preset) => (
              <button
                key={preset.name}
                onClick={() => applyPreset(preset)}
                className="p-3 text-left bg-muted hover:bg-muted rounded-lg border border-border hover:border-blue-300 transition-all"
              >
                <div className="font-medium text-sm text-foreground">{preset.name}</div>
                <div className="text-xs text-muted-foreground mt-0.5">{preset.description}</div>
              </button>
            ))}
          </div>
        )}
      </div>

      {/* Segment List */}
      <div className="space-y-2">
        {segments.map((segment, index) => {
          const preset = SEGMENT_PRESETS[segment.type];
          return (
            <div
              key={segment.id}
              className={`flex items-center gap-3 p-3 rounded-lg border ${preset.lightBg} ${preset.border}`}
            >
              <GripVertical className="w-4 h-4 text-muted-foreground cursor-grab" />

              <span className={`px-2 py-0.5 rounded text-xs font-semibold text-white ${preset.color}`}>
                {index + 1}
              </span>

              {/* Segment Type Selector */}
              <select
                value={segment.type}
                onChange={(e) => {
                  const newType = e.target.value as SegmentType;
                  const updates: Partial<Segment> = { type: newType };
                  if (newType === 'denovo') {
                    updates.minLength = 50;
                    updates.maxLength = 50;
                  } else if (newType === 'fixed') {
                    updates.chain = 'A';
                    updates.startResidue = 1;
                    updates.endResidue = 50;
                  } else if (newType === 'linker') {
                    updates.gapMin = 0;
                    updates.gapMax = 0;
                  }
                  updateSegment(segment.id, updates);
                }}
                className="bg-card border border-border rounded-lg px-2.5 py-1.5 text-sm text-foreground focus:border-blue-500 focus:ring-1 focus:ring-blue-500 outline-none"
              >
                <option value="denovo">De Novo</option>
                <option value="fixed">Fixed</option>
                <option value="linker">Linker</option>
              </select>

              {/* Segment-specific inputs */}
              {segment.type === 'denovo' && (
                <div className="flex items-center gap-1.5 text-sm">
                  <input
                    type="number"
                    value={segment.minLength || 0}
                    onChange={(e) => updateSegment(segment.id, { minLength: parseInt(e.target.value) || 0 })}
                    className="w-16 bg-card border border-border rounded-lg px-2 py-1.5 text-foreground focus:border-blue-500 focus:ring-1 focus:ring-blue-500 outline-none"
                    min={1}
                    placeholder="Min"
                  />
                  <span className="text-muted-foreground">-</span>
                  <input
                    type="number"
                    value={segment.maxLength || 0}
                    onChange={(e) => updateSegment(segment.id, { maxLength: parseInt(e.target.value) || 0 })}
                    className="w-16 bg-card border border-border rounded-lg px-2 py-1.5 text-foreground focus:border-blue-500 focus:ring-1 focus:ring-blue-500 outline-none"
                    min={1}
                    placeholder="Max"
                  />
                  <span className="text-muted-foreground text-xs">residues</span>
                </div>
              )}

              {segment.type === 'fixed' && (
                <div className="flex items-center gap-1.5 text-sm">
                  <span className="text-muted-foreground">Chain</span>
                  <input
                    type="text"
                    value={segment.chain || 'A'}
                    onChange={(e) => updateSegment(segment.id, { chain: e.target.value.toUpperCase().slice(0, 1) })}
                    className="w-10 bg-card border border-border rounded-lg px-2 py-1.5 text-center text-foreground focus:border-blue-500 focus:ring-1 focus:ring-blue-500 outline-none"
                    maxLength={1}
                  />
                  <span className="text-muted-foreground">:</span>
                  <input
                    type="number"
                    value={segment.startResidue || 1}
                    onChange={(e) => updateSegment(segment.id, { startResidue: parseInt(e.target.value) || 1 })}
                    className="w-16 bg-card border border-border rounded-lg px-2 py-1.5 text-foreground focus:border-blue-500 focus:ring-1 focus:ring-blue-500 outline-none"
                    min={1}
                  />
                  <span className="text-muted-foreground">-</span>
                  <input
                    type="number"
                    value={segment.endResidue || 1}
                    onChange={(e) => updateSegment(segment.id, { endResidue: parseInt(e.target.value) || 1 })}
                    className="w-16 bg-card border border-border rounded-lg px-2 py-1.5 text-foreground focus:border-blue-500 focus:ring-1 focus:ring-blue-500 outline-none"
                    min={1}
                  />
                </div>
              )}

              {segment.type === 'linker' && (
                <div className="flex items-center gap-1.5 text-sm">
                  <span className="text-muted-foreground">Gap</span>
                  <input
                    type="number"
                    value={segment.gapMin ?? 0}
                    onChange={(e) => updateSegment(segment.id, { gapMin: parseInt(e.target.value) || 0 })}
                    className="w-14 bg-card border border-border rounded-lg px-2 py-1.5 text-foreground focus:border-blue-500 focus:ring-1 focus:ring-blue-500 outline-none"
                    min={0}
                  />
                  <span className="text-muted-foreground">-</span>
                  <input
                    type="number"
                    value={segment.gapMax ?? 0}
                    onChange={(e) => updateSegment(segment.id, { gapMax: parseInt(e.target.value) || 0 })}
                    className="w-14 bg-card border border-border rounded-lg px-2 py-1.5 text-foreground focus:border-blue-500 focus:ring-1 focus:ring-blue-500 outline-none"
                    min={0}
                  />
                </div>
              )}

              <div className="flex-1" />

              <button
                onClick={() => removeSegment(segment.id)}
                disabled={segments.length <= 1}
                className="p-1.5 hover:bg-card/80 rounded-lg disabled:opacity-30 transition-colors"
                title="Remove segment"
              >
                <Trash2 className="w-5 h-5 text-red-500" />
              </button>
            </div>
          );
        })}
      </div>

      {/* Add Segment Buttons */}
      <div className="flex gap-2">
        <button
          onClick={() => addSegment('denovo')}
          className="flex-1 py-2.5 bg-blue-50 hover:bg-blue-100 border border-blue-200 hover:border-blue-300 rounded-lg text-sm font-medium text-blue-700 flex items-center justify-center gap-1.5 transition-all"
        >
          <Plus className="w-4 h-4" /> De Novo
        </button>
        <button
          onClick={() => addSegment('fixed')}
          className="flex-1 py-2.5 bg-emerald-50 hover:bg-emerald-100 border border-emerald-200 hover:border-emerald-300 rounded-lg text-sm font-medium text-emerald-700 flex items-center justify-center gap-1.5 transition-all"
        >
          <Plus className="w-4 h-4" /> Fixed
        </button>
        <button
          onClick={() => addSegment('linker')}
          className="flex-1 py-2.5 bg-violet-50 hover:bg-violet-100 border border-violet-200 hover:border-violet-300 rounded-lg text-sm font-medium text-violet-700 flex items-center justify-center gap-1.5 transition-all"
        >
          <Plus className="w-4 h-4" /> Linker
        </button>
      </div>

      {/* Preview */}
      <div className="bg-muted rounded-lg p-4 border border-border">
        <div className="flex items-center justify-between mb-3">
          <span className="text-xs font-semibold text-muted-foreground uppercase tracking-wider">Generated Contig</span>
          <button
            onClick={applyContig}
            className="px-4 py-1.5 bg-blue-600 hover:bg-blue-700 rounded-lg text-xs font-semibold text-white shadow-sm transition-colors"
          >
            Apply to Design
          </button>
        </div>
        <code className="text-sm text-blue-700 font-mono bg-card px-3 py-2 rounded-lg border border-border block">
          {contigString || '(empty)'}
        </code>
      </div>
    </div>
  );
}
