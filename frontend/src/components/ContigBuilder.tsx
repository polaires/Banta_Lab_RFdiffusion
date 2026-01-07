'use client';

import { useState } from 'react';
import { Plus, Trash2, GripVertical, Wand2, HelpCircle, Zap } from 'lucide-react';

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
  denovo: { label: 'De Novo Region', description: 'Design new residues from scratch', color: 'bg-blue-600' },
  fixed: { label: 'Fixed Region', description: 'Keep existing residues from input PDB', color: 'bg-green-600' },
  linker: { label: 'Flexible Linker', description: 'Variable-length connection between segments', color: 'bg-purple-600' },
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
    // Generate new IDs for segments
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

  // Apply to parent
  const applyContig = () => {
    onContigChange(contigString);
  };

  return (
    <div className="bg-gray-800/50 rounded-lg p-4 space-y-4 border border-gray-700">
      <div className="flex items-center justify-between">
        <div className="flex items-center gap-2">
          <Wand2 className="w-5 h-5 text-blue-400" />
          <h3 className="font-semibold text-gray-200">Contig Builder</h3>
        </div>
        <button
          onClick={() => setShowHelp(!showHelp)}
          className="p-1 hover:bg-gray-700 rounded"
          title="Help"
        >
          <HelpCircle className="w-4 h-4 text-gray-400" />
        </button>
      </div>

      {showHelp && (
        <div className="bg-gray-900/50 p-3 rounded text-sm text-gray-300 space-y-2">
          <p><strong>De Novo Region:</strong> Generates new amino acids. Specify length or range.</p>
          <p><strong>Fixed Region:</strong> Keeps residues from input PDB. Specify chain and residue numbers.</p>
          <p><strong>Flexible Linker:</strong> Gap between segments. Use 0 for direct connection.</p>
          <p className="pt-2 border-t border-gray-700 mt-2">
            <strong>Binder Design:</strong> Upload target PDB, add Fixed region (target), Linker (/0), De Novo region (binder).
          </p>
        </div>
      )}

      {/* Design Presets */}
      <div className="space-y-2">
        <button
          onClick={() => setShowPresets(!showPresets)}
          className="flex items-center gap-2 text-sm text-blue-400 hover:text-blue-300"
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
                className="p-2 text-left bg-gray-700/50 hover:bg-gray-700 rounded border border-gray-600 transition"
              >
                <div className="font-medium text-sm text-gray-200">{preset.name}</div>
                <div className="text-xs text-gray-400">{preset.description}</div>
              </button>
            ))}
          </div>
        )}
      </div>

      {/* Segment List */}
      <div className="space-y-2">
        {segments.map((segment, index) => (
          <div
            key={segment.id}
            className={`flex items-center gap-2 p-3 rounded border ${
              segment.type === 'denovo' ? 'border-blue-600/50 bg-blue-900/20' :
              segment.type === 'fixed' ? 'border-green-600/50 bg-green-900/20' :
              'border-purple-600/50 bg-purple-900/20'
            }`}
          >
            <GripVertical className="w-4 h-4 text-gray-500 cursor-grab" />

            <span className={`px-2 py-0.5 rounded text-xs font-medium ${SEGMENT_PRESETS[segment.type].color}`}>
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
              className="bg-gray-700 border border-gray-600 rounded px-2 py-1 text-sm"
            >
              <option value="denovo">De Novo</option>
              <option value="fixed">Fixed</option>
              <option value="linker">Linker</option>
            </select>

            {/* Segment-specific inputs */}
            {segment.type === 'denovo' && (
              <div className="flex items-center gap-1 text-sm">
                <input
                  type="number"
                  value={segment.minLength || 0}
                  onChange={(e) => updateSegment(segment.id, { minLength: parseInt(e.target.value) || 0 })}
                  className="w-16 bg-gray-700 border border-gray-600 rounded px-2 py-1"
                  min={1}
                  placeholder="Min"
                />
                <span className="text-gray-500">-</span>
                <input
                  type="number"
                  value={segment.maxLength || 0}
                  onChange={(e) => updateSegment(segment.id, { maxLength: parseInt(e.target.value) || 0 })}
                  className="w-16 bg-gray-700 border border-gray-600 rounded px-2 py-1"
                  min={1}
                  placeholder="Max"
                />
                <span className="text-gray-400 text-xs">residues</span>
              </div>
            )}

            {segment.type === 'fixed' && (
              <div className="flex items-center gap-1 text-sm">
                <span className="text-gray-400">Chain</span>
                <input
                  type="text"
                  value={segment.chain || 'A'}
                  onChange={(e) => updateSegment(segment.id, { chain: e.target.value.toUpperCase().slice(0, 1) })}
                  className="w-10 bg-gray-700 border border-gray-600 rounded px-2 py-1 text-center"
                  maxLength={1}
                />
                <span className="text-gray-400">:</span>
                <input
                  type="number"
                  value={segment.startResidue || 1}
                  onChange={(e) => updateSegment(segment.id, { startResidue: parseInt(e.target.value) || 1 })}
                  className="w-16 bg-gray-700 border border-gray-600 rounded px-2 py-1"
                  min={1}
                />
                <span className="text-gray-500">-</span>
                <input
                  type="number"
                  value={segment.endResidue || 1}
                  onChange={(e) => updateSegment(segment.id, { endResidue: parseInt(e.target.value) || 1 })}
                  className="w-16 bg-gray-700 border border-gray-600 rounded px-2 py-1"
                  min={1}
                />
              </div>
            )}

            {segment.type === 'linker' && (
              <div className="flex items-center gap-1 text-sm">
                <span className="text-gray-400">Gap</span>
                <input
                  type="number"
                  value={segment.gapMin ?? 0}
                  onChange={(e) => updateSegment(segment.id, { gapMin: parseInt(e.target.value) || 0 })}
                  className="w-14 bg-gray-700 border border-gray-600 rounded px-2 py-1"
                  min={0}
                />
                <span className="text-gray-500">-</span>
                <input
                  type="number"
                  value={segment.gapMax ?? 0}
                  onChange={(e) => updateSegment(segment.id, { gapMax: parseInt(e.target.value) || 0 })}
                  className="w-14 bg-gray-700 border border-gray-600 rounded px-2 py-1"
                  min={0}
                />
              </div>
            )}

            <div className="flex-1" />

            <button
              onClick={() => removeSegment(segment.id)}
              disabled={segments.length <= 1}
              className="p-1 hover:bg-gray-700 rounded disabled:opacity-30"
              title="Remove segment"
            >
              <Trash2 className="w-4 h-4 text-red-400" />
            </button>
          </div>
        ))}
      </div>

      {/* Add Segment Buttons */}
      <div className="flex gap-2">
        <button
          onClick={() => addSegment('denovo')}
          className="flex-1 py-2 bg-blue-600/20 hover:bg-blue-600/30 border border-blue-600/50 rounded text-sm flex items-center justify-center gap-1"
        >
          <Plus className="w-4 h-4" /> De Novo
        </button>
        <button
          onClick={() => addSegment('fixed')}
          className="flex-1 py-2 bg-green-600/20 hover:bg-green-600/30 border border-green-600/50 rounded text-sm flex items-center justify-center gap-1"
        >
          <Plus className="w-4 h-4" /> Fixed
        </button>
        <button
          onClick={() => addSegment('linker')}
          className="flex-1 py-2 bg-purple-600/20 hover:bg-purple-600/30 border border-purple-600/50 rounded text-sm flex items-center justify-center gap-1"
        >
          <Plus className="w-4 h-4" /> Linker
        </button>
      </div>

      {/* Preview */}
      <div className="bg-gray-900 rounded p-3 border border-gray-700">
        <div className="flex items-center justify-between mb-2">
          <span className="text-xs font-medium text-gray-400">Generated Contig</span>
          <button
            onClick={applyContig}
            className="px-3 py-1 bg-blue-600 hover:bg-blue-700 rounded text-xs font-medium"
          >
            Use in Design
          </button>
        </div>
        <code className="text-sm text-green-400 font-mono break-all">
          {contigString || '(empty)'}
        </code>
      </div>
    </div>
  );
}
