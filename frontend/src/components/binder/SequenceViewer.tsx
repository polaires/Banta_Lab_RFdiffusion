'use client';

import { useState, useMemo } from 'react';
import { Dna, Copy, Check } from 'lucide-react';

interface SequenceViewerProps {
  sequence: string;
  interfaceResidues?: number[];  // 1-indexed residue numbers at the interface
  highlightedResidue?: number;   // Currently highlighted residue
  onResidueClick?: (residue: number, aminoAcid: string) => void;
  onResidueHover?: (residue: number | null, aminoAcid: string | null) => void;
  title?: string;
  showNumbers?: boolean;
  chunkSize?: number;
}

// Amino acid properties for coloring
const AA_PROPERTIES: Record<string, { color: string; type: string; name: string }> = {
  // Hydrophobic
  A: { color: 'bg-amber-100 text-amber-800', type: 'hydrophobic', name: 'Alanine' },
  V: { color: 'bg-amber-100 text-amber-800', type: 'hydrophobic', name: 'Valine' },
  L: { color: 'bg-amber-100 text-amber-800', type: 'hydrophobic', name: 'Leucine' },
  I: { color: 'bg-amber-100 text-amber-800', type: 'hydrophobic', name: 'Isoleucine' },
  M: { color: 'bg-amber-100 text-amber-800', type: 'hydrophobic', name: 'Methionine' },
  F: { color: 'bg-amber-200 text-amber-900', type: 'aromatic', name: 'Phenylalanine' },
  W: { color: 'bg-amber-200 text-amber-900', type: 'aromatic', name: 'Tryptophan' },
  Y: { color: 'bg-amber-200 text-amber-900', type: 'aromatic', name: 'Tyrosine' },
  P: { color: 'bg-amber-100 text-amber-800', type: 'hydrophobic', name: 'Proline' },
  G: { color: 'bg-slate-100 text-slate-700', type: 'special', name: 'Glycine' },
  // Polar
  S: { color: 'bg-green-100 text-green-800', type: 'polar', name: 'Serine' },
  T: { color: 'bg-green-100 text-green-800', type: 'polar', name: 'Threonine' },
  N: { color: 'bg-green-100 text-green-800', type: 'polar', name: 'Asparagine' },
  Q: { color: 'bg-green-100 text-green-800', type: 'polar', name: 'Glutamine' },
  C: { color: 'bg-yellow-100 text-yellow-800', type: 'special', name: 'Cysteine' },
  // Positive
  K: { color: 'bg-blue-100 text-blue-800', type: 'positive', name: 'Lysine' },
  R: { color: 'bg-blue-100 text-blue-800', type: 'positive', name: 'Arginine' },
  H: { color: 'bg-blue-100 text-blue-800', type: 'positive', name: 'Histidine' },
  // Negative
  D: { color: 'bg-red-100 text-red-800', type: 'negative', name: 'Aspartate' },
  E: { color: 'bg-red-100 text-red-800', type: 'negative', name: 'Glutamate' },
};

export function SequenceViewer({
  sequence,
  interfaceResidues = [],
  highlightedResidue,
  onResidueClick,
  onResidueHover,
  title = 'Binder Sequence',
  showNumbers = true,
  chunkSize = 10,
}: SequenceViewerProps) {
  const [hoveredResidue, setHoveredResidue] = useState<number | null>(null);
  const [copied, setCopied] = useState(false);

  const interfaceSet = useMemo(() => new Set(interfaceResidues), [interfaceResidues]);

  const handleCopy = () => {
    navigator.clipboard.writeText(sequence);
    setCopied(true);
    setTimeout(() => setCopied(false), 2000);
  };

  const handleResidueHover = (idx: number | null) => {
    setHoveredResidue(idx);
    if (idx !== null) {
      onResidueHover?.(idx + 1, sequence[idx]);
    } else {
      onResidueHover?.(null, null);
    }
  };

  // Split sequence into chunks for better readability
  const chunks = useMemo(() => {
    const result: { start: number; residues: string[] }[] = [];
    for (let i = 0; i < sequence.length; i += chunkSize) {
      result.push({
        start: i,
        residues: sequence.slice(i, i + chunkSize).split(''),
      });
    }
    return result;
  }, [sequence, chunkSize]);

  // Calculate statistics
  const stats = useMemo(() => {
    const counts: Record<string, number> = {};
    for (const aa of sequence) {
      counts[aa] = (counts[aa] || 0) + 1;
    }

    const types = {
      hydrophobic: 0,
      polar: 0,
      positive: 0,
      negative: 0,
      aromatic: 0,
      special: 0,
    };

    for (const [aa, count] of Object.entries(counts)) {
      const prop = AA_PROPERTIES[aa];
      if (prop) {
        types[prop.type as keyof typeof types] += count;
      }
    }

    return {
      length: sequence.length,
      interfaceCount: interfaceResidues.length,
      types,
    };
  }, [sequence, interfaceResidues]);

  return (
    <div className="bg-white rounded-xl border border-slate-200 overflow-hidden">
      {/* Header */}
      <div className="bg-gradient-to-r from-emerald-50 to-teal-50 px-4 py-3 border-b border-emerald-100">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-2">
            <Dna className="h-5 w-5 text-emerald-600" />
            <h4 className="font-semibold text-slate-900 text-sm">{title}</h4>
            <span className="text-xs text-slate-500 bg-slate-200 px-2 py-0.5 rounded-full">
              {sequence.length} aa
            </span>
          </div>
          <button
            onClick={handleCopy}
            className="flex items-center gap-1 px-3 py-1 text-xs font-medium text-emerald-600 hover:bg-emerald-100 rounded-lg transition-colors"
          >
            {copied ? (
              <Check className="h-4 w-4" />
            ) : (
              <Copy className="h-4 w-4" />
            )}
            {copied ? 'Copied!' : 'Copy'}
          </button>
        </div>
      </div>

      {/* Sequence Display */}
      <div className="p-4 font-mono text-sm overflow-x-auto">
        <div className="space-y-1">
          {chunks.map((chunk, chunkIdx) => (
            <div key={chunkIdx} className="flex items-center gap-2">
              {/* Line number */}
              {showNumbers && (
                <span className="text-xs text-slate-400 w-8 text-right select-none">
                  {chunk.start + 1}
                </span>
              )}

              {/* Residues */}
              <div className="flex gap-px">
                {chunk.residues.map((aa, i) => {
                  const globalIdx = chunk.start + i;
                  const residueNum = globalIdx + 1;
                  const isInterface = interfaceSet.has(residueNum);
                  const isHighlighted = highlightedResidue === residueNum;
                  const isHovered = hoveredResidue === globalIdx;
                  const props = AA_PROPERTIES[aa] || { color: 'bg-gray-100 text-gray-700', type: 'unknown', name: 'Unknown' };

                  return (
                    <button
                      key={globalIdx}
                      onClick={() => onResidueClick?.(residueNum, aa)}
                      onMouseEnter={() => handleResidueHover(globalIdx)}
                      onMouseLeave={() => handleResidueHover(null)}
                      className={`
                        w-5 h-6 flex items-center justify-center rounded text-xs font-medium transition-all
                        ${props.color}
                        ${isInterface ? 'ring-2 ring-violet-500 ring-offset-1' : ''}
                        ${isHighlighted ? 'ring-2 ring-orange-500 ring-offset-1 scale-110' : ''}
                        ${isHovered ? 'scale-110 shadow-md z-10' : ''}
                        hover:scale-110 hover:shadow-md cursor-pointer
                      `}
                      title={`${residueNum}: ${props.name} (${aa})${isInterface ? ' - Interface' : ''}`}
                    >
                      {aa}
                    </button>
                  );
                })}
              </div>

              {/* Chunk separator */}
              {chunkIdx < chunks.length - 1 && (
                <span className="text-slate-300 select-none">|</span>
              )}
            </div>
          ))}
        </div>
      </div>

      {/* Hovered residue info */}
      {hoveredResidue !== null && (
        <div className="px-4 pb-2">
          <div className="bg-slate-100 rounded-lg px-3 py-2 text-xs">
            <span className="font-medium text-slate-900">
              Position {hoveredResidue + 1}:
            </span>
            <span className="ml-2 text-slate-600">
              {AA_PROPERTIES[sequence[hoveredResidue]]?.name || 'Unknown'} ({sequence[hoveredResidue]})
            </span>
            {interfaceSet.has(hoveredResidue + 1) && (
              <span className="ml-2 px-2 py-0.5 rounded-full bg-violet-100 text-violet-700 text-xs">
                Interface
              </span>
            )}
          </div>
        </div>
      )}

      {/* Legend */}
      <div className="px-4 py-3 border-t border-slate-100 bg-slate-50">
        <div className="flex flex-wrap items-center gap-3 text-xs">
          <span className="text-slate-500 font-medium">Residue types:</span>
          <div className="flex items-center gap-1">
            <div className="w-3 h-3 rounded bg-amber-100" />
            <span className="text-slate-600">Hydrophobic</span>
          </div>
          <div className="flex items-center gap-1">
            <div className="w-3 h-3 rounded bg-green-100" />
            <span className="text-slate-600">Polar</span>
          </div>
          <div className="flex items-center gap-1">
            <div className="w-3 h-3 rounded bg-blue-100" />
            <span className="text-slate-600">Positive</span>
          </div>
          <div className="flex items-center gap-1">
            <div className="w-3 h-3 rounded bg-red-100" />
            <span className="text-slate-600">Negative</span>
          </div>
          <div className="flex items-center gap-1">
            <div className="w-3 h-3 rounded ring-2 ring-violet-500" />
            <span className="text-slate-600">Interface</span>
          </div>
        </div>
      </div>

      {/* Statistics */}
      <div className="px-4 py-3 border-t border-slate-200 bg-white">
        <div className="grid grid-cols-4 gap-4 text-xs">
          <div>
            <div className="text-slate-400 uppercase tracking-wide">Length</div>
            <div className="font-semibold text-slate-900">{stats.length} aa</div>
          </div>
          <div>
            <div className="text-slate-400 uppercase tracking-wide">Interface</div>
            <div className="font-semibold text-violet-600">{stats.interfaceCount} residues</div>
          </div>
          <div>
            <div className="text-slate-400 uppercase tracking-wide">Hydrophobic</div>
            <div className="font-semibold text-amber-600">
              {((stats.types.hydrophobic + stats.types.aromatic) / stats.length * 100).toFixed(0)}%
            </div>
          </div>
          <div>
            <div className="text-slate-400 uppercase tracking-wide">Charged</div>
            <div className="font-semibold text-blue-600">
              {((stats.types.positive + stats.types.negative) / stats.length * 100).toFixed(0)}%
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}
