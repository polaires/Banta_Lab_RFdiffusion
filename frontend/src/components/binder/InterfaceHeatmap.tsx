'use client';

import { useRef, useEffect, useState } from 'react';
import * as d3 from 'd3';
import { LayoutGrid, Info } from 'lucide-react';

interface Contact {
  target_res: number;
  binder_res: number;
  distance: number;
  contact_type?: 'hbond' | 'hydrophobic' | 'ionic' | 'other';
}

interface InterfaceHeatmapProps {
  contacts: Contact[];
  targetResidueRange: [number, number];  // [start, end]
  binderResidueRange: [number, number];  // [start, end]
  targetChainLabel?: string;
  binderChainLabel?: string;
  onCellClick?: (targetRes: number, binderRes: number) => void;
  onCellHover?: (targetRes: number, binderRes: number, distance: number) => void;
  width?: number;
  height?: number;
}

export function InterfaceHeatmap({
  contacts,
  targetResidueRange,
  binderResidueRange,
  targetChainLabel = 'Target',
  binderChainLabel = 'Binder',
  onCellClick,
  onCellHover,
  width = 400,
  height = 300,
}: InterfaceHeatmapProps) {
  const svgRef = useRef<SVGSVGElement>(null);
  const tooltipRef = useRef<HTMLDivElement>(null);
  const [hoveredCell, setHoveredCell] = useState<{ target: number; binder: number; distance: number } | null>(null);

  useEffect(() => {
    if (!svgRef.current || contacts.length === 0) return;

    const svg = d3.select(svgRef.current);
    svg.selectAll('*').remove();

    const margin = { top: 40, right: 20, bottom: 60, left: 60 };
    const innerWidth = width - margin.left - margin.right;
    const innerHeight = height - margin.top - margin.bottom;

    // Create scales
    const targetResidues = d3.range(targetResidueRange[0], targetResidueRange[1] + 1);
    const binderResidues = d3.range(binderResidueRange[0], binderResidueRange[1] + 1);

    const xScale = d3.scaleBand()
      .domain(targetResidues.map(String))
      .range([0, innerWidth])
      .padding(0.05);

    const yScale = d3.scaleBand()
      .domain(binderResidues.map(String))
      .range([0, innerHeight])
      .padding(0.05);

    // Color scale (closer = more intense)
    const colorScale = d3.scaleSequential(d3.interpolateBlues)
      .domain([8, 2]);  // 8Å = light, 2Å = dark

    // Create container group
    const g = svg.append('g')
      .attr('transform', `translate(${margin.left},${margin.top})`);

    // Create contact matrix
    const contactMatrix: Map<string, number> = new Map();
    contacts.forEach(c => {
      const key = `${c.target_res}-${c.binder_res}`;
      const existing = contactMatrix.get(key);
      if (!existing || c.distance < existing) {
        contactMatrix.set(key, c.distance);
      }
    });

    // Draw cells
    const cells = g.selectAll('rect.cell')
      .data(contacts)
      .enter()
      .append('rect')
      .attr('class', 'cell')
      .attr('x', d => xScale(String(d.target_res)) || 0)
      .attr('y', d => yScale(String(d.binder_res)) || 0)
      .attr('width', xScale.bandwidth())
      .attr('height', yScale.bandwidth())
      .attr('fill', d => colorScale(d.distance))
      .attr('rx', 1)
      .attr('stroke', 'none')
      .style('cursor', 'pointer')
      .on('mouseover', function(event, d) {
        d3.select(this)
          .attr('stroke', '#7c3aed')
          .attr('stroke-width', 2);
        setHoveredCell({ target: d.target_res, binder: d.binder_res, distance: d.distance });
        onCellHover?.(d.target_res, d.binder_res, d.distance);
      })
      .on('mouseout', function() {
        d3.select(this)
          .attr('stroke', 'none');
        setHoveredCell(null);
      })
      .on('click', (event, d) => {
        onCellClick?.(d.target_res, d.binder_res);
      });

    // X-axis (target residues)
    const xAxisG = g.append('g')
      .attr('transform', `translate(0,${innerHeight})`)
      .call(d3.axisBottom(xScale)
        .tickValues(xScale.domain().filter((_, i) => i % Math.ceil(targetResidues.length / 10) === 0))
      );

    xAxisG.selectAll('text')
      .attr('font-size', '10px')
      .attr('transform', 'rotate(-45)')
      .attr('text-anchor', 'end');

    // Y-axis (binder residues)
    const yAxisG = g.append('g')
      .call(d3.axisLeft(yScale)
        .tickValues(yScale.domain().filter((_, i) => i % Math.ceil(binderResidues.length / 10) === 0))
      );

    yAxisG.selectAll('text')
      .attr('font-size', '10px');

    // X-axis label
    svg.append('text')
      .attr('x', margin.left + innerWidth / 2)
      .attr('y', height - 10)
      .attr('text-anchor', 'middle')
      .attr('font-size', '12px')
      .attr('fill', '#64748b')
      .text(`${targetChainLabel} Residues`);

    // Y-axis label
    svg.append('text')
      .attr('transform', 'rotate(-90)')
      .attr('x', -(margin.top + innerHeight / 2))
      .attr('y', 15)
      .attr('text-anchor', 'middle')
      .attr('font-size', '12px')
      .attr('fill', '#64748b')
      .text(`${binderChainLabel} Residues`);

    // Title
    svg.append('text')
      .attr('x', width / 2)
      .attr('y', 20)
      .attr('text-anchor', 'middle')
      .attr('font-size', '14px')
      .attr('font-weight', '600')
      .attr('fill', '#1e293b')
      .text('Interface Contact Map');

  }, [contacts, targetResidueRange, binderResidueRange, width, height, targetChainLabel, binderChainLabel, onCellClick, onCellHover]);

  return (
    <div className="bg-white rounded-xl border border-slate-200 overflow-hidden">
      {/* Header */}
      <div className="bg-gradient-to-r from-indigo-50 to-purple-50 px-4 py-3 border-b border-indigo-100">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-2">
            <LayoutGrid className="h-5 w-5 text-indigo-600" />
            <h4 className="font-semibold text-slate-900 text-sm">Interface Heatmap</h4>
          </div>
          <div className="text-xs text-slate-500">
            {contacts.length} contacts
          </div>
        </div>
      </div>

      {/* Heatmap */}
      <div className="p-4 relative">
        <svg
          ref={svgRef}
          width={width}
          height={height}
          className="mx-auto"
        />

        {/* Tooltip */}
        {hoveredCell && (
          <div
            ref={tooltipRef}
            className="absolute bg-slate-900 text-white px-3 py-2 rounded-lg text-xs shadow-lg pointer-events-none z-10"
            style={{
              left: '50%',
              top: '50%',
              transform: 'translate(-50%, -50%)',
            }}
          >
            <div className="font-medium">
              {targetChainLabel}:{hoveredCell.target} ↔ {binderChainLabel}:{hoveredCell.binder}
            </div>
            <div className="text-slate-300 mt-1">
              Distance: {hoveredCell.distance.toFixed(2)} Å
            </div>
          </div>
        )}
      </div>

      {/* Legend */}
      <div className="px-4 pb-4">
        <div className="flex items-center justify-center gap-4 text-xs">
          <span className="text-slate-500">Distance:</span>
          <div className="flex items-center gap-1">
            <div className="w-4 h-3 rounded" style={{ backgroundColor: d3.interpolateBlues(0.9) }} />
            <span className="text-slate-600">Close (2Å)</span>
          </div>
          <div className="flex items-center gap-1">
            <div className="w-4 h-3 rounded" style={{ backgroundColor: d3.interpolateBlues(0.5) }} />
            <span className="text-slate-600">Medium</span>
          </div>
          <div className="flex items-center gap-1">
            <div className="w-4 h-3 rounded" style={{ backgroundColor: d3.interpolateBlues(0.2) }} />
            <span className="text-slate-600">Far (8Å)</span>
          </div>
        </div>
      </div>

      {/* Instructions */}
      <div className="bg-slate-50 px-4 py-2 border-t border-slate-200 text-xs text-slate-500 flex items-center gap-1">
        <Info className="h-3 w-3" />
        Hover over cells to see contact details. Click to highlight in 3D viewer.
      </div>
    </div>
  );
}
