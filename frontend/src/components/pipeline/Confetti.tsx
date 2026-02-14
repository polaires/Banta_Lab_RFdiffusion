'use client';

import { useEffect, useState } from 'react';

const CONFETTI_COLORS = [
  '#D4922A', // amber/gold
  '#339966', // sage green
  '#5588AA', // dusty blue
  '#C45544', // warm rose
  '#8B6FC0', // soft purple
  '#E09520', // rich amber
];

interface Particle {
  id: number;
  x: number;
  y: number;
  color: string;
  rotation: number;
  size: number;
  delay: number;
}

/**
 * Lightweight confetti burst for pipeline completion celebrations.
 * Pure CSS animations — no canvas, no dependencies.
 */
export function Confetti() {
  const [particles, setParticles] = useState<Particle[]>([]);

  useEffect(() => {
    const items: Particle[] = Array.from({ length: 24 }, (_, i) => ({
      id: i,
      x: 40 + Math.random() * 20, // cluster around center
      y: 30 + Math.random() * 10,
      color: CONFETTI_COLORS[i % CONFETTI_COLORS.length],
      rotation: Math.random() * 360,
      size: 4 + Math.random() * 4,
      delay: Math.random() * 0.3,
    }));
    setParticles(items);
  }, []);

  if (particles.length === 0) return null;

  return (
    <div
      className="absolute inset-0 pointer-events-none overflow-hidden z-10"
      aria-hidden="true"
    >
      {particles.map((p) => (
        <div
          key={p.id}
          className="absolute rounded-sm"
          style={{
            left: `${p.x}%`,
            top: `${p.y}%`,
            width: p.size,
            height: p.size,
            backgroundColor: p.color,
            transform: `rotate(${p.rotation}deg)`,
            animation: `confetti-burst 1.8s ease-out ${p.delay}s forwards`,
            opacity: 0,
          }}
        />
      ))}
    </div>
  );
}
