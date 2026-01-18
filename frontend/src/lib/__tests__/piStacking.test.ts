import { describe, it, expect } from 'vitest';
import {
  calculateRingGeometry,
  analyzePiStacking,
  type RingGeometry,
  type PiStackingResult,
} from '../ligandAnalysis';

describe('calculateRingGeometry', () => {
  it('should calculate centroid of a hexagonal ring', () => {
    // Regular hexagon in XY plane centered at origin
    const positions: [number, number, number][] = [
      [1, 0, 0],
      [0.5, 0.866, 0],
      [-0.5, 0.866, 0],
      [-1, 0, 0],
      [-0.5, -0.866, 0],
      [0.5, -0.866, 0],
    ];

    const result = calculateRingGeometry(positions);

    expect(result).not.toBeNull();
    expect(result!.centroid[0]).toBeCloseTo(0, 1);
    expect(result!.centroid[1]).toBeCloseTo(0, 1);
    expect(result!.centroid[2]).toBeCloseTo(0, 1);
  });

  it('should calculate normal vector perpendicular to ring plane', () => {
    // Ring in XY plane should have normal along Z
    const positions: [number, number, number][] = [
      [1, 0, 0],
      [0.5, 0.866, 0],
      [-0.5, 0.866, 0],
      [-1, 0, 0],
      [-0.5, -0.866, 0],
      [0.5, -0.866, 0],
    ];

    const result = calculateRingGeometry(positions);

    expect(result).not.toBeNull();
    // Normal should be [0, 0, 1] or [0, 0, -1]
    expect(Math.abs(result!.normal[2])).toBeCloseTo(1, 1);
  });

  it('should return null for fewer than 5 atoms', () => {
    const positions: [number, number, number][] = [
      [0, 0, 0],
      [1, 0, 0],
      [1, 1, 0],
    ];

    const result = calculateRingGeometry(positions);

    expect(result).toBeNull();
  });
});

describe('analyzePiStacking', () => {
  it('should detect parallel stacking for rings with parallel normals', () => {
    const ring1: RingGeometry = {
      centroid: [0, 0, 0],
      normal: [0, 0, 1],
    };
    const ring2: RingGeometry = {
      centroid: [0, 0, 3.5],
      normal: [0, 0, 1],
    };

    const result = analyzePiStacking(ring1, ring2);

    expect(result.type).toBe('parallel');
    expect(result.isStacking).toBe(true);
    expect(result.distance).toBeCloseTo(3.5, 1);
    expect(result.angle).toBeCloseTo(0, 5);
  });

  it('should detect t-shaped stacking for perpendicular rings', () => {
    const ring1: RingGeometry = {
      centroid: [0, 0, 0],
      normal: [0, 0, 1],
    };
    const ring2: RingGeometry = {
      centroid: [3, 0, 0],
      normal: [1, 0, 0],
    };

    const result = analyzePiStacking(ring1, ring2);

    expect(result.type).toBe('t-shaped');
    expect(result.isStacking).toBe(true);
    expect(result.angle).toBeCloseTo(90, 5);
  });

  it('should detect offset-parallel for displaced parallel rings', () => {
    const ring1: RingGeometry = {
      centroid: [0, 0, 0],
      normal: [0, 0, 1],
    };
    const ring2: RingGeometry = {
      centroid: [2.5, 0, 3.5],
      normal: [0, 0, 1],
    };

    const result = analyzePiStacking(ring1, ring2);

    expect(result.type).toBe('offset-parallel');
    expect(result.isStacking).toBe(true);
  });

  it('should return none for distant rings', () => {
    const ring1: RingGeometry = {
      centroid: [0, 0, 0],
      normal: [0, 0, 1],
    };
    const ring2: RingGeometry = {
      centroid: [10, 0, 0],
      normal: [0, 0, 1],
    };

    const result = analyzePiStacking(ring1, ring2);

    expect(result.type).toBe('none');
    expect(result.isStacking).toBe(false);
  });
});
