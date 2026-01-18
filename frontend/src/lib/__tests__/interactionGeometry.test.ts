import { describe, it, expect } from 'vitest';
import {
  computeInteractionLines,
  type InteractionLine,
} from '../interactionGeometry';

describe('computeInteractionLines', () => {
  const samplePdb = `
ATOM      1  N   SER A  10       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  SER A  10       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  OG  SER A  10       2.000   1.000   0.000  1.00  0.00           O
HETATM   10  C1  LIG B   1       4.000   1.500   0.000  1.00  0.00           C
HETATM   11  O1  LIG B   1       5.000   1.000   0.000  1.00  0.00           O
END
`;

  it('should return empty array when no ligand contacts', () => {
    const result = computeInteractionLines('', []);
    expect(result).toEqual([]);
  });

  it('should compute line endpoints for hydrogen bond contacts', () => {
    const contacts = [
      {
        ligandChain: 'B',
        ligandResSeq: 1,
        ligandAtom: 'O1',
        proteinResidue: 'SER10',
        proteinChain: 'A',
        proteinAtom: 'OG',
        interactionType: 'hydrogen_bond' as const,
        distance: 3.0,
      },
    ];

    const result = computeInteractionLines(samplePdb, contacts);

    expect(result).toHaveLength(1);
    expect(result[0].type).toBe('hydrogen_bond');
    expect(result[0].start).toBeDefined();
    expect(result[0].end).toBeDefined();
    expect(result[0].color).toBeDefined();
  });

  it('should assign correct colors based on interaction type', () => {
    const contacts = [
      {
        ligandChain: 'B',
        ligandResSeq: 1,
        ligandAtom: 'O1',
        proteinResidue: 'SER10',
        proteinChain: 'A',
        proteinAtom: 'OG',
        interactionType: 'hydrogen_bond' as const,
        distance: 3.0,
      },
    ];

    const result = computeInteractionLines(samplePdb, contacts);

    expect(result[0].color).toBe(0x3B82F6); // Blue for H-bonds
  });
});
