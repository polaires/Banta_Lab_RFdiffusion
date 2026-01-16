import { describe, it, expect } from 'vitest';
import {
  extractPharmacophoreFeatures,
  type LigandContact,
  type PharmacophoreFeature,
  type PharmacophoreType,
} from '../ligandAnalysis';

describe('extractPharmacophoreFeatures', () => {
  describe('H-bond donors', () => {
    it('should classify N atoms in hydrogen bonds as donors', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'SER123',
          chain: 'A',
          atom: 'N',
          distance: 3.0,
          interactionType: 'hydrogen_bond',
        },
      ];

      const features = extractPharmacophoreFeatures(contacts);

      expect(features).toHaveLength(1);
      expect(features[0].type).toBe('donor');
      expect(features[0].residue).toBe('SER123');
      expect(features[0].chain).toBe('A');
      expect(features[0].atom).toBe('N');
      expect(features[0].distance).toBe(3.0);
    });

    it('should classify NH atoms (like NE2, ND2) as donors in hydrogen bonds', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'GLN45',
          chain: 'B',
          atom: 'NE2',
          distance: 2.8,
          interactionType: 'hydrogen_bond',
        },
      ];

      const features = extractPharmacophoreFeatures(contacts);

      expect(features).toHaveLength(1);
      expect(features[0].type).toBe('donor');
    });
  });

  describe('H-bond acceptors', () => {
    it('should classify O atoms in hydrogen bonds as acceptors', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'ASP67',
          chain: 'A',
          atom: 'O',
          distance: 3.2,
          interactionType: 'hydrogen_bond',
        },
      ];

      const features = extractPharmacophoreFeatures(contacts);

      expect(features).toHaveLength(1);
      expect(features[0].type).toBe('acceptor');
      expect(features[0].atom).toBe('O');
    });

    it('should classify OD1, OE1 etc. as acceptors in hydrogen bonds', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'GLU89',
          chain: 'C',
          atom: 'OE1',
          distance: 2.9,
          interactionType: 'hydrogen_bond',
        },
      ];

      const features = extractPharmacophoreFeatures(contacts);

      expect(features).toHaveLength(1);
      expect(features[0].type).toBe('acceptor');
    });
  });

  describe('Aromatic features', () => {
    it('should classify PHE residues as aromatic', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'PHE101',
          chain: 'A',
          atom: 'CG',
          distance: 3.8,
          interactionType: 'pi_stacking',
        },
      ];

      const features = extractPharmacophoreFeatures(contacts);

      expect(features).toHaveLength(1);
      expect(features[0].type).toBe('aromatic');
    });

    it('should classify TYR residues as aromatic', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'TYR55',
          chain: 'A',
          atom: 'CE1',
          distance: 3.6,
          interactionType: 'pi_stacking',
        },
      ];

      const features = extractPharmacophoreFeatures(contacts);

      expect(features).toHaveLength(1);
      expect(features[0].type).toBe('aromatic');
    });

    it('should classify TRP residues as aromatic', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'TRP200',
          chain: 'B',
          atom: 'CD2',
          distance: 3.9,
          interactionType: 'hydrophobic',
        },
      ];

      const features = extractPharmacophoreFeatures(contacts);

      expect(features).toHaveLength(1);
      expect(features[0].type).toBe('aromatic');
    });

    it('should classify HIS residues with pi_stacking as aromatic', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'HIS77',
          chain: 'A',
          atom: 'CD2',
          distance: 3.5,
          interactionType: 'pi_stacking',
        },
      ];

      const features = extractPharmacophoreFeatures(contacts);

      expect(features).toHaveLength(1);
      expect(features[0].type).toBe('aromatic');
    });

    it('should classify pi_stacking interactions as aromatic regardless of residue', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'XXX999',
          chain: 'Z',
          atom: 'C1',
          distance: 4.0,
          interactionType: 'pi_stacking',
        },
      ];

      const features = extractPharmacophoreFeatures(contacts);

      expect(features).toHaveLength(1);
      expect(features[0].type).toBe('aromatic');
    });
  });

  describe('Hydrophobic features', () => {
    it('should classify C atoms with hydrophobic interaction as hydrophobic', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'LEU44',
          chain: 'A',
          atom: 'CD1',
          distance: 3.7,
          interactionType: 'hydrophobic',
        },
      ];

      const features = extractPharmacophoreFeatures(contacts);

      expect(features).toHaveLength(1);
      expect(features[0].type).toBe('hydrophobic');
      expect(features[0].residue).toBe('LEU44');
    });

    it('should classify carbon contacts in non-aromatic residues as hydrophobic', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'VAL88',
          chain: 'A',
          atom: 'CG1',
          distance: 3.5,
          interactionType: 'hydrophobic',
        },
        {
          residue: 'ILE92',
          chain: 'A',
          atom: 'CD1',
          distance: 3.9,
          interactionType: 'hydrophobic',
        },
      ];

      const features = extractPharmacophoreFeatures(contacts);

      expect(features).toHaveLength(2);
      expect(features.every((f) => f.type === 'hydrophobic')).toBe(true);
    });
  });

  describe('Charged features (salt bridges)', () => {
    it('should classify ARG with salt_bridge as positive', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'ARG150',
          chain: 'A',
          atom: 'NH1',
          distance: 3.2,
          interactionType: 'salt_bridge',
        },
      ];

      const features = extractPharmacophoreFeatures(contacts);

      expect(features).toHaveLength(1);
      expect(features[0].type).toBe('positive');
    });

    it('should classify LYS with salt_bridge as positive', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'LYS55',
          chain: 'B',
          atom: 'NZ',
          distance: 3.0,
          interactionType: 'salt_bridge',
        },
      ];

      const features = extractPharmacophoreFeatures(contacts);

      expect(features).toHaveLength(1);
      expect(features[0].type).toBe('positive');
    });

    it('should classify HIS with salt_bridge as positive', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'HIS33',
          chain: 'A',
          atom: 'NE2',
          distance: 3.4,
          interactionType: 'salt_bridge',
        },
      ];

      const features = extractPharmacophoreFeatures(contacts);

      expect(features).toHaveLength(1);
      expect(features[0].type).toBe('positive');
    });

    it('should classify ASP with salt_bridge as negative', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'ASP120',
          chain: 'A',
          atom: 'OD1',
          distance: 3.1,
          interactionType: 'salt_bridge',
        },
      ];

      const features = extractPharmacophoreFeatures(contacts);

      expect(features).toHaveLength(1);
      expect(features[0].type).toBe('negative');
    });

    it('should classify GLU with salt_bridge as negative', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'GLU88',
          chain: 'C',
          atom: 'OE2',
          distance: 3.3,
          interactionType: 'salt_bridge',
        },
      ];

      const features = extractPharmacophoreFeatures(contacts);

      expect(features).toHaveLength(1);
      expect(features[0].type).toBe('negative');
    });
  });

  describe('Mixed contacts', () => {
    it('should correctly classify multiple different contact types', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'PHE50',
          chain: 'A',
          atom: 'CG',
          distance: 3.8,
          interactionType: 'pi_stacking',
        },
        {
          residue: 'SER60',
          chain: 'A',
          atom: 'OG',
          distance: 3.0,
          interactionType: 'hydrogen_bond',
        },
        {
          residue: 'LEU70',
          chain: 'A',
          atom: 'CD1',
          distance: 3.6,
          interactionType: 'hydrophobic',
        },
        {
          residue: 'ARG80',
          chain: 'A',
          atom: 'NH2',
          distance: 3.2,
          interactionType: 'salt_bridge',
        },
      ];

      const features = extractPharmacophoreFeatures(contacts);

      expect(features).toHaveLength(4);

      const types = features.map((f) => f.type);
      expect(types).toContain('aromatic');
      expect(types).toContain('acceptor');
      expect(types).toContain('hydrophobic');
      expect(types).toContain('positive');
    });
  });

  describe('Edge cases', () => {
    it('should return empty array for empty contacts', () => {
      const features = extractPharmacophoreFeatures([]);
      expect(features).toHaveLength(0);
    });

    it('should handle "other" interaction types appropriately', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'ALA10',
          chain: 'A',
          atom: 'CB',
          distance: 4.5,
          interactionType: 'other',
        },
      ];

      const features = extractPharmacophoreFeatures(contacts);

      // Should still produce a feature based on atom type (C = hydrophobic)
      expect(features).toHaveLength(1);
      expect(features[0].type).toBe('hydrophobic');
    });
  });
});
