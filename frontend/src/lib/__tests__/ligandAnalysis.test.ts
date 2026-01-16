import { describe, it, expect } from 'vitest';
import {
  extractPharmacophoreFeatures,
  classifyShellResidues,
  type LigandContact,
  type PharmacophoreFeature,
  type PharmacophoreType,
  type ShellAnalysis,
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

describe('classifyShellResidues', () => {
  describe('Primary shell classification', () => {
    it('should classify residues within 4A as primary shell', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'SER10',
          chain: 'A',
          atom: 'OG',
          distance: 2.5,
          interactionType: 'hydrogen_bond',
        },
        {
          residue: 'LEU20',
          chain: 'A',
          atom: 'CD1',
          distance: 3.5,
          interactionType: 'hydrophobic',
        },
        {
          residue: 'PHE30',
          chain: 'A',
          atom: 'CG',
          distance: 3.9,
          interactionType: 'pi_stacking',
        },
      ];

      const result = classifyShellResidues(contacts);

      expect(result.primary).toHaveLength(3);
      expect(result.secondary).toHaveLength(0);
      expect(result.stats.primaryCount).toBe(3);
      expect(result.stats.secondaryCount).toBe(0);
    });

    it('should not include 4A exactly in primary shell', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'ALA15',
          chain: 'A',
          atom: 'CB',
          distance: 4.0,
          interactionType: 'hydrophobic',
        },
      ];

      const result = classifyShellResidues(contacts);

      expect(result.primary).toHaveLength(0);
      expect(result.secondary).toHaveLength(1);
    });
  });

  describe('Secondary shell classification', () => {
    it('should classify residues 4-6A as secondary shell', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'VAL40',
          chain: 'A',
          atom: 'CG1',
          distance: 4.0,
          interactionType: 'hydrophobic',
        },
        {
          residue: 'ILE50',
          chain: 'A',
          atom: 'CD1',
          distance: 5.0,
          interactionType: 'other',
        },
        {
          residue: 'MET60',
          chain: 'A',
          atom: 'CE',
          distance: 6.0,
          interactionType: 'other',
        },
      ];

      const result = classifyShellResidues(contacts);

      expect(result.primary).toHaveLength(0);
      expect(result.secondary).toHaveLength(3);
      expect(result.stats.secondaryCount).toBe(3);
    });

    it('should exclude residues beyond 6A', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'GLY70',
          chain: 'A',
          atom: 'CA',
          distance: 6.1,
          interactionType: 'other',
        },
        {
          residue: 'ALA80',
          chain: 'A',
          atom: 'CB',
          distance: 8.0,
          interactionType: 'other',
        },
      ];

      const result = classifyShellResidues(contacts);

      expect(result.primary).toHaveLength(0);
      expect(result.secondary).toHaveLength(0);
      expect(result.stats.totalCount).toBe(0);
    });
  });

  describe('Shell statistics', () => {
    it('should calculate average distances correctly', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'SER10',
          chain: 'A',
          atom: 'OG',
          distance: 2.0,
          interactionType: 'hydrogen_bond',
        },
        {
          residue: 'LEU20',
          chain: 'A',
          atom: 'CD1',
          distance: 3.0,
          interactionType: 'hydrophobic',
        },
        {
          residue: 'VAL40',
          chain: 'A',
          atom: 'CG1',
          distance: 4.5,
          interactionType: 'other',
        },
        {
          residue: 'ILE50',
          chain: 'A',
          atom: 'CD1',
          distance: 5.5,
          interactionType: 'other',
        },
      ];

      const result = classifyShellResidues(contacts);

      expect(result.stats.primaryCount).toBe(2);
      expect(result.stats.secondaryCount).toBe(2);
      expect(result.stats.totalCount).toBe(4);
      expect(result.stats.avgPrimaryDistance).toBe(2.5); // (2.0 + 3.0) / 2
      expect(result.stats.avgSecondaryDistance).toBe(5.0); // (4.5 + 5.5) / 2
    });

    it('should return correct total count', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'SER10',
          chain: 'A',
          atom: 'OG',
          distance: 3.0,
          interactionType: 'hydrogen_bond',
        },
        {
          residue: 'VAL40',
          chain: 'A',
          atom: 'CG1',
          distance: 5.0,
          interactionType: 'other',
        },
        {
          residue: 'GLY70',
          chain: 'A',
          atom: 'CA',
          distance: 7.0,
          interactionType: 'other',
        },
      ];

      const result = classifyShellResidues(contacts);

      // Only 2 contacts within cutoff (1 primary, 1 secondary)
      expect(result.stats.totalCount).toBe(2);
    });
  });

  describe('Edge cases', () => {
    it('should return empty arrays for empty contacts', () => {
      const result = classifyShellResidues([]);

      expect(result.primary).toHaveLength(0);
      expect(result.secondary).toHaveLength(0);
      expect(result.stats.primaryCount).toBe(0);
      expect(result.stats.secondaryCount).toBe(0);
      expect(result.stats.totalCount).toBe(0);
      expect(result.stats.avgPrimaryDistance).toBe(0);
      expect(result.stats.avgSecondaryDistance).toBe(0);
    });

    it('should return 0 for avgPrimaryDistance when primary shell is empty', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'VAL40',
          chain: 'A',
          atom: 'CG1',
          distance: 5.0,
          interactionType: 'other',
        },
      ];

      const result = classifyShellResidues(contacts);

      expect(result.primary).toHaveLength(0);
      expect(result.stats.avgPrimaryDistance).toBe(0);
      expect(result.stats.avgSecondaryDistance).toBe(5.0);
    });

    it('should return 0 for avgSecondaryDistance when secondary shell is empty', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'SER10',
          chain: 'A',
          atom: 'OG',
          distance: 3.0,
          interactionType: 'hydrogen_bond',
        },
      ];

      const result = classifyShellResidues(contacts);

      expect(result.secondary).toHaveLength(0);
      expect(result.stats.avgSecondaryDistance).toBe(0);
      expect(result.stats.avgPrimaryDistance).toBe(3.0);
    });

    it('should handle all contacts in one shell', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'SER10',
          chain: 'A',
          atom: 'OG',
          distance: 2.0,
          interactionType: 'hydrogen_bond',
        },
        {
          residue: 'LEU20',
          chain: 'A',
          atom: 'CD1',
          distance: 3.0,
          interactionType: 'hydrophobic',
        },
        {
          residue: 'PHE30',
          chain: 'A',
          atom: 'CG',
          distance: 3.5,
          interactionType: 'pi_stacking',
        },
      ];

      const result = classifyShellResidues(contacts);

      expect(result.primary).toHaveLength(3);
      expect(result.secondary).toHaveLength(0);
      expect(result.stats.totalCount).toBe(3);
      expect(result.stats.avgSecondaryDistance).toBe(0);
    });

    it('should preserve contact data in classified arrays', () => {
      const contacts: LigandContact[] = [
        {
          residue: 'ARG100',
          chain: 'B',
          atom: 'NH1',
          distance: 3.2,
          interactionType: 'salt_bridge',
        },
      ];

      const result = classifyShellResidues(contacts);

      expect(result.primary).toHaveLength(1);
      expect(result.primary[0].residue).toBe('ARG100');
      expect(result.primary[0].chain).toBe('B');
      expect(result.primary[0].atom).toBe('NH1');
      expect(result.primary[0].distance).toBe(3.2);
      expect(result.primary[0].interactionType).toBe('salt_bridge');
    });
  });
});
