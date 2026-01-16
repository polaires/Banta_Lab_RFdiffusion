import { describe, it, expect } from 'vitest';
import {
  extractPharmacophoreFeatures,
  classifyShellResidues,
  calculateBinderQualityScore,
  type LigandContact,
  type PharmacophoreFeature,
  type PharmacophoreType,
  type ShellAnalysis,
  type BinderQualityResult,
  type QualityRating,
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

describe('calculateBinderQualityScore', () => {
  describe('High quality binding site (EXCELLENT)', () => {
    it('should return EXCELLENT rating for well-formed binding site with diverse interactions', () => {
      const contacts: LigandContact[] = [
        // 3 H-bonds (30 points)
        { residue: 'SER10', chain: 'A', atom: 'OG', distance: 2.8, interactionType: 'hydrogen_bond' },
        { residue: 'ASN20', chain: 'A', atom: 'ND2', distance: 3.0, interactionType: 'hydrogen_bond' },
        { residue: 'THR30', chain: 'A', atom: 'OG1', distance: 2.9, interactionType: 'hydrogen_bond' },
        // 5 hydrophobic contacts (25 points)
        { residue: 'LEU40', chain: 'A', atom: 'CD1', distance: 3.5, interactionType: 'hydrophobic' },
        { residue: 'ILE50', chain: 'A', atom: 'CD1', distance: 3.6, interactionType: 'hydrophobic' },
        { residue: 'VAL60', chain: 'A', atom: 'CG1', distance: 3.7, interactionType: 'hydrophobic' },
        { residue: 'ALA70', chain: 'A', atom: 'CB', distance: 3.8, interactionType: 'hydrophobic' },
        { residue: 'MET80', chain: 'A', atom: 'CE', distance: 3.9, interactionType: 'hydrophobic' },
        // 2 pi-stacking (15 points)
        { residue: 'PHE90', chain: 'A', atom: 'CG', distance: 3.8, interactionType: 'pi_stacking' },
        { residue: 'TYR100', chain: 'A', atom: 'CE1', distance: 3.9, interactionType: 'pi_stacking' },
        // 2 salt bridges (15 points)
        { residue: 'ARG110', chain: 'A', atom: 'NH1', distance: 3.2, interactionType: 'salt_bridge' },
        { residue: 'GLU120', chain: 'A', atom: 'OE1', distance: 3.3, interactionType: 'salt_bridge' },
      ];

      const result = calculateBinderQualityScore(contacts);

      expect(result.score).toBe(100); // Max score
      expect(result.rating).toBe('EXCELLENT');
      expect(result.breakdown.hbondScore).toBe(30);
      expect(result.breakdown.hydrophobicScore).toBe(25);
      expect(result.breakdown.piStackScore).toBe(15);
      expect(result.breakdown.saltBridgeScore).toBe(15);
      expect(result.breakdown.burialnessScore).toBe(15);
      expect(result.suggestions).toHaveLength(0);
    });

    it('should cap individual scores at their maximums', () => {
      const contacts: LigandContact[] = [
        // 5 H-bonds (should cap at 30)
        { residue: 'SER1', chain: 'A', atom: 'OG', distance: 2.8, interactionType: 'hydrogen_bond' },
        { residue: 'SER2', chain: 'A', atom: 'OG', distance: 2.9, interactionType: 'hydrogen_bond' },
        { residue: 'SER3', chain: 'A', atom: 'OG', distance: 2.8, interactionType: 'hydrogen_bond' },
        { residue: 'SER4', chain: 'A', atom: 'OG', distance: 3.0, interactionType: 'hydrogen_bond' },
        { residue: 'SER5', chain: 'A', atom: 'OG', distance: 2.7, interactionType: 'hydrogen_bond' },
        // 10 hydrophobic (should cap at 25)
        { residue: 'LEU1', chain: 'A', atom: 'CD1', distance: 3.5, interactionType: 'hydrophobic' },
        { residue: 'LEU2', chain: 'A', atom: 'CD1', distance: 3.6, interactionType: 'hydrophobic' },
        { residue: 'LEU3', chain: 'A', atom: 'CD1', distance: 3.7, interactionType: 'hydrophobic' },
        { residue: 'LEU4', chain: 'A', atom: 'CD1', distance: 3.8, interactionType: 'hydrophobic' },
        { residue: 'LEU5', chain: 'A', atom: 'CD1', distance: 3.9, interactionType: 'hydrophobic' },
        { residue: 'LEU6', chain: 'A', atom: 'CD1', distance: 3.5, interactionType: 'hydrophobic' },
        { residue: 'LEU7', chain: 'A', atom: 'CD1', distance: 3.6, interactionType: 'hydrophobic' },
        { residue: 'LEU8', chain: 'A', atom: 'CD1', distance: 3.7, interactionType: 'hydrophobic' },
        { residue: 'LEU9', chain: 'A', atom: 'CD1', distance: 3.8, interactionType: 'hydrophobic' },
        { residue: 'LEU10', chain: 'A', atom: 'CD1', distance: 3.9, interactionType: 'hydrophobic' },
      ];

      const result = calculateBinderQualityScore(contacts);

      expect(result.breakdown.hbondScore).toBe(30);
      expect(result.breakdown.hydrophobicScore).toBe(25);
      expect(result.breakdown.burialnessScore).toBe(15); // 15 contacts * 1.5 = 22.5, capped at 15
      expect(result.score).toBe(70); // 30 + 25 + 0 + 0 + 15
      expect(result.rating).toBe('EXCELLENT');
    });
  });

  describe('Low quality binding site (POOR)', () => {
    it('should return POOR rating for sparse contacts', () => {
      const contacts: LigandContact[] = [
        { residue: 'ALA10', chain: 'A', atom: 'CB', distance: 4.5, interactionType: 'other' },
      ];

      const result = calculateBinderQualityScore(contacts);

      expect(result.score).toBeLessThan(30);
      expect(result.rating).toBe('POOR');
      expect(result.suggestions.length).toBeGreaterThan(0);
    });

    it('should return zero scores for empty contacts', () => {
      const result = calculateBinderQualityScore([]);

      expect(result.score).toBe(0);
      expect(result.rating).toBe('POOR');
      expect(result.breakdown.hbondScore).toBe(0);
      expect(result.breakdown.hydrophobicScore).toBe(0);
      expect(result.breakdown.piStackScore).toBe(0);
      expect(result.breakdown.saltBridgeScore).toBe(0);
      expect(result.breakdown.burialnessScore).toBe(0);
    });
  });

  describe('Breakdown values', () => {
    it('should correctly calculate individual H-bond scores', () => {
      const contacts: LigandContact[] = [
        { residue: 'SER10', chain: 'A', atom: 'OG', distance: 2.8, interactionType: 'hydrogen_bond' },
        { residue: 'SER20', chain: 'A', atom: 'OG', distance: 2.9, interactionType: 'hydrogen_bond' },
      ];

      const result = calculateBinderQualityScore(contacts);

      expect(result.breakdown.hbondScore).toBe(20); // 2 * 10 = 20
      expect(result.breakdown.hydrophobicScore).toBe(0);
      expect(result.breakdown.piStackScore).toBe(0);
      expect(result.breakdown.saltBridgeScore).toBe(0);
      expect(result.breakdown.burialnessScore).toBe(3); // 2 * 1.5 = 3
    });

    it('should correctly calculate individual hydrophobic scores', () => {
      const contacts: LigandContact[] = [
        { residue: 'LEU10', chain: 'A', atom: 'CD1', distance: 3.5, interactionType: 'hydrophobic' },
        { residue: 'LEU20', chain: 'A', atom: 'CD1', distance: 3.6, interactionType: 'hydrophobic' },
        { residue: 'LEU30', chain: 'A', atom: 'CD1', distance: 3.7, interactionType: 'hydrophobic' },
      ];

      const result = calculateBinderQualityScore(contacts);

      expect(result.breakdown.hydrophobicScore).toBe(15); // 3 * 5 = 15
      expect(result.breakdown.burialnessScore).toBe(4.5); // 3 * 1.5 = 4.5
    });

    it('should correctly calculate pi-stacking scores with 7.5 multiplier', () => {
      const contacts: LigandContact[] = [
        { residue: 'PHE10', chain: 'A', atom: 'CG', distance: 3.8, interactionType: 'pi_stacking' },
      ];

      const result = calculateBinderQualityScore(contacts);

      expect(result.breakdown.piStackScore).toBe(7.5); // 1 * 7.5 = 7.5
    });

    it('should correctly calculate salt bridge scores with 7.5 multiplier', () => {
      const contacts: LigandContact[] = [
        { residue: 'ARG10', chain: 'A', atom: 'NH1', distance: 3.2, interactionType: 'salt_bridge' },
      ];

      const result = calculateBinderQualityScore(contacts);

      expect(result.breakdown.saltBridgeScore).toBe(7.5); // 1 * 7.5 = 7.5
    });

    it('should count "other" interaction type toward burialness only', () => {
      const contacts: LigandContact[] = [
        { residue: 'ALA10', chain: 'A', atom: 'CB', distance: 4.0, interactionType: 'other' },
        { residue: 'ALA20', chain: 'A', atom: 'CB', distance: 4.1, interactionType: 'other' },
        { residue: 'ALA30', chain: 'A', atom: 'CB', distance: 4.2, interactionType: 'other' },
        { residue: 'ALA40', chain: 'A', atom: 'CB', distance: 4.3, interactionType: 'other' },
      ];

      const result = calculateBinderQualityScore(contacts);

      expect(result.breakdown.hbondScore).toBe(0);
      expect(result.breakdown.hydrophobicScore).toBe(0);
      expect(result.breakdown.piStackScore).toBe(0);
      expect(result.breakdown.saltBridgeScore).toBe(0);
      expect(result.breakdown.burialnessScore).toBe(6); // 4 * 1.5 = 6
      expect(result.score).toBe(6);
    });
  });

  describe('Rating thresholds', () => {
    it('should return EXCELLENT for score >= 70', () => {
      // Create contacts for exactly 70 points
      const contacts: LigandContact[] = [
        // 3 H-bonds = 30 points
        { residue: 'SER1', chain: 'A', atom: 'OG', distance: 2.8, interactionType: 'hydrogen_bond' },
        { residue: 'SER2', chain: 'A', atom: 'OG', distance: 2.9, interactionType: 'hydrogen_bond' },
        { residue: 'SER3', chain: 'A', atom: 'OG', distance: 3.0, interactionType: 'hydrogen_bond' },
        // 5 hydrophobic = 25 points
        { residue: 'LEU1', chain: 'A', atom: 'CD1', distance: 3.5, interactionType: 'hydrophobic' },
        { residue: 'LEU2', chain: 'A', atom: 'CD1', distance: 3.6, interactionType: 'hydrophobic' },
        { residue: 'LEU3', chain: 'A', atom: 'CD1', distance: 3.7, interactionType: 'hydrophobic' },
        { residue: 'LEU4', chain: 'A', atom: 'CD1', distance: 3.8, interactionType: 'hydrophobic' },
        { residue: 'LEU5', chain: 'A', atom: 'CD1', distance: 3.9, interactionType: 'hydrophobic' },
        // 2 more for burialness to reach 10 total = 15 points
        { residue: 'GLY1', chain: 'A', atom: 'CA', distance: 4.0, interactionType: 'other' },
        { residue: 'GLY2', chain: 'A', atom: 'CA', distance: 4.1, interactionType: 'other' },
      ];

      const result = calculateBinderQualityScore(contacts);

      expect(result.score).toBe(70); // 30 + 25 + 0 + 0 + 15 = 70
      expect(result.rating).toBe('EXCELLENT');
    });

    it('should return GOOD for score >= 50 and < 70', () => {
      const contacts: LigandContact[] = [
        // 2 H-bonds = 20 points
        { residue: 'SER1', chain: 'A', atom: 'OG', distance: 2.8, interactionType: 'hydrogen_bond' },
        { residue: 'SER2', chain: 'A', atom: 'OG', distance: 2.9, interactionType: 'hydrogen_bond' },
        // 4 hydrophobic = 20 points
        { residue: 'LEU1', chain: 'A', atom: 'CD1', distance: 3.5, interactionType: 'hydrophobic' },
        { residue: 'LEU2', chain: 'A', atom: 'CD1', distance: 3.6, interactionType: 'hydrophobic' },
        { residue: 'LEU3', chain: 'A', atom: 'CD1', distance: 3.7, interactionType: 'hydrophobic' },
        { residue: 'LEU4', chain: 'A', atom: 'CD1', distance: 3.8, interactionType: 'hydrophobic' },
        // 4 more for burialness = 10 total
        { residue: 'GLY1', chain: 'A', atom: 'CA', distance: 4.0, interactionType: 'other' },
        { residue: 'GLY2', chain: 'A', atom: 'CA', distance: 4.1, interactionType: 'other' },
        { residue: 'GLY3', chain: 'A', atom: 'CA', distance: 4.2, interactionType: 'other' },
        { residue: 'GLY4', chain: 'A', atom: 'CA', distance: 4.3, interactionType: 'other' },
      ];

      const result = calculateBinderQualityScore(contacts);

      expect(result.score).toBe(55); // 20 + 20 + 0 + 0 + 15 = 55
      expect(result.rating).toBe('GOOD');
    });

    it('should return MODERATE for score >= 30 and < 50', () => {
      const contacts: LigandContact[] = [
        // 1 H-bond = 10 points
        { residue: 'SER1', chain: 'A', atom: 'OG', distance: 2.8, interactionType: 'hydrogen_bond' },
        // 2 hydrophobic = 10 points
        { residue: 'LEU1', chain: 'A', atom: 'CD1', distance: 3.5, interactionType: 'hydrophobic' },
        { residue: 'LEU2', chain: 'A', atom: 'CD1', distance: 3.6, interactionType: 'hydrophobic' },
        // 7 more for burialness = 10 total = 15 points (but only 10 total contacts)
        { residue: 'GLY1', chain: 'A', atom: 'CA', distance: 4.0, interactionType: 'other' },
        { residue: 'GLY2', chain: 'A', atom: 'CA', distance: 4.1, interactionType: 'other' },
        { residue: 'GLY3', chain: 'A', atom: 'CA', distance: 4.2, interactionType: 'other' },
        { residue: 'GLY4', chain: 'A', atom: 'CA', distance: 4.3, interactionType: 'other' },
        { residue: 'GLY5', chain: 'A', atom: 'CA', distance: 4.4, interactionType: 'other' },
        { residue: 'GLY6', chain: 'A', atom: 'CA', distance: 4.5, interactionType: 'other' },
        { residue: 'GLY7', chain: 'A', atom: 'CA', distance: 4.6, interactionType: 'other' },
      ];

      const result = calculateBinderQualityScore(contacts);

      expect(result.score).toBe(35); // 10 + 10 + 0 + 0 + 15 = 35
      expect(result.rating).toBe('MODERATE');
    });

    it('should return POOR for score < 30', () => {
      const contacts: LigandContact[] = [
        // 1 H-bond = 10 points
        { residue: 'SER1', chain: 'A', atom: 'OG', distance: 2.8, interactionType: 'hydrogen_bond' },
        // 2 hydrophobic = 10 points
        { residue: 'LEU1', chain: 'A', atom: 'CD1', distance: 3.5, interactionType: 'hydrophobic' },
        { residue: 'LEU2', chain: 'A', atom: 'CD1', distance: 3.6, interactionType: 'hydrophobic' },
      ];

      const result = calculateBinderQualityScore(contacts);

      expect(result.score).toBe(24.5); // 10 + 10 + 0 + 0 + 4.5 = 24.5
      expect(result.rating).toBe('POOR');
    });
  });

  describe('Suggestions', () => {
    it('should suggest adding H-bonds when count is low', () => {
      const contacts: LigandContact[] = [
        { residue: 'SER1', chain: 'A', atom: 'OG', distance: 2.8, interactionType: 'hydrogen_bond' },
      ];

      const result = calculateBinderQualityScore(contacts);

      expect(result.suggestions.some((s) => s.includes('H-bond'))).toBe(true);
    });

    it('should suggest increasing hydrophobic contacts when count is low', () => {
      const contacts: LigandContact[] = [
        { residue: 'LEU1', chain: 'A', atom: 'CD1', distance: 3.5, interactionType: 'hydrophobic' },
      ];

      const result = calculateBinderQualityScore(contacts);

      expect(result.suggestions.some((s) => s.includes('hydrophobic'))).toBe(true);
    });

    it('should suggest adding aromatic residues when no pi-stacking', () => {
      const contacts: LigandContact[] = [
        { residue: 'SER1', chain: 'A', atom: 'OG', distance: 2.8, interactionType: 'hydrogen_bond' },
      ];

      const result = calculateBinderQualityScore(contacts);

      expect(result.suggestions.some((s) => s.includes('aromatic'))).toBe(true);
    });

    it('should suggest adding charged residues when no salt bridges', () => {
      const contacts: LigandContact[] = [
        { residue: 'SER1', chain: 'A', atom: 'OG', distance: 2.8, interactionType: 'hydrogen_bond' },
      ];

      const result = calculateBinderQualityScore(contacts);

      expect(result.suggestions.some((s) => s.includes('charged') || s.includes('salt bridge'))).toBe(true);
    });

    it('should suggest deeper pocket when total contacts are low', () => {
      const contacts: LigandContact[] = [
        { residue: 'SER1', chain: 'A', atom: 'OG', distance: 2.8, interactionType: 'hydrogen_bond' },
      ];

      const result = calculateBinderQualityScore(contacts);

      expect(result.suggestions.some((s) => s.includes('shallow') || s.includes('10+'))).toBe(true);
    });

    it('should not generate suggestions for perfect binding site', () => {
      const contacts: LigandContact[] = [
        // 3 H-bonds
        { residue: 'SER1', chain: 'A', atom: 'OG', distance: 2.8, interactionType: 'hydrogen_bond' },
        { residue: 'SER2', chain: 'A', atom: 'OG', distance: 2.9, interactionType: 'hydrogen_bond' },
        { residue: 'SER3', chain: 'A', atom: 'OG', distance: 3.0, interactionType: 'hydrogen_bond' },
        // 5 hydrophobic
        { residue: 'LEU1', chain: 'A', atom: 'CD1', distance: 3.5, interactionType: 'hydrophobic' },
        { residue: 'LEU2', chain: 'A', atom: 'CD1', distance: 3.6, interactionType: 'hydrophobic' },
        { residue: 'LEU3', chain: 'A', atom: 'CD1', distance: 3.7, interactionType: 'hydrophobic' },
        { residue: 'LEU4', chain: 'A', atom: 'CD1', distance: 3.8, interactionType: 'hydrophobic' },
        { residue: 'LEU5', chain: 'A', atom: 'CD1', distance: 3.9, interactionType: 'hydrophobic' },
        // 2 pi-stacking
        { residue: 'PHE1', chain: 'A', atom: 'CG', distance: 3.8, interactionType: 'pi_stacking' },
        { residue: 'PHE2', chain: 'A', atom: 'CG', distance: 3.9, interactionType: 'pi_stacking' },
        // 1 salt bridge
        { residue: 'ARG1', chain: 'A', atom: 'NH1', distance: 3.2, interactionType: 'salt_bridge' },
      ];

      const result = calculateBinderQualityScore(contacts);

      expect(result.suggestions).toHaveLength(0);
    });
  });
});
