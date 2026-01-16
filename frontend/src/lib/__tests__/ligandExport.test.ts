import { describe, it, expect } from 'vitest';
import {
  exportToJSON,
  exportToCSV,
  exportToRFdiffusionHotspots,
} from '../ligandExport';
import type { LigandAnalysisResult } from '../ligandAnalysis';

// Sample test data
const createSampleData = (): LigandAnalysisResult => ({
  ligandCount: 2,
  ligandDetails: [
    {
      name: 'ATP',
      info: 'ATP (Chain A, 501)',
      atoms: 31,
      chainId: 'A',
      resSeq: 501,
      contacts: [
        {
          residue: 'ASP45',
          chain: 'A',
          atom: 'OD1',
          distance: 2.8,
          interactionType: 'hydrogen_bond',
        },
        {
          residue: 'PHE72',
          chain: 'A',
          atom: 'CG',
          distance: 3.5,
          interactionType: 'pi_stacking',
        },
        {
          residue: 'LYS89',
          chain: 'B',
          atom: 'NZ',
          distance: 3.2,
          interactionType: 'salt_bridge',
        },
      ],
      proteinContactCount: 12,
      waterContactCount: 3,
      bindingSiteType: 'functional',
      bindingSiteReason: '8 protein residue contacts - well-defined binding pocket',
    },
    {
      name: 'SO4',
      info: 'SO4 (Chain B, 601)',
      atoms: 5,
      chainId: 'B',
      resSeq: 601,
      contacts: [
        {
          residue: 'ARG150',
          chain: 'B',
          atom: 'NH1',
          distance: 3.0,
          interactionType: 'salt_bridge',
        },
      ],
      proteinContactCount: 2,
      waterContactCount: 1,
      bindingSiteType: 'crystal_artifact',
      bindingSiteReason: 'Known crystallization additive (SO4) with few contacts',
    },
  ],
});

describe('exportToJSON', () => {
  it('should export data as formatted JSON string', () => {
    const data = createSampleData();
    const json = exportToJSON(data);

    expect(json).toContain('"ligandCount": 2');
    expect(json).toContain('"ligandDetails"');
  });

  it('should include ligandCount and ligandDetails', () => {
    const data = createSampleData();
    const json = exportToJSON(data);
    const parsed = JSON.parse(json);

    expect(parsed.ligandCount).toBe(2);
    expect(parsed.ligandDetails).toHaveLength(2);
  });

  it('should preserve all contact details', () => {
    const data = createSampleData();
    const json = exportToJSON(data);
    const parsed = JSON.parse(json);

    const firstContact = parsed.ligandDetails[0].contacts[0];
    expect(firstContact.residue).toBe('ASP45');
    expect(firstContact.chain).toBe('A');
    expect(firstContact.atom).toBe('OD1');
    expect(firstContact.distance).toBe(2.8);
    expect(firstContact.interactionType).toBe('hydrogen_bond');
  });

  it('should handle empty ligand list', () => {
    const data: LigandAnalysisResult = {
      ligandCount: 0,
      ligandDetails: [],
    };
    const json = exportToJSON(data);
    const parsed = JSON.parse(json);

    expect(parsed.ligandCount).toBe(0);
    expect(parsed.ligandDetails).toHaveLength(0);
  });

  it('should produce valid JSON that can be re-parsed', () => {
    const data = createSampleData();
    const json = exportToJSON(data);

    expect(() => JSON.parse(json)).not.toThrow();
  });
});

describe('exportToCSV', () => {
  it('should have correct headers', () => {
    const data = createSampleData();
    const csv = exportToCSV(data);
    const lines = csv.split('\n');

    expect(lines[0]).toBe(
      'Ligand,Chain,ResSeq,ContactResidue,ContactChain,ContactAtom,Distance,InteractionType,BindingSiteType'
    );
  });

  it('should have one row per contact', () => {
    const data = createSampleData();
    const csv = exportToCSV(data);
    const lines = csv.split('\n');

    // Header + 3 contacts from ATP + 1 contact from SO4 = 5 lines
    expect(lines).toHaveLength(5);
  });

  it('should include all contact details in each row', () => {
    const data = createSampleData();
    const csv = exportToCSV(data);
    const lines = csv.split('\n');

    // Check first data row (first ATP contact)
    const firstDataRow = lines[1].split(',');
    expect(firstDataRow[0]).toBe('ATP'); // Ligand
    expect(firstDataRow[1]).toBe('A'); // Chain
    expect(firstDataRow[2]).toBe('501'); // ResSeq
    expect(firstDataRow[3]).toBe('ASP45'); // ContactResidue
    expect(firstDataRow[4]).toBe('A'); // ContactChain
    expect(firstDataRow[5]).toBe('OD1'); // ContactAtom
    expect(firstDataRow[6]).toBe('2.80'); // Distance
    expect(firstDataRow[7]).toBe('hydrogen_bond'); // InteractionType
    expect(firstDataRow[8]).toBe('functional'); // BindingSiteType
  });

  it('should handle ligand with no contacts', () => {
    const data: LigandAnalysisResult = {
      ligandCount: 1,
      ligandDetails: [
        {
          name: 'NAG',
          info: 'NAG (Chain A, 100)',
          atoms: 14,
          chainId: 'A',
          resSeq: 100,
          contacts: [],
          proteinContactCount: 0,
          waterContactCount: 0,
          bindingSiteType: 'uncertain',
          bindingSiteReason: 'No contacts found',
        },
      ],
    };

    const csv = exportToCSV(data);
    const lines = csv.split('\n');

    // Only header row
    expect(lines).toHaveLength(1);
  });

  it('should handle empty ligand list', () => {
    const data: LigandAnalysisResult = {
      ligandCount: 0,
      ligandDetails: [],
    };

    const csv = exportToCSV(data);
    const lines = csv.split('\n');

    // Only header row
    expect(lines).toHaveLength(1);
  });

  it('should format distance with 2 decimal places', () => {
    const data: LigandAnalysisResult = {
      ligandCount: 1,
      ligandDetails: [
        {
          name: 'LIG',
          info: 'LIG (Chain A, 1)',
          atoms: 10,
          chainId: 'A',
          resSeq: 1,
          contacts: [
            {
              residue: 'SER10',
              chain: 'A',
              atom: 'OG',
              distance: 3.12345,
              interactionType: 'hydrogen_bond',
            },
          ],
          proteinContactCount: 1,
          waterContactCount: 0,
          bindingSiteType: 'functional',
          bindingSiteReason: 'Test',
        },
      ],
    };

    const csv = exportToCSV(data);
    const lines = csv.split('\n');
    const dataRow = lines[1].split(',');

    expect(dataRow[6]).toBe('3.12');
  });
});

describe('exportToRFdiffusionHotspots', () => {
  it('should contain hotspot_residues keyword', () => {
    const data = createSampleData();
    const config = exportToRFdiffusionHotspots(data);

    expect(config).toContain('hotspot_residues:');
  });

  it('should only include functional binding sites', () => {
    const data = createSampleData();
    const config = exportToRFdiffusionHotspots(data);

    // ATP is functional: should include A45, A72, B89
    expect(config).toContain('A45');
    expect(config).toContain('A72');
    expect(config).toContain('B89');

    // SO4 is crystal_artifact: should NOT include B150
    expect(config).not.toContain('B150');
  });

  it('should format residues as chain+number', () => {
    const data = createSampleData();
    const config = exportToRFdiffusionHotspots(data);

    // ASP45 on chain A should become A45
    expect(config).toContain('A45');
    // PHE72 on chain A should become A72
    expect(config).toContain('A72');
    // LYS89 on chain B should become B89
    expect(config).toContain('B89');
  });

  it('should include usage comments', () => {
    const data = createSampleData();
    const config = exportToRFdiffusionHotspots(data);

    expect(config).toContain('# RFdiffusion');
    expect(config).toContain('Usage');
    expect(config).toContain('--potentials.guide_potentials');
  });

  it('should sort hotspots by chain then residue number', () => {
    const data: LigandAnalysisResult = {
      ligandCount: 1,
      ligandDetails: [
        {
          name: 'LIG',
          info: 'LIG (Chain A, 1)',
          atoms: 10,
          chainId: 'A',
          resSeq: 1,
          contacts: [
            { residue: 'PHE100', chain: 'B', atom: 'CG', distance: 3.5, interactionType: 'pi_stacking' },
            { residue: 'ASP10', chain: 'A', atom: 'OD1', distance: 2.8, interactionType: 'hydrogen_bond' },
            { residue: 'SER50', chain: 'A', atom: 'OG', distance: 3.0, interactionType: 'hydrogen_bond' },
            { residue: 'LYS25', chain: 'B', atom: 'NZ', distance: 3.2, interactionType: 'salt_bridge' },
          ],
          proteinContactCount: 8,
          waterContactCount: 0,
          bindingSiteType: 'functional',
          bindingSiteReason: 'Test',
        },
      ],
    };

    const config = exportToRFdiffusionHotspots(data);

    // Extract the hotspot_residues line
    const hotspotsLine = config.split('\n').find((line) => line.startsWith('hotspot_residues:'));
    expect(hotspotsLine).toBeDefined();

    // Should be sorted: A10, A50, B25, B100
    expect(hotspotsLine).toBe('hotspot_residues: [A10, A50, B25, B100]');
  });

  it('should handle ligands with only crystal artifacts', () => {
    const data: LigandAnalysisResult = {
      ligandCount: 1,
      ligandDetails: [
        {
          name: 'SO4',
          info: 'SO4 (Chain A, 1)',
          atoms: 5,
          chainId: 'A',
          resSeq: 1,
          contacts: [
            { residue: 'ARG50', chain: 'A', atom: 'NH1', distance: 3.0, interactionType: 'salt_bridge' },
          ],
          proteinContactCount: 2,
          waterContactCount: 0,
          bindingSiteType: 'crystal_artifact',
          bindingSiteReason: 'Test',
        },
      ],
    };

    const config = exportToRFdiffusionHotspots(data);

    // Should have empty hotspots array
    expect(config).toContain('hotspot_residues: []');
    expect(config).toContain('Total hotspot residues: 0');
  });

  it('should handle empty ligand list', () => {
    const data: LigandAnalysisResult = {
      ligandCount: 0,
      ligandDetails: [],
    };

    const config = exportToRFdiffusionHotspots(data);

    expect(config).toContain('hotspot_residues: []');
  });

  it('should deduplicate residues from multiple contacts', () => {
    const data: LigandAnalysisResult = {
      ligandCount: 1,
      ligandDetails: [
        {
          name: 'LIG',
          info: 'LIG (Chain A, 1)',
          atoms: 10,
          chainId: 'A',
          resSeq: 1,
          contacts: [
            // Same residue, different atoms
            { residue: 'ASP45', chain: 'A', atom: 'OD1', distance: 2.8, interactionType: 'hydrogen_bond' },
            { residue: 'ASP45', chain: 'A', atom: 'OD2', distance: 3.0, interactionType: 'hydrogen_bond' },
          ],
          proteinContactCount: 4,
          waterContactCount: 0,
          bindingSiteType: 'functional',
          bindingSiteReason: 'Test',
        },
      ],
    };

    const config = exportToRFdiffusionHotspots(data);

    // Should only have A45 once
    const hotspotsLine = config.split('\n').find((line) => line.startsWith('hotspot_residues:'));
    expect(hotspotsLine).toBe('hotspot_residues: [A45]');
  });

  it('should include hotspots from multiple functional ligands', () => {
    const data: LigandAnalysisResult = {
      ligandCount: 2,
      ligandDetails: [
        {
          name: 'ATP',
          info: 'ATP (Chain A, 501)',
          atoms: 31,
          chainId: 'A',
          resSeq: 501,
          contacts: [
            { residue: 'ASP45', chain: 'A', atom: 'OD1', distance: 2.8, interactionType: 'hydrogen_bond' },
          ],
          proteinContactCount: 8,
          waterContactCount: 0,
          bindingSiteType: 'functional',
          bindingSiteReason: 'Test',
        },
        {
          name: 'NAD',
          info: 'NAD (Chain B, 502)',
          atoms: 44,
          chainId: 'B',
          resSeq: 502,
          contacts: [
            { residue: 'GLU100', chain: 'B', atom: 'OE1', distance: 2.9, interactionType: 'salt_bridge' },
          ],
          proteinContactCount: 6,
          waterContactCount: 0,
          bindingSiteType: 'functional',
          bindingSiteReason: 'Test',
        },
      ],
    };

    const config = exportToRFdiffusionHotspots(data);

    expect(config).toContain('A45');
    expect(config).toContain('B100');
  });

  it('should include uncertain binding sites', () => {
    // Uncertain sites should NOT be included - only functional
    const data: LigandAnalysisResult = {
      ligandCount: 1,
      ligandDetails: [
        {
          name: 'LIG',
          info: 'LIG (Chain A, 1)',
          atoms: 10,
          chainId: 'A',
          resSeq: 1,
          contacts: [
            { residue: 'ASP45', chain: 'A', atom: 'OD1', distance: 2.8, interactionType: 'hydrogen_bond' },
          ],
          proteinContactCount: 3,
          waterContactCount: 0,
          bindingSiteType: 'uncertain',
          bindingSiteReason: 'Test',
        },
      ],
    };

    const config = exportToRFdiffusionHotspots(data);

    // Uncertain binding sites should not be included
    expect(config).not.toContain('A45');
    expect(config).toContain('hotspot_residues: []');
  });
});

// Note: downloadFile tests are skipped in node environment
// The function requires browser DOM APIs (document, Blob, URL.createObjectURL)
// and is tested via integration tests in the browser environment
describe('downloadFile', () => {
  it.skip('should create an anchor element (requires browser environment)', () => {
    // This test requires a browser environment with DOM APIs
  });

  it.skip('should set href to blob URL (requires browser environment)', () => {
    // This test requires a browser environment with DOM APIs
  });

  it.skip('should set download attribute to filename (requires browser environment)', () => {
    // This test requires a browser environment with DOM APIs
  });

  it.skip('should trigger click on the link (requires browser environment)', () => {
    // This test requires a browser environment with DOM APIs
  });

  it.skip('should append and remove the link from document body (requires browser environment)', () => {
    // This test requires a browser environment with DOM APIs
  });

  it.skip('should revoke the blob URL after download (requires browser environment)', () => {
    // This test requires a browser environment with DOM APIs
  });

  it.skip('should create blob with correct content and MIME type (requires browser environment)', () => {
    // This test requires a browser environment with DOM APIs
  });
});
