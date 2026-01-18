/**
 * Test script for metal and ligand analysis functions
 * Run with: node test-analysis.mjs
 */

import fs from 'fs';
import path from 'path';

// Read the compiled analysis functions
// Since these are TypeScript, we'll inline a simplified version for testing

// Metal elements
const METAL_ELEMENTS = new Set([
  'LI', 'NA', 'K', 'RB', 'CS', 'FR',
  'BE', 'MG', 'CA', 'SR', 'BA', 'RA',
  'SC', 'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN',
  'Y', 'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD',
  'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG',
  'AL', 'GA', 'IN', 'SN', 'TL', 'PB', 'BI', 'PO',
  'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU',
  'AC', 'TH', 'PA', 'U', 'NP', 'PU', 'AM', 'CM', 'BK', 'CF', 'ES', 'FM', 'MD', 'NO', 'LR'
]);

function distance(a, b) {
  const dx = a[0] - b[0];
  const dy = a[1] - b[1];
  const dz = a[2] - b[2];
  return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

function findMetalCoordinationFromPDB(pdbContent, coordinationRadius = 3.0) {
  const lines = pdbContent.split('\n');
  const atoms = [];

  for (const line of lines) {
    if (line.startsWith('ATOM') || line.startsWith('HETATM')) {
      const isHETATM = line.startsWith('HETATM');
      try {
        atoms.push({
          serial: parseInt(line.slice(6, 11).trim()),
          name: line.slice(12, 16).trim(),
          resName: line.slice(17, 20).trim(),
          chainId: line.slice(21, 22).trim(),
          resSeq: parseInt(line.slice(22, 26).trim()),
          x: parseFloat(line.slice(30, 38).trim()),
          y: parseFloat(line.slice(38, 46).trim()),
          z: parseFloat(line.slice(46, 54).trim()),
          element: line.slice(76, 78).trim().toUpperCase() || line.slice(12, 14).trim().toUpperCase(),
          isHETATM
        });
      } catch (e) {
        // Skip malformed lines
      }
    }
  }

  // Find metal ions
  const metals = atoms.filter(a => METAL_ELEMENTS.has(a.element));
  const results = [];

  for (const metal of metals) {
    const metalPos = [metal.x, metal.y, metal.z];
    const coordinating = [];

    // Find coordinating atoms
    for (const atom of atoms) {
      if (atom.serial === metal.serial) continue;

      const atomPos = [atom.x, atom.y, atom.z];
      const dist = distance(metalPos, atomPos);

      if (dist <= coordinationRadius) {
        const isWater = atom.resName === 'HOH' || atom.resName === 'WAT';
        coordinating.push({
          atom: atom.name,
          residue: `${atom.resName}${atom.resSeq}`,
          chain: atom.chainId,
          distance: dist,
          isWater,
          position: atomPos
        });
      }
    }

    // Sort by distance
    coordinating.sort((a, b) => a.distance - b.distance);

    const proteinContacts = coordinating.filter(c => !c.isWater).length;
    const waterContacts = coordinating.filter(c => c.isWater).length;

    results.push({
      element: metal.element,
      info: `${metal.element} (${metal.resName}${metal.resSeq}, Chain ${metal.chainId})`,
      chainId: metal.chainId,
      resSeq: metal.resSeq,
      resName: metal.resName,
      pos: metalPos,
      coordinating,
      coordinationNumber: coordinating.length,
      proteinContacts,
      waterContacts
    });
  }

  return results;
}

// Test with calmodulin
console.log('=== Testing Metal Analysis with Calmodulin (1CLL) ===\n');

// Try to read from WSL path via Windows
const wslPath = '\\\\wsl$\\Ubuntu\\tmp\\1cll.pdb';
let pdbContent;

try {
  pdbContent = fs.readFileSync(wslPath, 'utf8');
  console.log(`Loaded PDB: ${pdbContent.length} bytes\n`);
} catch (e) {
  console.log('Could not read from WSL path, trying alternate...');
  // If running from WSL, try direct path
  try {
    pdbContent = fs.readFileSync('/tmp/1cll.pdb', 'utf8');
    console.log(`Loaded PDB: ${pdbContent.length} bytes\n`);
  } catch (e2) {
    console.error('Could not load PDB file:', e2.message);
    process.exit(1);
  }
}

// Run analysis
const metals = findMetalCoordinationFromPDB(pdbContent, 3.0);

console.log(`Found ${metals.length} metal ions:\n`);

for (const metal of metals) {
  console.log(`--- ${metal.info} ---`);
  console.log(`  Position: (${metal.pos[0].toFixed(2)}, ${metal.pos[1].toFixed(2)}, ${metal.pos[2].toFixed(2)})`);
  console.log(`  Coordination number: ${metal.coordinationNumber}`);
  console.log(`  Protein contacts: ${metal.proteinContacts}`);
  console.log(`  Water contacts: ${metal.waterContacts}`);
  console.log(`  Coordinating atoms:`);

  for (const coord of metal.coordinating.slice(0, 8)) {
    const waterNote = coord.isWater ? ' (water)' : '';
    console.log(`    ${coord.residue} ${coord.atom}: ${coord.distance.toFixed(2)} A${waterNote}`);
  }
  if (metal.coordinating.length > 8) {
    console.log(`    ... and ${metal.coordinating.length - 8} more`);
  }
  console.log('');
}

// Test with TB-substituted calmodulin
console.log('\n=== Testing with TB-substituted Calmodulin ===\n');

try {
  const tbPdb = fs.readFileSync(wslPath.replace('1cll.pdb', '1cll_tb.pdb'), 'utf8');
  const tbMetals = findMetalCoordinationFromPDB(tbPdb, 3.0);

  console.log(`Found ${tbMetals.length} metal ions:\n`);

  for (const metal of tbMetals) {
    console.log(`${metal.info}: CN=${metal.coordinationNumber}, Protein=${metal.proteinContacts}, Water=${metal.waterContacts}`);
  }
} catch (e) {
  console.log('Could not test TB-substituted structure');
}

// Test with RFdiffusion output
console.log('\n=== Testing with RFdiffusion Output ===\n');

try {
  const outputPdb = fs.readFileSync(wslPath.replace('1cll.pdb', 'rfd3_output.pdb'), 'utf8');
  const outputMetals = findMetalCoordinationFromPDB(outputPdb, 3.0);

  console.log(`Found ${outputMetals.length} metal ions:\n`);

  for (const metal of outputMetals) {
    console.log(`--- ${metal.info} ---`);
    console.log(`  Coordination number: ${metal.coordinationNumber}`);
    console.log(`  Protein contacts: ${metal.proteinContacts}`);
    console.log(`  Water contacts: ${metal.waterContacts}`);
    console.log(`  Top coordinating atoms:`);
    for (const coord of metal.coordinating.slice(0, 6)) {
      console.log(`    ${coord.residue} ${coord.atom}: ${coord.distance.toFixed(2)} A`);
    }
  }
} catch (e) {
  console.log('Could not test RFdiffusion output:', e.message);
}

console.log('\n=== Test Complete ===');
