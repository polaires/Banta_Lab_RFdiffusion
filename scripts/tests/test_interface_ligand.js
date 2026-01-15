/**
 * Test: Symmetric Ligand-at-Interface Design
 *
 * This tests the NEW interface_ligand parameter which:
 * 1. Orients the ligand along the C2 symmetry axis (N=N bond along Z)
 * 2. Combines symmetry + ligand + RASA conditioning
 * 3. Validates that both chains contact the ligand
 */

const http = require('http');
const fs = require('fs');

const AZOBENZENE_SMILES = 'c1ccc(cc1)N=Nc2ccccc2';

async function submitJob(input) {
  return new Promise((resolve, reject) => {
    const postData = JSON.stringify({ input });
    const options = {
      hostname: 'localhost',
      port: 8000,
      path: '/runsync',
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Content-Length': Buffer.byteLength(postData)
      },
      timeout: 600000  // 10 minute timeout for long jobs
    };
    const req = http.request(options, (res) => {
      let data = '';
      res.on('data', chunk => data += chunk);
      res.on('end', () => {
        try { resolve(JSON.parse(data)); }
        catch (e) { resolve({ error: data }); }
      });
    });
    req.on('error', reject);
    req.on('timeout', () => {
      req.destroy();
      reject(new Error('Request timeout'));
    });
    req.write(postData);
    req.end();
  });
}

/**
 * Parse structure and count contacts between ligand and each chain
 */
function analyzeInterfaceContacts(pdbContent) {
  const lines = pdbContent.split('\n');
  const atomLines = lines.filter(l => l.startsWith('ATOM'));
  const hetatmLines = lines.filter(l => l.startsWith('HETATM'));

  // Get chain info
  const chainAtoms = {};
  for (const line of atomLines) {
    const chain = line[21];
    const x = parseFloat(line.substring(30, 38));
    const y = parseFloat(line.substring(38, 46));
    const z = parseFloat(line.substring(46, 54));
    if (!chainAtoms[chain]) chainAtoms[chain] = [];
    chainAtoms[chain].push({x, y, z});
  }

  // Get ligand atoms
  const ligandAtoms = [];
  for (const line of hetatmLines) {
    const x = parseFloat(line.substring(30, 38));
    const y = parseFloat(line.substring(38, 46));
    const z = parseFloat(line.substring(46, 54));
    ligandAtoms.push({x, y, z});
  }

  // Count contacts (within 4.0 A)
  const contactCutoff = 4.0;
  const contacts = {};

  for (const [chain, atoms] of Object.entries(chainAtoms)) {
    contacts[chain] = 0;
    for (const ligAtom of ligandAtoms) {
      for (const protAtom of atoms) {
        const dx = ligAtom.x - protAtom.x;
        const dy = ligAtom.y - protAtom.y;
        const dz = ligAtom.z - protAtom.z;
        const dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
        if (dist < contactCutoff) {
          contacts[chain]++;
          break;  // Count each ligand atom once per chain
        }
      }
    }
  }

  // Get ligand center (should be near origin for proper symmetry)
  const ligandCenter = ligandAtoms.length > 0 ? {
    x: ligandAtoms.reduce((s, a) => s + a.x, 0) / ligandAtoms.length,
    y: ligandAtoms.reduce((s, a) => s + a.y, 0) / ligandAtoms.length,
    z: ligandAtoms.reduce((s, a) => s + a.z, 0) / ligandAtoms.length
  } : null;

  // Calculate chain centroids
  const chainCentroids = {};
  for (const [chain, atoms] of Object.entries(chainAtoms)) {
    if (atoms.length > 0) {
      chainCentroids[chain] = {
        x: atoms.reduce((s, a) => s + a.x, 0) / atoms.length,
        y: atoms.reduce((s, a) => s + a.y, 0) / atoms.length,
        z: atoms.reduce((s, a) => s + a.z, 0) / atoms.length
      };
    }
  }

  // Count clashes between chain A and chain B (atoms within 2.0 A)
  const clashCutoff = 2.0;
  let interchainClashes = 0;
  const chainA = chainAtoms['A'] || [];
  const chainB = chainAtoms['B'] || [];

  for (const atomA of chainA) {
    for (const atomB of chainB) {
      const dx = atomA.x - atomB.x;
      const dy = atomA.y - atomB.y;
      const dz = atomA.z - atomB.z;
      const dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
      if (dist < clashCutoff) {
        interchainClashes++;
      }
    }
  }

  // Calculate min distance between chains
  let minInterchainDist = Infinity;
  for (const atomA of chainA) {
    for (const atomB of chainB) {
      const dx = atomA.x - atomB.x;
      const dy = atomA.y - atomB.y;
      const dz = atomA.z - atomB.z;
      const dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
      if (dist < minInterchainDist) {
        minInterchainDist = dist;
      }
    }
  }

  return {
    chains: Object.keys(chainAtoms),
    contacts,
    ligandAtomCount: ligandAtoms.length,
    ligandCenter,
    chainAtomCounts: Object.fromEntries(
      Object.entries(chainAtoms).map(([k, v]) => [k, v.length])
    ),
    chainCentroids,
    interchainClashes,
    minInterchainDist: minInterchainDist === Infinity ? null : minInterchainDist
  };
}

async function main() {
  console.log('=== Testing Symmetric Ligand-at-Interface Design ===\n');
  console.log('This tests the NEW interface_ligand parameter that:');
  console.log('1. Orients ligand along C2 symmetry axis');
  console.log('2. Uses RASA conditioning (select_partially_buried)');
  console.log('3. Ensures both dimer chains contact the ligand\n');

  // Test 1: Generate C2 dimer with ligand at interface
  console.log('=== Test 1: C2 Dimer with Azobenzene at Interface ===\n');
  console.log('Submitting job with:');
  console.log('  - task: rfd3');
  console.log('  - ligand_smiles: c1ccc(cc1)N=Nc2ccccc2 (azobenzene)');
  console.log('  - conformer_method: xtb (quantum-accurate)');
  console.log('  - symmetry: { id: "C2" }');
  console.log('  - interface_ligand: true  <-- NEW PARAMETER');
  console.log('  - length: 50-70 (per chain)');
  console.log('');

  const startTime = Date.now();

  const result = await submitJob({
    task: 'rfd3',
    ligand_smiles: AZOBENZENE_SMILES,
    conformer_method: 'xtb',  // Use quantum conformer for accuracy
    symmetry: { id: 'C2' },
    interface_ligand: true,   // NEW: Place ligand at symmetric interface
    length: '50-70',
    num_timesteps: 200,
    step_scale: 1.5,
    num_designs: 1,
    seed: 42
  });

  const elapsed = ((Date.now() - startTime) / 1000).toFixed(1);
  console.log(`Job completed in ${elapsed}s\n`);

  // Check result
  if (result.status !== 'COMPLETED' && result.status !== 'completed') {
    console.log('ERROR: Job failed!');
    console.log('Status:', result.status);
    console.log('Error:', result.error || JSON.stringify(result).substring(0, 500));
    return;
  }

  console.log('Job completed successfully!\n');

  const design = result.output?.result?.designs?.[0] || result.result?.designs?.[0];
  if (!design) {
    console.log('ERROR: No design in result');
    console.log('Result structure:', JSON.stringify(result, null, 2).substring(0, 1000));
    return;
  }

  // Analyze the structure
  console.log('=== Analyzing Generated Structure ===\n');
  const analysis = analyzeInterfaceContacts(design.content);

  console.log(`Chains present: ${analysis.chains.join(', ')}`);
  console.log(`Chain atom counts: ${JSON.stringify(analysis.chainAtomCounts)}`);
  console.log(`Ligand atoms: ${analysis.ligandAtomCount}`);
  if (analysis.ligandCenter) {
    const c = analysis.ligandCenter;
    console.log(`Ligand center: (${c.x.toFixed(2)}, ${c.y.toFixed(2)}, ${c.z.toFixed(2)})`);
    const distFromOrigin = Math.sqrt(c.x*c.x + c.y*c.y + c.z*c.z);
    console.log(`Distance from origin: ${distFromOrigin.toFixed(2)} A`);
  }

  // Show chain positions
  console.log('\n=== Chain Geometry ===\n');
  if (analysis.chainCentroids['A']) {
    const c = analysis.chainCentroids['A'];
    console.log(`Chain A centroid: (${c.x.toFixed(2)}, ${c.y.toFixed(2)}, ${c.z.toFixed(2)})`);
  }
  if (analysis.chainCentroids['B']) {
    const c = analysis.chainCentroids['B'];
    console.log(`Chain B centroid: (${c.x.toFixed(2)}, ${c.y.toFixed(2)}, ${c.z.toFixed(2)})`);
  }
  console.log(`Interchain clashes (atoms < 2.0Å): ${analysis.interchainClashes}`);
  console.log(`Min interchain distance: ${analysis.minInterchainDist?.toFixed(2) || 'N/A'} A`);

  if (analysis.interchainClashes === 0) {
    console.log('NO CLASHES - Chains are properly separated');
  } else if (analysis.interchainClashes < 10) {
    console.log('MINOR CLASHES - Some overlap detected');
  } else {
    console.log('MAJOR CLASHES - Significant overlap between chains');
  }

  console.log('\n=== Interface Contact Analysis ===\n');
  for (const [chain, count] of Object.entries(analysis.contacts)) {
    if (chain !== 'L') {  // Skip ligand chain
      const status = count >= 3 ? 'GOOD' : 'POOR';
      console.log(`Chain ${chain} contacts with ligand: ${count} [${status}]`);
    }
  }

  // Validate interface placement
  const contactsA = analysis.contacts['A'] || 0;
  const contactsB = analysis.contacts['B'] || 0;
  const isSymmetric = contactsA >= 3 && contactsB >= 3;
  const symmetryScore = (contactsA + contactsB) > 0
    ? 1.0 - Math.abs(contactsA - contactsB) / (contactsA + contactsB)
    : 0;

  console.log('\n=== Validation Results ===\n');
  console.log(`Interface valid: ${isSymmetric ? 'YES' : 'NO'}`);
  console.log(`Symmetry score: ${(symmetryScore * 100).toFixed(1)}%`);

  if (isSymmetric) {
    console.log('\nSUCCESS: Ligand is at the symmetric interface!');
    console.log('Both chains have sufficient contacts with the ligand.');
  } else {
    console.log('\nWARNING: Ligand may not be properly at interface.');
    if (contactsA < 3) console.log('  - Chain A has too few contacts');
    if (contactsB < 3) console.log('  - Chain B has too few contacts');
  }

  // Run binding evaluation
  console.log('\n=== Running Binding Evaluation ===\n');
  const bindingResult = await submitJob({
    task: 'binding_eval',
    pdb_content: design.content,
    chain_a: 'A',
    chain_b: 'B',
    ligand_smiles: AZOBENZENE_SMILES,
    run_gnina: true,
    whole_protein_search: true
  });

  const binding = bindingResult.output || bindingResult;

  let bindingQuality = 'unknown';
  let affinityValue = null;

  if (binding.gnina_scoring?.status === 'completed') {
    const gnina = binding.gnina_scoring.result || {};
    affinityValue = gnina.best_affinity;
    console.log('GNINA Scoring:');
    console.log(`  Best affinity: ${gnina.best_affinity?.toFixed(2) || 'N/A'} kcal/mol`);
    console.log(`  Best CNN score: ${gnina.best_cnn_score?.toFixed(3) || 'N/A'}`);

    // Binding quality thresholds
    if (affinityValue !== null && affinityValue !== undefined) {
      if (affinityValue < -8) {
        bindingQuality = 'strong';
        console.log('  --> STRONG BINDER (drug-like affinity)');
      } else if (affinityValue < -5) {
        bindingQuality = 'moderate';
        console.log('  --> MODERATE BINDER (acceptable affinity)');
      } else if (affinityValue < -2) {
        bindingQuality = 'weak';
        console.log('  --> WEAK BINDER (low affinity)');
      } else if (affinityValue < 0) {
        bindingQuality = 'very_weak';
        console.log('  --> VERY WEAK BINDER (essentially no specific binding)');
      } else {
        bindingQuality = 'none';
        console.log('  --> NO BINDING (unfavorable, steric clashes)');
      }
    }
  }

  if (binding.interface_analysis?.status === 'completed') {
    const iface = binding.interface_analysis.metrics || {};
    console.log('\nInterface Analysis:');
    console.log(`  Interface residues: ${iface.nres_int || 0}`);
    console.log(`  Atomic contacts: ${iface.contacts || 0}`);
    console.log(`  Hydrogen bonds: ${iface.hbonds_int || 0}`);
  }

  // Save the design
  const outputPath = 'G:\\Github_local_repo\\Banta_Lab_RFdiffusion\\azobenzene_interface_design.pdb';
  fs.writeFileSync(outputPath, design.content);
  console.log(`\nDesign saved to: ${outputPath}`);

  // Final summary
  console.log('\n=== Test Summary ===\n');
  console.log('GEOMETRY:');
  console.log(`  - C2 symmetry: ${analysis.chains.includes('A') && analysis.chains.includes('B') ? 'YES' : 'NO'}`);
  console.log(`  - Ligand embedded: ${analysis.ligandAtomCount > 0 ? 'YES' : 'NO'}`);
  console.log(`  - Interface placement: ${isSymmetric ? 'VALID' : 'INVALID'}`);
  console.log(`  - Interchain clashes: ${analysis.interchainClashes}`);
  console.log(`  - Chain A contacts: ${contactsA}`);
  console.log(`  - Chain B contacts: ${contactsB}`);
  console.log(`  - Symmetry score: ${(symmetryScore * 100).toFixed(1)}%`);

  console.log('\nBINDING QUALITY:');
  console.log(`  - GNINA affinity: ${affinityValue !== null ? affinityValue.toFixed(2) + ' kcal/mol' : 'N/A'}`);
  console.log(`  - Binding quality: ${bindingQuality.toUpperCase()}`);

  // Overall assessment considering BOTH geometry AND binding
  const geometryValid = isSymmetric && analysis.ligandAtomCount > 0 && analysis.interchainClashes === 0;
  const bindingValid = bindingQuality === 'strong' || bindingQuality === 'moderate';
  const bindingAcceptable = bindingValid || bindingQuality === 'weak';

  console.log('\n=== OVERALL ASSESSMENT ===\n');
  if (geometryValid && bindingValid) {
    console.log('SUCCESS: Geometry valid AND binding affinity good!');
    console.log('This design is ready for experimental validation.');
  } else if (geometryValid && bindingAcceptable) {
    console.log('PARTIAL SUCCESS: Geometry valid but binding affinity is weak.');
    console.log('Consider generating more designs or optimizing the binding site.');
  } else if (geometryValid) {
    console.log('GEOMETRY ONLY: Structure is valid but binding is poor.');
    console.log('The binding site may not be properly designed for this ligand.');
  } else if (analysis.ligandAtomCount > 0) {
    console.log('NEEDS WORK: Ligand present but interface placement invalid.');
    console.log('Check for clashes or insufficient contacts.');
  } else {
    console.log('FAILED: Ligand not embedded in design.');
  }
}

main().catch(err => {
  console.error('Test failed with error:', err);
  process.exit(1);
});
