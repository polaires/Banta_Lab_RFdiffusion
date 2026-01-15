/**
 * Test: Generate Multiple Designs and Filter by Binding Affinity
 *
 * This script generates 5 designs with interface_ligand mode,
 * evaluates each for binding affinity using GNINA,
 * and reports the best one.
 */

const http = require('http');
const fs = require('fs');

const AZOBENZENE_SMILES = 'c1ccc(cc1)N=Nc2ccccc2';
const NUM_DESIGNS = 10;

async function submitJob(input, timeout = 600000) {
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
      timeout: timeout
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

function analyzeDesign(pdbContent) {
  const lines = pdbContent.split('\n');
  const atomLines = lines.filter(l => l.startsWith('ATOM'));
  const hetatmLines = lines.filter(l => l.startsWith('HETATM'));

  const chainAtoms = {};
  for (const line of atomLines) {
    const chain = line[21];
    const x = parseFloat(line.substring(30, 38));
    const y = parseFloat(line.substring(38, 46));
    const z = parseFloat(line.substring(46, 54));
    if (!chainAtoms[chain]) chainAtoms[chain] = [];
    chainAtoms[chain].push({x, y, z});
  }

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
          break;
        }
      }
    }
  }

  // Check clashes between chains
  const clashCutoff = 2.0;
  let clashes = 0;
  const chainA = chainAtoms['A'] || [];
  const chainB = chainAtoms['B'] || [];
  for (const atomA of chainA) {
    for (const atomB of chainB) {
      const dx = atomA.x - atomB.x;
      const dy = atomA.y - atomB.y;
      const dz = atomA.z - atomB.z;
      const dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
      if (dist < clashCutoff) clashes++;
    }
  }

  const contactsA = contacts['A'] || 0;
  const contactsB = contacts['B'] || 0;
  const symmetryScore = (contactsA + contactsB) > 0
    ? 1.0 - Math.abs(contactsA - contactsB) / (contactsA + contactsB)
    : 0;

  return {
    chains: Object.keys(chainAtoms),
    contactsA,
    contactsB,
    totalContacts: contactsA + contactsB,
    symmetryScore,
    clashes,
    ligandAtoms: ligandAtoms.length,
    isValid: contactsA >= 3 && contactsB >= 3 && clashes === 0
  };
}

async function main() {
  console.log('=== Multi-Design Generation and Filtering ===\n');
  console.log(`Generating ${NUM_DESIGNS} designs with interface_ligand mode...`);
  console.log('Will evaluate each with GNINA and select the best.\n');

  const startTime = Date.now();

  // Generate designs
  console.log('Step 1: Generating designs...\n');
  const rfdResult = await submitJob({
    task: 'rfd3',
    ligand_smiles: AZOBENZENE_SMILES,
    conformer_method: 'xtb',
    symmetry: { id: 'C2' },
    interface_ligand: true,
    length: '80-100',
    num_timesteps: 200,
    step_scale: 1.5,
    num_designs: NUM_DESIGNS,
    seed: 456
  });

  if (rfdResult.status !== 'COMPLETED' && rfdResult.status !== 'completed') {
    console.log('ERROR: RFD3 job failed!');
    console.log('Status:', rfdResult.status);
    console.log('Error:', rfdResult.error || JSON.stringify(rfdResult).substring(0, 500));
    return;
  }

  const designs = rfdResult.output?.result?.designs || rfdResult.result?.designs || [];
  console.log(`Generated ${designs.length} designs in ${((Date.now() - startTime) / 1000).toFixed(1)}s\n`);

  if (designs.length === 0) {
    console.log('ERROR: No designs generated');
    return;
  }

  // Analyze and evaluate each design
  console.log('Step 2: Analyzing geometry and running binding evaluation...\n');
  const results = [];

  for (let i = 0; i < designs.length; i++) {
    const design = designs[i];
    console.log(`--- Design ${i + 1}/${designs.length} ---`);

    // Geometric analysis
    const geom = analyzeDesign(design.content);
    console.log(`  Contacts: A=${geom.contactsA}, B=${geom.contactsB}, Total=${geom.totalContacts}`);
    console.log(`  Symmetry: ${(geom.symmetryScore * 100).toFixed(0)}%`);
    console.log(`  Clashes: ${geom.clashes}`);
    console.log(`  Valid: ${geom.isValid ? 'YES' : 'NO'}`);

    // Binding evaluation
    let bindingResult = null;
    let affinity = null;
    let cnnScore = null;

    try {
      bindingResult = await submitJob({
        task: 'binding_eval',
        pdb_content: design.content,
        chain_a: 'A',
        chain_b: 'B',
        ligand_smiles: AZOBENZENE_SMILES,
        run_gnina: true,
        whole_protein_search: true
      }, 120000);

      if (bindingResult.output?.gnina_scoring?.status === 'completed') {
        const gnina = bindingResult.output.gnina_scoring.result || {};
        affinity = gnina.best_affinity;
        cnnScore = gnina.best_cnn_score;
        console.log(`  GNINA affinity: ${affinity?.toFixed(2) || 'N/A'} kcal/mol`);
        console.log(`  CNN score: ${cnnScore?.toFixed(3) || 'N/A'}`);
      } else {
        console.log('  GNINA: Not available');
      }
    } catch (err) {
      console.log(`  Binding eval error: ${err.message}`);
    }

    results.push({
      index: i,
      design: design,
      geometry: geom,
      affinity: affinity,
      cnnScore: cnnScore,
      // Composite score: prioritize affinity, then contacts, then symmetry
      score: (affinity !== null ? -affinity * 10 : 0) +
             geom.totalContacts +
             (geom.symmetryScore * 5) +
             (geom.isValid ? 10 : 0)
    });

    console.log('');
  }

  // Sort by score (higher is better)
  results.sort((a, b) => b.score - a.score);

  // Summary
  console.log('=== Results Summary ===\n');
  console.log('Rank | Design | Affinity | Contacts | Symmetry | Valid | Score');
  console.log('-----|--------|----------|----------|----------|-------|------');

  for (let i = 0; i < results.length; i++) {
    const r = results[i];
    const affStr = r.affinity !== null ? `${r.affinity.toFixed(2)}` : 'N/A';
    console.log(
      `  ${i + 1}  |   ${r.index + 1}    |  ${affStr.padStart(6)}  |    ${r.geometry.totalContacts.toString().padStart(2)}    |   ${(r.geometry.symmetryScore * 100).toFixed(0).padStart(3)}%   |  ${r.geometry.isValid ? 'YES' : 'NO '}  | ${r.score.toFixed(1)}`
    );
  }

  // Best design
  const best = results[0];
  console.log('\n=== Best Design ===\n');
  console.log(`Design #${best.index + 1}`);
  console.log(`  GNINA affinity: ${best.affinity?.toFixed(2) || 'N/A'} kcal/mol`);
  console.log(`  CNN score: ${best.cnnScore?.toFixed(3) || 'N/A'}`);
  console.log(`  Contacts: A=${best.geometry.contactsA}, B=${best.geometry.contactsB}`);
  console.log(`  Symmetry score: ${(best.geometry.symmetryScore * 100).toFixed(1)}%`);
  console.log(`  Interface valid: ${best.geometry.isValid ? 'YES' : 'NO'}`);

  // Binding quality assessment
  if (best.affinity !== null) {
    let quality = 'none';
    if (best.affinity < -8) quality = 'STRONG (drug-like)';
    else if (best.affinity < -5) quality = 'MODERATE (useful)';
    else if (best.affinity < -2) quality = 'WEAK (marginal)';
    else if (best.affinity < 0) quality = 'VERY WEAK (thermal noise)';
    else quality = 'NO BINDING';
    console.log(`  Binding quality: ${quality}`);
  }

  // Save best design
  const outputPath = 'G:\\Github_local_repo\\Banta_Lab_RFdiffusion\\best_interface_design.pdb';
  fs.writeFileSync(outputPath, best.design.content);
  console.log(`\nBest design saved to: ${outputPath}`);

  // Also save all designs for comparison
  for (let i = 0; i < results.length; i++) {
    const path = `G:\\Github_local_repo\\Banta_Lab_RFdiffusion\\interface_design_${i + 1}.pdb`;
    fs.writeFileSync(path, results[i].design.content);
  }
  console.log(`All ${results.length} designs saved for comparison.`);

  const totalTime = ((Date.now() - startTime) / 1000).toFixed(1);
  console.log(`\nTotal time: ${totalTime}s`);
}

main().catch(err => {
  console.error('Test failed with error:', err);
  process.exit(1);
});
