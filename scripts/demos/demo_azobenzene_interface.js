/**
 * DEMO: Azobenzene Dimerization Binder
 *
 * Creates a C2 symmetric protein dimer with azobenzene ligand at the interface.
 *
 * Usage:
 *   node demo_azobenzene_interface.js [mode] [num_designs]
 *
 *   mode: "natural" (default) or "binder"
 *     - natural: Basic design with post-processing symmetry
 *     - binder: Uses H-bond and burial conditioning
 *
 *   num_designs: Number of designs to generate (default: 10)
 *
 * Expected Results:
 *   - Natural mode: Best affinity around -3.78 kcal/mol
 *   - Binder mode: Best affinity around -3.30 kcal/mol, better symmetry
 *
 * Target: -5 to -8 kcal/mol for useful binding (current gap: ~1-4 kcal/mol)
 *
 * See docs/INTERFACE_LIGAND_WORKFLOW.md for complete methodology.
 */

const http = require('http');
const fs = require('fs');
const path = require('path');

// Configuration
const CONFIG = {
  AZOBENZENE_SMILES: 'c1ccc(cc1)N=Nc2ccccc2',
  CHAIN_LENGTH: '80-100',
  SYMMETRY: { id: 'C2' },
  NUM_TIMESTEPS: 200,
  STEP_SCALE: 1.5,
  SERVER_HOST: 'localhost',
  SERVER_PORT: 8000
};

// Parse command line arguments
const args = process.argv.slice(2);
const MODE = args[0] || 'natural';
const NUM_DESIGNS = parseInt(args[1]) || 10;
const SEED = parseInt(args[2]) || Math.floor(Math.random() * 10000);

/**
 * Submit job to RFD3 server
 */
async function submitJob(input, timeout = 600000) {
  return new Promise((resolve, reject) => {
    const postData = JSON.stringify({ input });
    const options = {
      hostname: CONFIG.SERVER_HOST,
      port: CONFIG.SERVER_PORT,
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
        try {
          resolve(JSON.parse(data));
        } catch (e) {
          resolve({ error: data });
        }
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
 * Analyze PDB structure for interface validity
 */
function analyzeDesign(pdbContent) {
  const lines = pdbContent.split('\n');
  const atomLines = lines.filter(l => l.startsWith('ATOM'));
  const hetatmLines = lines.filter(l => l.startsWith('HETATM'));

  // Parse chain atoms
  const chainAtoms = {};
  for (const line of atomLines) {
    const chain = line[21];
    const x = parseFloat(line.substring(30, 38));
    const y = parseFloat(line.substring(38, 46));
    const z = parseFloat(line.substring(46, 54));
    if (!chainAtoms[chain]) chainAtoms[chain] = [];
    chainAtoms[chain].push({ x, y, z });
  }

  // Parse ligand atoms
  const ligandAtoms = [];
  for (const line of hetatmLines) {
    const x = parseFloat(line.substring(30, 38));
    const y = parseFloat(line.substring(38, 46));
    const z = parseFloat(line.substring(46, 54));
    ligandAtoms.push({ x, y, z });
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
        const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (dist < contactCutoff) {
          contacts[chain]++;
          break;
        }
      }
    }
  }

  // Check interchain clashes
  const clashCutoff = 2.0;
  let clashes = 0;
  const chainA = chainAtoms['A'] || [];
  const chainB = chainAtoms['B'] || [];
  for (const atomA of chainA) {
    for (const atomB of chainB) {
      const dx = atomA.x - atomB.x;
      const dy = atomA.y - atomB.y;
      const dz = atomA.z - atomB.z;
      const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
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

/**
 * Get binding quality label
 */
function getBindingQuality(affinity) {
  if (affinity === null) return 'N/A';
  if (affinity < -8) return 'STRONG (drug-like)';
  if (affinity < -5) return 'MODERATE (useful)';
  if (affinity < -2) return 'WEAK (marginal)';
  if (affinity < 0) return 'VERY WEAK';
  return 'NO BINDING (clashing)';
}

/**
 * Main demo function
 */
async function runDemo() {
  console.log('='.repeat(60));
  console.log('  DEMO: Azobenzene Dimerization Binder');
  console.log('='.repeat(60));
  console.log();
  console.log(`Mode: ${MODE.toUpperCase()}`);
  console.log(`Designs: ${NUM_DESIGNS}`);
  console.log(`Seed: ${SEED}`);
  console.log(`Chain length: ${CONFIG.CHAIN_LENGTH} residues`);
  console.log(`Ligand: Azobenzene (${CONFIG.AZOBENZENE_SMILES})`);
  console.log();

  const startTime = Date.now();

  // Step 1: Generate designs
  console.log('Step 1: Generating designs...');
  console.log();

  const interfaceLigand = MODE === 'binder' ? 'binder' : true;

  const rfdResult = await submitJob({
    task: 'rfd3',
    ligand_smiles: CONFIG.AZOBENZENE_SMILES,
    conformer_method: 'xtb',
    symmetry: CONFIG.SYMMETRY,
    interface_ligand: interfaceLigand,
    length: CONFIG.CHAIN_LENGTH,
    num_timesteps: CONFIG.NUM_TIMESTEPS,
    step_scale: CONFIG.STEP_SCALE,
    num_designs: NUM_DESIGNS,
    seed: SEED
  });

  if (rfdResult.status !== 'COMPLETED' && rfdResult.status !== 'completed') {
    console.log('ERROR: RFD3 job failed!');
    console.log('Status:', rfdResult.status);
    console.log('Error:', rfdResult.error || JSON.stringify(rfdResult).substring(0, 500));
    process.exit(1);
  }

  const designs = rfdResult.output?.result?.designs || rfdResult.result?.designs || [];
  const genTime = ((Date.now() - startTime) / 1000).toFixed(1);
  console.log(`Generated ${designs.length} designs in ${genTime}s`);
  console.log();

  if (designs.length === 0) {
    console.log('ERROR: No designs generated');
    process.exit(1);
  }

  // Step 2: Analyze and evaluate
  console.log('Step 2: Analyzing geometry and binding...');
  console.log();

  const results = [];

  for (let i = 0; i < designs.length; i++) {
    const design = designs[i];
    process.stdout.write(`  Design ${i + 1}/${designs.length}... `);

    const geom = analyzeDesign(design.content);
    let affinity = null;
    let cnnScore = null;
    let hbonds = 0;

    try {
      const bindingResult = await submitJob({
        task: 'binding_eval',
        pdb_content: design.content,
        chain_a: 'A',
        chain_b: 'B',
        ligand_smiles: CONFIG.AZOBENZENE_SMILES,
        run_gnina: true,
        whole_protein_search: true
      }, 120000);

      if (bindingResult.output?.gnina_scoring?.status === 'completed') {
        const gnina = bindingResult.output.gnina_scoring.result || {};
        affinity = gnina.best_affinity;
        cnnScore = gnina.best_cnn_score;
      }

      if (bindingResult.output?.interface_analysis?.metrics) {
        hbonds = bindingResult.output.interface_analysis.metrics.hbonds_int || 0;
      }
    } catch (err) {
      // Binding eval failed, continue with geometry only
    }

    const affStr = affinity !== null ? `${affinity.toFixed(2)} kcal/mol` : 'N/A';
    const validStr = geom.isValid ? 'VALID' : 'invalid';
    console.log(`${affStr}, ${geom.totalContacts} contacts, ${validStr}`);

    results.push({
      index: i,
      design: design,
      geometry: geom,
      affinity: affinity,
      cnnScore: cnnScore,
      hbonds: hbonds,
      score: (affinity !== null ? -affinity * 10 : 0) +
             (hbonds * 5) +
             geom.totalContacts +
             (geom.symmetryScore * 5) +
             (geom.isValid ? 10 : 0)
    });
  }

  // Sort by score
  results.sort((a, b) => b.score - a.score);

  // Step 3: Results summary
  console.log();
  console.log('='.repeat(60));
  console.log('  RESULTS SUMMARY');
  console.log('='.repeat(60));
  console.log();
  console.log('Rank | Design | Affinity   | Contacts | Symmetry | Valid');
  console.log('-----|--------|------------|----------|----------|------');

  for (let i = 0; i < Math.min(5, results.length); i++) {
    const r = results[i];
    const affStr = r.affinity !== null ? `${r.affinity.toFixed(2).padStart(6)}` : '   N/A';
    const symStr = `${(r.geometry.symmetryScore * 100).toFixed(0)}%`.padStart(4);
    console.log(
      `  ${(i + 1).toString().padStart(2)}  |   ${(r.index + 1).toString().padStart(2)}   | ${affStr} | ` +
      `   ${r.geometry.totalContacts.toString().padStart(2)}     |   ${symStr}   |  ${r.geometry.isValid ? 'YES' : 'NO '}`
    );
  }

  // Best design details
  const best = results[0];
  console.log();
  console.log('-'.repeat(60));
  console.log('  BEST DESIGN');
  console.log('-'.repeat(60));
  console.log();
  console.log(`  Design #${best.index + 1}`);
  console.log(`  GNINA Affinity: ${best.affinity?.toFixed(2) || 'N/A'} kcal/mol`);
  console.log(`  CNN Score: ${best.cnnScore?.toFixed(3) || 'N/A'}`);
  console.log(`  H-bonds: ${best.hbonds}`);
  console.log(`  Contacts: Chain A = ${best.geometry.contactsA}, Chain B = ${best.geometry.contactsB}`);
  console.log(`  Symmetry: ${(best.geometry.symmetryScore * 100).toFixed(1)}%`);
  console.log(`  Interface Valid: ${best.geometry.isValid ? 'YES' : 'NO'}`);
  console.log(`  Binding Quality: ${getBindingQuality(best.affinity)}`);

  // Save best design
  const outputPath = path.join(__dirname, `azobenzene_demo_${MODE}_best.pdb`);
  fs.writeFileSync(outputPath, best.design.content);
  console.log();
  console.log(`  Saved to: ${outputPath}`);

  // Statistics
  const validCount = results.filter(r => r.geometry.isValid).length;
  const negativeCount = results.filter(r => r.affinity !== null && r.affinity < 0).length;

  console.log();
  console.log('-'.repeat(60));
  console.log('  STATISTICS');
  console.log('-'.repeat(60));
  console.log();
  console.log(`  Valid interfaces: ${validCount}/${results.length} (${(validCount/results.length*100).toFixed(0)}%)`);
  console.log(`  Negative affinity: ${negativeCount}/${results.length} (${(negativeCount/results.length*100).toFixed(0)}%)`);

  const totalTime = ((Date.now() - startTime) / 1000).toFixed(1);
  console.log();
  console.log(`  Total time: ${totalTime}s`);
  console.log();
  console.log('='.repeat(60));

  // Mode comparison reference
  if (MODE === 'natural') {
    console.log('  TIP: Try binder mode for H-bond conditioning:');
    console.log('       node demo_azobenzene_interface.js binder');
  } else {
    console.log('  TIP: Compare with natural mode:');
    console.log('       node demo_azobenzene_interface.js natural');
  }
  console.log('='.repeat(60));
}

// Run the demo
runDemo().catch(err => {
  console.error('Demo failed:', err.message);
  process.exit(1);
});
