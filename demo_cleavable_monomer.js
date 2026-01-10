/**
 * DEMO: Cleavable Monomer for Interface Ligand Design
 *
 * This algorithm achieves stronger ligand binding at dimer interfaces by:
 * 1. Designing a monomer with ligand fully buried in center
 * 2. Finding loop regions where backbone crosses near the ligand
 * 3. Cleaving at optimal site to create dimer with ligand at interface
 *
 * Usage:
 *   node demo_cleavable_monomer.js [ligand] [num_designs] [create_homo]
 *
 *   ligand: "azobenzene" (default), "caffeine", or custom SMILES
 *   num_designs: Number of monomer designs to try (default: 3)
 *   create_homo: "true" to create homo-dimer via C2 symmetry (default: false)
 *
 * Expected Results:
 *   - Affinity: -5 to -8 kcal/mol (improvement over direct symmetric design)
 *   - Complete binding pocket (vs. split pocket from post-processing symmetry)
 *
 * See docs/CLEAVABLE_MONOMER_ALGORITHM.md for complete methodology.
 */

const http = require('http');
const fs = require('fs');
const path = require('path');

// Configuration
const CONFIG = {
  SERVER_HOST: 'localhost',
  SERVER_PORT: 8000,
  TARGET_LENGTH: '80-100',
  CLEAVAGE_STRATEGY: 'balanced'
};

// Common ligands
const LIGANDS = {
  'azobenzene': 'c1ccc(cc1)N=Nc2ccccc2',
  'caffeine': 'Cn1cnc2c1c(=O)n(c(=O)n2C)C',
  'benzene': 'c1ccccc1',
  'naphthalene': 'c1ccc2ccccc2c1'
};

// Parse command line arguments
const args = process.argv.slice(2);
const LIGAND_NAME = args[0] || 'azobenzene';
const LIGAND_SMILES = LIGANDS[LIGAND_NAME.toLowerCase()] || LIGAND_NAME;
const NUM_DESIGNS = parseInt(args[1]) || 3;
const CREATE_HOMO = args[2] === 'true';
const SEED = parseInt(args[3]) || Math.floor(Math.random() * 10000);

/**
 * Submit job to RFD3 server
 */
async function submitJob(input, timeout = 900000) {  // 15 min default
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
 * Get binding quality label
 */
function getBindingQuality(affinity) {
  if (affinity === null || affinity === undefined) return 'N/A';
  if (affinity < -8) return 'STRONG (drug-like)';
  if (affinity < -5) return 'MODERATE (useful)';
  if (affinity < -2) return 'WEAK (marginal)';
  if (affinity < 0) return 'VERY WEAK';
  return 'NO BINDING (clashing)';
}

/**
 * Format cleavage site info
 */
function formatCleavageSite(site) {
  if (!site) return 'N/A';
  return `Res ${site.residue_id} (${site.residue_name}), ` +
         `${site.distance_to_ligand.toFixed(1)}A from ligand, ` +
         `SS=${site.secondary_structure}, ` +
         `contacts: A=${site.contacts_chain_a}/B=${site.contacts_chain_b}`;
}

/**
 * Main demo function
 */
async function runDemo() {
  console.log('='.repeat(64));
  console.log('  DEMO: Cleavable Monomer for Interface Ligand Design');
  console.log('='.repeat(64));
  console.log();
  console.log(`Ligand: ${LIGAND_NAME} (${LIGAND_SMILES.substring(0, 40)}${LIGAND_SMILES.length > 40 ? '...' : ''})`);
  console.log(`Monomer designs: ${NUM_DESIGNS}`);
  console.log(`Chain length: ${CONFIG.TARGET_LENGTH} residues`);
  console.log(`Cleavage strategy: ${CONFIG.CLEAVAGE_STRATEGY}`);
  console.log(`Create homo-dimer: ${CREATE_HOMO}`);
  console.log(`Seed: ${SEED}`);
  console.log();

  const startTime = Date.now();

  // Run the cleavable monomer algorithm
  console.log('Running cleavable monomer algorithm...');
  console.log('  Step 1: Generate monomers with buried ligand');
  console.log('  Step 2: Find cleavage sites in loops near ligand');
  console.log('  Step 3: Score and select best cleavage site');
  console.log('  Step 4: Cleave to create dimer');
  console.log('  Step 5: Validate contacts, clashes, affinity');
  if (CREATE_HOMO) {
    console.log('  Step 6: Create homo-dimer with C2 symmetry');
  }
  console.log();

  const result = await submitJob({
    task: 'cleavable_monomer',
    ligand_smiles: LIGAND_SMILES,
    target_length: CONFIG.TARGET_LENGTH,
    num_designs: NUM_DESIGNS,
    cleavage_strategy: CONFIG.CLEAVAGE_STRATEGY,
    create_homo_dimer: CREATE_HOMO,
    seed: SEED,
    // Advanced options
    min_chain_length: 25,
    min_contacts_per_chain: 3,
    distance_range: [4.0, 12.0],
    remove_loop: true,
    loop_window: 3
  });

  const elapsed = ((Date.now() - startTime) / 1000).toFixed(1);
  console.log(`Completed in ${elapsed}s`);
  console.log();

  // Check result status
  const status = result.status || result.output?.status;
  if (status !== 'COMPLETED' && status !== 'completed') {
    console.log('='.repeat(64));
    console.log('  RESULT: FAILED');
    console.log('='.repeat(64));
    console.log();

    const error = result.error || result.output?.error || 'Unknown error';
    console.log(`Error: ${error}`);

    // Show diagnostics if available
    const diagnostics = result.diagnostics || result.output?.diagnostics;
    if (diagnostics) {
      console.log();
      console.log('Diagnostics:');
      console.log(`  Designs tried: ${diagnostics.num_designs_tried}`);
      if (diagnostics.suggestions) {
        console.log('  Suggestions:');
        for (const s of diagnostics.suggestions) {
          console.log(`    - ${s}`);
        }
      }
    }

    console.log();
    console.log('='.repeat(64));
    process.exit(1);
  }

  // Extract result data
  const data = result.output?.result || result.result;

  console.log('='.repeat(64));
  console.log('  RESULT: SUCCESS');
  console.log('='.repeat(64));
  console.log();

  // Dimer info
  const dimer = data.dimer;
  console.log('Dimer Structure:');
  console.log(`  Chain A: ${dimer.chain_a_length} residues`);
  console.log(`  Chain B: ${dimer.chain_b_length} residues`);
  console.log(`  Residues removed: ${dimer.residues_removed}`);
  console.log();

  // Cleavage site info
  const site = data.cleavage_site;
  console.log('Cleavage Site:');
  console.log(`  Residue: ${site.residue_id} (${site.residue_name})`);
  console.log(`  Distance to ligand: ${site.distance_to_ligand.toFixed(2)} A`);
  console.log(`  Secondary structure: ${site.secondary_structure} (C=coil/loop)`);
  console.log(`  Contacts before: ${site.contacts_chain_a}`);
  console.log(`  Contacts after: ${site.contacts_chain_b}`);
  console.log(`  Cleavage score: ${site.score.toFixed(2)}`);
  console.log();

  // Validation results
  const validation = data.validation;
  console.log('Validation:');
  console.log(`  Valid: ${validation.is_valid ? 'YES' : 'NO'}`);

  if (validation.chain_lengths) {
    console.log(`  Chains OK: A=${validation.chain_lengths.chain_a_ok}, B=${validation.chain_lengths.chain_b_ok}`);
  }

  if (validation.contacts) {
    console.log(`  Contacts: A=${validation.contacts.chain_a_contacts}, B=${validation.contacts.chain_b_contacts}`);
    console.log(`  Symmetry: ${(validation.contacts.symmetry_score * 100).toFixed(0)}%`);
  }

  if (validation.clashes) {
    console.log(`  Clashes: ${validation.clashes.has_clashes ? 'YES (bad)' : 'NO (good)'}`);
    if (validation.clashes.min_distance) {
      console.log(`  Min distance: ${validation.clashes.min_distance.toFixed(2)} A`);
    }
  }

  if (validation.gnina) {
    const gnina = validation.gnina;
    console.log(`  GNINA Affinity: ${gnina.affinity?.toFixed(2) || 'N/A'} kcal/mol`);
    console.log(`  Binding Quality: ${getBindingQuality(gnina.affinity)}`);
    if (gnina.cnn_score) {
      console.log(`  CNN Score: ${gnina.cnn_score.toFixed(3)}`);
    }
  }
  console.log();

  // Monomer info
  const monomer = data.monomer;
  console.log('Source Monomer:');
  console.log(`  Design index: ${monomer.design_index}`);
  console.log();

  // Homo-dimer info if created
  if (data.homo_dimer) {
    const homo = data.homo_dimer;
    console.log('Homo-dimer (C2 Symmetry):');
    console.log(`  Chain C: ${homo.chain_c_length} residues`);
    console.log(`  Chain D: ${homo.chain_d_length} residues`);
    if (homo.rmsd) {
      console.log(`  Symmetry RMSD: ${homo.rmsd.toFixed(2)} A`);
    }
    console.log();
  }

  // All results summary
  if (data.all_results && data.all_results.length > 0) {
    console.log('-'.repeat(64));
    console.log('  ALL DESIGNS SUMMARY');
    console.log('-'.repeat(64));
    console.log();
    console.log('Design | Score  | Valid | Affinity    | Contacts');
    console.log('-------|--------|-------|-------------|----------');

    for (const r of data.all_results) {
      const affStr = r.validation?.gnina?.affinity
        ? `${r.validation.gnina.affinity.toFixed(2).padStart(6)}`.padEnd(12)
        : 'N/A'.padEnd(12);
      const contacts = r.cleavage_site
        ? `A=${r.cleavage_site.contacts_chain_a}/B=${r.cleavage_site.contacts_chain_b}`
        : 'N/A';
      const validStr = r.validation?.is_valid ? 'YES' : 'NO ';
      console.log(
        `   ${r.design_index.toString().padStart(2)}  | ${r.overall_score.toFixed(1).padStart(6)} |  ${validStr}  | ${affStr}| ${contacts}`
      );
    }
    console.log();
  }

  // Save output files
  const baseFilename = `cleavable_${LIGAND_NAME.toLowerCase()}`;

  // Save dimer PDB
  const dimerPath = path.join(__dirname, `${baseFilename}_dimer.pdb`);
  fs.writeFileSync(dimerPath, dimer.pdb_content);
  console.log(`Saved dimer:   ${dimerPath}`);

  // Save monomer PDB
  const monomerPath = path.join(__dirname, `${baseFilename}_monomer.pdb`);
  fs.writeFileSync(monomerPath, monomer.pdb_content);
  console.log(`Saved monomer: ${monomerPath}`);

  // Save homo-dimer if created
  if (data.homo_dimer?.pdb_content) {
    const homoPath = path.join(__dirname, `${baseFilename}_homodimer.pdb`);
    fs.writeFileSync(homoPath, data.homo_dimer.pdb_content);
    console.log(`Saved homo-dimer: ${homoPath}`);
  }

  console.log();
  console.log('='.repeat(64));

  // Compare with direct symmetric design
  const affinity = validation.gnina?.affinity;
  console.log('  COMPARISON WITH DIRECT SYMMETRIC DESIGN');
  console.log('-'.repeat(64));
  console.log();
  console.log('  Method                  | Best Affinity | Pocket');
  console.log('  ------------------------|---------------|----------');
  console.log('  Direct symmetric (prev) | -3.78 kcal/mol| Split (2 halves)');
  console.log(`  Cleavable monomer       | ${affinity ? affinity.toFixed(2).padStart(5) : '  N/A'} kcal/mol| Complete`);
  console.log();

  if (affinity !== null && affinity !== undefined) {
    const improvement = -3.78 - affinity;
    if (improvement > 0) {
      console.log(`  Improvement: ${improvement.toFixed(2)} kcal/mol better`);
    } else if (improvement < 0) {
      console.log(`  Note: ${(-improvement).toFixed(2)} kcal/mol worse (try more designs)`);
    } else {
      console.log('  No change in affinity');
    }

    // Check if target reached
    if (affinity < -5) {
      console.log('  TARGET REACHED: < -5 kcal/mol (useful binding)');
    } else {
      console.log(`  Gap to target: ${(affinity - (-5)).toFixed(2)} kcal/mol`);
    }
  }

  console.log();
  console.log('='.repeat(64));
  console.log('  Tips:');
  console.log('  - Try more designs: node demo_cleavable_monomer.js azobenzene 10');
  console.log('  - Create homo-dimer: node demo_cleavable_monomer.js azobenzene 3 true');
  console.log('  - Longer chains: modify TARGET_LENGTH to "100-140"');
  console.log('='.repeat(64));
}

// Run the demo
runDemo().catch(err => {
  console.error('Demo failed:', err.message);
  console.error(err.stack);
  process.exit(1);
});
