/**
 * Test script for RFD3 Fe→Tb conversion
 * Uses the same infrastructure as the frontend
 */

const fs = require('fs');
const https = require('https');
const http = require('http');

// Fetch PDB from RCSB
async function fetchPdb(pdbId) {
  return new Promise((resolve, reject) => {
    const url = `https://files.rcsb.org/download/${pdbId.toUpperCase()}.pdb`;
    https.get(url, (res) => {
      let data = '';
      res.on('data', chunk => data += chunk);
      res.on('end', () => resolve(data));
      res.on('error', reject);
    }).on('error', reject);
  });
}

// Submit job to backend
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
      }
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
    req.write(postData);
    req.end();
  });
}

async function main() {
  console.log('=== RFD3 Fe→Tb Conversion Test ===\n');

  // Step 1: Fetch PDB
  console.log('1. Fetching 1BRF from RCSB...');
  let pdbContent = await fetchPdb('1BRF');
  console.log(`   Downloaded ${pdbContent.length} bytes`);

  // Step 1b: Replace FE with TB in the PDB
  console.log('   Replacing FE with TB for metal conversion...');
  pdbContent = pdbContent.replace(/HETATM(.{12})FE(.{4}) FE /g, 'HETATM$1TB$2 TB ');
  pdbContent = pdbContent.replace(/ FE\s*$/gm, ' TB');
  console.log('   Metal replaced: FE → TB\n');

  // Step 2: Analyze structure first
  console.log('2. Analyzing structure...');
  const analysis = await submitJob({
    task: 'analyze',
    pdb_content: pdbContent
  });

  console.log('   Analysis result:');
  if (analysis.status === 'completed') {
    const result = analysis.output?.result || analysis.result;
    console.log(`   - Residues: ${result?.num_residues}`);
    console.log(`   - Ligands found: ${JSON.stringify(result?.ligands?.map(l => l.name))}`);
    console.log(`   - Binding sites: ${result?.binding_sites?.length}`);
    if (result?.suggestions?.length > 0) {
      console.log(`   - Suggestions: ${result.suggestions[0].description}`);
    }
  } else {
    console.log('   Error:', analysis.error);
  }

  // Step 3: Try RFD3 with proper format
  // The contig format for partial diffusion needs to include the metal
  // For metal replacement, we need to modify the PDB first or use specific format
  console.log('\n3. Submitting RFD3 design job...');

  // Try without ligand parameter first - just partial diffusion of the pocket
  const result = await submitJob({
    task: 'rfd3',
    contig: 'A1-53',
    partial_t: 12.0,
    num_timesteps: 100,
    step_scale: 1.3,
    gamma_0: 0.8,
    num_designs: 1,
    seed: 42,
    pdb_content: pdbContent
  });

  console.log('4. RFD3 Result:');
  console.log(JSON.stringify(result, null, 2));

  if (result.status === 'COMPLETED' || result.status === 'completed') {
    console.log('\n=== RFD3 SUCCESS ===');

    // Extract designed PDB
    const designedPdb = result.output?.result?.designs?.[0]?.content;
    if (designedPdb) {
      console.log(`   Generated structure with ${designedPdb.split('ATOM').length - 1} atoms`);

      // Step 5: Run MPNN to design sequences
      console.log('\n5. Running MPNN sequence design...');
      const mpnnResult = await submitJob({
        task: 'mpnn',
        pdb_content: designedPdb,
        num_sequences: 3,
        temperature: 0.1,
        model_type: 'ligand_mpnn'
      });

      console.log('6. MPNN Result:');
      if (mpnnResult.status === 'COMPLETED' || mpnnResult.status === 'completed') {
        console.log('   Status: SUCCESS');
        const sequences = mpnnResult.output?.result?.sequences || [];
        console.log(`   Generated ${sequences.length} sequences`);
        if (sequences.length > 0) {
          console.log('   Top sequence:');
          console.log(`   ${sequences[0]?.sequence?.substring(0, 60)}...`);
        }
      } else {
        console.log('   Status: ' + mpnnResult.status);
        console.log('   Error:', mpnnResult.error);
      }
    }

    console.log('\n=== PIPELINE COMPLETE ===');
    console.log('Next steps for Fe→Tb conversion:');
    console.log('1. Modify the designed structure to place TB at the metal site');
    console.log('2. Run RF3 to validate the designed sequences fold correctly');
    console.log('3. Evaluate coordination geometry for TB');

  } else if (result.error) {
    console.log('\n=== ERROR ===');
    console.log('Error:', result.error);
  }
}

main().catch(console.error);
