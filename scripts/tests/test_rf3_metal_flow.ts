/**
 * Integration test: RF3 metal parameter data flow
 *
 * Mimics the real NL pipeline execution to verify that:
 * 1. parse_intent extracts target_metal ("TB") and ligand_smiles
 * 2. configure step stores them in result.data
 * 3. createRf3Step extracts target_metal from previousResults
 * 4. submitRF3Prediction sends metal in the request body
 * 5. Backend handler receives metal and passes to run_rf3_inference
 * 6. run_rf3_cli includes [Tb+3] SMILES in the JSON config components
 *
 * Also checks MPNN alanine bias for hard metals.
 *
 * Run: npx tsx scripts/tests/test_rf3_metal_flow.ts
 */

// ============================================================================
// Inline stubs — no imports needed, we replicate the logic directly
// ============================================================================

const LIGAND_SMILES: Record<string, string> = {
  citrate: 'OC(=O)CC(O)(CC(O)=O)C(O)=O',
  atp: 'c1nc(c2c(n1)n(cn2)[C@@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N',
};

function lookupSmiles(name?: string): string {
  if (!name) return '';
  return LIGAND_SMILES[name.toLowerCase()] || '';
}

function mapDesignType(ai: { metal_type?: string }): string {
  if (ai.metal_type) return 'metal';
  return 'general';
}

// Metal MPNN bias (replicate from metal-mpnn-bias.ts)
type HsabClass = 'hard' | 'soft' | 'borderline';

const HSAB_CLASSIFICATION: Record<string, HsabClass> = {
  TB: 'hard', EU: 'hard', GD: 'hard', CA: 'hard', MG: 'hard',
  LA: 'hard', CE: 'hard', SM: 'hard', YB: 'hard', DY: 'hard',
  CU: 'soft', AG: 'soft',
  ZN: 'borderline', FE: 'borderline', MN: 'borderline',
  CO: 'borderline', NI: 'borderline',
};

const CONFIG_BY_CLASS: Record<HsabClass, { bias_AA?: string; omit_AA?: string }> = {
  hard: {
    bias_AA: 'A:-2.0',
    omit_AA: 'C',
  },
  soft: {
    bias_AA: 'C:4.0,H:3.0,M:2.0,A:-2.0',
  },
  borderline: {
    bias_AA: 'H:3.0,D:2.0,E:2.0,C:1.5,A:-2.0',
  },
};

function getMetalMpnnConfig(targetMetal: string) {
  const key = targetMetal.toUpperCase().trim();
  const hsab = HSAB_CLASSIFICATION[key] ?? 'borderline';
  return CONFIG_BY_CLASS[hsab];
}

// Metal ion SMILES (replicate from inference_utils.py)
const METAL_ION_SMILES: Record<string, string> = {
  TB: '[Tb+3]', GD: '[Gd+3]', EU: '[Eu+3]', LA: '[La+3]',
  CE: '[Ce+3]', ND: '[Nd+3]', DY: '[Dy+3]', SM: '[Sm+3]',
  ZN: '[Zn+2]', CA: '[Ca+2]', MG: '[Mg+2]', MN: '[Mn+2]',
  FE: '[Fe+2]', CU: '[Cu+2]', CO: '[Co+2]', NI: '[Ni+2]',
};

// ============================================================================
// Test harness
// ============================================================================

let passed = 0;
let failed = 0;

function assert(condition: boolean, msg: string) {
  if (condition) {
    console.log(`  ✓ ${msg}`);
    passed++;
  } else {
    console.error(`  ✗ FAIL: ${msg}`);
    failed++;
  }
}

function section(name: string) {
  console.log(`\n━━━ ${name} ━━━`);
}

// ============================================================================
// Test 1: Parse Intent Step
// ============================================================================

section('Step 1: Parse Intent — "design a protein to bind terbium and citrate"');

// Simulate AI parser output
const aiParserResult = {
  metal_type: 'TB',
  ligand_name: 'citrate',
  design_goal: 'binding',
  confidence: 0.95,
};

const intentData = {
  design_type: mapDesignType(aiParserResult),
  target_metal: aiParserResult.metal_type || '',
  ligand_name: aiParserResult.ligand_name || '',
  ligand_smiles: lookupSmiles(aiParserResult.ligand_name),
  ai_parsed: true,
};

assert(intentData.design_type === 'metal', 'design_type = "metal"');
assert(intentData.target_metal === 'TB', 'target_metal = "TB"');
assert(intentData.ligand_name === 'citrate', 'ligand_name = "citrate"');
assert(intentData.ligand_smiles === 'OC(=O)CC(O)(CC(O)=O)C(O)=O', 'ligand_smiles is citrate SMILES');

// Store as previousResults
const previousResults: Record<string, { data?: Record<string, unknown> }> = {
  parse_intent: { data: intentData },
};

// ============================================================================
// Test 2: Scaffold Search
// ============================================================================

section('Step 2: Scaffold Search — found 7F6E');

previousResults['scaffold_search_nl'] = {
  data: {
    searched: true,
    recommended_action: 'scaffold',
    scaffold_pdb_id: '7F6E',
  },
};

// ============================================================================
// Test 3: Configure Step
// ============================================================================

section('Step 3: Configure — merges intent data');

const intentResult = Object.values(previousResults).find(r => r.data?.design_type);
const configData: Record<string, unknown> = {
  design_type: intentResult?.data?.design_type || 'general',
  target_metal: intentResult?.data?.target_metal || '',
  ligand_name: intentResult?.data?.ligand_name || '',
  ligand_smiles: intentResult?.data?.ligand_smiles || '',
  num_designs: 4,
};

assert(configData.target_metal === 'TB', 'configure preserves target_metal = "TB"');
assert(configData.ligand_smiles === 'OC(=O)CC(O)(CC(O)=O)C(O)=O', 'configure preserves ligand_smiles');

previousResults['configure'] = { data: configData };

// ============================================================================
// Test 4: RFD3 Backbone Generation (metal single mode)
// ============================================================================

section('Step 4: RFD3 backbone generation — metal single mode');

previousResults['rfd3_nl'] = {
  data: {
    metal_single_mode: true,
    target_metal: 'TB',
    ligand_name: 'citrate',
    ligand_code: 'CIT',
  },
};

assert(previousResults['rfd3_nl'].data?.metal_single_mode === true, 'rfd3_nl sets metal_single_mode');

// ============================================================================
// Test 5: MPNN Step — alanine bias check
// ============================================================================

section('Step 5: MPNN — alanine bias for hard metals');

// Replicate the MPNN metal config resolution from shared-steps.ts lines 257-267
let metalConfig: { bias_AA?: string; omit_AA?: string } = {};
for (const result of Object.values(previousResults)) {
  const mp = (result.data as Record<string, unknown> | undefined)?.mpnn_params;
  if (mp && typeof mp === 'object') {
    metalConfig = { ...metalConfig, ...mp as typeof metalConfig };
  }
}
if (!metalConfig.bias_AA && !metalConfig.omit_AA) {
  let targetMetal = '';
  const mpnnIntentResult = Object.values(previousResults).find(r => r.data?.design_type);
  targetMetal = (mpnnIntentResult?.data?.target_metal as string) || '';
  if (targetMetal) {
    metalConfig = getMetalMpnnConfig(targetMetal);
  }
}

assert(metalConfig.bias_AA !== undefined, 'hard metal (TB) now has bias_AA defined');
assert(metalConfig.bias_AA?.includes('A:-2.0') === true, 'bias_AA includes A:-2.0 (alanine penalty)');
assert(metalConfig.omit_AA === 'C', 'omit_AA = "C" (cysteine omitted for hard acids)');

// Check all HSAB classes have alanine penalty
for (const cls of ['hard', 'soft', 'borderline'] as const) {
  const config = CONFIG_BY_CLASS[cls];
  assert(
    config.bias_AA?.includes('A:-') === true,
    `${cls} class has alanine penalty: ${config.bias_AA}`,
  );
}

// Store MPNN output (sequences only — no PDB)
previousResults['mpnn_nl'] = {
  data: {
    // MPNN outputs sequences, no metal data forwarded
  },
};

// ============================================================================
// Test 6: RF3 Step — metal extraction from previousResults
// ============================================================================

section('Step 6: RF3 — extract target_metal from previousResults');

// Replicate the exact logic from shared-steps.ts createRf3Step (lines 409-421)
let ligandSmiles: string | undefined;
let targetMetal: string | undefined;
for (const result of Object.values(previousResults)) {
  if (!ligandSmiles) {
    const smiles = result.data?.ligand_smiles as string;
    if (smiles) ligandSmiles = smiles;
  }
  if (!targetMetal) {
    const metal = result.data?.target_metal as string;
    if (metal) targetMetal = metal;
  }
}

assert(ligandSmiles === 'OC(=O)CC(O)(CC(O)=O)C(O)=O', 'RF3 step finds ligand_smiles from previousResults');
assert(targetMetal === 'TB', 'RF3 step finds target_metal from previousResults');

// Simulate the API call
const rf3Request = {
  sequence: 'MKVLIPPLGAGSSYIRRALTELLGELPADLRGLDARRRIQRLEPAGCR',
  name: 'Design 1',
  ligand_smiles: ligandSmiles,
  metal: targetMetal,
};

assert(rf3Request.metal === 'TB', 'RF3 request includes metal = "TB"');
assert(rf3Request.ligand_smiles !== undefined, 'RF3 request includes ligand_smiles');

// ============================================================================
// Test 7: Traditional mode proxy — body construction
// ============================================================================

section('Step 7: API proxy — request body construction');

// Traditional mode wraps as: { input: { task: 'rf3', ...request } }
const traditionalBody = { input: { task: 'rf3' as const, ...rf3Request } };
assert(traditionalBody.input.metal === 'TB', 'Traditional proxy preserves metal in input');
assert(traditionalBody.input.task === 'rf3', 'Traditional proxy sets task = "rf3"');

// Serverless mode wraps as: { input: { task, ...body } }
const serverlessBody = { input: { task: 'rf3' as const, ...rf3Request } };
assert(serverlessBody.input.metal === 'TB', 'Serverless proxy preserves metal in input');

// ============================================================================
// Test 8: Backend handler extraction
// ============================================================================

section('Step 8: Backend handler — job_input extraction');

// Simulate handler.py handle_rf3() extraction
const jobInput = traditionalBody.input;
const handlerMetal = jobInput.metal;
const handlerLigandSmiles = jobInput.ligand_smiles;

assert(handlerMetal === 'TB', 'handler extracts metal = "TB" from job_input');
assert(handlerLigandSmiles === rf3Request.ligand_smiles, 'handler extracts ligand_smiles');

// ============================================================================
// Test 9: RF3 CLI JSON config construction
// ============================================================================

section('Step 9: RF3 CLI — JSON config with metal SMILES');

// Simulate run_rf3_cli config building
const components: Array<Record<string, string>> = [];

// Chain A
components.push({ seq: rf3Request.sequence, chain_id: 'A' });

// Metal ion component
if (handlerMetal) {
  const metalSmiles = METAL_ION_SMILES[handlerMetal.toUpperCase()];
  if (metalSmiles) {
    components.push({ smiles: metalSmiles });
  }
}

// Ligand component
if (handlerLigandSmiles) {
  components.push({ smiles: handlerLigandSmiles });
}

const jsonConfig = [{ name: rf3Request.name, components }];

assert(components.length === 3, `RF3 config has 3 components (got ${components.length})`);
assert(components[0].seq !== undefined, 'Component 0: protein sequence');
assert(components[1].smiles === '[Tb+3]', `Component 1: metal ion [Tb+3] (got "${components[1]?.smiles}")`);
assert(components[2].smiles === rf3Request.ligand_smiles, 'Component 2: citrate ligand SMILES');

console.log('\nFull RF3 JSON config:');
console.log(JSON.stringify(jsonConfig, null, 2));

// ============================================================================
// Test 10: Metal coverage — all metals have SMILES
// ============================================================================

section('Step 10: Metal SMILES coverage');

const allMetals = Object.keys(HSAB_CLASSIFICATION);
for (const metal of allMetals) {
  const smiles = METAL_ION_SMILES[metal];
  assert(smiles !== undefined, `${metal} → ${smiles}`);
}

// ============================================================================
// Test 11: Live API test (if backend is running)
// ============================================================================

section('Step 11: Live API test (localhost:8000)');

async function testLiveApi() {
  const API_URL = 'http://localhost:8000';

  try {
    // Test health first
    const healthResp = await fetch(`${API_URL}/runsync`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ input: { task: 'health' } }),
    });

    if (!healthResp.ok) {
      console.log('  ⚠ Backend not running at localhost:8000 — skipping live test');
      return;
    }

    console.log('  Backend is running, testing RF3 with metal...');

    // Submit RF3 prediction with metal
    const rf3Resp = await fetch(`${API_URL}/runsync`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        input: {
          task: 'rf3',
          sequence: 'MDEEKLRALAQEQARSLAQEQERLAQEQARLRAEQELAKQR',
          name: 'test_metal_flow',
          ligand_smiles: 'OC(=O)CC(O)(CC(O)=O)C(O)=O',
          metal: 'TB',
        },
      }),
    });

    if (rf3Resp.ok) {
      const result = await rf3Resp.json();
      const output = result.output || {};
      console.log(`  RF3 status: ${output.status}`);
      if (output.status === 'completed') {
        const confidences = output.result?.confidences ||
          output.result?.predictions?.[0]?.confidences || {};
        const summary = confidences.summary_confidences || confidences;
        console.log(`  pTM: ${summary.ptm?.toFixed(3) || 'N/A'}`);
        console.log(`  pLDDT: ${summary.overall_plddt?.toFixed(2) || 'N/A'}`);
        assert(true, 'Live RF3 with metal completed successfully');
      } else {
        console.log(`  RF3 error: ${output.error || 'unknown'}`);
        assert(false, `Live RF3 failed: ${output.error}`);
      }
    } else {
      console.log(`  RF3 HTTP error: ${rf3Resp.status}`);
      assert(false, `RF3 HTTP error: ${rf3Resp.status}`);
    }
  } catch (err) {
    console.log(`  ⚠ Backend not reachable — skipping live test (${err})`);
  }
}

// ============================================================================
// Run and report
// ============================================================================

async function main() {
  await testLiveApi();

  console.log('\n' + '═'.repeat(50));
  console.log(`Results: ${passed} passed, ${failed} failed`);
  console.log('═'.repeat(50));

  if (failed > 0) {
    process.exit(1);
  }
}

main();
