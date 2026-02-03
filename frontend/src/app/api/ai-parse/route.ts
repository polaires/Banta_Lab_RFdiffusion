/**
 * Vercel API route for AI intent parsing.
 *
 * Calls an OpenAI-compatible chat completions endpoint so that natural-language
 * design queries are parsed without touching RunPod / GPU infrastructure.
 *
 * Uses gpt-5-mini by default — fast and cheap enough for structured extraction.
 *
 * Env vars (set in Vercel dashboard):
 *   LLM_API_KEY   – required
 *   LLM_API_URL   – optional, defaults to https://yinli.one
 *   LLM_MODEL     – optional, defaults to gpt-5-mini
 */

import { NextRequest, NextResponse } from 'next/server';

// ---------------------------------------------------------------------------
// Config
// ---------------------------------------------------------------------------

const LLM_API_KEY = process.env.LLM_API_KEY;
const LLM_API_URL = process.env.LLM_API_URL || 'https://yinli.one';
const LLM_MODEL = process.env.LLM_MODEL || 'gpt-5-mini';

// ---------------------------------------------------------------------------
// System prompt & user template — ported from nl_design_parser.py
// ---------------------------------------------------------------------------

const PARSER_SYSTEM_PROMPT = `You are a protein design assistant that extracts structured information from natural language design requests.

Your task is to parse user queries about protein design and extract key parameters.

For each query, extract:
1. metal_type: The metal ion involved (use standard symbols: ZN, FE, CA, TB, EU, etc.)
2. ligand_name: The ligand/cofactor name (citrate, PQQ, heme, etc.)
3. design_goal: One of: binding, catalysis, sensing, structural
4. target_topology: One of: monomer, dimer, symmetric, custom
5. chain_length: Suggested protein size (small: 60-80, medium: 80-120, large: 120-160)
6. bury_ligand: Whether to bury the metal/ligand inside the protein (true/false).
   - true (default): creates a deep binding pocket — use for binding, catalysis, most designs
   - false: leaves metal/ligand partially exposed — use for sensing, surface sites, structural roles
7. isomer_specification: If the user specifies cis/trans or E/Z isomer, extract it ("cis", "trans", or null)
8. confidence: Your confidence in the parsing (0.0 to 1.0)

Common metal-ligand combinations:
- Citrate with lanthanides (Tb, Eu, Gd) - luminescent biosensors
- PQQ with calcium - quinoprotein dehydrogenases
- Heme with iron - electron transfer, catalysis
- Zinc fingers - DNA binding

Return ONLY valid JSON with these fields. No explanation text.`;

const PARSER_USER_TEMPLATE = (query: string) =>
  `Parse this protein design request:\n\n"${query}"\n\nReturn JSON with fields: metal_type, ligand_name, isomer_specification, design_goal, target_topology, chain_length_min, chain_length_max, bury_ligand, confidence, reasoning`;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/** Strip markdown code fences if the model wraps its output in them. */
function extractJson(text: string): string {
  const fenced = text.match(/```(?:json)?\s*([\s\S]*?)```/);
  if (fenced) return fenced[1].trim();
  return text.trim();
}

// ---------------------------------------------------------------------------
// Route handler
// ---------------------------------------------------------------------------

export async function POST(request: NextRequest) {
  if (!LLM_API_KEY) {
    return NextResponse.json(
      { error: 'LLM_API_KEY is not configured' },
      { status: 503 },
    );
  }

  let query: string;
  try {
    const body = await request.json();
    // Accept the same shape the frontend already sends via runsync:
    //   { input: { query: "..." } }
    query = body?.input?.query ?? body?.query;
    if (!query || typeof query !== 'string') {
      return NextResponse.json(
        { error: 'Missing `query` in request body' },
        { status: 400 },
      );
    }
  } catch {
    return NextResponse.json(
      { error: 'Invalid JSON body' },
      { status: 400 },
    );
  }

  try {
    // OpenAI-compatible chat completions format
    const llmResponse = await fetch(`${LLM_API_URL}/v1/chat/completions`, {
      method: 'POST',
      headers: {
        'Authorization': `Bearer ${LLM_API_KEY}`,
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        model: LLM_MODEL,
        max_tokens: 4096,
        messages: [
          { role: 'system', content: PARSER_SYSTEM_PROMPT },
          { role: 'user', content: PARSER_USER_TEMPLATE(query) },
        ],
      }),
    });

    if (!llmResponse.ok) {
      const errText = await llmResponse.text().catch(() => '');
      console.error('[ai-parse] LLM API error:', llmResponse.status, errText);
      return NextResponse.json(
        { error: `LLM API error (${llmResponse.status}): ${errText}` },
        { status: 502 },
      );
    }

    const llmData = await llmResponse.json();

    // OpenAI chat completions format:
    //   { choices: [{ message: { content: "..." } }] }
    const rawText = llmData.choices?.[0]?.message?.content ?? '';

    const jsonStr = extractJson(rawText);
    let parsed: Record<string, unknown>;
    try {
      parsed = JSON.parse(jsonStr);
    } catch {
      console.error('[ai-parse] Failed to parse LLM JSON:', jsonStr.substring(0, 300));
      return NextResponse.json(
        {
          error: 'Failed to parse LLM response as JSON',
          raw: jsonStr.substring(0, 500),
        },
        { status: 502 },
      );
    }

    // Augment with metadata the frontend expects
    const result = {
      ...parsed,
      raw_query: query,
      corrected_query: query,
      parser_type: 'ai' as const,
      warnings: (parsed.warnings as string[]) ?? [],
      suggestions: (parsed.suggestions as string[]) ?? [],
      parsed_entities: (parsed.parsed_entities as Record<string, unknown>) ?? {},
      typo_corrections: (parsed.typo_corrections as string[]) ?? [],
    };

    // Wrap in the same envelope the frontend already unpacks:
    //   { output: { status: "completed", result: { ... } } }
    return NextResponse.json({
      output: {
        status: 'completed',
        result,
      },
    });
  } catch (err) {
    console.error('[ai-parse] Unexpected error:', err);
    return NextResponse.json(
      { error: err instanceof Error ? err.message : 'Unknown error' },
      { status: 500 },
    );
  }
}

export const runtime = 'edge';
