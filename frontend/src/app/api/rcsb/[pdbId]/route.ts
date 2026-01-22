/**
 * RCSB PDB Proxy API Route
 *
 * This route proxies requests to RCSB to avoid CORS issues.
 * The RCSB server doesn't include CORS headers, so we need to
 * fetch from the server side and return the content.
 *
 * Endpoint:
 *   GET /api/rcsb/[pdbId] -> Fetches PDB file from RCSB
 */

import { NextRequest, NextResponse } from 'next/server';

export async function GET(request: NextRequest) {
  // Extract PDB ID from pathname (same pattern as runpod route)
  const pathname = request.nextUrl.pathname;
  const pdbIdMatch = pathname.match(/\/api\/rcsb\/([^/]+)/);
  const pdbId = pdbIdMatch ? pdbIdMatch[1] : null;

  // Validate PDB ID format (4 characters, alphanumeric)
  if (!pdbId || !/^[A-Za-z0-9]{4}$/.test(pdbId)) {
    return NextResponse.json(
      { error: 'Invalid PDB ID. Must be exactly 4 alphanumeric characters.' },
      { status: 400 }
    );
  }

  const normalizedPdbId = pdbId.toUpperCase();

  try {
    const rcsbUrl = `https://files.rcsb.org/download/${normalizedPdbId}.pdb`;

    const response = await fetch(rcsbUrl, {
      headers: {
        'User-Agent': 'Mozilla/5.0 (compatible; RFdiffusion/1.0)',
        'Accept': 'text/plain, chemical/x-pdb, */*',
      },
    });

    if (!response.ok) {
      if (response.status === 404) {
        return NextResponse.json(
          { error: `PDB ID "${normalizedPdbId}" not found in RCSB database` },
          { status: 404 }
        );
      }
      return NextResponse.json(
        { error: `Failed to fetch PDB: ${response.statusText}` },
        { status: response.status }
      );
    }

    const pdbContent = await response.text();

    // Return the PDB content as plain text
    return new NextResponse(pdbContent, {
      status: 200,
      headers: {
        'Content-Type': 'text/plain',
        'Cache-Control': 'public, max-age=86400', // Cache for 24 hours
      },
    });
  } catch (error) {
    console.error('[RCSB Proxy] Error fetching PDB:', error);
    return NextResponse.json(
      { error: error instanceof Error ? error.message : 'Failed to fetch PDB from RCSB' },
      { status: 500 }
    );
  }
}

// Use Node.js runtime instead of Edge - Edge has issues with some external fetches
export const runtime = 'nodejs';
