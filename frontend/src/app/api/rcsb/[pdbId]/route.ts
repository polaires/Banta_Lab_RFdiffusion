/**
 * RCSB PDB Proxy API Route
 *
 * This route proxies requests to RCSB models.rcsb.org to avoid CORS issues
 * in environments where direct access might be blocked.
 *
 * Endpoint:
 *   GET /api/rcsb/[pdbId] -> Fetches mmCIF file from models.rcsb.org
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
    // Use models.rcsb.org v1 API which is reliable and has CORS support
    const rcsbUrl = `https://models.rcsb.org/v1/${normalizedPdbId}/full`;

    const response = await fetch(rcsbUrl, {
      headers: {
        'Accept': 'text/plain, */*',
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
        { error: `Failed to fetch structure: ${response.statusText}` },
        { status: response.status }
      );
    }

    const cifContent = await response.text();

    // Return the CIF content as plain text
    return new NextResponse(cifContent, {
      status: 200,
      headers: {
        'Content-Type': 'text/plain',
        'Cache-Control': 'public, max-age=86400', // Cache for 24 hours
      },
    });
  } catch (error) {
    console.error('[RCSB Proxy] Error fetching structure:', error);
    return NextResponse.json(
      { error: error instanceof Error ? error.message : 'Failed to fetch structure from RCSB' },
      { status: 500 }
    );
  }
}

// Use Node.js runtime for reliability
export const runtime = 'nodejs';
