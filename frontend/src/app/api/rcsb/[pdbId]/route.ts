/**
 * RCSB PDB Proxy API Route
 *
 * This route proxies requests to fetch PDB files, avoiding CORS issues.
 * Uses the wwPDB archive which is more reliable than files.rcsb.org/download/
 *
 * Endpoint:
 *   GET /api/rcsb/[pdbId] -> Fetches PDB file from wwPDB archive
 */

import { NextRequest, NextResponse } from 'next/server';
import { gunzipSync } from 'zlib';

export async function GET(request: NextRequest) {
  // Extract PDB ID from pathname
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

  const normalizedPdbId = pdbId.toLowerCase();
  // wwPDB archive organizes by middle 2 characters of PDB ID
  const middleChars = normalizedPdbId.substring(1, 3);

  try {
    // Use wwPDB archive - more reliable than files.rcsb.org/download/
    // Format: /pub/pdb/data/structures/divided/pdb/{middle2}/pdb{pdbid}.ent.gz
    const archiveUrl = `https://files.wwpdb.org/pub/pdb/data/structures/divided/pdb/${middleChars}/pdb${normalizedPdbId}.ent.gz`;

    const response = await fetch(archiveUrl, {
      headers: {
        'Accept': 'application/gzip, */*',
      },
    });

    if (!response.ok) {
      if (response.status === 404) {
        return NextResponse.json(
          { error: `PDB ID "${pdbId.toUpperCase()}" not found in RCSB database` },
          { status: 404 }
        );
      }
      if (response.status === 502 || response.status === 503) {
        return NextResponse.json(
          { error: `RCSB servers are temporarily unavailable. Please try again later.` },
          { status: 503 }
        );
      }
      return NextResponse.json(
        { error: `Failed to fetch PDB: ${response.statusText}` },
        { status: response.status }
      );
    }

    // Decompress the gzipped content
    const compressedBuffer = await response.arrayBuffer();
    const decompressed = gunzipSync(Buffer.from(compressedBuffer));
    const pdbContent = decompressed.toString('utf-8');

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

// Use Node.js runtime (required for zlib)
export const runtime = 'nodejs';
