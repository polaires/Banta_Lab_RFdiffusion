#!/usr/bin/env python3
"""Clean a PDB file: remove hydrogens, handle alternate conformations."""

import sys

seen = set()
for line in sys.stdin:
    if line.startswith('ATOM'):
        chain = line[21]
        if chain != 'A':
            continue
        # Skip hydrogens
        atom_name = line[12:16].strip()
        if atom_name.startswith('H') or (len(atom_name) > 1 and atom_name[0].isdigit() and atom_name[1] == 'H'):
            continue
        # Handle alternate conformations - take first one only
        alt_loc = line[16]
        res_num = line[22:26].strip()
        key = (res_num, atom_name)
        if key in seen:
            continue
        seen.add(key)
        # Remove alternate location indicator
        if alt_loc != ' ':
            line = line[:16] + ' ' + line[17:]
        print(line, end='')
print('END')
