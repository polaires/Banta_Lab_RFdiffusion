"""
Extract PDB content from Docker logs for completed r7v5 jobs.
This is a workaround to recover the designs that were generated.
"""

import subprocess
import json
import re
from pathlib import Path

OUTPUT_DIR = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\experiments\ln_citrate_scaffold\outputs\round_07_scaffold")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Get docker logs
cmd = 'wsl -d Ubuntu -e bash -c "docker compose -f /mnt/g/Github_local_repo/Banta_Lab_RFdiffusion/backend/serverless/docker-compose.local.yml logs --tail=2000"'
result = subprocess.run(cmd, shell=True, capture_output=True, text=True, encoding='utf-8', errors='replace')

logs = result.stdout + result.stderr

# Find completed RFD3 jobs with PDB content
# Pattern: run_job return: {'output': {'status': 'completed', 'result': {'designs': [...]
completed_pattern = r"run_job return: \{\'output\': \{\'status\': \'completed\', \'result\': \{\'designs\': \[(\{.*?\})\]"

# Simpler approach: Find lines with PDB ATOM content
pdb_sections = []
current_pdb = []
in_pdb = False

for line in logs.split('\n'):
    if "'content': 'ATOM" in line or "'content': 'HETATM" in line:
        # Start of PDB content
        start_idx = line.find("'content': '") + len("'content': '")
        # Extract content - might be truncated
        content_part = line[start_idx:]
        if content_part.endswith("'"):
            content_part = content_part[:-1]
        pdb_sections.append(content_part)

print(f"Found {len(pdb_sections)} PDB sections in logs")

# The logs are truncated, so we can't extract full PDBs from logs
# Instead, let's make a simple test API call to get ONE design
print("\nThe completed designs were generated but the save function had a bug.")
print("The status is 'COMPLETED' (uppercase) in outer response, 'completed' in inner output.")
print("\nLet me create a corrected save script...")
