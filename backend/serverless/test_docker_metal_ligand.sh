#!/bin/bash
# Test script for metal-ligand complex dimer design in Docker
#
# Usage:
#   1. Start Docker container: docker-compose -f docker-compose.local.yml up --build
#   2. Run this script: bash test_docker_metal_ligand.sh
#
# Prerequisites:
#   - Docker container running on port 8000
#   - curl installed

set -e

API_URL="http://localhost:8000/runsync"
CONTENT_TYPE="Content-Type: application/json"

echo "============================================"
echo "Metal-Ligand Complex Dimer Design Tests"
echo "============================================"
echo ""

# Test 1: Health check
echo "Test 1: Health check..."
curl -s -X POST "$API_URL" \
  -H "$CONTENT_TYPE" \
  -d '{"input": {"task": "health"}}' | python -m json.tool
echo ""

# Test 2: List available templates
echo "Test 2: List available metal-ligand templates..."
curl -s -X POST "$API_URL" \
  -H "$CONTENT_TYPE" \
  -d '{"input": {"task": "interface_metal_ligand_design", "list_templates": true}}' | python -m json.tool
echo ""

# Test 3: PQQ-Ca template info
echo "Test 3: Get PQQ-Ca template info..."
curl -s -X POST "$API_URL" \
  -H "$CONTENT_TYPE" \
  -d '{"input": {"task": "interface_metal_ligand_design", "template_name": "pqq_ca", "info_only": true}}' | python -m json.tool
echo ""

# Test 4: Citrate-Tb template info
echo "Test 4: Get Citrate-Tb template info..."
curl -s -X POST "$API_URL" \
  -H "$CONTENT_TYPE" \
  -d '{"input": {"task": "interface_metal_ligand_design", "template_name": "citrate_tb", "info_only": true}}' | python -m json.tool
echo ""

# Test 5: PQQ-Ca mock design (no GPU required)
echo "Test 5: PQQ-Ca mock design..."
curl -s -X POST "$API_URL" \
  -H "$CONTENT_TYPE" \
  -d '{
    "input": {
      "task": "interface_metal_ligand_design",
      "template_name": "pqq_ca",
      "approach": "joint",
      "chain_length": "60-80",
      "num_designs": 1,
      "use_mock": true
    }
  }' | python -m json.tool
echo ""

# Test 6: Citrate-Tb mock design (no GPU required)
echo "Test 6: Citrate-Tb mock design..."
curl -s -X POST "$API_URL" \
  -H "$CONTENT_TYPE" \
  -d '{
    "input": {
      "task": "interface_metal_ligand_design",
      "template_name": "citrate_tb",
      "approach": "joint",
      "chain_length": "60-80",
      "num_designs": 1,
      "chain_a_donors": ["Glu", "Asp"],
      "chain_b_donors": ["Glu", "Asp", "Asn"],
      "use_mock": true
    }
  }' | python -m json.tool
echo ""

# Test 7: Error handling - unknown template
echo "Test 7: Error handling - unknown template..."
curl -s -X POST "$API_URL" \
  -H "$CONTENT_TYPE" \
  -d '{
    "input": {
      "task": "interface_metal_ligand_design",
      "template_name": "unknown_template"
    }
  }' | python -m json.tool
echo ""

# Test 8: Error handling - missing required input
echo "Test 8: Error handling - missing required input..."
curl -s -X POST "$API_URL" \
  -H "$CONTENT_TYPE" \
  -d '{
    "input": {
      "task": "interface_metal_ligand_design"
    }
  }' | python -m json.tool
echo ""

# Test 9: Custom complex PDB (mock)
echo "Test 9: Custom complex PDB mock design..."
curl -s -X POST "$API_URL" \
  -H "$CONTENT_TYPE" \
  -d '{
    "input": {
      "task": "interface_metal_ligand_design",
      "complex_pdb": "HETATM    1  C1  PQQ L   1      48.000  50.000  48.000  1.00  0.00           C\nHETATM    2  O5  PQQ L   1      50.500  51.000  50.000  1.00  0.00           O\nHETATM    3  N6  PQQ L   1      49.500  50.500  50.500  1.00  0.00           N\nHETATM    4  O7  PQQ L   1      50.000  49.000  50.000  1.00  0.00           O\nHETATM   10 CA   CA  M   2      50.000  50.000  52.400  1.00  0.00          CA\nEND",
      "approach": "joint",
      "chain_length": "60-80",
      "num_designs": 1,
      "use_mock": true
    }
  }' | python -m json.tool
echo ""

echo "============================================"
echo "All tests completed!"
echo "============================================"

# Optional: Run pytest inside container
echo ""
echo "To run unit tests inside Docker, execute:"
echo "  docker exec -it \$(docker ps -q -f name=rfdiffusion) pytest /app/test_metal_ligand_complex.py -v"
