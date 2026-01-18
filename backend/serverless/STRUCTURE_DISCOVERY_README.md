# AI-Driven Structure Discovery Infrastructure

A comprehensive infrastructure for AI-assisted protein structure discovery and metal coordination chemistry queries. This system provides natural language interfaces for querying structural databases, extracting coordination information, and generating design recommendations for protein engineering.

## Overview

### Purpose

The AI-Driven Structure Discovery Infrastructure enables:

- **Natural language queries** about metal coordination chemistry
- **Multi-database aggregation** from MetalPDB, RCSB PDB, UniProt, and AlphaFold
- **Intelligent template selection** for protein design workflows
- **Design recommendations** with amino acid biasing for sequence design

### Architecture

```
+------------------------------------------+
|           AI Interface Layer              |
|   (ai_structure_interface.py)             |
|   - answer_structure_question()           |
|   - plan_structure_search()               |
|   - get_design_recommendations()          |
+------------------------------------------+
                    |
                    v
+------------------------------------------+
|      Structure Discovery Service          |
|      (structure_discovery.py)             |
|   - search_by_intent()                    |
|   - search_by_metal()                     |
|   - get_best_template()                   |
|   - extract_coordination_info()           |
+------------------------------------------+
                    |
                    v
+------------------------------------------+
|         Database Adapters                 |
|   (database_adapters/)                    |
|   - MetalPDBAdapter   - UniProtAdapter    |
|   - PubChemAdapter    - AlphaFoldAdapter  |
+------------------------------------------+
```

### Components

| Component | Location | Purpose |
|-----------|----------|---------|
| AI Interface | `ai_structure_interface.py` | Natural language Q&A entry point |
| Structure Discovery | `structure_discovery.py` | Multi-database search aggregation |
| MetalPDB Adapter | `database_adapters/metalpdb_adapter.py` | Metal coordination site database |
| UniProt Adapter | `database_adapters/uniprot_adapter.py` | Protein function and sequence data |
| PubChem Adapter | `database_adapters/pubchem_adapter.py` | Ligand SMILES and chemical data |
| AlphaFold Adapter | `database_adapters/alphafold_adapter.py` | Predicted structure access |

---

## Quick Start

### Basic Usage

```python
from structure_discovery import StructureDiscovery
from ai_structure_interface import (
    answer_structure_question,
    plan_structure_search,
    get_design_recommendations,
)

# Natural language query about metal coordination
result = answer_structure_question("What is the coordination geometry of terbium?")
print(result["answer"])
# Output: "Terbium (TB3+) typically has coordination number 8-9 with
#          tricapped_trigonal_prismatic geometry..."

# Get design recommendations for a biosensor
recs = get_design_recommendations(metal="TB", task="biosensor")
print(recs["donors"])           # ['Glu', 'Asp', 'Asn', 'Gln', 'His']
print(recs["bias_aa"])          # E:3.0,D:3.0,N:2.4,Q:2.4,H:0.5,C:-5.0,A:-2.0
print(recs["coordination_info"]) # {'coordination_number_range': (8, 9), ...}

# Search for structures by natural language
discovery = StructureDiscovery()
results = discovery.search_by_intent("terbium binding protein", limit=10)
for r in results[:3]:
    print(f"{r.pdb_id}: {r.metal} - {r.geometry}")
```

### Workflow Example

```python
# Complete workflow: Plan -> Search -> Design

# 1. Plan the search
plan = plan_structure_search("I need to design a terbium citrate biosensor")
print(f"Metal: {plan['metal']}")           # TB
print(f"Ligand: {plan['ligand']}")         # CIT
print(f"Strategy: {plan['search_strategy']}")  # lanthanide_specialized

# 2. Search for structures
discovery = StructureDiscovery()
results = discovery.search_by_metal("TB", limit=5)

# 3. Get best template
template = discovery.get_best_template("TB", ligand="CIT")
print(f"Template source: {template['source']}")
print(f"Coordination: {template['coordination_number']}")

# 4. Extract coordination info from a structure
coord_info = discovery.extract_coordination_info("6MI5", "TB")
print(f"Geometry: {coord_info['geometry']}")
print(f"Coordinating residues: {coord_info['coordinating_residues']}")

# 5. Get design recommendations
recs = get_design_recommendations("TB", task="biosensor", ligand="CIT")
print(f"Use bias string: {recs['bias_aa']}")
```

---

## Database Adapters

### MetalPDBAdapter

Queries metal coordination sites from the MetalPDB database with automatic fallback to RCSB PDB.

```python
from database_adapters import MetalPDBAdapter

adapter = MetalPDBAdapter(timeout=30, enable_fallback=True)

# Search by metal type
results = adapter.search_by_metal(
    metal="ZN",
    coordination_number=4,  # Optional filter
    geometry="tetrahedral", # Optional filter
    limit=20
)

# Get site details
details = adapter.get_site_details("1CA2", "ZN", site_index=0)
print(f"Coordination: {details.coordination_number}")
print(f"Geometry: {details.geometry}")
print(f"Residues: {details.coordinating_residues}")
```

**Returns:** `MetalSiteResult` dataclass with:
- `pdb_id`: PDB identifier
- `metal`: Metal symbol
- `coordination_number`: Number of coordinating atoms
- `geometry`: Coordination geometry (tetrahedral, octahedral, etc.)
- `resolution`: Structure resolution
- `coordinating_residues`: List of residue strings
- `distances`: Metal-ligand distances

### UniProtAdapter

Queries protein function and metal binding annotations from UniProt.

```python
from database_adapters import UniProtAdapter

adapter = UniProtAdapter(timeout=30)

# Search by protein function
proteins = adapter.search_by_function(
    function_query="zinc finger",
    organism="human",
    reviewed_only=True,
    limit=10
)

# Search by metal binding
binding_proteins = adapter.search_by_metal_binding(metal="calcium", limit=20)

# Get PDB mappings for a protein
pdb_mappings = adapter.get_pdb_mappings("P00918")  # Carbonic anhydrase 2
for pdb in pdb_mappings:
    print(f"PDB: {pdb['pdb_id']}, Chains: {pdb['chains']}")

# Get protein sequence
sequence = adapter.get_sequence("P00918")
```

### PubChemAdapter

Queries ligand chemical information including SMILES structures.

```python
from database_adapters import PubChemAdapter

adapter = PubChemAdapter(timeout=30)

# Search compounds by name
compounds = adapter.search_by_name("citric acid", limit=5)
for c in compounds:
    print(f"CID: {c['cid']}, SMILES: {c['smiles']}")

# Get SMILES for common ligands (fast lookup)
atp_smiles = adapter.get_smiles("ATP")
citrate_smiles = adapter.get_smiles("citrate")

# Get compound by CID
compound = adapter.get_compound_by_cid(311)  # Citrate
print(f"Formula: {compound['formula']}")
print(f"Weight: {compound['molecular_weight']}")
```

### AlphaFoldAdapter

Accesses predicted protein structures from the AlphaFold Database.

```python
from database_adapters import AlphaFoldAdapter

adapter = AlphaFoldAdapter(timeout=30)

# Check if prediction exists
prediction = adapter.get_prediction("P00918")
if prediction:
    print(f"Model version: {prediction['model_version']}")
    print(f"Mean pLDDT: {prediction['plddt']}")

# Download structure
structure_pdb = adapter.download_structure("P00918", format="pdb")
structure_cif = adapter.download_structure("P00918", format="cif")

# Get confidence scores
confidence = adapter.get_confidence_scores("P00918")
```

---

## StructureDiscovery Service

The main service class that aggregates results from all database adapters.

### Methods

#### `search_by_intent(query, limit=20)`

Search using natural language intent parsing.

```python
discovery = StructureDiscovery()

# Parses query to extract metal, ligand, and protein type
results = discovery.search_by_intent("PQQ dehydrogenase", limit=10)
results = discovery.search_by_intent("terbium binding protein", limit=10)
results = discovery.search_by_intent("zinc finger transcription factor", limit=10)
```

**Recognized patterns:**
- Metal names: zinc, terbium, calcium, iron, copper, etc.
- Ligand names: PQQ, citrate, ATP, NAD, heme, etc.
- Protein types: dehydrogenase, zinc finger, EF-hand, carbonic anhydrase, etc.

#### `search_by_metal(metal, coordination_number=None, geometry=None, limit=20)`

Search specifically by metal type.

```python
# Basic metal search
results = discovery.search_by_metal("ZN", limit=20)

# With filters
results = discovery.search_by_metal(
    metal="ZN",
    coordination_number=4,
    geometry="tetrahedral",
    limit=10
)
```

**Priority order:**
1. Curated reference structures (high-quality templates)
2. MetalPDB (specialized coordination database)
3. RCSB PDB (general structure database)
4. UniProt (metal binding annotations)

#### `get_best_template(metal, ligand=None, use_architector=True)`

Get the best template for a metal-ligand combination.

```python
# Get template for terbium-citrate
template = discovery.get_best_template("TB", ligand="CIT")

# Template contains:
# - pdb_id: Source structure
# - pdb_content: PDB file content (if available)
# - coordination_number: Expected CN
# - geometry: Expected geometry
# - coordinating_residues: Residue info
# - source: "curated", "metalpdb", "rcsb", or "architector"
```

#### `extract_coordination_info(pdb_id, metal, site_index=0)`

Extract detailed coordination information from a specific structure.

```python
info = discovery.extract_coordination_info("1CA2", "ZN")

print(f"PDB: {info['pdb_id']}")
print(f"Metal: {info['metal']}")
print(f"CN: {info['coordination_number']}")
print(f"Geometry: {info['geometry']}")
print(f"Residues: {info['coordinating_residues']}")
print(f"Distances: {info['distances']}")
print(f"Coords: {info['metal_coords']}")
```

---

## AI Interface Functions

### `answer_structure_question(question)`

Main entry point for natural language questions about metal coordination.

```python
from ai_structure_interface import answer_structure_question

# Coordination questions
result = answer_structure_question("What is the coordination geometry of terbium?")
print(result["answer"])
print(result["coordination_number"])  # [8, 9]
print(result["geometries"])           # ['tricapped_trigonal_prismatic']

# Find structure questions
result = answer_structure_question("Find me structures with calcium and PQQ")
print(result["structures"])           # List of structure dicts
print(result["count"])                # Number found

# Donor questions
result = answer_structure_question("What are the preferred donors for zinc?")
print(result["donor_weights"])        # {'N': 3.0, 'O': 2.5, 'S': 2.0}
print(result["common_residues"])      # ['His', 'Cys', 'Asp', 'Glu']

# Template questions
result = answer_structure_question("Get a template for terbium biosensor design")
print(result["template_name"])
print(result["template_info"])
```

**Question types detected:**
- Coordination: "coordination geometry", "coordination number", "CN of"
- Find structure: "find", "search", "show me structures"
- Template: "template", "starting point", "scaffold"
- Donors: "preferred donors", "which residues coordinate"

### `plan_structure_search(task_description)`

Plan a search strategy based on a task description.

```python
from ai_structure_interface import plan_structure_search

plan = plan_structure_search("I need to design a terbium biosensor")

print(plan["steps"])           # ['Search for terbium binding proteins', ...]
print(plan["metal"])           # 'TB'
print(plan["ligand"])          # None or detected ligand
print(plan["task_type"])       # 'biosensor'
print(plan["search_strategy"]) # 'lanthanide_specialized'
print(plan["databases"])       # ['metalpdb', 'rcsb', 'uniprot']
```

### `get_design_recommendations(metal, task=None, ligand=None)`

Get comprehensive design recommendations for a metal binding site.

```python
from ai_structure_interface import get_design_recommendations

recs = get_design_recommendations(metal="TB", task="biosensor")

# Key outputs:
print(recs["donors"])              # ['Glu', 'Asp', 'Asn', 'Gln', 'His']
print(recs["bias_aa"])             # LigandMPNN bias string
print(recs["avoid_residues"])      # ['Cys', 'Met'] for hard acids
print(recs["hsab_class"])          # 'hard' / 'soft' / 'borderline'
print(recs["coordination_info"])   # {'coordination_number_range': (8, 9), ...}
print(recs["templates"])           # Available templates
print(recs["task_recommendations"]) # Task-specific advice
```

**Bias string format:** `"E:3.0,D:3.0,N:2.4,Q:2.4,H:0.5,C:-5.0,A:-2.0"`

---

## Integration with Task Planner

The structure discovery system integrates with the AI Task Planner for automated workflows.

### DISCOVER_STRUCTURES Step Type

```python
from backend.ai.task_planner import TaskStep, TaskStepType

# Create a discovery step
step = TaskStep(
    id="discover_1",
    type=TaskStepType.DISCOVER_STRUCTURES,
    name="Find terbium templates",
    description="Search for terbium binding sites",
    params={
        "query": "terbium binding protein",
        "metal": "TB",
        "ligand": "CIT",
        "limit": 10,
    }
)
```

### Enhanced FETCH_PDB with Auto-Discovery

When a PDB ID is not specified, FETCH_PDB can auto-discover structures:

```python
step = TaskStep(
    id="fetch_1",
    type=TaskStepType.FETCH_PDB,
    name="Get terbium template",
    description="Auto-discover and fetch terbium structure",
    params={
        "metal": "TB",
        "auto_discover": True,  # Enable auto-discovery
        "prefer_curated": True, # Prefer curated templates
    }
)
```

### Example Task Plan

```python
from backend.ai.task_planner import TaskPlan, TaskStep, TaskStepType

plan = TaskPlan(
    id="tb_biosensor_001",
    name="Terbium Biosensor Design",
    goal="Design a terbium-citrate biosensor",
    description="Multi-step protein design for Tb(III) sensing",
    steps=[
        TaskStep(
            id="step1",
            type=TaskStepType.DISCOVER_STRUCTURES,
            name="Discover templates",
            description="Find terbium binding protein templates",
            params={"metal": "TB", "ligand": "CIT", "limit": 5}
        ),
        TaskStep(
            id="step2",
            type=TaskStepType.FETCH_PDB,
            name="Fetch best template",
            description="Download the best discovered structure",
            params={"use_discovery_result": True},
            depends_on=["step1"]
        ),
        TaskStep(
            id="step3",
            type=TaskStepType.RFD3_DESIGN,
            name="Scaffold generation",
            description="Generate scaffold around metal site",
            params={"use_template_coords": True},
            depends_on=["step2"]
        ),
        TaskStep(
            id="step4",
            type=TaskStepType.MPNN_DESIGN,
            name="Sequence design",
            description="Design sequences with metal-biased AA",
            params={"use_metal_bias": True},
            depends_on=["step3"]
        ),
    ]
)
```

---

## Testing

### Run All Tests

```bash
# Navigate to serverless directory
cd backend/serverless

# Run all unit tests
pytest tests/ -v

# Run with coverage
pytest tests/ -v --cov=. --cov-report=html
```

### Run Integration Tests

Integration tests require network access to external databases.

```bash
# Run integration tests only
pytest tests/test_integration.py -v -m integration

# Skip network-dependent tests
pytest tests/ -v -m "not integration"
```

### Test Categories

| Test File | Purpose |
|-----------|---------|
| `tests/test_metalpdb_adapter.py` | MetalPDB adapter unit tests |
| `tests/test_uniprot_adapter.py` | UniProt adapter unit tests |
| `tests/test_pubchem_adapter.py` | PubChem adapter unit tests |
| `tests/test_alphafold_adapter.py` | AlphaFold adapter unit tests |
| `tests/test_structure_discovery.py` | StructureDiscovery service tests |
| `tests/test_ai_interface.py` | AI interface function tests |
| `tests/test_coordination_knowledge.py` | Coordination chemistry tests |
| `tests/test_integration.py` | End-to-end integration tests |

### Example Test Run

```bash
# Test a specific adapter
pytest tests/test_metalpdb_adapter.py -v

# Test AI interface
pytest tests/test_ai_interface.py -v

# Run with verbose output and stop on first failure
pytest tests/ -v -x

# Run specific test class
pytest tests/test_structure_discovery.py::TestSearchByIntent -v
```

---

## Configuration

### Environment Variables

| Variable | Default | Purpose |
|----------|---------|---------|
| `METALPDB_TIMEOUT` | 30 | MetalPDB request timeout (seconds) |
| `UNIPROT_TIMEOUT` | 30 | UniProt request timeout (seconds) |
| `PUBCHEM_TIMEOUT` | 30 | PubChem request timeout (seconds) |
| `ALPHAFOLD_TIMEOUT` | 30 | AlphaFold request timeout (seconds) |
| `ENABLE_FALLBACK` | True | Enable fallback between databases |

### Initialization Options

```python
from structure_discovery import StructureDiscovery

# Custom configuration
discovery = StructureDiscovery(
    timeout=60,           # Longer timeout for slow connections
    enable_fallback=True, # Enable RCSB fallback when MetalPDB fails
)
```

---

## Supported Metals

The system includes comprehensive data for:

**Transition Metals:** ZN, FE, CU, MN, CO, NI, MO, W, V, CR, CD, HG, PT, AU, AG

**Alkaline Earth:** CA, MG, BA, SR

**Lanthanides:** TB, EU, GD, LA, CE, PR, ND, SM, DY, HO, ER, TM, YB, LU

**Alkali:** NA, K, LI

Each metal includes:
- HSAB classification (hard/soft/borderline)
- Coordination number ranges by oxidation state
- Preferred donor atoms
- Common coordinating residues
- Example PDB structures

---

## Error Handling

```python
from structure_discovery import StructureDiscovery
from ai_structure_interface import answer_structure_question

# Handle unknown metals
result = answer_structure_question("What is the coordination of unobtainium?")
if result.get("error") == "unknown_metal":
    print("Metal not in database")
    print(f"Valid metals: {result['valid_metals']}")

# Handle network errors gracefully
discovery = StructureDiscovery(enable_fallback=True)
try:
    results = discovery.search_by_metal("ZN")
except Exception as e:
    print(f"Search failed: {e}")
    # Fallback behavior is automatic if enabled
```

---

## Dependencies

**Required:**
- Python 3.8+
- `requests` (for HTTP API calls)

**Optional:**
- `pytest` (for testing)
- `pytest-cov` (for coverage)

**Related Modules:**
- `metal_chemistry.py`: Metal coordination database
- `metal_ligand_templates.py`: Template library
- `metal_site_fetcher.py`: RCSB PDB queries

---

## Contributing

1. Add tests for any new functionality
2. Follow existing code style (PEP 8)
3. Update documentation for API changes
4. Run `pytest tests/ -v` before submitting

---

## License

Part of the Banta Lab RFdiffusion project. Internal use only.
