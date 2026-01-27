# Utils package
from .sasa import calculate_sasa, classify_burial
from .dssp import (
    assign_secondary_structure,
    get_ss_codes_by_residue,
    extract_ss_segments,
    get_ss_context,
    suggest_ss_constraints,
)

# Conservation analysis utilities (ConSurf methodology)
# These may not be available if dependencies are missing
try:
    from .blast_client import (
        NCBIBlastClient,
        BlastHit,
        BlastResult,
        extract_sequence_from_pdb,
        fetch_sequences,
    )
    BLAST_CLIENT_AVAILABLE = True
except ImportError:
    BLAST_CLIENT_AVAILABLE = False

try:
    from .msa_utils import (
        filter_blast_hits,
        cluster_sequences_mmseqs,
        run_muscle_alignment,
        prepare_msa_from_blast,
        MSAResult,
    )
    MSA_UTILS_AVAILABLE = True
except ImportError:
    MSA_UTILS_AVAILABLE = False

try:
    from .rate4site_wrapper import (
        run_rate4site,
        parse_rate4site_output,
        assign_consurf_grades,
        calculate_simple_conservation,
        Rate4SiteScore,
        Rate4SiteResult,
        RATE4SITE_AVAILABLE,
    )
except ImportError:
    RATE4SITE_AVAILABLE = False
