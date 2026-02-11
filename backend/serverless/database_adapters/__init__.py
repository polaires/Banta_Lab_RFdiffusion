"""Database adapters for external structure databases.

Imports are lazy to prevent one broken adapter from cascading failures.
Production code imports directly: `from database_adapters.pubchem_adapter import PubChemAdapter`
"""


def __getattr__(name):
    """Lazy import adapters on first access."""
    if name == "MetalPDBAdapter":
        from .metalpdb_adapter import MetalPDBAdapter
        return MetalPDBAdapter
    elif name == "UniProtAdapter":
        from .uniprot_adapter import UniProtAdapter
        return UniProtAdapter
    elif name == "PubChemAdapter":
        from .pubchem_adapter import PubChemAdapter
        return PubChemAdapter
    elif name == "AlphaFoldAdapter":
        from .alphafold_adapter import AlphaFoldAdapter
        return AlphaFoldAdapter
    raise AttributeError(f"module 'database_adapters' has no attribute {name!r}")


__all__ = ["MetalPDBAdapter", "UniProtAdapter", "PubChemAdapter", "AlphaFoldAdapter"]
