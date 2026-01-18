"""Database adapters for external structure databases."""
from .metalpdb_adapter import MetalPDBAdapter
from .uniprot_adapter import UniProtAdapter

__all__ = ["MetalPDBAdapter", "UniProtAdapter"]
