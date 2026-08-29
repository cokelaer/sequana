# Import base class from sequana_report (Phase 2: refactor)
from sequana_report import SequanaBaseModule  # noqa: F401

from .kegg_enrichment import ModuleKEGGEnrichment  # noqa: F401
from .panther_enrichment import ModulePantherEnrichment  # noqa: F401

__all__ = [
    "SequanaBaseModule",
    "ModulePantherEnrichment",
    "ModuleKEGGEnrichment",
]
