# Import base class from sequana_report (Phase 2: refactor)
from sequana_report import SequanaBaseModule  # noqa: F401

from .panther_enrichment import ModulePantherEnrichment
from .kegg_enrichment import ModuleKEGGEnrichment
