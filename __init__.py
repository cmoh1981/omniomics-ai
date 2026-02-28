"""
OmniOmics AI - Multi-omics Analysis Platform
============================================
LLM-powered proteomics and multiomics analysis.
"""

__version__ = "0.1.0"
__author__ = "OmniOmics Team"

from omniomics_ai.modules.llm.orchestrator import OmniLLM
from omniomics_ai.modules.proteomics.analysis import ProteomicsPipeline, ProteomicsData
from omniomics_ai.modules.integration.multiomics import MultiOmicsIntegrator

__all__ = [
    "OmniLLM",
    "ProteomicsPipeline",
    "ProteomicsData",
    "MultiOmicsIntegrator",
]
