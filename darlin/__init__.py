"""
DARLIN Python - Python implementation of CARLIN sequence analysis tools

Main features:
- CRISPR-Cas9 sequence alignment
- Multiple amplicon templates support (Col1a1, Rosa, Tigre)
- Allele calling
- Mutation annotation
- Advanced sequence analysis API

Author: Jarning Gau
"""

__version__ = "0.1.5"
__author__ = "Jarning Gau"

# Export main API interfaces
from .api import analyze_sequences
from .config.amplicon_configs import AmpliconConfig
from .utils.build_config import build_carlin_config

__all__ = [
    '__version__',
    '__author__',
    # Main API
    'analyze_sequences',
    # Create new template config
    'AmpliconConfig',
    'build_carlin_config'
] 