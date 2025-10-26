"""
DARLIN Python - Python implementation of CARLIN sequence analysis tools

Main features:
- CRISPR-Cas9 sequence alignment
- Allele calling
- Mutation annotation
- Advanced sequence analysis API

Author: DARLIN-toolkits Team
"""

__version__ = "0.1.0"
__author__ = "Jarning Gau"

# Export main API interfaces
from .api import analyze_sequences, quick_analyze, AnalysisResult

# Export core functionality modules
from .alignment.cas9_align import cas9_align, nt2int, int2nt
from .mutations import Mutation, MutationType, annotate_mutations

# Export configuration
from .config.amplicon_configs import ORIGINAL_CARLIN

__all__ = [
    '__version__',
    '__author__',
    # Main API
    'analyze_sequences',
    'quick_analyze', 
    'AnalysisResult',
    # Core functionality
    'cas9_align',
    'nt2int',
    'int2nt',
    'Mutation',
    'MutationType',
    'annotate_mutations',
    # Configuration
    'ORIGINAL_CARLIN'
] 