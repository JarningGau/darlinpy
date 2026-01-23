"""
CARLIN allele calling module

This module provides complete allele calling functionality, including:
- Exact and coarse-grained calling algorithms
- Allele calling result data structures
- Statistical analysis functions

@deprecated Some classes in this module (BulkAlleleCallResult, AlleleCallStatistics) 
are not used in the main API and may be removed in future versions.
"""

from .allele_caller import AlleleCaller
from .allele_data import AlleleCallResult

__all__ = [
    'AlleleCaller',
    'AlleleCallResult'
] 