"""
CARLIN allele calling module

This module provides complete allele calling functionality, including:
- Exact and coarse-grained calling algorithms
- Allele calling result data structures
- Statistical analysis functions

@deprecated Some classes in this module (BulkAlleleCallResult, AlleleCallStatistics) 
are not used in the main API and may be removed in future versions.
"""

import warnings

from .allele_caller import AlleleCaller
from .allele_data import (
    AlleleCallResult, 
    BulkAlleleCallResult, 
    AlleleCallStatistics
)

# Warn about deprecated classes
warnings.warn(
    "BulkAlleleCallResult and AlleleCallStatistics are deprecated and not used in the main API. "
    "They may be removed in future versions.",
    DeprecationWarning,
    stacklevel=2
)

__all__ = [
    'AlleleCaller',
    'AlleleCallResult',
    'BulkAlleleCallResult', 
    'AlleleCallStatistics'
] 