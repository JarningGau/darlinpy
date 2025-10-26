"""
CARLIN allele calling module

This module provides complete allele calling functionality, including:
- Exact and coarse-grained calling algorithms
- Allele calling result data structures
- Statistical analysis functions
"""

from .allele_caller import AlleleCaller
from .allele_data import (
    AlleleCallResult, 
    BulkAlleleCallResult, 
    AlleleCallStatistics
)

__all__ = [
    'AlleleCaller',
    'AlleleCallResult',
    'BulkAlleleCallResult', 
    'AlleleCallStatistics'
] 