"""
Utility module - Provides basic sequence processing functions required for CARLIN analysis

@deprecated This module is not used in the main API (analyze_sequences, AmpliconConfig, build_carlin_config).
It may be removed in future versions.
"""

import warnings

warnings.warn(
    "darlin.utils module is deprecated and not used in the main API. "
    "It may be removed in future versions.",
    DeprecationWarning,
    stacklevel=2
)

from .sequence_utils import (
    reverse_complement,
    is_valid_dna,
    clean_sequence,
    calculate_gc_content,
    count_nucleotides,
    hamming_distance,
    format_sequence,
    parse_fasta_string,
    to_fasta_string
)

__all__ = [
    'reverse_complement',
    'is_valid_dna',
    'clean_sequence',
    'calculate_gc_content',
    'count_nucleotides',
    'hamming_distance',
    'format_sequence',
    'parse_fasta_string',
    'to_fasta_string'
] 