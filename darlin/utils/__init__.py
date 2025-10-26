"""
Utility module - Provides basic sequence processing functions required for CARLIN analysis
"""

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