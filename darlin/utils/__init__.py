"""
实用工具模块 - 提供CARLIN分析所需的基础序列处理功能
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