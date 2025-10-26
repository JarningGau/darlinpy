#!/usr/bin/env python3
"""
Sequence processing utility functions

Provides basic DNA sequence processing functionality required for CARLIN analysis
"""

from typing import List, Tuple, Dict
from collections import Counter


# DNA sequence related constants
DNA_BASES = set('ATCG')
DNA_BASES_AMBIGUOUS = set('ATCGRYSWKMBDHVN')
COMPLEMENT_MAP = {
    'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
    'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
    'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
    'H': 'D', 'V': 'B', 'N': 'N'
}


def reverse_complement(seq: str) -> str:
    """
    Calculate reverse complement of DNA sequence
    
    Args:
        seq: DNA sequence string
        
    Returns:
        str: Reverse complement sequence
        
    Examples:
        >>> reverse_complement("ATCG")
        'CGAT'
    """
    if not seq:
        return ""
    
    # Convert to uppercase and validate
    seq = seq.upper()
    if not is_valid_dna(seq, allow_ambiguous=True):
        raise ValueError(f"Invalid DNA sequence: {seq}")
    
    # Generate reverse complement sequence
    complement = ''.join(COMPLEMENT_MAP.get(base, 'N') for base in seq)
    return complement[::-1]


def is_valid_dna(seq: str, allow_ambiguous: bool = False) -> bool:
    """
    Validate if sequence is valid DNA
    
    Args:
        seq: Sequence to validate
        allow_ambiguous: Whether to allow ambiguous bases
        
    Returns:
        bool: Whether sequence is valid DNA
    """
    if not seq:
        return True
    
    seq = seq.upper()
    allowed_bases = DNA_BASES_AMBIGUOUS if allow_ambiguous else DNA_BASES
    
    return all(base in allowed_bases for base in seq)


def clean_sequence(seq: str, remove_chars: str = "- \n\t") -> str:
    """
    Clean sequence by removing specified characters
    
    Args:
        seq: Original sequence
        remove_chars: Characters to remove
        
    Returns:
        str: Cleaned sequence
    """
    if not seq:
        return ""
    
    cleaned = seq.upper()
    for char in remove_chars:
        cleaned = cleaned.replace(char, '')
    
    return cleaned


def calculate_gc_content(seq: str) -> float:
    """
    Calculate GC content
    
    Args:
        seq: DNA sequence
        
    Returns:
        float: GC content percentage (0-100)
    """
    if not seq:
        return 0.0
    
    seq = clean_sequence(seq)
    if not is_valid_dna(seq):
        raise ValueError(f"Invalid DNA sequence: {seq}")
    
    gc_count = seq.count('G') + seq.count('C')
    return (gc_count / len(seq)) * 100 if len(seq) > 0 else 0.0


def count_nucleotides(seq: str) -> Dict[str, int]:
    """
    Count nucleotide frequencies
    
    Args:
        seq: DNA sequence
        
    Returns:
        Dict[str, int]: Count of each nucleotide
    """
    if not seq:
        return {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    
    seq = clean_sequence(seq)
    if not is_valid_dna(seq):
        raise ValueError(f"Invalid DNA sequence: {seq}")
    
    counts = Counter(seq)
    return {
        'A': counts.get('A', 0),
        'T': counts.get('T', 0),
        'G': counts.get('G', 0),
        'C': counts.get('C', 0)
    }


def hamming_distance(seq1: str, seq2: str) -> int:
    """
    Calculate Hamming distance between two equal-length sequences
    
    Args:
        seq1: First sequence
        seq2: Second sequence
        
    Returns:
        int: Hamming distance
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must have equal length")
    
    return sum(c1 != c2 for c1, c2 in zip(seq1.upper(), seq2.upper()))


def format_sequence(seq: str, line_length: int = 80, add_numbers: bool = False) -> str:
    """
    Format sequence output
    
    Args:
        seq: DNA sequence
        line_length: Length per line
        add_numbers: Whether to add line numbers
        
    Returns:
        str: Formatted sequence
    """
    if not seq:
        return ""
    
    seq = clean_sequence(seq)
    lines = []
    
    for i in range(0, len(seq), line_length):
        line = seq[i:i+line_length]
        if add_numbers:
            line_num = i + 1
            line = f"{line_num:>8} {line}"
        lines.append(line)
    
    return '\n'.join(lines)


# Keep simple FASTA processing functionality as it may be needed for input file processing
def parse_fasta_string(fasta_string: str) -> List[Tuple[str, str]]:
    """
    Parse FASTA format string
    
    Args:
        fasta_string: FASTA format string
        
    Returns:
        List[Tuple[str, str]]: List of (sequence name, sequence) tuples
    """
    sequences = []
    lines = fasta_string.strip().split('\n')
    
    current_name = None
    current_seq = []
    
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            # Save previous sequence
            if current_name is not None:
                sequences.append((current_name, ''.join(current_seq)))
            
            # Start new sequence
            current_name = line[1:]  # Remove '>'
            current_seq = []
        elif line and current_name is not None:
            current_seq.append(line)
    
    # Save last sequence
    if current_name is not None:
        sequences.append((current_name, ''.join(current_seq)))
    
    return sequences


def to_fasta_string(sequences: List[Tuple[str, str]], line_length: int = 80) -> str:
    """
    Convert sequence list to FASTA format string
    
    Args:
        sequences: List of (sequence name, sequence) tuples
        line_length: Sequence length per line
        
    Returns:
        str: FASTA format string
    """
    fasta_lines = []
    
    for name, seq in sequences:
        fasta_lines.append(f">{name}")
        
        # Output sequence in lines
        for i in range(0, len(seq), line_length):
            fasta_lines.append(seq[i:i+line_length])
    
    return '\n'.join(fasta_lines) 