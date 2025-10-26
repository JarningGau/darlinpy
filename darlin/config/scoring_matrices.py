#!/usr/bin/env python3
"""
Scoring matrices module

Implements NUC44 and other nucleic acid sequence alignment scoring matrices
"""

import numpy as np
from typing import Dict, Tuple


def create_nuc44_matrix() -> np.ndarray:
    """
    Create NUC44 nucleic acid substitution scoring matrix
    
    This is the standard nucleic acid scoring matrix used by NCBI BLAST
    Index order: 0=gap, 1=A, 2=C, 3=G, 4=T
    
    Returns:
        np.ndarray: 5x5 scoring matrix flattened to 25-element 1D array
    """
    # NUC44 matrix (5x5)
    # Row/column order: gap(0), A(1), C(2), G(3), T(4)
    nuc44_2d = np.array([
        #    gap   A    C    G    T
        [   0,   0,   0,   0,   0],  # gap
        [   0,   5,  -4,  -4,  -4],  # A
        [   0,  -4,   5,  -4,  -4],  # C  
        [   0,  -4,  -4,   5,  -4],  # G
        [   0,  -4,  -4,  -4,   5],  # T
    ], dtype=np.float64)
    
    # Flatten to 1D array, consistent with cas9_align expected format
    return nuc44_2d.flatten()


def create_simple_scoring_matrix(match_score: float = 5.0, mismatch_score: float = -4.0) -> np.ndarray:
    """
    Create simple match/mismatch scoring matrix
    
    Args:
        match_score: Match score
        mismatch_score: Mismatch score
        
    Returns:
        np.ndarray: Flattened 5x5 scoring matrix
    """
    matrix = np.zeros(25, dtype=np.float64)
    
    # A=1, C=2, G=3, T=4, gap=0
    for i in range(1, 5):
        for j in range(1, 5):
            idx = i * 5 + j
            if i == j:
                matrix[idx] = match_score    # match
            else:
                matrix[idx] = mismatch_score # mismatch
    
    return matrix


def create_transition_transversion_matrix(
    match_score: float = 5.0,
    transition_score: float = -1.0,  # A<->G, C<->T
    transversion_score: float = -4.0  # other substitutions
) -> np.ndarray:
    """
    Create scoring matrix considering transitions/transversions
    
    Transition: A<->G, C<->T (usually more common)
    Transversion: A<->C, A<->T, G<->C, G<->T
    
    Args:
        match_score: Perfect match score
        transition_score: Transition score (lighter penalty)
        transversion_score: Transversion score (heavier penalty)
        
    Returns:
        np.ndarray: Flattened 5x5 scoring matrix
    """
    matrix = np.zeros(25, dtype=np.float64)
    
    # Define transition pairs
    transitions = {(1, 3), (3, 1), (2, 4), (4, 2)}  # A<->G, C<->T
    
    for i in range(1, 5):
        for j in range(1, 5):
            idx = i * 5 + j
            if i == j:
                matrix[idx] = match_score
            elif (i, j) in transitions:
                matrix[idx] = transition_score
            else:
                matrix[idx] = transversion_score
    
    return matrix


class ScoringConfig:
    """Scoring configuration class"""
    
    def __init__(self, matrix_type: str = "nuc44", **kwargs):
        """
        Initialize scoring configuration
        
        Args:
            matrix_type: Matrix type ("nuc44", "simple", "transition_transversion")
            **kwargs: Parameters passed to matrix creation functions
        """
        self.matrix_type = matrix_type
        self.kwargs = kwargs
        
        if matrix_type == "nuc44":
            self.substitution_matrix = create_nuc44_matrix()
        elif matrix_type == "simple":
            self.substitution_matrix = create_simple_scoring_matrix(**kwargs)
        elif matrix_type == "transition_transversion":
            self.substitution_matrix = create_transition_transversion_matrix(**kwargs)
        else:
            raise ValueError(f"Unknown matrix type: {matrix_type}")
    
    def get_score(self, base1: int, base2: int) -> float:
        """
        Get score between two bases
        
        Args:
            base1, base2: Base encoding (0=gap, 1=A, 2=C, 3=G, 4=T)
            
        Returns:
            float: Substitution score
        """
        return self.substitution_matrix[base1 * 5 + base2]
    
    def summary(self) -> str:
        """Return scoring configuration summary"""
        lines = [
            f"=== Scoring Matrix Configuration ({self.matrix_type}) ===",
            "",
            "Scoring matrix (row: query sequence, column: reference sequence):",
            "        gap     A     C     G     T"
        ]
        
        base_names = ["gap", "A  ", "C  ", "G  ", "T  "]
        for i in range(5):
            row_scores = []
            for j in range(5):
                score = self.substitution_matrix[i * 5 + j]
                row_scores.append(f"{score:5.1f}")
            lines.append(f"{base_names[i]:>3}  " + "  ".join(row_scores))
        
        return "\n".join(lines)


# Predefined scoring configurations
DEFAULT_NUC44 = ScoringConfig("nuc44")
DEFAULT_SIMPLE = ScoringConfig("simple", match_score=5.0, mismatch_score=-4.0)


def get_default_scoring_config() -> ScoringConfig:
    """Get default scoring configuration (NUC44)"""
    return DEFAULT_NUC44 