"""
CARLIN allele calling data structures

This module defines the core data structures used in the allele calling process.
"""

from typing import List, Optional, Dict, Any, Union
from dataclasses import dataclass
import numpy as np
from ..alignment.aligned_seq import AlignedSEQ


@dataclass
class AlleleCallResult:
    """
    Result of a single allele call
    
    Attributes:
        allele: Called allele sequence (AlignedSEQ object)
        constituents: List of sequence indices supporting this allele
        weight_contribution: Weight of each supporting sequence
        confidence: Calling confidence (0-1)
        calling_method: Calling method used ('exact' or 'coarse_grain')
        dominant_fraction: Proportion of dominant allele
    """
    allele: Optional[AlignedSEQ]
    constituents: List[int]
    weight_contribution: List[float]
    confidence: float = 0.0
    calling_method: str = 'coarse_grain'
    dominant_fraction: float = 0.0
    
    def __post_init__(self):
        """Validate data consistency"""
        if self.allele is not None:
            if len(self.constituents) != len(self.weight_contribution):
                raise ValueError("constituents and weight_contribution must have the same length")
            
    @property
    def total_weight(self) -> float:
        """Return total weight supporting this allele"""
        return sum(self.weight_contribution)
        
    @property
    def num_supporting_sequences(self) -> int:
        """Return number of sequences supporting this allele"""
        return len(self.constituents)
        
    @property
    def event_structure(self) -> Optional[str]:
        """Return event structure of the allele"""
        if self.allele is None:
            return None
        return ''.join(self.allele.get_event_structure())
        
    @property
    def sequence(self) -> Optional[str]:
        """Return allele sequence (with gaps removed)"""
        if self.allele is None:
            return None
        return self.allele.get_seq().replace('-', '')
        
    def is_callable(self) -> bool:
        """Check if allele was successfully called"""
        return self.allele is not None
