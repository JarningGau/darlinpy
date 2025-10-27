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


@dataclass 
class BulkAlleleCallResult:
    """
    Result of bulk allele calling
    
    Attributes:
        individual_results: List of calling results for each UMI/group
        summary_alleles: List of summary alleles
        allele_frequencies: Frequency of each allele
        total_callable_sequences: Total number of callable sequences
        calling_parameters: Calling parameters
        
    @deprecated This class is not used in the main API (analyze_sequences, AmpliconConfig, build_carlin_config).
    It may be removed in future versions.
    """
    individual_results: List[AlleleCallResult]
    summary_alleles: List[AlignedSEQ]
    allele_frequencies: List[float]
    total_callable_sequences: int
    calling_parameters: Dict[str, Any]
    
    def __post_init__(self):
        """Validate data consistency"""
        if len(self.summary_alleles) != len(self.allele_frequencies):
            raise ValueError("summary_alleles and allele_frequencies must have the same length")
            
    @property
    def num_called_alleles(self) -> int:
        """Return number of successfully called alleles"""
        return len([r for r in self.individual_results if r.is_callable()])
        
    @property
    def calling_success_rate(self) -> float:
        """Return calling success rate"""
        if not self.individual_results:
            return 0.0
        return self.num_called_alleles / len(self.individual_results)
        
    @property
    def dominant_allele(self) -> Optional[AlignedSEQ]:
        """Return the most dominant allele"""
        if not self.summary_alleles:
            return None
        return self.summary_alleles[0]
        
    @property 
    def dominant_frequency(self) -> float:
        """Return frequency of the most dominant allele"""
        if not self.allele_frequencies:
            return 0.0
        return self.allele_frequencies[0]
        
    def get_allele_by_rank(self, rank: int) -> Optional[AlignedSEQ]:
        """Get allele by frequency rank"""
        if rank < 0 or rank >= len(self.summary_alleles):
            return None
        return self.summary_alleles[rank]
        
    def get_frequency_by_rank(self, rank: int) -> float:
        """Get allele frequency by rank"""
        if rank < 0 or rank >= len(self.allele_frequencies):
            return 0.0
        return self.allele_frequencies[rank]


class AlleleCallStatistics:
    """
    Allele calling statistics
    
    @deprecated This class is not used in the main API (analyze_sequences, AmpliconConfig, build_carlin_config).
    It may be removed in future versions.
    """
    
    def __init__(self, results: Union[AlleleCallResult, BulkAlleleCallResult, List[AlleleCallResult]]):
        """
        Initialize statistics calculator
        
        Args:
            results: Allele calling results
        """
        if isinstance(results, AlleleCallResult):
            self.individual_results = [results]
        elif isinstance(results, BulkAlleleCallResult):
            self.individual_results = results.individual_results
        elif isinstance(results, list):
            self.individual_results = results
        else:
            raise ValueError("Unsupported result type")
            
    def compute_calling_metrics(self) -> Dict[str, float]:
        """
        Compute calling metrics
        
        Returns:
            Dictionary containing various statistical metrics
        """
        callable_results = [r for r in self.individual_results if r.is_callable()]
        
        if not self.individual_results:
            return {
                'success_rate': 0.0,
                'mean_confidence': 0.0,
                'mean_weight': 0.0,
                'mean_supporting_sequences': 0.0
            }
            
        metrics = {
            'success_rate': len(callable_results) / len(self.individual_results),
            'total_attempts': len(self.individual_results),
            'successful_calls': len(callable_results)
        }
        
        if callable_results:
            metrics.update({
                'mean_confidence': np.mean([r.confidence for r in callable_results]),
                'median_confidence': np.median([r.confidence for r in callable_results]),
                'mean_weight': np.mean([r.total_weight for r in callable_results]),
                'median_weight': np.median([r.total_weight for r in callable_results]),
                'mean_supporting_sequences': np.mean([r.num_supporting_sequences for r in callable_results]),
                'median_supporting_sequences': np.median([r.num_supporting_sequences for r in callable_results])
            })
        else:
            metrics.update({
                'mean_confidence': 0.0,
                'median_confidence': 0.0,
                'mean_weight': 0.0,
                'median_weight': 0.0,
                'mean_supporting_sequences': 0.0,
                'median_supporting_sequences': 0.0
            })
            
        return metrics
        
    def compute_event_distribution(self) -> Dict[str, int]:
        """
        Compute event type distribution
        
        Returns:
            Mapping from event type to count
        """
        event_counts = {}
        for result in self.individual_results:
            if result.is_callable():
                event_structure = result.event_structure
                if event_structure:
                    event_counts[event_structure] = event_counts.get(event_structure, 0) + 1
                    
        return event_counts
        
    def compute_method_distribution(self) -> Dict[str, int]:
        """
        Compute calling method distribution
        
        Returns:
            Mapping from calling method to count
        """
        method_counts = {}
        for result in self.individual_results:
            if result.is_callable():
                method = result.calling_method
                method_counts[method] = method_counts.get(method, 0) + 1
                
        return method_counts 