"""
CARLIN allele calling algorithm

This module implements the core CARLIN allele calling algorithm, including both exact calling and coarse-grained calling methods.
"""

from typing import List, Optional, Dict, Any, Tuple, Union
import numpy as np
from collections import Counter
from scipy import stats

from ..alignment.aligned_seq import AlignedSEQ, desemble_sequence
from ..config.amplicon_configs import AmpliconConfig
from .allele_data import AlleleCallResult, BulkAlleleCallResult


class AlleleCaller:
    """
    CARLIN allele caller
    
    Implements two calling strategies:
    1. call_alleles_exact: Exact calling based on complete sequences
    2. call_alleles_coarse_grain: Coarse-grained calling based on event structure
    """
    
    def __init__(self, 
             locus: str = "Col1a1",
             amplicon_config: Optional[AmpliconConfig] = None, 
             dominant_threshold: float = 0.5):
        """
        Initialize allele caller
        
        Args:
            locus: Locus name for selecting corresponding JSON configuration template, defaults to "Col1a1"
            amplicon_config: CARLIN amplicon configuration, takes priority if provided
            dominant_threshold: Minimum proportion threshold for dominant allele
        """
        # Load configuration based on locus
        if amplicon_config is None:
            from ..config.amplicon_configs import load_carlin_config_by_locus
            self.amplicon_config = load_carlin_config_by_locus(locus)
        else:
            self.amplicon_config = amplicon_config
        
        self.dominant_threshold = dominant_threshold
        
    def unique_by_frequency(self, 
                           values: List[Any], 
                           weights: Optional[List[float]] = None) -> Tuple[List[Any], List[int], List[int]]:
        """
        Get unique values sorted by frequency
        
        Args:
            values: Input value list
            weights: Weight for each value, defaults to 1
            
        Returns:
            (unique_values, backtrack, unique_indices)
            - unique_values: Unique values sorted by frequency in descending order
            - backtrack: Mapping from original index to unique value index
            - unique_indices: Unique value index corresponding to each original value
        """
        if weights is None:
            weights = [1.0] * len(values)
            
        if len(values) != len(weights):
            raise ValueError("values and weights must have the same length")
            
        # Create value to index mapping
        value_to_indices = {}
        for i, value in enumerate(values):
            # For AlignedSEQ objects, use sequence as key
            if isinstance(value, AlignedSEQ):
                key = value.get_seq()
            elif isinstance(value, str):
                key = value
            else:
                key = str(value)
                
            if key not in value_to_indices:
                value_to_indices[key] = []
            value_to_indices[key].append(i)
            
        # Calculate weight for each unique value
        unique_keys = []
        unique_weights = []
        unique_values = []
        
        for key, indices in value_to_indices.items():
            total_weight = sum(weights[i] for i in indices)
            unique_keys.append(key)
            unique_weights.append(total_weight)
            unique_values.append(values[indices[0]])  # Use first encountered value as representative
            
        # Sort by weight in descending order
        sorted_indices = sorted(range(len(unique_weights)), 
                              key=lambda i: unique_weights[i], 
                              reverse=True)
        
        sorted_unique_values = [unique_values[i] for i in sorted_indices]
        
        # Create mapping
        key_to_sorted_index = {}
        for new_idx, old_idx in enumerate(sorted_indices):
            key_to_sorted_index[unique_keys[old_idx]] = new_idx
            
        # Create backtrack and unique_indices
        backtrack = sorted_indices
        unique_indices = []
        
        for value in values:
            if isinstance(value, AlignedSEQ):
                key = value.get_seq()
            elif isinstance(value, str):
                key = value
            else:
                key = str(value)
            unique_indices.append(key_to_sorted_index[key])
            
        return sorted_unique_values, backtrack, unique_indices
        
    def call_alleles_exact(self, 
                          aligned_seqs: List[AlignedSEQ],
                          seq_weights: Optional[List[float]] = None,
                          dominant_only: bool = True) -> AlleleCallResult:
        """
        Exact allele calling based on complete sequences
        
        Args:
            aligned_seqs: List of aligned sequences
            seq_weights: Weight for each sequence
            dominant_only: Whether to return only dominant allele
            
        Returns:
            Allele calling result
        """
        if seq_weights is None:
            seq_weights = [1.0] * len(aligned_seqs)
            
        # Filter empty sequences
        valid_mask = [i for i, seq in enumerate(aligned_seqs) if seq is not None]
        if not valid_mask:
            return AlleleCallResult(
                allele=None,
                constituents=[],
                weight_contribution=[],
                calling_method='exact'
            )
            
        valid_seqs = [aligned_seqs[i] for i in valid_mask]
        valid_weights = [seq_weights[i] for i in valid_mask]
        
        # Extract sequence strings
        seq_strings = [seq.get_seq() for seq in valid_seqs]
        
        # Get unique sequences by frequency
        unique_seqs, _, which_seq = self.unique_by_frequency(seq_strings, valid_weights)
        
        # Calculate weight for each unique sequence
        seq_weight_sums = []
        for i in range(len(unique_seqs)):
            weight_sum = sum(valid_weights[j] for j, seq_idx in enumerate(which_seq) if seq_idx == i)
            seq_weight_sums.append(weight_sum)
            
        total_weight = sum(seq_weight_sums)
        
        # Check dominance
        dominant_fraction = seq_weight_sums[0] / total_weight if total_weight > 0 else 0
        
        if dominant_only and dominant_fraction < self.dominant_threshold:
            return AlleleCallResult(
                allele=None,
                constituents=[],
                weight_contribution=[],
                calling_method='exact',
                confidence=dominant_fraction
            )
            
        # Find representative of highest weight sequence
        target_seq_indices = [i for i, seq_idx in enumerate(which_seq) if seq_idx == 0]
        ref_idx = max(target_seq_indices, key=lambda i: valid_weights[i])
        
        allele = valid_seqs[ref_idx]
        constituents = [valid_mask[i] for i in target_seq_indices]
        weight_contribution = [valid_weights[i] for i in target_seq_indices]
        
        return AlleleCallResult(
            allele=allele,
            constituents=constituents,
            weight_contribution=weight_contribution,
            confidence=dominant_fraction,
            calling_method='exact',
            dominant_fraction=dominant_fraction
        )
        
    def call_alleles_coarse_grain(self, 
                                 aligned_seqs: List[AlignedSEQ],
                                 seq_weights: Optional[List[float]] = None,
                                 dominant_only: bool = True) -> AlleleCallResult:
        """
        Coarse-grained allele calling based on event structure
        
        Args:
            aligned_seqs: List of aligned sequences
            seq_weights: Weight for each sequence
            dominant_only: Whether to return only dominant allele
            
        Returns:
            Allele calling result
        """
        if seq_weights is None:
            seq_weights = [1.0] * len(aligned_seqs)
            
        # Filter empty sequences
        valid_mask = [i for i, seq in enumerate(aligned_seqs) if seq is not None]
        if not valid_mask:
            return AlleleCallResult(
                allele=None,
                constituents=[],
                weight_contribution=[],
                calling_method='coarse_grain'
            )
            
        valid_seqs = [aligned_seqs[i] for i in valid_mask]
        valid_weights = [seq_weights[i] for i in valid_mask]
        
        # Extract event structures
        event_structures = [seq.get_event_structure() for seq in valid_seqs]
        event_strings = [''.join(events) for events in event_structures]
        
        # Get unique event structures by frequency
        unique_events, _, which_event = self.unique_by_frequency(event_strings, valid_weights)
        
        # Calculate weight for each unique event
        event_weight_sums = []
        for i in range(len(unique_events)):
            weight_sum = sum(valid_weights[j] for j, event_idx in enumerate(which_event) if event_idx == i)
            event_weight_sums.append(weight_sum)
            
        total_weight = sum(event_weight_sums)
        
        # Check dominance
        dominant_fraction = event_weight_sums[0] / total_weight if total_weight > 0 else 0
        
        if dominant_only and dominant_fraction < self.dominant_threshold:
            return AlleleCallResult(
                allele=None,
                constituents=[],
                weight_contribution=[],
                calling_method='coarse_grain',
                confidence=dominant_fraction
            )
            
        # Collect sequences with dominant event structure
        target_event_indices = [i for i, event_idx in enumerate(which_event) if event_idx == 0]
        event_seqs = [valid_seqs[i] for i in target_event_indices]
        event_weights = [valid_weights[i] for i in target_event_indices]
        
        # Build consensus sequence
        consensus_allele = self._build_consensus_sequence(event_seqs, event_weights)
        
        constituents = [valid_mask[i] for i in target_event_indices]
        weight_contribution = event_weights
        
        return AlleleCallResult(
            allele=consensus_allele,
            constituents=constituents,
            weight_contribution=weight_contribution,
            confidence=dominant_fraction,
            calling_method='coarse_grain',
            dominant_fraction=dominant_fraction
        )
        
    def _build_consensus_sequence(self, 
                                 seqs: List[AlignedSEQ], 
                                 weights: List[float]) -> AlignedSEQ:
        """
        Build consensus sequence from multiple sequences
        
        Args:
            seqs: Sequence list
            weights: Weight list
            
        Returns:
            Consensus AlignedSEQ object
        """
        if not seqs:
            raise ValueError("Sequence list cannot be empty")
            
        # Extract sequence strings
        seq_strings = [seq.get_seq() for seq in seqs]
        ref_strings = [seq.get_ref() for seq in seqs]
        
        # Check if all sequences have the same length
        seq_lengths = [len(s) for s in seq_strings]
        if len(set(seq_lengths)) == 1:
            # Same length, use mode statistics
            consensus_seq = self._compute_mode_sequence(seq_strings, weights)
            # Use first sequence's reference as reference
            consensus_ref = ref_strings[0]
        else:
            # Different lengths, use most common length sequences
            length_weights = {}
            for i, length in enumerate(seq_lengths):
                if length not in length_weights:
                    length_weights[length] = 0
                length_weights[length] += weights[i]
                
            mode_length = max(length_weights.keys(), key=lambda x: length_weights[x])
            
            # Filter to mode length sequences
            mode_indices = [i for i, length in enumerate(seq_lengths) if length == mode_length]
            mode_seqs = [seq_strings[i] for i in mode_indices]
            mode_weights = [weights[i] for i in mode_indices]
            
            consensus_seq = self._compute_mode_sequence(mode_seqs, mode_weights)
            consensus_ref = ref_strings[mode_indices[0]]
            
        # Rebuild AlignedSEQ object
        try:
            # Use desemble_sequence to rebuild
            from ..alignment.aligned_seq import calculate_motif_boundaries
            motif_boundaries = calculate_motif_boundaries(consensus_ref, self.amplicon_config)
            return desemble_sequence(consensus_seq, consensus_ref, motif_boundaries)
        except Exception as e:
            # If failed, create simple single motif AlignedSEQ
            return AlignedSEQ([consensus_seq], [consensus_ref])
            
    def _compute_mode_sequence(self, 
                              sequences: List[str], 
                              weights: List[float]) -> str:
        """
        Compute weighted mode sequence
        
        Args:
            sequences: List of sequence strings
            weights: Weight list
            
        Returns:
            Consensus sequence string
        """
        if not sequences:
            return ""
            
        # Ensure all sequences have the same length
        seq_length = len(sequences[0])
        if not all(len(seq) == seq_length for seq in sequences):
            raise ValueError("All sequences must have the same length")
            
        consensus = []
        
        for pos in range(seq_length):
            # Collect characters and weights at this position
            char_weights = {}
            for seq, weight in zip(sequences, weights):
                char = seq[pos]
                if char not in char_weights:
                    char_weights[char] = 0
                char_weights[char] += weight
                
            # Choose character with highest weight
            mode_char = max(char_weights.keys(), key=lambda x: char_weights[x])
            consensus.append(mode_char)
            
        return ''.join(consensus)
        
    def call_bulk_alleles(self, 
                         sequence_groups: List[List[AlignedSEQ]],
                         group_weights: Optional[List[List[float]]] = None,
                         method: str = 'coarse_grain',
                         dominant_only: bool = True) -> BulkAlleleCallResult:
        """
        Bulk allele calling
        
        Args:
            sequence_groups: List of sequence groups, each group represents a UMI or cell
            group_weights: Weights for sequences in each group
            method: Calling method ('exact' or 'coarse_grain')
            dominant_only: Whether to return only dominant allele
            
        Returns:
            Bulk allele calling result
        """
        if group_weights is None:
            group_weights = [[1.0] * len(group) for group in sequence_groups]
            
        if len(sequence_groups) != len(group_weights):
            raise ValueError("sequence_groups and group_weights must have the same length")
            
        # Call alleles for each group
        individual_results = []
        calling_method = getattr(self, f'call_alleles_{method}')
        
        for group_seqs, group_w in zip(sequence_groups, group_weights):
            result = calling_method(group_seqs, group_w, dominant_only)
            individual_results.append(result)
            
        # Aggregate all successfully called alleles
        successful_alleles = [r.allele for r in individual_results if r.is_callable()]
        successful_weights = [r.total_weight for r in individual_results if r.is_callable()]
        
        if not successful_alleles:
            return BulkAlleleCallResult(
                individual_results=individual_results,
                summary_alleles=[],
                allele_frequencies=[],
                total_callable_sequences=0,
                calling_parameters={
                    'method': method,
                    'dominant_only': dominant_only,
                    'dominant_threshold': self.dominant_threshold
                }
            )
            
        # Aggregate alleles by frequency
        allele_seqs = [allele.get_seq() for allele in successful_alleles]
        unique_alleles, _, _ = self.unique_by_frequency(allele_seqs, successful_weights)
        
        # Calculate frequencies
        allele_counts = Counter()
        total_weight = sum(successful_weights)
        
        for allele_seq, weight in zip(allele_seqs, successful_weights):
            allele_counts[allele_seq] += weight
            
        # Sort and calculate frequencies
        sorted_alleles = sorted(allele_counts.items(), key=lambda x: x[1], reverse=True)
        summary_alleles = []
        allele_frequencies = []
        
        for allele_seq, count in sorted_alleles:
            # Find corresponding AlignedSEQ object
            for allele in successful_alleles:
                if allele.get_seq() == allele_seq:
                    summary_alleles.append(allele)
                    break
            allele_frequencies.append(count / total_weight)
            
        return BulkAlleleCallResult(
            individual_results=individual_results,
            summary_alleles=summary_alleles,
            allele_frequencies=allele_frequencies,
            total_callable_sequences=len(successful_alleles),
            calling_parameters={
                'method': method,
                'dominant_only': dominant_only,
                'dominant_threshold': self.dominant_threshold,
                'total_groups': len(sequence_groups),
                'successful_groups': len(successful_alleles)
            }
        ) 