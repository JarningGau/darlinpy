#!/usr/bin/env python3
"""
DARLIN Main API Module

Provides simplified high-level interface for CARLIN sequence analysis
"""

from typing import List, Optional, Dict, Any, Union, Tuple
from dataclasses import dataclass, field
import time
import pandas as pd

from .config.amplicon_configs import AmpliconConfig
from .alignment.carlin_aligner import CARLINAligner
from .calling.allele_caller import AlleleCaller, AlleleCallResult
from .mutations.mutation import Mutation, annotate_mutations


@dataclass
class AnalysisResult:
    """
    CARLIN sequence analysis results
    
    Contains complete analysis results including allele calling, mutation annotation and statistics
    
    Attributes:
        called_alleles: List of called alleles
        mutations: List of mutations corresponding to each allele
        alignment_scores: Alignment scores for each sequence
        aligned_query: Aligned input sequences
        aligned_reference: Aligned reference sequences
        valid_sequences: List of valid sequences after filtering
        summary_stats: Summary statistics
        processing_time: Processing time in seconds
        config_used: Configuration name used
        method_used: Calling method used
    """
    called_alleles: List[AlleleCallResult]
    mutations: List[List[Mutation]]
    alignment_scores: List[float]
    summary_stats: Dict[str, Any] = field(default_factory=dict)  # Provide default to avoid dataclass parameter order errors
    aligned_query: List[str] = field(default_factory=list)
    aligned_reference: List[str] = field(default_factory=list)
    valid_sequences: List[str] = field(default_factory=list)
    processing_time: float = 0.0
    config_used: str = "OriginalCARLIN"
    method_used: str = "coarse_grain"

    def __post_init__(self):
        """Validate result data consistency"""
        if len(self.called_alleles) != len(self.mutations):
            raise ValueError("called_alleles and mutations must have the same length")
    
    @property
    def num_sequences(self) -> int:
        """Return total number of analyzed sequences"""
        return len(self.alignment_scores)
    
    @property
    def num_called_alleles(self) -> int:
        """Return number of successfully called alleles"""
        return len([a for a in self.called_alleles if a.is_callable()])
    
    @property
    def calling_success_rate(self) -> float:
        """Return allele calling success rate"""
        if not self.called_alleles:
            return 0.0
        return self.num_called_alleles / len(self.called_alleles)
    
    @property
    def total_mutations(self) -> int:
        """Return total number of detected mutations"""
        return sum(len(muts) for muts in self.mutations)
    
    @property
    def average_alignment_score(self) -> float:
        """Return average alignment score"""
        if not self.alignment_scores:
            return 0.0
        return sum(self.alignment_scores) / len(self.alignment_scores)
    
    def get_mutation_summary(self) -> Dict[str, int]:
        """Get mutation type statistics"""
        mutation_counts = {}
        for mut_list in self.mutations:
            for mut in mut_list:
                mut_type = mut.type.value
                mutation_counts[mut_type] = mutation_counts.get(mut_type, 0) + 1
        return mutation_counts

    def to_df(self) -> pd.DataFrame:
        """Convert analysis results to pandas DataFrame
        
        Returns:
            pd.DataFrame: DataFrame containing analysis results with columns:
                - query: Original input sequences
                - query_len: Length of each sequence
                - aligned_query: Aligned query sequences
                - aligned_ref: Aligned reference sequences
                - scores: Alignment scores
                - mutations: Comma-separated list of mutations in HGVS format
                - confidence: Comma-separated list of mutation confidence labels
        """
        mutations = []
        confidence = []
        for mut in self.mutations:
            _mut = ','.join([m.to_hgvs() for m in mut]) if len(mut) > 0 else []
            _conf = ','.join([m.confidence_label for m in mut]) if len(mut) > 0 else []
            mutations.append(_mut)
            confidence.append(_conf)
        
        results_df = pd.DataFrame({
            'query': self.valid_sequences,
            'query_len': [len(s) for s in self.valid_sequences], 
            'aligned_query': self.aligned_query, 
            'aligned_ref': self.aligned_reference,
            'scores': self.alignment_scores, 
            'mutations': mutations, 
            'confidence': confidence
        })
        return results_df
    
    def print_summary(self):
        """Print analysis results summary"""
        print(f"CARLIN Sequence Analysis Results Summary")
        print(f"=" * 40)
        print(f"Configuration: {self.config_used}")
        print(f"Method: {self.method_used}")
        print(f"Processing time: {self.processing_time:.2f} seconds")
        print(f"")
        print(f"Sequence statistics:")
        print(f"  Total sequences: {self.num_sequences}")
        print(f"  Successfully called alleles: {self.num_called_alleles}")
        print(f"  Calling success rate: {self.calling_success_rate:.1%}")
        print(f"  Average alignment score: {self.average_alignment_score:.2f}")
        print(f"")
        print(f"Mutation statistics:")
        print(f"  Total mutations: {self.total_mutations}")
        mut_summary = self.get_mutation_summary()
        for mut_type, count in mut_summary.items():
            print(f"  {mut_type} type: {count}")


def analyze_sequences(
    sequences: List[str],
    config: Union[str, AmpliconConfig] = 'Col1a1',  # Changed to default Col1a1
    method: str = 'coarse_grain',
    min_reads: int = 1,
    dominant_threshold: float = 0.5,
    annotate_mutations_flag: bool = True,
    merge_adjacent_mutations: bool = True,
    verbose: bool = False,
    min_sequence_length: int = 50
) -> AnalysisResult:
    """
    Analyze CARLIN sequences
    
    Args:
        sequences: List of sequences to analyze
        config: Configuration, can be locus name string or AmpliconConfig object, defaults to 'Col1a1', options: "Col1a1", "Rosa", "Tigre"
        method: Calling method ('exact' or 'coarse_grain')
        min_reads: Minimum number of reads
        dominant_threshold: Dominant allele threshold
        annotate_mutations_flag: Whether to annotate mutations
        merge_adjacent_mutations: Whether to merge adjacent mutations
        verbose: Whether to show detailed information
        min_sequence_length: Minimum sequence length threshold, sequences shorter than this will be filtered
        
    Returns:
        Analysis results
    """
    start_time = time.time()
    
    if verbose:
        print(f"Starting analysis of {len(sequences)} sequences...")
    
    # Process configuration parameters
    if isinstance(config, str):
        # If string, treat as locus
        if config == 'OriginalCARLIN':
            # Backward compatibility: if user explicitly specifies OriginalCARLIN, use the original way
            from .config.amplicon_configs import get_original_carlin_config
            amplicon_config = get_original_carlin_config()
        else:
            # Treat as locus
            from .config.amplicon_configs import load_carlin_config_by_locus
            amplicon_config = load_carlin_config_by_locus(config)
    else:
        # Directly use the provided AmpliconConfig
        amplicon_config = config
    
    if verbose:
        print(f"Using configuration: {config}")
    
    if method not in ['exact', 'coarse_grain']:
        raise ValueError(f"Unsupported calling method: {method}")
    
    if not sequences:
        raise ValueError("Sequence list cannot be empty")
    
    # Filter sequences (remove sequences that are too short or empty)
    valid_sequences = [seq for seq in sequences if seq and len(seq) >= min_sequence_length]
    if len(valid_sequences) < len(sequences):
        if verbose:
            print(f"Filtered {len(sequences) - len(valid_sequences)} short sequences (< {min_sequence_length}bp)")
    
    if not valid_sequences:
        raise ValueError("No valid sequences available for analysis")
    
    try:
        # 2. Sequence alignment
        if verbose:
            print("Performing sequence alignment...")
        
        aligner = CARLINAligner(amplicon_config=amplicon_config)
        alignment_results = aligner.align_sequences(valid_sequences)
        
        # Extract alignment scores and sequences
        alignment_scores = [result['alignment_score'] for result in alignment_results]
        aligned_sequences = [result.get('aligned_seq_obj') for result in alignment_results]
        
        if verbose:
            avg_score = sum(alignment_scores) / len(alignment_scores)
            print(f"Alignment completed, average score: {avg_score:.2f}")
        
        # 3. Allele calling
        if verbose:
            print("Performing allele calling...")
        
        # Initialize AlleleCaller
        caller = AlleleCaller(amplicon_config=amplicon_config, dominant_threshold=dominant_threshold)
        
        # Choose calling strategy based on method
        called_alleles = []
        for aligned_seq in aligned_sequences:
            if method == 'exact':
                allele_result = caller.call_alleles_exact([aligned_seq])
            else:  # coarse_grain
                allele_result = caller.call_alleles_coarse_grain([aligned_seq])
            
            called_alleles.append(allele_result)
        
        if verbose:
            callable_count = len([a for a in called_alleles if a.is_callable()])
            print(f"Allele calling completed, successfully called: {callable_count}/{len(called_alleles)}")
        
        # 4. Mutation annotation
        mutations_list = []
        if annotate_mutations_flag:
            if verbose:
                print("Performing mutation annotation...")
            
            # Get CARLIN cut site information
            cut_sites = _get_cut_sites(amplicon_config)
            
            for allele_result in called_alleles:
                if allele_result.is_callable() and allele_result.allele:
                    mutations = annotate_mutations(
                        allele_result.allele,
                        cas9_mode=True,
                        cut_sites=cut_sites,
                        merge_adjacent=merge_adjacent_mutations
                    )
                else:
                    mutations = []
                mutations_list.append(mutations)
            
            if verbose:
                total_mutations = sum(len(muts) for muts in mutations_list)
                print(f"Mutation annotation completed, detected {total_mutations} mutations")
        else:
            mutations_list = [[] for _ in called_alleles]
        
        # 5. Generate statistics
        summary_stats = _generate_summary_stats(
            valid_sequences, alignment_results, called_alleles, mutations_list
        )
        
        # 6. Build result object
        processing_time = time.time() - start_time
        
        result = AnalysisResult(
            called_alleles=called_alleles,
            mutations=mutations_list,
            alignment_scores=alignment_scores,
            summary_stats=summary_stats,
            aligned_query=[r["aligned_seq_obj"].get_seq() for r in alignment_results],
            aligned_reference=[r["aligned_seq_obj"].get_ref() for r in alignment_results],
            valid_sequences=valid_sequences,
            processing_time=processing_time,
            config_used=str(config),
            method_used=method
        )
        
        if verbose:
            print(f"Analysis completed, time taken: {processing_time:.2f} seconds")
            result.print_summary()
        
        return result
        
    except Exception as e:
        raise RuntimeError(f"Error occurred during analysis: {str(e)}") from e


def _get_cut_sites(amplicon_config: AmpliconConfig) -> List[int]:
    """Get CARLIN Cas9 cut sites
    
    Uses the cutsite intervals (last 7bp) of each segment in the configuration, takes the center of each interval and converts to 1-based coordinates.
    """
    # positions['cutsites'] are 0-based half-open intervals (start, end) relative to CARLIN
    # Here we take the center point of each interval and convert to 1-based coordinates to match mutation positioning
    cutsite_intervals = amplicon_config.positions.get('cutsites', [])
    centers_1_based: List[int] = []
    for start_0_based, end_0_based in cutsite_intervals:
        # Interval is [start, end), last base index is end-1
        last_index = end_0_based - 1
        center_0_based = start_0_based + (last_index - start_0_based) // 2
        centers_1_based.append(center_0_based + 1)
    return centers_1_based


def _generate_summary_stats(
    sequences: List[str],
    alignment_results: List,
    called_alleles: List[AlleleCallResult],
    mutations_list: List[List[Mutation]]
) -> Dict[str, Any]:
    """Generate summary statistics"""
    
    # Basic statistics
    stats = {
        'total_sequences': len(sequences),
        'avg_sequence_length': sum(len(seq) for seq in sequences) / len(sequences),
        'called_alleles_count': len([a for a in called_alleles if a.is_callable()]),
        'calling_success_rate': len([a for a in called_alleles if a.is_callable()]) / len(called_alleles),
    }
    
    # Alignment statistics
    scores = [r['alignment_score'] for r in alignment_results]
    stats.update({
        'avg_alignment_score': sum(scores) / len(scores),
        'min_alignment_score': min(scores),
        'max_alignment_score': max(scores),
    })
    
    # Mutation statistics
    mutation_counts = {}
    total_mutations = 0
    for mut_list in mutations_list:
        total_mutations += len(mut_list)
        for mut in mut_list:
            mut_type = mut.type.value
            mutation_counts[mut_type] = mutation_counts.get(mut_type, 0) + 1
    
    stats.update({
        'total_mutations': total_mutations,
        'avg_mutations_per_allele': total_mutations / len(called_alleles) if called_alleles else 0,
        'mutation_type_distribution': mutation_counts
    })
    
    return stats
