#!/usr/bin/env python3
"""
CARLIN-specific sequence aligner

Integrates CARLIN configuration and cas9_align algorithm to provide advanced alignment interface
"""

import numpy as np
from typing import Tuple, Optional, Dict, List
from pathlib import Path

from .cas9_align import cas9_align, nt2int, int2nt, print_cas9_alignment
from .aligned_seq import AlignedSEQ, AlignedSEQMotif, SequenceSanitizer, desemble_sequence, calculate_motif_boundaries
from ..config import AmpliconConfig, get_original_carlin_config, ScoringConfig, get_default_scoring_config


class CARLINAligner:
    """
    CARLIN-specific sequence aligner
    
    Provides advanced alignment functionality for CARLIN amplicons with automatic configuration and parameter handling
    """
    
    def __init__(self, 
                 locus: str = "Col1a1",
                 amplicon_config: Optional[AmpliconConfig] = None,
                 scoring_config: Optional[ScoringConfig] = None):
        """
        Initialize aligner
        
        Args:
            locus: Locus name for selecting corresponding JSON configuration template, defaults to "Col1a1"
            amplicon_config: CARLIN amplicon configuration, takes priority if provided
            scoring_config: Scoring configuration, defaults to NUC44
        """
        # Load configuration based on locus
        if amplicon_config is None:
            from ..config.amplicon_configs import load_carlin_config_by_locus
            self.amplicon_config = load_carlin_config_by_locus(locus)
        else:
            self.amplicon_config = amplicon_config
        
        self.scoring_config = scoring_config or get_default_scoring_config()
        
        # Get precomputed parameters
        self.reference_sequence = self.amplicon_config.get_reference_sequence()
        self.open_penalty_array, self.close_penalty_array = self.amplicon_config.get_penalty_arrays()
        self.substitution_matrix = self.scoring_config.substitution_matrix
        
        # Encode reference sequence
        self.reference_encoded = nt2int(self.reference_sequence)
        
        print(f"✅ CARLIN aligner initialized successfully")
        print(f"   - Reference sequence length: {len(self.reference_sequence)} bp")
        print(f"   - Scoring matrix: {self.scoring_config.matrix_type}")
        print(f"   - Gap penalty range: {self.open_penalty_array.min():.1f}-{self.open_penalty_array.max():.1f}")
    
    def align_sequence(self, query_sequence: str, 
                      verbose: bool = False,
                      sanitize: bool = True) -> Dict:
        """
        Align a single sequence to CARLIN reference
        
        Args:
            query_sequence: Query sequence
            verbose: Whether to show detailed information
            sanitize: Whether to perform sequence normalization
            
        Returns:
            Dict: Alignment result dictionary
        """
        # Input validation
        if not query_sequence or not all(c in 'ACGTN-' for c in query_sequence.upper()):
            raise ValueError("Sequence can only contain ACGTN- characters")
        
        # Encode query sequence
        query_encoded = nt2int(query_sequence.upper().replace('N', 'A'))  # Replace N with A
        
        # Execute alignment
        score, aligned_query, aligned_ref = cas9_align(
            query_encoded, 
            self.reference_encoded,
            self.open_penalty_array,
            self.close_penalty_array,
            self.substitution_matrix
        )
        
        # Decode alignment results
        aligned_query_str = int2nt(aligned_query)
        aligned_ref_str = int2nt(aligned_ref)
        
        # Sequence normalization (optional)
        sanitized_aligned_seq = None
        if sanitize:
            try:
                # Decompose alignment results into motifs
                motif_boundaries = calculate_motif_boundaries(aligned_ref_str, self.amplicon_config)
                aligned_seq_obj = desemble_sequence(aligned_query_str, aligned_ref_str, motif_boundaries)
                
                # Perform normalization
                sanitized_seq = self._perform_sanitization(aligned_seq_obj)
                
                # Update alignment results with normalized sequence
                sanitized_aligned_query_str = sanitized_seq.get_seq()
                sanitized_aligned_ref_str = sanitized_seq.get_ref()
                sanitized_aligned_seq = sanitized_seq
                
                # Use normalized sequence for subsequent analysis
                aligned_query_str = sanitized_aligned_query_str
                aligned_ref_str = sanitized_aligned_ref_str
                
            except Exception as e:
                if verbose:
                    print(f"⚠️  Sequence normalization failed: {e}")
                    print("   Continuing with original alignment results")
        
        # Re-encode for statistical calculation
        if sanitize and sanitized_aligned_seq:
            # Use normalized sequence for statistics
            aligned_query_for_stats = nt2int(aligned_query_str)
            aligned_ref_for_stats = nt2int(aligned_ref_str) 
        else:
            aligned_query_for_stats = aligned_query
            aligned_ref_for_stats = aligned_ref
        
        # Calculate alignment statistics
        stats = self._calculate_alignment_stats(aligned_query_for_stats, aligned_ref_for_stats)
        
        # Build result
        result = {
            'query_sequence': query_sequence,
            'reference_sequence': self.reference_sequence,
            'aligned_query': aligned_query_str,
            'aligned_reference': aligned_ref_str,
            'alignment_score': score,
            'statistics': stats,
            'motif_analysis': self._analyze_motifs(aligned_query_for_stats, aligned_ref_for_stats),
            'sanitized': sanitize and sanitized_aligned_seq is not None,
            'aligned_seq_obj': sanitized_aligned_seq if sanitized_aligned_seq else None
        }
        
        if verbose:
            self._print_alignment_result(result)
        
        return result
    
    def align_sequences(self, sequences: List[str], 
                       verbose: bool = False) -> List[Dict]:
        """
        Batch align multiple sequences
        
        Args:
            sequences: List of sequences
            verbose: Whether to show detailed information
            
        Returns:
            List[Dict]: List of alignment results
        """
        results = []
        
        for i, seq in enumerate(sequences):
            if verbose:
                print(f"\n=== Aligning sequence {i+1}/{len(sequences)} ===")
            
            try:
                result = self.align_sequence(seq, verbose=verbose)
                results.append(result)
            except Exception as e:
                print(f"⚠️  Sequence {i+1} alignment failed: {e}")
                results.append({
                    'query_sequence': seq,
                    'error': str(e),
                    'alignment_score': float('-inf')
                })
        
        return results
    
    def _calculate_alignment_stats(self, aligned_query: np.ndarray, 
                                  aligned_ref: np.ndarray) -> Dict:
        """Calculate alignment statistics"""
        matches = 0
        mismatches = 0
        query_gaps = 0
        ref_gaps = 0
        
        for i in range(len(aligned_query)):
            q_base = aligned_query[i]
            r_base = aligned_ref[i]
            
            if q_base == 0:  # Query sequence gap
                query_gaps += 1
            elif r_base == 0:  # Reference sequence gap
                ref_gaps += 1
            elif q_base == r_base:  # Match
                matches += 1
            else:  # Mismatch
                mismatches += 1
        
        total_aligned = len(aligned_query)
        identity = matches / total_aligned if total_aligned > 0 else 0.0
        
        return {
            'aligned_length': total_aligned,
            'matches': matches,
            'mismatches': mismatches,
            'query_gaps': query_gaps,
            'reference_gaps': ref_gaps,
            'identity': identity,
            'similarity': (matches - mismatches) / total_aligned if total_aligned > 0 else 0.0
        }
    
    def _analyze_motifs(self, aligned_query: np.ndarray, 
                       aligned_ref: np.ndarray) -> Dict:
        """Analyze alignment status for each motif"""
        motif_stats = {
            'prefix': {'matches': 0, 'mismatches': 0, 'gaps': 0},
            'segments': [{} for _ in range(10)],
            'pams': [{} for _ in range(9)],
            'postfix': {'matches': 0, 'mismatches': 0, 'gaps': 0}
        }
        
        # Analyze motif information for each position
        ref_pos = 0  # Actual position in reference sequence
        
        for i in range(len(aligned_query)):
            q_base = aligned_query[i]
            r_base = aligned_ref[i]
            
            # Only analyze motif when reference sequence is not a gap
            if r_base != 0:
                motif_info = self.amplicon_config.get_motif_info(ref_pos)
                motif_type = motif_info['type']
                motif_id = motif_info['motif_id']
                
                # Count matches/mismatches/gaps
                if q_base == 0:
                    event = 'gaps'
                elif q_base == r_base:
                    event = 'matches'
                else:
                    event = 'mismatches'
                
                # Update statistics for corresponding motif
                if motif_type in ['prefix', 'postfix']:
                    motif_stats[motif_type][event] += 1
                elif motif_type in ['consite', 'cutsite']:
                    seg_id = motif_id
                    if seg_id not in motif_stats['segments'][seg_id]:
                        motif_stats['segments'][seg_id] = {'matches': 0, 'mismatches': 0, 'gaps': 0,
                                                         'consite': {'matches': 0, 'mismatches': 0, 'gaps': 0},
                                                         'cutsite': {'matches': 0, 'mismatches': 0, 'gaps': 0}}
                    
                    motif_stats['segments'][seg_id][event] += 1
                    motif_stats['segments'][seg_id][motif_type][event] += 1
                elif motif_type == 'pam':
                    pam_id = motif_id
                    if pam_id not in motif_stats['pams'][pam_id]:
                        motif_stats['pams'][pam_id] = {'matches': 0, 'mismatches': 0, 'gaps': 0}
                    motif_stats['pams'][pam_id][event] += 1
                
                ref_pos += 1
        
        return motif_stats
    
    def _print_alignment_result(self, result: Dict):
        """Print alignment results"""
        print(f"Query sequence length: {len(result['query_sequence'])} bp")
        print(f"Alignment score: {result['alignment_score']:.2f}")
        
        stats = result['statistics']
        print(f"Alignment statistics:")
        print(f"  - Aligned length: {stats['aligned_length']}")
        print(f"  - Matches: {stats['matches']} ({stats['identity']*100:.1f}%)")
        print(f"  - Mismatches: {stats['mismatches']}")
        print(f"  - Query gaps: {stats['query_gaps']}")
        print(f"  - Reference gaps: {stats['reference_gaps']}")
        
        # Show alignment results
        print("\nAlignment results:")
        aligned_query = result['aligned_query']
        aligned_ref = result['aligned_reference']
        
        # Display in segments, 60 characters per line
        line_length = 60
        for start in range(0, len(aligned_query), line_length):
            end = min(start + line_length, len(aligned_query))
            
            query_line = aligned_query[start:end]
            ref_line = aligned_ref[start:end]
            
            # Build match indicator line
            match_line = ""
            for i in range(len(query_line)):
                if query_line[i] == ref_line[i] and query_line[i] != '-':
                    match_line += "|"
                elif query_line[i] == '-' or ref_line[i] == '-':
                    match_line += " "
                else:
                    match_line += "."
            
            print(f"Query {start+1:>3}: {query_line}")
            print(f"           {match_line}")
            print(f"Ref   {start+1:>3}: {ref_line}")
            print()
    
    def get_config_summary(self) -> str:
        """Get configuration summary"""
        lines = [
            "=== CARLIN Aligner Configuration ===",
            f"Reference sequence: {self.reference_sequence[:50]}...",
            f"Sequence length: {len(self.reference_sequence)} bp",
            f"Scoring matrix: {self.scoring_config.matrix_type}",
            f"Match score: {self.scoring_config.get_score(1, 1)}",
            f"Mismatch score: {self.scoring_config.get_score(1, 2)}",
            f"Gap penalty range: {self.open_penalty_array.min():.1f} - {self.open_penalty_array.max():.1f}",
            f"Number of segments: {len(self.amplicon_config.sequence.segments)}",
            f"PAM sequence: {self.amplicon_config.sequence.pam}"
        ]
        return "\n".join(lines)
    
    def _perform_sanitization(self, aligned_seq: AlignedSEQ) -> AlignedSEQ:
        """
        Perform sequence normalization
        
        Args:
            aligned_seq: Decomposed alignment sequence object
            
        Returns:
            Normalized AlignedSEQ object
        """
        # Step 1: Prefix/Postfix normalization
        sanitized_seq = SequenceSanitizer.sanitize_prefix_postfix(aligned_seq)
        
        # Step 2: Conserved region normalization
        # Need to determine cutsite motif indices
        cutsite_indices = self._get_cutsite_motif_indices()
        sanitized_seq = SequenceSanitizer.sanitize_conserved_regions(sanitized_seq, cutsite_indices)
        
        return sanitized_seq
    
    def _get_cutsite_motif_indices(self) -> List[int]:
        """
        Get cutsite motif indices in motif list
        
        According to CARLIN configuration structure, motif order is:
        prefix, consite1, cutsite1, pam1, consite2, cutsite2, pam2, ..., consite10, cutsite10, postfix
        
        Returns:
            List of cutsite motif indices
        """
        cutsite_indices = []
        motif_index = 0
        
        # Prefix (index 0)
        motif_index += 1
        
        # 10 segments, each containing consite and cutsite
        for i in range(10):
            # Consite
            motif_index += 1
            
            # Cutsite (this is what we want to mark)
            cutsite_indices.append(motif_index)
            motif_index += 1
            
            # PAM (first 9 segments have PAM after them)
            if i < 9:
                motif_index += 1
        
        # Postfix is already at the end
        
        return cutsite_indices


def create_default_aligner() -> CARLINAligner:
    """Create default CARLIN aligner"""
    return CARLINAligner()


def align_to_carlin(sequence: str, verbose: bool = True) -> Dict:
    """
    Convenience function: align a single sequence to CARLIN reference sequence
    
    Args:
        sequence: Query sequence
        verbose: Whether to show detailed information
        
    Returns:
        Dict: Alignment result
    """
    aligner = create_default_aligner()
    return aligner.align_sequence(sequence, verbose=verbose) 