#!/usr/bin/env python3
"""
CARLIN alignment sequence data structures and normalization functionality

Implements AlignedSEQ, AlignedSEQMotif classes and normalization algorithms for prefix/postfix and conserved regions
"""

from typing import List, Tuple, Optional, Union
import numpy as np


class AlignedSEQMotif:
    """
    Aligned sequence motif fragment
    
    Represents a motif fragment in an aligned sequence, containing sequence, reference sequence and event type
    """
    
    def __init__(self, seq: str, ref: str):
        """
        Initialize motif
        
        Args:
            seq: Query sequence fragment
            ref: Reference sequence fragment
        """
        if len(seq) != len(ref):
            raise ValueError(f"Sequence length ({len(seq)}) does not match reference sequence length ({len(ref)})")
        
        self.seq = seq
        self.ref = ref
        self.event = self._classify_motif_event(seq, ref)
    
    def _classify_motif_event(self, seq: str, ref: str) -> str:
        """
        Classify motif event type
        
        Args:
            seq: Query sequence
            ref: Reference sequence
            
        Returns:
            str: Event type
                'N' - No change
                'M' - Mismatch
                'D' - Deletion
                'I' - Insertion
                'E' - Empty
        """
        if seq == ref:
            return 'N'
        elif all(c == '-' for c in seq):
            return 'E'
        elif all(c != '-' for c in ref):
            if any(c == '-' for c in seq):
                return 'D'
            else:
                return 'M'
        else:
            return 'I'
    
    def __repr__(self):
        return f"AlignedSEQMotif(seq='{self.seq}', ref='{self.ref}', event='{self.event}')"


class AlignedSEQ:
    """
    Complete aligned sequence object
    
    Composed of multiple AlignedSEQMotif objects, representing a complete aligned sequence
    """
    
    def __init__(self, seq_segments: List[str], ref_segments: List[str]):
        """
        Initialize aligned sequence
        
        Args:
            seq_segments: List of motif fragments for query sequence
            ref_segments: List of motif fragments for reference sequence
        """
        if len(seq_segments) != len(ref_segments):
            raise ValueError(f"Number of sequence segments ({len(seq_segments)}) does not match number of reference segments ({len(ref_segments)})")
        
        self.motifs = [AlignedSEQMotif(seq, ref) for seq, ref in zip(seq_segments, ref_segments)]
    
    def get_seq(self) -> str:
        """Get complete query sequence"""
        return ''.join(motif.seq for motif in self.motifs)
    
    def get_ref(self) -> str:
        """Get complete reference sequence"""
        return ''.join(motif.ref for motif in self.motifs)
    
    def get_event_structure(self) -> List[str]:
        """Get event structure"""
        return [motif.event for motif in self.motifs]
    
    def copy(self) -> 'AlignedSEQ':
        """Create deep copy"""
        seq_segments = [motif.seq for motif in self.motifs]
        ref_segments = [motif.ref for motif in self.motifs]
        return AlignedSEQ(seq_segments, ref_segments)
    
    def __repr__(self):
        events = ''.join(self.get_event_structure())
        return f"AlignedSEQ(motifs={len(self.motifs)}, events='{events}')"


class SequenceSanitizer:
    """
    Sequence normalizer
    
    Implements normalization for CARLIN sequence prefix/postfix and conserved regions
    """
    
    @staticmethod
    def sanitize_prefix_postfix(aligned_seqs: Union[AlignedSEQ, List[AlignedSEQ]]) -> Union[AlignedSEQ, List[AlignedSEQ]]:
        """
        Clean up non-functional insertions in prefix and suffix
        
        Remove insertions in CARLIN prefix and suffix regions, which are usually PCR artifacts or sequencing errors
        
        Args:
            aligned_seqs: Single AlignedSEQ or list of AlignedSEQ
            
        Returns:
            Normalized AlignedSEQ or list
        """
        # Handle single object case
        is_single = False
        if isinstance(aligned_seqs, AlignedSEQ):
            is_single = True
            aligned_seqs = [aligned_seqs]
        
        if not aligned_seqs:
            return aligned_seqs
        
        result_seqs = []
        
        for aligned_seq in aligned_seqs:
            event_structure = aligned_seq.get_event_structure()
            motifs = aligned_seq.motifs.copy()
            
            # Handle head insertions (prefix)
            if event_structure and event_structure[0] == 'I':
                first_motif = motifs[0]
                ref_seq = first_motif.ref
                
                # Find position of first non-gap character in reference sequence
                first_non_gap = None
                for i, char in enumerate(ref_seq):
                    if char != '-':
                        first_non_gap = i
                        break
                
                # If there are gaps at the beginning, trim them
                if first_non_gap is not None and first_non_gap > 0:
                    trimmed_seq = first_motif.seq[first_non_gap:]
                    trimmed_ref = first_motif.ref[first_non_gap:]
                    motifs[0] = AlignedSEQMotif(trimmed_seq, trimmed_ref)
            
            # Handle tail insertions (postfix)
            if event_structure and event_structure[-1] == 'I':
                last_motif = motifs[-1]
                ref_seq = last_motif.ref
                
                # Find position of last non-gap character in reference sequence
                last_non_gap = None
                for i in range(len(ref_seq) - 1, -1, -1):
                    if ref_seq[i] != '-':
                        last_non_gap = i
                        break
                
                # If there are gaps at the end, trim them
                if last_non_gap is not None and last_non_gap < len(ref_seq) - 1:
                    trimmed_seq = last_motif.seq[:last_non_gap + 1]
                    trimmed_ref = last_motif.ref[:last_non_gap + 1]
                    motifs[-1] = AlignedSEQMotif(trimmed_seq, trimmed_ref)
            
            # Rebuild AlignedSEQ
            seq_segments = [motif.seq for motif in motifs]
            ref_segments = [motif.ref for motif in motifs]
            result_seqs.append(AlignedSEQ(seq_segments, ref_segments))
        
        return result_seqs[0] if is_single else result_seqs
    
    @staticmethod
    def sanitize_conserved_regions(aligned_seqs: Union[AlignedSEQ, List[AlignedSEQ]], 
                                 cutsite_motif_indices: List[int]) -> Union[AlignedSEQ, List[AlignedSEQ]]:
        """
        Clean up sequencing errors in conserved regions
        
        Restore mismatches in non-cutsite motifs to reference sequence, as these regions are biologically conserved
        
        Args:
            aligned_seqs: Single AlignedSEQ or list of AlignedSEQ
            cutsite_motif_indices: List of cutsite motif indices
            
        Returns:
            Normalized AlignedSEQ or list
        """
        # Handle single object case
        is_single = False
        if isinstance(aligned_seqs, AlignedSEQ):
            is_single = True
            aligned_seqs = [aligned_seqs]
        
        if not aligned_seqs:
            return aligned_seqs
        
        result_seqs = []
        
        for aligned_seq in aligned_seqs:
            event_structure = aligned_seq.get_event_structure()
            motifs = aligned_seq.motifs.copy()
            
            # Original sequence length (for validation)
            original_length = len(aligned_seq.get_seq())
            
            # Identify motifs that need cleaning
            # Only clean mismatches ('M') that are not cutsites
            motifs_to_clean = []
            for i, (event, motif) in enumerate(zip(event_structure, motifs)):
                if event == 'M' and i not in cutsite_motif_indices:
                    # Ensure this motif has no gaps
                    if '-' not in motif.seq and '-' not in motif.ref:
                        motifs_to_clean.append(i)
            
            # Execute cleaning: restore mismatched sequences to reference sequence
            for i in motifs_to_clean:
                old_motif = motifs[i]
                # Replace query sequence with reference sequence
                motifs[i] = AlignedSEQMotif(old_motif.ref, old_motif.ref)
            
            # Rebuild AlignedSEQ
            seq_segments = [motif.seq for motif in motifs]
            ref_segments = [motif.ref for motif in motifs]
            cleaned_seq = AlignedSEQ(seq_segments, ref_segments)
            
            # Validate that sequence length has not changed
            new_length = len(cleaned_seq.get_seq())
            if new_length != original_length:
                raise RuntimeError(f"Sequence length changed after normalization: {original_length} -> {new_length}")
            
            result_seqs.append(cleaned_seq)
        
        return result_seqs[0] if is_single else result_seqs


def desemble_sequence(aligned_query: str, aligned_ref: str, 
                     motif_boundaries: List[Tuple[int, int]]) -> AlignedSEQ:
    """
    Decompose alignment results into motifs
    
    Args:
        aligned_query: Aligned query sequence
        aligned_ref: Aligned reference sequence
        motif_boundaries: List of motif boundaries, each element is (start, end)
        
    Returns:
        AlignedSEQ object
    """
    if len(aligned_query) != len(aligned_ref):
        raise ValueError(f"Aligned sequence lengths do not match: {len(aligned_query)} vs {len(aligned_ref)}")
    
    seq_segments = []
    ref_segments = []
    
    for start, end in motif_boundaries:
        if start < 0 or end > len(aligned_query) or start >= end:
            raise ValueError(f"Invalid motif boundary: ({start}, {end})")
        
        seq_segments.append(aligned_query[start:end])
        ref_segments.append(aligned_ref[start:end])
    
    result = AlignedSEQ(seq_segments, ref_segments)
    
    # Validate decomposition results
    if result.get_seq() != aligned_query:
        raise RuntimeError("Decomposed sequence does not match original sequence")
    if result.get_ref() != aligned_ref:
        raise RuntimeError("Decomposed reference sequence does not match original sequence")
    
    return result


def calculate_motif_boundaries(aligned_ref: str, carlin_config) -> List[Tuple[int, int]]:
    """
    Calculate motif boundaries
    
    Calculate boundaries for each motif based on aligned reference sequence and CARLIN configuration
    
    Args:
        aligned_ref: Aligned reference sequence
        carlin_config: CARLIN configuration object
        
    Returns:
        List of motif boundaries
    """
    # Find non-gap positions in reference sequence
    ref_positions = [i for i, char in enumerate(aligned_ref) if char != '-']
    
    if len(ref_positions) != len(carlin_config.carlin_sequence):
        raise ValueError(f"Number of non-gap positions in reference sequence ({len(ref_positions)}) does not match CARLIN sequence length ({len(carlin_config.carlin_sequence)})")
    
    # Get CARLIN internal motif boundaries (relative to original CARLIN sequence)
    carlin_boundaries = []
    
    # Prefix
    prefix_start, prefix_end = carlin_config.positions['prefix']
    carlin_boundaries.append((prefix_start, prefix_end))
    
    # Segments and PAMs
    for i in range(10):
        # Consite
        consite_start, consite_end = carlin_config.positions['consites'][i]
        carlin_boundaries.append((consite_start, consite_end))
        
        # Cutsite
        cutsite_start, cutsite_end = carlin_config.positions['cutsites'][i]
        carlin_boundaries.append((cutsite_start, cutsite_end))
        
        # PAM (first 9 segments have PAM after them)
        if i < 9:
            pam_start, pam_end = carlin_config.positions['pams'][i]
            carlin_boundaries.append((pam_start, pam_end))
    
    # Postfix
    postfix_start, postfix_end = carlin_config.positions['postfix']
    carlin_boundaries.append((postfix_start, postfix_end))
    
    # Map CARLIN boundaries to alignment sequence boundaries
    aligned_boundaries = []
    for carlin_start, carlin_end in carlin_boundaries:
        aligned_start = ref_positions[carlin_start]
        aligned_end = ref_positions[carlin_end - 1] + 1  # end is exclusive
        aligned_boundaries.append((aligned_start, aligned_end))
    
    # Handle gaps: ensure boundaries cover all aligned sequences
    # Adjust start of first boundary and end of last boundary
    if aligned_boundaries:
        aligned_boundaries[0] = (0, aligned_boundaries[0][1])
        aligned_boundaries[-1] = (aligned_boundaries[-1][0], len(aligned_ref))
        
        # Handle gaps in between
        for i in range(len(aligned_boundaries) - 1):
            current_end = aligned_boundaries[i][1]
            next_start = aligned_boundaries[i + 1][0]
            
            if current_end < next_start:
                # There is a gap, assign gap to next motif
                aligned_boundaries[i + 1] = (current_end, aligned_boundaries[i + 1][1])
    
    return aligned_boundaries 