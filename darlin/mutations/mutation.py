#!/usr/bin/env python3
"""
CARLIN mutation annotation module

Implements mutation event identification, classification and HGVS format annotation functionality
"""

from typing import List, Optional, Dict, Tuple, Union
from dataclasses import dataclass
from enum import Enum
import re

from ..alignment.aligned_seq import AlignedSEQ, AlignedSEQMotif


class MutationType(Enum):
    """Mutation type enumeration"""
    SUBSTITUTION = "M"  # Substitution/mismatch
    INSERTION = "I"     # Insertion
    DELETION = "D"      # Deletion 
    COMPLEX = "C"       # Complex mutation (contains multiple types)
    INDEL = "DI"        # Deletion-insertion
    CONSERVED = "N"     # Conserved (no change)


@dataclass
class Mutation:
    """
    Mutation event class
    
    Represents a single mutation event with HGVS format annotation support
    
    Attributes:
        type: Mutation type (MutationType)
        loc_start: Mutation start position (1-based)
        loc_end: Mutation end position (1-based, inclusive)
        seq_old: Original sequence fragment
        seq_new: Mutated sequence fragment
        motif_index: Index of motif where mutation occurred (0-based)
        confidence: Confidence of mutation call (0-1)
        confidence_label: Discrete confidence label ('High' | 'Low')
    """
    type: MutationType
    loc_start: int
    loc_end: int
    seq_old: str
    seq_new: str
    motif_index: int = -1
    confidence: float = 1.0
    confidence_label: str = "Low"
    
    def __post_init__(self):
        """Validate mutation data"""
        if self.loc_start < 1:
            raise ValueError("Mutation position must be 1-based positive integer")
        if self.loc_end < self.loc_start:
            raise ValueError("End position cannot be less than start position")
        if self.confidence < 0 or self.confidence > 1:
            raise ValueError("Confidence must be between 0-1")
    
    @property
    def length_change(self) -> int:
        """Return length change caused by mutation"""
        return len(self.seq_new) - len(self.seq_old)
    
    @property
    def is_indel(self) -> bool:
        """Check if mutation is an indel"""
        return self.type in [MutationType.INSERTION, MutationType.DELETION, MutationType.INDEL]
    
    @property
    def affected_length(self) -> int:
        """Return affected sequence length"""
        return self.loc_end - self.loc_start + 1
    
    def to_hgvs(self, reference_name: str = "CARLIN") -> str:
        """
        Convert to HGVS format annotation
        
        Args:
            reference_name: Reference sequence name
            
        Returns:
            str: HGVS format mutation annotation
            
        Examples:
            - Substitution: g.100A>T
            - Deletion: g.100_102del
            - Insertion: g.100_101insACG
            - Complex: g.100_102delinsATG
        """
        # prefix = f"{reference_name}:g."
        prefix = ""
        
        if self.type == MutationType.SUBSTITUTION:
            if len(self.seq_old) == 1 and len(self.seq_new) == 1:
                # Single base substitution
                return f"{prefix}{self.loc_start}{self.seq_old}>{self.seq_new}"
            else:
                # Multi-base substitution 
                return f"{prefix}{self.loc_start}_{self.loc_end}{self.seq_old}>{self.seq_new}"
        
        elif self.type == MutationType.DELETION:
            # Always use range format for deletions (e.g., 2_2del for single base, 2_5del for multi-base)
            return f"{prefix}{self.loc_start}_{self.loc_end}del"
        
        elif self.type == MutationType.INSERTION:
            # Insertion mutation
            if self.loc_start == self.loc_end:
                return f"{prefix}{self.loc_start}_{self.loc_start + 1}ins{self.seq_new}"
            else:
                return f"{prefix}{self.loc_start}_{self.loc_end}ins{self.seq_new}"
        
        elif self.type in [MutationType.COMPLEX, MutationType.INDEL]:
            # Complex mutation or indel
            if self.seq_new:
                return f"{prefix}{self.loc_start}_{self.loc_end}delins{self.seq_new}"
            else:
                return f"{prefix}{self.loc_start}_{self.loc_end}del"
        
        elif self.type == MutationType.CONSERVED:
            return f"{prefix}{self.loc_start}_{self.loc_end}="
        
        else:
            return f"{prefix}{self.loc_start}_{self.loc_end}unknown"
    
    def __str__(self) -> str:
        """String representation"""
        return self.to_hgvs()
    
    def __repr__(self) -> str:
        """Detailed representation"""
        return (f"Mutation(type={self.type.value}, "
                f"pos={self.loc_start}-{self.loc_end}, "
                f"'{self.seq_old}'->'{self.seq_new}', "
                f"motif={self.motif_index})")


class MutationIdentifier:
    """
    Mutation identifier
    
    Identifies and annotates mutation events from AlignedSEQ objects
    """
    
    def __init__(self, min_confidence: float = 0.8):
        """
        Initialize mutation identifier
        
        Args:
            min_confidence: Minimum confidence threshold
        """
        self.min_confidence = min_confidence
    
    def identify_sequence_events(self, aligned_seq: AlignedSEQ) -> List[Mutation]:
        """
        Identify basic mutation events in sequence
        
        Args:
            aligned_seq: Aligned sequence object
            
        Returns:
            List[Mutation]: List of identified mutations
        """
        mutations = []
        position = 1  # 1-based position counting
        
        for motif_idx, motif in enumerate(aligned_seq.motifs):
            mutation = self._identify_motif_mutation(motif, position, motif_idx)
            if mutation:
                mutations.append(mutation)
            
            # Update position count (based on reference sequence)
            position += len(motif.ref.replace('-', ''))
        
        return mutations
    
    def _identify_motif_mutation(self, motif: AlignedSEQMotif, start_pos: int, motif_idx: int) -> Optional[Mutation]:
        """
        Identify mutations in a single motif
        
        Args:
            motif: Motif object
            start_pos: Start position of motif in reference sequence
            motif_idx: Motif index
            
        Returns:
            Optional[Mutation]: Identified mutation, returns None if no mutation
        """
        seq = motif.seq
        ref = motif.ref
        
        # First perform gap-based indel identification to avoid misclassification of insertions/deletions after gap removal
        seq_nogap = seq.replace('-', '')
        ref_nogap = ref.replace('-', '')

        if seq == ref:
            # Perfect match, no mutation
            return None

        # Check for complex mutation: large length change with mismatches
        seq_clean = seq.replace('-', '')
        ref_clean = ref.replace('-', '')
        length_change = len(seq_clean) - len(ref_clean)
        is_large_length_change = abs(length_change) > 3
        
        # Check for mismatches by comparing aligned positions
        has_mismatch = False
        min_len = min(len(seq), len(ref))
        for i in range(min_len):
            if seq[i] != ref[i] and seq[i] != '-' and ref[i] != '-':
                has_mismatch = True
                break
        
        # Complex mutation: large length change with mismatches
        if is_large_length_change and has_mismatch:
            mutation_type = MutationType.COMPLEX
            loc_start = start_pos
            loc_end = start_pos + len(ref_clean) - 1
            ref_nogap = ref_clean
            seq_nogap = seq_clean
        elif ('-' in ref) and ('-' not in seq):
            # Insertion: at positions where reference has gaps, query sequence provides inserted bases
            mutation_type = MutationType.INSERTION
            # Find exact inserted fragment and insertion position (between reference bases before and after gap)
            inserted_bases = []
            ref_bases_before_gap = 0
            seen_gap = False
            for i in range(len(ref)):
                if ref[i] == '-':
                    inserted_bases.append(seq[i])
                    seen_gap = True
                else:
                    if not seen_gap:
                        ref_bases_before_gap += 1
            # Insertion occurs between the reference base before the first gap and the reference base after it
            insertion_after = start_pos + max(ref_bases_before_gap - 1, 0)
            loc_start = insertion_after
            loc_end = insertion_after  # Insertion uses start==end representation
            ref_nogap = ''
            seq_nogap = ''.join(inserted_bases)
        elif ('-' in seq) and ('-' not in ref):
            # Deletion or delins: query has gaps at reference base positions
            # Find continuous mutation region: deletions and/or adjacent substitutions
            ref_bases_in_region = []
            seq_bases_in_region = []
            first_mut_index = None
            last_mut_index = None
            
            # Step 1: Find the start of mutation region (first gap)
            for i in range(len(seq)):
                if seq[i] == '-':
                    if first_mut_index is None:
                        first_mut_index = i
                    last_mut_index = i
            
            # Step 2: Extend the region to include adjacent substitutions
            # Check positions after the last gap for substitutions
            # Continue until we find a continuous match region (at least 2 consecutive matches)
            if last_mut_index is not None:
                consecutive_matches = 0
                for i in range(last_mut_index + 1, len(seq)):
                    if ref[i] != '-' and seq[i] != '-':
                        if ref[i] == seq[i]:
                            consecutive_matches += 1
                            # If we have 2 consecutive matches, stop extending
                            if consecutive_matches >= 2:
                                break
                        else:
                            # Found a substitution, reset match counter and extend region
                            consecutive_matches = 0
                            last_mut_index = i
                    else:
                        # Gap or other case, reset match counter
                        consecutive_matches = 0
            
            # Step 3: Also check for substitutions within the gap region
            # (handles cases like ref='ACGA', seq='--GT' where gaps and substitutions are interleaved)
            if first_mut_index is not None:
                for i in range(first_mut_index, min(last_mut_index + 1 if last_mut_index is not None else len(seq), len(seq))):
                    if ref[i] != '-' and seq[i] != '-' and ref[i] != seq[i]:
                        # Found substitution within or adjacent to deletion region
                        if last_mut_index is None or i > last_mut_index:
                            last_mut_index = i
            
            # Step 4: Collect all reference and query bases in the mutation region
            if first_mut_index is not None:
                for i in range(first_mut_index, (last_mut_index + 1) if last_mut_index is not None else len(seq)):
                    if ref[i] != '-':
                        ref_bases_in_region.append(ref[i])
                    if seq[i] != '-':
                        seq_bases_in_region.append(seq[i])
            
            # Step 5: Determine mutation type
            has_substitution = False
            if first_mut_index is not None:
                for i in range(first_mut_index, (last_mut_index + 1) if last_mut_index is not None else len(seq)):
                    if ref[i] != '-' and seq[i] != '-' and ref[i] != seq[i]:
                        has_substitution = True
                        break
            
            # Check if ref and seq sequences differ (indicating delins)
            ref_str = ''.join(ref_bases_in_region)
            seq_str = ''.join(seq_bases_in_region)
            
            if has_substitution or (ref_str and seq_str and ref_str != seq_str):
                # This is a delins (deletion-insertion) event
                mutation_type = MutationType.INDEL
                ref_nogap = ref_str
                seq_nogap = seq_str
            else:
                # Pure deletion
                mutation_type = MutationType.DELETION
                ref_nogap = ref_str
                seq_nogap = ''
            
            # Step 6: Calculate start and end positions
            if first_mut_index is not None:
                ref_non_gap_before = sum(1 for i in range(first_mut_index) if ref[i] != '-')
                loc_start = start_pos + ref_non_gap_before
                if ref_nogap:
                    loc_end = loc_start + len(ref_nogap) - 1
                else:
                    loc_end = loc_start
            else:
                # Fallback (should not happen)
                loc_start = start_pos
                loc_end = start_pos
        elif ('-' in ref) and ('-' in seq):
            # Both ref and seq have gaps - this is a delins (deletion-insertion) event
            # Find the mutation region (where ref and seq differ or have gaps)
            first_mut_index = None
            last_mut_index = None
            
            # Find the start and end of mutation region
            for i in range(len(ref)):
                if ref[i] != seq[i] or ref[i] == '-' or seq[i] == '-':
                    if first_mut_index is None:
                        first_mut_index = i
                    last_mut_index = i
            
            # Collect reference and query bases in the mutation region only
            ref_bases = []
            seq_bases = []
            
            if first_mut_index is not None:
                for i in range(first_mut_index, last_mut_index + 1):
                    if ref[i] != '-':
                        ref_bases.append(ref[i])
                    if seq[i] != '-':
                        seq_bases.append(seq[i])
            
            ref_nogap = ''.join(ref_bases)
            seq_nogap = ''.join(seq_bases)
            
            # This is a delins event
            mutation_type = MutationType.INDEL
            
            # Calculate position: count ref bases before mutation region
            if first_mut_index is not None:
                ref_non_gap_before = sum(1 for i in range(first_mut_index) if ref[i] != '-')
                loc_start = start_pos + ref_non_gap_before
                if ref_nogap:
                    loc_end = loc_start + len(ref_nogap) - 1
                else:
                    loc_end = loc_start
            else:
                loc_start = start_pos
                loc_end = start_pos
        else:
            # Other cases fall back to length/content comparison classification
            # But first check if this should be complex mutation based on both length change and content change
            seq_clean = seq.replace('-', '')
            ref_clean = ref.replace('-', '')
            length_change = len(seq_clean) - len(ref_clean)
            
            # Check if there's both length change > 3 AND content change (not pure insertion/deletion)
            has_content_change = seq_clean != ref_clean
            is_large_length_change = abs(length_change) > 3
            
            if is_large_length_change and has_content_change:
                # Complex mutation: large length change with content change
                mutation_type = MutationType.COMPLEX
                loc_start = start_pos
                loc_end = start_pos + len(ref_clean) - 1
                ref_nogap = ref_clean
                seq_nogap = seq_clean
            else:
                # Use standard classification
                mutation_type, seq_old, seq_new = self._classify_mutation(seq, ref)
                # Use _classify_mutation results to override no-gap sequences
                ref_nogap = seq_old if isinstance(seq_old, str) else ref_clean
                seq_nogap = seq_new if isinstance(seq_new, str) else seq_clean
        
        # Calculate precise position range
        if 'loc_start' not in locals():
            # For non-specialized branches (or complex/substitution) position by motif overall reference length
            loc_start = start_pos
            if len(ref_nogap) > 0:
                loc_end = start_pos + len(ref_nogap) - 1
            else:
                # For insertion case, end position equals start position
                loc_end = start_pos

        # HGVS-style normalization: shrink shared prefix/suffix only for
        # complex/indel (delins) events so that we always use the minimal
        # representation, but keep pure substitutions unchanged.
        if mutation_type in (MutationType.COMPLEX, MutationType.INDEL) and ref_nogap and seq_nogap:
            # Trim common prefix
            prefix_len = 0
            max_pref = min(len(ref_nogap), len(seq_nogap))
            while prefix_len < max_pref and ref_nogap[prefix_len] == seq_nogap[prefix_len]:
                prefix_len += 1
            if prefix_len > 0:
                ref_nogap = ref_nogap[prefix_len:]
                seq_nogap = seq_nogap[prefix_len:]
                loc_start += prefix_len

            # Trim common suffix, but only if both sequences have remaining content
            # This prevents incorrectly removing suffixes when one sequence is a substring of the other
            if ref_nogap and seq_nogap:
                suffix_len = 0
                max_suf = min(len(ref_nogap), len(seq_nogap))
                while (
                    suffix_len < max_suf
                    and ref_nogap[-(suffix_len + 1)] == seq_nogap[-(suffix_len + 1)]
                ):
                    suffix_len += 1
                # Only trim suffix if it doesn't completely consume the shorter sequence
                # This prevents cases like 'GA' vs 'A' where 'A' is a suffix of 'GA'
                if suffix_len > 0 and suffix_len < min(len(ref_nogap), len(seq_nogap)):
                    if suffix_len < len(ref_nogap):
                        ref_nogap = ref_nogap[:-suffix_len]
                    else:
                        ref_nogap = ""
                    if suffix_len < len(seq_nogap):
                        seq_nogap = seq_nogap[:-suffix_len]
                    else:
                        seq_nogap = ""
                    loc_end -= suffix_len

            # Ensure loc_end is not before loc_start (can happen if everything matched)
            if loc_end < loc_start:
                loc_end = loc_start
        
        return Mutation(
            type=mutation_type,
            loc_start=loc_start,
            loc_end=loc_end,
            seq_old=ref_nogap,
            seq_new=seq_nogap,
            motif_index=motif_idx,
            confidence=self.min_confidence,
            confidence_label="Low"
        )
    
    def _classify_mutation(self, seq: str, ref: str) -> Tuple[MutationType, str, str]:
        """
        Classify mutation type
        
        Args:
            seq: Query sequence (with gaps)
            ref: Reference sequence (with gaps)
            
        Returns:
            Tuple[MutationType, str, str]: (mutation type, original sequence, new sequence)
        """
        seq_clean = seq.replace('-', '')
        ref_clean = ref.replace('-', '')
        
        if not seq_clean and ref_clean:
            # Pure deletion
            return MutationType.DELETION, ref_clean, ''
        elif seq_clean and not ref_clean:
            # Pure insertion  
            return MutationType.INSERTION, '', seq_clean
        elif len(seq_clean) == len(ref_clean):
            # Equal length, substitution
            return MutationType.SUBSTITUTION, ref_clean, seq_clean
        else:
            # Unequal length, complex mutation
            if abs(len(seq_clean) - len(ref_clean)) <= 3:
                return MutationType.INDEL, ref_clean, seq_clean
            else:
                return MutationType.COMPLEX, ref_clean, seq_clean
    
    def identify_cas9_events(self, aligned_seq: AlignedSEQ, cut_sites: Optional[List[int]] = None) -> List[Mutation]:
        """
        Identify Cas9-specific editing events
        
        Args:
            aligned_seq: Aligned sequence object
            cut_sites: Cas9 cut site list (1-based positions)
            
        Returns:
            List[Mutation]: Cas9-specific mutation list
        """
        # First get basic mutation events
        basic_mutations = self.identify_sequence_events(aligned_seq)
        
        # If no cut sites provided, return basic mutations
        if not cut_sites:
            return basic_mutations
        
        # Filter and enhance Cas9-specific events
        cas9_mutations = []
        window = 4

        # Support cut_sites as 1-based integer list or interval list (inclusive endpoints)
        def _normalize_intervals(cuts: List[Union[int, Tuple[int, int]]]) -> List[Tuple[int, int]]:
            intervals: List[Tuple[int, int]] = []
            for cs in cuts:
                if isinstance(cs, int):
                    intervals.append((cs, cs))
                else:
                    start, end = cs
                    if start > end:
                        start, end = end, start
                    intervals.append((start, end))
            return intervals

        cut_intervals: List[Tuple[int, int]] = _normalize_intervals(cut_sites) if cut_sites else []
        for mutation in basic_mutations:
            # Use interval overlap to determine near: expand mutation interval by window and check intersection with cutsite intervals
            mut_start = mutation.loc_start
            mut_end = mutation.loc_start if mutation.type == MutationType.INSERTION else mutation.loc_end
            expanded_start = mut_start - window
            expanded_end = mut_end + window
            near_cutsite = any(not (expanded_end < cs_start or expanded_start > cs_end)
                               for cs_start, cs_end in cut_intervals)
            
            if near_cutsite:
                # Increase confidence
                mutation.confidence = min(1.0, mutation.confidence + 0.1)
                mutation.confidence_label = "High"
                cas9_mutations.append(mutation)
            elif mutation.type in [MutationType.INSERTION, MutationType.DELETION, MutationType.INDEL, MutationType.COMPLEX]:
                # Even if not near cut sites, indels and complex mutations could be Cas9-induced
                mutation.confidence = max(0.5, mutation.confidence - 0.2)
                mutation.confidence_label = "Low"
                cas9_mutations.append(mutation)
        
        return cas9_mutations
    
    def merge_adjacent_mutations(self, mutations: List[Mutation], max_distance: int = 3) -> List[Mutation]:
        """
        Merge adjacent mutation events
        
        Args:
            mutations: Mutation list
            max_distance: Maximum merge distance
            
        Returns:
            List[Mutation]: Merged mutation list
        """
        if not mutations:
            return []
        
        # Sort by position
        sorted_mutations = sorted(mutations, key=lambda m: m.loc_start)
        merged = [sorted_mutations[0]]
        
        for current in sorted_mutations[1:]:
            last = merged[-1]
            
            # Check if should merge
            if current.loc_start - last.loc_end <= max_distance:
                # Merge mutations
                merged_mutation = self._merge_two_mutations(last, current)
                merged[-1] = merged_mutation
            else:
                merged.append(current)
        
        return merged
    
    def _merge_two_mutations(self, mut1: Mutation, mut2: Mutation) -> Mutation:
        """
        Merge two mutation events
        
        Args:
            mut1: First mutation
            mut2: Second mutation
            
        Returns:
            Mutation: Merged mutation
        """
        # Determine merged position range
        loc_start = min(mut1.loc_start, mut2.loc_start)
        loc_end = max(mut1.loc_end, mut2.loc_end)
        
        # Merge sequence information
        seq_old = mut1.seq_old + mut2.seq_old
        seq_new = mut1.seq_new + mut2.seq_new
        
        # Determine merged type
        if mut1.type == mut2.type:
            merged_type = mut1.type
        else:
            merged_type = MutationType.COMPLEX
        
        # Average confidence
        avg_confidence = (mut1.confidence + mut2.confidence) / 2
        
        return Mutation(
            type=merged_type,
            loc_start=loc_start,
            loc_end=loc_end,
            seq_old=seq_old,
            seq_new=seq_new,
            motif_index=min(mut1.motif_index, mut2.motif_index),
            confidence=avg_confidence
        )


# Convenience function
def annotate_mutations(aligned_seq: AlignedSEQ, 
                      cas9_mode: bool = True,
                      cut_sites: Optional[List[int]] = None,
                      merge_adjacent: bool = True) -> List[Mutation]:
    """
    Annotate mutations for aligned sequence
    
    Args:
        aligned_seq: Aligned sequence object
        cas9_mode: Whether to use Cas9 mode
        cut_sites: Cas9 cut site list
        merge_adjacent: Whether to merge adjacent mutations
        
    Returns:
        List[Mutation]: Annotated mutation list
    """
    identifier = MutationIdentifier()
    
    if cas9_mode:
        mutations = identifier.identify_cas9_events(aligned_seq, cut_sites)
    else:
        mutations = identifier.identify_sequence_events(aligned_seq)
    
    if merge_adjacent:
        mutations = identifier.merge_adjacent_mutations(mutations)
    
    # After merging, recalculate near based on final intervals and set labels/scores (window=4), support cutsite as points or intervals
    if cas9_mode and cut_sites:
        window = 4

        def _normalize_intervals(cuts: List[Union[int, Tuple[int, int]]]) -> List[Tuple[int, int]]:
            intervals: List[Tuple[int, int]] = []
            for cs in cuts:
                if isinstance(cs, int):
                    intervals.append((cs, cs))
                else:
                    start, end = cs
                    if start > end:
                        start, end = end, start
                    intervals.append((start, end))
            return intervals

        cut_intervals: List[Tuple[int, int]] = _normalize_intervals(cut_sites)
        for m in mutations:
            mut_start = m.loc_start
            mut_end = m.loc_start if m.type == MutationType.INSERTION else m.loc_end
            expanded_start = mut_start - window
            expanded_end = mut_end + window
            near = any(not (expanded_end < cs_start or expanded_start > cs_end)
                       for cs_start, cs_end in cut_intervals)
            if near:
                m.confidence_label = "High"
                m.confidence = min(1.0, max(m.confidence, 0.8) + 0.1)
            elif m.type in [MutationType.INSERTION, MutationType.DELETION]:
                m.confidence_label = "Low"
                m.confidence = max(0.5, min(m.confidence, 1.0) - 0.2)

    return mutations 