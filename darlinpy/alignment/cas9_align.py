import numpy as np
from typing import Tuple, List
from enum import IntEnum

try:
    from . import _cas9_align as _cas9_align_module
except ImportError:  # pragma: no cover - optional acceleration
    _cas9_align_module = None

HAS_CPP_IMPL = _cas9_align_module is not None


class Mutation(IntEnum):
    S = 0  # Substitution
    D = 1  # Deletion
    I = 2  # Insertion
    N_MUTS = 3

def max3(vals: np.ndarray) -> Tuple[float, np.ndarray]:
    """
    Find the maximum value among three values and return which positions reach the maximum
    
    Args:
        vals: Array of length 3 [S_score, D_score, I_score]
    
    Returns:
        maxval: Maximum value
        argmax3: Boolean array of length 3 indicating which positions reach the maximum
    """
    maxval = np.max(vals)
    argmax3 = (vals == maxval)
    return maxval, argmax3

def cas9_align_py(seq: np.ndarray, ref: np.ndarray,
                  open_penalty: np.ndarray, close_penalty: np.ndarray,
                  sub_score: np.ndarray) -> Tuple[float, List[int], List[int]]:
    """
    CRISPR-Cas9 specific sequence alignment algorithm
    
    Args:
        seq: Query sequence encoded as integers (A=1, C=2, G=3, T=4, gap=0)
        ref: Reference sequence with same encoding format
        open_penalty: Position-specific gap opening penalty array, size Lref+1
        close_penalty: Position-specific gap closing penalty array, size Lref+1  
        sub_score: 25x1 substitution score matrix (for calculating match/mismatch scores)
                   Index calculation: score = sub_score[seq_nt * 5 + ref_nt]
    
    Returns:
        best_score: Optimal alignment score
        al_seq: Aligned sequence
        al_ref: Aligned reference sequence
    """
    
    Lseq = len(seq)
    Lref = len(ref)
    
    # Initialize three-state score matrix (S, D, I)
    score = np.full((Mutation.N_MUTS, Lseq + 1, Lref + 1), 
                    -np.inf, dtype=np.float64)
    
    # Initialize backtrack matrix [from_state][to_state][i][j]
    backtrack = np.zeros((Mutation.N_MUTS, Mutation.N_MUTS, Lseq + 1, Lref + 1), 
                        dtype=bool)
    
    # Initialize boundary conditions
    score[Mutation.S, 0, 0] = 0.0
    
    # Initialize first row (insertion sequence)
    for j in range(1, Lseq + 1):
        score[Mutation.I, j, 0] = -open_penalty[0]
        backtrack[Mutation.I, Mutation.I, j, 0] = True
    
    # Initialize first column (deletion sequence) 
    for k in range(1, Lref + 1):
        score[Mutation.D, 0, k] = -open_penalty[0] 
        backtrack[Mutation.D, Mutation.D, 0, k] = True
    
    # Fill dynamic programming matrix
    for j in range(1, Lseq + 1):
        for k in range(1, Lref + 1):
            # Calculate substitution score (note: sequences are 1-indexed encoded)
            ss = sub_score[seq[j-1] * 5 + ref[k-1]]
            
            # Build state transition score matrix [to_state][from_state]
            score_jk = np.array([
                # Transition to S (substitution)
                [score[Mutation.S, j-1, k-1] + ss,
                 score[Mutation.D, j-1, k-1] - close_penalty[k-1] + ss,
                 score[Mutation.I, j-1, k-1] - close_penalty[k] + ss],
                
                # Transition to D (deletion)  
                [score[Mutation.S, j, k-1] - open_penalty[k],
                 score[Mutation.D, j, k-1],
                 score[Mutation.I, j, k-1]],
                
                # Transition to I (insertion)
                [score[Mutation.S, j-1, k] - open_penalty[k],
                 score[Mutation.D, j-1, k],
                 score[Mutation.I, j-1, k]]
            ])
            
            # Find optimal predecessor state for each state
            for m in range(Mutation.N_MUTS):
                score[m, j, k], out = max3(score_jk[m])
                backtrack[Mutation.S, m, j, k] = out[Mutation.S]
                backtrack[Mutation.D, m, j, k] = out[Mutation.D]  
                backtrack[Mutation.I, m, j, k] = out[Mutation.I]
    
    # Find final optimal state
    score_LL = np.array([score[Mutation.S, Lseq, Lref],
                        score[Mutation.D, Lseq, Lref], 
                        score[Mutation.I, Lseq, Lref]])
    
    best_score, state_LL = max3(score_LL)
    
    # Determine optimal end state
    if state_LL[Mutation.S]:
        cur_state = Mutation.S
    elif state_LL[Mutation.D]:
        cur_state = Mutation.D
    elif state_LL[Mutation.I]:
        cur_state = Mutation.I
    else:
        cur_state = -1
    
    # Backtrack to reconstruct alignment
    al_seq = []
    al_ref = []
    cur_j = Lseq
    cur_k = Lref
    
    while cur_j > 0 or cur_k > 0:
        assert 0 <= cur_state < Mutation.N_MUTS
        
        # Get backtrack information for current state
        cur_bt = np.array([
            backtrack[Mutation.S, cur_state, cur_j, cur_k],
            backtrack[Mutation.D, cur_state, cur_j, cur_k], 
            backtrack[Mutation.I, cur_state, cur_j, cur_k]
        ])
        
        if cur_state == Mutation.S:
            # Substitution: consume one character from both sequence and reference
            al_seq.append(seq[cur_j - 1])
            al_ref.append(ref[cur_k - 1])
            cur_j -= 1
            cur_k -= 1
            # Choose predecessor state
            if cur_bt[Mutation.S]:
                cur_state = Mutation.S
            elif cur_bt[Mutation.D]:
                cur_state = Mutation.D
            elif cur_bt[Mutation.I]:
                cur_state = Mutation.I
            else:
                cur_state = -1
                
        elif cur_state == Mutation.D:
            # Deletion: only consume one character from reference sequence
            al_seq.append(0)  # gap
            al_ref.append(ref[cur_k - 1])
            cur_k -= 1
            # Choose predecessor state (priority: D > S > I)
            if cur_bt[Mutation.D]:
                cur_state = Mutation.D
            elif cur_bt[Mutation.S]:
                cur_state = Mutation.S
            elif cur_bt[Mutation.I]:
                cur_state = Mutation.I
            else:
                cur_state = -1
                
        elif cur_state == Mutation.I:
            # Insertion: only consume one character from sequence
            al_seq.append(seq[cur_j - 1])
            al_ref.append(0)  # gap
            cur_j -= 1
            # Choose predecessor state (priority: I > S > D)
            if cur_bt[Mutation.I]:
                cur_state = Mutation.I
            elif cur_bt[Mutation.S]:
                cur_state = Mutation.S
            elif cur_bt[Mutation.D]:
                cur_state = Mutation.D
            else:
                cur_state = -1
        else:
            cur_state = -1
    
    # Reverse to get correctly ordered alignment
    al_seq.reverse()
    al_ref.reverse()
    
    return best_score, al_seq, al_ref


def cas9_align(seq: np.ndarray, ref: np.ndarray,
               open_penalty: np.ndarray, close_penalty: np.ndarray,
               sub_score: np.ndarray) -> Tuple[float, List[int], List[int]]:
    """Dispatch to the fastest available cas9 alignment implementation."""

    if HAS_CPP_IMPL:
        return _cas9_align_module.cas9_align(seq, ref, open_penalty, close_penalty, sub_score)
    return cas9_align_py(seq, ref, open_penalty, close_penalty, sub_score)

def nt2int(sequence: str) -> np.ndarray:
    """Convert DNA sequence to integer encoding"""
    mapping = {'A': 1, 'C': 2, 'G': 3, 'T': 4}
    return np.array([mapping.get(nt.upper(), 0) for nt in sequence], dtype=np.uint8)

def int2nt(int_seq: List[int]) -> str:
    """Convert integer encoding back to DNA sequence"""
    mapping = {0: '-', 1: 'A', 2: 'C', 3: 'G', 4: 'T'}
    return ''.join(mapping.get(nt, 'N') for nt in int_seq)

def print_cas9_alignment(al_seq: List[int], al_ref: List[int], score: float):
    """Print alignment results"""
    seq_str = int2nt(al_seq)
    ref_str = int2nt(al_ref)
    
    print(f"Sequence: {seq_str}")
    print(f"Reference: {ref_str}")
    print(f"Score: {score:.2f}")
    
    # Print match information
    match_line = ""
    for i in range(len(al_seq)):
        if al_seq[i] == al_ref[i] and al_seq[i] != 0:
            match_line += "|"
        elif al_seq[i] == 0 or al_ref[i] == 0:
            match_line += " "
        else:
            match_line += "."
    print(f"Match: {match_line}")

# Example usage
if __name__ == "__main__":
    # Example sequences
    seq_str = "ACGTACGT"
    ref_str = "ACGTGCGT"
    
    # Convert to integer encoding
    seq = nt2int(seq_str)
    ref = nt2int(ref_str)
    
    # Set simple penalty parameters (should be adjusted based on CRISPR cut sites in actual use)
    open_penalty = np.full(len(ref_str) + 1, 2.0)  # Uniform opening penalty
    close_penalty = np.full(len(ref_str) + 1, 1.0)     # Uniform closing penalty
    
    # Set simple NUC44 scoring matrix (should use complete 25x25 matrix in practice)
    sub_score = np.zeros(25, dtype=np.float64)
    # Simplified match/mismatch scores
    for i in range(1, 5):
        for j in range(1, 5):
            if i == j:
                sub_score[i * 5 + j] = 2.0  # Match
            else:
                sub_score[i * 5 + j] = -1.0  # Mismatch
    
    # Execute alignment
    score, al_seq, al_ref = cas9_align(seq, ref, open_penalty, close_penalty, sub_score)
    
    # Print results
    print("=== Example 1: Simple sequence alignment ===")
    print_cas9_alignment(al_seq, al_ref, score)
    
    # Example 2: Sequence with insertions and deletions
    print("\n=== Example 2: Sequence with insertions and deletions ===")
    seq_str2 = "ACGTGCGT"
    ref_str2 = "ACGTACCCGT"
    
    seq2 = nt2int(seq_str2)
    ref2 = nt2int(ref_str2)
    
    # Adjust penalty parameter sizes
    open_penalty2 = np.full(len(ref_str2) + 1, 2.0)
    close_penalty2 = np.full(len(ref_str2) + 1, 1.0)
    
    # Set lower gap penalties around position 5 to simulate CRISPR cut sites
    open_penalty2[4:7] = 0.5
    close_penalty2[4:7] = 0.2
    
    score2, al_seq2, al_ref2 = cas9_align(seq2, ref2, open_penalty2, close_penalty2, sub_score)
    print_cas9_alignment(al_seq2, al_ref2, score2)
    
    print(f"\nNote: Lower gap penalties set at positions 4-6 ({open_penalty2[5]:.1f}),")
    print(f"simulating CRISPR cut sites, with gap penalties of {open_penalty2[0]:.1f} at other positions")