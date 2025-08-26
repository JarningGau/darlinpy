import numpy as np
from typing import Tuple, List
from enum import IntEnum

class Mutation(IntEnum):
    S = 0  # Substitution
    D = 1  # Deletion
    I = 2  # Insertion
    N_MUTS = 3

def max3(vals: np.ndarray) -> Tuple[float, np.ndarray]:
    """
    找到三个值中的最大值，并返回哪些位置达到最大值
    
    Args:
        vals: 长度为3的数组 [S_score, D_score, I_score]
    
    Returns:
        maxval: 最大值
        argmax3: 长度为3的布尔数组，指示哪些位置达到最大值
    """
    maxval = np.max(vals)
    argmax3 = (vals == maxval)
    return maxval, argmax3

def cas9_align(seq: np.ndarray, ref: np.ndarray, 
               open_penalty: np.ndarray, close_penalty: np.ndarray,
               sub_score: np.ndarray) -> Tuple[float, List[int], List[int]]:
    """
    CRISPR-Cas9特异性序列比对算法
    
    Args:
        seq: 待比对序列，编码为整数 (A=1, C=2, G=3, T=4, gap=0)
        ref: 参考序列，编码格式同上
        open_penalty: 位置特异性的gap开启惩罚数组，大小为Lref+1
        close_penalty: 位置特异性的gap关闭惩罚数组，大小为Lref+1  
        sub_score: 25x1替换得分矩阵 (用于计算匹配/不匹配得分)
                  索引计算：score = sub_score[seq_nt * 5 + ref_nt]
    
    Returns:
        best_score: 最优比对得分
        al_seq: 比对后的序列
        al_ref: 比对后的参考序列
    """
    
    Lseq = len(seq)
    Lref = len(ref)
    
    # 初始化三个状态的得分矩阵 (S, D, I)
    score = np.full((Mutation.N_MUTS, Lseq + 1, Lref + 1), 
                    -np.inf, dtype=np.float64)
    
    # 初始化回溯矩阵 [from_state][to_state][i][j]
    backtrack = np.zeros((Mutation.N_MUTS, Mutation.N_MUTS, Lseq + 1, Lref + 1), 
                        dtype=bool)
    
    # 初始化边界条件
    score[Mutation.S, 0, 0] = 0.0
    
    # 初始化第一行 (插入序列)
    for j in range(1, Lseq + 1):
        score[Mutation.I, j, 0] = -open_penalty[0]
        backtrack[Mutation.I, Mutation.I, j, 0] = True
    
    # 初始化第一列 (删除序列) 
    for k in range(1, Lref + 1):
        score[Mutation.D, 0, k] = -open_penalty[0] 
        backtrack[Mutation.D, Mutation.D, 0, k] = True
    
    # 填充动态规划矩阵
    for j in range(1, Lseq + 1):
        for k in range(1, Lref + 1):
            # 计算替换得分 (注意：序列是1-indexed编码的)
            ss = sub_score[seq[j-1] * 5 + ref[k-1]]
            
            # 构建状态转移得分矩阵 [to_state][from_state]
            score_jk = np.array([
                # 转移到 S (substitution)
                [score[Mutation.S, j-1, k-1] + ss,
                 score[Mutation.D, j-1, k-1] - close_penalty[k-1] + ss,
                 score[Mutation.I, j-1, k-1] - close_penalty[k] + ss],
                
                # 转移到 D (deletion)  
                [score[Mutation.S, j, k-1] - open_penalty[k],
                 score[Mutation.D, j, k-1],
                 score[Mutation.I, j, k-1]],
                
                # 转移到 I (insertion)
                [score[Mutation.S, j-1, k] - open_penalty[k],
                 score[Mutation.D, j-1, k],
                 score[Mutation.I, j-1, k]]
            ])
            
            # 为每个状态找到最优前驱状态
            for m in range(Mutation.N_MUTS):
                score[m, j, k], out = max3(score_jk[m])
                backtrack[Mutation.S, m, j, k] = out[Mutation.S]
                backtrack[Mutation.D, m, j, k] = out[Mutation.D]  
                backtrack[Mutation.I, m, j, k] = out[Mutation.I]
    
    # 找到最终的最优状态
    score_LL = np.array([score[Mutation.S, Lseq, Lref],
                        score[Mutation.D, Lseq, Lref], 
                        score[Mutation.I, Lseq, Lref]])
    
    best_score, state_LL = max3(score_LL)
    
    # 确定最优结束状态
    if state_LL[Mutation.S]:
        cur_state = Mutation.S
    elif state_LL[Mutation.D]:
        cur_state = Mutation.D
    elif state_LL[Mutation.I]:
        cur_state = Mutation.I
    else:
        cur_state = -1
    
    # 回溯重构比对
    al_seq = []
    al_ref = []
    cur_j = Lseq
    cur_k = Lref
    
    while cur_j > 0 or cur_k > 0:
        assert 0 <= cur_state < Mutation.N_MUTS
        
        # 获取当前状态的回溯信息
        cur_bt = np.array([
            backtrack[Mutation.S, cur_state, cur_j, cur_k],
            backtrack[Mutation.D, cur_state, cur_j, cur_k], 
            backtrack[Mutation.I, cur_state, cur_j, cur_k]
        ])
        
        if cur_state == Mutation.S:
            # 替换：消耗序列和参考的一个字符
            al_seq.append(seq[cur_j - 1])
            al_ref.append(ref[cur_k - 1])
            cur_j -= 1
            cur_k -= 1
            # 选择前驱状态
            if cur_bt[Mutation.S]:
                cur_state = Mutation.S
            elif cur_bt[Mutation.D]:
                cur_state = Mutation.D
            elif cur_bt[Mutation.I]:
                cur_state = Mutation.I
            else:
                cur_state = -1
                
        elif cur_state == Mutation.D:
            # 删除：只消耗参考序列的一个字符
            al_seq.append(0)  # gap
            al_ref.append(ref[cur_k - 1])
            cur_k -= 1
            # 选择前驱状态 (优先级：D > S > I)
            if cur_bt[Mutation.D]:
                cur_state = Mutation.D
            elif cur_bt[Mutation.S]:
                cur_state = Mutation.S
            elif cur_bt[Mutation.I]:
                cur_state = Mutation.I
            else:
                cur_state = -1
                
        elif cur_state == Mutation.I:
            # 插入：只消耗序列的一个字符
            al_seq.append(seq[cur_j - 1])
            al_ref.append(0)  # gap
            cur_j -= 1
            # 选择前驱状态 (优先级：I > S > D)
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
    
    # 反转得到正确顺序的比对
    al_seq.reverse()
    al_ref.reverse()
    
    return best_score, al_seq, al_ref

def nt2int(sequence: str) -> np.ndarray:
    """将DNA序列转换为整数编码"""
    mapping = {'A': 1, 'C': 2, 'G': 3, 'T': 4}
    return np.array([mapping.get(nt.upper(), 0) for nt in sequence], dtype=np.uint8)

def int2nt(int_seq: List[int]) -> str:
    """将整数编码转换回DNA序列"""
    mapping = {0: '-', 1: 'A', 2: 'C', 3: 'G', 4: 'T'}
    return ''.join(mapping.get(nt, 'N') for nt in int_seq)

def print_cas9_alignment(al_seq: List[int], al_ref: List[int], score: float):
    """打印比对结果"""
    seq_str = int2nt(al_seq)
    ref_str = int2nt(al_ref)
    
    print(f"序列: {seq_str}")
    print(f"参考: {ref_str}")
    print(f"得分: {score:.2f}")
    
    # 打印匹配情况
    match_line = ""
    for i in range(len(al_seq)):
        if al_seq[i] == al_ref[i] and al_seq[i] != 0:
            match_line += "|"
        elif al_seq[i] == 0 or al_ref[i] == 0:
            match_line += " "
        else:
            match_line += "."
    print(f"匹配: {match_line}")

# 示例用法
if __name__ == "__main__":
    # 示例序列
    seq_str = "ACGTACGT"
    ref_str = "ACGTGCGT"
    
    # 转换为整数编码
    seq = nt2int(seq_str)
    ref = nt2int(ref_str)
    
    # 设置简单的惩罚参数 (实际使用中应该根据CRISPR切割位点调整)
    open_penalty = np.full(len(ref_str) + 1, 2.0)  # 统一的开启惩罚
    close_penalty = np.full(len(ref_str) + 1, 1.0)     # 统一的关闭惩罚
    
    # 设置简单的NUC44评分矩阵 (实际应该使用完整的25x25矩阵)
    sub_score = np.zeros(25, dtype=np.float64)
    # 简化的匹配/不匹配得分
    for i in range(1, 5):
        for j in range(1, 5):
            if i == j:
                sub_score[i * 5 + j] = 2.0  # 匹配
            else:
                sub_score[i * 5 + j] = -1.0  # 不匹配
    
    # 执行比对
    score, al_seq, al_ref = cas9_align(seq, ref, open_penalty, close_penalty, sub_score)
    
    # 打印结果
    print("=== 示例1: 简单序列比对 ===")
    print_cas9_alignment(al_seq, al_ref, score)
    
    # 示例2：包含插入删除的序列
    print("\n=== 示例2: 包含插入删除的序列 ===")
    seq_str2 = "ACGTGCGT"
    ref_str2 = "ACGTACCCGT"
    
    seq2 = nt2int(seq_str2)
    ref2 = nt2int(ref_str2)
    
    # 调整惩罚参数大小
    open_penalty2 = np.full(len(ref_str2) + 1, 2.0)
    close_penalty2 = np.full(len(ref_str2) + 1, 1.0)
    
    # 在位置5附近设置较低的gap惩罚，模拟CRISPR切割位点
    open_penalty2[4:7] = 0.5
    close_penalty2[4:7] = 0.2
    
    score2, al_seq2, al_ref2 = cas9_align(seq2, ref2, open_penalty2, close_penalty2, sub_score)
    print_cas9_alignment(al_seq2, al_ref2, score2)
    
    print(f"\n注意：在位置4-6设置了较低的gap惩罚({open_penalty2[5]:.1f})，")
    print(f"模拟CRISPR切割位点，其他位置gap惩罚为{open_penalty2[0]:.1f}")