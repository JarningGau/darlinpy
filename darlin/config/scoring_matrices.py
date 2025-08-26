#!/usr/bin/env python3
"""
评分矩阵模块

实现NUC44和其他核酸序列比对评分矩阵
"""

import numpy as np
from typing import Dict, Tuple


def create_nuc44_matrix() -> np.ndarray:
    """
    创建NUC44核酸替换评分矩阵
    
    这是NCBI BLAST使用的标准核酸评分矩阵
    索引顺序: 0=gap, 1=A, 2=C, 3=G, 4=T
    
    Returns:
        np.ndarray: 5x5的评分矩阵，打平为25元素的一维数组
    """
    # NUC44矩阵 (5x5)
    # 行/列顺序: gap(0), A(1), C(2), G(3), T(4)
    nuc44_2d = np.array([
        #    gap   A    C    G    T
        [   0,   0,   0,   0,   0],  # gap
        [   0,   5,  -4,  -4,  -4],  # A
        [   0,  -4,   5,  -4,  -4],  # C  
        [   0,  -4,  -4,   5,  -4],  # G
        [   0,  -4,  -4,  -4,   5],  # T
    ], dtype=np.float64)
    
    # 打平为一维数组，与cas9_align期望的格式一致
    return nuc44_2d.flatten()


def create_simple_scoring_matrix(match_score: float = 5.0, mismatch_score: float = -4.0) -> np.ndarray:
    """
    创建简单的匹配/不匹配评分矩阵
    
    Args:
        match_score: 匹配得分
        mismatch_score: 不匹配得分
        
    Returns:
        np.ndarray: 打平的5x5评分矩阵
    """
    matrix = np.zeros(25, dtype=np.float64)
    
    # A=1, C=2, G=3, T=4, gap=0
    for i in range(1, 5):
        for j in range(1, 5):
            idx = i * 5 + j
            if i == j:
                matrix[idx] = match_score    # 匹配
            else:
                matrix[idx] = mismatch_score # 不匹配
    
    return matrix


def create_transition_transversion_matrix(
    match_score: float = 5.0,
    transition_score: float = -1.0,  # A<->G, C<->T
    transversion_score: float = -4.0  # 其他替换
) -> np.ndarray:
    """
    创建考虑转换/颠换的评分矩阵
    
    转换 (transition): A<->G, C<->T (通常更常见)
    颠换 (transversion): A<->C, A<->T, G<->C, G<->T
    
    Args:
        match_score: 完全匹配得分
        transition_score: 转换得分 (惩罚较轻)
        transversion_score: 颠换得分 (惩罚较重)
        
    Returns:
        np.ndarray: 打平的5x5评分矩阵
    """
    matrix = np.zeros(25, dtype=np.float64)
    
    # 定义转换对
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
    """评分配置类"""
    
    def __init__(self, matrix_type: str = "nuc44", **kwargs):
        """
        初始化评分配置
        
        Args:
            matrix_type: 矩阵类型 ("nuc44", "simple", "transition_transversion")
            **kwargs: 传递给矩阵创建函数的参数
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
            raise ValueError(f"未知的矩阵类型: {matrix_type}")
    
    def get_score(self, base1: int, base2: int) -> float:
        """
        获取两个碱基之间的得分
        
        Args:
            base1, base2: 碱基编码 (0=gap, 1=A, 2=C, 3=G, 4=T)
            
        Returns:
            float: 替换得分
        """
        return self.substitution_matrix[base1 * 5 + base2]
    
    def summary(self) -> str:
        """返回评分配置摘要"""
        lines = [
            f"=== 评分矩阵配置 ({self.matrix_type}) ===",
            "",
            "评分矩阵 (行:查询序列, 列:参考序列):",
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


# 预定义的评分配置
DEFAULT_NUC44 = ScoringConfig("nuc44")
DEFAULT_SIMPLE = ScoringConfig("simple", match_score=5.0, mismatch_score=-4.0)


def get_default_scoring_config() -> ScoringConfig:
    """获取默认的评分配置 (NUC44)"""
    return DEFAULT_NUC44 