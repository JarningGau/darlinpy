#!/usr/bin/env python3
"""
序列处理工具函数

提供CARLIN分析所需的基础DNA序列处理功能
"""

from typing import List, Tuple, Dict
from collections import Counter


# DNA序列相关常量
DNA_BASES = set('ATCG')
DNA_BASES_AMBIGUOUS = set('ATCGRYSWKMBDHVN')
COMPLEMENT_MAP = {
    'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
    'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
    'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
    'H': 'D', 'V': 'B', 'N': 'N'
}


def reverse_complement(seq: str) -> str:
    """
    计算DNA序列的反向互补序列
    
    Args:
        seq: DNA序列字符串
        
    Returns:
        str: 反向互补序列
        
    Examples:
        >>> reverse_complement("ATCG")
        'CGAT'
    """
    if not seq:
        return ""
    
    # 转换为大写并验证
    seq = seq.upper()
    if not is_valid_dna(seq, allow_ambiguous=True):
        raise ValueError(f"无效的DNA序列: {seq}")
    
    # 生成反向互补序列
    complement = ''.join(COMPLEMENT_MAP.get(base, 'N') for base in seq)
    return complement[::-1]


def is_valid_dna(seq: str, allow_ambiguous: bool = False) -> bool:
    """
    验证是否为有效的DNA序列
    
    Args:
        seq: 待验证的序列
        allow_ambiguous: 是否允许模糊碱基
        
    Returns:
        bool: 是否为有效DNA序列
    """
    if not seq:
        return True
    
    seq = seq.upper()
    allowed_bases = DNA_BASES_AMBIGUOUS if allow_ambiguous else DNA_BASES
    
    return all(base in allowed_bases for base in seq)


def clean_sequence(seq: str, remove_chars: str = "- \n\t") -> str:
    """
    清理序列，移除指定字符
    
    Args:
        seq: 原始序列
        remove_chars: 要移除的字符
        
    Returns:
        str: 清理后的序列
    """
    if not seq:
        return ""
    
    cleaned = seq.upper()
    for char in remove_chars:
        cleaned = cleaned.replace(char, '')
    
    return cleaned


def calculate_gc_content(seq: str) -> float:
    """
    计算GC含量
    
    Args:
        seq: DNA序列
        
    Returns:
        float: GC含量百分比 (0-100)
    """
    if not seq:
        return 0.0
    
    seq = clean_sequence(seq)
    if not is_valid_dna(seq):
        raise ValueError(f"无效的DNA序列: {seq}")
    
    gc_count = seq.count('G') + seq.count('C')
    return (gc_count / len(seq)) * 100 if len(seq) > 0 else 0.0


def count_nucleotides(seq: str) -> Dict[str, int]:
    """
    统计核苷酸数量
    
    Args:
        seq: DNA序列
        
    Returns:
        Dict[str, int]: 各核苷酸的数量
    """
    if not seq:
        return {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    
    seq = clean_sequence(seq)
    if not is_valid_dna(seq):
        raise ValueError(f"无效的DNA序列: {seq}")
    
    counts = Counter(seq)
    return {
        'A': counts.get('A', 0),
        'T': counts.get('T', 0),
        'G': counts.get('G', 0),
        'C': counts.get('C', 0)
    }


def hamming_distance(seq1: str, seq2: str) -> int:
    """
    计算两个等长序列的汉明距离
    
    Args:
        seq1: 第一个序列
        seq2: 第二个序列
        
    Returns:
        int: 汉明距离
    """
    if len(seq1) != len(seq2):
        raise ValueError("序列长度必须相等")
    
    return sum(c1 != c2 for c1, c2 in zip(seq1.upper(), seq2.upper()))


def format_sequence(seq: str, line_length: int = 80, add_numbers: bool = False) -> str:
    """
    格式化序列输出
    
    Args:
        seq: DNA序列
        line_length: 每行长度
        add_numbers: 是否添加行号
        
    Returns:
        str: 格式化的序列
    """
    if not seq:
        return ""
    
    seq = clean_sequence(seq)
    lines = []
    
    for i in range(0, len(seq), line_length):
        line = seq[i:i+line_length]
        if add_numbers:
            line_num = i + 1
            line = f"{line_num:>8} {line}"
        lines.append(line)
    
    return '\n'.join(lines)


# 保留简单的FASTA处理功能，因为可能需要处理输入文件
def parse_fasta_string(fasta_string: str) -> List[Tuple[str, str]]:
    """
    解析FASTA格式字符串
    
    Args:
        fasta_string: FASTA格式的字符串
        
    Returns:
        List[Tuple[str, str]]: (序列名, 序列)的列表
    """
    sequences = []
    lines = fasta_string.strip().split('\n')
    
    current_name = None
    current_seq = []
    
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            # 保存前一个序列
            if current_name is not None:
                sequences.append((current_name, ''.join(current_seq)))
            
            # 开始新序列
            current_name = line[1:]  # 移除'>'
            current_seq = []
        elif line and current_name is not None:
            current_seq.append(line)
    
    # 保存最后一个序列
    if current_name is not None:
        sequences.append((current_name, ''.join(current_seq)))
    
    return sequences


def to_fasta_string(sequences: List[Tuple[str, str]], line_length: int = 80) -> str:
    """
    将序列列表转换为FASTA格式字符串
    
    Args:
        sequences: (序列名, 序列)的列表
        line_length: 每行序列长度
        
    Returns:
        str: FASTA格式字符串
    """
    fasta_lines = []
    
    for name, seq in sequences:
        fasta_lines.append(f">{name}")
        
        # 分行输出序列
        for i in range(0, len(seq), line_length):
            fasta_lines.append(seq[i:i+line_length])
    
    return '\n'.join(fasta_lines) 