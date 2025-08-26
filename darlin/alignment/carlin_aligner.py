#!/usr/bin/env python3
"""
CARLIN特异性序列比对器

集成CARLIN配置和cas9_align算法，提供高级比对接口
"""

import numpy as np
from typing import Tuple, Optional, Dict, List

from .cas9_align import cas9_align, nt2int, int2nt, print_cas9_alignment
from .aligned_seq import AlignedSEQ, AlignedSEQMotif, SequenceSanitizer, desemble_sequence, calculate_motif_boundaries
from ..config import AmpliconConfig, get_original_carlin_config, ScoringConfig, get_default_scoring_config


class CARLINAligner:
    """
    CARLIN特异性序列比对器
    
    提供针对CARLIN扩增子的高级比对功能，自动处理配置和参数
    """
    
    def __init__(self, 
                 amplicon_config: Optional[AmpliconConfig] = None,
                 scoring_config: Optional[ScoringConfig] = None):
        """
        初始化比对器
        
        Args:
            amplicon_config: CARLIN扩增子配置，默认使用OriginalCARLIN
            scoring_config: 评分配置，默认使用NUC44
        """
        # 使用默认配置
        self.amplicon_config = amplicon_config or get_original_carlin_config()
        self.scoring_config = scoring_config or get_default_scoring_config()
        
        # 获取预计算的参数
        self.reference_sequence = self.amplicon_config.get_reference_sequence()
        self.open_penalty_array, self.close_penalty_array = self.amplicon_config.get_penalty_arrays()
        self.substitution_matrix = self.scoring_config.substitution_matrix
        
        # 编码参考序列
        self.reference_encoded = nt2int(self.reference_sequence)
        
        print(f"✅ CARLIN比对器初始化成功")
        print(f"   - 参考序列长度: {len(self.reference_sequence)} bp")
        print(f"   - 评分矩阵: {self.scoring_config.matrix_type}")
        print(f"   - Gap惩罚范围: {self.open_penalty_array.min():.1f}-{self.open_penalty_array.max():.1f}")
    
    def align_sequence(self, query_sequence: str, 
                      verbose: bool = False,
                      sanitize: bool = True) -> Dict:
        """
        比对单个序列到CARLIN参考序列
        
        Args:
            query_sequence: 查询序列
            verbose: 是否显示详细信息
            sanitize: 是否执行序列标准化
            
        Returns:
            Dict: 比对结果字典
        """
        # 输入验证
        if not query_sequence or not all(c in 'ACGTN-' for c in query_sequence.upper()):
            raise ValueError("序列只能包含ACGTN-字符")
        
        # 编码查询序列
        query_encoded = nt2int(query_sequence.upper().replace('N', 'A'))  # 将N替换为A
        
        # 执行比对
        score, aligned_query, aligned_ref = cas9_align(
            query_encoded, 
            self.reference_encoded,
            self.open_penalty_array,
            self.close_penalty_array,
            self.substitution_matrix
        )
        
        # 解码比对结果
        aligned_query_str = int2nt(aligned_query)
        aligned_ref_str = int2nt(aligned_ref)
        
        # 序列标准化（可选）
        sanitized_aligned_seq = None
        if sanitize:
            try:
                # 将比对结果分解为motifs
                motif_boundaries = calculate_motif_boundaries(aligned_ref_str, self.amplicon_config)
                aligned_seq_obj = desemble_sequence(aligned_query_str, aligned_ref_str, motif_boundaries)
                
                # 执行标准化
                sanitized_seq = self._perform_sanitization(aligned_seq_obj)
                
                # 更新比对结果为标准化后的序列
                sanitized_aligned_query_str = sanitized_seq.get_seq()
                sanitized_aligned_ref_str = sanitized_seq.get_ref()
                sanitized_aligned_seq = sanitized_seq
                
                # 使用标准化后的序列进行后续分析
                aligned_query_str = sanitized_aligned_query_str
                aligned_ref_str = sanitized_aligned_ref_str
                
            except Exception as e:
                if verbose:
                    print(f"⚠️  序列标准化失败: {e}")
                    print("   继续使用原始比对结果")
        
        # 重新编码用于统计计算
        if sanitize and sanitized_aligned_seq:
            # 使用标准化后的序列计算统计
            aligned_query_for_stats = nt2int(aligned_query_str)
            aligned_ref_for_stats = nt2int(aligned_ref_str) 
        else:
            aligned_query_for_stats = aligned_query
            aligned_ref_for_stats = aligned_ref
        
        # 计算比对统计
        stats = self._calculate_alignment_stats(aligned_query_for_stats, aligned_ref_for_stats)
        
        # 构建结果
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
        批量比对多个序列
        
        Args:
            sequences: 序列列表
            verbose: 是否显示详细信息
            
        Returns:
            List[Dict]: 比对结果列表
        """
        results = []
        
        for i, seq in enumerate(sequences):
            if verbose:
                print(f"\n=== 比对序列 {i+1}/{len(sequences)} ===")
            
            try:
                result = self.align_sequence(seq, verbose=verbose)
                results.append(result)
            except Exception as e:
                print(f"⚠️  序列 {i+1} 比对失败: {e}")
                results.append({
                    'query_sequence': seq,
                    'error': str(e),
                    'alignment_score': float('-inf')
                })
        
        return results
    
    def _calculate_alignment_stats(self, aligned_query: np.ndarray, 
                                  aligned_ref: np.ndarray) -> Dict:
        """计算比对统计信息"""
        matches = 0
        mismatches = 0
        query_gaps = 0
        ref_gaps = 0
        
        for i in range(len(aligned_query)):
            q_base = aligned_query[i]
            r_base = aligned_ref[i]
            
            if q_base == 0:  # 查询序列gap
                query_gaps += 1
            elif r_base == 0:  # 参考序列gap
                ref_gaps += 1
            elif q_base == r_base:  # 匹配
                matches += 1
            else:  # 不匹配
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
        """分析各motif的比对情况"""
        motif_stats = {
            'prefix': {'matches': 0, 'mismatches': 0, 'gaps': 0},
            'segments': [{} for _ in range(10)],
            'pams': [{} for _ in range(9)],
            'postfix': {'matches': 0, 'mismatches': 0, 'gaps': 0}
        }
        
        # 分析每个位置的motif信息
        ref_pos = 0  # 参考序列的实际位置
        
        for i in range(len(aligned_query)):
            q_base = aligned_query[i]
            r_base = aligned_ref[i]
            
            # 只有当参考序列不是gap时才分析motif
            if r_base != 0:
                motif_info = self.amplicon_config.get_motif_info(ref_pos)
                motif_type = motif_info['type']
                motif_id = motif_info['motif_id']
                
                # 统计匹配/不匹配/gap
                if q_base == 0:
                    event = 'gaps'
                elif q_base == r_base:
                    event = 'matches'
                else:
                    event = 'mismatches'
                
                # 更新相应motif的统计
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
        """打印比对结果"""
        print(f"查询序列长度: {len(result['query_sequence'])} bp")
        print(f"比对得分: {result['alignment_score']:.2f}")
        
        stats = result['statistics']
        print(f"比对统计:")
        print(f"  - 比对长度: {stats['aligned_length']}")
        print(f"  - 匹配: {stats['matches']} ({stats['identity']*100:.1f}%)")
        print(f"  - 不匹配: {stats['mismatches']}")
        print(f"  - 查询gaps: {stats['query_gaps']}")
        print(f"  - 参考gaps: {stats['reference_gaps']}")
        
        # 显示比对结果
        print("\n比对结果:")
        aligned_query = result['aligned_query']
        aligned_ref = result['aligned_reference']
        
        # 分段显示，每行60个字符
        line_length = 60
        for start in range(0, len(aligned_query), line_length):
            end = min(start + line_length, len(aligned_query))
            
            query_line = aligned_query[start:end]
            ref_line = aligned_ref[start:end]
            
            # 构建匹配标记行
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
        """获取配置摘要"""
        lines = [
            "=== CARLIN比对器配置 ===",
            f"参考序列: {self.reference_sequence[:50]}...",
            f"序列长度: {len(self.reference_sequence)} bp",
            f"评分矩阵: {self.scoring_config.matrix_type}",
            f"匹配得分: {self.scoring_config.get_score(1, 1)}",
            f"不匹配得分: {self.scoring_config.get_score(1, 2)}",
            f"Gap惩罚范围: {self.open_penalty_array.min():.1f} - {self.open_penalty_array.max():.1f}",
            f"Segments数量: {len(self.amplicon_config.sequence.segments)}",
            f"PAM序列: {self.amplicon_config.sequence.pam}"
        ]
        return "\n".join(lines)
    
    def _perform_sanitization(self, aligned_seq: AlignedSEQ) -> AlignedSEQ:
        """
        执行序列标准化
        
        Args:
            aligned_seq: 分解后的比对序列对象
            
        Returns:
            标准化后的AlignedSEQ对象
        """
        # 步骤1: Prefix/Postfix标准化
        sanitized_seq = SequenceSanitizer.sanitize_prefix_postfix(aligned_seq)
        
        # 步骤2: 保守区域标准化
        # 需要确定cutsite motif的索引
        cutsite_indices = self._get_cutsite_motif_indices()
        sanitized_seq = SequenceSanitizer.sanitize_conserved_regions(sanitized_seq, cutsite_indices)
        
        return sanitized_seq
    
    def _get_cutsite_motif_indices(self) -> List[int]:
        """
        获取cutsite motif在motif列表中的索引
        
        根据CARLIN配置结构，motif顺序为：
        prefix, consite1, cutsite1, pam1, consite2, cutsite2, pam2, ..., consite10, cutsite10, postfix
        
        Returns:
            cutsite motif的索引列表
        """
        cutsite_indices = []
        motif_index = 0
        
        # Prefix (index 0)
        motif_index += 1
        
        # 10个segments，每个包含consite和cutsite
        for i in range(10):
            # Consite
            motif_index += 1
            
            # Cutsite (这是我们要标记的)
            cutsite_indices.append(motif_index)
            motif_index += 1
            
            # PAM (前9个segment后面有PAM)
            if i < 9:
                motif_index += 1
        
        # Postfix已经在最后
        
        return cutsite_indices


def create_default_aligner() -> CARLINAligner:
    """创建默认的CARLIN比对器"""
    return CARLINAligner()


def align_to_carlin(sequence: str, verbose: bool = True) -> Dict:
    """
    便捷函数：将单个序列比对到CARLIN参考序列
    
    Args:
        sequence: 查询序列
        verbose: 是否显示详细信息
        
    Returns:
        Dict: 比对结果
    """
    aligner = create_default_aligner()
    return aligner.align_sequence(sequence, verbose=verbose) 