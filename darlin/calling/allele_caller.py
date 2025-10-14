"""
CARLIN等位基因调用算法

这个模块实现了CARLIN的等位基因调用核心算法，包括精确调用和粗粒度调用两种方法。
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
    CARLIN等位基因调用器
    
    实现了两种调用策略：
    1. call_alleles_exact: 基于完整序列的精确调用
    2. call_alleles_coarse_grain: 基于事件结构的粗粒度调用
    """
    
    def __init__(self, 
             locus: str = "Col1a1",
             amplicon_config: Optional[AmpliconConfig] = None, 
             dominant_threshold: float = 0.5):
        """
        初始化等位基因调用器
        
        Args:
            locus: 位点名称，用于选择对应的JSON配置模板，默认为"Col1a1"
            amplicon_config: CARLIN扩增子配置，如果提供则优先使用
            dominant_threshold: 主导等位基因的最小比例阈值
        """
        # 根据locus加载配置
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
        按频率排序的唯一值获取
        
        Args:
            values: 输入值列表
            weights: 每个值的权重，默认为1
            
        Returns:
            (unique_values, backtrack, unique_indices)
            - unique_values: 按频率降序排列的唯一值
            - backtrack: 原始索引到唯一值索引的映射
            - unique_indices: 每个原始值对应的唯一值索引
        """
        if weights is None:
            weights = [1.0] * len(values)
            
        if len(values) != len(weights):
            raise ValueError("values和weights长度必须一致")
            
        # 创建值到索引的映射
        value_to_indices = {}
        for i, value in enumerate(values):
            # 对于AlignedSEQ对象，使用序列作为键
            if isinstance(value, AlignedSEQ):
                key = value.get_seq()
            elif isinstance(value, str):
                key = value
            else:
                key = str(value)
                
            if key not in value_to_indices:
                value_to_indices[key] = []
            value_to_indices[key].append(i)
            
        # 计算每个唯一值的权重
        unique_keys = []
        unique_weights = []
        unique_values = []
        
        for key, indices in value_to_indices.items():
            total_weight = sum(weights[i] for i in indices)
            unique_keys.append(key)
            unique_weights.append(total_weight)
            unique_values.append(values[indices[0]])  # 使用第一个遇到的值作为代表
            
        # 按权重降序排序
        sorted_indices = sorted(range(len(unique_weights)), 
                              key=lambda i: unique_weights[i], 
                              reverse=True)
        
        sorted_unique_values = [unique_values[i] for i in sorted_indices]
        
        # 创建映射
        key_to_sorted_index = {}
        for new_idx, old_idx in enumerate(sorted_indices):
            key_to_sorted_index[unique_keys[old_idx]] = new_idx
            
        # 创建backtrack和unique_indices
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
        基于完整序列的精确等位基因调用
        
        Args:
            aligned_seqs: 已对齐的序列列表
            seq_weights: 每个序列的权重
            dominant_only: 是否只返回主导等位基因
            
        Returns:
            等位基因调用结果
        """
        if seq_weights is None:
            seq_weights = [1.0] * len(aligned_seqs)
            
        # 过滤空序列
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
        
        # 提取序列字符串
        seq_strings = [seq.get_seq() for seq in valid_seqs]
        
        # 按频率获取唯一序列
        unique_seqs, _, which_seq = self.unique_by_frequency(seq_strings, valid_weights)
        
        # 计算每个唯一序列的权重
        seq_weight_sums = []
        for i in range(len(unique_seqs)):
            weight_sum = sum(valid_weights[j] for j, seq_idx in enumerate(which_seq) if seq_idx == i)
            seq_weight_sums.append(weight_sum)
            
        total_weight = sum(seq_weight_sums)
        
        # 检查主导性
        dominant_fraction = seq_weight_sums[0] / total_weight if total_weight > 0 else 0
        
        if dominant_only and dominant_fraction < self.dominant_threshold:
            return AlleleCallResult(
                allele=None,
                constituents=[],
                weight_contribution=[],
                calling_method='exact',
                confidence=dominant_fraction
            )
            
        # 找到最高权重序列的代表
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
        基于事件结构的粗粒度等位基因调用
        
        Args:
            aligned_seqs: 已对齐的序列列表
            seq_weights: 每个序列的权重
            dominant_only: 是否只返回主导等位基因
            
        Returns:
            等位基因调用结果
        """
        if seq_weights is None:
            seq_weights = [1.0] * len(aligned_seqs)
            
        # 过滤空序列
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
        
        # 提取事件结构
        event_structures = [seq.get_event_structure() for seq in valid_seqs]
        event_strings = [''.join(events) for events in event_structures]
        
        # 按频率获取唯一事件结构
        unique_events, _, which_event = self.unique_by_frequency(event_strings, valid_weights)
        
        # 计算每个唯一事件的权重
        event_weight_sums = []
        for i in range(len(unique_events)):
            weight_sum = sum(valid_weights[j] for j, event_idx in enumerate(which_event) if event_idx == i)
            event_weight_sums.append(weight_sum)
            
        total_weight = sum(event_weight_sums)
        
        # 检查主导性
        dominant_fraction = event_weight_sums[0] / total_weight if total_weight > 0 else 0
        
        if dominant_only and dominant_fraction < self.dominant_threshold:
            return AlleleCallResult(
                allele=None,
                constituents=[],
                weight_contribution=[],
                calling_method='coarse_grain',
                confidence=dominant_fraction
            )
            
        # 收集具有主导事件结构的序列
        target_event_indices = [i for i, event_idx in enumerate(which_event) if event_idx == 0]
        event_seqs = [valid_seqs[i] for i in target_event_indices]
        event_weights = [valid_weights[i] for i in target_event_indices]
        
        # 构建共识序列
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
        从多个序列构建共识序列
        
        Args:
            seqs: 序列列表
            weights: 权重列表
            
        Returns:
            共识AlignedSEQ对象
        """
        if not seqs:
            raise ValueError("序列列表不能为空")
            
        # 提取序列字符串
        seq_strings = [seq.get_seq() for seq in seqs]
        ref_strings = [seq.get_ref() for seq in seqs]
        
        # 检查所有序列长度是否相同
        seq_lengths = [len(s) for s in seq_strings]
        if len(set(seq_lengths)) == 1:
            # 长度相同，使用模式统计
            consensus_seq = self._compute_mode_sequence(seq_strings, weights)
            # 使用第一个序列的参考作为参考
            consensus_ref = ref_strings[0]
        else:
            # 长度不同，使用最常见长度的序列
            length_weights = {}
            for i, length in enumerate(seq_lengths):
                if length not in length_weights:
                    length_weights[length] = 0
                length_weights[length] += weights[i]
                
            mode_length = max(length_weights.keys(), key=lambda x: length_weights[x])
            
            # 过滤到模式长度的序列
            mode_indices = [i for i, length in enumerate(seq_lengths) if length == mode_length]
            mode_seqs = [seq_strings[i] for i in mode_indices]
            mode_weights = [weights[i] for i in mode_indices]
            
            consensus_seq = self._compute_mode_sequence(mode_seqs, mode_weights)
            consensus_ref = ref_strings[mode_indices[0]]
            
        # 重新构建AlignedSEQ对象
        try:
            # 使用desemble_sequence重新构建
            from ..alignment.aligned_seq import calculate_motif_boundaries
            motif_boundaries = calculate_motif_boundaries(consensus_ref, self.amplicon_config)
            return desemble_sequence(consensus_seq, consensus_ref, motif_boundaries)
        except Exception as e:
            # 如果失败，创建简单的单motif AlignedSEQ
            return AlignedSEQ([consensus_seq], [consensus_ref])
            
    def _compute_mode_sequence(self, 
                              sequences: List[str], 
                              weights: List[float]) -> str:
        """
        计算加权模式序列
        
        Args:
            sequences: 序列字符串列表
            weights: 权重列表
            
        Returns:
            共识序列字符串
        """
        if not sequences:
            return ""
            
        # 确保所有序列长度相同
        seq_length = len(sequences[0])
        if not all(len(seq) == seq_length for seq in sequences):
            raise ValueError("所有序列必须具有相同长度")
            
        consensus = []
        
        for pos in range(seq_length):
            # 收集该位置的字符和权重
            char_weights = {}
            for seq, weight in zip(sequences, weights):
                char = seq[pos]
                if char not in char_weights:
                    char_weights[char] = 0
                char_weights[char] += weight
                
            # 选择权重最高的字符
            mode_char = max(char_weights.keys(), key=lambda x: char_weights[x])
            consensus.append(mode_char)
            
        return ''.join(consensus)
        
    def call_bulk_alleles(self, 
                         sequence_groups: List[List[AlignedSEQ]],
                         group_weights: Optional[List[List[float]]] = None,
                         method: str = 'coarse_grain',
                         dominant_only: bool = True) -> BulkAlleleCallResult:
        """
        批量等位基因调用
        
        Args:
            sequence_groups: 序列组列表，每组代表一个UMI或细胞
            group_weights: 每组中序列的权重
            method: 调用方法 ('exact' 或 'coarse_grain')
            dominant_only: 是否只返回主导等位基因
            
        Returns:
            批量等位基因调用结果
        """
        if group_weights is None:
            group_weights = [[1.0] * len(group) for group in sequence_groups]
            
        if len(sequence_groups) != len(group_weights):
            raise ValueError("sequence_groups和group_weights长度必须一致")
            
        # 调用每个组的等位基因
        individual_results = []
        calling_method = getattr(self, f'call_alleles_{method}')
        
        for group_seqs, group_w in zip(sequence_groups, group_weights):
            result = calling_method(group_seqs, group_w, dominant_only)
            individual_results.append(result)
            
        # 汇总所有成功调用的等位基因
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
            
        # 按频率汇总等位基因
        allele_seqs = [allele.get_seq() for allele in successful_alleles]
        unique_alleles, _, _ = self.unique_by_frequency(allele_seqs, successful_weights)
        
        # 计算频率
        allele_counts = Counter()
        total_weight = sum(successful_weights)
        
        for allele_seq, weight in zip(allele_seqs, successful_weights):
            allele_counts[allele_seq] += weight
            
        # 排序并计算频率
        sorted_alleles = sorted(allele_counts.items(), key=lambda x: x[1], reverse=True)
        summary_alleles = []
        allele_frequencies = []
        
        for allele_seq, count in sorted_alleles:
            # 找到对应的AlignedSEQ对象
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