#!/usr/bin/env python3
"""
CARLIN比对序列数据结构和标准化功能

实现AlignedSEQ、AlignedSEQMotif类以及prefix/postfix和保守区域的标准化算法
"""

from typing import List, Tuple, Optional, Union
import numpy as np


class AlignedSEQMotif:
    """
    比对的序列motif片段
    
    表示比对序列中的一个motif片段，包含序列、参考序列和事件类型
    """
    
    def __init__(self, seq: str, ref: str):
        """
        初始化motif
        
        Args:
            seq: 查询序列片段
            ref: 参考序列片段
        """
        if len(seq) != len(ref):
            raise ValueError(f"序列长度({len(seq)})与参考序列长度({len(ref)})不匹配")
        
        self.seq = seq
        self.ref = ref
        self.event = self._classify_motif_event(seq, ref)
    
    def _classify_motif_event(self, seq: str, ref: str) -> str:
        """
        分类motif事件类型
        
        Args:
            seq: 查询序列
            ref: 参考序列
            
        Returns:
            str: 事件类型
                'N' - No change (无变化)
                'M' - Mismatch (错配)
                'D' - Deletion (删除)
                'I' - Insertion (插入)
                'E' - Empty (空)
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
    完整的比对序列对象
    
    由多个AlignedSEQMotif组成，表示完整的比对序列
    """
    
    def __init__(self, seq_segments: List[str], ref_segments: List[str]):
        """
        初始化比对序列
        
        Args:
            seq_segments: 查询序列的motif片段列表
            ref_segments: 参考序列的motif片段列表
        """
        if len(seq_segments) != len(ref_segments):
            raise ValueError(f"序列segments数量({len(seq_segments)})与参考segments数量({len(ref_segments)})不匹配")
        
        self.motifs = [AlignedSEQMotif(seq, ref) for seq, ref in zip(seq_segments, ref_segments)]
    
    def get_seq(self) -> str:
        """获取完整的查询序列"""
        return ''.join(motif.seq for motif in self.motifs)
    
    def get_ref(self) -> str:
        """获取完整的参考序列"""
        return ''.join(motif.ref for motif in self.motifs)
    
    def get_event_structure(self) -> List[str]:
        """获取事件结构"""
        return [motif.event for motif in self.motifs]
    
    def copy(self) -> 'AlignedSEQ':
        """创建深拷贝"""
        seq_segments = [motif.seq for motif in self.motifs]
        ref_segments = [motif.ref for motif in self.motifs]
        return AlignedSEQ(seq_segments, ref_segments)
    
    def __repr__(self):
        events = ''.join(self.get_event_structure())
        return f"AlignedSEQ(motifs={len(self.motifs)}, events='{events}')"


class SequenceSanitizer:
    """
    序列标准化器
    
    实现CARLIN序列的prefix/postfix和保守区域标准化
    """
    
    @staticmethod
    def sanitize_prefix_postfix(aligned_seqs: Union[AlignedSEQ, List[AlignedSEQ]]) -> Union[AlignedSEQ, List[AlignedSEQ]]:
        """
        清理前缀和后缀的非功能性插入
        
        移除在CARLIN前缀和后缀区域的插入，这些通常是PCR人工制品或测序错误
        
        Args:
            aligned_seqs: 单个AlignedSEQ或AlignedSEQ列表
            
        Returns:
            标准化后的AlignedSEQ或列表
        """
        # 处理单个对象的情况
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
            
            # 处理头部插入 (prefix)
            if event_structure and event_structure[0] == 'I':
                first_motif = motifs[0]
                ref_seq = first_motif.ref
                
                # 找到参考序列中第一个非gap字符的位置
                first_non_gap = None
                for i, char in enumerate(ref_seq):
                    if char != '-':
                        first_non_gap = i
                        break
                
                # 如果开头有gap，则修剪
                if first_non_gap is not None and first_non_gap > 0:
                    trimmed_seq = first_motif.seq[first_non_gap:]
                    trimmed_ref = first_motif.ref[first_non_gap:]
                    motifs[0] = AlignedSEQMotif(trimmed_seq, trimmed_ref)
            
            # 处理尾部插入 (postfix)
            if event_structure and event_structure[-1] == 'I':
                last_motif = motifs[-1]
                ref_seq = last_motif.ref
                
                # 找到参考序列中最后一个非gap字符的位置
                last_non_gap = None
                for i in range(len(ref_seq) - 1, -1, -1):
                    if ref_seq[i] != '-':
                        last_non_gap = i
                        break
                
                # 如果结尾有gap，则修剪
                if last_non_gap is not None and last_non_gap < len(ref_seq) - 1:
                    trimmed_seq = last_motif.seq[:last_non_gap + 1]
                    trimmed_ref = last_motif.ref[:last_non_gap + 1]
                    motifs[-1] = AlignedSEQMotif(trimmed_seq, trimmed_ref)
            
            # 重新构建AlignedSEQ
            seq_segments = [motif.seq for motif in motifs]
            ref_segments = [motif.ref for motif in motifs]
            result_seqs.append(AlignedSEQ(seq_segments, ref_segments))
        
        return result_seqs[0] if is_single else result_seqs
    
    @staticmethod
    def sanitize_conserved_regions(aligned_seqs: Union[AlignedSEQ, List[AlignedSEQ]], 
                                 cutsite_motif_indices: List[int]) -> Union[AlignedSEQ, List[AlignedSEQ]]:
        """
        清理保守区域的测序错误
        
        将非cutsite motif中的错配恢复为参考序列，因为这些区域在生物学上是保守的
        
        Args:
            aligned_seqs: 单个AlignedSEQ或AlignedSEQ列表
            cutsite_motif_indices: cutsite motif的索引列表
            
        Returns:
            标准化后的AlignedSEQ或列表
        """
        # 处理单个对象的情况
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
            
            # 原始序列长度（用于验证）
            original_length = len(aligned_seq.get_seq())
            
            # 标识需要清理的motifs
            # 只清理错配('M')且不是cutsite的motifs
            motifs_to_clean = []
            for i, (event, motif) in enumerate(zip(event_structure, motifs)):
                if event == 'M' and i not in cutsite_motif_indices:
                    # 确保该motif没有gaps
                    if '-' not in motif.seq and '-' not in motif.ref:
                        motifs_to_clean.append(i)
            
            # 执行清理：将错配的序列恢复为参考序列
            for i in motifs_to_clean:
                old_motif = motifs[i]
                # 将查询序列替换为参考序列
                motifs[i] = AlignedSEQMotif(old_motif.ref, old_motif.ref)
            
            # 重新构建AlignedSEQ
            seq_segments = [motif.seq for motif in motifs]
            ref_segments = [motif.ref for motif in motifs]
            cleaned_seq = AlignedSEQ(seq_segments, ref_segments)
            
            # 验证序列长度没有改变
            new_length = len(cleaned_seq.get_seq())
            if new_length != original_length:
                raise RuntimeError(f"标准化后序列长度改变: {original_length} -> {new_length}")
            
            result_seqs.append(cleaned_seq)
        
        return result_seqs[0] if is_single else result_seqs


def desemble_sequence(aligned_query: str, aligned_ref: str, 
                     motif_boundaries: List[Tuple[int, int]]) -> AlignedSEQ:
    """
    将比对结果分解为motifs
    
    Args:
        aligned_query: 比对后的查询序列
        aligned_ref: 比对后的参考序列
        motif_boundaries: motif边界列表，每个元素为(start, end)
        
    Returns:
        AlignedSEQ对象
    """
    if len(aligned_query) != len(aligned_ref):
        raise ValueError(f"比对序列长度不匹配: {len(aligned_query)} vs {len(aligned_ref)}")
    
    seq_segments = []
    ref_segments = []
    
    for start, end in motif_boundaries:
        if start < 0 or end > len(aligned_query) or start >= end:
            raise ValueError(f"无效的motif边界: ({start}, {end})")
        
        seq_segments.append(aligned_query[start:end])
        ref_segments.append(aligned_ref[start:end])
    
    result = AlignedSEQ(seq_segments, ref_segments)
    
    # 验证分解结果
    if result.get_seq() != aligned_query:
        raise RuntimeError("分解后序列与原始序列不匹配")
    if result.get_ref() != aligned_ref:
        raise RuntimeError("分解后参考序列与原始序列不匹配")
    
    return result


def calculate_motif_boundaries(aligned_ref: str, carlin_config) -> List[Tuple[int, int]]:
    """
    计算motif边界
    
    基于比对后的参考序列和CARLIN配置计算各motif的边界
    
    Args:
        aligned_ref: 比对后的参考序列
        carlin_config: CARLIN配置对象
        
    Returns:
        motif边界列表
    """
    # 找到参考序列中非gap位置
    ref_positions = [i for i, char in enumerate(aligned_ref) if char != '-']
    
    if len(ref_positions) != len(carlin_config.carlin_sequence):
        raise ValueError(f"参考序列非gap位置数量({len(ref_positions)})与CARLIN序列长度({len(carlin_config.carlin_sequence)})不匹配")
    
    # 获取CARLIN内部的motif边界（相对于原始CARLIN序列）
    carlin_boundaries = []
    
    # Prefix
    prefix_start, prefix_end = carlin_config.positions['prefix']
    carlin_boundaries.append((prefix_start, prefix_end))
    
    # Segments和PAMs
    for i in range(10):
        # Consite
        consite_start, consite_end = carlin_config.positions['consites'][i]
        carlin_boundaries.append((consite_start, consite_end))
        
        # Cutsite
        cutsite_start, cutsite_end = carlin_config.positions['cutsites'][i]
        carlin_boundaries.append((cutsite_start, cutsite_end))
        
        # PAM (前9个segment后面有PAM)
        if i < 9:
            pam_start, pam_end = carlin_config.positions['pams'][i]
            carlin_boundaries.append((pam_start, pam_end))
    
    # Postfix
    postfix_start, postfix_end = carlin_config.positions['postfix']
    carlin_boundaries.append((postfix_start, postfix_end))
    
    # 将CARLIN边界映射到比对序列边界
    aligned_boundaries = []
    for carlin_start, carlin_end in carlin_boundaries:
        aligned_start = ref_positions[carlin_start]
        aligned_end = ref_positions[carlin_end - 1] + 1  # end是exclusive的
        aligned_boundaries.append((aligned_start, aligned_end))
    
    # 处理gaps：确保边界覆盖所有比对序列
    # 调整第一个边界的开始和最后一个边界的结束
    if aligned_boundaries:
        aligned_boundaries[0] = (0, aligned_boundaries[0][1])
        aligned_boundaries[-1] = (aligned_boundaries[-1][0], len(aligned_ref))
        
        # 处理中间的gaps
        for i in range(len(aligned_boundaries) - 1):
            current_end = aligned_boundaries[i][1]
            next_start = aligned_boundaries[i + 1][0]
            
            if current_end < next_start:
                # 有gap，将gap分配给下一个motif
                aligned_boundaries[i + 1] = (current_end, aligned_boundaries[i + 1][1])
    
    return aligned_boundaries 