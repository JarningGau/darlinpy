#!/usr/bin/env python3
"""
CARLIN突变注释模块

实现突变事件的识别、分类和HGVS格式注释功能
"""

from typing import List, Optional, Dict, Tuple, Union
from dataclasses import dataclass
from enum import Enum
import re

from ..alignment.aligned_seq import AlignedSEQ, AlignedSEQMotif


class MutationType(Enum):
    """突变类型枚举"""
    SUBSTITUTION = "M"  # 替换/错配
    INSERTION = "I"     # 插入
    DELETION = "D"      # 删除 
    COMPLEX = "C"       # 复合突变(包含多种类型)
    INDEL = "DI"        # 缺失-插入
    CONSERVED = "N"     # 保守(无变化)


@dataclass
class Mutation:
    """
    突变事件类
    
    表示单个突变事件，支持HGVS格式注释
    
    Attributes:
        type: 突变类型 (MutationType)
        loc_start: 突变起始位置 (1-based)
        loc_end: 突变结束位置 (1-based, inclusive)
        seq_old: 原始序列片段
        seq_new: 突变后序列片段
        motif_index: 发生突变的motif索引 (0-based)
        confidence: 突变调用的置信度 (0-1)
    """
    type: MutationType
    loc_start: int
    loc_end: int
    seq_old: str
    seq_new: str
    motif_index: int = -1
    confidence: float = 1.0
    
    def __post_init__(self):
        """验证突变数据"""
        if self.loc_start < 1:
            raise ValueError("突变位置必须为1-based正整数")
        if self.loc_end < self.loc_start:
            raise ValueError("结束位置不能小于起始位置")
        if self.confidence < 0 or self.confidence > 1:
            raise ValueError("置信度必须在0-1之间")
    
    @property
    def length_change(self) -> int:
        """返回突变导致的长度变化"""
        return len(self.seq_new) - len(self.seq_old)
    
    @property
    def is_indel(self) -> bool:
        """检查是否为indel突变"""
        return self.type in [MutationType.INSERTION, MutationType.DELETION, MutationType.INDEL]
    
    @property
    def affected_length(self) -> int:
        """返回受影响的序列长度"""
        return self.loc_end - self.loc_start + 1
    
    def to_hgvs(self, reference_name: str = "CARLIN") -> str:
        """
        转换为HGVS格式注释
        
        Args:
            reference_name: 参考序列名称
            
        Returns:
            str: HGVS格式的突变注释
            
        Examples:
            - 替换: g.100A>T
            - 删除: g.100_102del
            - 插入: g.100_101insACG
            - 复合: g.100_102delinsATG
        """
        prefix = f"{reference_name}:g."
        
        if self.type == MutationType.SUBSTITUTION:
            if len(self.seq_old) == 1 and len(self.seq_new) == 1:
                # 单碱基替换
                return f"{prefix}{self.loc_start}{self.seq_old}>{self.seq_new}"
            else:
                # 多碱基替换 
                return f"{prefix}{self.loc_start}_{self.loc_end}{self.seq_old}>{self.seq_new}"
        
        elif self.type == MutationType.DELETION:
            if self.loc_start == self.loc_end:
                # 单碱基删除
                return f"{prefix}{self.loc_start}del"
            else:
                # 多碱基删除
                return f"{prefix}{self.loc_start}_{self.loc_end}del"
        
        elif self.type == MutationType.INSERTION:
            # 插入突变
            if self.loc_start == self.loc_end:
                return f"{prefix}{self.loc_start}_{self.loc_start + 1}ins{self.seq_new}"
            else:
                return f"{prefix}{self.loc_start}_{self.loc_end}ins{self.seq_new}"
        
        elif self.type in [MutationType.COMPLEX, MutationType.INDEL]:
            # 复合突变或indel
            if self.seq_new:
                return f"{prefix}{self.loc_start}_{self.loc_end}delins{self.seq_new}"
            else:
                return f"{prefix}{self.loc_start}_{self.loc_end}del"
        
        elif self.type == MutationType.CONSERVED:
            return f"{prefix}{self.loc_start}_{self.loc_end}="
        
        else:
            return f"{prefix}{self.loc_start}_{self.loc_end}unknown"
    
    def __str__(self) -> str:
        """字符串表示"""
        return self.to_hgvs()
    
    def __repr__(self) -> str:
        """详细表示"""
        return (f"Mutation(type={self.type.value}, "
                f"pos={self.loc_start}-{self.loc_end}, "
                f"'{self.seq_old}'->'{self.seq_new}', "
                f"motif={self.motif_index})")


class MutationIdentifier:
    """
    突变识别器
    
    从AlignedSEQ对象中识别和注释突变事件
    """
    
    def __init__(self, min_confidence: float = 0.8):
        """
        初始化突变识别器
        
        Args:
            min_confidence: 最小置信度阈值
        """
        self.min_confidence = min_confidence
    
    def identify_sequence_events(self, aligned_seq: AlignedSEQ) -> List[Mutation]:
        """
        识别序列中的基础突变事件
        
        Args:
            aligned_seq: 比对后的序列对象
            
        Returns:
            List[Mutation]: 识别出的突变列表
        """
        mutations = []
        position = 1  # 1-based位置计数
        
        for motif_idx, motif in enumerate(aligned_seq.motifs):
            mutation = self._identify_motif_mutation(motif, position, motif_idx)
            if mutation:
                mutations.append(mutation)
            
            # 更新位置计数(基于参考序列)
            position += len(motif.ref.replace('-', ''))
        
        return mutations
    
    def _identify_motif_mutation(self, motif: AlignedSEQMotif, start_pos: int, motif_idx: int) -> Optional[Mutation]:
        """
        识别单个motif中的突变
        
        Args:
            motif: motif对象
            start_pos: motif在参考序列中的起始位置
            motif_idx: motif索引
            
        Returns:
            Optional[Mutation]: 识别出的突变，如果没有突变则返回None
        """
        seq = motif.seq
        ref = motif.ref
        
        # 移除gap并分析
        seq_nogap = seq.replace('-', '')
        ref_nogap = ref.replace('-', '')
        
        if seq == ref:
            # 完全匹配，无突变
            return None
        
        # 确定突变类型和位置
        mutation_type, seq_old, seq_new = self._classify_mutation(seq, ref)
        
        # 计算精确的位置范围
        loc_start = start_pos
        if len(ref_nogap) > 0:
            loc_end = start_pos + len(ref_nogap) - 1
        else:
            # 插入的情况，结束位置等于起始位置
            loc_end = start_pos
        
        return Mutation(
            type=mutation_type,
            loc_start=loc_start,
            loc_end=loc_end,
            seq_old=ref_nogap,
            seq_new=seq_nogap,
            motif_index=motif_idx,
            confidence=self.min_confidence
        )
    
    def _classify_mutation(self, seq: str, ref: str) -> Tuple[MutationType, str, str]:
        """
        分类突变类型
        
        Args:
            seq: 查询序列(含gap)
            ref: 参考序列(含gap)
            
        Returns:
            Tuple[MutationType, str, str]: (突变类型, 原始序列, 新序列)
        """
        seq_clean = seq.replace('-', '')
        ref_clean = ref.replace('-', '')
        
        if not seq_clean and ref_clean:
            # 纯删除
            return MutationType.DELETION, ref_clean, ''
        elif seq_clean and not ref_clean:
            # 纯插入  
            return MutationType.INSERTION, '', seq_clean
        elif len(seq_clean) == len(ref_clean):
            # 长度相等，替换
            return MutationType.SUBSTITUTION, ref_clean, seq_clean
        else:
            # 长度不等，复合突变
            if abs(len(seq_clean) - len(ref_clean)) <= 3:
                return MutationType.INDEL, ref_clean, seq_clean
            else:
                return MutationType.COMPLEX, ref_clean, seq_clean
    
    def identify_cas9_events(self, aligned_seq: AlignedSEQ, cut_sites: Optional[List[int]] = None) -> List[Mutation]:
        """
        识别Cas9特异性编辑事件
        
        Args:
            aligned_seq: 比对后的序列对象
            cut_sites: Cas9切割位点列表(1-based位置)
            
        Returns:
            List[Mutation]: Cas9特异性突变列表
        """
        # 首先获取基础突变事件
        basic_mutations = self.identify_sequence_events(aligned_seq)
        
        # 如果没有提供切割位点，返回基础突变
        if not cut_sites:
            return basic_mutations
        
        # 过滤和增强Cas9特异性事件
        cas9_mutations = []
        for mutation in basic_mutations:
            # 检查突变是否靠近切割位点
            near_cutsite = any(
                abs(mutation.loc_start - cutsite) <= 10 or
                abs(mutation.loc_end - cutsite) <= 10
                for cutsite in cut_sites
            )
            
            if near_cutsite:
                # 增加置信度
                mutation.confidence = min(1.0, mutation.confidence + 0.1)
                cas9_mutations.append(mutation)
            elif mutation.type in [MutationType.INSERTION, MutationType.DELETION]:
                # 即使不在切割位点附近，indel也可能是Cas9导致的
                mutation.confidence = max(0.5, mutation.confidence - 0.2)
                cas9_mutations.append(mutation)
        
        return cas9_mutations
    
    def merge_adjacent_mutations(self, mutations: List[Mutation], max_distance: int = 3) -> List[Mutation]:
        """
        合并相邻的突变事件
        
        Args:
            mutations: 突变列表
            max_distance: 最大合并距离
            
        Returns:
            List[Mutation]: 合并后的突变列表
        """
        if not mutations:
            return []
        
        # 按位置排序
        sorted_mutations = sorted(mutations, key=lambda m: m.loc_start)
        merged = [sorted_mutations[0]]
        
        for current in sorted_mutations[1:]:
            last = merged[-1]
            
            # 检查是否应该合并
            if current.loc_start - last.loc_end <= max_distance:
                # 合并突变
                merged_mutation = self._merge_two_mutations(last, current)
                merged[-1] = merged_mutation
            else:
                merged.append(current)
        
        return merged
    
    def _merge_two_mutations(self, mut1: Mutation, mut2: Mutation) -> Mutation:
        """
        合并两个突变事件
        
        Args:
            mut1: 第一个突变
            mut2: 第二个突变
            
        Returns:
            Mutation: 合并后的突变
        """
        # 确定合并后的位置范围
        loc_start = min(mut1.loc_start, mut2.loc_start)
        loc_end = max(mut1.loc_end, mut2.loc_end)
        
        # 合并序列信息
        seq_old = mut1.seq_old + mut2.seq_old
        seq_new = mut1.seq_new + mut2.seq_new
        
        # 确定合并后的类型
        if mut1.type == mut2.type:
            merged_type = mut1.type
        else:
            merged_type = MutationType.COMPLEX
        
        # 平均置信度
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


# 便利函数
def annotate_mutations(aligned_seq: AlignedSEQ, 
                      cas9_mode: bool = True,
                      cut_sites: Optional[List[int]] = None,
                      merge_adjacent: bool = True) -> List[Mutation]:
    """
    对比对序列进行突变注释
    
    Args:
        aligned_seq: 比对后的序列对象
        cas9_mode: 是否使用Cas9模式
        cut_sites: Cas9切割位点列表
        merge_adjacent: 是否合并相邻突变
        
    Returns:
        List[Mutation]: 注释的突变列表
    """
    identifier = MutationIdentifier()
    
    if cas9_mode:
        mutations = identifier.identify_cas9_events(aligned_seq, cut_sites)
    else:
        mutations = identifier.identify_sequence_events(aligned_seq)
    
    if merge_adjacent:
        mutations = identifier.merge_adjacent_mutations(mutations)
    
    return mutations 