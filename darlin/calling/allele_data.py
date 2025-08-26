"""
CARLIN等位基因调用数据结构

这个模块定义了等位基因调用过程中使用的核心数据结构。
"""

from typing import List, Optional, Dict, Any, Union
from dataclasses import dataclass
import numpy as np
from ..alignment.aligned_seq import AlignedSEQ


@dataclass
class AlleleCallResult:
    """
    单个等位基因调用的结果
    
    Attributes:
        allele: 调用的等位基因序列 (AlignedSEQ对象)
        constituents: 支持这个等位基因的序列索引列表
        weight_contribution: 每个支持序列的权重
        confidence: 调用的置信度 (0-1)
        calling_method: 使用的调用方法 ('exact' 或 'coarse_grain')
        dominant_fraction: 主导等位基因的比例
    """
    allele: Optional[AlignedSEQ]
    constituents: List[int]
    weight_contribution: List[float]
    confidence: float = 0.0
    calling_method: str = 'coarse_grain'
    dominant_fraction: float = 0.0
    
    def __post_init__(self):
        """验证数据一致性"""
        if self.allele is not None:
            if len(self.constituents) != len(self.weight_contribution):
                raise ValueError("constituents和weight_contribution长度必须一致")
            
    @property
    def total_weight(self) -> float:
        """返回支持这个等位基因的总权重"""
        return sum(self.weight_contribution)
        
    @property
    def num_supporting_sequences(self) -> int:
        """返回支持这个等位基因的序列数量"""
        return len(self.constituents)
        
    @property
    def event_structure(self) -> Optional[str]:
        """返回等位基因的事件结构"""
        if self.allele is None:
            return None
        return ''.join(self.allele.get_event_structure())
        
    @property
    def sequence(self) -> Optional[str]:
        """返回等位基因序列（去除gap）"""
        if self.allele is None:
            return None
        return self.allele.get_seq().replace('-', '')
        
    def is_callable(self) -> bool:
        """检查是否成功调用了等位基因"""
        return self.allele is not None


@dataclass 
class BulkAlleleCallResult:
    """
    批量等位基因调用的结果
    
    Attributes:
        individual_results: 每个UMI/组的调用结果列表
        summary_alleles: 汇总的等位基因列表
        allele_frequencies: 每个等位基因的频率
        total_callable_sequences: 可调用的序列总数
        calling_parameters: 调用参数
    """
    individual_results: List[AlleleCallResult]
    summary_alleles: List[AlignedSEQ]
    allele_frequencies: List[float]
    total_callable_sequences: int
    calling_parameters: Dict[str, Any]
    
    def __post_init__(self):
        """验证数据一致性"""
        if len(self.summary_alleles) != len(self.allele_frequencies):
            raise ValueError("summary_alleles和allele_frequencies长度必须一致")
            
    @property
    def num_called_alleles(self) -> int:
        """返回成功调用的等位基因数量"""
        return len([r for r in self.individual_results if r.is_callable()])
        
    @property
    def calling_success_rate(self) -> float:
        """返回调用成功率"""
        if not self.individual_results:
            return 0.0
        return self.num_called_alleles / len(self.individual_results)
        
    @property
    def dominant_allele(self) -> Optional[AlignedSEQ]:
        """返回最主导的等位基因"""
        if not self.summary_alleles:
            return None
        return self.summary_alleles[0]
        
    @property 
    def dominant_frequency(self) -> float:
        """返回最主导等位基因的频率"""
        if not self.allele_frequencies:
            return 0.0
        return self.allele_frequencies[0]
        
    def get_allele_by_rank(self, rank: int) -> Optional[AlignedSEQ]:
        """根据频率排名获取等位基因"""
        if rank < 0 or rank >= len(self.summary_alleles):
            return None
        return self.summary_alleles[rank]
        
    def get_frequency_by_rank(self, rank: int) -> float:
        """根据频率排名获取等位基因频率"""
        if rank < 0 or rank >= len(self.allele_frequencies):
            return 0.0
        return self.allele_frequencies[rank]


class AlleleCallStatistics:
    """
    等位基因调用统计信息
    """
    
    def __init__(self, results: Union[AlleleCallResult, BulkAlleleCallResult, List[AlleleCallResult]]):
        """
        初始化统计计算器
        
        Args:
            results: 等位基因调用结果
        """
        if isinstance(results, AlleleCallResult):
            self.individual_results = [results]
        elif isinstance(results, BulkAlleleCallResult):
            self.individual_results = results.individual_results
        elif isinstance(results, list):
            self.individual_results = results
        else:
            raise ValueError("不支持的结果类型")
            
    def compute_calling_metrics(self) -> Dict[str, float]:
        """
        计算调用指标
        
        Returns:
            包含各种统计指标的字典
        """
        callable_results = [r for r in self.individual_results if r.is_callable()]
        
        if not self.individual_results:
            return {
                'success_rate': 0.0,
                'mean_confidence': 0.0,
                'mean_weight': 0.0,
                'mean_supporting_sequences': 0.0
            }
            
        metrics = {
            'success_rate': len(callable_results) / len(self.individual_results),
            'total_attempts': len(self.individual_results),
            'successful_calls': len(callable_results)
        }
        
        if callable_results:
            metrics.update({
                'mean_confidence': np.mean([r.confidence for r in callable_results]),
                'median_confidence': np.median([r.confidence for r in callable_results]),
                'mean_weight': np.mean([r.total_weight for r in callable_results]),
                'median_weight': np.median([r.total_weight for r in callable_results]),
                'mean_supporting_sequences': np.mean([r.num_supporting_sequences for r in callable_results]),
                'median_supporting_sequences': np.median([r.num_supporting_sequences for r in callable_results])
            })
        else:
            metrics.update({
                'mean_confidence': 0.0,
                'median_confidence': 0.0,
                'mean_weight': 0.0,
                'median_weight': 0.0,
                'mean_supporting_sequences': 0.0,
                'median_supporting_sequences': 0.0
            })
            
        return metrics
        
    def compute_event_distribution(self) -> Dict[str, int]:
        """
        计算事件类型分布
        
        Returns:
            事件类型到计数的映射
        """
        event_counts = {}
        for result in self.individual_results:
            if result.is_callable():
                event_structure = result.event_structure
                if event_structure:
                    event_counts[event_structure] = event_counts.get(event_structure, 0) + 1
                    
        return event_counts
        
    def compute_method_distribution(self) -> Dict[str, int]:
        """
        计算调用方法分布
        
        Returns:
            调用方法到计数的映射
        """
        method_counts = {}
        for result in self.individual_results:
            if result.is_callable():
                method = result.calling_method
                method_counts[method] = method_counts.get(method, 0) + 1
                
        return method_counts 