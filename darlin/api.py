#!/usr/bin/env python3
"""
DARLIN主API模块

提供简化的高级接口来进行CARLIN序列分析
"""

from typing import List, Optional, Dict, Any, Union, Tuple
from dataclasses import dataclass, field
import time
from pathlib import Path

from .config.amplicon_configs import ORIGINAL_CARLIN, AmpliconConfig
from .alignment.carlin_aligner import CARLINAligner
from .calling.allele_caller import AlleleCaller, AlleleCallResult, BulkAlleleCallResult
from .mutations.mutation import Mutation, annotate_mutations


@dataclass
class AnalysisResult:
    """
    CARLIN序列分析结果
    
    包含完整的分析结果，包括等位基因调用、突变注释和统计信息
    
    Attributes:
        called_alleles: 调用的等位基因列表
        mutations: 每个等位基因对应的突变列表
        alignment_scores: 每个序列的比对得分
        summary_stats: 汇总统计信息
        processing_time: 处理耗时(秒)
        config_used: 使用的配置名称
        method_used: 使用的调用方法
    """
    called_alleles: List[AlleleCallResult]
    mutations: List[List[Mutation]]
    alignment_scores: List[float]
    summary_stats: Dict[str, Any]
    processing_time: float = 0.0
    config_used: str = "OriginalCARLIN"
    method_used: str = "coarse_grain"
    
    def __post_init__(self):
        """验证结果数据一致性"""
        if len(self.called_alleles) != len(self.mutations):
            raise ValueError("called_alleles和mutations长度必须一致")
    
    @property
    def num_sequences(self) -> int:
        """返回分析的序列总数"""
        return len(self.alignment_scores)
    
    @property
    def num_called_alleles(self) -> int:
        """返回成功调用的等位基因数量"""
        return len([a for a in self.called_alleles if a.is_callable()])
    
    @property
    def calling_success_rate(self) -> float:
        """返回等位基因调用成功率"""
        if not self.called_alleles:
            return 0.0
        return self.num_called_alleles / len(self.called_alleles)
    
    @property
    def total_mutations(self) -> int:
        """返回检测到的突变总数"""
        return sum(len(muts) for muts in self.mutations)
    
    @property
    def average_alignment_score(self) -> float:
        """返回平均比对得分"""
        if not self.alignment_scores:
            return 0.0
        return sum(self.alignment_scores) / len(self.alignment_scores)
    
    def get_mutation_summary(self) -> Dict[str, int]:
        """获取突变类型统计"""
        mutation_counts = {}
        for mut_list in self.mutations:
            for mut in mut_list:
                mut_type = mut.type.value
                mutation_counts[mut_type] = mutation_counts.get(mut_type, 0) + 1
        return mutation_counts
    
    def print_summary(self):
        """打印分析结果摘要"""
        print(f"CARLIN序列分析结果摘要")
        print(f"=" * 40)
        print(f"配置: {self.config_used}")
        print(f"方法: {self.method_used}")
        print(f"处理时间: {self.processing_time:.2f}秒")
        print(f"")
        print(f"序列统计:")
        print(f"  总序列数: {self.num_sequences}")
        print(f"  成功调用等位基因: {self.num_called_alleles}")
        print(f"  调用成功率: {self.calling_success_rate:.1%}")
        print(f"  平均比对得分: {self.average_alignment_score:.2f}")
        print(f"")
        print(f"突变统计:")
        print(f"  突变总数: {self.total_mutations}")
        mut_summary = self.get_mutation_summary()
        for mut_type, count in mut_summary.items():
            print(f"  {mut_type}类型: {count}")


def analyze_sequences(
    sequences: List[str],
    config: Union[str, AmpliconConfig] = 'Col1a1',  # 改为默认Col1a1
    method: str = 'coarse_grain',
    min_reads: int = 1,
    dominant_threshold: float = 0.5,
    annotate_mutations_flag: bool = True,
    merge_adjacent_mutations: bool = True,
    verbose: bool = False
) -> AnalysisResult:
    """
    分析CARLIN序列
    
    Args:
        sequences: 待分析的序列列表
        config: 配置，可以是locus名称字符串或AmpliconConfig对象，默认为'Col1a1'
        method: 调用方法 ('exact' 或 'coarse_grain')
        min_reads: 最小读取数
        dominant_threshold: 主导等位基因阈值
        annotate_mutations_flag: 是否注释突变
        merge_adjacent_mutations: 是否合并相邻突变
        verbose: 是否显示详细信息
        
    Returns:
        分析结果
    """
    start_time = time.time()
    
    if verbose:
        print(f"开始分析 {len(sequences)} 条序列...")
    
    # 处理配置参数
    if isinstance(config, str):
        # 如果是字符串，作为locus处理
        if config == 'OriginalCARLIN':
            # 向后兼容性：如果用户明确指定OriginalCARLIN，使用原来的方式
            from .config.amplicon_configs import get_original_carlin_config
            amplicon_config = get_original_carlin_config()
        else:
            # 作为locus处理
            from .config.amplicon_configs import load_carlin_config_by_locus
            amplicon_config = load_carlin_config_by_locus(config)
    else:
        # 直接使用提供的AmpliconConfig
        amplicon_config = config
    
    if verbose:
        print(f"使用配置: {config}")
    
    if method not in ['exact', 'coarse_grain']:
        raise ValueError(f"不支持的调用方法: {method}")
    
    if not sequences:
        raise ValueError("序列列表不能为空")
    
    # 过滤序列
    valid_sequences = [seq for seq in sequences if seq and len(seq) >= 50]
    if len(valid_sequences) < len(sequences):
        if verbose:
            print(f"过滤了 {len(sequences) - len(valid_sequences)} 条短序列")
    
    if not valid_sequences:
        raise ValueError("没有有效的序列可供分析")
    
    try:
        # 2. 序列比对
        if verbose:
            print("正在进行序列比对...")
        
        aligner = CARLINAligner(amplicon_config)
        alignment_results = aligner.align_sequences(valid_sequences)
        
        # 提取比对得分和序列
        alignment_scores = [result['alignment_score'] for result in alignment_results]
        aligned_sequences = [result.get('aligned_seq_obj') for result in alignment_results]
        
        if verbose:
            avg_score = sum(alignment_scores) / len(alignment_scores)
            print(f"比对完成，平均得分: {avg_score:.2f}")
        
        # 3. 等位基因调用
        if verbose:
            print("正在进行等位基因调用...")
        
        # 修改AlleleCaller的初始化
        caller = AlleleCaller(amplicon_config=amplicon_config, dominant_threshold=dominant_threshold)
        
        # 根据方法选择调用策略
        called_alleles = []
        for aligned_seq in aligned_sequences:
            if method == 'exact':
                allele_result = caller.call_alleles_exact([aligned_seq])
            else:  # coarse_grain
                allele_result = caller.call_alleles_coarse_grain([aligned_seq])
            
            called_alleles.append(allele_result)
        
        if verbose:
            callable_count = len([a for a in called_alleles if a.is_callable()])
            print(f"等位基因调用完成，成功调用: {callable_count}/{len(called_alleles)}")
        
        # 4. 突变注释
        mutations_list = []
        if annotate_mutations_flag:
            if verbose:
                print("正在进行突变注释...")
            
            # 获取CARLIN的切割位点信息
            cut_sites = _get_cut_sites(amplicon_config)
            
            for allele_result in called_alleles:
                if allele_result.is_callable() and allele_result.allele:
                    mutations = annotate_mutations(
                        allele_result.allele,
                        cas9_mode=True,
                        cut_sites=cut_sites,
                        merge_adjacent=merge_adjacent_mutations
                    )
                else:
                    mutations = []
                mutations_list.append(mutations)
            
            if verbose:
                total_mutations = sum(len(muts) for muts in mutations_list)
                print(f"突变注释完成，共检测到 {total_mutations} 个突变")
        else:
            mutations_list = [[] for _ in called_alleles]
        
        # 5. 生成统计信息
        summary_stats = _generate_summary_stats(
            valid_sequences, alignment_results, called_alleles, mutations_list
        )
        
        # 6. 构建结果对象
        processing_time = time.time() - start_time
        
        result = AnalysisResult(
            called_alleles=called_alleles,
            mutations=mutations_list,
            alignment_scores=alignment_scores,
            summary_stats=summary_stats,
            processing_time=processing_time,
            config_used=str(config),
            method_used=method
        )
        
        if verbose:
            print(f"分析完成，耗时: {processing_time:.2f}秒")
            result.print_summary()
        
        return result
        
    except Exception as e:
        raise RuntimeError(f"分析过程中发生错误: {str(e)}") from e


def _load_config(config: Union[str, AmpliconConfig]) -> AmpliconConfig:
    """加载扩增子配置"""
    if isinstance(config, AmpliconConfig):
        return config
    
    config_map = {
        'OriginalCARLIN': ORIGINAL_CARLIN,
        'original': ORIGINAL_CARLIN
    }
    
    if config not in config_map:
        available = ', '.join(config_map.keys())
        raise ValueError(f"不支持的配置名称: {config}. 可用配置: {available}")
    
    return config_map[config]


def _get_cut_sites(amplicon_config: AmpliconConfig) -> List[int]:
    """获取CARLIN的Cas9切割位点"""
    # 基于CARLIN的设计，切割位点通常在segments之间
    # 这里提供一个简化的实现，实际位点可能需要根据具体配置调整
    cut_sites = []
    position = len(amplicon_config.sequence.prefix)
    
    for i, segment in enumerate(amplicon_config.sequence.segments):
        position += len(segment)
        if i < len(amplicon_config.sequence.segments) - 1:
            # 在segments之间添加潜在的切割位点
            cut_sites.append(position)
    
    return cut_sites


def _generate_summary_stats(
    sequences: List[str],
    alignment_results: List,
    called_alleles: List[AlleleCallResult],
    mutations_list: List[List[Mutation]]
) -> Dict[str, Any]:
    """生成汇总统计信息"""
    
    # 基础统计
    stats = {
        'total_sequences': len(sequences),
        'avg_sequence_length': sum(len(seq) for seq in sequences) / len(sequences),
        'called_alleles_count': len([a for a in called_alleles if a.is_callable()]),
        'calling_success_rate': len([a for a in called_alleles if a.is_callable()]) / len(called_alleles),
    }
    
    # 比对统计
    scores = [r['alignment_score'] for r in alignment_results]
    stats.update({
        'avg_alignment_score': sum(scores) / len(scores),
        'min_alignment_score': min(scores),
        'max_alignment_score': max(scores),
    })
    
    # 突变统计
    mutation_counts = {}
    total_mutations = 0
    for mut_list in mutations_list:
        total_mutations += len(mut_list)
        for mut in mut_list:
            mut_type = mut.type.value
            mutation_counts[mut_type] = mutation_counts.get(mut_type, 0) + 1
    
    stats.update({
        'total_mutations': total_mutations,
        'avg_mutations_per_allele': total_mutations / len(called_alleles) if called_alleles else 0,
        'mutation_type_distribution': mutation_counts
    })
    
    return stats


# 便利函数
def quick_analyze(sequences: List[str], config: str = 'OriginalCARLIN') -> AnalysisResult:
    """
    快速分析函数，使用默认参数
    
    Args:
        sequences: 序列列表
        config: 配置名称
        
    Returns:
        AnalysisResult: 分析结果
    """
    return analyze_sequences(sequences, config=config, verbose=True)


def batch_analyze_files(file_paths: List[Union[str, Path]], **kwargs) -> List[AnalysisResult]:
    """
    批量分析多个文件
    
    Args:
        file_paths: 文件路径列表
        **kwargs: 传递给analyze_sequences的其他参数
        
    Returns:
        List[AnalysisResult]: 每个文件的分析结果列表
    """
    results = []
    
    for file_path in file_paths:
        # 这里需要实现文件读取逻辑
        # 暂时留作占位符，后续可以扩展支持FASTA/FASTQ等格式
        pass
    
    return results 