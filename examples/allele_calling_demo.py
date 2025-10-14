#!/usr/bin/env python3
"""
CARLIN等位基因调用演示

展示CARLIN等位基因调用的核心功能，包括：
1. 精确调用和粗粒度调用
2. 批量调用处理
3. 统计分析
"""

import numpy as np
from typing import List, Tuple, Optional, Union

# CARLIN imports
from darlin.config import get_original_carlin_config
from darlin.alignment import create_default_aligner, AlignedSEQ, AlignedSEQMotif
from darlin.calling import AlleleCaller, AlleleCallResult, BulkAlleleCallResult, AlleleCallStatistics
from darlin.config.amplicon_configs import load_carlin_config_by_locus


def create_sample_sequences() -> List[AlignedSEQ]:
    """创建一些示例AlignedSEQ对象用于演示"""
    
    # 创建不同类型的序列变异
    sequences = []
    
    # 1. 无变化序列 (野生型)
    seq_segments = ["ATCGTCGA", "CGATCGAT", "TAGCTAGC"]
    ref_segments = ["ATCGTCGA", "CGATCGAT", "TAGCTAGC"]
    sequences.append(AlignedSEQ(seq_segments, ref_segments))
    
    # 2. 单点突变
    seq_segments = ["ATCGTCGA", "CGATGGAT", "TAGCTAGC"]  # 中间段有突变
    ref_segments = ["ATCGTCGA", "CGATCGAT", "TAGCTAGC"]
    sequences.append(AlignedSEQ(seq_segments, ref_segments))
    
    # 3. 删除突变
    seq_segments = ["ATCGTCGA", "CGAT-GAT", "TAGCTAGC"]  # 中间段有删除
    ref_segments = ["ATCGTCGA", "CGATCGAT", "TAGCTAGC"]
    sequences.append(AlignedSEQ(seq_segments, ref_segments))
    
    # 4. 插入突变
    seq_segments = ["ATCGTCGA", "CGATTCGAT", "TAGCTAGC"]  # 中间段有插入
    ref_segments = ["ATCGTCGA", "CGAT-CGAT", "TAGCTAGC"]
    sequences.append(AlignedSEQ(seq_segments, ref_segments))
    
    # 5. 复合突变
    seq_segments = ["ATCGTCGA", "CGAA-GAT", "TAGCAAAC"]  # 多处变异
    ref_segments = ["ATCGTCGA", "CGATCGAT", "TAGCTAGC"]
    sequences.append(AlignedSEQ(seq_segments, ref_segments))
    
    return sequences


def demonstrate_basic_calling():
    """演示基本的等位基因调用"""
    print("=== 基本等位基因调用演示 ===")
    
    # 使用Col1a1配置（默认）
    caller = AlleleCaller(locus="Col1a1", dominant_threshold=0.5)
    
    # 或者使用Rosa配置
    # caller = AlleleCaller(locus="Rosa", dominant_threshold=0.5)
    
    # 创建示例序列
    sequences = create_sample_sequences()
    weights = [5.0, 2.0, 1.5, 1.0, 0.5]  # 野生型权重最高
    
    print("输入序列:")
    for i, seq in enumerate(sequences):
        events = ''.join(seq.get_event_structure())
        print(f"  序列 {i+1}: 事件结构 {events}, 权重 {weights[i]}")
    
    print(f"\n主导阈值: {caller.dominant_threshold}")
    print(f"总权重: {sum(weights)}")
    print(f"野生型比例: {weights[0]/sum(weights):.2f}")
    
    # 精确调用
    print("\n--- 精确调用 ---")
    exact_result = caller.call_alleles_exact(sequences, weights, dominant_only=True)
    
    if exact_result.is_callable():
        print(f"✅ 调用成功!")
        print(f"   置信度: {exact_result.confidence:.3f}")
        print(f"   支持序列数: {exact_result.num_supporting_sequences}")
        print(f"   总权重: {exact_result.total_weight:.1f}")
        print(f"   事件结构: {exact_result.event_structure}")
    else:
        print(f"❌ 调用失败 (置信度: {exact_result.confidence:.3f})")
    
    # 粗粒度调用
    print("\n--- 粗粒度调用 ---")
    coarse_result = caller.call_alleles_coarse_grain(sequences, weights, dominant_only=True)
    
    if coarse_result.is_callable():
        print(f"✅ 调用成功!")
        print(f"   置信度: {coarse_result.confidence:.3f}")
        print(f"   支持序列数: {coarse_result.num_supporting_sequences}")
        print(f"   总权重: {coarse_result.total_weight:.1f}")
        print(f"   事件结构: {coarse_result.event_structure}")
        print(f"   共识序列长度: {len(coarse_result.allele.get_seq())}")
    else:
        print(f"❌ 调用失败 (置信度: {coarse_result.confidence:.3f})")


def simulate_umi_groups(aligned_seqs: List[AlignedSEQ], num_groups: int = 10) -> tuple:
    """模拟UMI组数据"""
    np.random.seed(42)  # 确保可重现性
    
    sequence_groups = []
    group_weights = []
    
    for i in range(num_groups):
        # 每组随机选择1-4个序列
        group_size = np.random.randint(1, 5)
        group_seqs = []
        group_w = []
        
        for j in range(group_size):
            # 从原始序列中随机选择
            seq_idx = np.random.randint(0, len(aligned_seqs))
            group_seqs.append(aligned_seqs[seq_idx])
            # 随机权重
            weight = np.random.exponential(1.0)
            group_w.append(weight)
            
        sequence_groups.append(group_seqs)
        group_weights.append(group_w)
        
    return sequence_groups, group_weights


def demonstrate_bulk_calling():
    """演示批量等位基因调用"""
    print("=== 批量等位基因调用演示 ===")
    
    # 使用Tigre配置
    caller = AlleleCaller(locus="Tigre", dominant_threshold=0.6)
    
    # 创建示例序列和模拟UMI组
    sequences = create_sample_sequences()
    sequence_groups, group_weights = simulate_umi_groups(sequences, num_groups=15)
    
    print(f"模拟了 {len(sequence_groups)} 个UMI组")
    print(f"每组序列数分布: {[len(group) for group in sequence_groups]}")
    
    # 批量调用 - 精确方法
    print("\n--- 批量精确调用 ---")
    bulk_exact = caller.call_bulk_alleles(
        sequence_groups, 
        group_weights, 
        method='exact',
        dominant_only=True
    )
    
    print(f"总尝试次数: {len(bulk_exact.individual_results)}")
    print(f"成功调用数: {bulk_exact.num_called_alleles}")
    print(f"调用成功率: {bulk_exact.calling_success_rate:.2f}")
    print(f"可调用序列总数: {bulk_exact.total_callable_sequences}")
    
    if bulk_exact.summary_alleles:
        print(f"\n发现的等位基因数: {len(bulk_exact.summary_alleles)}")
        for i, (allele, freq) in enumerate(zip(bulk_exact.summary_alleles[:3], bulk_exact.allele_frequencies[:3])):
            events = ''.join(allele.get_event_structure())
            print(f"  等位基因 {i+1}: 事件结构 {events}, 频率 {freq:.3f}")
    
    # 批量调用 - 粗粒度方法
    print("\n--- 批量粗粒度调用 ---")
    bulk_coarse = caller.call_bulk_alleles(
        sequence_groups, 
        group_weights, 
        method='coarse_grain',
        dominant_only=False  # 允许非主导等位基因
    )
    
    print(f"成功调用数: {bulk_coarse.num_called_alleles}")
    print(f"调用成功率: {bulk_coarse.calling_success_rate:.2f}")
    
    if bulk_coarse.summary_alleles:
        print(f"\n发现的等位基因数: {len(bulk_coarse.summary_alleles)}")
        for i, (allele, freq) in enumerate(zip(bulk_coarse.summary_alleles, bulk_coarse.allele_frequencies)):
            events = ''.join(allele.get_event_structure())
            print(f"  等位基因 {i+1}: 事件结构 {events}, 频率 {freq:.3f}")
    
    return bulk_exact, bulk_coarse


def demonstrate_statistical_analysis(bulk_results: List[BulkAlleleCallResult]):
    """演示统计分析功能"""
    print("\n\n=== 统计分析演示 ===\n")
    
    for i, bulk_result in enumerate(bulk_results):
        method = "精确" if i == 0 else "粗粒度"
        print(f"--- {method}调用统计 ---")
        
        # 创建统计对象
        stats = AlleleCallStatistics(bulk_result.individual_results)
        
        # 计算调用指标
        metrics = stats.compute_calling_metrics()
        print(f"调用成功率: {metrics['success_rate']:.2f}")
        print(f"平均置信度: {metrics['mean_confidence']:.3f}")
        print(f"中位数置信度: {metrics['median_confidence']:.3f}")
        print(f"平均权重: {metrics['mean_weight']:.2f}")
        print(f"平均支持序列数: {metrics['mean_supporting_sequences']:.1f}")
        
        # 事件分布
        event_dist = stats.compute_event_distribution()
        print(f"\n事件分布:")
        for event_type, count in sorted(event_dist.items()):
            print(f"  {event_type}: {count}")
        
        # 方法分布
        method_dist = stats.compute_method_distribution()
        print(f"\n方法分布:")
        for method_name, count in sorted(method_dist.items()):
            print(f"  {method_name}: {count}")
        
        print()


def main():
    """主函数"""
    print("🧬 CARLIN等位基因调用演示")
    print("=" * 50)
    
    # 基本调用演示
    demonstrate_basic_calling()
    
    # 批量调用演示
    bulk_results = demonstrate_bulk_calling()
    
    # 统计分析演示
    demonstrate_statistical_analysis(bulk_results)
    
    print("\n🎉 演示完成!")
    print("\n这个演示展示了CARLIN等位基因调用的主要功能:")
    print("  ✅ 精确调用和粗粒度调用")
    print("  ✅ 批量处理UMI组数据")
    print("  ✅ 全面的统计分析")
    print("  ✅ 灵活的阈值配置")


if __name__ == "__main__":
    main() 