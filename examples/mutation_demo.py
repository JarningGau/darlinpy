#!/usr/bin/env python3
"""
CARLIN突变注释功能演示

展示如何使用突变注释模块识别和分析CARLIN序列中的CRISPR编辑事件
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from darlinpy.mutations import (
    Mutation, 
    MutationType, 
    MutationIdentifier, 
    annotate_mutations
)
from darlinpy.alignment.aligned_seq import AlignedSEQ
from darlinpy.config.amplicon_configs import ORIGINAL_CARLIN


def demo_basic_mutations():
    """演示基础突变类型的识别"""
    print("=== 基础突变类型演示 ===\n")
    
    # 1. 替换突变
    print("1. 替换突变示例:")
    substitution_seq = AlignedSEQ(
        seq_segments=["ACGTACGTACGTACGTACGT", "TGCATGCATGCATGCATGCA"],
        ref_segments=["ACGTACGTACGTACGTACGT", "AGCATGCATGCATGCATGCA"]  # T->A替换
    )
    
    mutations = annotate_mutations(substitution_seq, cas9_mode=False)
    for mut in mutations:
        print(f"  突变类型: {mut.type.value}")
        print(f"  HGVS格式: {mut.to_hgvs()}")
        print(f"  位置: {mut.loc_start}-{mut.loc_end}")
        print(f"  序列变化: '{mut.seq_old}' -> '{mut.seq_new}'")
        print()
    
    # 2. 删除突变
    print("2. 删除突变示例:")
    deletion_seq = AlignedSEQ(
        seq_segments=["ACGTACGTACGTACGTACGT", "--------------------"],  # 完全删除
        ref_segments=["ACGTACGTACGTACGTACGT", "TGCATGCATGCATGCATGCA"]
    )
    
    mutations = annotate_mutations(deletion_seq, cas9_mode=False)
    for mut in mutations:
        print(f"  突变类型: {mut.type.value}")
        print(f"  HGVS格式: {mut.to_hgvs()}")
        print(f"  长度变化: {mut.length_change} bp")
        print(f"  是否为indel: {mut.is_indel}")
        print()
    
    # 3. 插入突变
    print("3. 插入突变示例:")
    insertion_seq = AlignedSEQ(
        seq_segments=["ACGTACGTACGTACGTACGT", "TGCATGCATGCATGCATGCA"],
        ref_segments=["ACGTACGTACGTACGTACGT", "--------------------"]  # 插入
    )
    
    mutations = annotate_mutations(insertion_seq, cas9_mode=False)
    for mut in mutations:
        print(f"  突变类型: {mut.type.value}")
        print(f"  HGVS格式: {mut.to_hgvs()}")
        print(f"  长度变化: {mut.length_change} bp")
        print()


def demo_cas9_specific_analysis():
    """演示Cas9特异性分析"""
    print("=== Cas9特异性分析演示 ===\n")
    
    # 模拟CARLIN序列中的Cas9编辑事件
    print("模拟CARLIN序列的Cas9编辑分析:")
    
    # 创建包含多种编辑事件的序列
    cas9_edited_seq = AlignedSEQ(
        seq_segments=[
            "ACGTACGTACGTACGTACGT",  # segment 1: 无变化
            "--------------------",  # segment 2: 完全删除
            "TGCATGCATGCATGCATGCA",  # segment 3: 无变化
            "GGCCGGCCGGCCGGCCGGCC",  # segment 4: 替换
            "ATATAT--GCGCGCGCGCGC",  # segment 5: 部分删除
        ],
        ref_segments=[
            "ACGTACGTACGTACGTACGT",  # segment 1
            "TTTTTTTTTTTTTTTTTTTT",  # segment 2
            "TGCATGCATGCATGCATGCA",  # segment 3
            "AAAAAAAAAAAAAAAAAAAA",  # segment 4 
            "ATATATATGCGCGCGCGCGC",  # segment 5
        ]
    )
    
    # 定义CARLIN的切割位点(示例位置)
    cut_sites = [25, 45, 85]  # 在segment边界附近
    
    # 使用Cas9模式分析
    mutations = annotate_mutations(
        cas9_edited_seq, 
        cas9_mode=True, 
        cut_sites=cut_sites,
        merge_adjacent=True
    )
    
    print(f"检测到 {len(mutations)} 个突变事件:\n")
    
    for i, mut in enumerate(mutations, 1):
        print(f"突变 {i}:")
        print(f"  类型: {mut.type.value}")
        print(f"  位置: {mut.loc_start}-{mut.loc_end}")
        print(f"  HGVS: {mut.to_hgvs()}")
        print(f"  置信度: {mut.confidence:.2f}")
        print(f"  Motif索引: {mut.motif_index}")
        
        # 检查是否靠近切割位点
        near_cutsite = any(
            abs(mut.loc_start - cs) <= 10 or abs(mut.loc_end - cs) <= 10 
            for cs in cut_sites
        )
        print(f"  靠近切割位点: {'是' if near_cutsite else '否'}")
        print()


def demo_complex_scenarios():
    """演示复杂的突变场景"""
    print("=== 复杂突变场景演示 ===\n")
    
    # 1. 复合突变
    print("1. 复合突变 (删除+插入):")
    complex_seq = AlignedSEQ(
        seq_segments=[
            "ACGTACGTACGTACGTACGT",
            "ATGGCCTTAAGGCCTT----"  # 部分替换+删除
        ],
        ref_segments=[
            "ACGTACGTACGTACGTACGT", 
            "TGCATGCATGCATGCATGCA"
        ]
    )
    
    mutations = annotate_mutations(complex_seq, cas9_mode=False)
    for mut in mutations:
        print(f"  类型: {mut.type.value}")
        print(f"  HGVS: {mut.to_hgvs()}")
        print(f"  长度变化: {mut.length_change} bp")
        print()
    
    # 2. 相邻突变的合并
    print("2. 相邻突变合并演示:")
    
    # 创建相邻的突变
    mut1 = Mutation(
        type=MutationType.DELETION,
        loc_start=10,
        loc_end=12,
        seq_old="ATG",
        seq_new="",
        confidence=0.95
    )
    
    mut2 = Mutation(
        type=MutationType.INSERTION,
        loc_start=14,
        loc_end=14,
        seq_old="",
        seq_new="GCC",
        confidence=0.90
    )
    
    identifier = MutationIdentifier()
    
    print("  合并前:")
    print(f"    突变1: {mut1.to_hgvs()}")
    print(f"    突变2: {mut2.to_hgvs()}")
    
    merged = identifier.merge_adjacent_mutations([mut1, mut2], max_distance=3)
    
    print("  合并后:")
    for mut in merged:
        print(f"    合并突变: {mut.to_hgvs()}")
        print(f"    类型: {mut.type.value}")
        print(f"    置信度: {mut.confidence:.2f}")
    print()


def demo_hgvs_formats():
    """演示各种HGVS格式"""
    print("=== HGVS格式演示 ===\n")
    
    # 不同类型突变的HGVS格式示例
    examples = [
        # (type, start, end, old, new, description)
        (MutationType.SUBSTITUTION, 100, 100, "A", "T", "单碱基替换"),
        (MutationType.SUBSTITUTION, 100, 102, "ATG", "GCC", "多碱基替换"),
        (MutationType.DELETION, 100, 100, "A", "", "单碱基删除"),
        (MutationType.DELETION, 100, 105, "ATGGCC", "", "多碱基删除"),
        (MutationType.INSERTION, 100, 100, "", "GCC", "插入"),
        (MutationType.COMPLEX, 100, 102, "ATG", "GCCAATT", "复合突变"),
        (MutationType.INDEL, 100, 101, "AT", "GCCG", "Indel"),
    ]
    
    for mut_type, start, end, old, new, desc in examples:
        mutation = Mutation(
            type=mut_type,
            loc_start=start,
            loc_end=end,
            seq_old=old,
            seq_new=new
        )
        
        print(f"{desc}:")
        print(f"  HGVS格式: {mutation.to_hgvs()}")
        print(f"  长度变化: {mutation.length_change} bp")
        print(f"  受影响长度: {mutation.affected_length} bp")
        print()


def main():
    """主演示函数"""
    print("CARLIN突变注释功能演示")
    print("=" * 50)
    print()
    
    try:
        demo_basic_mutations()
        demo_cas9_specific_analysis()
        demo_complex_scenarios()
        demo_hgvs_formats()
        
        print("演示完成！")
        print("\n突变注释模块功能包括：")
        print("✓ 多种突变类型识别 (替换、插入、删除、复合)")
        print("✓ HGVS标准格式注释")
        print("✓ Cas9特异性事件分析") 
        print("✓ 相邻突变智能合并")
        print("✓ 置信度评估")
        print("✓ 位置和长度变化统计")
        
    except Exception as e:
        print(f"演示过程中出现错误: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main() 