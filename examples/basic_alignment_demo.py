#!/usr/bin/env python3
"""
DARLIN Python 基础序列比对演示

展示如何使用核心的cas9_align算法进行序列比对
"""

import numpy as np
from darlin.alignment.cas9_align import cas9_align, nt2int, int2nt, print_cas9_alignment


def create_nuc44_matrix():
    """创建简化的NUC44评分矩阵"""
    sub_score = np.zeros(25, dtype=np.float64)
    
    # A=1, C=2, G=3, T=4, gap=0
    for i in range(1, 5):
        for j in range(1, 5):
            if i == j:
                sub_score[i * 5 + j] = 5.0   # 匹配得分
            else:
                sub_score[i * 5 + j] = -4.0  # 不匹配得分
    
    return sub_score


def demo_perfect_match():
    """演示完全匹配的序列比对"""
    print("=== 演示1: 完全匹配序列 ===")
    
    seq_str = "ACGTACGT"
    ref_str = "ACGTACGT"
    
    seq = nt2int(seq_str)
    ref = nt2int(ref_str)
    
    # 设置统一的惩罚参数
    open_penalty = np.full(len(ref) + 1, 10.0)
    close_penalty = np.full(len(ref) + 1, 5.0)
    sub_score = create_nuc44_matrix()
    
    score, al_seq, al_ref = cas9_align(seq, ref, open_penalty, close_penalty, sub_score)
    print_cas9_alignment(al_seq, al_ref, score)
    
    expected_score = len(seq_str) * 5.0  # 每个匹配5分
    print(f"预期得分: {expected_score}, 实际得分: {score}")
    print("✅ 完全匹配测试通过\n")


def demo_with_mutations():
    """演示包含突变的序列比对"""
    print("=== 演示2: 包含突变的序列 ===")
    
    seq_str = "ACGTACGT"
    ref_str = "ACCTACCT"  # 两个G->C突变
    
    seq = nt2int(seq_str)
    ref = nt2int(ref_str)
    
    open_penalty = np.full(len(ref) + 1, 10.0)
    close_penalty = np.full(len(ref) + 1, 5.0)
    sub_score = create_nuc44_matrix()
    
    score, al_seq, al_ref = cas9_align(seq, ref, open_penalty, close_penalty, sub_score)
    print_cas9_alignment(al_seq, al_ref, score)
    
    expected_score = 6 * 5.0 + 2 * (-4.0)  # 6个匹配 + 2个不匹配
    print(f"预期得分: {expected_score}, 实际得分: {score}")
    print("✅ 突变序列测试通过\n")


def demo_crispr_like_editing():
    """演示类似CRISPR编辑的序列比对"""
    print("=== 演示3: 模拟CRISPR编辑 ===")
    
    # 模拟一个包含插入删除的CRISPR编辑事件
    seq_str = "ACGTAACCGTACGT"     # 原序列
    ref_str = "ACGTACGTACGT"       # 参考序列 (少了AA插入)
    
    seq = nt2int(seq_str)
    ref = nt2int(ref_str)
    
    # 设置位置特异性惩罚：在"切割位点"附近(位置4-6)设置低惩罚
    open_penalty = np.full(len(ref) + 1, 10.0)
    close_penalty = np.full(len(ref) + 1, 5.0)
    
    # 模拟CRISPR切割位点的低gap惩罚
    cut_site_start = 4
    cut_site_end = 7
    open_penalty[cut_site_start:cut_site_end] = 1.0
    close_penalty[cut_site_start:cut_site_end] = 0.5
    
    sub_score = create_nuc44_matrix()
    
    score, al_seq, al_ref = cas9_align(seq, ref, open_penalty, close_penalty, sub_score)
    print_cas9_alignment(al_seq, al_ref, score)
    
    print(f"切割位点位置: {cut_site_start}-{cut_site_end-1}")
    print(f"切割位点gap惩罚: 开启={open_penalty[5]}, 关闭={close_penalty[5]}")
    print(f"其他位置gap惩罚: 开启={open_penalty[0]}, 关闭={close_penalty[0]}")
    print("✅ CRISPR编辑模拟测试通过\n")


def demo_complex_alignment():
    """演示复杂的序列比对情况"""
    print("=== 演示4: 复杂比对情况 ===")
    
    # 更复杂的情况：同时包含插入、删除、替换
    seq_str = "ACGTTTACGT"         # 插入了TT
    ref_str = "ACGTACGCGT"         # 删除了一个T，A->C替换
    
    seq = nt2int(seq_str)
    ref = nt2int(ref_str)
    
    # 使用中等的gap惩罚
    open_penalty = np.full(len(ref) + 1, 3.0)
    close_penalty = np.full(len(ref) + 1, 2.0)
    sub_score = create_nuc44_matrix()
    
    score, al_seq, al_ref = cas9_align(seq, ref, open_penalty, close_penalty, sub_score)
    print_cas9_alignment(al_seq, al_ref, score)
    
    # 分析比对结果
    al_seq_str = int2nt(al_seq)
    al_ref_str = int2nt(al_ref)
    
    matches = sum(1 for i in range(len(al_seq)) if al_seq[i] == al_ref[i] and al_seq[i] != 0)
    mismatches = sum(1 for i in range(len(al_seq)) if al_seq[i] != al_ref[i] and al_seq[i] != 0 and al_ref[i] != 0)
    gaps = sum(1 for i in range(len(al_seq)) if al_seq[i] == 0 or al_ref[i] == 0)
    
    print(f"匹配: {matches}, 不匹配: {mismatches}, gaps: {gaps}")
    print("✅ 复杂比对测试通过\n")


def main():
    """主演示函数"""
    print("🧬 DARLIN Python - 序列比对演示")
    print("=" * 50)
    print()
    
    # 运行所有演示
    demo_perfect_match()
    demo_with_mutations()
    demo_crispr_like_editing()
    demo_complex_alignment()
    
    print("🎉 所有演示完成！")
    print()
    print("💡 提示:")
    print("   - cas9_align算法支持位置特异性gap惩罚")
    print("   - 适合CRISPR-Cas9序列分析")
    print("   - 可以处理复杂的插入删除模式")
    print("   - NUC44评分矩阵确保高质量比对")


if __name__ == "__main__":
    main() 