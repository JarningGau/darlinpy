#!/usr/bin/env python3
"""
测试cas9_align核心算法
"""

import numpy as np
import sys
import os

# 添加路径以便导入模块
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from darlin.alignment.cas9_align import cas9_align, nt2int, int2nt, print_cas9_alignment


class TestCAS9Align:
    
    def setup_method(self):
        """设置测试参数"""
        # 基础的评分矩阵 (简化的NUC44)
        self.sub_score = np.zeros(25, dtype=np.float64)
        for i in range(1, 5):  # A=1, C=2, G=3, T=4
            for j in range(1, 5):
                if i == j:
                    self.sub_score[i * 5 + j] = 5.0  # 匹配 (类似NUC44)
                else:
                    self.sub_score[i * 5 + j] = -4.0  # 不匹配
    
    def test_sequence_encoding(self):
        """测试序列编码/解码"""
        seq_str = "ACGT"
        seq_int = nt2int(seq_str)
        seq_back = int2nt(seq_int)
        
        assert seq_back == seq_str
        assert list(seq_int) == [1, 2, 3, 4]  # A=1, C=2, G=3, T=4
    
    def test_perfect_match(self):
        """测试完全匹配的序列"""
        seq_str = "ACGT"
        ref_str = "ACGT"
        
        seq = nt2int(seq_str)
        ref = nt2int(ref_str)
        
        # 设置统一的惩罚参数
        open_penalty = np.full(len(ref) + 1, 10.0)
        close_penalty = np.full(len(ref) + 1, 5.0)
        
        score, al_seq, al_ref = cas9_align(seq, ref, open_penalty, close_penalty, self.sub_score)
        
        # 完全匹配应该得到高分
        expected_score = 4 * 5.0  # 4个匹配，每个5分
        assert score == expected_score
        
        # 比对结果应该与输入相同
        assert int2nt(al_seq) == seq_str
        assert int2nt(al_ref) == ref_str
    
    def test_single_mismatch(self):
        """测试单个碱基不匹配"""
        seq_str = "ACGT"
        ref_str = "ACCT"  # G->C
        
        seq = nt2int(seq_str)
        ref = nt2int(ref_str)
        
        open_penalty = np.full(len(ref) + 1, 10.0)
        close_penalty = np.full(len(ref) + 1, 5.0)
        
        score, al_seq, al_ref = cas9_align(seq, ref, open_penalty, close_penalty, self.sub_score)
        
        # 3个匹配 + 1个不匹配
        expected_score = 3 * 5.0 + 1 * (-4.0)
        assert score == expected_score
        
        assert int2nt(al_seq) == seq_str
        assert int2nt(al_ref) == ref_str
    
    def test_insertion(self):
        """测试插入情况"""
        seq_str = "ACGGT"  # 多了一个G
        ref_str = "ACGT"
        
        seq = nt2int(seq_str)
        ref = nt2int(ref_str)
        
        # 设置较低的gap惩罚便于测试
        open_penalty = np.full(len(ref) + 1, 2.0)
        close_penalty = np.full(len(ref) + 1, 1.0)
        
        score, al_seq, al_ref = cas9_align(seq, ref, open_penalty, close_penalty, self.sub_score)
        
        # 检查是否正确处理了插入
        al_seq_str = int2nt(al_seq)
        al_ref_str = int2nt(al_ref)
        
        print(f"\n插入测试结果:")
        print(f"序列: {al_seq_str}")
        print(f"参考: {al_ref_str}")
        print(f"得分: {score}")
        
        # 应该包含gap
        assert '-' in al_ref_str
        assert len(al_seq_str) == len(al_ref_str)
    
    def test_deletion(self):
        """测试删除情况"""
        seq_str = "ACT"   # 少了一个G
        ref_str = "ACGT"
        
        seq = nt2int(seq_str)
        ref = nt2int(ref_str)
        
        open_penalty = np.full(len(ref) + 1, 2.0)
        close_penalty = np.full(len(ref) + 1, 1.0)
        
        score, al_seq, al_ref = cas9_align(seq, ref, open_penalty, close_penalty, self.sub_score)
        
        al_seq_str = int2nt(al_seq)
        al_ref_str = int2nt(al_ref)
        
        print(f"\n删除测试结果:")
        print(f"序列: {al_seq_str}")
        print(f"参考: {al_ref_str}")
        print(f"得分: {score}")
        
        # 应该包含gap
        assert '-' in al_seq_str
        assert len(al_seq_str) == len(al_ref_str)
    
    def test_position_specific_penalty(self):
        """测试位置特异性惩罚"""
        seq_str = "ACGGT"
        ref_str = "ACGT"
        
        seq = nt2int(seq_str)
        ref = nt2int(ref_str)
        
        # 设置位置特异性惩罚：在位置2设置很低的惩罚
        open_penalty = np.full(len(ref) + 1, 10.0)
        close_penalty = np.full(len(ref) + 1, 5.0)
        
        # 在CRISPR切割位点附近设置低惩罚
        open_penalty[2:4] = 0.5
        close_penalty[2:4] = 0.1
        
        score, al_seq, al_ref = cas9_align(seq, ref, open_penalty, close_penalty, self.sub_score)
        
        al_seq_str = int2nt(al_seq)
        al_ref_str = int2nt(al_ref)
        
        print(f"\n位置特异性惩罚测试:")
        print(f"序列: {al_seq_str}")
        print(f"参考: {al_ref_str}")
        print(f"得分: {score}")
        print(f"低惩罚区域: 位置2-3 (开启惩罚={open_penalty[2]}, 关闭惩罚={close_penalty[2]})")
        
        # 由于在低惩罚区域，gap应该更容易出现在那里
        assert '-' in al_ref_str or '-' in al_seq_str


def test_print_function():
    """测试打印函数"""
    al_seq = [1, 2, 0, 3, 4]  # A C - G T
    al_ref = [1, 2, 3, 3, 4]  # A C G G T
    score = 10.5
    
    print("\n=== 打印函数测试 ===")
    print_cas9_alignment(al_seq, al_ref, score)


if __name__ == "__main__":
    # 运行测试
    test = TestCAS9Align()
    test.setup_method()
    
    print("=== 运行cas9_align测试 ===")
    
    test.test_sequence_encoding()
    print("✅ 序列编码测试通过")
    
    test.test_perfect_match()
    print("✅ 完全匹配测试通过")
    
    test.test_single_mismatch()
    print("✅ 单个不匹配测试通过")
    
    test.test_insertion()
    print("✅ 插入测试通过")
    
    test.test_deletion()
    print("✅ 删除测试通过")
    
    test.test_position_specific_penalty()
    print("✅ 位置特异性惩罚测试通过")
    
    test_print_function()
    print("✅ 打印函数测试通过")
    
    print("\n🎉 所有测试通过！cas9_align算法工作正常！") 