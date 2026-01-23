"""
等位基因调用模块测试

测试CARLIN等位基因调用的各种功能和算法。
"""

import unittest
import numpy as np
from darlinpy.alignment import AlignedSEQ, AlignedSEQMotif, CARLINAligner, create_default_aligner
from darlinpy.calling import AlleleCaller, AlleleCallResult
from darlinpy.config import get_original_carlin_config
from typing import Optional
from darlinpy.config.amplicon_configs import load_carlin_config_by_locus


class TestAlleleData(unittest.TestCase):
    """测试等位基因数据结构"""
    
    def setUp(self):
        """设置测试数据"""
        # 创建简单的AlignedSEQ对象用于测试
        seq_segments = ["ATCG", "A-CG", "ATCGT"]
        ref_segments = ["ATCG", "ATCG", "ATC-G"]
        self.test_allele = AlignedSEQ(seq_segments, ref_segments)
        
    def test_allele_call_result_creation(self):
        """测试AlleleCallResult创建"""
        result = AlleleCallResult(
            allele=self.test_allele,
            constituents=[0, 1, 2],
            weight_contribution=[1.0, 2.0, 1.5],
            confidence=0.8,
            calling_method='coarse_grain'
        )
        
        self.assertIsNotNone(result.allele)
        self.assertEqual(result.total_weight, 4.5)
        self.assertEqual(result.num_supporting_sequences, 3)
        self.assertTrue(result.is_callable())
        self.assertEqual(result.event_structure, 'NDI')
        
    def test_allele_call_result_validation(self):
        """测试AlleleCallResult数据验证"""
        with self.assertRaises(ValueError):
            AlleleCallResult(
                allele=self.test_allele,
                constituents=[0, 1],  # 长度不匹配
                weight_contribution=[1.0, 2.0, 1.5],  # 长度不匹配
            )
            
    def test_allele_call_result_properties(self):
        """测试AlleleCallResult属性"""
        result = AlleleCallResult(
            allele=self.test_allele,
            constituents=[0, 1],
            weight_contribution=[1.0, 2.0],
            confidence=0.8,
            calling_method='coarse_grain'
        )
        
        self.assertEqual(result.total_weight, 3.0)
        self.assertEqual(result.num_supporting_sequences, 2)
        self.assertTrue(result.is_callable())
        self.assertEqual(result.event_structure, 'NDI')
        self.assertEqual(result.sequence, 'ATCGACGATCGT')


class TestAlleleCaller(unittest.TestCase):
    """测试等位基因调用算法"""
    
    def setUp(self):
        """设置测试环境"""
        # 使用Col1a1配置（默认）
        self.caller = AlleleCaller(locus="Col1a1")
        self.aligner = create_default_aligner()
        
        # 创建测试序列
        self.reference_seq = self.caller.amplicon_config.full_sequence
        
        # 创建一些测试序列（无突变、单个突变、多个突变）
        self.test_sequences = [
            self.reference_seq,  # 完美匹配
            self.reference_seq[:100] + 'A' + self.reference_seq[101:],  # 单个替换
            self.reference_seq[:100] + self.reference_seq[101:],  # 单个删除
            self.reference_seq[:100] + 'ATCG' + self.reference_seq[100:],  # 插入
        ]
        
    def test_unique_by_frequency(self):
        """测试按频率排序的唯一值获取"""
        values = ['A', 'B', 'A', 'C', 'B', 'A']
        weights = [1.0, 2.0, 1.0, 0.5, 2.0, 1.0]
        
        unique_vals, backtrack, unique_indices = self.caller.unique_by_frequency(values, weights)
        
        # 检查结果
        self.assertEqual(len(unique_vals), 3)  # A, B, C
        # B应该是第一个（权重=4.0），A第二个（权重=3.0），C第三个（权重=0.5）
        self.assertEqual(unique_vals[0], 'B')
        self.assertEqual(unique_vals[1], 'A')
        self.assertEqual(unique_vals[2], 'C')
        
    def test_call_alleles_exact_simple(self):
        """测试精确等位基因调用 - 简单情况"""
        # 创建aligned sequences
        aligned_seqs = []
        for seq in self.test_sequences[:2]:  # 只使用前两个序列
            result = self.aligner.align_sequence(seq, sanitize=True)
            if result['sanitized']:
                aligned_seqs.append(result['aligned_seq_obj'])
            else:
                # 如果标准化失败，手动创建简单的AlignedSEQ
                aligned_seqs.append(AlignedSEQ([seq], [self.reference_seq]))
                
        weights = [3.0, 1.0]  # 第一个序列权重更高
        
        result = self.caller.call_alleles_exact(aligned_seqs, weights, dominant_only=True)
        
        self.assertTrue(result.is_callable())
        self.assertEqual(result.calling_method, 'exact')
        self.assertGreater(result.confidence, 0.5)  # 应该是主导的
        
    def test_call_alleles_coarse_grain_simple(self):
        """测试粗粒度等位基因调用 - 简单情况"""
        # 创建aligned sequences
        aligned_seqs = []
        for seq in self.test_sequences:
            result = self.aligner.align_sequence(seq, sanitize=True)
            if result['sanitized']:
                aligned_seqs.append(result['aligned_seq_obj'])
            else:
                # 如果标准化失败，手动创建简单的AlignedSEQ
                aligned_seqs.append(AlignedSEQ([seq], [self.reference_seq]))
                
        weights = [5.0, 1.0, 1.0, 1.0]  # 第一个（参考）序列权重更高
        
        result = self.caller.call_alleles_coarse_grain(aligned_seqs, weights, dominant_only=True)
        
        self.assertTrue(result.is_callable())
        self.assertEqual(result.calling_method, 'coarse_grain')
        self.assertGreater(result.confidence, 0.5)
        
    def test_call_alleles_no_dominant(self):
        """测试没有主导等位基因的情况"""
        # 直接创建三个不同的AlignedSEQ对象，确保没有主导等位基因
        seq_segments1 = ["ATCG"]
        ref_segments1 = ["ATCG"]
        seq_segments2 = ["TTTT"]
        ref_segments2 = ["ATCG"]
        seq_segments3 = ["GGGG"]
        ref_segments3 = ["ATCG"]
        
        test_seqs = [
            AlignedSEQ(seq_segments1, ref_segments1),
            AlignedSEQ(seq_segments2, ref_segments2),
            AlignedSEQ(seq_segments3, ref_segments3)
        ]
        test_weights = [0.2, 0.2, 0.2]  # 权重相等，每个序列占0.2/0.6 = 0.33 < 0.5
        
        test_result = self.caller.call_alleles_exact(test_seqs, test_weights, dominant_only=True)
        
        # 应该调用失败，因为没有主导等位基因
        self.assertFalse(test_result.is_callable())
        self.assertLess(test_result.confidence, self.caller.dominant_threshold)
        
    def test_consensus_sequence_building(self):
        """测试共识序列构建"""
        # 创建几个相似的AlignedSEQ对象
        seq_segments = ["ATCG"]
        ref_segments = ["ATCG"]
        
        seqs = [
            AlignedSEQ(seq_segments, ref_segments), 
            AlignedSEQ(seq_segments, ref_segments), 
            AlignedSEQ(seq_segments, ref_segments)
        ]
        weights = [2.0, 1.0, 1.0]
        
        # 测试内部方法
        consensus = self.caller._build_consensus_sequence(seqs, weights)
        self.assertIsInstance(consensus, AlignedSEQ)


class TestAlleleCallingEdgeCases(unittest.TestCase):
    """测试等位基因调用边界情况"""
    
    def setUp(self):
        """设置测试环境"""
        self.caller = AlleleCaller(locus="Col1a1")
        
    def test_empty_sequence_handling(self):
        """测试空序列处理"""
        # 测试空序列列表
        result = self.caller.call_alleles_exact([], [], dominant_only=True)
        self.assertFalse(result.is_callable())
        
        # 测试包含None的序列列表
        seqs = [None, AlignedSEQ(["ATCG"], ["ATCG"]), None]
        weights = [1.0, 2.0, 1.0]
        
        result = self.caller.call_alleles_exact(seqs, weights, dominant_only=True)
        self.assertTrue(result.is_callable())
        self.assertEqual(len(result.constituents), 1)  # 只有一个有效序列
        
    def test_threshold_boundary_conditions(self):
        """测试阈值边界条件"""
        seqs = [
            AlignedSEQ(["ATCG"], ["ATCG"]),
            AlignedSEQ(["TTCG"], ["ATCG"])
        ]
        
        # 测试正好达到阈值的情况 (0.5)
        weights = [5.0, 5.0]  # 5.0/10.0 = 0.5，正好等于阈值
        result = self.caller.call_alleles_exact(seqs, weights, dominant_only=True)
        self.assertTrue(result.is_callable())  # 应该通过，因为是 >= 阈值

class TestAdvancedAlleleCalling(unittest.TestCase):
    """高级等位基因调用测试"""
    
    def setUp(self):
        """设置测试环境"""
        self.amplicon_config = get_original_carlin_config()
        # 使用关键字参数明确指定amplicon_config
        self.caller = AlleleCaller(amplicon_config=self.amplicon_config, dominant_threshold=0.6)  # 更高的阈值
        self.aligner = create_default_aligner()
        
    def test_complex_weight_scenarios(self):
        """测试复杂权重场景"""
        # 创建多个不同的序列
        seqs = [
            AlignedSEQ(["ATCG"], ["ATCG"]),  # 无变化
            AlignedSEQ(["A-CG"], ["ATCG"]),  # 删除
            AlignedSEQ(["ATCG"], ["ATCG"]),  # 无变化（重复）
            AlignedSEQ(["TTCG"], ["ATCG"]),  # 替换
        ]
        
        # 权重分布：第一种序列应该主导
        weights = [3.0, 1.0, 2.0, 0.5]  # 总共6.5，第一种序列类型5.0/6.5 ≈ 0.77 > 0.6
        
        result = self.caller.call_alleles_exact(seqs, weights, dominant_only=True)
        
        self.assertTrue(result.is_callable())
        self.assertGreater(result.confidence, 0.6)
        self.assertEqual(len(result.constituents), 2)  # 两个"ATCG"序列
        
    def test_coarse_grain_vs_exact_consistency(self):
        """测试粗粒度调用与精确调用的一致性"""
        # 创建具有相同事件结构但不同具体序列的数据
        seqs = [
            AlignedSEQ(["ATCG"], ["ATCG"]),  # N事件
            AlignedSEQ(["ATCG"], ["ATCG"]),  # N事件（相同）
            AlignedSEQ(["ATGG"], ["ATCG"]),  # M事件（替换）
        ]
        weights = [2.0, 1.0, 0.5]
        
        exact_result = self.caller.call_alleles_exact(seqs, weights, dominant_only=False)
        coarse_result = self.caller.call_alleles_coarse_grain(seqs, weights, dominant_only=False)
        
        # 两种方法都应该成功调用
        self.assertTrue(exact_result.is_callable())
        self.assertTrue(coarse_result.is_callable())
        
        # 检查事件结构一致性
        self.assertEqual(exact_result.event_structure, 'N')
        self.assertEqual(coarse_result.event_structure, 'N')
        
    def test_empty_and_null_handling(self):
        """测试空输入和空序列处理"""
        # 测试空序列列表
        result = self.caller.call_alleles_exact([], [], dominant_only=True)
        self.assertFalse(result.is_callable())
        
        # 测试包含None的序列列表
        seqs = [None, AlignedSEQ(["ATCG"], ["ATCG"]), None]
        weights = [1.0, 2.0, 1.0]
        
        result = self.caller.call_alleles_exact(seqs, weights, dominant_only=True)
        self.assertTrue(result.is_callable())
        self.assertEqual(len(result.constituents), 1)  # 只有一个有效序列
        
    def test_threshold_boundary_conditions(self):
        """测试阈值边界条件"""
        seqs = [
            AlignedSEQ(["ATCG"], ["ATCG"]),
            AlignedSEQ(["TTCG"], ["ATCG"])
        ]
        
        # 测试正好达到阈值的情况 (0.6)
        weights = [6.0, 4.0]  # 6.0/10.0 = 0.6，正好等于阈值
        result = self.caller.call_alleles_exact(seqs, weights, dominant_only=True)
        self.assertTrue(result.is_callable())  # 应该通过，因为是 >= 阈值
        
        # 测试刚好低于阈值的情况
        weights = [5.9, 4.1]  # 5.9/10.0 = 0.59 < 0.6
        result = self.caller.call_alleles_exact(seqs, weights, dominant_only=True)
        self.assertFalse(result.is_callable())
        
    def test_complex_weight_scenarios(self):
        """测试复杂权重场景"""
        # 创建多个不同的序列
        seqs = [
            AlignedSEQ(["ATCG"], ["ATCG"]),  # 无变化
            AlignedSEQ(["A-CG"], ["ATCG"]),  # 删除
            AlignedSEQ(["ATCG"], ["ATCG"]),  # 无变化（重复）
            AlignedSEQ(["TTCG"], ["ATCG"]),  # 替换
        ]
        
        # 权重分布：第一种序列应该主导
        weights = [3.0, 1.0, 2.0, 0.5]  # 总共6.5，第一种序列类型5.0/6.5 ≈ 0.77 > 0.5
        
        result = self.caller.call_alleles_exact(seqs, weights, dominant_only=True)
        
        self.assertTrue(result.is_callable())
        self.assertGreater(result.confidence, 0.5)
        self.assertEqual(len(result.constituents), 2)  # 两个"ATCG"序列
        
    def test_coarse_grain_vs_exact_consistency(self):
        """测试粗粒度调用与精确调用的一致性"""
        # 创建具有相同事件结构但不同具体序列的数据
        seqs = [
            AlignedSEQ(["ATCG"], ["ATCG"]),  # N事件
            AlignedSEQ(["ATCG"], ["ATCG"]),  # N事件（相同）
            AlignedSEQ(["ATGG"], ["ATCG"]),  # M事件（替换）
        ]
        weights = [2.0, 1.0, 0.5]
        
        exact_result = self.caller.call_alleles_exact(seqs, weights, dominant_only=False)
        coarse_result = self.caller.call_alleles_coarse_grain(seqs, weights, dominant_only=False)
        
        # 两种方法都应该成功调用
        self.assertTrue(exact_result.is_callable())
        self.assertTrue(coarse_result.is_callable())
        
        # 检查事件结构一致性
        self.assertEqual(exact_result.event_structure, 'N')
        self.assertEqual(coarse_result.event_structure, 'N')


def run_tests():
    """运行所有测试"""
    unittest.main(verbosity=2)


if __name__ == '__main__':
    run_tests() 