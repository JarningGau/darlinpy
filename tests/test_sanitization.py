#!/usr/bin/env python3
"""
测试序列标准化功能
"""

import sys
import os

# 添加路径以便导入模块
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from darlin.alignment import AlignedSEQ, AlignedSEQMotif, SequenceSanitizer, desemble_sequence, create_default_aligner
from darlin.config import get_original_carlin_config


class TestAlignedSEQMotif:
    
    def test_motif_classification(self):
        """测试motif事件分类"""
        # 测试无变化 (N)
        motif = AlignedSEQMotif("ACGT", "ACGT")
        assert motif.event == 'N'
        
        # 测试错配 (M)
        motif = AlignedSEQMotif("ACGT", "ACCT")
        assert motif.event == 'M'
        
        # 测试删除 (D)
        motif = AlignedSEQMotif("AC-T", "ACGT")
        assert motif.event == 'D'
        
        # 测试插入 (I)
        motif = AlignedSEQMotif("ACGGT", "AC--T")
        assert motif.event == 'I'
        
        # 测试空 (E)
        motif = AlignedSEQMotif("----", "ACGT")
        assert motif.event == 'E'
        
        print("✅ Motif事件分类测试通过")


class TestAlignedSEQ:
    
    def test_aligned_seq_construction(self):
        """测试AlignedSEQ构建"""
        seq_segments = ["ACGT", "GG-T", "AAAA"]
        ref_segments = ["ACGT", "GGTT", "AAAA"]
        
        aligned_seq = AlignedSEQ(seq_segments, ref_segments)
        
        assert aligned_seq.get_seq() == "ACGTGG-TAAAA"
        assert aligned_seq.get_ref() == "ACGTGGTTAAAA"
        assert aligned_seq.get_event_structure() == ['N', 'D', 'N']
        
        print("✅ AlignedSEQ构建测试通过")
    
    def test_aligned_seq_copy(self):
        """测试AlignedSEQ拷贝"""
        seq_segments = ["ACGT", "GGGT"]
        ref_segments = ["ACCT", "GGGT"]
        
        original = AlignedSEQ(seq_segments, ref_segments)
        copied = original.copy()
        
        assert copied.get_seq() == original.get_seq()
        assert copied.get_ref() == original.get_ref()
        assert copied.get_event_structure() == original.get_event_structure()
        
        # 验证是深拷贝
        assert copied is not original
        assert copied.motifs is not original.motifs
        
        print("✅ AlignedSEQ拷贝测试通过")


class TestSequenceSanitizer:
    
    def test_sanitize_prefix_postfix(self):
        """测试prefix/postfix标准化"""
        # 测试头部插入修剪
        seq_segments = ["--ACGT", "GGTT", "AAAA"]
        ref_segments = ["--ACGT", "GGTT", "AAAA"] 
        
        # 模拟头部插入
        seq_segments[0] = "XXACGT"
        ref_segments[0] = "--ACGT"
        
        aligned_seq = AlignedSEQ(seq_segments, ref_segments)
        sanitized = SequenceSanitizer.sanitize_prefix_postfix(aligned_seq)
        
        # 应该修剪掉开头的XX
        assert sanitized.motifs[0].seq == "ACGT"
        assert sanitized.motifs[0].ref == "ACGT"
        
        print("✅ Prefix/Postfix标准化测试通过")
        print(f"   - 原始: {aligned_seq.get_seq()}")
        print(f"   - 标准化: {sanitized.get_seq()}")
    
    def test_sanitize_conserved_regions(self):
        """测试保守区域标准化"""
        # 创建带错配的序列，其中一些是cutsite，一些不是
        seq_segments = ["ACGT", "CCCC", "TTTT", "AAAA"]  # 第2和4个motif有错配
        ref_segments = ["ACGT", "GGGG", "TTTT", "CCCC"]
        
        aligned_seq = AlignedSEQ(seq_segments, ref_segments)
        
        # 假设第2个motif(index=1)是cutsite，第4个motif(index=3)不是
        cutsite_indices = [1]  # 只有index=1是cutsite
        
        sanitized = SequenceSanitizer.sanitize_conserved_regions(aligned_seq, cutsite_indices)
        
        # cutsite motif (index=1) 应该保持不变
        assert sanitized.motifs[1].seq == "CCCC"  # 保持原样
        assert sanitized.motifs[1].ref == "GGGG"
        
        # 非cutsite motif (index=3) 应该被修正为参考序列
        assert sanitized.motifs[3].seq == "CCCC"  # 修正为参考序列
        assert sanitized.motifs[3].ref == "CCCC"
        
        print("✅ 保守区域标准化测试通过")
        print(f"   - 原始事件: {aligned_seq.get_event_structure()}")
        print(f"   - 标准化事件: {sanitized.get_event_structure()}")
    
    def test_batch_sanitization(self):
        """测试批量标准化"""
        # 创建多个序列
        seqs = []
        for i in range(3):
            seq_segments = [f"ACG{i}", "GGTT"]
            ref_segments = ["ACGT", "GGTT"]
            seqs.append(AlignedSEQ(seq_segments, ref_segments))
        
        # 批量prefix/postfix标准化
        sanitized_seqs = SequenceSanitizer.sanitize_prefix_postfix(seqs)
        
        assert len(sanitized_seqs) == 3
        assert all(isinstance(seq, AlignedSEQ) for seq in sanitized_seqs)
        
        print("✅ 批量标准化测试通过")


class TestDesembleSequence:
    
    def test_desemble_sequence(self):
        """测试序列分解"""
        aligned_query = "ACGTGGTTAAAA"
        aligned_ref = "ACGTGGTTAAAA"
        motif_boundaries = [(0, 4), (4, 8), (8, 12)]
        
        aligned_seq = desemble_sequence(aligned_query, aligned_ref, motif_boundaries)
        
        assert len(aligned_seq.motifs) == 3
        assert aligned_seq.motifs[0].seq == "ACGT"
        assert aligned_seq.motifs[1].seq == "GGTT"
        assert aligned_seq.motifs[2].seq == "AAAA"
        
        # 验证重建的序列
        assert aligned_seq.get_seq() == aligned_query
        assert aligned_seq.get_ref() == aligned_ref
        
        print("✅ 序列分解测试通过")


class TestIntegratedSanitization:
    
    def test_carlin_aligner_with_sanitization(self):
        """测试集成的CARLIN比对器标准化功能"""
        config = get_original_carlin_config()
        perfect_seq = config.get_reference_sequence()
        
        # 创建一个带有小变化的序列
        mutated_seq = perfect_seq[:50] + "A" + perfect_seq[51:]  # 单个替换
        
        aligner = create_default_aligner()
        
        # 不使用标准化
        result_no_sanitize = aligner.align_sequence(mutated_seq, verbose=False, sanitize=False)
        
        # 使用标准化
        result_with_sanitize = aligner.align_sequence(mutated_seq, verbose=False, sanitize=True)
        
        # 验证标准化标志
        assert not result_no_sanitize['sanitized']
        assert result_with_sanitize['sanitized']
        
        # 验证标准化后的对象
        assert result_no_sanitize['aligned_seq_obj'] is None
        assert result_with_sanitize['aligned_seq_obj'] is not None
        assert isinstance(result_with_sanitize['aligned_seq_obj'], AlignedSEQ)
        
        print("✅ 集成标准化测试通过")
        print(f"   - 未标准化得分: {result_no_sanitize['alignment_score']:.1f}")
        print(f"   - 标准化得分: {result_with_sanitize['alignment_score']:.1f}")
        print(f"   - 标准化对象motifs数: {len(result_with_sanitize['aligned_seq_obj'].motifs)}")
    
    def test_sanitization_with_real_carlin_mutations(self):
        """测试真实CARLIN突变的标准化"""
        config = get_original_carlin_config()
        ref_seq = config.get_reference_sequence()
        
        # 在cutsite区域引入真实编辑
        cutsite_start, cutsite_end = config.positions['cutsites'][0]
        
        # 在cutsite中间插入
        edited_seq = ref_seq[:cutsite_start+3] + "AAA" + ref_seq[cutsite_start+3:]
        
        aligner = create_default_aligner()
        result = aligner.align_sequence(edited_seq, verbose=True, sanitize=True)
        
        # 验证标准化成功
        assert result['sanitized']
        
        # 获取标准化后的对象
        sanitized_obj = result['aligned_seq_obj']
        assert sanitized_obj is not None
        
        # 分析事件结构
        events = sanitized_obj.get_event_structure()
        print(f"   - 事件结构: {''.join(events)}")
        
        # 应该有插入事件
        assert 'I' in events
        
        print("✅ 真实CARLIN突变标准化测试通过")


def main():
    """主测试函数"""
    print("=== 序列标准化测试 ===")
    print()
    
    # 基础测试
    motif_test = TestAlignedSEQMotif()
    motif_test.test_motif_classification()
    
    seq_test = TestAlignedSEQ()
    seq_test.test_aligned_seq_construction()
    seq_test.test_aligned_seq_copy()
    
    sanitizer_test = TestSequenceSanitizer()
    sanitizer_test.test_sanitize_prefix_postfix()
    sanitizer_test.test_sanitize_conserved_regions()
    sanitizer_test.test_batch_sanitization()
    
    desemble_test = TestDesembleSequence()
    desemble_test.test_desemble_sequence()
    
    print()
    
    # 集成测试
    integrated_test = TestIntegratedSanitization()
    integrated_test.test_carlin_aligner_with_sanitization()
    integrated_test.test_sanitization_with_real_carlin_mutations()
    
    print()
    print("🎉 所有序列标准化测试通过！")
    print()
    print("✅ 验证的功能:")
    print("   - AlignedSEQMotif事件分类")
    print("   - AlignedSEQ序列重建")
    print("   - Prefix/Postfix标准化")
    print("   - 保守区域标准化") 
    print("   - 序列分解和边界计算")
    print("   - CARLIN比对器集成")
    print("   - 真实突变序列处理")


if __name__ == "__main__":
    main() 