#!/usr/bin/env python3
"""
CARLIN突变注释模块测试
"""

import pytest
import numpy as np
from darlin.mutations import (
    Mutation, 
    MutationType, 
    MutationIdentifier, 
    annotate_mutations
)
from darlin.alignment.aligned_seq import AlignedSEQ, AlignedSEQMotif


class TestMutation:
    """测试Mutation类"""
    
    def test_mutation_creation(self):
        """测试突变对象创建"""
        mutation = Mutation(
            type=MutationType.SUBSTITUTION,
            loc_start=10,
            loc_end=10,
            seq_old="A",
            seq_new="T",
            motif_index=0
        )
        
        assert mutation.type == MutationType.SUBSTITUTION
        assert mutation.loc_start == 10
        assert mutation.loc_end == 10
        assert mutation.seq_old == "A"
        assert mutation.seq_new == "T"
        assert mutation.motif_index == 0
        assert mutation.confidence == 1.0
    
    def test_mutation_validation(self):
        """测试突变对象验证"""
        with pytest.raises(ValueError):
            # 位置不能小于1
            Mutation(
                type=MutationType.SUBSTITUTION,
                loc_start=0,
                loc_end=5,
                seq_old="A",
                seq_new="T"
            )
        
        with pytest.raises(ValueError):
            # 结束位置不能小于起始位置
            Mutation(
                type=MutationType.SUBSTITUTION,
                loc_start=10,
                loc_end=5,
                seq_old="A",
                seq_new="T"
            )
        
        with pytest.raises(ValueError):
            # 置信度必须在0-1之间
            Mutation(
                type=MutationType.SUBSTITUTION,
                loc_start=10,
                loc_end=10,
                seq_old="A",
                seq_new="T",
                confidence=1.5
            )
    
    def test_mutation_properties(self):
        """测试突变属性"""
        # 测试删除
        deletion = Mutation(
            type=MutationType.DELETION,
            loc_start=10,
            loc_end=12,
            seq_old="ATG",
            seq_new=""
        )
        assert deletion.length_change == -3
        assert deletion.is_indel == True
        assert deletion.affected_length == 3
        
        # 测试插入
        insertion = Mutation(
            type=MutationType.INSERTION,
            loc_start=10,
            loc_end=10,
            seq_old="",
            seq_new="GCC"
        )
        assert insertion.length_change == 3
        assert insertion.is_indel == True
        assert insertion.affected_length == 1
        
        # 测试替换
        substitution = Mutation(
            type=MutationType.SUBSTITUTION,
            loc_start=10,
            loc_end=10,
            seq_old="A",
            seq_new="T"
        )
        assert substitution.length_change == 0
        assert substitution.is_indel == False
        assert substitution.affected_length == 1
    
    def test_hgvs_annotation(self):
        """测试HGVS格式注释"""
        # 单碱基替换
        sub = Mutation(
            type=MutationType.SUBSTITUTION,
            loc_start=100,
            loc_end=100,
            seq_old="A",
            seq_new="T"
        )
        assert sub.to_hgvs() == "CARLIN:g.100A>T"
        
        # 多碱基替换
        multi_sub = Mutation(
            type=MutationType.SUBSTITUTION,
            loc_start=100,
            loc_end=102,
            seq_old="ATG",
            seq_new="GCC"
        )
        assert multi_sub.to_hgvs() == "CARLIN:g.100_102ATG>GCC"
        
        # 单碱基删除
        del_single = Mutation(
            type=MutationType.DELETION,
            loc_start=100,
            loc_end=100,
            seq_old="A",
            seq_new=""
        )
        assert del_single.to_hgvs() == "CARLIN:g.100del"
        
        # 多碱基删除
        del_multi = Mutation(
            type=MutationType.DELETION,
            loc_start=100,
            loc_end=102,
            seq_old="ATG",
            seq_new=""
        )
        assert del_multi.to_hgvs() == "CARLIN:g.100_102del"
        
        # 插入
        insertion = Mutation(
            type=MutationType.INSERTION,
            loc_start=100,
            loc_end=100,
            seq_old="",
            seq_new="GCC"
        )
        assert insertion.to_hgvs() == "CARLIN:g.100_101insGCC"
        
        # 复合突变
        complex_mut = Mutation(
            type=MutationType.COMPLEX,
            loc_start=100,
            loc_end=102,
            seq_old="ATG",
            seq_new="GCCAA"
        )
        assert complex_mut.to_hgvs() == "CARLIN:g.100_102delinsGCCAA"
    
    def test_string_representation(self):
        """测试字符串表示"""
        mutation = Mutation(
            type=MutationType.SUBSTITUTION,
            loc_start=100,
            loc_end=100,
            seq_old="A",
            seq_new="T",
            motif_index=3
        )
        
        # 测试__str__
        assert str(mutation) == "CARLIN:g.100A>T"
        
        # 测试__repr__
        repr_str = repr(mutation)
        assert "Mutation(type=M" in repr_str
        assert "pos=100-100" in repr_str
        assert "'A'->'T'" in repr_str
        assert "motif=3" in repr_str


class TestMutationIdentifier:
    """测试MutationIdentifier类"""
    
    def test_identify_sequence_events_simple(self):
        """测试简单的序列事件识别"""
        # 创建一个简单的比对序列：替换
        aligned_seq = AlignedSEQ(
            seq_segments=["ACGT", "TGCA"],
            ref_segments=["ACGT", "AGCA"]  # 第二个motif有T->A替换
        )
        
        identifier = MutationIdentifier()
        mutations = identifier.identify_sequence_events(aligned_seq)
        
        # 应该识别出一个替换突变
        assert len(mutations) == 1
        mutation = mutations[0]
        assert mutation.type == MutationType.SUBSTITUTION
        assert mutation.loc_start == 5  # 第二个motif起始位置
        assert mutation.loc_end == 8    # 第二个motif结束位置
        assert mutation.seq_old == "AGCA"
        assert mutation.seq_new == "TGCA"
        assert mutation.motif_index == 1
    
    def test_identify_deletion(self):
        """测试删除事件识别"""
        # 创建包含删除的比对序列
        aligned_seq = AlignedSEQ(
            seq_segments=["ACGT", "----"],  # 第二个motif完全删除
            ref_segments=["ACGT", "TTGG"]
        )
        
        identifier = MutationIdentifier()
        mutations = identifier.identify_sequence_events(aligned_seq)
        
        assert len(mutations) == 1
        mutation = mutations[0]
        assert mutation.type == MutationType.DELETION
        assert mutation.loc_start == 5
        assert mutation.loc_end == 8
        assert mutation.seq_old == "TTGG"
        assert mutation.seq_new == ""
        assert mutation.motif_index == 1
    
    def test_identify_insertion(self):
        """测试插入事件识别"""
        # 创建包含插入的比对序列
        aligned_seq = AlignedSEQ(
            seq_segments=["ACGT", "TTGG"],
            ref_segments=["ACGT", "----"]  # 第二个motif是插入
        )
        
        identifier = MutationIdentifier()
        mutations = identifier.identify_sequence_events(aligned_seq)
        
        assert len(mutations) == 1
        mutation = mutations[0]
        assert mutation.type == MutationType.INSERTION
        assert mutation.loc_start == 5
        assert mutation.loc_end == 5  # 插入的结束位置(与起始位置相同)
        assert mutation.seq_old == ""
        assert mutation.seq_new == "TTGG"
        assert mutation.motif_index == 1
    
    def test_identify_complex_mutation(self):
        """测试复合突变识别"""
        # 创建包含复合突变的比对序列(长度变化大于3)
        aligned_seq = AlignedSEQ(
            seq_segments=["ACGT", "ATGGCCAA"],  # 8bp替换4bp
            ref_segments=["ACGT", "TTGG----"]
        )
        
        identifier = MutationIdentifier()
        mutations = identifier.identify_sequence_events(aligned_seq)
        
        assert len(mutations) == 1
        mutation = mutations[0]
        assert mutation.type == MutationType.COMPLEX
        assert mutation.seq_old == "TTGG"
        assert mutation.seq_new == "ATGGCCAA"
    
    def test_cas9_event_identification(self):
        """测试Cas9特异性事件识别"""
        # 创建在切割位点附近的突变
        aligned_seq = AlignedSEQ(
            seq_segments=["ACGT", "----"],  # 删除在位置5-8
            ref_segments=["ACGT", "TTGG"]
        )
        
        identifier = MutationIdentifier()
        cut_sites = [7]  # 切割位点在第7位
        
        cas9_mutations = identifier.identify_cas9_events(aligned_seq, cut_sites)
        
        assert len(cas9_mutations) == 1
        mutation = cas9_mutations[0]
        assert mutation.type == MutationType.DELETION
        # 置信度应该增加(靠近切割位点)
        assert mutation.confidence > 0.8
    
    def test_merge_adjacent_mutations(self):
        """测试相邻突变合并"""
        # 创建两个相邻的突变
        mut1 = Mutation(
            type=MutationType.DELETION,
            loc_start=10,
            loc_end=12,
            seq_old="ATG",
            seq_new="",
            confidence=0.9
        )
        
        mut2 = Mutation(
            type=MutationType.DELETION,
            loc_start=15,
            loc_end=17,
            seq_old="GCC",
            seq_new="",
            confidence=0.8
        )
        
        identifier = MutationIdentifier()
        merged = identifier.merge_adjacent_mutations([mut1, mut2], max_distance=3)
        
        # 应该合并为一个突变
        assert len(merged) == 1
        merged_mut = merged[0]
        assert merged_mut.loc_start == 10
        assert merged_mut.loc_end == 17
        assert abs(merged_mut.confidence - 0.85) < 1e-10  # 平均置信度(考虑浮点精度)
    
    def test_no_merge_distant_mutations(self):
        """测试距离较远的突变不合并"""
        mut1 = Mutation(
            type=MutationType.DELETION,
            loc_start=10,
            loc_end=12,
            seq_old="ATG",
            seq_new=""
        )
        
        mut2 = Mutation(
            type=MutationType.DELETION,
            loc_start=20,
            loc_end=22,
            seq_old="GCC",
            seq_new=""
        )
        
        identifier = MutationIdentifier()
        merged = identifier.merge_adjacent_mutations([mut1, mut2], max_distance=3)
        
        # 距离太远，不应该合并
        assert len(merged) == 2


class TestAnnotateMutations:
    """测试便利函数annotate_mutations"""
    
    def test_annotate_mutations_basic(self):
        """测试基础突变注释"""
        aligned_seq = AlignedSEQ(
            seq_segments=["ACGT", "TGCA"],
            ref_segments=["ACGT", "AGCA"]
        )
        
        mutations = annotate_mutations(aligned_seq, cas9_mode=False)
        
        assert len(mutations) == 1
        assert mutations[0].type == MutationType.SUBSTITUTION
    
    def test_annotate_mutations_cas9_mode(self):
        """测试Cas9模式突变注释"""
        aligned_seq = AlignedSEQ(
            seq_segments=["ACGT", "----"],
            ref_segments=["ACGT", "TTGG"]
        )
        
        mutations = annotate_mutations(
            aligned_seq, 
            cas9_mode=True, 
            cut_sites=[7],
            merge_adjacent=True
        )
        
        assert len(mutations) == 1
        assert mutations[0].type == MutationType.DELETION
        assert mutations[0].confidence > 0.8  # Cas9模式下置信度应该较高
    
    def test_annotate_mutations_no_mutations(self):
        """测试无突变序列"""
        aligned_seq = AlignedSEQ(
            seq_segments=["ACGT", "TTGG"],
            ref_segments=["ACGT", "TTGG"]  # 完全匹配
        )
        
        mutations = annotate_mutations(aligned_seq)
        
        assert len(mutations) == 0


if __name__ == "__main__":
    pytest.main([__file__]) 