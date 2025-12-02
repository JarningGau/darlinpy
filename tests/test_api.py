#!/usr/bin/env python3
"""
DARLIN主API模块测试
"""

import pytest
import time
from darlin.api import analyze_sequences, AnalysisResult


class TestAnalysisResult:
    """测试AnalysisResult数据类"""
    
    def test_analysis_result_creation(self):
        """测试AnalysisResult创建"""
        # 创建最小的有效结果
        result = AnalysisResult(
            called_alleles=[],
            mutations=[],
            alignment_scores=[85.5, 92.1, 78.3],
            summary_stats={'total_sequences': 3}
        )
        
        assert result.num_sequences == 3
        assert result.num_called_alleles == 0
        assert result.calling_success_rate == 0.0
        assert result.total_mutations == 0
        assert abs(result.average_alignment_score - 85.3) < 0.1
    
    def test_analysis_result_validation(self):
        """测试AnalysisResult数据验证"""
        with pytest.raises(ValueError):
            # called_alleles和mutations长度不一致
            AnalysisResult(
                called_alleles=[None],  # 1个元素
                mutations=[[], []],     # 2个元素
                alignment_scores=[85.5],
                summary_stats={}
            )
    
    def test_analysis_result_properties(self):
        """测试AnalysisResult属性方法"""
        from darlin.calling.allele_data import AlleleCallResult
        from darlin.mutations.mutation import Mutation, MutationType
        from darlin.alignment.aligned_seq import AlignedSEQ
        
        # 创建模拟数据
        aligned_seq = AlignedSEQ(
            seq_segments=["ACGT", "TTGG"],
            ref_segments=["ACGT", "TTGG"]
        )
        
        allele_result = AlleleCallResult(
            allele=aligned_seq,
            constituents=[0],
            weight_contribution=[1.0],
            confidence=0.95
        )
        
        mutation = Mutation(
            type=MutationType.SUBSTITUTION,
            loc_start=5,
            loc_end=5,
            seq_old="A",
            seq_new="T"
        )
        
        result = AnalysisResult(
            called_alleles=[allele_result],
            mutations=[[mutation]],
            alignment_scores=[85.5, 92.1],
            summary_stats={'total_sequences': 2}
        )
        
        assert result.num_sequences == 2
        assert result.num_called_alleles == 1
        assert result.calling_success_rate == 1.0
        assert result.total_mutations == 1
        
        # 测试突变摘要
        mut_summary = result.get_mutation_summary()
        assert mut_summary['M'] == 1
    
    def test_to_df_method(self):
        """测试to_df方法"""
        from darlin.calling.allele_data import AlleleCallResult
        from darlin.mutations.mutation import Mutation, MutationType
        from darlin.alignment.aligned_seq import AlignedSEQ
        
        # 创建模拟数据
        aligned_seq = AlignedSEQ(
            seq_segments=["ACGT", "TTGG"],
            ref_segments=["ACGT", "TTGG"]
        )
        
        allele_result = AlleleCallResult(
            allele=aligned_seq,
            constituents=[0],
            weight_contribution=[1.0],
            confidence=0.95
        )
        
        mutation = Mutation(
            type=MutationType.SUBSTITUTION,
            loc_start=5,
            loc_end=5,
            seq_old="A",
            seq_new="T"
        )
        
        result = AnalysisResult(
            called_alleles=[allele_result, None],  # 2个元素，第二个为None
            mutations=[[mutation], []],  # 2个元素，第二个为空列表
            alignment_scores=[85.5, 92.1],
            aligned_query=["ACGTACGT", "TTGGTTGG"],
            aligned_reference=["ACGTACGT", "TTGGTTGG"],
            valid_sequences=["ACGTACGT", "TTGGTTGG"],
            summary_stats={'total_sequences': 2}
        )
        
        # 测试to_df方法
        df = result.to_df()
        
        # 验证DataFrame结构
        assert len(df) == 2
        assert 'query' in df.columns
        assert 'query_len' in df.columns
        assert 'aligned_query' in df.columns
        assert 'aligned_ref' in df.columns
        assert 'scores' in df.columns
        assert 'mutations' in df.columns
        assert 'confidence' in df.columns
        
        # 验证数据内容
        assert df['query'].iloc[0] == "ACGTACGT"
        assert df['query_len'].iloc[0] == 8
        assert df['scores'].iloc[0] == 85.5
        assert '5A>T' in df['mutations'].iloc[0]  # HGVS格式的突变
        
        # 验证第二个序列（无突变）
        assert df['query'].iloc[1] == "TTGGTTGG"
        assert df['query_len'].iloc[1] == 8
        assert df['scores'].iloc[1] == 92.1
        assert df['mutations'].iloc[1] == []  # 无突变
        assert df['confidence'].iloc[1] == []  # 无置信度
    
    def test_print_summary(self, capsys):
        """测试打印摘要功能"""
        result = AnalysisResult(
            called_alleles=[],
            mutations=[],
            alignment_scores=[85.5],
            summary_stats={},
            config_used="Col1a1",
            method_used="coarse_grain"
        )
        
        result.print_summary()
        captured = capsys.readouterr()
        
        assert "CARLIN Sequence Analysis Results Summary" in captured.out
        assert "Col1a1" in captured.out
        assert "coarse_grain" in captured.out


class TestAnalyzeSequences:
    """测试analyze_sequences主函数"""
    
    def test_analyze_empty_sequences(self):
        """测试空序列列表"""
        with pytest.raises(ValueError, match="Sequence list cannot be empty"):
            analyze_sequences([])
    
    def test_analyze_invalid_config(self):
        """测试无效配置"""
        sequences = ["ACGTACGTACGTACGTACGT" * 10]  # 200bp序列
        
        with pytest.raises(ValueError, match="Unsupported locus"):
            analyze_sequences(sequences, config="InvalidConfig")
    
    def test_analyze_invalid_method(self):
        """测试无效方法"""
        sequences = ["ACGTACGTACGTACGTACGT" * 10]
        
        with pytest.raises(ValueError, match="Unsupported calling method"):
            analyze_sequences(sequences, method="invalid_method")
    
    def test_analyze_short_sequences(self):
        """测试过短的序列"""
        short_sequences = ["ACGT", "TTGG", "AAAA"]  # 都少于50bp
        
        with pytest.raises(ValueError, match="No valid sequences available for analysis"):
            analyze_sequences(short_sequences)
    
    def test_analyze_valid_sequences_basic(self):
        """测试基本的有效序列分析"""
        # 创建基于CARLIN结构的测试序列
        from darlin.config.amplicon_configs import load_carlin_config_by_locus
        carlin_config = load_carlin_config_by_locus('Col1a1')
        carlin_ref = carlin_config.get_full_reference_sequence()
        
        # 创建一些变异序列用于测试
        sequences = [
            carlin_ref[:200],  # 前200bp
            carlin_ref[50:250],  # 中间200bp  
            carlin_ref[-200:],  # 后200bp
        ]
        
        result = analyze_sequences(
            sequences, 
            config='Col1a1',
            method='coarse_grain',
            verbose=False
        )
        
        # 验证基本结果结构
        assert isinstance(result, AnalysisResult)
        assert len(result.called_alleles) == len(sequences)
        assert len(result.mutations) == len(sequences)
        assert len(result.alignment_scores) == len(sequences)
        assert result.config_used == 'Col1a1'
        assert result.method_used == 'coarse_grain'
        assert result.processing_time > 0
    
    def test_analyze_sequences_exact_method(self):
        """测试精确方法"""
        from darlin.config.amplicon_configs import load_carlin_config_by_locus
        carlin_config = load_carlin_config_by_locus('Col1a1')
        carlin_ref = carlin_config.get_full_reference_sequence()
        sequences = [carlin_ref[:200]]
        
        result = analyze_sequences(
            sequences,
            method='exact',
            verbose=False
        )
        
        assert result.method_used == 'exact'

    def test_analyze_sequences_expected_delins_allele(self):
        """测试特定等位基因的delins注释是否为23_265delinsG"""
        from darlin import analyze_sequences as analyze_sequences_public

        sequences = [
            "CGCCGGACTGCACGACAGTCGAGCGATGGGAGCT",
            "CGCCGGACTGCACGACAGTCGACGATGGAGTCGACACGACTCGCGCATACGATGGAGTCGACTACAGTCGCTACGACGATGGAGTCGCGAGCGCTATGAGCGACTATGGAGTCGATACGATACGCGCACGCTATGGAGTCGAGAGCGCGCTCGTCGACTATGGAGTCGCGACTGTACGCACACGCGATGGAGTCGATAGTATGCGTACACGCGATGGAGTCGAGTCGAGACGCTGACGATATGGAGTCGATACGTAGCACGCAGACGATGGGAGCT",
        ]

        results = analyze_sequences_public(
            sequences,
            config="Col1a1",
            method="exact",
            min_sequence_length=20,
            verbose=False,
        )
        df = results.to_df()
        assert df["mutations"].tolist() == ["23_265delinsG", []]
    
    def test_analyze_sequences_no_mutations(self):
        """测试不进行突变注释"""
        from darlin.config.amplicon_configs import load_carlin_config_by_locus
        carlin_config = load_carlin_config_by_locus('Col1a1')
        carlin_ref = carlin_config.get_full_reference_sequence()
        sequences = [carlin_ref[:200]]
        
        result = analyze_sequences(
            sequences,
            annotate_mutations_flag=False,
            verbose=False
        )
        
        # 所有突变列表应该为空
        assert all(len(muts) == 0 for muts in result.mutations)
    
    def test_analyze_sequences_verbose(self, capsys):
        """测试详细输出模式"""
        from darlin.config.amplicon_configs import load_carlin_config_by_locus
        carlin_config = load_carlin_config_by_locus('Col1a1')
        carlin_ref = carlin_config.get_full_reference_sequence()
        sequences = [carlin_ref[:200]]
        
        analyze_sequences(sequences, verbose=True)
        
        captured = capsys.readouterr()
        assert "Starting analysis" in captured.out
        assert "Alignment completed" in captured.out
        assert "Allele calling completed" in captured.out
    
    def test_analyze_different_configs(self):
        """测试不同配置"""
        from darlin.config.amplicon_configs import load_carlin_config_by_locus
        carlin_config = load_carlin_config_by_locus('Col1a1')
        carlin_ref = carlin_config.get_full_reference_sequence()
        sequences = [carlin_ref[:200]]
        
        # 测试不同的配置名称格式
        configs_to_test = ['Col1a1', 'Rosa']
        
        for config in configs_to_test:
            result = analyze_sequences(sequences, config=config, verbose=False)
            assert result.config_used == config
    
    def test_analyze_custom_parameters(self):
        """测试自定义参数"""
        from darlin.config.amplicon_configs import load_carlin_config_by_locus
        carlin_config = load_carlin_config_by_locus('Col1a1')
        carlin_ref = carlin_config.get_full_reference_sequence()
        sequences = [carlin_ref[:200]]
        
        result = analyze_sequences(
            sequences,
            dominant_threshold=0.7,
            merge_adjacent_mutations=False,
            verbose=False
        )
        
        assert isinstance(result, AnalysisResult)


class TestConfigLoading:
    """测试配置加载功能"""
    
    def test_load_config_by_locus(self):
        """测试通过locus名称加载配置"""
        from darlin.config.amplicon_configs import load_carlin_config_by_locus
        
        config = load_carlin_config_by_locus('Col1a1')
        assert config is not None
        
        config = load_carlin_config_by_locus('Rosa')
        assert config is not None
        
        config = load_carlin_config_by_locus('Tigre')
        assert config is not None


class TestStatistics:
    """测试统计功能"""
    
    def test_generate_summary_stats(self):
        """测试统计信息生成"""
        from darlin.api import _generate_summary_stats
        from darlin.alignment.aligned_seq import AlignedSEQ
        from darlin.calling.allele_data import AlleleCallResult
        from darlin.mutations.mutation import Mutation, MutationType
        
        # 创建模拟数据
        sequences = ["ACGTACGT" * 25, "TGCATGCA" * 25]  # 200bp序列
        
        # 模拟比对结果
        alignment_results = [
            {'alignment_score': 85.5},
            {'alignment_score': 92.1}
        ]
        
        # 模拟等位基因调用结果
        aligned_seq = AlignedSEQ(
            seq_segments=["ACGT", "TTGG"],
            ref_segments=["ACGT", "TTGG"]
        )
        
        allele_results = [
            AlleleCallResult(
                allele=aligned_seq,
                constituents=[0],
                weight_contribution=[1.0]
            )
        ]
        
        # 模拟突变
        mutations_list = [[
            Mutation(
                type=MutationType.SUBSTITUTION,
                loc_start=5,
                loc_end=5,
                seq_old="A",
                seq_new="T"
            )
        ]]
        
        stats = _generate_summary_stats(
            sequences, alignment_results, allele_results, mutations_list
        )
        
        # 验证统计结果
        assert stats['total_sequences'] == 2
        assert stats['avg_sequence_length'] == 200.0
        assert stats['called_alleles_count'] == 1
        assert stats['avg_alignment_score'] == (85.5 + 92.1) / 2
        assert stats['total_mutations'] == 1
        assert stats['mutation_type_distribution']['M'] == 1


class TestErrorHandling:
    """测试错误处理"""
    
    def test_analyze_sequences_runtime_error(self):
        """测试运行时错误处理"""
        # 这个测试可能需要更复杂的设置来模拟运行时错误
        # 目前我们测试基本的错误传播
        sequences = ["ACGT"]  # 太短的序列
        
        with pytest.raises(ValueError):
            analyze_sequences(sequences)


class TestComplexMutationCases:
    """测试复杂突变案例（来自issues.tmp）"""
    
    def test_case_1_delins_identification(self):
        """测试案例1：23_76delinsGT的正确识别"""
        # Query sequence from issues.tmp line 102 (removing dashes)
        query = "CGCCGGACTGCACGACAGTCGAGTCGATGGAGTCGCGAGCGCTATGAGCGACGATGGAGTCGAGTCGAGACGCTGACGAAATATGGAGTCGATACGTAGCACGCAGAACGATGGGAGCT"
        
        results = analyze_sequences(
            [query],
            config='Col1a1',
            method='exact',
            min_sequence_length=20,
            verbose=False,
            merge_adjacent_mutations=True
        )
        
        # Extract mutation HGVS strings
        mutations = [m.to_hgvs() for m in results.mutations[0]]
        expected = ["23_76delinsGT", "104_211del", "238_239insAA", "265_266insA"]
        
        assert mutations == expected, f"Expected {expected}, got {mutations}"
    
    def test_case_2_mixed_delins_identification(self):
        """测试案例2：51_211delinsA的正确识别（混合插入和删除）"""
        # Query sequence from issues.tmp line 110 (removing dashes)
        query = "CGCCGGACTGCACGACAGTCGAGATGGAGTCGACACGACTCGCGCATACACGATGGAGTCGAGTCGAGACGCTGACGCATATGGAGTCGATACGTAGCACGCAGAGGCGATGGGAGCT"
        
        results = analyze_sequences(
            [query],
            config='Col1a1',
            method='exact',
            min_sequence_length=20,
            verbose=False,
            merge_adjacent_mutations=True
        )
        
        # Extract mutation HGVS strings
        mutations = [m.to_hgvs() for m in results.mutations[0]]
        expected = ["23_23del", "51_211delinsA", "238_239insC", "265_266insGG"]
        
        assert mutations == expected, f"Expected {expected}, got {mutations}"
    
    def test_case_3_cross_motif_delins_identification(self):
        """测试案例3：3_97delinsA的正确识别（跨motif的delins事件）
        
        这个测试验证了修复后的代码能够正确识别跨motif边界的delins事件。
        之前的问题是INDEL类型的突变在identify_cas9_events中被过滤掉了。
        修复后，INDEL和COMPLEX类型的突变都会被正确保留并合并。
        """
        # Query sequence from issues.tmp line 119 (removing dashes)
        # Original: CG-----------------------------------------------------------------------------------------------AGAGCG--------------------------CGCTATGGAGTCGAGAGCGCGCTCGTCGA-----------------------------------------------------------------------------------------------------------ACGATGGGAGCT
        query = "CGAGAGCGCGCTATGGAGTCGAGAGCGCGCTCGTCGAACGATGGGAGCT"
        
        results = analyze_sequences(
            [query],
            config='Col1a1',
            method='exact',
            min_sequence_length=20,
            verbose=False,
            merge_adjacent_mutations=True
        )
        
        # Extract mutation HGVS strings
        mutations = [m.to_hgvs() for m in results.mutations[0]]
        expected = ["3_97delinsA", "103_128del", "158_264del"]
        
        assert mutations == expected, f"Expected {expected}, got {mutations}"
        
        # 验证第一个突变是delins类型（合并后变成COMPLEX类型）
        first_mutation = results.mutations[0][0]
        assert first_mutation.loc_start == 3
        assert first_mutation.loc_end == 97
        assert first_mutation.seq_new == "A", f"Expected inserted sequence 'A', got '{first_mutation.seq_new}'"
        assert first_mutation.type.value in ["DI", "C"], f"Expected INDEL or COMPLEX type, got {first_mutation.type.value}"
    
    def test_case_4_single_base_deletion_format(self):
        """测试案例4：单碱基缺失的HGVS格式（2_2del而不是2del）
        
        这个测试验证了修复后的代码能够正确格式化单碱基缺失。
        根据HGVS标准，即使是单碱基缺失也应该使用范围格式（如2_2del）。
        """
        # 直接测试Mutation对象的to_hgvs方法，验证单碱基缺失的格式
        from darlin.mutations.mutation import Mutation, MutationType
        
        # 创建单碱基缺失突变（位置2）
        single_base_del = Mutation(
            type=MutationType.DELETION,
            loc_start=2,
            loc_end=2,
            seq_old="G",
            seq_new=""
        )
        
        # 验证格式为2_2del而不是2del
        hgvs_str = single_base_del.to_hgvs()
        assert hgvs_str == "2_2del", f"Expected '2_2del', got '{hgvs_str}'"
        
        # 测试另一个单碱基缺失
        another_del = Mutation(
            type=MutationType.DELETION,
            loc_start=23,
            loc_end=23,
            seq_old="A",
            seq_new=""
        )
        assert another_del.to_hgvs() == "23_23del", f"Expected '23_23del', got '{another_del.to_hgvs()}'"
        
        # 验证多碱基缺失仍然使用范围格式
        multi_base_del = Mutation(
            type=MutationType.DELETION,
            loc_start=24,
            loc_end=35,
            seq_old="CGATGGAGTCG",
            seq_new=""
        )
        assert multi_base_del.to_hgvs() == "24_35del", f"Expected '24_35del', got '{multi_base_del.to_hgvs()}'"
    
    def test_case_5_large_delins_with_insertion(self):
        """测试案例5：大片段删除插入的正确识别（14_265delinsAGT）
        
        这个测试验证了修复后的代码能够正确识别大片段删除后跟小片段插入的delins事件。
        之前的问题是代码在找到2个连续匹配后停止检查，错过了后续的替换/插入。
        修复后，代码要求至少3个连续匹配才停止，并且会检查后续位置是否有替换。
        """
        # Query sequence from issues.tmp line 136
        query = "CGCCGGACTGCACAGTCGATGGGAGCT"
        
        results = analyze_sequences(
            [query],
            config='Col1a1',
            method='exact',
            min_sequence_length=20,
            verbose=False,
            merge_adjacent_mutations=True
        )
        
        # Extract mutation HGVS strings
        mutations = [m.to_hgvs() for m in results.mutations[0]]
        expected = ["14_265delinsAGT"]
        
        assert mutations == expected, f"Expected {expected}, got {mutations}"
        
        # Verify the mutation details
        first_mutation = results.mutations[0][0]
        assert first_mutation.loc_start == 14
        assert first_mutation.loc_end == 265
        assert first_mutation.seq_new == "AGT", f"Expected inserted sequence 'AGT', got '{first_mutation.seq_new}'"
        assert first_mutation.type.value in ["DI", "C"], f"Expected INDEL or COMPLEX type, got {first_mutation.type.value}"
    
    def test_case_6_delins_with_insertion_before_gap(self):
        """测试案例6：gap之前有插入的delins事件正确识别（19_270delinsGGGA）
        
        这个测试验证了修复后的代码能够正确识别在gap之前有替换/插入的delins事件。
        之前的问题是：
        1. 代码在找到gap之前没有检查替换，导致gap之前的插入被遗漏
        2. HGVS规范化错误地trim了不对齐的suffix，导致插入序列被截断
        修复后，代码会：
        1. 在找到gap之前检查替换，将突变区域扩展到替换的开始位置
        2. 改进HGVS规范化逻辑，只有当序列长度差<=1时才trim suffix（除非suffix只有1个字符）
        """
        # Query sequence from issues.tmp line 144
        query = "CGCCGGACTGCACGACAGGGGAGGAGCT"
        
        results = analyze_sequences(
            [query],
            config='Col1a1',
            method='exact',
            min_sequence_length=20,
            verbose=False,
            merge_adjacent_mutations=True
        )
        
        # Extract mutation HGVS strings
        mutations = [m.to_hgvs() for m in results.mutations[0]]
        expected = ["19_270delinsGGGA"]
        
        assert mutations == expected, f"Expected {expected}, got {mutations}"
        
        # Verify the mutation details
        first_mutation = results.mutations[0][0]
        assert first_mutation.loc_start == 19
        assert first_mutation.loc_end == 270
        assert first_mutation.seq_new == "GGGA", f"Expected inserted sequence 'GGGA', got '{first_mutation.seq_new}'"
        assert first_mutation.type.value in ["DI", "C"], f"Expected INDEL or COMPLEX type, got {first_mutation.type.value}"
    
    def test_case_7_complex_delins_with_multiple_insertions(self):
        """测试案例7：复杂delins事件中多个插入的正确识别（97_244delinsGCGCGTATCGTATCGACTCCAT）
        
        这个测试验证了修复后的代码能够正确识别包含多个插入片段的复杂delins事件。
        之前的问题是HGVS规范化错误地trim了不对齐的prefix，导致插入序列被截断。
        修复后，代码会：
        1. 检查序列长度差异，只有当长度差<=1时才trim prefix（除非prefix只有1个字符且长度差<=3）
        2. 避免trim不对齐的prefix，确保插入序列完整
        """
        # Query sequence from issues.tmp line 152
        query = "CGCCGGACTGCACGACAG----------------------------------------------------------CGATGGAGTCGCGAGCGCTA------------------------------------------------GCGCG----TATCGTATCGACTCCAT-------------------------------------------------------------------------------------AGTCGATACGTAGCACGCAGAGGCGATGGGAGCT"
        query_clean = query.replace('-', '')
        
        results = analyze_sequences(
            [query_clean],
            config='Col1a1',
            method='exact',
            min_sequence_length=20,
            verbose=False,
            merge_adjacent_mutations=True
        )
        
        # Extract mutation HGVS strings
        mutations = [m.to_hgvs() for m in results.mutations[0]]
        expected = ["19_76del", "97_244delinsGCGCGTATCGTATCGACTCCAT", "265_266insGG"]
        
        assert mutations == expected, f"Expected {expected}, got {mutations}"
        
        # Verify the second mutation (delins) details
        if len(results.mutations[0]) >= 2:
            second_mutation = results.mutations[0][1]
            assert second_mutation.loc_start == 97
            assert second_mutation.loc_end == 244
            assert second_mutation.seq_new == "GCGCGTATCGTATCGACTCCAT", f"Expected inserted sequence 'GCGCGTATCGTATCGACTCCAT', got '{second_mutation.seq_new}'"
            assert second_mutation.type.value in ["DI", "C"], f"Expected INDEL or COMPLEX type, got {second_mutation.type.value}"


if __name__ == "__main__":
    pytest.main([__file__]) 