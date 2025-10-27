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


if __name__ == "__main__":
    pytest.main([__file__]) 