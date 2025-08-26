#!/usr/bin/env python3
"""
DARLIN主API功能演示

展示如何使用analyze_sequences函数进行完整的CARLIN序列分析
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from darlin import analyze_sequences, quick_analyze, ORIGINAL_CARLIN
from darlin.utils import reverse_complement, calculate_gc_content


def demo_basic_analysis():
    """演示基本的序列分析"""
    print("=== 基本序列分析演示 ===\n")
    
    # 获取CARLIN参考序列并创建测试数据
    carlin_ref = ORIGINAL_CARLIN.get_full_reference_sequence()
    print(f"CARLIN参考序列长度: {len(carlin_ref)} bp")
    print(f"GC含量: {calculate_gc_content(carlin_ref):.1f}%")
    print()
    
    # 创建测试序列（包含一些变异）
    sequences = [
        carlin_ref[:250],                           # 前250bp
        carlin_ref[50:300],                         # 中间250bp
        carlin_ref[-250:],                          # 后250bp
        reverse_complement(carlin_ref[:250]),       # 反向互补
        carlin_ref[:200] + "N" * 50,               # 含有N的序列
    ]
    
    print(f"准备分析 {len(sequences)} 条测试序列...")
    
    # 进行分析
    result = analyze_sequences(
        sequences,
        config='OriginalCARLIN',
        method='coarse_grain',
        verbose=True
    )
    
    print("\n" + "="*50)
    print("分析结果详情:")
    print("="*50)
    
    # 打印详细结果
    for i, (allele, mutations) in enumerate(zip(result.called_alleles, result.mutations)):
        print(f"\n序列 {i+1}:")
        print(f"  比对得分: {result.alignment_scores[i]:.2f}")
        
        if allele.is_callable():
            print(f"  等位基因调用: 成功")
            print(f"  置信度: {allele.confidence:.2f}")
            print(f"  事件结构: {allele.event_structure}")
            print(f"  检测到突变: {len(mutations)} 个")
            
            for j, mut in enumerate(mutations):
                print(f"    突变 {j+1}: {mut.to_hgvs()}")
        else:
            print(f"  等位基因调用: 失败")


def demo_different_methods():
    """演示不同的分析方法"""
    print("\n=== 不同分析方法对比 ===\n")
    
    carlin_ref = ORIGINAL_CARLIN.get_full_reference_sequence()
    test_sequences = [carlin_ref[:250], carlin_ref[100:350]]
    
    methods = ['coarse_grain', 'exact']
    
    for method in methods:
        print(f"使用 {method} 方法:")
        result = analyze_sequences(
            test_sequences,
            method=method,
            verbose=False
        )
        
        print(f"  处理时间: {result.processing_time:.3f}秒")
        print(f"  调用成功率: {result.calling_success_rate:.1%}")
        print(f"  平均比对得分: {result.average_alignment_score:.2f}")
        print(f"  检测突变总数: {result.total_mutations}")
        print()


def demo_custom_parameters():
    """演示自定义参数的使用"""
    print("=== 自定义参数演示 ===\n")
    
    carlin_ref = ORIGINAL_CARLIN.get_full_reference_sequence()
    sequences = [carlin_ref[:250]]
    
    # 测试不同的参数组合
    param_sets = [
        {
            'dominant_threshold': 0.3,
            'annotate_mutations_flag': True,
            'merge_adjacent_mutations': True,
            'description': '低阈值，完整注释'
        },
        {
            'dominant_threshold': 0.8,
            'annotate_mutations_flag': False,
            'merge_adjacent_mutations': False,
            'description': '高阈值，无突变注释'
        }
    ]
    
    for params in param_sets:
        desc = params.pop('description')
        print(f"{desc}:")
        
        result = analyze_sequences(sequences, verbose=False, **params)
        
        print(f"  调用成功率: {result.calling_success_rate:.1%}")
        print(f"  检测突变数: {result.total_mutations}")
        print(f"  处理时间: {result.processing_time:.3f}秒")
        print()


def demo_quick_analyze():
    """演示quick_analyze便利函数"""
    print("=== 快速分析演示 ===\n")
    
    # 使用CARLIN参考序列的片段作为测试序列
    carlin_ref = ORIGINAL_CARLIN.get_full_reference_sequence()
    test_sequences = [
        carlin_ref[:180],
        carlin_ref[100:280], 
        carlin_ref[200:380],
    ]
    
    print("使用quick_analyze进行快速分析...")
    result = quick_analyze(test_sequences)
    
    print(f"\n快速分析摘要:")
    print(f"配置: {result.config_used}")
    print(f"方法: {result.method_used}")
    print(f"序列数: {result.num_sequences}")
    print(f"调用成功: {result.num_called_alleles}")


def demo_utils_integration():
    """演示工具函数与主API的集成使用"""
    print("\n=== 工具函数集成演示 ===\n")
    
    from darlin.utils import (
        reverse_complement, 
        format_sequence,
        to_fasta_string
    )
    
    # 创建一个示例序列
    carlin_ref = ORIGINAL_CARLIN.get_full_reference_sequence()
    test_seq = carlin_ref[:300]
    
    print("原始序列分析:")
    print(f"长度: {len(test_seq)} bp")
    print(f"GC含量: {calculate_gc_content(test_seq):.1f}%")
    
    # 分析原始序列
    result_forward = analyze_sequences([test_seq], verbose=False)
    
    # 分析反向互补序列
    rc_seq = reverse_complement(test_seq)
    result_rc = analyze_sequences([rc_seq], verbose=False)
    
    print(f"\n正向序列比对得分: {result_forward.alignment_scores[0]:.2f}")
    print(f"反向互补序列比对得分: {result_rc.alignment_scores[0]:.2f}")
    
    # 格式化输出
    print("\n格式化序列输出:")
    formatted = format_sequence(test_seq[:80], line_length=20, add_numbers=True)
    print(formatted)


def demo_error_handling():
    """演示错误处理"""
    print("\n=== 错误处理演示 ===\n")
    
    # 演示各种错误情况
    error_cases = [
        {
            'name': '空序列列表',
            'sequences': [],
            'params': {}
        },
        {
            'name': '过短序列',
            'sequences': ['ATCG', 'TTGG'],
            'params': {}
        },
        {
            'name': '无效配置',
            'sequences': ['A' * 100],
            'params': {'config': 'InvalidConfig'}
        },
        {
            'name': '无效方法',
            'sequences': ['A' * 100],
            'params': {'method': 'invalid_method'}
        }
    ]
    
    for case in error_cases:
        print(f"测试错误情况: {case['name']}")
        try:
            analyze_sequences(case['sequences'], **case['params'])
            print("  ❌ 未捕获到预期错误")
        except Exception as e:
            print(f"  ✅ 正确捕获错误: {type(e).__name__}: {str(e)[:50]}...")
        print()


def main():
    """主演示函数"""
    print("DARLIN主API功能演示")
    print("=" * 50)
    print()
    
    try:
        demo_basic_analysis()
        demo_different_methods()
        demo_custom_parameters()
        demo_quick_analyze()
        demo_utils_integration()
        demo_error_handling()
        
        print("\n" + "="*50)
        print("演示完成！")
        print("\n主API功能特性：")
        print("✓ 一站式序列分析接口 (analyze_sequences)")
        print("✓ 灵活的配置和参数选项")
        print("✓ 完整的结果数据结构 (AnalysisResult)")
        print("✓ 详细的统计信息和摘要")
        print("✓ 快速分析便利函数 (quick_analyze)")
        print("✓ 丰富的序列处理工具 (utils)")
        print("✓ 健壮的错误处理机制")
        print("✓ 与CARLIN生物学特性的深度集成")
        
    except Exception as e:
        print(f"演示过程中出现错误: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main() 