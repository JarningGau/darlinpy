#!/usr/bin/env python3
"""
基准测试：比较cas9_align_py和_cas9_align_module.cas9_align的性能
"""

import numpy as np
import time
import sys
import os
from typing import List, Tuple
import random

# 添加路径以便导入模块
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from darlinpy.alignment.cas9_align import (
    cas9_align_py,
    nt2int,
    HAS_CPP_IMPL,
    _cas9_align_module,
)


class CAS9AlignBenchmark:
    """CAS9对齐算法基准测试类"""
    
    def __init__(self):
        """初始化测试参数"""
        # 设置随机种子以确保可重复性
        random.seed(42)
        np.random.seed(42)
        
        # 基础的评分矩阵 (简化的NUC44)
        self.sub_score = np.zeros(25, dtype=np.float64)
        for i in range(1, 5):  # A=1, C=2, G=3, T=4
            for j in range(1, 5):
                if i == j:
                    self.sub_score[i * 5 + j] = 5.0  # 匹配
                else:
                    self.sub_score[i * 5 + j] = -4.0  # 不匹配
        
        # 检查C++实现是否可用
        self.has_cpp = HAS_CPP_IMPL
        if self.has_cpp:
            print("✓ C++扩展模块可用")
        else:
            print("⚠ C++扩展模块不可用，只测试Python实现")
    
    def generate_test_sequences(self, length: int, num_pairs: int = 10) -> List[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
        """生成测试序列对"""
        test_cases = []
        
        for _ in range(num_pairs):
            # 生成随机DNA序列
            seq_str = ''.join(random.choices('ACGT', k=length))
            ref_str = ''.join(random.choices('ACGT', k=length))
            
            # 转换为整数编码
            seq_int = nt2int(seq_str)
            ref_int = nt2int(ref_str)
            
            # 生成位置特定的gap penalty
            open_penalty = np.random.uniform(0.5, 2.0, len(ref_int) + 1)
            close_penalty = np.random.uniform(0.1, 0.5, len(ref_int) + 1)
            
            test_cases.append((seq_int, ref_int, open_penalty, close_penalty))
        
        return test_cases
    
    def benchmark_function(self, func, test_cases: List[Tuple], name: str, num_runs: int = 3):
        """基准测试单个函数"""
        print(f"\n=== 基准测试: {name} ===")
        
        times = []
        for run in range(num_runs):
            start_time = time.perf_counter()
            
            for seq, ref, open_pen, close_pen in test_cases:
                try:
                    score, al_seq, al_ref = func(seq, ref, open_pen, close_pen, self.sub_score)
                except Exception as e:
                    print(f"错误: {e}")
                    return None
            
            end_time = time.perf_counter()
            elapsed = end_time - start_time
            times.append(elapsed)
            print(f"运行 {run + 1}: {elapsed:.4f} 秒")
        
        avg_time = np.mean(times)
        std_time = np.std(times)
        min_time = np.min(times)
        
        print(f"平均时间: {avg_time:.4f} ± {std_time:.4f} 秒")
        print(f"最快时间: {min_time:.4f} 秒")
        print(f"总测试用例: {len(test_cases)}")
        print(f"每用例平均时间: {avg_time/len(test_cases)*1000:.2f} 毫秒")
        
        return {
            'name': name,
            'avg_time': avg_time,
            'std_time': std_time,
            'min_time': min_time,
            'num_cases': len(test_cases),
            'time_per_case': avg_time / len(test_cases)
        }
    
    def run_benchmark(self, sequence_lengths: List[int] = [50, 100, 200, 500], 
                     num_pairs: int = 10, num_runs: int = 3):
        """运行完整的基准测试"""
        print("=" * 60)
        print("CAS9对齐算法性能基准测试")
        print("=" * 60)
        
        results = []
        
        for length in sequence_lengths:
            print(f"\n{'='*20} 序列长度: {length} {'='*20}")
            
            # 生成测试数据
            test_cases = self.generate_test_sequences(length, num_pairs)
            
            # 测试Python实现
            py_result = self.benchmark_function(
                cas9_align_py, test_cases, 
                f"Python实现 (长度={length})", num_runs
            )
            if py_result:
                results.append(py_result)
            
            # 测试C++实现（如果可用）
            if self.has_cpp:
                cpp_result = self.benchmark_function(
                    _cas9_align_module.cas9_align, test_cases,
                    f"C++实现 (长度={length})", num_runs
                )
                if cpp_result:
                    results.append(cpp_result)
                    
                    # 计算加速比
                    if py_result:
                        speedup = py_result['time_per_case'] / cpp_result['time_per_case']
                        print(f"🚀 C++加速比: {speedup:.2f}x")
        
        # 生成总结报告
        self.print_summary(results)
        
        return results
    
    def print_summary(self, results: List[dict]):
        """打印基准测试总结"""
        print("\n" + "=" * 60)
        print("基准测试总结")
        print("=" * 60)
        
        # 按实现类型分组
        py_results = [r for r in results if 'Python' in r['name']]
        cpp_results = [r for r in results if 'C++' in r['name']]
        
        if py_results:
            print("\nPython实现性能:")
            for result in py_results:
                length = result['name'].split('长度=')[1].split(')')[0]
                print(f"  长度 {length}: {result['time_per_case']*1000:.2f} ms/用例")
        
        if cpp_results:
            print("\nC++实现性能:")
            for result in cpp_results:
                length = result['name'].split('长度=')[1].split(')')[0]
                print(f"  长度 {length}: {result['time_per_case']*1000:.2f} ms/用例")
        
        # 计算平均加速比
        if py_results and cpp_results:
            print("\n平均加速比:")
            for py_result in py_results:
                py_length = py_result['name'].split('长度=')[1].split(')')[0]
                for cpp_result in cpp_results:
                    cpp_length = cpp_result['name'].split('长度=')[1].split(')')[0]
                    if py_length == cpp_length:
                        speedup = py_result['time_per_case'] / cpp_result['time_per_case']
                        print(f"  长度 {py_length}: {speedup:.2f}x")


def main():
    """主函数"""
    benchmark = CAS9AlignBenchmark()
    
    # 运行基准测试
    # 可以调整参数：序列长度、测试用例数量、运行次数
    results = benchmark.run_benchmark(
        sequence_lengths=[50, 100, 200],  # 测试不同长度的序列
        num_pairs=5,                      # 每个长度的测试用例数
        num_runs=3                        # 每个测试运行次数
    )
    
    print(f"\n基准测试完成！共测试了 {len(results)} 个配置。")


if __name__ == "__main__":
    main()
