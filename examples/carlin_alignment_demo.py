#!/usr/bin/env python3
"""
CARLIN比对系统完整演示

展示如何使用集成的CARLIN比对器进行真实的序列分析
"""

from darlinpy.alignment import CARLINAligner, align_to_carlin, create_default_aligner
from darlinpy.config import get_original_carlin_config


def demo_perfect_carlin_sequence():
    """演示完美的CARLIN序列比对"""
    print("=== 演示1: 完美CARLIN序列 ===")
    
    # 获取原始CARLIN参考序列
    config = get_original_carlin_config()
    perfect_sequence = config.get_reference_sequence()
    
    # 进行比对
    result = align_to_carlin(perfect_sequence, verbose=True)
    
    print(f"✅ 完美序列比对完成")
    print(f"   - 比对得分: {result['alignment_score']:.2f}")
    print(f"   - 序列一致性: {result['statistics']['identity']*100:.1f}%")
    print()


def demo_carlin_with_mutations():
    """演示带突变的CARLIN序列"""
    print("=== 演示2: 带突变的CARLIN序列 ===")
    
    # 获取原始序列并引入一些突变
    config = get_original_carlin_config()
    original_seq = config.get_reference_sequence()
    
    # 模拟一些真实的CRISPR编辑事件
    mutations = [
        (50, "TTTT"),    # 在位置50插入TTTT
        (120, ""),       # 在位置120删除3个碱基
        (123, ""),
        (124, ""),
        (200, "A")       # 位置200的C->A替换
    ]
    
    # 应用突变
    mutated_seq = list(original_seq)
    for pos, insert_seq in sorted(mutations, reverse=True):  # 从后往前应用
        if insert_seq:  # 插入
            mutated_seq.insert(pos, insert_seq)
        else:  # 删除
            if pos < len(mutated_seq):
                del mutated_seq[pos]
    
    # 单独处理替换
    mutated_seq = ''.join(mutated_seq)
    if len(mutated_seq) > 200:
        mutated_seq = mutated_seq[:200] + 'A' + mutated_seq[201:]
    
    print(f"原始序列长度: {len(original_seq)} bp")
    print(f"突变序列长度: {len(mutated_seq)} bp")
    print("引入的突变: 插入(pos 50), 删除(pos 120-122), 替换(pos 200)")
    
    # 进行比对
    result = align_to_carlin(mutated_seq, verbose=True)
    
    print(f"✅ 突变序列比对完成")
    print()


def demo_batch_alignment():
    """演示批量序列比对"""
    print("=== 演示3: 批量序列比对 ===")
    
    config = get_original_carlin_config()
    ref_seq = config.get_reference_sequence()
    
    # 创建一批测试序列
    test_sequences = [
        ref_seq,  # 完美匹配
        ref_seq[:200],  # 截断序列
        ref_seq[:100] + "AAAAA" + ref_seq[100:200],  # 中间插入
        ref_seq.replace(ref_seq[50:55], "NNNNN"),  # 含N的序列
        ref_seq[:150] + ref_seq[160:],  # 删除片段
    ]
    
    sequence_names = [
        "完美序列",
        "截断序列", 
        "插入序列",
        "含N序列",
        "删除序列"
    ]
    
    # 创建比对器
    aligner = create_default_aligner()
    
    print(f"批量比对 {len(test_sequences)} 个序列...")
    
    # 批量比对
    results = aligner.align_sequences(test_sequences, verbose=False)
    
    # 汇总结果
    print("\n=== 批量比对结果汇总 ===")
    print(f"{'序列名称':<12} {'长度':<8} {'得分':<10} {'一致性':<10} {'状态'}")
    print("-" * 55)
    
    for i, (name, result) in enumerate(zip(sequence_names, results)):
        if 'error' in result:
            print(f"{name:<12} {len(test_sequences[i]):<8} {'ERROR':<10} {'N/A':<10} {'失败'}")
        else:
            seq_len = len(result['query_sequence'])
            score = result['alignment_score']
            identity = result['statistics']['identity'] * 100
            status = "成功"
            
            print(f"{name:<12} {seq_len:<8} {score:<10.1f} {identity:<10.1f}% {status}")
    
    print()


def demo_cutsite_analysis():
    """演示CRISPR切割位点分析"""
    print("=== 演示4: CRISPR切割位点分析 ===")
    
    config = get_original_carlin_config()
    ref_seq = config.get_reference_sequence()
    
    # 模拟在特定切割位点的编辑
    # 获取第一个segment的cutsite位置
    cutsite_start, cutsite_end = config.positions['cutsites'][0]
    
    print(f"目标切割位点: Segment 1, 位置 {cutsite_start}-{cutsite_end}")
    print(f"切割位点序列: {ref_seq[cutsite_start:cutsite_end]}")
    
    # 模拟典型的CRISPR编辑：在切割位点插入/删除
    edited_sequences = [
        # 在切割位点插入
        ref_seq[:cutsite_start+3] + "AGCT" + ref_seq[cutsite_start+3:],
        # 在切割位点删除
        ref_seq[:cutsite_start+2] + ref_seq[cutsite_start+5:],
        # 切割位点替换
        ref_seq[:cutsite_start] + "AAAAAAA" + ref_seq[cutsite_end:],
    ]
    
    edit_types = ["插入AGCT", "删除3bp", "替换cutsite"]
    
    aligner = create_default_aligner()
    
    for i, (edit_seq, edit_type) in enumerate(zip(edited_sequences, edit_types)):
        print(f"\n--- 编辑类型: {edit_type} ---")
        result = aligner.align_sequence(edit_seq, verbose=False)
        
        stats = result['statistics']
        print(f"比对得分: {result['alignment_score']:.2f}")
        print(f"序列一致性: {stats['identity']*100:.1f}%")
        print(f"Gap数量: 查询={stats['query_gaps']}, 参考={stats['reference_gaps']}")
        
        # 显示切割位点附近的比对情况
        aligned_query = result['aligned_query']
        aligned_ref = result['aligned_reference']
        
        # 找到切割位点在比对中的大致位置
        region_start = max(0, cutsite_start - 10)
        region_end = min(len(aligned_ref), cutsite_end + 10)
        
        print(f"切割位点附近比对:")
        if region_end <= len(aligned_query):
            print(f"Query: {aligned_query[region_start:region_end]}")
            print(f"Ref  : {aligned_ref[region_start:region_end]}")
    
    print()


def demo_position_specific_penalties():
    """演示位置特异性gap惩罚的效果"""
    print("=== 演示5: 位置特异性Gap惩罚效果 ===")
    
    config = get_original_carlin_config()
    ref_seq = config.get_reference_sequence()
    
    # 展示不同区域的gap惩罚
    print("Gap惩罚配置:")
    print(f"Prefix区域: {config.open_penalty.prefix}")
    print(f"Cutsite区域: {config.open_penalty.cutsites}")
    print(f"Consite区域: {config.open_penalty.consites[:5]}... (前5个值)")
    print(f"PAM区域: {config.open_penalty.pam}")
    print(f"Postfix区域: {config.open_penalty.postfix}")
    
    # 在不同区域插入相同的序列，观察gap惩罚的影响
    insertion = "AAAA"
    test_positions = [
        (2, "Prefix区域"),      # prefix内
        (50, "Cutsite区域"),    # 第一个cutsite内
        (30, "Consite区域"),    # 第一个consite内
        (270, "Postfix区域"),   # postfix内
    ]
    
    aligner = create_default_aligner()
    
    print(f"\n在不同区域插入相同序列'{insertion}'的比对结果:")
    print(f"{'区域':<12} {'插入位置':<8} {'比对得分':<10} {'Gap处理'}")
    print("-" * 45)
    
    for pos, region_name in test_positions:
        if pos < len(ref_seq):
            test_seq = ref_seq[:pos] + insertion + ref_seq[pos:]
            result = aligner.align_sequence(test_seq, verbose=False)
            
            score = result['alignment_score']
            stats = result['statistics']
            gap_info = f"Q:{stats['query_gaps']},R:{stats['reference_gaps']}"
            
            print(f"{region_name:<12} {pos:<8} {score:<10.1f} {gap_info}")
    
    print()


def main():
    """主演示函数"""
    print("🧬 CARLIN比对系统完整演示")
    print("=" * 60)
    print()
    
    # 显示配置信息
    aligner = create_default_aligner()
    print(aligner.get_config_summary())
    print()
    
    # 运行所有演示
    demo_perfect_carlin_sequence()
    demo_carlin_with_mutations() 
    demo_batch_alignment()
    demo_cutsite_analysis()
    demo_position_specific_penalties()
    
    print("🎉 CARLIN比对演示完成！")
    print()
    print("💡 关键特性总结:")
    print("   ✅ 完整的CARLIN配置管理")
    print("   ✅ 位置特异性gap惩罚")
    print("   ✅ 详细的比对统计分析")
    print("   ✅ CRISPR编辑事件检测")
    print("   ✅ 批量序列处理")
    print("   ✅ Motif级别的分析")


if __name__ == "__main__":
    main() 