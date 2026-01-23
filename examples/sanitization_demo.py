#!/usr/bin/env python3
"""
CARLIN序列标准化演示

展示prefix/postfix和保守区域标准化的实际效果
"""

from darlinpy.alignment import CARLINAligner, create_default_aligner, AlignedSEQ, SequenceSanitizer
from darlinpy.config import get_original_carlin_config


def demo_prefix_postfix_sanitization():
    """演示prefix/postfix标准化"""
    print("=== 演示1: Prefix/Postfix标准化 ===")
    
    config = get_original_carlin_config()
    ref_seq = config.get_reference_sequence()
    
    # 在prefix和postfix区域添加人工制品
    # 这些通常是PCR扩增或测序过程中的技术错误
    contaminated_seq = (
        "NNNNN" +           # 前缀污染
        ref_seq +           # 核心CARLIN序列
        "TTTTT"             # 后缀污染
    )
    
    print(f"原始CARLIN长度: {len(ref_seq)} bp")
    print(f"污染序列长度: {len(contaminated_seq)} bp")
    print(f"前缀污染: {contaminated_seq[:5]}")
    print(f"后缀污染: {contaminated_seq[-5:]}")
    
    aligner = create_default_aligner()
    
    # 不使用标准化
    result_raw = aligner.align_sequence(contaminated_seq, verbose=False, sanitize=False)
    
    # 使用标准化
    result_clean = aligner.align_sequence(contaminated_seq, verbose=False, sanitize=True)
    
    print(f"\n比对结果比较:")
    print(f"未标准化:")
    print(f"  - 比对得分: {result_raw['alignment_score']:.1f}")
    print(f"  - 序列一致性: {result_raw['statistics']['identity']*100:.1f}%")
    print(f"  - 比对长度: {result_raw['statistics']['aligned_length']}")
    
    print(f"标准化后:")
    print(f"  - 比对得分: {result_clean['alignment_score']:.1f}")
    print(f"  - 序列一致性: {result_clean['statistics']['identity']*100:.1f}%")
    print(f"  - 比对长度: {result_clean['statistics']['aligned_length']}")
    print(f"  - 标准化成功: {result_clean['sanitized']}")
    
    print("✅ Prefix/Postfix标准化演示完成\n")


def demo_conserved_region_sanitization():
    """演示保守区域标准化"""
    print("=== 演示2: 保守区域标准化 ===")
    
    config = get_original_carlin_config()
    ref_seq = config.get_reference_sequence()
    
    # 在保守区域(PAM、consite)引入"测序错误"
    # 这些错误在生物学上不太可能，更可能是技术问题
    
    # 在第一个PAM序列中引入错误
    pam_start, pam_end = config.positions['pams'][0]
    pam_sequence = ref_seq[pam_start:pam_end]
    
    # 在第一个consite中引入错误
    consite_start, consite_end = config.positions['consites'][0]
    consite_sequence = ref_seq[consite_start:consite_end]
    
    print(f"原始PAM序列: {pam_sequence}")
    print(f"原始Consite序列: {consite_sequence}")
    
    # 创建带错误的序列
    error_seq = list(ref_seq)
    
    # PAM错误: TGGAGTC -> TGTAGTC (测序错误)
    if pam_start + 2 < len(error_seq):
        error_seq[pam_start + 2] = 'T'  # G->T
    
    # Consite错误 (在保守区域引入单个错配)
    if consite_start + 5 < len(error_seq):
        error_seq[consite_start + 5] = 'T' if error_seq[consite_start + 5] != 'T' else 'A'
    
    error_seq_str = ''.join(error_seq)
    
    print(f"错误PAM序列: {error_seq_str[pam_start:pam_end]}")
    print(f"错误Consite序列: {error_seq_str[consite_start:consite_end]}")
    
    aligner = create_default_aligner()
    
    # 不使用标准化
    result_with_errors = aligner.align_sequence(error_seq_str, verbose=False, sanitize=False)
    
    # 使用标准化
    result_corrected = aligner.align_sequence(error_seq_str, verbose=False, sanitize=True)
    
    print(f"\n保守区域错误修正结果:")
    print(f"有错误的序列:")
    print(f"  - 比对得分: {result_with_errors['alignment_score']:.1f}")
    print(f"  - 序列一致性: {result_with_errors['statistics']['identity']*100:.1f}%")
    print(f"  - 不匹配数: {result_with_errors['statistics']['mismatches']}")
    
    print(f"标准化修正后:")
    print(f"  - 比对得分: {result_corrected['alignment_score']:.1f}")
    print(f"  - 序列一致性: {result_corrected['statistics']['identity']*100:.1f}%")
    print(f"  - 不匹配数: {result_corrected['statistics']['mismatches']}")
    print(f"  - 标准化成功: {result_corrected['sanitized']}")
    
    if result_corrected['sanitized']:
        sanitized_obj = result_corrected['aligned_seq_obj']
        events = sanitized_obj.get_event_structure()
        print(f"  - 事件结构: {''.join(events)}")
        print(f"  - 错配事件数: {events.count('M')}")
    
    print("✅ 保守区域标准化演示完成\n")


def demo_cutsite_preservation():
    """演示cutsite区域保留真实编辑"""
    print("=== 演示3: Cutsite编辑保留 ===")
    
    config = get_original_carlin_config()
    ref_seq = config.get_reference_sequence()
    
    # 在cutsite区域引入真实的CRISPR编辑
    cutsite_start, cutsite_end = config.positions['cutsites'][0]
    cutsite_sequence = ref_seq[cutsite_start:cutsite_end]
    
    print(f"目标Cutsite: {cutsite_sequence} (位置 {cutsite_start}-{cutsite_end})")
    
    # 模拟典型的CRISPR编辑：在cutsite中间删除2bp
    edited_seq = ref_seq[:cutsite_start+3] + ref_seq[cutsite_start+5:]
    
    print(f"编辑类型: 在cutsite位置{cutsite_start+3}-{cutsite_start+5}删除2bp")
    print(f"原始序列长度: {len(ref_seq)} bp")
    print(f"编辑序列长度: {len(edited_seq)} bp")
    
    aligner = create_default_aligner()
    result = aligner.align_sequence(edited_seq, verbose=False, sanitize=True)
    
    print(f"\nCutsite编辑分析:")
    print(f"  - 比对得分: {result['alignment_score']:.1f}")
    print(f"  - 序列一致性: {result['statistics']['identity']*100:.1f}%")
    print(f"  - 标准化成功: {result['sanitized']}")
    
    if result['sanitized']:
        sanitized_obj = result['aligned_seq_obj']
        events = sanitized_obj.get_event_structure()
        print(f"  - 事件结构: {''.join(events)}")
        print(f"  - 删除事件数: {events.count('D')}")
        print(f"  - 插入事件数: {events.count('I')}")
        print(f"  - 错配事件数: {events.count('M')}")
        
        # 分析具体哪些motif有编辑
        edited_motifs = []
        for i, event in enumerate(events):
            if event in ['D', 'I', 'M']:
                motif = sanitized_obj.motifs[i]
                edited_motifs.append((i, event, motif.seq, motif.ref))
        
        print(f"  - 检测到编辑的motifs:")
        for i, event, seq, ref in edited_motifs:
            print(f"    Motif {i}: {event} - '{seq}' vs '{ref}'")
    
    print("✅ Cutsite编辑保留演示完成\n")


def demo_batch_sanitization():
    """演示批量序列标准化"""
    print("=== 演示4: 批量序列标准化 ===")
    
    config = get_original_carlin_config()
    ref_seq = config.get_reference_sequence()
    
    # 创建一批带不同问题的序列
    test_sequences = [
        ref_seq,                                           # 完美序列
        "NNNNN" + ref_seq + "TTTTT",                      # 前后缀污染
        ref_seq[:50] + "A" + ref_seq[51:],                # 保守区域错误 (C->A替换)
        ref_seq[:100] + ref_seq[105:],                    # Cutsite删除
        "AAA" + ref_seq[:-3],                             # 前缀插入+后缀截断
    ]
    
    sequence_names = [
        "完美序列",
        "前后缀污染",
        "保守区域错误", 
        "Cutsite删除",
        "混合问题"
    ]
    
    aligner = create_default_aligner()
    
    print("批量标准化结果:")
    print(f"{'序列类型':<12} {'原长度':<8} {'标准化':<8} {'得分':<8} {'一致性':<8}")
    print("-" * 60)
    
    for i, (name, seq) in enumerate(zip(sequence_names, test_sequences)):
        result = aligner.align_sequence(seq, verbose=False, sanitize=True)
        
        sanitized = "✅" if result['sanitized'] else "❌"
        score = result['alignment_score']
        identity = result['statistics']['identity'] * 100
        
        print(f"{name:<12} {len(seq):<8} {sanitized:<8} {score:<8.1f} {identity:<8.1f}%")
    
    print("✅ 批量序列标准化演示完成\n")


def demo_sanitization_effects():
    """演示标准化对下游分析的影响"""
    print("=== 演示5: 标准化对分析的影响 ===")
    
    config = get_original_carlin_config()
    ref_seq = config.get_reference_sequence()
    
    # 创建一个复杂的问题序列
    complex_seq = (
        "NNNNN" +                    # 前缀污染
        ref_seq[:30] +               # 保留开头
        "A" +                        # 保守区域错误 (替换)
        ref_seq[31:100] +            # 中间部分
        ref_seq[105:200] +           # 删除5bp (模拟真实编辑)
        "T" +                        # 另一个保守区域错误 (替换)
        ref_seq[201:] +              # 结尾部分
        "TTTTT"                      # 后缀污染
    )
    
    print(f"复杂问题序列长度: {len(complex_seq)} bp")
    print("包含问题:")
    print("  - 前缀污染 (NNNNN)")
    print("  - 后缀污染 (TTTTT)")
    print("  - 保守区域错误 (替换)")
    print("  - 真实编辑 (5bp删除)")
    
    aligner = create_default_aligner()
    
    # 分别测试不使用和使用标准化
    result_raw = aligner.align_sequence(complex_seq, verbose=False, sanitize=False)
    result_clean = aligner.align_sequence(complex_seq, verbose=False, sanitize=True)
    
    print(f"\n标准化效果对比:")
    print(f"{'指标':<20} {'未标准化':<12} {'标准化后':<12} {'改善'}")
    print("-" * 60)
    
    # 比对得分
    score_diff = result_clean['alignment_score'] - result_raw['alignment_score']
    print(f"{'比对得分':<20} {result_raw['alignment_score']:<12.1f} {result_clean['alignment_score']:<12.1f} {score_diff:+.1f}")
    
    # 序列一致性
    identity_raw = result_raw['statistics']['identity'] * 100
    identity_clean = result_clean['statistics']['identity'] * 100
    identity_diff = identity_clean - identity_raw
    print(f"{'序列一致性(%)':<20} {identity_raw:<12.1f} {identity_clean:<12.1f} {identity_diff:+.1f}")
    
    # 不匹配数
    mismatch_diff = result_clean['statistics']['mismatches'] - result_raw['statistics']['mismatches']
    print(f"{'不匹配数':<20} {result_raw['statistics']['mismatches']:<12} {result_clean['statistics']['mismatches']:<12} {mismatch_diff:+d}")
    
    # Gap数量
    gaps_raw = result_raw['statistics']['query_gaps'] + result_raw['statistics']['reference_gaps']
    gaps_clean = result_clean['statistics']['query_gaps'] + result_clean['statistics']['reference_gaps']
    gaps_diff = gaps_clean - gaps_raw
    print(f"{'Gap总数':<20} {gaps_raw:<12} {gaps_clean:<12} {gaps_diff:+d}")
    
    if result_clean['sanitized'] and result_clean['aligned_seq_obj']:
        events = result_clean['aligned_seq_obj'].get_event_structure()
        print(f"\n标准化后的事件结构:")
        print(f"  - 事件字符串: {''.join(events)}")
        print(f"  - 无变化(N): {events.count('N')}")
        print(f"  - 错配(M): {events.count('M')}")
        print(f"  - 删除(D): {events.count('D')}")
        print(f"  - 插入(I): {events.count('I')}")
        print(f"  - 空(E): {events.count('E')}")
    
    print("✅ 标准化效果演示完成\n")


def main():
    """主演示函数"""
    print("🧬 CARLIN序列标准化功能演示")
    print("=" * 60)
    print()
    
    demo_prefix_postfix_sanitization()
    demo_conserved_region_sanitization() 
    demo_cutsite_preservation()
    demo_batch_sanitization()
    demo_sanitization_effects()
    
    print("🎉 CARLIN序列标准化演示完成！")
    print()
    print("💡 标准化的关键作用:")
    print("   ✅ 移除技术人工制品 (前后缀污染)")
    print("   ✅ 修正保守区域测序错误")
    print("   ✅ 保留真实的生物学编辑")
    print("   ✅ 提高下游分析准确性")
    print("   ✅ 标准化序列表示")
    print()
    print("🔬 生物学意义:")
    print("   - PAM和consite区域在生物学上高度保守")
    print("   - 这些区域的变化更可能是技术错误")
    print("   - Cutsite区域是CRISPR的真实编辑目标")
    print("   - 标准化确保分析专注于功能性变化")


if __name__ == "__main__":
    main() 