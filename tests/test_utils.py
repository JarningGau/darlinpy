#!/usr/bin/env python3
"""
序列处理工具函数测试
"""

import pytest
from darlin.utils.sequence_utils import (
    reverse_complement,
    is_valid_dna,
    clean_sequence,
    calculate_gc_content,
    count_nucleotides,
    hamming_distance,
    format_sequence,
    parse_fasta_string,
    to_fasta_string
)


class TestReverseComplement:
    """测试反向互补功能"""
    
    def test_basic_reverse_complement(self):
        """测试基本反向互补"""
        assert reverse_complement("ATCG") == "CGAT"
        assert reverse_complement("AAATTTGGGCCC") == "GGGCCCAAATTT"
        assert reverse_complement("GCTAGCTA") == "TAGCTAGC"
    
    def test_empty_sequence(self):
        """测试空序列"""
        assert reverse_complement("") == ""
    
    def test_ambiguous_bases(self):
        """测试模糊碱基"""
        assert reverse_complement("ATCGRYSWKM") == "KMWSRYCGAT"
    
    def test_invalid_sequence(self):
        """测试无效序列"""
        with pytest.raises(ValueError):
            reverse_complement("ATCGXYZ")


class TestSequenceValidation:
    """测试序列验证功能"""
    
    def test_valid_dna(self):
        """测试有效DNA序列"""
        assert is_valid_dna("ATCG") == True
        assert is_valid_dna("AAATTTGGGCCC") == True
        assert is_valid_dna("") == True  # 空序列被认为是有效的
    
    def test_invalid_dna(self):
        """测试无效DNA序列"""
        assert is_valid_dna("ATCGXYZ") == False
        assert is_valid_dna("ATCG123") == False
    
    def test_ambiguous_bases(self):
        """测试模糊碱基"""
        assert is_valid_dna("ATCGRYSWKM", allow_ambiguous=True) == True
        assert is_valid_dna("ATCGRYSWKM", allow_ambiguous=False) == False


class TestSequenceCleaning:
    """测试序列清理功能"""
    
    def test_clean_sequence(self):
        """测试序列清理"""
        assert clean_sequence("AT-CG") == "ATCG"
        assert clean_sequence("at cg\n") == "ATCG"
        assert clean_sequence("A\tT\nC\tG") == "ATCG"
    
    def test_custom_remove_chars(self):
        """测试自定义移除字符"""
        assert clean_sequence("A@T#C$G", remove_chars="@#$") == "ATCG"


class TestGCContent:
    """测试GC含量计算"""
    
    def test_calculate_gc_content(self):
        """测试GC含量计算"""
        assert calculate_gc_content("ATCG") == 50.0  # 2 GC out of 4
        assert calculate_gc_content("AAAA") == 0.0   # No GC
        assert calculate_gc_content("GGGG") == 100.0 # All GC
    
    def test_empty_sequence(self):
        """测试空序列"""
        assert calculate_gc_content("") == 0.0


class TestNucleotideCounting:
    """测试核苷酸计数"""
    
    def test_count_nucleotides(self):
        """测试核苷酸计数"""
        counts = count_nucleotides("ATCGATCG")
        assert counts['A'] == 2
        assert counts['T'] == 2  
        assert counts['C'] == 2
        assert counts['G'] == 2
    
    def test_empty_sequence(self):
        """测试空序列"""
        counts = count_nucleotides("")
        assert all(count == 0 for count in counts.values())


class TestHammingDistance:
    """测试汉明距离计算"""
    
    def test_hamming_distance(self):
        """测试汉明距离"""
        assert hamming_distance("ATCG", "ATCG") == 0  # 完全相同
        assert hamming_distance("ATCG", "TTCG") == 1  # 1个差异
        assert hamming_distance("ATCG", "TTGG") == 2  # 2个差异
    
    def test_unequal_length(self):
        """测试不等长序列"""
        with pytest.raises(ValueError):
            hamming_distance("ATCG", "ATCGG")


class TestSequenceFormatting:
    """测试序列格式化"""
    
    def test_format_sequence(self):
        """测试序列格式化"""
        seq = "ATCGATCGATCGATCG"
        formatted = format_sequence(seq, line_length=8)
        
        lines = formatted.split('\n')
        assert len(lines) == 2
        assert lines[0] == "ATCGATCG"
        assert lines[1] == "ATCGATCG"
    
    def test_format_with_numbers(self):
        """测试带行号的格式化"""
        seq = "ATCGATCGATCGATCG"
        formatted = format_sequence(seq, line_length=8, add_numbers=True)
        
        lines = formatted.split('\n')
        assert "1 ATCGATCG" in lines[0]
        assert "9 ATCGATCG" in lines[1]


class TestFASTAHandling:
    """测试FASTA格式处理"""
    
    def test_parse_fasta_string(self):
        """测试FASTA字符串解析"""
        fasta = """>seq1
ATCGATCG
TTGGCCAA
>seq2
GGGGCCCC"""
        
        sequences = parse_fasta_string(fasta)
        
        assert len(sequences) == 2
        assert sequences[0] == ("seq1", "ATCGATCGTTGGCCAA")
        assert sequences[1] == ("seq2", "GGGGCCCC")
    
    def test_to_fasta_string(self):
        """测试FASTA字符串生成"""
        sequences = [
            ("seq1", "ATCGATCGATCGATCG"),
            ("seq2", "GGGGCCCCTTTTAAAA")
        ]
        
        fasta = to_fasta_string(sequences, line_length=8)
        
        assert ">seq1" in fasta
        assert ">seq2" in fasta
        assert "ATCGATCG" in fasta
        assert "GGGGCCCC" in fasta
    
    def test_empty_fasta(self):
        """测试空FASTA"""
        assert parse_fasta_string("") == []


if __name__ == "__main__":
    pytest.main([__file__]) 