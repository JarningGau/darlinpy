"""
DARLIN Python - CARLIN序列分析工具的Python实现

主要功能：
- CRISPR-Cas9序列比对
- 等位基因调用
- 突变注释
- 高级序列分析API

Author: DARLIN-toolkits Team
"""

__version__ = "0.1.0"
__author__ = "DARLIN-toolkits Team"

# 导出主要API接口
from .api import analyze_sequences, quick_analyze, AnalysisResult

# 导出核心功能模块
from .alignment.cas9_align import cas9_align, nt2int, int2nt
from .mutations import Mutation, MutationType, annotate_mutations

# 导出配置
from .config.amplicon_configs import ORIGINAL_CARLIN

__all__ = [
    '__version__',
    '__author__',
    # 主API
    'analyze_sequences',
    'quick_analyze', 
    'AnalysisResult',
    # 核心功能
    'cas9_align',
    'nt2int',
    'int2nt',
    'Mutation',
    'MutationType',
    'annotate_mutations',
    # 配置
    'ORIGINAL_CARLIN'
] 