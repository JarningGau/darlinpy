"""
CARLIN等位基因调用模块

这个模块提供了完整的等位基因调用功能，包括：
- 精确和粗粒度调用算法
- 等位基因调用结果数据结构
- 统计分析功能
"""

from .allele_caller import AlleleCaller
from .allele_data import (
    AlleleCallResult, 
    BulkAlleleCallResult, 
    AlleleCallStatistics
)

__all__ = [
    'AlleleCaller',
    'AlleleCallResult',
    'BulkAlleleCallResult', 
    'AlleleCallStatistics'
] 