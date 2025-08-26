"""
突变注释模块 - 实现CARLIN突变事件的识别和注释
"""

from .mutation import (
    Mutation,
    MutationType,
    MutationIdentifier,
    annotate_mutations
)

__all__ = [
    'Mutation',
    'MutationType', 
    'MutationIdentifier',
    'annotate_mutations'
] 