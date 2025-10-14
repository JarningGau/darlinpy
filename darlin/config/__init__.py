"""
配置模块 - 包含CARLIN扩增子配置和相关常量
"""

from .amplicon_configs import AmpliconConfig, get_original_carlin_config, load_carlin_config_by_locus
from .scoring_matrices import ScoringConfig, get_default_scoring_config, create_nuc44_matrix

__all__ = [
    'AmpliconConfig',
    'get_original_carlin_config',
    'load_carlin_config_by_locus',
    'ScoringConfig', 
    'get_default_scoring_config',
    'create_nuc44_matrix'
] 