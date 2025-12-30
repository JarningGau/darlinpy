"""
Utility modules
"""

from .build_config import build_carlin_config
from .misc import open_fastq_file, iter_fastq_single
__all__ = [
    'build_carlin_config',
    'open_fastq_file',
    'iter_fastq_single'
]