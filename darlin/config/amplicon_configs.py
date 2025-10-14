#!/usr/bin/env python3
"""
CARLIN扩增子配置模块

实现CARLIN特异性的序列配置和惩罚参数管理
"""

import json
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
from dataclasses import dataclass


@dataclass
class SequenceConfig:
    """序列配置数据类"""
    segments: List[str]
    pam: str
    prefix: str
    postfix: str
    primer5: str
    primer3: str
    secondary_sequence: str
    
    def __post_init__(self):
        """验证配置参数"""
        if len(self.segments) != 10:
            raise ValueError(f"CARLIN必须有10个segments，但得到了{len(self.segments)}个")
        
        # 验证每个segment长度为20bp
        for i, seg in enumerate(self.segments):
            if len(seg) != 20:
                raise ValueError(f"Segment {i+1} 长度应为20bp，但得到了{len(seg)}bp")


@dataclass 
class PenaltyConfig:
    """惩罚参数配置数据类"""
    init: float
    cutsites: List[float]
    consites: List[float]
    pam: List[float]
    prefix: List[float]
    postfix: List[float]
    
    def __post_init__(self):
        """验证惩罚参数"""
        if len(self.cutsites) != 7:
            raise ValueError(f"cutsites惩罚应有7个值，得到{len(self.cutsites)}个")
        if len(self.consites) != 13:
            raise ValueError(f"consites惩罚应有13个值，得到{len(self.consites)}个")
        if len(self.pam) != 7:
            raise ValueError(f"pam惩罚应有7个值，得到{len(self.pam)}个")


class AmpliconConfig:
    """
    CARLIN扩增子配置类
    
    管理CARLIN序列结构、评分矩阵和位置特异性惩罚参数
    """
    
    def __init__(self, config_data: Optional[Dict] = None, config_file: Optional[str] = None):
        """
        初始化配置
        
        Args:
            config_data: 配置字典数据
            config_file: JSON配置文件路径
        """
        if config_file:
            self._load_from_file(config_file)
        elif config_data:
            self._load_from_dict(config_data)
        else:
            raise ValueError("必须提供config_data或config_file")
        
        # 构建完整的CARLIN序列
        self._build_carlin_sequence()
        
        # 构建位置特异性惩罚数组
        self._build_penalty_arrays()
    
    def _load_from_file(self, config_file: str):
        """从JSON文件加载配置"""
        config_path = Path(config_file)
        if not config_path.exists():
            raise FileNotFoundError(f"配置文件不存在: {config_file}")
        
        with open(config_path, 'r', encoding='utf-8') as f:
            config_data = json.load(f)
        
        self._load_from_dict(config_data)
    
    def _load_from_dict(self, config_data: Dict):
        """从字典加载配置"""
        # 解析序列配置
        seq_data = config_data['sequence']
        self.sequence = SequenceConfig(
            segments=seq_data['segments'],
            pam=seq_data['pam'],
            prefix=seq_data['prefix'],
            postfix=seq_data['postfix'],
            primer5=seq_data['Primer5'],
            primer3=seq_data['Primer3'],
            secondary_sequence=seq_data['SecondarySequence']
        )
        
        # 解析匹配得分
        self.match_scores = config_data.get('match_score', {})
        
        # 解析惩罚配置
        open_data = config_data['open_penalty']
        self.open_penalty = PenaltyConfig(
            init=open_data['init'],
            cutsites=open_data['cutsites'],
            consites=open_data['consites'],
            pam=open_data['pam'],
            prefix=open_data['prefix'],
            postfix=open_data['postfix']
        )
        
        close_data = config_data['close_penalty']
        self.close_penalty = PenaltyConfig(
            init=close_data.get('init', open_data['init']),  # 如果没有指定，使用open_penalty的init
            cutsites=close_data['cutsites'],
            consites=close_data['consites'],
            pam=close_data['pam'],
            prefix=close_data['prefix'],
            postfix=close_data['postfix']
        )
    
    def _build_carlin_sequence(self):
        """构建完整的CARLIN序列结构"""
        # CARLIN内部结构: prefix + (segment + pam) x 10 - 最后一个pam + postfix
        carlin_parts = [self.sequence.prefix]
        
        # 添加10个 (segment + pam)，最后一个segment不加pam
        for i, segment in enumerate(self.sequence.segments):
            carlin_parts.append(segment)
            if i < 9:  # 前9个segment后面加pam
                carlin_parts.append(self.sequence.pam)
        
        # 添加postfix
        carlin_parts.append(self.sequence.postfix)
        
        # 组装完整序列
        self.carlin_sequence = ''.join(carlin_parts)
        
        # 完整的扩增子序列: Primer5 + CARLIN + SecondarySequence + Primer3
        self.full_sequence = (
            self.sequence.primer5 + 
            self.carlin_sequence + 
            self.sequence.secondary_sequence + 
            self.sequence.primer3
        )
        
        # 计算各部分的位置信息
        self._calculate_positions()
    
    def _calculate_positions(self):
        """计算序列各部分的位置信息"""
        self.positions = {}
        
        # Primer5位置
        primer5_start = 0
        primer5_end = len(self.sequence.primer5)
        self.positions['primer5'] = (primer5_start, primer5_end)
        
        # CARLIN内部位置 (相对于CARLIN序列)
        carlin_start = primer5_end
        
        # Prefix位置
        prefix_start = 0
        prefix_end = len(self.sequence.prefix)
        self.positions['prefix'] = (prefix_start, prefix_end)
        
        # Segments和PAMs的位置
        self.positions['segments'] = []
        self.positions['pams'] = []
        self.positions['cutsites'] = []
        self.positions['consites'] = []
        
        current_pos = prefix_end
        
        for i, segment in enumerate(self.sequence.segments):
            # Segment位置
            segment_start = current_pos
            segment_end = current_pos + len(segment)
            self.positions['segments'].append((segment_start, segment_end))
            
            # Consite (前13bp) 和 Cutsite (后7bp)
            consite_start = segment_start
            consite_end = segment_start + 13
            cutsite_start = consite_end
            cutsite_end = segment_end
            
            self.positions['consites'].append((consite_start, consite_end))
            self.positions['cutsites'].append((cutsite_start, cutsite_end))
            
            current_pos = segment_end
            
            # PAM位置 (前9个segment后面有PAM)
            if i < 9:
                pam_start = current_pos
                pam_end = current_pos + len(self.sequence.pam)
                self.positions['pams'].append((pam_start, pam_end))
                current_pos = pam_end
        
        # Postfix位置
        postfix_start = current_pos
        postfix_end = current_pos + len(self.sequence.postfix)
        self.positions['postfix'] = (postfix_start, postfix_end)
        
        # 相对于完整序列的位置
        self.positions['carlin'] = (carlin_start, carlin_start + len(self.carlin_sequence))
        
        # Secondary sequence位置
        secondary_start = carlin_start + len(self.carlin_sequence)
        secondary_end = secondary_start + len(self.sequence.secondary_sequence)
        self.positions['secondary'] = (secondary_start, secondary_end)
        
        # Primer3位置
        primer3_start = secondary_end
        primer3_end = primer3_start + len(self.sequence.primer3)
        self.positions['primer3'] = (primer3_start, primer3_end)
    
    def _build_penalty_arrays(self):
        """构建位置特异性的惩罚数组"""
        carlin_length = len(self.carlin_sequence)
        
        # 初始化惩罚数组
        self.open_penalty_array = np.full(carlin_length + 1, self.open_penalty.init, dtype=np.float64)
        self.close_penalty_array = np.full(carlin_length + 1, self.close_penalty.init, dtype=np.float64)
        
        # 设置prefix惩罚
        prefix_start, prefix_end = self.positions['prefix']
        self._set_penalty_values(self.open_penalty_array, prefix_start, prefix_end, self.open_penalty.prefix)
        self._set_penalty_values(self.close_penalty_array, prefix_start, prefix_end, self.close_penalty.prefix)
        
        # 设置各segment的consite和cutsite惩罚
        for i in range(10):
            # Consite惩罚
            consite_start, consite_end = self.positions['consites'][i]
            self._set_penalty_values(self.open_penalty_array, consite_start, consite_end, self.open_penalty.consites)
            self._set_penalty_values(self.close_penalty_array, consite_start, consite_end, self.close_penalty.consites)
            
            # Cutsite惩罚
            cutsite_start, cutsite_end = self.positions['cutsites'][i]
            self._set_penalty_values(self.open_penalty_array, cutsite_start, cutsite_end, self.open_penalty.cutsites)
            self._set_penalty_values(self.close_penalty_array, cutsite_start, cutsite_end, self.close_penalty.cutsites)
            
            # PAM惩罚 (前9个)
            if i < 9:
                pam_start, pam_end = self.positions['pams'][i]
                self._set_penalty_values(self.open_penalty_array, pam_start, pam_end, self.open_penalty.pam)
                self._set_penalty_values(self.close_penalty_array, pam_start, pam_end, self.close_penalty.pam)
        
        # 设置postfix惩罚
        postfix_start, postfix_end = self.positions['postfix']
        self._set_penalty_values(self.open_penalty_array, postfix_start, postfix_end, self.open_penalty.postfix)
        self._set_penalty_values(self.close_penalty_array, postfix_start, postfix_end, self.close_penalty.postfix)
    
    def _set_penalty_values(self, penalty_array: np.ndarray, start: int, end: int, values: List[float]):
        """设置指定区域的惩罚值"""
        region_length = end - start
        if len(values) != region_length:
            raise ValueError(f"惩罚值数量({len(values)})与区域长度({region_length})不匹配")
        
        for i, value in enumerate(values):
            penalty_array[start + i] = value
    
    def get_reference_sequence(self) -> str:
        """获取参考序列（仅CARLIN部分）"""
        return self.carlin_sequence
    
    def get_full_reference_sequence(self) -> str:
        """获取完整的参考序列"""
        return self.full_sequence
    
    def get_penalty_arrays(self) -> Tuple[np.ndarray, np.ndarray]:
        """获取位置特异性惩罚数组"""
        return self.open_penalty_array.copy(), self.close_penalty_array.copy()
    
    def get_motif_info(self, position: int) -> Dict[str, any]:
        """获取指定位置的motif信息"""
        motif_info = {
            'type': 'unknown',
            'motif_id': None,
            'position_in_motif': None
        }
        
        # 检查各个区域
        prefix_start, prefix_end = self.positions['prefix']
        if prefix_start <= position < prefix_end:
            motif_info.update({
                'type': 'prefix',
                'motif_id': 0,
                'position_in_motif': position - prefix_start
            })
            return motif_info
        
        # 检查segments
        for i, (start, end) in enumerate(self.positions['segments']):
            if start <= position < end:
                # 判断是consite还是cutsite
                consite_start, consite_end = self.positions['consites'][i]
                if consite_start <= position < consite_end:
                    motif_info.update({
                        'type': 'consite',
                        'motif_id': i,
                        'position_in_motif': position - consite_start
                    })
                else:
                    cutsite_start, _ = self.positions['cutsites'][i]
                    motif_info.update({
                        'type': 'cutsite',
                        'motif_id': i,
                        'position_in_motif': position - cutsite_start
                    })
                return motif_info
        
        # 检查PAMs
        for i, (start, end) in enumerate(self.positions['pams']):
            if start <= position < end:
                motif_info.update({
                    'type': 'pam',
                    'motif_id': i,
                    'position_in_motif': position - start
                })
                return motif_info
        
        # 检查postfix
        postfix_start, postfix_end = self.positions['postfix']
        if postfix_start <= position < postfix_end:
            motif_info.update({
                'type': 'postfix',
                'motif_id': 0,
                'position_in_motif': position - postfix_start
            })
        
        return motif_info
    
    def summary(self) -> str:
        """返回配置摘要信息"""
        lines = [
            "=== CARLIN扩增子配置摘要 ===",
            f"完整序列长度: {len(self.full_sequence)} bp",
            f"CARLIN序列长度: {len(self.carlin_sequence)} bp",
            f"Segments数量: {len(self.sequence.segments)}",
            f"PAM序列: {self.sequence.pam}",
            f"Prefix: {self.sequence.prefix}",
            f"Postfix: {self.sequence.postfix}",
            "",
            "=== 序列结构 ===",
            f"Primer5 ({len(self.sequence.primer5)} bp): {self.sequence.primer5}",
            f"CARLIN ({len(self.carlin_sequence)} bp): {self.carlin_sequence[:50]}...",
            f"Secondary ({len(self.sequence.secondary_sequence)} bp): {self.sequence.secondary_sequence}",
            f"Primer3 ({len(self.sequence.primer3)} bp): {self.sequence.primer3}",
            "",
            "=== 惩罚参数范围 ===",
            f"Open penalty: {self.open_penalty_array.min():.2f} - {self.open_penalty_array.max():.2f}",
            f"Close penalty: {self.close_penalty_array.min():.2f} - {self.close_penalty_array.max():.2f}"
        ]
        return "\n".join(lines)


def load_original_carlin_config() -> AmpliconConfig:
    """加载原始CARLIN配置"""
    # 自动查找配置文件路径
    possible_paths = [
        "../allele_calling/carlin-master/cfg/amplicon/OriginalCARLIN.json",
        "../../allele_calling/carlin-master/cfg/amplicon/OriginalCARLIN.json", 
        "/home/jarning/project-core/DARLIN-toolkits/allele_calling/carlin-master/cfg/amplicon/OriginalCARLIN.json"
    ]
    
    for path in possible_paths:
        if Path(path).exists():
            return AmpliconConfig(config_file=path)
    
    raise FileNotFoundError("找不到OriginalCARLIN.json配置文件")


# 预定义的配置
ORIGINAL_CARLIN = None  # 将在首次访问时延迟加载

def get_original_carlin_config() -> AmpliconConfig:
    """获取原始CARLIN配置（延迟加载）"""
    global ORIGINAL_CARLIN
    if ORIGINAL_CARLIN is None:
        ORIGINAL_CARLIN = load_original_carlin_config()
    return ORIGINAL_CARLIN

# 立即初始化配置
try:
    ORIGINAL_CARLIN = get_original_carlin_config()
except FileNotFoundError:
    # 如果找不到配置文件，使用硬编码的默认配置
    ORIGINAL_CARLIN = AmpliconConfig() 

def load_carlin_config_by_locus(locus: str = "Col1a1") -> AmpliconConfig:
    """
    根据locus名称加载对应的CARLIN配置
    
    Args:
        locus: 位点名称，支持 "Col1a1", "Rosa", "Tigre"
    
    Returns:
        AmpliconConfig: 对应的配置对象
    """
    # 验证locus参数
    valid_loci = ["Col1a1", "Rosa", "Tigre"]
    if locus not in valid_loci:
        raise ValueError(f"不支持的locus: {locus}。支持的locus: {valid_loci}")
    
    # 构建配置文件路径
    config_dir = Path(__file__).parent / "data"
    config_file = config_dir / f"array_{locus}.json"
    
    if not config_file.exists():
        raise FileNotFoundError(f"找不到配置文件: {config_file}")
    
    return AmpliconConfig(config_file=str(config_file)) 