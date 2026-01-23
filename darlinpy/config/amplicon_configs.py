#!/usr/bin/env python3
"""
CARLIN amplicon configuration module

Implements CARLIN-specific sequence configuration and penalty parameter management
"""

import json
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
from dataclasses import dataclass


@dataclass
class SequenceConfig:
    """Sequence configuration data class"""
    segments: List[str]
    pam: str
    prefix: str
    postfix: str
    primer5: str
    primer3: str
    secondary_sequence: str
    
    def __post_init__(self):
        """Validate configuration parameters"""
        if len(self.segments) != 10:
            raise ValueError(f"CARLIN must have 10 segments, but got {len(self.segments)}")
        
        # Validate each segment length is 20bp
        for i, seg in enumerate(self.segments):
            if len(seg) != 20:
                raise ValueError(f"Segment {i+1} should be 20bp long, but got {len(seg)}bp")


@dataclass 
class PenaltyConfig:
    """Penalty parameter configuration data class"""
    init: float
    cutsites: List[float]
    consites: List[float]
    pam: List[float]
    prefix: List[float]
    postfix: List[float]
    
    def __post_init__(self):
        """Validate penalty parameters"""
        if len(self.cutsites) != 7:
            raise ValueError(f"cutsites penalties should have 7 values, got {len(self.cutsites)}")
        if len(self.consites) != 13:
            raise ValueError(f"consites penalties should have 13 values, got {len(self.consites)}")
        if len(self.pam) != 7:
            raise ValueError(f"pam penalties should have 7 values, got {len(self.pam)}")


class AmpliconConfig:
    """
    CARLIN amplicon configuration class
    
    Manages CARLIN sequence structure, scoring matrices and position-specific penalty parameters
    """
    
    def __init__(self, config_data: Optional[Dict] = None, config_file: Optional[str] = None):
        """
        Initialize configuration
        
        Args:
            config_data: Configuration dictionary data
            config_file: JSON configuration file path
        """
        if config_file:
            self._load_from_file(config_file)
        elif config_data:
            self._load_from_dict(config_data)
        else:
            raise ValueError("Must provide config_data or config_file")
        
        # Build complete CARLIN sequence
        self._build_carlin_sequence()
        
        # Build position-specific penalty arrays
        self._build_penalty_arrays()
    
    def _load_from_file(self, config_file: str):
        """Load configuration from JSON file"""
        config_path = Path(config_file)
        if not config_path.exists():
            raise FileNotFoundError(f"Configuration file does not exist: {config_file}")
        
        with open(config_path, 'r', encoding='utf-8') as f:
            config_data = json.load(f)
        
        self._load_from_dict(config_data)
    
    def _load_from_dict(self, config_data: Dict):
        """Load configuration from dictionary"""
        # Parse sequence configuration
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
        
        # Parse match scores
        self.match_scores = config_data.get('match_score', {})
        
        # Parse penalty configuration
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
            init=close_data.get('init', open_data['init']),  # If not specified, use open_penalty's init
            cutsites=close_data['cutsites'],
            consites=close_data['consites'],
            pam=close_data['pam'],
            prefix=close_data['prefix'],
            postfix=close_data['postfix']
        )
    
    def _build_carlin_sequence(self):
        """Build complete CARLIN sequence structure"""
        # CARLIN internal structure: prefix + (segment + pam) x 10 - last pam + postfix
        carlin_parts = [self.sequence.prefix]
        
        # Add 10 (segment + pam), last segment without pam
        for i, segment in enumerate(self.sequence.segments):
            carlin_parts.append(segment)
            if i < 9:  # First 9 segments followed by pam
                carlin_parts.append(self.sequence.pam)
        
        # Add postfix
        carlin_parts.append(self.sequence.postfix)
        
        # Assemble complete sequence
        self.carlin_sequence = ''.join(carlin_parts)
        
        # Complete amplicon sequence: Primer5 + CARLIN + SecondarySequence + Primer3
        self.full_sequence = (
            self.sequence.primer5 + 
            self.carlin_sequence + 
            self.sequence.secondary_sequence + 
            self.sequence.primer3
        )
        
        # Calculate position information for each part
        self._calculate_positions()
    
    def _calculate_positions(self):
        """Calculate position information for sequence parts"""
        self.positions = {}
        
        # Primer5 position
        primer5_start = 0
        primer5_end = len(self.sequence.primer5)
        self.positions['primer5'] = (primer5_start, primer5_end)
        
        # CARLIN internal positions (relative to CARLIN sequence)
        carlin_start = primer5_end
        
        # Prefix position
        prefix_start = 0
        prefix_end = len(self.sequence.prefix)
        self.positions['prefix'] = (prefix_start, prefix_end)
        
        # Segments and PAMs positions
        self.positions['segments'] = []
        self.positions['pams'] = []
        self.positions['cutsites'] = []
        self.positions['consites'] = []
        
        current_pos = prefix_end
        
        for i, segment in enumerate(self.sequence.segments):
            # Segment position
            segment_start = current_pos
            segment_end = current_pos + len(segment)
            self.positions['segments'].append((segment_start, segment_end))
            
            # Consite (first 13bp) and Cutsite (last 7bp)
            consite_start = segment_start
            consite_end = segment_start + 13
            cutsite_start = consite_end
            cutsite_end = segment_end
            
            self.positions['consites'].append((consite_start, consite_end))
            self.positions['cutsites'].append((cutsite_start, cutsite_end))
            
            current_pos = segment_end
            
            # PAM position (first 9 segments followed by PAM)
            if i < 9:
                pam_start = current_pos
                pam_end = current_pos + len(self.sequence.pam)
                self.positions['pams'].append((pam_start, pam_end))
                current_pos = pam_end
        
        # Postfix position
        postfix_start = current_pos
        postfix_end = current_pos + len(self.sequence.postfix)
        self.positions['postfix'] = (postfix_start, postfix_end)
        
        # Positions relative to complete sequence
        self.positions['carlin'] = (carlin_start, carlin_start + len(self.carlin_sequence))
        
        # Secondary sequence position
        secondary_start = carlin_start + len(self.carlin_sequence)
        secondary_end = secondary_start + len(self.sequence.secondary_sequence)
        self.positions['secondary'] = (secondary_start, secondary_end)
        
        # Primer3 position
        primer3_start = secondary_end
        primer3_end = primer3_start + len(self.sequence.primer3)
        self.positions['primer3'] = (primer3_start, primer3_end)
    
    def _build_penalty_arrays(self):
        """Build position-specific penalty arrays"""
        carlin_length = len(self.carlin_sequence)
        
        # Initialize penalty arrays
        self.open_penalty_array = np.full(carlin_length + 1, self.open_penalty.init, dtype=np.float64)
        self.close_penalty_array = np.full(carlin_length + 1, self.close_penalty.init, dtype=np.float64)
        
        # Set prefix penalties
        prefix_start, prefix_end = self.positions['prefix']
        self._set_penalty_values(self.open_penalty_array, prefix_start, prefix_end, self.open_penalty.prefix)
        self._set_penalty_values(self.close_penalty_array, prefix_start, prefix_end, self.close_penalty.prefix)
        
        # Set consite and cutsite penalties for each segment
        for i in range(10):
            # Consite penalties
            consite_start, consite_end = self.positions['consites'][i]
            self._set_penalty_values(self.open_penalty_array, consite_start, consite_end, self.open_penalty.consites)
            self._set_penalty_values(self.close_penalty_array, consite_start, consite_end, self.close_penalty.consites)
            
            # Cutsite penalties
            cutsite_start, cutsite_end = self.positions['cutsites'][i]
            self._set_penalty_values(self.open_penalty_array, cutsite_start, cutsite_end, self.open_penalty.cutsites)
            self._set_penalty_values(self.close_penalty_array, cutsite_start, cutsite_end, self.close_penalty.cutsites)
            
            # PAM penalties (first 9)
            if i < 9:
                pam_start, pam_end = self.positions['pams'][i]
                self._set_penalty_values(self.open_penalty_array, pam_start, pam_end, self.open_penalty.pam)
                self._set_penalty_values(self.close_penalty_array, pam_start, pam_end, self.close_penalty.pam)
        
        # Set postfix penalties
        postfix_start, postfix_end = self.positions['postfix']
        self._set_penalty_values(self.open_penalty_array, postfix_start, postfix_end, self.open_penalty.postfix)
        self._set_penalty_values(self.close_penalty_array, postfix_start, postfix_end, self.close_penalty.postfix)
    
    def _set_penalty_values(self, penalty_array: np.ndarray, start: int, end: int, values: List[float]):
        """Set penalty values for specified region"""
        region_length = end - start
        if len(values) != region_length:
            raise ValueError(f"Number of penalty values ({len(values)}) does not match region length ({region_length})")
        
        for i, value in enumerate(values):
            penalty_array[start + i] = value
    
    def get_reference_sequence(self) -> str:
        """Get reference sequence (CARLIN part only)"""
        return self.carlin_sequence
    
    def get_full_reference_sequence(self) -> str:
        """Get complete reference sequence"""
        return self.full_sequence
    
    def get_penalty_arrays(self) -> Tuple[np.ndarray, np.ndarray]:
        """Get position-specific penalty arrays"""
        return self.open_penalty_array.copy(), self.close_penalty_array.copy()
    
    def get_motif_info(self, position: int) -> Dict[str, any]:
        """Get motif information for specified position"""
        motif_info = {
            'type': 'unknown',
            'motif_id': None,
            'position_in_motif': None
        }
        
        # Check each region
        prefix_start, prefix_end = self.positions['prefix']
        if prefix_start <= position < prefix_end:
            motif_info.update({
                'type': 'prefix',
                'motif_id': 0,
                'position_in_motif': position - prefix_start
            })
            return motif_info
        
        # Check segments
        for i, (start, end) in enumerate(self.positions['segments']):
            if start <= position < end:
                # Determine if consite or cutsite
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
        
        # Check PAMs
        for i, (start, end) in enumerate(self.positions['pams']):
            if start <= position < end:
                motif_info.update({
                    'type': 'pam',
                    'motif_id': i,
                    'position_in_motif': position - start
                })
                return motif_info
        
        # Check postfix
        postfix_start, postfix_end = self.positions['postfix']
        if postfix_start <= position < postfix_end:
            motif_info.update({
                'type': 'postfix',
                'motif_id': 0,
                'position_in_motif': position - postfix_start
            })
        
        return motif_info
    
    def summary(self) -> str:
        """Return configuration summary information"""
        lines = [
            "=== CARLIN Amplicon Configuration Summary ===",
            f"Complete sequence length: {len(self.full_sequence)} bp",
            f"CARLIN sequence length: {len(self.carlin_sequence)} bp",
            f"Number of segments: {len(self.sequence.segments)}",
            f"PAM sequence: {self.sequence.pam}",
            f"Prefix: {self.sequence.prefix}",
            f"Postfix: {self.sequence.postfix}",
            "",
            "=== Sequence Structure ===",
            f"Primer5 ({len(self.sequence.primer5)} bp): {self.sequence.primer5}",
            f"CARLIN ({len(self.carlin_sequence)} bp): {self.carlin_sequence[:50]}...",
            f"Secondary ({len(self.sequence.secondary_sequence)} bp): {self.sequence.secondary_sequence}",
            f"Primer3 ({len(self.sequence.primer3)} bp): {self.sequence.primer3}",
            "",
            "=== Penalty Parameter Range ===",
            f"Open penalty: {self.open_penalty_array.min():.2f} - {self.open_penalty_array.max():.2f}",
            f"Close penalty: {self.close_penalty_array.min():.2f} - {self.close_penalty_array.max():.2f}"
        ]
        return "\n".join(lines)


def load_carlin_config_by_locus(locus: str = "Col1a1") -> AmpliconConfig:
    """
    Load corresponding CARLIN configuration by locus name
    
    Args:
        locus: Locus name, supports "Col1a1", "Rosa", "Tigre"
    
    Returns:
        AmpliconConfig: Corresponding configuration object
    """
    # Validate locus parameter
    valid_loci = ["Col1a1", "Rosa", "Tigre"]
    if locus not in valid_loci:
        raise ValueError(f"Unsupported locus: {locus}. Supported loci: {valid_loci}")
    
    # Build configuration file path
    config_dir = Path(__file__).parent / "data"
    config_file = config_dir / f"array_{locus}.json"
    
    if not config_file.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_file}")
    
    return AmpliconConfig(config_file=str(config_file)) 


# Predefined configurations
ORIGINAL_CARLIN = None  # Will be lazily loaded on first access

def get_original_carlin_config() -> AmpliconConfig:
    """Get original CARLIN configuration (lazy loading)"""
    global ORIGINAL_CARLIN
    if ORIGINAL_CARLIN is None:
        ORIGINAL_CARLIN = load_carlin_config_by_locus("Col1a1")
    return ORIGINAL_CARLIN

# Initialize configuration immediately
try:
    ORIGINAL_CARLIN = get_original_carlin_config()
except FileNotFoundError:
    # If configuration file not found, use hardcoded default configuration
    ORIGINAL_CARLIN = AmpliconConfig()
