# 更新日志

## [Unreleased] - 待发布

### 新增功能
- **突变注释模块** (`mutations/`)
  - 全新的`Mutation`类，支持多种突变类型和HGVS格式注释
  - `MutationIdentifier`类，实现智能突变识别和分析
  - Cas9特异性事件识别和置信度评估
  - 相邻突变智能合并算法
  - 15个全面的单元测试用例
  - 交互式演示脚本 (`examples/mutation_demo.py`)

### 改进
- 完善了AlignedSEQ与突变注释系统的集成
- 增强了错误处理和数据验证

---

## [0.1.0] - 2024-12-19

### 主要功能

#### 核心比对算法
- **cas9_align**: 实现了三状态动态规划算法，支持位置特异性gap惩罚
- **序列编码**: `nt2int`和`int2nt`函数，支持高效的数值计算
- **CARLINAligner**: 高级比对器，提供批量处理和统计分析功能

#### 配置管理

### Added
- Initial release of DARLIN Python
- Core CARLIN sequence alignment functionality
- Position-specific gap penalty system
- NUC44 and multiple scoring matrices support
- CARLIN amplicon configuration system
- Sequence sanitization and standardization
- AlignedSEQ and AlignedSEQMotif data structures
- Comprehensive test suite
- Batch sequence processing
- Detailed alignment statistics and analysis

### Core Features
- **cas9_align algorithm**: High-precision sequence alignment with position-specific penalties
- **AmpliconConfig**: Complete CARLIN amplicon configuration management
- **ScoringConfig**: Multiple scoring matrices (NUC44, simple, transition/transversion)
- **CARLINAligner**: Integrated high-level alignment interface
- **SequenceSanitizer**: Two-stage sequence standardization
  - Prefix/postfix contamination removal
  - Conserved region error correction
- **AlignedSEQ system**: Motif-level sequence analysis and event classification

### Technical Implementation
- Position-specific gap penalties based on biological motifs
- Automatic motif boundary calculation
- Robust error handling and input validation
- Comprehensive documentation and examples
- Development environment setup scripts

### Testing and Quality
- Unit tests for all core components
- Integration tests for end-to-end workflows
- Sanitization functionality tests
- Performance benchmarks
- Code examples and demonstrations

### Documentation
- Complete API documentation
- Usage examples and tutorials
- Development setup guide
- Biological background and methodology

### Compatibility
- Python 3.8+
- NumPy 1.20.0+
- SciPy 1.7.0+
- BioPython 1.79+

### Performance
- Efficient NumPy-based algorithms
- Optimized motif boundary calculations
- Batch processing capabilities
- Memory-efficient sequence handling 